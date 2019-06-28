#! /usr/bin/env python3

"""Main script and class to compute BAM QC metrics"""

import argparse, csv, json, os, re, pybedtools, pysam, sys, tempfile

DEFAULT_INSERT_MAX = 1500

class bam_qc:

    DOWNSAMPLE_WARNING_THRESHOLD = 1000
    METADATA_KEYS = [
        'barcode',
        'instrument',
        'lane',
        'library',
        'run name',
        'sample'
    ]
    PRECISION = 1 # number of decimal places for rounded output
    
    def __init__(self, bam_path, target_path, metadata_path=None, mark_duplicates_path=None,
                 trim_quality=None, expected_insert_max=None, sample_rate=None, tmpdir=None):
        self.target_path = target_path
        self.trim_quality = trim_quality
        self.tmp_object = None
        if tmpdir==None:
            self.tmp_object = tempfile.TemporaryDirectory(prefix='bam_qc_')
            self.tmpdir = self.tmp_object.name
        else:
            self.tmpdir = tmpdir
        # apply downsampling (if any)
        unfiltered_bam_path = None
        if sample_rate != None:
            self.sample_rate = sample_rate
            unfiltered_bam_path = self.generate_downsampled_bam(bam_path)
        else:
            self.sample_rate = 1
            unfiltered_bam_path = bam_path
        # apply quality filter (if any)
        excluded_reads_path = os.path.join(self.tmpdir, 'excluded.bam')
        if self.trim_quality:
            self.filtered_bam_path = os.path.join(self.tmpdir, 'filtered.bam')
            pysam.view(unfiltered_bam_path,
                       '-b',
                       '-q', str(self.trim_quality),
                       '-o', self.filtered_bam_path,
                       '-U', excluded_reads_path,
                       catch_stdout=False)
            self.qual_fail_reads = int(pysam.view('-c', excluded_reads_path).strip())
        else:
            self.filtered_bam_path = unfiltered_bam_path
            self.qual_fail_reads = 0
        # read required keys from metadata (if any)
        if metadata_path != None:
            with open(metadata_path) as f: raw_metadata = json.loads(f.read())
            self.metadata = {key: raw_metadata.get(key) for key in self.METADATA_KEYS}
        else:
            sys.stderr.write("WARNING: Metadata file not given, using empty defaults\n")
            self.metadata = {key: None for key in self.METADATA_KEYS}
        # find metrics
        self.mark_duplicates_metrics = {
            "ESTIMATED_LIBRARY_SIZE": None,
            "HISTOGRAM": {},
            "LIBRARY": None,
            "PERCENT_DUPLICATION": None,
            "READ_PAIRS_EXAMINED": None,
            "READ_PAIR_DUPLICATES": None,
            "READ_PAIR_OPTICAL_DUPLICATES": None,
            "UNMAPPED_READS": None,
            "UNPAIRED_READS_EXAMINED": None,
            "UNPAIRED_READ_DUPLICATES": None
        }
        if mark_duplicates_path != None:
            self.mark_duplicates_metrics = self.read_mark_duplicates_metrics(mark_duplicates_path)
        if expected_insert_max != None:
            self.expected_insert_max = expected_insert_max
        else:
            self.expected_insert_max = DEFAULT_INSERT_MAX
        self.bedtools_metrics = self.evaluate_bedtools_metrics()
        self.samtools_metrics = self.evaluate_samtools_metrics()
        read1_length = max(self.samtools_metrics['read 1 length histogram'].keys())
        read2_length = max(self.samtools_metrics['read 2 length histogram'].keys())
        self.custom_metrics = self.evaluate_custom_metrics(read1_length, read2_length)

    def cleanup(self):
        '''
        Temporary directory object (if any) will be automatically deleted when bam_qc object
        is out of scope. This method allows explicit cleanup, eg. to avoid warnings in the
        test harness.
        '''
        if self.tmp_object != None:
            self.tmp_object.cleanup()
        else:
            sys.stderr.write("Omitting cleanup for user-specified temporary directory %s\n" % self.tmpdir)

    def count_mapped_abnormally_far(self, insert_size_histogram):
        count = 0
        for key in insert_size_histogram.keys():
            if key >= self.expected_insert_max:
                count += insert_size_histogram[key]
        return count

    def evaluate_bedtools_metrics(self):
        metrics = {}
        bamBedTool = pybedtools.BedTool(self.filtered_bam_path)
        targetBedTool = pybedtools.BedTool(self.target_path)
        metrics['number of targets'] = targetBedTool.count()
        metrics['reads on target'] = len(bamBedTool.intersect(self.target_path))
        size = 0
        with open(self.target_path, newline='') as bedfile:
            reader = csv.reader(bedfile, delimiter="\t")
            for row in reader:
                size += int(row[2]) - int(row[1])
        metrics['total target size'] = size
        # TODO add bedtools coverage metrics?
        #coverage = targetBedTool.coverage(self.bam_path)
        #print(coverage)
        return metrics
        
    def evaluate_custom_metrics(self, read1_length, read2_length):
        '''
        Iterate over the BAM file to compute custom metrics
        Processes CIGAR strings; see p. 7 of https://samtools.github.io/hts-specs/SAMv1.pdf
        '''
        # Relevant CIGAR operations
        op_names = {
            0: 'aligned',
            1: 'insertion',
            2: 'deletion',
            4: 'soft clip',
            5: 'hard clip',
            8: 'mismatch'
        }
        # initialize the metrics data structure
        metrics = {
            'hard clip bases': 0,
            'soft clip bases': 0,
            'readsMissingMDtags': 0,
        }
        read_names = ['1', '2', '?']
        read_lengths = [read1_length, read2_length, max(read1_length, read2_length)]
        for op_name in op_names.values():
            for i in range(len(read_names)):
                key = 'read %s %s by cycle' % (read_names[i], op_name)
                metrics[key] = {j:0 for j in range(1, read_lengths[i]+1) }
        # metrics for unknown reads -- equivalents for reads 1 and 2 are derived from samtools
        ur_count = 0
        ur_length_total = 0
        ur_length_histogram = {}
        ur_quality_by_cycle = {i : 0 for i in range(1, read_lengths[2]+1)}
        ur_total_by_cycle = {i : 0 for i in range(1, read_lengths[2]+1)}
        ur_quality_histogram = {}
        # iterate over the BAM file
        consumes_query = set([0,1,4,7,8]) # CIGAR op indices which increment the query cycle
        for read in pysam.AlignmentFile(self.filtered_bam_path, 'rb').fetch(until_eof=True):
            if not read.has_tag('MD'):
                metrics['readsMissingMDtags'] += 1
            cycle = 1
            read_index = None
            if read.is_read1: read_index = 0
            elif read.is_read2: read_index = 1
            else: read_index = 2
            for (op, length) in read.cigartuples:
                if op in op_names:
                    if op == 4: metrics['soft clip bases'] += length
                    elif op == 5: metrics['hard clip bases'] += length
                    for i in range(length):
                        key = 'read %s %s by cycle' % (read_names[read_index], op_names[op])
                        metrics[key][cycle] += 1
                        if op in consumes_query: cycle += 1
                elif op in consumes_query:
                    cycle += length
            if read_index == 2:
                ur_count += 1
                ur_length = read.query_length
                ur_length_total += ur_length
                if ur_length in ur_length_histogram:
                    ur_length_histogram[ur_length] += 1
                else:
                    ur_length_histogram[ur_length] = 1
                for i in range(read.query_length):
                    q = read.query_qualities[i]
                    ur_quality_by_cycle[i+1] += q
                    ur_total_by_cycle[i+1] += 1
                    if q in ur_quality_histogram:
                        ur_quality_histogram[q] += 1
                    else:
                        ur_quality_histogram[q] = 1
        metrics['read ? average length'] = round(float(ur_length_total) / ur_count, self.PRECISION) if ur_count > 0 else None
        metrics['read ? length histogram'] = ur_length_histogram
        for cycle in ur_quality_by_cycle.keys():
            quality = ur_quality_by_cycle[cycle]
            total = ur_total_by_cycle[cycle]
            ur_quality_by_cycle[cycle] = round(float(quality)/total, self.PRECISION) if total > 0 else 0
        metrics['read ? quality by cycle'] = ur_quality_by_cycle
        metrics['read ? quality histogram'] = ur_quality_histogram
        return metrics

    def evaluate_samtools_metrics(self):
        '''Process metrics derived from samtools output'''
        # summary numbers (SN) fields denoted in float_keys are floats; integers otherwise
        float_keys = set([
            'error rate',
            'average quality',
            'insert size average',
            'insert size standard deviation',
            'percentage of properly paired reads (%)'
        ])
        # map from SN field names to output keys
        key_map = {
            'bases mapped (cigar)': 'bases mapped',
            'average length': 'average read length',
            'insert size average': 'insert size average',
            'insert size standard deviation': 'insert size standard deviation',
            'reads mapped': 'mapped reads',
            'reads mapped and paired': 'reads mapped and paired',
            'mismatches': 'mismatched bases',
            # 'reads paired' > 0 implies 'paired end' == True
            'reads paired': 'paired reads',
            'pairs on different chromosomes': 'pairsMappedToDifferentChr',
            'reads properly paired': 'properly paired reads',
            'raw total sequences': 'total reads',
            'reads unmapped': 'unmapped reads',
            'non-primary alignments': 'non primary reads',
        }
        metrics = {}
        metrics['inserted bases'] = 0
        metrics['deleted bases'] = 0
        metrics['insert size histogram'] = {}
        labels_to_store = set(['FFQ', 'FRL', 'LFQ', 'LRL', 'RL'])
        stored = {} # store selected rows for later processing
        for label in labels_to_store: stored[label] = []
        result = pysam.stats(self.filtered_bam_path)
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", result)),
            delimiter="\t"
        )
        for row in reader:
            if row[0] in labels_to_store:
                stored[row[0]].append(row)
            elif row[0] == 'ID':
                metrics['inserted bases'] += int(row[1]) * int(row[2])
                metrics['deleted bases'] += int(row[1]) * int(row[3])
            elif row[0] == 'IS':
                metrics['insert size histogram'][int(row[1])] = int(row[2])
            elif row[0] == 'SN':
                samtools_key = re.sub(':$', '', row[1])
                if samtools_key not in key_map: continue
                if samtools_key in float_keys: val = float(row[2])
                else: val = int(row[2])
                metrics[key_map[samtools_key]] = val
        metrics['average read length'] = self.mean_read_length(stored['RL'])
        metrics['paired end'] = metrics['paired reads'] > 0
        metrics['read 1 average length'] = self.mean_read_length(stored['FRL'])
        metrics['read 2 average length'] = self.mean_read_length(stored['LRL'])
        metrics['read 1 length histogram'] = self.read_length_histogram(stored['FRL'])
        metrics['read 2 length histogram'] = self.read_length_histogram(stored['LRL'])
        (ffq_mean_by_cycle, ffq_histogram) = self.fq_stats(stored['FFQ'])
        (lfq_mean_by_cycle, lfq_histogram) = self.fq_stats(stored['LFQ'])
        metrics['read 1 quality by cycle'] = ffq_mean_by_cycle
        metrics['read 1 quality histogram'] = ffq_histogram
        metrics['read 2 quality by cycle'] = lfq_mean_by_cycle
        metrics['read 2 quality histogram'] = lfq_histogram
        metrics['pairsMappedAbnormallyFar'] = self.count_mapped_abnormally_far(metrics['insert size histogram'])
        return metrics

    def fq_stats(self, rows):
        '''
        Compute quality metrics from either FFQ or LFQ entries in samtools stats:
            - Mean quality by cycle
            - Quality histogram
        '''
        meanByCyc = {}
        max_width = max([len(row) for row in rows])
        histogram = {q: 0 for q in range(max_width-2)}
        for row in rows:
            cycle = int(row[1])
            counts = [int(n) for n in row[2:]]
            total = 0
            count = 0
            for qscore in range(len(counts)):
                total += counts[qscore]*qscore
                count += counts[qscore]
                histogram[qscore] += counts[qscore]
            meanByCyc[cycle] = round(float(total) / count, self.PRECISION) if count > 0 else 0
        return (meanByCyc, histogram)

    def generate_downsampled_bam(self, bam_path):
        '''
        Write a temporary downsampled BAM file for all subsequent input

        This is fully deterministic -- for sample rate N, takes every (N*2)th pair of reads
        '''
        if self.sample_rate == 1:
            sys.stderr.write("Sample rate = 1, omitting down sampling\n")
            return bam_path
        sorted_bam_path = os.path.join(self.tmpdir, 'sorted.bam')
        sampled_bam_path = os.path.join(self.tmpdir, 'downsampled.bam')
        # ensure file is sorted by name, so pairs are together
        pysam.sort('-n', '-o', sorted_bam_path, bam_path)
        bam_in = pysam.AlignmentFile(sorted_bam_path)
        bam_out = pysam.AlignmentFile(sampled_bam_path, 'wb', template=bam_in)
        count = 0
        # sample two reads at a time -- should be read 1 and read 2, if read is paired
        # first read is odd-numbered (should be read 1), second is even-numbered (read 2)
        interval = self.sample_rate * 2
        sample_next = False
        sampled = 0
        for read in bam_in:
            count += 1
            if (count + 1) % interval == 0:
                bam_out.write(read)
                sampled += 1
                sample_next = True
            elif sample_next:
                bam_out.write(read)
                sample_next = False
                sampled += 1
        bam_in.close()
        bam_out.close()
        if sampled < self.DOWNSAMPLE_WARNING_THRESHOLD:
            sys.stderr.write("WARNING: Only %i reads remain after downsampling\n" % sampled)
        return sampled_bam_path

    def mean_read_length(self, rows):
        '''Process RL (read length), FRL (first RL) or LRL (last RL) rows for mean read length'''
        total = 0
        count = 0
        for row in rows:
            length = int(row[1])
            length_count = int(row[2])
            total += length * length_count
            count += length_count
        mean_rl = round(float(total) / count, self.PRECISION) if count > 0 else 0
        return mean_rl

    def read_length_histogram(self, rows):
        '''Process RL (read length), FRL (first RL) or LRL (last RL) rows for read length histogram'''
        histogram = {}
        for row in rows:
            histogram[int(row[1])] = int(row[2])
        return histogram
    
    def read_mark_duplicates_metrics(self, input_path):
        section = 0
        line_count = 0
        with open(input_path) as f: lines = f.readlines()
        keys = []
        values = []
        hist = {}
        for line in lines:
            line_count += 1
            line = line.strip()
            if re.match('## METRICS CLASS\s+net\.sf\.picard\.sam\.DuplicationMetrics', line):
                section += 1
            elif section == 1:
                keys = re.split("\t", line)
                section += 1
            elif section == 2:
                values = re.split("\t", line)
                section += 1
            elif section == 3 and re.match('## HISTOGRAM', line):
                section += 1
            elif section == 4 and re.match("BIN\tVALUE", line):
                section += 1
            elif section == 5 and re.match("[0-9]+\.[0-9]+\t[0-9]+\.{0,1}[0-9]*$", line):
                terms = re.split("\t", line)
                # JSON doesn't allow numeric dictionary keys, so hist_bin is stringified in output
                # but rounding removes the trailing '.0'
                hist_bin = round(float(terms[0])) if re.search('\.0$', terms[0]) else terms[0]
                hist[hist_bin] = float(terms[1])
            elif re.match('#', line) and section < 4 or line == '':
                continue
            else:
                params = (input_path, section, line)
                msg = "Failed to parse duplicate metrics path %s, section %d, line %d" % params
                raise ValueError(msg)
        if len(keys) != len(values):
            raise ValueError("Key and value lists from %s are of unequal length" % input_path)
        metrics = {}
        for i in range(len(keys)):
            if keys[i] == 'PERCENT_DUPLICATION':
                metrics[keys[i]] = float(values[i])
            elif keys[i] == 'LIBRARY':
                metrics[keys[i]] = values[i]
            else:
                metrics[keys[i]] = int(values[i])
        metrics['HISTOGRAM'] = hist
        return metrics
        
    def write_output(self, out_path):
        output = {}
        for key in self.METADATA_KEYS:
            output[key] = self.metadata.get(key)
        for key in self.bedtools_metrics.keys():
            output[key] = self.bedtools_metrics.get(key)
        for key in self.samtools_metrics.keys():
            output[key] = self.samtools_metrics.get(key)
        for key in self.custom_metrics.keys():
            output[key] = self.custom_metrics.get(key)
        output['qual cut'] = self.trim_quality
        output['qual fail reads'] = self.qual_fail_reads
        output['mark duplicates'] = self.mark_duplicates_metrics
        output['target file'] = os.path.abspath(self.target_path)
        output['sample rate'] = self.sample_rate
        output['insert max'] = self.expected_insert_max
        if out_path != '-':
            out_file = open(out_path, 'w')
        else:
            out_file = sys.stdout
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()


def validate_input_file(path_arg):
    valid = True
    if not os.path.exists(path_arg):
        sys.stderr.write("ERROR: Path %s does not exist.\n" % path_arg)
        valid = False
    elif not os.path.isfile(path_arg):
        sys.stderr.write("ERROR: Path %s is not a file.\n" % path_arg)
        valid = False
    elif not os.access(path_arg, os.R_OK):
        sys.stderr.write("ERROR: Path %s is not readable.\n" % path_arg)
        valid = False
    return valid

def validate_output_dir(dir_path):
    valid = True
    if not os.path.exists(dir_path):
        sys.stderr.write("ERROR: Directory %s does not exist.\n" % dir_path)
        valid = False
    elif not os.path.isdir(dir_path):
        sys.stderr.write("ERROR: Path %s is not a directory.\n" % dir_path)
        valid = False
    elif not os.access(dir_path, os.W_OK):
        sys.stderr.write("ERROR: Directory %s is not writable.\n" % dir_path)
        valid = False
    return valid

def validate_positive_integer(arg, name):
    param = None
    valid = True
    try:
        param = int(arg)
    except ValueError:
        sys.stderr.write("ERROR: %s must be an integer.\n" % name)
        valid = False
    if param != None and param < 0:
        sys.stderr.write("ERROR: %s cannot be negative.\n" % name)
        valid = False
    return valid

def validate_args(args):
    valid = True
    # flip valid from True to False if a check is failed; never flip back to True
    if args.trim_quality != None:
        valid = validate_positive_integer(args.trim_quality, 'Quality score')
    if args.insert_max != None:
        valid = valid and validate_positive_integer(args.insert_max, 'Max insert size')
    if args.sample_rate != None:
        valid = valid and validate_positive_integer(args.sample_rate, 'Downsampling rate')
    for path_arg in (args.bam, args.target, args.metadata, args.mark_duplicates):
        if path_arg != None:
            valid = valid and validate_input_file(path_arg)
    if args.out != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        valid = valid and validate_output_dir(parent_path)
    if args.temp_dir != None:
        valid = valid and validate_output_dir(args.temp_dir)
    return valid

def main():
    parser = argparse.ArgumentParser(description='QC for BAM files.')
    parser.add_argument('-b', '--bam', metavar='PATH', required=True,
                        help='Path to input BAM file. Required.')
    parser.add_argument('-d', '--mark-duplicates', metavar='PATH',
                        help='Path to text file output by Picard MarkDuplicates. Optional.')
    parser.add_argument('-i', '--insert-max', metavar='INT',
                        help='Maximum expected value for insert size; higher values will be counted as abnormal. '+\
                        'Optional; default = %i.' % DEFAULT_INSERT_MAX)
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata. Optional.')
    parser.add_argument('-o', '--out', metavar='PATH', required=True,
                        help='Path for JSON output, or - for STDOUT. Required.')
    parser.add_argument('-q', '--trim-quality', metavar='QSCORE',
                        help='Samtools threshold for trimming alignment quality. Optional.')
    parser.add_argument('-s', '--sample-rate', metavar='INT',
                        help='Sample every Nth read, where N is the argument. Optional, defaults to 1 (no sampling).')
    parser.add_argument('-t', '--target', metavar='PATH',
                        help='Path to target BED file, containing targets to calculate coverage '+\
                        'against. Optional; if given, must be sorted in same order as BAM file.')
    parser.add_argument('-T', '--temp-dir', metavar='PATH', help='Directory for temporary output files; optional, defaults to /tmp')
    args = parser.parse_args()
    if not validate_args(args): exit(1)
    trim_quality = None if args.trim_quality == None else int(args.trim_quality)
    insert_max = None if args.insert_max == None else int(args.insert_max)
    sample_rate = None if args.sample_rate == None else int(args.sample_rate)
    qc = bam_qc(args.bam, args.target, args.metadata, args.mark_duplicates, trim_quality, insert_max, sample_rate, args.temp_dir)
    qc.write_output(args.out)

if __name__ == "__main__":
    main()
