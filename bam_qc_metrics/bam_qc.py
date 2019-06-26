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
    
    def __init__(self, bam_path, target_path, metadata_path=None, mark_duplicates_path=None,
                 trim_quality=None, expected_insert_max=None, sample_rate=None):
        with open(metadata_path) as f: self.metadata = json.loads(f.read())
        if sample_rate != None:
            self.sample_rate = sample_rate
            self.tmp = tempfile.TemporaryDirectory(prefix='bam_qc_downsampled_')
            self.bam_path = self.generate_downsampled_bam(bam_path)
        else:
            self.sample_rate = 1
            self.tmp_dir = None
            self.bam_path = bam_path
        self.target_path = target_path
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
        self.trim_quality = trim_quality
        if expected_insert_max != None:
            self.expected_insert_max = expected_insert_max
        else:
            self.expected_insert_max = DEFAULT_INSERT_MAX
        self.bedtools_metrics = self.evaluate_bedtools_metrics()
        self.samtools_metrics = self.evaluate_samtools_metrics()
        read1_length = max(self.samtools_metrics['read 1 length histogram'].keys())
        read2_length = max(self.samtools_metrics['read 2 length histogram'].keys())
        self.custom_metrics = self.evaluate_custom_metrics(read1_length, read2_length)

    def count_mapped_abnormally_far(self, insert_size_histogram):
        count = 0
        for key in insert_size_histogram.keys():
            if key >= self.expected_insert_max:
                count += insert_size_histogram[key]
        return count

    def evaluate_bedtools_metrics(self):
        metrics = {}
        bamBedTool = pybedtools.BedTool(self.bam_path)
        targetBedTool = pybedtools.BedTool(self.target_path)
        metrics['number of targets'] = targetBedTool.count()
        metrics['reads on target'] = len(bamBedTool.intersect(self.target_path))
        size = 0
        with open(self.target_path, newline='') as bedfile:
            reader = csv.reader(bedfile, delimiter="\t")
            for row in reader:
                size += int(row[2]) - int(row[1])
        metrics['target size'] = size
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
            'qual cut': self.trim_quality,
            'qual fail reads': 0
        }
        read_names = ['1', '2', '?']
        read_lengths = [read1_length, read2_length, max(read1_length, read2_length)]
        for op_name in op_names.values():
            for i in range(len(read_names)):
                key = 'read %s %s by cycle' % (read_names[i], op_name)
                metrics[key] = {j:0 for j in range(1, read_lengths[i]+1) }
        # iterate over the BAM file
        consumes_query = set([0,1,4,7,8]) # CIGAR op indices which increment the query cycle
        for read in pysam.AlignmentFile(self.bam_path, 'rb').fetch(until_eof=True):
            if self.trim_quality != None and read.mapping_quality < self.trim_quality:
                metrics['qual fail reads'] += 1
                continue
            if not read.has_tag('MD'):
                metrics['readsMissingMDtags'] += 1
            cycle = 0
            for (op, length) in read.cigartuples:
                if op in op_names:
                    for i in range(length):
                        if op in consumes_query: cycle += 1
                        if op == 4: metrics['soft clip bases'] += length
                        elif op == 5: metrics['hard clip bases'] += length
                        if read.is_read1: read_index = 0
                        elif read.is_read2: read_index = 1
                        else: read_index = 2
                        key = 'read %s %s by cycle' % (read_names[read_index], op_names[op])
                        metrics[key][cycle] += 1
                elif op in consumes_query:
                    cycle += length
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
            'insert size average': 'insert mean',
            'insert size standard deviation': 'insert stdev',
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
        if self.trim_quality == None:
            result = pysam.stats(self.bam_path)
        else:
            result = pysam.stats("-q", str(self.trim_quality), self.bam_path)
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

    def fq_stats(self, rows, precision=6):
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
            meanByCyc[cycle] = round(float(total) / count, precision) if count > 0 else 0
        return (meanByCyc, histogram)

    def generate_downsampled_bam(self, bam_path):
        '''
        Write a temporary downsampled BAM file for all subsequent input
        '''
        if self.sample_rate == 1:
            sys.stderr.write("Sample rate = 1, omitting down sampling\n")
            return bam_path
        sorted_bam_path = os.path.join(self.tmp.name, 'sorted.bam')
        sampled_bam_path = os.path.join(self.tmp.name, 'downsampled.bam')
        # ensure file is sorted by name, so pairs are together
        pysam.sort('-n', '-o', sorted_bam_path, bam_path)
        bam_in = pysam.AlignmentFile(sorted_bam_path)
        bam_out = pysam.AlignmentFile(sampled_bam_path, 'wb', template=bam_in)
        count = 0
        # sample two reads at a time -- should be read 1 and read 2, if read is paired
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
        mean_rl = float(total) / count if count > 0 else 0
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
        output['mark duplicates'] = self.mark_duplicates_metrics
        output['target'] = os.path.abspath(self.target_path)
        output['sample rate'] = self.sample_rate
        output['quality cutoff'] = self.trim_quality
        output['insert max'] = self.expected_insert_max
        if out_path != '-':
            out_file = open(out_path, 'w')
        else:
            out_file = sys.stdout
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()

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
    if args.trim_quality != None:
        valid = validate_positive_integer(args.trim_quality, 'Quality score')
    if args.insert_max != None:
        valid = validate_positive_integer(args.insert_max, 'Max insert size')
    if args.sample_rate != None:
        valid = validate_positive_integer(args.sample_rate, 'Downsampling rate')
    for path_arg in (args.bam, args.target, args.metadata, args.mark_duplicates):
        if path_arg == None:
            continue
        if not os.path.exists(path_arg):
            sys.stderr.write("ERROR: Path %s does not exist.\n" % path_arg)
            valid = False
        elif not os.path.isfile(path_arg):
            sys.stderr.write("ERROR: Path %s is not a file.\n" % path_arg)
            valid = False
        elif not os.access(path_arg, os.R_OK):
            sys.stderr.write("ERROR: Path %s is not readable.\n" % path_arg)
            valid = False
    if args.out != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        if not os.path.exists(parent_path):
            sys.stderr.write("ERROR: Parent directory of %s does not exist.\n" % args.out)
            valid = False
        elif not os.path.isdir(parent_path):
            sys.stderr.write("ERROR: Parent of %s is not a directory.\n" % args.out)
            valid = False
        elif not os.access(parent_path, os.W_OK):
            sys.stderr.write("ERROR: Parent directory of %s is not writable.\n" % args.out)
            valid = False
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
    args = parser.parse_args()
    if not validate_args(args): exit(1)
    trim_quality = None if args.trim_quality == None else int(args.trim_quality)
    insert_max = None if args.insert_max == None else int(args.insert_max)
    sample_rate = None if args.sample_rate == None else int(args.sample_rate)
    qc = bam_qc(args.bam, args.target, args.metadata, args.mark_duplicates, trim_quality, insert_max, sample_rate)
    qc.write_output(args.out)

if __name__ == "__main__":
    main()
