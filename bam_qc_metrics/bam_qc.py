#! /usr/bin/env python3

"""Main class to compute BAM QC metrics"""

import csv, json, os, re, pybedtools, pysam, sys, tempfile

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
    READ_1_INDEX = 0
    READ_2_INDEX = 1
    READ_UNKNOWN_INDEX = 2

    def __init__(self, bam_path, target_path, expected_insert_max, metadata_path=None,
                 mark_duplicates_path=None, trim_quality=None, sample_rate=None, tmpdir=None):
        # define instance variables
        self.bedtools_metrics = None
        self.custom_metrics = None
        self.expected_insert_max = expected_insert_max
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
        self.metadata = None
        self.qc_input_bam_path = None
        self.qual_fail_reads = None
        self.sample_rate = 1
        self.samtools_metrics = None
        self.target_path = target_path
        self.trim_quality = trim_quality
        self.tmp_object = None
        self.unmapped_excluded_reads = None
        # set up temporary directory
        if tmpdir==None:
            self.tmp_object = tempfile.TemporaryDirectory(prefix='bam_qc_')
            self.tmpdir = self.tmp_object.name
        else:
            self.tmpdir = tmpdir
        # apply downsampling (if any)
        if sample_rate != None and sample_rate > 1:
            self.sample_rate = sample_rate
            unfiltered_bam_path = self.generate_downsampled_bam(bam_path, self.sample_rate)
        else:
            unfiltered_bam_path = bam_path
        # apply quality filter (if any)
        excluded_reads_path = os.path.join(self.tmpdir, 'excluded.bam')
        if self.trim_quality != None and self.trim_quality > 0:
            self.qc_input_bam_path = os.path.join(self.tmpdir, 'filtered.bam')
            pysam.view(unfiltered_bam_path,
                       '-b',
                       '-q', str(self.trim_quality),
                       '-o', self.qc_input_bam_path,
                       '-U', excluded_reads_path,
                       catch_stdout=False)
            self.qual_fail_reads = int(pysam.view('-c', excluded_reads_path).strip())
            # unmapped reads will fail the mapping quality filter, by definition
            # so if the quality filter is applied, find unmapped total from the excluded reads
            self.unmapped_excluded_reads = self.find_unmapped_reads(excluded_reads_path)
            if self.sample_rate != 1:
                os.remove(unfiltered_bam_path) # downsampled, unfiltered file is no longer needed
        else:
            self.qc_input_bam_path = unfiltered_bam_path
            self.qual_fail_reads = 0
        # read required keys from metadata (if any)
        if metadata_path != None:
            with open(metadata_path) as f: raw_metadata = json.loads(f.read())
            self.metadata = {key: raw_metadata.get(key) for key in self.METADATA_KEYS}
        else:
            sys.stderr.write("Metadata file not given, using empty defaults\n")
            self.metadata = {key: None for key in self.METADATA_KEYS}
        # find metrics
        if mark_duplicates_path != None:
            self.mark_duplicates_metrics = self.read_mark_duplicates_metrics(mark_duplicates_path)
        self.bedtools_metrics = self.evaluate_bedtools_metrics()
        self.samtools_metrics = self.evaluate_samtools_metrics()
        read1_hist = self.samtools_metrics['read 1 length histogram']
        read2_hist = self.samtools_metrics['read 1 length histogram']
        read1_length = max(read1_hist.keys()) if len(read1_hist) > 0 else 0
        read2_length = max(read2_hist.keys()) if len(read2_hist) > 0 else 0
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
        bamBedTool = pybedtools.BedTool(self.qc_input_bam_path)
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
        for read in pysam.AlignmentFile(self.qc_input_bam_path, 'rb').fetch(until_eof=True):
            if not read.has_tag('MD'):
                metrics['readsMissingMDtags'] += 1
            cycle = 1
            read_index = None
            if read.is_read1: read_index = self.READ_1_INDEX
            elif read.is_read2: read_index = self.READ_2_INDEX
            else: read_index = self.READ_UNKNOWN_INDEX
            if read.cigartuples != None:
                if read.is_reverse: cigar_list = reversed(read.cigartuples)
                else: cigar_list = read.cigartuples
                for (op, length) in cigar_list:
                    if op in op_names:
                        if op == 4: metrics['soft clip bases'] += length
                        elif op == 5: metrics['hard clip bases'] += length
                        for i in range(length):
                            key = 'read %s %s by cycle' % (read_names[read_index], op_names[op])
                            metrics[key][cycle] += 1
                            if op in consumes_query: cycle += 1
                    elif op in consumes_query:
                        cycle += length
            if read_index == self.READ_UNKNOWN_INDEX:
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
        result = pysam.stats(self.qc_input_bam_path)
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
                if samtools_key not in key_map:
                    continue
                elif samtools_key == 'reads unmapped':
                    if self.unmapped_excluded_reads != None:
                        # filtering in effect; unmapped reads excluded by alignment quality filter
                        val = self.unmapped_excluded_reads
                        if int(row[2]) != 0:
                            # This *should* never happen, but just in case it does
                            msg = "WARNING: %s reads were unmapped AND had alignment scores "+\
                                  "> %d; not counted in 'unmapped reads' total. Inconsistent "+\
                                  "data in BAM input?\n" % (row[2], self.trim_quality)
                            sys.stderr.write(msg)
                    else:
                        # no quality filtering; use unmapped reads count from main input
                        val = int(row[2])
                elif samtools_key in float_keys:
                    val = float(row[2])
                else:
                    val = int(row[2])
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

    def find_unmapped_reads(self, bam_path):
        result = pysam.stats(bam_path)
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", result)),
            delimiter="\t"
        )
        unmapped = 0
        for row in reader:
            if row[0]=='SN' and row[1]=='reads unmapped:':
                unmapped = int(row[2])
                break
        return unmapped

    def fq_stats(self, rows):
        '''
        Compute quality metrics from either FFQ or LFQ entries in samtools stats:
            - Mean quality by cycle
            - Quality histogram
        '''
        meanByCyc = {}
        histogram = {}
        if len(rows) > 0:
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

    def generate_downsampled_bam(self, bam_path, sample_rate):
        '''
        Write a temporary downsampled BAM file for all subsequent input

        This is fully deterministic -- for sample rate N, takes every (N*2)th pair of reads
        '''
        if sample_rate == 1:
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
        interval = sample_rate * 2
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
            if re.match('## METRICS CLASS\s.*picard\.sam\.DuplicationMetrics$', line):
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
                params = (input_path, section, line_count)
                msg = "Failed to parse duplicate metrics path %s, section %d, line %d" % params
                raise ValueError(msg)
        if len(keys) == len(values) + 1 and keys[-1] == 'ESTIMATED_LIBRARY_SIZE':
            # field is empty (no trailing \t) for low coverage; append a default value
            values.append(0)
        elif len(keys) != len(values):
            # otherwise, mismatched key/value totals are an error
            raise ValueError("Numbers of keys and values in %s do not match" % input_path)
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
        output['target file'] = os.path.split(self.target_path)[-1]
        output['sample rate'] = self.sample_rate
        output['insert max'] = self.expected_insert_max
        if out_path != '-':
            out_file = open(out_path, 'w')
        else:
            out_file = sys.stdout
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()


