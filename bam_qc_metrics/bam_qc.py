#! /usr/bin/env python3

"""Main class to compute BAM QC metrics"""

import csv, json, os, re, pybedtools, pysam, sys, tempfile

class base:

    """
    Class for shared constants
    """

    PRECISION = 1 # number of decimal places for rounded output
    FINE_PRECISION = 3 # finer precision, eg. for reads_per_start_point output
    READ_1_LENGTH_KEY = 'read 1'
    READ_2_LENGTH_KEY = 'read 2'
    MAX_READ_LENGTH_KEY = 'max_read_length'
    UNMAPPED_READS_KEY = 'unmapped reads'

class bam_qc(base):

    DEFAULT_MARK_DUPLICATES_METRICS =  {
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
    DOWNSAMPLE_WARNING_THRESHOLD = 1000
    METADATA_KEYS = [
        'barcode',
        'instrument',
        'lane',
        'library',
        'run name',
        'sample'
    ]
    RANDOM_SEED = 42

    def __init__(self, bam_path, target_path, expected_insert_max, metadata_path=None,
                 mark_duplicates_path=None, n_as_mismatch=False, skip_below_mapq=None,
                 reference=None, sample_rate=None, tmpdir=None, verbose=True):
        # define instance variables
        self.bedtools_metrics = None
        self.custom_metrics = None
        self.expected_insert_max = expected_insert_max
        self.fast_metrics = None
        self.mark_duplicates_metrics = self.DEFAULT_MARK_DUPLICATES_METRICS
        self.metadata = None
        self.mismatches_by_read = None
        self.qual_fail_reads = None
        self.reference = reference
        self.sample_rate = None
        self.skip_below_mapq = skip_below_mapq
        self.target_path = target_path
        self.tmp_object = None
        self.tmpdir = None
        self.unmapped_excluded_reads = None
        self.verbose = verbose # if False, suppress non-critical messages to stderr
        (self.tmpdir, self.tmp_object) = self.setup_tmpdir(tmpdir)
        # read metadata and MarkDuplicates metrics
        self.metadata = self.read_metadata(metadata_path)
        if mark_duplicates_path != None:
            self.mark_duplicates_metrics = self.read_mark_duplicates_metrics(mark_duplicates_path)
        # apply quality filter (if any); get input path for downsampling
        result = self.apply_mapq_filter(bam_path)
        (fast_finder_input_path, self.qual_fail_reads, self.unmapped_excluded_reads) = result
        # find 'fast' metrics on full dataset -- after filtering, before downsampling
        fast_finder = fast_metric_finder(fast_finder_input_path,
                                         self.reference,
                                         self.expected_insert_max)
        self.fast_metrics = self.update_unmapped_count(fast_finder.metrics)
        # apply downsampling (if any)
        if sample_rate != None and sample_rate > 1:
            self.sample_rate = sample_rate
            slow_finder_input_path = self.generate_downsampled_bam(fast_finder_input_path,
                                                                   self.sample_rate)
        else:
            self.sample_rate = 1
            slow_finder_input_path = fast_finder_input_path
        # find 'slow' metrics on (maybe) downsampled dataset
        slow_finder = slow_metric_finder(slow_finder_input_path,
                                         self.target_path,
                                         fast_finder.read_length_summary())
        self.slow_metrics = slow_finder.metrics

    def apply_mapq_filter(self, bam_path):
        """
        Write a BAM file filtered by alignment quality
        Also report the number of reads failing quality filters, and number of unmapped reads
        """
        filtered_bam_path = None
        if self.mapq_filter_is_active():
            excluded_by_mapq_path = os.path.join(self.tmpdir, 'excluded.bam')
            included_by_mapq_path = os.path.join(self.tmpdir, 'included.bam')
            pysam.view(bam_path,
                       '-b',
                       '-q', str(self.skip_below_mapq),
                       '-o', included_by_mapq_path,
                       '-U', excluded_by_mapq_path,
                       catch_stdout=False)
            total_failed_reads = int(pysam.view('-c', excluded_by_mapq_path).strip())
            # unmapped reads will fail the mapping quality filter, by definition
            # so if the quality filter is applied, find unmapped total from the excluded reads
            unmapped_raw = pysam.view('-c', '-f', '4', excluded_by_mapq_path)
            total_unmapped_and_excluded = int(unmapped_raw.strip())
            filtered_bam_path = included_by_mapq_path
        else:
            total_failed_reads = 0
            total_unmapped_and_excluded = 0
            filtered_bam_path = bam_path
        return (filtered_bam_path, total_failed_reads, total_unmapped_and_excluded)

    def cleanup(self):
        """
        Temporary directory object (if any) will be automatically deleted when bam_qc object
        is out of scope. This method allows explicit cleanup, eg. to avoid warnings in the
        test harness.
        """
        if self.tmp_object != None:
            self.tmp_object.cleanup()
        elif self.verbose:
            sys.stderr.write("Omitting cleanup for user-specified "+\
                             "temporary directory %s\n" % self.tmpdir)

    def generate_downsampled_bam(self, bam_path, sample_rate):
        """
        Write a temporary downsampled BAM file
        """
        downsampled_path = None
        if sample_rate == 1:
            if self.verbose: sys.stderr.write("Sample rate = 1, omitting down sampling\n")
            downsampled_path = bam_path
        elif sample_rate < 1:
            raise ValueError("Sample rate cannot be less than 1")
        else:
            # We are sampling every Nth read
            # Argument to `samtools -s` is of the form RANDOM_SEED.DECIMAL_RATE
            # Eg. for random seed 42 and sample rate 4, 42 + (1/4) = 42.25
            downsampled_path = os.path.join(self.tmpdir, 'downsampled.bam')
            sample_decimal = round(1.0/sample_rate, self.FINE_PRECISION)
            sample_arg = str(self.RANDOM_SEED + sample_decimal)
            pysam.view('-u', '-s', sample_arg, '-o', downsampled_path, bam_path, catch_stdout=False)
            # sanity check on the downsampled file
            # TODO Report size of downsampled file in JSON output?
            sampled = int(pysam.view('-c', downsampled_path).strip())
            if sampled < self.DOWNSAMPLE_WARNING_THRESHOLD:
                sys.stderr.write("WARNING: Only %i reads remain after downsampling\n" % sampled)
        return downsampled_path

    def mapq_filter_is_active(self):
        return self.skip_below_mapq != None and self.skip_below_mapq > 0

    def read_metadata(self, metadata_path):
        metadata = None
        if metadata_path != None:
            with open(metadata_path) as f: raw_metadata = json.loads(f.read())
            metadata = {key: raw_metadata.get(key) for key in self.METADATA_KEYS}
        else:
            if self.verbose: sys.stderr.write("Metadata file not given, using empty defaults\n")
            metadata = {key: None for key in self.METADATA_KEYS}
        return metadata
    
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

    def setup_tmpdir(self, tmpdir):
        if tmpdir==None:
            tmp_object = tempfile.TemporaryDirectory(prefix='bam_qc_')
            tmpdir = self.tmp_object.name
        else:
            tmp_object = None
            tmpdir = tmpdir
        return (tmpdir, tmp_object)

    def update_unmapped_count(self, metrics):
        """
        Input 'metrics' was calculated after mapping quality filter was applied
        If mapping quality filter is in effect, all unmapped reads should have been filtered out
        Check this is so, and update the 'unmapped reads' metric entry
        """
        if self.mapq_filter_is_active():
            if metrics[self.UNMAPPED_READS_KEY] > 0:
                raise ValueError("Mapping quality filter is in effect, so all unmapped reads should "+\
                                 "have been removed from BAM input; but 'reads unmapped' field from "+\
                                 "samtools stats is non-zero.")
            else:
                metrics[self.UNMAPPED_READS_KEY] = self.unmapped_excluded_reads
        return metrics

    def write_output(self, out_path):
        output = {}
        for key in self.METADATA_KEYS:
            output[key] = self.metadata.get(key)
        for key in self.fast_metrics.keys():
            output[key] = self.fast_metrics.get(key)
        for key in self.slow_metrics.keys():
            output[key] = self.slow_metrics.get(key)
        output['alignment reference'] = self.reference
        output['insert max'] = self.expected_insert_max
        output['mark duplicates'] = self.mark_duplicates_metrics
        output['qual cut'] = self.skip_below_mapq
        output['qual fail reads'] = self.qual_fail_reads
        output['sample rate'] = self.sample_rate
        output['target file'] = os.path.split(self.target_path)[-1]
        if out_path != '-':
            out_file = open(out_path, 'w')
        else:
            out_file = sys.stdout
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()


class fast_metric_finder(base):

    """
    Find "fast" metric types which can be evaluated before downsampling
    Mostly uses native samtools stats, rather than custom metrics
    """
    
    # summary numbers (SN) fields denoted in float_keys are floats; integers otherwise
    FLOAT_KEYS = set([
        'error rate',
        'average quality',
        'insert size average',
        'insert size standard deviation',
        'percentage of properly paired reads (%)'
    ])
    
    def __init__(self, bam_path, reference, expected_insert_max):
        self.bam_path = bam_path
        self.reference = reference
        self.expected_insert_max = expected_insert_max
        # map from SN field names to output keys
        self.sn_key_map = {
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
            'reads unmapped': self.UNMAPPED_READS_KEY,
            'non-primary alignments': 'non primary reads',
        }
        self.samtools_stats = pysam.stats(self.bam_path)
        self.metrics = self.evaluate_all_metrics()
        # find secondary stats -- not for JSON export, but used to calculate subsequent metrics
        self.max_read_length = self.evaluate_max_read_length()
        read1_hist = self.metrics['read 1 length histogram']
        read2_hist = self.metrics['read 2 length histogram']
        self.read_1_length = max(read1_hist.keys()) if len(read1_hist) > 0 else 0
        self.read_2_length = max(read2_hist.keys()) if len(read2_hist) > 0 else 0

    def count_mapped_abnormally_far(self, insert_size_histogram):
        count = 0
        for key in insert_size_histogram.keys():
            if key >= self.expected_insert_max:
                count += insert_size_histogram[key]
        return count
        
    def evaluate_all_metrics(self):
        metrics = {}
        metric_subsets = [
            self.evaluate_mismatches(),
            self.read_length_and_quality_metrics(),
            self.misc_stats_metrics()
        ]   
        for metric_subset in metric_subsets:
            for key in metric_subset.keys():
                metrics[key] = metric_subset[key]
        return metrics

    def evaluate_max_read_length(self):
        """
        Find max read length from samtools stats output.
        Not part of JSON output, but needed to evaluate custom metrics.
        """
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", self.samtools_stats)),
            delimiter="\t"
        )
        max_read_length = 0
        for row in reader:
            if row[0] == 'SN' and row[1] == 'maximum length:':
                max_read_length = int(row[2])
                break
        return max_read_length

    def evaluate_mismatches(self):
        """
        Find mismatches by read using samtools stats with:
        -r to specify a reference and get mismatch counts
        -f/F to specify read1/read2/unknown
        Use 'samtools view' to find unknown reads
        Inefficient, but simpler than custom processing of MD tags

        Reference may be None -- if so, return 3 empty dictionaries
        """
        mismatches = {}
        # find read using flags; unknown read is neither R1 nor R2
        if self.reference != None:
            r1_mismatch = self.parse_mismatch(pysam.stats('-r', reference, '-f', '64', self.bam_path))
            r2_mismatch = self.parse_mismatch(pysam.stats('-r', reference, '-f', '128', self.bam_path))
            ur_mismatch = self.parse_mismatch(pysam.stats('-r', reference, '-F', '192', self.bam_path))
        else:
            r1_mismatch = {}
            r2_mismatch = {}
            ur_mismatch = {}
        mismatches['read 1 mismatch by cycle'] = r1_mismatch
        mismatches['read 2 mismatch by cycle'] = r2_mismatch
        mismatches['read ? mismatch by cycle'] = ur_mismatch
        return mismatches

    def fq_stats(self, rows):
        """
        Compute quality metrics from either FFQ or LFQ entries in samtools stats:
            - Mean quality by cycle
            - Quality histogram
        """
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

    def mean_read_length(self, rows):
        """Process RL (read length), FRL (first RL) or LRL (last RL) rows for mean read length"""
        total = 0
        count = 0
        for row in rows:
            length = int(row[1])
            length_count = int(row[2])
            total += length * length_count
            count += length_count
        mean_rl = round(float(total) / count, self.PRECISION) if count > 0 else 0
        return mean_rl
        
    def misc_stats_metrics(self):
        """
        Process the output from 'samtools stats' to derive miscellaneous metrics
        """
        metrics = {}
        metrics['inserted bases'] = 0
        metrics['deleted bases'] = 0
        metrics['insert size histogram'] = {}
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", self.samtools_stats)),
            delimiter="\t"
        )
        for row in reader:
            if row[0] == 'ID':
                metrics['inserted bases'] += int(row[1]) * int(row[2])
                metrics['deleted bases'] += int(row[1]) * int(row[3])
            elif row[0] == 'IS':
                metrics['insert size histogram'][int(row[1])] = int(row[2])
            elif row[0] == 'SN':
                samtools_key = re.sub(':$', '', row[1])
                if samtools_key not in self.sn_key_map:
                    continue
                elif samtools_key in self.FLOAT_KEYS:
                    val = float(row[2])
                else:
                    val = int(row[2])
                metrics[self.sn_key_map[samtools_key]] = val
        metrics['paired end'] = metrics['paired reads'] > 0
        abnormal_count = self.count_mapped_abnormally_far(metrics['insert size histogram'])
        metrics['pairsMappedAbnormallyFar'] = abnormal_count
        return metrics
    
    def parse_mismatch(self, samtools_stats):
        """
        Input is the string returned by pysam.stats
        Parse the mismatches by cycle; return an empty dictionary if no data found
        """
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", samtools_stats)),
            delimiter="\t"
        )
        mismatch_by_cycle = {}
        # parse rows with the 'MPC' identifier for mismatch-per-cycle
        # column 1 is cycle; 2 is number of N's; 3 and subsequent are mismatches by quality
        empty = False
        for row in reader:
            if row[0] == 'MPC':
                cycle = int(row[1])
                mismatches = 0
                if self.n_as_mismatch: start = 2
                else: start = 3
                for m in row[start:]:
                    mismatches += int(m)
                mismatch_by_cycle[cycle] = mismatches
            elif row[0] == 'SN' and row[1] == 'sequences:' and int(row[2]) == 0:
                empty = True
        if empty:
            # special case -- no reads found (eg. no unknown reads)
            # MPC returns 0 mismatches for cycle 1, but we want an empty dictionary
            mismatch_by_cycle = {}
        return mismatch_by_cycle

    def read_length_histogram(self, rows):
        """Process RL (read length), FRL (first RL) or LRL (last RL) rows for read length histogram"""
        histogram = {}
        for row in rows:
            histogram[int(row[1])] = int(row[2])
        return histogram

    def read_length_and_quality_metrics(self):
        """
        Process the output from 'samtools stats' to derive read length and quality metrics
        """
        metrics = {}
        stored = {label: [] for label in ['FFQ', 'FRL', 'LFQ', 'LRL', 'RL']}
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", self.samtools_stats)),
            delimiter="\t"
        )
        for row in reader:
            if row[0] in stored:
                stored[row[0]].append(row)
        metrics['average read length'] = self.mean_read_length(stored['RL'])
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
        return metrics

    def read_length_summary(self):
        summary = {
            self.READ_1_LENGTH_KEY: self.read_1_length,
            self.READ_2_LENGTH_KEY: self.read_2_length,
            self.MAX_READ_LENGTH_KEY: self.max_read_length,
        }
        return summary
    
class slow_metric_finder(base):

    """
    Find "slow" metric types which should be evaluated after downsampling (if any)

    Includes:
    - bedtools
    - custom iteration over BAM reads, eg. to process CIGAR strings
    """

    READ_1_INDEX = 0
    READ_2_INDEX = 1
    READ_UNKNOWN_INDEX = 2
    READ_NAMES =  ['1', '2', '?']

    # Relevant CIGAR operations
    CIGAR_OP_NAMES = {
        0: 'aligned',
        1: 'insertion',
        2: 'deletion',
        4: 'soft clip',
        5: 'hard clip',
    }
    # CIGAR op indices which increment the query cycle
    CONSUMES_QUERY = set([0,1,4,7,8])

    # metric field names
    HARD_CLIP_KEY = 'hard clip bases'
    SOFT_CLIP_KEY = 'soft clip bases'
    MISSING_MD_KEY = 'readsMissingMDtags'

    # fields for unknown read stats
    READ_COUNT_KEY = 'read count'
    LENGTH_TOTAL_KEY = 'length total'
    LENGTH_HISTOGRAM_KEY = 'length histogram'
    QUALITY_BY_CYCLE_KEY = 'quality by cycle'
    TOTAL_BY_CYCLE_KEY = 'total by cycle'
    QUALITY_HISTOGRAM_KEY = 'quality histogram'

    def __init__(self, bam_path, target_path, read_lengths):
        self.bam_path = bam_path
        self.target_path = target_path
        self.read_lengths = read_lengths
        self.metrics = self.evaluate_all_metrics()

    def evaluate_all_metrics(self):
        metrics = {}
        metric_subsets = [
            self.evaluate_bedtools_metrics(),
            self.evaluate_custom_metrics()
        ]
        for metric_subset in metric_subsets:
            for key in metric_subset.keys():
                metrics[key] = metric_subset[key]
        return metrics

    def evaluate_bedtools_metrics(self):
        metrics = {}
        bamBedTool = pybedtools.BedTool(self.bam_path)
        targetBedTool = pybedtools.BedTool(self.target_path)
        metrics['number of targets'] = targetBedTool.count()
        metrics['total target size'] = sum(len(f) for f in targetBedTool.features())
        metrics['reads on target'] = len(bamBedTool.intersect(self.target_path))
        # TODO add bedtools coverage metrics?
        #coverage = targetBedTool.coverage(self.bam_path)
        #print(coverage)
        return metrics

    def evaluate_custom_metrics(self):
        [metrics, ur_stats, start_points] = self.process_bam_file()
        metrics['reads per start point'] = self.evaluate_reads_per_start_point(start_points)
        ur_count = ur_stats[self.READ_COUNT_KEY]
        if ur_count > 0:
            total = float(ur_stats[self.LENGTH_TOTAL_KEY])
            metrics['read ? average length'] = round(total/ur_count, self.PRECISION)
        else:
            metrics['read ? average length'] = None
        metrics['read ? length histogram'] = ur_stats[self.LENGTH_HISTOGRAM_KEY]
        ur_quality_by_cycle = ur_stats[self.QUALITY_BY_CYCLE_KEY]
        ur_total_by_cycle = ur_stats[self.TOTAL_BY_CYCLE_KEY]
        for cycle in ur_quality_by_cycle.keys():
            quality = ur_quality_by_cycle[cycle]
            total = ur_total_by_cycle[cycle]
            q = round(float(quality)/total, self.PRECISION) if total > 0 else 0
            ur_quality_by_cycle[cycle] = q
        metrics['read ? quality by cycle'] = ur_quality_by_cycle
        metrics['read ? quality histogram'] = ur_stats[self.QUALITY_HISTOGRAM_KEY]
        return metrics

    def evaluate_reads_per_start_point(self, start_points):
        """
        Find reads per start point, defined as (unique reads)/(unique start points)
        """
        reads_per_sp = None
        # count the reads, excluding secondary alignments (but including unmapped reads)
        unique_reads = int(pysam.view('-c', '-F', '256', self.bam_path).strip())
        if start_points > 0:
            reads_per_sp = round(float(unique_reads)/start_points, self.FINE_PRECISION)
        else:
            reads_per_sp = 0.0
        return reads_per_sp

    def initialize_bam_metrics(self):
        metrics = {
            self.HARD_CLIP_KEY: 0,
            self.SOFT_CLIP_KEY: 0,
            self.MISSING_MD_KEY: 0,
        }
        read_length_keys = [
            self.READ_1_LENGTH_KEY,
            self.READ_2_LENGTH_KEY,
            self.MAX_READ_LENGTH_KEY,
        ]
        for op_name in self.CIGAR_OP_NAMES.values():
            for i in range(len(self.READ_NAMES)):
                key = 'read %s %s by cycle' % (self.READ_NAMES[i], op_name)
                read_length = self.read_lengths[read_length_keys[i]]
                if read_length == 0:
                    metrics[key] = {} # placeholder
                else:
                    metrics[key] = {j:0 for j in range(1, read_length+1) }
        return metrics

    def initialize_unknown_read_stats(self):
        max_read_length = self.read_lengths[self.MAX_READ_LENGTH_KEY]
        stats = {
            self.READ_COUNT_KEY: 0,
            self.LENGTH_TOTAL_KEY: 0,
            self.LENGTH_HISTOGRAM_KEY: {},
            self.QUALITY_BY_CYCLE_KEY: {i : 0 for i in range(1, max_read_length+1)},
            self.TOTAL_BY_CYCLE_KEY: {i : 0 for i in range(1, max_read_length+1)},
            self.QUALITY_HISTOGRAM_KEY: {}
        }
        return stats

    def process_bam_file(self):
        """
        Iterate through the BAM file and process:
        - CIGAR strings
        - Missing MD tags
        - Hard clip counts
        - Unknown read stats (corresponding stats for R1/R2 are in fast_metric_finder)

        Data will be further processed to get metrics for export in JSON
        """
        metrics = self.initialize_bam_metrics()
        ur_stats = self.initialize_unknown_read_stats()
        start_point_set = set()
        for read in pysam.AlignmentFile(self.bam_path, 'rb').fetch(until_eof=True):
            read_index = None
            if read.is_read1: read_index = self.READ_1_INDEX
            elif read.is_read2: read_index = self.READ_2_INDEX
            else: read_index = self.READ_UNKNOWN_INDEX
            start_point_set = self.update_start_point_set(start_point_set, read)
            metrics = self.update_metrics(metrics, read, read_index)
            if read_index == self.READ_UNKNOWN_INDEX:
                ur_stats[self.READ_COUNT_KEY] += 1
                ur_len = read.query_length
                ur_stats[self.LENGTH_TOTAL_KEY] += ur_len
                if ur_len in ur_stats[self.LENGTH_HISTOGRAM_KEY]:
                    ur_stats[self.LENGTH_HISTOGRAM_KEY][ur_len] += 1
                else:
                    ur_stats[self.LENGTH_HISTOGRAM_KEY][ur_len] += 1
                for i in range(read.query_length):
                    q = read.query_qualities[i]
                    ur_stats[self.QUALITY_BY_CYCLE_KEY][i+1] += q
                    ur_stats[self.TOTAL_BY_CYCLE_KEY][i+1] += 1
                    if q in ur_stats[self.QUALITY_HISTOGRAM_KEY]:
                        ur_stats[self.QUALITY_HISTOGRAM_KEY][q] += 1
                    else:
                        ur_stats[self.QUALITY_HISTOGRAM_KEY][q] = 1
        start_points = len(start_point_set)
        return (metrics, ur_stats, start_points)

    def update_metrics(self, metrics, read, read_index):
        """ Update the metrics dictionary for a pysam 'read' object """
        if not read.has_tag('MD'):
            metrics[self.MISSING_MD_KEY] += 1
        if read.query_length == 0: # all bases are hard clipped
            metrics[self.HARD_CLIP_KEY] += read.infer_read_length()
            return metrics
        cycle = 1
        if read.cigartuples != None:
            if read.is_reverse: cigar_list = reversed(read.cigartuples)
            else: cigar_list = read.cigartuples
            for (op, length) in cigar_list:
                if op in self.CIGAR_OP_NAMES:
                    if op == 4: metrics['soft clip bases'] += length
                    elif op == 5: metrics['hard clip bases'] += length
                    for i in range(length):
                        key = 'read %s %s by cycle' % (self.READ_NAMES[read_index],
                                                       self.CIGAR_OP_NAMES[op])
                        metrics[key][cycle] += 1
                        if op in self.CONSUMES_QUERY: cycle += 1
                elif op in self.CONSUMES_QUERY:
                    cycle += length
        return metrics

    def update_start_point_set(self, start_point_set, read):
        """ Update the start point set for a pysam 'read' object """
        if not read.is_unmapped:
            start_point_set.add((read.reference_name, read.reference_start))
        return start_point_set
