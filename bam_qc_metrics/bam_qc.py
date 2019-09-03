#! /usr/bin/env python3

"""Main class to compute BAM QC metrics"""

import bam_qc_metrics, csv, json, logging, math, os, re, pybedtools, pysam, sys, tempfile

class base(object):

    """
    Class for shared constants
    """

    PRECISION = 1 # number of decimal places for rounded output
    FINE_PRECISION = 3 # finer precision, eg. for reads_per_start_point output
    READ_1_LENGTH_KEY = 'read 1'
    READ_2_LENGTH_KEY = 'read 2'
    MAX_READ_LENGTH_KEY = 'max_read_length'
    UNMAPPED_READS_KEY = 'unmapped reads'
    PACKAGE_VERSION_KEY = 'package version'


class version_updater(base):

    def __init__(self, data_dir=None, version=None):
        if data_dir == None:
            self.data_dir = bam_qc_metrics.get_data_dir_path()
        else:
            self.data_dir = data_dir
        if version==None:
            self.package_version = bam_qc_metrics.read_package_version()
        else:
            self.package_version = version

    def update_files(self):
        filenames = ['expected.json', 'expected_downsampled.json', 'expected_no_target.json']
        for name in filenames:
            json_path = os.path.join(self.data_dir, name)
            with open(json_path) as f:
                data = json.loads(f.read())
            data[self.PACKAGE_VERSION_KEY] = self.package_version
            out = open(json_path, 'w')
            out.write(json.dumps(data, sort_keys=True, indent=4))
            out.close()


class bam_qc(base):

    CONFIG_KEY_BAM = 'bam'
    CONFIG_KEY_DEBUG = 'debug'
    CONFIG_KEY_INSERT_MAX = 'insert max'
    CONFIG_KEY_LOG = 'log path'
    CONFIG_KEY_MARK_DUPLICATES = 'mark duplicates'
    CONFIG_KEY_METADATA = 'metadata'
    CONFIG_KEY_N_AS_MISMATCH = 'n as mismatch'
    CONFIG_KEY_RANDOM_SEED = 'random_seed'
    CONFIG_KEY_REFERENCE = 'reference'
    CONFIG_KEY_SAMPLE_RATE = 'sample rate'
    CONFIG_KEY_SKIP_BELOW_MAPQ = 'skip below mapq'
    CONFIG_KEY_TARGET = 'target'
    CONFIG_KEY_TEMP_DIR = 'temp dir'
    CONFIG_KEY_VERBOSE = 'verbose'
    CONFIG_KEY_WORKFLOW_VERSION = 'workflow version'
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
    DEFAULT_RANDOM_SEED = 42

    def __init__(self, config):
        self.validate_config_fields(config)
        # read instance variables from config
        self.debug = config[self.CONFIG_KEY_DEBUG] # enable logging first in case of warnings
        self.verbose = config[self.CONFIG_KEY_VERBOSE]
        self.logger = self.configure_logger(config[self.CONFIG_KEY_LOG])
        self.expected_insert_max = config[self.CONFIG_KEY_INSERT_MAX]
        self.mark_duplicates_metrics = self.read_mark_dup(config[self.CONFIG_KEY_MARK_DUPLICATES])
        self.metadata = self.read_metadata(config[self.CONFIG_KEY_METADATA])
        self.n_as_mismatch = config[self.CONFIG_KEY_N_AS_MISMATCH]
        seed = config[self.CONFIG_KEY_RANDOM_SEED]
        self.random_seed = seed if seed != None else self.DEFAULT_RANDOM_SEED
        self.reference = config[self.CONFIG_KEY_REFERENCE]
        self.sample_rate = config[self.CONFIG_KEY_SAMPLE_RATE]
        self.skip_below_mapq = config[self.CONFIG_KEY_SKIP_BELOW_MAPQ]
        self.target_path = config[self.CONFIG_KEY_TARGET]
        self.workflow_version = config[self.CONFIG_KEY_WORKFLOW_VERSION]
        # define other instance variables
        self.fast_metrics = None
        self.package_version = bam_qc_metrics.read_package_version()
        self.qual_fail_reads = None
        self.slow_metrics = None
        (self.tmpdir, self.tmp_object) = self.setup_tmpdir(config[self.CONFIG_KEY_TEMP_DIR])
        self.unmapped_excluded_reads = None
        # apply quality filter (if any); get input path for downsampling
        self.logger.info("Started bam_qc processing")
        result = self.apply_mapq_filter(config[self.CONFIG_KEY_BAM])
        (fast_finder_input_path, self.qual_fail_reads, self.unmapped_excluded_reads) = result
        # find 'fast' metrics on full dataset -- after filtering, before downsampling
        self.logger.info("Started computing fast bam_qc metrics")
        fast_finder = fast_metric_finder(fast_finder_input_path,
                                         self.reference,
                                         self.expected_insert_max,
                                         self.n_as_mismatch,
                                         self.logger)
        self.fast_metrics = self.update_unmapped_count(fast_finder.metrics)
        self.logger.info("Finished computing fast bam_qc metrics")
        # apply downsampling (if any)
        if self.sample_rate != None and self.sample_rate > 1:
            slow_finder_input_path = self.generate_downsampled_bam(fast_finder_input_path,
                                                                   self.sample_rate)
        else:
            self.logger.info("Downsampling is not in effect")
            self.sample_rate = 1
            slow_finder_input_path = fast_finder_input_path
        # find 'slow' metrics on (maybe) downsampled dataset
        self.logger.info("Started computing slow bam_qc metrics")
        slow_finder = slow_metric_finder(slow_finder_input_path,
                                         self.target_path,
                                         fast_finder.read_length_summary(),
                                         self.logger)
        self.slow_metrics = slow_finder.metrics
        self.logger.info("Finished computing slow bam_qc metrics")
        self.logger.info("Finished computing all bam_qc metrics")

    def apply_mapq_filter(self, bam_path):
        """
        Write a BAM file filtered by alignment quality
        Also report the number of reads failing quality filters, and number of unmapped reads
        Samtools steps may be slow on large files; debug logging can be used to track progress
        """
        filtered_bam_path = None
        if self.mapq_filter_is_active():
            excluded_by_mapq_path = os.path.join(self.tmpdir, 'excluded.bam')
            included_by_mapq_path = os.path.join(self.tmpdir, 'included.bam')
            self.logger.info("Started filtering input BAM file by mapping quality")
            self.logger.debug("Started writing included/excluded bam files")
            pysam.view(bam_path,
                       '-b',
                       '-q', str(self.skip_below_mapq),
                       '-o', included_by_mapq_path,
                       '-U', excluded_by_mapq_path,
                       catch_stdout=False)
            self.logger.debug("Finished writing included/excluded bam files")
            self.logger.debug("Started counting excluded reads")
            total_failed_reads = int(pysam.view('-c', excluded_by_mapq_path).strip())
            self.logger.debug("Finished counting excluded reads")
            self.logger.info("Total reads excluded by mapq filter: %d", total_failed_reads)
            # unmapped reads will fail the mapping quality filter, by definition
            # so if the quality filter is applied, find unmapped total from the excluded reads
            self.logger.debug("Started counting unmapped reads")
            unmapped_raw = pysam.view('-c', '-f', '4', excluded_by_mapq_path)
            self.logger.debug("Finished counting unmapped reads")
            total_unmapped_and_excluded = int(unmapped_raw.strip())
            self.logger.debug("Total unmapped and excluded reads: %d", total_unmapped_and_excluded)
            filtered_bam_path = included_by_mapq_path
            self.logger.info("Finished filtering input BAM file by mapping quality")
        else:
            self.logger.info("Mapping quality filter is not in effect")
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
            self.logger.info("Cleaned up default temporary directory %s" % self.tmpdir)
        else:
            self.logger.info("Omitting cleanup for user-specified "+\
                             "temporary directory %s" % self.tmpdir)

    def configure_logger(self, log_path=None):
        logger = logging.getLogger(__name__)
        if self.debug:
            log_level = getattr(logging, 'DEBUG')
        elif self.verbose:
            log_level = getattr(logging, 'INFO')
        else:
            log_level = getattr(logging, 'WARN')
        logger.setLevel(log_level)
        handler = None
        if log_path==None:
            handler = logging.StreamHandler()
        else:
            handler = logging.FileHandler(log_path)
        handler.setLevel(log_level)
        formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s: %(message)s',
                                      datefmt='%Y-%m-%d_%H:%M:%S')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        return logger
            
    def generate_downsampled_bam(self, bam_path, sample_rate):
        """
        Write a temporary downsampled BAM file
        """
        downsampled_path = None
        if sample_rate == 1:
            self.logger.info("Sample rate = 1, omitting down sampling")
            downsampled_path = bam_path
        elif sample_rate < 1:
            msg = "Sample rate cannot be less than 1"
            self.logger.error(msg)
            raise ValueError(msg)
        else:
            # We are sampling every Nth read
            # Argument to `samtools -s` is of the form RANDOM_SEED.DECIMAL_RATE
            # Eg. for random seed 42 and sample rate 4, 42 + (1/4) = 42.25
            self.logger.info("Starting downsampling with sample rate %i", sample_rate)
            downsampled_path = os.path.join(self.tmpdir, 'downsampled.bam')
            x = 1.0/sample_rate
            sf = 5 # significant figures for rounding the sample rate
            sample_decimal = round(x, sf - int(math.floor(math.log10(x)))-1)
            sample_arg = str(self.random_seed + sample_decimal)
            self.logger.debug("Sampling argument to 'samtools view' is %s", sample_arg)
            pysam.view('-u', '-s', sample_arg, '-o', downsampled_path, bam_path, catch_stdout=False)
            # sanity check on the downsampled file
            self.logger.info("Finished downsampling with sample rate %i", sample_rate)
            sampled = int(pysam.view('-c', downsampled_path).strip())
            self.logger.info("%i reads in downsampled set", sampled)
            if sampled < self.DOWNSAMPLE_WARNING_THRESHOLD:
                self.logger.warning("WARNING: Only %i reads remain after downsampling", sampled)
        return downsampled_path

    def mapq_filter_is_active(self):
        return self.skip_below_mapq != None and self.skip_below_mapq > 0

    def read_metadata(self, metadata_path):
        metadata = None
        if metadata_path != None:
            with open(metadata_path) as f: raw_metadata = json.loads(f.read())
            metadata = {key: raw_metadata.get(key) for key in self.METADATA_KEYS}
        else:
            self.logger.info("Metadata file not given, using empty defaults")
            metadata = {key: None for key in self.METADATA_KEYS}
        return metadata
    
    def read_mark_dup(self, input_path):
        if input_path == None:
            return self.DEFAULT_MARK_DUPLICATES_METRICS
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
                self.logger.error(msg)
                raise ValueError(msg)
        if len(keys) == len(values) + 1 and keys[-1] == 'ESTIMATED_LIBRARY_SIZE':
            # field is empty (no trailing \t) for low coverage; append a default value
            values.append(0)
        elif len(keys) != len(values):
            # otherwise, mismatched key/value totals are an error
            msg = "Numbers of keys and values in %s do not match" % input_path
            self.logger.error(msg)
            raise ValueError(msg)
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
            tmpdir = tmp_object.name
        else:
            tmp_object = None
        return (tmpdir, tmp_object)

    def update_unmapped_count(self, metrics):
        """
        Input 'metrics' was calculated after mapping quality filter was applied
        If mapping quality filter is in effect, all unmapped reads should have been filtered out
        Check this is so, and update the 'unmapped reads' metric entry
        """
        if self.mapq_filter_is_active():
            if metrics[self.UNMAPPED_READS_KEY] > 0:
                msg = "Mapping quality filter is in effect, so all unmapped reads should have been "+\
                      "removed from BAM input; but 'reads unmapped' field from samtools stats is %d." \
                      % metrics[self.UNMAPPED_READS_KEY]
                self.logger.error(msg)
                raise ValueError(msg)
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
        output[self.PACKAGE_VERSION_KEY] = self.package_version
        output['sample rate'] = self.sample_rate
        output['target file'] = os.path.split(self.target_path)[-1] if self.target_path else None
        output['workflow version'] = self.workflow_version
        if out_path != '-':
            out_file = open(out_path, 'w')
            dest = out_path
        else:
            out_file = sys.stdout
            dest = 'STDOUT'
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()
        self.logger.debug("Wrote JSON output to %s" % dest)

    def validate_config_fields(self, config):
        """ Validate keys of the config dictionary for __init__ """
        expected = set([
            self.CONFIG_KEY_BAM,
            self.CONFIG_KEY_DEBUG,
            self.CONFIG_KEY_INSERT_MAX,
            self.CONFIG_KEY_LOG,
            self.CONFIG_KEY_MARK_DUPLICATES,
            self.CONFIG_KEY_METADATA,
            self.CONFIG_KEY_N_AS_MISMATCH,
            self.CONFIG_KEY_RANDOM_SEED,
            self.CONFIG_KEY_REFERENCE,
            self.CONFIG_KEY_SAMPLE_RATE,
            self.CONFIG_KEY_SKIP_BELOW_MAPQ,
            self.CONFIG_KEY_TARGET,
            self.CONFIG_KEY_TEMP_DIR,
            self.CONFIG_KEY_VERBOSE,
            self.CONFIG_KEY_WORKFLOW_VERSION
        ])
        found = set(config.keys())
        if not expected == found:
            not_found_set = expected - found
            not_found = not_found_set if len(not_found_set) > 0 else None
            not_expected_set = found - expected
            not_expected = not_expected_set if len(not_expected_set) > 0 else None
            msg = "Config fields are not valid\n"
            msg = msg+"Fields expected and not found: "+str(not_found)+"\n"
            msg = msg+"Fields found and not expected: "+str(not_expected)+"\n"
            # do not log this message; logger not yet initialized
            raise ValueError(msg)

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
    
    def __init__(self, bam_path, reference, expected_insert_max, n_as_mismatch, logger):
        self.bam_path = bam_path
        self.reference = reference
        self.expected_insert_max = expected_insert_max
        self.n_as_mismatch = n_as_mismatch
        self.logger = logger
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
            r1_mismatch = self.parse_mismatch(pysam.stats('-r', self.reference,
                                                          '-f', '64',
                                                          self.bam_path))
            r2_mismatch = self.parse_mismatch(pysam.stats('-r', self.reference,
                                                          '-f', '128',
                                                          self.bam_path))
            ur_mismatch = self.parse_mismatch(pysam.stats('-r', self.reference,
                                                          '-F', '192',
                                                          self.bam_path))
        else:
            r1_mismatch = {}
            r2_mismatch = {}
            ur_mismatch = {}
        mismatches['read 1 mismatch by cycle'] = r1_mismatch
        mismatches['read 2 mismatch by cycle'] = r2_mismatch
        mismatches['read ? mismatch by cycle'] = ur_mismatch
        self.logger.debug("Found mismatches by cycle")
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
        self.logger.debug("Found miscellaneous metrics from samtools stats")
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
        self.logger.debug("Found read length and quality metrics")
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

    def __init__(self, bam_path, target_path, read_lengths, logger):
        self.bam_path = bam_path
        self.target_path = target_path
        self.read_lengths = read_lengths
        self.logger = logger
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
        if self.target_path:
            bamBedTool = pybedtools.BedTool(self.bam_path)
            targetBedTool = pybedtools.BedTool(self.target_path)
            metrics['number of targets'] = targetBedTool.count()
            metrics['total target size'] = sum(len(f) for f in targetBedTool.features())
            metrics['reads on target'] = len(bamBedTool.intersect(self.target_path))
            self.logger.debug("Found bedtools metrics")
        else:
            metrics['number of targets'] = None
            metrics['total target size'] = None
            metrics['reads on target'] = None
            self.logger.warning("No target path given, omitting bedtools metrics")
        # placeholders; TODO implement these metrics
        metrics['total coverage'] = None
        metrics['coverage per target'] = None
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
        self.logger.debug("Found custom metrics from CIGAR strings etc.")
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
               ur_stats = self.update_ur_stats(ur_stats, read)
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

    def update_ur_stats(self, ur_stats, read):
        """ Update the unknown-read stats for a pysam 'read' object """
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
        return ur_stats
