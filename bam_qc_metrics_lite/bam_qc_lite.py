#! /usr/bin/env python3

"""Classes to compute BAM QC metrics"""

import csv
import json
import logging
import os
import re
import sys
import tempfile
import pysam

import bam_qc_metrics_lite


class base_constants(object):
    """
    Class for shared constants
    """
    PRECISION = 1  # number of decimal places for rounded output
    FINE_PRECISION = 3  # finer precision, e.g. for reads_per_start_point output
    MAX_REPORTED_COVERAGE = 1000
    ALIGNMENT_REF_KEY = 'alignment reference'
    INSERT_MAX_KEY = 'insert max'
    READ_1_LENGTH_KEY = 'read 1'
    READ_2_LENGTH_KEY = 'read 2'
    MAX_READ_LENGTH_KEY = 'max_read_length'
    TOTAL_READS_KEY = 'total reads'
    UNMAPPED_READS_KEY = 'unmapped reads'
    PACKAGE_VERSION_KEY = 'package version'
    # shared keys for config dictionary
    CONFIG_KEY_DEBUG = 'debug'
    CONFIG_KEY_INSERT_MAX = 'insert max'
    CONFIG_KEY_LOG = 'log path'
    CONFIG_KEY_REFERENCE = 'reference'
    CONFIG_KEY_VERBOSE = 'verbose'
    # mean coverage keys
    MEAN_TARGET_COVERAGE = 'mean target coverage'
    MEAN_GENOMIC_COVERAGE = 'mean genomic coverage'


class base(base_constants):
    """
    Class for methods shared between bam_qc and fast_metric_writer
    """

    @staticmethod
    def configure_logger(log_path=None, debug=False, verbose=False):
        logger = logging.getLogger(__name__)
        log_level = logging.WARN
        if debug:
            log_level = logging.DEBUG
        elif verbose:
            log_level = logging.INFO
        logger.setLevel(log_level)
        handler = logging.StreamHandler() if log_path is None else logging.FileHandler(log_path)
        handler.setLevel(log_level)
        formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s: %(message)s',
                                      datefmt='%Y-%m-%d_%H:%M:%S')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        return logger

    @staticmethod
    def validate_key_sets(found, expected):
        """Compare set objects containing config keys; raise an informative error on mismatch"""
        if expected != found:
            not_found_set = expected - found
            not_found = not_found_set if len(not_found_set) > 0 else None
            not_expected_set = found - expected
            not_expected = not_expected_set if len(not_expected_set) > 0 else None
            msg = "Config fields are not valid\n"
            msg = msg + "Fields expected and not found: " + str(not_found) + "\n"
            msg = msg + "Fields found and not expected: " + str(not_expected) + "\n"
            # do not log this message; logger not yet initialized
            raise ValueError(msg)


class validator:
    """Utility functions for validating arguments to command-line scripts"""

    MINIMUM_SAMPLE_LEVEL = 1000

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def validate_positive_integer(arg, name):
        param = None
        valid = True
        try:
            param = int(arg)
        except ValueError:
            sys.stderr.write("ERROR: %s must be an integer.\n" % name)
            valid = False
        if param is not None and param < 0:
            sys.stderr.write("ERROR: %s cannot be negative.\n" % name)
            valid = False
        return valid


class version_updater(base_constants):
    FILENAMES = ['expected.json',
                 'expected_downsampled.json',
                 'expected_fast_metrics.json',
                 'expected_no_target.json',
                 'expected_downsampled_rs88.json']

    def __init__(self):
        self.package_version = bam_qc_metrics_lite.read_package_version()

    def get_filenames(self):
        return self.FILENAMES

    def update_files(self, data_dir, version):
        for name in self.FILENAMES:
            json_path = os.path.join(data_dir, name)
            with open(json_path) as f:
                data = json.loads(f.read())
            data[self.PACKAGE_VERSION_KEY] = version
            out = open(json_path, 'w')
            out.write(json.dumps(data, sort_keys=True, indent=4))
            out.close()


class bam_qc_lite(base):
    CONFIG_KEY_READS_ON_TARGET = 'reads on target'
    CONFIG_KEY_BAM = 'bam file'
    CONFIG_KEY_TARGET_PATH = 'target path'
    CONFIG_KEY_SAMSTATS = 'samstats'
    CONFIG_KEY_COVERAGE = 'coverage histogram'
    CONFIG_KEY_METADATA = 'metadata'
    CONFIG_KEY_MARK_DUPLICATES = 'mark duplicates'
    CONFIG_KEY_TARGET_COVERAGE = 'targeted coverage'
    CONFIG_KEY_UNIQUE_READS = 'unique reads'
    CONFIG_KEY_CIGAR_METRICS = 'cigar metrics'
    CONFIG_KEY_DOWNSAMPLED_READS = 'downsampled reads'
    CONFIG_KEY_RANDOM_SEED = 'random seed'
    CONFIG_KEY_SAMPLE = 'sample'
    CONFIG_KEY_TARGET = 'target'
    CONFIG_KEY_TEMP_DIR = 'temp dir'
    CONFIG_KEY_WORKFLOW_VERSION = 'workflow version'
    DEFAULT_MARK_DUPLICATES_METRICS = {
        "LIBRARY": None,
        "PERCENT_DUPLICATION": None,
        "READ_PAIRS_EXAMINED": None,
        "READ_PAIR_DUPLICATES": None,
        "UNPAIRED_READS_EXAMINED": None,
        "UNPAIRED_READ_DUPLICATES": None
    }
    DEFAULT_TARGET_COVERAGE_METRICS = {
        "bases per target": None,
        "coverage per target": None,
        "coverage histogram": None,
        "number of targets": None,
        "reads on target": None,
        "target sizes": None,
        "total bases on target": None,
        "total target size": None
    }
    METADATA_KEYS = [
        # for un-merged BAM input
        'barcode',
        'instrument',
        'lane',
        'library',
        'run name',
        'sample',
        # for merged BAM input
        'donor',
        'group id',
        'library design',
        'tissue origin',
        'tissue type'
        # input/filtered read counts (for merged and un-merged input)
        # 'meta' distinguishes from similarly named fields generated in this package, eg:
        # - 'unmapped reads' from the in-package quality filter
        # - 'total reads' from samtools stats output
    ]
    DEFAULT_RANDOM_SEED = 42

    def __init__(self, config):
        self.validate_config_fields(config)
        # read instance variables from config
        self.logger = self.configure_logger(
            config[self.CONFIG_KEY_LOG],
            config[self.CONFIG_KEY_DEBUG],
            config[self.CONFIG_KEY_VERBOSE]
        )
        self.expected_insert_max = config[self.CONFIG_KEY_INSERT_MAX]
        '''
           Next few calls for reading metrics from external tools, 
           - samblaster
           - mosdepth
           - samtools
        '''
        self.bam_file = config[self.CONFIG_KEY_BAM]
        self.mark_duplicates_metrics = self.read_mark_dup(config[self.CONFIG_KEY_MARK_DUPLICATES])
        ''' If we have a bam file, retrieve CIGAR metrics'''
        self.coverage_metrics = self.read_coverage_histogram(config[self.CONFIG_KEY_COVERAGE])
        self.target_coverage_metrics = None
        #if config[self.CONFIG_KEY_TARGET_PATH] is not None:
        #    if validator.validate_input_file(config[self.CONFIG_KEY_TARGET_PATH]):
        self.target_coverage_metrics = self.read_targeted_coverage_summary(config[self.CONFIG_KEY_TARGET_COVERAGE])
        self.target_path = config[self.CONFIG_KEY_TARGET_PATH]
        self.downsample_to = config[self.CONFIG_KEY_DOWNSAMPLED_READS]
        '''Fetching external metrics ends here'''
        self.metadata = self.read_metadata(config[self.CONFIG_KEY_METADATA])
        seed = config[self.CONFIG_KEY_RANDOM_SEED]
        self.random_seed = seed if seed is not None else self.DEFAULT_RANDOM_SEED
        self.reference = config[self.CONFIG_KEY_REFERENCE]
        self.targeted_coverage = config[self.CONFIG_KEY_TARGET_COVERAGE]
        self.workflow_version = config[self.CONFIG_KEY_WORKFLOW_VERSION]
        self.reads_on_target = config[self.CONFIG_KEY_READS_ON_TARGET]
        u_reads = config[self.CONFIG_KEY_UNIQUE_READS]
        self.unique_reads = u_reads if u_reads is not None else None
        # define other instance variables
        self.fast_metrics = None
        self.package_version = bam_qc_metrics_lite.read_package_version()
        self.qual_fail_reads = None
        self.slow_metrics = None
        (self.tmpdir, self.tmp_object) = self.setup_tmpdir(config[self.CONFIG_KEY_TEMP_DIR])
        self.unmapped_excluded_reads = None
        # find metrics; if an error occurs, do logging and cleanup before exit
        try:
            self._find_metrics(config[self.CONFIG_KEY_SAMSTATS])
        except Exception as e:
            self.logger.exception("Unexpected error: {0}".format(e))
            self.cleanup()
            raise

    def _find_metrics(self, samstats_input_path):
        self.logger.info("Started bam_qc_lite processing")
        self.logger.info("Started computing fast bam_qc_lite metrics")
        fast_finder = fast_metric_finder(samstats_input_path,
                                         self.reference,
                                         self.expected_insert_max,
                                         self.logger)
        self.fast_metrics = fast_finder.get_metrics()
        if self.bam_file is not None and validator.validate_input_file(self.bam_file):
            slow_finder = slow_metric_finder(self.bam_file, self.downsample_to, fast_finder.read_length_summary(), self.unique_reads, self.logger)
            self.slow_metrics = slow_finder.metrics
        self.logger.info("Finished computing fast bam_qc_lite metrics")

    def cleanup(self):
        """
        Cleanup temporary files and directory created by bam-qc-metrics

        Temporary directory object (if any) will be automatically deleted when bam_qc object
        is out of scope. This method allows explicit cleanup for:
           - Avoiding warnings in the test harness
           - Ensuring temporary file cleanup on error
        """
        self.tmp_object.cleanup()
        self.logger.info("Cleaned up temporary directory %s" % self.tmpdir)

    def read_metadata(self, metadata_path):
        if metadata_path is not None:
            with open(metadata_path) as f:
                raw_metadata = json.loads(f.read())
            metadata = {key: raw_metadata.get(key) for key in self.METADATA_KEYS}
        else:
            self.logger.info("Metadata file not given, using empty defaults")
            metadata = {key: None for key in self.METADATA_KEYS}
        return metadata

    def read_mark_dup(self, input_path):
        # scrape metrics file
        if input_path is None:
            return self.DEFAULT_MARK_DUPLICATES_METRICS
        reading_in = False
        with open(input_path) as f:
            lines = f.readlines()
        metrics = {}
        for line in lines:
            line = line.strip()
            if re.match('samblaster: -*$', line):
                # DuplicationMetrics section start
                reading_in = True
                continue
            if not reading_in:
                continue

            temp_values = line.split()
            if len(temp_values) > 5:
                if temp_values[1] == 'Orphan/Singleton':
                    metrics['UNPAIRED_READS_EXAMINED'] = int(temp_values[2])
                    metrics['UNPAIRED_READ_DUPLICATES'] = int(temp_values[4])
                elif temp_values[1] == 'Total':
                    metrics['READ_PAIRS_EXAMINED'] = int(temp_values[2])
                    metrics['READ_PAIR_DUPLICATES'] = int(temp_values[4])
                    metrics['PERCENT_DUPLICATION'] = float(temp_values[5])

        return metrics

    def read_coverage_histogram(self, c_histo_path):
        """ Do not report coverage higher than MAX_REPORTED_COVERAGE """
        try:
            if validator.validate_input_file(c_histo_path):
                metrics = {}
                with open(c_histo_path, 'r') as ch:
                    my_dict = json.load(ch)
                    sorted_dict = dict(sorted(my_dict.items(), key=lambda item: int(item[0])))
                    # Exclude keys greater than 10 and values equal 0:
                    metrics = dict(filter(lambda item: int(item[0]) <= self.MAX_REPORTED_COVERAGE and item[1] != 0,
                                          sorted_dict.items()))
                return metrics
            return None
        except Exception as e:
            self.logger.exception("Error reading coverage data: {0}".format(e))
            return None

    def read_targeted_coverage_summary(self, c_summary_path):
        """ scrape metrics file """
        number_of_targets_key = 'number of targets'
        total_target_size_key = 'total target size'
        coverage_per_target_key = 'coverage per target'
        bases_per_target_key = 'bases per target'
        total_bases_on_target_key = 'total bases on target'
        target_sizes_key = 'target sizes'
        if c_summary_path is None:
            return self.DEFAULT_TARGET_COVERAGE_METRICS
        with open(c_summary_path) as f:
            lines = f.readlines()

        region_lines = []
        genome_lines = []

        lines.pop(0)  # remove the header line
        for line in lines:
            line.strip()
            fields = line.split("\t")
            if re.search(r'_region$', fields[0]):
                region_lines.append(line)
            else:
                genome_lines.append(line)

        metrics = {}
        data = self.evaluate_coverage(region_lines) if len(region_lines) > 0 else self.evaluate_coverage(genome_lines)

        if data is not None and isinstance(data, dict):
            metrics[number_of_targets_key] = len(data[target_sizes_key])
            metrics[total_target_size_key] = sum(data[target_sizes_key][f] for f in data[target_sizes_key].keys())
            metrics[total_bases_on_target_key] = sum(data[bases_per_target_key][f] for f in data[bases_per_target_key].keys())

            metrics[coverage_per_target_key] = data[coverage_per_target_key]
            metrics[bases_per_target_key] = data[bases_per_target_key]
            metrics[target_sizes_key] = data[target_sizes_key]

            if self.MEAN_TARGET_COVERAGE in data.keys():
                metrics[self.MEAN_TARGET_COVERAGE] = data[self.MEAN_TARGET_COVERAGE]
            if self.MEAN_GENOMIC_COVERAGE in data.keys():
                metrics[self.MEAN_GENOMIC_COVERAGE] = data[self.MEAN_GENOMIC_COVERAGE]

        self.logger.info("Processed mosdepth metrics")
        return metrics

    def evaluate_coverage(self, summary_lines: list):
        if summary_lines is None or len(summary_lines) == 0:
            return None
        """
           We expect the following fields:
           chrom - target
           length - bases per target, target size
           bases - cumulative bases covered, 
           mean - use for coverage per target
           min - not used
           max - not used
        """
        results = {'bases per target': {},
                   'coverage per target': {},
                   'target sizes': {}}
        coverage_buffer = 0
        total_length = 0
        coverage_type = self.MEAN_GENOMIC_COVERAGE
        for line in summary_lines:
            line = line.strip()
            fields = line.split("\t")
            if fields[0].startswith("total"):
                continue
            try:
                target = fields[0]
                target = target.rstrip('_region')
                coverage_type = self.MEAN_TARGET_COVERAGE if fields[0] != target else self.MEAN_GENOMIC_COVERAGE
                [length, bases] = [int(x) for x in fields[1:3]]
                mean_coverage = float(fields[3])
                total_length += length
                coverage_buffer += mean_coverage*length
                bases_covered = int(bases/mean_coverage) if mean_coverage > 0 else 0
                results['bases per target'][target] = bases_covered
                results['coverage per target'][target] = round(mean_coverage, self.PRECISION)
                results['target sizes'][target] = length
            except ValueError:
                msg = "Cannot parse coverage in mosdepth output: '" + str(line) + "'"
                self.logger.error(msg)
                raise
        results[coverage_type] = round(coverage_buffer/total_length, self.PRECISION) if total_length > 0 else 0
        return results

    @staticmethod
    def setup_tmpdir(tmpdir_base=None):
        tmp_object = tempfile.TemporaryDirectory(prefix='bam_qc_', dir=tmpdir_base)
        tmpdir = tmp_object.name
        return tmpdir, tmp_object

    def validate_config_fields(self, config):
        """ Validate keys of the config dictionary for __init__ """
        expected = {self.CONFIG_KEY_BAM,
                    self.CONFIG_KEY_SAMSTATS,
                    self.CONFIG_KEY_COVERAGE,
                    self.CONFIG_KEY_METADATA,
                    self.CONFIG_KEY_MARK_DUPLICATES,
                    self.CONFIG_KEY_DOWNSAMPLED_READS,
                    self.CONFIG_KEY_TARGET_COVERAGE,
                    self.CONFIG_KEY_READS_ON_TARGET,
                    self.CONFIG_KEY_UNIQUE_READS,
                    self.CONFIG_KEY_CIGAR_METRICS,
                    self.CONFIG_KEY_TARGET_PATH,
                    self.CONFIG_KEY_REFERENCE,
                    self.CONFIG_KEY_DEBUG,
                    self.CONFIG_KEY_INSERT_MAX,
                    self.CONFIG_KEY_LOG,
                    self.CONFIG_KEY_RANDOM_SEED,
                    self.CONFIG_KEY_TEMP_DIR,
                    self.CONFIG_KEY_VERBOSE,
                    self.CONFIG_KEY_WORKFLOW_VERSION}
        found = set(config.keys())
        self.validate_key_sets(found, expected)

    def write_output(self, out_path):
        output = {}
        for key in self.METADATA_KEYS:
            output[key] = self.metadata.get(key)
        for key in self.fast_metrics.keys():
            output[key] = self.fast_metrics.get(key)
        if self.slow_metrics is not None and isinstance(self.slow_metrics, dict):
            for key in self.slow_metrics.keys():
                output[key] = self.slow_metrics.get(key)
        output[self.ALIGNMENT_REF_KEY] = self.reference
        output[self.INSERT_MAX_KEY] = self.expected_insert_max
        if self.mark_duplicates_metrics is not None and isinstance(self.mark_duplicates_metrics, dict):
            output['mark duplicates'] = self.mark_duplicates_metrics
        if self.coverage_metrics is not None and isinstance(self.coverage_metrics, dict):
            output['coverage_histogram'] = self.coverage_metrics
        if self.target_coverage_metrics is not None and isinstance(self.target_coverage_metrics, dict):
            output.update(self.target_coverage_metrics)
        if self.reads_on_target is not None:
            output['reads on target'] = int(self.reads_on_target)
        output['qual fail reads'] = self.qual_fail_reads
        output[self.PACKAGE_VERSION_KEY] = self.package_version
        output['target file'] = os.path.split(self.target_path)[-1] if self.target_path else None
        output['workflow version'] = self.workflow_version
        if out_path != '-':
            out_file = open(out_path, 'w')
            dest = out_path
        else:
            out_file = sys.stdout
            dest = 'STDOUT'
        print(json.dumps(output, sort_keys=True), file=out_file)
        if out_path != '-':
            out_file.close()
        self.logger.debug("Wrote JSON output to %s" % dest)


class fast_metric_finder(base_constants):
    """
       uses native samtools stats, rather than custom metrics
    """

    # summary numbers (SN) fields denoted in float_keys are floats; integers otherwise
    FLOAT_KEYS = {'error rate',
                  'average quality',
                  'insert size average',
                  'insert size standard deviation',
                  'percentage of properly paired reads (%)'}

    def __init__(self, samstats_path, reference, expected_insert_max, logger):
        self.samstats_path = samstats_path
        self.reference = reference
        self.expected_insert_max = expected_insert_max
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
            'raw total sequences': self.TOTAL_READS_KEY,
            'reads unmapped': self.UNMAPPED_READS_KEY,
            'non-primary alignments': 'non primary reads',
        }
        self.samtools_stats = self.load_samtools_stats(self.samstats_path)
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
            self.read_length_and_quality_metrics(),
            {"mismatch by cycle": self.parse_mismatch()},
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
            filter(lambda line: line != "" and line[0] != '#', self.samtools_stats),
            delimiter="\t"
        )
        max_read_length = 0
        for row in reader:
            if row[0] == 'SN' and row[1] == 'maximum length:':
                max_read_length = int(row[2])
                break
        return max_read_length

    def get_metrics(self):
        return self.metrics

    def fq_stats(self, rows):
        """
        Compute quality metrics from either FFQ or LFQ entries in samtools stats:
            - Mean quality by cycle
            - Quality histogram
        """
        mean_by_cyc = {}
        histogram = {}
        if len(rows) > 0:
            max_width = max([len(row) for row in rows])
            histogram = {q: 0 for q in range(max_width - 2)}
            for row in rows:
                cycle = int(row[1])
                counts = [int(n) for n in row[2:]]
                total = 0
                count = 0
                for qscore in range(len(counts)):
                    total += counts[qscore] * qscore
                    count += counts[qscore]
                    histogram[qscore] += counts[qscore]
                mean_by_cyc[cycle] = round(float(total) / count, self.PRECISION) if count > 0 else 0
        return mean_by_cyc, histogram

    def mean_read_length(self, rows):
        """
           Process RL (read length), FRL (first RL) or LRL (last RL) rows for mean read length
        """
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
        metrics = {'inserted bases': 0, 'deleted bases': 0, 'insert size histogram': {}}
        reader = csv.reader(
            filter(lambda line: line != "" and line[0] != '#', self.samtools_stats),
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

    def parse_mismatch(self):
        """
           Input is the string returned by 'samtools stats'
           Parse the mismatches by cycle; return an empty dictionary if no data found
        """
        if self.samtools_stats is None:
            self.logger.debug("None input to mismatch parser; returning empty dictionary")
            return {}
        reader = csv.reader(
            filter(lambda line: line != "" and line[0] != '#', self.samtools_stats),
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
                start = 3
                for m in row[start:]:
                    mismatches += int(m)
                mismatch_by_cycle[cycle] = mismatches
            elif row[0] == 'SN' and row[1] == 'sequences:' and int(row[2]) == 0:
                empty = True
        if empty:
            # special case -- no reads found (e.g. no unknown reads)
            # MPC returns 0 mismatches for cycle 1, but we want an empty dictionary
            mismatch_by_cycle = {}
        return mismatch_by_cycle

    @staticmethod
    def read_length_histogram(rows):
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
            filter(lambda line: line != "" and line[0] != '#', self.samtools_stats),
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

    @staticmethod
    def load_samtools_stats(stats_path):
        """
        Re-written function, does not run samtools stats but load results of external samtools stats
        process.
        """
        stats_str = None
        with open(stats_path, 'r') as stats:
            stats_str = stats.readlines()
        return stats_str


class slow_metric_finder(base_constants):
    """
    Find "slow" metric types which should be evaluated if requested

    Includes:
    - custom iteration over BAM reads, e.g. to process CIGAR strings
    """

    READ_1_INDEX = 0
    READ_2_INDEX = 1
    READ_UNKNOWN_INDEX = 2
    READ_NAMES = ['1', '2', '?']

    # Relevant CIGAR operations
    CIGAR_OP_NAMES = {
        0: 'aligned',
        1: 'insertion',
        2: 'deletion',
        4: 'soft clip',
        5: 'hard clip',
    }
    # CIGAR op indices which increment the query cycle
    CONSUMES_QUERY = {0, 1, 4, 7, 8}

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

    def __init__(self, bam_path, downsample_to, read_lengths, unique_reads, logger):
        self.bam_path = bam_path
        self.read_lengths = read_lengths
        self.logger = logger
        self.downsample_to = downsample_to
        self.unique_reads = unique_reads
        self.metrics = self.evaluate_all_metrics()

    def evaluate_all_metrics(self):
        metrics = {}
        metric_subsets = [
            self.evaluate_custom_metrics()
        ]
        for metric_subset in metric_subsets:
            for key in metric_subset.keys():
                metrics[key] = metric_subset[key]
        return metrics

    def evaluate_custom_metrics(self):
        [metrics, ur_stats, start_points] = self.process_bam_file(self.downsample_to)
        # Note that the below instruction needs to run on unique reads for downsampled bam
        if self.unique_reads is not None and start_points > 0:
            metrics['reads per start point'] = round(self.unique_reads/start_points, self.PRECISION)
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
            q = round(float(quality) / total, self.PRECISION) if total > 0 else 0
            ur_quality_by_cycle[cycle] = q
        metrics['read ? quality by cycle'] = ur_quality_by_cycle
        metrics['read ? quality histogram'] = ur_stats[self.QUALITY_HISTOGRAM_KEY]
        self.logger.debug("Found custom metrics from CIGAR strings etc.")
        return metrics

    def process_bam_file(self, downsample_to):
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
        for read in pysam.AlignmentFile(self.bam_path, 'rb').head(downsample_to):
            if read.is_read1:
                read_index = self.READ_1_INDEX
            elif read.is_read2:
                read_index = self.READ_2_INDEX
            else:
                read_index = self.READ_UNKNOWN_INDEX
            start_point_set = self.update_start_point_set(start_point_set, read)
            metrics = self.update_metrics(metrics, read, read_index)
            if read_index == self.READ_UNKNOWN_INDEX:
                ur_stats = self.update_ur_stats(ur_stats, read)
        start_points = len(start_point_set)
        metrics["downsampled total"] = downsample_to
        return metrics, ur_stats, start_points

    def initialize_unknown_read_stats(self):
        max_read_length = self.read_lengths[self.MAX_READ_LENGTH_KEY]
        stats = {
            self.READ_COUNT_KEY: 0,
            self.LENGTH_TOTAL_KEY: 0,
            self.LENGTH_HISTOGRAM_KEY: {},
            self.QUALITY_BY_CYCLE_KEY: {i: 0 for i in range(1, max_read_length + 1)},
            self.TOTAL_BY_CYCLE_KEY: {i: 0 for i in range(1, max_read_length + 1)},
            self.QUALITY_HISTOGRAM_KEY: {}
        }
        return stats

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
                    metrics[key] = {}  # placeholder
                else:
                    metrics[key] = {j: 0 for j in range(1, read_length + 1)}
        return metrics

    def update_metrics(self, metrics, read, read_index):
        """ Update the metrics dictionary for a pysam 'read' object """
        if not read.has_tag('MD'):
            metrics[self.MISSING_MD_KEY] += 1
        if read.query_length == 0:  # all bases are hard clipped
            metrics[self.HARD_CLIP_KEY] += read.infer_read_length()
            return metrics
        cycle = 1
        if read.cigartuples != None:
            if read.is_reverse:
                cigar_list = reversed(read.cigartuples)
            else:
                cigar_list = read.cigartuples
            for (op, length) in cigar_list:
                if op in self.CIGAR_OP_NAMES:
                    if op == 4:
                        metrics['soft clip bases'] += length
                    elif op == 5:
                        metrics['hard clip bases'] += length
                    for _ in range(length):
                        key = 'read %s %s by cycle' % (self.READ_NAMES[read_index],
                                                       self.CIGAR_OP_NAMES[op])
                        metrics[key][cycle] += 1
                        if op in self.CONSUMES_QUERY:
                            cycle += 1
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
        lhk = self.LENGTH_HISTOGRAM_KEY
        qhk = self.QUALITY_HISTOGRAM_KEY
        ur_stats[self.READ_COUNT_KEY] += 1
        ur_len = read.query_length
        ur_stats[self.LENGTH_TOTAL_KEY] += ur_len
        ur_stats[lhk][ur_len] = ur_stats[lhk].get(ur_len, 0) + 1
        for i in range(read.query_length):
            q = read.query_qualities[i]
            ur_stats[self.QUALITY_BY_CYCLE_KEY][i + 1] += q
            ur_stats[self.TOTAL_BY_CYCLE_KEY][i + 1] += 1
            ur_stats[qhk][q] = ur_stats[qhk].get(q, 0) + 1
        return ur_stats
