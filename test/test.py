#! /usr/bin/env python3

import json, os, re, shutil, subprocess, sys, tempfile, unittest

from bam_qc_metrics import bam_qc, fast_metric_finder, fast_metric_writer, version_updater

class test(unittest.TestCase):

    OICR_REF = '/oicr/data/reference/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19/fasta/hg19.fa'
    LOCAL_REF = os.path.join(os.path.expanduser('~'), 'data/reference/hg19.fa')

    def setUp(self):
        self.quality = 30
        self.insert_max = 1500
        self.sample_level = 10000
        self.sample_default = 1100000
        self.tmp = tempfile.TemporaryDirectory(prefix='bam_qc_test_')
        self.tmpdir = self.tmp.name
        self.testdir = os.path.dirname(os.path.realpath(__file__))
        self.datadir = os.path.realpath(os.path.join(self.testdir, '..', 'data'))
        self.metadata_path = os.path.join(self.datadir, 'metadata.json')
        self.metadata_path_alternate = os.path.join(self.datadir, 'metadata_alternate.json')
        self.bam_path = os.path.join(self.datadir, 'neat_5x_EX_hg19_chr21.bam')
        self.markdup_path = os.path.join(self.datadir, 'marked_dup_metrics.txt')
        self.markdup_path_low_cover = os.path.join(self.datadir, 'marked_dup_metrics_low_cover.txt')
        self.markdup_path_picard2 = os.path.join(self.datadir, 'marked_dup_metrics_picard2.txt')
        self.markdup_path_picard2_no_histogram = os.path.join(self.datadir, 'marked_dup_metrics_picard2_no_histogram.txt')
        self.markdup_path_picard2_multiple_libraries = os.path.join(self.datadir, 'marked_dup_metrics_picard2_multiple_libraries.txt')
        self.target_path = os.path.join(self.datadir,'SureSelect_All_Exon_V4_Covered_Sorted_chr21.bed')
        self.expected_path = os.path.join(self.datadir, 'expected.json')
        self.expected_alternate_meta = os.path.join(self.datadir, 'expected_alternate_metadata.json')
        self.expected_fast_metrics = os.path.join(self.datadir, 'expected_fast_metrics.json')
        self.expected_no_target = os.path.join(self.datadir, 'expected_no_target.json')
        self.expected_path_downsampled = os.path.join(self.datadir, 'expected_downsampled.json')
        self.expected_path_rs88 = os.path.join(self.datadir, 'expected_downsampled_rs88.json')
        self.expected_metrics_low_cover = os.path.join(self.datadir, 'expected_metrics_low_cover.json')
        self.expected_picard2 = os.path.join(self.datadir, 'expected_picard2.json')
        self.expected_picard2_no_histogram = os.path.join(self.datadir, 'expected_picard2_no_histogram.json')
        self.expected_picard2_multiple_libraries = os.path.join(self.datadir, 'expected_picard2_multiple_libraries.json')
        if os.path.exists(self.OICR_REF):
            self.reference = self.OICR_REF
        elif os.path.exists(self.LOCAL_REF):
            self.reference = self.LOCAL_REF
        else:
            msg = "WARNING: Default reference file %s or %s not found" % \
                  (self.OICR_REF, self.LOCAL_REF)
            print(msg, file=sys.stderr)
            self.reference = None
        self.n_as_mismatch = False
        self.log_path = os.path.join(self.tmpdir, 'test.log')
        self.debug = False
        self.verbose = False
        self.dummy_version = "0.0.0_TEST"
        self.workflow_version = self.dummy_version
        #self.maxDiff = None # uncomment to show the (very long) full output diff

    def assert_default_output_ok(self, actual_path, expected_path):
        expected_variables = {
            "deleted bases": 458,
            "insert size average": 249.8,
            "inserted bases": 315,
            "readsMissingMDtags": 80020,
            "sample level": self.sample_default,
            "total reads": 80020,
            "total target size": 527189,
        }
        self.assert_output_ok(actual_path, expected_path, expected_variables)

    def assert_fast_metric_output_ok(self, actual_path, expected_path):
        expected_variables = {
            "deleted bases": 458,
            "insert size average": 249.8,
            "inserted bases": 315,
            "total reads": 80020,
        }
        self.assert_output_ok(actual_path, expected_path, expected_variables)

    def assert_output_ok(self, actual_path, expected_path, expected_variables):
        with (open(actual_path)) as f: output = json.loads(f.read())
        # do individual sanity checks on some variables
        # helps validate results if expected output JSON file has been changed
        for key in expected_variables.keys():
            expected = expected_variables[key]
            got = output[key]
            try:
                self.assertEqual(expected, got)
            except AssertionError:
                print("\nFailed on metric '"+key+"': Expected", expected, ", got", got,
                      file=sys.stderr)
                raise
        # reference path output depends on local filesystem
        # make test portable by just checking the filename
        self.assertTrue(re.search('/hg19.fa$', output['alignment reference']))
        # now check all output data (aside from the reference path)
        with (open(expected_path)) as f: expected = json.loads(f.read())
        del output['alignment reference']
        del expected['alignment reference']
        self.assertEqual(output, expected)

    def test_alternate_metadata(self):
        # alternate metadata for merged BAM files
        config =  {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path_alternate,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_default,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assert_default_output_ok(out_path, self.expected_alternate_meta)
        qc.cleanup()
        
    def test_default_analysis(self):
        config =  {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_default,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assert_default_output_ok(out_path, self.expected_path)
        qc.cleanup()

    def test_default_analysis_picard2(self):
        config =  {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path_picard2,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_default,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assert_default_output_ok(out_path, self.expected_picard2)
        qc.cleanup()

    def test_default_analysis_picard2_no_histogram(self):
        config = {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path_picard2_no_histogram,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_default,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assert_default_output_ok(out_path, self.expected_picard2_no_histogram)
        qc.cleanup()

    def test_default_analysis_picard2_multiple_libraries(self):
        config = {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path_picard2_multiple_libraries,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_default,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assert_default_output_ok(out_path, self.expected_picard2_multiple_libraries)
        qc.cleanup()

    def test_downsampled_analysis(self):
        config =  {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_level,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out_downsampled.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        # do individual sanity checks on some variables
        # helps validate results if expected output JSON file has been changed
        expected_variables = {
            "inserted bases": 315,
            "reads per start point": 1.003, # downsampled
            "readsMissingMDtags": 9762, # downsampled
            "sample level": self.sample_level,
            "total reads": 80020,
            "total target size": 527189,
        }
        for key in expected_variables.keys():
            expected = expected_variables[key]
            got = output[key]
            try:
                self.assertEqual(expected, got)
            except AssertionError:
                print("\nFailed on metric '"+key+"': Expected", expected, ", got", got,
                      file=sys.stderr)
                raise
        # reference path output depends on local filesystem
        # make test portable by just checking the filename
        self.assertTrue(re.search('/hg19.fa$', output['alignment reference']))
        # now check all output data (aside from the reference)
        with (open(self.expected_path_downsampled)) as f: expected = json.loads(f.read())
        del output['alignment reference']
        self.assertEqual(output, expected)
        qc.cleanup()

    def test_fast_metric_script(self):
        # test the fast metric command-line script
        filename = 'write_fast_metrics.py'
        relative_path = os.path.join(os.path.dirname(__file__), os.pardir, 'bin', filename)
        script = os.path.realpath(relative_path)
        out_path = os.path.join(self.tmpdir, 'fast_metric_script_out.json')
        args = [
            script,
            '--bam', self.bam_path,
            '--insert-max', str(self.insert_max),
            '--out', out_path,
            '--reference', self.reference,
        ]
        if self.n_as_mismatch:
            args.append('--n-as-mismatch')
        if self.debug:
            args.append('--debug')
        if self.verbose:
            args.append('--verbose')
        result = subprocess.run(args)
        try:
            result.check_returncode()
        except subprocess.CalledProcessError:
            print("STANDARD OUTPUT:", result.stdout, file=sys.stderr)
            print("STANDARD ERROR:", result.stderr, file=sys.stderr)
            raise
        self.assertTrue(os.path.exists(out_path))
        self.assert_fast_metric_output_ok(out_path, self.expected_fast_metrics)

    def test_fast_metric_writer(self):
        # test the standalone fast_metric_writer class
        config =  {
            fast_metric_writer.CONFIG_KEY_BAM: self.bam_path,
            fast_metric_writer.CONFIG_KEY_DEBUG: self.debug,
            fast_metric_writer.CONFIG_KEY_INSERT_MAX: self.insert_max,
            fast_metric_writer.CONFIG_KEY_LOG: self.log_path,
            fast_metric_writer.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            fast_metric_writer.CONFIG_KEY_REFERENCE: self.reference,
            fast_metric_writer.CONFIG_KEY_VERBOSE: self.verbose,
        }
        fmw = fast_metric_writer(config)
        out_path = os.path.join(self.tmpdir, 'out_fast_metrics.json')
        fmw.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        self.assert_fast_metric_output_ok(out_path, self.expected_fast_metrics)

    def test_missing_inputs(self):
        # test possible missing inputs:
        # - ESTIMATED_LIBRARY_SIZE in mark duplicates text
        # - FFQ/LFQ in samtools stats
        config =  {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_default,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        # for low-coverage runs, ESTIMATED_LIBRARY_SIZE value is missing from mark duplicates text
        # test input file also has variant '## METRICS CLASS ...' line
        metrics_found = qc.read_mark_dup(self.markdup_path_low_cover)
        with (open(self.expected_metrics_low_cover)) as f: metrics_expected = json.loads(f.read())
        # Found/expected HISTOGRAM keys are integers and strings, respectively.
        # (Annoyingly, JSON format insists dictionary keys must be strings)
        histogram_found = metrics_found['HISTOGRAM']
        histogram_expected = metrics_expected['HISTOGRAM']
        self.assertEqual(len(histogram_found), len(histogram_expected))
        for histogram_type in histogram_found.keys():
            for key in histogram_found[histogram_type]:
                self.assertEqual(histogram_found[histogram_type][key], histogram_expected[histogram_type][str(key)])
        del metrics_found['HISTOGRAM']
        del metrics_expected['HISTOGRAM']
        self.assertEqual(metrics_found, metrics_expected)
        # test empty FFQ/LFQ result from samtools stats; may occur for small input datasets
        # requires a fast_metric_finder object
        fast_finder = fast_metric_finder(self.bam_path,
                                         self.reference,
                                         self.insert_max,
                                         self.n_as_mismatch,
                                         qc.logger)
        fq_result = fast_finder.fq_stats([])
        fq_expected = ({},{})
        self.assertEqual(fq_expected, fq_result)
        qc.cleanup()

    def test_random_seed(self):
        # test downsampling with a different random seed from the default
        config =  {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: self.target_path,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: 88,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_level,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out_downsampled_88.json')
        qc.write_output(out_path)
        with (open(out_path)) as f: output = json.loads(f.read())
        with (open(self.expected_path_rs88)) as f: expected = json.loads(f.read())
        # do not test the alignment reference local path
        del expected['alignment reference']
        del output['alignment reference']
        self.assertEqual(output, expected)
        qc.cleanup()

    def test_main_script(self):
        relative_path = os.path.join(os.path.dirname(__file__), os.pardir, 'bin', 'run_bam_qc.py')
        script = os.path.realpath(relative_path)
        out_path = os.path.join(self.tmpdir, 'script_out.json')
        args = [
            script,
            '--bam', self.bam_path,
            '--target', self.target_path,
            '--insert-max', str(self.insert_max),
            '--metadata', self.metadata_path,
            '--mark-duplicates', self.markdup_path,
            '--out', out_path,
            '--skip-below-mapq', str(self.quality),
            '--reference', self.reference,
            '--temp-dir', self.tmpdir,
            '--workflow-version', self.workflow_version
        ]
        if self.n_as_mismatch:
            args.append('--n-as-mismatch')
        if self.debug:
            args.append('--debug')
        if self.verbose:
            args.append('--verbose')
        result = subprocess.run(args)
        try:
            result.check_returncode()
        except subprocess.CalledProcessError:
            print("STANDARD OUTPUT:", result.stdout, file=sys.stderr)
            print("STANDARD ERROR:", result.stderr, file=sys.stderr)
            raise
        self.assert_default_output_ok(out_path, self.expected_path)

    def test_version_updater(self):
        # 'update' copies of the test files with a dummy package version
        # then check the dummy version has been written correctly
        data_dir = self.datadir
        updater = version_updater()
        temp_paths = []
        for name in updater.get_filenames():
            dest = os.path.join(self.tmpdir, name)
            temp_paths.append(dest)
            shutil.copyfile(os.path.join(data_dir, name), dest)
        updater.update_files(self.tmpdir, self.dummy_version)
        for temp_path in temp_paths:
            with (open(temp_path)) as f: data = json.loads(f.read())
            self.assertEqual(data[updater.PACKAGE_VERSION_KEY], self.dummy_version)

    def test_without_target(self):
        config =  {
            bam_qc.CONFIG_KEY_BAM: self.bam_path,
            bam_qc.CONFIG_KEY_DEBUG: self.debug,
            bam_qc.CONFIG_KEY_TARGET: None,
            bam_qc.CONFIG_KEY_INSERT_MAX: self.insert_max,
            bam_qc.CONFIG_KEY_LOG: self.log_path,
            bam_qc.CONFIG_KEY_METADATA: self.metadata_path,
            bam_qc.CONFIG_KEY_MARK_DUPLICATES: self.markdup_path,
            bam_qc.CONFIG_KEY_N_AS_MISMATCH: self.n_as_mismatch,
            bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: self.quality,
            bam_qc.CONFIG_KEY_RANDOM_SEED: None,
            bam_qc.CONFIG_KEY_REFERENCE: self.reference,
            bam_qc.CONFIG_KEY_SAMPLE: self.sample_default,
            bam_qc.CONFIG_KEY_TEMP_DIR: self.tmpdir,
            bam_qc.CONFIG_KEY_VERBOSE: self.verbose,
            bam_qc.CONFIG_KEY_WORKFLOW_VERSION: self.workflow_version
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        with (open(out_path)) as f: output = json.loads(f.read())
        # do individual sanity checks on some variables
        # helps validate results if expected output JSON file has been changed
        expected_variables = {
            "inserted bases": 315,
            "reads per start point": 1.031,
            "readsMissingMDtags": 80020,
            "sample level": self.sample_default,
            "total reads": 80020,
            "total target size": None,
        }
        for key in expected_variables.keys():
            expected = expected_variables[key]
            got = output[key]
            try:
                self.assertEqual(expected, got)
            except AssertionError:
                print("\nFailed on metric '"+key+"': Expected", expected, ", got", got,
                      file=sys.stderr)
                raise
        # reference path output depends on local filesystem
        # make test portable by just checking the filename
        self.assertTrue(re.search('/hg19.fa$', output['alignment reference']))
        # now check all output data (aside from the reference path)
        with (open(self.expected_no_target)) as f: expected = json.loads(f.read())
        del output['alignment reference']
        self.assertEqual(output, expected)
        qc.cleanup()

    def tearDown(self):
        self.tmp.cleanup()

if __name__ == '__main__':
    unittest.main()
