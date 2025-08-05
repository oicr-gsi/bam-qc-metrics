#!/usr/bin/env python3

"""
   Main script to compute BAM QC metrics
   This is a modified version which accepts inputs for multiple lanes and
   does not produce all metrics which original bamQC produces. We mostly
   interested in CIGAR analysis implemented in original code, other metrics come from
   external sources and get parsed by bam-qc-metrics parsing functions
"""

import argparse
import os
import re
import sys
import tempfile

from bam_qc_metrics_lite import bam_qc_lite, read_package_version, validator

DEFAULT_INSERT_MAX = 1500
DEFAULT_DOWNSAMPLED_READS = 500000

def validate_args(args):
    valid = True

    # flip valid from True to False if a check is failed; never flip back to True
    if args.insert_max is not None:
        valid = valid and validator.validate_positive_integer(args.insert_max, 'Max insert size')
    for path_arg in (args.samstats, args.mark_duplicates, args.coverage, args.metadata):
        if path_arg is not None:
            valid = valid and validator.validate_input_file(path_arg)
    if args.targeted is not None:
        valid = valid and validator.validate_input_file(args.targeted)
    if args.out != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        valid = valid and validator.validate_output_dir(parent_path)
    if args.log_path is not None:
        parent_path = os.path.abspath(os.path.join(args.log_path, os.pardir))
        valid = valid and validator.validate_output_dir(parent_path)
    if args.temp_dir is not None:
        valid = valid and validator.validate_output_dir(args.temp_dir)
    if args.workflow_version is not None:
        # version string should be of the form 0.12.34 or (for instance) 0.12.34_alpha
        pattern = r'^\d+\.\d+\.\d+\w*$'
        if not re.match(pattern, args.workflow_version):
            sys.stderr.write("ERROR: Workflow version does not match pattern "+pattern+"\n")
            valid = False
    return valid

def main():
    parser = argparse.ArgumentParser(description='QC for BAM files.')
    parser.add_argument('-b', '--bam', metavar='PATH',
                        help='Path to input BAM file. Optional.')
    parser.add_argument('-s', '--samstats', metavar='PATH', required=True,
                        help='Path to text file output by samtools stats.')
    parser.add_argument('-S', '--downsample-to', metavar='INT', default=DEFAULT_DOWNSAMPLED_READS,
                        help='Reads to use for calculation of metrics, if requested.')
    parser.add_argument('-d', '--mark-duplicates', metavar='PATH',
                        help='Path to text file output by Samblaster. Optional.')
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata. Optional.')
    parser.add_argument('-c', '--coverage', metavar='PATH',
                        help='Path to JSON file containing coverage histogram. Optional.')
    parser.add_argument('-t', '--targeted', metavar='PATH',
                        help='Path to optional regions.gz file containing targeted coverage data. Optional.')
    parser.add_argument('-tf', '--target-path', metavar='PATH', required=False,
                        help='Path to target .bed file. Optional.')
    parser.add_argument('-ot', '--reads-on-target', required=False,
                        help='integer showing on target reads int. Optional.')
    parser.add_argument('-D', '--debug', action='store_true',
                        help='Most verbose; write messages of priority DEBUG and higher to log')
    parser.add_argument('-i', '--insert-max', metavar='INT', default=DEFAULT_INSERT_MAX,
                        help='Maximum expected value for insert size; higher values will be ' +
                        'counted as abnormal. Optional; default = %i.' % DEFAULT_INSERT_MAX)
    parser.add_argument('-g', '--cigar_metrics', metavar='PATH',
                        help='Path to optional .json file with lane-level CIGAR metrics.')
    parser.add_argument('-l', '--log-path', metavar='PATH', help='Path of file where log output ' +
                        'will be appended. Optional, defaults to STDERR.')
    parser.add_argument('-o', '--out', metavar='PATH', required=True,
                        help='Path for JSON output, or - for STDOUT. Required.')
    parser.add_argument('-p', '--profile', action='store_true', help='Write runtime profile to ' +
                        'STDOUT. For development use only. Should not be combined with writing ' +
                        'JSON to STDOUT.')
    parser.add_argument('-r', '--reference', metavar='PATH',
                        help='Path to FASTA reference used to align the BAM file. Used for annotation only.')
    parser.add_argument('-u', '--unique-reads', metavar='INT',
                        help='Set unique reads calculated externally, this needs to come from downsampled bam file.')
    parser.add_argument('-R', '--random-seed', metavar='INT', help='Set sampling random seed to ' +
                        'INT. Has no effect if --sample-rate not specified. Optional; if not ' +
                        'given, a default seed will be used.')
    parser.add_argument('-T', '--temp-dir', metavar='PATH', help='Directory for temporary output ' +
                        'files; optional, defaults to %s (the current system tempdir).'
                        % tempfile.gettempdir())
    parser.add_argument('-v', '--version', action='version',
                        version=read_package_version(),
                        help='Print the version number of bam-qc-metrics and exit')
    parser.add_argument('-V', '--verbose', action='store_true',
                        help='More verbose; write messages of priority INFO and higher to log')
    parser.add_argument('-w', '--workflow-version', metavar='VERSION',
                        help='Version of the workflow being used to run bam-qc-metrics. ' +
                        'Optional. If given, will be recorded in JSON output.')
    args = parser.parse_args()
    if not validate_args(args):
        print("For usage, run with -h or --help")
        exit(1)
    insert_max = None if args.insert_max is None else int(args.insert_max)
    random_seed = None if args.random_seed is None else int(args.random_seed)
    reads_on_target = None if args.reads_on_target is None else int(args.reads_on_target)
    unique_reads = None if args.unique_reads is None else int(args.unique_reads)
    downsample_to = None if args.downsample_to is None else int(args.downsample_to)
    bam_file = args.bam
    bam_file = bam_file if bam_file is not None and str(bam_file).endswith(".bam") else None

    config = {
        bam_qc_lite.CONFIG_KEY_BAM: bam_file,
        bam_qc_lite.CONFIG_KEY_SAMSTATS: args.samstats,
        bam_qc_lite.CONFIG_KEY_COVERAGE: args.coverage,
        bam_qc_lite.CONFIG_KEY_METADATA: args.metadata,
        bam_qc_lite.CONFIG_KEY_MARK_DUPLICATES: args.mark_duplicates,
        bam_qc_lite.CONFIG_KEY_TARGET_COVERAGE: args.targeted,
        bam_qc_lite.CONFIG_KEY_DOWNSAMPLED_READS: downsample_to,
        bam_qc_lite.CONFIG_KEY_CIGAR_METRICS: args.cigar_metrics,
        bam_qc_lite.CONFIG_KEY_TARGET_PATH: args.target_path,
        bam_qc_lite.CONFIG_KEY_READS_ON_TARGET: reads_on_target,
        bam_qc_lite.CONFIG_KEY_UNIQUE_READS: unique_reads,
        bam_qc_lite.CONFIG_KEY_REFERENCE: args.reference,
        bam_qc_lite.CONFIG_KEY_DEBUG: args.debug,
        bam_qc_lite.CONFIG_KEY_INSERT_MAX: insert_max,
        bam_qc_lite.CONFIG_KEY_LOG: args.log_path,
        bam_qc_lite.CONFIG_KEY_RANDOM_SEED: random_seed,
        bam_qc_lite.CONFIG_KEY_TEMP_DIR: args.temp_dir,
        bam_qc_lite.CONFIG_KEY_VERBOSE: args.verbose,
        bam_qc_lite.CONFIG_KEY_WORKFLOW_VERSION: args.workflow_version
    }

    '''Before running, validate that we have all arguments in a good shape'''

    '''Assume we have inputs from multiple lanes, use threading to process inputs separately by lane'''
    qc = bam_qc_lite(config)

    '''Merge if we have > 1 lane, otherwise return report for a single set of inputs'''

    qc.write_output(args.out)

if __name__ == "__main__":
    main()
