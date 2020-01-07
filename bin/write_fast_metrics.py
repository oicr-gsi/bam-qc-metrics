#!/usr/bin/env python3

"""Script to compute faster BAM QC metrics derived from samtools stats"""

import argparse, cProfile, os, re, sys, tempfile
from bam_qc_metrics import fast_metric_writer, read_package_version, validator

DEFAULT_INSERT_MAX = 1500

def validate_args(args):
    valid = True
    # flip valid from True to False if a check is failed; never flip back to True
    if args.insert_max != None:
        valid = valid and validator.validate_positive_integer(args.insert_max, 'Max insert size')
    if args.bam == None:
        valid = False
        sys.stderr.write("ERROR: -b/--bam argument is required\n")
    for path_arg in (args.bam, args.reference):
        if path_arg != None:
            valid = valid and validator.validate_input_file(path_arg)
    if args.out != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        valid = valid and validator.validate_output_dir(parent_path)
    if args.log_path != None:
        parent_path = os.path.abspath(os.path.join(args.log_path, os.pardir))
        valid = valid and validator.validate_output_dir(parent_path)
    return valid

def main():
    parser = argparse.ArgumentParser(description='Compute low-overhead QC metrics for BAM files.')
    parser.add_argument('-b', '--bam', metavar='PATH', required=True,
                        help='Path to input BAM file. Required.')
    parser.add_argument('-D', '--debug', action='store_true',
                        help='Most verbose; write messages of priority DEBUG and higher to log')
    parser.add_argument('-i', '--insert-max', metavar='INT', default=DEFAULT_INSERT_MAX,
                        help='Maximum expected value for insert size; higher values will be '+\
                        'counted as abnormal. Optional; default = %i.' % DEFAULT_INSERT_MAX)
    parser.add_argument('-l', '--log-path', metavar='PATH', help='Path of file where log output '+\
                        'will be appended. Optional, defaults to STDERR.')
    parser.add_argument('-n', '--n-as-mismatch', action='store_true',
                        help='Record N calls as mismatches in mismatch-per-cycle counts. '+\
                        'Only relevant if a reference is given with -r.')
    parser.add_argument('-o', '--out', metavar='PATH', required=True,
                        help='Path for JSON output, or - for STDOUT. Required.')
    parser.add_argument('-r', '--reference', metavar='PATH',
                        help='Path to FASTA reference used to align the BAM file. Used to find '+\
                        'mismatches by cycle using samtools. Optional; if not supplied, '+\
                        'mismatches by cycle will be empty.')
    parser.add_argument('-v', '--version', action='version',
                        version=read_package_version(),
                        help='Print the version number of bam-qc-metrics and exit')
    parser.add_argument('-V', '--verbose', action='store_true',
                        help='More verbose; write messages of priority INFO and higher to log')
    args = parser.parse_args()
    if not validate_args(args):
        print("For usage, run with -h or --help")
        exit(1)
    insert_max = None if args.insert_max == None else int(args.insert_max)
    config = {
        fast_metric_writer.CONFIG_KEY_BAM: args.bam,
        fast_metric_writer.CONFIG_KEY_DEBUG: args.debug,
        fast_metric_writer.CONFIG_KEY_INSERT_MAX: insert_max,
        fast_metric_writer.CONFIG_KEY_LOG: args.log_path,
        fast_metric_writer.CONFIG_KEY_N_AS_MISMATCH: args.n_as_mismatch,
        fast_metric_writer.CONFIG_KEY_REFERENCE: args.reference,
        fast_metric_writer.CONFIG_KEY_VERBOSE: args.verbose,
    }
    fast_metric_writer(config).write_output(args.out)

if __name__ == "__main__":
    main()
