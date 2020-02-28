#!/usr/bin/env python3

"""Main script to compute BAM QC metrics"""

import argparse, cProfile, os, re, sys, tempfile
from bam_qc_metrics import bam_qc, read_package_version, validator

DEFAULT_INSERT_MAX = 1500
DEFAULT_SAMPLE_LEVEL = 1100000

def validate_args(args):
    valid = True

    # flip valid from True to False if a check is failed; never flip back to True
    if args.skip_below_mapq != None:
        valid = validator.validate_positive_integer(args.skip_below_mapq, 'Quality score')
    if args.insert_max != None:
        valid = valid and validator.validate_positive_integer(args.insert_max, 'Max insert size')
    if args.sample != None:
        valid = valid and validator.validate_positive_integer(args.sample, 'Downsampling level')
        valid = valid and validator.validate_sample_level(args.all_reads, args.sample)
    if args.downsampled_bam != None:
        if args.all_reads != None or args.sample != None:
            valid = False
            sys.stderr.write("ERROR: -S/--downsampled-bam argument incompatible with "+\
                             "-a/--all-reads or -s/--sample")
        valid = valid and validator.validate_input_file(args.downsampled_bam)
    if args.bam == None:
        valid = False
        sys.stderr.write("ERROR: -b/--bam argument is required\n")
    for path_arg in (args.bam, args.target, args.metadata, args.mark_duplicates, args.reference):
        if path_arg != None:
            valid = valid and validator.validate_input_file(path_arg)
    if args.out != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        valid = valid and validator.validate_output_dir(parent_path)
    if args.log_path != None:
        parent_path = os.path.abspath(os.path.join(args.log_path, os.pardir))
        valid = valid and validator.validate_output_dir(parent_path)
    if args.temp_dir != None:
        valid = valid and validator.validate_output_dir(args.temp_dir)
    if args.workflow_version != None:
        # version string should be of the form 0.12.34 or (for instance) 0.12.34_alpha
        pattern = r'^[0-9]+\.[0-9]+\.[0-9]+\w*$'
        if not re.match(pattern, args.workflow_version):
            sys.stderr.write("ERROR: Workflow version does not match pattern "+pattern+"\n")
            valid = False
    return valid

def main():
    parser = argparse.ArgumentParser(description='QC for BAM files.')
    parser.add_argument('-a', '--all-reads', action='store_true', help='Do not apply downsampling; '+\
                        'use all reads as input to all QC metrics. Incompatible with --sample.')
    parser.add_argument('-b', '--bam', metavar='PATH', required=True,
                        help='Path to input BAM file. Required.')
    parser.add_argument('-d', '--mark-duplicates', metavar='PATH',
                        help='Path to text file output by Picard MarkDuplicates. Optional.')
    parser.add_argument('-D', '--debug', action='store_true',
                        help='Most verbose; write messages of priority DEBUG and higher to log')
    parser.add_argument('-i', '--insert-max', metavar='INT', default=DEFAULT_INSERT_MAX,
                        help='Maximum expected value for insert size; higher values will be '+\
                        'counted as abnormal. Optional; default = %i.' % DEFAULT_INSERT_MAX)
    parser.add_argument('-l', '--log-path', metavar='PATH', help='Path of file where log output '+\
                        'will be appended. Optional, defaults to STDERR.')
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata. Optional.')
    parser.add_argument('-n', '--n-as-mismatch', action='store_true',
                        help='Record N calls as mismatches in mismatch-per-cycle counts. '+\
                        'Only relevant if a reference is given with -r.')
    parser.add_argument('-o', '--out', metavar='PATH', required=True,
                        help='Path for JSON output, or - for STDOUT. Required.')
    parser.add_argument('-p', '--profile', action='store_true', help='Write runtime profile to '+\
                        'STDOUT. For development use only. Should not be combined with writing '+\
                        'JSON to STDOUT.')
    parser.add_argument('-q', '--skip-below-mapq', metavar='QSCORE',
                        help='Threshold to skip reads with low alignment quality. Optional.')
    parser.add_argument('-r', '--reference', metavar='PATH',
                        help='Path to FASTA reference used to align the BAM file. Used to find '+\
                        'mismatches by cycle using samtools. Optional; if not supplied, '+\
                        'mismatches by cycle will be empty.')
    parser.add_argument('-R', '--random-seed', metavar='INT', help='Set sampling random seed to '+\
                        'INT. Has no effect if --sample-rate not specified. Optional; if not '+\
                        'given, a default seed will be used.')
    parser.add_argument('-s', '--sample', metavar='INT',
                        help='Sample a total of INT reads from the BAM file, for input to slower '+\
                        'QC metrics. Defaults to 1.1 million. Incompatible with --all-reads.')
    parser.add_argument('-S', '--downsampled-bam', metavar='PATH',
                        help='Downsampled BAM file for input to slow QC metrics. Incompatible with '+\
                        '--all-reads and --sample.')
    parser.add_argument('-t', '--target', metavar='PATH',
                        help='Path to target BED file, containing targets to calculate coverage '+\
                        'against. Optional. If given, must be sorted in same order as BAM file. '+\
                        'If not given, bedtools coverage metrics will be omitted.')
    parser.add_argument('-T', '--temp-dir', metavar='PATH', help='Directory for temporary output '+\
                        'files; optional, defaults to %s (the current system tempdir).' \
                        % tempfile.gettempdir())
    parser.add_argument('-v', '--version', action='version',
                        version=read_package_version(),
                        help='Print the version number of bam-qc-metrics and exit')
    parser.add_argument('-V', '--verbose', action='store_true',
                        help='More verbose; write messages of priority INFO and higher to log')
    parser.add_argument('-w', '--workflow-version', metavar='VERSION',
                        help='Version of the workflow being used to run bam-qc-metrics. '+\
                        'Optional. If given, will be recorded in JSON output.')
    args = parser.parse_args()
    if not validate_args(args):
        print("For usage, run with -h or --help")
        exit(1)
    skip_below_mapq = None if args.skip_below_mapq == None else int(args.skip_below_mapq)
    insert_max = None if args.insert_max == None else int(args.insert_max)
    random_seed = None if args.random_seed == None else int(args.random_seed)
    if args.all_reads:
        sample = None
    else:
        sample = DEFAULT_SAMPLE_LEVEL if args.sample == None else int(args.sample)
    config = {
        bam_qc.CONFIG_KEY_BAM: args.bam,
        bam_qc.CONFIG_KEY_DEBUG: args.debug,
        bam_qc.CONFIG_KEY_DOWNSAMPLED_BAM: args.downsampled_bam,
        bam_qc.CONFIG_KEY_TARGET: args.target,
        bam_qc.CONFIG_KEY_INSERT_MAX: insert_max,
        bam_qc.CONFIG_KEY_LOG: args.log_path,
        bam_qc.CONFIG_KEY_METADATA: args.metadata,
        bam_qc.CONFIG_KEY_MARK_DUPLICATES: args.mark_duplicates,
        bam_qc.CONFIG_KEY_N_AS_MISMATCH: args.n_as_mismatch,
        bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: skip_below_mapq,
        bam_qc.CONFIG_KEY_RANDOM_SEED: random_seed,
        bam_qc.CONFIG_KEY_REFERENCE: args.reference,
        bam_qc.CONFIG_KEY_SAMPLE: sample,
        bam_qc.CONFIG_KEY_TEMP_DIR: args.temp_dir,
        bam_qc.CONFIG_KEY_VERBOSE: args.verbose,
        bam_qc.CONFIG_KEY_WORKFLOW_VERSION: args.workflow_version
    }
    if args.profile:
        # sort order = 2, sorts profile by cumulative time
        cProfile.runctx('bam_qc(config).write_output(out_path)',
                        {'bam_qc': bam_qc, 'config': config, 'out_path': args.out},
                        {},
                        None,
                        2)
    else:
        qc = bam_qc(config)
        qc.write_output(args.out)

if __name__ == "__main__":
    main()
