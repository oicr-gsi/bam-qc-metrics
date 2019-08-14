#! /usr/bin/env python3

"""Main script to compute BAM QC metrics"""

import argparse, os, re, sys, tempfile
from bam_qc_metrics import bam_qc, read_package_version

DEFAULT_INSERT_MAX = 1500

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
    if args.skip_below_mapq != None:
        valid = validate_positive_integer(args.skip_below_mapq, 'Quality score')
    if args.insert_max != None:
        valid = valid and validate_positive_integer(args.insert_max, 'Max insert size')
    if args.sample_rate != None:
        valid = valid and validate_positive_integer(args.sample_rate, 'Downsampling rate')
    if args.bam == None:
        valid = False
        sys.stderr.write("ERROR: -b/--bam argument is required\n")
    for path_arg in (args.bam, args.target, args.metadata, args.mark_duplicates, args.reference):
        if path_arg != None:
            valid = valid and validate_input_file(path_arg)
    if args.out != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        valid = valid and validate_output_dir(parent_path)
    if args.temp_dir != None:
        valid = valid and validate_output_dir(args.temp_dir)
    if args.workflow_version != None:
        # version string should be of the form 0.12.34 or (for instance) 0.12.34_alpha
        pattern = r'^[0-9]+\.[0-9]+\.[0-9]+\w*$'
        if not re.match(pattern, args.workflow_version):
            sys.stderr.write("ERROR: Workflow version does not match pattern "+pattern+"\n")
            valid = False
    return valid

def main():
    parser = argparse.ArgumentParser(description='QC for BAM files.')
    parser.add_argument('-b', '--bam', metavar='PATH', required=True,
                        help='Path to input BAM file. Required.')
    parser.add_argument('-d', '--mark-duplicates', metavar='PATH',
                        help='Path to text file output by Picard MarkDuplicates. Optional.')
    parser.add_argument('-i', '--insert-max', metavar='INT', default=DEFAULT_INSERT_MAX,
                        help='Maximum expected value for insert size; higher values will be '+\
                        'counted as abnormal. Optional; default = %i.' % DEFAULT_INSERT_MAX)
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata. Optional.')
    parser.add_argument('-n', '--n-as-mismatch', action='store_true',
                        help='Record N calls as mismatches in mismatch-per-cycle counts. '+\
                        'Only relevant if a reference is given with -r.')
    parser.add_argument('-o', '--out', metavar='PATH', required=True,
                        help='Path for JSON output, or - for STDOUT. Required.')
    parser.add_argument('-q', '--skip-below-mapq', metavar='QSCORE',
                        help='Threshold to skip reads with low alignment quality. Optional.')
    parser.add_argument('-r', '--reference', metavar='PATH',
                        help='Path to FASTA reference used to align the BAM file. Used to find '+\
                        'mismatches by cycle using samtools. Optional; if not supplied, '+\
                        'mismatches by cycle will be empty.')
    parser.add_argument('-s', '--sample-rate', metavar='INT',
                        help='Sample every Nth read, where N is the argument. Optional, defaults to 1 (no sampling).')
    parser.add_argument('-t', '--target', metavar='PATH',
                        help='Path to target BED file, containing targets to calculate coverage '+\
                        'against. Optional; if given, must be sorted in same order as BAM file.')
    parser.add_argument('-T', '--temp-dir', metavar='PATH', help='Directory for temporary output files; optional, defaults to %s (the current system tempdir).' % tempfile.gettempdir())
    parser.add_argument('-v', '--version', action='version',
                        version=read_package_version(),
                        help='Print the version number of bam-qc-metrics and exit')
    parser.add_argument('-V', '--verbose', action='store_true', help='Print additional messages to STDERR')
    parser.add_argument('-w', '--workflow-version', metavar='VERSION',
                        help='Version of the workflow being used to run bam-qc-metrics. '+\
                        'Optional. If given, will be recorded in JSON output.')
    args = parser.parse_args()
    if not validate_args(args):
        print("For usage, run with -h or --help")
        exit(1)
    skip_below_mapq = None if args.skip_below_mapq == None else int(args.skip_below_mapq)
    insert_max = None if args.insert_max == None else int(args.insert_max)
    sample_rate = None if args.sample_rate == None else int(args.sample_rate)
    config = {
        bam_qc.CONFIG_KEY_BAM: args.bam,
        bam_qc.CONFIG_KEY_TARGET: args.target,
        bam_qc.CONFIG_KEY_INSERT_MAX: insert_max,
        bam_qc.CONFIG_KEY_METADATA: args.metadata,
        bam_qc.CONFIG_KEY_MARK_DUPLICATES: args.mark_duplicates,
        bam_qc.CONFIG_KEY_N_AS_MISMATCH: args.n_as_mismatch,
        bam_qc.CONFIG_KEY_SKIP_BELOW_MAPQ: skip_below_mapq,
        bam_qc.CONFIG_KEY_REFERENCE: args.reference,
        bam_qc.CONFIG_KEY_SAMPLE_RATE: sample_rate,
        bam_qc.CONFIG_KEY_TEMP_DIR: args.temp_dir,
        bam_qc.CONFIG_KEY_VERBOSE: args.verbose,
        bam_qc.CONFIG_KEY_WORKFLOW_VERSION: args.workflow_version
    }
    qc = bam_qc(config)
    qc.write_output(args.out)

if __name__ == "__main__":
    main()
