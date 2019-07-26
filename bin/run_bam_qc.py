#! /usr/bin/env python3

"""Main script to compute BAM QC metrics"""

import argparse, cProfile, os, sys
from bam_qc_metrics import bam_qc

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
    if args.trim_quality != None:
        valid = validate_positive_integer(args.trim_quality, 'Quality score')
    if args.insert_max != None:
        valid = valid and validate_positive_integer(args.insert_max, 'Max insert size')
    if args.sample_rate != None:
        valid = valid and validate_positive_integer(args.sample_rate, 'Downsampling rate')
    for path_arg in (args.bam, args.target, args.metadata, args.mark_duplicates):
        if path_arg != None:
            valid = valid and validate_input_file(path_arg)
    if args.profile != None and args.profile != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        valid = valid and validate_output_dir(parent_path)
    if args.out != '-':
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        valid = valid and validate_output_dir(parent_path)
    if args.temp_dir != None:
        valid = valid and validate_output_dir(args.temp_dir)
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
    parser.add_argument('-o', '--out', metavar='PATH', required=True,
                        help='Path for JSON output, or - for STDOUT. Required.')
    parser.add_argument('-p', '--profile', metavar='PATH',
                        help='Write profile stats to PATH, or - for STDOUT. WARNING: This will significantly increase runtime.')
    parser.add_argument('-q', '--trim-quality', metavar='QSCORE',
                        help='Samtools threshold for trimming alignment quality. Optional.')
    parser.add_argument('-s', '--sample-rate', metavar='INT',
                        help='Sample every Nth read, where N is the argument. Optional, defaults to 1 (no sampling).')
    parser.add_argument('-t', '--target', metavar='PATH',
                        help='Path to target BED file, containing targets to calculate coverage '+\
                        'against. Optional; if given, must be sorted in same order as BAM file.')
    parser.add_argument('-T', '--temp-dir', metavar='PATH', help='Directory for temporary output files; optional, defaults to /tmp')
    args = parser.parse_args()
    if not validate_args(args): exit(1)
    trim_quality = None if args.trim_quality == None else int(args.trim_quality)
    insert_max = None if args.insert_max == None else int(args.insert_max)
    sample_rate = None if args.sample_rate == None else int(args.sample_rate)
    if args.profile != None:
        method_args = (args.bam, args.target, insert_max, args.metadata, args.mark_duplicates, trim_quality, sample_rate, args.temp_dir)
        if args.profile == '-':
            params = [args.bam, args.target, insert_max, args.metadata, args.mark_duplicates, trim_quality, sample_rate, args.temp_dir]
            method_args = tuple([str(p) for p in params])
            cmd = "bam_qc('%s', '%s', %s, %s, %s, %s, %s, %s)" % method_args
            cProfile.run(cmd)
        else:
            sys.stderr.write(args.profile+"\n")
            profile_path = str(args.profile)
            qc = cProfile.run(bam_qc(*method_args), profile_path)
            qc.write_output(args.out)
    else:
        qc = bam_qc(args.bam, args.target, insert_max, args.metadata, args.mark_duplicates, trim_quality, sample_rate, args.temp_dir)
        qc.write_output(args.out)

if __name__ == "__main__":
    main()
