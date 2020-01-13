#! /usr/bin/env python3

import argparse, json, os
from bam_qc_metrics import version_updater

def main():
    updater = version_updater()
    parser = argparse.ArgumentParser(description='Convenience script to update package version '+\
                                     'in JSON test data')
    parser.add_argument('--dir', metavar="PATH", help="Input directory, defaults to './data/'",
                        default='./data/', dest='input_dir')
    parser.add_argument('--run', action='store_true', help='Run update on JSON files')
    parser.add_argument('--show-version', action='version', version=updater.package_version,
                        help="Show default package version and exit")
    parser.add_argument('--set-version', metavar='VERSION', default=updater.package_version,
                        help="Version string to use; defaults to the package version as shown "+\
                        "by --show-version")
    args = parser.parse_args()
    if args.run:
        updater.update_files(args.input_dir, args.set_version)
    else:
        print("Run with -h/--help for usage")

if __name__ == "__main__":
    main()
