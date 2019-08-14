#! /usr/bin/env python3

import argparse, json, os
from bam_qc_metrics import version_updater

def main():
    updater = version_updater()
    parser = argparse.ArgumentParser(description='Convenience script to package version '+\
                                     'in JSON test data')
    parser.add_argument('--run', action='store_true', help='Run update on JSON files')
    parser.add_argument('--version', action='version', version=updater.package_version)
    args = parser.parse_args()
    if args.run:
        updater.update_files()
    else:
        print("Run with -h/--help for usage")

if __name__ == "__main__":
    main()
