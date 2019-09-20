#! /usr/bin/env python3

"""Script to collate metadata into JSON format"""

import argparse, json

def main():
    desc ='Collate metadata for BAM QC and write in JSON format.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--barcode', help="barcode metadata value")
    parser.add_argument('--instrument', help="instrument metadata value")
    parser.add_argument('--lane', help="lane metadata value")
    parser.add_argument('--library', help="library metadata value")
    parser.add_argument('--run-name', help="run_name metadata value")
    parser.add_argument('--sample', help="sample metadata value")
    parser.add_argument('--out', help="Path for JSON output")
    args = parser.parse_args()
    data = {
        'barcode': args.barcode,
        'instrument': args.instrument,
        'lane': args.lane,
        'library': args.library,
        'run name': args.run_name,
        'sample': args.sample
    }
    out = open(args.out, 'w')
    print(json.dumps(data), file=out)
    out.close()

if __name__ == "__main__":
    main()
