#! /usr/bin/env python3

"""
Script to diff old and new BAMQC JSON outputs, generated by bamqc.pl and bam_qc.py respectively

JSON data structure has changed; maps old fields to new fields, and disregards obsolete ones
"""

import argparse, json, sys


class diff:

    OLD_KEYS = set(['aligned bases', 'alignedCount', 'average read length', 'barcode', 'collapsed bases covered', 'deleted bases', 'deletionCount', 'hard clip bases', 'hardClipCount', 'insert histogram', 'insert mean', 'insert stdev', 'insertCount', 'inserted bases', 'instrument', 'lane', 'library', 'mapped reads', 'mark duplicates', 'mate unmapped reads', 'meanInsert', 'mismatch bases', 'mismatchCount', 'non collapsed bases covered', 'non primary reads', 'number of ends', 'number of targets', 'paired reads', 'pairsMappedAbnormallyFar', 'pairsMappedToDifferentChr', 'properly paired reads', 'qual cut', 'qual fail reads', 'read 1 aligned by cycle', 'read 1 average length', 'read 1 deletion by cycle', 'read 1 hard clip by cycle', 'read 1 insertion by cycle', 'read 1 length histogram', 'read 1 mismatch by cycle', 'read 1 quality by cycle', 'read 1 quality histogram', 'read 1 soft clip by cycle', 'read 2 aligned by cycle', 'read 2 average length', 'read 2 deletion by cycle', 'read 2 hard clip by cycle', 'read 2 insertion by cycle', 'read 2 length histogram', 'read 2 mismatch by cycle', 'read 2 quality by cycle', 'read 2 quality histogram', 'read 2 soft clip by cycle', 'read ? aligned by cycle', 'read ? average length', 'read ? deletion by cycle', 'read ? hard clip by cycle', 'read ? insertion by cycle', 'read ? length histogram', 'read ? mismatch by cycle', 'read ? quality by cycle', 'read ? quality histogram', 'read ? soft clip by cycle', 'reads on target', 'reads per start point', 'readsMissingMDtags', 'run name', 'sample', 'soft clip bases', 'softClipCount', 'stdevInsert', 'target file', 'target size', 'total reads', 'unmapped reads'])
    NEW_KEYS = set(['average read length', 'barcode', 'bases mapped', 'deleted bases', 'hard clip bases', 'insert max', 'insert size average', 'insert size histogram', 'insert size standard deviation', 'inserted bases', 'instrument', 'lane', 'library', 'mapped reads', 'mark duplicates', 'mismatched bases', 'non primary reads', 'number of targets', 'paired end', 'paired reads', 'pairsMappedAbnormallyFar', 'pairsMappedToDifferentChr', 'properly paired reads', 'qual cut', 'qual fail reads', 'read 1 aligned by cycle', 'read 1 average length', 'read 1 deletion by cycle', 'read 1 hard clip by cycle', 'read 1 insertion by cycle', 'read 1 length histogram', 'read 1 mismatch by cycle', 'read 1 quality by cycle', 'read 1 quality histogram', 'read 1 soft clip by cycle', 'read 2 aligned by cycle', 'read 2 average length', 'read 2 deletion by cycle', 'read 2 hard clip by cycle', 'read 2 insertion by cycle', 'read 2 length histogram', 'read 2 mismatch by cycle', 'read 2 quality by cycle', 'read 2 quality histogram', 'read 2 soft clip by cycle', 'read ? aligned by cycle', 'read ? average length', 'read ? deletion by cycle', 'read ? hard clip by cycle', 'read ? insertion by cycle', 'read ? length histogram', 'read ? mismatch by cycle', 'read ? quality by cycle', 'read ? quality histogram', 'read ? soft clip by cycle', 'reads mapped and paired', 'reads on target', 'readsMissingMDtags', 'run name', 'sample', 'sample rate', 'soft clip bases', 'target file', 'total reads', 'total target size', 'unmapped reads'])

    OBSOLETE_KEYS = set([
        'alignedCount',
        'collapsed bases covered',
        'deletionCount',
        'hardClipCount',
        'insertCount',
        'meanInsert',
        'mismatchCount',
        'non collapsed bases covered',
        'reads per start point',
        'softClipCount',
        'stdevInsert',
        'workflow version',
        'workflow name',
        'external name',
        'mate unmaped reads',
        'sample group',
        'group id',
        'last modified',
        'sequencing type'
    ])
    
    def __init__(self, old_path, new_path):
        with open(old_path) as f: self.old_data = json.loads(f.read())
        with open(new_path) as f: self.new_data = json.loads(f.read())

        self.old_to_new_keys = {
            'aligned bases': 'bases mapped',
            'insert histogram': 'insert size histogram',
            'insert mean': 'insert size average',
            'insert stdev': 'insert size standard deviation',
            'mate unmapped reads': 'reads mapped and paired',
            'mate unmaped reads': 'reads mapped and paired', # typo in Perl script
            'mismatch bases': 'mismatched bases',
            'number of ends': 'paired end',
            'target size': 'total target size'
        }
        for key in self.OLD_KEYS:
            if key not in self.old_to_new_keys:
                self.old_to_new_keys[key] = key


        #print(sorted(self.old_data.keys()))
        #print(sorted(self.new_data.keys()))
        

    def run(self):        
        if not(set(self.old_data.keys()) == self.OLD_KEYS):
            sys.stderr.write("Unexpected key set in old output\n")
            print("Added:", set(self.old_data.keys()) - self.OLD_KEYS, file=sys.stderr)
            print("Removed:", self.OLD_KEYS - set(self.old_data.keys()), file=sys.stderr)

        if not(set(self.new_data.keys()) == self.NEW_KEYS):
            sys.stderr.write("Unexpected key set in new output\n")
            print("Added:", set(self.new_data.keys()) - self.NEW_KEYS, file=sys.stderr)
            print("Removed:", self.NEW_KEYS - set(self.new_data.keys()), file=sys.stderr)

        mismatch = 0
        total = 0
        for old_key in self.old_data.keys():
            if old_key in self.OBSOLETE_KEYS:
                continue
            total += 1
            new_key = self.old_to_new_keys[old_key]
            old_output = self.old_data[old_key]
            new_output = self.new_data[new_key]
            if self.old_data[old_key] != self.new_data[new_key]:
                # old data records some numbers as strings; ignore
                if str(self.old_data[old_key]) == str(self.new_data[new_key]):
                    continue
                elif type(self.old_data[old_key]) is dict:
                    if (self.old_data=={} or self.isAllZero(self.old_data)) and self.isAllZero(self.new_data):
                        continue
                    mismatch_old_to_new = 0
                    mismatch_new_to_old = 0
                    old_subdict = self.old_data[old_key]
                    new_subdict = self.new_data[new_key]
                    old_example = []
                    new_example = []
                    max_examples = 4
                    for sub_key in sorted(old_subdict.keys()):
                        if old_subdict.get(sub_key) != new_subdict.get(sub_key):
                            mismatch_old_to_new += 1
                            if len(old_example) < max_examples:
                                old_example.append([sub_key, old_subdict.get(sub_key), new_subdict.get(sub_key)])
                    for sub_key in sorted(new_subdict.keys()):
                        if new_subdict.get(sub_key) != old_subdict.get(sub_key):
                            mismatch_new_to_old += 1
                            if len(new_example) < max_examples:
                                new_example.append([sub_key, old_subdict.get(sub_key), new_subdict.get(sub_key)])
                    old_output = "%i of %i fields in old data do not match\nExamples: %s" % (mismatch_old_to_new, len(old_subdict), str(old_example))
                    new_output = "%i of %i fields in new data do not match\nExamples: %s" % (mismatch_new_to_old, len(new_subdict), str(new_example))
                print("###\nMismatched data for keys (%s, %s)" % (old_key, new_key))
                print(old_output)
                print(new_output)
                mismatch += 1
        print("Mismatches found for %i of %i top-level fields" % (mismatch, total))

    def isAllZero(self, dictionary):
        allZero = True
        for val in dictionary.values():
            if val != 0:
                allZero = False
                break
        return allZero
        
def main():
    parser = argparse.ArgumentParser(description='Diff for old and new BAMQC JSON files.')
    parser.add_argument('-n', '--new', metavar='PATH', required=True,
                        help='Path to new-format input JSON file. Required.')
    parser.add_argument('-o', '--old', metavar='PATH', required=True,
                        help='Path to old-format input JSON file. Required.')
    args = parser.parse_args()
    diff(args.old, args.new).run()
    
if __name__ == "__main__":
    main()
