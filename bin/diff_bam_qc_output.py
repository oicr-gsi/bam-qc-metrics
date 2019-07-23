#! /usr/bin/env python3

"""
Script to diff old and new BAMQC JSON outputs, generated by bamqc.pl and bam_qc.py respectively

JSON data structure has changed; maps old fields to new fields, and disregards obsolete ones
"""

import argparse, json, re, sys

class diff:

    ALL_NEW_KEYS = [
        'average read length',
        'barcode',
        'bases mapped',
        'deleted bases',
        'hard clip bases',
        'insert max',
        'insert size average',
        'insert size histogram',
        'insert size standard deviation',
        'inserted bases',
        'instrument',
        'lane',
        'library',
        'mapped reads',
        'mark duplicates',
        'mismatched bases',
        'non primary reads',
        'number of targets',
        'paired end',
        'paired reads',
        'pairsMappedAbnormallyFar',
        'pairsMappedToDifferentChr',
        'properly paired reads',
        'qual cut',
        'qual fail reads',
        'read 1 aligned by cycle',
        'read 1 average length',
        'read 1 deletion by cycle',
        'read 1 hard clip by cycle',
        'read 1 insertion by cycle',
        'read 1 length histogram',
        'read 1 mismatch by cycle',
        'read 1 quality by cycle',
        'read 1 quality histogram',
        'read 1 soft clip by cycle',
        'read 2 aligned by cycle',
        'read 2 average length',
        'read 2 deletion by cycle',
        'read 2 hard clip by cycle',
        'read 2 insertion by cycle',
        'read 2 length histogram',
        'read 2 mismatch by cycle',
        'read 2 quality by cycle',
        'read 2 quality histogram',
        'read 2 soft clip by cycle',
        'read ? aligned by cycle',
        'read ? average length',
        'read ? deletion by cycle',
        'read ? hard clip by cycle',
        'read ? insertion by cycle',
        'read ? length histogram',
        'read ? mismatch by cycle',
        'read ? quality by cycle',
        'read ? quality histogram',
        'read ? soft clip by cycle',
        'reads mapped and paired',
        'reads on target',
        'readsMissingMDtags',
        'run name',
        'sample',
        'sample rate',
        'soft clip bases',
        'target file',
        'total reads',
        'total target size',
        'unmapped reads'
    ]

    NEW_KEYS_TO_COMPARE = [
        'average read length',
        # 'barcode', # metadata
        'bases mapped',
        'deleted bases',
        'hard clip bases',
        # 'insert max', # new metric
        'insert size average',
        'insert size histogram',
        'insert size standard deviation',
        'inserted bases',
        # 'instrument', # metadata
        # 'lane', # metadata
        # 'library', # metadata
        'mapped reads',
        # 'mark duplicates', # new metric
        'mismatched bases',
        'non primary reads',
        'number of targets',
        # 'paired end', # format change
        'paired reads',
        # 'pairsMappedAbnormallyFar', # new metric
        # 'pairsMappedToDifferentChr', # new metric
        'properly paired reads',
        'qual cut',
        'qual fail reads',
        'read 1 aligned by cycle',
        'read 1 average length',
        'read 1 deletion by cycle',
        'read 1 hard clip by cycle',
        'read 1 insertion by cycle',
        'read 1 length histogram',
        'read 1 mismatch by cycle',
        'read 1 quality by cycle',
        'read 1 quality histogram',
        'read 1 soft clip by cycle',
        'read 2 aligned by cycle',
        'read 2 average length',
        'read 2 deletion by cycle',
        'read 2 hard clip by cycle',
        'read 2 insertion by cycle',
        'read 2 length histogram',
        'read 2 mismatch by cycle',
        'read 2 quality by cycle',
        'read 2 quality histogram',
        'read 2 soft clip by cycle',
        # 'read ? aligned by cycle', # format change for empty fields
        # 'read ? average length',
        # 'read ? deletion by cycle',
        # 'read ? hard clip by cycle',
        # 'read ? insertion by cycle',
        # 'read ? length histogram',
        # 'read ? mismatch by cycle',
        # 'read ? quality by cycle',
        # 'read ? quality histogram',
        # 'read ? soft clip by cycle',
        'reads mapped and paired',
        'reads on target',
        # 'readsMissingMDtags', # not present in old data
        # 'run name', # metadata
        # 'sample', # metadata
        # 'sample rate', # metadata
        'soft clip bases',
        # 'target file', # metadata
        'total reads',
        'total target size',
        'unmapped reads'
    ]

    NEW_KEYS_TO_OLD = {
        'bases mapped': 'aligned bases',
        'insert size histogram': 'insert histogram',
        'insert size average': 'insert mean',
        'insert size standard deviation': 'insert stdev',
        'mismatched bases': 'mismatch bases',
        'reads mapped and paired' : 'mate unmaped reads', # typo in Perl script
        'total target size': 'target size',
    }

    HASHMAP_KEYS = [
        'insert size histogram',
        'read 1 aligned by cycle',
        'read 1 deletion by cycle',
        'read 1 hard clip by cycle',
        'read 1 insertion by cycle',
        'read 1 length histogram',
        'read 1 mismatch by cycle',
        'read 1 quality by cycle',
        'read 1 quality histogram',
        'read 1 soft clip by cycle',
        'read 2 aligned by cycle',
        'read 2 deletion by cycle',
        'read 2 hard clip by cycle',
        'read 2 insertion by cycle',
        'read 2 length histogram',
        'read 2 mismatch by cycle',
        'read 2 quality by cycle',
        'read 2 quality histogram',
        'read 2 soft clip by cycle',
    ]

    def __init__(self, old_path, new_path):
        with open(old_path) as f: self.old_data = json.loads(f.read())
        with open(new_path) as f: self.new_data = json.loads(f.read())

    def find_diff(self):
        # account for string conversion, changed key names
        all_diff = {}
        non_equivalent_diff = {}
        compared = 0
        non_equivalent = 0
        for key in self.NEW_KEYS_TO_COMPARE:
            compared += 1
            new_val = self.new_data[key]
            old_key = self.NEW_KEYS_TO_OLD.get(key, key)
            old_val = self.old_data[old_key]
            diff = None
            equivalent = True
            if key in self.HASHMAP_KEYS:
                diff = {}
                for subkey in new_val.keys():
                    gain = new_val[subkey] - old_val.get(subkey, 0)
                    if gain != 0:
                        equivalent = False
                    diff[subkey] = [
                        new_val[subkey],
                        old_val.get(subkey, 0),
                        gain
                    ]
            elif key in ['average read length', 'insert size average', 'insert size standard deviation', 'mapped reads',
                         'non primary reads', 'number of targets', 'qual cut', 'qual fail reads', 'unmapped reads']:
                gain = new_val - float(old_val)
                if gain != 0:
                    equivalent = False
                diff = [
                    new_val,
                    float(old_val),
                    gain
                ]
            else:
                gain = new_val - old_val
                if gain != 0:
                    equivalent = False
                diff = [
                    new_val,
                    old_val,
                    gain
                ]
            all_diff[key] = diff
            if not equivalent:
                non_equivalent += 1
                non_equivalent_diff[key] = diff
        print("COMPARED:", compared, file=sys.stderr)
        print("NON-EQUIVALENT:", non_equivalent, file=sys.stderr)
        for key in sorted(non_equivalent_diff.keys()):
            print("\t", key, file=sys.stderr)
        return (all_diff, non_equivalent_diff)

    def find_filtered_diff(self):
        # filter out 'solved' issues from the non-equivalent diff structure
        return self.filter_known_diff(self.find_diff()[1])
    
    def filter_known_diff(self, diff):
        # check if difference is 'OK', eg. a known rounding error
        # return the remaining "non-OK" diff
        updated = {}
        for key in diff.keys():
            diff_ok = False
            if key == 'total reads':
                # new 'total reads' excludes quality failures
                [new, old, gain] = diff[key]
                [new_qf, old_qf, gain_qf] = diff['qual fail reads']
                if new + new_qf - old == 0:
                    diff_ok = True
                else:
                    extra_key = 'old definition total reads'
                    updated[extra_key] = [new + new_qf, old, new + new_qf - old]
            elif key == 'hard clip bases':
                # old script didn't detect hard clip bases
                [new, old, gain] = diff[key]
                if old == 0:
                    diff_ok = True
            elif re.search('hard clip by cycle', key):
                # old script didn't detect hard clip bases
                total = 0
                for [new, old, gain] in diff[key].values():
                    total += old
                if total == 0:
                    diff_ok = True
            elif re.search('mismatch by cycle', key):
                # new script doesn't detect mismatches in CIGAR string
                total = 0
                for [new, old, gain] in diff[key].values():
                    total += new
                if total == 0:
                    diff_ok = True
            elif re.search('quality by cycle', key):
                # rounding artifact
                max_gain = 0
                for [new, old, gain] in diff[key].values():
                    if abs(gain) > max_gain: max_gain = abs(gain)
                if max_gain < 1:
                    diff_ok = True
            if not diff_ok:
                updated[key] = diff[key]
        return updated
            

class diff_old:

    # was the diff class
    # kept for now, for convenience referring back

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
    
    def __init__(self, old_path, new_path, brief, tsv_path):
        with open(old_path) as f: self.old_data = json.loads(f.read())
        with open(new_path) as f: self.new_data = json.loads(f.read())
        self.brief = brief
        self.tsv_path = tsv_path
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
        if self.tsv_path:
            tsv = open(self.tsv_path, 'w')
        else:
            tsv = None
        for old_key in sorted(self.old_data.keys()):
            if old_key in self.OBSOLETE_KEYS:
                continue
            total += 1
            new_key = self.old_to_new_keys[old_key]
            old_output = self.old_data[old_key]
            new_output = self.new_data[new_key]
            tsv_out = [str(x) for x in [old_key, new_key, self.old_data[old_key], self.new_data[new_key]]]
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
                    if self.brief:
                        old_example.append('...')
                        new_example.append('...')
                        old_diff = str(old_example)
                        new_diff = str(new_example)
                    else:
                        old_diff = str(old_subdict)
                        new_diff = str(new_subdict)
                    old_output = "%i of %i fields in old data do not match: %s" % (mismatch_old_to_new, len(old_subdict), old_diff)
                    new_output = "%i of %i fields in new data do not match: %s" % (mismatch_new_to_old, len(new_subdict), new_diff)
                print("###\nMismatched data for keys (%s, %s)" % (old_key, new_key))
                print(old_output)
                print(new_output)
                mismatch += 1
                tsv_out.append('Not equal')
            else:
                tsv_out.append('Equal')
            if tsv != None:
                print("\t".join(tsv_out), file=tsv)
        print("Mismatches found for %i of %i top-level fields" % (mismatch, total))
        if tsv != None:
            tsv.close()

    def isAllZero(self, dictionary):
        allZero = True
        for val in dictionary.values():
            if val != 0:
                allZero = False
                break
        return allZero
        
def main():
    parser = argparse.ArgumentParser(description='Diff for old and new BAMQC JSON files.')
    parser.add_argument('-b', '--brief', action="store_true", help="Show only examples instead of the full diff by field")
    parser.add_argument('-n', '--new', metavar='PATH', required=True,
                        help='Path to new-format input JSON file. Required.')
    parser.add_argument('-o', '--old', metavar='PATH', required=True,
                        help='Path to old-format input JSON file. Required.')
    parser.add_argument('-t', '--tsv', metavar='PATH',
                        help='Path for TSV output. Optional.')
    args = parser.parse_args()
    result = diff(args.old, args.new).find_filtered_diff()
    print(json.dumps(result))
    
if __name__ == "__main__":
    main()
