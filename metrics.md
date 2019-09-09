# Metric specifications for BamQC

BamQC metrics are computed using a number of methods, including third-party software tools, and output in JSON format. The output also contains metadata, such as the instrument and lane names; and input parameters, such as the target file.

This document:
- Summarizes how metrics are computed
- Describes read filtering options
- Details some special cases of interest

Some metrics have been renamed in the move from the old Perl to new Python implementation.

## Summary of fields

| Name                            | Former name             | Source                  | DS |
| --------------------------------|-------------------------|-------------------------|----|
| alignment reference             |                         | input parameter         | N  |
| average read length             |                         | samtools stats: RL      | N  |
| barcode                         |                         | metadata                | N  |
| bases mapped                    | aligned bases           | samtools stats: SN      | N  |
| coverage per target             |                         | bedtools                | Y  |
| coverage histogram              |                         | bedtools                | Y  |
| deleted bases                   |                         | samtools stats: ID      | N  |
| hard clip bases                 |                         | CIGAR                   | Y  |
| insert max                      |                         | input parameter         | N  |
| insert size average             | insert mean             | samtools stats: SN      | N  |
| insert size histogram           | insert histogram        | samtools stats: IS      | N  |
| insert size standard deviation  | insert stdev            | samtools stats: SN      | N  |
| inserted bases                  |                         | samtools stats: ID      | N  |
| instrument                      |                         | metadata                | N  |
| lane                            |                         | metadata                | N  |
| library                         |                         | metadata                | N  |
| mapped reads                    |                         | samtools stats: SN      | N  |
| mark duplicates                 |                         | Picard MarkDuplicates   | N  |
| mismatched bases                | mismatch bases          | samtools stats: SN      | N  |
| non primary reads               |                         | samtools stats: SN      | N  |
| number of targets               |                         | bedtools                | N  |
| package version                 |                         | package version file    | N  |
| paired end                      | number of ends          | samtools stats: SN      | N  |
| paired reads                    |                         | samtools stats: SN      | N  |
| pairsMappedAbnormallyFar        |                         | insert size histogram   | Y  |
| pairsMappedToDifferentChr       |                         | samtools stats: SN      | N  |
| properly paired reads           |                         | samtools stats: SN      | N  |
| qual cut                        |                         | input parameter         | N  |
| qual fail reads                 |                         | samtools view           | N  |
| read 1 aligned by cycle         |                         | CIGAR                   | Y  |
| read 1 average length           |                         | samtools stats: FRL     | N  |
| read 1 deletion by cycle        |                         | CIGAR                   | Y  |
| read 1 hard clip by cycle       |                         | CIGAR                   | Y  |
| read 1 insertion by cycle       |                         | CIGAR                   | Y  |
| read 1 length histogram         |                         | samtools stats: FRL     | N  |
| read 1 mismatch by cycle        |                         | samtools stats: MPC     | N  |
| read 1 quality by cycle         |                         | samtools stats: FFQ     | N  |
| read 1 quality histogram        |                         | samtools stats: FFQ     | N  |
| read 1 soft clip by cycle       |                         | CIGAR                   | Y  |
| read 2 aligned by cycle         |                         | CIGAR                   | Y  |
| read 2 average length           |                         | samtools stats: LRL     | N  |
| read 2 deletion by cycle        |                         | CIGAR                   | Y  |
| read 2 hard clip by cycle       |                         | CIGAR                   | Y  |
| read 2 insertion by cycle       |                         | CIGAR                   | Y  |
| read 2 length histogram         |                         | samtools stats: LRL     | N  |
| read 2 mismatch by cycle        |                         | samtools stats: MPC     | N  |
| read 2 quality by cycle         |                         | samtools stats: LFQ     | N  |
| read 2 quality histogram        |                         | samtools stats: LFQ     | N  |
| read 2 soft clip by cycle       |                         | CIGAR                   | Y  |
| read ? aligned by cycle         |                         | CIGAR                   | Y  |
| read ? average length           |                         | CIGAR                   | Y  |
| read ? deletion by cycle        |                         | CIGAR                   | Y  |
| read ? hard clip by cycle       |                         | CIGAR                   | Y  |
| read ? insertion by cycle       |                         | CIGAR                   | Y  |
| read ? length histogram         |                         | CIGAR                   | Y  |
| read ? mismatch by cycle        |                         | samtools stats: MPC     | N  |
| read ? quality by cycle         |                         | CIGAR                   | Y  |
| read ? quality histogram        |                         | CIGAR                   | Y  |
| read ? soft clip by cycle       |                         | CIGAR                   | Y  |
| reads mapped and paired         | mate unmapped reads     | samtools stats: SN      | N  |
| reads on target                 |                         | bedtools                | Y  |
| reads per start point           |                         | RNAME/POS fields        | Y  |
| readsMissingMDtags              |                         | CIGAR                   | Y  |
| run name                        |                         | metadata                | N  |
| sample                          |                         | metadata                | N  |
| sample level                    |                         | input parameter         | N  |
| sample total                    |                         | samtools view -s        | N  |
| soft clip bases                 |                         | CIGAR                   | Y  |
| target file                     |                         | input parameter         | N  |
| target sizes                    |                         | bedtools                | Y  |
| total coverage                  |                         | bedtools                | Y  |
| total reads                     |                         | samtools stats: SN      | N  |
| total target size               | target size             | bedtools                | N  |
| unmapped reads                  |                         | samtools stats: SN      | N  |
| workflow version                |                         | input parameter         | N  |

- The DS column indicates which metrics are affected by downsampling (if any).

## General notes on data sources

- `samtools` and `bedtools` are wrapped by the Python packages `pysam` and `pybedtools`, respectively. (Exception: For coverage metrics, `bedtools` is called in a subprocess to circumvent a bug in `pybedtools`.)
- `pysam` is also used for processing CIGAR strings.
- `pairsMappedAbnormallyFar` is a secondary metric, computed from the insert size histogram, which in turn is derived from `samtools stats`. In this instance, 'abnormal' is defined by the `insert max` parameter.

## Filtering

There are two possible filter mechanisms in BamQC: _downsampling_ and _quality filtering_.

If both are in effect, quality filtering is applied first.

### Speed considerations

On large BAM files, quality filtering is slow, as it writes a large temporary file containing the reads which pass filters. Omitting quality filtering may result in a substantial speed gain. Downsampling is considerably faster, as the output file written is small.

For example, time taken to process a NovaSeq BAM file (22 GB, 545M reads) was as follows:
- Quality filtering on, downsampling off: 10 hours
- Quality filtering on, downsampling on: 67 minutes
- Quality filtering off, downsampling off: 22 minutes

### Downsampling

Downsampling, in which metrics are computed on a randomly selected subset of reads, may be applied for faster and more efficient data processing.

If downsampling is in effect, metrics derived from CIGAR strings and bedtools are evaluated on the downsampled reads. The estimated value of the metric on the full dataset must be scaled accordingly. For example, if 500 million reads have been downsampled to 1 million, and the reported 'hard clip bases' is 250, there will be approximately 125,000 hard clip bases in the entire BAM file. Metrics derived directly from samtools will always be calculated on the full dataset, without downsampling. The `DS` column in the summary table records which metrics are affected by downsampling.

Downsampling is carried out using `samtools view -s`. This method is fast and preserves read pairs, but does not allow sampling an exact number of reads, as it uses a probability of sampling rather than a fixed sample size. This is a [known issue](https://github.com/samtools/samtools/issues/955) with samtools and may be fixed at a later date. Alternatively, a workaround may be implemented in a subsequent version of `bam-qc-metrics`.

The default is downsampling to approximately 1.1 million reads, if the input BAM contains more than this number; and no downsampling otherwise. The default level was chosen so the downsampled file should contain at least 1 million reads. Downsampling can be deactivated, or a custom rate specified, using command-line options. In JSON output, `sample_level` is the target number of reads to downsample, and `sample_total` is the actual number of reads returned by samtools.

A fixed random seed is used for downsampling, so results will be consistent on repeated runs of `bam-qc-metrics`.

### Quality filtering

BamQC has an option to exclude reads with alignment quality below a given threshold. The filter is implemented using `samtools view -q`.

If quality filtering is in effect, the `unmapped reads` metric is computed only on reads which _failed_ the filter -- because by definition, an unmapped read does not have an alignment quality. All other metrics are computed only on reads which _passed_ the filter.

## Special case: Sequence mismatches

The `mismatched bases` field may not be consistent with the mismatch-by-cycle fields for each read.

The former is derived from `samtools stats SN`; the latter from `samtools stats MPC` with a given alignment reference.

Note also that `samtools stats MPC` reports a number of cycles 1 greater than the actual number in the BAM file. For example, a 101-base read will have an entry for cycle 102 with zero mismatches. This value is retained in JSON output for consistency with the samtools results.
