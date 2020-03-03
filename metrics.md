# Metric specifications for BamQC

BamQC metrics are computed using a number of methods, including third-party software tools, and output in JSON format. The output also contains metadata, such as the instrument and lane names; and input parameters, such as the target file.

This document:
- Summarizes how metrics are computed
- Describes read filtering options
- Details some special cases of interest

Some metrics have been renamed in the move from the old Perl to new Python implementation.

## Summary of fields

- The 'DS' column indicates which metrics are affected by downsampling (if any).
- The 'WFM' column indicates which fields are in the subset output by `write_fast_metrics.py`.

| Name                            | Former name          | Source                  | DS | WFM |
| --------------------------------|----------------------|-------------------------|----|-----|
| alignment reference             |                      | input parameter         | N  |  Y  |
| average read length             |                      | samtools stats: RL      | N  |  Y  |
| bases per target                |                      | bedtools                | Y  |  N  |
| barcode                         |                      | metadata                | N  |  N  |
| bases mapped                    | aligned bases        | samtools stats: SN      | N  |  Y  |
| coverage per target             |                      | bedtools                | Y  |  N  |
| coverage histogram              |                      | bedtools                | Y  |  N  |
| deleted bases                   |                      | samtools stats: ID      | N  |  Y  |
| donor                           |                      | metadata                | N  |  N  |
| group id                        |                      | metadata                | N  |  N  |
| hard clip bases                 |                      | CIGAR                   | Y  |  N  |
| insert max                      |                      | input parameter         | N  |  Y  |
| insert size average             | insert mean          | samtools stats: SN      | N  |  Y  |
| insert size histogram           | insert histogram     | samtools stats: IS      | N  |  Y  |
| insert size standard deviation  | insert stdev         | samtools stats: SN      | N  |  Y  |
| inserted bases                  |                      | samtools stats: ID      | N  |  Y  |
| instrument                      |                      | metadata                | N  |  N  |
| lane                            |                      | metadata                | N  |  N  |
| library                         |                      | metadata                | N  |  N  |
| library design                  |                      | metadata                | N  |  N  |
| low-quality reads meta          |                      | metadata                | N  |  N  |
| mapped reads                    |                      | samtools stats: SN      | N  |  Y  |
| mark duplicates                 |                      | Picard MarkDuplicates   | N  |  N  |
| mismatched bases                | mismatch bases       | samtools stats: SN      | N  |  Y  |
| non primary reads               |                      | samtools stats: SN      | N  |  Y  |
| non-primary reads meta          |                      | metadata                | N  |  N  |
| number of targets               |                      | bedtools                | N  |  N  |
| package version                 |                      | package version file    | N  |  Y  |
| paired end                      | number of ends       | samtools stats: SN      | N  |  Y  |
| paired reads                    |                      | samtools stats: SN      | N  |  Y  |
| pairsMappedAbnormallyFar        |                      | insert size histogram   | Y  |  Y  |
| pairsMappedToDifferentChr       |                      | samtools stats: SN      | N  |  Y  |
| properly paired reads           |                      | samtools stats: SN      | N  |  Y  |
| qual cut                        |                      | input parameter         | N  |  N  |
| qual fail reads                 |                      | samtools view           | N  |  N  |
| read 1 aligned by cycle         |                      | CIGAR                   | Y  |  N  |
| read 1 average length           |                      | samtools stats: FRL     | N  |  Y  |
| read 1 deletion by cycle        |                      | CIGAR                   | Y  |  N  |
| read 1 hard clip by cycle       |                      | CIGAR                   | Y  |  N  |
| read 1 insertion by cycle       |                      | CIGAR                   | Y  |  N  |
| read 1 length histogram         |                      | samtools stats: FRL     | N  |  Y  |
| read 1 mismatch by cycle        |                      | samtools stats: MPC     | N  |  Y  |
| read 1 quality by cycle         |                      | samtools stats: FFQ     | N  |  Y  |
| read 1 quality histogram        |                      | samtools stats: FFQ     | N  |  Y  |
| read 1 soft clip by cycle       |                      | CIGAR                   | Y  |  N  |
| read 2 aligned by cycle         |                      | CIGAR                   | Y  |  N  |
| read 2 average length           |                      | samtools stats: LRL     | N  |  Y  |
| read 2 deletion by cycle        |                      | CIGAR                   | Y  |  N  |
| read 2 hard clip by cycle       |                      | CIGAR                   | Y  |  N  |
| read 2 insertion by cycle       |                      | CIGAR                   | Y  |  N  |
| read 2 length histogram         |                      | samtools stats: LRL     | N  |  Y  |
| read 2 mismatch by cycle        |                      | samtools stats: MPC     | N  |  Y  |
| read 2 quality by cycle         |                      | samtools stats: LFQ     | N  |  Y  |
| read 2 quality histogram        |                      | samtools stats: LFQ     | N  |  Y  |
| read 2 soft clip by cycle       |                      | CIGAR                   | Y  |  N  |
| read ? aligned by cycle         |                      | CIGAR                   | Y  |  N  |
| read ? average length           |                      | CIGAR                   | Y  |  N  |
| read ? deletion by cycle        |                      | CIGAR                   | Y  |  N  |
| read ? hard clip by cycle       |                      | CIGAR                   | Y  |  N  |
| read ? insertion by cycle       |                      | CIGAR                   | Y  |  N  |
| read ? length histogram         |                      | CIGAR                   | Y  |  N  |
| read ? mismatch by cycle        |                      | samtools stats: MPC     | N  |  Y  |
| read ? quality by cycle         |                      | CIGAR                   | Y  |  N  |
| read ? quality histogram        |                      | CIGAR                   | Y  |  N  |
| read ? soft clip by cycle       |                      | CIGAR                   | Y  |  N  |
| reads mapped and paired         | mate unmapped reads  | samtools stats: SN      | N  |  Y  |
| reads on target                 |                      | bedtools                | Y  |  N  |
| RPSP (DEPRECATED)               |                      | RNAME/POS fields        | Y  |  N  |
| readsMissingMDtags              |                      | CIGAR                   | Y  |  N  |
| run name                        |                      | metadata                | N  |  N  |
| sample                          |                      | metadata                | N  |  N  |
| sample level                    |                      | input parameter         | N  |  N  |
| sample total                    |                      | samtools view -s        | N  |  N  |
| soft clip bases                 |                      | CIGAR                   | Y  |  N  |
| target file                     |                      | input parameter         | N  |  N  |
| target sizes                    |                      | bedtools                | Y  |  N  |
| tissue origin                   |                      | metadata                | N  |  N  |
| tissue type                     |                      | metadata                | N  |  N  |
| total bases on target           |                      | bedtools                | Y  |  N  |
| total input reads meta          |                      | metadata                | N  |  N  |
| total reads                     |                      | samtools stats: SN      | N  |  Y  |
| total target size               | target size          | bedtools                | N  |  N  |
| unmapped reads                  |                      | samtools stats: SN      | N  |  Y  |
| unmapped reads meta             |                      | metadata                | N  |  N  |
| workflow version                |                      | input parameter         | N  |  N  |

## General notes on data sources

### samtools and bedtools

- `samtools` and `bedtools` are wrapped by the Python packages `pysam` and `pybedtools`, respectively.
- Exception: For coverage metrics, `bedtools` is called in a subprocess to circumvent a bug in `pybedtools`.)
- `pysam` is also used for processing CIGAR strings.

### pairsMappedAbnormallyFar

`pairsMappedAbnormallyFar` is a secondary metric, computed from the insert size histogram, which in turn is derived from `samtools stats`. In this instance, 'abnormal' is defined by the `insert max` parameter.

### Read counts: Meta and Otherwise

The metrics `non-primary reads meta`, `low-quality reads meta`, `total input reads meta`, and `unmapped reads meta` are supplied as metadata, to report filtering carried out by upstream workflow tasks. They respectively correspond to `non primary reads`, `qual fail reads`, `total reads`, and `unmapped reads`, which are generated by the `bam-qc-metrics` package using `samtools`.

The latter group of metrics are computed on the BAM file input to `bam-qc-metrics` (before downsampling, if any). For example, we might filter out unmapped reads at a previous workflow step and report the total in `unmapped reads meta`. The filtered BAM file input to `bam-qc-metrics` has no unmapped reads, so the `unmapped reads` metric will be zero. Similar considerations apply to the other metrics.

The read counts computed by `bam-qc-metrics` will be redundant, and potentially confusing, when BAM QC is fully migrated to a workflow which does filtering upstream. Therefore, they [will be removed](https://github.com/oicr-gsi/bam-qc-metrics/issues/33) in a future release.

## Metadata

The metadata fields are supplied as JSON input to the `run_bam_qc.py` script.

Two sets of metadata fields are in use:
- `['barcode', 'instrument', 'lane', 'library', 'run name', 'sample']` for lane-level BAM input
- `['donor', 'group id', 'library design', 'tissue origin', 'tissue type']` for merged BAM input

Metadata fields not appropriate for the given BAM input will appear in JSON output with a null value.

## Filtering

There are two possible filter mechanisms in BamQC: _downsampling_ and _quality filtering_.

If both are in effect, quality filtering is applied first.

Neither filter is supported by the `write_fast_metrics.py` script.

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
