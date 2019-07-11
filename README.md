
# bam-qc-metrics

Quality control metrics for BAM alignment files.

Process a BAM file and write metrics and metadata in JSON format.

## Prerequisites

The `bamtools` executable must be available on the current `PATH`. Note that version 2.28.0 of the statically compiled binary on Github has an inconsistent version number, which means it does not work with `pybedtools`. Instead it is necessary to download the tarball and compile the `bamtools` binary locally.

## Usage

Run the script `bin/run_bam_qc.py` with `--help` for instructions.

## Tests

Run `test/test.py` for Python tests.

## Note on sequence mismatches

The `mismatched bases` field may not be consistent with the mismatch-by-cycle fields for each read.

The former is derived from `samtools stats`; the latter from CIGAR strings. The CIGAR operation M for 'alignment match' may represent a sequence match or mismatch. So, the CIGAR string does not necessarily record all mismatches. In the event of inconsistency, `mismatched bases` should be taken as correct.

## Conventions

- Changelog: See CHANGELOG.md and https://keepachangelog.com/en/1.0.0/
- Semantic Versioning: See https://semver.org/
