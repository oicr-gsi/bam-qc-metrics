
# bam-qc-metrics

Quality control metrics for BAM alignment files.

Process a BAM file and write metrics and metadata in JSON format.

## Prerequisites

The `bamtools` executable must be available on the current `PATH`. Note that version 2.28.0 of the statically compiled binary on Github has an inconsistent version number, which means it does not work with `pybedtools`. Instead it is necessary to download the tarball and compile the `bamtools` binary locally.

## Usage

Run the script `bam_qc_metrics/bam_qc.py` with `--help` for instructions.

## Conventions

- Changelog: See CHANGELOG.md and https://keepachangelog.com/en/1.0.0/
- Semantic Versioning: See https://semver.org/
