
# bam-qc-metrics

Quality control metrics for BAM alignment files.

Process a BAM file and write metrics and metadata in JSON format.

## Prerequisites

The `bedtools` executable must be available on the current `PATH`. Note that version 2.28.0 of the statically compiled binary on Github has an inconsistent version number, which means it does not work with `pybedtools`. Instead it is necessary to download the tarball and compile the `bedtools` binary locally.

## Usage

There are two command-line scripts:
- `bin/run_bam_qc.py` to run all metrics
- `bin/write_fast_metrics.py` to run a faster subset of metrics, eg. for RNASeqQC

Run either script with `--help` for instructions.

## Tests

Run `test/test.py` for Python tests.

## Release procedure

Before a release, the package version number in `./etc/versions/bam_qc_metrics/VERSION` must
be updated.

The package version number also appears in JSON output. The script
`bin/update_test_data_version.py` updates the expected test data. Any new test output files
must be added to the list for the update script, which is located in class
`version_updater` in `bam-qc-metrics/bam_qc.py`.

## Metric notes

See `metrics.md` for a detailed account of the metrics in use.

## Conventions

- Changelog: See CHANGELOG.md and https://keepachangelog.com/en/1.0.0/
- Semantic Versioning: See https://semver.org/
