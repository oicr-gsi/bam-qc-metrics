"""
Setup script for bam-qc-metrics
"""

import os
from setuptools import setup, find_packages

package_identifier = "bam_qc_metrics"
version_dir = os.path.join('etc', 'versions', package_identifier)
package_identifier_lite = "bam_qc_metrics_lite"
version_dir_lite = os.path.join('etc', 'versions', package_identifier_lite)
version_filename = 'VERSION'
with open(os.path.join(os.path.dirname(__file__), version_dir_lite, version_filename)) as version_file:
    package_version = version_file.read().strip()

setup(
    name='bam-qc-metrics',
    version=package_version,
    scripts=['bin/run_bam_qc.py', 'bin/write_fast_metrics.py', 'bin/bam_qc_merger.py', 'bin/run_bam_qc_lite.py'],
    packages=find_packages(),
    install_requires=['attrs', 'jsonschema', 'pybedtools', 'pyrsistent', 'pysam', 'six'],
    data_files=[(version_dir, [os.path.join(version_dir, version_filename)]),(version_dir_lite, [os.path.join(version_dir_lite, version_filename)])],
    python_requires='>=3.12',
    author="Iain Bancarz",
    author_email="ibancarz@oicr.on.ca",
    description="BAM QC metrics",
    long_description="Python implementation of BAM QC metrics in use at the Ontario Institute for Cancer Research, https://oicr.on.ca",
    url="https://github.com/oicr-gsi/bam-qc-metrics",
)
