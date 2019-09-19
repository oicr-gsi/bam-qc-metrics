"""
Setup script for bam-qc-metrics
"""

import os
from setuptools import setup

package_identifier = "bam_qc_metrics"
version_dir = os.path.join('etc', 'versions', package_identifier)
version_filename = 'VERSION'
with open(os.path.join(os.path.dirname(__file__), version_dir, version_filename)) as version_file:
    package_version = version_file.read().strip()

setup(
    name='bam-qc-metrics',
    version=package_version,
    scripts=['bin/run_bam_qc.py', ],
    packages=[package_identifier],
    install_requires=['attrs', 'jsonschema', 'pybedtools', 'pyrsistent', 'pysam', 'six'],
    data_files=[(version_dir, [os.path.join(version_dir, version_filename)])],
    python_requires='>=3.5',
    author="Iain Bancarz",
    author_email="ibancarz@oicr.on.ca",
    description="BAM QC metrics",
    long_description="Python implementation of BAM QC metrics in use at the Ontario Institute for Cancer Research, https://oicr.on.ca",
    url="https://github.com/oicr-gsi/bam-qc-metrics",
)
