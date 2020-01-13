
from .bam_qc import bam_qc, fast_metric_finder, fast_metric_writer, validator, version_updater

# 'base' and 'fast_metric_finder' classes are exported for tests
# 'version_updater' class is exported for update_test_data_version.py

import os

def read_package_version():
    """
    This method depends on relative path to the VERSION file
    So it has been placed in the package __init__ file, whose location will never change

    VERSION file is in 'etc/versions/bam-qc-metrics'
    'etc' directory may be in one of two places relative to __init__.py:
    - Parent directory, if testing the source code
    - 4 directories up, following install, eg:
      - bam-qc-metrics-0.1.8/etc/versions/bam_qc_metrics/VERSION
      - bam-qc-metrics-0.1.8/lib/python3.6/site-packages/bam_qc_metrics/__init__.py
    """
    in_path = None
    ver_path = os.path.join('etc', 'versions', 'bam_qc_metrics', 'VERSION')
    test_path = os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir, ver_path))
    install_path = os.path.realpath(os.path.join(os.path.dirname(__file__), *[os.pardir]*4, ver_path))
    if os.path.exists(test_path):
        in_path = test_path
    elif os.path.exists(install_path):
        in_path = install_path
    else:
        raise FileNotFoundError("Cannot find VERSION file; bad installation?")
    with open(in_path) as version_file:
        package_version = version_file.read().strip()
    return package_version
