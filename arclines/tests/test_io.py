# Module to run tests on I/O


import numpy as np
import os
import pytest

from astropy.table import Table

from arclines import io as arcl_io

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_load_sources():
    sources = arcl_io.load_sources()
    # Test
    assert isinstance(sources, Table)

