# Module to run tests on I/O


import numpy as np
import os
import pytest

from astropy.table import Table

from arclines import io as arcl_io

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_load_line_lists():
    line_lists = arcl_io.load_line_lists(['HgI','ZnI'])
    pytest.set_trace()

def test_load_line_list():
    import arclines
    llst = arcl_io.load_line_list(arclines.__path__[0]+'/data/lists/HgI_lines.dat')


def test_load_sources():
    sources = arcl_io.load_source_table()
    # Test
    assert isinstance(sources, Table)

