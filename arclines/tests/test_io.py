# Module to run tests on I/O


import numpy as np
import os
import pytest

from astropy.table import Table

import arclines
from arclines import io as arcl_io


#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def test_load_line_lists():
    line_lists = arcl_io.load_line_lists(['HgI','ZnI'])
    # Unknown
    line_lists = arcl_io.load_line_lists(['HgI','ZnI'], unknown=True)


def test_load_line_list():
    import arclines
    llst = arcl_io.load_line_list(arclines.__path__[0]+'/data/lists/HgI_lines.dat')


def test_load_sources():
    sources = arcl_io.load_source_table()
    # Test
    assert isinstance(sources, Table)

def test_load_nist():
    from astropy.table import vstack

    NIST = True
    nist_path = arclines.__path__[0]+'/data/NIST/'

    lines = ['ArI', 'CdI', 'CuI', 'HgI', 'KrI', 'NeI', 'XeI', 'ZnI']
    lists = []
    for line in lines:
        line_file = nist_path+'{:s}_vacuum.ascii'.format(line)
        tbl = arclines.io.load_line_list(line_file, NIST=NIST)
        lists.append(tbl)
    # Stack
    line_lists = vstack(lists, join_type='exact')
    # Test
    for key in ['wave','Aki','Rel.','Ion','NIST']:
        assert key in line_lists.keys()

