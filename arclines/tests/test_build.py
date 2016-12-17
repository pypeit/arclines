# Module to run tests on I/O


import numpy as np
import os
import pytest

from astropy.table import Table

from arclines import build_lists

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_init_line_list():
    init_tbl = build_lists.init_line_list()
    #
    assert len(init_tbl) == 1

