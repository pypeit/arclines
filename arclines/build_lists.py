""" Module for building arcline lists
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import pdb

from collections import OrderedDict

from astropy.table import Table

import arclines  # For path
from arclines import io as arcl_io
from arclines import load_source


def init_line_list(len_src=30):
    """ Initialize a Table for a linelist
    Rigidly enforces table column formats
    Strings are the most annoying

    Returns
    -------
    init_tbl : Table
      One dummy row
    """
    # Load sources to check
    sources = arcl_io.load_source_table()
    src_files = sources['File'].data
    if src_files.dtype.itemsize > len_src:
        raise ValueError("Sources now exceeds table.  Should fix source")
    dummy_src = str('#')*len_src
    # Ion in Lamp
    dummy_ion = str('#')*5
    #

    # Dict for Table
    idict = OrderedDict()
    idict['ion']=dummy_ion
    idict['wave'] = 0.
    idict['NIST']=1
    idict['amplitude']=1.
    idict['Source']=dummy_src

    # Table
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    init_tbl = Table(lst, names=tkeys)

    # Return
    return init_tbl


def master_build(check_only=True):
    """ Master loop to build the line lists

    Parameters
    ----------
    check_only : bool, optional
      Only check to see what would be created, i.e. do not execute
    """
    # Load sources
    sources = arcl_io.load_source_table()

    # Loop on sources
    for source in sources:
        # Load line table
        if source['Format'] == 'PYPIT1':
            load_source.pypit(1,)
        else:
            raise IOError("Format {:s} for source {:s} is not supporrted".format(
                source['Format'], source['File']))
