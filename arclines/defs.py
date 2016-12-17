""" Module for definitions
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from astropy.table import Table


def instruments():
    """ Dict to convert instrument to bitwise flag
    WARNING: Modifying any of the following is a *bad* idea
      Adding is ok

    Returns
    -------
    instr_dict : Table

    """
    instr_dict = {}
    #
    instr_dict['LRISr'] = 2**0
    #
    return instr_dict

def lamps():
    """ Dict of lamps included in this database
    WARNING: Modifying any of the following is a *bad* idea
      Adding is ok

    Returns
    -------
    lamp_dict : dict

    """
    lamp_dict = {}
    #
    lamp_dict['ArI'] = 2**0
    lamp_dict['HgI'] = 2**1
    lamp_dict['KrI'] = 2**2
    lamp_dict['NeI'] = 2**3
    lamp_dict['XeI'] = 2**4
    #
    return lamp_dict


def str_len():
    """ Hard-codes length of strings in the database
    WARNING: Modifying any of the following is a *bad* idea

    Returns
    -------
    strlen_dict : dict

    """
    strlen_dict = {}
    # Length of ion name
    strlen_dict['ion'] = 6
    # Length of data file name for line source
    strlen_dict['Source'] = 30
    # Return
    return strlen_dict
