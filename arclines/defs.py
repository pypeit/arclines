""" Module for definitions
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from astropy.table import Table


def instruments():
    """ Dict to convert instrument to bitwise flag

    Returns
    -------
    instr_dict : Table

    """
    instr_dict = {}
    #
    instr_dict['LRISr'] = 2**0
    #
    return instr_dict

