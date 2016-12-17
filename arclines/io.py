""" Module for I/O in arclines
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import datetime
import pdb

from astropy.table import Table, vstack

from arclines import defs

def load_line_list(line_file, add_path=False):
    """ Simple load of a line list

    Parameters
    ----------
    line_file : str
    add_path : bool, optional
      Not yet implemented

    Returns
    -------
    line_list : Table

    """
    line_list = Table.read(line_file, format='ascii.fixed_width', comment='#')
    # Return
    return line_list


def load_line_lists(lines, unknown=False):
    """ Loads a series of line list files

    Parameters
    ----------
    lamps : list
    unknown : bool, optional

    Returns
    -------
    line_list : Table

    """
    import arclines # For path
    line_path = arclines.__path__[0]+'/data/lists/'

    # Read standard files
    lists = []
    for line in lines:
        line_file = line_path+'{:s}_lines.dat'.format(line)
        if not os.path.isfile(line_file):
            raise IOError("Input line {:s} is not included in arclines".format(line))
        else:
            lists.append(load_line_list(line_file))
    # Stack
    line_lists = vstack(lists, join_type='exact')

    # Unknown
    if unknown:
        unkn_lines = load_unknown_list(lines)
        unkn_lines.remove_column('line_flag')  # may wish to have this info
        # Stack
        line_lists = vstack([line_lists, unkn_lines])

    # Return
    return line_lists


def load_source_table():
    """ Load table of arcline sources

    Returns
    -------
    sources : Table

    """
    import arclines # For path
    src_file = arclines.__path__[0]+'/data/sources/arcline_sources.ascii'
    # Load
    sources = Table.read(src_file, format='ascii.fixed_width', comment='#')
    # Return
    return sources


def load_unknown_list(lines, unknwn_file=None):
    """
    Parameters
    ----------
    lines : list
    unknwn_file : str, optional

    Returns
    -------
    unknwn_lines : Table

    """
    import arclines # For path
    line_dict = defs.lines()
    # Load
    line_path = arclines.__path__[0]+'/data/lists/'
    if unknwn_file is None:
        unknwn_file = line_path+'UNKNWN_lines.dat'
    # Cut on input lamps
    line_list = load_line_list(unknwn_file)
    msk = np.array([False]*len(line_list))
    for line in lines:
        line_flag = line_dict[line]
        match = line_list['line_flag'] % 2*line_flag >= line_flag
        msk[match] = True
        pdb.set_trace()
    # Finish
    return line_list[msk]


def write_line_list(tbl, outfile):
    """
    Parameters
    ----------
    tbl
    outfile
    """
    # Format
    tbl['wave'].format = '10.4f'
    # Write
    with open(outfile,'w') as f:
        f.write('#Creation Date: {:s}\n'.format(str(datetime.date.today().strftime('%Y-%b-%d'))))
        tbl.write(f, format='ascii.fixed_width')
