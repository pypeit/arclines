""" Module for I/O in arclines
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import datetime

from astropy.table import Table


def load_line_list(line_file, add_path=False):
    """
    Parameters
    ----------
    line_file
    add_path

    Returns
    -------
    line_list : Table

    """
    line_list = Table.read(line_file, format='ascii.fixed_width', comment='#')
    # Return
    return line_list

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
