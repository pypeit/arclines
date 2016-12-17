""" Module for loading source files
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from astropy.table import Table


def load_sources():
    """ Load table of arcline sources

    Returns
    -------
    sources : Table

    """
    import arclines # For path
    src_file = arclines.__path__[0]+'/data/sources/arcline_sources.ascii'
    # Load
    sources = Table.read(src_file, format='ascii.fixed_width')
    # Return
    return sources
