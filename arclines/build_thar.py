#!/usr/bin/env python
"""
Generate a ThAr linelist based on the NIST data and Michael Murphy's list.
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str


import arclines  # for path
line_path = arclines.__path__[0] + '/data/lists/'
nist_path = arclines.__path__[0] + '/data/NIST/'


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Build the arclines line lists from scratch')
    parser.add_argument("-w", "--write", default=False, action='store_true', help="Actually write files?")
    parser.add_argument("--skip_stop", default=False, action='store_true', help="Skip warning stop?")
    parser.add_argument("--unknowns", default=False, action='store_true', help="Create Unknown list?")
    parser.add_argument("--plots", default=False, action='store_true', help="Create plots?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    """ This script converts all Th and Ar lines in NIST
    to an astropy Table file that can be read by PYPIT

    Parameters

    ----------
    args

    Returns
    -------

    """
    # Grab arguments
    pargs = parser(options=args)

    import numpy as np
    from arclines import build_lists
    from arclines import arcio as arcl_io
    from arclines import load_source
    from arclines import plots as arcl_plots
    from arclines import utils as arcl_utils


    print("=============================================================")
    print("This script is for EXPERTS ONLY")
    print("Continue only if you know what you are doing")
    print("Otherwise exit")
    print("p.s.  You need to remove the files you wish to re-build")
    print("=============================================================")
    if not pargs.skip_stop:
        pdb.set_trace()

    # Load sources
    sources = arcl_io.load_source_table()






    # Load the NIST ThAr list
    llist = arcl_io.load_line_lists(["ThAr"], NIST=True)
    pdb.set_trace()

    linelist = build_lists.init_line_list()

    # Delete this, it's just for reference
    # idict = OrderedDict()
    # idict['ion'] = dummy_line
    # idict['wave'] = 0.
    # idict['NIST'] = 0
    # idict['Instr'] = 0  # Flag for instrument
    # idict['amplitude'] = 0
    # idict['Source'] = dummy_src

    # Add all lines, and flag if they meet the quality described by
    # Murphy et al. (2007), MNRAS, 378, 221
    # http://www.astronomy.swin.edu.au/~mmurphy/thar/thar.html



    # Write?
    if not pargs.write:
        print("=============================================================")
        print("Rerun with --write if you are happy with what you see.")
        print("=============================================================")
        return

    # Write the table to disk
    thar_file = line_path + '{:s}_lines.dat'.format("ThAr")
    arcl_io.write_line_list(linelist, thar_file)
    return

if __name__ == '__main__':
    main()
