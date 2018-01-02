#!/usr/bin/env python
"""
Generate a ThAr linelist based on the NIST data and Michael Murphy's list.
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

from collections import OrderedDict
from astropy.table import Table, vstack
import sys
import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str


import arclines  # for path
line_path = arclines.__path__[0] + '/data/lists/'
src_path  = arclines.__path__[0] + '/data/sources/'
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

    # Set the tolerance for matching lines (in A)
    # For reference, at R ~ 50k, lambda = 5000A, 1 pixel is roughly 0.04A
    toler = 0.001

    import numpy as np
    from arclines import build_lists
    from arclines import arcio as arcl_io

    print("=============================================================")
    print("This script is for EXPERTS ONLY")
    print("Continue only if you know what you are doing")
    print("Otherwise exit")
    print("p.s.  You need to remove the files you wish to re-build")
    print("=============================================================")
    if not pargs.skip_stop:
        pass
        #pdb.set_trace()

    # Load the NIST ThAr list
    llist = arcl_io.load_line_lists(["ThAr"], NIST=True)
    llist['Ion'] = llist['Spectrum']
    # ['Spectrum', 'wave', 'Aki', 'Rel.', 'Ion', 'NIST']

    # Generate a table
    linelist = build_lists.init_line_list()

    # Check that all lines in the following list are included. If not, add them.
    # Murphy et al. (2007), MNRAS, 378, 221
    # http://www.astronomy.swin.edu.au/~mmurphy/thar/thar.html
    vwn = np.loadtxt(src_path+"ThAr_Murphy2007.dat", unpack=True, usecols=(0,))
    vwv = 1.0E8 / vwn  # Convert vacuum wavenumber (in cm^-1) to vacuum wavelength (in Angstroms)
    elm, ion = np.loadtxt(src_path+"ThAr_Murphy2007.dat", unpack=True, dtype='|S4', usecols=(3, 4))
    elm = np.core.defchararray.add(elm, np.array([' ']*elm.size))
    elmion = np.core.defchararray.add(elm, ion)
    wxx = np.where(elmion == 'XX 0')
    elmion[wxx] = 'UNKNWN'
    flagnist = np.zeros(vwv.size)

    # now add all NIST lines
    nlines = llist['Ion'].size
    for ll in range(nlines):
        # Check if this line is in the Murphy list
        wm = np.where((elmion == llist['Ion'][ll]) & (np.abs(vwv-llist['wave'][ll]) < toler))
        if wm[0].size == 1:
            if flagnist[wm[0][0]] != 0:
                print("WARNING :: line already flagged in Murphy catalog.")
                print("Possibly a duplicate NIST entry?", llist['Ion'][ll], llist['wave'][ll])
                print("Ignoring this line\n---")
                continue
            linelist.add_row([llist['Ion'][ll], llist['wave'][ll], 1, 0, llist['Rel.'][ll], 'NIST+MURPHY'])
            flagnist[wm[0][0]] = 1
        elif wm[0].size == 0:
            linelist.add_row([llist['Ion'][ll], llist['wave'][ll], 1, 0, llist['Rel.'][ll], 'NIST'])
        else:
            print("ERROR :: too many lines match. Reduce tolerance")
            sys.exit()
        if ll % 1000 == 0:
            print(ll, '/', nlines)
    # Remove the first dummy row
    linelist.remove_row(0)

    # Now include all Murphy lines, which are not in NIST
    wm = np.where(flagnist == 0)[0]
    if wm.size != 0:
        print("{0:d} lines not in Murphy list, but not in NIST".format(wm.size))
        for ll in range(wm.size):
            linelist.add_row([elmion[wm[ll]], vwv[wm[ll]], 0, 0, 1.0, 'MURPHY'])  # Relative Intensity is just a dummy number here

    # Finally, sort the list by increasing wavelength
    linelist.sort('wave')

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
