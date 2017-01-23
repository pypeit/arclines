#!/usr/bin/env python
"""
Add a new source to the line lists
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Add a new source to the database')
    parser.add_argument("-w", "--write", default=False, action='store_true', help="Actually write files?")
    parser.add_argument("-s", "--skip_stop", default=False, action='store_true', help="Skip warning stop?")
    parser.add_argument("--no_unknowns", default=False, action='store_true', help="Skip UNKNOWNS in source")
    parser.add_argument("--plots", default=False, action='store_true', help="Create plots?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    """ Run
    Parameters
    ----------
    args

    Returns
    -------

    """
    # Grab arguments
    pargs = parser(options=args)

    import numpy as np
    from astropy.table import vstack
    from arclines import build_lists
    from arclines import io as arcl_io
    from arclines import load_source
    from arclines import plots as arcl_plots
    from arclines import utils as arcl_utils


    print("=============================================================")
    print("This script is for EXPERTS ONLY")
    print("Continue only if you know what you are doing")
    print("This script adds only the *last* source listed")
    print("Otherwise exit")
    print("=============================================================")
    if not pargs.skip_stop:
        pdb.set_trace()

    # Load sources
    sources = arcl_io.load_source_table()
    source = sources[-1]

    # Load line lists
    uions = arcl_utils.unique_ions(source)
    llist_dict = {}
    for ion in uions:
        try:
            llist_dict[ion] = arcl_io.load_line_list(ion, use_ion=True)
        except IOError:
            print("No linelist found for ion={:s}.  Will create one".format(ion))


    # IDs
    print("=============================================================")
    print("Working on adding IDs from source {:s}".format(source['File']))
    print("=============================================================")
    llist_dict = build_lists.source_to_line_lists(source, write=pargs.write,
                                     llist_dict=llist_dict)
    # Create line list
    ll = []
    for key in llist_dict.keys():
        ll.append(llist_dict[key])
    line_lists = vstack(ll)

    # Purge unknowns
    print("=============================================================")
    print("Purging UNKNOWNs")
    print("=============================================================")
    build_lists.purge_unknowns(line_lists, write=pargs.write)

    print("=============================================================")
    print("Working on adding Unknowns from source {:s}".format(source['File']))
    print("=============================================================")
    if not pargs.no_unknowns:
        unknwns = build_lists.source_to_unknowns(source, write=pargs.write)
        unknwns.remove_column('line_flag')
    else:
        unknwns = None

    if unknwns is not None:
        ll.append(unknwns)
    # Loop to my loop
    print("=============================================================")
    print("Working on plot for {:s}".format(source['File']))
    print("=============================================================")
    # Load
    src_dict = load_source.load(source)
    uions = arcl_utils.unique_ions(source, src_dict=src_dict)
    src_dict['uions'] = uions
    # Outfile
    iext = source['File'].rfind('.')
    outfile = source['File'][:iext]+'.pdf'
    title = source['File'][:iext]
    #
    arcl_plots.show_source(src_dict, line_lists, outfile, title=title, clobber=True)


if __name__ == '__main__':
    main()
