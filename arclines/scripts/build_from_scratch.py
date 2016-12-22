#!/usr/bin/env python
"""
Run a build of the line lists
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
        description='Build the arclines line lists from scratch')
    parser.add_argument("-w", "--write", default=False, action='store_true', help="Actually write files?")
    parser.add_argument("--skip_stop", default=False, action='store_true', help="Actually write files?")
    parser.add_argument("--unknowns", default=False, action='store_true', help="Actually write files?")
    #parser.add_argument("-s", "--step_by_step", default=False, action='store_true', help="Step by step build?")

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
    from arclines import build_lists
    from arclines import io as arcl_io


    print("=============================================================")
    print("This script is for EXPERTS ONLY")
    print("Continue if you know what you are doing")
    print("Otherwise exit")
    print("p.s.  You need to remove the files you wish to re-build")
    print("=============================================================")
    if not pargs.skip_stop:
        pdb.set_trace()

    # Load sources
    sources = arcl_io.load_source_table()

    # IDs
    llist_dict = {}
    if not pargs.unknowns:
        for kk,source in enumerate(sources):
            print("=============================================================")
            print("Working on adding IDs from source {:s}".format(source['File']))
            print("=============================================================")
            llist_dict = build_lists.source_to_line_lists(source, write=pargs.write,
                                             llist_dict=llist_dict)
        # By-hand
        build_lists.by_hand(llist_dict, write=pargs.write)

        # Write?
        if not pargs.write:
            print("=============================================================")
            print("Rerun with --write if you are happy with what you see.")
            print("=============================================================")
        return

if __name__ == '__main__':
    main()
