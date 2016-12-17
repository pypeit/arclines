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
        description='Build the arclines line lists')
    parser.add_argument("-v", "--version", help="DB version to generate")
    parser.add_argument("-t", "--test", default=False, action='store_true', help="Test?")
    parser.add_argument("--write", default=False, action='store_true', help="Actually write files?")

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
    from arclines import build_lists

    # Grab arguments
    pargs = parser(options=args)

    build_lists.master_build(write=pargs.write)
    if not pargs.write:
        print("Ran script without writing.  Use --write to write")

if __name__ == '__main__':
    main()
