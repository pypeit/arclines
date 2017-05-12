#!/usr/bin/env python
"""
Match input spectrum to ID lines
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
        description='Match input spectrum to arclines line lists')
    parser.add_argument("spectrum", type=str, help="Spectrum file (.ascii, .fits, .json)")
    parser.add_argument("wvcen", type=float, help="Guess at central wavelength (within 1000A)")
    parser.add_argument("disp", type=float, help="Accurate dispersion (Ang/pix)")
    parser.add_argument("lines", type=str, help="Comma separated list of lamps")
    parser.add_argument("--outroot", type=str, help="Root filename for plot, IDs")
    parser.add_argument("--min_ampl", type=float, help="Minimum amplitude for line analysis [default: 100.]")
    parser.add_argument("--debug", default=False, action='store_true', help="Debug")
    parser.add_argument("--fit", default=False, action='store_true', help="Fit the lines?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(pargs=None):
    """ Run
    Parameters
    ----------
    args

    Returns
    -------

    """
    import numpy as np

    from linetools import utils as ltu

    from arclines import io as arcl_io
    from arclines.holy import utils as arch_utils
    from arclines.holy.grail import semi_brute
    from arclines.holy import patterns as arch_patt
    from arclines.holy import fitting as arch_fit


    if pargs.outroot is None:
        pargs.outroot = 'tmp_matches'
    # Defaults
    min_ampl = (pargs.min_ampl if (pargs.min_ampl is not None) else 100.)

    # Load spectrum
    spec = arcl_io.load_spectrum(pargs.spectrum)
    # Arc lines
    lines = pargs.lines.split(',')

    # Call brute
    best_dict, final_fit = semi_brute(spec, lines, pargs.wvcen, pargs.disp, min_ampl=min_ampl,
               debug=pargs.debug, outroot=pargs.outroot, do_fit=pargs.fit,
               verbose=True)
    if pargs.debug:
        pdb.set_trace()

    if pargs.fit:
        ltu.savejson(pargs.outroot+'_fit.json', ltu.jsonify(final_fit), easy_to_read=True, overwrite=True)


