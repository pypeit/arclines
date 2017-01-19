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
    parser.add_argument("spectrum", type=str, help="Spectrum file (.ascii, .fits)")
    parser.add_argument("wvcen", type=float, help="Guess at central wavelength (within 1000A)")
    parser.add_argument("disp", type=float, help="Accurate dispersion (Ang/pix)")
    parser.add_argument("lines", type=str, help="Comma separated list of lamps")
    parser.add_argument("--outroot", default='tmp_matches', action='store_true', help="Root filename for plot, IDs")
    parser.add_argument("--min_ampl", type=float, help="Minimum amplitude for line analysis [default: 100.]")
    parser.add_argument("--debug", default=False, action='store_true', help="Debug")

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
    from astropy.table import vstack

    from linetools import utils as ltu

    from arclines import io as arcl_io
    from arclines import plots as arcl_plots
    from arclines.holy import utils as arch_utils
    from arclines.holy import patterns as arch_patt

    # Defaults
    min_ampl = (pargs.min_ampl if (pargs.min_ampl is not None) else 100.)

    # Load spectrum
    spec = arcl_io.load_spectrum(pargs.spectrum)
    npix = spec.size

    # Lines
    all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=min_ampl)

    # Load line lists
    lines = pargs.lines.split(',')
    line_lists = arcl_io.load_line_lists(lines)
    unknwns = arcl_io.load_unknown_list(lines)

    #delta_wv = npix * pargs.disp
    #in_spec = (wvdata > pargs.wvcen-delta_wv/2) & (wvdata < pargs.wvcen+delta_wv/2)
    #nline_in = np.sum(in_spec)

    #print("{:d} lines detected in the spectrum.".format(cut_tcent.size))
    #print("Approximately {:d} lines to match against in your nominal wavelength range".format(nline_in))

    # 3 things to fiddle:
    #  pix_tol -- higher for fewer lines  1/2
    #  unknowns -- on for fewer lines  off/on
    #  scoring -- weaken for more lines ??


    # Best
    best_dict = dict(nmatch=0, ibest=-1, bwv=0.)

    # Loop on unknowns
    for unknown in [False, True]:
        if unknown:
            tot_list = vstack([line_lists,unknwns])
        else:
            tot_list = line_lists
        wvdata = np.array(tot_list['wave'].data) # Removes mask if any
        wvdata.sort()
        sav_nmatch = best_dict['nmatch']

        # Loop on pix_tol
        for pix_tol in [1.,2.]:
            # Scan on wavelengths
            arch_patt.scan_for_matches(pargs.wvcen, pargs.disp, npix, cut_tcent, wvdata,
                               best_dict=best_dict, pix_tol=pix_tol)
        # Save linelist?
        if best_dict['nmatch'] > sav_nmatch:
            best_dict['line_list'] = tot_list
            best_dict['unknown'] = unknown

    # Report
    print('---------------------------------------------------')
    print('Report:')
    print('::   Number of lines recovered = {:d}'.format(all_tcent.size))
    print('::   Number of lines analyzed = {:d}'.format(cut_tcent.size))
    print('::   Number of Perf/Good/Ok matches = {:d}'.format(best_dict['nmatch']))
    print('::   Best central wavelength = {:g}A'.format(best_dict['bwv']))
    print('::   Best solution used pix_tol = {}'.format(best_dict['pix_tol']))
    print('::   Best solution had unknown = {}'.format(best_dict['unknown']))
    print('---------------------------------------------------')

    if pargs.debug:
        match_idx = best_dict['midx']
        for kk in match_idx.keys():
            uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
            print('kk={}, {}, {}, {}'.format(kk, uni, counts, np.sum(counts)))

    # Write scores
    #out_dict = best_dict['scores']
    #jdict = ltu.jsonify(out_dict)
    #ltu.savejson(pargs.outroot+'.scores', jdict, easy_to_read=True, overwrite=True)
    #print("Wrote: {:s}".format(pargs.outroot+'.scores'))

    # Write IDs
    out_dict = dict(pix=cut_tcent, IDs=best_dict['IDs'])
    jdict = ltu.jsonify(out_dict)
    ltu.savejson(pargs.outroot+'.json', jdict, easy_to_read=True, overwrite=True)
    print("Wrote: {:s}".format(pargs.outroot+'.json'))

    # Plot
    arcl_plots.match_qa(spec, cut_tcent, best_dict['line_list'],
                        best_dict['IDs'], best_dict['scores'], pargs.outroot+'.pdf')
    print("Wrote: {:s}".format(pargs.outroot+'.pdf'))
    if pargs.debug:
        pdb.set_trace()
