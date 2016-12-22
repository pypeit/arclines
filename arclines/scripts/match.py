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
    #parser.add_argument("--plots", default=False, action='store_true', help="Create plots?")

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
    from arclines import plots as arcl_plots
    from arclines.holy import utils as arch_utils
    from arclines.holy import patterns as arch_patt

    # Load spectrum
    spec = arcl_io.load_spectrum(pargs.spectrum)
    npix = spec.size

    # Load line lists
    lines = pargs.lines.split(',')
    line_lists = arcl_io.load_line_lists(lines, unknown=True)
    wvdata = line_lists['wave'].data  # NIST + Extra
    isrt = np.argsort(wvdata)
    wvdata = wvdata[isrt]

    # Lines
    all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec)#, siglev=siglev, min_ampl=min_ampl)

    # Scan on wv_cen
    swv_uncertainty=350.
    wvoff=1000.
    dcen = swv_uncertainty*0.8
    wvcens = np.arange(pargs.wvcen-wvoff, pargs.wvcen+wvoff+dcen, dcen)
    # Best
    best_dict = dict(nmatch=0, ibest=-1, bwv=0.)
    for ss,iwv_cen in enumerate(wvcens):
        # Wavelength array
        wave = iwv_cen + (np.arange(npix) - npix/2.)*pargs.disp
        match_idx, scores = arch_patt.run_quad_match(cut_tcent, wave, wvdata,
                                                 pargs.disp, swv_uncertainty=swv_uncertainty,
                                                 pix_tol=2.)
        # Score
        mask = np.array([False]*len(all_tcent))
        IDs = []
        for kk,score in enumerate(scores):
            if score in ['Perf', 'Good', 'Ok']:
                mask[icut[kk]] = True
                uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
                imx = np.argmax(counts)
                IDs.append(wvdata[uni[imx]])
            else:
                IDs.append(0.)
        ngd_match = np.sum(mask)
        if ngd_match > best_dict['nmatch']:
            best_dict['nmatch'] = ngd_match
            best_dict['midx'] = match_idx
            best_dict['scores'] = scores
            best_dict['ibest'] = ss
            best_dict['bwv'] = iwv_cen
            best_dict['IDs'] = IDs
    # Report
    print('---------------------------------------------------')
    print('Report:')
    print('::   Number of lines analyzed = {:d}'.format(cut_tcent.size))
    print('::   Number of Perf/Good/Ok matches = {:d}'.format(best_dict['nmatch']))
    print('::   Best central wavelength = {:g}A'.format(best_dict['bwv']))
    print('---------------------------------------------------')

    # Write IDs
    out_dict = dict(pix=cut_tcent, IDs=best_dict['IDs'])
    jdict = ltu.jsonify(out_dict)
    ltu.savejson(pargs.outroot+'.json', jdict, easy_to_read=True)
    print("Wrote: {:s}".format(pargs.outroot+'.json'))

    # Plot
    arcl_plots.match_qa(spec, cut_tcent, line_lists[isrt],
                        best_dict['IDs'], best_dict['scores'], pargs.outroot+'.pdf')
    print("Wrote: {:s}".format(pargs.outroot+'.pdf'))
