""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import pdb

from arclines import io as arcl_io
from arclines.holy import patterns as arch_patt
from arclines.holy import fitting as arch_fit
from arclines.pypit_utils import find_peaks


def basic(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
          swv_uncertainty=350., pix_tol=2, plot_fil=None, min_match=5):
    """ Basic holy grail algorithm

    Parameters
    ----------
    spec
    lines
    wv_cen
    disp
    siglev
    min_ampl
    swv_uncertainty
    pix_tol
    plot_fil

    Returns
    -------
    status : int

    """

    # Init line-lists and wavelength 'guess'
    npix = spec.size
    wave = wv_cen + (np.arange(npix) - npix/2.)*disp

    line_lists = arcl_io.load_line_lists(lines, unknown=True)
    wvdata = line_lists['wave'].data  # NIST + Extra
    isrt = np.argsort(wvdata)
    wvdata = wvdata[isrt]

    # Find peaks
    tampl, tcent, twid, w, yprep = find_peaks(spec, siglev=siglev)
    all_tcent = tcent[w]
    all_tampl = tampl[w]

    # Cut on Amplitude??
    cut_amp = all_tampl > min_ampl
    cut_tcent = all_tcent[cut_amp]
    icut = np.where(cut_amp)[0]
    nlin = cut_tcent.size

    # Matching
    match_idx, scores = arch_patt.run_quad_match(cut_tcent, wave, wvdata,
                                                 disp, swv_uncertainty=swv_uncertainty,
                                                 pix_tol=pix_tol)

    # Check quadrants
    xquad = npix//4 + 1
    print("================================================================")
    print("Checking quadrants:")
    print("----------------------------------------------------------------")
    for jj in range(4):
        tc_in_q = (cut_tcent >= jj*xquad) & (cut_tcent < (jj+1)*xquad)
        cstat = 'quad {:d}: ndet={:d}'.format(jj, np.sum(tc_in_q))
        # Stats
        for key in ['Perf', 'Good', 'OK', 'Amb']:
            in_stat = scores[tc_in_q] == key
            cstat += ' {:s}={:d}'.format(key, np.sum(in_stat))
        # Print
        print(cstat)
    print("----------------------------------------------------------------")

    # Go for it!?
    mask = np.array([False]*len(all_tcent))
    IDs = []
    for kk,score in enumerate(scores):
        if score in ['Perf', 'Good', 'Ok']:
            mask[icut[kk]] = True
            uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
            imx = np.argmax(counts)
            IDs.append(wvdata[uni[imx]])
    ngd_match = np.sum(mask)
    if ngd_match < min_match:
        print("Insufficient matches to continue")
        status = -1
        return status, ngd_match, match_idx, scores, None

    # Fit
    NIST_lines = line_lists['NIST'] > 0
    ifit = np.where(mask)[0]
    final_fit = arch_fit.iterative_fitting(spec, all_tcent, ifit,
                               IDs, line_lists[NIST_lines], disp, plot_fil=plot_fil)
    # Return
    status = 1
    return status, ngd_match, match_idx, scores, final_fit