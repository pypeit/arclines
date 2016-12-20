""" Module for running tests on spectrum to solution
Yes, the whole enchilada
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import pdb

from arclines import io as arcl_io
from arclines.holy import patterns as arch_patt
from arclines.holy import fitting as arch_fit
from arclines.pypit_utils import find_peaks

from xastropy.xutils import xdebug as xdb


def test_enchilada(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
                   swv_uncertainty=350., pix_tol=2, plot_fil=None):
    """
    Returns
    -------

    """
    # Init
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
    # Init
    match_idx = {}
    for ii in range(nlin):
        match_idx[ii] = {}
        match_idx[ii]['matches'] = []
    # Run -- quad
    for idx in range(nlin-4):
        for jj in range(4):
            sub_idx = idx + np.arange(5).astype(int)
            msk = np.array([True]*5)
            msk[jj+1] = False
            # Setup
            sidx = sub_idx[msk]
            spec_lines = np.array(cut_tcent)[sidx]
            #
            widx = int(np.round(cut_tcent[idx]))
            wvmnx = [wave[widx]-swv_uncertainty, wave[widx]+swv_uncertainty]
            if idx == 0:
                twv_min = wave[widx]
            # Run
            matches = arch_patt.match_quad_to_list(spec_lines, wvdata, wvmnx, disp, tol=pix_tol)
            # Save
            for match in matches:
                for ii in range(4):
                    match_idx[sidx[ii]]['matches'].append(match[ii])
    # Score
    scores = arch_patt.score_quad_matches(match_idx)
    scores = np.array(scores)

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
    NIST_lines = line_lists['NIST'] > 0
    #pdb.set_trace()
    ifit = np.where(mask)[0]
    arch_fit.iterative_fitting(spec, all_tcent, ifit,
                               IDs, line_lists[NIST_lines], disp, plot_fil=plot_fil)

def main(flg_tst):
    import json
    import arclines

    # Test on LRISb_600 from PYPIT
    if (flg_tst % 2**1) >= 2**0:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/sources/'
        src_file = 'lrisb_600_4000_PYPIT.json'
        with open(test_arc_path+src_file,'r') as f:
            pypit_fit = json.load(f)
        # Run
        test_enchilada(np.array(pypit_fit['spec']), ['CdI','HgI','ZnI'],
                       4400., 1.26, plot_fil='lrisb_fit.pdf')

    # Test on Kastb from PYPIT
    if (flg_tst % 2**2) >= 2**1:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/sources/'
        src_file = 'kastb_600_PYPIT.json'
        with open(test_arc_path+src_file,'r') as f:
            pypit_fit = json.load(f)
        # Run
        test_enchilada(np.array(pypit_fit['spec']), ['CdI','HeI','HgI'],
                       4400., 1.02, plot_fil='kastb_fit.pdf')

# Test
if __name__ == '__main__':
    flg_tst = 0
    #flg_tst += 2**0   # LRISb 600
    flg_tst += 2**1   # Kastb

    main(flg_tst)