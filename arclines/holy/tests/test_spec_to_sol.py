""" Module for running tests on spectrum to solution
Yes, the whole enchilada
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import pdb

from arclines.holy import grail

from xastropy.xutils import xdebug as xdb


def test_basic(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
               swv_uncertainty=350., pix_tol=2, plot_fil=None):
    """ As named
    """
    stuff = grail.basic(spec, lines, wv_cen, disp, siglev=siglev, min_ampl=min_ampl,
                   swv_uncertainty=swv_uncertainty, pix_tol=pix_tol, plot_fil=plot_fil)


def test_unknwn_wvcen(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
               swv_uncertainty=350., pix_tol=2, plot_fil=None):
    """ As named
    """
    dcen = swv_uncertainty*0.8
    wvcens = np.arange(wv_cen-500., wv_cen+500+dcen, dcen)
    # Loop
    for iwv_cen in wvcens:
        print('wv_cen guess = {:g}'.format(iwv_cen))
        stuff = grail.basic(spec, lines, iwv_cen, disp, siglev=siglev, min_ampl=min_ampl,
                swv_uncertainty=swv_uncertainty, pix_tol=pix_tol, plot_fil=None)
        if stuff[0] == -1:
            print("Solution failed for iwv_cen={:g}. Nmatch={:d}".format(iwv_cen, stuff[1]))
        elif stuff[0] == 1:
            # Unpack
            status, nmatch, match_idx, scores, final_fit = stuff
            print("Fit finished for iwv_cen={:g}".format(iwv_cen))
            # Metrics -- nmatch, nID, RMS
            print("   Nmatch={:d}, nID={:d}, RMS={:g}".format(nmatch, len(final_fit['xfit']),
                                                              final_fit['rms']))

def main(flg_tst):
    import json
    import arclines

    # Basic test on LRISb_600 from PYPIT
    if (flg_tst % 2**1) >= 2**0:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/sources/'
        src_file = 'lrisb_600_4000_PYPIT.json'
        with open(test_arc_path+src_file,'r') as f:
            pypit_fit = json.load(f)
        # Run
        test_basic(np.array(pypit_fit['spec']), ['CdI','HgI','ZnI'],
                       4400., 1.26, plot_fil='lrisb_fit.pdf')

    # Test on Kastb from PYPIT
    if (flg_tst % 2**2) >= 2**1:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/sources/'
        src_file = 'kastb_600_PYPIT.json'
        with open(test_arc_path+src_file,'r') as f:
            pypit_fit = json.load(f)
        # Run
        test_basic(np.array(pypit_fit['spec']), ['CdI','HeI','HgI'],
                       4400., 1.02, plot_fil='kastb_fit.pdf')

    # Basic test on LRISb_600 from PYPIT
    if (flg_tst % 2**3) >= 2**2:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/sources/'
        src_file = 'lrisb_600_4000_PYPIT.json'
        with open(test_arc_path+src_file,'r') as f:
            pypit_fit = json.load(f)
        # Run
        test_unknwn_wvcen(np.array(pypit_fit['spec']), ['CdI','HgI','ZnI'],
                   4400., 1.26, plot_fil='lrisb_fit.pdf')

# Test
if __name__ == '__main__':
    flg_tst = 0
    #flg_tst += 2**0   # LRISb 600
    #flg_tst += 2**1   # Kastb
    flg_tst += 2**2   # LRISb with unknown wv_cen

    main(flg_tst)