""" Module for running tests on spectrum to solution
Yes, the whole enchilada
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import pdb

from arclines.holy import grail

#from xastropy.xutils import xdebug as xdb


def test_basic(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
               swv_uncertainty=350., pix_tol=2, plot_fil=None):
    """ As named
    """
    stuff = grail.basic(spec, lines, wv_cen, disp, siglev=siglev, min_ampl=min_ampl,
                   swv_uncertainty=swv_uncertainty, pix_tol=pix_tol, plot_fil=plot_fil)


def test_unknwn_wvcen(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
               swv_uncertainty=350., pix_tol=2, plot_fil=None,
                      wvoff=1000.):
    """ As named
    """
    dcen = swv_uncertainty*0.8
    wvcens = np.arange(wv_cen-wvoff, wv_cen+wvoff+dcen, dcen)
    # Best
    best_dict = dict(nmatch=0, nID=0, rms=0., ibest=-1, bwv=0.)
    # Loop
    for ss,iwv_cen in enumerate(wvcens):
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
            nID = len(final_fit['xfit'])
            print("   Nmatch={:d}, nID={:d}, RMS={:g}".format(nmatch, nID,
                                                              final_fit['rms']))
            # Best?
            if nID > best_dict['nID']:
                best_dict['nmatch'] = nmatch
                best_dict['nID'] = nID
                best_dict['rms'] = final_fit['rms']
                best_dict['ibest'] = ss
                best_dict['bwv'] = iwv_cen
                best_dict['fit'] = final_fit
    # Report
    print("Here is the best:")
    print(best_dict)
    # QA
    if plot_fil is not None:
        from pypit import pyputils
        msgs = pyputils.get_dummy_logger()
        from pypit import arqa
        arqa.arc_fit_qa(None, best_dict['fit'], outfil=plot_fil)


def main(flg_tst):
    import json
    import h5py
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

    # Varying wvcen on LRISb
    if (flg_tst % 2**3) >=2**2:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/sources/'
        src_file = 'lrisb_600_4000_PYPIT.json'
        with open(test_arc_path+src_file,'r') as f:
            pypit_fit = json.load(f)
        # Run
        test_unknwn_wvcen(np.array(pypit_fit['spec']), ['CdI','HgI','ZnI'],
                   4400., 1.26, plot_fil='lrisb_scan_fit.pdf')

    # Varying wvcen on LRISb, LowRedux and far off center
    if (flg_tst % 2**4) >= 2**3:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/test_arcs/'
        hdf_file = test_arc_path+'LRISb_600_LRX.hdf5'  # Create with low_redux.py if needed
        hdf = h5py.File(hdf_file,'r')
        spec = hdf['arcs/18/spec'].value
        # Run
        test_unknwn_wvcen(spec, ['CdI','HgI','ZnI'],
                   5000., 1.26, plot_fil='lrisb_off_fit.pdf')

    # Varying wvcen on LRISr
    if (flg_tst % 2**5) >=2**4:
        # Load spectrum
        test_arc_path = arclines.__path__[0]+'/data/sources/'
        src_file = 'lrisr_600_7500_PYPIT.json'
        with open(test_arc_path+src_file,'r') as f:
            pypit_fit = json.load(f)
        # Run
        test_unknwn_wvcen(np.array(pypit_fit['spec']), ['ArI','HgI','KrI','NeI','XeI'],
                   7000., 1.6, plot_fil='lrisr_fit.pdf')

# Test
if __name__ == '__main__':
    flg_tst = 0
    #flg_tst += 2**0   # LRISb 600
    #flg_tst += 2**1   # Kastb
    #flg_tst += 2**2   # LRISb with unknown wv_cen
    #flg_tst += 2**3   # LRISb off to red with unknown wv_cen
    flg_tst += 2**4   # LRISr nominal

    main(flg_tst)