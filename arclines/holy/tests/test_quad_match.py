""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import h5py
import pdb

from arclines import io as arcl_io
from arclines.holy import patterns as arch_patt

from xastropy.xutils import xdebug as xdb


def grade_fidx_results(fidx):
    """  Grades the quad_match results
    Parameters
    ----------
    fidx

    Returns
    -------

    """
    grades = {}
    grades['ndetect'] = len(fidx)
    for key in ['nOK', 'nRisk', 'nPerf', 'nFail', 'nAmb', 'nGood']:
        grades[key] = 0
    # Loop on indices
    for key in fidx.keys():
        matches = np.array(fidx[key]['matches'])
        nmatch = matches.size
        uni, counts = np.unique(matches, return_counts=True)
        nanswers = len(uni)
        # Truth
        truth = uni == fidx[key]['truth']
        itruth = np.where(truth)[0]
        if len(itruth) == 1:
            ntruth = counts[itruth[0]]
        else:
            ntruth = 0
        # Bad
        if np.sum(~truth) > 0:
            max_bad = max(counts[~truth])
        else:
            max_bad = 0
        # Score
        if (ntruth == nmatch) & (nmatch >= 4):
            grades['nPerf'] += 1
        elif nmatch == 0:
            grades['nAmb'] += 1
        elif (ntruth == 4) & (nmatch == 5):
            grades['nGood'] += 1
        elif (ntruth/nmatch >= 2./3) & (nmatch >= 6):
            grades['nGood'] += 1
        elif (ntruth == 3) & (nmatch == 3):
            grades['nGood'] += 1
        elif (ntruth == 3) & (nmatch == 4):
            grades['nOK'] += 1
        elif (ntruth == 2) & (nmatch == 2):
            grades['nRisk'] += 1
        elif (nmatch >= 3) & (max_bad >= 3) & (max_bad > ntruth):
            grades['nFail'] += 1
            #print(fidx[key])
            #pdb.set_trace()
        elif (nmatch in [2,3]) & (max_bad >= 2):
            grades['nFail'] += 1
        else:
            grades['nAmb'] += 1
    # Return
    return grades


def test_quad_match_with_lowredux(low_redux_hdf, instr, swv_uncertainty=500.):
    """
    Returns
    -------

    """
    # Load for instrument
    if instr == 'LRISb':
        llist = arcl_io.load_line_lists(['CdI','HgI','ZnI'], unknown=True)
        # Cut
        gdwv = (llist['wave'] > 3000.) & (llist['wave'] < 5600.)
        cut_llist = llist[gdwv]
        #
        wvdata = cut_llist['wave'].data
        # Add a key missing line
        #wvdata = np.append(wvdata, [3404.6978])  
        # Sort
        wvdata.sort()
        #
        disp = 1.26  # Ang/binned pix
        pix_tol = 2
        cut_choice = 2  # 0 might be a touch better..
    elif instr == 'LRISr':
        #llist = arcl_io.load_line_lists(['ArI','HgI','KrI','NeI','XeI'], unknown=True)
        llist = arcl_io.load_line_lists(['ArI','HgI','NeI'], unknown=True)  # No Kr or Xe in LowRedux
        # Cut
        gdwv = (llist['wave'] > 5600.) & (llist['wave'] < 9000.)  # WILL WANT TO EXTEND
        cut_llist = llist[gdwv]
        #
        wvdata = cut_llist['wave'].data
        # Add a key missing line
        #wvdata = np.append(wvdata, [3404.6978])
        # Sort
        wvdata.sort()
        #
        disp = 1.6  # Ang/binned pix
        pix_tol = 1
        cut_choice = 2  # ARBITRARY
    else:
        raise IOError("Not ready for this instrument")

    # Open
    hdf = h5py.File(low_redux_hdf,'r')
    mdict = {}
    for key in hdf['meta'].keys():
        mdict[key] = hdf['meta'][key].value

    # Loop on spec
    extras = []
    for ispec in range(mdict['nspec']):
        all_tcent = hdf['arcs/'+str(ispec)+'/pixpk'].value

        spec = hdf['arcs/'+str(ispec)+'/spec'].value
        wave = hdf['arcs/'+str(ispec)+'/wave'].value # vacuum
        npix = wave.size
        #
        if False:
            from scipy.interpolate import interp1d
            fwv = interp1d(np.arange(npix), wave, kind='cubic')
            #fwv(all_tcent[3])
            xdb.set_trace()
        amps = []
        for itc in all_tcent:
            pix = int(np.round(itc))
            amps.append(spec[pix])
        amps = np.array(amps)

        # Trim tcent on amplitude
        if cut_choice == 0:
            cut_amp = amps > 1000.
        elif cut_choice == 1:
            mxa = np.max(amps)
            cut_amp = amps > 0.2*mxa
        elif cut_choice == 2:
            cut_amp = amps > 500.
        tcent = all_tcent[cut_amp]
        nlin = tcent.size

        # init with Truth
        final_idx = {}
        for ii in range(nlin):
            final_idx[ii] = {}
            final_idx[ii]['matches'] = []
            # Truth (if any)
            widx = int(np.round(tcent[ii]))
            mtw = np.where(np.abs(wvdata-wave[widx]) < 2*disp)[0]  # Catches bad LRISb line
            if len(mtw) == 0:
                if (instr=='LRISb') & (wave[widx] < 5600): # LRISb only
                    extras.append(wave[widx])
                #    print("No match for index={:d}, wave={:g}, amp={:g}".format( ii,wave[widx],spec[widx]))
                final_idx[ii]['truth'] = -1
            elif len(mtw) == 1:
                final_idx[ii]['truth'] = mtw[0]
            else:
                if instr == 'LRISb':
                    final_idx[ii]['truth'] = mtw[0] # Might have had this wrong! Not 4681.45
                else:
                    pdb.set_trace()
        #
        for idx in range(nlin-4):
            for jj in range(4):
                sub_idx = idx + np.arange(5).astype(int)
                msk = np.array([True]*5)
                msk[jj+1] = False
                # Setup
                sidx = sub_idx[msk]
                spec_lines = np.array(tcent)[sidx]
                #
                widx = int(np.round(tcent[idx]))
                wvmnx = [wave[widx]-swv_uncertainty, wave[widx]+swv_uncertainty]
                if idx == 0:
                    twv_min = wave[widx]
                # Run
                matches = arch_patt.match_quad_to_list(spec_lines, wvdata, wvmnx, disp, tol=pix_tol)
                # Save
                for match in matches:
                    for ii in range(4):
                        final_idx[sidx[ii]]['matches'].append(match[ii])
        '''
        for idx in range(nlin-4):
            # Setup
            spec_lines = np.array(tcent[idx:idx+4])
            #
            widx = int(np.round(tcent[idx]))
            wvmnx = [wave[widx]-swv_uncertainty, wave[widx]+swv_uncertainty]
            if idx == 0:
                twv_min = wave[widx]
            # Run
            matches = arch_patt.match_quad_to_list(spec_lines, wvdata, wvmnx, disp, tol=pix_tol)
            # Save
            for match in matches:
                for ii in range(4):
                    final_idx[idx+ii]['matches'].append(match[ii])
        '''
        # Grade
        grades = grade_fidx_results(final_idx)
        # PRINT
        if ispec == 0:
            print("II nDet nPerf nGood nOK nRisk nAmb nFail wvmin")
        print("{:2d} {:3d}  {:3d}  {:3d}   {:3d}   {:3d}  {:3d}  {:3d}   {:g}".format(
            ispec, grades['ndetect'], grades['nPerf'], grades['nGood'], grades['nOK'],
            grades['nRisk'], grades['nAmb'], grades['nFail'],
                twv_min))
        if ispec == 94:
            xdb.set_trace()
    extras = np.array(extras)
    extras.sort()
    xdb.set_trace()


def main(flg_tst):
    import arclines
    test_arc_path = arclines.__path__[0]+'/data/test_arcs/'

    # Test on LRISb_600
    if (flg_tst % 2**1) >= 2**0:
        hdf_file = test_arc_path+'LRISb_600_LRX.hdf5'  # Create with low_redux.py if needed
        test_quad_match_with_lowredux(hdf_file, 'LRISb')
    # Test on LRISr_600
    if (flg_tst % 2**2) >= 2**1:
        hdf_file = test_arc_path+'LRISr_600_7500.hdf5'  # Create with low_redux.py if needed
        test_quad_match_with_lowredux(hdf_file, 'LRISr')


# Test
if __name__ == '__main__':
    flg_tst = 0
    flg_tst += 2**0   # LRISb 600
    #flg_tst += 2**1   # LRISr 600

    main(flg_tst)