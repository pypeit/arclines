""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)


import numpy as np
import pdb


def match_quad_to_list(spec_lines, line_list, wv_guess, dwv_guess,
                  tol=2., dwv_uncertainty=0.2, min_ftol=0.005):
    """
    Parameters
    ----------
    spec_lines : ndarray
      pixel space
    line_list
    tol
    min_ftol : float, optional
      Minimum tolerance for matching 

    Returns
    -------
    possible_matches : list
      list of indices of matching quads

    """
    # Setup spec
    npix = spec_lines[-1]-spec_lines[0]
    spec_values = (spec_lines[1:-1]-spec_lines[0])/(
        spec_lines[-1]-spec_lines[0])
    ftol = max(tol/npix, min_ftol)
    #
    possible_start = np.where((line_list > wv_guess[0]) & (line_list < wv_guess[1]))[0]
    possible_matches = []
    for start in possible_start:
        #print("Trying {:g}".format(line_list[start]))
        # Find possible ends
        possible_ends = np.where( (line_list > line_list[start] + npix*dwv_guess*(1-dwv_uncertainty)) &
                                     (line_list < line_list[start] + npix*dwv_guess*(1+dwv_uncertainty)))[0]
        # Loop on ends
        for end in possible_ends:
            values = (line_list[start+1:end]-line_list[start]) / (
                line_list[end]-line_list[start])
            # Test
            diff0 = np.abs(values-spec_values[0])
            tst0 = diff0 < ftol
            diff1 = np.abs(values-spec_values[1])
            tst1 = diff1 < ftol
            #if np.abs(line_list[start]-6097.8) < 0.2:
            #    debugger.set_trace()
            if np.any(tst0) & np.any(tst1):
                i0 = np.argmin(diff0)
                i1 = np.argmin(diff1)
                #if np.sum(tst0) > 1:
                #    pdb.set_trace()
                #possible_matches.append([start, start+1+np.where(tst0)[0][0],
                #                         start+1+np.where(tst1)[0][0], end])
                possible_matches.append([start, start+1+i0, start+1+i1, end])
    # Return
    return possible_matches


def run_quad_match(tcent, twave, llist_wv, disp, swv_uncertainty=250.,
                   pix_tol=1.):
    """
    Parameters
    ----------
    tcent : ndarray
      Pixel positions of arc lines
    twave : ndarray
      Crude guess at wavelength solution, e.g. from wvcen, disp
    llist_wv : ndarray
      Lines to match against (from a line list)
    pix_tol : float
      Tolerance in units of pixels to match to

    Returns
    -------
    match_idx : dict
      Record of matches
    scores : ndarray
      str array of scores
    """

    # Init
    nlin = tcent.size
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
            spec_lines = np.array(tcent)[sidx]
            #
            widx = int(np.round(tcent[idx]))
            wvmnx = [twave[widx]-swv_uncertainty, twave[widx]+swv_uncertainty]
            # Run
            matches = match_quad_to_list(spec_lines, llist_wv, wvmnx, disp, tol=pix_tol)
            # Save
            for match in matches:
                for ii in range(4):
                    match_idx[sidx[ii]]['matches'].append(match[ii])
            #pdb.set_trace()
    # Score
    scores = score_quad_matches(match_idx)
    scores = np.array(scores)

    # Return
    return match_idx, scores


def scan_for_matches(wvcen, disp, npix, cut_tcent, wvdata, best_dict=None,
                     swv_uncertainty=350., wvoff=1000., pix_tol=2., ampl=None):
    """
    Parameters
    ----------
    wvcen : float
      Guess at central wavelength
    disp : float
    npix
    cut_tcent
    wvdata
    best_dict
    swv_uncertainty
    wvoff
    pix_tol

    Returns
    -------
    best_dict is updated in place
    """

    # Setup
    #wvoff=10.
    #pdb.set_trace()
    dcen = swv_uncertainty*0.8
    wvcens = np.arange(wvcen-wvoff, wvcen+wvoff+dcen, dcen)
    # Best
    if best_dict is None:
        best_dict = dict(nmatch=0, ibest=-1, bwv=0.)

    # Scan on wv_cen
    for ss,iwv_cen in enumerate(wvcens):
        # Wavelength array
        wave = iwv_cen + (np.arange(npix) - npix/2.)*disp
        match_idx, scores = run_quad_match(cut_tcent, wave, wvdata, disp,
                                           swv_uncertainty=swv_uncertainty,
                                           pix_tol=pix_tol)
        # Score
        mask = np.array([False]*len(cut_tcent))
        IDs = []
        for kk,score in enumerate(scores):
            if score in ['Perf', 'Good', 'Ok']:
                mask[kk] = True
                uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
                imx = np.argmax(counts)
                IDs.append(wvdata[uni[imx]])
            else:
                IDs.append(0.)
        ngd_match = np.sum(mask)
        # Update in place
        if ngd_match > best_dict['nmatch']:
            best_dict['nmatch'] = ngd_match
            best_dict['midx'] = match_idx
            best_dict['mask'] = mask
            best_dict['scores'] = scores
            best_dict['ibest'] = ss
            best_dict['bwv'] = iwv_cen
            best_dict['IDs'] = IDs
            # Search parameters
            best_dict['swv_uncertainty'] = swv_uncertainty
            best_dict['wvoff'] = wvoff
            best_dict['pix_tol'] = pix_tol
            best_dict['ampl'] = ampl

def score_quad_matches(fidx):
    """  Grades quad_match results
    Parameters
    ----------
    fidx

    Returns
    -------
    scores : list

    """
    # Loop on indices
    scores = []
    for key in fidx.keys():
        if len(fidx[key]['matches']) == 0:
            scores.append('None')
            continue
        matches = np.array(fidx[key]['matches'])
        nmatch = matches.size
        uni, counts = np.unique(matches, return_counts=True)
        nuni = len(uni)
        max_counts = max(counts)
        # Score
        if (nuni==1) & (nmatch >= 4):
            scores.append('Perf')
        elif nmatch == 0:
            scores.append('None')
        elif (max_counts == 4) & (nmatch == 5):
            scores.append('Good')
        elif (max_counts/nmatch >= 2./3) & (nmatch >= 6):
            scores.append('Good')
        elif (nuni == 1) & (nmatch == 3):
            scores.append('Good')
        elif (max_counts == 3) & (nmatch == 4):
            scores.append('OK')
        elif (nuni == 1) & (nmatch == 2):
            scores.append('Risk')
        else:
            scores.append('Amb')
    # Return
    return scores
