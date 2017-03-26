""" Utilities for building arcline lists
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import pdb


def unique_ions(source, src_dict=None):
    """ Unique ions from src_dict and source
    of just the source

    Parameters
    ----------
    src_dict
    source

    Returns
    -------
    uions : ndarray
      str array of unique ions

    """
    # Lines
    src_lines = source['Lines'].split(',')
    if src_dict is None:
        return src_lines
    # Check, as applicable
    if src_dict['ID_lines'] is not None:
        uions = np.unique(src_dict['ID_lines']['ion'].data)
        for src_line in src_lines:
            if src_line not in uions.tolist():
                raise ValueError("Line {:s} not found in ID_lines".format(src_line))
        return uions
    else:
        return src_lines


def vette_unkwn_against_lists(U_lines, uions, tol_NIST=0.2, NIST_only=False,
                              tol_llist=2., verbose=False):
    """ Query unknown lines against NIST database

    Parameters
    ----------
    U_lines : Table
    uions : list or ndarray
      list of lines to check against
    tol_NIST : float, optional
      Tolerance for a match with NIST
    tol_llist : float, optional
      Tolerance for a match with arclines line lists

    Returns
    -------
    mask : int array
      2 = NIST (and add)
      1 = Add these
      0 = Do not add these
    wv_match : ndarray
      str array

    """
    from arclines import io as arcl_io

    mask = np.ones(len(U_lines)).astype(int)
    wv_match = np.array(['XXI   12233.2312']*len(U_lines))
    # Loop on NIST
    for ion in uions:
        # Load
        nist = arcl_io.load_nist(ion)
        # Try to match
        for ss,row in enumerate(U_lines):
            dwv = np.abs(nist['wave']-row['wave'])
            imin = np.argmin(np.abs(dwv))
            #if verbose:
            #    print("Closest match to ion={:s} for {:g} is".format(ion,row['wave']))
            #    print(nist[['Ion','wave','RelInt']][imin])
            # Match?
            if dwv[imin] < tol_NIST:
                wv_match[ss] = '{:s} {:.4f}'.format(ion,nist['wave'][imin])
                mask[ss] = 2
                if verbose:
                    print("UNKNWN Matched to NIST: ion={:s} {:g} with {:g}".format(
                        ion,nist['wave'][imin], row['wave']))
                #print(nist[['Ion','wave','RelInt','Aki']][imin])
    if NIST_only:
        return mask, wv_match

    # Our line lists
    line_list = arcl_io.load_line_lists(uions, skip=True)
    if line_list is None:
        return mask, wv_match
    for ss,row in enumerate(U_lines):
        dwv = np.abs(line_list['wave']-row['wave'])
        imin = np.argmin(np.abs(dwv))
        # Match?
        if dwv[imin] < tol_llist:
            mask[ss] = 0
            if verbose:
                print("UNKNWN Matched to arclines: ion={:s} {:g} with {:g}".format(
                        line_list['ion'][imin], line_list['wave'][imin], row['wave']))
                print("  ---- Will not add it")
    return mask, wv_match

