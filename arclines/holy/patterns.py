""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np


def match_quad_to_list(spec_lines, line_list, wv_guess, dwv_guess,
                  tol=2., dwv_uncertainty=0.2):
    """
    Parameters
    ----------
    spec_lines : ndarray
      pixel space
    line_list
    tol

    Returns
    -------
    possible_matchtes : list
      list of indices of matching quads

    """
    # Setup spec
    npix = spec_lines[-1]-spec_lines[0]
    spec_values = (spec_lines[1:-1]-spec_lines[0])/(
        spec_lines[-1]-spec_lines[0])
    ftol = tol/npix
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
            tst0 = np.abs(values-spec_values[0]) < ftol
            tst1 = np.abs(values-spec_values[1]) < ftol
            #if np.abs(line_list[start]-6097.8) < 0.2:
            #    debugger.set_trace()
            if np.any(tst0) & np.any(tst1):
                possible_matches.append([start, start+1+np.where(tst0)[0][0],
                                       start+1+np.where(tst1)[0][0], end])
    # Return
    return possible_matches
