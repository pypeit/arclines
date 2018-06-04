""" Hiding calls to PYPIT used in arclines
Avoid double dependency if possible
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import pdb

from arclines import utils as al_utils


def find_peaks(censpec):
    """
    Parameters
    ----------
    censpec
    siglev
    bpfit : int, optional
      Order for background continuum

    Returns
    -------
    tampl, tcent, twid, w, detns

    """
    fitp = 7  #slf._argflag['arc']['calibrate']['nfitpix']
    if len(censpec.shape) == 3:
        detns = censpec[:, 0].flatten()
    else:
        detns = censpec.copy()
    xrng = np.arange(float(detns.size))

    # Find all significant detections
    pixt = np.where((detns > 0.0) &
                    (detns > np.roll(detns, 1)) & (detns >= np.roll(detns, -1)) &
                    (np.roll(detns, 1) > np.roll(detns, 2)) & (np.roll(detns, -1) > np.roll(detns, -2)) &
                    (np.roll(detns, 2) > np.roll(detns, 3)) & (np.roll(detns, -2) > np.roll(detns, -3)))[0]
    tampl, tcent, twid = fit_arcspec(xrng, detns, pixt, fitp)
    w = np.where((np.isnan(twid) == False) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent > 0.0) & (tcent < xrng[-1]))
    # Return
    return tampl, tcent, twid, w, detns


def fit_arcspec(xarray, yarray, pixt, fitp):

    # Setup the arrays with fit parameters
    sz_p = pixt.size
    sz_a = yarray.size
    ampl, cent, widt = -1.0*np.ones(sz_p, dtype=np.float),\
                       -1.0*np.ones(sz_p, dtype=np.float),\
                       -1.0*np.ones(sz_p, dtype=np.float)

    for p in range(sz_p):
        pmin = pixt[p]-(fitp-1)//2
        pmax = pixt[p]-(fitp-1)//2 + fitp
        if pmin < 0:
            pmin = 0
        if pmax > sz_a:
            pmax = sz_a
        if pmin == pmax:
            continue
        if pixt[p]-pmin <= 1 or pmax-pixt[p] <= 1:
            continue  # Probably won't be a good solution
        # Fit the gaussian
        try:
            popt = al_utils.func_fit(xarray[pmin:pmax], yarray[pmin:pmax], "gaussian", 3)
            ampl[p] = popt[0]
            cent[p] = popt[1]
            widt[p] = popt[2]
        except RuntimeError:
            pass
    return ampl, cent, widt
