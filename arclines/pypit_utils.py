""" Hiding calls to PYPIT used in arclines
Avoid double dependency if possible
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import pdb


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
    from pypit import ararc
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
    tampl, tcent, twid = ararc.fit_arcspec(xrng, detns, pixt, fitp)
    w = np.where((np.isnan(twid) == False) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent > 0.0) & (tcent < xrng[-1]))
    # Return
    return tampl, tcent, twid, w, detns

