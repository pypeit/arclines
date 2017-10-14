""" Hiding calls to PYPIT used in arclines
Avoid double dependency if possible
"""
from __future__ import absolute_import, division, print_function

import numpy as np


def find_peaks(censpec, siglev=6., bpfit=5):
    """
    Parameters
    ----------
    censpec
    siglev
    bpfit : int, optional
      Order for background continuum

    Returns
    -------
    tampl, tcent, twid, w, yprep

    """
    from pypit import arcyarc
    fitp = 7 #slf._argflag['arc']['calibrate']['nfitpix']
    if len(censpec.shape) == 3: detns = censpec[:, 0].flatten()
    else: detns = censpec.copy()
    xrng = np.arange(float(detns.size))
    yrng = np.zeros(detns.size)
    mask = np.zeros(detns.size, dtype=np.int)
    mskcnt = 0
    # Continuum
    while True:
        w = np.where(mask == 0)
        xfit = xrng[w]
        yfit = detns[w]
        ct = np.polyfit(xfit, yfit, bpfit)
        yrng = np.polyval(ct, xrng)
        sigmed = 1.4826*np.median(np.abs(detns[w]-yrng[w]))
        w = np.where(detns > yrng+1.5*sigmed)
        mask[w] = 1
        if mskcnt == np.sum(mask):
            break  # No new values have been included in the mask
        mskcnt = np.sum(mask)
    #
    w = np.where(mask == 0)
    xfit = xrng[w]
    yprep = detns - yrng
    sfit = 1.4826*np.abs(detns[w]-yrng[w])
    ct = np.polyfit(xfit, sfit, bpfit)
    yerr = np.polyval(ct, xrng)
    myerr = np.median(np.sort(yerr)[:yerr.size/2])
    yerr[np.where(yerr < myerr)] = myerr
    # Find all significant detections
    # The last argument is the overall minimum significance level of an arc line detection and the second
    # last argument is the level required by an individual pixel before the neighbourhood of this pixel is searched.
    pixt = np.where((yprep/yerr>0.0) & (yprep>np.roll(yprep,1)) & (yprep>=np.roll(yprep,-1))
                                      & (np.roll(yprep,1)>np.roll(yprep,2)) & (np.roll(yprep,-1)>np.roll(yprep,-2))
                                      & (np.roll(yprep,2)>np.roll(yprep,3)) & (np.roll(yprep,-2)>np.roll(yprep,-3)))[0]
    tampl, tcent, twid, ngood = arcyarc.fit_arcorder(xrng, yprep, pixt, fitp)
    w = np.where((np.isnan(twid) == False) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent > 0.0) & (tcent < xrng[-1]))
    # Return
    return tampl, tcent, twid, w, yprep

