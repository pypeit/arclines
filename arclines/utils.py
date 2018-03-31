""" Utilities for building arcline lists
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import warnings
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


def robust_polyfit(xarray, yarray, order, weights=None, maxone=True, sigma=3.0,
                   function="polynomial", initialmask=None, forceimask=False,
                   minv=None, maxv=None, guesses=None, **kwargs):
    """ Taken from PYPIT
    A robust (equally weighted) polynomial fit is performed to the xarray, yarray pairs
    mask[i] = 1 are masked values

    :param xarray: independent variable values
    :param yarray: dependent variable values
    :param order: the order of the polynomial to be used in the fitting
    :param weights: weights to be used in the fitting (weights = 1/sigma)
    :param maxone: If True, only the most deviant point in a given iteration will be removed
    :param sigma: confidence interval for rejection
    :param function: which function should be used in the fitting (valid inputs: 'polynomial', 'legendre', 'chebyshev', 'bspline')
    :param initialmask: a mask can be supplied as input, these values will be masked for the first iteration. 1 = value masked
    :param forceimask: if True, the initialmask will be forced for all iterations
    :param minv: minimum value in the array (or the left limit for a legendre/chebyshev polynomial)
    :param maxv: maximum value in the array (or the right limit for a legendre/chebyshev polynomial)
    :return: mask, ct -- mask is an array of the masked values, ct is the coefficients of the robust polyfit.
    """
    # Setup the initial mask
    if initialmask is None:
        mask = np.zeros(xarray.size, dtype=np.int)
        if forceimask:
            warnings.warn("Initial mask cannot be enforced -- no initital mask supplied")
            forceimask = False
    else:
        mask = initialmask.copy()
    mskcnt = np.sum(mask)
    # Iterate, and mask out new values on each iteration
    ct = guesses
    while True:
        w = np.where(mask == 0)
        xfit = xarray[w]
        yfit = yarray[w]
        if weights is not None:
            wfit = weights[w]
        else:
            wfit = None
        ct = func_fit(xfit, yfit, function, order, w=wfit,
                      guesses=ct, minv=minv, maxv=maxv, **kwargs)
        yrng = func_val(ct, xarray, function, minv=minv, maxv=maxv)
        sigmed = 1.4826*np.median(np.abs(yfit-yrng[w]))
        if xarray.size-np.sum(mask) <= order+2:
            warnings.warn("More parameters than data points - fit might be undesirable")
            break  # More data was masked than allowed by order
        if maxone:  # Only remove the most deviant point
            tst = np.abs(yarray[w]-yrng[w])
            m = np.argmax(tst)
            if tst[m] > sigma*sigmed:
                mask[w[0][m]] = 1
        else:
            if forceimask:
                w = np.where((np.abs(yarray-yrng) > sigma*sigmed) | (initialmask==1))
            else:
                w = np.where(np.abs(yarray-yrng) > sigma*sigmed)
            mask[w] = 1
        if mskcnt == np.sum(mask): break  # No new values have been included in the mask
        mskcnt = np.sum(mask)
    # Final fit
    w = np.where(mask == 0)
    xfit = xarray[w]
    yfit = yarray[w]
    if weights is not None:
        wfit = weights[w]
    else:
        wfit = None
    ct = func_fit(xfit, yfit, function, order, w=wfit, minv=minv, maxv=maxv, **kwargs)
    return mask, ct


def calc_fit_rms(xfit, yfit, fit, func, minv=None, maxv=None):
    """ Simple RMS calculation

    Parameters
    ----------
    xfit : ndarray
    yfit : ndarray
    fit : coefficients
    func : str
    minv : float, optional
    maxv : float, optional

    Returns
    -------
    rms : float

    """
    values = func_val(fit, xfit, func, minv=minv, maxv=maxv)
    rms = np.std(yfit-values)
    # Return
    return rms


def func_fit(x, y, func, deg, minv=None, maxv=None, w=None, guesses=None,
             **kwargs):
    """ General routine to fit a function to a given set of x,y points

    Parameters
    ----------
    x : ndarray
    y : ndarray
    func : str
      polynomial, legendre, chebyshev, bspline, gaussian
    deg : int
      degree of the fit
    minv : float, optional
    maxv
    w
    guesses : tuple
    kwargs

    Returns
    -------
    coeff : ndarray or tuple
      ndarray for standard function fits
      tuple for bspline

    """
    if func == "polynomial":
        return np.polynomial.polynomial.polyfit(x, y, deg, w=w)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legfit(xv, y, deg, w=w)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebfit(xv, y, deg, w=w)
    else:
        raise IOError("Fitting function '{0:s}' is not implemented yet\n"+"Please choose from 'polynomial', 'legendre', 'chebyshev','bspline'")


def func_val(c, x, func, minv=None, maxv=None):
    """ Generic routine to return an evaluated function
    Functional forms include:
      polynomial, legendre, chebyshev, bspline, gauss

    Parameters
    ----------
    c : ndarray
      coefficients
    x
    func
    minv
    maxv

    Returns
    -------
    values : ndarray

    """
    if func == "polynomial":
        return np.polynomial.polynomial.polyval(x, c)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legval(xv, c)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebval(xv, c)
    else:
        raise ValueError("Fitting function '{0:s}' is not implemented yet\n"+"Please choose from 'polynomial', 'legendre', 'chebyshev', 'bspline'")
