""" Generate a series of fake data and test how well
the calibration algorithm recovers the correct result
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import warnings
import pdb

from arclines.holy import grail
from astropy.table import Table


def tst_holy(spec, lines, solution, test='general', tol=0.1):
    """
    Run a test with a wavelength calibration algorithm

    Parameters
    ----------
    spec : ndarray
      Fake Spectrum
    lines : Table
      Astropy table containing the fake linelist
    solution : ndarray
      The true wavelength solution of every pixel
    test : str
      Which calibration algorithm should be tested
    tol : float
      maximum allowed pixel deviation of the wavelength solution

    Returns
    -------
    grade : bool
      PASS or FAIL
    best_dict : dict
      Initial ID solution
    final_fit : dict
      Final fit solution
    """
    from pypit import arutils

    # Run
    if test == 'semi_brute':
        print("Not implemented yet!")
    elif test == 'general':
        best_dict, final_fit = grail.general(spec, lines, islinelist=True, outroot="test")
    else:
        pdb.set_trace()

    # Generate the estimated solution
    xsol = np.linspace(0.0, 1.0, spec.size)
    ysol = arutils.func_val(final_fit['fitc'], xsol, final_fit['function'],
                            minv=final_fit['fmin'], maxv=final_fit['fmax'])
    ysolp = arutils.func_val(final_fit['fitc'], xsol+tol/spec.size, final_fit['function'],
                             minv=final_fit['fmin'], maxv=final_fit['fmax'])
    ysolm = arutils.func_val(final_fit['fitc'], xsol-tol/spec.size, final_fit['function'],
                             minv=final_fit['fmin'], maxv=final_fit['fmax'])

    # Compare this with the true solution and ensure that the
    # maximum deviation is less than the allowed tolerance
    if ysolp[0] > ysolm[0]:
        passed = np.all((ysolm < solution) & (ysolp > solution))
    else:
        passed = np.all((ysolp < solution) & (ysolm > solution))

    if not passed:
        from matplotlib import pyplot as plt
        plt.plot(xsol, (solution-ysol)/ysol, 'r-')
        plt.plot(xsol, (ysolp-ysol)/ysol, 'g-')
        plt.plot(xsol, (ysolm-ysol)/ysol, 'b-')
        plt.show()
        pdb.set_trace()

    # Return the result
    return passed, best_dict, final_fit


def gen_fakedata(wavecen, disp, nonlinear, ndet=40, nlines=100, nspurious=5, rms=0.1, npixels=2048):
    """
    Generate a fake linelist and fake data.

    Parameters
    ----------
    wavecen : float
      Central wavelength
    disp : float
      Spectral dispersion (Angstroms/pixel)
    nonlinear : float
      As a percentage, how non-linear is the pixel-wavelength solution
    ndet : int
      Number of detected lines used for the fake data
    nlines : int
      Number of lines to be used in the linelist
    nspurious : int
      Number of lines detected in the data but are not in the linelist
    rms : float
      Perturb each line detection by a fraction of a pixel
    npixels : int
      Number of pixels in the dispersion direction

    Returns
    -------
    detlines : ndarray
      Fake data (including spurious)
    linelist : ndarray
      Fake linelist (including spurious)
    idxlines : ndarray
      mask containing 0's (spurious) or the index in the linelist of each detection
    solution : ndarray
      The correct wavelength of every pixel
    """

    # Randomly choose whether pixels correlate or anticorrelate with wavelength
    sign = np.random.choice([-1.0, 1.0])

    # Setup a linear coefficient array
    wavcoeff = np.array([sign*disp*npixels/2.0, wavecen])

    # Generate random pixels on a detector
    pixlist = np.random.uniform(-1.0, 1.0, ndet)
    pixspur = np.random.uniform(-1.0, 1.0, nspurious)
    pixsoln = np.linspace(-1.0, 1.0, npixels)

    # Perturb the linear values onto a non-linear solution
    npert = 4
    while True:
        nlcff = np.random.uniform(-nonlinear, +nonlinear, npert-2)
        nlncoeff = np.append(nlcff, np.array([0.0, 1.0]))
        if np.all(np.abs(np.polyval(nlncoeff, pixsoln)-1.0) < nonlinear):
            break

    # Determine their wavelengths
    waves = np.polyval(wavcoeff, pixlist) * np.polyval(nlncoeff, pixlist)
    linelist = np.append(waves, np.random.uniform(3000.0, 10000.0, nlines))
    linelist.sort()
    solution = np.polyval(wavcoeff, pixsoln) * np.polyval(nlncoeff, pixsoln)

    # Get true list
    truwaves = np.polyval(wavcoeff, pixlist) * np.polyval(nlncoeff, pixlist)

    # Pick ndet pixels that are detected and add in nspurious
    detlines = npixels * (np.append(pixlist, pixspur) + 1.0) / 2.0
    idxlines = np.append(np.searchsorted(linelist, truwaves), np.zeros(nspurious, dtype=np.int))
    srt = np.argsort(detlines)
    detlines += np.random.normal(0.0, rms, detlines.size)

    # Plot the solution?
    plotsol = False
    if plotsol:
        from matplotlib import pyplot as plt
        plt.plot(detlines[:ndet], linelist[idxlines[:ndet]], 'bx')
        plt.show()
    return detlines[srt], linelist, idxlines[srt], solution, sign


def gen_linelist(lines):
    linetable = Table()
    linetable['ion'] = ['FAKELINE']*lines.size
    linetable['wave'] = lines
    linetable['NIST'] = [1]*lines.size
    linetable['Instr'] = [1]*lines.size
    linetable['amplitude'] = [1.0]*lines.size
    linetable['Source'] = ['FAKELINE']*lines.size
    return linetable


def gen_spectrum(detlines, npixels=2048):
    xspec = np.arange(npixels, dtype=np.float)
    yspec = np.zeros(npixels, dtype=np.float)
    ampl = np.random.uniform(0.0, 20000.0, detlines.size)
    sigma = 1.0  # pixels
    for i in range(detlines.size):
        yspec += ampl[i] * np.exp(-(xspec-detlines[i])**2/(2.0*sigma)**2)
    return np.random.normal(yspec, np.sqrt(yspec+100.0))


def main(flg_tst, nsample=1000):

    if flg_tst in [1]:
        test = 'semi_brute'
    elif flg_tst in [2]:
        test = 'general'
    else:
        print("Test not implemented")
        pdb.set_trace()

    # Note, in most cases, the deviation from linear is < 0.5% (i.e. nonlinear = 0.005)
    # so a one per cent deviation from linear (i.e. nonlinear = 0.01) is very reasonable.
    wavecen, disp, nonlinear = 5000.0, 1.0, 0.01
    npixels = 2048

    # Run it
    sv_grade = np.zeros(nsample, dtype=np.bool)  # for the end, just in case
    for i in range(nsample):
        # Generate a new set of fake data
        detlines, linelist, idxlines, solution, sign = gen_fakedata(wavecen, disp, nonlinear, npixels=npixels)
        spec = gen_spectrum(detlines, npixels=npixels)
        lltable = gen_linelist(linelist)
        passed, best_dict, final_fit = tst_holy(spec, lltable, solution, test=test, tol=1.0)
        sv_grade[i] = passed
        if not passed:
            print(sign)

    # Report results
    print('==============================================================')
    print(np.sum(sv_grade), "/", nsample)


# Test
if __name__ == '__main__':
    #flg_tst = 1   # Run em all with semi-brute
    flg_tst = 2   # Run em all with general

    main(flg_tst)
