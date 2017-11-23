""" Generate a series of fake data and test how well
the calibration algorithm recovers the correct result
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import warnings
import pdb

from arclines.holy import grail
from astropy.table import Table


def tst_holy(spec, lines, solution, test='general', tol=0.01):
    from pypit import arutils

    # Run
    if test == 'semi_brute':
        print("Not implemented yet!")
    elif test == 'general':
        best_dict, final_fit = grail.general(spec, lines, islinelist=True)
    else:
        pdb.set_trace()

    # Generate the estimated solution
    xsol = np.linspace(-1.0, 1.0, spec.size)
    ysol = arutils.func_val(final_fit['fitc'], xsol, final_fit['function'],
                            minv=final_fit['fmin'], maxv=final_fit['fmax'])

    # Compare this with the true solution and ensure that the
    # maximum deviation is less than the allowed tolerance
    MAKE THIS A PIXEL DEVIATION BY MULTIPLYING BY THE DERIVATIVE?
    max_diff = np.max(np.abs((ysol-solution)/solution))
    if max_diff < tol:
        grade = 'PASSED'
    else:
        grade = 'FAILED'

    # Return the result
    return grade, best_dict, final_fit


def gen_fakedata(wavecen, disp, nonlinear, ndet=25, nlines=100, nspurious=5, rms=0.1, npixels=2048):
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
    nlncoeff = np.array([np.random.uniform(-nonlinear, +nonlinear), 0.0, 1.0])

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
    plotsol = True
    if plotsol:
        from matplotlib import pyplot as plt
        plt.plot(detlines[:ndet], linelist[idxlines[:ndet]], 'bx')
        plt.show()
    return detlines[srt], linelist, idxlines[srt], solution


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

    wavecen, disp, nonlinear = 5000.0, 1.0, 0.02
    npixels = 2048

    # Run it
    sv_grade = []  # for the end, just in case
    for i in range(nsample):
        # Generate a new set of fake data
        detlines, linelist, idxlines, solution = gen_fakedata(wavecen, disp, nonlinear, npixels=npixels)
        spec = gen_spectrum(detlines, npixels=npixels)
        lltable = gen_linelist(linelist)
        grade, best_dict, final_fit = tst_holy(spec, lltable, solution, test=test, tol=0.01)
        sv_grade.append(grade)

    # Report it
    print('==============================================================')
    for grade in sv_grade:
        print("{:s}".format(grade))


# Test
if __name__ == '__main__':
    #flg_tst = 1   # Run em all with semi-brute
    flg_tst = 2   # Run em all with general

    main(flg_tst)
