""" Generate a series of fake data and test how well
the calibration algorithm recovers the correct result
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import warnings
import pdb

from arclines.holy import grail


def tst_holy(spec, lines, idx, npixels=2048, test='semi_brute'):

    # Run
    if test == 'semi_brute':
        print("Not implemented yet!")
    elif test == 'general':
        best_dict, final_fit = grail.general(spec, lines, npix=npixels, isspec=False, islinelist=True)
    else:
        pdb.set_trace()

    pdb.set_trace()

    # Score
    grade = 'PASSED'
    if final_fit['rms'] > score['rms']:
        grade = 'FAILED'
        warnings.warn("Solution for {:s} failed RMS!!".format(name))
    if len(final_fit['xfit']) < score['nxfit']:
        grade = 'FAILED'
        warnings.warn("Solution for {:s} failed N xfit!!".format(name))
    if best_dict['nmatch'] < score['nmatch']:
        grade = 'FAILED'
        warnings.warn("Solution for {:s} failed N match!!".format(name))

    # Warn
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
    detidx : ndarray
      mask containing 0's (spurious) or the index in the linelist of each detection
    """

    # Randomly choose whether pixels correlate or anticorrelate with wavelength
    sign = np.random.choice([-1.0, 1.0])

    # Setup a linear coefficient array
    wavcoeff = np.array([sign*disp*npixels/2.0, wavecen])

    # Generate random pixels on a detector
    pixlist = np.random.uniform(-1.0, 1.0, ndet)
    pixspur = np.random.uniform(-1.0, 1.0, nspurious)

    # Perturb the linear values onto a non-linear solution
    nlncoeff = np.array([np.random.uniform(-nonlinear, +nonlinear), 0.0, 1.0])

    # Determine their wavelengths
    waves = np.polyval(wavcoeff, pixlist) * np.polyval(nlncoeff, pixlist)
    linelist = np.append(waves, np.random.uniform(3000.0, 10000.0, nlines))
    linelist.sort()

    # Get true list
    truwaves = np.polyval(wavcoeff, pixlist) * np.polyval(nlncoeff, pixlist)

    # Pick ndet pixels that are detected and add in nspurious
    detlines = npixels * (np.append(pixlist, pixspur) + 1.0) / 2.0
    idxlines = np.append(np.searchsorted(linelist, truwaves), np.zeros(nspurious, dtype=np.int))
    srt = np.argsort(detlines)
    detlines += np.random.normal(0.0, rms, detlines.size)
    return detlines[srt], linelist, idxlines[srt]


def main(flg_tst, nsample=1000):

    if flg_tst in [1]:
        test = 'semi_brute'
    elif flg_tst in [2]:
        test = 'general'

    wavecen, disp, nonlinear = 5000.0, 1.0, 0.1
    npixels = 2048

    # Run it
    sv_grade = []  # for the end, just in case
    for i in range(nsample):
        # Generate a new set of fake data
        detlines, linelist, idxlines = gen_fakedata(wavecen, disp, nonlinear, npixels=npixels)
        grade, best_dict, final_fit = tst_holy(detlines, linelist, idxlines, npixels=npixels, test=test)
        sv_grade.append(grade)

    # Report it
    print('==============================================================')
    for name, grade in zip(names,sv_grade):
        print("{:s} {:s}".format(name,grade))


# Test
if __name__ == '__main__':
    #flg_tst = 1   # Run em all with semi-brute
    flg_tst = 2   # Run em all with general

    main(flg_tst)
