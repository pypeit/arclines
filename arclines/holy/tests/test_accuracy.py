""" Generate a series of fake data and test how well
the calibration algorithm recovers the correct result
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import warnings
import pdb

from arclines.holy import grail

import arclines
test_arc_path = arclines.__path__[0]+'/data/test_arcs/'
outdir = 'TEST_SUITE_OUTPUT/'


def tst_holy(name, spec_file, lines, wv_cen, disp, score, fidx, test='semi_brute'):

    # Favored parameters (should match those in the defaults)
    siglev=20.
    min_ampl=1000.
    min_match = 10

    # Load spectrum
    exten = spec_file.split('.')[-1]
    if exten == 'json':
        with open(test_arc_path+spec_file,'r') as f:
            pypit_fit = json.load(f)
        spec = np.array(pypit_fit['spec'])
    elif exten == 'hdf5':
        hdf = h5py.File(test_arc_path+spec_file,'r')
        spec = hdf['arcs/{:d}/spec'.format(fidx)].value
    else:
        pdb.set_trace()

    # Run
    outroot = outdir+name
    if test == 'semi_brute':
        best_dict, final_fit = grail.semi_brute(spec, lines, wv_cen, disp, siglev=siglev,
                                                min_ampl=min_ampl, min_nmatch=10, outroot=outroot)
    elif test == 'general':
        best_dict, final_fit = grail.general(spec, lines, siglev=siglev,
                                             min_ampl=min_ampl, min_nmatch=10, outroot=outroot)
    else:
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

    # Run it
    sv_grade = [] # for the end, just in case
    for name,src_file,lines,wvcen,disp,score,fidx in zip(
            names,src_files,all_lines,all_wvcen,all_disp,scores,fidxs):
        #if '8500' not in name:
        #    continue
        grade, best_dict, final_fit = tst_holy(name, src_file, lines, wvcen, disp, score, fidx,
                                               test=test)
        sv_grade.append(grade)
        #if '900' in name:
        #    pdb.set_trace()

    # Report it
    print('==============================================================')
    for name,grade in zip(names,sv_grade):
        print("{:s} {:s}".format(name,grade))


# Test
if __name__ == '__main__':
    #flg_tst = 1   # Run em all with semi-brute
    flg_tst = 2   # Run em all with general

    main(flg_tst)
