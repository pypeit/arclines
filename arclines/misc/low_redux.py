""" Generate hdf5 files from LowRedux save files
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import numpy as np
import pdb
import h5py
from scipy.io.idl import readsav

from astropy.table import Table
from astropy import units as u

import arclines
out_path = arclines.__path__[0]+'/data/test_arcs/'


def fcheby(xnrm,order):
    leg = np.zeros((len(xnrm),order))
    leg[:,0] = 1.
    if order >= 2:
        leg[:,1] = xnrm
    # For loop
    for j in range(2,order):
        leg[:,j] = 2.0 * xnrm * leg[:,j-1] - leg[:,j-2]
    # Return
    return leg


def cheby_val(coeff, x, nrm, order):
    #
    xnrm = 2. * (x - nrm[0])/nrm[1]
    # Matrix first
    leg = fcheby(xnrm, order)
    # Dot
    return np.dot(leg, coeff)


def poly_val(coeff, x, nrm):
    #
    xnrm = 2. * (x - nrm[0])/nrm[1]
    #
    n = len(coeff)-1
    y = coeff[n]
    #for i=n-1,0,-1 do y = TEMPORARY(y) * x + c[i]
    for ii in range(n-1,-1,-1):
        y = y*xnrm + coeff[ii]
    return y


def generate_hdf(sav_file, instr, lamps, outfil, dtoler=0.6):
    """ Given an input LR IDL save file, generate an hdf5
    IDs arc lines too

    Parameters
    ----------
    sav_file : str
      Root name of the IDL save file from LowRedux, e.g. lris_blue_600.sav
    lamps
    outfil

    Returns
    -------

    """
    from pypit import pyputils
    msgs = pyputils.get_dummy_logger()

    from pypit import arwave
    from pypit import arutils
    arutils.dummy_settings()
    #
    from arclines.pypit_utils import find_peaks
    from arclines.io import load_line_lists
    #

    # Read IDL save file
    sav_file = os.getenv('LONGSLIT_DIR')+'calib/linelists/'+sav_file
    s = readsav(sav_file)
    ctbl = Table(s['calib'])  # For writing later

    # Line list
    alist = load_line_lists(lamps)

    # One spectrum?
    ashape = s['archive_arc'].shape
    if len(ashape) == 1:
        nspec = 1
        npix = ashape[0]
    else:
        nspec = s['archive_arc'].shape[0]
        npix = ashape[1]

    # Meta data
    mdict = dict(npix=npix, instr=instr,
                 lamps=[str(ilamp) for ilamp in lamps],  # For writing to hdf5
                 nspec=nspec, infil=sav_file, IDairvac='vac')
    print("Processing {:d} spectra in {:s}".format(mdict['nspec'], sav_file))

    # Start output
    outh5 = h5py.File(out_path+outfil, 'w')
    outh5.create_group('arcs')

    # Loop on spectra
    for ss in range(mdict['nspec']):
        sss = str(ss)
        # Parse
        if nspec == 1:
            spec = s['archive_arc']
        else:
            spec = s['archive_arc'][ss]
        calib = s['calib'][ss]
        # Peaks
        tampl, tcent, twid, w, yprep = find_peaks(spec)
        pixpk = tcent[w]
        pixampl = tampl[w]

        # Wavelength solution
        if calib['func'] == 'CHEBY':
            wv_air = cheby_val(calib['ffit'], np.arange(mdict['npix']),
                       calib['nrm'], calib['nord'])
        elif calib['func'] == 'POLY':
            wv_air = poly_val(calib['ffit'], np.arange(mdict['npix']),
                               calib['nrm'])
        # Check blue->red or vice-versa
        if ss == 0:
            if wv_air[0] > wv_air[-1]:
                mdict['bluered'] = False
            else:
                mdict['bluered'] = True

        # Peak waves
        if calib['func'] == 'CHEBY':
            twave_air = cheby_val(calib['ffit'], pixpk,
                              calib['nrm'], calib['nord'])
        else:
            twave_air = poly_val(calib['ffit'], pixpk, calib['nrm'])
        # Air to Vac
        twave_vac = arwave.airtovac(twave_air*u.AA)
        wave_vac = arwave.airtovac(wv_air*u.AA)
        if ss == 0:
            disp = np.median(np.abs(wave_vac-np.roll(wave_vac,1)))
            print("Average dispersion = {:g}".format(disp))
        # IDs
        idwv = np.zeros_like(pixpk)
        idsion = np.array([str('12345')]*len(pixpk))
        for kk,twv in enumerate(twave_vac.value):
            # diff
            diff = np.abs(twv-alist['wave'])
            if np.min(diff) < dtoler:
                imin  = np.argmin(diff)
                idwv[kk] = alist['wave'][imin]
                #idsion[kk] = alist['Ion'][imin]  NIST
                idsion[kk] = alist['ion'][imin]
        # Red to blue?
        if mdict['bluered'] is False:
            pixpk = mdict['npix']-1 - pixpk
            # Re-sort
            asrt = np.argsort(pixpk)
            pixpk = pixpk[asrt]
            idwv = idwv[asrt]
            # Reverse
            spec = spec[::-1]
            wave_vac = wave_vac[::-1]
        # Output
        outh5['arcs'].create_group(sss)
        # Datasets
        outh5['arcs'][sss]['wave'] = wave_vac
        outh5['arcs'][sss]['wave'].attrs['airvac'] = 'vac'
        outh5['arcs'][sss]['spec'] = spec
        outh5['arcs'][sss]['spec'].attrs['flux'] = 'counts'
        outh5['arcs'][sss]['pixpk'] = pixpk
        outh5['arcs'][sss]['ID'] = idwv
        outh5['arcs'][sss]['ID'].attrs['airvac'] = 'vac'
        outh5['arcs'][sss]['Ion'] = idsion
        # LR wavelengths
        outh5['arcs'][sss]['LR_wave'] = wv_air
        outh5['arcs'][sss]['LR_wave'].attrs['airvac'] = 'air'
        # LR Fit
        outh5['arcs'][sss].create_group('LR_fit')
        for key in ctbl.keys():
            outh5['arcs'][sss]['LR_fit'][key] = ctbl[ss][key]

    # Meta data
    outh5.create_group('meta')
    for key in mdict.keys():
        try:
            outh5['meta'][key] = mdict[key]
        except TypeError:  # Probably a unicode thing
            pdb.set_trace()
    # Close
    outh5.close()
    print('Wrote {:s}'.format(out_path+outfil))


# Command line execution
def main(flg_tst):

    # LRISb 600
    if (flg_tst % 2**1) >= 2**0:
        generate_hdf('lris_blue_600.sav', 'LRISb_600',
                     ['ZnI', 'CdI', 'HgI', 'NeI', 'ArI'],
                     'LRISb_600_LRX.hdf5')
    # LRISr 600
    if (flg_tst % 2**2) >= 2**1:
        generate_hdf('lris_red_600_7500.sav', 'LRISr_600', ['ArI', 'HgI', 'KrI', 'NeI', 'XeI'], 'LRISr_600_7500.hdf5')

    # Kastb 600
    if (flg_tst % 2**3) >= 2**2:
        generate_hdf('kast_600_4310.sav', 'Kastb_600', ['ArI', 'HgI', 'KrI', 'NeI', 'XeI'], 'LRISr_600_7500.hdf5')

    # MMT RCS
    if (flg_tst % 2**4) >= 2**3:
        generate_hdf('mmt_rcs_600_6310.sav', 'MMT_RCS', ['ArI', 'NeI'], 'MMT_RCS_600_6310.hdf5')

    # MODS
    if (flg_tst % 2**5) >= 2**4:
        generate_hdf('mods_blue_400ms.sav', 'MODSb', ['XeI', 'KrI'], 'MODS_blue_400.hdf5')
        generate_hdf('mods_red_670.sav', 'MODSr', ['NeI', 'ArI'], 'MODS_red_670.hdf5')


# Test
if __name__ == '__main__':
    flg_tst = 0
    #flg_tst += 2**0   # LRISb 600
    #flg_tst += 2**1   # LRISr 600
    #flg_tst += 2**2   # Kastb 600
    #flg_tst += 2**3   # MMT RCS 600_6310
    flg_tst += 2**4   # MODS

    main(flg_tst)
