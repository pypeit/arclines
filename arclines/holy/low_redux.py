""" Generate hdf5 files from LowRedux save files
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, json
import numpy as np
import pdb
import h5py
from scipy.io.idl import readsav

from matplotlib import pyplot as plt

from astropy.table import Table
from astropy import units as u

from arclines.holy.utils import find_peaks

from pypit import pyputils
msgs = pyputils.get_dummy_logger()

from pypit import ararclines
from pypit import arwave
from pypit import arutils
arutils.dummy_settings()


try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def fcheby(xnrm,order):
    leg = np.zeros((len(xnrm),order))
    leg[:,0] = 1.
    if order >= 2:
        leg[:,1] = xnrm
    # For looop
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
    import arclines
    out_path = arclines.__path__[0]+'/data/test_arcs/'
    # Read IDL save file
    sav_file = os.getenv('LONGSLIT_DIR')+'calib/linelists/'+sav_file
    s = readsav(sav_file)
    ctbl = Table(s['calib'])  # For writing later

    # Line list
    alist = ararclines.load_arcline_list(None, None, lamps, None)

    # Meta data
    mdict = dict(npix=len(s['archive_arc'][0]), instr=instr,
                 lamps=[str(ilamp) for ilamp in lamps],  # For writing to hdf5
                 nspec=len(s['archive_arc']), infil=sav_file, IDairvac='vac')
    print("Processing {:d} spectra in {:s}".format(mdict['nspec'], sav_file))

    # Start output
    outh5 = h5py.File(out_path+outfil, 'w')
    outh5.create_group('arcs')

    # Loop on spectra
    for ss in range(mdict['nspec']):
        sss = str(ss)
        # Spec
        spec = s['archive_arc'][ss]
        # Peaks
        tampl, tcent, twid, w, yprep = find_peaks(spec)
        pixpk = tcent[w]
        pixampl = tampl[w]

        # Wavelength solution
        wv_air = cheby_val(s['calib'][ss]['ffit'], np.arange(mdict['npix']),
                   s['calib'][ss]['nrm'],s['calib'][ss]['nord'])
        # Check blue->red or vice-versa
        if ss == 0:
            if wv_air[0] > wv_air[-1]:
                mdict['bluered'] = False
            else:
                mdict['bluered'] = True

        # Peak waves
        twave_air = cheby_val(s['calib'][ss]['ffit'], pixpk,
                              s['calib'][ss]['nrm'],s['calib'][ss]['nord'])
        # Air to Vac
        twave_vac = arwave.airtovac(twave_air*u.AA)
        wave_vac = arwave.airtovac(wv_air*u.AA)
        # IDs
        idwv = np.zeros_like(pixpk)
        idsion = np.array([str('12345')]*len(pixpk))
        for kk,twv in enumerate(twave_vac.value):
            # diff
            diff = np.abs(twv-alist['wave'])
            if np.min(diff) < dtoler:
                imin  = np.argmin(diff)
                idwv[kk] = alist['wave'][imin]
                idsion[kk] = alist['Ion'][imin]
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
            debugger.set_trace()
    # Close
    outh5.close()
    print('Wrote {:s}'.format(out_path+outfil))


# Command line execution
def main(flg_tst):

    # LRISb 600
    if (flg_tst % 2**1) >= 2**0:
        generate_hdf('lris_blue_600.sav', 'LRISb_600', ['ZnI', 'CdI', 'HgI', 'NeI', 'ArI'], 'LRISb_600.hdf5')
    # LRISr 600
    if (flg_tst % 2**2) >= 2**1:
        generate_hdf('lris_red_600_7500.sav', 'LRISr_600', ['ArI', 'HgI', 'KrI', 'NeI', 'XeI'], 'LRISr_600_7500.hdf5')
    # Kastb 600
    if (flg_tst % 2**3) >= 2**2:
        generate_hdf('kast_600_4310.sav', 'Kastb_600', ['ArI', 'HgI', 'KrI', 'NeI', 'XeI'], 'LRISr_600_7500.hdf5')

# Test
if __name__ == '__main__':
    flg_tst = 0
    #flg_tst += 2**0   # LRISb 600
    #flg_tst += 2**1   # LRISr 600
    flg_tst += 2**2   # Kastb 600

    main(flg_tst)
