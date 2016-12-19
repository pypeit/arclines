""" Module for loading source files
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import json
import pdb

from astropy.table import Table

import arclines # For path
src_path = arclines.__path__[0]+'/data/sources/'

from arclines import plots as arcl_plots

# Hard-coded string lengths, and more
from arclines import defs
str_len_dict = defs.str_len()
instr_dict = defs.instruments()
line_dict = defs.lines()

def load(src_file, format, **kwargs):
    """
    Parameters
    ----------
    src_file
    format

    Returns
    -------

    """
    # Load
    if format == 'PYPIT1':
        ID_lines, U_lines = load_pypit(1, src_file, **kwargs)
    elif format == 'LRDX1':
        ID_lines, U_lines = load_low_redux(1, src_file, **kwargs)
    else:
        raise IOError("Format {:s} for source {:s} is not supported".format(
                format, src_file))
    # Reject lines

    # Return
    return ID_lines, U_lines


def load_pypit(version, src_file, plot=False, **kwargs):
    """ Load from PYPIT output

    Parameters
    ----------
    version : int
      Flag indicating version of PYPIT output
      1 : JSON format from PYPIT v1
    src_file : str
    plot : bool, optional
      Generate a plot?

    Returns
    -------
    ID_lines : Table
      Table of arc lines with IDs
    U_lines : Table or None
      Additional lines

    """
    from pypit import arutils
    # Load
    if version != 1:
        raise IOError("Unimplemented version!")
    with open(src_path+src_file,'r') as f:
        pypit_fit = json.load(f)

    npix = len(pypit_fit['spec'])
    # ID lines -- Assumed in NIST if from PYPIT
    ions = np.array(pypit_fit['ions'],
                    dtype='S{:d}'.format(str_len_dict['ion']))
    ID_lines = Table()
    ID_lines['ion'] = ions
    ID_lines['wave'] = pypit_fit['yfit']
    ID_lines['NIST'] = 1
    amps = []
    for jj,xfit in enumerate(pypit_fit['xfit']):
        pix = int(np.round(xfit*(npix-1)))
        amps.append(int(pypit_fit['spec'][pix]))
    ID_lines['amplitude'] = amps

    mn_ID, mx_ID = min(ID_lines['wave']), max(ID_lines['wave'])
    # Unknown
    # Use PYPIT to decode
    wave = arutils.func_val(pypit_fit['fitc'],
                            np.array(pypit_fit['tcent'])/(npix-1),
                            pypit_fit['function'],
                            minv=pypit_fit['fmin'], maxv=pypit_fit['fmax'])
    eamps, extras, epix = [], [], []
    for kk,iwave in enumerate(wave):
        if (np.min(np.abs(ID_lines['wave'].data-iwave)) > 0.5) and (
                iwave > mn_ID) and (iwave < mx_ID):  # NO EXTRAPOLATION
            extras.append(iwave)
            pix = int(np.round(pypit_fit['tcent'][kk]))
            epix.append(pix) # For plotting
            eamps.append(int(pypit_fit['spec'][pix]))
    if len(extras) == 0:
        U_lines = None
    else:
        U_lines = Table()
        U_lines['wave'] = extras
        U_lines['ion'] = str('UNKNWN').rjust(str_len_dict['ion'])
        U_lines['NIST'] = 0
        U_lines['amplitude'] = eamps

    # Plot
    if plot:
        # Generate IDs
        IDs = []
        for row in ID_lines:
            IDs.append('{:s} {:.4f}'.format(row['ion'], row['wave']))
        # Extras
        if U_lines is not None:
            pextras = dict(x=epix, IDs=[])
            for row in U_lines:
                pextras['IDs'].append('{:s} {:.4f}'.format(row['ion'], row['wave']))
        else:
            pextras = None
        #
        arcl_plots.arc_ids(np.array(pypit_fit['spec']),
                           np.array(pypit_fit['xfit'])*(npix-1),
                           IDs, src_file.replace('.json', '.pdf'),
                           title=src_file.replace('.json', ''),
                           extras=pextras)

    # Return
    return ID_lines, U_lines


def load_low_redux(version, src_file, plot=False, min_hist=10,
                   cut_amp_val=400., wvmnx=[0., 1e9], ions=None):
    """
    Parameters
    ----------
    version : int
      1 : Find extras from set of LowRedux fits
    src_file : str
    plot

    Returns
    -------

    """
    if version != 1:
        raise IOError("Unimplemented version!")

    import h5py
    from scipy.interpolate import interp1d
    from arclines.pypit_utils import find_peaks
    from arclines.io import load_line_lists

    # Load existing line lists
    line_list = load_line_lists(ions)
    wvdata = line_list['wave'].data

    # Open
    hdf = h5py.File(src_path+src_file,'r')
    mdict = {}
    for key in hdf['meta'].keys():
        mdict[key] = hdf['meta'][key].value

    # Loop on spec
    extras = []
    eamps = []
    for ispec in range(mdict['nspec']):
        spec = hdf['arcs/'+str(ispec)+'/spec'].value
        wave = hdf['arcs/'+str(ispec)+'/wave'].value # vacuum
        disp = np.median(np.abs(wave-np.roll(wave,1)))
        npix = wave.size
        # Find peaks for extras
        tampl, tcent, twid, w, yprep = find_peaks(spec, siglev=10.)
        all_tcent = tcent[w]

        # Function for more precise wavelengths
        fwv = interp1d(np.arange(npix), wave)#, kind='cubic')
        #
        amps = []
        for itc in all_tcent:
            pix = int(np.round(itc))
            amps.append(np.max(spec[pix-1:pix+2]))
        amps = np.array(amps)

        # Trim tcent on amplitude
        cut_amp = amps > cut_amp_val  # 500.
        tcent = all_tcent[cut_amp]
        nlin = tcent.size

        # init with Truth
        for ii in range(nlin):
            wvt = float(fwv(tcent[ii]))
            mtw = np.where(np.abs(wvdata-wvt) < 4*disp)[0]  # Deals with bad wavelength solutions
            if len(mtw) == 0:
                if (wvt > wvmnx[0]) & (wvt < wvmnx[1]): # LRISb only
                    extras.append(wvt)
                    eamps.append(amps[ii])
    # Repackage
    extras = np.array(extras)
    isort = np.argsort(extras)
    extras = extras[isort]
    eamps = np.array(eamps)[isort]
    # Group -- Super-crude friends of friends
    final_extras = []
    final_amps = []
    cnt = 0
    while cnt <= extras.size:
        # First try
        ingroup = np.abs(extras-extras[cnt]) < 2*disp
        for ii in range(3):  # For some convergence
            # Median
            mngroup = np.median(extras[ingroup])
            ingroup = np.abs(extras-mngroup) < 2*disp
        if np.sum(ingroup) > min_hist:
            final_extras.append(np.median(extras[ingroup]))
            final_amps.append(np.median(eamps[ingroup]))
        # Step
        newe = mngroup + 5*disp
        gdcnt = np.where(extras > newe)[0]
        if len(gdcnt) == 0:
            break
        else:
            cnt = gdcnt[0]

    # Plot??
    if plot:
        # Find the best spectrum
        max_nex = 0
        for ispec in range(mdict['nspec']):
            spec = hdf['arcs/'+str(ispec)+'/spec'].value
            wave = hdf['arcs/'+str(ispec)+'/wave'].value # vacuum
            minwv, maxwv = np.min(wave), np.max(wave)
            nex = np.sum((final_extras>minwv) & (final_extras<maxwv))
            if nex > max_nex:
                svi = ispec
                max_nex = nex
        # Find pixel values
        spec = hdf['arcs/'+str(svi)+'/spec'].value
        wave = hdf['arcs/'+str(svi)+'/wave'].value # vacuum
        npix = wave.size
        fpix = interp1d(wave, np.arange(npix))#, kind='cubic')
        pextras = dict(x=fpix(final_extras), IDs=[])
        for fex in final_extras:
            pextras['IDs'].append('UNKNWN {:.4f}'.format(fex))
        # Plot
        arcl_plots.arc_ids(spec, [], [],
                           src_file.replace('.hdf5', '.pdf'),
                           title=src_file.replace('.hdf5', ''),
                           extras=pextras)
    # Table
    U_lines = Table()
    U_lines['wave'] = final_extras
    U_lines['ion'] = str('UNKNWN').rjust(str_len_dict['ion'])
    U_lines['NIST'] = 0
    U_lines['amplitude'] = final_amps

    # Return
    return None, U_lines


