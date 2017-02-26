""" Module for finding patterns in arc line spectra
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import pdb

from arclines import io as arcl_io
from arclines.holy import patterns as arch_patt
from arclines.holy import fitting as arch_fit
from arclines.holy import utils as arch_utils
from arclines.pypit_utils import find_peaks


def basic(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
          swv_uncertainty=350., pix_tol=2, plot_fil=None, min_match=5,
          **kwargs):
    """ Basic holy grail algorithm

    Parameters
    ----------
    spec : spectrum
    lines : list
      List of arc lamps on
    wv_cen : float
      Guess at central wavelength
    disp : float
      Dispersion A/pix
    siglev
    min_ampl
    swv_uncertainty
    pix_tol
    plot_fil

    Returns
    -------
    status : int

    """

    # Init line-lists and wavelength 'guess'
    npix = spec.size
    wave = wv_cen + (np.arange(npix) - npix/2.)*disp

    line_lists = arcl_io.load_line_lists(lines, unknown=True)
    wvdata = line_lists['wave'].data  # NIST + Extra
    isrt = np.argsort(wvdata)
    wvdata = wvdata[isrt]

    # Find peaks
    all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, siglev=siglev, min_ampl=min_ampl)

    # Matching
    match_idx, scores = arch_patt.run_quad_match(cut_tcent, wave, wvdata,
                                                 disp, swv_uncertainty=swv_uncertainty,
                                                 pix_tol=pix_tol)

    # Check quadrants
    xquad = npix//4 + 1
    print("================================================================")
    print("Checking quadrants:")
    print("----------------------------------------------------------------")
    for jj in range(4):
        tc_in_q = (cut_tcent >= jj*xquad) & (cut_tcent < (jj+1)*xquad)
        cstat = 'quad {:d}: ndet={:d}'.format(jj, np.sum(tc_in_q))
        # Stats
        for key in ['Perf', 'Good', 'OK', 'Amb']:
            in_stat = scores[tc_in_q] == key
            cstat += ' {:s}={:d}'.format(key, np.sum(in_stat))
        # Print
        print(cstat)
    print("----------------------------------------------------------------")

    # Go for it!?
    mask = np.array([False]*len(all_tcent))
    IDs = []
    for kk,score in enumerate(scores):
        if score in ['Perf', 'Good', 'Ok']:
            mask[icut[kk]] = True
            uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
            imx = np.argmax(counts)
            IDs.append(wvdata[uni[imx]])
    ngd_match = np.sum(mask)
    if ngd_match < min_match:
        print("Insufficient matches to continue")
        status = -1
        return status, ngd_match, match_idx, scores, None

    # Fit
    NIST_lines = line_lists['NIST'] > 0
    ifit = np.where(mask)[0]
    final_fit = arch_fit.iterative_fitting(spec, all_tcent, ifit,
                               IDs, line_lists[NIST_lines], disp, plot_fil=plot_fil)
    # Return
    status = 1
    return status, ngd_match, match_idx, scores, final_fit


def semi_brute(spec, lines, wv_cen, disp, siglev=20., min_ampl=300.,
               outroot=None, debug=False, do_fit=True, verbose=False,
               fit_parm=None, min_nmatch=0, lowest_ampl=200.):
    """
    Parameters
    ----------
    spec
    lines
    wv_cen
    disp
    siglev
    min_ampl
    outroot
    debug
    do_fit
    verbose
    fit_parm
    min_nmatch
    lowest_ampl

    Returns
    -------
    best_dict : dict
    final_fit : dict

    """
    # imports
    from astropy.table import vstack
    from linetools import utils as ltu
    from arclines import plots as arcl_plots
    # Load line lists
    line_lists = arcl_io.load_line_lists(lines)
    unknwns = arcl_io.load_unknown_list(lines)

    npix = spec.size

    # Lines
    all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=min_ampl)

    # Best
    best_dict = dict(nmatch=0, ibest=-1, bwv=0., min_ampl=min_ampl)

    # 3 things to fiddle:
    #  pix_tol -- higher for fewer lines  1/2
    #  unknowns -- on for fewer lines  off/on
    #  scoring -- weaken for more lines ??

    # Loop on unknowns
    for unknown in [False, True]:
        if unknown:
            tot_list = vstack([line_lists,unknwns])
        else:
            tot_list = line_lists
        wvdata = np.array(tot_list['wave'].data) # Removes mask if any
        wvdata.sort()
        sav_nmatch = best_dict['nmatch']

        # Loop on pix_tol
        for pix_tol in [1.,2.]:
            # Scan on wavelengths
            arch_patt.scan_for_matches(wv_cen, disp, npix, cut_tcent, wvdata,
                                      best_dict=best_dict, pix_tol=pix_tol)
            # Lower minimum amplitude
            ampl = min_ampl
            while(best_dict['nmatch'] < min_nmatch):
                ampl /= 2.
                if ampl < lowest_ampl:
                    break
                all_tcent, cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=ampl)
                arch_patt.scan_for_matches(wv_cen, disp, npix, cut_tcent, wvdata,
                                       best_dict=best_dict, pix_tol=pix_tol, ampl=ampl)

        # Save linelist?
        if best_dict['nmatch'] > sav_nmatch:
            best_dict['line_list'] = tot_list
            best_dict['unknown'] = unknown
            best_dict['ampl'] = unknown

    if best_dict['nmatch'] == 0:
        print('---------------------------------------------------')
        print('Report:')
        print('::   No matches!  Could be you input a bad wvcen or disp value')
        print('---------------------------------------------------')
        return

    # Report
    print('---------------------------------------------------')
    print('Report:')
    print('::   Number of lines recovered = {:d}'.format(all_tcent.size))
    print('::   Number of lines analyzed = {:d}'.format(cut_tcent.size))
    print('::   Number of Perf/Good/Ok matches = {:d}'.format(best_dict['nmatch']))
    print('::   Best central wavelength = {:g}A'.format(best_dict['bwv']))
    print('::   Best solution used pix_tol = {}'.format(best_dict['pix_tol']))
    print('::   Best solution had unknown = {}'.format(best_dict['unknown']))
    print('---------------------------------------------------')

    if debug:
        match_idx = best_dict['midx']
        for kk in match_idx.keys():
            uni, counts = np.unique(match_idx[kk]['matches'], return_counts=True)
            print('kk={}, {}, {}, {}'.format(kk, uni, counts, np.sum(counts)))

    # Write scores
    #out_dict = best_dict['scores']
    #jdict = ltu.jsonify(out_dict)
    #ltu.savejson(pargs.outroot+'.scores', jdict, easy_to_read=True, overwrite=True)
    #print("Wrote: {:s}".format(pargs.outroot+'.scores'))

    # Write IDs
    if outroot is not None:
        out_dict = dict(pix=cut_tcent, IDs=best_dict['IDs'])
        jdict = ltu.jsonify(out_dict)
        ltu.savejson(outroot+'.json', jdict, easy_to_read=True, overwrite=True)
        print("Wrote: {:s}".format(outroot+'.json'))

    # Plot
    if outroot is not None:
        arcl_plots.match_qa(spec, cut_tcent, best_dict['line_list'],
                            best_dict['IDs'], best_dict['scores'], outroot+'.pdf')
        print("Wrote: {:s}".format(outroot+'.pdf'))

    # Fit
    final_fit = None
    if do_fit:
        #
        NIST_lines = line_lists['NIST'] > 0
        ifit = np.where(best_dict['mask'])[0]
        if outroot is not None:
            plot_fil = outroot+'_fit.pdf'
        else:
            plot_fil = None
        # Purge UNKNOWNS from ifit
        imsk = np.array([True]*len(ifit))
        for kk, idwv in enumerate(np.array(best_dict['IDs'])[ifit]):
            if np.min(np.abs(line_lists['wave'][NIST_lines]-idwv)) > 0.01:
                imsk[kk] = False
        ifit = ifit[imsk]
        # Allow for weaker lines in the fit
        all_tcent, weak_cut_tcent, icut = arch_utils.arc_lines_from_spec(spec, min_ampl=lowest_ampl)
        add_weak = []
        for weak in weak_cut_tcent:
            if np.min(np.abs(cut_tcent-weak)) > 5.:
                add_weak += [weak]
        if len(add_weak) > 0:
            cut_tcent = np.concatenate([cut_tcent, np.array(add_weak)])
        # Fit
        final_fit = arch_fit.iterative_fitting(spec, cut_tcent, ifit,
                                               np.array(best_dict['IDs'])[ifit], line_lists[NIST_lines],
                                               disp, plot_fil=plot_fil, verbose=verbose, aparm=fit_parm)
        if plot_fil is not None:
            print("Wrote: {:s}".format(plot_fil))

    # Return
    return best_dict, final_fit
