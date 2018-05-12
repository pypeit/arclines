""" Module for arc line plots
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import pdb

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

from arclines import utils as arcl_utils

# Default Path
import arclines
plot_path = arclines.__path__[0]+'/data/plots/'


def show_source(src_dict, line_lists, outfile, title=None, path=None, clobber=False,
                min_unk_ampl=0.):
    """ Plot of an input source for arclines

    Color code:
      blue = Good line
      gray = previously used unknown
      orange = new unknown
      red = Unknown in good line list (this shouldn't happen)
      green = ??
      brown = Unused UNKNOWN

    Parameters
    ----------
    line_lists : Table
      Existing line lists with the lamps of interest
    outfile : str
      Name of output file
    min_unk_ampl : float (optional)
      Cut on Amplitude for UNKNOWNs
    """
    # Path
    if path is None:
        path = plot_path
    iext = outfile.rfind('.')

    # Begin
    if os.path.isfile(plot_path+outfile) and (not clobber):
        print("Plot {:s} exists.  Remove if you wish to remake it".format(outfile))
        return

    # Parse input
    arc_spec = src_dict['spec']

    # In line_list?
    def chk_line_list(wave):
        mtw = np.where(np.abs(line_lists['wave']-wave) < 1e-3)[0]
        if len(mtw) == 0: # Not in the line list
            return 1
        elif len(mtw) != 1:
            pdb.set_trace()
        # Source used?
        if outfile[:iext] in line_lists['Source'][mtw[0]]:
            return 0
        else:
            return 2

    # IDs
    clrs = ['green', 'red', 'blue']  # Used, not used
    if src_dict['xIDs'] is not None:
        xIDs = src_dict['xIDs']
        IDlbls, IDclrs = [], []
        for row in src_dict['ID_lines']:
            IDlbls.append('{:s} {:.4f}'.format(row['ion'], row['wave']))
            # Color
            in_llist = chk_line_list(row['wave'])
            IDclrs.append(clrs[in_llist])
    else:
        xIDs = []

    # Unknowns
    U_lines = src_dict['U_lines']
    clrs = ['red', 'orange', 'green']  # Used, not used, ??
    if U_lines is not None:
        # Match to NIST
        mask, wv_match = arcl_utils.vette_unkwn_against_lists(U_lines, src_dict['uions'])
        xepix = src_dict['epix']
        extras = dict(x=[], IDs=[], clrs=[])
        for ss,row in enumerate(U_lines):
            # Check against minimum amplitude
            if row['amplitude'] < min_unk_ampl:
                continue
            # x
            extras['x'].append(xepix[ss])
            # Lbl
            if mask[ss] == 2:  # Matched to NIST
                lbl = '{:.4f}'.format(row['wave']) + ' [{:s}]'.format(wv_match[ss])
            else:
                lbl = 'UNKNWN {:.4f}'.format(row['wave'])
            extras['IDs'].append(lbl)
            # Color me
            if mask[ss] == 0: # In current line list
                extras['clrs'].append('gray')
            elif mask[ss] == -1: # Assuming Unknowns skipped
                extras['clrs'].append('brown')
            else:
                in_llist = chk_line_list(row['wave'])
                extras['clrs'].append(clrs[in_llist])
    else:
        extras = None

    # Plot
    pp = PdfPages(plot_path+outfile)
    plt.figure(figsize=(11, 8.5))
    plt.clf()
    gs = gridspec.GridSpec(2, 1)
    idfont = 'small'

    # Simple spectrum plot
    for qq in range(2):
        ax_spec = plt.subplot(gs[qq])
        ax_spec.plot(np.arange(len(arc_spec)), arc_spec, 'k')
        ymin, ymax = 0., np.max(arc_spec)
        ysep = ymax*0.03
        mn_yline = 1e9
        # Standard IDs
        for kk, x in enumerate(xIDs):
            yline = np.max(arc_spec[int(x)-2:int(x)+2])
            mn_yline = min(mn_yline, yline)
            # Tick mark
            ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], '-', color=IDclrs[kk])
            # label
            ax_spec.text(x, yline+ysep*1.3, '{:s}'.format(IDlbls[kk]), ha='center', va='bottom',
                size=idfont, rotation=90., color=IDclrs[kk])
        # Extras?
        if extras is not None:
            for kk, x in enumerate(extras['x']):
                yline = np.max(arc_spec[int(x)-2:int(x)+2])
                mn_yline = min(mn_yline, yline)
                # Tick mark
                ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], '-', color=extras['clrs'][kk])
                # label
                ax_spec.text(x, yline+ysep*1.3, '{:s}'.format(extras['IDs'][kk]), ha='center', va='bottom',
                    size=idfont, rotation=90., color=extras['clrs'][kk])
        # Axes
        ax_spec.set_xlim(0., len(arc_spec))
        if qq==1:
            ax_spec.set_yscale("log", nonposy='clip')
            ax_spec.set_ylim(mn_yline/10., 5*ymax)
        else:
            ax_spec.set_ylim(ymin, ymax*1.3)
        if qq == 0:
            ax_spec.set_xlabel('Pixel')
        ax_spec.minorticks_on()
        ax_spec.set_ylabel('Counts')
        if title is not None:
            ax_spec.text(0.04, 0.93, title, transform=ax_spec.transAxes,
                         size='x-large', ha='left')#, bbox={'facecolor':'white'})
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    pp.savefig(bbox_inches='tight')
    pp.close()
    plt.close()
    print("Wrote {:s}".format(outfile))
    return


def match_qa(arc_spec, tcent, line_list, IDs, scores, outfile, title=None, path=None):
    """
    Parameters
    ----------
    arc_spec
    tcent
    line_list
    IDs
    scores
    outfile
    title
    path

    Returns
    -------

    """


    # Plot
    pp = PdfPages(outfile)
    plt.figure(figsize=(11, 8.5))
    plt.clf()
    gs = gridspec.GridSpec(1, 1)
    idfont = 'small'

    # Simple spectrum plot
    ax_spec = plt.subplot(gs[0])
    ax_spec.plot(np.arange(len(arc_spec)), arc_spec, 'k')
    ymin, ymax = 0., np.max(arc_spec)
    ysep = ymax*0.03
    mn_yline = 1e9

    # Standard IDs
    clrs = dict(Perf='green', Good='blue', Ok='orange',
                Perfect='green')
    clrs['Very Good'] = 'blue'
    for kk, score in enumerate(scores):
        x = tcent[kk]
        # Color
        try:
            clr = clrs[score]
        except KeyError:
            clr = 'gray'
        yline = np.max(arc_spec[int(x)-2:int(x)+2])
        mn_yline = min(mn_yline, yline)
        # Tick mark
        ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], '-', color=clr)
        if score in ['Perf', 'Good', 'Ok', 'Perfect', 'Very Good']:
            # Label
            imin = np.argmin(np.abs(line_list['wave']-IDs[kk]))
            row = line_list[imin]
            lbl = '{:s} {:.4f}'.format(row['ion'], row['wave'])
            # label
            ax_spec.text(x, yline+ysep*1.3, '{:s}'.format(lbl), ha='center', va='bottom',
                size=idfont, rotation=90., color=clr)
    # Axes
    ax_spec.set_xlim(0., len(arc_spec))
    ax_spec.set_ylim(ymin, ymax*1.3)
    ax_spec.set_xlabel('Pixel')
    ax_spec.minorticks_on()
    ax_spec.set_ylabel('Counts')
    if title is not None:
        ax_spec.text(0.04, 0.93, title, transform=ax_spec.transAxes,
                     size='x-large', ha='left')#, bbox={'facecolor':'white'})
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    pp.savefig(bbox_inches='tight')
    pp.close()
    plt.close()
    return
