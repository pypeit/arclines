""" Module for arc line plots
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import pdb

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

# Default Path
import arclines
plot_path = arclines.__path__[0]+'/data/plots/'


def arc_ids(arc_spec, xIDs, IDs, outfile, title=None, path=None,
            clobber=False, extras=None):
    """
    QA for Arc spectrum

    Parameters
    ----------
    arc_spec : ndarray
      Arc spectrum
    xIDs : ndarray or list
      Pixel values of ID'd lines
    IDs : ndarray or list
      str array of ID labels
    extras : dict, optional
      x: list of x values
      IDs: list of str labels
    outfile : str
      Name of output file
    """
    # Path
    if path is None:
        path = plot_path

    # Begin
    if os.path.isfile(plot_path+outfile):
        print("Plot {:s} exists.  Remove if you wish to remake it".format(outfile))

    pp = PdfPages(plot_path+outfile)
    plt.figure(figsize=(11, 8.5))
    plt.clf()
    gs = gridspec.GridSpec(2, 1)
    idfont = 'small'

    # Simple spectrum plot
    for qq in range(2):
        ax_spec = plt.subplot(gs[qq])
        ax_spec.plot(np.arange(len(arc_spec)), arc_spec)
        ymin, ymax = 0., np.max(arc_spec)
        ysep = ymax*0.03
        mn_yline = 1e9
        # Standard IDs
        for kk, x in enumerate(xIDs):
            yline = np.max(arc_spec[int(x)-2:int(x)+2])
            mn_yline = min(mn_yline, yline)
            # Tick mark
            ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], 'g-')
            # label
            ax_spec.text(x, yline+ysep*1.3, '{:s}'.format(IDs[kk]), ha='center', va='bottom',
                size=idfont, rotation=90., color='green')
        # Extras?
        if extras is not None:
            for kk, x in enumerate(extras['x']):
                yline = np.max(arc_spec[int(x)-2:int(x)+2])
                mn_yline = min(mn_yline, yline)
                # Tick mark
                ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], 'r-')
                # label
                ax_spec.text(x, yline+ysep*1.3, '{:s}'.format(extras['IDs'][kk]), ha='center', va='bottom',
                    size=idfont, rotation=90., color='red')
        # Axes
        ax_spec.set_xlim(0., len(arc_spec))
        if qq==1:
            ax_spec.set_yscale("log", nonposy='clip')
            ax_spec.set_ylim(mn_yline/2., 5*ymax)
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
    return

