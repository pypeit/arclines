.. highlight:: rest

**************
Match to Lines
**************

This document describes a script that will
match an input spectrum to the arclines database.

Basic Methodology
=================

The script will scan through 3 regions of parameter space:

 1. linelist -- it will try using the unknown line list and without;
 2. pix_tol -- values of 1 and 2 are attempted
 3. wv_cen  -- solutions are scanned from +/- 1000. from the input value

Script
======

The script is named arclines_match and has the following usage::

    usage: arclines_match [-h] [--outroot] [--min_ampl MIN_AMPL]
                      spectrum wvcen disp lines

    Match input spectrum to arclines line lists

    positional arguments:
      spectrum             Spectrum file (.ascii, .fits)
      wvcen                Guess at central wavelength (within 1000A)
      disp                 Accurate dispersion (Ang/pix)
      lines                Comma separated list of lamps

    optional arguments:
      -h, --help           show this help message and exit
      --outroot            Root filename for plot, IDs
      --min_ampl MIN_AMPL  Minimum amplitude for line analysis [default: 100]

Here is a typical call::

    arclines_match tests/files/LRISb_600_spec.ascii 4000. 1.26 CdI,HgI,ZnI

where I've used an example spectrum in the test suite.

Spectrum Formats
================

The following arc spectrum formats are allowed:

======= =====================================================================
Type    Description
======= =====================================================================
ascii   ASCII file which must be read by Table.read(filename,format='ascii').
        Program will assume that the first column contains the arc spectrum
fits    FITS file where the first extension contains an ndarray of the spectrum
hdf5    Currently only reads the LowRedux conversion file (see misc.low_redux)
======= =====================================================================
