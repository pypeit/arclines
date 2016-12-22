""" Module for building arcline lists
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import pdb

from collections import OrderedDict

from astropy.table import Table, vstack

from arclines import io as arcl_io
from arclines import utils as arcl_utils
from arclines import load_source
from arclines import defs

import arclines # For path
llist_path = arclines.__path__[0]+'/data/lists/'
src_path = arclines.__path__[0]+'/data/sources/'


def by_hand(llist_dict, write=False):
    """
    Returns
    -------

    """
    # By-hand
    print("=============================================================")
    print("Adding lines by hand!")
    print("=============================================================")
    handIDs = arcl_io.load_by_hand()
    uhions = np.unique(handIDs['ion'].data)
    # Loop on ID ions
    for ion in uhions:
        # Parse
        idx = handIDs['ion'] == ion
        sub_tbl = handIDs[idx]
        # Update only
        llist_dict[ion], updated = update_line_list(llist_dict[ion], sub_tbl, None, None)
        if write and updated:
            ion_file = llist_path+'{:s}_lines.dat'.format(ion)
            arcl_io.write_line_list(llist_dict[ion], ion_file)


def init_line_list():
    """ Initialize a Table for a linelist
    Rigidly enforces table column formats
    Strings are the most annoying

    Returns
    -------
    init_tbl : Table
      One dummy row
    """
    # Get str lengths from defs
    len_line = defs.str_len()['ion']
    len_src = defs.str_len()['Source']
    # Load sources to check
    sources = arcl_io.load_source_table()
    src_files = sources['File'].data
    if src_files.dtype.itemsize > len_src:
        raise ValueError("Sources now exceeds table.  Should fix source name")
    dummy_src = str('#')*len_src
    # Arc Line name
    dummy_line = str('#')*len_line
    #

    # Dict for Table
    idict = OrderedDict()
    idict['ion'] = dummy_line
    idict['wave'] = 0.
    idict['NIST'] = 0
    idict['Instr'] = 0  # Flag for instrument
    idict['amplitude'] = 0
    idict['Source'] = dummy_src

    # Table
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    init_tbl = Table(lst, names=tkeys)

    # Return
    return init_tbl


def add_instr_source(line_tbl, instr, source_file):
    """
    Parameters
    ----------
    tbl
    instr
    source_file

    Returns
    -------
    Modified in place

    """
    str_len_dict = defs.str_len()
    # Add instrument flag
    line_tbl['Instr'] = defs.instruments()[instr]
    # Add source
    source = np.array([source_file]*len(line_tbl),
                      dtype='S{:d}'.format(str_len_dict['Source']))
    line_tbl['Source'] = source


def create_line_list(line_tbl, source_file, instr,
                     unknown=False, ions=None):
    """
    Parameters
    ----------
    line_tbl
    source_file : str
    instr : str
      Converted to a flag

    """
    # Init -- Mainly to insure formatting
    ini_tbl = init_line_list()
    # Add columns (in place)
    add_instr_source(line_tbl, instr, source_file)
    # Stack
    new_tbl = vstack([ini_tbl, line_tbl], join_type='exact')
    # Cut off first line
    cut_tbl = new_tbl[1:]

    # Unknown list??
    if unknown:
        if ions is None:
            raise IOError("Must provide list of possible ions if unknown")
        line_flag = get_line_flag(ions)
        cut_tbl['line_flag'] = line_flag

    # Return
    return cut_tbl


def get_line_flag(ions):
    """
    Parameters
    ----------
    ions

    Returns
    -------

    """
    line_dict = defs.lines()
    line_flag = 0
    for uion in ions:
        line_flag += line_dict[uion]
    # Return
    return line_flag


def source_to_line_lists(source, write=False, llist_dict=None):
    """
    Parameters
    ----------
    source
    write
    scratch : bool, optional
      Build from scratch -- EXPERTS ONLY

    Returns
    -------

    """
    if llist_dict is None:  # For building without writing
        llist_dict = {}
    # Parse
    src_dict = load_source.load(source)
    # Check
    if src_dict['ID_lines'] is None:
        print("No IDs in source: {:s}".format(source['File']))
        return llist_dict

    # Unique ions
    uions = arcl_utils.unique_ions(source, src_dict=src_dict)

    # Loop on ID ions
    for ion in uions:
        ion_file = llist_path+'{:s}_lines.dat'.format(ion)
        # Parse
        idx = src_dict['ID_lines']['ion'] == ion
        sub_tbl = src_dict['ID_lines'][idx]
        # Generate?
        if (not os.path.isfile(ion_file)) and (ion not in llist_dict.keys()):
            # New
            llist_dict[ion] = create_line_list(sub_tbl, source['File'], source['Instr'])
            if not write:
                print("Would generate line list:\n   {:s}".format(ion_file))
            else:
                print("Generating line list:\n   {:s}".format(ion_file))
                arcl_io.write_line_list(llist_dict[ion], ion_file)
        else:
            try:
                llist_dict[ion], updated = update_line_list(llist_dict[ion], sub_tbl, source['File'], source['Instr'])
            except KeyError:
                raise KeyError("You are trying to build from scratch but didn't remove {:s}".format(ion_file))
            # Write
            if write and updated:
                arcl_io.write_line_list(llist_dict[ion], ion_file)
    # Return
    return llist_dict


def update_line_list(line_list, new_lines, source_file, instr, tol_wave=0.1, NIST_tol=0.0001):
    """ Update/add to lines in line list as applicable
    Not for use on UNKNWN lines

    Parameters
    ----------
    line_list : Table
    new_lines : Table
    source_file : str
    instr : str
      Converted to a flag
    tol_wave : float, optional
      Matching tolerance in wavelength
      Anything closer than this, even if real, is trouble

    """
    start = "\x1B["
    end = "\x1B[" + "0m"

    # Add columns (in place)
    if 'Instr' not in new_lines.keys():
        add_instr_source(new_lines, instr, source_file)

    # Loop to my loop
    updated = False
    for line in new_lines:
        # NIST
        if line['NIST'] != 1:
            print("Not ready for this")
            pdb.set_trace()
        # Search for wavelength match within tolerance
        mtch_wave = np.where(np.abs(line_list['wave']-line['wave']) < tol_wave)[0]
        if len(mtch_wave) == 0:
            print(start+"1:34m"+"ADDING "+end+"the following line to {:s} line list".format(line['ion']))
            print(line)
            line_list = vstack([line_list, line]) # Insures columns are matched
            updated = True
        elif len(mtch_wave) == 1:
            idx = mtch_wave[0]
            if np.abs(line_list['wave'][idx]-line['wave']) > NIST_tol:
                print("Bad match for a NIST line")
                print(line)
                print(line_list[idx])
                pdb.set_trace()
            else:  # Check instrument
                if (line_list['Instr'][idx] % (2*line['Instr'])) >= line['Instr']:
                    pass
                else:
                    line_list['Instr'][idx] += line['Instr']
                    print("Updating INSTRUMENT in this line:")
                    print(line_list[idx])
                    updated = True
    # Sort
    line_list.sort('wave')
    # Return
    return line_list, updated


def update_uline_list(new_lines, source_file, instr, line_file,
                      ions, tol_wave=0.5, write=False):
    """ Update/add to UNKNWN line list as applicable

    Parameters
    ----------
    line_tbl
    source_file : str
    instr : str
      Converted to a flag
    line_file : str
    ions : list

    tol_wave : float, optional
      Matching tolerance in wavelength
      Anything closer than this, even if real, is trouble

    """
    # Load
    line_list = arcl_io.load_line_list(line_file)
    # Add columns (in place)
    add_instr_source(new_lines, instr, source_file)

    line_flag = get_line_flag(ions)
    new_lines['line_flag'] = line_flag

    # Loop to my loop
    updated = False
    for line in new_lines:
        # Search for wavelength match within tolerance
        mtch_wave = np.where(np.abs(line_list['wave']-line['wave']) < tol_wave)[0]
        if len(mtch_wave) == 0:
            if write is False:
                print("Would add the following line to {:s}".format(line_file))
                print(line)
            line_list = vstack([line_list, line]) # Insures columns are matched
            updated = True
        elif len(mtch_wave) == 1:
            idx = mtch_wave[0]
            # Check instrument
            if (line_list['Instr'][idx] % (2*line['Instr'])) >= line['Instr']:
                pass
            else:
                line_list['Instr'][idx] += line['Instr']
                if write is False:
                    print("Would update instrument in this line:")
                    print(line_list[idx])
                updated = True
    # Sort
    line_list.sort('wave')
    # Write
    if write and updated:
        arcl_io.write_line_list(line_list, line_file)


def master_build(write=False, nsources=None, plots=True, verbose=True):
    """ Master loop to build the line lists

    Parameters
    ----------
    check_only : bool, optional
      Only check to see what would be created, i.e. do not execute
    nsources : int, optional
      Mainly for step-by-step build
    plots : bool, optional
      Generate plots, if write=True too
    """
    # Load sources
    sources = arcl_io.load_source_table()
    if nsources is not None:
        sources = sources[0:nsources]

    # Loop on sources
    for source in sources:
        print("=============================================================")
        print("Working on source {:s}".format(source['File']))
        print("=============================================================")
        # Load line table
        ID_lines, U_lines = load_source.load(source['File'], source['Format'],
                                             source['Lines'].split(','),
                                             plot=plots, wvmnx=[source['wvmin'], source['wvmax']])
        # Lines (Double check)
        src_lines = source['Lines'].split(',')
        if ID_lines is None:
            uions = src_lines
        else:
            # Check
            uions = np.unique(ID_lines['ion'].data)
            for src_line in src_lines:
                if src_line not in uions.tolist():
                    raise ValueError("Line {:s} not found in ID_lines".format(src_line))
            # Loop on ID ions
            for ion in uions:
                # Parse
                idx = ID_lines['ion'] == ion
                sub_tbl = ID_lines[idx]
                # Generate?
                ion_file = llist_path+'{:s}_lines.dat'.format(ion)
                if not os.path.isfile(ion_file):
                    if not write:
                        print("Would generate line list:\n   {:s}".format(ion_file))
                    else:
                        print("Generating line list:\n   {:s}".format(ion_file))
                        create_line_list(sub_tbl, source['File'],
                                     source['Instr'], ion_file)
                else:
                    update_line_list(sub_tbl, source['File'],
                                     source['Instr'], ion_file, write=write)
        # UNKNWN lines
        if U_lines is None:
            continue
        unk_file = llist_path+'UNKNWN_lines.dat'
        # Check against 'complete' NIST and our line lists
        mask, _ = arcl_utils.vette_unkwn_against_lists(U_lines, uions, verbose=verbose)
        if np.sum(mask) == 0:
            continue
        if not os.path.isfile(unk_file): # Generate?
            if write:
                print("Generating line list:\n   {:s}".format(unk_file))
                create_line_list(U_lines[mask>0], source['File'], source['Instr'],
                             unk_file, unknown=True, ions=uions)
        else: # Update
            update_uline_list(U_lines[mask>0], source['File'], source['Instr'],
                              unk_file, uions, write=write)


