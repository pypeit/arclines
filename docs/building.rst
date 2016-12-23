.. highlight:: rest

*******************
Building Line Lists
*******************

This document describes how to run scripts to build the
line lists.

**WARNING:**  This is for Expert's only.

From Scratch
============

This is truly for experts.  This generally occurs in 3
stages.

Standard Line Lists
-------------------

Here is the recommended sequence::

    rm previous_existing_line_lists
    ./scripts/build_from_scratch.py  # Scan for complete rebuild
    ./scripts/build_from_scratch.py --write   # Write standard line list files

Unknowns
--------

Here is the standard sequence::

    ./scripts/build_from_scratch.py --unknowns  # Check the output
    ./scripts/build_from_scratch.py --unknowns --write   # Write standard line list files

Plots
-----

Generate PDF summaries of each source::

    ./scripts/build_from_scratch.py --plots

