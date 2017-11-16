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

Here is a summary of the color codes:

======== =====================================================================
Color    Description
======== =====================================================================
Green    New line (to be) added to ID database
Blue     Line previously exists in ID database
Gray     UNKNOWN line previously exists in UNKNOWN database
======== =====================================================================

Add Source
==========

The following steps will add a new arcline source to the database.
Do them in the order below.

First and foremost:  Choose wisely.  Stick with sources that
are well vetted.

Code Edits
----------

If you are adding a new instrument, add to the list in defs.instruments

If you will be generating a new line list (i.e. a new lamp), add to
defs.lines


Source List
-----------

Add the new source entry to the bottom of the arcline_sources.ascii file.

Run Add Setup
-------------

Inspection
++++++++++

Run the add setup script without --write::

    ./scripts/add_source.py

Inspect the output to the terminal.

  - Check the lines that would be added to the main line lists
  - Check the UNKNOWN lines that would be purged (should only match to new ones)
  - Check the UNKNOWN lines that would be added

Inspect the output plot.

  - Green are new lines to be added
  - Blue were present previously
  - Orange are new UNKNOWN lines
  - Grey are old UNKNOWN lines

Write
+++++

Write to the disk (and github)::

    ./scripts/add_source --write

Decide whether to include UNKNOWN lines.  This is only
recommended when there is a large set or for a very sparse
arc spectrum (i.e. where every line matters).::

    ./scripts/add_source --write --no_unknowns

Remake plots?
-------------

Consider remaking the plots for the entire database.
