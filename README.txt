/* Andrey Elagin, September 5, 2014
 * LightReco is a light-weight standalone vertex and (in the near future) directionality
 * reconsruction code for 0vbb-decay events in liquid scintillator. The code is based on
 * quadruplet-based vertex-finding method by Michael Smy. Many lines are directly copied
 * from WCSimAnalysis package.
 *
 * This file contains instructions on how to use and notes on significant updates.
 */

- LightReco.C is the main file
- help_func.C contains a few filtering functions for initial reconstruction of
              the off-center events
- run_reco.C simple run interface

--------------------
HOWTO RUN

1) The code is designed to run from a command line to be compatible with batch jobs.
Assuming the ROOT has been installed type:

$ root -b -l -q run_reco.C\(\"f_input.root\",\"f_output\",12\)

2) To run from ROOT interactively type:
$ root
[0] .L LightReco.C+
[1] LightReco([specify input parameters])
=========================

Input files are located at /u/nobackup/lwinslow/elagin/data/sph_out/

The last parameter is explained in run_reco.C

Control total number of events by changing global variable EVT_NUM in
LightReco.C
