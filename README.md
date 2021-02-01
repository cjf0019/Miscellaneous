Optizing atomic structure for heavy atoms proves tough, particularly for neutral species. The directory "/AtomicPhysics/ParallelizedStructureOpt" provides code for optimizing atomic structure, through a parallelized implementation of the program AUTOSTRUCTURE (cit. https://ui.adsabs.harvard.edu/abs/2016ascl.soft12014B/abstract, written by Nigel Badnell). Structure optimization of lambda parameters, max one per orbital, can be done either with L-BFGS-B or from a brute force calculation across a parameter grid. The cost function contour, for both energy and A-value optimization separately, can be further graphed (graph.py) for a more comprehensive view of the optimization process. The user can specify the number of workers as the number of concurrent AUTOSTRUCTURE runs at one time and can be run on HPC clusters. 


In addition, the following codes represent the most useful developed Python codes from the dissertation work "The role of excited state electron-impact ionization and configuration mixing in the
population modeling of neutral Ne and W in fusion-capable plasmas using large-scale R-matrix collisional data". Broadly, they represent a suite of programs for handling large-scale R-matrix and generalized collisional radiative modeling calculations. They can be found in the "/AtomicPhysics" subdirectory:

1. pystgsig.py: Extract cross sections from an OMEGA file, given initial and final transition
numbers. Can optionally return cross sections as collision strengths.

2. pydstgsig.py: Same as ‘pystgsig.py’ but takes a dstgsig file (same format as original
stgsig code) and can be used to extract specific transitions through a list of 0’s and 1’s
after the namelist: one value on each succeeding line, with 0 and 1 designating not to
include and to include respectively. Can use “generatestgsiglist.py” to generate this list.

3. omegautility.py: Includes OMEGA and XSEC classes, the former a wrapper for
OMEGA files, storing their basic information, and the latter for cross sections. The
OMEGA class can be used to extract cross sections (see pystgsig.py and pydstgsig.py).
Rost-Pattard and Younger fits as well as both fit and raw plots can be performed using
the XSEC class (see ionsnfitandplot.py for such usage). Cross sections reflecting an n4ˆ
scaling (for ionization of high n states) can also be generated.

4. fitecip winfile.py: Generate scaled ECIP fits to raw cross sections from a file
‘scale ecip input’ that includes sets of configuration, term, or level-specific ionization
potentials, the occupation number of the ionized shell, the cross section file name, and
the output file names for the fit, graph, and raw ECIP cross section.

5. fitecip.py: the base code for generating scaled ECIP fits to raw cross sections. A least
squares fit is applied to a raw cross section to determine a multiplying constant, ’a’ for ’a
unscaled ecip.’

6. ionsnfitandplot.py: Generate Rost-Pattard and Younger fits to raw cross sections in
lower and higher incident energies respectively.

7. generatestgsiglist.py: See the “pydstgsig.py” description. Reads in a TERMS file from
AUTOSTRUCTURE.

8. xsec.py: Contains the base class XSEC for reading, writing and fitting cross sections.

9. adf04utility.py: Contains various functions related to processing adf04 files (which store
the rate coefficients of various atomic processes like electron-impact and recombination.)

10. gettransnums.py: Get the transition numbers for OMEGA extraction of ionization cross
sections. NOTE: Assumes elastic collisions are not included. Needs modification if they
are present (one extra transition per level).

11. convertatomicunits.py: Contains a list of functions for converting between atomic units
(eV, Ryd, cm-1).
