import numpy as np
import ecip
from xsec import read_xsec, XSEC
import fitecip

"""
Used in conjunction with 'fitecip.py' to fit multiple cross sections at once. 

Dependencies include:
1) ecip_fxn.for (need to compile with f2py as the executable "ecip.x"). This fortran
    function actually calculates the raw ECIP cross section.
2) xsec.py: for reading in the raw cross section and for setting up the ECIP class.
3) fitecip.py: the base code for generating scaled ECIP fits to raw cross sections.
    A least squares fit is applied to a raw cross section to determine a 
    multiplying constant, 'a' for 'a * unscaled_ecip'.

The format for the cross sections in 'scale_ecip_input' is, per line/cross sectin:

term-specific_ion._pot.   occ_#   xsec_file   output_fit_file   output_graph_file   output_scaled_ecip_xsec_file

"""


fitecip.fit_ecip_winfile('scale_ecip_input')
