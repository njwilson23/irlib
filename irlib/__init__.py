"""
   ICE-PENETRATING RADAR LIBRARY

   Provides the underlying classes for ice radar utilities.
   Depends on:
       h5py                HDF file interface
       numpy/scipy         Numerics and math libraries
       matplotlib          Plotting

   For some specific functionality, depends on:
       mig_fk              FK migration module
       agc_cy              Cython-accelerated AGC
       pywavelet           Wavelet transform function
"""

import survey
from survey import *
import gather
from gather import *
import recordlist
from recordlist import *
import filehandler
from filehandler import *
import irlib
from irlib import *

from filter_defs import ApplyFilter

try:
	import itools
except:
	print "Warning: itools not imported"
