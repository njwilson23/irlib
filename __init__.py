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

import irlib
from irlib.survey import Survey
from irlib.gather import Gather, CommonOffsetGather, CommonMidpointGather
from irlib.gather import LineGather
from irlib.recordlist import RecordList
from irlib.filehandler import FileHandler

from irlib.filter_defs import ApplyFilter

try:
	import irlib.itools
	import irlib.itools as itools
except:
	print "Warning: itools not imported"
