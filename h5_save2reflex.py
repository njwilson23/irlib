#! /usr/bin/env python
#
#
"""
Created on Mon Jul 25 17:44:11 2016

@author: dmueller
"""

import os
import irlib
from irlib.itools import *
from pylab import *
import numpy as np


args = sys.argv[1:]
if len(args) != 1:
    print("""
    SYNTAX: h5_save2reflex INFILE

        Saves each line in an h5 file as an ASCII matrix with each row as
        an individual trace so it can be imported into REFLEX 
    """)
    sys.exit(1)

INFILE = sys.argv[1]
    
# Open INFILE as an HDF5 dataset
if not os.path.exists(INFILE):
    print("No such file: {0}".format(INFILE))
    sys.exit(1)

p = os.path.dirname(INFILE)
filename = os.path.basename(INFILE)
basename, ext = os.path.splitext(filename)

S = irlib.Survey(os.path.join(p,filename))
#S.GetChannelsInLine(1)

print()
print("Found {0} lines to convert".format(len(S.GetLines())))

for line in S.GetLines():
    lnum = int(line.split('_')[1])
    try:
        L = S.ExtractLine(lnum)
        np.savetxt(os.path.join(p,basename+"_"+str(lnum)+'.ipr'), transpose(L.data), fmt='%.4e', delimiter=' ') 
        print('saved line {0}'.format(lnum))
    except:
        print('problem with line {0}'.format(lnum))
        
