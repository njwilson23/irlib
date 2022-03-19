#! /usr/bin/env python
'''
h5_export.py - export a line from HDF5 to an ASCII, REFLEX or binary file

Formerly known as h52a.py
Nat Wilson
Derek Mueller/Adriana Caswell


TODO - combine with h52mat.py eventually (it can all happen here including the cleaning steps)
TODO - also add export back to h5 format.  Could be useful. 
    

'''

import sys
import os.path
import argparse
import numpy as np

import irlib

def clobbercheck(f,clobber):
    if (not clobber) and os.path.exists(f):
        print("File exists.  To overwrite use --clobber")
        sys.exit(1)

def line2ascii(file, L):
    """ Write a numpy array arr as an ascii grid to file f. """
    with open(file, 'w') as f:
        for row in L.data:
            f.write(str(row).strip('[').strip(']').replace('\n','') + '\n')


        
prog_description = """h52a.py - export a line from HDF5 to an ASCII, REFLEX or BINARY file
"""
parser = argparse.ArgumentParser(description = prog_description)
parser.add_argument('outformat', choices=['ascii', 'binary', 'reflex'], help="Select which format to export to - either ascii, binary or reflex")
parser.add_argument('infile', help="input HDF (.h5) filename, with or without path")
parser.add_argument('-o','--outfile', help="output filename, basename only NO extension; defaults to infile")
parser.add_argument('-l','--line', help='line number to export - defaults to all')
parser.add_argument('--clobber', help = 'overwrite existing files', default=False, action='store_true')

args = parser.parse_args()

# check input file 
infile = args.infile
if not os.path.exists(infile):
    print("No such file: {}".format(infile))
    sys.exit(1)

p = os.path.dirname(infile)
filename = os.path.basename(infile)
basename = os.path.splitext(filename)[0]

if args.outfile == None:
   outfile = os.path.join(p,basename)
else: 
    outfile = os.path.join(p,args.outfile)

    
# Open INFILE as an HDF5 dataset
S = irlib.Survey(infile)
lns = [int(line.split('_')[1]) for line in S.GetLines()]


if args.line == None:  # default to all lines
    lines = lns
else:
    lines = [int(l) for l in args.line] # convert to list of int

for l in lines:
    try:
        L = S.ExtractLine(l)
    except:
        sys.stderr.write("h52a: error extracting line {} data \n".format(l))
        continue
    if args.outformat == "reflex":
        fout = '{}_{:02d}.ipr'.format(outfile,l)
        clobbercheck(fout, args.clobber)
        np.savetxt(fout, np.transpose(L.data), fmt='%.4e', delimiter=' ') 
        
    elif args.outformat == 'ascii':
        fout = '{}_{:02d}.txt'.format(outfile,l)
        clobbercheck(fout, args.clobber)
        line2ascii(fout, L)
    else: 
        fout = '{}_{:02d}.bin'.format(outfile,l)
        clobbercheck(fout, args.clobber)
        L.data.tofile(fout)
        
    print('...saved file {}'.format(fout))