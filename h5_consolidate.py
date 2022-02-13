#! /usr/bin/env python
#
#   Combines existing HDF datasets into a single file
#

import h5py
import argparse

prog_description = """
    SYNTAX: h5_consolidate INFILE1 INFILE2 [...] -o OUTFILE

    Combines multiple datasets (>1) into a single concatenated dataset.
    """

parser = argparse.ArgumentParser(description = prog_description)
parser.add_argument('infile', action = 'append', nargs = '*')
parser.add_argument('-o', '--outfile')

args = parser.parse_args()


fout = args.outfile
fins = args.infile[0]

h5out = h5py.File(fout, 'w')


# Read each input HDF in the order provided, and copy each line to h5out
i = 0
for fin in fins:
    print("Copying " + fin)
    h5in = h5py.File(fin, 'r')
    # Get each line in the input HDF
    groups = [item[0] for item in h5in.items()]
    if 'LabVIEW Boolean' in groups:
        groups.pop(groups.index('LabVIEW Boolean'))
    groups.sort(key=lambda s: int(s.split('_')[1]))
    for group in groups:
        destgroup = 'line_'+str(i)
        print("\t{0} named {1}".format(group, destgroup))
        h5in.copy(h5in[group], h5out, name=destgroup)
        i += 1
    h5in.close()

h5out.close()
