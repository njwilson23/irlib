#! /usr/bin/env python
#
#   Combines existing HDF datasets into a single file
#

import irlib
import h5py
import sys, getopt

def print_syntax():
    print """
    SYNTAX: h5_consolidate INFILE1 INFILE2 [...] -o OUTFILE
    """

optlist, fins = getopt.gnu_getopt(sys.argv[1:], 'o:')
optdict = dict(optlist)

if '-o' in optdict.keys():
    fout = optdict['-o']
else:
    print_syntax()
    sys.exit()

h5out = h5py.File(fout, 'w')

# Read each input HDF in the order provided, and copy each line to h5out
i = 0
for fin in fins:
    print "Copying " + fin
    h5in = h5py.File(fin, 'r')
    # Get each line in the input HDF
    groups = [item[0] for item in h5in.items()]
    if 'LabVIEW Boolean' in groups:
        groups.pop(groups.index('LabVIEW Boolean'))
    groups.sort(key=lambda s: int(s.split('_')[1]))
    for group in groups:
        destgroup = 'line_'+str(i)
        print "\t{0} named {1}".format(group, destgroup)
        h5in.copy(h5in[group], h5out, name=destgroup)
        i += 1
    h5in.close()

h5out.close()

