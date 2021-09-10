#! /usr/bin/env python
#
#   irgram - viewing tool for survey lines
#

import sys, getopt, os.path
from irlib.misc import ExtractLine, PlotLine
import traceback

def print_syntax():
    print( """
    STNTAX: plotline infile Nline [options]

    Options:
        -G[file]    save segment to a file
    """)

optlist, args = getopt.gnu_getopt(sys.argv[1:], 'G:', ['clobber'])
optdict = dict(optlist)

# Test the syntax - display error if it doesn't make sense
if len(args) != 2:
    print_syntax()
    sys.exit(1)
else:
    infile = args[0]
    line = int(args[1])

if '-G' in optdict.keys():
    outfile = optdict['-G']
    if os.path.isfile(outfile) and not '--clobber' in optdict.keys():
        sys.stderr.write("Output file already exists\n")
        sys.exit(1)
else:
    outfile = None

# Pull out the requested values
bounds = (None, None)   # TEMPORARY
try:
    line_data = ExtractLine(infile, line, bounds)
except:
    #traceback.print_exc()
    sys.stderr.write("Error extracting line data\n")
    sys.exit(1)

# Determine whether output goes to the display or to a file
if '-G' in optdict.keys():
    outfile = optdict['-G']
else:
    outfile = None

# Plot it
title = "{0} line {1}".format(infile, line)
try:
    PlotLine(line_data, outfile, title=title)
except:
    traceback.print_exc()
    sys.stderr.write("Error while plotting\n")
    sys.exit(1)
