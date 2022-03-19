#! /usr/bin/env python
#
#   irtrace - viewing tool for ice radar traces
#
#   Draws a radar trace, along with neighbouring traces. Can either
#   show the trace immediately or save it to a file as a PNG image.
#
#   Takes a survey file, a line number, and a location number as
#   arguments.
#

import sys, getopt, os.path
from irlib.misc import ExtractTrace, PlotTrace
import traceback

def print_syntax():
    print ("""
    SYNTAX: plottrace infile Nline Nlocation [options]

    Options:
        -s          show single trace, without nearby traces
        -f          apply linear gain control filtering
        -r          specify sampling rate (default = 2.5e-9)
        -d          specify vertical spacing (default = 0.1)
        -G[file]    save trace to a file
        --clobber   overwrite existing files
    """)


optlist, args = getopt.gnu_getopt(sys.argv[1:], 'fsr:d:G:', ['clobber'])
optdict = dict(optlist)

if len(args) != 3:
    print_syntax()
    sys.exit(1)
else:
    infile = args[0]
    nline = int(args[1])
    nlocation = int(args[2])

if '-G' in optdict.keys():
    outfile = optdict['-G']
    if os.path.isfile(outfile) and not '--clobber' in optlist:
        sys.stderr.write("Output file already exists\n")
        sys.exit(1)
else:
    outfile = None

if '-r' in optdict.keys():
    rate = optdict['-r']
else:
    rate = 4e-9

if '-d' in optdict.keys():
    spacing = float(optdict['-d'])
else:
    spacing = 0.1

# NOT YET IMPLEMENTED
# To make this work nicely, it would be good to update the structure of the program
# to the OO side of irlib
if '-f' in optdict.keys():
    do_gc = True
else:
    do_gc = False

try:
    D = ExtractTrace(infile, nline, nlocation)
except:
    traceback.print_exc()
    sys.exit(1)

Dn = None
Dnn = None
Dp = None
Dpp = None

if '-s' not in optdict.keys():
    try:
        Dp = ExtractTrace(infile, nline, nlocation-1)
        try:
            Dpp = ExtractTrace(infile, nline, nlocation-2)
        except KeyError:
            pass
    except KeyError:
        pass

    try:
        Dn = ExtractTrace(infile, nline, nlocation+1)
        try:
            Dnn = ExtractTrace(infile, nline, nlocation+2)
        except KeyError:
            pass
    except KeyError:
        pass

PlotTrace(D, Dp=Dp, Dn=Dn, Dpp=Dpp, Dnn=Dnn, outfile=outfile, rate=rate, spacing=spacing)

