#! /usr/bin/env python
#
#   h5dumpmeta - HDF5 metadata dump
#
#   Utility to write the metadata from raw radar data files to ASCII.
#   Data is sent to stdout, and may be redirected or piped at the
#   terminal.
#

import StringIO
import sys, getopt, os.path
import traceback

def syntax():
    sys.stderr.write("""
    SYNTAX: h5_dumpmeta infile [-f] [--clobber] > outfile

    Options:
        -f          automatically name the output file
        --clobber   overwrite existing files
    \n""")

optlist, args = getopt.gnu_getopt(sys.argv[1:], 'f', ['clobber'])

if len(args) == 0:
    sys.stderr.write("No input given\n")
    syntax()
else:

    import irlib.misc       # Delay import for speed

    for infile in args:

        sys.stderr.write(infile + '\n')

        stringbuffer = StringIO.StringIO()

        try:
            records = irlib.misc.ExtractAttrs(infile,
                                              fout=stringbuffer, 
                                              eastern_hemisphere=False)
        except:
            traceback.print_exc()
            sys.stderr.write("Error reading radar data\n")
            sys.exit(1)

        stringbuffer.seek(0)

        if '-f' in dict(optlist).keys():
            # Print it to an auto-generated file
            outfile = infile.rsplit('.',1)[0] + '.csv'
            if not os.path.isfile(outfile):
                with open(outfile, 'w') as f:
                    f.write(stringbuffer.read())
            else:
                if '--clobber' in dict(optlist).keys():
                    with open(outfile, 'w') as f:
                        f.write(stringbuffer.read())
                    sys.stderr.write(
                        "\t{fnm} overwritten\n".format(
                            fnm=os.path.basename(outfile)))
                else:
                    sys.stderr.write(
                        "\t{fnm} already exists\n".format(
                            fnm=os.path.basename(outfile)))
        else:
            # Print to stdout
            sys.stdout.write(stringbuffer.read())
