""" Script that converts FIDs calculated from Line/Location indices to FIDs
read from Gather.metadata.

Historical explanations:

Prior to August 2013, the FIDs used in picking databases were calculated based
on the line and location number within a particular gather instance,
essentially replicating the functionality in recordlist.py. This proved
problematic, because processed cache data might not have the same number and
order of traces as the originals, and so picking files depended on a particular
filtering regime.

Following commit 30458624b9417fc7b7d8cdd4f1b07d058d344d46, new picking files
are created with FIDs based on the original raw file FID, i.e. the FID retained
by Gather.metadata. This script converts the FIDs in old picking files to the
new FIDs, so that old picking files may be read properly.
"""

import os
import argparse
import irlib


def parse_command_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--hdf5", required=True, help="HDF5 database")
    parser.add_argument("-i", "--indir", required=True, help="input picking files")
    parser.add_argument("-o", "--outdir", required=True, help="output picking files")
    parser.add_argument("-c", "--cachedir", help="radar caches")
    parser.add_argument("--clobber", action="store_true", help="overwrite existing pickfiles in the output directory")
    return parser.parse_args()


def get_pick_fnm(infile, lino):
    """ Autogenerate a filename for pickfiles. """
    fnm = os.path.join('picking',
            '{0}_line{1}.csv'.format(os.path.basename(infile).split('.')[0], lino))
    return fnm


if __name__ == "__main__":

    # infile shouldn't be necessary - get the proper names from the input HDF5

    args = parse_command_args()

    surv = irlib.Survey(args.hdf5)

    for line in surv.GetLines():

        lino = line.split("_")[1]
        pickfile = get_pick_fnm(args.hdf5, lino)

        if args.cachedir:
            use_cache = True
            cachedir = args.cachedir
        else:
            use_cache = False
            cachedir = ""

        lg = surv.ExtractLine(lino, fromcache=use_cache, cache_dir=cachedir)
        fh = irlib.FileHandler(pickfile, lino)

        old_fids = fh.fids
        new_fids = lg.metadata.fids

        if len(old_fids) == len(new_fids):
            fh.fids = new_fids

            pickfilename = os.path.split(fh.fnm)[1]
            fh.fnm = os.path.join(args.outdir, pickfilename)

            if (not os.path.isfile(fh.fnm)) or args.clobber:
                fh.Write()

        else:
            print pickfile
            print "\tError - old pick file had {0} recs but line had {1} recs".format(len(old_fids), len(new_fids))




