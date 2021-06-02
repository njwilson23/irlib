#! /usr/bin/env python
#   Generate pickled line caches for a survey, so that viewing lines later
#   is much faster.
#

from irlib import Survey
import os, sys, argparse
import traceback

#replacing getopt and def syntax() with argparse
prog_description = """
    SYNTAX: h5_generate_caches HDF_SURVEY [OPTIONS]
         
   **REMEMBER TO WRITE DOWN YOUR SETTINGS, SO YOU CAN RECREATE THESE FILES**
    """
prog_epilog = " Example: h5_generate_caches.py survey.h5  "

parser = argparse.ArgumentParser(description = prog_description, epilog=prog_epilog)
parser.add_argument("infile", help="input HDF (.h5) filename")
parser.add_argument("-d", "--dir_cache", help="cache directory (default: cache/)", default='cache/')
parser.add_argument("-r", "--remove_within", help="remove stationary traces by averaging all traces within # m (defaults to 0 m or off), recommend 3 for L1 GPS", 
                    default=0,type=float)
parser.add_argument("--dc", help="specify datacapture (default: 0)", default=0,type=int)
parser.add_argument("-n","--remove_nans", help="remove traces with NaN coordinates", action="store_true")
parser.add_argument("-i","--interp_nans", help="interpolate over NaN coordinates (overrides --remove_nans)", action="store_true")
parser.add_argument("-s", "--smoothen_coords", help="smoothen coordinates (overrides --interp_nans)", action="store_true")
parser.add_argument("-b", "--remove_blanks", help="remove blank traces caused by triggering failure", action="store_true")
parser.add_argument("-g", "--fix_gps", help="fix static GPS issues", action="store_true")
parser.add_argument("-f", "--force_cache", help="force regeneration of existing caches", action="store_true")
parser.add_argument("-q", "--quiet", help="silence standard output", action="store_true")
parser.add_argument("-v", "--verbose", help="print failed datacaptures", action="store_true") 

args = parser.parse_args()

if not os.path.isdir(args.dir_cache):
    os.makedirs(args.dir_cache)

S = Survey(args.infile)

if not args.quiet:
    print("Working on {0}".format(args.infile))

lines = S.GetLines()

for line in lines:
    line_no = line.split('_')[1]
    cache_fnm = S.GetLineCacheName(line_no, dc=args.dc, cache_dir=args.dir_cache)
    if os.path.isfile(cache_fnm) and not args.force_cache:
        pass
    else:
        if not args.quiet:
            print("\tCaching line {0}, datacapture {1}...".format(str(line),
                                                                  str(args.dc)))
        try:
            L = S.ExtractLine(line_no, datacapture=args.dc, verbose=args.verbose)
        except (AttributeError, IndexError):
            print("\t\tData invalid")
            L = None
        except KeyboardInterrupt:
            sys.exit(0)

        if L is not None:
            try:
                L.RemoveBadLocations()
                if args.fix_gps:
                    L.FixStaticGPS()
                try:
                    if args.smoothen_coords:
                        L.SmoothenGPS()
                    elif args.interp_nans:
                        L.InterpolateGPSNaNs()
                    elif args.remove_nans:
                        L.RemoveGPSNaNs()
                except ValueError:
                    print("\t\tNo valid UTM coordinates")
                if args.remove_blanks:
                    L.RemoveBlankTraces()
                if args.remove_within > 0.0:
                    L.RemoveStationary(args.remove_within)
                L.Dump(cache_fnm)
                del L
            except:
                traceback.print_exc()

