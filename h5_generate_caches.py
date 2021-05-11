#! /usr/bin/python
#
#   Generate pickled line caches for a survey, so that viewing lines later
#   is much faster.
#


from irlib import Survey
import os, sys, getopt
import traceback

def print_syntax():
    print("""
    SYNTAX: h5_generate_caches HDF_SURVEY [OPTIONS]

        -d [DIR]            cache directory (default: cache/)
        -g                  fix static GPS issues
        --remove-nans       remove traces with NaN coordinates
        --interp-nans       interpolate over NaN coordinates (overrides --remove-nans)
        -s
        --smoothen-coords   smoothen coordinates (overrides --interp-nans)
        -b                  remove blank traces caused by triggering failure
        -r                  remove stationary traces
        -f                  force regeneration of existing caches
        -q                  silence standard output
        -e                  print failed datacaptures
        --dc=[#]            specify datacapture (default: 0)
    """)

optlist, fins = getopt.gnu_getopt(sys.argv[1:], 'd:gsbrfqe',
                    ['remove-nans', 'interp-nans', 'smoothen-cords', 'dc='])
optdict = dict(optlist)


# Parse input switches
cache_dir = optdict.get('-d', 'cache')
fix_gps = '-g' in optdict
remove_stationary = '-r' in optdict
remove_blanks = '-b' in optdict
force_cache = '-f' in optdict
be_quiet = '-q' in optdict
verbose = '-e' in optdict

smoothen_gps = False
interpolate_nans = False
remove_nans = False
if '-s' in optdict or '--smoothen-coords' in optdict:
    smoothen_gps = True
elif '--interpolate_nans' in optdict:
    interpolate_nans = True
elif '--remove-nans' in optdict:
    remove_nans = True

try:
    dc = int(optdict.get('--dc', 0))
except ValueError:
    print("key for --dc must be an integer")
    sys.exit(1)

try:
    survey_fnm = fins[0]
except IndexError:
    print_syntax()
    sys.exit(1)

if not os.path.isdir(cache_dir):
    os.makedirs(cache_dir)

S = Survey(survey_fnm)

if not be_quiet:
    print("Working on {0}".format(survey_fnm))

lines = S.GetLines()

for line in lines:
    line_no = line.split('_')[1]
    cache_fnm = S.GetLineCacheName(line_no, dc=dc, cache_dir=cache_dir)
    if os.path.isfile(cache_fnm) and not force_cache:
        pass
    else:
        if not be_quiet:
            print("\tCaching line {0}, datacapture {1}...".format(str(line),
                                                                  str(dc)))
        try:
            L = S.ExtractLine(line_no, datacapture=dc, verbose=verbose)
        except (AttributeError, IndexError):
            print("\t\tData invalid")
            L = None
        except KeyboardInterrupt:
            sys.exit(0)

        if L is not None:
            try:
                L.RemoveBadLocations()
                if fix_gps:
                    L.FixStaticGPS()
                try:
                    if smoothen_gps:
                        L.SmoothenGPS()
                    elif interpolate_nans:
                        L.InterpolateGPSNaNs()
                    elif remove_nans:
                        L.RemoveGPSNaNs()
                except ValueError:
                    print("\t\tNo valid UTM coordinates")
                if remove_blanks:
                    L.RemoveBlankTraces()
                if remove_stationary:
                    L.RemoveStationary(3.0)
                L.Dump(cache_fnm)
                del L
            except:
                traceback.print_exc()

