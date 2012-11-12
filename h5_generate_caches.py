#! /usr/bin/python
#
#   Generate pickled line caches for a survey, so that viewing lines later
#   is much faster.
#


import irlib
import os, sys, getopt

def print_syntax():
    print """
    SYNTAX: gen_cache HDF_SURVEY [OPTIONS]

        d [DIR] cache directory (default: ./cache/)
        g       fix static GPS issues
        s       smoothen coordinates
        b       remove blank traces (trigger failure)
        r       remove stationary traces
        f       force regeneration of existing caches
        q       silence standard output
    """

optlist, fins = getopt.gnu_getopt(sys.argv[1:], 'd:gsbrfq')
optdict = dict(optlist)


if '-d' in optdict.keys():
    cache_dir = optdict['-d']
else:
    cache_dir = "cache"

if '-g' in optdict.keys():
    fix_gps = True
else:
    fix_gps = False

if '-s' in optdict.keys():
    smoothen_gps = True
else:
    smoothen_gps = False

if '-r' in optdict.keys():
    remove_stationary = True
else:
    remove_stationary = False

if '-b' in optdict.keys():
    remove_blanks = True
else:
    remove_blanks = False

if '-f' in optdict.keys():
    force_cache = True
else:
    force_cache = False

if '-q' in optdict.keys():
    be_quiet = True
else:
    be_quiet = False

try:
    survey_fnm = sys.argv[1]
except IndexError:
    print_syntax()
    sys.exit()

S = irlib.Survey(survey_fnm)

lines = S.GetLines()

for line in lines:
    line_no = line.split('_')[1]
    cache_fnm = S.GetLineCacheName(line_no, cache_dir=cache_dir)
    if os.path.isfile(cache_fnm) and not force_cache:
        pass
    else:
        if not be_quiet:
            print "Generating data for " + str(line) + " in " + cache_dir + "..."
        try:
            L = S.ExtractLine(line_no)
            L.RemoveBadLocations()
            if fix_gps:
                L.FixStaticGPS()
            if smoothen_gps:
                L.SmoothenGPS()
            if remove_blanks:
                L.RemoveBlankTraces()
            if remove_stationary:
                L.RemoveStationary(3.0)
            L.Dump(cache_fnm)
            del L
        except AttributeError:
            print "\tfailed"
        except KeyboardInterrupt:
            sys.exit(0)
