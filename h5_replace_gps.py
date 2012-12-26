#! /usr/bin/python
#
#   Replace onboard GPS coordinates on a line with Magellan data
#
#   Should be adapted to work when some onboard GPS data are good
#   Use '-f' flag to control overwriting
#

import datetime, time
import getopt
import sys
import shutil
import re
sys.path.append("/home/njw23/python")
import pdb

def print_syntax():
    print """
    SYNTAX: h5_replace_gps -h HDF_FILE -g GPX_FILE -o OUTFILE [OPTIONS]

    Required:
        -h s    (s) is an HDF dataset for which GPS timestamps exist
        -g s    (s) is a GPS interchange file downloaded from a hand-held GPS
        -o s    (s) is the name of the output; if this file exists, it will be
                overwritten

    Valid options:
        -l n    Work only on line (n); default works on all lines
        -t n    Set the max time delta permissible for matching locations
                to (n) seconds; default is 15 seconds
        -i n    Interpolate GPS times to (n) second interval; default is,
                no interpolation
                Interpolation is done without projection under the
                assumption that over small distances the difference between
                a small circle and a great circle is negligible.
                For reasons mysterious, interpolation is very CPU-intensive.
    """
    sys.exit(0)

def get_time_lon_lat(xml, savetimestamp):
    """ Return a tuple with time, lon, and lat for a location. """
    dd,mm,yy = savetimestamp.split("_")[0].split("/")
    m = re.search(r'<Name>GPS_timestamp_UTC</Name>\r\n<Val>[0-9]+.[0-9]+</Val>', xml)
    try:
        ts = re.search('[0-9]+.[0-9]+', m.group()).group()
        t = datetime.datetime(int(yy), int(mm), int(dd), int(ts[0:2]), int(ts[2:4]), int(ts[4:6]))
    except ValueError:      # There's a bad loc
        t = None
    except AttributeError:  # re.search returned None
        t = None

    m = re.search(r'<Name>Long_ W</Name>\r\n<Val>[0-9]+.[0-9]+</Val>', xml)
    if m is not None:
        lon = dm2dec(re.search('[0-9]+.[0-9]+', m.group()).group())
        if lon is not None:
            lon = -lon
    else:
        lon = None

    m = re.search(r'<Name>Lat_N</Name>\r\n<Val>[0-9]+.[0-9]+</Val>', xml)
    if m is not None:
        lat = dm2dec(re.search('[0-9]+.[0-9]+', m.group()).group())
    else:
        lat = None

    return (t, lon, lat)

def dm2dec(dmstr):
    """ Convert the degree - decimal minute codes in radar data
    to a decimal degree coordinate. dmstr is expected to a string.
    """
    if dmstr == '': return
    try:
        return round(float(dmstr.split('.')[0][:-2])
        +   float(dmstr[-7:]) / 60., 6)
    except ValueError:
        return None

def dec2dm(dec_flt):
    """ Convert a decimal degree coordinate to a degree decimal minute.
    Make the coordinate positive. Return a string.
    """
    deg = int(abs(dec_flt) // 1)
    mm = round(abs(dec_flt) % 1 * 60.0, 4)
    return str(deg) + str(int(mm//1)).rjust(2, '0') + "." + str(mm%1)[2:]

def gpx_str2datetime(s):
    """ Take a GPX timestamp and return a datetime.datetime() object. """
    hh, mn, ss = s.split('T')[1].strip('Z').split(':')
    yy, mm, dd = s.split('T')[0].split('-')
    return datetime.datetime(int(yy), int(mm), int(dd), int(hh), int(mn), int(round(float(ss))))

def sortby(A, B):
    """ Sort items in list A using items in B as a sorting key. """
    C = zip(B, A)
    C.sort()
    return [i[1] for i in C]

optlist, fins = getopt.gnu_getopt(sys.argv[1:], 'i:h:g:o:l:t:')
optdict = dict(optlist)

try:
    hdf_file = optdict["-h"]
    gpx_file = optdict["-g"]
    outfile = optdict["-o"]
    # Delay loading modules so it responds faster to bad arguments
    import numpy as np
    import irlib
    import h5py
    from geo_tools.vector.gpxparser import GPXParser
except KeyError:
    print_syntax()

if "-l" in optdict.keys():
    lineno = optdict["-l"]
else:
    lineno = None

if "-t" in optdict.keys():
    max_dt = datetime.timedelta(0, int(optdict["-t"]))
else:
    max_dt = datetime.timedelta(0, 15)      # 15 seconds default

if "-i" in optdict.keys():
    interp_dt = int(optdict["-i"])
    import numpy as np
else:
    interp_dt = None

print "Replacing coordinates in {infile}\n".format(infile=hdf_file) + \
      "with those from {gpx_file}\n".format(gpx_file=gpx_file) + \
      "where the difference is less than {max_dt}\n".format(max_dt=max_dt) + \
      "and saving as {outfile}".format(outfile=outfile)

# Copy the HDF to the new file location, and then modify in-place
shutil.copyfile(hdf_file, outfile)

# Load GPX file
gpx = GPXParser(gpx_file)

# Load hdf file
hdf = h5py.File(outfile, "r+")

# Read in the GPX file
gpx_times = []
gpx_lons = []
gpx_lats = []
for name in gpx.getnames():
    t = gpx.tracks[name].keys()
    _ = [gpx_times.append(gpx_str2datetime(i)) for i in t]
    _ = [gpx_lons.append(gpx.tracks[name][i]['lon']) for i in t]
    _ = [gpx_lats.append(gpx.tracks[name][i]['lat']) for i in t]

if interp_dt is not None:
    # Find the min and max times, and set interpolation points between them
    gpx_times = map(lambda t: time.mktime(t.timetuple()), gpx_times)
    min_t = min(gpx_times)
    max_t = max(gpx_times)
    gpx_lons = sortby(gpx_lons, gpx_times)
    gpx_lats = sortby(gpx_lats, gpx_times)
    gpx_times.sort()
    gpx_times = np.array(gpx_times)
    interp_t = np.arange(min_t, max_t, interp_dt)
    all_t = np.hstack([interp_t, gpx_times])
    gpx_lons = np.interp(all_t, gpx_times, gpx_lons)
    gpx_lats = np.interp(all_t, gpx_times, gpx_lats)
    gpx_times = np.array(map(lambda ts: datetime.datetime.fromtimestamp(ts), all_t))

else:
    gpx_times = np.array(gpx_times)

# For each location in the hdf file, read the timestamp, and figure out
# whether there's an appropriate position in GPX
if lineno is None:
    lines = hdf.keys()
else:
    lines = ["line_"+lineno]

nans = 0
replaces = 0
DT_rec = []

for line in lines:

    for loc in hdf[line].keys():

        dataset = hdf[line][loc]["datacapture_0/echogram_0"]
        xml = dataset.attrs['GPS Cluster- MetaData_xml']
        savetimestamp = dataset.attrs['PCSavetimestamp']

        # Read the timestamp, lon, and lat from xml
        hdf_t, hdf_lon, hdf_lat = get_time_lon_lat(xml, savetimestamp)

        if None not in (hdf_t, hdf_lon, hdf_lat):
            # Find the nearest observation in the GPX data
            dt = gpx_times - hdf_t
            nearest_i = np.nonzero(abs(dt)==abs(dt).min())[0][0]
            dt_i = abs(dt).min()
            DT_rec.append(dt_i)

            if dt_i < max_dt:
                # Insert the GPX coordinates
                gpx_lon = gpx_lons[nearest_i]
                gpx_lat = gpx_lats[nearest_i]
                aa = float(dec2dm(gpx_lon))
                bb = float(dec2dm(gpx_lat))
                if (aa < 13900) or (aa > 14000) or (bb < 6000) or (bb > 6100):
                    print gpx_lon, gpx_lat, aa, bb
                xml = re.sub(
                    r'<Name>Long_ W</Name>\r\n<Val>([0-9][0-9]+\.[0-9]+?)</Val>',
                    r'<Name>Long_ W</Name>\r\n<Val>{0}</Val>'.format(dec2dm(gpx_lon)),
                    xml)
                xml = re.sub(
                    r'<Name>Lat_N</Name>\r\n<Val>([0-9][0-9]+\.[0-9]+?)</Val>',
                    r'<Name>Lat_N</Name>\r\n<Val>{0}</Val>'.format(dec2dm(gpx_lat)),
                    xml)
                replaces += 1

            else:
                # Insert NaN
                xml = re.sub(
                    r'<Name>Long_ W</Name>\r\n<Val>([0-9][0-9]+\.[0-9]+?)</Val>',
                    r'<Name>Long_ W</Name>\r\n<Val>{0}</Val>'.format("NaN"),
                    xml)
                xml = re.sub(
                    r'<Name>Lat_N</Name>\r\n<Val>([0-9][0-9]+\.[0-9]+?)</Val>',
                    r'<Name>Lat_N</Name>\r\n<Val>{0}</Val>'.format("NaN"),
                    xml)
                nans += 1

        dataset.attrs.modify('GPS Cluster- MetaData_xml', xml)

hdf.close()
print "{0} matches found out of {1} points".format(replaces, replaces+nans)
DT_rec = [i.seconds for i in DT_rec]
print "Mean offset: {0} sec".format(float(sum(DT_rec)) / float(len(DT_rec)))
DT_rec.sort()
print "Median offset: {0} sec".format(DT_rec[len(DT_rec)//2])


