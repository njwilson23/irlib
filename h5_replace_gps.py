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

def print_syntax():
    print """
    SYNTAX: h5_replace_gps -h HDF_FILE -g GPX_FILE -o OUTFILE [OPTIONS]

    This tool replaces the existing geographical data in a ice radar HDF
    database with data taken from a GPX file, e.g. obtained from a handheld or
    external GPS unit.

    Required:

        -h s    (s) is an HDF dataset for which GPS timestamps exist

        -g s    (s) is a GPS interchange (GPX)

        -o s    (s) is the name of the output; if this file exists, it will be
                overwritten

        --tz n  (n) is the hour offset of the GPR computer from UTC

    Optional:

        -l n    Work only on line (n); default works on all lines

        -t n    Set the max time delta permissible for matching locations
                to (n) seconds; default is 15 seconds

        -n      Replace coordinates in HDF with no appropriate GPX counterpart
                with 'NaN'. By default, the original coordinates are retained.
    """
    sys.exit(0)

def get_time(gps_timestamp, pctime, tzoffset):
    """ Figure out the time and date from the HDF XML. 
    
    This is trickier than it sounds, because the IceRadar records the PC time
    as well as the UTC hour, but not the UTC date. The correct UTC date is
    first calculated based on the time zone offset (*tzoffset*), and then the
    correct UTC time is substituted in.

    Returns a datetime.datetime object.
    """
    mm,dd,yy = pctime.split("_")[0].split("/")
    ts = gps_timestamp
    utcdate = datetime.datetime(int(yy), int(mm), int(dd)) \
            - datetime.timedelta(0, 3600*tzoffset)
    return datetime.datetime(utcdate.year, utcdate.month, utcdate.day,
                             int(ts[0:2]), int(ts[2:4]), int(ts[4:6]))

def dec2dm(dec_flt):
    """ Convert a decimal degree coordinate to a degree decimal minute.
    Make the coordinate positive. Return a string. """
    deg = int(abs(dec_flt) // 1)
    mm = round(abs(dec_flt) % 1 * 60.0, 4)
    return str(deg) + str(int(mm//1)).rjust(2, '0') + "." + str(mm%1)[2:]

def gpxtime2dt(s):
    """ Take a GPX timestamp and return a datetime.datetime object. """
    hh, mn, ss = s.split('T')[1].strip('Z').split(':')
    YY, MM, DD = s.split('T')[0].split('-')
    return datetime.datetime(int(YY), int(MM), int(DD), int(hh), int(mn), int(ss))

def dateseconds(dt):
    """ Return a date in seconds from a datetime. """
    try:
        return (datetime.datetime(1990, 1, 1) - dt).seconds
    except TypeError:
        return np.nan

optlist, fins = getopt.gnu_getopt(sys.argv[1:], 'i:h:g:o:l:t:n', ['tz='])
optdict = dict(optlist)

try:
    hdf_file = optdict["-h"]
    gpx_file = optdict["-g"]
    outfile = optdict["-o"]
    tzoffset = int(optdict["--tz"])
except KeyError:
    print_syntax()

import numpy as np
import h5py
import irlib.gpx

lineno = optdict.get("l", None)
max_dt = int(optdict.get("-t", 15))     # 15 second default
insert_nans = ("-n" in optdict)

print ("Replacing coordinates in {infile} (UTC {tz:+d})\n"
        "\twith those from {gpx_file}\n"
        "\twhere the time delta is < {max_dt} seconds\n"
        "\tand saving as {outfile}".format(infile=hdf_file, tz=tzoffset,
        gpx_file=gpx_file, max_dt=max_dt, outfile=outfile))

# Read in the GPX file
gpxtimes = []
gpxlons = []
gpxlats = []
gpxeles = []

trackfile = irlib.gpx.GPX(gpx_file)
for trk in trackfile.tracks.values():
    for trkseg in trk.segments:
        for xy, z, t in zip(trkseg.vertices, trkseg.data['ele'], trkseg.data['time']):
            gpxtimes.append(gpxtime2dt(t))
            gpxlons.append(float(xy[0]))
            gpxlats.append(float(xy[1]))
            gpxeles.append(float(z))

# Copy the HDF to the new file location, and then modify in-place
shutil.copyfile(hdf_file, outfile)

# Load HDF file
hdf = h5py.File(outfile, "r+")

if lineno is None:
    lines = hdf.keys()
else:
    lines = ["line_"+lineno]

# Read out all acquisition times in the HDF file
hdfaddrs = []
hdftimes = []
for line in lines:
    for loc in hdf[line]:
        for dc in hdf[line][loc]:


            dataset = hdf[line][loc][dc]["echogram_0"]
            hdfaddrs.append(dataset)

            try:

                gpscluster = dataset.attrs['GPS Cluster- MetaData_xml']
                m = re.search(r'<Name>GPS_timestamp_UTC</Name>\r\n<Val>[0-9]{6}</Val>',
                              gpscluster)
                gpstimestamp = re.search('[0-9]{6}', m.group()).group()
                hdftimes.append(get_time(gpstimestamp,
                                         dataset.attrs["PCSavetimestamp"],
                                         tzoffset))

            except (ValueError, AttributeError):
                hdftimes.append(np.nan)

# Interpolate the GPX positions
hdfseconds = np.array([dateseconds(a) for a in hdftimes])
gpxseconds = np.array([dateseconds(a) for a in gpxtimes])
interp_lons = np.interp(hdfseconds, gpxseconds, gpxlons)
interp_lats = np.interp(hdfseconds, gpxseconds, gpxlats)
interp_eles = np.interp(hdfseconds, gpxseconds, gpxeles)

# Create a mask indicating where to use the interpolants
dts = np.array([np.min(np.abs(t-gpxseconds)) for t in hdfseconds])
dt_mask = dts < max_dt


# Replace the GPS data where the mask is true
for i, dataset in enumerate(hdfaddrs):

    if insert_nans or dt_mask[i]:
        xml = dataset.attrs["GPS Cluster- MetaData_xml"]

        if dt_mask[i]:
            xml = re.sub(r"<Name>Long_ W</Name>\r\n<Val>[0-9.]+?</Val>",
                         r"<Name>Long_ W</Name>\r\n<Val>{0}</Val>"
                         .format(dec2dm(interp_lons[i])), xml)

            xml = re.sub(r"<Name>Lat_N</Name>\r\n<Val>[0-9.]+?</Val>",
                         r"<Name>Lat_N</Name>\r\n<Val>{0}</Val>"
                         .format(dec2dm(interp_lats[i])), xml)
            
            xml = re.sub(r"<Name>Alt_asl_m</Name>\r\n<Val>[0-9.]+?</Val>",
                         r"<Name>Alt_asl_m</Name>\r\n<Val>{0}</Val>"
                         .format(interp_eles[i]), xml)

        elif insert_nans:
            xml = re.sub(r"<Name>Long_ W</Name>\r\n<Val>[0-9.]+?</Val>",
                         r"<Name>Long_ W</Name>\r\n<Val>NaN</Val>", xml)

            xml = re.sub(r"<Name>Lat_N</Name>\r\n<Val>[0-9.]+?</Val>",
                         r"<Name>Lat_N</Name>\r\n<Val>NaN</Val>", xml)
            
            xml = re.sub(r"<Name>Alt_asl_m</Name>\r\n<Val>[0-9.]+?</Val>",
                         r"<Name>Alt_asl_m</Name>\r\n<Val>NaN</Val>", xml)

        dataset.attrs.modify("GPS Cluster- MetaData_xml", xml)

    else:
        pass

hdf.close()

print "Mean/median time deltas:"
print "\t{0} / {1}".format(np.mean(dts), np.median(dts))

