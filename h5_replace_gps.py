#! /usr/bin/python
#
#   Replace onboard GPS coordinates with data from a CSV file  - from NRCan PPP
#### modified from irlib h5_replace_gps.py
####  Derek Mueller April-June 2016
#

import sys
import getopt
import datetime
import shutil
import re

def print_syntax():
    print("""
    SYNTAX: h5_replace_ppp -h HDF_FILE -g CSV_FILE -o OUTFILE [OPTIONS]

    This tool replaces the existing geographical data in a ice radar HDF
    database with data taken from a CSV file, e.g. obtained from a PPP output
    of GPS data

    Required:

        -h s    (s) is an HDF dataset for which GPS timestamps exist

        -g s    (s) is a CSV file from NRCan PPP output

        -o s    (s) is the name of the output; if this file exists, it will be
                overwritten

        --tz n  (n) is the hour offset of the GPR computer from UTC

    Optional:

        -l n    Work only on line (n); default works on all lines

        -t n    Set the max time delta permissible for matching locations
                to (n) seconds; default is 15 seconds

        -n      Replace coordinates in HDF with no appropriate CSV counterpart
                with 'NaN'. By default, the original coordinates are retained.
    """)
    sys.exit(0)

def get_time(gps_timestamp, pcdatetime, tzoffset):
    """ Figure out the time and date from the HDF XML.

    This is trickier than it sounds, because the IceRadar records the PC time
    as well as the UTC hour, but not the UTC date. The correct UTC date is
    first calculated based on the time zone offset (*tzoffset*), and then the
    correct UTC time is substituted in.

    Returns a datetime.datetime object.
    """
    pcdate, pctime = pcdatetime.split("_")
    mm,dd,yy = pcdate.split("/")
    hms, ampm = pctime.split()
    h,m,s = hms.split(":")
    utcdate = datetime.datetime(int(yy), int(mm), int(dd), int(h), int(m), int(s)) \
            - datetime.timedelta(0, 3600*tzoffset)
    if ampm.lower() == "pm" and int(h) < 12:
        utcdate += datetime.timedelta(0, 3600*12)
    ts = gps_timestamp
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
        return (dt - datetime.datetime(1990, 1, 1)).total_seconds()
    except TypeError:
        return np.nan

def substituteXMLval(name, newval, xml):
    """ Replace the floating point "value" text in one of Blue System's IPR
    metadata fragments with *newval*.
    """
    newxml = re.sub(r'<Name>{0}</Name>[\r]?\n<Val>-?[0-9.]+?</Val>'.format(
                        name.replace(' ', '\s')),
                    r'<Name>{0}</Name>\r\n<Val>{1}</Val>'.format(
                        name, newval),
                    xml, flags=re.IGNORECASE)
    return newxml

optlist, fins = getopt.gnu_getopt(sys.argv[1:], 'i:h:g:o:l:t:n', ['tz='])
optdict = dict(optlist)

try:
    hdf_file = optdict["-h"]
    csv_file = optdict["-g"]
    outfile = optdict["-o"]
    tzoffset = int(optdict["--tz"])
except KeyError:
    print_syntax()

import numpy as np
import h5py
#import irlib.gpx

lineno = optdict.get("l", None)
max_dt = int(optdict.get("-t", 15))     # 15 second default
insert_nans = ("-n" in optdict)

print("Performing coordinate replacement")
print()
print("\t======== PARAMETERS =========\n"
      "\tSRC dataset:      {infile}\n"
      "\tDST dataset:      {outfile}\n"
      "\tCSV file:         {csv_file}\n\n"
      "\tMAX timedelta:    {max_dt} sec\n"
      "\tTZ offset:        {tz:+} hr\n"
      "\tINSERT NaNs:      {insert_nans}\n".format(infile=hdf_file, tz=tzoffset,
       csv_file=csv_file, max_dt=max_dt, outfile=outfile, insert_nans=insert_nans))

# Read in the CSV file

with open(csv_file) as f:
    lines = f.readlines()

ppplats = [float(a.split(",")[0]) for a in lines[1:]]
ppplons = [float(a.split(",")[1]) for a in lines[1:]]
pppeles = [float(a.split(",")[6]) for a in lines[1:]]
pppyear = [float(a.split(",")[5]) for a in lines[1:]]
pppdoy =  [float(a.split(",")[4]) for a in lines[1:]]
ppphr =  [float(a.split(",")[3]) for a in lines[1:]]

ppptimes = []
for i in range(len(pppyear)):
   ppptimes.append (datetime.datetime(int(pppyear[i]),1,1) + datetime.timedelta(int(pppdoy[i])-1) + datetime.timedelta(hours=ppphr[i]) )

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

# Interpolate the CSV positions
hdfseconds = np.array([dateseconds(a) for a in hdftimes])
pppseconds = np.array([dateseconds(a) for a in ppptimes])

def sortby(a, b):
    return [a_ for (b_, a_) in sorted(zip(b,a))]

ppplons = sortby(ppplons, pppseconds)
ppplats = sortby(ppplats, pppseconds)
pppeles = sortby(pppeles, pppseconds)
pppseconds.sort()

interp_lons = np.interp(hdfseconds, pppseconds, ppplons)
interp_lats = np.interp(hdfseconds, pppseconds, ppplats)
interp_eles = np.interp(hdfseconds, pppseconds, pppeles)

# Create a mask indicating where to use the interpolants
dts = np.array([np.min(np.abs(t-pppseconds)) for t in hdfseconds])
dt_mask = np.isfinite(dts) & (dts < max_dt)

# Replace the GPS data where the mask is true, otherwise put NaN
irep = 0
for i, dataset in enumerate(hdfaddrs):

    if insert_nans or dt_mask[i]:  # replace with better gps or NAN
        xml = dataset.attrs["GPS Cluster- MetaData_xml"]

        if dt_mask[i]:  # better gps

            xml = substituteXMLval("Long_ W", dec2dm(interp_lons[i]), xml)
            xml = substituteXMLval("Lat_N", dec2dm(interp_lats[i]), xml)
            xml = substituteXMLval("Alt_asl_m", str(interp_eles[i]), xml)

            irep += 1

        elif insert_nans:  # exceed max_dt so gps not good, replace with NAN

            xml = substituteXMLval("Long_ W", "NaN", xml)
            xml = substituteXMLval("Lat_N", "NaN", xml)
            xml = substituteXMLval("Alt_asl_m", "NaN", xml)

        dataset.attrs.modify("GPS Cluster- MetaData_xml", xml)

    else:
        pass

hdf.close()

print("\t========== RESULTS ==========\n"
      "\tTOTAL traces:     {ntraces}\n"
      "\tMODIFIED traces:  {irep}\n"
      "\tMEAN timedelta:   {mn:.2f} sec\n"
      "\tMEDIAN timedelta: {md:g} sec".format(irep=irep, ntraces=i+1,
      mn=np.mean(dts[np.isnan(dts) == False]), md=np.median(dts)))

