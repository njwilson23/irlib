#! /usr/bin/python
#
#   Replace onboard GPS coordinates with data from a GPX file
#
#   Updated by Derek Mueller April-Jun 2016 to handle CSV files
#   Further updated for more direct PPP support
#
#

from __future__ import print_function
import sys
import os
import glob
import getopt
import datetime
import shutil
import re
import numpy as np
import h5py
import pandas as pd
import irlib
import pdb
pdb.set_trace()

def print_syntax():
    print("""
    SYNTAX: h5_replace_gps -h HDF_FILE [-g GPX_FILE] [-c CSV_FILE] -o OUTFILE --tz OFFSET [OPTIONS]

    This tool replaces the existing geographical data in a ice radar HDF
    database with data taken from a GPX file, e.g. obtained from a handheld or
    external GPS unit or from a CSV file, e.g. obtained from a PPP output of
    GPS data

    Required:

        -h s    (s) is an HDF dataset for which GPS timestamps exist

        -g s    (s) is a GPS eXchange (GPX) file

        -c s    (s) is a CSV file, e.g. obtained from a PPP output of GPS data
                The file should have columns:
                latitude_decimal_degree,longitude_decimal_degree,ellipsoidal_height_m,decimal_hour,day_of_year,year,rcvr_clk_ns
                    
        -o s    (s) is the name of the output; if this file exists, it will be
                overwritten

        --tz n  (n) is the hour offset of the GPR computer from UTC

    Optional:

        -l n    Work only on line (n); default works on all lines

        -t n    Set the max time delta permissible for matching locations
                to (n) seconds; default is 15 seconds

        -n      Replace coordinates in HDF with no appropriate supplementary GPS 
                counterpart with 'NaN'. By default, the original coordinates 
                are retained.
                        
        -p      Keep all coordinates positive (use with old h5 format where Lat_N
                and Long_W)    
    """)
    sys.exit(0)

def get_time(gps_timestamp, timestamp, tzoffset, gpsmissing=False):
    """ 
    Figure out the time and date from the HDF XML.  First read the PC timestamp
       
    If gpsmissing is True, the pcdatetime is the only timestamp to use.  
    If gpsmissing is False, the gps_timestamp refines the pcdatetime (but only
    to the nearest second).
    
    Returns a datetime.datetime object.
    """
    utcdate = timestamp - datetime.timedelta(0, 3600*tzoffset) 
    if gpsmissing:
        return utcdate
    ts = gps_timestamp
    return datetime.datetime(utcdate.year, utcdate.month, utcdate.day,
                             int(ts[0:2]), int(ts[2:4]), int(ts[4:6]))

def dec2dm(dec_flt, signed=True):
    """ Convert a decimal degree coordinate to a degree decimal minute.
    Return a string. 
    Force all coordinates positive, only if signed is False    
    """
    
    if dec_flt<0:
        hem = -1
    else: 
        hem = 1        
    deg = int(abs(dec_flt) // 1)
    mm = round(abs(dec_flt) % 1 * 60.0, 4)
    if not signed:
        return str(deg) + str(int(mm//1)).rjust(2, '0') + "." + str(mm%1)[2:]
    else:
        return str(hem*deg) + str(int(mm//1)).rjust(2, '0') + "." + str(mm%1)[2:]

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
    newxml = re.sub(r'<Name>{0}</Name>[\r]?\n<Val>.*</Val>'.format(
                    name.replace(' ', '\s')),
                    r'<Name>{0}</Name>\r\n<Val>{1}</Val>'.format(
                    name, newval), xml, flags=re.IGNORECASE)
    
    return newxml    
    
def readppp(csv_file):
    '''
    This function reads in an NRCan PPP file (or files)

    Parameters
    ----------
    csv_file : str or list of strings
        file name or a pattern to match.  Example:
        ppp_out.csv  or *.csv

    Returns
    -------
    dataframe with gpstimes, gpslats, gpslons, gpseles as columns

    '''
    # Make single file into a list
    if type(csv_file) == str:
        csv_file = [csv_file]
    gps_dfs = []
    for file in csv_file:
        print('Reading {}'.format(file))
        try:
            gps_dfs.append(pd.read_csv(file))
            gps = pd.concat(gps_dfs, ignore_index=True)
        except:
            print("Could not read file {}".format(file))            
    
        # It may be that there are more columns - ortho ht (second last)
    gps.columns = ['gpslats','gpslons','gpseles','hr','doy','year','rcvr']
    
    gpstimes = []
    for i in range(len(gps.year)):
        gpstimes.append(datetime.datetime(int(gps.year[i]),1,1) +
                       datetime.timedelta(int(gps.doy[i])-1) +
                       datetime.timedelta(hours=gps.hr[i]))

    gps.insert(0,"gpstimes",gpstimes)
    gps = gps.drop(['hr','doy','year','rcvr'], axis=1)
    
    # makes sure there are no duplicate entries (takes first one) and sorts by time
    gps = gps.drop_duplicates(subset = ['gpstimes'])
    gps = gps.sort_values(by='gpstimes')

    return gps

def readgpx(gpx_file):
    """
    reads supplementary gps data from a GPX track file

    Parameters
    ----------
    gpx_file : str
        GPX track file

    Returns
    -------
    None.

    """

    import irlib.gpx  #leaving this here for now... 
    
    # Read in the GPX file
    gpstimes = []
    gpslons = []
    gpslats = []
    gpseles = []

    trackfile = irlib.gpx.GPX(gpx_file)
    for trk in trackfile.tracks:
        for trkseg in trk.trksegs:
            for pt in trkseg.trkpts:
                gpstimes.append(gpxtime2dt(pt.properties["time"]))
                gpslons.append(float(pt.lonlat[0]))
                gpslats.append(float(pt.lonlat[1]))
                gpseles.append(float(pt.properties["ele"]))
    
    gps = pd.DataFrame(list(zip(gpstimes,gpslats,gpslons,gpseles)),columns = 
                       ['gpstimes','gpslats','gpslons','gpseles'])
                       
    # makes sure there are no duplicate entries (takes first one) and sorts by time
    gps = gps.drop_duplicates(subset = ['gpstimes'])
    gps = gps.sort_values(by='gpstimes')
    
    return gps


optlist, fins = getopt.gnu_getopt(sys.argv[1:], 'i:h:g:c:o:l:t:n:p', ['tz='])
optdict = dict(optlist)

try:
    hdf_file = optdict["-h"]
    if "-g" in optdict:
        gpx_file = optdict["-g"]
        coordinate_source = "gpx"
    elif "-c" in optdict:
        if len(fins) != 0:
            csv_file = fins
        else:
            csv_file = optdict["-c"]
        coordinate_source = "csv"
    else:
        print("missing GPX or CSV file with GPS coordinates")
        raise KeyError()
    outfile = optdict["-o"]
    tzoffset = int(optdict["--tz"])
except KeyError:
    print_syntax()


lineno = optdict.get("l", None)
max_dt = int(optdict.get("-t", 15))     # 15 second default
insert_nans = ("-n" in optdict)
positive = ("-p" in optdict)

if coordinate_source == "csv":
    _gps_fnm = csv_file
else:
    _gps_fnm = gpx_file
    
print("Performing coordinate replacement")
print()
print("\t======== PARAMETERS =========\n"
      "\tSRC dataset:      {infile}\n"
      "\tDST dataset:      {outfile}\n"
      "\tGPS source:       {coordinate_source}\n"
      "\tGPS file:         {gps_fnm}\n\n"
      "\tMAX timedelta:    {max_dt} sec\n"
      "\tTZ offset:        {tz:+} hr\n"
      "\tINSERT NaNs:      {insert_nans}\n".format(infile=hdf_file, tz=tzoffset,
          coordinate_source=coordinate_source, gps_fnm=_gps_fnm, max_dt=max_dt,
          outfile=outfile, insert_nans=insert_nans))

if coordinate_source == "csv":
    gps = readppp(csv_file)

elif coordinate_source == "gpx":
    gps = readgpx(gpx_file)

### Start the interpolation and substitution.. 

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
            pcdatetime = dataset.attrs["PCSavetimestamp"]
            if len(pcdatetime.split(",")) == 4:
                timestamp, startbuf,buftime,pps = pcdatetime.split(",")
                timestamp = irlib.recordlist.pcdateconvert(timestamp, datefmt='ddmm')
            else:
                timestamp = irlib.recordlist.pcdateconvert(pcdatetime, datefmt='mmdd') # guessing the format      
            try:
                gpscluster = dataset.attrs['GPS Cluster- MetaData_xml']
                m = re.search(r'<Name>GPS_timestamp_UTC</Name>\r\n<Val>[0-9]{6}</Val>',
                              gpscluster)
                if m == None:
                    hdftimes.append(get_time(None,timestamp,tzoffset,gpsmissing=True))
                else:
                    gpstimestamp = re.search('[0-9]{6}', m.group()).group()
                    hdftimes.append(get_time(gpstimestamp,
                                         timestamp, tzoffset=0))  # GPS time IS in UTC
            except (ValueError, AttributeError):
                hdftimes.append(np.nan)

# Interpolate the supplementary GPS positions
hdfseconds = np.array([dateseconds(a) for a in hdftimes])
gpsseconds = np.array([dateseconds(a) for a in gps.gpstimes])

interp_lons = np.interp(hdfseconds, gpsseconds, gps.gpslons)
interp_lats = np.interp(hdfseconds, gpsseconds, gps.gpslats)
interp_eles = np.interp(hdfseconds, gpsseconds, gps.gpseles)

# Create a mask indicating where to use the interpolants
dts = np.array([np.min(np.abs(t-gpsseconds)) for t in hdfseconds])
dt_mask = np.isfinite(dts) & (dts < max_dt)

# Replace the GPS data where the mask is true, otherwise put NaN
irep = 0
for i, dataset in enumerate(hdfaddrs):

    if insert_nans or dt_mask[i]:  # replace with better gps or NAN
        xml = dataset.attrs["GPS Cluster- MetaData_xml"]

        if dt_mask[i]:  # better gps
            if positive:
                xml = substituteXMLval("Long_ W", dec2dm(interp_lons[i], signed=False), xml)
                xml = substituteXMLval("Lat_N", dec2dm(interp_lats[i], signed=False), xml)
                xml = substituteXMLval("Alt_asl_m", str(interp_eles[i], signed=False), xml)
            else:
                xml = substituteXMLval("Long", dec2dm(interp_lons[i]), xml)
                xml = substituteXMLval("Lat", dec2dm(interp_lats[i]), xml)
                xml = substituteXMLval("Alt_asl_m", str(interp_eles[i]), xml)

            irep += 1

        elif insert_nans:  # exceed max_dt so gps not good, replace with NAN
            if positive:
                xml = substituteXMLval("Long_ W", "NaN", xml)
                xml = substituteXMLval("Lat_N", "NaN", xml)
            else:
                xml = substituteXMLval("Long", "NaN", xml)
                xml = substituteXMLval("Lat", "NaN", xml)
                
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

