#! /usr/bin/env python
#
#   Replace onboard GPS coordinates with data from a GPX file
#
#   Updated by Derek Mueller April-Jun 2016 to handle CSV files
#   Further updated for more direct PPP support
#
#

from __future__ import print_function
import sys
import glob
import argparse
import datetime
import shutil
import re
import numpy as np
import h5py
import pandas as pd
import irlib

def get_time(gps_timestamp, timestamp, tzoffset, gpsmissing=False):
    """ 
    Figure out the time and date from the HDF XML.  First read the PC timestamp
    
    gps_timestamp is 'HHMMSS' or 'HHMMSS.SS' as a string
       
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
    csv_file : str
        file name or a pattern to match.  Example:
        ppp_out.csv  or *.csv

    Returns
    -------
    dataframe with gpstimes, gpslats, gpslons, gpseles as columns

    '''
    # Make single file into a list
    csv_files = glob.glob(csv_file)    
    gps_dfs = []
    for file in csv_files:
        print('Reading {}'.format(file))
        try:
            gps_dfs.append(pd.read_csv(file))
            gps = pd.concat(gps_dfs, ignore_index=True)
        except Exception as e:
            print("Could not read file {} due to {}".format(file, e))            
    print('\n')
        # It may be that there are more columns - ortho ht (second last)
    if gps.shape[1] == 7:
        gps.columns = ['gpslats','gpslons','gpseles','hr','doy','year','rcvr']
        print('Elevation is height above ellipsoid\n')
    if gps.shape[1] == 8:
        gps.columns = ['gpslats','gpslons','hae','hr','doy','year','gpseles','rcvr']
        print('Elevation is height above geoid \n')
              
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


## MAIN PROGRAM

prog_description = 'This tool replaces the existing geographical data in a ice radar HDF \
    database with data taken from a GPX file, e.g. obtained from a handheld or \
    external GPS unit or from a CSV file, e.g. obtained from a PPP output of GPS data'
prog_epilog = " Example: h5_replace_gps.py -n survey.h5 survey_ppp.h5 ppp*.csv ppp "
parser = argparse.ArgumentParser(description = prog_description, epilog=prog_epilog)
parser.add_argument("infile", help="input HDF (.h5) filename, with or without path, for which GPS or PC timestamps exist")
parser.add_argument("outfile", help="output HDF (.h5) filename, with or without path, if this file exists, it will be overwritten")
parser.add_argument("gpsfile", help="GPS filename(s), with enhanced location, with or without path / wildcards")
parser.add_argument('gpssource', choices=['gpx', 'ppp'], help="Select which format the gps file is in - either gpx or ppp")
parser.add_argument('timesource', choices=['iprgps', 'iprpc', 'both'], help="Select which timestamp to match gps timestamps to - iprgps (recommended), iprpc (if iprgps not available) or both (use caution)")
parser.add_argument("-t", "--tzoffset", help="is the hour offset (hh) of the GPR computer from UTC (default = 0)", default=0, type=float)
parser.add_argument("-l", "--line", help="Work only on line (n); default works on all lines", type=int)
parser.add_argument("-d", "--deltatimemax", help="Set the max time delta permissible for matching locations to (n) seconds; default is 15 seconds", 
                    default=15, type=int)
parser.add_argument("-o", "--offsetElev", help="Adds an offset to the elevations to account for the height of GPS off the ice or different geoid, use a neg. number to subtract.", 
                    default=0, type=float)
parser.add_argument("-n", "--replaceNaN", help="Replace coordinates in HDF with no appropriate supplementary GPS counterpart with 'NaN'. By default, the original coordinates are retained.", action="store_true")
parser.add_argument("-p", "--positivecoords", help="Keep all coordinates positive (use with old h5 format where Lat_N and Long_W)", action="store_true")

args = parser.parse_args()
  
print("Performing coordinate replacement\n\n")
print("\t======== PARAMETERS =========\n"
      "\tSRC dataset:      {infile}\n"
      "\tDST dataset:      {outfile}\n"
      "\tGPS source:       {gpssource}\n"
      "\tTime source:      {timesource}\n"
      "\tGPS file:         {gpsfile}\n\n"
      "\tDelta time max:   {max_dt} sec\n"
      "\tTZ offset:        {tzoffset:+} hr\n"
      "\tINSERT NaNs:      {insert_nans}\n"
      "\tOFFSET Elev:      {offsetElev}\n"
      "\tLine:             {line}\n"
      "\tPositive Coords:  {positive}\n".format(
          infile=args.infile, tzoffset=args.tzoffset,
          gpssource=args.gpssource, timesource=args.timesource,
          gpsfile=args.gpsfile, max_dt=args.deltatimemax,
          outfile=args.outfile, insert_nans=args.replaceNaN,line=args.line,
          offsetElev=args.offsetElev, positive=args.positivecoords))

max_dt = int(args.deltatimemax)     # 15 second default
insert_nans = args.replaceNaN
positive = args.positivecoords

if args.gpssource == "ppp":
    gps = readppp(args.gpsfile)

elif args.gpssource == "gpx":
    gps = readgpx(args.gpsfile)

if args.offsetElev:  # adds offset to GPS data 
    gps.gpseles = gps.gpseles + args.offsetElev

### Start the interpolation and substitution.. 

# Copy the HDF to the new file location, and then modify in-place
shutil.copyfile(args.infile, args.outfile)

# Load HDF file
hdf = h5py.File(args.outfile, "r+")

if args.line is None:
    lines = hdf.keys()
else:
    lines = ["line_"+args.line]

# Read out all acquisition times in the HDF file
hdfaddrs = []
hdfpctimes = []
hdfgpstimes = []
for line in lines:   # for every line... 
    if line == 'LabVIEW Boolean':
            continue  # this seems to be present in some files - not a line!
    for loc in hdf[line]:    #for every trace... 
        for dc in hdf[line][loc]: 

            dataset = hdf[line][loc][dc]["echogram_0"]
            hdfaddrs.append(dataset)
            # get the timestamp from the EPU computer
            pcdatetime = dataset.attrs["PCSavetimestamp"]
            if not type(pcdatetime) is str:   #  I lost track of what this may be if not a string! 
                pcdatetime = pcdatetime.astype(str) # make sure it is a string
            if len(pcdatetime.split(",")) == 4:    # this is a new (> 2016 file format)
                timestamp, startbuf,buftime,pps = pcdatetime.split(",")
                timestamp = irlib.recordlist.pcdateconvert(timestamp, datefmt='ddmm')
            elif len(pcdatetime.split(",")) == 3:    # this is a format where the date is in a comment field
                timestamp = irlib.recordlist.TimeFromComment(args.outfile, line, loc)
            else:   # this is an older file format  
                timestamp = irlib.recordlist.pcdateconvert(pcdatetime, datefmt='mmdd') # guessing the format      
            
            # add this to list
            hdfpctimes.append(get_time(None,timestamp,args.tzoffset,gpsmissing=True))
            
            try:
                gpscluster = dataset.attrs['GPS Cluster- MetaData_xml']
                # search for a timestamp from the onboard GPS
                # These vary according to the IPR hdf version - HHMMSS or HHMMSS.SS
                m = re.search(r'<Name>GPS_timestamp_UTC</Name>\r\n<Val>[0-9]{6}\.?\d*</Val>',
                              gpscluster)
                if m == None:  
                    # if there is no GPS timestamp - use the one from the EPU computer
                    hdfgpstimes.append(np.nan)
                else:
                    # if there is a GPS timestamp use that instead
                    gpstimestamp = re.search('[0-9]{6}\.?\d*', m.group()).group()
                    hdfgpstimes.append(get_time(gpstimestamp,
                                         timestamp, 0))  # GPS time IS in UTC so tzoffset must be 0
            except (ValueError, AttributeError):
                hdfgpstimes.append(np.nan)

# convert all datetimes to seconds
hdfgpsseconds = np.array([dateseconds(a) for a in hdfgpstimes])  # represents the timestamps from IPR GPS
hdfpcseconds = np.array([dateseconds(a) for a in hdfpctimes])  # represents the timestamps from  EPU computer
gpsseconds = np.array([dateseconds(a) for a in gps.gpstimes])  # represents the timestamps from 'better' GPS file provided

#DM compared hdfgpstimes and hdfpctimes for some files: 
# What is the offset?  Can you use pc time instead of gps time when it is missing? 
# I would say no.  The difference between these is variable and there will be an offset
# between gps and pc time unless the pc was just set to UTC... 
#np.nanmean(hdfgpsseconds-hdfpcseconds)
#np.nanmax(hdfgpsseconds-hdfpcseconds)
#np.nanmin(hdfgpsseconds-hdfpcseconds)
#np.nanmedian(hdfgpsseconds-hdfpcseconds)
#np.nanstd(hdfgpsseconds-hdfpcseconds)
# You CAN use pctimes instead of the gps times but this is not recommended (see timesource) 
# This may be improved if the pc times are lagged accordingly

if args.timesource == 'iprgps':
    hdfseconds = hdfgpsseconds
elif args.timesource == 'pcgps':    
    hdfseconds = hdfpcseconds
else:
    hdfseconds = hdfgpsseconds
    # if the hdfgps data is nan, then use the other source
    hdfseconds[np.isnan(hdfseconds)] = hdfpcseconds[np.isnan(hdfgpsseconds)]


# Check here to see if the timestamps worked out: 
if sum(np.isnan(np.array(hdfseconds))) == len(hdfseconds):
    print(''''All IPR onboard timestamps are missing. Check that you have the correct timesource setting.  
          If you expect them to be there, there might be an issue with the code.''')
    sys.exit(1)
     
# Interpolate the supplementary GPS positions
#One-dimensional linear interpolation for monotonically increasing sample points.
#Returns the one-dimensional piecewise linear interpolant to a function with given discrete data points 
interp_lons = np.interp(hdfseconds, gpsseconds, gps.gpslons)
interp_lats = np.interp(hdfseconds, gpsseconds, gps.gpslats)
interp_eles = np.interp(hdfseconds, gpsseconds, gps.gpseles)

# Determine the time difference between supplementary gps and hdf timestamps for each trace
dts = np.array([np.min(np.abs(t-gpsseconds)) for t in hdfseconds])
# Create a mask indicating where to use the interpolants 
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
                xml = substituteXMLval("Alt_asl_m", str(interp_eles[i]), xml)
            else:
                xml = substituteXMLval("Long", dec2dm(interp_lons[i]), xml)
                xml = substituteXMLval("Lat", dec2dm(interp_lats[i]), xml)
                xml = substituteXMLval("Alt_asl_m", str(interp_eles[i]), xml)

            irep += 1

        elif insert_nans:  # exceed max_dt so gps not good, replace with NAN
            if positive: 
                xml = substituteXMLval("Long_ W", "NaN", xml)
                xml = substituteXMLval("Lat_N", "NaN", xml)
                xml = substituteXMLval("Alt_asl_m", "NaN", xml)
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
      "\t MEAN timedelta:   {mn:.2f} sec\n"
      "\t MEDIAN timedelta: {md:g} sec\n"
      "\t MAX timedelta: {ma:g} sec\n"
      "\tMODIFIED traces:  {irep}\n".format(ntraces=i+1,irep=irep, 
        mn=np.nanmean(dts), md=np.nanmedian(dts), ma=np.nanmax(dts)))
if not len(dts[dt_mask]) == 0:
    print("\t MEAN timedelta (modified):   {mnm:.2f} sec\n"
      "\t MEDIAN timedelta  (modified): {mdm:g} sec\n"
      "\t MAX timedelta  (modified): {mam:g} sec\n".format(
          mnm=np.nanmean(dts[dt_mask]),
          mdm=np.nanmedian(dts[dt_mask]),
          mam=np.nanmax(dts[dt_mask])))    
