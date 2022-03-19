#! /usr/bin/env python
"""
   h5dumpmeta - HDF5 metadata dump

   Utility to write the metadata from raw radar data files to ASCII.
   Data is sent to stdout by default or to csv/shapefile depending on args
   
   Note that the csv file will output metadata for all traces.  
   The shapefile output (by necessity) will be only traces with lat/lon/elevation

"""
# standard libraries
import sys
import os.path
import glob
import traceback
import argparse

try:
    import StringIO ## for Python 2
except ImportError:
    import io as StringIO ## for Python 3

import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString

def meta2pd(infile):
    '''
    Takes a single h5 file and dumps the metadata to a pandas dataframe

    Parameters
    ----------
    infile : str
        An h5 file name, with extension (with path or in current directory)

    Returns
    -------
    a pandas dataframe with slightly different col names (10 char)

    '''
    stringbuffer = StringIO.StringIO()
    
    try:
        irlib.misc.ExtractAttrs(infile, fout=stringbuffer, eastern_hemisphere=False)
    except:
        sys.stderr.write("Error reading radar data\n")
        raise

    stringbuffer.seek(0)
    meta = pd.read_csv(stringbuffer,header=0)
    # change column headings to be 10 chars or less
    meta.columns = ['FID', 'filename', 'line', 'location', 'datacapt', 'echogram', \
                    'timestamp', 'lat', 'lon', 'gps_time', 'fix_qual', 'num_sat',
                    'dilution', 'alt_asl', 'geoid_ht', 'gps_fix', 'gps_ok',
                    'vrange', 'samplerate', 'stacking', 'trig_level', 'rec_len',
                    'startbuf', 'buftime', 'pps', 'datum', 'easting', 'northing',
                    'elevation', 'zone']
    meta = meta.sort_values('FID')
    
    #format FID
    fid = ['{:0>16}'.format(fid) for fid in meta.FID]
    meta.FID = fid
    
    return meta


## MAIN HERE  ###

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="input HDF (*.h5) filename, with or without path, if you use wildcards in linux, put this in quotes")
parser.add_argument("-o", "--outfile", help="output file BASENAME [if missing, will be automatically generated]")
parser.add_argument("-c", "--csv", help="create csv metadata file", action="store_true")
parser.add_argument("-w", "--wpt", help="create a waypoint metadata shapefile", action="store_true")
parser.add_argument("-l", "--line", help="create a line metadata shapefile", action="store_true")
parser.add_argument("--clobber", help="overwrite existing files", action="store_true")
args = parser.parse_args()

try:
    import irlib.misc       # Delay import for speed
except ImportError:
    print("Cannot find irlib libraries, check your path/environment")


infiles = glob.glob(args.infile)
infiles.sort()     
outfiles = [f.rsplit('.')[0] for f in infiles]

if len(infiles) > 1 and args.outfile:
    sys.stderr.write('You cannot name metadata from these {} input files {}\n'.format(
        len(infiles),args.outfile))
    sys.exit(2)
if len(infiles) == 1 and args.outfile:
    outfiles[0] = args.outfile


for i, infile in enumerate(infiles):
    sys.stderr.write('\n\n---\t'+infile+'\t---\n')
    
    if not os.path.isfile(infile):
        sys.stderr.write(infile + 'does not exist \n')
        sys.exit(1)
        
    meta = meta2pd(infile)
    # print to standard out - I commented this out since it is tedious
    #sys.stdout.write(meta.to_csv(sep=',', index=False))
    print('\n')
    print(meta)  # this will give users a sample of what is there without the long scroll
    print('\n')
    if args.csv: 
        if not os.path.isfile(outfiles[i]+'.csv'):
            meta.to_csv(outfiles[i]+'.csv', sep=',',index=False)
            sys.stderr.write(
                    "\t{fnm} written\n".format(
                    fnm=os.path.basename(outfiles[i]+".csv")))
        else:
            if args.clobber:
                meta.to_csv(outfiles[i]+'.csv', sep=',', index=False)
                sys.stderr.write(
                    "\t{fnm} overwritten\n".format(
                        fnm=os.path.basename(outfiles[i]+'.csv')))
            else:
                sys.stderr.write(
                    "\t{fnm} already exists\n".format(
                        fnm=os.path.basename(outfiles[i]+'.csv')))
       

    if args.wpt:
        # Create a shapefile for the metadata
        meta = meta[(meta.lon != "None")]
        meta = meta[(meta.lat != "None")]
        meta = meta[(meta.alt_asl != "None")]
        proj='EPSG:4326'  #Assuming WGS84
        #Creating a points while zipping 3 coordinates(3 dimension)
        if meta.shape[0]==0: 
            print('No valid location data found - cannot generate shapefile(s)')
            continue
        pts_gd=gpd.GeoDataFrame(meta,geometry=gpd.points_from_xy(
                meta.lon.astype(float), meta.lat.astype(float),z=meta.alt_asl, crs=proj))
                            
        if not os.path.isfile(outfiles[i]+"_wpt.shp"):
            pts_gd.to_file(outfiles[i]+"_wpt.shp")
            sys.stderr.write(
                    "\t{fnm} written\n".format(
                    fnm=os.path.basename(outfiles[i]+"_wpt.shp")))
        else:    
            if args.clobber:
                pts_gd.to_file(outfiles[i]+"_wpt.shp")
                sys.stderr.write(
                    "\t{fnm} overwritten\n".format(
                    fnm=os.path.basename(outfiles[i]+"_wpt.shp")))
            else:
                sys.stderr.write(
                   "\t{fnm} exists and was not overwritten\n".format(
                   fnm=os.path.basename(outfiles[i]+"_wpt.shp")))
                        
    if args.line:
        # Create a shapefile for the metadata - lines only
        meta = meta[(meta.lon != "None")]
        meta = meta[(meta.lat != "None")]
        meta = meta[(meta.alt_asl != "None")]
        if meta.shape[0] == 0: 
            print('No valid location data found - cannot generate shapefile(s)')
            continue       
        proj='EPSG:4326'  #Assuming WGS84
        #Creating a line while zipping 3 coordinates(3 dimension)
        lines=meta.line.unique() #finding all the unique object ids
        index_list=[]
        geometry_list=[]
        for ln in lines: #looping through lines
            ln_df= meta.loc[meta.line == ln] #subsetting pandas dataframe to the specific object
            line=LineString(zip(ln_df.lon.astype(float), ln_df.lat.astype(float), ln_df.alt_asl)) 
            line.wkt
            index_list.append(ln)
            geometry_list.append(line)
        line_gd=gpd.GeoDataFrame(index=index_list,crs=proj,geometry=geometry_list)
                    
        if not os.path.isfile(outfiles[i]+"_ln.shp"):
            line_gd.to_file(outfiles[i]+"_ln.shp")
            sys.stderr.write(
                "\t{fnm} written\n".format(
                    fnm=os.path.basename(outfiles[i]+"_ln.shp")))
        else:            
            if args.clobber:
                line_gd.to_file(outfiles[i]+"_ln.shp")
                sys.stderr.write(
                    "\t{fnm} overwritten\n".format(
                        fnm=os.path.basename(outfiles[i]+"_ln.shp")))
            else:
                sys.stderr.write(
                    "\t{fnm} exists and was not overwritten\n".format(
                    fnm=os.path.basename(outfiles[i]+"_ln.shp")))
                
                
