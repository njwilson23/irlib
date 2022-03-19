#! /usr/bin/env python
#
#   Add UTM coordinates to HDF5 radar data
#   Assumes that lat and lon can be positive or negative

import sys
import os
import shutil
import argparse
import traceback
import h5py
import irlib
import pyproj

def calculate_centroid(X, Y):
    """ Return the planar centroid of a matched pair of X and Y coordinates.
    """
    count_good = lambda L: sum((1 for a in L if a is not None))
    Xc = (a for a in X if a is not None)
    Yc = (a for a in Y if a is not None)
    xm = sum(Xc) / float(count_good(X))
    ym = sum(Yc) / float(count_good(Y))
    return (xm, ym)

def calculate_utm_zone(xll, yll):
    """ Determine the UTM zone that a lon-lat point falls in. Returns and
    integer and a string, either ('N') or ('S'). """
    if yll >= 0: hemi = 'N'
    else: hemi = 'S'
    zone = int((180 + xll) // 6) + 1
    return (zone, hemi)

# Command-line parsing
prog_description = """
    SYNTAX: h5_add_utm --swap_lon --swap_lat INFILE OUTFILE

        Replaces geographical coordinates in INFILE with UTM coordinates
        in OUTFILE. Does not perform any datum shift. 

	Works with 2 formats from BSI HDF files: 
  	  Old format - Latitude and longitude data in BSI HDF files are unsigned. It is assumed
        	to be in the western hemisphere by default. Passing the --swap_lon key
        	forces longitudes to be interpretted from the eastern hemisphere.
		UTM projection is calculated assuming that the data from neither from western Norway nor Svalbard.
	  New format - Latitude and longigude data in BSI HDF files are signed to indicate 
		hemisphere. If any lat or lon values are negative, the --swap_lon key and --swap_lat key is disabled
        """

parser = argparse.ArgumentParser(description = prog_description)
parser.add_argument("infile", help="input HDF (*.h5) filename, with or without path")
parser.add_argument("outfile")
parser.add_argument("--swap_lon",action='store_true')
parser.add_argument("--swap_lat",action='store_true')
args = parser.parse_args()

INFILE = args.infile
OUTFILE = args.outfile


print('operating on {0}'.format(INFILE))

# Open INFILE as an HDF5 dataset
if os.path.exists(INFILE):
    fin = h5py.File(INFILE, 'r')
else:
    print("No such file: {0}".format(INFILE))
    sys.exit(1)

# Get a list of all datasets and grab all metadata
print("querying input dataset...")
names = []
fin.visit(names.append)
datasets = [name for name in names if (isinstance(fin[name], h5py.Dataset)
                                       and 'picked' not in name)]

metadata = irlib.RecordList(fin)
print("reading metadata...")
failed = []
for i, dataset in enumerate(datasets):
    try:
        metadata.AddDataset(fin[dataset])
    except Exception as e:
        sys.stderr.write('Failed to read {0} due to {1}\n'.format(dataset, e))
        failed.append(i)
fin.close()
for i in failed[::-1]:
    del datasets[i]
print("\tdone")

# Based on centroid of survey, determine if this is old format or new format
xlm, ylm = calculate_centroid(metadata.lons, metadata.lats)
if xlm>0 and ylm>0:
    # This is Northern and Eastern Hemisphere... See if you should swap_lon
    if args.swap_lon:
        print("swapping sign on longitudes for eastern hemisphere")
        lons = metadata.lons
    else:
        lons = [-lon if lon is not None else None for lon in metadata.lons]
else:
    lons = metadata.lons
if args.swap_lat:
    lats = [-lat if lat is not None else None for lat in metadata.lats]
else:
    lats = metadata.lats    

num_sat = metadata.num_sat
fix_qual = metadata.fix_qual
gps_fix_valid = metadata.gps_fix_valid
gps_message_ok = metadata.gps_message_ok

# Project these to UTM using USGS proj.4
eastings = []
northings = []

xlm, ylm = calculate_centroid(lons, lats)
zone, hemi = calculate_utm_zone(xlm, ylm)

#projector = pyproj.Proj(proj='utm', zone=7, north=True)    # St Elias Range
#projector = pyproj.Proj(proj='utm', zone=16, north=True)   # Milne Ice Shelf

if hemi == 'N':
    projector = pyproj.Proj(proj='utm', zone=zone, north=True, datum="WGS84") # Auto-determined
if hemi == 'S':
    projector = pyproj.Proj(proj='utm', zone=zone, south=True, datum="WGS84") # Auto-determined

print("Projecting to UTM zone {0}{1}".format(zone, hemi))

for i, (lon, lat) in enumerate(zip(lons, lats)):
    if lon is not None and lat is not None:
        x, y = projector(lon, lat)
    else:
        x = 'NaN'
        y = 'NaN'
    eastings.append(x)
    northings.append(y)
print("\t{0} coordinate pairs projected".format(len([i for i in eastings if i != 'NaN'])))

print("writing new HDF5...")
# Copy INFILE to OUTFILE, and open in read/write mode
shutil.copyfile(INFILE, OUTFILE)
fout = h5py.File(OUTFILE, 'r+')

# For each dataset in OUTFILE, modify the UTM attribute cluster in place
for i, dataset in enumerate(datasets):
    try:
        try:   # This is the oldParseError way h5py library decodes based on data type specified
            xml = fout[dataset].attrs['GPS Cluster_UTM-MetaData_xml'].decode("utf-8")
        except AttributeError:  # This is the newer way, should work h5py >= 3.0 
            xml = fout[dataset].attrs['GPS Cluster_UTM-MetaData_xml']
        new_xml = (
                xml.replace('<Name>Datum</Name>\r\n<Val>NaN</Val>',
                            '<Name>Datum</Name>\r\n<Val>WGS84</Val>')
                .replace('<Name>Easting_m</Name>\r\n<Val></Val>',
                         '<Name>Easting_m</Name>\r\n<Val>{0}</Val>'
                            .format(eastings[i]))
                .replace('<Name>Northing_m</Name>\r\n<Val>NaN</Val>',
                         '<Name>Northing_m</Name>\r\n<Val>{0}</Val>'
                            .format(northings[i]))
                .replace('<Name>Zone</Name>\r\n<Val>NaN</Val>',
                         '<Name>Zone</Name>\r\n<Val>{0}</Val>'
                            .format(zone))
                .replace('<Name>GPS Fix Valid (dup)</Name>\r\n<Val></Val>',
                         '<Name>GPS Fix Valid (dup)</Name>\r\n<Val>{0}</Val>'
                            .format(gps_fix_valid[i]))
                .replace('<Name>GPS Message ok (dup)</Name>\r\n<Val></Val>',
                         '<Name>GPS Message ok (dup)</Name>\r\n<Val>{0}</Val>'
                            .format(gps_message_ok[i]))
                )       # num_sats appears to already be there
        fout[dataset].attrs.modify('GPS Cluster_UTM-MetaData_xml', new_xml.encode("utf-8"))
    except KeyError:
        # If the dataset doesn't have a UTM cluster, then add one
        utm_string = """<Cluster>
<Name>GPS_UTM Cluster</Name>
<NumElts>10</NumElts>
<String>
<Name>Datum</Name>
<Val>WGS84</Val>
</String>
<String>
<Name>Easting_m</Name>
<Val>{x}</Val>
</String>
<String>
<Name>Northing_m</Name>
<Val>{y}</Val>
</String>
<String>
<Name>Elevation</Name>
<Val>NaN</Val>
</String>
<String>
<Name>Zone</Name>
<Val>{zone}</Val>
</String>
<String>
<Name>Satellites (dup)</Name>
<Val>07</Val>
</String>
<Boolean>
<Name>GPS Fix Valid (dup)</Name>
<Val>{fix}</Val>
</Boolean>
<Boolean>
<Name>GPS Message ok (dup)</Name>
<Val>{ok}</Val>
</Boolean>
<Boolean>
<Name>Flag_1</Name>
<Val>0</Val>
</Boolean>
<Boolean>
<Name>Flag_2</Name>
<Val>0</Val>
</Boolean>
</Cluster>""".format(x=eastings[i], y=northings[i], zone=zone,
                     fix=gps_fix_valid[i], ok=gps_message_ok[i])
        fout[dataset].attrs['GPS Cluster_UTM-MetaData_xml'] = utm_string
    except:
        # Something else went wrong
        sys.stderr.write('problem creating {0}\n'.format(fout[dataset].name))
        traceback.print_exc()
fout.close()
