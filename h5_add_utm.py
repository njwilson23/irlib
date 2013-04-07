#! /usr/bin/env python
#
#   Add UTM coordinates to HDF5 radar data
#

import sys
import os
import shutil
import traceback

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

if len(sys.argv) < 3:
    print """
    SYNTAX: h5_add_utm INFILE OUTFILE

        Replaces geographical coordinates in INFILE with UTM coordinates
        in OUTFILE. Does not perform any datum shift. Projection is calculated
        assuming that the data from neither from western Norway nor Svalbard.
    """
    sys.exit(0)
else:
    import irlib
    from irlib.recordlist import ParseError
    import h5py
    import pyproj
    INFILE = sys.argv[1]
    OUTFILE = sys.argv[2]

print 'operating on {0}'.format(INFILE)

# Open INFILE as an HDF5 dataset
if os.path.exists(INFILE):
    fin = h5py.File(INFILE, 'r')
else:
    print "No such file: {0}".format(INFILE)
    sys.exit(0)

# Get a list of all datasets and grab all metadata
print "querying input dataset..."
names = []
fin.visit(names.append)
datasets = [name for name in names if (isinstance(fin[name], h5py.Dataset)
                                       and 'picked' not in name)]

metadata = irlib.RecordList(fin)
print "reading metadata..."
failed = []
for i, dataset in enumerate(datasets):
    try:
        metadata.AddDataset(fin[dataset])
    except ParseError:
        sys.stderr.write('Failed to read {0}\n'.format(dataset))
        failed.append(i)
fin.close()
for i in failed[::-1]:
    del datasets[i]
print "\tdone"

# Pull out the geographical data for convenience
lons = metadata.lons
lats = metadata.lats
num_sat = metadata.num_sat
fix_qual = metadata.fix_qual
gps_fix_valid = metadata.gps_fix_valid
gps_message_ok = metadata.gps_message_ok

# Project these to UTM using USGS proj.4
eastings = []
northings = []

xlm, ylm = calculate_centroid(lons, lats)
zone, hemi = calculate_utm_zone(-xlm, ylm)
if hemi == 'N': north = True
else: north = False
#projector = pyproj.Proj(proj='utm', zone=7, north=True)    # St Elias Range
#projector = pyproj.Proj(proj='utm', zone=16, north=True)   # Milne Ice Shelf
projector = pyproj.Proj(proj='utm', zone=zone, north=north) # Auto-determined
print "Projecting to UTM zone {0}{1}".format(zone, hemi)

for i, (lon, lat) in enumerate(zip(lons, lats)):
    if lon is not None and lat is not None:
        x, y = projector(-lon, lat)
        #if (x < 600000) or (x > 700000) or (y < 6000000) or (y > 7000000):
        #    print "Warning:",i,'\t',-lon,'\t',lat,'\t',x,'\t',y
    else:
        x = 'NaN'
        y = 'NaN'
    eastings.append(x)
    northings.append(y)
print "{0} coordinate pairs projected".format(sum([1 for i in eastings if i != 'NaN']))

print "writing new HDF5..."
# Copy INFILE to OUTFILE, and open in read/write mode
shutil.copyfile(INFILE, OUTFILE)
fout = h5py.File(OUTFILE, 'r+')

# For each dataset in OUTFILE, modify the UTM attribute cluster in place
for i, dataset in enumerate(datasets):
    try:
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
                         '<Name>Zone</Name>\r\n<Val>7</Val>')
                .replace('<Name>GPS Fix Valid (dup)</Name>\r\n<Val></Val>',
                         '<Name>GPS Fix Valid (dup)</Name>\r\n<Val>{0}</Val>'
                            .format(gps_fix_valid[i]))
                .replace('<Name>GPS Message ok (dup)</Name>\r\n<Val></Val>',
                         '<Name>GPS Message ok (dup)</Name>\r\n<Val>{0}</Val>'
                            .format(gps_message_ok[i]))
                )       # num_sats appears to already be there
        fout[dataset].attrs.modify('GPS Cluster_UTM-MetaData_xml', new_xml)
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
<Val>7</Val>
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
</Cluster>""".format(x=eastings[i], y=northings[i],
                     fix=gps_fix_valid[i], ok=gps_message_ok[i])
        fout[dataset].attrs['GPS Cluster_UTM-MetaData_xml'] = utm_string
    except:
        # Something else went wrong
        sys.stderr.write('problem creating {0}\n'.format(fout[dataset].name))
        traceback.print_exc()
fout.close()
