#! /usr/bin/env python
#
#   Add UTM coordinates to HDF5 radar data
#

import sys, shutil
import traceback

if len(sys.argv) < 3:
    print """
    SYNTAX: h5_add_utm INFILE OUTFILE

        Replaces geographical coordinates in INFILE with UTM coordinates
        in OUTFILE. Does not perform any datum shift. Projection assumed
        to be UTM zone 7N.
    """
    sys.exit(0)
else:
    import irlib
    import h5py
    import pyproj
    INFILE = sys.argv[1]
    OUTFILE = sys.argv[2]

print 'operating on {0}'.format(INFILE)

# Open INFILE as an HDF5 dataset
fin = h5py.File(INFILE, 'r')

# Get a list of all datasets and grab all metadata
print "querying input dataset..."
names = []
fin.visit(names.append)
datasets = [name for name in names if (isinstance(fin[name], h5py.Dataset) and 'picked' not in name)]

metadata = irlib.RecordList(fin)
print "reading metadata..."
for dataset in datasets:
    metadata.AddDataset(fin[dataset])
fin.close()
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
projector = pyproj.Proj(proj='utm', zone=7, north=True)
for i, (lon, lat) in enumerate(zip(lons, lats)):
    if lon is not None and lat is not None:
        x, y = projector(-lon, lat)
        if (x < 600000) or (x > 700000) or (y < 6000000) or (y > 7000000):
            print "Warning:",i,'\t',-lon,'\t',lat,'\t',x,'\t',y
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
                xml.replace('<Name>Datum</Name>\r\n<Val>NaN</Val>', '<Name>Datum</Name>\r\n<Val>WGS84</Val>')
                .replace('<Name>Easting_m</Name>\r\n<Val></Val>', '<Name>Easting_m</Name>\r\n<Val>{0}</Val>'.format(eastings[i]))
                .replace('<Name>Northing_m</Name>\r\n<Val>NaN</Val>', '<Name>Northing_m</Name>\r\n<Val>{0}</Val>'.format(northings[i]))
                .replace('<Name>Zone</Name>\r\n<Val>NaN</Val>', '<Name>Zone</Name>\r\n<Val>7</Val>')
                .replace('<Name>GPS Fix Valid (dup)</Name>\r\n<Val></Val>', '<Name>GPS Fix Valid (dup)</Name>\r\n<Val>{0}</Val>'.format(gps_fix_valid[i]))
                .replace('<Name>GPS Message ok (dup)</Name>\r\n<Val></Val>', '<Name>GPS Message Valid (dup)</Name>\r\n<Val>{0}</Val>'.format(gps_message_ok[i]))
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
        sys.stderr.write('problem getting attributes from {0}\n'.format(f[dataset].name))
        traceback.print_exc()
fout.close()
