"""
 Contains the RecordList() class, which in addition to being useful for
 functions that try to directly read XML metadata from HDF datasets, is also
 used as the metadata container for Gather() objects.
"""

import sys
import os
import re
import traceback, pdb

class RecordList:
    """ Class to simplify the extraction of metadata from HDF5 radar
    datasets.

    Usage:
        initialize a RecordList instance with a filename (arbitrary,
        but should be the HDF filename)

        add datasets by passing h5 dataset objects to self.AddDataset()

        cut out parts of lines with self.Cut()

        export to CSV with self.Write()
    """
    def __init__(self, filename=None):
        self.filename = filename
        self.fids = []
        self.filenames = []
        self.lines = []
        self.locations = []
        self.datacaptures = []
        self.echograms = []
        self.timestamps = []
        self.lats = []
        self.lons = []
        self.fix_qual = []
        self.num_sat = []
        self.dilution = []
        self.alt_asl = []
        self.geoid_height = []
        self.gps_fix_valid = []
        self.gps_message_ok = []
        self.datums = []
        self.eastings = []
        self.northings = []
        self.elevations = []
        self.zones = []
        self.vrange = []
        self.sample_rate = []

        self.hasUTM = False

    def _xmlGetValF(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>([0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return float(m.group().split('<Val>')[1].split('</Val>')[0])
        else:
            return

    def _xmlGetValI(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>([0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return int(float(m.group().split('<Val>')[1].split('</Val>')[0]))
        else:
            return

    def _xmlGetValS(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>([0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return m.group().split('<Val>')[1].split('</Val>')[0]
        else:
            return

    def _xmlGetValF0(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        xmlsplit = xml.splitlines()
        for i, line in enumerate(xmlsplit):
            if line == '<Name>'+name+'</Name>':
                valstr = xmlsplit[i+1]
                try:
                    return float(valstr.split('<Val>')[1].split('</Val>')[0])
                except ValueError:
                    return
        return

    def _xmlGetValI0(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        xmlsplit = xml.splitlines()
        for i, line in enumerate(xmlsplit):
            if line == '<Name>'+name+'</Name>':
                valstr = xmlsplit[i+1]
                try:
                    return int(valstr.split('<Val>')[1].split('</Val>')[0])
                except ValueError:
                    return
        return

    def _xmlGetValS0(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        xmlsplit = xml.splitlines()
        for i, line in enumerate(xmlsplit):
            if line == '<Name>'+name+'</Name>':
                valstr = xmlsplit[i+1]
                try:
                    return valstr.split('<Val>')[1].split('</Val>')[0]
                except ValueError:
                    return
        return

    def _dm2dec(self, dmstr):
        """ Convert the degree - decimal minute codes in radar data
        to a decimal degree coordinate. dmstr is expected to a string.
        """
        if dmstr == '': return
        try:
            a,b = dmstr.split(".")
            return round(float(a[:-2]) +
                         float(a[-2:])/60. + float("." + b)/60.,6)
        except ValueError:
            return
        except AttributeError:
            # *dmstr* is not a string or string-like
            return

    def AddDataset(self, dataset, fid='9'*16):
        """ Add metadata from a new dataset to the RecordList instance. Updates
        the RecordList internal lists with data parsed from the radar xml.

                dataset is an actual h5py dataset (fh5[path])

        Does not read pick data, and returns None if attempted.
        """
        if 'picked' in dataset.name:
            sys.stderr.write('RecordList: did not attempt to parse {0}\n'
                            .format(dataset.name))
            return 1

        error = 0
        self.filenames.append(self.filename)
        self.fids.append(fid)

        try:
            splitname = dataset.name.split('/')
            self.lines.append(int(splitname[1].split('_')[1]))
            self.locations.append(int(splitname[2].split('_')[1]))
            self.datacaptures.append(int(splitname[3].split('_')[1]))
            self.echograms.append(int(splitname[4].split('_')[1]))
        except:
            traceback.print_exc()
            sys.stderr.write("\tirlib: Could not properly parse dataset name\n")
            sys.stderr.write('\t' + dataset.name + '\n')
            error += 1
            # Make sure that all list lengths are still equal
            n_items = min(len(self.lines), len(self.locations),
                        len(self.datacaptures), len(self.echograms))
            for lst in [self.lines, self.locations, self.datacaptures, self.echograms]:
                while len(lst) > n_items:
                    lst.pop()

        try:
            # Different attribute labels have been used over time
            if 'Save timestamp' in dataset.attrs:
                # 2008
                self.timestamps.append(dataset.attrs['Save timestamp'])
            elif 'PCSavetimestamp' in dataset.attrs:
                # 2009 and later
                self.timestamps.append(dataset.attrs['PCSavetimestamp'])
            else:
                error += 1
                sys.stderr.write("\t{d}: timestamp could not be found\n".format(d=dataset.name))
        except:
            traceback.print_exc()
            sys.stderr.write("\tirlib: Could not save timestamp field\n")
            sys.stderr.write('\t' + dataset.name + '\n')
            self.timestamps.append(None)
            error += 1

        # XML parsing code (unused categories set to None for speed)

        # Parse main cluster
        try:
            xml = dataset.attrs['GPS Cluster- MetaData_xml']
            self.lats.append(self._dm2dec(self._xmlGetValS(xml, 'Lat_N')))
            self.lons.append(self._dm2dec(self._xmlGetValS(xml, 'Long_ W')))
            self.fix_qual.append(self._xmlGetValI(xml, 'Fix_Quality'))
            self.num_sat.append(self._xmlGetValI(xml, 'Num _Sat'))
            #self.dilution.append(self._xmlGetValF(xml, 'Dilution'))
            self.dilution.append(None)
            #self.alt_asl.append(self._xmlGetValF(xml, 'Alt_asl_m'))
            self.alt_asl.append(None)       # GPS alt_asl is very inaccurate
            #self.geoid_height.append(self._xmlGetValF(xml, 'Geoid_Heigh_m'))
            self.geoid_height.append(None)
            self.gps_fix_valid.append(self._xmlGetValI(xml, 'GPS Fix valid'))
            self.gps_message_ok.append(self._xmlGetValI(xml, 'GPS Message ok'))
        except:
            traceback.print_exc()
            sys.stderr.write("\tirlib: Could not save GPS metadata\n")
            sys.stderr.write('\t' + dataset.name + '\n')
            sys.stderr.write(xml)
            sys.exit(1)
            error += 1
            return error

        # Parse digitizer cluster
        try:
            xml = dataset.attrs['Digitizer-MetaData_xml']
            self.vrange.append(self._xmlGetValF(xml, 'vertical range'))
            self.sample_rate.append(self._xmlGetValF(xml, ' sample rate'))
        except:
            sys.stderr.write("\tirlib: Could not save digitizer metadata\n")
            sys.stderr.write('\t' + dataset.name + '\n')
            error += 1
            return error


        # Parse UTM cluster if available (2009 and later?)
        if 'GPS Cluster_UTM-MetaData_xml' in dataset.attrs:
            self.hasUTM = True
            try:
                xml = dataset.attrs['GPS Cluster_UTM-MetaData_xml']
                self.datums.append(self._xmlGetValS(xml, 'Datum'))
                self.eastings.append(self._xmlGetValF(xml, 'Easting_m'))
                self.northings.append(self._xmlGetValF(xml, 'Northing_m'))
                #self.elevations.append(self._xmlGetValF(xml, 'Elevation'))
                self.elevations.append(None)
                self.zones.append(self._xmlGetValI(xml, 'Zone'))
            except:
                sys.stderr.write("\tCould not save GPS UTM metadata\n")
                sys.stderr.write('\t' + dataset.name + '\n')
                error += 1
                return error
        return error

    def Write(self, f, flip_lon=True):
        """ Write out the data stored internally in CSV format to a file
        object f. """
        error = 0

        if flip_lon:
            # Invert longitudes to be in the right hemisphere
            # Be careful not to accidentally do this more than once!
            self.lons = [-i if i!=None else None for i in self.lons]

        header = (
            "FID," +
            "filename," +
            "line," +
            "location," +
            "datacapture," +
            "echogram," +
            "timestamp," +
            "lat," +
            "lon," +
            "fix_qual," +
            "num_sat," +
            "dilution," +
            "alt_asl," +
            "geoid_ht," +
            "gps_fix," +
            "gps_ok," +
            "vertical_range," +
            "sample_rate")
        if self.hasUTM:
            header += (
                ",datums," +
                "eastings," +
                "northings," +
                "elevations," +
                "zones")
        header += "\n"
        f.write(header)

        for i in range(len(self.filenames)):
            try:
                sout = (
                    "\"" + self.fids[i] + "\"" + "," +
                    os.path.basename(self.filenames[i]) + "," +
                    str(self.lines[i]) + "," +
                    str(self.locations[i]) + "," +
                    str(self.datacaptures[i]) + "," +
                    str(self.echograms[i]) + "," +
                    "\"" + self.timestamps[i] + "\"" + "," +
                    str(self.lats[i]) + "," +
                    str(self.lons[i]) + "," +
                    str(self.fix_qual[i]) + "," +
                    str(self.num_sat[i]) + "," +
                    str(self.dilution[i]) + "," +
                    str(self.alt_asl[i]) + "," +
                    str(self.geoid_height[i]) + "," +
                    str(self.gps_fix_valid[i]) + "," +
                    str(self.gps_message_ok[i]) + "," +
                    str(self.vrange[i]) + "," +
                    str(self.sample_rate[i]))
                if self.hasUTM:
                    sout += (
                        "," + str(self.datums[i]) + "," +
                        str(self.eastings[i]) + "," +
                        str(self.northings[i]) + "," +
                        str(self.elevations[i]) + "," +
                        str(self.zones[i]))
                sout += "\n"
                f.write(sout)
            except:
                traceback.print_exc()
                sys.stderr.write("\tError writing record to file ({0})\n".format(i))
                error += 1
        return error

    def Reverse(self):
        """ Reverse data in place. """
        self.fids.reverse()
        self.filenames.reverse()
        self.lines.reverse()
        self.locations.reverse()
        self.datacaptures.reverse()
        self.echograms.reverse()
        self.timestamps.reverse()
        self.lats.reverse()
        self.lons.reverse()
        self.fix_qual.reverse()
        self.num_sat.reverse()
        self.dilution.reverse()
        self.alt_asl.reverse()
        self.geoid_height.reverse()
        self.gps_fix_valid.reverse()
        self.gps_message_ok.reverse()
        self.vrange.reverse()
        self.sample_rate.reverse()
        self.datums.reverse()
        self.eastings.reverse()
        self.northings.reverse()
        self.elevations.reverse()
        return

    def Cut(self, start, end):
        """ Drop section out of all attribute lists in place. """
        del self.fids[start:end]
        del self.filenames[start:end]
        del self.lines[start:end]
        del self.locations[start:end]
        del self.datacaptures[start:end]
        del self.echograms[start:end]
        del self.timestamps[start:end]
        del self.lats[start:end]
        del self.lons[start:end]
        del self.fix_qual[start:end]
        del self.num_sat[start:end]
        del self.dilution[start:end]
        del self.alt_asl[start:end]
        del self.geoid_height[start:end]
        del self.gps_fix_valid[start:end]
        del self.gps_message_ok[start:end]
        del self.vrange[start:end]
        del self.sample_rate[start:end]
        del self.datums[start:end]
        del self.eastings[start:end]
        del self.northings[start:end]
        del self.elevations[start:end]
        return
