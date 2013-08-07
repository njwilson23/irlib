""" Contains the `RecordList` class, which in addition to being useful for
functions that try to directly read XML metadata from HDF datasets, is also
used as the metadata container for `Gather` objects. """

import sys
import os
import re
import numpy as np
import traceback

class RecordList:
    """ Class to simplify the extraction of metadata from HDF5 radar
    datasets.

    Usage:

    - initialize a RecordList instance with a filename (arbitrary, but should
      be the HDF filename)

    - add datasets by passing h5 dataset objects to `self.AddDataset()`
    """
    def __init__(self, filename=None):
        self.filename = filename

        self.attrs = ['fids', 'filenames', 'lines', 'locations',
                      'datacaptures', 'echograms', 'timestamps', 'lats',
                      'lons', 'fix_qual', 'num_sat', 'dilution', 'alt_asl',
                      'geoid_height', 'gps_fix_valid', 'gps_message_ok',
                      'datums', 'eastings', 'northings', 'elevations', 'zones',
                      'vrange', 'sample_rate', 'comments']

        for attr in self.attrs:
            setattr(self, attr, [])
        self.hasUTM = False
        return

    def _xmlGetValF(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>(-?[0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return float(m.group().split('<Val>')[1].split('</Val>')[0])
        else:
            return np.nan

    def _xmlGetValI(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>([0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return int(float(m.group().split('<Val>')[1].split('</Val>')[0]))
        else:
            return None

    def _xmlGetValS(self, xml, name):
        """ Look up a value in an XML fragment. Return None if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>([0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return m.group().split('<Val>')[1].split('</Val>')[0]
        else:
            return ''

    def _dm2dec(self, dmstr):
        """ Convert the degree - decimal minute codes in radar data
        to a decimal degree coordinate. dmstr is expected to a string.
        """
        if dmstr == '': return
        try:
            a,b = dmstr.split(".")
            return round(float(a[:-2]) +
                         float(a[-2:])/60. + float("." + b)/60.,6)
        except AttributeError:
            return None
        except ValueError:
            return None
        #except AttributeError:
        #    # *dmstr* is not a string or string-like
        #    return

    def AddDataset(self, dataset, fid):
        """ Add metadata from a new dataset to the RecordList instance. Updates
        the RecordList internal lists with data parsed from the radar xml.

        Does not read pick data.

        Parameters
        ----------
        dataset : a h5py dataset at the `echogram` level
                  (fh5[line][location][datacapture][echogram])
        fid : pre-defined FID for the dataset

        Returns None
        """
        # Is this really a good way? Seems inelegant... -njw
        if 'picked' in dataset.name:
            sys.stderr.write('RecordList: did not attempt to parse {0}\n'
                            .format(dataset.name))
            return

        self.filenames.append(self.filename)
        self.fids.append(fid)


        # Parse dataset name
        splitname = dataset.name.split('/')
        self.lines.append(int(splitname[1].split('_')[1]))
        self.locations.append(int(splitname[2].split('_')[1]))
        self.datacaptures.append(int(splitname[3].split('_')[1]))
        self.echograms.append(int(splitname[4].split('_')[1]))


        # Timestamps
        if 'Save timestamp' in dataset.attrs:
            # 2008
            self.timestamps.append(dataset.attrs['Save timestamp'])
        elif 'PCSavetimestamp' in dataset.attrs:
            # 2009 and later
            self.timestamps.append(dataset.attrs['PCSavetimestamp'])
        else:
            raise ParseError('Timestamp read failure', dataset.name)



        # XML parsing code (unused categories set to None for speed)

        # Parse main cluster
        try:
            xml = dataset.attrs['GPS Cluster- MetaData_xml']
            self.lats.append(self._dm2dec(self._xmlGetValS(xml, 'Lat_N')))
            self.lons.append(self._dm2dec(self._xmlGetValS(xml, 'Long_ W')))
            self.fix_qual.append(self._xmlGetValI(xml, 'Fix_Quality'))
            self.num_sat.append(self._xmlGetValI(xml, 'Num _Sat'))
            self.dilution.append(self._xmlGetValF(xml, 'Dilution'))
            self.alt_asl.append(self._xmlGetValF(xml, 'Alt_asl_m'))
            self.geoid_height.append(self._xmlGetValF(xml, 'Geoid_Heigh_m'))
            self.gps_fix_valid.append(self._xmlGetValI(xml, 'GPS Fix valid'))
            self.gps_message_ok.append(self._xmlGetValI(xml, 'GPS Message ok'))
        except:
            with open('error.log', 'w') as f:
                traceback.print_exc(file=f)
            raise ParseError('GPS cluster read failure', dataset.name)

        # Parse digitizer cluster
        try:
            xml = dataset.attrs['Digitizer-MetaData_xml']
            self.vrange.append(self._xmlGetValF(xml, 'vertical range'))
            self.sample_rate.append(self._xmlGetValF(xml, ' sample rate'))
        except:
            with open('error.log', 'w') as f:
                traceback.print_exc(file=f)
            raise ParseError('Digitizer cluster read failure', dataset.name)

        # Parse UTM cluster if available (2009 and later?)
        if 'GPS Cluster_UTM-MetaData_xml' in dataset.attrs:
            self.hasUTM = True
            try:
                xml = dataset.attrs['GPS Cluster_UTM-MetaData_xml']
                self.datums.append(self._xmlGetValS(xml, 'Datum'))
                self.eastings.append(self._xmlGetValF(xml, 'Easting_m'))
                self.northings.append(self._xmlGetValF(xml, 'Northing_m'))
                self.elevations.append(self._xmlGetValF(xml, 'Elevation'))
                self.zones.append(self._xmlGetValI(xml, 'Zone'))
            except:
                with open('error.log', 'w') as f:
                    traceback.print_exc(file=f)
                raise ParseError('Digitizer cluster read failure', dataset.name)

        # Parse comment
        try:
            self.comments.append(dataset.parent.id.get_comment('.'))
        except:
            with open('error.log', 'w') as f:
                traceback.print_exc(file=f)
            raise ParseError('HDF Group comment read failure')

        return

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

    def CropRecords(self):
        """ Ensure that all records are the same length. This should be called
        if adding a dataset fails, potentially leaving dangling records. """
        nrecs = min([len(getattr(self, attr)) for attr in self.attrs])
        for attr in self.attrs:
            data = getattr(self, attr)
            while len(data) > nrecs:
                data.pop(-1)
        return

    def Reverse(self):
        """ Reverse data in place. """
        for attr in self.attrs:
            data = getattr(self, attr)
            data.reverse()
        return

    def Cut(self, start, end):
        """ Drop section out of all attribute lists in place. """
        for attr in self.attrs:
            data = getattr(self, attr)
            del data[start:end]
        return

class ParseError(Exception):
    def __init__(self, message='', fnm=''):
        self.message = message + ": {0}".format(fnm)
    def __str__(self):
        return self.message

