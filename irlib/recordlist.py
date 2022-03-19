""" Contains the `RecordList` class, which in addition to being useful for
functions that try to directly read XML metadata from HDF datasets, is also
used as the metadata container for `Gather` objects. """

import sys
import os
import re
import numpy as np
import traceback
import datetime
import h5py 
import pdb


def pcdateconvert(pcdatetime, datefmt='ddmm'):
    '''
    This will take a pcdatetime in BSI format - either 
    dd/mm/yy_hh:mm:ss.sss AM/PM (default) or  mm/dd/yy_.... (if requested)
    and make it into a datetime object
    
    datefmt can only be 'ddmm' (since fall 2016) or 'mmdd' <2016-Sept
    '''
    
    assert (datefmt == 'ddmm' or datefmt == 'mmdd' ), "Bad date format specified"
     
    pcdate, pctime = pcdatetime.split("_")
    if datefmt == 'ddmm':
        dd,mm,yy = pcdate.split("/")
    else:        
        mm,dd,yy = pcdate.split("/")  # It USED to be this way 
    hms, ampm = pctime.split()
    h,m,s = hms.split(":")
    if len(yy) == 2:
        yy = "20"+yy
    convdate = datetime.datetime(int(yy), int(mm), int(dd), int(h), int(m),\
                                int(float(s)), int(float(s)%1*1e6) )
    if ampm.lower() == "pm" and int(h) < 12:
        convdate += datetime.timedelta(0, 3600*12)
    if ampm.lower() == "am" and int(h) == 12:
        convdate = convdate - datetime.timedelta(0, 3600*12)
        
    return convdate

def TimeFromComment(infile, line, loc):
        '''
        This function finds the comment field from the datacapture_0 group using 
        the low-level hdf 5 api called h5g.  
        
        This is the only place where PC timestamp is located in some versions of the 
        hdf data format
        
        Note that when a dataset is opened as an h5l object, the group comment 
        is not available, hence the reason for reopening the file (maybe could be 
        done only once for efficiency but it doesn't seem to be an issue so far)
            
        infile - hdf file name
        line  - the line as str 'line_0'
        loc  - the loc as str 'loc_0'
        returns a datetime object
         
        '''
        
        h = h5py.File(infile)
        dt = h[line][loc].id.get_comment(b'datacapture_0').decode()
        return datetime.datetime.strptime(dt, "%m/%d/%Y %I:%M:%S %p")
    
    
def isodate(dt):
    ''' Formats datetime object dt to iso date
    '''
    return dt.strftime("%Y-%m-%d %H:%M:%S.%f")




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

        #now with more metadata fields
        self.attrs = ['fids', 'filenames', 'lines', 'locations',
                      'datacaptures', 'echograms', 'timestamps', 'lats',
                      'lons','gps_time', 'fix_qual', 'num_sat', 'dilution', 'alt_asl',
                      'geoid_height', 'gps_fix_valid', 'gps_message_ok',
                      'datums', 'eastings', 'northings', 'elevations', 'zones',
                      'vrange', 'sample_rate','stacking','trig_level','rec_len', 
                      'startbuf', 'buftime', 'pps','comments']

        for attr in self.attrs:
            setattr(self, attr, [])
        self.hasUTM = False
        return

    @staticmethod
    def _xmlGetValF(xml, name):
        """ Look up a float value in an XML fragment. Return NaN if not found.
        """
        # Modified search to be less restrictive
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>(-?[0-9.E-]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)        
        if m is not None:
            return float(m.group().split('<Val>')[1].split('</Val>')[0])
        else:
            return np.nan

    @staticmethod
    def _xmlGetValI(xml, name):
        """ Look up an integer value in an XML fragment. Return None if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>([0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return int(float(m.group().split('<Val>')[1].split('</Val>')[0]))
        else:
            return None

    @staticmethod
    def _xmlGetValS(xml, name):
        """ Look up a string value in an XML fragment. Return an empty string if not found.
        """
        m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>(-?[0-9.]+?)</Val>'.format(
                        name.replace(' ', '\s')), xml, flags=re.IGNORECASE)
        if m is not None:
            return m.group().split('<Val>')[1].split('</Val>')[0]
        else:
            return ''

    @staticmethod
    def _dm2dec(dmstr):
        """ Convert the degree - decimal minute codes in radar data
        to a decimal degree coordinate. dmstr is expected to a string.
        """
        if dmstr == '': return
        try:
            hem = 1
            a,b = dmstr.split(".")
            if float(a[:-2]) <0:
                hem = -1
                a = a[1:]
            return hem*round(float(a[:-2]) +
                         float(a[-2:])/60. + float("." + b)/60.,6)
        except (AttributeError, ValueError):
            return None

    def AddDataset(self, dataset, fid=None):
        """ Add metadata from a new dataset to the RecordList instance. Updates
        the RecordList internal lists with data parsed from the radar xml.

        Does not read pick data.

        Parameters
        ----------
        dataset : an h5py dataset at the `echogram` level
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
        line, loc, dc, eg = [int(s.split("_")[1]) for s in splitname[1:5]]
        self.lines.append(line)
        self.locations.append(loc)
        self.datacaptures.append(dc)
        self.echograms.append(eg)

        if fid is None:
            fid = "{0:0>4}{1:0>4}{2:0>4}{3:0>4}".format(line, loc, dc, eg)

        # Timestamps
        if 'Save timestamp' in dataset.attrs:
            # 2008
            self.timestamps.append(dataset.attrs['Save timestamp'])
        elif 'PCSavetimestamp' in dataset.attrs:
            # 2009 and later
            pcdatetime = dataset.attrs["PCSavetimestamp"]
            # there are various formats.  Decide which is which by splitting the string
            # and manipulating it. 
            
            if not type(pcdatetime) == str:
                pcdatetime = pcdatetime.astype(str) #convert to string, this converts byte-encoded data
            if len(pcdatetime.split(",")) == 4:
                timestamp, startbuf,buftime,pps = pcdatetime.split(",")
                self.timestamps.append(isodate(pcdateconvert(timestamp, datefmt='ddmm')))
                self.startbuf.append(startbuf.split(":")[1])
                self.buftime.append(buftime.split(":")[1])
                self.pps.append(pps)
            elif len(pcdatetime.split(",")) == 3:
                startbuf,buftime,pps = pcdatetime.split(",")
                self.startbuf.append(startbuf.split(":")[1])
                self.buftime.append(buftime.split(":")[1])
                self.pps.append(pps)
                #timestamp for this is in a completely different place.  
                self.timestamps.append(isodate(TimeFromComment(self.filename, 
                                                              splitname[1], splitname[2])))
            else:
                self.timestamps.append(isodate(pcdateconvert(pcdatetime, datefmt='mmdd'))) # guessing the format                
                self.startbuf.append("")
                self.buftime.append("")
                self.pps.append("")
        else:
            raise ParseError('Timestamp read failure', dataset.name)

        # XML parsing code (unused categories set to None for speed)
        # Parse main cluster
        try:

            try:   # This is the old way h5py library decodes based on data type specified
                xml = dataset.attrs['GPS Cluster- MetaData_xml'].decode("utf-8")
            except:  # This is the newer way, should work h5py >= 3.0 
                xml = dataset.attrs['GPS Cluster- MetaData_xml'] 
            
            if self._xmlGetValS(xml, 'Lat') == '' :   # old format            
                self.lats.append(self._dm2dec(self._xmlGetValS(xml, 'Lat_N'))) 
                # the Long_ W space here is important since this IS the variable name (will change in ver 6.2 IceRadar)
                self.lons.append(self._dm2dec(self._xmlGetValS(xml, 'Long_ W'))) 
            else:
                self.lats.append(self._dm2dec(self._xmlGetValS(xml, 'Lat')))  #Changed from Lat_N
                self.lons.append(self._dm2dec(self._xmlGetValS(xml, 'Long')))  #Changed from Long_W
            
            self.gps_time.append(self._xmlGetValS(xml, 'GPS_timestamp_UTC'))            

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

            try:   # This is the old way h5py library decodes based on data type specified
                xml = dataset.attrs['Digitizer-MetaData_xml'].decode("utf-8")
            except:  # This is the newer way, should work h5py >= 3.0 
                xml = dataset.attrs['Digitizer-MetaData_xml'] 
            self.vrange.append(self._xmlGetValF(xml, 'vertical range'))
            self.sample_rate.append(self._xmlGetValF(xml, ' sample rate'))
            self.stacking.append(self._xmlGetValI(xml, 'Stacking'))
            self.trig_level.append(self._xmlGetValF(xml, 'trigger level'))
            self.rec_len.append(self._xmlGetValI(xml, 'Record Length'))
            
        except:
            with open('error.log', 'w') as f:
                traceback.print_exc(file=f)
            raise ParseError('Digitizer cluster read failure', dataset.name)

        # Parse UTM cluster if available (2009 and later?)
        if 'GPS Cluster_UTM-MetaData_xml' in dataset.attrs:
            self.hasUTM = True
            try:

                try:   # This is the old way h5py library decodes based on data type specified
                    xml = dataset.attrs['GPS Cluster_UTM-MetaData_xml'].decode("utf-8")
                except:  # This is the newer way, should work h5py >= 3.0 
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
            self.comments.append(dataset.parent.id.get_comment('.'.encode("utf-8")))
        except:
            with open('error.log', 'w') as f:
                traceback.print_exc(file=f)
            raise ParseError('HDF Group comment read failure')

        return

    def Write(self, f, eastern_hemisphere=False):
        """ Write out the data stored internally in CSV format to a file
        object f. 
        
        IF lat and lon are both al
        
        """
        error = 0
        
        # If this is true, then either we are in the Eastern & Northern Hemisphere
        # Or this is the old format where the lat/lon were unsigned
        if len([lon for lon in self.lons if (lon is not None and lon < 0)]) + \
                len([lat for lat in self.lats if (lat is not None and lat < 0)]) == 0:
            if not eastern_hemisphere:
                # Invert longitudes to be in the western hemisphere
                self.lons = [-i if i is not None else i for i in self.lons]
        

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
            "gps_time," +
            "fix_qual," +
            "num_sat," +
            "dilution," +
            "alt_asl," +
            "geoid_ht," +
            "gps_fix," +
            "gps_ok," +
            "vertical_range," +
            "sample_rate," +
            "stacking," +
            "trig_level," +
            "rec_len," +
            "startbuf," +
            "buftime," +
            "pps")
        if self.hasUTM:
            header += (
                ",datum," +
                "easting," +
                "northing," +
                "elevation," +
                "zone")
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
                    str(self.gps_time[i]) + "," +
                    str(self.fix_qual[i]) + "," +
                    str(self.num_sat[i]) + "," +
                    str(self.dilution[i]) + "," +
                    str(self.alt_asl[i]) + "," +
                    str(self.geoid_height[i]) + "," +
                    str(self.gps_fix_valid[i]) + "," +
                    str(self.gps_message_ok[i]) + "," +
                    str(self.vrange[i]) + "," +
                    str(self.sample_rate[i]) + "," +
                    str(self.stacking[i]) + "," +
                    str(self.trig_level[i]) + "," +
                    str(self.rec_len[i]) + "," +
                    str(self.startbuf[i]) + "," +
                    str(self.buftime[i]) + "," +
                    str(self.pps[i])                    
                    )
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

