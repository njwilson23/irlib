""" Contains the `Survey` class, which is the overarching `irlib` structure.
`Survey` classes handle interaction with the raw HDF datasets, and spawn off
`Gather` classes for analysis. Each HDF dataset can be opened as a `Survey`,
which stores file references and collects metadata in the form of a
`FileHandler`. Radar lines can be created from a `Survey` using the
`ExtractLine` method, which returns a `Gather`. """


import h5py
import numpy as np
import os
import cPickle
import traceback

from gather import *
from filehandler import *
from recordlist import *

from autovivification import AutoVivification

class Survey:
    """ Surveys can be broken down into **Gathers** and *traces*. To create a
    survey and spawn a gather, do something like::

        # Create the survey
        S = Survey("mysurvey.h5")

        # Create the gather (Line)
        linenumber = 0      # This can be any nonzero integer
        datacapture = 0     # This corresponds to the channel frequency in
                            # dualdar setups
        L = S.ExtractLine(linenumber, dc=0)

        # To see what can be done with *line gather* **L**, refer to the
        # documentation in irlib/gather.py
    """

    def __init__(self, datafile):
        """ Instantiate a **Survey** object. A survey encompasses one HDF5 file
        generated from Blue System Inc. IceRadar software.         """
        self.datafile = datafile
        try:
            self.f = h5py.File(self.datafile, 'r')
            self.status = 'open'
            # Create 2-level boolean map of the dataset
            self.retain = AutoVivification()
            for line in self.f:
                if isinstance(self.f[line], h5py.Group):
                    for location in list(self.f[line]):
                        self.retain[line][location] = True
        except IOError:
            sys.stdout.write("No survey exists with the" +
                             " filename:\n\t{0}\n".format(datafile))
            self.close()
        return

    def __del__(self):
        try:
            self.f.close()
        except AttributeError:
            pass
        return

    def __repr__(self):
        return self.status + " survey object: " + self.datafile

    def _getdatasets(self, line=None):
        """ Return a list of datasets. Optionally, specify an integer
        line. """
        if isinstance(line, int):
            path = 'line_{0}/'.format(line)
        else:
            path = '/'
        names = []
        self.f[path].visit(names.append)
        datasets = []
        for name in names:
            if (type(f[name]) == h5py.Dataset) and \
              ('picked' not in self.f[path][name].name):
                datasets.append(path+name)
        return datasets

    def _path2fid(self, path, linloc_only = False):
        """ Based on a path, return a unique FID for database
        relations. """
        try:
            # Index from [1:] to cut out any '/' that might be present
            # Line number
            lin = int(path[1:].split('/',1)[0].split('_',1)[1])
            # Location number
            loc = int(path[1:].split('/',2)[1].split('_',1)[1])
            if not linloc_only:
                # Datacapture number
                dc = int(path[1:].split('/',3)[2].split('_',1)[1])
                # Echogram number
                eg = int(path[1:].split('/',3)[2].split('_',1)[1])
            else:
                dc = 0
                eg = 0
            fid = str(lin).rjust(4,'0') + str(loc).rjust(4,'0') \
                + str(dc).rjust(4,'0') + str(eg).rjust(4,'0')
            return fid
        except:
            traceback.print_exc()
            sys.stderr.write('survey: failed at path2fid')
            return None

    def Lines(self):
        return self.GetLines()

    def GetLines(self):
        """ Return the lines (str) contained within the survey. """
        lines = [name for name in self.f.keys() if name[:4] == 'line']
        lines.sort(key=(lambda s: int(s.split('_')[1])))
        return lines

    def GetChannelsInLine(self, lineno):
        """ Return the number of channels (datacaptures per location) in a
        line. If the number is not constant throughout the line, then return
        the maximum number. """
        try:
            line = self.GetLines()[lineno]
        except IndexError:
            sys.stderr.write("lineno out of range ({0} not in {1}:{2})\n"
                    .format(lineno, 0, len(self.GetLines)))
        dclist = [self.f[line][loc].keys() for loc in self.f[line].keys()]
        return max([len(a) for a in dclist])

    def ExtractTrace(self, line, location, datacapture=0, echogram=0):
        """ Extract the values for a trace and return as a vector. """
        path = ('line_{lin}/location_{loc}/datacapture_{dc}/'
                'echogram_{eg}'.format(lin=line, loc=location, dc=datacapture,
                                       eg=echogram))
        vec = self.f[path].value
        return vec

    def ExtractLine(self, line, bounds=(None,None), datacapture=0,
                    fromcache=False, cache_dir="cache", print_fnm=False,
                    verbose=False):
        """ Extract every trace on a line. If bounds are supplied
        (min, max), limit extraction to only the range specified.
        Return a LineGather instance.

        Parameters
        ----------
        line : line number to extract (int)
        bounds : return a specific data slice (iter x2)
        datacapture : datacapture subset to load (int or iter<ints)
        fromcache : attempt to load from a cached file (bool)
        cache_dir : specify a cache directory (str)
        print_fnm : print the cache search path (bool)
        """

        if fromcache:
            fnm = self.GetLineCacheName(line, dc=datacapture, cache_dir=cache_dir)
            if print_fnm:
                print fnm
            if os.path.isfile(fnm):
                try:
                    with open(fnm, 'r') as f:
                        unpickler = cPickle.Unpickler(f)
                        gatherdata = unpickler.load()
                except:
                    traceback.print_exc()
                    gatherdata = None
                return gatherdata
            else:
                print "Cached file not available; loading from HDF"

        path = 'line_{lin}'.format(lin=line)

        # Separate out all datasets on the correct line
        names = []
        self.f[path].visit(names.append)

        # Filter out the datasets, then the correct datacaptures
        # The following is a nested filter that first keeps elements of type
        # *h5py.Dataset*, next discards picked data, and finally restricts the
        # names to the datacaptures specified by *datacapture*.
        try:
            allowed_datacaptures = ["datacapture_{0}".format(dc) for dc in datacapture]
        except TypeError:
            allowed_datacaptures = ["datacapture_{0}".format(datacapture)]
        datasets = filter(lambda c: c.split('/')[-2] in allowed_datacaptures,
                    filter(lambda b: 'picked' not in self.f[path][b].name,
                    filter(lambda a: isinstance(self.f[path][a], h5py.Dataset),
                    names)))
        if len(datasets) == 0:
            sys.stderr.write("no datasets match the specified channel(s)\n")

        # Sort the datasets by location number
        try:
            datasets.sort(key=(lambda s: int(s.split('/')[0].split('_')[1])))
        except:
            traceback.print_exc()
            sys.stderr.write("Error sorting datasets by " + \
                        "location number in ExtractLine()\n")

        # If bounds are specified, slice out the excess locations
        try:
            if bounds[1] != None: datasets = datasets[:bounds[1]]
            if bounds[0] != None: datasets = datasets[bounds[0]:]
        except TypeError:
            sys.stderr.write("bounds kwarg in ExtractLine() " + \
                            "must be a two element list or tuple\n")

        # Grab XML metadata
        metadata = RecordList(self.datafile)
        for trace in datasets:
            full_path = path + '/' + trace
            try:
                metadata.AddDataset(self.f[path][trace], fid=self._path2fid(full_path))
            except ParseError as e:
                if verbose:
                    sys.stderr.write(e.message + '\n')
                metadata.CropRecords()

        # Pull out all the data, concatenate into a single array
        try:
            arr = np.array((self.f[path][datasets[0]].value,)).T
        except IndexError:
            sys.stderr.write("error indexing {0} - it might be empty\n".format(path))
            return

        for dataset in datasets[1:]:
            newrow = np.array((self.f[path][dataset].value,)).T
            # Test to see if the sample length is the same.
            # Resize to match the largest sample langth.
            # This will cause some funny steps in the radargram
            if newrow.shape[0] < arr.shape[0]:
                newrow = np.vstack((newrow, np.zeros((arr.shape[0]-newrow.shape[0], 1))))
            elif newrow.shape[0] > arr.shape[0]:
                arr = np.vstack((arr, np.zeros((newrow.shape[0]-arr.shape[0], arr.shape[1]))))
            arr = np.concatenate((arr, newrow), axis=1)

        return LineGather(arr, infile=self.datafile, line=line,
                metadata=metadata, retain=self.retain['line_{0}'.format(line)],
                dc=datacapture)

    def GetLineCacheName(self, line, dc=0, cache_dir="cache"):
        """ Return a standard cache name for a given *line* (int) and *dc*
        (int). """
        cnm = os.path.join(cache_dir,
                os.path.splitext(os.path.basename(self.datafile))[0] + \
                '_line' + str(line) + '_' + str(dc) + '.ird')
        return cnm

    def GetWaveletPolarity(self, line=None, aw_threshold=3., bw_threshold=5.,
        window=31, bw_early_lim=60, outfile='polarity_data.csv'):
        """ Attempt to identify wavelet polarities. Return a dictionary
        with trace IDs as keys and a tuple with polarity tag and
        estimates of certainty for the air and basal waves.

        Parameters
        ----------

        line : line to be queried; if None, will do all lines
        aw_threshold : factor of power increase required for a
                       positive identification of air wave
        bw_threshold : factor of power increase required for a
                       positive identification of reflected wave
        window : window for determining local background power;
                 must be odd

        Returns
        -------

        Polarity tag : `1` if air and basal reflection wave has same polarity,
            `-1` if air and basal reflection wave have opposite polarity, and
            `0` if uncertain

        Certainty estimate : difference between instantaneous power and
                             windowed power
        """
        makeplot = True
        aw_val = 0
        bw_val = 0

        datasets = self._getdatasets(line)
        outdict = {}
        for path in datasets:
            vec = self.f[path].value
            dc_phase = 0
            bw_polarity = 0

            # Calculate initial window rms
            square_sum = (vec[0:window]**2).sum()

            # Try to get the airwave polarity
            for val in range((window-1)/2, bw_early_lim):
                if vec[val]**2 > aw_threshold * square_sum/window:
                    # Airwave found
                    dc_phase = int(vec[val] / abs(vec[val]))
                    aw_cert = abs(vec[val] - math.sqrt(square_sum/window))
                    ### DEBUG ###
                    if makeplot:
                        aw_val = val
                    #############
                    break
                else:
                    # Go to the next val after updating square_sum
                    square_sum -= (vec[val - (window-1)/2]) ** 2
                    square_sum += (vec[val + (window+1)/2]) ** 2

            if dc_phase == 0:
                aw_cert = 0.
                bw_cert = 0.
                aw_val = 0
                # Skip looking for the reflection wave - no point

            else:
                # Calculate initial window rms
                square_sum = (vec[bw_early_lim-(window-1)/2:bw_early_lim+(window+1)/2]**2).sum()

                # Try to get the reflection polarity
                for val in range(bw_early_lim, vec.shape[0]-(window+1)/2):
                    if vec[val]**2 > bw_threshold * square_sum/window:
                        # Bed reflection wave found
                        # I make the polarity negative here, because it finds
                        # the second cycle of a Ricker wavelet
                        bw_polarity = -int(vec[val] / abs(vec[val]))
                        bw_cert = abs(vec[val] - math.sqrt(square_sum/window))
                        ### DEBUG ###
                        if makeplot:
                            bw_val = val
                            print path,'\t',bw_val,'\t',bw_polarity
                        #############
                        break
                    else:
                        # Go to the next val after updating square_sum
                        square_sum -= (vec[val - (window-1)/2]) ** 2
                        square_sum += (vec[val + (window+1)/2]) ** 2

                if bw_polarity == 0:
                    bw_cert = 0.
                    bw_val = 0

            # Compare polarities
            if dc_phase == 0:
                polbool = 0
            elif dc_phase == bw_polarity:
                polbool = 1
            else:
                polbool = -1

            if makeplot:
                plt.plot(vec, '-k')
                plt.plot([aw_val,aw_val],[-.05,.05],'-b')
                plt.plot([bw_val,bw_val],[-.05,.05],'-r')
                plt.savefig('pick_plots/'+self._path2fid(path, linloc_only=True))
                plt.clf()

            # Make a note
            outdict[self._path2fid(path, linloc_only=True)] = \
                (polbool, aw_cert, bw_cert)

        # Write outdict to a CSV
        keys = outdict.keys()
        keys.sort()
        with open(outfile, 'w') as f:
            f.write('"FID","TAG","AW_CERT","BW_CERT"\n')
            for key in keys:
                dval = outdict[key]
                f.write('\"{key}\",{tag},{aw_cert},{bw_cert}\n'.format(
                    key=key, tag=dval[0], aw_cert=dval[1], bw_cert=dval[2]))
        return


    def WriteHDF5(self, fnm, overwrite=False):
        """ Given a filename, write the contents of the original file to the
        new HDF5, whereever self.retain == True. The usage case for this is
        when bad data have been identified in the original file.

        Note that for now, this does not preserve HDF5 object comments.
        """
        if os.path.exists(fnm) and not overwrite:
            print 'already exists'
            return

        fout = h5py.File(fnm, 'w')
        print "working"

        for line in self.f:
            if isinstance(self.f[line], h5py.Group):
                try:
                    fout.create_group(line)
                except ValueError:
                    print "somehow, {0} already existed in Survey.WriteHDF5()". format(line)
                    print "this might be a problem, and you should look into it."
                print "\t{0}".format(line)
                for location in list(self.f[line]):
                    if self.retain[line][location]:
                        self.f.copy('{0}/{1}'.format(line, location), fout[line])
        fout.close()
        return

    def close(self):
        try:
            self.f.close()
        except AttributeError:
            pass
        self.status = 'closed'
        return
