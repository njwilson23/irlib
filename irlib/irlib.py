#! /usr/bin/env python
#
#   ICE-PENETRATING RADAR LIBRARY
#
#   Provides the underlying classes for ice radar utilities.
#   Depends on:
#       h5py                HDF file interface
#       numpy/scipy         Numerics and math libraries
#       matplotlib          Plotting
#       aaigrid_driver      ASCII grid class for raster IO
#
#   Works best with:
#       mig_fk              FK migration module
#       agc_cy              Cython-accelerated AGC
#       pywavelet           Wavelet transform function
#

import h5py
import numpy as np
import aaigrid_driver as aai
import scipy.signal as sig
import scipy.spatial as spatial
import matplotlib.pyplot as plt
import os, os.path, sys, math, cPickle
import copy
import traceback, pdb




def path2fid(path, linloc_only = False):
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


def LoadCoords(line, filename, L):
    """ Load the coordinates from filename. """
    # Get the paths for the traces within the line
    with h5py.File(filename, 'r') as f:
        allnames = []
        f.visit(allnames.append)
        paths = []
        for name in allnames:
            if (isinstance(f[name], h5py.Dataset)
                and name.split('/')[0] == 'line_'+str(line)):
                paths.append(name)

        # Sort the paths by location number
        paths.sort(key=(lambda s: int(s.split('/')[1].split('_')[1])))

        # Add the traces to a RecordList object
        R = RecordList(filename)
        fids = []
        for path in paths:
            fid = path2fid(path)
            fids.append(fid)
            R.AddDataset(f[path], fid=fid)

    # Extract the coordinates
    L.lats = R.lats
    L.lons = R.lons
    if R.hasUTM:
        L.eastings = R.eastings
        L.northings = R.northings

    return (R.lats, R.lons), (R.eastings, R.northings), fids


def ExtractAttrs(h5file, outfile=None, fout=None, flip_lon=True):
    """ Extract the metadata for each trace in a radar archive.
    Optionally write this data out to a comma-delimited file that can be
    imported into a GIS.
    """
    # Open the file object
    f = h5py.File(h5file, 'r')

    # Get all of the names of datasets in the file
    names = []
    f.visit(names.append)
    datasets = []
    for name in names:
        if (type(f[name]) == h5py.Dataset) and \
          ('picked' not in f[name].name):
            datasets.append(name)
    sys.stderr.write("\t{n} datasets found\n".format(n=len(names)))

    # Record the ID params, UTC time, and coordinates for each dataset
    records = RecordList(h5file)
    e = 0
    for name in datasets:
        fid = path2fid(name)
        e += records.AddDataset(f[name], fid=fid)

    sys.stderr.write("\t{f} read: {e} problems\n".format(
            f=os.path.basename(h5file), e=e))

    eout = 0
    if outfile:
        # Export as a CSV, creating a file with name outfile
        with open(outfile, 'w') as fout:
            eout = records.Write(fout, flip_lon=flip_lon)
        sys.stderr.write("\t{f} written: {e} problems\n".format(
                f=os.path.basename(output), e=eout))
    elif fout:
        # Export as a CSV, using the provided file object
        eout = records.Write(fout, flip_lon=flip_lon)
        sys.stderr.write("\tstream generated: {e} problems\n".format(e=eout))

    f.close()

    return e+eout, records


def ExtractTrace(h5file, line, location, datacapture=0, echogram=0):
    """ Extract the values for a trace and return as a vector. """
    with h5py.File(h5file, 'r') as f:
        path = 'line_{lin}/location_{loc}/datacapture_{dc}/echogram_{eg}'.format(
            lin=line, loc=location, dc=datacapture, eg=echogram)
        vec = f[path].value
    return vec


def ExtractLine(h5file, line, bounds=(None,None)):
    """ Extract the values from every trace on a line. If bounds are
    supplied (min, max), limit extraction to only the range specified.
    Return a Numpy array.
    """
    with h5py.File(h5file, 'r') as f:
        path = 'line_{lin}'.format(lin=line)

        # Separate out all datasets on the correct line
        names = []
        f[path].visit(names.append)
        datasets = []
        for name in names:
            if (type(f[path][name]) == h5py.Dataset) and \
              ('picked' not in f[path][name].name):
                datasets.append(name)

        # Sort the datasets by location number
        try:
            datasets.sort(key=(
                    lambda s: int(s.split('/')[0].split('_')[1])
                    ))
        except:
            traceback.print_exc()
            sys.stderr.write("Error sorting datasets by location " + \
                            "number in ExtractLine()\n")

        # If bounds are specified, slice out the excess locations
        try:
            if bounds[1] != None: datasets = datasets[:bounds[1]]
            if bounds[0] != None: datasets = datasets[bounds[0]:]
        except TypeError:
            sys.stderr.write("bounds kwarg in ExtractLine() must be" + \
                            " a two element list or tuple")

        # Pull them all out, concatenate into a single array
        try:
            arr = np.array((f[path][datasets[0]].value,)).T
        except IndexError:
            sys.stderr.write("error indexing {0} - it might be empty\n".format(path))
            return

        for dataset in datasets[1:]:
            newrow = np.array((f[path][dataset].value,)).T
            # Test to see if the sample length is the same.
            # Resize to match the largest sample langth.
            # This will cause some funny steps in the radargram -
            # hopefully these will disappear when setting zero-time
            if newrow.shape[0] < arr.shape[0]:
                newrow = np.vstack((newrow, np.zeros((arr.shape[0]-newrow.shape[0], 1))))
            elif newrow.shape[0] > arr.shape[0]:
                arr = np.vstack((arr, np.zeros((newrow.shape[0]-arr.shape[0], arr.shape[1]))))
            arr = np.concatenate((arr, newrow), axis=1)

    return arr


#def ApplyFilter(GatherInstance, cmd):
    ### DEPRECATED IN FAVOUR OF filter_defs.py ###
    #if isinstance(cmd, list):
        #args = [i for i in cmd[1:]]
        #cmd = cmd[0]
    #else:
        #args = None

    #try:
        #if cmd == 'lowpass':
            #GatherInstance.DoWindowedSinc(cutoff=25.e6, bandwidth=15.e6, mode='lowpass')
        #elif cmd == 'highpass':
            #GatherInstance.DoWindowedSinc(cutoff=5.e6, bandwidth=10.e6, mode='highpass')
        #elif cmd == 'lowpass_ma':
            #GatherInstance.DoMoveAvg(21, kind='blackman', mode='lowpass')
        #elif cmd == 'highpass_ma':
            #GatherInstance.DoMoveAvg(31, kind='blackman', mode='highpass')
        #elif cmd == 'dewow':
            #GatherInstance.Dewow()
        #elif cmd == 'bed10':
            #GatherInstance.DoWindowedSinc(cutoff=8.e6, bandwidth=5.e6, mode='highpass')
            #GatherInstance.DoWindowedSinc(cutoff=25.e6, bandwidth=10.e6, mode='lowpass')
        #elif cmd == 'bed35':
            #GatherInstance.DoWindowedSinc(cutoff=15.e6, bandwidth=10.e6, mode='highpass')
            #GatherInstance.DoWindowedSinc(cutoff=45.e6, bandwidth=10.e6, mode='lowpass')
        #elif cmd == 'bed50':
            #GatherInstance.DoWindowedSinc(cutoff=35.e6, bandwidth=10.e6, mode='highpass')
            #GatherInstance.DoWindowedSinc(cutoff=45.e6, bandwidth=10.e6, mode='lowpass')
        #elif cmd == 'bed':
            #GatherInstance.DoWindowedSinc(cutoff=20.e6, bandwidth=10.e6, mode='highpass')
            #GatherInstance.DoWindowedSinc(cutoff=30.e6, bandwidth=10.e6, mode='lowpass')
        #elif cmd == 'mult':
            #if len(args) >= 2:
                #m = float(args[1])
            #else:
                #m = 2.0
            #GatherInstance.MultiplyAmplitude(m)
        #elif cmd == 'gc':
            #GatherInstance.DoTimeGainControl(ncoef=2.e6, npow=1.)
            #GatherInstance.Dewow()
        #elif cmd == 'gc2':
            #GatherInstance.DoTimeGainControl(ncoef=2.e6, npow=2.)
            #GatherInstance.Dewow()
        #elif cmd == 'agc':
            #GatherInstance.DoAutoGainControl(25e-8)
            #GatherInstance.DoWindowedSinc(cutoff=5.e6, bandwidth=10.e6, mode='highpass')
        #elif cmd == 'fkmig':
            #GatherInstance.MigrateFK()
        #elif cmd == 'kirmig':
            #print "not yet implemented"
        #elif cmd == 'eng10':
            ## Englacial scatter enhancer that works well for 10 MHz data
            #GatherInstance.Reset()
            #GatherInstance.Dewow()
            #GatherInstance.DoTimeGainControl(ncoef=2.e6, npow=1.)
            #GatherInstance.DoWindowedSinc(40.e6, bandwidth=4.e6, mode='highpass')
            #D_masked = np.ma.masked_where(GatherInstance.data==0, np.abs(GatherInstance.data))
            #D_transformed = np.abs(D_masked**0.33) * np.sign(D_masked)
            #D_transformed.mask = np.ma.nomask
            #GatherInstance.data = D_transformed
            #GatherInstance.history.append(('cubic_transformation'))
            #GatherInstance.history.append(('absolute_value'))
            #GatherInstance.DoMoveAvg(7, kind='blackman', mode='lowpass')
        #elif cmd == 'eng35':
            ## Englacial scatter enhancer that works better for 35 and 50 MHz
            #GatherInstance.DoWindowedSinc(cutoff=40e6, bandwidth=5e6, mode='highpass')
            #GatherInstance.DoWindowedSinc(cutoff=100e6, bandwidth=20e6, mode='lowpass')
            #GatherInstance.data = np.abs(GatherInstance.data)
        #elif cmd == 'engd':
            ## Difference based englacial scatter filter. Not particularly impressive.
            #A = GatherInstance.data.copy()
            #GatherInstance.DoWindowedSinc(cutoff=40e6, bandwidth=5e6, mode='highpass')
            #GatherInstance.DoWindowedSinc(cutoff=100e6, bandwidth=20e6, mode='lowpass')
            #B = GatherInstance.data.copy()
            #GatherInstance.data = A-B
            #GatherInstance.Dewow()
            #r0 = A.max() - A.min()
            #r = GatherInstance.data.max() - GatherInstance.data.min()
            #GatherInstance.data = GatherInstance.data * (r0 / r)
        #else:
            #print "Filter type not recognized"
    #except:
        #traceback.print_exc()
    #return


def PlotTrace(D, Dp=None, Dpp=None, Dn=None, Dnn=None, outfile=None, rate=1e-8, spacing=0.1, title=None):
    """ Plot a trace (vector) and its neighbours if provided. Vectors
    are expected to be numpy arrays. If 'outfile' is not None, then
    save image rather than showing.
    """
    n = D.size
    T = np.arange(0, n*rate, rate)      # time axis

    plt.figure(figsize=(8,4))
    plt.axes(axisbg='black')
    plt.grid(b=True, color='white')
    plt.xlabel('Time (ns)')
    plt.ylabel('Voltage (V)')

    plt.plot(T, D, '-w')

    if Dp is not None:
        plt.plot(T, Dp-spacing, '-y')

        if Dpp is not None:
            plt.plot(T, Dpp-2*spacing, '-y')

    if Dn is not None:
        plt.plot(T, Dn+spacing, '-y')

        if Dnn is not None:
            plt.plot(T, Dnn+2*spacing, '-y')

    locs, labels = plt.xticks()
    plt.xticks(locs, locs*10/rate)

    if title is not None:
        plt.title(title)

    if outfile == None:
        plt.show()
    else:
        plt.savefig(outfile)

    return


def PlotLine(arr, outfile=None, cmap='gray', title='Radar Line',
            rate=1e-8, c=1.68e8):
    """ Plot an array containing radar data along a line. """

    n = arr.size
    T = np.arange(0, n*rate, rate)      # time axis

    # Find the luminescence range so that plot intensity is symmetric
    lum_bound = max((abs(arr.max()), abs(arr.min())))

    # Draw it
    img = plt.imshow(arr, aspect='auto', cmap=cmap, vmin=-lum_bound, vmax=lum_bound)
    plt.title(title)
    plt.xlabel("Location Number")
    plt.ylabel("Time (ns)")
    locs, labels = plt.yticks()
    plt.yticks(locs, locs*1.e9*rate)

    if outfile == None:
        plt.show()
    else:
        plt.savefig(outfile)

    return img


def TryCache(fnm):
    """ Tests whether a cached dataset with the given filename exists. If it
    does, then it unpickles it, and returns True, along with a reference to the
    dataset. It if does not, it returns False and a None value.
    """
    if os.path.isfile(fnm):
        try:
            with open(fnm, 'r') as f:
                unpickler = cPickle.Unpickler(f)
                dataset = unpickler.load()
        except:
            traceback.print_exc()
            dataset = None
        return True, dataset
    else:
        return False, None

