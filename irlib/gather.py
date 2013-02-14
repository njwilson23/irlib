"""
 Defines various kinds of Gather() classes, which are perhaps the most
 important classes in the irlib library. These contain the data from
 individual radar lines, and make it easy to apply various processing steps to
 the data. The Gather() class is a parent class for the CommonOffsetGather()
 and CommonMidpointGather() subclasses. The LineGather() object is deprecated
 and is now just an alias for the CommonOffsetGather(), kept for backwards
 compatibility.
"""

import scipy.signal as sig
import scipy.spatial as spatial
import numpy as np
import matplotlib.pyplot as plt
import math, copy, cPickle
import os, sys
import traceback, pdb

import aaigrid_driver as aai
from filehandler import *
from autovivification import AutoVivification

try:
    import agc_cy
except ImportError:
    sys.stderr.write("Warning: falling back to non-accellerated AGC\n")

try:
    from mig_fk import fkmig
except ImportError:
    sys.stderr.write("Warning: F-K migration module not loaded\n")

try:
    import pywavelet
except ImportError:
    sys.stderr.write("Warning: pywavelet module not loaded\n")


class Gather:
    """ Base class for common offset and common midpoint lines. """

    def __init__(self, arr, infile=None, line=None, metadata=None, retain=None,
        dc=0):
        """
        arr         data held within a numpy array
        infile      the original HDF5 dataset
        line        the line enumeration
        metadata    RecordList instance
        """
        self.raw_data = arr.copy()
        self.data = self.raw_data.copy()
        self.line = line
        self.datacapture = dc
        self.infile = infile
        self.metadata = metadata            # reference to a RecordList
        if self.metadata is not None:
            self.rate = 1./self.metadata.sample_rate[0]
            self.metadata_copy = copy.deepcopy(self.metadata)
        else:
            self.rate = None
            self.metadata_copy = None

        self.ny = self.data.shape[0]
        self.nx = self.data.shape[1]

        if retain is None:
            self.retain = AutoVivification()
            for col in range(self.nx):
                self.retain['location_{0}'.format(col)] = True
        else:
            self.retain = retain

        # Initialize FID list and pick data
        gen_fid = lambda i: (str(self.line).rjust(4,'0')
                            + str(i).rjust(4,'0') + 8*'0')
        self.fids = [gen_fid(i) for i in range(self.nx)]
        self.fids_copy = copy.copy(self.fids)
        self.bed_picks = 999 * np.ones(self.nx)
        self.bed_phase = 999 * np.ones(self.nx)
        self.dc_picks = 999 * np.ones(self.nx)
        self.dc_phase = 999 * np.ones(self.nx)

        self.history = [('init')]

    def __repr__(self):
        return ("Gather instance for {fnm} at line {line}, "
                "channel{chan}".format(fnm=self.infile, line=self.line,
                                       chan=self.datacapture))

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
            sys.stderr.write('LineGather: failed at path2fid\n')
            return None

    def _getkernel(self, width, kind):
        """ Generate a lowpass convolution kernel. """
        if width % 2 == 0:
            sys.stderr.write('convolution kernel width must be odd\n')
            return

        if kind == 'boxcar':
            kernel = 1./width * np.ones((width,))

        elif kind == 'gaussian':
            sigma = .25
            n = np.arange(0, width, 1)
            # Compute the gaussian kernel, normalizing by the integral
            kernel = (math.e ** (-.5 * ((n-(width-1)/2.) \
                                  / (sigma*(width-1)/2.))**2))
            kernel /= kernel.sum()

        elif kind == 'blackman':
            def integrate_blackman(w, a0, a1, a2):
                # Compute the integral of a blackman window
                val = a0*w \
                    - a1*(w-1)/(2*math.pi)*math.sin(2*math.pi*w/(w-2)) \
                    + a2*(w-1)/(2*math.pi)*math.sin(4*math.pi*w/(w-1))
                return val
            alpha = 0.16
            a0 = (1-alpha)/2.
            a1 = 1./2.
            a2 = alpha/2.
            n = np.arange(0, width, 1)
            # Compute the blackman kernel, normalizing by the integral
            kernel = (a0 - a1*np.cos(2*np.pi*n/(width-1)) \
                        + a2*np.cos(4*np.pi*n/(width-1)))
            kernel /= kernel.sum()

        else:
            sys.stderr.write("filter kind not recognized: {0}\n".format(
                            kind))

        return kernel

    def _lowpassma(self, width, kind='boxcar'):
        """ Internal low-pass convolution filter implementation.
        Width is the filter window in discrete samples, and must be odd.
        """
        sys.stderr.write("Method '_lowpassma' deprecated; use 'DoMoveAvg' instead\n")

        # Get the kernel in lowpass form
        kernel = self._getkernel(width, kind=kind)

        # Perform convolution on data
        for i in range(self.data.shape[1]):
            _trace = np.convolve(kernel, self.data[:,i], mode='same')
            _trace = (_trace / np.sqrt(np.mean(_trace**2))
                    * np.sqrt(np.mean(self.data[:,i]**2)))
            self.data[:,i] = _trace.copy()

        self.history.append(('convolution_A', width, kind, 'lowpass'))

    def _highpassma(self, width, kind='boxcar'):
        """ Internal high-pass convolution filter implementation.
        Uses spectral inversion to generate filter. Width is the filter
        window in discrete samples, and must be odd.
        """
        sys.stderr.write("Method '_highpassma' deprecated; use 'DoMoveAvg' instead\n")

        # Get the kernel in lowpass form
        kernel = self._getkernel(width, kind=kind)

        # Spectral inversion to generate high-pass filter
        kernel = -kernel
        kernel[(width-1)/2] += 1.

        # Perform convolution on data
        for i in range(self.data.shape[1]):
            _trace = np.convolve(kernel, self.data[:,i], mode='same')
            _trace = (_trace / np.sqrt(np.mean(_trace**2))
                    * np.sqrt(np.mean(self.data[:,i]**2)))
            self.data[:,i] = _trace.copy()

        self.history.append(('convolution_A', width, kind, 'highpass'))

    def _svd(self):
        """ Perform SVD on trace data. Results are:
            U: horizontal space domain eigenvectors
            S: singular values
            V: vertical time domain eigenvectors
        """
        U, S, V = np.linalg.svd(self.data)
        return U, S, V

    def _svd_reconstruct(self, U, S, V, i):
        """ Return an eigenimage given the SVD products. """
        return S[i] * np.dot(np.atleast_2d(U[:,i]).T, np.atleast_2d(V[i,:]))


    def GetFID(self, loc):
        """ Given a location (array index), return the location's
        FID. This solves the problem of knowing the FID of an array
        that has been sliced.
        """
        return self.fids[loc]

    def FindFID(self, fid):
        """ Given a FID, find the index within the data array. Takes
        either a single FID, or a list. """
        if hasattr(fid, "__iter__"):
            return [self.fids.index(i) for i in fid]
        else:
            return self.fids.index(fid)

    def GetCacheName(self, cache_dir="cache"):
        """ Returns a logical name for a pickled cache of self. """
        cnm = os.path.join(cache_dir,
                os.path.splitext(os.path.basename(self.infile))[0] + \
                '_line' + str(self.line) + '_' + str(self.datacapture) + \
                '.ird')
        return cnm

    def GetDigitizerFilename(self):
        """ Returns the likely name of a file containing digitized line data. """
        fnm = 'englacial/' + os.path.basename(self.infile).split('.')[0] + \
            "_line" + str(self.line) + ".txt"
        return fnm

    def LoadTopography(self, topofnm, smooth=True):
        """ Load topography along line gather, reading from an ASC file.
        Obviously, this requires the Gather to have a valid metadata attribute.

        If smooth=True, then apply a boxcar filter to soften the effects of the
        DEM discretization.
        """
        G = aai.Grid()
        if os.path.isfile(topofnm):
            G.FromFile(topofnm)
            try:
                self.topography = np.array(map(lambda x,y: G.Sample(x,y)[0],
                        self.metadata.eastings, self.metadata.northings))
            except:
                pdb.set_trace()
            self.history.append(('load_topo', topofnm))
            if smooth:
                self.SmoothenTopography()
            self.topography_copy = self.topography.copy()
        else:
            print "{0} is not a file".format(topofnm)
        return

    def SmoothenTopography(self):
        """ Load topography along line gather, reading from an ASC file.
        Obviously, this requires the Gather to have a valid metadata attribute.

        If smooth=True, then apply a boxcar filter to soften the effects of the
        DEM discretization.
        """
        def interpolate_nans(V):
            Vreal = np.nonzero(np.isnan(V) == False)[0]
            Vnan = np.nonzero(np.isnan(V))[0]
            return np.interp(range(len(V)), Vreal, V[Vreal])

        try:
            # Fill in NaNs
            self.topography = interpolate_nans(self.topography)

            # Pad and then convolve with boxcar
            self.topography = np.convolve(
                np.hstack([np.ones(5)*self.topography[0], self.topography,
                np.ones(5)*self.topography[-1]]),
                np.ones(11) / 11.0, mode='valid')
            self.history.append(('smooth_topo'))
        except AttributeError:
            print "topography does not exist"
        return

    def SmoothenGPS(self, win=5):
        """ Run a moving average over the GPS northings and eastings. """
        def interpolate_nans(V):
            Vreal = np.nonzero(np.isnan(V) == False)[0]
            Vnan = np.nonzero(np.isnan(V))[0]
            return np.interp(range(len(V)), Vreal, V[Vreal])

        try:
            eastings0 = np.array(self.metadata.eastings)
            northings0 = np.array(self.metadata.northings)
            # iInterpolate over NaNs, unless the entire array is NaN
            n_valid_eastings = len(np.nonzero(np.isnan(eastings0)==False)[0])
            n_valid_northings = len(np.nonzero(np.isnan(northings0)==False)[0])
            if (n_valid_eastings > 0) and (n_valid_northings > 0):
                eastings0 = interpolate_nans(eastings0)
                northings0 = interpolate_nans(northings0)
            # Pad and then convolve with boxcar
            eastings1 = np.convolve(
                np.hstack([np.ones(win//2)*eastings0[0], eastings0,
                np.ones(win//2)*eastings0[-1]]),
                np.ones(win) / float(win), mode='valid')
            northings1 = np.convolve(
                np.hstack([np.ones(win//2)*northings0[0], northings0,
                np.ones(win//2)*northings0[-1]]),
                np.ones(win) / float(win), mode='valid')
            self.metadata.eastings = eastings1.tolist()
            self.metadata.northings = northings1.tolist()
            self.history.append(('smooth_gps'))
        except AttributeError:
            print "UTM coordinates do not exist"
        return

    def DoMoveAvg(self, width, kind='blackman', mode='lowpass'):
        """ Time domain convolution filter implementation over A scan.
        Width is the filter window in discrete samples, and must be odd.
        """
        # Get the kernel in lowpass form
        kernel = self._getkernel(width, kind=kind)

        if mode == 'highpass':
            # Spectral inversion to generate high-pass filter
            kernel = -kernel
            kernel[(width-1)/2] += 1.

        # Perform convolution on data
        for i in range(self.data.shape[1]):
            _trace = np.convolve(kernel, self.data[:,i], mode='same')
            _trace = (_trace / np.sqrt(np.mean(_trace**2))
                    * np.sqrt(np.mean(self.data[:,i]**2)))
            self.data[:,i] = _trace.copy()

        self.data[np.isnan(self.data)] = 0.0
        self.history.append(('convolution_A', width, kind, mode))
        return

    def DoMoveAvgB(self, width, kind='blackman', mode='lowpass'):
        """ Time domain convolution filter implementation over B scan.
        Width is the filter window in discrete samples, and must be odd.
        """
        # Get the kernel in lowpass form
        kernel = self._getkernel(width, kind=kind)

        if mode == 'highpass':
            # Spectral inversion to generate high-pass filter
            kernel = -kernel
            kernel[(width-1)/2] += 1.

        # Perform convolution on data
        for i in range(self.data.shape[0]):
            _row = np.convolve(kernel, self.data[i,:], mode='same')
            _row = (_row / np.sqrt(np.mean(_row**2))
                    * np.sqrt(np.mean(self.data[i,:]**2)))
            self.data[i,:] = _row.copy()

        self.data[np.isnan(self.data)] = 0.0
        self.history.append(('convolution_B', width, kind, mode))
        return

    def Dewow(self, cutoff=4e6):
        """ Apply a dewow (highpass) filter to remove very low
        frequency signals. This is a step up from using a simple
        demean operation.

            cutoff: cutoff minimum frequency, in Hz

        This is just an interface to _highpassma.
        """
        # Compute filter width as a function of the cutoff frequency
        width = 41      # SEEMS TO WORK
        self.DoMoveAvg(width, kind='blackman', mode='highpass')
        return

    def MultiplyAmplitude(self, multiplier):
        """ Just makes viewing easier. """
        self.data *= multiplier
        self.history.append(("MultiplyAmplitude", multiplier))
        return

    def DoTimeGainControl(self, ncoef=1., npow=1., nexp=0., gamma=1., bias=0.):
        """ Apply a gain enhancement as a function of time.

            Transform F: f(t) -> F(f(t))

            F: (f(t) * t^npow * exp(nexp*t)) ^ gamma + bias

            Note: ncoef functionality has been removed and average
            power is preserved instead.

            Claerbout (1985) suggests:
                npow = 2.
                gamma = 0.5
        """
        if self.rate is not None:
            dt = self.rate
        else:
            dt = 1e-8

        t = dt * np.arange(1, self.data.shape[0]+1)

        for i in range(self.data.shape[1]):
            if self.data[:,i].mean() != 0:
                A = (self.data[:,i] * (t)**npow * np.exp(nexp*t))
                B = A/abs(A) * abs(A) ** gamma + bias
                B[np.isnan(B)] = 0.

                # Scale to preserve average amplitude
                self.data[:,i] = (B / np.sqrt(np.mean(B**2))
                                * np.sqrt(np.mean(self.data[:,i]**2)))

        self.history.append(('time_gain_control',ncoef,npow,nexp,gamma,bias,dt))
        return

    def DoMurrayGainControl(self, npow=2., tswitch=100.):
        """ Apply a time-power gain enhancement up to a time limit,
        after which gain is constant (e.g. Murray et al, 1997).

            tswitch is in samples
        """
        t = np.arange(1, self.data.shape[0]+1)
        tf = min((tswitch, self.data.shape[0]))
        for i in range(self.data.shape[1]):
            self.data[0:tf,i] *= t[0:tf]**npow
        if self.data.shape[0] > tf:
            for i in range(self.data.shape[1]):
                self.data[tf:,i] *= tswitch**npow
        self.history.append(('murray_gain_control',npow,tswitch))
        return

    def DoAutoGainControl(self, timewin=20e-8):
        """ Apply the RMS-based AGC algorithm from Seismic Unix*.

            timewin: time half-window in seconds
            dt: sample interval in seconds

        Try to use a fast Cython-accellerated version. If that fails, fall
        back to the Numpy version.

        *Cohen, J.K. and Stockwell, J.W. CWP/SU: Seismic Unix Release
        42. Colorado School of Mines, Center for Wave Phenomena. (1996)
        """

        if self.rate is not None:
            dt = self.rate
        else:
            dt = 1e-8

        try:
            agc_cy.DoAutoGainControl(self, timewin)

        except NameError:
            agcdata = np.zeros(self.data.shape)
            iwagc = int(round(timewin/dt))     # Half-window in samples
            nwin = iwagc*2 + 1

            for i in range(self.data.shape[1]):

                if sum(self.data[:,i]**2) != 0:
                    trace = self.data[:,i]
                    agctrace = np.zeros(trace.shape)
                    nt = len(trace)
                    blocksum = np.sum(trace[:iwagc]**2)

                    # Compute the first half-window
                    for j in range(iwagc):
                        blocksum = blocksum + trace[j+iwagc]**2
                        agctrace[j] = trace[j] / math.sqrt(blocksum / (j+iwagc))

                    # Compute the middle, full-rms window
                    for j in range(iwagc, nt-iwagc):
                        blocksum = blocksum - trace[j-iwagc]**2
                        blocksum = blocksum + trace[j+iwagc]**2
                        agctrace[j] = trace[j] / math.sqrt(blocksum / (2*iwagc))

                    # Compute the last half-window
                    for j in range(nt-iwagc, nt):
                        blocksum = blocksum - trace[j-iwagc]**2
                        try:
                            agctrace[j] = trace[j] / math.sqrt(blocksum / (nt-j+iwagc))
                        except ValueError:
                            # Sometimes there's a domain error here - should look into it
                            traceback.print_exc()

                    # Scale to preserve average amplitude
                    agctrace[np.isnan(agctrace)] = 0.
                    agctrace = (agctrace / np.sqrt(np.mean(agctrace**2))
                                * np.sqrt(np.mean(trace**2)))

                    # Copy array into new line
                    self.data[:,i] = agctrace.copy()

        self.history.append(('auto_gain_control',timewin,dt))
        return

    def DoWindowedSinc(self, cutoff, bandwidth=2.e6, mode='lowpass'):
        """ Implement a windowed sinc frequency-domain filter. This has
        better performance characteristics than a Chebyschev filter, at
        the expense of execution speed.

            cutoff: cutoff frequency in Hz
            bandwidth: transition bandwidth in Hz
            mode: 'lowpass' or 'highpass'
        """
        if self.rate is not None:
            r = 1.0 / self.rate
        else:
            r = 1e8

        fc = cutoff / r
        bw = bandwidth / r
        if max((fc, bw)) > 0.5:
            raise LineGatherError("cutoff or bandwidth exceed Nyquist "
                                  "frequency")
        M = int(4 // bw)
        if M % 2 == 1:
            M += 1

        X = np.arange(M+1) - M / 2
        kernel = sig.blackman(M+1) * np.sinc(2*fc*X)
        kernel /= np.sum(kernel)
        if mode == 'highpass':
            kernel = -kernel
            kernel[M/2] += 1.0

        for i in range(self.data.shape[1]):
            self.data[:,i] = np.convolve(kernel, self.data[:,i], mode='same')

        self.history.append(('windowed_sinc', cutoff, bandwidth, mode))
        return

    def DoWindowedSinc_old(self, cutoff, bandwidth=2.e6, mode='lowpass'):
        """ Implement a windowed sinc frequency-domain filter. This has
        better performance characteristics than a Chebyschev filter, at
        the expense of execution speed.

            cutoff: cutoff frequency in Hz
            bandwidth: transition bandwidth in Hz
        """
        if self.rate is not None:
            rate = self.rate
        else:
            rate = 1e-8

        bw = bandwidth * rate
        fc = cutoff * rate
        if bw > 0.5 or fc > 0.5:
            print ("cutoff and bandwidth must be below the Nyquist "
                   "frequency ({0:e}".format(0.5/rate))
            return
        N = int(round(4./bw))
        if N % 2 == 0:
            N += 1
        kernel = sig.firwin(N, 2*fc, window='blackman')
        if mode == 'highpass':
            kernel = -kernel
            kernel[(N-1)/2] += 1.0
        # Convolve every trace
        for i in range(self.data.shape[1]):
            self.data[:,i] = np.convolve(kernel, self.data[:,i], mode='same')
        # Record
        self.history.append(('windowed_sinc', cutoff, bandwidth, mode))
        return

    def DoRecursiveFilter(self, wp, ws, gp=0.001, gs=20.0, ftype='cheby1'):
        """ Perform filtering with an IIR-type recursive filter. Shallow
        wrapper around scipy.iirdesign. Performs zero-phase filtering,
        otherwise picking times might not be accurate."""
        try:
            b, a = sig.iirdesign(2.0*wp*self.rate, 2.0*ws*self.rate,
                                gpass=gp, gstop=gs, ftype=ftype, output='ba')
            for i in range(self.data.shape[1]):
                self.data[:,i] = sig.filtfilt(b, a, self.data[:,i])
            self.history.append(('iir_filter', ftype, (wp, ws, gp, gs)))
        except OverflowError:
            sys.stderr.write('failed to create IIR filter - bad parameters\n')
        return

    def DoWienerFilter(self, window=5, noise=None):
        """ Noise removal using a Wiener statistical filter. *window* must be
        an odd positive integer. """
        assert window > 0
        assert window % 2 == 1
        for i in range(self.data.shape[1]):
            self.data[:,i] = sig.wiener(self.data[:,i], mysize=window,
                                        noise=noise)
        self.history.append(('wiener_filter', window, noise))
        return

    def ConstructEigenimage(self, i):
        """ Return the i-th eigenimage of data. """
        assert i < min((self.nx, self.ny))
        U, S, V = self._svd()
        return self._svd_reconstruct(U, S, V, i)

    def RetainEigenimageRange(self, rng):
        """ Replace data with a slice of its eigenimages. Takes a slice object
        `rng` as an argument. """
        U, S, V = self._svd()
        E = np.dstack([self._svd_reconstruct(U, S, V, i) for i in range(len(S))])
        self.data = np.sum(E[:,:,rng], axis=2)
        self.history.append(('eigenimage_range', rng.start, rng.stop, rng.step))
        return

    def RemoveRinging(self, n=2):
        """ Attempt to filter out ringing clutter by substracting the first `n`
        eigenimages. """
        self.RetainEigenimageRange(slice(5, None))
        return

    def RemoveHorizontal(self):
        """ Remove horizontal features (clutter) by de-meaning each row in the
        data array. For example, this will attenuate the air wave. """
        self.data -= np.tile(self.data.mean(axis=1), [self.data.shape[1], 1]).T
        return

    def WaveletTransform(self, trno, m=0):
        """ Perform a wavelet transform of a single trace.
        Arguments:
            trno:       trace number to transform (< self.data.shape[1])
            mother:     mother wavelet (0==Morlet (default), 1==Paul, 2==DOG)
        Returns
            wva:        complex 2-D array, where amplitude is np.abs(wva)
            scales:     scales used
            period:     fourier periods of the scales used
            coi:        e-folding factor used for cone-of-influence
        """
        y = self.data[:,trno].copy()
        dt = self.rate
        mother = 0              # Mother wavelet (0==Morlet, 1==Paul, 2==DOG)
        param = -1              # Use default wavelet parameter
        if mother == 0:         # Smallest wavelet scale; heuristic based on type
            s0 = dt
        elif mother == 1:
            s0 = dt/4
        else:
            s0 = 2*dt
        dj = 0.05               # Spacing between discrete scales (<= 0.25)
        jtot = 125               # Total number of scales
        npad = 2**np.ceil(np.log(len(y)) / np.log(2)).astype(int)
                                # Total number of points including padding
        try:
            wva, scales, periods, coi \
                = pywavelet.pywavelet(y, dt, mother, param, s0, dj, jtot, npad)
        except NameError:
            sys.stderr.write("pywavelet module required for wavelet transform\n")
        return wva, scales, periods, coi


    def PickBed(self, sbracket=(60,200), bounds=(None,None), phase=1):
        """ Attempt to pick bed reflections along the line. Return pick
        data and estimated polarity in vectors (also stored internally).

            sbracket: tuple defining minimum and maximum times (by
                sample number) during which the event can be picked

            bounds: location number limits for picking
        """
        if self.rate is not None:
            rate = self.rate
        else:
            rate = 4e-9

        if bounds[0] == None:
            istart = 0
        else:
            istart = int(bounds[0])

        if bounds[1] == None:
            iend = self.data.shape[1]
        else:
            iend = int(bounds[1]) + 1

        def first_break_bed(A, prewin=30):
            """ Find the most negative point, and then find the first positive
            dV/dt before, within a buffer **prewin**. May require dewow for
            good performance. """
            try:
                # Most negative V
                ineg = (A == A[prewin:].min()).nonzero()[0][0]
                # Most positive dV/dt prior
                dA = np.diff(A[ineg-prewin:ineg])
                ipos1 = np.argmax(dA[10:]) + ineg - prewin + 10
                # Closest point with V ~ 0 prior
                prewin = min((ipos1, prewin))
                abs_win = np.abs(A[ipos1-prewin:ipos1])
                try:
                    izero = (abs_win < 1e-3).nonzero()[0][-1] + ipos1 - prewin
                except IndexError:
                    izero = np.argmin(abs_win) + ipos1 - prewin
            except:
                # Except what?
                izero = np.nan
            return izero

        # Apply function over all traces
        try:
            picks = map(first_break_bed, self.data[sbracket[0]:sbracket[1],istart:iend].T)
            picks = np.where(np.isnan(picks), 999, picks)
            self.bed_picks[istart:iend] = picks + sbracket[0]
            self.bed_phase[istart:iend] = 1
        except:
            traceback.print_exc()

        self.history.append(('pick_bed',sbracket,bounds))
        return self.bed_picks, self.bed_phase

    def PickDC(self, sbracket=(20,50), bounds=(None,None)):
        """ Attempt to pick direct coupling waves along the line.
        Return pick data in a vector (also stored internally).

            sbracket: tuple defining minimum and maximum times (by
                sample number) during which the event can be picked

            bounds: location number limits for picking
        """
        if bounds[0] == None:
            istart = 0
        else:
            istart = int(bounds[0])

        if bounds[1] == None:
            iend = self.data.shape[1]
        else:
            iend = int(bounds[1]) + 1

        def first_break_dc(A, prewin=10, dmax=1e-3):
            """ Find the first break of a voltage spike that is 50% of
            the highest in the window. The first break is where dV/dt
            is less than dmax.
            """
            try:
                # High energy indices
                he_bool = abs(A) > 0.5 * abs(A).max()
                ihalfpwr = he_bool.nonzero()[0][0]
                # Indices of sufficiently low power change
                le_bool = np.diff(A[ihalfpwr-prewin:ihalfpwr]) < dmax
                return le_bool.nonzero()[0][-1] + ihalfpwr - prewin
            except:
                return 999

        # Apply function over all traces
        try:
            picks = map(first_break_dc, self.data[sbracket[0]:sbracket[1],istart:iend].T)
            self.dc_picks[istart:iend] = picks
            self.dc_picks[istart:iend] += sbracket[0]
            self.dc_phase[istart:iend] = 1
        except:
            traceback.print_exc()

        self.history.append(('pick_dc',sbracket,bounds))
        return self.dc_picks, self.dc_phase

    def SavePicks(self, outfile, picks, mode='bed'):
        """ Save picks to a text file. """
        FH = FileHandler(outfile, self.line, fids=self.metadata.fids)
        if mode is 'bed':
            FH.AddBedPicks(picks)
        elif mode is 'dc':
            FH.AddDCPicks(picks)
        FH.ComputeTravelTimes()
        FH.Write()
        del FH
        return

    def _loadpicks(self, infile):
        """ Load picks from a text file.
        Future: this should probably be done with a FileHandler() instead.
        """
        with open(infile, 'r') as f:
            loadvec = 999 * np.ones(self.bed_picks.shape)
            _ = f.readline()
            for i in range(len(loadvec)):
                fline = f.readline()
                loadvec[i] = int(float(fline.split(',')[1]))
        return loadvec

    def LoadPicks(self, infile):
        """ Load bed picks from a text file. Employs a FileHandler. """
        try:
            F = FileHandler(infile, self.line)
            dc_points, bed_points = F.GetEventVals()
            self.dc_picks = np.array(dc_points)
            self.bed_picks = np.array(bed_points)
        except:
            traceback.print_exc()
            sys.stderr.write("Load failed\n")
        return

    def LoadLineFeatures(self, infile):
        """ Load digitized line features, such as those generated by irview's
        dexport command. Returns a dictionary with the feature data.
        """
        features = {}
        i = 0
        fidlist = self.fids.tolist()

        with open(infile, 'r') as f:
            while True:
                # Read a feature
                pnt_list = []

                while True:
                    # Read a point
                    s = f.readline()
                    if s in ('\n', ''):
                        break
                    else:
                        try:
                            slist = s.split()
                            fid = slist[0]
                            # x needs to refer to an actual trace, so find the
                            # one with the matching FID
                            try:
                                x = float(fidlist.index(fid))
                            except TypeError:
                                x = float(fidlist.index(fid + 4*'0'))
                            # y is just the depths in samples
                            y = float(slist[3])
                            pnt_list.append((x, y))
                        except:
                            traceback.print_exc()
                            sys.stdout.write("Failed to read record:\n\t{0}".format(s))

                if len(pnt_list) == 0:
                    break
                else:
                    features[i] = ['', pnt_list]
                    i += 1

        return features

    def RemoveBlankTraces(self, nsmp=100, threshold=4e-5):
        """ Attempt to identify and remove traces that did not trigger properly
        based on the cumulative energy in the first nsmp (int) samples. """
        demean = lambda A: A - A.mean()
        sum_energy = lambda A: (A**2).sum()
        kill_list = []
        for i, trace in enumerate(self.data.T):
            if sum_energy(demean(trace)) < threshold:
                kill_list.append(i)
        self.RemoveTraces(kill_list)

        self.metadata_copy = copy.deepcopy(self.metadata)
        self.raw_data = self.data.copy()
        self.fids_copy = copy.copy(self.fids)
        self.history.append(('remove_blank', nsmp, threshold))
        return

    def Reset(self):
        """ Reset data and metadata using the internal *_copy attribute
        variables. Does not undo the effects of operations that overwrite
        these attribute (such as most preprocessing routines).
        """
        self.data = self.raw_data.copy()
        self.nx = self.data.shape[1]
        if self.metadata_copy is not None:
            self.metadata = copy.deepcopy(self.metadata_copy)
        try:
            self.topography = self.topography_copy.copy()
        except AttributeError:
            pass
        self.fids = copy.copy(self.fids_copy)
        self.bed_picks = 999 * np.ones(self.data.shape[1])
        self.bed_phase = 999 * np.ones(self.data.shape[1])
        self.dc_picks = 999 * np.ones(self.data.shape[1])
        self.dc_phase = 999 * np.ones(self.data.shape[1])
        self.history = [('init')]
        return

    def RemoveTraces(self, kill_list):
        """ Remove the traces indicated by indices in kill_list (iterable).
        Remove metadata as well. This is intended to be easier to use than
        the CutRegion() and CutSingle() methods. """
        kill_list.sort()
        kill_list.reverse()
        keep_list = list(set(range(
                len(self.metadata.locations))).difference(set(kill_list)))

        if len(keep_list) > 0:

            try:
                self.fids = [self.fids[i] for i in keep_list]
                self.data = np.vstack([self.data[:,i] for i in keep_list]).T
                self.bed_picks = np.hstack([self.bed_picks[i] for i in keep_list])
                self.bed_phase = np.hstack([self.bed_phase[i] for i in keep_list])
                self.dc_picks = np.hstack([self.dc_picks[i] for i in keep_list])
                self.dc_phase = np.hstack([self.dc_phase[i] for i in keep_list])
            except IndexError:
                print "Inconsistent data lengths in {0}".format(repr(self))
                traceback.print_exc()

            if hasattr(self, 'topography'):
                self.topography = np.hstack([self.topography[i] for i in keep_list])

            for i in kill_list:
                self.metadata.Cut(i, i+1)
                self.retain['location_{0}'.format(i)] = False

            self.nx = self.data.shape[1]

        else:
            self.fids = []
            self.data = np.array([])
            self.bed_picks = np.array([])
            self.bed_phase = np.array([])
            self.dc_picks = np.array([])
            self.dc_phase = np.array([])

            if hasattr(self, 'topography'):
                self.topography = np.array([])

            for i in kill_list:
                self.metadata.Cut(i, i+1)
                self.retain['location_{0}'.format(i)] = False

            self.nx = 0

        return

    def RemoveMetadata(self, kill_list, update_registers=True):
        """ Remove the metadata corresponding to traces indicated by indices
        in kill_list (iterable). This is intended to be easier to use than
        the CutRegion() and CutSingle() methods. """
        kill_list.sort()
        kill_list.reverse()
        all_locs = range(len(self.metadata.locations))
        keep_list = list(set(all_locs).difference(set(kill_list)))

        try:
            self.fids = [self.fids[i] for i in keep_list]
            self.bed_picks = np.hstack([self.bed_picks[i] for i in keep_list])
            self.bed_phase = np.hstack([self.bed_phase[i] for i in keep_list])
            self.dc_picks = np.hstack([self.dc_picks[i] for i in keep_list])
            self.dc_phase = np.hstack([self.dc_phase[i] for i in keep_list])
        except IndexError:
            print "Inconsistent data lengths in {0}".format(repr(self))
            traceback.print_exc()

        if hasattr(self, 'topography'):
            self.topography = np.hstack([self.topography[i] for i in keep_list])

        for i in kill_list:
            self.metadata.Cut(i, i+1)
            self.retain['location_{0}'.format(i)] = False
        return

    def CutRegion(self, start=None, end=None, cut_data=False):
        """ Cleanly cut a region specified start, end from working data and
        internal pick registers. If either start or end is None, then cut
        to the limit of the array. If both are None, raise an exception.

        Does not adjust metadata; this should be done manually with the
        RecordList.Cut() method.

        Try RemoveTraces(). It's better.
        """
        if (start is None) and (end is None):
            raise LineGatherError('must specify limits in Cut()')
        elif start is None:
            if cut_data:
                self.data = self.data[:,end+1:]
            self.fids = self.fids[end+1:]
            self.bed_picks = self.bed_picks[end+1:]
            self.bed_phase = self.bed_phase[end+1:]
            self.dc_picks = self.dc_picks[end+1:]
            self.dc_phase = self.dc_phase[end+1:]
        elif end is None:
            if cut_data:
                self.data = self.data[:,:start]
            self.fids = self.fids[:start]
            self.bed_picks = self.bed_picks[:start]
            self.bed_phase = self.bed_phase[:start]
            self.dc_picks = self.dc_picks[:start]
            self.dc_phase = self.dc_phase[:start]
        else:
            if cut_data:
                self.data = np.hstack(
                    [self.data[:,:start], self.data[:,end+1:]])
            del self.fids[start:end+1]
            self.bed_picks = np.hstack(
                    [self.bed_picks[:start], self.bed_picks[end+1:]])
            self.bed_phase = np.hstack(
                    [self.bed_phase[:start], self.bed_phase[end+1:]])
            self.dc_picks = np.hstack(
                    [self.dc_picks[:start], self.dc_phase[end+1:]])
            self.dc_phase = np.hstack(
                    [self.dc_phase[:start], self.dc_phase[end+1:]])
        return

    def CutSingle(self, loc, cut_data=False):
        """ Cleanly cut a single location from working data and
        internal pick registers. If cut_data is True, then cut the data.
        Otherwise, just cut the pick registers.

        Does not adjust metadata; this should be done manually with the
        RecordList.Cut() method.

        Try RemoveTraces(). It's better.
        """
        if loc != len(self.fids):
            if cut_data:
                self.data = np.hstack(
                    [self.data[:,:loc], self.data[:,loc+1:]])
            del self.fids[loc]
            self.bed_picks = np.hstack(
                    [self.bed_picks[:loc], self.bed_picks[loc+1:]])
            self.bed_phase = np.hstack(
                    [self.bed_phase[:loc], self.bed_phase[loc+1:]])
            self.dc_picks = np.hstack(
                    [self.dc_picks[:loc], self.dc_phase[loc+1:]])
            self.dc_phase = np.hstack(
                    [self.dc_phase[:loc], self.dc_phase[loc+1:]])
        else:
            if cut_data:
                self.data = self.data[:,:loc]
            del self.fids[loc]
            self.bed_picks = self.bed_picks[:loc]
            self.bed_phase = self.bed_phase[:loc]
            self.dc_picks = self.dc_phase[:loc]
            self.dc_phase = self.dc_phase[:loc]
        return

    def Dump(self, fnm=None):
        """ Dumps self into a cache with the given filename using
        cPickle. Returns boolean on exit indicating success or
        failure.
        """
        if fnm is None:
            fnm = self.GetCacheName()
        try:
            with open(fnm, 'w') as f:
                pickler = cPickle.Pickler(f, cPickle.HIGHEST_PROTOCOL)
                pickler.dump(self)
            return True
        except:
            traceback.print_exc()
            return False


class CommonOffsetGather(Gather):
    """ Subclass defining common-offset specific data and
    operations. """

    def Reverse(self):
        """ Flip gather data. """
        self.data = self.data[:,::-1]
        self.fids.reverse()
        self.metadata.Reverse()
        self.bed_picks = self.bed_picks[::-1]
        self.bed_phase = self.bed_picks[::-1]
        self.dc_picks = self.bed_picks[::-1]
        self.dc_phase = self.bed_picks[::-1]
        try:
            self.topography = self.topography[::-1]
        except AttributeError:
            # no topography
            pass
        return

    def GetTopoCorrectedData(self):
        """ Stop-gap method for getting topographically corrected
        radargrams. Eventually a more graceful way should be built
        into the API, but for now, this returns a new array without
        messing with self.data (the repercussions of which I haven't
        fully considered yet).
        """
        try:
            assert self.topography is not None
        except AttributeError:
            raise LineGatherError("no topography")
            return

        # Pad data so that there's room for the whole vertical range
        topo_range = max(self.topography) - min(self.topography)
        max_topo = int(math.ceil(max(self.topography)))
        time_range = topo_range / 1.68e8 * 2
        if self.rate is not None:
            smp_range = int(math.ceil(time_range / self.rate))
        else:
            # Assume Pico digitizer
            smp_range = int(math.ceil(time_range / 4e-9))
        ntr = self.data.shape[1]
        padded_data = np.vstack([self.data, np.nan * np.ones([smp_range, ntr])])
        data = np.zeros_like(padded_data)

        # Shift all of the data by an amount according to topography
        for itr, trace in enumerate(padded_data.T):
            shift = int(round((max_topo - self.topography[itr]) / 1.68e8 / self.rate)) * 2
            data[:,itr] = np.roll(trace, shift)

        return data

    def FindLineBreaks(self, threshold=0.35):
        """ Find vertices along a scattered line. """
        try:
            x = np.array(self.metadata.eastings)
            y = np.array(self.metadata.northings)

            # Find the unit direction vector between every point
            dx = np.diff(x)
            dy = np.diff(y)
            dxn = dx / (np.sqrt(dx**2 + dy**2))
            dyn = dy / (np.sqrt(dx**2 + dy**2))

            # Smooth over the vectors
            kernel = np.ones([20]) * 0.05
            dxnf = np.convolve(kernel, dxn, mode='same')
            dynf = np.convolve(kernel, dyn, mode='same')

            breaks = [0]
            # Start at first 20, and check running dot product magnitude from there
            i=20
            while i < len(dxnf)-1:
                u = np.array([dxnf[breaks[-1]:i].mean(),
                              dynf[breaks[-1]:i].mean()])
                v = np.array([dxnf[i+1], dynf[i+1]])
                # This needs to normalize the vectors, stupid!
                if abs(np.dot(u,v)) < threshold:
                    breaks.append(i)
                i += 1
            breaks.append(i+1)
        except:
            traceback.print_exc()
            breaks = None
        return breaks

    def ProjNearestNeighbour(self, dx):
        """ Determine an equally-spaced best fit straight line, and
        place the nearest trace onto each node.
        Spacing in meters is given by dx.
        DEPRECATED - SEE self.LineProject_Nearest()
        """
        Xmesh, Ymesh, proj_arr, sum_dist = self.LineProject_Nearest(dp=dx)
        return proj_arr

    def LineProjectXY(self, bounds=None, eastings=None, northings=None, sane=True):
        """ Project coordinates onto a best-fit line. Return the
        projected eastings and northings and the line fit information.

        If sane is False, sanity checks on coordinates will be skipped
        """
        # Fit a line, using metadata if no coordinates are provided
        if (eastings is None) or (northings is None):
            eastings = np.array(self.metadata.eastings)
            northings = np.array(self.metadata.northings)
        else:
            eastings = np.array(eastings)
            northings = np.array(northings)
            if len(eastings) != len(northings):
                raise LineGatherError("eastings and northings have unequal length")
                return
        if bounds is not None:
            eastings = eastings[bounds[0]:bounds[1]]
            northings = northings[bounds[0]:bounds[1]]

        if sane:
            # Do sanity checks; these impose requirements on allowable
            # metadata that are not in practise expected to be too restrictive
            if False in (eastings > -7e6):    # No straddling more than 1 UTM zone
                raise LineGatherError("unallowable eastings in metadata")
            if False in (northings > 0):    # No equatorial straddling
                raise LineGatherError("northings must be non-negative")
            if False in (eastings < 7e6):  # No straddling more than 1 UTM zone
                raise LineGatherError("unallowable eastings in metadata")
            if False in (northings < 11e7): # No extra-polar coordinates
                raise LineGatherError("northings must not cross geographical poles")

        # If the vertical range is much larger than the horizontal
        # range, polyfit can give bad results, so swap them here
        extremes = lambda a: max(a) - min(a)
        if extremes(northings) / extremes(eastings) > 1.0:
            reversed_xy = True
            _eastings = eastings.copy()
            eastings = northings.copy()
            northings = _eastings
        else:
            reversed_xy = False

        p = np.polyfit(eastings, northings, 1)

        # Project data vectors onto best-fit line vector
        m = p[0]
        c = p[1]
        X = (m*northings + eastings - c*m) / (m**2 + 1.0)   # Projection formula
        Y = (m*X + c)

        if reversed_xy:
            _X = X.copy()
            X = Y.copy()
            Y = _X

        return X, Y, p

    def LineProject_Nearest(self, bounds=None, eastings=None, northings=None, dp=4.0):
        """ Populate a best-fit line with traces nearest to equally-spaced
        points. Return the selected traces in a 2D array, as well as the
        coordinates of the projection.

        argument dp is point spacing in meters
        """
        # Get the best-fit line p and the projections of all traces onto p
        X, Y, p = self.LineProjectXY(bounds=bounds, eastings=eastings,
                                    northings=northings)
        distances = np.sqrt((X-X[0])**2 + (Y-Y[0])**2)

        # Step along line from here, creating regular points
        dx = math.sqrt(dp**2 / (p[0]**2+1))
        dy = math.sqrt(dp**2 - dx**2)
        Xmesh = np.arange(min(X), max(X)+dx, dx)
        if X[0] > X[-1]:
            Xmesh = Xmesh[::-1]
        Ymesh = p[0] * Xmesh + p[1]

        # Find the nearest trace to every point on the line
        self.proj_arr = np.zeros((self.data.shape[0], len(Xmesh)))
        kdtree = spatial.KDTree(zip(X, Y))
        sum_dist = 0.
        for i,pt in enumerate(zip(Xmesh, Ymesh)):
            dist, ii = kdtree.query((pt[0], pt[1]))
            self.proj_arr[:,i] = self.data[:,ii]
            sum_dist += dist
        return Xmesh, Ymesh, self.proj_arr, sum_dist

    def FixStaticGPS(self):
        """ Attempt to approximate locations with data from when the GPS was
        operating on static mode.
        """
        # Find regions where displacement is 0 for a while, and then large
        eastings = np.array(self.metadata.eastings)
        northings = np.array(self.metadata.northings)
        dx = eastings[1:] - eastings[:-1]
        dy = northings[1:] - northings[:-1]
        displacement = np.sqrt(dx**2 + dy**2)
        izreg = np.nonzero(displacement==0.0)[0]    # Indices of zero displacement
        d_izreg = izreg[1:] - izreg[:-1]            # Differences in index
        first_izreg = izreg[2:][(d_izreg[:-1]>1) * (d_izreg[1:]==1)]
                    # This sets beginnings of static regions to be where d_izreg
                    # is 1 but the previous value was larger
        static_regions = []
        for i in first_izreg:
            try:
                iend = np.nonzero(displacement[i:] > 0.0)[0][0]
                if displacement[i+iend] > 10.0:
                    # Create a list of indices delineating static regions
                    static_regions.append((i, i+iend+1))
            except IndexError:
                # Data ended during the last region (no way to fix this)
                # Best to discard this data
                self.CutRegion(start=i, cut_data=True)
                self.metadata.Cut(i, self.nx)
                self.metadata.Cut(i, self.nx)
                for iloc in range(i, self.nx):
                    self.retain['location_{0}'.format(iloc)] = False

        # Interpolate linearly, fixing the metadata as we go
        if len(static_regions) > 0:
            for region in static_regions:
                x0 = self.metadata.eastings[region[0]]
                xf = self.metadata.eastings[region[1]]
                y0 = self.metadata.northings[region[0]]
                yf = self.metadata.northings[region[1]]
                di = region[1] - region[0] - 1
                dx = (xf-x0) / float(di)
                dy = (yf-y0) / float(di)
                for i in range(1, di+1):
                    new_x = x0 + dx*float(i)
                    new_y = y0 + dy*float(i)
                    self.metadata.eastings[i+region[0]] = new_x
                    self.metadata.northings[i+region[0]] = new_y

        self.raw_data = self.data.copy()        # I hope I don't regret this
        self.metadata_copy = copy.deepcopy(self.metadata)
        self.fids_copy = copy.copy(self.fids)
        self.nx = self.data.shape[1]
        self.history.append(('fix_static_gps'))
        return

    def RemoveBadLocations(self, bbox=None):
        """ Remove traces where the location is None or outside of an
        optional bounding box.

            bbox    [east, west, south, north] bounding box (optional)
        """
        kill_list = []
        for i, (x, y) in enumerate(zip(self.metadata.eastings,
                                       self.metadata.northings)):
            if None in (x,y):
                kill_list.append(i)
            elif bbox is not None and \
                ((x<bbox[0]) or (x>bbox[1]) or (y<bbox[2]) or (y>bbox[3])):
                kill_list.append(i)
        self.RemoveTraces(kill_list)

        self.metadata_copy = copy.deepcopy(self.metadata)
        self.raw_data = self.data.copy()
        self.fids_copy = copy.copy(self.fids)
        self.history.append(('remove_bad_locs', bbox))
        return

    def RemoveStationary(self, threshold=3.0, debug=False):
        """ Remove consecutive points with very similar GPS locations. Do this
        by finding points within a minimum distance of each other, and averaging
        them.
        Beware incorrect GPS readings, which need to be rectified first.

            threshold is the minimum distance to be considered to be moving
        """
        if self.metadata.hasUTM is False:
            raise LineGatherError('RemoveStationary: no UTM coordinates available')
            return
        eastings = np.array(self.metadata.eastings)
        northings = np.array(self.metadata.northings)

        dbg_traces_deleted = 0

        # Get rid of traces with no location data
        kill_list = list(np.union1d(
                np.nonzero(np.isnan(eastings))[0],
                np.nonzero(np.isnan(northings))[0]
                ))
        if len(kill_list) > 0:
            self.RemoveTraces(kill_list)
            dbg_traces_deleted += len(kill_list)

        if debug:
            print 'null location traces deleted:', dbg_traces_deleted

        # Loop through the rest and filter
        kill_list = []
        i = 0
        while i < self.nx:
            x = eastings[i]
            y = northings[i]

            # Find subsequent points within threshold
            dx = eastings[i:] - x
            dy = northings[i:] - y
            displacement = np.sqrt(dx**2 + dy**2)
            try:
                iend = np.nonzero(displacement > threshold)[0][0]
            except IndexError:
                iend = self.nx-i

            # Average them all together
            if iend > 1:
                if debug:
                    print i, iend, displacement[iend-1:iend+1]
                avg_trace = self.data[:,i:i+iend].sum(axis=1) / float(iend)
                dbg_traces_deleted += (iend - 1)
                # Mark spaces held by former traces for deletion
                _ = [kill_list.append(j) for j in range(i+1, i+iend)]
                for trace in range(i+1,i+iend):
                    self.retain['location_{0}'.format(trace)] = False
            else:
                avg_trace = self.data[:,i]

            # Update the working data array
            try:
                move_data = np.hstack([move_data, avg_trace.reshape([-1, 1])])
            except NameError:
                move_data = avg_trace.reshape([-1, 1])

            i += iend

        if self.nx > 0:
            self.data = move_data
            self.nx = self.data.shape[1]
            self.RemoveMetadata(kill_list)

        if debug:
            print 'stationary location traces deleted:', dbg_traces_deleted

        # I won't claim this bit isn't risky
        self.metadata_copy = copy.deepcopy(self.metadata)
        self.raw_data = self.data.copy()
        self.fids_copy = copy.copy(self.fids)

        self.history.append(('remove_stationary', threshold))
        return

    def Interpolate(self, X_int, X, arr=None):
        """ Returns array of data interpolated over space.
                X_int       locations to interpolate to
                X           current x locations
        """
        if arr is None: arr = self.data
        D_int = np.zeros([self.data.shape[0], X_int.shape[0]])
        for i, row in enumerate(arr):
            D_int[i,:] = np.interp(X_int, X, row)
        return D_int

    def LineProjectMultiSegment(self, dx=4.0, threshold=0.35, verbose=False):
        """ Projects data to a sequence of approximating line
        segments with even specing.

            dx          :   point spacing in meters
            threshold   :   threshold for segment bending
            verbose     :   print out status message

        Returns a tuple containing:
            segments
            Xdbg
            Ydbg
            Pdbg
            Pgriddbg

        THIS MAY BREAK THINGS because picks and metadata aside from
        coordinates aren't handled. Use with caution.
        """
        # Find places where line direction changed
        breaks = self.FindLineBreaks(threshold=threshold)
        sections = [(breaks[i], breaks[i+1]) for i in range(len(breaks)-1)]

        # Integrate short sections into the previous section
        fdiff = lambda a: a[1]-a[0]
        section_lengths = map(fdiff, sections)
        while min(section_lengths) < 10:
            i = section_lengths.index(min(section_lengths))
            sections[i-1] = (sections[i-1][0], sections[i][1])
            sections.pop(i)
            section_lengths = map(fdiff, sections)

        # Add 1 to the last bound so that the full array is included
        sections[-1] = (sections[-1][0], sections[-1][1]+1)
        if verbose:
            print "projecting across {0} line segments".format(len(sections))

        # Project each of the sections to a least-squares line and grid
        proj_data = np.array([])
        proj_topo = np.array([])
        proj_eastings = np.array([])
        proj_northings = np.array([])
        Xdbg, Ydbg, Pdbg, Pgriddbg = [], [], [], []

        for bounds in sections:
            # Project positions to a line
            X, Y, p = self.LineProjectXY(bounds=bounds)
            data_part = self.data[:,bounds[0]:bounds[1]]
            try:
                topo_part = self.topography[bounds[0]:bounds[1]]
            except AttributeError:      # No topography
                pass

            # Interpolate the bounded region to regular spacing
            P = np.sqrt((X-X[0])**2 + (Y-Y[0])**2)
            Pgrid = np.arange(P[0], P[-1], dx)
            proj_data_part = self.Interpolate(Pgrid, P, arr=data_part)
            try:
                proj_topo_part = np.interp(Pgrid, P, topo_part)
            except UnboundLocalError:      # No topography
                pass

            proj_x_part = np.interp(Pgrid, P, X)
            proj_y_part = np.interp(Pgrid, P, Y)

            # Append region chunk to data array
            if proj_data.size == 0:
                proj_data = proj_data_part
            else:
                proj_data = np.hstack([proj_data, proj_data_part])

            try:
                topo_part = self.topography[bounds[0]:bounds[1]]
                proj_topo_part = np.interp(Pgrid, P, topo_part)
                if len(proj_topo) == 0:
                    proj_topo = proj_topo_part
                else:
                    proj_topo = np.hstack([proj_topo, proj_topo_part])
            except AttributeError:      # No topography
                pass

            proj_eastings = np.hstack([proj_eastings, proj_x_part])
            proj_northings = np.hstack([proj_northings, proj_y_part])

            # Append debugging data (happens for each segment)
            Xdbg.extend(X.tolist())
            Ydbg.extend(Y.tolist())
            Pdbg.extend(P.tolist())
            Pgriddbg.extend(Pgrid.tolist())

        # Update state values
        self.data = proj_data
        self.metadata.eastings = proj_eastings.tolist()
        self.metadata.northings = proj_northings.tolist()
        try:
            self.topography = proj_topo
        except AttributeError:
            pass
        self.history.append(('projection_as_segments', bounds, dx))
        return (sections, Xdbg, Ydbg, Pdbg, Pgriddbg)

    def MigrateFK(self, dx=4.0, t0_adjust=0):
        """ Perform Stolt migration over multiple sections.

            dx          :   gridding interval

            t0_adjust   :   zero-time offset from top of self.data, in samples
                            e.g. if the radargram contains data from before t0,
                            this corrects for that.

        Returns a dictionary of statistics about the transformation.
        """
        full_Dmig = np.zeros(self.data.shape)

        # Find places where line direction changed
        breaks = self.FindLineBreaks(threshold = 0.35)
        migsections = [(breaks[i], breaks[i+1]) for i in range(len(breaks)-1)]

        # Integrate short sections into the previous section
        section_length = map(lambda a: a[1]-a[0], migsections)
        while min(section_length) < 10:
            i = section_length.index(min(section_length))
            migsections[i-1] = (migsections[i-1][0], migsections[i][1])
            migsections.pop(i)
            section_length = map(lambda a: a[1]-a[0], migsections)

        # Add 1 to the last bound so that the full array is included
        migsections[-1] = (migsections[-1][0], migsections[-1][1]+1)
        print "migrating in {0} sections".format(len(migsections))

        # Migrate each of the sections separately
        for bounds in migsections:
            # Project positions to a line
            X, Y, p = self.LineProjectXY(bounds=bounds)
            if t0_adjust > 0:
                arr = np.roll(self.data[:,bounds[0]:bounds[1]],
                              -t0_adjust, axis=0)
            else:
                arr = self.data[:,bounds[0]:bounds[1]]

            # Interpolate the bounded region to regular spacing
            P = np.sqrt((X-X[0])**2 + (Y-Y[0])**2)
            Pmesh = np.arange(P[0], P[-1], dx)

            if len(Pmesh) > 10:
                proj_arr = self.Interpolate(Pmesh, P, arr=arr)

                # Perform migration
                #pdb.set_trace()
                Dmig, tmig, xmig = fkmig(proj_arr, self.rate, dx, 1.68e8)

                # Interpolate back
                Dmig = self.Interpolate(P, Pmesh, arr=Dmig)

                # Place migrated segment in data array
                full_Dmig[:,bounds[0]:bounds[1]] = Dmig

            else:
                # Section does not have enough horizontal displacement
                # for migration to be meaningful
                full_Dmig[:,bounds[0]:bounds[1]] = 0.0

        if t0_adjust > 0:
            self.data = np.roll(full_Dmig, t0_adjust, axis=0)
        else:
            self.data = full_Dmig
        self.history.append(('fk_migration', bounds, dx))
        return migsections

    def Plot(self, outfile=None, cmap='gray', title='Radar Line',
        rate=1e-8, c=1.68e8):
        """ Plot an array containing radar data along a line. """
        n = self.data.shape[0]
        T = np.arange(0, n*rate, rate)      # time axis

        # Find the luminescence range so that plot intensity is symmetric
        lum_bound = max((abs(self.data.max()), abs(self.data.min())))

        # Draw it
        img = plt.imshow(self.data, aspect='auto', cmap=cmap,
                         vmin=-lum_bound, vmax=lum_bound)
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


class CommonMidpointGather(Gather):
    """ Subclass defining common-midpoint specific data and operations. """

    def ReadIndex(index):
        """ Read CMP index file and return the offsets as a 1-d array and
        location key as a 2-d array. """
        with open(index, 'r') as f:
            lines = f.readlines()
        # Remove title row
        lines.pop(0)
        xr = []
        xt = []
        loc_0 = []
        loc_f = []
        for line in lines:
            try:
                s = line.split(',')
                xr.append(float(s[0]))
                xt.append(float(s[1]))
                loc_0.append(int(s[2]))
                loc_f.append(int(s[3]))
            except:
                traceback.print_exc()
                return
        XR = np.array(xr)
        XT = np.array(xt)
        KEY = np.array(zip(loc_0, loc_f))
        OFFSETS = np.abs(XR - XT)
        return OFFSETS, KEY

    def CalcAveragePicks(key, picks):
        """ Average picks by location, based on a key. The key is a list of
        paired iterables, containing first and last locations for every shot
        gather.
        """
        avg_picks = []
        avg = lambda L: sum(L) / len(L)
        try:
            for pair in key:
                avg_picks.append(avg(picks[pair[0]:pair[1]+1]))
        except IndexError:
            sys.stderr.write("key and pick arguments to CalcAveragePicks() are"
                             " inconsistent\n")
            avg_picks = None
        return np.array(avg_picks)

    def Plot(offsets, events, curve=None):
        """ Plot picked events. If a set of points defining a reference NMO
        hyperbola is provided, plot that too. """
        plt.ion()
        plt.plot(offsets, -events, 'xk')
        if analytical is not None:
            try:
                plt.plot(offsets, -curve, '-r')
            except:
                pass
        return


class LineGather(CommonOffsetGather):
    """ Class for handling individual survey lines.

        arr         data held within a numpy array
        infile      the original HDF5 dataset
        line        the line enumeration
        metadata    RecordList instance containing metadata

        DEPRECATED
    """
    pass


class LineGatherError(Exception):
    def __init__(self, message="No message"):
        self.message = message
    def __str__(self):
        return self.message
