""" Convenience functions for dealing with radar data. This module contains
bits of code that aren't good enough to be in `irlib` proper, but that I got
tired to rewriting to solve actual problems. """

import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import scipy.io
try:
    import pywavelet
except ImportError:
    pass
from gather import Gather, CommonOffsetGather, LineGatherError
from recordlist import RecordList

#from irlib.pEKKOdriver import read_pulseEKKO
from pEKKOdriver import read_pulseEKKO


font = FontProperties(family='sans-serif', weight='normal', size=11)
defaultcm = matplotlib.cm.gray

def plotax(ax, L, gain=5, annotate=True, font=None, nan_fill=None,
           cmap=matplotlib.cm.gray):
    """ Make a radargram plot. Replaces `plotl()` and `plotlt()` functions.

    Parameters
    ----------
    ax : a matplotlib.Axes
    L : an irlib.Gather
    gain : specifies the display contrast
    annotate : turns on axes labels
    font : annotation FontProperties
    mask_topo : if not None, then must be a value to replace nans with

    Returns
    -------
    ax : same Axes instance as passed in
    """
    rate = L.rate
    lbnd = max([-L.data.min() / gain, L.data.max() / gain])

    if font is None:
        font = matplotlib.font_manager.FontProperties()

    try:
        data = L.GetTopoCorrectedData()
        if nan_fill is not None:
            data[np.isnan(data)] = nan_fill
    except LineGatherError:
        data = L.data

    X, Y = np.meshgrid(np.arange(L.nx), L.rate * np.arange(L.ny))
    ax.pcolormesh(X, Y, data, cmap=cmap, vmin=-lbnd, vmax=lbnd, rasterized=True)
    ax.set_ylabel('Time (ns)', fontproperties=font)
    ax.set_xlabel('Trace number', fontproperties=font)
    ax.set_xticklabels([int(a) for a in ax.get_xticks()], fontproperties=font)
    ax.set_ylim(sorted(ax.get_ylim(), reverse=True))
    ax.set_yticklabels((ax.get_yticks()*1e9).astype(int), fontproperties=font)
    ax.axis('tight')

    if annotate:
        ax2 = ax.twinx()
        ax2.set_ylabel('Sample number', rotation=270, fontproperties=font)
        locs2 = ax2.get_yticks()
        labels2 = [int(round(i)) for i in locs2*L.data.shape[0]]
        labels2.reverse()
        ax2.set_yticklabels(labels2, fontproperties=font)

    return ax


def plotwv(wva, scales, vlim=None, fnm=None, for_pub=False):
    """ Plot a wavelet transform. """
    # Calculate fourier wavelength and frequency from scales
    # Use w0=6 (Torrence and Compo, 1998)
    radicand = np.sqrt(2+6**2)
    flamda = 4*np.pi*scales / (6+radicand)
    freq = 1.0/flamda

    logwva = np.log(abs(wva.T))

    if vlim is not None:
        vmax = np.log(vlim)
    else:
        vmax = logwva.max()

    if for_pub:
        fig = plt.figure(figsize=[4,2])
    else:
        fig = plt.figure(figsize=[6,3])
    ax = fig.add_axes([0.13, 0.1, 0.86, 0.89])
    ax.contour(logwva[::-1,:], vmax=vmax, colors="black", linestyles="solid")
    ax.contourf(logwva[::-1,:], vmax=vmax, cmap=defaultcm)

    locs = np.linspace(10, len(freq)-11, 6).astype(int)

    ax.set_yticks(locs)
    ax.set_yticklabels((freq[locs[::-1]] / 1e6).astype(int), fontproperties=font)
    ax.set_ylabel('Frequency (Mhz)', fontproperties=font)
    ax.set_xticklabels([int(a) for a in ax.get_xticks()])
    ax.set_xlabel('Sample number', fontproperties=font)
    if fnm is None:
        plt.show()
    else:
        fig.savefig(fnm)
    return fig

def wvlt_trans(L, traceno, show=False, gain=1.0):
    """ Do a wavelet transform of a given trace, and plot it.
        wvlt_trans(Gather, traceno, show=False)
    """
    # Wavelet transform
    Z1 = L.data[:,traceno].copy()
    wva, scales, periods, coi = L.WaveletTransform(traceno)
    # Make pretty
    plt.subplot(2,1,1)
    T = np.arange(Z1.size)*4e-9
    plt.plot(T, Z1*1000, '-k')
    plt.xlim(0, T[-1])
    plt.ylim(-max(abs(Z1*1000)), max(abs(Z1*1000)))
    plt.ylabel('Voltage (mV)')
    plt.xticks([-10])
    plt.subplot(2,1,2)
    wvpower = np.abs(wva.T)
    plt.imshow(wvpower, vmin=wvpower.min(), vmax=wvpower.max()/gain)
    plt.axis('tight')
    freq = 1.0 / (scales)
    locs = range(0, len(freq), 25)
    plt.yticks(locs, [int(round(i)) for i in freq[locs] / 1e6])
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Sample number')
    if show:
        plt.show()

dt2m = lambda x1, x2: abs(x1-x2) * 1e-9 * 1.68e8 / 2.0

ds2m = lambda x1, x2: abs(x1-x2) * 4e-9 * 1.68e8 / 2.0

find_nearest1d = lambda L,x: np.nonzero(abs(L-x)==min(abs(L-x)))

def find_nearest2d(x, y, X, Y):
    D = np.sqrt((X-x)**2 + (Y-y)**2)
    return np.nonzero(D==D.min())

def get_n_indices(n, lst):
    """ Return n (int) indices chosen at random from lst (list-like). Disallow
    repeats. """
    nmax = len(lst)
    ns = [-1]
    for i in range(n):
        j=-1
        while j in ns:
            j = np.random.randint(0, nmax-1)
        ns.append(j)
    ns.pop(0)
    return ns

def split_sample(fin, fout1, fout2, n):
    """ Split a file object into two others. Extract `n` (int) lines from `fin`
    (file) and write to `fout2` (file). Write the remaining lines to `fout1`
    (file). This was useful, once. """
    lines = fin.readlines()
    I = sorted(get_n_indices(0, len(lines)-1))
    smp = [lines[i] for i in I]
    fout2.writelines(smp)
    I.reverse()
    for i in I: lines.pop(i)
    fout1.writelines(lines)
    return fout1, fout2

def find_projected_nearest_to(x, y, D):
    """ This is the worst kind of side-effect-ridden code. Sorry. """
    pX = np.array(D[1])
    pY = np.array(D[2])
    P = np.array(D[3])
    Pmesh = np.array(D[4])

    def remove_breaks(X):
        for i in range(1, len(X)):
            if X[i] < X[i-1]:
                X[i:] = X[i:] + X[i-1] - X[i]

        return X

    P = remove_breaks(P)
    Pmesh = remove_breaks(Pmesh)

    ind = find_nearest2d(x, y, pX, pY)
    print (ind)
    return find_nearest1d(Pmesh, P[ind])

def read_cresis_mat(matfile):
    """ Read in a CRESIS Antarctica/Greenland MATLAB datafile as a
    `CommonOffsetGather`. """
    C = scipy.io.loadmat(matfile)
    R = RecordList(None)
    R.lats = C['Latitude']
    R.lons = C['Longitude']
    R.timestamps = C['GPS_time']
    R.sample_rate = 1.0 / np.diff(C['Time'][0]).mean() * np.ones(C['Data'].shape[0])
    return CommonOffsetGather(C['Data'], line=0, metadata=R)

def read_pEKKO_as_gather(fstem):
    lmeta, tmeta, darray = read_pulseEKKO(fstem)
    metadata = RecordList(None)
    smple_rate = float(lmeta.get("NOMINAL FREQUENCY", np.nan))
    metadata.sample_rate = [smple_rate for i in range(len(tmeta))]
    return CommonOffsetGather(darray, line=0, metadata=metadata)
