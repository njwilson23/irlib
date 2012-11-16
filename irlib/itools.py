#! /usr/bin/python

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import scipy.io
import pywavelet
from gather import Gather, CommonOffsetGather, LineGatherError
from recordlist import RecordList


#font = FontProperties(fname='Verdana.ttf', weight='normal', size=12)
font = FontProperties(family='sans-serif', weight='normal', size=11)

def plotax(ax, L, gain=5, annotate=True, font=None, nan_fill=None):
    """ Make a radargram plot.

    plotax(ax, L, gain=5, annotate=True)
    where *ax* is a matplotlib.Axes, *L* is an irlib.Gather, *gain* specifies
    the contrast, and *annotate* turns on axes labels.

    *font*          : annotation FontProperties

    *mask_topo*     : if not None, then must be a value to replace nans with

    Replaces plotl() and plotlt() functions.
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

    ax.imshow(data, aspect='auto', cmap='gray', vmin=-lbnd, vmax=lbnd)

    if annotate:
        ax2 = ax.twinx()
        ax.set_ylabel('Time (ns)', fontproperties=font)
        ax2.set_ylabel('Sample number', rotation=270, fontproperties=font)
        ax.set_xlabel('Trace number', fontproperties=font)
        ax.set_xticklabels([int(a) for a in ax.get_xticks()], fontproperties=font)
        timevec = ax.get_yticks() * L.rate * 1e9
        ax.set_yticklabels(timevec.astype(int), fontproperties=font)
        locs2 = ax2.get_yticks()
        labels2 = [int(round(i)) for i in locs2*L.data.shape[0]]
        labels2.reverse()
        ax2.set_yticklabels(labels2, fontproperties=font)

    return ax


def plotl(D, gain=5, fnm=None, for_pub=False):
    """ Make a radargram plot.
        plotl(Gather, gain=1.25)
        DEPRECATED - use *plotax* instead
    """
    if for_pub:
        fig = plt.figure(figsize=[7,2])
        ax1 = fig.add_axes([0.15,0.22,0.83,0.7])
    else:
        fig = plt.figure(figsize=[9,9])
        ax1 = fig.add_axes([0.1,0.1,0.85,0.8])
        ax2 = ax1.twinx()

    if isinstance(D, np.ndarray):
        data = D
        rate = 4e-9
    elif isinstance(D, Gather):
        data = D.data
        rate = D.rate

    lbnd = max([-data.min() / gain, data.max() / gain])

    ax1.imshow(data, aspect='auto', cmap='gray', vmin=-lbnd, vmax=lbnd)
    ax1.set_ylabel('Time (ns)')
    ax1.set_xticklabels([int(a) for a in ax1.get_xticks()])
    ax1.set_xlabel('Trace number')

    locs1 = ax1.get_yticks()
    ax1.set_yticklabels([locs1[i] * rate * 1e9 for i in range(len(locs1))], fontproperties=font)
    if not for_pub:
        ax2.set_ylabel('Sample number')
        locs2 = ax2.get_yticks()
        labels2 = [int(round(i)) for i in locs2*data.shape[0]]
        labels2.reverse()
        ax2.set_yticklabels(labels2, fontproperties=font)

    if fnm is not None:
        fig.savefig(fnm)
    return fig

def plotlt(L, gain=5, fnm=None, for_pub=False, small=False):
    """ Make a radargram plot. Use topography.
        plotlt(Gather, gain=1.25)
        DEPRECATED - use *plotax* instead
    """
    if for_pub:
        figsize=[6,2]
        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_axes([0.02,0.05,0.96,0.9])
    else:
        figsize=[8,3]
        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_axes([0.05,0.05,0.95,0.92])

    ax2 = ax1.twinx()

    lbnd = max([-L.data.min() / gain, L.data.max() / gain])
    data = L.GetTopoCorrectedData()
    ax1.imshow(data, aspect='auto', cmap='gray', vmin=-lbnd, vmax=lbnd)

    if small:
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax2.set_yticks([])

    else:
        ax1.set_ylabel('Time (ns)')
        ax2.set_ylabel('Sample number')
        ax1.set_xlabel('Trace number')
        ax1.set_xticklabels([int(a) for a in ax1.get_xticks()])
        locs1 = ax1.get_yticks()
        ax1.set_yticklabels([locs1[i] * L.rate * 1e9 for i in range(len(locs1))], fontproperties=font)
        locs2 = ax2.get_yticks()
        labels2 = [int(round(i)) for i in locs2*L.data.shape[0]]
        labels2.reverse()
        ax2.set_yticklabels(labels2, fontproperties=font)

    if fnm is not None:
        fig.savefig(fnm)
    return fig

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
    ax.contourf(logwva[::-1,:], vmax=vmax, cmap=plt.cm.gray)

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
            j = random.randint(0, nmax-1)
        ns.append(j)
    ns.pop(0)
    return ns

def split_sample(fin, fout1, fout2, n):
    """ Extract n (int) lines from fin (file) and write to
    fout2 (file). Write the remaining lines to fout1 (file). """
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
    print ind
    return find_nearest1d(Pmesh, P[ind])

def read_cresis_mat(matfile):
    """ Read in a CRESIS Antarctica/Greenland MATLAB datafile as a
    CommonOffsetGather. """
    C = scipy.io.loadmat(matfile)
    R = RecordList(None)
    R.lats = C['Latitude']
    R.lons = C['Longitude']
    R.timestamps = C['GPS_time']
    R.sample_rate = 1.0 / np.diff(C['Time'][0]).mean() * np.ones(C['Data'].shape[0])
    return CommonOffsetGather(C['Data'], line=np.nan, metadata=R)
