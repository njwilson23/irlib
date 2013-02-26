"""
Functions relevant to the analysis of bed reflection power (BRP) in radar
traces.
"""

import os
import numpy as np

def get_pickfnm(L, directory):
    filenamestem = os.path.splitext(os.path.basename(L.infile))[0]
    loadfile = os.path.join(directory, filenamestem + '_line' + str(L.line) \
                                       + '.csv')
    return loadfile

def extract_wavelet(tr, i, relwin):
    """ Extract a window from vector `tr` at entry `i` of size tuple `relwin`
    """
    try:
        out = tr[i-relwin[0]:i+relwin[1]]
    except TypeError:
        out = np.nan
    return out

def brp_windowed(L, relwin=[None, None]):
    """ Calculate the bed reflection power within a window relative to the
    picks in Gather `L`. """
    if relwin[0] is None:
        relwin[0] = -2e8*L.rate
    if relwin[1] is None:
        relwin[1] = 1e9*L.rate
    picks = L.bed_picks
    power = lambda a: np.sum(a**2)
    bed_wavelets = map(lambda tr, i: extract_wavelet(tr, i, relwin),
                       L.data.T, L.bed_picks)
    brp = map(power, bed_wavelets)
    return brp

