"""
Functions relevant to the analysis of bed reflection power (BRP) in radar
traces.
"""

import os
import numpy as np

def get_pickfnm(L, directory):
    """ Return a pickfile name given a Gather `L` and a base directory. """
    filenamestem = os.path.splitext(os.path.basename(L.infile))[0]
    loadfile = os.path.join(directory, filenamestem + '_line' + str(L.line) \
                                       + '.csv')
    return loadfile

def extract_window(tr, i, relwin):
    """ Extract a window from vector `tr` at entry `i` of size tuple `relwin`
    """
    try:
        out = tr[i-relwin[0]:i+relwin[1]]
    except TypeError:
        out = np.nan
    return out

def brp_windowed(L, relwin=(None, None)):
    """ Calculate the bed reflection power within a window relative to the
    picks in Gather `L`. """
    if None in relwin:
        relwin = (int(round(-2e8*L.rate)), int(round(1e9*L.rate)))
    power = lambda a: np.sum(a**2) * L.rate
    bed_wavelets = map(lambda tr, i: extract_window(tr, i, relwin),
                       L.data.T, L.bed_picks)
    brp = map(power, bed_wavelets)
    return brp

