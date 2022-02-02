#! /usr/bin/env python

"""
Functions relevant to the analysis of bed reflection power (BRP) in radar
traces.

Author Nat Wilson
Improvements by Yulia Antropova 
TODO - continue to implement and test


"""
import irlib
import os
import numpy as np
import pandas as pd
import argparse

def get_pickfnm(L, directory):
    """ Return a pickfile name given a Gather `L` and a base directory. """
    filenamestem = os.path.splitext(os.path.basename(L.infile))[0]
    loadfile = os.path.join(directory, filenamestem + '_line' + str(L.line) \
                                       + '.csv')
    return loadfile

def extract_window_around(tr, i, relwin):
    """ Extract a window from vector `tr` at entry `i` of size tuple `relwin`
    """
    try:
        out = tr[i+relwin[0]:i+relwin[1]]
    except TypeError:
        out = np.nan
    return out

def get_brp_windowed(L, relwin=(None, None)):
    """ Calculate the bed reflection power within a window (specified in
    samples) relative to the picks in Gather `L`. """
    if None in relwin:
        #relwin = (int(round(-1e9*L.rate)), int(round(5e9*L.rate)))
        # Use Copland and Sharp's window (100 ns prior, 400 ns after)
        # or modify as (100 ns prior, 100 ns after): 
        relwin = (int(round(-100e-9/L.rate)), int(round(100e-9/L.rate)))
    power = lambda a: np.sum(a**2) * L.rate
    bed_wavelets = map(lambda tr, i: extract_window_around(tr, i, relwin),
                       L.data.T, L.bed_picks)
    brp = map(power, bed_wavelets)
    return brp

def get_irp(L):
    """ Calculate the bed reflection power in line `L` between the airwave and
    bed wave. 

    Note: IRP will be calculated for traces where no picks were done, to correct
    for that delete all IRP values that correspond to blank BRP values.    
    
    """
    power = lambda a: np.sum(a**2) * L.rate
    # Yulia had this as F.dcvals and F.bedvals - not sure where F comes from?! so used L
    L.dcvals =[x if str(x)!='nan' else 0 for x in L.dcvals]
    air = L.dcvals
    L.bedvals =[x if str(x)!='nan' else 0 for x in L.bedvals]
    bed = L.bedvals

    # Amount to avoid airwave/bedwave by, in samples:
    pad = int(round(100e-9 / L.rate)) # 100 ns
    internal_windows = [L.data[a+pad:b-pad] for a,b in zip(air,bed)]
    irp = map(power, internal_windows)
    # irp = map(power, internal_windows if {a!='nan' or b!='nan'} else None)
 
    return irp


def makeBRPdf(L):
    '''
    This function may be useful for making a dataframe but it needs work
    need to call irp and brp functions

    Parameters
    ----------
    L : TYPE
        DESCRIPTION.

    Returns
    -------
    brp_df : TYPE
        DESCRIPTION.

    '''
    
    traces = [int(s[4:8]) for s in L.metadata.fids]
    fids = []
    lons = []
    lats = []
    alt = []
    
    for fid in F.fids:
        n = int(fid[4:8])
        fids.append (fid)
        lons.append (L.metadata.lons[traces.index(n)])
        lats.append (L.metadata.lats[traces.index(n)])
        alt.append (L.metadata.alt_asl[traces.index(n)])
        
    brp_df = pd.DataFrame(list(map(list, zip(fids,lons,lats,alts, brp, irp))))    
    
    return brp_df





if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description = "BRP Analysis")
    parser.add_argument("infile", help="input HDF (*.h5) filename, with or without path")
    parser.add_argument("brpfile", help="output BRP filename")
    parser.add_argument("pickdir", help="path to picking directory")
    parser.add_argument("brpdir", help="path to BRP output directory")
        
    args = parser.parse_args()

    pick_directory = args.pick_directory
    fname = args.infile
    brp_fname = args.brpfile
    brp_path = args.brpdir
    
    S = irlib.Survey(fname)
    brp = pd.DataFrame()

    for line_number in range(len(S.GetLines())):

        L = S.ExtractLine(line_number, fromcache=True)
    
        if os.path.exists(get_pickfnm(L, pick_directory)):

            F = irlib.FileHandler(get_pickfnm(L, pick_directory), line_number)
            F.bedvals
        
            ##To account for values, where no picking were done 
            #(the resulting BRP will a blank cell):
            F.bedvals =[x if str(x)!='nan' else 'nan' for x in F.bedvals]
    
            # Open radar line to retrieve UTM coords
            L = S.ExtractLine(line_number)
        
            sample_rate = L.rate
    
            brp_new = get_brp_windowed(L)
            print("BRP Complete for line {}\n".format(line_number))
            #print (brp)
            brp_all = pd.concat([brp, brp_new])
            brp = brp_all
        
    brp_out = "BRP_IRP_centred_100_100" + str(brp_fname) + ".csv"
    brp_path = os.path.join(brp_path, brp_out)
    brp_all.columns = ['fid','longitude', 'latitude', 'brp', 'irp']
    brp_all.to_csv(brp_path, sep=',',index=False)
    
    # could also make a shapefile quite easily here... 
    