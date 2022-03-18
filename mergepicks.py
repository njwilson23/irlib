#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mergepicks.py

This commandline script will take an UPDATED h5 file/cache/picking directory 
(maybe it was reprocessed with better GPS data or different filters, or somebody
started repicking it, etc) and generate picking files in a NEW picking directory 
that are merged with some OLDER picking files found in a separate directory

When merging the OLD picking files with the UPDATED h5/cache/picking files, 
you can select to retain the OLD picks or the UPDATED picks in case a trace 
has been picked in both versions. 

The script looks for picking files in the picking folder and will use them as UPDATED picks, if found
If not, it will generate picking files from the cache (if found), otherwise, 
picking files will be generated from the h5 files. 

Created on Sat May 29 23:12:34 2021
@author: dmueller
"""
import irlib
import os, sys, argparse
import traceback
import pandas as pd
import os
import numpy as np
import pdb

def get_pick_fnm(h5,lineno,dirpicking):
    """ Autogenerate a filename for pickfiles. (modified from components.py)"""
    fnm = os.path.join(dirpicking, '{0}_line{1}.csv'.format(
        os.path.basename(h5).split('.')[0],lineno))
    return fnm

def save_picks(pickableline, lineno, fnm):
    """ Save picks. (modified from components.py)
    pickableline - a pickablegather object
    lineno - a line number (string ok ('0'), but not 'line_0'
    fnm - picking file path and name 
    
    """
    fh = irlib.FileHandler(fnm, lineno, fids=pickableline.fids)
    fh.AddBedPicks(pickableline.fids, pickableline.bed_picks)
    fh.AddDCPicks(pickableline.fids, pickableline.dc_picks)
    fh.ComputeTravelTimes()
    fh.Write()
    return

#replacing getopt and def syntax() with argparse
prog_description = """
    SYNTAX: mergepicks infile outdir oldpicks [OPTIONS]
    
    Note this script will overwrite picking files in the 'outdir' and 
      'outdir' cannot be named 'picking'
         
    """
parser = argparse.ArgumentParser(description = prog_description)
parser.add_argument("infile", help="input HDF (.h5) filename")
parser.add_argument("outdir", help="subfolder where new picking files will be written")
parser.add_argument("oldpicks", help="folder where old picking files are found")
parser.add_argument("-d", "--dir_cache", help="cache directory (default: cache/)", default='cache/')
parser.add_argument("-n","--newpick_priority", help="will priviledge new picks over old picks in case of conflict", action="store_true")
parser.add_argument("--dc", help="specify datacapture (default: 0)", default=0,type=int)
#parser.add_argument("--clobber", help="Allow files to be overwritten", action="store_true")
#TODO = add the clobber logic here eventually so only runs with clobber set overwrite files

args = parser.parse_args()

h5 = args.infile 
wrkdir = os.path.dirname(h5)
dircache = args.dir_cache
oldpickdir = args.oldpicks
if os.path.isdir(oldpickdir):
    # assume file pattern is all the same
    oldh5 = os.listdir(oldpickdir)[0].split('_line')[0]+'.h5'
else:
    print("Cannot find folder with old picking files")
    sys.exit(1)
if args.outdir == 'picking':
    print("You cannot specify 'picking' as an outdir... quitting")
    sys.exit(1)
    
if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

# open h5 file.   
S = irlib.Survey(h5)#args.infile)
# list all the lines in the survey...     
linenames = S.GetLines()
linenos = [L.split(sep="_")[1] for L in linenames]

for line in linenos: 
    print('Line {}....'.format(line))
    cfile = S.GetLineCacheName(line, cache_dir=dircache)
        
    #cache directory and file separated:
    cdir = os.path.join(wrkdir, os.path.dirname(cfile))
    cfile = os.path.basename(cfile)    
    
    # extract line from survey
    if os.path.isfile(os.path.join(cdir, cfile)):
        L = S.ExtractLine(line,fromcache=True,cache_dir=cdir, datacapture=args.dc)
    else:
        L = S.ExtractLine(line,fromcache=False, datacapture=args.dc)
    
    if L == None:
        print(' Line {} in current survey version could not be opened'.format(line))
        continue
    
    #create a pickable p_updated.bed.count() gather object
    p = irlib.PickableGather(L)
    
    if os.path.isfile(get_pick_fnm(h5, line, os.path.join(wrkdir,'picking'))):
        p.LoadPicks(get_pick_fnm(h5, line, os.path.join(wrkdir,'picking')))
    
    p_updated = pd.DataFrame(list(zip(p.fids, list(p.dc_picks), list(p.bed_picks))), 
                                 columns =['FID', 'dc', 'bed'])
    
    print(' Current survey version has {} traces, with {} dc picks and {} bed picks'.format(
            p_updated.FID.count(), p_updated.dc.count(), p_updated.bed.count()))    

    try: 
        #read in from csv
        p_old = pd.read_csv(get_pick_fnm(oldh5, line, oldpickdir))
        #recreate the fid left zero padding
        p_old.FID = ['{:0>16}'.format(fid) for fid in p_old.FID]
        print(' Older survey version has {} traces, with {} dc picks and {} bed picks'.format(
            p_old.FID.count(), p_old.dc.count(), p_old.bed.count()))    

        #merge dataframes
        p_merge = pd.merge(p_updated,p_old, on="FID", how='left')
        if args.newpick_priority:
            # when you want all the newpicks and only oldpicks when new is NaN
            p_merge.dc_x[p_merge.dc_x.issnull()] = p_merge.dc_y[p_merge.dc_x.isnull()]
            p_merge.bed_x[p_merge.bed_x.issnull()] = p_merge.bed_y[p_merge.bed_x.isnull()]
        else:
            # only take old picks
            p_merge.dc_x = p_merge.dc_y
            p_merge.bed_x = p_merge.bed_y
            
        print(' Merged survey version has {} traces, with {} dc picks and {} bed picks'.format(
            p_merge.FID.count(), p_merge.dc_x.count(), p_merge.bed_x.count()))    
            
        p.bed_picks = np.array(p_merge.bed_x)
        p.dc_picks = np.array(p_merge.dc_x)
        save_picks(p, line, get_pick_fnm(h5, line, args.outdir))
    except:
        save_picks(p, line, get_pick_fnm(h5, line, args.outdir))
        print(' Could not find/read file {}. New [unmerged] picking file saved'.format(get_pick_fnm(oldh5, line, oldpickdir)))     
            
    print('\n')
        