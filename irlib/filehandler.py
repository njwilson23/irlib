""" Contains the `FileHandler` class, which abstracts the process of reading
and writing picking files. """

import os
import numpy as np

class FileHandler:
    """ Class for reading and writing pick and rating files cleanly.
    If no list of FIDs is provided, then they will be guessed
    assuming locations are in order and complete. If a list of FIDs
    is provided, it will be used instead.
    """
    def __init__(self, fnm, line, fids=None):
        self.fnm = fnm
        self.line = line

        if os.path.isfile(fnm):
            try:
                with open(fnm) as f:
                    _ = f.readline()            # Bypass the header
                    self.prevrecs = f.readlines()
                    self.Parse(self.prevrecs)
            except:
                traceback.print_exc()
        elif fids is None:
            raise IOError("Filehandler must either recieve an existing file "
                          "or a list of FIDs\n"
                          "  {0} is not a valid file".format(fnm))
        else:
            self.prevlines = None
            self.fids = fids
            self.dcvals = [np.nan for a in self.fids]
            self.bedvals = [np.nan for a in self.fids]
            self.traveltimes = [np.nan for a in self.fids]

    @property
    def nrecs(self):
        return len(self.fids)

    def sort(self):
        """ Sort bedvals and dcvals by FID in-place. """
        def sortby(a, b):
            return [a_ for (b_, a_) in sorted(zip(b,a))]
        self.bedvals = sortby(self.bedvals, [int(a) for a in self.fids])
        self.dcvals = sortby(self.dcvals, [int(a) for a in self.fids])
        self.fids.sort(key=lambda a: int(a))
        return

    def Parse(self, recs):
        """ Read pick file records, and split into fields. """
        self.fids = []
        self.dcvals = []
        self.bedvals = []
        self.traveltimes = []
        for row in recs:
            if row != '':
                row = row.rstrip('\n')
                self.fids.append(row.split(',')[0])
                self.dcvals.append(float(row.split(',')[1]))
                self.bedvals.append(float(row.split(',')[2]))
                self.traveltimes.append(float(row.split(',')[3]))

    def GetEventVals(self):
        """ Return the airwave and bed reflection values (lists). """
        if self.nrecs == 0:
            raise FileHandlerError('Event file does not exist: {0}'.format(self.fnm))
        return self.dcvals, self.bedvals

    def GetEventValsByFID(self, fids):
        """ Return the airwave and bed reflection picks for a list of FIDs.
        """
        if self.nrecs == 0:
            raise FileHandlerError('Event file does not exist: {0}'.format(self.fnm))
        dcvals, bedvals = [], []
        if hasattr(fids, '__iter__'):
            for fid in fids:
                dcvals.append(searchbylist(fid, self.fids, self.dcvals))
                bedvals.append(searchbylist(fid, self.fids, self.bedvals))
        else:
            dcvals.append(searchbylist(fids, self.fids, self.dcvals))
            bedvals.append(searchbylist(fids, self.fids, self.bedvals))
        return dcvals, bedvals

    def GetEventVals_Interpolated(self, max_fid=None):
        """ Similar to GetEventVals, however returns a value for
        every FID. If max_fid is given, it is taken to be the ending
        FID, otherwise the last value in self.fids will be used as
        the ending FID.

        self.fids must not be None, otherwise this won't work.

        This method is only valid if missing records can be assumed
        to be linearly related to existing records. A less naive
        approach is to use locations from a line's metadata to do
        this interpolation.
        """
        try:
            assert self.fids is not None
        except AssertionError:
            raise FileHandlerError("Cannot interpolate without defined FIDs")

        if max_fid is None:
            max_fid = self.fids[-1]

        locs = [int(s[4:8]) for s in self.fids]
        max_loc = int(max_fid[4:8])
        all_locs = range(max_loc+1)
        bedvals_interpolated = np.interp(all_locs, locs, self.bedvals)
        dcvals_interpolated = np.interp(all_locs, locs, self.dcvals)
        return dcvals_interpolated, bedvals_interpolated

    def AddBedPicks(self, fids, vals):
        """ Add reflection picks at locations given by FIDs """
        for fid, val in zip(fids, vals):
            if fid in self.fids:
                self.bedvals[self.fids.index(fid)] = val
            else:
                self.fids.append(fid)
                self.bedvals.append(val)
        self.sort()
        return

    def AddDCPicks(self, fids, vals):
        """ Add direct wave picks at locations given by FIDs """
        for fid, val in zip(fids, vals):
            if fid in self.fids:
                self.dcvals[self.fids.index(fid)] = val
            else:
                self.fids.append(fid)
                self.dcvals.append(val)
        self.sort()
        return

    def ComputeTravelTimes(self):
        """ Where possible, subtract dc times from bed times. """
        try:

            if self.nrecs == 0:
                return
		# changed nodata=999 to -999, since 999 is a legitimate value (DM)
            self.traveltimes = [-999 for i in range(self.nrecs)]  

            for i in range(self.nrecs):
                if self.bedvals[i] != -999 and self.dcvals[i] != -999:
                    self.traveltimes[i] = self.bedvals[i] - self.dcvals[i]
        except:
            traceback.print_exc()

    def Write(self):
        """ Write to file. """
        savedir = os.path.split(self.fnm)[0]
        if not os.path.isdir(savedir):
            os.makedirs(savedir)
        with open(self.fnm, 'w') as f:
            f.write('"FID","dc","bed","trav_time"\n')
            data = zip(self.dcvals, self.bedvals, self.traveltimes)

            if self.fids is None:
                sys.stderr.write("Warning: guessing FIDS in FileHandler.Write()\n")

            for i,row in enumerate(data):
                fid = self.fids[i]
                f.write('{fid},{dcval},{bedval},{tt}\n'.format(
                    fid=fid, dcval=row[0], bedval=row[1], tt=row[2]))


def searchbylist(key, keylist, vallist, notfound=np.nan):
    """ Search a `keylist` for a `key`, and return the corresponding value from
    `vallist`. If `key` is not in `keylist`, return `notfound`.
    """
    out = notfound
    for k,v in zip(keylist, vallist):
        if k == key:
            out = v
            break
    return out


class FileHandlerError(Exception):
    def __init__(self, message="No message"):
        self.message = message
    def __str__(self):
        return self.message
