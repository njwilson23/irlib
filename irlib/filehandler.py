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
                    self.nrecs = len(self.prevrecs)
                    self.Parse(self.prevrecs)
            except:
                traceback.print_exc()
        else:
            self.prevlines = None
            self.nrecs = None
            self.fids = fids
            self.dcvals = None
            self.bedvals = None
            self.traveltimes = None

    def _linloc2fid(self, lin, loc):
        """ Based on a path, return a unique FID for database
        relations. """
        try:
            dc = 0
            eg = 0
            fid = str(lin).rjust(4,'0') + str(loc).rjust(4,'0') \
                + str(dc).rjust(4,'0') + str(eg).rjust(4,'0')
            return fid
        except:
            traceback.print_exc()
            sys.stderr.write('SaveFile: failed at linloc\n')
            return None

    def Parse(self, recs):
        """ Read pick file records, and split into fields. """
        self.fids = []
        self.dcvals = []
        self.bedvals = []
        self.traveltimes = []
        for row in recs:
            if row is not '':
                row = row.rstrip('\n')
                self.fids.append(row.split(',')[0])
                self.dcvals.append(float(row.split(',')[1]))
                self.bedvals.append(float(row.split(',')[2]))
                self.traveltimes.append(float(row.split(',')[3]))

    def GetEventVals(self):
        """ Return the airwave and bed reflection values (lists). """
        if self.nrecs is None:
            raise FileHandlerError('Event file does not exist: {0}'.format(self.fnm))
        return self.dcvals, self.bedvals

    def GetEventValsByFID(self, fids):
        """ Return the airwave and bed reflection picks for a list of FIDs.
        """
        if self.nrecs is None:
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

    def AddBedPicks(self, vec):
        # SHOULD BE FIXED TO CHECK IF FIDS ARE THE SAME
        if self.nrecs is None:
            self.nrecs = len(vec)
            self.dcvals = [999 for i in range(self.nrecs)]
        if len(vec) == self.nrecs:
            self.bedvals = list(vec)
        else:
            raise FileHandlerError('SaveFile: data lengths inconsistent (bed)\n')

    def AddDCPicks(self, vec):
        # SHOULD BE FIXED TO CHECK IF FIDS ARE THE SAME
        if self.nrecs is None:
            self.nrecs = len(vec)
            self.bedvals = [999 for i in range(self.nrecs)]
        if len(vec) == self.nrecs:
            self.dcvals = list(vec)
        else:
            raise FileHandlerError('SaveFile: data lengths inconsistent (dc)\n')

    def ComputeTravelTimes(self):
        """ Where possible, subtract dc times from bed times. """
        try:

            if self.nrecs is None:
                return

            self.traveltimes = [999 for i in range(self.nrecs)]

            for i in range(self.nrecs):
                if self.bedvals[i] != 999 and self.dcvals[i] != 999:
                    self.traveltimes[i] = self.bedvals[i] - self.dcvals[i]
        except:
            traceback.print_exc()

    def Write(self):
        """ Write to file. """
        with open(self.fnm, 'w') as f:
            f.write('"FID","dc","bed","trav_time"\n')
            data = zip(self.dcvals, self.bedvals, self.traveltimes)

            if self.fids is None:
                sys.stderr.write("Warning: guessing FIDS in FileHandler.Write()\n")

            for i,row in enumerate(data):
                if self.fids is not None:
                    fid = self.fids[i]
                else:
                    fid = self._linloc2fid(self.line, i)
                f.write('{fid},{dcval},{bedval},{tt}\n'.format(
                    fid=fid, dcval=row[0], bedval=row[1], tt=row[2]))


def searchbylist(key, keylist, vallist, notfound=999):
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
