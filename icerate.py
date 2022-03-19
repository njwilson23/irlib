#! /usr/bin/env python
'''
   IceRate - Graphical trace quality rating tool using matplotlib and irlib


TODO: replace getopt with argparse

'''

import irlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import sys, os, getopt
import traceback, pdb
try:
  import readline
except ImportError:
  import pyreadline as readline
from random import shuffle

np.seterr(invalid='ignore')

class RatingWindow(object):
    def __init__(self, L, picks, ratings=None, rate=1.e-8, bias=.02):
        self.isopen = True
        self.bias = bias
        self.rate = rate
        self.picks = picks
        self.arr = L.data
        self.line = L.line
        self.spacing = 0.25
        self.interval = 5
        if ratings is None:
            self.ratings = -9 * np.ones((len(picks)))
        elif len(ratings) == len(picks):
            self.ratings = ratings
        else:
            self.ratings = -9 * np.ones((len(picks)))
            sys.stderr.write("ratings passed to RatingWindow have inconsistent length\n")

        ny = self.arr.shape[0]
        # made list of this range since it could not be sorted otherwise
        self.traces = list(range(0, self.arr.shape[1], self.interval))  
        #shuffle(self.traces)
        self.n = 0
        self.cur_trace = self.traces[self.n]

        self.T = np.arange(0, ny) * self.rate

        # Set up the plotting window
        plt.ion()
        self.fig1 = plt.figure(1, figsize=(12,5))
     
        self.fig1.canvas.manager.set_window_title("Pick rater")

        # Axes 1 is the radargram; axes 2 is the trace being rated
        self.ax1 = self.fig1.add_axes([0.1, 0.05, 0.7, 0.9])
        self.ax2 = self.fig1.add_axes([0.8, 0.05, 0.1, 0.9])
        self.ax2.set_xticklabels('')
        self.ax2.set_yticklabels('')
        
        # Turn off default shortcuts

        key_press_cids = self.fig1.canvas.callbacks.callbacks.get('key_press_event', {}).copy()
        for cid in key_press_cids.keys():

            self.fig1.canvas.mpl_disconnect(cid)

        # Connect event handlers
        #self.cid_click = self.fig1.canvas.mpl_connect( \
                            #'button_press_event', self._onclick)
        self.cid_key = self.fig1.canvas.mpl_connect( \
                            'key_press_event', self._onkeypress)
        self.fig1.canvas.mpl_connect('close_event', self._onclose)

        self.ShowRadargram()
        self.ShowTraces()

    def _drawpick(self, trace, yi):
        """ Draw a pick mark on the given trace at yi. """
        
        trace = trace / trace.max() * self.spacing / 3.0
        self.ax2.plot(trace[int(yi)], -self.T[int(yi)], 'ob', alpha=0.4)
        return self.ax2.lines

    def _onclick(self, event):
        """ Event handler for mouse clicks."""
        pass

    def _onclose(self, event):
        print("exit\n")
        self.isopen = False

    def _rate(self, rating):
        """ Assign a rating to a specific event pick. """
        self.ratings[self.cur_trace:self.cur_trace+self.interval] = rating

    def _next_trace(self):
        """ Iterate to the next trace and reload the radargram and trace
        window. """
        self.n += 1
        self.n = min(self.n, len(self.traces)-1)
        self.cur_trace = self.traces[self.n]
        self.ShowRadargram()
        self.ShowTraces()


    def _prev_trace(self):
        """ Jump to the previous trace. """
        self.n -= 1
        self.n = max(self.n, 0)
        self.cur_trace = self.traces[self.n]
        self.ShowRadargram()
        self.ShowTraces()

    def _onkeypress(self, event):
        """ Event handler for keystrokes. """
        # Record the rating, or handle a command
        try:
            if int(event.key) >= 1 and int(event.key) <= 5:
                self._rate(int(event.key))
                self._next_trace()
        except ValueError:
            if event.key == 'n':
                self._next_trace()
            elif event.key == 'p':
                self._prev_trace()

    def ShowTraces(self, view=(None, None)):
        """ Show traces as echograms in PickerWindow axes 1. """

        self.ax2.clear()

        # Plot the trace
        trace = self.arr[:,self.cur_trace]
        trace = trace / trace.max() * self.spacing / 3.0
        self.ax2.plot(trace, -self.T, '-k')
        
        # Plot the picked time
        if self.picks[self.cur_trace] != -999:  #DM changed from 999 to -999
            self.mode = 'bed'
            if not np.isnan(self.picks[self.cur_trace]):
                self._drawpick(trace, self.picks[self.cur_trace])
        else:
            # should skip automatically
            pass

        plt.draw()

        return self.ax2.lines

    def ShowRadargram(self, cmap='gray', rate=1e-8, c=1.68e8, repaint=False):
        """ Display a radargram on axes2. Paints in background, and
        all subsequent calls update lines. Passing repaint as True
        forces the background to be redrawn (for example, after a
        filter opperation).
        """
        #self.ax1.lines = []

        # Find the luminescense range so that plot intensity is symmetric
        lum_bound = max((abs(self.arr.max()), abs(self.arr.min())))

        # Paint on the background
        if len(self.ax1.images) == 0 or repaint is True:
            self.ax1.imshow(self.arr, aspect='auto', cmap=cmap, vmin=-lum_bound, vmax=lum_bound)
            locs = self.ax1.get_yticks().tolist()
            self.ax1.yaxis.set_major_locator(mticker.FixedLocator(locs))
            self.ax1.set_yticklabels(["{:0.0f}".format((x*10)) for x in locs])

        # Draw location indicator
        self.ax1.plot((self.cur_trace, self.cur_trace),
                        (0,self.arr.shape[0]), '-y', alpha=0.6, linewidth=2.)
        self.ax1.plot((self.cur_trace+self.interval, self.cur_trace+self.interval),
                        (0,self.arr.shape[0]), '-y', alpha=0.6, linewidth=2.)

        # Draw the picked event region
        try:
            xa = self.cur_trace
            xb = min(self.cur_trace + self.interval + 1, self.arr.shape[1])
            upper_line = map(lambda n: (n==1004 and np.nan or n),
                            [i+5 for i in self.picks[xa:xb]])
            lower_line = map(lambda n: (n==994 and np.nan or n),
                            [i-5 for i in self.picks[xa:xb]])
            self.ax1.plot(list(range(xa, xb)), list(upper_line), '-y', alpha=0.6, linewidth=2.)
            self.ax1.plot(list(range(xa, xb)), list(lower_line), '-y', alpha=0.6, linewidth=2.)
        except:
            traceback.print_exc()

        self.ax1.set_xlim([0, self.arr.shape[1]-1])
        self.ax1.set_ylim([self.arr.shape[0]-1, 0])
        self.ax1.set_ylabel("Time (ns)")

        plt.draw()

        return


def LoadRatings(infile):
    """ Load ratings from file. """
    try:
        with open(infile, 'r') as f:
            ratings = f.readlines()
        valid = True
        err = 0

        for line in ratings:
            try:
                int(line.split('\t')[0])
            except ValueError:
                valid = False
                err += 1
            except IndexError:
                valid = False
                err += 1

        if valid:
            ratings = [int(float(line.split('\t')[1].rstrip('\n'))) for line in ratings]
        else:
            ratings = None

    except:
        traceback.print_exc()
    return ratings, err


def OpenLine(infile, line, pickfile, fromcache=True):
    """ Start up a RatingWindow object. """
    print(pickfile)
    if os.path.isfile(pickfile) == False :
        sys.stderr.write("\n\tNo pick file found for line " + str(line) + ".\n\tPress 'up' key and change line number.\n\n")
        exit(1)
    S = irlib.Survey(infile)
    
    L = S.ExtractLine(line, fromcache=fromcache)
    if not fromcache:
        try:
            L.FixStaticGPS()
            L.RemoveStationary(threshold=3.0)
        except irlib.LineGatherError:
            pass
    try:
        F = irlib.FileHandler(pickfile, L.line)
        dc_points, bed_points = F.GetEventVals()
        R = RatingWindow(L, bed_points)

    except irlib.FileHandlerError as err_message:
        print(err_message)
        dc_points = None
        bed_points = None
        err = 1
        R = None

    return R, L, S


def Autosave(L, R):
    """ Do an autosave to avoid 'oh shit' moments. """
    print("Temporary backup saved to rating/autosave.txt")
    if os.path.isfile('rating/autosave.txt'):
        os.remove('rating/autosave.txt')
    SaveRatings('rating/autosave.txt', R, L)


def SaveRatings(outfile, R, L):
    """ Save ratings to file, with FID and rating. """
    with open(outfile, 'w') as f:
        for fid, rating in zip(L.metadata.fids, list(R.ratings)):
            f.write(fid + "\t" + str(int(rating)) + "\n")
    return


def linloc2fid(lin, loc):
    """ Based on a line and location, return a unique FID for database
    relations. """
    try:
        dc = 0
        eg = 0
        fid = str(lin).rjust(4,'0') + str(loc).rjust(4,'0') \
            + str(dc).rjust(4,'0') + str(eg).rjust(4,'0')
        return fid
    except:
        traceback.print_exc()
        sys.stderr.write('icerate.py: failed at linloc2fid\n')
        return None


def HandleCommand(s, infile, R, L, S):

    # List args. Handle empty input.
    args = s.split(' ')
    if s == '':
        return R, L

    # Remove double spaces
    while 1:
        try:
            args.remove('')
        except ValueError:
            break

    if args[0] in ('exit', 'q', 'quit'):    # EXIT
        Autosave(L, R)
        sys.exit(0)

    elif args[0] == 'info':                 # INFO
        print(infile)
        print('line: ' + str(R.line))
        print('location: ' + str(R.cur_trace))
        print('nx: ' + str(R.arr.shape[1]))
        print('nz: ' + str(R.arr.shape[0]))
        print('rated: ' + str(sum([1 for r in R.ratings if r != -9])))
        print('unrated: ' + str(sum([1 for r in R.ratings if r == -9])))

    elif args[0] == 'ls':                   # LS
        print(S.GetLines())

    elif args[0] == 'save':                 # SAVE
        outfile = 'rating/' + \
                os.path.basename(infile).split('.')[0] + \
                "_line" + str(R.line) + ".txt"
        try:
            SaveRatings(outfile, R, L)
            print("Saved to {0}".format(outfile))
        except:
            traceback.print_exc()

    elif args[0] == 'load':                 # LOAD
        # Should probably do this through irlib call like 'save'
        loadfile = 'rating/' + \
                os.path.basename(infile).split('.')[0] + \
                "_line" + str(R.line) + ".txt"
        try:
            ratings, e = LoadRatings(loadfile)
            if e == 0:
                R.ratings = ratings
                R.ShowTraces()
                R.ShowRadargram()
                print("Loaded from {0}".format(loadfile))
        except:
            traceback.print_exc()

    elif args[0] == 'open':                 # OPEN
        if len(args) < 2:
            print("Not enough arguments")
        line = args[1]

        if 'line_{0}'.format(line) in S.GetLines():
            print("Opening line {0}".format(line))

            if len(args) >= 3:
                pickfile = args[2]
            else:
                pickfile = None

            if pickfile and not os.path.exists(pickfile):
                print("Pick file {0} does not exist".format(pickfile))
                pickfile = None

            if pickfile is None:
                pickfile = get_pickfnm(infile, line)

            print(pickfile)

            R.fig1.clf()
            del R
            R,L,S = OpenLine(infile, line, pickfile)
        else:
            print("Line {0} does not exist".format(line))

    elif args[0] in ('filter', 'f'):        # FILTER
        try:
            irlib.filter_defs.ApplyFilter(L, args[1:])
            R.arr = L.data
            R.ShowRadargram(repaint=True)
        except IndexError:
            print(L.history)

    elif args[0] in ('nofilter', 'nf'):     # NOFILTER
        L.Reset()
        R.arr = L.data
        R.ShowTraces()
        R.ShowRadargram(repaint=True)

    elif args[0] == 'order':                # ORDER
        print('location: ' + str(R.cur_trace))
        for n in R.traces:
            sys.stdout.write(str(n) + '\t')
        sys.stdout.write('\n')

    elif args[0] == 'sort':                  # SORT
        R.traces.sort()
        R.cur_trace = R.traces[0]
        R.ShowRadargram()
        R.ShowTraces()

    elif args[0] == 'randomize':            # RANDOMIZE
        shuffle(R.traces)
        R.ShowRadargram()
        R.ShowTraces()

    elif args[0] == 'ratings':              # RATINGS
        for n in R.ratings:
            sys.stdout.write(str(n) + '\t')
        sys.stdout.write('\n')

    elif args[0] == 'debug':
        pdb.set_trace()

    elif args[0] == 'help':                 # [empty input]
        print("""
        info
        save
        load
        open [line#]
        filter [lowpass|highpass|lowpass_ma|highpass_ma|dewow|
                gc|agc|
                bed|englacial
                fkmig]
        nofilter
        order
        sort
        randomize
        ratings
        debug
        exit
        """)

    else:
        print("Command not recognized")

    return R, L

def get_pickfnm(infile, line):
    return os.path.join("picking",
            os.path.basename(infile).split(".")[0]
                + "_line" + str(line) + ".csv")

# Cold start interface
def main():

    def print_syntax():
        print("\t icerate -f file_name [-L line_number] [--pick pick_filename]")
        return

    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], 'f:L:', ['pick='])
    except getopt.GetoptError:
        print("Error collecting arguments - check syntax.")
        print_syntax()
        sys.exit(1)
    optdict = dict(optlist)

    try:
        infile = optdict["-f"]
    except IndexError:
        print("A survey filename must be supplied:")
        print_syntax()
        sys.exit(0)

    line = int(optdict.get('-L', 0))
    pickfile = optdict.get('--pick', get_pickfnm(infile, line))
    R,L,S = OpenLine(infile, line, pickfile)

    if not os.path.isdir('rating'):
        os.mkdir('rating')

    # Begin main loop
    print("IceRate")

    while R.isopen:
        if sys.version_info[0]  == 2:
            s = raw_input('>> ')
        else:
            s = input(">> ")
        R, L = HandleCommand(s, infile, R, L, S)

if __name__ == '__main__':
    main()
