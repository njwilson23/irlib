#! /usr/bin/env python
#
#   icepick - graphical picking tool using matplotlib and irlib
#

import irlib
import numpy as np
import matplotlib.pyplot as plt
import sys, os
import getopt
import readline
import traceback, pdb

np.seterr(invalid='ignore')

class PickWindow:
    def __init__(self, arr=None, rate=1.e-8, bias=.02):
        self.rate = rate
        self.ntraces = 8
        self.bias = bias
        self.spacing = 0.1
        #self.spacing = 0.25
        self.intended_trace = None
        self.mode = 'bed'

        plt.ion()
        self.fig1 = plt.figure(1, figsize=(12,10))

        # Turn off default shortcuts
        for i in self.fig1.canvas.callbacks.callbacks:
            if i == 'key_press_event':
                self.fig1.canvas.mpl_disconnect(self.fig1.canvas.callbacks.callbacks[i].keys()[0])

        # Axes 1 is the trace display; axes 2 is the radargram; axes 3 permits
        # axes 2 to have yticks on the right side
        self.ax1 = self.fig1.add_axes([0.125, 0.05, 0.8, 0.6])
        self.ax1.set_autoscale_on(False)
        self.ax2 = self.fig1.add_axes([0.125, 0.65, 0.8, 0.3])
        self.lum_scale = 0.8
        self.ax3 = self.ax2.twinx()
        self.ax3.invert_yaxis()

        # Connect event handlers
        self.cid_click = self.fig1.canvas.mpl_connect( \
                            'button_press_event', self._onclick)
        self.cid_key = self.fig1.canvas.mpl_connect( \
                            'key_press_event', self._onkeypress)

        if arr is not None:
            self.Open(arr)
            self.ShowTraces()
        else:
            self.points = None
            self.T = None
            self.pick_reg = None

    def _shiftx(self, n):
        """ Shift an x-coordinate to the proper trace position on a
        multi-trace plot. """
        return n*self.spacing - self.ntraces/2*self.spacing

    def _clearlastpick(self):
        """ Remove the last pick point drawn, according to the
        registry. """
        if self.intended_trace is not None:
            self.ax1.lines.remove(self.pick_reg[self.intended_trace])
            self.pick_reg[self.intended_trace] = None
        return

    def _drawpick(self, trace, yi, intended_trace):
        """ Draw a pick mark on the given trace at yi, and update the
        registry. """
        trace = trace / trace.max() * self.spacing / 3.0
        if self.mode == 'bed':
            xshift = self._shiftx(intended_trace)
            self.ax1.plot(trace[yi] + self._shiftx(intended_trace),
                -self.T[yi], 'ob', alpha=0.4)
            self.bed_pick_reg[intended_trace] = self.ax1.lines[-1]
        elif self.mode == 'dc':
            self.ax1.plot(trace[yi] + self._shiftx(intended_trace),
                -self.T[yi], 'or', alpha=0.4)
            self.dc_pick_reg[intended_trace] = self.ax1.lines[-1]
        return self.ax1.lines

    def _clear_registers(self):
        try:
            mode = self.mode
            self.bed_pick_reg = [None for i in self.bed_pick_reg]
            self.dc_pick_reg = [None for i in self.dc_pick_reg]
            self.ChangeMode(mode)
        except AttributeError:
            pass

    def _onclick(self, event):
        """ Event handler for mouse clicks. Attempts to place a pick
        point where the mouse was clicked. """

        # Identify the trace that was aimed for
        if event.button == 2:
            self.ShowRadargram()
            return

        if event.xdata == None: return
        intended_trace = int(round(event.xdata/self.spacing + self.ntraces/2))
        try:
            trace = self.arr[:,self.cur_view[0]+intended_trace]
        except IndexError:
            return

        # Determine if trace already has a point, and if so, remove it
        if self.pick_reg[intended_trace] is not None:
            self.ax1.lines.remove(self.pick_reg[intended_trace])
            self.pick_reg[intended_trace] = None

        if event.button == 1:
            # Now put a point there
            yi = round(-event.ydata/self.rate)
            self._drawpick(trace, yi, intended_trace)

            # Save the picked point
            self.points[intended_trace+self.cur_view[0]] = yi

            self.yi = yi
            self.intended_trace = intended_trace

        else:
            self.pick_reg[intended_trace] = None
            self.points[intended_trace+self.cur_view[0]] = 999
            self.intended_trace = None

        plt.draw()

    def _onkeypress(self, event):
        """ Event handler for keystrokes. Moves pick locations or loads
        different traces.
        """

        # See about moving a placed pick
        if self.intended_trace is not None:
            try:
                trace = self.arr[:,self.cur_view[0]+self.intended_trace]
            except:
                pass

            if event.key == 'j':
                self._clearlastpick()
                self.yi += 1
                self._drawpick(trace, self.yi, self.intended_trace)
                self.points[self.intended_trace+self.cur_view[0]] = self.yi

            elif event.key == 'k':
                self._clearlastpick()
                self.yi -= 1
                self._drawpick(trace, self.yi, self.intended_trace)
                self.points[self.intended_trace+self.cur_view[0]] = self.yi

            plt.draw()

        if event.key == 'h':
            # Try moving to the left
            if self.cur_view[0] > 0:
                self.intended_trace = None
                self._clear_registers()
                self.cur_view[0] -= self.ntraces
                self.cur_view[1] -= self.ntraces
                self.ShowTraces()
                self.ShowRadargram()

        elif event.key == 'l':
            # Try moving to the right
            if self.cur_view[1] < self.arr.shape[1]:
                self.intended_trace = None
                self._clear_registers()
                self.cur_view[0] += self.ntraces
                self.cur_view[1] += self.ntraces
                self.ShowTraces()
                self.ShowRadargram()

    def Open(self, arr):
        """ Opens a line and initializes picking. """
        n = arr.shape[0]
        self.arr = arr
        self.bed_points = 999 * np.ones(arr.shape[1])
        self.dc_points = 999 * np.ones(arr.shape[1])
        self.points = self.bed_points
        self.bed_pick_reg = [None for i in range(self.ntraces)]
        self.dc_pick_reg = [None for i in range(self.ntraces)]
        self.pick_reg = self.bed_pick_reg
        self.cur_view = [0, self.ntraces]
        self.T = np.arange(0, n) * self.rate

        locs = self.ax3.get_yticks()
        self.ax3.set_yticklabels([int(round(i)) for i in locs*self.arr.shape[0]])

        self.ShowRadargram(repaint=True)
        self.ShowTraces()
        return

    def ShowTraces(self):
        """ Show traces as echograms in PickWindow axes 1. """

        self.ax1.lines = []
        self.ax1.set_autoscale_on(False)
        self.ax1.set_xlim(-self.spacing * (self.ntraces+2) / 2, self.spacing * self.ntraces / 2)
        self.ax1.cla()

        for i,j in enumerate(range(self.cur_view[0], self.cur_view[1])):
            try:
                # Plot the trace
                trace = self.arr[:,j]
                trace = trace / max(abs(trace)) * self.spacing / 3.0
                self.ax1.plot(trace + self._shiftx(i), -self.T, '-k')

                # Plot any existing picks
                oldmode = self.mode
                if self.bed_points[j] != 999:
                    self.mode = 'bed'
                    self._drawpick(trace, self.bed_points[j], i)
                if self.dc_points[j] != 999:
                    self.mode = 'dc'
                    self._drawpick(trace, self.dc_points[j], i)
                self.mode = oldmode

            except IndexError:
                pass

        try:
            self.ax2.set_title('Trace: {0}-{1} ~ Mode: {2}'
                .format(self.cur_view[0], self.cur_view[1]-1, self.mode))
            locs = self.ax1.get_yticks()
            self.ax1.set_yticklabels(locs*-1e9)
        except:
            traceback.print_exc()

        plt.draw()
        return self.cur_view

    def ShowRadargram(self, cmap='gray', c=1.68e8, repaint=False, lum_scale=None):
        """ Display a radargram on axes2. Paints in background, and
        all subsequent calls update lines. Passing repaint as True
        forces the background to be redrawn (for example, after a
        filter opperation).
        """
        n = self.arr.shape[0]
        T = np.arange(0, n*self.rate, self.rate)      # time axis

        # Clear any existing lines
        self.ax2.lines = []

        # Find the luminescense range so that plot intensity is symmetric
        if lum_scale is None:
            lum_scale = self.lum_scale
        else:
            self.lum_scale = lum_scale
        lum_bound = max((abs(self.arr.max()), abs(self.arr.min()))) * lum_scale

        # Paint on the background
        if len(self.ax2.images) == 0 or repaint is True:
            self.ax2.cla()
            self.ax2.imshow(self.arr, aspect='auto', cmap=cmap,
                            vmin=-lum_bound, vmax=lum_bound)
            locs = self.ax2.get_yticks()
            self.ax2.set_yticklabels(locs*self.rate*1e9)

        # Draw boundary brackets
        self.ax2.plot((self.cur_view[0]-1, self.cur_view[0]-1),
                        (0,self.arr.shape[0]), '-y', alpha=0.5)
        self.ax2.plot((self.cur_view[1]-1, self.cur_view[1]-1),
                        (0,self.arr.shape[0]), '-y', alpha=0.5)

        # Draw the picked event lines
        points = np.ma.array(self.bed_points)
        self.ax2.plot(np.ma.masked_where(points==999, points), '-b', alpha=0.5)

        points = np.ma.array(self.dc_points)
        self.ax2.plot(np.ma.masked_where(points==999, points), '-r', alpha=0.7)

        self.ax2.set_xlim([0, self.arr.shape[1]-1])
        self.ax2.set_ylim([self.arr.shape[0]-1, 0])
        self.ax2.set_ylabel("Time (ns)")

        self.ax2.set_title('Trace: {0}-{1} ~ Mode: {2}'
            .format(self.cur_view[0], self.cur_view[1]-1, self.mode))

        plt.draw()

        return

    def ChangeMode(self, mode):
        """ Change picking mode between bed and direct coupling. """
        self.mode = mode
        if mode == 'bed':
            self.points = self.bed_points
            self.pick_reg = self.bed_pick_reg
        elif mode == 'dc':
            self.points = self.dc_points
            self.pick_reg = self.dc_pick_reg

        self.ax2.set_title('Trace: {0}-{1} ~ Mode: {2}'
            .format(self.cur_view[0], self.cur_view[1]-1,
                    self.mode))
        plt.draw()
        return

    def Clear(self):
        """ Preserve memory by clearing lines before opening new ones. """
        self.ax1.clear()
        self.ax2.clear()


def LoadPicks(pick_window, line_gather, infile, event='bed'):
    """ Load picks from file. """
    try:
        F = irlib.FileHandler(infile, line_gather.line)
        dc_points, bed_points = F.GetEventVals()
        pick_window.bed_points = bed_points
        pick_window.dc_points = dc_points
        err = 0
    except irlib.FileHandlerError as err_message:
        print err_message
        dc_points = None
        bed_points = None
        err = 1
    return err


def OpenLine(P, infile, line, init_filters=False, fromcache=True, tocache=True):
    """ Update a PickWindow object with a new file. """
    S = irlib.Survey(infile)
    try:
        loaded = False
        cnm = S.GetLineCacheName(line)
        if fromcache:
            loaded, L = irlib.misc.TryCache(cnm)
        if loaded == False:
            L = S.ExtractLine(line)
            try:
                L.FixStaticGPS()
                L.RemoveStationary(threshold=3.0)
            except irlib.LineGatherError:
                pass
            if tocache:
                L.Dump(cnm)
        #L.Dewow()
        P.Clear()
        P.Open(L.data)
    except SystemExit:
        L = None
        P = None
    return S, L, P


def Autosave(L, P):
    """ Do an autosave to avoid 'oh shit' moments. """
    if os.path.isfile('picking/autosave.csv'):
        os.remove('picking/autosave.csv')
    try:
        L.SavePicks('picking/autosave.csv', P.bed_points, mode='bed')
        L.SavePicks('picking/autosave.csv', P.dc_points, mode='dc')
    except AttributeError:
        pass
    except irlib.FileHandlerError as s:
        print s
        print 'autosave failed'


def HandleCommand(s, infile, line, S, L, P):
    # List args. Handle empty input.
    args = s.split(' ')
    if s == '':
        return infile, line, S, L, P

    # Remove double spaces
    while 1:
        try:
            args.remove('')
        except ValueError:
            break

    if args[0] in ('exit', 'q', 'quit'):    # EXIT
        Autosave(L, P)
        sys.exit(0)

    elif args[0] == 'info':                 # INFO
        print infile
        print 'line: ' + str(line)
        try:
            print 'channel: ' + str(P.datacapture)
        except AttributeError:
            pass
        print '# traces: ' + str(P.arr.shape[1])
        print '# samples: ' + str(P.arr.shape[0])
        print 'sample interval: ' + str(P.rate) + ' s'
        print 'depth resolution: ' + str(P.rate * 1.68e8 / 2.0) + ' m'
        print 'vertical range: ' + str(1.68e8*P.arr.shape[0]*P.rate / 2.0) + ' m'
        print 'pick-mode: ' + P.mode

    elif args[0] == 'ls':                   # LS
        print [str(lnstr) for lnstr in S.GetLines()]

    elif args[0] == 'save':                 # SAVE
        outfile = 'picking/' + \
                os.path.basename(infile).split('.')[0] + \
                "_line" + str(line) + ".csv"
        try:
            L.SavePicks(outfile, P.bed_points, mode='bed')
            L.SavePicks(outfile, P.dc_points, mode='dc')
            print "Saved to {0}".format(outfile)
        except irlib.FileHandlerError as err_message:
            print err_message
        except:
            traceback.print_exc()

    elif args[0] == 'load':                 # LOAD
        # Should probably do this through irlib call like 'save' does
        loadfile = 'picking/' + \
                os.path.basename(infile).split('.')[0] + \
                "_line" + str(line) + ".csv"
        try:
            mode = P.mode
            e = LoadPicks(P, L, loadfile)   # The mode gets overridden
            P.ChangeMode(mode)              # so it needs to be reset
            P.ShowTraces()
            P.ShowRadargram()
            if e == 0:
                print "Loaded from {0}".format(loadfile)
        except:
            traceback.print_exc()

    elif args[0] == 'mode':                 # MODE
        try:
            if args[1] == 'bed':
                P.ChangeMode('bed')
                print "Picking mode set to 'bed'"
            elif args[1] == 'dc':
                P.ChangeMode('dc')
                print "Picking mode set to 'dc'"
            else:
                print "Mode option not recognized"
        except:
            print "Mode is currently '{0}'".format(P.mode)

    elif args[0] == 'open':                 # OPEN
        try:
            line = args[1]
        except:
            print "Bad arguments"
        try:
            if 'line_{0}'.format(line) in S.GetLines():
                print "Opening line {0}".format(line)
                P._clear_registers()
                S,L,P = OpenLine(P, infile, line, init_filters=False)
            else:
                print "Line {0} does not exist".format(line)
        except SystemExit:
            pass
        except:
            traceback.print_exc()

    elif args[0] == 'autobed':              # AUTOBED
        Autosave(L, P)
        lbound = None
        rbound = None
        try:
            tstart = int(args[1])
            tend = int(args[2])
            try:
                lbound = int(args[3])
                rbound = int(args[4])
            except IndexError:
                pass
        except Exception as err:
            if isinstance(err, ValueError):
                print "autobed arguments not understood", args[1:]
            elif isinstance(err, IndexError):
                pass
            else:
                print type(err), err.message
            tstart = 150
            tend = 511

        L.PickBed(sbracket=(tstart,tend), bounds=(lbound,rbound), phase=1)
        mode = P.mode
        P.bed_points = L.bed_picks  # The mode gets overridden
        P.ChangeMode(mode)          # so it needs to be reset
        P.ShowTraces()
        P.ShowRadargram()

    elif args[0] == 'autodc':               # AUTODC
        Autosave(L, P)
        try:
            tstart = int(args[1])
            tend = int(args[2])
        except Exception as err:
            if isinstance(err, ValueError):
                print "autodc arguments not understood", args[1:]
            elif isinstance(err, IndexError):
                pass
            else:
                print type(err), err.message
            tstart = 50
            tend = 150

        L.PickDC(sbracket=(tstart,tend))
        mode = P.mode
        P.dc_points = L.dc_picks    # The mode gets overridden
        P.ChangeMode(mode)          # so it needs to be reset
        P.ShowTraces()
        P.ShowRadargram()

    elif args[0] == 'clear':                # CLEAR
        Autosave(L, P)
        try:
            mode = P.mode
            if args[1] == 'bed':
                P.bed_points = 999 * np.ones(P.arr.shape[1])
                P.bed_pick_reg = [None for i in P.bed_pick_reg]
            elif args[1] == 'dc':
                P.dc_points = 999 * np.ones(P.arr.shape[1])
                P.dc_pick_reg = [None for i in P.dc_pick_reg]
            P.ChangeMode(mode)          # so it needs to be reset
            P.ShowTraces()
            P.ShowRadargram()
        except IndexError:
            pass

    elif args[0] in ('filter', 'f'):        # FILTER
        try:
            irlib.filter_defs.ApplyFilter(L, args[1:])
            P.arr = L.data
            P.ShowTraces()
            P.ShowRadargram(repaint=True)
        except IndexError:
            print L.history

    elif args[0] == 'gain':                 # GAIN
        try:
            gain = float(args[1])
            P.ShowRadargram(repaint=True, lum_scale=1.0/gain)
        except IndexError:
            sys.stdout.write('gain = ' + str(1.0 / P.lum_scale) + '\n')
        except:
            traceback.print_exc()

    elif args[0] in ('nofilter', 'nf'):     # NOFILTER
        L.Reset()
        L.Dewow()
        P.arr = L.data
        P.ShowTraces()
        P.ShowRadargram(repaint=True)

    elif args[0] == 'debug':                # DEBUG
        pdb.set_trace()

    elif args[0] == 'pick_reg':
        print P.pick_reg

    elif args[0] == 'help':                 # [empty input]
        print """
        info
        ls
        save
        load
        mode [bed|dc]
        ls
        open [line#]
        autobed [min [max]]
        autodc [min [max]]
        clear bed|dc
        filter [lowpass|highpass|lowpass_ma|highpass_ma|dewow|
                gc|agc|
                bed|englacial
                fkmig]
        nofilter
        exit
        """

    else:
        print "Command not recognized"

    return infile, line, S, L, P


# Cold start interface
def main():

    def print_syntax():
        print "\t icepick -f file_name [-L line_number]"
        return

    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], 'f:L:')
    except getopt.GetoptError:
        print "Error collecting arguments - check syntax."
        print_syntax()
        sys.exit(1)
    optdict = dict(optlist)

    try:
        infile = optdict['-f']
    except KeyError:
        print "A survey filename must be supplied:"
        print_syntax()
        sys.exit(0)

    line = int(optdict.get('-L', 0))

    # Initialize a PickWindow instance and load it with data
    P = PickWindow(rate=4e-9)
    S,L,P = OpenLine(P, infile, line, init_filters=False)

    if not os.path.isdir('picking'):
        os.mkdir('picking')

    # Begin main loop
    print "IcePick"

    while True:
        s = raw_input('>> ')
        infile,line,S,L,P = HandleCommand(s, infile, line, S, L, P)

if __name__ == '__main__':
    main()
