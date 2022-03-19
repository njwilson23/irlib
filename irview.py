#! /usr/bin/env python
#
#   irview - simple radar line viewer and feature digitizer
#

import irlib
from irlib.misc import TryCache
from irlib import app
import numpy as np
import matplotlib.pyplot as plt
import sys, os
import getopt
try:
  import readline
except ImportError:
  import pyreadline as readline
import traceback, pdb

np.seterr(invalid='ignore')

class ImageWindow:
    def __init__(self, L, fh5=None, rate=1.e-8):
        try:
            # Assume sample rate is const.
            self.rate = 1./L.metadata.sample_rate[0]
        except:
            traceback.print_exc()
            sys.stderr.write('something went wrong reading digitizer rate\n')
            self.rate = rate
        self.L = L
        self.arr = L.data
        self.line = L.line
        self.datafile = L.infile
        self.fh5 = fh5
        self.interval = 5
        self.radar_fids = [L.GetFID(i) for i in range(L.data.shape[1])]

        # Feature digitizing
        self.fid = 0
        self.active_coords = []
        self.features = {}
        self.digitize = False

        nx = self.arr.shape[1]
        ny = self.arr.shape[0]
        self.T = np.arange(0, ny*self.rate, self.rate)

        # Set up the plotting window
        plt.ion()
        self.fig1 = plt.figure(1, figsize=(12,5))

        # Axes 1 is the radargram
        self.ax1 = self.fig1.add_axes([0.1, 0.1, 0.85, 0.8])
        self.lum_scale = 0.25

        # Turn off default shortcuts

        key_press_cids = self.fig1.canvas.callbacks.callbacks.get('key_press_event', {}).copy()
        for cid in key_press_cids.keys():

            self.fig1.canvas.mpl_disconnect(cid)

        # Connect event handlers
        self.cid_click = self.fig1.canvas.mpl_connect( \
                            'button_press_event', self._onclick)
        self.cid_key = self.fig1.canvas.mpl_connect( \
                            'key_press_event', self._onkeypress)

        self.ShowRadargram()

    def _onclick(self, event):
        """ Event handler for mouse clicks."""
        if self.digitize:

            if event.button == 1:
                self.AddPoint(event)

            elif event.button == 2:
                self.RemoveLastPoint()

            elif event.button == 3:
                self.AddPoint(event)
                self.EndFeature()

        else:

            try:
                x = int(round(event.xdata))

                if event.button == 1:
                    sys.stdout.write(
                        "\n\tFID: {0}\n\tx: {1}\t\ty:{2}\t\t\tt: {3} ns\n>> "
                        .format(self.radar_fids[x],
                                int(round(event.xdata)),
                                int(round(event.ydata)),
                                round(event.ydata*self.rate*1e9,2)))

                elif event.button == 2:
                    self.ShowRadargram(repaint=True)

                else:
                    sys.stdout.write(
                        "\n\teasting: {0:.1f}\tnorthing: {1:.1f}\n>> ".format(
                            self.L.metadata.eastings[x],
                            self.L.metadata.northings[x]))

            except TypeError:
                pass

    def _onkeypress(self, event):
        """ Event handler for keystrokes. """
        if event.key == 'N':
            self.AddFeature('')

        elif event.key == 'E':
            self.EndFeature()

    def _linloc2fid(self, loc):
        """ Based on a line and location, return a unique FID for
        database relations. """
        dc = 0
        eg = 0
        fid = str(self.line).rjust(4,'0') + str(loc).rjust(4,'0') \
            + str(dc).rjust(4,'0') + str(eg).rjust(4,'0')
        return fid

    def AddFeature(self, s):
        if self.digitize:
            self.EndFeature()
        else:
            self.digitize = True
        self.active_feature_name = s
        self.active_coords = []
        self.ShowRadargram()
        return self.fid

    def EndFeature(self):
        self.features[self.fid] = [self.active_feature_name,
                            [xy for xy in self.active_coords]]
        self.digitize = False
        self.fid += 1
        self.ShowRadargram()
        return self.fid-1

    def RemoveFeature(self, fid):
        try:
            self.features.pop(fid)
            self.active_coords = []
            self.ShowRadargram(repaint=True)
        except KeyError:
            pass

    def AddPoint(self, event):
        """ Record a vertex. """
        self.active_coords.append((event.xdata, event.ydata))
        self.ShowRadargram()

    def RemoveLastPoint(self):
        try:
            self.active_coords.pop()
            self.ShowRadargram()
        except IndexError:
            pass

    def ShowRadargram(self, cmap='gray', c=1.68e8, repaint=False, lum_scale=None):
        """ Display a radargram on axes. Paints in background, and
        all subsequent calls update lines. Passing repaint as True
        forces the background to be redrawn (for example, after a
        filter opperation).
        """
        n = self.arr.shape[0]
        self.ax1.lines = []

        # Find the luminescense range so that plot intensity is symmetric
        if lum_scale is None:
            lum_scale = self.lum_scale
        else:
            self.lum_scale = lum_scale
        lum_bound = max((abs(self.arr.max()), abs(self.arr.min()))) * lum_scale

        # Paint on the background
        if len(self.ax1.images) == 0 or repaint is True:
            self.ax1.cla()
            self.ax1.imshow(self.arr, aspect='auto', cmap=cmap,
                            vmin=-lum_bound, vmax=lum_bound)
            locs = self.ax1.get_yticks()
            self.ax1.set_yticklabels(locs*self.rate*1e9)

        # Draw nodes
        drawxy = lambda xy: self.ax1.plot(xy[0], xy[1], 'or', markersize=5.0,
                                          markeredgewidth=0.0, alpha=0.5)
        points_ = map(drawxy, self.active_coords)

        # Draw previous features
        if len(self.features) > 0:

            drawline = lambda lsxy: self.ax1.plot(
                [i[0] for i in lsxy[1]], [i[1] for i in lsxy[1]],
                '--r')
            lines_ = map(drawline, self.features.values())

            labelfeature = lambda key, lsxy: self.ax1.text(lsxy[1][-1][0],
                lsxy[1][-1][1]-20, str(key), fontsize=12, color='r')
            text_ = map(labelfeature, self.features.keys(), self.features.values())

        # Force tight bounding
        self.ax1.set_xlim([0, self.arr.shape[1]-1])
        self.ax1.set_ylim([self.arr.shape[0]-1, 0])

        # Decorate and draw
        self.ax1.set_ylabel("Time (ns)")
        self.ax1.set_xlabel("Location number")
        if self.digitize:
            self.ax1.set_title("Line {0} [feature {1}]".format(self.line, self.fid))
        else:
            self.ax1.set_title("Line {0} [viewing]".format(self.line))
        plt.draw()
        return

    def GetDigitizerFilename(self):
        fnm = 'englacial/' + os.path.basename(self.datafile).split('.')[0] + \
            "_line" + str(self.line) + ".txt"
        return fnm

    def Import(self, f):
        """ Parse a digitizer file and return a dictionary with list
        entries that can be dropped directly into an ImageWindow.
        """
        self.digitize = False
        features = {}
        i = 0

        while True:
            # Read a feature
            pnt_list = []

            while True:
                # Read a point
                s = f.readline()
                if s in ('\n', ''):
                    break
                else:
                    try:
                        slist = s.split()
                        fid = slist[0]
                        try:
                            x = float(self.radar_fids.index(fid))
                        except:
                            x = float(self.radar_fids.index(fid + 4*'0'))
                        y = float(slist[3])
                        pnt_list.append((x, y))
                    except:
                        traceback.print_exc()
                        sys.stdout.write("Failed to read record:\n\t{0}".format(s))

            if len(pnt_list) == 0:
                break
            else:
                features[i] = ['', pnt_list]
                i += 1

        self.features = features
        self.fid = i
        self.ShowRadargram()

        return

    def Export(self):
        """ Returns a list of dictionaries containing the latitude, longitude,
        and the y-axis value of each vertex. Dictionary keys are standard
        linloc FIDs.
        """

        if self.datafile is None: return 1

        dict_list = []

        for fnum in self.features:
            coords = self.features[fnum]

            # Look up the data to be exported
            locs = [int(round(xy[0])) for xy in coords[1]]
            depths = [xy[1] for xy in coords[1]]
            fids = [self.radar_fids[loc] for loc in locs]
            RL = irlib.RecordList(self.datafile)
            h5addrs = ['line_{0}/location_{1}/datacapture_0/echogram_0'.format(self.line, loc) for loc in locs]
            map(lambda h5addr: RL.AddDataset(self.fh5[h5addr]), h5addrs)
            lons = [-lon for lon in RL.lons]
            lats = RL.lats

            # Combine it into a dictionary
            feature_dict = dict(zip(fids, zip(lons, lats, depths)))
            feature_dict['fnum'] = fnum
            dict_list.append(feature_dict)

        return dict_list

    def Close(self):
        self.fig1.clf()
        self.fig1.canvas.mpl_disconnect(self.cid_click)
        self.fig1.canvas.mpl_disconnect(self.cid_key)


def OpenLine(S, line, fh5, fromcache=True, tocache=False, datacapture=0):
    """ Start up an AppWindow object. Optionally attempts to read Gather
    instance from a cache. If not, or if this fails, it can write a cached
    version for future use.

            S (irlib.Survey) is a survey instance
            line (int) is the number of the line to open
            fh5 (file-like) is a file object
    """
    try:
        loaded = False
        cnm = S.GetLineCacheName(line, dc=datacapture)
        if fromcache:
            loaded, L = TryCache(cnm)
        if loaded == False:
            assert datacapture < S.GetChannelsInLine(int(line))

            L = S.ExtractLine(line, datacapture=datacapture)
            try:
                L.RemoveBadLocations()
                L.FixStaticGPS()
                L.RemoveBlankTraces()
                L.SmoothenGPS()
                L.RemoveStationary(threshold=3.0)
            except irlib.LineGatherError:
                pass
            if tocache:
                L.Dump(cnm)
        IW = ImageWindow(L, fh5=fh5)
    except (IndexError, SystemExit):
        L = None
        IW = None
    except AssertionError:
        print ("no channel {0} found".format(datacapture))
        L = None
        IW = None
    return IW, L


def StrFilterHistory(L):
    """ Return a printable string summarizing the filter history. """
    s_out = ""
    for operation in L.history:
        if hasattr(operation, "__iter__"):
            s_out += "\t" + operation[0] + " [ "
            for option in operation[1:]:
                s_out += str(option) + ", "
            s_out += " ] \n"
        else:
            s_out += "\t" + operation + " [ ] \n"
    return s_out

def WriteFeatures(IW, outfnm):
    """ Write digitized features from `IW` to `outfnm`, creating folder if
    necessary. """
    dict_list = IW.Export()

    if not os.path.isdir(os.path.dirname(outfnm)):
        os.makedirs(os.path.dirname(outfnm))

    with open(outfnm, 'w') as fout:
        for fdict in dict_list:
            keys = fdict.keys()
            keys.sort()
            for key in keys:
                vals = fdict[key]
                if isinstance(vals, __builtins__.tuple):
                    fout.write('{0}\t{1}\t{2}\t{3}\n'.format(
                        key, vals[0], vals[1], vals[2]))
            fout.write('\n')
    return

def HandleCommand(s, S, IW, L):

    # List args. Handle empty input.
    args = s.split(' ')
    if s == '':
        return IW, L

    # Remove double spaces
    while 1:
        try:
            args.remove('')
        except ValueError:
            break

    if args[0] in ('exit', 'q', 'quit'):    # EXIT
        sys.exit(0)

    elif args[0] == 'info':                 # INFO
        print (S.datafile)
        print ('line: ' + str(IW.line))
        try:
            print ('channel: ' + str(L.datacapture))
        except AttributeError:
            pass
        print ('# traces: ' + str(IW.arr.shape[1]))
        print ('# samples: ' + str(IW.arr.shape[0]))
        print ('sample interval: ' + str(IW.rate) + ' s')
        print ('depth resolution: ' + str(IW.rate * 1.68e8 / 2.0) + ' m')
        print ('vertical range: ' + str(1.68e8*IW.arr.shape[0]*IW.rate / 2.0) + ' m')
        print ('available channels: ' + str(S.GetChannelsInLine(int(L.line))))

    elif args[0] == 'ls':                   # LS
        print (reduce(lambda a,b: a + "{0:>9}".format(b), S.GetLines()))
        print (9*" "*IW.line + 4*" " + "^")

    elif args[0] == 'open':                 # OPEN
        try:
            line = args[1]
        except IndexError:
            print ("Must supply at least a line number")
            return IW, L

        try:
            datacapture = int(args[2])
        except IndexError:
            datacapture = 0
        except ValueError:
            print ("Bad channel number: {0}".format(args[2]))
            return IW, L

        try:
            if 'line_{0}'.format(line) in S.GetLines() and \
                datacapture < S.GetChannelsInLine(int(line)):
                print ("Opening line {0}, channel {1}".format(line, datacapture))
                if IW:
                    IW.Close()
                    del IW
                    del L
                IW, L = OpenLine(S, line, S.f, datacapture=datacapture)
            else:
                print ("Line {0} channel {1} does "
                       "not exist".format(line, datacapture))
        except:
            traceback.print_exc()

    elif args[0] in ('filter', 'f'):        # FILTER
        try:
            #irlib.ApplyFilter(L, args[1:])
            app.command_parser.apply_filter(args[1:], L)
            IW.arr = L.data
            IW.ShowRadargram(repaint=True)
        except app.command_parser.CommandSearchError as e:
            print (e.message)
        except IndexError:
            print (StrFilterHistory(L))

    elif args[0] in ('nofilter', 'nf'):     # NOFILTER
        L.Reset()
        IW.arr = L.data
        IW.ShowRadargram(repaint=True)

    elif args[0] == 'hist':                 # HIST
        plt.figure(2)
        plt.clf()
        plt.hist(L.data.flatten(), 100, cumulative=False, normed=True, histtype='step')
        plt.title('Diagnostic histogram')
        plt.figure(1)

    elif args[0] == 'dnew':                 # DNEW
        if len(args) == 1:
            args.append('')
        try:
            IW.AddFeature(args[1])
        except:
            traceback.print_exc()

    elif args[0] == 'drm':                  # DRM
        try:
            IW.RemoveFeature(int(args[1]))
        except IndexError:
            print ("No such feature")

    elif args[0] == 'dls':                  # DLS
        print ("Feature List")
        for key in IW.features.keys():
            print ("{0}: {1} vertices".format(key, len(IW.features[key][1])))

    elif args[0] == 'dsave':                # DSAVE
        try:
            dict_list = IW.Export()
            outfnm = IW.GetDigitizerFilename()
            WriteFeatures(IW, outfnm)
            print ("Features saved to " + outfnm)
        except:
            traceback.print_exc()

    elif args[0] == 'dimport':              # DIMPORT
        infnm = IW.GetDigitizerFilename()
        try:
            with open(infnm, 'r') as fin:
                IW.Import(fin)
        except IOError:
            sys.stdout.write("Digitizer file ({0}) not found\n".format(infnm))
        except:
            traceback.print_exc()

    elif args[0] == 'imsave':               # IMSAVE
        try:
            fnm = args[1]
            if os.path.splitext(fnm)[1] != '.png':
                fnm += '.png'
            IW.fig1.savefig(fnm)
        except IndexError:
            sys.stdout.write("Please supply a save location/name\n")
        except:
            traceback.print_exc()

    elif args[0] == 'gain':                 # GAIN
        try:
            gain = float(args[1])
            IW.ShowRadargram(repaint=True, lum_scale=1.0/gain)
        except IndexError:
            sys.stdout.write('gain = ' + str(1.0 / IW.lum_scale) + '\n')
        except:
            traceback.print_exc()

    elif args[0] == 'debug':                # DEBUG
        pdb.set_trace()

    elif args[0] == 'help':                 # HELP
        if len(args) == 1:
            print ("""\tApplication commands:

            info                print line metadata
            ls                  list lines in survey
            open [line#]        open a different line
            gain [#]            adjust display contrast

            filter *name*       apply a filter (see below)
            nofilter            remove all filters

            dnew                start digitizing new feature
            drm [#]             remove a feature
            dls                 list features
            dsave               export features to text

            imsave [file]       save the radargram as an image
            hist                show a brightness histogram
            exit                exit irview
            debug
            """)

            print ("\tAvailable filter commands:\n")
            for name in app.command_parser.list_filters():
                print ("\t\t{0}".format(name))

        else:
            print (app.command_parser.help_filter(args[1]))

    else:
        print ("Command not recognized")

    return IW, L


# Cold start
def main():

    def print_syntax():
        print ("\t irview -f file_name [-L line_number]")
        return

    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], 'f:L:')
    except getopt.GetoptError:
        print ("Error collecting arguments - check syntax.")
        print_syntax()
        sys.exit(1)
    optdict = dict(optlist)

    try:
        infile = optdict['-f']
    except KeyError:
        print ("A survey filename must be supplied:")
        print_syntax()
        sys.exit(0)

    line = int(optdict.get('-L', 0))

    S = irlib.Survey(infile)
    fh5 = S.f
    IW,L = OpenLine(S, line, fh5)

    # Begin main loop
    print ("Ice Radar Viewing Tool (IRView)")
    print ("Press [shift+n] to digitize a feature.")
    print ("Right-click or press [shift+e] to end the feature.")
    print ("Type 'help' for additional commands.")
    while 1:
        s = raw_input('>> ')
        IW, L = HandleCommand(s, S, IW, L)

if __name__ == '__main__':
    main()
