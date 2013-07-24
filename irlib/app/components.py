""" Define application components """

import sys
import os
import irlib
from irlib.misc import TryCache
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

class AppWindow(object):
    """ This is the generic application window class, and contains an axes
    instance and event-handling and update machinery.

    Subclasses must overload _onclick() and _onkeypress() methods.
    """

    def __init__(self, winsize):
        self.fig = plt.figure(figsize=winsize)
        self.ax = self.fig.add_subplot(1,1,1)

        # Turn off default shortcuts
        for i in self.fig.canvas.callbacks.callbacks:
            try:
                if i == 'key_press_event':
                    self.fig.canvas.mpl_disconnect(
                            self.fig.canvas.callbacks.callbacks[i].keys()[0])
            except IndexError:
                # Already disconnected
                pass

        # Connect event handlers
        self.cid_click = self.fig.canvas.mpl_connect('button_press_event', self._onclick)
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event', self._onkeypress)
        return

    def __del__(self):
        self.fig.clf()
        self.fig.canvas.mpl_disconnect(self.cid_click)
        self.fig.canvas.mpl_disconnect(self.cid_key)
        del self.fig
        return

    def _onclick(self, event):
        pass

    def _onkeypress(self, event):
        pass

    def update(self):
        """ Redraw the axes """
        self.ax.draw()
        return


class Radargram(AppWindow):
    """ Shows a radargram, and has methods for displaying annotations and
    collecting digitized features. """

    fid = 0
    active_coords = []
    features = {}
    digitize_mode = False

    lum_scale = 0.25
    cmap = plt.cm.binary

    def __init__(self, L):
        super(Radargram, self).__init__((8, 4))
        self._newline(L)
        return

    def _newline(self, L):
        """ Replace internal state with a new radar line, discarding all
        digitizing and feature information and redrawing. """
        self.rate = L.metadata.sample_rate[0]
        self.L = L
        self.data = L.data
        self.repaint()
        self.update()
        return


    def _onclick(self, event):
        """ Event handler for mouse clicks."""
        if self.digitize_mode:

            if event.button == 1:
                self.add_point(event)

            elif event.button == 2:
                self.remove_last_point()

            elif event.button == 3:
                self.add_point(event)
                self.end_feature()

        else:

            try:
                x = int(round(event.xdata))
                y = int(round(event.ydata))

                if event.button == 1:

                    xr = self.ax.get_xlim()
                    yr = self.ax.get_ylim()
                    self.ax.plot(xr, (y, y), "-k")
                    self.ax.plot((x, x), yr, "-k")
                    print x, y
                    print xr, yr
                    plt.draw()
                    #self.ax.update()



                    print ("\n\tFID: {0}\n\tx: {1}\t\ty:{2}\t\t\tt: {3} ns "
                            .format(self.L.fids[x],
                                int(round(event.xdata)), int(round(event.ydata)),
                                round(event.ydata*self.rate*1e-9,2)))

                elif event.button == 2:
                    self.repaint()
                    self.update()

                else:
                    print ("\n\teasting: {0:.1f}\tnorthing: {1:.1f}"
                            .format(self.L.metadata.eastings[x],
                                    self.L.metadata.northings[x]))

            except TypeError:
                pass
        return

    def _onkeypress(self, event):
        """ Event handler for keystrokes. """
        if event.key == 'N':
            self.add_feature('')

        elif event.key == 'E':
            self.end_feature()

    def _linloc2fid(self, loc, dc=0):
        """ Based on a line and location, return a unique FID for
        database relations. """
        eg = 0
        fid = str(self.L.line).rjust(4,'0') + str(loc).rjust(4,'0') \
            + str(dc).rjust(4,'0') + str(eg).rjust(4,'0')
        return fid

    def add_feature(self, s):
        if self.digitize_mode:
            self.end_feature()
        else:
            self.digitize_mode = True
        self.active_feature_name = s
        self.active_coords = []
        self.update()
        return self.fid

    def end_feature(self):
        self.features[self.fid] = [self.active_feature_name,
                            [xy for xy in self.active_coords]]
        self.digitize_mode = False
        self.fid += 1
        self.update()
        return self.fid-1

    def remove_feature(self, fid):
        try:
            self.features.pop(fid)
            self.active_coords = []
            self.update()
        except KeyError:
            pass

    def add_point(self, event):
        """ Record a vertex. """
        self.active_coords.append((event.xdata, event.ydata))
        self.update()

    def remove_last_point(self):
        try:
            self.active_coords.pop()
            self.update()
        except IndexError:
            pass

    def repaint(self, lum_scale=None):
        if lum_scale is None:
            lum_scale = self.lum_scale
        else:
            self.lum_scale = lum_scale

        lum_bound = max((abs(self.data.max()), abs(self.data.min()))) * lum_scale

        self.ax.cla()
        self.ax.imshow(self.data, aspect='auto', cmap=self.cmap, vmin=-lum_bound, vmax=lum_bound)
        locs = self.ax.get_yticks()
        self.ax.set_yticklabels(locs*self.rate*1e9)
        return

    def update(self, cmap='gray', c=1.68e8, repaint=False, lum_scale=None):
        """ Display a radargram on axes. Paints in background, and
        all subsequent calls update lines. Passing repaint as True
        forces the background to be redrawn (for example, after a
        filter opperation).
        """
        n = self.data.shape[0]
        self.ax.lines = []

        # Draw nodes
        drawxy = lambda xy: self.ax.plot(xy[0], xy[1], 'or', markersize=5.0,
                                          markeredgewidth=0.0, alpha=0.5)
        points_ = map(drawxy, self.active_coords)

        # Draw previous features
        if len(self.features) > 0:

            drawline = lambda lsxy: self.ax.plot(
                [i[0] for i in lsxy[1]], [i[1] for i in lsxy[1]],
                '--r')
            lines_ = map(drawline, self.features.values())

            labelfeature = lambda key, lsxy: self.ax.text(lsxy[1][-1][0],
                lsxy[1][-1][1]-20, str(key), fontsize=12, color='r')
            text_ = map(labelfeature, self.features.keys(), self.features.values())

        # Force tight bounding
        self.ax.set_xlim([0, self.data.shape[1]-1])
        self.ax.set_ylim([self.data.shape[0]-1, 0])

        # Decorate and draw
        self.ax.set_ylabel("Time (ns)")
        self.ax.set_xlabel("Location number")
        if self.digitize_mode:
            self.ax.set_title("Line {0} [feature {1}]".format(self.line, self.fid))
        else:
            self.ax.set_title("Line {0} [viewing]".format(self.L.line))
        plt.draw()
        return

    def get_digitizer_filename(self):
        fnm = os.path.join("englacial",
            os.path.basename(self.L.infile).split(".")[0] + "_line" + str(self.L.line) + ".txt")
        return fnm

    def load(self, f):
        """ Parse a digitizer file and return a dictionary with list
        entries that can be dropped directly into an ImageWindow.
        """
        self.digitize_mode = False
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
        self.update()

        return

    def save(self):
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
            h5addrs = ['line_{0}/location_{1}/datacapture_0/echogram_0'.format(self.line, loc) 
                        for loc in locs]
            map(lambda h5addr: RL.AddDataset(self.fh5[h5addr]), h5addrs)
            lons = [-lon for lon in RL.lons]
            lats = RL.lats

            # Combine it into a dictionary
            feature_dict = dict(zip(fids, zip(lons, lats, depths)))
            feature_dict['fnum'] = fnum
            dict_list.append(feature_dict)

        return dict_list

class PickWindow(AppWindow):
    pass

class MapWindow(AppWindow):
    pass

