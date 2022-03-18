""" Define application components """

import sys
import os
import irlib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from irlib.misc import TryCache

plt.ion()
matplotlib.rcParams["toolbar"] = "none"
matplotlib.rcParams["keymap.home"] = ""
matplotlib.rcParams["keymap.back"] = ""
matplotlib.rcParams["keymap.pan"] = ""
matplotlib.rcParams["keymap.zoom"] = ""
matplotlib.rcParams["keymap.forward"] = ""
matplotlib.rcParams["keymap.xscale"] = ""
matplotlib.rcParams["keymap.yscale"] = ""
matplotlib.rcParams["keymap.quit"] = ""

class AppWindow(object):
    """ This is the generic application window class, and contains an axes
    instance and event-handling and update machinery.

    Subclasses must overload _onclick() and _onkeypress() methods.
    """

    def __init__(self, winsize):
        """Parameters
        ----------
        winsize : size of the window, in inches
        """
        self.fig = plt.figure(figsize=winsize)
        self.ax = self.fig.add_subplot(1,1,1)

        # Turn off default shortcuts

        key_press_cids = self.fig.canvas.callbacks.callbacks.get('key_press_event', {}).copy()
        for cid in key_press_cids.keys():
            self.fig.canvas.mpl_disconnect(cid)

        # Connect event handlers
        self.cid_click = self.fig.canvas.mpl_connect('button_press_event', self._onclick)
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event', self._onkeypress)
        return

    def __del__(self):
        #self.fig.clf()
        self.fig.canvas.mpl_disconnect(self.cid_click)
        self.fig.canvas.mpl_disconnect(self.cid_key)
        plt.close(self.fig)
        del self.fig
        return

    def _onclick(self, event):
        pass

    def _onkeypress(self, event):
        pass

    def _newline(self, line):
        """ This method must be overloaded to define how the AppWindow reacts
        when the line changes. """
        return

    def update(self):
        """ Redraw the axes """
        self.fig.canvas.draw()
        return


class Radargram(AppWindow):
    """ Shows a radargram, and has methods for displaying annotations and
    collecting digitized features. """

    fid = 0
    active_coords = []
    annotations = {}
    #features = {}

    lum_scale = 0.25
    cmap = plt.cm.binary

    bbox = [None, None, None, None]

    def __init__(self, L):
        """Parameters
        ----------
        L : Gather instance containing the radar data to be displayed
        """
        super(Radargram, self).__init__((10, 4))
        self.digitize_mode = False
        self._newline(L)
        self.fig.canvas.manager.set_window_title("Radargram") 
        self.fig.tight_layout()
        return

    def _newline(self, L):
        """ Replace internal state with a new radar line, discarding all
        digitizing and feature information and redrawing. """
        self.rate = L.metadata.sample_rate[0]
        self.L = L
        self.data = L.data
        self.annotations = {}
        self.repaint()
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

                    #self.remove_annotation("x-hair")
                    #self.remove_annotation("x-hair-text")
                    if "x-hair" in self.annotations:
                        anns = self.annotations.get("x-hair", []) + \
                                self.annotations.get("x-hair-text", [])
                        for item in anns:
                            item.remove()

                    xr = self.ax.get_xlim()
                    yr = self.ax.get_ylim()
                    lines = []
                    lines.extend(self.ax.plot(xr, (y, y), ":k"))
                    lines.extend(self.ax.plot((x, x), yr, ":k"))
                    self.annotations["x-hair"] = lines


                    s = "{fid}\nsamp {samp}\ntime {time} ns".format(
                                        fid=self.L.fids[x],
                                        samp=y,
                                        time=round(y/self.rate*1e9,2))
                    
                    ha = "left" if (x < self.L.data.shape[1]//2) else "right"
                    va = "top" if (y < self.L.data.shape[0]//2) else "bottom"
                    xoff = 1 if ha == "left" else -1
                    yoff = 10 if va == "top" else -10
                    txtbbox = dict(facecolor='k', alpha=0.2, pad=3)

                    txt = self.ax.text(x+xoff, y+yoff, s, size=10, color="w",
                                        weight="bold", ha=ha, va=va, bbox=txtbbox)
                    self.annotations["x-hair-text"] = [txt]
                    self.fig.canvas.draw()

                elif event.button == 2:
                    self.repaint()

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

    #def add_feature(self, s):
    #    if self.digitize_mode:
    #        self.end_feature()
    #    else:
    #        self.digitize_mode = True
    #    self.active_feature_name = s
    #    self.active_coords = []
    #    self.update()
    #    return self.fid

    #def end_feature(self):
    #    self.features[self.fid] = [self.active_feature_name,
    #                        [xy for xy in self.active_coords]]
    #    self.digitize_mode = False
    #    self.fid += 1
    #    self.update()
    #    return self.fid-1

    #def remove_feature(self, fid):
    #    try:
    #        self.features.pop(fid)
    #        self.active_coords = []
    #        self.update()
    #    except KeyError:
    #        pass

    #def add_point(self, event):
    #    """ Record a vertex. """
    #    self.active_coords.append((event.xdata, event.ydata))
    #    self.update()

    #def remove_last_point(self):
    #    try:
    #        self.active_coords.pop()
    #        self.update()
    #    except IndexError:
    #        pass

    def remove_annotation(self, name):
        """ Properly remove an annotation from the Radargram axes. """
        if name in self.annotations:
            self.annotations[name][0].remove()
            del self.annotations[name]
            return True
        else:
            return False

    def repaint(self, lum_scale=None, **kwargs):
        """ Redraw the radargram raster """
        if lum_scale is None:
            lum_scale = self.lum_scale
        else:
            self.lum_scale = lum_scale

        lum_bound = max((abs(self.data.max()), abs(self.data.min()))) * lum_scale
        ybnds = self.bbox[2:]
        data_ybnds = [self.data.shape[0]-1, 0]

        self.ax.cla()
        if ybnds[0] is None:
            y0 = data_ybnds[0]
        else:
            y0 = ybnds[0]
        if ybnds[1] is None:
            y1 = data_ybnds[1]
        else:
            y1 = ybnds[1]
        self.ax.set_ylim((y0, y1))

        self.ax.imshow(self.data, aspect='auto', cmap=self.cmap, vmin=-lum_bound, vmax=lum_bound)
        locs = np.arange(self.ax.get_ylim()[1], self.ax.get_ylim()[0], 50)
        self.ax.set_yticks(locs)
        self.ax.set_yticklabels(np.round(locs / self.rate * 1e9).astype(int))
        self.fig.canvas.draw()
        self.update(**kwargs)
        return

    def update(self):
        """ Display a radargram on axes. Paints in background, and
        all subsequent calls update lines. Passing repaint as True
        forces the background to be redrawn (for example, after a
        filter opperation).
        """
        n = self.data.shape[0]
        # These next 2 lines of code cause problems in matplotlib 3.5 - can't set them
        #https://matplotlib.org/stable/api/prev_api_changes/api_changes_3.5.0.html?highlight=axes.lines#behaviour-changes
        #self.ax.lines = []   ##TODO need workaround as this throws an error (can't set)
        #self.ax.texts = []   ##TODO need workaround as this throws an error (can't set)

        # Draw nodes
        drawxy = lambda xy: self.ax.plot(xy[0], xy[1], 'or', markersize=5.0,
                                          markeredgewidth=0.0, alpha=0.5)
        points_ = map(drawxy, self.active_coords)

        # Draw annotations
        for annotation in self.annotations:
            for item in self.annotations[annotation]:
                self.ax.add_artist(item)

        ## Draw previous features
        #if len(self.features) > 0:

        #    drawline = lambda lsxy: self.ax.plot(
        #        [i[0] for i in lsxy[1]], [i[1] for i in lsxy[1]],
        #        '--r')
        #    lines_ = map(drawline, self.features.values())

        #    labelfeature = lambda key, lsxy: self.ax.text(lsxy[1][-1][0],
        #        lsxy[1][-1][1]-20, str(key), fontsize=12, color='r')
        #    text_ = map(labelfeature, self.features.keys(), self.features.values())

        # Force tight bounding
        if self.data.shape[1] > 1:
            self.ax.set_xlim([0, self.data.shape[1]-1])
        else:
            self.ax.set_xlim([-0.5, 0.5])

        # Decorate and draw
        self.ax.set_ylabel("Time (ns)")
        self.ax.set_xlabel("Trace number")
        if self.digitize_mode:
            self.ax.set_title("Line {0} [feature {1}]".format(self.L.line, self.fid))
        else:
            self.ax.set_title("Line {0} [viewing]".format(self.L.line))
        self.fig.canvas.draw()
        return

    def get_digitizer_filename(self):
        """ Return an automatically-generated filename for digitized features. """
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
    """ Show a series of A-scan traces and allow picking """

    annotations = {}
    spacing = 0.1
    rg = None

    def __init__(self, L, ntraces=8):
        """Parameters
        ----------
        L : Gather instance containing the radar data to be displayed
        ntraces : number of consecutive traces to display at once
        """
        super(PickWindow, self).__init__((10, 8))
        self.fig.canvas.manager.set_window_title("Picker")
        self.ax.set_autoscale_on(False)

        # Set up scale slider widget
        self.ax_slider = self.fig.add_axes([0.2, 0.02, 0.4, 0.03])
        self.trace_scale = 0.4
        self.slider = matplotlib.widgets.Slider(self.ax_slider, "Amplitude Scale",
                                                0.1, 2.0,
                                                valinit=self.trace_scale)
        self.slider.on_changed(self._set_trace_scale)

        # Set up change mode widget
        self.ax_mode = self.fig.add_axes([0.8, 0.01, 0.05, 0.06])
        self.modebuttons = matplotlib.widgets.RadioButtons(self.ax_mode,
                                                    ["Bed", "DC"], 0)
        def change_mode(txt):
            self.change_mode(txt)
            self.update()
        self.modebuttons.on_clicked(change_mode)

        self.mode = "bed"
        self.trace0 = 0
        self.activetrace = None
        self.ntraces = ntraces
        self._newline(L)
        return

    #def __del__(self):
    #    super(PickWindow, self).__del__()

    def _newline(self, L):
        """ Replace internal state with a new radar line, discarding all
        digitizing and feature information and redrawing. """
        self.rate = L.metadata.sample_rate[0]
        self.L = irlib.PickableGather(L)
        self.data = L.data
        self.time = np.arange(L.data.shape[0]) / self.rate
        self.ylim = [-self.time[-1], -self.time[0]]

        self.bed_points = np.nan * np.zeros(self.data.shape[1])
        self.dc_points = np.nan * np.zeros(self.data.shape[1])
        self.bed_pick_reg = [None for i in range(self.ntraces)]
        self.dc_pick_reg = [None for i in range(self.ntraces)]

        self.points = self.bed_points

        self.trace0 = 0
        self.annotations = {}
        self.update()
        return

    def _onclick(self, event):
        """ Event handler for mouse clicks. Attempts to place a pick
        point where the mouse was clicked. """

        if event.button == 2:

            newmode = "bed" if (self.mode == "dc") else "dc"
            self.change_mode(newmode)

        else:

            # Identify the trace that was aimed for
            if event.xdata == None or event.inaxes is not self.ax:
                return
            activetrace = int(round(event.xdata/self.spacing + self.ntraces/2))
            try:
                trace = self.data[:,self.trace0 + activetrace]
            except IndexError:
                return

            # Determine if trace already has a point, and if so, remove it
            pr = self._active_reg()
            if pr[activetrace] is not None:
                self.ax.lines.remove(pr[activetrace][0])
                pr[activetrace] = None

            if event.button == 1:
                # Now put a point there
                yi = round(-event.ydata * self.rate)
                self._drawpick(trace, yi, activetrace)

                # Save the picked point
                self.points[activetrace + self.trace0] = yi

                self.yi = yi
                self.activetrace = activetrace

            elif event.button == 3:
                # Remove the point
                pr[activetrace] = None
                self.points[activetrace + self.trace0] = np.nan
                self.activetrace = None

        self.update()
        return

    def _set_trace_scale(self, val):
        self.trace_scale = val
        self.update()
        return

    def _active_reg(self):
        if self.mode == "bed":
            return self.bed_pick_reg
        else:
            return self.dc_pick_reg

    def _onkeypress(self, event):
        """ Event handler for keystrokes. Moves pick locations or loads
        different traces.
        """

        # See about moving a placed pick
        if self.activetrace is not None:

            trace = self.data[:,self.trace0 + self.activetrace]

            if event.key in ('j', 'down'):
                self._clearlastpick()
                self.yi += 1
                self._drawpick(trace, self.yi, self.activetrace)
                self.points[self.activetrace + self.trace0] = self.yi

            elif event.key in ('k', 'up'):
                self._clearlastpick()
                self.yi -= 1
                self._drawpick(trace, self.yi, self.activetrace)
                self.points[self.activetrace + self.trace0] = self.yi

            self.fig.canvas.draw()
            self.update_radargram()

        if event.key in ('h', 'left'):
            # Try moving to the left
            if self.trace0 > 0:
                self.trace0 -= self.ntraces
                self.activetrace = None
                self._clear_registers()
                self.update()

        elif event.key in ('l', 'right'):
            # Try moving to the right
            if self.trace0 < self.data.shape[1] - self.ntraces:
                self.trace0 += self.ntraces
                self.activetrace = None
                self._clear_registers()
                self.update()

        return

    def _shiftx(self, n):
        """ Shift an x-coordinate to the proper trace position on a
        multi-trace plot. """
        return n*self.spacing - self.ntraces/2*self.spacing

    def _clear_registers(self):
        self.bed_pick_reg = [None for i in self.bed_pick_reg]
        self.dc_pick_reg = [None for i in self.dc_pick_reg]
        return

    def _clearlastpick(self):
        """ Remove the last pick point drawn, according to the
        registry. """
        pr = self._active_reg()
        if self.activetrace is not None:
            self.ax.lines.remove(pr[self.activetrace][0])
            pr[self.activetrace] = None
        return

    def _drawpick(self, trace, yi, activetrace):
        """ Draw a pick mark on the given trace at yi, and update the
        registry. """
        trace = trace / abs(trace).max() * self.spacing * self.trace_scale

        if self.mode == 'bed':
            pr = self.bed_pick_reg
            sty = "ob"
        elif self.mode == 'dc':
            pr = self.dc_pick_reg
            sty = "or"

        pr[activetrace] = self.ax.plot(trace[int(yi)] + self._shiftx(activetrace),
                                       -self.time[int(yi)], sty, alpha=0.4)

        return self.ax.lines

    def _get_pick_fnm(self):
        """ Autogenerate a filename for pickfiles. """
        fnm = os.path.join('picking', '{0}_line{1}.csv'.format(
                                os.path.basename(self.L.infile).split('.')[0],
                                self.L.line))
        return fnm

    def update(self):
        """ Redraw axes and data """
        
        self.ax.cla()
        self.ax.set_xlim(-self.spacing * (self.ntraces+2) / 2,
                          self.spacing * self.ntraces / 2)
        self.ax.set_ylim(self.ylim)

        for i in range(self.ntraces):
            trno = self.trace0 + i

            if trno < self.data.shape[1]:
                # Plot the trace
                trace = self.data[:,trno]
                trace = trace / abs(trace).max() * self.spacing * self.trace_scale
                self.ax.plot(trace + self._shiftx(i), -self.time, '-k')

                # Plot any existing picks
                oldmode = self.mode
                if not np.isnan(self.bed_points[trno]):
                    self.mode = 'bed'
                    self._drawpick(trace, self.bed_points[trno], i)
                if not np.isnan(self.dc_points[trno]):
                    self.mode = 'dc'
                    self._drawpick(trace, self.dc_points[trno], i)
                self.mode = oldmode

        # using fixed locator for labels
        locs = self.ax.get_yticks().tolist()
        self.ax.yaxis.set_major_locator(mticker.FixedLocator(locs))
        self.ax.set_yticklabels(["{:0d}".format(np.round(x*-1e9).astype(int)) for x in locs])

        self.ax.set_title("Line {0}, mode: {1}".format(self.L.line, self.mode))

        self.fig.canvas.draw()
        try:
            self.update_radargram()
        except ConnectionError:
            pass

        return

    def update_radargram(self):
        """ Send picking annotations to a connected Radargram. """
        if self.rg is not None:
            for name in ("picklim0", "picklim1", "bedpick", "dcpick"):
                self.rg.remove_annotation(name)
            yl = self.rg.ax.get_ylim()

            self.rg.annotations["picklim0"] = \
                self.rg.ax.plot((self.trace0, self.trace0), yl, "-y")
            self.rg.annotations["picklim1"] = \
                self.rg.ax.plot((self.trace0+self.ntraces, self.trace0+self.ntraces), yl, "-y")

            self.rg.annotations["bedpick"] = \
                self.rg.ax.plot(self.bed_points, "-b", alpha=0.5)
            self.rg.annotations["dcpick"] = \
                self.rg.ax.plot(self.dc_points, "-r", alpha=0.7)

            self.rg.fig.canvas.draw()
        else:
            raise ConnectionError("No Radargram instance attached")
        return

    def change_mode(self, mode):
        """ Change picking mode between bed and direct coupling. """
        self.mode = mode.lower()
        if mode.lower() == 'bed':
            self.points = self.bed_points
        elif mode.lower() == 'dc':
            self.points = self.dc_points
        return

    def connect_radargram(self, rg):
        """ Connect a Radargram instance so that the PickWindow can modify it. """
        self.rg = rg
        self.update()
        return

    def autopick_dc(self, t0=10, tf=150):
        """ Attempt to pick the first break of the direct-coupling wave.
        Optional constraints on start and end time can be passed to improve
        results.

        Parameters
        ----------
        t0 : start time, in samples
        tf : end time, in samples
        """
        self.L.data = self.data
        self.L.PickDC(sbracket=(t0, tf))
        self.dc_points = self.L.dc_picks
        if self.mode == "dc":
            self.points = self.dc_points
        self.update()
        if self.rg is not None:
            self.rg.update()
        return

    def autopick_bed(self, t0=150, tf=10000, lbnd=None, rbnd=None):
        """ Attempt to pick the first break of the direct-coupling wave.
        Optional constraints on start and end time can be passed to improve
        results.

        Parameters
        ----------
        t0 : start time, in samples
        tf : end time, in samples
        """
        self.L.data = self.data
        self.L.PickBed(sbracket=(t0, tf), bounds=(lbnd, rbnd), phase=1)
        self.bed_points = self.L.bed_picks
        if self.mode == "bed":
            self.points = self.bed_points
        self.update()
        if self.rg is not None:
            self.rg.update()
        return

    def save_picks(self, fnm=None):
        """ Save picks. If no fnm is provided, generate one based on self.L. """
        if fnm is None:
            fnm = self._get_pick_fnm()
        fh = irlib.FileHandler(fnm, self.L.line, fids=self.L.fids)
        fh.AddBedPicks(self.L.fids, self.bed_points)
        fh.AddDCPicks(self.L.fids, self.dc_points)
        fh.ComputeTravelTimes()
        fh.Write()
        return

    def load_picks(self, fnm=None):
        """ Load picks. If no *fnm* is provided, attempt to load from an
        autogenerated location. """
        if fnm is None:
            fnm = self._get_pick_fnm()
        if os.path.isfile(fnm):
            fh = irlib.FileHandler(fnm, self.L.line)
            dc_points, bed_points = fh.GetEventValsByFID(self.L.fids)
            self.dc_points = dc_points
            self.bed_points = bed_points
            self.change_mode(self.mode)
            self.update()
            if self.rg is not None:
                self.rg.update()
        else:
            sys.stderr.write("pickfile ({0}) does not exist\n".format(fnm))
        return


class MapWindow(AppWindow):
    """ Displays a simple map of trace locations """

    def __init__(self, L):
        super(MapWindow, self).__init__((4, 4))
        self._newline(L)
        self.fig.canvas.manager.set_window_title("Map") 
        self.fig.tight_layout()
        return

    def _newline(self, L):
        """ Replace internal state with a new radar line, discarding all
        digitizing and feature information and redrawing. """
        self.L = L
        self.x = L.metadata.eastings
        self.y = L.metadata.northings
        self.ax.cla()
        self.update()
        return

    def update(self):
        """ Redraw the map. """
        self.ax.set_xlabel("Eastings (m)")
        self.ax.set_ylabel("Northings (m)")
        self.ax.plot(self.x, self.y, ".k")
        self.ax.axis("equal")
        self.fig.canvas.draw()

class ConnectionError(Exception):
    def __init__(self, message="No message"):
        self.message = message
    def __str__(self):
        return self.message

