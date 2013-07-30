""" Define a Console instance, which is the central app controller as of irlib
v0.4. """

import sys
import getopt
import readline
import matplotlib.pyplot as plt
import irlib
import command_parser
from .components import Radargram, MapWindow, PickWindow
import traceback
import pdb

class Console(object):

    survey = None
    line = None
    appwindows = []

    def __init__(self, progname, bannertext=""):

        self.progname = progname

        try:
            optlist, args = getopt.gnu_getopt(sys.argv[1:], 'f:L:')
        except getopt.GetoptError:
            print "Error collecting arguments - check syntax."
            self.print_syntax()
            sys.exit(1)
        optdict = dict(optlist)

        try:
            self.infile = optdict['-f']
        except KeyError:
            print "A survey filename must be supplied:"
            self.print_syntax()
            sys.exit(0)

        lineno = int(optdict.get('-L', 0))

        self.open_line(lineno)
        self.appwindows.append(Radargram(self.line))

        # Begin main loop
        print bannertext

        while True:
            cmd = self.get_command()
            self.handle_command(cmd)
        return

    def print_syntax(self):
        print "\t {0} -f file_name [-L line_number]".format(self.progname)
        return

    def open_line(self, lineno, dcno=0, fromcache=True):
        """ Open a line from a survey """
        self.survey = irlib.Survey(self.infile)
        loaded = False
        if fromcache:
            cachename = self.survey.GetLineCacheName(lineno, dcno)
            loaded, line = irlib.misc.TryCache(cachename)
        if loaded == False:
            line = self.survey.ExtractLine(lineno, datacapture=dcno)
            try:
                line.RemoveBadLocations()
                line.FixStaticGPS()
                line.RemoveBlankTraces()
                line.SmoothenGPS()
                line.RemoveStationary(threshold=3.0)
            except irlib.LineGatherError:
                pass
        self.line = line
        return

    def get_command(self):
        """ Get a command from console input. """
        cmd = raw_input('>> ')
        return cmd

    def get_appwindows(self, t=None):
        """ Get all windows of a particular type from the window list. """
        if t is None:
            return self.appwindows
        elif not hasattr(t, "__iter__"):
            return [a for a in self.appwindows if type(a) == t]
        else:
            return reduce(lambda a,b: a+b, [self.get_appwindows(a) for a in t])

    def remove_appwindow(self, ref):
        """ Remove a window from the window list. """
        for i, win in enumerate(self.appwindows):
            if win is ref:
                break
        self.appwindows.pop(i)
        return

    def handle_command(self, cmd):
        # List args. Handle empty input.
        args = cmd.split(' ')

        if cmd == '':
            return

        # Remove double spaces
        while 1:
            try:
                args.remove('')
            except ValueError:
                break

        if args[0] in ('exit', 'q', 'quit'):    # EXIT
            sys.exit(0)

        elif args[0] == 'info':                 # INFO
            print self.survey.datafile
            print 'line: ' + str(self.line.line)
            try:
                print 'channel: ' + str(self.line.datacapture)
            except AttributeError:
                pass
            print '# traces: ' + str(self.line.data.shape[1])
            print '# samples: ' + str(self.line.data.shape[0])
            rate = 1.0 / self.line.metadata.sample_rate[0]
            print 'sample interval: ' + str(rate) + ' s'
            print 'depth resolution: ' + str(rate * 1.68e8 / 2.0) + ' m'
            print 'vertical range: ' + str(1.68e8*self.line.data.shape[0]*rate / 2.0) + ' m'
            print 'available channels: ' + str(self.survey.GetChannelsInLine(int(self.line.line)))

        elif args[0] == 'ls':                   # LS
            cursor = 2
            print "Lines:"
            print "  ",
            for linestr in self.survey.GetLines():
                lineno = int(linestr.split("_")[1])
                if self.line.line == lineno:
                    print "<{0}>".format(lineno),
                else:
                    print " {0} ".format(lineno),
                cursor += (2 + len(str(lineno)))
                if cursor > 70:
                    print
                    print "  ",
            print

        elif args[0] == 'open':                 # OPEN
            try:
                lineno = int(args[1])
            except (IndexError, ValueError):
                print "Must supply an integer line number"
                return

            try:
                dcno = int(args[2])
            except IndexError:
                dcno = 0
            except ValueError:
                print "Bad channel number: {0}".format(args[2])
                return

            try:
                if 'line_{0}'.format(lineno) in self.survey.GetLines() and \
                    dcno < self.survey.GetChannelsInLine(lineno):
                    print "Opening line {0}, channel {1}".format(lineno, dcno)
                    del self.line
                    self.open_line(lineno, dcno=dcno)

                    for w in self.get_appwindows((Radargram, MapWindow, PickWindow)):
                        w._newline(self.line)

                else:
                    print ("Line {0} channel {1} does "
                           "not exist".format(lineno, dcno))
            except:
                traceback.print_exc()

        elif args[0] in ('filter', 'f'):        # FILTER
            if len(args) == 0:
                print StrFilterHistory(self.line)
            else:
                try:
                    command_parser.apply_filter(args[1:], self.line)
                    for rg in self.get_appwindows(Radargram):
                        rg.data = self.line.data
                        rg.repaint()
                    for pw in self.get_appwindows(PickWindow):
                        pw.data = self.line.data
                        pw.update()
                except command_parser.CommandSearchError as e:
                    print e.message

        elif args[0] in ('nofilter', 'nf'):     # NOFILTER
            self.line.Reset()
            for rg in self.get_appwindows(Radargram):
                rg.data = self.line.data
                rg.repaint()
            for pw in self.get_appwindows(PickWindow):
                pw.data = self.line.data
                pw.update()

        #elif args[0] == 'dnew':                 # DNEW
        #    if len(args) == 1:
        #        args.append('')
        #    try:
        #        IW.AddFeature(args[1])
        #    except:
        #        traceback.print_exc()

        #elif args[0] == 'drm':                  # DRM
        #    try:
        #        IW.RemoveFeature(int(args[1]))
        #    except IndexError:
        #        print "No such feature"

        #elif args[0] == 'dls':                  # DLS
        #    print "Feature List"
        #    for key in IW.features.keys():
        #        print "{0}: {1} vertices".format(key, len(IW.features[key][1]))

        #elif args[0] == 'dsave':                # DSAVE
        #    try:
        #        dict_list = IW.Export()
        #        outfnm = IW.GetDigitizerFilename()
        #        WriteFeatures(IW, outfnm)
        #        print "Features saved to " + outfnm
        #    except:
        #        traceback.print_exc()

        #elif args[0] == 'dimport':              # DIMPORT
        #    infnm = IW.GetDigitizerFilename()
        #    try:
        #        with open(infnm, 'r') as fin:
        #            IW.Import(fin)
        #    except IOError:
        #        sys.stdout.write("Digitizer file ({0}) not found\n".format(infnm))
        #    except:
        #        traceback.print_exc()

        elif args[0] == 'imsave':               # IMSAVE
            try:
                fnm = args[1]
                if os.path.splitext(fnm)[1] != '.png':
                    fnm += '.png'
                self.get_appwindows(Radargram)[0].fig.savefig(fnm)
            except IndexError:
                sys.stdout.write("Please supply a save location/name\n")
            except:
                traceback.print_exc()

        elif args[0] == 'gain':                 # GAIN
            try:
                gain = float(args[1])
                for rg in self.get_appwindows(Radargram):
                    rg.repaint(lum_scale=1.0/gain)
            except IndexError:
                print "gain: " + str(1.0 / self.get_appwindows(Radargram)[0].lum_scale)

        elif args[0] == 'map':                  # MAP
            if len(args) < 2:
                if len(self.get_appwindows(MapWindow)) < 1:
                    print "Map window: off"
                else:
                    print "Map window: on"
            else:
                if args[1] == "on":
                    w = MapWindow(self.line)
                    self.appwindows.append(w)

                elif args[1] == "off":
                    for w in self.get_appwindows(MapWindow):
                        plt.close(w.fig)
                        self.remove_appwindow(w)
                        del w
                else:
                    print "Command not recognized"

        elif args[0] == 'pick':                 # PICK
            if len(args) < 2:
                if len(self.get_appwindows(MapWindow)) < 1:
                    print "Picking window: off"
                else:
                    print "Picking window: on"
            else:
                if args[1] == "on":
                    w = PickWindow(self.line)
                    w.connect_radargram(self.get_appwindows(Radargram)[0])
                    self.appwindows.append(w)

                elif args[1] == "off":
                    for w in self.get_appwindows(PickWindow):
                        plt.close(w.fig)
                        self.remove_appwindow(w)
                        del w
                elif args[1] == "save":
                    for w in self.get_appwindows(PickWindow):
                        w.save_picks()
                elif args[1] == "load":
                    for w in self.get_appwindows(PickWindow):
                        w.load_picks()
                elif args[1] == "bed":
                    for w in self.get_appwindows(PickWindow):
                        if len(args) > 2:
                            try:
                                pickargs = [int(a) for a in args[2:]]
                            except ValueError:
                                print "arguments must be integer"
                                return
                        else:
                            pickargs = []
                        w.autopick_bed(*pickargs)
                elif args[1] == "dc":
                    for w in self.get_appwindows(PickWindow):
                        if len(args) > 2:
                            try:
                                pickargs = [int(a) for a in args[2:]]
                            except ValueError:
                                print "arguments must be integer"
                                return
                        else:
                            pickargs = []
                        w.autopick_dc(*pickargs)
                else:
                    print "Command not recognized"

        elif args[0] == 'debug':                # DEBUG
            pdb.set_trace()

        elif args[0] == 'help':                 # HELP
            if len(args) == 1:
                print """\tApplication commands:

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

                map on|off          open a map of the displayed line

                exit                exit irview
                debug
                """

                print "\tAvailable filter commands:\n"
                for name in command_parser.list_filters():
                    print "\t\t{0}".format(name)

            else:
                print command_parser.help_filter(args[1])

        else:
            print "Command not recognized"

        return


