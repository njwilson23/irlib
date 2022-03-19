""" The Console class is the controller for irlib-based apps. Windows can be
attached and detached from the Console, which handles user input and passes
directives to its windows. """

from __future__ import print_function

import sys
import os
import getopt
try:
  import readline
except ImportError:
  import pyreadline as readline
import atexit
import matplotlib.pyplot as plt
import six
from .. import gather
from .. import survey
from .. import misc
from . import command_parser
from . import commands
from .components import Radargram, MapWindow, PickWindow

class Console(object):
    """ App-controller with a user input readline loop. """

    survey = None
    line = None
    appwindows = []
    readline_hist = ".history"

    def __init__(self, progname, bannertext=""):
        """ Initiate a new Console instance.

        Parameters
        ----------
        progname : Name of the application [string]
        bannertext : optional lines to display when application starts [string]
        """

        self.progname = progname
        self.bannertext = bannertext

        try:
            optlist, args = getopt.gnu_getopt(sys.argv[1:], 'f:L:')
        except getopt.GetoptError:
            print("Error collecting arguments - check syntax.")
            self.print_syntax()
            sys.exit(1)
        optdict = dict(optlist)

        try:
            self.infile = optdict['-f']
        except KeyError:
            if len(args) > 0:
                self.infile = args[0]
            else:
                self.print_syntax()
                sys.exit(0)

        lineno = int(optdict.get('-L', 0))

        try:
            self.survey = survey.Survey(self.infile)
        except IOError as e:
            print(e)
            raise e
        self.open_line(lineno)
        self.appwindows.append(Radargram(self.line))

        # Register basic commands
        self.command_registry = {}
        self.register(commands)

        # Setup readline
        try:
            readline.read_history_file(self.readline_hist)
        except IOError:
            pass
        atexit.register(readline.write_history_file, self.readline_hist)

        return

    def start(self):
        """ Begin input-output loop """
        print(self.bannertext)

        while True:
            cmd = self.get_command()
            self.handle_command(cmd)
        return

    def register(self, module):
        """ Load commands from a module a register them for use. """
        for item in module.__dict__.values():
            if isinstance(item, type) and commands.Command in item.mro():
                if item.cmd not in self.command_registry:
                    self.command_registry[item.cmd] = {}
                self.command_registry[item.cmd][item._type] = item
        return

    def print_syntax(self):
        """ Print start-up syntax for the forgetful. """
        print("\tUSAGE: {0} <HDF_survey> [-L line_number]".format(self.progname))
        return

    def open_line(self, lineno, dcno=0, fromcache=True):
        """ Open a line from a survey """
        loaded = False
        if fromcache:
            cachename = self.survey.GetLineCacheName(lineno, dcno)
            loaded, line = misc.TryCache(cachename)
            self.line = line
        if loaded is False:
            line = self.survey.ExtractLine(lineno, datacapture=dcno)
            if line.nx >= 2:
                try:
                    line.RemoveBadLocations()
                    line.RemoveGPSNaNs()
                    line.FixStaticGPS()
                    line.RemoveBlankTraces()
                    line.RemoveStationary(threshold=3.0)
                    self.line = line
                except gather.LineGatherError:
                    pass
            elif line.nx == 1:
                self.line = line
            else:
                print("line {0}:{0} contains no data".format(lineno, dcno))
                self.line = None
        return

    def get_command(self):
        """ Get a command from console input. """
        cmd = six.moves.input('>> ')
        return cmd

    def get_appwindows(self, t=None):
        """ Get all windows of a particular type from the window list. """
        if t is None:
            return self.appwindows
        elif not hasattr(t, "__iter__"):
            return [a for a in self.appwindows if type(a) == t]
        else:
            windows = []
            for _t in t:
                windows.extend(self.get_appwindows(_t))
            return windows

    def add_appwindow(self, ref):
        """ Add a window to be managed by the Console. """
        self.appwindows.append(ref)
        if isinstance(ref, PickWindow):
            rgs = self.get_appwindows(Radargram)
            if len(rgs) > 0:
                ref.connect_radargram(rgs[0])
        return

    def remove_appwindow(self, ref):
        """ Remove a window from Console management. """
        for i, win in enumerate(self.appwindows):
            if win is ref:
                break
        self.appwindows.pop(i)
        return

    def handle_command(self, cmd):
        """ Parse and handle user input. This will be deprecated int he future
        for a more modular approach that allows commands to be added and
        removed dynamically. """
        # List args. Handle empty input.
        if cmd == '':
            return

        args = cmd.split(' ')

        # Remove double spaces
        while 1:
            try:
                args.remove('')
            except ValueError:
                break

        try:
            command_parser.apply_command(self.command_registry, args, self,
                                         "General")
        except KeyError:
            print("No command '{0}' exists".format(args[0]))
        return

