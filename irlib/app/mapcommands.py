
from __future__ import print_function
import matplotlib.pyplot
from . import commands
from . import command_parser as cp
from .components import MapWindow

class MapCall(commands.Command):
    """ This is a "General" type command that dispatches to "Map" type commands.
    """

    _type = "General"
    cmd = "map"
    helpstr = """Manage map window

    map <cmd>

    Commands related to the map window. Type 'map', and then the specific
    command. Type 'help' for a list of possibilities.
    """
    def apply(self, app, args):
        if len(args) == 0:
            print("Type 'help' or 'help map' for instructions.")
        try:
            cp.apply_command(app.command_registry, args, app, "Map")
        except KeyError:
            print("No mapping command '{0}' exists".format(args[0]))
        return

class MapCommandBase(commands.Command):
    _type = "Map"
    cmd = None
    helpstr = ""

class MapOn(MapCommandBase):

    cmd = "on"
    helpstr = "\tOpen a map.\n"

    def apply(self, app, args):
        w = MapWindow(app.line)
        app.appwindows.append(w)
        return

class MapOff(MapCommandBase):

    cmd = "off"
    helpstr = "\tClose map.\n"

    def apply(self, app, args):
        for w in app.get_appwindows(MapWindow):
            matplotlib.pyplot.close(w.fig)
            app.remove_appwindow(w)
            del w


