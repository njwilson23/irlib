
import matplotlib.pyplot
import commands
import command_parser as cp
from .components import MapWindow

class MapCommand(commands.Command):

    _type = "Map"
    cmd = "map"
    helpstr = """Manage map window

    map <cmd>

    Commands related to the map window. Type 'map', and then the specific
    command. Type 'help' for a list of possibilities.
    """
    def apply(self, app, args):
        if len(args) == 0:
            print "Type 'help' or 'help map' for instructions."
        try:
            cp.apply_command(app.command_registry, args, app)
        except KeyError:
            print "No mapping command '{0}' exists".format(args[0])
        return

class MapOn(MapCommand):

    cmd = "on"
    helpstr = "\tOpen a map.\n"

    def apply(self, app, args):
        w = MapWindow(app.line)
        app.appwindows.append(w)
        return

class MapOff(MapCommand):

    cmd = "off"
    helpstr = "\tClose map.\n"

    def apply(self, app, args):
        for w in app.get_appwindows(MapWindow):
            matplotlib.pyplot.close(w.fig)
            app.remove_appwindow(w)
            del w


