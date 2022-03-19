""" Define a Command class for interactive apps. """

from __future__ import print_function

import sys
import itertools
import traceback

from . import command_parser as cp
from .components import Radargram, MapWindow, PickWindow
from ..survey import EmptyLineError

class Command(object):

    _type = "General"
    cmd = "__commandclass"
    helpstr = "This is a base class for Commands. It doesn't do anything."

    def apply(self, app, args):
        """ Overload this method to manipulate an Console instance. """
        raise Exception("apply() is undefined in the Command base class")

class Exit(Command):

    cmd = "exit"
    helpstr = "Exit"

    def apply(self, app, args):
        sys.exit(0)

class ExitAlt(Exit):
    cmd = "q"

class PrintInfo(Command):

    cmd = "info"
    helpstr = """Print Information

    info

    Print information about the currently open survey and radar line.
    """

    def apply(self, app, args):
        print(app.survey.datafile)
        print('line: ' + str(app.line.line))
        try:
            print('channel: ' + str(app.line.datacapture))
        except AttributeError:
            pass
        print('# traces: ' + str(app.line.data.shape[1]))
        print('# samples: ' + str(app.line.data.shape[0]))
        rate = 1.0 / app.line.metadata.sample_rate[0]
        print('sample interval: ' + str(rate) + ' s')
        print('depth resolution: ' + str(rate * 1.68e8 / 2.0) + ' m')
        print('vertical range: ' + str(1.68e8*app.line.data.shape[0]*rate / 2.0) + ' m')
        print('available channels: ' + str(app.survey.GetChannelsInLine(int(app.line.line))))
        return

class ListLines(Command):

    cmd = "ls"
    helpstr = "\tList available lines.\n"

    def apply(self, app, args):
        cursor = 2
        print("Lines:")
        print("  ", end="")
        for linestr in app.survey.GetLines():
            lineno = int(linestr.split("_")[1])
            if app.line.line == lineno:
                print("<{0}>".format(lineno), end="")
            else:
                print(" {0} ".format(lineno), end="")
            cursor += (2 + len(str(lineno)))
            if cursor > 50:
                print("\n  ", end="")
                cursor = 2
        print("\n")

class OpenLine(Command):

    cmd = "open"
    helpstr = """Open Line

    open <lineno> [dcno=0]

    Attempts to open the line with number *lineno* and the channel with number
    *dcno* (default 0).
    """

    def apply(self, app, args):
        try:
            lineno = int(args[0])
        except (IndexError, ValueError):
            print("Must supply an integer line number")
            return

        try:
            dcno = int(args[1])
        except IndexError:
            dcno = 0
        except ValueError:
            print("Bad channel number: {0}".format(args[1]))
            return

        try:
            if 'line_{0}'.format(lineno) in app.survey.GetLines() and \
                dcno < app.survey.GetChannelsInLine(lineno):
                print("Opening line {0}, channel {1}".format(lineno, dcno))
                del app.line
                app.open_line(lineno, dcno=dcno)

                for w in app.get_appwindows((Radargram, MapWindow, PickWindow)):
                    w._newline(app.line)

            else:
                print("Line {0} channel {1} does "
                      "not exist".format(lineno, dcno))
        except EmptyLineError:
            print("Line {0} channel {1} could not be opened because it "
                  "contains no data".format(lineno, dcno))
        except:
            traceback.print_exc()

class ApplyFilter(Command):
    """ This is a "General" command that dispatched to "Filter" commands. """

    cmd = "filter"
    helpstr = """Apply Filter

    filter <filtername> [args]

    Applies a filter to the data. For a list of filters, type "help". For help
    on individual filters, type "help <filtername>".

    Can be abbreviated as "f".
    """

    def apply(self, app, args):
        if len(args) == 0:
            print(app.line.PprintHistory())
        else:
            try:
                # args are the inputs -> 'filter'
                cp.apply_command(app.command_registry, args, app.line, "Filter")
                for rg in app.get_appwindows(Radargram):
                    rg.data = app.line.data
                    rg.repaint()
                for pw in app.get_appwindows(PickWindow):
                    pw.data = app.line.data
                    pw.update()
            except KeyError:
                print("No filter type '{0}' exists".format(args[0]))

class ApplyFilterAlt(ApplyFilter):
    cmd = "f"

class NoFilter(Command):

    cmd = "nofilter"
    helpstr = "\tRemove filters, restoring the radar data to it's original state.\n"

    def apply(self, app, args):
        app.line.Reset()
        for rg in app.get_appwindows(Radargram):
            rg.data = app.line.data
            rg.repaint()
        for pw in app.get_appwindows(PickWindow):
            pw.data = app.line.data
            pw.update()
        return

class NoFilterAlt(NoFilter):
    cmd = "nf"

class GainAdjuster(Command):

    cmd = "gain"
    helpstr = """Adjust radargram contrast

    gain [value]

    Apply a contrast enhancement determined by *value*. If no argument is
    supplied, the current constrast enhancement is returned."""

    def apply(self, app, args):
        try:
            gain = float(args[0])
            for rg in app.get_appwindows(Radargram):
                rg.repaint(lum_scale=1.0/gain)
        except IndexError:
            print("\tgain: {0}".format(1.0 / app.get_appwindows(Radargram)[0].lum_scale))

class YLimAdjuster(Command):

    cmd = "ylim"
    helpstr = """Adjust the vertical display limits

    ylim [t0 t1]

    Print or modify the vertical limits in all applicable windows. If no argument is
    supplied, the current limits are returned.
    """

    def apply(self, app, args):
        rate_ns = app.line.rate * 1e9
        if len(args) == 0:
            rg = app.get_appwindows(Radargram)[0]
            ylim_ns = [a*rate_ns for a in rg.ax.get_ylim()][::-1]
            print("Vertical range: {0} - {1} ns".format(*ylim_ns))
        elif len(args) == 2:
            try:
                for rg in app.get_appwindows(Radargram):
                    rg.bbox[2:] = [float(a)/rate_ns for a in args][::-1]
                    rg.repaint()
                for w in app.get_appwindows(PickWindow):
                    w.ylim = [-float(a)*1e-9 for a in args][::-1]
                    w.update()
            except ValueError:
                print("Could not understand '{0}'".format(args))
        else:
            print("Incorrect expression, type 'help ylim'")
        return

class SaveImage(Command):
    cmd = "saveimage"
    helpstr = """Save image

    saveimg <filename>

    Save the current radargram image to a PNG file.
    """

    def apply(self, app, args):
        if len(args) > 0:
            fnm = args[0]
        else:
            raise ValueError("Command '{0}' must by followed by a "
                             "filename".format(self.cmd))
        rg = app.get_appwindows(Radargram)[0]
        rg.fig.savefig(fnm, dpi=300)
        return

class Debug(Command):

    cmd = "debug"
    helpstr = "\tIf this is useful to you, you already know what it does.\n"

    def apply(self, app, args):
        import pdb
        pdb.set_trace()
        return

class HelpPrinter(Command):
    """ This is pretty gross, but it works. """

    cmd = "help"
    helpstr = "\tLet's be serious...\n"

    def apply(self, app, args):
        if len(args) == 0:
            print("""\tBasic application commands:

            info                print line metadata
            ls                  list lines in survey
            open [line#]        open a different line
            gain [#]            adjust radargram display contrast
            ylim [t0 tf]        adjust vertical radargram range in nanoseconds

            filter *name*       apply a filter (see below)
            nofilter            remove all filters

            pick on|off         open a picking window
            map on|off          open a map of the displayed line

            exit                exit irview
            debug
            """)

            cmdobjs = [b for b in itertools.chain(*[a.values() for a in app.command_registry.values()])]
            commandtypes = set([a._type for a in cmdobjs])

            for ct in filter(lambda a: a != "General", commandtypes):
                print("\n\tAvailable {0} commands\n".format(ct))
                for cmdobj in filter(lambda a: a._type == ct, cmdobjs):
                    if cmdobj.cmd != None and not cmdobj.cmd.startswith("_"):
                        print("\t{0}".format(cmdobj.cmd))

        else:
            if len(app.command_registry.get(args[-1], [])) > 1:
                print("The following help topics were found:")
                multiple_topics = True
            else:
                multiple_topics = False

            try:
                for cmdobj in app.command_registry[args[-1]].values():
                    if multiple_topics:
                        print("{0} {1}:".format(cmdobj._type.lower(), cmdobj.cmd))
                    print(cmdobj.helpstr)
            except KeyError:
                raise KeyError("No command definition '{0}' found".format(args[-1]))

        return


