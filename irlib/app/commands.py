""" Define a Command class for interactive apps. """

import sys

class Command(object):

    _type = "Command"
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

    Print information about the currently open survey and radar line. """

    def apply(self, app, args):
        print app.survey.datafile
        print 'line: ' + str(app.line.line)
        try:
            print 'channel: ' + str(app.line.datacapture)
        except AttributeError:
            pass
        print '# traces: ' + str(app.line.data.shape[1])
        print '# samples: ' + str(app.line.data.shape[0])
        rate = 1.0 / app.line.metadata.sample_rate[0]
        print 'sample interval: ' + str(rate) + ' s'
        print 'depth resolution: ' + str(rate * 1.68e8 / 2.0) + ' m'
        print 'vertical range: ' + \
                str(1.68e8*app.line.data.shape[0]*rate / 2.0) + ' m'
        print 'available channels: ' + \
                str(app.survey.GetChannelsInLine(int(app.line.line)))
        return

class ListLines(Command):

    cmd = "ls"
    helpstr = "List available lines."

    def apply(self, app, args):
        cursor = 2
        print "Lines:"
        print "  ",
        for linestr in app.survey.GetLines():
            lineno = int(linestr.split("_")[1])
            if app.line.line == lineno:
                print "<{0}>".format(lineno),
            else:
                print " {0} ".format(lineno),
            cursor += (2 + len(str(lineno)))
            if cursor > 60:
                print
                print "  ",
                cursor = 2
        print

class OpenLine(Command):

    cmd = "open"
    helpstr = """Open Line

    open <lineno> [dcno=0]

    Attempts to open the line with number *lineno* and the channel with number
    *dcno* (default 0)."""

    def apply(self, app, args):
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
            if 'line_{0}'.format(lineno) in app.survey.GetLines() and \
                dcno < app.survey.GetChannelsInLine(lineno):
                print "Opening line {0}, channel {1}".format(lineno, dcno)
                del app.line
                app.open_line(lineno, dcno=dcno)

                for w in app.get_appwindows((Radargram, MapWindow, PickWindow)):
                    w._newline(app.line)

            else:
                print ("Line {0} channel {1} does "
                       "not exist".format(lineno, dcno))
        except irlib.EmptyLineError:
            print ("Line {0} channel {1} could not be opened because it "
                   "contains no data".format(lineno, dcno))
        except:
            traceback.print_exc()

class ApplyFilter(Command):

    cmd = "filter"
    helpstr = """Apply Filter

    filter <filtername> [args]

    Applies a filter to the data. For a list of filters, type "help". For help
    on individual filters, type "help <filtername>".

    Can be abbreviated as "f"."""

    def apply(self, app, args):
        if len(args) == 1:
            print app.line.PprintHistory()
        else:
            try:
                command_parser.apply_filter(args[1:], app.line)
                for rg in app.get_appwindows(Radargram):
                    rg.data = app.line.data
                    rg.repaint()
                for pw in app.get_appwindows(PickWindow):
                    pw.data = app.line.data
                    pw.update()
            except command_parser.CommandSearchError as e:
                print e.message

class ApplyFilterAlt(ApplyFilter):
    cmd = "f"

class NoFilter(Command):

    cmd = "nofilter"
    helpstr = "Remove filters, restoring the radar data to it's original state."

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
            gain = float(args[1])
            for rg in app.get_appwindows(Radargram):
                rg.repaint(lum_scale=1.0/gain)
        except IndexError:
            print "gain: " + str(1.0 / app.get_appwindows(Radargram)[0].lum_scale)

class YLimAdjuster(Command):

    cmd = "ylim"
    helpstr = """Adjust the vertical display limits

    ylim [value]

    Modify the vertical limits in all applicable windows. If no argument is
    supplied, the current limits are returned."""

    def apply(self, app, args):
        rate_ns = app.line.rate * 1e9
        if len(args) == 0:
            rg = app.get_appwindows(Radargram)[0]
            ylim_ns = [a*rate_ns for a in rg.ax.get_ylim()][::-1]
            print "Vertical range: {0} - {1} ns".format(*ylim_ns)
        elif len(args) == 2:
            try:
                for rg in app.get_appwindows(Radargram):
                    rg.bbox[2:] = [float(a)*rate_ns for a in args][::-1]
                    rg.repaint()
            except ValueError:
                print "Could not understand '{0}'".format(args)
        else:
            print "'ylim' must be followed by times in nanoseconds"
        return

class Debug(Command):

    cmd = "debug"
    helpstr = "If this is useful to you, you already know what it does."

    def apply(self, app, args):
        import pdb
        pdb.set_trace()
        return

class HelpPrinter(Command):

    cmd = "help"
    helpstr = "Let's be serious..."

    def apply(self, app, args):
        if len(args) == 1:
            print """\tBasic application commands:

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
            """

            commandtypes = set([a._type for a in app.command_registry.values()])

            for ct in filter(lambda a: a is not "Command", commandtypes):
                print"\n\tAvailable {0} commands\n".format(ct)
                for name in (a for a in app.command_registry if a._type == ct):
                    print "\t\t{0}".format(name)

        else:
            print command_parser.help_filter(args[1])
        return


