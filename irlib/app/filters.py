""" This file defines the filter commands available from irlib-based apps. """

from commands import Command

def handle_no_args():
    # intended to be used as a decorator - if there are no arguments, then print filter history as exit
    pass

class FilterCommandBase(Command):
    """ A FilterCommand is implemented as a class with a command-line signature and
    an *apply* method that takes an appropriate Gather object as an argument.
    """

    _type = "Filter"
    cmd = "__filtercommandclass"
    helpstr = ("This is the base class for FilterCommands, which can be used "
               "to modify Gather data in-place.")

    def __init__(self):
        return

    def apply(self, G, args):
        """ Overload this method to perform operations on Gather object *G*. """
        raise Exception("apply() is undefined in FilterCommand baseclass")

class LinearGainControl(FilterCommandBase):
    cmd = "gc"
    helpstr = """Linear gain control

    gc [n]

    Apply a time-dependent linear gain control (t^n) to each trace, with the
    exponent `n` taking the default value of 1.0. """

    def apply(self, G, args):
        if len(args) > 0:
            npow = float(args[0])
        else:
            npow = 1.0
        G.DoTimeGainControl(npow=npow)
        return

class AutoGainControl(FilterCommandBase):
    cmd = "agc"
    helpstr = """Automatic gain control

    agc

    Apply automatic gain control, normalizing the power in a finite window to a
    constant value."""

    def apply(self, G, args):
        G.Dewow()
        G.DoAutoGainControl(5e-8)
        return

class ReflectionPower(FilterCommandBase):
    cmd = "power"
    helpstr = """Reflection power

    power

    Square the trace amplitude.
    """
    def apply(self, G, args):
        G.data = G.data**2
        return

class Lowpass_FD(FilterCommandBase):
    cmd = "lowpass"
    helpstr = """Frequency-domain lowpass

    lowpass [cut [bw]]

    Apply a frequency domain lowpass filter implemented using a windowed sinc
    kernel. Optional arguments are cutoff frequency and bandwidth. """
    def apply(self, G, args):
        if len(args) > 0:
            co = float(args[0])
        else:
            co = 25e6
        if len(args) > 1:
            bw = float(args[1])
        else:
            bw = 5e6
        G. DoWindowedSinc(cutoff=co, bandwidth=bw, mode="lowpass")
        return

class Highpass_FD(FilterCommandBase):
    cmd = "highpass"
    helpstr = """Frequency-domain highpass

    highpass [cut [bw]]

    Apply a frequency domain highpass filter implemented using a windowed sinc
    kernel. Optional arguments are cutoff frequency and bandwidth. """
    def apply(self, G, args):
        if len(args) > 0:
            co = float(args[0])
        else:
            co = 25e6
        if len(args) > 1:
            bw = float(args[1])
        else:
            bw = 5e6
        G.DoWindowedSinc(cutoff=co, bandwidth=bw, mode='highpass')
        return

class Reverse(FilterCommandBase):
    cmd = "reverse"
    helpstr = """Reverse

    reverse

    Flip the traces in the radar line such that they are ordered from last to
    first.
    """
    def apply(self, G, args):
        G.Reverse()
        return

class Lowpass_TD(FilterCommandBase):
    cmd = "lowpass_td"
    helpstr = """Time-domain lowpass

    lowpass [nsamp]

    Apply a frequency domain lowpass filter implemented as a moving average
    with a blackman window. Optional argument is sample width of the filter. """
    def apply(self, G, args):
        if len(args) > 0:
            ns = int(args[0])
            if ns % 2 == 0:
                ns += 1
        else:
            ns = 21
        G. DoMoveAvg(ns, kind="blackman", mode="lowpass")
        return

class Highpass_TD(FilterCommandBase):
    cmd = "highpass_td"
    helpstr = """Time-domain highpass

    highpass [nsamp]

    Apply a frequency domain highpass filter implemented as a
    spectrally-inverted moving average with a blackman window. Optional
    argument is sample width of the filter. """
    def apply(self, G, args):
        if len(args) > 0:
            ns = int(args[0])
            if ns % 2 == 0:
                ns += 1
        else:
            ns = 7
        G.DoMoveAvg(ns, kind="blackman", mode='highpass')
        return

class Dewow(FilterCommandBase):
    cmd = "dewow"
    helpstr = """Dewow

    dewow

    Apply a "dewowing" filter to remove instrument drift."""
    def apply(self, G, args):
        G.Dewow()
        return

class RemoveRinging(FilterCommandBase):
    helpstr = """De-ringing filter

    ringing

    Remove constant horizontal reflectors (e.g. instrument ringing) through
    eigenimage decomposition."""
    cmd = "ringing"
    def apply(self, G, args):
        G.RemoveRinging()
        return

class MigrateFK(FilterCommandBase):
    cmd = "migfk"
    helpstr = """Frequency-wavenumber migration

    migfk [t0_offset]

    Inversion to correct for line-parallel off-nadir reflections. Optionally
    takes `t0_offset`, which is an amount to vertically adjust the radargram to
    align zero-time."""

    def apply(self, G, args):
        if len(args) > 0:
            t0_offset = int(round(float(args[0])))
        else:
            t0_offset = 0
        G.MigrateFK(t0_adjust=t0_offset)
        return

