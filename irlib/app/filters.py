""" This file defines the filter commands available from irlib-based apps. """

class Command(object):
    """ A Command is implemented as a class with a command-line signature and
    an *apply* method that takes an appropriate Gather object as an argument.
    """

    helpstr = "No help"

    def __init__(self):
        return

    def apply(self, G):
        """ Overload this method to perform operations on Gather object *G*. """
        raise Exception("apply() is undefined in Command baseclass")

class LinearGainControl(Command):
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

class AutoGainControl(Command):
    cmd = "agc"
    helpstr = """Automatic gain control

    agc

    Apply automatic gain control, normalizing the power in a finite window to a
    constant value."""

    def apply(self, G, args):
        G.DoAutoGainControl(5e-8)
        return

class ReflectionPower(Command):
    cmd = "power"
    def apply(self, G, args):
        G.data = G.data**2
        return

class Lowpass_FD(Command):
    cmd = "lowpass"
    helpstr = """Frequency-domain lowpass

    lowpass

    Apply a frequency domain lowpass filter implemented using a windowed sinc
    kernel."""
    def apply(self, G, args):
        G. DoWindowedSinc(cutoff=25e6, bandwidth=5e6, mode="lowpass")
        return

class Highpass_FD(Command):
    cmd = "highpass"
    helpstr = """Frequency-domain highpass

    highpass

    Apply a frequency domain highpass filter implemented using a windowed sinc
    kernel."""
    def apply(self, G, args):
        G.DoWindowedSinc(cutoff=25.e6, bandwidth=5.e6, mode='highpass')
        return

class Lowpass_TD(Command):
    cmd = "lowpass_td"
    helpstr = """Time-domain lowpass

    lowpass

    Apply a frequency domain lowpass filter implemented as a
    moving average with a blackman window."""
    def apply(self, G, args):
        G. DoMoveAvg(21, kind="blackman", mode="lowpass")
        return

class Highpass_TD(Command):
    cmd = "highpass_td"
    helpstr = """Time-domain highpass

    highpass

    Apply a frequency domain highpass filter implemented as a
    spectrally-inverted moving average with a blackman window."""
    def apply(self, G, args):
        G.DoMoveAvg(7, kind="blackman", mode='highpass')
        return

class Dewow(Command):
    cmd = "dewow"
    helpstr = """Dewow

    dewow

    Apply a "dewowing" filter to remove instrument drift."""
    def apply(self, G, args):
        G.Dewow()
        return

class RemoveRinging(Command):
    helpstr = """De-ringing filter

    ringing

    Remove constant horizontal reflectors (e.g. instrument ringing) through
    eigenimage decomposition."""
    cmd = "ringing"
    def apply(self, G, args):
        G.RemoveRinging()
        return

class MigrateFK(Command):
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

