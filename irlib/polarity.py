"""
Functions relevant to the analysis of wavelet polarity in radar traces.
"""

from math import pi
import numpy as np
import scipy.signal as sig


def coherence_power(tr, wavelet):
    """ Return the power obtained by convolving a model wavelet over a trace.
    """
    raise NotImplementedError


def phase_angle(tr, pick, res=pi/16.0, **kwargs):
    """ Calculate the polarization angle in a picked trace as the wave rotation
    angle at which the power is maximized. `kwargs` are passed to
    `polarization_spectrum`. """
    parray = phase_spectrum(tr, **kwargs)
    pickslice = parray[:,pick]
    return np.argmax(pickslice) * res


def phase_spectrum(tr, **kwargs):
    """ Return the polarization angle dividing a model waveform and the
    observed signal. A polarization angle of 0 indicated constructive
    interferences, while a polarization of pi (180d) is perfect destructive
    interference.

    Valid kwargs are:
        wavelength : the target wavelength (integer, in samples)
        res : angular resolution (float, in radians)
    """
    wavelength = int(kwargs.get('wavelength', 10))
    res = kwargs.get('res', pi/16.0)
    n_omega = int(round(2*pi/res))
    ns = len(tr)
    parray = np.empty((n_omega, ns))

    t = np.linspace(0, 2*pi, wavelength) * np.ones((n_omega, wavelength))
    omega = np.arange(0, 2*pi, res)
    mwave = np.sin((t + np.atleast_2d(omega).T))
    pad = wavelength / 2
    parray = sig.convolve(mwave, np.atleast_2d(tr), mode='full')[:,pad:-pad+1]

    return parray





