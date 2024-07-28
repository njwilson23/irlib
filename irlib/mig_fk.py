#! /usr/bin/env python
#
#   Implements an FK (Stolt) migration routine.
#
#   Dmig, tmig, xmig = fkmig(D, dt, dx, v, params)
#
#       D           data array
#       dt          temporal sample rate
#       dx          spatial sample rate
#       v           constant migration velocity
#       params      migration parameters (not yet implemented)
#
#   Code translated from CREWES MATLAB algorithm, and similar restrictions
#   presumably apply. The original license terms are included below.
#
#   BEGIN TERMS OF USE LICENSE
#
#   This SOFTWARE is maintained by the CREWES Project at the Department
#   of Geology and Geophysics of the University of Calgary, Calgary,
#   Alberta, Canada.  The copyright and ownership is jointly held by
#   its author (identified above) and the CREWES Project.  The CREWES
#   project may be contacted via email at:  crewesinfo@crewes.org
#
#   The term 'SOFTWARE' refers to the Matlab source code, translations to
#   any other computer language, or object code
#
#   Terms of use of this SOFTWARE
#
#   1) Use of this SOFTWARE by any for-profit commercial organization is
#      expressly forbidden unless said organization is a CREWES Project
#      Sponsor.
#
#   2) A CREWES Project sponsor may use this SOFTWARE under the terms of the
#      CREWES Project Sponsorship agreement.
#
#   3) A student or employee of a non-profit educational institution may
#      use this SOFTWARE subject to the following terms and conditions:
#      - this SOFTWARE is for teaching or research purposes only.
#      - this SOFTWARE may be distributed to other students or researchers
#        provided that these license terms are included.
#      - reselling the SOFTWARE, or including it or any portion of it, in any
#        software that will be resold is expressly forbidden.
#      - transfering the SOFTWARE in any form to a commercial firm or any
#        other for-profit organization is expressly forbidden.
#
#

from __future__ import print_function
import math
import numpy as np
import pdb, traceback

def csinci():
    """ Complex valued sinc function interpolation.

    trout = csinci(trin, t, tout, sizetable)
    """

def fktran(D, t, x, ntpad=None, nxpad=None, percent=0., ishift=1):
    """ F-K transform using fft on time domain and ifft on space domain. """
    nsamp = D.shape[0]
    ntr = D.shape[1]

    if len(t) != nsamp:
        raise Exception('Time domain length is inconsistent in input')
    if len(x) != ntr:
        raise Exception('Space domain length is inconsistent in input')

    if ntpad is None:
        ntpad = 2.0**nextpow2(t)
    if nxpad is None:
        nxpad = 2.0**nextpow2(x)

    # Get real values of transform with fftrl
    specfx, f = fftrl(D, t, percent, ntpad)

    # Taper and pad in space domain
    if percent > 0.:
        mw = np.tile(mwindow(ntr, percent), (ntr, 1))
        specfx = specfx * mw
    if ntr < nxpad:
        ntr = nxpad                     # this causes ifft to apply the x padding

    spec = np.fft.ifft(specfx.T, n=ntr, axis=0).T
    # Compute kx
    kxnyq = 1. / (2. * (x[1] - x[0]))
    dkx = 2. * kxnyq / ntr
    kx = np.hstack([np.arange(0, kxnyq, dkx), np.arange(-kxnyq, 0, dkx)])

    if ishift:
        tmp = zip(kx, spec)
        tmp.sort()
        kx = [i[0] for i in tmp]
        spec = [i[1] for i in tmp]
    return spec, f, kx


def fftrl(s, t, percent=0.0, n=None):
    """ Returns the real part of the forward Fourier transform. """
    # Determine the number of traces in ensemble
    l = s.shape[0]
    m = s.shape[1]
    ntraces = 1
    itr = 0                             # transpose flag
    if l == 1:
        nsamps = m
        itr = 1
        s = s.T
    elif m == 1:
        nsamps = l
    else:
        nsamps = l
        ntraces = m
    if nsamps != len(t):
        t = t[0] + (t[1] - t[0]) * np.arange(0, nsamps)
    if n is None:
        n = len(t)

    # Apply the taper
    if percent > 0.0:
        mw = np.tile(mwindow(nsamps, percent), (ntraces, 1))
        s = s * mw
    # Pad s if needed
    if nsamps < n:
        s = np.vstack([s, np.zeros([n-nsamps, ntraces])])
        nsamps = n

    # Do the transformation
    spec = np.fft.fft(s, n=nsamps, axis=0)
    spec = spec[:int(n/2)+1, :]              # save only positive frequencies

    # Build the frequency vector
    fnyq = 1. / (2 * (t[1] - t[0]))
    nf = spec.shape[0]
    df = 2.0 * fnyq / n
    f = df * np.arange(0,nf).T
    if itr:
        f = f.T
        spec = spec.T
    return spec, f

def ifktran(spec, f, kx, nfpad=None, nkpad=None, percent=0.0):
    """ Inverse f-k transform.
        Arguments:
            spec    complex valued f-k series
            f       frequency components for rows of spec
            kx      wavenumber components for columns of spec
            nfpad   size to pad spec rows to
            nkpad   size to pad spec columns to
            percent controls cosine taper

        Returns:
            D       2-d array; one trace per column
            t       time coordinates for D
            x       space coordinates for D
    """
    nf,nkx = spec.shape

    if len(f) != nf:
        raise Exception('frequency coordinate vector is wrong size')
    elif len(kx) != nkx:
        raise Exception('wavenumber coordinate vector is wrong size')

    if nfpad is None:
        nfpad = 2.0**nextpow2(len(f))
    if nkpad is None:
        nkpad = 2.0**nextpow2(len(kx))

    # Determine if kx needs to be wrapped
    if kx[0] < 0.0:
        # Looks unwrapped (is this wise?)
        ind = kx >= 0.0
        kx = np.hstack([kx[ind], kx[np.arange(ind[0])]])
        spec = np.hstack([spec[:,ind], spec[:,np.arange(ind[0])]])
    else:
        ind = False

    # Taper and pad in kx
    if percent > 0.0:
        mw = mwindow(nkx, percent)
        if ind.any():
            mw = np.hstack([mw[ind], mw[np.arange(ind[0])]])
        mw = mw.repeat(nkz, axis=0)
        spec = spec * mw
    if nkx < nkpad:
        nkx = nkpad

    # Performs the transforms
    specfx = np.fft.fft(spec, nkx)
    D, t = ifftrl(specfx, f)

    # Compute x
    dkx = kx[1] - kx[0]
    xmax = 1.0 / dkx
    dx = xmax / nkx
    x = np.arange(0, xmax, dx)
    return D, t, x

def ifftrl(spec, f):
    """ Inverse Fourier transform for real-valued series.
        Arguments:
            spec    input spectrum
            f       input frequency coordinates
        Returns:
            r       output trace
            t       output time vector
    """
    m,n = spec.shape            # Will be a problem if spec is 1-dimensional
    itr = 0
    if (m == 1) or (n == 1):
        if m == 1:
            spec = spec.T
            itr = 1
        nsamp = len(spec)
        ntr = 1
    else:
        nsamp = m
        ntr = n

    # Form the conjugate symmetric complex spectrum expected by ifft
    # Test for nyquist
    nyq = 0
    if (spec[-1] == np.real(spec[-1])).all():
        nyq = 1
    if nyq:
        L1 = np.arange(nsamp)
        L2 = L1[-2:0:-1]
    else:
        L1 = np.arange(nsamp)
        L2 = L1[-2:0:-1]            # WTF? -njw
    symspec = np.vstack([spec[L1,:], np.conj(spec[L2,:])])
    # Transform the array
    r = (np.fft.ifft(symspec.T)).real.T
    # Build the time vector
    n = len(r)
    df = f[1] - f[0]
    dt = 1.0 / (n*df)
    t = dt * np.arange(n).T
    if itr == 1:
        r = r.T
        t = t.T
    return r, t

def mwindow(n, percent=10.):
    """ Creates a boxcar window with raised-cosine tapers. """
    if type(n) is not int and type(n) is not float:
        n = len(n)
    # Compute the hanning function
    if percent > 50. or percent < 0.:
        raise Exception('Invalid percent in function mwindow (={0})'.format(percent))
    m = 2.0 * math.floor(percent * n / 100.)
    h = np.hanning(m)
    return np.hstack([h[:m/2], np.ones([n-m]), h[m/2:]])

def mwhalf(n, percent=10.):
    """ Half mwindow. """
    if type(n) is not int and type(n) is not float:
        n = len(n)
    # Compute the hanning function
    if percent > 100. or percent < 0.:
        raise Exception('Invalid percent in function mwhalf (={0})'.format(percent))
    m = int(math.floor(percent * n / 100.))
    h = np.hanning(2*m)
    return np.hstack([np.ones([n-m]), h[m:0:-1]])

def nextpow2(a):
    """ Gives the next power of 2 larger than a. """
    return np.ceil(np.log(a) / np.log(2)).astype(int)

def fkmig(D, dt, dx, v, params=None):

    nsamp = D.shape[0]
    ntr = D.shape[1]
    t = np.arange(0, nsamp) * dt
    x = np.arange(0, ntr) * dx
    interpolated = True

    fnyq = 1.0 / (2.0*dt)
    knyq = 1.0 / (2.0*dx)
    tmax = t[-1]
    xmax = abs(x[-1]-x[0])

    # Deal with parameters
    if params == None:
        fmax = 0.6 * fnyq
        fwid = 0.2 * (fnyq - fmax)
        dipmax = 85.0
        dipwid = 90.0 - dipmax
        tpad = min([0.5 * tmax, abs(tmax / math.cos(math.pi*dipmax / 180.0))])
        xpad = min([0.5 * xmax, xmax / math.sin(math.pi*dipmax / 180.0)])
        padflag = 1
        intflag = 3
        cosflag = 1
        lsinc = 1
        ntable = 25
        mcflag = 0      # Faster, less memory-efficient transform (not implemented)
        kpflag = 50.0

    # Apply padding
    # tpad
    nsampnew = int(2.0**nextpow2( round((tmax+tpad) / dt + 1.0) ))
    tmaxnew = (nsampnew-1)*dt
    tnew = np.arange(t[0], tmaxnew+dt, dt)
    ntpad = nsampnew-nsamp
    D = np.vstack([D,np.zeros([ntpad,ntr])])

    # xpad
    ntrnew = 2 ** nextpow2(round((xmax+xpad) / dx + 1))
    xmaxnew = (ntrnew-1)*dx + x[0]
    xnew = np.arange(x[0], xmaxnew+dx, dx)
    nxpad = ntrnew-ntr
    D = np.hstack([D, np.zeros([nsampnew,nxpad])])

    # Forward f-k transform
    fkspec, f, kx = fktran(D, tnew, xnew, nsampnew, ntrnew, 0, 0)
    df = f[1] - f[0]
    nf = len(f)

    # Compute frequency mask
    ifmaxmig = int(round((fmax+fwid) / df + 1.0))
    pct = 100.0 * (fwid / (fmax+fwid))
    fmask = np.hstack([mwhalf(ifmaxmig,pct), np.zeros([nf-ifmaxmig])])
    fmaxmig = (ifmaxmig-1)*df       # i.e. fmax+fwid to nearest sample

    # Now loop over wavenumbers
    ve = v / 2.0                    # exploding reflector velocity
    dkz = df / ve
    kz = (np.arange(0,len(f)) * dkz).T
    kz2 = kz ** 2

    th1 = dipmax * math.pi / 180.0
    th2 = (dipmax+dipwid) * math.pi / 180.0
    if th1 == th2:
        print("No dip filtering")

    for j,kxi in enumerate(kx):
        # Evanescent cut-off
        fmin = abs(kxi) * ve
        ifmin = int(math.ceil(fmin / df)) + 1

        # Compute dip mask
        if th1 != th2:
            # First physical frequency excluding dc
            ifbeg = max([ifmin, 1])+1
            # Frequencies to migrate
            ifuse = np.arange(ifbeg, ifmaxmig+1)
            if len(ifuse) == 1:
                # Special case
                dipmask = np.zeros(f.shape)
                dipmask[ifuse-1] = 1
            else:
                # Physical dips for each frequency
                theta = np.arcsin(fmin / f[ifuse])
                # Sample number to begin ramp
                if1 = round(fmin / (math.sin(th1) * df))
                if1 = max([if1, ifbeg])
                # sample number to end ramp
                if2 = round(fmin / (math.sin(th2) * df))
                if2 = max([if2, ifbeg])
                # Initialize mask to zeros
                dipmask = np.zeros(f.shape)
                # Pass these dips
                dipmask[if1:nf-1] = 1
                dipmask[if2:if1] = 0.5 + 0.5 * np.cos(
                        (theta[np.arange(if2, if1, -1) - ifbeg] - th1)
                        * math.pi / float(th2-th1))
        else:
            dipmask = np.ones(f.shape)

        # Apply masks
        tmp = fkspec[:, j] * fmask * dipmask

        # Compute f that map to kz
        fmap = ve * np.sqrt(kx[j]**2 + kz2)
        # Contains one value for each kz giving the frequency
        # that maps there to migrate the data
        # Many of these frequencies will be far too high
        ind = np.vstack(np.nonzero(fmap <= fmaxmig)).T
        # ind is an array of indicies of fmap which will always start at 1
        # and end at the highest f to be migrated

        # Now map samples by interpolation
        fkspec[:, j] *= 0.0             # initialize output spectrum to zero
        if len(ind) != 0:
            # Compute cosine scale factor
            if cosflag:
                if fmap[ind].all() == 0:
                    scl = np.ones(ind.shape[0])
                    li = ind.shape[0]
                    scl[1:li] = (ve * kz[ind[1:li]] / fmap[ind[1:li]])[:,0]
                else:
                    scl = ve * kz[ind] / fmap[ind]
            else:
                scl = np.ones(ind.shape[0])
            if intflag == 0:
                # Nearest neighbour interpolation
                ifmap = (fmap[ind] / df).astype(int)
                fkspec[ind, j] = (scl.squeeze() \
                    * tmp[ifmap.squeeze()]).reshape([-1,1])
            elif intflag == 1:
                # Complex sinc interpolation
                fkspec[ind, j] = scl \
                        * csinci(tmp, f, fmap[ind], np.hstack([lsinc,ntable]))
            elif intflag == 2:
                # Spline interpolation
                # Not implemented
                pass
            elif intflag == 3:
                # Linear interpolation
                r_interp = scl.squeeze() \
                    * np.interp(fmap[ind], f, tmp.real).squeeze()
                j_interp = scl.squeeze() \
                    * np.interp(fmap[ind], f, tmp.imag).squeeze()
                fkspec[ind, j] = (r_interp + j_interp * 1j).reshape([-1,1])

    # Inverse transform
    Dmig, tmig, xmig = ifktran(fkspec, f, kx)

    # Remove padding, if desired
    if padflag:
        Dmig = Dmig[:nsamp, :ntr]
        tmig = tmig[:nsamp]
        xmig = xmig[:ntr]

    return Dmig, tmig, xmig
