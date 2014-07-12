#! /usr/bin/python
#
#   An implementation of Kirchoff migration
#

import math
import numpy as np

def calc_toffset(X, T, xoffset=12.0, v=1.68e8):
    """ Perform normal move-out and calculate a diffraction hyperbola.
        X           vector of horizontal offsets
        T           vector of time offsets
    Return an array of one-sided vertical offsets in terms of time.
    """
    # Define arrays of time and distance as starting points
    TT = T.reshape([-1,1]).repeat(len(X), axis=1)
    XX = X.reshape([1,-1]).repeat(len(T), axis=0)

    # Perform normal moveout for the given xoffset to get depth
    TT = np.ma.masked_less(TT, math.sqrt(xoffset**2 / v**2))
    DD = np.sqrt(TT**2 * v**2 / 4 - xoffset**2 / 4)

    # Compute the diffraction hyperbolas in terms of time
    DHYP = np.sqrt(DD**2 + XX**2) - DD

    # Reverse moveout to get back into time domain
    THYP = 2 * np.sqrt(DHYP**2 + xoffset**2 / 4) / v

    return THYP

def warp_array_linear(A, T, i):
    """ Warp an array A by shifting the values up by elements in array T. """
    if A.shape != T.shape:
        raise Exception('warp_array: A and T must have identical dimensions')
        return None
    x = np.arange(A.shape[1])
    y = np.arange(A.shape[0])
    X, Y = np.meshgrid(x, y)
    Y_warp = Y + T
    W = np.zeros(A.shape)

    # Interpolate the warped values and put into W
    for icol in range(A.shape[1]):
        W[:,icol] = np.interp(Y_warp[:,icol], Y[:,icol], A[:,icol])
    return W

def mirror(arr):
    return np.stack([arr[:0:-1], arr])

def mig_kirchoff(D, X, T, v=1.68e8, xoffset=12.0, xwindow=50.0, pad_spacing=1.0):
    """ Perform Kerchoff migration on a 2D line.
    Requires
        D           2D array of traces
        X           vector of horizontal offsets
        T           vector of time offsets
    Optionally also takes
        v           wave velocity (m/s, assumed to be for ice)
        xoffset     transmitter-reciever spacing (m)
        xwindow     half the horizontal window to be considered for migration
                    larger is slower, but important when ice is deep

    Kirchoff migration works by summing all of the wavelets along a
    pre-calculated hyperbola that come from the centre reflector. This
    performed in a very long loop on a trace-by-trace basis through the entire
    line. Compared the Stolt (FK) migration, this has the advantage of being
    amenable to irregular horizontal spacing.
    """
    # Create a set of horizontal coordinates
    X_hyp_half = np.arange(0, xwindow, 1.0)
    X_hyp = np.hstack([-X_hyp_half[:0:-1], X_hyp_half])

    # Create a look-up table of theoretical hyperbolae for all travel times
    THYP = calc_toffset(X_hyp_half, T, xoffset=xoffset, v=v) \
            / np.hstack([1, np.diff(T)]).reshape([-1,1]).repeat(int(xwindow),
            axis=1)
    THYP_sym = np.hstack([THYP[:,:0:-1], THYP])

    # Pad the data with zeros
    halfwin = math.ceil(xwindow / pad_spacing)
    pad = np.zeros([len(T), halfwin])
    DPAD = np.hstack([pad, D, pad])
    XPAD = np.hstack([np.arange(-xwindow, 0, pad_spacing),
                        X, np.arange(X[-1], X[-1]+xwindow, pad_spacing)])

    # Loop through traces in groups of size (xwindow*2)+1
    DMIG = np.zeros(D.shape)
    for i in range(D.shape[1]):

        if i % 10 == 0: print 'calculating trace ', i

        ipad = i + halfwin

        # Modify lookup table to handle irregular spacing
        THYP_mod = np.zeros(THYP_sym.shape)
        #for irow in range(THYP_sym.shape[0]):
            #THYP_mod[irow,:] = np.interp(XPAD[ipad-halfwin+1:ipad+halfwin], X_hyp, THYP_sym[irow,:])
        THYP_mod = THYP_sym

        # Build a warped trace group with linear interpolation
        TR_GRP = warp_array_linear(DPAD[:,ipad-halfwin:ipad+halfwin-1], THYP_sym, i)

        # For every point on the centre trace, sum the wavelet contributions
        DMIG[:,i] = TR_GRP.sum(axis=1)

    # Smear the vertical

    return DMIG

def test1():
    # Test function

    X = np.arange(0.0, 300.0, 4.0)
    dt = 1e-8
    T = np.arange(0.0, 256*dt, dt)
    THYP = calc_toffset(X, T)

    print THYP.max(), THYP.min(), THYP.mean()

    plt.figure()
    plt.imshow(THYP)
    plt.axis('tight')
    plt.show()

def test2():
    import irlib

    INFILE = '../field/radar_200905/Glacier2/Glacier2_May09_utm.h5'
    S = irlib.Survey(INFILE)
    L = S.ExtractLine(1)
    L.DoTimeGainControl()
    D = L.data.copy()
    eastings, northings, p = L.LineProject()
    P = np.sqrt((eastings-eastings[0])**2 + (northings-northings[0])**2)
    T = np.arange(L.data.shape[0]) * 1e-8
    Dmig = mig_kirchoff(D, P, T, xoffset=0.0, xwindow=10.0)

    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(D, cmap='gray')
    plt.subplot(1,2,2)
    plt.imshow(Dmig, cmap='gray')
    plt.show()

def test3():
    # Test on synthetic planes
    D = np.zeros([256, 256])
    D[63,5:251] = 1
    D[127,5:251] = 1
    D[191,5:251] = 1
    P = np.arange(0, 512, 2)
    T = np.arange(256) * 1e-8
    Dmig = mig_kirchoff(D, P, T, xoffset=15.0)

    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(D, cmap='gray')
    plt.subplot(1,2,2)
    plt.imshow(Dmig, cmap='gray')
    plt.show()

def test4():
    # Test on synthetic slopes
    D = np.zeros([256, 256])
    for r in range(5,251):
        for c in range(60,81):
            if r == round(0.1*c + 60):
                D[r,c] = 1
    for r in range(5,251):
        for c in range(160,181):
            if r == round(0.1*c + 160):
                D[r,c] = 1
    for r in range(5,251):
        for c in range(90,150):
            if r == round(-0.3*c + 150):
                D[r,c] = 1
    for r in range(5,251):
        for c in range(150,250):
            if r == round(-0.7*c + 200):
                D[r,c] = 1
    P = np.arange(0, 512, 2)
    T = np.arange(256) * 1e-8
    Dmig = mig_kirchoff(D, P, T, xoffset=15.0)

    plt.figure()
    plt.subplot(1,2,1)
    plt.imshow(D, cmap='gray')
    plt.subplot(1,2,2)
    plt.imshow(Dmig, cmap='gray')
    plt.show()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test1()
    test2()
    test3()
    test4()






