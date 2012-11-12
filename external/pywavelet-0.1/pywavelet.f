C Routines that can be exposed to python, using f2py
C Author: Arnold Moene (arnold.moene@wur.nl)
C         Wageningen University, The Netherlands
C Date:   December 4, 2002
C Version: 0.1
C $Id: pywavelet.f,v 1.0 2002-12-04 10:19:54+01 arnold Exp $
C Note:   To be used with wavelet software of Torrence and Compo:
C         Wavelet software was provided by C. Torrence and G. Compo,
C         and is available at URL: http://paos.colorado.edu/research/wavelets
C License: GPL (this file and the docs for the exension module only !)


C****************************************************************************
C PYWAVELET: computes the wavelet transform of a time series,
C	     with appropriate parameters.
C
C
C INPUTS:
C
C  n [INT] = the number of points in "y".
C
C  y [DOUBLE PRECISION] = the time series of length "n".
C
C  dt [DOUBLE PRECISION] = amount of time between each Y value, i.e. the sampling time.
C
C  mother [INT] = An integer giving the mother wavelet to use.
C		      0='Morlet'
C		      1='Paul'
C		      2='DOG' (derivative of Gaussian)
C		  If (mother<0 or >2) then default is 'Morlet'.
C
C  param [DOUBLE PRECISION] = mother wavelet parameter. If <0 then default is used.
C	     For 'Morlet' this is k0 (wavenumber), default is 6.
C	     For 'Paul' this is m (order), default is 4.
C	     For 'DOG' this is m (m-th derivative), default is 2.
C
C
C  s0 [DOUBLE PRECISION] = the smallest scale of the wavelet.  Typically = 2*dt.
C	  Note: for accurate reconstruction and variance computation
C	      set s0=dt for Morlet; s0=dt/4 for Paul
C
C  dj [DOUBLE PRECISION] = the spacing between discrete scales. Typically = 0.25.
C	  A smaller # will give better scale resolution, but be slower.
C
C  jtot [INT] = the # of scales.
C		Scales range from s0 up to s0*2^[(jtot-1)*dj],
C		Typically jtot=1+(LOG2(n dt/s0))/dj
C
C  npad [INT] = the total number of points (including padding) to
C		use for the wavelet transform. Typically this is some
C		power of 2. It must be greater or equal to "n".
C		If npad>n, then zeroes are padded onto the end
C		of the time series.
C
C
C OUTPUTS:
C
C  wave [DCMPLX(n,jtot)] = 2D array of the real & imaginary parts
C		    of the wavelet transform, versus time & scale.
C		    CABS(wave) gives the WAVELET amplitude,
C		    ATAN2(AIMAG(wave),DBLE(wave)) gives WAVELET phase.
C		    The WAVELET power spectrum is CABS(wave)**2.
C
C  scale [DOUBLE PRECISION(jtot)] = the wavelet scales that were used.
C
C  period [DOUBLE PRECISION(jtot)] = the "Fourier" periods (in time units) corresponding
C	     to "scale".
C
C  coi [DOUBLE PRECISION(n)] = the e-folding factor used for the cone of influence.
C
C
C REQUIRES:   WAVE_FUNCTION, CFFTPACK
C

      SUBROUTINE PyWAVELET (n,y,dt,mother,param,s0,dj,jtot,npad,
     &	  		    wave,scale,period,coi)
      IMPLICIT none

      INTEGER n,mother,jtot,npad
      DOUBLE PRECISION y(n),dt,param,s0,dj,scale(jtot),period(jtot),
     &	coi(n)
      DOUBLE COMPLEX wave(n,jtot)
      INTEGER i,j,k,nk
      DOUBLE PRECISION ymean,freq1,pi,period1,coi1
Cf2py intent(in)  jtot
Cf2py intent(out) wave
Cf2py intent(out) scale
Cf2py intent(out) period
Cf2py intent(out)  coi
Cf2py depend(jtot) wave
Cf2py depend(n) wave
Cf2py depend(jtot) scale
Cf2py depend(jtot) period
Cf2py depend(n) coi

      CALL WAVELET (n,y,dt,mother,param,s0,dj,jtot,npad,
     &			  wave,scale,period,coi)
      END


C****************************************************************************
C WAVE_FUNCTION: computes the daughter wavelets for a particular
C		 wavelet function, with appropriate parameters.
C
C
C INPUTS:
C
C  nk [INT] = the number of points in "kwave"
C
C  dt [DOUBLE PRECISION] = amount of time between each Y value, i.e. the sampling time.
C
C  mother [INT] = An integer giving the mother wavelet to use.
C		      0='Morlet'
C		      1='Paul'
C		      2='DOG' (derivative of Gaussian)
C
C  param [DOUBLE PRECISION] = mother wavelet parameter. If <0 then default is used.
C	     For 'Morlet' this is k0 (wavenumber), default is 6.
C	     For 'Paul' this is m (order), default is 4.
C	     For 'DOG' this is m (m-th derivative), default is 2.
C
C  scale1 [DOUBLE PRECISION] = the wavelet scale used to construct the daughter.
C
C  kwave [DOUBLE PRECISION(n)] = vector of wavenumbers, used to construct daughter.
C
C
C OUTPUTS:
C
C  period1 [DOUBLE PRECISION] = the "Fourier" period (in time units) that corresponds
C	     to "scale1".
C
C  coi1 [DOUBLE PRECISION] = the e-folding factor used for the cone of influence.
C
C  daughter [DCMPLX(nk)] = real & imaginary parts of the wavelet function
C			  at "scale1" and "kwave".
C
C
C REQUIRES:   FACTORIAL, CHISQR
C
C
C Reference: Tables 1 & 2 in
C	     Torrence, C. and G. P. Compo, 1998: A Practical Guide to
C	     Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
C
C  Copyright (C) 1998, Christopher Torrence and Gilbert P. Compo
C This software may be used, copied, or redistributed as long as it is not
C sold and this copyright notice is reproduced on each copy made.  This
C routine is provided as is without any express or implied warranties
C whatsoever.

      SUBROUTINE pyWAVE_FUNCTION (nk,dt,mother,param,scale1,
     &				kwave,period1,coi1,daughter)
      IMPLICIT none

      INTEGER nk,mother
      DOUBLE PRECISION dt,kwave(nk),param,scale1,period1,coi1
      DOUBLE COMPLEX daughter(nk)

      DOUBLE PRECISION expnt,sk,pi,norm,fourier_factor
      INTEGER k,m,factorial
      DOUBLE PRECISION gamma
Cf2py intent(out) period1
Cf2py intent(out) daughter
Cf2py intent(out) coi1

      CALL WAVE_FUNCTION (nk,dt,mother,param,scale1,
     &	  		kwave,period1,coi1,daughter)

      END





C****************************************************************************
C WAVE_SIGNIF: computes the significance levels for a wavelet transform.
C
C
C INPUTS:
C
C  isigtest [INT] = 0, 1, or 2.
C
C	   If 0, then just do a regular chi-square test,
C		 i.e. Eqn (18) from Torrence & Compo.
C	   If 1, then do a "time-average" test, i.e. Eqn (23).
C		 In this case, DOF(j) should be set to NA, the number
C		 of local wavelet spectra that were averaged together
C		 at each scale. For the Global Wavelet Spectrum,
C		 this would be dof(j)=N-scale(j),
C		 where N is the number of points in your time series.
C	   If 2, then do a "scale-average" test, i.e. Eqns (25)-(28).
C		 In this case, "dof(1)" and "dof(2)" should be set to the
C		 smallest (S1) and largest (S2) scales that were averaged
C		 together, respectively.
C		 e.g. if you scale-averaged scales between 2 and 8,
C		      then dof(1)=2.0 and dof(2)=8.0
C
C
C  n [INT] = the number of points in "y".
C
C  y [DOUBLE PRECISION] = the time series of length "n".
C
C  dt [DOUBLE PRECISION] = amount of time between each Y value, i.e. the sampling time.
C
C  mother [INT] = An integer giving the mother wavelet to use.
C		      0='Morlet'
C		      1='Paul'
C		      2='DOG' (derivative of Gaussian)
C
C  param [DOUBLE PRECISION] = mother wavelet parameter.
C
C  s0 [DOUBLE PRECISION] = the smallest scale of the wavelet.
C
C  dj [DOUBLE PRECISION] = the spacing between discrete scales.
C
C  jtot [INT] = the # of scales.
C
C  scale [DOUBLE PRECISION(jtot)] = the wavelet scales that were used.
C
C  period [DOUBLE PRECISION(jtot)] = the "Fourier" periods corresponding to "scale".
C
C  lag1 [DOUBLE PRECISION] = lag 1 Autocorrelation, used for signif levels.
C		 Default is 0.0, which corresponds to white-noise.
C
C  siglvl [DOUBLE PRECISION] = significance level to use. Default is 0.05 (the "5%" level)
C
C  dof [DOUBLE PRECISION(jtot)] = degrees-of-freedom for signif test.
C     IF SIGTEST=0, then the input dof is ignored.
C     IF SIGTEST=1, then dof(j) = NA, the number of times averaged together.
C     IF SIGTEST=2, then dof(1)=S1, dof(2)=S2, the range of scales averaged.
C
C
C OUTPUTS:
C
C  dof [DOUBLE PRECISION(jtot)] = degrees-of-freedom that were actually used.
C     IF SIGTEST=0, then dof(j) = 2 (or 1 for the 'DOG')
C     IF SIGTEST=1, then dof(j) = degrees-of-freedom versus scale.
C     IF SIGTEST=2, then dof(1)=degrees-of-freedom, dof(2...jtot)=0.0
C
C  fft_theor [DOUBLE PRECISION(jtot)] = theoretical red-noise spectrum vs scale.
C     IF SIGTEST=2, then fft_theor(1) = the average spectrum from S1-->S2
C			 fft_theor(2...jtot) = 0.0
C
C  signif [DOUBLE PRECISION(jtot)] = significance levels vs scale.
C     IF SIGTEST=2, then signif(1) = the significance level
C			 signif(2...jtot) = 0.0
C
C  ymean [DOUBLE PRECISION] = the mean of the time series.
C
C  variance [DOUBLE PRECISION] = the variance of the time series.
C
C  Cdelta [DOUBLE PRECISION] = the constant "Cdelta" for the mother wavelet (Table 2).
C
C  psi0[DOUBLE PRECISION] = the constant 'psi(0)' for the mother wavelet (Table 2)
C
C REQUIRES:   CHISQR
C
C  Copyright (C) 1998, Christopher Torrence and Gilbert P. Compo
C This software may be used, copied, or redistributed as long as it is not
C sold and this copyright notice is reproduced on each copy made.  This
C routine is provided as is without any express or implied warranties
C whatsoever.

      SUBROUTINE pyWAVE_SIGNIF (isigtest,n,y,dt,mother,param,dj,jtot,
     &	     scale,period,lag1,siglvl,dof,fft_theor,signif,
     &	     ymean,variance,Cdelta,psi0)
      IMPLICIT none

      INTEGER isigtest,n,mother,jtot
      DOUBLE PRECISION y(n),dt,param,dj,scale(jtot),period(jtot)
      DOUBLE PRECISION lag1,siglvl,dof(jtot),fft_theor(jtot),
     &	signif(jtot)
      DOUBLE PRECISION ymean,variance,Cdelta,psi0

      INTEGER i,j,m,status,javg1,javg2,navg
      DOUBLE PRECISION pi,freq1,dofmin,gammafac,dj0,Savg,Smid
      DOUBLE PRECISION fft_theor1
      DOUBLE PRECISION chisqr,p,q,bound
Cf2py intent(out) dof
Cf2py intent(out) fft_theor
Cf2py intent(out) signif
Cf2py intent(out) ymean
Cf2py intent(out) variance
Cf2py intent(out) Cdelta
Cf2py intent(out) psi0

C
      CALL WAVE_SIGNIF (isigtest,n,y,dt,mother,param,dj,jtot,
     &	     scale,period,lag1,siglvl,dof,fft_theor,signif,
     &	     ymean,variance,Cdelta,psi0)


      END

