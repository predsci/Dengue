      subroutine nserovechost(ns, nweeks, Nh, muH, muV, 
     $     sigmaH, sigmaV,gamma, beta, alpha, rho, k, Va, 
     $     eta, delta, chi, omega,fracS0, sEp, sIp,
     $     Ve, Vi, iseed, rtn)

      implicit none

      integer iday_per_week, day_per_year
      parameter (iday_per_week=7)
      parameter (day_per_year = 365)
! If this is a subroutine
      integer nstep
      parameter(nstep = 50)
      integer ns, nweeks, iseed

      real*8 Nh, muH, sigmaH, muV, sigmaV, k, Va
      real*8 alpha(ns), gamma(ns), eta(ns), delta(ns)
      real*8 rho(ns), chi(ns), omega(ns), beta(ns)
      real*8 dsdt(ns, (nweeks+2)*nstep*iday_per_week)
      real*8 fracS0, sEp(ns), sIp(ns), Ve(ns), Vi(ns)
      real*8 dt,tps(nweeks) !time-series
      real*8 pC, e_nondengue
      real*8 rtn(ns, nweeks)
      integer i

      pC = 1.0d0
      e_nondengue = 0.0d0
 
      dt = 1.0d0 / dble(nstep - 1)

! Comment this section if this is a program

      do i = 1, nweeks
         tps(i) = dble(i)
      enddo


      dsdt = 0.0d0

      call RK4OneD(ns, Nh, muH, sigmaH, muV, sigmaV, k, Va,
     $     alpha, gamma, eta, delta, rho, chi, omega, beta,
     $     sEp, sIp, Ve, Vi,nweeks,nstep,dt,tps,dsdt, fracS0)


      call weekly1D(ns, nweeks, nstep, dsdt, pC,
     $     e_nondengue, rtn)
 

c$$$      do i = 1, nweeks
c$$$         write(6,999) tps(i), rtn(1:ns,i)
c$$$      enddo
c$$$ 999  format(f10.1,1x,10(f12.5,1x))

!      return
      end subroutine nserovechost


!
! -------------------------------------------------------------------
!


      subroutine RK4OneD(ns, Nh, muH, sigmaH, muV, sigmaV, k,
     $     Va, alpha, gamma, eta, delta, rho, chi, omega, beta,
     $     sEp, sIp, Ve, Vi, nweeks,nstep,dt,tps,dsdt,fracS0)


      implicit none
      integer iday_per_week
      parameter (iday_per_week=7)
      integer iseed
      integer nweeks, nstep
      integer ns
      real*8 nH, muH, sigmaH, muV, sigmaV, k, Va
      real*8 alpha(ns), gamma(ns), eta(ns), delta(ns)
      real*8 rho(ns), chi(ns), omega(ns), beta(ns)
      real*8 dsdt(ns, (nweeks+2)*nstep*iday_per_week)
      real*8 fracS0
      real*8 dt,tps(nweeks) !time-series
      real*8 tps2(0:(nweeks+1))
      real*8 t_cur, p5, p3, p6
      real*8 ran1
      integer i, iweek, istep, iday
      integer icount
      logical debug

! States/Compartments
      real*8 S0, Vs0, sRtot
      real*8 sEp(ns), sIp(ns), sC(ns), sA(ns), sR(ns)
      real*8 sEs(ns), sIs(ns), sRTN(ns)
      real*8 Ve(ns), Vi(ns)
! tmp States/Compartments
      real*8 tmpS0, tmpVs0, tmpRtot
      real*8 tmpEp(ns), tmpIp(ns), tmpC(ns), tmpA(ns), tmpR(ns)
      real*8 tmpEs(ns), tmpIs(ns), tmpRTN(ns)
      real*8 tmpVe(ns), tmpVi(ns)

! Derivative of States/Compartments
      real*8 dS01, dVs01, dRtot1
      real*8 dEp1(ns), dIp1(ns), dC1(ns), dA1(ns), dR1(ns)
      real*8 dEs1(ns), dIs1(ns), dRTN1(ns)
      real*8 dVe1(ns), dVi1(ns)

      real*8 dS02, dVs02, dRtot2
      real*8 dEp2(ns), dIp2(ns), dC2(ns), dA2(ns), dR2(ns)
      real*8 dEs2(ns), dIs2(ns), dRTN2(ns)
      real*8 dVe2(ns), dVi2(ns)

      real*8 dS03, dVs03, dRtot3
      real*8 dEp3(ns), dIp3(ns), dC3(ns), dA3(ns), dR3(ns)
      real*8 dEs3(ns), dIs3(ns), dRTN3(ns)
      real*8 dVe3(ns), dVi3(ns)

      real*8 dS04, dVs04, dRtot4
      real*8 dEp4(ns), dIp4(ns), dC4(ns), dA4(ns), dR4(ns)
      real*8 dEs4(ns), dIs4(ns), dRTN4(ns)
      real*8 dVe4(ns), dVi4(ns)



      debug = .false.

      p5 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0
!
! ran1 works with -iseed
!
      iseed = -abs(iseed)


! initialize the vectors

      S0 = fracS0 * Nh

      sC  = 0.0d0
      sA  = 0.0d0
      sR  = 0.0d0
      sEs = 0.0d0
      sIs = 0.0d0
      sRTN = 0.0d0

      sRtot = Nh - (s0 + sum(sEp + sIp + sEs + sIs + sC + sA + sR))
     
      Vs0 = k * Nh - sum(Ve+Vi)

      dsdt = 0.0d0

      tps2 = 0.0d0
      tps2(1:nweeks) = tps
      tps2(0) = tps2(1) + (tps(2) - tps(1)) 
      tps2((nweeks+1)) = tps2(nweeks)  + (tps(2) - tps(1)) 

      t_cur = tps2(0)

      icount = 0

      do iweek=0,nweeks+1

         if (debug) Then
            write(6,999) tps2(iweek), S0 + sRtot 
     $           + sum(sEp + sEs +sIp + sIs + sC + sA + sR)
     $           ,Vs0+sum(Ve+Vi)
         endif
 999        format(f12.1,1x,10(f15.4,1x))
! Eventually we want to change this to work on days and then integrate each day-done

         do iday=1,iday_per_week
            do istep=1,nstep
               t_cur = t_cur+dt

               call diff1D(ns, Nh, muH, sigmaH, muV, sigmaV, k, Va, 
     $              alpha, gamma, eta, delta, rho, chi, omega, beta,
     $              S0, sEp, sIp, sC, sA, sR, sEs, sIs, sRtot,
     $              sRTN,vS0, Ve, Vi,
     $              dS01, dEp1, dIp1, dC1, dA1, dR1, dEs1, dIs1, dRtot1,
     $              dRTN1, dVs01, dVe1, dVi1, t_cur)


               tmpS0 = S0  + dt * dS01 * p5
               tmpEp  = sEp  + dt * dEp1 * p5
               tmpIp = sIp + dt * dIp1 * p5
               tmpC  = sC  + dt * dC1 * p5
               tmpA  = sA  + dt * dA1 * p5
               tmpR  = sR  + dt * dR1 * p5
               tmpEs = sEs + dt * dEs1 * p5
               tmpIs = sIs + dt * dIs1 * p5
               tmpRtot = sRtot + dt * dRtot1 * p5
               tmpRTN = sRTN + dt * dRTN1 * p5
               tmpVs0 = Vs0 + dt * dVs01 * p5
               tmpVe  = Ve + dt * dVe1 * p5
               tmpVi  = Vi + dt * dVi1 * p5


               call diff1D(ns, Nh, muH, sigmaH, muV, sigmaV, k, Va, 
     $              alpha, gamma, eta, delta, rho, chi, omega, beta,
     $              tmpS0, tmpEp, tmpIp, tmpC, tmpA, tmpR, tmpEs, tmpIs,
     $              tmpRtot, tmpRTN, tmpvS0, tmpVe, tmpVi,
     $              dS02, dEp2, dIp2, dC2, dA2, dR2, dEs2, dIs2, dRtot2,
     $              dRTN2, dVs02, dVe2, dVi2, t_cur)

               tmpS0 = S0 + dt * dS02 * p5
               tmpEp  = sEp + dt * dEp2 * p5
               tmpIp  = sIp + dt * dIp2 * p5
               tmpC  = sC + dt * dC2 * p5
               tmpA  = sA + dt * dA2 * p5
               tmpR  = sR + dt * dR2 * p5
               tmpEs = sEs + dt * dEs2 * p5
               tmpIs = sIs + dt * dIs2 * p5
               tmpRtot = sRtot + dt * dRtot2 * p5
               tmpRTN = sRTN + dt * dRTN2 * p5
               tmpVs0 = Vs0 + dt * dVs02 * p5
               tmpVe  = Ve + dt * dVe2 * p5
               tmpVi  = Vi + dt * dVi2 * p5

               call diff1D(ns, Nh, muH, sigmaH, muV, sigmaV, k, Va, 
     $              alpha, gamma, eta, delta, rho, chi, omega, beta,
     $              tmpS0, tmpEp, tmpIp, tmpC, tmpA, tmpR, tmpEs, tmpIs,
     $              tmpRtot, tmpRTN, tmpvS0, tmpVe, tmpVi,
     $              dS03, dEp3, dIp3, dC3, dA3, dR3, dEs3, dIs3, dRtot3,
     $              dRTN3, dVs03, dVe3, dVi3, t_cur)


               tmpS0 = S0 + dt * dS03 * p5
               tmpEp  = sEp + dt * dEp3 * p5
               tmpIp  = sIp + dt * dIp3 * p5
               tmpC  = sC + dt * dC3 * p5
               tmpA  = sA + dt * dA3 * p5
               tmpR  = sR + dt * dR3 * p5
               tmpEs = sEs + dt * dEs3 * p5
               tmpIs = sIs + dt * dIs3 * p5
               tmpRtot = sRtot + dt * dRtot3 * p5
               tmpRTN = sRTN + dt * dRTN3 * p5
               tmpVs0 = Vs0 + dt * dVs03 * p5
               tmpVe  = Ve + dt * dVe3 * p5
               tmpVi  = Vi + dt * dVi3 * p5


               call diff1D(ns, Nh, muH, sigmaH, muV, sigmaV, k, Va, 
     $              alpha, gamma, eta, delta, rho, chi, omega, beta,
     $              tmpS0, tmpEp, tmpIp, tmpC, tmpA, tmpR, tmpEs, tmpIs,
     $              tmpRtot, tmpRTN, tmpvS0, tmpVe, tmpVi,
     $              dS04, dEp4, dIp4, dC4, dA4, dR4, dEs4, dIs4, dRtot4,
     $              dRTN4, dVs04, dVe4, dVi4, t_cur)


               S0 = S0 + 
     $              dt *(dS01 * p6 + dS02 * p3 + dS03 * p3 + dS04 * p6)

               sEp = sEp + 
     $              dt *(dEp1 * p6 + dEp2 * p3 + dEp3 * p3 + dEp4 * p6)

               sIp = sIp + 
     $              dt *(dIp1 * p6 + dIp2 * p3 + dIp3 * p3 + dIp4 * p6)

               sC = sC + 
     $              dt *(dC1 * p6 + dC2 * p3 + dC3 * p3 + dC4 * p6)

               sA = sA + 
     $              dt *(dA1 * p6 + dA2 * p3 + dA3 * p3 + dA4 * p6)

               sR = sR + 
     $              dt *(dR1 * p6 + dR2 * p3 + dR3 * p3 + dR4 * p6)

               sEs = sEs + 
     $              dt *(dEs1 * p6 + dEs2 * p3 + dEs3 * p3 + dEs4 * p6)

               sIs = sIs + 
     $              dt *(dIs1 * p6 + dIs2 * p3 + dIs3 * p3 + dIs4 * p6)

               sRtot = sRtot + dt *(dRtot1 * p6 + dRtot2 * p3 +
     $              dRtot3 * p3 + dRtot4 * p6)


               Vs0 = vS0 + dt *(dvS01 * p6 + dvS02 * p3 + 
     $              dvS03 * p3 + dvS04 * p6)


               Ve = Ve + 
     $              dt *(dVe1 * p6 + dVe2 * p3 + dVe3 * p3 + dVe4 * p6)

               Vi = Vi + 
     $              dt *(dVi1 * p6 + dVi2 * p3 + dVi3 * p3 + dVi4 * p6)

               sRTN = sRTN + 
     $              dt *(dRTN1 * p6 + dRTN2 * p3 + dRTN3 * p3 + 
     $              dRTN4 * p6)
!     Update -dS/dt       
               icount= icount + 1
               
               dsdt(1:ns,icount) = sRTN

            enddo ! End of a single day         
         enddo                  ! End of days_per_week
         
      enddo                     ! End of loop on weeks

      return
      end subroutine RK4OneD

!
! -------------------------------------------------------------------
!

!! List of states/compartments that we are integrating
!! Human: S0, sEp, sIp, sC, sA, sR, sEs, sIs, sRtot
!! Vector: Vs0 Ve Vi
!! Most are vectors of length ns
!! These are numbers: S0, Rtot and Vs0

!! Parameters
!! Numbers: NH, muH, sigmaH, muV, sigmaV, k, Va
!! vectors: alpha, gamma, eta, delta, rho, chi, omega, beta


      SUBROUTINE Diff1D(ns, Nh, muH, sigmaH, muV, sigmaV, k, Va, 
     $     alpha, gamma, eta, delta, rho, chi, omega, beta,
     $     S0, sEp, sIp, sC, sA, sR, sEs, sIs, sRtot,
     $     sRTN, vS0, Ve, Vi,
     $     dS0, dEp, dIp, dC, dA, dR, dEs, dIs, dRtot, dRTN,
     $     dVs0, dVe, dVi, time)

      implicit none

      integer ns
      real*8 nH, muH, sigmaH, muV, sigmaV, k, Va
      real*8 alpha(ns), gamma(ns), eta(ns), delta(ns)
      real*8 rho(ns), chi(ns), omega(ns), beta(ns)
! States/Compartments
      real*8 S0, Vs0, sRtot
      real*8 sEp(ns), sIp(ns), sC(ns), sA(ns), sR(ns)
      real*8 sEs(ns), sIs(ns), sRTN(ns)
      real*8 Ve(ns), Vi(ns)
! Derivative of States/Compartments
      real*8 dS0, dVs0, dRtot
      real*8 dEp(ns), dIp(ns), dC(ns), dA(ns), dR(ns)
      real*8 dEs(ns), dIs(ns), dRTN(ns)
      real*8 dVe(ns), dVi(ns)
! Current time
      real*8 time
! Auxilary Matrix 
      real*8 fMat(ns,ns), fvec(ns)
! two pi and vector birth rate
      real*8 twopi, scale, vBirth, Nv
      integer j

      fvec = 1.0
      fMat = 1.0d0
      do j = 1, ns
         fMat(j,j) = 0.0d0
      enddo


      twopi = 2.0d0 * acos(-1.0d0)

! This gives the birth rate a periodicity of one year
! We will need to add a phase 

      scale = twopi / 365.0d0
      
      vBirth = k * Nh * (1.0d0 - Va * cos(scale*time))
      Nv = k * Nh

! Integrate the equations 

      dS0 = (Nh - S0) * muH - S0 /Nh * sum(alpha * Vi)

      dEp = S0 / Nh * (alpha * Vi) - (sigmaH + muH) * sEp

      dIp = sigmaH * sEp - gamma * sIp - muH * sIp

      dC = gamma * sIp 
     $     - sC/Nh * matmul(fMat, eta * alpha * Vi) 
     $     - delta * sC - muH * sC

      dA = delta * sC * (fvec - rho) 
     $     - sA / Nh * matmul(fMat, chi * alpha * Vi)
     $     - omega * sA - muH * sA

      dR = omega * sA 
     $     - sR / Nh * matmul(fMat, alpha * Vi)
     $     - muH * sR

      dEs = alpha * Vi / Nh * (
     $     + eta * matmul(fMat, sC)
     $     + chi * matmul(fMat, sA) + matmul(fMat, sR)) 
     $     - (sigmaH + muH) * sEs

      dIs = sigmaH * sEs - gamma * sIs - muH * sIs

      dRtot = sum(gamma * sIs) - muH * sRtot

      dVs0 = (vBirth - Vs0) * muV - Vs0 / Nh * sum( beta * (sIp + sIs))

      dVe = Vs0 / nH * beta * (sIp + sIs) - (sigmaV + muV) * Ve

      dVi = sigmaV * Ve - muV * Vi 

!! Cumulative infectious - coming from primary and secondary exposed
      dRTN = sigmaH * (sEp + sEs)

      RETURN
      END subroutine Diff1D


!
! -------------------------------------------------------------------
!

           
        subroutine weekly1D(ns,iweek,nstep,dsdt,pC,e_nondengue,x)
        
        implicit none
	integer iday_per_week
        parameter (iday_per_week=7)  
        integer ns, iweek,nstep
        real*8 dsdt(ns,iweek*nstep*iday_per_week)
	real*8 x(ns, iweek), pC, e_nondengue
        integer ishift, i, j

        ishift = nstep * iday_per_week * 0.50d0

        do j = 1, ns
           do i=1,iweek
c$$$           x(i) = dsdt(i*nstep*iday_per_week) - 
c$$$     $     dsdt(1+(i-1)*nstep*iday_per_week)

              x(j,i) = dsdt(j,i*nstep*iday_per_week+ishift) - 
     $             dsdt(j,i*nstep*iday_per_week-ishift)

           enddo

        enddo
        x= x* pC +e_nondengue
        return
        end subroutine weekly1D


!
! -------------------------------------------------------------------
!



      FUNCTION ran1(idum)  
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV  
      REAL*8 ran1,AM,EPS,RNMX  
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,  
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)  
      INTEGER j,k,iv(NTAB),iy  
      SAVE iv,iy  
      DATA iv /NTAB*0/, iy /0/  
      if (idum.le.0.or.iy.eq.0) then  
        idum=max(-idum,1)  
        do 11 j=NTAB+8,1,-1  
          k=idum/IQ  
          idum=IA*(idum-k*IQ)-IR*k  
          if (idum.lt.0) idum=idum+IM  
          if (j.le.NTAB) iv(j)=idum  
11      continue  
        iy=iv(1)  
      endif  
      k=idum/IQ  
      idum=IA*(idum-k*IQ)-IR*k  
      if (idum.lt.0) idum=idum+IM  
      j=1+iy/NDIV  
      iy=iv(j)  
      iv(j)=idum  
      ran1=min(AM*iy,RNMX)  
      return  
      END  

      FUNCTION gammln(xx)  
      REAL*8 gammln,xx  
      INTEGER j  
      REAL*8 ser,stp,tmp,x,y,cof(6)  
      SAVE cof,stp  
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,  
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  
     *-.5395239384953d-5,2.5066282746310005d0/  
      x=xx  
      y=x  
      tmp=x+5.5d0  
      tmp=(x+0.5d0)*log(tmp)-tmp  
      ser=1.000000000190015d0  
      do 11 j=1,6  
        y=y+1.d0  
        ser=ser+cof(j)/y  
11    continue  
      gammln=tmp+log(stp*ser/x)  
      return  
      END  
