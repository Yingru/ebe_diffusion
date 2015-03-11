ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   This file contains the functions and subroutines for the calculation
c   of heavy quark energy loss due to single gluon radiation
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function alphas(kT0,temp0)

c This is a function to calculate the strong coupling costant, temp is the 
c medium temperature, and alphas=alphas(pi*temp) if kT<pi*temp due to 
c screening effect. It's a first order calculation.

      implicit none

      include 'ucoms.f'
      double precision alphas,kT0,kT1,temp0,nflavor,lambdas,error_para

      error_para=1.d0

      if(kT0.lt.PI*temp0*error_para) then
         kT1=PI*temp0*error_para
      else
         kT1=kT0
      endif

      alphas=4.d0*PI/(11.d0-2.d0*nflavor(kT1)/3.d0)/2.d0/
     &       Log(kT1/lambdas(kT1))

c      write(*,*) alphas

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function nflavor(kT0)

      implicit none

      include 'df_coms.f'
      double precision kT0,nflavor

      if (kT0.lt.cMass) then
         nflavor=3.d0
      elseif (kT0.lt.bMass) then
         nflavor=4.d0
      else
         nflavor=5.d0
      endif

c      write(*,*) nflavor

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function lambdas(kT0)

      implicit none

      include 'df_coms.f'
      double precision kT0,lambdas

      if (kT0.lt.cMass) then
         lambdas=0.2d0
      elseif (kT0.lt.bMass) then
         lambdas=0.172508d0
      else
         lambdas=0.130719d0
      endif

c      write(*,*) lambdas

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function splittingP(z0)

      implicit none

      include 'df_coms.f'
      double precision splittingP,z0

      splittingP = (2.d0-2.d0*z0+z0*z0)*C_F/z0

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function tau_f(x0g,y0g)

      implicit none
      
      include 'df_coms.f'
      double precision tau_f,x0g,y0g

      tau_f = 2.d0*HQenergy*x0g*(1.d0-x0g)/((x0g*y0g*HQenergy)**2
     &        +x0g*x0g*HQmass*HQmass)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function dNg_over_dxdydt(x0g,y0g)

c single gluon radiation formula from Wang's paper

      implicit none

      include 'df_coms.f'
      include 'ucoms.f'

      double precision dNg_over_dxdydt,alphas,splittingP,tau_f
      double precision x0g,y0g

      qhat = 4.d0*alpha*temp_med**3/C_F

      if(x0g*HQenergy.lt.PI*temp_med) then
         dNg_over_dxdydt = 0.d0
         goto 111
      endif

      if(tau_f(x0g,y0g).lt.1.d0/PI/temp_med) then
         dNg_over_dxdydt = 0.d0
         goto 111
      endif   

c      if(x0g*y0g*HQenergy.lt.sqrt(2d0/3d0)*PI*temp_med) then
c         dNg_over_dxdydt = 0.d0
c         goto 111
c      endif   

      dNg_over_dxdydt = 4.d0/PI*N_c*alphas(x0g*y0g*HQenergy,temp_med)*
     &                 splittingP(x0g)*qhat*y0g**5*(sin(time_gluon
     &                 /2.d0/tau_f(x0g,y0g)/inv_fm_to_GeV))**2*
     &                 (HQenergy*HQenergy/(y0g*y0g*HQenergy*HQenergy+
     &                 HQmass*HQmass))**4/x0g/x0g/HQenergy/HQenergy/
     &                 inv_fm_to_GeV

c      write(6,*) dNg_over_dxdydt

 111   continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
