      subroutine parameterRead

      implicit none

      character*10 flag
      character*78 inputstr
      integer open_status
      include 'df_coms.f'
      include 'ucoms.f'

c set default values
      NUMSAMP=1
      flag_rad=2 ! 0-collsional; 1-radiative; 2-both.

      iflav=4     ! flavor of parton to select
      exp_setup=1 ! 1 for LHC 2.76~TeV and 2 for RHIC 200~GeV
      out_skip=0  ! if 1, only output the last time step
      reweight=2  ! 0: no re-weighting
                  ! 1: re-sample pT distribution according to wt table
                  ! 2: use uniform distribution
                  ! if 0, no evolution, output init. condition
      corr_flag=0 ! 1 for calculation of correlation
      HQ_input=1  ! 1:read in oscar file, 2:read in OSU file
                  ! 3:read in a list of x-y positions only
      ebe_flag=2  ! 0:smooth initial condition, 1:e-by-e study
                  ! 2:no rotation -- plain 
                  ! (use 2 when HQ_input is not 2)

c diffusion/noise parameters:
      D2piT=6d0
      alpha=6.2832d0/D2piT
      eta_cut=0.5d0   ! rapidity cut, symmetric inteval
      pT_init_min=0.25d0
      pT_init_max=70.25d0

c other paramters:
      static=2 ! 0-Chiho's hydro, 1-Static, 2-OSU hydro, 3-CCNU
      tsteps_cut=140
      static_cool=0
      T_static=0.30d0
      temp_cut=1 ! 0 for no cut, 1 for cut
      Tcut_critical=0.163d0
      wt_num=140 ! only meaningful for reweight=1
      wt_int=0.5d0 ! bin size of the weight table
      wt_Tab_min=0.5d0

      num_binary=0 ! 0: use number of xy points as npart
                   ! <0: use -num_binary as npart
                   ! >0: calculate npart with cross section
      sigma_c0=0.68432d0
      sigma_b0=0.056263d0

      OPEN(UNIT=15,FILE="parameters_df.dat",
     &     STATUS='OLD',FORM='FORMATTED',IOSTAT=open_status)

      if (open_status.ne.0) then
         write(6,*) "Error occurs during reading parameters_df.dat."
         write(6,*) "Use default parameters."
         goto 2
      endif

c read input lines
 1    continue
      read(15,99) flag,inputstr
 99   format(1A10,1A78)

c # : treat line as a comment
      if(flag(1:1).eq.'#') goto 1
      if(flag(1:4).eq."    ") goto 1
c xx: treat line as end of input marker
      if(flag(1:2).eq.'xx') goto 2

      SELECT CASE (flag)
         CASE ("NUMSAMP...")
            read(inputstr,fmt=*,err=88,end=88) NUMSAMP
         CASE ("iflav.....")
            read(inputstr,fmt=*,err=88,end=88) iflav
         CASE ("out_skip..")
            read(inputstr,fmt=*,err=88,end=88) out_skip
         CASE ("reweight..")
            read(inputstr,fmt=*,err=88,end=88) reweight
         CASE ("eta_cut...")
            read(inputstr,fmt=*,err=88,end=88) eta_cut
         CASE ("corr_flag.")
            read(inputstr,fmt=*,err=88,end=88) corr_flag
         CASE ("HQ_input..")
            read(inputstr,fmt=*,err=88,end=88) HQ_input
         CASE ("ebe_flag..")
            read(inputstr,fmt=*,err=88,end=88) ebe_flag
         CASE ("exp_setup.")
            read(inputstr,fmt=*,err=88,end=88) exp_setup
         CASE ("D2piT.....")
            read(inputstr,fmt=*,err=88,end=88) D2piT
            alpha=6.2832d0/D2piT
         CASE ("pT_min....")
            read(inputstr,fmt=*,err=88,end=88) pT_init_min
         CASE ("pT_max....")
            read(inputstr,fmt=*,err=88,end=88) pT_init_max
         CASE ("flag_rad..")
            read(inputstr,fmt=*,err=88,end=88) flag_rad
         CASE ("static....")
            read(inputstr,fmt=*,err=88,end=88) static
         CASE ("tsteps_cut")
            read(inputstr,fmt=*,err=88,end=88) tsteps_cut
         CASE ("stat_cool.")
            read(inputstr,fmt=*,err=88,end=88) static_cool
         CASE ("T_static..")
            read(inputstr,fmt=*,err=88,end=88) T_static
         CASE ("temp_cut..")
            read(inputstr,fmt=*,err=88,end=88) temp_cut
         CASE ("Tcut......")
            read(inputstr,fmt=*,err=88,end=88) Tcut_critical
         CASE ("wt_num....")
            read(inputstr,fmt=*,err=88,end=88) wt_num
         CASE ("wt_int....")
            read(inputstr,fmt=*,err=88,end=88) wt_int   
         CASE ("wt_Tab_min")
            read(inputstr,fmt=*,err=88,end=88) wt_Tab_min
         CASE ("sigma_c0..")
            read(inputstr,fmt=*,err=88,end=88) sigma_c0
         CASE ("sigma_b0..")
            read(inputstr,fmt=*,err=88,end=88) sigma_b0
         CASE ("num_binary")
            read(inputstr,fmt=*,err=88,end=88) num_binary
         CASE DEFAULT
            write(6,*) flag,"NOT a valid parameter!"
            write(6,*) "Terminating ..."
            stop
      END SELECT

      goto 1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 2    continue
      if(open_status.eq.0)
     &   write(6,*) "Parameters for Langevin has been read in."
      if(.TRUE.) then
         write(6,*) "flag_rad: ",flag_rad
         write(6,*) "iflav: ",iflav
         write(6,*) "D(2piT): ",D2piT
         write(6,*) "HQ_input: ",HQ_input
         write(6,*) "NUMSAMP: ",NUMSAMP
         write(6,*) "reweight: ",reweight
         write(6,*) "eta_cut: ",eta_cut
         write(6,*) "pT_min: ",pT_init_min
         write(6,*) "pT_max: ",pT_init_max
         write(6,*) "wt_num: ",wt_num
         write(6,*) "wt_int: ",wt_int
         write(6,*) "wt_Tab_min: ",wt_Tab_min
         write(6,*) "out_skip: ",out_skip
         write(6,*) "exp_setup: ",exp_setup
         write(6,*) "ebe_flag: ",ebe_flag
         write(6,*) "corr_flag: ",corr_flag
         write(6,*) "static: ",static
         write(6,*) "static_cool: ",static_cool
         write(6,*) "T_static: ",T_static
         write(6,*) "tsteps_cut: ",tsteps_cut 
         write(6,*) "temp_cut: ",temp_cut
         write(6,*) "Tcut_critical: ",Tcut_critical
         write(6,*) "sigma_c0: ",sigma_c0
         write(6,*) "sigma_b0: ",sigma_b0
         write(6,*) "num_binary: ",num_binary
      endif 
      return

c error-exit
 88   write(6,*) 'Syntax-error in the parameters_df.dat ...'
      write(6,*) 'Terminating ...'
      stop
      end

