ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     common-blocks:

      integer AMAXX,AMAXY,AMAXZ,ntsteps,ndiss,ncoeff,ncharge
      INTEGER MAXX,MAXY,MAXZ,MAXN,static,static_cool,
     &        rotation,temp_cut
      
      double precision initt,T_static,alpha,w,Tcut_critical,D2piT
      double precision tmax,x0,xmax,y0,ymax,h0,hmax,tau
      integer pdep_flag,iflav,ref_frame,exp_setup,tsteps_cut
      double precision eta_cut
      integer ebe_flag

      integer mxpart,npt,evsamp,des_cntr,out_skip,HQ_input
      parameter (mxpart=2000000)
!      parameter (mxpart=20000)
      PARAMETER (AMAXX =169,AMAXY=169,AMAXZ=169)

c set number of Langevin steps within hydro step 
      integer nlang
      integer Moore,reweight,wt_num,corr_flag,wt_MAXnum,num_binary
      double precision wt_int,wt_Tab_min
      parameter(nlang=10) ! eqiv to dt_lang=0.01 for MABIKT=10
      parameter(evsamp=1) 
      parameter(pdep_flag=0)
      parameter(rotation=0) ! several choices, for checking purpose
      parameter(ref_frame=1) ! ref_frame: 1-c.m frame, 2-local rest frame
      parameter(Moore=0) ! re-sample pT distribution by power law
      parameter(w=1.0d0) ! weight between white noise and memory
      parameter(wt_MAXnum=300)

      double precision cpT_wt(wt_MAXnum),bpT_wt(wt_MAXnum)
      double precision pT_init_min,pT_init_max

      common/weight/cpT_wt,bpT_wt
      common/pTinit/pT_init_min,pT_init_max,eta_cut
      common/userFlag/reweight,iflav,out_skip,corr_flag,HQ_input,
     &                ebe_flag,exp_setup
      common/diffusionConst/D2piT,alpha
      common/otherFlag1/static,static_cool,T_static,
     &                  temp_cut,wt_num,tsteps_cut,num_binary
      common/otherFlag2/Tcut_critical,wt_int,wt_Tab_min

      integer tstep,lstep,flag_stop
      double precision h_rx(0:amaxx,0:amaxy,0:amaxz),
     &                 h_ry(0:amaxx,0:amaxy,0:amaxz)
      double precision h_reta(0:amaxx,0:amaxy,0:amaxz)
      double precision h_vx(0:amaxx,0:amaxy,0:amaxz)
      double precision h_vy(0:amaxx,0:amaxy,0:amaxz)
      double precision h_veta(0:amaxx,0:amaxy,0:amaxz)
      double precision h_rqgp(0:amaxx,0:amaxy,0:amaxz)
      double precision temp(0:amaxx,0:amaxy,0:amaxz)
      double precision ener(0:amaxx,0:amaxy,0:amaxz)
      double precision press(0:amaxx,0:amaxy,0:amaxz)
      double precision nb(0:amaxx,0:amaxy,0:amaxz),
     &                 mu(0:amaxx,0:amaxy,0:amaxz)


      double precision par_plane_phi
      double precision sigma_pptot,sigma_ctot,sigma_btot
      double precision sigma_c0,sigma_b0

      common/ihydro/tstep,lstep,maxx,maxy,maxz,maxn,flag_stop
      common/ihydro2/ntsteps,ndiss,ncoeff,ncharge
      common/rhydro/h_rx,h_ry,h_reta
      common/rhydro2/initt,tmax,x0,xmax,y0,ymax,h0,hmax,tau
      common/vhydro/h_vx,h_vy,h_veta
      common/thydro/temp,mu,h_rqgp,ener,press,nb
      common/cross/sigma_ctot,sigma_btot,sigma_pptot,
     &             sigma_c0,sigma_b0
      common/ebeAngle/par_plane_phi


      integer pid(mxpart,evsamp)
      double precision p_rx(mxpart,evsamp),p_ry(mxpart,evsamp),
     &                 p_rz(mxpart,evsamp)
      double precision p_px(mxpart,evsamp),p_py(mxpart,evsamp),
     &       p_pz(mxpart,evsamp),p_p0(mxpart,evsamp)
      double precision p_mass(mxpart,evsamp),p_r0(mxpart,evsamp),
     &       p_reta(mxpart,evsamp),Thydro(mxpart,evsamp),
     &       c_vx(mxpart,evsamp),c_vy(mxpart,evsamp),
     &       c_vz(mxpart,evsamp),p_wt(mxpart,evsamp),
     &       edensity(mxpart,evsamp),sdensity(mxpart,evsamp)
      double precision xi_x(mxpart,evsamp),xi_y(mxpart,evsamp),
     &                 xi_z(mxpart,evsamp)
      double precision p_ipx(mxpart,evsamp),p_ipy(mxpart,evsamp),
     &       p_ipz(mxpart,evsamp),p_ip0(mxpart,evsamp),
     &       p_ipT(mxpart,evsamp)
      double precision p_vbx(mxpart,evsamp),p_vby(mxpart,evsamp),
     &       p_vbz(mxpart,evsamp),P_kvtx(mxpart,evsamp),
     &       p_kvty(mxpart,evsamp),p_kvtz(mxpart,evsamp),dEsum


      common/ipart/pid,npt,des_cntr
      common/rpart/p_rx,p_ry,p_rz,p_px,p_py,p_pz,p_p0,p_mass,
     &     p_r0,p_reta,xi_x,xi_y,xi_z,
     &     p_ipx,p_ipy,p_ipz,p_ip0,p_ipT,p_vbx,p_vby,p_vbz,
     &     p_kvtx,p_kvty,p_kvtz,dEsum,Thydro,c_vx,c_vy,c_vz,p_wt,
     &     edensity,sdensity

      double precision rrlu
      integer mrlu
      COMMON/VNIDAR/MRLU(6),RRLU(100)

      double precision inv_fm_to_GeV
      parameter(inv_fm_to_GeV=0.1973d0)

      integer numXY,NUMSAMP
      double precision initX(mxpart),initY(mxpart) 
      common/initXY1/numXY,NUMSAMP
      common/initXY2/initX,initY      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c parameters/variables for gluon radiation

      double precision C_F,N_c
      parameter (N_c=3.d0)
      parameter (C_F=4.d0/3.d0)

      integer flag_rad,narrow
      parameter(narrow=4)
      common/radflag1/flag_rad

      integer fgluon(mxpart,evsamp),init_int(mxpart,evsamp)
      double precision time_lim(mxpart,evsamp)
      double precision previous_kT(mxpart,evsamp)
      common/radflag2/fgluon,time_lim,init_int,previous_kT

      integer HQener_gn,t_gn,temp_gn
      parameter (HQener_gn=200)
      parameter (t_gn=70)
      parameter (temp_gn=100)

      double precision dNg_over_dt(t_gn+1,temp_gn,HQener_gn)
      double precision max_dNgfnc(t_gn+1,temp_gn,HQener_gn)

      double precision HQener_max,t_max,temp_max,temp_min
      parameter(HQener_max=80.d0)
      parameter(t_max=14.d0)
      parameter(temp_max=0.55d0)
      parameter(temp_min=0.15d0)

      double precision delta_HQener,delta_tg,delta_temp

      common/int_table/dNg_over_dt,delta_HQener,delta_tg,delta_temp
     &                 ,max_dNgfnc

      double precision cMass,bMass
      parameter(cMass=1.27d0,bMass=4.19d0)

      double precision HQenergy,HQmass,time,time_init,temp_med,time_tg
      double precision qhat,time_gluon
      double precision t_init(mxpart,evsamp),Tinteval_lrf(mxpart,evsamp)

      common/rad1/HQenergy,HQmass,time,time_init,temp_med,time_tg
      common/rad2/qhat,t_init,Tinteval_lrf,time_gluon

      integer time_num,HQenergy_num,temp_num
      integer cprob_gt1,ctemp_gtmax,ct_gtmax,num_gluon,cHQE_gtmax
      integer count_Ecut
      double precision delta_Ng,Etot_gluon

      common/rad_search1/time_num,HQenergy_num,temp_num
      common/rad_search2/delta_Ng
      common/rad_debug/cprob_gt1,ctemp_gtmax,ct_gtmax,num_gluon,
     &                 cHQE_gtmax,count_Ecut,Etot_gluon
