cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_osc_event(iret)

      implicit none
      include 'ucoms.f'

      character*12 oscar_tag, file_tag
      character*1 cdummy1
      character*3 cdummy3
      character*170 inputline
      integer ntestp,iret,aap,aat,i,j

      double precision embeam,emtarget,eb,galab,belab,gaeq,beeq,ppeq
      double precision pteq,pcms,dummy


      read(unit=5,fmt=306,err=199,end=199)inputline

 306  format(a170)


      if(inputline(1:3).eq.'OSC') then
c     file header:

         read (unit=inputline,fmt=901,err=198,end=199) oscar_tag
         
         read (unit=5,fmt=901,err=198,end=199) file_tag

 901     format (a12)

         read (unit=5,fmt=306,err=198,end=199) inputline

         read (unit=inputline,fmt=902,err=198,end=199) 
     &        model_tag, version_tag, cdummy1, Ap, cdummy1, 
     &        Zp, cdummy3, At, cdummy1, Zt, cdummy1,
     &        reffram, ecm, ntestp

 902     format (2(a8,2x),a1,i3,a1,i6,a3,i3,a1,i6,a1,2x,a4,2x,
     &        e10.4,2x,i8)

         iret=1

         aap=abs(Ap)
         aat=abs(At)

c assume nuclear projectile and target for the time being
         embeam=EMNUC
         emtarget=EMNUC
c         
         pcm=pcms(ecm,embeam,emtarget)
         ebeam=sqrt(embeam**2 + (pcm*ecm/emtarget)**2) - embeam
         pbeam=pcm*ecm/emtarget
c     now revert to full quantities
         ebeam=AAp*ebeam
         pbeam=AAp*pbeam
         embeam=embeam*AAp
         emtarget=emtarget*AAt
         eb=sqrt(embeam**2+pbeam**2)


c now do the calculation of equal_speed quantities

         galab=eb/embeam        ! gamma_lab
         belab=pbeam/eb         ! beta_lab
         gaeq=sqrt(0.5*(1+galab)) ! gamma_equal_speed
         beeq=belab*galab/(1+galab) ! beta_equal_speed
         ppeq=gaeq*beeq*embeam  ! p_projectile(eq)
         pteq=-(gaeq*beeq*emtarget) ! p_target(eq)

c reduce to per particle quantities
         ppeq=ppeq/dble(AAp)
         if(AAt.ne.0) then
            pteq=pteq/dble(AAt)
            emtarget=emtarget/dble(AAt)
         endif
         embeam=embeam/dble(AAp)
         pbeam=pbeam/dble(AAp)
         ebeam=ebeam/dble(AAp)
c     the following is the NN sqrt(s)
         ecm=sqrt(embeam**2+2*eb/dble(AAp)*emtarget+emtarget**2) 
ccccccccccccccccccccccccccccccccccccccccccccc

c     compute transformation betas for output

         pcm=max(1d-10,pbeam*emtarget/ecm)

c     write(6,*) ' srt, ebeam, pbeam,pcm',ecm,ebeam,pbeam,pcm


         betann=0.d0
         betatar=pcm/sqrt(pcm**2+emtarget**2)
         betapro=-(1.*pcm/sqrt(pcm**2+embeam**2))

cccccccccccccccccccccccccccccccccccccccc

c read next line for event header
         read(unit=5,fmt=306,err=199,end=199)inputline

      endif

c event header

         read (unit=inputline,fmt=903,err=198,end=199) 
     &        event, npart, bimp, dummy

 903     format (i10,2x,i10,2x,f8.3,2x,f8.3)

c  read particles

         do 99 i=1,npart
            read(unit=5,fmt=904,err=198,end=199) j, ityp(i), 
     .           px(i), py(i), pz(i), p0(i), fmass(i),     
     .           rx(i), ry(i), rz(i), r0(i), lstcoll(i),
     .           frrx(i),frry(i),frrz(i),frr0(i),origin(i)
 99      continue

c     904  format (i10,2x,i10,2x,9(e12.6,2x))
 904     FORMAT(I10,2X,I10,14(2X,D12.6),I5)


         return


 198  continue
      iret=-1
      write(6,*) 'ERROR encountered while reading OSCAR file'
      write(6,*) '   ... terminating  '
      return



 199  continue
      iret=0
      write(6,*) 'EOF encountered while reading OSCAR file'
      write(6,*) '   ... terminating  '
      return

      end


      double precision function pcms(srt,m1,m2)
c calculates the CM-momentum in a 2-body decay/coll. depending on ecm
      implicit none
      double precision srt,m1,m2,s

      if (srt.lt.m1+m2) then
         pcms = 0.0
         return
      endif

      s = srt*srt
      pcms = sqrt( (s-(m1+m2)**2)*(s-(m1-m2)**2)/(4.0*s))
      return
      end
