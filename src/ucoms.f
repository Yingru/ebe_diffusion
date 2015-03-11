cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c standard common block for uout2nt
c     Unit     : all modules
c     Author   : Steffen A. Bass, mofied by L.A.Winckelmann
c     Date     : 05/04/94
c     Revision : 1.0 beta - untested
c
c
      double precision EMNUC
      parameter (EMNUC = 0.938d0)

      integer nmax
c      parameter (nmax = 2000)
      parameter (nmax = 2000000)

c     this parameter MUST match MSTV(3) in BMS
      integer mstv3
      parameter (mstv3 = 1000)

      integer npart, ap,at,zp,zt
      integer event
      double precision  ebeam, bimp,ecm,pbeam,pcm
      double precision betann,betapro,betatar
c 7 integer

      character*8 model_tag, version_tag
      character*4 reffram

      common /names/model_tag,version_tag,reffram

      integer seed, refsys
      common /upara/ seed,refsys

      common /sys/ npart,event,ap,zp,at,zt

      common /rsys/bimp,ebeam,ecm,pbeam,pcm,
     +             betapro,betatar,betann
c 2*nmax*nmax logical

      integer ityp(nmax),origin(nmax)
c 6*nmax integer

      double precision eps, er0, pi, rho0, eps1
      parameter (eps  = 1.0E-12,
     +           eps1 = 1.0E-5,
     +           er0  = 1.128379167d0,
     +           pi   = 3.1415926535d0,
     +           rho0 = 0.16d0)


      double precision r0(nmax),rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax),px(nmax), py(nmax), pz(nmax),
     +     frr0(nmax),frrx(nmax), frry(nmax), frrz(nmax),
     +     fmass(nmax),lstcoll(nmax),pweight(nmax)
      

      common /isys/ ityp, lstcoll, origin
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass,
     &              frr0,frrx,frry,frrz,pweight





