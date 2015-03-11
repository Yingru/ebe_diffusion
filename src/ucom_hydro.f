cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  standard common block 
c 


      integer  NT, NX, NY, NZ
      parameter (NT = 101)
      parameter (NX = 81)
      parameter (NY = 81)
      parameter (NZ = 51)
!      parameter (max_Ncol=20000)

      double precision tau_min
      parameter (tau_min = 0.6)

      double precision posi_x(1:NT, -NX:NX, -NY:NY, -NZ:NZ), 
     &                 posi_y(1:NT, -NX:NX, -NY:NY, -NZ:NZ),
     &                 Rap(1:NT, -NX:NX, -NY:NY, -NZ:NZ), 
     &                 Vx(1:NT, -NX:NX, -NY:NY, -NZ:NZ),
     &                 Vy(1:NT, -NX:NX, -NY:NY, -NZ:NZ),
     &                 Vz(1:NT, -NX:NX, -NY:NY, -NZ:NZ),
     &                 Energy(1:NT, -NX:NX, -NY:NY, -NZ:NZ),  
     &                 TJT(1:NT, -NX:NX, -NY:NY, -NZ:NZ), 
     &                 TJAT(1:NT, -NX:NX, -NY:NY, -NZ:NZ),
     &                 CCNU_T(1:NT, -NX:NX, -NY:NY, -NZ:NZ)
      
! make sure the initial matrix are zero (number=NT*(2NX+1)*(2NY+1)*(2NZ+1))
! it turns out fortran will initially set value about 0
!       data posi_x /98517852*0/
!       data posi_y /98517852*0/
!       data Rap /98517852*0/
!       data Vx /98517852*0/
!       data Vy /98517852*0/
!       data Vz /98517852*0/
!       data Energy /98517852*0/
!       data TJT /98517852*0/
!       data TJAT /98517852*0/
!       data CCNU_T /98517852*0/

!      double precision, dimension(1:NT, -NX:NX, -NY:NY, -NZ:NZ) :: posi_x, posi_y, Rap, Vx, Vy, Vz, Energy, TJT, TJAT
      common /hydroInfo/posi_x, posi_y, Rap, Vx, Vy, Vz, 
     &                  Energy, TJT, TJAT, CCNU_T
      

