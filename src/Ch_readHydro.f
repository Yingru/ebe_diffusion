       subroutine hydro_read_CCNU
!       program hydro_read_CCNU
       implicit none
       include 'ucom_hydro.f'
     

       double precision :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a0
       integer :: t,i,j,k
       integer :: valid_NT, ios
       double precision :: dum_t

! determine the valid timestep of hydro    
       open(unit=20, file='Hydro_steps.dat')
       
       do j=1,NT 
           read(20, *, iostat=ios) dum_t
           if (ios/=0) exit
           if (j==NT) then
               write (6,*) "Error: maximum timestep for hydro record exceeds"
               write (6,*) "Exiting program now..."
               stop
           endif
           valid_NT=valid_NT+1
       enddo
       rewind(1)

       write(6,*) "CCNU Hydro ends at timesteps", valid_NT
! read initial hydro information       
       open(unit=10, file='CCNU_hydro.dat')

       do t=1, valid_NT
           read(10, *)
           read(10, *)
           do i=-NX, NX
               do j=-NY, NY
                   do k=-NZ, NZ
                       read(10, *)  a1, a2, a3, a4, a5, a6, a7, a0, a8, a9
                       posi_x(t, i, j, k)=a1
                       posi_y(t, i, j, k)=a2
                       Rap(t, i, j, k)=a3
                       Vx(t, i, j, k)=a4
                       Vy(t, i, j, k)=a5
                       Vz(t, i, j, k)=a6
                       Energy(t, i, j, k)=a7
                       CCNU_T(t, i, j, k)=a0
                       TJT(t, i, j, k)=a8
                       TJAT(t, i, j, k)=a9

                  enddo
               enddo
            enddo
       enddo
       close(10)

!      do i=-NZ, NZ
!         print *, posi_x(8,-NX,-NY,i), posi_y(8, -NX, NY, i),  Rap(8, -NX, NY, i),  Energy(8, -NX, NY, i)
!      enddo

      end

