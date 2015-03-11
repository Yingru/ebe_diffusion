
      subroutine read_radTable
      
      implicit none
      include 'df_coms.f'

      integer i,j,k,rad_read
      character*10 dummy_c
      integer dummy_n
c      double precision time_tg
      character*77 file10
      double precision time_tab,temp_tab,HQenergy_tab

      file10='       '
      call getenv('ftn10',file10)
      if (file10(1:4).ne.'    ') then
         OPEN(UNIT=10,FILE=file10,STATUS='OLD',FORM='FORMATTED')
      else 
         write(6,*) "Error in loading dNg_over_dt table :("
         stop
      endif
         
      k=1
      do while(k.lt.t_gn+1)
         i=1

         read(unit=10,fmt=2111,err=2098,end=2097) 
     &        dummy_c,dummy_n,dummy_c,time_tg

c read in the integral table for one time step

         do while(i.lt.temp_gn+1)
            j=1
            do while(j.lt.HQener_gn+1)
               read (unit=10,fmt=2112,err=2098,end=2097) 
     &              time_tab,temp_tab,HQenergy_tab,
     &              dNg_over_dt(k+1,i,j),max_dNgfnc(k+1,i,j)
               j=j+1
            enddo
            i=i+1
         enddo
         k=k+1
      enddo


 2111 format(A10,2X,I3,2X,A6,2X,F8.3)
 2112 format(e10.4,4(2X,e10.4))

c      write(6,*) dNg_over_dt(71,100,200)
c      write(6,*) max_dNgfnc(71,100,200)
      write(6,*) "Read integral table successfully :)"

      close(10)
      return

 2097   continue
        rad_read=1
        write(6,*) 'ERROR: EOF reached in dNg_over_dt table'
        write(6,*) 'terminating ...'
        return

 2098   continue
        rad_read=-1
        write(6,*) 'READ-ERROR in dNg_over_dt table'
        write(6,*) 'terminating ...'
        return

      end
         
