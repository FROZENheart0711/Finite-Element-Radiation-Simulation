        subroutine matset( mvoledge, mvolstore, iqvol, ipvol, volmat,
     &  ia, ja, a)
        implicit none
        integer mvoledge,mvolstore
        integer iqvol(*), ipvol(*)
        complex volmat(*)
        integer ia(*), ja(*)
        complex a(*)
        integer i,j
        integer ii
        ii=0
        do i=1,mvoledge
           do j=iqvol(i),iqvol(i+1)-1
              ii=ii+1
              ia(ii)=i
              ja(ii)=ipvol(j)
              a(ii)=volmat(j)
           enddo
         enddo

         return
         end
