      subroutine assemsur_bd(ele,kip,freq,msurele,
     &   lp,s1hfemi,s1hfemj,schfemmatrix,hfemnum,
     &   hfemstore)

      implicit none

c.....Input Data
      integer kip,msurele,hfemstore
      integer lp(3,*)
      real    freq
      complex ele(3,3)

c.....Output Data
      integer hfemnum
      integer s1hfemi(*),s1hfemj(*)
      complex schfemmatrix(*)

c.....Working Variables
      real    pi
      integer i, j, k1, k2, indexo, indexe, ii, jj
      complex xj
      integer cofmatrix(3)
      pi=3.1415926
      xj=(0.,1.)
c        xj=(1.0,0.)
      do i=1,3
        do j=1,3
        ele(i,j)=20*pi*xj*freq*ele(i,j)/3.
        enddo
      enddo

      do i=1,3
          cofmatrix(i)=0
      end do

      do i=1,3
           cofmatrix(i)=(lp(i,kip))
      end do

      do i=1,3
         do j=1,3
            hfemnum=hfemnum+1
            s1hfemi(hfemnum)=cofmatrix(i)
            s1hfemj(hfemnum)=cofmatrix(j)
            schfemmatrix(hfemnum)=ele(i,j)
         end do
      end do
      RETURN
      END
