      subroutine fdposi
     & ( mvoladj,
     &   ket, klt, ke, kl,
     &   lv,  le,  ln, lp, ip, iq,
     &   rp,  iqof,iqfo,   ipele,
     &   mvolstore, msurstore)

      implicit none

      integer i, j, n
      integer mvoladj, mvolstore,msurstore

      integer lv(9, *),le(6, *)
      integer ln(4, *),lp(3, *)
      integer ip(*), iq(*)
      integer ipele(mvoladj, *)

      integer iqof(*), iqfo(*)
	    real    rp(*)

      integer ii,jj,kk,k,mtl,num,ket,klt,ke,kl

      do i=1,ket
         ip(i)=0
      enddo

      do i=1,klt
      do j=1,6
         n=le(j,i) 
         ip(n)=ip(n)+1
         if(ip(n).gt.mvoladj)then
           print*,'error1!',ip(n),mvoladj,n,i,j,ket
           exit
         endif
         ipele(ip(n),n)=i
      end do
      end do

c.......find the number of neighbour edges for every edge in volume

      do 210 i=1,ket
         mtl=0
      do 220 j=1,ip(i)
         jj=ipele(j,i)
      do 220 k=1,6
         kk=le(k,jj)
         if (i.gt.kk) goto 220
            mtl=mtl+1
            rp(mtl)=kk
 220     continue
      call hpsort(mtl,rp,iqof,iqfo)
      num=1
      do 230 j=1,mtl-1
         jj=j+1
         if (rp(jj).ne.rp(j)) then
           num=num+1
         else
         endif
 230  continue
      ip(i)=num
 210  continue

      mvolstore=1
      do i=1,ket
         mvolstore=mvolstore+ip(i)
      enddo
c.......find the connecting elements for every edge in surface

      do i=1,ke
         iq(i)=0
      enddo
      do i=1,kl
      do j=1,3
         n=lp(j,i)
         if(n.gt.0)then
            iq(n)=iq(n)+1
            if(iq(n).gt.mvoladj)then
                print*,'error 3!',iq(n),mvoladj,n
                exit
            endif
            ipele(iq(n),n)=i
         endif
      enddo
      enddo

c.......find the number of neighbour edges for every edge in surface

      do 110 i=1,ke
         mtl=0
      do 120 j=1,iq(i)
         jj=ipele(j,i)
      do 120 k=1,3
         kk=lp(k,jj)
         if (i.gt.kk) goto 120
            mtl=mtl+1
            rp(mtl)=lp(k,jj)
 120     continue
      call hpsort(mtl,rp,iqof,iqfo)
      num=1
      do 130 j=1,mtl-1
         jj=j+1
         if (rp(jj).ne.rp(j)) then
         num=num+1
      else
      endif
 130  continue
      iq(i)=num
 110  continue

      msurstore=1
      do i=1,ke
         msurstore=msurstore+iq(i)
      enddo
      return
      end
