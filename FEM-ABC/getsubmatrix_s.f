      subroutine getsubmatrix_s(mvolnode,mvolele, bfface,blpt,bsurele,
     & xyz,freq,
     & hfemstore,hfemrealnum,submatrixi,submatrixj,submatrix, 
     & judgee,facetotal,kein)
      implicit none

c..... Input Data
      integer mvolnode, mvolele, bsurele
      integer facetotal,kein
      real    xyz(3,mvolnode),freq
      integer bhfemstore,hfemstore
      integer bfface(4,facetotal),blpt(3,facetotal)
      integer judgee(*)
c.....Output Data
!      integer bmatrixi(bhfemstore),bmatrixj(bhfemstore)
!      complex bmatrix(bhfemstore)

      integer hfemrealnum,bhfemnum

      integer submatrixi(hfemstore),submatrixj(hfemstore)
      complex submatrix(hfemstore)

c......Working Variable
      integer i, j, k,iele
      integer,allocatable::bmatrixi(:), bmatrixj(:)
      complex,allocatable::bmatrix(:)
      real,allocatable::rdis(:)
      integer,allocatable::ipof(:),ipfo(:),ip(:)
      complex elep(3,3)
      integer,allocatable:: s1truehfemi(:),s1truehfemj(:)
      complex,allocatable:: s1truehfemmatrix(:)
      integer num1i,num1j,num2i,num2j,num3i,num3j
      real    small,dxx,dyy
      integer nt,hfemnum

      bhfemnum=0

      bhfemstore=3*3*(bsurele)
      allocate(bmatrixi(bhfemstore))
      allocate(bmatrixj(bhfemstore))
      allocate(bmatrix(bhfemstore))


      do i=1,bhfemstore
         bmatrixi(i)=0
         bmatrixj(i)=0
         bmatrix(i)=0.0
      end do

cccccccccccccccccccccccccccccccc get bmatrix cccccccccccccccccccccccccccccccccccccccccccc

         hfemnum=0
         if(bsurele.gt.0)then

         do i=1,bsurele
           iele=i

           call elesur_md(elep,iele,mvolele,
     &              mvolnode,xyz,bfface)

           call assemsur_bd(elep,iele,freq,bsurele,
     &  blpt,bmatrixi,bmatrixj,bmatrix,hfemnum,
     &  bhfemstore)
         enddo
ccccccccccccccccccccccc compress unknownsccccccccccccccccccccccccccccccccccccccccccccccc
        allocate (rdis(hfemnum),ipof(hfemnum),ipfo(hfemnum),ip(hfemnum))

        do i=1,hfemnum
           ip(i)=0
           rdis(i)=0
           ipof(i)=0
           ipfo(i)=0
        enddo

        small = 0.001

        do i=1,hfemnum
           rdis(i)=abs(bmatrixi(i)+bmatrixj(i))
        end do
        call hpsort(hfemnum,rdis,ipof,ipfo)

        nt = 1
        ip(nt) = ipfo(1)
        ipof(ipfo(1)) = nt

        do i= 2,hfemnum
            j=i-1
13      if (abs(rdis(j)-rdis(i)).gt.small) goto 23
          num1i = bmatrixi(ipfo(i))
          num1j = bmatrixj(ipfo(i))

          num2i = bmatrixi(ipfo(j))
          num2j = bmatrixj(ipfo(j))

           dxx = num1i-num2i
           dyy = num1j-num2j

        if ((abs(dxx).le.small).and.(abs(dyy).le.small)) then

             ipof(ipfo(i))=ipof(ipfo(j))
              goto 43
        else
            j=j-1
           if (j.le.0) goto 23
            goto 13
        endif
23        nt=nt+1
           ipof(ipfo(i))=nt
           ip(nt)=ipfo(i)
43      enddo

        bhfemnum = nt
        allocate (s1truehfemi(nt),s1truehfemj(nt),s1truehfemmatrix(nt))

        do i=1,bhfemnum
           s1truehfemi(i)=0
           s1truehfemj(i)=0
           s1truehfemmatrix(i)=0.
        end do

        do i=1,bhfemnum
           s1truehfemi(i) = bmatrixi(ip(i))
           s1truehfemj(i) = bmatrixj(ip(i))
        end do

        do i=1,hfemnum
            s1truehfemmatrix(ipof(i))=s1truehfemmatrix(ipof(i))
     &              +bmatrix(i)
        end do

        do i=1,3*3*bsurele
           bmatrixi(i)= 0
           bmatrixj(i)= 0
           bmatrix(i)= 0.0
        enddo

        do i=1,bhfemnum
           bmatrixi(i)= s1truehfemi(i)
           bmatrixj(i)= s1truehfemj(i)
           bmatrix(i)= s1truehfemmatrix(i)
        enddo

        deallocate (s1truehfemi,s1truehfemj,s1truehfemmatrix)
        deallocate (rdis,ipof,ipfo,ip)
        end if
cccccccccccccccccccccccccccccccc add mmatrix to vol cccccccccccccccccccccccccccccccccccc

      hfemnum=bhfemnum+hfemrealnum

      allocate(s1truehfemi(hfemnum),
     &       s1truehfemj(hfemnum),
     &       s1truehfemmatrix(hfemnum))
      k=0
      do i=1,bhfemnum
         k=k+1
         s1truehfemi(k) = bmatrixi(i)
         s1truehfemj(k) = bmatrixj(i)
       s1truehfemmatrix(k)=bmatrix(i)
      enddo

      do i=1,hfemrealnum
         k=k+1
         s1truehfemi(k) = submatrixi(i)
         s1truehfemj(k) = submatrixj(i)
       s1truehfemmatrix(k)=submatrix(i)
      enddo

      allocate (rdis(hfemnum),ipof(hfemnum),ipfo(hfemnum),
     &          ip(hfemnum))

      do i=1,hfemnum
        ip(i)=0
        rdis(i)=0
        ipof(i)=0
        ipfo(i)=0
      enddo

      small = 0.001

      do i=1,hfemnum
        rdis(i)=abs(s1truehfemi(i)+s1truehfemj(i))
      end do

      call hpsort(hfemnum,rdis,ipof,ipfo)

      nt = 1
      ip(nt) = ipfo(1)
      ipof(ipfo(1)) = nt

      do i= 2,hfemnum
         j=i-1
12     if (abs(rdis(j)-rdis(i)).gt.small) goto 22
          num1i = s1truehfemi(ipfo(i))
          num1j = s1truehfemj(ipfo(i))

          num2i = s1truehfemi(ipfo(j))
          num2j = s1truehfemj(ipfo(j))

           dxx = num1i-num2i
           dyy = num1j-num2j

        if ((abs(dxx).le.small).and.(abs(dyy).le.small)) then

             ipof(ipfo(i))=ipof(ipfo(j))
              goto 42
        else
            j=j-1
           if (j.le.0) goto 22
            goto 12
        endif
22        nt=nt+1
           ipof(ipfo(i))=nt
           ip(nt)=ipfo(i)
42      enddo

      hfemrealnum = nt

      do i=1,hfemstore
         submatrixi(i)=0
         submatrixj(i)=0
         submatrix(i)=0.
      end do

      do i=1,hfemrealnum
         submatrixi(i) = s1truehfemi(ip(i))
         submatrixj(i) = s1truehfemj(ip(i))
      end do

      do i=1,hfemnum
       submatrix(ipof(i))=submatrix(ipof(i))
     &              +s1truehfemmatrix(i)
      end do
cccccccccccccccccccccccccccc enfbnd cccccccccccccccccccccccccccccccccccc
      if(kein.gt.0)then
      do i=1,hfemrealnum
        if ((judgee(submatrixi(i)).eq.1)
     &  .or.(judgee(submatrixj(i))).eq.1)then

          if (submatrixi(i).eq.submatrixj(i)) then
            submatrix(i) = 1.0
          else
            submatrix(i) = 0.0
          end if
        end if
      end do
      end if
        hfemnum=0
      do i=1,hfemrealnum
       if(submatrix(i).ne.0.)then
         hfemnum=hfemnum+1
         s1truehfemi(hfemnum) = submatrixi(i)
         s1truehfemj(hfemnum) = submatrixj(i)
         s1truehfemmatrix(hfemnum)= submatrix(i)
       endif
      enddo
      hfemrealnum=hfemnum

      do i=1,hfemstore
         submatrixi(i)=0
         submatrixj(i)=0
         submatrix(i)=0.
      end do

      do i=1,hfemrealnum
       submatrixi(i)=s1truehfemi(i)
       submatrixj(i)=s1truehfemj(i)
       submatrix(i)=s1truehfemmatrix(i)
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       deallocate (s1truehfemi,s1truehfemj,s1truehfemmatrix)
       deallocate (rdis,ipof,ipfo,ip)
 
      RETURN
      END
