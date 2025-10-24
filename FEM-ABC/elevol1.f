      subroutine ELEVOL1(ele,kip,freq,
     &                    mvolnode,mvolele, 
     &                    xyz,lv,
     &                    ep,qm,xi3)

      implicit none
        
c.....Input Data

      real    freq
      INTEGER kip, mvolnode, mvolele
      real    xyz(3,mvolnode)
      INTEGER lv(5,mvolele)
      COMPLEX ep(3,3), qm(3,3), xi3(3)

c.......Output Data
      COMPLEX ele(6,6)

c.....Working Variables

      real    length(4,4)
      INTEGER ne(4)
      COMPLEX elea(6,6),eleb(6,6), elec(6,6)
      real co(4),xx(4),yy(4),zz(4)
      real abcd(4,4)
      real rot(6,3)
      real    pi, k0, eta0, ck, v, det, xxx, yyy, zzz
      INTEGER i, ii, j, i1, i2, j1, j2, k1, k2, ktest
      real    co1, co2, co3, co4


      pi=3.1415926
      k0=20.*pi*freq/3.
      ck=k0*k0
      eta0=120*pi

      do 10 i=1,4
         ne(i)=lv(i,kip)
         co(i)=1
         xx(i)=xyz(1,ne(i))
         yy(i)=xyz(2,ne(i))
10       zz(i)=xyz(3,ne(i))
      v=0.

      do 20 i=1,4
         ii=mod(i,2)
         call detm(xx,yy,zz,i,det)
         if (ii.eq.1) then
            abcd(i,1)=det
         else
            abcd(i,1)=-det
         endif

         v=v+abcd(i,1)

         call detm(co,yy,zz,i,det)
         if (ii.eq.1) then
            abcd(i,2)=-det
         else
            abcd(i,2)=det
         endif

         call detm(co,xx,zz,i,det)
         if (ii.eq.1) then
            abcd(i,3)=det
         else
            abcd(i,3)=-det
         endif

         call detm(co,xx,yy,i,det)
         if (ii.eq.1) then
            abcd(i,4)=-det
         else
            abcd(i,4)= det
         endif

20    continue
c        print*,'v in elevol1',v
c        pause
      do 30 i=1,4
      do 30 j=1,4

         xxx=xx(i)-xx(j)
         yyy=yy(i)-yy(j)
         zzz=zz(i)-zz(j)

         length(i,j)=sqrt(xxx*xxx+yyy*yyy+zzz*zzz)
30    continue

      do i1=1,3
      do i2=i1+1,4
         call labeltran(i1,i2,i)
      do j=1,3
        j1=mod(j,3)+2
        j2=mod(j+1,3)+2
        rot(i,j)=abcd(i1,j1)*abcd(i2,j2)-abcd(i2,j1)*abcd(i1,j2)
      enddo
      enddo
      enddo

      do 40 i1=1,3
      do 40 i2=i1+1,4
      do 40 j1=1,3
      do 40 j2=j1+1,4

         call labeltran(i1,i2,k1)
         call labeltran(j1,j2,k2)
         elea(k1,k2)=0.0
         eleb(k1,k2)=0.0
         elec(k1,k2)=0.0

         do i=1,3
         do j=1,3

c.......discretization of rot E*rot E

        elea(k1,k2)=elea(k1,k2)+qm(i,j)*rot(k1,i)*rot(k2,j)

C.......discretization of E*E

        if (i1.eq.j1) then
           co1=2
        else
           co1=1
        endif
        if (i1.eq.j2) then
           co2=2
        else
           co2=1
        endif
        if (i2.eq.j1) then
           co3=2
        else
           co3=1
        endif
        if (i2.eq.j2) then
           co4=2
        else
           co4=1
        endif

        eleb(k1,k2)=eleb(k1,k2)+ep(i,j)*(co1*abcd(i2,i+1)*
     &                    abcd(j2,j+1)-co2*abcd(i2,i+1)*abcd(j1,j+1)-
     &                    co3*abcd(i1,i+1)*abcd(j2,j+1)+co4*abcd(i1,i+1)
     &                    *abcd(j1,j+1))              

        enddo

c......discretization of E * rot E

        elec(k1,k2)=elec(k1,k2)+(abcd(i2,i+1)-abcd(i1,i+1))*
     &              rot(k2,i)*xi3(i)

       enddo

        call signcheck(ne(i1),ne(i2),ne(j1),ne(j2),ktest)

       if (ktest.gt.0) then
         elea(k1,k2)= 4.0*elea(k1,k2)*length(i1,i2)*length(j1,j2)/v/v/v
         eleb(k1,k2)= eleb(k1,k2)*length(i1,i2)*length(j1,j2)/20.0/v
         elec(k1,k2)= elec(k1,k2)*length(i1,i2)*length(j1,j2)/v/v/4.0
       else
         elea(k1,k2)=-4.0*elea(k1,k2)*length(i1,i2)*length(j1,j2)/v/v/v
         eleb(k1,k2)=-eleb(k1,k2)*length(i1,i2)*length(j1,j2)/20.0/v
         elec(k1,k2)=-elec(k1,k2)*length(i1,i2)*length(j1,j2)/v/v/4.0
       endif

40     continue

      do 45 i=1,6
      do 45 j=1,6
45       ele(i,j)=(elea(i,j)-eleb(i,j)*ck)/6.0
!45       ele(i,j)=(elea(i,j)-2.*eta0*k0*elec(i,j)-eleb(i,j)*ck)/6.0
 
c        if (kip.eq.2) then
c           do i=1,6
c             do j=1,6
c               print*,ele(i,j)
c             end do
c           end do
c        else
c        end if
c        pause


      RETURN
      END
