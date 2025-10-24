      subroutine elesur_md(elep,kip,mvolele,mvolnode,xyz,lc)

      implicit none

c.....Input Data
      integer kip, mvolnode, mvolele
      real    xyz(3,mvolnode)
      integer lc(4,4*mvolele)

c.....Output Data
      complex elep(3,3)

c.....Working Variables
      real    length(3),xn(3),yn(3)
      integer ne(3)
      real    xx(3),yy(3),zz(3)
      real    a(3),b(3),c(3)
      integer i, k1, k2, j1, j2, l1, l2, ktest
      real    area


      do 10 i=1,3
10       ne(i)=lc(i,kip)

      do 20 i=1,3
         xx(i)=xyz(1,ne(i))
         yy(i)=xyz(2,ne(i))
20       zz(i)=xyz(3,ne(i))

      length(1)=sqrt((xx(1)-xx(2))**2+(yy(1)-yy(2))**2
     &            +(zz(1)-zz(2))**2)
      length(2)=sqrt((xx(2)-xx(3))**2+(yy(2)-yy(3))**2
     &            +(zz(2)-zz(3))**2)
      length(3)=sqrt((xx(1)-xx(3))**2+(yy(1)-yy(3))**2
     &            +(zz(1)-zz(3))**2)

      xn(1)=0.
      yn(1)=0.
      xn(2)=length(1)
      yn(2)=0.
      xn(3)=(length(1)**2+length(3)**2-length(2)**2)/
     &        length(1)/2.
      yn(3)=sqrt(length(3)**2-xn(3)**2)

      do 90 i=1,3
         k1=mod(i,3)+1
         k2=mod(i+1,3)+1
         a(i)=xn(k1)*yn(k2)-xn(k2)*yn(k1)
         b(i)=yn(k1)-yn(k2)
90       c(i)=xn(k2)-xn(k1)

         area=b(1)*c(2)-b(2)*c(1)
c        print*,'area',area
      do 100 j1=1,3
      do 100 j2=1,3

         if (j1.eq.3) then
            k1=3
            k2=1
         else
            k1=j1
            k2=j1+1
         endif

         if (j2.eq.3) then
            l1=3
            l2=1
         else
            l1=j2
            l2=j2+1
         endif

         elep(j1,j2)=0.

         if (k1.eq.l1) then
            elep(j1,j2)=elep(j1,j2)+(b(k2)*b(l2)+c(l2)*c(k2))/12.
         else
            elep(j1,j2)=elep(j1,j2)+(b(k2)*b(l2)+c(l2)*c(k2))/24.
         endif

         if (k1.eq.l2) then
            elep(j1,j2)=elep(j1,j2)-(b(k2)*b(l1)+c(l1)*c(k2))/12.
         else
            elep(j1,j2)=elep(j1,j2)-(b(k2)*b(l1)+c(l1)*c(k2))/24.
         endif

         if (k2.eq.l1) then
            elep(j1,j2)=elep(j1,j2)-(b(k1)*b(l2)+c(l2)*c(k1))/12.
         else
            elep(j1,j2)=elep(j1,j2)-(b(k1)*b(l2)+c(l2)*c(k1))/24.
         endif

         if (k2.eq.l2) then
            elep(j1,j2)=elep(j1,j2)+(b(k1)*b(l1)+c(l1)*c(k1))/12.
         else
            elep(j1,j2)=elep(j1,j2)+(b(k1)*b(l1)+c(l1)*c(k1))/24.
         endif

        call signcheck(lc(k1,kip),lc(k2,kip),lc(l1,kip),
     &                 lc(l2,kip),ktest)

        if (ktest.gt.0) then
           elep(j1,j2)=elep(j1,j2)*length(j1)*length(j2)/area
        else
           elep(j1,j2)=-elep(j1,j2)*length(j1)*length(j2)/area
        endif
100   continue

      RETURN
      END
