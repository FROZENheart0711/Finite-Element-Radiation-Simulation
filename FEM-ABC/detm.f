      subroutine detm(x,y,z,i,det)

      implicit none

c.....Input Data
      integer i
      real    x(4),y(4),z(4)

c.....Output Data
      real    det

c.....Working Variables
      integer i1, i2, i3
      i1=mod(i,4)+1
      i2=mod(i+1,4)+1
      i3=mod(i+2,4)+1

      det=x(i1)*y(i2)*z(i3)+x(i2)*y(i3)*z(i1)+x(i3)*y(i1)*z(i2)
     &      -x(i1)*y(i3)*z(i2)-x(i2)*y(i1)*z(i3)-x(i3)*y(i2)*z(i1)
     
      RETURN
      END
