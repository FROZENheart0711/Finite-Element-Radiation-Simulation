      subroutine dis2o(xyz1,xyz2,r)

      real xyz1(3), xyz2(3),r

      r=0
      do i=1,3
         xyzc=(xyz1(i)+xyz2(i))/2.
         r=r+xyzc*xyzc
      enddo

      return
      end
