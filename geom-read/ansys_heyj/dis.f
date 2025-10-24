      subroutine  dis(p1, p2, q1, q2, d)

      real p1(3), p2(3), q1(3), q2(3)

      d=0.0

      do i=1, 3
         d=d+(p1(i)-q1(i))*(p1(i)-q1(i))
         d=d+(p2(i)-q2(i))*(p2(i)-q2(i))
      enddo

      return
      end
