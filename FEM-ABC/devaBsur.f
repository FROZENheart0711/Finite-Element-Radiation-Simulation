        subroutine devaBsur(u,v,B)

        real*8 u,v
        complex B(2,8)

      B(1,1) = -v
      B(2,1) = u

      B(1,3) = - ( v - 1 )
      B(2,3) = u
  
      B(1,2) = - v
      B(2,2) = u - 1

      B(1,4) = v
      B(2,4) = u

      B(1,6) = 1 - 2*u - v
      B(2,6) = - u

      B(1,5) = -v
      B(2,5) = 1 - u - 2*v
      
      B(1,7) = v + u*v - v**2
      B(2,7) = u + u*v - u**2

      B(1,8) = v - u*v - v**2
      B(2,8) = -(u - u**2 - u*v)
c        do i=4,8
c          do j=1,2
c            B(j,i)=0
c          enddo
c        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        return
        end


