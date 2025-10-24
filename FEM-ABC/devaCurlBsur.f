        subroutine devaCurlBsur(u,v,CurlB)

        real*8 u,v
        complex CurlB(3,8)

      CurlB(1,1) = 0.
      CurlB(2,1) = 0.
      CurlB(3,1) = 2.
   
      CurlB(1,2) = 0.
      CurlB(2,2) = 0.
      CurlB(3,2) = 2.

      CurlB(1,3) = 0.
      CurlB(2,3) = 0.
      CurlB(3,3) = 2.

      CurlB(1,4) = 0.
      CurlB(2,4) = 0.
      CurlB(3,4) = 0.

      CurlB(1,5) = 0.
      CurlB(2,5) = 0.
      CurlB(3,5) = 0.

      CurlB(1,6) = 0.
      CurlB(2,6) = 0.
      CurlB(3,6) = 0.

      CurlB(1,7) = 0.
      CurlB(2,7) = 0.
      CurlB(3,7) = 3*v - 3*u

      CurlB(1,8) = 0.
      CurlB(2,8) = 0.
      CurlB(3,8) = 3*u + 3*v - 2.



      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        return
        end


