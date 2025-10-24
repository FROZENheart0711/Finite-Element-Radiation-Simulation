      function cdotmul(x,y)
      implicit none
      complex cdotmul,y(3)
      real x(3)
      cdotmul=cmplx(x(1),0)*y(1)+cmplx(x(2),0)*y(2)+cmplx(x(3),0)*y(3)
      return
      end
