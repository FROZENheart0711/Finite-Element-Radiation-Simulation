        subroutine devaNsur(u,v,N)

        real*8 u,v,w
        complex N(6)
        integer i
      w=1.-u-v
    
      N(1)=(2*u-1)*u;

      N(2)=(2*v-1)*v;

      N(3)=(2*w-1)*w;

      N(4)=4*u*v;

      N(5)=4*v*w;

      N(6)=4*w*u

      i=1

      if(i.eq.0)then

      N(1)=u;

      N(2)=v;

      N(3)=w;

      endif

        return
        end


