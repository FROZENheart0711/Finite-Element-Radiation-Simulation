        subroutine devagradNsur(u,v,gradN)

        real*8 u,v,w
        complex gradN(2,6)
        integer i

      w=1.-u-v
    
c      N(1)=(2*u-1)*u;
      gradN(1,1)=4*u-1.
      gradN(2,1)=0.

c      N(2)=(2*v-1)*v;
      gradN(1,2)=0
      gradN(2,2)=4*v-1.

c      N(3)=(2*w-1)*w;
      gradN(1,3)=-4.+4*u+4*v+1.
      gradN(2,3)=-4.+4*u+4*v+1.

c      N(4)=4*u*v;
      gradN(1,4)=4*v
      gradN(2,4)=4*u
      
c      N(5)=4*v*w;
      gradN(1,5)=-4*v
      gradN(2,5)=4.-4*u-8*v

c      N(6)=4*w*u
      gradN(1,6)=4.-8*u-4*v
      gradN(2,6)=-4*u

      i=1

      if(i.eq.0)then

      gradN(1,1)=1.
      gradN(2,1)=0.

c      N(2)=(2*v-1)*v;
      gradN(1,2)=0
      gradN(2,2)=1.

c      N(3)=(2*w-1)*w;
      gradN(1,3)=-1.
      gradN(2,3)=-1.

      endif


        return
        end


