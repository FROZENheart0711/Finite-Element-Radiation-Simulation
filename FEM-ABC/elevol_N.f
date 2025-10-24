        subroutine elevol_N(xx,yy,zz,xe,ye,ze,
     &                      vol,length,N,b,c,d)

        IMPLICIT NONE

c.......Input Data

        real     xx(4),yy(4),zz(4),xe,ye,ze 
        real     vol,length(4,4),b(4),c(4),d(4)

c.......Output Data

        real    N(6,3)

c.......Working Variables

        INTEGER  i,ii,i1,i2
        real   x(4),y(4),z(4),a(4),v(4)
        real   L(4),gradL(4,3),det


        do i=1,4
          v(i)=0.
        enddo

ccccccccc..............v(1)..................cccccccccccc

        do i=1,4
          x(i)=xx(i)
          y(i)=yy(i)
          z(i)=zz(i)
        enddo
        
        x(1)=xe
        y(1)=ye
        z(1)=ze

        do 10 i=1,4
        ii=mod(i,2)
        call detm(x,y,z,i,det)       
        if (ii.eq.1) then
          a(i)=det
        else
          a(i)=-det
        endif
10      v(1)=v(1)+a(i)  
       
cccccccc................v(2).................cccccccccccc

        do i=1,4
          x(i)=xx(i)
          y(i)=yy(i)
          z(i)=zz(i)
        enddo

        x(2)=xe
        y(2)=ye
        z(2)=ze

        do 20 i=1,4
        ii=mod(i,2)
        call detm(x,y,z,i,det)       
        if (ii.eq.1) then
          a(i)=det
        else
          a(i)=-det
        endif
20      v(2)=v(2)+a(i)

cccccccc................v(3)................cccccccccccc

        do i=1,4
          x(i)=xx(i)
          y(i)=yy(i)
          z(i)=zz(i)
        enddo

        x(3)=xe
        y(3)=ye
        z(3)=ze

        do 30 i=1,4
        ii=mod(i,2)
        call detm(x,y,z,i,det)       
        if (ii.eq.1) then
          a(i)=det
        else
          a(i)=-det
        endif
30      v(3)=v(3)+a(i)

cccccccc...............v(4).................cccccccccccc

        do i=1,4
          x(i)=xx(i)
          y(i)=yy(i)
          z(i)=zz(i)
        enddo

        x(4)=xe
        y(4)=ye
        z(4)=ze

        do 40 i=1,4
        ii=mod(i,2)
        call detm(x,y,z,i,det)       
        if (ii.eq.1) then
          a(i)=det
        else
          a(i)=-det
        endif
40      v(4)=v(4)+a(i)

ccccccccccccc...............L(i).................cccccccccccc

        do i=1,4
          L(i)=v(i)/vol
        enddo

ccccccccccccc...............gradL(i).............cccccccccccc

        do i=1,4
          gradL(i,1)=b(i)/vol
          gradL(i,2)=c(i)/vol
          gradL(i,3)=d(i)/vol
        enddo

ccccccccccccc...............N(i).................cccccccccccc


        do i1=1,3
        do i2=i1+1,4

          call labeltran(i1,i2,i)

          do ii=1,3
           N(i,ii)=(L(i1)*gradL(i2,ii)-L(i2)*gradL(i1,ii))*length(i1,i2)
          enddo
      
        enddo
        enddo


        RETURN
        END
