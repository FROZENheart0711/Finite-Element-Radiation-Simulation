        subroutine cal_et(xyz,xx,yy,zz,xct,yct,zct,
     &   kip,nsrc,mvoledge,wpsp,wpse,et,
     &   exx,eyy,ezz,area,
     &   totalsrc) ! newadded
     

      implicit none
!.....Input Data
      real xx(3),yy(3),zz(3),xct,yct,zct
      integer kip
      integer mvoledge
      integer nsrc
      integer wpsp(4,*)
      integer wpse(3,*)
      real xyz(3,*)
      integer totalsrc
      real et(*)
!.....Output Data
      real exx,eyy,ezz
      real area
!.....Working Variables
      integer i,j,k1,k2,k3,ktest1,ktest2,ktest3
      integer ii,jj
      real*8 xx1(3),yy1(3),zz1(3)
      real*8 length(3),length1(3)
      real*8 area1,areal(3)
      real*8 xn(3),yn(3)
      real*8 w_xn(3),w_yn(3)
      real*8 ex,ey
      real*8 xx_xyz(3),temp31(3),tempzz(3),yy_xyz(3)
      real*8 xxlength,yylength
      double precision u,v,weigh,detJ,detJ1
      real*8 x1,x2,x3,y1,y2,y3,z1,z2,z3
      complex*8 Jmatrix(2,3),IJmatrix(3,2),TIJmatrix(2,3)
      complex*8 Jmatrix1(3,3),IJmatrix1(3,3),TJmatrix1(3,3)
      complex*8 B1(2,8)
      complex*8 Eu,Ev
c.....begin
      call areapatch(xx,yy,zz,length,area)

        do i=1,3
          
          xx1=xx
          yy1=yy
          zz1=zz
          xx1(i)=xct
          yy1(i)=yct
          zz1(i)=zct
c          print*,"xx1",i,xx1(i),yy1(i),zz1(i)

          call dareapatch(xx1,yy1,zz1,length1,area1)

          areal(i)=area1/area

        enddo

        u=areal(1)
        v=areal(2)

      call devaBsur(u,v,B1)

      x1=xx(1)
      x2=xx(2)
      x3=xx(3)
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      z1=zz(1)
      z2=zz(2)
      z3=zz(3)
        
      detJ1=sqrt( ((y1-y3)*(z2-z3)-(z1-z3)*(y2-y3))**2+
     &         ((z1-z3)*(x2-x3)-(x1-x3)*(z2-z3))**2+
     &         ((x1-x3)*(y2-y3)-(y1-y3)*(x2-x3))**2)

        Jmatrix1(1,1)=x1-x3
        Jmatrix1(1,2)=y1-y3
        Jmatrix1(1,3)=z1-z3
        Jmatrix1(2,1)=x2-x3
        Jmatrix1(2,2)=y2-y3
        Jmatrix1(2,3)=z2-z3
        Jmatrix1(3,1)=1/detJ1*((y1-y3)*(z2-z3)-
     &           (z1-z3)*(y2-y3) )
        Jmatrix1(3,2)=1/detJ1*((z1-z3)*(x2-x3)-
     &           (x1-x3)*(z2-z3) )
        Jmatrix1(3,3)=1/detJ1*((x1-x3)*(y2-y3)-
     &           (y1-y3)*(x2-x3) )
         do ii=1,2
            do jj=1,3
               Jmatrix(ii,jj)=Jmatrix1(ii,jj)
            enddo
         enddo

         call invJmatrixf(Jmatrix1,IJmatrix1,detJ)

         detJ=detJ1

         do ii=1,3
             do jj=1,2
                IJmatrix(ii,jj)=IJmatrix1(ii,jj)
             enddo
         enddo

        do ii=1,2
         do jj=1,3
            TIJmatrix(ii,jj)=IJmatrix(jj,ii)
         enddo
        enddo

        do ii=1,3
         do jj=1,3
            TJmatrix1(ii,jj)=Jmatrix1(jj,ii)
         enddo
        enddo

        call signcheck(wpsp(1,kip),wpsp(2,kip),
     &                 1,2,ktest1)
        call signcheck(wpsp(2,kip),wpsp(3,kip),
     &                 1,2,ktest2)
        call signcheck(wpsp(3,kip),wpsp(1,kip),
     &                 1,2,ktest3)

         k1=wpse(1,kip)
         k2=wpse(2,kip)
         k3=wpse(3,kip)

         Eu=et(k1)*B1(1,1)*ktest1+et(k2)*B1(1,2)*ktest2
     &         +et(k3)*B1(1,3)*ktest3
         Ev=et(k1)*B1(2,1)*ktest1+et(k2)*B1(2,2)*ktest2
     &         +et(k3)*B1(2,3)*ktest3

         Eu=Eu+et(k1+mvoledge)*B1(1,4)+et(k2+mvoledge)*B1(1,5)
     &         +et(k3+mvoledge)*B1(1,6)
         Ev=Ev+et(k1+mvoledge)*B1(2,4)+et(k2+mvoledge)*B1(2,5)
     &         +et(k3+mvoledge)*B1(2,6)

         Eu=Eu+et(kip+mvoledge*2)*B1(1,7)+
     &         et(kip+mvoledge*2+totalsrc)*B1(1,8)
         Ev=Ev+et(kip+mvoledge*2)*B1(2,7)+
     &         et(kip+mvoledge*2+totalsrc)*B1(2,8)
c         print*,'edge:',k1,k2,k3
c         print*,'et:',et(k1),et(k2),et(k3)
         exx=IJmatrix(1,1)*Eu+IJmatrix(1,2)*Ev
         eyy=IJmatrix(2,1)*Eu+IJmatrix(2,2)*Ev
         ezz=IJmatrix(3,1)*Eu+IJmatrix(3,2)*Ev
c         print*,'e:',exx,eyy,ezz
         

        
        

        return
        end

        subroutine cal_et2(xyz,xx,yy,zz,xct,yct,zct,
     &   kip,nsrc,mvoledge,wpsp,wpse,et,
     &   exx,eyy,ezz,area,p_srcspfrt,
     &   p_mvolnodefrt,totalsrc) ! newadded
     

      implicit none
!.....Input Data
      real*8 xx(3),yy(3),zz(3),xct,yct,zct
      integer kip
      integer mvoledge
      integer nsrc
      integer wpsp(4,*)
      integer wpse(3,*)
      real xyz(3,*)
      integer totalsrc
      real et(*)
!.....Output Data
      real*8 exx,eyy,ezz
      real*8 area
!.....Working Variables
      integer i,j,k1,k2,k3,ktest1,ktest2,ktest3
      integer ii,jj
      real*8 xx1(3),yy1(3),zz1(3)
      real*8 length(3),length1(3)
      real*8 area1,areal(3)
      real*8 xn(3),yn(3)
      real*8 w_xn(3),w_yn(3)
      real*8 ex,ey
      real*8 xx_xyz(3),temp31(3),tempzz(3),yy_xyz(3)
      real*8 xxlength,yylength
      integer*8 p_srcspfrt
      integer*8 p_mvolnodefrt
      double precision u,v,weigh,detJ,detJ1
      real*8 x1,x2,x3,y1,y2,y3,z1,z2,z3
      complex*8 Jmatrix(2,3),IJmatrix(3,2),TIJmatrix(2,3)
      complex*8 Jmatrix1(3,3),IJmatrix1(3,3),TJmatrix1(3,3)
      complex*8 B1(2,8)
      complex*8 Eu,Ev
c.....begin

      call dareapatch(xx,yy,zz,length,area)

        do i=1,3

          xx1=xx
          yy1=yy
          zz1=zz
          xx1(i)=xct
          yy1(i)=yct
          zz1(i)=zct

          call dareapatch(xx1,yy1,zz1,length1,area1)

          areal(i)=area1/area

        enddo

        u=areal(1)
        v=areal(2)


      call devaBsur(u,v,B1)

      x1=xx(1)
      x2=xx(2)
      x3=xx(3)
      y1=yy(1)
      y2=yy(2)
      y3=yy(3)
      z1=zz(1)
      z2=zz(2)
      z3=zz(3)
        
      detJ1=sqrt( ((y1-y3)*(z2-z3)-(z1-z3)*(y2-y3))**2+
     &         ((z1-z3)*(x2-x3)-(x1-x3)*(z2-z3))**2+
     &         ((x1-x3)*(y2-y3)-(y1-y3)*(x2-x3))**2)

        Jmatrix1(1,1)=x1-x3
        Jmatrix1(1,2)=y1-y3
        Jmatrix1(1,3)=z1-z3
        Jmatrix1(2,1)=x2-x3
        Jmatrix1(2,2)=y2-y3
        Jmatrix1(2,3)=z2-z3
        Jmatrix1(3,1)=1/detJ1*((y1-y3)*(z2-z3)-
     &           (z1-z3)*(y2-y3) )
        Jmatrix1(3,2)=1/detJ1*((z1-z3)*(x2-x3)-
     &           (x1-x3)*(z2-z3) )
        Jmatrix1(3,3)=1/detJ1*((x1-x3)*(y2-y3)-
     &           (y1-y3)*(x2-x3) )
         do ii=1,2
            do jj=1,3
               Jmatrix(ii,jj)=Jmatrix1(ii,jj)
            enddo
         enddo

         call invJmatrixf(Jmatrix1,IJmatrix1,detJ)

         detJ=detJ1

         do ii=1,3
             do jj=1,2
                IJmatrix(ii,jj)=IJmatrix1(ii,jj)
             enddo
         enddo

        do ii=1,2
         do jj=1,3
            TIJmatrix(ii,jj)=IJmatrix(jj,ii)
         enddo
        enddo

        do ii=1,3
         do jj=1,3
            TJmatrix1(ii,jj)=Jmatrix1(jj,ii)
         enddo
        enddo

        call signcheck(wpsp(1,kip),wpsp(2,kip),
     &                 1,2,ktest1)
        call signcheck(wpsp(2,kip),wpsp(3,kip),
     &                 1,2,ktest2)
        call signcheck(wpsp(3,kip),wpsp(1,kip),
     &                 1,2,ktest3)

         k1=wpse(1,kip)
         k2=wpse(2,kip)
         k3=wpse(3,kip)

         Eu=et(k1)*B1(1,1)*ktest1+et(k2)*B1(1,2)*ktest2
     &         +et(k3)*B1(1,3)*ktest3
         Ev=et(k1)*B1(2,1)*ktest1+et(k2)*B1(2,2)*ktest2
     &         +et(k3)*B1(2,3)*ktest3

         Eu=Eu+et(k1+mvoledge)*B1(1,4)+et(k2+mvoledge)*B1(1,5)
     &         +et(k3+mvoledge)*B1(1,6)
         Ev=Ev+et(k1+mvoledge)*B1(2,4)+et(k2+mvoledge)*B1(2,5)
     &         +et(k3+mvoledge)*B1(2,6)

         Eu=Eu+et(kip+mvoledge*2)*B1(1,7)+
     &         et(kip+mvoledge*2+totalsrc)*B1(1,8)
         Ev=Ev+et(kip+mvoledge*2)*B1(2,7)+
     &         et(kip+mvoledge*2+totalsrc)*B1(2,8)
c         print*,'edge:',k1,k2,k3
c         print*,'et:',et(k1),et(k2),et(k3)
         exx=IJmatrix(1,1)*Eu+IJmatrix(1,2)*Ev
         eyy=IJmatrix(2,1)*Eu+IJmatrix(2,2)*Ev
         ezz=IJmatrix(3,1)*Eu+IJmatrix(3,2)*Ev
c         print*,'e:',exx,eyy,ezz
         

        
        

        return
        end
