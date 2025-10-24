      subroutine elesursrcm(kip,mvolnode,mvoledge,mvolele,xyz,nsurele,
     &           wpsp,wpse,le,et,ia,a)

      implicit none
!.....Input Data
      integer kip
      integer mvolnode,mvoledge,mvolele,nsurele
      real  xyz(3,mvolnode)
      integer wpsp(4,nsurele),wpse(3,nsurele)
      integer le(6,mvolele)
      real et(*)

!.....Output Data
      integer ia(3)
      complex a(3)
      
!.....Working Variables
      integer maxrule,maxgrid,irule
      integer ngrid(5)
      real  vt1(7,5),vt2(7,5),
     &        vt3(7,5),wt(7,5)

      integer ne1(3)
      real  length1(3)
      real  xx1(3),yy1(3),zz1(3)
      real  area1
      real  ns1(3)
      integer i,j,k,ii,i1,i2,i3,j1,j2,j3,k1,k2,k3
      real  xct,yct,zct
      real dotmul,temp1(3),temp2(3)
      real  temp,temp3(3),temp4(3)
      integer ktest
      real  pole(3)
      complex matrix(3)




      maxrule=5
      irule=5
      maxgrid=7
      call gausspoints(maxrule,maxgrid,ngrid,vt1,vt2,vt3,wt)
      
         
         

        do i=1,3
            matrix(i)=(0.0,0.0)
        enddo
        do i=1,3
            ne1(i)=wpsp(i,kip)
            xx1(i)=xyz(1,ne1(i))
            yy1(i)=xyz(2,ne1(i))
            zz1(i)=xyz(3,ne1(i))
        end do
        call getns(xx1,yy1,zz1,ns1)
        call areapatch(xx1,yy1,zz1,length1,area1)

        do j=1,ngrid(irule)

            xct=vt1(j,irule)*xx1(1)+vt2(j,irule)*xx1(2)
     &               +vt3(j,irule)*xx1(3)                 
            yct=vt1(j,irule)*yy1(1)+vt2(j,irule)*yy1(2)
     &               +vt3(j,irule)*yy1(3)
            zct=vt1(j,irule)*zz1(1)+vt2(j,irule)*zz1(2)
     &               +vt3(j,irule)*zz1(3)

            do k1=1,3
                i1=k1
                i2=mod(k1,3)+1
                i3=6-i1-i2
                call signcheck(ne1(i1),ne1(i2),2,1,ktest)
                temp=0.0

                temp1(1)=(xx1(i3)-xct)*length1(k1)/area1
                temp1(2)=(yy1(i3)-yct)*length1(k1)/area1
                temp1(3)=(zz1(i3)-zct)*length1(k1)/area1

c                print*,'temp1:',temp1
c                print*,'ns1:',ns1
        
                call xmul(temp1,ns1,temp2)

                call cal_et(xyz,xx1,yy1,zz1,
     &                      xct,yct,zct,
     &                      kip,nsurele,mvoledge,wpsp,
     &                      wpse,et,
     &                      pole(1),pole(2),pole(3),area1,
     &                      nsurele)
c                  print*,'pole:',pole
c                  print*,'temp2:',temp2

                temp=dotmul(temp2,pole)*real(ktest)


                matrix(k1)=matrix(k1)+temp*area1*wt(j,irule)*0.5
            end do
        enddo 

        do i=1,3
            ia(i)=wpse(i,kip)
            a(i)=matrix(i)
            write(9988,*)a(i),ia(i)
        enddo
                 

      
         
      
      return
      end

      subroutine getns(xx,yy,zz,norm)
        implicit none
        real xx(3),yy(3),zz(3)
        real norm(3)
        real vector21(3),vector31(3)
        real absn
        integer i
  
        vector21(1)=xx(2)-xx(1)
        vector21(2)=yy(2)-yy(1)
        vector21(3)=zz(2)-zz(1)
  
        vector31(1)=xx(3)-xx(1)
        vector31(2)=yy(3)-yy(1)
        vector31(3)=zz(3)-zz(1)
  
        call xmul(vector21,vector31,norm)
        absn=norm(1)**2+norm(2)**2+norm(3)**2
        absn=sqrt(absn)
        do i=1,3
           norm(i)=norm(i)/absn
        end do
        return
        end

      subroutine xmul(x,y,xy)
c         ! cross multiplication of two vectors.
      implicit none
      real x(3),y(3),xy(3)
         xy(1)=x(2)*y(3)-x(3)*y(2)
         xy(2)=x(3)*y(1)-x(1)*y(3)
         xy(3)=x(1)*y(2)-x(2)*y(1)
      return
      end

      subroutine dxmul(x,y,xy)
c         ! cross multiplication of two vectors.
      implicit none
      real*8 x(3),y(3),xy(3)
         xy(1)=x(2)*y(3)-x(3)*y(2)
         xy(2)=x(3)*y(1)-x(1)*y(3)
         xy(3)=x(1)*y(2)-x(2)*y(1)
      return
      end
