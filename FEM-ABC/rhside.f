      subroutine rhside(mvolnode,mvoledge,mvolele,xyz,nsurele,
     &           nfface,nlpt,le,et,crhs,kz, ep, mu, mxmat, lv,
     &           amp,phase)

      implicit none
!.....Input Data
      
      integer mvolnode,mvoledge,mvolele,nsurele
      integer mxmat
      real  xyz(3,mvolnode)
      integer nfface(4,*),nlpt(3,*)
      integer le(6,mvolele)
      real et(*)
      complex kz
      complex ep(3,3,mxmat), mu(3,3,mxmat)
      integer lv(5,mvolele)
      real amp,phase
      

!.....Output Data
      complex crhs(mvoledge)
      
!.....Working Variables
      integer kip
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
      complex kcmplx
      complex eptt,qmtt
      complex xj
      real pi
      real E0

      read(15,*) amp,phase
      close(15)

      xj=(0.0,1.0)
      pi=3.14159265359
      E0=sqrt(120*pi)

      eptt=ep(1,1,lv(5,nfface(4,1)))
      qmtt=mu(1,1,lv(5,nfface(4,1))) 

      kcmplx=sqrt(eptt*qmtt)
      print*,'kcmplx:',kcmplx

      maxrule=5
      irule=5
      maxgrid=7
      call gausspoints(maxrule,maxgrid,ngrid,vt1,vt2,vt3,wt)
      
         
         

      do kip=1,nsurele
        do i=1,3
            ne1(i)=nfface(i,kip)
            xx1(i)=xyz(1,ne1(i))
            yy1(i)=xyz(2,ne1(i))
            zz1(i)=xyz(3,ne1(i))
        end do
        call getns(xx1,yy1,zz1,ns1)
        call areapatch(xx1,yy1,zz1,length1,area1)
        do i=1,3
            matrix(i)=(0.0,0.0)
        enddo

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
        
                call xmul(temp1,ns1,temp2)

                call cal_et(xyz,xx1,yy1,zz1,
     &                      xct,yct,zct,
     &                      kip,nsurele,mvoledge,nfface,
     &                      nlpt,et,
     &                      pole(1),pole(2),pole(3),area1,
     &                      nsurele)

                temp=dotmul(temp2,pole)*real(ktest)


                matrix(k1)=matrix(k1)+temp*area1*wt(j,irule)*0.5

            end do
        enddo 

        do i=1,3
            crhs(nlpt(i,kip))=
     &       crhs(nlpt(i,kip))+
     &          matrix(i)*2.0*xj*kz*amp*exp(xj*phase/180*pi)*
     &          E0/sqrt(kcmplx)/qmtt

        enddo
       enddo
                 

      
         
      
      return
      end
