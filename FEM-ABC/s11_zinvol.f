      subroutine s11_zinvol(mvolnode,mvolele,mvoledge,
     &           xyz,lv,le,nfface,nlpt,nsurele,
     &           crhs,et1,mxmat,ep,mu,amp,phase)

      implicit none
!.....Input Data
      integer mvolnode,mvolele,mvoledge
      integer nsurele
      real  xyz(3,mvolnode)
      integer lv(5,mvolele),le(6,mvolele)
      integer nfface(4,nsurele),nlpt(3,nsurele)
      complex crhs(mvoledge)
      real et1(*)
      real amp,phase
      integer mxmat
      complex ep(3,3,mxmat), mu(3,3,mxmat)

!.....Output Data
      real s11

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
      integer isub,kip
      integer i,j,k,ii,i1,i2,i3,j1,j2,j3,k1,k2,k3
      real  xct,yct,zct
      real  r_r,rrt
      real  ddotmul,temp1(3),temp2(3)
      real  temp,temp3(3),temp4(3)
      integer ktest
      real  thesin,thecos,phisin,phicos
      real  rki(3),pole(3),n_E(3),n_rki(3)
      real EMod(3)
      real  pi,E0,eTEM
      complex zin,reflt

      real  length(4,4)
      integer ne(4)
      real  co(4),xx(4),yy(4),zz(4)
      real  v, det, xxx, yyy, zzz
      real  a(4),b(4),c(4),d(4)
      integer iele
     
      real  delNv(3),n_delNv(3),rn_delNv(3),rrn_delNv(3)
      real  Nv(6,3)
      real  n_Nv(3)
      real  r_in,r_out,x_feed,y_feed,z_feed,kz
      integer iport,iiport
      complex enRlt,enInc,Et(3),Einc(3),Er(3)
      complex xj
      complex eptt,qmtt
      complex kcmplx

c.......Begin
      open(111,file='portE.txt',status='unknown')

      pi=3.14159265359
      E0=sqrt(120.0*pi) !1.0
      xj=(0.0,1.0)
      eptt=ep(1,1,lv(5,nfface(4,1)))
      qmtt=mu(1,1,lv(5,nfface(4,1))) 
      print*,'ep:',eptt
      print*,'mu:',qmtt
      print*,'amp:',amp
      print*,'phase:',phase

      kcmplx=sqrt(eptt*qmtt)


      maxrule=5
      irule=5
      maxgrid=7
      call gausspoints(maxrule,maxgrid,ngrid,vt1,vt2,vt3,wt)

c      do i=1,ngrid(irule)
c         print*,'vt wt:',vt1(i,irule),vt2(i,irule),vt3(i,irule),
c     &            wt(i,irule)
c      enddo



            reflt=0.0
            enRlt=0.0
            enInc=0.0

           do 2000 ii=1,nsurele
              kip=ii
              do i=1,3
                 ne1(i)=nfface(i,kip)
                 xx1(i)=xyz(1,ne1(i))
                 yy1(i)=xyz(2,ne1(i))
                 zz1(i)=xyz(3,ne1(i))
              end do
              call getns(xx1,yy1,zz1,ns1)
              call areapatch(xx1,yy1,zz1,length1,area1)
              iele=nfface(4,kip)
              do i=1,4
                 ne(i)=lv(i,iele)
                 co(i)=1
                 xx(i)=xyz(1,ne(i))
                 yy(i)=xyz(2,ne(i))
                 zz(i)=xyz(3,ne(i))

              end do
              v=0.
              do 20 i=1,4
                 j=mod(i,2)
                 call detm(xx,yy,zz,i,det)
                 if (j.eq.1) then
                   a(i)=det
                 else
                   a(i)=-det
                 endif
            
                 v=v+a(i)

                 call detm(co,yy,zz,i,det)
                 if (j.eq.1) then
                    b(i)=-det
                 else 
                 b(i)=det
                 endif

                 call detm(co,xx,zz,i,det)
                 if (j.eq.1) then
                    c(i)=det
                 else
                    c(i)=-det
                 endif

                 call detm(co,xx,yy,i,det)
                 if (j.eq.1) then
                    d(i)=-det
                 else
                    d(i)= det
                 endif
20            continue

              do 30 i=1,4
              do 30 j=1,4
                 xxx=xx(i)-xx(j)
                 yyy=yy(i)-yy(j)
                 zzz=zz(i)-zz(j)
                 length(i,j)=sqrt(xxx*xxx+yyy*yyy+zzz*zzz)
30            continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              do j=1,ngrid(irule)
                 xct=vt1(j,irule)*xx1(1)+vt2(j,irule)*xx1(2)
     &              +vt3(j,irule)*xx1(3)
                 yct=vt1(j,irule)*yy1(1)+vt2(j,irule)*yy1(2)
     &              +vt3(j,irule)*yy1(3)
                 zct=vt1(j,irule)*zz1(1)+vt2(j,irule)*zz1(2)
     &              +vt3(j,irule)*zz1(3)
                 
               call cal_et(xyz,xx1,yy1,zz1,
     &                      xct,yct,zct,
     &                      kip,nsurele,mvoledge,nfface,
     &                      nlpt,et1,
     &                      pole(1),pole(2),pole(3),area1,
     &                      nsurele)
               
c               write(111,*) xct,yct,zct
c               write(111,*) abs(pole(1)),abs(pole(2)),abs(pole(3))
                 call elevol_N(xx,yy,zz,xct,yct,zct,v,length,Nv,b,c,d)
 
                 Et(1)=0.0
                 Et(2)=0.0
                 Et(3)=0.0
                 
                 do 45 i1=1,3
                 do 45 i2=i1+1,4
                    call labeltran(i1,i2,i)
                    do k=1,3
                       temp1(k)=Nv(i,k)
                    end do
c                    call cro_mul(ns1,temp2,temp1)
                    call signcheck(lv(i1,iele),lv(i2,iele),2,1,ktest)
                    if(ktest.gt.0)then
                    Et(1)=Et(1)+temp1(1)*crhs(le(i,iele))
                    Et(2)=Et(2)+temp1(2)*crhs(le(i,iele))
                    Et(3)=Et(3)+temp1(3)*crhs(le(i,iele))
                    else
                    Et(1)=Et(1)-temp1(1)*crhs(le(i,iele))
                    Et(2)=Et(2)-temp1(2)*crhs(le(i,iele))
                    Et(3)=Et(3)-temp1(3)*crhs(le(i,iele))
                    endif
45               continue

               write(111,*) xct,yct,zct
               write(111,*) abs(Et(1)),abs(Et(2)),abs(Et(3))

                  Einc(1)=E0*amp*exp(xj*phase/180*pi)*
     &                  pole(1)/sqrt(kcmplx)

                  Einc(2)=E0*amp*exp(xj*phase/180*pi)*
     &                  pole(2)/sqrt(kcmplx)

                  Einc(3)=E0*amp*exp(xj*phase/180*pi)*
     &                  pole(3)/sqrt(kcmplx)

                  EMod(1)=pole(1)
c*exp(-xj*realkz(imode)*z_feed)
                  EMod(2)=pole(2)
c*exp(-xj*realkz(imode)*z_feed)
                  EMod(3)=pole(3)

                  Er(1)=Et(1)-Einc(1)
                  Er(2)=Et(2)-Einc(2)
                  Er(3)=Et(3)-Einc(3)

                  enRlt=enRlt+(Er(1)*EMod(1)
     &                 +Er(2)*EMod(2)+
     &                 Er(3)*EMod(3))*area1*wt(j,irule)*0.5

                  enInc=enInc+(Einc(1)*EMod(1)+Einc(2)*EMod(2)+
     &                   Einc(3)*EMod(3))*area1*wt(j,irule)*0.5

              enddo
2000       end do
            
           reflt=enRlt/enInc

           s11=20*log10(abs(reflt))

           print*,'s11:',s11

           close(111)

      return
      end

