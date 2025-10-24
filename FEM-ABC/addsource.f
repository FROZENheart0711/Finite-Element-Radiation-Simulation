      subroutine add_source(thetai,phii,ipol,freq,mvolnode,mvolele,
     &           mvoledge,xyz,lv,le,bfface,blpt,bsurele,
     &           crhs)


      IMPLICIT NONE

C.....Input Data
      real    thetai,phii,freq
      integer ipol,mvolnode,mvolele,mvoledge
      real    xyz(3,mvolnode)
      integer lv(5,mvolele),le(6,mvolele)
      integer bfface(4,4*mvolele),blpt(3,4*mvolele)
      integer bsurele
c.....Output Data
      complex crhs(mvoledge)
 
c.....Working Data
      real    k0,r_r
      complex ci,xj
      real    pi
      parameter(pi=3.141592653)
c      parameter (rad      = 0.017453293)
      integer ktest
      integer i,j,k,     n1,     n2,     n3
      real    rki(3), pole(3)
      real    thesin, thecos, phisin, phicos
      integer isub,kip,iface,iele
      real    length(4,4)
      integer ne(4)
      real    co(4),xx(4),yy(4),zz(4)
      real    v, det, xxx, yyy, zzz,length_sur(3),area
      integer ii, i1, i2, j1, j2, k1, k2
      real    a(4),b(4),c(4),d(4),e23(3),e21(3),ns(3)
      real    nsabs
      real    as(3),bs(3),cs(3)
      real    xn(3),yn(3),points(3,3),pc(3,3)
      real    Nv(6,3)
      real    temp1(3),temp2(3)
      real    n_E(3),n_rki(3),dotmul
      complex temp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      thesin=sin(thetai)
      thecos=cos(thetai)
      phisin=sin(phii)
      phicos=cos(phii)
      rki(1)=-thesin*phicos
      rki(2)=-thesin*phisin
      rki(3)=-thecos
      print*,'rki',rki

      if(ipol.eq.1) then
         pole(1)= thecos*phicos
         pole(2)= thecos*phisin
         pole(3)=-thesin
      end if
      if(ipol.eq.2) then
         pole(1)= -phisin
         pole(2)=  phicos
         pole(3)=  0.0
      end if
      print*,'pole',pole
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccexp(-xj*k0*(kx*xc+ky*yc+kz*zc))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ci=(0.0,1.0)
      xj=ci
      k0=20.0*pi*freq/3.0

         do i=1,mvoledge
            crhs(i)=0.
         end do

         do kip=1,bsurele
            iface=kip
            iele=bfface(4,iface)
            do i=1,3
            do j=1,3
               points(i,j)=xyz(j,bfface(i,iface))
            enddo
            enddo

            do 25 i=1,3
            do 25 j=1,3
25             pc(i,j)=(points(i,j)+points((mod(i,3)+1),j))/2.0 

            do 10 i=1,4
               ne(i)=lv(i,iele)
               co(i)=1
               xx(i)=xyz(1,ne(i))
               yy(i)=xyz(2,ne(i))
10             zz(i)=xyz(3,ne(i))

            v=0.

            do 20 i=1,4
               ii=mod(i,2)
               call detm(xx,yy,zz,i,det)
               if (ii.eq.1) then
                  a(i)=det
               else
                  a(i)=-det
               endif
               v=v+a(i)

               call detm(co,yy,zz,i,det)
               if (ii.eq.1) then
                   b(i)=-det
               else
                   b(i)=det
               endif

               call detm(co,xx,zz,i,det)
               if (ii.eq.1) then
                  c(i)=det
               else
                  c(i)=-det
               endif

               call detm(co,xx,yy,i,det)
               if (ii.eq.1) then
                  d(i)=-det
               else
                  d(i)= det
               endif

20      continue

            do 30 i=1,4
            do 30 j=1,4
               xxx=xx(i)-xx(j)
               yyy=yy(i)-yy(j)
               zzz=zz(i)-zz(j)
               length(i,j)=sqrt(xxx*xxx+yyy*yyy+zzz*zzz)
30      continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            do i=1,3
               e23(i)=xyz(i,bfface(3,iface))-xyz(i,bfface(2,iface))
               e21(i)=xyz(i,bfface(1,iface))-xyz(i,bfface(2,iface))
            enddo

            call cro_mul(e23,e21,ns)
            nsabs=sqrt(ns(1)**2+ns(2)**2+ns(3)**2)
            do i=1,3
               ns(i)=ns(i)/nsabs
            enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            call refer_trans2(points,xn,yn,length_sur)
            call generate_abc(xn,yn,as,bs,cs,area)
            area = abs(area)
cccccccccccccccccccccc compute Ni ccccccccccccccccccccccccccccccccccccc
            do  j=1,3
              r_r=0.
              do k=1,3
                r_r=r_r+rki(k)*pc(j,k)
              end do
              call elevol_N(xx,yy,zz,pc(j,1),pc(j,2),pc(j,3),
     &                      v,length,Nv,b,c,d)          
            do  i1=1,3
            do  i2=i1+1,4
              temp=0.
              call labeltran(i1,i2,i)
              call signcheck(lv(i1,iele),lv(i2,iele),2,1,ktest)

              do k=1,3
                temp1(k)=Nv(i,k)
              enddo
ccccccccccccccccccc cro_mul ccccccccccccccccccccccccccccccccccccccccccc
              call cro_mul(ns,pole,n_E)
              call cro_mul(n_E,ns,temp2)
              temp=temp+dotmul(temp1,temp2)*real(ktest)

              call cro_mul(rki,pole,n_rki)
              call cro_mul(ns,n_rki,temp2)
ccccccccccccccccc dot mul cccccccccccccccccccccccccccccccccccccccccccc
              temp=temp+dotmul(temp1,temp2)*real(ktest)

             crhs((le(i,iele)))=
     &       crhs((le(i,iele)))+
     &       temp*xj*k0*exp(-xj*k0*r_r)*area/3.0
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            enddo
            enddo
            enddo
         enddo

      RETURN
      END

