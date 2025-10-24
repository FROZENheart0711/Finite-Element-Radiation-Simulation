      subroutine farfld(thetas,phis,crhs,freq,
     &           mvolnode,mvolele,mvoledge,
     &           xyz,lv,le,
     &           sfface,ssurele,totalp)

      implicit none

C.....Input Data
      integer dd,mvolnode,mvolele,mvoledge,dedgetotal
      real    freq,thetas,phis
      real    xyz(3,mvolnode)
      integer lv(5,mvolele),le(6,mvolele)
      integer sfface(4,4*mvolele),ssurele
      complex crhs(mvoledge)
      real    totalp
    
c.....Working Variables
      real,parameter::pi=3.1415926
      real    k0,rad
      complex ci
      real    thesin,thecos,phisin,phicos
      integer i, j, k, ii,i1, i2
      real    rks(3), the(3), phi(3), rc(3), bn(3)
      real    length(4,4)
      integer ne(4)
      real    co(4),xx(4),yy(4),zz(4)
      real    v, det, xxx, yyy, zzz,length_sur(3),area
      real    a(4),b(4),c(4),d(4),e23(3),e21(3),ns(3)
      real    nsabs
      real    as(3),bs(3),cs(3)
      real    xn(3),yn(3),points(3,3),pc(3,3)
      real    Nv(6,3),n_Nv(3)
      real    delNv(3),n_delNv(3),temp1(3),temp2(3)
      complex xj,Ig(3),cdotmul
      integer kip,iface,iele,ktest
      real    r_r
      real    rn_delNv(3),rrn_delNv(3)
      real    sigma,sigma_dB,sigma2,sigma_dB2
ccccccccccccccc begin ccccccccccccccccccccccccccccccccccccccccccccccccc
      thesin=sin(thetas)
      thecos=cos(thetas)
      phisin=sin(phis)
      phicos=cos(phis)

      rks(1)=thesin*phicos
      rks(2)=thesin*phisin
      rks(3)=thecos

      the(1)= thecos*phicos
      the(2)= thecos*phisin
      the(3)=-thesin

      phi(1)=-phisin
      phi(2)= phicos
      phi(3)= 0.0

      ci=(0.,1.)
      xj=ci
      k0=20.*pi*freq/3.
      rad=pi/180.
      do i=1,3
         Ig(i)=0.
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do kip=1,ssurele

            iface=kip
            do i=1,3
            do j=1,3
               points(i,j)=xyz(j,sfface(i,iface))
            enddo
            enddo

            do 25 i=1,3
            do 25 j=1,3
25             pc(i,j)=(points(i,j)+points((mod(i,3)+1),j))/2.0 

            iele=sfface(4,iface)
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
               e23(i)=xyz(i,sfface(3,iface))-xyz(i,sfface(2,iface))
               e21(i)=xyz(i,sfface(1,iface))-xyz(i,sfface(2,iface))
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
cccccccccccc...................cccccccccccccccc

            do 40 j=1,3
               r_r=0.
               do k=1,3
                 r_r=r_r+rks(k)*pc(j,k)
               enddo
            do 40 i1=1,3
            do 40 i2=i1+1,4

              delNv(1)=2*length(i1,i2)*(c(i1)*d(i2)-c(i2)*d(i1))/v**2
              delNv(2)=-2*length(i1,i2)*(b(i1)*d(i2)-b(i2)*d(i1))/v**2
              delNv(3)=2*length(i1,i2)*(b(i1)*c(i2)-b(i2)*c(i1))/v**2
              call cro_mul(ns,delNv,n_delNv)
              call cro_mul(n_delNv,rks,rn_delNv)
              call cro_mul(rks,rn_delNv,rrn_delNv)
              call labeltran(i1,i2,i)
              call signcheck(lv(i1,iele),lv(i2,iele),2,1,ktest)

              if (ktest.gt.0) then

              do k=1,3
               Ig(k)=Ig(k)+rrn_delNv(k)*crhs(le(i,iele))
     &        *exp(xj*k0*r_r)*area/3.0
              enddo
              else
              
              do k=1,3
                Ig(k)=Ig(k)-rrn_delNv(k)*crhs(le(i,iele))
     &        *exp(xj*k0*r_r)*area/3.0
              enddo
              endif
40          continue
ccccccccc


            do 45 j=1,3                             !.....n*E*r
               r_r=0.
               do k=1,3
                  r_r=r_r+rks(k)*pc(j,k)
               end do
               call elevol_N(xx,yy,zz,pc(j,1),pc(j,2),pc(j,3),
     &                      v,length,Nv,b,c,d)
            do 45 i1=1,3
            do 45 i2=i1+1,4
               call labeltran(i1,i2,i)
               do k=1,3
                  temp1(k)=Nv(i,k)
               enddo

               call cro_mul(ns,temp1,n_Nv)
               call cro_mul(n_Nv,rks,temp2)

c              do k=1,3
c                n_Nv_r(i,k)=temp2(k)
c              enddo
              call signcheck(lv(i1,iele),lv(i2,iele),2,1,ktest)

              if (ktest.gt.0) then
         
              do k=1,3
                Ig(k)=Ig(k)+xj*k0*temp2(k)*crhs(le(i,iele))
     &         *exp(xj*k0*r_r)*area/3.0
              enddo
              else

              do k=1,3
                 Ig(k)=Ig(k)-xj*k0*temp2(k)*crhs(le(i,iele))
     &         *exp(xj*k0*r_r)*area/3.0
              enddo
            endif

45       continue

         enddo
        sigma=cabs(cdotmul(the,ig))**2!/4.0/pi
        sigma2=cabs(cdotmul(phi,ig))**2!/4.0/pi
        sigma_dB=10*log10(sigma/totalp)
        sigma_dB2=10*log10(sigma2/totalp)

      write(10,*) thetas/rad,sigma_dB,sigma_dB2

      RETURN
      END

      subroutine farfld1(thetas,phis,crhs,freq,
     &           mvolnode,mvolele,mvoledge,
     &           xyz,lv,le,
     &           sfface,ssurele,
     &           totalp,cweigh)

      implicit none

C.....Input Data
      integer dd,mvolnode,mvolele,mvoledge,dedgetotal
      real    freq,thetas,phis
      real    xyz(3,mvolnode)
      integer lv(5,mvolele),le(6,mvolele)
      integer sfface(4,4*mvolele),ssurele
      complex crhs(mvoledge)
      real    totalp,cweigh
    
c.....Working Variables
      real,parameter::pi=3.1415926
      real    k0,rad
      complex ci
      real    thesin,thecos,phisin,phicos
      integer i, j, k, ii,i1, i2
      real    rks(3), the(3), phi(3), rc(3), bn(3)
      real    length(4,4)
      integer ne(4)
      real    co(4),xx(4),yy(4),zz(4)
      real    v, det, xxx, yyy, zzz,length_sur(3),area
      real    a(4),b(4),c(4),d(4),e23(3),e21(3),ns(3)
      real    nsabs
      real    as(3),bs(3),cs(3)
      real    xn(3),yn(3),points(3,3),pc(3,3)
      real    Nv(6,3),n_Nv(3)
      real    delNv(3),n_delNv(3),temp1(3),temp2(3)
      complex xj,Ig(3),cdotmul
      integer kip,iface,iele,ktest
      real    r_r
      real    rn_delNv(3),rrn_delNv(3)
      real    sigma,sigma_dB,sigma2,sigma_dB2
ccccccccccccccc begin ccccccccccccccccccccccccccccccccccccccccccccccccc
      thesin=sin(thetas)
      thecos=cos(thetas)
      phisin=sin(phis)
      phicos=cos(phis)

      rks(1)=thesin*phicos
      rks(2)=thesin*phisin
      rks(3)=thecos

      the(1)= thecos*phicos
      the(2)= thecos*phisin
      the(3)=-thesin

      phi(1)=-phisin
      phi(2)= phicos
      phi(3)= 0.0

      ci=(0.,1.)
      xj=ci
      k0=20.*pi*freq/3.
      rad=pi/180.
      do i=1,3
         Ig(i)=0.
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do kip=1,ssurele

            iface=kip
            do i=1,3
            do j=1,3
               points(i,j)=xyz(j,sfface(i,iface))
            enddo
            enddo

            do 25 i=1,3
            do 25 j=1,3
25             pc(i,j)=(points(i,j)+points((mod(i,3)+1),j))/2.0 

            iele=sfface(4,iface)
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
               e23(i)=xyz(i,sfface(3,iface))-xyz(i,sfface(2,iface))
               e21(i)=xyz(i,sfface(1,iface))-xyz(i,sfface(2,iface))
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
cccccccccccc...................cccccccccccccccc

            do 40 j=1,3
               r_r=0.
               do k=1,3
                 r_r=r_r+rks(k)*pc(j,k)
               enddo
            do 40 i1=1,3
            do 40 i2=i1+1,4

              delNv(1)=2*length(i1,i2)*(c(i1)*d(i2)-c(i2)*d(i1))/v**2
              delNv(2)=-2*length(i1,i2)*(b(i1)*d(i2)-b(i2)*d(i1))/v**2
              delNv(3)=2*length(i1,i2)*(b(i1)*c(i2)-b(i2)*c(i1))/v**2
              call cro_mul(ns,delNv,n_delNv)
              call cro_mul(n_delNv,rks,rn_delNv)
              call cro_mul(rks,rn_delNv,rrn_delNv)
              call labeltran(i1,i2,i)
              call signcheck(lv(i1,iele),lv(i2,iele),2,1,ktest)

              if (ktest.gt.0) then

              do k=1,3
               Ig(k)=Ig(k)+rrn_delNv(k)*crhs(le(i,iele))
     &        *exp(xj*k0*r_r)*area/3.0
              enddo
              else
              
              do k=1,3
                Ig(k)=Ig(k)-rrn_delNv(k)*crhs(le(i,iele))
     &        *exp(xj*k0*r_r)*area/3.0
              enddo
              endif
40          continue
ccccccccc


            do 45 j=1,3                             !.....n*E*r
               r_r=0.
               do k=1,3
                  r_r=r_r+rks(k)*pc(j,k)
               end do
               call elevol_N(xx,yy,zz,pc(j,1),pc(j,2),pc(j,3),
     &                      v,length,Nv,b,c,d)
            do 45 i1=1,3
            do 45 i2=i1+1,4
               call labeltran(i1,i2,i)
               do k=1,3
                  temp1(k)=Nv(i,k)
               enddo

               call cro_mul(ns,temp1,n_Nv)
               call cro_mul(n_Nv,rks,temp2)

c              do k=1,3
c                n_Nv_r(i,k)=temp2(k)
c              enddo
              call signcheck(lv(i1,iele),lv(i2,iele),2,1,ktest)

              if (ktest.gt.0) then
         
              do k=1,3
                Ig(k)=Ig(k)+xj*k0*temp2(k)*crhs(le(i,iele))
     &         *exp(xj*k0*r_r)*area/3.0
              enddo
              else

              do k=1,3
                 Ig(k)=Ig(k)-xj*k0*temp2(k)*crhs(le(i,iele))
     &         *exp(xj*k0*r_r)*area/3.0
              enddo
            endif

45       continue

         enddo
        sigma=cabs(cdotmul(the,ig))**2/4.0/pi
        sigma2=cabs(cdotmul(phi,ig))**2/4.0/pi
         totalp=totalp+(sigma+sigma2)*cweigh

      RETURN
      END

        subroutine cro_mul(a,b,c)
           
        IMPLICIT NONE

c.......Input Data

        real   a(3),b(3)

c.......Output Data

        real   c(3)


        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)


        return
        end

      subroutine refer_trans2(points,xn,yn,length)

      implicit none

c.....Input Data
      real points(3,3)

c.....Output Data
      real xn(3),yn(3),length(3)

        length(1)=sqrt((points(1,1)-points(2,1))**2
     &                +(points(1,2)-points(2,2))**2
     &                +(points(1,3)-points(2,3))**2)
        length(2)=sqrt((points(2,1)-points(3,1))**2
     &                +(points(2,2)-points(3,2))**2
     &                +(points(2,3)-points(3,3))**2)
        length(3)=sqrt((points(1,1)-points(3,1))**2
     &                +(points(1,2)-points(3,2))**2
     &                +(points(1,3)-points(3,3))**2)

        xn(1)=0.
        yn(1)=0.
        xn(2)=length(1)
        yn(2)=0.
        xn(3)=(length(1)**2+length(3)**2-length(2)**2)/
     &        length(1)/2.
        yn(3)=sqrt(length(3)**2-xn(3)**2)


        RETURN
        END

      subroutine generate_abc(xn,yn,a,b,c,area)

      implicit none

c.....Input Data
      real    xn(3),yn(3)

c.....Output Data

      real    a(3),b(3),c(3),area

c.....Working Variables
      integer i,k1,k2
 


        do 10 i=1,3
        k1=mod(i,3)+1
        k2=mod(i+1,3)+1
        a(i)=xn(k1)*yn(k2)-xn(k2)*yn(k1)
        b(i)=yn(k1)-yn(k2)
10      c(i)=xn(k2)-xn(k1)

        area=0.5*(b(1)*c(2)-b(2)*c(1))

        RETURN
        END

        subroutine comi(k1,k2)
        if (k1.gt.k2) then
        k=k1
        k1=k2
        k2=k
        else
        endif
        return
        end


