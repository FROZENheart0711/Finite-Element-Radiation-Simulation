      subroutine numint(ls,lsmax,xgl,wgl,ftheta,fphi)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  this subroutine used for define cof and wt for spherical integral  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
!      include 'mlfma_const.inc'
!      include 'mlfma_param.inc'
      integer ls, lsmax
      real    xgl(lsmax), wgl(lsmax), ftheta(2,*),fphi(2,*)

      integer it, iph, kp
      real    cst, snt, phi, csp, snp,dphi
      real    pi
      pi=3.1415926
      call dgauleg2(-1e0,1e0,xgl,wgl,ls,lsmax)

      dphi=2e0*pi/real(2*ls)

      do it=1,ls
        cst=xgl(it)
        snt=sqrt((1e0-cst)*(1e0+cst))
        ftheta(1,it)=cst
        ftheta(2,it)=snt
        do iph=1,2*ls
          kp=(it-1)*2*ls+iph
          phi=real(iph-1)*dphi
          csp=cos(phi)
          snp=sin(phi)
          fphi(1,iph)=csp
          fphi(2,iph)=snp
        enddo
      enddo

      return
      end

      subroutine dgauleg2(x1,x2,x,w,n,nd)
        implicit real*4 (a-h,o-z)
        integer n
        real*4 x1,x2,x(nd),w(nd)
        double precision eps
        parameter (eps=3.d-14)
        integer i,j,m
        double precision p1,p2,p3,pp,xl,xm,z,z1
        m=(n+1)/2
        xm=0.5d0*(x2+x1)
        xl=0.5d0*(x2-x1)
        do 12 i=1,m
          z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
  1       continue
            p1=1.d0
            p2=0.d0
            do 11 j=1,n
              p3=p2
              p2=p1
              p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
  11        continue
            pp=n*(z*p1-p2)/(z*z-1.d0)
            z1=z
            z=z1-p1/pp
          if(abs(z-z1).gt.eps)goto 1
          x(i)=xm-xl*z
          x(n+1-i)=xm+xl*z
          w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
          w(n+1-i)=w(i)
  12    continue
        return
        end