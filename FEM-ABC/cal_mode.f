      subroutine cal_mod(freq, mvolnode, mvoledge,
     &           xyz,nsrc,wpsp,wpse,etout,
     &            ep, ur, kz)
!,
!     &           p_mnsrcsfrt_new,irank) ! newadded


      implicit none
!.....Input Data
      integer mvolnode, mvoledge
      real  xyz(3,*), freq
      integer nsrc
      integer wpsp(4,*)
      integer wpse(3,*) 
      complex ep,ur, kz
!.....Output Data
      real etout(mvoledge*2+nsrc*2)
!.....Working Variables
      real*8 k0
      integer ne(3)
      real*8  length(3)
      real*8  xx(3),yy(3),zz(3)
      real*8  area
      integer isub,kip
      integer i,j,k,ii,i1,i2,i3,j1,j2,j3,k1,k2,k3,kk,jj
      integer l1,l2
      real*8  xct,yct,zct
      real*8  xct1,yct1,zct1
      real*8  r_r,rrt
      integer ktest,ktest1,ktest2,ktest3
      real*8 pi,eTEM
      real*8  r_in,r_out,x_feed,y_feed,z_feed
      real*8 totalp
!     new added
      integer nsnode, nsedge
      integer,allocatable::ipoe(:), ipfe(:), ipon(:), ipfn(:)
      integer,allocatable::lc(:,:), lp(:,:)
      integer,allocatable::edgesign(:),nodesign(:)
      integer,allocatable::ebnd(:), nbnd(:) ! boundary
      integer,allocatable::pecnode(:)
      complex*16 att(8,8), btt(8,8), btz(8,6), bzt(6,8), bzz(6,6)
      complex*16,allocatable::at(:,:), bt(:,:)
      integer mdim
      character jobvl,jobvr
      integer n,nrhs,lda,ldb,ldvl,ldvr,lwork,info, ierr
      complex*16,allocatable :: alphar(:),alphai(:),alpha(:),
     &   beta(:),vr(:,:),work(:)
      complex*16,allocatable :: vl(:,:)
      double precision,allocatable:: rwork(:)
      complex*16,allocatable :: kz1(:)
      integer,allocatable :: evindex(:)
      complex*16,allocatable :: kz2(:)
      real*8,allocatable:: kz2_r(:)
      real*8,allocatable:: et(:)
      complex*16,allocatable :: att_f(:,:)
      complex*16,allocatable :: btt_1(:,:),btz_1(:,:),bzt_1(:,:)
      complex*16,allocatable :: bzz_1(:,:)
      integer,allocatable :: ipiv(:)
      complex*16,allocatable :: idm(:,:)
      complex*16,allocatable :: btt_2(:,:),btt_3(:,:),btt_f(:,:)
      integer,allocatable :: ipfo(:),ipof(:)

      real*8,allocatable :: exx(:),eyy(:),ezz(:)
      integer level
      real*8 exx_temp,eyy_temp,ezz_temp
      integer*8 p_srcspfrt
      integer*8 p_mvolnodefrt
      real*8 tempp,maxp
      real*8 tempmul
      integer order
       real ep0, u0, rk0
       double precision var(7), varxyz(3,7)
       real*8 tan
       integer totalsrc
       real kzreal,kzaimag
       pi=atan2(1.0,1.0)*4.0
       ep0=8.8541853e-4;
       u0=pi*40.0;
       rk0=2.0*pi*freq*sqrt(ep0*u0)*10.0;
       k0=rk0

!      integer p_mnsrcsfrt_new(*)
!      integer irank
c.......Begin
!     renumber edge and nodes
c      k0=9
      level=1
c      level=1000
      
      totalsrc=nsrc
      
c      print*,"nsrc=",nsrc
c      print*,"mvoledge=",mvoledge 
c      print*,"mvolnode=",mvolnode     

      allocate(edgesign(mvoledge))
      allocate(ipoe(mvoledge), ipfe(mvoledge))
c      print*,"11"
      edgesign=0
      ipoe=0
      ipfe=0
c      print*,"p_srcspfrt=",p_srcspfrt
      do i=1, nsrc
        do j=1,3
c         print*,"i,j,p_srcspfrt=",i,j,p_srcspfrt
c         if(irank.eq.1)then 
c         print*,"wpse=",irank,wpse(j,i)
c         endif
         edgesign(wpse(j,i))=
     &           edgesign(wpse(j,i))+1
        enddo
      enddo

      nsedge=0
      do i=1,mvoledge
       if(edgesign(i).gt.0)then
         nsedge=nsedge+1
         ipoe(i)=nsedge
         ipfe(nsedge)=i
       endif
      enddo
c      print*,"nsedge=",nsedge
      allocate(ebnd(nsedge))
      ebnd=0
c      print*,"11"
      do i=1,nsedge
        if(edgesign(ipfe(i)).eq.1)then
          ebnd(i)=1
        endif
      enddo
c      print*,"22" 
      

      deallocate(edgesign)
      allocate(nodesign(mvolnode))
      allocate(ipon(mvolnode), ipfn(mvolnode))
      nodesign=0
      ipon=0
      ipfn=0
c      print*,"nsrc=",nsrc
      do i=1, nsrc
        do j=1,3
c         if(irank.eq.1)then
c          print*,"wpsp",i,p_srcspfrt
c,wpsp(j,i)
c         endif
         nodesign(wpsp(j,i))=1
        enddo
      enddo
      
      nsnode=0
      do i=1,mvolnode
       if(nodesign(i).eq.1)then
         nsnode=nsnode+1
         ipon(i)=nsnode
         ipfn(nsnode)=i
       endif
      enddo
c      print*,"nsnode=",nsnode      
      allocate(nbnd(nsnode))
      nbnd=0
c      print*,"irank=",irank

      do ii=1,nsrc
        do i=1,3
          k1=i
          k2=mod(i,3)+1
          if(ebnd(ipoe(wpse(i,ii))).eq.1)then
            nbnd(ipon(wpsp(k1,ii)))=1
            nbnd(ipon(wpsp(k2,ii)))=1
          endif
        enddo
      enddo

c      print*,'ebnd:'
      do i=1,nsedge
        if(ebnd(i).gt.0)then
c        print*,ebnd(i),i
        endif
      enddo

      deallocate(nodesign)
c      print*,"irank=",irank
      mdim=nsedge+nsnode
c      print*,"mdim=",mdim,nsedge,nsnode
c      do i=1,nsnode
c        if(nbnd(i).eq.1)then
c          print*,xyz(1,ipfn(i)),xyz(2,ipfn(i)),xyz(3,ipfn(i))
c        endif
c      enddo
c      stop

      allocate(att_f(nsedge*2+nsrc*2,nsedge*2+nsrc*2))
      allocate(btt_1(nsedge*2+nsrc*2,nsedge*2+nsrc*2))
      allocate(btz_1(nsedge*2+nsrc*2,nsnode+nsedge))
      allocate(bzt_1(nsnode+nsedge,nsedge*2+nsrc*2))
      allocate(bzz_1(nsnode+nsedge,nsnode+nsedge))
      att_f=0.
      btt_1=0.
      btz_1=0.
      bzt_1=0.
      bzz_1=0.


c      stop

!     begin
c      print*,'calculating begin'
      do ii=1,nsrc
        call elesrc(ii,xyz,p_srcspfrt,p_mvolnodefrt,wpsp,
     &     k0,ur,ep,
     &     att,btt,btz,bzt,bzz)
        

        
         call assemsrc(ii,att,btt,btz,bzt,bzz,att_f,btt_1,
     &          btz_1,bzt_1,bzz_1,nsrc,nsedge,nsnode,p_mvolnodefrt,
     &          wpse,wpsp,p_srcspfrt,mvoledge,mvolnode,ipoe,ipon,
     &          ebnd,nbnd)
      enddo


c      print*,'k0:',k0
c      print*,'ur:',ur
c      print*,'ep:',ep

c      print*,'assem matrix finished!'

      order=1
      if(order.eq.0)then
      do i=nsedge+1,nsedge*2+nsrc*2
        do j=nsedge+1,nsedge*2+nsrc*2
          if(i.eq.j)then
            att_f(i,j)=1.
            btt_1(i,j)=1.
          else
            att_f(i,j)=0.
            btt_1(i,j)=0.
          endif
        enddo
      enddo

      do i=nsedge+1,nsedge*2+nsrc*2
        do j=1,nsedge          
            att_f(i,j)=0.
            btt_1(i,j)=0.
            att_f(j,i)=0.
            btt_1(j,i)=0.          
        enddo
      enddo

      do i=nsnode+1,nsnode+nsedge
        do j=1,nsedge*2+nsrc*2
            bzt_1(i,j)=0.
            btz_1(j,i)=0.
        enddo
      enddo

      do i=1,nsnode
        do j=nsedge+1,nsedge*2+nsrc*2
          bzt_1(i,j)=0.
          btz_1(j,i)=0.
        enddo
      enddo

      do i=nsnode+1,nsnode+nsedge
        do j=nsnode+1,nsnode+nsedge
            if(i.eq.j)then
              bzz_1(i,j)=1.
            else
              bzz_1(i,j)=0.
            endif
        enddo
      enddo

      do i=nsnode+1,nsnode+nsedge
        do j=1,nsnode
            bzz_1(i,j)=0.
            bzz_1(j,i)=0.
        enddo
      enddo
      endif

      n=nsnode+nsedge
      nrhs=nsnode+nsedge
      lda=nsnode+nsedge
      ldb=nsnode+nsedge
      allocate(ipiv(n))
      allocate(idm(nsnode+nsedge,nsnode+nsedge))
      idm=0
      do i=1,nsnode+nsedge
        idm(i,i)=1.
      enddo



      call zgesv(n,nrhs,bzz_1,lda,ipiv,idm,ldb,info)

      allocate(btt_2(nsedge*2+nsrc*2,nsnode+nsedge))

      do i=1,nsedge*2+nsrc*2
           do j=1,nsnode+nsedge
               btt_2(i,j)=0
               do k=1,nsnode+nsedge
                   btt_2(i,j)=btt_2(i,j)+btz_1(i,k)*idm(k,j)
               enddo
           enddo
      enddo 
      
      allocate(btt_3(nsedge*2+nsrc*2,nsedge*2+nsrc*2))

      do i=1,nsedge*2+nsrc*2
           do j=1,nsedge*2+nsrc*2
               btt_3(i,j)=0
               do k=1,nsnode+nsedge
                   btt_3(i,j)=btt_3(i,j)+btt_2(i,k)*bzt_1(k,j)
               enddo
           enddo
      enddo

      allocate(btt_f(nsedge*2+nsrc*2,nsedge*2+nsrc*2))

      do i=1,nsedge*2+nsrc*2
        do j=1,nsedge*2+nsrc*2
          btt_f(i,j)=btt_3(i,j)-btt_1(i,j)
        enddo
      enddo


      jobvl='n'
      jobvr='v'
!n=maxsys2
      n=nsedge*2+nsrc*2
!lda=maxsys2
      lda=nsedge*2+nsrc*2
!ldb=maxsys2
      ldb=nsedge*2+nsrc*2
      ldvl=1
!ldvr=maxsys2
      ldvr=nsedge*2+nsrc*2
!lwork=160*maxsys2
      lwork=160*(nsedge*2+nsrc*2)

!allocate(alphar(maxsys2),alphai(maxsys2))
!allocate(beta(maxsys2),vr(maxsys2,maxsys2),work(lwork))
      allocate(alphar(nsedge*2+nsrc*2),alphai(nsedge*2+nsrc*2))
      allocate(alpha(n))
      allocate(beta(nsedge*2+nsrc*2),vr(nsedge*2+nsrc*2,
     &  nsedge*2+nsrc*2),work(lwork)) !
      allocate(rwork(8*n))
      allocate(vl(ldvl,n))
      alphar=0
      alphai=0
      alpha=0
      beta=0
      vr=0
      work=0
      vl=0
    
!call sggev(jobvl,jobvr,n,a3,lda,b3,ldb,alphar,alphai,beta,vl,ldvl,vr,ldvr,work,lwork,info)
c      call dggev(jobvl,jobvr,n,att_f,lda,btt_f,ldb,alphar,
c     &  alphai,beta,vl,ldvl,vr,ldvr,work,lwork,info)


      call zggev(jobvl, jobvr, n, att_f, lda, btt_f, ldb, alpha, 
     &    beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info)

c      print*,'info:',info

c      do j=1,nsedge*2+nsrc*2
c        do i=1,nsedge*2+nsrc*2
c          if(vr(i,j).ne.0)then
c            print*,'vr(',i,j,')',vr(i,j)
c          endif
c        enddo
c      enddo
c      stop

      allocate(kz1(nsedge*2+nsrc*2),evindex(nsedge*2+nsrc*2))
      kz1=0
      do i=1,nsedge*2+nsrc*2
        if(zabs(beta(i)).ge.1e-12)then
c          if(alphai(i).eq.0)then
          kz1(i)=alpha(i)/beta(i)

c          endif
        endif
c        print*,alpha(i),beta(i),i
c        print*,''
        evindex(i)=i
      enddo
     
      allocate(kz2(nsedge*2+nsrc*2))
      allocate(kz2_r(nsedge*2+nsrc*2))

      do i=1,nsedge*2+nsrc*2
        kz2(i)=cdsqrt(kz1(i))
        kz2_r(i)=real(kz2(i))
c        print*,kz2(i)
      enddo


      allocate(ipfo(nsedge*2+nsrc*2),ipof(nsedge*2+nsrc*2))

      call dhpsort(nsedge*2+nsrc*2,kz2_r,ipof,ipfo)

      
      kzreal=real(kz2(ipfo(nsedge*2+nsrc*2)))
      kzaimag=aimag(kz2(ipfo(nsedge*2+nsrc*2)))
      print*,'kzreal and aimag:',kzreal,kzaimag
      kz=cmplx(kzreal,kzaimag)
      print*,'kz and k0:',kz,k0
  

      allocate(et(nsedge*2+nsrc*2))



      do i=1,nsedge*2+nsrc*2
        if(real(vr(i,ipfo(nsedge*2+nsrc*2))).ne.0)then
          tan=aimag(vr(i,ipfo(nsedge*2+nsrc*2)))/
     &    real(vr(i,ipfo(nsedge*2+nsrc*2)))
          et(i)=real(vr(i,ipfo(nsedge*2+nsrc*2)))*
     &      dsqrt(1+tan**2)
        else
          et(i)=0.
        endif
c        print*,''
c        print*,real(vr(i,ipfo(nsedge*2+nsrc*2))),
c     &         aimag(vr(i,ipfo(nsedge*2+nsrc*2)))
c        print*,vr(i,ipfo(nsedge*2+nsrc*2)),i
c        print*,et(i),i
      enddo
c      print*,vr(570,ipfo(nsedge*2+nsrc*2))
c      stop
      etout=0

      do i=1,nsedge
        etout(ipfe(i))=et(i)
      enddo

      do i=1,nsedge
        etout(ipfe(i)+mvoledge)=et(i+nsedge)
      enddo

      do i=1,nsrc
        etout(i+mvoledge*2)=et(i+nsedge*2)
      enddo

      do i=1,nsrc
        etout(i+mvoledge*2+nsrc)=et(i+nsedge*2+nsrc)
      enddo


      allocate(exx(nsrc),eyy(nsrc),ezz(nsrc))

      
      totalp=0.
      maxp=0.
      kk=0
      xct1=0.
      yct1=0.
      zct1=0.
      jj=0


      open(111,file='portE.txt',status='unknown')  
      do ii=1,nsrc

        do i=1,3
            ne(i)=wpsp(i,ii)
            xx(i)=xyz(1,ne(i))*level
            yy(i)=xyz(2,ne(i))*level
            zz(i)=xyz(3,ne(i))*level
        end do



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       var(1) = 9./20.
       var(2) = 1./20.
       var(3) = 1./20.
       var(4) = 1./20.
       var(5) = 2./15.
       var(6) = 2./15.
       var(7) = 2./15.

         varxyz(1,1) = (xx(1)+xx(2)+xx(3))/3.
         varxyz(1,2) =  xx(1)
         varxyz(1,3) =  xx(2)
         varxyz(1,4) =  xx(3)

       varxyz(1,5) = (xx(1)+xx(2))/2.
       varxyz(1,6) = (xx(2)+xx(3))/2.
       varxyz(1,7) = (xx(3)+xx(1))/2.

         varxyz(2,1) = (yy(1)+yy(2)+yy(3))/3.
         varxyz(2,2) =  yy(1)
         varxyz(2,3) =  yy(2)
         varxyz(2,4) =  yy(3)

       varxyz(2,5) = (yy(1)+yy(2))/2.
       varxyz(2,6) = (yy(2)+yy(3))/2.
       varxyz(2,7) = (yy(3)+yy(1))/2.

         varxyz(3,1) = (zz(1)+zz(2)+zz(3))/3.
         varxyz(3,2) =  zz(1)
         varxyz(3,3) =  zz(2)
         varxyz(3,4) =  zz(3)

       varxyz(3,5) = (zz(1)+zz(2))/2.
       varxyz(3,6) = (zz(2)+zz(3))/2.
       varxyz(3,7) = (zz(3)+zz(1))/2.
      xct1=xct1+varxyz(1,1)
      yct1=yct1+varxyz(2,1)
      zct1=zct1+varxyz(3,1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do j=1,7
        
        xct= varxyz(1,j) 
        yct= varxyz(2,j) 
        zct= varxyz(3,j) 

        call dareapatch(xx,yy,zz,length,area)

        area=area/2.
        

        call cal_et2(xyz,xx,yy,zz,xct,yct,zct,
     &   ii,nsrc,mvoledge,wpsp,wpse,etout,
     &   exx_temp,eyy_temp,ezz_temp,area,
     &   p_srcspfrt,p_mvolnodefrt, totalsrc)
c        print*,"finish cal_et2"
        exx(ii)=exx_temp
        eyy(ii)=eyy_temp
        ezz(ii)=ezz_temp
c        print*,exx(ii),eyy(ii),ezz(ii)
        write(111,*) xct,yct,zct
        write(111,*) exx_temp,eyy_temp,ezz_temp

        call dareapatch(xx,yy,zz,length,area)
        area=area/2.


        totalp=totalp+(exx(ii)**2+eyy(ii)**2+ezz(ii)**2)*area*var(j)
        

        tempp=exx(ii)**2+eyy(ii)**2+ezz(ii)**2
        if(tempp.gt.maxp)then
c         print*,'jj,j',jj,j
          maxp=tempp
          kk=ii
          jj=j
        endif
         enddo
      enddo
      close(111)

      xct1=xct1/nsrc
      yct1=yct1/nsrc
      zct1=zct1/nsrc


      totalp=dsqrt(totalp)
      exx=exx/totalp
      eyy=eyy/totalp
      ezz=ezz/totalp

      et=et/totalp

      etout=etout/totalp
      
      tempmul=0.
    
        do i=1,3
            ne(i)=wpsp(i,kk)
            xx(i)=xyz(1,ne(i))*level
            yy(i)=xyz(2,ne(i))*level
            zz(i)=xyz(3,ne(i))*level
        end do
         varxyz(1,1) = (xx(1)+xx(2)+xx(3))/3.
         varxyz(1,2) =  xx(1)
         varxyz(1,3) =  xx(2)
         varxyz(1,4) =  xx(3)

       varxyz(1,5) = (xx(1)+xx(2))/2.
       varxyz(1,6) = (xx(2)+xx(3))/2.
       varxyz(1,7) = (xx(3)+xx(1))/2.

         varxyz(2,1) = (yy(1)+yy(2)+yy(3))/3.
         varxyz(2,2) =  yy(1)
         varxyz(2,3) =  yy(2)
         varxyz(2,4) =  yy(3)

       varxyz(2,5) = (yy(1)+yy(2))/2.
       varxyz(2,6) = (yy(2)+yy(3))/2.
       varxyz(2,7) = (yy(3)+yy(1))/2.

         varxyz(3,1) = (zz(1)+zz(2)+zz(3))/3.
         varxyz(3,2) =  zz(1)
         varxyz(3,3) =  zz(2)
         varxyz(3,4) =  zz(3)

       varxyz(3,5) = (zz(1)+zz(2))/2.
       varxyz(3,6) = (zz(2)+zz(3))/2.
       varxyz(3,7) = (zz(3)+zz(1))/2.

      xct= varxyz(1,jj)
      yct= varxyz(2,jj)
      zct= varxyz(3,jj)

      tempmul=exx(kk)*(xct-xct1)+
     &        eyy(kk)*(yct-yct1)+
     &        ezz(kk)*(zct-zct1)
   
 
      if(tempmul.gt.0)then
        ktest=1
      else
        ktest=-1
      endif
     
      etout=etout*ktest
c      open(9999,file='Etout.txt',status='unknown')  
c      do i=1,mvoledge*2+totalsrc*2
c        if(etout(i).ne.0)then
c          print*,etout(i),i
c          write(9999,*)etout(i),i
      
c        endif
c       enddo
c      close(9999)
c      stop

      exx=exx*ktest
      eyy=eyy*ktest
      ezz=ezz*ktest







      j=0
      if(j.eq.1)then

       do ii=1,nsrc
        do i=1,3
          ne(i)=wpsp(i,ii)
          xx(i)=xyz(1,ne(i))*level
          yy(i)=xyz(2,ne(i))*level
          zz(i)=xyz(3,ne(i))*level
        end do

           xct=(xx(1)+xx(2)+xx(3))/3.0
           yct=(yy(1)+yy(2)+yy(3))/3.0
           zct=(zz(1)+zz(2)+zz(3))/3.0

            r_r=dsqrt((xct-xct1)**2+(yct-yct1)**2
     &               +(zct-zct1)**2)
           

            eTEM=1/dsqrt(2.0*pi*log(3.49e-3/6e-4))/r_r

            exx_temp=(xct-xct1)/r_r*eTEM
            eyy_temp=(yct-yct1)/r_r*eTEM
            ezz_temp=(zct-zct1)/r_r*eTEM
            print*,xct,yct,zct
            print*,exx(ii),eyy(ii),ezz(ii)
            print*,exx_temp,eyy_temp,ezz_temp
            print*,''

        

       enddo
      endif


      deallocate(ipoe,ipfe,ebnd,ipon,ipfn,nbnd,
     &           att_f,btt_1,btz_1,bzt_1,bzz_1,ipiv,
     &           idm,btt_2,btt_3,btt_f,alphar,alphai,
     &           beta,vr,work,kz1,evindex,ipfo,ipof,
     &           et,exx,eyy,ezz,alpha,vl,rwork,kz2)

      print*,'finish cal_mode'

!.....test
     


      return
      end

