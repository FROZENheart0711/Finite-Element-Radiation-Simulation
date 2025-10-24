      subroutine vol_edge_label
     &( kl, ke, kefemm, nbnd, kei,
     &  lv, le, lnc, lpc, ldv, nc,
     &  xyz, xyzc,
     &  rp, iqof, iqfo , ip1, ip2, ip3, iq, ip)
       implicit none
      integer lv(9, *),   le(6, *), 
     &        lnc(3, *), lpc(3, *), ldv(2,*), nc(*)
       real   xyz(3, *),   xyzc(3, *)
      integer iqof(*), iqfo(*), ip1(*), ip2(*), ip3(*), ip(*), iq(*)
       real   rp(*)
      integer ni(3),type1,type2,type3
      integer mtl,kl,ke,nbnd,kei,i,j,k,k1,k2,k3,i1,i2,i3
      integer kk,rt1,rt2,rt3,kls,mvolele,l1,l2,l3,kr,j1,j2,r1,r2
      integer kefemm,kefemo,ks,keis,klc,km,kn(3)
      integer,allocatable::ipedge(:),iqb(:)
      integer,allocatable::lec(:,:),ldvc(:,:)
      integer ie1,ie2,iflag1,lo,iflag2,kf1,kf2,kff1,kff2,k4
       REAL nx,nx1,nx2,ny,ny1,ny2,nz,nz1,nz2
      integer tt,kk1,kk2,kk3
       real xi(3),norm(3),xyzc1(3)
       real xc,yc,zc,test
	mtl=4*kl

c.......define the quantity for every triangular

	do 15 i=1,kl
	do 15 j=1,4

	kk=4*(i-1)+j
        k1=mod(j,4)+1
        k2=mod(j+1,4)+1
        k3=mod(j+2,4)+1
        k1=lv(k1,i)
        k2=lv(k2,i)
        k3=lv(k3,i)
        call comi(k1,k2) 
        call comi(k1,k3) 
        call comi(k2,k3) 
	rt1=k1
	rt2=k2
	rt3=k3

	rp(kk)=real(rt1)*log(real(rt2))*log10(real(rt3))

        ip1(kk)=5

15	continue

	call hpsort(mtl,rp,iqof,iqfo)
        rp(mtl+1)=-100
c.......find the triangular on the surface

	i=1
	kls=0
        do 20 i=1, mtl
	i1=i+1
24      if (rp(i).ne.rp(i1)) goto 25
           call fcnd(mvolele, lv, iqfo(i), k1,k2,k3)
           call fcnd(mvolele, lv, iqfo(i1),l1,l2,l3)
	if ((k1.eq.l1).and.(k2.eq.l2).and.(k3.eq.l3)) then
	  ip1(iqfo(i))=0
	  ip1(iqfo(i1))=0
	else
	  i1=i1+1
	  goto 24
25        if (ip1(iqfo(i)).eq.0) goto 20
          kr=iqfo(i)-((iqfo(i)-1)/4)*4
          km=(iqfo(i)-1)/4+1
          if (lv(5+kr,km).eq.1) then
             ip1(iqfo(i))=1
          else
             ip1(iqfo(i))=2
          endif
	  kls=kls+1
  	endif
20      continue

c.......label edges on cavity

        do i=1, 6*kl
           ip(i)=0
        enddo

        mtl    = 0
 
        do 30 i=1,kl
	do 30 j=1,6
           mtl = mtl + 1
           call labelinv(j, j1, j2)
           j1=lv(j1,i)
           j2=lv(j2,i)
           call comi(j1,j2)
           r1=j1
           r2=j2
           rp(mtl) = r1 * log(real(r2))
 30     continue

	call hpsort(mtl,rp,iqof,iqfo)
        ke     = 0
        nbnd   = 0
        ke=1
        ip(1)=1
        do 201 i=2,mtl
        i1=i-1
301        if (rp(i).ne.rp(i1)) goto 401
           k=iqfo(i)-((iqfo(i)-1)/6)*6
          ie1=(iqfo(i)-1)/6+1
          call labelinv(k,k1,k2)
          k1=lv(k1,ie1)
          k2=lv(k2,ie1)
          call com(k1,k2,iflag1)
          lo=iqfo(i1)-((iqfo(i1)-1)/6)*6
          ie2=(iqfo(i1)-1)/6+1
          call labelinv(lo,l1,l2)
          l1=lv(l1,ie2)
          l2=lv(l2,ie2)
          call com(l1,l2,iflag2)
        if ((k1.eq.l1).and.(k2.eq.l2)) then
        ip(i)=ip(i1)
        else
         i1=i1-1
         if(i1.le.0) goto 401
        goto 301
401      ke=ke+1
        ip(i)=ke
        endif
201     continue
        do 50 i=1,kl
        do 50 j=1,6
        k=6*(i-1)+j
        le(j,i)=ip(iqof(k))
50      continue
        do i=1,kl
        do j=1,6
         if(le(j,i).gt.ke)then
        print*,'before error!',i,j,le(j,i)
        pause
         endif
        enddo
        enddo
ccccccccccccccccccccccccccfind edge typecccccccccccccc
       allocate(ipedge(ke))
         nbnd=0
        do i=1,ke
        ipedge(i)=0
        enddo
        do i=1,ke
        ip(i)=0
        enddo
        do i=1,ke
        ip3(i)=0
        enddo
        do i=1,kl
        do j=1,6
        call edgeface(j,kf1,kf2)
        call labelinv(j,k1,k2)
        call comi(k1,k2)
        kff1 = 4*(i-1)+kf1
        kff2 = 4*(i-1)+kf2
       if ((ip1(kff1).eq.2).or.(ip1(kff2).eq.2)) then
       ip(le(j,i))=2
       else
       endif
       enddo
       enddo
        do i=1,kl
        do j=1,6
        call edgeface(j,kf1,kf2)
        call labelinv(j,k1,k2)
        call comi(k1,k2)
c        ldv(1,le(j,i)) = lv(k1,i)
c        ldv(2,le(j,i)) =  lv(k2,i)
        kff1 = 4*(i-1)+kf1
        kff2 = 4*(i-1)+kf2
       if ((lv(5+kf1,i).eq.1).or.(lv(5+kf2,i).eq.1)) then
       if(ip(le(j,i)).eq.2)then
        nbnd=nbnd+1
        ip3(le(j,i))=nbnd
        endif
       ip(le(j,i))=1
       else
       endif
       enddo
       enddo
ccccccccccccccccccccccccresort edgescccccccccccccccccccccccccccccc
       type1=0
       type2=0
       type3=0
       do i=1,ke
       if(ip(i).eq.1)then
       type1=type1+1
        endif
       if(ip(i).eq.2)then
       type2=type2+1
        endif
       if(ip(i).eq.0)then
       type3=type3+1
       endif
       enddo
        print*,'type',type1,type2,type3
        do i=1,ke
        ipedge(i)=0
        enddo
       k1=0
       k2=0
       k3=0
       k4=0     
       do i=1,ke
       if(ip(i).eq.2)then
        k2=k2+1
        ipedge(i)=k2
        endif
       if(ip(i).eq.1)then
        if(ip3(le(j,i)).gt.0)then
        k4=k4+1
        ipedge(i)=k4+type2+type3
        iq(k4)=k4+type2+type3
        else
        k1=k1+1
        ipedge(i)=k1+type2+type3
       endif
        endif
       if(ip(i).eq.0)then
        k3=k3+1
        ipedge(i)=type2+k3
       endif
       enddo
        print*,'k1 k2 k3 k4',k1+k4,k2,k3,type1,type2,type3
        print*,'total',k1+k2+k3+k4,ke
        do i=1,ke
         if((ipedge(i).gt.ke).or.(ipedge(i).eq.0))then
        print*,'ipedge error!',ipedge(i),i,ip(i),ke
        pause
         endif
        enddo
       allocate(lec(6,kl))
       do i=1,kl
       do j=1,6
       lec(j,i)=ipedge(le(j,i))
       enddo
       enddo
        do i=1,kl
        do j=1,6
         le(j,i)=lec(j,i)
         enddo
         enddo
         deallocate(lec,ipedge)
        do i=1,kl
        do j=1,6
         if(le(j,i).gt.ke)then
        print*,'error!',i,j,le(j,i)
        pause
         endif
        enddo
        enddo
       open(32,file='le',status='unknown')
       do i=1,kl
       do j=1,6
       write(32,*) le(j,i),j,i
       enddo
       enddo
       close(32)
ccccccccccccccccccccccccget edges on cavity facecccccccccccccccccc
ccccccccccccccccccccccccget edges on PEC boundcccccccccccccccccccc
ccccccccccccccccccccccccget edges inner partccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       kefemm=type1
       ks = type2+type3
       keis= type2
       kei = type2+nbnd
       print*,'ks keis and kei',ks,keis,kei
        do i=1, kei
           if (i.le.keis) then
              ip(i)=i
           else
              ip(i)=type2+type3+i-keis
           endif
        enddo
        do i=1,kl
        do j=1,6
        if(le(j,i).gt.ks)then
        call labelinv(j,k1,k2)
        call comi(k1,k2)
        ldv(1,le(j,i)-ks)=lv(k1,i)
        ldv(2,le(j,i)-ks)=lv(k2,i)
        endif
        enddo
        enddo
        klc = 0
        mtl = 0
    	do 76 i=1,kl
	do 76 j=1,4
           if (lv(5+j,i).ne.1) goto 76
           klc = klc +1
           do k=1,3
              ni(k)=mod(j+k-1,4)+1
           enddo
ccccccccccccccctest normcccccccccccccccccccccccccc
          nx1=xyz(1,lv(ni(1),i))-xyz(1,lv(ni(2),i))
          ny1=xyz(2,lv(ni(1),i))-xyz(2,lv(ni(2),i))
          nz1=xyz(3,lv(ni(1),i))-xyz(3,lv(ni(2),i))
          nx2=xyz(1,lv(ni(3),i))-xyz(1,lv(ni(2),i))
          ny2=xyz(2,lv(ni(3),i))-xyz(2,lv(ni(2),i))
          nz2=xyz(3,lv(ni(3),i))-xyz(3,lv(ni(2),i))
          nx=ny1*nz2-nz1*ny2
          ny=nz1*nx2-nx1*nz2
          nz=nx1*ny2-ny1*nx2
          kk=lv(j,i)
          do k=1,3
          xyzc1(k)=xyz(k,kk)
          enddo
          kk1=lv(ni(1),i)
          kk2=lv(ni(2),i)
          kk3=lv(ni(3),i)
          xc=(xyz(1,kk1)+xyz(1,kk2)+xyz(1,kk3))/3.
          yc=(xyz(2,kk1)+xyz(2,kk2)+xyz(2,kk3))/3.
          zc=(xyz(3,kk1)+xyz(3,kk2)+xyz(3,kk3))/3.
          norm(1)=xc-xyzc1(1)
          norm(2)=yc-xyzc1(2)
          norm(3)=zc-xyzc1(3)
         test =nx*norm(1)+ny*norm(2)+nz*norm(3)
          if(test.gt.0)then
          tt=ni(2)
          ni(2)=ni(3)
          ni(3)=tt
          else
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccc
           do k=1,3
              mtl = mtl + 1
              k1=mod(k-1,3)+1
              k2=mod(k,  3)+1
              k1=ni(k1)
              k2=ni(k2)
              lnc(k,klc)=lv(k1,i)
              call comi(k1,k2)
              call labeltran(k1, k2, kk)
              lpc(k, klc)=le(kk,i)-ks
            enddo
76	continue
       print*,'klc============',klc,ks
 	return
	end

      subroutine area(xyz1,xyz2,xyz3,det)

      real xyz1(3),xyz2(3),xyz3(3)

      det=xyz1(1)*xyz2(2)*xyz3(3)+xyz2(1)*xyz3(2)*xyz1(3)+
     &    xyz1(2)*xyz2(3)*xyz3(1)-xyz3(1)*xyz2(2)*xyz1(3)-
     &    xyz2(1)*xyz1(2)*xyz3(3)-xyz3(2)*xyz2(3)*xyz1(1)

      return
      end
      subroutine edgeface
     &( in, kf1, kf2)
       integer in,kf1,kf2

      if (in.eq.1) then
         kf1=3
         kf2=4
      else
      endif

      if (in.eq.2) then
         kf1=2
         kf2=4
      else
      endif

      if (in.eq.3) then
         kf1=2
         kf2=3
      else
      endif

      if (in.eq.4) then
         kf1=1
         kf2=4
      else
      endif

      if (in.eq.5) then
         kf1=1
         kf2=3
      else
      endif

      if (in.eq.6) then
         kf1=1
         kf2=2
      else
      endif

      return
      end

