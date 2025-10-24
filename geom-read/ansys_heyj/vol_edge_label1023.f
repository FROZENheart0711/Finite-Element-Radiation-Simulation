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
      integer ni(3)
      integer mtl,kl,ke,nbnd,kei,i,j,k,k1,k2,k3,i1,i2,i3
      integer kk,rt1,rt2,rt3,kls,mvolele,l1,l2,l3,kr,j1,j2,r1,r2
      integer kefemm,kefemo,ks,keis,klc,km
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
           rp(mtl) = r1 * sqrt(real(r2))
 30     continue

	call hpsort(mtl,rp,iqof,iqfo)

        ke     = 0
        nbnd   = 0

        call labelfemm( mvolele, mtl, ke, nbnd, lv, rp, iqof, iqfo, 
     &                  ip1, ip, iq, ldv , xyz)

        kefemm=ke

        call labelfemo( mvolele, mtl, ke, lv, rp, iqof, iqfo, 
     &                  ip, ip1 )
        kefemo=ke-kefemm+nbnd

        call labelfemi( mvolele, mtl, ke, lv, rp, iqof, iqfo, 
     &                  ip )

       ks = ke - kefemm
       keis= kefemo-nbnd
       kei = kefemo
       
	do 8 i=1,kl
	do 8 j=1,6
           k = 6*(i-1)+j
           if (ip(k).le.kefemm) then
              le(j,i)=ip(k)+ks
           else
              le(j,i)=ip(k)-kefemm
           endif
8	continue

        do i=1, kei
           if (i.le.keis) then
              ip(i)=i
           else
              ip(i)=iq(i-keis)+ks
           endif
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
