      subroutine vol_edge_label
     &( kl, ke, klc,kefemm, nbnd, kei,
     &  lv, le, lnc, lpc, ldv, nc,
     &  xyz, xyzc,
     &  rp, iqof, iqfo , ip1, ip2, ip3, iq, ip,
     &  kefemabc, kefemsrc, labelabc, labelsrc)

      integer lv(9, *),   le(6, *), 
     &        lnc(4, *), lpc(3, *), ldv(2,*), nc(*)
       real   xyz(3, *),   xyzc(3, *)
      integer iqof(*), iqfo(*), ip1(*), ip2(*), ip3(*), ip(*), iq(*)
       real   rp(*)
      integer ni(3),tt
       real xi(3),norm(3),xyzc1(3)
       REAL nx,nx1,nx2,ny,ny1,ny2,nz,nz1,nz2
       mtl=4*kl

c.......define the quantity for every triangular

       do 15 i=1,kl
       do 15 j=1,4

        kk=4*(i-1)+j
!        k1=mod(j,4)+1
!        k2=mod(j+1,4)+1
!        k3=mod(j+2,4)+1
        if (j.eq.1) then
           k1=1
           k2=3
           k3=2
        else if (j.eq.2) then
           k1=1
           k2=2
           k3=4
        else if (j.eq.3) then
           k1=1
           k2=4
           k3=3
        else if (j.eq.4) then
           k1=2
           k2=3
           k3=4
        end if
           
        k1=lv(k1,i)
        k2=lv(k2,i)
        k3=lv(k3,i)
        call comi(k1,k2) 
        call comi(k1,k3) 
        call comi(k2,k3) 
        rt1=k1
        rt2=k2
        rt3=k3

        rp(kk)=rt1*log(rt2)*log10(rt3)

        ip1(kk)=5

15      continue

        call hpsort(mtl,rp,iqof,iqfo)

c.......find the triangular on the surface

        i=1
        kls=0
        do 20 i=1, mtl
           i1=i+1
24         if (rp(i).ne.rp(i1)) goto 25
           call fcnd(mvolele, lv, iqfo(i), k1,k2,k3)
           call fcnd(mvolele, lv, iqfo(i1),l1,l2,l3)
           if ((k1.eq.l1).and.(k2.eq.l2).and.(k3.eq.l3)) then
              ip1(iqfo(i))=0
              ip1(iqfo(i1))=0
              goto 25
           else
              i1=i1+1
              goto 24
!25            if (ip1(iqfo(i)).eq.0) goto 20
25            kr=iqfo(i)-((iqfo(i)-1)/4)*4
              km=(iqfo(i)-1)/4+1
              if (lv(5+kr,km).eq.labelabc) then
                ip1(iqfo(i))=1
              elseif(lv(5+kr,km).eq.labelsrc) then  !new added (src face)
                ip1(iqfo(i))=11
              elseif(lv(5+kr,km).gt.0) then
                ip1(iqfo(i))=2
              endif
              kls=kls+1
           endif
20      continue

        do i=1,mtl
            if(ip1(i).eq.5)then
               print*,'ip1 error!!!'
               stop
            endif
        enddo

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
           rp(mtl) = r1 * sqrt(r2)
 30     continue

        call hpsort(mtl,rp,iqof,iqfo)

        ke     = 0
        nbnd   = 0

        call labelfemm( mvolele, mtl, ke, nbnd, lv, rp, iqof, iqfo, 
     &                  ip1, ip, iq, ldv , xyz, labelabc)
        
        kefemabc=ke !kefemabc is the number of edges on cavity surface



        nbnd=0
        

        call labelfemm( mvolele, mtl, ke, nbnd, lv, rp, iqof, iqfo, 
     &                  ip1, ip, iq, ldv , xyz, labelsrc)

        kefemsrc=ke-kefemabc  !kefemsrc is the number of edges on port surface include PEC boudary

        kefemm=ke !kefemm is the number of edges on cavity surface and port surface include PEC boudary

        call labelfemo( mvolele, mtl, ke, lv, rp, iqof, iqfo, 
     &                  ip, ip1 )
        kefemo=ke-kefemm+nbnd !kefemo is the number fo edges on PEC surface 

        call labelfemi( mvolele, mtl, ke, lv, rp, iqof, iqfo, 
     &                  ip )

       ks = ke - kefemm
       keis= kefemo-nbnd
       kei = kefemo
!       print*,'ks keis and kei',ks,keis,kei
       do 8 i=1,kl
       do 8 j=1,6
           k = 6*(i-1)+j
           if (ip(k).le.kefemm) then
              le(j,i)=ip(k)+ks
           else
              le(j,i)=ip(k)-kefemm
           endif
8      continue

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
!           do k=1,3
!              ni(k)=mod(j+k-1,4)+1
!           enddo
        if (j.eq.1) then
           ni(1)=1
           ni(2)=3
           ni(3)=2
        else if (j.eq.2) then
           ni(1)=1
           ni(2)=2
           ni(3)=4
        else if (j.eq.3) then
           ni(1)=1
           ni(2)=4
           ni(3)=3
        else if (j.eq.4) then
           ni(1)=2
           ni(2)=3
           ni(3)=4
        end if
        lnc(4,klc)=i
 
ccccccccccccccctest normcccccccccccccccccccccccccc
          nx1=xyz(1,lv(ni(1),i))-xyz(1,lv(ni(2),i))
          ny1=xyz(2,lv(ni(1),i))-xyz(2,lv(ni(2),i))
          nz1=xyz(3,lv(ni(1),i))-xyz(3,lv(ni(2),i))
          nx2=xyz(1,lv(ni(3),i))-xyz(1,lv(ni(2),i))
          ny2=xyz(2,lv(ni(3),i))-xyz(2,lv(ni(2),i))
          nz2=xyz(3,lv(ni(3),i))-xyz(3,lv(ni(2),i))
          nx=ny2*nz1-nz2*ny1
          ny=nz2*nx1-nx2*nz1
          nz=nx2*ny1-ny2*nx1
          if(j.eq.1)then
            kk=lv(4,i)
          else if(j.eq.2)then
            kk=lv(3,i)
          else if(j.eq.3)then
            kk=lv(2,i)
          else if(j.eq.4)then
            kk=lv(1,i)
          end if
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
          if(test.le.0)then
          print*,"norm error in lv:",j,i
!          tt=ni(2)
!          ni(2)=ni(3)
!          ni(3)=tt
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
76    continue

      return
      end

      subroutine area(xyz1,xyz2,xyz3,det)

      real xyz1(3),xyz2(3),xyz3(3)

      det=xyz1(1)*xyz2(2)*xyz3(3)+xyz2(1)*xyz3(2)*xyz1(3)+
     &    xyz1(2)*xyz2(3)*xyz3(1)-xyz3(1)*xyz2(2)*xyz1(3)-
     &    xyz2(1)*xyz1(2)*xyz3(3)-xyz3(2)*xyz2(3)*xyz1(1)

      return
      end
