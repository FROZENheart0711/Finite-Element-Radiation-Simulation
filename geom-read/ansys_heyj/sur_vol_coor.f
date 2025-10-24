      subroutine sur_vol_coor
     &( mn, kls, klv,
     &  nsnode, nvnode,
     &  ln, lv,
     &  xyz, xyzs, xyzv,
     &  rp,iqof,iqfo,ip )
      implicit none
      real    xyz(3, *),xyzs(3, *),xyzv(3, *)
      integer ln(4, *), lv(9, *)
      real    rp(*)
      integer iqof(*), iqfo(*), ip(*)
      integer i,j,k,kls,mn,klv,nsnode,nvnode,mtl
      integer i1,i2
c....information for surface
      do i=1, kls
      do j=1, 3
         k=3*(i-1)+j
         rp(k)=ln(j,i)
      enddo
      enddo
      mtl = 3*kls
      call hpsort(mtl, rp, iqof, iqfo)

      nsnode = 1
      ip(1)  = nsnode
      do j=1,3
         xyzs(j,nsnode)=xyz(j,int(rp(1)))
      enddo


      do i=2, mtl
         i1=i-1
         if (rp(i).ne.rp(i1)) then
            nsnode = nsnode+1
            do j=1,3
               xyzs(j,nsnode)=xyz(j,int(rp(i)))
            enddo
         else
         endif
            ip(i)  = nsnode
      enddo

      do i=1, mtl
         i1=(iqfo(i)-1)/3+1
         i2=iqfo(i)-(i1-1)*3
         ln(i2,i1)=ip(i)
      enddo

c....Information for volume
      do i=1, klv
      do j=1, 4
         k=4*(i-1)+j
         rp(k)=lv(j, i)
      enddo
      enddo
      mtl = 4*klv
      call hpsort(mtl, rp, iqof, iqfo)
      nvnode = 1
      ip(1)  = nvnode
      do j=1,3
         xyzv(j,nvnode) = xyz(j,int(rp(1)))
      enddo

      do i=2, mtl
         if (rp(i).ne.rp(i-1)) then
            nvnode  = nvnode + 1
            do j=1,3
               xyzv(j,nvnode) = xyz(j,int(rp(i)))
            enddo
         else
         endif
            ip(i) = nvnode
      enddo

      do i=1, mtl
         i1=(iqfo(i)-1)/4+1
         i2=iqfo(i)-(i1-1)*4
         lv(i2,i1)=ip(i)
      enddo
      return
      end

      
