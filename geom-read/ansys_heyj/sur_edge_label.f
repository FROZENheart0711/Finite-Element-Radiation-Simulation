        subroutine sur_edge_label
     & ( nnode, npatch, nedge,
     &   ncnode, ncele, ncedge, nbnd,  
     &   ipat, iedge, 
     &   ln, lp, lds, le,
     &   ip, iq, nc,
     &   xyz, xyzc,
     &   rp, iqof, iqfo, labelabc, labelsrc )

c.....nnode, npatch, nedge are the number of nodes...on the whole surface
c.....ncnode, ncele, ncedge are the number of nodes..on the cavity surface
c.....nbnd is the number of nodes on the boundary between cavity  and non-cavity
c.....ipat, iedge are for MoM
c.....ln, lp are information for cavity
c.....ip is the relation of edge label between local label on cavity and global
c.....xyzc is coordinate of nodes on cavity
c.....nc is the global label of element on cavity
   
        integer ipat(4, *),iedge(5, *)
        integer ln(3, *), le(3, *), lp(3, *),lds(2, *)
        integer ip(*), iq(*), nc(*)
        real    xyz(3, *), xyzc(3, *)
        integer iqof(*),iqfo(*)
        real    rp(*)

c.....define all edges and sort

        mtl  = 3*npatch
        ncele= 0
        nedge= 0
        nbnd = 0

        do i=1,mtl
           iq(i)=0
           nc(i)=0
        enddo

        do i=1,npatch
        do j=1,3
           call edgenod(j,k1,k2)
           kk=3*(i-1)+j
           rt1=ipat(k1,i)
           rt2=ipat(k2,i)
           if (rt1.gt.rt2) then
              rts=rt1
              rt1=rt2
              rt2=rts
           else
           endif
           rp(kk)=rt2*log(rt1+5)
        enddo
        enddo
        call hpsort(mtl,rp,iqof,iqfo)

c.....label edges on conducting body

      do 20 i=1,mtl


         ie1=(iqfo(i)-1)/3+1
         if ((ipat(4,ie1).ne.labelabc).and.
     &     (ipat(4,ie1).ne.labelsrc)) then
            testi=1
         else
            testi=0
         endif
      
         if (testi.eq.1) then
      
             i1=i-1

           if (i1.le.0) goto 40

30         if (rp(i).ne.rp(i1)) goto 40

             k=iqfo(i)-((iqfo(i)-1)/3)*3
             call edgenod(k,k1,k2)
             k1=ipat(k1,ie1)
             k2=ipat(k2,ie1)
             call com(k1,k2,iflag1)

             lo=iqfo(i1)-((iqfo(i1)-1)/3)*3
             ie2=(iqfo(i1)-1)/3+1
             call edgenod(lo,l1,l2)
             l1=ipat(l1,ie2)
             l2=ipat(l2,ie2)
             call com(l1,l2,iflag2)

             if ((k1.eq.l1).and.(k2.eq.l2)) then
               if ((ipat(4,ie2).ne.labelabc) .and.
     &              (ipat(4,ie2).ne.labelsrc)) then
                  testj=1
               else
                  testj=0
               endif

               if (testj.eq.1) then
                   iq(iqfo(i))=iq(iqfo(i1))
               else
                  nbnd  = nbnd+1
                  nedge = nedge+1
                  iq(iqfo(i)) = nedge
                  iq(iqfo(i1))= nedge
               endif

               if (iflag1.eq.1) then
                  iedge(1,iq(iqfo(i)))=k1
                  iedge(2,iq(iqfo(i)))=k2
                  iedge(3,iq(iqfo(i)))=ie1
                  iedge(4,iq(iqfo(i)))=ie2
               else
                  iedge(1,iq(iqfo(i)))=k1
                  iedge(2,iq(iqfo(i)))=k2
                  iedge(3,iq(iqfo(i)))=ie2
                  iedge(4,iq(iqfo(i)))=ie1
               endif

            else
               i1=i1-1
               goto 30
 40            nedge=nedge+1
               iq(iqfo(i))=nedge
            endif

      else

c.......case for element on cavity surface
        ie1=(iqfo(i)-1)/3+1
        do ii=1,ncele
           if (nc(ii).eq.ie1) goto 45
        enddo
        ncele=ncele+1
        nc(ncele)=ie1

45      i1=i-1
50      if (rp(i).ne.rp(i1)) goto 20
            k=iqfo(i)-((iqfo(i)-1)/3)*3
            call edgenod(k,k1,k2)
            k1=ipat(k1,ie1)
            k2=ipat(k2,ie1)
            call com(k1,k2,iflag1)

            lo=iqfo(i1)-((iqfo(i1)-1)/3)*3
            ie2=(iqfo(i1)-1)/3+1
            call edgenod(lo,l1,l2)
            l1=ipat(l1,ie2)
            l2=ipat(l2,ie2)
            call com(l1,l2,iflag2)

            if ((k1.eq.l1).and.(k2.eq.l2)) then

              if ((ipat(4,ie2).ne.labelabc).and.
     &            (ipat(4,ie2).ne.labelsrc)) then
                  testj=1
              else
                  testj=0
              endif
      
              if (testj.ne.1) then
                 i1=i1-1
                 goto 50
              else
                nbnd=nbnd+1
                iq(iqfo(i))=iq(iqfo(i1))
              endif

             if (iflag1.eq.1) then
                iedge(1,iq(iqfo(i)))=k1
                iedge(2,iq(iqfo(i)))=k2
                iedge(3,iq(iqfo(i)))=ie1
                iedge(4,iq(iqfo(i)))=ie2
              else
                iedge(1,iq(iqfo(i)))=k1
                iedge(2,iq(iqfo(i)))=k2
                iedge(3,iq(iqfo(i)))=ie2
                iedge(4,iq(iqfo(i)))=ie1
              endif
            else
              i1=i1-1
              goto 50
            endif
        endif
 20     continue
        nledge=nedge !nledge is the number of PEC edges include boundary
c.....define edges on cavity surface and sort

        mtl=3*ncele  

        do i=1,ncele
        do j=1,3
           ii=nc(i)
           call edgenod(j,k1,k2)
           kk=3*(i-1)+j
           rt1=ipat(k1,ii)
           rt2=ipat(k2,ii)
           if (rt1.gt.rt2) then
              rts=rt1
              rt1=rt2
              rt2=rts
           else
           endif
           rp(kk)=rt2*log(rt1+5)
        enddo
        enddo

        call hpsort(mtl,rp,iqof,iqfo)

c.....label edges on cavity surface

        do 80 i=1,mtl
           ki=int((iqfo(i)-1)/3)+1
           kj=iqfo(i)-3*(ki-1)
           kc=3*(nc(ki)-1)+kj

           if (iq(kc).ne.0) goto 80
           i1=i-1
90         if (rp(i).ne.rp(i1)) goto 100

            k=iqfo(i)-((iqfo(i)-1)/3)*3
            ie1=(iqfo(i)-1)/3+1
            ie1=nc(ie1)
            call edgenod(k,k1,k2)
            k1=ipat(k1,ie1)
            k2=ipat(k2,ie1)
            call com(k1,k2,iflag1)

            lo=iqfo(i1)-((iqfo(i1)-1)/3)*3
            ie2=(iqfo(i1)-1)/3+1
            ie2=nc(ie2)
            call edgenod(lo,l1,l2)
            l1=ipat(l1,ie2)
            l2=ipat(l2,ie2)
            call com(l1,l2,iflag2)

            if ((k1.eq.l1).and.(k2.eq.l2)) then
            kd=3*(ie2-1)+lo
            iq(kc)=iq(kd)
            if (iflag1.eq.1) then
               iedge(1,iq(kc))=k1
               iedge(2,iq(kc))=k2
               iedge(3,iq(kc))=ie1
               iedge(4,iq(kc))=ie2
            else
               iedge(1,iq(kc))=k1
               iedge(2,iq(kc))=k2
               iedge(3,iq(kc))=ie2
               iedge(4,iq(kc))=ie1
            endif
            else
              i1=i1-1
              goto 90
100           nedge=nedge+1
              iq(kc)=nedge
            endif
 80      continue
!        nedge is the number of all edges: pec+cavity
         do i = 1, npatch
         do j = 1, 3
            le(j,i) = iq(3*(i-1)+j)
         enddo
         enddo
c.....For FEM nodes and coordinate

      if (ncele.eq.0) then
         ncnode=0
         ncedge=0
         goto 200
      else
      endif

      mtl=3*ncele
      do i=1, ncele
      do j=1, 3
         ii=nc(i)
         kk=3*(i-1)+j
         rp(kk)=ipat(j,ii)
      enddo
      enddo

      call  hpsort(mtl,rp,iqof,iqfo)

      k=1
      iq(int(rp(k)))=k
      do i=1,3
         xyzc(i,k)=xyz(i,int(rp(k)))
      enddo

      do i=2,mtl
         i1=i-1
         if (rp(i).ne.rp(i1)) then
           k=k+1
           iq(int(rp(i)))=k 
           do j=1,3
              xyzc(j, k)=xyz(j, int(rp(i)))
           enddo
         else
           iq(int(rp(i)))=iq(int(rp(i1)))
         endif
      enddo
      ncnode=k
      do i=1,ncele
      do j=1,3
         ii=nc(i)
         ln(j,i)=iq(ipat(j,ii))
      enddo
      enddo

c.....For FEM edges
      do i=1,ncele
      do j=1,3
         ii=nc(i)
         kk=3*(i-1)+j
         rp(kk)=le(j,ii)
      enddo
      enddo

      call  hpsort(mtl,rp,iqof,iqfo)

      k=1
      iq(int(rp(k)))=k
      ip(k)=int(rp(k))

      do i=2,mtl
         i1=i-1
         if (rp(i).ne.rp(i1)) then
           k=k+1
           iq(int(rp(i)))=k
           ip(k)=int(rp(i))      
         else
           iq(int(rp(i)))=iq(int(rp(i1)))      
         endif
      enddo

      ncedge=k

      do i=1,ncele
      do j=1,3
         ii=nc(i)
         lp(j,i)=iq(le(j,ii))
      enddo
      enddo

      do i=1,ncedge
         k1=iedge(1,ip(i))
         k2=iedge(2,ip(i))
         call comi(k1,k2)
         lds(1,i)=k1
         lds(2,i)=k2
      enddo
 200  continue
      return
      end
