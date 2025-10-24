      subroutine findpecedge(klv,mvoledge,npecpat,lv,
     &           le,pecpat,pecjudge,mpecedge,xyzv,
     &           jedgen,jedgeco,jejudge)
      implicit none
      integer klv,npecpat,mpecedge,mvoledge,jedgen
      integer lv(9,klv),le(6,klv),pecpat(3,npecpat)
      integer jedgeco(2,jedgen),jejudge(jedgen)
      real    xyzv(3,*)
      integer nt,i,j,k,k1,k2
      integer j1,j2,numi1,numi2,numj1,numj2
      real    r1,r2
      integer pecjudge(*)
      integer,allocatable::edgecot(:,:),pecedgeco(:,:)
      integer,allocatable::edgeco(:,:)
      real,allocatable::rp(:)
      integer,allocatable::iqof(:),iqfo(:),ip(:),iq(:)
      real,allocatable::middlenode(:,:)
      !==========================================
      !  Get the edges on PEC
      !==========================================
      print*,klv,mvoledge,npecpat,mpecedge
      allocate(edgecot(2,6*klv),edgeco(2,mvoledge))
      allocate(pecedgeco(2,3*npecpat))
      allocate(rp(6*klv),iqof(6*klv),iqfo(6*klv),
     &         ip(6*klv),iq(6*klv))
      allocate(middlenode(3,mvoledge))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,klv
         do j=1,6
           call labelinv(j,k1,k2)
           k1=lv(k1,i)
           k2=lv(k2,i)
           call comi(k1,k2)
           edgeco(1,le(j,i))=k1
           edgeco(2,le(j,i))=k2
           if(le(j,i).gt.mvoledge)then
             print*,"error in le:",i,j
           end if
         end do
      end do
      do i=1,mvoledge
         do j=1,3
         middlenode(j,i)=(xyzv(j,edgeco(1,i))+xyzv(j,edgeco(2,i)))/2.0
         end do
      end do 
      nt=0
      do i=1,npecpat
         do j=1,3
           nt=nt+1
           if(j.eq.1)then
              j1=1
              j2=2
           else if(j.eq.2)then
              j1=2
              j2=3
           else if(j.eq.3)then
              j1=1
              j2=3
           end if     
           j1=pecpat(j1,i)
           j2=pecpat(j2,i)
           call comi(j1,j2)
           edgecot(1,nt)=j1
           edgecot(2,nt)=j2
           r1=j1
           r2=j2
           rp(nt)=r1*sqrt(r2)
         end do
       end do
       call hpsort(nt,rp,iqof,iqfo)
       nt=1
       ip(iqfo(1))=1
       iq(1)=iqfo(1)

       do i=2,3*npecpat
          j=i-1
101       if(abs(rp(j)-rp(i)).gt.0.01) goto 201

          numi1=edgecot(1,iqfo(i))
          numi2=edgecot(2,iqfo(i))

          numj1=edgecot(1,iqfo(j))
          numj2=edgecot(2,iqfo(j))
          
          if((numi1.eq.numj1).and.(numi2.eq.numj2)) then
             ip(iqfo(i))=ip(iqfo(j))
             goto 401
          else
             j=j-1
             if(j.le.0) goto 201
             goto 101
          end if
201       nt=nt+1
          ip(iqfo(i))=nt
          iq(nt)=iqfo(i)
401   end do
      mpecedge=nt
      print*,"edge on PEC:",mpecedge 
      do i=1,mpecedge
         pecedgeco(1,i)=edgecot(1,iq(i))
         pecedgeco(2,i)=edgecot(2,iq(i))
      end do
 
      nt=0
      do i=1,mpecedge
         nt=nt+1
         edgecot(1,nt)=pecedgeco(1,i)
         edgecot(2,nt)=pecedgeco(2,i)
         r1=pecedgeco(1,i)
         r2=pecedgeco(2,i)
         rp(nt)=r1*sqrt(r2)
      end do
      do i=1,mvoledge
         nt=nt+1
         edgecot(1,nt)=edgeco(1,i)
         edgecot(2,nt)=edgeco(2,i)
         r1=edgeco(1,i)
         r2=edgeco(2,i)
         rp(nt)=r1*sqrt(r2)
      end do
      do i=1,nt
         ip(i)=-100
      end do
      call hpsort(nt,rp,iqof,iqfo)
      do i=2,nt
         j=i-1
103      if(abs(rp(j)-rp(i)).gt.0.01) goto 203
         numi1=edgecot(1,iqfo(i))
         numi2=edgecot(2,iqfo(i))

         numj1=edgecot(1,iqfo(j))
         numj2=edgecot(2,iqfo(j))
         if((numi1.eq.numj1).and.(numi2.eq.numj2)) then
            ip(iqfo(i))=iqfo(j)
            ip(iqfo(j))=iqfo(i)
            goto 203
         else 
           j=j-1
            if(j.le.0) goto 203
            goto 103
         end if
203   end do
      
      do i=1,mpecedge
         if(ip(i).le.0)then
          write(*,*),"error:",i,ip(i)
          if(ip(i).le.mpecedge)then
            write(*,*),"error:",ip(i),mpecedge
          end if
          pause
         end if
         pecjudge(i)=ip(i)-mpecedge
!         print*,i,pecjudge(i),(middlenode(j,pecjudge(i)),j=1,3)
      enddo

      do i=1,jedgen
         j1=jedgeco(1,i)
         j2=jedgeco(2,i)
         call comi(j1,j2)
         jedgeco(1,i)=j1
         jedgeco(2,i)=j2
      end do

      nt=0
      do i=1,jedgen
         nt=nt+1
         edgecot(1,nt)=jedgeco(1,i)
         edgecot(2,nt)=jedgeco(2,i)
         r1=jedgeco(1,i)
         r2=jedgeco(2,i)
         rp(nt)=r1*sqrt(r2)
      end do
      do i=1,mvoledge
         nt=nt+1
         edgecot(1,nt)=edgeco(1,i)
         edgecot(2,nt)=edgeco(2,i)
         r1=edgeco(1,i)
         r2=edgeco(2,i)
         rp(nt)=r1*sqrt(r2)
      end do
      do i=1,nt
         ip(i)=-100
      end do
      call hpsort(nt,rp,iqof,iqfo)
      do i=2,nt
         j=i-1
104      if(abs(rp(j)-rp(i)).gt.0.01) goto 204
         numi1=edgecot(1,iqfo(i))
         numi2=edgecot(2,iqfo(i))

         numj1=edgecot(1,iqfo(j))
         numj2=edgecot(2,iqfo(j))
         if((numi1.eq.numj1).and.(numi2.eq.numj2)) then
            ip(iqfo(i))=iqfo(j)
            ip(iqfo(j))=iqfo(i)
            goto 204
         else 
           j=j-1
            if(j.le.0) goto 204
            goto 104
         end if
204   end do
      
      do i=1,jedgen
         if(ip(i).le.0)then
          write(*,*),"error:",i,ip(i)
          if(ip(i).le.jedgen)then
            write(*,*),"error:",ip(i),jedgen
          end if
          pause
         end if
         jejudge(i)=ip(i)-jedgen
         print*,i,jejudge(i),(middlenode(j,jejudge(i)),j=1,3)
      enddo


      return
      end
 
