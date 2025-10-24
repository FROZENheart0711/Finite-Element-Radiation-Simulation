      Implicit none
      integer maxnode,patch(10),maxpatch,tetra(10),maxtetra,maxelement
      integer nsur,npec,nmom,nfem,nsurfem,iFlagpec,iTrpec
      
      character*256 Path,Dataname
      integer Uinp,Uinp_para,Uout

      character*256 tmp_filename1,tmp_filename2
      integer tmp1,tmp2,tmp3,iNodStr,nnodeq,iPthStr,nspatchq(10)
      integer iTetStr,nvtetraq(10) 

      integer i,j,k,ii,jj,kk,i1,i2,iFlag,iNum(2)
      integer m,n
      real    temp1,temp2,temp3,temp4

      character*256 tmp_char,tmp_char1,tmp_char2,tmp_char3
      integer nargs, iargc

      integer nt,j1,j2,j3
      integer numi1,numi2,numi3,numj1,numj2,numj3
      real    r1,r2,r3,sf,xt,yt,zt
      real    ri1,ri2,ri3,rj1,rj2,rj3
      integer facet
      integer ie1,ie2,iflag1,iflag2
      integer x1,x2,y1,y2,z1,z2
      real,allocatable::    xyz_n(:,:)
      integer,allocatable:: lva(:,:)
      real,allocatable::    rp(:)
      integer,allocatable:: iqof(:),iqfo(:),ip(:),iq(:)
      integer,allocatable:: fface(:,:),arr1_tmp(:)
      integer,allocatable:: ipat(:,:)
      integer,allocatable:: ele_all(:,:)
      real nx1,ny1,nz1,nx2,ny2,nz2,nx,ny,nz
      real norm(3),xc,yc,zc,test
      real xe,ye,ze
      real small
!c.....get commend line 
      Uinp=101
      Uout=102
      Uinp_para=103
      nargs = iargc()
      if(nargs.ne.1) then
         write(*,*)" USAGE: *.exe [filename]"
         write(*,*)" For example(Windows): "
         write(*,*)" a.exe sph0 "
         write(*,*)" For example(Linux): "
         write(*,*)" a.x sph0 "
         stop
      endif

      call getarg(1, tmp_char)
      tmp_char2=tmp_char
      Dataname=trim(tmp_char2);
      sf=0.001

      write(*,*) "transforming the *.dat to *.xyz,*.ipat,*.lv";
      write(*,*) "waiting....."

!      tmp_filename1=trim(Path)//trim(Dataname);
      tmp_filename1=trim(Dataname)
      tmp_filename2='.dat'
      tmp_filename1=trim(tmp_filename1)//trim(tmp_filename2)

      open( Uinp,file=tmp_filename1 )
      read(Uinp,*) tmp_char1
      nnodeq = 0
      iNodStr = 0           ! number of lines before nodcrd
      do while(trim(tmp_char1) .ne. 'GRID*') 
         read(Uinp,*) tmp_char1
         iNodStr = iNodStr+1
      enddo
      nnodeq = iNodStr   
       
      iFlag = 0
      tmp1 = 1
      tmp2 = 0
      do while(iFlag .ne. 1)
         read(Uinp,*) tmp_char1
         if(trim(tmp_char1) .eq. '*') then
            tmp2 = tmp2+1
            read(Uinp,*) tmp_char1
            if(trim(tmp_char1) .eq. 'GRID*') then
                tmp1 = tmp1+1
            else
                 iFlag = 1
            endif
         endif
      enddo
      if (tmp1 .eq. tmp2 ) then
        maxnode=tmp1
      else
        write(*,*) "Maxnode is wrong after seek!"
        stop
      endif
       
      do i=1,2
         iNum(i)= 0
      enddo
      do i=1,10
         nspatchq(i) = 0
         patch(i) = 0
         nvtetraq(i) = 0
         tetra(i) = 0
      enddo
      maxpatch = 0
      maxtetra = 0

      tmp_filename1=trim(Dataname)
      tmp_filename2='.para_inp'
      tmp_filename1=trim(tmp_filename1)//trim(tmp_filename2)      
      open( Uinp_para,file=tmp_filename1 )
      read(Uinp_para,*) iFlagpec
      iTrpec = 0
      if(iFlagpec .ne. 0) then
           read(Uinp_para,*) iTrpec
      endif
      
      iPthStr = 0            !number of lines before patch
      do while(trim(tmp_char1) .ne. 'CTRIA3')
         read(Uinp,*) tmp_char1
         iPthStr = iPthStr+1
      enddo
      iNum(1) = iNum(1) + 1
      nspatchq(iNum(1)) = iPthStr

      iFlag = 0
      tmp1 = 1
      do while(iFlag .ne. 1)
         read(Uinp,*) tmp_char1
         if(trim(tmp_char1) .eq. 'CTRIA3') then
            tmp1 = tmp1+1
         else
            iFlag = 1
         endif
      enddo
      if (tmp1 .gt. 0 ) then
        patch(iNum(1)) =  tmp1
      else
        write(*,*) "Maxpatch(die1) is wrong after seek!"
        stop
      endif      
      
800   tmp1 = 0
      do while(iand((trim(tmp_char1).ne.'CTRIA3'),
     &               (trim(tmp_char1).ne.'CTETRA')))
             read(Uinp,*) tmp_char1
             tmp1 = tmp1+1
      enddo
      if(trim(tmp_char1) .eq. 'CTRIA3') then
         iNum(1) = iNum(1) + 1
         nspatchq(iNum(1)) = tmp1
      elseif(trim(tmp_char1).eq.'CTETRA') then
         iNum(2) = iNum(2) + 1
         nvtetraq(iNum(2)) = tmp1
      else  
         write(*,*) "Element is wrong in seeking!"
         stop
      endif
       
      iTetStr = 0           !number of lines before tetra
      iFlag = 0
      tmp1 = 1
      if (iNum(2) .gt. 0) then
          do while(iFlag .ne. 1)
             read(Uinp,*) tmp_char1
             if(trim(tmp_char1) .eq. 'CTETRA') then
                 tmp1 = tmp1+1
             else
                 iFlag = 1
             endif
          enddo
          if (tmp1 .gt. 0 ) then
              tetra(iNum(2)) =  tmp1
          else
              write(*,*) "Maxtetra(1) is wrong after seek!"
              stop
          endif
          goto 400
      else 
          do while(iFlag .ne. 1)
             read(Uinp,*) tmp_char1
             if(trim(tmp_char1) .eq. 'CTRIA3') then
                 tmp1 = tmp1+1
             else
                  iFlag = 1
             endif
          enddo
          if (tmp1 .gt. 0 ) then
             patch(iNum(1)) =  tmp1
          else
             write(*,*) "Maxpatch is wrong with number :",iNum(1)
             stop
          endif
          goto 800
      endif

400   if (iNum(1) .lt. iTrpec) then
        write(*,*) "Patch kinds is wrong: iNum is less than iTrcpec!"
        stop
      endif
      do i=1,iNum(1)
         maxpatch = maxpatch + patch(i)
      enddo      

600   read(Uinp,*) tmp_char1
      if(trim(tmp_char1) .ne. 'ENDDATA') then
         iTetStr = 1
      do while(trim(tmp_char1) .ne. 'CTETRA')
         read(Uinp,*) tmp_char1
         iTetStr = iTetStr+1
      enddo
      iNum(2) = iNum(2) + 1      
      nvtetraq(iNum(2)) = iTetStr

      iFlag = 0
      tmp1 = 1
      do while(iFlag .ne. 1)
         read(Uinp,*) tmp_char1
         if(trim(tmp_char1) .eq. 'CTETRA') then
            tmp1 = tmp1+1
         else
            iFlag = 1
         endif
      enddo
      if (tmp1 .gt. 0 ) then
       tetra(iNum(2)) = tmp1
      else
        write(*,*) "Maxtetra is wrong after seek!"
        stop
      endif
         goto 600
      endif

      do i=1,iNum(2)
         maxtetra = maxtetra + tetra(i)
      enddo

      write(*,*),"********************************"
      write(*,*),"Maxnode by reading =",maxnode
      write(*,*),"NUMBER:  patch type=",iNum(1)
      if(iFlagpec .ne. 0) then
         write(*,*),"NUMBER:  pec1 patch=",iTrpec
      endif
      do i=1,iNum(1)
         write(*,*),"type and max:patch =",i,nspatchq(i),patch(i)
      enddo
      write(*,*),"Maxpatch by reading=",maxpatch
      write(*,*),"NUMBER:  tetra type=",iNum(2)
      do i=1,iNum(2)
        write(*,*),"type and max:tetra =",i+iNum(1),nvtetraq(i),tetra(i)
      enddo
      write(*,*),"Maxtetra by reading=",maxtetra
      rewind(Uinp) 
      maxelement = maxpatch + maxtetra       
       
      !============================================
      !  Read  xyz of nodes
      !============================================
      allocate(xyz_n(3,maxnode))
      do ii=1,nnodeq
        read(Uinp,*)
      enddo
      do ii=1,maxnode
        read(Uinp,1201) tmp_char1,tmp1,tmp2,xt,yt
        read(Uinp,1210) tmp_char2,tmp3,zt       
        if(tmp1 .eq. tmp3)  then
             xyz_n(1,ii)=xt
             xyz_n(2,ii)=yt
             xyz_n(3,ii)=zt 
        else
             write(*,*) "The read of Node Coordinate is Error!" 
             stop
        endif
      enddo
      !============================================
      !  Read nodes of element
      !============================================
      allocate(ele_all(5,maxelement))
      allocate(arr1_tmp(6) )
      npec=0
      nsur=0
      nfem=0
      do i=1,6
         arr1_tmp(i)=0
      end do
      kk = 0
      do jj=1,iNum(1)
      do ii=1,nspatchq(jj)
        read(Uinp,*)
      enddo
      do ii=1,patch(jj)
         kk = kk + 1 
        read(Uinp,1211) tmp_char1,arr1_tmp(1:5)
        tmp2=arr1_tmp(2)
        if(iFlagpec .ne. 0) then
          if(tmp2 .lt. iTrpec)then
             nsur=nsur+1 
             ele_all(5,kk)=1
!          else if(tmp2 .eq. iTrpec)then
          else
             npec=npec+1
             ele_all(5,kk)=2
!          else
!            write(*,*)  "The surpatch mesh is ERROR! PEC not one type!"
!            stop
          end if
        else
           nsur=nsur+1
           ele_all(5,kk)=1
        endif
        ele_all(1,kk)=arr1_tmp(3)
        ele_all(2,kk)=arr1_tmp(4)
        ele_all(3,kk)=arr1_tmp(5)
!        ele_all(5,kk)=tmp2
      end do
      end do
      if(nsur+npec .ne. maxpatch) then
         write(*,*)  "ERROR in reading CTRIA3!"
         stop
      endif
      if(kk .ne. maxpatch) then
         write(*,*)  "ERROR in writing CTRIA3!"
         stop
      endif
      
      do jj=1,iNum(2)
      do ii=1,nvtetraq(jj)
        read(Uinp,*)
      enddo
      do ii=1,tetra(jj)
        kk = kk + 1
        read(Uinp,1212) tmp_char1,arr1_tmp
        nfem=nfem+1
        ele_all(1,kk)=arr1_tmp(3)
        ele_all(2,kk)=arr1_tmp(4)
        ele_all(3,kk)=arr1_tmp(5)
        ele_all(4,kk)=arr1_tmp(6)
        ele_all(5,kk)=arr1_tmp(2)
      end do
      end do
      if(nfem .ne. maxtetra) then
         write(*,*)  "ERROR in reading CTETRA!"
         stop
      endif
      if(kk .ne. maxelement) then
         write(*,*)  "ERROR in writing CTETRA!"
         stop
      endif
      if(ele_all(5,maxpatch+1) .ge. iNum(1)+1) then
         if(ele_all(5,maxpatch+1) .gt. iNum(1)+1) then
         write(*,*) "Warning:there are blank seri-number before CTETRA!"
         endif
            kk=maxpatch
            tmp1=ele_all(5,maxpatch+1)-1 
         do ii=1,maxtetra
            kk=kk+1
            ele_all(5,kk)=ele_all(5,kk)-tmp1
         enddo
      else
         write(*,*)  "ERROR!"
         stop
      endif
      write(*,*),"********************************"
      write(*,*),"triangle on mat BI =",nsur
      write(*,*),"triangle on pec BI =",npec
      write(*,*),"tetrahedron in FEM =",nfem
      close(Uinp)  
      deallocate( arr1_tmp )
      nmom=nsur+npec
      !===================================================
      !   find outside surface patch of FEM
      !===================================================
      allocate(rp(4*nfem),iqof(4*nfem),iqfo(4*nfem))
      allocate(ip(4*nfem),iq(4*nfem))
      allocate(fface(4,4*nfem))
      nt=0
      do i=1,nfem
         ii=i+nmom
         do j=1,4
            if (j.eq.1) then
               j1 = ele_all(1,ii)
               j2 = ele_all(3,ii)
               j3 = ele_all(2,ii)
            else if (j.eq.2) then
               j1 = ele_all(1,ii)
               j2 = ele_all(2,ii)
               j3 = ele_all(4,ii)
            else if (j.eq.3) then
               j1 = ele_all(1,ii)
               j2 = ele_all(4,ii)
               j3 = ele_all(3,ii)
            else if(j.eq.4) then
               j1 = ele_all(2,ii)
               j2 = ele_all(3,ii)
               j3 = ele_all(4,ii)
            end if
            nt=nt+1
            fface(1,nt)=j1
            fface(2,nt)=j2
            fface(3,nt)=j3
            fface(4,nt)=100
            if(j.eq.1)then
               k=4
            else if(j.eq.2)then
               k=3
            else if(j.eq.3)then
               k=2
            else if(j.eq.4)then
               k=1
            end if
            nx1=xyz_n(1,j1)-xyz_n(1,j2)
            ny1=xyz_n(2,j1)-xyz_n(2,j2)
            nz1=xyz_n(3,j1)-xyz_n(3,j2)
            
            nx2=xyz_n(1,j3)-xyz_n(1,j2)
            ny2=xyz_n(2,j3)-xyz_n(2,j2)
            nz2=xyz_n(3,j3)-xyz_n(3,j2)

            nx=ny2*nz1-nz2*ny1
            ny=nz2*nx1-nx2*nz1
            nz=nx2*ny1-ny2*nx1
            
            xc=(xyz_n(1,j1)+xyz_n(1,j2)+xyz_n(1,j3))/3.0
            yc=(xyz_n(2,j1)+xyz_n(2,j2)+xyz_n(2,j3))/3.0
            zc=(xyz_n(3,j1)+xyz_n(3,j2)+xyz_n(3,j3))/3.0

            norm(1)=xc-xyz_n(1,ele_all(k,ii))
            norm(2)=yc-xyz_n(2,ele_all(k,ii))
            norm(3)=zc-xyz_n(3,ele_all(k,ii))
            test=nx*norm(1)+ny*norm(2)+nz*norm(3)
            if(test.le.0)then
              print*,"error in norm of lv:",i,j
!              print*,nx,ny,nz,norm(1),norm(2),norm(3)
!              print*,j1,xyz_n(1,j1),xyz_n(2,j1),xyz_n(3,j1)
!              print*,j2,xyz_n(1,j2),xyz_n(2,j2),xyz_n(3,j2)
!              print*,j3,xyz_n(1,j3),xyz_n(2,j3),xyz_n(3,j3)
              pause
            end if
         end do
      end do
        
      facet=nt
      i=0
      ii=0
      do jj=1,patch(4)
         j=jj+patch(1)+patch(2)+patch(3)
            j1 = ele_all(1,j)
            j2 = ele_all(2,j)
            j3 = ele_all(3,j)     
   
            nx1=xyz_n(1,j1)-xyz_n(1,j2)
            ny1=xyz_n(2,j1)-xyz_n(2,j2)
            nz1=xyz_n(3,j1)-xyz_n(3,j2)

            nx2=xyz_n(1,j3)-xyz_n(1,j2)
            ny2=xyz_n(2,j3)-xyz_n(2,j2)
            nz2=xyz_n(3,j3)-xyz_n(3,j2)

            nx=ny2*nz1-nz2*ny1
            ny=nz2*nx1-nx2*nz1
            nz=nx2*ny1-ny2*nx1

            xc=(xyz_n(1,j1)+xyz_n(1,j2)+xyz_n(1,j3))/3.0
            yc=(xyz_n(2,j1)+xyz_n(2,j2)+xyz_n(2,j3))/3.0
            zc=(xyz_n(3,j1)+xyz_n(3,j2)+xyz_n(3,j3))/3.0

            norm(1)=xc
            norm(2)=yc
            norm(3)=zc+4060
            test=nx*norm(1)+ny*norm(2)+nz*norm(3)
            if(test.le.0)then
!               print*,"error in norm of patch:",j
!              print*,nx,ny,nz,norm(1),norm(2),norm(3)
!              print*,j1,xyz_n(1,j1),xyz_n(2,j1),xyz_n(3,j1)
!              print*,j2,xyz_n(1,j2),xyz_n(2,j2),xyz_n(3,j2)
!              print*,j3,xyz_n(1,j3),xyz_n(2,j3),xyz_n(3,j3)
!               pause
               i = i+1
            else
               ii = ii +1
            end if
       enddo
       write(*,*) "the norm of patch is in-vector",i
       write(*,*) "the norm of patch is out-vector",ii

      do i=1,facet
         j1 = fface(1,i)
         j2 = fface(2,i)
         j3 = fface(3,i)
           
         call comi(j1,j2)
         call comi(j1,j3)
         call comi(j2,j3)

         r1=j1
         r2=j2
         r3=j3
         rp(i)=r1*log(r2)*log10(r3)
      end do
      call hpsort(facet,rp,iqof,iqfo)

      nt = 1
      ip(nt) = iqfo(1)
      iq(iqfo(1)) = nt

      do i=2,facet
         j=i-1
104      if (abs(rp(j)-rp(i)).gt.0.01) goto 204
         numi1 = fface(1,iqfo(i))
         numi2 = fface(2,iqfo(i))
         numi3 = fface(3,iqfo(i))
         call comi(numi1,numi2)
         call comi(numi1,numi3)
         call comi(numi2,numi3)

         numj1 = fface(1,iqfo(j))
         numj2 = fface(2,iqfo(j))
         numj3 = fface(3,iqfo(j))
         call comi(numj1,numj2)
         call comi(numj1,numj3)
         call comi(numj2,numj3)

         if((numi1.eq.numj1).and.(numi2.eq.numj2)
     &     .and.(numi3.eq.numj3))then
           ip(iqfo(i))=ip(iqfo(j))
           fface(4,iqfo(i)) = 0
           fface(4,iqfo(j)) = 0
           goto 404
         else
            j=j-1
            if (j.le.0) goto 204
          goto 104
        endif

 204    nt=nt+1
        ip(nt)=iqfo(i)
        iq(iqfo(i))=nt
 404  enddo
      nt=0
      do i=1,facet
         if(fface(4,i).eq.100)then
            nt=nt+1
         end if
      end do
      nsurfem=nt
      write(*,*),"********************************"
      write(*,*),"surpatch of FEM    =",nsurfem
      allocate(ipat(5,nsurfem))
      allocate(lva(4,nfem))
      do i=1,nfem
      do j=1,4
         lva(j,i)=0
      end do
      end do
 
      nt=0
      do i=1,facet
         if(fface(4,i).ne.0)then
            nt=nt+1
            ipat(1,nt)=fface(1,i)
            ipat(2,nt)=fface(2,i)
            ipat(3,nt)=fface(3,i)
            j=(i-1)/4+1
            ipat(4,nt)=j
            k=i-((i-1)/4)*4
            ipat(5,nt)=k
         end if
      end do
      deallocate(fface)
      deallocate(rp,iqof,iqfo,ip,iq)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      allocate(rp(nsurfem+nsur),iqof(nsurfem+nsur),iqfo(nsurfem+nsur))
      allocate(ip(nsurfem+nsur),iq(nsurfem+nsur))
      allocate(fface(4,nsurfem+nsur))
      nt=0
      do i=1,nsurfem
         nt=nt+1
         fface(1,nt)=ipat(1,i)
         fface(2,nt)=ipat(2,i)
         fface(3,nt)=ipat(3,i)
         fface(4,nt)=100
      end do
      do i=1,nsur
         nt=nt+1
         fface(1,nt)=ele_all(1,i)
         fface(2,nt)=ele_all(2,i)
         fface(3,nt)=ele_all(3,i)
         fface(4,nt)=100
      end do 
      facet=nt
      do i=1,facet
         j1 = fface(1,i)
         j2 = fface(2,i)
         j3 = fface(3,i)
           
         call comi(j1,j2)
         call comi(j1,j3)
         call comi(j2,j3)

         r1=j1
         r2=j2
         r3=j3
         rp(i)=r1*log(r2)*log10(r3)
      end do
      call hpsort(facet,rp,iqof,iqfo)
      nt = 1
      ip(nt) = iqfo(1)
      iq(iqfo(1)) = nt

      do i=2,facet
         j=i-1
105      if (abs(rp(j)-rp(i)).gt.0.01) goto 205
         numi1 = fface(1,iqfo(i))
         numi2 = fface(2,iqfo(i))
         numi3 = fface(3,iqfo(i))
         call comi(numi1,numi2)
         call comi(numi1,numi3)
         call comi(numi2,numi3)

         numj1 = fface(1,iqfo(j))
         numj2 = fface(2,iqfo(j))
         numj3 = fface(3,iqfo(j))
         call comi(numj1,numj2)
         call comi(numj1,numj3)
         call comi(numj2,numj3)

         if((numi1.eq.numj1).and.(numi2.eq.numj2)
     &     .and.(numi3.eq.numj3))then
           ip(iqfo(i))=ip(iqfo(j))
           fface(4,iqfo(i)) = 0
           fface(4,iqfo(j)) = 0
           goto 405
         else
            j=j-1
            if (j.le.0) goto 205
          goto 105
        endif

 205    nt=nt+1
        ip(nt)=iqfo(i)
        iq(iqfo(i))=nt
 405  enddo
      nt=0
      do i=1,nsurfem
         if(fface(4,i).eq.0)then
           nt=nt+1
           lva(ipat(5,i),ipat(4,i))=1
         end if
      end do
      if(nt.ne.nsur)then
        write(*,*),"nt .ne. nsur:",nt,nsur
!        pause
      end if
      deallocate(fface)
      deallocate(rp,iqof,iqfo,ip,iq)
      write(*,*),"********************************"
!cccc   Record the *.xyz
!      tmp_filename1=trim(Path)//trim(Dataname);
      tmp_filename1=trim(Dataname)
      tmp_filename2='.xyz'
      tmp_filename1=trim(tmp_filename1)//trim(tmp_filename2)
      open(Uout,file=tmp_filename1)
      write(Uout,*) maxnode
      do i=1,maxnode
         do j=1,3
            xyz_n(j,i)=xyz_n(j,i)*sf
         end do         
         write(Uout,1202) (xyz_n(j,i),j=1,3)
      end do
      close(Uout)
      write(*,*) "finish creating ",trim(tmp_filename1) 
!cccc   Record the *.ipat
!      tmp_filename1=trim(Path)//trim(Dataname);
      tmp_filename1=trim(Dataname)
      tmp_filename2='.ipat'
      tmp_filename1=trim(tmp_filename1)//trim(tmp_filename2)
      open(Uout,file=tmp_filename1)
      write(Uout,1204) nmom
      do i=1,nmom
         write(Uout,1203) (ele_all(j,i),j=1,3),ele_all(5,i)
      end do
      close(Uout)
      write(*,*) "finish creating ",trim(tmp_filename1)
!cccc  Record the *.lv     
!      tmp_filename1=trim(Path)//trim(Dataname);
      tmp_filename1=trim(Dataname)
      tmp_filename2='.lv'
      tmp_filename1=trim(tmp_filename1)//trim(tmp_filename2)
      open(Uout,file=tmp_filename1)
      write(Uout,1205) nfem
      do i=1,nfem
         write(Uout,1205) (ele_all(j,i+nmom),j=1,4),
     &            ele_all(5,i+nmom),(lva(j,i),j=1,4)
      end do
      close(Uout)
      write(*,*) "finish creating ",trim(tmp_filename1)
1201  format (1A16,1I8,3E16.9)
1210  format (1A1,1I7,1E16.9)
1211  format (1A8,5I8)
1212  format (1A8,6I8)
1202  format (3E20.7)
1203  format (19I8)
1204  format (3I16)
1205  format (9I12)
1206  format (6I12)
1207  format (1I12)
1208  format (4I12)
      end
