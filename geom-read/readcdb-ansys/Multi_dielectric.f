      Implicit none
      integer maxnode,maxpatch
      integer nsur,npec,nmom,nfem,nsurfem
      
      character*256 Path,Dataname
      integer Uinp,Uout

      character*256 tmp_filename1,tmp_filename2
      integer tmp1,tmp2,tmp3

      integer i,j,k,ii,jj,kk,i1,i2
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
      nargs = iargc()
      if(nargs.ne.1) then
         write(*,*)" USAGE: *.exe [filename]"
         write(*,*)" For example(Windows): "
         write(*,*)" a.exe sph0 "
         write(*,*)" For example(Linux): "
         write(*,*)" a.x sph0 "
         stop
      endif

!      tmp_char='./ttt'
!      inquire(directory=trim(tmp_char),exist=aaa)
!      print*,aaa
!      pause
!      tmp_char='mkdir ttt'
!      call system(tmp_char) 
  
      call getarg(1, tmp_char)
      tmp_char2=tmp_char
!      call getarg(2, tmp_char)
!      tmp_char3=tmp_char
         
!      Path=trim(tmp_char3);
      Dataname=trim(tmp_char2);

      write(*,*) "transforming the *.cdb to *.xyz,*.ipat,*.lv";
      write(*,*) "waiting....."

!      tmp_filename1=trim(Path)//trim(Dataname);
      tmp_filename1=trim(Dataname)
      tmp_filename2='.cdb'
      tmp_filename1=trim(tmp_filename1)//trim(tmp_filename2)

      open( Uinp,file=tmp_filename1 )
      do ii=1,5
        read(Uinp,*) 
      enddo
      read(Uinp,*) tmp_char1,tmp_char2,tmp1
   
      if (tmp1>0 ) then
        maxnode=tmp1
      else
        write(*,*) "Maxnode is wrong at line 6."
        stop
      endif

      read(Uinp,*) tmp_char1,tmp_char2,tmp1
   
      if (tmp1>0 ) then
         maxpatch=tmp1
      else
         write(*,*) "Maxpatch is wrong at line 7."
         stop
      endif
               
      do ii=1,2
        read(Uinp,*) 
      enddo

      read(Uinp,*) tmp_char1,tmp_char2

      tmp_char3=trim(tmp_char1)//','//trim(tmp_char2);
      if (trim(tmp_char3)/='(3i8,6e20.13)') then
      write(*,*) " The file's format dosen't conform..."
      write(*,*) " The program's format is (3i8,6e16.9), but
     &     the format of the *.cdb (line 10) is: ",trim(tmp_char3)
      write(*,*) " Please check the *.cdb again or 
     &             change the format of this transform program "
      stop
      endif
               
      !============================================
      !  Read  xyz of nodes
      !============================================
      allocate(xyz_n(3,maxnode))
      do ii=1,maxnode
        read(Uinp,1201) tmp1,tmp2,tmp3,xt,yt,zt
        xyz_n(1,ii)=xt
        xyz_n(2,ii)=yt
        xyz_n(3,ii)=zt 
      enddo
      !============================================
      !  Read nodes of element
      !============================================
      read(Uinp,*)
      read(Uinp,*)
      read(Uinp,*) tmp_char1

      if (trim(tmp_char1)/='(19i8)') then
        write(*,*) " The file's format dosen't conform..."
        write(*,*) " The program's format is (19i8), but 
     &     the format of the *.cdb (line 849) is: ",trim(tmp_char1)
        write(*,*) " Please check the *.cdb again or change the format
     &               of this transform program "
        stop
      endif
      allocate(ele_all(5,maxpatch))
      allocate(arr1_tmp(19) )
      do i=1,19
         arr1_tmp(i)=0
      end do
      npec=0
      nsur=0
      nfem=0
      do ii=1,maxpatch
        read(Uinp,1203) arr1_tmp
        tmp2=arr1_tmp(2)-3
        if(tmp2.eq.1)then
          nsur=nsur+1
        else if(tmp2.eq.2)then
          npec=npec+1
        else
          nfem=nfem+1
        end if
        ele_all(1,ii)=arr1_tmp(12)
        ele_all(2,ii)=arr1_tmp(13)
        ele_all(3,ii)=arr1_tmp(14)
        ele_all(4,ii)=arr1_tmp(15)
        ele_all(5,ii)=tmp2
      end do
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
!              print*,xyz_n(1,j1),xyz_n(2,j1),xyz_n(3,j1)
!              print*,xyz_n(1,j2),xyz_n(2,j2),xyz_n(3,j2)
!              print*,xyz_n(1,j3),xyz_n(2,j3),xyz_n(3,j3)
              pause
            end if
         end do
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
      write(*,*),"surpatch of FEM=",nsurfem
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
        pause
      end if
      deallocate(fface)
      deallocate(rp,iqof,iqfo,ip,iq)
!cccc   Record the *.xyz
!      tmp_filename1=trim(Path)//trim(Dataname);
      tmp_filename1=trim(Dataname)
      tmp_filename2='.xyz'
      tmp_filename1=trim(tmp_filename1)//trim(tmp_filename2)
      open(Uout,file=tmp_filename1)
      write(Uout,*) maxnode
      do i=1,maxnode
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
     &            ele_all(5,i+nmom)+3,(lva(j,i),j=1,4)
      end do
      close(Uout)
      write(*,*) "finish creating ",trim(tmp_filename1)
1201  format (3I8,6E20.13)
1202  format (3E20.7)
1203  format (19I8)
1204  format (3I16)
1205  format (9I12)
1206  format (6I12)
1207  format (1I12)
1208  format (4I12)
      end
