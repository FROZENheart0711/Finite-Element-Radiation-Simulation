      Implicit none
      integer maxnode,patch(20),maxpatch,tetra(20),maxtetra,maxelement
      integer nsur,npec,nmom,nfem,nsurfem,iFlagpec,iTrpec
      
      character*256 Path,Dataname
      integer Uinp,Uinp_para,Uout

      character*256 tmp_filename1,tmp_filename2
      integer tmp1,tmp2,tmp3,iNodStr,nnodeq,iPthStr,nspatchq(20)
      integer iTetStr,nvtetraq(20) 

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
      do i=1,20
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

      write(*,*),"********************************"
      write(*,*),"Maxnode by reading =",maxnode
      write(*,*),"**********Patch-Element*********"
      write(*,*),"NUMBER:  patch type=",iNum(1)
      if(iFlagpec .ne. 0) then
      write(*,*),"NUMBER:N1.pec patch=",iTrpec
      endif
      do i=1,iNum(1)
         write(*,*),"type and max:patch =",i,patch(i)
      enddo
      write(*,*),"Maxpatch by reading=",maxpatch
      rewind(Uinp) 
      maxelement = maxpatch
       
      !============================================
      !  Read  xyz of nodes
      !============================================
      allocate(xyz_n(3,maxnode))
      do ii=1,nnodeq
        read(Uinp,*)
      enddo
      do ii=1,maxnode
        read(Uinp,1201) tmp_char1,tmp1,tmp2,xt,yt
        if(tmp_filename1.eq.'ship-3000w.dat')then
          read(Uinp,1209) tmp_char2,tmp3,zt       
        else
          read(Uinp,1210) tmp_char2,zt
        endif
c        if(tmp1 .eq. tmp3)  then
        if(tmp1.eq.ii) then
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
      maxpatch=patch(1)
      allocate(ele_all(5,maxpatch))
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
      
      write(*,*),"********************************"
      write(*,*),"triangle on mat BI =",nsur
      write(*,*),"triangle on pec BI =",npec
      write(*,*),"tetrahedron in FEM =",nfem
      close(Uinp)  
      deallocate( arr1_tmp )
      nmom=nsur+npec
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
1201  format (1A16,1I8,3E16.9)
1209  format (1A1,1I7,1E16.9)
1210  format (1A8,1E16.9)
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
