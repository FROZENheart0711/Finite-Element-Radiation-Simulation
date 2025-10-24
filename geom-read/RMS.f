      implicit none
      integer i,j
      real freq, theta, phi, rcst, rcsc, phase

      real rcs1(3601), rcs2(3601)    
      real rcsv, rcsh, phasev,phaseh  

      real RMS
      character*256 runname, tailname, scomb
      integer nargs, iargc
      integer nangle

      call getarg(1, runname)
 
      print*,'Please input number of data:'
      read(*,*) nangle

      tailname = '.mumps-rcs'
      scomb = trim(runname)//trim(tailname)
      open(33,file=scomb,status='unknown')


      tailname = '.HLU-rcs'
      scomb = trim(runname)//trim(tailname)
      open(44,file=scomb,status='unknown')

      do i=1,nangle
      read(33,*) theta, rcst, rcsh
      rcs1(i)=rcst
      read(44,*) theta, rcsv, rcsh
      rcs2(i)=rcsv
      enddo
     
      RMS=0
      do i=1,nangle
        print*,'dbrcs',abs(rcs1(i)-rcs2(i)),rcs1(i), rcs2(i),i
        RMS=(rcs1(i)-rcs2(i))**2+RMS
      enddo

      RMS=sqrt(RMS/real(nangle))

      print*,'RMS is:',RMS
      end
