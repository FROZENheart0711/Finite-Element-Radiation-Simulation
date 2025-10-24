      subroutine read_fem_inp(mvolnode,mvolele,mvoledge,msurnode, 
     &           msurele, msuredge,msurstore,maxnode, maxpatch, 
     &           maxedge,kein, ipol, tyrcs, freq, tol,
     &           thetai, phii, thetas, phis, mxmat, ep, mu, xi,
     &           mdef1, mdef2)
      implicit none
      integer msurnode, msurele, msuredge,
     &        msurstore,mvoledge,mvolnode,mvolele,
     &        maxnode, maxpatch, maxedge,kein
      integer i,j,k
      integer nargs,iargc
      character*256 scomb,tailname

      integer mxmat,ipol,tyrcs
      real    freq,tol
      integer thetai(3),phii(3),thetas(3),phis(3)
      complex ep(3,3,mxmat),mu(3,3,mxmat),xi(mxmat)
      integer mdef1, mdef2
cccccccccccccccccccccccc Start ccccccccccccccccccccccccccccccccccccc

       read(2,*) mvolnode, mvolele, mvoledge
       read(2,*) msurnode, msurele, msuredge,kein
       read(2,*) msurstore
       read(2,*) maxnode, maxpatch, maxedge
       close(2)

       read(1,*)
       read(1,*) freq
       read(1,*)
       read(1,*) ipol
       read(1,*)
       read(1,*) thetai(1),thetai(2),thetai(3)
       read(1,*)
       read(1,*) phii(1),phii(2),phii(3)
       read(1,*)
       read(1,*) 
       read(1,*) tyrcs
       read(1,*)
       read(1,*) thetas(1),thetas(2),thetas(3)
       read(1,*)
       read(1,*) phis(1),phis(2),phis(3)
       read(1,*)
       read(1,*) mdef1, mdef2
       read(1,*)

!       allocate(ep(3,3,mxmat),mu(3,3,mxmat),xi(mxmat))

       do k=1,mxmat
          read(1,*)
          do i=1,3
           read(1,*)(ep(i,j,k),j=1,3)
          enddo
          read(1,*)
          do i=1,3
           read(1,*)(mu(i,j,k),j=1,3)
          enddo
          read(1,*)
          read(1,*) xi(k)
       end do
       read(1,*)
       read(1,*) tol
       close(1)


       
       return
       END

