      subroutine fem(mvolnode,mvolele,mvoledge,msuredge,msurele,
     &               bsurele, hfemstore, dhfemnum, kein, ipin, freq,
     &               mxmat, ep, mu, xi,
     &               xyz, lv, le, lp, lc, bfface, blpt, edgesign2,
     &               submatrixi, submatrixj, submatrix,
     &               kefemabc, kefemsrc,
     &               nsurele,nfface,nlpt,et,kz)
      implicit none

      integer mvolnode, mvolele, mvoledge,bsurele
      integer msuredge, msurele,kein
      integer mxmat
      complex ep(3,3,mxmat), mu(3,3,mxmat), xi(mxmat)
      real    xyz(3,mvolnode)
      integer lv(5,mvolele),le(6,mvolele)
      integer lp(3,msurele),lc(3,msurele), ipin(kein)
      integer hfemstore,dhfemnum
      integer submatrixi(*),submatrixj(*)
      complex submatrix(*)
      integer bfface(4,*),blpt(3,*)
      integer edgesign2(*)
      real freq
      real et(*)
      complex kz

      integer facetotal, delemax
      integer i,j
      integer kefemabc, kefemsrc
      integer nsurele
      integer nfface(4,*),nlpt(3,*)
cccccccccccccccccc memory for preddm ccccccccccccccccccccccccccccccccccc
      delemax=mvolele
      hfemstore=delemax*36+msurele**2*9
      facetotal=4*mvolele

      if (kein.gt.0) then
      do i=1,kein
         edgesign2(ipin(i))=1
      end do
      end if
ccccccccccccccc ccccccccccccccccccccccccccccccccccccccccccccc
!!!! get volume matrix 
      call getsubmatrix_v(mvolnode, bfface,blpt,
     &     mvolele,mvoledge,msurele,msuredge,xyz,lv,le,lc,
     &     lp,bsurele,freq, mxmat, ep, mu, xi,
     &     hfemstore,delemax,dhfemnum,
     &     submatrixi,submatrixj,submatrix,
     &     facetotal, kefemabc, kefemsrc,
     &     nsurele,nfface,nlpt,et,kz)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!! get surface matrix
      call getsubmatrix_s(mvolnode,mvolele,bfface,blpt,
     &     bsurele,xyz,freq,
     &     hfemstore,dhfemnum,
     &     submatrixi,submatrixj,submatrix,edgesign2,
     &     facetotal,kein)

      return
      end


