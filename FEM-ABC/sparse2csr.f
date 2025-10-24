      subroutine sparse2csr(mvoledge, nzp, mati, matj, mat,
     & iqvol, ipvol, volmat)

      implicit none
ccc input 
      integer mvoledge
      integer nzp
      integer mati(*), matj(*), iqvol(*), ipvol(*)
      complex mat(*), volmat(*)

      integer i,j, nii, maxsys, ii, k, kk

      real rp(800)
      integer iqfo(800),iqof(800)
      complex a(800)
      integer mtl


      call zcoocsr(mvoledge,nzp,mat,mati,matj,volmat,ipvol,iqvol)

ccccccccccccc resort ccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,mvoledge
        mtl=0
        do j=iqvol(i),iqvol(i+1)-1
          mtl=mtl+1
          rp(mtl)=ipvol(j)
          a(mtl)=volmat(j)
        enddo

        if(mtl.gt.100)then
c         print*,'Size 20 is not enough! increase!',mtl
c         stop
        endif

        call hpsort(mtl,rp,iqof,iqfo)

        mtl=0
        do j=iqvol(i),iqvol(i+1)-1
          mtl=mtl+1
          volmat(j)=a(iqfo(mtl))
          ipvol(j)=int(rp((mtl)))
        enddo
      enddo
110   continue

      return
      end


      subroutine zcoocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
c-----------------------------------------------------------------------
      complex a(nnz),ao(nnz),x
      integer ir(nnz),jc(nnz),jao(nnz),iao(nrow+1)
c------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
c------------- end of coocsr -------------------------------------------
c-----------------------------------------------------------------------
      end

