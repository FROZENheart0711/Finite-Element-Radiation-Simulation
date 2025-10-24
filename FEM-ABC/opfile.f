      subroutine OPFILE(runname,clength,mxmat)
      implicit none
      integer clength,mxmat
      character(len=clength)::runname
      character*256 tailname, scomb

       tailname='.amp'
       scomb=trim(runname)//trim(tailname)
       open(15,file=scomb,status='unknown')
       
       

       tailname = '.fem'
       scomb = trim(runname)//trim(tailname)
       open(5, file=scomb, status='unknown')

       tailname='.rcs'

       scomb = trim(runname)//trim(tailname)
       open(10, file=scomb, status='unknown')

       tailname='.fmm_max'
       scomb = trim(runname)//trim(tailname)
       open(2, file=scomb, status='unknown')

       tailname='.inp'
       scomb=trim(runname)//trim(tailname)
       open(1,file=scomb,status='unknown')
       read(1,*)
       read(1,*) mxmat

       return
      end

      subroutine CASUM(N,DX,INCX,csum)
*     .. Scalar Arguments ..
      INTEGER INCX,N
      real csum
*     ..
*     .. Array Arguments ..
      COMPLEX DX(*)
      COMPLEX DTEMP
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CABS,MOD
*     ..
      CSUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DTEMP = DTEMP + CABS(DX(I))
   10 CONTINUE
      CSUM = DTEMP
      RETURN
   20 M = MOD(N,6)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + CABS(DX(I))
   30 CONTINUE
      IF (N.LT.6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
          DTEMP = DTEMP + CABS(DX(I)) + CABS(DX(I+1)) + CABS(DX(I+2)) +
     +            CABS(DX(I+3)) + CABS(DX(I+4)) + CABS(DX(I+5))
   50 CONTINUE
   60 CSUM = DTEMP
      RETURN
      END


