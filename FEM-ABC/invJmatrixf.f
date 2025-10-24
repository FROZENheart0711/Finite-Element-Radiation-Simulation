       subroutine invJmatrixf(Jmatrix,IJmatrix,detJ)
	 implicit none
      INTEGER IERR

       integer i,j,k,comm_nodes,ipiv(3),info
       complex Jmatrix(3,3),sjmatrix(3,3)
       complex IJmatrix(3,3)
       real*8 detJ
       real*8 temp1(3),temp2(3)
       complex Imatrix(3,3)
       do i=1,3
       do j=1,3
       Imatrix(i,j)=0
       sjmatrix(i,j)=jmatrix(i,j)
       enddo
       enddo
       do i=1,3
       Imatrix(i,i)=1.0
       enddo
       detJ = Jmatrix(1,1)*(Jmatrix(2,2)*Jmatrix(3,3)
     &                     -Jmatrix(2,3)*Jmatrix(3,2))
     &       -Jmatrix(1,2)*(Jmatrix(2,1)*Jmatrix(3,3)
     &                     -Jmatrix(2,3)*Jmatrix(3,1))
     &       +Jmatrix(1,3)*(Jmatrix(2,1)*Jmatrix(3,2)
     &                     -Jmatrix(2,2)*Jmatrix(3,1))

       call cgesv(3,3,jmatrix,3,ipiv,imatrix,3,info)

       do i=1,3
       do j=1,3
       ijmatrix(i,j)=imatrix(i,j)
       jmatrix(i,j)=sjmatrix(i,j)
       enddo
       enddo

           return
           end










