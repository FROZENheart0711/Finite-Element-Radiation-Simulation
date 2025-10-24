       subroutine assemfem(matrix,kip,
     &           mvolele,mvolnode,mvoledge,msuredge,
     &           hfemstore,le,
     &           hfemi,hfemj,hfemmatrix,hfemnum)

       complex  matrix(6,6)

        INTEGER mvolnode, mvolele, mvoledge,
     &          mvolstore, msurstore
        Integer hfemstore
        Integer msuredge

        INTEGER le(6,mvolele)

        integer hfemnum
        integer hfemi(hfemstore),hfemj(hfemstore)
        complex hfemmatrix(hfemstore)

        integer cofmatrix(6),kip
        integer i,j,k


        do i=1,6
          cofmatrix(i)=0
        end do

       do i=1,6
           k=i
           cofmatrix(k)=le(i,kip)
       end do

       do i=1,6
          do j=1,6
            hfemnum=hfemnum+1
            hfemi(hfemnum)=cofmatrix(i)
            hfemj(hfemnum)=cofmatrix(j)
            hfemmatrix(hfemnum)=matrix(i,j)
          end do
       end do


        return
        end
