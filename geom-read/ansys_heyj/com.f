        subroutine com(i,j,iflag)

        if (i.le.j) then
        k=i
        i=j
        j=k
        iflag=-1
        else
        iflag=1
        endif

        return
        end

