        subroutine fcnd(mvolele, lv, i,k1,k2,k3)

        integer lv(9,mvolele)

        k=i-((i-1)/4)*4
!        k1=mod(k,4)+1
!        k2=mod(k+1,4)+1
!        k3=mod(k+2,4)+1
        if (k.eq.1) then
           k1=1
           k2=3
           k3=2
        else if (k.eq.2) then
           k1=1
           k2=2
           k3=4
        else if (k.eq.3) then
           k1=1
           k2=4
           k3=3
        else if (k.eq.4) then
           k1=2
           k2=3
           k3=4
        end if
 
        k=(i-1)/4+1
        k1=lv(k1,k)
        k2=lv(k2,k)
        k3=lv(k3,k)

        call comi(k1,k2)
        call comi(k1,k3)
        call comi(k2,k3)

        return
        end
