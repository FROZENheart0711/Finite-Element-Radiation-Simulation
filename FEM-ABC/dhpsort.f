      subroutine hpsort(n,ra,iqof,iqfo)
        integer n
        real ra(n)
        integer i,ir,j,l,iqof(n),iqfo(n)
       real rra
        do 5 i=1,n
        iqof(i)=i
5       iqfo(i)=i
        if (n.lt.2) return
        l=n/2+1
        ir=n
10      continue
        if (l.gt.1) then
        l=l-1
        rra=ra(l)
        else
        rra=ra(ir)
        ra(ir)=ra(1)
        it1=iqfo(1)
        it2=iqfo(ir)
        iqfo(1)=it2
        iqfo(ir)=it1
        iqof(it2)=1
        iqof(it1)=ir
        ir=ir-1
        if (ir.eq.1) then
        ra(1)=rra
        return
        endif
        endif
        i=l
        j=l+l
        it1=iqfo(i)
20      if (j.le.ir) then
           if (j.lt.ir) then
             if (ra(j).lt.ra(j+1)) j=j+1
             endif
             if (rra.lt.ra(j)) then
             ra(i)=ra(j)
             it2=iqfo(j)
             iqof(it2)=i
             iqfo(i)=it2
             i=j
             j=j+j
             else
             j=ir+1
             endif
        goto 20
        endif
        ra(i)=rra
             iqfo(i)=it1
             iqof(it1)=i
        goto 10
        end


        
      subroutine dhpsort(n,ra,iqof,iqfo)
        integer n
        real*8 ra(n)
        integer i,ir,j,l,iqof(n),iqfo(n)
        real*8 rra
        do 5 i=1,n
        iqof(i)=i
5       iqfo(i)=i
        if (n.lt.2) return
        l=n/2+1
        ir=n
10      continue
        if (l.gt.1) then
        l=l-1
        rra=ra(l)
        else
        rra=ra(ir)
        ra(ir)=ra(1)
        it1=iqfo(1)
        it2=iqfo(ir)
        iqfo(1)=it2
        iqfo(ir)=it1
        iqof(it2)=1
        iqof(it1)=ir
        ir=ir-1
        if (ir.eq.1) then
        ra(1)=rra
        return
        endif
        endif
        i=l
        j=l+l
        it1=iqfo(i)
20      if (j.le.ir) then
           if (j.lt.ir) then
             if (ra(j).lt.ra(j+1)) j=j+1
             endif
             if (rra.lt.ra(j)) then
             ra(i)=ra(j)
             it2=iqfo(j)
             iqof(it2)=i
             iqfo(i)=it2
             i=j
             j=j+j
             else
             j=ir+1
             endif
        goto 20
        endif
        ra(i)=rra
             iqfo(i)=it1
             iqof(it1)=i
        goto 10
        end
