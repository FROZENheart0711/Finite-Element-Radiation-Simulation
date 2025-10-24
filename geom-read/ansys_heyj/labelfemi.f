      subroutine labelfemi
     &( mvolele, mtl, ke, lv, rp, iqof, iqfo,
     &  ip )

      integer lv(9, mvolele)
       real   rp(*)
      integer iqof(*), iqfo(*), ip(*)

      k=iqfo(1)
      if (ip(k).eq.0) then
        ke = ke + 1
        ip(iqfo(1)) = ke
      else
      endif

      do 10 i=2, mtl
 
         if (ip(iqfo(i)).ne.0) goto 10
         i1=i
 30      i1=i1-1
         if(i1.le.0)then
          print*,'in edge not found'
          goto 10
          endif
         if (rp(i).ne.rp(i1)) then
             ke = ke + 1
             ip(iqfo(i))=ke
         else
         endif
         if (rp(i).eq.rp(i1)) then
            call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                ie, je, k1, k2, l1, l2 )
            if ((k1.eq.l1).and.(k2.eq.l2)) then
               ip(iqfo(i))=ip(iqfo(i1))
            else
               goto 30
            endif
         else
         endif

10      continue

        return
        end
 
                 
                   
                   
