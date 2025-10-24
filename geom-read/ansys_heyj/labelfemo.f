      subroutine labelfemo
     &( mvolele, mtl, ke, lv, rp, iqof, iqfo,
     &  ip, ipt )

      integer lv(9, mvolele)
       real   rp(*)
      integer iqof(*), iqfo(*), ip(*), ipt(*)

      k=iqfo(1)
      call idsame( mvolele, lv, k, k, 
     &             ie, je, k1, k2, l1, l2 )
      call nodeface(k, kf1, kf2)
      if (ip(k).eq.0) then
         kf1  = 4*(ie-1)+kf1
         kf2  = 4*(ie-1)+kf2
         if ((ipt(kf1).eq.2).or.(ipt(kf2).eq.2)) then
            ke = ke + 1
            ip(iqfo(1)) = ke
         else
         endif
      else
      endif

      do 10 i=2, mtl
 
         if (ip(iqfo(i)).ne.0) goto 10
         i1=i
 30      i1=i1-1
 
         call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                ie, je, k1, k2, l1, l2 )
         call nodeface(iqfo(i), kf1, kf2)
         call nodeface(iqfo(i1),lf1, lf2)
         kf1=4*(ie-1)+kf1
         kf2=4*(ie-1)+kf2
         lf1=4*(je-1)+lf1
         lf2=4*(je-1)+lf2
         if ((ipt(kf1).eq.2).or.(ipt(kf2).eq.2)) then
            if (rp(i).ne.rp(i1)) then
               ke = ke + 1
               ip(iqfo(i))=ke
            else
            endif
            if (rp(i).eq.rp(i1)) then
               if ((k1.ne.l1).or.(k2.ne.l2)) goto 30
 25            if ((ipt(lf1).eq.2).or.(ipt(lf2).eq.2)) then
                  ip(iqfo(i))=ip(iqfo(i1))
               else
 35               i1=i1-1
                  if (rp(i).eq.rp(i1)) then
                      call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                             ie, je, k1, k2, l1, l2 )
                      call nodeface(iqfo(i1),lf1, lf2) 
                      lf1=4*(je-1)+lf1
                      lf2=4*(je-1)+lf2                    
                      if ((k1.eq.l1).and.(k2.eq.l2)) then
                         goto 25
                      else
                         goto 35
                      endif
                  else
                  endif
                  i1=i-1
                  ke = ke + 1
                  ip(iqfo(i)) = ke 
50                if (rp(i).ne.rp(i1)) goto 10
                   call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                       ie, je, k1, k2, l1, l2 )
                  if ((k1.eq.l1).and.(k2.eq.l2)) then
                     ip(iqfo(i1))= ke
                  else
                  endif
                  i1=i1-1
                  goto 50
               endif
            else
            endif
         else
 60         if (rp(i).ne.rp(i1)) goto 10

            call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                   ie, je, k1, k2, l1, l2 )
            call nodeface(iqfo(i), kf1, kf2)
            call nodeface(iqfo(i1),lf1, lf2)

            if ((k1.eq.l1).and.(k2.eq.l2)) then
               if ((ipt(lf1).eq.2).or.(ipt(lf2).eq.2)) then
                  ip(iqfo(i))=ip(iqfo(i1))
               else
                  i1=i1-1
                  goto 60
               endif
            else
               i1=i1-1
               goto 60
           endif
        endif
10      continue

        return
        end
 
                 
                   
                   
