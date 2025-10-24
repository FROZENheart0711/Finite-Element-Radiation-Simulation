      subroutine labelfemm
     &( mvolele, mtl, kefemm, nbnd, lv, rp, iqof, iqfo,
     &  ip, ipt, ldv )

      integer lv(6, mvolele), ldv(2, *)
       real   rp(*)
      integer iqof(*), iqfo(*), ip(*), ipt(*)

       ip(iqfo(1)) = 1
       kefemm      = 1

      do 10 i=2, mtl
 
         i1=i
 30      i1=i1-1
 
         call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                ie, je, k1, k2, l1, l2 )
         call nodeface(iqfo(i), kf1, kf2)
         call nodeface(iqfo(i1),lf1, lf2)

         if ((kf1.eq.lv(6,ie)).or.(kf2.eq.lv(6,ie))) then
            if (rp(i).ne.rp(i1)) then
               kefemm=kefemm+1
               ip(iqfo(i))=kefemm
            else
            endif
            if (rp(i).eq.rp(i1)) then
               if ((k1.ne.l1).or.(k2.ne.l2)) goto 30
               if ((lf1.eq.lv(6,je)).or.(lf2.eq.lv(6,je))) then
                  ip(iqfo(i))=ip(iqfo(i1))
               else
                  kefemm = kefemm + 1
                  ip(iqfo(i)) = kefemm
                  ip(iqfo(i1))= kefemm
 40               lf1=4*(je-1)+lf1
                  lf2=4*(je-1)+lf2
                  if ((ipt(lf1).eq.1).or.(ipt(lf2).eq.1)) then
                     nbnd=nbnd+1
                     ldv(1,nbnd)=k1
                     ldv(2,nbnd)=k2
                  else
                  endif
 50               i1=i1-1
                  if (rp(i).ne.rp(i1)) goto 10
                  call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                       ie, je, k1, k2, l1, l2 )
                  if ((k1.eq.l1).and.(k2.eq.l2)) then
                     ip(iqfo(i1))= kefemm
                     call nodeface(iqfo(i1),lf1, lf2)
                     goto 40
                  else
                  endif
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
               if ((lf1.eq.lv(6,je)).or.(lf2.eq.lv(6,je))) then
                  ip(iqfo(i))=ip(iqfo(i1))
                  kf1=4*(ie-1)+kf1
                  kf2=4*(ie-1)+kf2
                  if ((ipt(kf1).eq.1).or.(ipt(kf2).eq.1)) then
                     nbnd=nbnd+1
                     ldv(1,nbnd)=k1
                     ldv(2,nbnd)=k2
                  else
                  endif
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
        endif
 
                 
                   
                   
