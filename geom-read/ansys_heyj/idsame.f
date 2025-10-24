      subroutine idsame
     &( mvolele, lv, i, i1, ie, je, k1, k2, j1, j2)
      
       integer lv(9, mvolele)

         ie=(i-1)/6 + 1
         in=i-(ie-1)*6
         call labelinv(in, kk1, kk2)
         k1=lv(kk1,ie)
         k2=lv(kk2,ie)
         call comi(k1,k2)

         je=(i1-1)/6 + 1
         jn=i1-(je-1)*6
         call labelinv(jn, jj1, jj2)
         j1=lv(jj1,je)
         j2=lv(jj2,je)
         call comi(j1,j2)

       return
       end
