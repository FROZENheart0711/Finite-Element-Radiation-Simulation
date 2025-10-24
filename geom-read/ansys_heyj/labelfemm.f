      subroutine labelfemm
     &( mvolele, mtl, ke, nbnd, lv, rp, iqof, iqfo,
     &  ipt, ip, ipb, ldv, xyz, label )

      integer lv(9, mvolele), ldv(2, *)
       real   rp(*), xyz(3,*)
      integer iqof(*), iqfo(*), ipt(*), ip(*), ipb(*)
      logical t, t1, t2
      real,allocatable:: rp1(:)
      integer,allocatable:: iqof1(:), iqfo1(:)
      integer nbndt
      integer,allocatable:: ipbt(:)

        call idsame( mvolele, lv, iqfo(1), iqfo(1), 
     &               ie, je, k1, k2, l1, l2 )
        call nodeface(iqfo(1), kf1, kf2)

       if ((lv(5+kf1,ie).eq.label).or.(lv(5+kf2,ie).eq.label)) then
          ke = ke + 1
          ip(iqfo(1)) = ke
          ldv(1,ke) = k1
          ldv(2,ke) = k2
       else
       endif
      
       kff1 = 4*(ie-1)+kf1
       kff2 = 4*(ie-1)+kf2
       t1 = (lv(5+kf1,ie).eq.label).and.(ipt(kff2).eq.2)
       t2 = (lv(5+kf2,ie).eq.label).and.(ipt(kff1).eq.2)
       t  = t1.or.t2
       if (t.eqv..true.) then
          nbnd=nbnd+1
          ipb(nbnd)  = ke
       else
       endif   

      do 10 i=2, mtl
 
         i1=i
         i2=0
 30      i1=i1-1
 
         call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                ie, je, k1, k2, l1, l2 )
         call nodeface(iqfo(i), kf1, kf2)
         call nodeface(iqfo(i1),lf1, lf2)

         if ((lv(5+kf1,ie).eq.label).or.(lv(5+kf2,ie).eq.label)) then
             if (rp(i).ne.rp(i1)) then
               ke = ke + 1
               ip(iqfo(i)) = ke
               ldv(1,ke)=k1
               ldv(2,ke)=k2
               kff1=4*(ie-1)+kf1
               kff2=4*(ie-1)+kf2
               t1 = (lv(5+kf1,ie).eq.label).and.(ipt(kff2).eq.2)
               t2 = (lv(5+kf2,ie).eq.label).and.(ipt(kff1).eq.2)
               t  = t1.or.t2
               if (t.eqv..true.) then
                  nbnd=nbnd+1
                  ipb(nbnd)  =ke
               else
               endif   
            else
            endif
            if (rp(i).eq.rp(i1)) then
               if ((k1.ne.l1).or.(k2.ne.l2)) goto 30
 25            if ((lv(5+lf1,je).eq.label).or.
     &             (lv(5+lf2,je).eq.label)) then
                  ip(iqfo(i))=ip(iqfo(i1))
               else
 35               i1=i1-1
                  if (rp(i).eq.rp(i1)) then
                      call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                             ie, je, k1, k2, l1, l2 )
                      call nodeface(iqfo(i1),lf1, lf2) 
                      if ((k1.eq.l1).and.(k2.eq.l2)) then
                         goto 25
                      else
                         goto 35
                      endif
                  else
                  endif
                  i1=i-1
                  call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                         ie, je, k1, k2, l1, l2 )
                  call nodeface(iqfo(i1),lf1, lf2)
                  ke = ke + 1
                  ip(iqfo(i)) = ke
                  ldv(1,ke)=k1
                  ldv(2,ke)=k2
                  kff1=4*(ie-1)+kf1
                  kff2=4*(ie-1)+kf2
                  t1 = (lv(5+kf1,ie).eq.label).and.(ipt(kff2).eq.2)
                  t2 = (lv(5+kf2,ie).eq.label).and.(ipt(kff1).eq.2)
                  t  = t1.or.t2
                  if (t.eqv..true.) then
                     nbnd=nbnd+1
                     ipb(nbnd)  =ke
                  else
                  endif   
                  if ((k1.eq.l1).and.(k2.eq.l2)) then
                  ip(iqfo(i1)) = ke
 40               lf1=4*(je-1)+lf1
                  lf2=4*(je-1)+lf2
                     if ((ipt(lf1).eq.2).or.(ipt(lf2).eq.2)) then
                        nbnd=nbnd+1
                        ipb(nbnd)  =ke
                     else
                     endif
                  else
                  endif
 50               i1=i1-1
                  if (rp(i).ne.rp(i1)) goto 10
                  call idsame( mvolele, lv, iqfo(i), iqfo(i1), 
     &                         ie, je, k1, k2, l1, l2 )
                  if ((k1.eq.l1).and.(k2.eq.l2)) then
                     ip(iqfo(i1))= ke
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
               if ((lv(5+lf1,je).eq.label).or.
     &             (lv(5+lf2,je).eq.label)) then
                  ip(iqfo(i))=ip(iqfo(i1))
                  kf1=4*(ie-1)+kf1
                  kf2=4*(ie-1)+kf2
                 if ((ipt(kf1).eq.2).or.(ipt(kf2).eq.2)) then
                     nbnd=nbnd+1
                     ipb(nbnd)=ip(iqfo(i1))
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
ccccccccccccccccccc......new added

        allocate(rp1(nbnd))
        allocate(iqfo1(nbnd))
        allocate(iqof1(nbnd))
        allocate(ipbt(nbnd))

        do i=1,nbnd
            rp1(i)=ipb(i)
        enddo

        call hpsort(nbnd,rp1,iqof1,iqfo1)
        
        nbndt=1
        ipbt(1)=rp1(1)

        do i=2,nbnd
            i1=i
            i1=i1-1
            if(rp1(i).ne.rp1(i1))then
               nbndt=nbndt+1
               ipbt(nbndt)=rp1(i)
            endif         
        enddo

        do i=1,nbnd
            ipb(i)=0
        enddo

        nbnd=nbndt

        do i=1,nbnd
           ipb(i)=ipbt(i)
        enddo

        deallocate(rp1)
        deallocate(iqfo1)
        deallocate(iqof1)
        deallocate(ipbt)

cccccccccccccccccccccccccccccccccccccccccccc      

        return
        end
 
                 
                   
                   
 
