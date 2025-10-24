      subroutine read_geom_fem
     & (mvolnode, mvolele, mvoledge,
     &  msurnode, msurele, msuredge, kein, kc,
     &  mvolstore, msurstore,
     &  xyz, lv, le, lp, lc, ipin,
     &  kefemabc,kefemsrc)
      implicit none
c.....Input Data
      integer mvolnode, mvolele, mvoledge, 
     &        msurnode, msuredge, msurele,
     &        mvolstore, msurstore,kein

c......Output Data       
      real    xyz(3,mvolnode)
      integer lv(5,mvolele),le(6,mvolele)
      integer lp(3,msurele),lc(3,msurele)
      integer ipin(kein)

c......Working Variable
      integer mn, kv, ke, klout, keout, kc,keint
      integer kefemabc,kefemsrc
      integer  i, j, k,ii

c.......'number of total nodes: mn'
c.......'number of total elements: kv'
c.......'number of total edges: ke'
c.......'number of cavity surface elements: klout'
c.......'number of edges on inner FEM boundary: kein'
c.......'number of edges on outer FEM boundary: keout'
c.......'number of nodes on intersection boundary: kc'
        
        read(5,*) 
        read(5,*) 
        read(5,*) 
        read(5,*) 
        read(5,*) 
        read(5,*) 
        read(5,*) 
        read(5,*) mn,kv,ke,klout,keint,keout,kc,kefemabc,kefemsrc

        if((kefemabc+kefemsrc).ne.msuredge)then
         write(*,*) 'abc+src.NE.msuredge'
         print*, kefemabc,kefemsrc,msuredge
         stop
        endif

        if (mn.ne.mvolnode) then
         write(*,*) 'mn.NE.mvolnode'
         print*, mn,mvolnode
         stop
        else
        endif

        if (kv.ne.mvolele) then
          write(*,*) 'kv.NE.mvolele'
          print*, kv,mvolele
          stop
        else
        endif

        if (ke.ne.mvoledge) then
          write(*,*) 'ke.NE.mvoledge'
          stop
        else
        endif
        if(keint.ne.kein) then
          write(*,*) 'keint.NE.kein'
          stop
        else
        end if

        read(5,*)
        do i=1,mn
           read(5,*) (xyz(j,i),j=1,3)
        enddo
        read(5,*) 
        do i=1,kv
         read(5,*) (lv(j,i),j=1,5)
        enddo

        read(5,*) 
        do i=1,kv
           read(5,*) (le(j,i),j=1,3)
           read(5,*) (le(j,i),j=4,6)
        end do

        read(5,*) 
        do i=1,klout
           read(5,*) (lc(j,i),j=1,3)
        enddo

        read(5,*)
        do i=1,klout
           read(5,*) (lp(j,i),j=1,3)
        end do

        read(5,*) 
        do i=1,keout
           read(5,*) 
        enddo

        read(5,*)
        do i=1,kein
           read(5,*) ipin(i)
        enddo

        close(5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        write(*,*)  'cccccccccccccccccccccccccccccccccccccc'
        write(*,*)  'Information for FEM' 
        write(*,*)  'cccccccccccccccccccccccccccccccccccccc'
        write(*,*) 'mvoledge, msuredge, kein, kc =',
     &              mvoledge, msuredge, kein, kc
        RETURN
        END
