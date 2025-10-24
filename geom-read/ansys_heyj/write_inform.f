        subroutine write_inform
     & ( nnode, npatch, nedge,
     &   xyzs,
     &   ipat, iedge,
     &   mnt, klt, ket, mn, kl, keout, kein, kc,                       
     &   xyzv,
     &   lv, le,
     &   ln, lp, ip, iq, ipposi, iqposi,
     &   mvolstore, msurstore,runname, 
     &   kefemabc, kefemsrc)

        integer ipat(4, *),iedge(5, *)
        integer lv(9, *),le(6, *)
        integer ln(4, *),lp(3, *)
	integer iq(*), ip(*)
         real   xyzs(3,*), xyzv(3,*)
        integer ipposi(*), iqposi(*)
        real ml,nl,klx,xx,yy,zz
        integer mnx,nn,kn
        integer isub
        integer tm,tn,tl
        integer writeFlg
        character*256 runname
        integer,allocatable::ipat2(:,:),iedge2(:,:)
        real,allocatable::middlenode(:,:)

        do i=1,klt
!        read(15,*) isub
        lv(6,i)=1 !isub+1
        enddo
!        close(15)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 do i=1,nedge
	   iedge(5,i)=0
	 enddo

        write(11,*) 'number of nodes and elements on surface' 
        write(11,*) nnode, npatch    

        write(11,*) 'coordiantes of nodes'
        do i=1,nnode
        write(11,*) (xyzs(j,i),j=1,3)
        end do

        write(11,*) 'relation of global label and local label for nodes'
        do i=1,npatch
        write(11,*) (ipat(j,i),j=1,3)
        end do

        write(12,*) 'number of edges on surface'
        write(12,*) nedge

        write(12,*) 'relation of global label and local label for edges'
        do i=1,nedge
        write(12,*) (iedge(j,i),j=1,5)
        end do
        
	  write(13,*) 'number of total nodes: mnt'
	  write(13,*) 'number of total elements: klt'
	  write(13,*) 'number of total edges: ket'
	  write(13,*) 'number of cavity and surface elements: kl'
	  write(13,*) 'number of edges on inner fem boundary: kein'
	  write(13,*) 'number of edges on outer fem boundary: keout'
	  write(13,*) 'number of nodes on intersection boundary: kc'

	  write(13,*) mnt,klt,ket,kl,kein,keout,kc,kefemabc,kefemsrc

	  write(13,*) 'coordinates of fem nodes'
	  do i=1,mnt
	  write(13,*) (xyzv(j,i),j=1,3)
	  enddo
 
	  write(13,*) 'rlt of global and local label for fem nodes'
	  do i=1,klt
		write(13,*) (lv(j,i),j=1,6)
	  enddo

	  write(13,*) 'rlt of global and local label for fem edges'
	  do i=1,klt
		write(13,*) (le(j,i),j=1,3)
		write(13,*) (le(j,i),j=4,6)
	  enddo

	  write(13,*) 'rlt of global and local label for nodes onsurface'
	  do i=1,kl
		write(13,*) (ln(j,i),j=1,3)
	  end do

	  write(13,*) 'rlt of global and local label for edges on surface'
	  do i=1,kl
	  	write(13,*) (lp(j,i),j=1,3)
	  end do

	  write(13,*) 'rlt of whole and cavity surface label for edges'
	  do i=1,keout
	  	write(13,*) iq(i)
	  enddo 	   

	  write(13,*) 'label of edges on inner fem surface'
	  do i=1,kein
	  	write(13,*) ip(i)
	  enddo    
 
	  write(13,*) 'position of row-index storage for total surface edges' 
	  do i=1,ket
		 write(13,*) ipposi(i)
	  enddo

	  write(13,*) 'position of row-index storage for cavity surface edges' 
	  do i=1,keout
		 write(13,*) iqposi(i)
	  enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccwrite .b* files for parallel mlfmacccccccccccccccccccc
       call WriteMomMax(nnode,npatch,nedge,1,runname);
       writeFlg=0
       call WriteNodCrd(writeFlg,runname,nnode,xyzs);
       allocate(ipat2(3,npatch))
       do i=1,npatch
           do j=1,3
            ipat2(j,i)=ipat(j,i)
           enddo
       enddo
       call WritePthNod(writeFlg,runname,npatch,ipat2);
       allocate(iedge2(4,nedge))
       do i=1,nedge
       do j=1,4
       iedge2(j,i)=iedge(j,i)
       enddo
       enddo
       call WriteEdgPth(writeFlg,runname,nedge, iedge2);
       allocate(middlenode(3,nedge))

      do i=1,nedge
        do j=1,3
          middlenode(j,i)=(xyzs(j,iedge(1,i))+xyzs(j,iedge(2,i)))/2.0
        end do
      end do
      open(44,file='middlenode',status='unknown')
      do i=1,nedge
      write(44,*),'middlenode',(middlenode(j,i),j=1,3)
      enddo
      close(44)
       call WriteEdgCtr(writeFlg,runname,nedge,middlenode);
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         print*,'mnt',mnt
	  write(14,*) mnt,klt,ket
	  write(14,*) mn, kl, keout,kein
	  write(14,*) msurstore
	  write(14,*) nnode, npatch, nedge
        return
        end

        subroutine getm(xx,yy,tm)
        real xx,yy
        integer mn,tm

        real tdepth

        if((xx.gt.0).and.(yy.gt.0))then
        tm=1
        else if((xx.gt.0).and.(yy.lt.0))then
        tm=4
        else if((xx.lt.0).and.(yy.gt.0))then
        tm=2
        else
        tm=3
        endif

        return
        end

        subroutine getn(xx,ml,mn,tm)
        real xx,ml
        integer mn,tm

        real tdepth

        if(xx.lt.0)then
        tm=1
        else
        tm=1
        endif

        return
        end







