       subroutine getsubmatrix_v(mvolnode,bfface,blpt,
     &            mvolele,mvoledge,msurele,msuredge,xyz,lv,le,lc,
     &            lp,bsurele,freq, mxmat, ep, mu, xi,
     &            hfemstore,delemax,hfemrealnum,
     &            submatrixi,submatrixj,submatrix,facetotal,
     &            kefemabc, kefemsrc,
     &            nsurele,nfface,nlpt,et,kz)
      implicit none

c.....Input Data
      integer mvolnode, mvolele, mvoledge, 
     &        msuredge, msurele
      integer facetotal, mxmat
      real freq
      complex ep(3,3,mxmat), mu(3,3,mxmat), xi(mxmat)
      integer kefemabc, kefemsrc

c......Output Data       
      real    xyz(3,mvolnode)
      integer lv(5,mvolele),le(6,mvolele)
      integer lp(3,msurele),lc(3,msurele)
      integer bsurele, nsurele
      integer hfemstore,delemax,hfemrealnum
      integer submatrixi(hfemstore),submatrixj(hfemstore)
      complex submatrix(hfemstore)
      integer bfface(4,facetotal),blpt(3,facetotal)
      integer nfface(4,facetotal),nlpt(3,facetotal)
      real et(*)
      complex kz

c......Working Variable
      integer i, j, k,iele,layerele,l
      integer,allocatable::fface(:,:),ipof(:),ipfo(:),ip(:)
      real,allocatable::rdis(:)
      integer facetop,num
      complex det,pivot(3)
      integer ipvot(3),index(3,2)
      complex mut(3,3),ept(3,3),qmt(3,3),xit(3)
      complex eptt,qmtt
      integer nm,isym,p1,p2,p3
      complex ele(6,6)
      integer,allocatable:: s1truehfemi(:),s1truehfemj(:)
      complex,allocatable:: s1truehfemmatrix(:)
      integer,allocatable::lps(:,:)
      integer num1i,num1j,num2i,num2j,num3i,num3j
      real small,dxx,dyy
      integer k1,k2,k3,nt,hfemnum
      integer facenum,bfacenum
      real rt1,rt2,rt3
      integer ii,jj      
      integer,allocatable:: ia(:)
      complex,allocatable:: a(:)
      complex kcmplx
ccccccccccccccccccccccc sort face, find boundary elements for each subdomain ccccccc

      do i=1,3
         xit(i)=0.
      end do

         bsurele=0
         nsurele=0
  
        layerele=mvolele
        allocate(fface(5,4*mvolele),rdis(4*mvolele),
     & lps(3,4*mvolele),
     & ipof(4*mvolele),ipfo(4*mvolele),ip(4*mvolele))

       do i=1,4*mvolele
         do j=1,4
          fface(j,i)=0
         enddo
         fface(5,i)=100
         do j=1,3
          lps(j,i)=0
         enddo
         rdis(i)=0
         ipof(i)=0
         ipfo(i)=0
         ip(i)=0
       enddo

       facetop=0
       do iele=1,layerele
          do j=1,4
             if (j.eq.1) then
               k1 = lv(1,iele)
               k2 = lv(3,iele)
               k3 = lv(2,iele)
             else if (j.eq.2) then
               k1 = lv(1,iele)
               k2 = lv(2,iele)
               k3 = lv(4,iele)
             else if (j.eq.3) then
               k1 = lv(1,iele)
               k2 = lv(4,iele)
               k3 = lv(3,iele)
             else if(j.eq.4) then
               k1 = lv(2,iele)
               k2 = lv(3,iele)
               k3 = lv(4,iele)
             end if
          call facetoedge(j,p1,p2,p3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         facetop = facetop+1

         fface(1,facetop) = k1
         fface(2,facetop) = k2
         fface(3,facetop) = k3
         fface(4,facetop) = iele

         lps(1,facetop)=le(p1,iele)
         lps(2,facetop)=le(p2,iele)
         lps(3,facetop)=le(p3,iele)
         enddo
       enddo
       num = 0
       do i=1,facetop
           num = num+1
           k1 = fface(1,i)
           k2 = fface(2,i)
           k3 = fface(3,i)
           call comi(k1,k2)
           call comi(k1,k3)
           call comi(k2,k3)

           rt1=k1
           rt2=k2
           rt3=k3

           rdis(num)=rt1*log(rt2)*log10(rt3)
       end do

       call hpsort(num,rdis,ipof,ipfo)

       nt = 1
       ip(nt) = ipfo(1)
       ipof(ipfo(1)) = nt

       do i = 2, num
          j=i-1
100            if (rdis(j).ne.rdis(i)) goto 200
                  num1i = fface(1,ipfo(i))
                  num2i = fface(2,ipfo(i))
                  num3i = fface(3,ipfo(i))

                  call comi(num1i,num2i)
                  call comi(num1i,num3i)
                  call comi(num2i,num3i)

                  num1j = fface(1,ipfo(j))
                  num2j = fface(2,ipfo(j))
                  num3j = fface(3,ipfo(j))
                  call comi(num1j,num2j)
                  call comi(num1j,num3j)
                  call comi(num2j,num3j)

        if((num1i.eq.num1j).and.(num2i.eq.num2j)
     &                     .and.(num3i.eq.num3j))then
          ipof(ipfo(i))=ipof(ipfo(j))
          fface(5,ipfo(i)) = 0
          fface(5,ipfo(j)) = 0
          goto 400
        else
          j=j-1
          if (j.le.0) goto 200
          goto 100
        endif

 200       nt=nt+1
           ipof(ipfo(i))=nt
           ip(nt)=ipfo(i)
 400      enddo

        facenum = nt
        bfacenum=0
        do i=1,4*layerele
         if(fface(5,i).ne.0)then
           if((lps(1,i).gt.(mvoledge-msuredge)).and.
     &        (lps(1,i).le.(mvoledge-kefemsrc)).and.
     &        (lps(2,i).gt.(mvoledge-msuredge)).and.
     &        (lps(2,i).le.(mvoledge-kefemsrc)).and.
     &        (lps(3,i).gt.(mvoledge-msuredge)).and.
     &        (lps(3,i).le.(mvoledge-kefemsrc))) then
              bfacenum=bfacenum+1
           else
           endif
         endif
        enddo
        bsurele=bfacenum

         bfacenum=0
         do i=1,4*layerele
         if(fface(5,i).ne.0)then
           if((lps(1,i).gt.(mvoledge-msuredge)).and.
     &        (lps(1,i).le.(mvoledge-kefemsrc)).and.
     &        (lps(2,i).gt.(mvoledge-msuredge)).and.
     &        (lps(2,i).le.(mvoledge-kefemsrc)).and.
     &        (lps(3,i).gt.(mvoledge-msuredge)).and.
     &        (lps(3,i).le.(mvoledge-kefemsrc))) then
              bfacenum=bfacenum+1
          do j=1,4
          bfface(j,bfacenum)=fface(j,i)
          enddo

          do j=1,3
          blpt(j,bfacenum)=lps(j,i)
          enddo
          else
         end if
       
         endif
         enddo


         bfacenum=0
        do i=1,4*layerele
         if(fface(5,i).ne.0)then
           if((lps(1,i).gt.(mvoledge-kefemsrc)).and.
     &        (lps(2,i).gt.(mvoledge-kefemsrc)).and.
     &        (lps(3,i).gt.(mvoledge-kefemsrc))) then
              bfacenum=bfacenum+1
           else
           endif
         endif
        enddo
        nsurele=bfacenum

         bfacenum=0
         do i=1,4*layerele
         if(fface(5,i).ne.0)then
           if((lps(1,i).gt.(mvoledge-kefemsrc)).and.
     &        (lps(2,i).gt.(mvoledge-kefemsrc)).and.
     &        (lps(3,i).gt.(mvoledge-kefemsrc))) then
              bfacenum=bfacenum+1
          do j=1,4
          nfface(j,bfacenum)=fface(j,i)
          enddo

          do j=1,3
          nlpt(j,bfacenum)=lps(j,i)
          enddo
          else
         end if
       
         endif
         enddo
      deallocate(fface,rdis,ipof,ipfo,ip,lps)
      print*,'bsurele',bsurele
      print*,'nsurele',nsurele
      print*,'msurele',msurele
      print*,'mvoledge:',mvoledge
      print*,'kefemabc:',kefemabc
      print*,'kefemsrc:',kefemsrc
      if ((bsurele+nsurele).ne.msurele)then
        print*,'(bsurele+nsurele).ne.msurele !!!!'
        stop
      endif

cccccccccccccccccccccccccc get submatrix for each subdomain ccccccccccccccccccccccccccc
      nm=3
      isym=0

      do i=1, mxmat
        do j=1,2
          do k=j+1,3
           if ((real(ep(j,k,i)).ne.real(ep(k,j,i))).or.
     &         (aimag(ep(j,k,i)).ne.aimag(ep(k,j,i))).or.
     &          (real(mu(j,k,i)).ne.real(mu(k,j,i))).or.
     &           (aimag(mu(j,k,i)).ne.aimag(mu(k,j,i)))) then
              isym=1
              print*,'The FEM Matrix are asymmetric'
              goto 91
           endif
          enddo
        enddo
      enddo
91    do i=1, mxmat
        call cmatin(nm,mu(1,1,i), nm, mut, nm, det,ipvot, index, pivot)
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      allocate (rdis(hfemstore),ipof(hfemstore),ipfo(hfemstore),
     &          ip(hfemstore))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1,hfemstore
            rdis(i)=0
            ipof(i)=0
            ipfo(i)=0
            ip(i)=0
         enddo
         hfemnum=0

         do iele=1,mvolele

           do j=1,3
           do k=1,3
              ept(j,k)=ep(j,k,lv(5,iele))
              qmt(j,k)=mu(j,k,lv(5,iele))
           enddo
           enddo
         call elevol1(ele,iele,freq,
     &              mvolnode,mvolele,
     &              xyz,lv,
     &              ept,qmt,xit)
         call assemfem(ele,iele,
     &                mvolele,mvolnode,mvoledge,msuredge,
     &                hfemstore,le,
     &                submatrixi,submatrixj,
     &                submatrix,hfemnum)
         enddo

         print*,'finish elevol and assem'

         eptt=ep(1,1,lv(5,nfface(4,1)))
         qmtt=mu(1,1,lv(5,nfface(4,1)))
         kcmplx=sqrt(eptt*qmtt)

         call cal_mod(freq, mvolnode, 
     &        mvoledge,xyz,nsurele,
     &        nfface,nlpt,
     &        et,
     &        eptt, qmtt, kz)

         allocate(ia(3*nsurele))
         allocate(a(3*nsurele))
         open(9988,file='elesursrcmportE.txt',status='unknown')  
         do iele=1,nsurele
            call elesursrcm(iele,mvolnode,mvoledge,mvolele,xyz,nsurele,
     &           nfface,nlpt,le,et,ia(3*(iele-1)+1),a(3*(iele-1)+1))
         enddo
         close(9988)

         do i=1,nsurele
          do l=1,nsurele
          do j=1,3
            do k=1,3
              hfemnum=hfemnum+1
              submatrixi(hfemnum)=ia(3*(i-1)+j)
              submatrixj(hfemnum)=ia(3*(l-1)+k)
              submatrix(hfemnum)=a(3*(i-1)+j)*a(3*(l-1)+k)
     &              *kz*cmplx(0.0,1.0)/qmtt!/sqrt(kcmplx)
              if((submatrixi(hfemnum).gt.mvoledge).or.
     &          (submatrixi(hfemnum).le.(mvoledge-kefemsrc)))then
                print*,'elesursrcm error!'
                stop
              endif
              if((submatrixj(hfemnum).gt.mvoledge).or.
     &          (submatrixj(hfemnum).le.(mvoledge-kefemsrc)))then
                print*,'elesursrcm error!'
                stop
              endif
            enddo
          enddo
         enddo
        enddo
         deallocate(ia)
         deallocate(a)
ccccccccccccccccccccccccccccccccc compress unknowns cccccccccccccccccccccccccccccccccccc
      small = 0.001

      do i=1,hfemnum
        rdis(i)=abs(submatrixi(i)+submatrixj(i))
      end do

      call hpsort(hfemnum,rdis,ipof,ipfo)

      nt = 1
      ip(nt) = ipfo(1)
      ipof(ipfo(1)) = nt

      do i= 2,hfemnum
         j=i-1
10     if (abs(rdis(j)-rdis(i)).gt.small) goto 20
          num1i = submatrixi(ipfo(i))
          num1j = submatrixj(ipfo(i))

          num2i = submatrixi(ipfo(j))
          num2j = submatrixj(ipfo(j))

           dxx = num1i-num2i
           dyy = num1j-num2j

        if ((abs(dxx).le.small).and.(abs(dyy).le.small)) then

             ipof(ipfo(i))=ipof(ipfo(j))
              goto 40
        else
            j=j-1
           if (j.le.0) goto 20
            goto 10
        endif
20        nt=nt+1
           ipof(ipfo(i))=nt
           ip(nt)=ipfo(i)
40      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      hfemrealnum= nt
      allocate (s1truehfemi(nt),s1truehfemj(nt),s1truehfemmatrix(nt))

      do i=1,hfemrealnum
         s1truehfemi(i)=0
         s1truehfemj(i)=0
         s1truehfemmatrix(i)=0.
      end do

      do i=1,hfemrealnum
         s1truehfemi(i) = submatrixi(ip(i))
         s1truehfemj(i) = submatrixj(ip(i))
      end do

      do i=1,hfemnum
       s1truehfemmatrix(ipof(i))=s1truehfemmatrix(ipof(i))
     &              +submatrix(i)
      end do

      do i=1,hfemstore
         submatrixi(i)= 0
         submatrixj(i)= 0
         submatrix(i)= 0.0
      enddo
     
      do i=1,hfemrealnum
         submatrixi(i)= s1truehfemi(i)
         submatrixj(i)= s1truehfemj(i)
         submatrix(i)= s1truehfemmatrix(i)
      enddo
 
      deallocate (s1truehfemi,s1truehfemj,s1truehfemmatrix)

      deallocate (rdis,ipof,ipfo,ip)



      RETURN
      END

        subroutine facetoedge(j,p1,p2,p3)
        implicit none
        integer j,p1,p2,p3
        if(j.eq.1)then
         p1=2
         p2=4
         p3=1
        else if(j.eq.2)then
         p1=1
         p2=5
         p3=3
        else if(j.eq.3)then
         p1=3
         p2=6
         p3=2
        else if(j.eq.4)then
         p1=4
         p2=6
         p3=5
        else
        print*,'error J!should <=4',j
        stop
        endif

        return
        end

