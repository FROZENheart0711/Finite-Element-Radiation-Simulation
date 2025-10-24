       subroutine getfaceinfo(mvolele,mvoledge,
     & mvolface,lv,lf,le,facef,fedge)
       implicit none

        INTEGER mvolele, mvoledge,mvolface
        integer lv(5,mvolele),le(6,mvolele)
        integer facef(5,*),fedge(3,*),lf(4,mvolele)

c......Working Variable

        INTEGER i, j, k,ii,iele,layerele,id
        integer,allocatable::fface(:,:),ipof(:),lps(:,:),
     &  ipfo(:),ip(:)
        real,allocatable::rdis(:)
        integer facetop,num
        integer num1i,num1j,num2i,num2j,num3i,num3j
        real small,dxx,dyy
        integer k1,k2,k3,k4,nt,hfemnum,facenum
        real rt1,rt2,rt3
        integer p1,p2,p3
cccccccccccccccccccccccccccccc read material parameter ccccccccccccccccccccccccccccc
        layerele=mvolele
        allocate(fface(5,4*mvolele),rdis(4*mvolele),
     & lps(3,4*mvolele),
     & ipof(4*mvolele),ipfo(4*mvolele),ip(4*mvolele))

       do i=1,4*mvolele
         do j=1,5
          fface(j,i)=0
         enddo
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
               k2 = lv(2,iele)
               k3 = lv(3,iele)
             else if (j.eq.2) then
               k1 = lv(1,iele)
               k2 = lv(2,iele)
               k3 = lv(4,iele)
             else if (j.eq.3) then
               k1 = lv(1,iele)
               k2 = lv(3,iele)
               k3 = lv(4,iele)
             else if(j.eq.4) then
               k1 = lv(2,iele)
               k2 = lv(3,iele)
               k3 = lv(4,iele)
             end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         facetop = facetop+1
         fface(1,facetop) = k1
         fface(2,facetop) = k2
         fface(3,facetop) = k3
         fface(4,facetop) = iele
         call facetoedge(j,p1,p2,p3)
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
          fface(5,ipfo(i)) =  fface(4,ipfo(j))
          fface(5,ipfo(j)) =  fface(4,ipfo(i))
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

        mvolface=nt

        do i=1,nt
          do j=1,5
           facef(j,i) = fface(j,ip(i))
          end do

          do j=1,3
           fedge(j,i)=lps(j,ip(i))
          enddo

        end do

        do i=1,mvolele
          do j=1,4
            ii=j+(i-1)*4
            lf(j,i)=ipof(ii) !!!! face number of each element
          enddo
        enddo

      deallocate(fface,rdis,ipof,ipfo,ip,lps)

        RETURN
        END

