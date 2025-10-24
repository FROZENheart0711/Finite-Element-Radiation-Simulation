       subroutine find_face(mvolnode,mvolele,xyz,lv,
     &            sfface,ssurele, facetotal,mdef1, mdef2)

      implicit none

c.....Input 
      integer mvolnode, mvolele
      real    xyz(3,mvolnode)
      integer lv(5,mvolele)
      real    x_central,y_central,z_central
      real    radius_output1,radius_output2,radius_output3
      integer facetotal
c.....Output Data
      integer sfface(4,facetotal),ssurele
c......Working Variable
      integer i, j,iele
      integer mdef1, mdef2
      integer,allocatable::ipof(:),ipfo(:),ip(:)
      integer,allocatable::fface(:,:)
      real, allocatable:: rdis(:)
      real rt1,rt2,rt3
      integer nt,num1i,num1j,num2i,num2j,num3i,num3j
      integer facenum, k1, k2, k3, num
      integer ele1, ele2 

      print*,'mdef', mdef1, mdef2
      ele1=0
      ele2=0
      ssurele=0
      allocate(fface(6,4*mvolele))
      allocate(rdis(4*mvolele),ipof(4*mvolele),ipfo(4*mvolele),
     &ip(4*mvolele))

      do i=1,4*mvolele
       do j=1,6
          fface(j,i)=0
       enddo
       ip(i)=0
       ipof(i)=0
       ipfo(i)=0
       rdis(i)=0.0
      enddo

          facenum=0

          do i=1,mvolele

            do j=1,4

             if (j.eq.1) then
               k1 = lv(1,i)
               k2 = lv(2,i)
               k3 = lv(3,i)
             else if (j.eq.2) then
               k1 = lv(1,i)
               k2 = lv(2,i)
               k3 = lv(4,i)
             else if (j.eq.3) then
               k1 = lv(1,i)
               k2 = lv(3,i)
               k3 = lv(4,i)
             else if(j.eq.4) then
               k1 = lv(2,i)
               k2 = lv(3,i)
               k3 = lv(4,i)
             end if

              if((lv(5,i).eq.mdef1).or.(lv(5,i).eq.mdef2))then
               facenum=facenum+1
              if(lv(5,i).eq.mdef1)then
               if((j.eq.2).or.(j.eq.4))then
               fface(1,facenum) = k1
               fface(2,facenum) = k2
               fface(3,facenum) = k3
               else
               fface(1,facenum) = k1
               fface(2,facenum) = k3
               fface(3,facenum) = k2
               endif
              else
               if((j.eq.1).or.(j.eq.3))then
               fface(1,facenum) = k1
               fface(2,facenum) = k2
               fface(3,facenum) = k3
               else
               fface(1,facenum) = k1
               fface(2,facenum) = k3
               fface(3,facenum) = k2
               endif
              endif
               fface(4,facenum)=i
               fface(5,facenum)=lv(5,i)
              endif
           enddo
        enddo


        num = 0
        do i=1,facenum
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

10            if (rdis(j).ne.rdis(i)) goto 20
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
     &                    .and.(num3i.eq.num3j))then
          fface(6,ipfo(i))=fface(5,ipfo(j))
          fface(6,ipfo(j))=fface(5,ipfo(i))

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
        nt=0
        do i=1,facenum
        if((fface(5,(i)).eq.mdef1).and.(fface(6,(i)).eq.mdef2))then
          nt=nt+1
          do j=1,4
           sfface(j,nt)=fface(j,i)           
          enddo
          endif
       end do
       ssurele=nt
       print*,'facenum===', ssurele,facenum  
      deallocate(fface)
      deallocate(rdis,ipof,ipfo,ip)


      RETURN
      END


