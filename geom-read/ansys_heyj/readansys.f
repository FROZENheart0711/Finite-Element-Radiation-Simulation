      subroutine readansys
     &( mn, kls, klv,  
     &  xyz, ln, lv )

      real    xyz(3, *)
      integer ln(4, *), lv(9, *)
      integer tlv(9)
    
      read(1,*) mn
      do i=1, mn
         read(1,*) (xyz(j,i),j=1,3)
      enddo

      read(2,*) kls
      do i=1,kls
         read(2,*) (ln(j,i),j=1,4)
      enddo

      read(3,*) klv
      do i=1,klv
         read(3,*) (lv(j,i),j=1,9)
      enddo

      do i=1,klv
	
	if (lv(7,i).eq.1) then
	tlv(1) = lv(2,i)
	tlv(2) = lv(4,i)
	tlv(3) = lv(3,i)
	tlv(4) = lv(1,i)
	tlv(5) = lv(5,i)
	tlv(6) = 1
	tlv(7) = 0
	tlv(8) = 0
	tlv(9) = 0
	
	do j=1,9
	 lv(j,i) = tlv(j)
	end do

	else if (lv(8,i).eq.1) then
	tlv(1) = lv(3,i)
	tlv(2) = lv(4,i)
	tlv(3) = lv(1,i)
	tlv(4) = lv(2,i)
	tlv(5) = lv(5,i)
	tlv(6) = 1
	tlv(7) = 0
	tlv(8) = 0
	tlv(9) = 0
	
	do j=1,9
	 lv(j,i) = tlv(j)
	end do

	else if (lv(9,i).eq.1) then
	tlv(1) = lv(4,i)
	tlv(2) = lv(1,i)
	tlv(3) = lv(3,i)
	tlv(4) = lv(2,i)
	tlv(5) = lv(5,i)
	tlv(6) = 1
	tlv(7) = 0
	tlv(8) = 0
	tlv(9) = 0
	
	do j=1,9
	 lv(j,i) = tlv(j)
	end do
   
      else
	end if
	end do

      do i=1,kls
      do j=1,3
         k1=mod(j-1,3)+1
         k2=mod(j,3)+1
         k1=ln(k1,i)
         k2=ln(k2,i)
         d1=xyz(1,k1)-xyz(1,k2)
         d2=xyz(2,k1)-xyz(2,k2)
         d3=xyz(3,k1)-xyz(3,k2)
         d=sqrt(d1*d1+d2*d2+d3*d3)
         if (d.ge.0.3) then
            print*,'the length of the edge is larger than 0.3' 
            print*,i,k1,k2
         else
         endif
       enddo
       enddo

      return
      end

      
