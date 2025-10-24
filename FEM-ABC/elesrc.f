       subroutine elesrc(kip,xyz,p_srcspfrt,p_mvolnodefrt,wpsp,
     &     k0,ur,ep,
     &     att,btt,btz,bzt,bzz)
        
        
        

        implicit none
!.......Input Data
        integer kip
        real  xyz(3,*)
        integer*8 p_srcspfrt
        integer*8 p_mvolnodefrt
        integer wpsp(4,*)
        real*8 k0
        complex ep,ur
!.......Output Data
        complex*16 att(8,8)
        complex*16 btt(8,8)
        complex*16 btz(8,6)
        complex*16 bzt(6,8)
        complex*16 bzz(6,6)
!.......Working Variables
        complex*16 ele(8,8),ele2(8,8)
        integer ne1(3),ne2(3)
        real*8    length1(3),length2(3)
        real*8    xx1(3),yy1(3),zz1(3)
        real*8    xx2(3),yy2(3),zz2(3)
        real*8    xn1(3),yn1(3)
        real*8 x1,x2,x3,y1,y2,y3,z1,z2,z3
        real*8    area1,area2
        integer i,j,k,i1,i2,i3,j1,j2,j3,k1,k2,k3
        integer ii,jj
        integer ktest,ktest1,ktest2
c        real    rcurn1,rcurn2
        integer id1
        real rdis(3)
        integer ipof(3),ipfo(3)
        integer cofvertex(3),cofedge(3)
        integer orginedge(3,3)
        real*8 uu(7),vv(7),nweigh(7)
        real*8 u,v,weigh,detJ,detJ1
        real volume
        complex B(2,8),CurlB(3,8)
        complex N(6),gradN(2,6)
        complex Jmatrix(2,3),IJmatrix(3,2),TIJmatrix(2,3)
        complex Jmatrix1(3,3),IJmatrix1(3,3),TJmatrix1(3,3)
        complex*16 attmtx1(8,8),attmtx2(8,8)
        complex*16 bzzmtx1(6,6),bzzmtx2(6,6)
        complex*16 btzmtx(8,6)
        complex temp1(2),temp2(2),temp3(3),temp4(3),temp
        complex temp5(3)
        complex temp11(3),temp22(3),temp33(3),temp44(3)
        complex dotmul_c
        integer ksur,kk
        integer cofmatrix1(8),cofmatrix2(8),cofmatrix(8)
c.......pre
         uu(1)= 1./3.
         vv(1)= 1./3.
         
         uu(2)= 1.
         vv(2)= 0.
         
         uu(3)= 0.
         vv(3)= 1.

         uu(4)= 0.
         vv(4)= 0.
         
         uu(5)= 1./2.
         vv(5)= 1./2.

         uu(6)= 0.
         vv(6)= 1./2.

         uu(7)= 1./2.
         vv(7)= 0.

         nweigh(1) = 9./20.
         nweigh(2) = 1./20.
         nweigh(3) = 1./20.
         nweigh(4) = 1./20.
         nweigh(5) = 2./15.
         nweigh(6) = 2./15.
         nweigh(7) = 2./15.

        

         do i=1,8
            do j=1,2
            B(j,i)=0.
            enddo
         enddo

         do i=1,8
            do j=1,3
            CurlB(j,i)=0.
            enddo
         enddo

        do i=1,3
           ne1(i)=wpsp(i,kip)
           xx1(i)=xyz(1,ne1(i))
           yy1(i)=xyz(2,ne1(i))
           zz1(i)=xyz(3,ne1(i))           
        enddo     
         
        volume=1.0/2.0

        attmtx1=0.
        attmtx2=0.
        bzzmtx1=0.
        bzzmtx2=0.
        btzmtx=0.

c..............begin! 

        do i=1,7

        u=uu(i)
        v=vv(i)
        weigh=nweigh(i)
        x1=xx1(1)
        x2=xx1(2)
        x3=xx1(3)
        y1=yy1(1)
        y2=yy1(2)
        y3=yy1(3)
        z1=zz1(1)
        z2=zz1(2)
        z3=zz1(3)

        detJ1=sqrt( ((y1-y3)*(z2-z3)-(z1-z3)*(y2-y3))**2+
     &         ((z1-z3)*(x2-x3)-(x1-x3)*(z2-z3))**2+
     &         ((x1-x3)*(y2-y3)-(y1-y3)*(x2-x3))**2)

        Jmatrix1(1,1)=x1-x3
        Jmatrix1(1,2)=y1-y3
        Jmatrix1(1,3)=z1-z3
        Jmatrix1(2,1)=x2-x3
        Jmatrix1(2,2)=y2-y3
        Jmatrix1(2,3)=z2-z3
        Jmatrix1(3,1)=1/detJ1*((y1-y3)*(z2-z3)-
     &           (z1-z3)*(y2-y3) )
        Jmatrix1(3,2)=1/detJ1*((z1-z3)*(x2-x3)-
     &           (x1-x3)*(z2-z3) )
        Jmatrix1(3,3)=1/detJ1*((x1-x3)*(y2-y3)-
     &           (y1-y3)*(x2-x3) )
         do ii=1,2
            do jj=1,3
               Jmatrix(ii,jj)=Jmatrix1(ii,jj)
            enddo
         enddo

         call invJmatrixf(Jmatrix1,IJmatrix1,detJ)

         detJ=detJ1

         do ii=1,3
             do jj=1,2
                IJmatrix(ii,jj)=IJmatrix1(ii,jj)
             enddo
         enddo

        do ii=1,2
         do jj=1,3
            TIJmatrix(ii,jj)=IJmatrix(jj,ii)
         enddo
        enddo

        do ii=1,3
         do jj=1,3
            TJmatrix1(ii,jj)=Jmatrix1(jj,ii)
         enddo
        enddo

        call devaBsur(u,v,B)

        call devaCurlBsur(u,v,CurlB)

        call devaNsur(u,v,N)

        call devagradNsur(u,v,gradN)

c........att curlE dot curlE!

        do ii=1,8
         do jj=1,8
            do k=1,3
               temp11(k)=CurlB(k,ii)
               temp22(k)=CurlB(k,jj)
            enddo       

            do i1=1,3 
               temp33(i1) = temp11(1)*Jmatrix1(1,i1)+
     &                      temp11(2)*Jmatrix1(2,i1)+
     &                      temp11(3)*Jmatrix1(3,i1)
            end do        

            do i1=1,3   
               temp44(i1) = temp22(1)*TJmatrix1(i1,1)+
     &                      temp22(2)*TJmatrix1(i1,2)+
     &                      temp22(3)*TJmatrix1(i1,3)
            end do
                
            attmtx2(ii,jj)=attmtx2(ii,jj)+weigh*volume*
     &            (temp33(1)*temp44(1)+temp33(2)*temp44(2)+
     &            temp33(3)*temp44(3))/dabs(detJ)
         enddo
        enddo

c........att E dot E!
        do ii=1,8
         do jj=1,8
            do k=1,2
               temp1(k)=B(k,ii)
               temp2(k)=B(k,jj)
            enddo


            do i1=1,3 
               temp3(i1) = temp1(1)*IJmatrix(i1,1)+
     &                     temp1(2)*IJmatrix(i1,2)
            end do

            do i1=1,3   
               temp4(i1) = temp2(1)*IJmatrix(i1,1)+
     &                     temp2(2)*IJmatrix(i1,2)
            end do
              
            attmtx1(ii,jj)=attmtx1(ii,jj)+weigh*volume*abs(detJ)*
     &            (temp3(1)*temp4(1)+temp3(2)*temp4(2)+
     &            temp3(3)*temp4(3))
         enddo
        enddo

c........bzz N dot N!
        do ii=1,6
         do jj=1,6
            bzzmtx1(ii,jj)=bzzmtx1(ii,jj)+weigh*volume*
     &      N(ii)*N(jj)*dabs(detJ)
         enddo
        enddo

c........bzz gradN dot gradN!
        do ii=1,6
         do jj=1,6
            do k=1,2
               temp1(k)=gradN(k,ii)
               temp2(k)=gradN(k,jj)
            enddo


            do i1=1,3 
               temp3(i1) = temp1(1)*TIJmatrix(1,i1)+
     &                     temp1(2)*TIJmatrix(2,i1)
            end do

            do i1=1,3   
               temp4(i1) = temp2(1)*IJmatrix(i1,1)+
     &                     temp2(2)*IJmatrix(i1,2)
            end do
              
            bzzmtx2(ii,jj)=bzzmtx2(ii,jj)+weigh*volume*
     &            (temp3(1)*temp4(1)+temp3(2)*temp4(2)+
     &            temp3(3)*temp4(3))*dabs(detJ)
         enddo
        enddo

c........btz E dot gradN!
        do ii=1,8
         do jj=1,6
            do k=1,2
               temp1(k)=B(k,ii)
               temp2(k)=gradN(k,jj)
            enddo


            do i1=1,3 
               temp3(i1) = temp1(1)*TIJmatrix(1,i1)+
     &                     temp1(2)*TIJmatrix(2,i1)
            end do

            do i1=1,3   
               temp4(i1) = temp2(1)*IJmatrix(i1,1)+
     &                     temp2(2)*IJmatrix(i1,2)
            end do
              
            btzmtx(ii,jj)=btzmtx(ii,jj)+weigh*volume*
     &            (temp3(1)*temp4(1)+temp3(2)*temp4(2)+
     &            temp3(3)*temp4(3))*dabs(detJ)
         enddo
        enddo


       enddo

       



ccccccccccccccccccccccccccccccccc

c........fill matrix att and btt!

       do i=1,3  ! edge
         do j=1,3  ! edge
            i1=i
            i2=mod(i,3)+1
            j1=j
            j2=mod(j,3)+1

            call signcheck(wpsp(i1,kip),
     &                     wpsp(i2,kip),
     &                     wpsp(j1,kip),
     &               wpsp(j2,kip),ktest)
            call signcheck(wpsp(i1,kip),
     &                     wpsp(i2,kip),
     &                     1,2,ktest1)
            call signcheck(1,2,
     &                     wpsp(j1,kip),
     &               wpsp(j2,kip),ktest2)

            att(i,j)=(attmtx1(i,j)*(-k0*k0)*ep+
     &                attmtx2(i,j)/ur)*ktest
            att(i+3,j)=(attmtx1(i+3,j)*(-k0*k0)*ep+
     &                attmtx2(i+3,j)/ur)*ktest2
            att(i,j+3)=(attmtx1(i,j+3)*(-k0*k0)*ep+
     &                attmtx2(i,j+3)/ur)*ktest1
            att(i+3,j+3)=(attmtx1(i+3,j+3)*(-k0*k0)*ep+
     &                attmtx2(i+3,j+3)/ur)
            
            btt(i,j)=(attmtx1(i,j)/ur)*ktest
            btt(i+3,j)=(attmtx1(i+3,j)/ur)*ktest2
            btt(i,j+3)=(attmtx1(i,j+3)/ur)*ktest1
            btt(i+3,j+3)=(attmtx1(i+3,j+3)/ur)
            
         enddo
      enddo

      do i=1,3  ! edge
         do j=1,2  ! face
            i1=i
            i2=mod(i,3)+1

            call signcheck(wpsp(i1,kip),
     &                     wpsp(i2,kip),
     &                     1,2,ktest)

            att(i,j+6)=(attmtx1(i,j+6)*(-k0*k0)*ep+
     &                attmtx2(i,j+6)/ur)*ktest
            att(i+3,j+6)=(attmtx1(i+3,j+6)*(-k0*k0)*ep+
     &                attmtx2(i+3,j+6)/ur)
            att(j+6,i)=(attmtx1(j+6,i)*(-k0*k0)*ep+
     &                attmtx2(j+6,i)/ur)*ktest
            att(j+6,i+3)=(attmtx1(j+6,i+3)*(-k0*k0)*ep+
     &                attmtx2(j+6,i+3)/ur)

            btt(i,j+6)=(attmtx1(i,j+6)/ur)*ktest
            btt(i+3,j+6)=(attmtx1(i+3,j+6)/ur)
            btt(j+6,i)=(attmtx1(j+6,i)/ur)*ktest
            btt(j+6,i+3)=(attmtx1(j+6,i+3)/ur)
     
         enddo
      enddo

      do i=1,2  ! face
         do j=1,2  ! face

            att(i+6,j+6)=(attmtx1(j+6,i+6)*(-k0*k0)*ep+
     &                attmtx2(j+6,i+6)/ur)

            btt(i+6,j+6)=(attmtx1(j+6,i+6)/ur)

         enddo
      enddo

c........fill matrix bzz!

      do i=1,6  ! node&edge
         do j=1,6  ! node&edge

            bzz(i,j)=bzzmtx1(i,j)*(-k0*k0)*ep+
     &            bzzmtx2(i,j)/ur

         enddo
      enddo

c........fill matrix btz!
      do i=1,3  ! edge
         do j=1,6  ! node&edge
            i1=i
            i2=mod(i,3)+1

            call signcheck(wpsp(i1,kip),
     &                     wpsp(i2,kip),
     &                     1,2,ktest)

            btz(i,j)=btzmtx(i,j)/ur*ktest

            btz(i+3,j)=btzmtx(i+3,j)/ur

         enddo
      enddo

      do i=1,2  ! face
         do j=1,6  ! node&edge

            btz(i+6,j)=btzmtx(i+6,j)/ur

         enddo
      enddo

c........fill matrix bzt!
      do i=1,6  ! node&edge
         do j=1,8  ! edge&edge&face

            bzt(i,j)=btz(j,i)

         enddo
      enddo
      i=0
      if(i.eq.1)then      
      print*,'att(',kip,'):'
      do i=1,3
      print*,att(i,1),att(i,2),att(i,3)
      enddo
      print*,''

      print*,'btt(',kip,'):'
      do i=1,3
      print*,btt(i,1),btt(i,2),btt(i,3)
      enddo
      print*,''

      print*,'bzz(',kip,'):'
      do i=1,3
      print*,bzz(i,1),bzz(i,2),bzz(i,3)
      enddo
      print*,''

      print*,'btz(',kip,'):'
      do i=1,3
      print*,btz(i,1),btz(i,2),btz(i,3)
      enddo
      print*,''

      print*,'attmtx1(',kip,'):'
      do i=1,3
      print*,attmtx1(i,1),attmtx1(i,2),attmtx1(i,3)
      enddo
      print*,''
      print*,'attmtx2(',kip,'):'
      do i=1,3
      print*,attmtx2(i,1),attmtx2(i,2),attmtx2(i,3)
      enddo
      print*,''
      print*,'detJ:',detJ
      print*,''
      print*,'cccccccccccccccccccccccccccccccccccccc'
      print*,''
      endif

     



      return
      end

