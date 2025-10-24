       subroutine assemsrc(kip,att,btt,btz,bzt,bzz,att_f,btt_1,
     &          btz_1,bzt_1,bzz_1,nsrc,nsedge,nsnode,p_mvolnodefrt,
     &          wpse,wpsp,p_srcspfrt,mvoledge,mvolnode,ipoe,ipon,
     &          ebnd,nbnd)
        
        
        implicit none

!.......Input Data
        integer kip
        complex*16 att(8,8), btt(8,8)
        complex*16 btz(8,6), bzt(6,8), bzz(6,6)
        integer nsrc
        integer nsedge,nsnode
        integer wpse(3,*)
        integer wpsp(4,*)
        integer*8 p_srcspfrt
        integer*8 p_mvolnodefrt
        integer mvolnode, mvoledge
        integer ipoe(mvoledge),ipon(mvolnode)
        integer ebnd(nsedge),nbnd(nsnode)

!.......Output Data
        complex*16 att_f(nsedge*2+nsrc*2,nsedge*2+nsrc*2)
        complex*16 btt_1(nsedge*2+nsrc*2,nsedge*2+nsrc*2)
        complex*16 btz_1(nsedge*2+nsrc*2,nsnode+nsedge)
        complex*16 bzt_1(nsnode+nsedge,nsedge*2+nsrc*2)
        complex*16 bzz_1(nsnode+nsedge,nsnode+nsedge)

!.......Working Variables
        integer i,j,k1,k2


c.......begin!

c........assem Att!
        do i=1,3  ! edge
            do j=1,3  ! edge
                k1=ipoe(wpse(i,kip))
                k2=ipoe(wpse(j,kip))

                if((ebnd(k1).le.0).and.(ebnd(k2).le.0))then
                    att_f(k1,k2)=att_f(k1,k2)+att(i,j)
                    att_f(k1,k2+nsedge)=att_f(k1,k2+nsedge)
     &                                  +att(i,j+3)
                    att_f(k1+nsedge,k2)=att_f(k1+nsedge,k2)
     &                                  +att(i+3,j)
                    att_f(k1+nsedge,k2+nsedge)=att_f(k1+nsedge,
     &                                  k2+nsedge)+att(i+3,j+3)
                elseif(k1.eq.k2)then
                    att_f(k1,k2)=1.
                    att_f(k1+nsedge,k2+nsedge)=1.                 
                endif
            enddo
        enddo

        do i=1,3  ! edge
            do j=1,2  ! face
                k1=ipoe(wpse(i,kip))
                k2=kip+nsrc*(j-1)+nsedge*2

                if((ebnd(k1).le.0))then
                    att_f(k1,k2)=att_f(k1,k2)+att(i,6+j)

                    att_f(k1+nsedge,k2)=att_f(k1+nsedge,k2)+att(i+3,6+j)
                    
                    att_f(k2,k1)=att_f(k2,k1)+att(6+j,i)
     
                    att_f(k2,k1+nsedge)=att_f(k2,k1+nsedge)+att(6+j,i+3)

                endif
            enddo
        enddo

        do i=1,2  ! face
            do j=1,2  ! face
                k1=kip+nsrc*(i-1)+nsedge*2
                k2=kip+nsrc*(j-1)+nsedge*2

                att_f(k1,k2)=att_f(k1,k2)+att(i+6,j+6)
            enddo
        enddo

c........assem Btt!
        do i=1,3  ! edge
            do j=1,3  ! edge
                k1=ipoe(wpse(i,kip))
                k2=ipoe(wpse(j,kip))

                if((ebnd(k1).le.0).and.(ebnd(k2).le.0))then
                    btt_1(k1,k2)=btt_1(k1,k2)+btt(i,j)
                    btt_1(k1,k2+nsedge)=btt_1(k1,k2+nsedge)
     &                                  +btt(i,j+3)
                    btt_1(k1+nsedge,k2)=btt_1(k1+nsedge,k2)
     &                                  +btt(i+3,j)
                    btt_1(k1+nsedge,k2+nsedge)=btt_1(k1+nsedge,
     &                                  k2+nsedge)+btt(i+3,j+3)
                elseif(k1.eq.k2)then
                    btt_1(k1,k2)=1.
                    btt_1(k1+nsedge,k2+nsedge)=1.                 
                endif
            enddo
        enddo

        do i=1,3  ! edge
            do j=1,2  ! face
                k1=ipoe(wpse(i,kip))
                k2=kip+nsrc*(j-1)+nsedge*2

                if((ebnd(k1).le.0))then
                    btt_1(k1,k2)=btt_1(k1,k2)+btt(i,6+j)

                    btt_1(k1+nsedge,k2)=btt_1(k1+nsedge,k2)+btt(i+3,6+j)
                    
                    btt_1(k2,k1)=btt_1(k2,k1)+btt(6+j,i)
     
                    btt_1(k2,k1+nsedge)=btt_1(k2,k1+nsedge)+btt(6+j,i+3)

                endif
            enddo
        enddo

        do i=1,2  ! face
            do j=1,2  ! face
                k1=kip+nsrc*(i-1)+nsedge*2
                k2=kip+nsrc*(j-1)+nsedge*2

                btt_1(k1,k2)=btt_1(k1,k2)+btt(i+6,j+6)
            enddo
        enddo

c........assem Bzz!
        do i=1,3  ! node
            do j=1,3  ! node
                k1=ipon(wpsp(i,kip))
                k2=ipon(wpsp(j,kip))

                if((nbnd(k1).le.0).and.(nbnd(k2).le.0))then
                    bzz_1(k1,k2)=bzz_1(k1,k2)+bzz(i,j)
                elseif(k1.eq.k2)then
                    bzz_1(k1,k2)=1.
                endif
            enddo
        enddo

        do i=1,3  ! node
            do j=1,3  ! edge
                k1=ipon(wpsp(i,kip))
                k2=ipoe(wpse(j,kip))+nsnode

                if((nbnd(k1).le.0).and.(ebnd(k2-nsnode).le.0))then
                    bzz_1(k1,k2)=bzz_1(k1,k2)+bzz(i,j+3)
                    bzz_1(k2,k1)=bzz_1(k2,k1)+bzz(j+3,i)
                endif
            enddo
        enddo

        do i=1,3  ! edge
            do j=1,3  ! edge
                k1=ipoe(wpse(i,kip))
                k2=ipoe(wpse(j,kip))

                if((ebnd(k1).le.0).and.(ebnd(k2).le.0))then
                    bzz_1(k1+nsnode,k2+nsnode)=
     &                      bzz_1(k1+nsnode,k2+nsnode)+bzz(i+3,j+3)
                elseif(k1.eq.k2)then
                    bzz_1(k1+nsnode,k2+nsnode)=1.
                endif
            enddo
        enddo

c........assem Btz!
        do i=1,3  ! edge
            do j=1,3  ! node
                k1=ipoe(wpse(i,kip))
                k2=ipon(wpsp(j,kip))

                if((ebnd(k1).le.0).and.(nbnd(k2).le.0))then
                    btz_1(k1,k2)=btz_1(k1,k2)+btz(i,j)
                    btz_1(k1+nsedge,k2)=btz_1(k1+nsedge,k2)+btz(i+3,j)
                endif   
            enddo
        enddo

        do i=1,3  ! edge
            do j=1,3  ! edge
                k1=ipoe(wpse(i,kip))
                k2=ipoe(wpse(j,kip))

                if((ebnd(k1).le.0).and.(ebnd(k2).le.0))then
                    btz_1(k1,k2+nsnode)=btz_1(k1,k2+nsnode)+btz(i,j+3)
                    btz_1(k1+nsedge,k2+nsnode)=
     &                      btz_1(k1+nsedge,k2+nsnode)+btz(i+3,j+3)
                endif   
            enddo
        enddo

        do i=1,2  ! face
            do j=1,3  ! node
                k1=kip+nsrc*(i-1)+nsedge*2
                k2=ipon(wpsp(j,kip))

                if(nbnd(k2).le.0)then
                    btz_1(k1,k2)=btz_1(k1,k2)+btz(i+6,j)
                endif
            enddo
        enddo

        do i=1,2  ! face
            do j=1,3  ! edge
                k1=kip+nsrc*(i-1)+nsedge*2
                k2=ipoe(wpse(j,kip))

                if(ebnd(k2).le.0)then
                    btz_1(k1,k2+nsnode)=btz_1(k1,k2+nsnode)+btz(i+6,j+3)  
                endif
            enddo
        enddo

c........assem Bzt!
        do i=1,nsnode+nsedge
            do j=1,nsedge*2+nsrc*2
                bzt_1(i,j)=btz_1(j,i)
            enddo
        enddo

        return
        end
