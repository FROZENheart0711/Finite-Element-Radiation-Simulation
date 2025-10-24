      subroutine gausspoints(maxrulef,maxgridf,ngridf,
     &                       vtf1,vtf2,vtf3,wtf)
      implicit none
      integer maxrulef,maxgridf
      integer ngridf(maxrulef)
      real    vtf1(maxgridf,maxrulef),vtf2(maxgridf,maxrulef)
      real    vtf3(maxgridf,maxrulef),wtf(maxgridf,maxrulef)
      integer i
      do i=1,maxrulef
           if(i.eq.1)then
              ngridf(i)=1
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =1.0
           else if(i.eq.2)then
              ngridf(i)=3
              vtf1(1,i)=0.5
              vtf2(1,i)=0.5
              vtf3(1,i)=0.0
              wtf(1,i) =0.33333333
 
              vtf1(2,i)=0.0
              vtf2(2,i)=0.5
              vtf3(2,i)=0.5
              wtf(2,i) =0.33333333
             
              vtf1(3,i)=0.5
              vtf2(3,i)=0.0
              vtf3(3,i)=0.5
              wtf(3,i) =0.33333333
          else if(i.eq.3)then
              ngridf(i)=4
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =-0.5625000
 
              vtf1(2,i)=0.73333333
              vtf2(2,i)=0.13333333
              vtf3(2,i)=0.13333333
              wtf(2,i) =0.52083333

              vtf1(3,i)=0.13333333
              vtf2(3,i)=0.73333333
              vtf3(3,i)=0.13333333
              wtf(3,i) =0.52083333
          
              vtf1(4,i)=0.13333333
              vtf2(4,i)=0.13333333
              vtf3(4,i)=0.73333333
              wtf(4,i) =0.52083333
          else if(i.eq.4)then
              ngridf(i)=7
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =0.45000000
 
              vtf1(2,i)=0.5
              vtf2(2,i)=0.5
              vtf3(2,i)=0.0
              wtf(2,i) =0.13333333

              vtf1(3,i)=0.0
              vtf2(3,i)=0.5
              vtf3(3,i)=0.5
              wtf(3,i) =0.13333333
              
              vtf1(4,i)=0.5
              vtf2(4,i)=0.0
              vtf3(4,i)=0.5
              wtf(4,i) =0.13333333

              vtf1(5,i)=1.0
              vtf2(5,i)=0.0
              vtf3(5,i)=0.0
              wtf(5,i) =0.05000000
             
              vtf1(6,i)=0.0
              vtf2(6,i)=1.0
              vtf3(6,i)=0.0
              wtf(6,i) =0.05000000
         
              vtf1(7,i)=0.0
              vtf2(7,i)=0.0
              vtf3(7,i)=1.0
              wtf(7,i) =0.05000000
          else if(i.eq.5)then
              ngridf(i)=7
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =0.22500000
              
 
              vtf1(2,i)=0.05961587
              vtf2(2,i)=0.47014206
              vtf3(2,i)=0.47014206
              wtf(2,i) =0.13239415

              vtf1(3,i)=0.47014206
              vtf2(3,i)=0.05961587
              vtf3(3,i)=0.47014206
              wtf(3,i) =0.13239415
              
              vtf1(4,i)=0.47014206
              vtf2(4,i)=0.47014206
              vtf3(4,i)=0.05961587
              wtf(4,i) =0.13239415

              vtf1(5,i)=0.79742699
              vtf2(5,i)=0.10128651
              vtf3(5,i)=0.10128651
              wtf(5,i) =0.12593918
             
              vtf1(6,i)=0.10128651
              vtf2(6,i)=0.79742699
              vtf3(6,i)=0.10128651
              wtf(6,i) =0.12593918
         
              vtf1(7,i)=0.10128651
              vtf2(7,i)=0.10128651
              vtf3(7,i)=0.79742699
              wtf(7,i) =0.12593918
          end if
        end do
      return 
      end

      subroutine dgausspoints(maxrulef,maxgridf,ngridf,
     &                       vtf1,vtf2,vtf3,wtf)
      implicit none
      integer maxrulef,maxgridf
      integer ngridf(maxrulef)
      real*8  vtf1(maxgridf,maxrulef),vtf2(maxgridf,maxrulef)
      real*8  vtf3(maxgridf,maxrulef),wtf(maxgridf,maxrulef)
      integer i
      do i=1,maxrulef
           if(i.eq.1)then
              ngridf(i)=1
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =1.0
           else if(i.eq.2)then
              ngridf(i)=3
              vtf1(1,i)=0.5
              vtf2(1,i)=0.5
              vtf3(1,i)=0.0
              wtf(1,i) =0.33333333
 
              vtf1(2,i)=0.0
              vtf2(2,i)=0.5
              vtf3(2,i)=0.5
              wtf(2,i) =0.33333333
             
              vtf1(3,i)=0.5
              vtf2(3,i)=0.0
              vtf3(3,i)=0.5
              wtf(3,i) =0.33333333
          else if(i.eq.3)then
              ngridf(i)=4
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =-0.5625000
 
              vtf1(2,i)=0.73333333
              vtf2(2,i)=0.13333333
              vtf3(2,i)=0.13333333
              wtf(2,i) =0.52083333

              vtf1(3,i)=0.13333333
              vtf2(3,i)=0.73333333
              vtf3(3,i)=0.13333333
              wtf(3,i) =0.52083333
          
              vtf1(4,i)=0.13333333
              vtf2(4,i)=0.13333333
              vtf3(4,i)=0.73333333
              wtf(4,i) =0.52083333
          else if(i.eq.4)then
              ngridf(i)=7
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =0.45000000
 
              vtf1(2,i)=0.5
              vtf2(2,i)=0.5
              vtf3(2,i)=0.0
              wtf(2,i) =0.13333333

              vtf1(3,i)=0.0
              vtf2(3,i)=0.5
              vtf3(3,i)=0.5
              wtf(3,i) =0.13333333
              
              vtf1(4,i)=0.5
              vtf2(4,i)=0.0
              vtf3(4,i)=0.5
              wtf(4,i) =0.13333333

              vtf1(5,i)=1.0
              vtf2(5,i)=0.0
              vtf3(5,i)=0.0
              wtf(5,i) =0.05000000
             
              vtf1(6,i)=0.0
              vtf2(6,i)=1.0
              vtf3(6,i)=0.0
              wtf(6,i) =0.05000000
         
              vtf1(7,i)=0.0
              vtf2(7,i)=0.0
              vtf3(7,i)=1.0
              wtf(7,i) =0.05000000
          else if(i.eq.5)then
              ngridf(i)=7
              vtf1(1,i)=0.33333333
              vtf2(1,i)=0.33333333
              vtf3(1,i)=0.33333333
              wtf(1,i) =0.22500000
              
 
              vtf1(2,i)=0.05961587
              vtf2(2,i)=0.47014206
              vtf3(2,i)=0.47014206
              wtf(2,i) =0.13239415

              vtf1(3,i)=0.47014206
              vtf2(3,i)=0.05961587
              vtf3(3,i)=0.47014206
              wtf(3,i) =0.13239415
              
              vtf1(4,i)=0.47014206
              vtf2(4,i)=0.47014206
              vtf3(4,i)=0.05961587
              wtf(4,i) =0.13239415

              vtf1(5,i)=0.79742699
              vtf2(5,i)=0.10128651
              vtf3(5,i)=0.10128651
              wtf(5,i) =0.12593918
             
              vtf1(6,i)=0.10128651
              vtf2(6,i)=0.79742699
              vtf3(6,i)=0.10128651
              wtf(6,i) =0.12593918
         
              vtf1(7,i)=0.10128651
              vtf2(7,i)=0.10128651
              vtf3(7,i)=0.79742699
              wtf(7,i) =0.12593918
          end if
        end do
      return 
      end 
