        subroutine SIGNCHECK(m1,m2,n1,n2,ktest)

        IMPLICIT NONE

c.......Input Data

        INTEGER m1,m2,n1,n2

c.......Output Data

        INTEGER ktest

c.......Working Variables

        INTEGER itest1, itest2

         if ((m1-m2).gt.0) then
           itest1=1
         else
           itest1=-1
         endif
         if ((n1-n2).gt.0) then
           itest2=1
         else
           itest2=-1
         endif

         ktest=itest1*itest2

        RETURN
        END
