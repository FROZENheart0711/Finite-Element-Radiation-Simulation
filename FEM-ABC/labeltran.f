      subroutine LABELTRAN(i,j,mpresent)

      IMPLICIT NONE

c.....Input Data
      INTEGER i, j

c.....Output Data

      INTEGER mpresent

      if (i.eq.1) then
         mpresent=j-1
      else
      endif
      if (i.eq.2) then
         mpresent=j+1
      else
      endif
      if (i.eq.3) then
         mpresent=6
      else
      endif

      RETURN
      END
