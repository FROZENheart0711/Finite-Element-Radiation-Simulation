      subroutine nodeface
     &( i, kf1, kf2)

      ie=(i-1)/6 + 1
      in=i-(ie-1)*6

      if (in.eq.1) then
          kf1=1
          kf2=2
      else
      endif
      
      if (in.eq.2) then
          kf1=1
          kf2=3
      else
      endif

      if (in.eq.3) then
          kf1=2
          kf2=3
      else
      endif

      if (in.eq.4) then
          kf1=1
          kf2=4
      else
      endif

      if (in.eq.5) then
          kf1=2
          kf2=4
      else
      endif

      if (in.eq.6) then
          kf1=3
          kf2=4
      else
      endif

      return
      end
 
