      Subroutine Operate_four(tmp1,tmp2,tmp3,
     &                        A1,A2,A3,A4)
     
      implicit none
      integer A1,A2,A3,A4
      integer tmp1,tmp2,tmp3,tmp4

      A1=tmp1+tmp2+tmp3;
      A2=tmp1*tmp2*tmp3;

      A3=tmp1*tmp1+tmp2*tmp2+tmp3*tmp3;
      A4=tmp1*tmp1*tmp1+tmp2*tmp2*tmp2+tmp3*tmp3*tmp3;
      return
      end
