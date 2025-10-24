      subroutine sur_vol_edge
     &( kec,
     &  ld1, ip1, ld2, ip2, ip,
     &  rp1, iqof1, iqfo1, rp2, iqof2, iqfo2, 
     &  xyz1, xyz2)
      implicit none
      integer ld1(2, *), ip1(*),
     &        ld2(2, *), ip2(*), ip(*),
     &        iqof1(*), iqfo1(*), iqof2(*), iqfo2(*)

       real   xyz1(3, *), xyz2(3, *), rp1(*), rp2(*)
       integer kec,i,k1,k2,j,ij,l1,l2
       real d1,d2
       do i=1, kec
          k1=ld1(1,i)
          k2=ld1(2,i)
          call dis2o(xyz1(1,k1),xyz1(1,k2),rp1(i))
          k1=ld2(1,i)
          k2=ld2(2,i)
          call dis2o(xyz2(1,k1),xyz2(1,k2),rp2(i))
       enddo
  
          call hpsort(kec, rp1, iqof1, iqfo1)
          call hpsort(kec, rp2, iqof2, iqfo2)
       ij=0
       j=1
       do 10 i=1, kec
          j=j-ij
          ij=0
 20       if (rp1(i).eq.rp2(j)) then
             k1=ld1(1,iqfo1(i))
             k2=ld1(2,iqfo1(i))
             l1=ld2(1,iqfo2(j))
             l2=ld2(2,iqfo2(j))
             call dis(xyz1(1,k1),xyz1(1,k2),xyz2(1,l1),xyz2(1,l2), d1)
             call dis(xyz1(1,k1),xyz1(1,k2),xyz2(1,l2),xyz2(1,l1), d2)
             if ((d1.le.0.00001).or.(d2.le.0.00001)) then
                ip(iqfo1(i))=ip2(iqfo2(j))
             else
                ij=ij+1
                j=j+1
                if((j.gt.kec).or.(ij.gt.kec)) goto 10
                goto 20
             endif
          else
             j=j+1
             if(j.gt.kec) goto 10
             goto 20
          endif
 10    continue

       return
       end
