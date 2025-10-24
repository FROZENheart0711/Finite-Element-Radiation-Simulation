      subroutine areapatch(xx,yy,zz,length,area)
      implicit none
      real area,length(3)
      real xx(3),yy(3),zz(3)
      real vetac(3),vetbc(3),vett(3)
      integer i,k1,k2
      length(1)=sqrt((xx(1)-xx(2))**2+(yy(1)-yy(2))**2
     &            +(zz(1)-zz(2))**2)
      length(2)=sqrt((xx(2)-xx(3))**2+(yy(2)-yy(3))**2
     &            +(zz(2)-zz(3))**2)
      length(3)=sqrt((xx(1)-xx(3))**2+(yy(1)-yy(3))**2
     &            +(zz(1)-zz(3))**2)
      vetac(1)=xx(1)-xx(3)
      vetac(2)=yy(1)-yy(3)
      vetac(3)=zz(1)-zz(3)

      vetbc(1)=xx(2)-xx(3)
      vetbc(2)=yy(2)-yy(3)
      vetbc(3)=zz(2)-zz(3)
      call xmul(vetac,vetbc,vett)
      area=sqrt(vett(1)**2+vett(2)**2+vett(3)**2)
      return
      end

      subroutine dareapatch(xx,yy,zz,length,area)
      implicit none
      real*8 area,length(3)
      real*8 xx(3),yy(3),zz(3)
      real*8 vetac(3),vetbc(3),vett(3)
      integer i,k1,k2
      length(1)=sqrt((xx(1)-xx(2))**2+(yy(1)-yy(2))**2
     &            +(zz(1)-zz(2))**2)
      length(2)=sqrt((xx(2)-xx(3))**2+(yy(2)-yy(3))**2
     &            +(zz(2)-zz(3))**2)
      length(3)=sqrt((xx(1)-xx(3))**2+(yy(1)-yy(3))**2
     &            +(zz(1)-zz(3))**2)
      vetac(1)=xx(1)-xx(3)
      vetac(2)=yy(1)-yy(3)
      vetac(3)=zz(1)-zz(3)

      vetbc(1)=xx(2)-xx(3)
      vetbc(2)=yy(2)-yy(3)
      vetbc(3)=zz(2)-zz(3)
c      print*,"vet",vetbc(1),vetbc(2),vetbc(3)
      vett(1)=vetac(2)*vetbc(3)-vetac(3)*vetbc(2)
c      print*,"vett1",vett(1)
      vett(2)=vetac(3)*vetbc(1)-vetac(1)*vetbc(3)
      vett(3)=vetac(1)*vetbc(2)-vetac(2)*vetbc(1)
c      call dxmul(vetac,vetbc,vett)
c      print*,"vett",vett(1),vett(2),vett(3)
      area=sqrt(vett(1)**2+vett(2)**2+vett(3)**2)
c      print*,"area",area
      return
      end


