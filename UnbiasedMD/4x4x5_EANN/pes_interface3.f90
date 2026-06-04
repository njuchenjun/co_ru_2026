module pes_interface3

contains

subroutine pes_nCO_RuHD(natom,x,e,dedx) ! angstrom , ev

  implicit none
  integer,intent(in) :: natom
  integer :: nco
  integer,parameter :: nsurf=32
  real*8,intent(in) :: x(3,natom)
  real*8,intent(out) :: e,dedx(3,natom)

  real*8 :: csurf(3,nsurf),vsurf,dvsurf(3,nsurf)

!   real*8 :: cvdw(3,nsurf+2),vvdw,dvvdw(3,nsurf+2)
  real*8 :: cpbe(3,nsurf+2),vpbe,dvpbe(3,nsurf+2) !!!2025/06/16 hje

  real*8 :: ctri(3,nsurf+4),vtri,dvtri(3,nsurf+4) !!!2025/08/18 hje intents for tribody potentials

  real*8 :: x0(3,natom),xn(3,natom)

  !! for co-co interactions
  real*8 :: ococ(3,4),vint,dococ(3,4)

  integer :: i,j,k,m1,m2

  real*8,parameter :: a0=2.740296191462d0
  integer,parameter :: ncell=4
  real*8,parameter :: a(3)=[a0*sqrt(3.d0)/2.d0,-a0/2.d0,0.d0]*ncell
  real*8,parameter :: b(3)=[0.d0,a0,0.d0]*ncell

  e=0.d0
  dedx=0.d0

  nco=(natom-nsurf)/2
  if(nco.le.0) stop "wrong INPUT!"

  !! Ru surface

  csurf=x(1:3,1:nsurf)
  vsurf=0.d0; dvsurf=0.d0
  call Ru_surf_pes(nsurf,csurf,vsurf,dvsurf)
  e=e+vsurf
  dedx(1:3,1:nsurf)=dvsurf

 !print*,vsurf,e

  !! CO - Ru HD interaction (together with CO 1D high-level PES)
!!!2025/06/16 hje for Ru-CO interactions
   do i=1,nco
     cpbe(1:3,1:nsurf)=csurf
     cpbe(1:3,(nsurf+1):(nsurf+2))=x(1:3,(nsurf+i*2-1):(nsurf+i*2))
     vpbe=0.d0; dvpbe=0.d0
     call CORu_int_pes(nsurf+2,cpbe,vpbe,dvpbe) ! contains CO 1D PES
     e=e+vpbe
    !print*,vpbe,e
     dedx(1:3,1:nsurf)=dedx(1:3,1:nsurf)+dvpbe(1:3,1:nsurf)
     dedx(1:3,(nsurf+i*2-1):(nsurf+i*2))=dedx(1:3,(nsurf+i*2-1):(nsurf+i*2))+dvpbe(1:3,nsurf+1:nsurf+2)
  enddo

  !! to calculate co-co interactions, move all CO to the 4x4 cell first
  call move_nCOatRu_to_4x4_cell(natom,nsurf,x,x0)
 !print*,'here'
  !! CO-CO interaction at the original supercell 
  !! And CO-CO-surface interation at the original supercell 
  do i=1,nco-1
  do j=i+1,nco
     ococ=x0(1:3,[nsurf+i*2,nsurf+i*2-1,nsurf+j*2,nsurf+j*2-1])
     !ctri(1:3,1:nsurf) = csurf !!! coordination for tribody system
     ctri(1:3,1:nsurf) = x0(1:3,1:nsurf) !!! coordination for tribody system
     ctri(1:3,(nsurf+1):(nsurf+4)) = x0(1:3,[nsurf+i*2-1,nsurf+i*2,nsurf+j*2-1,nsurf+j*2])
     vint=0.d0; dococ=0.d0
     vtri=0.d0; dvtri=0.d0 !!! initialize tribody potentials 
     call energy1(ococ,vint,dococ)
     call tribody_int_pes(nsurf+4,ctri,vtri,dvtri)
     e=e+vint
     e=e+vtri
    !print*,vint,e
     dedx(1:3,[nsurf+i*2,nsurf+i*2-1,nsurf+j*2,nsurf+j*2-1])=dedx(1:3,[nsurf+i*2,nsurf+i*2-1,nsurf+j*2,nsurf+j*2-1])-dococ
     dedx(1:3,1:nsurf)=dedx(1:3,1:nsurf)+dvtri(1:3,1:nsurf)
     dedx(1:3,[nsurf+i*2-1,nsurf+i*2,nsurf+j*2-1,nsurf+j*2])=dedx(1:3,[nsurf+i*2-1,nsurf+i*2,nsurf+j*2-1,nsurf+j*2])+dvtri(1:3,(nsurf+1):(nsurf+4))
  enddo
  enddo

 !print*,'here here'
  !! CO-CO interaction with the neighbouring supercell
  do m1=-1,0
  do m2=-1,1
     if(m1.eq.0 .and. m2.ne.1) cycle !! only 4-pairs are required
     do k=1,natom
        xn(1:3,k)=x0(1:3,k)+m1*a+m2*b
     enddo
     do i=1,nco
     do j=1,nco
        ococ(1:3,1:2)=x0(1:3,[nsurf+i*2,nsurf+i*2-1])
        ococ(1:3,3:4)=xn(1:3,[nsurf+j*2,nsurf+j*2-1])
        vint=0.d0; dococ=0.d0
        call energy1(ococ,vint,dococ)
        e=e+vint
       !print*,vint,vtri,e,m1,m2,i,j
       !print*,dsqrt(dot_product(ococ(1:3,2)-ococ(1:3,4),ococ(1:3,2)-ococ(1:3,4))),vint
        dedx(1:3,[nsurf+i*2,nsurf+i*2-1])=dedx(1:3,[nsurf+i*2,nsurf+i*2-1])-dococ(:,1:2)
        dedx(1:3,[nsurf+j*2,nsurf+j*2-1])=dedx(1:3,[nsurf+j*2,nsurf+j*2-1])-dococ(:,3:4)
     enddo
     enddo
  enddo
  enddo

  !! CO-CO-SURF interaction

  return
end subroutine

subroutine move_nCOatRu_to_4x4_cell(natom,nsurf,x,x0)
  implicit none
  real*8,parameter :: a0=2.740296191462d0
  integer,parameter :: ncell=4
  real*8,parameter :: a(3)=[a0*sqrt(3.d0)/2.d0,-a0/2.d0,0.d0]*ncell
  real*8,parameter :: b(3)=[0.d0,a0,0.d0]*ncell
  integer,intent(in) :: natom,nsurf
  real*8,intent(in) :: x(3,natom)
  real*8,intent(out) :: x0(3,natom)
  integer :: i,j,k,nco

  if(mod(natom,2).eq.1 .or. mod(nsurf,2).eq.1) stop "nsurf is an odd number!"

  nco=(natom-nsurf)/2
  x0=x

  do i=(nsurf/2+1),(natom/2)
     !! left limit: 0
     do while(x0(1,i*2-1).lt.0.0)
        x0(:,i*2-1)=x0(:,i*2-1)+a
        x0(:,i*2)=x0(:,i*2)+a
     enddo
     !! right limit: a(1)
     do while(x0(1,i*2-1).gt.a(1))
        x0(:,i*2-1)=x0(:,i*2-1)-a
        x0(:,i*2)=x0(:,i*2)-a
     enddo
     !! down limit: -x/sqrt(3)
     do while(x0(2,i*2-1).lt.-x0(1,i*2-1)/sqrt(3.d0))
        x0(:,i*2-1)=x0(:,i*2-1)+b
        x0(:,i*2)=x0(:,i*2)+b
     enddo
     !! up limit: b(2)-x/sqrt(3)
     do while(x0(2,i*2-1).gt.b(2)-x0(1,i*2-1)/sqrt(3.d0))
        x0(:,i*2-1)=x0(:,i*2-1)-b
        x0(:,i*2)=x0(:,i*2)-b
     enddo
  enddo

  return
end subroutine

end module pes_interface3
