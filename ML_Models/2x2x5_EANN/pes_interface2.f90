subroutine pes_nCO_Pt6x6(natom,x,e,dedx)
  implicit none
  integer,intent(in) :: natom
  real*8,intent(in) :: x(3,natom)
  real*8,intent(out) :: e,dedx(3,natom)

  integer :: nco
  real*8 :: c(3,2),v,dv(3,2)

  real*8 :: x0(3,natom),xn(3,natom)
  real*8 :: ococ(3,4),vint,dococ(3,4)

  integer :: i,j,k,m1,m2

  real*8,parameter :: a0=2.812256750888d0
  integer,parameter :: ncell=6
  real*8,parameter :: a(3)=[a0*sqrt(3.d0)/2.d0,-a0/2.d0,0.d0]*ncell
  real*8,parameter :: b(3)=[0.d0,a0,0.d0]*ncell

  e=0.d0
  dedx=0.d0
  nco=natom/2
  !! CO - Pt 6D
  do i=1,nco
     c=x(1:3,(i*2-1):(i*2))
     call co_pt_pes(c,v,1,dv) !! Angstrom, eV
     e=e+v
     dedx(1:3,(i*2-1):(i*2))=-dv
  enddo
  !! to calculate co-co interactions, move all CO to the 6x6 cell first
  call move_nco_to_6x6_cell(natom,x,x0)

  !! CO-CO interaction at the original supercell
  do i=1,nco-1
  do j=i+1,nco
     ococ=x0(1:3,[i*2,i*2-1,j*2,j*2-1])
     call energy1(ococ,vint,dococ)
     e=e+vint
     dedx(1:3,[i*2,i*2-1,j*2,j*2-1])=dedx(1:3,[i*2,i*2-1,j*2,j*2-1])-dococ
  enddo
  enddo

  !! CO-CO interaction with the neighbouring supercell
  do m1=-1,0
  do m2=-1,1
     if(m1.eq.0 .and. m2.ne.1) cycle !! only 4-pairs are required
     do k=1,natom
        xn(1:3,k)=x0(1:3,k)+m1*a+m2*b
     enddo
     do i=1,nco
     do j=1,nco
        ococ(1:3,1:2)=x0(1:3,[i*2,i*2-1])
        ococ(1:3,3:4)=xn(1:3,[j*2,j*2-1])
        call energy1(ococ,vint,dococ)
        e=e+vint
       !print*,dsqrt(dot_product(ococ(1:3,2)-ococ(1:3,4),ococ(1:3,2)-ococ(1:3,4))),vint
        dedx(1:3,[i*2,i*2-1])=dedx(1:3,[i*2,i*2-1])-dococ(:,1:2)
        dedx(1:3,[j*2,j*2-1])=dedx(1:3,[j*2,j*2-1])-dococ(:,3:4)
     enddo
     enddo
  enddo
  enddo

  return
end subroutine

subroutine move_nCO_to_6x6_cell(natom,x,x0)
  implicit none
  real*8,parameter :: a0=2.812256750888d0
  integer,parameter :: ncell=6
  real*8,parameter :: a(3)=[a0*sqrt(3.d0)/2.d0,-a0/2.d0,0.d0]*ncell
  real*8,parameter :: b(3)=[0.d0,a0,0.d0]*ncell
  integer,intent(in) :: natom
  real*8,intent(in) :: x(3,natom)
  real*8,intent(out) :: x0(3,natom)
  integer :: i,j,k,nco
  x0=x

  nco=natom/2
  do i=1,nco
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
