program display
  implicit none
  integer :: natom
  integer*8 :: istep
  real*8 :: pot
  real*8,allocatable :: cart(:,:),cart0(:,:)
  character(kind=1,len=1) :: atom

  real*8,parameter :: a0=2.812256750888d0
  real*8,parameter :: a(3)=[a0*sqrt(3.d0)/2.d0,-a0/2.d0,0.d0]
  real*8,parameter :: b(3)=[0.d0,a0,0.d0]
  real*8,parameter :: c(3)=[0.d0,0.d0,29.184792087162d0]
  integer :: i,j,k,m

  real*8 :: surf(3,49)
  integer :: iframe

  !! define the surface atoms (1-layer)
  m=0
  do i=1,7
  do j=1,7
     m=m+1
     surf(1,m)=(i-1)*a(1)
     surf(2,m)=(i-1)*a(2)+(j-1)*b(2)
     surf(3,m)=0.d0
  enddo
  enddo

  !! check natom
  open(991,file='x',status='old',action='read')
  read(991,*) natom
  allocate(cart(3,natom))
  allocate(cart0(3,natom))
  rewind(991)

  !! read and write
  iframe=0
  do while(.true.)
     read(991,*,end=98) natom
     read(991,*,end=98) istep,pot
     do m=1,natom
        read(991,*,end=98) atom,cart(1:3,m)
     enddo

     call move_nCO_to_6x6_cell(natom,cart,cart0)

     iframe=iframe+1

     if (mod(iframe,100).ne.1) cycle

     write(995,*) natom+49
     write(995,'(i12,f16.5)') istep,pot
     do m=1,natom
        if(mod(m,2).eq.1) write(995,'(a3,3f10.5)') " C ",cart0(:,m)
        if(mod(m,2).eq.0) write(995,'(a3,3f10.5)') " O ",cart0(:,m)
     enddo
     do m=1,49
        write(995,'(a3,3f10.5)') " Pt",surf(:,m)
     enddo

  enddo
  98 close(991)
  deallocate(cart)
  deallocate(cart0)
end

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

