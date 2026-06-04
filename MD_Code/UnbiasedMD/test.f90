program main
 !call test_6x6_top
  call test_6x6_hollow
end

subroutine test_6x6_hollow
  implicit none
  real*8,parameter :: a0=2.812256750888d0
  real*8,parameter :: a(3)=[a0*sqrt(3.d0)/2.d0,-a0/2.d0,0.d0]
  real*8,parameter :: b(3)=[0.d0,a0,0.d0]
  real*8,parameter :: c(3)=[0.d0,0.d0,29.184792087162d0]
  integer :: i,j,k,m
  real*8 :: cart(3,72),v,dvdc(3,72),rd(3)

  call random_seed()
  m=0
  do i=1,6
  do j=1,6
     do k=1,2
        m=m+1
        call random_number(rd); rd=rd*0.15d0
        cart(1,m)=rd(1)+(i-1)*a(1)             + a0/2.d0/sqrt(3.d0)
        cart(2,m)=rd(2)+(i-1)*a(2)+(j-1)*b(2)  + a0/2.d0
        cart(3,m)=rd(3)+(k-1)*1.402d0+7.85d0
     enddo
  enddo
  enddo

  do m=1,72
     if(mod(m,2).eq.1) write(*,'(a12,3f10.5)') " C  12.00000",cart(:,m)
     if(mod(m,2).eq.0) write(*,'(a12,3f10.5)') " O  15.99491",cart(:,m)
  enddo

  call pes_co36_pt6x6(cart,v,dvdc)

  print*,v
end


subroutine test_6x6_top
  implicit none
  real*8,parameter :: a0=2.812256750888d0
  real*8,parameter :: a(3)=[a0*sqrt(3.d0)/2.d0,-a0/2.d0,0.d0]
  real*8,parameter :: b(3)=[0.d0,a0,0.d0]
  real*8,parameter :: c(3)=[0.d0,0.d0,29.184792087162d0]
  integer :: i,j,k,m
  real*8 :: cart(3,72),v,dvdc(3,72),rd(3)

  call random_seed()
  m=0
  do i=1,6
  do j=1,6
     do k=1,2
        m=m+1
        call random_number(rd); rd=rd*0.15d0
        cart(1,m)=rd(1)+(i-1)*a(1)
        cart(2,m)=rd(2)+(i-1)*a(2)+(j-1)*b(2)
        cart(3,m)=rd(3)+(k-1)*1.402d0+1.85d0
     enddo
  enddo
  enddo

  do m=1,72
     if(mod(m,2).eq.1) write(*,'(a12,3f10.5)') " C  12.00000",cart(:,m)
     if(mod(m,2).eq.0) write(*,'(a12,3f10.5)') " O  15.99491",cart(:,m)
  enddo

  call pes_co36_pt6x6(cart,v,dvdc)

  print*,v
end

subroutine test_disorption_top_1d
  implicit none
  real*8 :: c(3,2),v,dv(3,2)
  integer :: i,j,k

  ! minima on top: z(C-surf)=1.85 Angstrom
  do i=14,100
     c=0.d0
     c(3,1)=i*0.1d0
     c(3,2)=c(3,1)+1.402d0
     call co_pt_pes(c,v,1,dv)
     write(*,*) c(3,1),v
  enddo
end subroutine

