subroutine co_pt_pes(c,v,idv,dv)
  implicit none

  real*8,intent(in)  :: c(6,1)!! Cartesian of C O on Pt(111) surface, Angstrom
  real*8,intent(out) :: v(1)  !! Potential Energy in eV
  integer,intent(in) :: idv   !! 0: energy only; 1: energy and force
  real*8,intent(out) :: dv(6,1)
  real*8 :: ctrans(6,1),dvtrans(6,1),frac(6,1),zmat(6,1),dc,cx(6,1),cy(6,1),vx(1),vy(1),dvt(6,1)
  integer :: j

  real*8 :: v1,v2,v3,v4,v5
  real*8,parameter :: vinf=-130.706d0

  !! calculate potential energy
  ctrans=c
  call cart2frac(ctrans,frac)
  call frac2unit(frac,ctrans,zmat,2,0)

  call nnfit_b1(6,1,ctrans,v1,0,dvtrans);
  v=v1-vinf

  ! idv=1: calculate analytical derivatives, but directions should be rotated back

  !! calculate derivatives numerically
  if(idv.eq.1) then
    dv=0.d0
    dc=1.d-3
    do j=1,6
      cx=c; cx(j,1)=cx(j,1)-dc
      cy=c; cy(j,1)=cy(j,1)+dc

      ctrans=cx
      call cart2frac(ctrans,frac)
      call frac2unit(frac,ctrans,zmat,2,0)
      call nnfit_b1(6,1,ctrans,vx,0,dvt); vx=vx-vinf

      ctrans=cy
      call cart2frac(ctrans,frac)
      call frac2unit(frac,ctrans,zmat,2,0)
      call nnfit_b1(6,1,ctrans,vy,0,dvt); vy=vy-vinf

      dv(j,1)=(vy(1)-vx(1))/2.d0/dc
    enddo
  endif

  return
end subroutine
! --------------------------------------------------------------------
subroutine nnfit_b1(ndim,ntot,r0,vx,idv,dv)
  implicit none
  integer,intent(in) :: ndim,ntot
  real*8,intent(in) :: r0(ndim,ntot)
  real*8,intent(out) :: vx(ntot)

  integer,parameter :: s0=6,s1=40,s2=40
  character*99,save :: wfile="pes_copt6d_b48.txt"
  real*8,save :: w1(s0,s1),b1(s1),w2(s1,s2),b2(s2),w3(s2),b3,rg(2,s0),vg(2)

  integer :: idv
  real*8 :: dv(ndim,ntot)
  integer :: i
  integer,save :: init=0
  integer,save :: fid
  if(ndim.ne.s0) stop "ndim .ne. s0"
  if (init.eq.0) then
    init=1
    fid=123456
    open(fid,file=trim(wfile),status='old',action='read')
    read(fid,*) w1,b1,w2,b2,w3,b3,rg,vg
    close(fid)
  endif
  if(ntot.ge.24) then
    call nsimx(ntot,r0,vx,idv,dv,s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
  else
    do i=1,ntot
    call nsim(r0(1,i),vx(i),idv,dv(:,i),s0,s1,s2,rg,w1,b1,w2,b2,w3,b3,vg)
    enddo
  endif
  return
end subroutine
! --------------------------------------------------------------------
subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
  ! ---- simulate a neural network with two hidden layers
  ! ----      n0-n1-n2-1
  ! ---- blas routines used in this subroutine: dgemv, ddot
  implicit none
  integer,intent(in) :: n0,n1,n2,idv
  real*8,intent(in)  :: r0(n0),rg(2,n0),vg(2)
  real*8,intent(in)  :: w1(n0,n1),b1(n1)
  real*8,intent(in)  :: w2(n1,n2),b2(n2)
  real*8,intent(in)  :: w3(n2),b3
  real*8,intent(out) :: v,dv(n0)
  integer :: i,j,k
  real*8  :: r(n0),rgg(n0),vgg,ax(n1),bx(n2)
  real*8  :: dvtm,rtmp,rt1(n1),rt2(n2)
  real*8,external :: ddot
  v=0.d0
  r=r0
  vgg=vg(2)-vg(1)
  ! mapminmax [-1,1]
  do i=1,n0
    rgg(i)=rg(2,i)-rg(1,i)
    r(i)=2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
  end do
  ! 1st layer
  rt1=b1
  call dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1)
  ax=dtanh(rt1)
  ! 2nd layer
  rt2=b2
  call dgemv('t',n1,n2,1.d0,w2,n1,ax,1,1.d0,rt2,1)
  bx=dtanh(rt2)
  ! output layer
  v=b3+ddot(n2,w3,1,bx,1)
  !reverse map
  v=vgg*(v+1.d0)/2+vg(1)
  if(idv.ne.1) return
  ! calculate first derivatives, dv(i)=dv/dr(i)
  dv=0.d0
  do i=1,n0
    do k=1,n2
      dvtm=0.d0
      do j=1,n1
        dvtm=dvtm+w2(j,k)*w1(i,j)*(1-ax(j)**2)
      enddo
      dv(i)=dv(i)+w3(k)*dvtm*(1-bx(k)**2)
    enddo
    dv(i)=dv(i)*vgg/rgg(i)
  enddo
  return
end subroutine nsim
! --------------------------------------------------------------------
subroutine nsimx(nt,r,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
  ! ---- simulate a neural network with two hidden layers
  ! ----    n0-n1(tansig)-n2(tansig)-1(pureline)
  ! ---- blas routines used in this subroutine: dgemm, dgemv
  ! ---- the syntax keeps the same with nsim, except input nt at first
  ! ---- Chen Jun, 2013/11/26
  !
  ! ---- description of parameters
  ! ---- nt  : total number of points in r,v and dv
  ! ---- r   : input of nn, r(n0,nt), r(1:n0,i) is the i-th input
  ! ---- v   : output of nn, v(nt), v(i) is the i-th input
  ! ---- idv : if idv==0, calc. v only, else calc. v and derivatives dv
  ! ---- dv  : derivatives of dv/dr, dv(n0,nt), dv(j,i)=d v(i) / d r(j,i)
  ! ---- n0  : neurons of the input layer, dimension of r
  ! ---- n1  : neurons of the first hidden layer
  ! ---- n2  : neurons of the second hidden layer
  ! ---- rg  : input ranges of the training set, rg(2,n0)
  ! ---- w1  : weights connect input layer with 1-st hidden layer
  ! ---- b1  : biases of 1-st hidden layer
  ! ---- w2  : weights connect 1-st hidden layer with 2-nd hidden layer
  ! ---- b2  : biases of 2-nd hidden layer
  ! ---- w3  : weights connect 2-nd hidden layer with output layer
  ! ---- b3  : biase of output layer
  ! ---- vg  : output range of the training set, vg(2)
  implicit none
  integer,intent(in) :: nt,idv,n0,n1,n2
  real*8,intent(in) :: r(n0,nt)
  real*8,intent(out) :: v(nt),dv(n0,nt)
  real*8,intent(in) :: w1(n0,n1),b1(n1),w2(n1,n2),b2(n2),w3(n2),b3
  real*8,intent(in) :: rg(2,n0),vg(2)
  real*8 :: x0(n0,nt),x1(nt,n1),x2(nt,n2),tmp,rgg(n0),vgg
  integer :: i,j,k,n

  v=0.d0
  do i=1,n0
    rgg(i)=rg(2,i)-rg(1,i)
  enddo
  vgg=vg(2)-vg(1)

  x0=r

  do i=1,n0
    x0(i,1:nt)=(x0(i,1:nt)-rg(1,i))/rgg(i)*2.d0-1.d0
  enddo

  do i=1,n1
   x1(1:nt,i)=b1(i)
  enddo

  ! ---- x1=trans(x0(n0,nt))*w1(n0,n1)+x1(nt,n1)
  !call gpu_dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)
  call dgemm('t','n',nt,n1,n0,1.d0,x0,n0,w1,n0,1.d0,x1,nt)

  x1=dtanh(x1)

  do i=1,n2
   x2(1:nt,i)=b2(i)
  enddo

  ! ---- x2=x1(nt,n1)*w2(n1,n2)+x2(nt,n2)
  !call gpu_dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)
  call dgemm('n','n',nt,n2,n1,1.d0,x1,nt,w2,n1,1.d0,x2,nt)

  x2=dtanh(x2)

  v(1:nt)=b3

  ! ---- v=x2(nt,n2)*w3(n2)+v(nt)
  call dgemv('n',nt,n2,1.d0,x2,nt,w3,1,1.d0,v,1)

  v=(v+1.d0)/2.d0*vgg+vg(1)

  if (idv.ne.1) return
  ! ---- calculate first derivatives
  dv=0.d0
  do n=1,nt
    do i=1,n0
      do k=1,n2
        tmp=0.d0
        do j=1,n1
          tmp=tmp+w2(j,k)*w1(i,j)*(1-x1(n,j)**2)
        enddo
        dv(i,n)=dv(i,n)+w3(k)*tmp*(1-x2(n,k)**2)
      enddo
      dv(i,n)=dv(i,n)*vgg/(rgg(i))
    enddo
  enddo
  return
end subroutine nsimx
! --------------------------------------------------------------------

subroutine frac2unit(frac0,cart,zmat,level,vasp)
  !! level 0:  move C and O atoms into a same 2x2 cell, then trans to Cartesian only.
  !! level 1:  move CO into the 1x1 cell, then trans to Cartesian.
  !! level 2:  symmetrize the Cartesian to the minimum triangle

  implicit none
  real*8,parameter :: pi=dacos(-1.d0)
  real*8,intent(in) :: frac0(3,2)
  real*8,intent(out) :: cart(3,2),zmat(6)
  integer,intent(in) :: level,vasp
  real*8 :: frac(3,2),com(3),dist2d,comc(3),tmpx,tmpy,tmpa,phi,alpha

  real*8 :: a(3),b(3),c(3),a0
  integer :: i,j

  !! must be a 2x2 super cell !!
  data a/4.870971576466d0, -2.812256750888d0,  0.000000000000d0/
  data b/0.000000000000d0,  5.624513501775d0,  0.000000000000d0/
  data c/0.000000000000d0,  0.000000000000d0, 29.184792087162d0/
  data a0/5.624513501775d0/

 !data a/9.7419431529316274d0,-5.6245135017751027d0, 0.0000000000000000d0/
 !data b/0.0000000000000000d0,11.2490270035502071d0, 0.0000000000000000d0/
 !data c/0.0000000000000000d0, 0.0000000000000000d0,29.1847920871624140d0/
 !data a0/11.2490270035502071d0/

  frac=frac0

 !!! move C and O to the same 2x2 cell, only for vasp trajectory results
 !if(vasp.eq.1) then
 !dist2d=dsqrt(dot_product(frac(1:2,1)-frac(1:2,2),frac(1:2,1)-frac(1:2,2)))
 !if(dist2d.gt.0.5d0) then
 !  if(frac(1,1)-frac(1,2).gt. 0.5d0) frac(1,2)=frac(1,2)+1.d0
 !  if(frac(1,1)-frac(1,2).lt.-0.5d0) frac(1,1)=frac(1,1)+1.d0
 !  if(frac(2,1)-frac(2,2).gt. 0.5d0) frac(2,2)=frac(2,2)+1.d0
 !  if(frac(2,1)-frac(2,2).lt.-0.5d0) frac(2,1)=frac(2,1)+1.d0
 !endif
 !endif

  com(1:3)=(frac(1:3,1)+frac(1:3,2))*0.5d0

  if(level.eq.0) goto 990 !! translate to cartesian directly

  !! move all poins to a unit cell
  !! move com(1) to [0.0,0.5)
  !! move com(2) to [0.0,0.5)
  do while(com(1).ge.0.5d0)
    com(1)=com(1)-0.5d0
    frac(1,:)=frac(1,:)-0.5d0
  enddo
  do while(com(1).lt.0.0d0)
    com(1)=com(1)+0.5d0
    frac(1,:)=frac(1,:)+0.5d0
  enddo
  do while(com(2).ge.0.5d0)
    com(2)=com(2)-0.5d0
    frac(2,:)=frac(2,:)-0.5d0
  enddo
  do while(com(2).lt.0.0d0)
    com(2)=com(2)+0.5d0
    frac(2,:)=frac(2,:)+0.5d0
  enddo

  990 continue

  !! translate to cartesian coordinate, assuming a0=1
  do i=1,2
  cart(1,i)=sqrt(3.d0)/2.d0*frac(1,i)
  cart(2,i)=-0.5d0*frac(1,i)+frac(2,i)
  cart(3,i)=frac(3,i)*c(3)
  enddo
  comc(1)=sqrt(3.d0)/2.d0*com(1)
  comc(2)=-0.5d0*com(1)+com(2)
  comc(3)=com(3)*c(3)

  if(level.eq.0 .or. level.eq.1)  goto 991 !! only move to a unit cell

  !! move CO to a minimum triangle through Cs symmetries
  ! Cs-1: y=0.d0
  if(comc(2).lt.0.d0) then
    cart(2,:)=-cart(2,:)
    comc(2)=-comc(2)
  endif

  ! CS-2: y=sqrt(3.d0)*x
  if(comc(2).gt.sqrt(3.d0)*comc(1)) then
   !write(91,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
    tmpx=comc(1); tmpy=comc(2)
    comc(1)=-0.5d0*tmpx + 0.5d0*sqrt(3.d0)*tmpy
    comc(2)= 0.5d0*sqrt(3.d0)*tmpx + 0.5d0*tmpy
    do i=1,2
      tmpx=cart(1,i); tmpy=cart(2,i)
      cart(1,i)=-0.5d0*tmpx + 0.5d0*sqrt(3.d0)*tmpy
      cart(2,i)= 0.5d0*sqrt(3.d0)*tmpx + 0.5d0*tmpy
    enddo
   !write(92,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
  endif

  ! CS-3: y=0.25d0 !! remind that the vasp use 2x2 supercell
  if(comc(2).gt.0.25d0) then
   !write(91,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
    comc(2)=0.5d0-comc(2)
    cart(2,:)=0.5d0-cart(2,:)
   !write(92,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
  endif

  ! CS-4: y=sqrt(3.d0)*x-0.5d0 !! remind that the vasp use 2x2 supercell
  if(comc(2).lt.(sqrt(3.d0)*comc(1)-0.5d0)) then
   !write(91,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
    tmpx=comc(1); tmpy=comc(2); !tmpa=sqrt(3.d0)*tmpx-tmpy-0.5d0
    !! x' = x - sqrt(3)*0.5*a; y' = y + 0.5*a
    comc(1)=-0.5d0*tmpx + 0.5d0*sqrt(3.d0)*tmpy + 0.25d0*sqrt(3.d0)
    comc(2)= 0.5d0*sqrt(3.d0)*tmpx + 0.5d0*tmpy - 0.25d0
    do i=1,2
      tmpx=cart(1,i); tmpy=cart(2,i); !tmpa=sqrt(3.d0)*tmpx-tmpy-0.5d0
      cart(1,i)=-0.5d0*tmpx + 0.5d0*sqrt(3.d0)*tmpy + 0.25d0*sqrt(3.d0)
      cart(2,i)= 0.5d0*sqrt(3.d0)*tmpx + 0.5d0*tmpy - 0.25d0
    enddo
   !write(92,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
  endif

  ! CS-5: y=-sqrt(3.d0)*x+0.5d0 !! remind that the vasp use 2x2 supercell
  if(comc(2).gt.(-sqrt(3.d0)*comc(1)+0.5d0)) then
   !write(91,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
    tmpx=comc(1); tmpy=comc(2); tmpa=sqrt(3.d0)*tmpx+tmpy-0.5d0
    !! x' = x - sqrt(3)*0.5*a; y' = y - 0.5*a
    comc(1)=tmpx-0.5d0*sqrt(3.d0)*tmpa
    comc(2)=tmpy-0.5d0*tmpa
    do i=1,2
      tmpx=cart(1,i); tmpy=cart(2,i); tmpa=sqrt(3.d0)*tmpx+tmpy-0.5d0
      cart(1,i)=tmpx-0.5d0*sqrt(3.d0)*tmpa
      cart(2,i)=tmpy-0.5d0*tmpa
    enddo
   !write(92,'(3f10.5)') comc(1)*a0,comc(2)*a0,comc(3)
  endif

  991 continue
  cart(1,:)=cart(1,:)*a0
  cart(2,:)=cart(2,:)*a0
  comc(1)=comc(1)*a0
  comc(2)=comc(2)*a0

  !! translate cartesian to zmat
  zmat(1:3)=comc
  zmat(4)=dsqrt(dot_product(cart(:,1)-cart(:,2),cart(:,1)-cart(:,2)))

  ! phi: angle between CO and z axis
  phi=acos(max(-1.d0,min(1.d0,((cart(3,2)-cart(3,1))/zmat(4)))))/pi*180.d0
  zmat(5)=phi

  ! alpha: angle between (CO's projection on xy plane) and y axis
  dist2d=dsqrt(dot_product(cart(1:2,1)-cart(1:2,2),cart(1:2,1)-cart(1:2,2)))
  alpha=acos(max(-1.d0,min(1.d0,((cart(2,2)-cart(2,1))/dist2d))))/pi*180.d0
  if(cart(1,1).gt.cart(1,2)) alpha=-alpha
  zmat(6)=alpha

  return
end subroutine

subroutine trans2cart(n,frac,cart)
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: frac(3,n)
  real*8,intent(out) :: cart(3,n)

  real*8 :: a(3),b(3),c(3),a0
  integer :: i,j

  !! this is a 2x2 super cell !!
  data a/4.870971576466d0, -2.812256750888d0,  0.000000000000d0/
  data b/0.000000000000d0,  5.624513501775d0,  0.000000000000d0/
  data c/0.000000000000d0,  0.000000000000d0, 29.184792087162d0/
  data a0/5.624513501775d0/

  do i=1,n
    cart(1:3,i)=frac(1,i)*a&
               +frac(2,i)*b&
               +frac(3,i)*c
  enddo
  return
end subroutine

subroutine rotate(xyz,axis,direction,degree)
  !! examples:
  !! call rotate(x(1:3,2),"y","a",45.d0)
  !! call rotate(x(1:3,4),"z","c",120.d0)
  implicit none
  real(kind=8),intent(inout) :: xyz(3)
  character(kind=1,len=1),intent(in) :: axis
  character(kind=1,len=1),intent(in) :: direction
  real(kind=8),intent(in) :: degree
  real(kind=8),parameter :: pi=3.14159265358979323846d0
  real*8 :: tx,ty,tz,ta
  tx=xyz(1)
  ty=xyz(2)
  tz=xyz(3)
  ta=degree/180.d0*pi
  if(direction.ne."a" .and. direction.ne."c") stop "direction = a or c"
  if(axis.ne."x" .and. axis.ne."y" .and. axis.ne."z") stop "axis = x y z"
  if(direction.eq."c") ta=-ta
  if(axis.eq."x") then
    xyz(2)=cos(ta)*ty-sin(ta)*tz
    xyz(3)=sin(ta)*ty+cos(ta)*tz
  elseif(axis.eq."y") then
    xyz(3)=cos(ta)*tz-sin(ta)*tx
    xyz(1)=sin(ta)*tz+cos(ta)*tx
  elseif(axis.eq."z") then
    xyz(1)=cos(ta)*tx-sin(ta)*ty
    xyz(2)=sin(ta)*tx+cos(ta)*ty
  endif
  return
end subroutine

subroutine cart2frac(cart,frac)
  implicit none
  real*8,intent(in) :: cart(3,2)
  real*8,intent(out) :: frac(3,2)
  real*8 :: a(3),b(3),c(3),a0
  integer :: i,j

  !! this is a 2x2 super cell !!
  data a/4.870971576466d0, -2.812256750888d0,  0.000000000000d0/
  data b/0.000000000000d0,  5.624513501775d0,  0.000000000000d0/
  data c/0.000000000000d0,  0.000000000000d0, 29.184792087162d0/
  data a0/5.624513501775d0/

  frac(1,:)=cart(1,:)/a0*2.d0/sqrt(3.d0)
  frac(2,:)=cart(2,:)/a0+0.5d0*frac(1,:)
  frac(3,:)=cart(3,:)/c(3)

  return
end subroutine

subroutine zmat2cart(zmat,cart)
  implicit none
  real*8,intent(in) :: zmat(6)
  real*8,intent(out) :: cart(3,2)
  real*8 :: comc(3),phi,dist2d,alpha,distco
  comc=zmat(1:3)
  distco=zmat(4)
  phi=zmat(5)
  alpha=zmat(6)
  ! put CO along y axis
  cart=0.d0
  cart(2,1)=-distco/2.d0
  cart(2,2)=distco/2.d0
  ! rotate around x axis for (90-phi) anticlockwisely
  call rotate(cart(:,1),'x','a',90.d0-phi)
  call rotate(cart(:,2),'x','a',90.d0-phi)
  ! rotate around z axis for alpha clockwisely
  call rotate(cart(:,1),'z','c',alpha)
  call rotate(cart(:,2),'z','c',alpha)
  ! move the COM
  cart(:,1)=cart(:,1)+comc
  cart(:,2)=cart(:,2)+comc
  return
end subroutine

subroutine cart2zmat(cart,zmat)
  implicit none
  real*8,intent(in) :: cart(3,2)
  real*8,intent(out) :: zmat(6)
  real*8 :: comc(3),phi,dist2d,alpha
  real*8,parameter :: pi=dacos(-1.d0)

  comc=(cart(:,1)+cart(:,2))/2.d0
  !! translate cartesian to zmat
  zmat(1:3)=comc
  zmat(4)=dsqrt(dot_product(cart(:,1)-cart(:,2),cart(:,1)-cart(:,2)))

  ! phi: angle between CO and z axis
  phi=acos(max(-1.d0,min(1.d0,((cart(3,2)-cart(3,1))/zmat(4)))))/pi*180.d0
  zmat(5)=phi

  ! alpha: angle between (CO's projection on xy plane) and y axis
  dist2d=dsqrt(dot_product(cart(1:2,1)-cart(1:2,2),cart(1:2,1)-cart(1:2,2)))
  alpha=acos(max(-1.d0,min(1.d0,((cart(2,2)-cart(2,1))/dist2d))))/pi*180.d0
  if(cart(1,1).gt.cart(1,2)) alpha=-alpha
  zmat(6)=alpha
  return
end subroutine
