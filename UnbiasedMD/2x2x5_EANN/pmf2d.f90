program main
  implicit none
  integer,parameter :: nbatch=10000

  real(kind=8),parameter :: kb=0.0019872041d0 ! kb*T result in kcal/mol
  real(kind=8) :: TK ! Temperature / Kelvin
  character(kind=1,len=99) :: fname
  character(kind=1,len=4) :: ptype
  integer :: ndim
  real(kind=8) :: xmin,xmax
  integer :: xbin
  real(kind=8) :: ymin,ymax
  integer :: ybin

  integer :: i,j,k
  character(kind=1,len=99) :: arg(0:99)
  integer :: narg

  integer :: alive
  real(kind=8) :: tid
  character(kind=1,len=6) :: aid
  integer(kind=8) :: fsize,nline ! nline=fsize/4/6, nline=ns*nt+nr
  integer(kind=8) :: ns,nt,nr    ! read ns one time, repeat nt times, remaining nr lines.
  real(kind=4),allocatable :: cart(:,:),frac(:,:),car6(:,:),fra6(:,:)
  integer(kind=8),allocatable :: site(:,:),hist(:,:),sit6(:,:)
  integer(kind=8) :: ntotal,nmax,nmin
  real(kind=8),allocatable :: pmf(:,:)

  !! read setup and parameters
  i=0; arg="";
  do while(.true.)
     call get_command_argument(i, arg(i))
     if(len_trim(arg(i)) == 0) exit
     i=i+1
     if(i.gt.99) stop "argv too much"
  enddo
  narg=i-1
  if(narg.ne.9) stop "argv not enough. usage:./exec temperature(K) inputfie(binary) pmf2 0 1 100 0 1 100"

  read(arg(1),*) TK
  read(arg(2),*) fname
  read(arg(3),*) ptype
  read(arg(4),*) xmin
  read(arg(5),*) xmax
  read(arg(6),*) xbin
  read(arg(7),*) ymin
  read(arg(8),*) ymax
  read(arg(9),*) ybin

  if(ptype.eq."PMF1" .or. ptype.eq."pmf1") then
      ndim=1
      stop "PMF1 unsupported"
  elseif(ptype.eq."PMF2" .or. ptype.eq."pmf2") then
      ndim=2
  else
      stop "ptype = PMF1 or PMF2"
  endif
  print*,"Temperature / Kelvin: ",TK
  print*,"Read data from: ",trim(fname)
  print*,"PMF type: ",ptype

  if(ndim.eq.1) then

  else
      print*,"hist_minx hist_maxx num_binsx bin_widthx"
      write(*,'(2f10.6,i5,f10.6)') xmin,xmax,xbin,(xmax-xmin)/xbin
      print*,"hist_miny hist_maxy num_binsy bin_widthy"
      write(*,'(2f10.6,i5,f10.6)') ymin,ymax,ybin,(ymax-ymin)/ybin
  endif
  allocate(hist(xbin,ybin),pmf(xbin,ybin))

  !! read data
  call random_seed()
  call random_number(tid)
  tid=tid*999999
  write(aid,'(i6.6)') floor(tid)
  print*,"file size write to: ",'tmp-size-'//aid//'.txt'
  call system('du -b "'//trim(fname)//'" | tee tmp-size-'//aid//'.txt')
  open(101,file='tmp-size-'//aid//'.txt',status='old',action='read')
  read(101,*) fsize
  close(101)
  nline=fsize/4/6
  print*," file size (byte) and lines: ",fsize,nline

  if(nline.le.1000000) then
     ns=nline
     nt=1
     nr=0
  else
     ns=1000000
     nt=nline/ns
     nr=nline-ns*nt
  endif
  print*,"number of lines in a batch: ",ns
  print*,"number of batches: ",nt
  print*,"remaining lines: ",nr

  allocate(cart(6,ns),frac(2,ns),site(2,ns))
  allocate(car6(2,ns*6),fra6(2,ns*6),sit6(2,ns*6))
  hist=0
  open(501,file=trim(fname),status='old',action='read',form='unformatted',access='stream')
  do i=1,nt ! @@@@ debug
     print*,"reading data: ",i," of",nt," million"
     read(501) cart

     ! let top-site be (0,0)
     cart(1,:)=cart(1,:)-0.79117d0 ! from pes_interface.f90
     cart(2,:)=cart(2,:)+0.30431d0 !

    !! version 1: use original MD traj points
    !!call MapSix_Ru0001(ns,cart,car6)
    !frac(1,:)=cart(1,:)/2.3731661157d0
    !frac(2,:)=(cart(2,:)+frac(1,:)*1.37014809573109455d0)/2.7402961914621904d0
    !frac=frac-floor(frac)
    !site(1,:)=ceiling(frac(1,:)*(xbin-2.d-6)+1d-6)
    !site(2,:)=ceiling(frac(2,:)*(ybin-2.d-6)+1d-6)
    !do j=1,ns
    !   hist(site(1,j),site(2,j))=hist(site(1,j),site(2,j))+1
    !enddo

     ! version 2: copy one point to 6 points using symmetry
     call MapSix_Ru0001(ns,cart,car6)
     fra6(1,:)=car6(1,:)/2.3731661157d0
     fra6(2,:)=(car6(2,:)+fra6(1,:)*1.37014809573109455d0)/2.7402961914621904d0
     fra6=fra6-floor(fra6)
     sit6(1,:)=ceiling(fra6(1,:)*(xbin-2.d-6)+1d-6)
     sit6(2,:)=ceiling(fra6(2,:)*(ybin-2.d-6)+1d-6)
     do j=1,ns*6
        hist(sit6(1,j),sit6(2,j))=hist(sit6(1,j),sit6(2,j))+1
     enddo


  enddo

  if(nr.ne.0) then
     print*,nr," lines not used"
  endif

  close(501)
  deallocate(cart,frac,site)
  deallocate(car6,fra6,sit6)

  !! calculate pmf
  ntotal=sum(hist)
  nmax=maxval(hist)
  nmin=minval(hist)
  print*,"total number of lines used: ",ntotal
  print*,"maximum occupation number: ",nmax
  print*,"minimum occupation number: ",nmin
  pmf=-kb*TK*log(1.d0*hist/nmax)

  do i=1,xbin
     do j=1,ybin
        write(602,'(3f10.6,1x,i11)') (i-0.5)/xbin*(xmax-xmin),(j-0.5)/ybin*(ymax-ymin),pmf(i,j),hist(i,j)
     enddo
     write(602,*)
  enddo
  deallocate(hist)
  stop
end program

subroutine MapSix_Ru0001(ns,cart,car6)
  implicit none
  integer,intent(in) :: ns
  real(kind=4),intent(inout) :: cart(6,ns)
  real(kind=4),intent(out) :: car6(2,ns*6)

  real(kind=8) :: x(2),carbon(2),carbo6(2,6)

  integer :: i,j,k

  do i=1,ns
     x=cart(1:2,i) ! x , y of Carbon
     call one_cell_carbon(x,carbon)
     cart(1:2,i)=carbon
     call map_six(carbon(1),carbon(2),carbo6)
     car6(:,(i*6-5):(i*6))=carbo6
  enddo

  return
end subroutine

subroutine one_cell_carbon(x,carbon)
  implicit none
  real*8,intent(in) :: x(2)
  real*8,intent(out) :: carbon(2)
  real*8 :: dax,day,dby
  dax=4.7463322314000056d0/2.d0
  day=-2.7402961914621891d0/2.d0
  dby=5.4805923829243808d0/2.d0

  carbon=x

  do while(carbon(1).gt.0.d0)
     carbon(1)=carbon(1)-dax
     carbon(2)=carbon(2)-day
  enddo

  do while(carbon(1).lt.0.791170d0)
     carbon(1)=carbon(1)+dax
     carbon(2)=carbon(2)+day
  enddo

  do while(carbon(2).gt.-2.d0)
     carbon(2)=carbon(2)-dby
  enddo

  do while(carbon(2).lt.((-1.d0/sqrt(3.d0)*(carbon(1)-0.791170d0))-0.304310d0))
     carbon(2)=carbon(2)+dby
  enddo

  return
end subroutine

subroutine one_cell(x,co)
  implicit none
  real*8,intent(in) :: x(3,2)
  real*8,intent(out) :: co(3,2)
  real*8 :: dax,day,dby
  dax=4.7463322314000056d0/2.d0
  day=-2.7402961914621891d0/2.d0
  dby=5.4805923829243808d0/2.d0

  co=x

  do while(co(1,1).gt.0.d0)
     co(1,1:2)=co(1,1:2)-dax
     co(2,1:2)=co(2,1:2)-day
  enddo

  do while(co(1,1).lt.0.791170d0)
     co(1,1:2)=co(1,1:2)+dax
     co(2,1:2)=co(2,1:2)+day
  enddo

  do while(co(2,1).gt.-2.d0)
     co(2,1:2)=co(2,1:2)-dby
  enddo

  do while(co(2,1).lt.((-1.d0/sqrt(3.d0)*(co(1,1)-0.791170d0))-0.304310d0))
     co(2,1:2)=co(2,1:2)+dby
  enddo

  return
end subroutine

subroutine axsym(p1,p2,xold,xnew)
  implicit none
  real*8 :: p1(2),p2(2),xold(2),xnew(2)
  real*8 :: x1,y1,x2,y2,x0,y0,A,B,C
  x0=xold(1)
  y0=xold(2)

  x1=p1(1)
  y1=p1(2)
  x2=p2(1)
  y2=p2(2)

  ! line: Ax+By+C=0
  A=y1-y2
  B=x2-x1
  C=x1*y2-x2*y1

  xnew(1)=x0-A*2.d0*(A*x0+B*y0+C)/(A**2+B**2)
  xnew(2)=y0-B*2.d0*(A*x0+B*y0+C)/(A**2+B**2)
  return
end subroutine

subroutine map_six(x,y,xy)
  implicit none
  real*8,parameter :: a0=5.4805923829243808d0/2.d0
  real*8,intent(in) :: x,y
  real*8,intent(out) :: xy(2,6)
  real*8 :: orange(2),xt,yt

  xt=x; yt=y
  if(yt.gt.-sqrt(3.d0)*xt+a0) then
      call axsym([0.d0,a0],[a0/sqrt(3.d0),0.d0],[xt,yt],orange)
      xt=orange(1); yt=orange(2)
  endif
  if(yt.gt.sqrt(3.d0)*xt) then
      call axsym([0.d0,0.d0],[1.d0,sqrt(3.d0)],[xt,yt],orange)
      xt=orange(1); yt=orange(2)
  endif
  if(yt.lt.0.d0) then
      yt=-yt
  endif
  if(yt.gt.-sqrt(3.d0)*xt+a0) then
      call axsym([0.d0,a0],[a0/sqrt(3.d0),0.d0],[xt,yt],orange)
      xt=orange(1); yt=orange(2)
  endif

  xy(:,1)=[xt,yt]
  if(yt.gt.xt/sqrt(3.d0)) then
    call axsym([0.d0,0.d0],[1.d0,sqrt(3.d0)],xy(:,1),xy(:,2))
    call axsym([0.d0,a0/2.d0],[1.d0,a0/2.d0],xy(:,2),xy(:,3))
  else
    call axsym([0.d0,0.d0],[1.d0,0.d0],xy(:,1),xy(:,2))
    call axsym([a0/sqrt(3.d0),0.d0],[a0*sqrt(3.d0)/2.d0,a0/2.d0],xy(:,2),xy(:,3))
  endif
  call axsym([0.d0,a0],[a0/sqrt(3.d0),0.d0],xy(:,1),xy(:,4))
  call axsym([0.d0,a0],[a0/sqrt(3.d0),0.d0],xy(:,2),xy(:,5))
  call axsym([0.d0,a0],[a0/sqrt(3.d0),0.d0],xy(:,3),xy(:,6))
  return
end subroutine

