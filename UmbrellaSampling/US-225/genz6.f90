program main
  implicit none
  real*8 :: temp
  character(kind=1,len=3) :: atemp
  character(kind=1,len=4) :: aid
  character(kind=1,len=99) :: ainp(40),binp
  character(kind=1,len=10) :: cinp,dinp,einp
  integer :: id,j,i,k,l,outstep,step,m,n
  double precision :: invm,invn
  real*8 :: czcoor,ozcoor,x_o,y_o,x_r,y_r
  real*8 :: a(3),b(3)

  !call random_seed()
  data a/4.7463322314000056d0,-2.7402961914621891d0, 0.0000000000000000d0/
  data b/0.0000000000000000d0, 5.4805923829243808d0, 0.0000000000000000d0/

 !temp=300.d0
  write(*,*) 'Temperature:'
  read*,temp
  print*,temp
  write(*,*) 'step:'
  read*,step
  print*,step
  write(*,*) 'outstep:'
  read*,outstep
  print*,outstep
  write(*,*) 'm:'
  read*,m
  print*,m
  write(*,*) 'n:'
  read*,n
  print*,n
  write(atemp,'(i3.3)') int(temp)
  call system("mkdir run-"//atemp)
  open(101,file='./INPUT',status='old',action='read')
  do l=1,40
     read(101,'(a99)') ainp(l)
  enddo
  close(101)

  call random_seed()
  invm=1.0d0/m
  write(*,*) invm
  invn=1.0d0/n
  write(*,*) invn
  do i=1,m
      do j=1,n
         do k=1,64
         id = k+(j-1)*64+(i-1)*n*64
         write(aid,'(i4.4)') id
     write(ainp(6),*) outstep
     write(ainp(4),*) step
     !call random_number(xrand); xrand=xrand*5.d0
     !call random_number(yrand); yrand=yrand*5.d0
     call random_number(x_r)
     x_r = x_r * invm + (i-1)*invm
     call random_number(y_r)
     y_r = y_r * invn + (j-1)*invn
     x_o = x_r * a(1)/2 + y_r * b(1)/2 + 0.79117
     y_o = x_r * a(2)/2 + y_r * b(2)/2 - 0.30431

     czcoor = 8.47202 + 1.60d0 + (k-1)*0.05000
     binp=ainp(28)!change the coordination of C,32 mean the C line in the INPUT text 
     write(cinp,'(f10.5)') czcoor!change z value 20:29 mean z value range
     write(dinp,'(f10.5)') x_o
     write(einp,'(f10.5)') y_o
     binp(20:29)=dinp
     binp(30:39)=einp
     binp(40:49)=cinp
     ainp(28)=binp

     ozcoor = czcoor + 1.20000
     binp=ainp(29) !for O
     write(cinp,'(f10.5)') ozcoor
     write(dinp,'(f10.5)') x_o
     write(einp,'(f10.5)') y_o
     binp(20:29)=dinp
     binp(30:39)=einp
     binp(40:49)=cinp
     ainp(29)=binp

     binp=ainp(10)
     write(cinp,'(f10.2)') temp
     binp(1:10)=cinp
     ainp(10)=binp

     binp=ainp(36)!the window index
     write(cinp,'(I4)') id
     binp(1:4)=cinp
     ainp(36)=binp

     binp=ainp(38)!the z0
     write(cinp,'(f10.5)') czcoor
     binp(1:10)=cinp
     ainp(38)=binp

     call system("mkdir run-"//atemp//"/run-"//aid)
     open(200+id,file="run-"//atemp//"/run-"//aid//"/INPUT",status='unknown')
     do l=1,40!lines of the input file
        write(200+id,*) trim(ainp(l))
     enddo
     close(200+id)
     !call system("ln -s /home/chenjun/CO-Ru0001/3-xMD/2-xMD-HD/src-3layers-binary/pes8.dat run-"//atemp//"/run-"//aid//"/pes8.dat")
     !call system("ln -s /home/hejinen/usdata/with6dpes/xmd.x run-"//atemp//"/run-"//aid//"/xmd")
   enddo
   enddo
  enddo

end program
