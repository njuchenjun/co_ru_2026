program main
  implicit none
  real*8 :: temp
  character(kind=1,len=3) :: atemp
  character(kind=1,len=4) :: aid
  character(kind=1,len=99) :: ainp(40),binp
  character(kind=1,len=10) :: cinp,dinp,einp
  integer :: id,j,i,k,l,outstep,step
  real*8 :: x_or,y_or,x_o,y_o,x_r,y_r,x1,y1,x2,y2
  real*8 :: a(3),b(3)

  !call random_seed()
  data a/4.7463322314000056d0,-2.7402961914621891d0, 0.0000000000000000d0/
  data b/0.0000000000000000d0, 5.4805923829243808d0, 0.0000000000000000d0/

 !temp=300.d0
 ! read*,temp
 ! print*,temp
 ! read*,step
 ! print*,step
 ! read*,outstep
 ! print*,outstep
 ! write(atemp,'(i3.3)') int(temp)
 ! call system("mkdir run-"//atemp)
  open(101,file='./INPUT',status='old',action='read')
  do l=1,40
     read(101,'(a99)') ainp(l)
  enddo
  close(101)

  call random_seed()

     !call random_number(xrand); xrand=xrand*5.d0
     !call random_number(yrand); yrand=yrand*5.d0
     call random_number(x_r)
     x_r = (x_r*2.0-1.0) * 0.25
     call random_number(y_r)
     y_r = (y_r*2.0-1.0) * 0.15
     call random_number(x_or)
     call random_number(y_or)
     x_o = x_r+(x_or*2.0-1.0)*0.001
     y_o = x_r+(y_or*2.0-1.0)*0.001
     write(*,*) x_r,y_r
     write(*,*) x_o,y_o

     binp=ainp(28)!change the coordination of C,32 mean the C line in the INPUT text 
     dinp=binp(21:30)
     einp=binp(31:40)
     read(dinp,'(f10.5)') x1
     read(einp,'(f10.5)') y1
     write(*,*) x1,y1
     x1=x1+x_r
     y1=y1+y_r
     write(*,*) x1,y1
     write(dinp,'(f10.5)') x1
     write(einp,'(f10.5)') y1
     binp(21:30)=dinp
     binp(31:40)=einp
     ainp(28)=binp

     binp=ainp(29) !for O
     dinp=binp(21:30)
     einp=binp(31:40)
     read(dinp,'(f10.5)') x2
     read(einp,'(f10.5)') y2
     write(*,*) x2,y2
     x2=x2+x_o
     y2=y2+y_o
     write(*,*) x2,y2
     write(dinp,'(f10.5)') x2
     write(einp,'(f10.5)') y2
     binp(21:30)=dinp
     binp(31:40)=einp
     ainp(29)=binp

    ! binp=ainp(36)!the window index
    ! write(cinp,'(I4)') id
    ! binp(1:4)=cinp
    ! ainp(36)=binp


     open(114,file="./INPUT",status='old',action='write')
     do l=1,40!lines of the input file
        write(114,'(a)') trim(ainp(l))
     enddo
     close(114)
     !call system("ln -s /home/chenjun/CO-Ru0001/3-xMD/2-xMD-HD/src-3layers-binary/pes8.dat run-"//atemp//"/run-"//aid//"/pes8.dat")
     !call system("ln -s /home/hejinen/usdata/with6dpes/xmd.x run-"//atemp//"/run-"//aid//"/xmd")

end program
