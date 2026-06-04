

program main
  implicit none
  integer :: natom,na,ng,j,nsurf
  integer,parameter :: nmax=200
  real*4 :: cmat(3,nmax) !! single precision
  character(kind=1,len=2) :: elem(nmax)
  character(kind=1,len=40) :: infile,outfile,surfa
  integer*8 :: i_step
  real*8 :: vpot

  integer :: utrj=9005
  call getarg(1,infile)
  call getarg(2,outfile)
  call getarg(3,surfa)
  read(surfa,*) nsurf
  ! usage: ./exec inputfile outfile surfatomsnum

  open(unit=utrj,file=outfile,status='unknown',form='unformatted',access='stream') ! full traj of CO adsorbate

  open(1001,file=infile,status='old',action='read')
  read(1001,*,end=991) natom
  rewind(1001)

  do while(.true.)
     read(1001,*,end=991) na
     read(1001,*,end=991) i_step,vpot
     do j=1,natom
        read(1001,*,end=991) elem(j),cmat(1:3,j)
     enddo

     write(utrj) cmat(:,nsurf+1:natom) !!@@@@ for HD pes, 4*4 2 layers of Ru atoms

  enddo
  991 continue
  close(1001)
  close(utrj)
end program

