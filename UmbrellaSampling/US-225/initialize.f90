  Subroutine Initialize()

  Use ctrl

  Implicit None

  Character(len=80) :: blankline
  Integer :: iii,jjj,kkk
  Double precision :: cscale
  Double precision :: gaussian
  Double precision :: aaa,bbb,lwin,start

! read control
  Open(unit=uinp,file='INPUT')
  read(uinp,*) blankline
  read(uinp,*) ensemble
  read(uinp,*) blankline
  read(uinp,*) nstep
  read(uinp,*) blankline
  read(uinp,*) outstep
  read(uinp,*) blankline
  read(uinp,*) timestep
  read(uinp,*) blankline
  read(uinp,*) Temperature
  read(uinp,*) blankline
  read(uinp,*) fcpl
  read(uinp,*) blankline
  read(uinp,*) lconst
  read(uinp,*) blankline
  read(uinp,*) supercell(1),supercell(2)
  read(uinp,*) blankline
  read(uinp,*) natom

  Allocate(element(natom))
  Allocate(Mass(natom))
  Allocate(Cmat(3,natom))
  Allocate(Cmatsingle(3,natom))
  Allocate(Vmat(3,natom))
  Allocate(Fmat(3,natom))

  read(uinp,*) blankline
  do iii = 1, natom
    read(uinp,*) element(iii),Mass(iii),(Cmat(jjj,iii),jjj=1,3)
  end do

  read(uinp,*) blankline
  read(uinp,*) zupper
  read(uinp,*) blankline
  if (blankline /= '') then
    read(uinp,*) Vlower
  end if

  read(uinp,*) blankline
  if (blankline(1:3) /= 'END') then
    flag_umbrella = .true.
    read(uinp,*) blankline
    read(uinp,*) us%ius
    read(uinp,*) blankline
    read(uinp,*) us%z0
    read(uinp,*) blankline
    read(uinp,*) us%fc
  end if

  Close(uinp)

!! updated 2022/2/17 @@@@ chenjun
! open(unit=uout,file='OUTPUT')
! open(unit=utrj,file='TRAJECTORY')
! open(unit=uvel,file='VELOCITY')
  open(unit=uout,file='traj.log',status='unknown') ! a short output
  open(unit=utrj,file='traj_co.dat',status='unknown',form='unformatted',access='stream') ! full traj of CO adsorbate
  open(unit=utrjsurf,file='traj_surf.dat',status='unknown',form='unformatted',access='stream') ! full traj of surface
  open(unit=utrjtxt,file='traj.txt',status='unknown') ! partial movie molden


  if (debug == 1) then
    open(unit=udeb,file='DEBUG')
  end if

  if (ensemble == 'NVE') then
    entype = 0
  else if (ensemble == 'NVT') then
    entype = 1
  end if

! Output information
  write(uout,"(2x,A)") 'MD simulation for CO diffusion'
  write(uout,*) ''
  write(uout,"(2x,A,2x,A)") 'Ensemble:',ensemble
  write(uout,*) ''
  write(uout,"(2x,A,I14)") 'Number of steps:',nstep
  write(uout,*) ''
  write(uout,"(2x,A,I14)") 'Output interval:',outstep
  write(uout,*) ''
  write(uout,"(2x,A,E14.6)") 'Time step:',timestep
  write(uout,*) ''
  write(uout,"(2x,A,F12.3)") 'Temperature:',Temperature
  write(uout,*) ''
  write(uout,"(2x,A,E14.6)") 'Lattice constant:',lconst
  write(uout,*) ''
  write(uout,"(2x,A,2x,I8,2x,I8)") 'Supercell:',supercell(1),supercell(2)
  write(uout,*) ''
  write(uout,"(2x,A,I14)") 'Number of atoms',natom
  write(uout,*) ''
  write(uout,"(2x,A)") 'Coordinates'
  do iii = 1, natom
    write(uout,"(2x,A,2x,3E14.6)") element(iii),Cmat(1,iii), &
                                   Cmat(2,iii),Cmat(3,iii)
  end do
  write(uout,*) ''
  write(uout,"(2x,A,F16.8)") 'Upper bound of z direction:',zupper
  write(uout,*) ''

  write(uout,"(2x,A,F16.8,A)") 'Lower bound of potential energy:',Vlower,' eV'
  Vlower = Vlower * eeVtoint
  write(uout,"(2x,A,E18.8,A)") 'Lower bound of potential energy:',Vlower
  write(uout,*) ''


! Number of freedom
  nfr = natom * 3
  write(uout,"(2x,A,I8)") 'Number of freedom:',nfr
  write(uout,*) ''

  if (debug == 0) then
    write(uout,"(2x,A)") 'Debug:          disable'
  else
    write(uout,"(2x,A)") 'Debug:           enable'
  end if
  write(uout,*) ''

! Lattice vector
  lattice_vector(1,1) = sqrt(3.0D0) / 2.0D0
  lattice_vector(2,1) = -0.5D0
  lattice_vector(1,2) = 0.0D0
  lattice_vector(2,2) = 1.0D0

! Surface lattice
  lsc(1) = dble(supercell(1)) * lconst
  lsc(2) = dble(supercell(2)) * lconst

! Randomly assign velocity
  do iii = 1, natom
    do jjj = 1, 3
      Vmat(jjj,iii) = gaussian() / sqrt(Mass(iii))
    end do
  end do

  Call Calctemp()

  cscale = sqrt(Temperature / tempnow)

! Rescale velocity
  Call RescaleV(cscale)

! Initialize NHC
  Call iniNHC()

! Initialize umbrella sampling if necessary
  if (flag_umbrella) then

!   Find the index of C/O atom
    do iii = 1, natom
      if (element(iii) == 'C') then
        us%iC = iii
      end if

      if (element(iii) == 'O') then
        us%iO = iii
      end if
    end do

!   Calculate the ratio of mass
    us%rmC = Mass(us%iC) / (Mass(us%iC) + Mass(us%iO))
    us%rmO = Mass(us%iO) / (Mass(us%iC) + Mass(us%iO))

!   Transfer fc from kcal/mol to eV
    us%fc = us%fc / eeVtokcal

  end if

  End Subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Function gaussian()

!    *******************************************************************
!    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         **
!    **                                                               **
!    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
!    **    ADDISON-WESLEY), 1978                                      **
!    **                                                               **
!    *******************************************************************

  Implicit None

  Double precision :: gaussian

  Double precision, Parameter :: a1 = 3.949846138
  Double precision, Parameter :: a3 = 0.252408784
  Double precision, Parameter :: a5 = 0.076542912
  Double precision, Parameter :: a7 = 0.008355968
  Double precision, Parameter :: a9 = 0.029899776

  Integer :: iii,jjj,kkk
  Double precision :: rrr,rr2
  Double precision , dimension(12) :: uni

  uni = 0.0D0
  do iii = 1, 12
    Call Random_number(uni(iii))
  end do

  rrr = (sum(uni) - 6.0D0) / 4.0D0
  rr2 = rrr * rrr
  gaussian = rrr * (a1 + rr2 * (a3 + rr2 * (a5 + rr2 * (a7 + rr2 * a9))))

  End Function gaussian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



