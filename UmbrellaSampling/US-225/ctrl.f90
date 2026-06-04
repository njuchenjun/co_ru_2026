Module ctrl

  Use const
  Use thermostat

  Implicit None

  Character (len=3) :: ensemble
  Integer (kind=8) :: nstep      
  Integer (kind=8) :: outstep    
  Double precision :: timestep     ! ps
  Double precision :: Temperature  ! K
  Double precision :: lconst       ! A
  Integer , dimension(2) :: supercell
  Integer :: natom
  Double precision :: zupper
  Double precision :: Vlower = 0.0D0
  Character (len=2) , dimension(:) , allocatable :: element
  Double precision , dimension(:) , allocatable :: Mass     ! au
  Double precision , dimension(:,:) , allocatable :: Cmat   ! A
  real*4 , dimension(:,:) , allocatable :: Cmatsingle   ! A
  Double precision , dimension(:,:) , allocatable :: Vmat   ! A/ps
  Double precision , dimension(:,:) , allocatable :: Fmat   ! 10J/mol/A

  Double precision , dimension(2,2) :: lattice_vector       ! lattice vector
  Double precision , dimension(2) :: lsc 

  Integer :: nfr ! Number of freedom

  Integer :: debug = 0
  Integer :: entype = 0

  Double precision :: Ek,Ev,Etot,tempnow

  Logical :: flag_umbrella = .false.

  Type :: umbrella_sampling

    Integer :: ius
!   Window index

    Double precision :: z0
!   The center point of each window

    Double precision :: fc
!   Force constant: biased function W(z) = 0.5 * fc * (z - z0)

    Integer :: iC, iO
!   The index of C and O atom

    Double precision :: rmC,rmO
!   rmC = mass(iC) / (mass(iC) + mass(iO)); rmO = mass(iO) / (mass(iC) + mass(iO))

  End Type umbrella_sampling

  Type(umbrella_sampling) :: us

  Double precision :: Ebias
! Biased term of potential energy

End Module ctrl
