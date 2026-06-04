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

End Module ctrl
