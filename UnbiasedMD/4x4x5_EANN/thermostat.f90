Module thermostat

  Implicit None

  Integer , Parameter :: NNOS = 6
! NNOS  : Number of Nose Hoover thermostats.

  Integer , Parameter :: NYOSH = 5
! Nys in ref

  Integer , Parameter :: NRESN = 2
! Nc in ref

  Double precision , dimension(NYOSH) :: WDTI2,WDTI4,WDTI8
! WDTI2  : wj*tstep/2/nc
! WDTI4  : wj*tstep/4/nc
! WDTI8  : wj*tstep/8/nc

  Double precision  :: GNKT,GKT
! GNKT   : GNKT = Nfr*kB*T
! GKT    : GKT = kB*T

  Double precision , dimension(NNOS) :: Qmass,Glogs,Vlogs,Xlogs
! Qmass : Mass of NHC thermostats

  Double precision :: fcpl
! fcpl : couple frequency

End Module thermostat

