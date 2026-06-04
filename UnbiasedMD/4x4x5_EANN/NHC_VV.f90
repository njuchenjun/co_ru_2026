  Subroutine NHC_VV(job)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                   !
! Velocity verlet with Nose-Hoover Chain thermostat !
!                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use ctrl

  Implicit None

  Integer , intent(in) :: job

  if (job == 1) then
! First stage of Nose Hoover chain thermostat
! Update thermostat velocity 
    call NHCINT()

! First stage of velocity verlet algorithm
! Update velocity and position
    call MD_VV_half()
 
  else if (job == 2) then

! Second stage of velocity verlet algorithm
! Update velocity and position
    call MD_VV_half()

! Second stage of Nose Hoover chain thermostat
! Update thermostat velocity
    call NHCINT()

  end if
  
  End Subroutine NHC_VV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine NHCINT()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
! Integrate Nose Hoover chain thermostat from t=0 to t=tstep/2 !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use ctrl

  Implicit None

  Integer ::  i,j,k,INOS,IRESN,IYOSH
  Double precision :: cscale,hstep,AA,AKIN

  cscale = 1.0D0

! Calculate kinetic energy
  Call CalcTemp()
  AKIN = 2.0D0 * Ek

! Update forces
!  GLOGS(1)=(2.0_wp*Ek-GNKT)/QMASS(1)

  if (debug.eq.1) then
    write(udeb,*) 'NHCINT'
    write(udeb,*) 'GLOGS'
    write(udeb,*) (GLOGS(j),j=1,NNOS)
    write(udeb,*) 'QMASS'
    write(udeb,*) (QMASS(j),j=1,NNOS)
    write(udeb,*) 'XLOGS'
    write(udeb,*) (XLOGS(j),j=1,NNOS)
    write(udeb,*) 'VLOGS'
    write(udeb,*) (VLOGS(j),j=1,NNOS)
  end if

! Start the multiple time step procedure
  do IRESN = 1, NRESN
    do IYOSH = 1, NYOSH

! Update the thermostat forces
      GLOGS(NNOS) = (QMASS(NNOS-1)*VLOGS(NNOS-1)*VLOGS(NNOS-1)-GKT)/QMASS(NNOS)
      do INOS = 1, NNOS - 2
        i = NNOS - INOS
        GLOGS(i) = (QMASS(i-1)*VLOGS(i-1)*VLOGS(i-1)-GKT)/QMASS(i)
      end do
      GLOGS(1) = (AKIN-GNKT) / QMASS(1)

! Update the thermostat velocities
      VLOGS(NNOS) = VLOGS(NNOS) + GLOGS(NNOS) * WDTI4(IYOSH)
      do INOS = 1, NNOS - 1
        AA = Exp(-WDTI8(IYOSH)*VLOGS(NNOS+1-INOS))
        VLOGS(NNOS-INOS)=VLOGS(NNOS-INOS)*AA*AA+WDTI4(IYOSH)*GLOGS(NNOS-INOS)*AA
      end do

! Update the particle velolcities and kinetic energy
      AA = Exp(-WDTI2(IYOSH)*VLOGS(1))
      cscale = cscale * AA
      AKIN = AKIN * AA * AA

      if (debug.eq.1) then
        write(udeb,*) 'AA=',AA
      end if

! Update the thermostat positions
      do INOS = 1, NNOS
        XLOGS(INOS) = XLOGS(INOS) + VLOGS(INOS) * WDTI2(IYOSH)
      end do

! Update forces
      GLOGS(1) = (AKIN - GNKT) / QMASS(1)

! Update the thermostat forces and velocities
      do INOS = 1, NNOS - 1
        AA = Exp(-WDTI8(IYOSH)*VLOGS(INOS+1))
        VLOGS(INOS) = VLOGS(INOS)*AA*AA+WDTI4(IYOSH)*GLOGS(INOS)*AA
        GLOGS(INOS+1) = (QMASS(INOS)*VLOGS(INOS)*VLOGS(INOS)-GKT)/QMASS(INOS+1)
      end do
      VLOGS(NNOS) = VLOGS(NNOS) + GLOGS(NNOS) * WDTI4(IYOSH)

      if (debug.eq.1) then
        write(udeb,*) 'Scale=',cscale
      end if
    end do
  end do

  if (debug.eq.1) then
    write(udeb,*) 'NHCINT Rescale factor',cscale
  end if

  call RescaleV(cscale)

  End Subroutine NHCINT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine iniNHC()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                          !
! Initialize Nose Hoover chain thermostat  !
!                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use ctrl

  Implicit None

  Integer :: iii,jjj,kkk
  Double precision , dimension(NYOSH) :: tmp
  Double precision :: Nfreq,ttt
  Double precision :: gaussian

  GNKT = dble(nfr) * Temperature * ektoint
  GKT = Temperature * ektoint

  Nfreq = fcmtoint * fcpl * 2.0D0 * pi

! Set QMASS array
  Qmass(1) = GNKT / Nfreq / Nfreq
  do iii = 2, NNOS
    Qmass(iii) = GKT / Nfreq / Nfreq
  end do

! Set WDTI arrays
  do iii = 1, NYOSH
    tmp(iii) = 1.0D0 / (4.0D0 - 1.5874010519682)
  end do
  tmp(3) = 1.0D0 - 4.0D0 * tmp(1)

  do iii = 1, NYOSH
    WDTI2(iii) = timestep * tmp(iii) / 2.0D0 / dble(NRESN)
    WDTI4(iii) = timestep * tmp(iii) / 4.0D0 / dble(NRESN)
    WDTI8(iii) = timestep * tmp(iii) / 8.0D0 / dble(NRESN)
  end do

  if (debug.eq.1) then
    write(udeb,*) 'Initialize NHC thermostat'
    write(udeb,*) 'GNKT=',GNKT
    write(udeb,*) 'GKT=',GKT
    write(udeb,*) 'WDTI2'
    write(udeb,*) (WDTI2(iii),iii=1,NYOSH)
    write(udeb,*) 'WDTI4'
    write(udeb,*) (WDTI4(iii),iii=1,NYOSH)
    write(udeb,*) 'WDTI8'
    write(udeb,*) (WDTI8(iii),iii=1,NYOSH)
  end if

! Initialize variables of thermostat
  Glogs(1) = 0.0D0
  Xlogs(1) = 1.0D0
  do iii = 2, NNOS
    Glogs(iii) = 0.0D0
    Xlogs(iii) = 1.0D0
  end do

  Vlogs(1) = sqrt(2.0D0 * GNKT / Qmass(1))
  do iii = 2, NNOS
    VLOGS(iii) = gaussian() / sqrt(QMASS(iii))
  end do

  ttt=0.0D0
  do iii = 2, NNOS
    ttt = ttt + VLOGS(iii) * VLOGS(iii) * QMASS(iii)
  end do

  if (ttt.gt.0.000000001) then
    do iii = 2, NNOS
      VLOGS(iii) = VLOGS(iii) * sqrt(GKT * Dble(NNOS-1) / ttt)
    end do
  end if

  End Subroutine iniNHC  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

