  Subroutine simulation()

  Use ctrl

  Implicit None

  Integer (kind=8) :: i_step
  Integer :: iii,jjj,kkk

  Double precision :: avg_T,var_T,expt_var_T

  write(uout,"(12x,A,6x,A,6x,A,6x,A)") 'Step:','Temperature (K)', &
                                      'Potential (eV)','Total energy ( eV)'

  write(uout,"(80A)") ('-',iii=1,80)

  avg_T = 0.0D0
  var_T = 0.0D0

  do i_step = 1, nstep

!    print*,'i_step',i_step

    Call MD_VV(i_step,0.0D0)
!   Call constrain()

    avg_T = avg_T + tempnow
    var_T = var_T + (tempnow - Temperature) ** 2

    if (Mod(i_step,outstep) == 0) then
      write(uout,"(2x,I14,6x,F12.3,2ES18.8)") i_step,tempnow,Ev/1.d2/96.4853d0,(Ev+Ek)/1.d2/96.4853d0
                                                           !! from 10J/mol to eV
      if (mod(i_step,outstep).eq.0) then
         write(utrjtxt,*) natom
         write(utrjtxt,"(2x,I14,ES18.8)") i_step,Ev/1.d2/96.4853d0
         do iii = 1, natom
            write(utrjtxt,"(2x,A,2x,3E16.6)") element(iii),(Cmat(jjj,iii),jjj=1,3)
         end do
      endif

      Cmatsingle=Cmat
      write(utrj) Cmatsingle(:,(natom-1):natom)
      write(utrjsurf) Cmatsingle(:,1:(natom-2))

     !write(uvel,"(2x,A,2x,I14)") 'Step:',i_step
     !do iii = 1, natom
     !  write(uvel,"(2x,A,2x,3E16.6)") element(iii),(Vmat(jjj,iii),jjj=1,3)
     !end do

      if ((Ev < Vlower).and.(Vlower < 0.0D0)) then
        write(uout,"(2x,A)") "Error: potential energy too low!"
        exit
      end if
    end if
  end do

  write(uout,"(80A)") ('-',iii=1,80)
  write(uout,"(2x,A)") 'End of simulation.'

  write(uout,"(A)") ''

  avg_T = avg_T / dble(nstep)
  var_T = var_T / dble(nstep)
  expt_var_T = 2.0D0 * temperature ** 2 / dble(nfr)
  write(uout,"(2x,A,F14.3)") 'Expected temperature: ',Temperature
  write(uout,"(2x,A,F14.3)") 'Averaged temperature: ',avg_T
  write(uout,"(A)") ''

  write(uout,"(2x,A,F14.3)") 'Expected variance of temperature: ',expt_var_T
  write(uout,"(2x,A,9x,F14.3)") 'Variance of temperature: ',var_T
  write(uout,"(A)") ''

  End Subroutine simulation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine constrain()

  Use ctrl

  Implicit None

  Integer :: iii,jjj,kkk
  Double precision :: rrr,xxx,yyy

  do iii = 1, natom

    xxx = Cmat(1,iii) / lattice_vector(1,1)

    if (xxx < 0.0D0) then

      xxx = xxx + lsc(1)
      Cmat(1,iii) = Cmat(1,iii) + lsc(1) * lattice_vector(1,1)
      Cmat(2,iii) = Cmat(2,iii) + lsc(1) * lattice_vector(2,1)

    else if (xxx > lsc(1)) then

      xxx = xxx - lsc(1)
      Cmat(1,iii) = Cmat(1,iii) - lsc(1) * lattice_vector(1,1)
      Cmat(2,iii) = Cmat(2,iii) - lsc(1) * lattice_vector(2,1)

    end if

    yyy = Cmat(2,iii) - xxx * lattice_vector(2,1)

    if (yyy < 0.0D0) then

      yyy = yyy + lsc(2)
      Cmat(1,iii) = Cmat(1,iii) + lsc(2) * lattice_vector(1,2)
      Cmat(2,iii) = Cmat(2,iii) + lsc(2) * lattice_vector(2,2)

    end if

    if (yyy > lsc(2)) then

      yyy = yyy - lsc(2)
      Cmat(1,iii) = Cmat(1,iii) - lsc(2) * lattice_vector(1,2)
      Cmat(2,iii) = Cmat(2,iii) - lsc(2) * lattice_vector(2,2)

    end if

!    print*,'atom',iii
!    print*,'xxx',xxx,'yyy',yyy

  end do

  End Subroutine constrain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine umbrella()

  Use ctrl

  Implicit None

  Integer (kind=8) :: i_step
  Integer :: iii,jjj,kkk,ius

  Double precision :: avg_T,var_T,expt_var_T,zmc

  Character (len=20) :: WINNAME,WINC

  write(uout,"(2x,A)") 'Umbrella sampling starts...'

  write(uout,"(2x,A,2x,I6)") 'Sampling simulation for window ',us%ius
  write(uout,"(2x,A,F14.5,A,F14.4,A)") 'Equilibrium position: ', us%z0,' A,   Force constant: ',us%fc,' eV/A/A'

  write(uout,"(12x,A,2x,A18,2x,A18,2x,A18,2x,A18)") 'Step:','Temperature (K)', &
                                      'Potential (eV)','Biased V (eV)','Total energy (eV)'

  write(uout,"(80A)") ('-',iii=1,80)

  avg_T = 0.0D0
  var_T = 0.0D0

  write(WINNAME,"(I4)") us%ius
  write(WINC,"(I4)") us%ius
  do iii = 1, 4
    if (WINNAME(iii:iii) == ' ') then
      WINNAME(iii:iii) = '0'
    end if

    if (WINC(iii:iii) == ' ') then
      WINC(iii:iii) = '0'
    end if

  end do

  WINNAME = 'US-a'//WINNAME
  WINC = 'tra-a'//WINC

  open(unit=uwin,file=WINNAME)
  open(unit=uwc,file=WINC)

  do i_step = 1, nstep

!    print*,'i_step',i_step

    Call MD_VV(i_step,us%z0,zmc)
!   Call constrain()

    avg_T = avg_T + tempnow
    var_T = var_T + (tempnow - Temperature) ** 2

    if (Mod(i_step,outstep) == 0) then
      write(uout,"(2x,I14,2x,F20.6,3ES20.8)") i_step,tempnow,Ev/eeVtoint,Ebias,(Ev+Ek)/eeVtoint ! eV
      if (mod(i_step,outstep).eq.0) then
         write(utrjtxt,*) natom
         write(utrjtxt,"(2x,I14,2ES18.8)") i_step,Ev/eeVtoint,Ebias ! eV
         do iii = 1, natom
            write(utrjtxt,"(2x,A,2x,3E16.6)") element(iii),(Cmat(jjj,iii),jjj=1,3)
         end do
      endif

      Cmatsingle=Cmat
      write(utrj) Cmatsingle(:,(natom-1):natom)
      write(utrjsurf) Cmatsingle(:,1:(natom-2))

     !write(uvel,"(2x,A,2x,I14)") 'Step:',i_step
     !do iii = 1, natom
      ! write(uvel,"(2x,A,2x,3E16.6)") element(iii),(Vmat(jjj,iii),jjj=1,3)
     !end do

!    Output umbrella sampling
    if(natom.eq.2) then !input contral
      write(uwc,*) natom+4
      write(uwc,"(1x,f10.5,A7,2x,I15)") Ev/100.d0/96.4853d0 ! 'step:' i_step - 1
      do iii = 1,natom
        write(uwc,"(A2,3F12.4)") element(iii),(Cmat(jjj,iii),jjj = 1,3)
      end do
      write(uwc,"(A2,3F12.4)") "Ru",0.79117,   2.43599,  0.00000
      write(uwc,"(A2,3F12.4)") "Ru",3.16433,   1.06584,  0.00000
      write(uwc,"(A2,3F12.4)") "Ru",0.79117,  -0.30431,  0.00000
      write(uwc,"(A2,3F12.4)") "Ru",3.16433,  -1.67446,  0.00000
    else
      write(uwc,*) natom
      write(uwc,"(1x,f10.5,A7,2x,I15)") Ev/100.d0/96.4853d0
      do iii = 1,natom
        write(uwc,"(A2,3F12.4)") element(iii),(Cmat(jjj,iii),jjj = 1,3)
      end do
    end if
    ! write(uwc,"(A,2x,I15)") 'Step: ', i_step - 1
     !do iii = 1, natom
      ! write(uwc,"(A6,3F12.4)") element(iii),(Cmat(jjj,iii),jjj = 1,3)
     !end do !the original code

     write(uwin,"(I15,2x,F16.8,2x,F16.8)") i_step - 1,zmc,us%z0

      if ((Ev < Vlower).and.(Vlower < 0.0D0)) then
        write(uout,"(2x,A)") "Error: potential energy too low!"
        exit
      end if
    end if
  end do

  close(uwin)

  write(uout,"(80A)") ('-',iii=1,80)
  write(uout,"(2x,A)") 'End of simulation.'

  write(uout,"(A)") ''

  avg_T = avg_T / dble(nstep)
  var_T = var_T / dble(nstep)
  expt_var_T = 2.0D0 * temperature ** 2 / dble(nfr)
  write(uout,"(2x,A,F14.3)") 'Expected temperature: ',Temperature
  write(uout,"(2x,A,F14.3)") 'Averaged temperature: ',avg_T
  write(uout,"(A)") ''

  write(uout,"(2x,A,F14.3)") 'Expected variance of temperature: ',expt_var_T
  write(uout,"(2x,A,9x,F14.3)") 'Variance of temperature: ',var_T
  write(uout,"(A)") ''


  End Subroutine umbrella

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


