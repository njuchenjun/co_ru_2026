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

!   print*,'i_step',i_step

    Call MD_VV(i_step)

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
     !write(utrj) Cmatsingle(:,1:natom)
     !write(utrj) Cmatsingle(:,(natom-1):natom)
     !write(utrjsurf) Cmatsingle(:,1:(natom-2))

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


