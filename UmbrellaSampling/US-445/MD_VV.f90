  Subroutine MD_VV(i_step,zeq,zmc)

  Use ctrl
  use pes_interface3, only : pes_nCO_RuHD
  Implicit None
  integer,intent(in) :: i_step
  Double precision , intent(in) :: zeq
  Integer :: iii,jjj,kkk
  Double precision :: rrr
  Double precision , dimension(nfr,1) :: cvec,fvec
  Double precision :: dedx(3,natom)
  Double precision :: egy
  Double precision , intent(out) :: zmc   !z coordinate of the mass center
  Double precision :: dz
  Logical :: flag_reverse

! First stage
! Update velocity
  if (entype == 1) then
    Call NHC_VV(1)
  else if (entype == 0) then
    Call MD_VV_half()
  end if

! Update coordinates
  do iii = 1, natom
    do jjj = 1, 3
      Cmat(jjj,iii) = Cmat(jjj,iii) + Vmat(jjj,iii) * timestep
    end do
  end do


!!--------------------> interface to PES <----------------------
! Update force !! updated 2024/6/10 @@@@
! call pes_nCO_pt6x6(natom,Cmat,egy,dedx)

  !! CO molecules on relaxed Pt surface (6x6 2 layers)
  !!@@@@ 2024/6/28, 1-72 surf atoms, 73-natom CO molecules
! if(natom.le.72) stop "wrong INPUT type!"
  egy=0.d0; dedx=0.d0
  call pes_nCO_RuHD(natom,Cmat,egy,dedx)

! print*,"pot energy",egy,"ev"
!!--------------------> end link to PES  <----------------------

  if (flag_umbrella) then

!   Calculate the mass center of CO
!    zmc = (Cmat(3,us%iC) * Mass(us%iC) + Cmat(3,us%iO) * Mass(us%iO)) /
!    (Mass(us%iC) + Mass(us%iO))
    zmc = Cmat(3,us%iC) * us%rmC + Cmat(3,us%iO) * us%rmO

!   Calculate delta z and biased potential
    dz = zmc - zeq
    Ebias = 0.5 * us%fc * dz ** 2
    egy = egy + Ebias

!   The force from biased potential
    dedx(3,us%iC) = dedx(3,us%iC) - us%fc * dz * us%rmC
    dedx(3,us%iO) = dedx(3,us%iO) - us%fc * dz * us%rmO

  end if

  Ev=egy*eeVtoint

  do iii = 1,natom
  do jjj = 1,3
  Fmat(jjj,iii)=dedx(jjj,iii)*eeVtoint
  enddo
  enddo

! Second stage
! Update velocity
  if (entype == 1) then
    Call NHC_VV(2)
  else if (entype == 0) then
    Call MD_VV_half()
  end if

  Call CalcTemp()

!! -----> for CO molecules on relaxed Pt surface (6x6 2 layers)
!!@@@@ 2024/6/28, 1-72 surf atoms, 73-natom CO molecules
!! reverse the velocity in z-direction when CO reach to the upper limit

    do iii=1,(natom-72)/2
  
    if ( Cmat(3,72+iii*2-1) .gt. zupper .and. Vmat(3,72+iii*2-1).gt.0.d0) then
      Vmat(3,72+iii*2-1) = -Vmat(3,72+iii*2-1) ! C
      Vmat(3,72+iii*2) = -Vmat(3,72+iii*2)     ! O
      write(uout,*) "Desorption (ps) / i-th CO ",i_step*timestep,iii
    endif
  
    enddo

! print*,'an integral step finished'
  return
  End Subroutine MD_VV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine MD_VV_half()

  Use ctrl

  Implicit None

  Integer :: iii,jjj,kkk
  Double precision :: hstep

  hstep = 0.5D0 * timestep

  do iii = 1, natom
    do jjj = 1, 3
      Vmat(jjj,iii) = Vmat(jjj,iii) + Fmat(jjj,iii) / mass(iii) * hstep
    end do
  end do

  End Subroutine MD_VV_half

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine Calctemp()

  Use ctrl

  Implicit None

  Integer :: iii,jjj,kkk
  Double precision :: rrr

  Ek = 0.0D0
  do iii = 1, natom
    do jjj = 1, 3
      Ek = Ek + Mass(iii) * Vmat(jjj,iii) ** 2
    end do
  end do

  tempnow = Ek / ektoint / Dble(nfr)
  Ek = Ek * 0.5D0

  End Subroutine Calctemp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine RescaleV(cscale)

  Use ctrl

  Implicit None

  Double precision , intent(in) :: cscale

  Integer :: iii,jjj,kkk

  do iii = 1, natom
    do jjj = 1, 3
      Vmat(jjj,iii) = Vmat(jjj,iii) * cscale
    end do
  end do

  End Subroutine RescaleV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


