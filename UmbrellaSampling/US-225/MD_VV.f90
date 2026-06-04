  Subroutine MD_VV(i_step,zeq,zmc)

  Use ctrl

  Implicit None
  Integer , intent(in) :: i_step
  Double precision , intent(in) :: zeq
  Integer :: iii,jjj,kkk
  Double precision :: rrr
  Double precision , dimension(nfr,1) :: cvec,fvec
  Double precision , dimension(3,natom) :: dedx
!  Double precision , dimension(1) :: egy
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

! Update force !! updated 2022/2/17 @@@@
  if(natom.eq.10) then
    Call pes_ru8co(Cmat,egy,dedx)
  elseif(natom.eq.14) then
    Call pes_ru12co(Cmat,egy,dedx)
  elseif(natom.eq.2) then
    call pes_ruFixco(Cmat,egy,dedx)!Cmat --> the martix?
  else
    stop "natom = 10 or 14"
  endif

  if (flag_umbrella) then

!   Calculate the mass center of CO
!    zmc = (Cmat(3,us%iC) * Mass(us%iC) + Cmat(3,us%iO) * Mass(us%iO)) / (Mass(us%iC) + Mass(us%iO))
    zmc = Cmat(3,us%iC) * us%rmC + Cmat(3,us%iO) * us%rmO

!   Calculate delta z and biased potential
    dz = zmc - zeq
    Ebias = 0.5 * us%fc * dz ** 2
    egy = egy + Ebias

!   The force from biased potential
    dedx(3,us%iC) = dedx(3,us%iC) - us%fc * dz * us%rmC
    dedx(3,us%iO) = dedx(3,us%iO) - us%fc * dz * us%rmO

  end if

  Ev = egy * eeVtoint

  do iii = 1, natom
    do jjj = 1, 3
      Fmat(jjj,iii) = dedx(jjj,iii) * eeVtoint
    end do
  end do

! Second stage
! Update velocity
  if (entype == 1) then
    Call NHC_VV(2)
  else if (entype == 0) then
    Call MD_VV_half()
  end if

  Call CalcTemp()

! Constrain !! updated 2022/2/17 @@@@
! !  flag_reverse = .false.
!   do iii = 1, natom
!     if ((Cmat(3,iii) > zupper).and.(Vmat(3,iii) > 0.0D0)) then
! !      flag_reverse = .true.
!       Vmat(3,iii) = -Vmat(3,iii)
!     end if
!   end do
!  if ( Cmat(3,natom-1) .gt. zupper .and. Vmat(3,natom-1).gt.0.d0) then
!  ! flag_reverse = .true.
!    Vmat(3,natom-1) = -Vmat(3,natom-1) ! C
!    Vmat(3,natom) = -Vmat(3,natom)     ! O
!    write(uout,*) "Desorption (ps) ",i_step*timestep
!  endif

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


