  Subroutine MD_VV(i_step)

  Use ctrl
  use pes_interface3, only : pes_nCO_RuHD
  Implicit None
  integer,intent(in) :: i_step
  Integer :: iii,jjj,kkk
  Double precision :: rrr
  Double precision , dimension(nfr,1) :: cvec,fvec
  Double precision :: dedx(3,natom)
  Double precision :: egy

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
! call pes_nCO_ru6x6(natom,Cmat,egy,dedx)

  !! CO molecules on relaxed ru surface (4x4 2 layers)
  !!@@@@ 2024/6/28, 1-32 surf atoms, 33-natom CO molecules
! if(natom.le.32) stop "wrong INPUT type!"
  egy=0.d0; dedx=0.d0
  call pes_nCO_RuHD(natom,Cmat,egy,dedx)

! print*,"pot energy",egy,"ev"
!!--------------------> end link to PES  <----------------------

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

!! -----> for CO molecules on relaxed ru surface (4x4 3 layers)
!!@@@@ 2024/6/28, 1-32 surf atoms, 33-natom CO molecules
!! reverse the velocity in z-direction when CO reach to the upper limit

    do iii=1,(natom-32)/2
  
    if ( Cmat(3,32+iii*2-1) .gt. zupper .and. Vmat(3,32+iii*2-1).gt.0.d0) then
      Vmat(3,32+iii*2-1) = -Vmat(3,32+iii*2-1) ! C
      Vmat(3,32+iii*2) = -Vmat(3,32+iii*2)     ! O
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


