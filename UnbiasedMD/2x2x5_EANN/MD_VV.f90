  Subroutine MD_VV(i_step)

  Use ctrl

  Implicit None
  integer,intent(in) :: i_step
  Integer :: iii,jjj,kkk
  Double precision :: rrr
  Double precision , dimension(nfr,1) :: cvec,fvec
  Double precision , dimension(3,natom) :: dedx
!  Double precision , dimension(1) :: egy
  Double precision :: egy
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

! Update force !! updated 2024/6/10 @@@@
 !if(natom.eq.72) then
 !  call pes_co36_pt6x6(Cmat,egy,dedx)
 !elseif(natom.eq.54 .or. natom.eq.36 .or. natom.eq.24 .or. natom.eq.18 .or. natom.eq.144 .or. natom.eq.56) then
 !  call pes_nCO_pt6x6(natom,Cmat,egy,dedx)
 !else
 !  stop "natom = 72/54/36/24/18"
 !endif
 if (natom.eq.10) then
   call pes_ru8co(Cmat,egy,dedx)
 else if(natom .eq. 14) then
   call pes_ru12co(Cmat,egy,dedx)
 else if(natom .eq. 2) then
   call pes_ruFIXco(Cmat,egy,dedx)
 else
   stop "MD for Ru-CO , natom = 2 or 10 or 14"
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

!! constrain original, by thshen @@ last two atoms are CO
!! reverse the velocity in z-direction when CO reach to the upper limit

! if ( Cmat(3,natom-1) .gt. zupper .and. Vmat(3,natom-1).gt.0.d0) then
! ! flag_reverse = .true.
!   Vmat(3,natom-1) = -Vmat(3,natom-1) ! C
!   Vmat(3,natom) = -Vmat(3,natom)     ! O
!   write(uout,*) "Desorption (ps) ",i_step*timestep
! endif

!! constrain original, by thshen @@ last two atoms are CO
!! reverse the velocity in z-direction when CO reach to the upper limit
!! updated 2024/6/10, chenjun
!! check each CO (36 times)

  !! 1 ML    = 36CO = 72 atoms
  !! 3/4 ML  = 27CO = 54 atoms
  !! 1/2 ML  = 18CO = 36 atoms
  !! 1/3 ML  = 12CO = 24 atoms
  !! 1/4 ML  =  9CO = 18 atoms
 !if(natom.eq.72 &
 !   .or. natom.eq.54 .or. natom.eq.36 .or. natom.eq.24 .or. natom.eq.18 &
 !   .or. natom.eq.144 .or. natom.eq.56) then
  do iii=1,natom/2

  if ( Cmat(3,iii*2-1) .gt. zupper .and. Vmat(3,iii*2-1).gt.0.d0) then
  ! flag_reverse = .true.
    Vmat(3,iii*2-1) = -Vmat(3,iii*2-1) ! C
    Vmat(3,iii*2) = -Vmat(3,iii*2)     ! O
    write(uout,*) "Desorption (ps) / i-th CO ",i_step*timestep,iii
  endif

  enddo
 !endif

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


