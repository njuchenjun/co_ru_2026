! --------------------------------------------------------------------
!     Q:    CARTESIAN IN ANGSTROM
!     VOUT: ENERGY IN eV (Interaction energy only)
!     DVDQ: DERIVATIVES IN eV/Angstrom
! --------------------------------------------------------------------
#     DEFINE NATOM 4
#     DEFINE NBOND (NATOM*(NATOM-1)/2)
! --------------------------------------------------------------------
      SUBROUTINE ENERGY1(Q,VOUT,DVDQ)
      IMPLICIT NONE
      REAL*8,INTENT(IN) :: Q(3,NATOM)
      REAL*8,INTENT(OUT) :: VOUT,DVDQ(3,NATOM)

      INTEGER I,J,K
      REAL*8 :: R(NBOND),DR(3,NATOM,NBOND),DVDR(NBOND)
      INTEGER,SAVE :: MB(NBOND),NB(NBOND)
      INTEGER,SAVE :: INIT=0

      REAL*8,PARAMETER :: EV=27.21138505D0
      REAL*8,PARAMETER :: TOEE=27.21138505D0*23.0605D0*0.04184D0
      REAL*8,PARAMETER :: TOANG=0.5291772083D0

      real*8 :: v12,v01,v02,dvco1,dvco2
      real*8 :: dist,weight

      IF(INIT.EQ.0) THEN
        INIT=1
        K=1
        DO I=1,NATOM-1
        DO J=I+1,NATOM
          MB(K)=I
          NB(K)=J
          K=K+1
        ENDDO
        ENDDO
      ENDIF

      !! cut-off at 10 to 12 Angstrom, using distances of C--C
      !! 2024/6/10
      dist=dsqrt(dot_product(Q(1:3,2)-Q(1:3,4),Q(1:3,2)-Q(1:3,4)))
      weight=0.d0

      if(dist.gt.12.d0) then
        vout=0.d0
        dvdq=0.d0
        weight=0.d0
        return
      elseif(dist.le.10.d0) then
        weight=1.d0
      else
        weight=dcos((dist-10.d0)/2.d0*3.14159265358979d0)/2.d0+0.5d0
      endif
      !! end of cut-off definition

      CALL INTERN(Q,MB,NB,R,DR)

      !! R: O1 C1 O2 C2 Angstrom
      call pes_ococ(r/toang,v12,1,dvdr)
!     call potco(r(1)/toang,v01,dvco1)
!     call potco(r(6)/toang,v02,dvco2)
!     vout=(v12+v01+v02)/219474.5d0*toee
      vout=v12/219474.5d0*EV*weight

      ! dvdr: cm-1/bohr
!     dvdr(1)=dvdr(1)+dvco1
!     dvdr(6)=dvdr(6)+dvco2
!     dvdr=dvdr/219474.5d0*toee/TOANG

      dvdr=dvdr/219474.5d0*EV/TOANG*weight ! now in eV/Angstrom

      CALL DVDX(MB,NB,DR,DVDR,DVDQ)

!     print*,"test"
!     print*,Q
!     write(*,'(6f10.5)') r
!     write(*,'(6f10.5)') r/toang
!     write(*,*) vout
!     stop

      RETURN
      END
! --------------------------------------------------------------------
      SUBROUTINE INTERN(Q,MB,NB,R,DR)
! ---- CALCULATES BONDLENGTHS AND THEIR DERIVATIVES
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8,INTENT(IN) :: Q(3,NATOM)
      INTEGER,INTENT(IN) :: MB(NBOND),NB(NBOND)
      REAL*8,INTENT(OUT) :: R(NBOND),DR(3,NATOM,NBOND)
! ---- BOND LENGTHS
      DO I=1,NBOND
        R(I)=DSQRT(
     .   (Q(1,MB(I))-Q(1,NB(I)))*(Q(1,MB(I))-Q(1,NB(I)))
     .  +(Q(2,MB(I))-Q(2,NB(I)))*(Q(2,MB(I))-Q(2,NB(I)))
     .  +(Q(3,MB(I))-Q(3,NB(I)))*(Q(3,MB(I))-Q(3,NB(I))))
      ENDDO
! ---- BOND LENGTH DERIVATIVES
      DO I=1,NBOND
         DR(1,MB(I),I)=(Q(1,MB(I))-Q(1,NB(I)))/R(I)
         DR(2,MB(I),I)=(Q(2,MB(I))-Q(2,NB(I)))/R(I)
         DR(3,MB(I),I)=(Q(3,MB(I))-Q(3,NB(I)))/R(I)
         DR(1,NB(I),I)=-DR(1,MB(I),I)
         DR(2,NB(I),I)=-DR(2,MB(I),I)
         DR(3,NB(I),I)=-DR(3,MB(I),I)
      ENDDO
      RETURN
      END SUBROUTINE INTERN
! --------------------------------------------------------------------
      SUBROUTINE DVDX(MB,NB,DR,DVR,DV)
! ---- CALCULATES DV/DX FROM DV/DR*DR/DX
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MB(NBOND),NB(NBOND)
      REAL*8,INTENT(IN) :: DR(3,NATOM,NBOND),DVR(NBOND)
      REAL*8,INTENT(OUT) :: DV(3,NATOM)
      INTEGER :: I
      DV=0.D0
      DO I=1,NBOND
         DV(1,MB(I))=DV(1,MB(I))+DVR(I)*DR(1,MB(I),I)
         DV(1,NB(I))=DV(1,NB(I))+DVR(I)*DR(1,NB(I),I)
         DV(2,MB(I))=DV(2,MB(I))+DVR(I)*DR(2,MB(I),I)
         DV(2,NB(I))=DV(2,NB(I))+DVR(I)*DR(2,NB(I),I)
         DV(3,MB(I))=DV(3,MB(I))+DVR(I)*DR(3,MB(I),I)
         DV(3,NB(I))=DV(3,NB(I))+DVR(I)*DR(3,NB(I),I)
      ENDDO
      RETURN
      END SUBROUTINE DVDX
! --------------------------------------------------------------------
