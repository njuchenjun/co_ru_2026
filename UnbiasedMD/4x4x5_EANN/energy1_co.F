! --------------------------------------------------------------------
!     Q:    CARTESIAN IN ANGSTROM
!     VOUT: ENERGY IN EV
!     DVDQ: DERIVATIVES IN EV/ANGSTROM
! --------------------------------------------------------------------
#     DEFINE NATOM 2
#     DEFINE NBOND (NATOM*(NATOM-1)/2)
! --------------------------------------------------------------------
      SUBROUTINE ENERGY1_CO(Q,VOUT,DVDQ)
      IMPLICIT NONE
      REAL*8,INTENT(IN) :: Q(3,NATOM)
      REAL*8,INTENT(OUT) :: VOUT,DVDQ(3,NATOM)

      INTEGER I,J,K
      REAL*8 :: R(NBOND),DR(3,NATOM,NBOND),DVDR(NBOND)
      INTEGER,SAVE :: MB(NBOND),NB(NBOND)
      INTEGER,SAVE :: INIT=0

      REAL*8,PARAMETER :: TOANG=0.5291772083D0
      REAL*8,PARAMETER :: TOEV=27.21138505D0
      real*8 :: v01,dvco

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

      CALL INTERN_CO(Q/TOANG,MB,NB,R,DR) ! original Q in Angstrom, R in bohr

      !! R: bohr
      call potco(r(1),v01,dvdr) ! v01 in cm-1, dvdr in cm-1/bohr

      ! dvdr: cm-1/bohr
      CALL DVDX_CO(MB,NB,DR,DVDR,DVDQ)

      VOUT=v01/219474.5D0*toev ! now VOUT in eV
      DVDQ=DVDQ/219474.5D0*toev/TOANG ! now in eV/Angstrom

      RETURN
      END
! --------------------------------------------------------------------
      SUBROUTINE INTERN_CO(Q,MB,NB,R,DR)
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
      END SUBROUTINE
! --------------------------------------------------------------------
      SUBROUTINE DVDX_CO(MB,NB,DR,DVR,DV)
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
      END SUBROUTINE
! --------------------------------------------------------------------
