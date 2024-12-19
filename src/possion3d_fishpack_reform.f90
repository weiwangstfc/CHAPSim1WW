MODULE FISHPACK_POIS3D
  USE WPRECISION
  
  INTEGER(4), save :: LPEROD, MPEROD, NPEROD
  INTEGER(4), save :: LP, MP, NP
  INTEGER(4), save :: L, M, NG, NL, ML
  INTEGER(4), save :: LDIMF, MDIMF
  INTEGER(4), save :: WSZ

  REAL(WP), save :: C1, C2
  REAL(WP), save :: SCALX, SCALY
  REAL(WP), ALLOCATABLE, save :: A(:), B(:), C(:), D(:), BB(:)
  REAL(WP), ALLOCATABLE, save :: XRT(:), YRT(:), WX(:), WY(:)
  REAL(WP), ALLOCATABLE, save :: FR(:, :, :), FK(:, :, :), T(:)
  REAL(WP), ALLOCATABLE, save :: F_io   (:, :, :)

  private
  private :: FFTPACK_ROOT
  private :: TRID0
  private :: FFTPACK_XZ
  public :: FISHPACK_POIS3D_INIT
  public :: FISHPACK_POIS3D_SIMPLE

CONTAINS
!==========================================================================================================
!==========================================================================================================
  SUBROUTINE TRID0(NG, A, BB, C, T)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NG
    REAL(WP), DIMENSION(NG), INTENT(IN)    :: A
    REAL(WP), DIMENSION(NG), INTENT(IN)    :: BB
    REAL(WP), DIMENSION(NG), INTENT(IN)    :: C
    REAL(WP), DIMENSION(NG), INTENT(INOUT) :: T
    

    INTEGER :: NR, MM1, I, IP
    REAL(WP) :: Z
    REAL(WP), DIMENSION(NG) :: D

    NR = NG
    MM1 = NR - 1
    Z = 1.0_WP / BB(1)
    D(1) = C(1) * Z
    T(1) = T(1) * Z
    DO I = 2, MM1
        Z = 1.0_WP / (BB(I) - A(I) * D(I - 1))
        D(I) = C(I) * Z
        T(I) = (T(I) - A(I) * T(I - 1)) * Z
    END DO
    Z = BB(NR) - A(NR) * D(MM1)
    IF (DABS(Z) > 1.0E-14_WP) THEN
      T(NR) = (T(NR) - A(NR) * T(MM1)) /Z
    ELSE
      T(NR) = 0.0_WP
    END IF
    DO IP = 1, MM1
      I = NR-IP
      T(I) = T(I) - D(I) * T(I + 1)
    END DO

    RETURN
  END SUBROUTINE

!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FFTPACK_ROOT_1D(L, LP, C1, SCALX, WX, XRT)
  ! Generate FFT transform roots for 1-D direction
  ! Arguments:
  !   L    (INTEGER(4), IN): Grid dimensions in X direction
  !   LP  (INTEGER(4), IN): Parameters determining transform type
  !   C1  (REAL(WP), IN): Scaling factors for X transform
  !   SCALX(REAL(WP), IN): scaling factor 
  !   WX  (REAL(WP), DIMENSION(:), INOUT): Work arrays for FFT
  !   XRT (REAL(WP), DIMENSION(:), OUT): Transform roots for X

    IMPLICIT NONE

    ! Arguments
    INTEGER(4), INTENT(IN) :: L, LP
    REAL(WP), INTENT(IN) :: C1
    REAL(WP), INTENT(OUT) :: SCALX
    REAL(WP), DIMENSION(:), INTENT(OUT) :: WX
    REAL(WP), DIMENSION(:), INTENT(OUT) :: XRT
    

    ! Local variables
    REAL(WP) :: PI, DX, DI
    INTEGER(4) :: LR, LRDEL, I

    ! Compute PI
    PI = 2.0_WP * DASIN(1.0_WP)
    !-----------------------------------------------------------
    ! X direction transform roots
    !-----------------------------------------------------------
    LR = L
    LRDEL = ((LP - 1) * (LP - 3) * (LP - 5)) / 3
    SCALX = DBLE(LR + LRDEL)
    DX = PI / (2.0_WP * SCALX)

    SELECT CASE (LP)
    CASE (1)
      ! RFFTI     INITIALIZE  RFFTF AND RFFTB
      ! RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
      ! RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
      XRT(1) = 0.0_WP
      XRT(LR) = -4.0_WP * C1
      DO I = 3, LR, 2
        XRT(I - 1) = -4.0_WP * C1 * (DSIN(DBLE((I - 1)) * DX))**2
        XRT(I) = XRT(I - 1)
      END DO
      CALL RFFTI(LR, WX)

    CASE (2)
      ! SINTI     INITIALIZE SINT
      ! SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
      DI = 0.00_WP
      DO I = 1, LR
        XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
      END DO
      SCALX = 2.0_WP * SCALX
      CALL SINTI (LR, WX)

    CASE (3)
    ! SINQI     INITIALIZE SINQF AND SINQB
    ! SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
    ! SINQB     UNNORMALIZED INVERSE OF SINQF
      DI = 0.50_WP
      SCALX = 2.0_WP * SCALX
      DO I = 1, LR
        XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
      ENDDO
      SCALX = 2.0_WP * SCALX
      CALL SINQI (LR, WX)

    CASE (4)
    ! COSTI     INITIALIZE COST
    ! COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
      DI = 1.00_WP
      DO I = 1, LR
          XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
      END DO
      SCALX = 2.0_WP * SCALX
      CALL COSTI (LR, WX)

    CASE (5)
    ! COSQI     INITIALIZE COSQF AND COSQB
    ! COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
    ! COSQB     UNNORMALIZED INVERSE OF COSQF
      DI = 0.50_WP
      SCALX = 2.0_WP * SCALX
      DO I = 1, LR
        XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
      END DO
      SCALX = 2.0_WP * SCALX
      CALL COSQI (LR, WX)

    END SELECT

    RETURN
  END SUBROUTINE 

!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FFTPACK_X1D(IFWRD, MR, NR, LR, LP, WX, FR)
  ! FOR FFT TRANSFORM FROM SPACE TO WAVENUMBER, IFWRD = 1, IS = 1
  ! FOR FFT TRANSFORM FROM WAVENUMBER TO SPACE, IFWRD = 2, IS = -1
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: IFWRD
    INTEGER(4), INTENT(IN) :: MR, NR, LR, LP
    REAL(WP), DIMENSION(:),       INTENT(IN)    :: WX
    REAL(WP), DIMENSION(:, :, :), INTENT(INOUT) :: FR

    REAL(WP) :: T(LR)
    INTEGER(4) :: I, J, K
    
    IF(IFWRD == 1) then
      DO J = 1, MR
        DO k = 1, NR
          DO I = 1, LR  
              T(I) = FR(I, J, K)
          END DO
          SELECT CASE (LP)
          CASE (1)
            CALL RFFTF (LR, T, WX)
          CASE (2)
            CALL SINT (LR, T, WX)
          CASE (3)
            CALL SINQF (LR, T, WX)
          CASE (4)
            CALL COST (LR, T, WX)
          CASE (5)
            CALL COSQF (LR, T, WX)
          END SELECT
          DO I = 1, LR
              FR(I, J, K) = T(I)
          END DO
        ENDDO
      END DO

    ELSE IF (IFWRD == 2) then

      DO J = 1, MR
        DO k = 1, NR
          DO I = 1, LR  
              T(I) = FR(I, J, K)
          END DO
          SELECT CASE (LP)
          CASE (1)
            CALL RFFTB (LR, T, WX)
          CASE (2)
            CALL SINT (LR, T, WX)
          CASE (3)
            CALL SINQB (LR, T, WX)
          CASE (4)
            CALL COST (LR, T, WX)
          CASE (5)
            CALL COSQB (LR, T, WX)
          END SELECT
          DO I = 1, LR
              FR(I, J, K) = T(I)
          END DO
        ENDDO
      END DO

    END IF

  END SUBROUTINE
!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FFTPACK_Z1D(IFWRD, MR, NR, LR, MP, WY, FR)
    ! FOR FFT TRANSFORM FROM SPACE TO WAVENUMBER, IFWRD = 1, IS = 1
  ! FOR FFT TRANSFORM FROM WAVENUMBER TO SPACE, IFWRD = 2, IS = -1
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: IFWRD
    INTEGER(4), INTENT(IN) :: MR, NR, LR, MP
    REAL(WP), DIMENSION(:),       INTENT(IN)    :: WY
    REAL(WP), DIMENSION(:, :, :), INTENT(INOUT) :: FR

    REAL(WP) :: T(MR)
    INTEGER(4) :: I, J, K
  
    IF(IFWRD == 1) then
      DO I = 1, LR
          DO k = 1, NR
            DO  J = 1, MR
                T(J) = FR(I, J, K)
            END DO
            SELECT CASE (MP)
            CASE (1)
              CALL RFFTF (MR, T, WY)
            CASE (2)
              CALL SINT (MR, T, WY)
            CASE (3)
              CALL SINQF (MR, T, WY)
            CASE (4)
              CALL COST (MR, T, WY)
            CASE (5)
              CALL COSQF (MR, T, WY)
            END SELECT
            DO J = 1, MR
                FR(I, J, K) = T(J)
            END DO
          END DO
        END DO

      ELSE IF (IFWRD == 2) THEN

        DO I = 1, LR
          DO k = 1, NR
            DO  J = 1, MR
                T(J) = FR(I, J, K)
            END DO
            SELECT CASE (MP)
            CASE (1)
              CALL RFFTB (MR, T, WY)
            CASE (2)
              CALL SINT (MR, T, WY)
            CASE (3)
              CALL SINQB (MR, T, WY)
            CASE (4)
              CALL COST (MR, T, WY)
            CASE (5)
              CALL COSQB (MR, T, WY)
            END SELECT
            
            DO J = 1, MR
                FR(I, J, K) = T(J)
            END DO
          END DO
        END DO
      END IF
  
  
    END SUBROUTINE
!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FFTPACK_XZ2D(IFWRD, MR, NR, LR, LP, MP, WX, WY, FR)
  ! FOR FFT TRANSFORM FROM SPACE TO WAVENUMBER, IFWRD = 1, IS = 1
  ! FOR FFT TRANSFORM FROM WAVENUMBER TO SPACE, IFWRD = 2, IS = -1
    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: IFWRD
    INTEGER(4), INTENT(IN) :: MR, NR, LR, MP, LP
    REAL(WP), DIMENSION(:),       INTENT(IN)    :: WX, WY
    REAL(WP), DIMENSION(:, :, :), INTENT(INOUT) :: FR


    IF(IFWRD == 1) then
      !     TRANSFORM X FORWARD: PHY -> SPEC
      CALL FFTPACK_X1D(IFWRD, MR, NR, LR, LP, WX, FR)
      !     TRANSFORM Z FORWARD: PHY -> SPEC
      CALL FFTPACK_Z1D(IFWRD, MR, NR, LR, MP, WY, FR)

    ELSE IF (IFWRD == 2) then
      !     TRANSFORM Z BACKWARD: SPEC -> PHY
      CALL FFTPACK_Z1D(IFWRD, MR, NR, LR, MP, WY, FR)
      !     TRANSFORM X BACKWARD: SPEC -> PHY
      CALL FFTPACK_X1D(IFWRD, MR, NR, LR, LP, WX, FR)
    END IF

  END SUBROUTINE


!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FISHPACK_POIS3D_INIT(BCX_io, BCZ, NCL1_io, NCL2, NCL3, N2DO, N3DO, DXQI, DZQI)
    !    SHOULD BE CALLED IN SLAVES BEFORE FIRSTLY
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: BCX_io(2)
    INTEGER, INTENT(IN) :: BCZ(2)
    INTEGER, INTENT(IN) :: NCL1_io, NCL2, NCL3, N2DO, N3DO
    REAL(WP), INTENT(IN) :: DXQI, DZQI
    
    INTEGER, INTENT(OUT) :: LP, MP, NP
    INTEGER, INTENT(OUT) :: L, M, NG, NL, ML
    REAL(WP), INTENT(OUT) :: C1, C2


    INTEGER(4) :: K
    INTEGER :: LPEROD, MPEROD, NPEROD

    IF(BCX_io(1) == 3 .AND. BCX_io(2) == 3 )  LPEROD = 0
    IF(BCX_io(1) == 1 .AND. BCX_io(2) == 1 )  LPEROD = 1
    IF(BCX_io(1) == 1 .AND. BCX_io(2) == 2 )  LPEROD = 2
    IF(BCX_io(1) == 2 .AND. BCX_io(2) == 2 )  LPEROD = 3
    IF(BCX_io(1) == 2 .AND. BCX_io(2) == 1 )  LPEROD = 4

    IF(BCZ(1) == 3 .AND. BCZ(2) == 3 )  MPEROD = 0
    IF(BCZ(1) == 1 .AND. BCZ(2) == 1 )  MPEROD = 1
    IF(BCZ(1) == 1 .AND. BCZ(2) == 2 )  MPEROD = 2
    IF(BCZ(1) == 2 .AND. BCZ(2) == 2 )  MPEROD = 3
    IF(BCZ(1) == 2 .AND. BCZ(2) == 1 )  MPEROD = 4

    NPEROD = 1

    LP = LPEROD + 1
    MP = MPEROD + 1
    NP = NPEROD + 1

    L  = NCL1_io
    M  = NCL3
    NG = NCL2
    NL = N2DO!(MYID)
    ML = N3DO!(MYID)

    IF( (L <= 3) .OR. (M <= 3) .OR. (NG <= 3) ) &
    CALL ERRHDL('Dimensions in poISson solver should be lARger than 3!', MYID)

    C1 = DXQI
    C2 = DZQI

    WSZ = 30 + L + M + 2 * NG + MAX(L, M, NG) + &
                7*(INT((L+1)/2) + INT((M+1)/2)) + 128

    ALLOCATE (A (NG) ) ; A = 0.0_WP
    ALLOCATE (B (NG) ) ; B = 0.0_WP
    ALLOCATE (C (NG) ) ; C = 0.0_WP
    ALLOCATE (D (NG) ) ; D = 0.0_WP
    ALLOCATE (BB(NG) ) ; BB = 0.0_WP

    MEMPC_Byte = MEMPC_Byte + NG*10 *8

    ALLOCATE (XRT(L))  ; XRT = 0.0_WP
    ALLOCATE (YRT(M))  ; YRT = 0.0_WP
    ALLOCATE (WX(WSZ)) ; WX = 0.0_WP
    ALLOCATE (WY(WSZ)) ; WY = 0.0_WP

    ALLOCATE (FR(L, M, NL) ) ; FR = 0.0_WP
    ALLOCATE (FK(L, ML, NG) ) ; FK = 0.0_WP
    ALLOCATE (T(MAX0(L, M, NG)) ) ; T = 0.0_WP

    ALLOCATE ( F_io   (NCL1_io, NCL2, N3DO(0) )     )       ;  F_io = 0.0_WP

    MEMPC_Byte = MEMPC_Byte + (L+M + WSZ*2 +L * M * NL+L * ML * NG+MAX0(L, M, NG)) *8

    DO K = 1, NG
        A(K) = AMPH(K)
        B(K) = ACPH(K)
        C(K) = APPH(K)
    END DO

    ! Build up fft root for x
    CALL FFTPACK_ROOT_1D(L, LP, C1, SCALX, WX, XRT)
    ! Build up fft root for y (z in this code)
    CALL FFTPACK_ROOT_1D(M, MP, C2, SCALY, WY, YRT)

    RETURN

  END SUBROUTINE


!==========================================================================================================
!==========================================================================================================
  SUBROUTINE FISHPACK_POIS3D_SIMPLE(RHS, VAR)
    !    CALLED IN SLAVES EVERY RK STAGE TO CALCULATE pressure CORRECTION TERMS.
    USE mesh_info
    USE init_info
    IMPLICIT NONE
    REAL(WP), INTENT(INOUT) :: RHS (NCL1_io, N2DO(0), NCL3 )
    REAL(WP), INTENT(OUT) :: VAR (NCL1_io, N2DO(0), NCL3 )
    INTEGER(4) :: IFWRD
    INTEGER(4) :: I, J, K, JJ

    !   ======== RECONSTRUCT RHS TO FIT POIS3D.========
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
                FR(I, K, J) = RHS(I, J, K)
            END DO
        END DO
    END DO

    !   ========forward FFT IN X AND Z direction ========
    IFWRD = 1
    CALL FFTPACK_XZ(IFWRD)

    !   ======== RESTORE RHSLLPHI ========
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
              RHS(I, J, K) = FR(I, K, J)
            END DO
        END DO
    END DO

    !   ======== TRANSPORT Y- DECOMP TO K - DECOMP ========
    !CALL TRASP_Y2Z_RHSLLPHI_io
    CALL TRASP23_Y2Z(NCL1_io, 1, N2DO(0), RHS, F_io)

    !   ========CONSTRUCT DATA FOR TDMA========
    DO I = 1, L
        DO K = 1, N3DO(MYID)
            DO J = 1, NG
                FK(I, K, J) = F_io(I, J, K)
            END DO
        END DO
    END DO

    !   ======== TDMA IN Y direction FOR PART OF FR(:,PART, :) ========
    DO I = 1, L
        DO J = 1, N3DO(MYID)
            JJ = KCL2G(J)
            DO K = 1, NG
                BB(K) = B(K) + XRT(I) / RCCI2(K) + YRT(JJ)
                T(K) = FK(I, J, K)
            END DO

            CALL TRID0

            DO K = 1, NG
                FK(I, J, K) = T(K)
            END DO

        END DO
    END DO

    !      ======== RE -CONSTRUCT DATA BACK ========
    DO I = 1, L
        DO K = 1, N3DO(MYID)
            DO J = 1, NG
                F_io(I, J, K) = FK(I, K, J)
            END DO
        END DO
    END DO

    !   ======== TRANSPORT Z - DECOMP TO Y- DECOMP ========
    !CALL TRASP_Z2Y_RHSLLPHI_io
    CALL TRASP23_Z2Y(NCL1_io, 1, N2DO(0), RHS, F_io)

    !   ======== RESTORE RHSLLPHI ========
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
                FR(I, K, J) = RHS(I, J, K)
            END DO
        END DO
    END DO

    !   ========backward FFT IN Z AND X directionS ========
    IFWRD = 2
    CALL FFTPACK_XZ(IFWRD)

    !   ======== SCALE THE CALCULATED VALUE ========
    DO I = 1, L
        DO J = 1, M
            DO K = 1, NL
                FR(I, J, K) = FR(I, J, K) / (SCALX * SCALY)
            END DO
        END DO
    END DO


    !   ======== RE -STORE AND ASSIGN DATA TO DPH ========
    VAR = 0.0_WP
    DO I = 1, L
        DO K = 1, M
            DO J = 1, NL
              VAR(I, J, K) = FR(I, K, J)
            END DO
        END DO
    END DO

    RETURN

  END SUBROUTINE

END MODULE

































SUBROUTINE FFTPACK_ROOT(L, M, LP, MP, C1, C2, WX, WY, XRT, YRT)
  ! Generate FFT transform roots for X and Y directions
  ! Arguments:
  !   L, M    (INTEGER(4), IN): Grid dimensions in X and Y directions
  !   LP, MP  (INTEGER(4), IN): Parameters determining transform type
  !   C1, C2  (REAL(WP), IN): Scaling factors for X and Y transforms
  !   WX, WY  (REAL(WP), DIMENSION(:), INOUT): Work arrays for FFT
  !   XRT, YRT (REAL(WP), DIMENSION(:), OUT): Transform roots for X and Y

  IMPLICIT NONE

  ! Arguments
  INTEGER(4), INTENT(IN) :: L, M, LP, MP
  REAL(WP), INTENT(IN) :: C1, C2
  REAL(WP), DIMENSION(:), INTENT(INOUT) :: WX, WY
  REAL(WP), DIMENSION(:), INTENT(OUT) :: XRT, YRT

  ! Local variables
  REAL(WP) :: PI, DX, DI, DY, DJ, SCALX, SCALY
  INTEGER(4) :: LR, MR, LRDEL, MRDEL, I, J

  ! Compute PI
  PI = 2.0_WP * DASIN(1.0_WP)

  !-----------------------------------------------------------
  ! X direction transform roots
  !-----------------------------------------------------------
  LR = L
  LRDEL = ((LP - 1) * (LP - 3) * (LP - 5)) / 3
  SCALX = DBLE(LR + LRDEL)
  DX = PI / (2.0_WP * SCALX)

  SELECT CASE (LP)
  CASE (1)
    ! RFFTI     INITIALIZE  RFFTF AND RFFTB
    ! RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
    ! RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
    XRT(1) = 0.0_WP
    XRT(LR) = -4.0_WP * C1
    DO I = 3, LR, 2
      XRT(I - 1) = -4.0_WP * C1 * (DSIN(DBLE((I - 1)) * DX))**2
      XRT(I) = XRT(I - 1)
    END DO
    CALL RFFTI(LR, WX)

  CASE (2)
    ! SINTI     INITIALIZE SINT
    ! SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
    DI = 0.00_WP
    DO I = 1, LR
      XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
    END DO
    SCALX = 2.0_WP * SCALX
    CALL SINTI (LR, WX)

  CASE (3)
   ! SINQI     INITIALIZE SINQF AND SINQB
   ! SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
   ! SINQB     UNNORMALIZED INVERSE OF SINQF
    DI = 0.50_WP
    SCALX = 2.0_WP * SCALX
    DO I = 1, LR
      XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
    ENDDO
    SCALX = 2.0_WP * SCALX
    CALL SINQI (LR, WX)

  CASE (4)
   ! COSTI     INITIALIZE COST
   ! COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
    DI = 1.00_WP
    DO I = 1, LR
        XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
    END DO
    SCALX = 2.0_WP * SCALX
    CALL COSTI (LR, WX)

  CASE (5)
   ! COSQI     INITIALIZE COSQF AND COSQB
   ! COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
   ! COSQB     UNNORMALIZED INVERSE OF COSQF
    DI = 0.50_WP
    SCALX = 2.0_WP * SCALX
    DO I = 1, LR
      XRT(I) = -4.0_WP * C1 * (DSIN((DBLE(I) - DI) * DX))**2
    END DO
    SCALX = 2.0_WP * SCALX
    CALL COSQI (LR, WX)

  END SELECT

  !-----------------------------------------------------------
  ! Y direction transform roots
  !-----------------------------------------------------------
  MR = M
  MRDEL = ((MP - 1) * (MP - 3) * (MP - 5)) / 3
  SCALY = DBLE(MR + MRDEL)
  DY = PI / (2.0_WP * SCALY)

  SELECT CASE (MP)
  CASE (1)
    YRT(1) = 0.0_WP
    YRT(MR) = -4.0_WP * C2
    DO J = 3, MR, 2
      YRT(J - 1) = -4.0_WP * C2 * (DSIN(DBLE((J - 1)) * DY))**2
      YRT(J) = YRT(J - 1)
    END DO
    CALL RFFTI (MR, WY)

  CASE (2)
    DJ = 0.00_WP
    DO J = 1, MR
      YRT(J) = -4.0_WP * C2 * (DSIN((DBLE(J) - DJ) * DY))**2
    ENDDO
    SCALY = 2.0_WP * SCALY
    CALL SINTI (MR, WY)

  CASE (3)
    DJ = 0.50_WP
    SCALY = 2.0_WP * SCALY
    DO J = 1, MR
      YRT(J) = -4.0_WP * C2 * (DSIN((DBLE(J) - DJ) * DY))**2
    ENDDO
    SCALY = 2.0_WP * SCALY
    CALL SINQI (MR, WY)

  CASE (4)
    DJ = 1.00_WP
    DO J = 1, MR
      YRT(J) = -4.0_WP * C2 * (DSIN((DBLE(J) - DJ) * DY))**2
    ENDDO
    SCALY = 2.0_WP * SCALY
    CALL COSTI (MR, WY)

  CASE (5)
    DJ = 0.50_WP
    SCALY = 2.0_WP * SCALY
    DO J = 1, MR
      YRT(J) = -4.0_WP * C2 * (DSIN((DBLE(J) - DJ) * DY))**2
    ENDDO
    SCALY = 2.0_WP * SCALY
    CALL COSQI (MR, WY)

  END SELECT


  RETURN
END SUBROUTINE FFTPACK_ROOT