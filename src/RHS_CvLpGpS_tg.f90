!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Kui He, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> to calculate RHS of Eq.(A6)
!
!>@param      QC  [in]  H11 etc - \partial{u_1 u_j}/\partial{x_j}
!>@param      AL  [in]  alpha coefficients in Eq.A4a of Mehdi paper.
!>@param     GAL  [in]  gamma coefficients in Eq.A4a of Mehdi paper.
!>@param     ROL  [in]  zeta  coefficients in Eq.A4a of Mehdi paper.
!>@param      PR  [in]  pressure
!>@param      NQ  [in]  =IVEL, velocity direction 1=x, 2=y, 3=z. 
!>@param     NIN  [in]  =initial index in the direction of x,y,z which are 1,2,1
!>@param     RHS  [out] RHS of Eq.(A6) 
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! ??/??/201?- Add pipe treatment, by Kui He.
! 06/02/2014- Subroutine structured is optimized by Wei Wang
! 14/02/2014- Remove the du/dy at low wall for pipe and commented the wall
!             shear stress calcualation. by Wei Wang
!**********************************************************************
    SUBROUTINE RHS_CvLpGpS_tg(NS,IDR)
        use init_info
        use mesh_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE

        INTEGER(4),INTENT(IN) :: NS
        INTEGER(4),INTENT(IN) :: IDR

        REAL(WP) :: INTGRHSY
        REAL(WP) :: INTGRHSY_WORK
        !REAL(WP) :: INTGU
        !REAL(WP) :: INTGU_WORK
        REAL(WP) :: DPGRNS
        !REAL(WP) :: INTGVOL
        !REAL(WP) :: INTGVOL_WORK       
        !REAL(WP) :: RHSCCC
        REAL(WP) :: SUCACJ
        REAL(WP) :: DERQ
        REAL(WP) :: RHSC, RHSL
        !REAL(WP) :: VOLtmp
        REAL(WP) :: PGM      
        INTEGER(4) :: NII
        INTEGER(4) :: I, IC, IM, IP
        INTEGER(4) :: J, JC, JM, JP, JJ
        INTEGER(4) :: K, KC, KM, KP
        REAL(WP)    :: RMC2(N2DO(0))
        REAL(WP)    :: CONVH_tg(NCL1_tg,N2DO(0),NCL3)  
        REAL(WP)    :: COE1,COE2
 
!>      @note Setup the initial y values 
        NII=1
        IF((IDR.EQ.2) .and. (MYID.EQ.0)) NII=2

        RMC2 = 0.0_WP
        IF(IDR.EQ.2)THEN
            DO JC=NII, N2DO(MYID)
                JJ=JCL2G(JC)
                RMC2(JC)=RNDI2(JJ)
            END DO
        ELSE
            DO JC=NII, N2DO(MYID)
                JJ=JCL2G(JC)
                RMC2(JC)=RCCI2(JJ)
            END DO
        END IF  
      
        !======================construct the convectiont term===============================
        CONVH_tg = 0.0_WP
        IF(IDR.EQ.1) THEN
      
            DO I=1,NCL1_tg
                DO J=NII,N2DO(MYID)
                    DO K=1,NCL3
                        CONVH_tg(I,J,K) = Qtmp_tg(I,J,K)
                    END DO
                END DO
            END DO
      
        ELSE IF(IDR .EQ. 2) THEN
      
            DO I=1,NCL1_tg
                DO J=NII,N2DO(MYID)
                    DO K=1,NCL3
                        CONVH_tg(I,J,K) = DPH_tg(I,J,K)
                    END DO
                END DO
            END DO

        ELSE IF(IDR .EQ. 3) THEN
      
            DO I=1,NCL1_tg
                DO J=NII,N2DO(MYID)
                    DO K=1,NCL3
                        CONVH_tg(I,J,K) = RHSLLPHI_tg(I,J,K)
                    END DO
                END DO
            END DO
      
        ELSE
        END IF
             
        !======================construct viscous term and summing up the convection/viscous terms===============================
        COE1 = TALP(NS)*CVISC
        RHSC = 0.0_WP
        RHSL = 0.0_WP
        
        DO KC=1,NCL3
            KM=KMV(KC)
            KP=KPV(KC)
            DO JC=NII,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                DO IC=1,NCL1_tg
                    IP=IPV_tg(IC)
                    IM=IMV_tg(IC)
                    RHSC  = TGAM(NS)*CONVH_tg(IC,JC,KC)+TROH(NS)*CONVH0_tg(IC,JC,KC,IDR)
                    RHSL  = COE1*( (Q_tg(IP,JC,KC,IDR)-2.0_WP*Q_tg(IC,JC,KC,IDR)+Q_tg(IM,JC,KC,IDR))*DXQI +              &
                                   (Q_tg(IC,JC,KP,IDR)-2.0_WP*Q_tg(IC,JC,KC,IDR)+Q_tg(IC,JC,KM,IDR))*DZQI*RMC2(JC)+      &
                                   (Q_tg(IC,JP,KC,IDR)*APVR(JJ,IDR)+             &
                                    Q_tg(IC,JC,KC,IDR)*ACVR(JJ,IDR)+             &
                                    Q_tg(IC,JM,KC,IDR)*AMVR(JJ,IDR) )            &
                                 )
                    CONVH0_tg(IC,JC,KC,IDR)=CONVH_tg(IC,JC,KC)
                    RHS_tg(IC,JC,KC)=(RHSC+RHSL)*DT
                ENDDO
            ENDDO
        ENDDO
            
        !=========================pressure gradient terms==================================================
        PGM=0.0_WP
        COE2 = TALP(NS)*DT
        IF (IDR.EQ.1) THEN 
            DO KC=1,NCL3
                DO JC=NII,N2DO(MYID)
                    DO IC=1,NCL1_tg
                        IM=IMV_tg(IC)
                        PGM = (PR_tg(IC,JC,KC)-PR_tg(IM,JC,KC))*DXI*COE2 
                        RHS_tg(IC,JC,KC)=RHS_tg(IC,JC,KC)-PGM
                    ENDDO
                ENDDO
            ENDDO
        ELSE IF (IDR.EQ.2) THEN 
            DO KC=1,NCL3
                DO JC=NII,N2DO(MYID)
                    JM=JLMV(JC)
                    JJ=JCL2G(JC)
                    SUCACJ=DYCI(JJ)/RNDI1(jj)
                    DO IC=1,NCL1_tg
                        PGM = (PR_tg(IC,JC,KC)-PR_tg(IC,JM,KC))*SUCACJ*COE2
                        RHS_tg(IC,JC,KC)=RHS_tg(IC,JC,KC)-PGM
                    ENDDO
                ENDDO
            ENDDO      
        ELSE IF (IDR.EQ.3) THEN
            DO KC=1,NCL3
                KM=KMV(KC)
                DO JC=NII,N2DO(MYID)
                    DO IC=1,NCL1_tg
                        PGM = (PR_tg(IC,JC,KC)-PR_tg(IC,JC,KM))*DZI*COE2
                        RHS_tg(IC,JC,KC)=RHS_tg(IC,JC,KC)-PGM
                    ENDDO
                ENDDO
            ENDDO  
        ELSE       
        ENDIF
    
        !====================flow drive terms (source terms) in periodic streamwise flow===========================
        DPGRNS=0.0_WP
        DERQ=0.0_WP
        IF (IDR.EQ.NFLOW) THEN
            !=============constant mass flow rate=============================
            IF(FLOWDRVTP==1) THEN         
                INTGRHSY=0.0_WP

                DO IC=1,NCL1_tg
                    DO JC=1,N2DO(MYID)
                        DO KC=1,NCL3
                            JJ=JCL2G(JC)
                            INTGRHSY=INTGRHSY+RHS_tg(IC,JC,KC)/DYFI(JJ)/RCCI1(JJ)
                        ENDDO
                    ENDDO
                ENDDO
                
                CALL MPI_BARRIER(ICOMM,IERROR)
                CALL MPI_ALLREDUCE(INTGRHSY,INTGRHSY_WORK,1,MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
                IF(MYID.EQ.0) DPGRNS=INTGRHSY_WORK/VOLM_tg   
                
                
            ELSE IF(FLOWDRVTP==2) THEN 
                IF(MYID.EQ.0) DPGRNS = -0.50_WP*CFGV*COE2 !dimensionless based on \delta and U_m   
                !constant pressure gradient, dimensionless based on \delta and U_m
                !DPGRNS = -2.0_WP      ! DIMENSIONLESS BASED ON U_TAU
            ELSE
            END IF   
            
            CALL MPI_BCAST( DPGRNS, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

            DO K=1,NCL3
                DO I=1,NCL1_tg
                    DO J=NII,N2DO(MYID)
                        RHS_tg(I,J,K)=RHS_tg(I,J,K) - DPGRNS
                    ENDDO
                ENDDO
            ENDDO

        END IF 

        RETURN
     END SUBROUTINE
