!***********************************************************************
!> @author 
!> Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Factorization Approximation of the momentume equation.
!> @param AL  [in] alpha coefficients in Eq.A4a of Mehdi paper.
!> @param IDR [in] velocity IDRection. 
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! 06/02/2014- Initial Version, by Wei Wang
!**********************************************************************    
    SUBROUTINE MOMFA_io(NS,IDR)
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
 
        INTEGER(4),INTENT(IN)  :: NS
        INTEGER(4),INTENT(IN)  :: IDR
      
        INTEGER(4)  :: N2I
        INTEGER(4)  :: I, J, K
      
        FACOE_io=0.50_WP*TALP(NS)*DT*CVISC     !=alpha*dt/(2Re)
    
        N2I=1
        IF ((MYID.EQ.0).AND.(IDR.EQ.2).AND.(ICASE.NE.IBOX3P))  N2I=2

        CALL MOMFA1_X_io(N2I,IDR)
        CALL MOMFA2_Z_io(N2I,IDR)
        CALL MOMFA3_Y_io(IDR)
      
     
        DO K=1,NCL3
            DO J=N2I,N2DO(MYID)
                DO I=1,NCL1_io
                    G_io(I,J,K,IDR)=RHS_io(I,J,K) + G_io(I,J,K,IDR)
                ENDDO
            ENDDO
        ENDDO

        RETURN
    END SUBROUTINE
    
!********************************************************************************    
    SUBROUTINE MOMFA1_X_io(N2I,IDR)
        use mesh_info
        use flow_info
        use thermal_info
        use init_info
        IMPLICIT NONE
                  
        INTEGER(4),INTENT(IN)    :: N2I
        INTEGER(4),INTENT(IN)    :: IDR
      
        REAL(WP)  :: AMIV(NCL1_io,N2DO(0))
        REAL(WP)  :: ACIV(NCL1_io,N2DO(0))
        REAL(WP)  :: APIV(NCL1_io,N2DO(0))
        REAL(WP)  :: FI  (NCL1_io,N2DO(0))
        REAL(WP)  :: BCI (N2I:N2DO(MYID),2)
        REAL(WP)  :: COE1, COE2, AMU, APU
        INTEGER(4) :: IC, JC, KC, IM, IP, JSZ
       
        
        COE1 = FACOE_io * 4.0_wp/3.0_wp *DXQI    
        COE2 = FACOE_io *DXQI
        DO KC=1,NCL3
            AMIV = 0.0_WP
            ACIV = 0.0_WP
            APIV = 0.0_WP
            FI  = 0.0_WP
        
            DO IC=1, NCL1_io
                IM = IMV_io(IC)
                IP = IPV_io(IC)
                DO JC=N2I,N2DO(MYID)
                    IF(IDR==1) THEN ! U in X-IDRection
                        AMU = -COE1 * VISCOUSITY(IM,JC,KC)
                        APU = -COE1 * VISCOUSITY(IC,JC,KC)
                        AMIV(IC,JC) = AMU/RHO_STG(IM,JC,KC,IDR)
                        APIV(IC,JC) = APU/RHO_STG(IP,JC,KC,IDR)
                        ACIV(IC,JC) = 1.0_WP-(AMU+APU)/RHO_STG(IC,JC,KC,IDR)
                    END IF
                    
                    IF(IDR==2) THEN ! V in X-IDRection
                        AMU = -COE2 * MU_STG(IC,JC,KC,1)
                        APU = -COE2 * MU_STG(IP,JC,KC,1)
                        AMIV(IC,JC) = AMU/RHO_STG(IM,JC,KC,IDR)
                        APIV(IC,JC) = APU/RHO_STG(IP,JC,KC,IDR)
                        ACIV(IC,JC) = 1.0_WP-(AMU+APU)/RHO_STG(IC,JC,KC,IDR)
                    END IF
                    
                    IF(IDR==3) THEN ! W in X-IDRection
                        AMU = -COE2 * MU_STG(IC,JC,KC,2)
                        APU = -COE2 * MU_STG(IP,JC,KC,2)
                        AMIV(IC,JC) = AMU/RHO_STG(IM,JC,KC,IDR)
                        APIV(IC,JC) = APU/RHO_STG(IP,JC,KC,IDR)
                        ACIV(IC,JC) = 1.0_WP-(AMU+APU)/RHO_STG(IC,JC,KC,IDR)
                    END IF
                    FI(IC,JC)  = RHS_io(IC,JC,KC)       !RHSi 
                END DO
            END DO
            
            IF(TGFLOWflg) THEN ! to do
                DO JC=N2I,N2DO(MYID)
                    BCI(JC,1) = BC_U_SSTAR(JC,KC,IDR)
                    BCI(JC,2) = BC_TDMA(1, JC, KC, IDR)
                END DO
            END IF
            
            JSZ = N2DO(MYID)-N2I+1
            
            IF(TGFLOWFLG) THEN ! for inlet/outlet
                CALL TDMAIJI_nonCYC ( &
                        AMIV(1:NCL1_io,N2I:N2DO(MYID)),  &
                        ACIV(1:NCL1_io,N2I:N2DO(MYID)),  &
                        APIV(1:NCL1_io,N2I:N2DO(MYID)),  &
                        FI  (1:NCL1_io,N2I:N2DO(MYID)),  &
                        BCI (          N2I:N2DO(MYID) ,1:2),1,NCL1_io,N2I,JSZ )
            ELSE    ! for periodic x IDRection
                CALL TDMAIJI_CYC( &
                        AMIV(1:NCL1_io,N2I:N2DO(MYID)), &
                        ACIV(1:NCL1_io,N2I:N2DO(MYID)), &
                        APIV(1:NCL1_io,N2I:N2DO(MYID)), &
                        FI  (1:NCL1_io,N2I:N2DO(MYID)), &
                        1,NCL1_io,N2I,JSZ)
            END IF
            
            DO IC=1,NCL1_io
                DO JC=N2I,N2DO(MYID)
                    RHS_io(IC,JC,KC) = FI(IC,JC)
                END DO
            END DO

        ENDDO
      
        RETURN
    END SUBROUTINE
    
!******************************************************************************** 
    SUBROUTINE MOMFA2_Z_io(N2I,IDR)
        use mesh_info
        use flow_info
        use thermal_info
        use init_info
        IMPLICIT NONE
                  
        INTEGER(4),INTENT(IN)    :: N2I
        INTEGER(4),INTENT(IN)    :: IDR
      
        REAL(WP)      :: AMKV(NCL3,N2DO(MYID))
        REAL(WP)      :: ACKV(NCL3,N2DO(MYID))
        REAL(WP)      :: APKV(NCL3,N2DO(MYID))
        REAL(WP)      :: FK  (NCL3,N2DO(MYID))
        REAL(WP)      :: RMC2(N2DO(MYID))
        REAL(WP)      :: COE1, COE2, AMU, APU
        INTEGER(4)    :: IC, JC, KC, KM,KP,JJ, JSZ
       
        RMC2= 0.0_WP
        DO JC=N2I, N2DO(MYID)
            JJ = JCL2G(JC)
            IF (n2i.eq.2) THEN
                RMC2(JC)=RNDI2(JJ)
            ELSE
                RMC2(JC)=RCCI2(JJ)
            END IF
        END DO

        COE1 = FACOE_io * 4.0_wp/3.0_wp *DZQI    
        COE2 = FACOE_io *DZQI
        DO IC=1,NCL1_io
            AMKV = 0.0_WP
            ACKV = 0.0_WP
            APKV = 0.0_WP
            FK  = 0.0_WP
      
            DO KC=1,NCL3
                KM = KMV(KC)
                KP = KPV(KC)
                DO JC=N2I,N2DO(MYID)
                                    
                    IF(IDR==1) THEN ! U in Z-IDRection
                        AMU = -COE2 * MU_STG(IC,JC,KC,2)*RMC2(JC)
                        APU = -COE2 * MU_STG(IC,JC,KP,2)*RMC2(JC)
                        AMKV(KC,JC) = AMU/RHO_STG(IC,JC,KM,IDR)
                        APKV(KC,JC) = APU/RHO_STG(IC,JC,KP,IDR)
                        ACKV(KC,JC) = 1.0_WP-(AMU+APU)/RHO_STG(IC,JC,KC,IDR)
                    END IF
                    
                    IF(IDR==2) THEN ! V in Z-IDRection
                        AMU = -COE2 * MU_STG(IC,JC,KC,3)*RMC2(JC)
                        APU = -COE2 * MU_STG(IC,JC,KP,3)*RMC2(JC)
                        AMKV(KC,JC) = AMU/RHO_STG(IC,JC,KM,IDR)
                        APKV(KC,JC) = APU/RHO_STG(IC,JC,KP,IDR)
                        ACKV(KC,JC) = 1.0_WP-(AMU+APU)/RHO_STG(IC,JC,KC,IDR)
                    END IF
                    
                    IF(IDR==3) THEN ! W in Z-IDRection
                        AMU = -COE1 * VISCOUSITY(IC,JC,KM)*RMC2(JC)
                        APU = -COE1 * VISCOUSITY(IC,JC,KC)*RMC2(JC)
                        AMKV(KC,JC) = AMU/RHO_STG(IC,JC,KM,IDR)
                        APKV(KC,JC) = APU/RHO_STG(IC,JC,KP,IDR)
                        ACKV(KC,JC) = 1.0_WP-(AMU+APU)/RHO_STG(IC,JC,KC,IDR)
                    END IF
                    
                    FK(KC,JC) = RHS_io(IC,JC,KC)     
                ENDDO
            ENDDO

            JSZ = N2DO(MYID)-N2I+1
            CALL TDMAIJI_CYC( &
                    AMKV(1:NCL3,N2I:N2DO(MYID)), &
                    ACKV(1:NCL3,N2I:N2DO(MYID)), &
                    APKV(1:NCL3,N2I:N2DO(MYID)), &
                    FK  (1:NCL3,N2I:N2DO(MYID)), &
                    1,NCL3,N2I,JSZ)
                    
            DO KC=1,NCL3
                DO JC=N2I,N2DO(MYID)
                    RHS_io(IC,JC,KC)=FK(KC,JC)
                ENDDO
            ENDDO
         
        ENDDO
      
        RETURN
    END SUBROUTINE
    
!********************************************************************************
    SUBROUTINE MOMFA3_Y_io(IDR)
        use mesh_info
        use flow_info
        use thermal_info
        use init_info
        IMPLICIT NONE      
      
        INTEGER(4),INTENT(IN) :: IDR
      
        INTEGER(4)  :: IC, JC, KC, JM, JP , JJ
        REAL(WP)     :: FJJ (NCL1_io,NCL2) 
        REAL(WP)     :: AMJV(NCL1_io,NCL2)
        REAL(WP)     :: ACJV(NCL1_io,NCL2)
        REAL(WP)     :: APJV(NCL1_io,NCL2)
        REAL(WP)     :: BCJJ(NCL1_io,2)
        INTEGER(4)   :: NYS, NYE, JSZ
        REAL(WP)      :: F_io(NCL1_io,NCL2,N3DO(0) )
        REAL(WP)      :: FMU0(NCL1_io,NCL2,N3DO(0) )
        REAL(WP)      :: FDEN(NCL1_io,NCL2,N3DO(0) )
        REAL(WP)      :: VISN(NCL1_io,NCL3,2)
        REAL(WP)      :: COE1, COE2, AMU, APU, MUQT

        
        !IF(IDR.EQ.2) THEN
            !NYS=1
            !NYE=NND2
        !ELSE
            NYS=1
            NYE=NCL2
        !END IF 
        
        
        F_io    = 0.0_WP
        CALL TRASP23_Y2Z(NCL1_io, 1, N2DO(0), RHS_io, F_io)
       
        COE1 = FACOE_io * 4.0_wp/3.0_wp
        COE2 = FACOE_io
        
        IF(IDR==2) THEN
            CALL TRASP23_Y2Z(NCL1_io, 0, N2DO(0)+1, VISCOUSITY(:,:,:), FMU0)
            CALL TRASP23_Y2Z(NCL1_io, 0, N2DO(0)+1, RHO_STG(:,:,:,2),  FDEN)
        ELSE
            IF(MYID==NPSLV) THEN
                DO IC=1,NCL1_io
                    DO KC=1,NCL3
                        VISN(IC,KC,1) = VISCOUSITY(IC,N2DO(MYID),KC)
                        VISN(IC,KC,2) = VISCOUSITY(IC,N2DO(MYID)+1,KC)
                    END DO
                END DO
                CALL MPI_BCAST( VISN,       NCL1_io*NCL3*2, MPI_DOUBLE_PRECISION, NPSLV, ICOMM, IERROR )
            ELSE
                CALL MPI_BCAST( VISN,       NCL1_io*NCL3*2, MPI_DOUBLE_PRECISION, NPSLV, ICOMM, IERROR )
            END IF
            
            IF(IDR==1) THEN
                CALL TRASP23_Y2Z(NCL1_io, 0, N2DO(0)+1, MU_STG (:,:,:,IDR), FMU0)
                CALL TRASP23_Y2Z(NCL1_io, 0, N2DO(0)+1, RHO_STG(:,:,:,IDR), FDEN)
            END IF
        
            IF(IDR==3) THEN
                CALL TRASP23_Y2Z(NCL1_io, 0, N2DO(0)+1, MU_STG (:,:,:,IDR), FMU0)
                CALL TRASP23_Y2Z(NCL1_io, 0, N2DO(0)+1, RHO_STG(:,:,:,IDR), FDEN)
            END IF
        END IF
        
        
        DO KC=1,N3DO(MYID)
            FJJ  = 0.0_WP
            AMJV = 0.0_WP
            APJV = 0.0_WP
            ACJV = 0.0_WP
            BCJJ = 0.0_WP
!>         @note Prepare the coefficient for (1-A12)dU=dUStar         
            DO IC=1,NCL1_io 
                DO JC=NYS, NYE
                    JM = JC-1
                    JP = JC+1
                    IF(IDR==1 .OR. IDR==3) THEN ! U in Y-IDRection
                        ! ACVR IS NOT USED, BUT RE-CALCUATED....
                        IF(JC==1) THEN
                            MUQT= (VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JM,KC))*0.5_WP ! AT 1/4 ABOVE WALL...
                            AMU = -COE2 * AMVR1*MUQT
                            APU = -COE2 * APVR(JC,IDR)*FMU0(IC,JP,KC)
                            AMJV(IC,JC) = 0.0_wp
                            APJV(IC,JC) = APU/FDEN(IC,JP,KC)
                        ELSE IF (JC==NCL2) THEN
                            MUQT= (VISN(IC,KC,1) + VISN(IC,KC,2))*0.5_WP ! AT 1/4 ABOVE WALL...) 
                            AMU = -COE2 * AMVR(JC,IDR)*FMU0(IC,JC,KC)
                            APU = -COE2 * APVRN*MUQT
                            AMJV(IC,JC) = AMU/FDEN(IC,JM,KC)
                            APJV(IC,JC) = 0.0_wp
                        ELSE
                            AMU = -COE2 * AMVR(JC,IDR)*FMU0(IC,JC,KC)
                            APU = -COE2 * APVR(JC,IDR)*FMU0(IC,JP,KC)
                            AMJV(IC,JC) = AMU/FDEN(IC,JM,KC)
                            APJV(IC,JC) = APU/FDEN(IC,JP,KC)
                        END IF
                        ACJV(IC,JC) = 1.0_WP-(AMU+APU)/FDEN(IC,JC,KC)
                    END IF
                    
                    IF(IDR==2) THEN ! V in Y-IDRection
                        IF(JC==1) THEN
                            AMU = -COE1 * AMVR(JC,IDR) * VISCOUSITY(IC,JM,KC)
                            APU = -COE1 * APVR(JC,IDR) * VISCOUSITY(IC,JC,KC) 
                            AMJV(IC,JC) = AMU/RHO_STG(IC,JM,KC,IDR)
                            APJV(IC,JC) = APU/RHO_STG(IC,JP,KC,IDR)
                            ACJV(IC,JC) = 1.0_WP-(AMU+APU)/RHO_STG(IC,JC,KC,IDR)
                        ELSE
                            AMU = -COE1 * AMVR(JC,IDR) * FMU0(IC,JM,KC)
                            APU = -COE1 * APVR(JC,IDR) * FMU0(IC,JC,KC)
                            AMJV(IC,JC) = AMU/FDEN(IC,JM,KC)
                            APJV(IC,JC) = APU/FDEN(IC,JC,KC) 
                            ACJV(IC,JC) = 1.0_WP-(AMU+APU)/FDEN(IC,JC,KC)
                        END IF
                        
                    END IF
                    
                    FJJ(IC,JC) =  F_io(IC,JC,KC)
                ENDDO
                BCJJ(IC,:) = 0.0_WP
            ENDDO
         
            JSZ = NYE-NYS+1
            CALL TDMAIJJ_nonCYC (&
                    AMJV(1:NCL1_io,NYS:NYE), &
                    ACJV(1:NCL1_io,NYS:NYE), &
                    APJV(1:NCL1_io,NYS:NYE), &
                    FJJ (1:NCL1_io,NYS:NYE),  &
                    BCJJ(1:NCL1_io,1:2),   &
                    NYS,JSZ,1,NCL1_io)
         
            DO IC=1,NCL1_io
                DO JC=NYS, NYE
                    F_io(IC,JC,KC)=FJJ(IC,JC)
                ENDDO
            ENDDO
        ENDDO
      
        !CALL TRASP23G2L_RHS_io
        CALL TRASP23_Z2Y(NCL1_IO, 1, N2DO(0), RHS_IO, F_io)
      
        RETURN
    END SUBROUTINE MOMFA3_Y_io
    
    
    
 !***************************************************************************************************************   
    SUBROUTINE MOMFA_imcomp_io(NS,IDR)
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
 
        INTEGER(4),INTENT(IN)  :: NS
        INTEGER(4),INTENT(IN)  :: IDR
      
        INTEGER(4)  :: N2I
        INTEGER(4)  :: I, J, K
      
    
        FACOE_io=0.50_WP*TALP(NS)*DT*CVISC     !=alpha*dt/(2Re)
 
!>      @note N2I is setup for B.C. or not treatment.      
        N2I=1
        IF ((MYID.EQ.0).AND.(IDR.EQ.2).AND.(ICASE.NE.IBOX3P)) THEN
            N2I=2
        ENDIF

!>      @note Solving Eq.(A8a) to get dUstarstar
        CALL MOMFA1_imcomp_X_io(N2I,IDR)
!>      @note Solving Eq.(A8b) to get dUstar
        CALL MOMFA2_imcomp_Z_io(N2I)
!>      @note Solving Eq.(A8c) to get du
        CALL MOMFA3_imcomp_Y_io(IDR)
      
!>       @note update u using u_star-u=du       
        DO K=1,NCL3
            DO J=N2I,N2DO(MYID)
                DO I=1,NCL1_io
                    G_io(I,J,K,IDR)=RHS_io(I,J,K)+G_io(I,J,K,IDR)
                ENDDO
            ENDDO
        ENDDO

        RETURN
    END SUBROUTINE


    SUBROUTINE MOMFA1_imcomp_X_io(N2I,IDR)
        use mesh_info
        use flow_info
        IMPLICIT NONE
                  
        INTEGER(4),INTENT(IN)    :: N2I
        INTEGER(4),INTENT(IN)    :: IDR
      
        REAL(WP)  :: AMIV(NCL1_io,N2DO(MYID))
        REAL(WP)  :: ACIV(NCL1_io,N2DO(MYID))
        REAL(WP)  :: APIV(NCL1_io,N2DO(MYID))
        REAL(WP)  :: FI  (NCL1_io,N2DO(MYID))
        REAL(WP)  :: BCI(N2I:N2DO(MYID),2)
      
        INTEGER(4) :: I, J, K, JSZ
      
 
!>      @note To solve A8a, TDMA in x direction.         
        DO K=1,NCL3
      
            AMIV = 0.0_WP
            ACIV = 0.0_WP
            APIV = 0.0_WP
            FI   = 0.0_WP
            BCI  = 0.0_WP
         
            DO J=N2I,N2DO(MYID)
                DO I=1, NCL1_io
                    ACIV(I,J)= 1.0_WP-FACOE_io*(-2.0_WP*DXQI)             !bi
                    APIV(I,J)=       -FACOE_io*DXQI      !ci
                    AMIV(I,J)=       -FACOE_io*DXQI      !ai
                    FI(I,J)=       RHS_io(I,J,K)       !RHSi 
                END DO
                IF(TGFLOWflg) THEN
                    BCI(J,1) = BC_U_SSTAR(J,K,IDR)
                    BCI(J,2) = BC_TDMA(1, J, K, IDR)
                END IF
            END DO

            JSZ = N2DO(MYID)-N2I+1
            IF(TGFLOWFLG) THEN ! for inlet/outlet
                CALL TDMAIJI_nonCYC ( &
                        AMIV(1:NCL1_io,N2I:N2DO(MYID)),  &
                        ACIV(1:NCL1_io,N2I:N2DO(MYID)),  &
                        APIV(1:NCL1_io,N2I:N2DO(MYID)),  &
                        FI  (1:NCL1_io,N2I:N2DO(MYID)),  &
                        BCI (          N2I:N2DO(MYID) ,1:2),1,NCL1_io,N2I,JSZ )
            ELSE    ! for periodic x IDRection
                CALL TDMAIJI_CYC( &
                        AMIV(1:NCL1_io,N2I:N2DO(MYID)), &
                        ACIV(1:NCL1_io,N2I:N2DO(MYID)), &
                        APIV(1:NCL1_io,N2I:N2DO(MYID)), &
                        FI  (1:NCL1_io,N2I:N2DO(MYID)), &
                        1,NCL1_io,N2I,JSZ)
            END IF
                                       
            DO J=N2I,N2DO(MYID)
                DO I=1,NCL1_io
                    RHS_io(I,J,K) = FI(I,J)
                END DO
            END DO
         
        ENDDO
      
        RETURN
    END SUBROUTINE
    
!*******************************************************************     
    SUBROUTINE MOMFA2_imcomp_Z_io(N2I)
        use mesh_info
        use flow_info
        IMPLICIT NONE
                  
        INTEGER(4),INTENT(IN)    :: N2I
      
        REAL(WP)      :: AMKV(NCL3,N2DO(MYID))
        REAL(WP)      :: ACKV(NCL3,N2DO(MYID))
        REAL(WP)      :: APKV(NCL3,N2DO(MYID))
        REAL(WP)      :: FK  (NCL3,N2DO(MYID))
      
        REAL(WP)      :: RMC2(N2DO(MYID))
      
        INTEGER(4) :: I, J, K, JJ, JSZ
      

        AMKV = 0.0_WP
        ACKV = 0.0_WP
        APKV = 0.0_WP
        FK  = 0.0_WP
      
        RMC2= 0.0_WP

!>    @note To solve A8b, TDMA in z direction.    
        DO J=N2I, N2DO(MYID)
            JJ = JCL2G(J)
            IF (n2i.eq.2) THEN
                RMC2(J)=RNDI2(JJ)
            ELSE
                RMC2(J)=RCCI2(JJ)
            END IF
        END DO

        DO I=1,NCL1_io
      
            DO K=1,NCL3
                DO J=N2I,N2DO(MYID)
                    ACKV(K,J)= 1.0_WP-FACOE_io*(-2.0_WP*DZQI)*RMC2(J)
                    APKV(K,J)=       -FACOE_io*DZQI*RMC2(J)    !@     
                    AMKV(K,J)=       -FACOE_io*DZQI*RMC2(J)    !@
                    FK(K,J) = RHS_io(I,J,K)     
                ENDDO
            ENDDO

            JSZ = N2DO(MYID)-N2I+1
            CALL TDMAIJI_CYC(AMKV(1:NCL3,N2I:N2DO(MYID)),&
                             ACKV(1:NCL3,N2I:N2DO(MYID)),&
                             APKV(1:NCL3,N2I:N2DO(MYID)),&
                             FK(1:NCL3,N2I:N2DO(MYID)),&
                             1,NCL3,N2I,JSZ)
            DO K=1,NCL3
                DO J=N2I,N2DO(MYID)
                    RHS_io(I,J,K)=FK(K,J)
                ENDDO
            ENDDO
         
        ENDDO
      
        RETURN
    END SUBROUTINE
    
!**********************************************************************    
    SUBROUTINE MOMFA3_imcomp_Y_io(IDR)
        use mesh_info
        use flow_info
        IMPLICIT NONE      
      
        INTEGER(4),INTENT(IN) :: IDR
      
        INTEGER(4)  :: I, J, K 
        REAL(WP)     :: FJ  (NCL1_io,NCL2) 
        REAL(WP)     :: AMJV(NCL1_io,NCL2)
        REAL(WP)     :: ACJV(NCL1_io,NCL2)
        REAL(WP)     :: APJV(NCL1_io,NCL2)
        REAL(WP)     :: BCJ (NCL1_io,2)
        REAL(WP)     :: F_io(NCL1_io,NCL2,N3DO(0) )
        INTEGER(4)   :: NF 


        !IF(IDR.EQ.2) THEN
            !NF = NND2
        !ELSE
            NF = NCL2 
        !END IF 
        FJ  = 0.0_WP
        AMJV = 0.0_WP
        APJV = 0.0_WP
        ACJV = 0.0_WP
        BCJ = 0.0_WP
        F_IO = 0.0_WP
!>       @note SENDRECV globally QK=dUStar    
        CALL TRASP23_Y2Z(NCL1_io, 1, N2DO(0), RHS_io, F_io)
       
        DO K=1,N3DO(MYID)
!>        @note Prepare the coefficient for (1-A12)dU=dUStar         
            DO I=1,NCL1_io 
                DO J=1,NF
                    FJ(I,J) = F_io(I,J,K)
                    ACJV(I,J)= 1.0_WP-FACOE_io*ACVR(J,IDR)
                    APJV(I,J)=       -FACOE_io*APVR(J,IDR)
                    AMJV(I,J)=       -FACOE_io*AMVR(J,IDR)
                ENDDO
                BCJ(I,:) = 0.0_WP
            ENDDO
         
            CALL TDMAIJJ_nonCYC (AMJV(1:NCL1_io,1:NF),&
                              ACJV(1:NCL1_io,1:NF),&
                              APJV(1:NCL1_io,1:NF),&
                              FJ (1:NCL1_io,1:NF),&
                              BCJ(1:NCL1_io,1:2),&
                              1,NF,1,NCL1_io)
         
            DO I=1,NCL1_io
                DO J=1,NF
                    F_io(I,J,K)=FJ(I,J)
                ENDDO
            ENDDO
        ENDDO
      
      CALL TRASP23_Z2Y(NCL1_IO, 1, N2DO(0), RHS_IO, F_io)
      
      RETURN
      END SUBROUTINE
