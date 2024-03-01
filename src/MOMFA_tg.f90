!***********************************************************************
!> @author 
!> Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Factorization Approximation of the momentume equation.
!> @param AL  [in] alpha coefficients in Eq.A4a of Mehdi paper.
!> @param IDR [in] velocity direction. 
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! 06/02/2014- Initial Version, by Wei Wang
!**********************************************************************    
    SUBROUTINE MOMFA_tg(NS,IDR)
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
 
        INTEGER(4),INTENT(IN)  :: NS
        INTEGER(4),INTENT(IN)  :: IDR
      
        INTEGER(4)  :: N2I
        INTEGER(4)  :: I, J, K
      
      
        FACOE_TG=0.50_WP*TALP(NS)*DT*CVISC     !=alpha*dt/(2Re)
   
        N2I=1
        IF ((MYID.EQ.0).AND.(IDR.EQ.2))  N2I=2

        CALL MOMFA1_X_tg(N2I,IDR)
        CALL MOMFA2_Z_tg(N2I)
        CALL MOMFA3_Y_tg(IDR)
        
        DO K=1,NCL3
            DO J=N2I,N2DO(MYID)
                DO I=1,NCL1_tg
                    Q_tg(I,J,K,IDR)=RHS_tg(I,J,K)+Q_tg(I,J,K,IDR)
                ENDDO
            ENDDO
        ENDDO

        RETURN
    END SUBROUTINE


!*************************************************************************************
    SUBROUTINE MOMFA1_X_tg(N2I,IDR)
        use mesh_info
        use flow_info
        IMPLICIT NONE
                  
        INTEGER(4),INTENT(IN)    :: N2I
        INTEGER(4),INTENT(IN)    :: IDR
      
        REAL(WP)  :: AMIV(NCL1_tg,N2DO(0))
        REAL(WP)  :: ACIV(NCL1_tg,N2DO(0))
        REAL(WP)  :: APIV(NCL1_tg,N2DO(0))
        REAL(WP)  :: FI  (NCL1_tg,N2DO(0))
        REAL(WP)  :: COE
      
        INTEGER(4) :: I, J, K, JSZ
       
        AMIV = 0.0_WP
        ACIV = 0.0_WP
        APIV = 0.0_WP
        FI  = 0.0_WP
            
        COE = -FACOE_TG*DXQI
        DO K=1,NCL3
        
        
            DO J=N2I,N2DO(MYID)
                DO I=1, NCL1_tg
                    APIV(I,J)= COE     !ci
                    AMIV(I,J)= COE      !ai
                    ACIV(I,J)= 1.0_WP-APIV(I,J)-AMIV(I,J)             !bi
                    FI(I,J)  = RHS_tg(I,J,K)       !RHSi 
                END DO
            END DO
!             CALL TRIPVI(AMI,ACI,API,FI,1,N1M,N2I,N2DO)
            JSZ = N2DO(MYID)-N2I+1
            CALL TDMAIJI_CYC(AMIV(1:NCL1_tg,N2I:N2DO(MYID)),&
                             ACIV(1:NCL1_tg,N2I:N2DO(MYID)),&
                             APIV(1:NCL1_tg,N2I:N2DO(MYID)),&
                             FI(1:NCL1_tg,N2I:N2DO(MYID)),&
                            1,NCL1_tg,N2I,JSZ)
            DO J=N2I,N2DO(MYID)
                DO I=1,NCL1_tg
                    RHS_tg(I,J,K) = FI(I,J)
                END DO
            END DO
         
            IF(IOFLOWflg) THEN
                DO J=N2I,N2DO(MYID)
                    BC_U_SSTAR(J,K,IDR) = RHS_tg(NCL1_tg,J,K)
                END DO
            END IF
         
        ENDDO
      
        RETURN
    END SUBROUTINE
    
    
!*************************************************************************************
    SUBROUTINE MOMFA2_Z_tg(N2I)
        use mesh_info
        use flow_info
        IMPLICIT NONE
                  
        INTEGER(4),INTENT(IN)    :: N2I
      
        REAL(WP)      :: AMKV(NCL3,N2DO(0))
        REAL(WP)      :: ACKV(NCL3,N2DO(0))
        REAL(WP)      :: APKV(NCL3,N2DO(0))
        REAL(WP)      :: FK  (NCL3,N2DO(0))
      
        REAL(WP)      :: RMC2(N2DO(MYID))
      
        INTEGER(4) :: I, J, K, JJ, JSZ
      
        AMKV = 0.0_WP
        ACKV = 0.0_WP
        APKV = 0.0_WP
        FK  = 0.0_WP
      
        RMC2= 0.0_WP
     
        DO J=N2I, N2DO(MYID)
            JJ = JCL2G(J)
            IF (n2i.eq.2) THEN
                RMC2(J)=RNDI2(jj) * (-FACOE_TG*DZQI)
            ELSE
                RMC2(J)=RCCI2(jj) * (-FACOE_TG*DZQI)
            END IF
        END DO

        DO I=1,NCL1_tg
      
            DO K=1,NCL3
                DO J=N2I,N2DO(MYID)
                    APKV(K,J)=RMC2(J)    !@     
                    AMKV(K,J)=RMC2(J)    !@
                    ACKV(K,J)= 1.0_WP-APKV(K,J)-AMKV(K,J)
                    FK(K,J)  = RHS_tg(I,J,K)     
                ENDDO
            ENDDO

!         CALL  TRVPJK(AMK,ACK,APK,FK,1,N3M,N2I,N2DO )
            JSZ = N2DO(MYID)-N2I+1
            CALL TDMAIJI_CYC(AMKV(1:NCL3,N2I:N2DO(MYID)),&
                             ACKV(1:NCL3,N2I:N2DO(MYID)),&
                             APKV(1:NCL3,N2I:N2DO(MYID)),&
                             FK(1:NCL3,N2I:N2DO(MYID)),&
                             1,NCL3,N2I,JSZ)
            DO K=1,NCL3
                DO J=N2I,N2DO(MYID)
                    RHS_tg(I,J,K)=FK(K,J)
                ENDDO
            ENDDO
         
        ENDDO
      
        RETURN
    END SUBROUTINE
    
!*************************************************************************************
    SUBROUTINE MOMFA3_Y_tg(IDR)
        use mesh_info
        use flow_info
        IMPLICIT NONE      
      
        INTEGER(4),INTENT(IN) :: IDR
      
        INTEGER(4)  :: I, J, K 
        REAL(WP)     :: FJJ (NCL1_tg,NND2) 
        REAL(WP)     :: AMJV(NCL1_tg,NND2)
        REAL(WP)     :: ACJV(NCL1_tg,NND2)
        REAL(WP)     :: APJV(NCL1_tg,NND2)
        REAL(WP)     :: BCJJ(NCL1_tg,2)
        
        REAL(WP)     :: F(NCL1_TG,NCL2,N3DO(0) )
        
        INTEGER(4)  :: NYS, NYE, JSZ


        !IF(IDR.EQ.2) THEN
            !NYS=1
            !NYE=NND2
        !ELSE
            NYS=1
            NYE=NCL2
        !END IF 
        FJJ = 0.0_WP
        AMJV = 0.0_WP
        APJV = 0.0_WP
        ACJV = 0.0_WP
        BCJJ = 0.0_WP
        F    = 0.0_WP
  
        !CALL TRASP23L2G_RHS
        CALL TRASP23_Y2Z(NCL1_tg, 1, N2DO(0), RHS_TG, F)
       
        DO K=1,N3DO(MYID)
     
            DO I=1,NCL1_tg 
                DO J=NYS, NYE
                    FJJ(I,J) =  F(I,J,K)
                    ACJV(I,J)= 1.0_WP-FACOE_TG*ACVR(J,IDR)
                    APJV(I,J)=-FACOE_TG*APVR(J,IDR)
                    AMJV(I,J)=-FACOE_TG*AMVR(J,IDR)
                ENDDO
                BCJJ(I,:) = 0.0_WP
            ENDDO
         
            JSZ = NYE-NYS+1
            CALL TDMAIJJ_nonCYC (AMJV(1:NCL1_tg,NYS:NYE),&
                                 ACJV(1:NCL1_tg,NYS:NYE),&
                                 APJV(1:NCL1_tg,NYS:NYE),&
                                 FJJ(1:NCL1_tg,NYS:NYE),&
                                 BCJJ(1:NCL1_tg,1:2),&
                                 NYS,JSZ,1,NCL1_tg)
         
            DO I=1,NCL1_tg
                DO J=NYS,NYE
                    F(I,J,K)=FJJ(I,J)
                ENDDO
            ENDDO
        ENDDO
      
        !CALL TRASP23G2L_RHS
        CALL TRASP23_Z2Y(NCL1_TG, 1, N2DO(0), RHS_TG, F)
      
        RETURN
    END SUBROUTINE
    
