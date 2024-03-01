
    SUBROUTINE RANDOM_FL_FLD_TG
        use init_info
        use mesh_info
        use flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE      
             
        INTEGER(4) :: I, L
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        real(wp)    :: Umean1_tg
 
 
!>===========Generate scaled random Q and PR=1==========================
!       Turblence Generator***
        IF(MYID.EQ.0) CALL CHKHDL('   (1-1) TG: Generating random velocity',myid)
        CALL INIFIELD_FLOW_tg
        CALL VMAV_tg                  
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('         TG: VMV(1:3) in Generated random velocity field ',myid)      
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Max. U, V, W', (VMAX_tg(L),L=1,3)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Min. U, V, W', (VMIN_tg(L),L=1,3)
        END IF
       
!>================ Q=Q-Qxzmean=========================================
!        Turblence Generator***
        CALL INTFC_VARS3(1,NCL1_tg,1,NCL1_tg,Q_tg)
        CALL INITIAL_MEAN_TG
        DO L=1,3
            DO J=1,N2DO(MYID)
                Q_tg(:,J,:,L)= Q_tg(:,J,:,L) - UU(J,L,1)
            ENDDO
        END DO   
        CALL VMAV_tg
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (2-1) TG: VMV(1:3) in random velocity field after subtracting meanXZ ',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Max. U, V, W', (VMAX_tg(L),L=1,3)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Min. U, V, W', (VMIN_tg(L),L=1,3)
        END IF
       
!>===========update Q in the incoming flow direction====================
!        Turblence Generator*** 
        DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            DO I=1,NCL1_tg
                DO K=1,NCL3
                    Q_tg(I,J,K,NFLOW)= Q_tg(I,J,K,NFLOW) + Vini(JJ)  
                ENDDO
            ENDDO
        ENDDO         
            
        !IF(ICASE.NE.1) CALL VELO2RVELO_tg  
        
!>       @note Set Q(I,J,K,2)=v at Wall b.c. be zero.   
        DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            IF ( (JJ.EQ.1) .OR. (JJ.EQ.(NND2)) ) THEN
                Q_tg(:,J,:,2)= 0.0_WP
            ENDIF
        ENDDO   

!==============CORRECT THE MASS FLOW RATE=1=========================================
        CALL BULK_VELOCITY_TG(Umean1_tg)
        IF(MYID.EQ.0) call CHKRLHDL  ('        TG: The bulk velocity (original)=        ',MYID,Umean1_tg)
        
        DO J=1,N2DO(MYID)
            Q_tg(:,J,:,NFLOW)= Q_tg(:,J,:,NFLOW)/Umean1_tg
        END DO
        
        CALL BULK_VELOCITY_TG(Umean1_tg)
        IF(MYID.EQ.0) call CHKRLHDL  ('        TG: The bulk velocity (corrected)=       ',MYID,Umean1_tg)
       
!>============CHECK DIVERGENCE======================================
        CALL INTFC_VARS3(1,NCL1_tg,1,NCL1_tg,Q_tg)
        CALL DIVGCK_tg    
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (3-1) TG: Max divergence of initial flow field',myid) 
            WRITE(*,'(A,25X,A,1ES15.7)') '#','==>', MAXDIVGV_tg(2)
        END IF
        
        CALL CALL_TEC360
      
        RETURN
      
    END SUBROUTINE
      
!*********************************************************************************************************************
    SUBROUTINE CALC_INITIALIZATION_tg
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE      
       
        INTEGER(4) :: L
      
        !======Calculate the QP fields =========================================
!        Turblence Generator***
        IF(MYID.EQ.0) CALL CHKHDL('   (4-1) TG: Initial flow field, calc QP...',myid)
        
        CALL SOLVERRK3_MOM_tg(0)
       
        !==============CHECK INFO========================================
        CALL VMAV_tg
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (5-1) TG: Initial flow field with first calculated velocity, VMAX(1:3)',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Max. U, V, W', (VMAX_tg(L),L=1,3)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Min. U, V, W', (VMIN_tg(L),L=1,3)
        END IF
       
        CALL DIVGCK_tg
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (6-1) TG: MAX DIVERGENCE OF VELOCITY in final initial flow field',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==>', MAXDIVGV_tg(2)
        END IF
        MAXDIVGV_tg(2) = MAXDIVGV_tg(1)
        
        RETURN
    END SUBROUTINE
!**************************************************************************************************
    SUBROUTINE BULK_VELOCITY_TG(UUMEAN_work)
        use init_info
        use mesh_info
        use flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE      
             
        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        REAL(WP)    :: UUMEAN_work, UUMEAN
    
        UUMEAN = 0.0_wp
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            DO I=1,NCL1_tg
                DO K=1,NCL3
                    UUMEAN=UUMEAN + Q_TG(I,J,K,NFLOW)/DYFI(JJ)/RCCI1(JJ)
                END DO
            END DO
        ENDDO
        UUMEAN =  UUMEAN/DZI
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(UUMEAN, UUMEAN_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        
        UUMEAN_work = UUMEAN_work/(Area_INLET*DBLE(NCL1_tg))
        
        RETURN
    END SUBROUTINE

!**************************************************************************************************
    SUBROUTINE INIFIELD_FLOW_tg
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE      
       
       
        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        INTEGER(4) :: RDCOUNT, seed
!       INTEGER(4) :: IDUM
        REAL(WP)    :: VPERQ
        REAL(WP)    :: XRDM(3)
        REAL(WP)    :: XRD(3)
       
        INTEGER(4)  :: RANDOMTYPE=2  !=1 for real random, =2 for fixed random
       
   
       
        IF (RANDOMTYPE==1) THEN
            XRDM = 0.0_WP  
            !@note Generate the random flow field u,v,w.      
            DO J=1,N2DO(MYID)           ! only local y_cell no.      
                JJ=JCL2G(J)              
            !>@param VPERQ=scaled magnitude of perturbation based on Up (parabolic max. velo.) 
            !>@note for 1/4 near wall region, perturbation ratio is decreased by 25%. 
                VPERQ=VPERG                               
                IF((1.0_WP-DABS(YCC(JJ))).LT.0.250_WP) VPERQ=VPERG*SVPERG
              
                DO I=1,NCL1_tg
                    DO K=1,NCL3             
    
                        CALL RANDOM_NUMBER(XRDM)
                        !>@note Transfer random real no. 0~1 to -1~1                
                        XRD(1)=(-1.0_WP+2.0_WP*XRDM(1))
                        XRD(2)=(-1.0_WP+2.0_WP*XRDM(2))
                        XRD(3)=(-1.0_WP+2.0_WP*XRDM(3))
                    
                        !>@note the scaled random perturbations, in three directions.                
                        Q_tg(I,J,K,1)=VPERQ*XRD(1)
                        Q_tg(I,J,K,2)=VPERQ*XRD(2) /RNDI1(JJ)
                        Q_tg(I,J,K,3)=VPERQ*XRD(3) /RCCI1(JJ)
                                  
                    ENDDO
                ENDDO
                IF(JJ==1) Q_tg(:,J,:,2) = 0.0_WP
            ENDDO
        END IF
       
       
        IF (RANDOMTYPE==2) THEN
            XRDM = 0.0_WP  
            RDCOUNT = 0
            !>@note Generate the random flow field u,v,w.      
            DO J=1,N2DO(MYID)           ! only local y_cell no.      
                JJ=JCL2G(J)              
                !>@param VPERQ=scaled magnitude of perturbation based on Up (parabolic max. velo.) 
                !>@note for 1/4 near wall region, perturbation ratio is decreased by 25%. 
                VPERQ=VPERG                               
                IF((1.0_WP-DABS(YCC(JJ))).LT.0.250_WP) VPERQ=VPERG*SVPERG
              
                DO I=1,NCL1_tg
                    DO K=1,NCL3
                        !RDCOUNT = RDCOUNT + 1
                        RDCOUNT = I + K + JJ  ! make it independent of processor number.
                        seed = 1973+RDCOUNT          
                        call random_initialize ( seed )
                        call rvec_random ( -1.0_WP, 1.0_WP, 3, XRD )
                        !>@note the scaled random perturbations, in three directions.                
                        Q_tg(I,J,K,1)=VPERQ*XRD(1)
                        Q_tg(I,J,K,2)=VPERQ*XRD(2) /RNDI1(JJ)
                        Q_tg(I,J,K,3)=VPERQ*XRD(3) /RCCI1(JJ)
                        !write(*,'(3I4.2,3ES13.5)') I,J,K,XRD(1),XRD(2),XRD(3) !test             
                                       
                    ENDDO
                ENDDO
                IF(JJ==1) Q_tg(:,J,:,2) = 0.0_WP
            ENDDO
        END IF
        
        PR_tg(:,:,:)=0.0_WP 
        
 
        RETURN
    END SUBROUTINE
    
!***********************************************************************************    
!    SUBROUTINE VELO2RVELO_TG
!        use init_info
!        use mesh_info
!        use flow_info
!        IMPLICIT NONE    
        
!        INTEGER(4) :: I, J, K, NYI, JJ
        
!        NYI = 1
!        IF(MYID==0) THEN
!            J = 1
!            DO I=1, NCL1_TG
!                DO K=1,NCL3
!                    Q_TG(I,J,K,1) = Q_TG(I,J,K,1) 
!                    Q_TG(I,J,K,2) = 0.0_WP
!                    Q_TG(I,J,K,3) = Q_TG(I,J,K,3) / RCCI1(JJ)
!                END DO
!            END DO
!        NYI = 2
!        END IF
        
!        DO J=NYI,N2DO(MYID)
!            JJ = JCL2G(J)
!            DO I=1, NCL1_TG
!                DO K=1,NCL3
!                    Q_TG(I,J,K,1) = Q_TG(I,J,K,1) 
!                    Q_TG(I,J,K,2) = Q_TG(I,J,K,2) / RNDI1(JJ)
!                    Q_TG(I,J,K,3) = Q_TG(I,J,K,3) / RCCI1(JJ)
!                END DO
!            END DO
!        END DO
        
!        RETURN
!    END SUBROUTINE
