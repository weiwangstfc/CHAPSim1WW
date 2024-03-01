!************************************************************************************************************************
    SUBROUTINE BC_PIPE_CENTRE_Q(IDOMAIN)
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: IDOMAIN
        INTEGER(4) :: I, K, KS
        
        IF(ICASE.NE.IPIPEC) RETURN
        IF(MYID.NE.0) RETURN
        
        IF(IDOMAIN==ITG) THEN
        
            DO I=1,NCL1_tg
                DO K=1,NCL3
                    KS=KSYM(K)                
                    Q_tg(I,0,K,1)=Q_tg(I,1,KS,1)       
                    Q_tg(I,0,K,3)=Q_tg(I,1,KS,3)  
                    Q_tg(I,0,K,2)=Q_tg(I,2,KS,2)               
                    Q_tg(I,1,K,2)=0.0_WP
                 END DO
            ENDDO 
            
        ELSE IF(IDOMAIN==IIO) THEN
        
            DO I=NCL1S,NCL1E
                DO K=1,NCL3
                    KS=KSYM(K)                
                    Q_io(I,0,K,1)=Q_io(I,1,KS,1)       
                    Q_io(I,0,K,3)=Q_io(I,1,KS,3)  
                    Q_io(I,0,K,2)=Q_io(I,2,KS,2)               
                    Q_io(I,1,K,2)=0.0_WP
                 END DO
            ENDDO 
            
        ELSE
        
        END IF
        
        
        RETURN
    END SUBROUTINE
    
!************************************************************************************************************************
    SUBROUTINE BC_PIPE_CENTRE_PR(IDOMAIN)
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: IDOMAIN
        INTEGER(4) :: I, K, KS
        
        IF(ICASE.NE.IPIPEC) RETURN
        IF(MYID.NE.0) RETURN
        
        IF(IDOMAIN==ITG) THEN
        
            DO I=1,NCL1_tg
                DO K=1,NCL3
                    KS=KSYM(K)                
                    PR_TG(I,0, K) = PR_TG(I,1,KS)
                 END DO
            ENDDO  
            
        ELSE IF(IDOMAIN==IIO) THEN
        
            DO I=NCL1S,NCL1E
                DO K=1,NCL3
                    KS=KSYM(K)                
                    PR_IO(I,0, K) = PR_IO(I,1,KS)
                 END DO
            ENDDO
            
        ELSE
        
        END IF

        RETURN
    END SUBROUTINE
    
!************************************************************************************************************************
    SUBROUTINE BC_PIPE_CENTRE_DPH(IDOMAIN)
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: IDOMAIN
        INTEGER(4) :: I, K, KS
        
        IF(ICASE.NE.IPIPEC) RETURN
        IF(MYID.NE.0) RETURN
        
        IF(IDOMAIN==ITG) THEN
        
            DO I=1,NCL1_tg
                DO K=1,NCL3
                    KS=KSYM(K)                
                    DPH_TG(I,0, K) = DPH_TG(I,1,KS)
                 END DO
            ENDDO  
            
        ELSE IF(IDOMAIN==IIO) THEN
        
            DO I=NCL1S,NCL1E
                DO K=1,NCL3
                    KS=KSYM(K)                
                    DPH_IO(I,0, K) = DPH_IO(I,1,KS)
                 END DO
            ENDDO
            
        ELSE
        
        END IF

        RETURN
    END SUBROUTINE 
    
!************************************************************************************************************************    
    SUBROUTINE BC_PIPE_CENTRE_G(IDOMAIN)
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: IDOMAIN
        INTEGER(4) :: I, K, KS
        
        IF(ICASE.NE.IPIPEC) RETURN
        IF(IDOMAIN.NE.IIO) RETURN
        IF(MYID.NE.0) RETURN
        
        DO I=NCL1S,NCL1E
            DO K=1,NCL3
                KS=KSYM(K)                
                G_io(I,0,K,1)=G_io(I,1,KS,1)       
                G_io(I,0,K,3)=G_io(I,1,KS,3)  
                G_io(I,0,K,2)=G_io(I,2,KS,2)               
                G_io(I,1,K,2)=0.0_WP
             END DO
        ENDDO 

        RETURN
    END SUBROUTINE 
    
!************************************************************************************************************************    
    SUBROUTINE BC_PIPE_CENTRE_THERMAL
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        INTEGER(4) :: I, K, KS
        
        IF(ICASE.NE.IPIPEC) RETURN
        IF(MYID.NE.0) RETURN
        
        DO I=NCL1S,NCL1E
            DO K=1,NCL3
                KS=KSYM(K)                
                ENTHALPY  (I,0,K) = ENTHALPY  (I,1,KS)
                DENSITY   (I,0,K) = DENSITY   (I,1,KS)
                TEMPERATURE(I,0,K) = TEMPERATURE(I,1,KS)
                VISCOUSITY(I,0,K) = VISCOUSITY(I,1,KS)
                THERMCONDT(I,0,K) = THERMCONDT(I,1,KS)
                RHOH      (I,0,K) = ENTHALPY  (I,1,KS) *  &
                                    DENSITY   (I,1,KS)
             END DO
        ENDDO 
        
        RETURN
    END SUBROUTINE 
    
!************************************************************************************************************************       
!************************************************************************************************************************       
!************************************************************************************************************************
    SUBROUTINE BC_WALL_Q(N,IDOMAIN)
    
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        
        INTEGER(4),INTENT(IN) :: N
        INTEGER(4),INTENT(IN) :: IDOMAIN
        
        INTEGER(4) :: I, K
        
        IF(N.EQ.Ibotwall .AND. ICASE.EQ.IPIPEC) RETURN
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        IF(IDOMAIN==ITG) THEN
        
            IF (MYID.EQ.0 .AND. N.EQ.Ibotwall) THEN !=========bottom wall==================
            
                DO I=1,NCL1_tg
                    DO K=1,NCL3                
                        Q_tg(I,0,K,1)=0.0_WP!-Q_tg(I,1,K,N)        
                        Q_tg(I,0,K,3)=0.0_WP!-Q_tg(I,3,K,N)                     
                        Q_tg(I,1,K,2)=0.0_WP
                        Q_tg(I,0,K,2)=0.0_WP
                    ENDDO
                ENDDO  
            
            END IF
            
            IF (MYID.EQ.NPSLV .AND. N.EQ.Itopwall) THEN !=========top wall==================
        
                DO I=1,NCL1_tg
                    DO K=1,NCL3
                        Q_tg(I,N2DO(MYID)+1,K,1) = 0.0_WP
                        Q_tg(I,N2DO(MYID)+1,K,3) = 0.0_WP
                        Q_tg(I,N2DO(MYID)+1,K,2) = 0.0_WP    
                    END DO
                END DO
        
            END IF
            
        ELSE IF(IDOMAIN==IIO) THEN
        
            IF (MYID.EQ.0 .AND. N.EQ.Ibotwall) THEN !=========bottom wall==================
            
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3                
                        Q_io(I,0,K,1)=0.0_WP!-Q_io(I,1,K,N)        
                        Q_io(I,0,K,3)=0.0_WP!-Q_io(I,3,K,N)                     
                        Q_io(I,1,K,2)=0.0_WP
                        Q_io(I,0,K,2)=0.0_WP
                    ENDDO
                ENDDO  
            
            END IF
            
            IF (MYID.EQ.NPSLV .AND. N.EQ.Itopwall) THEN !=========top wall==================
        
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3
                        Q_io(I,N2DO(MYID)+1,K,1) = 0.0_WP
                        Q_io(I,N2DO(MYID)+1,K,3) = 0.0_WP
                        Q_io(I,N2DO(MYID)+1,K,2) = 0.0_WP    
                    END DO
                END DO
        
            END IF
        
        ELSE
        
        END IF

        RETURN
    END SUBROUTINE
    
    
    
!************************************************************************************************************************
    SUBROUTINE BC_WALL_PR(N,IDOMAIN)
    
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        
        INTEGER(4),INTENT(IN) :: N
        INTEGER(4),INTENT(IN) :: IDOMAIN
        
        INTEGER(4) :: I, K
        
        IF(N.EQ.Ibotwall .AND. ICASE.EQ.IPIPEC) RETURN
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        IF(IDOMAIN==ITG) THEN
        
            IF (MYID.EQ.0 .AND. N.EQ.Ibotwall) THEN !=========bottom wall==================
            
                DO I=1,NCL1_tg
                    DO K=1,NCL3                
                       PR_tg(I,0,K) = PR_tg(I,1,K)
                    ENDDO
                ENDDO  
            
            END IF
            
            IF (MYID.EQ.NPSLV .AND. N.EQ.Itopwall) THEN !=========top wall==================
        
                DO I=1,NCL1_tg
                    DO K=1,NCL3
                        PR_tg(I,N2DO(MYID)+1,K) = PR_tg(I,N2DO(MYID),K)  
                    END DO
                END DO
        
            END IF
            
        ELSE IF(IDOMAIN==IIO) THEN
        
            IF (MYID.EQ.0 .AND. N.EQ.Ibotwall) THEN !=========bottom wall==================
            
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3                
                        PR_io(I,0,K) = PR_io(I,1,K)
                    ENDDO
                ENDDO  
            
            END IF
            
            IF (MYID.EQ.NPSLV .AND. N.EQ.Itopwall) THEN !=========top wall==================
        
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3
                        PR_io(I,N2DO(MYID)+1,K) = PR_io(I,N2DO(MYID),K)  
                    END DO
                END DO
        
            END IF
        
        ELSE
        
        END IF

        RETURN
    END SUBROUTINE
    
    
!************************************************************************************************************************
    SUBROUTINE BC_WALL_DPH(N,IDOMAIN)
    
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        
        INTEGER(4),INTENT(IN) :: N
        INTEGER(4),INTENT(IN) :: IDOMAIN
        
        INTEGER(4) :: I, K
        
        IF(N.EQ.Ibotwall .AND. ICASE.EQ.IPIPEC) RETURN
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        IF(IDOMAIN==ITG) THEN
        
            IF (MYID.EQ.0 .AND. N.EQ.Ibotwall) THEN !=========bottom wall==================
            
                DO I=1,NCL1_tg
                    DO K=1,NCL3                
                       DPH_tg(I,0,K) = DPH_tg(I,1,K)
                    ENDDO
                ENDDO  
            
            END IF
            
            IF (MYID.EQ.NPSLV .AND. N.EQ.Itopwall) THEN !=========top wall==================
        
                DO I=1,NCL1_tg
                    DO K=1,NCL3
                        DPH_tg(I,N2DO(MYID)+1,K) = DPH_tg(I,N2DO(MYID),K)  
                    END DO
                END DO
        
            END IF
            
        ELSE IF(IDOMAIN==IIO) THEN
        
            IF (MYID.EQ.0 .AND. N.EQ.Ibotwall) THEN !=========bottom wall==================
            
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3                
                        DPH_io(I,0,K) = DPH_io(I,1,K)
                    ENDDO
                ENDDO  
            
            END IF
            
            IF (MYID.EQ.NPSLV .AND. N.EQ.Itopwall) THEN !=========top wall==================
        
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3
                        DPH_io(I,N2DO(MYID)+1,K) = DPH_io(I,N2DO(MYID),K)  
                    END DO
                END DO
        
            END IF
        
        ELSE
        
        END IF

        RETURN
    END SUBROUTINE
    
  
!************************************************************************************************************************
    SUBROUTINE BC_WALL_G(N,IDOMAIN)
    
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        
        INTEGER(4),INTENT(IN) :: N
        INTEGER(4),INTENT(IN) :: IDOMAIN
        
        INTEGER(4) :: I, K
        
        IF(IDOMAIN.NE.IIO) RETURN
        IF(N.EQ.Ibotwall .AND. ICASE.EQ.IPIPEC) RETURN
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        IF (MYID.EQ.0 .AND. N.EQ.Ibotwall) THEN !=========bottom wall==================
            
            DO I=NCL1S,NCL1E
                DO K=1,NCL3                
                    G_io(I,0,K,1)=0.0_WP!-Q_io(I,1,K,N)        
                    G_io(I,0,K,3)=0.0_WP!-Q_io(I,3,K,N)                     
                    G_io(I,1,K,2)=0.0_WP
                    G_io(I,0,K,2)=0.0_WP
                ENDDO
            ENDDO  
        
        END IF
        
        IF (MYID.EQ.NPSLV .AND. N.EQ.Itopwall) THEN !=========top wall==================
    
            DO I=NCL1S,NCL1E
                DO K=1,NCL3
                    G_io(I,N2DO(MYID)+1,K,1) = 0.0_WP
                    G_io(I,N2DO(MYID)+1,K,3) = 0.0_WP
                    G_io(I,N2DO(MYID)+1,K,2) = 0.0_WP    
                END DO
            END DO
    
        END IF

        RETURN
    END SUBROUTINE
    
!************************************************************************************************************************   
    
    SUBROUTINE BC_WALL_ISOTHERMAL(N)
!>      @note
!>      All thermal variables at wall are located on the wall,
!>      Rather than the ghost cells in exterior
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE init_info
        IMPLICIT NONE
        
        INTEGER(4), INTENT(IN) :: N
        INTEGER(4) :: I,J, K
        
        IF(N.EQ.Ibotwall .AND. ICASE.EQ.IPIPEC) RETURN
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        !write(*,*) 'myid', myid, N !test
        IF(MYID.eq.0 .and. N==Ibotwall) THEN 
        !write(*,*) 'myidsss',  myid, N !test
            J=0
            DO I=NCL1S,NCL1E
                DO K=1,NCL3
                    TEMPERATURE(I,J,K) = T_WAL_GV(I,N)
                    ENTHALPY  (I,J,K) = H_WAL_GV(I,N)
                    DENSITY   (I,J,K) = D_WAL_GV(I,N)
                    VISCOUSITY(I,J,K) = M_WAL_GV(I,N)
                    THERMCONDT(I,J,K) = K_WAL_GV(I,N)
                    RHOH      (I,J,K) = ENTHALPY(I,J,K) * DENSITY(I,J,K)
                END DO
            END DO
            !WRITE(*,*) 'N,J,T',N, J, T_WAL_GV(1,N) , T_WAL_GV(NCL1_io/2,N)!test
        END IF
        
        IF(MYID.eq.NPSLV .and. N==Itopwall) THEN
        !write(*,*) 'myidsss',  myid, N !test 
            J=N2DO(MYID)+1
            DO I=NCL1S,NCL1E
                DO K=1,NCL3
                    TEMPERATURE(I,J,K) = T_WAL_GV(I,N)
                    ENTHALPY  (I,J,K) = H_WAL_GV(I,N)
                    DENSITY   (I,J,K) = D_WAL_GV(I,N)
                    VISCOUSITY(I,J,K) = M_WAL_GV(I,N)
                    THERMCONDT(I,J,K) = K_WAL_GV(I,N)
                    RHOH      (I,J,K) = ENTHALPY(I,J,K) * DENSITY(I,J,K)
                END DO
            END DO
            !WRITE(*,*) 'N,J,T',N, J, T_WAL_GV(1,N) , T_WAL_GV(NCL1_io/2,N)!test
        END IF
        
        
        
        RETURN
    END SUBROUTINE
   
    
    
    
!************************************************************************************************************************        
    SUBROUTINE BC_WALL_ISOFLUX(N,IREGION)
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        USE INIT_INFO
        IMPLICIT NONE
    
        INTEGER(4),INTENT(IN) :: IREGION, N
        INTEGER(4) :: J, K, I1, I2, I, JS, KS
        REAL(WP)    :: H_TEMP, M_TEMP, K_TEMP, D_TEMP, T_temp, P_temp
        
        IF(N.EQ.Ibotwall .AND. ICASE.EQ.IPIPEC) RETURN
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        IF(IREGION==IOULET) THEN  ! outlet
            I1 = NCL1_IO+1
            I2 = NCL1_IO+1
        ELSE IF(IREGION==IINLET) THEN !INLET
            I1 = 0
            I2 = 0
        ELSE IF(IREGION==IMAIND) THEN
            I1 = 1
            I2 = NCL1_IO
        ELSE IF(IREGION==IALLDM) THEN
            I1 = NCL1S
            I2 = NCL1E
        ELSE
            I1 = 1
            I2 = NCL1_IO
        END IF
        
        
        IF(MYID.eq.0 .or. MYID.eq.NPSLV) THEN
        
            IF(MYID.EQ.0 .AND. N==Ibotwall)     THEN
                J=0
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3
                        ENTHALPY(I,J,K) = 2.0_wp*ENTHALPY(I,J+1,K) &
                                        - YCL2ND_WFF(J+2)*ENTHALPY(I,J+2,K) &
                                        - YCL2ND_WFB(J+2)*ENTHALPY(I,J+1,K)
                        H_TEMP = ENTHALPY(I,J,K)
                        IF(thermoStat==search_table) Then
                            call NIST_SLEVAL_HT(H_temp,T_temp)
                            call NIST_SLEVAL_TM(T_temp,M_temp)
                            call NIST_SLEVAL_TK(T_temp,K_temp)
                            CALL NIST_SLEVAL_TD(T_temp,D_TEMP)
                        ELSE IF(thermoStat==idealgas_law) Then
                            P_temp = PR_IO(I,J,K)
                            !call idealgas_HT(H_temp,T_temp)
                            !call idealgas_TM(T_temp,M_temp)
                            !call idealgas_TK(T_temp,K_temp)
                            !call idealgas_TD(T_temp,D_TEMP, P_temp)
                        ELSE
                        END IF
                        
                        DENSITY(I,J,K)    = D_TEMP
                        TEMPERATURE(I,J,K) = T_temp
                        VISCOUSITY(I,J,K) = M_temp
                        THERMCONDT(I,J,K) = K_temp
                        
                        RHOH(I,J,K) = H_TEMP * D_TEMP
                    END DO
                END DO
                
            END IF
            
            IF(MYID.EQ.NPSLV .AND. N==Itopwall) THEN
                J = N2DO(MYID) + 1
                DO I=NCL1S,NCL1E
                    DO K=1,NCL3
                        ENTHALPY(I,J,K) = 2.0_wp*ENTHALPY(I,J-1,K) &
                                        - YCL2ND_WFF(J-1)*ENTHALPY(I,J-1,K) &
                                        - YCL2ND_WFB(J-1)*ENTHALPY(I,J-2,K)
                        H_TEMP = ENTHALPY(I,J,K)
                        IF(thermoStat==search_table) Then
                            call NIST_SLEVAL_HT(H_temp,T_temp)
                            call NIST_SLEVAL_TM(T_temp,M_temp)
                            call NIST_SLEVAL_TK(T_temp,K_temp)
                            CALL NIST_SLEVAL_TD(T_temp,D_TEMP)
                        ELSE IF(thermoStat==idealgas_law) Then
                            P_temp = PR_IO(I,J,K)
                            !call idealgas_HT(H_temp,T_temp)
                            !call idealgas_TM(T_temp,M_temp)
                            !call idealgas_TK(T_temp,K_temp)
                            !call idealgas_TD(T_temp,D_TEMP, P_temp)
                        ELSE
                        END IF
                        
                        DENSITY(I,J,K)    = D_TEMP
                        TEMPERATURE(I,J,K) = T_temp
                        VISCOUSITY(I,J,K) = M_temp
                        THERMCONDT(I,J,K) = K_temp
                        
                        RHOH(I,J,K) = H_TEMP * D_TEMP
                    END DO
                END DO
                
            END IF
        END IF
            

        RETURN
    END SUBROUTINE
    
    !***********************************************************************************************************************
    SUBROUTINE BC_WALL_THERMAL(IREGION)
        USE init_info
        USE thermal_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: IREGION
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        IF(BCWALLHEAT(itopwall)==isoFluxWall)    CALL BC_WALL_ISOFLUX   (itopwall, IREGION)
        IF(BCWALLHEAT(itopwall)==isoThermalWall) CALL BC_WALL_ISOTHERMAL(itopwall)
        
        IF(ICASE.ne.Ipipec) then
          if(BCWALLHEAT(ibotwall)==isoFluxWall)    CALL BC_WALL_ISOFLUX   (ibotwall, IREGION)
          IF(BCWALLHEAT(ibotwall)==isoThermalWall) CALL BC_WALL_ISOTHERMAL(ibotwall)
        ELSE
            CALL BC_PIPE_CENTRE_THERMAL
        END IF
    
        RETURN
    END SUBROUTINE
    !************************************************************************************************************************
    SUBROUTINE BC_WALL_Q_tg
        USE init_info
        IMPLICIT NONE
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        CALL BC_WALL_Q(Itopwall,ITG)
        IF(ICASE.ne.Ipipec) CALL BC_WALL_Q(Ibotwall,ITG)
        IF(ICASE.eq.Ipipec) CALL BC_PIPE_CENTRE_Q(ITG)
        RETURN
    END SUBROUTINE
    
    !************************************************************************************************************************
    SUBROUTINE BC_WALL_Q_io
        USE init_info
        IMPLICIT NONE
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        CALL BC_WALL_Q(Itopwall,IIO)
        IF(ICASE.ne.Ipipec) CALL BC_WALL_Q(Ibotwall,IIO)
        IF(ICASE.eq.Ipipec) CALL BC_PIPE_CENTRE_Q(IIO)
    
        RETURN
    END SUBROUTINE
    
    !************************************************************************************************************************
    SUBROUTINE BC_WALL_G_io
        USE init_info
        IMPLICIT NONE
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        CALL BC_WALL_G(Itopwall,IIO)
        IF(ICASE.ne.Ipipec) CALL BC_WALL_G(Ibotwall,IIO)
        IF(ICASE.eq.Ipipec) CALL BC_PIPE_CENTRE_G(IIO)
    
        RETURN
    END SUBROUTINE
    
    !************************************************************************************************************************
    SUBROUTINE BC_WALL_PR_TG
        USE init_info
        IMPLICIT NONE
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        CALL BC_WALL_PR(Itopwall,ITG)
        IF(ICASE.ne.Ipipec) CALL BC_WALL_PR(Ibotwall,ITG)
        IF(ICASE.eq.Ipipec) CALL BC_PIPE_CENTRE_PR(ITG)
    
        RETURN
    END SUBROUTINE
    
    !************************************************************************************************************************
    SUBROUTINE BC_WALL_PR_io
        USE init_info
        IMPLICIT NONE
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        CALL BC_WALL_PR(Itopwall,IIO)
        IF(ICASE.ne.Ipipec) CALL BC_WALL_PR(Ibotwall,IIO)
        IF(ICASE.eq.Ipipec) CALL BC_PIPE_CENTRE_PR(IIO)
    
        RETURN
    END SUBROUTINE
    
    !************************************************************************************************************************
    SUBROUTINE BC_WALL_DPH_TG
        USE init_info
        IMPLICIT NONE
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        CALL BC_WALL_DPH(Itopwall,ITG)
        IF(ICASE.ne.Ipipec) CALL BC_WALL_DPH(Ibotwall,ITG)
        IF(ICASE.eq.Ipipec) CALL BC_PIPE_CENTRE_DPH(ITG)
    
        RETURN
    END SUBROUTINE
    
    !************************************************************************************************************************
    SUBROUTINE BC_WALL_DPH_io
        USE init_info
        IMPLICIT NONE
        
        IF(BCY(1).NE.IBCWALL) RETURN
        IF(BCY(2).NE.IBCWALL) RETURN
        
        CALL BC_WALL_DPH(Itopwall,IIO)
        IF(ICASE.ne.Ipipec) CALL BC_WALL_DPH(Ibotwall,IIO)
        IF(ICASE.eq.Ipipec) CALL BC_PIPE_CENTRE_DPH(IIO)
    
        RETURN
    END SUBROUTINE

