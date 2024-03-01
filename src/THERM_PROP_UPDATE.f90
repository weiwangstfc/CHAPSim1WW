
    SUBROUTINE ENTHALPY_UPDATE(IREGION)
        USE THERMAL_INFO
        USE MESH_INFO
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: IREGION
        INTEGER(4) :: I,J,K, I1, I2
      
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
            I1 = 0
            I2 = NCL1_IO+1
        END IF
         
        DO J=1,N2DO(MYID)           
            DO I=I1,I2  
                DO K=1,NCL3
                    ENTHALPY(I,J,K) = RHOH(I,J,K)/DENSITY(I,J,K)
                END DO
            END DO
        END DO
      RETURN     
      END SUBROUTINE
      
      
!=================================================================================
    SUBROUTINE RHOH_UPDATE(IREGION)
        USE THERMAL_INFO
        USE MESH_INFO
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: IREGION
        INTEGER(4) :: I,J,K, I1, I2
      
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
            I1 = 0
            I2 = NCL1_IO+1
        END IF
         
        DO J=1,N2DO(MYID)           
            DO I=I1,I2  
                DO K=1,NCL3
                    RHOH(I,J,K) = ENTHALPY(I,J,K) * DENSITY(I,J,K)
                END DO
            END DO
        END DO
        RETURN     
    END SUBROUTINE


!*******************************************************************************
    SUBROUTINE DENSITY_UPDATE(IREGION)
        USE THERMAL_INFO
        USE FLOW_INFO
        USE MESH_INFO
        USE INIT_INFO
        IMPLICIT NONE
      
        INTEGER(4),INTENT(IN) :: IREGION
      
        INTEGER(4) :: I,J,K
        INTEGER(4) :: I1, I2
        REAL(WP)    :: H_TEMP, D_TEMP, T_TEMP, P_TEMP
      
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
            I1 = 0
            I2 = NCL1_IO+1
        END IF
        
      
        DO J=1,N2DO(MYID)            ! INCLUDE WALL B.C.
            DO I=I1,I2
                DO K=1,NCL3
                    H_TEMP = ENTHALPY(I,J,K)
                    IF(thermoStat==search_table) Then
                        CALL NIST_SLEVAL_HT(H_TEMP,T_TEMP)
                        CALL NIST_SLEVAL_TD(T_TEMP,D_TEMP)
                    else if (thermoStat==idealgas_law) then
                        P_TEMP = PR_IO(I,J,K)
                        !call idealgas_HT(H_TEMP,T_TEMP)
                        !call idealgas_TD(T_TEMP,D_TEMP, P_TEMP)
                    else
                    end if
                    DENSITY(I,J,K) = D_TEMP
                END DO
            END DO
        END DO
        
        

        RETURN
    END SUBROUTINE 
      
!*******************************************************************************
    SUBROUTINE THERM_PROP_UPDATE(IREGION)
!>    @NOTE: WHY TO USE "H" TO UPDATE OTHER PROPERTIES, NOT "H*D"
!>           Answer: FOR "HD", FOR EACH VALUE OF "HD", THERE IS MORE THAN 
!>                   ONE NUMBER FOUND FOR OTHER PROPERTIES.
!>                   AS "H" AND "D" HAS OPPOPSITE CHANGING TREND AS "T" CHANGES.
        use init_info
        USE thermal_info
        USE FLOW_INFO
        use mesh_info
        IMPLICIT NONE
      
        INTEGER(4),INTENT(IN) :: IREGION
        INTEGER(4) :: I,J,K
        INTEGER(4) :: I1, I2
        REAL(WP)    :: H_TEMP, T_TEMP, M_TEMP, K_TEMP, D_TEMP, P_TEMP

      
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
            I1 = 0
            I2 = NCL1_IO+1
        END IF
      
        DO J=1,N2DO(MYID)      
            DO I=I1,I2
                DO K=1,NCL3
                    H_TEMP = ENTHALPY(I,J,K)
               
                    IF(DABS(H_TEMP).LT.1.0E-10_wp) THEN
                        T_temp = 1.0_wp
                        M_temp = 1.0_wp
                        K_temp = 1.0_wp
                        D_temp = 1.0_wp
                        
                    ELSE
                        IF(thermoStat==search_table) Then
                            call NIST_SLEVAL_HT(H_temp,T_temp)
                            call NIST_SLEVAL_TM(T_temp,M_temp)
                            call NIST_SLEVAL_TK(T_temp,K_temp)
                            call NIST_SLEVAL_TD(T_temp,D_temp)
                        else if (thermoStat==idealgas_law) then
                            P_TEMP = PR_IO(I,J,K)
                            !call idealgas_HT(H_temp,T_temp)
                            !call idealgas_TM(T_temp,M_temp)
                            !call idealgas_TK(T_temp,K_temp)
                            !call idealgas_TD(T_temp,D_temp,P_TEMP)
                        else
                        end if
                        
                   
                        TEMPERATURE(I,J,K) = T_temp
                        VISCOUSITY(I,J,K) = M_temp
                        THERMCONDT(I,J,K) = K_temp
                        DENSITY(I,J,K)    = D_temp
                    END IF
                    
                    !RHOH(I,J,K) = DENSITY(I,J,K)*ENTHALPY(I,J,K)
                    
                    !if(myid==0) write(*,*) H_TEMP, T_temp, D_temp, M_temp, K_temp !test
                END DO
            END DO
        END DO
        
      
        RETURN
    END SUBROUTINE 
    
      
!**********************************************************************************
!    SUBROUTINE DENSITY_BKUP
!        USE thermal_info
!        USE mesh_info
!        IMPLICIT NONE
      
!        INTEGER(4) :: I,J,K, IS, IE
        
!        DENSITY1(:,:,:)   = DENSITY0(:,:,:)
!        DENSITY0(:,:,:)   = DENSITY(:,:,:)
        

!        RETURN
!    END SUBROUTINE
      
      
!********************************************************************************
    SUBROUTINE DENSITY_TIME_DERIVATION(NS)
        USE thermal_info
        USE MESH_INFO
        use flow_info
        use init_info
        IMPLICIT NONE
     
        INTEGER(4) :: IC,JC,KC, JP, JJ, IP, KP, NS
        real(wp)   :: DIVX, DIVY, DIVZ, COE2, COE3, DIVX0, DIVY0, DIVZ0
        
        DO JC=1,N2DO(MYID)   
            JP=JLPV(JC)
            JJ=JCL2G(JC)   
            COE2 =   DYFI(JJ) * RCCI1(JJ)
            COE3 =   DZI * RCCI2(JJ)
            DO IC=1,NCL1_IO  
                IP=IPV_io(IC) 
                DO KC=1,NCL3
                    KP=KPV(KC)
                    
!                    !METHOD 1==================================
!                    !==========d(\rho u)/dx at (i,j,k)===
!                    DIVX  = ( G_IO(IP,JC,KC,1) - G_IO(IC,JC,KC,1) ) * DXI
               
!                    !==========d(\rho v)/dy at (i,j,k)===
!                    DIVY  = ( G_io(IC,JP,KC,2) - G_io(IC,JC,KC,2) ) * COE2
               
!                    !==========d(\rho w)/dz at (i,j,k)===
!                    DIVZ  = ( G_io(IC,JC,KP,3) - G_io(IC,JC,KC,3) ) * COE3
!                    !(d-dp)/dt+gv=0  dp-d=gv*dt
!                    !DENSITYP(IC,JC,KC) = DENSITY(IC,JC,KC)+DT*(DIVX+DIVY+DIVZ)
!                    DrhoDtP(IC,JC,KC)= -1.0_wp*(DIVX+DIVY+DIVZ)

                    !METHOD 2====================================
                    !DrhoDtP(IC,JC,KC)= (DENSITY(IC,JC,KC)-DENSITY0(IC,JC,KC))*0.5_WP/DT
                    
                    !METHOD 4======Adam-Bashforth method===========
                    !DrhoDtP(IC,JC,KC)= (3.0_wp*DENSITY(IC,JC,KC)-4.0_WP*DENSITY0(IC,JC,KC)+DENSITY1(IC,JC,KC) )*0.5_WP/DT
                    
                    
                                        
                    !METHOD 3===========RK==(divegence free)=GOOD==================
                    !==========d(\rho u)/dx at (i,j,k)===
                    DIVX  = ( G_IO(IP,JC,KC,1) - G_IO(IC,JC,KC,1) )  * DXI
                    DIVX0 = ( G0_IO(IP,JC,KC,1)- G0_IO(IC,JC,KC,1) ) * DXI
               
                    !==========d(\rho v)/dy at (i,j,k)===
                    DIVY  = ( G_io(IC,JP,KC,2) - G_io(IC,JC,KC,2) )  * COE2
                    DIVY0 = ( G0_io(IC,JP,KC,2)- G0_io(IC,JC,KC,2) ) * COE2
               
                    !==========d(\rho w)/dz at (i,j,k)===
                    DIVZ  = ( G_io(IC,JC,KP,3) - G_io(IC,JC,KC,3) )  * COE3
                    DIVZ0 = ( G0_io(IC,JC,KP,3)- G0_io(IC,JC,KC,3) ) * COE3

                    DrhoDtP(IC,JC,KC)= ( (TGAM(NS)+TROH(NS))*(DIVX +DIVY +DIVZ) &
                                                  + TROH(NS)*(DIVX0+DIVY0+DIVZ0) )/(TGAM(NS)-2.0_wp)
                    
                END DO
            END DO
        END DO
        
        !THIS DENSITY IS USED ONLY IN CONTINUITY EQUATION.
     
     
        RETURN
     
    END SUBROUTINE
    
    
