!**********************************************************************
    SUBROUTINE BC_TINLET_FLOW
!>     @note
!>     Inlet bc will use the data from the adjacent cells in the turbulence 
!>     generator.
!>     For TG, energy equation may be solved with no heat flux fed in, therefore
!>     It is still like incompressible flow with constant thermal properties.
!>     Thus, the inlet thermal properties use constant values.      
    
        use mesh_info
        use flow_info
        use thermal_info
        IMPLICIT NONE
      
        INTEGER(4) :: J, K
        
        IF(.NOT. TGFLOWFLG) RETURN
      
        !===================FOR FLOW FIELD=========================
        DO J = 0, N2DO(MYID)+1
            DO K=1, NCL3
            
                !=====U=========
                Q_io  (0, J, K, 1) = Q_tg(NCL1_tg, J, K , 1) 
                Q_io  (1, J, K, 1) = Q_tg(1,       J, K , 1) 
            
                G_io  (0, J, K, 1) = Q_io  (0, J, K, 1) *   DENSITY(0,J,K)
                G_io  (1, J, K, 1) = Q_io  (1, J, K, 1) * ( DENSITY(0,J,K) + DENSITY(1,J,K) ) * XND2CL
                !G_io  (1, J, K, 1) = Q_io  (1, J, K, 1) *   DENSITY(0,J,K)  !CHECK!!!
            
                !=====V=========
                Q_io  (0, J, K, 2) = Q_tg(NCL1_tg, J, K , 2) 
                G_io  (0, J, K, 2) = Q_io  (0, J, K, 2) * DENSITY(0,J,K)
            
                !=====W=========
                Q_io  (0, J, K, 3) = Q_tg(NCL1_tg, J, K , 3) 
                G_io  (0, J, K, 3) = Q_io(0, J, K, 3) * DENSITY(0,J,K)
            
                !==========P===========
                PR_io (0, J, K)    = PR_tg(NCL1_tg, J, K)
            END DO
        END DO
        DPH_io(0, :, :)    = 0.0_WP
      
        !==as j=0~n2do+1, the intfc is not necessary=====
        !CALL INTFC_INL_G_io
        !CALL INTFC_INL_Q_io
        !CALL INTFC_INL_P_io

        RETURN
    END SUBROUTINE

!**********************************************************************
!**********************************************************************
    SUBROUTINE BC_TINLET_THERML
        use mesh_info
        use flow_info
        use thermal_info
        IMPLICIT NONE
      
        INTEGER(4) :: J, K
        
        IF(.NOT. TGFLOWFLG) RETURN
      
        DO J = 0, N2DO(MYID)+1, 1
            DO K=1, NCL3
                TEMPERATURE(0,J,K) = 1.0_WP
                THERMCONDT(0,J,K) = 1.0_WP
                DENSITY   (0,J,K) = 1.0_WP
                !DENSITY0  (0,J,K) = 1.0_WP
                VISCOUSITY(0,J,K) = 1.0_WP
                
                RHOH      (0,J,K) = 0.0_WP
                ENTHALPY  (0,J,K) = 0.0_WP
                
            END DO
        END DO
      
!        !================FOR THERMAL FIELD============================
!        !====ONLY TO BE CALLED BEFORE MAIN SOLVER ONCE=======
!        IF(BCWALLHEAT==isoFluxWall) THEN
!            DO J = 0, N2DO(MYID)+1, 1
!                DO K=1, NCL3
!                    TEMPRATURE(0,J,K) = 1.0_WP
!                    THERMCONDT(0,J,K) = 1.0_WP
!                    DENSITY   (0,J,K) = 1.0_WP
!                    !DENSITY0  (0,J,K) = 1.0_WP
!                    VISCOUSITY(0,J,K) = 1.0_WP
                    
!                    RHOH      (0,J,K) = 0.0_WP
!                    ENTHALPY  (0,J,K) = 0.0_WP
                    
!                END DO
!            END DO
!        END IF
        
!        IF(BCWALLHEAT==isoThermalWall) THEN
!            DO J = 0, N2DO(MYID)+1, 1
!                IF(MYID==0     .and. J==0             ) cycle
!                IF(MYID==NPSLV .and. J==(N2DO(MYID)+1)) cycle
!                DO K=1,NCL3
!                    TEMPRATURE(0,J,K) = 2.0_wp*TEMPRATURE(1,J,K)-TEMPRATURE(2,J,K)
!                    THERMCONDT(0,J,K) = 2.0_wp*THERMCONDT(1,J,K)-THERMCONDT(2,J,K)
!                    DENSITY   (0,J,K) = 2.0_wp*DENSITY   (1,J,K)-DENSITY   (2,J,K)
!                    !DENSITY0  (0,J,K) = 2.0_wp*DENSITY0  (1,J,K)-DENSITY0  (2,J,K)
!                    VISCOUSITY(0,J,K) = 2.0_wp*VISCOUSITY(1,J,K)-VISCOUSITY(2,J,K)
!                    ENTHALPY  (0,J,K) = 2.0_wp*ENTHALPY  (1,J,K)-ENTHALPY  (2,J,K)
!                    RHOH      (0,J,K) = 2.0_wp*RHOH      (1,J,K)-RHOH      (2,J,K)
!                END DO
!            END DO
!            CALL INTFC_INL_THERMAL_io
!        END IF
      
        RETURN
    END SUBROUTINE
