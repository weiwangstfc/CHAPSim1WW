 !**********************************************************************
    SUBROUTINE BC_COUTLET_ENEG_RK3(NS)
        use mesh_info
        use flow_info
        use init_info
        use thermal_info
        IMPLICIT NONE
      
        INTEGER(4),INTENT(IN) :: NS
      
        INTEGER(4) :: J, K
        REAL(WP)    :: BC_RHS
        REAL(WP)    :: BC_CONV
        
        
        IF(.NOT. TGFLOWFLG) RETURN

        CALL CONVCTION_OUTLET_U
      
        DO K=1, NCL3
            DO J=1, N2DO(myid)
                !U_OUTLET=      !Q_io(NCL1_io+1,J,K,1) which may be negative.
                BC_CONV=-DXI*U_OUTLET * &
                        ( RHOH(NCL1_io+1,J,K) - RHOH(NCL1_io,  J,K) ) 
                BC_RHS = TGAM(NS)*BC_CONV + TROH(NS)*BC_CONV0_ENG(J,K)
                BC_CONV0_ENG(J,K)=BC_CONV

                RHOH(NCL1_io+1,J,K) = RHOH(NCL1_io+1,J,K) + DT * BC_RHS
            END DO
        END DO
        
        CALL ENTHALPY_UPDATE(IOULET)
        CALL DENSITY_UPDATE(IOULET)
        CALL ENTHALPY_UPDATE(IOULET)
      
        CALL THERM_PROP_UPDATE(IOULET)
        CALL INTFC_OUL_THERMAL_io
        CALL BC_WALL_THERMAL(IOULET)
        
        RETURN
    END SUBROUTINE 
