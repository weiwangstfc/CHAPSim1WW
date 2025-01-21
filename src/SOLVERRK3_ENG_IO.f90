    SUBROUTINE SOLVERRK3_ENG_IO(NS)
        USE init_info
        use mesh_info
        USE thermal_info
        USE flow_info
        IMPLICIT NONE
      
        INTEGER(4),INTENT(IN) :: NS
        INTEGER(4) :: I,J, K
        REAL(WP)     :: COE1, COE0
        
        
        VISCOUSITY0(:,:,:) = VISCOUSITY(:,:,:) ! for time stagger in momentum equation.
        
        DENSITY1(:,:,:)   = DENSITY0(:,:,:)
        DENSITY0(:,:,:)   = DENSITY(:,:,:)

        !========STEP 2: MAIN DOMAIN  ENERGY EQUATION ==\RHOH=================
        COE1 = TGAM(NS)*DT
        COE0 = TROH(NS)*DT

        CALL RHS_ENERGY_EXPLICIT
#ifdef DEBUG
        call wrt_3d_pt_debug(RHS_ENERGY(1:NCL1_io, 1:N2DO(myid), 1:NCL3), '', 'energy_rhs@bf st', ITERG, NS) ! debug4chapsim2
#endif
        RHOH0(:,:,:)      = RHOH(:,:,:)
        DO J=1,N2DO(MYID)
            DO I=1,NCL1_IO
                DO K=1,NCL3
                    RHOH(I,J,K) = RHOH(I,J,K) + COE1*RHS_ENERGY (I,J,K) + COE0*RHS_ENERGY0(I,J,K)
                    IF(thermoStat==search_table) RHOH(I,J,K) = DMIN1(RHOHmax+REALMIN, RHOH(I,J,K)) ! limiter!
                    
                    RHS_ENERGY0(I,J,K) = RHS_ENERGY(I,J,K)
                END DO
            END DO 
        END DO
#ifdef DEBUG
        call wrt_3d_pt_debug(RHOH(1:NCL1_io, 1:N2DO(myid), 1:NCL3), '', 'rhoh@af st', ITERG, NS) ! debug4chapsim2
        call wrt_3d_pt_debug(TEMPERATURE(1:NCL1_io, 1:N2DO(myid), 1:NCL3), '', 'T@bf st', ITERG, NS) ! debug4chapsim2
#endif
        !==================check=============
        !IF(NS==3) CALL CHK_EnegConsv_io
        !IF(MYID==0) WRITE(*,*) 'eneg conservation(inter): NS ', NS, CHK_ENEG_CONSV0
      
        !========STEP 3:  PROVISONAL ESTIMATE FOR H================
        CALL ENTHALPY_UPDATE(IMAIND) ! rhoh/rho
      
        !========STEP 4:  DENSITY UPDATE===========================
        CALL DENSITY_UPDATE(IMAIND)  ! h->T->rho
      
        !========STEP 5:  ENTHALPY UPDATE==========================
        CALL ENTHALPY_UPDATE(IMAIND) ! rhoh/rho
      
        !========STEP 6:  THERMAL_PROPERTY_UPDATE==================
        CALL THERM_PROP_UPDATE(IMAIND) ! no wall b.c.
      
        !========STEP 7: SET UP THERMAL PROPERTIES IN INTERFACES============
        CALL INTFC_MFD_THERMAL_io
        
        !========STEP 8:  b.c. ENERGY EQUATION ==\RHOH=================
        IF(TGFLOWFLG) THEN
            CALL BC_TINLET_THERML
            CALL BC_COUTLET_ENEG_RK3(NS)
        END IF
        
        CALL BC_WALL_THERMAL(IMAIND)
        
        !CALL DEBUG_WRT_THERMAL
        !==================check=============
        !CALL CHK_EnegConsv_io
        !IF(MYID==0) WRITE(*,*) 'eneg conservation(final): NS ', NS, CHK_ENEG_CONSV0
#ifdef DEBUG
        call wrt_3d_pt_debug(TEMPERATURE(1:NCL1_io, 1:N2DO(myid), 1:NCL3), '', 'T@af st', ITERG, NS) ! debug4chapsim2
        write(*,*) 'T-e', TEMPERATURE(1, 1:4, 1)
#endif
        RETURN
    END SUBROUTINE 
