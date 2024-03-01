    SUBROUTINE SOLVERRK3_MOM_IO(NS)
        use cparam
        use thermal_info 
        use mesh_info
        use flow_info
        use init_info
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: NS
        INTEGER(4) :: IDR
      
        !integer(4) :: JC, JJ, KC, IC
        
        !============for time stagereed in energy equation=====================
        G0_io(:,:,:,:) = G_io(:,:,:,:)
        
        
        !=======STEP1: CALCULATE OUTLET B.C.=================================
        IF(TGFLOWFLG) THEN
            CALL BC_TINLET_FLOW
            CALL BC_COUTLET_MOM_RK3(NS)  ! OUTLET
            
            IF (visthemflg == visimplicit) CALL BC_CBC_TDMA(NS) !!to do...
        END IF
    

        !=======STEP1: CALCULATE CONVECTION TERMS=============================== 
        CALL CONVECTION_X_io 
        CALL CONVECTION_Y_io 
        CALL CONVECTION_Z_io 
      
        !=======STEP2: CALCULATE THE WHOLE RHS FOR X,Y,Z==========================
        IF(visthemflg == visexplicit) THEN
            CALL DIVG_U_io
          
            IDR = 1
            CALL VISCOUS_ALL_EXPLT_X_io
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
          
            IDR = 2
            CALL VISCOUS_ALL_EXPLT_Y_io
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
          
            IDR = 3
            CALL VISCOUS_ALL_EXPLT_Z_io
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
    
            CALL MASSFLUX_CALC_IO
            
        ELSE IF (visthemflg == visimplicit) THEN
            IDR = 1
            CALL VISCOUS_PAR_EXPLT_X_io
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
            CALL MOMFA_io(NS,IDR)
            
            IDR = 2
            CALL VISCOUS_PAR_EXPLT_Y_io
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
            CALL MOMFA_io(NS,IDR)
            
            IDR = 3
            CALL VISCOUS_PAR_EXPLT_Z_io
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
            CALL MOMFA_io(NS,IDR)

        ELSE
            CALL ERRHDL(' No Such Scheme for viscous term',MYID) 
        END IF

        CALL INTFC_VARS3(1,NCL1_io,NCL1S,NCL1E,G_io)
        CALL BC_WALL_G_io
        
        !=======STEP4: CONSTRUCTING AND SOLVING POISSION EQ.=====================
        CALL DIVG_io(NS)
        
        !CALL DEBUG_WRT_LOCAL(RHSLLPHI_io,1,N2DO(MYID),'divg') !test
        
        IF(TGFLOWFLG) THEN
            CALL FISHPACK_POIS3D_SIMPLE
        ELSE
            CALL FFT99_POIS3D_periodicxz(IIO) !Method One 
            !CALL FISHPACK_POIS3D_SIMPLE ! Method Two,  good
        END IF
        
        CALL INTFC_VARS1(1,NCL1_io,NCL1S,NCL1E,DPH_io)
        CALL BC_WALL_DPH_io
        
        !CALL CHECK_FFT_SOLVER!test
        !CALL DEBUG_WRT_LOCAL(DPH_io,0,N2DO(MYID)+1,'dphi') !test
    
        !=======STEP5: CALCULATE PRESSURE=================
        CALL PRCALC_io(NS)
        
        !=======STEP6: CALCULATE mass flux=================
        CALL MASSFLUX_UPDATE_IO(NS)

        !=======STEP6: CALCULATE VELOCITY=================
        !CALL VELOCITY_CALC_io
        
        !==========check=============
        !CALL CHK_MassConsv_io
        !IF(MYID==0) WRITE(*,*) 'mass conservation(final): NS ', NS, CHK_MASS_CONSV0
        !IF(MYID==0)     WRITE(*,*) 'pressure1', &
        !PR_io(4,0,4),         PR_io(4,1,4),           PR_io(4,N2Do(myid),4),PR_io(4,N2DO(myid)+1,4)
        !IF(MYID==NPSLV) WRITE(*,*) 'pressure2', &
        !PR_io(4,N2DO(myid),4),PR_io(4,N2DO(myid)+1,4),PR_io(4,0,4),         PR_io(4,1,4)
        
        RETURN
    END SUBROUTINE
