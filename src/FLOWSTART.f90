
    SUBROUTINE FLOWSTART
        use init_info
        use flow_info 
        use thermal_info
        use mesh_info
        IMPLICIT NONE
     
        !========INITIAL FLOW FIELD and CONDITIONS===FOR TURBULENCE GENERATOR==========================================
        IF(TGFLOWflg) THEN
        
            IF (RSTflg_tg.EQ.0) THEN  
                IF (MYID.EQ.0) &
                CALL CHKHDL('14.TG: Flow initialization from random velocity field', myid)
                CALL RANDOM_FL_FLD_TG
                CALL CALC_INITIALIZATION_tg
                
            ELSE IF (RSTflg_tg.EQ.1) THEN 
                IF (MYID.EQ.0) &
                CALL CHKHDL('14.TG: Flow initialization from previous coarse mesh to a fine mesh using extrapolation.', myid)
                CALL INITIAL_INTERP_tg
                CALL WRT_INSTANT_VARS_TG
              
            ELSE IF (RSTflg_tg.EQ.2) THEN     
                IF (MYID.EQ.0)  &
                CALL CHKHDL('14.TG: Flow initialization from last step using the same mesh.', myid) 
                CALL RESTART_INSTANT_VARS_TG(RSTtim_tg)
            ENDIF
            
            CALL INTFC_VARS1(1,NCL1_tg,1,NCL1_tg,PR_tg)
            CALL INTFC_VARS3(1,NCL1_tg,1,NCL1_tg,Q_tg)
            CALL BC_WALL_Q_tg
            CALL BC_WALL_PR_TG
            CALL DIVGCK_tg
            IF(MYID.EQ.0) THEN
                CALL CHKHDL('    TG: MAX DIVERGENCE OF VELOCITY in restarted flow field',myid) 
                WRITE(*,'(A,20X,1ES15.7)') '#',MAXDIVGV_tg(2)
            END IF
            MAXDIVGV_tg(2) = MAXDIVGV_tg(1)
            CALL PP_TMEAN_INI_TG
            
        END IF
        
       
        !========INITIAL FLOW FIELD and CONDITIONS===FOR THE MAIN DOMAIN============================================
        IF(IOFLOWflg) THEN
            
            !================initial flow field===========================
            IF (RSTflg_io.EQ.0) THEN  ! 
                IF (MYID.EQ.0) CALL CHKHDL('14.IO: Flow initialization from random velocity field', myid)  
                
                CALL RANDOM_FL_THEML_FLD_io
                !CALL CALC_INITIALIZATION_io ! debug4chapsim2
            ELSE IF (RSTflg_io.EQ.1) THEN  !
                !====NREAD=1 : from coarse mesh to fine mesh using extrapolation======
                IF (MYID.EQ.0) &
                CALL CHKHDL('14.IO: Flow initialization from previous coarse mesh to a fine mesh using extrapolation.', myid)
                CALL INITIAL_INTERP_io
                CALL WRT_INSTANT_VARS_io
              
            ELSE IF (RSTflg_io.EQ.2) THEN     
                !====NREAD=2 restarting from the same mesh======
                IF (MYID.EQ.0) THEN
                    CALL CHKHDL('14.IO: Flow initialization from last step using the same mesh.', myid) 
                END IF
                
                CALL RESTART_INSTANT_VARS_IO(RSTtim_io)
                
                CALL DIVGCK_io    
                IF(MYID.EQ.0) THEN
                    CALL CHKHDL('   IO: Max divergence of initial flow field. Main Domain, Inlet and Outlet',myid) 
                    WRITE(*,'(A,25X,A,3ES23.15)') '#','==>', MAXDIVGV_io(1:3)
                END IF
                IF(weightedpressure==1) THEN
                    PR0_io(:,:,:,1)  = PR_io(:,:,:)
                    PR0_io(:,:,:,2)  = PR_io(:,:,:)
                END IF
            ENDIF
            
            !===========after initial flow field==============================
            CALL PP_TMEAN_INI_IO
            
            !==========below is only a repeat for security====================
            IF(thermlflg==1) THEN 
                CALL BC_WALL_THERMAL(IALLDM)
                IF(TGFLOWFLG) CALL BC_TINLET_THERML
                CALL THERM_PROP_UPDATE(IALLDM)
                IF(TGFLOWFLG) THEN
                    CALL INTFC_OUL_THERMAL_io
                    CALL INTFC_INL_THERMAL_io
                END IF
                CALL INTFC_MFD_THERMAL_io
            END IF
            
            CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
            CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
            CALL INTFC_VARS1(NCL1S,NCL1E,NCL1S,NCL1E,PR_io)
            CALL BC_WALL_Q_io
            CALL BC_WALL_G_io
            CALL BC_WALL_PR_io
            
        END IF
        
        CALL CALL_TEC360

        call wrt_3d_pt_debug(Q_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 1), 'ux', '@bf solv', 0, 0) ! debug4chapsim2
        call wrt_3d_pt_debug(Q_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 2), 'uy', '@bf solv', 0, 0) ! debug4chapsim2
        call wrt_3d_pt_debug(Q_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 3), 'uz', '@bf solv', 0, 0) ! debug4chapsim2
        call wrt_3d_pt_debug(PR_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3),    'pr', '@bf solv', 0, 0) ! debug4chapsim2   
        
        IF(PPROCESSONLY.eq.1) THEN
            CALL MPI_BARRIER(ICOMM,IERROR)
            IF(MYID==0) CALL CHKHDL('<===Only postprocessed results, now the code stops...==>', myid)  
            STOP
        END IF

        !call TEST_POISSON  ! debug4chapsim2

        
        
       RETURN     
     END SUBROUTINE
     
     
