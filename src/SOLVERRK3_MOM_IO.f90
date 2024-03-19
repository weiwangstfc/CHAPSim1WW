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
            CALL VISCOUS_ALL_EXPLT_X_io !test
            call wrt_3d_pt_debug(Qtmp_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'ConVisX', '@bf stepping', ITERG, NS) ! debug4chapsim2
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
            call wrt_3d_pt_debug(Qtmp_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'ConVisX', '@total', ITERG, NS) ! debug4chapsim2
          
            IDR = 2
            CALL VISCOUS_ALL_EXPLT_Y_io
            call wrt_3d_pt_debug(DPH_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'ConVisY', '@bf stepping', ITERG, NS) ! debug4chapsim2
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
            call wrt_3d_pt_debug(DPH_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'ConVisY', '@total', ITERG, NS) ! debug4chapsim2
          
            IDR = 3
            CALL VISCOUS_ALL_EXPLT_Z_io
            call wrt_3d_pt_debug(RHSLLPHI_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'ConVisZ', '@bf stepping', ITERG, NS) ! debug4chapsim2
            CALL RHS_MOM_EXPLICIT_io(NS,IDR)
            call wrt_3d_pt_debug(RHSLLPHI_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'ConVisZ', '@total', ITERG, NS) ! debug4chapsim2
    
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

        call wrt_3d_pt_debug(G_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 1), 'ux', '@bf divg', ITERG, NS) ! debug4chapsim2
        call wrt_3d_pt_debug(G_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 2), 'uy', '@bf divg', ITERG, NS) ! debug4chapsim2
        call wrt_3d_pt_debug(G_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 3), 'uz', '@bf divg', ITERG, NS) ! debug4chapsim2
        !=======STEP4: CONSTRUCTING AND SOLVING POISSION EQ.=====================
        CALL DIVG_io(NS)
        call wrt_3d_pt_debug (RHSLLPHI_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'PhiRHS', '@RHS phi', ITERG, NS) ! debug4chapsim2
        call wrt_3d_all_debug(RHSLLPHI_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'PhiRHS', '@RHS phi', ITERG, NS)
        !CALL DEBUG_WRT_LOCAL(RHSLLPHI_io,1,N2DO(MYID),'divg') !test
        
        IF(TGFLOWFLG) THEN
            CALL FISHPACK_POIS3D_SIMPLE
        ELSE
            CALL FFT99_POIS3D_periodicxz(IIO) !Method One 
            !CALL FISHPACK_POIS3D_SIMPLE ! Method Two,  good
        END IF
        
        CALL INTFC_VARS1(1,NCL1_io,NCL1S,NCL1E,DPH_io)
        CALL BC_WALL_DPH_io

        !if(myid == 0) write(*,*) DPH_IO(1:NCL1_io, 1, 1)
        call wrt_3d_pt_debug (DPH_IO(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'phi', '@sol phi', ITERG, NS) ! debug4chapsim2
        call wrt_3d_all_debug(DPH_IO(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'phi', '@sol phi', ITERG, NS) ! debug4chapsim2
        
        !CALL CHECK_FFT_SOLVER!test
        !CALL DEBUG_WRT_LOCAL(DPH_io,0,N2DO(MYID)+1,'dphi') !test
    
        !=======STEP5: CALCULATE PRESSURE=================
        CALL PRCALC_io(NS)
        call wrt_3d_pt_debug(PR_IO(1:NCL1_io, 1:N2DO(myid), 1:NCL3), 'pr', '@updated', ITERG, NS) ! debug4chapsim2
        
        !=======STEP6: CALCULATE mass flux=================
        CALL MASSFLUX_UPDATE_IO(NS)
        call wrt_3d_pt_debug(G_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 1), 'ux', '@updated', ITERG, NS) ! debug4chapsim2
        call wrt_3d_pt_debug(G_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 2), 'uy', '@updated', ITERG, NS) ! debug4chapsim2
        call wrt_3d_pt_debug(G_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 3), 'uz', '@updated', ITERG, NS) ! debug4chapsim2

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


    SUBROUTINE TEST_POISSON
        use cparam
        use thermal_info 
        use mesh_info
        use flow_info
        use init_info
    implicit none
    integer :: i, j, k

    DO I = 1, NCL1_IO
      DO J = 1, N2DO(myid)
        DO K = 1, NCL3
            RHSLLPHI_io(i, j, k) = - 4._WP * dsin(xcc_io(i)) * dcos(xcc_io(i))
        END DO
      END DO
    END DO 
    IF(myid == 0) then
        do i = 1, NCL1_IO
            write(*,*) 'innput', i, RHSLLPHI_io(i, 8, 8)
        end do 
    end if

    CALL FFT99_POIS3D_periodicxz(IIO) !Method One 

    IF(myid == 0) then
        do i = 1, NCL1_IO
            write(*,*) 'output', i, DPH_IO(i, 8, 8)
        end do 
    end if


    END SUBROUTINE
