!***********************************************************************
    SUBROUTINE RANDOM_FL_THEML_FLD_io
        use init_info
        use mesh_info
        use flow_info
        use thermal_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE      
             
        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        INTEGER(4) :: L
        real(wp)    :: Umean1_io
 
        IF(.not.IOFLOWflg) RETURN

        !========Generate scaled random Q and PR=1==========================       
        !    INFLOW/OUTFLOW Domain***
        IF(MYID.EQ.0) CALL CHKHDL('   (1-2) IO: Generating random velocity',myid)
        IF(thermlflg==1) CALL INIFIELD_THERMAL_io  
        CALL INIFIELD_FLOW_io
        CALL VMAV_io
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('         IO: VMV(1:3) In random velocity field',myid)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Max. U, V, W', (VMAX_io(L),L=1,3)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Min. U, V, W', (VMIN_io(L),L=1,3)
        END IF

        !>================ Q=Q-Qxzmean=========================================
        !>    INFLOW/OUTFLOW Domain***
        IF(thermlflg==1) CALL INTFC_MFD_THERMAL_io
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
        CALL BC_WALL_Q_io
        CALL INITIAL_MEAN_IO
        DO L=1,3
            DO J=1,N2DO(MYID)
                Q_IO(:,J,:,L)= Q_IO(:,J,:,L) - UU(J,L,2)
            ENDDO
        END DO 
     
        CALL VMAV_io          
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (2-2) IO: VMV(1:3) in random velocity after subtracting meanXZ ',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Max. U, V, W', (VMAX_io(L),L=1,3)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Min. U, V, W', (VMIN_io(L),L=1,3)
        END IF
      
        IF(ICASE.NE.IBOX3P) THEN !!!added
        !>===========update Q in the incoming flow direction====================  
        DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            DO I=1,NCL1E
                DO K=1,NCL3
                    Q_io(I,J,K,NFLOW)= Q_io(I,J,K,NFLOW) + Vini(JJ)  
                ENDDO
            ENDDO
        ENDDO  
                
        !IF(ICASE.NE.ICHANL) CALL VELO2RVELO_IO
        
        !>    @note Set Q(I,J,K,2)=v at Wall b.c. be zero.    
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
        CALL BC_WALL_Q_io
        
        !==============CORRECT THE MASS FLOW RATE=1=========================================
        CALL BULK_VELOCITY_IO(Umean1_io)
        IF(MYID.EQ.0) call CHKRLHDL  ('        IO: The bulk velocity (original)=        ',MYID,Umean1_io)
       
        DO J=1,N2DO(MYID)
            Q_io(:,J,:,NFLOW)= Q_io(:,J,:,NFLOW)/Umean1_io
        END DO
        
        CALL BULK_VELOCITY_IO(Umean1_io)
        IF(MYID.EQ.0) call CHKRLHDL  ('        IO: The bulk velocity (corrected)=       ',MYID,Umean1_io)
        END IF
        
        
        !>============CHECK DIVERGENCE======================================          
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
        CALL BC_WALL_Q_io

        CALL VELO2MASSFLUX
        IF(TGFLOWFLG) CALL BC_TINLET_FLOW
       
        CALL Unified_massflux
        
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
        CALL BC_WALL_Q_io
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
        CALL BC_WALL_G_io
        CALL VMAV_io          
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (3-2) IO: VMV(1:3) in the real initial flow field ',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Max. U, V, W', (VMAX_io(L),L=1,3)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Min. U, V, W', (VMIN_io(L),L=1,3)
        END IF
       
        CALL DIVGCK_io    
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (3-2) IO: Max divergence of initial flow field. Main Domain, Inlet and Outlet',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==>', MAXDIVGV_io(1:3)
        END IF
        
        
                
                
        !CALL CALL_TEC360
      
        RETURN
      
    END SUBROUTINE
      
!*********************************************************************************************************************
    SUBROUTINE CALC_INITIALIZATION_io
        use init_info
        use mesh_info
        use flow_info
        use thermal_info
        IMPLICIT NONE      
             
        INTEGER(4) :: L
       
        IF(.not.IOFLOWflg) RETURN
        
        
        !======Calculate the QP fields ================================
        IF(MYID.EQ.0) CALL CHKHDL('   (4-2) IO: Initial flow field, calc QP...',myid) 
          
        IF(thermlflg==1) THEN
            CALL BC_WALL_THERMAL(IALLDM)
            CALL SOLVERRK3_ENG_IO(0) 
            CALL DENSITY_Staggered
            CALL MU_Staggered
        END IF
        CALL SOLVERRK3_MOM_IO(0) 
        
        !==============CHECK INFO========================================
        CALL VMAV_io
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (5-2) IO: Initial flow field in the first calcuted velocity, VMAX(1:3)',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Max. U, V, W', (VMAX_io(L),L=1,3)
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==> Min. U, V, W', (VMIN_io(L),L=1,3)
        END IF
        
        CALL DIVGCK_io
        IF(MYID.EQ.0) THEN
            CALL CHKHDL('   (6-2) IO: Max divergence of the flow field. Main Domain, Inlet and Outlet',myid) 
            WRITE(*,'(A,25X,A,3ES23.15)') '#','==>', MAXDIVGV_io(1:3)
        END IF
        !MAXDIVGV_io(1) = MAXDIVGV_io(2)
        
        RETURN
    END SUBROUTINE
!*********************************************************************************************************************       
    SUBROUTINE INIFIELD_FLOW_io  
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
        
        
        IF(.NOT.IOFLOWflg) RETURN
        
        IF(ICASE .EQ. IBOX3P) THEN
            RANDOMTYPE=100
        ELSE
            RANDOMTYPE=2
        END IF
       
        IF (RANDOMTYPE==1) THEN
            XRDM = 0.0_WP  
    !>       @note Generate the random flow field u,v,w.      
            DO J=1,N2DO(MYID)           ! only local y_cell no.      
                JJ=JCL2G(J)              
    !>           @param VPERQ=scaled magnitude of perturbation based on Up (parabolic max. velo.) 
    !>           @note for 1/4 near wall region, perturbation ratio is decreased by 25%. 
                VPERQ=VPERG                               
                IF((1.0_WP-DABS(YCC(JJ))).LT.0.250_WP) VPERQ=VPERG*SVPERG
              
                DO I=1,NCL1E
                    DO K=1,NCL3             
    
                        CALL RANDOM_NUMBER(XRDM)
    !>                  @note Transfer random real no. 0~1 to -1~1                
                        XRD(1)=(-1.0_WP+2.0_WP*XRDM(1))
                        XRD(2)=(-1.0_WP+2.0_WP*XRDM(2))
                        XRD(3)=(-1.0_WP+2.0_WP*XRDM(3))
    !>                  @note the scaled random perturbations, in three directions.                
                        Q_io(I,J,K,1)=VPERQ*XRD(1)
                        Q_io(I,J,K,2)=VPERQ*XRD(2) /RNDI1(JJ)
                        Q_io(I,J,K,3)=VPERQ*XRD(3) /RCCI1(JJ)
                    ENDDO
                ENDDO
                IF(JJ==1) Q_io(:,J,:,2) = 0.0_WP
            ENDDO   
        END IF
        
!******************************************************************
       
        IF (RANDOMTYPE==2) THEN
       
            XRD = 0.0_WP  
            RDCOUNT = 0
    !>     @note Generate the random flow field u,v,w.      
            DO J=1,N2DO(MYID)           ! only local y_cell no.      
                JJ=JCL2G(J)              
    !>          @param VPERQ=scaled magnitude of perturbation based on Up (parabolic max. velo.) 
    !>          @note for 1/4 near wall region, perturbation ratio is decreased by 25%. 
                VPERQ=VPERG                               
                IF((1.0_WP-DABS(YCC(JJ))).LT.0.250_WP) VPERQ=VPERG*SVPERG
              
                DO I=1,NCL1E
                    DO K=1,NCL3
                        !RDCOUNT = RDCOUNT + 1
                        RDCOUNT = I + K + JJ  ! make it independent of processor number.
                        seed = 1973+RDCOUNT+(NCL1_TG*NCL3*NCL2)+1024
                        call random_initialize ( seed )
                        call rvec_random ( -1.0_WP, 1.0_WP, 3, XRD )
    !>                   @note the scaled random perturbations, in three directions.                
                        Q_io(I,J,K,1)=VPERQ*XRD(1)
                        Q_io(I,J,K,2)=VPERQ*XRD(2) /RNDI1(JJ)
                        Q_io(I,J,K,3)=VPERQ*XRD(3) /RCCI1(JJ)
                        !write(*,'(3I4.2,3ES13.5)') I,J,K,XRD(1),XRD(2),XRD(3) !test
    !>               @note No perturbulation is added to pressure                
                    ENDDO
                ENDDO
                IF(JJ==1) Q_IO(:,J,:,2) = 0.0_WP
            ENDDO
        END IF
        PR_io(:,:,:)  = 0.0_WP
        
        !******************************************************************
        IF(ICASE .EQ. IBOX3P .and. RANDOMTYPE==100) THEN

            Q_io =0.0_WP
            PR_io=0.0_WP            
            
            DO I=1,NCL1E
                DO K=1,NCL3
                
                    DO J=1,N2DO(MYID)     
                        JJ=JCL2G(J) 
                        Q_io(I,J,K,1)= DSIN(XND_io(I))*DCOS(YCC(JJ))*DCOS(ZCC(K))
                        Q_io(I,J,K,2)=-DCOS(XCC_io(I))*DSIN(YND(JJ))*DCOS(ZCC(K))
                        Q_io(I,J,K,3)= 0.0_WP
                        PR_io(I,J,K) = (DCOS(2.0_WP*XCC_io(I))+DCOS(2.0_WP*YCC(JJ)))* &
                                       (DCOS(2.0_WP*ZCC(K))+2.0_WP)/16.0_WP
                                       
                        if((i == 8) .and. (jj == 8) .and. (k == 8)) write(*,*) Q_io(I,J,K,1:3), PR_io(I,J,K)
                        !if (myid==0) WRITE(*,*) I, K, JJ, Q_io(I,J,K,1:3),PR_io(I,J,K)
                    END DO
                    
                    !IF(MYID==NPSLV) THEN    
                    !    J=N2DO(MYID)+1
                    !    JJ=NND2
                    !    Q_io(I,J,K,2)=-DCOS(XCC_io(I))*DSIN(YND(JJ))*DCOS(ZCC(K))
                    !END IF
                    
                END DO
            END DO
                    ! TEST===========    
        ! call wrt_3d_pt_debug(Q_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 1), 'ux', '@af init') ! debug4chapsim2
        ! call wrt_3d_pt_debug(Q_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 2), 'uy', '@af init') ! debug4chapsim2
        ! call wrt_3d_pt_debug(Q_IO (1:NCL1_io, 1:N2DO(myid), 1:NCL3, 3), 'uz', '@af init') ! debug4chapsim2
        ! call wrt_3d_pt_debug(PR_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3),    'pr', '@af init') ! debug4chapsim2   
        END IF
        !*****************************************************************************
        
            

        IF(weightedpressure==1) THEN
            PR0_io(:,:,:,1)  = PR_io(:,:,:)
            PR0_io(:,:,:,2)  = PR_io(:,:,:)
        END IF
                   
        
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io) 
        CALL INTFC_VARS1(NCL1S,NCL1E,NCL1S,NCL1E,PR_io) 
                 
            
        CALL VELO2MASSFLUX
        
        !CALL DEBUG_WRT_QGP_io
 
        RETURN
    END SUBROUTINE
       
       
!*********************************************************************************************************************       
    SUBROUTINE INIFIELD_THERMAL_io
!>    @note
!>          Set up the initial thermal field.
!>    Known:  
!>           P_0, T_0 
!>           P_0 is only used to creat the NIST table, and not used in the code.
!>    To set up: 
!>           T, \rho, \mu, \kappa, h, h\rho 
        USE MESH_INFO
        USE thermal_info
        use init_info
        IMPLICIT NONE
        
        INTEGER(4) :: I, J, K, JJ, IE, IS
        REAL(WP)   :: T_temp, SFS, SFE, H_temp, D_temp, K_temp, M_temp
      
        NTHERMAL = 6
        
        
        RHOH(:,:,:)        = 0.0_WP
        ENTHALPY(:,:,:)    = 0.0_WP
        DENSITY(:,:,:)     = 1.0_WP
        TEMPERATURE(:,:,:) = 1.0_WP
        VISCOUSITY(:,:,:)  = 1.0_WP
        THERMCONDT(:,:,:)  = 1.0_WP
        
        RHOH0(:,:,:) = 0.0_wp
        
        RETURN
    END SUBROUTINE 
      
      
    !*********************************************************************************************************************   
    SUBROUTINE BULK_VELOCITY_IO(UUMEAN_work)
        use init_info
        use mesh_info
        use flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE      
             
        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        REAL(WP)    :: UUMEAN, UUMEAN_work
        
        UUMEAN = 0.0_wp
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            DO I=1,NCL1_IO
                DO K=1,NCL3
                    UUMEAN=UUMEAN + Q_IO(I,J,K,NFLOW)/DYFI(JJ)/RCCI1(JJ)
                END DO
            END DO
        ENDDO
        UUMEAN=UUMEAN/DZI
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(UUMEAN, UUMEAN_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        
        UUMEAN_work = UUMEAN_work/(Area_INLET*DBLE(NCL1_io))
    
        RETURN
    END SUBROUTINE
    
!*********************************************************************************************************************   
    SUBROUTINE BULK_MASSFLUX_IO(GMEAN_work)
        use init_info
        use mesh_info
        use flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE      
             
        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        REAL(WP)    :: GMEAN, GMEAN_work
        
        GMEAN = 0.0_wp
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            DO I=1,NCL1_IO
                DO K=1,NCL3
                    GMEAN=GMEAN + G_IO(I,J,K,NFLOW)/DYFI(JJ)/RCCI1(JJ)
                END DO
            END DO
        ENDDO
        GMEAN=GMEAN/DZI
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(GMEAN, GMEAN_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        
        GMEAN_work = GMEAN_work/(Area_INLET*DBLE(NCL1_io))
    
        RETURN
    END SUBROUTINE
    
!*********************************************************************************************************************   
    SUBROUTINE BULK_H_IO(HMEAN_work)
        use init_info
        use mesh_info
        use flow_info
        use thermal_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE      
             
        INTEGER(4) :: I
        INTEGER(4) :: J, JJ
        INTEGER(4) :: K
        REAL(WP)    :: HMEAN, HMEAN_work
        
        HMEAN = 0.0_wp
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            DO I=1,NCL1_IO
                DO K=1,NCL3
                    HMEAN=HMEAN + ENTHALPY(I,J,K)/DYFI(JJ)/RCCI1(JJ)
                END DO
            END DO
        ENDDO
        HMEAN=HMEAN/DZI
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(HMEAN, HMEAN_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        
        HMEAN_work = HMEAN_work/(Area_INLET*DBLE(NCL1_io))
    
        RETURN
    END SUBROUTINE    

    
    
!*********************************************************************************************************
    SUBROUTINE VELO2MASSFLUX  
        USE MESH_INFO
        USE THERMAL_INFO
        USE FLOW_INFO
        USE init_info
        IMPLICIT NONE
        INTEGER(4) :: IDR, NII, IC,JC,KC, IM, JJ, KM, JM, KS
        !REAL(wp)    :: RHO
      
      
!>      !=====================================X=============================
        IDR = 1
        DO IC=1,NCL1_io
            IM = IMV_io(IC)
            DO JC=1,N2DO(MYID) 
                DO KC=1,NCL3     
                    !RHO = ( DENSITY(IC,JC,KC) + DENSITY(IM,JC,KC) ) * XND2CL
                    G_io(IC,JC,KC,IDR)=Q_io(IC,JC,KC,IDR)*RHO_STG(IC,JC,KC,IDR)!*RHO
                ENDDO
            ENDDO
        ENDDO
      
!>      !=====================================Z=============================
        IDR = 3      
        DO KC=1,NCL3
            KM=KMV(KC)
            DO JC=1,N2DO(MYID)
                DO IC=1,NCL1_io
                    !RHO = ( DENSITY(IC,JC,KC) + DENSITY(IC,JC,KM) ) * ZND2CL
                    G_io(IC,JC,KC,IDR)=Q_io(IC,JC,KC,IDR)*RHO_STG(IC,JC,KC,IDR)!*RHO
                ENDDO
            ENDDO
        ENDDO
      
!>      !=====================================Y=============================
        IDR = 2  
        NII = 1
        IF(MYID==0 .and. icase.ne.IBOX3P) NII = 2 
        DO JC=NII,N2DO(MYID)
            JJ=JCL2G(JC)
            JM=JLMV(JC)   
            DO KC=1,NCL3
                DO IC=1,NCL1_io
                    !RHO = YCL2ND_WFF(JJ) * DENSITY(IC,JC,KC) + &
                    !      YCL2ND_WFB(JJ) * DENSITY(IC,JM,KC)
                    G_io(IC,JC,KC,IDR)=Q_io(IC,JC,KC,IDR)*RHO_STG(IC,JC,KC,IDR)!*RHO
                ENDDO
            ENDDO
        ENDDO 
                       
        IF (MYID.EQ.0 .and. icase.ne.IBOX3P) THEN
            DO KC=1,NCL3
                KS =KSYM(KC)
                DO IC=NCL1S,NCL1E
                    IF(ICASE.EQ.ICHANL) THEN
                        G_io(IC,1,KC,IDR)=0.0_WP
                        G_io(IC,0,KC,IDR)=0.0_WP ! not used!
                    END IF
                    IF(ICASE.EQ.IPIPEC) THEN
                        G_io(IC,1,KC,IDR)=0.0_WP
                        G_io(IC,0,KC,IDR)=G_io(IC,1,KS,IDR)
                    END IF
                ENDDO
            ENDDO
        ENDIF
      
        IF (MYID.EQ.NPSLV) THEN
            DO KC=1,NCL3
                DO IC=NCL1S,NCL1E
                    IF(ICASE.EQ.IBOX3P) THEN
                        G_io(IC,N2DO(MYID)+1,KC,IDR)=Q_io(IC,N2DO(MYID)+1,KC,IDR)*RHO_STG(IC,N2DO(MYID)+1,KC,IDR)
                    ELSE
                        G_io(IC,N2DO(MYID)+1,KC,IDR)=0.0_WP
                    END IF
                ENDDO
            ENDDO
        ENDIF
      
        RETURN
    END SUBROUTINE
    
    
!    SUBROUTINE INTFC_G_WALL
!        USE MESH_INFO
!        USE THERMAL_INFO
!        USE FLOW_INFO
!        USE init_info
!        IMPLICIT NONE
!        INTEGER(4) :: IDR, NII, IC,JC,KC, IM, JJ, KM, JM, KS,JP
!        REAL(wp)    :: RHO1,RHO2, RHO0
!        REAL(WP)    :: RHO_ND2, RHO_ND0, RHO_WAL, RHO_WL0
        
!        IF(MYID==0) THEN
!            JC=1
!            JP=0
!        ELSE IF(MYID==NPSLV) THEN
!            JC=N2DO(MYID)
!            JP=N2DO(MYID)+1
!        ELSE
!            RETURN
!        END IF
        
!        !=====================================X=============================
!        IDR = 1
!        DO IC=NCL1S,NCL1E
!            IF(IC==0) THEN
!                IM=0
!            ELSE
!                IM=IMV_io(IC)
!            END IF
           
!            DO KC=1,NCL3     
!                RHO1 = ( DENSITY(IC,JC,KC) + DENSITY(IM,JC,KC) ) * XND2CL
!                RHO2 = ( DENSITY(IC,JP,KC) + DENSITY(IM,JP,KC) ) * XND2CL
!                Q_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)/RHO1
!                G_io(IC,JP,KC,IDR)=Q_io(IC,JC,KC,IDR)*RHO2*(-1.0_wp)
!            ENDDO
            
!            !Q_io(:,JP,:,IDR)=-Q_io(:,JC,:,IDR) !test
!            !G_io(:,JP,:,IDR)=-G_io(:,JC,:,IDR) !test
            
!        ENDDO
      
!!>      !=====================================Z=============================
!        IDR = 3      
!        DO KC=1,NCL3
!            KM=KMV(KC)
!            DO IC=NCL1S,NCL1E
!                RHO1 = ( DENSITY(IC,JC,KC) + DENSITY(IC,JC,KM) ) * ZND2CL
!                RHO2 = ( DENSITY(IC,JP,KC) + DENSITY(IC,JP,KM) ) * ZND2CL
!                Q_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)/RHO1
!                G_io(IC,JP,KC,IDR)=Q_io(IC,JC,KC,IDR)*RHO2*(-1.0_wp)
!            ENDDO
            
!            !Q_io(:,JP,:,IDR)=-Q_io(:,JC,:,IDR) !test
!            !G_io(:,JP,:,IDR)=-G_io(:,JC,:,IDR) !test
            
!        ENDDO
        
!!>      !=====================================Y=============================
!        IF(MYID==0) THEN
!            IDR = 2      
!            DO KC=1,NCL3
!                KM=KMV(KC)
!                DO IC=NCL1S,NCL1E
                
!                    RHO_WAL = ( DENSITY (IC,0,KC) + DENSITY(IC,1,KC)  ) * 0.5_wp
!                    RHO_WL0 = ( DENSITY0(IC,0,KC) + DENSITY0(IC,1,KC) ) * 0.5_wp
!                    RHO_ND2 =   DENSITY(IC,1,KC)*YCL2ND_WFF(2) + DENSITY(IC,2,KC)*YCL2ND_WFF(2)
!					RHO_ND0 = 2.0_wp*RHO_WAL-RHO_ND2

!                    Q_io(IC,2,KC,IDR)=G_io(IC,2,KC,IDR)/RHO_ND2
!                    Q_io(IC,0,KC,IDR)=(-1.0_wp)*Q_io(IC,2,KC,IDR)!-2.0_wp/DT*(RHO_WAL-RHO_WL0)/RHO_WAL
!                    G_io(IC,0,KC,IDR)=Q_io(IC,0,KC,IDR)*RHO_ND0
!                ENDDO
!            ENDDO
            
!            !Q_io(:,0:1,:,IDR)=0.0_wp !test
!            !G_io(:,0:1,:,IDR)=0.0_wp !test
            
!        END IF
        
!        IF(MYID==NPSLV) THEN
!            IDR = 2      
!            Q_io(:,N2DO(MYID)+1,:,IDR)=0.0_wp
!            G_io(:,N2DO(MYID)+1,:,IDR)=0.0_wp
!        END IF
      
!        RETURN
!    END SUBROUTINE
      
!*********************************************************************      
    SUBROUTINE Unified_massflux
        USE MESH_INFO
        USE THERMAL_INFO
        USE FLOW_INFO
        USE init_info
        IMPLICIT NONE
        
        REAL(WP)  :: GMEAN1, GMEAN2, Umean1_io
        
        IF(ICASE.EQ.IBOX3P) RETURN
        
        GMEAN1 = 0.0_wp
        CALL BULK_MASSFLUX_IO(GMEAN1)
        IF(MYID==0) CALL CHKRLHDL('              The initial mass flux is  ',myid,GMEAN1)
        
        
        Q_IO(:,:,:,1:3) = Q_IO(:,:,:,1:3)/GMEAN1
        CALL BULK_VELOCITY_IO(Umean1_io)
        IF(MYID==0) CALL CHKRLHDL('              The unified Umean is      ',myid,Umean1_io)
        
        
        CALL VELO2MASSFLUX
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        GMEAN2 = 0.0_wp
        CALL BULK_MASSFLUX_IO(GMEAN2)
        IF(MYID==0) CALL CHKRLHDL('              The unified mass flux is  ',myid,GMEAN2)
    
        RETURN
    END SUBROUTINE

!!***********************************************************************************    
!    SUBROUTINE VELO2RVELO_IO
!        use init_info
!        use mesh_info
!        use flow_info
!        IMPLICIT NONE    
        
!        INTEGER(4) :: I, J, K, NYI, JJ
        
!        NYI = 1
!        IF(MYID==0) THEN
!            J = 1
!            DO I=1, NCL1_IO
!                DO K=1,NCL3
!                    Q_IO(I,J,K,1) = Q_IO(I,J,K,1) 
!                    Q_IO(I,J,K,2) = 0.0_WP
!                    Q_IO(I,J,K,3) = Q_IO(I,J,K,3) / RCCI1(JJ)
!                END DO
!            END DO
!        NYI = 2
!        END IF
        
!        DO J=NYI,N2DO(MYID)
!            JJ = JCL2G(J)
!            DO I=1, NCL1_IO
!                DO K=1,NCL3
!                    Q_IO(I,J,K,1) = Q_IO(I,J,K,1) 
!                    Q_IO(I,J,K,2) = Q_IO(I,J,K,2) / RNDI1(JJ)
!                    Q_IO(I,J,K,3) = Q_IO(I,J,K,3) / RCCI1(JJ)
!                END DO
!            END DO
!        END DO
        
!        RETURN
!    END SUBROUTINE
