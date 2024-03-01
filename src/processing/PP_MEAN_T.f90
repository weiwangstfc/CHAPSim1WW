!#########################
!  Notice: by Wei Wang on 9 April, 2005
!          consider below two methods to do a smooth average
!           -moving average
!           -Exponential smoothing

!*************************************************************************************************    
    SUBROUTINE PP_TMEAN_INI_TG
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE
        
        IF(RSTflg_tg.EQ.2 .and. phyTIME_TG .GT. TSTAV1) THEN ! RESTART
        
            CALL RESTART_AVERAGE_VARS_TG
            CALL WRT_AVERAGE_PPED_TG
            
        ELSE
        
            CALL PP_TMEAN_INI_ZERO_TG
            
        END IF
        
        RETURN
    END SUBROUTINE 
    
    SUBROUTINE PP_TMEAN_INI_ZERO_TG
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE
    
        NSTATIS_tg = 0
            
        U1xztL_tg = 0.0_WP
        UPxztL_tg = 0.0_WP
        
        U2xztL_tg = 0.0_WP
        U3xztL_tg = 0.0_WP
        
        DVDL1xztL_tg = 0.0_WP
        DVDLPxztL_tg = 0.0_WP
        DVDL2xztL_tg = 0.0_WP
        
        QUADUVxztL_tg = 0.0_WP
        QUADVzxztL_tg = 0.0_WP
        QUADTKxztL_tg = 0.0_WP
        QUADDRxztL_tg = 0.0_WP
            
        RETURN
    END SUBROUTINE 
!*************************************************************************************************     
    SUBROUTINE PP_TMEAN_INI_IO
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
        
        ! RSTflg_io==2 : restarting from the last step
        ! phyTIME_IO .GT. TSTAV1 : to do time averaging
        ! RST_type_flg==0 : Restart both, flow and thermal fields
        ! RST_time_set_flg==0 :  Restart following previous time
        IF( (RSTflg_io==2) .and. (phyTIME_IO .GT. TSTAV1) .and. (RST_type_flg==0) .and. (RST_time_set_flg==0) ) THEN ! RESTART
            IF(TGFLOWFLG) THEN
                !=============pp space averaged=======
                !CALL PP_MEAN_Z_FLOW_nonXperiodic_IO
                !IF(thermlflg==1) &
                !CALL PP_MEAN_Z_THEML_nonXperiodic_IO
                !============read in space and time averaged=========
                CALL RESTART_AVERAGE_VARS_nonXperiodic_IO
                !============write out space and time averaged=======
                CALL WRT_AVERAGE_PPED_nonXperiodic_IO
                CALL MPI_BARRIER(ICOMM,IERROR)
                !CALL PP_MONITOR_nonXperiodic_IO
            ELSE
                !=============pp space averaged=======
                !CALL PP_MEAN_ZX_FLOW_Xperiodic_IO    
                !(This will be replaced by below) 
                !============read in space and time averaged=========
                CALL RESTART_AVERAGE_VARS_Xperiodic_IO
                IF(TSTAV_RESET.GT.TSTAV1+1.0E-2_WP) THEN
                    phyTIME = RSTtim_io
                    CALL WRT_AVERAGE_VARS_Xperiodic_IO
                END IF
                !============write out space and time averaged=======
                CALL WRT_AVERAGE_PPED_Xperiodic_IO     
                IF(ppspectra==1) THEN
                    CALL PP_MEAN_ZX_FLOW_Xperiodic_IO
                    CALL PP_SPECOSPEC
                END IF
                CALL MPI_BARRIER(ICOMM,IERROR)
                !CALL PP_MONITOR_Xperiodic_IO
            END IF
            
        ELSE
        
            CALL PP_TMEAN_INI_ZERO_IO
            
        END IF
        
    RETURN
    END SUBROUTINE 
    
    
    SUBROUTINE PP_TMEAN_INI_ZERO_IO
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
        
        
        NSTATIS_IO = 0
            
        IF(TGFLOWFLG) THEN
            U1ztL_io = 0.0_WP
            G1ztL_io = 0.0_WP
                
            UPztL_io  = 0.0_WP
            
            U2ztL_IO  = 0.0_WP
            UGztL_IO  = 0.0_WP
            UGUztL_IO = 0.0_WP
            
            DVDL1ztL_io = 0.0_WP
            DVDLPztL_io = 0.0_WP
            DVDL2ztL_io = 0.0_WP
            
            IF(thermlflg==1) THEN    
                T1ztL_io  = 0.0_WP
                D1ztL_io  = 0.0_WP
                H1ztL_io  = 0.0_WP 
             
                T2ztL_io  = 0.0_WP
                D2ztL_io  = 0.0_WP
                H2ztL_io  = 0.0_WP
            
                UHztL_io  = 0.0_wp
                GHztL_io  = 0.0_wp
            END IF
        
        
        ELSE
        
            U1xztL_io =  0.0_WP
            G1xztL_io =  0.0_WP
            UPxztL_io =  0.0_WP
            
            U2xztL_IO =  0.0_WP
            UGxztL_IO =  0.0_WP
            UGUxztL_IO=  0.0_WP
            U3xztL_IO =  0.0_WP
            
            DVDL1xztL_io = 0.0_WP
            DVDLPxztL_io = 0.0_WP
            DVDL2xztL_io = 0.0_WP
            
            QUADUVxztL_io = 0.0_WP
            QUADVzxztL_io = 0.0_WP
            QUADTKxztL_io = 0.0_WP
            QUADDRxztL_io = 0.0_WP
            
            QUADDUV1xztL_io = 0.0_WP
            QUADDUV2xztL_io = 0.0_WP
            
            FUxztL_IO = 0.0_WP
            

            IF(thermlflg==1) THEN
                T1xztL_io  = 0.0_WP
                D1xztL_io  = 0.0_WP
                H1xztL_io  = 0.0_WP
                M1xztL_io  = 0.0_WP
                        
                T2xztL_io  = 0.0_WP
                D2xztL_io  = 0.0_WP
                H2xztL_io  = 0.0_WP
                        
                DHxztL_io  = 0.0_WP
                PHxztL_io  = 0.0_WP
                
                DVDL1MxztL_io = 0.0_WP
                DVDL1MHxztL_io= 0.0_WP
                DVDL1MUxztL_io= 0.0_WP
                DVDL2MxztL_io = 0.0_WP
                
                UHxztL_io  = 0.0_WP
                GHxztL_io  = 0.0_WP
                U2DHxztL_io= 0.0_WP
                
                DhDL1xztL_io  = 0.0_WP
                DhDLPxztL_io  = 0.0_WP
                DTDLKxztL_io  = 0.0_WP
                DTDLKUxztL_io = 0.0_WP
                DTDLKDVDLxztL_io = 0.0_WP 
                DHDLMDVDLxztL_io = 0.0_WP
            END IF

            if(myid==0) then
                !===========spectra======================
                !================Velocity =====================
                R11X1_xztLa=0.0_wp
                R22X1_xztLa=0.0_wp
                R33X1_xztLa=0.0_wp
                R12X1_xztLa=0.0_wp
                R13X1_xztLa=0.0_wp
                R23X1_xztLa=0.0_wp
                
                R11X3_xztLa=0.0_wp
                R22X3_xztLa=0.0_wp
                R33X3_xztLa=0.0_wp
                R12X3_xztLa=0.0_wp
                R13X3_xztLa=0.0_wp
                R23X3_xztLa=0.0_wp
                
                ENE11T_xztLa=0.0_wp
                ENE22T_xztLa=0.0_wp
                ENE33T_xztLa=0.0_wp
                ENE12T_xztLa=0.0_wp
                ENE13T_xztLa=0.0_wp
                ENE23T_xztLa=0.0_wp
                
                ENE11Z_xztLa=0.0_wp
                ENE22Z_xztLa=0.0_wp
                ENE33Z_xztLa=0.0_wp
                ENE12Z_xztLa=0.0_wp
                ENE13Z_xztLa=0.0_wp
                ENE23Z_xztLa=0.0_wp
                
                !================Voriticity ====================
                V11X1_xztLa=0.0_wp
                V22X1_xztLa=0.0_wp
                V33X1_xztLa=0.0_wp
                V12X1_xztLa=0.0_wp
                V13X1_xztLa=0.0_wp
                V23X1_xztLa=0.0_wp
                
                V11X3_xztLa=0.0_wp
                V22X3_xztLa=0.0_wp
                V33X3_xztLa=0.0_wp
                V12X3_xztLa=0.0_wp
                V13X3_xztLa=0.0_wp
                V23X3_xztLa=0.0_wp
                
                ENV11T_xztLa=0.0_wp
                ENV22T_xztLa=0.0_wp
                ENV33T_xztLa=0.0_wp
                ENV12T_xztLa=0.0_wp
                ENV13T_xztLa=0.0_wp
                ENV23T_xztLa=0.0_wp
                
                ENV11Z_xztLa=0.0_wp
                ENV22Z_xztLa=0.0_wp
                ENV33Z_xztLa=0.0_wp
                ENV12Z_xztLa=0.0_wp
                ENV13Z_xztLa=0.0_wp
                ENV23Z_xztLa=0.0_wp
                
                
                !===============Voriticity & Velocity=========================
                VO11X1_xztLa=0.0_wp
                VO12X1_xztLa=0.0_wp
                VO13X1_xztLa=0.0_wp
                
                VO21X1_xztLa=0.0_wp
                VO22X1_xztLa=0.0_wp
                VO23X1_xztLa=0.0_wp
                
                VO31X1_xztLa=0.0_wp
                VO32X1_xztLa=0.0_wp
                VO33X1_xztLa=0.0_wp
                
                VO11X3_xztLa=0.0_wp
                VO12X3_xztLa=0.0_wp
                VO13X3_xztLa=0.0_wp
                
                VO21X3_xztLa=0.0_wp
                VO22X3_xztLa=0.0_wp
                VO23X3_xztLa=0.0_wp
                
                VO31X3_xztLa=0.0_wp
                VO32X3_xztLa=0.0_wp
                VO33X3_xztLa=0.0_wp
                
                EVO11T_xztLa=0.0_wp
                EVO12T_xztLa=0.0_wp
                EVO13T_xztLa=0.0_wp
                
                EVO21T_xztLa=0.0_wp
                EVO22T_xztLa=0.0_wp
                EVO23T_xztLa=0.0_wp
                
                EVO31T_xztLa=0.0_wp
                EVO32T_xztLa=0.0_wp
                EVO33T_xztLa=0.0_wp
                
                EVO11Z_xztLa=0.0_wp
                EVO12Z_xztLa=0.0_wp
                EVO13Z_xztLa=0.0_wp
                
                EVO21Z_xztLa=0.0_wp
                EVO22Z_xztLa=0.0_wp
                EVO23Z_xztLa=0.0_wp
                
                EVO31Z_xztLa=0.0_wp
                EVO32Z_xztLa=0.0_wp
                EVO33Z_xztLa=0.0_wp
            
            end if
        END IF
        
    RETURN
    END SUBROUTINE     
    
!*************************************************************************************************     
    SUBROUTINE PP_MEAN_T_TG(TP)
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: TP
        REAL(WP) :: COE1, COE2, COE3
        
        NSTATIS_TG = NSTATIS_TG + 1
        
        IF(TP==T_Asymptotic_Average) THEN  ! gradually increasing the averaged accuracy
            !COE1 = DBLE(NSTATIS_IO - 1)
            COE2 = 1.0_WP/DBLE(NSTATIS_TG)
            COE3 = DBLE(NSTATIS_TG - 1)/DBLE(NSTATIS_TG)
        ELSE IF(TP==T_Summing_average) THEN ! just adding all them together
            !COE1 = 1.0_WP/
            COE2 = 1.0_WP/DBLE(pp_instn_sz)
            COE3 = 1.0_WP
        ELSE
        END IF
        
        
!        NSTATIS_TG = NSTATIS_TG + 1
        
!        COE1 = DBLE(NSTATIS_TG - 1)
!        COE2 = 1.0_WP/DBLE(NSTATIS_TG)
!        COE3 = DBLE(NSTATIS_TG - 1)/DBLE(NSTATIS_TG)
        
        
        U1xztL_tg =  COE3 * U1xztL_tg + COE2*U1xzL_tg
        UPxztL_tg =  COE3 * UPxztL_tg + COE2*UPxzL_tg
            
        U2xztL_tg =  COE3 * U2xztL_tg + COE2*U2xzL_tg
        U3xztL_tg =  COE3 * U3xztL_tg + COE2*U3xzL_tg
            
        DVDL1xztL_tg = COE3 * DVDL1xztL_tg + COE2*DVDL1xzL_tg
        DVDLPxztL_tg = COE3 * DVDLPxztL_tg + COE2*DVDLPxzL_tg
        DVDL2xztL_tg = COE3 * DVDL2xztL_tg + COE2*DVDL2xzL_tg
        
        QUADUVxztL_tg= COE3 * QUADUVxztL_tg+ COE2*QUADUVxzL_tg
        QUADVzxztL_tg= COE3 * QUADVzxztL_tg+ COE2*QUADVzxzL_tg
        QUADTKxztL_tg= COE3 * QUADTKxztL_tg+ COE2*QUADTKxzL_tg
        QUADDRxztL_tg= COE3 * QUADDRxztL_tg+ COE2*QUADDRxzL_tg
        
        
    RETURN
    END SUBROUTINE 
    
!*************************************************************************************************  
    SUBROUTINE PP_MEAN_T_nonXperiodic_IO(TP)
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: TP
        REAL(WP) :: COE1, COE2, COE3
        
        
        NSTATIS_IO = NSTATIS_IO + 1
        
        IF(TP==T_Asymptotic_Average) THEN  ! gradually increasing the averaged accuracy
            !COE1 = DBLE(NSTATIS_IO - 1)
            COE2 = 1.0_WP/DBLE(NSTATIS_IO)
            COE3 = DBLE(NSTATIS_IO - 1)/DBLE(NSTATIS_IO)
        ELSE IF(TP==T_Summing_average) THEN ! just adding all them together
            !COE1 = 1.0_WP/
            COE2 = 1.0_WP/DBLE(pp_instn_sz)
            COE3 = 1.0_WP
        ELSE
        END IF
        

        !COE1 = DBLE(NSTATIS_IO - 1)
        !COE2 = 1.0_WP/DBLE(NSTATIS_IO)
        !COE3 = DBLE(NSTATIS_IO - 1)/DBLE(NSTATIS_IO)
        
        U1ztL_io =  COE3* U1ztL_io + COE2*U1zL_io
        G1ztL_io =  COE3* G1ztL_io + COE2*G1zL_io
        
        UPztL_io =  COE3* UPztL_io + COE2*UPzL_io
        
        U2ztL_IO =  COE3* U2ztL_IO  + COE2*U2zL_IO
        UGztL_IO =  COE3* UGztL_IO  + COE2*UGzL_IO
        UGUztL_IO=  COE3* UGUztL_IO + COE2*UGUzL_IO
            
        DVDL1ztL_io = COE3* DVDL1ztL_io + COE2*DVDL1zL_io
        DVDLPztL_io = COE3* DVDLPztL_io + COE2*DVDLPzL_io
        DVDL2ztL_io = COE3* DVDL2ztL_io + COE2*DVDL2zL_io
        
        IF(thermlflg==1) THEN
            T1ztL_io  = COE3*T1ztL_io + COE2*T1zL_io
            D1ztL_io  = COE3*D1ztL_io + COE2*D1zL_io
            H1ztL_io  = COE3*H1ztL_io + COE2*H1zL_io
                    
            T2ztL_io  = COE3*T2ztL_io + COE2*T2zL_io
            D2ztL_io  = COE3*D2ztL_io + COE2*D2zL_io
            H2ztL_io  = COE3*H2ztL_io + COE2*H2zL_io
                    
            DHztL_io  = COE3*DHztL_io + COE2*DHzL_io
                    
            UHztL_io  = COE3*UHztL_io + COE2*UHzL_io
            GHztL_io  = COE3*GHztL_io + COE2*GHzL_io
        END IF
        
        
    RETURN
    END SUBROUTINE 
    
    
    SUBROUTINE PP_MEAN_T_Xperiodic_IO(TP)
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: TP
        REAL(WP) :: COE2, COE3
        
        
        NSTATIS_IO = NSTATIS_IO + 1
        
        IF(TP==T_Asymptotic_Average) THEN  ! gradually increasing the averaged accuracy
            !COE1 = DBLE(NSTATIS_IO - 1)
            COE2 = 1.0_WP/DBLE(NSTATIS_IO)
            COE3 = DBLE(NSTATIS_IO - 1)/DBLE(NSTATIS_IO)
        ELSE IF(TP==T_Summing_average) THEN ! just adding all them together
            !COE1 = 1.0_WP/
            COE2 = 1.0_WP/DBLE(pp_instn_sz)
            COE3 = 1.0_WP
        ELSE
        END IF
        
        U1xztL_io =  COE3* U1xztL_io + COE2* U1xzL_io
        G1xztL_io =  COE3* G1xztL_io + COE2* G1xzL_io
        UPxztL_io =  COE3* UPxztL_io + COE2* UPxzL_io
        
        U2xztL_IO =  COE3* U2xztL_IO + COE2* U2xzL_IO
        UGxztL_IO =  COE3* UGxztL_IO + COE2* UGxzL_IO
        
        U3xztL_IO =  COE3* U3xztL_IO + COE2* U3xzL_IO
        UGUxztL_IO=  COE3* UGUxztL_IO+ COE2* UGUxzL_IO
            
        DVDL1xztL_io = COE3* DVDL1xztL_io + COE2* DVDL1xzL_io
        DVDLPxztL_io = COE3* DVDLPxztL_io + COE2* DVDLPxzL_io
        DVDL2xztL_io = COE3* DVDL2xztL_io + COE2* DVDL2xzL_io
        
        QUADUVxztL_io= COE3 * QUADUVxztL_io+ COE2*QUADUVxzL_io
        QUADVzxztL_io= COE3 * QUADVzxztL_io+ COE2*QUADVzxzL_io
        QUADTKxztL_io= COE3 * QUADTKxztL_io+ COE2*QUADTKxzL_io
        QUADDRxztL_io= COE3 * QUADDRxztL_io+ COE2*QUADDRxzL_io
        QUADDUV1xztL_io= COE3 * QUADDUV1xztL_io+ COE2*QUADDUV1xzL_io
        QUADDUV2xztL_io= COE3 * QUADDUV2xztL_io+ COE2*QUADDUV2xzL_io
        
        OCTDUVxztL_io= COE3 * OCTDUVxztL_io+ COE2*OCTDUVxzL_io
        OCTDVzxztL_io= COE3 * OCTDVzxztL_io+ COE2*OCTDVzxzL_io
        OCTDTKxztL_io= COE3 * OCTDTKxztL_io+ COE2*OCTDTKxzL_io
        OCTDDRxztL_io= COE3 * OCTDDRxztL_io+ COE2*OCTDDRxzL_io
        OCTDDUV1xztL_io= COE3 * OCTDDUV1xztL_io+ COE2*OCTDDUV1xzL_io
        OCTDDUV2xztL_io= COE3 * OCTDDUV2xztL_io+ COE2*OCTDDUV2xzL_io
        
        OCTTUVxztL_io= COE3 * OCTTUVxztL_io+ COE2*OCTTUVxzL_io
        OCTTVzxztL_io= COE3 * OCTTVzxztL_io+ COE2*OCTTVzxzL_io
        OCTTTKxztL_io= COE3 * OCTTTKxztL_io+ COE2*OCTTTKxzL_io
        OCTTDRxztL_io= COE3 * OCTTDRxztL_io+ COE2*OCTTDRxzL_io
        OCTTDUV1xztL_io= COE3 * OCTTDUV1xztL_io+ COE2*OCTTDUV1xzL_io
        OCTTDUV2xztL_io= COE3 * OCTTDUV2xztL_io+ COE2*OCTTDUV2xzL_io
        
        FUxztL_IO      = COE3 * FUxztL_IO      + COE2*FUxzL_IO
        
        
        IF(thermlflg==1) THEN
            T1xztL_io  = COE3* T1xztL_io + COE2* T1xzL_io
            D1xztL_io  = COE3* D1xztL_io + COE2* D1xzL_io
            H1xztL_io  = COE3* H1xztL_io + COE2* H1xzL_io
            M1xztL_io  = COE3* M1xztL_io + COE2* M1xzL_io
                    
            T2xztL_io  = COE3* T2xztL_io + COE2* T2xzL_io
            D2xztL_io  = COE3* D2xztL_io + COE2* D2xzL_io
            H2xztL_io  = COE3* H2xztL_io + COE2* H2xzL_io
                    
            DHxztL_io  = COE3* DHxztL_io + COE2* DHxzL_io
            PHxztL_io  = COE3* PHxztL_io + COE2* PHxzL_io
            
            DVDL1MxztL_io = COE3* DVDL1MxztL_io + COE2* DVDL1MxzL_io
            DVDL1MHxztL_io= COE3* DVDL1MHxztL_io+ COE2* DVDL1MHxzL_io
            DVDL1MUxztL_io= COE3* DVDL1MUxztL_io+ COE2* DVDL1MUxzL_io
            DVDL2MxztL_io = COE3* DVDL2MxztL_io + COE2* DVDL2MxzL_io
            
            UHxztL_io  = COE3* UHxztL_io   + COE2* UHxzL_io
            GHxztL_io  = COE3* GHxztL_io   + COE2* GHxzL_io
            U2DHxztL_io= COE3* U2DHxztL_io + COE2* U2DHxzL_io
            
            DhDL1xztL_io  = COE3* DhDL1xztL_io + COE2* DhDL1xzL_io
            DhDLPxztL_io  = COE3* DhDLPxztL_io + COE2* DhDLPxzL_io
            DTDLKxztL_io  = COE3* DTDLKxztL_io + COE2* DTDLKxzL_io
            DTDLKUxztL_io = COE3* DTDLKUxztL_io+ COE2* DTDLKUxzL_io
            
            DTDLKDVDLxztL_io = COE3* DTDLKDVDLxztL_io+ COE2* DTDLKDVDLxzL_io   
            DHDLMDVDLxztL_io = COE3* DHDLMDVDLxztL_io+ COE2* DHDLMDVDLxzL_io    
        END IF
        
        if(myid==0) then
            !===========spectra======================
            !================Velocity =====================
            R11X1_xztLa = COE3 * R11X1_xztLa + COE2 *R11X1_xzLa
            R22X1_xztLa = COE3 * R22X1_xztLa + COE2 *R22X1_xzLa
            R33X1_xztLa = COE3 * R33X1_xztLa + COE2 *R33X1_xzLa
            R12X1_xztLa = COE3 * R12X1_xztLa + COE2 *R12X1_xzLa
            R13X1_xztLa = COE3 * R13X1_xztLa + COE2 *R13X1_xzLa
            R23X1_xztLa = COE3 * R23X1_xztLa + COE2 *R23X1_xzLa
            
            R11X3_xztLa = COE3 * R11X3_xztLa + COE2 *R11X3_xzLa
            R22X3_xztLa = COE3 * R22X3_xztLa + COE2 *R22X3_xzLa
            R33X3_xztLa = COE3 * R33X3_xztLa + COE2 *R33X3_xzLa
            R12X3_xztLa = COE3 * R12X3_xztLa + COE2 *R12X3_xzLa
            R13X3_xztLa = COE3 * R13X3_xztLa + COE2 *R13X3_xzLa
            R23X3_xztLa = COE3 * R23X3_xztLa + COE2 *R23X3_xzLa
            
            ENE11T_xztLa = COE3 * ENE11T_xztLa + COE2 *ENE11T_xzLa
            ENE22T_xztLa = COE3 * ENE22T_xztLa + COE2 *ENE22T_xzLa
            ENE33T_xztLa = COE3 * ENE33T_xztLa + COE2 *ENE33T_xzLa
            ENE12T_xztLa = COE3 * ENE12T_xztLa + COE2 *ENE12T_xzLa
            ENE13T_xztLa = COE3 * ENE13T_xztLa + COE2 *ENE13T_xzLa
            ENE23T_xztLa = COE3 * ENE23T_xztLa + COE2 *ENE23T_xzLa
            
            ENE11Z_xztLa = COE3 * ENE11Z_xztLa + COE2 *ENE11Z_xzLa
            ENE22Z_xztLa = COE3 * ENE22Z_xztLa + COE2 *ENE22Z_xzLa
            ENE33Z_xztLa = COE3 * ENE33Z_xztLa + COE2 *ENE33Z_xzLa
            ENE12Z_xztLa = COE3 * ENE12Z_xztLa + COE2 *ENE12Z_xzLa
            ENE13Z_xztLa = COE3 * ENE13Z_xztLa + COE2 *ENE13Z_xzLa
            ENE23Z_xztLa = COE3 * ENE23Z_xztLa + COE2 *ENE23Z_xzLa
            
            !================Voriticity ====================
            V11X1_xztLa = COE3 * V11X1_xztLa + COE2 *V11X1_xzLa
            V22X1_xztLa = COE3 * V22X1_xztLa + COE2 *V22X1_xzLa
            V33X1_xztLa = COE3 * V33X1_xztLa + COE2 *V33X1_xzLa
            V12X1_xztLa = COE3 * V12X1_xztLa + COE2 *V12X1_xzLa
            V13X1_xztLa = COE3 * V13X1_xztLa + COE2 *V13X1_xzLa
            V23X1_xztLa = COE3 * V23X1_xztLa + COE2 *V23X1_xzLa
            
            V11X3_xztLa = COE3 * V11X3_xztLa + COE2 *V11X3_xzLa
            V22X3_xztLa = COE3 * V22X3_xztLa + COE2 *V22X3_xzLa
            V33X3_xztLa = COE3 * V33X3_xztLa + COE2 *V33X3_xzLa
            V12X3_xztLa = COE3 * V12X3_xztLa + COE2 *V12X3_xzLa
            V13X3_xztLa = COE3 * V13X3_xztLa + COE2 *V13X3_xzLa
            V23X3_xztLa = COE3 * V23X3_xztLa + COE2 *V23X3_xzLa
            
            ENV11T_xztLa = COE3 * ENV11T_xztLa + COE2 *ENV11T_xzLa
            ENV22T_xztLa = COE3 * ENV22T_xztLa + COE2 *ENV22T_xzLa
            ENV33T_xztLa = COE3 * ENV33T_xztLa + COE2 *ENV33T_xzLa
            ENV12T_xztLa = COE3 * ENV12T_xztLa + COE2 *ENV12T_xzLa
            ENV13T_xztLa = COE3 * ENV13T_xztLa + COE2 *ENV13T_xzLa
            ENV23T_xztLa = COE3 * ENV23T_xztLa + COE2 *ENV23T_xzLa
            
            ENV11Z_xztLa = COE3 * ENV11Z_xztLa + COE2 *ENV11Z_xzLa
            ENV22Z_xztLa = COE3 * ENV22Z_xztLa + COE2 *ENV22Z_xzLa
            ENV33Z_xztLa = COE3 * ENV33Z_xztLa + COE2 *ENV33Z_xzLa
            ENV12Z_xztLa = COE3 * ENV12Z_xztLa + COE2 *ENV12Z_xzLa
            ENV13Z_xztLa = COE3 * ENV13Z_xztLa + COE2 *ENV13Z_xzLa
            ENV23Z_xztLa = COE3 * ENV23Z_xztLa + COE2 *ENV23Z_xzLa
            
            
            !===============Voriticity & Velocity=========================
            VO11X1_xztLa = COE3 * VO11X1_xztLa + COE2 *VO11X1_xzLa
            VO12X1_xztLa = COE3 * VO12X1_xztLa + COE2 *VO12X1_xzLa
            VO13X1_xztLa = COE3 * VO13X1_xztLa + COE2 *VO13X1_xzLa
            
            VO21X1_xztLa = COE3 * VO21X1_xztLa + COE2 *VO21X1_xzLa
            VO22X1_xztLa = COE3 * VO22X1_xztLa + COE2 *VO22X1_xzLa
            VO23X1_xztLa = COE3 * VO23X1_xztLa + COE2 *VO23X1_xzLa
            
            VO31X1_xztLa = COE3 * VO31X1_xztLa + COE2 *VO31X1_xzLa
            VO32X1_xztLa = COE3 * VO32X1_xztLa + COE2 *VO32X1_xzLa
            VO33X1_xztLa = COE3 * VO33X1_xztLa + COE2 *VO33X1_xzLa
            
            VO11X3_xztLa = COE3 * VO11X3_xztLa + COE2 *VO11X3_xzLa
            VO12X3_xztLa = COE3 * VO12X3_xztLa + COE2 *VO12X3_xzLa
            VO13X3_xztLa = COE3 * VO13X3_xztLa + COE2 *VO13X3_xzLa
            
            VO21X3_xztLa = COE3 * VO21X3_xztLa + COE2 *VO21X3_xzLa
            VO22X3_xztLa = COE3 * VO22X3_xztLa + COE2 *VO22X3_xzLa
            VO23X3_xztLa = COE3 * VO23X3_xztLa + COE2 *VO23X3_xzLa
            
            VO31X3_xztLa = COE3 * VO31X3_xztLa + COE2 *VO31X3_xzLa
            VO32X3_xztLa = COE3 * VO32X3_xztLa + COE2 *VO32X3_xzLa
            VO33X3_xztLa = COE3 * VO33X3_xztLa + COE2 *VO33X3_xzLa
            
            EVO11T_xztLa = COE3 * EVO11T_xztLa + COE2 *EVO11T_xzLa
            EVO12T_xztLa = COE3 * EVO12T_xztLa + COE2 *EVO12T_xzLa
            EVO13T_xztLa = COE3 * EVO13T_xztLa + COE2 *EVO13T_xzLa
            
            EVO21T_xztLa = COE3 * EVO21T_xztLa + COE2 *EVO21T_xzLa
            EVO22T_xztLa = COE3 * EVO22T_xztLa + COE2 *EVO22T_xzLa
            EVO23T_xztLa = COE3 * EVO23T_xztLa + COE2 *EVO23T_xzLa
            
            EVO31T_xztLa = COE3 * EVO31T_xztLa + COE2 *EVO31T_xzLa
            EVO32T_xztLa = COE3 * EVO32T_xztLa + COE2 *EVO32T_xzLa
            EVO33T_xztLa = COE3 * EVO33T_xztLa + COE2 *EVO33T_xzLa
            
            EVO11Z_xztLa = COE3 * EVO11Z_xztLa + COE2 *EVO11Z_xzLa
            EVO12Z_xztLa = COE3 * EVO12Z_xztLa + COE2 *EVO12Z_xzLa
            EVO13Z_xztLa = COE3 * EVO13Z_xztLa + COE2 *EVO13Z_xzLa
            
            EVO21Z_xztLa = COE3 * EVO21Z_xztLa + COE2 *EVO21Z_xzLa
            EVO22Z_xztLa = COE3 * EVO22Z_xztLa + COE2 *EVO22Z_xzLa
            EVO23Z_xztLa = COE3 * EVO23Z_xztLa + COE2 *EVO23Z_xzLa
            
            EVO31Z_xztLa = COE3 * EVO31Z_xztLa + COE2 *EVO31Z_xzLa
            EVO32Z_xztLa = COE3 * EVO32Z_xztLa + COE2 *EVO32Z_xzLa
            EVO33Z_xztLa = COE3 * EVO33Z_xztLa + COE2 *EVO33Z_xzLa
        end if
        
    RETURN
    END SUBROUTINE
