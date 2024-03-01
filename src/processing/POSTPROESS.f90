    SUBROUTINE POSTPROCESS_TG
        use init_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE
       
        
        !========================MONITORING THE CALCULATION RUNNING PROPERLY========================
        CALL DIVGCK_tg
        
        !=============== call average every time step after the setting tsav1=========================
        IF(phyTime .GT. TSTAV1) THEN
            CALL PP_MEAN_ZX_FLOW_TG
            CALL PP_MEAN_T_TG(T_Asymptotic_Average)
        END IF
        
        !===============screen printing for monitoring the code======================================
        IF (ITERG==1 .OR. MOD(ITERG,ITSCN)==0)  THEN
            IF(phyTime .LE. TSTAV1) THEN
                CALL PP_MEAN_ZX_FLOW_TG
            END IF
            CALL PP_MONITOR_TG
        END IF
        
        !=================write out data============================================================
        
        IF( DMOD(phyTIME,DTSAVE) .LT.DT) THEN
            !================write out basic variables=================
            CALL WRT_INSTANT_VARS_TG
            !================write out time averaged variables===========
            IF(phyTime .GT. TSTAV1) THEN
                CALL WRT_AVERAGE_VARS_TG
                IF(DMOD(phyTIME,DTSTATEC) .LT. DT) THEN
                    CALL WRT_AVERAGE_PPED_TG
                END IF
            END IF
            
        END IF
       
    RETURN
    END SUBROUTINE
      
!*******************************************************************

    SUBROUTINE POSTPROCESS_IO
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
       
        
        !============instantanous saving===================================
        IF( DMOD(phyTIME,DTSAVE) .LT.DT) CALL WRT_INSTANT_VARS_io
        
        !=============== call average every time step after the setting tsav1=========================
        IF(phyTime .GT. TSTAV1) THEN
        
            !===================instantanous processing for each step==========================
            IF(TGFLOWflg) THEN
                CALL PP_MEAN_Z_FLOW_nonXperiodic_IO
                IF(thermlflg==1) &
                CALL PP_MEAN_Z_THEML_nonXperiodic_IO
                CALL PP_MEAN_T_nonXperiodic_IO(T_Asymptotic_Average)
            ELSE
                CALL PP_MEAN_ZX_FLOW_Xperiodic_IO
                CALL PP_SPECOSPEC
                CALL PP_MEAN_T_Xperiodic_IO(T_Asymptotic_Average)
            END IF

            !==================writing out========================
            IF( DMOD(phyTIME,DTSAVE) .LT.DT) THEN
            
                IF(TGFLOWflg) THEN
                    CALL WRT_AVERAGE_VARS_nonXperiodic_IO
                    IF(DMOD(phyTIME,DTSTATEC) .LT. DT) THEN
                        CALL WRT_AVERAGE_PPED_nonXperiodic_IO
                    END IF
                ELSE
                    CALL WRT_AVERAGE_VARS_Xperiodic_IO
                    IF(DMOD(phyTIME,DTSTATEC) .LT. DT) THEN
                        CALL WRT_AVERAGE_PPED_Xperiodic_IO
                        IF(ppspectra==1) CALL PP_SPECOSPEC
                    END IF
                END IF
            END IF    
                
        END IF
        
        !===============screen printing for monitoring the code======================================
        IF ( (ITERG==1) .OR. MOD(ITERG,ITSCN)==0 )  THEN
            !========================MONITORING THE CALCULATION RUNNING PROPERLY==
            CALL DIVGCK_io
        
            IF(TGFLOWflg) THEN
            
                IF(phyTime .LE. TSTAV1) THEN
                    CALL PP_MEAN_Z_FLOW_nonXperiodic_IO
                    IF(thermlflg==1) &
                    CALL PP_MEAN_Z_THEML_nonXperiodic_IO
                END IF
                CALL PP_MONITOR_nonXperiodic_IO
                
            ELSE
            
                IF(phyTime .LE. TSTAV1) THEN
                    CALL PP_MEAN_ZX_FLOW_Xperiodic_IO
                END IF
                CALL PP_MONITOR_Xperiodic_IO
                
            END IF
                
        END IF
        
       RETURN
    END SUBROUTINE
      
!*******************************************************************
