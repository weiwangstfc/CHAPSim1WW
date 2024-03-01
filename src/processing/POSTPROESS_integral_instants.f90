    SUBROUTINE POSTPROCESS_INTEGRAL_INSTANS
        use init_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE
    
    
        INTEGER(4) :: N
        REAL(WP)   :: TI
    
        IF(pprocessonly.NE.2) RETURN
        
        IF(MYID.EQ.0) CALL CHKHDL('   IO: re-postprocessing all instantanous flows',myid)
        
        !===========FOR TG PART==========================
        IF(TGFLOWflg) THEN
            !==initialize the averaged values======
                CALL PP_TMEAN_INI_ZERO_TG
                
                !==read in all data=======
                DO N=1, pp_instn_sz
                    ti = pp_instn_tim(N)
                    CALL RESTART_INSTANT_VARS_TG(ti)
                    CALL PP_MEAN_ZX_FLOW_TG
                    CALL PP_MEAN_T_TG(T_Summing_average)
                END DO
                
                !====write data out===
                CALL WRT_AVERAGE_VARS_TG
                CALL WRT_AVERAGE_PPED_TG
                IF(ppspectra==1) CALL PP_SPECOSPEC
            
        
        
        END IF
        
        !==========FOR MAIN PART========================
        
        IF(IOFLOWflg) THEN
            !==initialize the averaged values======
            CALL PP_TMEAN_INI_ZERO_IO
            
            IF(TGFLOWflg) THEN
               ! to add....
            
            ELSE  !===periodic main io domain=======================
                !==read in all data=======
                DO N=1, pp_instn_sz
                    ti = pp_instn_tim(N)
                    CALL RESTART_INSTANT_VARS_IO(ti)
                    CALL PP_MEAN_ZX_FLOW_Xperiodic_IO
                    CALL PP_SPECOSPEC
                    CALL PP_MEAN_T_Xperiodic_IO(T_Summing_average)
                END DO
                
                phyTIME = phyTIME_io
                ITERG0  = ITERG0_io
                    
                !====write data out===
                CALL WRT_AVERAGE_VARS_Xperiodic_IO
                CALL WRT_AVERAGE_PPED_Xperiodic_IO
                !IF(ppspectra==1) CALL PP_SPECOSPEC
                
            END IF
        
        END IF
    

    RETURN
    END SUBROUTINE
    
    
    
    
    
    
    
    
