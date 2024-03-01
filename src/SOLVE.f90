!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Kui He, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> The main CFD solver
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! ??/??/201?- Add pipe treatment, by Kui He.
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!
!>        @warning : why from ITERG>1 to re-calcuate DT?
!**********************************************************************
    SUBROUTINE SOLVE
        use init_info
        use flow_info 
        use thermal_info
        IMPLICIT NONE      
   
        REAL(WP)    :: ENDTIME(2)
        REAL(WP)    :: STARTTIME
        REAL(WP)    :: REN_TEMP   
       
        INTEGER(4) :: NTSTF       
        INTEGER(4) :: NS
        INTEGER(4) :: ITERL
       
        !=========note for all ranks=================================================    
        !=========INITIAL TIME/STEP SETTING UP=======================================
        
        
        ITERL= 0
        IF(TGFLOWflg) THEN
            ITERG0_TG  = 0
            phyTIME_TG = 0.0_WP
            
            CONVH0_tg = 0.0_WP
        END IF
        
        IF(IOFLOWflg) THEN
            ITERG0_io  = 0
            phyTIME_io = 0.0_WP
            
            EXPLT0_io     = 0.0_WP
            IF(thermlflg==1) RHS_ENERGY0 = 0.0_WP
        END IF
        
        
        !========FLOW INITIALIZATION, EITHER FROM RANDOM OR FROM RESTARTING===========
        CALL FLOWSTART
        REN_TEMP=REN
        DT0 = DT
        IF(IOFLOWflg) THEN
            IF(RSTflg_io==2 .and. RST_time_set_flg==1) THEN
                ITERG0_io  = 0
                phyTIME_io = 0.0_WP  
            END IF
            phyTIME = phyTIME_io
            ITERG0  = ITERG0_io
        ELSE
            phyTIME = phyTIME_TG
            ITERG0  = ITERG0_TG
        END IF
        DTSAVE  = DTSAVE1
        
        !========SCREEN PRINGITNG=====================================================
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL PP_MONITOR_INI
        CALL MPI_BARRIER(ICOMM,IERROR)
        !=========MAIN STEPS ADVANCING=========MAIN SOLVE=============================
        !NTSTF = IDINT(TSTOP/DT)*10000
        NTSTF = 2000000000
        IF(MYID==0) THEN
            CALL CHKINTHDL(' The total iterations to calculation (given)=',myid,NTSTF)
            CALL CHKINTHDL(' The current iterations (restarted)=         ',myid,ITERG0)
        END IF
        DO ITERG=ITERG0+1,NTSTF
            STARTTIME=MPI_WTIME()
          
            !===========SET UP CFL===============
            IF(TGFLOWFLG) CALL CFL_tg
            IF(IOFLOWFLG) CALL CFL_io
            IF(IOFLOWFLG .AND. TGFLOWFLG) THEN
                CFLMM = DMAX1(CFLMM_tg,CFLMM_io)
            ELSE   
                IF(TGFLOWFLG)  CFLMM = CFLMM_tg
                IF(IOFLOWFLG)  CFLMM = CFLMM_io
            END IF
            
            !==========share the same common info for approaching=====
            IF ( phyTIME.LT.TLGRE ) THEN
                REN    = REINI
            ELSE
                REN    = REN_TEMP
            ENDIF
            CALL BCAST_COMM_STEP

            !============STEP AND TIME CONTROL=======================
            ITERL      = ITERL      + 1
            phyTIME    = phyTIME    + DT
            IF(TGFLOWFLG) phyTIME_tg = phyTIME
            IF(IOFLOWFLG) phyTIME_io = phyTIME
            
!            !===========SET UP RECORDING TIME INTERVAL===============
            IF(phyTIME.LT.TSTAV1) THEN
                DTSAVE  = DTSAVE1 * 2.0_WP
            ELSE
                DTSAVE  = DTSAVE1
            END IF
            
            
            !=============RK=====the main part===========================
            DO NS=1, NSST
                IF(TGFLOWFLG) CALL SOLVERRK3_MOM_tg(NS)
                IF(IOFLOWFLG .and. thermlflg==1) THEN
                    CALL SOLVERRK3_ENG_IO(NS)
                    CALL DENSITY_TIME_DERIVATION(NS)
                    CALL DENSITY_Staggered
                    CALL MU_Staggered
                    IF(visthemflg==visimplicit) THEN
                        CALL DENSITY_implicit_add
                    END IF
                END IF
                IF(IOFLOWFLG) THEN
                     !=======STEP6: CALCULATE VELOCITY=================
                    CALL VELOCITY_CALC_io
                    CALL SOLVERRK3_MOM_IO(NS)
                END IF
            END DO
            
            !===================================================
            ENDTIME(1)=MPI_WTIME()
            CPUTIME_TMP=ENDTIME(1)-STARTTIME

            
            !=======POSTPROCESS=============================
            IF( DMOD(phyTIME,DTTECCK) .LT.DT) & 
            CALL CALL_TEC360
            IF(TGFLOWflg) CALL POSTPROCESS_tg
            IF(IOFLOWflg) CALL POSTPROCESS_io

        END DO 

       
        RETURN
    END SUBROUTINE
    
 !*******************************************************************************************************************************
    SUBROUTINE BCAST_COMM_STEP
        use init_info
        use flow_info 
        use thermal_info
        IMPLICIT NONE      
    
        REAL(WP)    :: COMMINFO(5) 
        REAL(WP)    :: DT0HALF
    
        COMMINFO = 0.0_wp
        IF(MYID.EQ.0) THEN
            !============REYNOLDS NUMBER AND PRT0=====================
            
            CVISC  = 1.0_WP/REN
            IF(ioflowflg .and. (thermlflg==1)) CTHECD = 1.0_WP/REN/PRT0
          
            !==========CALCULATE DT FROM GIVEN DT AND CFL============
            IF(phyTIME.LT.TSTAV1) THEN
                DT0HALF = DT0
            ELSE
                DT0HALF = DT0 !* 0.5_wp
            END IF
            
            DT = DMIN1(DT0HALF,5.0_wp*DT,CFLGV/CFLMM) !DT can not exceed five times of the last step
            DT = DMAX1(DT, DTMIN)

            COMMINFO(1) = REN
            COMMINFO(2) = CVISC
            COMMINFO(3) = CTHECD
            COMMINFO(4) = CFLMM
            COMMINFO(5) = DT
        END IF
            
        CALL MPI_BCAST( COMMINFO, 5, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        IF(MYID.NE.0) THEN
            REN    = COMMINFO(1) 
            CVISC  = COMMINFO(2) 
            CTHECD = COMMINFO(3) 
            CFLMM  = COMMINFO(4) 
            DT     = COMMINFO(5) 
        END IF
    
    
        ! IN CASE OF ANY SMALL DT DUE TO DIVERGENCE
        IF((DT/DT0).lt.0.01_wp) THEN
        
            CALL CALL_TEC360
            IF(TGFLOWflg) CALL POSTPROCESS_tg
            IF(IOFLOWflg) CALL POSTPROCESS_io
            
            CALL MPI_BARRIER(ICOMM,IERROR)
            STOP
        END IF
    
    
        RETURN
    END SUBROUTINE
