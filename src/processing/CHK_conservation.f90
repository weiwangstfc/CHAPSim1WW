!**********************************************************************
    SUBROUTINE CHK_MassConsv_io
        use init_info
        use mesh_info
        use flow_info
        USE thermal_info
        IMPLICIT NONE
      
        REAL(WP) :: FluxINTG_INL,  FluxINTG_INL_WORK
        REAL(WP) :: FluxINTG_OUT,  FluxINTG_OUT_WORK
        REAL(WP) :: MASS_RATE_INTG,MASS_RATE_INTG_WORK
        REAL(WP) :: MASS_RATE, coe2
        REAL(WP) :: dummy(3), dummy_work(3)

        INTEGER(4) :: IC, JC, KC, JJ
               
        !========INLET/OUTLET===========================
        FluxINTG_INL   = 0.0_WP
        FluxINTG_OUT   = 0.0_WP   
        IF(TGFLOWflg) THEN
            DO JC=1,N2DO(myid)
                JJ=JCL2G(JC)
                COE2 = 1.0_WP/DYFI(JJ)/RCCI1(JJ)
                DO KC=1,NCL3
                    FluxINTG_INL   =FluxINTG_INL   + G_io(1,    JC,KC,NFLOW)*COE2
                    FluxINTG_OUT   =FluxINTG_OUT   + G_io(NCL1E,JC,KC,NFLOW)*COE2
                ENDDO
            ENDDO
            FluxINTG_INL   = FluxINTG_INL/DZI
            FluxINTG_OUT   = FluxINTG_OUT/DZI
        END IF
        
        !========INSIDE===========================
        MASS_RATE_INTG = 0.0_WP        
        DO JC=1,N2DO(myid)
            JJ=JCL2G(JC)
            COE2 = 1.0_WP/DYFI(JJ)/RCCI1(JJ)
            DO KC=1,NCL3
                DO IC = 1, NCL1_IO
                    !MASS_RATE      = ( DENSITY(IC,JC,KC) - DENSITYP(IC,JC,KC) )*COE2
                    MASS_RATE      =  DrhoDtP(IC,JC,KC)*COE2
                    MASS_RATE_INTG =  MASS_RATE + MASS_RATE_INTG
                END DO
            ENDDO
        ENDDO
        MASS_RATE_INTG = MASS_RATE_INTG/DZI/DXI
        
        !===========ADD TOGETHER FROM ALL RANKS===============================
        dummy(1) = FluxINTG_INL
        dummy(2) = FluxINTG_OUT
        dummy(3) = MASS_RATE_INTG
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(dummy,dummy_WORK,3,MPI_DOUBLE_PRECISION,MPI_SUM, ICOMM, IERROR)     
        FluxINTG_INL_WORK   = dummy_WORK(1)
        FluxINTG_OUT_WORK   = dummy_WORK(2)
        MASS_RATE_INTG_WORK = dummy_WORK(3)
                             
        !! method 1
        IF(TGFLOWflg) &
        CHK_Mass_CONSV0_SCALING = (FluxINTG_INL_WORK - MASS_RATE_INTG_WORK) / FluxINTG_OUT_WORK
        !CHK_MASS_CONSV0         = DABS(CHK_Mass_CONSV0_SCALING - 1.0_wp)
        
        !! method 2
        CHK_MASS_CONSV0 = FluxINTG_INL_WORK - MASS_RATE_INTG_WORK - FluxINTG_OUT_WORK
        
        !IF(MYID==0) write(*,'(A,I3.1, 5ES13.5)') 'mass coe=', myid, CHK_MASS_CONSV0,CHK_Mass_CONSV0_SCALING, &
        !FluxINTG_INL_WORK, MASS_RATE_INTG_WORK, FluxINTG_OUT_WORK
        
          
        RETURN
    END SUBROUTINE
    
!**********************************************************************
    

!**********************************************************************
    SUBROUTINE CHK_EnegConsv_io
        use init_info
        use mesh_info
        use flow_info
        use thermal_info
        IMPLICIT NONE
      
        REAL(WP) :: FluxINTG_INL,  FluxINTG_INL_WORK
        REAL(WP) :: FluxINTG_OUT,  FluxINTG_OUT_WORK
        REAL(WP) :: ENEG_RATE_INTG,ENEG_RATE_INTG_WORK, eneg_total, eneg_total_work
        
        REAL(WP) :: RHOHV_BTW, KDT_BTW
        REAL(WP) :: RHOHV_TPW, KDT_TPW
        REAL(WP) :: FluxINTG_BTW, FluxINTG_BTW_WORK
        REAL(WP) :: FluxINTG_TPW, FluxINTG_TPW_WORK
        REAL(WP) :: RHOHU_INL, RHOHU_OUT, KDT_INL, KDT_OUT, WHEATFLUX
        INTEGER(4) :: IC, JC, KC, JJ, JJP, JM, JP, IC1, IC2, IM1, IP2
        REAL(WP) :: COE1, COE2, COE221, COE222
        REAL(wp) :: dummy(6), dummy_work(6)
        
        
        
        !CALL CHK_EnegConsv_each_cell_io ! CHECK
        
        IF(thermlflg.ne.1) RETURN
        
        !==============INLET/OUTLET AT Y-Z PLANE=======================
        FluxINTG_INL = 0.0_wp
        FluxINTG_OUT = 0.0_wp
        IF(TGFLOWflg) THEN
            IC1 = 1
            IM1 = IMV_IO(IC1)
            
            IC2 = NCL1_io
            IP2 = IPV_IO(IC2)
            
            DO JC=1,N2DO(myid)
                JJ=JCL2G(JC)
                COE1 = XND2CL/DYFI(JJ)/RCCI1(JJ)
                COE2 = 1.0_WP/DYFI(JJ)/RCCI1(JJ)
                DO KC=1,NCL3
                    RHOHU_INL = ( ENTHALPY(IM1,JC,KC) + ENTHALPY(IC1,JC,KC) ) * G_io(IC1,JC,KC,1)
                    KDT_INL   = ( THERMCONDT(IC1,JC,KC) + THERMCONDT(IM1,JC,KC) ) * &
                                ( TEMPERATURE(IC1,JC,KC) - TEMPERATURE(IM1,JC,KC) ) * DXI * CTHECD
                                  
                    RHOHU_OUT = ( ENTHALPY(IC2,JC,KC) + ENTHALPY(IP2,JC,KC) ) * G_io(IP2,JC,KC,1)
                    KDT_OUT   = ( THERMCONDT(IP2,JC,KC) + THERMCONDT(IC2,JC,KC) ) * &
                                ( TEMPERATURE(IP2,JC,KC) - TEMPERATURE(IC2,JC,KC) ) * DXI * CTHECD
                                  
                    FluxINTG_INL  = FluxINTG_INL   + (RHOHU_INL-KDT_INL)*COE1
                    FluxINTG_OUT  = FluxINTG_OUT   + (RHOHU_OUT-KDT_OUT)*COE1
                                
                ENDDO
            ENDDO
            FluxINTG_INL =  FluxINTG_INL / DZI
            FluxINTG_OUT = -FluxINTG_OUT / DZI
        END IF
        
        !=================INSIDE=================================
        ENEG_RATE_INTG=0.0_wp
        eneg_total=0.0_wp
        DO JC=1,N2DO(myid)
            JJ=  JCL2G(JC)
            COE2 = 1.0_WP/DYFI(JJ)/RCCI1(JJ)
            DO KC=1,NCL3         
                DO IC = 1, NCL1_IO
                    ENEG_RATE_INTG =  (RHOH(IC,JC,KC)-RHOH0(IC,JC,KC))*COE2 + ENEG_RATE_INTG
                    eneg_total=eneg_total+RHOH(IC,JC,KC)*COE2
                END DO
            ENDDO
        ENDDO
        ENEG_RATE_INTG=ENEG_RATE_INTG/DZI/DXI/DT
        eneg_total=eneg_total/DZI/DXI
        
        !==============BOTTOM WALL AT X-Z PLANE=======================
        FluxINTG_BTW = 0.0_wp
        IF(MYID==0) THEN
            JC = 1
            JP = JLPV(JC)
            JM = JLMV(JC)
            JJ=  JCL2G(JC)
            JJP= JGPV(JJ)
            
            IF( BCWALLHEAT(ibotwall)==isoFluxWall    ) THEN
                IF(ICASE==IPIPEC) THEN
                    COE222= 0.0_wp
                ELSE
                    COE222= 1.0_wp /RNDI1(JJ)* CTHECD
                END IF
            END IF
            
            IF( BCWALLHEAT(ibotwall)==isoThermalWall    ) THEN
                IF(ICASE==IPIPEC) THEN
                    COE222= 0.0_wp
                ELSE
                    COE222= DYCI(JJ)/RNDI1(JJ)* CTHECD
                END IF
            END IF
                
            DO IC=1,NCL1_io
                DO KC=1,NCL3
                    RHOHV_BTW =  0.0_WP
                    !KDT_BTW   =  WHEATFLUX* CTHECD
                    
                    IF( BCWALLHEAT(ibotwall)==isoFluxWall    ) KDT_BTW =  -WALLFLUX(IC,1) * COE222
                    IF( BCWALLHEAT(ibotwall)==isoThermalWall ) THEN
                        KDT_BTW =  ( YCL2ND_WFF(JJ)  * THERMCONDT(IC,JC,KC) +   &
                                     YCL2ND_WFB(JJ)  * THERMCONDT(IC,JM,KC) ) * &
                                   ( TEMPERATURE(IC,JC,KC) - TEMPERATURE(IC,JM,KC) ) * COE222        !DYCI(JJ)/RNDI(JJ)
                    END IF

                    FluxINTG_BTW  = FluxINTG_BTW   + (RHOHV_BTW - KDT_BTW) 
                ENDDO
            ENDDO
            FluxINTG_BTW = FluxINTG_BTW/DXI/DZI
            IF(ICASE==IPIPEC) FluxINTG_BTW = 0.0_WP
        END IF
        
        !==============TOP WALL AT X-Z PLANE=======================
        FluxINTG_TPW = 0.0_wp
        IF(MYID==NPSLV) THEN
            JC = N2DO(MYID)
            JP = JLPV(JC)
            JM = JLMV(JC)
            JJ=  JCL2G(JC)
            JJP= JGPV(JJ)
            
            IF( BCWALLHEAT(itopwall)==isoFluxWall    ) THEN
                COE221= 1.0_wp /RNDI1(JJP)* CTHECD
            END IF
            IF( BCWALLHEAT(itopwall)==isoThermalWall    ) THEN
                COE221= DYCI(JJP) /RNDI1(JJP)* CTHECD
            END IF
                
            DO IC=1,NCL1_io
                
                DO KC=1,NCL3
                    RHOHV_TPW =  0.0_WP
                    IF( BCWALLHEAT(itopwall)==isoFluxWall    ) KDT_TPW =  -WALLFLUX(IC,2) * COE221
                    IF( BCWALLHEAT(itopwall)==isoThermalWall ) THEN
                        KDT_TPW =  ( YCL2ND_WFF(JJP) * THERMCONDT(IC,JP,KC) +   &
                                     YCL2ND_WFB(JJP) * THERMCONDT(IC,JC,KC) ) * & 
                                   ( TEMPERATURE(IC,JP,KC) - TEMPERATURE(IC,JC,KC) ) * COE221   !DYCI(JJP)/RNDI(JJP)
                    END IF
                        
                    FluxINTG_TPW  = FluxINTG_TPW   + (RHOHV_TPW + KDT_TPW) !positive
                ENDDO
            ENDDO
            FluxINTG_TPW = FluxINTG_TPW/DXI/DZI/RNDI1(NND2)
        END IF
        
        dummy(1) = FluxINTG_INL
        dummy(2) = FluxINTG_OUT
        dummy(3) = ENEG_RATE_INTG
        dummy(4) = FluxINTG_BTW
        dummy(5) = FluxINTG_TPW
        dummy(6) = ENEG_total
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(dummy,  dummy_work,  6, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR) 
        FluxINTG_INL_WORK    = dummy_WORK(1)
        FluxINTG_OUT_WORK    = dummy_WORK(2)
        ENEG_RATE_INTG_WORK  = dummy_WORK(3)
        FluxINTG_BTW_WORK    = dummy_WORK(4)
        FluxINTG_TPW_WORK    = dummy_WORK(5)
        ENEG_total_WORK      = dummy_WORK(6)
        
        !CALL MPI_ALLREDUCE(FluxINTG_INL,  FluxINTG_INL_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)    
        !CALL MPI_ALLREDUCE(FluxINTG_OUT,  FluxINTG_OUT_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        !CALL MPI_ALLREDUCE(ENEG_RATE_INTG,ENEG_RATE_INTG_WORK,1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)           
        !CALL MPI_ALLREDUCE(FluxINTG_BTW,  FluxINTG_BTW_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)   
        !CALL MPI_ALLREDUCE(FluxINTG_TPW,  FluxINTG_TPW_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        !CALL MPI_ALLREDUCE(ENEG_total,ENEG_total_WORK,1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)      
                             
        CHK_ENEG_CONSV0 = ENEG_RATE_INTG_WORK + FluxINTG_INL_WORK + FluxINTG_OUT_WORK + FluxINTG_BTW_WORK + FluxINTG_TPW_WORK
        CHK_ENEG_TOTAL  = ENEG_total_WORK
        !IF(MYID==0) write(*,*) 'eng coe=',CHK_ENEG_CONSV0,ENEG_RATE_INTG_WORK, &
        !-FluxINTG_BTW_WORK,FluxINTG_TPW_WORK, FluxINTG_TPW_WORK-FluxINTG_BTW_WORK
        
        RETURN
    END SUBROUTINE


!**********************************************************************
    SUBROUTINE CHK_EnegConsv_each_cell_io
        use init_info
        use mesh_info
        use flow_info
        use thermal_info
        IMPLICIT NONE
      
       
        
        INTEGER(4) :: IC, IP, IM
        INTEGER(4) :: KC, KP, KM
        INTEGER(4) :: JC, JP, JM, JJ, JJP
      
        REAL(WP) :: TOTLENG, ENEG_RATE_INTG, G1HP, G1HN, G2HP, G2HN, G3HP, G3HN, KXTP, KXTN, KZTP, KZTN, KYTP, KYTN
        REAL(WP) :: G1H, G2H, G3H, KXT, KYT, KZT
        
        INTEGER(4)            :: FLID
        CHARACTER(4) :: NNN
        
        
        IF(thermlflg.ne.1) RETURN
        
        WRITE(NNN,'(1I4.4)') myid
        FLID = MYID+100
        OPEN(FLID,FILE='OUT_VARS'//NNN//'.debug',position="append")
        if(myid==0) write(FLID,*) '#===check energy equation==='

        DO JC=1, N2DO(MYID)
            JP = JLPV(JC)
            JM = JLMV(JC)
            JJ=  JCL2G(JC)
            JJP= JGPV(JJ)
            DO IC=1, NCL1_io
                IP = IPV_io(IC)
                IM = IMV_io(IC)
                DO KC=1, NCL3
                    KP = KPV(KC)
                    KM = KMV(KC)
                    
                    !===RHOH=============================
                    ENEG_RATE_INTG =  (RHOH(IC,JC,KC)-RHOH0(IC,JC,KC))/DYFI(JJ)/DZI/DXI/DT
                    
                    !===G1H===============
                    G1HP = ( ENTHALPY(IP,JC,KC) + ENTHALPY(IC,JC,KC) ) * G_io(IP,JC,KC,1)/DYFI(JJ)/DZI
                    G1HN = ( ENTHALPY(IC,JC,KC) + ENTHALPY(IM,JC,KC) ) * G_io(IC,JC,KC,1)/DYFI(JJ)/DZI
              
                    !===G2H===============
                    G2HP = ( YCL2ND_WFF(JJP)*ENTHALPY(IC,JP,KC) + &
                             YCL2ND_WFB(JJP)*ENTHALPY(IC,JC,KC) ) * G_io(IC,JP,KC,2)/DZI/DXI
                    G2HN = ( YCL2ND_WFF(JJ) *ENTHALPY(IC,JC,KC) +                         &
                             YCL2ND_WFB(JJ) *ENTHALPY(IC,JM,KC) ) * G_io(IC,JC,KC,2)/DZI/DXI 
                
                    !===G3H===============
                    G3HP = ( ENTHALPY(IC,JC,KP) + ENTHALPY(IC,JC,KC) ) * G_io(IC,JC,KP,3)/DYFI(JJ)/DXI
                    G3HN = ( ENTHALPY(IC,JC,KC) + ENTHALPY(IC,JC,KM) ) * G_io(IC,JC,KC,3)/DYFI(JJ)/DXI
                    
                    !========KDTX=========================
                    KXTP = ( THERMCONDT(IP,JC,KC)+THERMCONDT(IC,JC,KC) ) * 0.5_WP* &
                           ( TEMPERATURE(IP,JC,KC)-TEMPERATURE(IC,JC,KC) )*DXI/DYFI(JJ)/DZI*CTHECD
                    KXTN = ( THERMCONDT(IC,JC,KC)+THERMCONDT(IM,JC,KC) ) * 0.5_WP* &
                           ( TEMPERATURE(IC,JC,KC)-TEMPERATURE(IM,JC,KC) )*DXI/DYFI(JJ)/DZI*CTHECD
                           
                    !========KDTZ=========================
                    KZTP = ( THERMCONDT(IC,JC,KP)+THERMCONDT(IC,JC,KC) )  * 0.5_WP* &
                           ( TEMPERATURE(IC,JC,KP)-TEMPERATURE(IC,JC,KC) )*DZI/DYFI(JJ)/DXI*CTHECD
                    KZTN = ( THERMCONDT(IC,JC,KM)+THERMCONDT(IC,JC,KC) ) * 0.5_WP* &
                           ( TEMPERATURE(IC,JC,KC)-TEMPERATURE(IC,JC,KM) )*DZI/DYFI(JJ)/DXI*CTHECD
                           
                    !========KDTY=========================
                    KYTP = ( YCL2ND_WFF(JJP) * THERMCONDT(IC,JP,KC) +   &
                             YCL2ND_WFB(JJP) * THERMCONDT(IC,JC,KC) ) * & 
                           ( TEMPERATURE(IC,JP,KC) - TEMPERATURE(IC,JC,KC) ) *DYCI(JJP)/DZI/DXI*CTHECD
                    KYTN = ( YCL2ND_WFF(JJ)  * THERMCONDT(IC,JC,KC) +   &
                             YCL2ND_WFB(JJ)  * THERMCONDT(IC,JM,KC) ) * &
                           ( TEMPERATURE(IC,JC,KC) - TEMPERATURE(IC,JM,KC) ) *DYCI(JJ)/DZI/DXI*CTHECD
                    IF(JJ==1 .and. BCWALLHEAT(Ibotwall)==isoFluxWall) THEN
                        KYTN =  -WALLFLUX(IC,Ibotwall) /DZI/DXI*CTHECD
                    ELSE IF(JJ==NCL2 .and. BCWALLHEAT(Itopwall)==isoFluxWall) THEN
                        KYTP =  -WALLFLUX(IC,Itopwall) /DZI/DXI*CTHECD
                    ELSE
                    END IF
                    
                    G1H = G1HN-G1HP
                    G2H = G2HN-G2HP
                    G3H = G3HN-G3HP
                    KXT = KXTP-KXTN
                    KYT = KYTP-KYTN 
                    KZT = KZTP-KZTN
                    
                    TOTLENG =-ENEG_RATE_INTG+G1H+G2H+G3H+KXT+KZT+KYT
                    IF( dabs(TOTLENG).gt.1.0E-5_wp) then
                        WRITE(FLID,*) JJ, IC, KC, &
                        TOTLENG, -ENEG_RATE_INTG,G1H+G2H+G3H,KXT+KyT+KzT
                        WRITE(FLID,*) G1HP, G1HN, G1H
                        WRITE(FLID,*) G2HP, G2HN, G2H
                        WRITE(FLID,*) G3HP, G3HN, G3H
                        WRITE(FLID,*) KXTP, KXTN, KXT
                        WRITE(FLID,*) KYTP, KYTN, KYT
                        WRITE(FLID,*) KZTP, KZTN, KZT
                    END IF
                END DO
            END DO
        END DO
        
       
        RETURN
    END SUBROUTINE
