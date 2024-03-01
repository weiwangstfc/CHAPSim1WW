    SUBROUTINE PP_MONITOR_INI
        use init_info
        use mesh_info
        use wrt_info
        USE thermal_info
        IMPLICIT NONE
      
        real(wp) :: dxplus, yMINplus, dzplus, L2, HY, dyMIN, dyMAX, yMAXplus
        
        IF(MYID.NE.0) RETURN
        
        !========================MESH INFO==================================================================
        IF(ICASE==ICHANL .or. ICASE==IBOX3P) THEN
            dyMIN = 1.0_wp/DYFI(1)
            dyMAX = 1.0_WP/DYFI(NCL2/2)
            HY    = 2.0_WP*ALX2
        ELSE
            dyMIN = 1.0_wp/DYFI(NCL2)
            dyMAX = 1.0_WP/DYFI(1)
            HY    = 2.0_WP*PI*ALX2
        END IF
        dxplus = DX * REN * DSQRT(0.5_wp*CFGV)
        dzplus = DZ * REN * DSQRT(0.5_wp*CFGV)
        yMINplus = DYMIN* REN * DSQRT(0.5_wp*CFGV)
        yMAXplus = DYMAX* REN * DSQRT(0.5_wp*CFGV)
    
  
        WRITE(*,'(A)') '#**********Mesh and flow details based on initial flow field******************************'
        IF(TGFLOWflg) &
        WRITE(*,'(A,I5,3(2X,A,F9.5))') '#   NCL1_TG= ',NCL1_tg, 'HX_TG=', ALX1(1),'DX=  ', DX,   'DX+= ', dxplus
        IF(IOFLOWflg) &
        WRITE(*,'(A,I5,3(2X,A,F9.5))') '#   NCL1_IO= ',NCL1_io, 'HX_IO=', ALX1(2),'DX=  ', DX,   'DX+= ', dxplus
        WRITE(*,'(A,I5,3(2X,A,F9.5))') '#   NCL3=    ',NCL3,    'HZ=   ', ALX3,   'DZ=  ', DZ,   'DZ+= ', dzplus
        WRITE(*,'(A,I5,5(2X,A,F9.5))') '#   NCL2=    ',NCL2,    'HY=   ', HY,     'DY1= ', DYMIN,'Y1+= ', yMINplus
        WRITE(*,'(A,I5,5(2X,A,F9.5))') '#   NCL2=    ',NCL2,    'HY=   ', HY,     'DYc= ', DYMAX,'Yc+= ', yMAXplus  
                                                                        
        IF(phyTime .GT. TSTAV1) THEN                                                       
            dxplus = DX * REN * DSQRT(0.5_wp*CFGV)
            dzplus = DZ * REN * DSQRT(0.5_wp*CFGV)
            yMINplus = DYMIN* REN * DSQRT(0.5_wp*CFGV)
            yMAXplus = DYMAX* REN * DSQRT(0.5_wp*CFGV)                                                                                             
        END IF
                                                                                    
        IF(TGFLOWflg) &
        WRITE(*,'(A,I12)'            ) '#   MESH_SIZE_TG   =',NCL1_tg*NCL2*NCL3
        IF(IOFLOWflg) &
        WRITE(*,'(A,I12)'            ) '#   MESH_SIZE_IO   =',NCL1_io*NCL2*NCL3
        IF(IOFLOWflg .AND. TGFLOWFLG) &
        WRITE(*,'(A,I12)'            ) '#   MESH_SIZE_total=',(NCL1_io+NCL1_tg)*NCL2*NCL3
        WRITE(*,'(A)') '#****************************************************************************************'
        WRITE(*,'(A,2X,F18.5)')        '#   CONSTANT MEMORY PER CORE (Mb) = ', DBLE(MEMPC_byte)/1024.0_WP/1024.0_WP
        WRITE(*,'(A)') '#****************************************************************************************'
        
        !==========================INDICATORS===============================================================
        IF(IOFLOWflg .AND. TGFLOWFLG)  THEN
            CALL date_and_time(DATE=date,TIME=time)
            fllog=date(1:4)//'.'//date(5:8)//'.'//time(1:4)//'.log'
            logflg_tg = 110
            OPEN(logflg_tg,FILE='history.periodicxz.'//fllog)
            logflg_io = 6
        ELSE
            IF(TGFLOWFLG) logflg_tg = 6
            IF(IOFLOWFLG) logflg_io = 6
        END IF
        
        
        IF(TGFLOWFLG)  THEN
            WRITE(logflg_tg,'(A, 8X,A, 5X,A, 2X,A,3X,A, 7X,A, 2X,A, 4(9X,A,8X,A),2(7X,A,6X,A), 2(10X,A))') &
                        '#','1STEP','2PhyT', '3CpuT', '4CFL', '5DT','6MaxDiv', &
                        '07UU','08UUt','09UV','10UVt','11VV','12VVt','13WW','14WWT', &
                        '15Cf_L','16Cf_Lt','17Cf_U','18Cf_Ut', '19Umean', '20Umaxx'
        END IF
        
        IF(IOFLOWFLG) THEN
            IF(TGFLOWFLG) THEN
                IF(thermlflg==1) then
                    WRITE(logflg_io, '(13A)')&
                        '# "1STEP","2PhyT", "3CpuT", "4CFL", "5DT", "6DivMFD", "7DivINL", "8DivOUL", ', &
                        '"09UU","10UUt","11UV","12UVt", "13VV","14VVt","15WW","16WWT", ', &
                        '"17DD","18DDt","19TT","20TTt", "21HH", "22HHt", ', &
                        '"23Cf_L","24Cf_Lt","25Cf_U","26Cf_Ut", ', &
                        '"27UU","28UUt","29UV","30UVt", "31VV","32VVt","33WW","34WWT", ', &
                        '"35DD","36DDt","37TT","38TTt", "39HH", "40HHt", ', &
                        '"41Cf_L","42Cf_Lt","43Cf_U","44Cf_Ut", ' ,&
                        '"45UU","46UUt","47UV","48UVt", "49VV","50VVt","51WW","52WWT", ', &
                        '"53DD","54DDt","55TT","56TTt", "57HH", "58HHt", ', &
                        '"59Cf_L","60Cf_Lt","61Cf_U","62Cf_Ut", ' ,&
                        '"63Gmean", "64Gmaxx","65MassConvs" , "66EnegConvs", "67EnegTotal", ', &
                        '"68total_energy_RA", "69total_energy_FA", "70total_enstrophy", ' ,&
                        '"71UUV","72UUVt","73UVV","74UVVt","75Fcdrv"'
                else
                    WRITE(logflg_io,'(6A)') &
                        '# "1STEP","2PhyT", "3CpuT", "4CFL", "5DT", "6DivMFD", "7DivINL", "8DivOUL", ',&
                        '"09UU","10UUt","11UV","12UVt", "13VV","14VVt","15WW","16WWT","17Cf_L","18Cf_Lt","19Cf_U","20Cf_Ut", ',&
                        '"21UU","22UUt","23UV","24UVt", "25VV","26VVt","27WW","28WWT","29Cf_L","30Cf_Lt","31Cf_U","32Cf_Ut", ',&
                        '"33UU","34UUt","35UV","36UVt", "37VV","38VVt","39WW","40WWT","41Cf_L","42Cf_Lt","43Cf_U","44Cf_Ut", ',&
                        '"45Gmean", "46Gmaxx","47MassConvs", "48total_energy_RA", "49total_energy_FA", "50total_enstrophy", ',&
                        '"51UUV","52UUVt","53UVV","54UVVt","55Fcdrv"'
                end if
            ELSE
                IF(thermlflg==1) then
                    WRITE(logflg_io,'(A)')'#TITLE = " instantanous xz-periodic with thermal "'
                    WRITE(logflg_io,'(10A)') &
                                        '#VARIABLES = ', &
                                        '"1STEP", "2PhyT", "3CpuT", "4CFL", "5DT", "6MaxDiv", ', &
                                        '"07Uc", "08Uct", "09UUc", "10UUct", "11UVC", "12UVCt", ', &
                                        '"13Cf_B", "14Cf_Bt", "15Cf_T", "16Cf_Tt", ', &
                                        '"17Gmean", "18Gmaxx", "19MassConvs", ', &
                                        '"20Tc", "21Tct", "22Trms_B", "23Trmst_Bt","24Trms_T", "25Trmst_Tt",', &
                                        '"26Qwflux_B", "27Qwflux_Bt", "28Qwflux_T", "29Qwflux_Tt", ', &
                                        '"30Tbulk", "31Tmaxx", "32EnegConvs", "33EnegTotal", "34QwRatio", "35QwRatio_t", ', &
                                        '"36total_energy_RA", "37total_energy_FA", "38total_enstrophy" ', &
                                        '"39UUV","40UUVt","41UVV","42UVVt","43Fcdrv"'
                                        
                    WRITE(logflg_io,'(A)')'#ZONE T="tracking"' 
                ELSE
                    WRITE(logflg_io,'(A)')'#TITLE = " instantanous xz-periodic without thermal "'
                    WRITE(logflg_io,'(7A)')'#VARIABLES = ', & 
                                        '"1STEP", "2PhyT", "3CpuT", "4CFL", "5DT", "6MaxDiv", ', &
                                        '"07Uc", "08Uct", "09UUc", "10UUct", "11UVC", "12UVCt", ', &
                                        '"13Cf_B", "14Cf_Bt", "15Cf_T", "16Cf_Tt", ', &
                                        '"17Gmean", "18Gmaxx", "19MassConvs", ', &
                                        '"20total_energy_RA","21tke_energy_ave",',&
                                        '"22total_enstrophy","23enstrophy_ave","24UUV","25UUVt","26UVV","27UVVt","28Fcdrv"'
                    WRITE(logflg_io,'(A)')'#ZONE T="tracking"' 
                    
                    
                END IF
            END IF
            
        END IF
    
        RETURN
    END SUBROUTINE


!******************************************************************************************************************************
    SUBROUTINE PP_MONITOR_TG
        use flow_info
        use init_info
        use mesh_info
        use postprocess_info
        USE WRT_INFO
        IMPLICIT NONE 

        REAL(WP) :: CF_XZT(4)  !1=LW, 2=T_LW, 3=UW,  4=T_UW
        REAL(WP) :: CF_XZT_WORK(4)
        REAL(WP) :: dudy_xz, dudy_xzt
        REAL(WP) :: UU_XZT(8), U1_XZT(6)
        REAL(WP) :: UU_XZT_WORK(8), U1_XZT_WORK(6)
        integer(4) :: JID, J, JJ, JJSCN
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CPUTIME_TMP, CPUTIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

        !========Cf============================================
        CF_XZT(:) = 0.0_WP
        IF(ICASE.NE.IPIPEC) THEN
            IF(MYID.EQ.0) THEN
                JID = 1
                JJ=JCL2G(JID)
                dudy_xz = (U1xzL_tg(JID,1)-0.0_wp)*DYFI(JJ)*2.0_wp
                !dudy_xz = (Q_tg(NCL1_tg/2,JID,NCL3/2,1)-0.0_wp)*DYFI(JJ)*2.0_wp
                CF_XZT(1) = dudy_xz  * 2.0_WP / REN
            
                dudy_xzt = (U1xztL_tg(JID,1)-0.0_wp)*DYFI(JID)*2.0_wp
                CF_XZT(2) = dudy_xzt * 2.0_WP / REN
            END IF
        END IF
        
        IF(MYID.EQ.NPSLV) THEN
            JID = N2DO(MYID)
            JJ=JCL2G(JID)
            dudy_xz = (U1xzL_tg(JID,1)-0.0_wp)*DYFI(JJ)*2.0_wp
            !dudy_xz = (Q_tg(NCL1_tg/2,JID,NCL3/2,1)-0.0_wp)*DYFI(JJ)*2.0_wp
            CF_XZT(3) = dudy_xz  * 2.0_WP / REN
        
            dudy_xzt = (U1xztL_tg(JID,1)-0.0_wp)*DYFI(JJ)*2.0_wp
            CF_XZT(4) = dudy_xzt * 2.0_WP / REN
        END IF
        
        !=========UU,VV,WW,UV==================================
        IF(ICASE==IPIPEC) THEN
            JJSCN = 1
        ELSE
            JJSCN = NCL2/2
        END IF
        
        U1_XZT(:) = 0.0_WP
        U1_XZT_WORK(:) = 0.0_WP
        UU_XZT(:) = 0.0_WP
        UU_XZT_WORK(:) = 0.0_WP
        
        DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            IF(JJ.EQ.JJSCN) THEN
                U1_XZT(1) = U1xzL_tg(J,1)
                U1_XZT(3) = U1xzL_tg(J,2)
                U1_XZT(5) = U1xzL_tg(J,3)
                
                U1_XZT(2) = U1xztL_tg(J,1)
                U1_XZT(4) = U1xztL_tg(J,2)
                U1_XZT(6) = U1xztL_tg(J,3)
                
                UU_XZT(1) = U2xzL_tg(J,1)  - U1xzL_tg(J,1) *U1xzL_tg(J,1)!UUXZ
                UU_XZT(2) = U2xztL_tg(J,1) - U1xztL_tg(J,1)*U1xztL_tg(J,1)
                
                UU_XZT(3) = U2xzL_tg(J,2)  - U1xzL_tg(J,1) *U1xzL_tg(J,2)!UVXZ
                UU_XZT(4) = U2xztL_tg(J,2) - U1xztL_tg(J,1)*U1xztL_tg(J,2)
                
                UU_XZT(5) = U2xzL_tg(J,4)  - U1xzL_tg(J,2) *U1xzL_tg(J,2)!VVXZ
                UU_XZT(6) = U2xztL_tg(J,4) - U1xztL_tg(J,2)*U1xztL_tg(J,2)
                
                UU_XZT(7) = U2xzL_tg(J,6)  - U1xzL_tg(J,3) *U1xzL_tg(J,3)!WWXZ
                UU_XZT(8) = U2xztL_tg(J,6) - U1xztL_tg(J,3)*U1xztL_tg(J,3)
            END IF
        END DO
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CF_XZT, CF_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU_XZT, UU_XZT_WORK, 8, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(U1_XZT, U1_XZT_WORK, 6, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        
        !================PRINT DATA ON SCREEN FOR MONITORING==============
        IF (MYID.EQ.0) THEN
            WRITE(logflg_tg,'(1I14, 1F10.3, 2F7.3, 1F10.5, 1ES10.2,18ES13.5,2F17.9)') &
                    ITERG, phyTIME, CPUTIME, CFLMM*DT, DT, MAXDIVGV_TG(2),&
                    U1_XZT_WORK(1),U1_XZT_WORK(2), &
                    U1_XZT_WORK(3),U1_XZT_WORK(4), &
                    U1_XZT_WORK(5),U1_XZT_WORK(6), &
                    UU_XZT_WORK(1),UU_XZT_WORK(2), &
                    UU_XZT_WORK(3),UU_XZT_WORK(4), &
                    UU_XZT_WORK(5),UU_XZT_WORK(6), &
                    UU_XZT_WORK(7),UU_XZT_WORK(8), &
                    CF_XZT_WORK(1),CF_XZT_WORK(2), &
                    CF_XZT_WORK(3),CF_XZT_WORK(4), &
                    U1MEAN_WORK_TG, U1MAXX_WORK_TG
        ENDIF
         
    RETURN
    END SUBROUTINE
    
    
!******************************************************************************************************************************
!******************************************************************************************************************************
    SUBROUTINE PP_MONITOR_Xperiodic_IO
        use flow_info
        use init_info
        use mesh_info
        use postprocess_info
        use thermal_info
        USE WRT_INFO
        IMPLICIT NONE 
        !===top and bottom==================
        !REAL(WP) :: CFw_XZT(4), CFw_XZT_WORK(4)  !1=LW, 2=T_LW, 3=UW,  4=T_UW
        REAL(WP) :: Cfw_io(4)
        REAL(WP) :: Trmsw_XZT(4), Trmsw_XZT_WORK(4)
        REAL(WP) :: WHw(4), Hw, Tw, WH_ratio(2)=0.0
        
        !========centre======================
        REAL(WP) :: UGc_XZT(4), UGc_XZT_WORK(4)
        REAL(WP) :: UGU_XZT(4),UGU_XZT_WORK(4)
        REAL(WP) :: Uc_XZT(2), Uc_XZT_WORK(2)
        REAL(WP) :: Tc_XZT(2), Tc_XZT_WORK(2)
        INTEGER(4) :: I, K, IP, KP, JP, IM, KM, JM  
        integer(4) :: JID, J, JJ, JJSCN, N, JJA, JJB, N1,N2
        INTEGER(4) :: H, M, LMNH
        REAL(wp)   :: tmp, tmp1, tmp2
        REAL(WP)   :: total_energy_FA, total_energy_FA_work
        REAL(WP)   :: total_energy_RA, total_energy_RA_work
        REAL(WP)   :: total_energy_M2, total_energy_M2_work
        REAL(WP)   :: total_enstrophy, total_enstrophy_work
        REAL(WP)   :: total_ensphy_M2, total_ensphy_M2_work
        
        !REAL(WP)   :: aveg_energy_RA, aveg_energy_RA_work
        !REAL(WP)   :: aveg_enstrophy, aveg_enstrophy_work
        

        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CPUTIME_TMP, CPUTIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        
        CALL CHK_MassConsv_io
        IF(thermlflg==1) CALL CHK_EnegConsv_io
        
        Cfw_io(:)        = 0.0_wP
        WHw(:)           = 0.0_wP
        WH_ratio(:)      = 0.0_WP
        !============================info on the wall=========================================
        CALL PP_wall_thermal_shear(flgxz)
        Cfw_io(1)  = 2.0_WP*Tauw_io(1)
        Cfw_io(3)  = 2.0_WP*Tauw_io(2) 
        
        IF(thermlflg==1) THEN
            IF( BCWALLHEAT(itopwall)==isoFluxWall .and. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
                WHw(1) = Twal(1)*T0 ! Tw_d
                WHw(3) = Twal(2)*T0 ! Tw_d
            END IF
            IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                WHw(1) = qw(1)*D0*U0*CP0*T0
                WHw(3) = qw(2)*D0*U0*CP0*T0
            END IF
        END IF
        
        IF(phyTime .GT. TSTAV1) THEN 
            CALL PP_wall_thermal_shear(flgxzt)
            Cfw_io(2)  = 2.0_WP*Tauw_io(1)
            Cfw_io(4)  = 2.0_WP*Tauw_io(2) 
            IF(thermlflg==1) THEN
                IF( BCWALLHEAT(itopwall)==isoFluxWall .and. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
                    WHw(2) = Twal(1)*T0 ! Tw_d
                    WHw(4) = Twal(2)*T0 ! Tw_d
                END IF
                IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                    WHw(2) = qw(1)*D0*U0*CP0*T0
                    WHw(4) = qw(2)*D0*U0*CP0*T0
                END IF
            END IF
        ELSE
            Cfw_io(2)  = Cfw_io(1)
            Cfw_io(4)  = Cfw_io(3)
            IF(thermlflg==1) THEN
                WHw(2) = WHw(1)
                WHw(4) = WHw(3)
            ENDIF
        END IF  
        
        !================================================================================
        Trmsw_XZT(:)     = 0.0_WP 
        Trmsw_XZT_WORK(:)= 0.0_WP 
        IF( (thermlflg==1) .and. ( (ICASE.NE.IPIPEC .and. MYID.EQ.0) .or. MYID.EQ.NPSLV ) ) THEN
            IF(ICASE.NE.IPIPEC .and. MYID.EQ.0) THEN
                JID = 1
                JJ  = JCL2G(JID)
                N1  = 1
                N2  = 2
            END IF
            IF(MYID.EQ.NPSLV) THEN
                JID = N2DO(MYID)
                JJ  = JCL2G(JID)
                N1  = 3
                N2  = 4
            END IF
            
            !==================Trms======================
            DO N=N1, N2
                tmp2 = T2xzL_io(JID)
                tmp1 = T1xzL_io(JID)
                IF(N==N2) THEN
                    IF(phyTime .GT. TSTAV1) THEN
                        tmp2 = T2xztL_io(JID)
                        tmp1 = T1xztL_io(JID)
                        
                    ELSE
                        Trmsw_XZT(N) = Trmsw_XZT(N1)
                        EXIT
                    END IF
                END IF
                Trmsw_XZT(N) = DSQRT(DABS(tmp2-tmp1*tmp1))
            END DO

        END IF
            
        !============================info on y-centre =========================================
        IF(ICASE==IPIPEC) THEN
            JJSCN = 1
        ELSE
            JJSCN = NCL2/2
        END IF
        
        UGc_XZT(:)    = 0.0_WP
        UGU_XZT(:)    = 0.0_WP
        Uc_XZT(:)     = 0.0_WP
        Tc_XZT(:)     = 0.0_WP
        
        UGc_XZT_WORK(:)    = 0.0_WP
        UGU_XZT_WORK(:)    = 0.0_WP
        Uc_XZT_WORK(:)     = 0.0_WP
        Tc_XZT_WORK(:)     = 0.0_WP
        DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            IF(JJ.EQ.JJSCN) THEN
                UGc_XZT(1) = UGxzL_io (J,1) - G1xzL_io (J,1) *G1xzL_io (J,1)/D1xzL_io (J)!UUXZ
                UGc_XZT(3) = UGxzL_io (J,2) - G1xzL_io (J,1) *G1xzL_io (J,2)/D1xzL_io (J)!UVXZ
                Uc_XZT(1)  = G1xzL_io (J,1)/D1xzL_io(J)
                
                M = 1
                N = 1
                H = 2
                LMNH = M*(6-M)+(N*(7-N))/2+H-8
                UGU_XZT(1) = UGUxzL_io(J,LMNH)/D1xzL_io(J) &
                            -G1xzL_io (J,M)/D1xzL_io(J) * UGxzL_io(J,(N*(7-N))/2+H-3)/D1xzL_io(J) &
                            -G1xzL_io (J,N)/D1xzL_io(J) * UGxzL_io(J,(M*(7-M))/2+H-3)/D1xzL_io(J) &
                            -G1xzL_io (J,H)/D1xzL_io(J) * UGxzL_io(J,(M*(7-M))/2+N-3)/D1xzL_io(J) &
                            +2.0_WP * G1xzL_io(J,M)/D1xzL_io(J) * &
                                      G1xzL_io(J,N)/D1xzL_io(J) * &
                                      G1xzL_io(J,H)/D1xzL_io(J)
                                      
                M = 1
                N = 2
                H = 2
                LMNH = M*(6-M)+(N*(7-N))/2+H-8
                UGU_XZT(3) = UGUxzL_io(J,LMNH)/D1xzL_io(J) &
                            -G1xzL_io (J,M)/D1xzL_io(J) * UGxzL_io(J,(N*(7-N))/2+H-3)/D1xzL_io(J) &
                            -G1xzL_io (J,N)/D1xzL_io(J) * UGxzL_io(J,(M*(7-M))/2+H-3)/D1xzL_io(J) &
                            -G1xzL_io (J,H)/D1xzL_io(J) * UGxzL_io(J,(M*(7-M))/2+N-3)/D1xzL_io(J) &
                            +2.0_WP * G1xzL_io(J,M)/D1xzL_io(J) * &
                                      G1xzL_io(J,N)/D1xzL_io(J) * &
                                      G1xzL_io(J,H)/D1xzL_io(J)
                                          
                IF(thermlflg==1) Tc_XZT(1) = T1xzL_io(J)
            END IF
        END DO
        
        IF(phyTime .GT. TSTAV1) THEN
            DO J=1,N2DO(MYID)
                JJ=JCL2G(J)
                IF(JJ.EQ.JJSCN) THEN
                    UGc_XZT(2) = UGxztL_io(J,1) - G1xztL_io(J,1) *G1xztL_io(J,1)/D1xztL_io(J)
                    UGc_XZT(4) = UGxztL_io(J,2) - G1xztL_io(J,1) *G1xztL_io(J,2)/D1xztL_io(J)
                    Uc_XZT(2)  = G1xztL_io (J,1)/D1xztL_io(J)
                    
                    M = 1
                    N = 1
                    H = 2
                    LMNH = M*(6-M)+(N*(7-N))/2+H-8
                    UGU_XZT(2) = UGUxztL_io(J,LMNH)/D1xztL_io(J) &
                                -G1xztL_io (J,M)/D1xztL_io(J) * UGxztL_io(J,(N*(7-N))/2+H-3)/D1xztL_io(J) &
                                -G1xztL_io (J,N)/D1xztL_io(J) * UGxztL_io(J,(M*(7-M))/2+H-3)/D1xztL_io(J) &
                                -G1xztL_io (J,H)/D1xztL_io(J) * UGxztL_io(J,(M*(7-M))/2+N-3)/D1xztL_io(J) &
                                +2.0_WP * G1xztL_io(J,M)/D1xztL_io(J) * &
                                          G1xztL_io(J,N)/D1xztL_io(J) * &
                                          G1xztL_io(J,H)/D1xztL_io(J)
                                          
                    M = 1
                    N = 2
                    H = 2
                    LMNH = M*(6-M)+(N*(7-N))/2+H-8
                    UGU_XZT(4) = UGUxztL_io(J,LMNH)/D1xztL_io(J) &
                                -G1xztL_io (J,M)/D1xztL_io(J) * UGxztL_io(J,(N*(7-N))/2+H-3)/D1xztL_io(J) &
                                -G1xztL_io (J,N)/D1xztL_io(J) * UGxztL_io(J,(M*(7-M))/2+H-3)/D1xztL_io(J) &
                                -G1xztL_io (J,H)/D1xztL_io(J) * UGxztL_io(J,(M*(7-M))/2+N-3)/D1xztL_io(J) &
                                +2.0_WP * G1xztL_io(J,M)/D1xztL_io(J) * &
                                          G1xztL_io(J,N)/D1xztL_io(J) * &
                                          G1xztL_io(J,H)/D1xztL_io(J)
                                      
                    IF(thermlflg==1) Tc_XZT(2) = T1xztL_io(J)
                END IF
            END DO
        ELSE
            UGc_XZT(2)=UGc_XZT(1)
            UGc_XZT(4)=UGc_XZT(3)
            
            UGU_XZT(2)=UGU_XZT(1)
            UGU_XZT(4)=UGU_XZT(3)
            
            Uc_XZT(2) =Uc_XZT(1)
            Tc_XZT(2) =Tc_XZT(1) 
        END IF
        
        !====total energy and enstrophy==============
        total_energy_RA = 0.0_wp
        total_energy_FA = 0.0_wp
        total_enstrophy = 0.0_wp
        DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            !total_energy_RA =   total_energy_RA + ( &
            !                    U2xzL_io (J,1) - U1xzL_io (J,1) * U1xzL_io (J,1) + &
            !                    U2xzL_io (J,2) - U1xzL_io (J,2) * U1xzL_io (J,2) + &
            !                    U2xzL_io (J,3) - U1xzL_io (J,3) * U1xzL_io (J,3) ) * 0.5_wp * (1.0W_P/DYFI(JJ))
                                
            total_energy_RA =   total_energy_RA + ( &
                                U2xzL_io (J,1)  + &
                                U2xzL_io (J,2)  + &
                                U2xzL_io (J,3)  ) * 0.5_wp * (1.0_WP/DYFI(JJ))
                                
            total_energy_FA =   total_energy_FA + ( &
                                UGxzL_io (J,1) - G1xzL_io (J,1) *G1xzL_io (J,1)/D1xzL_io (J) + &
                                UGxzL_io (J,2) - G1xzL_io (J,2) *G1xzL_io (J,2)/D1xzL_io (J) + &
                                UGxzL_io (J,3) - G1xzL_io (J,3) *G1xzL_io (J,3)/D1xzL_io (J) ) * 0.5_wp* (1.0_WP/DYFI(JJ))
                                
                
                                
            total_enstrophy =   total_enstrophy + ( &
                               (DVDL2xzL_io(J,(3-1)*NDV+2,(3-1)*NDV+2) - 2.0_wp* &
                                DVDL2xzL_io(J,(3-1)*NDV+2,(2-1)*NDV+3) +         &
                                DVDL2xzL_io(J,(2-1)*NDV+3,(2-1)*NDV+3) ) -       &
                               (DVDL1xzL_io(J,3,2)*DVDL1xzL_io(J,3,2) - 2.0_wp * &
                                DVDL1xzL_io(J,3,2)*DVDL1xzL_io(J,2,3)  +         &
                                DVDL1xzL_io(J,2,3)*DVDL1xzL_io(J,2,3) ) +        &
                               (DVDL2xzL_io(J,(1-1)*NDV+3,(1-1)*NDV+3) - 2.0_wp* &
                                DVDL2xzL_io(J,(1-1)*NDV+3,(3-1)*NDV+1) +         &
                                DVDL2xzL_io(J,(3-1)*NDV+1,(3-1)*NDV+1) ) -       &
                               (DVDL1xzL_io(J,1,3)*DVDL1xzL_io(J,1,3) - 2.0_wp * &
                                DVDL1xzL_io(J,1,3)*DVDL1xzL_io(J,3,1) +          &
                                DVDL1xzL_io(J,3,1)*DVDL1xzL_io(J,3,1) ) +        &
                               (DVDL2xzL_io(J,(2-1)*NDV+1,(2-1)*NDV+1) - 2.0_wp* &
                                DVDL2xzL_io(J,(2-1)*NDV+1,(1-1)*NDV+2) +         &
                                DVDL2xzL_io(J,(1-1)*NDV+2,(1-1)*NDV+2) ) -       &
                               (DVDL1xzL_io(J,2,1)*DVDL1xzL_io(J,2,1) - 2.0_wp * &
                                DVDL1xzL_io(J,2,1)*DVDL1xzL_io(J,1,2) +          &
                                DVDL1xzL_io(J,1,2)*DVDL1xzL_io(J,1,2) ) ) * 0.5_wp* (1.0_WP/DYFI(JJ))
        END DO
        
        
        ! Below is for TGV
        total_energy_M2 = 0.0_wp
        DO I=1, NCL1_io
            IP = IPV_io(I)
            DO K=1, NCL3
                KP = KPV(K)
                DO J=1, N2DO(MYID)
                    JP = JLPV(J)
                    JJ = JCL2G(J)
                    total_energy_M2 = total_energy_M2 + ( &
                    ((Q_io(I,J,K,1) + Q_io(IP,J,K,1))*0.5_WP)**2 + &
                    ((Q_io(I,J,K,2) + Q_io(I,JP,K,2))*0.5_WP)**2 + &
                    ((Q_io(I,J,K,3) + Q_io(I,J,KP,3))*0.5_WP)**2 )*0.5_WP*DX*DZ/DYFI(JJ)
                END DO
            END DO
        END DO
        
        
        total_ensphy_M2 = 0.0_WP
        DO I=1, NCL1_io
            IP = IPV_io(I)
            IM = IMV_io(I)
            DO K=1, NCL3
                KP = KPV(K)
                KM = KMV(K)
                DO J=1, N2DO(MYID)
                    JP = JLPV(J)
                    JM = JLMV(J)
                    JJ = JCL2G(J)
                    
                    total_ensphy_M2 = total_ensphy_M2 + ( &
                    ((( ( (Q_io(I,JP,K,3) + Q_io(I,JP,KP,3))*0.5_WP ) + &  
                        ( (Q_io(I,J, K,3) + Q_io(I,J, KP,3))*0.5_WP ) )*0.5_WP - &
                      ( ( (Q_io(I,J, K,3) + Q_io(I,J, KP,3))*0.5_WP ) + &  
                        ( (Q_io(I,JM,K,3) + Q_io(I,JM,KP,3))*0.5_WP ) )*0.5_WP)*DYFI(JJ) - &
                     (( ( (Q_io(I,J,KP,2) + Q_io(I,JP,KP,2))*0.5_WP ) + &  
                        ( (Q_io(I,J, K,2) + Q_io(I,JP, K,2))*0.5_WP ) )*0.5_WP - &
                      ( ( (Q_io(I,J, K,2) + Q_io(I,JP, K,2))*0.5_WP ) + &  
                        ( (Q_io(I,J,KM,2) + Q_io(I,JP,KM,2))*0.5_WP ) )*0.5_WP)*DZI )**2 +  &    
                    ((( ( (Q_io(I,J,KP,1) + Q_io(IP,J,KP,1))*0.5_WP ) + &  
                        ( (Q_io(I,J, K,1) + Q_io(IP,J, K,1))*0.5_WP ) )*0.5_WP - &
                      ( ( (Q_io(I,J, K,1) + Q_io(IP,J, K,1))*0.5_WP ) + &  
                        ( (Q_io(I,J,KM,1) + Q_io(IP,J,KM,1))*0.5_WP ) )*0.5_WP)*DZI - &
                     (( ( (Q_io(IP,J,K,3) + Q_io(IP,J,KP,3))*0.5_WP ) + &  
                        ( (Q_io(I,J, K,3) + Q_io(I,J, KP,3))*0.5_WP ) )*0.5_WP - &
                      ( ( (Q_io(I,J, K,3) + Q_io(I,J, KP,3))*0.5_WP ) + &  
                        ( (Q_io(IM,J,K,3) + Q_io(IM,J,KP,3))*0.5_WP ) )*0.5_WP)*DXI )**2 +  &     
                    ((( ( (Q_io(IP,J,K,2) + Q_io(IP,JP,K,2))*0.5_WP ) + &  
                        ( (Q_io(I,J, K,2) + Q_io(I,JP, K,2))*0.5_WP ) )*0.5_WP - &
                      ( ( (Q_io(I,J, K,2) + Q_io(I,JP, K,2))*0.5_WP ) + &  
                        ( (Q_io(IM,J,K,2) + Q_io(IM,JP,K,2))*0.5_WP ) )*0.5_WP)*DXI - &
                     (( ( (Q_io(I,JP,K,1) + Q_io(IP,JP,K,1))*0.5_WP ) + &  
                        ( (Q_io(I,J, K,1) + Q_io(IP,J, K,1))*0.5_WP ) )*0.5_WP - &
                      ( ( (Q_io(I,J, K,1) + Q_io(IP,J, K,1))*0.5_WP ) + &  
                        ( (Q_io(I,JM,K,1) + Q_io(IP,JM,K,1))*0.5_WP ) )*0.5_WP)*DYFI(JJ) )**2 )*DX*DZ/DYFI(JJ)
                END DO
            END DO
        END DO
        
        !============================info on y-centre =========================================
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(UGc_XZT, UGc_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UGU_XZT, UGU_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(Uc_XZT,  Uc_XZT_WORK,  2, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(total_energy_RA,  total_energy_RA_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(total_energy_FA,  total_energy_FA_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(total_enstrophy,  total_enstrophy_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(total_energy_M2,  total_energy_M2_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(total_ensphy_M2,  total_ensphy_M2_WORK,  1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        total_energy_RA_WORK  = total_energy_RA_WORK/((HYT-HYB))
        total_energy_FA_WORK  = total_energy_FA_WORK/((HYT-HYB))
        total_enstrophy_WORK  = total_enstrophy_WORK/((HYT-HYB))
        
        total_energy_M2_WORK  = total_energy_M2_WORK/(HX_io*HZ*(HYT-HYB))
        total_ensphy_M2_WORK  = total_ensphy_M2_WORK/(HX_io*HZ*(HYT-HYB))
        
        !IF(MYID==0) THEN
        !    WRITE(*,*) '##',phyTIME, total_energy_RA_WORK, total_energy_M2_WORK, total_enstrophy_WORK, total_ensphy_M2_WORK
       ! END IF
        IF(thermlflg==1) THEN
            CALL MPI_ALLREDUCE(Trmsw_XZT, Trmsw_XZT_WORK, 4, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
            CALL MPI_ALLREDUCE(Tc_XZT,    Tc_XZT_WORK,    2, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        END IF
        
        !=====================ON THE TOP WALL=================================
        IF(MYID.EQ.NPSLV) THEN
            IF(thermlflg==1) THEN
                IF(ICASE.NE.IPIPEC) THEN
                    WH_ratio(1) =  ( DABS(WHw(1)) - DABS(WHw(3)) )/ ( (DABS(WHw(1)) + DABS(WHw(3))) )
                    WH_ratio(2) =  ( DABS(WHw(2)) - DABS(WHw(4)) )/ ( (DABS(WHw(2)) + DABS(WHw(4))) )
                END IF
                    
                WRITE(logflg_io,'(1I14, 1F12.5, 3F9.5, 1ES11.3, 10ES13.5, 2F9.5, 1ES13.5, 10ES13.5, 2F9.5, 12ES13.5)') &
                    ITERG, phyTIME, CPUTIME, CFLMM*DT, DT, MAXDIVGV_IO(1), &
                    Uc_XZT_WORK(1:2), UGc_XZT_WORK(1:4),   Cfw_io(1:4), G1BULK_WORK_IO, G1MAXX_WORK_IO, CHK_MASS_CONSV0, &
                    Tc_XZT_WORK(1:2), Trmsw_XZT_WORK(1:4), WHw(1:4),    T1BULK_WORK_IO, T1MAXX_WORK_IO, CHK_ENEG_CONSV0, &
                    CHK_ENEG_TOTAL, WH_ratio(2), total_energy_RA_WORK, total_energy_FA_WORK, total_enstrophy_WORK, &
                    UGU_XZT_WORK(1:4), FcDrv_IO
                            
            ELSE
                WRITE(logflg_io,'(1I14, 1F12.5, 3F9.5, 1ES11.3, 10ES13.5, 2F9.5, 10ES13.5)') &
                    ITERG, phyTIME, CPUTIME, CFLMM*DT, DT, MAXDIVGV_IO(1), &
                    Uc_XZT_WORK(1:2), UGc_XZT_WORK(1:4),  Cfw_io(1:4), G1BULK_WORK_IO, G1MAXX_WORK_IO, CHK_MASS_CONSV0, &
                    total_energy_RA_WORK, total_energy_M2_WORK, total_enstrophy_WORK, total_ensphy_M2_WORK, &
                    UGU_XZT_WORK(1:4),FcDrv_IO
            
            END IF
        END IF
         
    RETURN
    END SUBROUTINE
      
!******************************************************************************************************************    
    SUBROUTINE PP_MONITOR_nonXperiodic_IO
!>  @note: three point at x direction.
        use flow_info
        use init_info
        use mesh_info
        use postprocess_info
        USE WRT_INFO
        use thermal_info
        IMPLICIT NONE 
        
        INTEGER(4) :: IID(NDV) ! NDV=3 ONLY SHOW NUMBERS, NOT MEANING.
        REAL(WP) :: CF_ZT(4,NDV)  !1=LW, 2=T_LW, 3=UW,  4=T_UW
        REAL(WP) :: CF_ZT_WORK(4,NDV)
        REAL(WP) :: UG_ZT(8,NDV)
        REAL(WP) :: UG_ZT_WORK(8,NDV)
        
        REAL(WP) :: D2_ZT(2,NDV)
        REAL(WP) :: H2_ZT(2,NDV)
        REAL(WP) :: T2_ZT(2,NDV)
        
        REAL(WP) :: dudy_z, dudy_zt
        INTEGER(4) :: N, JID, J, JJ, JJSCN
        
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CPUTIME_TMP, CPUTIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        
        !CALL CHK_MassConsv_io
        IF(thermlflg==1) CALL CHK_EnegConsv_io
        
        
        DO N=1,NDV
           IID(N) = NCL1_IO/(NDV+1)*N
        END DO
        
        !========Cf============================================
        CF_ZT(:,:) = 0.0_wp
        
        IF(ICASE.NE.IPIPEC) THEN
            IF(MYID.EQ.0) THEN
                JID = 1
                JJ=JCL2G(JID)
                DO N = 1, NDV
                    dudy_z = (U1zL_io(IID(N),JID,1)-0.0_wp)*DYFI(JID)*2.0_wp
                    !dudy_z = (Q_io(IID(N),JID,NCL3/2,1)-0.0_wp)*DYFI(JJ)*2.0_wp
                    CF_ZT(1,N) = dudy_z  * 2.0_WP / REN
                    
                    dudy_zt = (U1ztL_io(IID(N),JID,1)-0.0_wp)*DYFI(JJ)*2.0_wp
                    CF_ZT(2,N) = dudy_zt * 2.0_WP / REN
                END DO
            END IF
        END IF
        
        IF(MYID.EQ.NPSLV) THEN
            JID = N2DO(MYID)
            JJ=JCL2G(JID)
            DO N = 1, NDV
                dudy_z = (U1zL_io(IID(N),JID,1)-0.0_wp)*DYFI(jj)*2.0_wp
                !dudy_z = (Q_io(IID(N),JID,NCL3/2,1)-0.0_wp)*DYFI(JJ)*2.0_wp
                CF_ZT(3,N) = dudy_z  * 2.0_WP / REN
                dudy_zt = (U1ztL_io(IID(N),JID,1)-0.0_wp)*DYFI(JJ)*2.0_wp
                CF_ZT(4,N) = dudy_zt * 2.0_WP / REN
            END DO
        END IF
        
        !=========UU,VV,WW,UV==================================
        IF(ICASE==IPIPEC) THEN
            JJSCN = 1
        ELSE
            JJSCN = NCL2/2
        END IF
        
        UG_ZT(:,:) = 0.0_WP
        DO J=1,N2DO(MYID)
            JJ=JCL2G(J)
            IF(JJ.EQ.JJSCN) THEN
                DO N=1,NDV
                    UG_ZT(1,N) = UGzL_IO (IID(N),J,1) - G1zL_IO (IID(N),J,1) * G1zL_IO (IID(N),J,1)/D1zL_io (IID(N),J)!UUXZ
                    UG_ZT(3,N) = UGzL_IO (IID(N),J,2) - G1zL_IO (IID(N),J,1) * G1zL_IO (IID(N),J,2)/D1zL_io (IID(N),J)!UVXZ
                    UG_ZT(5,N) = UGzL_IO (IID(N),J,4) - G1zL_IO (IID(N),J,2) * G1zL_IO (IID(N),J,2)/D1zL_io (IID(N),J)!VVXZ
                    UG_ZT(7,N) = UGzL_IO (IID(N),J,6) - G1zL_IO (IID(N),J,3) * G1zL_IO (IID(N),J,3)/D1zL_io (IID(N),J)!WWXZ
                END DO
            END IF
        END DO
        
        IF(phyTime .GT. TSTAV1) THEN
            DO J=1,N2DO(MYID)
                JJ=JCL2G(J)
                IF(JJ.EQ.JJSCN) THEN
                    DO N=1,NDV
                        UG_ZT(2,N) = UGztL_IO(IID(N),J,1) - G1ztL_IO(IID(N),J,1) * G1ztL_IO(IID(N),J,1)/D1ztL_io(IID(N),J)
                        UG_ZT(4,N) = UGztL_IO(IID(N),J,2) - G1ztL_IO(IID(N),J,1) * G1ztL_IO(IID(N),J,2)/D1ztL_io(IID(N),J)
                        UG_ZT(6,N) = UGztL_IO(IID(N),J,4) - G1ztL_IO(IID(N),J,2) * G1ztL_IO(IID(N),J,2)/D1ztL_io(IID(N),J)
                        UG_ZT(8,N) = UGztL_IO(IID(N),J,6) - G1ztL_IO(IID(N),J,3) * G1ztL_IO(IID(N),J,3)/D1ztL_io(IID(N),J)
                    END DO
                END IF
            END DO
        END IF
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CF_ZT, CF_ZT_WORK, 4*NDV, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UG_ZT, UG_ZT_WORK, 8*NDV, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        !=====================ON THE TOP WALL=================================
        IF(MYID.EQ.NPSLV) THEN
        
        
            IF(thermlflg==1) THEN
        
                J = N2DO(MYID)
                DO N = 1, NDV
                    D2_ZT(1,N) = D2zL_io(IID(N),J) - D1zL_io(IID(N),J) * D1zL_io(IID(N),J)
                    T2_ZT(1,N) = T2zL_io(IID(N),J) - T1zL_io(IID(N),J) * T1zL_io(IID(N),J)
                    H2_ZT(1,N) = H2zL_io(IID(N),J) - H1zL_io(IID(N),J) * H1zL_io(IID(N),J)
                    
                    D2_ZT(2,N) = D2zTL_io(IID(N),J) - D1zTL_io(IID(N),J) * D1zTL_io(IID(N),J)
                    T2_ZT(2,N) = T2zTL_io(IID(N),J) - T1zTL_io(IID(N),J) * T1zTL_io(IID(N),J)
                    H2_ZT(2,N) = H2zTL_io(IID(N),J) - H1zTL_io(IID(N),J) * H1zTL_io(IID(N),J)
                    
                END DO
                WRITE(logflg_io,'(1I14, 1F10.3, 2F7.3, 1F10.5, 3ES10.2, 54ES13.5, 2F13.5, 3ES13.5)') &
                            ITERG, phyTIME, CPUTIME, CFLMM*DT,  DT, MAXDIVGV_IO(1:3), &
                          ( UG_ZT_WORK(1,N),UG_ZT_WORK(2,N),UG_ZT_WORK(3,N),UG_ZT_WORK(4,N), &
                            UG_ZT_WORK(5,N),UG_ZT_WORK(6,N),UG_ZT_WORK(7,N),UG_ZT_WORK(8,N), &
                            CF_ZT_WORK(1,N),CF_ZT_WORK(2,N),CF_ZT_WORK(3,N),CF_ZT_WORK(4,N), &
                            D2_ZT(1,N), D2_ZT(2,N), T2_ZT(1,N), T2_ZT(2,N), H2_ZT(1,N), H2_ZT(2,N), &
                            N=1,NDV), G1BULK_WORK_IO, G1MAXX_WORK_IO, CHK_MASS_CONSV0, CHK_ENEG_CONSV0, &
                            CHK_ENEG_TOTAL
                            
            ELSE
                WRITE(logflg_io,'(1I14, 1F10.3, 2F7.3, 1F10.5, 3ES10.2, 36ES13.5, 1F13.5, 2ES13.5)') &
                            ITERG, phyTIME, CPUTIME, CFLMM*DT,  DT, MAXDIVGV_IO(1:3), &
                          ( UG_ZT_WORK(1,N),UG_ZT_WORK(2,N),UG_ZT_WORK(3,N),UG_ZT_WORK(4,N), &
                            UG_ZT_WORK(5,N),UG_ZT_WORK(6,N),UG_ZT_WORK(7,N),UG_ZT_WORK(8,N), &
                            CF_ZT_WORK(1,N),CF_ZT_WORK(2,N),CF_ZT_WORK(3,N),CF_ZT_WORK(4,N), &
                            N=1,NDV), G1BULK_WORK_IO, G1MAXX_WORK_IO, CHK_MASS_CONSV0
            
            END IF
        END IF
         
        RETURN
    END SUBROUTINE
!****************************************************************************************************************** 

    SUBROUTINE PP_wall_thermal_shear(flg_xzt)
        use flow_info
        use init_info
        use mesh_info
        use postprocess_info
        use thermal_info
        USE WRT_INFO
        IMPLICIT NONE 
        
        INTEGER(4),INTENT(IN) :: FLG_XZT
        INTEGER(4) :: JID, JJ, JJA, N1, NSZ, J
        REAL(WP) :: DUMMY(7,2)
        REAL(WP) :: DUMMY_WORK(7,2)
        REAL(WP) :: DUDY
        REAL(WP) :: DEN, YL, VIS
        REAL(WP) :: Tau_diff, Tau_avag
        

        IF(thermlflg==1) CALL PP_Wall_thermal_properties(flg_xzt)

        !================arithmetic mean of density==================
        DEN = 0.0_WP
        VIS = 0.0_WP
        YL  = 0.0_WP
        DO J=1, N2DO(MYID)
            JJ = JCL2G(J)
            YL = YL + 1.0_WP/DYFI(JJ)
            IF(flg_xzt==flgxz) then
                DEN = DEN + D1xzL_io (J)/DYFI(JJ)
                VIS = VIS + M1xzL_io (J)/DYFI(JJ)
            ELSE IF (flg_xzt==flgxzt)  then
                DEN = DEN + D1xztL_io(J)/DYFI(JJ)
                VIS = VIS + M1xztL_io(J)/DYFI(JJ)
            ELSE
            END IF
        END DO
        
        !================arithmetic mean of density==================
        Utaw_io(1:2) = 0.0_wp
        Tauw_io(1:2) = 0.0_wp
        Ret_io(1:2)  = 0.0_wp
        DUDY = 0.0_WP
        IF((ICASE.NE.IPIPEC .and. MYID.EQ.0) .or. MYID.EQ.NPSLV ) THEN ! FOR TOP OR BOTTOM WALL RANKID
        
            
            IF(ICASE.NE.IPIPEC .and. MYID.EQ.0) THEN
                JID = 1
                JJ  = JCL2G(JID) !1
                JJA = 1
                N1  = 1
            END IF
            IF(MYID.EQ.NPSLV) THEN
                JID = N2DO(MYID) !LOCAL
                JJ  = JCL2G(JID) !GLOCAL
                JJA = NND2       !GLOBAL
                N1  = 2
            END IF
            
            !==================du/dy======================  
            IF(flg_xzt==flgxz) then
                DuDy =  ( U1xzL_io (JID,1)-0.0_wp)/(YCC(JJ)-YND(JJA))
            else if (flg_xzt==flgxzt)  then
                DuDy =  ( U1xztL_io(JID,1)-0.0_wp)/(YCC(JJ)-YND(JJA))
            else
            end if    
            
            !=================skin varaibles==================
            IF(thermlflg ==0) THEN
            
                Tauw_io(N1) = DUDy/REN               ! undim tau_w = tau_w/(\rho*u2)
                Utaw_io(N1) = DSQRT(dabs(Tauw_io(N1)))
                Ret_io(N1)  =  REN * Utaw_io(N1)
                
            ELSE IF (thermlflg ==1) THEN
            
                Tauw_io(N1) = DUDy * Mwal(N1) /REN
                Utaw_io(N1) = DSQRT( dabs(Tauw_io(N1)/Dwal(N1)) ) 
                Ret_io(N1)  = REN * Utaw_io(N1)*Dwal(N1)/Mwal(N1)
                
                Tauw_d_io(N1) = Tauw_io(N1) * U0 * U0 * D0 !dimensional
                Utaw_d_io(N1) = Utaw_io(N1) * U0           !dimensional
                
            ELSE
            END IF
            
        END IF

        DUMMY(1,1:2) =  Utaw_io  (1:2)
        DUMMY(2,1:2) =  Tauw_io  (1:2)
        DUMMY(3,1:2) =  Ret_io   (1:2)
        DUMMY(4,1:2) =  Utaw_d_io(1:2)
        DUMMY(5,1:2) =  Tauw_d_io(1:2)
        DUMMY(6,1)   =  YL
        DUMMY(6,2)   =  DEN 
        DUMMY(7,1)   =  VIS
        DUMMY(7,2)   =  0.0_WP
        NSZ = 7*2
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(DUMMY(1,1), DUMMY_WORK(1,1), NSZ, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        Utaw_io  (1:2) = DUMMY_WORK(1,1:2)
        Tauw_io  (1:2) = DUMMY_WORK(2,1:2)
        Ret_io   (1:2) = DUMMY_WORK(3,1:2)
        Utaw_d_io(1:2) = DUMMY_WORK(4,1:2) 
        Tauw_d_io(1:2) = DUMMY_WORK(5,1:2)
        
        YL = DUMMY_WORK(6,1)
        DEN= DUMMY_WORK(6,2)
        VIS= DUMMY_WORK(7,1)
        
        Tauw_io  (1:2) = DABS(Tauw_io  (1:2))
        
        IF(thermlflg ==1) THEN
            Tau_diff = Tauw_io(2)-Tauw_io(1)!;  write(*,*) 'Tau_diff',Tau_diff , Tauw_io(2),Tauw_io(1)
            Tau_avag = Tauw_io(2)+Tauw_io(1)!;  write(*,*) 'Tau_avag',Tau_avag
            Ldist_io(1) = DABS(-1.0_wp+dabs(Tau_diff/Tau_avag))
            Ldist_io(2) = DABS( 1.0_wp+dabs(Tau_diff/Tau_avag))
            !write(*,*) 'test',Ret_io(1:2),Ldist_io(1:2)
            Ret_io(1:2)= Ret_io(1:2)*Ldist_io(1:2)
        END IF
        
        !==============averaged tauw based variables====================
        Tauw_ave_io = 0.5_wp*(DABS(Tauw_io(1))+DABS(Tauw_io(2)))
        DenAvew = DEN/YL
        VisAvew = VIS/YL
        
        !WRITE(*,*) 'DEN,YL', DEN, YL
        Utaw_ave_io = DSQRT( Tauw_ave_io/ DenAvew)
        Ret_ave_io  =  REN * Utaw_ave_io * DenAvew / VisAvew
        
        Tauw_d_ave_io = Tauw_ave_io * U0 * U0 * D0
        Utaw_d_ave_io = Utaw_ave_io * U0
        
        !WRITE(*,*) '# Ret_io', Ret_io(1), Ret_io(2) !test
        
        RETURN
    END SUBROUTINE
    
!================================================================
    SUBROUTINE PP_Wall_thermal_properties(flg_xzt)
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: flg_xzt
        INTEGER(4) :: J, JP, JM, NSZ, JJ, JJP, JJM
        REAL(WP)   :: HTEMP, CPtemp, Dtemp, Ktemp, Mtemp, Ttemp, Ptemp
        REAL(WP)   :: DUMMY(8,2)
        REAL(WP)   :: DUMMY_WORK(8,2)
        
        
        Hwal_RA(1:2) = 0.0_wp
        Hwal_FA(1:2) = 0.0_wp
        Twal(1:2) = 0.0_wp
        Dwal(1:2) = 0.0_wp
        Mwal(1:2) = 0.0_wp
        Kwal(1:2) = 0.0_wp
        Cpwal(1:2) = 0.0_wp
        qw(1:2) = 0.0_wp

        !=============wall parameters==undim==================================
        IF(BCWALLHEAT(itopwall)==isoThermalWall .AND. MYID==NPSLV) THEN
            HWAL_RA(itopwall) = H_WAL_GV(NCL1_IO/2, itopwall)
            Dwal(itopwall)    = D_WAL_GV(NCL1_IO/2, itopwall)
            Twal(itopwall)    = T_WAL_GV(NCL1_IO/2, itopwall)
            Mwal(itopwall)    = M_WAL_GV(NCL1_IO/2, itopwall)
            Kwal(itopwall)    = K_WAL_GV(NCL1_IO/2, itopwall)
            Cpwal(itopwall)   = Cp_WAL_GV(NCL1_IO/2, itopwall)
            
            HWAL_FA(itopwall)  = HWAL_RA(itopwall)
            
            J  =N2DO(MYID)
            JJ = JCL2G(J)
            if(flg_xzt==flgxz) then
                qw(itopwall)    = -1.0_wp*Kwal(itopwall)*(Twal(itopwall)-T1xzL_io (J))/( YND(NND2)-YCC(NCL2) )*CTHECD
            else if (flg_xzt==flgxzt) then
                qw(itopwall)    = -1.0_wp*Kwal(itopwall)*(Twal(itopwall)-T1xztL_io(J))/( YND(NND2)-YCC(NCL2) )*CTHECD
            else
            end if
 
        END IF    
        
        IF(BCWALLHEAT(ibotwall)==isoThermalWall .AND. MYID==0) THEN
            HWAL_RA(ibotwall)  = H_WAL_GV(NCL1_IO/2, ibotwall)
            Dwal(ibotwall)  = D_WAL_GV(NCL1_IO/2, ibotwall)
            Twal(ibotwall)  = T_WAL_GV(NCL1_IO/2, ibotwall)
            Mwal(ibotwall)  = M_WAL_GV(NCL1_IO/2, ibotwall)
            Kwal(ibotwall)  = K_WAL_GV(NCL1_IO/2, ibotwall)
            Cpwal(ibotwall) = Cp_WAL_GV(NCL1_IO/2, ibotwall)
            
            HWAL_FA(ibotwall)  = HWAL_RA(ibotwall)
            
            if(flg_xzt==flgxz) then
                qw(ibotwall)    = -1.0_wp*Kwal(ibotwall)*(T1xzL_io(1)-Twal(ibotwall))/(YCC(1)-YND(1))*CTHECD
            else if (flg_xzt==flgxzt) then
                qw(ibotwall)    = -1.0_wp*Kwal(ibotwall)*(T1xztL_io(1)-Twal(ibotwall))/(YCC(1)-YND(1))*CTHECD
            else
            end if
            
        END IF        
        
        IF(BCWALLHEAT(itopwall)==isoFluxWall .AND. MYID==NPSLV) THEN

            J =N2DO(MYID)
            JP=J
            JM=N2DO(MYID)-1
            
            JJ = JCL2G(J)
            JJP= JJ
            JJM= JJ-1
            if(flg_xzt==flgxz) then
                HWAL_RA(itopwall) =  2.0_wp*H1xzL_io(J) - YCL2ND_WFF(JJP)*H1xzL_io(JP) - YCL2ND_WFB(JJP)*H1xzL_io(JM)
                HWAL_FA(itopwall) = (2.0_wp*DHxzL_io(J) - YCL2ND_WFF(JJP)*DHxzL_io(JP) - YCL2ND_WFB(JJP)*DHxzL_io(JM))/&
                                    (2.0_wp*D1xzL_io(J) - YCL2ND_WFF(JJP)*D1xzL_io(JP) - YCL2ND_WFB(JJP)*D1xzL_io(JM))
            ELSE if (flg_xzt==flgxzt) then
                HWAL_RA(itopwall) =  2.0_wp*H1xztL_io(J) - YCL2ND_WFF(JJP)*H1xztL_io(JP) - YCL2ND_WFB(JJP)*H1xztL_io(JM)
                HWAL_FA(itopwall) = (2.0_wp*DHxztL_io(J) - YCL2ND_WFF(JJP)*DHxztL_io(JP) - YCL2ND_WFB(JJP)*DHxztL_io(JM))/&
                                    (2.0_wp*D1xztL_io(J) - YCL2ND_WFF(JJP)*D1xztL_io(JP) - YCL2ND_WFB(JJP)*D1xztL_io(JM))
            ELSE
            END IF            
            HTEMP = HWAL_RA(itopwall)
            
            IF(thermoStat==search_table) Then
                call NIST_SLEVAL_HT (Htemp,Ttemp)
                call NIST_SLEVAL_TD (Ttemp,Dtemp)
                call NIST_SLEVAL_TM (Ttemp,Mtemp)
                call NIST_SLEVAL_TK (Ttemp,Ktemp)
                call NIST_SLEVAL_TCP(Ttemp,CPtemp)
            ELSE IF(thermoStat==idealgas_law) Then
                if(flg_xzt==flgxz) then
                    Ptemp = 2.0_wp*U1xzL_io(J,4) - YCL2ND_WFF(JJP)*U1xzL_io(JP,4) - YCL2ND_WFB(JJP)*U1xzL_io(JM,4)
                ELSE if (flg_xzt==flgxzt) then
                    Ptemp = 2.0_wp*U1xztL_io(J,4) - YCL2ND_WFF(JJP)*U1xztL_io(JP,4) - YCL2ND_WFB(JJP)*U1xztL_io(JM,4)
                ELSE
                END IF  
                
                !call idealgas_HT (Htemp,Ttemp)
                !call idealgas_TD (Ttemp,Dtemp,Ptemp)
                !call idealgas_TM (Ttemp,Mtemp)
                !call idealgas_TK (Ttemp,Ktemp)
                !call idealgas_TCP(Ttemp,CPtemp)
            ELSE
            END IF
            
            Dwal(itopwall)  = Dtemp
            Twal(itopwall)  = Ttemp
            Mwal(itopwall)  = Mtemp
            Kwal(itopwall)  = Ktemp
            Cpwal(itopwall) = Cptemp
            if(flg_xzt==flgxz) then
                qw(itopwall)    = -1.0_wp*Kwal(itopwall)*(Twal(itopwall)-T1xzL_io(J))/( YND(NND2)-YCC(NCL2) )*CTHECD
            ELSE if (flg_xzt==flgxzt) then
                qw(itopwall)    = -1.0_wp*Kwal(itopwall)*(Twal(itopwall)-T1xztL_io(J))/( YND(NND2)-YCC(NCL2) )*CTHECD
            ELSE
            END IF
            
        END IF
        
        
        IF((BCWALLHEAT(ibotwall)==isoFluxWall .or. icase==IPIPEC ) .AND. MYID==0) THEN

            J=1
            JP=2
            JM=1
            
            if(flg_xzt==flgxz) then    
                HWAL_RA(ibotwall) = 2.0_wp*H1xzL_io(J) - YCL2ND_WFF(JP)*H1xzL_io(JP) - YCL2ND_WFB(JP)*H1xzL_io(JM)
                HWAL_FA(ibotwall) =(2.0_wp*DHxzL_io(J) - YCL2ND_WFF(JP)*DHxzL_io(JP) - YCL2ND_WFB(JP)*DHxzL_io(JM))/ &
                                    2.0_wp*D1xzL_io(J) - YCL2ND_WFF(JP)*D1xzL_io(JP) - YCL2ND_WFB(JP)*D1xzL_io(JM)
            ELSE if (flg_xzt==flgxzt) then
                HWAL_RA(ibotwall) = 2.0_wp*H1xztL_io(J) - YCL2ND_WFF(JP)*H1xztL_io(JP) - YCL2ND_WFB(JP)*H1xztL_io(JM)
                HWAL_FA(ibotwall) =(2.0_wp*DHxztL_io(J) - YCL2ND_WFF(JP)*DHxztL_io(JP) - YCL2ND_WFB(JP)*DHxztL_io(JM))/ &
                                    2.0_wp*D1xztL_io(J) - YCL2ND_WFF(JP)*D1xztL_io(JP) - YCL2ND_WFB(JP)*D1xztL_io(JM)
            ELSE
            END IF
                        
            HTEMP = HWAL_RA(ibotwall)
            
            IF(thermoStat==search_table) Then
                call NIST_SLEVAL_HT (Htemp,Ttemp)
                call NIST_SLEVAL_TD (Ttemp,Dtemp)
                call NIST_SLEVAL_TM (Ttemp,Mtemp)
                call NIST_SLEVAL_TK (Ttemp,Ktemp)
                call NIST_SLEVAL_TCP(Ttemp,CPtemp)
            ELSE IF(thermoStat==idealgas_law) Then
                if(flg_xzt==flgxz) then    
                    Ptemp = 2.0_wp*U1xzL_io(J,4) - YCL2ND_WFF(JP)*U1xzL_io(JP,4) - YCL2ND_WFB(JP)*U1xzL_io(JM,4)
                ELSE if (flg_xzt==flgxzt) then
                    Ptemp = 2.0_wp*U1xztL_io(J,4) - YCL2ND_WFF(JP)*U1xztL_io(JP,4) - YCL2ND_WFB(JP)*U1xztL_io(JM,4)
                ELSE
                END IF
                !call idealgas_HT (Htemp,Ttemp)
                !call idealgas_TD (Ttemp,Dtemp,Ptemp)
                !call idealgas_TM (Ttemp,Mtemp)
                !call idealgas_TK (Ttemp,Ktemp)
                !call idealgas_TCP(Ttemp,CPtemp)
            ELSE
            END IF
            
            Dwal(ibotwall)  = Dtemp
            Twal(ibotwall)  = Ttemp
            Mwal(ibotwall)  = Mtemp
            Kwal(ibotwall)  = Ktemp
            Cpwal(ibotwall) = Cptemp
            
            if(flg_xzt==flgxz) then
                qw(ibotwall)    = -1.0_wp*Kwal(ibotwall)*(T1xzL_io(1)-Twal(ibotwall))/(YCC(1)-YND(1))*CTHECD
            else if (flg_xzt==flgxzt) then
                qw(ibotwall)    = -1.0_wp*Kwal(ibotwall)*(T1xztL_io(1)-Twal(ibotwall))/(YCC(1)-YND(1))*CTHECD
            else
            end if
        END IF


        DUMMY(1,1:2) = Hwal_RA(1:2)
        DUMMY(2,1:2) = Hwal_FA(1:2)
        DUMMY(3,1:2) = Twal(1:2)
        DUMMY(4,1:2) = Dwal(1:2)
        DUMMY(5,1:2) = Mwal(1:2)
        DUMMY(6,1:2) = Kwal(1:2) 
        DUMMY(7,1:2) = Cpwal(1:2) 
        DUMMY(8,1:2) = qw(1:2)
        NSZ = 8*2
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(DUMMY, DUMMY_WORK, NSZ, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        Hwal_RA(1:2) = DUMMY_WORK(1,1:2)
        Hwal_FA(1:2) = DUMMY_WORK(2,1:2)
        Twal(1:2)    = DUMMY_WORK(3,1:2)
        Dwal(1:2)    = DUMMY_WORK(4,1:2)
        Mwal(1:2)    = DUMMY_WORK(5,1:2)
        Kwal(1:2)    = DUMMY_WORK(6,1:2)
        Cpwal(1:2)   = DUMMY_WORK(7,1:2)
        qw(1:2)      = DUMMY_WORK(8,1:2)
        
        RETURN
    END SUBROUTINE 
