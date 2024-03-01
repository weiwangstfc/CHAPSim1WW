    MODULE VARS_AVERAGED_TG
        USE wrt_info
        
        REAL(WP) :: Cf_LW_TG, Cf_UW_TG, Cf_ave_TG
        REAL(WP) :: U_tau_LW_TG, Re_tau_LW_TG, U_tau_ave_TG
        REAL(WP) :: U_tau_UW_TG, Re_tau_UW_TG, Re_tau_ave_TG
        
        CHARACTER(15) :: PNTIM
        
        !===============GLOABL DATA===================
        REAL(WP),ALLOCATABLE :: U1xztL_F0_tg( :, : )
        REAL(WP),ALLOCATABLE :: UPxztL_F0_tg( :, : )
        REAL(WP),ALLOCATABLE :: U2xztL_F0_tg( :, : )
        REAL(WP),ALLOCATABLE :: U3xztL_F0_tg( :, : )

        REAL(WP),ALLOCATABLE :: DVDL1xztL_F0_tg( :, :, : )
        REAL(WP),ALLOCATABLE :: DVDLPxztL_F0_tg( :, :, : )
        REAL(WP),ALLOCATABLE :: DVDL2xztL_F0_tg( :, :, : )
        
        REAL(WP),ALLOCATABLE :: U2PER( :, :, : )
        REAL(WP),ALLOCATABLE :: U3PER( :, :, :, :)
        REAL(WP),ALLOCATABLE :: VORper2( :,:)
        REAL(WP),ALLOCATABLE :: DUDX1(:,:,:)
        REAL(WP),ALLOCATABLE :: Skewness(:,:)
        
        REAL(WP),ALLOCATABLE :: BUDG_productn( :, :)
        REAL(WP),ALLOCATABLE :: BUDG_dissipat( :, :)
        REAL(WP),ALLOCATABLE :: BUDG_Tpr_diff( :, :)
        REAL(WP),ALLOCATABLE :: BUDG_Vis_diff( :, :)
        REAL(WP),ALLOCATABLE :: BUDG_VPG_diff( :, :)
        REAL(WP),ALLOCATABLE :: BUDG_VPG_stra( :, :)
        

    END MODULE
    
    !****************************************************************************
    SUBROUTINE MEMO_ALLOCT_AVERAGE_TG
        USE VARS_AVERAGED_TG
        use mesh_info
        use init_info
        IMPLICIT NONE

        ALLOCATE( U1xztL_F0_tg( NCL2, NDV+1 ) )
        ALLOCATE( UPxztL_F0_tg( NCL2, NDV   )  )
        ALLOCATE( U2xztL_F0_tg( NCL2, NDV*(7-NDV)/2+NDV-3 )   )
        ALLOCATE( U3xztL_F0_tg( NCL2, NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) )

        ALLOCATE( DVDL1xztL_F0_tg( NCL2, NDV, NDV  ) )
        ALLOCATE( DVDLPxztL_F0_tg( NCL2, NDV, NDV  ) )
        ALLOCATE( DVDL2xztL_F0_tg( NCL2, NDV*NDV, NDV*NDV  ) )
        
        ALLOCATE( U2PER(NCL2,NDV,NDV)                     )
        ALLOCATE( U3PER(NCL2,NDV,NDV,NDV)            )
        ALLOCATE( VORper2(NCL2,NDV)                     )
        ALLOCATE( DUDX1(NCL2,NDV,NDV)                     )
        ALLOCATE( Skewness(NCL2,NDV)                    )
        
        
        ALLOCATE( BUDG_productn( NCL2, (NDV*(7-NDV))/2+NDV-3 ) )
        ALLOCATE( BUDG_Tpr_diff( NCL2, (NDV*(7-NDV))/2+NDV-3 ) )
        ALLOCATE( BUDG_Vis_diff( NCL2, (NDV*(7-NDV))/2+NDV-3 ) )
        ALLOCATE( BUDG_VPG_diff( NCL2, (NDV*(7-NDV))/2+NDV-3 ) )
        ALLOCATE( BUDG_VPG_stra( NCL2, (NDV*(7-NDV))/2+NDV-3 ) )
        ALLOCATE( BUDG_dissipat( NCL2, (NDV*(7-NDV))/2+NDV-3 ) )
        
    
        RETURN
    END SUBROUTINE
    
!***********************************************
    SUBROUTINE MEMO_DEALLT_AVERAGE_TG
        USE VARS_AVERAGED_TG
        use mesh_info
        use init_info
        IMPLICIT NONE
        
        DEALLOCATE(U1xztL_F0_tg)
        DEALLOCATE(UPxztL_F0_tg)
        DEALLOCATE(U2xztL_F0_tg)
        DEALLOCATE(U3xztL_F0_tg)
        
        DEALLOCATE(DVDL1xztL_F0_tg)
        DEALLOCATE(DVDLPxztL_F0_tg)
        DEALLOCATE(DVDL2xztL_F0_tg)
        
        DEALLOCATE( U2PER     )
        DEALLOCATE( U3PER     )
        DEALLOCATE( VORper2    )
        DEALLOCATE( DUDX1     )
        DEALLOCATE( Skewness  )
        
        DEALLOCATE( BUDG_productn )
        DEALLOCATE( BUDG_Tpr_diff )
        DEALLOCATE( BUDG_Vis_diff )
        DEALLOCATE( BUDG_VPG_diff )
        DEALLOCATE( BUDG_VPG_stra )
        DEALLOCATE( BUDG_dissipat )
    
        RETURN
    END SUBROUTINE
    
    !===========================================================================
    SUBROUTINE WRT_AVERAGE_PPED_TG
        USE VARS_AVERAGED_TG
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE

        CALL WRT_AVERAGE_PPED_TG_GATHER
        IF(MYID==0) THEN
        
            CALL Cf_Utau_TG
            CALL WRT_AVERAGE_PPED_TG_CALC_RSTE_BUDG
            
            CALL WRT_AVERAGE_PPED_TG_WRT_TEC
            
            CALL MEMO_DEALLT_AVERAGE_TG
        END IF

        RETURN
    END SUBROUTINE
    
!===========================================================================
    SUBROUTINE WRT_AVERAGE_PPED_TG_WRT_TEC
        USE VARS_AVERAGED_TG
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4)    :: TECFLG1, TECFLG2
        REAL(WP)       :: COE1, COE2, COE
        REAL(WP)       :: urms, vrms, wrms, uv, uw, vw, yplus
        
        INTEGER(4)  :: INN
        INTEGER(4)  :: J, JJ
        INTEGER(4)  :: L, IP, M
        INTEGER(4)  :: N2DOID

        
        
        IF(MYID.NE.0) RETURN 

        !================WRITE DATA OUT=======================
        
        TECFLG1=200
        
        COE1 = DSQRT(dabs(Cf_LW_tg*0.5_WP))
        COE2 = DSQRT(dabs(Cf_UW_tg*0.5_WP))
         
        WRITE(PNTIM,'(1ES15.9)') phyTIME_tg
    
        OPEN(TECFLG1,FILE=TRIM(filepath4)//'Result.TG.Reynolds.Averaged.Flow.'//TRIM(PNTIM)//'.Profile.tec')
        WRITE(TECFLG1,'(A)') 'TITLE = " '//TRIM(ADJUSTL(teczonename))//' " '
        WRITE(TECFLG1,'(A)') &
            'VARIABLES = "Y", "Y+", "Ut", "Ux", "Uy", "Uz", "P", "Ruu", "Rvv", "Rww", "Ruv", "Ruw", "Rvw", ', &
            '"vorxper2", "voryper2", "vorzper2", "Suuu", "Svvv", "Swww", ', &
            '"uu_Prodc", "uu_Tdiff", "uu_Vdiff", "uu_VPGdiff", "uu_VPG-stain", "uu_Dissip", ', &
            '"vv_Prodc", "vv_Tdiff", "vv_Vdiff", "vv_VPGdiff", "vv_VPG-stain", "vv_Dissip", ', &
            '"ww_Prodc", "ww_Tdiff", "ww_Vdiff", "ww_VPGdiff", "ww_VPG-stain", "ww_Dissip", ', &
            '"uv_Prodc", "uv_Tdiff", "uv_Vdiff", "uv_VPGdiff", "uv_VPG-stain", "uv_Dissip", ', &
            '"uw_Prodc", "uw_Tdiff", "uw_Vdiff", "uw_VPGdiff", "uw_VPG-stain", "uw_Dissip", ', &
            '"vw_Prodc", "vw_Tdiff", "vw_Vdiff", "vw_VPGdiff", "vw_VPG-stain", "vw_Dissip"  ', &
            '"dmeanUdy", "dmeanWdy"'
        WRITE(TECFLG1,'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        
        
        DO J=1,NCL2
            if(YND(J).LT.0.0_wp) then
                COE = COE1
            else
                COE = COE2
            end if
            urms = DSQRT( dabs(U2PER(J,1,1) ) )
            vrms = DSQRT( dabs(U2PER(J,2,2) ) )
            wrms = DSQRT( dabs(U2PER(J,3,3) ) )
            uv   = U2PER(J,1,2)
            uw   = U2PER(J,1,3)
            vw   = U2PER(J,2,3)
            yplus= (1.0_wp-abs(YCC(J)))*REN*COE
            
        
            WRITE(TECFLG1,'(57ES22.14)') &
                YCC(J), yplus, COE, &
                U1xztL_F0_tg(J,1:4), &
                urms, vrms, wrms, uv, uw, vw, &
                VORper2(J,1), VORper2(J,2), VORper2(J,3), &
                Skewness(J,1:3), &
                BUDG_productn(J,1),BUDG_Tpr_diff(J,1),BUDG_Vis_diff(J,1), &
                BUDG_VPG_diff(J,1),BUDG_VPG_stra(J,1),BUDG_dissipat(J,1), &
                BUDG_productn(J,4),BUDG_Tpr_diff(J,4),BUDG_Vis_diff(J,4), &
                BUDG_VPG_diff(J,4),BUDG_VPG_stra(J,4),BUDG_dissipat(J,4), &
                BUDG_productn(J,6),BUDG_Tpr_diff(J,6),BUDG_Vis_diff(J,6), &
                BUDG_VPG_diff(J,6),BUDG_VPG_stra(J,6),BUDG_dissipat(J,6), &
                BUDG_productn(J,2),BUDG_Tpr_diff(J,2),BUDG_Vis_diff(J,2), &
                BUDG_VPG_diff(J,2),BUDG_VPG_stra(J,2),BUDG_dissipat(J,2), &
                BUDG_productn(J,3),BUDG_Tpr_diff(J,3),BUDG_Vis_diff(J,3), &
                BUDG_VPG_diff(J,3),BUDG_VPG_stra(J,3),BUDG_dissipat(J,3), &
                BUDG_productn(J,5),BUDG_Tpr_diff(J,5),BUDG_Vis_diff(J,5), &
                BUDG_VPG_diff(J,5),BUDG_VPG_stra(J,5),BUDG_dissipat(J,5), &
                DUDX1(J,1,2), DUDX1(J,3,2)
                                        
            !WRITE(*,*) YCC(J), VORper2(J,1)/COE1/COE1, VORper2(J,2)/COE1/COE1, VORper2(J,3)/COE1/COE1

            
        END DO
        CLOSE(TECFLG1)
        
        RETURN
    END SUBROUTINE
    
    
    SUBROUTINE WRT_AVERAGE_PPED_TG_CALC_RSTE_BUDG
        USE VARS_AVERAGED_TG
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4) :: J
        INTEGER(4) :: M, N, H, LMN, LMH, LNH, LMNH, L
        INTEGER(4) :: H1(3), P1(3), L1(3), L2(3)
        

        DO J=1,NCL2
        
            !====u'u', u'v', u'w', v'w', w'w'============
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    LMN = (M*(7-M))/2+N-3 
                    U2PER(J,M,N)=U2xztL_F0_tg(J,LMN)-U1xztL_F0_tg(J,M)*U1xztL_F0_tg(J,N)
                END DO
            END DO
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                        U2PER(J,M,N)=U2PER(J,N,M)
                    END IF
                END DO
            END DO
            
            !==above tested OK======
            
            !=====u'u'u', u'u'v', u'u'w', u'v'v', u'v'w', u'w'w', v'v'v', v'v'w', v'w'w', w'w'w ===
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    DO H=1,NDV
                        IF(N.GT.H) CYCLE
                        LMNH = M*(6-M)+(N*(7-N))/2+H-8  
                        LMN = (M*(7-M))/2+N-3
                        LMH = (M*(7-M))/2+H-3
                        LNH = (N*(7-N))/2+H-3
                        
                        U3PER(J, M,N,H) =  U3xztL_F0_tg(J,LMNH)  &
                                          -U1xztL_F0_tg(J,M) * U2xztL_F0_tg(J,LNH) &
                                          -U1xztL_F0_tg(J,N) * U2xztL_F0_tg(J,LMH) &
                                          -U1xztL_F0_tg(J,H) * U2xztL_F0_tg(J,LMN)  &
                                          +2.0_WP*U1xztL_F0_tg(J,M)*U1xztL_F0_tg(J,N)*U1xztL_F0_tg(J,H)
                    END DO
                END DO
            END DO
            
            !{111} {112} {113}
            ![121] {122} {123}
            ![131] [132] {133}
            ![211] [212] [213]
            ![221] {222} {223}
            ![231] [232] {233}
            ![311] [312] [313]
            ![321] [322] [323]
            ![331] [332] {333}
            
           
            U3PER(J,1,2,1) = U3PER(J,1,1,2)
            U3PER(J,1,3,1) = U3PER(J,1,1,3)
            U3PER(J,1,3,2) = U3PER(J,1,2,3)
            
            U3PER(J,2,1,1) = U3PER(J,1,1,2)
            U3PER(J,2,1,2) = U3PER(J,1,2,2)
            U3PER(J,2,1,3) = U3PER(J,1,2,3)
            
            U3PER(J,2,2,1) = U3PER(J,1,2,2)
            U3PER(J,2,3,1) = U3PER(J,1,2,3)
            U3PER(J,2,3,2) = U3PER(J,2,2,3)
            
            U3PER(J,3,1,1) = U3PER(J,1,1,3)
            U3PER(J,3,1,2) = U3PER(J,1,2,3)
            U3PER(J,3,1,3) = U3PER(J,1,3,3)
            
            U3PER(J,3,2,1) = U3PER(J,1,2,3) 
            U3PER(J,3,2,2) = U3PER(J,2,2,3)
            U3PER(J,3,2,3) = U3PER(J,2,3,3)
            
            U3PER(J,3,3,1) = U3PER(J,1,3,3)
            U3PER(J,3,3,2) = U3PER(J,2,3,3)
            
            !==========
            DO M=1,NDV
                DO N=1,NDV
                    IF(N==2) THEN
                        IF(J==1) THEN
                            DUDX1(J,M,N)= ( ( YCL2ND_WFB(J+1)*U1xztL_F0_tg(J,M) + YCL2ND_WFF(J+1)*U1xztL_F0_tg(J+1,M) ) - &
                                            0.0_WP  ) * DYFI(J)
                        ELSE IF (J==NCL2) THEN
                            DUDX1(J,M,N)= (   0.0_WP - &
                                          ( YCL2ND_WFF(J)  *U1xztL_F0_tg(J,M) + YCL2ND_WFB(J)  *U1xztL_F0_tg(J-1,M) )  ) * DYFI(J)
                        ELSE
                            DUDX1(J,M,N)= ( ( YCL2ND_WFB(J+1)*U1xztL_F0_tg(J,M) + YCL2ND_WFF(J+1)*U1xztL_F0_tg(J+1,M) ) - &
                                            ( YCL2ND_WFF(J)  *U1xztL_F0_tg(J,M) + YCL2ND_WFB(J)  *U1xztL_F0_tg(J-1,M) ) ) * DYFI(J)
                        END IF
                    ELSE
                        DUDX1(J,M,N) = 0.0_WP
                    END IF
                   
                END DO
            END DO
            
        END DO
        
        
        DO J=1,NCL2
            DO M=1,NDV
                Skewness(J,M) = U3PER(J,M,M,M)/ ( U2PER(J,M,M)**(3.0_wp/2.0_wp) )
            END DO
        END DO
        
        
        !===refer to note on 19/11/2014
        DO J=1, NCL2
            VORper2(J,1)= DVDL2xztL_F0_tg(J,8,8) + DVDL2xztL_F0_tg(J,6,6) -2.0_wp * DVDL2xztL_F0_tg(J,8,6)  &
                         +DUDX1(J,3,2)**2 + DUDX1(J,2,3)**2  &
                         -2.0_wp * DVDL1xztL_F0_tg(J,3,2) * DUDX1(J,3,2)  &
                         -2.0_wp * DVDL1xztL_F0_tg(J,2,3) * DUDX1(J,2,3)  &
                         +2.0_wp * DVDL1xztL_F0_tg(J,3,2) * DUDX1(J,2,3)  &
                         +2.0_wp * DVDL1xztL_F0_tg(J,2,3) * DUDX1(J,3,2)  &
                         -2.0_wp * DUDX1(J,3,2)*DUDX1(J,2,3) 
            VORper2(J,2)= DVDL2xztL_F0_tg(J,3,3) + DVDL2xztL_F0_tg(J,7,7) -2.0_wp * DVDL2xztL_F0_tg(J,3,7)  &
                         +DUDX1(J,1,3)**2 + DUDX1(J,3,1)**2  &
                         -2.0_wp * DVDL1xztL_F0_tg(J,1,3) * DUDX1(J,1,3)  &
                         -2.0_wp * DVDL1xztL_F0_tg(J,3,1) * DUDX1(J,3,1)  &
                         +2.0_wp * DVDL1xztL_F0_tg(J,1,3) * DUDX1(J,3,1)  &
                         +2.0_wp * DVDL1xztL_F0_tg(J,3,1) * DUDX1(J,1,3)  &
                         -2.0_wp * DUDX1(J,1,3)*DUDX1(J,3,1) 
            VORper2(J,3)= DVDL2xztL_F0_tg(J,4,4) + DVDL2xztL_F0_tg(J,2,2) -2.0_wp * DVDL2xztL_F0_tg(J,4,2)  &
                         +DUDX1(J,2,1)**2 + DUDX1(J,1,2)**2  &
                         -2.0_wp * DVDL1xztL_F0_tg(J,2,1) * DUDX1(J,2,1)  &
                         -2.0_wp * DVDL1xztL_F0_tg(J,1,2) * DUDX1(J,1,2)  &
                         +2.0_wp * DVDL1xztL_F0_tg(J,2,1) * DUDX1(J,1,2)  &
                         +2.0_wp * DVDL1xztL_F0_tg(J,1,2) * DUDX1(J,2,1)  &
                         -2.0_wp * DUDX1(J,2,1)*DUDX1(J,1,2) 
            VORper2(J,1) = dsqrt(dabs(VORper2(J,1)))
            VORper2(J,2) = dsqrt(dabs(VORper2(J,2)))
            VORper2(J,3) = dsqrt(dabs(VORper2(J,3)))
            
        END DO
        
        !============================================================================================================
        DO J=1,NCL2    
            !============================BUDGET TERMS=========================================================
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    !================PRODUCTUON TERMS=======================================
                    BUDG_productn(J,L)= U2PER(J,1,M)*DUDX1(J,N,1) + &
                                        U2PER(J,2,M)*DUDX1(J,N,2) + &
                                        U2PER(J,3,M)*DUDX1(J,N,3) + &
                                        U2PER(J,1,N)*DUDX1(J,M,1) + &
                                        U2PER(J,2,N)*DUDX1(J,M,2) + &
                                        U2PER(J,3,N)*DUDX1(J,M,3)
                    BUDG_productn(J,L)= BUDG_productn(J,L)*(-1.0_WP)               

                    !================DIFFUSION TERMS========================================
                    
                    !=================TURBULENCE TRANSPORT RATE==(T)==========
                    H = 2
                    IF(J==1) THEN
                        BUDG_Tpr_diff(J,L)= ( ( YCL2ND_WFF(J+1)*U3PER(J+1,M,N,H) + &
                                                YCL2ND_WFB(J+1)*U3PER(J,  M,N,H) ) - &
                                                0.0_WP ) * DYFI(J)
                    ELSE IF (J==NCL2) THEN
                        BUDG_Tpr_diff(J,L)= ( 0.0_WP - &
                                            ( YCL2ND_WFF(J)  *U3PER(J,  M,N,H) + &
                                              YCL2ND_WFB(J)  *U3PER(J-1,M,N,H) ) ) * DYFI(J)
                    ELSE
                        BUDG_Tpr_diff(J,L)= ( ( YCL2ND_WFF(J+1)*U3PER(J+1,M,N,H) + &
                                                YCL2ND_WFB(J+1)*U3PER(J,  M,N,H) ) - &
                                              ( YCL2ND_WFF(J)  *U3PER(J,  M,N,H) + &
                                                YCL2ND_WFB(J)  *U3PER(J-1,M,N,H) ) ) * DYFI(J)
                    END IF
                    BUDG_Tpr_diff(J,L) = -1.0_wp*BUDG_Tpr_diff(J,L)
                    
                    
                    !===============VISCOUS DIFFUSION TERM==(D)===============
                    
                    IF(J==1) THEN
                        BUDG_Vis_diff(J,L)= U2PER(J+1,M,N)*APVR(J,1) + U2PER(J,M,N)*ACVR(J,1) + 0.0_WP*AMVR(J,1)
                    ELSE IF(J==NCL2) THEN
                        BUDG_Vis_diff(J,L)= 0.0_WP*APVR(J,1) + U2PER(J,M,N)*ACVR(J,1) + U2PER(J-1,M,N)*AMVR(J,1)
                    ELSE
                        BUDG_Vis_diff(J,L)= U2PER(J+1,M,N)*APVR(J,1) + U2PER(J,  M,N)*ACVR(J,1) + U2PER(J-1,M,N)*AMVR(J,1) 
                    END IF
                    BUDG_Vis_diff(J,L) = BUDG_Vis_diff(J,L) *CVISC

                    !=============VELOCITY-PRESSURE GRADIENT (DIFFUSION)======
                    IF(M==2 .AND. N==2) THEN
                        IF(J==1) THEN
                            BUDG_VPG_diff(J,L)= &
                                ((  YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,M)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,M) )   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 -( YCL2ND_WFF(J)  *0.0_WP   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                  ) *DYFI(J)  &
                                 +(( YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,N)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,N) )   &
                                    +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                    -( YCL2ND_WFF(J)  *0.0_WP   &
                                    +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                  ) *DYFI(J) 
                        ELSE IF (J==NCL2) THEN
                            BUDG_VPG_diff(J,L)= &
                                ((  YCL2ND_WFF(J+1)*0.0_WP   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,M)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,M) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 ) *DYFI(J)  &
                                +(( YCL2ND_WFF(J+1)*0.0_WP   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,N)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,N) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 ) *DYFI(J) 
                        ELSE
                            BUDG_VPG_diff(J,L)=  &
                                 (( YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,M)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,M) )   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,M)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,M) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 ) *DYFI(J)  &
                                +(( YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,N)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,N) )   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,N)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,N) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 ) *DYFI(J) 
                        END IF
                    END IF
                    
                    IF((n.ne.2) .and. (m.eq.2)) THEN
                        IF(J==1) THEN
                            BUDG_VPG_diff(J,L)= &
                                 (( YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,N)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,N) )   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 -( YCL2ND_WFF(J)  *0.0_WP   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 ) *DYFI(J) 
                        ELSE IF (J==NCL2) THEN
                            BUDG_VPG_diff(J,L)= &
                                 (( YCL2ND_WFF(J+1)*0.0_WP   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,N)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,N) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 ) *DYFI(J) 
                        ELSE
                            BUDG_VPG_diff(J,L)=  &
                                 (( YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,N)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,N) )   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,N)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,N) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  N)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  N) ) ) &
                                 ) *DYFI(J) 

                        END IF
                    END IF
                    
                    IF((M.ne.2) .and. (n.eq.2)) THEN
                        IF(J==1) THEN
                            BUDG_VPG_diff(J,L)= &
                                ((  YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,M)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,M) )   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 -( YCL2ND_WFF(J)  *0.0_WP   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 ) *DYFI(J)
                        ELSE IF (J==NCL2) THEN
                            BUDG_VPG_diff(J,L)= &
                                ((  YCL2ND_WFF(J+1)*0.0_WP   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,M)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,M) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 ) *DYFI(J)
                        ELSE
                            BUDG_VPG_diff(J,L)=  &
                                 (( YCL2ND_WFF(J+1)*( UPxztL_F0_tg(J+1,M)-U1xztL_F0_tg(J+1,4)*U1xztL_F0_tg(J+1,M) )   &
                                   +YCL2ND_WFB(J+1)*( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 -( YCL2ND_WFF(J)  *( UPxztL_F0_tg(J-1,M)-U1xztL_F0_tg(J-1,4)*U1xztL_F0_tg(J-1,M) )   &
                                   +YCL2ND_WFB(J)  *( UPxztL_F0_tg(J,  M)-U1xztL_F0_tg(J,  4)*U1xztL_F0_tg(J,  M) ) ) &
                                 ) *DYFI(J)
                        END IF
                    END IF
                    
                    IF (.not.(M.eq.2 .or. N.eq.2)) THEN
                        BUDG_VPG_diff(J,L) = 0.0_WP
                    END IF
                    
                    BUDG_VPG_diff(J,L) = -1.0_wp*BUDG_VPG_diff(J,L)
                    
                    
                    !================VELOCITY-PRESSURE GRADIENT (STRAIN)==============================
                    BUDG_VPG_stra(J,L)= DVDLPxztL_F0_tg(J,M,N) - U1xztL_F0_tg(J,4)*DVDL1xztL_F0_tg(J,M,N) + &
                                        DVDLPxztL_F0_tg(J,N,M) - U1xztL_F0_tg(J,4)*DVDL1xztL_F0_tg(J,N,M)
                                     
                                     
                    !================DISSIPATION RATE=================================================
                    H1(1)=1
                    H1(2)=2
                    H1(3)=3
                    P1(1)=1
                    P1(2)=2
                    P1(3)=3
                    L1(1:3)=(M-1)*3+H1(1:3)
                    L2(1:3)=(N-1)*3+P1(1:3)
                    BUDG_dissipat(J,L)=  DVDL2xztL_F0_tg(J,L1(1),L2(1)) &
                                      +DVDL2xztL_F0_tg(J,L1(2),L2(2)) &
                                      +DVDL2xztL_F0_tg(J,L1(3),L2(3)) &
                                      -DVDL1xztL_F0_tg(J,N,1)*DUDX1(J,M,1) &
                                      -DVDL1xztL_F0_tg(J,N,2)*DUDX1(J,M,2) &
                                      -DVDL1xztL_F0_tg(J,N,3)*DUDX1(J,M,3) &
                                      -DVDL1xztL_F0_tg(J,M,1)*DUDX1(J,N,1) &
                                      -DVDL1xztL_F0_tg(J,M,2)*DUDX1(J,N,2) &
                                      -DVDL1xztL_F0_tg(J,M,3)*DUDX1(J,N,3) &
                                      +DUDX1(J,M,1)*DUDX1(J,N,1) &
                                      +DUDX1(J,M,2)*DUDX1(J,N,2) &
                                      +DUDX1(J,M,3)*DUDX1(J,N,3) 
                    BUDG_dissipat(J,L)=BUDG_dissipat(J,L)*CVISC*2.0_wp
                END DO
            END DO

        END DO
        

        RETURN
    END SUBROUTINE
    

!**************************************************************************** 
    SUBROUTINE Cf_Utau_TG
        USE VARS_AVERAGED_TG
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4) :: TECFLG = 200
        INTEGER(4) :: J
        real(wp)   :: Mdot, area_tg
        
        ! method 1 and 2 give slightly different results.
        ! method 2
        !Cf_LW_TG =  2.0_WP*(DVDL1xztL_F0_tg(1,   1,2))/REN
        !Cf_UW_TG = -2.0_WP*(DVDL1xztL_F0_tg(NCL2,1,2))/REN
        
        ! method 1
        Cf_UW_TG = -2.0_WP*(U1xztL_F0_tg(NCL2,1)-0.0_wp)/(YCC(NCL2)-YND(NND2)) /REN
        IF(ICASE==IPIPEC) THEN
            Cf_LW_TG = Cf_UW_TG
        ELSE
            Cf_LW_TG =  2.0_WP*(U1xztL_F0_tg(1,   1)-0.0_wp)/(YCC(1)   -YND(1))    /REN
        END IF

        U_tau_LW_TG = DSQRT(dabs(Cf_LW_TG)*0.5_wp)
        U_tau_UW_TG = DSQRT(dabs(Cf_UW_TG)*0.5_wp)
        
        Re_tau_LW_TG= REN * U_tau_LW_TG
        Re_tau_UW_TG= REN * U_tau_UW_TG
        
        
        Cf_ave_TG= 0.5*(dabs(Cf_LW_TG) + dabs(Cf_UW_TG))
        U_tau_ave_TG = DSQRT(dabs(Cf_ave_TG)*0.5_wp)
        Re_tau_ave_TG= REN * U_tau_ave_TG
        
        
        Mdot = 0.0_wp
        Area_TG = 0.0_wp
        DO J =1, NCL2
            Mdot = Mdot + U1xztL_F0_tg(J,1) /DYFI(J) /RCCI1(J)
            Area_TG = Area_TG + 1.0_WP/DYFI(J)/RCCI1(J)
        END DO
        Mdot = Mdot/Area_TG
        
        
        
        OPEN(TECFLG,FILE=TRIM(filepath4)//'Result.TG.Reynolds.Averaged.Cf.Utau.table.tec')
        
        IF(icase.ne.IPIPEC)  &
        Write(TECFLG,'(A,1ES15.8)') 'Cf on the lower wall=     ', Cf_LW_TG
        Write(TECFLG,'(A,1ES15.8)') 'Cf on the upper wall=     ', Cf_UW_TG
        Write(TECFLG,'(A,1ES15.8)') 'Cf average of walls =     ', Cf_ave_TG
        
        Write(TECFLG,'(A)') '      '
        
        IF(icase.ne.IPIPEC)  &
        Write(TECFLG,'(A,1ES15.8)') 'U_tau on the lower wall=  ', U_tau_LW_TG
        Write(TECFLG,'(A,1ES15.8)') 'U_tau on the upper wall=  ', U_tau_UW_TG
        Write(TECFLG,'(A,1ES15.8)') 'U_tau average of walls =  ', U_tau_ave_TG
        
        Write(TECFLG,'(A)') '      '
        
        IF(icase.ne.IPIPEC) &
        Write(TECFLG,'(A,1ES15.8)') 'Re_tau on the lower wall= ', Re_tau_LW_TG
        Write(TECFLG,'(A,1ES15.8)') 'Re_tau on the upper wall= ', Re_tau_UW_TG
        Write(TECFLG,'(A,1ES15.8)') 'Re_tau average of walls = ', Re_tau_ave_TG
        
        Write(TECFLG,'(A)') '      '
        Write(TECFLG,'(A,1ES15.8)') 'Bulk Mass Flux = ', Mdot
        
        CLOSE(TECFLG)
        
        
        
        RETURN
    END SUBROUTINE
    
    
    SUBROUTINE WRT_AVERAGE_PPED_TG_GATHER
        USE VARS_AVERAGED_TG
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4)    :: TECFLG1, TECFLG2
        REAL(WP)       :: COE1, COE2, COE
        !REAL(WP)       :: urms, vrms, wrms, uv, uw, vw
        
        INTEGER(4)  :: INN
        INTEGER(4)  :: J, JJ
        INTEGER(4)  :: L, IP, M, N, H, P, L1, L2
        INTEGER(4)  :: N2DOID
        REAL(WP)     :: D1AUX (N2DO(MYID), NDV+1,                               1:NPTOT)
        REAL(WP)     :: D2AUX (N2DO(MYID), NDV,                                 1:NPTOT)
        REAL(WP)     :: D3AUX (N2DO(MYID),(NDV*(7-NDV)/2+NDV-3),                1:NPTOT)
        REAL(WP)     :: D4AUX (N2DO(MYID),(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8),  1:NPTOT)
        
        REAL(WP)     :: T1AUX (N2DO(MYID), NDV,                 NDV, 1:NPTOT)
        REAL(WP)     :: T2AUX (N2DO(MYID), NDV,                 NDV, 1:NPTOT)
        REAL(WP)     :: T3AUX (N2DO(MYID), NDV*NDV, NDV*NDV, 1:NPTOT)
        
        
    
        INN = N2DO(MYID)*(NDV+1)
        CALL MPI_GATHER( U1xztL_tg, INN, MPI_DOUBLE_PRECISION, D1AUX, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*NDV
        CALL MPI_GATHER( UPxztL_tg, INN, MPI_DOUBLE_PRECISION, D2AUX, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*(NDV*(7-NDV)/2+NDV-3)
        CALL MPI_GATHER( U2xztL_tg, INN, MPI_DOUBLE_PRECISION, D3AUX, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        
        INN = N2DO(MYID)*(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
        CALL MPI_GATHER( U3xztL_tg, INN, MPI_DOUBLE_PRECISION, D4AUX, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        
        INN = N2DO(MYID)*NDV*NDV
        CALL MPI_GATHER( DVDL1xztL_tg, INN, MPI_DOUBLE_PRECISION, T1AUX, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*NDV*NDV
        CALL MPI_GATHER( DVDLPxztL_tg, INN, MPI_DOUBLE_PRECISION, T2AUX, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*(NDV*NDV)*NDV*NDV
        CALL MPI_GATHER( DVDL2xztL_tg, INN, MPI_DOUBLE_PRECISION, T3AUX, INN,  &
               MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        
        IF(MYID==0) THEN
            !==================================================================
            CALL MEMO_ALLOCT_AVERAGE_TG
            DO IP = 0, NPSLV
                N2DOID=JDEWT(IP)-JDSWT(IP)+1
                DO J=1,N2DOID
                    JJ=JDSWT(IP)-1+J
                  
                    DO L=1,NDV+1
                        U1xztL_F0_tg(JJ,L)=D1AUX(J,L,IP+1)
                    ENDDO
                  
                    DO L=1, NDV
                        UPxztL_F0_tg(JJ,L)=D2AUX(J,L,IP+1)
                    END DO
                  
                    DO L=1, (NDV*(7-NDV)/2+NDV-3)
                        U2xztL_F0_tg(JJ,L)=D3AUX(J,L,IP+1)
                    END DO
                  
                    DO L=1,(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
                        U3xztL_F0_tg(JJ,L)=D4AUX(J,L,IP+1)
                    END DO
                  
                    DO L=1,NDV
                        DO M =1,NDV
                            DVDL1xztL_F0_tg(JJ,L,M)=T1AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,NDV
                        DO M =1,NDV
                            DVDLPxztL_F0_tg(JJ,L,M)=T2AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,(NDV*(7-NDV)/2+NDV-3)
                        DO M =1,NDV
                            DVDL2xztL_F0_tg(JJ,L,M)=T3AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO M=1,NDV
                        DO N=1,NDV
                            DO H=1,NDV
                                DO P=1, NDV
                                    L1=(M-1)*3+H
                                    L2=(N-1)*3+P
                                        DVDL2xztL_F0_tg(JJ,L1,L2)=T3AUX(J,L1,L2,IP+1)
                                END DO
                            END DO
                        END DO
                    END DO

                ENDDO
            ENDDO
            
        END IF
    
        RETURN
    END SUBROUTINE
    
