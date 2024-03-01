    SUBROUTINE PP_MEAN_ZX_FLOW_TG
!>  @NOTE :  I, J,  L
!>           1, 1,  1
!>           1, 2,  2
!>           1, 3,  3
!>           2, 2,  4
!>           2, 3,  5
!>           3, 3,  6
!>  @NOTE :  I, J, K,  L
!>           1, 1, 1,  1
!>           1, 1, 2,  2
!>           1, 1, 3,  3
!>           1, 2, 2,  4
!>           1, 2, 3,  5
!>           1, 3, 3,  6
!>           2, 2, 2,  7
!>           2, 2, 3,  8
!>           2, 3, 3,  9
!>           3, 3, 3,  10

        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4)  :: I, J, K, IP, JP, KP, IM, JM, KM, JJ, JJP, JJM, KS, KSM, KSP
        INTEGER(4)  :: M, N, H, L, P, L1, L2
        REAL(WP)    :: U_cct(NDV)
        REAL(WP)    :: U_cct_ip(NDV), U_cct_jp(NDV), U_cct_kp(NDV)
        REAL(WP)    :: U_cct_im(NDV), U_cct_jm(NDV), U_cct_km(NDV)
        REAL(WP)    :: DVDL_cct(NDV,NDV)
        REAL(WP)    :: COE0, COE1, COE2, COE3
        REAL(WP)    :: COG1(NDV,NDV), COG2(NDV*NDV,NDV*NDV)
        
        U_cct = 0.0_WP
        U_cct_ip = 0.0_WP
        U_cct_jp = 0.0_WP
        U_cct_kp = 0.0_WP
        U_cct_im = 0.0_WP
        U_cct_jm = 0.0_WP
        U_cct_km = 0.0_WP
        
        COE0 = VL1313_tg
        COE1 = VL1313_tg*0.5_WP
        COE2 = VL1313_tg*(0.5_wp**2)
        COE3 = VL1313_tg*(0.5_wp**3)
        
        
        COG1(1,1) = VL1313_tg * DXI
        COG1(1,2) = VL1313_tg * 0.25_WP
        COG1(1,3) = VL1313_tg * 0.25_WP * DZI
        
        COG1(2,1) = VL1313_tg * 0.25_WP * DXI
        COG1(2,2) = VL1313_tg
        COG1(2,3) = VL1313_tg * 0.25_WP * DZI
        
        COG1(3,1) = VL1313_tg * 0.25_WP * DXI
        COG1(3,2) = VL1313_tg * 0.25_WP
        COG1(3,3) = VL1313_tg * DZI
        

!        DO M=1,NDV
!            DO N=1,NDV
!                DO H=1,NDV
!                    DO P=1,NDV
!                        L1=(M-1)*3+H
!                        L2=(N-1)*3+P
!                        COG2(L1,L2) = COG1(M,H)*COG1(N,P)
!                    END DO
!                END DO
!            END DO
!        END DO
        DO M=1,NDV
            DO H=1, NDV
                L1   = (M-1)*NDV+H
                DO N=1, NDV
                    DO P=1, NDV
                        L2=(N-1)*NDV+P
                        COG2(L1,L2) = COG1(M,H)*COG1(N,P)
                    END DO
                END DO
            END DO
        END DO
                    
        
        DO J=1,N2DO(MYID)
            JJ = JCL2G(J)
            JM = JLMV(J)
            JP = JLPV(J)
            JJP= JGPV(JJ)
            JJM= JGMV(JJ)
            
            U1xzL_tg(J,:) = 0.0_WP
            UPxzL_tg(J,:) = 0.0_WP
            
            U2xzL_tg(J,:) = 0.0_WP
            U3xzL_tg(J,:) = 0.0_WP
            
            DVDL1xzL_tg(J,:,:) = 0.0_WP
            DVDLPxzL_tg(J,:,:) = 0.0_WP
            DVDL2xzL_tg(J,:,:) = 0.0_WP
            
            DO  I=1,NCL1_tg
                IP = IPV_TG(I)
                IM = IMV_TG(I)
                DO K=1,NCL3

                    KP = KPV(K)
                    KM = KMV(K)
                    KS = KSYM(K)
                    KSM= KSYM(KM)
                    KSP= KSYM(KP)
                    
                    U_cct(1)    = Q_tg(I,J,K,1)  + Q_tg(IP,J,K,1)            ! U at i, j, k 
                    U_cct_jp(1) = Q_tg(I,JP,K,1) + Q_tg(IP,JP,K,1)           ! U at i, j+1, k 
                    U_cct_jm(1) = Q_tg(I,JM,K,1) + Q_tg(IP,JM,K,1)           ! U at i, j-1, k 
                    U_cct_kp(1) = Q_tg(I,J,KP,1) + Q_tg(IP,J,KP,1)           ! U at i, j, k+1 
                    U_cct_km(1) = Q_tg(I,J,KM,1) + Q_tg(IP,J,KM,1)           ! U at i, j, k-1
                    
                    U_cct(3)    = (Q_tg(I,J,K,3)  + Q_tg(I,J,KP,3))*RCCI1(JJ)   ! W at i, j, k 
                    U_cct_ip(3) = (Q_tg(IP,J,K,3) + Q_tg(IP,J,KP,3))*RCCI1(JJ)  ! W at i+1, j, k 
                    U_cct_im(3) = (Q_tg(IM,J,K,3) + Q_tg(IM,J,KP,3))*RCCI1(JJ)  ! W at i-1, j, k 
                    U_cct_jp(3) = (Q_tg(I,JP,K,3) + Q_tg(I,JP,KP,3))*RCCI1(JJP) ! W at i, j+1, k 
                    U_cct_jm(3) = (Q_tg(I,JM,K,3) + Q_tg(I,JM,KP,3))*RCCI1(JJM) ! W at i, j-1, k 
                    
                    IF(JJ==1 .AND. ICASE==IPIPEC) THEN
                        U_cct(2)    = (Q_tg(I,JP,K,2) - Q_tg(I,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_tg(I,JP,K,2)*RNDI1(JJP)  ! V at i, j, k 
                        U_cct_ip(2) = (Q_tg(IP,JP,K,2) - Q_tg(IP,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_tg(IP,JP,K,2)*RNDI1(JJP)  ! V at i+1, j, k 
                        U_cct_im(2) = (Q_tg(IM,JP,K,2) - Q_tg(IM,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_tg(IM,JP,K,2)*RNDI1(JJP)  ! V at i-1, j, k 
                        U_cct_kp(2) = (Q_tg(I,JP,KP,2) - Q_tg(I,JP,KSP,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_tg(I,JP,KP,2)*RNDI1(JJP)  ! V at i, j, k+1 
                        U_cct_km(2) = (Q_tg(I,JP,KM,2) - Q_tg(I,JP,KSM,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_tg(I,JP,KM,2)*RNDI1(JJP)  ! V at i, j, k-1
                    ELSE
                        U_cct(2)    = (Q_tg(I,J,K,2)  + Q_tg(I,JP,K,2) )*RCCI1(JJ) ! V at i, j, k 
                        U_cct_ip(2) = (Q_tg(IP,J,K,2) + Q_tg(IP,JP,K,2))*RCCI1(JJ) ! V at i+1, j, k
                        U_cct_im(2) = (Q_tg(IM,J,K,2) + Q_tg(IM,JP,K,2))*RCCI1(JJ) ! V at i-1, j, k
                        U_cct_kp(2) = (Q_tg(I,J,KP,2) + Q_tg(I,JP,KP,2))*RCCI1(JJ) ! V at i, j, k+1
                        U_cct_km(2) = (Q_tg(I,J,KM,2) + Q_tg(I,JP,KM,2))*RCCI1(JJ) ! V at i, j, k-1
                    END IF
                    
                    !=================MEAN U,V,W,P=========================
                    DO M = 1, NDV
                        U1xzL_tg(J,M) = U1xzL_tg(J,M) + U_cct(M)
                        UPxzL_tg(J,M) = UPxzL_tg(J,M) + U_cct(M) * PR_tg(I,J,K)
                    END DO
                    U1xzL_tg(J,NDV+1) = U1xzL_tg(J,NDV+1) + PR_tg(I,J,K)
                    
                    !====MEAN UU,UV,UW,ETC...===============================
                    DO M=1,NDV
                        DO N=1,NDV
                            IF(M.GT.N) CYCLE
                            L = (M*(7-M))/2+N-3  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER TRI-MATRIX
                            U2xzL_tg(J,L) = U2xzL_tg(J,L) +  U_cct(M) * U_cct(N)
                        END DO
                    END DO
                    
                    !======MEAN UUU, UUV, UUW, UVV, UVW ETC...================
                    DO M=1,NDV
                        DO N=1,NDV
                            IF(M.GT.N) CYCLE
                            DO H=1,NDV
                                IF(N.GT.H) CYCLE
                                L = M*(6-M)+(N*(7-N))/2+H-8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER TRI-MATRIX
                                U3xzL_tg(J,L) = U3xzL_tg(J,L) +  U_cct(M) * U_cct(N) * U_cct(H)
                            END DO
                        END DO
                    END DO
                    
                    
                    !========================DU/DX.DY.DZ======================================
                    DVDL_cct(1,1) =   Q_tg(IP,J,K,1) - Q_tg(I,J,K,1)      ! *DXI
                    DVDL_cct(1,2) = ( U_cct_jp(1) - U_cct(1) ) * DYCI(JJP) + &
                                    ( U_cct(1) - U_cct_jm(1) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                    DVDL_cct(1,3) =   U_cct_kp(1) - U_cct_km(1)           ! * 0.5_WP * 0.5_WP * DZI
                    
                    !========================DV/DX.DY.DZ====================================== 
                    DVDL_cct(2,1) =   U_cct_ip(2) - U_cct_im(2)           ! * 0.5_WP * 0.5_WP * DXI
                    IF(JJ==1 .AND. ICASE.EQ.IPIPEC) THEN
                        DVDL_cct(2,2) = (   Q_tg(I,JP,K,2)*RNDI1(JJP) - &
                                          ( Q_tg(I,JP,K,2) - Q_tg(I,JP,KS,2) ) * 0.50_WP*RNDI1(JJP) ) * DYFI(JJ)
                    ELSE
                        DVDL_cct(2,2) = ( Q_tg(I,JP,K,2)*RNDI1(JJP) - Q_tg(I,J,K,2)*RNDI1(JJ) ) * DYFI(JJ)
                    END IF
                    DVDL_cct(2,3) =   U_cct_kp(2) - U_cct_km(2)           ! * 0.5_WP * 0.5_WP * DZI
                    
                    !========================DV/DX.DY.DZ====================================== 
                    DVDL_cct(3,1) =   U_cct_ip(3) - U_cct_im(3)           ! * 0.5_WP * 0.5_WP * DXI
                    DVDL_cct(3,2) = ( U_cct_jp(3) - U_cct(3) ) * DYCI(JJP) + &
                                    ( U_cct(3) - U_cct_jm(3) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                    DVDL_cct(3,3) =   Q_tg(I,J,KP,3)*RCCI1(JJ) - Q_tg(I,J,K,3)*RCCI1(JJ)      ! *DzI            
                    
                    !===========mean dU/dX  ,  P dU/dX =======================================
                    DO M=1,NDV
                        DO N=1,NDV
                            DVDL1xzL_tg(J,M,N) =  DVDL1xzL_tg(J,M,N) + DVDL_cct(M,N)
                            DVDLPxzL_tg(J,M,N) =  DVDLPxzL_tg(J,M,N) + DVDL_cct(M,N) * PR_tg(I,J,K)
                        END DO
                    END DO
                    
                    !===========mean dUi/dX * dUj/dX ==========================================
!                    DO M=1,NDV
!                        DO N=1,NDV
!                            DO H=1,NDV
!                                DO P=1,NDV
!                                    L1=(M-1)*3+H
!                                    L2=(N-1)*3+P
!                                    DVDL2xzL_tg(J,L1,L2) =  DVDL2xzL_tg(J,L1,L2) + DVDL_cct(M,H)*DVDL_cct(N,P)
!                                END DO
!                            END DO
!                        END DO
!                    END DO
                    
                    DO M=1,NDV
                        DO H=1, NDV
                            L1   = (M-1)*NDV+H
                            DO N=1, NDV
                                DO P=1, NDV
                                    L2=(N-1)*NDV+P
                                    DVDL2xzL_tg (J,L1,L2)=  DVDL2xzL_tg (J,L1,L2) + DVDL_cct(M,H)*DVDL_cct(N,P)
                                END DO
                            END DO
                        END DO
                    END DO
                    

                ENDDO
            ENDDO
            U1xzL_tg(J,NDV+1) = U1xzL_tg(J,NDV+1) * COE0
            U1xzL_tg(J,1:NDV) = U1xzL_tg(J,1:NDV) * COE1
            
            UPxzL_tg(J,:) = UPxzL_tg(J,:) * COE1
            U2xzL_tg(J,:) = U2xzL_tg(J,:) * COE2
            U3xzL_tg(J,:) = U3xzL_tg(J,:) * COE3
            
            
            DO M=1,NDV
                DO N=1,NDV
                    DVDL1xzL_tg(J,M,N) = DVDL1xzL_tg(J,M,N) * COG1(M,N)
                    DVDLPxzL_tg(J,M,N) = DVDLPxzL_tg(J,M,N) * COG1(M,N)
                END DO
            END DO
            
!            DO M=1,NDV
!                DO N=1,NDV
!                    DO H=1,NDV
!                        DO P=1, NDV
!                            L1=(M-1)*3+H
!                            L2=(N-1)*3+P
!                            DVDL2xzL_tg(J,L1,L2) = DVDL2xzL_tg(J,L1,L2) * COG2(L1,L2)
!                        END DO
!                    END DO
!                END DO
!            END DO
            
            DO M=1,NDV
                DO H=1, NDV
                    L1   = (M-1)*NDV+H
                    DO N=1, NDV
                        DO P=1, NDV
                            L2=(N-1)*NDV+P
                            DVDL2xzL_tg(J,L1,L2) = DVDL2xzL_tg(J,L1,L2) * COG2(L1,L2)
                        END DO
                    END DO
                END DO
            END DO
            
            
        ENDDO
        
        
        !CALL PP_QUADRANTANALYSIS_xzL(ITG) 
        
        !============================================
        
        U1MEAN_tg = 0.0_WP
        U1MAXX_tg = 0.0_WP
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            U1MEAN_tg=U1MEAN_tg + U1xzL_tg(J,1)/DYFI(JJ)/RCCI1(JJ)
            U1MAXX_tg=DMAX1(U1MAXX_tg,DABS(U1xzL_tg(J,1))) 
        ENDDO

        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(U1MEAN_tg, U1MEAN_WORK_tg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        CALL MPI_ALLREDUCE(U1MAXX_tg, U1MAXX_WORK_tg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        U1MEAN_WORK_tg = U1MEAN_WORK_tg/DZI*DBLE(NCL3)/AREA_INLET
        
        RETURN
    END SUBROUTINE     


!***********************************************************************************
!***********************************************************************************
    SUBROUTINE PP_MEAN_ZX_FLOW_Xperiodic_IO
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
        
        INTEGER(4)  :: I, J, K, IP, JP, KP, IM, JM, KM, JJ, JJP, JJM, M, N, H, L, KS, KSM, KSP,P, L1, L2
        REAL(WP)    :: U_cct(NDV), G_cct(NDV)
        REAL(WP)    :: U_cct_ip(NDV), U_cct_jp(NDV), U_cct_kp(NDV)
        REAL(WP)    :: U_cct_im(NDV), U_cct_jm(NDV), U_cct_km(NDV)
        REAL(WP)    :: DVDL_cct(NDV,NDV)
        REAL(WP)    :: dHdL_cct(NDV), dTdL_cct(NDV)
        REAL(WP)    :: rtmp

        !if (myid==0) CALL CHKHDL('==> CALL PP_MEAN_ZX_FLOW_Xperiodic_IO', myid)
        DO J=1,N2DO(MYID)
            JP=JLPV(J)
            JM=JLMV(J)
            JJ=JCL2G(J)
            JJP=JGPV(JJ)
            JJM=JGMV(JJ)
            
            U1xzL_io(J,:) = 0.0_WP
            G1xzL_io(J,:) = 0.0_WP
            UPxzL_io (J,:) = 0.0_WP
            
            U2xzL_IO (J,:) = 0.0_WP
            UGxzL_IO (J,:) = 0.0_WP
            
            U3xzL_io(J,:) = 0.0_WP
            UGUxzL_IO(J,:) = 0.0_WP
            
            DVDL1xzL_io  (J,:,:) = 0.0_WP
            DVDLPxzL_io  (J,:,:) = 0.0_WP
            DVDL2xzL_io  (J,:,:) = 0.0_WP
            
            
            IF(thermlflg==1) THEN
                T1xzL_io(J)  = 0.0_WP
                D1xzL_io(J)  = 0.0_WP
                H1xzL_io(J)  = 0.0_WP
                M1xzL_io(J)  = 0.0_WP 
                 
                T2xzL_io(J)  = 0.0_WP
                D2xzL_io(J)  = 0.0_WP
                H2xzL_io(J)  = 0.0_WP
                
                DHxzL_io(J)  = 0.0_WP
                PHxzL_io(J)  = 0.0_WP
                
                DVDL1MxzL_io (J,:,:  ) = 0.0_WP
                DVDL1MHxzL_io(J,:,:  ) = 0.0_WP
                DVDL1MUxzL_io(J,:,:,:) = 0.0_WP
                DVDL2MxzL_io (J,:,:  ) = 0.0_WP
                
                UHxzL_io(J,:)  = 0.0_wp
                GHxzL_io(J,:)  = 0.0_wp
                U2DHxzL_IO(J,:)= 0.0_wp
                
                DhDL1xzL_io(J,:) = 0.0_wp
                DhDLPxzL_io(J,:) = 0.0_wp
                DTDLKxzL_io(J,:) = 0.0_wp
                DTDLKUxzL_io(J,:,:) = 0.0_wp
                DTDLKDVDLxzL_io(J,:,:,:) = 0.0_wp
                DHDLMDVDLxzL_io(J,:,:,:) = 0.0_wp
            END IF
                
            
            DO  I=1,NCL1_io
                IP=IPV_io(I)
                IM=IMV_io(I)
                
                DO K=1,NCL3
                    KP = KPV(K)
                    KM = KMV(K)
                    KS = KSYM(K)
                    KSP= KSYM(KP)
                    KSM= KSYM(KM)
                    
                    G_cct(1)    = ( G_io(I,J,K,1)  + G_io(IP,J,K,1)  ) *0.5_WP ! RHO*U at i, j, k 
                    U_cct(1)    = ( Q_io(I,J,K,1)  + Q_io(IP,J,K,1)  ) *0.5_WP ! U at i, j, k 
                    U_cct_jp(1) = ( Q_io(I,JP,K,1) + Q_io(IP,JP,K,1) ) *0.5_WP ! U at i, j+1, k 
                    U_cct_jm(1) = ( Q_io(I,JM,K,1) + Q_io(IP,JM,K,1) ) *0.5_WP ! U at i, j-1, k 
                    U_cct_kp(1) = ( Q_io(I,J,KP,1) + Q_io(IP,J,KP,1) ) *0.5_WP ! U at i, j, k+1 
                    U_cct_km(1) = ( Q_io(I,J,KM,1) + Q_io(IP,J,KM,1) ) *0.5_WP ! U at i, j, k-1
                    
                    G_cct(3)    = ( (G_io(I,J,K,3)  + G_io(I,J,KP,3))*RCCI1(JJ)   ) *0.5_WP  ! RHO*U at i, j, k 
                    U_cct(3)    = ( (Q_io(I,J,K,3)  + Q_io(I,J,KP,3))*RCCI1(JJ)   ) *0.5_WP  ! W at i, j, k 
                    U_cct_ip(3) = ( (Q_io(IP,J,K,3) + Q_io(IP,J,KP,3))*RCCI1(JJ)  ) *0.5_WP  ! W at i+1, j, k 
                    U_cct_im(3) = ( (Q_io(IM,J,K,3) + Q_io(IM,J,KP,3))*RCCI1(JJ)  ) *0.5_WP  ! W at i-1, j, k 
                    U_cct_jp(3) = ( (Q_io(I,JP,K,3) + Q_io(I,JP,KP,3))*RCCI1(JJP) ) *0.5_WP  ! W at i, j+1, k 
                    U_cct_jm(3) = ( (Q_io(I,JM,K,3) + Q_io(I,JM,KP,3))*RCCI1(JJM) ) *0.5_WP  ! W at i, j-1, k 
                    
                    
                    IF(JJ==1 .AND. ICASE==IPIPEC) THEN
                        G_cct(2)    = ( (G_io(I,JP,K,2) - G_io(I,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                         G_io(I,JP,K,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i, j, k 
                        U_cct(2)    = ( (Q_io(I,JP,K,2) - Q_io(I,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                         Q_io(I,JP,K,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i, j, k 
                        
                        U_cct_ip(2) = ( (Q_io(IP,JP,K,2) - Q_io(IP,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                         Q_io(IP,JP,K,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i+1, j, k 
                        U_cct_im(2) = ( (Q_io(IM,JP,K,2) - Q_io(IM,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                         Q_io(IM,JP,K,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i-1, j, k 
                        U_cct_kp(2) = ( (Q_io(I,JP,KP,2) - Q_io(I,JP,KSP,2))*0.50_WP*RNDI1(JJP) + &
                                         Q_io(I,JP,KP,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i, j, k+1 
                        U_cct_km(2) = ( (Q_io(I,JP,KM,2) - Q_io(I,JP,KSM,2))*0.50_WP*RNDI1(JJP) + &
                                         Q_io(I,JP,KM,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i, j, k-1
                    ELSE
                        G_cct(2)    = ( (G_io(I,J,K,2) + G_io(I,JP,K,2))*RCCI1(JJ)   ) *0.5_WP  ! RHO*U at i, j, k 
                        U_cct(2)    = ( (Q_io(I,J,K,2) + Q_io(I,JP,K,2))*RCCI1(JJ)   ) *0.5_WP  ! V at i, j, k 
                        U_cct_ip(2) = ( (Q_io(IP,J,K,2) + Q_io(IP,JP,K,2))*RCCI1(JJ) ) *0.5_WP  ! V at i+1, j, k
                        U_cct_im(2) = ( (Q_io(IM,J,K,2) + Q_io(IM,JP,K,2))*RCCI1(JJ) ) *0.5_WP  ! V at i-1, j, k
                        U_cct_kp(2) = ( (Q_io(I,J,KP,2) + Q_io(I,JP,KP,2))*RCCI1(JJ) ) *0.5_WP  ! V at i, j, k+1
                        U_cct_km(2) = ( (Q_io(I,J,KM,2) + Q_io(I,JP,KM,2))*RCCI1(JJ) ) *0.5_WP  ! V at i, j, k-1
                    END IF
                    
                    
                    !================U,V,W,P, G1,G2,G3========================
                    DO M = 1, NDV
                        U1xzL_io(J,M) = U1xzL_io(J,M) + U_cct(M)
                        G1xzL_io(J,M) = G1xzL_io(J,M) + G_cct(M)                 !++Method1++
                        !G1xzL_io(J,M) = G1xzL_io(J,M) + U_cct(M) * DENSITY(I,J,K)!++Method2++  for re-pp will introduce errors
                        UPxzL_io(J,M) = UPxzL_io(J,M) + U_cct(M) * PR_io(I,J,K)
                    END DO
                    U1xzL_io(J,NDV+1) = U1xzL_io(J,NDV+1) + PR_io(I,J,K)
                    
                    !====U*FLUX===============================================
                    DO M=1,NDV
                        DO N=1,NDV
                            IF(M.GT.N) CYCLE
                            L = (M*(7-M))/2+N-3  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER TRI-MATRIX
                            UGxzL_IO(J,L) = UGxzL_IO(J,L) +  U_cct(M) * G_cct(N)                 !++Method1++
                            !UGxzL_IO(J,L) = UGxzL_IO(J,L) +  U_cct(M) * U_cct(N) * DENSITY(I,J,K)!++Method2++
                            U2xzL_IO(J,L) = U2xzL_IO(J,L) +  U_cct(M) * U_cct(N)
                        END DO
                    END DO
                    
                    !====U*FLUX*U============================================
                    DO M=1,NDV
                        DO N=1,NDV
                            IF(M.GT.N) CYCLE
                            DO H=1,NDV
                                IF(N.GT.H) CYCLE
                                L = M*(6-M)+(N*(7-N))/2+H-8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER TRI-MATRIX
                                UGUxzL_IO(J,L) = UGUxzL_IO(J,L) +  U_cct(M) * G_cct(N) * U_cct(H)                 !++Method1++
                                !UGUxzL_IO(J,L) = UGUxzL_IO(J,L) +  U_cct(M) * U_cct(N) * U_cct(H) * DENSITY(I,J,K)!++Method2++
                            END DO
                        END DO
                    END DO
                    
                    !======MEAN UUU, UUV, UUW, UVV, UVW ETC...================
                    DO M=1,NDV
                        DO N=1,NDV
                            IF(M.GT.N) CYCLE
                            DO H=1,NDV
                                IF(N.GT.H) CYCLE
                                L = M*(6-M)+(N*(7-N))/2+H-8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER TRI-MATRIX
                                U3xzL_io(J,L) = U3xzL_io(J,L) +  U_cct(M) * U_cct(N) * U_cct(H)
                            END DO
                        END DO
                    END DO
                    ! U3(1) = U1 U1 U1
                    ! U3(2) = U1 U1 U2
                    ! U3(3) = U1 U1 U3
                    ! U3(4) = U1 U2 U2
                    ! U3(5) = U1 U2 U3
                    ! U3(6) = U1 U3 U3
                    ! U3(7) = U2 U2 U2
                    ! U3(8) = U2 U2 U3
                    ! U3(9) = U2 U3 U3
                    ! U3(10)= U3 U3 U3
                    
                    
                    
                    !========================DU/DX.DY.DZ======================================
                    DVDL_cct(1,1) = ( Q_io(IP,J,K,1) - Q_io(I,J,K,1) )*DXI  
                    DVDL_cct(1,3) = ( U_cct_kp(1) - U_cct_km(1) ) * DZI * 0.5_WP
                    DVDL_cct(1,2) = ( ( U_cct_jp(1) - U_cct(1)    ) * DYCI(JJP) + &
                                      ( U_cct(1)    - U_cct_jm(1) ) * DYCI(JJ) ) * 0.5_WP
                    
                    
                    !========================DV/DX.DY.DZ====================================== 
                    DVDL_cct(2,1) = ( U_cct_ip(2) - U_cct_im(2) ) * DXI * 0.5_WP 
                    DVDL_cct(2,3) = ( U_cct_kp(2) - U_cct_km(2) ) * DZI * 0.5_WP 
                    IF(JJ==1 .AND. ICASE.EQ.IPIPEC) THEN
                        DVDL_cct(2,2) = ( Q_io(I,JP,K,2)*RNDI1(JJP) - &
                                        ( Q_io(I,JP,K,2) - Q_io(I,JP,KS,2) ) * 0.50_WP*RNDI1(JJP) ) * DYFI(JJ)
                    ELSE
                        DVDL_cct(2,2) = ( Q_io(I,JP,K,2)*RNDI1(JJP) - Q_io(I,J,K,2)*RNDI1(JJ) ) * DYFI(JJ)
                    END IF
                    

                    !========================DV/DX.DY.DZ====================================== 
                    DVDL_cct(3,3) = ( Q_io(I,J,KP,3)*RCCI1(JJ) - Q_io(I,J,K,3)*RCCI1(JJ) ) *DZI 
                    DVDL_cct(3,1) = ( U_cct_ip(3) - U_cct_im(3) ) * DXI * 0.5_WP
                    DVDL_cct(3,2) = ( ( U_cct_jp(3) - U_cct(3)    ) * DYCI(JJP) + &
                                      ( U_cct(3)    - U_cct_jm(3) ) * DYCI(JJ) ) * 0.5_WP

                    !===============d(u_m)/d(x_n)===============================
                    ! Eq. DVDL1xzL_io(J,M,N) = d(u_m)/d(x_n)
                    !     DVDLPxzL_io(J,M,N) = p*d(u_m)/d(x_n)
                    DO M=1,NDV
                        DO N=1,NDV
                            DVDL1xzL_io(J,M,N) =  DVDL1xzL_io(J,M,N) + DVDL_cct(M,N)
                            DVDLPxzL_io(J,M,N) =  DVDLPxzL_io(J,M,N) + DVDL_cct(M,N) * PR_io(I,J,K)
                        END DO
                    END DO
                    
                    !======d(u_m)/d(x_h) * d(u_n)/d(x_p)=========================
                    ! Eq. DVDL2xzL_IO (J,(M-1)*3+H,(N-1)*3+P) = d(u_m)/d(x_h) * d(u_n)/d(x_p)
!                    DO M=1,NDV
!                        DO N=1,NDV
!                            DO H=1,NDV
!                                DO P=1,NDV
!                                    L1=(M-1)*3+H ! 1,2,3,4,5,6,7,8,9
!                                    L2=(N-1)*3+P ! 1,2,3,4,5,6,7,8,9
!                                    DVDL2xzL_IO (J,L1,L2)=  DVDL2xzL_IO (J,L1,L2) + DVDL_cct(M,H)*DVDL_cct(N,P)
!                                END DO
!                            END DO
!                        END DO
!                    END DO
                    !du_m/dx_h * du_n/dx_p
                    DO M=1,NDV
                        DO H=1, NDV
                            L1   = (M-1)*NDV+H
                            DO N=1, NDV
                                DO P=1, NDV
                                    L2=(N-1)*NDV+P
                                    DVDL2xzL_IO (J,L1,L2)=  DVDL2xzL_IO (J,L1,L2) + DVDL_cct(M,H)*DVDL_cct(N,P)
                                    
                                    ! m=1, h=1, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,1,1) = du/dx * du/dx  (J, (1-C1)*NDV+1, (1-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,1,2) = du/dx * du/dy  (J, (1-C1)*NDV+1, (1-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,1,3) = du/dx * du/dz  (J, (1-C1)*NDV+1, (1-C1)*NDV+3 )
                                    ! m=1, h=1, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,1,4) = du/dx * dv/dx  (J, (1-C1)*NDV+1, (2-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,1,5) = du/dx * dv/dy  (J, (1-C1)*NDV+1, (2-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,1,6) = du/dx * dv/dz  (J, (1-C1)*NDV+1, (2-C1)*NDV+3 )
                                    ! m=1, h=1, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,1,7) = du/dx * dw/dx  (J, (1-C1)*NDV+1, (3-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,1,8) = du/dx * dw/dy  (J, (1-C1)*NDV+1, (3-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,1,9) = du/dx * dw/dz  (J, (1-C1)*NDV+1, (3-C1)*NDV+3 )
                                    
                                    ! m=1, h=2, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,2,1) = du/dy * du/dx  (J, (1-C1)*NDV+2, (1-C1)*NDV+1 ) = DVDL2xzL_IO (J,1,2)
                                    ! DVDL2xzL_IO (J,2,2) = du/dy * du/dy  (J, (1-C1)*NDV+2, (1-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,2,3) = du/dy * du/dz  (J, (1-C1)*NDV+2, (1-C1)*NDV+3 )
                                    ! m=1, h=2, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,2,4) = du/dy * dv/dx  (J, (1-C1)*NDV+2, (2-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,2,5) = du/dy * dv/dy  (J, (1-C1)*NDV+2, (2-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,2,6) = du/dy * dv/dz  (J, (1-C1)*NDV+2, (2-C1)*NDV+3 )
                                    ! m=1, h=2, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,2,7) = du/dy * dw/dx  (J, (1-C1)*NDV+2, (3-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,2,8) = du/dy * dw/dy  (J, (1-C1)*NDV+2, (3-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,2,9) = du/dy * dw/dz  (J, (1-C1)*NDV+2, (3-C1)*NDV+3 )
                                    
                                    ! m=1, h=3, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,3,1) = du/dz * du/dx  (J, (1-C1)*NDV+3, (1-C1)*NDV+1 ) = DVDL2xzL_IO (J,1,3)
                                    ! DVDL2xzL_IO (J,3,2) = du/dz * du/dy  (J, (1-C1)*NDV+3, (1-C1)*NDV+2 ) = DVDL2xzL_IO (J,2,3)
                                    ! DVDL2xzL_IO (J,3,3) = du/dz * du/dz  (J, (1-C1)*NDV+3, (1-C1)*NDV+3 )
                                    ! m=1, h=3, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,3,4) = du/dz * dv/dx  (J, (1-C1)*NDV+3, (2-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,3,5) = du/dz * dv/dy  (J, (1-C1)*NDV+3, (2-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,3,6) = du/dz * dv/dz  (J, (1-C1)*NDV+3, (2-C1)*NDV+3 )
                                    ! m=1, h=3, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,3,7) = du/dz * dw/dx  (J, (1-C1)*NDV+3, (3-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,3,8) = du/dz * dw/dy  (J, (1-C1)*NDV+3, (3-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,3,9) = du/dz * dw/dz  (J, (1-C1)*NDV+3, (3-C1)*NDV+3 )
                                    
                                    ! m=2, h=1, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,4,1) = dv/dx * du/dx  (J, (2-C1)*NDV+1, (1-C1)*NDV+1 ) = DVDL2xzL_IO (J,1,4) 
                                    ! DVDL2xzL_IO (J,4,2) = dv/dx * du/dy  (J, (2-C1)*NDV+1, (1-C1)*NDV+2 ) = DVDL2xzL_IO (J,2,4)
                                    ! DVDL2xzL_IO (J,4,3) = dv/dx * du/dz  (J, (2-C1)*NDV+1, (1-C1)*NDV+3 ) = DVDL2xzL_IO (J,3,4)
                                    ! m=2, h=1, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,4,4) = dv/dx * dv/dx  (J, (2-C1)*NDV+1, (2-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,4,5) = dv/dx * dv/dy  (J, (2-C1)*NDV+1, (2-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,4,6) = dv/dx * dv/dz  (J, (2-C1)*NDV+1, (2-C1)*NDV+3 )
                                    ! m=2, h=1, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,4,7) = dv/dx * dw/dx  (J, (2-C1)*NDV+1, (3-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,4,8) = dv/dx * dw/dy  (J, (2-C1)*NDV+1, (3-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,4,9) = dv/dx * dw/dz  (J, (2-C1)*NDV+1, (3-C1)*NDV+3 )
                                    
                                    ! m=2, h=2, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,5,1) = dv/dy * du/dx  (J, (2-C1)*NDV+2, (1-C1)*NDV+1 ) = DVDL2xzL_IO (J,1,5)
                                    ! DVDL2xzL_IO (J,5,2) = dv/dy * du/dy  (J, (2-C1)*NDV+2, (1-C1)*NDV+2 ) = DVDL2xzL_IO (J,2,5)
                                    ! DVDL2xzL_IO (J,5,3) = dv/dy * du/dz  (J, (2-C1)*NDV+2, (1-C1)*NDV+3 ) = DVDL2xzL_IO (J,3,5)
                                    ! m=2, h=2, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,5,4) = dv/dy * dv/dx  (J, (2-C1)*NDV+2, (2-C1)*NDV+1 ) = DVDL2xzL_IO (J,4,5)
                                    ! DVDL2xzL_IO (J,5,5) = dv/dy * dv/dy  (J, (2-C1)*NDV+2, (2-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,5,6) = dv/dy * dv/dz  (J, (2-C1)*NDV+2, (2-C1)*NDV+3 )
                                    ! m=2, h=2, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,5,7) = dv/dy * dw/dx  (J, (2-C1)*NDV+2, (3-C1)*NDV+1 )
                                    ! DVDL2xzL_IO (J,5,8) = dv/dy * dw/dy  (J, (2-C1)*NDV+2, (3-C1)*NDV+2 )
                                    ! DVDL2xzL_IO (J,5,9) = dv/dy * dw/dz  (J, (2-C1)*NDV+2, (3-C1)*NDV+3 )
                                    
                                    ! m=2, h=3, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,6,1) = dv/dz * du/dx = DVDL2xzL_IO (J,1,6)
                                    ! DVDL2xzL_IO (J,6,2) = dv/dz * du/dy = DVDL2xzL_IO (J,2,6)
                                    ! DVDL2xzL_IO (J,6,3) = dv/dz * du/dz = DVDL2xzL_IO (J,3,6)
                                    ! m=2, h=3, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,6,4) = dv/dz * dv/dx = DVDL2xzL_IO (J,4,6)
                                    ! DVDL2xzL_IO (J,6,5) = dv/dz * dv/dy = DVDL2xzL_IO (J,5,6)
                                    ! DVDL2xzL_IO (J,6,6) = dv/dz * dv/dz
                                    ! m=2, h=3, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,6,7) = dv/dz * dw/dx
                                    ! DVDL2xzL_IO (J,6,8) = dv/dz * dw/dy
                                    ! DVDL2xzL_IO (J,6,9) = dv/dz * dw/dz
                                    
                                    ! m=3, h=1, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,7,1) = dw/dx * du/dx = DVDL2xzL_IO (J,1,7)
                                    ! DVDL2xzL_IO (J,7,2) = dw/dx * du/dy = DVDL2xzL_IO (J,2,7)
                                    ! DVDL2xzL_IO (J,7,3) = dw/dx * du/dz = DVDL2xzL_IO (J,3,7)
                                    ! m=3, h=1, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,7,4) = dw/dx * dv/dx = DVDL2xzL_IO (J,4,7)
                                    ! DVDL2xzL_IO (J,7,5) = dw/dx * dv/dy = DVDL2xzL_IO (J,5,7)
                                    ! DVDL2xzL_IO (J,7,6) = dw/dx * dv/dz = DVDL2xzL_IO (J,6,7)
                                    ! m=3, h=1, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,7,7) = dw/dx * dw/dx
                                    ! DVDL2xzL_IO (J,7,8) = dw/dx * dw/dy
                                    ! DVDL2xzL_IO (J,7,9) = dw/dx * dw/dz
                                    
                                    ! m=3, h=2, n=1, m=1:3
                                    ! DVDL2xzL_IO (J,8,1) = dw/dy * du/dx = DVDL2xzL_IO (J,1,8)
                                    ! DVDL2xzL_IO (J,8,2) = dw/dy * du/dy = DVDL2xzL_IO (J,2,8)
                                    ! DVDL2xzL_IO (J,8,3) = dw/dy * du/dz = DVDL2xzL_IO (J,3,8)
                                    ! m=3, h=2, n=2, m=1:3
                                    ! DVDL2xzL_IO (J,8,4) = dw/dy * dv/dx = DVDL2xzL_IO (J,4,8)
                                    ! DVDL2xzL_IO (J,8,5) = dw/dy * dv/dy = DVDL2xzL_IO (J,5,8)
                                    ! DVDL2xzL_IO (J,8,6) = dw/dy * dv/dz = DVDL2xzL_IO (J,6,8)
                                    ! m=3, h=2, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,8,7) = dw/dy * dw/dx = DVDL2xzL_IO (J,7,8)
                                    ! DVDL2xzL_IO (J,8,8) = dw/dy * dw/dy
                                    ! DVDL2xzL_IO (J,8,9) = dw/dy * dw/dz
                                    
                                    ! m=3, h=3, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,9,1) = dw/dz * du/dx = DVDL2xzL_IO (J,1,9)
                                    ! DVDL2xzL_IO (J,9,2) = dw/dz * du/dy = DVDL2xzL_IO (J,2,9)
                                    ! DVDL2xzL_IO (J,9,3) = dw/dz * du/dz = DVDL2xzL_IO (J,3,9)
                                    ! m=3, h=3, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,9,4) = dw/dz * dv/dx = DVDL2xzL_IO (J,4,9)
                                    ! DVDL2xzL_IO (J,9,5) = dw/dz * dv/dy = DVDL2xzL_IO (J,5,9)
                                    ! DVDL2xzL_IO (J,9,6) = dw/dz * dv/dz = DVDL2xzL_IO (J,6,9)
                                    ! m=3, h=3, n=3, m=1:3
                                    ! DVDL2xzL_IO (J,9,7) = dw/dz * dw/dx = DVDL2xzL_IO (J,7,9)
                                    ! DVDL2xzL_IO (J,9,8) = dw/dz * dw/dy = DVDL2xzL_IO (J,8,9)
                                    ! DVDL2xzL_IO (J,9,9) = dw/dz * dw/dz
                                END DO
                            END DO
                        END DO
                    END DO
                    
                    
                    
                    IF(thermlflg==1) THEN
                    
                        T1xzL_io(J) = T1xzL_io(J) + TEMPERATURE(I,J,K)
                        D1xzL_io(J) = D1xzL_io(J) + DENSITY    (I,J,K)
                        H1xzL_io(J) = H1xzL_io(J) + ENTHALPY   (I,J,K)
                        M1xzL_io(J) = M1xzL_io(J) + VISCOUSITY (I,J,K)
                        
                        T2xzL_io(J) = T2xzL_io(J) + TEMPERATURE(I,J,K)*TEMPERATURE(I,J,K)
                        D2xzL_io(J) = D2xzL_io(J) + DENSITY    (I,J,K)*DENSITY    (I,J,K)
                        H2xzL_io(J) = H2xzL_io(J) + ENTHALPY   (I,J,K)*ENTHALPY   (I,J,K)
                        
                        DHxzL_io(J) = DHxzL_io(J) + ENTHALPY(I,J,K)*DENSITY(I,J,K)
                        PHxzL_io(J) = PHxzL_io(J) + ENTHALPY(I,J,K)*PR_io(I,J,K)
                        
                        !=====d(h)/d(x_i)==============================================
                        DhDL_cct(1) = ( ENTHALPY(IP,J,K) - ENTHALPY(IM,J,K) ) * DXI * 0.5_WP
                        DhDL_cct(3) = ( ENTHALPY(I,J,KP) - ENTHALPY(I,J,KM) ) * DZI * 0.5_WP
                        DhDL_cct(2) = ( ( ENTHALPY(I,JP,K) - ENTHALPY(I,J, K ) )*DYCI(JJP) + &
                                        ( ENTHALPY(I,J, K) - ENTHALPY(I,JM,K ) )*DYCI(JJ) )*0.5_WP
                                        
                        !=====d(T)/d(x_i)==============================================
                        DTDL_cct(1) = ( TEMPERATURE(IP,J,K) - TEMPERATURE(IM,J,K) ) * DXI * 0.5_WP
                        DTDL_cct(3) = ( TEMPERATURE(I,J,KP) - TEMPERATURE(I,J,KM) ) * DZI * 0.5_WP
                        DTDL_cct(2) = ( ( TEMPERATURE(I,JP,K) - TEMPERATURE(I,J, K ) )*DYCI(JJP) + &
                                        ( TEMPERATURE(I,J, K) - TEMPERATURE(I,JM,K ) )*DYCI(JJ) )*0.5_WP
                                        
                        !===============d(u_m)/d(x_n)===============================
                        ! Eq. 
                        !     DVDL1MxzL_io(J,M,N) = M*d(u_m)/d(x_n)
                        !     DVDL1MhxzL_io(J,M,N)= M*d(u_m)/d(x_n)*H
                        !     DVDL1MUxzL_io(J,M,N,H) = M*d(u_m)/d(x_n) * u_h
                        DO M=1,NDV
                            DO N=1,NDV
                                DVDL1MxzL_io (J,M,N)= DVDL1MxzL_io(J,M,N) + DVDL_cct(M,N) * VISCOUSITY(I,J,K)
                                DVDL1MHxzL_io(J,M,N)= DVDL1MHxzL_io(J,M,N)+ DVDL_cct(M,N) * VISCOUSITY(I,J,K) * ENTHALPY(I,J,K)
                                DO H=1,NDV
                                    DVDL1MUxzL_io(J,M,N,H) =  DVDL1MUxzL_io(J,M,N,H) + DVDL_cct(M,N) * VISCOUSITY(I,J,K) * U_cct(H)
                                END DO
                            END DO
                        END DO
                        
                        !======d(u_m)/d(x_h) * d(u_n)/d(x_p)=========================
                        ! Eq. 
                        !     DVDL2MxzL_IO(J,(M-1)*3+H,(N-1)*3+P) = d(u_m)/d(x_h) * d(u_n)/d(x_p) * M
!                        DO M=1,NDV
!                            DO N=1,NDV
!                                DO H=1,NDV
!                                    DO P=1,NDV
!                                        L1=(M-1)*3+H ! 1,2,3,4,5,6,7,8,9
!                                        L2=(N-1)*3+P ! 1,2,3,4,5,6,7,8,9
!                                        DVDL2MxzL_IO(J,L1,L2)=  DVDL2MxzL_IO(J,L1,L2) + &
!                                                                DVDL_cct(M,H)*DVDL_cct(N,P) * VISCOUSITY(I,J,K)
!                                    END DO
!                                END DO
!                            END DO
!                        END DO
                        
                        DO M=1,NDV
                            DO H=1, NDV
                                L1   = (M-1)*NDV+H
                                DO N=1, NDV
                                    DO P=1, NDV
                                        L2=(N-1)*NDV+P
                                        DVDL2MxzL_IO(J,L1,L2)=  DVDL2MxzL_IO(J,L1,L2) + &
                                                                DVDL_cct(M,H)*DVDL_cct(N,P) * VISCOUSITY(I,J,K)
                                    END DO
                                END DO
                            END DO
                        END DO
                    
                    
                        !=========================================
                        !Eq.  UHxzL_io(J,M) = u_m * h
                        !     GHxzL_io(J,M) = g_m * h
                        !     DhDL1xzL_io(J,M)        = d(h)/d(x_m)
                        !     DhDLPxzL_io(J,M)        = d(h)/d(x_m) * p
                        !     DhDLKxzL_io(J,M)        = d(h)/d(x_m) * k
                        !     DTDLKUxzL_io(J,M,N)     = d(h)/d(x_m) * k * u_n
                        !     DTDLKDVDLxzL_io(J,M,N,H)= d(h)/d(x_m) * k * d(u_n)/d(x_h)
                        !     DHDLMDVDLxzL_io(J,M,N,H)= d(h)/d(x_m) * mu* d(u_n)/d(x_h)
                        DO M=1,NDV
                        
                            UHxzL_io(J,M) = UHxzL_io(J,M) + ENTHALPY(I,J,K)*U_cct(M) 
                            GHxzL_io(J,M) = GHxzL_io(J,M) + ENTHALPY(I,J,K)*G_cct(M)  !++method1+++
                            !GHxzL_io(J,M) = GHxzL_io(J,M) + ENTHALPY(I,J,K)*U_cct(M)*DENSITY(I,J,K) ! ++method2+++
                            
                            DhDL1xzL_io(J,M) =  DhDL1xzL_io(J,M) + DhDL_cct(M)   
                            DhDLPxzL_io(J,M) =  DhDLPxzL_io(J,M) + DhDL_cct(M) * PR_io(I,J,K) 
                            DTDLKxzL_io(J,M) =  DTDLKxzL_io(J,M) + DTDL_cct(M) * THERMCONDT(I,J,K) 
                            
                            DO N=1, NDV
                                DTDLKUxzL_io(J,M,N) =  DTDLKUxzL_io(J,M,N) + DTDL_cct(M) * THERMCONDT(I,J,K) * U_cct(N)
                                DO H=1, NDV
                                    DTDLKDVDLxzL_io(J,M,N,H) = DTDLKDVDLxzL_io(J,M,N,H) + &
                                                               DTDL_cct(M) * THERMCONDT(I,J,K)*DVDL_cct(N,H)
                                    DHDLMDVDLxzL_io(J,M,N,H) = DHDLMDVDLxzL_io(J,M,N,H) + &
                                                               DHDL_cct(M) * VISCOUSITY(I,J,K)*DVDL_cct(N,H)
                                END DO
                            END DO
                        END DO

                        !====\rho h *U*U============================================
                        !U2DHxzL_IO(J,L) = \rho * h * u_m * u_n
                        DO M=1,NDV
                            DO N=1,NDV
                                IF(M.GT.N) CYCLE
                                L = (M*(7-M))/2+N-3
                                U2DHxzL_IO(J,L)= U2DHxzL_IO(J,L) + DENSITY(I,J,K)*ENTHALPY(I,J,K)*U_cct(M)*U_cct(N)
                            END DO
                        END DO
                        
                    END IF
                    
                    
                ENDDO
            ENDDO
            
            U1xzL_io(J,:)  = U1xzL_io(J,:)*VL1313_io
            G1xzL_io(J,:)  = G1xzL_io(J,:)*VL1313_io
            UPxzL_io (J,:) = UPxzL_io (J,:)*VL1313_io
            
            U2xzL_IO (J,:) = U2xzL_IO (J,:)*VL1313_io
            UGxzL_IO (J,:) = UGxzL_IO (J,:)*VL1313_io
            U3xzL_io(J,:)  = U3xzL_io(J,:) *VL1313_io
            UGUxzL_IO(J,:) = UGUxzL_IO(J,:)*VL1313_io
            
            DVDL1xzL_io  (J,:,:)  = DVDL1xzL_io  (J,:,:)*VL1313_io
            DVDLPxzL_io  (J,:,:)  = DVDLPxzL_io  (J,:,:)*VL1313_io
            DVDL2xzL_io  (J,:,:)  = DVDL2xzL_io  (J,:,:)*VL1313_io
            
            !WRITE(*,*) J, U1xzL_io(J,2), G1xzL_io(J,2) !test
            IF(thermlflg==1) THEN
                T1xzL_io(J)  = T1xzL_io(J)*VL1313_io
                D1xzL_io(J)  = D1xzL_io(J)*VL1313_io
                H1xzL_io(J)  = H1xzL_io(J)*VL1313_io
                M1xzL_io(J)  = M1xzL_io(J)*VL1313_io
                
               
                T2xzL_io(J)  = T2xzL_io(J)*VL1313_io
                D2xzL_io(J)  = D2xzL_io(J)*VL1313_io
                H2xzL_io(J)  = H2xzL_io(J)*VL1313_io
                
                DHxzL_io(J)  = DHxzL_io(J)*VL1313_io
                PHxzL_io(J)  = PHxzL_io(J)*VL1313_io
                
                DVDL1MxzL_io (J,:,:)  = DVDL1MxzL_io (J,:,:)  *VL1313_io
                DVDL1MHxzL_io(J,:,:)  = DVDL1MHxzL_io(J,:,:)  *VL1313_io
                DVDL1MUxzL_io(J,:,:,:)= DVDL1MUxzL_io(J,:,:,:)*VL1313_io
                DVDL2MxzL_io (J,:,:)  = DVDL2MxzL_io (J,:,:)  *VL1313_io
                
                !if (myid==0) THEN
!                DO L1=1, (3-1)*NDV+3
!                    DO L2=1, (3-1)*NDV+3
!                        WRITE(*,*) J, L1, L2, DVDL2xzL_IO(J,L1,L2), DVDL2MxzL_IO(J,L1,L2)
!                    END DO
!                END DO
                    !WRITE(*,*) J, DVDLPxzL_io(J,1,1),DVDL1xzL_io  (J,1,1), U1xzL_io(J,4)
                !END IF
                        
                UHxzL_io(J,:)  = UHxzL_io(J,:)*VL1313_io
                GHxzL_io(J,:)  = GHxzL_io(J,:)*VL1313_io
                U2DHxzL_IO(J,:)= U2DHxzL_IO(J,:)*VL1313_io
                
                DhDL1xzL_io(J,:) = DhDL1xzL_io(J,:)*VL1313_io
                DhDLPxzL_io(J,:) = DhDLPxzL_io(J,:)*VL1313_io
                DTDLKxzL_io(J,:) = DTDLKxzL_io(J,:)*VL1313_io
                DTDLKUxzL_io(J,:,:)      = DTDLKUxzL_io(J,:,:)*VL1313_io
                DTDLKDVDLxzL_io(J,:,:,:) = DTDLKDVDLxzL_io(J,:,:,:)*VL1313_io
                DHDLMDVDLxzL_io(J,:,:,:) = DHDLMDVDLxzL_io(J,:,:,:)*VL1313_io
            END IF
            
        ENDDO
        
        
        !CALL PP_QUADRANTANALYSIS_xzL(IIO) 
        
        CALL PP_DRIVEN_FORCE
        !==============================
        
        G1RATE_IO = 0.0_WP
        G1MAXX_IO = 0.0_WP
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            G1RATE_IO=G1RATE_IO + G1xzL_io(J,1)/DYFI(JJ)/RCCI1(JJ)
            G1MAXX_IO=DMAX1(G1MAXX_IO,DABS(G1xzL_io(J,1))) 
        ENDDO

        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(G1RATE_IO, G1RATE_WORK_IO, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        CALL MPI_ALLREDUCE(G1MAXX_IO, G1MAXX_WORK_IO, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        G1BULK_WORK_IO = G1RATE_WORK_IO/DZI*DBLE(NCL3)/AREA_INLET
        
        IF(thermlflg==1) THEN
            H1RATE_io = 0.0_WP
            T1MAXX_io = 0.0_WP
            DO J=1,N2DO(MYID)  !@
                JJ=JCL2G(J)
                H1RATE_io=H1RATE_io + GHxzL_io(J,1)/DYFI(JJ)/RCCI1(JJ)
                T1MAXX_io=DMAX1(T1MAXX_io,DABS(T1xzL_io(J))) 
            ENDDO
            CALL MPI_BARRIER(ICOMM,IERROR)
            CALL MPI_ALLREDUCE(H1RATE_io, H1RATE_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
            CALL MPI_ALLREDUCE(T1MAXX_io, T1MAXX_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
            
            H1bulk_WORK_io =  H1RATE_WORK_io/G1RATE_WORK_io
            
            IF(thermoStat==search_table) Then
                call NIST_SLEVAL_HT(H1bulk_WORK_io,T1bulk_WORK_io)
            ELSE IF(thermoStat==idealgas_law) Then
                !call idealgas_HT(H1bulk_WORK_io,T1bulk_WORK_io)
            ELSE
            END IF
        
        END IF
        
        
        RETURN
    END SUBROUTINE   
    
    
    SUBROUTINE INST_Qio_GRADIENT_CELLCENTRED(I,J,K,DVDL_cct)
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
        
        INTEGER(4), INTENT(IN) :: I,J,K
        REAL(WP),INTENT(OUT)   :: DVDL_cct(NDV,NDV)
        
        INTEGER(4) :: JJ,JM, JP, JJP, JJM, IP, IM, KP, KM, KS, KSP, KSM
        REAL(WP)    :: U_cct(NDV)
        REAL(WP)    :: U_cct_ip(NDV), U_cct_jp(NDV), U_cct_kp(NDV)
        REAL(WP)    :: U_cct_im(NDV), U_cct_jm(NDV), U_cct_km(NDV)
        

        JJ = JCL2G(J)
        JM = JLMV(J)
        JP = JLPV(J)
        JJP= JGPV(JJ)
        JJM= JGMV(JJ)
        
        IP=IPV_io(I)
        IM=IMV_io(I)
        KP = KPV(K)
        KM = KMV(K)
        KS = KSYM(K)
        KSP= KSYM(KP)
        KSM= KSYM(KM)
                    
        U_cct(1)    = ( Q_io(I,J,K,1)  + Q_io(IP,J,K,1)  ) *0.5_WP ! U at i, j, k 
        U_cct_jp(1) = ( Q_io(I,JP,K,1) + Q_io(IP,JP,K,1) ) *0.5_WP ! U at i, j+1, k 
        U_cct_jm(1) = ( Q_io(I,JM,K,1) + Q_io(IP,JM,K,1) ) *0.5_WP ! U at i, j-1, k 
        U_cct_kp(1) = ( Q_io(I,J,KP,1) + Q_io(IP,J,KP,1) ) *0.5_WP ! U at i, j, k+1 
        U_cct_km(1) = ( Q_io(I,J,KM,1) + Q_io(IP,J,KM,1) ) *0.5_WP ! U at i, j, k-1
        
        U_cct(3)    = ( (Q_io(I,J,K,3)  + Q_io(I,J,KP,3))*RCCI1(JJ)   ) *0.5_WP  ! W at i, j, k 
        U_cct_ip(3) = ( (Q_io(IP,J,K,3) + Q_io(IP,J,KP,3))*RCCI1(JJ)  ) *0.5_WP  ! W at i+1, j, k 
        U_cct_im(3) = ( (Q_io(IM,J,K,3) + Q_io(IM,J,KP,3))*RCCI1(JJ)  ) *0.5_WP  ! W at i-1, j, k 
        U_cct_jp(3) = ( (Q_io(I,JP,K,3) + Q_io(I,JP,KP,3))*RCCI1(JJP) ) *0.5_WP  ! W at i, j+1, k 
        U_cct_jm(3) = ( (Q_io(I,JM,K,3) + Q_io(I,JM,KP,3))*RCCI1(JJM) ) *0.5_WP  ! W at i, j-1, k 
        
        
        IF(JJ==1 .AND. ICASE==IPIPEC) THEN
            U_cct(2)    = ( (Q_io(I,JP,K,2) - Q_io(I,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                             Q_io(I,JP,K,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i, j, k 
            
            U_cct_ip(2) = ( (Q_io(IP,JP,K,2) - Q_io(IP,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                             Q_io(IP,JP,K,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i+1, j, k 
            U_cct_im(2) = ( (Q_io(IM,JP,K,2) - Q_io(IM,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                             Q_io(IM,JP,K,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i-1, j, k 
            U_cct_kp(2) = ( (Q_io(I,JP,KP,2) - Q_io(I,JP,KSP,2))*0.50_WP*RNDI1(JJP) + &
                             Q_io(I,JP,KP,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i, j, k+1 
            U_cct_km(2) = ( (Q_io(I,JP,KM,2) - Q_io(I,JP,KSM,2))*0.50_WP*RNDI1(JJP) + &
                             Q_io(I,JP,KM,2)*RNDI1(JJP)  ) *0.5_WP  ! V at i, j, k-1
        ELSE
            U_cct(2)    = ( (Q_io(I,J,K,2) + Q_io(I,JP,K,2))*RCCI1(JJ)   ) *0.5_WP  ! V at i, j, k 
            U_cct_ip(2) = ( (Q_io(IP,J,K,2) + Q_io(IP,JP,K,2))*RCCI1(JJ) ) *0.5_WP  ! V at i+1, j, k
            U_cct_im(2) = ( (Q_io(IM,J,K,2) + Q_io(IM,JP,K,2))*RCCI1(JJ) ) *0.5_WP  ! V at i-1, j, k
            U_cct_kp(2) = ( (Q_io(I,J,KP,2) + Q_io(I,JP,KP,2))*RCCI1(JJ) ) *0.5_WP  ! V at i, j, k+1
            U_cct_km(2) = ( (Q_io(I,J,KM,2) + Q_io(I,JP,KM,2))*RCCI1(JJ) ) *0.5_WP  ! V at i, j, k-1
        END IF


        !========================DU/DX.DY.DZ======================================
        DVDL_cct(1,1) = ( Q_io(IP,J,K,1) - Q_io(I,J,K,1) )*DXI  
        DVDL_cct(1,3) = ( U_cct_kp(1) - U_cct_km(1) ) * DZI * 0.5_WP
        DVDL_cct(1,2) = ( ( U_cct_jp(1) - U_cct(1)    ) * DYCI(JJP) + &
                          ( U_cct(1)    - U_cct_jm(1) ) * DYCI(JJ) ) * 0.5_WP
        
        
        !========================DV/DX.DY.DZ====================================== 
        DVDL_cct(2,1) = ( U_cct_ip(2) - U_cct_im(2) ) * DXI * 0.5_WP 
        DVDL_cct(2,3) = ( U_cct_kp(2) - U_cct_km(2) ) * DZI * 0.5_WP 
        IF(JJ==1 .AND. ICASE.EQ.IPIPEC) THEN
            DVDL_cct(2,2) = (   Q_io(I,JP,K,2)*RNDI1(JJP) - &
                            ( Q_io(I,JP,K,2) - Q_io(I,JP,KS,2) ) * 0.50_WP*RNDI1(JJP) ) * DYFI(JJ)
        ELSE
            DVDL_cct(2,2) = ( Q_io(I,JP,K,2)*RNDI1(JJP) - Q_io(I,J,K,2)*RNDI1(JJ) ) * DYFI(JJ)
        END IF
        

        !========================DW/DX.DY.DZ====================================== 
        DVDL_cct(3,3) = ( Q_io(I,J,KP,3)*RCCI1(JJ) - Q_io(I,J,K,3)*RCCI1(JJ) ) *DZI 
        DVDL_cct(3,1) = ( U_cct_ip(3) - U_cct_im(3) ) * DXI * 0.5_WP
        DVDL_cct(3,2) = ( ( U_cct_jp(3) - U_cct(3)    ) * DYCI(JJP) + &
                          ( U_cct(3)    - U_cct_jm(3) ) * DYCI(JJ) ) * 0.5_WP
        
        
        RETURN
    END SUBROUTINE
    
!##################################################################################################
    SUBROUTINE PP_DRIVEN_FORCE
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        USE THERMAL_INFO
        IMPLICIT NONE
        
        REAL(WP)   :: VisForcWall1, VisForcWall1_WORK
        REAL(WP)   :: VisForcWall2, VisForcWall2_WORK
        REAL(WP)   :: BuoFC, BuoFC_WORK
        INTEGER(4) :: IC, KC, M, JC, JJ
        
        IF(FLOWDRVTP.NE.1) RETURN
        
        
        ! VISCOUS FORCE on the bottom wall
        VisForcWall1 = 0.0_WP
        VisForcWall2 = 0.0_WP
        IF(MYID==0) THEN
            DO IC=1, NCL1_io
                DO KC=1, NCL3
                    IF(thermlflg==1) THEN
                        VisForcWall1 = VisForcWall1 + M_WAL_GV(1,ibotwall)*(Q_io(IC,1,KC,1)-0.0_WP)/(YCC(1)-YND(1))*CVISC*DX*DZ
                    ELSE
                        VisForcWall1 = VisForcWall1 + (Q_io(IC,1,KC,1)-0.0_WP)/(YCC(1)-YND(1))*CVISC*DX*DZ
                    END IF
                END DO
            END DO
            !WRITE(*,*) 'Wall1 Friction Force:', VisForcWall1
        END IF
        
        IF(MYID==NPSLV) THEN
            DO IC=1, NCL1_io
                DO KC=1, NCL3
                    IF(thermlflg==1) THEN
                        VisForcWall2 = VisForcWall2 + &
                                        M_WAL_GV(1,itopwall)*(Q_io(IC,N2DO(MYID),KC,1)-0.0_WP)/(YCC(NCL2)-YND(NND2))*CVISC*DX*DZ
                    ELSE
                        VisForcWall2 = VisForcWall2 + (Q_io(IC,N2DO(MYID),KC,1)-0.0_WP)/(YCC(NCL2)-YND(NND2))*CVISC*DX*DZ
                    END IF
                END DO
            END DO
            !WRITE(*,*) 'Wall2 Friction Force:', VisForcWall2
        END IF
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(VisForcWall1,VisForcWall1_WORK,1,MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VisForcWall2,VisForcWall2_WORK,1,MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        ! Gravity force 
        BuoFC = 0.0_WP
        BuoFC_WORK = 0.0_WP
        IF ( GRAVDIR==1 ) THEN
            DO JC=1, N2DO(MYID)
                JJ = JCL2G(JC)
                DO IC=1, NCL1_io
                    DO KC=1, NCL3
                        BuoFC = BuoFC + DENSITY(IC,JC,KC)*F_A/DYFI(JJ)*DX*DZ
                    END DO
                END DO
            END DO
            !WRITE(*,*) 'Gravity in myid', BuoFC, myid
            CALL MPI_BARRIER(ICOMM,IERROR)
            CALL MPI_ALLREDUCE(BuoFC,BuoFC_WORK,1,MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        END IF
        
        IF(MYID.EQ.0) THEN
            FcDrv_IO=( DABS(VisForcWall1_WORK)+DABS(VisForcWall2_WORK)-BuoFC_WORK)/2.0_WP/DBLE(NCL1_io*NCL3)/DX/DZ
            
            !WRITE(*,*) 'FRC', DABS(VisForcWall1_WORK)/DBLE(NCL1_io*NCL3)/DX/DZ, &
            !                  DABS(VisForcWall2_WORK)/DBLE(NCL1_io*NCL3)/DX/DZ, &
            !                  BuoFC_WORK/DBLE(NCL1_io*NCL3)/DX/DZ, FcDrv_IO
        END IF
        
        CALL MPI_BCAST( FcDrv_IO, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        DO JC=1, N2DO(MYID)
            DO M=1,NDV
                FUxzL_IO(JC,M) = FcDrv_IO * U1xzL_io(JC,M)
            END DO
            FUxzL_IO(JC,NDV+1) = FcDrv_IO
        END DO
        
        
        
        RETURN
    END SUBROUTINE
    
