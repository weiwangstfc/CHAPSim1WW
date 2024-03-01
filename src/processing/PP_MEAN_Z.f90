    SUBROUTINE PP_MEAN_Z_FLOW_nonXperiodic_IO
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        IMPLICIT NONE
        
        INTEGER(4)  :: I, J, K, IP, JP, KP, IM, JM, KM, JJ, JJP, JJM, M, N, H, L, KS, KSM, KSP, P, L1, L2
        REAL(WP)    :: U_cct(NDV), G_cct(NDV)
        REAL(WP)    :: U_cct_ip(NDV), U_cct_jp(NDV), U_cct_kp(NDV)
        REAL(WP)    :: U_cct_im(NDV), U_cct_jm(NDV), U_cct_km(NDV)
        REAL(WP)    :: DVDL_cct(NDV,NDV)
        REAL(WP)    :: COG1(NDV,NDV), COG2(NDV*NDV,NDV*NDV)
        REAL(WP)    :: COE0, COE1, COE2, COE3
        REAL(WP)    :: rtmp
        
        
        COE0 = 1.0_wp/DBLE(NCL3)
        COE1 = 1.0_wp/DBLE(NCL3)*0.5_WP
        COE2 = 1.0_wp/DBLE(NCL3)*(0.5_wp**2)
        COE3 = 1.0_wp/DBLE(NCL3)*(0.5_wp**3)
        
        COG1(1,1) = 1.0_wp/DBLE(NCL3) * DXI
        COG1(1,2) = 1.0_wp/DBLE(NCL3) * 0.25_WP
        COG1(1,3) = 1.0_wp/DBLE(NCL3) * 0.25_WP * DZI
        
        COG1(2,1) = 1.0_wp/DBLE(NCL3) * 0.25_WP * DXI
        COG1(2,2) = 1.0_wp/DBLE(NCL3)
        COG1(2,3) = 1.0_wp/DBLE(NCL3) * 0.25_WP * DZI
        
        COG1(3,1) = 1.0_wp/DBLE(NCL3) * 0.25_WP * DXI
        COG1(3,2) = 1.0_wp/DBLE(NCL3) * 0.25_WP
        COG1(3,3) = 1.0_wp/DBLE(NCL3) * DZI
        
        
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
            JP=JLPV(J)
            JM=JLMV(J)
            JJ=JCL2G(J)
            JJP=JGPV(JJ)
            JJM=JGMV(JJ)
            DO  I=1,NCL1_io
                IP=IPV_io(I)
                IM=IMV_io(I)
                U1zL_io(I,J,:) = 0.0_WP
                G1zL_io(I,J,:) = 0.0_WP
                UPzL_io (I,J,:) = 0.0_WP
                
                U2zL_IO (I,J,:) = 0.0_WP
                UGzL_IO (I,J,:) = 0.0_WP
                UGUzL_IO(I,J,:) = 0.0_WP
                
                DVDL1zL_io(I,J,:,:) = 0.0_WP
                DVDLPzL_io(I,J,:,:) = 0.0_WP
                DVDL2zL_io(I,J,:,:) = 0.0_WP
                
                DO K=1,NCL3
                    KP = KPV(K)
                    KM = KMV(K)
                    KS = KSYM(K)
                    KSP= KSYM(KP)
                    KSM= KSYM(KM)
                    U_cct(1)   =  Q_io(I,J,K,1) + Q_io(IP,J,K,1)              ! U at i, j, k 
                    G_cct(1)   =  G_io(I,J,K,1) + G_io(IP,J,K,1)             ! RHO*U at i, j, k 
                    U_cct_jp(1) = Q_io(I,JP,K,1) + Q_io(IP,JP,K,1)  ! U at i, j+1, k 
                    U_cct_jm(1) = Q_io(I,JM,K,1) + Q_io(IP,JM,K,1)  ! U at i, j-1, k 
                    U_cct_kp(1) = Q_io(I,J,KP,1) + Q_io(IP,J,KP,1)  ! U at i, j, k+1 
                    U_cct_km(1) = Q_io(I,J,KM,1) + Q_io(IP,J,KM,1)  ! U at i, j, k-1
                    
                    
                    U_cct(3)   = (Q_io(I,J,K,3) + Q_io(I,J,KP,3))*RCCI1(JJ)   ! W at i, j, k 
                    G_cct(3)   = (G_io(I,J,K,3) + G_io(I,J,KP,3))*RCCI1(JJ)  ! RHO*U at i, j, k 
                    U_cct_ip(3) = (Q_io(IP,J,K,3) + Q_io(IP,J,KP,3))*RCCI1(JJ)  ! W at i+1, j, k 
                    U_cct_im(3) = (Q_io(IM,J,K,3) + Q_io(IM,J,KP,3))*RCCI1(JJ)  ! W at i-1, j, k 
                    U_cct_jp(3) = (Q_io(I,JP,K,3) + Q_io(I,JP,KP,3))*RCCI1(JJP) ! W at i, j+1, k 
                    U_cct_jm(3) = (Q_io(I,JM,K,3) + Q_io(I,JM,KP,3))*RCCI1(JJM) ! W at i, j-1, k 
                    
                    
                    IF(JJ==1 .AND. ICASE==IPIPEC) THEN
                        U_cct(2)    = (Q_io(I,JP,K,2) - Q_io(I,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_io(I,JP,K,2)*RNDI1(JJP)  ! V at i, j, k 
                        G_cct(2)    = (G_io(I,JP,K,2) - G_io(I,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                       G_io(I,JP,K,2)*RNDI1(JJP)  ! V at i, j, k 
                        U_cct_ip(2) = (Q_io(IP,JP,K,2) - Q_io(IP,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_io(IP,JP,K,2)*RNDI1(JJP)  ! V at i+1, j, k 
                        U_cct_im(2) = (Q_io(IM,JP,K,2) - Q_io(IM,JP,KS,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_io(IM,JP,K,2)*RNDI1(JJP)  ! V at i-1, j, k 
                        U_cct_kp(2) = (Q_io(I,JP,KP,2) - Q_io(I,JP,KSP,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_io(I,JP,KP,2)*RNDI1(JJP)  ! V at i, j, k+1 
                        U_cct_km(2) = (Q_io(I,JP,KM,2) - Q_io(I,JP,KSM,2))*0.50_WP*RNDI1(JJP) + &
                                       Q_io(I,JP,KM,2)*RNDI1(JJP)  ! V at i, j, k-1
                    ELSE
                        U_cct(2)   = (Q_io(I,J,K,2) + Q_io(I,JP,K,2))*RCCI1(JJ)   ! V at i, j, k 
                        G_cct(2)   = (G_io(I,J,K,2) + G_io(I,JP,K,2))*RCCI1(JJ)  ! RHO*U at i, j, k 
                        U_cct_ip(2) = (Q_io(IP,J,K,2) + Q_io(IP,JP,K,2))*RCCI1(JJ)  ! V at i+1, j, k
                        U_cct_im(2) = (Q_io(IM,J,K,2) + Q_io(IM,JP,K,2))*RCCI1(JJ)  ! V at i-1, j, k
                        U_cct_kp(2) = (Q_io(I,J,KP,2) + Q_io(I,JP,KP,2))*RCCI1(JJ)  ! V at i, j, k+1
                        U_cct_km(2) = (Q_io(I,J,KM,2) + Q_io(I,JP,KM,2))*RCCI1(JJ)  ! V at i, j, k-1
                    END IF
                    
                    
                    !================U,V,W,P, G1,G2,G3========================
                    DO M = 1, NDV
                        U1zL_io(I,J,M) = U1zL_io(I,J,M) + U_cct(M)
                        G1zL_io(I,J,M) = G1zL_io(I,J,M) + G_cct(M)
                        UPzL_io(I,J,M) = UPzL_io(I,J,M) + U_cct(M) * PR_io(I,J,K)
                    END DO
                    U1zL_io(I,J,NDV+1) = U1zL_io(I,J,NDV+1) + PR_io(I,J,K)
                    
                    !====U*FLUX===============================================
                    DO M=1,NDV
                        DO N=1,NDV
                            IF(M.GT.N) CYCLE
                            L = (M*(7-M))/2+N-3  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER TRI-MATRIX
                            UGzL_IO(I,J,L) = UGzL_IO(I,J,L) +  U_cct(M) * G_cct(N)
                            U2zL_IO(I,J,L) = U2zL_IO(I,J,L) +  U_cct(M) * U_cct(N)
                        END DO
                    END DO
                    
                    !====U*FLUX*U============================================
                    DO M=1,NDV
                        DO N=1,NDV
                            IF(M.GT.N) CYCLE
                            DO H=1,NDV
                               IF(N.GT.H) CYCLE
                                L = M*(6-M)+(N*(7-N))/2+H-8  ! MATRIX INDEX TO STORAGE LINEAR INDEX FOR UPPER TRI-MATRIX
                                UGUzL_IO(I,J,L) = UGUzL_IO(I,J,L) +  U_cct(M) * G_cct(N) * U_cct(H)
                            END DO
                        END DO
                    END DO
                    
                    
                    !========================DU/DX.DY.DZ======================================
                    DVDL_cct(1,1) =   Q_io(IP,J,K,1) - Q_io(I,J,K,1)    ! *DXI
                    DVDL_cct(1,2) = ( U_cct_jp(1) - U_cct(1) ) * DYCI(JJP) + &
                                    ( U_cct(1) - U_cct_jm(1) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                    DVDL_cct(1,3) =   U_cct_kp(1) - U_cct_km(1)           ! * 0.5_WP * 0.5_WP * DZI
                    
                    !========================DV/DX.DY.DZ====================================== 
                    DVDL_cct(2,1) =   U_cct_ip(2) - U_cct_im(2)           ! * 0.5_WP * 0.5_WP * DXI
                    IF(JJ==1 .AND. ICASE.EQ.IPIPEC) THEN
                        DVDL_cct(2,2) = (   Q_io(I,JP,K,2)*RNDI1(JJP) - &
                                          ( Q_io(I,JP,K,2) - Q_io(I,JP,KS,2) ) * 0.50_WP*RNDI1(JJP) ) * DYFI(JJ)
                    ELSE
                        DVDL_cct(2,2) = ( Q_io(I,JP,K,2)*RNDI1(JJP) - Q_io(I,J,K,2)*RNDI1(JJ) ) * DYFI(JJ)
                    END IF
                    DVDL_cct(2,3) =   U_cct_kp(2) - U_cct_km(2)           ! * 0.5_WP * 0.5_WP * DZI

                    !========================DV/DX.DY.DZ====================================== 
                    DVDL_cct(3,1) =   U_cct_ip(3) - U_cct_im(3)           ! * 0.5_WP * 0.5_WP * DXI
                    DVDL_cct(3,2) = ( U_cct_jp(3) - U_cct(3) ) * DYCI(JJP) + &
                                    ( U_cct(3) - U_cct_jm(3) ) * DYCI(JJ) ! * 0.5_WP * 0.5_WP
                    DVDL_cct(3,3) =   Q_io(I,J,KP,3)*RCCI1(JJ) - Q_io(I,J,K,3)*RCCI1(JJ)      ! *DzI     
                    
                    
                    DO M=1,NDV
                        DO N=1,NDV
                            DVDL1zL_io(I,J,M,N) =  DVDL1zL_io(I,J,M,N) + DVDL_cct(M,N)
                            DVDLPzL_io(I,J,M,N) =  DVDLPzL_io(I,J,M,N) + DVDL_cct(M,N) * PR_io(I,J,K)
                        END DO
                    END DO
                    
!                     DO M=1,NDV
!                        DO N=1,NDV
!                            DO H=1,NDV
!                                DO P=1,NDV
!                                    L1=(M-1)*3+H
!                                    L2=(N-1)*3+P
!                                    DVDL2zL_IO(I,J,L1,L2) =  DVDL2zL_IO(I,J,L1,L2) + DVDL_cct(M,H)*DVDL_cct(N,P)
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
                                    DVDL2zL_IO (I,J,L1,L2)=  DVDL2zL_IO (I,J,L1,L2) + DVDL_cct(M,H)*DVDL_cct(N,P)
                                END DO
                            END DO
                        END DO
                    END DO
                    
                ENDDO
                
                U1zL_io (I,J,NDV+1) = U1zL_io(I,J,NDV+1) * COE0
                U1zL_io (I,J,1:NDV) = U1zL_io(I,J,1:NDV) * COE1
                G1zL_io (I,J,1:NDV) = G1zL_io(I,J,1:NDV) * COE1
                UPzL_io (I,J,:) = UPzL_io (I,J,:) * COE1
                
                U2zL_io (I,J,:) = U2zL_io (I,J,:) * COE2
                UGzL_io (I,J,:) = UGzL_io (I,J,:) * COE2
                
                UGUzL_io(I,J,:) = UGUzL_io(I,J,:) * COE3
                
                DO M=1,NDV
                    DO N=1,NDV
                        DVDL1zL_io(I,J,M,N) =  DVDL1zL_io(I,J,M,N) * COG1(M,N)
                        DVDLPzL_io(I,J,M,N) =  DVDLPzL_io(I,J,M,N) * COG1(M,N)
                    END DO
                END DO

!                DO M=1,NDV
!                    DO N=1,NDV
!                        DO H=1,NDV
!                            DO P=1, NDV
!                                L1=(M-1)*3+H
!                                L2=(N-1)*3+P
!                                DVDL2zL_IO(I,J,L1,L2) = DVDL2zL_IO(I,J,L1,L2) * COG2(L1,L2)
!                            END DO
!                        END DO
!                    END DO
!                END DO
                
                DO M=1,NDV
                    DO H=1, NDV
                        L1   = (M-1)*NDV+H
                        DO N=1, NDV
                            DO P=1, NDV
                                L2=(N-1)*NDV+P
                                DVDL2zL_IO(I,J,L1,L2) = DVDL2zL_IO(I,J,L1,L2) * COG2(L1,L2)
                            END DO
                        END DO
                    END DO
                END DO
                
                
            ENDDO
        ENDDO
        
        G1RATE_io = 0.0_WP
        G1MAXX_io = 0.0_WP
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            DO I=1,NCL1_io
                G1RATE_io=G1RATE_io + G1zL_io(I,J,1)/DYFI(JJ)/RCCI1(JJ)
                G1MAXX_io=DMAX1(G1MAXX_io,DABS(G1zL_io(I,J,1))) 
            END DO
        ENDDO
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(G1RATE_io, G1RATE_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
        CALL MPI_ALLREDUCE(G1MAXX_io, G1MAXX_WORK_io, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        G1BULK_WORK_io = G1RATE_WORK_io*DBLE(NCL3)/VOLM_IO

        
        RETURN
    END SUBROUTINE 
    
    

!***********************************************************************************
    SUBROUTINE PP_MEAN_Z_THEML_nonXperiodic_IO
        use mesh_info
        use init_info
        use flow_info
        use postprocess_info
        use thermal_info
        IMPLICIT NONE
        
        INTEGER(4)  :: I, J, K, IP, KP, JP, N, JJ
        REAL(WP)    :: COE0, COE1
        REAL(WP)    :: U_cct(NDV), G_cct(NDV)
        
        
        COE0 = 1.0_wp/DBLE(NCL3)
        COE1 = 1.0_wp/DBLE(NCL3)*0.5_WP
        
        DO J=1,N2DO(MYID)
            JP=JLPV(J)
            JJ=JCL2G(J)
            DO  I=1,NCL1_io
                IP=IPV_io(I)
                
                T1zL_io(I,J)  = 0.0_WP
                D1zL_io(I,J)  = 0.0_WP
                H1zL_io(I,J)  = 0.0_WP 
                 
                T2zL_io(I,J)  = 0.0_WP
                D2zL_io(I,J)  = 0.0_WP
                H2zL_io(I,J)  = 0.0_WP
                
                DHzL_io(I,J)  = 0.0_WP
                
                UHzL_io(I,J,:)  = 0.0_wp
                GHzL_io(I,J,:)  = 0.0_wp
                
                DO K=1,NCL3
                    KP = KPV(K)
                    
                    T1zL_io(I,J) = T1zL_io(I,J) + TEMPERATURE(I,J,K)
                    D1zL_io(I,J) = D1zL_io(I,J) + DENSITY(I,J,K)
                    H1zL_io(I,J) = H1zL_io(I,J) + ENTHALPY(I,J,K)
                    
                    T2zL_io(I,J) = T2zL_io(I,J) + TEMPERATURE(I,J,K)*TEMPERATURE(I,J,K)
                    D2zL_io(I,J) = D2zL_io(I,J) + DENSITY(I,J,K)   *DENSITY(I,J,K)
                    H2zL_io(I,J) = H2zL_io(I,J) + ENTHALPY(I,J,K)  *ENTHALPY(I,J,K)
                    
                    DHzL_io(I,J) = DHzL_io(I,J) + ENTHALPY(I,J,K)*DENSITY(I,J,K)
                    
                    U_cct(1)   =  Q_io(I,J,K,1) + Q_io(IP,J,K,1)              ! U at i, j, k 
                    U_cct(2)   = (Q_io(I,J,K,2) + Q_io(I,JP,K,2))*RCCI1(JJ)   ! V at i, j, k 
                    U_cct(3)   = (Q_io(I,J,K,3) + Q_io(I,J,KP,3))*RCCI1(JJ)   ! W at i, j, k 
                    
                    G_cct(1)   =  G_io(I,J,K,1) + G_io(IP,J,K,1)             ! RHO*U at i, j, k 
                    G_cct(2)   = (G_io(I,J,K,2) + G_io(I,JP,K,2))*RCCI1(JJ)  ! RHO*U at i, j, k 
                    G_cct(3)   = (G_io(I,J,K,3) + G_io(I,J,KP,3))*RCCI1(JJ)  ! RHO*U at i, j, k 
                    
                    
                    DO N=1,NDV
                        UHzL_io(I,J,N) = UHzL_io(I,J,N) + ENTHALPY(I,J,K)*U_cct(N)
                        GHzL_io(I,J,N) = GHzL_io(I,J,N) + ENTHALPY(I,J,K)*G_cct(N)
                    END DO
                    
                    
                END DO
                
                T1zL_io(I,J)  = T1zL_io(I,J) * COE0
                D1zL_io(I,J)  = D1zL_io(I,J) * COE0
                H1zL_io(I,J)  = H1zL_io(I,J) * COE0
                
                T2zL_io(I,J)  = T2zL_io(I,J) * COE0
                D2zL_io(I,J)  = D2zL_io(I,J) * COE0
                H2zL_io(I,J)  = H2zL_io(I,J) * COE0
                
                DHzL_io(I,J)  = DHzL_io(I,J) * COE0
                
                UHzL_io(I,J,:)  = UHzL_io(I,J,:) * COE1
                GHzL_io(I,J,:)  = GHzL_io(I,J,:) * COE1
                
            END DO
        END DO
        

        H1RATE_io = 0.0_WP
        T1MAXX_io = 0.0_WP
        DO J=1,N2DO(MYID)  !@
            JJ=JCL2G(J)
            DO I=1,NCL1_io
                H1RATE_io=H1RATE_io + GHzL_io(I,J,1)/DYFI(JJ)/RCCI1(JJ)
                T1MAXX_io=DMAX1(T1MAXX_io,DABS(T1zL_io(I,J))) 
            END DO
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

        
        RETURN
    END SUBROUTINE 
    
    

