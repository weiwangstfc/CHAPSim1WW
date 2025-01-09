    SUBROUTINE VISCOUS_ALL_EXPLT_Y_io
        USE init_info
        USE MESH_INFO
        USE FLOW_INFO
        USE THERMAL_INFO
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IM, IP
        INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
        INTEGER(4) :: KC, KM, KP, KS
        INTEGER(4) :: NYI
     
        REAL(WP)    :: VIS21C, VIS21P
        REAL(WP)    :: VIS23C, VIS23P
        REAL(WP)    :: DVDY, DVDX, DUDY, DVDZ, DWDY
        REAL(WP)    :: TAU21F, TAU21B, DTAU21DX
        REAL(WP)    :: TAU22F, TAU22B, DTAU22DY
        REAL(WP)    :: TAU23F, TAU23B, DTAU23DZ
        REAL(WP)    :: DTAU22DD, DWDZ, QR_R2, QT_R2, UR1, TAU24B, TAU24F
        REAL(WP)    :: COE1, COE2, COE3 !, TESTVAL,TESTVAL1,TESTVAL2,TESTVAL3
        
        !REAL(WP)     :: q1e, q1w    , rhsl 
        !REAL(WP)     :: d11q2e
        !integer(4)  :: idr
     
        IF(ICASE==ICHANL .or. ICASE==IBOX3P) THEN
            RHS_io = 0.0_WP
            NYI=1
            IF (MYID.EQ.0 .and. ICASE==ICHANL) NYI=2
     
            !COE1 = DXI*XND2CL*CVISC
            !COE3 = DZI*ZND2CL*CVISC
            
            COE1 = DXI*CVISC
            COE3 = DZI*CVISC
        
            DO JC=NYI,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                COE2 = 2.0_wp*DYCI(JJ)*CVISC
             
                DO KC=1,NCL3
                    KM=KMV(KC)
                    KP=KPV(KC)
                    
                    DO IC=1,NCL1_io
                       IP=IPV_io(IC)
                       IM=IMV_io(IC)
                       
                       !================DY_TAU_22=====================================
                       ! at (i,j,  k)
                       DVDY   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * DYFI(JJ)
                       TAU22F = 0.5_wp*( VISCOUSITY0(IC,JC,KC) + VISCOUSITY(IC,JC,KC) ) * ( DVDY - DivU_io(IC,JC,KC) )
                       ! at (i,j-1,k)
                       DVDY   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JM,KC,2) ) * DYFI(JJM)
                       TAU22B = 0.5_wp*( VISCOUSITY0(IC,JM,KC) + VISCOUSITY(IC,JM,KC) ) * ( DVDY - DivU_io(IC,JM,KC) )
                       ! at (i,j', k)
                       DTAU22DY = (TAU22F-TAU22B)*COE2
                       !================DX_TAU_21=====================================
                       ! at (i'+1,j', k)
                       !VIS21P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IP,JC,KC) ) * YCL2ND_WFF(JJ) + &
                       !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IP,JM,KC) ) * YCL2ND_WFB(JJ)
                       VIS21P = MU_STG(IP,JC,KC,1)
                       DVDX   = ( Q_io(IP,JC,KC,2) - Q_io(IC,JC,KC,2) ) * DXI 
                       DUDY   = ( Q_io(IP,JC,KC,1) - Q_io(IP,JM,KC,1) ) * DYCI(JJ)
                       TAU21F = VIS21P * ( DVDX + DUDY )
                       
                       ! at (i',  j', k)
                       !VIS21C = ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                       !         ( VISCOUSITY(IM,JM,KC) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                       VIS21C = MU_STG(IC,JC,KC,1)
                       DVDX   = ( Q_io(IC,JC,KC,2) - Q_io(IM,JC,KC,2) ) * DXI 
                       DUDY   = ( Q_io(IC,JC,KC,1) - Q_io(IC,JM,KC,1) ) * DYCI(JJ)
                       TAU21B = VIS21C * ( DVDX + DUDY ) 
                       
                       ! at (i,   j', k)
                       DTAU21DX = (TAU21F-TAU21B)*COE1
                       !write(*,*) 'G,j,i,k',jj,kc,ic,VISCOUSITY(IM,JM,KC), VISCOUSITY(IC,JM,KC)
                       
                       !================DZ_TAU_23=====================================
                       ! at (i, j', k'+1)
                       !VIS23P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KP) ) * YCL2ND_WFF(JJ) + &
                       !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KP) ) * YCL2ND_WFB(JJ)
                       VIS23P = MU_STG(IC,JC,KP,3)
                       DVDZ   = ( Q_io(IC,JC,KP,2) - Q_io(IC,JC,KC,2) ) * DZI 
                       DWDY   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JM,KP,3) ) * DYCI(JJ)
                       TAU23F = VIS23P * ( DVDZ + DWDY )
                       
                       ! at (i, j', k'  )
                       !VIS23C = ( VISCOUSITY(IC,JC,KM) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                       !         ( VISCOUSITY(IC,JM,KM) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                       VIS23C = MU_STG(IC,JC,KC,3)
                       DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI 
                       DWDY   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JM,KC,3) ) * DYCI(JJ)
                       TAU23B = VIS23C * ( DVDZ + DWDY )
                       
                       ! at (i, j', k   )
                       DTAU23DZ = (TAU23F-TAU23B)*COE3
                       !================D_TAU_Y DIRECTION=================================
                       DPH_io(IC,JC,KC)= DPH_io(IC,JC,KC) + DTAU21DX + DTAU22DY + DTAU23DZ
                       !if(myid==0 .and. IC == 1 .and. KC == 1 .and. JJ <=4) &
                       ! write(*,*) 'visy-123', DTAU21DX, DTAU22DY, DTAU23DZ, DTAU22DD
!
                    END DO
                END DO
            END DO
        END IF
        
        
        IF(ICASE==IPIPEC .OR. ICASE==IANNUL) THEN 
            RHS_io = 0.0_WP
            NYI=1
            IF (MYID.EQ.0) NYI=2
     
            !COE1 = DXI*XND2CL*CVISC
            !COE3 = DZI*ZND2CL*CVISC
            COE1 = DXI*CVISC
            COE3 = DZI*CVISC

            DO JC=NYI,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                COE2 = 2.0_wp*DYCI(JJ)*CVISC
             
                DO KC=1,NCL3
                    KM=KMV(KC)
                    KP=KPV(KC)
                    KS=KSYM(KC)
                    DO IC=1,NCL1_io
                        IP=IPV_io(IC)
                        IM=IMV_io(IC)
                        
                        !================DY_TAU_22=====================================
                        ! at (i,j,  k)
                        DVDY   = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) - Q_io(IC,JC,KC,2)*RNDI1(JJ) ) * DYFI(JJ)
                        TAU22F = 0.5_wp* ( VISCOUSITY0(IC,JC,KC) + VISCOUSITY(IC,JC,KC) ) * ( DVDY - DivU_io(IC,JC,KC) )/RCCI1(JJ)
                        ! at (i,j-1,k)
                        IF(icase==IPIPEC .and. jj==2) THEN
                            Ur1    = ( Q_IO(IC,JC,KC,2) - Q_IO(IC,JC,KS,2) )*0.50_WP*RNDI1(JJ)
                            DVDY   = ( Q_io(IC,JC,KC,2)*RNDI1(JJ)  - Ur1 ) * DYFI(JJM)
                        ELSE
                            DVDY   = ( Q_io(IC,JC,KC,2)*RNDI1(JJ)  - Q_io(IC,JM,KC,2)*RNDI1(JJM) ) * DYFI(JJM)
                        END IF
                        TAU22B = 0.5_wp*( VISCOUSITY0(IC,JM,KC) + VISCOUSITY(IC,JM,KC) )  * ( DVDY - DivU_io(IC,JM,KC) )/RCCI1(JJM)
                        ! at (i,j', k)
                        DTAU22DY = (TAU22F-TAU22B)*COE2
                        !WRITE(*,'(3I3.1,3ES13.5)') JJ, KC, IC,  TAU22F, TAU22B, DTAU22DY
                       
                        !================DX_TAU_21=====================================
                        ! at (i'+1,j', k)
                        !VIS21P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IP,JC,KC) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IP,JM,KC) ) * YCL2ND_WFB(JJ)
                        VIS21P = MU_STG(IP,JC,KC,1)
                        DVDX   = ( Q_io(IP,JC,KC,2) - Q_io(IC,JC,KC,2) ) * DXI * RNDI1(JJ)
                        DUDY   = ( Q_io(IP,JC,KC,1) - Q_io(IP,JM,KC,1) ) * DYCI(JJ)
                        TAU21F = VIS21P * ( DVDX + DUDY )
                        ! at (i',  j', k)
                        !VIS21C = ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IM,JM,KC) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                        VIS21C = MU_STG(IC,JC,KC,1)
                        DVDX   = ( Q_io(IC,JC,KC,2) - Q_io(IM,JC,KC,2) ) * DXI * RNDI1(JJ)
                        DUDY   = ( Q_io(IC,JC,KC,1) - Q_io(IC,JM,KC,1) ) * DYCI(JJ)
                        TAU21B = VIS21C * ( DVDX + DUDY ) 
                             
                        ! at (i,   j', k)
                        DTAU21DX = (TAU21F-TAU21B)*COE1 / RNDI1(JJ)  
                       
                        !================DZ_TAU_23=====================================
                        ! at (i, j', k'+1)
                        !VIS23P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KP) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KP) ) * YCL2ND_WFB(JJ)
                        VIS23P = MU_STG(IC,JC,KP,3)
                        DVDZ   = ( Q_io(IC,JC,KP,2) - Q_io(IC,JC,KC,2) ) * DZI * RNDI2(JJ)
                        DWDY   = ( Q_io(IC,JC,KP,3)*RCCI1(JJ) - Q_io(IC,JM,KP,3)*RCCI1(JJM) ) * DYCI(JJ)
                        Qt_R2  = ( Q_io(IC,JC,KP,3)*RCCI1(JJ) *YCL2ND_WFF(JJ) + &
                                   Q_io(IC,JM,KP,3)*RCCI1(JJM)*YCL2ND_WFB(JJ) )*RNDI1(JJ) 
                        TAU23F = VIS23P * ( DVDZ + DWDY - Qt_R2)
                       
                        ! at (i, j', k'  )
                        !VIS23C = ( VISCOUSITY(IC,JC,KM) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KM) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                        VIS23C = MU_STG(IC,JC,KC,3)
                        DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI * RNDI2(JJ)
                        DWDY   = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) - Q_io(IC,JM,KC,3)*RCCI1(JJM) ) * DYCI(JJ)
                        Qt_R2  = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) *YCL2ND_WFF(JJ) + &
                                   Q_io(IC,JM,KC,3)*RCCI1(JJM)*YCL2ND_WFB(JJ) )*RNDI1(JJ) 
                        TAU23B = VIS23C * ( DVDZ + DWDY - Qt_R2)
                       
                        ! at (i, j', k   )
                        DTAU23DZ = (TAU23F-TAU23B)*COE3

                        !=================ADDITIONAL======================================
                        ! at (i,j,  k)
                        DWDZ   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI * RCCI2(JJ)
                        Qr_R2  = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) + Q_io(IC,JC,KC,2)*RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
                        TAU24F = 0.5_wp* (VISCOUSITY0(IC,JC,KC)+VISCOUSITY(IC,JC,KC)) * ( DWDZ + Qr_R2 - DivU_io(IC,JC,KC) )
                       
                        ! at (i,j-1, k)
                        DWDZ   = ( Q_io(IC,JM,KP,3) - Q_io(IC,JM,KC,3) ) * DZI * RCCI2(JJM)
                        IF(icase==IPIPEC .and. jj==2) THEN
                            Ur1    = ( Q_IO(IC,JC,KC,2) - Q_IO(IC,JC,KS,2) )*0.50_WP*RNDI1(JJ)
                            Qr_R2  = ( Q_io(IC,JC,KC,2)*RNDI1(JJ) + Ur1 ) * YND2CL * RCCI1(JJM)
                        ELSE
                            Qr_R2  = ( Q_io(IC,JC,KC,2)*RNDI1(JJ) + Q_io(IC,JM,KC,2)*RNDI1(JJM) ) * YND2CL * RCCI1(JJM)
                        END IF
                        TAU24B = 0.5_wp*( VISCOUSITY0(IC,JM,KC) + VISCOUSITY(IC,JM,KC) ) * ( DWDZ + Qr_R2 - DivU_io(IC,JM,KC) )
                       
                        ! at (i,j',  k)
                        DTAU22DD = (TAU24F * YCL2ND_WFF(JJ) + TAU24B * YCL2ND_WFB(JJ))*2.0_wp*CVISC
                       
                        !================D_TAU_Y DIRECTION=================================
                        DPH_io(IC,JC,KC)= DPH_io(IC,JC,KC) + DTAU21DX + DTAU22DY + DTAU23DZ - DTAU22DD

                        if(myid==0 .and. IC == 1 .and. KC == 1 .and. JJ <=4) &
                        write(*,*) 'visy-123', DTAU21DX, DTAU22DY, DTAU23DZ, DTAU22DD
                    END DO
                END DO
            END DO

        END IF
        
        !CALL DEBUG_WRT_LOCAL(DPH_io,1,N2DO(MYID),'visy') ! test
        
!        IDR =2
!        DO KC=1,NCL3
!            KM=KMV(KC)
!            KP=KPV(KC)
!            DO JC=2,N2DO(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                DO IC=1,NCL1_io
!                    IP=IPV_io(IC)
!                    IM=IMV_io(IC)
!                    RHSL  = CVISC*( (Q_IO(IP,JC,KC,IDR)-2.0_WP*Q_IO(IC,JC,KC,IDR)+Q_IO(IM,JC,KC,IDR))*DXQI +              &
!                                    (Q_IO(IC,JC,KP,IDR)-2.0_WP*Q_IO(IC,JC,KC,IDR)+Q_IO(IC,JC,KM,IDR))*DZQI*RNDI2(JJ)+      &
!                                    (Q_IO(IC,JP,KC,IDR)*APVR(JJ,IDR)+             &
!                                     Q_IO(IC,JC,KC,IDR)*ACVR(JJ,IDR)+             &
!                                     Q_IO(IC,JM,KC,IDR)*AMVR(JJ,IDR) )            &
!                                 )
!                    q1e=q_IO(ic,jc,kp,3)*RCCI1(jj)+q_IO(ic,jm,kp,3)*RCCI1(jjm) 
!                    q1w=q_IO(ic,jc,kc,3)*RCCI1(jj)+q_IO(ic,jm,kc,3)*RCCI1(jjm)
!                    d11q2e= -(q1e-q1w)*DZI*RNDI1(jj)
!                    RHSL = d11q2e/ren + RHSL
!                    !if(myid==npslv) WRITE(*,'(a,3I3.1,3ES13.5)') 'vis_y',JJ, KC, IC, RHSL, RHS_io(IC,JC,KC), RHSL-RHS_io(IC,JC,KC)
!                ENDDO
!            ENDDO
!        ENDDO
     
     
        RETURN
    END SUBROUTINE 
    
    
!=================
    SUBROUTINE VISCOUS_PAR_EXPLT_Y_io
        USE init_info
        USE MESH_INFO
        USE FLOW_INFO
        USE THERMAL_INFO
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IM, IP
        INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
        INTEGER(4) :: KC, KM, KP, KS
        INTEGER(4) :: NYI
     
        REAL(WP)    :: VIS21C, VIS21P
        REAL(WP)    :: VIS23C, VIS23P
        REAL(WP)    :: DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ
        REAL(WP)    :: TAU21F_EX, TAU21B_EX, DTAU21DX_EX
        REAL(WP)    :: TAU22F_EX, TAU22B_EX, DTAU22DY_EX
        REAL(WP)    :: TAU23F_EX, TAU23B_EX, DTAU23DZ_EX
        REAL(WP)    :: TAU21F_IM, TAU21B_IM, DTAU21DX_IM
        REAL(WP)    :: TAU22F_IM, TAU22B_IM, DTAU22DY_IM
        REAL(WP)    :: TAU23F_IM, TAU23B_IM, DTAU23DZ_IM
        REAL(WP)    :: TAU21F_AD, TAU21B_AD, DTAU21DX_AD
        REAL(WP)    :: TAU22F_AD, TAU22B_AD, DTAU22DY_AD
        REAL(WP)    :: TAU23F_AD, TAU23B_AD, DTAU23DZ_AD
        REAL(WP)    :: DTAU22DD, QR_R2, QT_R2, UR1, TAU24B, TAU24F
        REAL(WP)    :: COE1, COE2, COE3 !, TESTVAL,TESTVAL1,TESTVAL2,TESTVAL3
        
        !REAL(WP)     :: q1e, q1w    , rhsl 
        !REAL(WP)     :: d11q2e
        !integer(4)  :: idr
     
        IF(ICASE==ICHANL) THEN
            RHS_io = 0.0_WP
            NYI=1
            IF (MYID.EQ.0) NYI=2
     
            !COE1 = DXI*XND2CL*CVISC
            !COE3 = DZI*ZND2CL*CVISC
            
            COE1 = DXI*CVISC
            COE3 = DZI*CVISC
        
            DO JC=NYI,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                COE2 = -2.0_wp/3.0_wp*DYCI(JJ)*CVISC
             
                DO KC=1,NCL3
                    KM=KMV(KC)
                    KP=KPV(KC)
                    
                    DO IC=1,NCL1_io
                        IP=IPV_io(IC)
                        IM=IMV_io(IC)
                       
                        !================DY_TAU_22===(PARTLY)================================
                        ! at (i,j,  k)
                        DUDX = ( Q_io(IP,JC,KC,1) - Q_io(IC,JC,KC,1) ) * DXI
                        DVDY = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * DYFI(JJ)
                        DWDZ = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI
                        TAU22F_EX = VISCOUSITY(IC,JC,KC) * (DUDX+DWDZ) 
                        TAU22F_IM = VISCOUSITY(IC,JC,KC) * (DVDY) 
                        
                        DVDY = ( G_io(IC,JP,KC,2)*DRHOI_STG(IC,JP,KC,2) - G_io(IC,JC,KC,2)*DRHOI_STG(IC,JC,KC,2) ) * DYFI(JJ)
                        TAU22F_AD = VISCOUSITY(IC,JC,KC) * (DVDY) 
                        
                        ! at (i,j-1,  k)
                        DUDX = ( Q_io(IP,JM,KC,1) - Q_io(IC,JM,KC,1) ) * DXI
                        DVDY = ( Q_io(IC,JC,KC,2) - Q_io(IC,JM,KC,2) ) * DYFI(JJM)
                        DWDZ = ( Q_io(IC,JM,KP,3) - Q_io(IC,JM,KC,3) ) * DZI
                        TAU22B_EX = VISCOUSITY(IC,JM,KC) * (DUDX+DWDZ) 
                        TAU22B_IM = VISCOUSITY(IC,JM,KC) * (DVDY) 
                        
                        DVDY = ( G_io(IC,JC,KC,2)*DRHOI_STG(IC,JC,KC,2) - G_io(IC,JM,KC,2)*DRHOI_STG(IC,JM,KC,2) ) * DYFI(JJM)
                        TAU22B_AD = VISCOUSITY(IC,JM,KC) * (DVDY) 
                        
                        ! at (i,j', k)
                        DTAU22DY_EX = (TAU22F_EX-TAU22B_EX)*COE2
                        DTAU22DY_IM = (TAU22F_IM-TAU22B_IM)*COE2*(-2.0_WP)
                        
                        DTAU22DY_AD = (TAU22F_AD-TAU22B_AD)*COE2*(-2.0_WP)
                        
                        !================DX_TAU_21=====================================
                        ! at (i'+1,j', k)
                        !VIS21P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IP,JC,KC) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IP,JM,KC) ) * YCL2ND_WFB(JJ)
                        VIS21P = MU_STG(IP,JC,KC,1)
                        DVDX   = ( Q_io(IP,JC,KC,2) - Q_io(IC,JC,KC,2) ) * DXI 
                        DUDY   = ( Q_io(IP,JC,KC,1) - Q_io(IP,JM,KC,1) ) * DYCI(JJ)
                        TAU21F_EX = VIS21P * ( DUDY )
                        TAU21F_IM = VIS21P * ( DVDX )
                        
                        DVDX = ( G_io(IP,JC,KC,2)*DRHOI_STG(IP,JC,KC,2) - G_io(IC,JC,KC,2)*DRHOI_STG(IC,JC,KC,2) ) * DXI 
                        TAU21F_AD = VIS21P * ( DVDX ) 
                       
                        ! at (i',  j', k)
                        !VIS21C = ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IM,JM,KC) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                        VIS21C = MU_STG(IC,JC,KC,1)
                        DVDX   = ( Q_io(IC,JC,KC,2) - Q_io(IM,JC,KC,2) ) * DXI 
                        DUDY   = ( Q_io(IC,JC,KC,1) - Q_io(IC,JM,KC,1) ) * DYCI(JJ)
                        TAU21B_EX = VIS21C * ( DUDY ) 
                        TAU21B_IM = VIS21C * ( DVDX ) 
                        
                        DVDX = ( G_io(IC,JC,KC,2)*DRHOI_STG(IC,JC,KC,2) - G_io(IM,JC,KC,2)*DRHOI_STG(IM,JC,KC,2) ) * DXI 
                        TAU21B_AD = VIS21C * ( DVDX ) 
                       
                        ! at (i,   j', k)
                        DTAU21DX_EX = (TAU21F_EX-TAU21B_EX)*COE1
                        DTAU21DX_IM = (TAU21F_IM-TAU21B_IM)*COE1
                        
                        DTAU21DX_AD = (TAU21F_AD-TAU21B_AD)*COE1
                       
                        !write(*,*) 'G,j,i,k',jj,kc,ic,VISCOUSITY(IM,JM,KC), VISCOUSITY(IC,JM,KC)
                       
                        !================DZ_TAU_23=====================================
                        ! at (i, j', k'+1)
                        !VIS23P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KP) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KP) ) * YCL2ND_WFB(JJ)
                        VIS23P = MU_STG(IC,JC,KP,3)
                        DVDZ   = ( Q_io(IC,JC,KP,2) - Q_io(IC,JC,KC,2) ) * DZI 
                        DWDY   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JM,KP,3) ) * DYCI(JJ)
                        TAU23F_EX = VIS23P * ( DWDY )
                        TAU23F_IM = VIS23P * ( DVDZ )
                        
                        DVDZ = ( G_io(IC,JC,KP,2)*DRHOI_STG(IC,JC,KP,2) - G_io(IC,JC,KC,2)*DRHOI_STG(IC,JC,KC,2) ) * DZI 
                        TAU23F_AD = VIS23P * ( DVDZ )
                        
                        ! at (i, j', k'  )
                        !VIS23C = ( VISCOUSITY(IC,JC,KM) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KM) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                        VIS23C = MU_STG(IC,JC,KC,3)
                        DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI 
                        DWDY   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JM,KC,3) ) * DYCI(JJ)
                        TAU23B_EX = VIS23C * ( DWDY )
                        TAU23B_IM = VIS23C * ( DVDZ )
                        
                        DVDZ = ( G_io(IC,JC,KC,2)*DRHOI_STG(IC,JC,KC,2) - G_io(IC,JC,KM,2)*DRHOI_STG(IC,JC,KM,2) ) * DZI 
                        TAU23B_AD = VIS23C * ( DVDZ )
                        
                        ! at (i, j', k   )
                        DTAU23DZ_EX = (TAU23F_EX-TAU23B_EX)*COE3
                        DTAU23DZ_IM = (TAU23F_IM-TAU23B_IM)*COE3
                        
                        DTAU23DZ_AD = (TAU23F_AD-TAU23B_AD)*COE3
                       
                        !================D_TAU_Y DIRECTION=================================
                        DPH_io (IC,JC,KC) = DTAU21DX_EX + DTAU22DY_EX + DTAU23DZ_EX + DPH_io(IC,JC,KC)
                        RHS_io (IC,JC,KC) = DTAU21DX_IM + DTAU22DY_IM + DTAU23DZ_IM
                        DIVU_io(IC,JC,KC) = DTAU21DX_AD + DTAU22DY_AD + DTAU23DZ_AD
                        !write(*,*) DTAU21DX_AD, DTAU22DY_AD, DTAU23DZ_AD
                       !IF(JJ==2 .and. IC==1 .and. KC==1) write(*,'(A,4ES13.5)') 'viscy',&
                        !DTAU21DX ,DTAU22DY,DTAU23DZ,RHS_io(IC,JC,KC)
                    END DO
                END DO
            END DO
        END IF
        
! for cylinderical coordinates, to do....
!        IF(ICASE==IPIPEC .OR. ICASE==IANNUL) THEN 
!            RHS_io = 0.0_WP
!            NYI=1
!            IF (MYID.EQ.0) NYI=2
     
!            COE1 = DXI*XND2CL*CVISC
!            COE3 = DZI*ZND2CL*CVISC

!            DO JC=NYI,N2DO(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                JJM=JGMV(JJ)
!                JJP=JGPV(JJ)
!                COE2 = 2.0_wp*DYCI(JJ)*CVISC
             
!                DO KC=1,NCL3
!                    KM=KMV(KC)
!                    KP=KPV(KC)
!                    KS=KSYM(KC)
!                    DO IC=1,NCL1_io
!                        IP=IPV_io(IC)
!                        IM=IMV_io(IC)
                        
!                        !================DY_TAU_22==(PARTLY)================================
!                        ! at (i,j,  k)
!                        DVDY   = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) - Q_io(IC,JC,KC,2)*RNDI1(JJ) ) * DYFI(JJ)
!                        TAU22F = VISCOUSITY(IC,JC,KC) * ( DVDY - DivU_io(IC,JC,KC) )/RCCI1(JJ)
!                        ! at (i,j-1,k)
!                        IF(icase==IPIPEC .and. jj==2) THEN
!                            Ur1    = ( Q_IO(IC,JC,KC,2) - Q_IO(IC,JC,KS,2) )*0.50_WP*RNDI1(JJ)
!                            DVDY   = ( Q_io(IC,JC,KC,2)*RNDI1(JJ)  - Ur1 ) * DYFI(JJM)
!                        ELSE
!                            DVDY   = ( Q_io(IC,JC,KC,2)*RNDI1(JJ)  - Q_io(IC,JM,KC,2)*RNDI1(JJM) ) * DYFI(JJM)
!                        END IF
!                        TAU22B = VISCOUSITY(IC,JM,KC) * ( DVDY - DivU_io(IC,JM,KC) )/RCCI1(JJM)
!                        ! at (i,j', k)
!                        DTAU22DY = (TAU22F-TAU22B)*COE2
!                        !WRITE(*,'(3I3.1,3ES13.5)') JJ, KC, IC,  TAU22F, TAU22B, DTAU22DY
                        
!                        ! at (i,j,  k)
!                        DUDX = ( Q_io(IP,JC,KC,1) - Q_io(IC,JC,KC,1) ) * DXI
!                        DWDZ = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI
!                        TAU22F = VISCOUSITY(IC,JC,KC) * (DUDX+DWDZ) 
                        
!                        ! at (i,j-1,  k)
!                        DUDX = ( Q_io(IP,JM,KC,1) - Q_io(IC,JM,KC,1) ) * DXI
!                        DWDZ = ( Q_io(IC,JM,KP,3) - Q_io(IC,JM,KC,3) ) * DZI
!                        TAU22B = VISCOUSITY(IC,JM,KC) * (DUDX+DWDZ) 
                        
!                        ! at (i,j', k)
!                        DTAU22DY = (TAU22F-TAU22B)*COE2
                       
!                        !================DX_TAU_21=====================================
!                        ! at (i'+1,j', k)
!                        VIS21P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IP,JC,KC) ) * YCL2ND_WFF(JJ) + &
!                                 ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IP,JM,KC) ) * YCL2ND_WFB(JJ)
!                        DVDX   = ( Q_io(IP,JC,KC,2) - Q_io(IC,JC,KC,2) ) * DXI * RNDI1(JJ)
!                        DUDY   = ( Q_io(IP,JC,KC,1) - Q_io(IP,JM,KC,1) ) * DYCI(JJ)
!                        TAU21F = VIS21P * ( DVDX + DUDY )
!                        ! at (i',  j', k)
!                        VIS21C = ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
!                                 ( VISCOUSITY(IM,JM,KC) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
!                        DVDX   = ( Q_io(IC,JC,KC,2) - Q_io(IM,JC,KC,2) ) * DXI * RNDI1(JJ)
!                        DUDY   = ( Q_io(IC,JC,KC,1) - Q_io(IC,JM,KC,1) ) * DYCI(JJ)
!                        TAU21B = VIS21C * ( DVDX + DUDY ) 
                             
!                        ! at (i,   j', k)
!                        DTAU21DX = (TAU21F-TAU21B)*COE1 / RNDI1(JJ)  
                       
!                        !================DZ_TAU_23=====================================
!                        ! at (i, j', k'+1)
!                        VIS23P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KP) ) * YCL2ND_WFF(JJ) + &
!                                 ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KP) ) * YCL2ND_WFB(JJ)
!                        DVDZ   = ( Q_io(IC,JC,KP,2) - Q_io(IC,JC,KC,2) ) * DZI * RNDI2(JJ)
!                        DWDY   = ( Q_io(IC,JC,KP,3)*RCCI1(JJ) - Q_io(IC,JM,KP,3)*RCCI1(JJM) ) * DYCI(JJ)
!                        Qt_R2  = ( Q_io(IC,JC,KP,3)*RCCI1(JJ) *YCL2ND_WFF(JJ) + &
!                                   Q_io(IC,JM,KP,3)*RCCI1(JJM)*YCL2ND_WFB(JJ) )*RNDI1(JJ) 
!                        TAU23F = VIS23P * ( DVDZ + DWDY - Qt_R2)
                       
!                        ! at (i, j', k'  )
!                        VIS23C = ( VISCOUSITY(IC,JC,KM) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
!                                 ( VISCOUSITY(IC,JM,KM) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
!                        DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI * RNDI2(JJ)
!                        DWDY   = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) - Q_io(IC,JM,KC,3)*RCCI1(JJM) ) * DYCI(JJ)
!                        Qt_R2  = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) *YCL2ND_WFF(JJ) + &
!                                   Q_io(IC,JM,KC,3)*RCCI1(JJM)*YCL2ND_WFB(JJ) )*RNDI1(JJ) 
!                        TAU23B = VIS23C * ( DVDZ + DWDY - Qt_R2)
                       
!                        ! at (i, j', k   )
!                        DTAU23DZ = (TAU23F-TAU23B)*COE3

!                        !=================ADDITIONAL======================================
!                        ! at (i,j,  k)
!                        DWDZ   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI * RCCI2(JJ)
!                        Qr_R2  = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) + Q_io(IC,JC,KC,2)*RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
!                        TAU24F = VISCOUSITY(IC,JC,KC) * ( DWDZ + Qr_R2 - DivU_io(IC,JC,KC) )
                       
!                        ! at (i,j-1, k)
!                        DWDZ   = ( Q_io(IC,JM,KP,3) - Q_io(IC,JM,KC,3) ) * DZI * RCCI2(JJM)
!                        IF(icase==IPIPEC .and. jj==2) THEN
!                            Ur1    = ( Q_IO(IC,JC,KC,2) - Q_IO(IC,JC,KS,2) )*0.50_WP*RNDI1(JJ)
!                            Qr_R2  = ( Q_io(IC,JC,KC,2)*RNDI1(JJ) + Ur1 ) * YND2CL * RCCI1(JJM)
!                        ELSE
!                            Qr_R2  = ( Q_io(IC,JC,KC,2)*RNDI1(JJ) + Q_io(IC,JM,KC,2)*RNDI1(JJM) ) * YND2CL * RCCI1(JJM)
!                        END IF
!                        TAU24B = VISCOUSITY(IC,JM,KC) * ( DWDZ + Qr_R2 - DivU_io(IC,JM,KC) )
                       
!                        ! at (i,j',  k)
!                        DTAU22DD = (TAU24F * YCL2ND_WFF(JJ) + TAU24B * YCL2ND_WFB(JJ))*2.0_wp*CVISC
                       
!                        !================D_TAU_Y DIRECTION=================================
!                        RHS_io(IC,JC,KC) = DTAU21DX + DTAU22DY + DTAU23DZ - DTAU22DD
!                        !RHS_io(IC,JC,KC) = DTAU22DY
!                        !WRITE(*,'(3I3.1,5ES13.5)') JJ, KC, IC,  RHS_io(IC,JC,KC), DTAU21DX, DTAU22DY, DTAU23DZ, -DTAU22DD
!                    END DO
!                END DO
!            END DO

!        END IF
        
!!        IDR =2
!!        DO KC=1,NCL3
!!            KM=KMV(KC)
!!            KP=KPV(KC)
!!            DO JC=2,N2DO(MYID)
!!                JM=JLMV(JC)
!!                JP=JLPV(JC)
!!                JJ=JCL2G(JC)
!!                DO IC=1,NCL1_io
!!                    IP=IPV_io(IC)
!!                    IM=IMV_io(IC)
!!                    RHSL  = CVISC*( (Q_IO(IP,JC,KC,IDR)-2.0_WP*Q_IO(IC,JC,KC,IDR)+Q_IO(IM,JC,KC,IDR))*DXQI +              &
!!                                    (Q_IO(IC,JC,KP,IDR)-2.0_WP*Q_IO(IC,JC,KC,IDR)+Q_IO(IC,JC,KM,IDR))*DZQI*RNDI2(JJ)+      &
!!                                    (Q_IO(IC,JP,KC,IDR)*APVR(JJ,IDR)+             &
!!                                     Q_IO(IC,JC,KC,IDR)*ACVR(JJ,IDR)+             &
!!                                     Q_IO(IC,JM,KC,IDR)*AMVR(JJ,IDR) )            &
!!                                 )
!!                    q1e=q_IO(ic,jc,kp,3)*RCCI1(jj)+q_IO(ic,jm,kp,3)*RCCI1(jjm) 
!!                    q1w=q_IO(ic,jc,kc,3)*RCCI1(jj)+q_IO(ic,jm,kc,3)*RCCI1(jjm)
!!                    d11q2e= -(q1e-q1w)*DZI*RNDI1(jj)
!!                    RHSL = d11q2e/ren + RHSL
!!                    !if(myid==npslv) WRITE(*,'(a,3I3.1,3ES13.5)') 'vis_y',JJ, KC, IC, RHSL, RHS_io(IC,JC,KC), RHSL-RHS_io(IC,JC,KC)
!!                ENDDO
!!            ENDDO
!!        ENDDO
     
     
        RETURN
    END SUBROUTINE 
    
    
    
    
    
    SUBROUTINE CHECK_GradP_ON_WALL
        USE init_info
        USE MESH_INFO
        USE FLOW_INFO
        USE THERMAL_INFO
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IM, IP
        INTEGER(4) :: JC, JM, JP, JJ, JJM, JJP
        INTEGER(4) :: KC, KM, KP, KS
        INTEGER(4) :: NYI
     
        REAL(WP)    :: VIS21C, VIS21P
        REAL(WP)    :: VIS23C, VIS23P
        REAL(WP)    :: DVDY, DVDX, DUDY, DVDZ, DWDY
        REAL(WP)    :: TAU21F, TAU21B, DTAU21DX
        REAL(WP)    :: TAU22F, TAU22B, DTAU22DY
        REAL(WP)    :: TAU23F, TAU23B, DTAU23DZ
        REAL(WP)    :: DTAU22DD, DWDZ, QR_R2, QT_R2, UR1, TAU24B, TAU24F, VISYALL
        REAL(WP)    :: COE1, COE2, COE3 !, TESTVAL,TESTVAL1,TESTVAL2,TESTVAL3
        
        !REAL(WP)     :: q1e, q1w    , rhsl 
        !REAL(WP)     :: d11q2e
        !integer(4)  :: idr
     
            IF(MYID.NE.0) RETURN
            
            COE1 = DXI*XND2CL*CVISC
            COE3 = DZI*ZND2CL*CVISC
            
            JC =1
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJM=JGMV(JJ)
            JJP=JGPV(JJ)
            COE2 = 2.0_wp*DYCI(JJ)*CVISC
         
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
                
                DO IC=1,NCL1_io
                   IP=IPV_io(IC)
                   IM=IMV_io(IC)
                   
                   !================DY_TAU_22=====================================
                   ! at (i,j,  k)
                   DVDY   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * DYFI(JJ)
                   TAU22F = VISCOUSITY(IC,JC,KC) * ( DVDY - DivU_io(IC,JC,KC) )
                   ! at (i,j-1,k)
                   DVDY   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JM,KC,2) ) * DYFI(JJM)
                   TAU22B = VISCOUSITY(IC,JM,KC) * ( DVDY - DivU_io(IC,JM,KC) )
                   ! at (i,j', k)
                   DTAU22DY = (TAU22F-TAU22B)*COE2
                   
                   !================DX_TAU_21=====================================
                   ! at (i'+1,j', k)
                   VIS21P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IP,JC,KC) ) * YCL2ND_WFF(JJ) + &
                            ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IP,JM,KC) ) * YCL2ND_WFB(JJ)
                   DVDX   = ( Q_io(IP,JC,KC,2) - Q_io(IC,JC,KC,2) ) * DXI 
                   DUDY   = ( Q_io(IP,JC,KC,1) - Q_io(IP,JM,KC,1) ) * DYCI(JJ)
                   TAU21F = VIS21P * ( DVDX + DUDY )
                   
                   ! at (i',  j', k)
                   VIS21C = ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                            ( VISCOUSITY(IM,JM,KC) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                   DVDX   = ( Q_io(IC,JC,KC,2) - Q_io(IM,JC,KC,2) ) * DXI 
                   DUDY   = ( Q_io(IC,JC,KC,1) - Q_io(IC,JM,KC,1) ) * DYCI(JJ)
                   TAU21B = VIS21C * ( DVDX + DUDY ) 
                   
                   ! at (i,   j', k)
                   DTAU21DX = (TAU21F-TAU21B)*COE1
                   
                   !write(*,*) 'G,j,i,k',jj,kc,ic,VISCOUSITY(IM,JM,KC), VISCOUSITY(IC,JM,KC)
                   
                   !================DZ_TAU_23=====================================
                   ! at (i, j', k'+1)
                   !VIS23P = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KP) ) * YCL2ND_WFF(JJ) + &
                   !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KP) ) * YCL2ND_WFB(JJ)
                   VIS23P = MU_STG(IC,JC,KP,3)
                   DVDZ   = ( Q_io(IC,JC,KP,2) - Q_io(IC,JC,KC,2) ) * DZI 
                   DWDY   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JM,KP,3) ) * DYCI(JJ)
                   TAU23F = VIS23P * ( DVDZ + DWDY )
                   
                   ! at (i, j', k'  )
                   !VIS23C = ( VISCOUSITY(IC,JC,KM) + VISCOUSITY(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                   !         ( VISCOUSITY(IC,JM,KM) + VISCOUSITY(IC,JM,KC) ) * YCL2ND_WFB(JJ)
                   VIS23C = MU_STG(IC,JC,KC,3)
                   DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI 
                   DWDY   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JM,KC,3) ) * DYCI(JJ)
                   TAU23B = VIS23C * ( DVDZ + DWDY )
                   
                   ! at (i, j', k   )
                   DTAU23DZ = (TAU23F-TAU23B)*COE3
                   
                   !================D_TAU_Y DIRECTION=================================
                   VISYALL = DTAU22DY + DTAU21DX + DTAU23DZ
                   WRITE(*,'(A, 2I3.1,4ES18.10)') 'G(P)_WALL', KC, IC, DTAU22DY, DTAU21DX , DTAU23DZ, VISYALL
                END DO
            END DO
 
        
        
     
        RETURN
    END SUBROUTINE 
