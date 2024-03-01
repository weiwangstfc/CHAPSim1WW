    SUBROUTINE VISCOUS_ALL_EXPLT_Z_io
        USE init_info
        USE MESH_INFO
        USE FLOW_INFO
        USE THERMAL_INFO
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IM, IP
        INTEGER(4) :: JC, JM, JP, JJ, JJP, JJM
        INTEGER(4) :: KC, KM, KP, KS, KSM, NYI
     
        REAL(WP)    :: VIS31C, VIS31P
        REAL(WP)    :: VIS32C, VIS32P
        REAL(WP)    :: DWDZ, DUDZ, DWDX, DVDZ, DWDY
        REAL(WP)    :: TAU31F, TAU31B, DTAU31DX
        REAL(WP)    :: TAU32F, TAU32B, DTAU32DY
        REAL(WP)    :: TAU33F, TAU33B, DTAU33DZ
        REAL(WP)    :: QR_R2, QT_R2, TAU34ADD, TAU34B, TAU34F
        REAL(WP)    :: UR1, UR2, UR21, UR22, UT21, UT22
        REAL(WP)    :: COE1, COE2, COE3  !,  TESTVAL,TESTVAL1,TESTVAL2,TESTVAL3
        
!        INTEGER(4)  :: imsy
!        REAL(WP)     :: QDX2 
!        REAL(WP)     :: q2s1, q2s2      
!        REAL(WP)     :: q2e, q2w      
!        REAL(WP)     :: d11q1e

        
        
        !=============================CARTISIAN COORINDATES========================================
        IF(ICASE==ICHANL .OR. ICASE==IBOX3P) THEN
            RHS_io = 0.0_WP
            !COE1 = DXI * XND2CL * ZND2CL *CVISC
            COE1 = DXI *CVISC
            COE3 = 2.0_WP*DZI*CVISC
            
            DO JC=1,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                !COE2 = DYFI(JJ) * ZND2CL*CVISC
                COE2 = DYFI(JJ) * CVISC
                DO KC=1,NCL3
                    KM=KMV(KC)
                    KP=KPV(KC)
                    
                    DO IC=1,NCL1_io
                       IP=IPV_io(IC)
                       IM=IMV_io(IC)
                       
                       !================DZ_TAU_33=====================================
                       ! at (i,j,k  )
                       DWDZ   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI 
                       TAU33F = 0.5_wp*( VISCOUSITY0(IC,JC,KC) + VISCOUSITY(IC,JC,KC) ) * ( DWDZ - DivU_io(IC,JC,KC) )
                       
                       ! at (i,j,k-1)
                       DWDZ   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JC,KM,3) ) * DZI 
                       TAU33B = 0.5_wp*( VISCOUSITY0(IC,JC,KM) + VISCOUSITY(IC,JC,KM) )* ( DWDZ - DivU_io(IC,JC,KM) )
                       
                       ! at (i,j,k')
                       DTAU33DZ = (TAU33F-TAU33B)*COE3
                       
                       !================DX_TAU_31=====================================
                       ! at (i'+1,j,k')
                       !VIS31P = ( VISCOUSITY(IP,JC,KC) + VISCOUSITY(IP,JC,KM) ) + &
                       !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) 
                       VIS31P = MU_STG(IP,JC,KC,2)
                       DWDX = ( Q_io(IP,JC,KC,3) - Q_io(IC,JC,KC,3) ) * DXI
                       DUDZ = ( Q_io(IP,JC,KC,1) - Q_io(IP,JC,KM,1) ) * DZI 
                       TAU31F = VIS31P * ( DWDX + DUDZ ) 
                       ! at (i',  j,k')
                       !VIS31C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) + &
                       !         ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IM,JC,KM) ) 
                       VIS31C = MU_STG(IC,JC,KC,2)
                       DWDX = ( Q_io(IC,JC,KC,3) - Q_io(IM,JC,KC,3) ) * DXI
                       DUDZ = ( Q_io(IC,JC,KC,1) - Q_io(IC,JC,KM,1) ) * DZI 
                       TAU31B = VIS31C * ( DWDX + DUDZ ) 
                       
                       ! at (i,   j,k')
                       DTAU31DX = (TAU31F-TAU31B)*COE1
                       
                       !================DY_TAU_32=====================================
                       ! at (i,j'+1,k')
                       !VIS32P = ( VISCOUSITY(IC,JP,KC) + VISCOUSITY(IC,JP,KM) ) * YCL2ND_WFF(JJP) + &
                       !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFB(JJP)
                       VIS32P = MU_STG(IC,JP,KC,3)
                       DWDY   = ( Q_io(IC,JP,KC,3) - Q_io(IC,JC,KC,3) ) * DYCI(JJP)
                       DVDZ   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JP,KM,2) ) * DZI 
                       TAU32F = VIS32P * ( DWDY + DVDZ )
                       
                       ! at (i,j',  k')
                       !VIS32C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFF(JJ) + &
                       !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KM) ) * YCL2ND_WFB(JJ)
                       VIS32C = MU_STG(IC,JC,KC,3)
                       DWDY   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JM,KC,3) ) * DYCI(JJ)
                       DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI 
                       TAU32B = VIS32C * ( DWDY + DVDZ )
                       
                       ! at (i,j,   k')
                       DTAU32DY = (TAU32F-TAU32B)*COE2
                       
                        !================D_TAU_Z DIRECTION=================================
                        RHSLLPHI_io(IC,JC,KC) = RHSLLPHI_io(IC,JC,KC) + DTAU31DX + DTAU32DY + DTAU33DZ
                       
                       !IF(JJ==2 .and. IC==1 .and. KC==1) write(*,'(A,4ES13.5)') 'viscz',&
                        !DTAU31DX ,DTAU32DY,DTAU33DZ,RHS_io(IC,JC,KC)
                    
                    END DO
                END DO
            END DO
        
        END IF
     
        !=============================CYLINDRICAL COORINDATES========================================
        IF(ICASE==IPIPEC .OR. ICASE==IANNUL) THEN 
            RHS_io = 0.0_WP
            !COE1 = DXI * XND2CL * ZND2CL *CVISC
            COE1 = DXI *CVISC
            COE3 = 2.0_WP*DZI*CVISC
            
            !=============================Main field========================================
            DO JC=1,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                !COE2 = DYFI(JJ) * ZND2CL*CVISC
                COE2 = DYFI(JJ) *CVISC
                DO KC=1,NCL3
                    KM=KMV(KC)
                    KP=KPV(KC)
                    KS=KSYM(KC)
                    KSM=KSYM(KM)
                    DO IC=1,NCL1_io
                        IP=IPV_io(IC)
                        IM=IMV_io(IC)
                       
                        !================DZ_TAU_33=====================================
                        ! at (i,j,k  )
                        DWDZ   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI    * RCCI2(JJ)
                        IF(ICASE==IPIPEC .and. JJ==1) THEN
                            Ur1    = ( Q_IO(IC,JP,KC,2) - Q_IO(IC,JP,KS,2) ) * 0.50_WP*RNDI1(JJP)
                            QR_R2  = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) +Ur1 ) * YND2CL * RCCI1(JJ)
                        ELSE
                            QR_R2  = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) + Q_io(IC,JC,KC,2)*RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
                        END IF
                        TAU33F = 0.5_wp*( VISCOUSITY0(IC,JC,KC) + VISCOUSITY(IC,JC,KC) ) * ( DWDZ - DivU_io(IC,JC,KC) + QR_R2)
                       
                        ! at (i,j,k-1)
                        DWDZ   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JC,KM,3) ) * DZI    * RCCI2(JJ)
                        IF(ICASE==IPIPEC .and. JJ==1) THEN
                            Ur1    = ( Q_IO(IC,JP,KM,2) - Q_IO(IC,JP,KSM,2) )*0.50_WP*RNDI1(JJP)
                            QR_R2  = ( Q_io(IC,JP,KM,2)*RNDI1(JJP) + Ur1 ) * YND2CL * RCCI1(JJ)
                        ELSE
                            QR_R2  = ( Q_io(IC,JP,KM,2)*RNDI1(JJP) + Q_io(IC,JC,KM,2)*RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
                        END IF
                        TAU33B = 0.5_wp*( VISCOUSITY0(IC,JC,KM) + VISCOUSITY(IC,JC,KM) )* ( DWDZ - DivU_io(IC,JC,KM) + QR_R2)
                       
                        ! at (i,j,k')
                        DTAU33DZ = (TAU33F-TAU33B)*COE3
                       
                        !================DX_TAU_31=====================================
                        ! at (i'+1,j,k')
                        !VIS31P = ( VISCOUSITY(IP,JC,KC) + VISCOUSITY(IP,JC,KM) ) + &
                        !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) 
                        VIS31P = MU_STG(IP,JC,KC,2)
                        DWDX = ( Q_io(IP,JC,KC,3) - Q_io(IC,JC,KC,3) ) * DXI
                        DUDZ = ( Q_io(IP,JC,KC,1) - Q_io(IP,JC,KM,1) ) * DZI 
                        TAU31F = VIS31P * ( DWDX + DUDZ ) 
                        ! at (i',  j,k')
                        !VIS31C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) + &
                        !         ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IM,JC,KM) ) 
                        VIS31C = MU_STG(IC,JC,KC,2)
                        DWDX = ( Q_io(IC,JC,KC,3) - Q_io(IM,JC,KC,3) ) * DXI
                        DUDZ = ( Q_io(IC,JC,KC,1) - Q_io(IC,JC,KM,1) ) * DZI 
                        TAU31B = VIS31C * ( DWDX + DUDZ ) 
                       
                        ! at (i,   j,k')
                        DTAU31DX = (TAU31F-TAU31B)*COE1
                        
                        !================DY_TAU_32=====================================
                        ! at (i,j'+1,k')
                        !VIS32P = ( VISCOUSITY(IC,JP,KC) + VISCOUSITY(IC,JP,KM) ) * YCL2ND_WFF(JJP) + &
                        !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFB(JJP)
                        VIS32P = MU_STG(IC,JP,KC,3)
                        DWDY   = ( Q_io(IC,JP,KC,3)*RCCI1(JJP) - Q_io(IC,JC,KC,3)*RCCI1(JJ)  ) * DYCI(JJP) / RNDI1(JJP)
                        DVDZ   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JP,KM,2) ) * DZI * RNDI1(JJP) 
                        QT_R2  =   Q_io(IC,JP,KC,3)*RCCI1(JJP) *YCL2ND_WFF(JJP) + &
                                   Q_io(IC,JC,KC,3)*RCCI1(JJ)  *YCL2ND_WFB(JJP)
                        TAU32F = VIS32P * ( DWDY + DVDZ - QT_R2)
                        TAU34F = VIS32P * ( DWDY + DVDZ - QT_R2) * RNDI1(JJP)
                       
                        ! at (i,j',  k')
                        !VIS32C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KM) ) * YCL2ND_WFB(JJ)
                        VIS32C = MU_STG(IC,JC,KC,3)
                        IF(ICASE==IPIPEC .and. JJ==1) THEN
                            DWDY   = 0.0_WP
                            UR1    = ( Q_IO(IC,JP,KC,2) - Q_IO(IC,JP,KS, 2) )*0.50_WP*RNDI1(JJP)
                            UR2    = ( Q_IO(IC,JP,KM,2) - Q_IO(IC,JP,KSM,2) )*0.50_WP*RNDI1(JJP)
                            DVDZ   = (UR1 - UR2) * DZI
                        ELSE
                            DWDY   = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) - Q_io(IC,JM,KC,3)*RCCI1(JJM)  ) * DYCI(JJ) / RNDI1(JJ)
                            DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI * RNDI1(JJ) 
                        END IF
                        QT_R2  =   Q_io(IC,JC,KC,3)*RCCI1(JJ)  *YCL2ND_WFF(JJ) + &
                                   Q_io(IC,JM,KC,3)*RCCI1(JJM) *YCL2ND_WFB(JJ)
                        TAU32B = VIS32C * ( DWDY + DVDZ - QT_R2)
                        
                        ! at (i,j,   k')
                        DTAU32DY = (TAU32F-TAU32B)*COE2  !DYFI(JJ) * ZND2CL*CVISC

                        !================DY_TAU_34=====================================
                        IF(ICASE==IPIPEC .and. JJ==1) THEN
                            DWDY   = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) - Q_io(IC,JM,KC,3)*RCCI1(JJM)  ) * DYCI(JJ)
                            UR21   = UR1*RNDI1(JJP)
                            UR22   = UR2*RNDI1(JJP)
                            DVDZ   = (UR21-UR22) *DZI
                        
                            Ut21   = ( Q_io(IC,JP,KC,3)*RCCI1(JJP)   *YCL2ND_WFF(JJP) +  &
                                       Q_io(IC,JC,KC,3)*RCCI1(JJ)    *YCL2ND_WFB(JJP) )
                            Ut22   = ( Q_io(IC,JP,KS,3)*RCCI1(JJP)   *YCL2ND_WFF(JJP) +  &
                                       Q_io(IC,JC,KS,3)*RCCI1(JJ)    *YCL2ND_WFB(JJP) )
                            QT_R2  = (Ut21 - Ut22)*0.50_WP*RNDI1(JJP)
                            TAU34B = VIS32C * ( DWDY + DVDZ - QT_R2)
                        ELSE
                            TAU34B = VIS32C * ( DWDY + DVDZ - QT_R2) * RNDI1(JJ) 
                        END IF
                       
                        
                       
                        TAU34ADD = (TAU34F + TAU34B)*YND2CL *CVISC
                       
                        !================D_TAU_Z DIRECTION=================================
                        RHSLLPHI_io(IC,JC,KC) = RHSLLPHI_io(IC,JC,KC) + DTAU31DX + DTAU32DY + DTAU33DZ + TAU34ADD
                    END DO
                END DO
            END DO
        
        END IF
        
        !CALL DEBUG_WRT_LOCAL(RHSLLPHI_io,1,N2DO(MYID),'visz') ! test
        
!        IDR =3
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
!                                    (Q_IO(IC,JC,KP,IDR)-2.0_WP*Q_IO(IC,JC,KC,IDR)+Q_IO(IC,JC,KM,IDR))*DZQI*RccI2(JJ)+      &
!                                    (Q_IO(IC,JP,KC,IDR)*APVR(JJ,IDR)+             &
!                                     Q_IO(IC,JC,KC,IDR)*ACVR(JJ,IDR)+             &
!                                     Q_IO(IC,JM,KC,IDR)*AMVR(JJ,IDR) )            &
!                                 )
!                        if (jc.eq.1.and.myid.eq.0) then    
!                            imsy = KSYM(kmv(kc))       
!                            q2s1= (q_tg(ic,jp,kc,2)-q_tg(ic,jp,KSYM(kc),2))*0.50_WP*RNDI1(jjp)
!                            q2s2= (q_tg(ic,jp,km,2)-q_tg(ic,jp,imsy,    2))*0.50_WP*RNDI1(jjp)
!                            q2e =  q_tg(ic,jp,kc,2)*RNDI1(jjp) + q2s1
!                            q2w =  q_tg(ic,jp,km,2)*RNDI1(jjp) + q2s2
!                        else
!                            q2e = q_tg(ic,jp,kc,2)*RNDI1(jjp)+q_tg(ic,jc,kc,2)*RNDI1(jj)
!                            q2w = q_tg(ic,jp,km,2)*RNDI1(jjp)+q_tg(ic,jc,km,2)*RNDI1(jj)
!                        end if
!                        d11q1e=(q2e-q2w)*DZI*RCCI1(jj)
!                        RHSL=RHSL+d11q1e/ren  
                        
!                    if(myid==npslv) WRITE(*,'(A,3I3.1,3ES13.5)') 'vis_z',JJ, KC, IC, RHSL, RHS_io(IC,JC,KC), RHSL-RHS_io(IC,JC,KC)
!                ENDDO
!            ENDDO
!        ENDDO
        
        
        
        RETURN
    END SUBROUTINE 
    
    
    SUBROUTINE VISCOUS_PAR_EXPLT_Z_io
        USE init_info
        USE MESH_INFO
        USE FLOW_INFO
        USE THERMAL_INFO
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IM, IP
        INTEGER(4) :: JC, JM, JP, JJ, JJP, JJM
        INTEGER(4) :: KC, KM, KP, KS, KSM, NYI
     
        REAL(WP)    :: VIS31C, VIS31P
        REAL(WP)    :: VIS32C, VIS32P
        REAL(WP)    :: DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ
        REAL(WP)    :: TAU31F_EX, TAU31B_EX, DTAU31DX_EX
        REAL(WP)    :: TAU32F_EX, TAU32B_EX, DTAU32DY_EX
        REAL(WP)    :: TAU33F_EX, TAU33B_EX, DTAU33DZ_EX
        REAL(WP)    :: TAU31F_IM, TAU31B_IM, DTAU31DX_IM
        REAL(WP)    :: TAU32F_IM, TAU32B_IM, DTAU32DY_IM
        REAL(WP)    :: TAU33F_IM, TAU33B_IM, DTAU33DZ_IM
        REAL(WP)    :: TAU31F_AD, TAU31B_AD, DTAU31DX_AD
        REAL(WP)    :: TAU32F_AD, TAU32B_AD, DTAU32DY_AD
        REAL(WP)    :: TAU33F_AD, TAU33B_AD, DTAU33DZ_AD
        REAL(WP)    :: QR_R2, QT_R2, TAU34ADD, TAU34B, TAU34F
        REAL(WP)    :: UR1, UR2, UR21, UR22, UT21, UT22
        REAL(WP)    :: COE1, COE2, COE3  !,  TESTVAL,TESTVAL1,TESTVAL2,TESTVAL3
        
!        INTEGER(4)  :: imsy
!        REAL(WP)     :: QDX2 
!        REAL(WP)     :: q2s1, q2s2      
!        REAL(WP)     :: q2e, q2w      
!        REAL(WP)     :: d11q1e

        
        
        !=============================CARTISIAN COORINDATES========================================
        IF(ICASE==ICHANL) THEN
            RHS_io = 0.0_WP
            !COE1 = DXI * XND2CL * ZND2CL *CVISC
            COE1 = DXI *CVISC
            COE3 = -2.0_WP/3.0_WP*DZI*CVISC
            DO JC=1,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                !COE2 = DYFI(JJ) * ZND2CL*CVISC
                COE2 = DYFI(JJ) * CVISC
                DO KC=1,NCL3
                    KM=KMV(KC)
                    KP=KPV(KC)
                    
                    DO IC=1,NCL1_io
                        IP=IPV_io(IC)
                        IM=IMV_io(IC)
                       
                        !================DZ_TAU_33===(PARTLY)==============================
                        ! at (i,j,k  )
                        DUDX   = ( Q_io(IP,JC,KC,1) - Q_io(IC,JC,KC,1) ) * DXI
                        DVDY   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * DYFI(JJ)
                        DWDZ   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI 
                        TAU33F_EX = VISCOUSITY(IC,JC,KC) * ( DUDX + DVDY )
                        TAU33F_IM = VISCOUSITY(IC,JC,KC) * ( DWDZ )
                       
                        DWDZ = ( G_io(IC,JC,KP,3)*DRHOI_STG(IC,JC,KP,3) - G_io(IC,JC,KC,3)*DRHOI_STG(IC,JC,KC,3) )  * DZI 
                        TAU33F_AD = VISCOUSITY(IC,JC,KC) * ( DWDZ )
                       
                        ! at (i,j,k-1)
                        DUDX   = ( Q_io(IP,JC,KM,1) - Q_io(IC,JC,KM,1) ) * DXI
                        DVDY   = ( Q_io(IC,JP,KM,2) - Q_io(IC,JC,KM,2) ) * DYFI(JJ)
                        DWDZ   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JC,KM,3) ) * DZI 
                        TAU33B_EX = VISCOUSITY(IC,JC,KM) * ( DUDX + DVDY )
                        TAU33B_IM = VISCOUSITY(IC,JC,KM) * ( DWDZ )
                        
                        DWDZ = ( G_io(IC,JC,KC,3)*DRHOI_STG(IC,JC,KC,3) - G_io(IC,JC,KM,3)*DRHOI_STG(IC,JC,KM,3) )  * DZI 
                        TAU33B_AD = VISCOUSITY(IC,JC,KM) * ( DWDZ )
                       
                        ! at (i,j,k')
                        DTAU33DZ_EX = (TAU33F_EX-TAU33B_EX)*COE3
                        DTAU33DZ_IM = (TAU33F_IM-TAU33B_IM)*COE3*(-2.0_WP)
                        
                        DTAU33DZ_AD = (TAU33F_AD-TAU33B_AD)*COE3*(-2.0_WP)
                       
                        !================DX_TAU_31=====================================
                        ! at (i'+1,j,k')
                        !VIS31P = ( VISCOUSITY(IP,JC,KC) + VISCOUSITY(IP,JC,KM) ) + &
                        !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) 
                        VIS31P = MU_STG(IP,JC,KC,2)
                        DWDX = ( Q_io(IP,JC,KC,3) - Q_io(IC,JC,KC,3) ) * DXI
                        DUDZ = ( Q_io(IP,JC,KC,1) - Q_io(IP,JC,KM,1) ) * DZI 
                        TAU31F_EX = VIS31P * ( DUDZ ) 
                        TAU31F_IM = VIS31P * ( DWDX ) 
                       
                        DWDX = ( G_io(IP,JC,KC,3)*DRHOI_STG(IP,JC,KC,3) - G_io(IC,JC,KC,3)*DRHOI_STG(IC,JC,KC,3) ) * DXI
                        TAU31F_AD = VIS31P * ( DWDX ) 
                         
                        ! at (i',  j,k')
                        !VIS31C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) + &
                        !         ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IM,JC,KM) ) 
                        VIS31C = MU_STG(IC,JC,KC,2)
                        DWDX = ( Q_io(IC,JC,KC,3) - Q_io(IM,JC,KC,3) ) * DXI
                        DUDZ = ( Q_io(IC,JC,KC,1) - Q_io(IC,JC,KM,1) ) * DZI 
                        TAU31B_EX = VIS31C * ( DUDZ ) 
                        TAU31B_IM = VIS31C * ( DWDX ) 
                        
                        DWDX = ( G_io(IC,JC,KC,3)*DRHOI_STG(IC,JC,KC,3) - G_io(IM,JC,KC,3)*DRHOI_STG(IM,JC,KC,3) ) * DXI
                        TAU31B_AD = VIS31C * ( DWDX ) 
                       
                        ! at (i,   j,k')
                        DTAU31DX_EX = (TAU31F_EX-TAU31B_EX)*COE1
                        DTAU31DX_IM = (TAU31F_IM-TAU31B_IM)*COE1
                        
                        DTAU31DX_AD = (TAU31F_AD-TAU31B_AD)*COE1
                       
                        !================DY_TAU_32=====================================
                        ! at (i,j'+1,k')
                        !VIS32P = ( VISCOUSITY(IC,JP,KC) + VISCOUSITY(IC,JP,KM) ) * YCL2ND_WFF(JJP) + &
                        !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFB(JJP)
                        VIS32P = MU_STG(IC,JP,KC,3)
                        DWDY   = ( Q_io(IC,JP,KC,3) - Q_io(IC,JC,KC,3) ) * DYCI(JJP)
                        DVDZ   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JP,KM,2) ) * DZI 
                        TAU32F_EX = VIS32P * ( DVDZ )
                        TAU32F_IM = VIS32P * ( DWDY )
                        
                        DWDY = ( G_io(IC,JP,KC,3)*DRHOI_STG(IC,JP,KC,3) - G_io(IC,JC,KC,3)*DRHOI_STG(IC,JC,KC,3) ) * DYCI(JJP)
                        TAU32F_AD= VIS32P * ( DWDY )
                       
                        ! at (i,j',  k')
                        !VIS32C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFF(JJ) + &
                        !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KM) ) * YCL2ND_WFB(JJ)
                        VIS32C = MU_STG(IC,JC,KC,3)
                        DWDY   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JM,KC,3) ) * DYCI(JJ)
                        DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI 
                        TAU32B_EX = VIS32C * ( DVDZ )
                        TAU32B_IM = VIS32C * ( DWDY )
                        
                        DWDY = ( G_io(IC,JC,KC,3)*DRHOI_STG(IC,JC,KC,3) - G_io(IC,JM,KC,3)*DRHOI_STG(IC,JM,KC,3) ) * DYCI(JJ)
                        TAU32B_AD= VIS32C * ( DWDY )
                       
                        ! at (i,j,   k')
                        DTAU32DY_EX = (TAU32F_EX-TAU32B_EX)*COE2
                        DTAU32DY_IM = (TAU32F_IM-TAU32B_IM)*COE2
                        
                        DTAU32DY_AD = (TAU32F_AD-TAU32B_AD)*COE2
                       
                       
                        !================D_TAU_Z DIRECTION=================================
                        RHSLLPHI_io(IC,JC,KC) = DTAU31DX_EX + DTAU32DY_EX + DTAU33DZ_EX + RHSLLPHI_io(IC,JC,KC)
                        RHS_io     (IC,JC,KC) = DTAU31DX_IM + DTAU32DY_IM + DTAU33DZ_IM
                        DIVU_io    (IC,JC,KC) = DTAU31DX_AD + DTAU32DY_AD + DTAU33DZ_AD
                        !WRITE(*,*) DTAU31DX_AD, DTAU32DY_AD, DTAU33DZ_AD
                       !IF(JJ==2 .and. IC==1 .and. KC==1) write(*,'(A,4ES13.5)') 'viscz',&
                        !DTAU31DX ,DTAU32DY,DTAU33DZ,RHS_io(IC,JC,KC)
                    
                    END DO
                END DO
            END DO
        
        END IF
     
!    BELOW TO DO ...
!        !=============================CYLINDRICAL COORINDATES========================================
!        IF(ICASE==IPIPEC .OR. ICASE==IANNUL) THEN 
!            RHS_io = 0.0_WP
!            COE1 = DXI * XND2CL * ZND2CL *CVISC
!            COE3 = 2.0_WP*DZI*CVISC
            
!            !=============================Main field========================================
!            DO JC=1,N2DO(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                JJM=JGMV(JJ)
!                JJP=JGPV(JJ)
!                COE2 = DYFI(JJ) * ZND2CL*CVISC
!                DO KC=1,NCL3
!                    KM=KMV(KC)
!                    KP=KPV(KC)
!                    KS=KSYM(KC)
!                    KSM=KSYM(KM)
!                    DO IC=1,NCL1_io
!                        IP=IPV_io(IC)
!                        IM=IMV_io(IC)
                       
!                        !================DZ_TAU_33=====================================
!                        ! at (i,j,k  )
!                        DWDZ   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI    * RCCI2(JJ)
!                        IF(ICASE==IPIPEC .and. JJ==1) THEN
!                            Ur1    = ( Q_IO(IC,JP,KC,2) - Q_IO(IC,JP,KS,2) ) * 0.50_WP*RNDI1(JJP)
!                            QR_R2  = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) +Ur1 ) * YND2CL * RCCI1(JJ)
!                        ELSE
!                            QR_R2  = ( Q_io(IC,JP,KC,2)*RNDI1(JJP) + Q_io(IC,JC,KC,2)*RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
!                        END IF
!                        TAU33F = VISCOUSITY(IC,JC,KC) * ( DWDZ - DivU_io(IC,JC,KC) + QR_R2)
                       
!                        ! at (i,j,k-1)
!                        DWDZ   = ( Q_io(IC,JC,KC,3) - Q_io(IC,JC,KM,3) ) * DZI    * RCCI2(JJ)
!                        IF(ICASE==IPIPEC .and. JJ==1) THEN
!                            Ur1    = ( Q_IO(IC,JP,KM,2) - Q_IO(IC,JP,KSM,2) )*0.50_WP*RNDI1(JJP)
!                            QR_R2  = ( Q_io(IC,JP,KM,2)*RNDI1(JJP) + Ur1 ) * YND2CL * RCCI1(JJ)
!                        ELSE
!                            QR_R2  = ( Q_io(IC,JP,KM,2)*RNDI1(JJP) + Q_io(IC,JC,KM,2)*RNDI1(JJ) ) * YND2CL * RCCI1(JJ)
!                        END IF
!                        TAU33B = VISCOUSITY(IC,JC,KM) * ( DWDZ - DivU_io(IC,JC,KM) + QR_R2)
                       
!                        ! at (i,j,k')
!                        DTAU33DZ = (TAU33F-TAU33B)*COE3
                       
!                        !================DX_TAU_31=====================================
!                        ! at (i'+1,j,k')
!                        VIS31P = ( VISCOUSITY(IP,JC,KC) + VISCOUSITY(IP,JC,KM) ) + &
!                                 ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) 
!                        DWDX = ( Q_io(IP,JC,KC,3) - Q_io(IC,JC,KC,3) ) * DXI
!                        DUDZ = ( Q_io(IP,JC,KC,1) - Q_io(IP,JC,KM,1) ) * DZI 
!                        TAU31F = VIS31P * ( DWDX + DUDZ ) 
!                        ! at (i',  j,k')
!                        VIS31C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) + &
!                                 ( VISCOUSITY(IM,JC,KC) + VISCOUSITY(IM,JC,KM) ) 
!                        DWDX = ( Q_io(IC,JC,KC,3) - Q_io(IM,JC,KC,3) ) * DXI
!                        DUDZ = ( Q_io(IC,JC,KC,1) - Q_io(IC,JC,KM,1) ) * DZI 
!                        TAU31B = VIS31C * ( DWDX + DUDZ ) 
                       
!                        ! at (i,   j,k')
!                        DTAU31DX = (TAU31F-TAU31B)*COE1
                        
!                        !================DY_TAU_32=====================================
!                        ! at (i,j'+1,k')
!                        VIS32P = ( VISCOUSITY(IC,JP,KC) + VISCOUSITY(IC,JP,KM) ) * YCL2ND_WFF(JJP) + &
!                                 ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFB(JJP)
!                        DWDY   = ( Q_io(IC,JP,KC,3)*RCCI1(JJP) - Q_io(IC,JC,KC,3)*RCCI1(JJ)  ) * DYCI(JJP) / RNDI1(JJP)
!                        DVDZ   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JP,KM,2) ) * DZI * RNDI1(JJP) 
!                        QT_R2  =   Q_io(IC,JP,KC,3)*RCCI1(JJP) *YCL2ND_WFF(JJP) + &
!                                   Q_io(IC,JC,KC,3)*RCCI1(JJ)  *YCL2ND_WFB(JJP)
!                        TAU32F = VIS32P * ( DWDY + DVDZ - QT_R2)
!                        TAU34F = VIS32P * ( DWDY + DVDZ - QT_R2) * RNDI1(JJP)
                       
!                        ! at (i,j',  k')
!                        VIS32C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IC,JC,KM) ) * YCL2ND_WFF(JJ) + &
!                                 ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IC,JM,KM) ) * YCL2ND_WFB(JJ)
!                        IF(ICASE==IPIPEC .and. JJ==1) THEN
!                            DWDY   = 0.0_WP
!                            UR1    = ( Q_IO(IC,JP,KC,2) - Q_IO(IC,JP,KS, 2) )*0.50_WP*RNDI1(JJP)
!                            UR2    = ( Q_IO(IC,JP,KM,2) - Q_IO(IC,JP,KSM,2) )*0.50_WP*RNDI1(JJP)
!                            DVDZ   = (UR1 - UR2) * DZI
!                        ELSE
!                            DWDY   = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) - Q_io(IC,JM,KC,3)*RCCI1(JJM)  ) * DYCI(JJ) / RNDI1(JJ)
!                            DVDZ   = ( Q_io(IC,JC,KC,2) - Q_io(IC,JC,KM,2) ) * DZI * RNDI1(JJ) 
!                        END IF
!                        QT_R2  =   Q_io(IC,JC,KC,3)*RCCI1(JJ)  *YCL2ND_WFF(JJ) + &
!                                   Q_io(IC,JM,KC,3)*RCCI1(JJM) *YCL2ND_WFB(JJ)
!                        TAU32B = VIS32C * ( DWDY + DVDZ - QT_R2)
                        
!                        ! at (i,j,   k')
!                        DTAU32DY = (TAU32F-TAU32B)*COE2  !DYFI(JJ) * ZND2CL*CVISC

!                        !================DY_TAU_34=====================================
!                        IF(ICASE==IPIPEC .and. JJ==1) THEN
!                            DWDY   = ( Q_io(IC,JC,KC,3)*RCCI1(JJ) - Q_io(IC,JM,KC,3)*RCCI1(JJM)  ) * DYCI(JJ)
!                            UR21   = UR1*RNDI1(JJP)
!                            UR22   = UR2*RNDI1(JJP)
!                            DVDZ   = (UR21-UR22) *DZI
                        
!                            Ut21   = ( Q_io(IC,JP,KC,3)*RCCI1(JJP)   *YCL2ND_WFF(JJP) +  &
!                                       Q_io(IC,JC,KC,3)*RCCI1(JJ)    *YCL2ND_WFB(JJP) )
!                            Ut22   = ( Q_io(IC,JP,KS,3)*RCCI1(JJP)   *YCL2ND_WFF(JJP) +  &
!                                       Q_io(IC,JC,KS,3)*RCCI1(JJ)    *YCL2ND_WFB(JJP) )
!                            QT_R2  = (Ut21 - Ut22)*0.50_WP*RNDI1(JJP)
!                            TAU34B = VIS32C * ( DWDY + DVDZ - QT_R2)
!                        ELSE
!                            TAU34B = VIS32C * ( DWDY + DVDZ - QT_R2) * RNDI1(JJ) 
!                        END IF
                       
                        
                       
!                        TAU34ADD = (TAU34F + TAU34B)*YND2CL * ZND2CL*CVISC
                       
!                        !================D_TAU_Z DIRECTION=================================
!                        RHS_io(IC,JC,KC) = DTAU31DX + DTAU32DY + DTAU33DZ + TAU34ADD
!                    END DO
!                END DO
!            END DO
        
!        END IF
        
!        IDR =3
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
!                                    (Q_IO(IC,JC,KP,IDR)-2.0_WP*Q_IO(IC,JC,KC,IDR)+Q_IO(IC,JC,KM,IDR))*DZQI*RccI2(JJ)+      &
!                                    (Q_IO(IC,JP,KC,IDR)*APVR(JJ,IDR)+             &
!                                     Q_IO(IC,JC,KC,IDR)*ACVR(JJ,IDR)+             &
!                                     Q_IO(IC,JM,KC,IDR)*AMVR(JJ,IDR) )            &
!                                 )
!                        if (jc.eq.1.and.myid.eq.0) then    
!                            imsy = KSYM(kmv(kc))       
!                            q2s1= (q_tg(ic,jp,kc,2)-q_tg(ic,jp,KSYM(kc),2))*0.50_WP*RNDI1(jjp)
!                            q2s2= (q_tg(ic,jp,km,2)-q_tg(ic,jp,imsy,    2))*0.50_WP*RNDI1(jjp)
!                            q2e =  q_tg(ic,jp,kc,2)*RNDI1(jjp) + q2s1
!                            q2w =  q_tg(ic,jp,km,2)*RNDI1(jjp) + q2s2
!                        else
!                            q2e = q_tg(ic,jp,kc,2)*RNDI1(jjp)+q_tg(ic,jc,kc,2)*RNDI1(jj)
!                            q2w = q_tg(ic,jp,km,2)*RNDI1(jjp)+q_tg(ic,jc,km,2)*RNDI1(jj)
!                        end if
!                        d11q1e=(q2e-q2w)*DZI*RCCI1(jj)
!                        RHSL=RHSL+d11q1e/ren  
                        
!                    if(myid==npslv) WRITE(*,'(A,3I3.1,3ES13.5)') 'vis_z',JJ, KC, IC, RHSL, RHS_io(IC,JC,KC), RHSL-RHS_io(IC,JC,KC)
!                ENDDO
!            ENDDO
!        ENDDO
        
        
        
        RETURN
    END SUBROUTINE
               
