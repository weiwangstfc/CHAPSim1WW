    SUBROUTINE VISCOUS_ALL_EXPLT_X_io
        USE MESH_INFO
        USE FLOW_INFO
        USE THERMAL_INFO
        use init_info
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IM, IP, NXI
        INTEGER(4) :: JC, JM, JP, JJ, JJP, JJM
        INTEGER(4) :: KC, KM, KP
     
        REAL(WP)    :: VIS12C, VIS12P
        REAL(WP)    :: VIS13C, VIS13P
        REAL(WP)    :: DUDX, DUDY, DVDX, DUDZ, DWDX
        REAL(WP)    :: TAU11F, TAU11B, DTAU11DX
        REAL(WP)    :: TAU12F, TAU12B, DTAU12DY
        REAL(WP)    :: TAU13F, TAU13B, DTAU13DZ
        REAL(WP)    :: COE1, COE2, COE3,  COE32 !, TESTVAL,TESTVAL1,TESTVAL2,TESTVAL3
     
        RHS_io = 0.0_WP
     
        COE1 = 2.0_WP*DXI*CVISC
        !COE3 = DZI * XND2CL * ZND2CL *CVISC
        COE3 = DZI *CVISC
        
        NXI =  1
        IF(TGFLOWFLG) NXI  =  2 
        
        DO JC=1,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJM=JGMV(JJ)
            JJP=JGPV(JJ)
            !COE2 = DYFI(JJ) * XND2CL*CVISC * RCCI1(JJ)
            COE2 = DYFI(JJ) * CVISC * RCCI1(JJ)
            COE32= COE3 * RCCI2(JJ)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
                
                DO IC=NXI,NCL1_io
                    IP=IPV_io(IC)
                    IM=IMV_io(IC)
                   
                    !================DX_TAU_11=====================================
                    ! at (i,  j,k)
                    DUDX   = ( Q_io(IP,JC,KC,1) - Q_io(IC,JC,KC,1) ) * DXI 
                    TAU11F = 0.5_wp*( VISCOUSITY0(IC,JC,KC) + VISCOUSITY(IC,JC,KC) ) * ( DUDX - DivU_io(IC,JC,KC) )
                   
                    ! at (i-1,j,k)
                    DUDX   = ( Q_io(IC,JC,KC,1) - Q_io(IM,JC,KC,1) ) * DXI
                    TAU11B = 0.5_wp*( VISCOUSITY0(IM,JC,KC) + VISCOUSITY(IM,JC,KC) ) * ( DUDX - DivU_io(IM,JC,KC) )
                   
                    ! at (i', j,k)
                    DTAU11DX = (TAU11F-TAU11B)*COE1
                   
                    !================DY_TAU_12=====================================
                    ! at (i',j'+1,k)
                    !VIS12P = ( VISCOUSITY(IC,JP,KC) + VISCOUSITY(IM,JP,KC) ) * YCL2ND_WFF(JJP) + &
                    !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) * YCL2ND_WFB(JJP)
                    VIS12P = MU_STG(IC,JP,KC,1)
                    DUDY   = ( Q_io(IC,JP,KC,1) - Q_io(IC,JC,KC,1) ) * DYCI(JJP) / RNDI1(JJP)
                    DVDX   = ( Q_io(IC,JP,KC,2) - Q_io(IM,JP,KC,2) ) * DXI 
                    TAU12F = VIS12P * ( DUDY + DVDX )
                   
                    ! at (i',j',k)
                    !VIS12C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IM,JM,KC) ) * YCL2ND_WFB(JJ)
                    VIS12C = MU_STG(IC,JC,KC,1)
                    IF(icase.eq.IPIPEC .and. jj.eq.1) THEN
                        DUDY   = 0.0_wp
                    ELSE
                        DUDY   = ( Q_io(IC,JC,KC,1) - Q_io(IC,JM,KC,1) ) * DYCI(JJ) / RNDI1(JJ)
                    END IF
                    DVDX   = ( Q_io(IC,JC,KC,2) - Q_io(IM,JC,KC,2) ) * DXI 
                    TAU12B = VIS12C * ( DUDY + DVDX )
                   
                    ! at (i', j,k)
                    DTAU12DY = (TAU12F-TAU12B)*COE2
                    !IF(myid==0) write(*,*) 'viscx', JC, KC, IC, &
                    !TAU12F, TAU12B, VIS12P, VIS12C
                   
                    !================DZ_TAU_13=====================================
                    ! at (i',j,k'+1)
                    !VIS13P = ( VISCOUSITY(IC,JC,KP) + VISCOUSITY(IM,JC,KP) ) + &
                    !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) 
                    VIS13P = MU_STG(IC,JC,KP,2)
                    DUDZ = ( Q_io(IC,JC,KP,1) - Q_io(IC,JC,KC,1) ) * DZI
                    DWDX = ( Q_io(IC,JC,KP,3) - Q_io(IM,JC,KP,3) ) * DXI 
                    TAU13F = VIS13P * ( DUDZ + DWDX ) 
                   
                    ! at (i',j,k')
                    !VIS13C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) + &
                    !         ( VISCOUSITY(IC,JC,KM) + VISCOUSITY(IM,JC,KM) )
                    VIS13C = MU_STG(IC,JC,KC,2)
                    DUDZ = ( Q_io(IC,JC,KC,1) - Q_io(IC,JC,KM,1) ) * DZI
                    DWDX = ( Q_io(IC,JC,KC,3) - Q_io(IM,JC,KC,3) ) * DXI 
                    TAU13B = VIS13C * ( DUDZ + DWDX ) 
                   
                    ! at (i', j,k)
                    DTAU13DZ = (TAU13F-TAU13B)*COE32
                   
                    !================D_TAU_X DIRECTION=================================
                    Qtmp_io(IC,JC,KC) = Qtmp_io(IC,JC,KC) + DTAU11DX + DTAU12DY + DTAU13DZ

                    !if(myid==0 .and. IC == 1 .and. KC == 1 .and. JJ <=4) &
                    !write(*,*) 'visx-123', DTAU11DX, DTAU12DY, DTAU13DZ
        
                END DO
            END DO
        END DO
        
        !CALL DEBUG_WRT_LOCAL(Qtmp_io,1,N2DO(MYID),'visx') ! test
        
        RETURN
    END SUBROUTINE 
    
    
!**************************************************************************    
    SUBROUTINE VISCOUS_PAR_EXPLT_X_io
    ! Refer to Nicoud (2000) for the split of viscous terms.
        USE MESH_INFO
        USE FLOW_INFO
        USE THERMAL_INFO
        use init_info
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IM, IP, NXI
        INTEGER(4) :: JC, JM, JP, JJ, JJP, JJM
        INTEGER(4) :: KC, KM, KP
     
        REAL(WP)    :: VIS12C, VIS12P
        REAL(WP)    :: VIS13C, VIS13P
        REAL(WP)    :: DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ
        REAL(WP)    :: TAU11F_EX, TAU11B_EX, DTAU11DX_EX
        REAL(WP)    :: TAU12F_EX, TAU12B_EX, DTAU12DY_EX
        REAL(WP)    :: TAU13F_EX, TAU13B_EX, DTAU13DZ_EX
        REAL(WP)    :: TAU11F_IM, TAU11B_IM, DTAU11DX_IM
        REAL(WP)    :: TAU12F_IM, TAU12B_IM, DTAU12DY_IM
        REAL(WP)    :: TAU13F_IM, TAU13B_IM, DTAU13DZ_IM
        REAL(WP)    :: TAU11F_AD, TAU11B_AD, DTAU11DX_AD
        REAL(WP)    :: TAU12F_AD, TAU12B_AD, DTAU12DY_AD
        REAL(WP)    :: TAU13F_AD, TAU13B_AD, DTAU13DZ_AD
        REAL(WP)    :: COE1, COE2, COE3,  COE32 !, TESTVAL,TESTVAL1,TESTVAL2,TESTVAL3
     
        RHS_io = 0.0_WP
     
        COE1 = -2.0_wp/3.0_wp*DXI*CVISC
        !COE3 = DZI * XND2CL * ZND2CL *CVISC
        COE3 = DZI *CVISC
        
        NXI =  1
        IF(TGFLOWFLG) NXI  =  2 
        
        DO JC=1,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJM=JGMV(JJ)
            JJP=JGPV(JJ)
            !COE2 = DYFI(JJ) * XND2CL*CVISC * RCCI1(JJ)
            COE2 = DYFI(JJ) * CVISC * RCCI1(JJ)
            COE32= COE3 * RCCI2(JJ)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
                
                DO IC=NXI,NCL1_io
                    IP=IPV_io(IC)
                    IM=IMV_io(IC)
                   
                    !================DX_TAU_11=(partly)====================================
                    ! at (i,  j,k)
                    DUDX   = ( Q_io(IP,JC,KC,1) - Q_io(IC,JC,KC,1) ) * DXI 
                    DVDY   = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * DYFI(JJ) * RCCI1(JJ)
                    DWDZ   = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * DZI * RCCI2(JJ)
                    TAU11F_EX = VISCOUSITY(IC,JC,KC) * ( DVDY + DWDZ )
                    TAU11F_IM = VISCOUSITY(IC,JC,KC) * DUDX
                    
                    DUDX = ( G_io(IP,JC,KC,1)*DRHOI_STG(IP,JC,KC,1) - &
                             G_io(IC,JC,KC,1)*DRHOI_STG(IC,JC,KC,1) ) * DXI 
                    TAU11F_AD = VISCOUSITY(IC,JC,KC) * DUDX
                   
                    ! at (i-1,j,k)
                    DUDX   = ( Q_io(IC,JC,KC,1) - Q_io(IM,JC,KC,1) ) * DXI 
                    DVDY   = ( Q_io(IM,JP,KC,2) - Q_io(IM,JC,KC,2) ) * DYFI(JJ) * RCCI1(JJ)
                    DWDZ   = ( Q_io(IM,JC,KP,3) - Q_io(IM,JC,KC,3) ) * DZI * RCCI2(JJ)
                    TAU11B_EX = VISCOUSITY(IM,JC,KC) * ( DVDY + DWDZ )
                    TAU11B_IM = VISCOUSITY(IM,JC,KC) * DUDX
                    
                    DUDX = ( G_io(IC,JC,KC,1)*DRHOI_STG(IC,JC,KC,1) - &
                             G_io(IM,JC,KC,1)*DRHOI_STG(IM,JC,KC,1) ) * DXI 
                    TAU11B_AD = VISCOUSITY(IM,JC,KC) * DUDX
                    
                    ! at (i', j,k)
                    DTAU11DX_EX = (TAU11F_EX-TAU11B_EX)*COE1
                    DTAU11DX_IM = (TAU11F_IM-TAU11B_IM)*COE1*(-2.0_WP)
                    
                    DTAU11DX_AD = (TAU11F_AD-TAU11B_AD)*COE1*(-2.0_WP)
                   
                    !================DY_TAU_12=====================================
                    ! at (i',j'+1,k)
                    !VIS12P = ( VISCOUSITY(IC,JP,KC) + VISCOUSITY(IM,JP,KC) ) * YCL2ND_WFF(JJP) + &
                    !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) * YCL2ND_WFB(JJP)
                    VIS12P = MU_STG(IC,JP,KC,1)
                    DVDX   = ( Q_io(IC,JP,KC,2) - Q_io(IM,JP,KC,2) ) * DXI 
                    DUDY   = ( Q_io(IC,JP,KC,1) - Q_io(IC,JC,KC,1) ) * DYCI(JJP) / RNDI1(JJP)
                    TAU12F_EX = VIS12P * ( DVDX )
                    TAU12F_IM = VIS12P * ( DUDY )
                    
                    DUDY = ( G_io(IC,JP,KC,1)*DRHOI_STG(IC,JP,KC,1) - &
                             G_io(IC,JC,KC,1)*DRHOI_STG(IC,JC,KC,1) ) * DYCI(JJP) / RNDI1(JJP)
                    TAU12F_AD = VIS12P * ( DUDY )
                    ! at (i',j',k)
                    !VIS12C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) * YCL2ND_WFF(JJ) + &
                    !         ( VISCOUSITY(IC,JM,KC) + VISCOUSITY(IM,JM,KC) ) * YCL2ND_WFB(JJ)
                    VIS12C = MU_STG(IC,JC,KC,1)
                    IF(icase.eq.IPIPEC .and. jj.eq.1) THEN
                        DUDY   = 0.0_wp
                    ELSE
                        DUDY   = ( Q_io(IC,JC,KC,1) - Q_io(IC,JM,KC,1) ) * DYCI(JJ) / RNDI1(JJ)
                    END IF
                    DVDX   = ( Q_io(IC,JC,KC,2) - Q_io(IM,JC,KC,2) ) * DXI 
                    TAU12B_EX = VIS12C * ( DVDX )
                    TAU12B_IM = VIS12C * ( DUDY )
                    
                    DUDY = ( G_io(IC,JC,KC,1)*DRHOI_STG(IC,JC,KC,1) - &
                             G_io(IC,JM,KC,1)*DRHOI_STG(IC,JM,KC,1) ) * DYCI(JJ) / RNDI1(JJ)
                    TAU12B_AD = VIS12C * ( DUDY )
                   
                    ! at (i', j,k)
                    DTAU12DY_EX = (TAU12F_EX-TAU12B_EX)*COE2
                    DTAU12DY_IM = (TAU12F_IM-TAU12B_IM)*COE2
                    
                    DTAU12DY_AD = (TAU12F_AD-TAU12B_AD)*COE2
                   
                    !================DZ_TAU_13=====================================
                    ! at (i',j,k'+1)
                    !VIS13P = ( VISCOUSITY(IC,JC,KP) + VISCOUSITY(IM,JC,KP) ) + &
                    !         ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) 
                    VIS13P = MU_STG(IC,JC,KP,2)
                    DUDZ = ( Q_io(IC,JC,KP,1) - Q_io(IC,JC,KC,1) ) * DZI
                    DWDX = ( Q_io(IC,JC,KP,3) - Q_io(IM,JC,KP,3) ) * DXI 
                    TAU13F_EX = VIS13P * ( DWDX ) 
                    TAU13F_IM = VIS13P * ( DUDZ ) 
                    
                    DUDZ = ( Q_io(IC,JC,KP,1)*DRHOI_STG(IC,JC,KP,1) - &
                             Q_io(IC,JC,KC,1)*DRHOI_STG(IC,JC,KC,1) ) * DZI
                    TAU13F_AD = VIS13P * ( DUDZ ) 
                   
                    ! at (i',j,k')
                    !VIS13C = ( VISCOUSITY(IC,JC,KC) + VISCOUSITY(IM,JC,KC) ) + &
                    !         ( VISCOUSITY(IC,JC,KM) + VISCOUSITY(IM,JC,KM) )
                    VIS13C = MU_STG(IC,JC,KC,2)
                    DUDZ = ( Q_io(IC,JC,KC,1) - Q_io(IC,JC,KM,1) ) * DZI
                    DWDX = ( Q_io(IC,JC,KC,3) - Q_io(IM,JC,KC,3) ) * DXI 
                    TAU13B_EX = VIS13C * ( DWDX ) 
                    TAU13B_IM = VIS13C * ( DUDZ ) 
                    
                    DUDZ = ( Q_io(IC,JC,KC,1)*DRHOI_STG(IC,JC,KC,1) - &
                             Q_io(IC,JC,KM,1)*DRHOI_STG(IC,JC,KM,1) ) * DZI
                    TAU13B_AD = VIS13C * ( DUDZ ) 
                   
                    ! at (i', j,k)
                    DTAU13DZ_EX = (TAU13F_EX-TAU13B_EX)*COE32
                    DTAU13DZ_IM = (TAU13F_IM-TAU13B_IM)*COE32
                    
                    DTAU13DZ_AD = (TAU13F_AD-TAU13B_AD)*COE32
                   
                    !================D_TAU_X DIRECTION=================================
                    Qtmp_io(IC,JC,KC) = DTAU11DX_EX + DTAU12DY_EX + DTAU13DZ_EX + Qtmp_io(IC,JC,KC)
                    RHS_io (IC,JC,KC) = DTAU11DX_IM + DTAU12DY_IM + DTAU13DZ_IM
                    DIVU_io(IC,JC,KC) = DTAU11DX_AD + DTAU12DY_AD + DTAU13DZ_AD
                    !write(*,*) DTAU11DX_AD,DTAU12DY_AD,DTAU13DZ_AD
                    !IF(JJ==2 .and. IC==1 .and. KC==1) write(*,'(A,4ES13.5)') 'viscx',&
                    !DTAU11DX ,DTAU12DY,DTAU13DZ,RHS_io(IC,JC,KC)
                END DO
            END DO
        END DO
        
        RETURN
    END SUBROUTINE 
               
               
               
