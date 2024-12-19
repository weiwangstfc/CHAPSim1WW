    SUBROUTINE RHS_ENERGY_EXPLICIT
!>   @note
!>    1) the calculation of heat flux is based on a constant velocity field. 
!>    2) \rho h u is treated as (\rho u)h, not (\rho h)u
!>    3) WALL BC IS NOT INCLUDED, WHICH WILL BE DEALED WITH SEPARATELY.  
   
        use thermal_info
        use mesh_info
        use flow_info
        USE INIT_INFO
        IMPLICIT NONE
     
        INTEGER(4) :: IC, IP, IM
        INTEGER(4) :: KC, KP, KM
        INTEGER(4) :: JC, JP, JM, JJ, JJP, JS,JE
        REAL(WP)     :: DX_RHOHU, DY_RHOHV, DZ_RHOHW 
        REAL(WP)     :: DX_KDT, DY_KDT, DZ_KDT, KDTy1, KDTy2
        REAL(WP)     :: COE1, COE2, COE3, COE11, COE22, COE221, COE222, COE33, COE331
     
        !logical,external :: isnan1, isinf1
     
        RHS_ENERGY = 0.0_wp
        
        COE1 = DXI  * XND2CL * 0.5_wp
        COE11= DXQI * XND2CL * CTHECD 
        DO JC=1,N2DO(MYID) 
            JP = JLPV(JC)
            JM = JLMV(JC)
            JJ=  JCL2G(JC)
            JJP= JGPV(JJ)
            
            !=============coefficients=============================
            !COEFFICIENT FOR CONVECTIONS
            COE2  = DYFI(JJ)     * RCCI1(JJ) * 0.5_wp
            COE3  = DZI * ZND2CL * RCCI2(JJ) * 0.5_wp
            
            !COEFFICIENTS FOR VISCOUS TERMS
            COE33 = DZQI * ZND2CL * RCCI2(JJ) * CTHECD 
            COE22 = DYFI(JJ) * RCCI1(JJ)* CTHECD 
            
            IF(JJ==1) THEN
                COE221= DYCI(JJP) /RNDI1(JJP)
                IF( BCWALLHEAT(Ibotwall)==isoFluxWall    ) THEN
                    IF(ICASE==IPIPEC) THEN
                        COE222= 0.0_wp
                    ELSE
                        COE222= 1.0_wp /RNDI1(JJ)
                    END IF
                END IF
                
                IF( BCWALLHEAT(Ibotwall)==isoThermalWall    ) THEN
                    IF(ICASE==IPIPEC) THEN
                        COE222= 0.0_wp
                    ELSE
                        COE222= DYCI(JJ)/RNDI1(JJ)
                    END IF
                END IF
                
            ELSE IF(JJ==NCL2) THEN
                COE222= DYCI(JJ) /RNDI1(JJ)
                IF( BCWALLHEAT(Itopwall)==isoFluxWall    ) THEN
                    COE221= 1.0_wp /RNDI1(JJP)
                END IF
                IF( BCWALLHEAT(Itopwall)==isoThermalWall    ) THEN
                    COE221= DYCI(JJP) /RNDI1(JJP)
                END IF
            ELSE
                COE221= DYCI(JJP) /RNDI1(JJP)
                COE222= DYCI(JJ)  /RNDI1(JJ)
            END IF
            !=============coefficients=============================
            
            DO IC=1,NCL1_IO
                IP = IPV_io(IC)
                IM = IMV_io(IC)
                DO KC=1,NCL3
                    KP = KPV(KC)
                    KM = KMV(KC)
                
                    !=====================CONVECTION===================
                    ! $ \delta_x [ (h)^x \cdot G_X ] $ 
                    ! (I, J, K) = (I'+1, J, K)  - (I', J, K)   DX
                    DX_RHOHU =  ( ENTHALPY(IP,JC,KC) + ENTHALPY(IC,JC,KC) ) * ( G0_io(IP,JC,KC,1) + G_io(IP,JC,KC,1) ) - &
                                ( ENTHALPY(IC,JC,KC) + ENTHALPY(IM,JC,KC) ) * ( G0_io(IC,JC,KC,1) + G_io(IC,JC,KC,1) )
                    DX_RHOHU =   DX_RHOHU * COE1  ! *XND2CL*DXI*0.5_wp

               
                    ! $ \delta_y [ (h)^y \cdot G_Y ] $
                    ! (I, J, K) = (I, J'+1, K)  - (I, J', K)   DY
                    DY_RHOHV =  ( YCL2ND_WFF(JJP)*ENTHALPY(IC,JP,KC) +                         &
                                  YCL2ND_WFB(JJP)*ENTHALPY(IC,JC,KC) ) * ( G0_io(IC,JP,KC,2) + G_io(IC,JP,KC,2) ) -    &
                                ( YCL2ND_WFF(JJ) *ENTHALPY(IC,JC,KC) +                         &
                                  YCL2ND_WFB(JJ) *ENTHALPY(IC,JM,KC) ) * ( G0_io(IC,JC,KC,2) + G_io(IC,JC,KC,2) )
                    DY_RHOHV = DY_RHOHV * COE2 !*RCCI1(JJ)*DYFI(JJ)*0.5_wp

                    ! $ \delta_z [ (\rho h)^z \cdot w ] $
                    ! (I, J, K) = (I, J, K'+1)  - (I, J, K')   DZ
                    DZ_RHOHW =  ( ENTHALPY(IC,JC,KP) + ENTHALPY(IC,JC,KC) ) * ( G0_io(IC,JC,KP,3) + G_io(IC,JC,KP,3) ) - &
                                ( ENTHALPY(IC,JC,KC) + ENTHALPY(IC,JC,KM) ) * ( G0_io(IC,JC,KC,3) + G_io(IC,JC,KC,3) )
                    DZ_RHOHW = DZ_RHOHW * COE3 !*ZND2CL*DZI*RCCI2(JJ)*0.5_wp
              
              
                    !==========CONDUCTION=================================
                    ! $ \delta_x [k^x (\delta_x T)] $
                    ! (I, J, K) = (I'+1, J, K)  - (I', J, K)   DX
                    DX_KDT = ( THERMCONDT(IP,JC,KC)+THERMCONDT(IC,JC,KC) ) * ( TEMPERATURE(IP,JC,KC)-TEMPERATURE(IC,JC,KC) ) - &
                             ( THERMCONDT(IC,JC,KC)+THERMCONDT(IM,JC,KC) ) * ( TEMPERATURE(IC,JC,KC)-TEMPERATURE(IM,JC,KC) )
                    DX_KDT = DX_KDT * COE11 !XND2CL*DXI*DXI*CTHECD
              
                    ! $ \delta_y [k^y (\delta_y T)] $   
                    ! (I, J, K) = (I, J'+1, K)  - (I, J', K)   DY  
                    KDTy1 =  ( YCL2ND_WFF(JJP) * THERMCONDT(IC,JP,KC) +   &
                               YCL2ND_WFB(JJP) * THERMCONDT(IC,JC,KC) ) * & 
                             ( TEMPERATURE(IC,JP,KC) - TEMPERATURE(IC,JC,KC) ) * COE221        !DYCI(JJP)/RNDI(JJP)
                    KDTy2 =  ( YCL2ND_WFF(JJ)  * THERMCONDT(IC,JC,KC) +   &
                               YCL2ND_WFB(JJ)  * THERMCONDT(IC,JM,KC) ) * &
                             ( TEMPERATURE(IC,JC,KC) - TEMPERATURE(IC,JM,KC) ) * COE222        !DYCI(JJ)/RNDI(JJ)
                    IF(JJ==1 .and. BCWALLHEAT(Ibotwall)==isoFluxWall) THEN
                        KDTy2 =  -WALLFLUX(IC,Ibotwall) * COE222
                    ELSE IF(JJ==NCL2 .and. BCWALLHEAT(Itopwall)==isoFluxWall) THEN
                        KDTy1 =  -WALLFLUX(IC,Itopwall) * COE221
                    ELSE
                    END IF
                    DY_KDT = (KDTy1-KDTy2)* COE22  !DYFI(JJ)*RCCI1(JJ)*CTHECD
              
                    ! $ \delta_z [k^z (\delta_z T)] $                
                    DZ_KDT = ( THERMCONDT(IC,JC,KP)+THERMCONDT(IC,JC,KC) ) * ( TEMPERATURE(IC,JC,KP)-TEMPERATURE(IC,JC,KC) ) - &
                             ( THERMCONDT(IC,JC,KM)+THERMCONDT(IC,JC,KC) ) * ( TEMPERATURE(IC,JC,KC)-TEMPERATURE(IC,JC,KM) )
                    DZ_KDT = DZ_KDT * COE33  !ZND2CL*DZI*DZI*RCCI2(JJ)*CTHECD
             
                    !=======EXPLICIT RHS OF ENERGY EQUATION====================
                    RHS_ENERGY(IC,JC,KC) = -DX_RHOHU - DY_RHOHV - DZ_RHOHW + DX_KDT + DY_KDT + DZ_KDT

                    if(IC==4 .and. JJ<=4 .and. KC==4 .and. myid==0) then
                        write(*,*)'conx-e, cony-e, conz-e, dif-x, dif-y, dif-z', &
                        DX_RHOHU, DY_RHOHV, DZ_RHOHW, DX_KDT, DY_KDT, DZ_KDT
                    end if
                    if(IC==4 .and. JC<=4 .and. KC==4 .and. myid==0) then
                        write(*,*)'diy-dT', 'dify-k', &
                        KDTy2, THERMCONDT(IC,JC,KC)
                    end if
                                           
                END DO
            END DO
        END DO 
     
        !CALL DEBUG_WRT_LOCAL(RHS_ENERGY,1,N2DO(MYID)) !test
     
        RETURN
    END SUBROUTINE
