!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Calculate the divergence of velocity 
!> \f$ \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}\f$ 
!> Ref: Eq. 4.37, 4.38 4.39 in Mehdi thesis.
!
!> @return RHSLLPHI The RHS of dp equation.  
!
!> @todo
!> Check the size of RHSLLPHI 
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! 20/01/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE DIVG_tg(NS)
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE

        INTEGER(4),INTENT(IN)  :: NS
      
        INTEGER(4)   :: KC, KP
        INTEGER(4)   :: JC, JP, JJ
        INTEGER(4)   :: IC, IP
        REAL(WP)      :: DVIGVELO
        REAL(WP)      :: COE0, COE1, COE2
       
        COE0 = 1.0_WP/(DT*TALP(NS))
        RHSLLPHI_tg = 0.0_WP   
        
        DO JC=1,N2DO(MYID)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            COE1 = DXI/RCCI2(jj)
            COE2 = DYFI(JJ)/RCCI1(jj)
            DO KC=1,NCL3
                KP=KPV(KC)
                DO IC=1,NCL1_tg
                    IP=IPV_tg(IC)
                    DVIGVELO=(Q_tg(IP,JC,KC,1)-Q_tg(IC,JC,KC,1))*COE1    &
                            +(Q_tg(IC,JP,KC,2)-Q_tg(IC,JC,KC,2))*COE2   & 
                            +(Q_tg(IC,JC,KP,3)-Q_tg(IC,JC,KC,3))*DZI
                    RHSLLPHI_tg(IC,JC,KC)=DVIGVELO*COE0
                ENDDO
            ENDDO
        ENDDO
      
        RETURN
      
    END SUBROUTINE
    
!*******************************************************************************************    
    SUBROUTINE DIVG_io(NS)
!>    @NOTE
!>    1) If the velocity field is calculated based on time-fixed density field, 
!>       the density time rate is excluded from the continuity equation.
        use thermal_info
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE

        INTEGER(4),INTENT(IN) :: NS
      
        INTEGER(4)   :: IC,JC,KC
        INTEGER(4)   :: IP,JP,KP,JJ
        REAL(WP)       :: DIVX
        REAL(WP)       :: DIVY
        REAL(WP)       :: DIVZ
        REAL(WP)       :: DTRHO, DQCAP
        REAL(WP)       :: COE0, COE1, COE2, COE4
      
        RHSLLPHI_io = 0.0_WP   
        DTRHO = 0.0_WP
        DQCAP = 0.0_WP

        IF(weightedpressure==1) THEN
            COE0 = 1.0_wp/DT/TALP(NS)/(0.25_wp+pres_epslon)
        ELSE
            COE0 = 1.0_wp/DT/TALP(NS)
        END IF
        
        DO JC=1,N2DO(MYID)
            JP = JLPV(JC)
            JJ = JCL2G(JC)
            COE1 = DXI/RCCI2(jj)
            COE2 = DYFI(JJ)/RCCI1(jj)
            COE4 = 1.0_WP/ DT /RCCI2(jj)
            DO KC=1,NCL3
                KP=KPV(KC)
                !==================MAX. DIV IN THE MAIN DOMAIN=================
                DO IC=1,NCL1_io
                    IP=IPV_io(IC)
               
                    !==========d(\rho u)/dx at (i,j,k)===========================
                    DIVX  = ( G_IO(IP,JC,KC,1) - G_IO(IC,JC,KC,1) ) * COE1
               
                    !==========d(\rho v)/dy at (i,j,k)===========================
                    DIVY  = ( G_io(IC,JP,KC,2) - G_io(IC,JC,KC,2) ) * COE2
               
                    !==========d(\rho w)/dz at (i,j,k)===========================
                    DIVZ  = ( G_io(IC,JC,KP,3) - G_io(IC,JC,KC,3) ) * DZI 
               
                    !==========D \RHO/ DT========================================
                    !DTRHO = ( DENSITY(IC,JC,KC) - DENSITYP(IC,JC,KC) ) * COE4
                    DTRHO = DrhoDtP(IC,JC,KC)/RCCI2(jj)
               
                    !===========TOTAL============================================
                    DQCAP = DIVX + DIVY + DIVZ + DTRHO 

                    RHSLLPHI_io(IC,JC,KC)=DQCAP*COE0

!                    IF(JJ==1)    THEN
!                        RHSLLPHI_io(IC,JC,KC)=RHSLLPHI_io(IC,JC,KC)+AMPH0*DPDYWAL(IC,KC,1)/DYFI(1)
!                        !write(*,*) JJ, DQCAP*COE0, AMPH0*DPDYWAL(IC,KC,1)*DYFI(1), RHSLLPHI_io(IC,JC,KC)
!                    END IF
!                    IF(JJ==NCL2) THEN
!                        RHSLLPHI_io(IC,JC,KC)=RHSLLPHI_io(IC,JC,KC)-APPH0*DPDYWAL(IC,KC,2)/DYFI(NCL2)
!                        !write(*,*) JJ, DQCAP*COE0, AMPH0*DPDYWAL(IC,KC,2)*DYFI(NCL2), RHSLLPHI_io(IC,JC,KC)
!                    END IF
                    
!                    IF(MYID.eq.0) &
!                    WRITE(*,'(A,4I3.1,7ES13.5)') 'divgoy', myid, JC, KC, IC, &
!                    DIVX/COE1,DIVY/COE2,DIVZ/DZI, &
!                    DENSITY(IC,JC,KC),DENSITY0(IC,JC,KC), &
!                    DTRHO/COE4,DQCAP
                    
                ENDDO
            ENDDO
        ENDDO
       !if(myid==0) write(*,*) 'rhs',  RHSLLPHI_io(8,8,8), RHSLLPHI_io(16,8,8), RHSLLPHI_io(32,8,8)
        
        
        RETURN
      
    END SUBROUTINE
      
!************************************************************************************************************      
    SUBROUTINE DIVG_U_io
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        INTEGER(4)   :: IC,JC,KC!, NYI
        INTEGER(4)   :: IP,JP,KP,JJ
        REAL(WP)       :: DIVX
        REAL(WP)       :: DIVY
        REAL(WP)       :: DIVZ
        REAL(WP)       :: COE2, COE3
      
        
        DivU_io = 0.0_wp
        
        
        DO JC=1,N2DO(MYID)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            COE2 = DYFI(JJ) * RCCI1(JJ)
            COE3 = DZI * RCCI2(JJ)
            DO KC=1,NCL3
                KP=KPV(KC)
                !==================MAX. DIV IN THE MAIN DOMAIN=================
                DO IC=1,NCL1_io
                    IP=IPV_io(IC)
               
                    !==========d(\rho u)/dx at (i,j,k)===========================
                    DIVX  = ( Q_IO(IP,JC,KC,1) - Q_IO(IC,JC,KC,1) ) * DXI
                    !==========d(\rho v)/dy at (i,j,k)===========================
                    DIVY  = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * COE2
                    !==========d(\rho w)/dz at (i,j,k)=========================== 
                    DIVZ  = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * COE3
               
                    !===========TOTAL============================================
                    DivU_io(IC,JC,KC)= (DIVX + DIVY + DIVZ)/3.0_wp
                    
                    !write(*,'(A,3I4.1,4ES13.5)') 'divxyz',JJ,IC,KC,DIVX, DIVY, DIVZ, DivU_io(IC,JC,KC)
                ENDDO
            ENDDO
        ENDDO
        
!        !===========================below bottom wall JC=0 is not used============================
!        IF(MYID==0) THEN
!            JC=0
!            JP=1
!            JJ=1
!            COE2 = DYFI(JJ) * RCCI1(JJ)
!            COE3 = DZI * RCCI2(JJ)
!            DO KC=1,NCL3
!                KP=KPV(KC)
!                !==================MAX. DIV IN THE MAIN DOMAIN=================
!                DO IC=1,NCL1_io
!                    IP=IPV_io(IC)
               
!                    !==========d(\rho u)/dx at (i,j,k)===========================
!                    DIVX  = ( Q_IO(IP,JC,KC,1) - Q_IO(IC,JC,KC,1) ) * DXI
!                    !==========d(\rho v)/dy at (i,j,k)===========================
!                    DIVY  = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * COE2
!                    !==========d(\rho w)/dz at (i,j,k)=========================== 
!                    DIVZ  = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * COE3
               
!                    !===========TOTAL============================================
!                    DivU_io(IC,JC,KC)= (DIVX + DIVY + DIVZ)/3.0_wp
                    
!                    !write(*,'(A,3I4.1,4ES13.5)') 'divxyz',JJ,IC,KC,DIVX, DIVY, DIVZ, DivU_io(IC,JC,KC)
!                ENDDO
!            ENDDO
!        END IF
        
!        !===========================below top wall JC=0 is not used============================
!        IF(MYID==NPSLV) THEN
!            JC=N2DO(MYID)+1
!            JP=N2DO(MYID)+2
!            JJ=NND2
!            COE2 = DYFI(NCL2) * RCCI1(JJ)
!            COE3 = DZI * RCCI2(JJ)
!            DO KC=1,NCL3
!                KP=KPV(KC)
!                !==================MAX. DIV IN THE MAIN DOMAIN=================
!                DO IC=1,NCL1_io
!                    IP=IPV_io(IC)
               
!                    !==========d(\rho u)/dx at (i,j,k)===========================
!                    DIVX  = ( Q_IO(IP,JC,KC,1) - Q_IO(IC,JC,KC,1) ) * DXI
!                    !==========d(\rho v)/dy at (i,j,k)===========================
!                    DIVY  = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * COE2
!                    !==========d(\rho w)/dz at (i,j,k)=========================== 
!                    DIVZ  = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * COE3
               
!                    !===========TOTAL============================================
!                    DivU_io(IC,JC,KC)= (DIVX + DIVY + DIVZ)/3.0_wp
                    
!                    !write(*,'(A,3I4.1,4ES13.5)') 'divxyz',JJ,IC,KC,DIVX, DIVY, DIVZ, DivU_io(IC,JC,KC)
!                ENDDO
!            ENDDO
!        END IF
        
        !===========================================================================================
        
        !CALL INTFC_ALL_DIVU_io
        CALL INTFC_VARS1(1,NCL1_io,NCL1S,NCL1E,DivU_io)
       
        RETURN
      
    END SUBROUTINE
