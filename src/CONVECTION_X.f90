!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Kui He, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Calculate  the convection term in x direction.  
!> \f$ Hx = - \partial{u_1 u_j}/\partial{x_j} \f$ 
!> Ref: Eq. 4.37, 4.38 4.39 in Mehdi thesis.
!
!> @return CONVH = \f$ = - \partial{u_1 u_j}/\partial{x_j} \f$
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! ??/??/201?- Add pipe treatment, by Kui He.
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************

    SUBROUTINE CONVECTION_X_tg
        use mesh_info
        use flow_info
        IMPLICIT NONE    
      
        INTEGER(4)  :: IC, IM, IP
        INTEGER(4)  :: JC, JM, JP, JJ, JJP
        INTEGER(4)  :: KC, KM, KP
        REAL(WP)     :: COE1, COE2, COE3
        REAL(WP)     :: H11, H12, H13
       
        Qtmp_tg(:,:,:) = 0.0_WP

        COE1 = XND2CL * XND2CL * DXI
        DO JC=1,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJP=JGPV(JJ)
            COE2 = DYFI(JJ) *XND2CL *RCCI1(jj)
            COE3 = XND2CL  * ZND2CL * DZI * RCCI2(jj)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
            
                DO IC=1,NCL1_tg
                    IP=IPV_tg(IC)
                    IM=IMV_tg(IC)
                    H11=( (Q_tg(IP,JC,KC,1)+Q_tg(IC,JC,KC,1)) *       &
                          (Q_tg(IP,JC,KC,1)+Q_tg(IC,JC,KC,1))     -   &
                          (Q_tg(IM,JC,KC,1)+Q_tg(IC,JC,KC,1)) *       &
                          (Q_tg(IM,JC,KC,1)+Q_tg(IC,JC,KC,1)) )*COE1
                          
                    H12=( (Q_tg(IC,JP,KC,2)+Q_tg(IM,JP,KC,2)) *                                &
                          (YCL2ND_WFF(JJP)*Q_tg(IC,JP,KC,1)+ &
                           YCL2ND_WFB(JJP)*Q_tg(IC,JC,KC,1) )- &
                          (Q_tg(IC,JC,KC,2)+Q_tg(IM,JC,KC,2)) *                                &
                          (YCL2ND_WFF(JJ) *Q_tg(IC,JC,KC,1)+ &
                           YCL2ND_WFB(JJ) *Q_tg(IC,JM,KC,1))  ) *COE2

                    H13=( (Q_tg(IC,JC,KP,3)+Q_tg(IM,JC,KP,3)) *      &
                          (Q_tg(IC,JC,KP,1)+Q_tg(IC,JC,KC,1))    -   &
                          (Q_tg(IC,JC,KC,3)+Q_tg(IM,JC,KC,3)) *      &
                          (Q_tg(IC,JC,KC,1)+Q_tg(IC,JC,KM,1)) )* COE3
                          
                    Qtmp_tg(IC,JC,KC)=-(H11+H12+H13)     !H for momentum x direction. 
                ENDDO
            ENDDO
        ENDDO
        
        RETURN
    END SUBROUTINE
    
!======================================================================================================================    
    SUBROUTINE CONVECTION_X_io
!>    The last time \rho is used. No combination of n+1 and n step of \rho      
      
        use thermal_info
        use mesh_info
        use flow_info
        IMPLICIT NONE    
      
        INTEGER(4)  :: IC, IM, IP
        INTEGER(4)  :: JC, JM, JP, JJ, JJP, NXI
        INTEGER(4)  :: KC, KM, KP
        REAL(WP)     :: H11, H12, H13
        REAL(WP)     :: H11F
        REAL(WP)     :: H11B
        REAL(WP)     :: H12F
        REAL(WP)     :: H12B
        REAL(WP)     :: H13F
        REAL(WP)     :: H13B
        REAL(WP)     :: COE1, COE2, COE3, COE32
       
        Qtmp_io(:,:,:) = 0.0_WP
        COE1 = DXI * XND2CL * XND2CL
        COE3 = DZI * XND2CL * ZND2CL
        
        NXI = 1
        IF(TGFLOWFLG) NXI = 2
        
        DO JC=1,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJP=JGPV(JJ)
            COE2 = DYFI(JJ) * XND2CL *RCCI1(jj)
            COE32= COE3*RCCI2(jj)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
            
                DO IC=NXI,NCL1_io
                    IP=IPV_io(IC)
                    IM=IMV_io(IC)
                   
                    ! \frac{\partial {\rho u u}}{\partial x}_{i',j,k}
                    !(I,J,K)
                    H11F = ( G_IO(IP,JC,KC,1) + G_IO(IC,JC,KC,1) ) *  &
                           ( Q_IO(IP,JC,KC,1) + Q_IO(IC,JC,KC,1) )     
                    !(I-1,J,K)         
                    H11B = ( G_IO(IC,JC,KC,1) + G_IO(IM,JC,KC,1) ) *  &
                           ( Q_IO(IC,JC,KC,1) + Q_IO(IM,JC,KC,1) )     
                    !(I', J,K)      
                    H11  = (H11F-H11B) * COE1                     
                    
                    ! \frac{\partial {\rho v u}}{\partial y}_{i',j,k}
                    !(I',J'+1,K)
                    H12F = ( G_IO(IC,JP,KC,2) +  G_IO(IM,JP,KC,2) ) *  &
                           ( YCL2ND_WFF(JJP) * Q_IO(IC,JP,KC,1) +      &
                             YCL2ND_WFB(JJP) * Q_IO(IC,JC,KC,1) )     
                    !(I',J',  K)         
                    H12B = ( G_IO(IC,JC,KC,2) +  G_IO(IM,JC,KC,2) ) *  &
                           ( YCL2ND_WFF(JJ)  * Q_IO(IC,JC,KC,1) +    &
                             YCL2ND_WFB(JJ)  * Q_IO(IC,JM,KC,1) )     
                    !(I',J,  K)         
                    H12  = ( H12F - H12B) * COE2
                             
                    ! \frac{\partial {\rho w u}}{\partial z}_{i',j,k} 
                    !(I',J,K'+1)
                    H13F = ( G_IO(IC,JC,KP,3) + G_IO(IM,JC,KP,3) ) *  &
                           ( Q_IO(IC,JC,KP,1) + Q_IO(IC,JC,KC,1) )     
                    !(I',J,K')       
                    H13B = ( G_IO(IC,JC,KC,3) + G_IO(IM,JC,KC,3) ) *  &
                           ( Q_IO(IC,JC,KC,1) + Q_IO(IC,JC,KM,1) )     
                    !(I',J,K)       
                    H13   = ( H13F - H13B) * COE32
                   
                    Qtmp_io(IC,JC,KC) = -(H11+H12+H13)     !H for momentum x direction. 
                    
                    !IF(myid==0 .and. JC==n2do(0)) write(*,'(A,3I4.1,5ES13.5)') 'convx',&
                    !JC,KC,IC,H12F,G_IO(IC,JP,KC,2),G_IO(IM,JP,KC,2),Q_IO(IC,JP,KC,1),Q_IO(IC,JC,KC,1)

#ifdef DEBUG
  if(myid==0 .and. ic==1 .and. kc==1 .and. jc<=4) write(*,*) 'conx-11,12,13',JC, -H11, -H12, -H13
#endif
                    
                END DO
            
            END DO
      
        END DO
        
        !CALL DEBUG_WRT_UWP_io ! test
        !CALL DEBUG_WRT_LOCAL(Qtmp_io,1,N2DO(MYID),'conx') ! test
        !  if(myid == 0) then        
        !      write(*,*) 'init1 ux', Q_io(:, 8, 8, 1), Q_io(:, 1, 8, 1)
        !      write(*,*) 'init1 uy', Q_io(:, 8, 8, 2), Q_io(:, 1, 8, 3)
        !      write(*,*) 'init1 uz', Q_io(:, 8, 8, 3), Q_io(:, 1, 8, 4)
        !      write(*,*) 'conx', Qtmp_io(:, 8, 8), Qtmp_io(:, 1, 8)
        !  end if

        RETURN
    END SUBROUTINE

