!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Kui He, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Calculate  the convection term in z direction.  
!> \f$ Hz = - \partial{u_3 u_j}/\partial{x_j} \f$ 
!> Ref: Eq. 4.37, 4.38 4.39 in Mehdi thesis.
!
!> @return CONVH = \f$ = - \partial{u_3 u_j}/\partial{x_j} \f$
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! ??/??/201?- Add pipe treatment, by Kui He.
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE CONVECTION_Z_tg
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE    
      
        INTEGER(4)  :: IC, IM, IP
        INTEGER(4)  :: JC, JM, JP, JJ, JJM, JJP
        INTEGER(4)  :: KC, KM, KP, KS, KSM
        REAL(WP)     :: COE1, COE2, COE3 
        REAL(WP)     :: H31, H32, H33
        REAL(WP)     :: q2s1, q2s2      
        REAL(WP)     :: h32n
        REAL(WP)     :: q2e, q2w      
        REAL(WP)     :: d11q1e
      
      
        RHSLLPHI_tg(:,:,:)=0.0_WP ! to store convection term in z direction.
        COE1 = XND2CL * ZND2CL * DXI 
        DO JC=1,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJM=JGMV(JJ)
            JJP=JGPV(JJ)
            COE2=DYFI(JJ)*ZND2CL
            COE3=ZND2CL*ZND2CL*DZI*RCCI2(jj)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
                DO IC=1,NCL1_tg
                    IP=IPV_tg(IC)
                    IM=IMV_tg(IC)           
                    
                    H31=((Q_tg(IP,JC,KC,1)+Q_tg(IP,JC,KM,1))*                  &
                         (Q_tg(IP,JC,KC,3)+Q_tg(IC,JC,KC,3))-                  &
                         (Q_tg(IC,JC,KC,1)+Q_tg(IC,JC,KM,1))*                  &
                         (Q_tg(IC,JC,KC,3)+Q_tg(IM,JC,KC,3)) )*COE1 
                               
                    H32=( ( Q_tg(IC,JP,KC,2)+Q_tg(IC,JP,KM,2) )*                &
                          ( YCL2ND_WFF(JJP) * Q_tg(IC,JP,KC,3)*RCCI1(jjp)+         &
                            YCL2ND_WFB(JJP) * Q_tg(IC,JC,KC,3)*RCCI1(jj)   )-      &       
                          ( Q_tg(IC,JC,KC,2)+Q_tg(IC,JC,KM,2) )*                &
                          ( YCL2ND_WFF(JJ)  * Q_tg(IC,JC,KC,3)*RCCI1(jj)+          &
                            YCL2ND_WFB(JJ)  * Q_tg(IC,JM,KC,3)*RCCI1(jjm) ) )*COE2   
                    
                    H33=((Q_tg(IC,JC,KP,3)+Q_tg(IC,JC,KC,3))*                   &
                         (Q_tg(IC,JC,KP,3)+Q_tg(IC,JC,KC,3))-                   &
                         (Q_tg(IC,JC,KC,3)+Q_tg(IC,JC,KM,3))*                   &
                         (Q_tg(IC,JC,KC,3)+Q_tg(IC,JC,KM,3)) )*COE3

                    RHSLLPHI_tg(IC,JC,KC)=-(H31+H32+H33)
                ENDDO
            ENDDO
        ENDDO
            
!>       @note below is only for pipe/coriolis force!!!!!!!!!!!!!!!      for athmut
        if (icase.eq.IPIPEC .or. icase.eq.IANNUL) then 
            DO jc=1,n2do(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)      !@
                do kc=1,NCL3
                    km=kmv(kc)
                    kp=kpv(kc)
                
                    do ic=1,NCL1_tg
                        if (jj.eq.1) then  
                            KS  = KSYM(KC)
                            KSM = KSYM(KM)
                            q2s1= (q_tg(ic,jp,kc,2)-q_tg(ic,jp,KS ,2))*0.50_WP*RNDI1(jjp)
                            q2s2= (q_tg(ic,jp,km,2)-q_tg(ic,jp,KSM,2))*0.50_WP*RNDI1(jjp)
                            q2e =  q_tg(ic,jp,kc,2)*RNDI1(jjp) + q2s1
                            q2w =  q_tg(ic,jp,km,2)*RNDI1(jjp) + q2s2
                        else
                            q2e = q_tg(ic,jp,kc,2)*RNDI1(jjp)+q_tg(ic,jc,kc,2)*RNDI1(jj)
                            q2w = q_tg(ic,jp,km,2)*RNDI1(jjp)+q_tg(ic,jc,km,2)*RNDI1(jj)
                        end if
                        h32n= q_tg(ic,jc,kc,3)*(q2e+q2w)*0.250_WP*RCCI1(jj)
                        d11q1e=(q2e-q2w)*DZI*RCCI1(jj)
                        RHSLLPHI_tg(ic,jc,kc)=-h32n+d11q1e/ren+RHSLLPHI_tg(ic,jc,kc)              
                    enddo
                enddo
            enddo 
     
        endif 
              
        RETURN
    END SUBROUTINE
    
!********************************************************************************************    
    SUBROUTINE CONVECTION_Z_io
        use thermal_info
        use mesh_info
        use flow_info
        use init_info
        IMPLICIT NONE    
      
        INTEGER(4)  :: IC, IM, IP
        INTEGER(4)  :: JC, JM, JP, JJ, JJP, JJM
        INTEGER(4)  :: KC, KM, KP, KS, KSM, NYI
      
        REAL(WP)     :: H31, H32, H33, H34
        REAL(WP)     :: H31F, H31B
        REAL(WP)     :: H32F, H32B
        REAL(WP)     :: H33F, H33B
        REAL(WP)     :: H34F, H34B, H34F1, H34F2, H34B1, H34B2
        REAL(WP)     :: COE1, COE2, COE30, COE31
        REAL(WP)     :: Ur1z1, Ur1z2, q2s1, q2s2, h32n, q2e, q2w
      
      
        RHSLLPHI_io(:,:,:)=0.0_WP ! to store convection term in z direction.
      
        COE1  = DXI * XND2CL * ZND2CL
        COE30 = DZI * ZND2CL * ZND2CL
        DO JC=1,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJP=JGPV(JJ)
            JJM=JGMV(jj)
            COE2 = DYFI(JJ) * ZND2CL
            COE31=COE30*RCCI2(JJ)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
                DO IC=1,NCL1_io
                    IP=IPV_io(IC)
                    IM=IMV_io(IC)
                   
                    ! \frac{\partial {\rho u w}}{\partial x}_{i,j,k'}
                    !{i'+1,j,k'}
                    H31F = ( G_IO(IP,JC,KC,1) + G_IO(IP,JC,KM,1) ) * &
                           ( Q_IO(IP,JC,KC,3) + Q_IO(IC,JC,KC,3) ) 
                    !{i',j,k'}             
                    H31B = ( G_IO(IC,JC,KC,1) + G_IO(IC,JC,KM,1) ) * &
                           ( Q_IO(IC,JC,KC,3) + Q_IO(IM,JC,KC,3) )  
                    !{i,j,k'}    
                    H31  = ( H31F - H31B) * COE1
                             
                    ! \frac{\partial {\rho v w}}{\partial y}_{i,j,k'}
                    !{i,j'+1,k'}
                    H32F = ( G_IO(IC,JP,KC,2) + G_IO(IC,JP,KM,2) ) * &
                           ( YCL2ND_WFF(JJP) * Q_IO(IC,JP,KC,3)*RCCI1(JJP) +   &
                             YCL2ND_WFB(JJP) * Q_IO(IC,JC,KC,3)*RCCI1(JJ) )        
                    !{i,j',k'}
                    H32B = ( G_IO(IC,JC,KC,2) + G_IO(IC,JC,KM,2) ) * &
                           ( YCL2ND_WFF(JJ)  * Q_IO(IC,JC,KC,3)*RCCI1(JJ) +   &
                             YCL2ND_WFB(JJ)  * Q_IO(IC,JM,KC,3)*RCCI1(JJM) ) 
                    !{i,j,k'}
                    H32   = ( H32F - H32B) * COE2
                 
                    ! \frac{\partial {\rho w w}}{\partial z}_{i,j,k'}
                    !(I,J,K)
                    H33F = ( G_IO(IC,JC,KP,3) + G_IO(IC,JC,KC,3) ) *  &
                           ( Q_IO(IC,JC,KP,3) + Q_IO(IC,JC,KC,3) )     
                    !(I,J,K-1)         
                    H33B = ( G_IO(IC,JC,KC,3) + G_IO(IC,JC,KM,3) ) *  &
                           ( Q_IO(IC,JC,KC,3) + Q_IO(IC,JC,KM,3) )     
                    !(I,J,K')
                    H33   = ( H33F - H33B) * COE31
                   
                    RHSLLPHI_io(IC,JC,KC) = -(H31+H32+H33)     !H for momentum Z direction. 
                    
                    !IF(JJ==1 .and. IC==1 .and. KC==1) write(*,'(A,4ES13.5)') 'convz',H31,H32,H33,RHSLLPHI_io(IC,JC,KC)
#ifdef DEBUG
  if(myid==0 .and. ic==1 .and. kc==1 .and. jc<=4) write(*,*) 'conz-31,32,33', JC, -H31, -H32, -H33
#endif
                END DO
            
            END DO
      
        END DO
      
      
        !below is only for pipe/coriolis force!!!!!!!!!!!!!!!
        if (icase.eq.IPIPEC .or. icase.eq.IANNUL ) then 
            DO jc=1,n2do(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)      !@
                do kc=1,NCL3
                    km=kmv(kc)
                    kp=kpv(kc)
                
                    do ic=1,NCL1_IO
                        if (jj.eq.1) then  
                            KS  = KSYM(KC)
                            KSM = KSYM(KM)
                            q2s1= (G_IO(ic,jp,kc,2)-G_IO(ic,jp,KS ,2))*0.50_WP*RNDI1(jjp)
                            q2s2= (G_IO(ic,jp,km,2)-G_IO(ic,jp,KSM,2))*0.50_WP*RNDI1(jjp)
                            q2e =  G_IO(ic,jp,kc,2)*RNDI1(jjp) + q2s1
                            q2w =  G_IO(ic,jp,km,2)*RNDI1(jjp) + q2s2
                        else
                            q2e = G_IO(ic,jp,kc,2)*RNDI1(jjp)+G_IO(ic,jc,kc,2)*RNDI1(jj)
                            q2w = G_IO(ic,jp,km,2)*RNDI1(jjp)+G_IO(ic,jc,km,2)*RNDI1(jj)
                        end if
                        h32n= Q_IO(ic,jc,kc,3)*(q2e+q2w)*YND2CL*ZND2CL*RCCI1(jj)
                        RHSLLPHI_IO(ic,jc,kc)=RHSLLPHI_IO(ic,jc,kc)-h32n            
                    enddo
                enddo
            enddo 
     
        endif
        
        
!        IF(ICASE==IPIPEC .OR. ICASE==IANNUL) THEN
!            NYI = 1
            
!            IF(MYID.EQ.0 .AND. ICASE.EQ.IPIPEC) THEN
!                !( Q_IO(IC,JC,KC,2)*RNDI1(JJ) ==>(Q_IO(ic,jp,kc,2)-Q_IO(ic,jp,KS, 2))*RNDI1(JJP)*0.5_WP
!                JC = 1
!                JJ = JCL2G(JC)
!                JM = JLMV(JC)
!                JP = JLPV(JC)
!                JJP= JGPV(JJ) !JPV(JJ)
!                JJM= JGMV(JJ) !JMV(JJ)
!                COE1 = ZND2CL*RNDI1(JJP)
            
!                DO KC=1,NCL3
!                    KP = KPV(KC)
!                    KM = KMV(KC)
!                    KS = KSYM(KC)
!                    KSM= KSYM(KM)
!                    DO IC=1,NCL1_io
!                        ! AT I, J'+1,K'
!                        H34F1=  ( G_IO(IC,JP,KC,2)+G_IO(IC,JP,KM,2) ) * COE1     !Gr average over theta
!                        H34F2=  ( YCL2ND_WFF(JJP) * Q_IO(IC,JP,KC,3)*RCCI1(JJP) + &
!                                  YCL2ND_WFB(JJP) * Q_IO(IC,JC,KC,3)*RCCI1(JJ) ) !G_theta average over r
!                        H34F = H34F1*H34F2
                        
!                        ! AT I, J', K'
!                        Ur1z1  = ( G_IO(ic,jp,kc,2)-G_IO(ic,jp,KS, 2))* COE1!*RNDI1(JJP)*0.5_WP
!                        Ur1z2  = ( G_IO(ic,jp,km,2)-G_IO(ic,jp,KSM,2))* COE1!*RNDI1(JJP)*0.5_WP
!                        H34B1=   (Ur1z1 + Ur1z2)*ZND2CL     !Gr average over theta
!                        H34B2=  ( YCL2ND_WFF(JJ) * Q_IO(IC,JC,KC,3)*RCCI1(JJ) +   &
!                                  YCL2ND_WFB(JJ) * Q_IO(IC,JM,KC,3)*RCCI1(JJM) ) !G_theta average over r
!                        H34B = H34B1*H34B2
                        
!                        ! AT I, J, K'
!                        H34 = (H34F + H34B)*YND2CL
                        
!                        RHSLLPHI_io(IC,JC,KC) = RHSLLPHI_io(IC,JC,KC) - H34
!                    END DO
!                END DO
                
!                NYI = 2
!            END IF
            
!            DO JC=NYI,N2DO(MYID)
!                JJ = JCL2G(JC)
!                JM = JLMV(JC)
!                JP = JLPV(JC)
!                JJP= JGPV(JJ) 
!                JJM= JGMV(JJ)
!                COE1 = ZND2CL*RNDI1(JJP)
!                COE2 = ZND2CL*RNDI1(JJ)
                
!                DO KC=1,NCL3
!                    KP = KPV(KC)
!                    KM = KMV(KC)
!                    DO IC=1,NCL1_io
!                        ! AT I, J'+1,K'
!                        H34F1=  ( G_io(IC,JP,KC,2)+G_io(IC,JP,KM,2) ) * COE1    !Gr average over theta
!                        H34F2=  ( YCL2ND_WFF(JJP) * Q_IO(IC,JP,KC,3)*RCCI1(JJP) + &
!                                  YCL2ND_WFB(JJP) * Q_IO(IC,JC,KC,3)*RCCI1(JJ) ) !G_theta average over r
!                        H34F = H34F1*H34F2
                        
!                        ! AT I, J', K'
!                        H34B1=  ( G_io(IC,JC,KC,2)+G_io(IC,JC,KM,2) ) * COE2     !Gr average over theta
!                        H34B2=  ( YCL2ND_WFF(JJ) * Q_IO(IC,JC,KC,3)*RCCI1(JJ) +   &
!                                  YCL2ND_WFB(JJ) * Q_IO(IC,JM,KC,3)*RCCI1(JJM) ) !G_theta average over r
!                        H34B = H34B1*H34B2
                        
!                        ! AT I, J, K'
!                        H34 = (H34F + H34B)*YND2CL
                        
!                        RHSLLPHI_io(IC,JC,KC) = RHSLLPHI_io(IC,JC,KC) - H34
!                    END DO
!                END DO
!            END DO
            
            
!        END IF
        
        !CALL DEBUG_WRT_LOCAL(RHSLLPHI_IO,1,N2DO(MYID),'conz') ! test
        
        RETURN
    END SUBROUTINE

