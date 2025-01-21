!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Kui He, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Calculate  the convection term in y direction.  
!> \f$ Hy = - \partial{u_2 u_j}/\partial{x_j} \f$ 
!> Ref: Eq. 4.37, 4.38 4.39 in Mehdi thesis.
!
!> @return CONVH = \f$ = - \partial{u_2 u_j}/\partial{x_j} \f$
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! ??/??/201?- Add pipe treatment, by Kui He.
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE CONVECTION_Y_tg
        use flow_info
        use mesh_info
        use init_info
        IMPLICIT NONE    
      
        INTEGER(4)  :: IC, IM, IP
        INTEGER(4)  :: JC, JM, JP, JJ, JJM, JJP
        INTEGER(4)  :: KC, KM, KP
        INTEGER(4)  :: KS
        INTEGER(4)  :: NYI
        REAL(WP)     :: H21, H22, H23    
        REAL(WP)     :: q2s1     
        REAL(WP)     :: h23n
        REAL(WP)     :: q1e, q1w     
        REAL(WP)     :: d11q2e
        REAL(WP)     :: COE1,COE2, COE30, COE31
      
        DPH_tg(:,:,:)=0.0_WP  ! to store the convection term in y direction.
     
        NYI=1
        IF (MYID.EQ.0) NYI=2
     
        COE1  = DXI * XND2CL
        COE30 = DZI * ZND2CL
        DO JC=NYI,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            jjm=JGMV(jj)
            jjp=JGPV(jj)
            COE2  = DYCI(JJ)*YND2CL*YND2CL
            COE31 = COE30 *RNDI1(JJ)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
                DO IC=1,NCL1_tg
       
                    IP=IPV_tg(IC)
                    IM=IMV_tg(IC)
                    ! I,J',K = (I'+1, J', K)  vs  (I', J', K) 
                    H21=( ( YCL2ND_WFF(JJ) * Q_tg(IP,JC,KC,1)+                      &
                            YCL2ND_WFB(JJ) * Q_tg(IP,JM,KC,1) )*                     &
                          ( Q_tg(IP,JC,KC,2)+Q_tg(IC,JC,KC,2) )-                     &
                          ( YCL2ND_WFF(JJ) * Q_tg(IC,JC,KC,1)+                      &
                            YCL2ND_WFB(JJ) * Q_tg(IC,JM,KC,1) ) *                    &
                          ( Q_tg(IC,JC,KC,2)+Q_tg(IM,JC,KC,2) ) )*COE1
                    if (icase.eq.IPIPEC .and. jj.eq.2) then
                        KS = KSYM(KC)
                        q2s1= (q_tg(ic,jc,kc,2) - q_tg(ic,jc,KS,2))*0.50_WP*RNDI1(jj)
                        h22=((q_tg(ic,jp,kc,2)*RNDI1(jjp)+q_tg(ic,jc,kc,2)*RNDI1(jj))*       &
                             (q_tg(ic,jp,kc,2)+q_tg(ic,jc,kc,2)) -                  &
                             (q_tg(ic,jc,kc,2)*RNDI1(jj)+q2s1)*                     &
                             (q_tg(ic,jc,kc,2)+q_tg(ic,jm,kc,2)) )*COE2  
                    else
                        ! I,J',K = (I, J, K)  vs  (I, J-1, K) 
                        H22=((Q_tg(IC,JP,KC,2)*RNDI1(jjp)+Q_tg(IC,JC,KC,2)*RNDI1(jj))*   &
                             (Q_tg(IC,JP,KC,2)+Q_tg(IC,JC,KC,2))-                  &
                             (Q_tg(IC,JC,KC,2)*RNDI1(jj)+Q_tg(IC,JM,KC,2)*RNDI1(jjm))*   &
                             (Q_tg(IC,JC,KC,2)+Q_tg(IC,JM,KC,2)) )*COE2
                    endif
                    ! I,J',K = (I, J', K'+1)  vs  (I, J', K') 
                    H23=( ( YCL2ND_WFF(JJ) * Q_tg(IC,JC,KP,3)*RCCI1(jj)+      &
                            YCL2ND_WFB(JJ) * Q_tg(IC,JM,KP,3)*RCCI1(jjm) )*   &
                          ( Q_tg(IC,JC,KP,2)+Q_tg(IC,JC,KC,2) )- &
                          ( YCL2ND_WFF(JJ) * Q_tg(IC,JC,KC,3)*RCCI1(jj)+      &
                            YCL2ND_WFB(JJ) * Q_tg(IC,JM,KC,3)*RCCI1(jjm) )*   &
                          ( Q_tg(IC,JC,KC,2)+Q_tg(IC,JC,KM,2) ) )*COE31        !@
                    DPH_tg(IC,JC,KC)=-(H21+H22+H23)   !H for momentum y direction.
                ENDDO
            ENDDO
        ENDDO


        if (icase.eq.IPIPEC .OR. icase.eq.IANNUL) then            !@   Centripetal force
            do jc=NYI,n2do(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                do kc=1,NCL3
                    km=kmv(kc)
                    kp=kpv(kc)
                    do ic=1,NCL1_tg
!                        q1e=q_tg(ic,jc,kp,3)*RCCI1(jj)+q_tg(ic,jm,kp,3)*RCCI1(jjm) 
!                        q1w=q_tg(ic,jc,kc,3)*RCCI1(jj)+q_tg(ic,jm,kc,3)*RCCI1(jjm)
!                        h23n=( (q1e+q1w)*0.250_WP )**2
!                        d11q2e= -(q1e-q1w)*DZI*RNDI1(jj)
!                        DPH_tg(ic,jc,kc)=DPH_tg(ic,jc,kc)+h23n+d11q2e/ren
                        q1e=q_tg(ic,jc,kp,3)*RCCI1(jj)*YCL2ND_WFF(JJ)+q_tg(ic,jm,kp,3)*RCCI1(jjm)*YCL2ND_WFB(JJ) 
                        q1w=q_tg(ic,jc,kc,3)*RCCI1(jj)*YCL2ND_WFF(JJ)+q_tg(ic,jm,kc,3)*RCCI1(jjm)*YCL2ND_WFB(JJ)
                        h23n=( (q1e+q1w)*ZND2CL )**2
                        d11q2e= -2.0_wp*(q1e-q1w)*DZI*RNDI1(jj)
                        DPH_tg(ic,jc,kc)=DPH_tg(ic,jc,kc)+h23n+d11q2e/ren
                    enddo
                enddo
            enddo
        endif
        
        RETURN
    END SUBROUTINE



!***********************************************************************
    SUBROUTINE CONVECTION_Y_io
        use thermal_info
        use mesh_info
        use flow_info
        use init_info
        IMPLICIT NONE    
        
        INTEGER(4)  :: IC, IM, IP
        INTEGER(4)  :: JC, JM, JP, JJ, JJM, JJP
        INTEGER(4)  :: KC, KM, KP, KS
        INTEGER(4)  :: NYI   
      
        REAL(WP)     :: H21, H22, H23, H24, h23n   
        REAL(WP)     :: H21F
        REAL(WP)     :: H21B
        REAL(WP)     :: H22F
        REAL(WP)     :: H22B
        REAL(WP)     :: H23F
        REAL(WP)     :: H23B
        REAL(WP)     :: H24F, H24B, H24F1, H24F2, H24B1, H24B2
        REAL(WP)     :: COE1, COE2, COE30, COE31, COE41, COE42
        REAL(WP)     :: Ur1, q1e, q1w, g1e, g1w
      
        DPH_io(:,:,:)=0.0_WP  ! to store the convection term in y direction.
     
        NYI=1
        IF (MYID.EQ.0 .and. ICASE.NE.IBOX3P) NYI=2
      
        COE1  = DXI * XND2CL
        COE30 = DZI * ZND2CL
        DO JC=NYI,N2DO(MYID)
            JM=JLMV(JC)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            JJM=JGMV(JJ)
            JJP=JGPV(JJ)
            COE2 = DYCI(JJ) * YND2CL * YND2CL
            COE31= COE30 *RNDI1(JJ)
            DO KC=1,NCL3
                KM=KMV(KC)
                KP=KPV(KC)
                
                DO IC=1,NCL1_io
                    IP=IPV_io(IC)
                    IM=IMV_io(IC)
                   
                    ! \frac{\partial {\rho u v}}{\partial x}_{i,j',k}
                    !{i'+1,j',k}
                    H21F = ( YCL2ND_WFF(JJ) * G_IO(IP,JC,KC,1) +   &
                             YCL2ND_WFB(JJ) * G_IO(IP,JM,KC,1) ) * &
                           ( Q_IO(IP,JC,KC,2) + Q_IO(IC,JC,KC,2) )         
                    !{i',j',k}         
                    H21B = ( YCL2ND_WFF(JJ) * G_IO(IC,JC,KC,1) +   &
                             YCL2ND_WFB(JJ) * G_IO(IC,JM,KC,1) ) * &
                           ( Q_IO(IC,JC,KC,2) + Q_IO(IM,JC,KC,2) )         
                    !{i,j',k}       
                    H21  = ( H21F - H21B) * COE1 
                   
                    ! \frac{\partial {\rho v v}}{\partial y}_{i,j',k}
                    !{i,j,k}
                    H22F = ( G_IO(IC,JP,KC,2) + G_IO(IC,JC,KC,2) ) * &
                           ( Q_IO(IC,JP,KC,2)*RNDI1(JJP) + Q_IO(IC,JC,KC,2)*RNDI1(JJ) )         
                    !{i,j-1,k}      
                    IF(icase.eq.IPIPEC .and. jj.eq.2) THEN
                        KS = KSYM(KC)
                        Ur1  = ( Q_IO(IC,JC,KC,2) - Q_IO(IC,JC,KS,2) )*0.50_WP*RNDI1(JJ)
                        H22B = ( G_IO(IC,JC,KC,2) + G_IO(IC,JM,KC,2) ) * &
                               ( Q_IO(IC,JC,KC,2)*RNDI1(JJ) + Ur1 )     
                    ELSE
                        H22B = ( G_IO(IC,JC,KC,2) + G_IO(IC,JM,KC,2) ) * &
                               ( Q_IO(IC,JC,KC,2)*RNDI1(JJ) + Q_IO(IC,JM,KC,2)*RNDI1(JJM) )   
                    END IF
                    !{i,j',k}       
                    H22   = ( H22F - H22B) * COE2
    
                    ! \frac{\partial {\rho w v}}{\partial x}_{i,j',k}
                    !{i,j',k'+1}
                    H23F = ( YCL2ND_WFF(JJ) * G_IO(IC,JC,KP,3)*RCCI1(JJ)   +   &
                             YCL2ND_WFB(JJ) * G_IO(IC,JM,KP,3)*RCCI1(JJM) ) * &
                           ( Q_IO(IC,JC,KP,2) + Q_IO(IC,JC,KC,2) ) !*RNDI1(JJ)     
                    !{i,j',k'}         
                    H23B = ( YCL2ND_WFF(JJ) * G_IO(IC,JC,KC,3)*RCCI1(JJ) +   &
                             YCL2ND_WFB(JJ) * G_IO(IC,JM,KC,3)*RCCI1(JJM) ) * &
                           ( Q_IO(IC,JC,KC,2) + Q_IO(IC,JC,KM,2) ) !*RNDI1(JJ)
                    !{i,j',k}       
                    H23   = ( H23F - H23B) * COE31
                   
                    DPH_io(IC,JC,KC) = -(H21+H22+H23)   !H for momentum y direction.
                    ! IF(myid==0 .and. JJ<=4 .and. IC==4 .and. KC==4) &
                    ! write(*,*) 'convy', JJ, H21,H22,H23,DPH_io(IC,JC,KC)
#ifdef DEBUG
  if(myid==0 .and. ic==1 .and. kc==1 .and. jc<=4) write(*,*) 'cony-21,22,23', JC, -H21, -H22, -H23
#endif
                END DO
            END DO
        END DO
        
        !@   Centripetal force
        if (icase.eq.IPIPEC .OR. icase.eq.IANNUL) then            !@   Centripetal force
            do jc=NYI,n2do(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                JJM=JGMV(JJ)
                JJP=JGPV(JJ)
                do kc=1,NCL3
                    km=kmv(kc)
                    kp=kpv(kc)
                    do ic=1,NCL1_io
                        q1e=q_IO(ic,jc,kp,3)*RCCI1(jj)*YCL2ND_WFF(JJ)+q_IO(ic,jm,kp,3)*RCCI1(jjm)*YCL2ND_WFB(JJ)
                        q1w=q_IO(ic,jc,kc,3)*RCCI1(jj)*YCL2ND_WFF(JJ)+q_IO(ic,jm,kc,3)*RCCI1(jjm)*YCL2ND_WFB(JJ)
                        g1e=g_IO(ic,jc,kp,3)*RCCI1(jj)*YCL2ND_WFF(JJ)+g_IO(ic,jm,kp,3)*RCCI1(jjm)*YCL2ND_WFB(JJ)
                        g1w=g_IO(ic,jc,kc,3)*RCCI1(jj)*YCL2ND_WFF(JJ)+g_IO(ic,jm,kc,3)*RCCI1(jjm)*YCL2ND_WFB(JJ)
                        
                        !h23n=( (q1e+q1w)*0.250_WP )**2
                        h23n=(q1e+q1w)*ZND2CL*(g1e+g1w)*ZND2CL
                        DPH_IO(ic,jc,kc)=DPH_IO(ic,jc,kc)+h23n
                    enddo
                enddo
            enddo
        endif
        
!        !@   Centripetal force
!        if (icase.eq.IPIPEC .OR. icase.eq.IANNUL) then   
         
!            do jc=NYI,n2do(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                JJM=JGMV(JJ)
!                JJP=JGPV(JJ)
!                COE41 = ZND2CL * ZND2CL * RCCI2(JJ)  
!                COE42 = ZND2CL * ZND2CL * RCCI2(JJM)   
!                do kc=1,NCL3
!                    km=kmv(kc)
!                    kp=kpv(kc)
                
!                    do ic=1,NCL1_io
!                        ! AT I,J,  K
!                        H24F1 = ( G_IO(IC,JC,KP,3) + G_IO(IC,JC,KC,3) ) !* ZND2CL *RCCI1(JJ)
!                        H24F2 = ( Q_IO(IC,JC,KP,3) + Q_IO(IC,JC,KC,3) ) !* ZND2CL *RCCI1(JJ)
!                        H24F  = H24F1 * H24F2 * COE41 
!                        ! AT I,J-1,K
!                        H24B1 = ( G_IO(IC,JM,KP,3) + G_IO(IC,JM,KC,3) ) !* ZND2CL * RCCI1(JJM)
!                        H24B2 = ( Q_IO(IC,JM,KP,3) + Q_IO(IC,JM,KC,3) ) !* ZND2CL * RCCI1(JJM)
!                        H24B  = H24B1 * H24B2 * COE42
                        
!                        ! AT I, J',K
!                        H24 = YCL2ND_WFF(JJ) * H24F + YCL2ND_WFB(JJ) * H24B

!                        DPH_IO(ic,jc,kc)=DPH_IO(ic,jc,kc)+H24
!                    enddo
!                enddo
!            enddo
!        endif
        
        !CALL DEBUG_WRT_LOCAL(DPH_IO,1,N2DO(MYID),'cony') ! test
        
        RETURN
    END SUBROUTINE
