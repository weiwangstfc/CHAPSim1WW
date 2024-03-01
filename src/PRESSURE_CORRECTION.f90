    SUBROUTINE PRCALC_tg(NS)
        use init_info
        use flow_info
        use mesh_info
        IMPLICIT NONE     

        INTEGER(4),INTENT(IN)  :: NS
        REAL(WP) :: COE1, BE, LLPHI !,  COE21, COE22, COE3
        INTEGER(4) :: IC, JC, KC !, KM, KP, IM, IP, JM, JP, JJ, JJP
           
        BE=0.50_WP*TALP(NS)*DT*CVISC ! \alpha_1 * dt /2/Re in Eq.(A1d) of Mehdi thesis. 
        
        
        DO KC=1,NCL3
            DO JC=1,N2DO(MYID)
                DO IC=1,NCL1_TG         
                    LLphi = RHSLLPHI_tg(IC,JC,KC)                                     
                    PR_TG(IC,JC,KC)=PR_TG(IC,JC,KC)+DPH_TG(IC,JC,KC)-BE*LLphi
                ENDDO
            ENDDO
        ENDDO
            
        
!        IF(ICASE==ICHANL) THEN
!            DO KC=1,NCL3
!                KP=KPV(KC)
!                KM=KMV(KC)
!                DO JC=1,N2DO(MYID)
!                    JM=JLMV(JC)
!                    JP=JLPV(JC)
!                    JJ=JCL2G(JC)
!                    DO IC=1,NCL1_TG
!                        IP=IPV_TG(IC)
!                        IM=IMV_TG(IC)               
!                        LLphi = (       DPH_TG(IP,JC,KC) - &
!                                 2.0_WP*DPH_TG(IC,JC,KC) + &
!                                        DPH_TG(IM,JC,KC) )*DXQI+ &
!                                ( DPH_TG(IC,JP,KC)*APPH(JJ)+ &
!                                  DPH_TG(IC,JC,KC)*ACPH(JJ)+ &
!                                  DPH_TG(IC,JM,KC)*AMPH(JJ) )+ &
!                                (       DPH_TG(IC,JC,KP) - &
!                                 2.0_WP*DPH_TG(IC,JC,KC) + &
!                                        DPH_TG(IC,JC,KM) )*DZQI
!                        !or
!                        !LLphi = RHSLLPHI_tg(IC,JC,KC)                                     
!                        PR_TG(IC,JC,KC)=PR_TG(IC,JC,KC)+DPH_TG(IC,JC,KC)-BE*LLphi
!                    ENDDO
!                ENDDO
!            ENDDO
!        END IF
        
!        IF(ICASE.NE.ICHANL) THEN
!            DO JC=1,N2DO(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                JJP=JGPV(JJ)
!                COE3=DZQI*RCCI2(JJ)
!                COE21 = 1.0_WP/RNDI1(JJP)*DYCI(JJP)*DYFI(JJ)*RCCI1(JJ)
!                COE22 = 1.0_WP/RNDI1(JJ) *DYCI(JJ) *DYFI(JJ)*RCCI1(JJ)
!                DO KC=1,NCL3
!                    KP=KPV(KC)
!                    KM=KMV(KC)
                
!                    DO IC=1,NCL1_TG
!                        IP=IPV_TG(IC)
!                        IM=IMV_TG(IC)               
!                        LLphi= (         DPH_TG(IP,JC,KC)-  &
!                                  2.0_WP*DPH_TG(IC,JC,KC)+  &
!                                         DPH_TG(IM,JC,KC) )*DXQI+  &
!                               ( (DPH_TG(IC,JP,KC)-DPH_TG(IC,JC,KC) )*COE21-      &
!                                 (DPH_TG(IC,JC,KC)-DPH_TG(IC,JM,KC) )*COE22  ) +  &       
!                               (         DPH_TG(IC,JC,KP)-                         &
!                                  2.0_WP*DPH_TG(IC,JC,KC)+                         &
!                                         DPH_TG(IC,JC,KM) )*COE3  
!                        PR_TG(IC,JC,KC)=PR_TG(IC,JC,KC)+DPH_TG(IC,JC,KC)-BE*LLphi
!                    ENDDO
!                ENDDO                                                                     
!            ENDDO
!        END IF
      
      
        RETURN
    END SUBROUTINE
    
    
!************************************************************************    
    SUBROUTINE PRCALC_io(NS)
        use init_info
        use flow_info
        use mesh_info
        IMPLICIT NONE     

        INTEGER(4),INTENT(IN) :: NS 
        INTEGER(4) :: IC
        INTEGER(4) :: JC
        INTEGER(4) :: KC
        REAL(WP)   :: BE,LLphi

        IF(visthemflg==visimplicit) THEN
        
            BE=0.50_WP*TALP(NS)*DT*CVISC ! \alpha_1 * dt /2/Re in Eq.(A1d) of Mehdi thesis. 
            DO KC=1,NCL3
                DO JC=1,N2DO(MYID)
                    DO IC=1,NCL1_io          
                        LLphi = RHSLLPHI_io(IC,JC,KC)                                     
                        PR_io(IC,JC,KC)=PR_io(IC,JC,KC)+DPH_io(IC,JC,KC)-BE*LLphi
                    ENDDO
                ENDDO
            ENDDO
            
        ELSE IF(visthemflg == visexplicit) THEN
        
            DO KC=1,NCL3
                DO JC=1,N2DO(MYID)
                    DO IC=1,NCL1_io          
                        PR_io(IC,JC,KC)=PR_io(IC,JC,KC)+DPH_io(IC,JC,KC)
                    ENDDO
                ENDDO
            ENDDO
            
        ELSE
        END IF
      
        !BELOW IS FOR IMPLICIT VISCOUS AND ICASE=2
!        DO KC=1,NCL3
!            KP=KPV(KC)
!            KM=KMV(KC)
!            DO JC=1,N2DO(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                JJP=JGPV(JJ)
!                COE2=DYFI(JJ)*RCCI1(JJ)
!                COE3=DZQI*RCCI2(JJ)
!                DO IC=1,NCL1_io
!                    IP=IPV(IC)
!                    IM=IMV(IC)               
!                    LLphi= (         DPH(IP,JC,KC)-  &
!                              2.0_WP*DPH(IC,JC,KC)+  &
!                                     DPH(IM,JC,KC) )*DXQI+  &
!                           ( (DPH(IC,JP,KC)-DPH(IC,JC,KC) )*rc(jpp)*DYCI(JJP)-      &
!                             (DPH(IC,JC,KC)-DPH(IC,JM,KC) )*rc(jj) *DYCI(jj)  )*COE2 + &       
!                           (         DPH(IC,JC,KP)-                         &
!                              2.0_WP*DPH(IC,JC,KC)+                         &
!                                     DPH(IC,JC,KM) )*COE3  
!                    PR(IC,JC,KC)=PR(IC,JC,KC)+DPH(IC,JC,KC)-BE*LLphi
!                ENDDO
!            ENDDO                                                                     
!        ENDDO

        CALL INTFC_VARS1(NCL1S,NCL1E,NCL1S,NCL1E,PR_io)
        CALL BC_WALL_PR_io
      
        RETURN
    END SUBROUTINE
    
    
    SUBROUTINE CHECK_FFT_SOLVER
        use init_info
        use flow_info
        use mesh_info
        IMPLICIT NONE     

        INTEGER(4) :: IC, JC, KC , KM, KP, IM, IP, JM, JP, JJ, JJP
        REAL(WP) :: LLphiX, LLphiY, LLphiZ, LLPHI, COE0
        
            
            DO JC=1,N2DO(MYID)
                JM=JLMV(JC)
                JP=JLPV(JC)
                JJ=JCL2G(JC)
                DO KC=1,NCL3
                    KP=KPV(KC)
                    KM=KMV(KC)
                    DO IC=1,NCL1_io
                        IP=IPV_io(IC)
                        IM=IMV_io(IC)     
                        LLphiX = (       DPH_io(IP,JC,KC) - &
                                  2.0_WP*DPH_io(IC,JC,KC) + &
                                         DPH_io(IM,JC,KC) )*DXQI
                                        
                        LLphiY = DPH_io(IC,JP,KC)*APPH(JJ)+ &
                                 DPH_io(IC,JC,KC)*ACPH(JJ)+ &
                                 DPH_io(IC,JM,KC)*AMPH(JJ)
                            
                        LLphiZ = (       DPH_io(IC,JC,KP) - &
                                  2.0_WP*DPH_io(IC,JC,KC) + &
                                         DPH_io(IC,JC,KM) )*DZQI 
                                  
                        LLphi = LLphiX + LLphiY + LLphiZ
                        COE0  = DABS( LLphi-RHSLLPhi_io(IC,JC,KC) )
                        
                        IF(COE0.GT.1.0E-8_wp ) &
                        WRITE(*,'(A,4I4.1,3ES22.14)') '# FFT-comp', myid, JC, KC, IC,LLphi,&
                        RHSLLPhi_io(IC,JC,KC), COE0   
                    ENDDO
                ENDDO
            ENDDO

        RETURN
    END SUBROUTINE
    
    
    
!    SUBROUTINE GRAD_P_WALL_NONZERO
!        use init_info
!        use flow_info
!        use mesh_info
!        USE thermal_info
!        IMPLICIT NONE  
        
!        INTEGER(4) :: IC, JC, JJ, KC, N
!        REAL(WP)   :: COE1, COE2, COE3


!        DPDYWAL(:,:,:) = 0.0_wp
!        IF(thermlflg==0) RETURN
!        IF(BCWALLHEAT.NE.isoThermalWall )  RETURN
!        IF(MYID.GT.0 .and. MYID.LT.NPSLV) RETURN
    
!        IF(MYID.EQ.0) THEN
!            JC = 1
!            N  = 1
!            JJ = 1
!            COE1 = 2.0_WP/3.0_wp*DYFI(JJ)/DT/REN
            
!            COE3 = AMPH0
!        END IF
        
!        IF(MYID.EQ.NPSLV) THEN
!            JC = N2DO(MYID)
!            N  = 2
!            JJ = NCL2
!            COE1 = -2.0_WP/3.0_wp*DYFI(JJ)/DT/REN
            
!            COE3 = APPH0
!        END IF
        
        
        
!        DO IC=1,NCL1_io
!            COE2 = M_WAL_GV(IC,N)*D_WAL_GV(IC,N)
!            DO KC=1,NCL3
!                DPDYWAL(IC,KC,N) =  COE1*COE2/DENSITY(IC,JC,KC)/(2.0_WP*D_WAL_GV(IC,N)-DENSITY(IC,JC,KC))* &
!                (DENSITY(IC,JC,KC)-DENSITY0(IC,JC,KC))
!                !WRITE(*,*) 'DPDY',myid,IC,COE1,COE2,DPDYWAL(IC,KC,N),DPDYWAL(IC,KC,N)*COE3/DYFI(JJ)
!            END DO
            
!        END DO
    
    
!        RETURN
!    END SUBROUTINE
    
    
    
!    SUBROUTINE REVERSE_DPH_GHOST
!        use init_info
!        use flow_info
!        use mesh_info
!        IMPLICIT NONE     

!        REAL(WP) :: COE1, BE, LLPHI,   COE21, COE22, COE3
!        INTEGER(4) :: IC, JC, KC , KM, KP, IM, IP, JM, JP, JJ, JJP
        
!        IF(MYID==0) THEN
!            DO JC=1,N2DO(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                DO KC=1,NCL3
!                KP=KPV(KC)
!                KM=KMV(KC)
!                    DO IC=1,NCL1_io
!                        IP=IPV_io(IC)
!                        IM=IMV_io(IC)               
!                        LLphi = (       DPH_io(IP,JC,KC) - &
!                                 2.0_WP*DPH_io(IC,JC,KC) + &
!                                        DPH_io(IM,JC,KC) )*DXQI+ &
!                                ( DPH_io(IC,JP,KC)*APPH(JJ)+ &
!                                  DPH_io(IC,JC,KC)*ACPH(JJ)+ &
!                                  0.0_wp                   )+ &
!                                (       DPH_io(IC,JC,KP) - &
!                                 2.0_WP*DPH_io(IC,JC,KC) + &
!                                        DPH_io(IC,JC,KM) )*DZQI
!                         DPH_io(IC,JM,KC) = (RHSLLPhi_io(IC,JC,KC)-LLphi)/AMPH(JJ) !INFINITY
!                         write(*,*) myid, JC, KC, IC, RHSLLPhi_io(IC,JC,KC), LLphi, (RHSLLPhi_io(IC,JC,KC)-LLphi), AMPH(JJ)
!                    ENDDO
!                ENDDO
!            ENDDO
!        END IF
        
!        IF(MYID==NPSLV) THEN
!            DO JC=1,N2DO(MYID)
!                JM=JLMV(JC)
!                JP=JLPV(JC)
!                JJ=JCL2G(JC)
!                DO KC=1,NCL3
!                KP=KPV(KC)
!                KM=KMV(KC)
!                    DO IC=1,NCL1_io
!                        IP=IPV_io(IC)
!                        IM=IMV_io(IC)               
!                        LLphi = (       DPH_io(IP,JC,KC) - &
!                                 2.0_WP*DPH_io(IC,JC,KC) + &
!                                        DPH_io(IM,JC,KC) )*DXQI+ &
!                                ( 0.0_wp+ &
!                                  DPH_io(IC,JC,KC)*ACPH(JJ)+ &
!                                  DPH_io(IC,JM,KC)*AMPH(JJ) )+ &
!                                (       DPH_io(IC,JC,KP) - &
!                                 2.0_WP*DPH_io(IC,JC,KC) + &
!                                        DPH_io(IC,JC,KM) )*DZQI
!                        DPH_io(IC,JP,KC)= (RHSLLPhi_io(IC,JC,KC)-LLphi)/APPH(JJ) !INFINITY
!                        write(*,*) myid, JC, KC, IC, RHSLLPhi_io(IC,JC,KC), LLphi, (RHSLLPhi_io(IC,JC,KC)-LLphi), APPH(JJ)
!                    ENDDO
!                ENDDO
!            ENDDO
!        END IF
    
    
!        RETURN
!    END SUBROUTINE
