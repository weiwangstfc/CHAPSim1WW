    SUBROUTINE RHS_MOM_EXPLICIT_io(NS,IDR) ! not using other Gs or Qs in the current step
        USE FLOW_INFO
        USE THERMAL_INFO
        USE MESH_INFO
        use init_info
        IMPLICIT NONE
     
        INTEGER(4),INTENT(IN) :: NS 
        INTEGER(4),INTENT(IN) :: IDR 
     
        REAL(WP)    :: EXPLT_io(NCL1_io,N2DO(0),NCL3) 
        REAL(WP)    :: DEN0
        REAL(WP)    :: COE2      
        REAL(WP)    :: RHSC, RHSL
        REAL(WP)    :: PGM     
        REAL(WP)    :: DPGRNS, INTGRHSY, INTGRHSY_WORK
        INTEGER(4) :: NXI, NYI
        INTEGER(4) :: IC, IM
        INTEGER(4) :: JC, JM, JJ
        INTEGER(4) :: KC, KM

     
        !===============SET UP INDEX WITHOUT B.C.===========================
        NXI=1
        IF(TGFLOWFLG .AND. (IDR.EQ.1) )  NXI = 2
        NYI=1
        IF(IDR.EQ.2 .AND. MYID.EQ.0 .AND. ICASE.NE.IBOX3P) NYI=2
      
        !===============SET UP THE CONVECTION TERM INTO ONE VARIABLE========
        EXPLT_io = 0.0_WP
        IF(IDR.EQ.1) THEN
      
            DO IC=NXI,NCL1_io
                DO JC=NYI,N2DO(MYID)
                    DO KC=1,NCL3
                        EXPLT_io(IC,JC,KC) = Qtmp_io(IC,JC,KC)
                    END DO
                END DO
            END DO
      
        ELSE IF(IDR .EQ. 2) THEN
      
            DO IC=NXI,NCL1_io
                DO JC=NYI,N2DO(MYID)
                    DO KC=1,NCL3
                        EXPLT_io(IC,JC,KC) = DPH_io(IC,JC,KC)
                    END DO
                END DO
            END DO
      
        ELSE IF(IDR .EQ. 3) THEN
      
            DO IC=NXI,NCL1_io
                DO JC=NYI,N2DO(MYID)
                    DO KC=1,NCL3
                        EXPLT_io(IC,JC,KC) = RHSLLPHI_io(IC,JC,KC)
                    END DO
                END DO
            END DO
      
        ELSE
        END IF
      
        !============SET UP THE RHS WITH CONVECTION AND VISCOUS TERMS ================
        IF(visthemflg == visexplicit) THEN
            DO IC = NXI, NCL1_io
                DO JC = NYI, N2DO(MYID)
                    DO KC = 1, NCL3
                        RHSC = EXPLT_io (IC,JC,KC)    !CURRENT 
                        RHSL = EXPLT0_io(IC,JC,KC,IDR) ! LAST
                        RHS_io(IC,JC,KC)= (TGAM(NS)*RHSC + TROH(NS)* RHSL) * DT
                        EXPLT0_io(IC,JC,KC,IDR) = RHSC
                        ! if(IC==1 .and. JC==1 .and. KC==1 .and. myid==0) then
                        !     write(*,*)'test1', TGAM(NS), TROH(NS), RHSC, RHSL, DT, RHS_io(IC,JC,KC)
                        ! end if
                    END DO
                END DO
            END DO   
        ELSE IF(visthemflg == visimplicit) THEN
            COE2 = TALP(NS) * DT
            DO IC = NXI, NCL1_io
                DO JC = NYI, N2DO(MYID)
                    DO KC = 1, NCL3
                        RHSC = EXPLT_io (IC,JC,KC)    !CURRENT 
                        RHSL = EXPLT0_io(IC,JC,KC,IDR) ! LAST
                        RHS_io(IC,JC,KC)= (TGAM(NS)*RHSC + TROH(NS)* RHSL) * DT + &
                                           COE2 * RHS_io(IC,JC,KC) + &
                                           0.5_WP*COE2 * DIVU_io(IC,JC,KC)
                        EXPLT0_io(IC,JC,KC,IDR) = RHSC
                        
                    END DO
                END DO
            END DO   
        ELSE
        END IF  

        !===========SET UP PRESSURE GRADIENT=========================================
        PGM=0.0_WP
        IF (IDR.EQ.1) THEN 
            COE2 = TALP(NS) * DT*DXI
            DO IC=NXI,NCL1_io
                IM=IMV_io(IC)
                DO KC=1,NCL3
                    DO JC=NYI,N2DO(MYID)
                        IF(weightedpressure==1) THEN
                            PGM = ( (PR0_io(IC,JC,KC,2)-PR0_io(IM,JC,KC,2)) * (0.25_wp-pres_epslon) + &
                                    (PR0_io(IC,JC,KC,1)-PR0_io(IM,JC,KC,1)) * 0.5_wp + &
                                    (PR_io(IC,JC,KC)   -PR_io(IM,JC,KC))    * (0.25_wp+pres_epslon))*COE2
                        ELSE
                            PGM = ( PR_io(IC,JC,KC)-PR_io(IM,JC,KC) )*COE2
                        END IF 
                        RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC)-PGM
                        ! if(IC==1 .and. JC==1 .and. KC==1 .and. myid==0) then
                        !     write(*,*)'test2', PGM/TALP(NS)/DT, RHS_io(IC,JC,KC)
                        ! end if
                    ENDDO
                ENDDO
            ENDDO
            call wrt_3d_pt_debug(RHS_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), '', 'ConVisPX@af st', ITERG, NS) ! debug4chapsim2
        ELSE IF (IDR.EQ.2) THEN 
            
            DO JC=NYI,N2DO(MYID)
                JM=JLMV(JC)
                JJ=JCL2G(JC)
                COE2 = TALP(NS) * DT*DYCI(JJ)/RNDI1(jj)
                DO IC=NXI,NCL1_io
                    DO KC=1,NCL3
                        IF(weightedpressure==1) THEN
                            PGM = ( (PR0_io(IC,JC,KC,2)-PR0_io(IC,JM,KC,2)) * (0.25_wp-pres_epslon) + &
                                    (PR0_io(IC,JC,KC,1)-PR0_io(IC,JM,KC,1)) * 0.5_wp + &
                                    (PR_io(IC,JC,KC)   -PR_io(IC,JM,KC))    * (0.25_wp+pres_epslon))*COE2
                        ELSE
                            PGM = (PR_io(IC,JC,KC)-PR_io(IC,JM,KC))*COE2
                        END IF
                        RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC)-PGM
                    ENDDO
                ENDDO
            ENDDO    
            call wrt_3d_pt_debug(RHS_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), '', 'ConVisPY@af st', ITERG, NS) ! debug4chapsim2  
        ELSE IF (IDR.EQ.3) THEN
            COE2 = TALP(NS)* DT*DZI
            DO KC=1,NCL3
                KM=KMV(KC)
                DO JC=NYI,N2DO(MYID)
                    DO IC=NXI,NCL1_io
                        IF(weightedpressure==1) THEN
                            PGM = ( (PR0_io(IC,JC,KC,2)-PR0_io(IC,JC,KM,2)) * (0.25_wp-pres_epslon) + &
                                    (PR0_io(IC,JC,KC,1)-PR0_io(IC,JC,KM,1)) * 0.5_wp + &
                                    (PR_io(IC,JC,KC)   -PR_io(IC,JC,KM))    * (0.25_wp+pres_epslon))*COE2
                        ELSE
                            PGM = (PR_io(IC,JC,KC)-PR_io(IC,JC,KM))*COE2
                        END IF
                        RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC)-PGM
                    ENDDO
                ENDDO
            ENDDO  
            call wrt_3d_pt_debug(RHS_io(1:NCL1_io, 1:N2DO(myid), 1:NCL3), '', 'ConVisPZ@af st', ITERG, NS) ! debug4chapsim2
        ELSE       
        ENDIF
        
        !=========SET UP THE BUOYANCY FORCE=======================
        IF(GRAVDIR==IDR) THEN
     
            IF(GRAVDIR==1) THEN  
                COE2 = TALP(NS) * DT*F_A
                DO IC=NXI,NCL1_IO
                    IM = IMV_IO(IC)
                    DO KC=1,NCL3
                        DO JC=NYI,N2DO(MYID)
                            !DEN0 =  (DENSITY(IC,JC,KC) + DENSITY(IM,JC,KC))*XND2CL
                            RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC) + RHO_STG(IC,JC,KC,IDR)*COE2
                            !WRITE(*,*) MYID, IC, KC, JC, DEN0, RHO_STG(IC,JC,KC,IDR), DEN0-RHO_STG(IC,JC,KC,IDR) !test
                        ENDDO
                    ENDDO
                ENDDO
            END IF
        
            IF(GRAVDIR==2) THEN 
                 DO JC=NYI,N2DO(MYID)
                    JJ=JCL2G(JC)
                    JM = JLMV(JC)
                    COE2 = TALP(NS) * DT*F_A/RNDI1(JJ)
                    DO KC=1,NCL3
                        DO IC=NXI,NCL1_IO
                            !DEN0 =  YCL2ND_WFF(JJ) * DENSITY(IC,JC,KC) + &
                            !        YCL2ND_WFB(JJ) * DENSITY(IC,JM,KC)
                            RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC) + RHO_STG(IC,JC,KC,IDR)*COE2
                        ENDDO
                    ENDDO
                ENDDO
            END IF
        
            IF(GRAVDIR==3) THEN 
                COE2 = TALP(NS) * DT
                
                DO JC=NYI,N2DO(MYID)
                    JJ=JCL2G(JC)
                    COE2 = TALP(NS) * DT*F_A/RCCI1(JJ)
                    DO KC=1,NCL3
                        KM = KMV(KC)
                        DO IC=NXI,NCL1_IO
                            !DEN0 =  DENSITY(IC,JC,KC) + DENSITY(IC,JC,KM)
                            RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC) + RHO_STG(IC,JC,KC,IDR)*COE2
                        ENDDO
                    ENDDO
                ENDDO
            END IF
        
        END IF
        
        !====================flow drive terms (source terms) in periodic streamwise flow===========================
        IF ( (.NOT.TGFLOWFLG) .and. (IDR.EQ.NFLOW) .and. (ICASE.NE.IBOX3P)) THEN
            
            DPGRNS=0.0_WP
            !=============constant mass flow rate=============================
            IF(FLOWDRVTP==1) THEN         
                INTGRHSY=0.0_WP

                DO JC=1,N2DO(MYID)
                    JJ=JCL2G(JC)
                    DO IC=1,NCL1_io
                        DO KC=1,NCL3  
                            INTGRHSY=INTGRHSY+RHS_io(IC,JC,KC)/DYFI(JJ)/RCCI1(JJ)
                        ENDDO
                    ENDDO
                ENDDO
                
                CALL MPI_BARRIER(ICOMM,IERROR)
                CALL MPI_ALLREDUCE(INTGRHSY,INTGRHSY_WORK,1,MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)      
                IF(MYID.EQ.0) DPGRNS=INTGRHSY_WORK/VOLM_IO
            ELSE IF(FLOWDRVTP==2) THEN 
                !constant pressure gradient, dimensionless based on \delta and U_m
                !DPGRNS = -2.0_WP      ! DIMENSIONLESS BASED ON U_TAU
                COE2 = TALP(NS) * DT
                IF(MYID.EQ.0) DPGRNS = -0.50_WP*CFGV*COE2 !dimensionless based on \delta and U_m
            ELSE
            END IF  
            if(myid==0) then
                write(*,*) 'compensition', DPGRNS
                !write(*,*) 'rhsx:', RHS_io(:, 1, 1), RHS_io(:, 4, 4)
            end if
            CALL MPI_BCAST( DPGRNS, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )

            DO KC=1,NCL3
                DO IC=NXI,NCL1_io
                    DO JC=NYI,N2DO(MYID)
                        RHS_io(IC,JC,KC)=RHS_io(IC,JC,KC) - DPGRNS
                    ENDDO
                ENDDO
            ENDDO
        END IF
        

        !=====SAVE DATA BACK===========================================
        IF(visthemflg == visexplicit) THEN
            IF(IDR.EQ.1) THEN
             
                DO IC=NXI,NCL1_io
                    DO JC=NYI,N2DO(MYID)
                        DO KC=1,NCL3
                            Qtmp_io(IC,JC,KC) = RHS_io(IC,JC,KC)
                        END DO
                    END DO
                END DO
                !CALL DEBUG_WRT_LOCAL(Qtmp_io,1,N2DO(MYID),'conx') ! test
            ELSE IF(IDR .EQ. 2) THEN
          
                DO IC=NXI,NCL1_io
                    DO JC=NYI,N2DO(MYID)
                        DO KC=1,NCL3
                            DPH_io(IC,JC,KC) = RHS_io(IC,JC,KC)
                        END DO
                    END DO
                END DO
                !CALL DEBUG_WRT_LOCAL(DPH_IO,1,N2DO(MYID),'cony') ! test
            ELSE IF(IDR .EQ. 3) THEN
          
                DO IC=NXI,NCL1_io
                    DO JC=NYI,N2DO(MYID)
                        DO KC=1,NCL3
                            RHSLLPHI_io(IC,JC,KC) = RHS_io(IC,JC,KC)
                        END DO
                    END DO
                END DO
                !CALL DEBUG_WRT_LOCAL(RHSLLPHI_IO,1,N2DO(MYID),'conz') ! test
            ELSE
            END IF
        END IF
      
        RETURN
    END SUBROUTINE 

