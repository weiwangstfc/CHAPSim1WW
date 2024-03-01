    SUBROUTINE MASSFLUX_CALC_IO
        USE FLOW_INFO
        USE THERMAL_INFO
        USE MESH_INFO
        use init_info
        IMPLICIT NONE
        INTEGER(4) :: NXI, NYI,IDR
        INTEGER(4) :: IC, JC, KC, KS
      
        !=====================================X=============================
        IDR = 1
        NXI = 1
        IF(TGFLOWFLG) NXI = 2
        NYI = 1
        DO IC=NXI,NCL1_io
            DO JC=NYI,N2DO(MYID)
                DO KC=1,NCL3
                    G_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR) + Qtmp_io(IC,JC,KC)
                END DO
            END DO
        END DO
      
        !=====================================Z=============================
        IDR = 3
        NXI = 1
        NYI = 1
        DO IC=NXI,NCL1_io
            DO JC=NYI,N2DO(MYID)
                DO KC=1,NCL3
                    G_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR) + RHSLLPHI_io(IC,JC,KC)
                END DO
            END DO
        END DO
      
        !=====================================Y=============================
        IDR = 2
        NXI = 1
        NYI = 1
        IF(MYID==0 .and. ICASE.ne.IBOX3P) NYI = 2
        DO IC=NXI,NCL1_io
            DO JC=NYI,N2DO(MYID)
                DO KC=1,NCL3
                    G_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR) + DPH_io(IC,JC,KC)
                END DO
            END DO
        END DO
      
!==========below part is repeating interfaces===========
!        IF(MYID==0) THEN
!            IF(ICASE.EQ.IPIPEC) THEN
!                DO KC=1,NCL3
!                    KS =KSYM(KC)
!                    DO IC=NXI,NCL1_io
!                        G_io(IC,1,KC,IDR)=0.0_WP
!                        G_io(IC,0,KC,IDR)=G_io(IC,2,KS,IDR)
!                    END DO
!                END DO
!            ELSE
!                DO KC=1,NCL3
!                    DO IC=NXI,NCL1_io
!                        G_io(IC,1,KC,IDR)=0.0_WP
!                        G_io(IC,0,KC,IDR)=0.0_WP
!                    END DO
!                END DO
!            END IF
!        END IF
      
!        IF(MYID==NPSLV) THEN
!            DO KC=1,NCL3
!                DO IC=NXI,NCL1_io
!                    G_io(IC,N2DO(MYID)+1,KC,IDR)=0.0_WP
!                END DO
!            END DO
!        END IF
      
        !CALL INTFC_G_WALL
      
        RETURN
    END SUBROUTINE
      
!*********************************************************************************************************
    SUBROUTINE VELOCITY_CALC_io
        USE MESH_INFO
        USE THERMAL_INFO
        USE FLOW_INFO
        use init_info
        IMPLICIT NONE
        INTEGER(4) :: IDR, NYI, IC,JC,KC, IM, JJ, KM, JM!, KS
        REAL(wp)    :: RHO
      
        !=====================================X=============================
        IDR = 1
        NYI = 1
        DO IC=1,NCL1_io
            IM=IMV_io(IC)
            DO JC=NYI,N2DO(MYID) 
                DO KC=1,NCL3     
                    !RHO = ( DENSITY(IC,JC,KC) + DENSITY(IM,JC,KC) ) * XND2CL
                    Q_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)/RHO_STG(IC,JC,KC,IDR)
                ENDDO
            ENDDO
        ENDDO
       
        !=====================================Z=============================
        IDR = 3   
        NYI = 1  
        DO KC=1,NCL3
            KM=KMV(KC)
            DO JC=NYI,N2DO(MYID)
                DO IC=1,NCL1_io
                    !RHO = ( DENSITY(IC,JC,KC) + DENSITY(IC,JC,KM) ) * ZND2CL
                    Q_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)/RHO_STG(IC,JC,KC,IDR)
                ENDDO
            ENDDO
        ENDDO
      
      
        !=====================================Y=============================
        IDR = 2  
        NYI = 1
        IF(MYID==0 .and. ICASE.NE.IBOX3P) NYI = 2 
        DO JC=NYI,N2DO(MYID)
            JJ=JCL2G(JC)
            JM=JLMV(JC)    
            DO KC=1,NCL3
                DO IC=1,NCL1_io
                    !RHO = YCL2ND_WFF(JJ) * DENSITY(IC,JC,KC) + &
                    !      YCL2ND_WFB(JJ) * DENSITY(IC,JM,KC)
                    Q_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)/RHO_STG(IC,JC,KC,IDR)
                ENDDO
            ENDDO
        ENDDO 
        
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
        CALL BC_WALL_Q_io
      
        RETURN
    END SUBROUTINE
    
    
    
    !*********************************************************************************************************
    SUBROUTINE MU_Staggered
        USE FLOW_INFO
        USE THERMAL_INFO
        USE MESH_INFO
        use init_info
        IMPLICIT NONE

        INTEGER(4) :: IDR, NYI, IC,JC,KC, IM, JM, KM, JJ
        
        
        IF(thermlflg.ne.1) THEN
            MU_STG = 1.0_WP
            RETURN
        END IF
        
        !=====================================X-Y=(i',j',k)============================
        IDR = 1
        NYI = 1
        IF(MYID==0) NYI = 2
        DO IC=1,NCL1_io
            IM=IMV_io(IC)
            DO JC=NYI,N2DO(MYID) 
                JJ=JCL2G(JC)
                JM=JLMV(JC)
                DO KC=1,NCL3     
                    MU_STG(IC,JC,KC,IDR) = &
                        ( ( VISCOUSITY(IC,JC,KC)  + VISCOUSITY(IM,JC,KC) + &
                            VISCOUSITY0(IC,JC,KC) + VISCOUSITY0(IM,JC,KC)) * YCL2ND_WFF(JJ) + &
                          ( VISCOUSITY(IC,JM,KC)  + VISCOUSITY(IM,JM,KC) + &
                            VISCOUSITY0(IC,JM,KC) + VISCOUSITY0(IM,JM,KC)) * YCL2ND_WFB(JJ) ) * XND2CL*0.5_wp
                ENDDO
            ENDDO
        ENDDO
       
        !=====================================X-Z=(i',j,k')=============================
        ! Notice: the brackets in below equation will introduce differences in the order of 1E-14.
        IDR = 2   
        NYI = 1
        DO KC=1,NCL3
            KM=KMV(KC)
            DO JC=NYI,N2DO(MYID)
                DO IC=1,NCL1_io
                    IM=IMV_io(IC)
                    MU_STG(IC,JC,KC,IDR) = &
                          ( VISCOUSITY(IC,JC,KC)  + VISCOUSITY(IM,JC,KC)  + &
                            VISCOUSITY0(IC,JC,KC) + VISCOUSITY0(IM,JC,KC) + &
                            VISCOUSITY(IC,JC,KM)  + VISCOUSITY(IM,JC,KM)  + &
                            VISCOUSITY0(IC,JC,KM) + VISCOUSITY0(IM,JC,KM)  ) * XND2CL * ZND2CL*0.5_wp
                ENDDO
            ENDDO
        ENDDO
      
      
        !=====================================Y-Z=(i,j',k')============================
        IDR = 3  
        NYI = 1
        IF(MYID==0) NYI = 2
        DO JC=NYI,N2DO(MYID)
            JJ=JCL2G(JC)
            JM=JLMV(JC)    
            DO KC=1,NCL3
                KM=KMV(KC)
                DO IC=1,NCL1_io
                    MU_STG(IC,JC,KC,IDR) = &
                        (( VISCOUSITY(IC,JC,KM)  + VISCOUSITY(IC,JC,KC) + &
                           VISCOUSITY0(IC,JC,KM) + VISCOUSITY0(IC,JC,KC) ) * YCL2ND_WFF(JJ) + &
                         ( VISCOUSITY(IC,JM,KM)  + VISCOUSITY(IC,JM,KC)  + &
                           VISCOUSITY0(IC,JM,KM) + VISCOUSITY0(IC,JM,KC) ) * YCL2ND_WFB(JJ) ) * ZND2CL * 0.5_wp
                ENDDO
            ENDDO
        ENDDO 
        
        !=======================================================================================
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,MU_STG)
        
        IF(MYID.eq.0) THEN
            MU_STG(:,1,:,1)  = 0.5_wp*( VISCOUSITY0(:,0,:) + VISCOUSITY(:,0,:) )
            MU_STG(:,1,:,3)  = 0.5_wp*( VISCOUSITY0(:,0,:) + VISCOUSITY(:,0,:) )
        END IF
        
        IF(MYID.eq.NPSLV) THEN
            MU_STG(:,N2DO(MYID)+1,:,1)  = 0.5_wp*( VISCOUSITY0(:,N2DO(MYID)+1,:) + VISCOUSITY(:,N2DO(MYID)+1,:) )
            MU_STG(:,N2DO(MYID)+1,:,3)  = 0.5_wp*( VISCOUSITY0(:,N2DO(MYID)+1,:) + VISCOUSITY(:,N2DO(MYID)+1,:) )
        END IF
        
        !CALL DEBUG_WRT_LOCAL(MU_STG(:,:,:,1),0,N2DO(MYID)+1,'MUG1')
        !CALL DEBUG_WRT_LOCAL(MU_STG(:,:,:,2),0,N2DO(MYID)+1,'MUG2')
        !CALL DEBUG_WRT_LOCAL(MU_STG(:,:,:,3),0,N2DO(MYID)+1,'MUG3')
      
        RETURN
    END SUBROUTINE
    
    
    
    
    SUBROUTINE DENSITY_Staggered
        USE FLOW_INFO
        USE THERMAL_INFO
        USE MESH_INFO
        use init_info
        IMPLICIT NONE
        INTEGER(4) :: IDR, NYI, IC,JC,KC, IM, JJ, KM, JM
      
        IF(thermlflg.ne.1) THEN
            RHO_STG = 1.0_WP
            RETURN
        END IF
        
        !=====================================X=(i',j,k)==========================
        IDR = 1
        NYI = 1
        DO IC=1,NCL1_io
            IM=IMV_io(IC)
            DO JC=NYI,N2DO(MYID) 
                DO KC=1,NCL3     
                    RHO_STG(IC,JC,KC,IDR) = &
                         ( DENSITY(IC,JC,KC)  + DENSITY(IM,JC,KC) + &
                           DENSITY0(IC,JC,KC) + DENSITY0(IM,JC,KC) ) * XND2CL * 0.5_wp
                ENDDO
            ENDDO
        ENDDO
       
        !=====================================Z===X=(i,j,k')=======================
        IDR = 3   
        NYI = 1  
        DO KC=1,NCL3
            KM=KMV(KC)
            DO JC=NYI,N2DO(MYID)
                DO IC=1,NCL1_io
                    RHO_STG(IC,JC,KC,IDR) =&
                         ( DENSITY(IC,JC,KC)  + DENSITY(IC,JC,KM) &
                         + DENSITY0(IC,JC,KC) + DENSITY0(IC,JC,KM)) * ZND2CL * 0.5_wp
                ENDDO
            ENDDO
        ENDDO
      
      
        !=====================================Y====X=(i,j',k)=======================
        IDR = 2  
        NYI = 1
        IF(MYID==0) NYI = 2 
        DO JC=NYI,N2DO(MYID)
            JJ=JCL2G(JC)
            JM=JLMV(JC)    
            DO KC=1,NCL3
                DO IC=1,NCL1_io
                    RHO_STG(IC,JC,KC,IDR) = &
                         ( YCL2ND_WFF(JJ) * ( DENSITY(IC,JC,KC) + DENSITY0(IC,JC,KC) ) + &
                           YCL2ND_WFB(JJ) * ( DENSITY(IC,JM,KC) + DENSITY0(IC,JM,KC) ) ) *0.5_wp
                ENDDO
            ENDDO
        ENDDO 
        
        CALL INTFC_VARS3(1,NCL1_io,NCL1S,NCL1E,RHO_STG)
        
        IF(MYID.eq.0) THEN
            RHO_STG(:,1,:,2)             = 0.5_wp*( DENSITY0(:,0,:) +DENSITY(:,0,:) )
        END IF
        
        IF(MYID.eq.NPSLV) THEN
            RHO_STG(:,N2DO(MYID)+1,:,2)  = 0.5_wp*( DENSITY0(:,N2DO(MYID)+1,:) + DENSITY(:,N2DO(MYID)+1,:) )
        END IF
      
        RETURN
    END SUBROUTINE
    
    
    
    
    SUBROUTINE DENSITY_implicit_add
        USE FLOW_INFO
        USE THERMAL_INFO
        USE MESH_INFO
        use init_info
        IMPLICIT NONE
        INTEGER(4) :: IDR, NYI, IC,JC,KC, IM, JJ, KM, JM
        REAL(WP)   :: DEN_STG
             
        IF(thermlflg.ne.1) THEN
            DRHOI_STG = 0.0_WP
            RETURN
        END IF
        !=====================================X=============================
        IDR = 1
        NYI = 1
        DO IC=1,NCL1_io
            IM=IMV_io(IC)
            DO JC=NYI,N2DO(MYID) 
                DO KC=1,NCL3     
                    DEN_STG = ( DENSITY0(IC,JC,KC) + DENSITY0(IM,JC,KC) ) * XND2CL
                    DRHOI_STG(IC,JC,KC,IDR) = 1.0_WP/RHO_STG(IC,JC,KC,IDR) - 1.0_WP/DEN_STG
                ENDDO
            ENDDO
        ENDDO
       
        !=====================================Z=============================
        IDR = 3   
        NYI = 1  
        DO KC=1,NCL3
            KM=KMV(KC)
            DO JC=NYI,N2DO(MYID)
                DO IC=1,NCL1_io
                    DEN_STG = ( DENSITY0(IC,JC,KC) + DENSITY0(IC,JC,KM) ) * ZND2CL
                    DRHOI_STG(IC,JC,KC,IDR) = 1.0_WP/RHO_STG(IC,JC,KC,IDR) - 1.0_WP/DEN_STG
                ENDDO
            ENDDO
        ENDDO
      
      
        !=====================================Y=============================
        IDR = 2  
        NYI = 1
        IF(MYID==0) NYI = 2 
        DO JC=NYI,N2DO(MYID)
            JJ=JCL2G(JC)
            JM=JLMV(JC)    
            DO KC=1,NCL3
                DO IC=1,NCL1_io
                    DEN_STG = YCL2ND_WFF(JJ) * DENSITY0(IC,JC,KC) + &
                              YCL2ND_WFB(JJ) * DENSITY0(IC,JM,KC)
                    DRHOI_STG(IC,JC,KC,IDR) = 1.0_WP/RHO_STG(IC,JC,KC,IDR) - 1.0_WP/DEN_STG
                ENDDO
            ENDDO
        ENDDO 
        
        CALL INTFC_VARS3(1,NCL1_io,NCL1S,NCL1E,DRHOI_STG)
        
        
        IF(MYID.eq.0) THEN
            DO IC=1,NCL1_io
                DRHOI_STG(IC,1,:,2)  = 0.0_WP
            END DO
        END IF
        
        IF(MYID.eq.NPSLV) THEN
            DO IC=1,NCL1_io
                DRHOI_STG(IC,N2DO(MYID)+1,:,2)  = 0.0_WP
            END DO
        END IF
        
        
        !test
!        do jc=1,n2do(myid)
!            do ic=1,ncl1_io
!                do kc=1,ncl3
!                    if(  dabs(DRHOI_STG(IC,jc,kc,1)) .gt. 1.0E-10_wp) &!
!                    write(*,*) myid,jc,ic,kc, DRHOI_STG(IC,jc,kc,1:ndv)
!                end do
!            end do
!        end do
        !CALL DEBUG_WRT_LOCAL(DRHOI_STG,0,N2DO(MYID)+1,'dstg')
      
        RETURN
    END SUBROUTINE
