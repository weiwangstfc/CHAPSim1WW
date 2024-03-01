 !**********************************************************************
    SUBROUTINE BC_COUTLET_MOM_RK3(NS)
        use mesh_info
        use flow_info
        use init_info
        use thermal_info
        IMPLICIT NONE
      
        INTEGER(4),INTENT(IN) :: NS
      
        INTEGER(4) :: J, K, IDR, N2I
        REAL(WP)    :: BC_RHS
        REAL(WP)    :: BC_CONV
        real(WP)    :: RHO, COE1
        
        IF(.NOT. TGFLOWFLG) RETURN
        
        !====method one=====================================================
        ! outlet in x direction. If in other directions, the cylindrical coordinates will introduce other terms, see eq.4.31 of Pierce thesis.
        CALL CONVCTION_OUTLET_U
      
        COE1 = DXI*U_OUTLET 
        DO IDR=1, NDV
      
            IF( (MYID.EQ.0) .AND. (IDR .EQ. 2)) then
                N2I = 2
            else 
                N2I = 1
            end if
          
            DO K=1, NCL3
                DO J=N2I, N2DO(myid)
                    BC_CONV=-( G_io(NCL1_io+1,J,K,IDR)- G_io(NCL1_io,  J,K,IDR) )*COE1   
                    BC_RHS = TGAM(NS)*BC_CONV + TROH(NS)*BC_CONV0(J,K,IDR)
                    BC_CONV0(J,K,IDR)=BC_CONV
                    
                    G_io(NCL1_io+1,J,K,IDR) = G_io(NCL1_io+1,J,K,IDR) + DT * BC_RHS
                                              
                    RHO = ( DENSITY(NCL1_io,J,K) + DENSITY(NCL1_io+1,J,K) ) * XND2CL
                    Q_io(NCL1_io+1,J,K,IDR) = G_io(NCL1_io+1,J,K,IDR)/RHO
                END DO
            END DO
          
        END DO
      
        
        CALL CHK_MassConsv_io
        ! add a limiter to the variations of the outlet velocity============
        !IF(dabs(CHK_Mass_CONSV0_SCALING-1.0_wp).lt.0.5_wp) then 
        DO J=1,N2DO(MYID) 
            DO K=1,NCL3
                G_io(NCL1_io+1,J,K,NFLOW) = G_io(NCL1_io+1,J,K,NFLOW) * (CHK_Mass_CONSV0_SCALING-1.0_wp) &
                                          + G_io(NCL1_io+1,J,K,NFLOW)
                RHO = ( DENSITY(NCL1_io,J,K) + DENSITY(NCL1_io+1,J,K) ) * XND2CL
                Q_io(NCL1_io+1,J,K,NFLOW) = G_io(NCL1_io+1,J,K,NFLOW)/RHO
            END DO
        END DO
        !end if
      
        !=========================update============================================
        CALL INTFC_VARS3(NCL1_io+1,NCL1_io+1,NCL1S,NCL1E,G_io)
        CALL BC_WALL_G_io
        CALL INTFC_VARS3(NCL1_io+1,NCL1_io+1,NCL1S,NCL1E,Q_io)
        CALL BC_WALL_Q_io
      
!        !=====method two, sometimes diverged==============================
!        ! use the convergence critria to calculate the Nflow streamwise velocity always introduces divergence after hundreds of steps.
!        CALL CBC_Gx
!        CALL CONVCTION_OUTLET_U
      
!        DO IDR=2, NDV
      
!            IF( (MYID.EQ.0) .AND. (IDR .EQ. 2)) then
!                N2I = 2
!            else 
!                N2I = 1
!            end if
          
!            DO K=1, NCL3
!                DO J=N2I, N2DO(myid)
!                    BC_CONV=-1.0_WP*( G_io(NCL1_io+1,J,K,IDR)- &
!                                      G_io(NCL1_io,  J,K,IDR) )*DXI*U_OUTLET 
!                    BC_RHS = TGAM(NS)*BC_CONV + TROH(NS)*BC_CONV0(J,K,IDR)
!                    BC_CONV0(J,K,IDR)=BC_CONV
                    
!                    G_io(NCL1_io+1,J,K,IDR) = G_io(NCL1_io+1,J,K,IDR) + &
!                                              DT * BC_RHS
!                    RHO = ( DENSITY(NCL1_io,J,K) + DENSITY(NCL1_io+1,J,K) ) * ZND2CL
!                    Q_io(NCL1_io+1,J,K,IDR) = G_io(NCL1_io+1,J,K,IDR)/RHO
!                END DO
!            END DO
          
!        END DO
      
      
        RETURN
    END SUBROUTINE 

!**********************************************************************
 !**********************************************************************
    SUBROUTINE CONVCTION_OUTLET_U
!     Below is based on constant density in the domian.
        use thermal_info
        use mpi_info
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        !REAL(WP) :: FluxINTG_INL, FluxINTG_INL_WORK
        !REAL(WP) :: FluxINTG_OUT, FluxINTG_OUT_WORK
        !REAL(WP) :: MASS_RATE_INTG, MASS_RATE
        !REAL(WP) :: MASS_RATE_INTG_WORK
        !INTEGER(4) :: IC, JC, KC, JJ, I
      
        INTEGER(4) ::  J, K
        REAL(WP) :: Uxmax, Uxmin, Uxmax1, Uxmin1
        REAL(WP) :: UCC
      
!>================Method One========================================    
!        FluxINTG_INL   = 0.0_WP
!        FluxINTG_OUT   = 0.0_WP   
!        MASS_RATE_INTG = 0.0_WP         
!        DO JC=1,N2DO(myid)
!            JJ=JCL2G(JC)
!            DO KC=1,NCL3
!                FluxINTG_INL   =FluxINTG_INL   + G_io(1,        JC,KC,NFLOW)/DYFI(JJ)
!                FluxINTG_OUT   =FluxINTG_OUT   + 0.5_wp*( DENSITY(NCL1_IO,JC,KC) + DENSITY(NCL1_IO+1,JC,KC) ) /DYFI(JJ)
                
!                DO IC = 1, NCL1_IO
!                    MASS_RATE      = ( DENSITY(IC,JC,KC) - DENSITY0(IC,JC,KC) )/DYFI(JJ)
!                    MASS_RATE_INTG =  MASS_RATE + MASS_RATE_INTG
!                END DO
!            ENDDO
!        ENDDO
!        FluxINTG_INL   = FluxINTG_INL/DZI
!        FluxINTG_OUT   = FluxINTG_OUT/DZI
!        MASS_RATE_INTG = MASS_RATE_INTG/DZI/DXI/DT
        
!        CALL MPI_ALLREDUCE(FluxINTG_INL,FluxINTG_INL_WORK,1,MPI_DOUBLE_PRECISION,&
!                             MPI_SUM, ICOMM, IERROR)  
!        CALL MPI_ALLREDUCE(FluxINTG_OUT,FluxINTG_OUT_WORK,1,MPI_DOUBLE_PRECISION,&
!                             MPI_SUM, ICOMM, IERROR)  
!        CALL MPI_ALLREDUCE(MASS_RATE_INTG,MASS_RATE_INTG_WORK,1,MPI_DOUBLE_PRECISION,&
!                             MPI_SUM, ICOMM, IERROR)  
                             
!        U_OUTLET = (FluxINTG_INL_WORK + MASS_RATE_INTG_WORK)/FluxINTG_OUT_WORK
          
          !write(*,*) 'U_OUTLET1',FluxINTG_INL_WORK, FluxINTG_OUT_WORK, U_OUTLET !test
     !END IF
!!!>=============Method Two==============================================
!     IF(FLOWDRVTP==2) THEN
        uxmax=-1.0E20_WP
        uxmin= 1.0E20_WP
        do k=1,NCL3
            do j=1,N2DO(MYID)
                UCC = ( Q_io(NCL1_io,J,K,NFLOW) + Q_io(NCL1_io+1,J,K,NFLOW) )*0.5_WP
                if (UCC.gt.uxmax) uxmax=UCC
                if (UCC.lt.uxmin) uxmin=UCC
            enddo
        enddo
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        call MPI_ALLREDUCE(uxmax,uxmax1,1,MPI_DOUBLE_PRECISION,MPI_MAX,ICOMM, IERROR)
        call MPI_ALLREDUCE(uxmin,uxmin1,1,MPI_DOUBLE_PRECISION,MPI_MIN,ICOMM, IERROR)

        U_OUTLET = 0.5_wp*(uxmax1+uxmin1)

        !write(*,*) 'U_OUTLET1',U_OUTLET !test
!      END IF
      
      
        RETURN
    END SUBROUTINE
      
      
!**********************************************************************
      SUBROUTINE BC_CBC_TDMA(NS)
!     Check N2I = 1 or 2
      use mesh_info
      use flow_info
      use init_info
      IMPLICIT NONE
      
      INTEGER(4),INTENT(IN) :: NS
      
      INTEGER(4) :: J, JJ, JM, JC, JP
      INTEGER(4) :: K,     KM, KC, KP
      INTEGER(4) :: IDR
      REAL(WP)    :: AMJV, ACJV, APJV
      REAL(WP)    :: AMKV, ACKV, APKV
      REAL(WP)    :: COE1
      
      
      COE1=0.50_WP*TALP(NS)*DT*CVISC
      
      CALL BC_INTFC_CBC(3)
!>    @note : With BC_INTFC_CBC(3) all known, and tri-matrix, 
!>            TRANSP23 is not necessary here.
      DO IDR = 1, NDV 
         DO K=1, NCL3
            DO JC=1,N2DO(myid)
               JM = JLMV(JC)
               JP = JLPV(JC)
               JJ = JCL2G(JC)
               ACJV = 1.0_WP-COE1*ACVR(JJ,IDR)
               APJV =       -COE1*APVR(JJ,IDR)
               AMJV =       -COE1*AMVR(JJ,IDR)
               BC_TDMA(2, JC, K, IDR) = AMJV*BC_TDMA(3,JM,K,IDR)+ &
                                        ACJV*BC_TDMA(3,JC,K,IDR)+ &
                                        APJV*BC_TDMA(3,JP,K,IDR)
            END DO
         END DO
      END DO
       
      !CALL BC_INTFC_CBC(2) not necessary

      DO IDR = 1, NDV 
         DO J=1,N2DO(myid)
            DO KC=1, NCL3
               KM = KMV(KC)
               KP = KPV(KC)
               ACKV = 1.0_WP-COE1*(-2.0_WP*DZQI)
               APKV =     -COE1*DZQI
               AMKV =     -COE1*DZQI
               BC_TDMA(1, J, KC, IDR) = AMKV*BC_TDMA(2,J,KM,IDR)+ &
                                        ACKV*BC_TDMA(2,J,KC,IDR)+ &
                                        APKV*BC_TDMA(2,J,KP,IDR)
            END DO
         END DO
      END DO
      
      !CALL BC_INTFC_CBC(1) not necessary
      

      RETURN
      END SUBROUTINE

!**********************************************************************

!************************************************************************
    SUBROUTINE BC_INTFC_CBC(FASTEP)
        use mpi_info
        use mesh_info
        use flow_info
        use init_info
        IMPLICIT NONE

        INTEGER(4),INTENT(IN)  :: FASTEP
       
        INTEGER(4)  :: K
        INTEGER(4)  :: IDR
        INTEGER(4)  :: ITAG
        INTEGER(4)  :: IDESF
        INTEGER(4)  :: IDESB
        INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
        INTEGER(4)  :: NSZ
       
        REAL(WP)    :: BSEN_F_io(NCL1_IO * NCL3,NDV)
        REAL(WP)    :: BSEN_L_io(NCL1_IO * NCL3,NDV)
        REAL(WP)    :: BREC_F_io(NCL1_IO * NCL3,NDV)
        REAL(WP)    :: BREC_L_io(NCL1_IO * NCL3,NDV)

        BSEN_F_io = 0.0_WP
        BSEN_L_io = 0.0_WP
        BREC_F_io = 0.0_WP
        BREC_L_io = 0.0_WP
        
        NSZ      = NCL1_IO * NCL3*NDV 
        
                       
        DO IDR=1,NDV
            DO K=1,NCL3
                BSEN_F_io(K,IDR)=BC_TDMA(FASTEP,1,         K,IDR) !y=local Floor
                BSEN_L_io(K,IDR)=BC_TDMA(FASTEP,N2DO(MYID),K,IDR) !y=local ROOF
            ENDDO
        ENDDO

       
        ITAG=0      
!>     @note rank=2,3,4...SIZE-1 (no b.c.)***************************************      
        IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
            IDESF = MYID+1    !next myid, forward
            IDESB = MYID-1    !last myid, backward


            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
        
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)
    

!>           @note Constructing p,u,v,w of ghost cells y_id=0 and N2NO+1, received from adjacent partitions.      
         
            DO IDR=1,NDV
                DO K=1,NCL3
                    BC_TDMA(FASTEP,0,           K,IDR)=BREC_L_io(K,IDR)
                    BC_TDMA(FASTEP,N2DO(MYID)+1,K,IDR)=BREC_F_io(K,IDR)
                ENDDO
            ENDDO
 
        ENDIF
       
!>     @note rank=SIZE-1 (top wall)***************************************
        IF(MYID.EQ.NPSLV) THEN
            IDESF = 0  
            IDESB = MYID-1   
          
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION, IDESB, ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION, IDESF, ITAG,ICOMM,TRC_STS,IERROR)
                
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                    MPI_DOUBLE_PRECISION, IDESF, ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION, IDESB, ITAG,ICOMM,TRC_STS,IERROR)
          
            DO IDR=1,3
                DO K=1,NCL3
                    BC_TDMA(FASTEP,0, K,IDR)=BREC_L_io(K,IDR)
                END DO
            ENDDO
          
          
            DO K=1,NCL3
                 BC_TDMA(FASTEP,N2DO(myid)+1, K,1) = -BC_TDMA(FASTEP,N2DO(myid), K,1)
                 BC_TDMA(FASTEP,N2DO(myid)+1, K,3) = -BC_TDMA(FASTEP,N2DO(myid), K,3)
                 BC_TDMA(FASTEP,N2DO(myid)+1, K,2) = 0.0_WP     
            ENDDO               
          
        ENDIF

!>       @note rank=0 (bottom wall)***************************************
        IF (MYID.EQ.0) THEN
            IDESF = MYID+1
            IDESB = NPSLV
         
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)
                
!>           @note bottom wall velocity B.C.
            DO IDR=1,3
                DO K=1,NCL3 
                    BC_TDMA(FASTEP,N2DO(myid)+1, K,IDR)=BREC_F_io(K,IDR)
                END DO
            ENDDO          
          
            DO K=1,NCL3                
                BC_TDMA(FASTEP,0,K,1)=-BC_TDMA(FASTEP,1,K,1)
                BC_TDMA(FASTEP,0,K,3)=-BC_TDMA(FASTEP,1,K,3)
                BC_TDMA(FASTEP,1,K,2)=0.0_WP
                BC_TDMA(FASTEP,0,K,2)=0.0_WP
            ENDDO
          
        ENDIF
       
        RETURN
      END SUBROUTINE
      
!***********************************************************************  

!************************
      SUBROUTINE VELOUPDT_CBC_U
      use init_info
      use mesh_info
      use flow_info
      IMPLICIT NONE 
      
      INTEGER(4) :: J,K,JP,JJ,KP
      
      DO J = 1, N2DO(MYID)
         JP=JLPV(J)
         JJ=JCL2G(J)
         DO K=1, NCL3
             KP=KPV(K)
             BC_TDMA(3,J,K,1) = Q_io(NCL1_io+1, J, K, 1)
             
             Q_io(NCL1_io+1, J, K, 1)= -( (Q_io(NCL1_io,JP,K,2)-Q_io(NCL1_io,J,K,2))*DYFI(JJ) + &
                                          (Q_io(NCL1_io,J,KP,3)-Q_io(NCL1_io,J,K,3))*DZI ) /DXI &
                                         + Q_io(NCL1_io,J,K,1)
         END DO
      END DO
      
      !CALL MassConservationCheck(coeffou)
      !write(*,*) 'coeffou',coeffou
      
      RETURN
      END SUBROUTINE


!************************
    SUBROUTINE CBC_Gx
        use init_info
        use mesh_info
        use flow_info
        use thermal_info
        IMPLICIT NONE 
      
        INTEGER(4) :: IC,JC,KC,JP,JJ,KP, IP
        REAL(WP)    :: RHO, DTRHO, DIVY, DIVZ
        
        IC = NCL1_IO
        IP = IC + 1
        DO JC = 1, N2DO(MYID)
            JP=JLPV(JC)
            JJ=JCL2G(JC)
            DO KC=1, NCL3
                KP=KPV(KC)
             
                !==========d(\rho u)/dx at (i,j,k)===========================
                !DIVX  = ( G_IO(IP,JC,KC,1) - G_IO(IC,JC,KC,1) ) * DXI
               
                !==========d(\rho v)/dy at (i,j,k)===========================
                DIVY  = ( G_io(IC,JP,KC,2) - G_io(IC,JC,KC,2) ) * DYFI(JJ)
               
                !==========d(\rho w)/dz at (i,j,k)===========================
                DIVZ  = ( G_io(IC,JC,KP,3) - G_io(IC,JC,KC,3) ) * DZI 
               
                !==========D \RHO/ DT========================================
                !DTRHO = ( DENSITY(IC,JC,KC) - DENSITYp(IC,JC,KC) ) / DT
                DTRHO = DrhoDtP(IC,JC,KC)
               
                G_IO(IP,JC,KC,1) = G_IO(IC,JC,KC,1)-(DTRHO+DIVY+DIVZ)*DX
                
                RHO = ( DENSITY(IC,JC,KC) + DENSITY(IP,JC,KC) ) * ZND2CL
                Q_IO(IP,JC,KC,1) = G_IO(IP,JC,KC,1)/RHO
         END DO
      END DO
      
      !CALL MassConservationCheck(coeffou)
      !write(*,*) 'coeffou',coeffou
      
      RETURN
      END SUBROUTINE
