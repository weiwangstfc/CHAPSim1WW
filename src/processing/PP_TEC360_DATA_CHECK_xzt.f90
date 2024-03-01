!=============================TG+IO==================================
    MODULE TEC360_INFO
        use init_info
        use flow_info
        use postprocess_info
        use mesh_info
        use wrt_info
        use thermal_info
        !========================TG=========================
        REAL(WP),ALLOCATABLE    :: U_F0_tg(:,:,:,:)
        REAL(WP),ALLOCATABLE    :: U1xztL_F0_tg(:,:)
        REAL(WP),ALLOCATABLE    :: U1xztL_INTP_tg(:,:)
        
        REAL(WP),ALLOCATABLE    :: U_INTP_tg(:, :, :, :)
        real(WP), allocatable   :: uprime_tg(:,:,:,:)
        
        real(WP), allocatable   :: Qcr_tg(:,:,:)
        real(WP), allocatable   :: vor_tg(:,:,:,:)
        real(WP), allocatable   :: delta_tg(:,:,:)
        real(WP), allocatable   :: lambda2_tg(:,:,:)
        real(WP), allocatable    :: swirlstrength_tg(:,:,:,:)
        
        !========================IO=========================
        REAL(WP),ALLOCATABLE    :: U_F0_io(:,:,:,:)
        REAL(WP),ALLOCATABLE    :: G_F0_io(:,:,:,:)
        REAL(WP),ALLOCATABLE    :: H_F0_io(:,:,:)
        REAL(WP),ALLOCATABLE    :: T_F0_io(:,:,:)
        REAL(WP),ALLOCATABLE    :: D_F0_io(:,:,:)
        REAL(WP),ALLOCATABLE    :: M_F0_io(:,:,:)
        

        REAL(WP),ALLOCATABLE    :: U1xztL_F0_io(:,:)
        REAL(WP),ALLOCATABLE    :: G1xztL_F0_io(:,:)
        REAL(WP),ALLOCATABLE    :: H1xztL_F0_io(:)
        REAL(WP),ALLOCATABLE    :: T1xztL_F0_io(:)
        REAL(WP),ALLOCATABLE    :: D1xztL_F0_io(:)
        REAL(WP),ALLOCATABLE    :: M1xztL_F0_io(:)
        
        REAL(WP),ALLOCATABLE    :: U1ztL_F0_io(:,:,:)
        REAL(WP),ALLOCATABLE    :: G1ztL_F0_io(:,:,:)
        REAL(WP),ALLOCATABLE    :: H1ztL_F0_io(:,:)
        REAL(WP),ALLOCATABLE    :: T1ztL_F0_io(:,:)
        REAL(WP),ALLOCATABLE    :: D1ztL_F0_io(:,:)
        REAL(WP),ALLOCATABLE    :: M1ztL_F0_io(:,:)
        
        
        REAL(WP),ALLOCATABLE     :: U_INTP_io(:, :, :, :)
        REAL(WP),ALLOCATABLE     :: G_INTP_io(:, :, :, :)
        REAL(WP),ALLOCATABLE     :: H_INTP_io(:, :, :)
        REAL(WP),ALLOCATABLE     :: T_INTP_io(:, :, :)
        REAL(WP),ALLOCATABLE     :: D_INTP_io(:, :, :)
        REAL(WP),ALLOCATABLE     :: M_INTP_io(:, :, :)
        
        
        REAL(WP),ALLOCATABLE    :: U1xztL_INTP_io(:,:)
        REAL(WP),ALLOCATABLE    :: G1xztL_INTP_io(:,:)
        REAL(WP),ALLOCATABLE    :: H1xztL_INTP_io(:)
        REAL(WP),ALLOCATABLE    :: T1xztL_INTP_io(:)
        REAL(WP),ALLOCATABLE    :: D1xztL_INTP_io(:)
        REAL(WP),ALLOCATABLE    :: M1xztL_INTP_io(:)
        
        REAL(WP),ALLOCATABLE    :: U1ztL_INTP_io(:,:,:)
        REAL(WP),ALLOCATABLE    :: G1ztL_INTP_io(:,:,:)
        REAL(WP),ALLOCATABLE    :: H1ztL_INTP_io(:,:)
        REAL(WP),ALLOCATABLE    :: T1ztL_INTP_io(:,:)
        REAL(WP),ALLOCATABLE    :: D1ztL_INTP_io(:,:)
        REAL(WP),ALLOCATABLE    :: M1ztL_INTP_io(:,:)
        
        real(WP), allocatable    :: uprime_io(:,:,:,:)
        real(WP), allocatable    :: gprime_io(:,:,:,:)
        real(WP), allocatable    :: Tprime_io(:,:,:)
        real(WP), allocatable    :: dprime_io(:,:,:)
        real(WP), allocatable    :: mprime_io(:,:,:)
        real(WP), allocatable    :: udprime_io(:,:,:,:)
        
        real(WP), allocatable    :: Qcr_io(:,:,:)
        real(WP), allocatable    :: vor_io(:,:,:,:)
        real(WP), allocatable    :: delta_io(:,:,:)
        real(WP), allocatable    :: lambda2_io(:,:,:)
        real(WP), allocatable    :: swirlstrength_io(:,:,:,:)
      
        INTEGER(4) :: NI
        !INTEGER(4) :: NJ
        INTEGER(4) :: NK 
        INTEGER(4),ALLOCATABLE :: IID(:)
        !INTEGER(4),ALLOCATABLE :: JID(:)
        INTEGER(4),ALLOCATABLE :: KID(:)
        INTEGER(4) :: NCOUNT(5) = 0
      
        
      
    END MODULE TEC360_INFO

    SUBROUTINE DEALLO_TEC
        use TEC360_INFO
        IMPLICIT NONE
        
        !========================TG=========================
        if ( allocated(U_F0_tg) )        deallocate(U_F0_tg)
        if ( allocated(U1xztL_F0_tg) )   deallocate(U1xztL_F0_tg)
        if ( allocated(U1xztL_INTP_tg) ) deallocate(U1xztL_INTP_tg)
        
        
        if ( allocated(uprime_tg) )  deallocate(uprime_tg)
        if ( allocated(U_INTP_tg) )  deallocate(U_INTP_tg)
        
        if ( allocated(Qcr_tg) )   deallocate(Qcr_tg)
        if ( allocated(vor_tg) )   deallocate(vor_tg)
        if ( allocated(delta_tg) ) deallocate(delta_tg)
        if ( allocated(lambda2_tg) ) deallocate(lambda2_tg)
        if ( allocated(swirlstrength_tg) )    deallocate(swirlstrength_tg)
        
        
        !========================IO=========================
        if ( allocated(U_F0_io) ) deallocate(U_F0_io)
        if ( allocated(G_F0_io) ) deallocate(G_F0_io)
        if ( allocated(H_F0_io) ) deallocate(H_F0_io)
        if ( allocated(T_F0_io) ) deallocate(T_F0_io)
        if ( allocated(D_F0_io) ) deallocate(D_F0_io)
        if ( allocated(D_F0_io) ) deallocate(M_F0_io)
        
        

        if ( allocated(U1xztL_F0_io) ) deallocate(U1xztL_F0_io)
        if ( allocated(G1xztL_F0_io) ) deallocate(G1xztL_F0_io)
        if ( allocated(T1xztL_F0_io) ) deallocate(H1xztL_F0_io)
        if ( allocated(T1xztL_F0_io) ) deallocate(T1xztL_F0_io)
        if ( allocated(D1xztL_F0_io) ) deallocate(D1xztL_F0_io)
        if ( allocated(M1xztL_F0_io) ) deallocate(M1xztL_F0_io)
        
        if ( allocated(U1ztL_F0_io) )  deallocate(U1ztL_F0_io)
        if ( allocated(G1ztL_F0_io) )  deallocate(G1ztL_F0_io)
        if ( allocated(T1ztL_F0_io) )  deallocate(H1ztL_F0_io)
        if ( allocated(T1ztL_F0_io) )  deallocate(T1ztL_F0_io)
        if ( allocated(D1ztL_F0_io) )  deallocate(D1ztL_F0_io)
        if ( allocated(M1ztL_F0_io) )  deallocate(M1ztL_F0_io)
        
        
        if ( allocated(U_INTP_io) ) deallocate(U_INTP_io)
        if ( allocated(G_INTP_io) ) deallocate(G_INTP_io)
        if ( allocated(H_INTP_io) ) deallocate(H_INTP_io)
        if ( allocated(T_INTP_io) ) deallocate(T_INTP_io)
        if ( allocated(D_INTP_io) ) deallocate(D_INTP_io)
        if ( allocated(D_INTP_io) ) deallocate(M_INTP_io)
        
        
        if ( allocated(U1xztL_INTP_io) ) deallocate(U1xztL_INTP_io)
        if ( allocated(G1xztL_INTP_io) ) deallocate(G1xztL_INTP_io)
        if ( allocated(T1xztL_INTP_io) ) deallocate(H1xztL_INTP_io)
        if ( allocated(T1xztL_INTP_io) ) deallocate(T1xztL_INTP_io)
        if ( allocated(D1xztL_INTP_io) ) deallocate(D1xztL_INTP_io)
        if ( allocated(M1xztL_INTP_io) ) deallocate(M1xztL_INTP_io)
        
        if ( allocated(U1ztL_INTP_io) )  deallocate(U1ztL_INTP_io)
        if ( allocated(G1ztL_INTP_io) )  deallocate(G1ztL_INTP_io)
        if ( allocated(T1ztL_INTP_io) )  deallocate(H1ztL_INTP_io)
        if ( allocated(T1ztL_INTP_io) )  deallocate(T1ztL_INTP_io)
        if ( allocated(D1ztL_INTP_io) )  deallocate(D1ztL_INTP_io)
        if ( allocated(M1ztL_INTP_io) )  deallocate(M1ztL_INTP_io)
        
        if ( allocated(uprime_io) )   deallocate(uprime_io)
        if ( allocated(gprime_io) )   deallocate(gprime_io)
        if ( allocated(Tprime_io) )   deallocate(Tprime_io)
        if ( allocated(dprime_io) )   deallocate(dprime_io)
        if ( allocated(mprime_io) )   deallocate(mprime_io)
        if ( allocated(udprime_io) )  deallocate(udprime_io)
        
        if ( allocated(Qcr_io) )        deallocate(Qcr_io)
        if ( allocated(vor_io) )        deallocate(vor_io)
        if ( allocated(delta_io) )      deallocate(delta_io)
        if ( allocated(lambda2_io) )    deallocate(lambda2_io)
        if ( allocated(swirlstrength_io) )    deallocate(swirlstrength_io)
        
        if ( allocated(IID) ) deallocate(IID)
        !if ( allocated(JID) ) deallocate(JID)
        if ( allocated(KID) ) deallocate(KID)
      

        RETURN
        
    END SUBROUTINE
    
    
!**********************WRITE INITIAL CALCUATED RESULTS OUT*******************************************

    SUBROUTINE CALL_TEC360
        use TEC360_INFO
        IMPLICIT NONE

        INTEGER(4) :: N
        
        !============WRITE INITIAL CALCUATED RESULTS OUT======================================
        IF(TGFLOWflg .AND. IOFLOWflg) THEN
            
            CALL GATHERING_TG
            CALL GATHERING_io
            
            IF(MYID==0) THEN
                NI = 4
                !NJ = 4
                NK = 2
                ALLOCATE( IID(NI) )
                !ALLOCATE( JID(NJ) )
                ALLOCATE( KID(NK) )
                
                CALL TEC360_EXTRPLAT2ND_MASTER_TG
                CALL TEC360_EXTRPLAT2ND_MASTER_IO
            END IF 

        ELSE
            IF(MYID==0) THEN
                NI = 1
                !NJ = 8
                NK = 1
                ALLOCATE( IID(NI) )
                !ALLOCATE( JID(NJ) )
                ALLOCATE( KID(NK) )
            END IF
            
            IF (TGFLOWflg) THEN
                CALL GATHERING_TG
                IF(MYID==0) CALL TEC360_EXTRPLAT2ND_MASTER_TG 
            END IF
             
            IF(IOFLOWflg) THEN
                CALL GATHERING_IO
                IF(MYID==0) CALL TEC360_EXTRPLAT2ND_MASTER_IO
            END IF
        END IF
        
        IF(MYID==0) THEN
            CALL TEC360_INSTANT_VORTEX_CRITERIA
            CALL TEC360_INSTANT_UPRIME
            
            IF(PPROCESSONLY==1 .and. ppinst == 1) THEN
                CALL TEC360_ALL_NODES
            END IF
    
            !===================WRITE X-SLICE================================== 
            
            DO N=1,NI
                IID(N) = N*(NND1_IO+1)/(2*NI) + 1
                CALL TEC360_XSLICE(N)
            END DO
            
            !===================WRITE Y-SLICE================================== 
            !DO N=1,NJ
                !JID(N) = N*(NND2+1)/(2*NJ) + 1
            DO N=1, MGRID
                CALL TEC360_YSLICE(N)
            END DO
            !===================WRITE Z-SLICE================================== 
            DO N=1,NK
                KID(N) = N*(NND3+1)/(2*NK) + 1
                CALL TEC360_ZSLICE(N)
            END DO
        END IF
    
        CALL DEALLO_TEC
        
        RETURN
    END SUBROUTINE
    
    
!**************************************************************************************************
    SUBROUTINE GATHERING_TG
        use TEC360_INFO
        IMPLICIT NONE
        real(wp),allocatable :: varloc1(:,:,:,:)
        real(wp),allocatable :: varglb1(:,:,:,:)
        real(wp),allocatable :: varloc2(:,:)
        real(wp),allocatable :: varglb2(:,:)
        integer(4) :: N, I, JJ, K, KS
        
        !=======================global =============================
        ALLOCATE ( U_F0_tg      (NCL1_tg,  NCL2,  NCL3, NDV+1) )
        ALLOCATE ( U1xztL_F0_tg (NCL2, NDV+1 ) )
        
        !=========================gather U,P==========================
        N = 4
        ALLOCATE ( varloc1(NCL1_tg,0:(N2DO(0)+1),NCL3,N))
        ALLOCATE ( varglb1(NCL1_tg,NCL2,         NCL3,N))
        
        varloc1(1:NCL1_tg,0:(N2DO(0)+1),1:NCL3,1:3) = Q_tg (1:NCL1_tg,0:(N2DO(0)+1),1:NCL3,1:3)
        varloc1(1:NCL1_tg,0:(N2DO(0)+1),1:NCL3,4)   = PR_tg(1:NCL1_tg,0:(N2DO(0)+1),1:NCL3    )
        
        CALL GATHERING_xyzn(varloc1, varglb1, NCL1_tg, N)
        
        IF(MYID==0) THEN
            DO I=1,NCL1_tg
                DO K=1,NCL3
                    DO JJ=1,NCL2
                        U_F0_tg(I, JJ, K,1) = varglb1(I, JJ, K, 1)
                        U_F0_tg(I, JJ, K,3) = varglb1(I, JJ, K, 3)*RCCI1(JJ)
                        U_F0_tg(I, JJ, K,4) = varglb1(I, JJ, K, 4)
                        IF(JJ==1 .and. icase==ipipec) THEN
                            KS = KSYM(K)
                            U_F0_tg(I, JJ, K, 2) = (varglb1(I, 2, K,2) - varglb1(I, 2, KS,2) ) * 0.50_WP*RNDI1(2)
                        ELSE
                            U_F0_tg(I, JJ, K, 2) = varglb1(I, JJ, K, 2)*RNDI1(JJ)
                        END IF
                    END DO
                END DO
            END DO
        END IF
        DEALLOCATE (varloc1)
        DEALLOCATE (varglb1)
        
        !=========================gather averaged==========================
        N = 4 
        ALLOCATE ( varloc2(1:N2DO(0),N))
        ALLOCATE ( varglb2(NCL2,     N))
        
        varloc2(1:N2DO(0),1:4) = U1xztL_tg (1:N2DO(0),1:4)   ! should space and time averaged? or just space averaged?
        
        CALL GATHERING_yn(varloc2, varglb2, N)
        
        IF(MYID==0) THEN
            U1xztL_F0_tg(1:NCL2, 1:4) = varglb2(1:NCL2, 1:4)
        END IF
        
        DEALLOCATE (varloc2)
        DEALLOCATE (varglb2)
        
        
        RETURN
    END SUBROUTINE
    
!**************************************************************************************************
    SUBROUTINE GATHERING_io
        use TEC360_INFO
        IMPLICIT NONE
        real(wp),allocatable :: varloc1(:,:,:,:)
        real(wp),allocatable :: varglb1(:,:,:,:)
        real(wp),allocatable :: varloc2(:,:)
        real(wp),allocatable :: varglb2(:,:)
        real(wp),allocatable :: varloc3(:,:,:)
        real(wp),allocatable :: varglb3(:,:,:)
        integer(4) :: N, I, K, JJ, KS
        
        !======================global===================================
        ALLOCATE ( U_F0_io(NCL1_IO, NCL2,  NCL3, NDV+1) )
        
        IF(thermlflg==1) THEN
            ALLOCATE ( G_F0_io(NCL1_IO, NCL2,  NCL3, NDV) )
            ALLOCATE ( T_F0_io(NCL1_IO, NCL2,  NCL3) )
            ALLOCATE ( D_F0_io(NCL1_IO, NCL2,  NCL3) )
            ALLOCATE ( M_F0_io(NCL1_IO, NCL2,  NCL3) )
            ALLOCATE ( H_F0_io(NCL1_IO, NCL2,  NCL3) )
        END IF
        
        IF(TGFLOWFLG) THEN
            ALLOCATE ( U1ztL_F0_io(NCL1_io,NCL2, NDV+1 ) )
            IF(thermlflg==1) THEN
                ALLOCATE ( G1ztL_F0_io(NCL1_io,NCL2, NDV ) )
                ALLOCATE ( H1ztL_F0_io(NCL1_io,NCL2      ) )
                ALLOCATE ( T1ztL_F0_io(NCL1_io,NCL2      ) )
                ALLOCATE ( D1ztL_F0_io(NCL1_io,NCL2      ) )
                ALLOCATE ( M1ztL_F0_io(NCL1_io,NCL2      ) )
                
            END IF
        ELSE
            ALLOCATE ( U1xztL_F0_io(NCL2, NDV+1 ) )
            IF(thermlflg==1) THEN
                ALLOCATE ( G1xztL_F0_io(NCL2, NDV ) )
                ALLOCATE ( H1xztL_F0_io(NCL2      ) )
                ALLOCATE ( T1xztL_F0_io(NCL2      ) )
                ALLOCATE ( D1xztL_F0_io(NCL2      ) )
                ALLOCATE ( M1xztL_F0_io(NCL2      ) )
            END IF
        END IF
        
        !=========================gather U,P==========================
        IF(thermlflg==1) then
            N = 11
        else
            N = 4
        end if
        ALLOCATE ( varloc1(NCL1_io,0:(N2DO(0)+1),NCL3,N))
        ALLOCATE ( varglb1(NCL1_io,NCL2,         NCL3,N))
        
        varloc1(1:NCL1_io,0:(N2DO(0)+1),1:NCL3,1:3) = Q_io (1:NCL1_io,0:(N2DO(0)+1),1:NCL3,1:3)
        varloc1(1:NCL1_io,0:(N2DO(0)+1),1:NCL3,4)   = PR_io(1:NCL1_io,0:(N2DO(0)+1),1:NCL3    )
        IF(thermlflg==1) then
            varloc1(1:NCL1_io,0:(N2DO(0)+1),1:NCL3,5:7) = G_io       (1:NCL1_io,0:(N2DO(0)+1),1:NCL3,1:3)
            varloc1(1:NCL1_io,0:(N2DO(0)+1),1:NCL3,8)   = ENTHALPY   (1:NCL1_io,0:(N2DO(0)+1),1:NCL3)
            varloc1(1:NCL1_io,0:(N2DO(0)+1),1:NCL3,9)   = TEMPERATURE(1:NCL1_io,0:(N2DO(0)+1),1:NCL3)
            varloc1(1:NCL1_io,0:(N2DO(0)+1),1:NCL3,10)  = DENSITY    (1:NCL1_io,0:(N2DO(0)+1),1:NCL3)
            varloc1(1:NCL1_io,0:(N2DO(0)+1),1:NCL3,11)  = VISCOUSITY (1:NCL1_io,0:(N2DO(0)+1),1:NCL3)
        END IF
        
        CALL GATHERING_xyzn(varloc1, varglb1, NCL1_io, N)
        
        IF(MYID==0) THEN

            DO I=1,NCL1_IO
                DO K=1,NCL3
                    DO JJ=1,NCL2
                        U_F0_IO(I, JJ, K, 1) = varglb1(I, JJ, K, 1)
                        U_F0_IO(I, JJ, K, 3) = varglb1(I, JJ, K, 3)*RCCI1(JJ)
                        U_F0_IO(I, JJ, K, 4) = varglb1(I, JJ, K, 4)
                        IF(JJ==1 .and. icase==ipipec) THEN
                            KS = KSYM(K)
                            U_F0_IO(I, JJ, K, 2) = (varglb1(I, 2, K,2) - varglb1(I, 2, KS,2) ) * 0.50_WP*RNDI1(2)
                        ELSE
                            U_F0_IO(I, JJ, K, 2) = varglb1(I, JJ, K, 2)*RNDI1(JJ)
                        END IF
                    END DO
                END DO
            END DO
        
            IF(thermlflg==1) then
                DO I=1,NCL1_IO
                    DO K=1,NCL3
                        DO JJ=1,NCL2
                            G_F0_IO(I, JJ, K, 1) = varglb1(I, JJ, K, 5)
                            G_F0_IO(I, JJ, K, 3) = varglb1(I, JJ, K, 7)*RCCI1(JJ)
                            IF(JJ==1 .and. icase==ipipec) THEN
                                KS = KSYM(K)
                                G_F0_IO(I, JJ, K, 2) = (varglb1(I, 2, K,6) - varglb1(I, 2, KS,6) ) * 0.50_WP*RNDI1(2)
                            ELSE
                                G_F0_IO(I, JJ, K, 2) = varglb1(I, JJ, K, 6)*RNDI1(JJ)
                            END IF
                        END DO
                    END DO
                END DO
                H_F0_io(1:NCL1_io, 1:NCL2, 1:NCL3) = varglb1(1:NCL1_io, 1:NCL2,  1:NCL3, 8)
                T_F0_io(1:NCL1_io, 1:NCL2, 1:NCL3) = varglb1(1:NCL1_io, 1:NCL2,  1:NCL3, 9)
                D_F0_io(1:NCL1_io, 1:NCL2, 1:NCL3) = varglb1(1:NCL1_io, 1:NCL2,  1:NCL3, 10)
                M_F0_io(1:NCL1_io, 1:NCL2, 1:NCL3) = varglb1(1:NCL1_io, 1:NCL2,  1:NCL3, 11)
            END IF
        END IF
        
        DEALLOCATE (varloc1)
        DEALLOCATE (varglb1)
        
        !=========================gather averaged==========================
        IF(TGFLOWFLG) THEN
            N = 4
            ALLOCATE ( varloc3(NCL1_io,1:N2DO(0),N))
            ALLOCATE ( varglb3(NCL1_io,NCL2,     N))
            
            varloc3(NCL1_io,1:N2DO(0),1:4) = U1ztL_io (NCL1_io,1:N2DO(0),1:4)
            
            CALL GATHERING_xyn(varloc3, varglb3, NCL1_io,N)
            
            IF(MYID==0) THEN
                U1ztL_F0_io(NCL1_io,1:NCL2, 1:4) = varglb3(NCL1_io,1:NCL2, 1:4)
            END IF
            DEALLOCATE (varloc3)
            DEALLOCATE (varglb3)
            
            IF(thermlflg==1) then
                N = 6
                
                ALLOCATE ( varloc3(NCL1_io,1:N2DO(0),N))
                ALLOCATE ( varglb3(NCL1_io,NCL2,     N))
                
                varloc3(NCL1_io,1:N2DO(0),1:3) = G1ztL_io (NCL1_io,1:N2DO(0),1:3)
                varloc3(NCL1_io,1:N2DO(0),4)   = H1ztL_io (NCL1_io,1:N2DO(0)    )
                varloc3(NCL1_io,1:N2DO(0),5)   = T1ztL_io (NCL1_io,1:N2DO(0)    )
                varloc3(NCL1_io,1:N2DO(0),6)   = D1ztL_io (NCL1_io,1:N2DO(0)    )
                !varloc3(NCL1_io,1:N2DO(0),7)   = M1ztL_io (NCL1_io,1:N2DO(0)    )
                
                CALL GATHERING_xyn(varloc3, varglb3, NCL1_io,N)
                
                IF(MYID==0) THEN
                    G1ztL_F0_io(NCL1_io,1:NCL2, 1:3) = varglb3(NCL1_io,1:NCL2, 1:3)
                    H1ztL_F0_io(NCL1_io,1:NCL2)      = varglb3(NCL1_io,1:NCL2, 4  )
                    T1ztL_F0_io(NCL1_io,1:NCL2)      = varglb3(NCL1_io,1:NCL2, 5  )
                    D1ztL_F0_io(NCL1_io,1:NCL2)      = varglb3(NCL1_io,1:NCL2, 6  )
                    !M1ztL_F0_io(NCL1_io,1:NCL2)      = varglb3(NCL1_io,1:NCL2, 7  )
                END IF
            
                DEALLOCATE (varloc3)
                DEALLOCATE (varglb3)
            END IF
            
        ELSE
            N = 4
            
            ALLOCATE ( varloc2(1:N2DO(0),N))
            ALLOCATE ( varglb2(NCL2,     N))
            
            varloc2(1:N2DO(0),1:4) = U1xztL_io (1:N2DO(0),1:4)
            
            CALL GATHERING_yn(varloc2, varglb2, N)
            
            IF(MYID==0) THEN
                U1xztL_F0_io(1:NCL2, 1:4) = varglb2(1:NCL2, 1:4)
            END IF
            
            DEALLOCATE (varloc2)
            DEALLOCATE (varglb2)
            
            IF(thermlflg==1) then
                N = 7
                
                ALLOCATE ( varloc2(1:N2DO(0),N))
                ALLOCATE ( varglb2(NCL2,     N))
                
                varloc2(1:N2DO(0),1:3) = G1xztL_io (1:N2DO(0),1:3)
                varloc2(1:N2DO(0),4)   = H1xztL_io (1:N2DO(0)    )
                varloc2(1:N2DO(0),5)   = T1xztL_io (1:N2DO(0)    )
                varloc2(1:N2DO(0),6)   = D1xztL_io (1:N2DO(0)    )
                varloc2(1:N2DO(0),7)   = M1xztL_io (1:N2DO(0)    )
                
                
                CALL GATHERING_yn(varloc2, varglb2, N)
                
                IF(MYID==0) THEN
                    G1xztL_F0_io(1:NCL2, 1:3) = varglb2(1:NCL2, 1:3)
                    H1xztL_F0_io(1:NCL2)      = varglb2(1:NCL2, 4  )
                    T1xztL_F0_io(1:NCL2)      = varglb2(1:NCL2, 5  )
                    D1xztL_F0_io(1:NCL2)      = varglb2(1:NCL2, 6  )
                    M1xztL_F0_io(1:NCL2)      = varglb2(1:NCL2, 7 )
                END IF
            
                DEALLOCATE (varloc2)
                DEALLOCATE (varglb2)
            END IF
            
            
        END IF
        
        RETURN
    END SUBROUTINE
    
!**************************************************************************************************
    SUBROUTINE GATHERING_xyzn(varloc, varglb, NSX, NSN)
        use TEC360_INFO
        implicit none
        INTEGER(4) :: NSX, NSN
        REAL(WP),INTENT(IN)  :: varloc(nsx,0:N2DO(0)+1,NCL3, NSN)
        REAL(WP),INTENT(OUT) :: varglb(nsx,1:NCL2,     NCL3, NSN)
        REAL(WP)              :: AUX   (nsx,0:N2DO(0)+1,NCL3, NSN,1:NPTOT)
        INTEGER(4)            :: INN, N, KK, J, N2DOID, I, JJ, K
        
        INN = nsx*(N2DO(0)+2)*NCL3*NSN
        CALL MPI_GATHER( varloc, INN, MPI_DOUBLE_PRECISION, AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)                
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !==============re-arrange gathered data======================== 
        IF (MYID.EQ.0) THEN
            DO N =1, nsn
                DO KK=0,NPSLV
                    N2DOID=JDEWT(KK)-JDSWT(KK)+1
                    DO  J=1,N2DOID
                        JJ=JDSWT(KK)-1+J
                        DO I=1,NSX
                            DO K=1,NCL3
                                varglb(I,JJ,K, N) = AUX(I,J,K,N, KK+1)
                            ENDDO
                        ENDDO
                    ENDDO
                END DO
            END DO
   
        END IF

        return
    END SUBROUTINE  
    
    
    SUBROUTINE GATHERING_yn(varloc, varglb,NSN)
        use tec360_info
        implicit none
        INTEGER(4) :: NSN
        REAL(WP),INTENT(IN)  :: varloc(N2DO(MYID),   NSN)
        REAL(WP),INTENT(OUT) :: varglb(NCL2,         NSN)
        REAL(WP)              :: AUX  (N2DO(MYID),   NSN,1:NPTOT)
        INTEGER(4)            :: INN, N, IP, J, N2DOID, I, JJ
        
        INN = N2DO(0)*NSN
        CALL MPI_GATHER( varloc, INN, MPI_DOUBLE_PRECISION, AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)                
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !==============re-arrange gathered data======================== 
        IF (MYID.EQ.0) THEN
            
            DO IP=0,NPSLV
                N2DOID=JDEWT(IP)-JDSWT(IP)+1
                DO  J=1,N2DOID
                    JJ=JDSWT(IP)-1+J
                    
                    DO N =1, nsn
                        varglb(JJ,N) = AUX(J,N, IP+1)
                    ENDDO
                    
                END DO
            END DO
   
        END IF

        return
    END SUBROUTINE
    
    SUBROUTINE GATHERING_xyn(varloc, varglb, NSX, NSN)
        use tec360_info
        implicit none
        INTEGER(4) :: NSX, NSN
        REAL(WP),INTENT(IN)  :: varloc(nsx,1:N2DO(0), NSN)
        REAL(WP),INTENT(OUT) :: varglb(nsx,1:NND2,    NSN)
        REAL(WP)              :: AUX   (nsx,1:N2DO(0), NSN,1:NPTOT)
        INTEGER(4)            :: INN, N, IP, J, N2DOID, I, JJ
        
        INN = nsx*(N2DO(0))*NSN
        CALL MPI_GATHER( varloc, INN, MPI_DOUBLE_PRECISION, AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)                
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !==============re-arrange gathered data======================== 
        IF (MYID.EQ.0) THEN
            
            DO IP=0,NPSLV
                N2DOID=JDEWT(IP)-JDSWT(IP)+1
                DO  J=1,N2DOID
                    JJ=JDSWT(IP)-1+J
                    DO I=1,NSX
                        DO N =1, nsn
                            varglb(I,JJ,N) = AUX(I,J,N,IP+1)
                        ENDDO
                    ENDDO
                END DO
            END DO
   
        END IF

        return
    END SUBROUTINE
    

!***************************************************************************************************
    SUBROUTINE TEC360_EXTRPLAT2ND_MASTER_TG
        use TEC360_INFO
        IMPLICIT NONE
      
        INTEGER(4)  :: I, J, K, IM, JM, KM, KS
      
       
        IF (MYID.NE.0) RETURN
        
        ALLOCATE ( U1xztL_INTP_tg(NND2, NDV+1 ) )
        
        DO J=2,NCL2
            JM = J-1
            U1xztL_INTP_tg(J,1:4) = (U1xztL_F0_tg(JM,1:4)+ U1xztL_F0_tg(J,1:4))*0.5_WP
        END DO
        U1xztL_INTP_tg(1,   1:3) = 0.0_WP
        U1xztL_INTP_tg(NND2,1:3) = 0.0_WP
        U1xztL_INTP_tg(1,   4) =  U1xztL_F0_tg(1,   4)
        U1xztL_INTP_tg(NND2,4) =  U1xztL_F0_tg(NCL2,4)

        
        ALLOCATE ( U_INTP_tg(NND1_tg,NND2,NND3, NDV+1) )

        !>=============INTERPOLATION ALL VALUES TO POINTS FOR X, Z PERIODIC B.C.==========
        
        DO I=1,NCL1_tg
            IM = IMV_tg(I)
            DO K=1,NCL3
                KM = KMV(K)
                KS = KSYM(K)
                DO J=2,NCL2
                    JM = JGMV(J)
                    U_INTP_tg(I,J,K,1) = 0.250_WP * ( U_F0_tg(I,JM,KM,1) + &
                                                      U_F0_tg(I,JM,K, 1)  + &
                                                      U_F0_tg(I,J, KM,1) + &
                                                      U_F0_tg(I,J,K,  1)    )
                    U_INTP_tg(I,J,K,3) = 0.250_WP * ( U_F0_tg(IM,JM,K,3) + &
                                                      U_F0_tg(IM,J,K, 3)  + &
                                                      U_F0_tg(I,JM,K, 3)  + &
                                                      U_F0_tg(I,J,K,  3) )
                    U_INTP_tg(I,J,K,2) = 0.250_WP * ( U_F0_tg(IM,J,KM,2) + &
                                                      U_F0_tg(IM,J,K ,2)  + &
                                                      U_F0_tg(I,J,KM ,2)  + &
                                                      U_F0_tg(I,J,K  ,2) )
                    U_INTP_tg(I,J,K,4) = 0.1250_WP * (U_F0_tg(IM,JM,KM,4) + U_F0_tg(IM,JM,K,4) + &
                                                      U_F0_tg(IM,J ,KM,4) + U_F0_tg(IM,J ,K,4) + &
                                                      U_F0_tg(I ,JM,KM,4) + U_F0_tg(I ,JM,K,4) + &
                                                      U_F0_tg(I ,J ,KM,4) + U_F0_tg(I ,J ,K,4) )
                END DO
                IF(icase==IPIPEC) THEN
                    U_INTP_tg(I,1,K,1:4) = 0.5_wp*( U_F0_tg(I,1,K,1:4) + U_F0_tg(I,1,KS,1:4))
                ELSE
                    U_INTP_tg(I,1,K,1:3) = 0.0_WP
                    U_INTP_tg(I,1,K, 4)  = U_F0_tg(I,1,K,4)
                END IF
                
            END DO
        END DO
        U_INTP_tg(NND1_tg,:,:,:) = U_INTP_tg(1,:,:,:)
        
        U_INTP_tg(:,:,NND3,:)    = U_INTP_tg(:,:,1,:)
        
        U_INTP_tg(:,NND2,:,1:3)    = 0.0_WP
        U_INTP_tg(:,NND2,:, 4)     = U_INTP_tg(:,NCL2,:,4)
        

        RETURN
    END SUBROUTINE

!********************************************************************************************************      
    SUBROUTINE TEC360_EXTRPLAT2ND_MASTER_IO
        use TEC360_INFO
        IMPLICIT NONE
      
        INTEGER(4)  :: I, J, K, IM, JM, KM, KS, IE
            
        IF(.not.IOFLOWflg) RETURN
        IF (MYID.NE.0) RETURN 
        
        ALLOCATE ( U_INTP_io(NND1_io,NND2,NND3, NDV+1) )
        
        IF(thermlflg==1) THEN
            ALLOCATE ( G_INTP_io(NND1_io,NND2,NND3, NDV) )
            ALLOCATE ( H_INTP_io(NND1_io,NND2,NND3) )
            ALLOCATE ( T_INTP_io(NND1_io,NND2,NND3) )
            ALLOCATE ( D_INTP_io(NND1_io,NND2,NND3) )
            ALLOCATE ( M_INTP_io(NND1_io,NND2,NND3) )
        END IF
        
        IF(TGFLOWFLG) THEN
            ALLOCATE ( U1ztL_INTP_io(NND1_io,NND2, NDV+1 ) )
            IF(thermlflg==1) THEN
                ALLOCATE ( G1ztL_INTP_io(NND1_io,NND2, NDV ) )
                ALLOCATE ( H1ztL_INTP_io(NND1_io,NND2) )
                ALLOCATE ( T1ztL_INTP_io(NND1_io,NND2) )
                ALLOCATE ( D1ztL_INTP_io(NND1_io,NND2) )
                ALLOCATE ( M1ztL_INTP_io(NND1_io,NND2) )
            END IF
        ELSE
            ALLOCATE ( U1xztL_INTP_io(NND2, NDV+1 ) )
            IF(thermlflg==1) THEN
                ALLOCATE ( G1xztL_INTP_io(NND2, NDV ) )
                ALLOCATE ( H1xztL_INTP_io(NND2) )
                ALLOCATE ( T1xztL_INTP_io(NND2) )
                ALLOCATE ( D1xztL_INTP_io(NND2) )
                ALLOCATE ( M1xztL_INTP_io(NND2) )
            END IF
        END IF
            
        !>============INTERPOLATION ALL VALUES TO POINTS FOR Z PERIODIC B.C. X INLET/OUTLET============================
        IE = NCL1_io
        IF(TGFLOWflg) IE = NND1_io
        
        
        IF(TGFLOWFLG) THEN
        
            DO I=1,IE
                IF(I==NND1_io) THEN 
                    IM = NCL1_io 
                ELSE
                    IM = IMV_io(I)
                END IF
                DO J=2,NCL2
                    JM = JGMV(J)
                    U1ztL_INTP_io(I, J,1:4) = ( U1ztL_F0_io(IM, JM,1:4)+ U1ztL_F0_io(IM, J,1:4) + &
                                                U1ztL_F0_io(I , JM,1:4)+ U1ztL_F0_io(I , J,1:4) )*0.25_WP
                END DO
                U1ztL_INTP_io(I, 1,   1:3) = 0.0_WP
                U1ztL_INTP_io(I, NND2,1:3) = 0.0_WP
                U1ztL_INTP_io(I, 1,   4)   =  U1ztL_F0_io(I, 1,   4)
                U1ztL_INTP_io(I, NND2,4)   =  U1ztL_F0_io(I, NCL2,4)
            END DO
            
            IF(thermlflg==1) THEN
                DO I=1,IE
                    IF(I==NND1_io) THEN 
                        IM = NCL1_io 
                    ELSE
                        IM = IMV_io(I)
                    END IF
                    DO J=2,NCL2
                        JM = JGMV(J)
                        G1ztL_INTP_io(I, J,1:3) = ( G1ztL_F0_io(IM, JM,1:3)+ G1ztL_F0_io(IM, J,1:3) + &
                                                    G1ztL_F0_io(I , JM,1:3)+ G1ztL_F0_io(I , J,1:3) )*0.25_WP
                        H1ztL_INTP_io(I, J    ) = ( H1ztL_F0_io(IM, JM    )+ H1ztL_F0_io(IM, J    ) + &
                                                    H1ztL_F0_io(I , JM    )+ H1ztL_F0_io(I , J    ) )*0.25_WP
                        T1ztL_INTP_io(I, J    ) = ( T1ztL_F0_io(IM, JM    )+ T1ztL_F0_io(IM, J    ) + &
                                                    T1ztL_F0_io(I , JM    )+ T1ztL_F0_io(I , J    ) )*0.25_WP
                        D1ztL_INTP_io(I, J    ) = ( D1ztL_F0_io(IM, JM    )+ D1ztL_F0_io(IM, J    ) + &
                                                    D1ztL_F0_io(I , JM    )+ D1ztL_F0_io(I , J    ) )*0.25_WP
                        M1ztL_INTP_io(I, J    ) = ( M1ztL_F0_io(IM, JM    )+ M1ztL_F0_io(IM, J    ) + &
                                                    M1ztL_F0_io(I , JM    )+ M1ztL_F0_io(I , J    ) )*0.25_WP
                        
                    END DO
                    G1ztL_INTP_io(I, 1,   1:3) = 0.0_WP
                    G1ztL_INTP_io(I, NND2,1:3) = 0.0_WP
                    H1ztL_INTP_io(I, 1)      =  H1ztL_F0_io(I, 1   )
                    H1ztL_INTP_io(I, NND2)   =  H1ztL_F0_io(I, NCL2)
                    T1ztL_INTP_io(I, 1)      =  T1ztL_F0_io(I, 1   )
                    T1ztL_INTP_io(I, NND2)   =  T1ztL_F0_io(I, NCL2)
                    D1ztL_INTP_io(I, 1)      =  D1ztL_F0_io(I, 1   )
                    D1ztL_INTP_io(I, NND2)   =  D1ztL_F0_io(I, NCL2)
                    M1ztL_INTP_io(I, 1)      =  M1ztL_F0_io(I, 1   )
                    M1ztL_INTP_io(I, NND2)   =  M1ztL_F0_io(I, NCL2)
                    
                END DO
            END IF
            
        ELSE
            DO J=2,NCL2
                JM = J-1
                U1xztL_INTP_io(J,1:4) = (U1xztL_F0_io(JM,1:4)+ U1xztL_F0_io(J,1:4))*0.5_WP
            END DO
            U1xztL_INTP_io(1,   1:3) = 0.0_WP
            U1xztL_INTP_io(NND2,1:3) = 0.0_WP
            U1xztL_INTP_io(1,   4) =  U1xztL_F0_io(1,   4)
            U1xztL_INTP_io(NND2,4) =  U1xztL_F0_io(NCL2,4)
            
            IF(thermlflg==1) THEN
                DO J=2,NCL2
                    JM = J-1
                    G1xztL_INTP_io(J,1:3) = (G1xztL_F0_io(JM,1:3)+ G1xztL_F0_io(J,1:3))*0.5_WP
                    H1xztL_INTP_io(J    ) = (H1xztL_F0_io(JM    )+ H1xztL_F0_io(J    ))*0.5_WP
                    T1xztL_INTP_io(J    ) = (T1xztL_F0_io(JM    )+ T1xztL_F0_io(J    ))*0.5_WP
                    D1xztL_INTP_io(J    ) = (D1xztL_F0_io(JM    )+ D1xztL_F0_io(J    ))*0.5_WP
                    M1xztL_INTP_io(J    ) = (M1xztL_F0_io(JM    )+ M1xztL_F0_io(J    ))*0.5_WP
                END DO
                G1xztL_INTP_io(1,   1:3) = 0.0_WP
                G1xztL_INTP_io(NND2,1:3) = 0.0_WP
                H1xztL_INTP_io(1)    =  HWAL_RA(ibotwall)
                H1xztL_INTP_io(NND2) =  HWAL_RA(itopwall)
                T1xztL_INTP_io(1)    =  TWAL(ibotwall)
                T1xztL_INTP_io(NND2) =  TWAL(itopwall)
                D1xztL_INTP_io(1)    =  DWAL(ibotwall)
                D1xztL_INTP_io(NND2) =  DWAL(itopwall)
                M1xztL_INTP_io(1)    =  MWAL(ibotwall)
                M1xztL_INTP_io(NND2) =  MWAL(itopwall)
            END IF
            
        END IF

        
        DO I=1,IE
            IF(I==NND1_io) THEN 
                IM = NCL1_io 
            ELSE
                IM = IMV_io(I)
            END IF
            DO K=1,NCL3
                KM = KMV(K)
                KS = KSYM(K)
                DO J=2,NCL2
                    JM = JGMV(J)
                    U_INTP_io(I,J,K,1) = 0.250_WP * ( U_F0_io(I,JM,KM,1) + &
                                                    U_F0_io(I,JM,K ,1)  + &
                                                    U_F0_io(I,J,KM ,1)  + &
                                                    U_F0_io(I,J,K  ,1) )
                    U_INTP_io(I,J,K,2) = 0.250_WP * ( U_F0_io(IM,J,KM,2) + &
                                                    U_F0_io(IM,J,K ,2)  + &
                                                    U_F0_io(I,J,KM ,2)  + &
                                                    U_F0_io(I,J,K  ,2) )
                    U_INTP_io(I,J,K,3) = 0.250_WP * ( U_F0_io(IM,JM,K,3) + &
                                                    U_F0_io(IM,J,K ,3)  + &
                                                    U_F0_io(I,JM,K ,3)  + &
                                                    U_F0_io(I,J,K  ,3) )
                    U_INTP_io(I,J,K,4) = 0.1250_WP * &
                                                  ( U_F0_io(IM,JM,KM,4) + U_F0_io(IM,JM,K,4) + &
                                                    U_F0_io(IM,J ,KM,4) + U_F0_io(IM,J ,K,4) + &
                                                    U_F0_io(I ,JM,KM,4) + U_F0_io(I ,JM,K,4) + &
                                                    U_F0_io(I ,J ,KM,4) + U_F0_io(I ,J ,K,4) )
                END DO
                IF(icase==IPIPEC) THEN
                    U_INTP_IO(I,1,K,:) = 0.5_wp*( U_F0_IO(I,1,K,:) + U_F0_IO(I,1,KS,:))
                ELSE
                    U_INTP_IO(I,1,K,1:3) = 0.0_WP
                    U_INTP_IO(I,1,K,4)   = U_F0_io(I,1,K,4)
                END IF
            END DO
        END DO
        
        IF(TGFLOWFLG) THEN
            U_INTP_IO(NND1_io,:,:,:) = U_F0_io(NCL1_io,:,:,:)
        ELSE
            U_INTP_IO(NND1_io,:,:,:) = U_INTP_IO(1,:,:,:)
        END IF
        
        U_INTP_io(:,:,NND3,:) = U_INTP_io(:,:,1,:)
        
        U_INTP_io(:,NND2,:,1:3) = 0.0_WP
        U_INTP_io(:,NND2,:,4)   = U_INTP_io(:,NCL2,:,4)
         
         
        IF(thermlflg==1) THEN
            DO I=1,NCL1_io
                IF(I==NND1_io) THEN 
                    IM = NCL1_io 
                ELSE
                    IM = IMV_io(I)
                END IF
                
                DO K=1,NCL3
                    KM = KMV(K)
                    KS = KSYM(K)
                    DO J=2,NCL2
                        JM = JGMV(J)
                        G_INTP_io(I,J,K,1) = 0.250_WP * ( G_F0_io(I,JM,KM,1) + &
                                                    G_F0_io(I,JM,K ,1)  + &
                                                    G_F0_io(I,J,KM ,1)  + &
                                                    G_F0_io(I,J,K  ,1) )
                        G_INTP_io(I,J,K,2) = 0.250_WP * ( G_F0_io(IM,J,KM,2) + &
                                                    G_F0_io(IM,J,K ,2)  + &
                                                    G_F0_io(I,J,KM ,2)  + &
                                                    G_F0_io(I,J,K  ,2) )
                        G_INTP_io(I,J,K,3) = 0.250_WP * ( G_F0_io(IM,JM,K,3) + &
                                                    G_F0_io(IM,J,K ,3)  + &
                                                    G_F0_io(I,JM,K ,3)  + &
                                                    G_F0_io(I,J,K  ,3) )
                                                    
                        H_INTP_io(I,J,K) = 0.1250_WP * &
                                         ( H_F0_io(IM,JM,KM) + H_F0_io(IM,JM,K) + &
                                           H_F0_io(IM,J ,KM) + H_F0_io(IM,J ,K) + &
                                           H_F0_io(I ,JM,KM) + H_F0_io(I ,JM,K) + &
                                           H_F0_io(I ,J ,KM) + H_F0_io(I ,J ,K) )
                                                                       
                        T_INTP_io(I,J,K) = 0.1250_WP * &
                                         ( T_F0_io(IM,JM,KM) + T_F0_io(IM,JM,K) + &
                                           T_F0_io(IM,J ,KM) + T_F0_io(IM,J ,K) + &
                                           T_F0_io(I ,JM,KM) + T_F0_io(I ,JM,K) + &
                                           T_F0_io(I ,J ,KM) + T_F0_io(I ,J ,K) )
                        D_INTP_io(I,J,K) = 0.1250_WP * &
                                         ( D_F0_io(IM,JM,KM) + D_F0_io(IM,JM,K) + &
                                           D_F0_io(IM,J ,KM) + D_F0_io(IM,J ,K) + &
                                           D_F0_io(I ,JM,KM) + D_F0_io(I ,JM,K) + &
                                           D_F0_io(I ,J ,KM) + D_F0_io(I ,J ,K) )
                        M_INTP_io(I,J,K) = 0.1250_WP * &
                                         ( M_F0_io(IM,JM,KM) + M_F0_io(IM,JM,K) + &
                                           M_F0_io(IM,J ,KM) + M_F0_io(IM,J ,K) + &
                                           M_F0_io(I ,JM,KM) + M_F0_io(I ,JM,K) + &
                                           M_F0_io(I ,J ,KM) + M_F0_io(I ,J ,K) )
                        
                    END DO
                    IF(icase==IPIPEC) THEN
                        G_INTP_IO(I,1,K,:) = 0.5_wp*( G_F0_IO(I,1,K,:) + G_F0_IO(I,1,KS,:))
                    ELSE
                        G_INTP_IO(I,1,K,:) = 0.0_WP
                    END IF
                END DO
            END DO

            IF(TGFLOWFLG) THEN
                G_INTP_IO(NND1_io,:,:,:) = G_F0_io(NCL1_io,:,:,:)
                H_INTP_io(NND1_io,1:NCL2,1:NCL3) = H_F0_io(NCL1_io,1:NCL2,1:NCL3)
                T_INTP_io(NND1_io,1:NCL2,1:NCL3) = T_F0_io(NCL1_io,1:NCL2,1:NCL3)
                D_INTP_io(NND1_io,1:NCL2,1:NCL3) = D_F0_io(NCL1_io,1:NCL2,1:NCL3)
                M_INTP_io(NND1_io,1:NCL2,1:NCL3) = M_F0_io(NCL1_io,1:NCL2,1:NCL3)
            ELSE
                G_INTP_IO(NND1_io,:,:,:) = G_INTP_IO(1,:,:,:)
                H_INTP_io(NND1_io,1:NCL2,1:NCL3) = H_INTP_io(1,1:NCL2,1:NCL3)
                T_INTP_io(NND1_io,1:NCL2,1:NCL3) = T_INTP_io(1,1:NCL2,1:NCL3)
                D_INTP_io(NND1_io,1:NCL2,1:NCL3) = D_INTP_io(1,1:NCL2,1:NCL3)
                M_INTP_io(NND1_io,1:NCL2,1:NCL3) = M_INTP_io(1,1:NCL2,1:NCL3)
            END IF
            G_INTP_io(:,:,NND3,:) = G_INTP_io(:,:,1,:)
            G_INTP_io(:,NND2,:,:) = 0.0_WP
            H_INTP_io(:,:,NND3) = H_INTP_io(:,:,1)
            T_INTP_io(:,:,NND3) = T_INTP_io(:,:,1)
            D_INTP_io(:,:,NND3) = D_INTP_io(:,:,1)
            M_INTP_io(:,:,NND3) = M_INTP_io(:,:,1)
            


            IF(BCWALLHEAT(ibotwall)==isoFluxWall .and. icase==IPIPEC) THEN
                H_INTP_io(:,1,:)    = 2.0_wp*H_INTP_io(:,2,:)-H_INTP_io(:,3,:)
                T_INTP_io(:,1,:)    = 2.0_wp*T_INTP_io(:,2,:)-T_INTP_io(:,3,:)
                D_INTP_io(:,1,:)    = 2.0_wp*D_INTP_io(:,2,:)-D_INTP_io(:,3,:)
                M_INTP_io(:,1,:)    = 2.0_wp*M_INTP_io(:,2,:)-M_INTP_io(:,3,:)
            ELSE IF(BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                H_INTP_io(:,1,:)    = H_WAL_GV(1,ibotwall) 
                T_INTP_io(:,1,:)    = T_WAL_GV(1,ibotwall) 
                D_INTP_io(:,1,:)    = D_WAL_GV(1,ibotwall)  
                M_INTP_io(:,1,:)    = M_WAL_GV(1,ibotwall)  
            ELSE
            END IF
            
            
            IF(BCWALLHEAT(itopwall)==isoFluxWall) THEN
                H_INTP_io(:,NND2,:) = 2.0_wp*H_INTP_io(:,NCL2,:)-H_INTP_io(:,NCL2-1,:)
                T_INTP_io(:,NND2,:) = 2.0_wp*T_INTP_io(:,NCL2,:)-T_INTP_io(:,NCL2-1,:)
                D_INTP_io(:,NND2,:) = 2.0_wp*D_INTP_io(:,NCL2,:)-D_INTP_io(:,NCL2-1,:)
                M_INTP_io(:,NND2,:) = 2.0_wp*M_INTP_io(:,NCL2,:)-M_INTP_io(:,NCL2-1,:)
            END IF
            IF(BCWALLHEAT(itopwall)==isoThermalWall) THEN
                H_INTP_io(:,NND2,:) = H_WAL_GV(1,itopwall) 
                T_INTP_io(:,NND2,:) = T_WAL_GV(1,itopwall) 
                D_INTP_io(:,NND2,:) = D_WAL_GV(1,itopwall) 
                M_INTP_io(:,NND2,:) = M_WAL_GV(1,itopwall) 
            END IF
            
        END IF
         
      
        RETURN
    END SUBROUTINE
!**********************************************************************************************************************************
    SUBROUTINE TEC360_ALL_NODES
        use TEC360_INFO
        IMPLICIT NONE
      
        INTEGER(4)      :: I, J, K, N1
        INTEGER(4)      :: TECFLG = 202
        character(128) :: FLNM
        logical        :: file_exists
        real(wp)       :: uprime(3)
           

        FLNM = TRIM(filepath5)//'RESULT.UVWP.ALLPOINTS.tec'
        
        INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
        
        IF(IOFLOWflg .AND. TGFLOWFLG) THEN 
            N1 = NND1_tg+NND1_io
        ELSE
            IF(TGFLOWflg) N1 = NND1_tg
            IF(IOFLOWFLG) N1 = NND1_io
        END IF
           
        IF(file_exists) then
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)), position='append')
            
            WRITE(TECFLG,'(A,1I13.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') &
                      'ZONE T=" ',ITERG,phyTIME,' ", I=', &
                       N1, ', J=',NND2,', K=',NND3,', F=POINT'
            write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'
            DO K=1,NND3
               DO J=1,NND2
                    IF(TGFLOWflg) THEN
                        DO I =1,NND1_tg   
                            WRITE(TECFLG,'(30ES15.7)') U_INTP_tg(I,J,K,1),U_INTP_tg(I,J,K,2),U_INTP_tg(I,J,K,3),&
                                U_INTP_tg(I,J,K,4), &
                                1.0_WP, 1.0_WP, 0.0_WP, &
                                Vor_tg(I,J,K,1:4), Qcr_tg(I,J,K), delta_tg(I,J,K), lambda2_tg(I,J,K), swirlstrength_tg(I,J,K,3), &
                                uprime_tg(I,J,K,1:4),0.0_WP,0.0_WP,0.0_WP,0.0_WP, &
                                0.0_WP,0.0_WP, udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END DO
                    END IF
                    IF(IOFLOWflg) THEN
                        DO I =1,NND1_io  
                            IF(thermlflg==1) THEN   
                                WRITE(TECFLG,'(30ES15.7)') &
                                U_INTP_io(I,J,K,1:4), &
                                T_INTP_io(I,J,K), D_INTP_io(I,J,K), H_INTP_io(I,J,K), &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3), &
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                                Dprime_io(I,J,K),Mprime_io(I,J,K),udprime_io(I,J,K,1:3)
                            ELSE
                                WRITE(TECFLG,'(30ES15.7)') &
                                U_INTP_io(I,J,K,1:4), &
                                1.0_WP, 1.0_WP, 0.0_WP,  &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3), &
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),0.0_WP, &
                                0.0_wp, 0.0_wp, udprime_io(I,J,K,1),0.0_WP,0.0_WP
                            END IF
                        END DO
                    END IF
                END DO
            END DO
            
        else
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)))
            WRITE(TECFLG,'(A)') 'TITLE = " '//TRIM(ADJUSTL(teczonename))//' " '
            WRITE(TECFLG,'(A)',advance="no") 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
            WRITE(TECFLG,'(A)',advance="no") '"Vorx","Vory","Vorz","VorM","Qcr","delta","lambda2",'
            WRITE(TECFLG,'(A)',advance="no") '"swirlstrengthR1","swirlstrengthR2","swirlstrengthI2",'
            WRITE(TECFLG,'(A)',advance="no") '"uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime",'
            WRITE(TECFLG,'(A)')              '"Tprime", "Dprime","Mprime","uvprime", "udprime", "vdprime"'
        
            WRITE(TECFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') &
                  'ZONE T=" ',ITERG,phyTIME,' ", I=', &
                   N1, ', J=',NND2,', K=',NND3,', F=POINT'
            DO K=1,NND3
               DO J=1,NND2
                    IF(TGFLOWflg) THEN
                        DO I =1,NND1_tg    
                            WRITE(TECFLG,'(33ES15.7)') XND_tg(I),YND(J),ZND(K), &
                            U_INTP_tg(I,J,K,1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_tg(I,J,K,1:4), Qcr_tg(I,J,K), delta_tg(I,J,K), lambda2_tg(I,J,K), swirlstrength_tg(I,J,K,1:3), &
                            uprime_tg(I,J,K,1:4),0.0_WP,0.0_WP,0.0_WP,0.0_WP, &
                            0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END DO
                    END IF
                    IF(IOFLOWflg) THEN
                        DO I =1,NND1_io  
                            IF(thermlflg==1) THEN 
                                WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                                U_INTP_io(I,J,K,1:4), &
                                T_INTP_io(I,J,K), D_INTP_io(I,J,K), H_INTP_io(I,J,K), &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3), &
                                uprime_io(I,J,K,1:4),gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                                Dprime_io(I,J,K),Mprime_io(I,J,K),udprime_io(I,J,K,1:3)
                            ELSE
                                WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                                U_INTP_io(I,J,K,1:4), &
                                1.0_WP, 1.0_WP, 0.0_WP, &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3), &
                                uprime_io(I,J,K,1:4),gprime_io(I,J,K,1:3),0.0_WP, &
                                0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                            END IF
                        END DO
                    END IF
               END DO
            END DO
            
        end if
        
        
        IF(TGFLOWFLG .and. IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(NND2)
        END IF
        
 
        CLOSE(TECFLG)
     
     
    END SUBROUTINE
     
!**********************************************************************************************************************************
    SUBROUTINE TEC360_XSLICE(N)
        use TEC360_INFO
        IMPLICIT NONE
      
        INTEGER(4) :: I, J, K, N
        CHARACTER(4) :: PNTIM
        INTEGER(4)  :: TECFLG = 202
        character(128) :: FLNM
        logical        :: file_exists
        real(wp)       :: uprime(3)
        
        
        
        IF(.NOT.IOFLOWflg) RETURN
        TECFLG = TECFLG + IID(N)
        
        WRITE(PNTIM,'(I4.4)') IID(N)
        
        FLNM = TRIM(filepath5)//'RESULT.UVWP.XSLICE.'//PNTIM//'.tec'
        
        INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
        
        IF(file_exists) THEN
        
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)), position='append')
            WRITE(TECFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                          ' ", I=', 1, ', J=',NND2,', K=',NND3,', F=POINT'
            write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'
    
            DO K=1,NND3
                DO J=1,NND2
                    DO I =IID(N),IID(N)      
                        IF(thermlflg==1) THEN    
                            WRITE(TECFLG,'(30ES15.7)') &
                                U_INTP_io(I,J,K,1:4), &
                                T_INTP_io(I,J,K),D_INTP_io(I,J,K), H_INTP_io(I,J,K), &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K),swirlstrength_io(I,J,K,1:3),&
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                                Dprime_io(I,J,K),Mprime_io(I,J,K), udprime_io(I,J,K,1:3)
                        ELSE
                            WRITE(TECFLG,'(30ES15.7)') &
                                U_INTP_io(I,J,K,1:4), &
                                1.0_WP, 1.0_WP, 0.0_WP, &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),0.0_WP, &
                                0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END IF
                    END DO
                END DO
            END DO
            
            
        ELSE 
        
        
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)) )
            WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW X-SLICE"'
            WRITE(TECFLG,'(A)',advance="no") 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
            WRITE(TECFLG,'(A)',advance="no") '"Vorx","Vory","Vorz","VorM","Qcr","delta","lambda2",'
            WRITE(TECFLG,'(A)',advance="no") '"swirlstrengthR1","swirlstrengthR2","swirlstrengthI2",'
            WRITE(TECFLG,'(A)',advance="no") '"uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime",'
            WRITE(TECFLG,'(A)')              '"Tprime", "Dprime","Mprime","uvprime", "udprime", "vdprime"'
            
        
            WRITE(TECFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG, phyTIME, &
                      ' ", I=', 1, ', J=',NND2,', K=',NND3,', F=POINT'
                      
            DO K=1,NND3
                DO J=1,NND2
                    DO I =IID(N),IID(N)   
                        IF(thermlflg==1) THEN            
                            WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                            U_INTP_io(I,J,K,1:4), &
                            T_INTP_io(I,J,K),D_INTP_io(I,J,K),H_INTP_io(I,J,K), &
                            Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                            uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                            Dprime_io(I,J,K),Mprime_io(I,J,K),udprime_io(I,J,K,1:3)
                        ELSE
                            WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                            U_INTP_io(I,J,K,1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                            uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),0.0_WP, &
                            0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END IF
                    END DO
                END DO
            END DO
        
        END IF
        
        IF(TGFLOWFLG .and. IOFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') ZND(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') ZND(1),YND(NND2)
        
        END IF

        CLOSE(TECFLG)
     
     
    END SUBROUTINE
!**********************************************************************************************************************************
    SUBROUTINE TEC360_YSLICE(N)
        use TEC360_INFO
        IMPLICIT NONE
      
        INTEGER(4) :: I, J, K, N1, N
        CHARACTER(4) :: PNTIM
        INTEGER(4)  :: TECFLG = 202
        character(128) :: FLNM
        logical        :: file_exists
        real(wp)       :: uprime(3)
           
        IF(N.LT.(MGRID/2)) THEN
            J= JGMOV(N)+1
        ELSE
            J= JGMOV(N)
        END IF
        TECFLG = TECFLG + J
        WRITE(PNTIM,'(I4.4)') J
        FLNM = TRIM(filepath5)//'RESULT.UVWP.YSLICE.'//PNTIM//'.tec'
        
        INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
        
        IF(IOFLOWflg .AND. TGFLOWFLG) THEN 
            N1 = NND1_tg+NND1_io
        ELSE
            IF(TGFLOWflg) N1 = NND1_tg
            IF(IOFLOWFLG) N1 = NND1_io
        END IF
        
        IF(file_exists) THEN
        
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)), position='append')
            WRITE(TECFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                          ' ", I=', N1, ', J=',1,', K=',NND3,', F=POINT'
            write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'
            
            DO K=1,NND3
                IF(TGFLOWflg) THEN
                    DO I =1,NND1_tg 
                        WRITE(TECFLG,'(30ES15.7)') U_INTP_tg(I,J,K,1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_tg(I,J,K,1:4), Qcr_tg(I,J,K), delta_tg(I,J,K), lambda2_tg(I,J,K), swirlstrength_tg(I,J,K,1:3),&
                            uprime_tg(I,J,K,1:4),0.0_WP,0.0_WP,0.0_WP,0.0_WP, &
                            0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                    END DO
                END IF
                IF(IOFLOWflg) THEN
                    DO I =1,NND1_io 
                        IF(thermlflg==1) THEN                   
                            WRITE(TECFLG,'(30ES15.7)') &
                            U_INTP_io(I,J,K,1:4), &
                            T_INTP_io(I,J,K),D_INTP_io(I,J,K),H_INTP_io(I,J,K), &
                            Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                            uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                            Dprime_io(I,J,K),Mprime_io(I,J,K),udprime_io(I,J,K,1:3)
                        ELSE
                            WRITE(TECFLG,'(30ES15.7)') &
                            U_INTP_io(I,J,K,1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                            uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),0.0_WP, &
                            0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END IF
                    END DO
                END IF
                !END DO
            END DO
        
        
        ELSE
        
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)) )
            WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW Y-SLICE"'
            WRITE(TECFLG,'(A)',advance="no") 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
            WRITE(TECFLG,'(A)',advance="no") '"Vorx","Vory","Vorz","VorM","Qcr","delta","lambda2",'
            WRITE(TECFLG,'(A)',advance="no") '"swirlstrengthR1","swirlstrengthR2","swirlstrengthI2",'
            WRITE(TECFLG,'(A)',advance="no") '"uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime",'
            WRITE(TECFLG,'(A)')              '"Tprime", "Dprime","Mprime","uvprime", "udprime", "vdprime"'
            
            WRITE(TECFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                          ' ", I=', N1, ', J=',1,', K=',NND3,', F=POINT'
                       
            DO K=1,NND3
                IF(TGFLOWflg) THEN
                    DO I=1,NND1_tg  
                        WRITE(TECFLG,'(33ES15.7)') XND_tg(I),YND(J),ZND(K), &
                        U_INTP_tg(I,J,K,1:4), &
                        1.0_WP, 1.0_WP, 0.0_WP, &
                        Vor_tg(I,J,K,1:4), Qcr_tg(I,J,K), delta_tg(I,J,K), lambda2_tg(I,J,K), swirlstrength_tg(I,J,K,1:3),&
                        uprime_tg(I,J,K,1:4), 0.0_WP,0.0_WP,0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                    END DO
                END IF
                IF(IOFLOWflg) THEN
                    DO I =1,NND1_io  
                        IF(thermlflg==1) THEN                   
                            WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                            U_INTP_io(I,J,K,1:4), &
                            T_INTP_io(I,J,K),D_INTP_io(I,J,K),H_INTP_io(I,J,K), &
                            Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                            uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                            Dprime_io(I,J,K),Mprime_io(I,J,K),udprime_io(I,J,K,1:3)
                        ELSE
                            WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                            U_INTP_io(I,J,K,1:4), &
                            1.0_WP, 1.0_WP, 0.0_WP, &
                            Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                            uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),0.0_WP, &
                            0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END IF
                    END DO
                END IF
                !END DO
            END DO
        
        
        END IF
        
        IF(IOFLOWflg .AND. TGFLOWFLG) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),ZND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),ZND(NND3)
        
        END IF

        CLOSE(TECFLG)
     
     
     END SUBROUTINE TEC360_YSLICE
         
!**********************************************************************************************************************************
    SUBROUTINE TEC360_ZSLICE(N)
        use TEC360_INFO
        IMPLICIT NONE
      
        INTEGER(4)     :: I, J, K, N1,N
        CHARACTER(4)   :: PNTIM
        INTEGER(4)     :: TECFLG = 202
        character(128) :: FLNM
        logical        :: file_exists
        real(wp)       :: uprime(3)
           
        TECFLG = TECFLG + KID(N)
        IF(IOFLOWflg .AND. TGFLOWFLG) THEN 
            N1 = NND1_tg+NND1_io
        ELSE
            IF(TGFLOWflg) N1 = NND1_tg
            IF(IOFLOWFLG) N1 = NND1_io
        END IF
        WRITE(PNTIM,'(I4.4)') KID(N)
        
        FLNM = TRIM(filepath5)//'RESULT.UVWP.ZSLICE.'//PNTIM//'.tec'
        
        INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
        
        IF(file_exists) THEN
        
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)), position='append')
            WRITE(TECFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                          ' ", I=', N1, ', J=',NND2,', K=',1,', F=POINT'
            write(TECFLG,'(A)') 'VARSHARELIST=([1-3]=1)'
            
            DO K=KID(N),KID(N)
                DO J=1,NND2
                    IF(TGFLOWflg) THEN
                        DO I =1,NND1_tg 
                            WRITE(TECFLG,'(30ES15.7)') U_INTP_tg(I,J,K,1),U_INTP_tg(I,J,K,2),U_INTP_tg(I,J,K,3),&
                                U_INTP_tg(I,J,K,4), &
                                1.0_WP, 1.0_WP, 0.0_WP, &
                                Vor_tg(I,J,K,1:4), Qcr_tg(I,J,K), delta_tg(I,J,K), lambda2_tg(I,J,K), swirlstrength_tg(I,J,K,1:3),&
                                uprime_tg(I,J,K,1:4), 0.0_WP,0.0_WP,0.0_WP,0.0_WP, &
                                0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END DO
                    END IF
                    IF(IOFLOWflg) THEN
                        DO I =1,NND1_io 
                        
                            IF(thermlflg==1) THEN            
                                WRITE(TECFLG,'(30ES15.7)') &
                                U_INTP_io(I,J,K,1:4), &
                                T_INTP_io(I,J,K),D_INTP_io(I,J,K),H_INTP_io(I,J,K), &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                                Dprime_io(I,J,K),Mprime_io(I,J,K),udprime_io(I,J,K,1:3)
                            ELSE
                                WRITE(TECFLG,'(30ES15.7)') &
                                U_INTP_io(I,J,K,1:4), &
                                1.0_WP, 1.0_WP, 0.0_WP, &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K), swirlstrength_io(I,J,K,1:3),&
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),0.0_WP, &
                                0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                            END IF
                        END DO
                    END IF
                END DO
            END DO
        
        
        ELSE
        
            OPEN(TECFLG,FILE=TRIM(ADJUSTL(FLNM)) )
            WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW Z-SLICE"'
            WRITE(TECFLG,'(A)',advance="no") 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "P", "T", "D", "H", '
            WRITE(TECFLG,'(A)',advance="no") '"Vorx","Vory","Vorz","VorM","Qcr","delta","lambda2",'
            WRITE(TECFLG,'(A)',advance="no") '"swirlstrengthR1","swirlstrengthR2","swirlstrengthI2",'
            WRITE(TECFLG,'(A)',advance="no") '"uprime","vprime","wprime","pprime","g1prime","g2prime","g3prime",'
            WRITE(TECFLG,'(A)')              '"Tprime", "Dprime","Mprime","uvprime", "udprime", "vdprime"'
            
            WRITE(TECFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T=" ',ITERG,phyTIME, &
                          ' ", I=', N1, ', J=',NND2,', K=',1,', F=POINT'
                       
    
            DO K=KID(N),KID(N)
                DO J=1,NND2
                    IF(TGFLOWflg) THEN
                        DO I=1,NND1_tg    
                            WRITE(TECFLG,'(33ES15.7)') XND_tg(I),YND(J),ZND(K), &
                                U_INTP_tg(I,J,K,1:4), &
                                1.0_WP, 1.0_WP, 0.0_WP, &
                                Vor_tg(I,J,K,1:4), Qcr_tg(I,J,K), delta_tg(I,J,K), lambda2_tg(I,J,K),swirlstrength_tg(I,J,K,1:3),&
                                uprime_tg(I,J,K,1:4), 0.0_WP,0.0_WP,0.0_WP,0.0_WP, &
                                0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                        END DO
                    END IF
                    IF(IOFLOWflg) THEN
                        DO I =1,NND1_io  
                            IF(thermlflg==1) THEN                
                                WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                                U_INTP_io(I,J,K,1:4), &
                                T_INTP_io(I,J,K),D_INTP_io(I,J,K),H_INTP_io(I,J,K), &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K),swirlstrength_io(I,J,K,1:3), &
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),Tprime_io(I,J,K), &
                                Dprime_io(I,J,K),Mprime_io(I,J,K),udprime_io(I,J,K,1:3)
                            ELSE
                                WRITE(TECFLG,'(33ES15.7)') XND_io(I),YND(J),ZND(K), &
                                U_INTP_io(I,J,K,1:4), &
                                1.0_WP, 1.0_WP, 0.0_WP, &
                                Vor_io(I,J,K,1:4), Qcr_io(I,J,K), delta_io(I,J,K), lambda2_io(I,J,K),swirlstrength_io(I,J,K,1:3),&
                                uprime_io(I,J,K,1:4), gprime_io(I,J,K,1:3),0.0_WP, &
                                0.0_WP,0.0_WP,udprime_io(I,J,K,1),0.0_WP,0.0_WP
                            END IF
                        END DO
                    END IF
               END DO
            END DO
        
        END IF
        
        IF(IOFLOWflg .AND. TGFLOWflg) THEN
            WRITE(TECFLG,'(A)') 'GEOMETRY X = 0.0, Y = 0.0, T = LINE, CS = GRID,'// &
                                 'L = SOLID, LT = 0.005, C = BLACK, FC = BLACK,'// &
                                 'F = POINT, S = GLOBAL'
            WRITE(TECFLG,'(I2.1)') 1
            WRITE(TECFLG,'(I2.1)') 2
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(1)
            WRITE(TECFLG,'(2ES15.7)') XND_io(1),YND(NND2)
        
        END IF

        CLOSE(TECFLG)
     
     
     END SUBROUTINE
     
     
    SUBROUTINE TEC360_INSTANT_UPRIME
        use TEC360_INFO
        IMPLICIT NONE
        
        INTEGER(4) I, J, K
        
        IF(TGFLOWflg) THEN
            ALLOCATE ( uprime_tg(NND1_io,NND2,NND3,NDV+1) )
            DO I=1,NND1_tg
                DO K=1, NND3
                    DO J=1, NND2
                        uprime_tg(I,J,K,1:4) = U_INTP_tg(I,J,K,1:4) - U1xztL_INTP_tg(J,1:4)
                    END DO
                END DO
            END DO
        END IF
        
        IF(IOFLOWFLG) THEN
            ALLOCATE ( uprime_io(NND1_io,NND2,NND3,NDV+1) )
            ALLOCATE ( udprime_io(NND1_io,NND2,NND3,NDV) )

            IF(TGFLOWflg) THEN
                DO I=1,NND1_io
                    DO K=1, NND3
                        DO J=1, NND2
                            uprime_io(I,J,K,1:4) = U_INTP_io(I,J,K,1:4) - U1ztL_INTP_io(I,J,1:4)
                        END DO
                    END DO
                END DO
                udprime_io(:,:,:,1)=uprime_io(:,:,:,1)*uprime_io(:,:,:,2) !u'v'
                
                IF(thermlflg==1) THEN
                    ALLOCATE ( gprime_io(NND1_io,NND2,NND3,NDV) )
                    ALLOCATE ( Tprime_io (NND1_io,NND2,NND3) )
                    ALLOCATE ( dprime_io (NND1_io,NND2,NND3) )
                    ALLOCATE ( mprime_io (NND1_io,NND2,NND3) )
                    
                    DO I=1,NND1_io
                        DO K=1, NND3
                            DO J=1, NND2
                                !gprime_io(I,J,K,1:3) = U_INTP_io(I,J,K,1:3) - G1ztL_INTP_io(I,J,1:3)/D1ztL_INTP_io(I,J)
                                gprime_io(I,J,K,1:3) = G_INTP_io(I,J,K,1:3) - G1ztL_INTP_io(I,J,1:3)
                                Tprime_io(I,J,K)      = T_INTP_io(I,J,K) -     T1ztL_INTP_io(I,J)
                                dprime_io(I,J,K)      = D_INTP_io(I,J,K) -     D1ztL_INTP_io(I,J)
                            END DO
                        END DO
                    END DO
                    udprime_io(:,:,:,2)=uprime_io(:,:,:,1)*dprime_io(:,:,:) !u'd'
                    udprime_io(:,:,:,3)=uprime_io(:,:,:,2)*dprime_io(:,:,:) !v'd'
                ELSE
                    gprime_io = uprime_io
                END IF
                
            ELSE  
                DO I=1,NND1_io
                    DO K=1, NND3
                        DO J=1, NND2
                            uprime_io(I,J,K,1:4) = U_INTP_io(I,J,K,1:4) - U1xztL_INTP_io(J,1:4)
                        END DO
                    END DO
                END DO
                udprime_io(:,:,:,1)=uprime_io(:,:,:,1)*uprime_io(:,:,:,2) !u'v'
                
                IF(thermlflg==1) THEN
                    ALLOCATE ( gprime_io(NND1_io,NND2,NND3,NDV) )
                    ALLOCATE ( Tprime_io (NND1_io,NND2,NND3) )
                    ALLOCATE ( Dprime_io (NND1_io,NND2,NND3) )
                    ALLOCATE ( Mprime_io (NND1_io,NND2,NND3) )
                    DO I=1,NND1_io
                        DO K=1, NND3
                            DO J=1, NND2
                                !gprime_io(I,J,K,1:3) = U_INTP_io(I,J,K,1:3) - G1xztL_INTP_io(J,1:3)/D1xztL_INTP_io(J)
                                gprime_io(I,J,K,1:3)  = G_INTP_io(I,J,K,1:3) - G1xztL_INTP_io(J,1:3)
                                Tprime_io(I,J,K)      = T_INTP_io(I,J,K) -     T1xztL_INTP_io(J)
                                dprime_io(I,J,K)      = D_INTP_io(I,J,K) -     D1xztL_INTP_io(J)
                                mprime_io(I,J,K)      = M_INTP_io(I,J,K) -     M1xztL_INTP_io(J)
                            END DO
                        END DO
                    END DO
                    udprime_io(:,:,:,2)=uprime_io(:,:,:,1)*dprime_io(:,:,:) !u'd'
                    udprime_io(:,:,:,3)=uprime_io(:,:,:,2)*dprime_io(:,:,:) !v'd'
                ELSE
                    gprime_io = uprime_io
                END IF
            END IF

        END IF
        
        
        
        RETURN
    END SUBROUTINE
     
    SUBROUTINE TEC360_INSTANT_VORTEX_CRITERIA
        use TEC360_INFO
        IMPLICIT NONE
        
        integer(4)  :: IC, JC, KC
        integer(4)  :: IM, JM, KM
        integer(4)  :: IP, JP, KP
        real(wp)    :: DUDX1, DUDX2, DUDX3, DUDX4, DUDX5, DUDX6, DUDX7, DUDX8 
        real(wp)    :: DVDY1, DVDY2, DVDY3, DVDY4, DVDY5, DVDY6, DVDY7, DVDY8 
        real(wp)    :: DWDZ1, DWDZ2, DWDZ3, DWDZ4, DWDZ5, DWDZ6, DWDZ7, DWDZ8
        real(wp)    :: DUDY1, DUDY2, DUDZ1, DUDZ2
        real(wp)    :: DWDY1, DWDY2, DVDZ1, DVDZ2
        real(wp)    :: DWDX1, DWDX2, DVDX1, DVDX2
        real(wp)    :: strSSYM(3,3)  !symmetric components of velcocity gradient 
        real(wp)    :: strASYM(3,3)  !antisymmetric components of velcocity gradient 
        real(wp)    :: gradU(3,3)
        real(wp)    :: S2plusO2(3,3)
        real(wp)    :: eig(3)
        real(wp)    :: detGradU 
        real(wp)    :: traceGradU
        real(wp)    :: STSSYM
        real(wp)    :: STASYM
        real(wp), allocatable :: dQ_io(:,:,:,:,:)
        real(wp), allocatable :: dQ_tg(:,:,:,:,:)
        
        real(wp)    :: ar(3,3),ai(3,3)
        real(wp)    :: wr(3), wi(3)
        real(wp)    :: zr(3,3), zi(3,3)
        integer(4)  :: nm, nn
        integer(4)  :: matz
        integer(4)  :: ierr
        real(wp)    :: fv1(3), fv2(3), fv3(3)
    
        
        IF(IOFLOWFLG) THEN
            allocate (Qcr_io    (NND1_io,NND2,NND3)       )
            allocate (vor_io    (NND1_io,NND2,NND3,4)     )
            allocate (delta_io  (NND1_io,NND2,NND3)       )
            allocate (lambda2_io(NND1_io,NND2,NND3)       )
            allocate (swirlstrength_io(NND1_io,NND2,NND3,3) )
            
            allocate (dQ_io    (NND1_io,NND2,NND3,3,3)   )
            
            !====================GET ALL du_i/dx_i on nodes===================================================
            DO IC=1,NCL1_io
                IM = IMV_IO(IC)
                IP = IPV_IO(IC)
                DO JC=2,NCL2
                    JM=JGMV(JC)
                    JP=JGPV(JC)
                    DO KC=1,NCL3
                        KP=KPV(KC)
                        KM=KMV(KC)
                        !===============================================================================
                        !===du/dx at i', j', k'============================
                        DUDX1 = (U_F0_io(IC,JC,KM,1)-U_F0_io(IM,JC,KM,1))*DXI       ! (i-1, j, k-1)
                        DUDX2 = (U_F0_io(IP,JC,KM,1)-U_F0_io(IC,JC,KM,1))*DXI       ! (i,   j, k-1)
                        
                        DUDX3 = (U_F0_io(IC,JC,KC,1)-U_F0_io(IM,JC,KC,1))*DXI       ! (i-1, j, k)
                        DUDX4 = (U_F0_io(IP,JC,KC,1)-U_F0_io(IC,JC,KC,1))*DXI       ! (i,   j, k)
                        
                        DUDX5 = (U_F0_io(IC,JM,KM,1)-U_F0_io(IM,JM,KM,1))*DXI       ! (i-1, j-1, k-1)
                        DUDX6 = (U_F0_io(IP,JM,KM,1)-U_F0_io(IC,JM,KM,1))*DXI       ! (i,   j-1, k-1)
                        
                        DUDX7 = (U_F0_io(IC,JM,KC,1)-U_F0_io(IM,JM,KC,1))*DXI       ! (i-1, j-1, k)
                        DUDX8 = (U_F0_io(IP,JM,KC,1)-U_F0_io(IC,JM,KC,1))*DXI       ! (i,   j-1, k)
                        
                        dQ_io(IC,JC,KC,1,1) = ( (DUDX1+DUDX2)*0.5_wp + (DUDX3+DUDX4)*0.5_wp )*0.5_wp*YCL2ND_WFF(JC) + &
                                              ( (DUDX3+DUDX4)*0.5_wp + (DUDX5+DUDX6)*0.5_wp )*0.5_wp*YCL2ND_WFB(JC) 
    
                        !===du/dy at i',j',k'============================
                        DUDY1 = ( (U_F0_io(IC,JC,KM,1)-U_F0_io(IC,JM,KM,1))- (U1xztL_F0_io(JC,1)-U1xztL_F0_io(JM,1)) )*DYCI(JC)  ! (i',j',k-1)
                        DUDY2 = ( (U_F0_io(IC,JC,KC,1)-U_F0_io(IC,JM,KC,1))- (U1xztL_F0_io(JC,1)-U1xztL_F0_io(JM,1)) )*DYCI(JC)  ! (i',j',k)
                        
                        dQ_io(IC,JC,KC,1,2) = 0.5_wp*(DUDY1+DUDY2)                 ! (i',j',k')
                        
                        !===du/dz at i',j',k'============================
                        DUDZ1 = (U_F0_io(IC,JM,KC,1)-U_F0_io(IC,JM,KM,1))*DZI       ! (i',j-1,k')
                        DUDZ2 = (U_F0_io(IC,JC,KC,1)-U_F0_io(IC,JC,KM,1))*DZI       ! (i',j,  k')
                        
                        dQ_io(IC,JC,KC,1,3) = DUDZ1*YCL2ND_WFB(JC)+DUDZ2*YCL2ND_WFF(JC)
                        
                        !===============================================================================
                        !===dv/dy at i', j', k'============================
                        IF(JC==NCL2) THEN
                            DVDY1 = ( (0.0_wp-U_F0_io(IM,JC,KM,2))- (0.0_WP-U1xztL_F0_io(JC,2)) )*DYFI(JC)       ! (i-1, j, k-1)
                            DVDY2 = ( (0.0_wp-U_F0_io(IC,JC,KM,2))- (0.0_WP-U1xztL_F0_io(JC,2)) )*DYFI(JC)       ! (i,   j, k-1)
                        
                            DVDY3 = (0.0_wp-U_F0_io(IM,JC,KC,2))*DYFI(JC)       ! (i-1, j, k)
                            DVDY4 = (0.0_wp-U_F0_io(IC,JC,KC,2))*DYFI(JC)       ! (i,   j, k)
                        ELSE
                            DVDY1 = ( (U_F0_io(IM,JP,KM,2)-U_F0_io(IM,JC,KM,2))- (U1xztL_F0_io(JP,2)-U1xztL_F0_io(JC,2)) )*DYFI(JC)       ! (i-1, j, k-1)
                            DVDY2 = ( (U_F0_io(IC,JP,KM,2)-U_F0_io(IC,JC,KM,2))- (U1xztL_F0_io(JP,2)-U1xztL_F0_io(JC,2)) )*DYFI(JC)       ! (i,   j, k-1)
                        
                            DVDY3 = ( (U_F0_io(IM,JP,KC,2)-U_F0_io(IM,JC,KC,2))- (U1xztL_F0_io(JP,2)-U1xztL_F0_io(JC,2)) )*DYFI(JC)       ! (i-1, j, k)
                            DVDY4 = ( (U_F0_io(IC,JP,KC,2)-U_F0_io(IC,JC,KC,2))- (U1xztL_F0_io(JP,2)-U1xztL_F0_io(JC,2)) )*DYFI(JC)       ! (i,   j, k)
                        END IF

                        DVDY5 = ( (U_F0_io(IM,JC,KM,2)-U_F0_io(IM,JM,KM,2))- (U1xztL_F0_io(JC,2)-U1xztL_F0_io(JM,2)) )*DYFI(JM)      ! (i-1, j-1, k-1)
                        DVDY6 = ( (U_F0_io(IC,JC,KM,2)-U_F0_io(IC,JM,KM,2))- (U1xztL_F0_io(JC,2)-U1xztL_F0_io(JM,2)) )*DYFI(JM)      ! (i,   j-1, k-1)
                        
                        DVDY7 = ( (U_F0_io(IM,JC,KC,2)-U_F0_io(IM,JM,KC,2))- (U1xztL_F0_io(JC,2)-U1xztL_F0_io(JM,2)) )*DYFI(JM)      ! (i-1, j-1, k)
                        DVDY8 = ( (U_F0_io(IC,JC,KC,2)-U_F0_io(IC,JM,KC,2))- (U1xztL_F0_io(JC,2)-U1xztL_F0_io(JM,2)) )*DYFI(JM)      ! (i,   j-1, k)
                        
                        dQ_io(IC,JC,KC,2,2) = ( (DVDY1+DVDY2)*0.5_wp + (DVDY3+DVDY4)*0.5_wp )*0.5_wp*YCL2ND_WFF(JC) + &
                                              ( (DVDY3+DVDY4)*0.5_wp + (DVDY5+DVDY6)*0.5_wp )*0.5_wp*YCL2ND_WFB(JC) 
    
                        !===dv/dx at i',j',k'============================
                        DVDX1 = (U_F0_io(IC,JC,KM,2)-U_F0_io(IM,JC,KM,2))*DXI  ! (i',j',k-1)
                        DVDX2 = (U_F0_io(IC,JC,KC,2)-U_F0_io(IM,JC,KC,2))*DXI  ! (i',j',k)
                        
                        dQ_io(IC,JC,KC,2,1) = 0.5_wp*(DVDX1+DVDX2)                 ! (i',j',k')
                        
                        !===dv/dz at i',j',k'============================
                        DVDZ1 = (U_F0_io(IM,JC,KC,2)-U_F0_io(IM,JC,KM,2))*DZI       ! (i-1,j',k')
                        DVDZ2 = (U_F0_io(IC,JC,KC,2)-U_F0_io(IC,JC,KM,2))*DZI       ! (i,  j',k')
                        
                        dQ_io(IC,JC,KC,2,3) = (DVDZ1+DVDZ2)*0.5_WP
                        
                        !===============================================================================
                        !===dw/dz at i', j', k'============================
                        DWDZ1 = (U_F0_io(IM,JC,KC,3)-U_F0_io(IM,JC,KM,3))*DXI       ! (i-1, j, k-1)
                        DWDZ2 = (U_F0_io(IC,JC,KC,3)-U_F0_io(IC,JC,KM,3))*DXI       ! (i,   j, k-1)
                        
                        DWDZ3 = (U_F0_io(IM,JC,KP,3)-U_F0_io(IM,JC,KC,3))*DXI       ! (i-1, j, k)
                        DWDZ4 = (U_F0_io(IC,JC,KP,3)-U_F0_io(IC,JC,KC,3))*DXI       ! (i,   j, k)
                        
                        DWDZ5 = (U_F0_io(IM,JM,KC,3)-U_F0_io(IM,JM,KM,3))*DXI       ! (i-1, j-1, k-1)
                        DWDZ6 = (U_F0_io(IC,JM,KC,3)-U_F0_io(IC,JM,KM,3))*DXI       ! (i,   j-1, k-1)
                        
                        DWDZ7 = (U_F0_io(IM,JM,KP,3)-U_F0_io(IM,JM,KC,3))*DXI       ! (i-1, j-1, k)
                        DWDZ8 = (U_F0_io(IC,JM,KP,3)-U_F0_io(IC,JM,KC,3))*DXI       ! (i,   j-1, k)
                        
                        dQ_io(IC,JC,KC,3,3) = ( (DWDZ1+DWDZ2)*0.5_wp + (DWDZ3+DWDZ4)*0.5_wp )*0.5_wp*YCL2ND_WFF(JC) + &
                                           ( (DWDZ3+DWDZ4)*0.5_wp + (DWDZ5+DWDZ6)*0.5_wp )*0.5_wp*YCL2ND_WFB(JC) 
    
                        !===dw/dy at i',j',k'============================
                        DWDY1 = ( (U_F0_io(IM,JC,KC,3)-U_F0_io(IM,JM,KC,3))- (U1xztL_F0_io(JC,3)-U1xztL_F0_io(JM,3)) )*DYCI(JC)  ! (i-1,j',k')
                        DWDY2 = ( (U_F0_io(IC,JC,KC,3)-U_F0_io(IC,JM,KC,3))- (U1xztL_F0_io(JC,3)-U1xztL_F0_io(JM,3)) )*DYCI(JC)  ! (i,  j',k)
                        
                        dQ_io(IC,JC,KC,3,2) = 0.5_wp*(DWDY1+DWDY2)                 ! (i',j',k')
                        
                        !===dw/dx at i',j',k'============================
                        DWDX1 = (U_F0_io(IC,JM,KC,3)-U_F0_io(IM,JM,KC,3))*DZI       ! (i',j-1,k')
                        DWDX2 = (U_F0_io(IC,JC,KC,3)-U_F0_io(IM,JC,KC,3))*DZI       ! (i',j,  k')
                        
                        dQ_io(IC,JC,KC,3,1) = DWDX1*YCL2ND_WFB(JC)+DWDX2*YCL2ND_WFF(JC)
                    END DO
                END DO
            END DO
            dQ_io(:,NND2,:,:,:)=2.0_wp*dQ_io(:,NCL2,:,:,:)-dQ_io(:,NCL2-1,:,:,:)
            dQ_io(:,1,   :,:,:)=2.0_wp*dQ_io(:,2,   :,:,:)-dQ_io(:,3,     :,:,:)
            IF(.NOT.TGFLOWFLG) dQ_io(NND1_io,:,:,   :,:)=dQ_io(1,:,:,:,:)
            dQ_io(:,      :,NND3,:,:)=dQ_io(:,:,1,:,:)
            
            
            
            DO IC=1,NND1_IO
                DO JC=1,NND2
                    DO KC=1,NND3
                        !===================vorticity===============================
                        VOR_io(IC,JC,KC,1) = dQ_io(IC,JC,KC,3,2) - dQ_io(IC,JC,KC,2,3) ! VOR_io-x 
                        VOR_io(IC,JC,KC,2) = dQ_io(IC,JC,KC,1,3) - dQ_io(IC,JC,KC,3,1) ! VOR_io-y
                        VOR_io(IC,JC,KC,3) = dQ_io(IC,JC,KC,2,1) - dQ_io(IC,JC,KC,1,2) ! VOR_io-z
                        VOR_io(IC,JC,KC,4) = DSQRT( DOT_PRODUCT( VOR_io(IC,JC,KC,1:3),VOR_io(IC,JC,KC,1:3) ) )!VOR_io-mag
                        
                        !==========================================================================================
                        !===================omega = antisymetric components of graduate U===============
                        strASYM(1,1) = dQ_io(IC,JC,KC,1,1) - dQ_io(IC,JC,KC,1,1)  !dU/dX - dU/dX
                        strASYM(1,2) = dQ_io(IC,JC,KC,2,1) - dQ_io(IC,JC,KC,1,2)  !dV/dX - dU/dY
                        strASYM(1,3) = dQ_io(IC,JC,KC,3,1) - dQ_io(IC,JC,KC,1,3)  !dW/dX - dU/dZ
                    
                        strASYM(2,1) = dQ_io(IC,JC,KC,1,2) - dQ_io(IC,JC,KC,2,1)  !dU/dY - dV/dX
                        strASYM(2,2) = dQ_io(IC,JC,KC,2,2) - dQ_io(IC,JC,KC,2,2)  !dV/dY - dV/dY
                        strASYM(2,3) = dQ_io(IC,JC,KC,3,2) - dQ_io(IC,JC,KC,2,3)  !dW/dY - dV/dZ
                    
                        strASYM(3,1) = dQ_io(IC,JC,KC,1,3) - dQ_io(IC,JC,KC,3,1)  !dU/dZ - dW/dX
                        strASYM(3,2) = dQ_io(IC,JC,KC,2,3) - dQ_io(IC,JC,KC,3,2)  !dV/dZ - dW/dY
                        strASYM(3,3) = dQ_io(IC,JC,KC,3,3) - dQ_io(IC,JC,KC,3,3)  !dW/dZ - dW/dZ
                        !===================omega = symetric components of graduate U===============
                        strSSYM(1,1) = dQ_io(IC,JC,KC,1,1) + dQ_io(IC,JC,KC,1,1)  !dU/dX + dU/dX
                        strSSYM(1,2) = dQ_io(IC,JC,KC,2,1) + dQ_io(IC,JC,KC,1,2)  !dV/dX + dU/dY
                        strSSYM(1,3) = dQ_io(IC,JC,KC,3,1) + dQ_io(IC,JC,KC,1,3)  !dW/dX + dU/dZ
                    
                        strSSYM(2,1) = dQ_io(IC,JC,KC,1,2) + dQ_io(IC,JC,KC,2,1)  !dU/dY + dV/dX
                        strSSYM(2,2) = dQ_io(IC,JC,KC,2,2) + dQ_io(IC,JC,KC,2,2)  !dV/dY + dV/dY
                        strSSYM(2,3) = dQ_io(IC,JC,KC,3,2) + dQ_io(IC,JC,KC,2,3)  !dW/dY + dV/dZ
                    
                        strSSYM(3,1) = dQ_io(IC,JC,KC,1,3) + dQ_io(IC,JC,KC,3,1)  !dU/dZ + dW/dX
                        strSSYM(3,2) = dQ_io(IC,JC,KC,2,3) + dQ_io(IC,JC,KC,3,2)  !dV/dZ + dW/dY
                        strSSYM(3,3) = dQ_io(IC,JC,KC,3,3) + dQ_io(IC,JC,KC,3,3)  !dW/dZ + dW/dZ
                        
                        strASYM  = strASYM  * 0.5_wp  ! Sij
                        strSSYM  = strSSYM  * 0.5_wp  ! wij
                        
                        ! trace of (STSYM*STSYM^T)
                        STSSYM   =  DOT_PRODUCT(strSSYM(:,1), strSSYM(:,1)) + &
                                    DOT_PRODUCT(strSSYM(:,2), strSSYM(:,2)) + &
                                    DOT_PRODUCT(strSSYM(:,3), strSSYM(:,3))    ! Sij*Sij
                        ! trace of (STASYM*STASYM^T)     
                        STASYM   = DOT_PRODUCT(strASYM(:,1), strASYM(:,1)) + &
                                   DOT_PRODUCT(strASYM(:,2), strASYM(:,2)) + &
                                   DOT_PRODUCT(strASYM(:,3), strASYM(:,3))    ! wij*wij
                                   
                        !================================Q cretiria===================================================
                        QCR_io(IC,JC,KC)=0.5_wp*(STASYM-STSSYM) ! Q
                        
                        
                        !================================delta_io criteria===============================================
                        ! grad U 
                        gradU   = strASYM + strSSYM
                        
                        !  negtive trace of velocity gradient  !P
                        traceGradU = (strASYM(1,1)+strASYM(2,2)+strASYM(3,3)+&
                                      strSSYM(1,1)+strSSYM(2,2)+strSSYM(3,3))*(-1.0_WP)
                        
                        ! det(gradU) !R
                        detGradU =(gradU(1,1) * ( gradU(2,2)*gradU(3,3) - gradU(2,3)*gradU(3,2) ) - &
                                   gradU(1,2) * ( gradU(2,1)*gradU(3,3) - gradU(3,1)*gradU(2,3) ) + &
                                   gradU(1,3) * ( gradU(2,1)*gradU(3,2) - gradU(2,2)*gradU(3,1) ))*(-1.0_WP)
                        ! delta for incompressible flow
                        !delta_io(IC,JC,KC) = (QCR_io(IC,JC,KC)/3.0_wp)**3 + (detGradU/2.0_wp)**2
                        
                        ! delta for compressible flow
                        delta_io(IC,JC,KC) = 27.0_WP*detGradU*detGradU + &
                                (4.0_wp*traceGradU*traceGradU*traceGradU-18.0_WP*traceGradU*QCR_io(IC,JC,KC))*detGradU + &
                                4.0_wp*QCR_io(IC,JC,KC)*QCR_io(IC,JC,KC)*QCR_io(IC,JC,KC)- &
                                traceGradU*traceGradU*QCR_io(IC,JC,KC)*QCR_io(IC,JC,KC)             
                        
                        !get rid of possible infinity 
                        if (delta_io(IC,JC,KC)> 1.0E30_wp) delta_io(IC,JC,KC)= 1.0E30_wp
                        if (delta_io(IC,JC,KC)<-1.0E30_wp) delta_io(IC,JC,KC)=-1.0E30_wp
             
                        !==============swirl strength ==================================
                        if (delta_io(IC,JC,KC)>0.0_WP) then
                        
                            matz = 0
                            nm   = 3
                            nn   = 3
                            ar(:,:)   = strASYM(:,:) + strSSYM(:,:)
                            ai(:,:)   = 0.0
                            
                            call cg(nm,nn,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
                            
                            if(dabs(wi(1))<(dabs(wi(2))+1.0E-14) .and. dabs(wi(1))<(dabs(wi(3))+1.0E-14) ) then
                                swirlstrength_io(IC,JC,KC,1) = wr(1)
                                swirlstrength_io(IC,JC,KC,2) = wr(2)
                                swirlstrength_io(IC,JC,KC,3) = wr(2)
                            else if (dabs(wi(2))<(dabs(wi(1))+1.0E-14) .and. dabs(wi(2))<(dabs(wi(3))+1.0E-14) ) then
                                swirlstrength_io(IC,JC,KC,1) = wr(2)
                                swirlstrength_io(IC,JC,KC,2) = wr(1)
                                swirlstrength_io(IC,JC,KC,3) = wr(1)
                            else if (dabs(wi(3))<(dabs(wi(1))+1.0E-14) .and. dabs(wi(3))<(dabs(wi(2))+1.0E-14) ) then
                                swirlstrength_io(IC,JC,KC,1) = wr(3)
                                swirlstrength_io(IC,JC,KC,2) = wr(1)
                                swirlstrength_io(IC,JC,KC,3) = wr(1)
                            else
                            end if
 
                        end if
             
             
                        !*************lamda2 criteria *************************
                        S2plusO2(1,1) = strSSYM(1,1)*strSSYM(1,1) + strSSYM(1,2)*strSSYM(2,1) + strSSYM(1,3)*strSSYM(3,1) + &
                                        strASYM(1,1)*strASYM(1,1) + strASYM(1,2)*strASYM(2,1) + strASYM(1,3)*strASYM(3,1)
                        S2plusO2(1,2) = strSSYM(1,1)*strSSYM(1,2) + strSSYM(1,2)*strSSYM(2,2) + strSSYM(1,3)*strSSYM(3,2) + &
                                        strASYM(1,1)*strASYM(1,2) + strASYM(1,2)*strASYM(2,2) + strASYM(1,3)*strASYM(3,2)
                        S2plusO2(1,3) = strSSYM(1,1)*strSSYM(1,3) + strSSYM(1,2)*strSSYM(2,3) + strSSYM(1,3)*strSSYM(3,3) + &
                                        strASYM(1,1)*strASYM(1,3) + strASYM(1,2)*strASYM(2,3) + strASYM(1,3)*strASYM(3,3)
                             
                        S2plusO2(2,1) = strSSYM(2,1)*strSSYM(1,1) + strSSYM(2,2)*strSSYM(2,1) + strSSYM(2,3)*strSSYM(3,1) + &
                                        strASYM(2,1)*strASYM(1,1) + strASYM(2,2)*strASYM(2,1) + strASYM(2,3)*strASYM(3,1)                 
                        S2plusO2(2,2) = strSSYM(2,1)*strSSYM(1,2) + strSSYM(2,2)*strSSYM(2,2) + strSSYM(2,3)*strSSYM(3,2) + &
                                        strASYM(2,1)*strASYM(1,2) + strASYM(2,2)*strASYM(2,2) + strASYM(2,3)*strASYM(3,2)           
                        S2plusO2(2,3) = strSSYM(2,1)*strSSYM(1,3) + strSSYM(2,2)*strSSYM(2,3) + strSSYM(2,3)*strSSYM(3,3) + &
                                        strASYM(2,1)*strASYM(1,3) + strASYM(2,2)*strASYM(2,3) + strASYM(2,3)*strASYM(3,3)
    
                        S2plusO2(3,1) = strSSYM(3,1)*strSSYM(1,1) + strSSYM(3,2)*strSSYM(2,1) + strSSYM(3,3)*strSSYM(3,1) + &
                                        strASYM(3,1)*strASYM(1,1) + strASYM(3,2)*strASYM(2,1) + strASYM(3,3)*strASYM(3,1)                        
                        S2plusO2(3,2) = strSSYM(3,1)*strSSYM(1,2) + strSSYM(3,2)*strSSYM(2,2) + strSSYM(3,3)*strSSYM(3,2) + &
                                        strASYM(3,1)*strASYM(1,2) + strASYM(3,2)*strASYM(2,2) + strASYM(3,3)*strASYM(3,2)                       
                        S2plusO2(3,3) = strSSYM(3,1)*strSSYM(1,3) + strSSYM(3,2)*strSSYM(2,3) + strSSYM(3,3)*strSSYM(3,3) + &
                                        strASYM(3,1)*strASYM(1,3) + strASYM(3,2)*strASYM(2,3) + strASYM(3,3)*strASYM(3,3)
    
                        call DSYEVC3(S2plusO2,eig)
                        if (((eig(1)-eig(2))*(eig(1)-eig(3)))<0.0_wp) lambda2_io(IC,JC,KC)=eig(1)
                        if (((eig(2)-eig(1))*(eig(2)-eig(3)))<0.0_wp) lambda2_io(IC,JC,KC)=eig(2)
                        if (((eig(3)-eig(1))*(eig(3)-eig(2)))<0.0_wp) lambda2_io(IC,JC,KC)=eig(3)
    
               
                        !get rid of possible infinity
                        if (lambda2_io(IC,JC,KC)> 1.0E30_wp) lambda2_io(IC,JC,KC)= 1.0E30_wp
                        if (lambda2_io(IC,JC,KC)<-1.0E30_wp) lambda2_io(IC,JC,KC)=-1.0E30_wp
                        !***********************************************************************************
                    END DO
                END DO
            END DO
            
            deallocate (dQ_io)
        END IF
        
        IF(TGFLOWFLG)THEN
            allocate (Qcr_tg    (NND1_tg,NND2,NND3)       )
            allocate (vor_tg    (NND1_tg,NND2,NND3,4)     )
            allocate (delta_tg  (NND1_tg,NND2,NND3)       )
            allocate (lambda2_tg(NND1_tg,NND2,NND3)       )
            allocate (swirlstrength_tg(NND1_tg,NND2,NND3,3) )
            allocate (dQ_tg     (NND1_tg,NND2,NND3,3,3)   )
            
           
            !====================GET ALL du_i/dx_i on nodes===================================================
            DO IC=1,NCL1_tg
                IM = IMV_tg(IC)
                IP = IPV_tg(IC)
                DO JC=2,NCL2
                    JM=JGMV(JC)
                    JP=JGPV(JC)
                    DO KC=1,NCL3
                        KP=KPV(KC)
                        KM=KMV(KC)
                        !===============================================================================
                        !===du/dx at i', j', k'============================
                        DUDX1 = (U_F0_tg(IC,JC,KM,1)-U_F0_tg(IM,JC,KM,1))*DXI       ! (i-1, j, k-1)
                        DUDX2 = (U_F0_tg(IP,JC,KM,1)-U_F0_tg(IC,JC,KM,1))*DXI       ! (i,   j, k-1)
                        
                        DUDX3 = (U_F0_tg(IC,JC,KC,1)-U_F0_tg(IM,JC,KC,1))*DXI       ! (i-1, j, k)
                        DUDX4 = (U_F0_tg(IP,JC,KC,1)-U_F0_tg(IC,JC,KC,1))*DXI       ! (i,   j, k)
                        
                        DUDX5 = (U_F0_tg(IC,JM,KM,1)-U_F0_tg(IM,JM,KM,1))*DXI       ! (i-1, j-1, k-1)
                        DUDX6 = (U_F0_tg(IP,JM,KM,1)-U_F0_tg(IC,JM,KM,1))*DXI       ! (i,   j-1, k-1)
                        
                        DUDX7 = (U_F0_tg(IC,JM,KC,1)-U_F0_tg(IM,JM,KC,1))*DXI       ! (i-1, j-1, k)
                        DUDX8 = (U_F0_tg(IP,JM,KC,1)-U_F0_tg(IC,JM,KC,1))*DXI       ! (i,   j-1, k)
                        
                        dQ_tg(IC,JC,KC,1,1) = ( (DUDX1+DUDX2)*0.5_wp + (DUDX3+DUDX4)*0.5_wp )*0.5_wp*YCL2ND_WFF(JC) + &
                                           ( (DUDX3+DUDX4)*0.5_wp + (DUDX5+DUDX6)*0.5_wp )*0.5_wp*YCL2ND_WFB(JC) 
    
                        !===du/dy at i',j',k'============================
                        DUDY1 = (U_F0_tg(IC,JC,KM,1)-U_F0_tg(IC,JM,KM,1))*DYCI(JC)  ! (i',j',k-1)
                        DUDY2 = (U_F0_tg(IC,JC,KC,1)-U_F0_tg(IC,JM,KC,1))*DYCI(JC)  ! (i',j',k)
                        
                        dQ_tg(IC,JC,KC,1,2) = 0.5_wp*(DUDY1+DUDY2)                 ! (i',j',k')
                        
                        !===du/dz at i',j',k'============================
                        DUDZ1 = (U_F0_tg(IC,JM,KC,1)-U_F0_tg(IC,JM,KM,1))*DZI       ! (i',j-1,k')
                        DUDZ2 = (U_F0_tg(IC,JC,KC,1)-U_F0_tg(IC,JC,KM,1))*DZI       ! (i',j,  k')
                        
                        dQ_tg(IC,JC,KC,1,3) = DUDZ1*YCL2ND_WFB(JC)+DUDZ2*YCL2ND_WFF(JC)
                        
                        !===============================================================================
                        !===dv/dy at i', j', k'============================
                        IF(JC==NCL2) THEN
                            DVDY1 = (0.0_wp-U_F0_tg(IM,JC,KM,2))*DYFI(JC)       ! (i-1, j, k-1)
                            DVDY2 = (0.0_wp-U_F0_tg(IC,JC,KM,2))*DYFI(JC)       ! (i,   j, k-1)
                        
                            DVDY3 = (0.0_wp-U_F0_tg(IM,JC,KC,2))*DYFI(JC)       ! (i-1, j, k)
                            DVDY4 = (0.0_wp-U_F0_tg(IC,JC,KC,2))*DYFI(JC)       ! (i,   j, k)
                        ELSE 
                            DVDY1 = (U_F0_tg(IM,JP,KM,2)-U_F0_tg(IM,JC,KM,2))*DYFI(JC)       ! (i-1, j, k-1)
                            DVDY2 = (U_F0_tg(IC,JP,KM,2)-U_F0_tg(IC,JC,KM,2))*DYFI(JC)       ! (i,   j, k-1)
                        
                            DVDY3 = (U_F0_tg(IM,JP,KC,2)-U_F0_tg(IM,JC,KC,2))*DYFI(JC)       ! (i-1, j, k)
                            DVDY4 = (U_F0_tg(IC,JP,KC,2)-U_F0_tg(IC,JC,KC,2))*DYFI(JC)       ! (i,   j, k)
                        END IF
                        
                        DVDY5 = (U_F0_tg(IM,JC,KM,2)-U_F0_tg(IM,JM,KM,2))*DYFI(JM)      ! (i-1, j-1, k-1)
                        DVDY6 = (U_F0_tg(IC,JC,KM,2)-U_F0_tg(IC,JM,KM,2))*DYFI(JM)      ! (i,   j-1, k-1)
                        
                        DVDY7 = (U_F0_tg(IM,JC,KC,2)-U_F0_tg(IM,JM,KC,2))*DYFI(JM)      ! (i-1, j-1, k)
                        DVDY8 = (U_F0_tg(IC,JC,KC,2)-U_F0_tg(IC,JM,KC,2))*DYFI(JM)      ! (i,   j-1, k)
                        
                        dQ_tg(IC,JC,KC,2,2) = ( (DVDY1+DVDY2)*0.5_wp + (DVDY3+DVDY4)*0.5_wp )*0.5_wp*YCL2ND_WFF(JC) + &
                                              ( (DVDY3+DVDY4)*0.5_wp + (DVDY5+DVDY6)*0.5_wp )*0.5_wp*YCL2ND_WFB(JC) 
    
                        !===dv/dx at i',j',k'============================
                        DVDX1 = (U_F0_tg(IC,JC,KM,2)-U_F0_tg(IM,JC,KM,2))*DXI  ! (i',j',k-1)
                        DVDX2 = (U_F0_tg(IC,JC,KC,2)-U_F0_tg(IM,JC,KC,2))*DXI  ! (i',j',k)
                        
                        dQ_tg(IC,JC,KC,2,1) = 0.5_wp*(DVDX1+DVDX2)                 ! (i',j',k')
                        
                        !===dv/dz at i',j',k'============================
                        DVDZ1 = (U_F0_tg(IM,JC,KC,2)-U_F0_tg(IM,JC,KM,2))*DZI       ! (i-1,j',k')
                        DVDZ2 = (U_F0_tg(IC,JC,KC,2)-U_F0_tg(IC,JC,KM,2))*DZI       ! (i,  j',k')
                        
                        dQ_tg(IC,JC,KC,2,3) = (DVDZ1+DVDZ2)*0.5_WP
                        
                        !===============================================================================
                        !===dw/dz at i', j', k'============================
                        DWDZ1 = (U_F0_tg(IM,JC,KC,3)-U_F0_tg(IM,JC,KM,3))*DXI       ! (i-1, j, k-1)
                        DWDZ2 = (U_F0_tg(IC,JC,KC,3)-U_F0_tg(IC,JC,KM,3))*DXI       ! (i,   j, k-1)
                        
                        DWDZ3 = (U_F0_tg(IM,JC,KP,3)-U_F0_tg(IM,JC,KC,3))*DXI       ! (i-1, j, k)
                        DWDZ4 = (U_F0_tg(IC,JC,KP,3)-U_F0_tg(IC,JC,KC,3))*DXI       ! (i,   j, k)
                        
                        DWDZ5 = (U_F0_tg(IM,JM,KC,3)-U_F0_tg(IM,JM,KM,3))*DXI       ! (i-1, j-1, k-1)
                        DWDZ6 = (U_F0_tg(IC,JM,KC,3)-U_F0_tg(IC,JM,KM,3))*DXI       ! (i,   j-1, k-1)
                        
                        DWDZ7 = (U_F0_tg(IM,JM,KP,3)-U_F0_tg(IM,JM,KC,3))*DXI       ! (i-1, j-1, k)
                        DWDZ8 = (U_F0_tg(IC,JM,KP,3)-U_F0_tg(IC,JM,KC,3))*DXI       ! (i,   j-1, k)
                        
                        dQ_tg(IC,JC,KC,3,3) = ( (DWDZ1+DWDZ2)*0.5_wp + (DWDZ3+DWDZ4)*0.5_wp )*0.5_wp*YCL2ND_WFF(JC) + &
                                           ( (DWDZ3+DWDZ4)*0.5_wp + (DWDZ5+DWDZ6)*0.5_wp )*0.5_wp*YCL2ND_WFB(JC) 
    
                        !===dw/dy at i',j',k'============================
                        DWDY1 = (U_F0_tg(IM,JC,KC,3)-U_F0_tg(IM,JM,KC,3))*DYCI(JC)  ! (i-1,j',k')
                        DWDY2 = (U_F0_tg(IC,JC,KC,3)-U_F0_tg(IC,JM,KC,3))*DYCI(JC)  ! (i,  j',k)
                        
                        dQ_tg(IC,JC,KC,3,2) = 0.5_wp*(DWDY1+DWDY2)                 ! (i',j',k')
                        
                        !===dw/dx at i',j',k'============================
                        DWDX1 = (U_F0_tg(IC,JM,KC,3)-U_F0_tg(IM,JM,KC,3))*DZI       ! (i',j-1,k')
                        DWDX2 = (U_F0_tg(IC,JC,KC,3)-U_F0_tg(IM,JC,KC,3))*DZI       ! (i',j,  k')
                        
                        dQ_tg(IC,JC,KC,3,1) = DWDX1*YCL2ND_WFB(JC)+DWDX2*YCL2ND_WFF(JC)
                    END DO
                END DO
            END DO
            
            dQ_tg(:,NND2,:,:,:)=2.0_wp*dQ_tg(:,NCL2,:,:,:)-dQ_tg(:,NCL2-1,:,:,:)
            dQ_tg(:,1,   :,:,:)=2.0_wp*dQ_tg(:,2,   :,:,:)-dQ_tg(:,3,     :,:,:)
            dQ_tg(NND1_tg,:,:,   :,:)=dQ_tg(1,:,:,:,:)
            dQ_tg(:,      :,NND3,:,:)=dQ_tg(:,:,1,:,:)
            
            
            DO IC=1,NND1_tg
                DO JC=1,NND2
                    DO KC=1,NND3
                        !===================vorticity===============================
                        VOR_tg(IC,JC,KC,1) = dQ_tg(IC,JC,KC,3,2) - dQ_tg(IC,JC,KC,2,3) ! VOR_tg-x 
                        VOR_tg(IC,JC,KC,2) = dQ_tg(IC,JC,KC,1,3) - dQ_tg(IC,JC,KC,3,1) ! VOR_tg-y
                        VOR_tg(IC,JC,KC,3) = dQ_tg(IC,JC,KC,2,1) - dQ_tg(IC,JC,KC,1,2) ! VOR_tg-z
                        VOR_tg(IC,JC,KC,4) = DSQRT( DOT_PRODUCT( VOR_tg(IC,JC,KC,1:3),VOR_tg(IC,JC,KC,1:3) ) )!VOR_tg-mag
                        
                        !==========================================================================================
                        !===================omega = antisymetric components of graduate U===============
                        strASYM(1,1) = dQ_tg(IC,JC,KC,1,1) - dQ_tg(IC,JC,KC,1,1)  !dU/dX - dU/dX
                        strASYM(1,2) = dQ_tg(IC,JC,KC,2,1) - dQ_tg(IC,JC,KC,1,2)  !dV/dX - dU/dY
                        strASYM(1,3) = dQ_tg(IC,JC,KC,3,1) - dQ_tg(IC,JC,KC,1,3)  !dW/dX - dU/dZ
                    
                        strASYM(2,1) = dQ_tg(IC,JC,KC,1,2) - dQ_tg(IC,JC,KC,2,1)  !dU/dY - dV/dX
                        strASYM(2,2) = dQ_tg(IC,JC,KC,2,2) - dQ_tg(IC,JC,KC,2,2)  !dV/dY - dV/dY
                        strASYM(2,3) = dQ_tg(IC,JC,KC,3,2) - dQ_tg(IC,JC,KC,2,3)  !dW/dY - dV/dZ
                    
                        strASYM(3,1) = dQ_tg(IC,JC,KC,1,3) - dQ_tg(IC,JC,KC,3,1)  !dU/dZ - dW/dX
                        strASYM(3,2) = dQ_tg(IC,JC,KC,2,3) - dQ_tg(IC,JC,KC,3,2)  !dV/dZ - dW/dY
                        strASYM(3,3) = dQ_tg(IC,JC,KC,3,3) - dQ_tg(IC,JC,KC,3,3)  !dW/dZ - dW/dZ
                        !===================omega = symetric components of graduate U===============
                        strSSYM(1,1) = dQ_tg(IC,JC,KC,1,1) + dQ_tg(IC,JC,KC,1,1)  !dU/dX + dU/dX
                        strSSYM(1,2) = dQ_tg(IC,JC,KC,2,1) + dQ_tg(IC,JC,KC,1,2)  !dV/dX + dU/dY
                        strSSYM(1,3) = dQ_tg(IC,JC,KC,3,1) + dQ_tg(IC,JC,KC,1,3)  !dW/dX + dU/dZ
                    
                        strSSYM(2,1) = dQ_tg(IC,JC,KC,1,2) + dQ_tg(IC,JC,KC,2,1)  !dU/dY + dV/dX
                        strSSYM(2,2) = dQ_tg(IC,JC,KC,2,2) + dQ_tg(IC,JC,KC,2,2)  !dV/dY + dV/dY
                        strSSYM(2,3) = dQ_tg(IC,JC,KC,3,2) + dQ_tg(IC,JC,KC,2,3)  !dW/dY + dV/dZ
                    
                        strSSYM(3,1) = dQ_tg(IC,JC,KC,1,3) + dQ_tg(IC,JC,KC,3,1)  !dU/dZ + dW/dX
                        strSSYM(3,2) = dQ_tg(IC,JC,KC,2,3) + dQ_tg(IC,JC,KC,3,2)  !dV/dZ + dW/dY
                        strSSYM(3,3) = dQ_tg(IC,JC,KC,3,3) + dQ_tg(IC,JC,KC,3,3)  !dW/dZ + dW/dZ
                        
                        strASYM  = strASYM  * 0.5_wp  ! Sij
                        strSSYM  = strSSYM  * 0.5_wp  ! wij
                        
                        ! trace of (STSYM*STSYM^T)
                        STSSYM   =  DOT_PRODUCT(strSSYM(:,1), strSSYM(:,1)) + &
                                    DOT_PRODUCT(strSSYM(:,2), strSSYM(:,2)) + &
                                    DOT_PRODUCT(strSSYM(:,3), strSSYM(:,3))    ! Sij*Sij
                        ! trace of (STASYM*STASYM^T)     
                        STASYM   = DOT_PRODUCT(strASYM(:,1), strASYM(:,1)) + &
                                   DOT_PRODUCT(strASYM(:,2), strASYM(:,2)) + &
                                   DOT_PRODUCT(strASYM(:,3), strASYM(:,3))    ! wij*wij
                                   
                        !================================Q cretiria===================================================
                        QCR_tg(IC,JC,KC)=0.5_wp*(STASYM-STSSYM)
                        
                        
                        !================================Q cretiria===================================================
                        QCR_io(IC,JC,KC)=0.5_wp*(STASYM-STSSYM) ! Q
                        
                        
                        !================================delta_io criteria===============================================
                        ! grad U 
                        gradU   = strASYM + strSSYM
                        
                        !  negtive trace of velocity gradient  !P
                        traceGradU = (strASYM(1,1)+strASYM(2,2)+strASYM(3,3)+&
                                      strSSYM(1,1)+strSSYM(2,2)+strSSYM(3,3))*(-1.0_WP)
                        
                        ! det(gradU) !R
                        detGradU =(gradU(1,1) * ( gradU(2,2)*gradU(3,3) - gradU(2,3)*gradU(3,2) ) - &
                                   gradU(1,2) * ( gradU(2,1)*gradU(3,3) - gradU(3,1)*gradU(2,3) ) + &
                                   gradU(1,3) * ( gradU(2,1)*gradU(3,2) - gradU(2,2)*gradU(3,1) ))*(-1.0_WP)
                        ! delta for incompressible flow
                        !delta_io(IC,JC,KC) = (QCR_io(IC,JC,KC)/3.0_wp)**3 + (detGradU/2.0_wp)**2
                        
                        ! delta for compressible flow
                        delta_tg(IC,JC,KC) = 27.0_WP*detGradU*detGradU + &
                                (4.0_wp*traceGradU*traceGradU*traceGradU-18.0_WP*traceGradU*QCR_tg(IC,JC,KC))*detGradU + &
                                 4.0_wp*QCR_tg(IC,JC,KC)*QCR_tg(IC,JC,KC)*QCR_tg(IC,JC,KC)- &
                                traceGradU*traceGradU*QCR_tg(IC,JC,KC)*QCR_tg(IC,JC,KC)             
                        
                        !get rid of possible infinity 
                        if (delta_tg(IC,JC,KC)> 1.0E30_wp) delta_tg(IC,JC,KC)= 1.0E30_wp
                        if (delta_tg(IC,JC,KC)<-1.0E30_wp) delta_tg(IC,JC,KC)=-1.0E30_wp
             
                        !==============swirl strength ==================================
                        if (delta_tg(IC,JC,KC)>0.0_WP) then
                        
                            matz = 0
                            nm   = 3
                            nn   = 3
                            ar(:,:)   = strASYM(:,:) + strSSYM(:,:)
                            ai(:,:)   = 0.0
                            
                            call cg(nm,nn,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
                            
                            if(dabs(wi(1))<(dabs(wi(2))+1.0E-14) .and. dabs(wi(1))<(dabs(wi(3))+1.0E-14) ) then
                                swirlstrength_tg(IC,JC,KC,1) = wr(1)
                                swirlstrength_tg(IC,JC,KC,2) = wr(2)
                                swirlstrength_tg(IC,JC,KC,3) = wr(2)
                            else if (dabs(wi(2))<(dabs(wi(1))+1.0E-14) .and. dabs(wi(2))<(dabs(wi(3))+1.0E-14) ) then
                                swirlstrength_tg(IC,JC,KC,1) = wr(2)
                                swirlstrength_tg(IC,JC,KC,2) = wr(1)
                                swirlstrength_tg(IC,JC,KC,3) = wr(1)
                            else if (dabs(wi(3))<(dabs(wi(1))+1.0E-14) .and. dabs(wi(3))<(dabs(wi(2))+1.0E-14) ) then
                                swirlstrength_tg(IC,JC,KC,1) = wr(3)
                                swirlstrength_tg(IC,JC,KC,2) = wr(1)
                                swirlstrength_tg(IC,JC,KC,3) = wr(1)
                            else
                            end if
 
                        end if
             
             
             
                        !*************lamda2 criteria *************************
                        S2plusO2(1,1) = strSSYM(1,1)*strSSYM(1,1) + strSSYM(1,2)*strSSYM(2,1) + strSSYM(1,3)*strSSYM(3,1) + &
                                        strASYM(1,1)*strASYM(1,1) + strASYM(1,2)*strASYM(2,1) + strASYM(1,3)*strASYM(3,1)
                        S2plusO2(1,2) = strSSYM(1,1)*strSSYM(1,2) + strSSYM(1,2)*strSSYM(2,2) + strSSYM(1,3)*strSSYM(3,2) + &
                                        strASYM(1,1)*strASYM(1,2) + strASYM(1,2)*strASYM(2,2) + strASYM(1,3)*strASYM(3,2)
                        S2plusO2(1,3) = strSSYM(1,1)*strSSYM(1,3) + strSSYM(1,2)*strSSYM(2,3) + strSSYM(1,3)*strSSYM(3,3) + &
                                        strASYM(1,1)*strASYM(1,3) + strASYM(1,2)*strASYM(2,3) + strASYM(1,3)*strASYM(3,3)
                             
                        S2plusO2(2,1) = strSSYM(2,1)*strSSYM(1,1) + strSSYM(2,2)*strSSYM(2,1) + strSSYM(2,3)*strSSYM(3,1) + &
                                        strASYM(2,1)*strASYM(1,1) + strASYM(2,2)*strASYM(2,1) + strASYM(2,3)*strASYM(3,1)                 
                        S2plusO2(2,2) = strSSYM(2,1)*strSSYM(1,2) + strSSYM(2,2)*strSSYM(2,2) + strSSYM(2,3)*strSSYM(3,2) + &
                                        strASYM(2,1)*strASYM(1,2) + strASYM(2,2)*strASYM(2,2) + strASYM(2,3)*strASYM(3,2)           
                        S2plusO2(2,3) = strSSYM(2,1)*strSSYM(1,3) + strSSYM(2,2)*strSSYM(2,3) + strSSYM(2,3)*strSSYM(3,3) + &
                                        strASYM(2,1)*strASYM(1,3) + strASYM(2,2)*strASYM(2,3) + strASYM(2,3)*strASYM(3,3)
    
                        S2plusO2(3,1) = strSSYM(3,1)*strSSYM(1,1) + strSSYM(3,2)*strSSYM(2,1) + strSSYM(3,3)*strSSYM(3,1) + &
                                        strASYM(3,1)*strASYM(1,1) + strASYM(3,2)*strASYM(2,1) + strASYM(3,3)*strASYM(3,1)                        
                        S2plusO2(3,2) = strSSYM(3,1)*strSSYM(1,2) + strSSYM(3,2)*strSSYM(2,2) + strSSYM(3,3)*strSSYM(3,2) + &
                                        strASYM(3,1)*strASYM(1,2) + strASYM(3,2)*strASYM(2,2) + strASYM(3,3)*strASYM(3,2)                       
                        S2plusO2(3,3) = strSSYM(3,1)*strSSYM(1,3) + strSSYM(3,2)*strSSYM(2,3) + strSSYM(3,3)*strSSYM(3,3) + &
                                        strASYM(3,1)*strASYM(1,3) + strASYM(3,2)*strASYM(2,3) + strASYM(3,3)*strASYM(3,3)
    
                        call DSYEVC3(S2plusO2,eig)
                        if (((eig(1)-eig(2))*(eig(1)-eig(3)))<0.0_wp) lambda2_tg(IC,JC,KC)=eig(1)
                        if (((eig(2)-eig(1))*(eig(2)-eig(3)))<0.0_wp) lambda2_tg(IC,JC,KC)=eig(2)
                        if (((eig(3)-eig(1))*(eig(3)-eig(2)))<0.0_wp) lambda2_tg(IC,JC,KC)=eig(3)
    
               
                        !get rid of possible infinity
                        if (lambda2_tg(IC,JC,KC)> 1.0E30_wp) lambda2_tg(IC,JC,KC)= 1.0E30_wp
                        if (lambda2_tg(IC,JC,KC)<-1.0E30_wp) lambda2_tg(IC,JC,KC)=-1.0E30_wp
                        !***********************************************************************************
                    END DO
                END DO
            END DO
            deallocate(dQ_tg)
        END IF
        
    
    
        RETURN
    END SUBROUTINE 
     
!**********************************************************************************************************************************
!    SUBROUTINE TEC360_PROFILE_TEST
!        USE QPGATHERING_INFO
!        use TEC360_INFO
!        use mesh_info
!        use flow_info
!        use init_info
!        use thermal_info
!        IMPLICIT NONE
      
!        INTEGER(4) :: I, J, K
!        INTEGER(4)  :: TECFLG = 202
      
!        LOGICAL :: file_exists
        
!        INQUIRE(FILE='RESULT.Centerline.tec', EXIST=file_exists) 
           
!        IF(file_exists) then
!            OPEN(TECFLG,FILE='RESULT.Centerline.tec', position='append')
!        else
!            OPEN(TECFLG,FILE='RESULT.Centerline.tec')
!        end if
        
!        WRITE(TECFLG,'(A)') 'TITLE = "DNS Centerline"'
!        WRITE(TECFLG,'(A)') 'VARIABLES = "X", "U", "V", "W", "P"'
        
!        WRITE(TECFLG,'(A,1I6.1,1ES13.5,A)') 'ZONE T=" ',ITERG,phyTIME, '"'
       
!        J=NCL2/2
!        K=NCL3/2
        
!        DO I=0,NCL1_io
!           write(TECFLG,'(2ES15.7)') 0.5_wp*(XND_io(I)+XND_io(I+1)), P_F0_io(I,J,K)
!        END DO
!        close(tecflg)
        
!        RETURN
!     END SUBROUTINE 
     
!    !***********************************************************************************************
!    SUBROUTINE BULK_PARAMETER_CHECK
!    use TEC360_INFO
!    use mesh_info
!    use flow_info
!    use init_info
!    IMPLICIT NONE
    
!    INTEGER(4) :: IM, IC, JC, KC
!    INTEGER(4)  :: TECFLG = 202
    
!    REAL(WP)   :: MASS_RATE(NCL1_IO)
!    REAL(WP)   :: MASS_FLUX(NCL1_IO)
    
!    REAL(WP)   :: ENTH_RATE
!    REAL(WP)   :: ENTH_FLUX(NCL1_IO)
    
!    REAL(WP)   :: RHOU, RHO, RHOU_Z
!    REAL(WP)   :: RHOH, RHOHU,RHOHU_Z
    
    
!    !=============BULK MASS FLUX========================
!    DO  IC=1,NCL1_IO
!        IM = IMV_IO(IC)
       
!        MASS_RATE(IC) = 0.0_WP
       
!        DO JC=1,NCL2
!            !====AVERAGE IN HOMOGENOUS DIRECTION
!            RHOU = 0.0_WP
!            DO KC=1,NCL3
!                RHO  = D_F0_io(IM,JC,KC) + D_F0_io(IC,JC,KC)
!                RHOU = RHOU+U_F0_io(IC,JC,KC)*RHO
!            END DO
!            RHOU_Z =  RHOU/NCL3/2.0_WP
          
!            MASS_RATE(IC) = MASS_RATE(IC) + RHOU_Z*DYFI(JC)/DZI
!        END DO
!        MASS_FLUX(IC) =  MASS_RATE(IC)/AREA_INLET
        
!    END DO
    
    
!    !=============BULK ENTHALPY========================
!    DO  IC=1,NCL1_IO
!        IM = IMV_IO(IC)
       
!        ENTH_RATE = 0.0_WP
       
!        DO JC=1,NCL2
!            !====AVERAGE IN HOMOGENOUS DIRECTION
!            RHOHU = 0.0_WP
!            DO KC=1,NCL3
!                RHOH  = D_F0_io(IM,JC,KC)*H_F0_io(IM,JC,KC) + &
!                        D_F0_io(IC,JC,KC)*H_F0_io(IC,JC,KC)
                        
!                RHOHU = RHOHU+U_F0_io(IC,JC,KC)*RHOH
!            END DO
!            RHOHU_Z =  RHOHU/NCL3/2.0_WP
          
!            ENTH_RATE = ENTH_RATE + RHOHU_Z*DYFI(JC)/DZI
!        END DO
!        ENTH_FLUX(IC) =  ENTH_RATE/MASS_RATE(IC)
!    END DO
    
!    !===============WRITE OUT================================
    
!    IF(NCOUNT(4)==0) THEN
!       OPEN(TECFLG,FILE='CHECK_CONSERVATION.tec')
!    ELSE
!       OPEN(TECFLG,FILE='CHECK_CONSERVATION.tec', position='append')
!    END IF
!    NCOUNT(4) = NCOUNT(4)+1
!    WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW CONSERVATION"'
!    WRITE(TECFLG,'(A)') 'VARIABLES = "X", "Gb", "Hb"'
!    WRITE(TECFLG,'(A,1I6.1,1ES13.5)') 'ZONE T=" ',ITERG,phyTIME
    
!    DO IC=1,NCL1
!       WRITE(*,*) XND_io(IC),MASS_FLUX(IC),ENTH_FLUX(IC)
!    END DO 
    
!    RETURN
!    END SUBROUTINE
     
     
     
              
         
        !===============CONSERVATION CHECK=============================================================================
!         !=============BULK MASS FLUX========================
!        DO  IC=1,NCL1_IO
!            IM = IMV_IO(IC)
           
!            MASS_RATE(IC) = 0.0_WP
           
!            DO JC=1,NCL2
!                !====AVERAGE IN HOMOGENOUS DIRECTION
!                RHOUt = 0.0_WP
!                DO KC=1,NCL3
!                    RHOt  = D_F0_io(IM,JC,KC) + D_F0_io(IC,JC,KC)
!                    RHOUt = RHOUt+U_F0_io(IC,JC,KC)*RHOt
!                END DO
!                RHOU_Z =  RHOUt/NCL3/2.0_WP
              
!                MASS_RATE(IC) = MASS_RATE(IC) + RHOU_Z/DYFI(JC)/DZI
!            END DO
!            MASS_FLUX(IC) =  MASS_RATE(IC)/AREA_INLET
            
!        END DO
        
        
!        !=============BULK ENTHALPY========================
!        DO  IC=1,NCL1_IO
!            IM = IMV_IO(IC)
           
!            ENTH_RATE = 0.0_WP
           
!            DO JC=1,NCL2
!                !====AVERAGE IN HOMOGENOUS DIRECTION
!                RHOHUt = 0.0_WP
!                DO KC=1,NCL3
!                    RHOHt  = D_F0_io(IM,JC,KC)*H_F0_io(IM,JC,KC) + &
!                             D_F0_io(IC,JC,KC)*H_F0_io(IC,JC,KC)
                            
!                    RHOHUt = RHOHUt+U_F0_io(IC,JC,KC)*RHOHt
!                END DO
!                RHOHU_Z =  RHOHUt/NCL3/2.0_WP
              
!                ENTH_RATE = ENTH_RATE + RHOHU_Z/DYFI(JC)/DZI
!            END DO
!            ENTH_FLUX(IC) =  ENTH_RATE/MASS_RATE(IC)
!        END DO
        
!        !===============WRITE OUT================================
        
!        IF(NCOUNT(5).EQ.0) THEN
!           OPEN(TECFLG,FILE='CHECK_CONSERVATION.tec')
!        ELSE
!           OPEN(TECFLG,FILE='CHECK_CONSERVATION.tec', position='append')
!        END IF
!        NCOUNT(5) = NCOUNT(5)+1
!        WRITE(TECFLG,'(A)') 'TITLE = "DNS FLOW CONSERVATION"'
!        WRITE(TECFLG,'(A)') 'VARIABLES = "X", "Gb", "Hb"'
!        WRITE(TECFLG,'(A,1I6.1,1ES13.5)') 'ZONE T=" ',ITERG,phyTIME
        
!        DO IC=1,NCL1_io
!           WRITE(TECFLG,'(3ES15.7)') XND_io(IC),MASS_FLUX(IC),ENTH_FLUX(IC)
!        END DO 
        
        
        !<=================================================================
     !===============================TG=========================

!=========================================================================================
!    SUBROUTINE QPGATHERING_TG 
!!>    @Warning: 1) MPI_GATHER requires the data account from each partition
!!>             sould be the same.  
!        USE QPGATHERING_INFO
!        USE MESH_INFO
!        USE FLOW_INFO
!        USE INIT_INFO
!        use postprocess_info
!        IMPLICIT NONE
        
!        INTEGER(4)  :: INN1, INN2, INN3
!        INTEGER(4)  :: I,L
!        INTEGER(4)  :: J, JJ, JJP, JP
!        INTEGER(4)  :: K, KK, KS
!        INTEGER(4)  :: N2DOID
        
        
!        IF( .NOT.TGFLOWflg) RETURN
        
!        !==================ALLOCATTE DATA=============================
!        ALLOCATE ( U_F0_tg(NCL1_tg,0:NND2,  NCL3) )
!        ALLOCATE ( U_F0_tg(NCL1_tg,0:NND2,NCL3) )
!        ALLOCATE ( U_F0_tg(NCL1_tg,0:NND2,  NCL3) )
!        ALLOCATE ( P_F0_tg(NCL1_tg,NCL2,  NCL3) )
        
!        ALLOCATE ( QAUX_tg(NCL1_tg,0:N2DO(0)+1,NCL3,NDV,1:NPTOT) )
!        ALLOCATE ( PAUX_tg(NCL1_tg,0:N2DO(0)+1,NCL3,    1:NPTOT) )
        
!        ALLOCATE ( D1AUX (N2DO(MYID), NDV+1,1:NPTOT) )
!        ALLOCATE ( U1xztL_F0_tg( NCL2, NDV+1 ) )
        
        
            
!        !=============gather data to master===========================
!        INN1 = NCL1_tg*(N2DO(0)+2)*NCL3*NDV
!        CALL MPI_GATHER( Q_tg, INN1, MPI_DOUBLE_PRECISION, QAUX_tg, INN1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)                
!        CALL MPI_BARRIER(ICOMM,IERROR)

!        INN2 = NCL1_tg*(N2DO(0)+2)*NCL3
!        CALL MPI_GATHER( PR_tg, INN2, MPI_DOUBLE_PRECISION, PAUX_tg, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
!        CALL MPI_BARRIER(ICOMM,IERROR)
        
!        INN3 = N2DO(MYID)*(NDV+1)                            
!        CALL MPI_GATHER( U1xztL_tg, INN3, MPI_DOUBLE_PRECISION, D1AUX, INN3,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR) 
!        CALL MPI_BARRIER(ICOMM,IERROR)

         
!        !==============re-arrange gathered data======================== 
!        IF (MYID.EQ.0) THEN
        
!            DO KK = 0, NPSLV
!                N2DOID=JDEWT(KK)-JDSWT(KK)+1
!                DO J=1,N2DOID
!                    JJ=JDSWT(KK)-1+J
!                    DO L=1,NDV+1
!                        U1xztL_F0_tg(JJ,L)=D1AUX(J,L,KK+1)
!                    ENDDO
!                ENDDO
!            ENDDO
      
!            DO KK=0,NPSLV
!                N2DOID=JDEWT(KK)-JDSWT(KK)+1
!                DO  J=1,N2DOID
!                    JJ=JDSWT(KK)-1+J
!                    JP = J + 1
!                    JJP=JJ+1
!                    DO I=1,NCL1_tg
!                        DO K=1,NCL3
!                            KS = KSYM(K)
                            
!                            U_F0_tg(I,JJ,K)=(QAUX_tg(I,J,K,1,KK+1))
!                            U_F0_tg(I,JJ,K)=(QAUX_tg(I,J,K,3,KK+1)*RCCI1(JJ))
!                            IF(JJ==1 .and. icase==IPIPEC) THEN
!                                U_F0_tg(I,JJ,K)=(QAUX_tg(I,JP,K,2,KK+1) - QAUX_tg(I,JP,KS,2,KK+1) ) * 0.50_WP*RNDI1(JJP)
!                            ELSE
!                                U_F0_tg(I,JJ,K)=(QAUX_tg(I,J,K,2,KK+1)*RNDI1(JJ))
!                            END IF
!                            P_F0_tg(I,JJ,K)=(PAUX_tg(I,J,K,KK+1))
!                        ENDDO
!                    ENDDO
!                ENDDO
!            END DO
   
!            DO I=1,NCL1_tg
!                DO K=1,NCL3
!                    U_F0_tg(I,NND2,K)=(QAUX_tg(I,N2DO(MYID)+1,K,2,NPTOT))  
!                ENDDO
!            ENDDO
            
!            DO I=1,NCL1_tg
!                DO K=1,NCL3
!                    U_F0_tg(I,0,   K)= U_F0_tg(I,1,   K)*(-1.0_wp)
!                    U_F0_tg(I,NND2,K)= U_F0_tg(I,NND2,K)*(-1.0_wp)
!                    U_F0_tg(I,0,   K)= U_F0_tg(I,1,   K)*(-1.0_wp)
!                    U_F0_tg(I,NND2,K)= U_F0_tg(I,NND2,K)*(-1.0_wp)
!                    U_F0_tg(I,0,   K)= 0.0_wp
!                ENDDO
!            ENDDO
    
!        END IF
        
!        RETURN
!    END SUBROUTINE
    
!!=========================IO===============================================    
!    SUBROUTINE QPGATHERING_IO 
!        USE MESH_INFO
!        USE FLOW_INFO
!        USE INIT_INFO
!        USE THERMAL_INFO
!        USE QPGATHERING_INFO
!        use postprocess_info
!        IMPLICIT NONE
        
!        INTEGER(4)  :: INN1, INN2, INN3
!        INTEGER(4)  :: I, L
!        INTEGER(4)  :: J, JJ, JJP, JP
!        INTEGER(4)  :: K, KK, KS
!        INTEGER(4)  :: N2DOID
        
        
        
!        IF( .NOT.IOFLOWflg) RETURN
        
!        !==================ALLOCATTE DATA=============================
        
!        ALLOCATE ( U_F0_io(NCL1S:NCL1E, 0:NND2,  NCL3) )
!        ALLOCATE ( U_F0_io(NCL1S:NCL1E, 0:NND2,  NCL3) )
!        ALLOCATE ( U_F0_io(NCL1S:NCL1E, 0:NND2,  NCL3) )
!        ALLOCATE ( QAUX_io(NCL1S:NCL1E, 0:N2DO(0)+1,NCL3,NDV,1:NPTOT) )
!        ALLOCATE ( PAUX_io(NCL1S:NCL1E, 0:N2DO(0)+1,NCL3,    1:NPTOT) )
!        ALLOCATE ( P_F0_io(NCL1S:NCL1E, NCL2,  NCL3) )
        
!        IF(thermlflg==1) THEN
!            ALLOCATE ( T_F0_io(NCL1S:NCL1E, NCL2,  NCL3) )
!            ALLOCATE ( D_F0_io(NCL1S:NCL1E, NCL2,  NCL3) )
!            ALLOCATE ( H_F0_io(NCL1S:NCL1E, NCL2,  NCL3) )
!            ALLOCATE ( K_F0_io(NCL1S:NCL1E, NCL2,  NCL3) )
!            ALLOCATE ( M_F0_io(NCL1S:NCL1E, NCL2,  NCL3) )
!            ALLOCATE ( TAUX_io(NCL1S:NCL1E, 0:N2DO(0)+1,NCL3,    1:NPTOT) )
!            ALLOCATE ( DAUX_io(NCL1S:NCL1E, 0:N2DO(0)+1,NCL3,    1:NPTOT) )
!            ALLOCATE ( HAUX_io(NCL1S:NCL1E, 0:N2DO(0)+1,NCL3,    1:NPTOT) )
!            ALLOCATE ( KAUX_io(NCL1S:NCL1E, 0:N2DO(0)+1,NCL3,    1:NPTOT) )
!            ALLOCATE ( MAUX_io(NCL1S:NCL1E, 0:N2DO(0)+1,NCL3,    1:NPTOT) )
!        END IF
        
!        ALLOCATE ( D1AUXX (NCL1_io,N2DO(MYID), NDV+1,1:NPTOT) )
!        ALLOCATE ( U1ztL_F0_io(NCL1_io,NCL2, NDV+1 ) )
        
!        ALLOCATE ( D1AUXZ (N2DO(MYID), NDV+1,1:NPTOT) )
!        ALLOCATE ( U1xztL_F0_io(NCL2, NDV+1 ) )
        
!        !if(myid==0) then
!          !write(*,*) 'U,myid',Q_io(ncl1_tg/2,:,ncl3/2,1)
!        !end if
        
!        !=============gather data to master===========================
!        INN1=(NCL1E-NCL1S+1)*(N2DO(0)+2)*NCL3*NDV
!        CALL MPI_GATHER( Q_io, INN1, MPI_DOUBLE_PRECISION, QAUX_io, INN1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR) 
!        CALL MPI_BARRIER(ICOMM,IERROR)
    
!        INN2=(NCL1E-NCL1S+1)*(N2DO(0)+2)*NCL3
!        CALL MPI_GATHER( PR_io, INN2, MPI_DOUBLE_PRECISION, PAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)          
!        CALL MPI_BARRIER(ICOMM,IERROR)
        
!        IF(thermlflg==1) THEN
!            CALL MPI_GATHER( TEMPERATURE, INN2, MPI_DOUBLE_PRECISION, TAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)            
!            CALL MPI_BARRIER(ICOMM,IERROR)
             
!            CALL MPI_GATHER( DENSITY, INN2, MPI_DOUBLE_PRECISION, DAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)           
!            CALL MPI_BARRIER(ICOMM,IERROR)
             
!            CALL MPI_GATHER( ENTHALPY, INN2, MPI_DOUBLE_PRECISION, HAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)             
!            CALL MPI_BARRIER(ICOMM,IERROR)
            
!            CALL MPI_GATHER( VISCOUSITY, INN2, MPI_DOUBLE_PRECISION, MAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)    
!            CALL MPI_BARRIER(ICOMM,IERROR)
             
!            CALL MPI_GATHER( THERMCONDT, INN2, MPI_DOUBLE_PRECISION, KAUX_io, INN2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)               
!            CALL MPI_BARRIER(ICOMM,IERROR)
!        END IF
        
!        IF(TGFLOWflg) THEN
!            INN3  =NCL1_io* N2DO(MYID)*(NDV+1)                            
!            CALL MPI_GATHER( U1ztL_io, INN3, MPI_DOUBLE_PRECISION, D1AUXX, INN3,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR) 
!            CALL MPI_BARRIER(ICOMM,IERROR)
!        ELSE
!            INN3 = N2DO(MYID)*(NDV+1)                            
!            CALL MPI_GATHER( U1xztL_io, INN3, MPI_DOUBLE_PRECISION, D1AUXZ, INN3,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR) 
!            CALL MPI_BARRIER(ICOMM,IERROR)
!        END IF
        
!        !==============re-arrange gathered data============================= 
!        IF (MYID.EQ.0) THEN
        
!            IF(TGFLOWflg) THEN
!                DO KK = 0, NPSLV
!                    N2DOID=JDEWT(KK)-JDSWT(KK)+1
!                    DO I=1,NCL1_io
!                        DO J=1,N2DOID
!                            JJ=JDSWT(KK)-1+J
!                            DO L=1,NDV+1
!                                !U1xtL_F0_io(I,JJ,L)=D1AUXX(I,J,L,KK+1)
!                            END DO
!                        END DO
!                    END DO
!                END DO
!            ELSE
!                DO KK = 0, NPSLV
!                    N2DOID=JDEWT(KK)-JDSWT(KK)+1
!                    DO J=1,N2DOID
!                        JJ=JDSWT(KK)-1+J
!                        DO L=1,NDV+1
!                            U1xztL_F0_io(JJ,L)=D1AUXZ(J,L,KK+1)
!                        END DO
!                    END DO
!                END DO
!            END IF
            
            

!            DO KK=0,NPSLV
!                N2DOID=JDEWT(KK)-JDSWT(KK)+1
!                DO J=1,N2DOID
!                    JJ=JDSWT(KK)-1+J
!                    JP = J+1
!                    JJP= JJ + 1
!                    DO I=NCL1S, NCL1E
!                        DO K=1,NCL3
!                            KS=KSYM(K)
!                            U_F0_io(I,JJ,K)= QAUX_io(I,J,K,1,KK+1)
!                            U_F0_io(I,JJ,K)= QAUX_io(I,J,K,3,KK+1)*RCCI1(JJ)
!                            IF(JJ==1 .and. icase==IPIPEC) THEN
!                                U_F0_io(I,JJ,K)=(QAUX_IO(I,JP,K,2,KK+1) - QAUX_IO(I,JP,KS,2,KK+1) ) * 0.50_WP*RNDI1(JJP)
!                            ELSE
!                                U_F0_io(I,JJ,K)=(QAUX_IO(I,J,K,2,KK+1)*RNDI1(JJ))
!                            END IF
!                            P_F0_io(I,JJ,K)= PAUX_io(I,J,K,KK+1)
!                            IF(thermlflg==1) THEN
!                                T_F0_io(I,JJ,K)= TAUX_io(I,J,K,KK+1)
!                                D_F0_io(I,JJ,K)= DAUX_io(I,J,K,KK+1)
!                                H_F0_io(I,JJ,K)= HAUX_io(I,J,K,KK+1)
!                                K_F0_io(I,JJ,K)= KAUX_io(I,J,K,KK+1)
!                                M_F0_io(I,JJ,K)= MAUX_io(I,J,K,KK+1)
!                            END IF
!                        ENDDO
!                    ENDDO
!                ENDDO
!            END DO
   
!            DO I=NCL1S, NCL1E
!                DO K=1,NCL3
!                    U_F0_io(I,NND2,K)= QAUX_io(I,N2DO(MYID)+1,K,2,NPTOT) 
!                ENDDO
!            ENDDO 
            
!            DO I=NCL1S, NCL1E
!                DO K=1,NCL3
!                    U_F0_io(I,0,   K)= U_F0_io(I,1,   K)*(-1.0_wp)
!                    U_F0_io(I,NND2,K)= U_F0_io(I,NND2,K)*(-1.0_wp)
!                    U_F0_io(I,0,   K)= U_F0_io(I,1,   K)*(-1.0_wp)
!                    U_F0_io(I,NND2,K)= U_F0_io(I,NND2,K)*(-1.0_wp)
!                    U_F0_io(I,0,   K)= 0.0_wp
!                ENDDO
!            ENDDO
            
!        END IF
        
        
!        RETURN    
!    END SUBROUTINE
