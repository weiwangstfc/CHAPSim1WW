    MODULE INTP_VARS
        USE WPRECISION
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        
        integer(4) :: NCLO1, NCLO2, NCLO3, NCL1
        integer(4) :: NCLOO1, NCLOO3
        
        REAL(WP)   :: HXO, HZO, HXOO, HZOO, HX
        integer(4),ALLOCATABLE :: N2DOO(:)
        integer(4),ALLOCATABLE :: N3DOO(:)
        REAL(WP),ALLOCATABLE :: XNDO(:)
        REAL(WP),ALLOCATABLE :: YNDO(:)
        REAL(WP),ALLOCATABLE :: ZNDO(:)
        REAL(WP),ALLOCATABLE :: XCCO(:)
        REAL(WP),ALLOCATABLE :: YCCO(:)
        REAL(WP),ALLOCATABLE :: ZCCO(:)
        !REAL(WP),ALLOCATABLE :: ZCC(:)
        REAL(WP),ALLOCATABLE :: XND(:)
        REAL(WP),ALLOCATABLE :: XCC(:)
    END MODULE
    
!=========================================================================================================
    SUBROUTINE INITIAL_INTERP_MESH
        USE INTP_VARS
        IMPLICIT NONE
        
        CHARACTER(128) :: sect          ! section names
        CHARACTER(128) :: skey          ! key names
        INTEGER(4) :: IOS=0
        INTEGER(4) :: INI=13             !file I/O id
        INTEGER(4) :: lens   
        real(wp)  :: rtmp
        character(256) :: stmp
        INTEGER(4)  :: I, J, K
        
        
        !*******Read in old y mesh*******************************************************
        if(myid==0) CALL CHKHDL('       Reading in intp y coordinates',myid) 
        
        OPEN(INI,FILE='INTP_COORD_YND.ini', status='old', iostat=ios)
        if(ios .NE. 0)  &
        call ERRHDL(' File INTP_COORD_YND.ini cannot be found.',MYID)
        
        read(ini,*) sect
        lens = len_trim(sect)
        IF (MYID.EQ.0) call CHKHDL('         Read in '//sect(1:lens),MYID)
        if(sect(1:lens) /= '[origeo]')  &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID)
        
        read(ini,*) skey, HXOO
        read(ini,*) skey, HZOO
        read(ini,*) skey, NCLOO1
        read(ini,*) skey, NCLO2
        read(ini,*) skey, NCLOO3
        
        IF(MYID==0) THEN
            IF((NCLO2/NPTOT).LT.1) CALL ERRHDL(' Interpolation from coarse mesh with a CPU number larger than NCL2',myid)
        END IF
        
        !*******allocate variables*************************************************
        
        ALLOCATE (YNDO(NCLO2+1))
        read(ini,*) sect
        lens = len_trim(sect)
        IF (MYID.EQ.0) call CHKHDL('         Read in '//sect(1:lens),MYID)
        if(sect(1:lens) /= '[ynd]')  &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID)
        
        
        DO J=1,NCLO2 + 1
            read(ini,*) lens, YNDO(J)
        END DO
        CLOSE(INI)
        
        !*******repeating (not stretching) to the real domain******************
        IF(TGFLOWflg .AND. IOFLOWflg) THEN
            WRITE(*,*) 'MODIFY CODE, PLEASE'
        ELSE
            IF(TGFLOWflg) THEN
                HX    = HX_TG
                NCL1  = NCL1_TG
            END IF
            IF(IOFLOWflg) THEN 
                HX    = HX_IO
                NCL1  = NCL1_IO
            END IF
            ALLOCATE (XND(NCL1+1))
            ALLOCATE (XCC(NCL1))
            IF(TGFLOWflg) THEN
                XND   = XND_tg
                XCC   = XCC_tg
            END IF
            IF(IOFLOWflg) THEN 
                XND   = XND_io
                XCC   = XCC_io
            END IF
        END IF
        
        NCLO1 = ceiling(HX/(HXOO/DBLE(NCLOO1)))
        NCLO3 = ceiling(HZ/(HZOO/DBLE(NCLOO3)))
        HXO = HX
        HZO = HZ
        
        IF(MYID==0) THEN
            call CHKHDL    ('         Repeating the originl domain to the real domain',MYID)
            call CHK2RLHDL ('            Lx00->Lx0',MYID,HXOO, HXO)
            call CHK2RLHDL ('            Lz00->Lz0',MYID,HZOO, HZO)
            call CHK2INTHDL('            Nx00->Nx0',MYID,NCLOO1, NCLO1)
            call CHK2INTHDL('            Nz00->Nz0',MYID,NCLOO3, NCLO3)
        END IF
        
        
        
        ALLOCATE (N3DOO(0:NPSLV))
        ALLOCATE (N2DOO(0:NPSLV))
        N2DOO(0:NPSLV) = INT(NCLO2/NPTOT)
        N3DOO(0:NPSLV) = INT(NCLO3/NPTOT)
        !*********OLD AND NEW COORDINATES******************************************************************
        ALLOCATE (YCCO(NCLO2) )
        ALLOCATE (XNDO(NCLO1+1))
        ALLOCATE (ZNDO(NCLO3+1))
        ALLOCATE (XCCO(NCLO1))
        ALLOCATE (ZCCO(NCLO3))
        !ALLOCATE (ZCC(NCL3))
        
        DO  I=1,(NCLO1+1)
            XNDO(I) = DBLE(I-1)*HXO/DBLE(NCLO1)
        ENDDO

        DO  K=1,(NCLO3+1)
            ZNDO(K) = DBLE(K-1)*HZO/DBLE(NCLO3)
        ENDDO
        
        DO I=1, NCLO1
            XCCO(I) = 0.5_WP*(XNDO(I)+XNDO(I+1))
        END DO
        
        DO K=1, NCLO3
            ZCCO(K) = 0.5_WP*(ZNDO(K)+ZNDO(K+1))
        END DO

        DO J=1, NCLO2
            YCCO(J) = 0.5_wp*( YNDO(J) + YNDO(J+1) )
        END DO
        
        !DO K=1, NCL3
            !ZCC(K) = 0.5_WP*( ZND(K) + ZND(K+1) )
        !END DO
    
        RETURN
    END SUBROUTINE
    
!=========================================================================================================    
    SUBROUTINE INITIAL_INTERP_tg
        USE INTP_VARS
        IMPLICIT NONE
        
        CHARACTER(LEN=256)  :: WRT_RST_FNM
        CHARACTER(15) :: PNTIM
        INTEGER(4) :: I,J,K,N, KK, II
        real(wp)    :: Umean1_tg
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:)
        REAL(WP),ALLOCATABLE :: QO   (:,:,:)
        REAL(WP),ALLOCATABLE :: QN   (:,:,:)
        
        if(myid==0) CALL CHKHDL('    TG: Interpolating from one mesh to another',myid) 
        
        !=========STEP 1, GENERATE MESH INFO==================
        CALL INITIAL_INTERP_MESH
        
        !=========STEP 2 ALLOCATE VARS===================
        ALLOCATE( DUMMY (NCLOO1,N2DOO(MYID),NCLOO3) )
        ALLOCATE( QO    (NCLO1, N2DOO(MYID),NCLO3)  )
        ALLOCATE( QN    (NCL1,  N2DO(MYID), NCL3 )  )
        
        !==========STEP 3 READ IN OLD DATA AND INTP INTO NEW ONES===============
        WRITE(PNTIM,'(1ES15.9)') RSTtim_tg
        DO N=1,NDV+1
            !==========(1) FILE NAME===========
            IF(N==1)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_U.D'
            IF(N==2)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_V.D'
            IF(N==3)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_W.D'
            IF(N==4)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_P.D'
            
            !==========(2) READ IN DATA===========
            CALL READ_3D_VARS(NCLOO1,NCLO2,NCLOO3,N2DOO(MYID),N2DOO(0)*MYID,ITERG0_TG,phyTIME_TG, DUMMY, WRT_RST_FNM)
            
            !==========(3) EXPAND READ-IN TO REAL LENGTH===========
            DO I=1,NCLO1
                DO J=1,N2DOO(MYID)
                    DO K=1,NCLO3
                    
                        II=MOD(I,NCLOO1)
                        IF(II==0) II=NCLOO1
                        
                        KK=MOD(K,NCLOO3)
                        IF(KK==0) KK=NCLOO3
                        
                        QO(I,J,K)=DUMMY(II,J,KK)
                    ENDDO
                ENDDO
            ENDDO
            
            !==========(4) INTEP DATA===========
            CALL INTP_3D_VARS(QO,QN,N)
            
            DO J=1,N2DO(MYID)
                DO I=1,NCL1
                    DO K=1,NCL3
                        IF(N==1) Q_TG(I,J,K,1)=QN(I,J,K) !U
                        IF(N==2) Q_TG(I,J,K,2)=QN(I,J,K) !V
                        IF(N==3) Q_TG(I,J,K,3)=QN(I,J,K) !W
                        IF(N==4) PR_TG(I,J,K) =QN(I,J,K) !P
                    ENDDO
                ENDDO
            ENDDO

        END DO
        
        DEALLOCATE (DUMMY)
        DEALLOCATE (QO)
        DEALLOCATE (QN)
        DEALLOCATE  (N2DOO)
        DEALLOCATE  (N3DOO)
        DEALLOCATE  (XNDO)
        DEALLOCATE  (YNDO)
        DEALLOCATE  (ZNDO)
        DEALLOCATE  (XCCO)
        DEALLOCATE  (YCCO)
        DEALLOCATE  (ZCCO)
        !DEALLOCATE  (ZCC)
        DEALLOCATE  (XCC)
        DEALLOCATE  (XND)
        
        CALL CALL_TEC360
        
        CALL BULK_VELOCITY_TG(Umean1_tg)
        IF(MYID.EQ.0) call CHKRLHDL  ('            TG: The bulk velocity (orignal)=',MYID, Umean1_tg)
        
        DO J=1,N2DO(MYID)
            Q_tg(:,J,:,NFLOW)= Q_tg(:,J,:,NFLOW)/Umean1_tg
        END DO
        
        CALL BULK_VELOCITY_TG(Umean1_tg)
        IF(MYID.EQ.0) call CHKRLHDL  ('            TG: The bulk velocity (corcted)=',MYID,Umean1_tg)
        

        CALL MPI_BARRIER(ICOMM,IERROR)
    
        
        if(myid==0) CALL CHKRLHDL('       TG: END of interplating from a coarse mesh at time',myid,phyTIME_TG) 
        
        
        
        RETURN
    END SUBROUTINE 
    
    
    !=========================================================================================================    
    SUBROUTINE INITIAL_INTERP_io
        USE INTP_VARS
        IMPLICIT NONE
        
        CHARACTER(LEN=256)  :: WRT_RST_FNM
        CHARACTER(15) :: PNTIM
        INTEGER(4) :: I,J,K,N, KK, II, JJ
        real(wp)    :: Umean1_io, Gmean1_io, HMEAN, HMEAN_WORK, HMEANO_WORK, VOLOO, VOLOO_work
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:)
        REAL(WP),ALLOCATABLE :: QO   (:,:,:)
        REAL(WP),ALLOCATABLE :: QN   (:,:,:)
        
        if(myid==0) CALL CHKHDL('    IO: Interpolating from one mesh to another',myid) 
        
        !=========STEP 1, GENERATE MESH INFO==================
        CALL INITIAL_INTERP_MESH
        
        !=========STEP 2 ALLOCATE VARS===================
        ALLOCATE( DUMMY (NCLOO1,N2DOO(MYID),NCLOO3) )
        ALLOCATE( QO    (NCLO1, N2DOO(MYID),NCLO3)  )
        ALLOCATE( QN    (NCL1,  N2DO(MYID), NCL3 )  )
        
        !==========STEP 3 READ IN OLD DATA AND INTP INTO NEW ONES===============
        WRITE(PNTIM,'(1ES15.9)') RSTtim_tg
        DO N=1,NDV+1
            !==========(1) FILE NAME===========
            IF(N==1)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_U.D'
            IF(N==2)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_V.D'
            IF(N==3)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_W.D'
            IF(N==4)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_P.D'
            
            !==========(2) READ IN DATA===========
            CALL READ_3D_VARS(NCLOO1,NCLO2,NCLOO3,N2DOO(MYID),N2DOO(0)*MYID,ITERG0_TG,phyTIME_TG, DUMMY, WRT_RST_FNM)
            
            !==========(3) EXPAND READ-IN TO REAL LENGTH===========
            DO I=1,NCLO1
                DO J=1,N2DOO(MYID)
                    DO K=1,NCLO3
                    
                        II=MOD(I,NCLOO1)
                        IF(II==0) II=NCLOO1
                        
                        KK=MOD(K,NCLOO3)
                        IF(KK==0) KK=NCLOO3
                        
                        QO(I,J,K)=DUMMY(II,J,KK)
                    ENDDO
                ENDDO
            ENDDO
            
            !==========(4) INTEP DATA===========
            CALL INTP_3D_VARS(QO,QN,N)

            IF(RST_type_flg.EQ.0) THEN
                DO I=1,NCL1E
                    DO J=1,N2DO(MYID)
                        DO K=1,NCL3
                            IF(N.EQ.1) G_IO(I,J,K,1)=QN(I,J,K) !RHO U 
                            IF(N.EQ.2) G_IO(I,J,K,2)=QN(I,J,K) !RHO V
                            IF(N.EQ.3) G_IO(I,J,K,3)=QN(I,J,K) !RHO W
                            IF(N.EQ.4) PR_IO(I,J,K) =QN(I,J,K) !P
                        ENDDO
                    ENDDO
                ENDDO
                
                IF(N==NFLOW) THEN
                    CALL BULK_MASSFLUX_IO( Gmean1_io )
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk G (original)=      ',MYID, Gmean1_io)
                    
                    
                    DO J=1,N2DO(MYID)
                        G_io(:,J,:,NFLOW)= G_io(:,J,:,NFLOW)/Gmean1_io
                    END DO
                    
                    CALL BULK_MASSFLUX_IO( Gmean1_io )
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk G (corrected)=     ',MYID,Gmean1_io)
                    
                END IF
        
            END IF
            
            IF(RST_type_flg.EQ.1) THEN
            
                DO I=1,NCL1E
                    DO J=1,N2DO(MYID)
                        DO K=1,NCL3
                            IF(N.EQ.1) Q_IO(I,J,K,1)=QN(I,J,K) !U 
                            IF(N.EQ.2) Q_IO(I,J,K,2)=QN(I,J,K) !V
                            IF(N.EQ.3) Q_IO(I,J,K,3)=QN(I,J,K) !W
                            IF(N.EQ.4) PR_IO(I,J,K) =QN(I,J,K) !P
                        ENDDO
                    ENDDO
                ENDDO

                IF(N==NFLOW) THEN
                    CALL BULK_VELOCITY_IO( Umean1_io )
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk velocity(original)=',MYID, Umean1_io)
                    
                
                    DO J=1,N2DO(MYID)
                        Q_io(:,J,:,NFLOW)= Q_io(:,J,:,NFLOW)/Umean1_io
                    END DO
                
                    CALL BULK_VELOCITY_IO( Umean1_io )
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk velocity (coreted)=',MYID,Umean1_io)
                    
                END IF
                
            END IF

        END DO
        
        IF(thermlflg==1 .AND. RST_type_flg.EQ.0 ) THEN
            DO N=1,3
                !==========(1) FILE NAME===========
                IF(N==1)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_T.D'
                IF(N==2)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_D.D'
                IF(N==3)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath3)//'DNS_perixz_INSTANT_T'//TRIM(PNTIM)//'_E.D'
                
                !==========(2) READ IN DATA===========
                CALL READ_3D_VARS(NCLOO1,NCLO2,NCLOO3,N2DOO(MYID),N2DOO(0)*MYID,ITERG0_TG,phyTIME_TG, DUMMY, WRT_RST_FNM)
                
                !==========(3) EXPAND READ-IN TO REAL LENGTH===========
                DO I=1,NCLO1
                    DO J=1,N2DOO(MYID)
                        DO K=1,NCLO3
                        
                            II=MOD(I,NCLOO1)
                            IF(II==0) II=NCLOO1
                            
                            KK=MOD(K,NCLOO3)
                            IF(KK==0) KK=NCLOO3
                            
                            QO(I,J,K)=DUMMY(II,J,KK)
                        ENDDO
                    ENDDO
                ENDDO
                
                !======================================
                IF(N==3) THEN
                    DO J=1,N2DOO(MYID)  
                        JJ=myid*N2DOO(0) + J
                        DO I=1,NCLO1
                            DO K=1,NCLO3
                                HMEAN=HMEAN + QO(I,J,K)*(YNDO(JJ+1)-YNDO(JJ))*(HXO/DBLE(NCLO1))*(HZO/DBLE(NCLO3))
                                VOLOO=VOLOO + (YNDO(JJ+1)-YNDO(JJ))*(HXO/DBLE(NCLO1))*(HZO/DBLE(NCLO3))
                            END DO
                        END DO
                    ENDDO
                    CALL MPI_BARRIER(ICOMM,IERROR)
                    CALL MPI_ALLREDUCE(HMEAN, HMEANO_work,1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
                    CALL MPI_ALLREDUCE(VOLOO, VOLOO_work, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)  
                    
                    HMEANO_work = HMEANO_work/VOLOO_work
    
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk H orig-coarse mesh=',MYID, HMEANO_work)
                    
                END IF
                
                !==========(4) INTEP DATA===========
                CALL INTP_3D_VARS(QO,QN,0)
                
                DO J=1,N2DO(MYID)
                    DO I=1,NCL1
                        DO K=1,NCL3
                            IF(N.EQ.1) TEMPERATURE(I, J, K)  =QN(I,J,K) !T
                            IF(N.EQ.2) DENSITY    (I, J, K)  =QN(I,J,K) !D
                            IF(N.EQ.3) ENTHALPY   (I, J, K)  =QN(I,J,K) !E
                        ENDDO
                    ENDDO
                ENDDO
                !================================
                IF(N==3) THEN
                    CALL BULK_H_IO(HMEAN_work)
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk H orig-fine mesh  =',MYID, HMEAN_work)
                    
                END IF
        
            END DO
            
            !=================scaled to original==================
            DO J=1,N2DO(MYID)
                ENTHALPY(:,J,:)= ENTHALPY(:,J,:)*HMEANO_work/HMEAN_work
            END DO
            
            IF(N==3) THEN
                CALL BULK_H_IO(HMEAN_work)
                IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk H-scaled fine mesh=        ',MYID, HMEAN_work)
            END IF
            
        
            !=============BUILD UP RHOH===========================================
            DO J=1,N2DO(MYID)
                DO I=1,NCL1E
                    DO K=1,NCL3
                        RHOH(I,J,K) = DENSITY(I,J,K) * ENTHALPY(I,J,K)
                    END DO
                END DO
            END DO
            
            !===========build up thermal BC and other thermal variables=======
            CALL BC_WALL_THERMAL(IALLDM)
            IF(TGFLOWFLG) CALL BC_TINLET_THERML
            
            CALL THERM_PROP_UPDATE(IALLDM)
            IF(TGFLOWFLG) THEN
                CALL INTFC_OUL_THERMAL_io
                CALL INTFC_INL_THERMAL_io
            END IF
            CALL INTFC_MFD_THERMAL_io
            
            
            !=============update velocity based on the given G and density===================
            CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
            CALL BC_WALL_G_io
            CALL INTFC_VARS1(NCL1S,NCL1E,NCL1S,NCL1E,PR_io)
            CALL BC_WALL_Pr_io
            IF(TGFLOWFLG) THEN
                CALL BC_TINLET_FLOW
            END IF
            CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
            CALL INTFC_VARS1(NCL1S,NCL1E,NCL1S,NCL1E,PR_io)
            CALL BC_WALL_G_io
            CALL BC_WALL_PR_io
            
            CALL DENSITY_Staggered
            CALL MU_Staggered
            CALL VELOCITY_CALC_io
        
        
        ELSE
            !=======build up thermal field=================
            !CALL INIFIELD_THERMAL_io 
            IF(thermlflg ==1) THEN
                CALL BC_WALL_THERMAL(IALLDM)
                IF(TGFLOWFLG) CALL BC_TINLET_THERML
                CALL THERM_PROP_UPDATE(IALLDM)
                IF(TGFLOWFLG) THEN
                    CALL INTFC_OUL_THERMAL_io
                    CALL INTFC_INL_THERMAL_io
                END IF
                CALL INTFC_MFD_THERMAL_io
            END IF
            
            !==========flow field===================================
            IF(RST_type_flg.EQ.1) THEN
                !====================
                IF(TGFLOWFLG) THEN
                    CALL BC_TINLET_FLOW
                END IF
                CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
                CALL INTFC_VARS1(NCL1S,NCL1E,NCL1S,NCL1E,PR_io)
                CALL BC_WALL_PR_io
                CALL BC_WALL_Q_io
            
                CALL VELO2MASSFLUX
                 
                CALL Unified_massflux
                 
                CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
                CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
                CALL BC_WALL_G_io
                CALL BC_WALL_Q_io
                
                phyTIME_io = 0.0_WP
                ITERG0_io  = 0
                
            END IF
            
            IF(RST_type_flg.EQ.0) THEN
                IF(TGFLOWFLG) THEN
                    CALL BC_TINLET_FLOW
                END IF
                CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
                CALL INTFC_VARS1(NCL1S,NCL1E,NCL1S,NCL1E,PR_io)
                CALL BC_WALL_G_io
                CALL BC_WALL_PR_io
                
                CALL DENSITY_Staggered
                CALL MU_Staggered
                CALL VELOCITY_CALC_io
                
                CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
                CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,Q_io)
                CALL BC_WALL_G_io
                CALL BC_WALL_Q_io
            END IF
        
        END IF
        
        DEALLOCATE (DUMMY)
        DEALLOCATE (QO)
        DEALLOCATE (QN)
        DEALLOCATE  (N2DOO)
        DEALLOCATE  (N3DOO)
        DEALLOCATE  (XNDO)
        DEALLOCATE  (YNDO)
        DEALLOCATE  (ZNDO)
        DEALLOCATE  (XCCO)
        DEALLOCATE  (YCCO)
        DEALLOCATE  (ZCCO)
        !DEALLOCATE  (ZCC)
        
        CALL CALL_TEC360
        
        
        if(myid==0) CALL CHKRLHDL('       IO: END of interplating from a coarse mesh at time',myid,phyTIME_IO) 
        
        
        
        RETURN
    END SUBROUTINE 
    
