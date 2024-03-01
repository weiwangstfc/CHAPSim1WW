    SUBROUTINE RESTART_INSTANT_VARS_IO(TFM)  ! FOR BOTH KINDS OF DOMAINS....
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE
        
        REAL(WP), INTENT(IN) :: TFM
        CHARACTER(LEN=256)  :: WRT_RST_FNM
        CHARACTER(LEN=64)   :: FLNAME
        CHARACTER(15) :: PNTIM
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:)
        INTEGER(4) :: I,J,K,N
        REAL(wp)   :: Gmean1_io, Umean1_io, Hmean1
        
        
        
        if(myid==0) CALL CHKHDL('IO: Restart instantanous flow field',myid) 
        
        ALLOCATE(DUMMY(NCL1E,N2DO(MYID),NCL3))
        WRITE(PNTIM,'(1ES15.9)') TFM !RSTtim_io
        
        IF(TGFLOWFLG .and. IOFLOWFLG) THEN
            FLNAME ='DNS_perioz_INSTANT_T'
        ELSE 
            FLNAME ='DNS_perixz_INSTANT_T'
        END IF
        
        DO N=1,NDV+1
            IF(N==1)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_U.D'
            IF(N==2)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_V.D'
            IF(N==3)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_W.D'
            IF(N==4)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_P.D'
            
            CALL READ_3D_VARS(NCL1E,NCL2,NCL3,N2DO(MYID),JCL2G(1)-1,ITERG0_IO,phyTIME_IO, DUMMY, WRT_RST_FNM)
            
            IF(RST_type_flg.EQ.0) THEN
                DO I=1,NCL1E
                    DO J=1,N2DO(MYID)
                        DO K=1,NCL3
                            IF(N.EQ.1) G_IO(I,J,K,1)=DUMMY(I,J,K) !RHO U 
                            IF(N.EQ.2) G_IO(I,J,K,2)=DUMMY(I,J,K) !RHO V
                            IF(N.EQ.3) G_IO(I,J,K,3)=DUMMY(I,J,K) !RHO W
                            IF(N.EQ.4) PR_IO(I,J,K) =DUMMY(I,J,K) !P
                        ENDDO
                    ENDDO
                ENDDO
                
                IF(N==NFLOW) THEN
                    CALL BULK_MASSFLUX_IO(Gmean1_io)
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk massflux (restart)=',MYID, Gmean1_io)
                END IF
                
            END IF
            
            IF(RST_type_flg.EQ.1) THEN
                DO I=1,NCL1E
                    DO J=1,N2DO(MYID)
                        DO K=1,NCL3
                            IF(N.EQ.1) Q_IO(I,J,K,1)=DUMMY(I,J,K) !U 
                            IF(N.EQ.2) Q_IO(I,J,K,2)=DUMMY(I,J,K) !V
                            IF(N.EQ.3) Q_IO(I,J,K,3)=DUMMY(I,J,K) !W
                            IF(N.EQ.4) PR_IO(I,J,K) =DUMMY(I,J,K) !P
                        ENDDO
                    ENDDO
                ENDDO
                IF(N==NFLOW) THEN
                    CALL BULK_VELOCITY_IO(Umean1_io)
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk velocity (restart)=',MYID, Umean1_io)
                END IF
            END IF

        END DO
            

        if(myid==0) CALL CHKRLHDL('       IO: End of restarting instantanous flow field',myid,phyTIME_IO)
        !============================thermal field==============================================================
        IF(thermlflg==1 .AND. RST_type_flg.EQ.0 ) THEN
            if(myid==0) CALL CHKHDL('    IO: Restart instantanous thermal field',myid) 
            DO N=1,3
                IF(N==1)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_T.D'
                IF(N==2)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_D.D'
                IF(N==3)  WRITE(WRT_RST_FNM,'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_E.D'
        
                CALL READ_3D_VARS(NCL1E,NCL2,NCL3,N2DO(MYID),JCL2G(1)-1,ITERG0_IO,phyTIME_IO, DUMMY, WRT_RST_FNM)
                
                DO I=1,NCL1E
                    DO J=1,N2DO(MYID)
                        DO K=1,NCL3
                            IF(N.EQ.1) TEMPERATURE(I, J, K)  =DUMMY(I,J,K)
                            IF(N.EQ.2) DENSITY    (I, J, K)  =DUMMY(I,J,K)
                            IF(N.EQ.3) ENTHALPY   (I, J, K)  =DUMMY(I,J,K)
                        ENDDO
                    ENDDO
                ENDDO
                
                IF(N==3) THEN
                    CALL BULK_H_IO(Hmean1)
                    IF(MYID.EQ.0) call CHKRLHDL  ('            IO: The bulk H (restart) =      ',MYID, Hmean1) 
                END IF
                
            END DO
            if(myid==0) CALL CHKRLHDL('       IO: End of restarting instantanous thermal field',myid,phyTIME_IO) 
            
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
                CALL INIFIELD_THERMAL_io
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
        
        if(myid==0) CALL CHKRLHDL('===========phyTIME_io==================',myid,phyTIME_IO) 
        RETURN
    END SUBROUTINE

!****************************************************************
    SUBROUTINE RESTART_AVERAGE_VARS_nonXperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE
        
        CHARACTER(15) :: PNTIM
        INTEGER(4) :: DFLG, NARRAY
        INTEGER(4) :: IERR,SIZES_ARRAY(4),NEWTYPE,SUBSIZES(4),STARTS(4) 
        INTEGER(4) :: INISIZE,IRLSIZE,INTMPI(4),INTBYTE,DBLBYTE
        INTEGER(4) :: BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND) OFFSET
        INTEGER(4) :: N1ML, N2ML, I, J, L, L1, L2, N, M, H, P
        REAL(WP) :: RLEMPI(3)
        REAL(WP) :: RENL, DTL
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:,:)
        logical  :: file_exists
        
        if(myid==0) CALL CHKHDL('    IO: Restart averaged flow field',myid) 
        
        NARRAY = 4
        ALLOCATE( DUMMY(1:NCL1_io, 1:N2DO(MYID), 1:(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8), 1:21) )
        DUMMY = 0.0_wp
        
        DFLG = 100
        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
        WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perioz_AVERAGD_T'//TRIM(PNTIM)//'_FLOW.D'
        
        INQUIRE(FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)), EXIST=file_exists) 
        if(.not.file_exists .and. myid==0) CALL ERRHDL('File '//WRT_AVE_FNM_io//' does not exist.',0)

        SIZES_ARRAY(1)=NCL1_io    
        SIZES_ARRAY(2)=NCL2
        SIZES_ARRAY(3)=NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8
        SIZES_ARRAY(4)=21
        
        SUBSIZES(1)=SIZES_ARRAY(1)
        SUBSIZES(2)=N2DO(MYID)
        SUBSIZES(3)=SIZES_ARRAY(3)
        SUBSIZES(4)=SIZES_ARRAY(4)
        
        STARTS(1)=0
        STARTS(2)=JCL2G(1)-1  !SUBSIZES(2) * MYID
        STARTS(3)=0
        STARTS(4)=0
        
        BUFSIZE=SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)* SUBSIZES(4)

        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY,MPI_INFO_NULL,MYFILE,IERR)
        
        INISIZE=4
        IRLSIZE=3
        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_READ(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
       
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        
        
        N1ML     =INTMPI(1)
        N2ML     =INTMPI(2)
        ITERG0_IO=INTMPI(3)
        NSTATIS_IO=INTMPI(4)
        
        phyTIME_IO=RLEMPI(1)
        RENL      =RLEMPI(2)  
        DTL       =RLEMPI(3)
        
        IF(MYID==0) THEN
            CALL CHKINTHDL('        I IN AVERAGED VARIABLES (I,J,L,M) ',MYID,N1ML)
            CALL CHKINTHDL('        J IN AVERAGED VARIABLES (I,J,L,M) ',MYID,N2ML)
            CALL CHKINTHDL('        ITERG0_IO IN AVERAGED VARIABLES ',MYID,ITERG0_IO)
            CALL CHKINTHDL('        NSTATIS_IO IN AVERAGED VARIABLES',MYID,NSTATIS_IO)
            
            CALL CHKRLHDL ('        phyTIME_IO IN AVERAGED VARIABLES',MYID,phyTIME_IO)
            CALL CHKRLHDL ('        RENL       IN AVERAGED VARIABLES',MYID,RENL)
            CALL CHKRLHDL ('        DTL        IN AVERAGED VARIABLES',MYID,DTL)
        END IF
        
        L1 = NDV+1
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    U1ztL_io(I,J,L)=DUMMY(I,J,L,1)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    G1ztL_io(I,J,L)=DUMMY(I,J,L,2)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    UPztL_io(I,J,L)=DUMMY(I,J,L,3)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV*(7-NDV)/2+NDV-3
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    U2ztL_io(I,J,L)=DUMMY(I,J,L,4)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV*(7-NDV)/2+NDV-3
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    UGztL_io(I,J,L)=DUMMY(I,J,L,5)
                ENDDO
            ENDDO
        END DO

        L1 = NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    UGUztL_io(I,J,L)=DUMMY(I,J,L,6)
                ENDDO
            ENDDO
        ENDDO
        
        
        L1 = NDV
        L2 = NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DO M=1,L2
                        DVDL1ztL_io(I, J, L, M)=DUMMY(I,J,L,M+6)
                    END DO
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV
        L2 = NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DO M=1,L2
                        DVDLPztL_io(I, J, L, M)=DUMMY(I,J,L,M+9)
                    END DO
                ENDDO
            ENDDO
        END DO
        
        DO M=1,NDV
            DO H=1,NDV
                DO N=1,NDV
                    DO P=1,NDV
                        L1=(M-1)*3+H
                        L2=(N-1)*3+P
                        DO J=1,N2DO(MYID)
                            DVDL2ztL_io(I,J,L1,L2)=DUMMY(I,J,L1,L2+12)
                        END DO
                    END DO
                END DO
            END DO
        END DO
        
        DEALLOCATE (DUMMY)
        
        IF(RST_type_flg.EQ.1) THEN
            phyTIME_io = 0.0_WP
            ITERG0_io  = 0
        END IF
        

        IF(thermlflg==1 .AND. RST_type_flg.EQ.0 ) call RESTART_AVERAGE_VARS_THERMAL_nonXperiodic_IO

        if(myid==0) CALL CHKRLHDL('       IO: End of restarting averaged flow field',myid,phyTIME_IO) 
        
        RETURN
    
    END SUBROUTINE
    
    
!****************************************************************
    SUBROUTINE RESTART_AVERAGE_VARS_THERMAL_nonXperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE
        
        CHARACTER(15) :: PNTIM
        INTEGER(4) :: DFLG, NARRAY
        INTEGER(4) :: IERR,SIZES_ARRAY(3),NEWTYPE,SUBSIZES(3),STARTS(3) 
        INTEGER(4) :: INISIZE,IRLSIZE,INTMPI(4),INTBYTE,DBLBYTE
        INTEGER(4) :: BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND) OFFSET
        INTEGER(4) :: N1ML, N2ML, I, J, L1, L,N
        REAL(WP) :: RLEMPI(3)
        REAL(WP) :: RENL, DTL
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:)
        logical  :: file_exists
        
        if(myid==0) CALL CHKHDL('    IO: Restart averaged thermal field',myid) 
        
        NARRAY = 3
        ALLOCATE( DUMMY(1:NCL1_io, 1:N2DO(MYID), 7+2*NDV ) )
        DUMMY = 0.0_wp
        
        DFLG = 100
        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
        WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perioz_AVERAGD_T'//TRIM(PNTIM)//'_THEL.D'
        
        INQUIRE(FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)), EXIST=file_exists) 
        if(.not.file_exists .and. myid==0) CALL ERRHDL('File '//WRT_AVE_FNM_io//' does not exist.',0)
        

        SIZES_ARRAY(1)=NCL1_io    
        SIZES_ARRAY(2)=NCL2
        SIZES_ARRAY(3)=7+2*NDV
        
        SUBSIZES(1)=SIZES_ARRAY(1)
        SUBSIZES(2)=N2DO(MYID)
        SUBSIZES(3)=SIZES_ARRAY(3)
        
        STARTS(1)=0
        STARTS(2)=JCL2G(1)-1  !SUBSIZES(2) * MYID
        STARTS(3)=0
        
        BUFSIZE=SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)

        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY,MPI_INFO_NULL,MYFILE,IERR)
        
        INISIZE=4
        IRLSIZE=3
        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_READ(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        
        
        N1ML     =INTMPI(1)
        N2ML     =INTMPI(2)
        ITERG0_IO=INTMPI(3)
        NSTATIS_IO=INTMPI(4)
        
        phyTIME_IO=RLEMPI(1)
        RENL      =RLEMPI(2)  
        DTL       =RLEMPI(3)
        
        N  = 1
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                T1ztL_io(I,J)=DUMMY(I,J,N)
            ENDDO
        END DO
        

        N  = 2
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                D1ztL_io(I,J)=DUMMY(I,J,N)
            ENDDO
        END DO
        
        N  = 3
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                H1ztL_io(I,J)=DUMMY(I,J,N)
            ENDDO
        END DO
        
        N  = 4
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                T2ztL_io(I,J)=DUMMY(I,J,N)
            ENDDO
        END DO
        

        N  = 5
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                D2ztL_io(I,J)=DUMMY(I,J,N)
            ENDDO
        END DO
        
        N  = 6
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                H2ztL_io(I,J)=DUMMY(I,J,N)
            ENDDO
        END DO
        
        N  = 7
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DHztL_io(I,J)=DUMMY(I,J,N)
            ENDDO
        END DO
        
        L1 = NDV
        N  = 7
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    UHztL_io(I, J, L)=DUMMY(I,J,L+N)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV
        N  = 7+NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    GHztL_io(I, J, L)=DUMMY(I,J,L+N)
                ENDDO
            ENDDO
        END DO
        
        DEALLOCATE (DUMMY)
        
        IF(RST_type_flg.EQ.1) THEN
            phyTIME_io = 0.0_WP
            ITERG0_io  = 0
        END IF
        
        if(myid==0) CALL CHKRLHDL('       IO: End of restarting averaged thermal field',myid, phyTIME_IO) 

        RETURN
    
    END SUBROUTINE    
    
    
    
    !****************************************************************
    SUBROUTINE RESTART_AVERAGE_VARS_Xperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE
        
        CHARACTER(15) :: PNTIM
        INTEGER(4) :: DFLG
        INTEGER(4) :: IERR,NEWTYPE
        INTEGER(4) :: INISIZE,IRLSIZE,INTBYTE,DBLBYTE
        INTEGER(4)    :: NARRAY, NSZ
        INTEGER(4),ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),STARTS(:)
        INTEGER(4),ALLOCATABLE :: INTMPI(:)
        REAL(WP),ALLOCATABLE :: RLEMPI(:)
        INTEGER(4) :: BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND) OFFSET
        INTEGER(4) :: N1ML, N2ML, J, L, L1, L2, N, NLLL, NNNL,M, K
        REAL(WP) :: RENL, DTL
        REAL(WP),ALLOCATABLE :: DUMMY(:,:), DUMMY1(:,:), DUMMY2(:,:)
        logical  :: file_exists
        INTEGER(4) :: NSTATIS_IO1, NSTATIS_IO2
        
        if(myid==0) CALL CHKHDL('    IO: Restart averaged flow field',myid) 
        
        !============================
        NARRAY=2
        ALLOCATE ( SIZES_ARRAY(NARRAY) )
        ALLOCATE ( SUBSIZES   (NARRAY) )
        ALLOCATE ( STARTS     (NARRAY) )
        INISIZE=4
        IRLSIZE=3
        ALLOCATE ( INTMPI(INISIZE)       )
        ALLOCATE ( RLEMPI(IRLSIZE)       )
        
        ! original 
        NSZ = NDV+1 + 2*NDV +2*(NDV*(7-NDV)/2+NDV-3) + (NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8) + &
              2*NDV*NDV + ((NDV-1)*3+NDV)*((NDV-1)*3+NDV)
              
        ! if with U3
        NSZ = NSZ + (NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
        
        ! if with quadrant
        NSZ = NSZ + 6*(4*QUADHN)
        
        ! if with driven force
        NSZ = NSZ + NDV+1
        
        ! if with octant
        !NSZ = NSZ + 12*(8*QUADHN)
        
              
        ALLOCATE( DUMMY (1:N2DO(MYID), NSZ ) )
        ALLOCATE( DUMMY1(1:N2DO(MYID), NSZ ) )
        DUMMY  = 0.0_wp
        DUMMY1 = 0.0_wp
        
        
        !=======prepare sizes and ranges===================================================
        SIZES_ARRAY(1)=NCL2
        SIZES_ARRAY(2)=NSZ
        
        SUBSIZES(1)=N2DO(MYID)
        SUBSIZES(2)=SIZES_ARRAY(2)
        
        STARTS(1)=JCL2G(1)-1  !SUBSIZES(2) * MYID
        STARTS(2)=0
        
        BUFSIZE=SUBSIZES(1) * SUBSIZES(2)
        
        !===========read data1==============================================================
        DFLG = 100
        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
        WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_FLOW.D'
        
        if(myid==0) CALL CHKHDL('    IO: Restart averaged flow field '//TRIM(WRT_AVE_FNM_io),myid) 
        
        INQUIRE(FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)), EXIST=file_exists) 
        if(.not.file_exists .and. myid==0) CALL ERRHDL('File '//WRT_AVE_FNM_io//' does not exist.',0)
        
        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY,MPI_INFO_NULL,MYFILE,IERR)

        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_READ(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
       
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ_ALL(MYFILE,DUMMY1,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        
        
        N2ML       =INTMPI(1)
        NNNL       =INTMPI(2)
        ITERG0_IO  =INTMPI(3)
        NSTATIS_IO1=INTMPI(4)
        
        phyTIME_IO =RLEMPI(1)
        RENL       =RLEMPI(2)  
        DTL        =RLEMPI(3)
        
        IF(MYID==0) THEN
            CALL CHK2INTHDL('        Expected and read-in Dummy Sizes are',MYID,NSZ,NNNL)
            IF(NSZ.NE.NNNL) CALL ERRHDL('Restarting averaged data formats are not consistent! ',myid)
            CALL CHKINTHDL('        J IN AVERAGED VARIABLES (J,M) ',MYID,N2ML)
            CALL CHKINTHDL('        M IN AVERAGED VARIABLES (J,M) ',MYID,NNNL)
            CALL CHKINTHDL('        ITERG0_IO IN AVERAGED VARIABLES ',MYID,ITERG0_IO)
            CALL CHKINTHDL('        NSTATIS_IO IN AVERAGED VARIABLES',MYID,NSTATIS_IO1)
            
            CALL CHKRLHDL ('        phyTIME_IO IN AVERAGED VARIABLES',MYID,phyTIME_IO)
            CALL CHKRLHDL ('        RENL       IN AVERAGED VARIABLES',MYID,RENL)
            CALL CHKRLHDL ('        DTL        IN AVERAGED VARIABLES',MYID,DTL)
        END IF
        
        NSTATIS_IO = NSTATIS_IO1
        DUMMY = DUMMY1
        
        !===========read data2==============================================================
        IF(TSTAV_RESET.GT.TSTAV1+1.0E-2_WP) THEN
            ALLOCATE( DUMMY2(1:N2DO(MYID), NSZ ) )
            DUMMY2 = 0.0_wp
            
            DFLG = 101
            WRITE(PNTIM,'(1ES15.9)') TSTAV_RESET
            WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_FLOW.D'
        
            if(myid==0) CALL CHKHDL('    IO: Resetting averaged flow field '//TRIM(WRT_AVE_FNM_io),myid) 
        
            INQUIRE(FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)), EXIST=file_exists) 
            if(.not.file_exists .and. myid==0) CALL ERRHDL('File '//WRT_AVE_FNM_io//' does not exist.',0)
        
            CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY,SUBSIZES,STARTS,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,NEWTYPE,IERR)
            CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
            CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY,MPI_INFO_NULL,MYFILE,IERR)

            OFFSET=0_MPI_OFFSET_KIND
            CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
            CALL MPI_FILE_READ(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
            CALL MPI_FILE_READ(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
       
            CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
            CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
            OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
            CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
            CALL MPI_FILE_READ_ALL(MYFILE,DUMMY2,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
            CALL MPI_FILE_CLOSE(MYFILE,IERR)
            CALL MPI_TYPE_FREE(NEWTYPE,IERR)
            

            NSTATIS_IO2=INTMPI(4)
            IF(MYID==0) THEN
                IF(NSZ.NE.INTMPI(2)) CALL ERRHDL('Restarting averaged data formats are not consistent! ',myid)
                CALL CHK2INTHDL('        Expected and read-in Dummy Sizes are',MYID,NSZ,INTMPI(2))
                CALL CHKINTHDL('        J IN AVERAGED VARIABLES (J,M) ',MYID,INTMPI(1))
                CALL CHKINTHDL('        M IN AVERAGED VARIABLES (J,M) ',MYID,INTMPI(2))
                CALL CHKINTHDL('        ITERG0_IO IN AVERAGED VARIABLES ',MYID,INTMPI(3))
                CALL CHKINTHDL('        NSTATIS_IO IN AVERAGED VARIABLES',MYID,NSTATIS_IO2)
                
                CALL CHKRLHDL ('        phyTIME_IO IN AVERAGED VARIABLES',MYID,RLEMPI(1))
                CALL CHKRLHDL ('        RENL       IN AVERAGED VARIABLES',MYID,RLEMPI(2))
                CALL CHKRLHDL ('        DTL        IN AVERAGED VARIABLES',MYID,RLEMPI(3))
            END IF
            
            DUMMY = (DUMMY1*DBLE(NSTATIS_IO1)-DUMMY2*DBLE(NSTATIS_IO2))/DBLE(NSTATIS_IO1-NSTATIS_IO2)
            NSTATIS_IO = NSTATIS_IO1-NSTATIS_IO2
            DEALLOCATE (DUMMY2)
        END IF
        
        !==========================================
        N=0
        L1 = NDV+1
        DO J=1,N2DO(MYID)
            DO L=1,L1
                U1xztL_io(J,L)=DUMMY(J,N+L) !1-4
            ENDDO
        ENDDO
        N = N + L1 !4
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in U1xztL_io',myid, L1, N) 
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1
                G1xztL_io(J,L)=DUMMY(J,N+L)!5-7
            ENDDO
        ENDDO
        N = N + L1 !7
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in G1xztL_io',myid, L1, N) 
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1
                UPxztL_io(J,L)=DUMMY(J,N+L) !8-10
            ENDDO
        ENDDO
        N = N + L1 !10
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in UPxztL_io',myid, L1, N) 
        
        !==========================================
        L1 = NDV*(7-NDV)/2+NDV-3
        DO J=1,N2DO(MYID)
            DO L=1,L1
                U2xztL_io(J,L)=DUMMY(J,N+L) !11-16
            ENDDO
        ENDDO
        N = N + L1 !16
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in U2xztL_io',myid, L1, N) 
        
        L1 = NDV*(7-NDV)/2+NDV-3
        DO J=1,N2DO(MYID)
            DO L=1,L1
                UGxztL_io(J,L)=DUMMY(J,N+L)
            ENDDO
        ENDDO
        N = N + L1 !22
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in UGxztL_io',myid, L1, N) 
        
        
        L1 = NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8 !10
        DO J=1,N2DO(MYID)
            DO L=1,L1
                UGUxztL_io(J,L)=DUMMY(J,N+L)
            ENDDO
        ENDDO
        N = N + L1 !32
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in UGUxztL_io',myid, L1, N) 
        
!!!!        ! Commented below is no U3 ===========================
        L1 = NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8 !10
        DO J=1,N2DO(MYID)
            DO L=1,L1
                U3xztL_io(J,L)=DUMMY(J,N+L)
            ENDDO
        ENDDO
        N = N + L1 !32
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in U3xztL_io',myid, L1, N) 
        
        !==========================================
        L1 = NDV
        L2 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DVDL1xztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in DVDL1xztL_io',myid, L1, N) 
        
        
        L1 = NDV
        L2 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DVDLPxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in DVDLPxztL_io',myid, L1, N) 
        
        L1= (NDV-1)*3+NDV !9
        L2= (NDV-1)*3+NDV !9
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DVDL2xztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in DVDL2xztL_io',myid, L1, N) 
        
!        !        ! Commented below is no quadrant
        !=============quadrant===============
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    QUADUVxztL_io(J, L, K)=DUMMY(J,M)
                    !write(*,*) QUADUVxztL_io(J, L, K) !!!test
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in QUADUVxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    QUADVzxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in QUADVzxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    QUADTKxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in QUADTKxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    QUADDRxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in QUADDRxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    QUADDUV1xztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in QUADDUV1xztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    QUADDUV2xztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in QUADDUV2xztL_io',myid, L1, N) 
        
        !=========================================================
        ! with driven force
        L1 = NDV+1
        DO J=1,N2DO(MYID)
            DO L=1,L1
                FUxztL_IO(J,L)=DUMMY(J,N+L) !1-4
            ENDDO
        ENDDO
        N = N + L1 !4
        !IF(MYID==0) CALL CHK2INTHDL('        Reading in FUxztL_io',myid, L1, N) 
        
!        !=======================================================
!        !=============octant==\rho=============
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTDUVxztL_io(J, L, K)=DUMMY(J,M)
!                    !write(*,*) OCTDUVxztL_io(J, L, K) !!!test
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTDUVxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTDVzxztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTDVzxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTDTKxztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTDTKxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTDDRxztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTDDRxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTDDUV1xztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTDDUV1xztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTDDUV2xztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTDDUV2xztL_io',myid, L1, N) 
        
!        !=============octant==T=============
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTTUVxztL_io(J, L, K)=DUMMY(J,M)
!                    !write(*,*) OCTTUVxztL_io(J, L, K) !!!test
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTTUVxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTTVzxztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTTVzxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTTTKxztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTTTKxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTTDRxztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTTDRxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTTDUV1xztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTTDUV1xztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    OCTTDUV2xztL_io(J, L, K)=DUMMY(J,M)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Reading in OCTTDUV2xztL_io',myid, L1, N) 
        
        
        
        IF(N.NE.NSZ) THEN
            WRITE(*,*) '#ERROR in RESTART_AVERAGE_VARS_Xperiodic_IO', N, NSZ
        END IF
        
        DEALLOCATE (DUMMY)
        DEALLOCATE (DUMMY1)
        if(myid==0) CALL CHKRLHDL('       IO: End of restarting averaged flow    field',myid,phyTIME_IO) 
        
        IF(thermlflg==1 .AND. RST_type_flg.EQ.0 ) call RESTART_AVERAGE_VARS_THERMAL_Xperiodic_IO

        IF(RST_type_flg.EQ.0 )  CALL RESTART_AVERAGE_SPECTRA
        
        RETURN
    
    END SUBROUTINE
    
    
!****************************************************************
    SUBROUTINE RESTART_AVERAGE_VARS_THERMAL_Xperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE
        
        CHARACTER(15) :: PNTIM
        INTEGER(4) :: DFLG
        INTEGER(4) :: IERR,NEWTYPE
        INTEGER(4) :: INISIZE,IRLSIZE,INTBYTE,DBLBYTE
        INTEGER(4)    :: NARRAY, NSZ
        INTEGER(4),ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),STARTS(:)
        INTEGER(4),ALLOCATABLE :: INTMPI(:)
        REAL(WP),ALLOCATABLE :: RLEMPI(:)
        INTEGER(4) :: BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND) OFFSET
        INTEGER(4) :: N1ML, N2ML, J, L1, L2, L3, L,N, NNNL, K, S, M
        REAL(WP) :: RENL, DTL
        REAL(WP),ALLOCATABLE :: DUMMY(:,:), DUMMY1(:,:), DUMMY2(:,:)
        logical :: file_exists
        INTEGER(4) :: NSTATIS_IO1, NSTATIS_IO2
        
        
        
        
        !============================
        NARRAY=2
        ALLOCATE ( SIZES_ARRAY(NARRAY) )
        ALLOCATE ( SUBSIZES   (NARRAY) )
        ALLOCATE ( STARTS     (NARRAY) )
        INISIZE=4
        IRLSIZE=3
        ALLOCATE ( INTMPI(INISIZE)       )
        ALLOCATE ( RLEMPI(IRLSIZE)       )


        !!NSZ = 9+5*NDV+((NDV*(7-NDV))/2+NDV-3)+4*NDV*NDV+3*NDV*NDV*NDV
        NSZ = 9+5*NDV+((NDV*(7-NDV))/2+NDV-3)+3*NDV*NDV+3*NDV*NDV*NDV+((NDV-1)*NDV+NDV)*((NDV-1)*NDV+NDV)
        ALLOCATE( DUMMY (1:N2DO(MYID), NSZ ) )
        ALLOCATE( DUMMY1(1:N2DO(MYID), NSZ ) )
        DUMMY  = 0.0_wp
        DUMMY1 = 0.0_wp
        
        !============================
        SIZES_ARRAY(1)=NCL2
        SIZES_ARRAY(2)=NSZ
        
        SUBSIZES(1)=N2DO(MYID)
        SUBSIZES(2)=SIZES_ARRAY(2)
        
        STARTS(1)=JCL2G(1)-1  !SUBSIZES(2) * MYID
        STARTS(2)=0
        
        BUFSIZE=SUBSIZES(1) * SUBSIZES(2)
        
        INISIZE=4
        IRLSIZE=3
        !============================
        
        !!===========read data1==============================================================
        DFLG = 100
        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
        WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_THEL.D'
        
        if(myid==0) CALL CHKHDL('    IO: Restart averaged thermal field'//WRT_AVE_FNM_io,myid) 
        
        INQUIRE(FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)), EXIST=file_exists) 
        if(.not.file_exists .and. myid==0) CALL ERRHDL('File '//WRT_AVE_FNM_io//' does not exist.',0)

        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY, SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE,IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY,MPI_INFO_NULL,MYFILE,IERR)
        
        
        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_READ(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_READ_ALL(MYFILE,DUMMY1,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        
        
        N2ML       =INTMPI(1)
        NNNL       =INTMPI(2)
        ITERG0_IO  =INTMPI(3)
        NSTATIS_IO1=INTMPI(4)
        
        phyTIME_IO=RLEMPI(1)
        RENL      =RLEMPI(2)  
        DTL       =RLEMPI(3)
        
        IF(MYID==0) THEN
            CALL CHK2INTHDL('        Expected and read-in Dummy Sizes are',MYID,NSZ,NNNL)
            CALL CHKINTHDL('        J IN AVERAGED VARIABLES (J,M) ',MYID,N2ML)
            CALL CHKINTHDL('        M IN AVERAGED VARIABLES (J,M) ',MYID,NNNL)
            CALL CHKINTHDL('        ITERG0_IO IN AVERAGED VARIABLES ',MYID,ITERG0_IO)
            CALL CHKINTHDL('        NSTATIS_IO IN AVERAGED VARIABLES',MYID,NSTATIS_IO1)
            
            CALL CHKRLHDL ('        phyTIME_IO IN AVERAGED VARIABLES',MYID,phyTIME_IO)
            CALL CHKRLHDL ('        RENL       IN AVERAGED VARIABLES',MYID,RENL)
            CALL CHKRLHDL ('        DTL        IN AVERAGED VARIABLES',MYID,DTL)
        END IF
        
        NSTATIS_IO = NSTATIS_IO1
        DUMMY = DUMMY1
        
        !!===========read data2==============================================================
        IF(TSTAV_RESET.GT.TSTAV1+1.0E-2_WP) THEN
            ALLOCATE( DUMMY2 (1:N2DO(MYID), NSZ ) )
            DUMMY2 = 0.0_wp
            
            DFLG = 101
            WRITE(PNTIM,'(1ES15.9)') TSTAV_RESET
            WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_THEL.D'
            
            if(myid==0) CALL CHKHDL('    IO: Resetting averaged thermal field'//WRT_AVE_FNM_io,myid) 
            
            INQUIRE(FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)), EXIST=file_exists) 
            if(.not.file_exists .and. myid==0) CALL ERRHDL('File '//WRT_AVE_FNM_io//' does not exist.',0)
    
            CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY,SUBSIZES,STARTS,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,NEWTYPE,IERR)
            CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
            CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io, MPI_MODE_RDONLY,MPI_INFO_NULL,MYFILE,IERR)
            
            
            OFFSET=0_MPI_OFFSET_KIND
            CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
            CALL MPI_FILE_READ(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
            CALL MPI_FILE_READ(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
            CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
            CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
            
            OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
            CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
            CALL MPI_FILE_READ_ALL(MYFILE,DUMMY2,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
            CALL MPI_FILE_CLOSE(MYFILE,IERR)
            CALL MPI_TYPE_FREE(NEWTYPE,IERR)
            
            NSTATIS_IO2=INTMPI(4)

            IF(MYID==0) THEN
                CALL CHK2INTHDL('        Expected and read-in Dummy Sizes are',MYID,NSZ,INTMPI(2))
                CALL CHKINTHDL('        J IN AVERAGED VARIABLES (J,M) ',MYID,INTMPI(1))
                CALL CHKINTHDL('        M IN AVERAGED VARIABLES (J,M) ',MYID,INTMPI(2))
                CALL CHKINTHDL('        ITERG0_IO IN AVERAGED VARIABLES ',MYID,INTMPI(3))
                CALL CHKINTHDL('        NSTATIS_IO IN AVERAGED VARIABLES',MYID,NSTATIS_IO2)
                
                CALL CHKRLHDL ('        phyTIME_IO IN AVERAGED VARIABLES',MYID,RLEMPI(1))
                CALL CHKRLHDL ('        RENL       IN AVERAGED VARIABLES',MYID,RLEMPI(2))
                CALL CHKRLHDL ('        DTL        IN AVERAGED VARIABLES',MYID,RLEMPI(3))
            END IF
            
            DUMMY = (DUMMY1*DBLE(NSTATIS_IO1)-DUMMY2*DBLE(NSTATIS_IO2))/DBLE(NSTATIS_IO1-NSTATIS_IO2)
            NSTATIS_IO = NSTATIS_IO1-NSTATIS_IO2
            DEALLOCATE (DUMMY2)
        END IF
        
        !==================================================
        
        N  = 1
        DO J=1,N2DO(MYID)
            T1xztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !2

        DO J=1,N2DO(MYID)
            D1xztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !3
        
        DO J=1,N2DO(MYID)
            H1xztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !4
        
        DO J=1,N2DO(MYID)
            M1xztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !5
        
        !===============================
        DO J=1,N2DO(MYID)
            T2xztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !6
        
        DO J=1,N2DO(MYID)
            D2xztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !7
        
        DO J=1,N2DO(MYID)
            H2xztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !8
        
        !=============================
        DO J=1,N2DO(MYID)
            DHxztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+1 !9
        
        DO J=1,N2DO(MYID)
            PHxztL_io(J)=DUMMY(J,N)
        ENDDO
        N  = N+0 !9
        
        !===============================
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1
                UHxztL_io(J, L)=DUMMY(J,L+N)
            ENDDO
        ENDDO
        N  = N+L1 ! 12
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1
                GHxztL_io(J, L)=DUMMY(J,L+N)
            ENDDO
        ENDDO
        N  = N+L1 ! 15
        
        L1 = (NDV*(7-NDV))/2+NDV-3 !6
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 16,17,18,19,20,21
                U2DHxztL_io(J, L)=DUMMY(J,L+N)
            ENDDO
        ENDDO
        N  = N+L1 ! 21
        
        !==========================================
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 22,23,24
                DhDL1xztL_io(J, L)=DUMMY(J,L+N)
            ENDDO
        ENDDO
        N  = N+L1 ! 24
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 25,26,27
                DhDLPxztL_io(J, L)=DUMMY(J,L+N)
            ENDDO
        ENDDO
        N  = N+L1 ! 27
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 28,29,30
                DTDLKxztL_io(J, L)=DUMMY(J,L+N)
            ENDDO
        ENDDO
        N  = N+L1 ! 30
        
        !==========================================
        L1 = NDV
        L2 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  !31-33,34-36,37-39
                    DVDL1MxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2 !39
        
        L1 = (NDV-1)*NDV+NDV
        L2 = (NDV-1)*NDV+NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  !40-42,43-45,46-48
                    DVDL2MxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2 !48
        
        L1 = NDV
        L2 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  !49-51,52-54,55-57
                    DVDL1MHxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2 !57
        
        L1 = NDV
        L2 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  !58-60,61-63,64-66
                    DTDLKUxztL_io(J, L, K)=DUMMY(J,M)
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2 !66
        
        !==========================================
        L1 = NDV
        L2 = NDV
        L3 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    DO S=1,L3
                        M = N+(L-1)*L2*L3+(K-1)*L3+S !67-75,76-84,85-93
                        DVDL1MUxztL_io(J, L, K, S)=DUMMY(J,M)
                    END DO
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2*L3 !93
        
        L1 = NDV
        L2 = NDV
        L3 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    DO S=1,L3
                        M = N+(L-1)*L2*L3+(K-1)*L3+S !
                        DTDLKDVDLxztL_io(J, L, K, S)=DUMMY(J,M)
                    END DO
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2*L3 !120
        
        L1 = NDV
        L2 = NDV
        L3 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    DO S=1,L3
                        M = N+(L-1)*L2*L3+(K-1)*L3+S !67-75,76-84,85-93
                        DHDLMDVDLxztL_io(J, L, K, S)=DUMMY(J,M)
                    END DO
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2*L3 !219
        
        IF(N.NE.NSZ) THEN
            WRITE(*,*) 'ERROR in RESTART_AVERAGE_VARS_THERMAL_Xperiodic_IO', N, NSZ
        END IF
        
        DEALLOCATE (DUMMY)
        DEALLOCATE (DUMMY1)
        
        IF(RST_type_flg.EQ.1) THEN
            phyTIME_io = 0.0_WP
            ITERG0_io  = 0
        END IF
        
        if(myid==0) CALL CHKRLHDL('       IO: End of restarting averaged thermal field',myid,phyTIME_IO) 

        RETURN
    
    END SUBROUTINE 
    

!*************************************************************************
    SUBROUTINE RESTART_AVERAGE_SPECTRA
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE
      
        CHARACTER(15) :: PNTIM
        CHARACTER(64) :: FIL1
        logical  :: file_exists
        INTEGER(4)   :: DFLG = 10
        INTEGER(4)   :: N1ML, N2DOL, N3ML
        
        
        REAL(WP)  :: RENL
        REAL(WP)  :: DTL
        INTEGER(4) :: ITEMP
        INTEGER(4) :: IOS

        !IF(ppspectra .NE. 1) RETURN
        IF(MYID .NE. 0)  RETURN
        

        WRITE(PNTIM,'(1ES15.9)') RSTtim_io    
        WRITE(WRT_AVE_FNM_io,'(A,I4.4,A4)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_SPEC.D'
        
        INQUIRE(FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)), EXIST=file_exists) 
        if(.not.file_exists .and. myid==0) CALL ERRHDL('File '//WRT_AVE_FNM_io//' does not exist.',0)
        
        !        !==============OPEN FILE=============================================
        OPEN(DFLG,FILE=TRIM(ADJUSTL(WRT_AVE_FNM_io)),FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IOS)

        READ(DFLG) ITEMP, ITEMP
        READ(DFLG) N1ML,N2DOL,N3ML,ITERG0_io
        READ(DFLG) NSTATIS_io
        READ(DFLG) phyTIME_io,RENL,DTL

        
        !================Velocity =====================
        READ(DFLG) R11X1_xztLa
        READ(DFLG) R22X1_xztLa
        READ(DFLG) R33X1_xztLa
        READ(DFLG) R12X1_xztLa
        READ(DFLG) R13X1_xztLa
        READ(DFLG) R23X1_xztLa
        
        READ(DFLG) R11X3_xztLa
        READ(DFLG) R22X3_xztLa
        READ(DFLG) R33X3_xztLa
        READ(DFLG) R12X3_xztLa
        READ(DFLG) R13X3_xztLa
        READ(DFLG) R23X3_xztLa
        
        READ(DFLG) ENE11T_xztLa
        READ(DFLG) ENE22T_xztLa
        READ(DFLG) ENE33T_xztLa
        READ(DFLG) ENE12T_xztLa
        READ(DFLG) ENE13T_xztLa
        READ(DFLG) ENE23T_xztLa
        
        READ(DFLG) ENE11Z_xztLa
        READ(DFLG) ENE22Z_xztLa
        READ(DFLG) ENE33Z_xztLa
        READ(DFLG) ENE12Z_xztLa
        READ(DFLG) ENE13Z_xztLa
        READ(DFLG) ENE23Z_xztLa
        
        !================Voriticity ====================
        READ(DFLG) V11X1_xztLa
        READ(DFLG) V22X1_xztLa
        READ(DFLG) V33X1_xztLa
        READ(DFLG) V12X1_xztLa
        READ(DFLG) V13X1_xztLa
        READ(DFLG) V23X1_xztLa
        
        READ(DFLG) V11X3_xztLa
        READ(DFLG) V22X3_xztLa
        READ(DFLG) V33X3_xztLa
        READ(DFLG) V12X3_xztLa
        READ(DFLG) V13X3_xztLa
        READ(DFLG) V23X3_xztLa
        
        READ(DFLG) ENV11T_xztLa
        READ(DFLG) ENV22T_xztLa
        READ(DFLG) ENV33T_xztLa
        READ(DFLG) ENV12T_xztLa
        READ(DFLG) ENV13T_xztLa
        READ(DFLG) ENV23T_xztLa
        
        READ(DFLG) ENV11Z_xztLa
        READ(DFLG) ENV22Z_xztLa
        READ(DFLG) ENV33Z_xztLa
        READ(DFLG) ENV12Z_xztLa
        READ(DFLG) ENV13Z_xztLa
        READ(DFLG) ENV23Z_xztLa
        
        
        !===============Voriticity & Velocity=========================
        READ(DFLG) VO11X1_xztLa
        READ(DFLG) VO12X1_xztLa
        READ(DFLG) VO13X1_xztLa
        
        READ(DFLG) VO21X1_xztLa
        READ(DFLG) VO22X1_xztLa
        READ(DFLG) VO23X1_xztLa
        
        READ(DFLG) VO31X1_xztLa
        READ(DFLG) VO32X1_xztLa
        READ(DFLG) VO33X1_xztLa
        
        READ(DFLG) VO11X3_xztLa
        READ(DFLG) VO12X3_xztLa
        READ(DFLG) VO13X3_xztLa
        
        READ(DFLG) VO21X3_xztLa
        READ(DFLG) VO22X3_xztLa
        READ(DFLG) VO23X3_xztLa
        
        READ(DFLG) VO31X3_xztLa
        READ(DFLG) VO32X3_xztLa
        READ(DFLG) VO33X3_xztLa
        
        READ(DFLG) EVO11T_xztLa
        READ(DFLG) EVO12T_xztLa
        READ(DFLG) EVO13T_xztLa
        
        READ(DFLG) EVO21T_xztLa
        READ(DFLG) EVO22T_xztLa
        READ(DFLG) EVO23T_xztLa
        
        READ(DFLG) EVO31T_xztLa
        READ(DFLG) EVO32T_xztLa
        READ(DFLG) EVO33T_xztLa
        
        READ(DFLG) EVO11Z_xztLa
        READ(DFLG) EVO12Z_xztLa
        READ(DFLG) EVO13Z_xztLa
        
        READ(DFLG) EVO21Z_xztLa
        READ(DFLG) EVO22Z_xztLa
        READ(DFLG) EVO23Z_xztLa
        
        READ(DFLG) EVO31Z_xztLa
        READ(DFLG) EVO32Z_xztLa
        READ(DFLG) EVO33Z_xztLa
        
        CLOSE(DFLG)
        if(myid==0) CALL CHKHDL( '       Finish reading '//TRIM(WRT_AVE_FNM_io),myid)
        
        
      
        RETURN
    END SUBROUTINE


!!***********************************************************************************
!    SUBROUTINE RESTART_INSTANT_GLOBAL_io
!        use init_info
!        use mesh_info
!        use flow_info
!        use thermal_info
!        use postprocess_info
!        IMPLICIT NONE
      
!        CHARACTER(15) :: PNTIM
!        CHARACTER(64) :: FIL1
      
!        !INTEGER(4)   :: ITIM
!        INTEGER(4)   :: I, J, K, JJ
!        INTEGER(4)   :: N1ML
!        INTEGER(4)   :: N2DOL
!        INTEGER(4)   :: N3ML
        
!        INTEGER(4)   :: IOS
      
!        REAL(WP),ALLOCATABLE  :: U_F0_io(:,:,:) 
!        REAL(WP),ALLOCATABLE  :: V_F0_io(:,:,:)
!        REAL(WP),ALLOCATABLE  :: W_F0_io(:,:,:)
!        REAL(WP),ALLOCATABLE  :: P_F0_io(:,:,:)
!        REAL(WP),ALLOCATABLE  :: T_F0_io(:,:,:)
!        REAL(WP),ALLOCATABLE  :: D_F0_io(:,:,:)
!        REAL(WP),ALLOCATABLE  :: H_F0_io(:,:,:)
!        REAL(WP),ALLOCATABLE  :: D0_F0_io(:,:,:)
!        REAL(WP)  :: RENL
!        REAL(WP)  :: DTL


!        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
      
!        ALLOCATE( U_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
!        ALLOCATE( W_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
!        ALLOCATE( V_F0_io(0:NCL1_io+1,1:NND2,1:NCL3) )
!        ALLOCATE( P_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
        
        
!        IF(RST_type_flg==0)THEN
!            ALLOCATE( T_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
!            ALLOCATE( D_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
!            ALLOCATE( H_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
!            ALLOCATE( D0_F0_io(0:NCL1_io+1,1:NCL2,1:NCL3) )
!        END IF
            
!        !=======================================U=============================================      
!        FIL1='DNS_inoudomain_RST_'//TRIM(PNTIM)//'_FLOW_INST.bin'
!        if(myid==0) CALL CHKHDL( '       Start reading data from '//FIL1, MYID)
        
!        OPEN(10,FILE=FIL1, FORM='UNFORMATTED', STATUS='OLD',IOSTAT=IOS)
!        IF(IOS/=0) THEN
!            CALL ERRHDL('# CANNOT OPEN FILE: '//TRIM(FIL1),MYID)
!        END IF
        
!        READ(10) N1ML,N2DOL,N3ML,ITERG0_io
!        READ(10) phyTIME_io,RENL,DTL
!        if(myid==0) then
!            call CHKRLHDL   ('           UVW field Time=     ',MYID,phyTIME_io)
!            call CHKINTHDL  ('           UVW field ITERG=    ',MYID,ITERG0_io)
!        end if
        
!        READ(10) U_F0_io
!        READ(10) V_F0_io
!        READ(10) W_F0_io
!        READ(10) P_F0_io
!        CLOSE(10)
!        if(myid==0) CALL CHKHDL( '       Finish reading '//TRIM(FIL1),myid)
!        IF(DABS(phyTIME_io-RSTtim_io)>TSAVE1) CALL CHKHDL('Warning: Restarting from an ealier time step!',myid)

!!====================================================================================
!        DO J=1,N2DO(MYID)
!            JJ=JCL2G(J)
!            DO I=0,NCL1_io+1,1
!                DO K=1,NCL3
!                    G_io(I,J,K,1)    = (U_F0_io(I,JJ,K))
!                    G_io(I,J,K,2)    = (V_F0_io(I,JJ,K))
!                    G_io(I,J,K,3)    = (W_F0_io(I,JJ,K))
!                    PR_io(I,J,K)     = (P_F0_io(I,JJ,K))
!                ENDDO
!            ENDDO
!        ENDDO
      
!        IF (MYID.EQ.NPSLV) THEN
!            DO I=0,NCL1_io+1,1
!                DO K=1,NCL3
!                    !Q_io(I,N2DO(MYID)+1,K,2)=(V_F0_io(I,NND2,K))
!                    G_io(I,N2DO(MYID)+1,K,2)=0.0_WP
!                ENDDO
!            ENDDO
!        ENDIF
      
!        IF (MYID.EQ.0) THEN
!            DO I=0,NCL1_io+1,1
!                DO K=1,NCL3
!                    G_io(I,1,K,2)=0.0_WP
!                    G_io(I,0,K,2)=0.0_WP
!                ENDDO
!            ENDDO
!        ENDIF
!        !====================================================================================
!        if( RST_type_flg==0 ) THEN  !restart thermal field as well
!            FIL1='DNS_inoudomain_RST_'//TRIM(PNTIM)//'_THEL_INST.bin'
!            if(myid==0) CALL CHKHDL( '       Start reading '//FIL1, MYID)
            
!            OPEN(10,FILE=FIL1, FORM='UNFORMATTED', STATUS='OLD',IOSTAT=IOS)
!            IF(IOS/=0) THEN
!                CALL ERRHDL('# CANNOT OPEN FILE: '//TRIM(FIL1),MYID)
!            END IF
            
!            READ(10) N1ML,N2DOL,N3ML,ITERG0_io
!            READ(10) phyTIME_io,RENL,DTL
!            if(myid==0) then
!                call CHKRLHDL   ('           TDH field Time=     ',MYID,phyTIME_io)
!                call CHKINTHDL  ('           TDH field ITERG=    ',MYID,ITERG0_io)
!            end if
!            READ(10) T_F0_io
!            READ(10) D_F0_io
!            READ(10) H_F0_io
!            READ(10) D0_F0_io
!            CLOSE(10)
!            if(myid==0) CALL CHKHDL( '       Finish reading '//TRIM(FIL1),myid)
!            IF(DABS(phyTIME_io-RSTtim_io)>TSAVE1) CALL CHKHDL('Warning: Restarting from an ealier time step!',myid)
          
!            !==============================================================================
!            DO J=1,N2DO(MYID)
!                JJ=JCL2G(J)
!                DO I=0,NCL1_io+1,1
!                    DO K=1,NCL3
!                        TEMPERATURE(I,J,K)= (T_F0_io(I,JJ,K))
!                        DENSITY(I,J,K)   = (D_F0_io(I,JJ,K))
!                        ENTHALPY(I,J,K)  = (H_F0_io(I,JJ,K))
!                        DENSITY0(I,J,K)  = (D0_F0_io(I,JJ,K))
!                    ENDDO
!                ENDDO
!            ENDDO
          
!            !=============BUILD UP RHOH============
!            DO J=1,N2DO(MYID)
!                DO I=0,NCL1_io+1,1
!                    DO K=1,NCL3
!                        RHOH(I,J,K) = DENSITY(I,J,K) * ENTHALPY(I,J,K)
!                    END DO
!                END DO
!            END DO
            
!            !===========build up other variables=======
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML
!            CALL THERM_PROP_UPDATE(IALLDM)
!            CALL INTFC_OUL_THERMAL_io
!            CALL INTFC_MFD_THERMAL_io
!            IF(BCWALLHEAT==isoFluxWall ) CALL BC_WALL_ISOFLUX(IALLDM)
!            CALL VELOCITY_CALC_io
!            !===========build up mass flux=========
!            !CALL Massflux_restart_io
!        END IF
        
!        IF(RST_type_flg==1) THEN
        
!            CALL INIFIELD_THERMAL_io  
!            Q_io = G_io   ! based on constant unit density
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML
        
!            IF(thermlflg==1)THEN
!                CALL SOLVERRK3_ENG_IO(0) 
!                phyTIME_io = 0.0_WP
!                ITERG0_io  = 0
!            END IF
        
!        END IF
        
!!        if(myid==0) then
!!           do j=1,N2DO(myid)
!!              do I=0,NCL1_io+1,1
!!                 do k=1,ncl3
!!                    write(*,'(3I5.1,7ES13.5)') J,I,K, U_F0_io(I,J,K),v_F0_io(I,J,K),w_F0_io(I,J,K),p_F0_io(I,J,K), &
!!                    T_F0_io(I,J,K),d_F0_io(I,J,K),h_F0_io(I,J,K)
!!                 end do
!!              end do
!!            end do
!!        end if

!!        DT=DTL
!!        REN=RENL
      
!        DEALLOCATE (U_F0_io,W_F0_io,P_F0_io,V_F0_io)
!        if( RST_type_flg==0 ) DEALLOCATE (T_F0_io, D_F0_io, H_F0_io)
        
!        CALL MPI_BARRIER(ICOMM,IERROR)
      
!        RETURN
!    END SUBROUTINE
    
!!***********************************************************************************
!    SUBROUTINE RESTART_INSTANT_LOCAL_io
!        use init_info
!        use mesh_info
!        use flow_info
!        use thermal_info
!        use postprocess_info
!        IMPLICIT NONE
      
!        CHARACTER(15) :: PNTIM
!        CHARACTER(15) :: STRIP
!        CHARACTER(64) :: FIL1
      
!        INTEGER(4)   :: I, J, K
!        INTEGER(4)   :: N1ML
!        INTEGER(4)   :: N2DOL
!        INTEGER(4)   :: N3ML
!        INTEGER(4)   :: DFLG = 10
!        INTEGER(4)   :: MYID_TMP, NPTOT_TMP
        
!        INTEGER(4)   :: IOS
      
!        REAL(WP)  :: RENL
!        REAL(WP)  :: DTL


!        !===============CREAT FILE NAME================================
!        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
!        WRITE(STRIP,'(1I3.3)') MYID
            
!        IF(MYID.LT.NPTOT-1) THEN
!            WRITE(FIL1,'(A,I4.4,A4)') &
!                'DNS_inoudomain_INSTANT_RST_'//TRIM(PNTIM)//'_ip'//TRIM(STRIP)//'C.bin'
!        ELSE
!            WRITE(FIL1,'(A,I4.4,A4)') &
!                'DNS_inoudomain_INSTANT_RST_'//TRIM(PNTIM)//'_ip'//TRIM(STRIP)//'E.bin'
!        END IF  
          
!        !==============OPEN FILE=============================================
!        OPEN(DFLG,FILE=TRIM(FIL1),FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IOS)
!        IF(IOS/=0) THEN
!            CALL ERRHDL('# CANNOT OPEN FILE: '//TRIM(FIL1),myid)
!        END IF
!        CALL CHKHDL( '#    READING FILE FROM: '//TRIM(FIL1), MYID)
        
!        !===============READ IN PROCESSOR INFO.====================
!        READ(DFLG) MYID_tmp, NPTOT_tmp
        
!        IF(MYID.NE.MYID_tmp) THEN
!            WRITE(*,*) 'MYID/=IP, WHICH ARE',MYID, MYID_tmp, ' IN FILE '//TRIM(FIL1)
!        END IF
!        IF(NPTOT_TMP.NE.NPTOT) THEN
!            WRITE(*,*) 'NPTOT_TMP/=NPTOT, WHICH ARE',NPTOT_TMP, NPTOT, ' IN FILE '//TRIM(FIL1)
!        END IF
        
!        !=====================READ IN DATA============================
        
!        READ(10) N1ML,N2DOL,N3ML,ITERG0_io
!        READ(10) phyTIME_io,RENL,DTL
        
!        READ(DFLG) G_io
!        READ(DFLG) PR_io
        
!        if( RST_type_flg==0 ) THEN
        
!            READ(DFLG) TEMPERATURE
!            READ(DFLG) DENSITY            
!            READ(DFLG) ENTHALPY            
!            READ(DFLG) DENSITY0

!            !=============BUILD UP RHOH============
!            DO J=1,N2DO(MYID)
!                DO I=0,NCL1_io+1,1
!                    DO K=1,NCL3
!                        RHOH(I,J,K) = DENSITY(I,J,K) * ENTHALPY(I,J,K)
!                    END DO
!                END DO
!            END DO
            
!            !===========build up other variables=======
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML
!            CALL THERM_PROP_UPDATE(IALLDM)
!            CALL INTFC_OUL_THERMAL_io
!            CALL INTFC_MFD_THERMAL_io
!            IF(BCWALLHEAT==isoFluxWall ) CALL BC_WALL_ISOFLUX(IALLDM)
!            CALL VELOCITY_CALC_io
!            !===========build up mass flux=========
!            !CALL Massflux_restart_io
            
!        END IF
!        CLOSE(DFLG)
        
!        IF(RST_type_flg==1) THEN
        
!            CALL INIFIELD_THERMAL_io  
!            Q_io = G_io   ! based on constant unit density
!            CALL BC_TINLET_FLOW
!            CALL BC_TINLET_THERML 
        
!            IF(thermlflg==1)THEN
!                CALL SOLVERRK3_ENG_IO(0) 
!                phyTIME_io = 0.0_WP
!                ITERG0_io  = 0
!            END IF
        
!        END IF
      
!        RETURN
!    END SUBROUTINE

!!*************************************************************************
!    SUBROUTINE RESTART_AVERAGE_GLOBAL_io
!        use init_info
!        use mesh_info
!        use flow_info
!        use thermal_info
!        use postprocess_info
!        IMPLICIT NONE
      
!        CHARACTER(15) :: PNTIM
!        CHARACTER(64) :: FIL1
      
!        INTEGER(4)   :: DFLG = 10
!        INTEGER(4)   :: I, J, JJ
!        INTEGER(4)   :: N1ML
!        INTEGER(4)   :: N2DOL
!        INTEGER(4)   :: N3ML
        
!        INTEGER(4)   :: IOS
      
!        REAL(WP),ALLOCATABLE     :: U1ztL_F0_io(:,:,:)!Ui
!        REAL(WP),ALLOCATABLE     :: G1ztL_F0_io(:,:,:)!Gi
!        REAL(WP),ALLOCATABLE     :: UPztL_F0_io(:,:,:)!UiP
!        REAL(WP),ALLOCATABLE     :: U2ztL_F0_io(:,:,:)!UiUj
!        REAL(WP),ALLOCATABLE     :: UGztL_F0_io(:,:,:)!UiGj
!        REAL(WP),ALLOCATABLE     :: UGUztL_F0_io(:,:,:)!UiGjUk
        
!        REAL(WP),ALLOCATABLE     :: DVDL1ztL_F0_IO(:,:,:,:) ! dUi/dXj
!        REAL(WP),ALLOCATABLE     :: DVDLPztL_F0_IO(:,:,:,:) !PdUi/dXj
!        REAL(WP),ALLOCATABLE     :: DVDL2ztL_F0_IO(:,:,:,:) ! dUi/dXj * dUi/dXj
        
!        REAL(WP),ALLOCATABLE     :: T1ztL_F0_io(:,:)
!        REAL(WP),ALLOCATABLE     :: D1ztL_F0_io(:,:)
!        REAL(WP),ALLOCATABLE     :: H1ztL_F0_io(:,:)
!        REAL(WP),ALLOCATABLE     :: T2ztL_F0_io(:,:)
!        REAL(WP),ALLOCATABLE     :: D2ztL_F0_io(:,:)
!        REAL(WP),ALLOCATABLE     :: H2ztL_F0_io(:,:)
        
!        REAL(WP),ALLOCATABLE     :: DHztL_F0_io(:,:)
        
!        REAL(WP),ALLOCATABLE     :: UHztL_F0_io(:,:,:)
!        REAL(WP),ALLOCATABLE     :: GHztL_F0_io(:,:,:)
        
!        REAL(WP)  :: RENL
!        REAL(WP)  :: DTL


!        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
        
!        ALLOCATE ( U1ztL_F0_io( NCL1_io,NCL2,NDV+1 ) )  
!        ALLOCATE ( G1ztL_F0_io( NCL1_io,NCL2,NDV   ) )
!        ALLOCATE ( UPztL_F0_io( NCL1_io,NCL2,NDV   ) ) 
!        ALLOCATE ( U2ztL_F0_io( NCL1_io,NCL2,NDV*(7-NDV)/2+NDV-3) )
!        ALLOCATE ( UGztL_F0_io( NCL1_io,NCL2,NDV*(7-NDV)/2+NDV-3) )
!        ALLOCATE ( UGUztL_F0_io(NCL1_io,NCL2,NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8) )
    
!        ALLOCATE ( DVDL1ztL_F0_io( NCL1_io,NCL2, NDV, NDV  ) )
!        ALLOCATE ( DVDLPztL_F0_io( NCL1_io,NCL2, NDV, NDV  ) ) 
!        ALLOCATE ( DVDL2ztL_F0_io( NCL1_io,NCL2, NDV*(7-NDV)/2+NDV-3, NDV  ) )
                    
!        IF(RST_type_flg==0)THEN            
!            ALLOCATE ( T1ztL_F0_io( NCL1_io,NCL2 ) )
!            ALLOCATE ( D1ztL_F0_io( NCL1_io,NCL2 ) )
!            ALLOCATE ( H1ztL_F0_io( NCL1_io,NCL2 ) )
                 
!            ALLOCATE ( T2ztL_F0_io( NCL1_io,NCL2 ) )
!            ALLOCATE ( D2ztL_F0_io( NCL1_io,NCL2 ) ) 
!            ALLOCATE ( H2ztL_F0_io( NCL1_io,NCL2 ) )
                
!            ALLOCATE ( UHztL_F0_io( NCL1_io,NCL2,NDV ) )
!            ALLOCATE ( GHztL_F0_io( NCL1_io,NCL2,NDV ) ) 
!        END IF
!        !=======================================U=============================================      
!        FIL1='DNS_inoudomain_RST_'//TRIM(PNTIM)//'_FLOW_AVEG.bin'
!        if(myid==0) CALL CHKHDL( '       Start reading '//FIL1, MYID)
        
!        OPEN(DFLG,FILE=FIL1, FORM='UNFORMATTED', STATUS='OLD',IOSTAT=IOS)
!        IF(IOS/=0) THEN
!            CALL ERRHDL('# CANNOT OPEN FILE: '//TRIM(FIL1),MYID)
!        END IF
        
!        READ(DFLG) N1ML,N2DOL,N3ML,ITERG0_io
!        READ(DFLG) NSTATIS_io
!        READ(DFLG) phyTIME_io,RENL,DTL
!        if(myid==0) then
!            call CHKRLHDL   ('           FLOW_AVEG field Time=     ',MYID,phyTIME_io)
!            call CHKINTHDL  ('           FLOW_AVEG field ITERG=    ',MYID,ITERG0_io)
!            call CHKINTHDL  ('           FLOW_AVEG field INSTATIS_io=    ',MYID,NSTATIS_io)
!        end if
        
!        READ(DFLG) U1ztL_F0_io
!        READ(DFLG) G1ztL_F0_io
!        READ(DFLG) UPztL_F0_io
        
!        READ(DFLG) U2ztL_F0_IO
!        READ(DFLG) UGztL_F0_IO
!        READ(DFLG) UGUztL_F0_IO
        
!        READ(DFLG) DVDL1ztL_F0_io
!        READ(DFLG) DVDLPztL_F0_io
!        READ(DFLG) DVDL2ztL_F0_io
        
!        CLOSE(DFLG)
!        if(myid==0) CALL CHKHDL( '       Finish reading '//TRIM(FIL1),myid)
        
!        DO J=1,N2DO(MYID)
!            JJ=JCL2G(J)
!            DO I=1,NCL1_io
!                U1ztL_io( I,J,: ) = U1ztL_F0_io( I,JJ,: ) 
!                G1ztL_io( I,J,: ) = G1ztL_F0_io( I,JJ,: ) 
!                UPztL_iO( I,J,: ) = UPztL_F0_iO( I,JJ,: )
                
!                U2ztL_io( I,J,: ) = U2ztL_F0_io( I,JJ,: )
!                UGztL_io( I,J,: ) = UGztL_F0_io( I,JJ,: )
                
!                UGUztL_io( I,J,: )= UGUztL_F0_io( I,JJ,: )
                
!                DVDL1ztL_io( I,J,:,: ) = DVDL1ztL_F0_io( I,JJ,:,: ) 
!                DVDLPztL_io( I,J,:,: ) = DVDLPztL_F0_io( I,JJ,:,: ) 
!                DVDL2ztL_io( I,J,:,: ) = DVDL2ztL_F0_io( I,JJ,:,: ) 
!            ENDDO
!        ENDDO
        
!        IF(RST_type_flg==0) THEN
        
!            FIL1='DNS_inoudomain_RST_'//TRIM(PNTIM)//'_THEL_AVEG.bin'
!            if(myid==0) CALL CHKHDL( '       Start reading data from '//FIL1, MYID)
            
!            OPEN(DFLG,FILE=FIL1, FORM='UNFORMATTED', STATUS='OLD',IOSTAT=IOS)
!            IF(IOS/=0) THEN
!                CALL ERRHDL('# CANNOT OPEN FILE: '//TRIM(FIL1),MYID)
!            END IF
            
!            READ(DFLG) N1ML,N2DOL,N3ML,ITERG0_io
!            READ(DFLG) NSTATIS_io
!            READ(DFLG) phyTIME_io,RENL,DTL
!            if(myid==0) then
!                call CHKRLHDL   ('           THEL_AVEG field Time=     ',MYID,phyTIME_io)
!                call CHKINTHDL  ('           THEL_AVEG field ITERG=    ',MYID,ITERG0_io)
!            end if
        
!            READ(DFLG) T1ztL_F0_io
!            READ(DFLG) D1ztL_F0_io
!            READ(DFLG) H1ztL_F0_io
            
!            READ(DFLG) T2ztL_F0_io
!            READ(DFLG) D2ztL_F0_io
!            READ(DFLG) H2ztL_F0_io
            
!            READ(DFLG) DHztL_F0_io
            
!            READ(DFLG) UHztL_F0_io
!            READ(DFLG) GHztL_F0_io
            
!            CLOSE(DFLG)
!            if(myid==0) CALL CHKHDL( '       Finish reading '//TRIM(FIL1),myid)
            
!            DO J=1,N2DO(MYID)
!                JJ=JCL2G(J)
!                DO I=1,NCL1_io
!                    T1ztL_io( I,J ) = T1ztL_F0_io( I,JJ )
!                    D1ztL_io( I,J ) = D1ztL_F0_io( I,JJ ) 
!                    H1ztL_io( I,J ) = H1ztL_F0_io( I,JJ ) 
                    
!                    T2ztL_io( I,J ) = T2ztL_F0_io( I,JJ ) 
!                    D2ztL_io( I,J ) = D2ztL_F0_io( I,JJ ) 
!                    H2ztL_io( I,J ) = H2ztL_F0_io( I,JJ ) 
                    
!                    UHztL_io( I,J,: ) = UHztL_F0_io( I,JJ,: ) 
!                    GHztL_io( I,J,: ) = GHztL_F0_io( I,JJ,: ) 
!                END DO
!            END DO
            
!        END IF
        
      
!        DEALLOCATE (U1ztL_F0_io, G1ztL_F0_io, UPztL_F0_iO, U2ztL_F0_io, UGztL_F0_io,&
!        UGUztL_F0_io, DVDL1ztL_F0_io, DVDLPztL_F0_io, DVDL2ztL_F0_io)
!        if( RST_type_flg==0 ) DEALLOCATE (T1ztL_F0_io, D1ztL_F0_io, H1ztL_F0_io, &
!        T2ztL_F0_io, D2ztL_F0_io, H2ztL_F0_io, UHztL_F0_io, GHztL_F0_io)
        
!        CALL MPI_BARRIER(ICOMM,IERROR)
      
!        RETURN
!    END SUBROUTINE
!!********************************************************************************
!    SUBROUTINE RESTART_AVERAGE_LOCAL_io
!        use init_info
!        use mesh_info
!        use flow_info
!        use thermal_info
!        use postprocess_info
!        IMPLICIT NONE
      
!        CHARACTER(15) :: PNTIM
!        CHARACTER(15) :: STRIP
!        CHARACTER(64) :: FIL1
      
!        INTEGER(4)   :: DFLG = 10
!        INTEGER(4)   :: N1ML
!        INTEGER(4)   :: N2DOL
!        INTEGER(4)   :: N3ML
!        INTEGER(4)   :: MYID_TMP, NPTOT_TMP
        
!        INTEGER(4)   :: IOS
      
!        REAL(WP)  :: RENL
!        REAL(WP)  :: DTL


!        !================CREAT FILE NAME====================================
!        WRITE(PNTIM,'(1ES15.9)') RSTtim_io
!        WRITE(STRIP,'(1I3.3)') MYID
!        IF(MYID.LT.NPTOT-1) THEN
!            WRITE(FIL1,'(A,I4.4,A4)') &
!                'DNS_inoudomain_AVERAGE_RST_'//TRIM(PNTIM)//'_ip'//TRIM(STRIP)//'C.bin'
!        ELSE
!            WRITE(FIL1,'(A,I4.4,A4)') &
!                'DNS_inoudomain_AVERAGE_RST_'//TRIM(PNTIM)//'_ip'//TRIM(STRIP)//'E.bin'
!        END IF
        
!        !==============OPEN FILE=============================================
!        OPEN(DFLG,FILE=TRIM(FIL1),FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IOS)
!        IF(IOS/=0) THEN
!            CALL ERRHDL('# CANNOT OPEN FILE: '//TRIM(FIL1),myid)
!        END IF
!        CALL CHKHDL( '#    READING FILE FROM: '//TRIM(FIL1), MYID)
        
!        !==============read in processor info.================================
!        READ(DFLG) MYID_tmp, NPTOT_TMP
        
!        IF(MYID.NE.MYID_tmp) THEN
!            WRITE(*,*) 'MYID/=IP, WHICH ARE',MYID, MYID_tmp, ' IN FILE '//TRIM(FIL1)
!        END IF
!        IF(NPTOT_TMP.NE.NPTOT) THEN
!            WRITE(*,*) 'NPTOT_TMP/=NPTOT, WHICH ARE',NPTOT_TMP, NPTOT, ' IN FILE '//TRIM(FIL1)
!        END IF
        
!        !=================read in data================================
!        READ(DFLG) N1ML,N2DOL,N3ML,ITERG0_io
!        READ(DFLG) NSTATIS_io
!        READ(DFLG) phyTIME_io,RENL,DTL
                
!        READ(DFLG) U1ztL_io
!        READ(DFLG) G1ztL_io
!        READ(DFLG) UPztL_io
        
!        READ(DFLG) U2ztL_IO
!        READ(DFLG) UGztL_IO
        
!        READ(DFLG) UGUztL_IO
        
!        READ(DFLG) DVDL1ztL_io
!        READ(DFLG) DVDLPztL_io
!        READ(DFLG) DVDL2ztL_io
        
!        IF(RST_type_flg==0) THEN
        
!            READ(DFLG) T1ztL_io
!            READ(DFLG) D1ztL_io
!            READ(DFLG) H1ztL_io
            
!            READ(DFLG) T2ztL_io
!            READ(DFLG) D2ztL_io
!            READ(DFLG) H2ztL_io
            
!            READ(DFLG) DHztL_io
            
!            READ(DFLG) UHztL_io
!            READ(DFLG) GHztL_io
            
!        END IF
        
!        CLOSE(DFLG)
      
!        RETURN
!    END SUBROUTINE
    

