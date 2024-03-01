    SUBROUTINE WRT_INSTANT_VARS_io
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE

        CHARACTER(15) :: PNTIM
        INTEGER(4)    :: DFLG(4), IVEL, I, J, K
        INTEGER(4)    :: IERR,SIZES_ARRAY(3)
        INTEGER(4)    :: NEWTYPE,SUBSIZES(3),STARTS(3)
        INTEGER(4)    :: IRLSIZE,INTMPI(4),INTBYTE,DBLBYTE,INISIZE
        INTEGER(4)    ::  BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND)::  OFFSET
        REAL(WP)               :: RLEMPI(3)
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:)
        character(64) :: FLNAME
        
        
        ALLOCATE(DUMMY(1:NCL1E,1:N2DO(MYID),1:NCL3))
        DUMMY = 0.0_wp
        IF(TGFLOWFLG .and. IOFLOWFLG) THEN
            FLNAME ='DNS_perioz_INSTANT_T'
        ELSE 
            FLNAME ='DNS_perixz_INSTANT_T'
        END IF
        
        DO IVEL=1,NDV+1
            !==================GENERATE FILE NAMES===============================================
            DFLG(IVEL) = 100+IVEL
            WRITE(PNTIM,'(1ES15.9)') phyTIME
            IF(IVEL.EQ.1) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_U.D'
            IF(IVEL.EQ.2) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_V.D'
            IF(IVEL.EQ.3) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_W.D'
            IF(IVEL.EQ.4) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_P.D'
            
            !========================RE-ASSIGN VARIABLES=======================================
            DO I=1,NCL1E
                DO J=1,N2DO(MYID)
                    DO K=1,NCL3
                        IF(IVEL.EQ.1) DUMMY(I,J,K)=G_io(I,J,K,1) !U 
                        IF(IVEL.EQ.2) DUMMY(I,J,K)=G_io(I,J,K,2) !V
                        IF(IVEL.EQ.3) DUMMY(I,J,K)=G_io(I,J,K,3) !W
                        IF(IVEL.EQ.4) DUMMY(I,J,K)=PR_io(I,J,K)  !P
                    ENDDO
                ENDDO
            ENDDO

            SIZES_ARRAY(1)=NCL1E
            SIZES_ARRAY(2)=NCL2
            SIZES_ARRAY(3)=NCL3
            
            SUBSIZES(1)=SIZES_ARRAY(1)
            SUBSIZES(2)=N2DO(MYID)
            SUBSIZES(3)=SIZES_ARRAY(3)
            
            STARTS(1)=0
            STARTS(2)=JCL2G(1)-1  !SUBSIZES(2) * MYID
            STARTS(3)=0
            BUFSIZE=SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)
            
            CALL MPI_TYPE_CREATE_SUBARRAY(3,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
            CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
            CALL MPI_FILE_OPEN(ICOMM,WRT_RST_FNM_io(IVEL),MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,MYFILE,IERR) 
                                

            INTMPI(1)=NCL1E
            INTMPI(2)=NCL2
            INTMPI(3)=NCL3
            INTMPI(4)=ITERG
            
            RLEMPI(1)=phyTIME
            RLEMPI(2)=REN 
            RLEMPI(3)=DT   
            
            INISIZE=4
            IRLSIZE=3

            OFFSET=0_MPI_OFFSET_KIND
            CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
            CALL MPI_FILE_WRITE(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
            CALL MPI_FILE_WRITE(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
            CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
            CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
            
            OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
            
            CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
            CALL MPI_FILE_WRITE_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
            CALL MPI_FILE_CLOSE(MYFILE,IERR)
            CALL MPI_TYPE_FREE(NEWTYPE,IERR)
            CALL MPI_BARRIER(ICOMM,IERROR)
        ENDDO
        
        DUMMY = 0.0_wp
        IF(thermlflg==1) THEN
            DO IVEL=1,3
                !==================GENERATE FILE NAMES===============================================
                DFLG(IVEL) = 120+IVEL
                WRITE(PNTIM,'(1ES15.9)') phyTIME
                IF(IVEL.EQ.1) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_T.D'
                IF(IVEL.EQ.2) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_D.D'
                IF(IVEL.EQ.3) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_E.D'
                !IF(IVEL.EQ.4) WRITE(WRT_RST_FNM_io(IVEL),'(A)') TRIM(filepath1)//TRIM(FLNAME)//TRIM(PNTIM)//'_D0.D'
                
                !========================RE-ASSIGN VARIABLES=======================================
                DO I=1,NCL1E
                    DO J=1,N2DO(MYID)
                        DO K=1,NCL3
                            IF(IVEL.EQ.1) DUMMY(I,J,K)=TEMPERATURE(I, J, K) 
                            IF(IVEL.EQ.2) DUMMY(I,J,K)=DENSITY(I, J, K) 
                            IF(IVEL.EQ.3) DUMMY(I,J,K)=ENTHALPY(I, J, K) 
                            !IF(IVEL.EQ.4) DUMMY(I,J,K)=DENSITY0(I, J, K) 
                        ENDDO
                    ENDDO
                ENDDO
    
                SIZES_ARRAY(1)=NCL1E
                SIZES_ARRAY(2)=NCL2
                SIZES_ARRAY(3)=NCL3
                
                SUBSIZES(1)=SIZES_ARRAY(1)
                SUBSIZES(2)=N2DO(MYID)
                SUBSIZES(3)=SIZES_ARRAY(3)
                
                STARTS(1)=0
                STARTS(2)=JCL2G(1)-1  !SUBSIZES(2) * MYID
                STARTS(3)=0
                
                BUFSIZE=SUBSIZES(1) * SUBSIZES(2) * SUBSIZES(3)
                
                CALL MPI_TYPE_CREATE_SUBARRAY(3,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, &
                                                MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
                CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
                CALL MPI_FILE_OPEN(ICOMM,WRT_RST_FNM_io(IVEL),MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,MYFILE,IERR) 
                                    
    
                INTMPI(1)=NCL1E
                INTMPI(2)=NCL2
                INTMPI(3)=NCL3
                INTMPI(4)=ITERG
                
                RLEMPI(1)=phyTIME
                RLEMPI(2)=REN 
                RLEMPI(3)=DT   
                
                INISIZE=4
                IRLSIZE=3
    
                OFFSET=0_MPI_OFFSET_KIND
                CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
                CALL MPI_FILE_WRITE(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
                CALL MPI_FILE_WRITE(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
                CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
                CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
                
                OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
                
                CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
                CALL MPI_FILE_WRITE_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
                CALL MPI_FILE_CLOSE(MYFILE,IERR)
                CALL MPI_TYPE_FREE(NEWTYPE,IERR)
                CALL MPI_BARRIER(ICOMM,IERROR)
            ENDDO

        END IF
        DEALLOCATE (DUMMY)
        
        
        
                
        RETURN
    END SUBROUTINE


!======================================================================================================================
    SUBROUTINE WRT_AVERAGE_VARS_nonXperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE

        CHARACTER(15) :: PNTIM
        INTEGER(4)    :: L1, L2
        INTEGER(4)    :: DFLG, I, J, L, M, N, H, P
        INTEGER(4)    :: IERR,SIZES_ARRAY(4)
        INTEGER(4)    :: NEWTYPE,SUBSIZES(4),STARTS(4)
        INTEGER(4)    :: IRLSIZE,INTMPI(4),INTBYTE,DBLBYTE,INISIZE
        INTEGER(4)    ::  BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND)::  OFFSET
        REAL(WP)               :: RLEMPI(3)
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:,:)

        ALLOCATE( DUMMY(1:NCL1_io, 1:N2DO(MYID), 1:(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8), 1:21 ) )
        DUMMY = 0.0_wp
        
        DFLG = 100
        WRITE(PNTIM,'(1ES15.9)') phyTIME
        WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perioz_AVERAGD_T'//TRIM(PNTIM)//'_FLOW.D'

        L1 = NDV+1
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L,1)=U1ztL_io(I,J,L)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L,2)=G1ztL_io(I,J,L)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L,3)=UPztL_io(I,J,L)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV*(7-NDV)/2+NDV-3
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L,4)=U2ztL_io(I,J,L)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV*(7-NDV)/2+NDV-3
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L,5)=UGztL_io(I,J,L)
                ENDDO
            ENDDO
        END DO

        L1 = NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L,6)=UGUztL_io(I,J,L)
                ENDDO
            ENDDO
        ENDDO
        
        
        L1 = NDV
        L2 = NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DO M=1,L2
                        DUMMY(I,J,L,M+6)=DVDL1ztL_io(I, J, L, M)
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
                        DUMMY(I,J,L,M+9)=DVDLPztL_io(I, J, L, M)
                    END DO
                ENDDO
            ENDDO
        END DO
        
        
!        DO M=1,NDV
!            DO N=1,NDV
!                DO H=1,NDV
!                    DO P=1, NDV
!                        L1=(M-1)*3+H
!                        L2=(N-1)*3+P
!                        DO J=1,N2DO(MYID)
!                            DUMMY(I,J,L1,L2+12)=DVDL2ztL_io(I, J, L1, L2)
!                        END DO
!                    END DO
!                END DO
!            END DO
!        END DO
        
        DO M=1,NDV
            DO H=1, NDV
                L1   = (M-1)*NDV+H
                DO N=1, NDV
                    DO P=1, NDV
                        L2=(N-1)*NDV+P
                        DO J=1,N2DO(MYID)
                            DUMMY(I,J,L1,L2+12)=DVDL2ztL_io(I, J, L1, L2)
                        END DO
                    END DO
                END DO
            END DO
        END DO
        
        
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
        
        CALL MPI_TYPE_CREATE_SUBARRAY(4,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io,MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,MYFILE,IERR)
         
        INISIZE=4
        IRLSIZE=3    
                        
        INTMPI(1)=NCL1_io
        INTMPI(2)=NCL2
        INTMPI(3)=ITERG
        INTMPI(4)=NSTATIS_IO
        
        RLEMPI(1)=PHYTIME
        RLEMPI(2)=REN 
        RLEMPI(3)=DT   
        

        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_WRITE(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        DEALLOCATE (DUMMY)
        
        IF(thermlflg==1) call WRT_AVERAGE_VARS_THERMAL_nonXperiodic_IO
       
        RETURN
    END SUBROUTINE



!***********************************************************************************************************
    SUBROUTINE WRT_AVERAGE_VARS_THERMAL_nonXperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE

        CHARACTER(15) :: PNTIM
        INTEGER(4)    :: L1
        INTEGER(4)    :: DFLG, I, J, L, N
        INTEGER(4)    :: IERR,SIZES_ARRAY(3)
        INTEGER(4)    :: NEWTYPE,SUBSIZES(3),STARTS(3)
        INTEGER(4)    :: IRLSIZE,INTMPI(4),INTBYTE,DBLBYTE,INISIZE
        INTEGER(4)    ::  BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND)::  OFFSET
        REAL(WP)               :: RLEMPI(3)
        REAL(WP),ALLOCATABLE :: DUMMY(:,:,:)

        
        ALLOCATE( DUMMY(1:NCL1_io, 1:N2DO(MYID), 7+2*NDV ) )
        DUMMY = 0.0_wp
        
        DFLG = 102
        WRITE(PNTIM,'(1ES15.9)') phyTIME
        WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perioz_AVERAGD_T'//TRIM(PNTIM)//'_THEL.D'


        N  = 1
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DUMMY(I,J,N)=T1ztL_io(I,J)
            ENDDO
        END DO
        

        N  = 2
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DUMMY(I,J,N)=D1ztL_io(I,J)
            ENDDO
        END DO
        
        N  = 3
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DUMMY(I,J,N)=H1ztL_io(I,J)
            ENDDO
        END DO
        
        N  = 4
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DUMMY(I,J,N)=T2ztL_io(I,J)
            ENDDO
        END DO
        

        N  = 5
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DUMMY(I,J,N)=D2ztL_io(I,J)
            ENDDO
        END DO
        
        N  = 6
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DUMMY(I,J,N)=H2ztL_io(I,J)
            ENDDO
        END DO
        
        N  = 7
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DUMMY(I,J,N)=DHztL_io(I,J)
            ENDDO
        END DO
        
        L1 = NDV
        N  = 7
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L+N)=UHztL_io(I, J, L)
                ENDDO
            ENDDO
        END DO
        
        L1 = NDV
        N  = 7+NDV
        DO I=1, NCL1_io
            DO J=1,N2DO(MYID)
                DO L=1,L1
                    DUMMY(I,J,L+N)=GHztL_io(I, J, L)
                ENDDO
            ENDDO
        END DO
        
        
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
        
        CALL MPI_TYPE_CREATE_SUBARRAY(3,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io,MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,MYFILE,IERR) 
        
        INISIZE=4
        IRLSIZE=3
        
        INTMPI(1)=NCL1_io
        INTMPI(2)=NCL2
        INTMPI(3)=ITERG
        INTMPI(4)=NSTATIS_io
        
        RLEMPI(1)=PHYTIME
        RLEMPI(2)=REN 
        RLEMPI(3)=DT   
        

        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_WRITE(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        CALL MPI_BARRIER(ICOMM,IERROR)
    

        DEALLOCATE (DUMMY)
        
                
        RETURN
    END SUBROUTINE


!*******************************************************************************************************************
!======================================================================================================================
    SUBROUTINE WRT_AVERAGE_VARS_Xperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE

        CHARACTER(15) :: PNTIM
        INTEGER(4)    :: L1, L2, NSZ
        INTEGER(4)    :: DFLG, J, L, M, N, K
        INTEGER(4)    :: IERR,NEWTYPE
        INTEGER(4)    :: NARRAY
        INTEGER(4),ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),STARTS(:)
        INTEGER(4),ALLOCATABLE :: INTMPI(:)
        REAL(WP),ALLOCATABLE :: RLEMPI(:)
        INTEGER(4)    :: IRLSIZE,INTBYTE,DBLBYTE,INISIZE
        INTEGER(4)    ::  BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND)::  OFFSET
        REAL(WP),ALLOCATABLE :: DUMMY(:,:)
        
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
        
              
        ALLOCATE( DUMMY(1:N2DO(MYID), NSZ ) )
        DUMMY = 0.0_wp
        !============================
        
        SIZES_ARRAY(1)=NCL2
        SIZES_ARRAY(2)=NSZ
        
        SUBSIZES(1)=N2DO(MYID)
        SUBSIZES(2)=SIZES_ARRAY(2)
        
        STARTS(1)=JCL2G(1)-1  !SUBSIZES(2) * MYID
        STARTS(2)=0
        BUFSIZE  =SUBSIZES(1) * SUBSIZES(2)
        
        
        INTMPI(1)=NCL2
        INTMPI(2)=NSZ
        INTMPI(3)=ITERG
        INTMPI(4)=NSTATIS_io
        
        RLEMPI(1)=PHYTIME
        RLEMPI(2)=REN 
        RLEMPI(3)=DT   
        

        !==========================================
        N=0
        L1 = NDV+1
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=U1xztL_io(J,L) !1-4
            ENDDO
        ENDDO
        N = N + L1 !4
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in U1xztL_io',myid, L1, N) 
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=G1xztL_io(J,L)!5-7
            ENDDO
        ENDDO
        N = N + L1 !7
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in G1xztL_io',myid, L1, N) 
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=UPxztL_io(J,L) !8-10
            ENDDO
        ENDDO
        N = N + L1 !10
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in UPxztL_io',myid, L1, N) 
        
        !==========================================
        L1 = NDV*(7-NDV)/2+NDV-3
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=U2xztL_io(J,L) !11-16
            ENDDO
        ENDDO
        N = N + L1 !16
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in U2xztL_io',myid, L1, N) 
        
        L1 = NDV*(7-NDV)/2+NDV-3
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=UGxztL_io(J,L)
            ENDDO
        ENDDO
        N = N + L1 !22
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in UGxztL_io',myid, L1, N) 
        
        
        L1 = NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8 !10
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=UGUxztL_io(J,L)
            ENDDO
        ENDDO
        N = N + L1 !32
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in UGUxztL_io',myid, L1, N) 
        
!!        ! Commented below is no U3=================================
        L1 = NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8 !10
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=U3xztL_io(J,L)
            ENDDO
        ENDDO
        N = N + L1 !32
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in U3xztL_io',myid, L1, N) 
        
        !==========================================
        L1 = NDV
        L2 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=DVDL1xztL_io(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in DVDL1xztL_io',myid, L1, N) 
        
        
        L1 = NDV
        L2 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=DVDLPxztL_io(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in DVDLPxztL_io',myid, L1, N) 
        
        L1= (NDV-1)*3+NDV
        L2= (NDV-1)*3+NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=DVDL2xztL_io(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in DVDL2xztL_io',myid, L1, N) 
        
        !=======for quadrant================
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=QUADUVxztL_IO(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in QUADUVxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=QUADVzxztL_IO(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in QUADVzxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=QUADTKxztL_IO(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in QUADTKxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=QUADDRxztL_IO(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in QUADDRxztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=QUADDUV1xztL_io(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in QUADDUV1xztL_io',myid, L1, N) 
        
        L1= 4
        L2= QUADHN
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 
                DO K=1,L2
                    M = N+(L-1)*L2+K  
                    DUMMY(J,M)=QUADDUV2xztL_io(J, L, K)
                END DO
            ENDDO
        ENDDO
        N = N + L1*L2
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in QUADDUV2xztL_io',myid, L1, N) 
        
        !=============with driven force=========================
        L1 = NDV+1
        DO J=1,N2DO(MYID)
            DO L=1,L1
                DUMMY(J,N+L)=FUxztL_IO(J,L) !1-4
            ENDDO
        ENDDO
        N = N + L1 !4
        !IF(MYID==0) CALL CHK2INTHDL('        Writing in FUxztL_io',myid, L1, N) 
        
!        !=======for octant===\rho ============
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTDUVxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTDUVxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTDVzxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTDVzxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTDTKxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTDTKxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTDDRxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTDDRxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTDDUV1xztL_io(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTDDUV1xztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTDDUV2xztL_io(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTDDUV2xztL_io',myid, L1, N) 
        
!        !=======for octant===T ============
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTTUVxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTTUVxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTTVzxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTTVzxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTTTKxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTTTKxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTTDRxztL_IO(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTTDRxztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTTDUV1xztL_io(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTTDUV1xztL_io',myid, L1, N) 
        
!        L1= 8
!        L2= QUADHN
!        DO J=1,N2DO(MYID)
!            DO L=1,L1 ! 
!                DO K=1,L2
!                    M = N+(L-1)*L2+K  
!                    DUMMY(J,M)=OCTTDUV2xztL_io(J, L, K)
!                END DO
!            ENDDO
!        ENDDO
!        N = N + L1*L2
!        !IF(MYID==0) CALL CHK2INTHDL('        Writing in OCTTDUV2xztL_io',myid, L1, N) 
        
        
        IF(N.NE.NSZ) THEN
            WRITE(*,*) 'ERROR in WRT_AVERAGE_VARS_Xperiodic_IO', N, NSZ
        END IF
        
        !!==========================================
            
        DFLG = 100
        WRITE(PNTIM,'(1ES15.9)') phyTIME
        WRITE(WRT_AVE_FNM_IO,'(A)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_FLOW.D'
        if(myid==0) CALL CHKHDL('    IO: WRITING averaged flow field '//TRIM(WRT_AVE_FNM_io),myid) 
        
        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io,MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,MYFILE,IERR) 

        IF(MYID==0) THEN
            CALL CHKINTHDL('        J IN AVERAGED VARIABLES (J,L)   ',MYID,INTMPI(1))
            CALL CHKINTHDL('        L IN AVERAGED VARIABLES (J,L)   ',MYID,INTMPI(2))
            CALL CHKINTHDL('        ITERG0_IO IN AVERAGED VARIABLES ',MYID,INTMPI(3))
            CALL CHKINTHDL('        NSTATIS_IO IN AVERAGED VARIABLES',MYID,INTMPI(4))
            
            CALL CHKRLHDL ('        phyTIME_TG IN AVERAGED VARIABLES',MYID,RLEMPI(1))
            CALL CHKRLHDL ('        RENL       IN AVERAGED VARIABLES',MYID,RLEMPI(2))
            CALL CHKRLHDL ('        DTL        IN AVERAGED VARIABLES',MYID,RLEMPI(3))
        END IF

        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_WRITE(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        
        DEALLOCATE (DUMMY)
        
        DEALLOCATE ( SIZES_ARRAY )
        DEALLOCATE ( SUBSIZES    )
        DEALLOCATE ( STARTS      )
        DEALLOCATE ( INTMPI      )
        DEALLOCATE ( RLEMPI      )
        
        IF(thermlflg==1) call WRT_AVERAGE_VARS_THERMAL_Xperiodic_IO
        CALL WRT_AVERAGE_SPECTRA
                
        RETURN
    END SUBROUTINE


!***********************************************************************************************************
    SUBROUTINE WRT_AVERAGE_VARS_THERMAL_Xperiodic_IO
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE POSTPROCESS_INFO
        use thermal_info
        IMPLICIT NONE

        CHARACTER(15) :: PNTIM
        INTEGER(4)    :: L1,L2,L3
        INTEGER(4)    :: DFLG, J, L, N, K,M, NSZ, S
        INTEGER(4)    :: IERR
        INTEGER(4)    :: NEWTYPE
        INTEGER(4)    :: NARRAY
        INTEGER(4),ALLOCATABLE :: SIZES_ARRAY(:),SUBSIZES(:),STARTS(:)
        INTEGER(4),ALLOCATABLE :: INTMPI(:)
        REAL(WP),ALLOCATABLE :: RLEMPI(:)
        INTEGER(4)    :: IRLSIZE,INTBYTE,DBLBYTE,INISIZE
        INTEGER(4)    ::  BUFSIZE,MYFILE,STATUS(MPI_STATUS_SIZE)
        INTEGER(KIND=MPI_OFFSET_KIND)::  OFFSET
        REAL(WP),ALLOCATABLE :: DUMMY(:,:)
        
        
        !============================
        NARRAY=2
        ALLOCATE ( SIZES_ARRAY(NARRAY) )
        ALLOCATE ( SUBSIZES   (NARRAY) )
        ALLOCATE ( STARTS     (NARRAY) )
        INISIZE=4
        IRLSIZE=3
        ALLOCATE ( INTMPI(INISIZE)       )
        ALLOCATE ( RLEMPI(IRLSIZE)       )


        NSZ = 9+5*NDV+((NDV*(7-NDV))/2+NDV-3)+3*NDV*NDV+3*NDV*NDV*NDV+((NDV-1)*NDV+NDV)*((NDV-1)*NDV+NDV)
        ALLOCATE( DUMMY(1:N2DO(MYID), NSZ ) )
        DUMMY = 0.0_wp
        
        !================================
        SIZES_ARRAY(1)=NCL2
        SIZES_ARRAY(2)=NSZ
        
        SUBSIZES(1)=N2DO(MYID)
        SUBSIZES(2)=SIZES_ARRAY(2)
        
        STARTS(1)=JCL2G(1)-1  !SUBSIZES(2) * MYID
        STARTS(2)=0
        
        BUFSIZE=SUBSIZES(1) * SUBSIZES(2)
        
        INTMPI(1)=NCL2
        INTMPI(2)=NSZ
        INTMPI(3)=ITERG
        INTMPI(4)=NSTATIS_io
        
        RLEMPI(1)=PHYTIME
        RLEMPI(2)=REN 
        RLEMPI(3)=DT   

        !==========================================
        N  = 1
        DO J=1,N2DO(MYID)
            DUMMY(J,N)=T1xztL_io(J)
        ENDDO
        N  = N+1 !2

        DO J=1,N2DO(MYID)
            DUMMY(J,N)=D1xztL_io(J)
        ENDDO
        N  = N+1 !3
        
        DO J=1,N2DO(MYID)
            DUMMY(J,N)=H1xztL_io(J)
        ENDDO
        N  = N+1 !4
        
        DO J=1,N2DO(MYID)
            DUMMY(J,N)=M1xztL_io(J)
        ENDDO
        N  = N+1 !5
        
        !==========================================
        DO J=1,N2DO(MYID)
            DUMMY(J,N)=T2xztL_io(J)
        ENDDO
        N  = N+1 !6

        DO J=1,N2DO(MYID)
            DUMMY(J,N)=D2xztL_io(J)
        ENDDO
        N  = N+1 !7
        
        DO J=1,N2DO(MYID)
            DUMMY(J,N)=H2xztL_io(J)
        ENDDO
        N  = N+1 !8
        
        !==========================================
        DO J=1,N2DO(MYID)
            DUMMY(J,N)=DHxztL_io(J)
        ENDDO
        N  = N+1 !9
        
        DO J=1,N2DO(MYID)
            DUMMY(J,N)=PHxztL_io(J)
        ENDDO
        N  = N+0 !9
        
        !==========================================
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 !10,11, 12
                DUMMY(J,L+N)=UHxztL_io(J, L)
            ENDDO
        ENDDO
        N  = N+L1 ! 12
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 13,14,15
                DUMMY(J,L+N)=GHxztL_io(J, L)
            ENDDO
        ENDDO
        N  = N+L1 ! 15
        
        L1 = (NDV*(7-NDV))/2+NDV-3 !6
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 16,17,18,19,20,21
                DUMMY(J,L+N)=U2DHxztL_io(J, L)
            ENDDO
        ENDDO
        N  = N+L1 ! 21
        
        !==========================================
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 22,23,24
                DUMMY(J,L+N)=DhDL1xztL_io(J, L)
            ENDDO
        ENDDO
        N  = N+L1 ! 24
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 25,26,27
                DUMMY(J,L+N)=DhDLPxztL_io(J, L)
            ENDDO
        ENDDO
        N  = N+L1 ! 27
        
        L1 = NDV
        DO J=1,N2DO(MYID)
            DO L=1,L1 ! 28,29,30
                DUMMY(J,L+N)=DTDLKxztL_io(J, L)
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
                    DUMMY(J,M)=DVDL1MxztL_io(J, L, K)
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
                    DUMMY(J,M)=DVDL2MxztL_io(J, L, K)
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
                    DUMMY(J,M)=DVDL1MHxztL_io(J, L, K)
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
                    DUMMY(J,M)=DTDLKUxztL_io(J, L, K)
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
                        DUMMY(J,M)=DVDL1MUxztL_io(J, L, K, S)
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
                        DUMMY(J,M)=DTDLKDVDLxztL_io(J, L, K, S)
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
                        DUMMY(J,M)=DHDLMDVDLxztL_io(J, L, K, S)
                    END DO
                END DO
            ENDDO
        ENDDO
        N = N+L1*L2*L3 !147
        
        IF(N.NE.NSZ) THEN
            WRITE(*,*) 'ERROR in WRT_AVERAGE_VARS_THERMAL_Xperiodic_IO', N, NSZ
        END IF
        
        !==================================
        DFLG = 102
        WRITE(PNTIM,'(1ES15.9)') phyTIME
        WRITE(WRT_AVE_FNM_io,'(A)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_THEL.D'
        if(myid==0) CALL CHKHDL('    IO: WRITING averaged thermal field '//TRIM(WRT_AVE_FNM_io),myid) 
        CALL MPI_TYPE_CREATE_SUBARRAY(NARRAY,SIZES_ARRAY, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, NEWTYPE, IERR)
        CALL MPI_TYPE_COMMIT(NEWTYPE,IERR)
        CALL MPI_FILE_OPEN(ICOMM,WRT_AVE_FNM_io,MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,MYFILE,IERR) 
        
        IF(MYID==0) THEN
            CALL CHKINTHDL('        J IN AVERAGED VARIABLES (J,L)   ',MYID,INTMPI(1))
            CALL CHKINTHDL('        L IN AVERAGED VARIABLES (J,L)   ',MYID,INTMPI(2))
            CALL CHKINTHDL('        ITERG0_IO IN AVERAGED VARIABLES ',MYID,INTMPI(3))
            CALL CHKINTHDL('        NSTATIS_IO IN AVERAGED VARIABLES',MYID,INTMPI(4))
            
            CALL CHKRLHDL ('        phyTIME_TG IN AVERAGED VARIABLES',MYID,RLEMPI(1))
            CALL CHKRLHDL ('        RENL       IN AVERAGED VARIABLES',MYID,RLEMPI(2))
            CALL CHKRLHDL ('        DTL        IN AVERAGED VARIABLES',MYID,RLEMPI(3))
        END IF
        

        OFFSET=0_MPI_OFFSET_KIND
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE(MYFILE,INTMPI,INISIZE,MPI_INTEGER4,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_WRITE(MYFILE,RLEMPI,IRLSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_TYPE_SIZE(MPI_INTEGER4, INTBYTE,IERR)
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, DBLBYTE,IERR)
        
        OFFSET=0_MPI_OFFSET_KIND + INISIZE*INTBYTE +IRLSIZE*DBLBYTE
        
        CALL MPI_FILE_SET_VIEW(MYFILE,OFFSET,MPI_DOUBLE_PRECISION,NEWTYPE,"NATIVE",MPI_INFO_NULL,IERR)
        CALL MPI_FILE_WRITE_ALL(MYFILE,DUMMY,BUFSIZE,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR)
        CALL MPI_FILE_CLOSE(MYFILE,IERR)
        CALL MPI_TYPE_FREE(NEWTYPE,IERR)
        CALL MPI_BARRIER(ICOMM,IERROR)
    

        DEALLOCATE (DUMMY)
        DEALLOCATE (SIZES_ARRAY)
        DEALLOCATE (SUBSIZES)
        DEALLOCATE (STARTS)
        DEALLOCATE (INTMPI)
        DEALLOCATE (RLEMPI)
                
        RETURN
    END SUBROUTINE









!!*******************************************************************************    
!    SUBROUTINE WRT_INSTANT_VARS_io
!        USE MESH_INFO
!        USE FLOW_INFO
!        USE INIT_INFO
!        USE WRT_INFO
!        USE THERMAL_INFO
!        USE POSTPROCESS_INFO
!        IMPLICIT NONE
        
!        CHARACTER(15) :: PNTIM
!        CHARACTER(15) :: STRRK
!        INTEGER(4)    :: DFLG = 10
        
!        IF(WRT_RST_FLG_io) THEN
!            WRITE(PNTIM,'(1ES15.9)') phyTIME
!            WRITE(STRRK,'(1I3.3)') MYID
            
!            IF(MYID.LT.NPSLV) THEN
!                WRITE(WRT_RST_FNM_io,'(A,I4.4,A4)') &
!                    'DNS_inoudomain_INSTANT_'//TRIM(PNTIM)//'_ip'//TRIM(STRRK)//'C.D'
!            ELSE
!                WRITE(WRT_RST_FNM_io,'(A,I4.4,A4)') &
!                    'DNS_inoudomain_INSTANT_'//TRIM(PNTIM)//'_ip'//TRIM(STRRK)//'E.D'
!            END IF
            
!            OPEN(DFLG,FILE=TRIM(WRT_RST_FNM_io),FORM='UNFORMATTED')
!            WRT_RST_FLG_io= .FALSE.
!        ELSE
!            OPEN(DFLG,FILE=TRIM(WRT_RST_FNM_io),FORM='UNFORMATTED',POSITION='APPEND')
!        END IF 
        
!        WRITE(DFLG) MYID, NPTOT
!        WRITE(DFLG) NCL1_io,N2DO(MYID),NCL3,ITERG
!        WRITE(DFLG) phyTIME,REN,DT
!        WRITE(DFLG) G_io
!        WRITE(DFLG) PR_io
!        IF(thermlflg==1) THEN
!            WRITE(DFLG) TEMPERATURE
!            WRITE(DFLG) DENSITY
!            WRITE(DFLG) ENTHALPY
!            WRITE(DFLG) DENSITY0
!        END IF

!        CLOSE(DFLG)
    
!    END SUBROUTINE
    
    
!******************************************************************************* 
    SUBROUTINE WRT_AVERAGE_SPECTRA
        USE MESH_INFO
        USE FLOW_INFO
        USE INIT_INFO
        USE WRT_INFO
        USE THERMAL_INFO
        USE POSTPROCESS_INFO
        IMPLICIT NONE
        
        CHARACTER(15) :: PNTIM
        CHARACTER(15) :: STRRK
        INTEGER(4)    :: DFLG = 11
        
        !IF(ppspectra .NE. 1) RETURN
        IF(MYID .NE. 0)  RETURN
        
        WRITE(PNTIM,'(1ES15.9)') phyTIME
        WRITE(WRT_AVE_FNM_io,'(A,I4.4,A4)') TRIM(filepath2)//'DNS_perixz_AVERAGD_T'//TRIM(PNTIM)//'_SPEC.D'
                    
            
        OPEN(DFLG,FILE=TRIM(WRT_AVE_FNM_io),FORM='UNFORMATTED')
        
        WRITE(DFLG) MYID, NPTOT
        WRITE(DFLG) NCL1_io,NCL2,NCL3,ITERG
        WRITE(DFLG) NSTATIS_io
        WRITE(DFLG) phyTIME,REN,DT
        
        !===========spectra======================
        !================Velocity =====================
        WRITE(DFLG) R11X1_xztLa
        WRITE(DFLG) R22X1_xztLa
        WRITE(DFLG) R33X1_xztLa
        WRITE(DFLG) R12X1_xztLa
        WRITE(DFLG) R13X1_xztLa
        WRITE(DFLG) R23X1_xztLa
        
        WRITE(DFLG) R11X3_xztLa
        WRITE(DFLG) R22X3_xztLa
        WRITE(DFLG) R33X3_xztLa
        WRITE(DFLG) R12X3_xztLa
        WRITE(DFLG) R13X3_xztLa
        WRITE(DFLG) R23X3_xztLa
        
        WRITE(DFLG) ENE11T_xztLa
        WRITE(DFLG) ENE22T_xztLa
        WRITE(DFLG) ENE33T_xztLa
        WRITE(DFLG) ENE12T_xztLa
        WRITE(DFLG) ENE13T_xztLa
        WRITE(DFLG) ENE23T_xztLa
        
        WRITE(DFLG) ENE11Z_xztLa
        WRITE(DFLG) ENE22Z_xztLa
        WRITE(DFLG) ENE33Z_xztLa
        WRITE(DFLG) ENE12Z_xztLa
        WRITE(DFLG) ENE13Z_xztLa
        WRITE(DFLG) ENE23Z_xztLa
        
        !================Voriticity ====================
        WRITE(DFLG) V11X1_xztLa
        WRITE(DFLG) V22X1_xztLa
        WRITE(DFLG) V33X1_xztLa
        WRITE(DFLG) V12X1_xztLa
        WRITE(DFLG) V13X1_xztLa
        WRITE(DFLG) V23X1_xztLa
        
        WRITE(DFLG) V11X3_xztLa
        WRITE(DFLG) V22X3_xztLa
        WRITE(DFLG) V33X3_xztLa
        WRITE(DFLG) V12X3_xztLa
        WRITE(DFLG) V13X3_xztLa
        WRITE(DFLG) V23X3_xztLa
        
        WRITE(DFLG) ENV11T_xztLa
        WRITE(DFLG) ENV22T_xztLa
        WRITE(DFLG) ENV33T_xztLa
        WRITE(DFLG) ENV12T_xztLa
        WRITE(DFLG) ENV13T_xztLa
        WRITE(DFLG) ENV23T_xztLa
        
        WRITE(DFLG) ENV11Z_xztLa
        WRITE(DFLG) ENV22Z_xztLa
        WRITE(DFLG) ENV33Z_xztLa
        WRITE(DFLG) ENV12Z_xztLa
        WRITE(DFLG) ENV13Z_xztLa
        WRITE(DFLG) ENV23Z_xztLa
        
        
        !===============Voriticity & Velocity=========================
        WRITE(DFLG) VO11X1_xztLa
        WRITE(DFLG) VO12X1_xztLa
        WRITE(DFLG) VO13X1_xztLa
        
        WRITE(DFLG) VO21X1_xztLa
        WRITE(DFLG) VO22X1_xztLa
        WRITE(DFLG) VO23X1_xztLa
        
        WRITE(DFLG) VO31X1_xztLa
        WRITE(DFLG) VO32X1_xztLa
        WRITE(DFLG) VO33X1_xztLa
        
        WRITE(DFLG) VO11X3_xztLa
        WRITE(DFLG) VO12X3_xztLa
        WRITE(DFLG) VO13X3_xztLa
        
        WRITE(DFLG) VO21X3_xztLa
        WRITE(DFLG) VO22X3_xztLa
        WRITE(DFLG) VO23X3_xztLa
        
        WRITE(DFLG) VO31X3_xztLa
        WRITE(DFLG) VO32X3_xztLa
        WRITE(DFLG) VO33X3_xztLa
        
        WRITE(DFLG) EVO11T_xztLa
        WRITE(DFLG) EVO12T_xztLa
        WRITE(DFLG) EVO13T_xztLa
        
        WRITE(DFLG) EVO21T_xztLa
        WRITE(DFLG) EVO22T_xztLa
        WRITE(DFLG) EVO23T_xztLa
        
        WRITE(DFLG) EVO31T_xztLa
        WRITE(DFLG) EVO32T_xztLa
        WRITE(DFLG) EVO33T_xztLa
        
        WRITE(DFLG) EVO11Z_xztLa
        WRITE(DFLG) EVO12Z_xztLa
        WRITE(DFLG) EVO13Z_xztLa
        
        WRITE(DFLG) EVO21Z_xztLa
        WRITE(DFLG) EVO22Z_xztLa
        WRITE(DFLG) EVO23Z_xztLa
        
        WRITE(DFLG) EVO31Z_xztLa
        WRITE(DFLG) EVO32Z_xztLa
        WRITE(DFLG) EVO33Z_xztLa
        
        CLOSE(DFLG)
    
    END SUBROUTINE
