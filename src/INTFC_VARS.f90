    SUBROUTINE INTFC_VARS3(L1,L2,NCL1ST,NCL1ED,VARS3)
        use mesh_info
        IMPLICIT NONE
        
        INTEGER(4), INTENT(IN)  :: L1, L2,NCL1ST,NCL1ED
        REAL(WP), INTENT(INOUT) :: VARS3(NCL1ST:NCL1ED,0:N2DO(MYID)+1,NCL3,NDV)
       
        INTEGER(4)  :: I, IK, JJ, JJP
        INTEGER(4)  :: K
        INTEGER(4)  :: N
        INTEGER(4)  :: ITAG
        INTEGER(4)  :: IDESF
        INTEGER(4)  :: IDESB
        INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
        INTEGER(4)  :: NSZ, N11
       
        REAL(WP)    :: BSEN_F((L2-L1+1) * NCL3,NDV)
        REAL(WP)    :: BSEN_L((L2-L1+1) * NCL3,NDV)
        REAL(WP)    :: BREC_F((L2-L1+1) * NCL3,NDV)
        REAL(WP)    :: BREC_L((L2-L1+1) * NCL3,NDV)

        BSEN_F = 0.0_WP
        BSEN_L = 0.0_WP
        BREC_F = 0.0_WP
        BREC_L = 0.0_WP
        NSZ =(L2-L1+1)*NCL3*NDV
        
        N11 = 0
        DO I=L1,L2
            N11 = N11+1
            DO K=1,NCL3
                IK=(N11-1)*NCL3+K
                DO N=1,NDV
                    BSEN_F(IK,N)=VARS3(I,1,         K,N)        !y=local Floor
                    BSEN_L(IK,N)=VARS3(I,N2DO(MYID),K,N)        !y=local ROOF
                ENDDO
            ENDDO
        ENDDO

       
        ITAG=0      
!>     @note rank=2,3,4...SIZE-1 (no b.c.)***************************************      
        IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
            IDESF = MYID+1    !next myid
            IDESB = MYID-1    !last myid
        ELSE IF(MYID.EQ.NPSLV) THEN
            IDESF = 0  
            IDESB = MYID-1   
        ELSE IF (MYID.EQ.0) THEN
            IDESF = MYID+1
            IDESB = NPSLV
        ELSE
        END IF

!>           @note SEND Floor B.C. plane p,u,v,w from current MYID to IDESB=MYID-1
!>              RECEIVE Floor B.C. plane p,u,v,w from IDESF=MYID+1 to current MYID
        CALL  MPI_SENDRECV( BSEN_F(1,1),NSZ,                &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F(1,1),NSZ,&
                MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
 
!>           @note    SEND TOP B.C. plane p,u,v,w from current MYID to IDESF=MYID+1
!>              RECEIVE TOP B.C. plane p,u,v,w from IDESB=MYID-1 to current MYID     
        CALL  MPI_SENDRECV( BSEN_L(1,1),NSZ,                   &
                MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


!>        @note Constructing p,u,v,w of ghost cells y_id=0 and N2NO+1, received from adjacent partitions.      
        N11 = 0
        DO I=L1,L2
            N11 = N11+1
            DO K=1,NCL3
                IK=(N11-1)*NCL3+K
                DO N=1,NDV
                    VARS3(I,0,           K,N)=BREC_L(IK,N)
                    VARS3(I,N2DO(MYID)+1,K,N)=BREC_F(IK,N)
                ENDDO
            ENDDO
        ENDDO
        
        ! Check below
        !IF(myid==NPSLV .or. myid==0) THEN
        !    IF(myid==0)     WRITE(*,*) 'Intf:Y1YNND',myid, VARS3(1,1,1,2)
        !    IF(MYID==NPSLV) WRITE(*,*) 'Intf:Y1YNND',myid, VARS3(1,N2DO(myid)+1,1,2)
        !END IF
        
        !
 
        RETURN
    END SUBROUTINE
    
!==============================================================================================    
    SUBROUTINE INTFC_VARS1(L1,L2,NCL1ST,NCL1ED,VARS1)
        use mesh_info
        IMPLICIT NONE
        
        INTEGER(4), INTENT(IN)  :: L1, L2,NCL1ST,NCL1ED
        REAL(WP), INTENT(INOUT) :: VARS1(NCL1ST:NCL1ED,0:N2DO(MYID)+1,NCL3)

        INTEGER(4)  :: I, IK
        INTEGER(4)  :: K
        INTEGER(4)  :: ITAG
        INTEGER(4)  :: IDESF
        INTEGER(4)  :: IDESB
        INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
        INTEGER(4)  :: NSZ, N11
        
        REAL(WP)    :: BSEN_F((L2-L1+1) * NCL3,NDV)
        REAL(WP)    :: BSEN_L((L2-L1+1) * NCL3,NDV)
        REAL(WP)    :: BREC_F((L2-L1+1) * NCL3,NDV)
        REAL(WP)    :: BREC_L((L2-L1+1) * NCL3,NDV)

        BSEN_F = 0.0_WP
        BSEN_L = 0.0_WP
        BREC_F = 0.0_WP
        BREC_L = 0.0_WP
        NSZ =(L2-L1+1)*NCL3
        
        N11 = 0   
        DO I=L1,L2
            N11 = N11+1
            DO K=1,NCL3
                IK=(N11-1)*NCL3+K
                BSEN_F(IK,1) = VARS1(I,1,         K)  !y=local Floor     
                BSEN_L(IK,1) = VARS1(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
        
        !===================no. wall domains======================================
        ITAG=0          
        IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
            IDESF = MYID+1    !next myid
            IDESB = MYID-1    !last myid
        ELSE IF(MYID.EQ.NPSLV) THEN
            IDESF = 0  
            IDESB = MYID-1   
        
        ELSE IF (MYID.EQ.0) THEN
            IDESF = MYID+1
            IDESB = NPSLV
        ELSE
        END IF

        !==========send/recv data=======================================
        CALL  MPI_SENDRECV( BSEN_F(1,1),NSZ,                     &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F(1,1),NSZ,  &
                MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
    
        CALL  MPI_SENDRECV( BSEN_L(1,1),NSZ,                   &
                MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L(1,1),NSZ, &
                MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


        !=======Constructing data from received===============================     
        N11 = 0   
        DO I=L1,L2
            N11 = N11+1
            DO K=1,NCL3
                IK=(N11-1)*NCL3+K
                VARS1(I,0,           K) = BREC_L(IK,1)
                VARS1(I,N2DO(MYID)+1,K) = BREC_F(IK,1)
            ENDDO
        ENDDO
        
        RETURN
    END SUBROUTINE
