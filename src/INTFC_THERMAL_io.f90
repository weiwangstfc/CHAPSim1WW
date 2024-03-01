    SUBROUTINE INTFC_MFD_THERMAL_io
!>      @NOTE: 
!>         TOP AND BOTTOM WALL HERE IS USING ZERO GRADIENT CONDITION.
!>        THIS TREATMENT IS NOT ACCURATE/PHYSICALLY REASONABLE. 
!>           A B.C. SUBROUTINE IS USED TO TREAT THE TOP AND BOTTOM WALL.
        use thermal_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        INTEGER(4)  :: I, IK
        INTEGER(4)  :: K
        INTEGER(4)  :: ITAG
        INTEGER(4)  :: IDESF
        INTEGER(4)  :: IDESB
        INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
        INTEGER(4)  :: NSZ
        REAL(WP)    :: BSEN_F_io(NCL1_IO * NCL3,NTHERMAL)
        REAL(WP)    :: BSEN_L_io(NCL1_IO * NCL3,NTHERMAL)
        REAL(WP)    :: BREC_F_io(NCL1_IO * NCL3,NTHERMAL)
        REAL(WP)    :: BREC_L_io(NCL1_IO * NCL3,NTHERMAL)
        
        BSEN_F_io = 0.0_WP
        BSEN_L_io = 0.0_WP
        BREC_F_io = 0.0_WP
        BREC_L_io = 0.0_WP
        NSZ = NCL1_IO * NCL3 * NTHERMAL

        !==================1 ENTHALPY=================================== 
        DO I=1,NCL1_io
            DO K=1,NCL3
                IK=(I-1)*NCL3+K
                BSEN_F_io(IK,1)=ENTHALPY(I,1,         K) !y=local Floor     
                BSEN_L_io(IK,1)=ENTHALPY(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
       
        !==================2 DENSITY================================== 
        DO I=1,NCL1_io
            DO K=1,NCL3
                IK=(I-1)*NCL3+K
                BSEN_F_io(IK,2)=DENSITY(I,1,         K)  !y=local Floor     
                BSEN_L_io(IK,2)=DENSITY(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
       
        !==================3 TEMPERATURE================================== 
        DO I=1,NCL1_io
            DO K=1,NCL3
                IK=(I-1)*NCL3+K
                BSEN_F_io(IK,3)=TEMPERATURE(I,1,         K)  !y=local Floor     
                BSEN_L_io(IK,3)=TEMPERATURE(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
       
       
        !==================4 VISCOUSITY================================== 
        DO I=1,NCL1_io
            DO K=1,NCL3
                IK=(I-1)*NCL3+K
                BSEN_F_io(IK,4)=VISCOUSITY(I,1,         K)  !y=local Floor     
                BSEN_L_io(IK,4)=VISCOUSITY(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
       
       !==================5 THERMAL CONDUCTIVITY========================= 
        DO I=1,NCL1_io
            DO K=1,NCL3
                IK=(I-1)*NCL3+K
                BSEN_F_io(IK,5)=THERMCONDT(I,1,         K)  !y=local Floor     
                BSEN_L_io(IK,5)=THERMCONDT(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
       
       !==================1 RHOH================================== 
        DO I=1,NCL1_io
            DO K=1,NCL3
                IK=(I-1)*NCL3+K
                BSEN_F_io(IK,6)=RHOH(I,1,         K)  !y=local Floor     
                BSEN_L_io(IK,6)=RHOH(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
       
       
        !=======TRANSFERING DATA IN THE MAIN PARTITIONS====================================
        ITAG=0            
        IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
            IDESF = MYID+1    !next myid
            IDESB = MYID-1    !last myid

            !SEND Floor B.C. plane p,u,v,w from current MYID to IDESB=MYID-1
            !RECEIVE Floor B.C. plane p,u,v,w from IDESF=MYID+1 to current MYID
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
     
            !SEND TOP B.C. plane p,u,v,w from current MYID to IDESF=MYID+1
            !RECEIVE TOP B.C. plane p,u,v,w from IDESB=MYID-1 to current MYID     
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


            !Constructing of ghost cells y_id=0 and N2NO+1, received from adjacent partitions.

            !==================1 ENTHALPY======================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    ENTHALPY(I,0,           K)=BREC_L_io(IK,1)
                    ENTHALPY(I,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
                ENDDO
            ENDDO
          
            !==================2 DENSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    DENSITY(I,0,           K)=BREC_L_io(IK,2)
                    DENSITY(I,N2DO(MYID)+1,K)=BREC_F_io(IK,2)
                ENDDO
            ENDDO
          
            !==================3 TEMPERATURE==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    TEMPERATURE(I,0,           K)=BREC_L_io(IK,3)
                    TEMPERATURE(I,N2DO(MYID)+1,K)=BREC_F_io(IK,3)
                ENDDO
            ENDDO
          
            !==================4 VISCOUSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    VISCOUSITY(I,0,           K)=BREC_L_io(IK,4)
                    VISCOUSITY(I,N2DO(MYID)+1,K)=BREC_F_io(IK,4)
                ENDDO
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    THERMCONDT(I,0,           K)=BREC_L_io(IK,5)
                    THERMCONDT(I,N2DO(MYID)+1,K)=BREC_F_io(IK,5)
                ENDDO
            ENDDO
          
           !==================1 RHOH==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    RHOH(I,0,K)           =BREC_L_io(IK,6)
                    RHOH(I,N2DO(MYID)+1,K)=BREC_F_io(IK,6)
                ENDDO
            ENDDO
          
 
        ENDIF

        !=======TRANSFERING DATA IN THE TOP WALL PARTITION========================================
        IF(MYID.EQ.NPSLV) THEN
            IDESF = 0  
            IDESB = MYID-1   
          
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)
          

            !@note var(I,N2DO(MYID)+1,K) not used actually.
            !==================1 ENTHALPY===========================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    ENTHALPY(I,0,K)           = BREC_L_io(IK,1)
                    !ENTHALPY(I,N2DO(MYID)+1,K)= 2.0_wp*ENTHALPY(I,N2DO(MYID),  K)- ENTHALPY(I,N2DO(MYID)-1,K)
                ENDDO
            ENDDO
          
            !==================2 DENSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    DENSITY(I,0,K)           =BREC_L_io(IK,2)
                    !DENSITY(I,N2DO(MYID)+1,K)=2.0_wp*DENSITY(I,N2DO(MYID),  K)- DENSITY(I,N2DO(MYID)-1,K)
                ENDDO
            ENDDO
          
            !==================3 TEMPERATURE==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    TEMPERATURE(I,0,K)           =BREC_L_io(IK,3)
                    !TEMPERATURE(I,N2DO(MYID)+1,K)=2.0_wp*TEMPERATURE(I,N2DO(MYID),  K)- TEMPERATURE(I,N2DO(MYID)-1,K)
                ENDDO
            ENDDO
          
            !==================4 VISCOUSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    VISCOUSITY(I,0,K)           =BREC_L_io(IK,4)
                    !VISCOUSITY(I,N2DO(MYID)+1,K)=2.0_wp*VISCOUSITY(I,N2DO(MYID),  K)- VISCOUSITY(I,N2DO(MYID)-1,K)
                ENDDO
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    THERMCONDT(I,0,K)           =BREC_L_io(IK,5)
                    !THERMCONDT(I,N2DO(MYID)+1,K)=2.0_wp*THERMCONDT(I,N2DO(MYID),  K)- THERMCONDT(I,N2DO(MYID)-1,K)
                ENDDO
            ENDDO
          
          
            !==================1 RHOH==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    RHOH(I,0,K)           =BREC_L_io(IK,6)
                    !RHOH(I,N2DO(MYID)+1,K)=2.0_wp*RHOH(I,N2DO(MYID),  K)- RHOH(I,N2DO(MYID)-1,K)
                ENDDO
            ENDDO
            
        ENDIF

        !=======TRANSFERING DATA IN THE BOTTOM WALL PARTITION========================================
        IF (MYID.EQ.0) THEN
            IDESF = MYID+1
            IDESB = NPSLV
         
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

            !==================1 ENTHALPY===========================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    !ENTHALPY(I,0,K)           = 2.0_wp*ENTHALPY(I,1,K) - ENTHALPY(I,2,K)
                    ENTHALPY(I,N2DO(MYID)+1,K)= BREC_F_io(IK,1)
                ENDDO
            ENDDO
          
          
            !==================2 DENSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    !DENSITY(I,0,K)           = 2.0_wp*DENSITY(I,1,K) - DENSITY(I,2,K)
                    DENSITY(I,N2DO(MYID)+1,K)= BREC_F_io(IK,2)
                ENDDO
            ENDDO
          
            !==================3 TEMPERATURE==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    !TEMPERATURE(I,0,K)           = 2.0_wp*TEMPERATURE(I,1,K) - TEMPERATURE(I,2,K)
                    TEMPERATURE(I,N2DO(MYID)+1,K)= BREC_F_io(IK,3)
                ENDDO
            ENDDO
          
            !==================4 VISCOUSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    !VISCOUSITY(I,0,K)           = 2.0_wp*VISCOUSITY(I,1,K) - VISCOUSITY(I,2,K)
                    VISCOUSITY(I,N2DO(MYID)+1,K)= BREC_F_io(IK,4)
                ENDDO
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    !THERMCONDT(I,0,K)           = 2.0_wp*THERMCONDT(I,1,K) - THERMCONDT(I,2,K)
                    THERMCONDT(I,N2DO(MYID)+1,K)= BREC_F_io(IK,5)
                ENDDO
            ENDDO
          
           !==================1 RHOH==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    !RHOH(I,0,K)           = 2.0_wp*RHOH(I,1,K) - RHOH(I,2,K)
                    RHOH(I,N2DO(MYID)+1,K)= BREC_F_io(IK,6)
                ENDDO
            ENDDO
            
        ENDIF

      RETURN
      END SUBROUTINE 
      
     
!******************************************************************************************************************      
    SUBROUTINE INTFC_OUL_THERMAL_io
        use thermal_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        INTEGER(4)  :: I
        INTEGER(4)  :: K
        INTEGER(4)  :: ITAG
        INTEGER(4)  :: IDESF
        INTEGER(4)  :: IDESB
        INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
        INTEGER(4)  :: NSZ
        REAL(WP)    :: BSEN_F_io(1* NCL3,NTHERMAL)
        REAL(WP)    :: BSEN_L_io(1* NCL3,NTHERMAL)
        REAL(WP)    :: BREC_F_io(1* NCL3,NTHERMAL)
        REAL(WP)    :: BREC_L_io(1* NCL3,NTHERMAL)
        
        !REAL(WP)    :: H_TEMP, M_TEMP, K_TEMP, D_TEMP, T_temp
        
        BSEN_F_io = 0.0_WP
        BSEN_L_io = 0.0_WP
        BREC_F_io = 0.0_WP
        BREC_L_io = 0.0_WP
        NSZ = 1* NCL3 * NTHERMAL

        !==================1 ENTHALPY=================================== 
        I=NCL1_io+1
        DO K=1,NCL3
            BSEN_F_io(K,1)=ENTHALPY(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,1)=ENTHALPY(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
        !==================2 DENSITY================================== 
        I=NCL1_io+1
        DO K=1,NCL3
            BSEN_F_io(K,2)=DENSITY(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,2)=DENSITY(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
        !==================3 TEMPERATURE================================== 
        I=NCL1_io+1
        DO K=1,NCL3
            BSEN_F_io(K,3)=TEMPERATURE(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,3)=TEMPERATURE(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       
        !==================4 VISCOUSITY================================== 
        I=NCL1_io+1
        DO K=1,NCL3
            BSEN_F_io(K,4)=VISCOUSITY(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,4)=VISCOUSITY(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       !==================5 THERMAL CONDUCTIVITY========================= 
        I=NCL1_io+1
        DO K=1,NCL3
            BSEN_F_io(K,5)=THERMCONDT(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,5)=THERMCONDT(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       !==================1 RHOH================================== 
        I=NCL1_io+1
        DO K=1,NCL3
            BSEN_F_io(K,6)=RHOH(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,6)=RHOH(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       
        !=======TRANSFERING DATA IN THE MAIN PARTITIONS====================================
        ITAG=0            
        IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
            IDESF = MYID+1    !next myid
            IDESB = MYID-1    !last myid

            !SEND Floor B.C. plane p,u,v,w from current MYID to IDESB=MYID-1
            !RECEIVE Floor B.C. plane p,u,v,w from IDESF=MYID+1 to current MYID
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
     
            !SEND TOP B.C. plane p,u,v,w from current MYID to IDESF=MYID+1
            !RECEIVE TOP B.C. plane p,u,v,w from IDESB=MYID-1 to current MYID     
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


            !Constructing of ghost cells y_id=0 and N2NO+1, received from adjacent partitions.

            !==================1 ENTHALPY======================================
            I=NCL1_io+1
            DO K=1,NCL3
                ENTHALPY(I,0,           K)=BREC_L_io(K,1)
                ENTHALPY(I,N2DO(MYID)+1,K)=BREC_F_io(K,1)
            ENDDO
          
            !==================2 DENSITY==================================
            I=NCL1_io+1
            DO K=1,NCL3
                DENSITY(I,0,           K)=BREC_L_io(K,2)
                DENSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,2)
            ENDDO
          
            !==================3 TEMPERATURE==================================
            I=NCL1_io+1
            DO K=1,NCL3
                TEMPERATURE(I,0,           K)=BREC_L_io(K,3)
                TEMPERATURE(I,N2DO(MYID)+1,K)=BREC_F_io(K,3)
            ENDDO
          
            !==================4 VISCOUSITY==================================
            I=NCL1_io+1
            DO K=1,NCL3
                VISCOUSITY(I,0,           K)=BREC_L_io(K,4)
                VISCOUSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,4)
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            I=NCL1_io+1
            DO K=1,NCL3
                THERMCONDT(I,0,           K)=BREC_L_io(K,5)
                THERMCONDT(I,N2DO(MYID)+1,K)=BREC_F_io(K,5)
            ENDDO
          
           !==================1 RHOH==================================
            I=NCL1_io+1
            DO K=1,NCL3
                RHOH(I,0,K)           =BREC_L_io(K,6)
                RHOH(I,N2DO(MYID)+1,K)=BREC_F_io(K,6)
            ENDDO
          
 
        ENDIF

        !=======TRANSFERING DATA IN THE TOP WALL PARTITION========================================
        IF(MYID.EQ.NPSLV) THEN
            IDESF = 0  
            IDESB = MYID-1   
          
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)
          

            !@note var(I,N2DO(MYID)+1,K) not used actually.
            !==================1 ENTHALPY===========================
            I=NCL1_io+1
            DO K=1,NCL3
                ENTHALPY(I,0,K)           =BREC_L_io(K,1)
                !ENTHALPY(I,N2DO(MYID)+1,K)=2.0_wp*ENTHALPY(I,N2DO(MYID),  K)- ENTHALPY(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================2 DENSITY==================================
            I=NCL1_io+1
            DO K=1,NCL3
                DENSITY(I,0,K)           =BREC_L_io(K,2)
                !DENSITY(I,N2DO(MYID)+1,K)=2.0_wp*DENSITY(I,N2DO(MYID),  K)- DENSITY(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================3 TEMPERATURE==================================
            I=NCL1_io+1
            DO K=1,NCL3
                TEMPERATURE(I,0,K)           =BREC_L_io(K,3)
                !TEMPERATURE(I,N2DO(MYID)+1,K)=2.0_wp*TEMPERATURE(I,N2DO(MYID),  K)- TEMPERATURE(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================4 VISCOUSITY==================================
            I=NCL1_io+1
            DO K=1,NCL3
                VISCOUSITY(I,0,K)           =BREC_L_io(K,4)
                !VISCOUSITY(I,N2DO(MYID)+1,K)=2.0_wp*VISCOUSITY(I,N2DO(MYID),  K)- VISCOUSITY(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            I=NCL1_io+1
            DO K=1,NCL3
                THERMCONDT(I,0,K)           =BREC_L_io(K,5)
                !THERMCONDT(I,N2DO(MYID)+1,K)=2.0_wp*THERMCONDT(I,N2DO(MYID),  K)- THERMCONDT(I,N2DO(MYID)-1,K)
            ENDDO
          
          
            !==================1 RHOH==================================
            I=NCL1_io+1
            DO K=1,NCL3
                RHOH(I,0,K)           =BREC_L_io(K,6)
                !RHOH(I,N2DO(MYID)+1,K)=2.0_wp*RHOH(I,N2DO(MYID),  K)- RHOH(I,N2DO(MYID)-1,K)
            ENDDO
            
        ENDIF

        !=======TRANSFERING DATA IN THE BOTTOM WALL PARTITION========================================
        IF (MYID.EQ.0) THEN
            IDESF = MYID+1
            IDESB = NPSLV
         
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

            !==================1 ENTHALPY===========================
            I=NCL1_io+1
            DO K=1,NCL3
                !ENTHALPY(I,0,K)           =2.0_wp*ENTHALPY(I,1,K) - ENTHALPY(I,2,K)
                ENTHALPY(I,N2DO(MYID)+1,K)=BREC_F_io(K,1)
            ENDDO
          
          
            !==================2 DENSITY==================================
            I=NCL1_io+1
            DO K=1,NCL3
                !DENSITY(I,0,K)           =2.0_wp*DENSITY(I,1,K) - DENSITY(I,2,K)
                DENSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,2)
            ENDDO
          
            !==================3 TEMPERATURE==================================
            I=NCL1_io+1
            DO K=1,NCL3
                !TEMPERATURE(I,0,K)           =2.0_wp*TEMPERATURE(I,1,K) - TEMPERATURE(I,2,K)
                TEMPERATURE(I,N2DO(MYID)+1,K)=BREC_F_io(K,3)
            ENDDO
          
            !==================4 VISCOUSITY==================================
            I=NCL1_io+1
            DO K=1,NCL3
                !VISCOUSITY(I,0,K)           =2.0_wp*VISCOUSITY(I,1,K) - VISCOUSITY(I,2,K)
                VISCOUSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,4)
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            I=NCL1_io+1
            DO K=1,NCL3
                !THERMCONDT(I,0,K)           =2.0_wp*THERMCONDT(I,1,K) - THERMCONDT(I,2,K)
                THERMCONDT(I,N2DO(MYID)+1,K)=BREC_F_io(K,5)
            ENDDO
          
           !==================1 RHOH==================================
            I=NCL1_io+1
            DO K=1,NCL3
                !RHOH(I,0,K)           =2.0_wp*RHOH(I,1,K) - RHOH(I,2,K)
                RHOH(I,N2DO(MYID)+1,K)=BREC_F_io(K,6)
            ENDDO
            
        ENDIF
      
      
      RETURN
      END SUBROUTINE 
      
!******************************************************************************************************************      
    SUBROUTINE INTFC_INL_THERMAL_io
        use thermal_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        INTEGER(4)  :: I
        INTEGER(4)  :: K
        INTEGER(4)  :: ITAG
        INTEGER(4)  :: IDESF
        INTEGER(4)  :: IDESB
        INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
        INTEGER(4)  :: NSZ
        REAL(WP)    :: BSEN_F_io(1* NCL3,NTHERMAL)
        REAL(WP)    :: BSEN_L_io(1* NCL3,NTHERMAL)
        REAL(WP)    :: BREC_F_io(1* NCL3,NTHERMAL)
        REAL(WP)    :: BREC_L_io(1* NCL3,NTHERMAL)
        
        !REAL(WP)    :: H_TEMP, M_TEMP, K_TEMP, D_TEMP, T_temp
        
        BSEN_F_io = 0.0_WP
        BSEN_L_io = 0.0_WP
        BREC_F_io = 0.0_WP
        BREC_L_io = 0.0_WP
        NSZ = 1* NCL3 * NTHERMAL

        !==================1 ENTHALPY=================================== 
        I=0
        DO K=1,NCL3
            BSEN_F_io(K,1)=ENTHALPY(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,1)=ENTHALPY(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
        !==================2 DENSITY================================== 
        I=0
        DO K=1,NCL3
            BSEN_F_io(K,2)=DENSITY(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,2)=DENSITY(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
        !==================3 TEMPERATURE================================== 
        I=0
        DO K=1,NCL3
            BSEN_F_io(K,3)=TEMPERATURE(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,3)=TEMPERATURE(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       
        !==================4 VISCOUSITY================================== 
        I=0
        DO K=1,NCL3
            BSEN_F_io(K,4)=VISCOUSITY(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,4)=VISCOUSITY(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       !==================5 THERMAL CONDUCTIVITY========================= 
        I=0
        DO K=1,NCL3
            BSEN_F_io(K,5)=THERMCONDT(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,5)=THERMCONDT(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       !==================1 RHOH================================== 
        I=0
        DO K=1,NCL3
            BSEN_F_io(K,6)=RHOH(I,1,         K)  !y=local Floor     
            BSEN_L_io(K,6)=RHOH(I,N2DO(MYID),K)  !y=local ROOF
        ENDDO
       
       
        !=======TRANSFERING DATA IN THE MAIN PARTITIONS====================================
        ITAG=0            
        IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
            IDESF = MYID+1    !next myid
            IDESB = MYID-1    !last myid

            !SEND Floor B.C. plane p,u,v,w from current MYID to IDESB=MYID-1
            !RECEIVE Floor B.C. plane p,u,v,w from IDESF=MYID+1 to current MYID
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
     
            !SEND TOP B.C. plane p,u,v,w from current MYID to IDESF=MYID+1
            !RECEIVE TOP B.C. plane p,u,v,w from IDESB=MYID-1 to current MYID     
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


            !Constructing of ghost cells y_id=0 and N2NO+1, received from adjacent partitions.

            !==================1 ENTHALPY======================================
            I=0
            DO K=1,NCL3
                ENTHALPY(I,0,           K)=BREC_L_io(K,1)
                ENTHALPY(I,N2DO(MYID)+1,K)=BREC_F_io(K,1)
            ENDDO
          
            !==================2 DENSITY==================================
            I=0
            DO K=1,NCL3
                DENSITY(I,0,           K)=BREC_L_io(K,2)
                DENSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,2)
            ENDDO
          
            !==================3 TEMPERATURE==================================
            I=0
            DO K=1,NCL3
                TEMPERATURE(I,0,           K)=BREC_L_io(K,3)
                TEMPERATURE(I,N2DO(MYID)+1,K)=BREC_F_io(K,3)
            ENDDO
          
            !==================4 VISCOUSITY==================================
            I=0
            DO K=1,NCL3
                VISCOUSITY(I,0,           K)=BREC_L_io(K,4)
                VISCOUSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,4)
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            I=0
            DO K=1,NCL3
                THERMCONDT(I,0,           K)=BREC_L_io(K,5)
                THERMCONDT(I,N2DO(MYID)+1,K)=BREC_F_io(K,5)
            ENDDO
          
           !==================1 RHOH==================================
            I=0
            DO K=1,NCL3
                RHOH(I,0,K)           =BREC_L_io(K,6)
                RHOH(I,N2DO(MYID)+1,K)=BREC_F_io(K,6)
            ENDDO
          
 
        ENDIF

        !=======TRANSFERING DATA IN THE TOP WALL PARTITION========================================
        IF(MYID.EQ.NPSLV) THEN
            IDESF = 0  
            IDESB = MYID-1   
          
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)
          

            !@note var(I,N2DO(MYID)+1,K) not used actually.
            !==================1 ENTHALPY===========================
            I=0
            DO K=1,NCL3
                ENTHALPY(I,0,K)           =BREC_L_io(K,1)
                !ENTHALPY(I,N2DO(MYID)+1,K)=2.0_wp*ENTHALPY(I,N2DO(MYID),  K)- ENTHALPY(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================2 DENSITY==================================
            I=0
            DO K=1,NCL3
                DENSITY(I,0,K)           =BREC_L_io(K,2)
                !DENSITY(I,N2DO(MYID)+1,K)=2.0_wp*DENSITY(I,N2DO(MYID),  K)- DENSITY(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================3 TEMPERATURE==================================
            I=0
            DO K=1,NCL3
                TEMPERATURE(I,0,K)           =BREC_L_io(K,3)
                !TEMPERATURE(I,N2DO(MYID)+1,K)=2.0_wp*TEMPERATURE(I,N2DO(MYID),  K)- TEMPERATURE(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================4 VISCOUSITY==================================
            I=0
            DO K=1,NCL3
                VISCOUSITY(I,0,K)           =BREC_L_io(K,4)
                !VISCOUSITY(I,N2DO(MYID)+1,K)=2.0_wp*VISCOUSITY(I,N2DO(MYID),  K)- VISCOUSITY(I,N2DO(MYID)-1,K)
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            I=0
            DO K=1,NCL3
                THERMCONDT(I,0,K)           =BREC_L_io(K,5)
                !THERMCONDT(I,N2DO(MYID)+1,K)=2.0_wp*THERMCONDT(I,N2DO(MYID),  K)- THERMCONDT(I,N2DO(MYID)-1,K)
            ENDDO
          
          
            !==================1 RHOH==================================
            I=0
            DO K=1,NCL3
                RHOH(I,0,K)           =BREC_L_io(K,6)
                !RHOH(I,N2DO(MYID)+1,K)=2.0_wp*RHOH(I,N2DO(MYID),  K)- RHOH(I,N2DO(MYID)-1,K)
            ENDDO
            
        ENDIF

        !=======TRANSFERING DATA IN THE BOTTOM WALL PARTITION========================================
        IF (MYID.EQ.0) THEN
            IDESF = MYID+1
            IDESB = NPSLV
         
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

            !==================1 ENTHALPY===========================
            I=0
            DO K=1,NCL3
                !ENTHALPY(I,0,K)           =2.0_wp*ENTHALPY(I,1,K) - ENTHALPY(I,2,K)
                ENTHALPY(I,N2DO(MYID)+1,K)=BREC_F_io(K,1)
            ENDDO
          
          
            !==================2 DENSITY==================================
            I=0
            DO K=1,NCL3
                !DENSITY(I,0,K)           =2.0_wp*DENSITY(I,1,K) - DENSITY(I,2,K)
                DENSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,2)
            ENDDO
          
            !==================3 TEMPERATURE==================================
            I=0
            DO K=1,NCL3
                !TEMPERATURE(I,0,K)           =2.0_wp*TEMPERATURE(I,1,K) - TEMPERATURE(I,2,K)
                TEMPERATURE(I,N2DO(MYID)+1,K)=BREC_F_io(K,3)
            ENDDO
          
            !==================4 VISCOUSITY==================================
            I=0
            DO K=1,NCL3
                !VISCOUSITY(I,0,K)           =2.0_wp*VISCOUSITY(I,1,K) - VISCOUSITY(I,2,K)
                VISCOUSITY(I,N2DO(MYID)+1,K)=BREC_F_io(K,4)
            ENDDO
          
            !==================5 THERMAL CONDUCTIVITY===========================
            I=0
            DO K=1,NCL3
                !THERMCONDT(I,0,K)           =2.0_wp*THERMCONDT(I,1,K) - THERMCONDT(I,2,K)
                THERMCONDT(I,N2DO(MYID)+1,K)=BREC_F_io(K,5)
            ENDDO
          
           !==================1 RHOH==================================
            I=0
            DO K=1,NCL3
                !RHOH(I,0,K)           =2.0_wp*RHOH(I,1,K) - RHOH(I,2,K)
                RHOH(I,N2DO(MYID)+1,K)=BREC_F_io(K,6)
            ENDDO
            
        ENDIF
      
      
      RETURN
      END SUBROUTINE 



!========================================
    SUBROUTINE INTFC_MFD_DENSITY_io
!>      @NOTE: 
!>         TOP AND BOTTOM WALL HERE IS USING ZERO GRADIENT CONDITION.
!>        THIS TREATMENT IS NOT ACCURATE/PHYSICALLY REASONABLE. 
!>           A B.C. SUBROUTINE IS USED TO TREAT THE TOP AND BOTTOM WALL.
        use thermal_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        INTEGER(4)  :: I, IK
        INTEGER(4)  :: K
        INTEGER(4)  :: ITAG
        INTEGER(4)  :: IDESF
        INTEGER(4)  :: IDESB
        INTEGER(4)  :: TRC_STS(MPI_STATUS_SIZE)
        INTEGER(4)  :: NSZ
        REAL(WP)    :: BSEN_F_io(NCL1_IO * NCL3,1)
        REAL(WP)    :: BSEN_L_io(NCL1_IO * NCL3,1)
        REAL(WP)    :: BREC_F_io(NCL1_IO * NCL3,1)
        REAL(WP)    :: BREC_L_io(NCL1_IO * NCL3,1)
        
        BSEN_F_io = 0.0_WP
        BSEN_L_io = 0.0_WP
        BREC_F_io = 0.0_WP
        BREC_L_io = 0.0_WP
        NSZ = NCL1_IO * NCL3 * 1

        
        !==================2 DENSITY================================== 
        DO I=1,NCL1_io
            DO K=1,NCL3
                IK=(I-1)*NCL3+K
                BSEN_F_io(IK,1)=DENSITY(I,1,         K)  !y=local Floor     
                BSEN_L_io(IK,1)=DENSITY(I,N2DO(MYID),K)  !y=local ROOF
            ENDDO
        ENDDO
       
        !=======TRANSFERING DATA IN THE MAIN PARTITIONS====================================
        ITAG=0            
        IF ( (MYID.GT.0) .AND. (MYID.LT.NPSLV) ) THEN
            IDESF = MYID+1    !next myid
            IDESB = MYID-1    !last myid

            !SEND Floor B.C. plane p,u,v,w from current MYID to IDESB=MYID-1
            !RECEIVE Floor B.C. plane p,u,v,w from IDESF=MYID+1 to current MYID
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,ICOMM,TRC_STS,IERROR)
     
            !SEND TOP B.C. plane p,u,v,w from current MYID to IDESF=MYID+1
            !RECEIVE TOP B.C. plane p,u,v,w from IDESB=MYID-1 to current MYID     
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                   &
                    MPI_DOUBLE_PRECISION,IDESF,ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB,ITAG,ICOMM,TRC_STS,IERROR)


            !Constructing of ghost cells y_id=0 and N2NO+1, received from adjacent partitions.
            !==================2 DENSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    DENSITY(I,0,           K)=BREC_L_io(IK,1)
                    DENSITY(I,N2DO(MYID)+1,K)=BREC_F_io(IK,1)
                ENDDO
            ENDDO
          
        ENDIF

        !=======TRANSFERING DATA IN THE TOP WALL PARTITION========================================
        IF(MYID.EQ.NPSLV) THEN
            IDESF = 0  
            IDESB = MYID-1   
          
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,&
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,               &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ, &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)
          

            !@note var(I,N2DO(MYID)+1,K) not used actually.
            !==================2 DENSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    DENSITY(I,0,K)           =BREC_L_io(IK,1)
                    DENSITY(I,N2DO(MYID)+1,K)=2.0_wp*DENSITY(I,N2DO(MYID),  K)- DENSITY(I,N2DO(MYID)-1,K)
                ENDDO
            ENDDO
          
        ENDIF

        !=======TRANSFERING DATA IN THE BOTTOM WALL PARTITION========================================
        IF (MYID.EQ.0) THEN
            IDESF = MYID+1
            IDESB = NPSLV
         
            CALL  MPI_SENDRECV( BSEN_F_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESB , ITAG,BREC_F_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESF, ITAG,ICOMM,TRC_STS,IERROR)
            CALL  MPI_SENDRECV( BSEN_L_io(1,1),NSZ,                 &
                    MPI_DOUBLE_PRECISION, IDESF , ITAG,BREC_L_io(1,1),NSZ,  &
                    MPI_DOUBLE_PRECISION,IDESB, ITAG,ICOMM,TRC_STS,IERROR)

        
            !==================2 DENSITY==================================
            DO I=1,NCL1_io
                DO K=1,NCL3
                    IK=(I-1)*NCL3+K
                    DENSITY(I,0,K)           = 2.0_wp*DENSITY(I,1,K) - DENSITY(I,2,K)
                    DENSITY(I,N2DO(MYID)+1,K)= BREC_F_io(IK,1)
                ENDDO
            ENDDO
          
        ENDIF
      
        RETURN
    END SUBROUTINE 
        
      
      
      
