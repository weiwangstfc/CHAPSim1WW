!***********************************************************************
!> @author 
!> Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Allocate and deallocate common variables.
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! 18/12/2013- Initial Version
!**********************************************************************
    SUBROUTINE MEM_ALLOCAT
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
         
         
        IF(MYID==0) CALL CHKHDL   ('       Allocating index...',MYID) 
        ALLOCATE ( JCL2G( N2DO(0) ) )        ; JCL2G = 0
        ALLOCATE ( KCL2G( N3DO(0) ) )        ; KCL2G = 0
        ALLOCATE ( JG2IP( NCL2)     )        ; JG2IP = 0
        ALLOCATE ( JG2LC( NCL2)     )        ; JG2LC = 0
        
        ALLOCATE ( JGMOV( MGRID)     )        ; JGMOV = 0
        
        MEMPC_byte = MEMPC_byte + ( N2DO(0) + N3DO(0) + NCL2 ) * 4
      
        ALLOCATE ( KPV(NCL3) )           ; KPV = 0
        ALLOCATE ( KMV(NCL3) )           ; KMV = 0
        ALLOCATE ( KSYM(NCL3) )          ; KSYM= 0
        IF(icase.eq.IBOX3P) THEN
            ALLOCATE ( JGPV(NCL2) )        ; JGPV = 0
            ALLOCATE ( JGMV(NCL2) )        ; JGMV = 0
        ELSE
            ALLOCATE ( JGPV(0:NCL2) )        ; JGPV = 0
            ALLOCATE ( JGMV(1:NND2) )        ; JGMV = 0
        END IF
        
        ALLOCATE ( JLPV(0:N2DO(0)) )     ; JLPV = 0
        ALLOCATE ( JLMV(1:N2DO(0)+1) )     ; JLMV = 0
        MEMPC_byte = MEMPC_byte + ( 3*NCL3 + 2*NND2 + 2*N2DO(0) ) * 4
      
        IF(MYID==0) CALL CHKHDL   ('       Allocating TDMA coefficients...',MYID) 
        ALLOCATE ( CFLVIS(N2DO(0)) )    ; CFLVIS  = 0.0_WP
        ALLOCATE ( AMPH(NCL2) )         ; AMPH  = 0.0_WP
        ALLOCATE ( ACPH(NCL2) )         ; ACPH  = 0.0_WP
        ALLOCATE ( APPH(NCL2) )         ; APPH  = 0.0_WP
        MEMPC_byte = MEMPC_byte + 4*NCL2 * 8  
        
        !ALLOCATE ( DPDYWAL(NCL1_io,NCL3,2) ) ; DPDYWAL = 0.0_WP
        !IF(TGFLOWflg .or. (visthemflg==visimplicit)) THEN
        ALLOCATE ( AMVR(NND2,3) )       ; AMVR  = 0.0_WP
        ALLOCATE ( ACVR(NND2,3) )       ; ACVR  = 0.0_WP
        ALLOCATE ( APVR(NND2,3) )       ; APVR  = 0.0_WP
        MEMPC_byte = MEMPC_byte + ( 3*3*NND2) * 8  
        !END IF
        
        IF(MYID==0) CALL CHKHDL   ('       Allocating Coordinates...',MYID) 
        ALLOCATE ( ZCC(NCL3) )          ; ZCC   = 0.0_WP
        ALLOCATE ( ZND(NND3) )          ; ZND   = 0.0_WP
        ALLOCATE ( YND(NND2) )          ; YND   = 0.0_WP
        ALLOCATE ( YCC(NCL2) )          ; YCC   = 0.0_WP
        ALLOCATE ( RCCI1(0:NND2) )      ; RCCI1    = 1.0_WP
        ALLOCATE ( RCCI2(0:NND2) )      ; RCCI2    = 1.0_WP
        ALLOCATE ( RNDI1(1:NND2) )      ; RNDI1    = 1.0_WP
        ALLOCATE ( RNDI2(1:NND2) )      ; RNDI2    = 1.0_WP
        
        ALLOCATE ( DYFI(NCL2) )         ; DYFI  = 0.0_WP
        ALLOCATE ( DYCI(NND2) )         ; DYCI  = 0.0_WP
        MEMPC_byte = MEMPC_byte + (NND3+NND2*2+NCL2*2+4*(NND2+1)) * 8  
        
        ALLOCATE ( YCL2ND_WFF(NND2) )   ; YCL2ND_WFF = 0.0_WP
        ALLOCATE ( YCL2ND_WFB(NND2) )   ; YCL2ND_WFB = 0.0_WP
        MEMPC_byte = MEMPC_byte + NND2*2*8  
        
        
        ALLOCATE ( Vini(NCL2)              ) ; Vini = 0.0_WP
        ALLOCATE ( UU (N2DO(0),NDV,2)      ) ; UU   = 0.0_WP
        MEMPC_byte = MEMPC_byte + ( NCL2 + N2DO(0)*NDV*2) * 8
        

        IF(TGFLOWflg) THEN
            IF(MYID==0) CALL CHKHDL   ('       Allocating TG instant variables...',MYID) 
            ALLOCATE ( IPV_tg( NCL1_tg ) )          ;  IPV_tg = 0
            ALLOCATE ( IMV_tg( NCL1_tg ) )          ;  IMV_tg = 0
            MEMPC_byte = MEMPC_byte + ( 2*NCL1_tg) * 4  
            
            ALLOCATE ( XND_tg( NND1_tg ) )          ;  XND_tg= 0.0_WP
            ALLOCATE ( XCC_tg( NCL1_tg ) )          ;  XCC_tg= 0.0_WP
            MEMPC_byte = MEMPC_byte + ( NCL1_tg + NND1_tg) * 8
            
            ALLOCATE ( Q_tg   (NCL1_tg,0:N2DO(0)+1,NCL3,NDV) )   ; Q_tg       = 0.0_WP
            ALLOCATE ( PR_tg  (NCL1_tg,0:N2DO(0)+1,NCL3    ) )   ; PR_tg      = 0.0_WP
            ALLOCATE ( Qtmp_tg(NCL1_tg,0:N2DO(0)+1,NCL3    ) )   ; Qtmp_tg    = 0.0_WP
            ALLOCATE ( DPH_tg (NCL1_tg,0:N2DO(0)+1,NCL3    ) )   ; DPH_tg     = 0.0_WP
        
            ALLOCATE ( CONVH0_tg   (NCL1_tg,N2DO(0),NCL3,NDV)  ) ; CONVH0_tg  = 0.0_WP
            ALLOCATE ( RHS_tg      (NCL1_tg,N2DO(0),NCL3)      ) ; RHS_tg= 0.0_WP
            ALLOCATE ( RHSLLPHI_tg (NCL1_tg,N2DO(0),NCL3)      ) ; RHSLLPHI_tg= 0.0_WP
      
            MEMPC_byte = MEMPC_byte + ( NCL1_tg *  (N2DO(0) + 2) * NCL3 * (2*NDV+5) ) * 8
            
             !====for postprocess=================================
             IF(MYID==0) CALL CHKHDL   ('       Allocating TG x_z averaged variables...',MYID) 
            ALLOCATE ( U1xzL_tg( N2DO(0), NDV+1 ) )                              ; U1xzL_tg = 0.0_wp
            ALLOCATE ( UPxzL_tg( N2DO(0), NDV   ) )                              ; UPxzL_tg = 0.0_wp
            ALLOCATE ( U2xzL_tg( N2DO(0), NDV*(7-NDV)/2+NDV-3                ) ) ; U2xzL_tg=0.0_WP
            ALLOCATE ( U3xzL_tg( N2DO(0), NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) ) ; U3xzL_tg=0.0_WP
    
            ALLOCATE ( DVDL1xzL_tg( N2DO(0), NDV, NDV  ) )                        ; DVDL1xzL_tg = 0.0_WP
            ALLOCATE ( DVDLPxzL_tg( N2DO(0), NDV, NDV  ) )                        ; DVDLPxzL_tg = 0.0_WP
            ALLOCATE ( DVDL2xzL_tg( N2DO(0), NDV*NDV, NDV*NDV  ) )        ; DVDL2xzL_tg = 0.0_WP
            
            IF(MYID==0) CALL CHKHDL   ('       Allocating TG x_z_t averaged variables...',MYID) 
            ALLOCATE ( U1xztL_tg( N2DO(0), NDV+1 ) )                              ; U1xztL_tg = 0.0_WP
            ALLOCATE ( UPxztL_tg( N2DO(0), NDV   ) )                              ; UPxztL_tg = 0.0_WP
            ALLOCATE ( U2xztL_tg( N2DO(0), NDV*(7-NDV)/2+NDV-3                ) ) ; U2xztL_tg=0.0_WP
            ALLOCATE ( U3xztL_tg( N2DO(0), NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) ) ; U3xztL_tg=0.0_WP
    
            ALLOCATE ( DVDL1xztL_tg( N2DO(0), NDV, NDV  ) )                       ; DVDL1xztL_tg = 0.0_WP
            ALLOCATE ( DVDLPxztL_tg( N2DO(0), NDV, NDV  ) )                       ; DVDLPxztL_tg = 0.0_WP
            ALLOCATE ( DVDL2xztL_tg( N2DO(0), NDV*NDV, NDV*NDV  ) )               ; DVDL2xztL_tg = 0.0_WP
            
            ALLOCATE ( QUADUVxzL_tg (N2DO(0),4,QUADHN) )                          ; QUADUVxzL_tg =0.0_wp
            ALLOCATE ( QUADVzxzL_tg (N2DO(0),4,QUADHN) )                          ; QUADVzxzL_tg =0.0_wp
            ALLOCATE ( QUADTKxzL_tg (N2DO(0),4,QUADHN) )                          ; QUADTKxzL_tg =0.0_wp
            ALLOCATE ( QUADDRxzL_tg (N2DO(0),4,QUADHN) )                          ; QUADDRxzL_tg =0.0_wp
            
            ALLOCATE ( QUADUVxztL_tg (N2DO(0),4,QUADHN) )                          ; QUADUVxztL_tg =0.0_wp
            ALLOCATE ( QUADVzxztL_tg (N2DO(0),4,QUADHN) )                          ; QUADVzxztL_tg =0.0_wp
            ALLOCATE ( QUADTKxztL_tg (N2DO(0),4,QUADHN) )                          ; QUADTKxztL_tg =0.0_wp
            ALLOCATE ( QUADDRxztL_tg (N2DO(0),4,QUADHN) )                          ; QUADDRxztL_tg =0.0_wp
                
    
            MEMPC_byte = MEMPC_byte + 2*N2DO(0)* ( &
                         NDV+1 + NDV + NDV + NDV*(7-NDV)/2+NDV-3 + NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8 + &
                         NDV*NDV*2 + (NDV*(7-NDV)/2+NDV-3)*NDV +4*4*QUADHN) * 8
        END IF
        
      
        IF(IOFLOWflg) THEN
            IF(MYID==0) CALL CHKHDL   ('       Allocating IO instant flow variables...',MYID) 
            ALLOCATE ( XND_io(1:NND1_io ) )                         ; XND_io= 0.0_WP
            ALLOCATE ( XCC_io(1:NCL1_io ) )                         ; XCC_io= 0.0_WP
            MEMPC_byte = MEMPC_byte + (NND1_io+NCL1_io)*8

            ALLOCATE ( Q_io   (NCL1S:NCL1E,0:N2DO(0)+1,NCL3,NDV) ) ; Q_io  = 0.0_WP
            ALLOCATE ( G_io   (NCL1S:NCL1E,0:N2DO(0)+1,NCL3,NDV) ) ; G_io  = 0.0_WP
            ALLOCATE ( G0_io   (NCL1S:NCL1E,0:N2DO(0)+1,NCL3,NDV) ) ; G0_io  = 0.0_WP
            ALLOCATE ( PR_io  (NCL1S:NCL1E,0:N2DO(0)+1,NCL3    ) ) ; PR_io = 0.0_WP
            ALLOCATE ( DPH_io (NCL1S:NCL1E,0:N2DO(0)+1,NCL3    ) ) ; DPH_io= 0.0_WP
            ALLOCATE ( Qtmp_io(NCL1S:NCL1E,0:N2DO(0)+1,NCL3    ) ) ; Qtmp_io = 0.0_WP
            
            IF(weightedpressure==1) THEN
                ALLOCATE ( PR0_io  (NCL1S:NCL1E,0:N2DO(0)+1,NCL3,2    ) ) ; PR0_io = 0.0_WP
            END IF
         
            ALLOCATE ( RHS_io      (NCL1_io,N2DO(0),NCL3        )  )   ; RHS_io      = 0.0_WP
            ALLOCATE ( RHSLLPHI_io (NCL1_io,N2DO(0),NCL3        )  )   ; RHSLLPHI_io = 0.0_WP
            ALLOCATE ( DivU_io     (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)  )   ; DivU_io     = 0.0_WP
            ALLOCATE ( EXPLT0_io   (NCL1_io,N2DO(0),NCL3,NDV    )  )   ; EXPLT0_io   = 0.0_WP
            MEMPC_byte = MEMPC_byte + ( (NCL1_io+2)*(N2DO(0)+2)*NCL3*(3*NDV+6) ) * 8
            
            IF(TGFLOWflg) THEN
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO Inlet/Outlet variables...',MYID) 
                ALLOCATE ( BC_CONV0    ( 0:N2DO(0)+1,NCL3,NDV    ) )   ; BC_CONV0     = 0.0_WP
                ALLOCATE ( BC_CONV0_ENG( 0:N2DO(0)+1,NCL3        ) )   ; BC_CONV0_ENG = 0.0_WP
                ALLOCATE ( BC_TDMA     ( 3,0:N2DO(0)+1,NCL3,NDV  ) )   ; BC_TDMA      = 0.0_WP
                ALLOCATE ( BC_U_SSTAR  (     N2DO(0),  NCL3,NDV  ) )   ; BC_U_SSTAR   = 0.0_WP
                MEMPC_byte = MEMPC_byte + (N2DO(0)+2)*NCL3*(NDV*3+1)*8
            END IF
          
            !================thermal info============================= 
            IF(thermlflg==1) THEN
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO instantanous thermal variables...',MYID) 
                ALLOCATE (RHOH       (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)  )   ; RHOH        = 0.0_WP
                ALLOCATE (RHOH0      (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)  )   ; RHOH0       = 0.0_WP
                ALLOCATE (ENTHALPY   (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)  )   ; ENTHALPY    = 0.0_WP
                ALLOCATE (TEMPERATURE(NCL1S:NCL1E,0:N2DO(0)+1,NCL3)  )   ; TEMPERATURE = 1.0_WP
                ALLOCATE (THERMCONDT (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)  )   ; THERMCONDT  = 1.0_WP
                MEMPC_byte = MEMPC_byte + (NCL1_io+2)*(N2DO(0)+2)*NCL3*5*8
                
                ALLOCATE ( RHS_ENERGY  (NCL1_io,N2DO(0),NCL3))          ; RHS_ENERGY = 0.0_WP
                ALLOCATE ( RHS_ENERGY0 (NCL1_io,N2DO(0),NCL3))          ; RHS_ENERGY0= 0.0_WP
                MEMPC_byte = MEMPC_byte + NCL1_io*N2DO(0)*NCL3*2*8
            END IF
            
            IF(MYID==0) CALL CHKHDL   ('       Allocating IO instantanous both thermal and flow shared variables...',MYID) 
            ALLOCATE (DENSITY    (NCL1S:NCL1E,0:N2DO(0)+1,NCL3    )  )   ; DENSITY    = 1.0_WP
            ALLOCATE (DrhoDtP    (NCL1S:NCL1E,0:N2DO(0)+1,NCL3    )  )   ; DrhoDtP    = 0.0_WP
            ALLOCATE (VISCOUSITY (NCL1S:NCL1E,0:N2DO(0)+1,NCL3    )  )   ; VISCOUSITY = 1.0_WP
            ALLOCATE (VISCOUSITY0(NCL1S:NCL1E,0:N2DO(0)+1,NCL3    )  )   ; VISCOUSITY0= 1.0_WP
            ALLOCATE (RHO_STG    (NCL1S:NCL1E,0:N2DO(0)+1,NCL3,NDV)  )   ; RHO_STG    = 1.0_WP
            ALLOCATE (MU_STG     (NCL1S:NCL1E,0:N2DO(0)+1,NCL3,NDV)  )   ; MU_STG     = 1.0_WP
            MEMPC_byte = MEMPC_byte + (NCL1_io+2)*(N2DO(0)+2)*NCL3*9*8
            
            ALLOCATE (DENSITY0  (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)      )   ; DENSITY0   = 1.0_WP
            ALLOCATE (DENSITY1  (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)      )   ; DENSITY0   = 1.0_WP
            
            IF(visthemflg==visimplicit) THEN
                !ALLOCATE (DENSITY0  (NCL1S:NCL1E,0:N2DO(0)+1,NCL3)      )   ; DENSITY0   = 1.0_WP
                ALLOCATE (DRHOI_STG (NCL1S:NCL1E,0:N2DO(0)+1,NCL3,NDV)  )   ; DRHOI_STG  = 0.0_WP
                MEMPC_byte = MEMPC_byte + (NCL1_io+2)*(N2DO(0)+2)*NCL3*4*8
            END IF

            IF(TGFLOWFLG) THEN
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO z averaged variables...',MYID) 
                ALLOCATE ( IPV_io(0:NCL1_io   ) )                       ; IPV_io = 0
                ALLOCATE ( IMV_io(1:NCL1_io+1 ) )                       ; IMV_io = 0
                MEMPC_byte = MEMPC_byte + 2*NND1_io*4
            
                !============FOR POSTPROCESS============================
                ALLOCATE ( U1zL_io( NCL1S:NCL1E,N2DO(0),NDV+1 ) )                               ; U1zL_io = 0.0_WP
                ALLOCATE ( G1zL_io( NCL1S:NCL1E,N2DO(0),NDV   ) )                               ; G1zL_io = 0.0_WP
                ALLOCATE ( UPzL_io( NCL1S:NCL1E,N2DO(0),NDV   ) )                               ; UPzL_io = 0.0_WP
                
                ALLOCATE ( U2zL_io( NCL1S:NCL1E,N2DO(0),NDV*(7-NDV)/2+NDV-3                ) )  ; U2zL_io = 0.0_WP
                ALLOCATE ( UGzL_io( NCL1S:NCL1E,N2DO(0),NDV*(7-NDV)/2+NDV-3                ) )  ; UGzL_io = 0.0_WP
                ALLOCATE ( UGUzL_io(NCL1S:NCL1E,N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) )  ; UGUzL_io= 0.0_WP
                ALLOCATE ( U3zL_io (NCL1S:NCL1E,N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) )  ; U3zL_io = 0.0_WP
            
                ALLOCATE ( DVDL1zL_io( NCL1S:NCL1E,N2DO(0), NDV,     NDV      ) )               ; DVDL1zL_io= 0.0_WP
                ALLOCATE ( DVDLPzL_io( NCL1S:NCL1E,N2DO(0), NDV,     NDV      ) )               ; DVDLPzL_io= 0.0_WP
                ALLOCATE ( DVDL2zL_io( NCL1S:NCL1E,N2DO(0), NDV*NDV, NDV*NDV  ) )               ; DVDL2zL_io= 0.0_WP
                
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO z-t averaged variables...',MYID) 
                ALLOCATE ( U1ztL_io( NCL1S:NCL1E,N2DO(0),NDV+1 ) )                              ; U1ztL_io = 0.0_WP
                ALLOCATE ( G1ztL_io( NCL1S:NCL1E,N2DO(0),NDV   ) )                              ; G1ztL_io = 0.0_WP
                ALLOCATE ( UPztL_io( NCL1S:NCL1E,N2DO(0),NDV  ) )                               ; UPztL_io = 0.0_WP
                
                ALLOCATE ( U2ztL_io( NCL1S:NCL1E,N2DO(0),NDV*(7-NDV)/2+NDV-3                ) ) ; U2ztL_io = 0.0_WP
                ALLOCATE ( UGztL_io( NCL1S:NCL1E,N2DO(0),NDV*(7-NDV)/2+NDV-3                ) ) ; UGztL_io = 0.0_WP
                ALLOCATE ( UGUztL_io(NCL1S:NCL1E,N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) ) ; UGUztL_io= 0.0_WP
                ALLOCATE ( U3ztL_io (NCL1S:NCL1E,N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) ) ; U3ztL_io = 0.0_WP
            
                ALLOCATE ( DVDL1ztL_io( NCL1S:NCL1E,N2DO(0), NDV, NDV  ) )                      ; DVDL1ztL_io= 0.0_WP
                ALLOCATE ( DVDLPztL_io( NCL1S:NCL1E,N2DO(0), NDV, NDV  ) )                      ; DVDLPztL_io= 0.0_WP
                ALLOCATE ( DVDL2ztL_io( NCL1S:NCL1E,N2DO(0), NDV*NDV, NDV*NDV ) )               ; DVDL2ztL_io= 0.0_WP
            
                MEMPC_byte = MEMPC_byte + NCL1_IO*N2DO(0)* ( &
                            NDV+1 + NDV + NDV + NDV*(7-NDV)/2+NDV-3 + NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8 + &
                            NDV*NDV*2 + (NDV*(7-NDV)/2+NDV-3)*NDV ) * 8
                            
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO z averaged quadrant variables...',MYID) 
                ALLOCATE ( QUADUVzL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADUVzL_IO =0.0_wp
                ALLOCATE ( QUADVzzL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADVzzL_IO =0.0_wp
                ALLOCATE ( QUADTKzL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADTKzL_IO =0.0_wp
                ALLOCATE ( QUADDRzL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADDRzL_IO =0.0_wp
                ALLOCATE ( QUADDUV1zL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                       ; QUADDUV1zL_IO =0.0_wp
                ALLOCATE ( QUADDUV2zL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                       ; QUADDUV2zL_IO =0.0_wp
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO z-t averaged quadrant variables...',MYID) 
                ALLOCATE ( QUADUVztL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADUVztL_IO =0.0_wp
                ALLOCATE ( QUADVzztL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADVzztL_IO =0.0_wp
                ALLOCATE ( QUADTKztL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADTKztL_IO =0.0_wp
                ALLOCATE ( QUADDRztL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                         ; QUADDRztL_IO =0.0_wp
                ALLOCATE ( QUADDUV1ztL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                       ; QUADDUV1ztL_IO =0.0_wp
                ALLOCATE ( QUADDUV2ztL_IO (NCL1S:NCL1E,N2DO(0),4,QUADHN) )                       ; QUADDUV2ztL_IO =0.0_wp
                            
                            
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO z averaged density and viscosity variables...',MYID)             
                ALLOCATE ( D1zL_io ( NCL1S:NCL1E,N2DO(0) ) )                 ; D1zL_io  = 1.0_WP
                ALLOCATE ( D1ztL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; D1ztL_io = 1.0_WP      
                ALLOCATE ( D2zL_io ( NCL1S:NCL1E,N2DO(0) ) )                 ; D2zL_io  = 1.0_WP   
                ALLOCATE ( D2ztL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; D2ztL_io = 1.0_WP   
                
                ALLOCATE ( M1zL_io ( NCL1S:NCL1E,N2DO(0) ) )                 ; M1zL_io  = 1.0_WP
                ALLOCATE ( M1ztL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; M1ztL_io = 1.0_WP      
                
                IF(thermlflg==1) THEN       
                    IF(MYID==0) CALL CHKHDL   ('       Allocating IO z averaged thermal variables...',MYID)            
                    ALLOCATE ( T1zL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; T1zL_io = 0.0_WP
                    ALLOCATE ( H1zL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; H1zL_io = 0.0_WP
                         
                    ALLOCATE ( T2zL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; T2zL_io = 0.0_WP
                    ALLOCATE ( H2zL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; H2zL_io = 0.0_WP
                    
                    ALLOCATE ( DHzL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; DHzL_io = 0.0_WP
                    ALLOCATE ( PHzL_io( NCL1S:NCL1E,N2DO(0) ) )                 ; PHzL_io = 0.0_WP
                    
                    ALLOCATE ( DVDL1MzL_io  (NCL1S:NCL1E,N2DO(0), NDV, NDV) )     ; DVDL1MzL_io  = 0.0_wp
                    ALLOCATE ( DVDL2MzL_io  (NCL1S:NCL1E,N2DO(0), (NDV-1)*3+NDV, (NDV-1)*3+NDV) )     ; DVDL2MzL_io  = 0.0_wp
                    ALLOCATE ( DVDL1MHzL_io (NCL1S:NCL1E,N2DO(0), NDV, NDV) )     ; DVDL1MHzL_io = 0.0_wp
                    ALLOCATE ( DVDL1MUzL_io (NCL1S:NCL1E,N2DO(0), NDV, NDV, NDV)) ; DVDL1MUzL_io = 0.0_wp
                    
                        
                    ALLOCATE ( UHzL_io( NCL1S:NCL1E,N2DO(0),NDV ) )             ; UHzL_io = 0.0_WP
                    ALLOCATE ( GHzL_io( NCL1S:NCL1E,N2DO(0),NDV ) )             ; GHzL_io = 0.0_WP
                    ALLOCATE ( U2DHzL_io( NCL1S:NCL1E,N2DO(0),(NDV*(7-NDV))/2+NDV-3 ) ) ; U2DHzL_io = 0.0_WP
                    
                    ALLOCATE ( DhDL1zL_io(NCL1S:NCL1E,N2DO(0), NDV) )               ; DhDL1zL_io = 0.0_wp
                    ALLOCATE ( DhDLPzL_io(NCL1S:NCL1E,N2DO(0), NDV) )               ; DhDLPzL_io = 0.0_wp
                    ALLOCATE ( DTDLKzL_io(NCL1S:NCL1E,N2DO(0), NDV) )               ; DTDLKzL_io = 0.0_wp
                    ALLOCATE ( DTDLKUzL_io(NCL1S:NCL1E,N2DO(0), NDV, NDV) )         ; DTDLKUzL_io = 0.0_wp
                    ALLOCATE ( DTDLKDVDLzL_io(NCL1S:NCL1E,N2DO(0), NDV, NDV, NDV) ) ; DTDLKDVDLzL_io = 0.0_wp
                    ALLOCATE ( DhDLMDVDLzL_io(NCL1S:NCL1E,N2DO(0), NDV, NDV, NDV) ) ; DhDLMDVDLzL_io = 0.0_wp
                    
                    
                    IF(MYID==0) CALL CHKHDL   ('       Allocating IO z-t averaged thermal variables...',MYID)            
                    ALLOCATE ( T1ztL_io( NCL1S:NCL1E,N2DO(0) ) )                ; T1ztL_io = 0.0_WP
                    ALLOCATE ( H1ztL_io( NCL1S:NCL1E,N2DO(0) ) )                ; H1ztL_io = 0.0_WP
                         
                    ALLOCATE ( T2ztL_io( NCL1S:NCL1E,N2DO(0) ) )                ; T2ztL_io = 0.0_WP
                    ALLOCATE ( H2ztL_io( NCL1S:NCL1E,N2DO(0) ) )                ; H2ztL_io = 0.0_WP
                    
                    ALLOCATE ( DHztL_io( NCL1S:NCL1E,N2DO(0) ) )                ; DHztL_io = 0.0_WP
                    ALLOCATE ( PHztL_io( NCL1S:NCL1E,N2DO(0) ) )                ; PHztL_io = 0.0_WP
                        
                    ALLOCATE ( DVDL1MztL_io  (NCL1S:NCL1E,N2DO(0), NDV, NDV) )     ; DVDL1MztL_io  = 0.0_wp
                    ALLOCATE ( DVDL2MztL_io  (NCL1S:NCL1E,N2DO(0), (NDV-1)*3+NDV, (NDV-1)*3+NDV))    ; DVDL2MztL_io  = 0.0_wp
                    ALLOCATE ( DVDL1MHztL_io (NCL1S:NCL1E,N2DO(0), NDV, NDV) )     ; DVDL1MHztL_io = 0.0_wp
                    ALLOCATE ( DVDL1MUztL_io (NCL1S:NCL1E,N2DO(0), NDV, NDV, NDV)) ; DVDL1MUztL_io = 0.0_wp
                        
                    ALLOCATE ( UHztL_io  ( NCL1S:NCL1E,N2DO(0),NDV ) )             ; UHztL_io   = 0.0_WP
                    ALLOCATE ( GHztL_io  ( NCL1S:NCL1E,N2DO(0),NDV ) )             ; GHztL_io   = 0.0_WP
                    ALLOCATE ( U2DHztL_io( NCL1S:NCL1E,N2DO(0),(NDV*(7-NDV))/2+NDV-3 ) ) ; U2DHztL_io = 0.0_WP
                    
                    ALLOCATE ( DhDL1ztL_io(NCL1S:NCL1E,N2DO(0), NDV) )               ; DhDL1ztL_io = 0.0_wp
                    ALLOCATE ( DhDLPztL_io(NCL1S:NCL1E,N2DO(0), NDV) )               ; DhDLPztL_io = 0.0_wp
                    ALLOCATE ( DTDLKztL_io(NCL1S:NCL1E,N2DO(0), NDV) )               ; DTDLKztL_io = 0.0_wp
                    ALLOCATE ( DTDLKUztL_io(NCL1S:NCL1E,N2DO(0), NDV, NDV) )         ; DTDLKUztL_io = 0.0_wp
                    ALLOCATE ( DTDLKDVDLztL_io(NCL1S:NCL1E,N2DO(0), NDV, NDV, NDV) ) ; DTDLKDVDLztL_io = 0.0_wp
                    ALLOCATE ( DhDLMDVDLztL_io(NCL1S:NCL1E,N2DO(0), NDV, NDV, NDV) ) ; DhDLMDVDLztL_io = 0.0_wp
                    
                END IF
                
            ELSE
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z averaged variables...',MYID) 
                ALLOCATE ( IPV_io(NCL1_io   ) )                       ; IPV_io = 0
                ALLOCATE ( IMV_io(NCL1_io   ) )                       ; IMV_io = 0
                MEMPC_byte = MEMPC_byte + 2*NCL1_io*4
                
                
                !============FOR POSTPROCESS============================
                ALLOCATE ( U1xzL_io( N2DO(0),NDV+1 ) )                               ; U1xzL_io = 0.0_WP
                ALLOCATE ( G1xzL_io( N2DO(0),NDV   ) )                               ; G1xzL_io = 0.0_WP
                ALLOCATE ( UPxzL_io( N2DO(0),NDV  ) )                                ; UPxzL_io = 0.0_WP
                
                ALLOCATE ( U2xzL_io( N2DO(0),NDV*(7-NDV)/2+NDV-3                ) )  ; U2xzL_io = 0.0_WP
                ALLOCATE ( UGxzL_io( N2DO(0),NDV*(7-NDV)/2+NDV-3                ) )  ; UGxzL_io = 0.0_WP
                ALLOCATE ( UGUxzL_io(N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) )  ; UGUxzL_io= 0.0_WP
                ALLOCATE ( U3xzL_io( N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) )  ; U3xzL_io = 0.0_WP
            
                ALLOCATE ( DVDL1xzL_io( N2DO(0), NDV, NDV  )         )               ; DVDL1xzL_io= 0.0_WP
                ALLOCATE ( DVDLPxzL_io( N2DO(0), NDV, NDV  )         )               ; DVDLPxzL_io= 0.0_WP
                ALLOCATE ( DVDL2xzL_io( N2DO(0), NDV*NDV, NDV*NDV  ) )               ; DVDL2xzL_io= 0.0_WP
                
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z-t averaged variables...',MYID) 
                ALLOCATE ( U1xztL_io( N2DO(0),NDV+1 ) )                              ; U1xztL_io = 0.0_WP
                ALLOCATE ( G1xztL_io( N2DO(0),NDV   ) )                              ; G1xztL_io = 0.0_WP
                ALLOCATE ( UPxztL_io( N2DO(0),NDV  ) )                               ; UPxztL_io = 0.0_WP
                
                ALLOCATE ( U2xztL_io( N2DO(0),NDV*(7-NDV)/2+NDV-3                ) ) ; U2xztL_io = 0.0_WP
                ALLOCATE ( UGxztL_io( N2DO(0),NDV*(7-NDV)/2+NDV-3                ) ) ; UGxztL_io = 0.0_WP
                ALLOCATE ( UGUxztL_io(N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) ) ; UGUxztL_io= 0.0_WP
                ALLOCATE ( U3xztL_io (N2DO(0),NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8  ) ) ; UGUxztL_io= 0.0_WP
            
                ALLOCATE ( DVDL1xztL_io( N2DO(0), NDV, NDV  ) )                      ; DVDL1xztL_io= 0.0_WP
                ALLOCATE ( DVDLPxztL_io( N2DO(0), NDV, NDV  ) )                      ; DVDLPxztL_io= 0.0_WP
                ALLOCATE ( DVDL2xztL_io( N2DO(0), NDV*NDV, NDV*NDV ) )               ; DVDL2xztL_io= 0.0_WP
                
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z averaged quadrant...',MYID) 
                ALLOCATE ( QUADUVxzL_IO (N2DO(0),4,QUADHN) )                         ; QUADUVxzL_IO =0.0_wp
                ALLOCATE ( QUADVzxzL_IO (N2DO(0),4,QUADHN) )                         ; QUADVzxzL_IO =0.0_wp
                ALLOCATE ( QUADTKxzL_IO (N2DO(0),4,QUADHN) )                         ; QUADTKxzL_IO =0.0_wp
                ALLOCATE ( QUADDRxzL_IO (N2DO(0),4,QUADHN) )                         ; QUADDRxzL_IO =0.0_wp
                ALLOCATE ( QUADDUV1xzL_IO (N2DO(0),4,QUADHN) )                       ; QUADDUV1xzL_IO =0.0_wp
                ALLOCATE ( QUADDUV2xzL_IO (N2DO(0),4,QUADHN) )                       ; QUADDUV2xzL_IO =0.0_wp
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z-t averaged quadrant...',MYID) 
                ALLOCATE ( QUADUVxztL_IO (N2DO(0),4,QUADHN) )                         ; QUADUVxztL_IO =0.0_wp
                ALLOCATE ( QUADVzxztL_IO (N2DO(0),4,QUADHN) )                         ; QUADVzxztL_IO =0.0_wp
                ALLOCATE ( QUADTKxztL_IO (N2DO(0),4,QUADHN) )                         ; QUADTKxztL_IO =0.0_wp
                ALLOCATE ( QUADDRxztL_IO (N2DO(0),4,QUADHN) )                         ; QUADDRxztL_IO =0.0_wp
                ALLOCATE ( QUADDUV1xztL_IO (N2DO(0),4,QUADHN) )                       ; QUADDUV1xztL_IO =0.0_wp
                ALLOCATE ( QUADDUV2xztL_IO (N2DO(0),4,QUADHN) )                       ; QUADDUV2xztL_IO =0.0_wp
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z averaged octant...',MYID) 
                ALLOCATE ( OCTTUVxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTTUVxzL_IO =0.0_wp
                ALLOCATE ( OCTTVzxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTTVzxzL_IO =0.0_wp
                ALLOCATE ( OCTTTKxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTTTKxzL_IO =0.0_wp
                ALLOCATE ( OCTTDRxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTTDRxzL_IO =0.0_wp
                ALLOCATE ( OCTTDUV1xzL_IO (N2DO(0),8,QUADHN) )                       ; OCTTDUV1xzL_IO =0.0_wp
                ALLOCATE ( OCTTDUV2xzL_IO (N2DO(0),8,QUADHN) )                       ; OCTTDUV2xzL_IO =0.0_wp
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z-t averaged octant...',MYID) 
                ALLOCATE ( OCTTUVxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTTUVxztL_IO =0.0_wp
                ALLOCATE ( OCTTVzxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTTVzxztL_IO =0.0_wp
                ALLOCATE ( OCTTTKxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTTTKxztL_IO =0.0_wp
                ALLOCATE ( OCTTDRxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTTDRxztL_IO =0.0_wp
                ALLOCATE ( OCTTDUV1xztL_IO (N2DO(0),8,QUADHN) )                       ; OCTTDUV1xztL_IO =0.0_wp
                ALLOCATE ( OCTTDUV2xztL_IO (N2DO(0),8,QUADHN) )                       ; OCTTDUV2xztL_IO =0.0_wp
                
                ALLOCATE ( OCTDUVxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTDUVxzL_IO =0.0_wp
                ALLOCATE ( OCTDVzxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTDVzxzL_IO =0.0_wp
                ALLOCATE ( OCTDTKxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTDTKxzL_IO =0.0_wp
                ALLOCATE ( OCTDDRxzL_IO (N2DO(0),8,QUADHN) )                         ; OCTDDRxzL_IO =0.0_wp
                ALLOCATE ( OCTDDUV1xzL_IO (N2DO(0),8,QUADHN) )                       ; OCTDDUV1xzL_IO =0.0_wp
                ALLOCATE ( OCTDDUV2xzL_IO (N2DO(0),8,QUADHN) )                       ; OCTDDUV2xzL_IO =0.0_wp
                
                ALLOCATE ( OCTDUVxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTDUVxztL_IO =0.0_wp
                ALLOCATE ( OCTDVzxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTDVzxztL_IO =0.0_wp
                ALLOCATE ( OCTDTKxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTDTKxztL_IO =0.0_wp
                ALLOCATE ( OCTDDRxztL_IO (N2DO(0),8,QUADHN) )                         ; OCTDDRxztL_IO =0.0_wp
                ALLOCATE ( OCTDDUV1xztL_IO (N2DO(0),8,QUADHN) )                       ; OCTDDUV1xztL_IO =0.0_wp
                ALLOCATE ( OCTDDUV2xztL_IO (N2DO(0),8,QUADHN) )                       ; OCTDDUV2xztL_IO =0.0_wp
                
                
                ALLOCATE ( FUxzL_IO(N2DO(0),NDV+1)  ) ; FUxzL_IO=0.0_WP
                ALLOCATE ( FUxztL_IO(N2DO(0),NDV+1) ) ; FUxztL_IO=0.0_WP
            
                MEMPC_byte = MEMPC_byte + N2DO(0)* ( &
                            NDV+1 + NDV + NDV + NDV*(7-NDV)/2+NDV-3 + NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8 + &
                            NDV*NDV*2 + (NDV*(7-NDV)/2+NDV-3)*NDV + 4*4*QUADHN) * 8
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z averaged density and viscosity variables...',MYID)                     
                ALLOCATE ( D1xzL_io ( N2DO(0) ) )                    ; D1xzL_io  = 1.0_WP
                ALLOCATE ( D1xztL_io( N2DO(0) ) )                    ; D1xztL_io = 1.0_WP
                ALLOCATE ( M1xzL_io ( N2DO(0) ) )                    ; M1xzL_io  = 1.0_WP
                ALLOCATE ( M1xztL_io( N2DO(0) ) )                    ; M1xztL_io = 1.0_WP
                
                !========CORRELATION===========================
                IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z averaged correlation and spectra...',MYID)         
                ! IN PHYSICAL SPACE
                ! CORRELATIONo R(IJ)X(K) = Velocity_i and Velocity_J along the direction K
                ALLOCATE ( R11X1_xzLa (NCL2,NCL1_io,2) ) ; R11X1_xzLa = 0.0_WP 
                ALLOCATE ( R22X1_xzLa (NCL2,NCL1_io,2) ) ; R22X1_xzLa = 0.0_WP 
                ALLOCATE ( R33X1_xzLa (NCL2,NCL1_io,2) ) ; R33X1_xzLa = 0.0_WP 
                ALLOCATE ( R12X1_xzLa (NCL2,NCL1_io,2) ) ; R12X1_xzLa = 0.0_WP 
                ALLOCATE ( R13X1_xzLa (NCL2,NCL1_io,2) ) ; R13X1_xzLa = 0.0_WP 
                ALLOCATE ( R23X1_xzLa (NCL2,NCL1_io,2) ) ; R23X1_xzLa = 0.0_WP 
                
                ALLOCATE ( R11X3_xzLa (NCL2,NCL3,2) ) ; R11X3_xzLa = 0.0_WP 
                ALLOCATE ( R22X3_xzLa (NCL2,NCL3,2) ) ; R22X3_xzLa = 0.0_WP 
                ALLOCATE ( R33X3_xzLa (NCL2,NCL3,2) ) ; R33X3_xzLa = 0.0_WP 
                ALLOCATE ( R12X3_xzLa (NCL2,NCL3,2) ) ; R12X3_xzLa = 0.0_WP 
                ALLOCATE ( R13X3_xzLa (NCL2,NCL3,2) ) ; R13X3_xzLa = 0.0_WP 
                ALLOCATE ( R23X3_xzLa (NCL2,NCL3,2) ) ; R23X3_xzLa = 0.0_WP 
                
                !===============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENE11T_xzLa (NCL2,NCL1_io,2) ) ; ENE11T_xzLa =0.0_WP 
                ALLOCATE ( ENE22T_xzLa (NCL2,NCL1_io,2) ) ; ENE22T_xzLa =0.0_WP 
                ALLOCATE ( ENE33T_xzLa (NCL2,NCL1_io,2) ) ; ENE33T_xzLa =0.0_WP 
                ALLOCATE ( ENE12T_xzLa (NCL2,NCL1_io,2) ) ; ENE12T_xzLa =0.0_WP 
                ALLOCATE ( ENE13T_xzLa (NCL2,NCL1_io,2) ) ; ENE13T_xzLa =0.0_WP 
                ALLOCATE ( ENE23T_xzLa (NCL2,NCL1_io,2) ) ; ENE23T_xzLa =0.0_WP 
                
                !=============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENE11Z_xzLa (NCL2,NCL3,2) ) ; ENE11Z_xzLa =0.0_WP 
                ALLOCATE ( ENE22Z_xzLa (NCL2,NCL3,2) ) ; ENE22Z_xzLa =0.0_WP 
                ALLOCATE ( ENE33Z_xzLa (NCL2,NCL3,2) ) ; ENE33Z_xzLa =0.0_WP  
                ALLOCATE ( ENE12Z_xzLa (NCL2,NCL3,2) ) ; ENE12Z_xzLa =0.0_WP 
                ALLOCATE ( ENE13Z_xzLa (NCL2,NCL3,2) ) ; ENE13Z_xzLa =0.0_WP 
                ALLOCATE ( ENE23Z_xzLa (NCL2,NCL3,2) ) ; ENE23Z_xzLa =0.0_WP 
        
                !========CORRELATION===========================
                ! IN PHYSICAL SPACE
                ! CORRELATIONo R(IJ)X(K) = voriticity_i and voriticity_J along the direction K
                ALLOCATE ( V11X1_xzLa (NCL2,NCL1_io,2) ) ; V11X1_xzLa = 0.0_WP 
                ALLOCATE ( V22X1_xzLa (NCL2,NCL1_io,2) ) ; V22X1_xzLa = 0.0_WP 
                ALLOCATE ( V33X1_xzLa (NCL2,NCL1_io,2) ) ; V33X1_xzLa = 0.0_WP 
                ALLOCATE ( V12X1_xzLa (NCL2,NCL1_io,2) ) ; V12X1_xzLa = 0.0_WP 
                ALLOCATE ( V13X1_xzLa (NCL2,NCL1_io,2) ) ; V13X1_xzLa = 0.0_WP 
                ALLOCATE ( V23X1_xzLa (NCL2,NCL1_io,2) ) ; V23X1_xzLa = 0.0_WP 
                
                ALLOCATE ( V11X3_xzLa (NCL2,NCL3,2) ) ; V11X3_xzLa = 0.0_WP 
                ALLOCATE ( V22X3_xzLa (NCL2,NCL3,2) ) ; V22X3_xzLa = 0.0_WP 
                ALLOCATE ( V33X3_xzLa (NCL2,NCL3,2) ) ; V33X3_xzLa = 0.0_WP 
                ALLOCATE ( V12X3_xzLa (NCL2,NCL3,2) ) ; V12X3_xzLa = 0.0_WP 
                ALLOCATE ( V13X3_xzLa (NCL2,NCL3,2) ) ; V13X3_xzLa = 0.0_WP 
                ALLOCATE ( V23X3_xzLa (NCL2,NCL3,2) ) ; V23X3_xzLa = 0.0_WP 
                
                !===============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENV11T_xzLa (NCL2,NCL1_io,2) ) ; ENV11T_xzLa =0.0_WP 
                ALLOCATE ( ENV22T_xzLa (NCL2,NCL1_io,2) ) ; ENV22T_xzLa =0.0_WP 
                ALLOCATE ( ENV33T_xzLa (NCL2,NCL1_io,2) ) ; ENV33T_xzLa =0.0_WP 
                ALLOCATE ( ENV12T_xzLa (NCL2,NCL1_io,2) ) ; ENV12T_xzLa =0.0_WP 
                ALLOCATE ( ENV13T_xzLa (NCL2,NCL1_io,2) ) ; ENV13T_xzLa =0.0_WP 
                ALLOCATE ( ENV23T_xzLa (NCL2,NCL1_io,2) ) ; ENV23T_xzLa =0.0_WP 
                
                !=============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENV11Z_xzLa (NCL2,NCL3,2) ) ; ENV11Z_xzLa =0.0_WP 
                ALLOCATE ( ENV22Z_xzLa (NCL2,NCL3,2) ) ; ENV22Z_xzLa =0.0_WP 
                ALLOCATE ( ENV33Z_xzLa (NCL2,NCL3,2) ) ; ENV33Z_xzLa =0.0_WP  
                ALLOCATE ( ENV12Z_xzLa (NCL2,NCL3,2) ) ; ENV12Z_xzLa =0.0_WP 
                ALLOCATE ( ENV13Z_xzLa (NCL2,NCL3,2) ) ; ENV13Z_xzLa =0.0_WP 
                ALLOCATE ( ENV23Z_xzLa (NCL2,NCL3,2) ) ; ENV23Z_xzLa =0.0_WP 
                
                !=============velocity and voriticity==============
                ALLOCATE ( VO11X1_xzLa (NCL2,NCL1_io,2) ) ; VO11X1_xzLa = 0.0_WP
                ALLOCATE ( VO12X1_xzLa (NCL2,NCL1_io,2) ) ; VO12X1_xzLa = 0.0_WP
                ALLOCATE ( VO13X1_xzLa (NCL2,NCL1_io,2) ) ; VO13X1_xzLa = 0.0_WP
                
                ALLOCATE ( VO21X1_xzLa (NCL2,NCL1_io,2) ) ; VO21X1_xzLa = 0.0_WP
                ALLOCATE ( VO22X1_xzLa (NCL2,NCL1_io,2) ) ; VO22X1_xzLa = 0.0_WP
                ALLOCATE ( VO23X1_xzLa (NCL2,NCL1_io,2) ) ; VO23X1_xzLa = 0.0_WP
                
                ALLOCATE ( VO31X1_xzLa (NCL2,NCL1_io,2) ) ; VO31X1_xzLa = 0.0_WP
                ALLOCATE ( VO32X1_xzLa (NCL2,NCL1_io,2) ) ; VO32X1_xzLa = 0.0_WP
                ALLOCATE ( VO33X1_xzLa (NCL2,NCL1_io,2) ) ; VO33X1_xzLa = 0.0_WP
                
                ALLOCATE ( VO11X3_xzLa (NCL2,NCL3,2) ) ; VO11X3_xzLa = 0.0_WP
                ALLOCATE ( VO12X3_xzLa (NCL2,NCL3,2) ) ; VO12X3_xzLa = 0.0_WP
                ALLOCATE ( VO13X3_xzLa (NCL2,NCL3,2) ) ; VO13X3_xzLa = 0.0_WP
                
                ALLOCATE ( VO21X3_xzLa (NCL2,NCL3,2) ) ; VO21X3_xzLa = 0.0_WP
                ALLOCATE ( VO22X3_xzLa (NCL2,NCL3,2) ) ; VO22X3_xzLa = 0.0_WP
                ALLOCATE ( VO23X3_xzLa (NCL2,NCL3,2) ) ; VO23X3_xzLa = 0.0_WP
                
                ALLOCATE ( VO31X3_xzLa (NCL2,NCL3,2) ) ; VO31X3_xzLa = 0.0_WP
                ALLOCATE ( VO32X3_xzLa (NCL2,NCL3,2) ) ; VO32X3_xzLa = 0.0_WP
                ALLOCATE ( VO33X3_xzLa (NCL2,NCL3,2) ) ; VO33X3_xzLa = 0.0_WP
                
                ALLOCATE ( EVO11T_xzLa (NCL2,NCL1_io,2) ) ; EVO11T_xzLa = 0.0_WP
                ALLOCATE ( EVO12T_xzLa (NCL2,NCL1_io,2) ) ; EVO12T_xzLa = 0.0_WP
                ALLOCATE ( EVO13T_xzLa (NCL2,NCL1_io,2) ) ; EVO13T_xzLa = 0.0_WP
                
                ALLOCATE ( EVO21T_xzLa (NCL2,NCL1_io,2) ) ; EVO21T_xzLa = 0.0_WP
                ALLOCATE ( EVO22T_xzLa (NCL2,NCL1_io,2) ) ; EVO22T_xzLa = 0.0_WP
                ALLOCATE ( EVO23T_xzLa (NCL2,NCL1_io,2) ) ; EVO23T_xzLa = 0.0_WP
                
                ALLOCATE ( EVO31T_xzLa (NCL2,NCL1_io,2) ) ; EVO31T_xzLa = 0.0_WP
                ALLOCATE ( EVO32T_xzLa (NCL2,NCL1_io,2) ) ; EVO32T_xzLa = 0.0_WP
                ALLOCATE ( EVO33T_xzLa (NCL2,NCL1_io,2) ) ; EVO33T_xzLa = 0.0_WP
                
                ALLOCATE ( EVO11Z_xzLa (NCL2,NCL3,2) ) ; EVO11Z_xzLa = 0.0_WP
                ALLOCATE ( EVO12Z_xzLa (NCL2,NCL3,2) ) ; EVO12Z_xzLa = 0.0_WP
                ALLOCATE ( EVO13Z_xzLa (NCL2,NCL3,2) ) ; EVO13Z_xzLa = 0.0_WP
                
                ALLOCATE ( EVO21Z_xzLa (NCL2,NCL3,2) ) ; EVO21Z_xzLa = 0.0_WP
                ALLOCATE ( EVO22Z_xzLa (NCL2,NCL3,2) ) ; EVO22Z_xzLa = 0.0_WP
                ALLOCATE ( EVO23Z_xzLa (NCL2,NCL3,2) ) ; EVO23Z_xzLa = 0.0_WP
                
                ALLOCATE ( EVO31Z_xzLa (NCL2,NCL3,2) ) ; EVO31Z_xzLa = 0.0_WP
                ALLOCATE ( EVO32Z_xzLa (NCL2,NCL3,2) ) ; EVO32Z_xzLa = 0.0_WP
                ALLOCATE ( EVO33Z_xzLa (NCL2,NCL3,2) ) ; EVO33Z_xzLa = 0.0_WP
                
                ! IN PHYSICAL SPACE
                ! CORRELATIONo R(IJ)X(K) = Velocity_i and Velocity_J along the direction K
                ALLOCATE ( R11X1_xztLa (NCL2,NCL1_io,2) ) ; R11X1_xztLa = 0.0_WP 
                ALLOCATE ( R22X1_xztLa (NCL2,NCL1_io,2) ) ; R22X1_xztLa = 0.0_WP 
                ALLOCATE ( R33X1_xztLa (NCL2,NCL1_io,2) ) ; R33X1_xztLa = 0.0_WP 
                ALLOCATE ( R12X1_xztLa (NCL2,NCL1_io,2) ) ; R12X1_xztLa = 0.0_WP 
                ALLOCATE ( R13X1_xztLa (NCL2,NCL1_io,2) ) ; R13X1_xztLa = 0.0_WP 
                ALLOCATE ( R23X1_xztLa (NCL2,NCL1_io,2) ) ; R23X1_xztLa = 0.0_WP 
                
                ALLOCATE ( R11X3_xztLa (NCL2,NCL3,2) ) ; R11X3_xztLa = 0.0_WP 
                ALLOCATE ( R22X3_xztLa (NCL2,NCL3,2) ) ; R22X3_xztLa = 0.0_WP 
                ALLOCATE ( R33X3_xztLa (NCL2,NCL3,2) ) ; R33X3_xztLa = 0.0_WP 
                ALLOCATE ( R12X3_xztLa (NCL2,NCL3,2) ) ; R12X3_xztLa = 0.0_WP 
                ALLOCATE ( R13X3_xztLa (NCL2,NCL3,2) ) ; R13X3_xztLa = 0.0_WP 
                ALLOCATE ( R23X3_xztLa (NCL2,NCL3,2) ) ; R23X3_xztLa = 0.0_WP 
                
                !===============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENE11T_xztLa (NCL2,NCL1_io,2) ) ; ENE11T_xztLa =0.0_WP 
                ALLOCATE ( ENE22T_xztLa (NCL2,NCL1_io,2) ) ; ENE22T_xztLa =0.0_WP 
                ALLOCATE ( ENE33T_xztLa (NCL2,NCL1_io,2) ) ; ENE33T_xztLa =0.0_WP 
                ALLOCATE ( ENE12T_xztLa (NCL2,NCL1_io,2) ) ; ENE12T_xztLa =0.0_WP 
                ALLOCATE ( ENE13T_xztLa (NCL2,NCL1_io,2) ) ; ENE13T_xztLa =0.0_WP 
                ALLOCATE ( ENE23T_xztLa (NCL2,NCL1_io,2) ) ; ENE23T_xztLa =0.0_WP 
                
                !=============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENE11Z_xztLa (NCL2,NCL3,2) ) ; ENE11Z_xztLa =0.0_WP 
                ALLOCATE ( ENE22Z_xztLa (NCL2,NCL3,2) ) ; ENE22Z_xztLa =0.0_WP 
                ALLOCATE ( ENE33Z_xztLa (NCL2,NCL3,2) ) ; ENE33Z_xztLa =0.0_WP  
                ALLOCATE ( ENE12Z_xztLa (NCL2,NCL3,2) ) ; ENE12Z_xztLa =0.0_WP 
                ALLOCATE ( ENE13Z_xztLa (NCL2,NCL3,2) ) ; ENE13Z_xztLa =0.0_WP 
                ALLOCATE ( ENE23Z_xztLa (NCL2,NCL3,2) ) ; ENE23Z_xztLa =0.0_WP 
        
                !========CORRELATION===========================
                ! IN PHYSICAL SPACE
                ! CORRELATIONo R(IJ)X(K) = voriticity_i and voriticity_J along the direction K
                ALLOCATE ( V11X1_xztLa (NCL2,NCL1_io,2) ) ; V11X1_xztLa = 0.0_WP 
                ALLOCATE ( V22X1_xztLa (NCL2,NCL1_io,2) ) ; V22X1_xztLa = 0.0_WP 
                ALLOCATE ( V33X1_xztLa (NCL2,NCL1_io,2) ) ; V33X1_xztLa = 0.0_WP 
                ALLOCATE ( V12X1_xztLa (NCL2,NCL1_io,2) ) ; V12X1_xztLa = 0.0_WP 
                ALLOCATE ( V13X1_xztLa (NCL2,NCL1_io,2) ) ; V13X1_xztLa = 0.0_WP 
                ALLOCATE ( V23X1_xztLa (NCL2,NCL1_io,2) ) ; V23X1_xztLa = 0.0_WP 
                
                ALLOCATE ( V11X3_xztLa (NCL2,NCL3,2) ) ; V11X3_xztLa = 0.0_WP 
                ALLOCATE ( V22X3_xztLa (NCL2,NCL3,2) ) ; V22X3_xztLa = 0.0_WP 
                ALLOCATE ( V33X3_xztLa (NCL2,NCL3,2) ) ; V33X3_xztLa = 0.0_WP 
                ALLOCATE ( V12X3_xztLa (NCL2,NCL3,2) ) ; V12X3_xztLa = 0.0_WP 
                ALLOCATE ( V13X3_xztLa (NCL2,NCL3,2) ) ; V13X3_xztLa = 0.0_WP 
                ALLOCATE ( V23X3_xztLa (NCL2,NCL3,2) ) ; V23X3_xztLa = 0.0_WP 
                
                !===============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENV11T_xztLa (NCL2,NCL1_io,2) ) ; ENV11T_xztLa =0.0_WP 
                ALLOCATE ( ENV22T_xztLa (NCL2,NCL1_io,2) ) ; ENV22T_xztLa =0.0_WP 
                ALLOCATE ( ENV33T_xztLa (NCL2,NCL1_io,2) ) ; ENV33T_xztLa =0.0_WP 
                ALLOCATE ( ENV12T_xztLa (NCL2,NCL1_io,2) ) ; ENV12T_xztLa =0.0_WP 
                ALLOCATE ( ENV13T_xztLa (NCL2,NCL1_io,2) ) ; ENV13T_xztLa =0.0_WP 
                ALLOCATE ( ENV23T_xztLa (NCL2,NCL1_io,2) ) ; ENV23T_xztLa =0.0_WP 
                
                !=============
                ! IN WAVE NUMBER SPACE
                ALLOCATE ( ENV11Z_xztLa (NCL2,NCL3,2) ) ; ENV11Z_xztLa =0.0_WP 
                ALLOCATE ( ENV22Z_xztLa (NCL2,NCL3,2) ) ; ENV22Z_xztLa =0.0_WP 
                ALLOCATE ( ENV33Z_xztLa (NCL2,NCL3,2) ) ; ENV33Z_xztLa =0.0_WP  
                ALLOCATE ( ENV12Z_xztLa (NCL2,NCL3,2) ) ; ENV12Z_xztLa =0.0_WP 
                ALLOCATE ( ENV13Z_xztLa (NCL2,NCL3,2) ) ; ENV13Z_xztLa =0.0_WP 
                ALLOCATE ( ENV23Z_xztLa (NCL2,NCL3,2) ) ; ENV23Z_xztLa =0.0_WP 
                
                !=============velocity and voriticity==============
                ALLOCATE ( VO11X1_xztLa (NCL2,NCL1_io,2) ) ; VO11X1_xztLa = 0.0_WP
                ALLOCATE ( VO12X1_xztLa (NCL2,NCL1_io,2) ) ; VO12X1_xztLa = 0.0_WP
                ALLOCATE ( VO13X1_xztLa (NCL2,NCL1_io,2) ) ; VO13X1_xztLa = 0.0_WP
                
                ALLOCATE ( VO21X1_xztLa (NCL2,NCL1_io,2) ) ; VO21X1_xztLa = 0.0_WP
                ALLOCATE ( VO22X1_xztLa (NCL2,NCL1_io,2) ) ; VO22X1_xztLa = 0.0_WP
                ALLOCATE ( VO23X1_xztLa (NCL2,NCL1_io,2) ) ; VO23X1_xztLa = 0.0_WP
                
                ALLOCATE ( VO31X1_xztLa (NCL2,NCL1_io,2) ) ; VO31X1_xztLa = 0.0_WP
                ALLOCATE ( VO32X1_xztLa (NCL2,NCL1_io,2) ) ; VO32X1_xztLa = 0.0_WP
                ALLOCATE ( VO33X1_xztLa (NCL2,NCL1_io,2) ) ; VO33X1_xztLa = 0.0_WP
                
                ALLOCATE ( VO11X3_xztLa (NCL2,NCL3,2) ) ; VO11X3_xztLa = 0.0_WP
                ALLOCATE ( VO12X3_xztLa (NCL2,NCL3,2) ) ; VO12X3_xztLa = 0.0_WP
                ALLOCATE ( VO13X3_xztLa (NCL2,NCL3,2) ) ; VO13X3_xztLa = 0.0_WP
                
                ALLOCATE ( VO21X3_xztLa (NCL2,NCL3,2) ) ; VO21X3_xztLa = 0.0_WP
                ALLOCATE ( VO22X3_xztLa (NCL2,NCL3,2) ) ; VO22X3_xztLa = 0.0_WP
                ALLOCATE ( VO23X3_xztLa (NCL2,NCL3,2) ) ; VO23X3_xztLa = 0.0_WP
                
                ALLOCATE ( VO31X3_xztLa (NCL2,NCL3,2) ) ; VO31X3_xztLa = 0.0_WP
                ALLOCATE ( VO32X3_xztLa (NCL2,NCL3,2) ) ; VO32X3_xztLa = 0.0_WP
                ALLOCATE ( VO33X3_xztLa (NCL2,NCL3,2) ) ; VO33X3_xztLa = 0.0_WP
                
                ALLOCATE ( EVO11T_xztLa (NCL2,NCL1_io,2) ) ; EVO11T_xztLa = 0.0_WP
                ALLOCATE ( EVO12T_xztLa (NCL2,NCL1_io,2) ) ; EVO12T_xztLa = 0.0_WP
                ALLOCATE ( EVO13T_xztLa (NCL2,NCL1_io,2) ) ; EVO13T_xztLa = 0.0_WP
                
                ALLOCATE ( EVO21T_xztLa (NCL2,NCL1_io,2) ) ; EVO21T_xztLa = 0.0_WP
                ALLOCATE ( EVO22T_xztLa (NCL2,NCL1_io,2) ) ; EVO22T_xztLa = 0.0_WP
                ALLOCATE ( EVO23T_xztLa (NCL2,NCL1_io,2) ) ; EVO23T_xztLa = 0.0_WP
                
                ALLOCATE ( EVO31T_xztLa (NCL2,NCL1_io,2) ) ; EVO31T_xztLa = 0.0_WP
                ALLOCATE ( EVO32T_xztLa (NCL2,NCL1_io,2) ) ; EVO32T_xztLa = 0.0_WP
                ALLOCATE ( EVO33T_xztLa (NCL2,NCL1_io,2) ) ; EVO33T_xztLa = 0.0_WP
                
                ALLOCATE ( EVO11Z_xztLa (NCL2,NCL3,2) ) ; EVO11Z_xztLa = 0.0_WP
                ALLOCATE ( EVO12Z_xztLa (NCL2,NCL3,2) ) ; EVO12Z_xztLa = 0.0_WP
                ALLOCATE ( EVO13Z_xztLa (NCL2,NCL3,2) ) ; EVO13Z_xztLa = 0.0_WP
                
                ALLOCATE ( EVO21Z_xztLa (NCL2,NCL3,2) ) ; EVO21Z_xztLa = 0.0_WP
                ALLOCATE ( EVO22Z_xztLa (NCL2,NCL3,2) ) ; EVO22Z_xztLa = 0.0_WP
                ALLOCATE ( EVO23Z_xztLa (NCL2,NCL3,2) ) ; EVO23Z_xztLa = 0.0_WP
                
                ALLOCATE ( EVO31Z_xztLa (NCL2,NCL3,2) ) ; EVO31Z_xztLa = 0.0_WP
                ALLOCATE ( EVO32Z_xztLa (NCL2,NCL3,2) ) ; EVO32Z_xztLa = 0.0_WP
                ALLOCATE ( EVO33Z_xztLa (NCL2,NCL3,2) ) ; EVO33Z_xztLa = 0.0_WP
        
                IF(thermlflg==1) THEN
                    IF(MYID==0) CALL CHKHDL   ('       Allocating IO x-z averaged thermal variables...',MYID)            
                    ALLOCATE ( T1xzL_io( N2DO(0) ) )                   ; T1xzL_io = 0.0_WP
                    ALLOCATE ( H1xzL_io( N2DO(0) ) )                   ; H1xzL_io = 0.0_WP
                         
                    ALLOCATE ( T2xzL_io( N2DO(0) ) )                   ; T2xzL_io = 0.0_WP
                    ALLOCATE ( D2xzL_io( N2DO(0) ) )                   ; D2xzL_io  = 1.0_WP
                    ALLOCATE ( H2xzL_io( N2DO(0) ) )                   ; H2xzL_io = 0.0_WP
                    
                    ALLOCATE ( DHxzL_io( N2DO(0) ) )                   ; DHxzL_io = 0.0_WP
                    ALLOCATE ( PHxzL_io( N2DO(0) ) )                   ; DHxzL_io = 0.0_WP
                    
                    ALLOCATE ( DVDL1MxzL_io  (N2DO(0), NDV, NDV) )     ; DVDL1MxzL_io  = 0.0_wp
                    ALLOCATE ( DVDL2MxzL_io  (N2DO(0), (NDV-1)*3+NDV, (NDV-1)*3+NDV) )     ; DVDL2MxzL_io  = 0.0_wp
                    ALLOCATE ( DVDL1MHxzL_io (N2DO(0), NDV, NDV) )     ; DVDL1MHxzL_io = 0.0_wp
                    ALLOCATE ( DVDL1MUxzL_io (N2DO(0), NDV, NDV, NDV)) ; DVDL1MUxzL_io = 0.0_wp
                        
                    ALLOCATE ( UHxzL_io  ( N2DO(0),NDV ) )             ; UHxzL_io   = 0.0_WP
                    ALLOCATE ( GHxzL_io  ( N2DO(0),NDV ) )             ; GHxzL_io   = 0.0_WP
                    ALLOCATE ( U2DHxzL_io( N2DO(0),(NDV*(7-NDV))/2+NDV-3 ) ) ; U2DHxzL_io = 0.0_WP
                    
                    ALLOCATE ( DhDL1xzL_io(N2DO(0), NDV) )               ; DhDL1xzL_io = 0.0_wp
                    ALLOCATE ( DhDLPxzL_io(N2DO(0), NDV) )               ; DhDLPxzL_io = 0.0_wp
                    ALLOCATE ( DTDLKxzL_io(N2DO(0), NDV) )               ; DTDLKxzL_io = 0.0_wp
                    ALLOCATE ( DTDLKUxzL_io(N2DO(0), NDV, NDV) )         ; DTDLKUxzL_io = 0.0_wp
                    ALLOCATE ( DTDLKDVDLxzL_io(N2DO(0), NDV, NDV, NDV) ) ; DTDLKDVDLxzL_io = 0.0_wp
                    ALLOCATE ( DhDLMDVDLxzL_io(N2DO(0), NDV, NDV, NDV) ) ; DhDLMDVDLxzL_io = 0.0_wp
                    
                    ALLOCATE ( T1xztL_io( N2DO(0) ) )                   ; T1xztL_io = 0.0_WP
                    ALLOCATE ( H1xztL_io( N2DO(0) ) )                   ; H1xztL_io = 0.0_WP
                         
                    ALLOCATE ( T2xztL_io( N2DO(0) ) )                   ; T2xztL_io = 0.0_WP
                    ALLOCATE ( D2xztL_io ( N2DO(0) ) )                  ; D2xztL_io  = 1.0_WP
                    ALLOCATE ( H2xztL_io( N2DO(0) ) )                   ; H2xztL_io = 0.0_WP
                    
                    ALLOCATE ( DHxztL_io( N2DO(0) ) )                   ; DHxztL_io = 0.0_WP
                    ALLOCATE ( PHxztL_io( N2DO(0) ) )                   ; DHxztL_io = 0.0_WP
                    
                    ALLOCATE ( DVDL1MxztL_io  (N2DO(0), NDV, NDV) )     ; DVDL1MxztL_io  = 0.0_wp
                    ALLOCATE ( DVDL2MxztL_io  (N2DO(0), (NDV-1)*3+NDV, (NDV-1)*3+NDV))    ; DVDL2MxztL_io  = 0.0_wp
                    ALLOCATE ( DVDL1MHxztL_io (N2DO(0), NDV, NDV) )     ; DVDL1MHxztL_io = 0.0_wp
                    ALLOCATE ( DVDL1MUxztL_io (N2DO(0), NDV, NDV, NDV)) ; DVDL1MUxztL_io = 0.0_wp
                        
                    ALLOCATE ( UHxztL_io  ( N2DO(0),NDV ) )             ; UHxztL_io   = 0.0_WP
                    ALLOCATE ( GHxztL_io  ( N2DO(0),NDV ) )             ; GHxztL_io   = 0.0_WP
                    ALLOCATE ( U2DHxztL_io( N2DO(0),(NDV*(7-NDV))/2+NDV-3 ) ) ; U2DHxztL_io = 0.0_WP
                    
                    ALLOCATE ( DhDL1xztL_io(N2DO(0), NDV) )               ; DhDL1xztL_io = 0.0_wp
                    ALLOCATE ( DhDLPxztL_io(N2DO(0), NDV) )               ; DhDLPxztL_io = 0.0_wp
                    ALLOCATE ( DTDLKxztL_io(N2DO(0), NDV) )               ; DTDLKxztL_io = 0.0_wp
                    ALLOCATE ( DTDLKUxztL_io(N2DO(0), NDV, NDV) )         ; DTDLKUxztL_io = 0.0_wp
                    ALLOCATE ( DTDLKDVDLxztL_io(N2DO(0), NDV, NDV, NDV) ) ; DTDLKDVDLxztL_io = 0.0_wp
                    ALLOCATE ( DhDLMDVDLxztL_io(N2DO(0), NDV, NDV, NDV) ) ; DhDLMDVDLxztL_io = 0.0_wp
                    
                    
                    MEMPC_byte = MEMPC_byte + 8*(N2DO(0)*(7+4*NDV+7))
                END IF

            END IF
            
            
      END IF
      
      RETURN
      END SUBROUTINE 
      
      
!**********************************************************************      
    SUBROUTINE MEM_DEALLOCAT
        use mesh_info
        use init_info
        use flow_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE   
        
        DEALLOCATE ( JCL2G )
        DEALLOCATE ( KCL2G )
        DEALLOCATE ( JG2IP )  
        DEALLOCATE ( JG2LC )
      
        DEALLOCATE ( KPV )
        DEALLOCATE ( KMV )
        DEALLOCATE ( KSYM )
        DEALLOCATE ( JGPV )
        DEALLOCATE ( JGMV )
        DEALLOCATE ( JLPV )
        DEALLOCATE ( JLMV )
      

        DEALLOCATE ( CFLVIS )
        
        DEALLOCATE ( AMPH )
        DEALLOCATE ( ACPH )
        DEALLOCATE ( APPH )
      
        DEALLOCATE ( AMVR )
        DEALLOCATE ( ACVR )
        DEALLOCATE ( APVR )
        
        DEALLOCATE ( ZND )
        DEALLOCATE ( ZCC )
        DEALLOCATE ( YND )
        DEALLOCATE ( YCC )
        DEALLOCATE ( RCCI1 )
        DEALLOCATE ( RNDI1 )
        DEALLOCATE ( RCCI2 )
        DEALLOCATE ( RNDI2 )
        DEALLOCATE ( DYFI )
        DEALLOCATE ( DYCI )
        
        DEALLOCATE ( YCL2ND_WFF )
        DEALLOCATE ( YCL2ND_WFB )
        
        
        DEALLOCATE ( Vini     )
        DEALLOCATE ( UU       )
        

        IF(TGFLOWflg) THEN
            DEALLOCATE ( IPV_tg )
            DEALLOCATE ( IMV_tg )
            
            DEALLOCATE ( XND_tg )
            DEALLOCATE ( XCC_tg )
            
            DEALLOCATE ( Q_tg    )
            DEALLOCATE ( PR_tg   )
            DEALLOCATE ( Qtmp_tg )
            DEALLOCATE ( DPH_tg  )
        
            DEALLOCATE ( CONVH0_tg    )
            DEALLOCATE ( RHS_tg       )
            DEALLOCATE ( RHSLLPHI_tg  )
            
             !====for postprocess=================================
            DEALLOCATE ( U1xzL_tg )
            DEALLOCATE ( UPxzL_tg )
            DEALLOCATE ( U2xzL_tg )
            DEALLOCATE ( U3xzL_tg )
    
            DEALLOCATE ( DVDL1xzL_tg )
            DEALLOCATE ( DVDLPxzL_tg )
            DEALLOCATE ( DVDL2xzL_tg )
            
            
            DEALLOCATE ( U1xztL_tg)
            DEALLOCATE ( UPxztL_tg )
            DEALLOCATE ( U2xztL_tg )
            DEALLOCATE ( U3xztL_tg )
    
            DEALLOCATE ( DVDL1xztL_tg )
            DEALLOCATE ( DVDLPxztL_tg )
            DEALLOCATE ( DVDL2xztL_tg )
            
        END IF
        
      
        IF(IOFLOWflg) THEN
        
            DEALLOCATE ( XND_io )
            DEALLOCATE ( XCC_io )

            DEALLOCATE ( Q_io    )
            DEALLOCATE ( G_io    )
            DEALLOCATE ( G0_io   )
            DEALLOCATE ( PR_io   )
            DEALLOCATE ( DPH_io  )
            DEALLOCATE ( Qtmp_io )
            
         
            DEALLOCATE ( RHS_io        )
            DEALLOCATE ( RHSLLPHI_io   )
            DEALLOCATE ( DivU_io       )
            DEALLOCATE ( EXPLT0_io     )
         
            DEALLOCATE ( BC_CONV0     )
            DEALLOCATE ( BC_CONV0_ENG )
            DEALLOCATE ( BC_TDMA      )
            DEALLOCATE ( BC_U_SSTAR   )
           
            !================thermal info============================= 

            DEALLOCATE (RHOH        )
            DEALLOCATE (ENTHALPY    )
            DEALLOCATE (TEMPERATURE )
            DEALLOCATE (THERMCONDT  )
            DEALLOCATE (DENSITY     )
            DEALLOCATE (DENSITY0    )
            DEALLOCATE (VISCOUSITY  )
         
         
            DEALLOCATE (WALLFLUX  )
        
            DEALLOCATE ( RHS_ENERGY  )
            DEALLOCATE ( RHS_ENERGY0 )
            
            
       
            IF(TGFLOWFLG) THEN
                DEALLOCATE ( IPV_io )
                DEALLOCATE ( IMV_io )
            
                !============FOR POSTPROCESS============================
                DEALLOCATE ( U1zL_io )
                DEALLOCATE ( G1zL_io )
                DEALLOCATE ( UPzL_io )
                DEALLOCATE ( U2zL_io )
                DEALLOCATE ( UGzL_io )
                DEALLOCATE ( UGUzL_io )
            
                DEALLOCATE ( DVDL1zL_io )
                DEALLOCATE ( DVDLPzL_io )
                DEALLOCATE ( DVDL2zL_io )
                
                DEALLOCATE ( U1ztL_io )
                DEALLOCATE ( G1ztL_io )
                DEALLOCATE ( UPztL_io )
                
                DEALLOCATE ( U2ztL_io )
                DEALLOCATE ( UGztL_io )
                
                DEALLOCATE ( UGUztL_io )
            
                DEALLOCATE ( DVDL1ztL_io )
                DEALLOCATE ( DVDLPztL_io )
                DEALLOCATE ( DVDL2ztL_io )
            
                            
                DEALLOCATE ( D1zL_io  )
                DEALLOCATE ( D1ztL_io ) 
                DEALLOCATE ( D2zL_io  )
                DEALLOCATE ( D2ztL_io )
                IF(thermlflg==1) THEN           
                    DEALLOCATE ( T1zL_io )
                    DEALLOCATE ( H1zL_io )
                         
                    DEALLOCATE ( T2zL_io )
                    
                    DEALLOCATE ( H2zL_io )
                    
                    DEALLOCATE ( DHzL_io )
                    
                        
                    DEALLOCATE ( UHzL_io )
                    DEALLOCATE ( GHzL_io )
                    
                    DEALLOCATE ( T1ztL_io )
                   
                    DEALLOCATE ( H1ztL_io )
                         
                    DEALLOCATE ( T2ztL_io )
                    
                    DEALLOCATE ( H2ztL_io )
                    
                    DEALLOCATE ( DHztL_io )
                        
                    DEALLOCATE ( UHztL_io )
                    DEALLOCATE ( GHztL_io )
                END IF
                
            ELSE
                DEALLOCATE ( IPV_io )
                DEALLOCATE ( IMV_io )
                
                
                !============FOR POSTPROCESS============================
                DEALLOCATE ( U1xzL_io )
                DEALLOCATE ( G1xzL_io )
                DEALLOCATE ( UPxzL_io )
                DEALLOCATE ( U2xzL_io )
                DEALLOCATE ( UGxzL_io )
                DEALLOCATE ( UGUxzL_io )
                DEALLOCATE ( U3xzL_io )
            
                DEALLOCATE ( DVDL1xzL_io )
                DEALLOCATE ( DVDLPxzL_io )
                DEALLOCATE ( DVDL2xzL_io )
                
                DEALLOCATE ( U1xztL_io )
                DEALLOCATE ( G1xztL_io )
                DEALLOCATE ( UPxztL_io )
                
                DEALLOCATE ( U2xztL_io )
                DEALLOCATE ( UGxztL_io )
                
                DEALLOCATE ( UGUxztL_io )
                DEALLOCATE ( U3xztL_io )
            
                DEALLOCATE ( DVDL1xztL_io )
                DEALLOCATE ( DVDLPxztL_io )
                DEALLOCATE ( DVDL2xztL_io )
            
                            
                DEALLOCATE ( D1xzL_io  )
                DEALLOCATE ( D1xztL_io )
                DEALLOCATE ( D2xztL_io )
                DEALLOCATE ( D2xzL_io  )
                IF(thermlflg==1) THEN
                    DEALLOCATE ( T1xzL_io )
                    DEALLOCATE ( H1xzL_io )
                         
                    DEALLOCATE ( T2xzL_io )
                    
                    DEALLOCATE ( H2xzL_io )
                    
                    DEALLOCATE ( DHxzL_io )
                    
                        
                    DEALLOCATE ( UHxzL_io )
                    DEALLOCATE ( GHxzL_io )
                    
                    
                    DEALLOCATE ( T1xztL_io )
                    DEALLOCATE ( H1xztL_io )
                         
                    DEALLOCATE ( T2xztL_io )
                    
                    DEALLOCATE ( H2xztL_io )
                    
                    DEALLOCATE ( DHxztL_io )
                        
                    DEALLOCATE ( UHxztL_io )
                    DEALLOCATE ( GHxztL_io )
                END IF

            END IF
            
            
      END IF
      
      
      
         
      RETURN
      END SUBROUTINE 
