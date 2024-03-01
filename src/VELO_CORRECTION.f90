!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> CULATION OF SOLENOIDAL VEL FIELD.
!> (Q(N+1)=QHAT-GRAD(DPH)*DT).               
!> Eq.(A1c) in Mehdi paper or Eq.(4.61) in Mehdi thesis.
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE VELOUPDT_tg(NS)
      
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE 

        INTEGER(4),INTENT(IN)  :: NS
      
        INTEGER(4) :: IC, IM
        INTEGER(4) :: JC, JM, JJ
        INTEGER(4) :: KC, KM
        INTEGER(4) :: NII, KS
        REAL(WP)    :: DFX11, DFX22,DFX33
        REAL(WP)    :: COE1
        
        
        !===================update u===================================      
        COE1 = DT*TALP(NS)*DXI
        DO KC=1,NCL3
            DO JC=1,N2DO(MYID)
                DO IC=1,NCL1_tg
                    IM=IMV_tg(IC)
                    DFX11 = DPH_tg(IC,JC,KC)-DPH_tg(IM,JC,KC)
                    Q_tg(IC,JC,KC,1)=Q_tg(IC,JC,KC,1)-DFX11*COE1
                ENDDO
            ENDDO
        ENDDO
      
        !===================update W=================================== 
        COE1 = DT*TALP(NS)*DZI     
        DO KC=1,NCL3
            KM=KMV(KC)
            DO JC=1,N2DO(MYID)
                DO IC=1,NCL1_tg
                    DFX33 = DPH_tg(IC,JC,KC)-DPH_tg(IC,JC,KM)
                    Q_tg(IC,JC,KC,3)=Q_tg(IC,JC,KC,3)-DFX33*COE1
                ENDDO
            ENDDO
        ENDDO
      
      
        !===================update V=================================== 
        NII=1            
        IF(MYID==0) NII = 2    
          
        DO JC=NII,N2DO(MYID) !======MAIN DOMAIN================
            JJ   = JCL2G(JC)
            JM   = JLMV(JC)
            COE1 = DT*TALP(NS)*DYCI(JJ)/RNDI1(jj)
            DO KC=1,NCL3   
                DO IC=1,NCL1_tg
                    DFX22=DPH_tg(IC,JC,KC)-DPH_tg(IC,JM,KC)
                    Q_tg(IC,JC,KC,2)=Q_tg(IC,JC,KC,2)-DFX22*COE1
                ENDDO
            ENDDO
        ENDDO
        
!        IF (MYID.EQ.0) THEN          !======BOTTOM==========repeat in the interface======
!            IF (ICASE==IPIPEC) THEN
!                DO KC=1,NCL3
!                    KS = KSYM(KC)
!                    DO IC=1,NCL1_tg
!                        Q_tg(IC,1,KC,2)=0.0_WP
!                        Q_tg(IC,0,KC,2)=Q_tg(IC,2,KS,2)
!                    ENDDO
!                ENDDO
!            ELSE
!                DO KC=1,NCL3
!                    DO IC=1,NCL1_tg
!                        Q_tg(IC,1,KC,2)=0.0_WP
!                        Q_tg(IC,0,KC,2)=0.0_WP
!                    ENDDO
!                ENDDO
!            END IF
!        ENDIF
      
!        IF (MYID.EQ.NPSLV) THEN      !======TOP================
!            DO KC=1,NCL3
!                DO IC=1,NCL1_tg
!                    Q_tg(IC,N2DO(MYID)+1,KC,2)=0.0_WP
!                ENDDO
!            ENDDO
!        ENDIF

        RETURN
    END SUBROUTINE
    
    
!**************************************************************************
    SUBROUTINE MASSFLUX_UPDATE_IO(NS)
        use thermal_info
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE 

        INTEGER(4),INTENT(IN) :: NS  
      
        INTEGER(4) :: IC, IM
        INTEGER(4) :: JC, JM, JJ
        INTEGER(4) :: KC, KM, KS
        INTEGER(4) :: NXI, NYI, IDR, NXE
        REAL(WP)    :: DFX11, DFX22,DFX33
        REAL(WP)    :: COE1
      
        
      
        !=====================================X=============================
        IDR = 1 
        NXI = 1
        NXE = NCL1_IO
        IF(TGFLOWFLG) THEN 
            NXI = 2
            NXE = NCL1_io+1
        END IF
        
        NYI = 1
        COE1 = DT*TALP(NS)*DXI
        DO IC=NXI,NXE
            IM=IMV_io(IC)
            DO JC=NYI,N2DO(MYID)
                DO KC=1,NCL3
                    DFX11= DPH_io(IC,JC,KC)-DPH_io(IM,JC,KC)    !d\phi / dx
                    G_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)-DFX11*COE1
                ENDDO
            ENDDO
        ENDDO
      
        !=====================================Z=============================
        IDR = 3  
        NXI = 1
        NYI = 1 
        COE1 = DT*TALP(NS)*DZI
        DO KC=1,NCL3
            KM=KMV(KC)
            DO JC=NYI,N2DO(MYID)
                DO IC=NXI,NCL1_io
                    DFX33=DPH_io(IC,JC,KC)-DPH_io(IC,JC,KM)
                    G_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)-DFX33*COE1
                ENDDO
            ENDDO
        ENDDO
      
        !=====================================Y=============================
        IDR = 2
        NXI = 1
        NYI = 1
        IF(MYID==0 .AND. ICASE.NE.IBOX3P) NYI = 2
        DO JC=NYI,N2DO(MYID)
            JJ=JCL2G(JC)
            JM=JLMV(JC)
            COE1 = DT*TALP(NS)*DYCI(JJ)/RNDI1(jj)
            DO IC=NXI,NCL1_io
                DO KC=1,NCL3
                    DFX22= DPH_io(IC,JC,KC)-DPH_io(IC,JM,KC)
                    G_io(IC,JC,KC,IDR)=G_io(IC,JC,KC,IDR)-DFX22*COE1
                ENDDO
            ENDDO
        ENDDO
        
        
        CALL INTFC_VARS3(NCL1S,NCL1E,NCL1S,NCL1E,G_io)
        CALL BC_WALL_G_io

        RETURN
    END SUBROUTINE
    
    
