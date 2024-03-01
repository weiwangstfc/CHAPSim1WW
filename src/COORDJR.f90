    SUBROUTINE COORDJR
    !> @note : plane channel  y=(-1, 1)
    !>         circular tube  y=(0, 1), where the polar centre is 0
        use mesh_info
        use init_info
        IMPLICIT NONE     

        INTEGER(4)  :: I
        INTEGER(4)  :: J, N
        INTEGER(4)  :: K
        INTEGER(4)  :: NFIL
        REAL(WP)      :: X2 
        REAL(WP)      :: COE1
        real(wp) :: dxplus, yMINplus, dzplus, L2, HY, dyMIN, dyMAX, yMAXplus
        integer(4) :: MGRIDHF
        REAL(8)    :: dup
        
        IF(MYID.NE.0) RETURN

        RNDI1(:)=1.0_WP
        RCCI1(:)=1.0_WP
        RNDI2(:)=1.0_WP
        RCCI2(:)=1.0_WP
        
        SELECT CASE (icase)
           CASE (ICHANL)
           
                COE1    = 0.5_wp   
                DO J=1,NND2             !!> @note: do j=1,NND2
                    X2=DBLE(J-1)/DBLE(NND2-1)   !!> @note: X2=(J-1)/(NND2-1), epslon, the relative computational domain 
                    YND(J)=DTANH(STR2*(X2-COE1))/DTANH(STR2*COE1)
                ENDDO    
                
                DO J=1,NCL2
                    YCC(J)=(YND(J)+YND(J+1))*0.50_WP
                ENDDO
                
           CASE (IPIPEC)
           
                DO J=1,NND2
                    X2=DBLE(J-1)/DBLE(NND2-1)  
                    YND(j)=DTANH(STR2*(X2))/DTANH(STR2) 
                    if(J.eq.1) then
                        RNDI1(j)= 1.0E20_wp ! A large no.
                        RNDI2(j)= 1.0E20_wp ! A large no.
                    else
                        RNDI1(j)=1.0_WP / YND(j)
                        RNDI2(j)=RNDI1(j)*RNDI1(j)
                    end if
                ENDDO 
                DO J=1,NCL2
                    YCC(J)=(YND(j)+YND(j+1))*0.50_WP
                    RCCI1(J)=1.0_WP/ YCC(J)
                    RCCI2(J)=RCCI1(J)*RCCI1(J)
                END DO
                RCCI1(0) = RCCI1(1)
                RCCI2(0) = RCCI2(1)
                RCCI1(NND2) = RNDI1(NND2) ! last variables on the wall ! 1.0_WP/( YCC(NCL2) + (YND(NND2)-YND(NCL2)) )
                RCCI2(NND2) = RNDI2(NND2) ! last variables on the wall ! RCCI1(NND2) * RCCI1(NND2)
                
           CASE (IANNUL)
           
                COE1    = 0.5_wp   
                DO J=1,NND2             !!> @note: do j=1,NND2
                    X2=DBLE(J-1)/DBLE(NND2-1)   !!> @note: X2=(J-1)/(NND2-1), epslon, the relative computational domain 
                    RNDI1(J)=1.0_WP/ ( (DTANH(STR2*(X2-COE1))/DTANH(STR2*COE1) + 1.0_WP)*0.5_WP*(HYT-HYB)+HYB)
                    RNDI2(J)=RNDI1(J)*RNDI1(J)
                    YND(J)=(DTANH(STR2*(X2-COE1))/DTANH(STR2*COE1) + 1.0_WP)*0.5_WP*(HYT-HYB)+HYB
                ENDDO        
                
                DO J=1,NCL2
                    RCCI1(J)=1.0_WP/( (1.0_WP/RNDI1(J)+1.0_WP/RNDI1(J+1))*0.50_WP )
                    RCCI2(J)=RCCI1(J)*RCCI1(J)
                    YCC(J)=(YND(J)+YND(J+1))*0.50_WP
                ENDDO
                
           CASE (IBOX3P)
           
                DO  J=1,NND2
                    YND(J) = HYB+DBLE(J-1)*((HYT-HYB)/DBLE(NCL2))
                ENDDO
                DO J=1,NCL2
                    YCC(J) = 0.5_WP*(YND(J)+YND(J+1))
                END DO
                
           CASE DEFAULT
           
                COE1    = 0.5_wp   
                DO J=1,NND2             !!> @note: do j=1,NND2
                    X2=DBLE(J-1)/DBLE(NND2-1)   !!> @note: X2=(J-1)/(NND2-1), epslon, the relative computational domain 
                    YND(J)=DTANH(STR2*(X2-COE1))/DTANH(STR2*COE1)
                ENDDO    
                
                DO J=1,NCL2
                    YCC(J)=(YND(J)+YND(J+1))*0.50_WP
                ENDDO
                
        END SELECT
        
        
        !===========================y distance==========================================
        DO J=1,NCL2
            DYFI(J)=1.0_WP/ ( YND(J+1)-YND(J) ) 
        ENDDO   
        
        DO J=2,NND2-1
            DYCI(J)= 1.0_WP/ ( ( YND(J+1)-YND(J-1) )*0.50_WP )  
        END DO
        
        SELECT CASE (icase)
            CASE (ICHANL)
                DYCI(1)   = 2.0_WP/( YND(2) -YND(1)    )   ! half length
                DYCI(NND2)= 2.0_WP/( YND(NND2)-YND(NCL2) ) ! half length
            CASE (IPIPEC)
                DYCI(1)   = 1.0_WP/( YND(2) -YND(1)    )   ! pipe centre, whole length
                DYCI(NND2)= 2.0_WP/( YND(NND2)-YND(NCL2) ) ! half length
            CASE (IANNUL)
                DYCI(1)   = 2.0_WP/( YND(2) -YND(1)    )   ! half length
                DYCI(NND2)= 2.0_WP/( YND(NND2)-YND(NCL2) ) ! half length
            CASE (IBOX3P)
                DYCI(1)   = 1.0_WP/(DABS( YCC(1) -YND(1))+DABS(YND(NND2)-YCC(NCL2)))   ! whole length
                DYCI(NND2)= 1.0_WP/(DABS( YCC(1) -YND(1))+DABS(YND(NND2)-YCC(NCL2)))   ! whole length
            CASE DEFAULT
                DYCI(1)   = 2.0_WP/( YND(2) -YND(1)    )   ! half length
                DYCI(NND2)= 2.0_WP/( YND(NND2)-YND(NCL2) ) ! half length
        END SELECT
        
        !================================x=============================================
        IF(TGFLOWflg .AND. IOFLOWflg) THEN
            DO  I=1,NND1_io
                XND_io(I) = DBLE(I-1)*DX
            ENDDO
            
            DO  I=1,NND1_tg
                XND_tg(I) = (-1.0_WP)*DBLE(NCL1_tg)*DX+DBLE(I-1)*DX
            ENDDO
            
            DO I=1, NCL1_tg
                XCC_TG(I) = 0.5_WP*(XND_TG(I)+XND_TG(I+1))
            END DO
            
            DO I=1, NCL1_IO
                XCC_IO(I) = 0.5_WP*(XND_IO(I)+XND_IO(I+1))
            END DO
        
        ELSE
            IF(TGFLOWflg) THEN
                DO  I=1,NND1_tg
                    XND_tg(I) = DBLE(I-1)*DX
                ENDDO
                
                DO I=1, NCL1_tg
                    XCC_TG(I) = 0.5_WP*(XND_TG(I)+XND_TG(I+1))
                END DO
            END IF
            
            IF(IOFLOWflg) THEN
                DO  I=1,NND1_io
                    XND_io(I) = DBLE(I-1)*DX
                ENDDO
                DO I=1, NCL1_IO
                    XCC_IO(I) = 0.5_WP*(XND_IO(I)+XND_IO(I+1))
                END DO
                
            END IF
        
        END IF
        
        !================================z=============================================
        DO K=1,NND3
            ZND(K)=DBLE(K-1)*DZ
        ENDDO
        DO K=1,NCL3
            ZCC(K) = 0.5_WP*(ZND(K)+ZND(K+1))
        END DO
        
        ! debug4chapsim2
        ! IF(ICASE==IBOX3P) THEN
        !     XND_io(:) = XND_io(:) - PI
        !     XCC_io(:) = XCC_io(:) - PI
        !     ZND(:) = ZND(:) - PI
        !     ZCC(:) = ZCC(:) - PI
        ! END IF

        !==============================shared variables================================
        CALL AREA_INOUTLET
        CALL MESH_YWEIGHT_FACTORS
        
        !===========identify the y locations for x-z statistics=====
        MGRIDHF = MGRID/2+1
        
        JGMOV(1)       = JINI
        JGMOV(MGRID)   = NCL2-JINI+1
        JGMOV(MGRIDHF) = NCL2/2
        dup = DBLE(JGMOV(MGRIDHF) - JGMOV(1) )/DBLE(MGRIDHF * (MGRIDHF-1)/2)

        DO J=2,MGRIDHF-1
            JGMOV(J) = JGMOV(J-1) + INT(DBLE(J-1)*dup)
        END DO

        DO J=MGRID-1, MGRIDHF+1, -1
            JGMOV(J) = JGMOV(MGRID)-( JGMOV(MGRID-J+1)-JGMOV(1) )
        END DO
        
    
        !===================writing mesh info into files===============================
        IF(myid==0) THEN
        
            NFIL=14  
            OPEN(NFIL,FILE='CHK_PROBE_for_spectra.dat')
            REWIND NFIL
            write(NFIL,'(A)') '#MGRID, JG, dJ, YCC(J), dYCC(J)'
            WRITE(NFIL,'(3I4.1,2ES13.5)') 1, JGMOV(1), 0, YCC(1), 0.0_wp
            DO N=2, MGRID
                J = JGMOV(N)
                WRITE(NFIL,'(3I4.1,2ES13.5)') N, JGMOV(N), JGMOV(N)-JGMOV(N-1), YCC(J), YCC(J)-YCC(J-1)  
            END DO
            CLOSE(NFIL)
            
            !===============X NODE COORNIDATES=========================================
            NFIL=NFIL+1  
            OPEN(NFIL,FILE='CHK_COORD_XND.dat')
            REWIND NFIL
            WRITE(NFIL,'(A)') '# INDEX, XND'  
            IF(TGFLOWflg .AND. IOFLOWflg) THEN
                WRITE(NFIL,*) '#', NND1_io + NND1_tg
                DO I=1,NND1_tg
                    WRITE(NFIL,*) I, XND_tg(I)
                END DO
                DO I=1,NND1_io
                    WRITE(NFIL,*) I, XND_io(I)
                END DO
            ELSE
                IF(TGFLOWflg) THEN
                    WRITE(NFIL,*) '#', NND1_tg
                    DO I=1,NND1_tg
                        WRITE(NFIL,*) I, XND_tg(I)
                    END DO
                END IF
                IF(IOFLOWflg) THEN
                    WRITE(NFIL,*) '#', NND1_IO
                    DO I=1,NND1_io
                        WRITE(NFIL,*) I, XND_io(I)
                    END DO
                END IF
            END IF 
            CLOSE(NFIL)
    
            NFIL=NFIL + 1
            OPEN(NFIL,FILE='CHK_COORD_YND.dat')
            REWIND NFIL
            IF(icase .EQ. ICHANL .or. icase .EQ. IBOX3P) THEN
                WRITE(NFIL,'(A)') '# JND, YND, DYF, YND+'   
                DO J=1,NND2
                    WRITE(NFIL,'(I10, 3ES15.7)') J, YND(J), 1.0_WP/DYCI(J), (ALX2-dabs(YND(J)))* REN * DSQRT(0.5_wp*CFGV)
                END DO  
                WRITE(NFIL,'(A)') '# JCL, YCC, DYC, DYC+'     
                DO J=1,NCL2
                    WRITE(NFIL,'(I10, 3ES15.7)') J, YCC(J), 1.0_WP/DYFI(J), 1.0_WP/DYFI(J)* REN * DSQRT(0.5_wp*CFGV)
                ENDDO
            ELSE IF (icase.EQ.IPIPEC .or. icase.EQ.IANNUL) THEN
                WRITE(NFIL,'(A)') '# JND, RND, RND2, RNDI1, RNDI2'  
                DO J=1,NND2
                    WRITE(NFIL,'(I10, 4ES15.7)') J, 1.0_WP/RNDI1(J), 1.0_WP/RNDI2(J), RNDI1(J), RNDI2(J)
                END DO
                WRITE(NFIL,'(A)') '# JCC, RCC, DYC, DYF, RCCI1, RCCI2'
                DO J=1,NCL2
                    WRITE(NFIL,'(I10, 5ES13.5)') J,1.0_WP/RCCI1(J),1.0_WP/DYCI(J),1.0_WP/DYFI(J), RCCI1(J), RCCI2(J)
                ENDDO
            ELSE
            END IF
            CLOSE(NFIL)
          
            NFIL=NFIL + 2
            OPEN(NFIL,FILE='CHK_COORD_ZND.dat')
            REWIND NFIL
            WRITE(NFIL,*) '#', NND3       
            DO K=1,NND3
                WRITE(NFIL,*) K, ZND(K)
            END DO  
            CLOSE(NFIL)
        END IF
        
        
        !========================MESH INFO==================================================================
        IF(ICASE==IPIPEC) THEN
            dyMIN = 1.0_wp/DYFI(NCL2)
            dyMAX = 1.0_WP/DYFI(1)
            HY    = HYT-HYB
        ELSE
            dyMIN = 1.0_wp/DYFI(1)
            dyMAX = 1.0_WP/DYFI(NCL2/2)
            HY    = HYT-HYB
        END IF
        dxplus = DX * REN * DSQRT(0.5_wp*CFGV)
        dzplus = DZ * REN * DSQRT(0.5_wp*CFGV)
        yMINplus = DYMIN* REN * DSQRT(0.5_wp*CFGV)
        yMAXplus = DYMAX* REN * DSQRT(0.5_wp*CFGV)
    
  
        WRITE(*,'(A)') '#**********Mesh and flow details based on initial flow field******************************'
        IF(TGFLOWflg) &
        WRITE(*,'(A,I5,3(2X,A,F9.5))') '#   NCL1_TG= ',NCL1_tg, 'HX_TG=', ALX1(1),'DX=  ', DX,   'DX+= ', dxplus
        IF(IOFLOWflg) &
        WRITE(*,'(A,I5,3(2X,A,F9.5))') '#   NCL1_IO= ',NCL1_io, 'HX_IO=', ALX1(2),'DX=  ', DX,   'DX+= ', dxplus
        WRITE(*,'(A,I5,3(2X,A,F9.5))') '#   NCL3=    ',NCL3,    'HZ=   ', ALX3,   'DZ=  ', DZ,   'DZ+= ', dzplus
        WRITE(*,'(A,I5,5(2X,A,F9.5))') '#   NCL2=    ',NCL2,    'HY=   ', HY,     'DY1= ', DYMIN,'Y1+= ', yMINplus
        WRITE(*,'(A,I5,5(2X,A,F9.5))') '#   NCL2=    ',NCL2,    'HY=   ', HY,     'DYc= ', DYMAX,'Yc+= ', yMAXplus  
                                                                                    
        IF(TGFLOWflg) &
        WRITE(*,'(A,I12)'            ) '#   MESH_SIZE_TG   =',NCL1_tg*NCL2*NCL3
        IF(IOFLOWflg) &
        WRITE(*,'(A,I12)'            ) '#   MESH_SIZE_IO   =',NCL1_io*NCL2*NCL3
        IF(IOFLOWflg .AND. TGFLOWFLG) &
        WRITE(*,'(A,I12)'            ) '#   MESH_SIZE_total=',(NCL1_io+NCL1_tg)*NCL2*NCL3
        WRITE(*,'(A)') '#****************************************************************************************'
        WRITE(*,'(A,2X,F18.5)')        '#   CONSTANT MEMORY PER CORE (Mb) = ', DBLE(MEMPC_byte)/1024.0_WP/1024.0_WP
        WRITE(*,'(A)') '#****************************************************************************************'
        
        !STOP !
      
        RETURN
    END SUBROUTINE
      
      
!****************************************************************************
    SUBROUTINE AREA_INOUTLET
        USE MESH_INFO
        USE init_info
        IMPLICIT NONE
     
        INTEGER(4) :: JC,KC,IC
        real(wp)   :: AREA_INLET_ideal, VOLM_tg_ideal, VOLM_IO_ideal
      
        IF(MYID.NE.0) RETURN
        
        
        SELECT CASE (icase)
            CASE (ICHANL)
                AREA_INLET_ideal = HZ * (HYT-HYB)
            CASE (IPIPEC)
                AREA_INLET_ideal = PI * (HYT-HYB) * (HYT-HYB)
            CASE (IANNUL)
                AREA_INLET_ideal = PI * HYT * HYT - PI * HYB * HYB
            CASE (IBOX3P)
                AREA_INLET_ideal = HZ * (HYT-HYB)
            CASE DEFAULT
                AREA_INLET_ideal = HZ * (HYT-HYB)
        END SELECT
        !===================common========================
        AREA_INLET=0.0_WP            
        DO JC=1,NCL2
            DO KC=1,NCL3
                AREA_INLET   = AREA_INLET + 1.0_WP/RCCI1(JC)/DYFI(JC)/DZI
            ENDDO
        ENDDO
        
        
        IF(TGFLOWflg) THEN
            VOLM_tg_ideal = AREA_INLET_ideal * HX_tg
            VOLM_tg=0.0_wp  
            DO KC=1,NCL3  
                DO JC=1,NCL2
                    DO IC=1,NCL1_tg
                        VOLM_tg =VOLM_tg + 1.0_WP/RCCI1(JC)/DYFI(JC) !/DZI/DXI
                    END DO
                END DO
            END DO
        END IF
        
        
        IF(IOFLOWflg) THEN
            VOLM_IO_ideal = AREA_INLET_ideal * HX_io
            VOLM_IO=0.0_wp  
            DO KC=1,NCL3  
                DO JC=1,NCL2
                    DO IC=1,NCL1_IO
                        VOLM_IO =VOLM_IO + 1.0_WP/RCCI1(JC)/DYFI(JC) !/DZI/DXI
                    END DO
                END DO
            END DO
        END IF
        

        IF (MYID.EQ.0) THEN
            CALL CHKRLHDL  ('        INLET/OUTLET SURFACE AREA(ideal)=    ',MYID,AREA_INLET_ideal)
            CALL CHKRLHDL  ('        INLET/OUTLET SURFACE AREA(calcu)=    ',MYID,AREA_INLET)
            IF(TGFLOWflg) CALL CHKRLHDL  ('        Volume of TG region(ideal)=          ',MYID,VOLM_tg_ideal)
            IF(TGFLOWflg) CALL CHKRLHDL  ('        Volume of TG region(calcu)=          ',MYID,VOLM_tg/DZI/DXI)
            IF(IOFLOWflg) CALL CHKRLHDL  ('        Volume of IO region(ideal)=          ',MYID,VOLM_IO_ideal)
            IF(IOFLOWflg) CALL CHKRLHDL  ('        Volume of IO region(calcu)=          ',MYID,VOLM_IO/DZI/DXI)
        END IF
        
        
        RETURN
    END SUBROUTINE      
      
!*****************************************************************************      
    SUBROUTINE MESH_YWEIGHT_FACTORS
        USE MESH_INFO
        USE thermal_info
        use init_info
        IMPLICIT NONE
        INTEGER(4) :: JC, JM, NFIL
     
        IF(MYID.NE.0) RETURN
        
        DO JC=2,NCL2
            JM = JGMV(JC)
            YCL2ND_WFF(JC) = ( YND(JC)-YCC(JM) ) / ( YCC(JC)-YCC(JM) ) 
            YCL2ND_WFB(JC) = ( YCC(JC)-YND(JC) ) / ( YCC(JC)-YCC(JM) ) 
        END DO


        SELECT CASE (icase)
            CASE (ICHANL)
                YCL2ND_WFB(1)    = 1.00_WP
                YCL2ND_WFF(1)    = 0.00_WP
                YCL2ND_WFB(NND2) = 0.00_WP
                YCL2ND_WFF(NND2) = 1.00_WP
            CASE (IPIPEC)
                YCL2ND_WFF(1)    = 0.50_WP
                YCL2ND_WFB(1)    = 0.50_WP
                YCL2ND_WFB(NND2) = 0.00_WP
                YCL2ND_WFF(NND2) = 1.00_WP
            CASE (IANNUL)
                YCL2ND_WFB(1)    = 1.00_WP
                YCL2ND_WFF(1)    = 0.00_WP
                YCL2ND_WFB(NND2) = 0.00_WP
                YCL2ND_WFF(NND2) = 1.00_WP
            CASE (IBOX3P)
                YCL2ND_WFB(1)    = 0.50_WP
                YCL2ND_WFF(1)    = 0.50_WP
                YCL2ND_WFB(NND2) = 0.50_WP
                YCL2ND_WFF(NND2) = 0.50_WP
            CASE DEFAULT
                YCL2ND_WFB(1)    = 1.00_WP
                YCL2ND_WFF(1)    = 0.00_WP
                YCL2ND_WFB(NND2) = 0.00_WP
                YCL2ND_WFF(NND2) = 1.00_WP
        END SELECT

        
        XND2CL = 0.50_WP
        YND2CL = 0.50_WP
        ZND2CL = 0.50_WP
        
        IF(myid==0) THEN
            NFIL=15  
          
            OPEN(NFIL,FILE='CHK_COEF_YCL2ND.dat')
            REWIND NFIL
            WRITE(NFIL,*) 'JND, DYCI, YCL2ND_WFF, YCL2ND_WFB'
            DO JC=1,NND2
                WRITE(NFIL,'(I4.1,3ES13.5)') JC, DYCI(JC), YCL2ND_WFF(JC),YCL2ND_WFB(JC)
            END DO
            CLOSE(NFIL)
        END IF
        
     
        RETURN
    END SUBROUTINE  
    
    
    
    SUBROUTINE BCAST_COORDJR
        use mesh_info
        use init_info
        IMPLICIT NONE
        
        
        CALL MPI_BCAST( YND(1),   NND2,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( YCC(1),   NCL2,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( RNDI1(1), NND2,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( RNDI2(1), NND2,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( RCCI1(0), NND2+1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( RCCI2(0), NND2+1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( DYFI(1),  NCL2,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DYCI(1),  NND2,   MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        IF(TGFLOWflg .AND. IOFLOWflg) THEN
            CALL MPI_BCAST( XND_io(1), NND1_io, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( XCC_io(1), NCL1_io, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( XND_tg(1), NND1_tg, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( XCC_tg(1), NCL1_tg, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        ELSE
            IF(TGFLOWflg) THEN
                CALL MPI_BCAST( XND_tg(1), NND1_tg, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
                CALL MPI_BCAST( XCC_tg(1), NCL1_tg, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            END IF
            
            IF(IOFLOWflg) THEN
                CALL MPI_BCAST( XND_io(1), NND1_io, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
                CALL MPI_BCAST( XCC_io(1), NCL1_io, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            END IF
        END IF
        
        CALL MPI_BCAST( ZND, NND3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ZCC, NCL3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        
        CALL MPI_BCAST( AREA_INLET, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( VOLM_tg,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( VOLM_io,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( YCL2ND_WFF, NND2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( YCL2ND_WFB, NND2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( XND2CL,     1,    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( YND2CL,     1,    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ZND2CL,     1,    MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( JGMOV, MGRID, MPI_INTEGER4, 0, ICOMM, IERROR )
        
        
        
        RETURN
    END SUBROUTINE
    
      
      
