!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Calculate the geometry, cell no, uniform mesh size info.
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
!**********************************************************************
    SUBROUTINE CONSPARA
        use init_info
        use mesh_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
        
        REAL(WP) :: DX1, DX2
        INTEGER(4) :: I, J
       
        IF(MYID.NE.0) RETURN
        
        !=================
        DO I=1,NDV
            DO J=1,NDV
                IF(I==J) THEN
                    Kronecker_delta(I,J) = 1
                ELSE
                    Kronecker_delta(I,J) = 0
                END IF
            END DO
        END DO
        !=================
       
        
       
        ALX1(1)=HX_tg
        ALX1(2)=HX_io
        
        IF(TGFLOWflg .AND. IOFLOWflg) THEN
            DX1=  ALX1(1)/DBLE(NCL1_TG)
            DX2=  ALX1(2)/DBLE(NCL1_IO)
            IF(DABS(DX1/DX2-1.0_WP) .GT. 1.0E-8_WP) THEN
                CALL ERRHDL('# DX_TG/=DX_IO',MYID)
            ELSE
                DX = DX1
            END IF
        ELSE
            IF(TGFLOWflg) THEN
                DX=  ALX1(1)/DBLE(NCL1_tg)  !check todo for clustering grids in x direction.
            END IF
            IF(IOFLOWflg) THEN
                DX=  ALX1(2)/DBLE(NCL1_io)
            END IF
        END IF

        DXI  = 1.0_WP/DX
        DXQI = DXI*DXI
        

        if(icase.eq.ICHANL .OR. ICASE==IBOX3P) then
            ALX3=HZ
        ELSE IF (icase.EQ.IPIPEC .or. icase.EQ.IANNUL) THEN
            ALX3=2.0_WP*PI
        ELSE
        endif
        DZ   = ALX3/DBLE(NCL3)
        DZI  = 1.0_WP/DZ
        DZQI = DZI*DZI
        
        ALX2   =(HYT-HYB)/2.0_WP   ! HALF CHANNEL HEIGHT
        
        CVISC = 1.0_WP/REN
       
        IF(TGFLOWflg) THEN
            VL1313_tg = 1.0_WP/DBLE(NCL1_tg*NCL3)
        END IF
        
        IF(IOFLOWflg) THEN
            VL1313_io = 1.0_WP/DBLE(NCL1_io*NCL3)
        END IF    
        
        !=========define threshold level H for quadrant analysis===============
        QUADHV(1)=0.00_WP
        QUADHV(2)=0.25_WP
        QUADHV(3)=0.50_WP
        QUADHV(4)=0.75_WP
        QUADHV(5)=1.00_WP
        QUADHV(6)=2.00_WP
        QUADHV(7)=3.00_WP
        QUADHV(8)=4.00_WP
        QUADHV(9)=5.00_WP    
        !
        
       
        RETURN
    END SUBROUTINE
    
    SUBROUTINE BCAST_CONSPARA
        use init_info
        use mesh_info
        use thermal_info
        use postprocess_info
        IMPLICIT NONE
        
        CALL MPI_BCAST( PI,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( ALX1, 2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ALX2, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ALX3, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( DX,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DXI,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DXQI, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( DZ,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DZI,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DZQI, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( CVISC, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( VL1313_tg,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( VL1313_io,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( Kronecker_delta,9, MPI_INTEGER4, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( QUADHV, QUADHN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        RETURN
    END SUBROUTINE
