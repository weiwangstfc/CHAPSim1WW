!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> decomposite the mesh in y-direction
!> @param SW11  the starting global cell-index in one rank   \n
!> @param EW11  the ending global cell-index in one rank     \n
!> @param IDSWT  IDSWT(0:SIZE-1) the starting y global cell index in this rank.
!>               IDSWT(i) means the starting global cell index in MYID=i.
!>               This variable exists in each ranks.
!> @param IDEWT  IDEWT(0:SIZE-1) the ending y global cell index in this rank 
!> @param NPSLV
!>     @param  DIMS(1:1)  number of processors, =SIZE, 
!                       only used in MPI_CART_CREATE 
!>     @param  NDIM11   dimention no. of Cartesian mesh, 
!>                      only used in MPI_CART_CREATE
!>     @param ICOORDS   only used in MPI_CART_GET
!>     @param PERIODS   only used in MPI_CART_CREAT/GET
!>     @param REORDER   only used in MPI_CART_CREAT
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE mesh_Ydecomp
        use mesh_info
        IMPLICIT NONE      

        INTEGER(4)  :: NTOT1
        INTEGER(4)  :: NLOCAL
        INTEGER(4)  :: DEFICIT
        INTEGER(4)  :: IP
         
!>     @note GET THE SIZE OF THE PROBLEM=NTOT1
!>           WE ARE CHOOSING 1-D DECOMPOSITION IN Y DIRECTION.
!>           Decomposition of mesh is only in Y (or r) direction

        ALLOCATE ( JDSWT(0:NPSLV) )
        ALLOCATE ( JDEWT(0:NPSLV) )
        ALLOCATE ( N2DO (0:NPSLV) )
        JDSWT=0
        JDEWT=0
!>       @warning: y velocity is stored on y cell surface, which are N2M+1
        NTOT1=NCL2  

        DO IP=0,NPSLV
          
            NLOCAL  = NTOT1 / NPTOT  
            DEFICIT = MOD(NTOT1,NPTOT)
          
            JDSWT(IP) = IP * NLOCAL + 1
            JDSWT(IP) = JDSWT(IP) + MIN(IP,DEFICIT)
!
            IF (IP .LT. DEFICIT) THEN
                NLOCAL = NLOCAL + 1
            ENDIF
!
            JDEWT(IP) = JDSWT(IP) + NLOCAL - 1
!
            IF ( (JDEWT(IP) .GT. NTOT1) .OR. ( IP .EQ. NPSLV ) ) THEN
                JDEWT(IP) = NTOT1 
            END IF
        END DO

        !call CHKHDL('   Ydecomp,      rank,  Start,  End,  Interval',MYID)
        DO IP=0,NPSLV
            N2DO(IP) = JDEWT(IP)-JDSWT(IP)+1
            !WRITE(*,'(A, 4I8.1)') '# MYID   0        Ydecomp,',ip, JDSWT(IP), JDEWT(IP),N2DO(IP)
        END DO
       
        RETURN
    end subroutine 
       
!**********************************************************************       
    SUBROUTINE mesh_Zdecomp
        use mesh_info
        IMPLICIT NONE      

        INTEGER(4)  :: NTOT1
        INTEGER(4)  :: NLOCAL
        INTEGER(4)  :: DEFICIT
        INTEGER(4)  :: IP
             
!>     @note GET THE NPTOT OF THE PROBLEM=NTOT1
!>           WE ARE CHOOSING 1-D DECOMPOSITION IN Y DIRECTION.
!>           Decomposition of mesh is only in Y (or r) direction
       
        ALLOCATE ( KDSWT(0:NPSLV) )
        ALLOCATE ( KDEWT(0:NPSLV) )
        ALLOCATE ( N3DO (0:NPSLV) )
        KDSWT=0
        KDEWT=0

!>       @warning: y velocity is stored on y cell surface, which are N2M+1
        NTOT1=NCL3   

        DO IP=0,NPSLV
          
            NLOCAL  = NTOT1 / NPTOT  
            DEFICIT = MOD(NTOT1,NPTOT)
          
            KDSWT(IP) = IP * NLOCAL + 1
            KDSWT(IP) = KDSWT(IP) + MIN(IP,DEFICIT)
!
            IF (IP .LT. DEFICIT) THEN
                NLOCAL = NLOCAL + 1
            ENDIF
!
            KDEWT(IP) = KDSWT(IP) + NLOCAL - 1
!
            IF ( (KDEWT(IP) .GT. NTOT1) .OR. ( IP .EQ. NPSLV ) ) THEN
                KDEWT(IP) = NTOT1 
            END IF
        END DO

        !call CHKHDL('   Zdecomp,      rank,  Start,  End,  Interval',MYID)
        DO IP=0,NPSLV
            N3DO(IP) = KDEWT(IP)-KDSWT(IP)+1
            !WRITE(*,'(A, 4I8.1)') '# MYID   0        Zdecomp,',ip, KDSWT(IP), KDEWT(IP),N3DO(IP)
        END DO
          
        RETURN
    end subroutine 
       
!***********************************************************************       
    SUBROUTINE mesh_decomp_bcast
        use mesh_info
        IMPLICIT NONE
       
        IF(MYID.NE.0) THEN
            ALLOCATE ( JDSWT(0:NPSLV) )
            ALLOCATE ( JDEWT(0:NPSLV) )
            ALLOCATE ( KDSWT(0:NPSLV) )
            ALLOCATE ( KDEWT(0:NPSLV) )
            ALLOCATE ( N2DO (0:NPSLV) )
            ALLOCATE ( N3DO (0:NPSLV) )
            KDSWT=0
            KDEWT=0
            JDSWT=0
            JDEWT=0
            N2DO =0
            N3DO =0
        END IF
           
        CALL MPI_BCAST(JDSWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
        CALL MPI_BCAST(JDEWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)

        CALL MPI_BCAST(KDSWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
        CALL MPI_BCAST(KDEWT(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
           
        CALL MPI_BCAST(N2DO(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
        CALL MPI_BCAST(N3DO(0),NPTOT,MPI_INTEGER4,0,ICOMM,IERROR)
           
        RETURN
    END SUBROUTINE 
       
