!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details
!> setup some constant parameter for RK method
!
!> @todo
!> setup coefficients for RK3 and Adams-Bashfort(notusedhere)           
!> RUNGE KUTTA & ADAMS-BASHFORT & SEMI IMPLICIT PARAMETERS
!> it's better to initialize NSST from ini file
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
!**********************************************************************
    SUBROUTINE CONS_RKCOEF
        use init_info
        IMPLICIT NONE

        IF(MYID.NE.0) RETURN
        
        NSST=3
        IF (NSST.EQ.3) THEN
            IF (MYID.EQ.0) call CHKHDL('         Time scheme: R-K 3rd order',MYID)
!             RUNGE KUTTA Coefficients
            TGAM(0)=   1.0_WP
            TGAM(1)=   8.0_WP/15.0_WP
            TGAM(2)=   5.0_WP/12.0_WP
            TGAM(3)=   3.0_WP/4.0_WP
           
            TROH(0)=   0.0_WP
            TROH(1)=   0.0_WP
            TROH(2)= -17.0_WP/60.0_WP
            TROH(3)=  -5.0_WP/12.0_WP
        ELSE
!            IF NSST=1, ADAMS-BASHFORT SCHEME
            IF (MYID.EQ.0) call CHKHDL('         Time scheme: ADAMS-BASHFORT 3rd order',MYID)
            TGAM(0)= 1.0_WP
            TGAM(1)= 1.5_WP
            TGAM(2)= 0.0_WP
            TGAM(3)= 0.0_WP
           
            TROH(0)= 0.0_WP
            TROH(1)=-0.5_WP
            TROH(2)= 0.0_WP
            TROH(3)= 0.0_WP
        ENDIF
        
        TALP(0)=  TGAM(0) + TROH(0)
        TALP(1)=  TGAM(1) + TROH(1)
        TALP(2)=  TGAM(2) + TROH(2)
        TALP(3)=  TGAM(3) + TROH(3)

        RETURN
    END SUBROUTINE
    
    
    SUBROUTINE BCAST_RKCOEF
        use init_info
        IMPLICIT NONE
        
        CALL MPI_BCAST( NSST,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( TGAM(0), 4, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( TROH(0), 4, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( TALP(0), 4, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        RETURN
    END SUBROUTINE
    
