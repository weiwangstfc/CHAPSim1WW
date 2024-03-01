!>**********************************************************************
!>@par SUBROUTINE start_mpi
!>     Initialize MPI
!>
!>@note OUTPUT:
!>     MYID: index of each CPU, from 0 to NPSLV
!>     SIZE: ammont of CPUs.
!>**********************************************************************
!***********************************************************************
!> @author 
!> Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Start the mpi lib
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! 06/12/2013- Initial Version, by Wei Wang
!**********************************************************************
    SUBROUTINE start_mpi
        use mpi_info
        use MPI
        IMPLICIT NONE   
         
        INTEGER  :: NDIM11
        INTEGER  :: DIMS(1)
        INTEGER  :: ICOORDS(1)
        INTEGER  :: NMINUS
        INTEGER  :: NPLUS
        INTEGER  :: IDIMS(1)
        LOGICAL        PERIODS(1),REORDER

        CALL MPI_INIT(IERROR)
        if(IERROR==1) then
            call ERRHDL('mpi_init fails!',0)
        end if
       
!>       @note FIND GLOBAL RANKING AND NUMBER OF PROCESSORS
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERROR)  
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,SIZE,IERROR)
       
        !*     PARAMETERS USED TO DEFINE PERIODIC CARTESIAN PROCESSOR:
        NDIM11=1
        DIMS(1)=SIZE
!
!*       ALLOWS TO FORCE PERIODIC NO. FOR PROCESSORS!!!(SO PROCESSOR AFTER
!*       (SIZE-1) IS 0, AND BEFORE 0 IS (SIZE-1)
        PERIODS(1)=.TRUE.   
!*       ALLOWS THE HARDWARE TO OPTIMIZE PROCESSOR ARRANGEMENT
        REORDER=.TRUE.    
!    
!>       @note HARDWARE OPTIMIZED PROCESSOR CONFIGURATOR AND NEW RANKING
        CALL MPI_CART_CREATE(MPI_COMM_WORLD,NDIM11,DIMS,PERIODS,  &
                             REORDER,ICOMM,IERROR)
!>       @note ABOVE replace MPI_COMM_WORLD by ICOMM
 
!>       GET MY POSITION IN THIS COMMUNICATOR AND MY NEIGHBORS
        CALL MPI_COMM_RANK(ICOMM,MYID,IERROR)
        CALL MPI_CART_SHIFT(ICOMM,0,1,NMINUS,NPLUS,IERROR)
        CALL MPI_CART_GET(ICOMM,1,IDIMS,PERIODS,ICOORDS,IERROR)
       
        NPTOT = SIZE
        NPSLV = SIZE-1 
       
        RETURN
    END SUBROUTINE
