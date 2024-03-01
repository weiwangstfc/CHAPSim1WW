!***********************************************************************
!> @author 
!> Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!>@details 
!> Convert the local z, y index to global index 
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! 18/12/2013- Initial Version, by Wei Wang
!**********************************************************************
    SUBROUTINE INDEXL2G
        use mesh_info
        IMPLICIT NONE     

        INTEGER(4) :: J, JJ, K
        INTEGER(4),ALLOCATABLE :: CIP(:)
        INTEGER(4),ALLOCATABLE :: CLC(:)
      
        ALLOCATE ( CIP(NCL2) )
        CIP    = 0 
      
        DO J=1,N2DO(MYID)
            JCL2G(J)=JDSWT(MYID)-1+J
        END DO
      
        DO K=1,N3DO(MYID)
            KCL2G(K)=KDSWT(MYID)-1+K
        END DO
       
        DO J=1,N2DO(MYID)
            JJ = JCL2G(J)
            CIP(JJ)=MYID
        END DO
       
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CIP(1),  JG2IP(1),  NCL2, MPI_INTEGER4,MPI_SUM, ICOMM, IERROR)
        DEALLOCATE (CIP)
        
        ALLOCATE ( CLC(NCL2) )
        CLC = 0
        DO JJ=1, NCL2
            IF(MYID==JG2IP(JJ)) THEN
                CLC(JJ) = JJ-JDSWT(MYID)+1
            END IF
        END DO
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CLC(1),  JG2LC(1),  NCL2, MPI_INTEGER4,MPI_SUM, ICOMM, IERROR)
        DEALLOCATE (CLC)
        
        IF(MYID==0) THEN
            OPEN(15,FILE='CHK_INDEX_J2G.dat')
            REWIND 15
            WRITE(15,'(A)') '# JG, JIP, JLC'  
            
            DO JJ=1, NCL2
                WRITE(15,'(3I8.1)') JJ,JG2IP(JJ),JG2LC(JJ) 
            END DO
            CLOSE(15)
        END IF
        
      
        RETURN
      
    END SUBROUTINE

      
      
      
                
!**********************************************************************
