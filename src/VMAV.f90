!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Calculte the max. velocity
!> @param Q(local)  [in]
!> @param VMV(1:3)  [out] Max. Q (only scaled perturbulation) in all directions and all locations.
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE VMAV_tg
      
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        INTEGER(4)  :: I, J, K, JJ, NYI
        REAL(WP)     :: VMA, VMI
        REAL(WP)     :: VMAX_WORK, VMIN_WORK
      
!>       @note Max. Q(i,j,k,1)       
        VMA = -1.0e20_WP
        VMI =  1.0e20_WP
        DO K=1,NCL3
            DO J=1,N2DO(MYID)
                DO I=1,NCL1_tg              
                    VMA=DMAX1(VMA,DABS(Q_tg(I,J,K,1)))
                    VMI=DMIN1(VMI,DABS(Q_tg(I,J,K,1)))
                ENDDO
            ENDDO
        ENDDO
     
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)   
        CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ICOMM, IERROR)    
        VMAX_tg(1)=VMAX_WORK
        VMIN_tg(1)=VMIN_WORK
       
       
!>       @note Max. Q(i,j,k,3)       
        VMA = -1.0e20_WP
        VMI =  1.0e20_WP
        DO K=1,NCL3
            DO J=1,N2DO(MYID)
                JJ=JCL2G(J) 
                DO I=1,NCL1_tg              
                    VMA=DMAX1( VMA,DABS( Q_tg(I,J,K,3)*RCCI1(jj) )   )
                    VMI=DMIN1( VMI,DABS( Q_tg(I,J,K,3)*RCCI1(jj) )   )
                ENDDO
            ENDDO
        ENDDO
     
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ICOMM, IERROR)
        VMAX_tg(3)=VMAX_WORK
        VMIN_tg(3)=VMIN_WORK
       
!>        @note Max. Q(i,j,k,2)       
        VMA = -1.0e20_WP
        VMI =  1.0e20_WP
        NYI = 1
        IF(MYID==0) NYI = 2
        DO K=1,NCL3
            DO J=NYI,N2DO(MYID)
                JJ=JCL2G(J)
                DO I=1,NCL1_tg              
                    VMA=DMAX1(VMA,DABS( Q_tg(I,J,K,2)*RNDI1(jj)  ) )
                    VMI=DMIN1(VMI,DABS( Q_tg(I,J,K,2)*RNDI1(jj)  ) )
                ENDDO
            ENDDO
        ENDDO
     
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)   
        CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ICOMM, IERROR)   
        VMAX_tg(2)=VMAX_WORK
        VMIN_tg(2)=VMIN_WORK     
       
        RETURN
      
    END SUBROUTINE
    
    !********************************************************************************
    SUBROUTINE VMAV_io
        use mesh_info
        use flow_info
        IMPLICIT NONE
      
        INTEGER(4)  :: I, J, K, JJ, NYI
        REAL(WP)     :: VMA, VMI
        REAL(WP)     :: VMAX_WORK, VMIN_WORK     
      
        
        VMA = -1.0e20_WP
        VMI =  1.0e20_WP
        DO K=1,NCL3
            DO J=1,N2DO(MYID)
                DO I=NCL1S, NCL1E
                    VMA=DMAX1(VMA,DABS(Q_io(I,J,K,1))) 
                    VMI=DMIN1(VMI,DABS(Q_io(I,J,K,1))) 
                ENDDO
            ENDDO
        ENDDO
     
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MAX, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MIN, ICOMM, IERROR)
        VMAX_io(1)=VMAX_WORK
        VMIN_io(1)=VMIN_WORK
       
       
!>      @note Max. Q(i,j,k,3)       
        VMA = -1.0e20_WP
        VMI =  1.0e20_WP
        DO K=1,NCL3
            DO J=1,N2DO(MYID)
                JJ=JCL2G(J) 
                DO I=NCL1S, NCL1E
                    VMA=DMAX1(  VMA,DABS(Q_io(I,J,K,3)*RCCI1(jj))  ) 
                    VMI=DMIN1(  VMI,DABS(Q_io(I,J,K,3)*RCCI1(jj))  ) 
                ENDDO
            ENDDO
        ENDDO
     
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MAX, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MIN, ICOMM, IERROR)
        VMAX_io(3)=VMAX_WORK
        VMIN_io(3)=VMIN_WORK
       
!>     @note Max. Q(i,j,k,2)       
        VMA = -1.0e20_WP
        VMI =  1.0e20_WP
        NYI = 1
        IF(MYID==0) NYI = 2
        DO K=1,NCL3
            DO J=NYI,N2DO(MYID)
                JJ=JCL2G(J)
                DO I=NCL1S, NCL1E
                    VMA=DMAX1(VMA,DABS(Q_io(I,J,K,2)*RNDI1(jj)))
                    VMI=DMIN1(VMA,DABS(Q_io(I,J,K,2)*RNDI1(jj)))
                ENDDO
            ENDDO
        ENDDO
        
     
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(VMA, VMAX_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MAX, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VMI, VMIN_WORK, 1, MPI_DOUBLE_PRECISION,  MPI_MIN, ICOMM, IERROR)
        VMAX_io(2)=VMAX_WORK
        VMIN_io(2)=VMIN_WORK  
       
        RETURN
      
      END SUBROUTINE
      
      
      
