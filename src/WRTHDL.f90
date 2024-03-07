!>**********************************************************************
!>@par SUBROUTINE  ERRORHANDLE
!>     Write out the error message
!***********************************************************************
!> @author 
!>  Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Write out message
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Wei Wang
!**********************************************************************
       SUBROUTINE ERRHDL(msg,rank)
       IMPLICIT NONE
       
       CHARACTER*(*) msg
       INTEGER       rank
       
   11  FORMAT(A,I3.1,5X,A)
       
       write(*,11) '# Error Msg in MYID ', rank, msg
       
       STOP
       
       END SUBROUTINE ERRHDL

!>**********************************************************************
!>@par SUBROUTINE  CHECKHANDLE
!>     Write out the check message
!>**********************************************************************           
       SUBROUTINE CHKHDL(msg,rank)
       IMPLICIT NONE
       
       CHARACTER*(*) msg
       INTEGER       rank
       
   11  FORMAT(A,I3.1,5X,A)
       
       write(*,11) '# MYID ', rank, msg
       
       END SUBROUTINE
       
!>**********************************************************************
!>@par SUBROUTINE  CHECKHANDLE
!>     Write out the check message
!>**********************************************************************            
       SUBROUTINE CHKINTHDL(msg,rank,n)
       IMPLICIT NONE
       
       CHARACTER*(*) msg
       INTEGER       rank
       INTEGER       n
       
   11  FORMAT(A,I3.1,5X,A,5X,I11.1)
       
       write(*,11) '# MYID ', rank, msg,n
       
       END SUBROUTINE
       
       
       SUBROUTINE CHK2INTHDL(msg,rank,n1,n2)
       IMPLICIT NONE
       
       CHARACTER*(*) msg
       INTEGER       rank
       INTEGER       n1,n2
       
   11  FORMAT(A,I3.1,5X,A,5X,2I11.1)
       
       write(*,11) '# MYID ', rank, msg,n1,n2
       
       END SUBROUTINE
       
!>**********************************************************************
!>@par SUBROUTINE  CHECKHANDLE
!>     Write out the check message
!>**********************************************************************            
       SUBROUTINE CHKRLHDL(msg,rank,a)
       IMPLICIT NONE
       
       CHARACTER*(*)        msg
       INTEGER              rank
       DOUBLE PRECISION    a
       
   11  FORMAT(A,I3.1,5X,A,5X,ES17.9)
       
       write(*,11) '# MYID ', rank, msg,a
       
       END SUBROUTINE
       
       SUBROUTINE CHK2RLHDL(msg,rank,a1,a2)
       IMPLICIT NONE
       
       CHARACTER*(*)        msg
       INTEGER              rank
       DOUBLE PRECISION    a1, a2
       
   11  FORMAT(A,I3.1,5X,A,5X,2ES17.9)
       
       write(*,11) '# MYID ', rank, msg,a1, a2
       
       END SUBROUTINE
!*******************************************************************************************************
        SUBROUTINE DEBUG_WRT_QP_tg
            USE FLOW_INFO
            USE INIT_INFO
            USE MESH_INFO
            IMPLICIT NONE
            
            
            INTEGER(4) :: I, J, K, JJ, FLID
            CHARACTER(4) :: NNN
            
            WRITE(NNN,'(1I4.4)') myid
            FLID = MYID + 20

            OPEN(FLID,FILE='OUT_UP'//NNN//'.debug',position="append")
            WRITE(FLID,'(A)') '#=====================check======================='
            DO J=0,N2DO(MYID)+1
                DO K=1,NCL3
                    DO I=1,NCL1_TG
                        WRITE(FLID,*) 'UGP',K, I, Q_tg(I,J, K,1:3), PR_tg(I,J,K)
                    END DO
                END DO
            END DO
            
            CLOSE(FLID)
            
            RETURN
        END SUBROUTINE
!***********************************************************************   
    SUBROUTINE DEBUG_WRT_UWP_io
            USE FLOW_INFO
            USE INIT_INFO
            USE MESH_INFO
            IMPLICIT NONE
            
            
            INTEGER(4) :: I, J, K, JJ, FLID
            CHARACTER(4) :: NNN
            
            CALL MPI_BARRIER(ICOMM,IERROR)
            WRITE(NNN,'(1I4.4)') myid
            
            FLID = MYID + 200
            OPEN(FLID,FILE='OUT_UGP'//NNN//'.debug',position="append")
            
            WRITE(FLID,'(A)') '#=====================check======================='
            IF(MYID.LT.(NPTOT/2)) THEN
                DO J=1,N2DO(MYID)
                    WRITE(FLID,'(A,I4.1)') 'J= ', J
                    DO K=1,NCL3
                        DO I=NCL1S,NCL1E
                            WRITE(FLID,*) 'UGP',K, I, Q_IO(I,J, K,1:3), G_IO(I,J, K,1:3),PR_IO(I,J,K)
                        END DO
                    END DO
                END DO
            ELSE
                DO J=N2DO(MYID),1
                    WRITE(FLID,'(A,I4.1)') 'J= ', J
                    DO K=1,NCL3
                        DO I=NCL1S,NCL1E
                            WRITE(FLID,*) 'UGP',K, I, Q_IO(I,J, K,1:3), G_IO(I,J, K,1:3),PR_IO(I,J,K)
                        END DO
                    END DO
                END DO

            END IF
            
            CALL MPI_BARRIER(ICOMM,IERROR)
            CLOSE(FLID)
            
            RETURN
        END SUBROUTINE
        

!***********************************************************************   
        SUBROUTINE DEBUG_WRT_LOCAL(VAR,NJS,NJE,STR)
            USE FLOW_INFO
            USE INIT_INFO
            USE MESH_INFO
            IMPLICIT NONE
            
            INTEGER(4),INTENT(IN) :: NJS,NJE
            REAL(WP),INTENT(IN)   :: VAR(NCL1S:NCL1E,NJS:NJE,1:NCL3)
            CHARACTER(4)          :: STR
            INTEGER(4)            :: FLID
            
            INTEGER(4) :: I, J, K, JJ
            CHARACTER(4) :: NNN
            
            CALL MPI_BARRIER(ICOMM,IERROR)
            WRITE(NNN,'(1I4.4)') myid
            
            FLID = MYID+100
            OPEN(FLID,FILE='OUT_VARS'//NNN//'.debug',position="append")
            
            WRITE(FLID,'(A)') '#=====================check======================='
            IF(MYID.LT.(NPTOT/2)) THEN
                DO J=NJS,NJE
                    WRITE(FLID,'(A,I4.1)') 'J= ', J
                    DO K=1,NCL3
                        DO I=NCL1S,NCL1E
                            !WRITE(MYID+100,'(A, 2I4.1, 1ES17.9)') 'VAR',K, I, VAR(I,J, K)
                            WRITE(FLID,*) STR,K, I, VAR(I,J, K)
                        END DO
                    END DO
                END DO
            ELSE
                DO J=NJE,NJS,-1
                    WRITE(FLID,'(A,I4.1)') 'J= ', J
                    DO K=1,NCL3
                        DO I=NCL1S,NCL1E
                            !WRITE(MYID+100,'(A, 2I4.1, 1ES17.9)') 'VAR', K, I, VAR(I,J, K)
                            WRITE(FLID,*) STR,K, I, VAR(I,J, K)
                        END DO
                    END DO
                END DO

            END IF
            
        
            CALL MPI_BARRIER(ICOMM,IERROR)
            CLOSE(FLID)
            
            RETURN
        END SUBROUTINE
        
!***********************************************************************           

        SUBROUTINE DEBUG_WRT_DENSITY
            USE FLOW_INFO
            USE INIT_INFO
            USE MESH_INFO
            USE thermal_info
            IMPLICIT NONE
            
            
            INTEGER(4) :: I, J, K, JJ
            CHARACTER(4) :: NNN
            
            CALL MPI_BARRIER(ICOMM,IERROR)
            WRITE(NNN,'(1I4.4)') myid
            OPEN(MYID+100,FILE='OUT_'//NNN//'.debug',position="append")
            
            WRITE(MYID+100,'(A,1ES13.5)') 'check: t=',phyTime
            DO J=0,N2DO(MYID)+1
                DO K=1,NCL3
                    DO I=1, NCL1_IO
                        WRITE(MYID+100,*) 'DENSITY', MYID, J, K, I, DENSITY(I,J, K), DENSITY0(I,J, K)
                    END DO
                END DO
            END DO
        
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            RETURN
        END SUBROUTINE
!***********************************************************************         
      !***********************************************************************   
    SUBROUTINE DEBUG_WRT_VG2_io
            USE FLOW_INFO
            USE INIT_INFO
            USE MESH_INFO
            IMPLICIT NONE
            
            
            INTEGER(4) :: I, J, K, JJ
            CHARACTER(4) :: NNN
            
            
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            WRITE(NNN,'(1I4.4)') myid
            OPEN(MYID+200,FILE='OUT_VG2A'//NNN//'.debug',position="append")

            WRITE(MYID+200,'(A)') '=====================check======================='
            
            
            IF(MYID.LT.(NPTOT/2)) THEN
                DO J=1,N2DO(MYID)+1,+1
                    WRITE(MYID+200,'(A,I4.1)') 'J= ', J
                    DO K=1,NCL3
                        DO I=NCL1S,NCL1E
                            WRITE(MYID+200,'(A, 2I4.1, 2ES15.7)') 'VG2',K, I, Q_IO(I,J, K,2),G_IO(I,J, K,2)
                        END DO
                    END DO
                END DO
            ELSE
                DO J=N2DO(MYID)+1,1,-1
                    WRITE(MYID+200,'(A,I4.1)') 'J= ', J
                    DO K=1,NCL3
                        DO I=NCL1S,NCL1E
                            WRITE(MYID+200,'(A, 2I4.1, 2ES15.7)') 'VG2',K, I, Q_IO(I,J, K,2),G_IO(I,J, K,2)
                        END DO
                    END DO
                END DO

            END IF
            
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            RETURN
    END SUBROUTINE
    
    
    
    SUBROUTINE DEBUG_WRT_THERMAL
        USE FLOW_INFO
        USE INIT_INFO
        USE MESH_INFO
        USE thermal_info
        IMPLICIT NONE
        
        
        INTEGER(4) :: I, J, K, JJ, FLID
        CHARACTER(4) :: NNN
        CHARACTER(5)          :: STR
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        WRITE(NNN,'(1I4.4)') myid
        STR='HDTMK'
        
        FLID = MYID+100
        OPEN(FLID,FILE='OUT_VARS'//NNN//'.debug',position="append")
        
        WRITE(FLID,'(A)') '#=====================check======================='
        IF(MYID.LT.(NPTOT/2)) THEN
            DO J=0,N2DO(MYID)+1
                WRITE(FLID,'(A,I4.1)') 'J= ', J
                DO K=1,NCL3
                    DO I=NCL1S,NCL1E
                        !WRITE(MYID+100,'(A, 2I4.1, 1ES17.9)') 'VAR',K, I, VAR(I,J, K)
                        WRITE(FLID,*) STR,K, I, ENTHALPY(I,J, K), DENSITY(I,J, K), TEMPERATURE(I,J, K), &
                                       VISCOUSITY(I,J,K), THERMCONDT(I,J,K)
                    END DO
                END DO
            END DO
        ELSE
            DO J=N2DO(MYID)+1,0,-1
                WRITE(FLID,'(A,I4.1)') 'J= ', J
                DO K=1,NCL3
                    DO I=NCL1S,NCL1E
                        !WRITE(MYID+100,'(A, 2I4.1, 1ES17.9)') 'VAR', K, I, VAR(I,J, K)
                        WRITE(FLID,*) STR,K, I, ENTHALPY(I,J, K), DENSITY(I,J, K), TEMPERATURE(I,J, K), &
                                       VISCOUSITY(I,J,K), THERMCONDT(I,J,K)
                    END DO
                END DO
            END DO

        END IF
        
    
        CALL MPI_BARRIER(ICOMM,IERROR)
        CLOSE(FLID)
        
        RETURN
    END SUBROUTINE
        
        
    subroutine mkdir
        use wrt_info
        implicit none
        
        call system('mkdir -p '//dir1)
        call system('mkdir -p '//dir2)
        call system('mkdir -p '//dir3)
        call system('mkdir -p '//dir4)
        call system('mkdir -p '//dir5)
        call system('mkdir -p '//dir6)
        call system('mkdir -p '//dir7)
        
        return
    end subroutine


    subroutine wrt_3d_pt_debug(var, str, loc, iter, irk)
        use init_info
        use flow_info
        use postprocess_info
        use mesh_info
        use wrt_info
        implicit none 
        
        real(WP), intent(in)     :: var(1:NCL1_io, 1:N2DO(myid), 1:NCL3)
        character(*), intent(in) :: str
        character(*), intent(in) :: loc
        integer, intent(in)      :: iter, irk

        integer, parameter :: NPT = 4
        integer, parameter :: nfil = 20
        INTEGER  :: nid(4, 3), a(12)

        CHARACTER(1) :: PNTIM
        character(128) :: FLNM
        logical :: file_exists
        integer :: n, i, j, k, jj

        a = (/8, 16, 32, 40, 8, 16, 32, 40, 8, 16, 32, 40/)
        nid = RESHAPE(a, (/4, 3/))


        !write(*, *) 'var sz', ubound(var, 1), ubound(var, 2), ubound(var, 3)

        do n = 1, NPT
            
            WRITE(PNTIM,'(I1.1)') n
            FLNM = 'code1ww_p'//PNTIM//'_'//trim(str)//'.dat'
            

            DO J = 1, N2DO(MYID)
                JJ = JCL2G(J)
                if(JJ == nid(n, 2)) then
                    DO I = 1, NCL1_IO
                        if(I == nid(n, 1)) then
                            DO K =1, NCL3
                                if(K == nid(n, 3)) then
                                   file_exists = .false.
                                   INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
                                   IF(file_exists) THEN
                                     OPEN(nfil+n,FILE=TRIM(ADJUSTL(FLNM)), position='append')
                                     !write(nfil+n,*) '# iter = ', ITERG
                                    ELSE
                                     OPEN(nfil+n,FILE=TRIM(ADJUSTL(FLNM)) )
                                     !write(nfil+n,*) '# iter = ', ITERG
                                    END IF
                                    !write(*, *) I, JJ, K, xnd_io(i), Xcc_io(i), ynd(jj), ycc(jj), znd(k), zcc(k)
                                    write(nfil+n, *) trim(str), &
                                    trim(loc), iter, irk, i, jj, k, var(i, j, k)
                                    close(nfil+n)
                                end if
                            end do
                        end if
                    END DO
                end if
            END DO
            
        end do

    return

    end subroutine

subroutine wrt_3d_all_debug(var, str, loc, iter, irk)
    use init_info
    use flow_info
    use postprocess_info
    use mesh_info
    use wrt_info
    implicit none 
    
    real(WP), intent(in)     :: var(1:NCL1_io, 1:N2DO(myid), 1:NCL3)
    character(*), intent(in) :: str
    character(*), intent(in) :: loc
    integer, intent(in)      :: iter, irk

    integer, parameter :: nfil = 20
    character(1) :: pntim
    character(128) :: FLNM
    logical :: file_exists
    integer :: n, i, j, k, jj


    write(pntim,'(i1.1)') myid
    FLNM = 'code1ww_'//trim(str)//'_'//pntim//'.dat'
    file_exists = .false.
    INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
    IF(file_exists) THEN
        OPEN(nfil+n,FILE=TRIM(ADJUSTL(FLNM)), position='append')
        !write(nfil+n,*) '# iter = ', ITERG
    ELSE
        OPEN(nfil+n,FILE=TRIM(ADJUSTL(FLNM)) )
        !write(nfil+n,*) '# iter = ', ITERG
    END IF
        
    DO J = 1, N2DO(MYID)
        JJ = JCL2G(J)
        DO K =1, NCL3
            DO I = 1, NCL1_IO    
                write(nfil, *) JJ, K, I, var(i, j, k)
            end do
        END DO
    END DO
    close(nfil)

return

end subroutine

