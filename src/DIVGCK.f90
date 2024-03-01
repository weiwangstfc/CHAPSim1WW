!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> CALCULATION OF DIVG(Q) AND CONTINUTY CHECK 
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
!**********************************************************************
    SUBROUTINE DIVGCK_tg
        use mesh_info
        use flow_info
        IMPLICIT NONE
                     
        REAL(WP) :: DIVX
        REAL(WP) :: QMA
        REAL(WP) :: QMAX_WORK
        REAL(WP) :: DQCAP
        INTEGER(4) :: K, KP
        INTEGER(4) :: J, JP, JJ
        INTEGER(4) :: I, IP
          
        DIVX      =0.0_WP
        QMA       =0.0_WP
        QMAX_WORK = 0.0_WP
        DQCAP     = 0.0_WP
        DO K=1,NCL3
            KP=KPV(K)
            DO J=1,N2DO(MYID)
                JP=JLPV(J)
                JJ=JCL2G(J)
                DO I=1,NCL1_tg
                    IP=IPV_tg(I)
                    DQCAP=(Q_tg(IP,J,K,1)-Q_tg(I,J,K,1))*DXI       &
                         +(Q_tg(I,JP,K,2)-Q_tg(I,J,K,2))*DYFI(JJ)*RCCI1(jj)     &
                         +(Q_tg(I,J,KP,3)-Q_tg(I,J,K,3))*DZI*RCCI2(jj)
                    DIVX=DMAX1(DABS(DQCAP),DIVX)
                ENDDO
            ENDDO
        ENDDO
      
        QMA=DIVX  

        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(QMA, QMAX_WORK, 1, MPI_DOUBLE_PRECISION,  &
             MPI_MAX, ICOMM, IERROR)
        
        MAXDIVGV_tg(1) =  MAXDIVGV_tg(2)
        MAXDIVGV_tg(2) =  QMAX_WORK
      
        RETURN
    END SUBROUTINE
      
!***************************************************************************      
    SUBROUTINE DIVGCK_io
        use thermal_info
        use mesh_info
        use flow_info
        use init_info !test
        IMPLICIT NONE
                     
        !REAL(WP) :: DIVX1, DIVX2, DIVX3
        REAL(WP) :: QMA1, QMA2, QMA3
        REAL(WP) :: QMAX_WORK1, QMAX_WORK2, QMAX_WORK3
      
        QMA1       = 0.0_WP
        QMAX_WORK1 = 0.0_WP
       
        QMA2       = 0.0_WP
        QMAX_WORK2 = 0.0_WP
       
        QMA3       = 0.0_WP
        QMAX_WORK3 = 0.0_WP
       
        IF(TGFLOWflg) THEN
            CALL DIVGCK_comm_io(2,NCL1_io-1,    QMA1)
            CALL DIVGCK_comm_io(0,1,            QMA2)
            CALL DIVGCK_comm_io(NCL1_io,NCL1_io,QMA3)
        ELSE
            CALL DIVGCK_comm_io(1,NCL1_io,    QMA1)
            QMA2 =  QMA1
            QMA3 =  QMA1
        END IF
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(QMA1, QMAX_WORK1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)  
        CALL MPI_ALLREDUCE(QMA2, QMAX_WORK2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(QMA3, QMAX_WORK3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)

        MAXDIVGV_io(1) =  QMAX_WORK1
        MAXDIVGV_io(2) =  QMAX_WORK2
        MAXDIVGV_io(3) =  QMAX_WORK3
        
        !IF(myid==0) WRITE(*,'(A,3ES15.7)') '#Continuity Eq.=', QMAX_WORK1, QMAX_WORK2, QMAX_WORK3
        
        !=============test below===================
        IF(TGFLOWflg) THEN
            CALL DIVGCK_Q_comm_io(2,NCL1_io-1,    QMA1)
            CALL DIVGCK_Q_comm_io(0,1,            QMA2)
            CALL DIVGCK_Q_comm_io(NCL1_io,NCL1_io,QMA3)
        ELSE
            CALL DIVGCK_Q_comm_io(1,NCL1_io,    QMA1)
            QMA2 =  QMA1
            QMA3 =  QMA1
        END IF
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(QMA1, QMAX_WORK1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)  
        CALL MPI_ALLREDUCE(QMA2, QMAX_WORK2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(QMA3, QMAX_WORK3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        
       !IF(myid==0) WRITE(*,'(A,3ES15.7)') '#Divergence of V=', QMAX_WORK1, QMAX_WORK2, QMAX_WORK3
        !=============test above==================
        
        RETURN
    END SUBROUTINE
    
    
!==================================================================================================    
    SUBROUTINE DIVGCK_comm_io(IS,IE,DIVMAX)
        use thermal_info
        use mesh_info
        use flow_info
        use init_info !test
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: IS
        INTEGER(4),INTENT(IN) :: IE
        REAL(WP),INTENT(OUT)  :: DIVMAX
                     
        INTEGER(4) :: KC, KP
        INTEGER(4) :: JC, JP, JJ
        INTEGER(4) :: IC, IP
        REAL(WP)    :: DIVX, DIVY, DIVZ
        REAL(WP)    :: DTRHO, DQCAP
        REAL(WP)    :: COE2, COE3, COE4
      
        DIVMAX     = 0.0_WP
!        IF(MYID==0) write(*,*) 'check:'
        
        COE4= 1.0_wp/ DT    
        DO JC=1,N2DO(MYID)
            JP  = JLPV(JC)
            JJ  = JCL2G(JC)
            COE2= DYFI(JJ)*RCCI1(jj)
            COE3= DZI*RCCI2(jj) 
            DO KC=1,NCL3
                KP=KPV(KC)
                !==================MAX. DIV IN THE MAIN DOMAIN=================
                DO IC=IS,IE
                    IP=IPV_io(IC)
                    !==========d(\rho u)/dx at (i,j,k)===========================
                    DIVX  = ( G_IO(IP,JC,KC,1) - G_IO(IC,JC,KC,1) ) * DXI
                    !==========d(\rho v)/dy at (i,j,k)===========================
                    DIVY  = ( G_io(IC,JP,KC,2) - G_io(IC,JC,KC,2) ) * COE2
                    !==========d(\rho w)/dz at (i,j,k)===========================
                    DIVZ  = ( G_io(IC,JC,KP,3) - G_io(IC,JC,KC,3) ) * COE3
                    !==========D \RHO/ DT========================================
                    !DTRHO = ( DENSITY(IC,JC,KC) - DENSITYP(IC,JC,KC) ) *COE4 
                    DTRHO = DrhoDtP(IC,JC,KC)
                    !===========TOTAL============================================
                    DQCAP = DIVX + DIVY + DIVZ + DTRHO 
                    DIVMAX=DMAX1(DABS(DQCAP),DIVMAX)
                    
!                    IF(DABS(DQCAP).GT.1.0E-7_wp) &
!                    !IF(JJ==1 .or. JJ==NCL2) &!
!                    WRITE(*,'(A,4I3.1,8ES13.5)') 'divgck', myid, JC, KC, IC, &
!                    DIVX/DXI,DIVY/COE2,DIVZ/COE3, &
!                    DIVX/DXI+DIVY/COE2+DIVZ/COE3, &
!                    DTRHO/COE4,DQCAP, &
!                    G_io(IC,JP,KC,2),G_io(IC,JC,KC,2)
                    
                     
                ENDDO
            ENDDO
        ENDDO
        
        RETURN
    END SUBROUTINE
    
    SUBROUTINE DIVGCK_Q_comm_io(IS,IE,DIVMAX)
        use thermal_info
        use mesh_info
        use flow_info
        use init_info !test
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: IS
        INTEGER(4),INTENT(IN) :: IE
        REAL(WP),INTENT(OUT)  :: DIVMAX
                     
        INTEGER(4) :: KC, KP
        INTEGER(4) :: JC, JP, JJ
        INTEGER(4) :: IC, IP
        REAL(WP)    :: DIVX, DIVY, DIVZ
        REAL(WP)    :: DTRHO, DQCAP
        REAL(WP)    :: COE2, COE3, COE4
      
        DIVMAX     = 0.0_WP
!        IF(MYID==0) write(*,*) 'check:'
        
        COE4= 1.0_wp/ DT    
        DO JC=1,N2DO(MYID)
            JP  = JLPV(JC)
            JJ  = JCL2G(JC)
            COE2= DYFI(JJ)*RCCI1(jj)
            COE3= DZI*RCCI2(jj) 
            DO KC=1,NCL3
                KP=KPV(KC)
                !==================MAX. DIV IN THE MAIN DOMAIN=================
                DO IC=IS,IE
                    IP=IPV_io(IC)
                    !==========d(\rho u)/dx at (i,j,k)===========================
                    DIVX  = ( Q_IO(IP,JC,KC,1) - Q_IO(IC,JC,KC,1) ) * DXI
                    !==========d(\rho v)/dy at (i,j,k)===========================
                    DIVY  = ( Q_io(IC,JP,KC,2) - Q_io(IC,JC,KC,2) ) * COE2
                    !==========d(\rho w)/dz at (i,j,k)===========================
                    DIVZ  = ( Q_io(IC,JC,KP,3) - Q_io(IC,JC,KC,3) ) * COE3
                   
                    DQCAP = DIVX + DIVY + DIVZ  
                    DIVMAX=DMAX1(DABS(DQCAP),DIVMAX)
                ENDDO
            ENDDO
        ENDDO
        
        RETURN
    END SUBROUTINE
    
    
    
