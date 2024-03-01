!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani,Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Calculate the CFL number. \f$ CFL=dt * Max. SUM_{i=1}^3{U_i*x_i} < CFL_max \f$ 
!> Details about CFL is given in 
!> Ref: http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
!
!> @return CFLM = \f$ Max. SUM_{i=1}^3{U_i*x_i} \f$
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! 06/02/2014- The original uniform dy is replaced by the real dy. by Wei Wang
!**********************************************************************
    SUBROUTINE CFL_tg
        use flow_info
        use init_info
        use mesh_info
      
        IMPLICIT NONE
      
        REAL(WP)      :: CFLMA
        REAL(WP)      :: CFLM_WORK
        REAL(WP)      :: QCF!, COE1 
        INTEGER(4)  :: I, IP
        INTEGER(4)  :: J, JP, JJ, JC
        INTEGER(4)  :: K, KP

        CFLMA=0.0_WP
        CFLM_WORK = 0.0_WP
        QCF = 0.0_WP
      
        DO K=1,NCL3
            KP=KPV(K)
            DO J=1,N2DO(MYID)
                JP=JLPV(J)
                JJ=JCL2G(J)
                JC=JJ !NCL2/2
                DO I=1,NCL1_tg
                    IP=IPV_tg(I)
                    QCF=( DABS( (Q_tg(I,J,K,1)+Q_tg(IP,J,K,1))*DXI                )+        &
                          DABS( (Q_tg(I,J,K,2)+Q_tg(I,JP,K,2))*DYFI(JC)*RCCI1(JC) )+        &
                          DABS( (Q_tg(I,J,K,3)+Q_tg(I,J,KP,3))*DZI*RCCI2(JC)      ) )*0.50_WP  
                    CFLMA=DMAX1(CFLMA,QCF)
                ENDDO
            ENDDO
        ENDDO
      
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CFLMA, CFLM_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        CFLMM_tg = CFLM_WORK
      
        RETURN
    END SUBROUTINE

!**********************************************************************
    SUBROUTINE CFL_io
!>    @note
!>    1) Energy equation has the same convection velocity as the 
!>       momentum equation, thus, no additional CFL is required 
!>       for energy equation.
!>    2) For conservative variables, the convection velocity are G

        use flow_info
        use init_info
        use mesh_info
      
        IMPLICIT NONE
      
        REAL(WP)      :: CFLMA
        REAL(WP)      :: CFLM_WORK
        REAL(WP)      :: QCF   ! , COE1
        INTEGER(4)  :: I, IP
        INTEGER(4)  :: J, JP, JJ!, JC
        INTEGER(4)  :: K, KP

        CFLMA=0.0_WP
        CFLM_WORK = 0.0_WP
        QCF = 0.0_WP
    
        DO K=1,NCL3
            KP=KPV(K)
            DO J=1,N2DO(MYID)
                JP=JLPV(J)
                JJ=JCL2G(J)
                DO I=1,NCL1_io
                    IP=IPV_io(I)
                    QCF=( DABS( (G_io(I,J,K,1)+G_io(IP,J,K,1))*DXI                 )+        &
                          DABS( (G_io(I,J,K,2)+G_io(I,JP,K,2))*DYFI(JJ)*RCCI1(JJ)  )+        &
                          DABS( (G_io(I,J,K,3)+G_io(I,J,KP,3))*DZI*RCCI2(JJ)       )   )*0.50_WP  
                    QCF=QCF + CFLVIS(J)
                    CFLMA=DMAX1(CFLMA,QCF)
!>                  @warning Why rm is involved in 3rd direction?
                ENDDO
            ENDDO
        ENDDO
      
        CALL MPI_BARRIER(ICOMM,IERROR)
        CALL MPI_ALLREDUCE(CFLMA, CFLM_WORK, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ICOMM, IERROR)
        CFLMM_io = CFLM_WORK
      
        RETURN
    END SUBROUTINE
    
!**********************************************************************    
    SUBROUTINE CFL_VISCOUS
        use flow_info
        use init_info
        use mesh_info
      
        IMPLICIT NONE
        

        INTEGER(4)  :: J, JJ

        CFLVIS = 0.0_WP
        IF(visthemflg == visexplicit) THEN
            DO J=1,N2DO(MYID)
                JJ=JCL2G(J)
                CFLVIS(J)= 2.0_wp*CVISC* (DXQI + DZQI*RCCI2(JJ) + DYFI(JJ)*DYFI(JJ) )
            ENDDO
        ELSE
            CFLVIS(:) = 0.0_WP
        END IF

        RETURN
    END SUBROUTINE
        
