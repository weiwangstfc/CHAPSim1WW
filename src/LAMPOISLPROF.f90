!**********************************************************************
! *********************************************************************
!> @par   SBR. INVPRO  :  CREATING LAMINAR POISEUILLE PROFILE           *
!> @param V30(1,N2M) [out]  cell velocity profile in 2nd direction
! *********************************************************************
!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Kui He @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> CREATING LAMINAR POISEUILLE PROFILE
!> The scale is to maintain Ubulk = 1.0.
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! ??/??/201?- Add pipe treatment, by Kui He.
!**********************************************************************
    SUBROUTINE LAMPOISLPROF
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE

        INTEGER(4) :: J, FLLG, k 
        REAL(WP)   :: Y1, Y2, PCOE, a, b, U1MEAN, U1MAXX
      
        IF(MYID.NE.0) RETURN
        IF(ICASE.EQ.IBOX3P) RETURN
        
        Y1 = HYB
        Y2 = HYT
        a = (HYT-HYB)/2.0_wp
        b = 0.0_wp
        PCOE= 3.0_WP/2.0_WP
        IF(icase .eq. Ipipec) THEN
            PCOE= 2.0_WP
            a = 2.0_WP*a
        ELSE IF(icase .EQ.IANNUL) THEN
            PCOE= 2.0_WP
            b  =  a + HYB
        ELSE
        END IF
        
        Vini(:)=0.0_WP
        DO J=1,NCL2
            IF ( YCC(J).LT.Y2 .AND. YCC(J).GT.Y1 ) THEN
                Vini(J)=(1.0_WP-((YCC(J)-b)**2)/a/a) * PCOE
            ELSE
                Vini(J)=0.0_WP
            ENDIF
        ENDDO
            
        U1MEAN = 0.0_wp
        U1MAXX = 0.0_WP
        DO J=1,NCL2
            DO K=1,NCL3
                U1MEAN=U1MEAN+Vini(j)/RCCI1(J)/DYFI(J)/DZI
                U1MAXX=DMAX1(U1MAXX,Vini(j))
            END DO
        ENDDO
        U1MEAN = U1MEAN/AREA_INLET
    
        IF(MYID==0) THEN   
            FLLG = 101
            OPEN(FLLG,FILE='CHK_INIL_POISEUILLE.dat')
            WRITE(FLLG,'(A)') '#  YCC, Uini'
            DO J=1,NCL2
                WRITE(FLLG,'(I8,2ES15.7)') J,YCC(J),Vini(J)
            END DO
            CLOSE(FLLG)
            
            call CHKRLHDL  ('   The given lamimar Poiseuille flow bulk velocity=        ',MYID,U1MEAN)
            call CHKRLHDL  ('   The given lamimar Poiseuille flow maxi velocity=        ',MYID,U1MAXX)
        END IF 
        
        RETURN
    
    END  SUBROUTINE
    
    SUBROUTINE BCAST_LAMPOISLPROF
        use init_info
        use mesh_info
        use flow_info
        IMPLICIT NONE
        
        CALL MPI_BCAST( Vini, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        RETURN
    END SUBROUTINE
    
    
