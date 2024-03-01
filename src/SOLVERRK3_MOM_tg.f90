!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Kui He, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Time marching with RK3 method
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! ??/??/201?- Add pipe treatment, by Kui He.
! 06/02/2014- Subroutine structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE SOLVERRK3_MOM_tg(NS) 
        use cparam
        USE flow_info
        use Mesh_info
        IMPLICIT NONE

        INTEGER(4),INTENT(IN)  :: NS
        INTEGER(4) :: IDR
      
        !=======CALCULATE CONVECTION TERMS==================
        CALL CONVECTION_X_tg
        CALL CONVECTION_Y_tg
        CALL CONVECTION_Z_tg
      
        !=======CALCULATE THE WHOLE RHS =====================
        DO IDR=1,3
            CALL RHS_CvLpGpS_tg(NS,IDR)
            CALL MOMFA_tg(NS,IDR)
        END DO
      
        !=======CONSTRUCTING B.C. OF U,V,W,P=================
        CALL INTFC_VARS3(1,NCL1_tg,1,NCL1_tg,Q_tg)
        CALL BC_WALL_Q_tg
        
        !=======CONSTRUCTING AND SOLVING POISSION EQ.======== 
        CALL DIVG_tg(NS)               
        CALL FFT99_POIS3D_periodicxz(ITG)
      
        !=======CONSTRUCTING B.C. OF PHI=====================
        CALL INTFC_VARS1(1,NCL1_tg,1,NCL1_tg,DPH_tg)
        CALL BC_WALL_DPH_tg
        
        !=======UPDATE U,V,W,P===============================
        CALL VELOUPDT_tg(NS)  
        CALL INTFC_VARS3(1,NCL1_tg,1,NCL1_tg,Q_tg)
        CALL BC_WALL_Q_tg
      
        CALL PRCALC_tg(NS)
        CALL INTFC_VARS1(1,NCL1_tg,1,NCL1_tg,PR_tg)
        CALL BC_WALL_PR_tg
      
      
        RETURN
    END SUBROUTINE 
