!/**
!>Document of the Code "IncompSim"

!>@section mp_sec_intro_tag Introduction
!>  "IncompSim" is a finite-different based Direct Numerical Simulation code.
  
!>@section mp_sec_devep_tag Code Development

!>@subsection mp_subsec_devep1_tag Main Developers
!>    - Dr Mehdi Seddighi
    
!>@subsection mp_subsec_devep2_tag Co-Developers
!>    - Dr Wei Wang
    
!>@subsection mp_subsec_devep3_tag Other-Developers
!>    - Akshat Mathur
!>    - Kui He
    
!>@subsection mp_subsec_hist_tag Code Development History

!>    - Mehdi Seddighi
!>        * Developped the main code of IncomSim_Main
!>        * Added the IBM method to tackle rought wall problems. (IncomSim_IBM)
        
!>    - Akshat Mathur (Based on IncomSim_Main)
!>        * Added Large Eddy Simulation Method to the code (IncomSim_LES)
        
!>    - Kui He (Based on IncomSim_Main)
!>        * Added Cylindrical Coordinates (IncomSim_pipe)
        
!>    - Wei Wang (Based on IncomSim_pipe)
!>        * Improved the code from F77 to F90 and optimized code structure/memory (IncomSim_F90)
!>        * Developped the developping flow with turbulence generator. (IncomSim_F90)
!>        * Developped the DNS with heat transfer, applied to supercritical fluids. (IncomSim_HT)


!>@subsection mp_subsec_grp_tag Group
!>    - Heat, Flow and Turbulence Research Group (http://www.sheffield.ac.uk/heft)
    
!>@subsection mp_subsec_pi_tag PI
!>    - Professor Shuisheng He
    
!>@subsection mp_subsec_cont_tag Contact
!>    - Dr Mehdi Seddighi, seddighi@sheffield.ac.uk
!>    - Dr Wei Wang, wei.wang@sheffield.ac.uk
!>    - Akshat Mathur, mep11am@sheffield.ac.uk 
!>    - Kui He, mep11kh@sheffield.ac.uk  
!>    - Professor Shuisheng He, s.he@sheffield.ac.uk 

!>@section mp_sec_tc1_tag Terms and Conditions

!>  Copyright (C) 2015  The University Of Sheffield
  
!>  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add information on how to contact you by electronic and paper mail.

!>  Head of Heat, Flow and Turbulence Research Group in the department of Mechanical Engineering, The University of sheffied, Shuisheng He, hereby disclaims all copyright interest in the CFD code 'IncompSim'.

!>   Shuisheng He (signature), 1 May, 2015
   
!>   Shuisheng He (printed name),  Head of Heat, Flow and Turbulence Research Group in Department of Mechanical Engineering, The University of sheffied
   
!>@section mp_sec_tc2_tag Terms and Conditions

!>   Copyright (C) 2015  The University Of Sheffield
   
!>   This is proprietary source code of Heat, Flow and Turbulence Research Group in the department of Mechanical Engineering, The University of sheffied. This code and any parts of it may not be used, copied or distributed by any means unless written permission is given by the Head of Heat, Flow and Turbulence Research Group in the department of Mechanical Engineering, The University of sheffied, UK. 
 
!>*/



!**********************************************************************
! DESCRIPTION: 
!> Brief description of the CFD solver 
!> @brief                                                         
!> Code Name: DNS OF UNSTEADY CHANNEL FLOW             
!> @author  
!> Mehdi Seddighi-Moornani, Kui He, Wei Wang @ HETF, ME, The University of Sheffield.
!> @version        VER CHD-2.5.0    ; 11Nov2012                
!>
!> @par        SIMULATION CODE(SINGLE PRECISION) 
!>             PARALLEL     
!>                    -  MPI                                          \n
!>             MESH                                                   \n
!>                    -  STAGGERED GRID                               \n
!>             NONLINEAR                                              \n
!>                    -  CONVECTION & DIVERGENCE FORM                 \n
!>                    -  2ND ORDER FINITE DIFFERENCE METHOD           \n
!>                    -  RUNGE KUTTA & ADAMS-BASHFORTH METHOD         \n
!>             VISCOUS                                                \n
!>                    -  CRANK-NICOLSON METHOD                        \n
!>             PRESSURE                                               \n
!>                    -  FFT                                          \n
!>                    -  FRACTIONAL STEP METHOD                       \n
!> REVISION HISTORY:
!> ??/11/2012- Initial Version, by Mehdi Seddighi
!> ??/??/201?- Add pipe treatment, by Kui He.
!> 06/02/2014- Code optimizaion by Wei Wang
! ********************************************************************

    PROGRAM DNS_THERMAL
        use init_info
        use mesh_info
        use FISHPACK_POIS3D
        IMPLICIT NONE  
                   
!>      !===============CODE STARTES=========================================================================
        CALL mkdir
        CALL START_MPI
        IF (MYID.EQ.0) CALL CHKHDL('****DNS solver starts******',myid)
        
        IF (MYID.EQ.0) CALL CHKHDL('1. Read ini file...',myid)   
        IF (MYID.EQ.0) &
        CALL READINI 
        CALL BCAST_READINI
        
        IF (MYID.EQ.0) CALL CHKHDL('2. Mesh decomposition (both y-dir and z-dir)...',myid)
        CALL mesh_Ydecomp
        CALL mesh_Zdecomp 
        
        IF (MYID.EQ.0) CALL CHKHDL('3. Allocate variables...',myid)
        CALL MEM_ALLOCAT
        
        IF (MYID.EQ.0) CALL CHKHDL('4. Set up global and local index...',myid)
        CALL INDEXL2G
        
        IF (MYID.EQ.0) CALL CHKHDL('5. Set up +1/-1 index (all three directions)...',myid)
        CALL INDEXSET 
        
        IF (MYID.EQ.0) CALL CHKHDL('6. Set up RK coefficients...',myid)  
        IF (MYID.EQ.0) &
        CALL CONS_RKCOEF 
        CALL BCAST_RKCOEF
        
        IF (MYID.EQ.0) CALL CHKHDL('7. Set up mesh information in x and z directions...',myid)    
        IF (MYID.EQ.0) &
        CALL CONSPARA   
        CALL BCAST_CONSPARA
        
        IF (MYID.EQ.0) CALL CHKHDL('8. Set up mesh distribution in y direction...',myid)
        IF (MYID.EQ.0) &
        CALL COORDJR
        CALL BCAST_COORDJR
        
        IF (MYID.EQ.0) CALL CHKHDL('9. Set up laminar (initial) velocity profile...',myid)   
        IF (MYID.EQ.0) &
        CALL LAMPOISLPROF
        CALL BCAST_LAMPOISLPROF
        
        IF (MYID.EQ.0) CALL CHKHDL('10. Set up mesh coefficient for TDMA...',myid)         
        IF (MYID.EQ.0) & 
        CALL TDMA_COEF
        CALL BCAST_TDMA_COEF
        CALL CFL_VISCOUS

        !====PREPARE POISSON SOLVER==========
        IF (MYID.EQ.0) CALL CHKHDL('11.Initialization of FFT solver...',myid) 
        IF(TGFLOWflg .AND. IOFLOWflg) THEN
            CALL FFT99_POIS3D_INIT(ITG)
            if(fishpack==1) then
              CALL FISHPACK_POIS3D_INIT(BCX_io, BCZ, NCL1_io, NCL2, NCL3, N2DO(MYID), N3DO(MYID), DXQI, DZQI, AMPH, ACPH, APPH)
            else
              CALL FISHPACK_POIS3D_INIT(BCX_io, BCZ, NCL1_io, NCL2, NCL3, N2DO(MYID), N3DO(MYID), DXQI, DZQI, AMPH, ACPH, APPH)
            end if
        ELSE
            IF(TGFLOWflg) CALL FFT99_POIS3D_INIT(ITG)
            IF(IOFLOWflg) THEN
                !CALL FFT99_POIS3D_INIT(IIO) !Method One
              if(fishpack==1) then
                CALL FISHPACK_POIS3D_INIT(BCX_io, BCZ, NCL1_io, NCL2, NCL3, N2DO(MYID), N3DO(MYID), DXQI, DZQI, AMPH, ACPH, APPH)  !Metthod Two, good
              else
                CALL FISHPACK_POIS3D_INIT(BCX_io, BCZ, NCL1_io, NCL2, NCL3, N2DO(MYID), N3DO(MYID), DXQI, DZQI, AMPH, ACPH, APPH)  !Metthod Two, good
                !CALL FFT99_POIS3D_INIT(IIO)
              end if
            END IF
        END IF
        
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !====PREPARE THERMAL INFO============
        IF(IOFLOWflg) THEN
            !=================thermal property setup=============================
            IF(thermlflg==1) THEN
                IF (MYID.EQ.0) CALL CHKHDL('12. Initialization of thermal properties...', myid) 
                IF (MYID.EQ.0) CALL thermal_init
                CALL BCAST_THERMAL_INIT
                
                CALL MPI_BARRIER(ICOMM,IERROR)
            END IF
        END IF
       
        IF(pprocessonly==2) THEN
            CALL POSTPROCESS_INTEGRAL_INSTANS
            CALL MPI_BARRIER(ICOMM,IERROR)
            IF(MYID==0) CALL CHKHDL('<===Only postprocessed given instantanous results, now the code stops...==>', myid)  
            STOP
        END IF
       
        !====The Main CFD Solver==========
        IF (MYID.EQ.0) CALL CHKHDL('13.The main CFD solver starts...',myid)       
        CALL SOLVE     
        
        CALL MEM_DEALLOCAT     
        CALL MPI_FINALIZE(IERROR)
        STOP
       
    END PROGRAM
