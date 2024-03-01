!***********************************************************************************
    SUBROUTINE thermal_init ! only in master
        USE THERMAL_INFO
        USE MESH_INFO
        USE INIT_INFO
        IMPLICIT NONE
        
        INTEGER(4)  :: IE, IS,  NIN1, NIN2, I
        REAL(wp)    :: massflux
        REAL(WP)   :: H_temp, M_temp, K_temp, D_temp, Cp_temp, P_temp, B_temp, T_temp
        
!=====================reference state===(dimensional)======different for perfect gas and supercriticla flow==========
        IF(thermoStat==idealgas_law) THEN ! for perfect gas
            ! The perfect gas is only thermal perfect gas, not calorically perfect gas.
            ! Thus, Cp and Cv changes with Time.
            
        !=========constants=================
            ! ref: http://en.wikipedia.org/wiki/Ideal_gas
            !knowns: T0, P0
            ! R=gas constant (J/Kg/K), p=rho R T
            R0 = 287.058_wp ! (J/Kg/K)
            
            !D0=density, equation of state for perfect gas
            D0 = P0/R0/T0
            
            !gamma0 = heat capacity ratio.
            !Cvhat0 = the dimensionless specific heat capacity at constant volume
            SELECT CASE (IdealGasType)
                CASE (monatomic_gas)      ! one atom in molecule
                    gamma0 = 5.0_wp/3.0_wp
                    Cvhat0 = 3.0_wp/2.0_wp
                    !PRT0   = 0.67_wp
                CASE (diatomic_gas)       ! two atoms in molecule (like, air)
                    gamma0 = 7.0_wp/5.0_wp
                    Cvhat0 = 5.0_wp/2.0_wp
                    !PRT0   = 0.73_wp
                CASE (trivalence_gas)
                    gamma0 = 1.3_wp
                    Cvhat0 = 3.0_wp
                    !PRT0   = 0.8_wp
                CASE DEFAULT
                    gamma0 = 7.0_wp/5.0_wp
                    Cvhat0 = 5.0_wp/2.0_wp
                    !PRT0   = 0.73_wp
            END SELECT
            
        !====calculation=================
            !===specific heat capacity (J/Kg/K)
            Cv0 = R0/(gamma0 - 1.0_wp)
            Cp0 = Cv0*gamma0
            !===internal energy/enthalpy (J/Kg)
            !e0 = Cvhat0 * R0 * T0
            H0 = (Cvhat0 + 1.0_wp)* R0 * T0
            
            !==thermal expansion ratio = 1/T====
            B0 = 1.0_wp/T0
            
            !molecular viscousity from the sutherland equation for perfect gas
            
            CALL Sutherland_Viscosity_dim(T0,M0)
            CALL Sutherland_Conductivity_dim(T0,K0)
            
            ! thermal conductivity from Prandlt no.
            PRT0 = M0*Cp0/K0
            
!            H0  = SFE*NIST_H(IE)  + SFS*NIST_H(IS)
!            D0  = SFE*NIST_D(IE)  + SFS*NIST_D(IS)
!            K0  = SFE*NIST_K(IE)  + SFS*NIST_K(IS)
!            M0  = SFE*NIST_M(IE)  + SFS*NIST_M(IS) 
!            CP0 = SFE*NIST_CP(IE) + SFS*NIST_CP(IS) ! updated
!            B0  = SFE*NIST_B(IE)  + SFS*NIST_B(IS)
            
            !========power lawer coefficients==========
            PL_MT  = 0.67_WP
            PL_CpT = 0.095_WP
            PL_HT  = PL_CpT+1.0_WP !1.095
            PL_KT  = 0.805_WP
        
        ELSE IF(thermoStat==search_table) THEN
            CALL NIST_TABLE_SCALING

            PRT0   = M0*CP0/K0
        ELSE
        END IF
        
!====================COMMON  FOR supercritical flow and perfect gas=====================
        CTHECD = 1.0_WP/REN/PRT0
        U0     = REN*M0/D0/L0
        
        RHOU20 = D0*U0*U0
        
        IF(GRAVDIR==0) THEN !not consider gravitational acceleration
            G_A = 0.0_WP           
        ELSE                !gravitational acceleration
            G_A = 9.80665_WP    
        END IF
        
        IBuoF = 0
        IF(FLOWDIR==0) THEN       ! HORIZONTAL FLOW
            F_A = -L0/U0/U0*G_A
            IBuoF(2)=1
        ELSE if(FLOWDIR==1) THEN  ! UPWORD FLOW
            F_A = -L0/U0/U0*G_A
            IBuoF(1)=1
        ELSE if(FLOWDIR==2) THEN  ! DOWNWARD FLOW
            F_A =  L0/U0/U0*G_A
            IBuoF(1)=1
        ELSE 
            F_A =  0.0_WP
        END IF
        
!===============Wall info=================================================
        !=====================isoflux========================
        IF( BCWALLHEAT(itopwall)==isoFluxWall .or. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
        
            ALLOCATE (WALLFLUX(NCL1S:NCL1E,2)  ); WALLFLUX = 0.0_wp
            
            !carefully in ini file 
            ! bottom wall heating (heat flux in ), qw = -kdT/dy = positive
            ! bottom wall cooling (heat flux out), qw = -kdT/dy = negtive
            ! top    wall heating (heat flux in ), qw = -kdT/dy = negtive
            ! top    wall cooling (heat flux out), qw = -kdT/dy = positive 
            
            IF(BCWALLHEAT(itopwall)==isoFluxWall) THEN
            
                WHEAT0_UNDIM(itopwall) = WHEAT0_DIM(itopwall)*L0/K0/T0
                NIN2 = 0
                DO I=NCL1S, NCL1E
                    !=====TOP WALL HEATING====
                    IF(I.LT.NIN2) THEN
                        WALLFLUX(I,itopwall) = 0.0_WP
                    ELSE
                        WALLFLUX(I,itopwall) = -WHEAT0_UNDIM(itopwall)
                    END IF
                END DO
                
            END IF
            IF(ICASE.ne.Ipipec .and. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
            
                WHEAT0_UNDIM(ibotwall) = WHEAT0_DIM(ibotwall)*L0/K0/T0
                NIN1 = 0
                DO I=NCL1S, NCL1E
                    !=====TOP WALL HEATING====
                    IF(I.LT.NIN1) THEN
                        WALLFLUX(I,ibotwall) = 0.0_WP
                    ELSE
                        WALLFLUX(I,ibotwall) = WHEAT0_UNDIM(ibotwall)
                    END IF
                END DO
            END IF
            
            IF( BCWALLHEAT(itopwall)==isoFluxWall .and. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
                Gr0(1:2)   = G_A*B0*WHEAT0_DIM(1:2)*(L0**4)/K0/((M0/D0)**2)
                BO0(1,1:2) = Gr0(1:2) / REN**(2.7_wp)
                BO0(2,1:2) = Gr0(1:2) / (REN**3.425_wp) / (PRT0**0.8_wp)
            END IF
        END IF
        
        
        !=====================isothermal========================
        IF( BCWALLHEAT(itopwall)==isoThermalWall .or. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
            ALLOCATE (T_WAL_GV(NCL1S:NCL1E,2)  ); T_WAL_GV = 0.0_wp
            ALLOCATE (H_WAL_GV(NCL1S:NCL1E,2)  ); H_WAL_GV = 0.0_wp
            ALLOCATE (D_WAL_GV(NCL1S:NCL1E,2)  ); D_WAL_GV = 0.0_wp
            ALLOCATE (M_WAL_GV(NCL1S:NCL1E,2)  ); M_WAL_GV = 0.0_wp
            ALLOCATE (K_WAL_GV(NCL1S:NCL1E,2)  ); K_WAL_GV = 0.0_wp
            ALLOCATE (Cp_WAL_GV(NCL1S:NCL1E,2) ); Cp_WAL_GV = 0.0_wp
            
            
            IF(BCWALLHEAT(itopwall)==isoThermalWall) THEN
                WHEAT0_UNDIM(itopwall) = WHEAT0_DIM(itopwall)/T0
                NIN2 = 0
                DO I=NCL1S, NCL1E
                    IF(I.LT.NIN2) THEN
                        T_WAL_GV(I,itopwall) = WHEAT0_UNDIM(itopwall)
                    ELSE
                        T_WAL_GV(I,itopwall) = WHEAT0_UNDIM(itopwall)
                    END IF
                    !WRITE(*,*) 'bottomT',I, WALLFLUX(I,1) !test
                END DO
                CALL WALL_ISOTHERMAL_GIVEN(itopwall)
            END IF
            IF(ICASE.ne.Ipipec .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                WHEAT0_UNDIM(ibotwall) = WHEAT0_DIM(ibotwall)/T0
                NIN1 = 0
                DO I=NCL1S, NCL1E
                    IF(I.LT.NIN1) THEN
                        T_WAL_GV(I,ibotwall) = WHEAT0_UNDIM(ibotwall)
                    ELSE
                        T_WAL_GV(I,ibotwall) = WHEAT0_UNDIM(ibotwall)
                    END IF
                END DO
                CALL WALL_ISOTHERMAL_GIVEN(ibotwall)
            END IF
        
            IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                Gr0(1:2)   = G_A*B0*DABS(WHEAT0_DIM(1)-WHEAT0_DIM(2))*((2.0_wp*L0)**3)/((M0/D0)**2)
                BO0(1,1:2) = Gr0(1:2) / (2.0*REN)**(2.7_wp)
                BO0(2,1:2) = Gr0(1:2) / ((2.0*REN)**3.425_wp) / (PRT0**0.8_wp)
            END IF
        END IF
        
        massflux = REN*M0/L0
      
        IF(MYID==0) THEN
        
            IF(thermoStat==search_table) THEN
                OPEN(10,FILE='CHK_NIST_TABLE_UNDIM.dat')
                WRITE(10,'(A,1ES13.5)') '# REF T0(K)=             ', T0
                WRITE(10,'(A,1ES13.5)') '# REF H0(J)=             ', H0
                WRITE(10,'(A,1ES13.5)') '# REF D0(Kg/m3)=         ', D0
                WRITE(10,'(A,1ES13.5)') '# REF M0(Pa-s)=          ', M0
                WRITE(10,'(A,1ES13.5)') '# REF K0(W/m-K)=         ', K0
                WRITE(10,'(A,1ES13.5)') '# REF B0(1/K)=           ', B0
                WRITE(10,'(A,1ES13.5)') '# REF massflux(Kg/m2s)=  ', B0
                WRITE(10,'(A,1ES13.5)') '# REF CP0=               ', CP0
                WRITE(10,'(A,1ES13.5)') '# REF H0/CP0/T0=         ', H0/CP0/T0
                
                IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                    WRITE(10,'(A,A,A)')     '#          ', 'BOTTOM WALL', ' TOP WALL'
                    WRITE(10,'(A,2ES13.5)') '# T=       ', T_WAL_GV(1,1), T_WAL_GV(1,2)
                    WRITE(10,'(A,2ES13.5)') '# H=       ', H_WAL_GV(1,1), H_WAL_GV(1,2)
                    WRITE(10,'(A,2ES13.5)') '# D=       ', D_WAL_GV(1,1), D_WAL_GV(1,2)
                    WRITE(10,'(A,2ES13.5)') '# M=       ', M_WAL_GV(1,1), M_WAL_GV(1,2)
                    WRITE(10,'(A,2ES13.5)') '# K=       ', K_WAL_GV(1,1), K_WAL_GV(1,2)
                END IF
            
                WRITE(10,'(A)') '#P(Mpa)      H     T    D    M     K      CP      BETA'
                WRITE(10,*) '# ', N_NIST
            
                DO I=1,N_NIST
                    WRITE(10,'(8ES13.5)') P0, &
                                NIST_H(I), NIST_T(I), NIST_D(I), &
                                NIST_M(I), NIST_K(I), NIST_CP(I),NIST_B(I)
                END DO
                CLOSE(10)
            END IF
            
            IF(thermoStat==idealgas_law)  THEN
                OPEN(10,FILE='CHK_StateEquation_UNDIM_T.dat')
                WRITE(10,'(A)')      '#1H     2T    3D    4M     5K      6CP      7BETA'
                WRITE(10,'(A,I3.1)') '# ', 256
                DO I=1,256
                    T_temp = T_WAL_GV(1,1)+(T_WAL_GV(1,2)-T_WAL_GV(1,1))/255.0_WP*(I-1)
                    P_temp = 1.0_wp
                    !call idealgas_TH(T_temp,H_temp)
                    !call idealgas_TD(T_temp,D_temp,P_temp)
                    !call idealgas_TM(T_temp,M_temp)
                    !call idealgas_TK(T_temp,K_temp)
                    !call idealgas_TCp(T_temp,Cp_temp)
                    !call idealgas_TB(T_temp,B_temp)
    
                    WRITE(10,'(7ES13.5)') H_temp, T_temp, D_temp, M_temp, K_temp, Cp_temp, B_temp
                END DO
                CLOSE(10)
            END IF
                
            call CHKHDL    ('   Wall States (dimensionaless) (y=-1) & y=(+1)',MYID)
            call CHK2RLHDL ('       Wall T=          ',MYID,T_WAL_GV(1,1),T_WAL_GV(1,2))    
            call CHK2RLHDL ('       Wall H=          ',MYID,H_WAL_GV(1,1),H_WAL_GV(1,2))    
            call CHK2RLHDL ('       Wall D=          ',MYID,D_WAL_GV(1,1),D_WAL_GV(1,2))    
            call CHK2RLHDL ('       Wall M=          ',MYID,M_WAL_GV(1,1),M_WAL_GV(1,2))    
            call CHK2RLHDL ('       Wall K=          ',MYID,K_WAL_GV(1,1),K_WAL_GV(1,2))    
            
            !========write ===========
            call CHKHDL    ('   Reference States (dimensional)',MYID)
            call CHKRLHDL  ('       REF P0(Pa)=      ',MYID,P0)
            call CHKRLHDL  ('       REF RHOU20(Pa)=  ',MYID,RHOU20)
            call CHKRLHDL  ('       REF L0(m)=       ',MYID,L0)
            call CHKRLHDL  ('       REF U0(m/s)=     ',MYID,U0)
            call CHKRLHDL  ('       REF Mdot(Kg/m2s)=',MYID,U0*D0)
            call CHKRLHDL  ('       REF Mdot(Kg/m2s)=',MYID,massflux)
            call CHKRLHDL  ('       REF T0(K)=       ',MYID,T0)
            call CHKRLHDL  ('       REF H0(J)=       ',MYID,H0)
            call CHKRLHDL  ('       REF D0(Kg/m3)=   ',MYID,D0)
            call CHKRLHDL  ('       REF Cp(J/Kg/K)=  ',MYID,CP0)
            call CHKRLHDL  ('       REF K0(W/m-K)=   ',MYID,K0)
            call CHKRLHDL  ('       REF M0(Pa-s)=    ',MYID,M0)
            call CHKRLHDL  ('       REF B0(1/K)=     ',MYID,B0)

            IF(thermoStat==idealgas_law) THEN
            call CHKRLHDL  ('       REF Cv(J/Kg/K)=  ',MYID,Cv0)
            call CHKRLHDL  ('       REF Gamma=       ',MYID,Gamma0)
            END IF
            

            
            IF( BCWALLHEAT(itopwall)==isoFluxWall .and. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
                call CHKRLHDL  ('       REF WALHFLX(W/m2)=',MYID,WHEAT0_DIM(1))
                call CHKRLHDL  ('       REF WALHFLX(W/m2)=',MYID,WHEAT0_DIM(2))
            END IF
            IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                call CHKRLHDL  ('       REF T1(K)=       ',MYID,WHEAT0_DIM(1))
                call CHKRLHDL  ('       REF T2(K)=       ',MYID,WHEAT0_DIM(2))
            END IF
            
            
            call CHKHDL    ('   Reference States (undimensional)',MYID)
            call CHKRLHDL  ('       REF CP0=         ',MYID,CP0)
            
            call CHKRLHDL  ('       REF Re=          ',MYID,REN)
            call CHKRLHDL  ('       REF PRT0=        ',MYID,PRT0)
            call CHKRLHDL  ('       REF 1/Fr2=        ',MYID,F_A)
                    
            call CHKRLHDL  ('       Wall HEAT (BOT) FROM X/L0=     ',MYID,XND_io(NIN1+1))
            call CHKRLHDL  ('       Wall HEAT (TOP) FROM X/L0=     ',MYID,XND_io(NIN2+1))
            
            IF( BCWALLHEAT(itopwall)==isoFluxWall .and. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
                call CHKRLHDL  ('       REF Gr0(Bottom Wall)=                 ',MYID,Gr0(1))
                call CHKRLHDL  ('       REF Gr0(Top Wall)=                    ',MYID,Gr0(2))
                
                call CHKRLHDL  ('       REF Gr0/RE^2.7(Bottom Wall)=          ',MYID,BO0(1,1))
                call CHKRLHDL  ('       REF Gr0/RE^2.7(Top Wall)=             ',MYID,BO0(1,2))
                
                call CHKRLHDL  ('       REF Gr0/RE^3.425/Pr^0.8(Bottom Wall)= ',MYID,BO0(2,1))
                call CHKRLHDL  ('       REF Gr0/RE^3.425/Pr^0.8(Top Wall)=    ',MYID,BO0(2,2))
                
                call CHKRLHDL  ('       REF Q+(WHFLUX0) (Bottom Wall)=        ',MYID,WHEAT0_UNDIM(1))
                call CHKRLHDL  ('       REF Q+(WHFLUX0) (Top Wall)=           ',MYID,WHEAT0_UNDIM(2))
            END IF
            
            IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                call CHKRLHDL  ('       REF Gr0(guess value for Bottom Wall)= ',MYID,Gr0(1))
                call CHKRLHDL  ('       REF Gr0(guess value for Top Wall)=    ',MYID,Gr0(2))
                
                call CHKRLHDL  ('       REF Gr0/RE^2.7(Bottom Wall)=          ',MYID,BO0(1,1))
                call CHKRLHDL  ('       REF Gr0/RE^2.7(Top Wall)=             ',MYID,BO0(1,2))
                
                call CHKRLHDL  ('       REF Gr0/RE^3.425/Pr^0.8(Bottom Wall)= ',MYID,BO0(2,1))
                call CHKRLHDL  ('       REF Gr0/RE^3.425/Pr^0.8(Top Wall)=    ',MYID,BO0(2,2))
                
                call CHKRLHDL  ('       REF T+(WALLTEMP) (Bottom Wall)=       ',MYID,WHEAT0_UNDIM(1))
                call CHKRLHDL  ('       REF T+(WALLTEMP) (Top Wall)=          ',MYID,WHEAT0_UNDIM(2))
            END IF
        
        END IF
        

        RETURN
    END SUBROUTINE 

!**********************************************************************
    SUBROUTINE NIST_TABLE_SCALING
!>    THE FORMAT OF NIST TABLE USED HERE SHUOLD BE
!>            # P H  T  \RHO \MU \KAPPA CP
!>            NUMBER OF LINES
!>            P  H  T  D  M  K CP
!>            ................
!>
!>    ALL H VALUSE ARE SORTED FROM SMALL TO BIG NUMBERS. 

!>    As H has been sorted, the Newton Bisection method can be used for searching.
!>    If the H data has not been sorted, binary tree method is recommended.      
      
        use thermal_info
        use mesh_info
        use flow_info
        USE INIT_INFO
        IMPLICIT NONE
      
        INTEGER(4)     :: IOS=0
        INTEGER(4)     :: FFLG=15           !file I/O id
        CHARACTER(64) :: STMP
        REAL(WP)        :: RTMP
        
        REAL(WP)    :: SFE, SFS
        INTEGER(4)  :: IE, IS,  NIN1, NIN2, I
        
        INTEGER(4):: JMAX=0, JMIN=0           ! how many peaks
        !=====max. 10 peaks max. and 10 peaks of mini.
        REAL(WP)  :: AMAXL(10)=0.0_WP, AMINL(10)=0.0_WP ! peak values
        INTEGER(4):: NMAXL(10)=0, NMINL(10)=0 ! peaks location
        
        
        IF(MYID.NE.0) RETURN
      
        !================READ IN NIST TABLE (DIMENSIONAL)=====================
        OPEN(FFLG,FILE=TRIM(NISTFLNM), status='old', iostat=ios)
        if(ios .NE. 0)  &
        call ERRHDL(' File '//TRIM(NISTFLNM)//' cannot be found.',MYID)
      
        READ(FFLG,*) STMP
        READ(FFLG,*) N_NIST
      
        ALLOCATE( NIST_H (N_NIST) ) ;  NIST_H  = 0.0_wp
        ALLOCATE( NIST_T (N_NIST) ) ;  NIST_T  = 0.0_wp
        ALLOCATE( NIST_D (N_NIST) ) ;  NIST_D  = 0.0_wp
        ALLOCATE( NIST_M (N_NIST) ) ;  NIST_M  = 0.0_wp
        ALLOCATE( NIST_K (N_NIST) ) ;  NIST_K  = 0.0_wp
        ALLOCATE( NIST_B (N_NIST) ) ;  NIST_B  = 0.0_wp
        ALLOCATE( NIST_CP(N_NIST) ) ;  NIST_CP = 0.0_wp
        MEMPC_byte = MEMPC_byte + N_NIST*7*8
      
        DO I=1,N_NIST
            READ(FFLG,*) RTMP, &
                        NIST_H(I), NIST_T(I), NIST_D(I), &
                        NIST_M(I), NIST_K(I), NIST_CP(I), NIST_B(I)
        END DO
        
        
        !==============reference statment====(dimensional)======================
        CALL BISECTION_NIST_T(T0,IE,IS)
        SFS = ( NIST_T(IE) - T0 )/( NIST_T(IE) - NIST_T(IS) )
        SFE = 1.0_WP - SFS
      
        H0  = SFE*NIST_H(IE)  + SFS*NIST_H(IS)
        D0  = SFE*NIST_D(IE)  + SFS*NIST_D(IS)
        K0  = SFE*NIST_K(IE)  + SFS*NIST_K(IS)
        M0  = SFE*NIST_M(IE)  + SFS*NIST_M(IS) 
        CP0 = SFE*NIST_CP(IE) + SFS*NIST_CP(IS) ! updated
        B0  = SFE*NIST_B(IE)  + SFS*NIST_B(IS)
        
        !================local peaks===========================================
        CALL FIND_LOCAL_MAXIMA(N_NIST,  NIST_K, JMAX, NMAXL, AMAXL, JMIN,  NMINL, AMINL)
        DO I = 1, JMAX
            call CHK2RLHDL  ('       Local Max. of K   = ', MYID, AMAXL(I)     ,AMAXL(I)/K0 )
            call CHK2RLHDL  ('             Locate at T = ', MYID, NIST_T(NMAXL(I)),NIST_T(NMAXL(I))/T0)
        END DO
        DO I = 1, JMIN
            call CHK2RLHDL  ('       Local Min. of K   = ', MYID, AMINL(I)     ,AMINL(I)/K0 )
            call CHK2RLHDL  ('             Locate at T = ', MYID, NIST_T(NMINL(I)),NIST_T(NMINL(I))/T0)
        END DO
        
        !===============scale the NIST table based on the reference state==(undimensionalization)===============
        DO I=1,N_NIST
            NIST_H(I)  = (NIST_H(I)-H0)/T0/CP0
            NIST_T(I)  =  NIST_T(I)/T0
            NIST_D(I)  =  NIST_D(I)/D0
            NIST_M(I)  =  NIST_M(I)/M0
            NIST_K(I)  =  NIST_K(I)/K0
            NIST_CP(I) =  NIST_CP(I)/CP0
            NIST_B(I)  =  NIST_B(I)/B0
        END DO
      
        call NIST_SPLINE_H
        call NIST_SPLINE_T
        
        
        RETURN
    END SUBROUTINE
    
!===============powerlaw for perfect gas====================    
    SUBROUTINE idealgas_HT(htmp,Ttmp)
        USE WPRECISION
        USE THERMAL_INFO
        IMPLICIT NONE
        REAL(WP),INTENT(IN)  :: Htmp
        REAL(WP),INTENT(OUT) :: Ttmp
        logical,external    :: isnan1


        !Ttmp = (PL_HT*(htmp+H0/CP0/T0))**(1.0_wp/PL_HT)
        Ttmp = (PL_HT*htmp+1.0_wp)**(1.0_wp/PL_HT)
        
        if(isnan1(Ttmp)) THEN
            write(*,*) 'idealgas_HT', htmp, ttmp
            STOP
        END IF
        if(isnan1(htmp)) THEN
            write(*,*) 'idealgas_HT',htmp, ttmp
            STOP
        END IF

        RETURN
    END SUBROUTINE
    
    SUBROUTINE idealgas_TH(Ttmp,Htmp)
        USE WPRECISION
        use thermal_info
        IMPLICIT NONE
        REAL(WP),INTENT(IN)  :: Ttmp
        REAL(WP),INTENT(OUT) :: Htmp
        
        !Htmp = (Ttmp**PL_HT)/PL_HT-H0/CP0/T0
        Htmp = (Ttmp**PL_HT - 1.0_wp)/PL_HT
        RETURN
    END SUBROUTINE
    
    SUBROUTINE idealgas_TD(Ttmp,dtmp, Ptmp)
        USE WPRECISION
        USE THERMAL_INFO
        IMPLICIT NONE
        REAL(WP),INTENT(IN)  :: Ttmp
        REAL(WP),INTENT(IN)  :: Ptmp
        REAL(WP),INTENT(OUT) :: dtmp
        REAL(WP) :: P_P0, Ptemp
        
        !Ptemp = Ptmp
        !Ptemp = 1.0_wp
        !P_P0 = (P0/D0/U0/U0 + Ptemp)/(P0/D0/U0/U0 + 1.0_wp)
        !WRITE(*,*) 'P/P0=', P_P0
        !dtmp = 1.0_WP/Ttmp
        !dtmp = 1.0_WP/Ttmp *  P_P0
        
        dtmp = 1.0_WP/Ttmp
        
        
        RETURN
    END SUBROUTINE
         
    SUBROUTINE idealgas_TM(Ttmp,Mtmp)  ! OK
        USE WPRECISION
        USE THERMAL_INFO
        IMPLICIT NONE
        REAL(WP),INTENT(IN)  :: Ttmp
        REAL(WP),INTENT(OUT) :: Mtmp
        
        REAL(WP) :: Ttmpd, Mtmpd
        
        IF(fstatetype==powerlawflg) THEN
            Mtmp = Ttmp**PL_MT
        ELSE IF(fstatetype==sutherlandflg) THEN
            Ttmpd = Ttmp * T0
            CALL Sutherland_Viscosity_dim(Ttmpd,Mtmpd)
            Mtmp = Mtmpd/M0
        ELSE
            Mtmp = Ttmp**PL_MT
        END IF
        
        RETURN
    END SUBROUTINE
    
    SUBROUTINE idealgas_TK(Ttmp,Ktmp)   ! OK
        USE WPRECISION
        USE THERMAL_INFO
        IMPLICIT NONE
        REAL(WP),INTENT(IN)  :: Ttmp
        REAL(WP),INTENT(OUT) :: Ktmp
    
        REAL(WP) :: Ttmpd, Ktmpd
    
        IF(fstatetype==powerlawflg) THEN
            Ktmp = Ttmp**PL_KT
        ELSE IF(fstatetype==sutherlandflg) THEN
            Ttmpd = Ttmp * T0
            CALL Sutherland_Conductivity_dim(Ttmpd,Ktmpd)
            Ktmp = Ktmpd/K0
        ELSE
            Ktmp = Ttmp**PL_KT
        END IF
        
        RETURN
    END SUBROUTINE
    
    SUBROUTINE idealgas_TCp(Ttmp,Cptmp) ! OK
        USE WPRECISION
        USE THERMAL_INFO
        IMPLICIT NONE
        REAL(WP),INTENT(IN)  :: Ttmp
        REAL(WP),INTENT(OUT) :: Cptmp
        
        Cptmp = Ttmp**PL_CpT
        RETURN
    END SUBROUTINE
    
    SUBROUTINE idealgas_TB(Ttmp,Btmp)  ! OK
        USE WPRECISION
        USE THERMAL_INFO
        IMPLICIT NONE
        REAL(WP),INTENT(IN)  :: Ttmp
        REAL(WP),INTENT(OUT) :: Btmp
        
        Btmp = 1.0_wp/Ttmp
        RETURN
    END SUBROUTINE
!===============powerlaw for perfect gas====================     


!=========Sutherland Law for dynamic viscousity================
    SUBROUTINE Sutherland_Viscosity_dim(TT,MM)
        USE WPRECISION
        USE Flow_State_Constant
        IMPLICIT NONE
        REAL(WP),INTENT(IN) :: TT
        REAL(WP),INTENT(OUT):: MM
        
        ! Ref : http://www.cfd-online.com/Wiki/Sutherland's_law
        
        MM = SL_MU0 * ( (TT/SL_T0)**(3.0_wp/2.0_wp) ) *(SL_T0+SL_S)/(TT+SL_S)
        !WRITE(*,*) MM
        !or
        
        !MM = SL_C1 * TT**(3.0_wp/2.0_wp) / (TT+SL_S)
        !WRITE(*,*) MM
         
        RETURN
    END SUBROUTINE 

!=========Sutherland Law for conductivity================    
    SUBROUTINE Sutherland_Conductivity_dim(TT,KK)
        USE WPRECISION
        USE Flow_State_Constant
        IMPLICIT NONE
        REAL(WP),INTENT(IN) :: TT
        REAL(WP),INTENT(OUT):: KK
        
        KK = SL_K0 * ( (TT/SL_T0)**(3.0_wp/2.0_wp) ) *(SL_T0+SL_S)/(TT+SL_S)
        !WRITE(*,*) KK
        
        RETURN
    END SUBROUTINE
    !*****************************************************************************    
    SUBROUTINE WALL_ISOTHERMAL_GIVEN(N)
        use flow_info
        USE MESH_INFO
        USE THERMAL_INFO
        use init_info
        IMPLICIT NONE
        
        INTEGER(4) :: I,N  !,I1, I2, IE, IS
        REAL(WP)   :: T_TEMP  !, SFS, SFE
        REAL(WP)   :: H_temp, M_temp, K_temp, D_temp, Cp_temp, P_temp
        
        !CALL DEBUG_WRT_LOCAL(TEMPERATURE,0,N2DO(MYID)+1)  !test
        
        IF(BCWALLHEAT(N).ne.isoThermalWall ) RETURN
        DO I=NCL1S,NCL1E
            T_TEMP =  T_WAL_GV(I,N)
            
            IF(thermoStat==search_table) Then
                call NIST_SLEVAL_TH(T_temp,H_temp)
                call NIST_SLEVAL_TD(T_temp,D_temp)
                call NIST_SLEVAL_TM(T_temp,M_temp)
                call NIST_SLEVAL_TK(T_temp,K_temp)
                call NIST_SLEVAL_TCp(T_temp,Cp_temp)

            ELSE IF(thermoStat==idealgas_law) Then
                P_temp = 1.0_wp
                !call idealgas_TH(T_temp,H_temp)
                !call idealgas_TD(T_temp,D_temp,P_temp)
                !call idealgas_TM(T_temp,M_temp)
                !call idealgas_TK(T_temp,K_temp)
                !call idealgas_TCp(T_temp,Cp_temp)
            ELSE
            END IF
            
            H_WAL_GV(I,N) = H_TEMP
            D_WAL_GV(I,N) = D_TEMP
            M_WAL_GV(I,N) = M_TEMP
            K_WAL_GV(I,N) = K_TEMP
            Cp_WAL_GV(I,N)= Cp_temp
        END DO
        
        
        RETURN
    END SUBROUTINE
    
    !*****************************************************************************    
    SUBROUTINE BCAST_THERMAL_INIT
        use thermal_info
        use mesh_info
        use flow_info
        USE INIT_INFO
        IMPLICIT NONE
        
        CALL MPI_BCAST( N_NIST, 1, MPI_INTEGER4,      0, ICOMM, IERROR )
        
        CALL MPI_BCAST( H0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( D0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( K0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( M0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( B0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( Cp0, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( PRT0,  1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( CTHECD,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( U0,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( RHOU20,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( G_A,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( F_A,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( IBuoF, 3, MPI_INTEGER4,         0, ICOMM, IERROR )
        
        CALL MPI_BCAST( RHOHmax,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( T4MaxRhoh, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( Cpmax,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( T4CpMax, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        CALL MPI_BCAST( WHEAT0_UNDIM, 2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        !CALL MPI_BCAST( WHFLUX_UND,   2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( Gr0,          2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( BO0,          4, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        
        IF(thermoStat==idealgas_law) THEN
            CALL MPI_BCAST( R0,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( Cv0,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( gamma0,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( Cvhat0,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( PL_MT, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( PL_CpT,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( PL_HT, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( PL_KT, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
        END IF
        
        IF(MYID.NE.0) THEN
        
            IF( BCWALLHEAT(itopwall)==isoFluxWall .or. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
                ALLOCATE (WALLFLUX(NCL1S:NCL1E,2)  )
            END IF
            
            IF( BCWALLHEAT(itopwall)==isoThermalWall .or. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
                ALLOCATE (T_WAL_GV(NCL1S:NCL1E,2)  )
                ALLOCATE (H_WAL_GV(NCL1S:NCL1E,2)  )
                ALLOCATE (D_WAL_GV(NCL1S:NCL1E,2)  )
                ALLOCATE (M_WAL_GV(NCL1S:NCL1E,2)  )
                ALLOCATE (K_WAL_GV(NCL1S:NCL1E,2)  )
                ALLOCATE (Cp_WAL_GV(NCL1S:NCL1E,2)  )
            END IF
            
            IF(thermoStat==search_table) THEN
                ALLOCATE( NIST_H (N_NIST) ) ;  NIST_H  = 0.0_wp
                ALLOCATE( NIST_T (N_NIST) ) ;  NIST_T  = 0.0_wp
                ALLOCATE( NIST_D (N_NIST) ) ;  NIST_D  = 0.0_wp
                ALLOCATE( NIST_M (N_NIST) ) ;  NIST_M  = 0.0_wp
                ALLOCATE( NIST_K (N_NIST) ) ;  NIST_K  = 0.0_wp
                ALLOCATE( NIST_B (N_NIST) ) ;  NIST_B  = 0.0_wp
                ALLOCATE( NIST_CP(N_NIST) ) ;  NIST_CP = 0.0_wp
                
                ALLOCATE(NIST_HT_B(N_NIST)); NIST_HT_B=0.0_WP
                ALLOCATE(NIST_HT_C(N_NIST)); NIST_HT_C=0.0_WP
                ALLOCATE(NIST_HT_D(N_NIST)); NIST_HT_D=0.0_WP
                
                ALLOCATE(NIST_HD_B(N_NIST)); NIST_HD_B=0.0_WP
                ALLOCATE(NIST_HD_C(N_NIST)); NIST_HD_C=0.0_WP
                ALLOCATE(NIST_HD_D(N_NIST)); NIST_HD_D=0.0_WP
                
                ALLOCATE(NIST_HM_B(N_NIST)); NIST_HM_B=0.0_WP
                ALLOCATE(NIST_HM_C(N_NIST)); NIST_HM_C=0.0_WP
                ALLOCATE(NIST_HM_D(N_NIST)); NIST_HM_D=0.0_WP
                
                ALLOCATE(NIST_HK_B(N_NIST)); NIST_HK_B=0.0_WP
                ALLOCATE(NIST_HK_C(N_NIST)); NIST_HK_C=0.0_WP
                ALLOCATE(NIST_HK_D(N_NIST)); NIST_HK_D=0.0_WP
                
                ALLOCATE(NIST_HCP_B(N_NIST)); NIST_HCP_B=0.0_WP
                ALLOCATE(NIST_HCP_C(N_NIST)); NIST_HCP_C=0.0_WP
                ALLOCATE(NIST_HCP_D(N_NIST)); NIST_HCP_D=0.0_WP
                
                ALLOCATE(NIST_HB_B(N_NIST)); NIST_HB_B=0.0_WP
                ALLOCATE(NIST_HB_C(N_NIST)); NIST_HB_C=0.0_WP
                ALLOCATE(NIST_HB_D(N_NIST)); NIST_HB_D=0.0_WP
                
                ALLOCATE(NIST_TH_B(N_NIST)); NIST_TH_B=0.0_WP
                ALLOCATE(NIST_TH_C(N_NIST)); NIST_TH_C=0.0_WP
                ALLOCATE(NIST_TH_D(N_NIST)); NIST_TH_D=0.0_WP
                
                ALLOCATE(NIST_TD_B(N_NIST)); NIST_TD_B=0.0_WP
                ALLOCATE(NIST_TD_C(N_NIST)); NIST_TD_C=0.0_WP
                ALLOCATE(NIST_TD_D(N_NIST)); NIST_TD_D=0.0_WP
                
                ALLOCATE(NIST_TM_B(N_NIST)); NIST_TM_B=0.0_WP
                ALLOCATE(NIST_TM_C(N_NIST)); NIST_TM_C=0.0_WP
                ALLOCATE(NIST_TM_D(N_NIST)); NIST_TM_D=0.0_WP
                
                ALLOCATE(NIST_TK_B(N_NIST)); NIST_TK_B=0.0_WP
                ALLOCATE(NIST_TK_C(N_NIST)); NIST_TK_C=0.0_WP
                ALLOCATE(NIST_TK_D(N_NIST)); NIST_TK_D=0.0_WP
                
                ALLOCATE(NIST_TCP_B(N_NIST)); NIST_TCP_B=0.0_WP
                ALLOCATE(NIST_TCP_C(N_NIST)); NIST_TCP_C=0.0_WP
                ALLOCATE(NIST_TCP_D(N_NIST)); NIST_TCP_D=0.0_WP
                
                ALLOCATE(NIST_TB_B(N_NIST)); NIST_TB_B=0.0_WP
                ALLOCATE(NIST_TB_C(N_NIST)); NIST_TB_C=0.0_WP
                ALLOCATE(NIST_TB_D(N_NIST)); NIST_TB_D=0.0_WP
            END IF
        
        END IF
        
        IF( BCWALLHEAT(itopwall)==isoFluxWall .or. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
            CALL MPI_BCAST( WALLFLUX,     2*(NCL1E-NCL1S+1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        END IF
        
        IF( BCWALLHEAT(itopwall)==isoThermalWall .or. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
            CALL MPI_BCAST( T_WAL_GV,     2*(NCL1E-NCL1S+1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( H_WAL_GV,     2*(NCL1E-NCL1S+1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( D_WAL_GV,     2*(NCL1E-NCL1S+1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( M_WAL_GV,     2*(NCL1E-NCL1S+1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( K_WAL_GV,     2*(NCL1E-NCL1S+1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST(Cp_WAL_GV,     2*(NCL1E-NCL1S+1), MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        END IF
        
        IF(thermoStat==search_table) THEN
            CALL MPI_BCAST( NIST_H, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_T, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_M, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_K, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_Cp,N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_HT_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HT_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HT_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_HD_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HD_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HD_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_HM_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HM_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HM_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_HK_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HK_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HK_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_HB_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HB_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HB_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_HCP_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HCP_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_HCP_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_TH_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TH_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TH_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_TD_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TD_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TD_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_TM_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TM_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TM_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_TK_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TK_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TK_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_TB_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TB_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TB_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            
            CALL MPI_BCAST( NIST_TCP_B, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TCP_C, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
            CALL MPI_BCAST( NIST_TCP_D, N_NIST, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        END IF
        
        
        RETURN
    END SUBROUTINE
    
!*********************************************************************************
    SUBROUTINE BISECTION_NIST_H(ENTH,IE,IS)
        USE thermal_info
        IMPLICIT NONE
      
        REAL(WP),INTENT(IN)      :: ENTH
        INTEGER(4),INTENT(OUT) :: IE,IS
        REAL(WP) :: VS, VM
        INTEGER(4) :: IM
      
!      IF( ( ENTH.LT.NIST_H(1) ).OR. ( ENTH .GT. NIST_H(N_NIST) ) ) THEN
!        WRITE(*,'(A,2ES13.5)') 'RANGE OF NIST_H =', NIST_H(1), NIST_H(N_NIST)
!        WRITE(*,'(A,1ES13.5)') 'CURRENT GIVEN ENTHALPY= ', ENTH
!        CALL ERRHDL(' ENTHOPHY IS OUT OF GIVEN NIST TABLE RANGE.',MYID)
!      END IF
      
        IS = 1
        IE = N_NIST
        DO WHILE ( (IE-IS).GT. 1 )
            IM = IS + (IE-IS)/2
            VS = NIST_H(IS) - ENTH
            VM = NIST_H(IM) - ENTH
            IF( VS*VM .GT. 0.0_WP) THEN
                IS = IM
            ELSE
                IE = IM
            END IF 
        END DO
      
        RETURN
    END SUBROUTINE 
      
!*********************************************************************************      
    SUBROUTINE BISECTION_NIST_T(TEMP,IE,IS)
!>    1) TABLE USED HERE IS DIMENSIONAL.
!>    2) AS H IS SORTED FROM SMALL TO LARGE NUMBERS, AND 
!>          T IS A MONOTONE FUNCTION OF H, THUS
!>          T IS ALSO FROM SMALL TO LARGE IN THE TABLE.

        USE thermal_info
        IMPLICIT NONE
      
        REAL(WP),INTENT(IN)      :: TEMP
        INTEGER(4),INTENT(OUT) :: IE,IS
        REAL(WP) :: VS, VM
        INTEGER(4) :: IM
      
        !IF( ( TEMP.LT.NIST_T(1) ).OR. ( TEMP .GT. NIST_T(N_NIST) ) ) &
        !CALL ERRHDL(' REFERENCE TEMPERATURE IS OUT OF GIVEN NIST TABLE RANGE.',MYID)
      
        IS = 1
        IE = N_NIST
        DO WHILE ( (IE-IS).GT. 1 )
            IM = IS + (IE-IS)/2
            VS = NIST_T(IS) - TEMP
            VM = NIST_T(IM) - TEMP
            IF( VS*VM .GT. 0.0_WP) THEN
                IS = IM
            ELSE
                IE = IM
            END IF 
        END DO
      
        RETURN
    END SUBROUTINE 
      
      
    !===================================================================
    SUBROUTINE NIST_SPLINE_H
        USE thermal_info
        IMPLICIT NONE
        INTEGER(4) :: I
        real(wp)   :: H_temp
        real(wp)   :: T_temp, D_temp, M_temp, K_temp, Cp_temp, B_temp, dDdH_temp
        real(wp)   :: D0_temp, H0_temp, dDdH0_temp
        
        IF(MYID.NE.0) RETURN
        
        ALLOCATE(NIST_HT_B(N_NIST)); NIST_HT_B=0.0_WP
        ALLOCATE(NIST_HT_C(N_NIST)); NIST_HT_C=0.0_WP
        ALLOCATE(NIST_HT_D(N_NIST)); NIST_HT_D=0.0_WP
        
        ALLOCATE(NIST_HD_B(N_NIST)); NIST_HD_B=0.0_WP
        ALLOCATE(NIST_HD_C(N_NIST)); NIST_HD_C=0.0_WP
        ALLOCATE(NIST_HD_D(N_NIST)); NIST_HD_D=0.0_WP
        
        ALLOCATE(NIST_HM_B(N_NIST)); NIST_HM_B=0.0_WP
        ALLOCATE(NIST_HM_C(N_NIST)); NIST_HM_C=0.0_WP
        ALLOCATE(NIST_HM_D(N_NIST)); NIST_HM_D=0.0_WP
        
        ALLOCATE(NIST_HK_B(N_NIST)); NIST_HK_B=0.0_WP
        ALLOCATE(NIST_HK_C(N_NIST)); NIST_HK_C=0.0_WP
        ALLOCATE(NIST_HK_D(N_NIST)); NIST_HK_D=0.0_WP
        
        ALLOCATE(NIST_HCP_B(N_NIST)); NIST_HCP_B=0.0_WP
        ALLOCATE(NIST_HCP_C(N_NIST)); NIST_HCP_C=0.0_WP
        ALLOCATE(NIST_HCP_D(N_NIST)); NIST_HCP_D=0.0_WP
        
        ALLOCATE(NIST_HB_B(N_NIST)); NIST_HB_B=0.0_WP
        ALLOCATE(NIST_HB_C(N_NIST)); NIST_HB_C=0.0_WP
        ALLOCATE(NIST_HB_D(N_NIST)); NIST_HB_D=0.0_WP

        CALL CUBIC_SPLINE (N_NIST, NIST_H, NIST_T, NIST_HT_B, NIST_HT_C, NIST_HT_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_H, NIST_D, NIST_HD_B, NIST_HD_C, NIST_HD_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_H, NIST_M, NIST_HM_B, NIST_HM_C, NIST_HM_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_H, NIST_K, NIST_HK_B, NIST_HK_C, NIST_HK_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_H, NIST_Cp,NIST_HCp_B,NIST_HCp_C,NIST_HCp_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_H, NIST_B, NIST_HB_B, NIST_HB_C, NIST_HB_D)
    
        
        IF(MYID==0) THEN
            OPEN(10,FILE='CHK_NIST_TABLE_UNDIM_SPLINE_H.dat')
            WRITE(10,'(A)') '#P(Mpa)      H     T    D    M     K      CP      BETA  dDdH'
            WRITE(10,*) '# ', N_NIST*3+1
            
            RHOHmax =-1E-14_wp
            CpMax=-1E-14_wp
            DO I=2,N_NIST*3

                D0_temp = D_temp
                H0_temp = H_temp

                H_temp = NIST_H(1)+DBLE(I-1)*(NIST_H(N_NIST)-NIST_H(1))/DBLE(N_NIST*3)
                call NIST_SLEVAL_HT(H_temp,T_temp)
                call NIST_SLEVAL_HD(H_temp,D_temp)
                call NIST_SLEVAL_HM(H_temp,M_temp)
                call NIST_SLEVAL_HK(H_temp,K_temp)
                call NIST_SLEVAL_HCp(H_temp,Cp_temp)
                call NIST_SLEVAL_HB(H_temp,B_temp)

                IF(Cp_Temp.GT.CpMax) THEN
                    CpMax   = Cp_temp
                    T4CpMax = T_temp 
                END IF

                IF((D_temp*H_temp).gt.RHOHmax) THEN
                    RHOHmax = D_temp*H_temp
                    T4MaxRhoh = T_temp
                END IF

                dDdH0_temp = (D_temp-D0_temp)/(H_temp-H0_temp)
                
                WRITE(10,'(9ES13.5)') P0, H_temp, T_temp, D_temp, M_temp, K_temp, Cp_temp, B_temp, dDdH0_temp
            END DO
            WRITE(10,*) '# Max RHOH:               ', RHOHmax
            WRITE(10,*) '# Corresponding Tempeture:', T4MaxRhoh
            CLOSE(10)
            
            call CHKRLHDL  ('       Max RHOH (undim):                = ',MYID,RHOHmax)
            call CHKRLHDL  ('           At Temperature(undim):       = ',MYID,T4MaxRhoh)
            call CHKRLHDL  ('           At Temperature(dim):         = ',MYID,T4MaxRhoh*T0)
            
            call CHKRLHDL  ('       Max Cp (undim):                  = ',MYID,CpMax)
            call CHKRLHDL  ('           (Tpc)At Temperature(undim):  = ',MYID,T4CpMax)
            call CHKRLHDL  ('           (Tpc)At Temperature(dim):    = ',MYID,T4CpMax*T0)
            
            
            
        END IF
        
      
    RETURN
    END SUBROUTINE
    
    
 
    !==============================================================================
    SUBROUTINE NIST_SLEVAL_HT(H,EVAL)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: H,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(H <= NIST_H(1)) then
          !EVAL = NIST_T(1)
          EVAL = ( NIST_T(1)-NIST_T(2) )/( NIST_H(1)-NIST_H(2) ) *(H-NIST_H(1))+NIST_T(1)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its minimum value in the given NIST table.', NIST_H(1), H 
          STOP
        end if
        if(H >= NIST_H(n)) then
          !EVAL = NIST_T(n)
          EVAL = ( NIST_T(n)-NIST_T(n-1) )/( NIST_H(n)-NIST_H(n-1) ) *(H-NIST_H(n))+NIST_T(n)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its maximum value in the given NIST table.', NIST_H(n), H 
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(h < NIST_H(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = h - NIST_H(i)
        EVAL = NIST_T(i) + dx*(NIST_HT_b(i) + dx*(NIST_HT_c(i) + dx*NIST_HT_d(i)))
        
        IF(DABS(H).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
        
    !==============================================================================
    SUBROUTINE NIST_SLEVAL_HD(H,EVAL)
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: H,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(H <= NIST_H(1)) then
          !EVAL = NIST_T(1)
          EVAL = ( NIST_D(1)-NIST_D(2) )/( NIST_H(1)-NIST_H(2) ) *(H-NIST_H(1))+NIST_D(1)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its minimum value in the given NIST table.', NIST_H(1), H 
          STOP
        end if
        if(H >= NIST_H(n)) then
          !EVAL = NIST_T(n)
          EVAL = ( NIST_D(n)-NIST_D(n-1) )/( NIST_H(n)-NIST_H(n-1) ) *(H-NIST_H(n))+NIST_D(n)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its maximum value in the given NIST table.', NIST_H(N), H 
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(h < NIST_H(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = h - NIST_H(i)
        EVAL = NIST_D(i) + dx*(NIST_HD_b(i) + dx*(NIST_HD_c(i) + dx*NIST_HD_d(i)))
        IF(DABS(H).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    
    !==============================================================================
    SUBROUTINE NIST_SLEVAL_HM(H,EVAL)
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: H,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(H <= NIST_H(1)) then
          !EVAL = NIST_T(1)
          EVAL = ( NIST_M(1)-NIST_M(2) )/( NIST_H(1)-NIST_H(2) ) *(H-NIST_H(1))+NIST_M(1)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its minimum value in the given NIST table.', NIST_H(1), H 
          STOP
        end if
        if(H >= NIST_H(n)) then
          !EVAL = NIST_T(n)
          EVAL = ( NIST_M(n)-NIST_M(n-1) )/( NIST_H(n)-NIST_H(n-1) ) *(H-NIST_H(n))+NIST_M(n)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its maximum value in the given NIST table.', NIST_H(N), H 
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(h < NIST_H(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = h - NIST_H(i)
        EVAL = NIST_M(i) + dx*(NIST_HM_b(i) + dx*(NIST_HM_c(i) + dx*NIST_HM_d(i)))
        IF(DABS(H).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    !==============================================================================
    SUBROUTINE NIST_SLEVAL_HK(H,EVAL)
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: H,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(H <= NIST_H(1)) then
          !EVAL = NIST_T(1)
          EVAL = ( NIST_K(1)-NIST_K(2) )/( NIST_H(1)-NIST_H(2) ) *(H-NIST_H(1))+NIST_K(1)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its minimum value in the given NIST table.', NIST_H(1), H 
          STOP
        end if
        if(H >= NIST_H(n)) then
          !EVAL = NIST_T(n)
          EVAL = ( NIST_K(n)-NIST_K(n-1) )/( NIST_H(n)-NIST_H(n-1) ) *(H-NIST_H(n))+NIST_K(n)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its maximum value in the given NIST table.', NIST_H(N), H 
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(h < NIST_H(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = h - NIST_H(i)
        EVAL = NIST_K(i) + dx*(NIST_HK_b(i) + dx*(NIST_HK_c(i) + dx*NIST_HK_d(i)))
        IF(DABS(H).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    !==============================================================================
    SUBROUTINE NIST_SLEVAL_HB(H,EVAL)
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: H,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(H <= NIST_H(1)) then
          !EVAL = NIST_T(1)
          EVAL = ( NIST_B(1)-NIST_B(2) )/( NIST_H(1)-NIST_H(2) ) *(H-NIST_H(1))+NIST_B(1)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its minimum value in the given NIST table.', NIST_H(1), H 
          STOP
        end if
        if(H >= NIST_H(n)) then
          !EVAL = NIST_T(n)
          EVAL = ( NIST_B(n)-NIST_B(n-1) )/( NIST_H(n)-NIST_H(n-1) ) *(H-NIST_H(n))+NIST_B(n)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its maximum value in the given NIST table.', NIST_H(N), H 
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(h < NIST_H(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = h - NIST_H(i)
        EVAL = NIST_B(i) + dx*(NIST_HB_b(i) + dx*(NIST_HB_c(i) + dx*NIST_HB_d(i)))
        IF(DABS(H).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    !==============================================================================
    SUBROUTINE NIST_SLEVAL_HCp(H,EVAL)
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: H,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(H <= NIST_H(1)) then
          !EVAL = NIST_T(1)
          EVAL = ( NIST_Cp(1)-NIST_Cp(2) )/( NIST_H(1)-NIST_H(2) ) *(H-NIST_H(1))+NIST_Cp(1)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its minimum value in the given NIST table.', NIST_H(1), H 
          STOP
        end if
        if(H >= NIST_H(n)) then
          !EVAL = NIST_T(n)
          EVAL = ( NIST_Cp(n)-NIST_Cp(n-1) )/( NIST_H(n)-NIST_H(n-1) ) *(H-NIST_H(n))+NIST_Cp(n)
          WRITE(*,'(A,2F12.5)') '# Warning: H reaches its maximum value in the given NIST table.', NIST_H(N), H 
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(h < NIST_H(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = h - NIST_H(i)
        EVAL = NIST_Cp(i) + dx*(NIST_HCp_b(i) + dx*(NIST_HCp_c(i) + dx*NIST_HCp_d(i)))
        IF(DABS(H).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    
!====================based on T=====================
!===================================================================
    SUBROUTINE NIST_SPLINE_T
        USE thermal_info
        IMPLICIT NONE
        INTEGER(4) :: I
        real(wp)   :: H_temp
        real(wp)   :: T_temp, D_temp, M_temp, K_temp, Cp_temp, B_temp
        
        IF(MYID.NE.0) RETURN
        
        
        ALLOCATE(NIST_TH_B(N_NIST)); NIST_TH_B=0.0_WP
        ALLOCATE(NIST_TH_C(N_NIST)); NIST_TH_C=0.0_WP
        ALLOCATE(NIST_TH_D(N_NIST)); NIST_TH_D=0.0_WP
        
        ALLOCATE(NIST_TD_B(N_NIST)); NIST_TD_B=0.0_WP
        ALLOCATE(NIST_TD_C(N_NIST)); NIST_TD_C=0.0_WP
        ALLOCATE(NIST_TD_D(N_NIST)); NIST_TD_D=0.0_WP
        
        ALLOCATE(NIST_TM_B(N_NIST)); NIST_TM_B=0.0_WP
        ALLOCATE(NIST_TM_C(N_NIST)); NIST_TM_C=0.0_WP
        ALLOCATE(NIST_TM_D(N_NIST)); NIST_TM_D=0.0_WP
        
        ALLOCATE(NIST_TK_B(N_NIST)); NIST_TK_B=0.0_WP
        ALLOCATE(NIST_TK_C(N_NIST)); NIST_TK_C=0.0_WP
        ALLOCATE(NIST_TK_D(N_NIST)); NIST_TK_D=0.0_WP
        
        ALLOCATE(NIST_TCP_B(N_NIST)); NIST_TCP_B=0.0_WP
        ALLOCATE(NIST_TCP_C(N_NIST)); NIST_TCP_C=0.0_WP
        ALLOCATE(NIST_TCP_D(N_NIST)); NIST_TCP_D=0.0_WP
        
        ALLOCATE(NIST_TB_B(N_NIST)); NIST_TB_B=0.0_WP
        ALLOCATE(NIST_TB_C(N_NIST)); NIST_TB_C=0.0_WP
        ALLOCATE(NIST_TB_D(N_NIST)); NIST_TB_D=0.0_WP

        CALL CUBIC_SPLINE (N_NIST, NIST_T, NIST_H, NIST_TH_B, NIST_TH_C, NIST_TH_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_T, NIST_D, NIST_TD_B, NIST_TD_C, NIST_TD_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_T, NIST_M, NIST_TM_B, NIST_TM_C, NIST_TM_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_T, NIST_K, NIST_TK_B, NIST_TK_C, NIST_TK_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_T, NIST_Cp,NIST_TCp_B,NIST_TCp_C,NIST_TCp_D)
        CALL CUBIC_SPLINE (N_NIST, NIST_T, NIST_B, NIST_TB_B, NIST_TB_C, NIST_TB_D)
    
        
        IF(MYID==0) THEN
            OPEN(10,FILE='CHK_NIST_TABLE_UNDIM_SPLINE_T.dat')
            WRITE(10,'(A)') '#P(Mpa)      H     T    D    M     K      CP      BETA'
            WRITE(10,*) '# ', N_NIST*3
            
            DO I=2,N_NIST*3
                T_temp = NIST_T(1)+DBLE(I-1)*(NIST_T(N_NIST)-NIST_T(1))/DBLE(N_NIST*3)
                call NIST_SLEVAL_TH(T_temp,H_temp)
                call NIST_SLEVAL_TD(T_temp,D_temp)
                call NIST_SLEVAL_TM(T_temp,M_temp)
                call NIST_SLEVAL_TK(T_temp,K_temp)
                call NIST_SLEVAL_TCp(T_temp,Cp_temp)
                call NIST_SLEVAL_TB(T_temp,B_temp)
                
                WRITE(10,'(8ES13.5)') P0, H_temp, T_temp, D_temp, M_temp, K_temp, Cp_temp, B_temp
            END DO
            CLOSE(10)
            
        END IF
        
      
    RETURN
    END SUBROUTINE
    
    
 
    !==============================================================================
    SUBROUTINE NIST_SLEVAL_TH(T,EVAL)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: T,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(T <= NIST_T(1)) then
          EVAL = ( NIST_H(1)-NIST_H(2) )/( NIST_T(1)-NIST_T(2) ) *(T-NIST_T(1))+NIST_H(1)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its minimum value in the given NIST table.', NIST_T(1), T
          STOP
        end if
        if(T >= NIST_T(n)) then
          EVAL = ( NIST_H(n)-NIST_H(n-1) )/( NIST_T(n)-NIST_T(n-1) ) *(T-NIST_T(n))+NIST_H(n)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its maximum value in the given NIST table.', NIST_T(N), T 
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(T < NIST_T(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = T - NIST_T(i)
        EVAL = NIST_H(i) + dx*(NIST_TH_b(i) + dx*(NIST_TH_c(i) + dx*NIST_TH_d(i)))
        IF(DABS(T-1.0_WP).LT.1.0E-10_WP) EVAL = 0.0_WP !test
        return
    end subroutine
    
    
    SUBROUTINE NIST_SLEVAL_TD(T,EVAL)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: T,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(T <= NIST_T(1)) then
          EVAL = ( NIST_D(1)-NIST_D(2) )/( NIST_T(1)-NIST_T(2) ) *(T-NIST_T(1))+NIST_D(1)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its minimum value in the given NIST table.', NIST_T(1), T
          STOP
        end if
        if(T >= NIST_T(n)) then
          EVAL = ( NIST_D(n)-NIST_D(n-1) )/( NIST_T(n)-NIST_T(n-1) ) *(T-NIST_T(n))+NIST_D(n)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its maximum value in the given NIST table.', NIST_T(N), T
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(T < NIST_T(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = T - NIST_T(i)
        EVAL = NIST_D(i) + dx*(NIST_TD_b(i) + dx*(NIST_TD_c(i) + dx*NIST_TD_d(i)))
        IF(DABS(T-1.0_WP).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    SUBROUTINE NIST_SLEVAL_TM(T,EVAL)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: T,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(T <= NIST_T(1)) then
          EVAL = ( NIST_M(1)-NIST_M(2) )/( NIST_T(1)-NIST_T(2) ) *(T-NIST_T(1))+NIST_M(1)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its minimum value in the given NIST table.', NIST_T(1), T
          STOP
        end if
        if(T >= NIST_T(n)) then
          EVAL = ( NIST_M(n)-NIST_M(n-1) )/( NIST_T(n)-NIST_T(n-1) ) *(T-NIST_T(n))+NIST_M(n)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its maximum value in the given NIST table.', NIST_T(N), T
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(T < NIST_T(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = T - NIST_T(i)
        EVAL = NIST_M(i) + dx*(NIST_TM_b(i) + dx*(NIST_TM_c(i) + dx*NIST_TM_d(i)))
        IF(DABS(T-1.0_WP).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    SUBROUTINE NIST_SLEVAL_TK(T,EVAL)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: T,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(T <= NIST_T(1)) then
          EVAL = ( NIST_K(1)-NIST_K(2) )/( NIST_T(1)-NIST_T(2) ) *(T-NIST_T(1))+NIST_K(1)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its minimum value in the given NIST table.', NIST_T(1), T
          STOP
        end if
        if(T >= NIST_T(n)) then
          EVAL = ( NIST_K(n)-NIST_K(n-1) )/( NIST_T(n)-NIST_T(n-1) ) *(T-NIST_T(n))+NIST_K(n)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its maximum value in the given NIST table.', NIST_T(N), T
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(T < NIST_T(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = T - NIST_T(i)
        EVAL = NIST_K(i) + dx*(NIST_TK_b(i) + dx*(NIST_TK_c(i) + dx*NIST_TK_d(i)))
        IF(DABS(T-1.0_WP).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    SUBROUTINE NIST_SLEVAL_TB(T,EVAL)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: T,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(T <= NIST_T(1)) then
          EVAL = ( NIST_B(1)-NIST_B(2) )/( NIST_T(1)-NIST_T(2) ) *(T-NIST_T(1))+NIST_B(1)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its minimum value in the given NIST table.', NIST_T(1), T
          STOP
        end if
        if(T >= NIST_T(n)) then
          EVAL = ( NIST_B(n)-NIST_B(n-1) )/( NIST_T(n)-NIST_T(n-1) ) *(T-NIST_T(n))+NIST_B(n)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its maximum value in the given NIST table.', NIST_T(N), T
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(T < NIST_T(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = T - NIST_T(i)
        EVAL = NIST_B(i) + dx*(NIST_TB_b(i) + dx*(NIST_TB_c(i) + dx*NIST_TB_d(i)))
        IF(DABS(T-1.0_WP).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
    SUBROUTINE NIST_SLEVAL_TCp(T,EVAL)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        use thermal_info
        IMPLICIT NONE
        REAL(WP)   :: T,DX,EVAL
        INTEGER(4) :: I,J,K,N
        
        N = N_NIST
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(T <= NIST_T(1)) then
          EVAL = ( NIST_Cp(1)-NIST_Cp(2) )/( NIST_T(1)-NIST_T(2) ) *(T-NIST_T(1))+NIST_Cp(1)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its minimum value in the given NIST table.', NIST_T(1), T
          STOP
        end if
        if(T >= NIST_T(n)) then
          EVAL = ( NIST_Cp(n)-NIST_Cp(n-1) )/( NIST_T(n)-NIST_T(n-1) ) *(T-NIST_T(n))+NIST_Cp(n)
          WRITE(*,'(A,2F12.5)') '# Warning: T reaches its maximum value in the given NIST table.', NIST_T(N), T
          STOP
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(T < NIST_T(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = T - NIST_T(i)
        EVAL = NIST_Cp(i) + dx*(NIST_TCp_b(i) + dx*(NIST_TCp_c(i) + dx*NIST_TCp_d(i)))
        IF(DABS(T-1.0_WP).LT.1.0E-10_WP) EVAL = 1.0_WP !test
        return
    end subroutine
    
!==================lib============================================================
    
    !==============================================================================
    SUBROUTINE CUBIC_SPLINE (N,X,Y,B,C,D)
!---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCREET FONCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!     REFERENCE:
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!---------------------------------------------------------------------
        USE WPRECISION
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: N
        REAL(WP),INTENT(IN)   :: X(N), Y(N)
        REAL(WP),INTENT(OUT)  :: B(N), C(N), D(N)
      
        INTEGER(4) :: NM1, I, L
        REAL(WP)    :: T
      
        NM1 = N-1
        IF (N.LT.2) RETURN
        IF (N.LT.3) GO TO 50

!       BUILD THE TRIDIAGONAL SYSTEM
!        B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)
    
        D(1) = X(2)-X(1)
        C(2) = (Y(2)-Y(1))/D(1)
        DO 10 I = 2,NM1
        D(I) = X(I+1)-X(I)
        B(I) = 2.0_WP*(D(I-1)+D(I))
        C(I+1) = (Y(I+1)-Y(I))/D(I)
        C(I) = C(I+1)-C(I)
   10   CONTINUE

!     CONDITIONS AT LIMITS
!     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

        B(1) = -D(1)
        B(N) = -D(N-1)
        C(1) = 0.0_wp
        C(N) = 0.0_wp
        IF (N.EQ.3) GO TO 15
        C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
        C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
        C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
        C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))

!     FORWARD ELIMINATION

   15   DO 20 I = 2,N
            T = D(I-1)/B(I-1)
            B(I) = B(I)-T*D(I-1)
            C(I) = C(I)-T*C(I-1)
   20   CONTINUE

!     BACK SUBSTITUTION

        C(N) = C(N)/B(N)
        DO 30 L = 1,NM1
            I = N-L
            C(I) = (C(I)-D(I)*C(I+1))/B(I)
   30   CONTINUE

!     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

        B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.0_WP*C(N))
        DO 40 I = 1,NM1
            B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.0_WP*C(I))
            D(I) = (C(I+1)-C(I))/D(I)
            C(I) = 3.0_WP*C(I)
   40   CONTINUE
        C(N) = 3.0_WP*C(N)
        D(N) = D(NM1)
        RETURN

!     CAS N = 2

   50   B(1) = (Y(2)-Y(1))/(X(2)-X(1))
        C(1) = 0.0_wp
        D(1) = 0.0_wp
        B(2) = B(1)
        C(2) = 0.0_wp
        D(2) = 0.0_wp
        RETURN
    
    END SUBROUTINE

!#################################################################
    SUBROUTINE FIND_LOCAL_MAXIMA(N, A, JMAX, NMAXL, AMAXL, JMIN,  NMINL, AMINL)
        USE WPRECISION
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: N
        REAL(WP), INTENT (IN) :: A(N)
        INTEGER(4),INTENT(OUT) :: JMAX, JMIN           ! how many peaks
        !=====max. 10 peaks max. and 10 peaks of mini.
        REAL(WP), INTENT (OUT) :: AMAXL(10), AMINL(10) ! peak values
        INTEGER(4),INTENT(OUT) :: NMAXL(10), NMINL(10) ! peaks location
        
        REAL(WP) :: MULTP, DL, DR
        INTEGER(4) :: I
        
        JMAX = 0
        JMIN = 0
        DO I = 2, N-1
            DL   = A(I)-A(I-1)
            DR   = A(I)-A(I+1)
            MULTP= DL * DR
            IF(MULTP.GT.0.0_WP) THEN
                IF(DL .GT. 0.0_WP) THEN
                    JMAX = JMAX + 1
                    AMAXL(JMAX) = A(I)
                    NMAXL(JMAX) = I
                END IF
                IF(DL .LT. 0.0_WP) THEN
                    JMIN = JMIN + 1
                    AMINL(JMIN) = A(I)
                    NMINL(JMIN) = I
                END IF
            END IF
        END DO
        
        RETURN
    END SUBROUTINE 
        
        
    
    
    
