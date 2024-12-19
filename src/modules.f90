!************************************************************************
    MODULE WPRECISION
        INTEGER,PARAMETER  :: WP=8 !KIND(0.0D0) !WORKING PRECESION
    END MODULE WPRECISION
!************************************************************************
    MODULE mpi_info
        !include 'mpif.h'
        use MPI
     
        INTEGER(4) ::  MYID
        INTEGER(4) ::  IERROR
        INTEGER(4) ::  SIZE
        INTEGER(4) ::  ICOMM
        integer(4) ::  STS(MPI_STATUS_SIZE)
     
        INTEGER(4) ::  NPSLV
        INTEGER(4) ::  NPTOT
       
      END MODULE mpi_info
      
!********************************MODULE**********************************
    MODULE cparam
        use WPRECISION
        use mpi_info
        INTEGER(4),PARAMETER ::  NDV=3 
!>      @param NDV is three directions

        REAL(WP), PARAMETER  :: REALMIN = 1.0E-14_WP 
        REAL(WP), PARAMETER  :: REALMAX = 1.0E+14_WP 
        REAL(WP), PARAMETER  :: Cmu0=0.09_wp
        REAL(WP), PARAMETER  :: TK2C=273.15_WP
        
        INTEGER(4),PARAMETER :: IINLET = -1 
        INTEGER(4),PARAMETER :: IOULET =  1
        INTEGER(4),PARAMETER :: IMAIND =  0
        INTEGER(4),PARAMETER :: IALLDM =  2
        
        INTEGER(4),PARAMETER :: IBCWALL = 1
        INTEGER(4),PARAMETER :: IBCPERI = 3
        INTEGER(4),PARAMETER :: IBCNEUM = 2
        
        
        INTEGER(4),PARAMETER :: Ibotwall =  1
        INTEGER(4),PARAMETER :: Itopwall =  2
        
        INTEGER(4),PARAMETER ::  ITG = 1
        INTEGER(4),PARAMETER ::  IIO = 2     
        
        
        INTEGER(4),PARAMETER ::  T_Asymptotic_Average=1
        INTEGER(4),PARAMETER ::  T_Summing_average=2
        
        INTEGER(4),PARAMETER ::  visimplicit = 1
        INTEGER(4),PARAMETER ::  visexplicit = 0

        INTEGER(4),PARAMETER ::  search_table = 1
        INTEGER(4),PARAMETER ::  idealgas_law = 2
        
        INTEGER(4),PARAMETER ::  flgxz  = 1
        INTEGER(4),PARAMETER ::  flgxzt = 2
        
        INTEGER(4),PARAMETER ::  monatomic_gas = 1
        INTEGER(4),PARAMETER ::  diatomic_gas  = 2
        INTEGER(4),PARAMETER ::  trivalence_gas= 3
        
        INTEGER(4),PARAMETER :: powerlawflg=1
        INTEGER(4),PARAMETER :: sutherlandflg=2
        
        INTEGER(4) :: icase
        INTEGER(4),PARAMETER :: ICHANL =  1
        INTEGER(4),PARAMETER :: IPIPEC =  2
        INTEGER(4),PARAMETER :: IANNUL =  3
        INTEGER(4),PARAMETER :: IBOX3P =  4
        
        integer(4)    :: thermlflg
        INTEGER(4)              ::  MGRID, JINI
        INTEGER(4),ALLOCATABLE ::  JGMOV(:)
        
        INTEGER(4) ::  NCL1_tg
        INTEGER(4) ::  NCL1_io
        INTEGER(4) ::  NCL2
        INTEGER(4) ::  NCL3
        
        INTEGER(4) :: NCL1S
        INTEGER(4) :: NCL1E

        INTEGER(4) ::  NND1_tg
        INTEGER(4) ::  NND1_io      
        INTEGER(4) ::  NND2
        INTEGER(4) ::  NND3 
             
        INTEGER(4) ::  Kronecker_delta(NDV,NDV)
        
      
        LOGICAL    ::  TGFLOWflg, IOFLOWflg
        
        INTEGER(4) :: MEMPC_byte = 0  ! memory per core in the unit of byte
      
      END MODULE cparam
!********************************MODULE**********************************
      MODULE wrt_info
        use cparam
        character(8)  :: date
        character(10) :: time
      
        character(18)   :: fllog
        integer(4)      :: logflg_tg = 200
        integer(4)      :: logflg_io = 6
        
        
        CHARACTER(LEN=256)  :: WRT_RST_FNM_TG(NDV+1)
        CHARACTER(LEN=256)  :: WRT_RST_FNM_IO(NDV+1)
        CHARACTER(LEN=256)  :: WRT_AVE_FNM_TG
        CHARACTER(LEN=256)  :: WRT_AVE_FNM_IO
        
        character(11) :: dir1='1_instant_D'
        character(11) :: dir2='2_averagd_D'
        character(11) :: dir3='3_intpcrs_D'
        character(14) :: dir4='4_rslt_ave_plt'
        character(14) :: dir5='5_rslt_ins_plt' 
        character(6)  :: dir6='6_logs' 
        character(6)  :: dir7='7_srcs' 
        
        CHARACTER(LEN=256)  :: filepath1='1_instant_D/'
        CHARACTER(LEN=256)  :: filepath2='2_averagd_D/'
        CHARACTER(LEN=256)  :: filepath3='3_intpcrs_D/'
        CHARACTER(LEN=256)  :: filepath4='4_rslt_ave_plt/'
        CHARACTER(LEN=256)  :: filepath5='5_rslt_ins_plt/'

      
      END MODULE wrt_info

!********************************MODULE**********************************
      MODULE mesh_info
        use cparam 
        
        INTEGER(4),ALLOCATABLE  :: N2DO(:)
        INTEGER(4),ALLOCATABLE  :: N3DO(:)
        
        INTEGER(4),ALLOCATABLE  :: JDSWT(:)
        INTEGER(4),ALLOCATABLE  :: JDEWT(:)
        INTEGER(4),ALLOCATABLE  :: KDSWT(:)
        INTEGER(4),ALLOCATABLE  :: KDEWT(:)
        INTEGER(4),ALLOCATABLE  :: JCL2G(:)
        INTEGER(4),ALLOCATABLE  :: KCL2G(:)
        INTEGER(4),ALLOCATABLE  :: JG2IP(:)
        INTEGER(4),ALLOCATABLE  :: JG2LC(:)
        
        INTEGER(4),ALLOCATABLE  :: IPV_tg(:)
        INTEGER(4),ALLOCATABLE  :: IMV_tg(:)
        INTEGER(4),ALLOCATABLE  :: IPV_io(:)
        INTEGER(4),ALLOCATABLE  :: IMV_io(:)
        INTEGER(4),ALLOCATABLE  :: KPV(:)
        INTEGER(4),ALLOCATABLE  :: KMV(:)
        INTEGER(4),ALLOCATABLE  :: JLPV(:)
        INTEGER(4),ALLOCATABLE  :: JLMV(:)
        INTEGER(4),ALLOCATABLE  :: JGPV(:)
        INTEGER(4),ALLOCATABLE  :: JGMV(:)
        INTEGER(4),ALLOCATABLE  :: KSYM(:)
        
        
        REAL(WP)      :: HX_tg, HX_io, HZ, HYB, HYT
        REAL(WP)      :: ALX1(2),ALX2,ALX3
        REAL(WP)      :: DX, DXI, DXQI
        REAL(WP)      :: DZ, DZI, DZQI
        REAL(WP)      :: VL1313_tg, VL1313_io
        
        REAL(WP),ALLOCATABLE  :: CFLVIS(:)
        
        REAL(WP),ALLOCATABLE  :: AMPH(:)
        REAL(WP),ALLOCATABLE  :: ACPH(:)
        REAL(WP),ALLOCATABLE  :: APPH(:)
        
        REAL(WP)  :: AMPH0
        REAL(WP)  :: APPH0
        !REAL(WP),ALLOCATABLE  :: DPDYWAL(:,:,:)
        
        REAL(WP)               :: AMVR1, APVRN
        
        REAL(WP),ALLOCATABLE  :: AMVR(:,:)
        REAL(WP),ALLOCATABLE  :: ACVR(:,:)
        REAL(WP),ALLOCATABLE  :: APVR(:,:)
        
        REAL(WP),ALLOCATABLE  :: XND_tg(:)
        REAL(WP),ALLOCATABLE  :: XND_io(:)
        REAL(WP),ALLOCATABLE  :: YND(:)
        REAL(WP),ALLOCATABLE  :: ZND(:)
        REAL(WP),ALLOCATABLE  :: ZCC(:)
            
        REAL(WP),ALLOCATABLE  :: RCCI1(:)
        REAL(WP),ALLOCATABLE  :: RNDI1(:) 
        REAL(WP),ALLOCATABLE  :: RCCI2(:)
        REAL(WP),ALLOCATABLE  :: RNDI2(:) 
      
        REAL(WP),ALLOCATABLE  :: YCC(:)
        
        REAL(WP),ALLOCATABLE  :: DYFI(:)
        REAL(WP),ALLOCATABLE  :: DYCI(:)
        
        !===============WEIGHTING FACTORS FOR X, Y, Z AVERAGE (NODE VS CELL CENTRE)==
        REAL(WP)      :: XND2CL
        REAL(WP)      :: YND2CL
        REAL(WP)      :: ZND2CL
        
        REAL(WP),ALLOCATABLE  :: Xcc_tg(:)
        REAL(WP),ALLOCATABLE  :: Xcc_io(:)
        
        REAL(WP),ALLOCATABLE  :: YCL2ND_WFF(:)
        REAL(WP),ALLOCATABLE  :: YCL2ND_WFB(:)
        
    END MODULE mesh_info


!********************************MODULE**********************************
    MODULE BC_INLOUT_info
        use cparam

        REAL(WP)  :: U_OUTLET
        REAL(WP), ALLOCATABLE :: BC_CONV0(:,:,:)
        REAL(WP), ALLOCATABLE :: BC_CONV0_ENG(:,:)
        REAL(WP), ALLOCATABLE :: BC_TDMA (:,:,:, :)
        REAL(WP), ALLOCATABLE :: BC_U_SSTAR(:,:,:)
      
    END MODULE BC_INLOUT_info


!********************************MODULE**********************************
    MODULE init_info
        use cparam
        
        INTEGER(4)  :: ITERG0_io
        INTEGER(4)  :: ITERG0_TG
        INTEGER(4)  :: ITERG0
        
        INTEGER(4)  :: ITERG
        
        REAL(WP)      :: phyTIME_io
        REAL(WP)      :: phyTIME_TG
        REAL(WP)      :: phyTIME
        

        
        INTEGER(4)  :: NSST
        INTEGER(4)  :: NFLOW
        INTEGER(4)  :: RSTflg_tg, RSTflg_io, RST_type_flg, RST_time_set_flg
        INTEGER(4)  :: pprocessonly, ppinst, ppspectra, ppdim, pp_instn_sz
        REAL(WP), ALLOCATABLE :: pp_instn_tim(:)
        
        INTEGER(4)  :: ISTR2
        character(64) :: teczonename
        
        INTEGER(4)   :: FLOWDRVTP
        
        INTEGER(4) ::  weightedpressure
        INTEGER(4) ::  visthemflg
        INTEGER(4) ::  thermoStat
        
        REAL(WP)      :: PI
        REAL(WP)      :: DT     ! real time step
        REAL(WP)      :: DT0    ! given time step
        REAL(WP)      :: DTMIN  ! the minimum time step
        
        REAL(WP)      :: CPUTIME, CPUTIME_TMP
        REAL(WP)      :: CFLGV
        
        
        REAL(WP)      :: RSTtim_tg
        REAL(WP)      :: RSTtim_io
        REAL(WP)      :: TSTOP
        REAL(WP)      :: TLGRE

        REAL(WP)      :: DTSAVE      
        REAL(WP)      :: DTSAVE1 
        REAL(WP)      :: DTSTATEC
        REAL(WP)      :: DTTECCK
        INTEGER(4)      :: ITSCN
        
        REAL(WP)      :: TSTAV1, TSTAV_RESET

        REAL(WP)      :: CFGV
   
        REAL(WP)      :: REINI
        REAL(WP)      :: REN
        REAL(WP)      :: PRT0
        
        REAL(WP)      :: AREA_INLET, VOLM_tg, VOLM_IO
        REAL(WP)      :: CVISC
 
        REAL(WP)      :: VPERG
        REAL(WP)      :: SVPERG
        
        REAL(WP)      :: STR2
        REAL(WP)      :: TGAM(0:3), TROH(0:3), TALP(0:3)
        
        REAL(WP),ALLOCATABLE :: Vini(:)
        
        
        INTEGER(4)    :: BCX_tg(2)
        INTEGER(4)    :: BCZ(2)
        INTEGER(4)    :: BCY(2)
        
        INTEGER(4)    :: BCX_io(2)
        
        
        
    END MODULE init_info
      
!********************************MODULE**********************************
    MODULE flow_info
        use cparam
        use BC_INLOUT_info
        
        REAL(WP)      :: FACOE_TG, FACOE_IO
        REAL(WP)      :: U1MEAN_tg,  U1MAXX_tg, U1MEAN_work_tg,  U1MAXX_work_tg
        
        REAL(WP)      :: G1MAXX_io, G1MAXX_work_io
        REAL(WP)      :: G1RATE_io, G1RATE_work_io, G1BULK_work_io
        REAL(WP)      :: T1MAXX_io, T1MAXX_work_io
        REAL(WP)      :: H1RATE_io, H1RATE_work_io, H1BULK_work_io, T1BULK_work_io
        
        
        REAL(WP)      :: CHK_Mass_CONSV0, CHK_Mass_CONSV0_SCALING
        REAL(WP)      :: CHK_ENEG_CONSV0, CHK_ENEG_TOTAL
        
        !===============TURBULENCE==GENERATOR=================================
        REAL(WP)      :: MAXDIVGV_tg(2)
        REAL(WP)      :: CFLMM_tg
        REAL(WP)      :: CFLMM
        REAL(WP)      :: VMAX_tg(3), VMIN_tg(3)
        
        !=========PRIMARY==VARIABLES==DIMENSIONLESS==========
        REAL(WP),ALLOCATABLE :: Q_tg   (:,:,:,:)
        REAL(WP),ALLOCATABLE :: PR_tg  (:,:,:)
        
        !=========VARIABLES==IN==CALCULATION=================
        REAL(WP),ALLOCATABLE :: RHS_tg     (:,:,:)
        REAL(WP),ALLOCATABLE :: CONVH0_tg  (:,:,:,:)
        REAL(WP),ALLOCATABLE :: DPH_tg     (:,:,:)
        REAL(WP),ALLOCATABLE :: RHSLLPHI_tg(:,:,:)
        
        !========ASSISTANT==VARIABLES==IN==CALCULATION=======
        !REAL(WP),ALLOCATABLE :: F   (:,:,:)
        REAL(WP),ALLOCATABLE :: Qtmp_tg(:,:,:)
        
        !===========MAIN==DOMAIN=================================================
        REAL(WP)      :: MAXDIVGV_io(3) ! 1formaindomain,2formaindomain,3forb.c.outlet.
        REAL(WP)      :: VMAX_io(3),VMIN_io(3)
        REAL(WP)      :: CFLMM_io
        
        !=========PRIMARY==VARIABLES==DIMENSIONLESS==========
        REAL(WP),ALLOCATABLE :: Q_io   (:,:,:,:)
        REAL(WP),ALLOCATABLE :: G_io   (:,:,:,:)
        REAL(WP),ALLOCATABLE :: G0_io   (:,:,:,:)
        REAL(WP),ALLOCATABLE :: PR_io  (:,:,:)
        REAL(WP),ALLOCATABLE :: PR0_io  (:,:,:,:)
        REAL(WP),parameter   :: pres_epslon = 0.001_wp
        
        !=========VARIABLES==IN==CALCULATION=================
        REAL(WP),ALLOCATABLE :: RHS_io     (:,:,:)

        !REAL(WP),ALLOCATABLE :: VISCS0_io  (:,:,:,:)
        REAL(WP),ALLOCATABLE :: EXPLT0_io  (:,:,:,:)
        REAL(WP),ALLOCATABLE :: DPH_io     (:,:,:)
        REAL(WP),ALLOCATABLE :: RHSLLPHI_io(:,:,:)
        REAL(WP),ALLOCATABLE :: DivU_io(:,:,:)
        
        
        !========ASSISTANT==VARIABLES==IN==CALCULATION=======
        REAL(WP),ALLOCATABLE :: Qtmp_io(:,:,:)
        
        REAL(WP),ALLOCATABLE :: UU(:,:,:)
        
    END MODULE flow_info
      
!*************************************************************************
    MODULE thermal_info
        use cparam
        !==============INI GIVEN ================
        INTEGER(4)    :: FLOWDIR
        INTEGER(4)    :: GRAVDIR
        CHARACTER(64) :: NISTFLNM
        integer(4),parameter    :: isoFluxWall=1
        integer(4),parameter    :: isoThermalWall=2
        integer(4)    :: BCWALLHEAT(2)
        
        REAL(WP)      :: CTHECD
        
        
        
        REAL(WP)      :: G_A
        REAL(WP)      :: F_A
        INTEGER(4)    :: IBuoF(3)
        
        REAL(WP)      :: Gr0(2)
        REAL(WP)      :: BO0(2,2)
        
        !===========for perfect gas==============
        INTEGER(4)    :: IdealGasType
        INTEGER(4)    :: fstatetype ! flow state equation.
        REAL(WP)      :: R0
        REAL(WP)      :: CV0
        REAL(WP)      :: gamma0
        REAL(WP)      :: Cvhat0
        REAL(WP)      :: RHOU20
        
        REAL(WP)      :: PL_MT  = 0.67_WP
        REAL(WP)      :: PL_CpT = 0.095_WP
        REAL(WP)      :: PL_HT  = 1.095_WP
        REAL(WP)      :: PL_KT  = 0.805_WP
         
        !==============INI FOR REFEENCE STATE==DIMENSIONAL====
        REAL(WP)      :: Ti
        REAL(WP)      :: L0  !  (m)
        REAL(WP)      :: U0  !  (m/s)
        REAL(WP)      :: P0  !  (Pa)
        REAL(WP)      :: T0  !  (K)
        REAL(WP)      :: H0  !  (Pa)
        REAL(WP)      :: D0  !  (J/Kg)
        REAL(WP)      :: M0  !  (Pa-s)
        REAL(WP)      :: K0  !  (W/m-K)
        REAL(WP)      :: CP0 ! (J/Kg-K)
        REAL(WP)      :: B0  !  (1/K)
        REAL(WP)      :: WHEAT0_DIM(2) ! (W/m2) given heat flux rate on the wall or ! (K) given constant wall temperature
        
        REAL(WP)      :: WHEAT0_UNDIM(2) 
        
        !REAL(WP)      :: WHFLUX_UND(2)
        REAL(WP)       :: RHOHmax, T4MaxRhoh
        REAL(WP)       :: CpMax, T4CpMax
        !===============VARIABLES====DIMENSIONALESS===============
        REAL(WP),ALLOCATABLE  :: RHOH(:,:,:)       
        REAL(WP),ALLOCATABLE  :: DENSITY(:,:,:)   
        REAL(WP),ALLOCATABLE  :: TEMPERATURE(:,:,:) 
        REAL(WP),ALLOCATABLE  :: VISCOUSITY(:,:,:) 
        REAL(WP),ALLOCATABLE  :: VISCOUSITY0(:,:,:) 
        REAL(WP),ALLOCATABLE  :: THERMCONDT(:,:,:) 
        REAL(WP),ALLOCATABLE  :: ENTHALPY(:,:,:)
        REAL(WP),ALLOCATABLE  :: RHO_STG(:,:,:,:)   
        REAL(WP),ALLOCATABLE  :: DRHOI_STG(:,:,:,:)   
        REAL(WP),ALLOCATABLE  :: MU_STG(:,:,:,:)   
        
        !==============VARIABLES====IN=CALCULATION================
        REAL(WP),ALLOCATABLE  :: RHS_ENERGY (:,:,:)
        REAL(WP),ALLOCATABLE  :: RHS_ENERGY0(:,:,:)
        REAL(WP),ALLOCATABLE  :: DENSITY0   (:,:,:)
        REAL(WP),ALLOCATABLE  :: DrhoDtP    (:,:,:)
        REAL(WP),ALLOCATABLE  :: RHOH0   (:,:,:)
        REAL(WP),ALLOCATABLE  :: WALLFLUX   (:,:)
        
        REAL(WP),ALLOCATABLE  :: DENSITY1   (:,:,:)
        
        !==============NIST TABLE INFO=====================
        INTEGER(4) :: N_NIST
        REAL(WP),ALLOCATABLE :: NIST_H(:)
        REAL(WP),ALLOCATABLE :: NIST_T(:)
        REAL(WP),ALLOCATABLE :: NIST_D(:)
        REAL(WP),ALLOCATABLE :: NIST_M(:)
        REAL(WP),ALLOCATABLE :: NIST_K(:)
        REAL(WP),ALLOCATABLE :: NIST_CP(:)
        REAL(WP),ALLOCATABLE :: NIST_B(:)
        
        
        REAL(WP),ALLOCATABLE :: NIST_HT_B(:)
        REAL(WP),ALLOCATABLE :: NIST_HD_B(:)
        REAL(WP),ALLOCATABLE :: NIST_HM_B(:)
        REAL(WP),ALLOCATABLE :: NIST_HK_B(:)
        REAL(WP),ALLOCATABLE :: NIST_HCP_B(:)
        REAL(WP),ALLOCATABLE :: NIST_HB_B(:)
        
        
        REAL(WP),ALLOCATABLE :: NIST_HT_C(:)
        REAL(WP),ALLOCATABLE :: NIST_HD_C(:)
        REAL(WP),ALLOCATABLE :: NIST_HM_C(:)
        REAL(WP),ALLOCATABLE :: NIST_HK_C(:)
        REAL(WP),ALLOCATABLE :: NIST_HCP_C(:)
        REAL(WP),ALLOCATABLE :: NIST_HB_C(:)
        
        
        REAL(WP),ALLOCATABLE :: NIST_HT_D(:)
        REAL(WP),ALLOCATABLE :: NIST_HD_D(:)
        REAL(WP),ALLOCATABLE :: NIST_HM_D(:)
        REAL(WP),ALLOCATABLE :: NIST_HK_D(:)
        REAL(WP),ALLOCATABLE :: NIST_HCP_D(:)
        REAL(WP),ALLOCATABLE :: NIST_HB_D(:)
        
        REAL(WP),ALLOCATABLE :: NIST_TH_B(:)
        REAL(WP),ALLOCATABLE :: NIST_TD_B(:)
        REAL(WP),ALLOCATABLE :: NIST_TM_B(:)
        REAL(WP),ALLOCATABLE :: NIST_TK_B(:)
        REAL(WP),ALLOCATABLE :: NIST_TCP_B(:)
        REAL(WP),ALLOCATABLE :: NIST_TB_B(:)
        
        REAL(WP),ALLOCATABLE :: NIST_TH_C(:)
        REAL(WP),ALLOCATABLE :: NIST_TD_C(:)
        REAL(WP),ALLOCATABLE :: NIST_TM_C(:)
        REAL(WP),ALLOCATABLE :: NIST_TK_C(:)
        REAL(WP),ALLOCATABLE :: NIST_TCP_C(:)
        REAL(WP),ALLOCATABLE :: NIST_TB_C(:)
        
        REAL(WP),ALLOCATABLE :: NIST_TH_D(:)
        REAL(WP),ALLOCATABLE :: NIST_TD_D(:)
        REAL(WP),ALLOCATABLE :: NIST_TM_D(:)
        REAL(WP),ALLOCATABLE :: NIST_TK_D(:)
        REAL(WP),ALLOCATABLE :: NIST_TCP_D(:)
        REAL(WP),ALLOCATABLE :: NIST_TB_D(:)
        
        
        !==============SOME CONSTANT==========================
        INTEGER(4) :: NTHERMAL = 6
        
        !========================
        REAL(WP),ALLOCATABLE :: T_WAL_GV(:,:)
        REAL(WP),ALLOCATABLE :: H_WAL_GV(:,:)
        REAL(WP),ALLOCATABLE :: D_WAL_GV(:,:)
        REAL(WP),ALLOCATABLE :: M_WAL_GV(:,:)
        REAL(WP),ALLOCATABLE :: K_WAL_GV(:,:)
        REAL(WP),ALLOCATABLE :: Cp_WAL_GV(:,:)
      
    END MODULE
    
!********************************MODULE**********************************      


!    MODULE WALLSTRESS_info
!        use cparam
      
!        REAL(WP),ALLOCATABLE     :: Cf_LW_io(:)
!        REAL(WP),ALLOCATABLE     :: Cf_UW_io(:)
!        REAL(WP),ALLOCATABLE     :: Cf_LW_tg(:)
!        REAL(WP),ALLOCATABLE     :: Cf_UW_tg(:)
!        INTEGER(4)               :: NCOUNT1 = 0
!        INTEGER(4)               :: TECFLG 
!    END MODULE WALLSTRESS_info
      
!********************************MODULE**********************************      
    MODULE postprocess_info
        use cparam
        
        INTEGER(4)  :: NSTATIS_tg
        INTEGER(4)  :: NSTATIS_IO
        
        !===============TG part====space averaged===========================
        REAL(WP),ALLOCATABLE     :: U1xzL_tg(:,:), U1xztL_tg(:,:) !Ui
        REAL(WP),ALLOCATABLE     :: UPxzL_tg(:,:), UPxztL_tg(:,:) !UiP
        REAL(WP),ALLOCATABLE     :: U2xzL_tg(:,:), U2xztL_tg(:,:) !UiUj
        REAL(WP),ALLOCATABLE     :: U3xzL_tg(:,:), U3xztL_tg(:,:) !UiUjUk
        
        REAL(WP),ALLOCATABLE     :: DVDL1xzL_tg(:,:,:), DVDL1xztL_tg(:,:,:) ! dUi/dXj
        REAL(WP),ALLOCATABLE     :: DVDLPxzL_tg(:,:,:), DVDLPxztL_tg(:,:,:) !PdUi/dXj
        REAL(WP),ALLOCATABLE     :: DVDL2xzL_tg(:,:,:), DVDL2xztL_tg(:,:,:) ! dUi/dXj * dUi/dXj
        
        !===============IO part (non-periodic)====space averaged================
        REAL(WP),ALLOCATABLE     :: U1zL_io(:,:,:),  U1ztL_io(:,:,:)!Ui
        REAL(WP),ALLOCATABLE     :: G1zL_io(:,:,:),  G1ztL_io(:,:,:)!Gi
        REAL(WP),ALLOCATABLE     :: UPzL_io(:,:,:),  UPztL_io(:,:,:)!UiP
        
        REAL(WP),ALLOCATABLE     :: U2zL_io(:,:,:),  U2ztL_io(:,:,:)!UiUj
        REAL(WP),ALLOCATABLE     :: UGzL_io(:,:,:),  UGztL_io(:,:,:)!UiGj
        
        REAL(WP),ALLOCATABLE     :: U3zL_io(:,:,:),  U3ztL_io(:,:,:)!UiGjUk
        REAL(WP),ALLOCATABLE     :: UGUzL_io(:,:,:), UGUztL_io(:,:,:)!UiGjUk
        
        REAL(WP),ALLOCATABLE     :: DVDL1zL_IO(:,:,:,:), DVDL1ztL_IO(:,:,:,:) ! dUi/dXj
        REAL(WP),ALLOCATABLE     :: DVDLPzL_IO(:,:,:,:), DVDLPztL_IO(:,:,:,:) !PdUi/dXj
        REAL(WP),ALLOCATABLE     :: DVDL2zL_IO(:,:,:,:), DVDL2ztL_IO(:,:,:,:) ! dUi/dXj * dUi/dXj
        
        REAL(WP),ALLOCATABLE     :: T1zL_io(:,:), T1ztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: D1zL_io(:,:), D1ztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: H1zL_io(:,:), H1ztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: M1zL_io(:,:), M1ztL_io(:,:)
        
        REAL(WP),ALLOCATABLE     :: T2zL_io(:,:), T2ztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: D2zL_io(:,:), D2ztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: H2zL_io(:,:), H2ztL_io(:,:)
        
        REAL(WP),ALLOCATABLE     :: DHzL_io(:,:), DHztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: PHzL_io(:,:), PHztL_io(:,:)
        
        REAL(WP),ALLOCATABLE     :: DVDL1MzL_io(:,:,:,:),    DVDL1MztL_io(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: DVDL1MHzL_io(:,:,:,:),   DVDL1MHztL_io(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: DVDL1MUzL_io(:,:,:,:,:), DVDL1MUztL_io(:,:,:,:,:)
        REAL(WP),ALLOCATABLE     :: DVDL2MzL_io(:,:,:,:),    DVDL2MztL_io(:,:,:,:)
        
        
        REAL(WP),ALLOCATABLE     :: UHzL_io(:,:,:), UHztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: GHzL_io(:,:,:), GHztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: U2DHzL_io(:,:,:), U2DHztL_io(:,:,:)
        
        REAL(WP),ALLOCATABLE     :: DhDL1zL_io(:,:,:),         DhDL1ztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: DhDLPzL_io(:,:,:),         DhDLPztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: DTDLKzL_io(:,:,:),         DTDLKztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: DTDLKUzL_io(:,:,:,:),      DTDLKUztL_io(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: DTDLKDVDLzL_io(:,:,:,:,:), DTDLKDVDLztL_io(:,:,:,:,:)
        REAL(WP),ALLOCATABLE     :: DhDLMDVDLzL_io(:,:,:,:,:), DhDLMDVDLztL_io(:,:,:,:,:)
        
        
        !===============IO part (periodic)====space averaged================
        REAL(WP),ALLOCATABLE     :: U1xzL_io(:,:),  U1xztL_io(:,:) !Ui
        REAL(WP),ALLOCATABLE     :: G1xzL_io(:,:),  G1xztL_io(:,:)!Gi
        REAL(WP),ALLOCATABLE     :: UPxzL_io(:,:),  UPxztL_io(:,:)!UiP

        REAL(WP),ALLOCATABLE     :: U2xzL_io(:,:),  U2xztL_io(:,:)!UiUj
        REAL(WP),ALLOCATABLE     :: UGxzL_io(:,:),  UGxztL_io(:,:)!UiGj
        
        REAL(WP),ALLOCATABLE     :: U3xzL_io(:,:),  U3xztL_io(:,:)!UiUjUk
        REAL(WP),ALLOCATABLE     :: UGUxzL_io(:,:), UGUxztL_io(:,:)!UiGjUk
        
        REAL(WP),ALLOCATABLE     :: DVDL1xzL_IO(:,:,:), DVDL1xztL_IO(:,:,:) ! dUi/dXj
        REAL(WP),ALLOCATABLE     :: DVDLPxzL_IO(:,:,:), DVDLPxztL_IO(:,:,:) ! PdUi/dXj
        REAL(WP),ALLOCATABLE     :: DVDL2xzL_IO(:,:,:), DVDL2xztL_IO(:,:,:) ! dUi/dXj * dUi/dXj
        
        REAL(WP),ALLOCATABLE     :: T1xzL_io(:), T1xztL_io(:)
        REAL(WP),ALLOCATABLE     :: D1xzL_io(:), D1xztL_io(:)
        REAL(WP),ALLOCATABLE     :: H1xzL_io(:), H1xztL_io(:)
        REAL(WP),ALLOCATABLE     :: M1xzL_io(:), M1xztL_io(:)

        REAL(WP),ALLOCATABLE     :: T2xzL_io(:), T2xztL_io(:)
        REAL(WP),ALLOCATABLE     :: D2xzL_io(:), D2xztL_io(:)
        REAL(WP),ALLOCATABLE     :: H2xzL_io(:), H2xztL_io(:)
        
        REAL(WP),ALLOCATABLE     :: DHxzL_io(:), DHxztL_io(:)
        REAL(WP),ALLOCATABLE     :: PHxzL_io(:), PHxztL_io(:)
        
        REAL(WP),ALLOCATABLE     :: DVDL1MxzL_io(:,:,:),    DVDL1MxztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: DVDL1MHxzL_io(:,:,:),   DVDL1MHxztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: DVDL1MUxzL_io(:,:,:,:), DVDL1MUxztL_io(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: DVDL2MxzL_io(:,:,:),    DVDL2MxztL_io(:,:,:)
        
        REAL(WP),ALLOCATABLE     :: UHxzL_io(:,:), UHxztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: GHxzL_io(:,:), GHxztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: U2DHxzL_io(:,:), U2DHxztL_io(:,:)
        
        REAL(WP),ALLOCATABLE     :: DhDL1xzL_io(:,:),         DhDL1xztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: DhDLPxzL_io(:,:),         DhDLPxztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: DTDLKxzL_io(:,:),         DTDLKxztL_io(:,:)
        REAL(WP),ALLOCATABLE     :: DTDLKUxzL_io(:,:,:),      DTDLKUxztL_io(:,:,:)
        REAL(WP),ALLOCATABLE     :: DTDLKDVDLxzL_io(:,:,:,:), DTDLKDVDLxztL_io(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: DhDLMDVDLxzL_io(:,:,:,:), DhDLMDVDLxztL_io(:,:,:,:)
        
        !================variables==================================================
        REAL(WP) :: Hwal_RA(2)
        REAL(WP) :: Hwal_FA(2)
        REAL(WP) :: Twal(2)
        REAL(WP) :: Dwal(2) = 1.0_WP
        REAL(WP) :: Mwal(2) = 1.0_WP
        REAL(WP) :: Kwal(2)
        REAL(WP) :: Cpwal(2)
        REAL(WP) :: qw(2)
        
        !======below values are non-dimensional===============
        
        REAL(WP)    :: Utaw_d_io(2), Utaw_io(2), Utaw_ave_io, Utaw_d_ave_io
        REAL(WP)    :: Tauw_d_io(2), Tauw_io(2), Tauw_ave_io, Tauw_d_ave_io
        REAL(WP)    :: Ret_io(2), Ldist_io(2)
        REAL(WP)    :: Ret_ave_io, DenAvew, VisAvew
        
        
        !===below for quadrant analysis=========
        INTEGER(4),PARAMETER :: QUADHN=9
        REAL(WP)             :: QUADHV(QUADHN)
        
        !==========tg=================================
        REAL(WP),ALLOCATABLE     :: QUADUVxzL_TG(:,:,:), QUADUVxztL_TG(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADVzxzL_TG(:,:,:), QUADVzxztL_TG(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADTKxzL_TG(:,:,:), QUADTKxztL_TG(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADDRxzL_TG(:,:,:), QUADDRxztL_TG(:,:,:)
        
        !==========developping io======================
        REAL(WP),ALLOCATABLE     :: QUADUVzL_IO(:,:,:,:), QUADUVztL_IO(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADVzzL_IO(:,:,:,:), QUADVzztL_IO(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADTKzL_IO(:,:,:,:), QUADTKztL_IO(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADDRzL_IO(:,:,:,:), QUADDRztL_IO(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADDUV1zL_IO(:,:,:,:), QUADDUV1ztL_IO(:,:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADDUV2zL_IO(:,:,:,:), QUADDUV2ztL_IO(:,:,:,:)
        
        !==========periodic io=========================
        REAL(WP),ALLOCATABLE     :: QUADUVxzL_IO(:,:,:), QUADUVxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADVzxzL_IO(:,:,:), QUADVzxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADTKxzL_IO(:,:,:), QUADTKxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADDRxzL_IO(:,:,:), QUADDRxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADDUV1xzL_IO(:,:,:), QUADDUV1xztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: QUADDUV2xzL_IO(:,:,:), QUADDUV2xztL_IO(:,:,:)

        REAL(WP),ALLOCATABLE     :: OCTDUVxzL_IO(:,:,:), OCTDUVxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTDVzxzL_IO(:,:,:), OCTDVzxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTDTKxzL_IO(:,:,:), OCTDTKxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTDDRxzL_IO(:,:,:),  OCTDDRxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTDDUV1xzL_IO(:,:,:), OCTDDUV1xztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTDDUV2xzL_IO(:,:,:), OCTDDUV2xztL_IO(:,:,:)
        
        REAL(WP),ALLOCATABLE     :: OCTTUVxzL_IO(:,:,:), OCTTUVxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTTVzxzL_IO(:,:,:), OCTTVzxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTTTKxzL_IO(:,:,:), OCTTTKxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTTDRxzL_IO(:,:,:), OCTTDRxztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTTDUV1xzL_IO(:,:,:), OCTTDUV1xztL_IO(:,:,:)
        REAL(WP),ALLOCATABLE     :: OCTTDUV2xzL_IO(:,:,:), OCTTDUV2xztL_IO(:,:,:)
        
        !=========DRIVEN FORCE============================
        REAL(WP)                 :: FcDrv_IO
        
        REAL(WP),ALLOCATABLE     :: FUxzL_IO(:,:)
        REAL(WP),ALLOCATABLE     :: FUxztL_IO(:,:)
        
        
        
        !===========spectra======================
        !================Velocity =====================
        REAL(WP), ALLOCATABLE :: R11X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R22X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R33X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R12X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R13X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R23X1_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: R11X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R22X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R33X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R12X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R13X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R23X3_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENE11T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE22T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE33T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE12T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE13T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE23T_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENE11Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE22Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE33Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE12Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE13Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE23Z_xzLa(:, :, :)
        
        !================Voriticity ====================
        REAL(WP), ALLOCATABLE :: V11X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V22X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V33X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V12X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V13X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V23X1_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: V11X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V22X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V33X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V12X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V13X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V23X3_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENV11T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV22T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV33T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV12T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV13T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV23T_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENV11Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV22Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV33Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV12Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV13Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV23Z_xzLa(:, :, :)
        
        
        !===============Voriticity & Velocity=========================
        REAL(WP), ALLOCATABLE :: VO11X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO12X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO13X1_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO21X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO22X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO23X1_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO31X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO32X1_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO33X1_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO11X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO12X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO13X3_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO21X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO22X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO23X3_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO31X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO32X3_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO33X3_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO11T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO12T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO13T_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO21T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO22T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO23T_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO31T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO32T_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO33T_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO11Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO12Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO13Z_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO21Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO22Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO23Z_xzLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO31Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO32Z_xzLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO33Z_xzLa(:, :, :)
        
        !===========spectra======================
        !================Velocity =====================
        REAL(WP), ALLOCATABLE :: R11X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R22X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R33X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R12X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R13X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R23X1_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: R11X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R22X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R33X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R12X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R13X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: R23X3_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENE11T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE22T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE33T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE12T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE13T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE23T_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENE11Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE22Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE33Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE12Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE13Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENE23Z_xztLa(:, :, :)
        
        !================Voriticity ====================
        REAL(WP), ALLOCATABLE :: V11X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V22X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V33X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V12X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V13X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V23X1_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: V11X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V22X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V33X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V12X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V13X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: V23X3_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENV11T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV22T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV33T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV12T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV13T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV23T_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENV11Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV22Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV33Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV12Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV13Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: ENV23Z_xztLa(:, :, :)
        
        
        !===============Voriticity & Velocity=========================
        REAL(WP), ALLOCATABLE :: VO11X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO12X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO13X1_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO21X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO22X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO23X1_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO31X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO32X1_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO33X1_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO11X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO12X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO13X3_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO21X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO22X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO23X3_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: VO31X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO32X3_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: VO33X3_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO11T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO12T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO13T_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO21T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO22T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO23T_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO31T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO32T_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO33T_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO11Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO12Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO13Z_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO21Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO22Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO23Z_xztLa(:, :, :)
        
        REAL(WP), ALLOCATABLE :: EVO31Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO32Z_xztLa(:, :, :)
        REAL(WP), ALLOCATABLE :: EVO33Z_xztLa(:, :, :)
            
    END MODULE postprocess_info 
      
!##############################################################
    MODULE Flow_State_Constant
        use cparam
        
        REAL(WP),Parameter :: SL_MU0 = 1.716E-5_WP  ! Kg/m/s
        REAL(WP),Parameter :: SL_T0  = 273.15_WP    ! K
        REAL(WP),Parameter :: SL_S   = 110.4_WP     ! K
        REAL(WP),Parameter :: SL_C1  = 1.458E-6_WP  ! kg/m/s/sqrt(K)
        
        REAL(WP),Parameter :: SL_K0  = 0.023961_WP  ! W/mK
        
    END MODULE Flow_State_Constant
    
module flatten_index_mod
 implicit none 
 
 interface flatten_index
   module procedure flatten_3d_to_1d
   module procedure flatten_2d_to_1d
 end interface
 
contains

 function flatten_3d_to_1d(i, j, k, Nx, Ny, Nz) result(n)
   integer, intent(in) :: i, j, k, Nx, Ny, Nz
   integer :: n
   n = i + Nx * (j - 1)  + Nx * Ny * (k - 1)
 end function
 
 function flatten_2d_to_1d(i, j, Nx, Ny) result(n)
   integer, intent(in) :: i, j, Nx, Ny
   integer :: n
   n = i + Nx * (j - 1)
 end function
 
end module