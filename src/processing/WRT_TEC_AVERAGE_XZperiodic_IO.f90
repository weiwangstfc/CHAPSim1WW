!================================PERIODIC BELOW========================================================
    MODULE VARS_AVERAGED_XZ_IO
        use thermal_info
        USE wrt_info
        use postprocess_info
        use mesh_info
        use init_info
        use flow_info
        CHARACTER(15) :: PNTIM
        
        INTEGER(4),PARAMETER  :: NX = 5
        
       !=============averaged global Dwta in each processor======================
        REAL(WP),ALLOCATABLE :: U1xztL_F0_io( :, : ) 
        REAL(WP),ALLOCATABLE :: G1xztL_F0_io( :, : )
        REAL(WP),ALLOCATABLE :: UPxztL_F0_io( :, : )
        
        REAL(WP),ALLOCATABLE :: U2xztL_F0_io ( :, : )
        REAL(WP),ALLOCATABLE :: UGxztL_F0_io ( :, : )
        REAL(WP),ALLOCATABLE :: UGUxztL_F0_io( :, : )
        REAL(WP),ALLOCATABLE :: U3xztL_F0_io( :, : )
        
        REAL(WP),ALLOCATABLE :: DVDL1xztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: DVDLPxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: DVDL2xztL_F0_io( :, : , :  )
        
        REAL(WP),ALLOCATABLE :: QuadUVxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: QuadVzxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: QuadTKxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: QuadDRxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: QuadDUV1xztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: QuadDUV2xztL_F0_io( :, : , :  )
        
        REAL(WP),ALLOCATABLE :: OctDUVxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctDVzxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctDTKxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctDDRxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctDDUV1xztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctDDUV2xztL_F0_io( :, : , :  )
        
        REAL(WP),ALLOCATABLE :: OctTUVxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctTVzxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctTTKxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctTDRxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctTDUV1xztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: OctTDUV2xztL_F0_io( :, : , :  )
        
        REAL(WP),ALLOCATABLE :: FUxztL_F0_io( :, :)
        
        !=============averaged global Dwta in each processor thermal======================
        REAL(WP),ALLOCATABLE :: T1xztL_F0_io( : )
        REAL(WP),ALLOCATABLE :: D1xztL_F0_io( : )
        REAL(WP),ALLOCATABLE :: H1xztL_F0_io( : )
        REAL(WP),ALLOCATABLE :: M1xztL_F0_io( : )
        
        REAL(WP),ALLOCATABLE :: T2xztL_F0_io( : )
        REAL(WP),ALLOCATABLE :: D2xztL_F0_io( : )
        REAL(WP),ALLOCATABLE :: H2xztL_F0_io( : )
        
        REAL(WP),ALLOCATABLE :: DHxztL_F0_io( : )
        REAL(WP),ALLOCATABLE :: PHxztL_F0_io( : )
        
        REAL(WP),ALLOCATABLE :: DVDL1MxztL_F0_io ( :, : , :  )
        REAL(WP),ALLOCATABLE :: DVDL1MHxztL_F0_io( :, : , :  )
        REAL(WP),ALLOCATABLE :: DVDL1MUxztL_F0_io( :, : , : , :)
        REAL(WP),ALLOCATABLE :: DVDL2MxztL_F0_io ( :, : , :  )
        
        REAL(WP),ALLOCATABLE :: UHxztL_F0_io( :, : )
        REAL(WP),ALLOCATABLE :: GHxztL_F0_io( :, : )
        REAL(WP),ALLOCATABLE :: U2DHxztL_F0_io( :, : )
        
        REAL(WP),ALLOCATABLE :: DhDL1xztL_F0_io( :, : )
        REAL(WP),ALLOCATABLE :: DhDLPxztL_F0_io( :, : )
        REAL(WP),ALLOCATABLE :: DTDLKxztL_F0_io( :, : )
        REAL(WP),ALLOCATABLE :: DTDLKUxztL_F0_io( :, : , :)
        REAL(WP),ALLOCATABLE :: DTDLKDVDLxztL_F0_io( :, : , :, :)
        REAL(WP),ALLOCATABLE :: DHDLMDVDLxztL_F0_io( :, : , :, :)
        
        !===========================================
        REAL(WP),ALLOCATABLE  :: TauwSD(:)
        REAL(WP),ALLOCATABLE  :: DensSD(:)
        REAL(WP),ALLOCATABLE  :: YWdiSD(:)
        REAL(WP),ALLOCATABLE  :: ViscSD(:)
        
        REAL(WP),ALLOCATABLE  :: CpSD(:)
        REAL(WP),ALLOCATABLE  :: QwSD(:)
        REAL(WP),ALLOCATABLE  :: TwSD(:)
        REAL(WP),ALLOCATABLE  :: HwSD(:)
        
        !========================pp intermediate variables ===============================
        !====RA related==============
        REAL(WP),ALLOCATABLE :: dPdX_RA   (:, :)          !d<p>/dx_i
        REAL(WP),ALLOCATABLE :: ufpf_RA   (:, :)          !<p' u'_i> = <p' u"_i>
        REAL(WP),ALLOCATABLE :: uf2_RA    (:, :, :)       !<u'_i u'_j>
        REAL(WP),ALLOCATABLE :: uf2d_RA   (:, :, :)       !<\rho> <u'_i u'_j>
        REAL(WP),ALLOCATABLE :: uf3_RA    (:, :, :, :)    !<u'_i u'_j u'_k>
        REAL(WP),ALLOCATABLE :: uf3d_RA   (:, :, :, :)    !<\rho> * <u'_i u'_j u'_k>
        REAL(WP),ALLOCATABLE :: UU_RA      (:, :, :)    !{u_i u_j}
        
        REAL(WP),ALLOCATABLE :: dUidXi    (:)             ! dUi/dxi
        REAL(WP),ALLOCATABLE :: StrainTensor(:, :, :)     ! 0.5*(dUi/dxj+dUj/dxi)
        REAL(WP),ALLOCATABLE :: VortcyTensor(:, :, :)     ! 0.5*(dUi/dxj-dUj/dxi)
        REAL(WP),ALLOCATABLE :: Skewness_RA(:,:)
        
        REAL(WP),ALLOCATABLE :: MKE_RA(:)
        REAL(WP),ALLOCATABLE :: TKE_RA(:)
        REAL(WP),ALLOCATABLE :: ufTKEfd_RA(:)
        REAL(WP),ALLOCATABLE :: ufMKEfd_RA(:)
        
        REAL(WP),ALLOCATABLE :: Omega2_RA(:,:)
        REAL(WP),ALLOCATABLE :: Omega_RA2(:,:)
        REAL(WP),ALLOCATABLE :: Omega_rms(:,:)
        
        REAL(WP),ALLOCATABLE :: ANISOTROPY_RA(:, :, :)
        REAL(WP),ALLOCATABLE :: ANISTPinva_RA(:, :)
        REAL(WP),ALLOCATABLE :: LumleyAxis_RA(:, :)
        
        !====FA related==============
        REAL(WP),ALLOCATABLE :: DrivenForce (:)         !
        
        REAL(WP),ALLOCATABLE :: U_FA       (:, :)       !{u_i}
        REAL(WP),ALLOCATABLE :: UU_FA      (:, :, :)    !{u_i u_j}
        REAL(WP),ALLOCATABLE :: dUdX_FA    (:, :, :)    !d{u_i}/dx_j
        
        REAL(WP),ALLOCATABLE :: uff_RA     (:, :)       ! <u"_i>
        REAL(WP),ALLOCATABLE :: uff2_FA    (:, :, :)    !{u"_i u"_j}
        REAL(WP),ALLOCATABLE :: uff2d_FA   (:, :, :)    !<\rho>*{u"_i u"_j}
        REAL(WP),ALLOCATABLE :: uff3_FA    (:, :, :, :) !{u"_i u"_j u"_k}
        REAL(WP),ALLOCATABLE :: uff3d_FA   (:, :, :, :) !<\rho>*{u"_i u"_j u"_k}
        REAL(WP),ALLOCATABLE :: TDIFU_FA   (:, :, :)    !
        REAL(WP),ALLOCATABLE :: dUidXiM      (:)        ! <\mu dUi/dxi>
        REAL(WP),ALLOCATABLE :: StrainTensorM(:, :, :)  ! <0.5*(dUi/dxj+dUj/dxi)*\mu>
        REAL(WP),ALLOCATABLE :: VortcyTensorM(:, :, :)  ! <0.5*(dUi/dxj-dUj/dxi)*\mu>
        REAL(WP),ALLOCATABLE :: Skewness_FA(:,:)
        
        REAL(WP),ALLOCATABLE :: MKE_FA(:)
        REAL(WP),ALLOCATABLE :: TKE_FA(:)
        REAL(WP),ALLOCATABLE :: uffTKEffd_FA(:)
        REAL(WP),ALLOCATABLE :: uffMKEffd_FA(:)
        
        REAL(WP),ALLOCATABLE :: Tau_Mean_RA (:, :, :)       !<tau_ij>
        REAL(WP),ALLOCATABLE :: Tau_meaU_RA (:, :, :)       !<tau_ij>(<S>,<\mu>)
        REAL(WP),ALLOCATABLE :: dTaudy_RA   (:, :, :)       !d<tau_ij>/dy
        REAL(WP),ALLOCATABLE :: dTSSdy_RA   (:, :, :)       !d<R_ij>/dy
        REAL(WP),ALLOCATABLE :: TauU_RA     (:, :, :, :)    !< u_h* \tau_mn >
        REAL(WP),ALLOCATABLE :: Taufuf_RA   (:, :, :, :)    !< u'_h* \tau'_mn >
        REAL(WP),ALLOCATABLE :: TauDvDL_RA  (:, :, :, :, :) !< d(u_m)/d(x_n)* \tau_hp >
        !REAL(WP),ALLOCATABLE :: Tau_ik_Du_jDx_i_RA  (:, :, :) !< d(u_m)/d(x_n)* \tau_hp >
        
        
        REAL(WP),ALLOCATABLE :: ANISOTROPY_FA(:, :, :)
        REAL(WP),ALLOCATABLE :: ANISTPinva_FA(:, :)
        REAL(WP),ALLOCATABLE :: LumleyAxis_FA(:, :)
        
        !==================Checking==============================================
        REAL(WP)   :: NSFbalt_RA(NDV)
        REAL(WP)   :: NSFbalt_FA(NDV)
        REAL(WP),ALLOCATABLE   :: NSFbal_RA (:,:)
        REAL(WP),ALLOCATABLE   :: NSFbal_FA (:,:)
        REAL(WP)   :: BuoyForceTT
        REAL(WP),ALLOCATABLE   :: ENEbal_FA (:)
        REAL(WP)   :: ENEbalt_FA
        
        REAL(WP),ALLOCATABLE   :: bdfcintg (:)
        REAL(WP),ALLOCATABLE   :: densintg (:)
        
        !===============budgets================================
        !==================Ruv=============================
        REAL(WP),ALLOCATABLE :: BUDG_prodc_stres_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_dissp_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_pdudx_stran_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_Turbu_diffu_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_dpudx_diffu_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_diffu_duiuj(:, :)
        
        REAL(WP),ALLOCATABLE :: BUDG_press_accl1_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_accl1_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc1_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_balance1_duiuj(:, :)
        
        REAL(WP),ALLOCATABLE :: BUDG_prodc_gvfc2_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc2_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_turss_accl2_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_balance2_duiuj(:, :)
        
        REAL(WP),ALLOCATABLE :: BUDG_pressure3_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_vistress3_duiuj(:, :)
        REAL(WP),ALLOCATABLE :: BUDG_balance3_duiuj(:, :)
        
        !==================Ruv==Sum along y=======================
        REAL(WP),ALLOCATABLE :: BUDG_prodc_stres_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_dissp_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_pdudx_stran_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_Turbu_diffu_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_dpudx_diffu_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_diffu_duiuj_ysum(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_press_accl1_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_accl1_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc1_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance1_duiuj_ysum(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_prodc_gvfc2_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc2_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_turss_accl2_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance2_duiuj_ysum(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_pressure3_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_vistress3_duiuj_ysum(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance3_duiuj_ysum(:)
        
        !==================TKE=============================
        REAL(WP),ALLOCATABLE :: BUDG_prodc_stres_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_dissp_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_pdudx_stran_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_Turbu_diffu_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_dpudx_diffu_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_diffu_TKE(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_press_accl1_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_accl1_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc1_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance1_TKE(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_prodc_gvfc2_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc2_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_turss_accl2_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance2_TKE(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_pressure3_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_vistress3_TKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance3_TKE(:)
        
        !==================TKE sum along y=============================
        REAL(WP) :: BUDG_prodc_stres_TKE_ysum
        REAL(WP) :: BUDG_viscs_dissp_TKE_ysum
        REAL(WP) :: BUDG_pdudx_stran_TKE_ysum
        REAL(WP) :: BUDG_Turbu_diffu_TKE_ysum
        REAL(WP) :: BUDG_dpudx_diffu_TKE_ysum
        REAL(WP) :: BUDG_viscs_diffu_TKE_ysum
        
        REAL(WP) :: BUDG_press_accl1_TKE_ysum
        REAL(WP) :: BUDG_viscs_accl1_TKE_ysum
        REAL(WP) :: BUDG_prodc_dvfc1_TKE_ysum
        REAL(WP) :: BUDG_balance1_TKE_ysum
        
        REAL(WP) :: BUDG_prodc_gvfc2_TKE_ysum
        REAL(WP) :: BUDG_prodc_dvfc2_TKE_ysum
        REAL(WP) :: BUDG_turss_accl2_TKE_ysum
        REAL(WP) :: BUDG_balance2_TKE_ysum
        
        REAL(WP) :: BUDG_pressure3_TKE_ysum
        REAL(WP) :: BUDG_vistress3_TKE_ysum
        REAL(WP) :: BUDG_balance3_TKE_ysum
        
        !==================MKE=============================
        REAL(WP),ALLOCATABLE :: BUDG_prodc_stres_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_dissp_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_pdudx_stran_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_Turbu_diffu_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_dpudx_diffu_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_diffu_MKE(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_press_accl1_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_accl1_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc1_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_gvfc1_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance1_MKE(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_prodc_gvfc2_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_dvfc2_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_turss_accl2_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance2_MKE(:)
        
        REAL(WP),ALLOCATABLE :: BUDG_pressure3_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_vistress3_MKE(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance3_MKE(:)
        
        !==================MKE sum along y=============================
        REAL(WP) :: BUDG_prodc_stres_MKE_ysum
        REAL(WP) :: BUDG_viscs_dissp_MKE_ysum
        REAL(WP) :: BUDG_pdudx_stran_MKE_ysum
        REAL(WP) :: BUDG_Turbu_diffu_MKE_ysum
        REAL(WP) :: BUDG_dpudx_diffu_MKE_ysum
        REAL(WP) :: BUDG_viscs_diffu_MKE_ysum
        
        REAL(WP) :: BUDG_press_accl1_MKE_ysum
        REAL(WP) :: BUDG_viscs_accl1_MKE_ysum
        REAL(WP) :: BUDG_prodc_dvfc1_MKE_ysum
        REAL(WP) :: BUDG_prodc_gvfc1_MKE_ysum
        REAL(WP) :: BUDG_balance1_MKE_ysum
        
        REAL(WP) :: BUDG_prodc_gvfc2_MKE_ysum
        REAL(WP) :: BUDG_prodc_dvfc2_MKE_ysum
        REAL(WP) :: BUDG_turss_accl2_MKE_ysum
        REAL(WP) :: BUDG_balance2_MKE_ysum
        
        REAL(WP) :: BUDG_pressure3_MKE_ysum
        REAL(WP) :: BUDG_vistress3_MKE_ysum
        REAL(WP) :: BUDG_balance3_MKE_ysum
        

        !==============For RANS=================================
        REAL(WP),ALLOCATABLE :: RANS_Mut(:)
        
        !========================pp variables 2===============================
        REAL(WP),ALLOCATABLE :: H_FA(:)
        REAL(WP),ALLOCATABLE :: hff_RA(:)
        REAL(WP),ALLOCATABLE :: hfpf_RA(:)
        REAL(WP),ALLOCATABLE :: dTdX(:,:)
        REAL(WP),ALLOCATABLE :: dDdX(:,:)
        REAL(WP),ALLOCATABLE :: dHdX_RA(:,:)
        REAL(WP),ALLOCATABLE :: dHdX_FA(:,:)
        REAL(WP),ALLOCATABLE :: UH_FA(:,:)
        REAL(WP),ALLOCATABLE :: uff2hffd_FA(:,:,:)
        
        REAL(WP),ALLOCATABLE :: uffhffd_FA(:,:)
        REAL(WP),ALLOCATABLE :: ufhfd_RA(:,:)
        REAL(WP),ALLOCATABLE :: ViscStressEnth_RA(:,:,:)
        REAL(WP),ALLOCATABLE :: ViscStressEnthGrad_RA(:,:,:,:)
        
        REAL(WP),ALLOCATABLE :: BUDG_prodc_stres_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_enthg_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_Turbu_diffu_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_press_accl1_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_dphdx_diffu_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_pdhdx_stran_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_ConHF_accel_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_ConHF_diffu_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_ConHF_dissp_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_accl1_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_diffu_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_dissp_thf(:,:)
        REAL(WP),ALLOCATABLE :: BUDG_balance1_thf(:,:)
        
        REAL(WP),ALLOCATABLE :: BUDG_prodc_stres_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_prodc_enthg_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_Turbu_diffu_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_press_accl1_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_dphdx_diffu_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_pdhdx_stran_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_ConHF_accel_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_ConHF_diffu_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_ConHF_dissp_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_accl1_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_diffu_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_viscs_dissp_IEN(:)
        REAL(WP),ALLOCATABLE :: BUDG_balance1_IEN(:)
        
        !========================body force=================================
        
        REAL(WP),ALLOCATABLE :: BUDG_prodc_gvfc2_thf(:)
        
        !=========================thermal===================================
        INTEGER(4) :: J4SS0, J4Tpc, J4Tbk, J4MaxU, J4TbukSsd(2)
        
        !======on the wall ==================
        REAL(WP) :: Hwal_RA_d(2)
        REAL(WP) :: Twal_d(2)
        REAL(WP) :: Dwal_d(2)
        REAL(WP) :: Mwal_d(2)
        REAL(WP) :: Kwal_d(2)
        REAL(WP) :: Cpwal_d(2)
        REAL(WP) :: Bowal_d(2)
        REAL(WP) :: qw_d(2), qw_ave, qw_d_ave
        
        !=======coeffecieint================
        REAL(WP) :: Cf0_io(2),        Cf0_ave_io
        REAL(WP) :: Cfbk_io(2),       Cfbk_ave_io
        REAL(WP) :: CfbkSsd_io(2),    CfbkSsd_ave_io
        
        REAL(WP) :: Rebk, RebkSsd(2), RebkTsd(2)
        REAL(WP) :: Prbk, PrbkSsd(2), PrbkTsd(2)
        REAL(WP) :: Grbk(2), Grbk_drho, GrbkSsd(2)
        REAL(wp) :: Bobk(2)
        REAL(WP) :: Nubk(2), Nubk_dTw(2), NubkSsd(2), NubkTsd(2), Nu_int
        
        REAL(WP) :: hc_d(2),  hc_dTw_d(2), hc_dTw_d_ave, hcsd_d(2) 

        
        REAL(WP) :: L4TbkSsd(2), L4TbkTsd(2), L4Tbk(2)
        
        
        REAL(WP) :: RichardsonNo(2), Ttau(2), Ttau_d(2)
        REAL(WP) :: Mdot
        REAL(WP) :: Hdot
        REAL(WP) :: Gbuk, Gbuk_d, GbukSsd(2), GbukTsd(2)
        REAL(WP) :: Hbuk, Hbuk_d, HbukSsd(2), HbukTsd(2)
        REAL(WP) :: Ubuk, Ubuk_d, UbukSsd(2), UbukTsd(2)
        REAL(WP) :: Tbuk, Tbuk_d, TbukSsd(2), TbukTsd(2)
        REAL(WP) :: Mbuk, Mbuk_d, MbukSsd(2), MbukTsd(2)
        REAL(WP) :: Kbuk, Kbuk_d, KbukSsd(2), KbukTsd(2)
        REAL(WP) :: Cpbk, Cpbk_d, CpbkSsd(2), CpbkTsd(2)
        REAL(WP) :: Bbuk, Bbuk_d, BbukSsd(2), BbukTsd(2)
        REAL(WP) :: Dbuk, Dbuk_d, DbukSsd(2), DbukTsd(2)
        REAL(WP) :: D_int, D_int_d, Dsd_int(2)
        
        
        !REAL(WP),ALLOCATABLE :: Nuy(:,:)
        
        
        
    END MODULE     
!*****************************************************************************************************
    SUBROUTINE MEMO_ALLOCT_AVERAGE_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        !========================================
        ALLOCATE( U1xztL_F0_io( NCL2,NDV+1 ) )
        ALLOCATE( G1xztL_F0_io( NCL2,NDV   ) )
        ALLOCATE( UPxztL_F0_io( NCL2,NDV   ) )
        
        ALLOCATE( U2xztL_F0_io( NCL2,NDV*(7-NDV)/2+NDV-3 ) )
        ALLOCATE( UGxztL_F0_io( NCL2,NDV*(7-NDV)/2+NDV-3 ) )
        ALLOCATE( UGUxztL_F0_io(NCL2,NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8) )
        ALLOCATE( U3xztL_F0_io (NCL2,NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8) )
        
        ALLOCATE( DVDL1xztL_F0_io( NCL2, NDV, NDV  ) )
        ALLOCATE( DVDLPxztL_F0_io( NCL2, NDV, NDV  ) )
        ALLOCATE( DVDL2xztL_F0_io( NCL2, (NDV-1)*3+NDV, (NDV-1)*3+NDV  ) )
        
        ALLOCATE( QuadUVxztL_F0_io  ( NCL2, 4, QUADHN  ) )
        ALLOCATE( QuadVzxztL_F0_io  ( NCL2, 4, QUADHN  ) )
        ALLOCATE( QuadTKxztL_F0_io  ( NCL2, 4, QUADHN  ) )
        ALLOCATE( QuadDRxztL_F0_io  ( NCL2, 4, QUADHN  ) )
        ALLOCATE( QuadDUV1xztL_F0_io( NCL2, 4, QUADHN  ) )
        ALLOCATE( QuadDUV2xztL_F0_io( NCL2, 4, QUADHN  ) )
        
        ALLOCATE( OctDUVxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctDVzxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctDTKxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctDDRxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctDDUV1xztL_F0_io( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctDDUV2xztL_F0_io( NCL2, 8, QUADHN  ) )
        
        ALLOCATE( OctTUVxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctTVzxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctTTKxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctTDRxztL_F0_io  ( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctTDUV1xztL_F0_io( NCL2, 8, QUADHN  ) )
        ALLOCATE( OctTDUV2xztL_F0_io( NCL2, 8, QUADHN  ) )
         
        ALLOCATE( FUxztL_F0_io( NCL2,NDV+1 ) )
        
        !=============RA =========================
        ALLOCATE( dPdX_RA   (NCL2, NDV)    )
        ALLOCATE( ufpf_RA   (NCL2, NDV)    )
        ALLOCATE( uf2_RA    (NCL2, NDV, NDV)    )
        ALLOCATE( uf2d_RA   (NCL2, NDV, NDV)    )
        ALLOCATE( uf3_RA    (NCL2, NDV, NDV, NDV) )
        ALLOCATE( uf3d_RA   (NCL2, NDV, NDV, NDV) )
        ALLOCATE( UU_RA         (NCL2, NDV, NDV) )
        
        ALLOCATE( dUidXi        (NCL2)            )
        ALLOCATE( StrainTensor  (NCL2, NDV, NDV)  )
        ALLOCATE( VortcyTensor  (NCL2, NDV, NDV)  )
        ALLOCATE( Skewness_RA   (NCL2, NDV)       )
        
        ALLOCATE( MKE_RA        (NCL2)            )
        ALLOCATE( TKE_RA        (NCL2)            )
        ALLOCATE( ufTKEfd_RA    (NCL2)            )
        ALLOCATE( ufMKEfd_RA    (NCL2)            )
        
        
        ALLOCATE( Omega2_RA(NCL2,NDV) )
        ALLOCATE( Omega_RA2(NCL2,NDV) )
        ALLOCATE( Omega_rms(NCL2,NDV) )
        
        ALLOCATE( ANISOTROPY_RA(NCL2,NDV,NDV) )
        ALLOCATE( ANISTPinva_RA(NCL2,NDV) )
        ALLOCATE( LumleyAxis_RA(NCL2,2) )
        
        
        
        !=============FA =========================
        
        ALLOCATE( DrivenForce   (NCL2)           )
        ALLOCATE( U_FA          (NCL2, NDV)      )
        ALLOCATE( UU_FA         (NCL2, NDV, NDV) )
        ALLOCATE( dUdX_FA       (NCL2, NDV, NDV) )
        
        
        ALLOCATE( uff_RA        (NCL2, NDV)      )
        ALLOCATE( uff2_FA       (NCL2, NDV, NDV) )
        ALLOCATE( uff2d_FA      (NCL2, NDV, NDV) )
        ALLOCATE( uff3_FA       (NCL2, NDV, NDV, NDV) )
        ALLOCATE( uff3d_FA      (NCL2, NDV, NDV, NDV) )
        ALLOCATE( TDIFU_FA      (NCL2, NDV, NDV) )

        ALLOCATE( dUidXiM       (NCL2)            )
        ALLOCATE( StrainTensorM (NCL2, NDV, NDV)  )
        ALLOCATE( VortcyTensorM (NCL2, NDV, NDV)  )
        ALLOCATE( Skewness_FA   (NCL2, NDV)       )
        
        ALLOCATE( MKE_FA        (NCL2)            )
        ALLOCATE( TKE_FA        (NCL2)            )
        ALLOCATE( uffTKEffd_FA  (NCL2)            )
        ALLOCATE( uffMKEffd_FA  (NCL2)            )
        
        ALLOCATE( Tau_Mean_RA   (NCL2, NDV, NDV) )
        ALLOCATE( Tau_meaU_RA   (NCL2, NDV, NDV) )
        ALLOCATE( dTaudy_RA     (NCL2, NDV, NDV) )
        ALLOCATE( dTSSdy_RA     (NCL2, NDV, NDV) )
        ALLOCATE( TauU_RA       (NCL2, NDV, NDV, NDV) )
        ALLOCATE( Taufuf_RA     (NCL2, NDV, NDV, NDV     ) )
        ALLOCATE( TauDvDL_RA    (NCL2, NDV, NDV, NDV, NDV) )
        !ALLOCATE( Tau_ik_Du_jDx_i_RA(NCL2, NDV, NDV))
        
        ALLOCATE( ANISOTROPY_FA(NCL2,NDV,NDV) )
        ALLOCATE( ANISTPinva_FA(NCL2,NDV) )
        ALLOCATE( LumleyAxis_FA(NCL2,2) )
        
        ALLOCATE( NSFbal_FA(NCL2,NDV))
        ALLOCATE( NSFbal_RA(NCL2,NDV))
        ALLOCATE( ENEbal_FA(NCL2))
        
        ALLOCATE( bdfcintg(NCL2)) ; bdfcintg=0.0_WP
        ALLOCATE( densintg(NCL2)) ; densintg=0.0_WP
        
        
        !========================================
        !==================Ruv=============================
        ALLOCATE ( BUDG_prodc_stres_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_viscs_dissp_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_pdudx_stran_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_Turbu_diffu_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_dpudx_diffu_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_viscs_diffu_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        
        ALLOCATE ( BUDG_press_accl1_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_viscs_accl1_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_prodc_dvfc1_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_balance1_duiuj   (NCL2, NDV*(7-NDV)/2+NDV-3))
        
        ALLOCATE ( BUDG_prodc_gvfc2_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_prodc_dvfc2_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_turss_accl2_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_balance2_duiuj   (NCL2, NDV*(7-NDV)/2+NDV-3))
        
        ALLOCATE ( BUDG_pressure3_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_vistress3_duiuj(NCL2, NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_balance3_duiuj (NCL2, NDV*(7-NDV)/2+NDV-3))
        
        !==================Ruv==Sum along y=======================
        ALLOCATE ( BUDG_prodc_stres_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_viscs_dissp_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_pdudx_stran_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_Turbu_diffu_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_dpudx_diffu_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_viscs_diffu_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        
        ALLOCATE ( BUDG_press_accl1_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_viscs_accl1_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_prodc_dvfc1_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_balance1_duiuj_ysum   (NDV*(7-NDV)/2+NDV-3))
        
        ALLOCATE ( BUDG_prodc_gvfc2_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_prodc_dvfc2_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_turss_accl2_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_balance2_duiuj_ysum   (NDV*(7-NDV)/2+NDV-3))
        
        ALLOCATE ( BUDG_pressure3_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_vistress3_duiuj_ysum(NDV*(7-NDV)/2+NDV-3))
        ALLOCATE ( BUDG_balance3_duiuj_ysum (NDV*(7-NDV)/2+NDV-3))
        
        !==================TKE=============================
        ALLOCATE ( BUDG_prodc_stres_TKE(NCL2) )
        ALLOCATE ( BUDG_viscs_dissp_TKE(NCL2) )
        ALLOCATE ( BUDG_pdudx_stran_TKE(NCL2) )
        ALLOCATE ( BUDG_Turbu_diffu_TKE(NCL2) )
        ALLOCATE ( BUDG_dpudx_diffu_TKE(NCL2) )
        ALLOCATE ( BUDG_viscs_diffu_TKE(NCL2) )
        
        ALLOCATE ( BUDG_press_accl1_TKE(NCL2) )
        ALLOCATE ( BUDG_viscs_accl1_TKE(NCL2) )
        ALLOCATE ( BUDG_prodc_dvfc1_TKE(NCL2) )
        ALLOCATE ( BUDG_balance1_TKE   (NCL2) )
        
        ALLOCATE ( BUDG_prodc_gvfc2_TKE(NCL2) )
        ALLOCATE ( BUDG_prodc_dvfc2_TKE(NCL2) )
        ALLOCATE ( BUDG_turss_accl2_TKE(NCL2) )
        ALLOCATE ( BUDG_balance2_TKE   (NCL2) )
        
        ALLOCATE ( BUDG_pressure3_TKE(NCL2) )
        ALLOCATE ( BUDG_vistress3_TKE(NCL2) )
        ALLOCATE ( BUDG_balance3_TKE (NCL2) )
        
        !==================MKE=============================
        ALLOCATE ( BUDG_prodc_stres_MKE(NCL2) )
        ALLOCATE ( BUDG_viscs_dissp_MKE(NCL2) )
        ALLOCATE ( BUDG_pdudx_stran_MKE(NCL2) )
        ALLOCATE ( BUDG_Turbu_diffu_MKE(NCL2) )
        ALLOCATE ( BUDG_dpudx_diffu_MKE(NCL2) )
        ALLOCATE ( BUDG_viscs_diffu_MKE(NCL2) )
        
        ALLOCATE ( BUDG_press_accl1_MKE(NCL2) )
        ALLOCATE ( BUDG_viscs_accl1_MKE(NCL2) )
        ALLOCATE ( BUDG_prodc_dvfc1_MKE(NCL2) )
        ALLOCATE ( BUDG_prodc_gvfc1_MKE(NCL2) )
        ALLOCATE ( BUDG_balance1_MKE   (NCL2) )
        
        ALLOCATE ( BUDG_prodc_gvfc2_MKE(NCL2) )
        ALLOCATE ( BUDG_prodc_dvfc2_MKE(NCL2) )
        ALLOCATE ( BUDG_turss_accl2_MKE(NCL2) )
        ALLOCATE ( BUDG_balance2_MKE   (NCL2) )
        
        ALLOCATE ( BUDG_pressure3_MKE(NCL2) )
        ALLOCATE ( BUDG_vistress3_MKE(NCL2) )
        ALLOCATE ( BUDG_balance3_MKE (NCL2) )
        
        
        !===========================================
        ALLOCATE ( RANS_Mut(NCL2) )
        !========================================
        ALLOCATE( TauwSD(NCL2)         )
        ALLOCATE( DensSD(NCL2)         )
        ALLOCATE( ViscSD(NCL2)         )
        ALLOCATE( YWdiSD(NCL2)         )

        ALLOCATE( D1xztL_F0_io( NCL2 ) ) ;  D1xztL_F0_io = 1.0_WP
        ALLOCATE( M1xztL_F0_io( NCL2 ) ) ;  M1xztL_F0_io = 1.0_WP
        
        ALLOCATE( DVDL1MxztL_F0_io ( NCL2, NDV, NDV  )                     )
        ALLOCATE( DVDL1MUxztL_F0_io( NCL2, NDV, NDV, NDV  )                )
        ALLOCATE( DVDL2MxztL_F0_io ( NCL2, (NDV-1)*3+NDV, (NDV-1)*3+NDV  ) )
        
        IF(thermlflg ==1) THEN
            !========================================
            ALLOCATE( T1xztL_F0_io( NCL2 ) )
            ALLOCATE( H1xztL_F0_io( NCL2 ) )
            
            ALLOCATE( DVDL1MHxztL_F0_io( NCL2, NDV, NDV  )                     )
            ALLOCATE( T2xztL_F0_io( NCL2 ) )
            ALLOCATE( D2xztL_F0_io( NCL2 ) )
            ALLOCATE( H2xztL_F0_io( NCL2 ) )
            
            ALLOCATE( DHxztL_F0_io( NCL2 ) )
            ALLOCATE( PHxztL_F0_io( NCL2 ) )
            
!            ALLOCATE( DVDL1MxztL_F0_io ( NCL2, NDV, NDV  )                     )
!            ALLOCATE( DVDL1MHxztL_F0_io( NCL2, NDV, NDV  )                     )
!            ALLOCATE( DVDL1MUxztL_F0_io( NCL2, NDV, NDV, NDV  )                )
!            ALLOCATE( DVDL2MxztL_F0_io ( NCL2, (NDV-1)*3+NDV, (NDV-1)*3+NDV  ) )
            
            ALLOCATE( UHxztL_F0_io  ( NCL2,NDV )                 )
            ALLOCATE( GHxztL_F0_io  ( NCL2,NDV )                 )
            ALLOCATE( U2DHxztL_F0_io( NCL2,NDV*(7-NDV)/2+NDV-3 ) )
            
            ALLOCATE( DhDL1xztL_F0_io ( NCL2,NDV )     )
            ALLOCATE( DhDLPxztL_F0_io ( NCL2,NDV )     )
            ALLOCATE( DTDLKxztL_F0_io ( NCL2,NDV )     )
            ALLOCATE( DTDLKUxztL_F0_io(NCL2,NDV, NDV ) )
            
            ALLOCATE( DTDLKDVDLxztL_F0_io(NCL2,NDV, NDV, NDV ) )
            ALLOCATE( DHDLMDVDLxztL_F0_io(NCL2,NDV, NDV, NDV ) )
            
            
            ALLOCATE( CpSD  (NCL2)         )
            ALLOCATE( QwSD  (NCL2)         )
            ALLOCATE( TwSD  (NCL2)         )
            ALLOCATE( HwSD  (NCL2)         )
        
            !============================
            !ALLOCATE( Nuy(NCL2,2))
        
            !========================================
            ALLOCATE( H_FA      (0:NND2)         )
            ALLOCATE( hff_RA    (NCL2)         )
            ALLOCATE( hfpf_RA (NCL2)         )
            ALLOCATE( dTdX      (NCL2,NDV)     )
            ALLOCATE( dDdX      (NCL2,NDV)     )
            ALLOCATE( dHdX_RA   (NCL2,NDV)     )
            ALLOCATE( dHdX_FA   (NCL2,NDV)     )
            ALLOCATE( UH_FA     (NCL2,NDV)     )
            ALLOCATE( uff2hffd_FA(NCL2,NDV,NDV) )
           
            ALLOCATE( uffhffd_FA        (NCL2, NDV)      )
            ALLOCATE( ufhfd_RA        (NCL2, NDV)      )
            ALLOCATE( ViscStressEnth_RA    (NCL2, NDV, NDV) )
            ALLOCATE( ViscStressEnthGrad_RA(NCL2, NDV, NDV, NDV) )
            !========================================
            
            ALLOCATE( BUDG_prodc_stres_thf(NCL2,NDV) ) ; BUDG_prodc_stres_thf=0.0_WP
            ALLOCATE( BUDG_prodc_enthg_thf(NCL2,NDV) ) ; BUDG_prodc_enthg_thf=0.0_WP
            ALLOCATE( BUDG_Turbu_diffu_thf(NCL2,NDV) ) ; BUDG_Turbu_diffu_thf=0.0_WP
            ALLOCATE( BUDG_press_accl1_thf(NCL2,NDV) ) ; BUDG_press_accl1_thf=0.0_WP
            ALLOCATE( BUDG_dphdx_diffu_thf(NCL2,NDV) ) ; BUDG_dphdx_diffu_thf=0.0_WP
            ALLOCATE( BUDG_pdhdx_stran_thf(NCL2,NDV) ) ; BUDG_pdhdx_stran_thf=0.0_WP
            ALLOCATE( BUDG_ConHF_accel_thf(NCL2,NDV) ) ; BUDG_ConHF_accel_thf=0.0_WP
            ALLOCATE( BUDG_ConHF_diffu_thf(NCL2,NDV) ) ; BUDG_ConHF_diffu_thf=0.0_WP
            ALLOCATE( BUDG_ConHF_dissp_thf(NCL2,NDV) ) ; BUDG_ConHF_dissp_thf=0.0_WP
            ALLOCATE( BUDG_viscs_accl1_thf(NCL2,NDV) ) ; BUDG_viscs_accl1_thf=0.0_WP
            ALLOCATE( BUDG_viscs_diffu_thf(NCL2,NDV) ) ; BUDG_viscs_diffu_thf=0.0_WP
            ALLOCATE( BUDG_viscs_dissp_thf(NCL2,NDV) ) ; BUDG_viscs_dissp_thf=0.0_WP
            ALLOCATE( BUDG_balance1_thf    (NCL2,NDV) ) ; BUDG_balance1_thf    =0.0_WP
            !========================body force=================================
            ALLOCATE( BUDG_prodc_gvfc2_thf(NCL2)     ) ; BUDG_prodc_gvfc2_thf=0.0_WP
            
            
            ALLOCATE( BUDG_prodc_stres_IEN(NCL2) ) ; BUDG_prodc_stres_IEN=0.0_WP
            ALLOCATE( BUDG_prodc_enthg_IEN(NCL2) ) ; BUDG_prodc_enthg_IEN=0.0_WP
            ALLOCATE( BUDG_Turbu_diffu_IEN(NCL2) ) ; BUDG_Turbu_diffu_IEN=0.0_WP
            ALLOCATE( BUDG_press_accl1_IEN(NCL2) ) ; BUDG_press_accl1_IEN=0.0_WP
            ALLOCATE( BUDG_dphdx_diffu_IEN(NCL2) ) ; BUDG_dphdx_diffu_IEN=0.0_WP
            ALLOCATE( BUDG_pdhdx_stran_IEN(NCL2) ) ; BUDG_pdhdx_stran_IEN=0.0_WP
            ALLOCATE( BUDG_ConHF_accel_IEN(NCL2) ) ; BUDG_ConHF_accel_IEN=0.0_WP
            ALLOCATE( BUDG_ConHF_diffu_IEN(NCL2) ) ; BUDG_ConHF_diffu_IEN=0.0_WP
            ALLOCATE( BUDG_ConHF_dissp_IEN(NCL2) ) ; BUDG_ConHF_dissp_IEN=0.0_WP
            ALLOCATE( BUDG_viscs_accl1_IEN(NCL2) ) ; BUDG_viscs_accl1_IEN=0.0_WP
            ALLOCATE( BUDG_viscs_diffu_IEN(NCL2) ) ; BUDG_viscs_diffu_IEN=0.0_WP
            ALLOCATE( BUDG_viscs_dissp_IEN(NCL2) ) ; BUDG_viscs_dissp_IEN=0.0_WP
            ALLOCATE( BUDG_balance1_IEN    (NCL2) ) ; BUDG_balance1_IEN    =0.0_WP
            
        END IF
    
        RETURN
    END SUBROUTINE 

    
    
!*****************************************************************************************************
    SUBROUTINE MEMO_DEALLT_AVERAGE_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        
        !========================================
        DEALLOCATE( U1xztL_F0_io )
        DEALLOCATE( G1xztL_F0_io )
        DEALLOCATE( UPxztL_F0_io )
        
        DEALLOCATE( U2xztL_F0_io  )
        DEALLOCATE( UGxztL_F0_io  )
        DEALLOCATE( UGUxztL_F0_io )
        DEALLOCATE( U3xztL_F0_io )
        
        DEALLOCATE( DVDL1xztL_F0_io )
        DEALLOCATE( DVDLPxztL_F0_io )
        DEALLOCATE( DVDL2xztL_F0_io )
         
        DEALLOCATE( QuadUVxztL_F0_io )
        DEALLOCATE( QuadVzxztL_F0_io )
        DEALLOCATE( QuadTKxztL_F0_io )
        DEALLOCATE( QuadDRxztL_F0_io )
        DEALLOCATE( QuadDUV1xztL_F0_io )
        DEALLOCATE( QuadDUV2xztL_F0_io )
        
        DEALLOCATE( OctTUVxztL_F0_io )
        DEALLOCATE( OctTVzxztL_F0_io )
        DEALLOCATE( OctTTKxztL_F0_io )
        DEALLOCATE( OctTDRxztL_F0_io )
        DEALLOCATE( OctTDUV1xztL_F0_io )
        DEALLOCATE( OctTDUV2xztL_F0_io )
        
        DEALLOCATE( OctDUVxztL_F0_io )
        DEALLOCATE( OctDVzxztL_F0_io )
        DEALLOCATE( OctDTKxztL_F0_io )
        DEALLOCATE( OctDDRxztL_F0_io )
        DEALLOCATE( OctDDUV1xztL_F0_io )
        DEALLOCATE( OctDDUV2xztL_F0_io )
        
        DEALLOCATE( FUxztL_F0_io )
        !========================================
        DEALLOCATE( dPdX_RA    )
        DEALLOCATE( ufpf_RA    )
        DEALLOCATE( uf2_RA     )
        DEALLOCATE( uf2d_RA    )
        DEALLOCATE( uf3_RA     )
        DEALLOCATE( uf3d_RA    )
        DEALLOCATE( UU_RA      )
        
        DEALLOCATE( dUidXi     )
        DEALLOCATE( StrainTensor)
        DEALLOCATE( VortcyTensor)
        DEALLOCATE( Skewness_RA)
        
        DEALLOCATE( MKE_RA     )
        DEALLOCATE( TKE_RA     )
        DEALLOCATE( ufTKEfd_RA )
        DEALLOCATE( ufMKEfd_RA )
        
        DEALLOCATE( Omega2_RA              )
        DEALLOCATE( Omega_RA2              )
        DEALLOCATE( Omega_rms              )
        
        DEALLOCATE( ANISOTROPY_RA          )
        DEALLOCATE( ANISTPinva_RA          )
        DEALLOCATE( LumleyAxis_RA          )
        
        !========================================
        DEALLOCATE( DrivenForce)
        DEALLOCATE( U_FA       )
        DEALLOCATE( UU_FA      )
        DEALLOCATE( dUdX_FA    )
        
        DEALLOCATE( uff_RA     )
        DEALLOCATE( uff2_FA    )
        DEALLOCATE( uff2d_FA   )
        DEALLOCATE( uff3_FA    ) 
        DEALLOCATE( uff3d_FA   )
        DEALLOCATE( TDIFU_FA   )
         
        DEALLOCATE( dUidXiM      )
        DEALLOCATE( StrainTensorM)
        DEALLOCATE( VortcyTensorM)
        DEALLOCATE( Skewness_FA  )
        
        DEALLOCATE( MKE_FA       )
        DEALLOCATE( TKE_FA       )
        DEALLOCATE( uffTKEffd_FA )
        DEALLOCATE( uffMKEffd_FA )
        
        DEALLOCATE( Tau_Mean_RA  )
        DEALLOCATE( Tau_meaU_RA  )
        DEALLOCATE( dTaudy_RA    )
        DEALLOCATE( dTSSdy_RA    )
        DEALLOCATE( TauU_RA      )
        DEALLOCATE( Taufuf_RA    )
        DEALLOCATE( TauDvDL_RA   )
        !DEALLOCATE( Tau_ik_Du_jDx_i_RA)
        
        DEALLOCATE( ANISOTROPY_FA          )
        DEALLOCATE( ANISTPinva_FA          )
        DEALLOCATE( LumleyAxis_FA          )
        
        DEALLOCATE( NSFbal_FA)
        DEALLOCATE( NSFbal_RA)
        DEALLOCATE( ENEbal_FA)
        
        DEALLOCATE( bdfcintg)
        DEALLOCATE( densintg)
        !========================================
        !==================Ruv=============================
        DEALLOCATE (  BUDG_prodc_stres_duiuj )
        DEALLOCATE (  BUDG_viscs_dissp_duiuj )
        DEALLOCATE (  BUDG_pdudx_stran_duiuj )
        DEALLOCATE (  BUDG_Turbu_diffu_duiuj )
        DEALLOCATE (  BUDG_dpudx_diffu_duiuj )
        DEALLOCATE (  BUDG_viscs_diffu_duiuj )
        
        DEALLOCATE (  BUDG_press_accl1_duiuj )
        DEALLOCATE (  BUDG_viscs_accl1_duiuj )
        DEALLOCATE (  BUDG_prodc_dvfc1_duiuj )
        DEALLOCATE (  BUDG_balance1_duiuj    )
        
        DEALLOCATE (  BUDG_prodc_gvfc2_duiuj )
        DEALLOCATE (  BUDG_prodc_dvfc2_duiuj )
        DEALLOCATE (  BUDG_turss_accl2_duiuj )
        DEALLOCATE (  BUDG_balance2_duiuj    )
        
        DEALLOCATE (  BUDG_pressure3_duiuj )
        DEALLOCATE (  BUDG_vistress3_duiuj )
        DEALLOCATE (  BUDG_balance3_duiuj  )
        
        !==================Ruv==Sum along y=======================
        DEALLOCATE (  BUDG_prodc_stres_duiuj_ysum )
        DEALLOCATE (  BUDG_viscs_dissp_duiuj_ysum )
        DEALLOCATE (  BUDG_pdudx_stran_duiuj_ysum )
        DEALLOCATE (  BUDG_Turbu_diffu_duiuj_ysum )
        DEALLOCATE (  BUDG_dpudx_diffu_duiuj_ysum )
        DEALLOCATE (  BUDG_viscs_diffu_duiuj_ysum )
        
        DEALLOCATE (  BUDG_press_accl1_duiuj_ysum )
        DEALLOCATE (  BUDG_viscs_accl1_duiuj_ysum )
        DEALLOCATE (  BUDG_prodc_dvfc1_duiuj_ysum )
        DEALLOCATE (  BUDG_balance1_duiuj_ysum    )
        
        DEALLOCATE (  BUDG_prodc_gvfc2_duiuj_ysum )
        DEALLOCATE (  BUDG_prodc_dvfc2_duiuj_ysum )
        DEALLOCATE (  BUDG_turss_accl2_duiuj_ysum )
        DEALLOCATE (  BUDG_balance2_duiuj_ysum    )
        
        DEALLOCATE (  BUDG_pressure3_duiuj_ysum )
        DEALLOCATE (  BUDG_vistress3_duiuj_ysum )
        DEALLOCATE (  BUDG_balance3_duiuj_ysum  )
        
        !==================TKE=============================
        DEALLOCATE (  BUDG_prodc_stres_TKE )
        DEALLOCATE (  BUDG_viscs_dissp_TKE )
        DEALLOCATE (  BUDG_pdudx_stran_TKE )
        DEALLOCATE (  BUDG_Turbu_diffu_TKE )
        DEALLOCATE (  BUDG_dpudx_diffu_TKE )
        DEALLOCATE (  BUDG_viscs_diffu_TKE )
        
        DEALLOCATE (  BUDG_press_accl1_TKE )
        DEALLOCATE (  BUDG_viscs_accl1_TKE )
        DEALLOCATE (  BUDG_prodc_dvfc1_TKE )
        DEALLOCATE (  BUDG_balance1_TKE    )
        
        DEALLOCATE (  BUDG_prodc_gvfc2_TKE )
        DEALLOCATE (  BUDG_prodc_dvfc2_TKE )
        DEALLOCATE (  BUDG_turss_accl2_TKE )
        DEALLOCATE (  BUDG_balance2_TKE    )
        
        DEALLOCATE (  BUDG_pressure3_TKE )
        DEALLOCATE (  BUDG_vistress3_TKE )
        DEALLOCATE (  BUDG_balance3_TKE  )
        
        !==================MKE=============================
        DEALLOCATE (  BUDG_prodc_stres_MKE )
        DEALLOCATE (  BUDG_viscs_dissp_MKE )
        DEALLOCATE (  BUDG_pdudx_stran_MKE )
        DEALLOCATE (  BUDG_Turbu_diffu_MKE )
        DEALLOCATE (  BUDG_dpudx_diffu_MKE )
        DEALLOCATE (  BUDG_viscs_diffu_MKE )
        
        DEALLOCATE (  BUDG_press_accl1_MKE )
        DEALLOCATE (  BUDG_viscs_accl1_MKE )
        DEALLOCATE (  BUDG_prodc_dvfc1_MKE )
        DEALLOCATE (  BUDG_prodc_gvfc1_MKE )
        DEALLOCATE (  BUDG_balance1_MKE    )
        
        DEALLOCATE (  BUDG_prodc_gvfc2_MKE )
        DEALLOCATE (  BUDG_prodc_dvfc2_MKE )
        DEALLOCATE (  BUDG_turss_accl2_MKE )
        DEALLOCATE (  BUDG_balance2_MKE    )
        
        DEALLOCATE (  BUDG_pressure3_MKE )
        DEALLOCATE (  BUDG_vistress3_MKE )
        DEALLOCATE (  BUDG_balance3_MKE  )

        
        !======================================
        DEALLOCATE( RANS_Mut)
        
        DEALLOCATE( YWdiSD )
        DEALLOCATE( TauwSD )
        DEALLOCATE( DensSD )
        DEALLOCATE( ViscSD )
            

        DEALLOCATE( D1xztL_F0_io ) 
        DEALLOCATE( M1xztL_F0_io ) 
        
        DEALLOCATE( DVDL1MxztL_F0_io  )
        DEALLOCATE( DVDL1MUxztL_F0_io )
        DEALLOCATE( DVDL2MxztL_F0_io  )
        
        IF(thermlflg ==1) THEN
            !========================================
            DEALLOCATE( T1xztL_F0_io )
            DEALLOCATE( H1xztL_F0_io )
            
            DEALLOCATE( DVDL1MHxztL_F0_io )
            DEALLOCATE( T2xztL_F0_io )
            DEALLOCATE( D2xztL_F0_io )
            DEALLOCATE( H2xztL_F0_io )
            
            DEALLOCATE( DHxztL_F0_io )
            DEALLOCATE( PHxztL_F0_io )
            
            DEALLOCATE( UHxztL_F0_io   )
            DEALLOCATE( GHxztL_F0_io   )
            DEALLOCATE( U2DHxztL_F0_io )
            
            DEALLOCATE( DhDL1xztL_F0_io      )
            DEALLOCATE( DhDLPxztL_F0_io      )
            DEALLOCATE( DTDLKxztL_F0_io      )
            DEALLOCATE( DTDLKUxztL_F0_io     )
            
            DEALLOCATE( DTDLKDVDLxztL_F0_io )
            DEALLOCATE( DHDLMDVDLxztL_F0_io )
            
            !==============================
            !DEALLOCATE( Nuy  )
            !==================
            
            DEALLOCATE( CpSD )
            DEALLOCATE( QwSD )
            DEALLOCATE( TwSD )
            DEALLOCATE( HwSD )
            
            !========================================
            DEALLOCATE( H_FA       )
            DEALLOCATE( hff_RA     )
            DEALLOCATE( hfpf_RA  )
            DEALLOCATE( dTdX       )
            DEALLOCATE( dDdX       )
            DEALLOCATE( dHdX_RA    )
            DEALLOCATE( dHdX_FA    )
            DEALLOCATE( UH_FA      )
            DEALLOCATE( uff2hffd_FA )
           
            DEALLOCATE( uffhffd_FA         )
            DEALLOCATE( ufhfd_RA         )
            DEALLOCATE( ViscStressEnth_RA     )
            DEALLOCATE( ViscStressEnthGrad_RA )
            !========================================
            
            DEALLOCATE( BUDG_prodc_stres_thf )
            DEALLOCATE( BUDG_prodc_enthg_thf )
            DEALLOCATE( BUDG_Turbu_diffu_thf )
            DEALLOCATE( BUDG_press_accl1_thf )
            DEALLOCATE( BUDG_dphdx_diffu_thf )
            DEALLOCATE( BUDG_pdhdx_stran_thf )
            DEALLOCATE( BUDG_ConHF_accel_thf )
            DEALLOCATE( BUDG_ConHF_diffu_thf )
            DEALLOCATE( BUDG_ConHF_dissp_thf )
            DEALLOCATE( BUDG_viscs_accl1_thf )
            DEALLOCATE( BUDG_viscs_diffu_thf )
            DEALLOCATE( BUDG_viscs_dissp_thf )
            DEALLOCATE( BUDG_balance1_thf    )
            
            DEALLOCATE( BUDG_prodc_stres_IEN )
            DEALLOCATE( BUDG_prodc_enthg_IEN )
            DEALLOCATE( BUDG_Turbu_diffu_IEN )
            DEALLOCATE( BUDG_press_accl1_IEN )
            DEALLOCATE( BUDG_dphdx_diffu_IEN )
            DEALLOCATE( BUDG_pdhdx_stran_IEN )
            DEALLOCATE( BUDG_ConHF_accel_IEN )
            DEALLOCATE( BUDG_ConHF_diffu_IEN )
            DEALLOCATE( BUDG_ConHF_dissp_IEN )
            DEALLOCATE( BUDG_viscs_accl1_IEN )
            DEALLOCATE( BUDG_viscs_diffu_IEN )
            DEALLOCATE( BUDG_viscs_dissp_IEN )
            DEALLOCATE( BUDG_balance1_IEN    )
            !========================body force=================================
            
            DEALLOCATE( BUDG_prodc_gvfc2_thf   )
    
        END IF

        RETURN
    END SUBROUTINE
    

!*****************************************************************************************************
    SUBROUTINE WRT_AVERAGE_PPED_Xperiodic_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        IF(MYID.EQ.0) CALL CHKRLHDL('15.IO: postprocessing general data at',myid,phyTIME_io)
         
        !==========Gather Dwta================
        CALL WRT_AVERAGE_PPED_XZ_IO_GATHER
        
        !==========Calcuate wall info================
        CALL PP_wall_thermal_shear(flgxzt)
        
        IF(MYID.NE.0) RETURN
        
        !=========Calculate data=======================
        CALL PP_FLOW_BASIC_VARS_XZ_IO
        CALL PP_FLOW_FA_RSTE_BUDG_XZ_IO
        CALL PP_SSzero_SIDED(1)
        !==========Write out wall info=================
        IF(thermlflg==1) THEN
            CALL WRT_HeatTransfer_Table_XZ_IO
        END IF
        CALL WRT_Cf_Table_XZ_IO ! must be after heat 
        !CALL PP_Umax_SIDED                          
        CALL WRT_FLOW_FA_Profile_XZ_IO
        IF(thermlflg==1) THEN
            CALL PP_HEAT_BASIC_VARS_XZ_IO
            CALL PP_HEAT_FA_RSTE_BUDG_XZ_IO
            CALL WRT_HEAT_FA_Profile_XZ_IO
        END IF
        
        CALL PP_FLOW_RA_noDen_RSTE_BUDG_XZ_IO
        if(thermlflg .ne. 1) CALL PP_SSzero_SIDED(2)
        CALL WRT_FLOW_RA_Profile_XZ_IO
        
        CALL WRT_Checking_TABLE_XZ_IO
        
        IF(ppspectra==1) THEN
            CALL WRITE_SPECO_AVE_PROFILE('FLOW')
            CALL WRITE_SPECO_AVE_Contour('FLOW')
            IF(thermlflg==1) THEN
            CALL WRITE_SPECO_AVE_PROFILE('HEAT')
            CALL WRITE_SPECO_AVE_Contour('HEAT')
            END IF
        END IF
        !CALL MEMO_DEALLT_INTP_XZ_IO
        CALL MEMO_DEALLT_AVERAGE_XZ_IO
        
        IF(MYID.EQ.0) CALL CHKRLHDL('   IO: finished postprocessing general data at',myid,phyTIME_io)
    
    END SUBROUTINE

!*****************************************************************************************************
    SUBROUTINE WRT_AVERAGE_PPED_XZ_IO_GATHER
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4)  :: INN
        INTEGER(4)  :: J, JJ, I
        INTEGER(4)  :: L, IP, M, H, N
        INTEGER(4)  :: N2DOID
        
        CHARACTER(15) :: NXSTR
        REAL(WP)   :: urms, vrms, wrms, uv, uw, vw
        REAL(WP)   ::  COE, Pwall
        INTEGER(4) :: TECFLG1, TECFLG2, TECFLG3
        
        REAL(WP)     :: FU_AUX (N2DO(MYID), NDV+1,                               1:NPTOT)
        REAL(WP)     :: U1_AUX (N2DO(MYID), NDV+1,                               1:NPTOT)
        REAL(WP)     :: G1_AUX (N2DO(MYID), NDV,                                 1:NPTOT)
        REAL(WP)     :: UP_AUX (N2DO(MYID), NDV,                                 1:NPTOT)
        
        REAL(WP)     :: U2_AUX (N2DO(MYID),(NDV*(7-NDV)/2+NDV-3),                1:NPTOT)
        REAL(WP)     :: UG_AUX (N2DO(MYID),(NDV*(7-NDV)/2+NDV-3),                1:NPTOT)
        REAL(WP)     :: UGU_AUX(N2DO(MYID),(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8),  1:NPTOT)
        REAL(WP)     :: U3_AUX (N2DO(MYID),(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8),  1:NPTOT)
        
        REAL(WP)     :: DVDL1_AUX (N2DO(MYID), NDV,                 NDV, 1:NPTOT)
        REAL(WP)     :: DVDLP_AUX (N2DO(MYID), NDV,                 NDV, 1:NPTOT)
        REAL(WP)     :: DVDL2_AUX (N2DO(MYID),(NDV-1)*3+NDV, (NDV-1)*3+NDV, 1:NPTOT)
        
        REAL(WP)     :: QuadUV_AUX (N2DO(MYID),   4,QUADHN, 1:NPTOT)
        REAL(WP)     :: QuadVz_AUX (N2DO(MYID),   4,QUADHN, 1:NPTOT)
        REAL(WP)     :: QuadTK_AUX (N2DO(MYID),   4,QUADHN, 1:NPTOT)
        REAL(WP)     :: QuadDR_AUX (N2DO(MYID),   4,QUADHN, 1:NPTOT)
        REAL(WP)     :: QuadDUV1_AUX (N2DO(MYID), 4,QUADHN, 1:NPTOT)
        REAL(WP)     :: QuadDUV2_AUX (N2DO(MYID), 4,QUADHN, 1:NPTOT)
        
        REAL(WP)     :: OctDUV_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctDVz_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctDTK_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctDDR_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctDDUV1_AUX (N2DO(MYID), 8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctDDUV2_AUX (N2DO(MYID), 8,QUADHN, 1:NPTOT)
        
        REAL(WP)     :: OctTUV_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctTVz_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctTTK_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctTDR_AUX (N2DO(MYID),   8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctTDUV1_AUX (N2DO(MYID), 8,QUADHN, 1:NPTOT)
        REAL(WP)     :: OctTDUV2_AUX (N2DO(MYID), 8,QUADHN, 1:NPTOT)
        
        REAL(WP)     :: T1_AUX (N2DO(MYID), 1:NPTOT)
        REAL(WP)     :: D1_AUX (N2DO(MYID), 1:NPTOT)
        REAL(WP)     :: H1_AUX (N2DO(MYID), 1:NPTOT)
        REAL(WP)     :: M1_AUX (N2DO(MYID), 1:NPTOT)
        
        REAL(WP)     :: T2_AUX (N2DO(MYID), 1:NPTOT)
        REAL(WP)     :: D2_AUX (N2DO(MYID), 1:NPTOT)
        REAL(WP)     :: H2_AUX (N2DO(MYID), 1:NPTOT)
        
        REAL(WP)     :: DH_AUX (N2DO(MYID), 1:NPTOT)
        REAL(WP)     :: PH_AUX (N2DO(MYID), 1:NPTOT)
        
        REAL(WP)     :: DVDL1M_AUX  (N2DO(MYID), NDV,                 NDV,     1:NPTOT)
        REAL(WP)     :: DVDL1MH_AUX (N2DO(MYID), NDV,                 NDV,     1:NPTOT)
        REAL(WP)     :: DVDL1MU_AUX (N2DO(MYID), NDV, NDV, NDV,                1:NPTOT)
        REAL(WP)     :: DVDL2M_AUX  (N2DO(MYID), (NDV-1)*3+NDV, (NDV-1)*3+NDV, 1:NPTOT)
        
        REAL(WP)     :: UH_AUX   (N2DO(MYID), NDV, 1:NPTOT)
        REAL(WP)     :: GH_AUX   (N2DO(MYID), NDV, 1:NPTOT)
        REAL(WP)     :: U2DH_AUX (N2DO(MYID), (NDV*(7-NDV)/2+NDV-3), 1:NPTOT)
        
        REAL(WP)     :: DhDL1_AUX (N2DO(MYID), NDV, 1:NPTOT)
        REAL(WP)     :: DhDLP_AUX (N2DO(MYID), NDV, 1:NPTOT)
        REAL(WP)     :: DTDLK_AUX     (N2DO(MYID), NDV, 1:NPTOT)
        REAL(WP)     :: DTDLKU_AUX    (N2DO(MYID), NDV, NDV, 1:NPTOT)
        REAL(WP)     :: DTDLKDVDL_AUX (N2DO(MYID), NDV, NDV, NDV, 1:NPTOT)
        REAL(WP)     :: DHDLMDVDL_AUX (N2DO(MYID), NDV, NDV, NDV, 1:NPTOT)
        
        !======================================
        INN = N2DO(MYID)*(NDV+1)
        CALL MPI_GATHER( U1xztL_io, INN, MPI_DOUBLE_PRECISION, U1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*NDV
        CALL MPI_GATHER( G1xztL_io, INN, MPI_DOUBLE_PRECISION, G1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)  
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*NDV
        CALL MPI_GATHER( UPxztL_io, INN, MPI_DOUBLE_PRECISION, UP_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !======================================
        INN = N2DO(MYID)*(NDV*(7-NDV)/2+NDV-3)
        CALL MPI_GATHER( U2xztL_io, INN, MPI_DOUBLE_PRECISION, U2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*(NDV*(7-NDV)/2+NDV-3)
        CALL MPI_GATHER( UGxztL_io, INN, MPI_DOUBLE_PRECISION, UG_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
        CALL MPI_GATHER( UGUxztL_io,INN, MPI_DOUBLE_PRECISION, UGU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
        CALL MPI_GATHER( U3xztL_io,INN, MPI_DOUBLE_PRECISION, U3_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !======================================
        INN = N2DO(MYID)*NDV*NDV
        CALL MPI_GATHER( DVDL1xztL_io, INN, MPI_DOUBLE_PRECISION, DVDL1_AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*NDV*NDV
        CALL MPI_GATHER( DVDLPxztL_io, INN, MPI_DOUBLE_PRECISION, DVDLP_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*((NDV-1)*3+NDV)*((NDV-1)*3+NDV)
        CALL MPI_GATHER( DVDL2xztL_io, INN, MPI_DOUBLE_PRECISION, DVDL2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)

        INN = N2DO(MYID)*4*QUADHN
        CALL MPI_GATHER( QuadUVxztL_io, INN, MPI_DOUBLE_PRECISION, QuadUV_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*4*QUADHN
        CALL MPI_GATHER( QuadVzxztL_io, INN, MPI_DOUBLE_PRECISION, QuadVz_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*4*QUADHN
        CALL MPI_GATHER( QuadTKxztL_io, INN, MPI_DOUBLE_PRECISION, QuadTK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*4*QUADHN
        CALL MPI_GATHER( QuadDRxztL_io, INN, MPI_DOUBLE_PRECISION, QuadDR_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*4*QUADHN
        CALL MPI_GATHER( QuadDUV1xztL_io, INN, MPI_DOUBLE_PRECISION, QuadDUV1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*4*QUADHN
        CALL MPI_GATHER( QuadDUV2xztL_io, INN, MPI_DOUBLE_PRECISION, QuadDUV2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !=======================================================
        INN = N2DO(MYID)*(NDV+1)
        CALL MPI_GATHER( FUxztL_io, INN, MPI_DOUBLE_PRECISION, FU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !=======================================================
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctDUVxztL_io, INN, MPI_DOUBLE_PRECISION, OctDUV_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctDVzxztL_io, INN, MPI_DOUBLE_PRECISION, OctDVz_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctDTKxztL_io, INN, MPI_DOUBLE_PRECISION, OctDTK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctDDRxztL_io, INN, MPI_DOUBLE_PRECISION, OctDDR_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctDDUV1xztL_io, INN, MPI_DOUBLE_PRECISION, OctDDUV1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctDDUV2xztL_io, INN, MPI_DOUBLE_PRECISION, OctDDUV2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        !=======================================================
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctTUVxztL_io, INN, MPI_DOUBLE_PRECISION, OctTUV_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctTVzxztL_io, INN, MPI_DOUBLE_PRECISION, OctTVz_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctTTKxztL_io, INN, MPI_DOUBLE_PRECISION, OctTTK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctTDRxztL_io, INN, MPI_DOUBLE_PRECISION, OctTDR_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctTDUV1xztL_io, INN, MPI_DOUBLE_PRECISION, OctTDUV1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        INN = N2DO(MYID)*8*QUADHN
        CALL MPI_GATHER( OctTDUV2xztL_io, INN, MPI_DOUBLE_PRECISION, OctTDUV2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
        CALL MPI_BARRIER(ICOMM,IERROR)
        
        IF(thermlflg ==1) THEN
            !===============================
            INN = N2DO(MYID)
            CALL MPI_GATHER( T1xztL_io, INN, MPI_DOUBLE_PRECISION, T1_AUX, INN, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)
            CALL MPI_GATHER( D1xztL_io, INN, MPI_DOUBLE_PRECISION, D1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)
            CALL MPI_GATHER( H1xztL_io, INN, MPI_DOUBLE_PRECISION, H1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)
            CALL MPI_GATHER( M1xztL_io, INN, MPI_DOUBLE_PRECISION, M1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            !===============================
            INN = N2DO(MYID)
            CALL MPI_GATHER( T2xztL_io, INN, MPI_DOUBLE_PRECISION, T2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)
            CALL MPI_GATHER( D2xztL_io, INN, MPI_DOUBLE_PRECISION, D2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)
            CALL MPI_GATHER( H2xztL_io, INN, MPI_DOUBLE_PRECISION, H2_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            !================================
            INN = N2DO(MYID)
            CALL MPI_GATHER( DHxztL_io, INN, MPI_DOUBLE_PRECISION, DH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)
            CALL MPI_GATHER( PHxztL_io, INN, MPI_DOUBLE_PRECISION, PH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            !================================
            INN = N2DO(MYID)*NDV*NDV
            CALL MPI_GATHER( DVDL1MxztL_io,  INN, MPI_DOUBLE_PRECISION, DVDL1M_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV*NDV
            CALL MPI_GATHER( DVDL1MHxztL_io, INN, MPI_DOUBLE_PRECISION, DVDL1MH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV*NDV*NDV
            CALL MPI_GATHER( DVDL1MUxztL_io, INN, MPI_DOUBLE_PRECISION, DVDL1MU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
        
            INN = N2DO(MYID)*((NDV-1)*3+NDV)*((NDV-1)*3+NDV)
            CALL MPI_GATHER( DVDL2MxztL_io,  INN, MPI_DOUBLE_PRECISION, DVDL2M_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
!            !======test=============
!            IF (myid==0) THEN
!            Do J=1, N2DO(myid)
!                DO L=1,(NDV-1)*3+NDV
!                    DO M =1,(NDV-1)*3+NDV
!                        !WRITE(*,*) L, M, DVDL2xztL_io(J,L,M)*M1xztL_io(J), DVDL2MxztL_io(J,L,M)
!                    END DO
!                END DO
!            END DO
!            END IF
!            !======test=============
            
            !============================
            INN = N2DO(MYID)*NDV
            CALL MPI_GATHER( UHxztL_io, INN, MPI_DOUBLE_PRECISION, UH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV
            CALL MPI_GATHER( GHxztL_io, INN, MPI_DOUBLE_PRECISION, GH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*(NDV*(7-NDV)/2+NDV-3)
            CALL MPI_GATHER( U2DHxztL_io, INN, MPI_DOUBLE_PRECISION, U2DH_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            !============================
            INN = N2DO(MYID)*NDV
            CALL MPI_GATHER( DhDL1xztL_io, INN, MPI_DOUBLE_PRECISION, DhDL1_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV
            CALL MPI_GATHER( DhDLPxztL_io, INN, MPI_DOUBLE_PRECISION, DhDLP_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV
            CALL MPI_GATHER( DTDLKxztL_io, INN, MPI_DOUBLE_PRECISION, DTDLK_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV*NDV
            CALL MPI_GATHER( DTDLKUxztL_io, INN, MPI_DOUBLE_PRECISION, DTDLKU_AUX, INN,  MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV*NDV*NDV
            CALL MPI_GATHER( DTDLKDVDLxztL_io, INN, MPI_DOUBLE_PRECISION,DTDLKDVDL_AUX,INN, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
            INN = N2DO(MYID)*NDV*NDV*NDV
            CALL MPI_GATHER( DHDLMDVDLxztL_io, INN, MPI_DOUBLE_PRECISION,DHDLMDVDL_AUX,INN, MPI_DOUBLE_PRECISION,0,ICOMM,IERROR)
            CALL MPI_BARRIER(ICOMM,IERROR)
            
        END IF
        
        
        IF(MYID==0) THEN
            !==================================================================
            CALL MEMO_ALLOCT_AVERAGE_XZ_IO
            
            DO IP = 0, NPSLV
                N2DOID=JDEWT(IP)-JDSWT(IP)+1
                DO J=1,N2DOID
                    JJ=JDSWT(IP)-1+J
                   
                    !===============================
                    DO L=1,NDV+1
                        U1xztL_F0_io(JJ,L)=U1_AUX(J,L,IP+1)
                    ENDDO
                    
                    DO L=1,NDV
                        G1xztL_F0_io(JJ,L)=G1_AUX(J,L,IP+1)
                    ENDDO
                    
                    DO L=1,NDV
                        UPxztL_F0_io(JJ,L)=UP_AUX(J,L,IP+1)
                    ENDDO
                    
                    !===============================
                    DO L=1,(NDV*(7-NDV)/2+NDV-3)
                        U2xztL_F0_io(JJ,L)=U2_AUX(J,L,IP+1)
                    ENDDO
                    
                    DO L=1,(NDV*(7-NDV)/2+NDV-3)
                        UGxztL_F0_io(JJ,L)=UG_AUX(J,L,IP+1)
                    ENDDO
                    
                    DO L=1,(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
                        UGUxztL_F0_io(JJ,L)=UGU_AUX(J,L,IP+1)
                    ENDDO
                    
                    DO L=1,(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
                        U3xztL_F0_io(JJ,L)=U3_AUX(J,L,IP+1)
                    ENDDO
                    
                    !===============================
                    DO L=1,NDV
                        DO M =1,NDV
                            DVDL1xztL_F0_io(JJ,L,M)=DVDL1_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,NDV
                        DO M =1,NDV
                            DVDLPxztL_F0_io(JJ,L,M)=DVDLP_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,(NDV-1)*3+NDV
                        DO M =1,(NDV-1)*3+NDV
                            DVDL2xztL_F0_io(JJ,L,M)=DVDL2_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    !=========quadrant======================
                    DO L=1,4
                        DO M =1,QUADHN
                            QuadUVxztL_F0_io(JJ,L,M)=QuadUV_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,4
                        DO M =1,QUADHN
                            QuadVzxztL_F0_io(JJ,L,M)=DSQRT(QuadVz_AUX(J,L,M,IP+1))
                        END DO
                    END DO
                    
                    DO L=1,4
                        DO M =1,QUADHN
                            QuadTKxztL_F0_io(JJ,L,M)=QuadTK_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,4
                        DO M =1,QUADHN
                            QuadDRxztL_F0_io(JJ,L,M)=QuadDR_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,4
                        DO M =1,QUADHN
                            QuadDUV1xztL_F0_io(JJ,L,M)=QuadDUV1_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,4
                        DO M =1,QUADHN
                            QuadDUV2xztL_F0_io(JJ,L,M)=QuadDUV2_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    !================================
                    DO L=1,NDV+1
                        FUxztL_F0_io(JJ,L)=FU_AUX(J,L,IP+1)
                    ENDDO
                    
                    !=========Octant==\rho==================
                    DO L=1,8
                        DO M =1,QUADHN
                            OctDUVxztL_F0_io(JJ,L,M)=OctDUV_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctDVzxztL_F0_io(JJ,L,M)=DSQRT(OctDVz_AUX(J,L,M,IP+1))
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctDTKxztL_F0_io(JJ,L,M)=OctDTK_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctDDRxztL_F0_io(JJ,L,M)=OctDDR_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctDDUV1xztL_F0_io(JJ,L,M)=OctDDUV1_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctDDUV2xztL_F0_io(JJ,L,M)=OctDDUV2_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    !=========Octant==T===================
                    DO L=1,8
                        DO M =1,QUADHN
                            OctTUVxztL_F0_io(JJ,L,M)=OctTUV_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctTVzxztL_F0_io(JJ,L,M)=DSQRT(OctTVz_AUX(J,L,M,IP+1))
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctTTKxztL_F0_io(JJ,L,M)=OctTTK_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctTDRxztL_F0_io(JJ,L,M)=OctTDR_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctTDUV1xztL_F0_io(JJ,L,M)=OctTDUV1_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    
                    DO L=1,8
                        DO M =1,QUADHN
                            OctTDUV2xztL_F0_io(JJ,L,M)=OctTDUV2_AUX(J,L,M,IP+1)
                        END DO
                    END DO
                    !============================
                    IF(thermlflg ==1) THEN
                        !===============================
                        T1xztL_F0_io(JJ)=T1_AUX(J,IP+1)
                        D1xztL_F0_io(JJ)=D1_AUX(J,IP+1)
                        H1xztL_F0_io(JJ)=H1_AUX(J,IP+1)
                        M1xztL_F0_io(JJ)=M1_AUX(J,IP+1)
                        
                        !===============================
                        T2xztL_F0_io(JJ)=T2_AUX(J,IP+1)
                        D2xztL_F0_io(JJ)=D2_AUX(J,IP+1)
                        H2xztL_F0_io(JJ)=H2_AUX(J,IP+1)
                        
                        !===============================
                        DHxztL_F0_io(JJ)=DH_AUX(J,IP+1)
                        PHxztL_F0_io(JJ)=PH_AUX(J,IP+1)
                        !===============================
                        
                        DO L=1,NDV
                            DO M =1,NDV
                                DVDL1MxztL_F0_io (JJ,L,M)=DVDL1M_AUX (J,L,M,IP+1)
                                DVDL1MHxztL_F0_io(JJ,L,M)=DVDL1MH_AUX(J,L,M,IP+1)
                                DO N=1, NDV
                                    DVDL1MUxztL_F0_io(JJ,L,M,N)=DVDL1MU_AUX(J,L,M,N,IP+1)
                                END DO
                            END DO
                        END DO
                        
                        DO L=1,(NDV-1)*3+NDV
                            DO M =1,(NDV-1)*3+NDV
                                DVDL2MxztL_F0_io(JJ,L,M)=DVDL2M_AUX(J,L,M,IP+1)
                            END DO
                        END DO
                        
                        DO L=1,NDV
                            UHxztL_F0_io(JJ,L)   =UH_AUX(J,L,IP+1)
                            GHxztL_F0_io(JJ,L)   =GH_AUX(J,L,IP+1)
                        END DO
                        
                        DO L=1,(NDV*(7-NDV)/2+NDV-3)
                            U2DHxztL_F0_io(JJ,L)=U2DH_AUX(J,L,IP+1)
                        ENDDO
                        
                        DO L=1, NDV
                            DhDL1xztL_F0_io(JJ,L)=DhDL1_AUX(J,L,IP+1)
                            DhDLPxztL_F0_io(JJ,L)=DhDLP_AUX(J,L,IP+1)
                            DTDLKxztL_F0_io(JJ,L)=DTDLK_AUX(J,L,IP+1)
                            DO M=1, NDV
                                DTDLKUxztL_F0_io(JJ,L,M)=DTDLKU_AUX(J,L,M,IP+1)
                                DO N=1, NDV
                                    DTDLKDVDLxztL_F0_io(JJ,L,M,N)=DTDLKDVDL_AUX(J,L,M,N,IP+1)
                                    DHDLMDVDLxztL_F0_io(JJ,L,M,N)=DHDLMDVDL_AUX(J,L,M,N,IP+1)
                                END DO
                            END DO
                        ENDDO

                    ELSE 
                        ! below is to test FA recovering to RA...
                        CALl CHKHDL('NOTICE: This is only for isothermal flow!', myid)
                        D1xztL_F0_io = 1.0_WP
                        M1xztL_F0_io = 1.0_WP
                        DO L=1,NDV
                            DO M =1,NDV
                                DVDL1MxztL_F0_io (JJ,L,M)=DVDL1xztL_F0_io (JJ,L,M)
                                DO N=1, NDV
                                    DVDL1MUxztL_F0_io(JJ,L,M,N)=0.0_WP ! to get !!!?
                                END DO
                            END DO
                        END DO
                        
                        
                        DO L=1,(NDV-1)*3+NDV
                            DO M =1,(NDV-1)*3+NDV
                                DVDL2MxztL_F0_io (JJ,L,M)=DVDL2xztL_F0_io (JJ,L,M)
                            END DO
                        END DO
                        
                        
                    END IF
                END DO
            END DO

            !================check data=============================
!            DO L=1,(NDV*(6-NDV)+(NDV*(7-NDV))/2+NDV-8)
!            DO J=1, NCL2
!                !WRITE(*,*) L, J, UGUxztL_F0_io(J,L)/D1xztL_F0_io(J), U3xztL_F0_io(J,L), &
!                                UGUxztL_F0_io(J,L)/D1xztL_F0_io(J)-U3xztL_F0_io(J,L) !test
!            END DO
!            END DO 
            
            !=================adjust any pressure term to make wall pressure zero=======================
            Pwall=U1xztL_F0_io(1,4)
            DO J=1, NCL2
            
                U1xztL_F0_io(J,4)  = U1xztL_F0_io(J,4)  - Pwall
                
                DO M =1, NDV
                    UPxztL_F0_io(J,M)= UPxztL_F0_io(J,M)- Pwall*U1xztL_F0_io(J,M)
                END DO
                
                DO L=1,NDV
                    DO M =1,NDV
                        DVDLPxztL_F0_io(J,L,M)=DVDLPxztL_F0_io(J,L,M)-Pwall*DVDL1xztL_F0_io(J,L,M)
                    END DO
                END DO
                
                IF(thermlflg ==1) THEN
                    DO L=1, NDV
                        DhDLPxztL_F0_io(J,L)=DhDLPxztL_F0_io(J,L)-Pwall*DhDL1xztL_F0_io(J,L)
                    END DO
                END IF
                !!WRITE(*,*) YCC(J), U1xztL_F0_io(J,4) !test
            END DO
            !==================adjust any pressure term to make wall pressure zero==============
!            CALL CHKHDL(' ==> In WRT..DVDL2 VS DVDL2M.', myid)
!            Do J=1, NCL2
!                DO L=1,(NDV-1)*3+NDV
!                    DO M =1,(NDV-1)*3+NDV
!                        !WRITE(*,*) J, L, M, DVDL2xztL_F0_io (J,L,M), DVDL2MxztL_F0_io (J,L,M), &
!                                            DVDL2MxztL_F0_io (J,L,M)/DVDL2xztL_F0_io (J,L,M)
!                    END DO
!                END DO
!            END DO
            
            
        END IF
        
        RETURN
    END SUBROUTINE 

!=====================================================================================================================
    SUBROUTINE PP_FLOW_BASIC_VARS_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: J
        INTEGER(4) :: M, N, H, P, LMNH, LMH, LMN, LNH
        REAL(WP)   :: Matrix(3,3)
        REAL(WP)   :: eig(3)
        REAL(WP)   :: FF, TT
        REAL(WP)   :: DrivenFCTT1, DrivenFCTT2
        REAL(WP)   :: DrivenFCTTU1, DrivenFCTTU2
        INTEGER(4) :: TECFLG=200
        REAL(WP)      :: DENtemp, ddenintg,denmintg
        
        !====================================================================  
        !=======introduction=================================================
        !====================================================================  
        ! {*} = FA = Favred Averaged Mean.
        ! ''  = ff = Favred Averaged Fluctuation 
        ! <>  = RA = Reynolds Averaged Mean
        ! '   = f  = Reynolds Averaged Fluctuation.   
             
        !====================================================================  
        !======================RA based======================================
        !====================================================================  

        CALL CHKHDL('      Calculating basic variables',myid)
        !==============dPdX_RA(CL,m)=d<p>/dx_m ========================
        dPdX_RA(:,:) = 0.0_WP
        N=2
        DO J=1,NCL2
            IF(J==1) THEN
                dPdX_RA(J,N)= ( ( YCL2ND_WFB(J+1)*U1xztL_F0_io(J,  4) + &
                                  YCL2ND_WFF(J+1)*U1xztL_F0_io(J+1,4) ) - U1xztL_F0_io(J,4)  ) * DYFI(J)
                                
            ELSE IF (J==NCL2) THEN
                dPdX_RA(J,N)= (  U1xztL_F0_io(J,4) - &
                                ( YCL2ND_WFF(J)  *U1xztL_F0_io(J,  4) + &
                                  YCL2ND_WFB(J)  *U1xztL_F0_io(J-1,4) )  ) * DYFI(J)
                              
            ELSE
                dPdX_RA(J,N)= ( ( YCL2ND_WFB(J+1)*U1xztL_F0_io(J,  4) + &
                                  YCL2ND_WFF(J+1)*U1xztL_F0_io(J+1,4) ) - &
                                ( YCL2ND_WFF(J)  *U1xztL_F0_io(J,  4) + &
                                  YCL2ND_WFB(J)  *U1xztL_F0_io(J-1,4) ) ) * DYFI(J)
            END IF
        END DO
        CALL CHKHDL('      ==>Calculated d<P>/dy',myid)
        
        !====<p' u'_i>=<p' u"_i>=<p u_i> - <p>*<u_i>==ufpf_RA(J,i)=============
        DO J=1,NCL2
            DO M=1,NDV
                ufpf_RA(J,M)=  UPxztL_F0_io(J,M) - U1xztL_F0_io(J,4) * U1xztL_F0_io(J,M)
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated <u`_i p`>',myid)
        
        !=======================================
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    LMN = (M*(7-M))/2+N-3 
                    UU_RA(J,M,N)= U2xztL_F0_io(J,LMN)
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                        UU_RA(J,M,N)=UU_RA(J,N,M)
                    END IF
                END DO
            END DO
            !WRITE(*,'(4ES13.5)') YCC(J), UU_FA(J,1:3,1:3) !test    
        END DO
        CALL CHKHDL('      ==>Calculated <UU>',myid)
        
        !==============<u'_i u'_j> ========================
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    !LMN = (M*(7-M))/2+N-3 
                    uf2_RA(J,M,N)  = UU_RA(J,M,N) - U1xztL_F0_io(J,M) * U1xztL_F0_io(J,N)
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                        uf2_RA(J,M,N)  =uf2_RA(J,N,M)
                    END IF
                END DO
            END DO
            !WRITE(*,'(10ES13.5)') YCC(J), uf2_RA(J,1,1), uf2_RA(J,2,2),uf2_RA(J,3,3), &
            !            uf2_RA(J,1,2), uf2_RA(J,2,1),  &
            !            uf2_RA(J,1,3), uf2_RA(J,3,1),  &
            !            uf2_RA(J,2,3), uf2_RA(J,3,2)  ! test
        END DO
        CALL CHKHDL('      ==>Calculated <u`_i u`_j>',myid)
        
        DO J=1, NCL2
            DO M=1, NDV
                DO N=1, NDV
                    uf2d_RA(J,M,N)= uf2_RA(J,M,N)*D1xztL_F0_io(J)  
                END DO
                !WRITE(*,'(10ES13.5)') YCC(J), uf2d_RA(J,1,1), uf2d_RA(J,2,2),uf2d_RA(J,3,3), &
                !        uf2d_RA(J,1,2), uf2d_RA(J,2,1),  &
                !        uf2d_RA(J,1,3), uf2d_RA(J,3,1),  &
                !        uf2d_RA(J,2,3), uf2d_RA(J,3,2)  ! test
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated <\rho>*<u`_i u`_j>',myid)
        
        !=====u'u'u', u'u'v', u'u'w', u'v'v', u'v'w', u'w'w', v'v'v', v'v'w', v'w'w', w'w'w ===
        !{111} {112} {113}
        ![121] {122} {123}
        ![131] [132] {133}
        ![211] [212] [213]
        ![221] {222} {223}
        ![231] [232] {233}
        ![311] [312] [313]
        ![321] [322] [323]
        ![331] [332] {333}
!        ! Below is for test
!        DO J=1,NCL2
!            DO M=1,NDV
!                DO N=1,NDV
!                    IF(M.GT.N) CYCLE
!                    DO H=1,NDV
!                        IF(N.GT.H) CYCLE
!                        LMNH = M*(6-M)+(N*(7-N))/2+H-8  
!                        LMN = (M*(7-M))/2+N-3
!                        LMH = (M*(7-M))/2+H-3
!                        LNH = (N*(7-N))/2+H-3
                        
!                        U3xztL_F0_io(J,LMNH) = ( UGUxztL_F0_io(J,LMNH)  &
!                            +3.0_WP*D1xztL_F0_io(J)*U1xztL_F0_io(J,M)*U1xztL_F0_io(J,N)*U1xztL_F0_io(J,H) &
!                            -U1xztL_F0_io(J,M)*UGxztL_F0_io(J,LNH) &
!                            -U1xztL_F0_io(J,N)*UGxztL_F0_io(J,LMH) &
!                            -U1xztL_F0_io(J,H)*UGxztL_F0_io(J,LMN) &
!                            +U1xztL_F0_io(J,M)*U1xztL_F0_io(J,N)*G1xztL_F0_io(J,H) &
!                            +U1xztL_F0_io(J,M)*U1xztL_F0_io(J,H)*G1xztL_F0_io(J,N) &
!                            +U1xztL_F0_io(J,N)*U1xztL_F0_io(J,H)*G1xztL_F0_io(J,M) &
!                            +D1xztL_F0_io(J)*U1xztL_F0_io(J,M)*U2xztL_F0_io(J,LNH) &
!                            +D1xztL_F0_io(J)*U1xztL_F0_io(J,N)*U2xztL_F0_io(J,LMH) &
!                            +D1xztL_F0_io(J)*U1xztL_F0_io(J,H)*U2xztL_F0_io(J,LMN) )/D1xztL_F0_io(J) 
!                    END DO
!                END DO
!            END DO
!        END DO
!        ! above treatment is an approximat
!   As third order introduces more differences than the second order, thus above approximation is too rough to be used.
        
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    DO H=1,NDV
                        IF(N.GT.H) CYCLE
                        LMNH = M*(6-M)+(N*(7-N))/2+H-8  
                        LMN = (M*(7-M))/2+N-3
                        LMH = (M*(7-M))/2+H-3
                        LNH = (N*(7-N))/2+H-3
                        
                        uf3_RA(J, M,N,H) =  U3xztL_F0_io(J,LMNH)  &
                                          -U1xztL_F0_io(J,M) * UU_RA(J,N,H) &
                                          -U1xztL_F0_io(J,N) * UU_RA(J,M,H) &
                                          -U1xztL_F0_io(J,H) * UU_RA(J,M,N)  &
                                          +2.0_WP*U1xztL_F0_io(J,M)*U1xztL_F0_io(J,N)*U1xztL_F0_io(J,H)
                    END DO
                END DO
            END DO
                                            
            uf3_RA(J,1,2,1) = uf3_RA(J,1,1,2)
            uf3_RA(J,1,3,1) = uf3_RA(J,1,1,3)
            uf3_RA(J,1,3,2) = uf3_RA(J,1,2,3)
            
            uf3_RA(J,2,1,1) = uf3_RA(J,1,1,2)
            uf3_RA(J,2,1,2) = uf3_RA(J,1,2,2)
            uf3_RA(J,2,1,3) = uf3_RA(J,1,2,3)
            
            uf3_RA(J,2,2,1) = uf3_RA(J,1,2,2)
            uf3_RA(J,2,3,1) = uf3_RA(J,1,2,3)
            uf3_RA(J,2,3,2) = uf3_RA(J,2,2,3)
            
            uf3_RA(J,3,1,1) = uf3_RA(J,1,1,3)
            uf3_RA(J,3,1,2) = uf3_RA(J,1,2,3)
            uf3_RA(J,3,1,3) = uf3_RA(J,1,3,3)
            
            uf3_RA(J,3,2,1) = uf3_RA(J,1,2,3) 
            uf3_RA(J,3,2,2) = uf3_RA(J,2,2,3)
            uf3_RA(J,3,2,3) = uf3_RA(J,2,3,3)
            
            uf3_RA(J,3,3,1) = uf3_RA(J,1,3,3)
            uf3_RA(J,3,3,2) = uf3_RA(J,2,3,3)
            
            !WRITE(*,'(28ES13.5)') YCC(J), uf3_RA(J,1:3,1:3,1:3)  ! test

        END DO
        CALL CHKHDL('      ==>Calculated <u`_i u`_j u`_k>',myid)
        
         DO J=1, NCL2
            DO M=1, NDV
                DO N=1, NDV
                    DO H=1,NDV
                        uf3d_RA(J, M,N,H)= uf3_RA(J, M,N,H)*D1xztL_F0_io(J)  
                    END DO
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated <\rho>*<u`_i u`_j u`_k>',myid)
        
        
        !=============dialation of each volume======================
        DO J=1, NCL2
            dUidXi(J) = DVDL1xztL_F0_io(J,1,1) + &
                        DVDL1xztL_F0_io(J,2,2) + &
                        DVDL1xztL_F0_io(J,3,3)
            !WRITE(*,'(2ES13.5)') YCC(J), dUidXi(J) !test
        END DO
        CALL CHKHDL('      ==>Calculated d<u>/dx+d<v>/dy+d<w>/dz ',myid)
        
        !=============mean strain rate and vorticity tensor ==============
        DO J=1, NCL2
            DO M=1,NDV
                DO N=1,NDV
                    StrainTensor(J, M, N) = 0.5_WP*( DVDL1xztL_F0_io(J,M,N) + DVDL1xztL_F0_io(J,N,M) )
                    VortcyTensor(J, M, N) = 0.5_WP*( DVDL1xztL_F0_io(J,M,N) - DVDL1xztL_F0_io(J,N,M) )
                END DO
            END DO
            !WRITE(*,'(7ES13.5)') YCC(J), StrainTensor(J, 1:3, 1:3) !test
        END DO
        CALL CHKHDL('      ==>Calculated semi and asymetric parts of d<u_i>/dx_j ',myid)
        
        !==============skewness========================
        DO J=1,NCL2
            DO M=1,NDV
                Skewness_RA(J,M) = uf3_RA(J,M,M,M)/ ( dabs(uf2_RA(J,M,M))**(3.0_wp/2.0_wp) + RealMin)
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated Skewness',myid)
        
        DO J=1, NCL2
            MKE_RA(J) = 0.5_WP*( U1xztL_F0_io(J,1)*U1xztL_F0_io(J,1)+ &
                                 U1xztL_F0_io(J,2)*U1xztL_F0_io(J,2)+ &
                                 U1xztL_F0_io(J,3)*U1xztL_F0_io(J,3) ) *D1xztL_F0_io(J)
            TKE_RA(J) = 0.5_WP*( uf2_RA(J,1,1)+ &
                                 uf2_RA(J,2,2)+ &
                                 uf2_RA(J,3,3) )      *D1xztL_F0_io(J)
        END DO
        CALL CHKHDL('      ==>Calculated MKE_RA=<\rho>*<u_i>*<u_i>/2 and TKE_RA=<\rho>*<u`_i*u`_i>/2',myid)
        
        
        !====Vorticity==================================
        ! Omega1= \partial u_3/ \partial x_2 -  \partial u_2/ \partial x_3
        ! Omega2= \partial u_1/ \partial x_3 -  \partial u_3/ \partial x_1
        ! Omega3= \partial u_2/ \partial x_1 -  \partial u_1/ \partial x_2
        ! Omega_i = <Omega_i> + Omega'_i (or omega)
        ! <Omega'_i * Omega'_i> = <Omega_i * Omega_i> - <Omega_i><Omega_i>
        DO J=1, NCL2
            ! Omega2_RA(J,1) = (dw/dy * dw/dy)-2*dw/dy*dv/dz + (dv/dz*dv/dz)
            Omega2_RA(J,1) = DVDL2xztL_F0_io(J,(3-1)*NDV+2,(3-1)*NDV+2) - 2.0_wp* &
                             DVDL2xztL_F0_io(J,(3-1)*NDV+2,(2-1)*NDV+3) + &
                             DVDL2xztL_F0_io(J,(2-1)*NDV+3,(2-1)*NDV+3)
            
            Omega2_RA(J,2) = DVDL2xztL_F0_io(J,(1-1)*NDV+3,(1-1)*NDV+3) - 2.0_wp* &
                             DVDL2xztL_F0_io(J,(1-1)*NDV+3,(3-1)*NDV+1) + &
                             DVDL2xztL_F0_io(J,(3-1)*NDV+1,(3-1)*NDV+1)
                             
            Omega2_RA(J,3) = DVDL2xztL_F0_io(J,(2-1)*NDV+1,(2-1)*NDV+1) - 2.0_wp* &
                             DVDL2xztL_F0_io(J,(2-1)*NDV+1,(1-1)*NDV+2) + &
                             DVDL2xztL_F0_io(J,(1-1)*NDV+2,(1-1)*NDV+2)
                             
            Omega_RA2(J,1) = DVDL1xztL_F0_io(J,3,2)*DVDL1xztL_F0_io(J,3,2) - 2.0_wp * &
                             DVDL1xztL_F0_io(J,3,2)*DVDL1xztL_F0_io(J,2,3) + &
                             DVDL1xztL_F0_io(J,2,3)*DVDL1xztL_F0_io(J,2,3)
                             
            Omega_RA2(J,2) = DVDL1xztL_F0_io(J,1,3)*DVDL1xztL_F0_io(J,1,3) - 2.0_wp * &
                             DVDL1xztL_F0_io(J,1,3)*DVDL1xztL_F0_io(J,3,1) + &
                             DVDL1xztL_F0_io(J,3,1)*DVDL1xztL_F0_io(J,3,1)
                             
            Omega_RA2(J,3) = DVDL1xztL_F0_io(J,2,1)*DVDL1xztL_F0_io(J,2,1) - 2.0_wp * &
                             DVDL1xztL_F0_io(J,2,1)*DVDL1xztL_F0_io(J,1,2) + &
                             DVDL1xztL_F0_io(J,1,2)*DVDL1xztL_F0_io(J,1,2)
             
            omega_rms(J,1) = DSQRT( dabs(Omega2_RA(J,1) - Omega_RA2(J,1)) )
            omega_rms(J,2) = DSQRT( dabs(Omega2_RA(J,2) - Omega_RA2(J,2)) )
            omega_rms(J,3) = DSQRT( dabs(Omega2_RA(J,3) - Omega_RA2(J,3)) )
        END DO
        CALL CHKHDL('      ==>Calculated omega_rms',myid)
        
        DO J=1, NCL2
            DO M=1,NDV
                DO N=1,NDV
                    ANISOTROPY_RA(J,M,N) = uf2_RA(J,M,N)/(uf2_RA(J,1,1)+uf2_RA(J,2,2)+uf2_RA(J,3,3)) - &
                                        DBLE(Kronecker_delta(M,N))/3.0_WP 
                END DO
            END DO
!            !WRITE(*,*) '  --  '
!            !WRITE(*,*) ANISOTROPY_RA(J,1,1:3)
!            !WRITE(*,*) ANISOTROPY_RA(J,2,1:3)
!            !WRITE(*,*) ANISOTROPY_RA(J,3,1:3)
!            !WRITE(*,*) '  --  '
        END DO
        
        DO J=1, NCL2
             ! below two methods are the same. Checked Good!
!            ANISTPinva_RA(J,1) = 0.0_WP
!            ANISTPinva_RA(J,2) = -( ANISOTROPY_RA(J,1,1)*ANISOTROPY_RA(J,1,1) + &
!                                    ANISOTROPY_RA(J,1,2)*ANISOTROPY_RA(J,1,2) + &
!                                    ANISOTROPY_RA(J,1,3)*ANISOTROPY_RA(J,1,3) + &
!                                    ANISOTROPY_RA(J,2,1)*ANISOTROPY_RA(J,2,1) + &
!                                    ANISOTROPY_RA(J,2,2)*ANISOTROPY_RA(J,2,2) + &
!                                    ANISOTROPY_RA(J,2,3)*ANISOTROPY_RA(J,2,3) + &
!                                    ANISOTROPY_RA(J,3,1)*ANISOTROPY_RA(J,3,1) + &
!                                    ANISOTROPY_RA(J,3,2)*ANISOTROPY_RA(J,3,2) + &
!                                    ANISOTROPY_RA(J,3,3)*ANISOTROPY_RA(J,3,3) )/2.0_WP
!            ANISTPinva_RA(J,3) = (  ANISOTROPY_RA(J,1,1)*ANISOTROPY_RA(J,1,1)*ANISOTROPY_RA(J,1,1) + &
!                                    ANISOTROPY_RA(J,1,2)*ANISOTROPY_RA(J,2,1)*ANISOTROPY_RA(J,1,1) + &
!                                    ANISOTROPY_RA(J,1,3)*ANISOTROPY_RA(J,3,1)*ANISOTROPY_RA(J,1,1) + &
!                                    ANISOTROPY_RA(J,2,1)*ANISOTROPY_RA(J,1,1)*ANISOTROPY_RA(J,1,2) + &
!                                    ANISOTROPY_RA(J,2,2)*ANISOTROPY_RA(J,2,1)*ANISOTROPY_RA(J,1,2) + &
!                                    ANISOTROPY_RA(J,2,3)*ANISOTROPY_RA(J,3,1)*ANISOTROPY_RA(J,1,2) + &
!                                    ANISOTROPY_RA(J,3,1)*ANISOTROPY_RA(J,1,1)*ANISOTROPY_RA(J,1,3) + &
!                                    ANISOTROPY_RA(J,3,2)*ANISOTROPY_RA(J,2,1)*ANISOTROPY_RA(J,1,3) + &
!                                    ANISOTROPY_RA(J,3,3)*ANISOTROPY_RA(J,3,1)*ANISOTROPY_RA(J,1,3) + &
!                                    ANISOTROPY_RA(J,1,1)*ANISOTROPY_RA(J,1,2)*ANISOTROPY_RA(J,2,1) + &
!                                    ANISOTROPY_RA(J,1,2)*ANISOTROPY_RA(J,2,2)*ANISOTROPY_RA(J,2,1) + &
!                                    ANISOTROPY_RA(J,1,3)*ANISOTROPY_RA(J,3,2)*ANISOTROPY_RA(J,2,1) + &
!                                    ANISOTROPY_RA(J,2,1)*ANISOTROPY_RA(J,1,2)*ANISOTROPY_RA(J,2,2) + &
!                                    ANISOTROPY_RA(J,2,2)*ANISOTROPY_RA(J,2,2)*ANISOTROPY_RA(J,2,2) + &
!                                    ANISOTROPY_RA(J,2,3)*ANISOTROPY_RA(J,3,2)*ANISOTROPY_RA(J,2,2) + &
!                                    ANISOTROPY_RA(J,3,1)*ANISOTROPY_RA(J,1,2)*ANISOTROPY_RA(J,2,3) + &
!                                    ANISOTROPY_RA(J,3,2)*ANISOTROPY_RA(J,2,2)*ANISOTROPY_RA(J,2,3) + &
!                                    ANISOTROPY_RA(J,3,3)*ANISOTROPY_RA(J,3,2)*ANISOTROPY_RA(J,2,3) + &
!                                    ANISOTROPY_RA(J,1,1)*ANISOTROPY_RA(J,1,3)*ANISOTROPY_RA(J,3,1) + &
!                                    ANISOTROPY_RA(J,1,2)*ANISOTROPY_RA(J,2,3)*ANISOTROPY_RA(J,3,1) + &
!                                    ANISOTROPY_RA(J,1,3)*ANISOTROPY_RA(J,3,3)*ANISOTROPY_RA(J,3,1) + &
!                                    ANISOTROPY_RA(J,2,1)*ANISOTROPY_RA(J,1,3)*ANISOTROPY_RA(J,3,2) + &
!                                    ANISOTROPY_RA(J,2,2)*ANISOTROPY_RA(J,2,3)*ANISOTROPY_RA(J,3,2) + &
!                                    ANISOTROPY_RA(J,2,3)*ANISOTROPY_RA(J,3,3)*ANISOTROPY_RA(J,3,2) + &
!                                    ANISOTROPY_RA(J,3,1)*ANISOTROPY_RA(J,1,3)*ANISOTROPY_RA(J,3,3) + &
!                                    ANISOTROPY_RA(J,3,2)*ANISOTROPY_RA(J,2,3)*ANISOTROPY_RA(J,3,3) + &
!                                    ANISOTROPY_RA(J,3,3)*ANISOTROPY_RA(J,3,3)*ANISOTROPY_RA(J,3,3) )/3.0_WP
                                    
!            LumleyAxis_RA(J,1) = DSQRT(-ANISTPinva_RA(J,2)/3.0_WP)
!            LumleyAxis_RA(J,2) = (ANISTPinva_RA(J,3)/2.0_WP)**(1.0_wp/3.0_wp)
!            !WRITE(*,*) 'invars', J, ANISTPinva_RA(J,1:3), LumleyAxis_RA(J,1:2)
            
            Matrix(1:3,1:3)=ANISOTROPY_RA(J,1:3,1:3)
            call DSYEVC3(Matrix,eig)
            ANISTPinva_RA(J,1) = 0.0_WP
            ANISTPinva_RA(J,2) = (eig(1)*eig(1)+eig(1)*eig(2)+eig(2)*eig(2))*(-1.0_WP)
            ANISTPinva_RA(J,3) = eig(1)*eig(2)*(eig(1)+eig(2))*(-1.0_WP)
            LumleyAxis_RA(J,1) = DSQRT(-ANISTPinva_RA(J,2)/3.0_WP)
            LumleyAxis_RA(J,2) = (ANISTPinva_RA(J,3)/2.0_WP)**(1.0_wp/3.0_wp)
            
            !!WRITE(*,*) 'invars', J, ANISTPinva_RA(J,1:3), LumleyAxis_RA(J,1:2)
            !!WRITE(*,*) 'eigenv', J, eig(1:3), -eig(1)-eig(2), eig(1)+eig(2)+eig(3)
            
        END DO
        
        
!====================================================================    
!======================FA based======================================  
!====================================================================  
  
        !===============U_FA(J,M)={u_M}==========================
        ! Eq. U_FA(J,M)={u_M} = <\rho u_M>/ <\rho>
        !     uff_RA(J,M) = <u"_M>=<u_M>-{u_M}
        U_FA   = 0.0_WP
        uff_RA = 0.0_WP
        DO J=1,NCL2
            DO M=1, NDV
                U_FA(J,M)  = G1xztL_F0_io(J,M)/D1xztL_F0_io(J)
                uff_RA(J,M)= U1xztL_F0_io(J,M) - U_FA(J,M)
            END DO
            !WRITE(*,'(3ES13.5)') YCC(J), U_FA(J,1:3), uff_RA(J,1:3) !test
        END DO  
        CALL CHKHDL('      ==>Calculated {U} and <u">',myid)
        
        
        !==============dUdX_FA(CL,m,n)=d{u_m}/dx_n========================
        dUdX_FA = 0.0_WP
        N=2
        DO M=1,NDV      ! u-components
            
            DO J=1,NCL2
                IF(J==1) THEN
                    dUdX_FA(J,M,N)= ( ( YCL2ND_WFB(J+1)*U_FA(J,  M) + &
                                        YCL2ND_WFF(J+1)*U_FA(J+1,M) ) - 0.0_WP  ) * DYFI(J)                
                ELSE IF (J==NCL2) THEN
                    dUdX_FA(J,M,N)= (  0.0_WP - &
                                    ( YCL2ND_WFF(J)*U_FA(J,  M) + &
                                      YCL2ND_WFB(J)*U_FA(J-1,M) )  ) * DYFI(J)              
                ELSE
                    dUdX_FA(J,M,N)= ( ( YCL2ND_WFB(J+1)*U_FA(J,  M) + &
                                        YCL2ND_WFF(J+1)*U_FA(J+1,M) ) - &
                                      ( YCL2ND_WFF(J)  *U_FA(J,  M) + &
                                        YCL2ND_WFB(J)  *U_FA(J-1,M) ) ) * DYFI(J)
                END IF
            END DO
            !WRITE(*,'(4ES13.5)') YCC(J), dUdX_FA(J,1:3,2) !test    
        END DO
        CALL CHKHDL('      ==>Calculated d{U}/dx',myid)
        
        !======UU_FA(CL,i,j) = {u_i u_j} = <\rho u_i u_j>/<\rho>  Not Equal to <u"_u u"_j>=======
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    LMN = (M*(7-M))/2+N-3 
                    UU_FA(J,M,N)= UGxztL_F0_io(J,LMN)/D1xztL_F0_io(J)
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                        UU_FA(J,M,N)=UU_FA(J,N,M)
                    END IF
                END DO
            END DO
            !WRITE(*,'(4ES13.5)') YCC(J), UU_FA(J,1:3,1:3) !test    
        END DO
        CALL CHKHDL('      ==>Calculated {UU}',myid)
        
        !=============={u"_i u"_j} ========================
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    uff2_FA(J,M,N)= UU_FA(J,M,N)  - U_FA(J,M)* U_FA(J,N)
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                        uff2_FA(J,M,N)=uff2_FA(J,N,M)
                    END IF
                END DO
            END DO
            !WRITE(*,'(4ES13.5)') YCC(J), uff2_FA(J,1:3,1:3)
        END DO
        CALL CHKHDL('      ==>Calculated {u"_i u"_j }',myid)
        
        DO J=1, NCL2
            DO M=1,NDV
                DO N=1,NDV
                    uff2d_FA(J,M,N) = uff2_FA(J,M,N) * D1xztL_F0_io(J)
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated <\rho>*{u"_i u"_j }',myid)
        
        !======{\rho u_2 u"_i u"_j}===for checking turb diffu, and exactly the same as the below one.
!        DO J=1,NCL2
        
!            DO M=1,NDV
!                DO N=1,NDV
!                    IF(M.GT.N) CYCLE   
!                    IF(M==1 .and. N==1)   LMNH =   2 
!                    IF(M==1 .and. N==2)   LMNH =   4 
!                    IF(M==1 .and. N==3)   LMNH =   5 
!                    IF(M==2 .and. N==2)   LMNH =   7 
!                    IF(M==2 .and. N==3)   LMNH =   8 
!                    IF(M==3 .and. N==3)   LMNH =   9 
                                 
!                    TDIFU_FA(J,M,N) =  UGUxztL_F0_io(J,LMNH) &
!                                        -U_FA(J,M) * UU_FA(J,N,2)*D1xztL_F0_io(J) &
!                                        -U_FA(J,N) * UU_FA(J,M,2)*D1xztL_F0_io(J) &
!                                        +U_FA(J,M) * U_FA(J,N) * U_FA(J,2)*D1xztL_F0_io(J)
!                END DO
!            END DO 
!            DO M=1,NDV
!                DO N=1,NDV
!                    IF(M.GT.N) THEN
!                        TDIFU_FA(J,M,N)=TDIFU_FA(J,N,M)
!                    END IF
!                END DO
!            END DO
!        END DO
        
        
        !==={u"_k u"_i u"_j}==============================
        !=====u"u"u", u"u"v'', u"u"w'', u"v''v'', u"v''w'' =
        !=====u"w''w'', v''v''v'', v''v''w'', v''w''w'', w''w''w'' ===
        !{111} {112} {113}
        ![121] {122} {123}
        ![131] [132] {133}
        ![211] [212] [213]
        ![221] {222} {223}
        ![231] [232] {233}
        ![311] [312] [313]
        ![321] [322] [323]
        ![331] [332] {333}
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    DO H=1,NDV
                        IF(N.GT.H) CYCLE
                        LMNH = M*(6-M)+(N*(7-N))/2+H-8                          
                        uff3_FA(J,M,N,H) =  UGUxztL_F0_io(J,LMNH)/D1xztL_F0_io(J) &
                                            -U_FA(J,M) * UU_FA(J,N,H) &
                                            -U_FA(J,N) * UU_FA(J,M,H) &
                                            -U_FA(J,H) * UU_FA(J,M,N) &
                                            +2.0_WP * U_FA(J,M) * U_FA(J,N) * U_FA(J,H)
                        !!WRITE(*,*) M,N, H, uff3_FA(J,M,N,H),uf3_RA(J,M,N,H),uff3_FA(J,M,N,H)-uf3_RA(J,M,N,H)
                    END DO
                END DO
            END DO
            
           !uff3_FA(J,1,1,1)  ! #1
           !uff3_FA(J,1,1,2)  ! #2
           !uff3_FA(J,1,1,3)  ! #3
            
            uff3_FA(J,1,2,1) = uff3_FA(J,1,1,2) !2
           !uff3_FA(J,1,2,2)  ! #4
           !uff3_FA(J,1,2,3)  ! #5
           
            uff3_FA(J,1,3,1) = uff3_FA(J,1,1,3) !3
            uff3_FA(J,1,3,2) = uff3_FA(J,1,2,3) !5
           !uff3_FA(J,1,3,3)  ! #6
           
            uff3_FA(J,2,1,1) = uff3_FA(J,1,1,2) !2
            uff3_FA(J,2,1,2) = uff3_FA(J,1,2,2) !4
            uff3_FA(J,2,1,3) = uff3_FA(J,1,2,3) !5
            
            uff3_FA(J,2,2,1) = uff3_FA(J,1,2,2) !4
           !uff3_FA(J,2,2,2)  ! #7
           !uff3_FA(J,2,2,3)  ! #8
           
            uff3_FA(J,2,3,1) = uff3_FA(J,1,2,3) !5
            uff3_FA(J,2,3,2) = uff3_FA(J,2,2,3) !8
           !uff3_FA(J,2,3,3)  ! #9
           
            uff3_FA(J,3,1,1) = uff3_FA(J,1,1,3) !3
            uff3_FA(J,3,1,2) = uff3_FA(J,1,2,3) !5
            uff3_FA(J,3,1,3) = uff3_FA(J,1,3,3) !6
            
            uff3_FA(J,3,2,1) = uff3_FA(J,1,2,3) !5
            uff3_FA(J,3,2,2) = uff3_FA(J,2,2,3) !8
            uff3_FA(J,3,2,3) = uff3_FA(J,2,3,3) !9
            
            uff3_FA(J,3,3,1) = uff3_FA(J,1,3,3) !6
            uff3_FA(J,3,3,2) = uff3_FA(J,2,3,3) !9
           !uff3_FA(J,3,3,3)  ! #10
            !WRITE(*,'(28ES13.5)') YCC(J), uff3_FA(J,1:3,1:3,1:3)  ! test
        END DO 
        
!        !WRITE(*,*) 'Check UUU'
!        DO J=1, NCL2
!!            !WRITE(*,*) TDIFU_FA(J,1,1), TDIFU_FA(J,1,2), TDIFU_FA(J,1,3), TDIFU_FA(J,2,2), TDIFU_FA(J,2,3), TDIFU_FA(J,3,3)
!            WRITE(*,'(9ES13.5)') YCC(J), uff3_FA (J,1,1,1),uf3_RA (J,1,1,1), &
!            UGUxztL_F0_io(J,1)/D1xztL_F0_io(J), U3xztL_F0_io(J,1), &
!            U_FA(J,1), U1xztL_F0_io(J,1), &
!            UU_FA(J,1,1), UU_RA(J,1,1)
!        END DO
! Conclusion: Above is correct, U3 indeed introduces large variations from UGU/D
        
        CALL CHKHDL('      ==>Calculated {u"_i u"_j u"_k}',myid)
        
        DO J=1, NCL2
            DO M=1,NDV
                DO N=1,NDV
                    DO H=1, NDV
                       uff3d_FA(J,M,N,H) = uff3_FA(J,M,N,H) * D1xztL_F0_io(J)
                    END DO
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated <\rho> * {u"_i u"_j u"_k}',myid)
        
        !=============dialation of each volume======================
        DO J=1, NCL2
            dUidXiM(J)= (   DVDL1MxztL_F0_io(J,1,1) + &
                            DVDL1MxztL_F0_io(J,2,2) + &
                            DVDL1MxztL_F0_io(J,3,3)    ) * CVISC
        END DO
        CALL CHKHDL('      ==>Calculated <\mu dU/dx> + <\mu dV/dy> + <\mu dW/dz>',myid)
        
        !=============mean strain rate and vorticity tensor ==============
        DO J=1, NCL2
            DO M=1,NDV
                DO N=1,NDV
                    StrainTensorM(J, M, N) = 0.5_WP*( DVDL1MxztL_F0_io(J,M,N) + DVDL1MxztL_F0_io(J,N,M) )* CVISC
                    VortcyTensorM(J, M, N) = 0.5_WP*( DVDL1MxztL_F0_io(J,M,N) - DVDL1MxztL_F0_io(J,N,M) )* CVISC
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated \mu S and \mu\Omega',myid)
        
        !==============skewness========================
        
        DO M=1,NDV
            DO J=1,NCL2
                Skewness_FA(J,M) = uff3_FA(J,M,M,M)/ ( dabs(uff2_FA(J,M,M))**(3.0_wp/2.0_wp) + RealMin)
                !WRITE(*,'(2I3.1, 3(3ES13.5, 2X))') M, J, Skewness_FA(J,M),Skewness_RA(J,M),&
                !                (Skewness_FA(J,M)-Skewness_RA(J,M))/Skewness_FA(J,M), &
                !                 uff3_FA(J,M,M,M),uf3_RA(J,M,M,M),(uff3_FA(J,M,M,M)-uf3_RA(J,M,M,M))/uff3_FA(J,M,M,M), &
                !                 uff2_FA(J,M,M),uf2_RA(J,M,M),(uff2_FA(J,M,M)-uf2_RA(J,M,M))/uff2_FA(J,M,M)
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated Skewness_FA',myid)
        
        
        DO J=1, NCL2
            MKE_FA(J) = 0.5_WP*( U_FA(J,1)*U_FA(J,1)+ &
                                 U_FA(J,2)*U_FA(J,2)+ &
                                 U_FA(J,3)*U_FA(J,3) ) *D1xztL_F0_io(J)
            TKE_FA(J) = 0.5_WP*( uff2_FA(J,1,1)+ &
                                 uff2_FA(J,2,2)+ &
                                 uff2_FA(J,3,3) )      *D1xztL_F0_io(J)
        END DO
        CALL CHKHDL('      ==>Calculated MKE_FA and TKE_FA',myid)
        
        
        !===========Viscous shear stress=====================
        !=========<tau_mn>(<u>) AND <tau_mn>(u')============
        ! Eq. Tau_Mean_RA(J,M,N) = ViscStress_Tau_Umea + ViscStress_Tau_Uper
        !     <tau_mn>(<u>)*REN= <mu> [ (\partial <u_m>)/(\partial x_n) + 
        !                               (\partial <u_n>)/(\partial x_m) ] -
        !                     2/3<mu> [ (\partial <u_l>)/(\partial x_l) \delta_mn]
        !     <tau_mn>(u')*REN= < mu' [ (\partial u'_m)/(\partial x_n) + 
        !                               (\partial u'_n)/(\partial x_m) ] > -
        !                     2/3<mu' [ (\partial u'_l)/(\partial x_l) \delta_mn]
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    !==Eq.2.20 in Huang1995====
                    Tau_Mean_RA(J,M,N) = 2.0_wp * &
                                                ( StrainTensorM(J,M,N)- dUidXiM(J)*DBLE(Kronecker_delta(M,N))/3.0_WP  )
                    Tau_meaU_RA(J,M,N) = 2.0_wp * M1xztL_F0_io(J) *CVISC * &
                                                ( StrainTensor (J,M,N)- dUidXi(J) *DBLE(Kronecker_delta(M,N))/3.0_WP  )
                END DO
            END DO
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                         Tau_Mean_RA(J,M,N) = Tau_Mean_RA(J,N,M)
                         Tau_meaU_RA(J,M,N) = Tau_meaU_RA(J,N,M)
                    END IF
                END DO
            END DO
            
        END DO
        CALL CHKHDL('      ==>Calculated <\mu S*>, <\mu>*<S*>, <\mu` S`>',myid)
        
        
        !========TauU_RA(J, M, N, H)=<u_h \tau_mn>======================
       !Eq.<u_h \tau_mn> = <\mu u_h \partial{u_m}/\partial{u_n}>+
        !                  <\mu u_h \partial{u_n}/\partial{u_m}>-
        !               2/3<\mu u_h \partial{u_l}/\partial{u_l}>delta_mn
        DO J=1, NCL2    
            DO M=1,NDV
                DO N=1, NDV
                    DO H=1, NDV
                        TauU_RA(J, M, N, H)= ( & 
                                            DVDL1MUxztL_F0_io(J,M,N,H) + &
                                            DVDL1MUxztL_F0_io(J,N,M,H) - &
                            2.0_WP/3.0_WP*( DVDL1MUxztL_F0_io(J,1,1,H) + &
                                            DVDL1MUxztL_F0_io(J,2,2,H) + &
                                            DVDL1MUxztL_F0_io(J,3,3,H) )*DBLE(Kronecker_delta(M,N)) ) * CVISC
                        Taufuf_RA(J, M, N, H) = TauU_RA(J, M, N, H) - &
                                            U1xztL_F0_io(J,H)*Tau_Mean_RA(J,M,N)
                    END DO
                END DO
            END DO
            
            !===below for test== test OK===
            !!WRITE(*,*) TauU_RA(J, 1, 2, 1:3)-TauU_RA(J, 2, 1, 1:3)
            !!WRITE(*,*) TauU_RA(J, 1, 3, 1:3)-TauU_RA(J, 3, 1, 1:3)
            !!WRITE(*,*) TauU_RA(J, 2, 3, 1:3)-TauU_RA(J, 3, 2, 1:3)
            !test
            
        END DO 
        CALL CHKHDL('      ==>Calculated TauU_RA',myid)
        
        !=========d<tau_mn>(<u>,u')/dy=dTaudy_RA(J,M,N)============
        DO J=1, NCL2
        
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    IF(J==1) THEN
                        if (M==1 .and. N==2) THEN
                            FF=DABS(Tauw_io(1))
                        ELSE
                            FF=0.0_WP
                        END IF
                        dTaudy_RA(J,M,N) = &
                            ( ( YCL2ND_WFF(J+1)*Tau_Mean_RA(J+1,M,N)+ &
                                YCL2ND_WFB(J+1)*Tau_Mean_RA(J,  M,N) ) -  &
                              FF )*DYFI(J)

                    ELSE IF(J==NCL2) THEN
                        if (M==1 .and. N==2) THEN
                            FF=-DABS(Tauw_io(2))
                        ELSE
                            FF=0.0_WP
                        END IF
                        dTaudy_RA(J,M,N) = &
                            ( FF -  &
                              ( YCL2ND_WFF(J)  *Tau_Mean_RA(J,  M,N)+ &
                                YCL2ND_WFB(J)  *Tau_Mean_RA(J-1,M,N) ) )*DYFI(J)

                    ELSE
                        dTaudy_RA(J,M,N) = &
                            ( ( YCL2ND_WFF(J+1)*Tau_Mean_RA(J+1,M,N)+ &
                                YCL2ND_WFB(J+1)*Tau_Mean_RA(J,  M,N) ) -  &
                              ( YCL2ND_WFF(J)  *Tau_Mean_RA(J,  M,N)+ &
                                YCL2ND_WFB(J)  *Tau_Mean_RA(J-1,M,N) ) )*DYFI(J)
                    END IF
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                         dTaudy_RA(J,M,N) = dTaudy_RA(J,N,M)
                    END IF
                END DO
            END DO
            
        END DO
        CALL CHKHDL('      ==>Calculated dTaudy_RA',myid)
        
        
        !=========d<R_mn>/dy=dTSSdy_RA(J,M,N)============
        !=========d<\rho u"_m u"_n>/dy = d (<\rho> {u"_m u"_n})/dy
        DO J=1, NCL2
        
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    IF(J==1) THEN
                        dTSSdy_RA(J,M,N) = &
                            ( ( YCL2ND_WFF(J+1)*uff2d_FA(J+1,M,N)+ &
                                YCL2ND_WFB(J+1)*uff2d_FA(J,  M,N) ) -  &
                              0.0_WP )*DYFI(J)

                    ELSE IF(J==NCL2) THEN
                        dTSSdy_RA(J,M,N) = &
                            ( 0.0_WP -  &
                              ( YCL2ND_WFF(J)  *uff2d_FA(J,  M,N)+ &
                                YCL2ND_WFB(J)  *uff2d_FA(J-1,M,N) ) )*DYFI(J)

                    ELSE
                        dTSSdy_RA(J,M,N) = &
                            ( ( YCL2ND_WFF(J+1)*uff2d_FA(J+1,M,N)+ &
                                YCL2ND_WFB(J+1)*uff2d_FA(J,  M,N) ) -  &
                              ( YCL2ND_WFF(J)  *uff2d_FA(J,  M,N)+ &
                                YCL2ND_WFB(J)  *uff2d_FA(J-1,M,N) ) )*DYFI(J)
                    END IF
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                         dTSSdy_RA(J,M,N) = dTSSdy_RA(J,N,M)
                    END IF
                END DO
            END DO
            
        END DO
        CALL CHKHDL('      ==>Calculated dTSSdy_RA',myid)
        
        
        !========<du_m/dx_n * \tau_hp >======================
       ! Eq. = <d(u_m)/d(x_n) * mu * d(u_h)/d(x_p) > +
       !       <d(u_m)/d(x_n) * mu * d(u_p)/d(x_h) > - 2/3*delta_hp
       !       <d(u_m)/d(x_n) * mu * d(u_l)/d(x_l) >
       ! (M-1)*NDV+H
        DO J=1, NCL2    
            DO M=1,NDV         
                DO N=1, NDV   
                    DO H=1, NDV           
                        DO P=1, NDV    
                            TauDvDL_RA(J, M, N, H, P) = (     &
                                DVDL2MxztL_F0_io(J, (M-1)*3+N, (H-1)*3+P ) + &
                                DVDL2MxztL_F0_io(J, (M-1)*3+N, (P-1)*3+H ) - &
                                2.0_wp/3.0_wp* DBLE(Kronecker_delta(H,P))* (& 
                                DVDL2MxztL_F0_io(J, (M-1)*3+N, (1-1)*3+1 ) + &
                                DVDL2MxztL_F0_io(J, (M-1)*3+N, (2-1)*3+2 ) + &
                                DVDL2MxztL_F0_io(J, (M-1)*3+N, (3-1)*3+3 ) ) ) * CVISC
                            
                        END DO
                    END DO
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated TauDvDL_RA',myid)
        
        
!        DO J=1, NCL2    
!            DO M=1,NDV         
!                DO N=1, NDV   
!                    Tau_ik_Du_jDx_i_RA(J,M,N) = ( &
!                            DVDL2MxztL_F0_io(J, (M-1)*3+1, (N-1)*3+1 ) + &
!                            DVDL2MxztL_F0_io(J, (M-1)*3+2, (N-1)*3+2 ) + &
!                            DVDL2MxztL_F0_io(J, (M-1)*3+3, (N-1)*3+3 ) + &
!                            DVDL2MxztL_F0_io(J, (1-1)*3+M, (N-1)*3+1 ) + &
!                            DVDL2MxztL_F0_io(J, (2-1)*3+M, (N-1)*3+2 ) + &
!                            DVDL2MxztL_F0_io(J, (3-1)*3+M, (N-1)*3+3 ) - &
!                            2.0_wp/3.0_wp* ( &
!                            DVDL2MxztL_F0_io(J, (1-1)*3+1, (N-1)*3+1 )*DBLE(Kronecker_delta(M,1)) + &
!                            DVDL2MxztL_F0_io(J, (1-1)*3+1, (N-1)*3+2 )*DBLE(Kronecker_delta(M,2)) + &
!                            DVDL2MxztL_F0_io(J, (1-1)*3+1, (N-1)*3+3 )*DBLE(Kronecker_delta(M,3)) + &
!                            DVDL2MxztL_F0_io(J, (2-1)*3+2, (N-1)*3+1 )*DBLE(Kronecker_delta(M,1)) + &
!                            DVDL2MxztL_F0_io(J, (2-1)*3+2, (N-1)*3+2 )*DBLE(Kronecker_delta(M,2)) + &
!                            DVDL2MxztL_F0_io(J, (2-1)*3+2, (N-1)*3+3 )*DBLE(Kronecker_delta(M,3)) + &
!                            DVDL2MxztL_F0_io(J, (3-1)*3+3, (N-1)*3+1 )*DBLE(Kronecker_delta(M,1)) + &
!                            DVDL2MxztL_F0_io(J, (3-1)*3+3, (N-1)*3+2 )*DBLE(Kronecker_delta(M,2)) + &
!                            DVDL2MxztL_F0_io(J, (3-1)*3+3, (N-1)*3+3 )*DBLE(Kronecker_delta(M,3)) ) &
!                            )* CVISC                    
!                END DO
!            END DO
!        END DO
        
        
        DO J=1, NCL2
            DO M=1,NDV
                DO N=1,NDV
                    ANISOTROPY_FA(J,M,N) = uff2_FA(J,M,N)/( uff2_FA(J,1,1) + uff2_FA(J,2,2) + uff2_FA(J,3,3) ) -&
                                        DBLE(Kronecker_delta(M,N))/3.0_WP 
                END DO
            END DO
!            !WRITE(*,*) '  --  '
!            !WRITE(*,*) ANISOTROPY_FA(J,1,1:3)
!            !WRITE(*,*) ANISOTROPY_FA(J,2,1:3)
!            !WRITE(*,*) ANISOTROPY_FA(J,3,1:3)
!            !WRITE(*,*) ANISOTROPY_FA(J,1,1) + ANISOTROPY_FA(J,2,2) + ANISOTROPY_FA(J,3,3)
!            !WRITE(*,*) '  --  '
        END DO
        CALL CHKHDL('      ==>ANISOTROPY_FA',myid)
        
        DO J=1, NCL2
             ! below two methods are the same, checked good!
!            ANISTPinva_FA(J,1) = ANISOTROPY_FA(J,1,1)+ANISOTROPY_FA(J,2,2)+ANISOTROPY_FA(J,3,3)
!            ANISTPinva_FA(J,2) = -( ANISOTROPY_FA(J,1,1)*ANISOTROPY_FA(J,1,1) + &
!                                    ANISOTROPY_FA(J,1,2)*ANISOTROPY_FA(J,1,2) + &
!                                    ANISOTROPY_FA(J,1,3)*ANISOTROPY_FA(J,1,3) + &
!                                    ANISOTROPY_FA(J,2,1)*ANISOTROPY_FA(J,2,1) + &
!                                    ANISOTROPY_FA(J,2,2)*ANISOTROPY_FA(J,2,2) + &
!                                    ANISOTROPY_FA(J,2,3)*ANISOTROPY_FA(J,2,3) + &
!                                    ANISOTROPY_FA(J,3,1)*ANISOTROPY_FA(J,3,1) + &
!                                    ANISOTROPY_FA(J,3,2)*ANISOTROPY_FA(J,3,2) + &
!                                    ANISOTROPY_FA(J,3,3)*ANISOTROPY_FA(J,3,3) )/2.0_WP
!            ANISTPinva_FA(J,3) = (  ANISOTROPY_FA(J,1,1)*ANISOTROPY_FA(J,1,1)*ANISOTROPY_FA(J,1,1) + &
!                                    ANISOTROPY_FA(J,1,2)*ANISOTROPY_FA(J,2,1)*ANISOTROPY_FA(J,1,1) + &
!                                    ANISOTROPY_FA(J,1,3)*ANISOTROPY_FA(J,3,1)*ANISOTROPY_FA(J,1,1) + &
!                                    ANISOTROPY_FA(J,2,1)*ANISOTROPY_FA(J,1,1)*ANISOTROPY_FA(J,1,2) + &
!                                    ANISOTROPY_FA(J,2,2)*ANISOTROPY_FA(J,2,1)*ANISOTROPY_FA(J,1,2) + &
!                                    ANISOTROPY_FA(J,2,3)*ANISOTROPY_FA(J,3,1)*ANISOTROPY_FA(J,1,2) + &
!                                    ANISOTROPY_FA(J,3,1)*ANISOTROPY_FA(J,1,1)*ANISOTROPY_FA(J,1,3) + &
!                                    ANISOTROPY_FA(J,3,2)*ANISOTROPY_FA(J,2,1)*ANISOTROPY_FA(J,1,3) + &
!                                    ANISOTROPY_FA(J,3,3)*ANISOTROPY_FA(J,3,1)*ANISOTROPY_FA(J,1,3) + &
!                                    ANISOTROPY_FA(J,1,1)*ANISOTROPY_FA(J,1,2)*ANISOTROPY_FA(J,2,1) + &
!                                    ANISOTROPY_FA(J,1,2)*ANISOTROPY_FA(J,2,2)*ANISOTROPY_FA(J,2,1) + &
!                                    ANISOTROPY_FA(J,1,3)*ANISOTROPY_FA(J,3,2)*ANISOTROPY_FA(J,2,1) + &
!                                    ANISOTROPY_FA(J,2,1)*ANISOTROPY_FA(J,1,2)*ANISOTROPY_FA(J,2,2) + &
!                                    ANISOTROPY_FA(J,2,2)*ANISOTROPY_FA(J,2,2)*ANISOTROPY_FA(J,2,2) + &
!                                    ANISOTROPY_FA(J,2,3)*ANISOTROPY_FA(J,3,2)*ANISOTROPY_FA(J,2,2) + &
!                                    ANISOTROPY_FA(J,3,1)*ANISOTROPY_FA(J,1,2)*ANISOTROPY_FA(J,2,3) + &
!                                    ANISOTROPY_FA(J,3,2)*ANISOTROPY_FA(J,2,2)*ANISOTROPY_FA(J,2,3) + &
!                                    ANISOTROPY_FA(J,3,3)*ANISOTROPY_FA(J,3,2)*ANISOTROPY_FA(J,2,3) + &
!                                    ANISOTROPY_FA(J,1,1)*ANISOTROPY_FA(J,1,3)*ANISOTROPY_FA(J,3,1) + &
!                                    ANISOTROPY_FA(J,1,2)*ANISOTROPY_FA(J,2,3)*ANISOTROPY_FA(J,3,1) + &
!                                    ANISOTROPY_FA(J,1,3)*ANISOTROPY_FA(J,3,3)*ANISOTROPY_FA(J,3,1) + &
!                                    ANISOTROPY_FA(J,2,1)*ANISOTROPY_FA(J,1,3)*ANISOTROPY_FA(J,3,2) + &
!                                    ANISOTROPY_FA(J,2,2)*ANISOTROPY_FA(J,2,3)*ANISOTROPY_FA(J,3,2) + &
!                                    ANISOTROPY_FA(J,2,3)*ANISOTROPY_FA(J,3,3)*ANISOTROPY_FA(J,3,2) + &
!                                    ANISOTROPY_FA(J,3,1)*ANISOTROPY_FA(J,1,3)*ANISOTROPY_FA(J,3,3) + &
!                                    ANISOTROPY_FA(J,3,2)*ANISOTROPY_FA(J,2,3)*ANISOTROPY_FA(J,3,3) + &
!                                    ANISOTROPY_FA(J,3,3)*ANISOTROPY_FA(J,3,3)*ANISOTROPY_FA(J,3,3) )/3.0_WP
                                    
!            LumleyAxis_FA(J,1) = DSQRT(-ANISTPinva_FA(J,2)/3.0_WP)
!            LumleyAxis_FA(J,2) = (ANISTPinva_FA(J,3)/2.0_WP)**(1.0_wp/3.0_wp)
            
            Matrix(1:3,1:3)=ANISOTROPY_FA(J,1:3,1:3)
            call DSYEVC3(Matrix,eig)
            
            ANISTPinva_FA(J,1) = eig(1) + eig(2) + eig(3)
            ANISTPinva_FA(J,2) = (eig(1)*eig(1)+eig(1)*eig(2)+eig(2)*eig(2))*(-1.0_WP)
            ANISTPinva_FA(J,3) = eig(1)*eig(2)*(eig(1)+eig(2))*(-1.0_WP)
            LumleyAxis_FA(J,1) = DSQRT(DABS(-ANISTPinva_FA(J,2)/3.0_WP))            ! \eta
            LumleyAxis_FA(J,2) = SIGN( DABS(ANISTPinva_FA(J,3)/2.0_WP)**(1.0_wp/3.0_WP),  ANISTPinva_FA(J,3)/2.0_WP) ! \xi
            
!            !WRITE(*,*) 'invars', J, ANISTPinva_FA(J,1:3), LumleyAxis_FA(J,1:2)
!            !WRITE(*,*) 'eigenv', J, eig(1:3), -eig(1)-eig(2), eig(1)+eig(2)+eig(3)
            
        END DO
        CALL CHKHDL('      ==>LumleyAxis_FA',myid)
         
        ! calculate DrivenForce in streamwise direction for periodic directions.
        DO J=1, NCL2
            
            IF(J==1) THEN
                FF = ( ( YCL2ND_WFF(J+1)*(Tau_Mean_RA(J+1,1,2)-uff2d_FA(J+1,1,2))+ &
                         YCL2ND_WFB(J+1)*(Tau_Mean_RA(J,  1,2)-uff2d_FA(J,  1,2)) ) -  &
                      DABS(Tauw_io(1)) )*DYFI(J)
                      

            ELSE IF(J==NCL2) THEN
                FF = ( -DABS(Tauw_io(2)) -  &
                      ( YCL2ND_WFF(J)  *(Tau_Mean_RA(J,  1,2)-uff2d_FA(J,  1,2))+ &
                        YCL2ND_WFB(J)  *(Tau_Mean_RA(J-1,1,2)-uff2d_FA(J-1,1,2)) ) )*DYFI(J)

                        
            ELSE
                FF = ( ( YCL2ND_WFF(J+1)*(Tau_Mean_RA(J+1,1,2)-uff2d_FA(J+1,1,2))+ &
                         YCL2ND_WFB(J+1)*(Tau_Mean_RA(J,  1,2)-uff2d_FA(J,  1,2)) ) -  &
                       ( YCL2ND_WFF(J)  *(Tau_Mean_RA(J,  1,2)-uff2d_FA(J,  1,2))+ &
                         YCL2ND_WFB(J)  *(Tau_Mean_RA(J-1,1,2)-uff2d_FA(J-1,1,2)) ) )*DYFI(J)
                    
            END IF
                    
            DrivenForce(J) =  (FF+F_A*IBuoF(1)*D1xztL_F0_io(J))*(-1.0_WP)   
            !!WRITE(*,*) YCC(J), DrivenForce(J), FF,F_A*IBuoF(1)*D1xztL_F0_io(J)
        END DO
        
        DrivenFCTT1  = 0.0_WP
        DrivenFCTT2  = 0.0_WP
        DrivenFCTTU1 = 0.0_WP
        DrivenFCTTU2 = 0.0_WP
        DO J=1, NCL2
            DrivenFCTT1 = DrivenFCTT1 +DrivenForce(J)/DYFI(J)
            DrivenFCTTU1= DrivenFCTTU1+DrivenForce(J)*U_FA(J,1)/DYFI(J)
            
            DrivenFCTT2 = DrivenFCTT2 +FUxztL_F0_io(J,4)/DYFI(J)
            DrivenFCTTU2= DrivenFCTTU2+FUxztL_F0_io(J,4)*U_FA(J,1)/DYFI(J)
        END DO
        
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        OPEN(TECFLG,FILE=TRIM(filepath4)//'Result.IO.Table.WallandBulk.'//TRIM(PNTIM)//'.tec', POSITION='APPEND')
        WRITE(TECFLG,'(A)') '===================================='
        WRITE(TECFLG,'(A,2ES20.7)') 'Driven Force (Constant) & its MKE production=', DrivenFCTT2,DrivenFCTTU2
        WRITE(TECFLG,'(A,2ES20.7)') 'Driven Force (inversed) & its MKE production=', DrivenFCTT1,DrivenFCTTU1
        CLOSE(TECFLG)

        CALL CHKHDL('      ==>DrivenForce',myid)
        
        
        !=============================================================
        ddenintg = 0.0_WP
        denmintg = 0.0_WP
        bdfcintg= 0.0_WP
        densintg= 0.0_WP
        DO J=1, NCL2
            ! second order intgeral 
!            IF(J==1) THEN
!                DENtemp =0.5_WP*( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J) + YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) +             &
!                                  Dwal(1) )
!            ELSE IF(J==NCL2) THEN
!                DENtemp =0.5_WP*( Dwal(2) +             &
!                                  ( YCL2ND_WFF(J)  *D1xztL_F0_io(J) + YCL2ND_WFB(J)  *D1xztL_F0_io(J-1) ) )
!            ELSE
!                DENtemp =0.5_WP*( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J) + YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) +             &
!                                  ( YCL2ND_WFF(J)  *D1xztL_F0_io(J) + YCL2ND_WFB(J)  *D1xztL_F0_io(J-1) ) )
!            END IF
            !first order integral
            DENtemp = D1xztL_F0_io(J)
            
            denmintg = denmintg+DENtemp/DYFI(J)
            densintg(J) = denmintg !; !WRITE(*,*) J, DENtemp, densintg(J)
            
            ddenintg = ddenintg + (DENtemp-DenAvew)/DYFI(J)
            bdfcintg(J) = F_A*ddenintg
        END DO
        CALL CHKHDL('      ==>Calculated bodyforce distribution',myid)
        
        
    RETURN
    END SUBROUTINE    
    
!=================================================================================================================
    SUBROUTINE PP_BUDG_INIT
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        !==================Ruv=============================
        BUDG_prodc_stres_duiuj = 0.0_WP
        BUDG_viscs_dissp_duiuj = 0.0_WP
        BUDG_pdudx_stran_duiuj = 0.0_WP
        BUDG_Turbu_diffu_duiuj = 0.0_WP
        BUDG_dpudx_diffu_duiuj = 0.0_WP
        BUDG_viscs_diffu_duiuj = 0.0_WP
        
        BUDG_press_accl1_duiuj = 0.0_WP
        BUDG_viscs_accl1_duiuj = 0.0_WP
        BUDG_prodc_dvfc1_duiuj = 0.0_WP
        BUDG_balance1_duiuj    = 0.0_WP
        
        BUDG_prodc_gvfc2_duiuj = 0.0_WP
        BUDG_prodc_dvfc2_duiuj = 0.0_WP
        BUDG_turss_accl2_duiuj = 0.0_WP
        BUDG_balance2_duiuj    = 0.0_WP
        
        BUDG_pressure3_duiuj = 0.0_WP
        BUDG_vistress3_duiuj = 0.0_WP
        BUDG_balance3_duiuj  = 0.0_WP
        
        !==================Ruv==Sum along y=======================
        BUDG_prodc_stres_duiuj_ysum = 0.0_WP
        BUDG_viscs_dissp_duiuj_ysum = 0.0_WP
        BUDG_pdudx_stran_duiuj_ysum = 0.0_WP
        BUDG_Turbu_diffu_duiuj_ysum = 0.0_WP
        BUDG_dpudx_diffu_duiuj_ysum = 0.0_WP
        BUDG_viscs_diffu_duiuj_ysum = 0.0_WP
        
        BUDG_press_accl1_duiuj_ysum = 0.0_WP
        BUDG_viscs_accl1_duiuj_ysum = 0.0_WP
        BUDG_prodc_dvfc1_duiuj_ysum = 0.0_WP
        BUDG_balance1_duiuj_ysum = 0.0_WP
        
        BUDG_prodc_gvfc2_duiuj_ysum = 0.0_WP
        BUDG_prodc_dvfc2_duiuj_ysum = 0.0_WP
        BUDG_turss_accl2_duiuj_ysum = 0.0_WP
        BUDG_balance2_duiuj_ysum    = 0.0_WP
        
        BUDG_pressure3_duiuj_ysum = 0.0_WP
        BUDG_vistress3_duiuj_ysum = 0.0_WP
        BUDG_balance3_duiuj_ysum  = 0.0_WP
        
        !==================TKE=============================
        BUDG_prodc_stres_TKE = 0.0_WP
        BUDG_viscs_dissp_TKE = 0.0_WP
        BUDG_pdudx_stran_TKE = 0.0_WP
        BUDG_Turbu_diffu_TKE = 0.0_WP
        BUDG_dpudx_diffu_TKE = 0.0_WP
        BUDG_viscs_diffu_TKE = 0.0_WP
        
        BUDG_press_accl1_TKE = 0.0_WP
        BUDG_viscs_accl1_TKE = 0.0_WP
        BUDG_prodc_dvfc1_TKE = 0.0_WP
        BUDG_balance1_TKE    = 0.0_WP
        
        BUDG_prodc_gvfc2_TKE = 0.0_WP
        BUDG_prodc_dvfc2_TKE = 0.0_WP
        BUDG_turss_accl2_TKE = 0.0_WP
        BUDG_balance2_TKE    = 0.0_WP
        
        BUDG_pressure3_TKE = 0.0_WP
        BUDG_vistress3_TKE = 0.0_WP
        BUDG_balance3_TKE  = 0.0_WP
        
        !==================TKE sum along y=============================
        BUDG_prodc_stres_TKE_ysum = 0.0_WP
        BUDG_viscs_dissp_TKE_ysum = 0.0_WP
        BUDG_pdudx_stran_TKE_ysum = 0.0_WP
        BUDG_Turbu_diffu_TKE_ysum = 0.0_WP
        BUDG_dpudx_diffu_TKE_ysum = 0.0_WP
        BUDG_viscs_diffu_TKE_ysum = 0.0_WP
        
        BUDG_press_accl1_TKE_ysum = 0.0_WP
        BUDG_viscs_accl1_TKE_ysum = 0.0_WP
        BUDG_prodc_dvfc1_TKE_ysum = 0.0_WP
        BUDG_balance1_TKE_ysum    = 0.0_WP
        
        BUDG_prodc_gvfc2_TKE_ysum = 0.0_WP
        BUDG_prodc_dvfc2_TKE_ysum = 0.0_WP
        BUDG_turss_accl2_TKE_ysum = 0.0_WP
        BUDG_balance2_TKE_ysum    = 0.0_WP
        
        BUDG_pressure3_TKE_ysum = 0.0_WP
        BUDG_vistress3_TKE_ysum = 0.0_WP
        BUDG_balance3_TKE_ysum  = 0.0_WP
        
        !==================MKE=============================
        BUDG_prodc_stres_MKE = 0.0_WP
        BUDG_viscs_dissp_MKE = 0.0_WP
        BUDG_pdudx_stran_MKE = 0.0_WP
        BUDG_Turbu_diffu_MKE = 0.0_WP
        BUDG_dpudx_diffu_MKE = 0.0_WP
        BUDG_viscs_diffu_MKE = 0.0_WP
        
        BUDG_press_accl1_MKE = 0.0_WP
        BUDG_viscs_accl1_MKE = 0.0_WP
        BUDG_prodc_dvfc1_MKE = 0.0_WP
        BUDG_balance1_MKE    = 0.0_WP
        
        BUDG_prodc_gvfc1_MKE = 0.0_WP
        BUDG_prodc_dvfc2_MKE = 0.0_WP
        BUDG_turss_accl2_MKE = 0.0_WP
        BUDG_balance2_MKE    = 0.0_WP
        
        BUDG_pressure3_MKE = 0.0_WP
        BUDG_vistress3_MKE = 0.0_WP
        BUDG_balance3_MKE  = 0.0_WP
        
        !==================MKE sum along y=============================
        BUDG_prodc_stres_MKE_ysum = 0.0_WP
        BUDG_viscs_dissp_MKE_ysum = 0.0_WP
        BUDG_pdudx_stran_MKE_ysum = 0.0_WP
        BUDG_Turbu_diffu_MKE_ysum = 0.0_WP
        BUDG_dpudx_diffu_MKE_ysum = 0.0_WP
        BUDG_viscs_diffu_MKE_ysum = 0.0_WP
        
        BUDG_press_accl1_MKE_ysum = 0.0_WP
        BUDG_viscs_accl1_MKE_ysum = 0.0_WP
        BUDG_prodc_dvfc1_MKE_ysum = 0.0_WP
        BUDG_prodc_gvfc1_MKE_ysum = 0.0_WP
        BUDG_balance1_MKE_ysum = 0.0_WP
        
        BUDG_prodc_gvfc2_MKE_ysum = 0.0_WP
        BUDG_prodc_dvfc2_MKE_ysum = 0.0_WP
        BUDG_turss_accl2_MKE_ysum = 0.0_WP
        BUDG_balance2_MKE_ysum = 0.0_WP
        
        BUDG_pressure3_MKE_ysum = 0.0_WP
        BUDG_vistress3_MKE_ysum = 0.0_WP
        BUDG_balance3_MKE_ysum = 0.0_WP
        
        RETURN
    END SUBROUTINE 

!###########################  FA   ############################################################## 
    SUBROUTINE PP_FLOW_FA_RSTE_BUDG_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: J
        INTEGER(4) :: M, N, L, k
        INTEGER(4) :: TECFLG=200
        REAL(WP)   :: FF, TT
        REAL(WP)   :: BUDG_prod_Pij1, BUDG_prod_Pij2, coe
        REAL(WP)   :: visc_dissipation_duiduj1, visc_dissipation_duiduj2
        REAL(WP)   :: sum_p_related_mke
        character(128) :: FLNM
        logical        :: file_exists
        
        
        CAll CHKHDL('=====Calculating FA Budegts=====',myid)
        CALL PP_BUDG_INIT
        
!===FA====PRODUCTION TERMS ==P ==================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       Eq: BUDG_prodc_stres_duiuj(J,L) = P*_{ij}=P_{ij}-2/3P\delta_{ij}
!       Eq: P_{ij} = -<\rho u"_i u"_k> (\partial {u_j} / \partial x_k) + 
!                    -<\rho u"_j u"_k> (\partial {u_i} / \partial x_k)
        DO J=1,NCL2
            !==============for TKE and MKE==================
            BUDG_prodc_stres_MKE(J)=uff2d_FA(J,1,2)*dUdX_FA(J,1,2) + &
                                    uff2d_FA(J,2,2)*dUdX_FA(J,2,2) + &
                                    uff2d_FA(J,3,2)*dUdX_FA(J,3,2) 
            BUDG_prodc_stres_TKE(J)=-1.0_wp*BUDG_prodc_stres_MKE(J)
            
            !==========for each Ruv=========================
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    BUDG_prod_Pij1= uff2d_FA(J,M,2)*dUdX_FA(J,N,2)
                    BUDG_prod_Pij2= uff2d_FA(J,N,2)*dUdX_FA(J,M,2)
                    BUDG_prodc_stres_duiuj(J,L) = -1.0_wp*(BUDG_prod_Pij1 + BUDG_prod_Pij2) 
                              
                END DO
            END DO
            
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_prodc_stres',myid)
        
!==FA=====Viscous ENERGY dissipation term ===========================================================
!       TauDvDL_RA(J,M,N,H,P) = <\partial(u_m)/\partial x_n \tau_hp>
        DO J=1, NCL2    
            !==============for TKE and MKE==================
            BUDG_viscs_dissp_TKE(J) = &
                TauDvDL_RA(J,1,1,1,1) - Tau_Mean_RA(J,1,1) * DVDL1xztL_F0_io(J,1,1) + &
                TauDvDL_RA(J,2,1,1,2) - Tau_Mean_RA(J,1,2) * DVDL1xztL_F0_io(J,2,1) + &
                TauDvDL_RA(J,3,1,1,3) - Tau_Mean_RA(J,1,3) * DVDL1xztL_F0_io(J,3,1) + &
                TauDvDL_RA(J,1,2,2,1) - Tau_Mean_RA(J,2,1) * DVDL1xztL_F0_io(J,1,2) + &
                TauDvDL_RA(J,2,2,2,2) - Tau_Mean_RA(J,2,2) * DVDL1xztL_F0_io(J,2,2) + &
                TauDvDL_RA(J,3,2,2,3) - Tau_Mean_RA(J,2,3) * DVDL1xztL_F0_io(J,3,2) + &
                TauDvDL_RA(J,1,3,3,1) - Tau_Mean_RA(J,3,1) * DVDL1xztL_F0_io(J,1,3) + &
                TauDvDL_RA(J,2,3,3,2) - Tau_Mean_RA(J,3,2) * DVDL1xztL_F0_io(J,2,3) + &
                TauDvDL_RA(J,3,3,3,3) - Tau_Mean_RA(J,3,3) * DVDL1xztL_F0_io(J,3,3) 
            BUDG_viscs_dissp_TKE(J) = -1.0_wp*BUDG_viscs_dissp_TKE(J) 
            
            BUDG_viscs_dissp_MKE(J) = &  
                Tau_Mean_RA(J,1,1) * DVDL1xztL_F0_io(J,1,1) + &
                Tau_Mean_RA(J,2,1) * DVDL1xztL_F0_io(J,2,1) + &
                Tau_Mean_RA(J,3,1) * DVDL1xztL_F0_io(J,3,1) + &
                Tau_Mean_RA(J,1,2) * DVDL1xztL_F0_io(J,1,2) + &
                Tau_Mean_RA(J,2,2) * DVDL1xztL_F0_io(J,2,2) + &
                Tau_Mean_RA(J,3,2) * DVDL1xztL_F0_io(J,3,2) + &
                Tau_Mean_RA(J,1,3) * DVDL1xztL_F0_io(J,1,3) + &
                Tau_Mean_RA(J,2,3) * DVDL1xztL_F0_io(J,2,3) + &
                Tau_Mean_RA(J,3,3) * DVDL1xztL_F0_io(J,3,3)
            BUDG_viscs_dissp_MKE(J) = -1.0_wp*BUDG_viscs_dissp_MKE(J) 
                
            !==========for each Ruv=========================
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3  !Tau_ik_Du_jDx_i_RA(J,M,N) + &Tau_ik_Du_jDx_i_RA(J,N,M) +  &
                    visc_dissipation_duiduj1= &
                        TauDvDL_RA(J,N,1,M,1) + &
                        TauDvDL_RA(J,N,2,M,2) + &
                        TauDvDL_RA(J,N,3,M,3) - &
                        Tau_Mean_RA(J,M,1) * DVDL1xztL_F0_io(J,N,1) - &
                        Tau_Mean_RA(J,M,2) * DVDL1xztL_F0_io(J,N,2) - &
                        Tau_Mean_RA(J,M,3) * DVDL1xztL_F0_io(J,N,3)
                    visc_dissipation_duiduj2= &
                        TauDvDL_RA(J,M,1,N,1) + &
                        TauDvDL_RA(J,M,2,N,2) + &
                        TauDvDL_RA(J,M,3,N,3) - &
                        Tau_Mean_RA(J,N,1) * DVDL1xztL_F0_io(J,M,1) - &
                        Tau_Mean_RA(J,N,2) * DVDL1xztL_F0_io(J,M,2) - &
                        Tau_Mean_RA(J,N,3) * DVDL1xztL_F0_io(J,M,3)
                    BUDG_viscs_dissp_duiuj(J,L) = -1.0_wp*(visc_dissipation_duiduj1 +visc_dissipation_duiduj2)
                    ! do not need times density, as it goes into the equation. See Huang1995.
                    !!WRITE(*,*) L, visc_dissipation_duiduj1, visc_dissipation_duiduj2
                END DO
            END DO
            
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_viscs_dissp_duiuj',myid) 
        
!==FA=====Velocity-Pressure gradient Pressure strain term=======(Strain )=======>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!       Eq = <p' \partial { u"_i}/\partial{x_j} > +
!            <p' \partial { u"_j}/\partial{x_i} >
        DO J=1,NCL2
            !==============for TKE and MKE==================
            BUDG_pdudx_stran_TKE(J) =   DVDLPxztL_F0_io(J,1,1) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,1,1) + &
                                        DVDLPxztL_F0_io(J,2,2) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,2,2) + &
                                        DVDLPxztL_F0_io(J,3,3) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,3,3)
            BUDG_pdudx_stran_MKE(J) =   U1xztL_F0_io(J,4)* dUidXi(J)
        
            !==========for each Ruv=========================
            DO M=1,NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    BUDG_pdudx_stran_duiuj(J,L) = DVDLPxztL_F0_io(J,M,N) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,M,N) + &
                                                  DVDLPxztL_F0_io(J,N,M) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,N,M)
                END DO
            END DO
            !!WRITE(*,*) J, YCC(J), DVDLPxztL_F0_io(J,1,1),U1xztL_F0_io(J,4), DVDL1xztL_F0_io(J,1,1), BUDG_pdudx_stran_duiuj(J,1) !test
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_pdudx_stran_duiuj',myid)
        
        
!==FA=====TURBUELCEN DIFFUSION TERMS=( Turb. Diffusion)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       Eq.  BUDG_Tdiff_duiuj(J,L) = ( \partial <\rho u"_i u"_j u"_k > ) / (\partial x_k )
        DO J=1,NCL2   
            uffMKEffd_FA(J)= U_FA(J,1)*uff2d_FA(J,1,2) + &
                             U_FA(J,2)*uff2d_FA(J,2,2) + &
                             U_FA(J,3)*uff2d_FA(J,3,2)
                                
            uffTKEffd_FA(J) = 0.5_WP*( uff3d_FA(J,1,1,2) + uff3d_FA(J,2,2,2) + uff3d_FA(J,3,3,2) )
        END DO
        
        DO J=1,NCL2 
            !==============for TKE and MKE==================
            IF(J==1) THEN
                BUDG_Turbu_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*uffTKEffd_FA(J+1) + &
                        YCL2ND_WFB(J+1)*uffTKEffd_FA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_wp)
                BUDG_Turbu_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*uffMKEffd_FA(J+1) + &
                        YCL2ND_WFB(J+1)*uffMKEffd_FA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_wp)
            ELSE IF (J==NCL2) THEN
                BUDG_Turbu_diffu_TKE(J)= &
                    ( 0.0_WP - &
                      ( YCL2ND_WFF(J)*uffTKEffd_FA(J  ) + &
                        YCL2ND_WFB(J)*uffTKEffd_FA(J-1) ) ) * DYFI(J) * (-1.0_wp)
                BUDG_Turbu_diffu_MKE(J)= &
                    ( 0.0_WP - &
                      ( YCL2ND_WFF(J)*uffMKEffd_FA(J  ) + &
                        YCL2ND_WFB(J)*uffMKEffd_FA(J-1) ) ) * DYFI(J) * (-1.0_wp)
                                  
            ELSE
                BUDG_Turbu_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*uffTKEffd_FA(J+1) + &
                        YCL2ND_WFB(J+1)*uffTKEffd_FA(J  ) ) -         &
                      ( YCL2ND_WFF(J)  *uffTKEffd_FA(J  ) + &
                        YCL2ND_WFB(J)  *uffTKEffd_FA(J-1) ) ) * DYFI(J) * (-1.0_wp)
                BUDG_Turbu_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*uffMKEffd_FA(J+1) + &
                        YCL2ND_WFB(J+1)*uffMKEffd_FA(J  ) ) -         &
                      ( YCL2ND_WFF(J)  *uffMKEffd_FA(J  ) + &
                        YCL2ND_WFB(J)  *uffMKEffd_FA(J-1) ) ) * DYFI(J) * (-1.0_wp)
            END IF
            
            !==========for each Ruv=========================
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    IF(J==1) THEN
                        BUDG_Turbu_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*uff3d_FA(J+1,M,N,2) + &
                                YCL2ND_WFB(J+1)*uff3d_FA(J,  M,N,2) ) - 0.0_WP ) * DYFI(J) * (-1.0_wp)
                    ELSE IF (J==NCL2) THEN
                        BUDG_Turbu_diffu_duiuj(J,L)= &
                            ( 0.0_WP - &
                              ( YCL2ND_WFF(J)*uff3d_FA(J,  M,N,2) + &
                                YCL2ND_WFB(J)*uff3d_FA(J-1,M,N,2) ) ) * DYFI(J) * (-1.0_wp)
                                          
                    ELSE
                        BUDG_Turbu_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*uff3d_FA(J+1,M,N,2) + &
                                YCL2ND_WFB(J+1)*uff3d_FA(J,  M,N,2) ) -         &
                              ( YCL2ND_WFF(J)  *uff3d_FA(J,  M,N,2) + &
                                YCL2ND_WFB(J)  *uff3d_FA(J-1,M,N,2) ) ) * DYFI(J) * (-1.0_wp)
                    END IF
                    
                    
!                    !=======test==checked the same as above==================
!                    FF=BUDG_Turbu_diffu_duiuj(J,L)
!                    IF(J==1) THEN
!                        BUDG_Turbu_diffu_duiuj(J,L)= &
!                            ( ( YCL2ND_WFF(J+1)*TDIFU_FA(J+1,M,N) + &
!                                YCL2ND_WFB(J+1)*TDIFU_FA(J,  M,N) ) - 0.0_WP ) * DYFI(J) * (-1.0_wp)
!                    ELSE IF (J==NCL2) THEN
!                        BUDG_Turbu_diffu_duiuj(J,L)= &
!                            ( 0.0_WP - &
!                              ( YCL2ND_WFF(J)*TDIFU_FA(J,  M,N) + &
!                                YCL2ND_WFB(J)*TDIFU_FA(J-1,M,N) ) ) * DYFI(J) * (-1.0_wp)
                                          
!                    ELSE
!                        BUDG_Turbu_diffu_duiuj(J,L)= &
!                            ( ( YCL2ND_WFF(J+1)*TDIFU_FA(J+1,M,N) + &
!                                YCL2ND_WFB(J+1)*TDIFU_FA(J,  M,N) ) -         &
!                              ( YCL2ND_WFF(J)  *TDIFU_FA(J,  M,N) + &
!                                YCL2ND_WFB(J)  *TDIFU_FA(J-1,M,N) ) ) * DYFI(J) * (-1.0_wp)
!                    END IF
!                    !WRITE(*,*) M, N, J, BUDG_Turbu_diffu_duiuj(J,L)-FF
!                    !========test
                    
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_Turbu_diffu_duiuj',myid)
        
            
!==FA=====Velocity-Pressure gradient diffusion term=======(pressure diffusion)==================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!       Eq. = - \partial (<p' u"_j>) /\partial (x_i) - \partial (<p' u"_i>) /\partial (x_j)
!           = - \partial (<p' u'_j>) /\partial (x_i) - \partial (<p' u'_i>) /\partial (x_j)
        DO J=1,NCL2
            !==============for TKE and MKE==================
            IF(J==1) THEN
                BUDG_dpudx_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufpf_RA(J+1,2) + &
                        YCL2ND_WFB(J+1)*ufpf_RA(J,  2) ) &
                         -  0.0_WP )*DYFI(J)* (-1.0_wp)
                         
                BUDG_dpudx_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*(U1xztL_F0_io(J+1,2)*U1xztL_F0_io(J+1,4)) + &
                        YCL2ND_WFB(J+1)*(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) ) &
                         -  0.0_WP )*DYFI(J)* (-1.0_wp)
                         
            ELSE IF(J==NCL2) THEN
                BUDG_dpudx_diffu_TKE(J)= &
                    ( 0.0_WP -  &
                    ( YCL2ND_WFF(J)  *ufpf_RA(J,  2) + &
                      YCL2ND_WFB(J)  *ufpf_RA(J-1,2) ) )*DYFI(J)* (-1.0_wp)
                      
                BUDG_dpudx_diffu_MKE(J)= &
                    ( 0.0_WP -  &
                    ( YCL2ND_WFF(J)  *(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) + &
                      YCL2ND_WFB(J)  *(U1xztL_F0_io(J-1,2)*U1xztL_F0_io(J-1,4)) ))*DYFI(J)* (-1.0_wp)
            ELSE
            
                BUDG_dpudx_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufpf_RA(J+1,2) + &
                        YCL2ND_WFB(J+1)*ufpf_RA(J,  2) ) -  &
                      ( YCL2ND_WFF(J)  *ufpf_RA(J,  2) + &
                        YCL2ND_WFB(J)  *ufpf_RA(J-1,2) ) )*DYFI(J)* (-1.0_wp)
                                    
                BUDG_dpudx_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*(U1xztL_F0_io(J+1,2)*U1xztL_F0_io(J+1,4)) + &
                        YCL2ND_WFB(J+1)*(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) ) -  &
                      ( YCL2ND_WFF(J)  *(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) + &
                        YCL2ND_WFB(J)  *(U1xztL_F0_io(J-1,2)*U1xztL_F0_io(J-1,4)) ) )*DYFI(J)* (-1.0_wp)
            END IF
        
        
            !==========for each Ruv=========================
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    IF(J==1) THEN
                        BUDG_dpudx_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*( ufpf_RA(J+1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J+1,N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J+1)*( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   ) )- &
                                0.0_WP )*DYFI(J)
                    ELSE IF(J==NCL2) THEN
                        BUDG_dpudx_diffu_duiuj(J,L)= &
                            (  0.0_WP- &
                              ( YCL2ND_WFF(J)  *( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J)  *( ufpf_RA(J-1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J-1,N)*DBLE(Kronecker_delta(M,2))   ) ) )*DYFI(J)
                    ELSE
                        BUDG_dpudx_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*( ufpf_RA(J+1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J+1,N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J+1)*( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   ) )- &
                              ( YCL2ND_WFF(J)  *( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J)  *( ufpf_RA(J-1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J-1,N)*DBLE(Kronecker_delta(M,2))   ) ) )*DYFI(J)
                    END IF
                    
                    
                    BUDG_dpudx_diffu_duiuj(J,L) = BUDG_dpudx_diffu_duiuj(J,L)* (-1.0_wp)
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_dpudx_diffu_duiuj',myid)
        

!===FA====Viscous diffusion term ===========================================================
!       Viscous stress term is based on RA decomposition like Huang, not FA. 
!       Eq. = \partial <u"_j tau'_ki> / partial (x_k) + 
!             \partial <u"_i tau'_kj> / partial (x_k)
        DO J=1, NCL2
            !==============for TKE and MKE==================
            IF(J==1) THEN
                BUDG_viscs_diffu_TKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J+1)* ( Taufuf_RA(J+1,1, 2, 1) + &
                                         Taufuf_RA(J+1,2, 2, 2) + &
                                         Taufuf_RA(J+1,3, 2, 3) ) )-&
                    0.0_WP &
                    )*DYFI(J)
                BUDG_viscs_diffu_MKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,  1,2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,  2,2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,  3,2) ) &
                     +YCL2ND_WFB(J+1)* ( U1xztL_F0_io(J+1,1)*Tau_Mean_RA(J+1,1,2) + &
                                         U1xztL_F0_io(J+1,2)*Tau_Mean_RA(J+1,2,2) + &
                                         U1xztL_F0_io(J+1,3)*Tau_Mean_RA(J+1,3,2) ) )-&
                    0.0_WP &
                    )*DYFI(J)
            ELSE IF (J==NCL2) THEN
                BUDG_viscs_diffu_TKE(J) =  ( &
                    0.0_WP-&
                    ( YCL2ND_WFF(J)  * ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J)  * ( Taufuf_RA(J-1,1, 2, 1) + &
                                         Taufuf_RA(J-1,2, 2, 2) + &
                                         Taufuf_RA(J-1,3, 2, 3) ) ) &
                    )*DYFI(J)
                    
                BUDG_viscs_diffu_MKE(J) =  ( &
                    0.0_WP-&
                    ( YCL2ND_WFF(J)  * ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,1,  2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,2,  2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,3,  2) ) &
                     +YCL2ND_WFB(J)  * ( U1xztL_F0_io(J-1,1)*Tau_Mean_RA(J-1,1,2) + &
                                         U1xztL_F0_io(J-1,2)*Tau_Mean_RA(J-1,2,2) + &
                                         U1xztL_F0_io(J-1,3)*Tau_Mean_RA(J-1,3,2) ) ) &
                    )*DYFI(J)
            ELSE
                BUDG_viscs_diffu_TKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J+1)* ( Taufuf_RA(J+1,1, 2, 1) + &
                                         Taufuf_RA(J+1,2, 2, 2) + &
                                         Taufuf_RA(J+1,3, 2, 3) ) )-&
                    ( YCL2ND_WFF(J)  * ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J)  * ( Taufuf_RA(J-1,1, 2, 1) + &
                                         Taufuf_RA(J-1,2, 2, 2) + &
                                         Taufuf_RA(J-1,3, 2, 3) ) ) &
                    )*DYFI(J)
                    
                BUDG_viscs_diffu_MKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,  1,2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,  2,2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,  3,2) ) &
                     +YCL2ND_WFB(J+1)* ( U1xztL_F0_io(J+1,1)*Tau_Mean_RA(J+1,1,2) + &
                                         U1xztL_F0_io(J+1,2)*Tau_Mean_RA(J+1,2,2) + &
                                         U1xztL_F0_io(J+1,3)*Tau_Mean_RA(J+1,3,2) ) )-&
                    ( YCL2ND_WFF(J)  * ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,  1,2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,  2,2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,  3,2) ) &
                     +YCL2ND_WFB(J)  * ( U1xztL_F0_io(J-1,1)*Tau_Mean_RA(J-1,1,2) + &
                                         U1xztL_F0_io(J-1,2)*Tau_Mean_RA(J-1,2,2) + &
                                         U1xztL_F0_io(J-1,3)*Tau_Mean_RA(J-1,3,2) ) ) &
                    )*DYFI(J)
            END IF
        
            !==========for each Ruv=========================
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    IF(J==1) THEN
                        BUDG_viscs_diffu_duiuj(J,L) =  ( &
                            ( YCL2ND_WFF(J+1)* Taufuf_RA(J,  2, M, N)   &
                             +YCL2ND_WFB(J+1)* Taufuf_RA(J+1,2, M, N) )-&
                            0.0_WP &
                            )*DYFI(J) + (&
                            ( YCL2ND_WFF(J+1)* Taufuf_RA(J,  2, N, M)   &
                             +YCL2ND_WFB(J+1)* Taufuf_RA(J+1,2, N, M) )-&
                            0.0_WP &
                            )*DYFI(J)
                    ELSE IF (J==NCL2) THEN
                        BUDG_viscs_diffu_duiuj(J,L) =  ( &
                            0.0_WP-&
                            ( YCL2ND_WFF(J)  * Taufuf_RA(J,  2, M, N)   &
                             +YCL2ND_WFB(J)  * Taufuf_RA(J-1,2, M, N) ) &
                            )*DYFI(J) + (&
                            0.0_WP-&
                            ( YCL2ND_WFF(J)  * Taufuf_RA(J,  2, N, M)   &
                             +YCL2ND_WFB(J)  * Taufuf_RA(J-1,2, N, M) ) &
                            )*DYFI(J)
                    ELSE
                        BUDG_viscs_diffu_duiuj(J,L) =  ( &
                            ( YCL2ND_WFF(J+1)* Taufuf_RA(J,  2, M, N)   &
                             +YCL2ND_WFB(J+1)* Taufuf_RA(J+1,2, M, N) )-&
                            ( YCL2ND_WFF(J)  * Taufuf_RA(J,  2, M, N)   &
                             +YCL2ND_WFB(J)  * Taufuf_RA(J-1,2, M, N) ) &
                            )*DYFI(J) + (&
                            ( YCL2ND_WFF(J+1)* Taufuf_RA(J,  2, N, M)   &
                             +YCL2ND_WFB(J+1)* Taufuf_RA(J+1,2, N, M) )-&
                            ( YCL2ND_WFF(J)  * Taufuf_RA(J,  2, N, M)   &
                             +YCL2ND_WFB(J)  * Taufuf_RA(J-1,2, N, M) ) &
                            )*DYFI(J)
                    END IF
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_viscs_diffu_duiuj',myid) 
        
!====FA==========below 1 first decomposition method ===========================================
        !==FA=====Pressure acceleration term, including all========(FA)===========>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        DO J=1,NCL2
            !==============for TKE and MKE==================
            BUDG_press_accl1_MKE(J) = uff_RA(J,1)*dPdX_RA(J,1) + uff_RA(J,2)*dPdX_RA(J,2) + uff_RA(J,3)*dPdX_RA(J,3)
            BUDG_press_accl1_TKE(J) = -1.0_WP*BUDG_press_accl1_MKE(J)
            
            !==========for each Ruv=========================
            DO M=1,NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    BUDG_press_accl1_duiuj(J,L) = -dPdX_RA(J,M)*uff_RA(J,N) &
                                                  -dPdX_RA(J,N)*uff_RA(J,M)  
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_press_accl1_duiuj',myid)      
                      
        !==FA=====Viscous acceleration, including all========(based on tau_RA)======>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       Eq. <u"_j> (\partial <tau_ki>)/ (\partial  x_k) +
!           <u"_i> (\partial <tau_kj>)/ (\partial  x_k)
        DO J=1, NCL2
            !==============for TKE and MKE==================
            BUDG_viscs_accl1_TKE(J) = dTaudy_RA(J,2,1)*uff_RA(J,1) + &
                                      dTaudy_RA(J,2,2)*uff_RA(J,2) + &
                                      dTaudy_RA(J,2,3)*uff_RA(J,3)
            BUDG_viscs_accl1_MKE(J) = -1.0_WP * BUDG_viscs_accl1_TKE(J)
            
            !==========for each Ruv=========================
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    BUDG_viscs_accl1_duiuj(J,L) = dTaudy_RA(J,2,M)*uff_RA(J,N) + &
                                                  dTaudy_RA(J,2,N)*uff_RA(J,M)
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_viscs_accl1_duiuj',myid)    
        
        !==FA=====Production by driven force method 1======>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        BUDG_prodc_dvfc1_duiuj=0.0_wp  
        DO J=1,NCL2
        
            !============exact method============================
            BUDG_prodc_dvfc1_duiuj(J,1) = ( FUxztL_F0_io(J,1)-FUxztL_F0_io(J,4)*U_FA(J,1)  )*2.0_wp
            BUDG_prodc_dvfc1_duiuj(J,2) =   FUxztL_F0_io(J,2)-FUxztL_F0_io(J,4)*U_FA(J,2)
            BUDG_prodc_dvfc1_duiuj(J,3) =   FUxztL_F0_io(J,3)-FUxztL_F0_io(J,4)*U_FA(J,3)
            BUDG_prodc_dvfc1_TKE(J)     = FUxztL_F0_io(J,1) - FUxztL_F0_io(J,4)*U_FA(J,1)
!            BUDG_prodc_dvfc1_MKE(J)    = FUxztL_F0_io(J,4) * U_FA(J,1)
            
            !==============rough method==================
!            BUDG_prodc_dvfc1_duiuj(J,1) =   DrivenForce(J)*uff_RA(J,1) *2.0_wp
!            BUDG_prodc_dvfc1_duiuj(J,2) =   DrivenForce(J)*uff_RA(J,2) 
!            BUDG_prodc_dvfc1_duiuj(J,3) =   DrivenForce(J)*uff_RA(J,3) 
!            BUDG_prodc_dvfc1_TKE(J)     =   DrivenForce(J)*uff_RA(J,1) 
            BUDG_prodc_dvfc1_MKE(J)     =   DrivenForce(J)*U_FA(J,1)   ! Good
            
        END DO
        
        !=========FA======method 1 decomposition, no gravity production for turb. =================
        DO J=1, NCL2
            BUDG_prodc_gvfc1_MKE(J) = F_A*D1xztL_F0_io(J)*U1xztL_F0_io(J,1)
        END DO
        
!====FA==========below 2 second decomposition method ===========================================
        !==FA=====Production by driven force method 2==
        DO J=1,NCL2
            BUDG_prodc_dvfc2_duiuj(J,1) = ( FUxztL_F0_io(J,1)-FUxztL_F0_io(J,4)*U1xztL_F0_io(J,1)  )*2.0_wp
            BUDG_prodc_dvfc2_duiuj(J,2) =   FUxztL_F0_io(J,2)-FUxztL_F0_io(J,4)*U1xztL_F0_io(J,2)
            BUDG_prodc_dvfc2_duiuj(J,3) =   FUxztL_F0_io(J,3)-FUxztL_F0_io(J,4)*U1xztL_F0_io(J,3)
            BUDG_prodc_dvfc2_TKE(J)     =   FUxztL_F0_io(J,1)-FUxztL_F0_io(J,4)*U1xztL_F0_io(J,1)
            BUDG_prodc_dvfc2_MKE(J)     =   DrivenForce(J)*U1xztL_F0_io(J,1) 
        END DO
        
        !==FA=====Production by gravity force method 2==
        BUDG_prodc_gvfc2_duiuj=0.0_wp  
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    IF(M.NE.GRAVDIR .AND. N.NE.GRAVDIR) CYCLE
                    
                    L = (M*(7-M))/2+N-3 
                       
                    IF(M.EQ.GRAVDIR .AND. N.EQ.GRAVDIR) THEN
                        COE = 2.0_WP 
                        K   = GRAVDIR    
                    ELSE IF (M.EQ.GRAVDIR) THEN
                        COE = 1.0_WP
                        K   = N
                    ELSE IF (N.EQ.GRAVDIR) THEN
                        COE = 1.0_WP
                        K   = M
                    ELSE
                        COE = 0.0_WP
                        K   = -1 !(WHICH WILL LEAD TO ERROR!)
                    END IF
                    ! F_A includes a postive or negtive sign
                    BUDG_prodc_gvfc2_duiuj(J,L)= F_A * COE * ( G1xztL_F0_io(J,K) - D1xztL_F0_io(J) * U1xztL_F0_io(J,K) )
                    
                    !du=( G1xztL_F0_io(J,1) - D1xztL_F0_io(J) * U1xztL_F0_io(J,1) )*F_A
                    !dv=( G1xztL_F0_io(J,2) - D1xztL_F0_io(J) * U1xztL_F0_io(J,2) )*F_A
                    !dw=( G1xztL_F0_io(J,3) - D1xztL_F0_io(J) * U1xztL_F0_io(J,3) )*F_A
            
                END DO
            END DO
            !==============for TKE and MKE==================
            BUDG_prodc_gvfc2_TKE(J) = 0.5_WP*(BUDG_prodc_gvfc2_duiuj(J,1)+BUDG_prodc_gvfc2_duiuj(J,4)+BUDG_prodc_gvfc2_duiuj(J,6))
            
            BUDG_prodc_gvfc2_MKE(J) = F_A*D1xztL_F0_io(J)*U_FA(J,1)
        END DO
        
        !==FA=====Production by gravity force method 2==
        BUDG_turss_accl2_duiuj=0.0_wp  
        DO J =1, NCL2
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    BUDG_turss_accl2_duiuj(J,L) = uff_RA(J,M) * dTSSdy_RA(J,2,N) + uff_RA(J,N) * dTSSdy_RA(J,2,M)
                END DO
            END DO
        
            BUDG_turss_accl2_TKE(J) =   uff_RA(J,1) * dTSSdy_RA(J,2,1) + &
                                        uff_RA(J,2) * dTSSdy_RA(J,2,2) + &
                                        uff_RA(J,3) * dTSSdy_RA(J,2,3)
            BUDG_turss_accl2_MKE(J) = -1.0_WP * BUDG_turss_accl2_TKE(J)
        END DO
        
!====FA==========below is the 3 third  sum of some terms===========================
        !========pressure related===================
        DO J=1, NCL2
            BUDG_pressure3_MKE(J) = -U_FA(J,2)*dPdX_RA(J,2)
            FF = 0.0_WP
            IF(J==1) THEN
                FF =  ( &
                    ( YCL2ND_WFF(J+1)* UPxztL_F0_io(J,  2)   &
                     +YCL2ND_WFB(J+1)* UPxztL_F0_io(J+1,2) )-&
                    0.0_WP)*DYFI(J) 
                    
            ELSE IF (J==NCL2) THEN
                FF =  ( &
                    0.0_WP-&
                    ( YCL2ND_WFF(J)  * UPxztL_F0_io(J,  2)   &
                     +YCL2ND_WFB(J)  * UPxztL_F0_io(J-1,2) ) &
                    )*DYFI(J)
            ELSE
                FF =  ( &
                    ( YCL2ND_WFF(J+1)* UPxztL_F0_io(J,  2)   &
                     +YCL2ND_WFB(J+1)* UPxztL_F0_io(J+1,2) )-&
                    ( YCL2ND_WFF(J)  * UPxztL_F0_io(J,  2)   &
                     +YCL2ND_WFB(J)  * UPxztL_F0_io(J-1,2) ) &
                    )*DYFI(J)
            END IF
            BUDG_pressure3_TKE(J) = FF - DVDLPxztL_F0_io(J,1,1)- DVDLPxztL_F0_io(J,2,2)- DVDLPxztL_F0_io(J,3,3) &
                                      - U_FA(J,2)*dPdX_RA(J,2)
            BUDG_pressure3_TKE(J) = -1.0_wp * BUDG_pressure3_TKE(J)
            
            
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    FF =0.0_WP
                    IF(N==2) THEN
                        IF(J==1) THEN
                            FF =  FF+ ( &
                                ( YCL2ND_WFF(J+1)* UPxztL_F0_io(J,  M)   &
                                 +YCL2ND_WFB(J+1)* UPxztL_F0_io(J+1,M) )-&
                                0.0_WP &
                                )*DYFI(J)
                        ELSE IF (J==NCL2) THEN
                            FF =  FF+ ( &
                                0.0_WP-&
                                ( YCL2ND_WFF(J)  * UPxztL_F0_io(J,  M)   &
                                 +YCL2ND_WFB(J)  * UPxztL_F0_io(J-1,M) ) &
                                )*DYFI(J)
                        ELSE
                            FF =  FF+ ( &
                                ( YCL2ND_WFF(J+1)* UPxztL_F0_io(J,  M)   &
                                 +YCL2ND_WFB(J+1)* UPxztL_F0_io(J+1,M) )-&
                                ( YCL2ND_WFF(J)  * UPxztL_F0_io(J,  M)   &
                                 +YCL2ND_WFB(J)  * UPxztL_F0_io(J-1,M) ) &
                                )*DYFI(J)
                        END IF
                    END IF
                    
                    IF(M==2) THEN
                        IF(J==1) THEN
                            FF = FF+  ( &
                                ( YCL2ND_WFF(J+1)* UPxztL_F0_io(J,  N)   &
                                 +YCL2ND_WFB(J+1)* UPxztL_F0_io(J+1,N) )-&
                                0.0_WP &
                                )*DYFI(J)
                        ELSE IF (J==NCL2) THEN
                            FF =  FF+ ( &
                                0.0_WP-&
                                ( YCL2ND_WFF(J)  * UPxztL_F0_io(J,  N)   &
                                 +YCL2ND_WFB(J)  * UPxztL_F0_io(J-1,N) ) &
                                )*DYFI(J)
                        ELSE
                            FF = FF+  ( &
                                ( YCL2ND_WFF(J+1)* UPxztL_F0_io(J,  N)   &
                                 +YCL2ND_WFB(J+1)* UPxztL_F0_io(J+1,N) )-&
                                ( YCL2ND_WFF(J)  * UPxztL_F0_io(J,  N)   &
                                 +YCL2ND_WFB(J)  * UPxztL_F0_io(J-1,N) ) &
                                )*DYFI(J)
                        END IF
                    END IF
                    
                    BUDG_pressure3_duiuj(J,L) = FF - DVDLPxztL_F0_io(J,M,N) - DVDLPxztL_F0_io(J,N,M) &
                                                  - U_FA(J,M)*dPdX_RA(J,N) - U_FA(J,N)*dPdX_RA(J,M)
                    BUDG_pressure3_duiuj(J,L) = -1.0_WP * BUDG_pressure3_duiuj(J,L)
                END DO
            END DO
        
        END DO
        
        
        !========stress related===================
        DO J=1, NCL2
            BUDG_vistress3_MKE(J) = U_FA(J,1)*dTaudy_RA(J,2,1) + &
                                   U_FA(J,2)*dTaudy_RA(J,2,2) + &
                                   U_FA(J,3)*dTaudy_RA(J,2,3)
            FF = 0.0_WP
            IF(J==1) THEN
                FF =  ( &
                    ( YCL2ND_WFF(J+1)* (TauU_RA(J,  1,2,1)+TauU_RA(J,  2,2,2)+TauU_RA(J,  3,2,3))   &
                     +YCL2ND_WFB(J+1)* (TauU_RA(J+1,1,2,1)+TauU_RA(J+1,2,2,2)+TauU_RA(J+1,3,2,3)) )-&
                    0.0_WP &
                    )*DYFI(J)
            ELSE IF (J==NCL2) THEN
                FF =  ( &
                    0.0_WP-&
                    ( YCL2ND_WFF(J)  * (TauU_RA(J,  1,2,1)+TauU_RA(J,  2,2,2)+TauU_RA(J,  3,2,3))   &
                     +YCL2ND_WFB(J)  * (TauU_RA(J-1,1,2,1)+TauU_RA(J-1,2,2,2)+TauU_RA(J-1,3,2,3)) ) &
                    )*DYFI(J)
            ELSE
                FF =  ( &
                    ( YCL2ND_WFF(J+1)* (TauU_RA(J,  1,2,1)+TauU_RA(J,  2,2,2)+TauU_RA(J,  3,2,3))   &
                     +YCL2ND_WFB(J+1)* (TauU_RA(J+1,1,2,1)+TauU_RA(J+1,2,2,2)+TauU_RA(J+1,3,2,3)) )-&
                    ( YCL2ND_WFF(J)  * (TauU_RA(J,  1,2,1)+TauU_RA(J,  2,2,2)+TauU_RA(J,  3,2,3))   &
                     +YCL2ND_WFB(J)  * (TauU_RA(J-1,1,2,1)+TauU_RA(J-1,2,2,2)+TauU_RA(J-1,3,2,3)) ) &
                    )*DYFI(J)
            END IF
            BUDG_vistress3_TKE(J) = FF &
                                    -TauDvDL_RA(J,1,1,1,1) - TauDvDL_RA(J,2,1,2,1) - TauDvDL_RA(J,3,1,3,1) &
                                    -TauDvDL_RA(J,1,2,1,2) - TauDvDL_RA(J,2,2,2,2) - TauDvDL_RA(J,3,2,3,2) &
                                    -TauDvDL_RA(J,1,3,1,3) - TauDvDL_RA(J,2,3,2,3) - TauDvDL_RA(J,3,3,3,3) &
                                    -U_FA(J,1)*dTaudy_RA(J,2,1) &
                                    -U_FA(J,2)*dTaudy_RA(J,2,2) &
                                    -U_FA(J,3)*dTaudy_RA(J,2,3)
                                      
                                      
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    FF =0.0_WP

                    IF(J==1) THEN
                        FF =  ( &
                            ( YCL2ND_WFF(J+1)* TauU_RA(J,  M,2,N)   &
                             +YCL2ND_WFB(J+1)* TauU_RA(J+1,M,2,N) )-&
                            0.0_WP &
                            )*DYFI(J) + ( &
                            ( YCL2ND_WFF(J+1)* TauU_RA(J,  N,2,M)   &
                             +YCL2ND_WFB(J+1)* TauU_RA(J+1,N,2,M) )-&
                            0.0_WP &
                            )*DYFI(J)
                    ELSE IF (J==NCL2) THEN
                        FF =  ( &
                            0.0_WP-&
                            ( YCL2ND_WFF(J)  * TauU_RA(J,  M,2,N)   &
                             +YCL2ND_WFB(J)  * TauU_RA(J-1,M,2,N) ) &
                            )*DYFI(J)+ ( &
                            0.0_WP-&
                            ( YCL2ND_WFF(J)  * TauU_RA(J,  N,2,M)   &
                             +YCL2ND_WFB(J)  * TauU_RA(J-1,N,2,M) ) &
                            )*DYFI(J)
                    ELSE
                        FF =  ( &
                            ( YCL2ND_WFF(J+1)* TauU_RA(J,  M,2,N)   &
                             +YCL2ND_WFB(J+1)* TauU_RA(J+1,M,2,N) )-&
                            ( YCL2ND_WFF(J)  * TauU_RA(J,  M,2,N)   &
                             +YCL2ND_WFB(J)  * TauU_RA(J-1,M,2,N) ) &
                            )*DYFI(J)+  ( &
                            ( YCL2ND_WFF(J+1)* TauU_RA(J,  N,2,M)   &
                             +YCL2ND_WFB(J+1)* TauU_RA(J+1,N,2,M) )-&
                            ( YCL2ND_WFF(J)  * TauU_RA(J,  N,2,M)   &
                             +YCL2ND_WFB(J)  * TauU_RA(J-1,N,2,M) ) &
                            )*DYFI(J)
                    END IF

                    

                    
                    BUDG_vistress3_duiuj(J,L) = FF &
                            - TauDvDL_RA(J,M,1,N,1) - TauDvDL_RA(J,M,2,N,2) - TauDvDL_RA(J,M,3,N,3) &
                            - TauDvDL_RA(J,N,1,M,1) - TauDvDL_RA(J,N,2,M,2) - TauDvDL_RA(J,N,3,M,3) &
                            - U_FA(J,M)*dTaudy_RA(J,2,N) - U_FA(J,N)*dTaudy_RA(J,2,M)
                            
                    !IF (L==1) THEN
                        !!WRITE(*,*) FF, - TauDvDL_RA(J,M,1,N,1) - TauDvDL_RA(J,M,2,N,2) - TauDvDL_RA(J,M,3,N,3), &
                        !- TauDvDL_RA(J,N,1,M,1) - TauDvDL_RA(J,N,2,M,2) - TauDvDL_RA(J,N,3,M,3), &
                        !- U_FA(J,M)*dTaudy_RA(J,2,N), - U_FA(J,N)*dTaudy_RA(J,2,M)
                        !!WRITE(*,*) - U_FA(J,M)*dTaudy_RA(J,2,N), - U_FA(J,N)*dTaudy_RA(J,2,M), U_FA(J,M), dTaudy_RA(J,2,M)
                    !END IF 
                    
                END DO
            END DO
        
        END DO
        
!====FA===========BALANCE===================
!==============buoyancy force direct production is zero======================
        DO J=1, NCL2
            !==============for TKE and MKE==================
            BUDG_balance1_TKE(J) =    BUDG_prodc_stres_TKE(J) + &
                                      BUDG_viscs_dissp_TKE(J) + &
                                      BUDG_pdudx_stran_TKE(J) + &
                                      BUDG_Turbu_diffu_TKE(J) + &
                                      BUDG_dpudx_diffu_TKE(J) + &
                                      BUDG_viscs_diffu_TKE(J) + &
                                      BUDG_press_accl1_TKE(J) + &
                                      BUDG_viscs_accl1_TKE(J) + &
                                      BUDG_prodc_dvfc1_TKE(J)
                                      
            BUDG_balance2_TKE(J) =    BUDG_prodc_stres_TKE(J) + &
                                      BUDG_viscs_dissp_TKE(J) + &
                                      BUDG_pdudx_stran_TKE(J) + &
                                      BUDG_Turbu_diffu_TKE(J) + &
                                      BUDG_dpudx_diffu_TKE(J) + &
                                      BUDG_viscs_diffu_TKE(J) + &
                                      BUDG_turss_accl2_TKE(J) + &
                                      BUDG_prodc_dvfc2_TKE(J) + &
                                      BUDG_prodc_gvfc2_TKE(J)
                                      
            BUDG_balance3_TKE(J) =    BUDG_prodc_stres_TKE(J) + &
                                      BUDG_Turbu_diffu_TKE(J) + &
                                      BUDG_prodc_dvfc1_TKE(J) + &
                                      BUDG_pressure3_TKE(J)   + &
                                      BUDG_vistress3_TKE(J)    
                                      
            BUDG_balance1_MKE(J) =    BUDG_prodc_stres_MKE(J) + &
                                      BUDG_viscs_dissp_MKE(J) + &
                                      BUDG_pdudx_stran_MKE(J) + &
                                      BUDG_Turbu_diffu_MKE(J) + &
                                      BUDG_dpudx_diffu_MKE(J) + &
                                      BUDG_viscs_diffu_MKE(J) + &
                                      BUDG_press_accl1_MKE(J) + &
                                      BUDG_viscs_accl1_MKE(J) + &
                                      BUDG_prodc_gvfc1_MKE(J) + &
                                      BUDG_prodc_dvfc1_MKE(J)
                                      
            BUDG_balance2_MKE(J) =    BUDG_prodc_stres_MKE(J) + &
                                      BUDG_viscs_dissp_MKE(J) + &
                                      BUDG_pdudx_stran_MKE(J) + &
                                      BUDG_Turbu_diffu_MKE(J) + &
                                      BUDG_dpudx_diffu_MKE(J) + &
                                      BUDG_viscs_diffu_MKE(J) + &
                                      BUDG_turss_accl2_MKE(J) + &
                                      BUDG_prodc_gvfc2_MKE(J) + &
                                      BUDG_prodc_dvfc2_MKE(J)
                                      
            BUDG_balance3_MKE(J) =    BUDG_prodc_stres_MKE(J) + &
                                      BUDG_Turbu_diffu_MKE(J) + &
                                      BUDG_prodc_dvfc1_MKE(J) + &
                                      BUDG_pressure3_MKE(J)    + &
                                      BUDG_vistress3_MKE(J)    
                                      
            !==========for each Ruv=========================                                  
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    BUDG_balance1_duiuj(J,L) =BUDG_prodc_stres_duiuj(J,L) + &
                                              BUDG_viscs_dissp_duiuj(J,L) + &
                                              BUDG_pdudx_stran_duiuj(J,L) + &
                                              BUDG_Turbu_diffu_duiuj(J,L) + &
                                              BUDG_dpudx_diffu_duiuj(J,L) + &
                                              BUDG_viscs_diffu_duiuj(J,L) + &
                                              BUDG_press_accl1_duiuj(J,L) + &
                                              BUDG_viscs_accl1_duiuj(J,L) + &
                                              BUDG_prodc_dvfc1_duiuj(J,L)
                                              
                    BUDG_balance2_duiuj(J,L) =BUDG_prodc_stres_duiuj(J,L) + &
                                              BUDG_viscs_dissp_duiuj(J,L) + &
                                              BUDG_pdudx_stran_duiuj(J,L) + &
                                              BUDG_Turbu_diffu_duiuj(J,L) + &
                                              BUDG_dpudx_diffu_duiuj(J,L) + &
                                              BUDG_viscs_diffu_duiuj(J,L) + &
                                              BUDG_turss_accl2_duiuj(J,L) + &
                                              BUDG_prodc_gvfc2_duiuj(J,L) + &
                                              BUDG_prodc_dvfc2_duiuj(J,L)
                                              
                    BUDG_balance3_duiuj(J,L) =BUDG_prodc_stres_duiuj(J,L) + &
                                              BUDG_Turbu_diffu_duiuj(J,L) + &
                                              BUDG_prodc_dvfc1_duiuj(J,L) + &
                                              BUDG_pressure3_duiuj(J,L) + &
                                              BUDG_vistress3_duiuj(J,L)
                                              
                    
                END DO
            END DO
        END DO
        
        
!================integral of each terms along Y==================================
        DO J=1, NCL2
            !============TKE======================================================================
            BUDG_prodc_stres_TKE_ysum = BUDG_prodc_stres_TKE_ysum + BUDG_prodc_stres_TKE(J)/DYFI(J)
            BUDG_viscs_dissp_TKE_ysum = BUDG_viscs_dissp_TKE_ysum + BUDG_viscs_dissp_TKE(J)/DYFI(J)
            BUDG_pdudx_stran_TKE_ysum = BUDG_pdudx_stran_TKE_ysum + BUDG_pdudx_stran_TKE(J)/DYFI(J)
            BUDG_Turbu_diffu_TKE_ysum = BUDG_Turbu_diffu_TKE_ysum + BUDG_Turbu_diffu_TKE(J)/DYFI(J)
            BUDG_dpudx_diffu_TKE_ysum = BUDG_dpudx_diffu_TKE_ysum + BUDG_dpudx_diffu_TKE(J)/DYFI(J)
            BUDG_viscs_diffu_TKE_ysum = BUDG_viscs_diffu_TKE_ysum + BUDG_viscs_diffu_TKE(J)/DYFI(J)
            
            BUDG_press_accl1_TKE_ysum = BUDG_press_accl1_TKE_ysum + BUDG_press_accl1_TKE(J)/DYFI(J)
            BUDG_viscs_accl1_TKE_ysum = BUDG_viscs_accl1_TKE_ysum + BUDG_viscs_accl1_TKE(J)/DYFI(J)
            BUDG_prodc_dvfc1_TKE_ysum = BUDG_prodc_dvfc1_TKE_ysum + BUDG_prodc_dvfc1_TKE(J)/DYFI(J)
            BUDG_balance1_TKE_ysum    = BUDG_balance1_TKE_ysum    + BUDG_balance1_TKE(J)/DYFI(J)
             
            BUDG_turss_accl2_TKE_ysum = BUDG_turss_accl2_TKE_ysum + BUDG_turss_accl2_TKE(J)/DYFI(J)
            BUDG_prodc_gvfc2_TKE_ysum = BUDG_prodc_gvfc2_TKE_ysum + BUDG_prodc_gvfc2_TKE(J)/DYFI(J)
            BUDG_prodc_dvfc2_TKE_ysum = BUDG_prodc_dvfc2_TKE_ysum + BUDG_prodc_dvfc2_TKE(J)/DYFI(J)
            BUDG_balance2_TKE_ysum    = BUDG_balance2_TKE_ysum    + BUDG_balance2_TKE(J)/DYFI(J)
        
            
            BUDG_pressure3_TKE_ysum   = BUDG_pressure3_TKE_ysum   + BUDG_pressure3_TKE(J)/DYFI(J)
            BUDG_vistress3_TKE_ysum   = BUDG_vistress3_TKE_ysum   + BUDG_vistress3_TKE(J)/DYFI(J)
            BUDG_balance3_TKE_ysum    = BUDG_balance3_TKE_ysum    + BUDG_balance3_TKE(J)/DYFI(J)
            
            
            !===========MKE=======================================================================
            BUDG_prodc_stres_MKE_ysum = BUDG_prodc_stres_MKE_ysum + BUDG_prodc_stres_MKE(J)/DYFI(J)
            BUDG_viscs_dissp_MKE_ysum = BUDG_viscs_dissp_MKE_ysum + BUDG_viscs_dissp_MKE(J)/DYFI(J)
            BUDG_pdudx_stran_MKE_ysum = BUDG_pdudx_stran_MKE_ysum + BUDG_pdudx_stran_MKE(J)/DYFI(J)
            BUDG_Turbu_diffu_MKE_ysum = BUDG_Turbu_diffu_MKE_ysum + BUDG_Turbu_diffu_MKE(J)/DYFI(J)
            BUDG_dpudx_diffu_MKE_ysum = BUDG_dpudx_diffu_MKE_ysum + BUDG_dpudx_diffu_MKE(J)/DYFI(J)
            BUDG_viscs_diffu_MKE_ysum = BUDG_viscs_diffu_MKE_ysum + BUDG_viscs_diffu_MKE(J)/DYFI(J)
            
            BUDG_press_accl1_MKE_ysum = BUDG_press_accl1_MKE_ysum + BUDG_press_accl1_MKE(J)/DYFI(J)
            BUDG_viscs_accl1_MKE_ysum = BUDG_viscs_accl1_MKE_ysum + BUDG_viscs_accl1_MKE(J)/DYFI(J)
            BUDG_prodc_dvfc1_MKE_ysum = BUDG_prodc_dvfc1_MKE_ysum + BUDG_prodc_dvfc1_MKE(J)/DYFI(J)
            BUDG_prodc_gvfc1_MKE_ysum = BUDG_prodc_gvfc1_MKE_ysum + BUDG_prodc_gvfc1_MKE(J)/DYFI(J)
            BUDG_balance1_MKE_ysum    = BUDG_balance1_MKE_ysum    + BUDG_balance1_MKE(J)/DYFI(J)
             
            BUDG_turss_accl2_MKE_ysum = BUDG_turss_accl2_MKE_ysum + BUDG_turss_accl2_MKE(J)/DYFI(J)
            BUDG_prodc_gvfc2_MKE_ysum = BUDG_prodc_gvfc2_MKE_ysum + BUDG_prodc_gvfc2_MKE(J)/DYFI(J)
            BUDG_prodc_dvfc2_MKE_ysum = BUDG_prodc_dvfc2_MKE_ysum + BUDG_prodc_dvfc2_MKE(J)/DYFI(J)
            BUDG_balance2_MKE_ysum    = BUDG_balance2_MKE_ysum    + BUDG_balance2_MKE(J)/DYFI(J)
        
            
            BUDG_pressure3_MKE_ysum   = BUDG_pressure3_MKE_ysum   + BUDG_pressure3_MKE(J)/DYFI(J)
            BUDG_vistress3_MKE_ysum   = BUDG_vistress3_MKE_ysum   + BUDG_vistress3_MKE(J)/DYFI(J)
            BUDG_balance3_MKE_ysum    = BUDG_balance3_MKE_ysum    + BUDG_balance3_MKE(J)/DYFI(J)
            
            
            DO L=1, (NDV*(7-NDV))/2+NDV-3
                BUDG_prodc_stres_duiuj_ysum(L) = BUDG_prodc_stres_duiuj_ysum(L) + BUDG_prodc_stres_duiuj(J,L)/DYFI(J)
                BUDG_viscs_dissp_duiuj_ysum(L) = BUDG_viscs_dissp_duiuj_ysum(L) + BUDG_viscs_dissp_duiuj(J,L)/DYFI(J)
                BUDG_pdudx_stran_duiuj_ysum(L) = BUDG_pdudx_stran_duiuj_ysum(L) + BUDG_pdudx_stran_duiuj(J,L)/DYFI(J)
                BUDG_Turbu_diffu_duiuj_ysum(L) = BUDG_Turbu_diffu_duiuj_ysum(L) + BUDG_Turbu_diffu_duiuj(J,L)/DYFI(J)
                BUDG_dpudx_diffu_duiuj_ysum(L) = BUDG_dpudx_diffu_duiuj_ysum(L) + BUDG_dpudx_diffu_duiuj(J,L)/DYFI(J)
                BUDG_viscs_diffu_duiuj_ysum(L) = BUDG_viscs_diffu_duiuj_ysum(L) + BUDG_viscs_diffu_duiuj(J,L)/DYFI(J)
                
                BUDG_press_accl1_duiuj_ysum(L) = BUDG_press_accl1_duiuj_ysum(L) + BUDG_press_accl1_duiuj(J,L)/DYFI(J)
                BUDG_viscs_accl1_duiuj_ysum(L) = BUDG_viscs_accl1_duiuj_ysum(L) + BUDG_viscs_accl1_duiuj(J,L)/DYFI(J)
                BUDG_prodc_dvfc1_duiuj_ysum(L) = BUDG_prodc_dvfc1_duiuj_ysum(L) + BUDG_prodc_dvfc1_duiuj(J,L)/DYFI(J)
                BUDG_balance1_duiuj_ysum(L)    = BUDG_balance1_duiuj_ysum(L)    + BUDG_balance1_duiuj(J,L)/DYFI(J)
                 
                BUDG_turss_accl2_duiuj_ysum(L) = BUDG_turss_accl2_duiuj_ysum(L) + BUDG_turss_accl2_duiuj(J,L)/DYFI(J)
                BUDG_prodc_gvfc2_duiuj_ysum(L) = BUDG_prodc_gvfc2_duiuj_ysum(L) + BUDG_prodc_gvfc2_duiuj(J,L)/DYFI(J)
                BUDG_prodc_dvfc2_duiuj_ysum(L) = BUDG_prodc_dvfc2_duiuj_ysum(L) + BUDG_prodc_dvfc2_duiuj(J,L)/DYFI(J)
                BUDG_balance2_duiuj_ysum(L)    = BUDG_balance2_duiuj_ysum(L)    + BUDG_balance2_duiuj(J,L)/DYFI(J)
            
                
                BUDG_pressure3_duiuj_ysum(L)   = BUDG_pressure3_duiuj_ysum(L)   + BUDG_pressure3_duiuj(J,L)/DYFI(J)
                BUDG_vistress3_duiuj_ysum(L)   = BUDG_vistress3_duiuj_ysum(L)   + BUDG_vistress3_duiuj(J,L)/DYFI(J)
                BUDG_balance3_duiuj_ysum(L)    = BUDG_balance3_duiuj_ysum(L)    + BUDG_balance3_duiuj(J,L)/DYFI(J)
            END DO
            
        END DO
        
!        !=============below for validation========(theoretically zero)===================
!        sum_p_related_mke=0.0_WP
!        DO J=1, NCL2
!            sum_p_related_mke = BUDG_pdudx_stran_MKE(J)+BUDG_dpudx_diffu_MKE(J)+BUDG_press_accl1_MKE(J)
!            !WRITE(*,*) '##Sum of P-MKE:', J, sum_p_related_mke
!        END DO
!        !================
        
        !========sum of the pressure term related in MKE is zero=========================
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        FLNM=TRIM(filepath4)//'Result.IO.budget.check.'//TRIM(PNTIM)//'.tec'
        INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
        if(file_exists) then
            OPEN(TECFLG,FILE=FLNM, POSITION='APPEND')
        else
            OPEN(TECFLG,FILE=FLNM)
        end if
        
        WRITE(TECFLG,'(A)') '=======FA=====uu,uv,uw,vv,vw,ww=========='
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_prodc_stres_duiuj_ysum(1:6) = ', BUDG_prodc_stres_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_viscs_dissp_duiuj_ysum(1:6) = ', BUDG_viscs_dissp_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_pdudx_stran_duiuj_ysum(1:6) = ', BUDG_pdudx_stran_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_Turbu_diffu_duiuj_ysum(1:6) = ', BUDG_Turbu_diffu_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_dpudx_diffu_duiuj_ysum(1:6) = ', BUDG_dpudx_diffu_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_viscs_diffu_duiuj_ysum(1:6) = ', BUDG_viscs_diffu_duiuj_ysum(1:6)
        
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_press_accl1_duiuj_ysum(1:6) = ', BUDG_press_accl1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_viscs_accl1_duiuj_ysum(1:6) = ', BUDG_viscs_accl1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_prodc_dvfc1_duiuj_ysum(1:6) = ', BUDG_prodc_dvfc1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_balance1_duiuj_ysum(1:6) =     ', BUDG_balance1_duiuj_ysum(1:6)
        
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_turss_accl2_duiuj_ysum(1:6) = ', BUDG_turss_accl2_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_prodc_gvfc2_duiuj_ysum(1:6) = ', BUDG_prodc_gvfc2_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_prodc_dvfc2_duiuj_ysum(1:6) = ', BUDG_prodc_dvfc2_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_balance2_duiuj_ysum(1:6)    = ', BUDG_balance2_duiuj_ysum(1:6)
        
        
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_pressure3_duiuj_ysum(1:6)(calc) = ',  BUDG_pressure3_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_pressure3_duiuj_ysum(1:6)(addt) = ',  BUDG_pdudx_stran_duiuj_ysum(1:6) + &
                                                                                BUDG_dpudx_diffu_duiuj_ysum(1:6) + &
                                                                                BUDG_press_accl1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_vistress3_duiuj_ysum(1:6)(calc) = ',  BUDG_vistress3_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_vistress3_duiuj_ysum(1:6)(addt) = ',  BUDG_viscs_dissp_duiuj_ysum(1:6) + &
                                                                                BUDG_viscs_diffu_duiuj_ysum(1:6) + &
                                                                                BUDG_viscs_accl1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_balance3_duiuj_ysum(1:6)        = ',  BUDG_balance3_duiuj_ysum(1:6)
        

        WRITE(TECFLG,'(A)') '=======FA=====TKE, MKE=========='
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_prodc_stres_ysum_TKE_and_MKE = ', BUDG_prodc_stres_TKE_ysum, BUDG_prodc_stres_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_viscs_dissp_ysum_TKE_and_MKE = ', BUDG_viscs_dissp_TKE_ysum, BUDG_viscs_dissp_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_pdudx_stran_ysum_TKE_and_MKE = ', BUDG_pdudx_stran_TKE_ysum, BUDG_pdudx_stran_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_Turbu_diffu_ysum_TKE_and_MKE = ', BUDG_Turbu_diffu_TKE_ysum, BUDG_Turbu_diffu_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_dpudx_diffu_ysum_TKE_and_MKE = ', BUDG_dpudx_diffu_TKE_ysum, BUDG_dpudx_diffu_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_viscs_diffu_ysum_TKE_and_MKE = ', BUDG_viscs_diffu_TKE_ysum, BUDG_viscs_diffu_MKE_ysum
        
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_press_accl1_ysum_TKE_and_MKE = ', BUDG_press_accl1_TKE_ysum, BUDG_press_accl1_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_viscs_accl1_ysum_TKE_and_MKE = ', BUDG_viscs_accl1_TKE_ysum, BUDG_viscs_accl1_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_prodc_dvfc1_ysum_TKE_and_MKE = ', BUDG_prodc_dvfc1_TKE_ysum, BUDG_prodc_dvfc1_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_balance1_ysum_TKE_and_MKE    = ', BUDG_balance1_TKE_ysum,    BUDG_balance1_MKE_ysum
        
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_turss_accl2_ysum_TKE_and_MKE = ', BUDG_turss_accl2_TKE_ysum, BUDG_turss_accl2_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_prodc_gvfc2_ysum_TKE_and_MKE = ', BUDG_prodc_gvfc2_TKE_ysum, BUDG_prodc_gvfc2_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_prodc_dvfc2_ysum_TKE_and_MKE = ', BUDG_prodc_dvfc2_TKE_ysum, BUDG_prodc_dvfc2_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_balance2_ysum_TKE_and_MKE    = ', BUDG_balance2_TKE_ysum,    BUDG_balance2_MKE_ysum
        
        
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_pressure3_ysum_TKE_and_MKE(calc) = ', BUDG_pressure3_TKE_ysum,BUDG_pressure3_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_pressure3_ysum_TKE_and_MKE(addt) = ', BUDG_pdudx_stran_TKE_ysum + &
                                                                                BUDG_dpudx_diffu_TKE_ysum + &
                                                                                BUDG_press_accl1_TKE_ysum,  &
                                                                                BUDG_pdudx_stran_MKE_ysum + &
                                                                                BUDG_dpudx_diffu_MKE_ysum + &
                                                                                BUDG_press_accl1_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_vistress3_ysum_TKE_and_MKE(calc) = ', BUDG_vistress3_TKE_ysum, BUDG_vistress3_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_vistress3_ysum_TKE_and_MKE(addt) = ', BUDG_viscs_dissp_TKE_ysum + &
                                                                                BUDG_viscs_diffu_TKE_ysum + &
                                                                                BUDG_viscs_accl1_TKE_ysum, &
                                                                                BUDG_viscs_dissp_MKE_ysum + &
                                                                                BUDG_viscs_diffu_MKE_ysum + &
                                                                                BUDG_viscs_accl1_MKE_ysum
        WRITE(TECFLG,'(A,2ES20.7)') 'BUDG_balance3_ysum_TKE_and_MKE        = ', BUDG_balance3_TKE_ysum, BUDG_balance3_MKE_ysum
                                                                    
        CLOSE(TECFLG)                                    
        
        
        
    RETURN
    END SUBROUTINE 

!###########################  RA   ##############################################################    
    SUBROUTINE PP_FLOW_RA_noDen_RSTE_BUDG_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: J
        INTEGER(4) :: M, N, H, L, K
        REAL(WP)   :: BUDG_prod_Pij1, BUDG_prod_Pij2, COE
        INTEGER(4) :: TECFLG=200
        character(128) :: FLNM
        logical        :: file_exists
        
        
        CAll CHKHDL('=====Calculating RA Budegts=====',myid)
        CALL PP_BUDG_INIT
    !=======PRODUCTION TERMS due to mean shear===(RS)=========================================================
!       Eq: BUDG_prodc_stres_duiuj(J,L) = P*_{ij}=P_{ij}-2/3P\delta_{ij}
!       Eq: P_{ij} = -<\rho u''_i u''_k> (\partial {u_j} / \partial x_k) + 
!                    -<\rho u''_j u''_k> (\partial {u_i} / \partial x_k)
        BUDG_prodc_stres_duiuj = 0.0_wp
        DO J=1,NCL2
            !==============for TKE and MKE==================
            BUDG_prodc_stres_MKE(J)=uf2_RA(J,1,1)*DVDL1xztL_F0_io(J,1,1) + &
                                    uf2_RA(J,1,2)*DVDL1xztL_F0_io(J,1,2) + &
                                    uf2_RA(J,1,3)*DVDL1xztL_F0_io(J,1,3) + &
                                    uf2_RA(J,2,1)*DVDL1xztL_F0_io(J,2,1) + &
                                    uf2_RA(J,2,2)*DVDL1xztL_F0_io(J,2,2) + &
                                    uf2_RA(J,2,3)*DVDL1xztL_F0_io(J,2,3) + &
                                    uf2_RA(J,3,1)*DVDL1xztL_F0_io(J,3,1) + &
                                    uf2_RA(J,3,2)*DVDL1xztL_F0_io(J,3,2) + &
                                    uf2_RA(J,3,3)*DVDL1xztL_F0_io(J,3,3) 
            BUDG_prodc_stres_TKE(J)=-1.0_wp*BUDG_prodc_stres_MKE(J)
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    BUDG_prod_Pij1= uf2_RA(J,1,M)*DVDL1xztL_F0_io(J,N,1) + &
                                    uf2_RA(J,2,M)*DVDL1xztL_F0_io(J,N,2) + &
                                    uf2_RA(J,3,M)*DVDL1xztL_F0_io(J,N,3) 
                    BUDG_prod_Pij2= uf2_RA(J,1,N)*DVDL1xztL_F0_io(J,M,1) + &
                                    uf2_RA(J,2,N)*DVDL1xztL_F0_io(J,M,2) + &
                                    uf2_RA(J,3,N)*DVDL1xztL_F0_io(J,M,3) 
                    BUDG_prodc_stres_duiuj(J,L) = -1.0_wp*(BUDG_prod_Pij1 + BUDG_prod_Pij2) 
                              
                END DO
            END DO
            
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_prodc_stres_duiuj',myid)
        
!==RA=====Viscous ENERGY dissipation term ===========================================================
!       ViscStressVeloGrad_RA(J,M,N,H,P) = <\partial(u_m)/\partial x_n \tau_hp>
        BUDG_viscs_dissp_TKE   = 0.0_WP
        BUDG_viscs_dissp_MKE   = 0.0_WP
        BUDG_viscs_dissp_duiuj = 0.0_WP
        DO J=1, NCL2    
            !==============for TKE and MKE==================
            BUDG_viscs_dissp_TKE(J) = &
                 DVDL2xztL_F0_io(J,(1-1)*3+1,(1-1)*3+1) - DVDL1xztL_F0_io(J,1,1)*DVDL1xztL_F0_io(J,1,1) & 
                +DVDL2xztL_F0_io(J,(1-1)*3+2,(1-1)*3+2) - DVDL1xztL_F0_io(J,1,2)*DVDL1xztL_F0_io(J,1,2) &  
                +DVDL2xztL_F0_io(J,(1-1)*3+3,(1-1)*3+3) - DVDL1xztL_F0_io(J,1,3)*DVDL1xztL_F0_io(J,1,3) &
                +DVDL2xztL_F0_io(J,(2-1)*3+1,(2-1)*3+1) - DVDL1xztL_F0_io(J,2,1)*DVDL1xztL_F0_io(J,2,1) & 
                +DVDL2xztL_F0_io(J,(2-1)*3+2,(2-1)*3+2) - DVDL1xztL_F0_io(J,2,2)*DVDL1xztL_F0_io(J,2,2) &  
                +DVDL2xztL_F0_io(J,(2-1)*3+3,(2-1)*3+3) - DVDL1xztL_F0_io(J,2,3)*DVDL1xztL_F0_io(J,2,3) &
                +DVDL2xztL_F0_io(J,(3-1)*3+1,(3-1)*3+1) - DVDL1xztL_F0_io(J,3,1)*DVDL1xztL_F0_io(J,3,1) & 
                +DVDL2xztL_F0_io(J,(3-1)*3+2,(3-1)*3+2) - DVDL1xztL_F0_io(J,3,2)*DVDL1xztL_F0_io(J,3,2) &  
                +DVDL2xztL_F0_io(J,(3-1)*3+3,(3-1)*3+3) - DVDL1xztL_F0_io(J,3,3)*DVDL1xztL_F0_io(J,3,3) 
            BUDG_viscs_dissp_TKE(J) = BUDG_viscs_dissp_TKE(J)*CVISC*(-1.0_wp)
            
            BUDG_viscs_dissp_MKE(J) = &  
                Tau_Mean_RA(J,1,1) * DVDL1xztL_F0_io(J,1,1) + &
                Tau_Mean_RA(J,2,1) * DVDL1xztL_F0_io(J,2,1) + &
                Tau_Mean_RA(J,3,1) * DVDL1xztL_F0_io(J,3,1) + &
                Tau_Mean_RA(J,1,2) * DVDL1xztL_F0_io(J,1,2) + &
                Tau_Mean_RA(J,2,2) * DVDL1xztL_F0_io(J,2,2) + &
                Tau_Mean_RA(J,3,2) * DVDL1xztL_F0_io(J,3,2) + &
                Tau_Mean_RA(J,1,3) * DVDL1xztL_F0_io(J,1,3) + &
                Tau_Mean_RA(J,2,3) * DVDL1xztL_F0_io(J,2,3) + &
                Tau_Mean_RA(J,3,3) * DVDL1xztL_F0_io(J,3,3)
            BUDG_viscs_dissp_MKE(J) = -1.0_wp*BUDG_viscs_dissp_MKE(J) 
            
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    BUDG_viscs_dissp_duiuj(J,L)=  &
                         DVDL2xztL_F0_io(J,(M-1)*3+1,(N-1)*3+1) - DVDL1xztL_F0_io(J,M,1)*DVDL1xztL_F0_io(J,N,1)             & 
                        +DVDL2xztL_F0_io(J,(M-1)*3+2,(N-1)*3+2) - DVDL1xztL_F0_io(J,M,2)*DVDL1xztL_F0_io(J,N,2)             &  
                        +DVDL2xztL_F0_io(J,(M-1)*3+3,(N-1)*3+3) - DVDL1xztL_F0_io(J,M,3)*DVDL1xztL_F0_io(J,N,3)             
                    
                    BUDG_viscs_dissp_duiuj(J,L)=BUDG_viscs_dissp_duiuj(J,L)*CVISC*(-2.0_wp)
                    !IF(dabs(dUdX_RA(J,M,N)-DVDL1xztL_F0_io(J,M,N)).GT.1.0E-12_wp) &
                    !WRITE(*,'(3I4.1,3ES13.5)') J,M,N, dUdX_RA(J,M,N),DVDL1xztL_F0_io(J,M,N),dUdX_RA(J,M,N)-DVDL1xztL_F0_io(J,M,N)
                    
                END DO
            END DO
        END DO
        
        CALL CHKHDL('      ==>Calculated BUDG_viscs_dissp_duiuj',myid) 
        

!==RA=====Velocity-Pressure gradient Pressure strain term=======(RS)============================== 
        DO J=1,NCL2
            !==============for TKE and MKE==================
            BUDG_pdudx_stran_TKE(J) =   DVDLPxztL_F0_io(J,1,1) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,1,1) + &
                                        DVDLPxztL_F0_io(J,2,2) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,2,2) + &
                                        DVDLPxztL_F0_io(J,3,3) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,3,3)
            BUDG_pdudx_stran_MKE(J) =   U1xztL_F0_io(J,4)* dUidXi(J)
            
            DO M=1,NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    BUDG_pdudx_stran_duiuj(J,L) = DVDLPxztL_F0_io(J,M,N) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,M,N) + &
                                                  DVDLPxztL_F0_io(J,N,M) - U1xztL_F0_io(J,4)*DVDL1xztL_F0_io(J,N,M)
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_pdudx_stran_duiuj',myid)
        
                
!==RA=====TURBUELCEN DIFFUSION TERMS=(turblence tranport rate)===(RS)========================
!       Eq.  BUDG_Tdiff_duiuj(J,L) = ( \partial <\rho u''_i u''_j u''_k > ) / (\partial x_k )
        DO J=1,NCL2   
            ufMKEfd_RA(J)=  0.5_WP*( U3xztL_F0_io(J,2) + U3xztL_F0_io(J,7) + U3xztL_F0_io(J,9) )
                                
            ufTKEfd_RA(J) = 0.5_WP*( uf3_RA(J,1,1,2) + uf3_RA(J,2,2,2) + uf3_RA(J,3,3,2) )
            
        END DO
        
        
        DO J=1,NCL2 
            !==============for TKE and MKE==================
            IF(J==1) THEN
                BUDG_Turbu_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufTKEfd_RA(J+1) + &
                        YCL2ND_WFB(J+1)*ufTKEfd_RA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_wp)
                BUDG_Turbu_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufMKEfd_RA(J+1) + &
                        YCL2ND_WFB(J+1)*ufMKEfd_RA(J  ) ) - 0.0_WP ) * DYFI(J) * (-1.0_wp)
            ELSE IF (J==NCL2) THEN
                BUDG_Turbu_diffu_TKE(J)= &
                    ( 0.0_WP - &
                      ( YCL2ND_WFF(J)*ufTKEfd_RA(J  ) + &
                        YCL2ND_WFB(J)*ufTKEfd_RA(J-1) ) ) * DYFI(J) * (-1.0_wp)
                BUDG_Turbu_diffu_MKE(J)= &
                    ( 0.0_WP - &
                      ( YCL2ND_WFF(J)*ufMKEfd_RA(J  ) + &
                        YCL2ND_WFB(J)*ufMKEfd_RA(J-1) ) ) * DYFI(J) * (-1.0_wp)
                                  
            ELSE
                BUDG_Turbu_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufTKEfd_RA(J+1) + &
                        YCL2ND_WFB(J+1)*ufTKEfd_RA(J  ) ) -         &
                      ( YCL2ND_WFF(J)  *ufTKEfd_RA(J  ) + &
                        YCL2ND_WFB(J)  *ufTKEfd_RA(J-1) ) ) * DYFI(J) * (-1.0_wp)
                BUDG_Turbu_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufMKEfd_RA(J+1) + &
                        YCL2ND_WFB(J+1)*ufMKEfd_RA(J  ) ) -         &
                      ( YCL2ND_WFF(J)  *ufMKEfd_RA(J  ) + &
                        YCL2ND_WFB(J)  *ufMKEfd_RA(J-1) ) ) * DYFI(J) * (-1.0_wp)
            END IF
            
 
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    IF(J==1) THEN
                        BUDG_Turbu_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*uf3_RA(J+1,M,N,2) + &
                                YCL2ND_WFB(J+1)*uf3_RA(J,  M,N,2) ) - 0.0_WP ) * DYFI(J) * (-1.0_wp)
                    ELSE IF (J==NCL2) THEN
                        BUDG_Turbu_diffu_duiuj(J,L)= &
                            ( 0.0_WP - &
                              ( YCL2ND_WFF(J)*uf3_RA(J,  M,N,2) + &
                                YCL2ND_WFB(J)*uf3_RA(J-1,M,N,2) ) ) * DYFI(J) * (-1.0_wp)
                                          
                    ELSE
                        BUDG_Turbu_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*uf3_RA(J+1,M,N,2) + &
                                YCL2ND_WFB(J+1)*uf3_RA(J,  M,N,2) ) -         &
                              ( YCL2ND_WFF(J)  *uf3_RA(J,  M,N,2) + &
                                YCL2ND_WFB(J)  *uf3_RA(J-1,M,N,2) ) ) * DYFI(J) * (-1.0_wp)
                    END IF
                    
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_Turbu_diffu_duiuj',myid)
        
!==RA=====Velocity-Pressure gradient diffusion term=======(RS)==============================        
!       Eq. = - \partial (<p' u''_j>) /\partial (x_i) - \partial (<p' u''_i>) /\partial (x_j)
        DO J=1,NCL2
            !==============for TKE and MKE==================
            IF(J==1) THEN
                BUDG_dpudx_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufpf_RA(J+1,2) + &
                        YCL2ND_WFB(J+1)*ufpf_RA(J,  2) ) &
                         -  0.0_WP )*DYFI(J)* (-1.0_wp)
                         
                BUDG_dpudx_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*(U1xztL_F0_io(J+1,2)*U1xztL_F0_io(J+1,4)) + &
                        YCL2ND_WFB(J+1)*(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) ) &
                         -  0.0_WP )*DYFI(J)* (-1.0_wp)
                         
            ELSE IF(J==NCL2) THEN
                BUDG_dpudx_diffu_TKE(J)= &
                    ( 0.0_WP -  &
                    ( YCL2ND_WFF(J)  *ufpf_RA(J,  2) + &
                      YCL2ND_WFB(J)  *ufpf_RA(J-1,2) ) )*DYFI(J)* (-1.0_wp)
                      
                BUDG_dpudx_diffu_MKE(J)= &
                    ( 0.0_WP -  &
                    ( YCL2ND_WFF(J)  *(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) + &
                      YCL2ND_WFB(J)  *(U1xztL_F0_io(J-1,2)*U1xztL_F0_io(J-1,4)) ))*DYFI(J)* (-1.0_wp)
            ELSE
            
                BUDG_dpudx_diffu_TKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*ufpf_RA(J+1,2) + &
                        YCL2ND_WFB(J+1)*ufpf_RA(J,  2) ) -  &
                      ( YCL2ND_WFF(J)  *ufpf_RA(J,  2) + &
                        YCL2ND_WFB(J)  *ufpf_RA(J-1,2) ) )*DYFI(J)* (-1.0_wp)
                                    
                BUDG_dpudx_diffu_MKE(J)= &
                    ( ( YCL2ND_WFF(J+1)*(U1xztL_F0_io(J+1,2)*U1xztL_F0_io(J+1,4)) + &
                        YCL2ND_WFB(J+1)*(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) ) -  &
                      ( YCL2ND_WFF(J)  *(U1xztL_F0_io(J,  2)*U1xztL_F0_io(J,  4)) + &
                        YCL2ND_WFB(J)  *(U1xztL_F0_io(J-1,2)*U1xztL_F0_io(J-1,4)) ) )*DYFI(J)* (-1.0_wp)
            END IF
            
            !==========for each Ruv=========================
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    IF(J==1) THEN
                        BUDG_dpudx_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*( ufpf_RA(J+1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J+1,N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J+1)*( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   ) )- &
                                0.0_WP )*DYFI(J)
                    ELSE IF(J==NCL2) THEN
                        BUDG_dpudx_diffu_duiuj(J,L)= &
                            (  0.0_WP- &
                              ( YCL2ND_WFF(J)  *( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J)  *( ufpf_RA(J-1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J-1,N)*DBLE(Kronecker_delta(M,2))   ) ) )*DYFI(J)
                    ELSE
                        BUDG_dpudx_diffu_duiuj(J,L)= &
                            ( ( YCL2ND_WFF(J+1)*( ufpf_RA(J+1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J+1,N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J+1)*( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   ) )- &
                              ( YCL2ND_WFF(J)  *( ufpf_RA(J,  M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J,  N)*DBLE(Kronecker_delta(M,2))   )+   &
                                YCL2ND_WFB(J)  *( ufpf_RA(J-1,M)*DBLE(Kronecker_delta(N,2)) +      &
                                                  ufpf_RA(J-1,N)*DBLE(Kronecker_delta(M,2))   ) ) )*DYFI(J)
                    END IF
                    
                    
                    BUDG_dpudx_diffu_duiuj(J,L) = BUDG_dpudx_diffu_duiuj(J,L)* (-1.0_wp)
                END DO
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_dpudx_diffu_duiuj',myid)
        

!==RA=====Viscous diffusion term ===========================================================
!       Eq. = \partial <u''_j tau'_ki> / partial (x_k) + 
!             \partial <u''_i tau'_kj> / partial (x_k)
        DO J=1, NCL2
            !==============for TKE and MKE==================
            IF(J==1) THEN
                BUDG_viscs_diffu_TKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J+1)* ( Taufuf_RA(J+1,1, 2, 1) + &
                                         Taufuf_RA(J+1,2, 2, 2) + &
                                         Taufuf_RA(J+1,3, 2, 3) ) )-&
                    0.0_WP &
                    )*DYFI(J)
                BUDG_viscs_diffu_MKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,  1,2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,  2,2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,  3,2) ) &
                     +YCL2ND_WFB(J+1)* ( U1xztL_F0_io(J+1,1)*Tau_Mean_RA(J+1,1,2) + &
                                         U1xztL_F0_io(J+1,2)*Tau_Mean_RA(J+1,2,2) + &
                                         U1xztL_F0_io(J+1,3)*Tau_Mean_RA(J+1,3,2) ) )-&
                    0.0_WP &
                    )*DYFI(J)
            ELSE IF (J==NCL2) THEN
                BUDG_viscs_diffu_TKE(J) =  ( &
                    0.0_WP-&
                    ( YCL2ND_WFF(J)  * ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J)  * ( Taufuf_RA(J-1,1, 2, 1) + &
                                         Taufuf_RA(J-1,2, 2, 2) + &
                                         Taufuf_RA(J-1,3, 2, 3) ) ) &
                    )*DYFI(J)
                    
                BUDG_viscs_diffu_MKE(J) =  ( &
                    0.0_WP-&
                    ( YCL2ND_WFF(J)  * ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,1,  2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,2,  2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,3,  2) ) &
                     +YCL2ND_WFB(J)  * ( U1xztL_F0_io(J-1,1)*Tau_Mean_RA(J-1,1,2) + &
                                         U1xztL_F0_io(J-1,2)*Tau_Mean_RA(J-1,2,2) + &
                                         U1xztL_F0_io(J-1,3)*Tau_Mean_RA(J-1,3,2) ) ) &
                    )*DYFI(J)
            ELSE
                BUDG_viscs_diffu_TKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J+1)* ( Taufuf_RA(J+1,1, 2, 1) + &
                                         Taufuf_RA(J+1,2, 2, 2) + &
                                         Taufuf_RA(J+1,3, 2, 3) ) )-&
                    ( YCL2ND_WFF(J)  * ( Taufuf_RA(J,  1, 2, 1) + &
                                         Taufuf_RA(J,  2, 2, 2) + &
                                         Taufuf_RA(J,  3, 2, 3) ) &
                     +YCL2ND_WFB(J)  * ( Taufuf_RA(J-1,1, 2, 1) + &
                                         Taufuf_RA(J-1,2, 2, 2) + &
                                         Taufuf_RA(J-1,3, 2, 3) ) ) &
                    )*DYFI(J)
                    
                BUDG_viscs_diffu_MKE(J) =  ( &
                    ( YCL2ND_WFF(J+1)* ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,  1,2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,  2,2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,  3,2) ) &
                     +YCL2ND_WFB(J+1)* ( U1xztL_F0_io(J+1,1)*Tau_Mean_RA(J+1,1,2) + &
                                         U1xztL_F0_io(J+1,2)*Tau_Mean_RA(J+1,2,2) + &
                                         U1xztL_F0_io(J+1,3)*Tau_Mean_RA(J+1,3,2) ) )-&
                    ( YCL2ND_WFF(J)  * ( U1xztL_F0_io(J,  1)*Tau_Mean_RA(J,  1,2) + &
                                         U1xztL_F0_io(J,  2)*Tau_Mean_RA(J,  2,2) + &
                                         U1xztL_F0_io(J,  3)*Tau_Mean_RA(J,  3,2) ) &
                     +YCL2ND_WFB(J)  * ( U1xztL_F0_io(J-1,1)*Tau_Mean_RA(J-1,1,2) + &
                                         U1xztL_F0_io(J-1,2)*Tau_Mean_RA(J-1,2,2) + &
                                         U1xztL_F0_io(J-1,3)*Tau_Mean_RA(J-1,3,2) ) ) &
                    )*DYFI(J)
            END IF
            
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    
                    IF(J==1) THEN
                        BUDG_viscs_diffu_duiuj(J,L)= uf2_RA(J+1,M,N)*APVR(J,1) + &
                                                     uf2_RA(J,  M,N)*ACVR(J,1) + 0.0_WP*AMVR(J,1)
                    ELSE IF(J==NCL2) THEN
                        BUDG_viscs_diffu_duiuj(J,L)= 0.0_WP*APVR(J,1) + &
                                                     uf2_RA(J,  M,N)*ACVR(J,1) + &
                                                     uf2_RA(J-1,M,N)*AMVR(J,1)
                    ELSE
                        BUDG_viscs_diffu_duiuj(J,L)= uf2_RA(J+1,M,N)*APVR(J,1) + &
                                                     uf2_RA(J,  M,N)*ACVR(J,1) + &
                                                     uf2_RA(J-1,M,N)*AMVR(J,1) 
                    END IF
                    BUDG_viscs_diffu_duiuj(J,L) = BUDG_viscs_diffu_duiuj(J,L) *CVISC
                    
                END DO
            END DO
            
            BUDG_viscs_diffu_TKE(J) = 0.5_wp* ( BUDG_viscs_diffu_duiuj(J,1) + &
                                                BUDG_viscs_diffu_duiuj(J,4) + &
                                                BUDG_viscs_diffu_duiuj(J,6) )
            
        END DO
        CALL CHKHDL('      ==>Calculated BUDG_viscs_diffu_duiuj',myid)    
        
        !================BODY FORCE/ BUOYANCY PRODUCTION, which is part of pressure acceleration term===========
        !       ! not an independe contribution, but belongs to part of BUDG_press_accl1_duiuj
        BUDG_prodc_gvfc2_duiuj=0.0_wp  
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    IF(M.NE.GRAVDIR .AND. N.NE.GRAVDIR) CYCLE
                    
                    L = (M*(7-M))/2+N-3 
                       
                    IF(M.EQ.GRAVDIR .AND. N.EQ.GRAVDIR) THEN
                        COE = 2.0_WP 
                        K   = GRAVDIR    
                    ELSE IF (M.EQ.GRAVDIR) THEN
                        COE = 1.0_WP
                        K   = N
                    ELSE IF (N.EQ.GRAVDIR) THEN
                        COE = 1.0_WP
                        K   = M
                    ELSE
                        COE = 0.0_WP
                        K   = -1 !(WHICH WILL LEAD TO ERROR!)
                    END IF
                    ! F_A includes a postive or negtive sign
                    BUDG_prodc_gvfc2_duiuj(J,L)= F_A * COE * ( G1xztL_F0_io(J,K) - D1xztL_F0_io(J) * U1xztL_F0_io(J,K) )
                    
                    !du=( G1xztL_F0_io(J,1) - D1xztL_F0_io(J) * U1xztL_F0_io(J,1) )*F_A
                    !dv=( G1xztL_F0_io(J,2) - D1xztL_F0_io(J) * U1xztL_F0_io(J,2) )*F_A
                    !dw=( G1xztL_F0_io(J,3) - D1xztL_F0_io(J) * U1xztL_F0_io(J,3) )*F_A
            
                END DO
            END DO
            
        END DO
        BUDG_press_accl1_duiuj(:,:)=0.0_wp
        BUDG_viscs_accl1_duiuj(:,:)=0.0_wp    
        
        
        BUDG_prodc_dvfc1_duiuj=0.0_wp  
        DO J=1,NCL2
            ! f_1 * u"_1 + f_1 * u"_1
            BUDG_prodc_dvfc1_duiuj(J,1) = ( FUxztL_F0_io(J,1)-FUxztL_F0_io(J,4)*U1xztL_F0_io(J,1)  )*2.0_wp
            ! f_1 * u"_2
            BUDG_prodc_dvfc1_duiuj(J,2) =   FUxztL_F0_io(J,2)-FUxztL_F0_io(J,4)*U1xztL_F0_io(J,2)
            ! f_1 * u"_3
            BUDG_prodc_dvfc1_duiuj(J,3) =   FUxztL_F0_io(J,3)-FUxztL_F0_io(J,4)*U1xztL_F0_io(J,3)
            
            !==============for TKE and MKE==================
            !BUDG_prodc_dvfc1_TKE(J) = DrivenForce(J)*uff_RA(J,1) 
            !BUDG_prodc_dvfc1_MKE(J) = DrivenForce(J)*U_FA(J,1)
            
            !!WRITE(*,*) 'Driven force', DrivenForce(J), FUxztL_F0_io(J,4)
            !!WRITE(*,*) 'FU', FUxztL_F0_io(J,1:4) !test
            
            BUDG_prodc_dvfc1_TKE(J) = FUxztL_F0_io(J,1) - FUxztL_F0_io(J,4)*U1xztL_F0_io(J,1)
            !BUDG_prodc_dvfc1_MKE(J) = FUxztL_F0_io(J,4) * U1xztL_F0_io(J,1)
            BUDG_prodc_dvfc1_MKE(J) = DrivenForce(J)*U_FA(J,1)
            
        END DO
        
        !====RA===========BALANCE===================
        DO J=1, NCL2
            !==============for TKE and MKE==================
            BUDG_balance1_TKE(J) =     BUDG_prodc_stres_TKE(J) + &
                                      BUDG_viscs_dissp_TKE(J) + &
                                      BUDG_pdudx_stran_TKE(J) + &
                                      BUDG_Turbu_diffu_TKE(J) + &
                                      BUDG_dpudx_diffu_TKE(J) + &
                                      BUDG_viscs_diffu_TKE(J) + &
                                      BUDG_press_accl1_TKE(J) + &
                                      BUDG_viscs_accl1_TKE(J) + &
                                      BUDG_prodc_gvfc2_TKE(J) + &
                                      BUDG_prodc_dvfc1_TKE(J)
                                      
            BUDG_balance1_MKE(J) =     BUDG_prodc_stres_MKE(J) + &
                                      BUDG_viscs_dissp_MKE(J) + &
                                      BUDG_pdudx_stran_MKE(J) + &
                                      BUDG_Turbu_diffu_MKE(J) + &
                                      BUDG_dpudx_diffu_MKE(J) + &
                                      BUDG_viscs_diffu_MKE(J) + &
                                      BUDG_press_accl1_MKE(J) + &
                                      BUDG_viscs_accl1_MKE(J) + &
                                      BUDG_prodc_gvfc1_MKE(J) + &
                                      BUDG_prodc_dvfc1_MKE(J)
                                      
            !==========for each Ruv=========================                                  
            DO M=1, NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    BUDG_balance1_duiuj(J,L) = BUDG_prodc_stres_duiuj(J,L) + &
                                              BUDG_viscs_dissp_duiuj(J,L) + &
                                              BUDG_pdudx_stran_duiuj(J,L) + &
                                              BUDG_Turbu_diffu_duiuj(J,L) + &
                                              BUDG_dpudx_diffu_duiuj(J,L) + &
                                              BUDG_pdudx_stran_duiuj(J,L) + &
                                              BUDG_viscs_diffu_duiuj(J,L) + &
                                              BUDG_press_accl1_duiuj(J,L) + &
                                              BUDG_viscs_accl1_duiuj(J,L) + &
                                              BUDG_prodc_gvfc2_duiuj(J,L) + &
                                              BUDG_prodc_dvfc1_duiuj(J,L)
                                              
                    
                END DO
            END DO
        END DO
        !================integral of each terms along Y==================================
        DO J=1, NCL2
            BUDG_prodc_stres_TKE_ysum = BUDG_prodc_stres_TKE_ysum + BUDG_prodc_stres_TKE(J)/DYFI(J)
            BUDG_viscs_dissp_TKE_ysum = BUDG_viscs_dissp_TKE_ysum + BUDG_viscs_dissp_TKE(J)/DYFI(J)
            BUDG_pdudx_stran_TKE_ysum = BUDG_pdudx_stran_TKE_ysum + BUDG_pdudx_stran_TKE(J)/DYFI(J)
            BUDG_Turbu_diffu_TKE_ysum = BUDG_Turbu_diffu_TKE_ysum + BUDG_Turbu_diffu_TKE(J)/DYFI(J)
            BUDG_dpudx_diffu_TKE_ysum = BUDG_dpudx_diffu_TKE_ysum + BUDG_dpudx_diffu_TKE(J)/DYFI(J)
            BUDG_viscs_diffu_TKE_ysum = BUDG_viscs_diffu_TKE_ysum + BUDG_viscs_diffu_TKE(J)/DYFI(J)
            BUDG_press_accl1_TKE_ysum = BUDG_press_accl1_TKE_ysum + BUDG_press_accl1_TKE(J)/DYFI(J)
            BUDG_viscs_accl1_TKE_ysum = BUDG_viscs_accl1_TKE_ysum + BUDG_viscs_accl1_TKE(J)/DYFI(J)
            BUDG_prodc_gvfc2_TKE_ysum = BUDG_prodc_gvfc2_TKE_ysum + BUDG_prodc_gvfc2_TKE(J)/DYFI(J)
            BUDG_prodc_dvfc1_TKE_ysum = BUDG_prodc_dvfc1_TKE_ysum + BUDG_prodc_dvfc1_TKE(J)/DYFI(J)
            BUDG_balance1_TKE_ysum     = BUDG_balance1_TKE_ysum     + BUDG_balance1_TKE(J)/DYFI(J)
        
            BUDG_prodc_stres_MKE_ysum = BUDG_prodc_stres_MKE_ysum + BUDG_prodc_stres_MKE(J)/DYFI(J)
            BUDG_viscs_dissp_MKE_ysum = BUDG_viscs_dissp_MKE_ysum + BUDG_viscs_dissp_MKE(J)/DYFI(J)
            BUDG_pdudx_stran_MKE_ysum = BUDG_pdudx_stran_MKE_ysum + BUDG_pdudx_stran_MKE(J)/DYFI(J)
            BUDG_Turbu_diffu_MKE_ysum = BUDG_Turbu_diffu_MKE_ysum + BUDG_Turbu_diffu_MKE(J)/DYFI(J)
            BUDG_dpudx_diffu_MKE_ysum = BUDG_dpudx_diffu_MKE_ysum + BUDG_dpudx_diffu_MKE(J)/DYFI(J)
            BUDG_viscs_diffu_MKE_ysum = BUDG_viscs_diffu_MKE_ysum + BUDG_viscs_diffu_MKE(J)/DYFI(J)
            BUDG_press_accl1_MKE_ysum = BUDG_press_accl1_MKE_ysum + BUDG_press_accl1_MKE(J)/DYFI(J)
            BUDG_viscs_accl1_MKE_ysum = BUDG_viscs_accl1_MKE_ysum + BUDG_viscs_accl1_MKE(J)/DYFI(J)
            BUDG_prodc_gvfc1_MKE_ysum = BUDG_prodc_gvfc1_MKE_ysum + BUDG_prodc_gvfc1_MKE(J)/DYFI(J)
            BUDG_prodc_dvfc1_MKE_ysum = BUDG_prodc_dvfc1_MKE_ysum + BUDG_prodc_dvfc1_MKE(J)/DYFI(J)
            BUDG_balance1_MKE_ysum     = BUDG_balance1_MKE_ysum     + BUDG_balance1_MKE(J)/DYFI(J)
            
            DO L=1, (NDV*(7-NDV))/2+NDV-3
                BUDG_prodc_stres_duiuj_ysum(L) = BUDG_prodc_stres_duiuj_ysum(L) + BUDG_prodc_stres_duiuj(J,L)/DYFI(J)
                BUDG_viscs_dissp_duiuj_ysum(L) = BUDG_viscs_dissp_duiuj_ysum(L) + BUDG_viscs_dissp_duiuj(J,L)/DYFI(J)
                BUDG_pdudx_stran_duiuj_ysum(L) = BUDG_pdudx_stran_duiuj_ysum(L) + BUDG_pdudx_stran_duiuj(J,L)/DYFI(J)
                BUDG_Turbu_diffu_duiuj_ysum(L) = BUDG_Turbu_diffu_duiuj_ysum(L) + BUDG_Turbu_diffu_duiuj(J,L)/DYFI(J)
                BUDG_dpudx_diffu_duiuj_ysum(L) = BUDG_dpudx_diffu_duiuj_ysum(L) + BUDG_dpudx_diffu_duiuj(J,L)/DYFI(J)
                BUDG_viscs_diffu_duiuj_ysum(L) = BUDG_viscs_diffu_duiuj_ysum(L) + BUDG_viscs_diffu_duiuj(J,L)/DYFI(J)
                BUDG_press_accl1_duiuj_ysum(L) = BUDG_press_accl1_duiuj_ysum(L) + BUDG_press_accl1_duiuj(J,L)/DYFI(J)
                BUDG_viscs_accl1_duiuj_ysum(L) = BUDG_viscs_accl1_duiuj_ysum(L) + BUDG_viscs_accl1_duiuj(J,L)/DYFI(J)
                BUDG_prodc_gvfc2_duiuj_ysum(L) = BUDG_prodc_gvfc2_duiuj_ysum(L) + BUDG_prodc_gvfc2_duiuj(J,L)/DYFI(J)
                BUDG_prodc_dvfc1_duiuj_ysum(L) = BUDG_prodc_dvfc1_duiuj_ysum(L) + BUDG_prodc_dvfc1_duiuj(J,L)/DYFI(J)
                BUDG_balance1_duiuj_ysum(L)     = BUDG_balance1_duiuj_ysum(L)     + BUDG_balance1_duiuj(J,L)/DYFI(J)
            END DO
            
        END DO
        
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        FLNM=TRIM(filepath4)//'Result.IO.budget.check.'//TRIM(PNTIM)//'.tec'
        INQUIRE(FILE=TRIM(ADJUSTL(FLNM)), EXIST=file_exists) 
        if(file_exists) then
            OPEN(TECFLG,FILE=FLNM, POSITION='APPEND')
        else
            OPEN(TECFLG,FILE=FLNM)
        end if
        WRITE(TECFLG,'(A)') '=======RA=====uu,uv,uw,vv,vw,ww=========='
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_prodc_stres_duiuj_ysum(1:6) = ', BUDG_prodc_stres_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_viscs_dissp_duiuj_ysum(1:6) = ', BUDG_viscs_dissp_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_pdudx_stran_duiuj_ysum(1:6) = ', BUDG_pdudx_stran_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_Turbu_diffu_duiuj_ysum(1:6) = ', BUDG_Turbu_diffu_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_dpudx_diffu_duiuj_ysum(1:6) = ', BUDG_dpudx_diffu_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_viscs_diffu_duiuj_ysum(1:6) = ', BUDG_viscs_diffu_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_press_accl1_duiuj_ysum(1:6) = ', BUDG_press_accl1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_viscs_accl1_duiuj_ysum(1:6) = ', BUDG_viscs_accl1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_prodc_gvfc2_duiuj_ysum(1:6) = ', BUDG_prodc_gvfc2_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_prodc_dvfc1_duiuj_ysum(1:6) = ', BUDG_prodc_dvfc1_duiuj_ysum(1:6)
        WRITE(TECFLG,'(A,6ES20.7)') 'BUDG_balance1_duiuj_ysum(1:6) =     ', BUDG_balance1_duiuj_ysum(1:6)
        
        WRITE(TECFLG,'(A)') '=======RA=====TKE=========='
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_prodc_stres_TKE_ysum = ', BUDG_prodc_stres_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_viscs_dissp_TKE_ysum = ', BUDG_viscs_dissp_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_pdudx_stran_TKE_ysum = ', BUDG_pdudx_stran_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_Turbu_diffu_TKE_ysum = ', BUDG_Turbu_diffu_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_dpudx_diffu_TKE_ysum = ', BUDG_dpudx_diffu_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_viscs_diffu_TKE_ysum = ', BUDG_viscs_diffu_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_press_accl1_TKE_ysum = ', BUDG_press_accl1_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_viscs_accl1_TKE_ysum = ', BUDG_viscs_accl1_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_prodc_gvfc2_TKE_ysum = ', BUDG_prodc_gvfc2_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_prodc_dvfc1_TKE_ysum = ', BUDG_prodc_dvfc1_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_balance1_TKE_ysum =     ', BUDG_balance1_TKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'Sum of all terms =          ', BUDG_prodc_stres_TKE_ysum + BUDG_viscs_dissp_TKE_ysum + &
                BUDG_pdudx_stran_TKE_ysum + BUDG_Turbu_diffu_TKE_ysum + BUDG_dpudx_diffu_TKE_ysum + BUDG_viscs_diffu_TKE_ysum &
              + BUDG_press_accl1_TKE_ysum + BUDG_viscs_accl1_TKE_ysum + BUDG_prodc_gvfc2_TKE_ysum + BUDG_prodc_dvfc1_TKE_ysum                 
        WRITE(TECFLG,'(A)') '======RA======MKE=========='
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_prodc_stres_MKE_ysum = ', BUDG_prodc_stres_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_viscs_dissp_MKE_ysum = ', BUDG_viscs_dissp_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_pdudx_stran_MKE_ysum = ', BUDG_pdudx_stran_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_Turbu_diffu_MKE_ysum = ', BUDG_Turbu_diffu_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_dpudx_diffu_MKE_ysum = ', BUDG_dpudx_diffu_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_viscs_diffu_MKE_ysum = ', BUDG_viscs_diffu_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_press_accl1_MKE_ysum = ', BUDG_press_accl1_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_viscs_accl1_MKE_ysum = ', BUDG_viscs_accl1_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_prodc_gvfc1_MKE_ysum = ', BUDG_prodc_gvfc1_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_prodc_dvfc1_MKE_ysum = ', BUDG_prodc_dvfc1_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'BUDG_balance1_MKE_ysum =     ', BUDG_balance1_MKE_ysum
        WRITE(TECFLG,'(A,1ES20.7)') 'Sum of all terms =          ', BUDG_prodc_stres_MKE_ysum + BUDG_viscs_dissp_MKE_ysum + &
                BUDG_pdudx_stran_MKE_ysum + BUDG_Turbu_diffu_MKE_ysum + BUDG_dpudx_diffu_MKE_ysum + BUDG_viscs_diffu_MKE_ysum &
              + BUDG_press_accl1_MKE_ysum + BUDG_viscs_accl1_MKE_ysum + BUDG_prodc_gvfc1_MKE_ysum + BUDG_prodc_dvfc1_MKE_ysum                                         
                                            
        CLOSE(TECFLG)

        RETURN
    END SUBROUTINE
    
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    SUBROUTINE PP_HEAT_BASIC_VARS_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: J
        INTEGER(4) :: M, N, H, P, L, K
        REAL(WP)    :: COE
        REAL(WP)    :: htaujc, htaujp
        
        !=============={h}===================================
        ! Eq. {h} =<\rho h>/<\rho>
        H_FA = 0.0_WP
        DO J=1,NCL2
             H_FA(J)   = DHxztL_F0_io(J)/D1xztL_F0_io(J)
             hff_RA(J) = H1xztL_F0_io(J) - H_FA(J)
        END DO
        H_FA(0)    = HWAL_FA(ibotwall)
        H_FA(NND2) = HWAL_FA(itopwall)
        
        CALL CHKHDL('      ==>Calculated {h} and <h">',myid)
        
        !====<p' h'>=<p' h''>=<ph> - <p>*<h>==hp_per_RA(J,i)=============
        DO J=1,NCL2
            hfpf_RA(J)=  PHxztL_F0_io(J) - U1xztL_F0_io(J,  4) * H1xztL_F0_io(J)
        END DO
        CALL CHKHDL('      ==>Calculated <h`p`>',myid)
        
        !==============d{h}/dx_m= dHdX_FA(CL,m)========================
        !==============d<h>/dx_m= dHdX_RA(CL,m)========================
        !==============d<T>/dx_m= dTdX(CL,m)========================
        M=2
        dHdX_FA = 0.0_WP
        dHdX_RA = 0.0_WP
        dTdX    = 0.0_WP
        dDdX    = 0.0_WP
        DO J=1,NCL2
          
            IF(J==1) THEN
                dHdX_RA(J,M)= ( ( YCL2ND_WFB(J+1)*H1xztL_F0_io(J)  + &
                                  YCL2ND_WFF(J+1)*H1xztL_F0_io(J+1) ) - &
                                  HWAL_RA(ibotwall)  ) * DYFI(J)
                                  
                dHdX_FA(J,M)= ( ( YCL2ND_WFB(J+1)*H_FA(J)  + &
                                  YCL2ND_WFF(J+1)*H_FA(J+1) ) - &
                                  H_FA(0)  ) * DYFI(J)
                                  
                dTdX(J,M)= ( ( YCL2ND_WFB(J+1)*T1xztL_F0_io(J)  + &
                                        YCL2ND_WFF(J+1)*T1xztL_F0_io(J+1) ) - &
                                        TWAL(ibotwall)  ) * DYFI(J)
                dDdX(J,M)= ( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J)  + &
                                        YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) - &
                                        DWAL(ibotwall)  ) * DYFI(J)
                                  
            ELSE IF (J==NCL2) THEN
                dHdX_RA(J,M)= (  HWAL_RA(itopwall) - &
                               ( YCL2ND_WFF(J)*H1xztL_F0_io(J) + &
                                 YCL2ND_WFB(J)*H1xztL_F0_io(J-1) )  ) * DYFI(J)
                                 
                dHdX_FA(J,M)= (  H_FA(NND2) - &
                               ( YCL2ND_WFF(J)*H_FA(J) + &
                                 YCL2ND_WFB(J)*H_FA(J-1) )  ) * DYFI(J)
                
                dTdX(J,M)= (  TWAL(itopwall) - &
                                    ( YCL2ND_WFF(J)*T1xztL_F0_io(J) + &
                                      YCL2ND_WFB(J)*T1xztL_F0_io(J-1) )  ) * DYFI(J)
                    
                dDdX(J,M)= (  DWAL(itopwall) - &
                                    ( YCL2ND_WFF(J)*D1xztL_F0_io(J) + &
                                      YCL2ND_WFB(J)*D1xztL_F0_io(J-1) )  ) * DYFI(J)
                              
            ELSE
                dHdX_RA(J,M)= ( ( YCL2ND_WFB(J+1)*H1xztL_F0_io(J  ) + &
                                  YCL2ND_WFF(J+1)*H1xztL_F0_io(J+1) ) - &
                                ( YCL2ND_WFF(J)  *H1xztL_F0_io(J  ) + &
                                  YCL2ND_WFB(J)  *H1xztL_F0_io(J-1) ) ) * DYFI(J)
                                  
                dHdX_FA(J,M)= ( ( YCL2ND_WFB(J+1)*H_FA(J  ) + &
                                  YCL2ND_WFF(J+1)*H_FA(J+1) ) - &
                                ( YCL2ND_WFF(J)  *H_FA(J  ) + &
                                  YCL2ND_WFB(J)  *H_FA(J-1) ) ) * DYFI(J)
                                  
                dTdX(J,M)= ( ( YCL2ND_WFB(J+1)*T1xztL_F0_io(J  ) + &
                               YCL2ND_WFF(J+1)*T1xztL_F0_io(J+1) ) - &
                             ( YCL2ND_WFF(J)  *T1xztL_F0_io(J  ) + &
                               YCL2ND_WFB(J)  *T1xztL_F0_io(J-1) ) ) * DYFI(J)
                                       
                dDdX(J,M)= ( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J  ) + &
                               YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) - &
                             ( YCL2ND_WFF(J)  *D1xztL_F0_io(J  ) + &
                               YCL2ND_WFB(J)  *D1xztL_F0_io(J-1) ) ) * DYFI(J)
                                       
            END IF
           
        END DO
        CALL CHKHDL('      ==>Calculated d<h>/dy, d<T>/dy and d{h}/dy',myid)
        

        !==============<\rho>{h'' u''_m} = <\rho h'' u''_m> = uffhffd_FA(CL,m)========================
        ! Eq: uffhffd_FA(J,M) = <\rho>{h'' u''_m} = <\rho u_m h> - <\rho u_m> * <\rho h>/ <\rho>
        !     UH_FA(J,M) = {u_m h} = <\rho u_m h>/ <\rho>
        DO J=1,NCL2
            DO M=1,NDV
                uffhffd_FA(J,M) = GHxztL_F0_io(J,M) - G1xztL_F0_io(J,M)*DHxztL_F0_io(J)/ D1xztL_F0_io(J)
                ufhfd_RA(J,M) =(UHxztL_F0_io(J,M) - U1xztL_F0_io(J,M)*H1xztL_F0_io(J))*D1xztL_F0_io(J)
                UH_FA(J,M)         = GHxztL_F0_io(J,M)/D1xztL_F0_io(J)
            END DO
        END DO
        CALL CHKHDL('      ==>Calculated <uh>, uffhffd_FA, ufhfd_RA',myid)
        
        !===<\rho u''_k u''_i h''> = <\rho> {u''_k u''_i h''}==============================
        ! Eq. : 
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    L = (M*(7-M))/2+N-3
                    uff2hffd_FA(J,M,N) = U2DHxztL_F0_io(J,L)                        &
                                        -D1xztL_F0_io(J) * H_FA(J)   * UU_FA(J,M,N) &
                                        -D1xztL_F0_io(J) * U_FA(J,N) * UH_FA(J,M) &
                                        -D1xztL_F0_io(J) * U_FA(J,M) * UH_FA(J,N) &
                                        +2.0_WP * D1xztL_F0_io(J) * H_FA(J) * U_FA(J,N) * U_FA(J,M)
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                        uff2hffd_FA(J,M,N) = uff2hffd_FA(J,N,M)
                    END IF
                END DO
            END DO
            
        END DO
        
        !========ViscStressEnth_RA(J, M, N)=<h \tau_mn>======================
        !Eq.
        DO J=1, NCL2    
            DO M=1,NDV
                DO N=1, NDV
                    IF(M.GT.N) CYCLE
                    ViscStressEnth_RA(J, M, N)= ( & 
                                        DVDL1MHxztL_F0_io(J,M,N) + &
                                        DVDL1MHxztL_F0_io(J,N,M) - &
                        2.0_WP/3.0_WP*( DVDL1MHxztL_F0_io(J,1,1) + &
                                        DVDL1MHxztL_F0_io(J,2,2) + &
                                        DVDL1MHxztL_F0_io(J,3,3) )*DBLE(Kronecker_delta(M,N)) )/REN
                END DO
            END DO
            
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) THEN
                        ViscStressEnth_RA(J, M, N) = ViscStressEnth_RA(J, N, M)
                    END IF
                END DO
            END DO
            
        END DO 
        
        !========<\partial(h)/\partial x_m \tau_hp>======================
        ! Eq. ViscStressEnthGrad_RA(J, M, H, P)
        !     = <mu * d(u_h)/d(x_p) * d(h)/d(x_m)> +
        !       <mu * d(u_p)/d(x_h) * d(h)/d(x_m)> - 2/3*
        !       <mu * d(u_l)/d(x_l) * d(h)/d(x_m)>*delta_hp
        DO J=1, NCL2    
            DO M=1,NDV      
               
                DO H=1, NDV           
                    DO P=1, NDV    
                        IF(H.GT.P) CYCLE
                        ViscStressEnthGrad_RA(J, M, H, P) = (     &
                            DHDLMDVDLxztL_F0_io(J, M, H, P ) + &
                            DHDLMDVDLxztL_F0_io(J, M, P, H ) - &
                            2.0_wp/3.0_wp* ( & 
                            DHDLMDVDLxztL_F0_io(J, M, 1, 1 ) + &
                            DHDLMDVDLxztL_F0_io(J, M, 2, 2 ) + &
                            DHDLMDVDLxztL_F0_io(J, M, 3, 3 )   &
                            )*DBLE(Kronecker_delta(H,P)) )/REN
                    END DO
                END DO
                
                DO H=1, NDV           
                    DO P=1, NDV    
                        IF(H.LT.P) CYCLE
                        ViscStressEnthGrad_RA(J, M, H, P) = ViscStressEnthGrad_RA(J, M, P, H)
                    END DO
                END DO
                
            END DO
        END DO 
    
    RETURN
    END SUBROUTINE 
        


    SUBROUTINE PP_HEAT_FA_RSTE_BUDG_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: J
        INTEGER(4) :: M, N, H, L, K
        REAL(WP) :: COE
        REAL(WP) :: htaujc, htaujp
        
    !==PRODUCTION TERMS due to mean shear and mean enthalpy gradient===(THF(turbulent heat flux))========    
       !Eq: BUDG_prod_thf_stres(J,L) = -<\rho h''   u''_k> (\partial {u_i} / \partial x_k)  
       !Eq: BUDG_prod_thf_entpg(J,L) = -<\rho u_i'' u''_k> (\partial {h}   / \partial x_k)  
        BUDG_prodc_stres_thf = 0.0_wp
        BUDG_prodc_enthg_thf = 0.0_WP
        DO J=1,NCL2
            DO M=1,NDV
                BUDG_prodc_stres_thf(J,M) =( uffhffd_FA(J,1) * dUdX_FA(J,M,1) + &
                                             uffhffd_FA(J,2) * dUdX_FA(J,M,2) + &
                                             uffhffd_FA(J,3) * dUdX_FA(J,M,3) )*(-1.0_wp)
                BUDG_prodc_enthg_thf(J,M) =( uff2d_FA(J,1,M) * dHdX_FA(J,1) + &
                                             uff2d_FA(J,2,M) * dHdX_FA(J,2) + &
                                             uff2d_FA(J,3,M) * dHdX_FA(J,3) )*(-1.0_wp)
            END DO
        END DO
        
       !==TURBUELCEN DIFFUSION TERMS=(turblence tranport rate)================================
       !Eq.  BUDG_Tdiff_duiuj(J,L) = ( \partial <\rho u''_i u''_j u''_k > ) / (\partial x_k )
        DO J=1,NCL2   
            DO M=1,NDV
                IF(J==1) THEN
                    BUDG_Turbu_diffu_thf(J,M)= &
                        ( ( YCL2ND_WFF(J+1)*uff2hffd_FA(J+1,M,2) + &
                            YCL2ND_WFB(J+1)*uff2hffd_FA(J,  M,2) ) - 0.0_WP ) * DYFI(J) *(-1.0_wp)
                ELSE IF (J==NCL2) THEN
                    BUDG_Turbu_diffu_thf(J,M)= &
                        ( 0.0_WP - &
                        ( YCL2ND_WFF(J) *uff2hffd_FA(J,  M,2) + &
                          YCL2ND_WFB(J) *uff2hffd_FA(J-1,M,2) ) ) * DYFI(J) *(-1.0_wp)
                                      
                ELSE
                    BUDG_Turbu_diffu_thf(J,M)= &
                        ( ( YCL2ND_WFF(J+1)*uff2hffd_FA(J+1,M,2) + &
                            YCL2ND_WFB(J+1)*uff2hffd_FA(J,  M,2) ) - &
                          ( YCL2ND_WFF(J)  *uff2hffd_FA(J,  M,2) + &
                            YCL2ND_WFB(J)  *uff2hffd_FA(J-1,M,2) ) ) * DYFI(J) *(-1.0_wp)
                END IF
            END DO
        END DO
        
        !==================pressure acceleration term===============
        DO J=1, NCL2
            DO M=1,NDV
                BUDG_press_accl1_thf(J,M) = -dPdX_RA(J,M)*hff_RA(J) 
            END DO
        END DO
        
        !==================enthalpy-pressure gradient diffusion===============
        DO J=1, NCL2
            DO M=1,NDV
                IF (M==2) THEN
                    IF(J==1) THEN
                        BUDG_dphdx_diffu_thf(J,M) = &
                            ( ( YCL2ND_WFF(J+1)*hfpf_RA(J+1) + &
                                YCL2ND_WFB(J+1)*hfpf_RA(J) ) -  &
                                0.0_wp )*DYFI(J)*(-1.0_wp)
                    ELSE IF(J==NCL2) THEN
                        BUDG_dphdx_diffu_thf(J,M) = &
                            (  0.0_wp -  &
                              ( YCL2ND_WFF(J)  *hfpf_RA(J) + &
                                YCL2ND_WFB(J)  *hfpf_RA(J-1) ) )*DYFI(J)*(-1.0_wp)
                    ELSE
                        BUDG_dphdx_diffu_thf(J,M) = &
                            ( ( YCL2ND_WFF(J+1)*hfpf_RA(J+1) + &
                                YCL2ND_WFB(J+1)*hfpf_RA(J) ) -  &
                              ( YCL2ND_WFF(J)  *hfpf_RA(J) + &
                                YCL2ND_WFB(J)  *hfpf_RA(J-1) ) )*DYFI(J)*(-1.0_wp)  
                    END IF
                ELSE
                    BUDG_dphdx_diffu_thf(J,M) = 0.0_wp
                END IF
            END DO
        END DO
            
        !==================enthalpy-pressure strain===============
        DO J=1,NCL2
            DO M=1, NDV
                BUDG_pdhdx_stran_thf(J,M) =  DHDLPxztL_F0_io(J,M) - U1xztL_F0_io(J,4)*DHDL1xztL_F0_io(J,M) 
            END DO
        END DO
    
        !================conductive heat flux acceleration=============
        DO J=1, NCL2
            DO M=1, NDV
                IF(J==1) THEN
                    BUDG_ConHF_accel_thf(J,M) = uff_RA(J,M)* &
                            ( ( YCL2ND_WFF(J+1)*DTDLKxztL_F0_io(J+1,2) + &
                                YCL2ND_WFB(J+1)*DTDLKxztL_F0_io(J,  2) ) -  &
                                qw(ibotwall) )*DYFI(J)*CTHECD
                
                ELSE IF(J==NCL2) THEN
                    BUDG_ConHF_accel_thf(J,M) = uff_RA(J,M)* &
                            (  qw(itopwall) -  &
                              ( YCL2ND_WFF(J)  *DTDLKxztL_F0_io(J,  2) + &
                                YCL2ND_WFB(J)  *DTDLKxztL_F0_io(J-1,2) ) )*DYFI(J)*CTHECD
                ELSE
                    BUDG_ConHF_accel_thf(J,M) = uff_RA(J,M)* &
                            ( ( YCL2ND_WFF(J+1)*DTDLKxztL_F0_io(J+1,2) + &
                                YCL2ND_WFB(J+1)*DTDLKxztL_F0_io(J,  2) ) -  &
                              ( YCL2ND_WFF(J)  *DTDLKxztL_F0_io(J,  2) + &
                                YCL2ND_WFB(J)  *DTDLKxztL_F0_io(J-1,2) ) )*DYFI(J)*CTHECD
                END IF
                
            END DO
        END DO
        
        !================conductive heat flux diffusion=============
        DO J=1, NCL2
            DO M=1, NDV
                IF(J==1) THEN
                    BUDG_ConHF_diffu_thf(J,M) = &
                        ( ( YCL2ND_WFF(J+1)*( DTDLKUxztL_F0_io(J+1,2,M) - U1xztL_F0_io(J+1,M)*DTDLKxztL_F0_io(J+1,2) ) + &
                            YCL2ND_WFB(J+1)*( DTDLKUxztL_F0_io(J  ,2,M) - U1xztL_F0_io(J  ,M)*DTDLKxztL_F0_io(J,  2) ) )-&
                            0.0_wp )&
                        *DYFI(J)*CTHECD
                ELSE IF(J==NCL2) THEN
                    BUDG_ConHF_diffu_thf(J,M) = &
                        ( 0.0_wp-&
                          ( YCL2ND_WFF(J)  *( DTDLKUxztL_F0_io(J  ,2,M) - U1xztL_F0_io(J  ,M)*DTDLKxztL_F0_io(J,  2) ) + &
                            YCL2ND_WFB(J)  *( DTDLKUxztL_F0_io(J-1,2,M) - U1xztL_F0_io(J-1,M)*DTDLKxztL_F0_io(J-1,2) ) ) )&
                        *DYFI(J)*CTHECD
                ELSE
                    BUDG_ConHF_diffu_thf(J,M) = &
                        ( ( YCL2ND_WFF(J+1)*( DTDLKUxztL_F0_io(J+1,2,M) - U1xztL_F0_io(J+1,M)*DTDLKxztL_F0_io(J+1,2) ) + &
                            YCL2ND_WFB(J+1)*( DTDLKUxztL_F0_io(J  ,2,M) - U1xztL_F0_io(J  ,M)*DTDLKxztL_F0_io(J,  2) ) )-&
                          ( YCL2ND_WFF(J)  *( DTDLKUxztL_F0_io(J  ,2,M) - U1xztL_F0_io(J  ,M)*DTDLKxztL_F0_io(J,  2) ) + &
                            YCL2ND_WFB(J)  *( DTDLKUxztL_F0_io(J-1,2,M) - U1xztL_F0_io(J-1,M)*DTDLKxztL_F0_io(J-1,2) ) ) )&
                        *DYFI(J)*CTHECD
                END IF
            END DO
        END DO
    
        !============conductive heat flux dissipation==================
        DO J=1, NCL2
            DO M=1, NDV
                 BUDG_ConHF_dissp_thf(J,M) = -1.0_WP* ( dTdLKdVdLxztL_F0_io(J,1,M,1) + &
                                                        dTdLKdVdLxztL_F0_io(J,2,M,2) + &
                                                        dTdLKdVdLxztL_F0_io(J,3,M,3) )*CTHECD + &
                                                ( DTDLKxztL_F0_io(J,1) * DVDL1xztL_F0_io(J,M,1) + &
                                                  DTDLKxztL_F0_io(J,2) * DVDL1xztL_F0_io(J,M,2) + &
                                                  DTDLKxztL_F0_io(J,3) * DVDL1xztL_F0_io(J,M,3) )*CTHECD
            END DO
        END DO
    
        !=======Viscous acceleration, including all========(based on tau_RA)===================================
        !       Eq. <h''> (\partial <tau_ki>)/ (\partial  x_k) 
        DO J=1, NCL2
            DO M=1, NDV
                BUDG_viscs_accl1_thf(J,M) = dTaudy_RA(J,2,M)*hff_RA(J)
            END DO
        END DO
        
        !=======Viscous diffusion term ===========================================================
        !       Eq. = \partial <h'' tau'_ki> / partial (x_k)
        DO J=1, NCL2
            DO M=1, NDV
                IF(J==1) THEN
                    htaujp = YCL2ND_WFF(J+1)*( ViscStressEnth_RA(J,  2, M)-H1xztL_F0_io(J  )*Tau_Mean_RA(J,  2,M))+&
                             YCL2ND_WFB(J+1)*( ViscStressEnth_RA(J+1,2, M)-H1xztL_F0_io(J+1)*Tau_Mean_RA(J+1,2,M))
                    htaujc = ViscStressEnth_RA(J,  2, M)-H1xztL_F0_io(J  )*Tau_Mean_RA(J,  2,M)
                    BUDG_viscs_diffu_thf(J,M) = (htaujp-htaujc)*DYFI(J)*2.0_wp
                ELSE IF(J==NCL2) THEN
                    htaujp =  ViscStressEnth_RA(J,  2, M)-H1xztL_F0_io(J  )*Tau_Mean_RA(J,  2,M)
                    htaujc = YCL2ND_WFF(J)  *( ViscStressEnth_RA(J,  2, M)-H1xztL_F0_io(J  )*Tau_Mean_RA(J,  2,M))+&
                             YCL2ND_WFB(J)  *( ViscStressEnth_RA(J-1,2, M)-H1xztL_F0_io(J-1)*Tau_Mean_RA(J-1,2,M))
                    BUDG_viscs_diffu_thf(J,M) = (htaujp-htaujc)*DYFI(J)*2.0_wp 
                ELSE
                    htaujp = YCL2ND_WFF(J+1)*( ViscStressEnth_RA(J,  2, M)-H1xztL_F0_io(J  )*Tau_Mean_RA(J,  2,M))+&
                             YCL2ND_WFB(J+1)*( ViscStressEnth_RA(J+1,2, M)-H1xztL_F0_io(J+1)*Tau_Mean_RA(J+1,2,M))
                    htaujc = YCL2ND_WFF(J)  *( ViscStressEnth_RA(J,  2, M)-H1xztL_F0_io(J  )*Tau_Mean_RA(J,  2,M))+&
                             YCL2ND_WFB(J)  *( ViscStressEnth_RA(J-1,2, M)-H1xztL_F0_io(J-1)*Tau_Mean_RA(J-1,2,M)) 
                    BUDG_viscs_diffu_thf(J,M) = (htaujp-htaujc)*DYFI(J) 
                END IF
                
            END DO
        END DO
    
        !=======Viscous ENERGY dissipation term ===========================================================
        DO J=1, NCL2    
            DO M=1, NDV
                BUDG_viscs_dissp_thf(J,M) = -1.0_wp * ( &
                    ViscStressEnthGrad_RA(J, 1, 1, M) + &
                    ViscStressEnthGrad_RA(J, 2, 2, M) + &
                    ViscStressEnthGrad_RA(J, 3, 3, M) - &
                    Tau_Mean_RA(J,1,M) * DHDX_RA(J,1) - &
                    Tau_Mean_RA(J,2,M) * DHDX_RA(J,2) - &
                    Tau_Mean_RA(J,3,M) * DHDX_RA(J,3) )
            END DO
        END DO
        
        DO J=1, NCL2
            DO M=1, NDV
                BUDG_balance1_thf(J,M) = BUDG_prodc_stres_thf(J,M)+BUDG_prodc_enthg_thf(J,M)+ &
                                        BUDG_Turbu_diffu_thf(J,M)+BUDG_press_accl1_thf(J,M)+ &
                                        BUDG_dphdx_diffu_thf(J,M)+BUDG_pdhdx_stran_thf(J,M)+ &
                                        BUDG_ConHF_accel_thf(J,M)+BUDG_ConHF_diffu_thf(J,M)+ &
                                        BUDG_ConHF_dissp_thf(J,M)+BUDG_viscs_accl1_thf(J,M)+ &
                                        BUDG_viscs_diffu_thf(J,M)+BUDG_viscs_dissp_thf(J,M)
            END DO
        END DO
    
        
        !================BODY FORCE/ BUOYANCY PRODUCTION, which is part of pressure acceleration term===========
        !       ! not an independe contribution, but belongs to part of BUDG_press_accl1_duiuj
        BUDG_prodc_gvfc2_duiuj=0.0_wp  
        DO J=1,NCL2
        
            DO M=1,NDV
                DO N=1,NDV
                    IF(M.GT.N) CYCLE
                    IF(M.NE.GRAVDIR .AND. N.NE.GRAVDIR) CYCLE
                    
                    L = (M*(7-M))/2+N-3 
                       
                    IF(M.EQ.GRAVDIR .AND. N.EQ.GRAVDIR) THEN
                        COE = 2.0_WP 
                        K   = GRAVDIR    
                    ELSE IF (M.EQ.GRAVDIR) THEN
                        COE = 1.0_WP
                        K   = N
                    ELSE IF (N.EQ.GRAVDIR) THEN
                        COE = 1.0_WP
                        K   = M
                    ELSE
                        COE = 0.0_WP
                        K   = -1 !(WHICH WILL LEAD TO ERROR!)
                    END IF
                    ! F_A includes a postive or negtive sign
                    BUDG_prodc_gvfc2_duiuj(J,L)= F_A * COE * ( G1xztL_F0_io(J,K) - D1xztL_F0_io(J) * U1xztL_F0_io(J,K) )
                    
                    !du=( G1xztL_F0_io(J,1) - D1xztL_F0_io(J) * U1xztL_F0_io(J,1) )*F_A
                    !dv=( G1xztL_F0_io(J,2) - D1xztL_F0_io(J) * U1xztL_F0_io(J,2) )*F_A
                    !dw=( G1xztL_F0_io(J,3) - D1xztL_F0_io(J) * U1xztL_F0_io(J,3) )*F_A
            
                END DO
            END DO
            
            IF(GRAVDIR==0) THEN
                BUDG_prodc_gvfc2_thf(J) = 0.0_wp
            ELSE
                BUDG_prodc_gvfc2_thf(J) = F_A * ( DHxztL_F0_io(J) - D1xztL_F0_io(J) * H1xztL_F0_io(J) )
            END IF
            
        END DO
            
    RETURN
    END SUBROUTINE
    
    

!**********************************************************************************************************************************
    SUBROUTINE WRT_FLOW_FA_Profile_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        CHARACTER(128) :: FLNM, FLNN
        CHARACTER(5)   :: STDIM(2)
        CHARACTER(4)   :: STFASTRESS(8)
        CHARACTER(2)   :: SJ2
        REAL(WP)   :: scaling1, scaling2, scaling3, scaling4,scaling5, scaling6, scaling7
        REAL(WP)   :: tke,Ruv_vis
        INTEGER(4)    :: I, J, N, L, TECFLG_FavAG(2), TECFLG_FavAG3, Nmax
        REAL(wp)   ::   COE, tem, Mtemp
        REAL(WP)   ::   TKE2, EPPSI, Fmu
        REAL(wp)   ::   Cmu, K2De
        
        REAL(WP)   :: FCT(0:NND2,NDV)
        REAL(WP)   :: COE1, COE2, DENtemp, Dintg, intgbdfc
        
!====================Favre Avderaged profiles===============================================================
        Nmax = 1
        STDIM(1) = 'undim'
        TECFLG_FavAG(1) = 101
        IF(ppdim==1) THEN
            Nmax = 2
            STDIM(2) = 'dimen'
            TECFLG_FavAG(2) = 102
        END IF
        
        
        !scaling1 = U0*U0*D0
        DO N=1, Nmax
            WRITE(PNTIM,'(1ES15.9)') phyTIME_io
            FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(N))//'.Profile.Flow.Favre.'//TRIM(PNTIM)//'.tec'
            OPEN (TECFLG_FavAG(N), FILE=TRIM(ADJUSTL(FLNM)))
            WRITE(TECFLG_FavAG(N),'(A)') 'TITLE = " Favre Averged Flow (35 variables)" '
            J=0;                          write(TECFLG_FavAG(N),'(A)',advance="no") 'VARIABLES = '
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Tauw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Dw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'DInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'D",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Muw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'MuInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'M",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'F_A",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'YCC",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ywd",' 
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ux",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Uy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Uz",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'P",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ruu",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Rvv",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Rww",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ruv",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ruw",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Rvw",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ruvvis",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dudy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dilation",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'densintg",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'bdfcintg",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'OmegRmsX",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'OmegRmsY",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'OmegRmsZ",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'RhoUU",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'RhoVV",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'RhoWW",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'MKE",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'TKE",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'UUUpp1_FA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'UUUpp2_FA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'UUUpp3_FA"'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'UUUpp112FA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'UUUpp332FA",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'skewness_U_FA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'skewness_V_FA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'skewness_W_FA"'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ANISTPinva2",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ANISTPinva3",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'LumleyX",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'LumleyY",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)'             ) '"'//SJ2//'drvFC" '
            
            WRITE(TECFLG_FavAG(N),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        END DO
        
        DO J=1,NCL2
        
            WRITE(TECFLG_FavAG(1),'(46ES20.12)') &
                TauwSD(J), &
                DensSD(J), DenAvew, D1xztL_F0_io(J), &
                ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                F_A, YCC(J), YWdiSD(J), &
                U_FA(J,1), U_FA(J,2), U_FA(J,3), U1xztL_F0_io(J,4),&
                uff2_FA(J,1,1), &
                uff2_FA(J,2,2), &
                uff2_FA(J,3,3), &
                uff2_FA(J,1,2), &
                uff2_FA(J,1,3), &
                uff2_FA(J,2,3), &
                Tau_Mean_RA(J,1,2), &
                dUdX_FA(J,1,2),&
                dUidXi(J),&
                densintg(J), &
                bdfcintg(J), &
                Omega_rms(J,1:3), &
                UGxztL_F0_io(J,1), &
                UGxztL_F0_io(J,4), &
                UGxztL_F0_io(J,6), &
                MKE_FA(J), &
                TKE_FA(J), &
                uff3_FA(J,1,1,1), &
                uff3_FA(J,2,2,2), &
                uff3_FA(J,3,3,3), &
                uff3_FA(J,1,1,2), &
                uff3_FA(J,3,3,2), &
                skewness_FA(J,1:3),&
                ANISTPinva_FA(J,2:3), &
                LumleyAxis_FA(J,1:2), &
                FUxztL_F0_io(J,4)

        END DO
        DO N=1,Nmax
            CLOSE(TECFLG_FavAG(N))
        END DO
        
        CALL WRT_FLOW_budgets_Profile_XZ_IO('Favre')
        
!======================================================================================================
        N=1
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        FLNM= TRIM(filepath4)//'Result.IO.VALIDATION.Profile.Flow.FANS.'//TRIM(PNTIM)//'.tec'
              
        OPEN (TECFLG_FavAG(N), FILE=TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_FavAG(N),'(A)') 'TITLE = " Favre Averged Flow (20 variables)" '
        J=0;                           write(TECFLG_FavAG(N),'(A)',advance="no") 'VARIABLES = '
        
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Tauw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Dw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'DInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'D",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Muw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'MuInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'M",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'F_A",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'YCC",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ywd",' 
        
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'X-VTau12-1",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'X-VTau12-2",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'X-TTauT12",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'X-bdfc",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'X-total",' 
        
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Y-pressure",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Y-VTau22-1",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Y-VTau22-2",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Y-TTauT22",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Y-bdfc",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Y-total",' 
        
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Z-VTau23-1",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Z-VTau23-2",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Z-TTauT23",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Z-bdfc",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)')              '"'//SJ2//'Z-total"' 
        
        WRITE(TECFLG_FavAG(N),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
       
        DO J=1, NCL2
            WRITE(TECFLG_FavAG(N),'(26ES20.12)')     &
            TauwSD(J), &
            DensSD(J), DenAvew, D1xztL_F0_io(J), &
            ViscSD(J), VisAvew, M1xztL_F0_io(J), &
            F_A, YCC(J), YWdiSD(J), &
            Tau_meaU_RA(J,1,2), Tau_Mean_RA(J,1,2)-Tau_meaU_RA(J,1,2), -uff2d_FA(J,1,2), F_A*IBuoF(1)*D1xztL_F0_io(J), &
            Tau_Mean_RA(J,1,2)-uff2d_FA(J,1,2)+F_A*IBuoF(1)*D1xztL_F0_io(J), &
            -U1xztL_F0_io(J,4), Tau_meaU_RA(J,2,2), Tau_Mean_RA(J,2,2)-Tau_meaU_RA(J,2,2), -uff2d_FA(J,2,2), &
            F_A*IBuoF(2)*D1xztL_F0_io(J), &
            -U1xztL_F0_io(J,4)+Tau_Mean_RA(J,2,2)-uff2d_FA(J,2,2)+F_A*IBuoF(2)*D1xztL_F0_io(J), &
            Tau_meaU_RA(J,2,3), &
            Tau_Mean_RA(J,2,3)-Tau_meaU_RA(J,2,3), &
           -uff2d_FA(J,2,3), &
            F_A*IBuoF(3)*D1xztL_F0_io(J), &
            Tau_Mean_RA(J,2,3)-uff2d_FA(J,2,3)+F_A*IBuoF(3)*D1xztL_F0_io(J)
        END DO
        
        !========================================================================================
    ! show RA FA relations
!        N=1
!        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
!        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(N))//'.Profile.Flow.FA_RA_Comp.'//TRIM(PNTIM)//'.tec'
!        OPEN (TECFLG_FavAG(N), FILE=TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(N),'(A)') 'TITLE = " Favre Averged Flow (29 variables)" '
!        J=0;                          write(TECFLG_FavAG(N),'(A)',advance="no") 'VARIABLES = '
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Utw",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Dw",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Muw",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Utave",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Dave",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'MuInt",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'YCC",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ywd",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'D",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'M",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'upp",' 
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'vpp",'
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'wpp",'
!        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)'             ) '"'//SJ2//'Mut"'
!        WRITE(TECFLG_FavAG(N),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        
!        DO J=1,NCL2
        
!            Ruv_vis=Tau_Mean_RA(J,1,2)!
!            tke   = 0.5_wp*(uff2d_FA(J,1,1)+uff2d_FA(J,2,2)+uff2d_FA(J,3,3))
!            EPPSI = 0.5_wp*(BUDG_viscs_dissp_duiuj(J,1)+BUDG_viscs_dissp_duiuj(J,4)+BUDG_viscs_dissp_duiuj(J,6))/D1xztL_F0_io(J)
!            !http://www.cfd-online.com/Wiki/Turbulence_dissipation_rate
            
!            WRITE(TECFLG_FavAG(N),'(29ES20.12)') &
!                UtauSD(J), DensSD(J), ViscSD(J), Utaw_ave_io, DenAvew, VisAvew, YCC(J), YWdiSD(J), &
!                D1xztL_F0_io(J), M1xztL_F0_io(J),&
!                U1xztL_F0_io(J,1)-U_FA(J,1), U1xztL_F0_io(J,2)-U_FA(J,2), U1xztL_F0_io(J,3)-U_FA(J,3)
!        END DO
!        CLOSE(TECFLG_FavAG(N))
        
        
        
        RETURN
    END SUBROUTINE
    
    
    !******************************************************************************************************
    SUBROUTINE WRT_FLOW_RA_Profile_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        CHARACTER(128) :: FLNM, FLNN
        CHARACTER(5)   :: STDIM(2)
        CHARACTER(4)   :: STFASTRESS(8)
        CHARACTER(2)   :: SJ2
        CHARACTER(3)   :: SJ3
        REAL(WP)       :: ux, uy, uz, p, &
                         tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, vorx, vory, vorz
        REAL(wp)    ::   COE, tem, Mtemp,Ruv_vis, dUdY
        REAL(WP)    :: scaling1 , scaling2
        INTEGER(4)  :: I, J, L
        INTEGER(4)  :: TECFLG_ReyAG(2), TECFLG_ReyAG3
        INTEGER(4)  :: N, Nmax, M
        REAL(WP)    :: TKE2, EPPSI, Fmu
        REAL(wp)    :: Cmu, K2De
        REAL(WP)    :: FCT(0:NCL2+1,NDV)
        REAL(WP)    :: COE1, COE2
        REAL(WP)    :: Yplus

        !==========================title=========================================
        Nmax = 1
        STDIM(1) = 'undim'
        TECFLG_ReyAG(1) = 101
        IF(ppdim==1) THEN
            Nmax = 2
            STDIM(2) = 'dimen'
            TECFLG_ReyAG(2) = 102
        END IF
        
        !scaling1 = U0*U0*D0
        
        DO N=1, Nmax
            WRITE(PNTIM,'(1ES15.9)') phyTIME_io
            FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(N))//'.Profile.Flow.Reynd.'//TRIM(PNTIM)//'.tec'
            OPEN (TECFLG_ReyAG(N), FILE=TRIM(ADJUSTL(FLNM)))
            WRITE(TECFLG_ReyAG(N),'(A)')     'TITLE = " Reynods Averaged Flow (41 variables)" '
            J=0;                          write(TECFLG_ReyAG(N),'(A)',advance="no") 'VARIABLES = '
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Tauw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Dw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'DInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'D",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Muw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'MuInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'M",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'F_A",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'YCC",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Ywd",' 
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Ux",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Uy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Uz",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'P",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Ruu",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Rvv",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Rww",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Ruv",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Ruw",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Rvw",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'Ruvvis",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'dudy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'dialation",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'densintg",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'bdfcintg",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'OmegRmsX",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'OmegRmsY",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'OmegRmsZ",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'RhoUU",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'RhoVV",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'RhoWW",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'MKE",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'TKE",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'UUUp1_RA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'UUUp2_RA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'UUUp3_RA"'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'UUUp112RA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'UUUp332RA"'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'skewness_U_RA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'skewness_V_RA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'skewness_W_RA"'
            
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'ANISTPinva2",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'ANISTPinva3",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'LumleyX",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)',advance="no") '"'//SJ2//'LumleyY",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(N),'(A)'             ) '"'//SJ2//'drvFC" '
            
            
            
            WRITE(TECFLG_ReyAG(N),'(A)')     'ZONE T=" '//TRIM(ADJUSTL(teczonename))//' " '
        END DO    
            
        !==========================write Flow info=========================================
        DO J=1,NCL2

            WRITE(TECFLG_ReyAG(1),'(46ES20.12)') &
                    TauwSD(J), &
                    DensSD(J), DenAvew, D1xztL_F0_io(J), &
                    ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                    F_A, YCC(J), YWdiSD(J), &
                    U1xztL_F0_io(J,1), &
                    U1xztL_F0_io(J,2), &
                    U1xztL_F0_io(J,3), &
                    U1xztL_F0_io(J,4), &
                    uf2_RA(J,1,1), &
                    uf2_RA(J,2,2), &
                    uf2_RA(J,3,3), &
                    uf2_RA(J,1,2), &
                    uf2_RA(J,1,3), &
                    uf2_RA(J,2,3), &
                    Tau_Mean_RA(J,1,2), &
                    DVDL1xztL_F0_io(J,1,2), &
                    dUidXi(J), &
                    densintg(J), &
                    bdfcintg(J), &
                    Omega_rms(J,1:3), &
                    UU_RA(J,1,1)*D1xztL_F0_io(J), &
                    UU_RA(J,2,2)*D1xztL_F0_io(J), &
                    UU_RA(J,3,3)*D1xztL_F0_io(J), &
                    MKE_RA(J), &
                    TKE_RA(J), &
                    uf3_RA(J,1,1,1), &
                    uf3_RA(J,2,2,2), &
                    uf3_RA(J,3,3,3), &
                    uf3_RA(J,1,1,2), &
                    uf3_RA(J,3,3,2), &
                    skewness_RA(J,1:3), &
                    ANISTPinva_RA(J,2:3), &
                    LumleyAxis_RA(J,1:2), &
                    FUxztL_F0_io(J,4)
            
        END DO
        DO N=1,Nmax
            CLOSE(TECFLG_ReyAG(N))
        END DO
            
        CALL WRT_FLOW_budgets_Profile_XZ_IO('Reynd')
        
!=======================quadrant terms=======================
        
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.Quad.Reynd.'//TRIM(PNTIM)//'.tec'
        
        OPEN (TECFLG_ReyAG(1), FILE=TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_ReyAG(1),'(A)') 'TITLE = " Reynolds Averged Flow (10+16*9 variables)" '
        J=0;                          write(TECFLG_ReyAG(1),'(A)',advance="no") 'VARIABLES = '
        
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Tauw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Dw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'DInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'D",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Muw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'MuInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'M",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'F_A",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'YCC",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Ywd",' 
        
        DO N=1, QUADHN
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADHV",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDR1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDR2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDR3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDR4",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADUV1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADUV2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADUV3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADUV4",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV11",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV12",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV13",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV14",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV21",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV22",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV23",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADDUV24",' 
        
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADVz1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADVz2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADVz3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADVz4",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADTK1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADTK2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADTK3",' 
            IF (N==QUADHN) THEN
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)'             ) '"'//SJ3//'QUADTK4",' 
            ELSE
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'QUADTK4",' 
            END IF
        END DO
        
        WRITE(TECFLG_ReyAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
    
        DO J=1,NCL2
            
            WRITE(TECFLG_ReyAG(1),'(10ES20.12)',advance="no") &
                TauwSD(J), &
                DensSD(J), DenAvew, D1xztL_F0_io(J), &
                ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                F_A, YCC(J), YWdiSD(J)
                
            DO N=1, QUADHN-1
                WRITE(TECFLG_ReyAG(1),'(25ES20.12)',advance="no") &
                    QUADHV(N), &
                    QUADDRxztL_F0_IO  (J,1:4,N),&
                    QUADUVxztL_F0_IO  (J,1:4,N),&
                    QUADDUV1xztL_F0_IO(J,1:4,N),&
                    QUADDUV2xztL_F0_IO(J,1:4,N),&
                    QUADVzxztL_F0_IO  (J,1:4,N),&
                    QUADTKxztL_F0_IO  (J,1:4,N)
            END DO
            N = QUADHN
            WRITE(TECFLG_ReyAG(1),'(25ES20.12)') &
                    QUADHV(N), &
                    QUADDRxztL_F0_IO  (J,1:4,N),&
                    QUADUVxztL_F0_IO  (J,1:4,N),&
                    QUADDUV1xztL_F0_IO(J,1:4,N),&
                    QUADDUV2xztL_F0_IO(J,1:4,N),&
                    QUADVzxztL_F0_IO  (J,1:4,N),&
                    QUADTKxztL_F0_IO  (J,1:4,N)
            

        END DO
        CLOSE(TECFLG_ReyAG(1))
        
        !=======================octant terms====\rho=================
        
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.OctD.Reynd.'//TRIM(PNTIM)//'.tec'
        
        OPEN (TECFLG_ReyAG(1), FILE=TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_ReyAG(1),'(A)') 'TITLE = " Reynolds Averged Flow (10+16*9 variables)" '
        J=0;                          write(TECFLG_ReyAG(1),'(A)',advance="no") 'VARIABLES = '
        
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Tauw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Dw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'DInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'D",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Muw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'MuInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'M",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'F_A",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'YCC",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Ywd",' 
        
        DO N=1, QUADHN
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDHV",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR7",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDR8",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV7",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDUV8",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV11",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV12",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV13",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV14",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV15",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV16",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV17",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV18",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV21",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV22",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV23",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV24",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV25",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV26",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV27",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDDUV28",' 
        
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz7",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDVz8",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK7",' 
            IF (N==QUADHN) THEN
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)'             ) '"'//SJ3//'OctDTK8",' 
            ELSE
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctDTK8",' 
            END IF
        END DO
        
        WRITE(TECFLG_ReyAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
    
        DO J=1,NCL2
            
            WRITE(TECFLG_ReyAG(1),'(10ES20.12)',advance="no") &
                TauwSD(J), &
                DensSD(J), DenAvew, D1xztL_F0_io(J), &
                ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                F_A, YCC(J), YWdiSD(J)
                
            DO N=1, QUADHN-1
                WRITE(TECFLG_ReyAG(1),'(49ES20.12)',advance="no") &
                    QuadHV(N), &
                    OctDDRxztL_F0_IO  (J,1:8,N),&
                    OctDUVxztL_F0_IO  (J,1:8,N),&
                    OctDDUV1xztL_F0_IO(J,1:8,N),&
                    OctDDUV2xztL_F0_IO(J,1:8,N),&
                    OctDVzxztL_F0_IO  (J,1:8,N),&
                    OctDTKxztL_F0_IO  (J,1:8,N)
            END DO
            N = QUADHN
            WRITE(TECFLG_ReyAG(1),'(49ES20.12)') &
                    QuadHV(N), &
                    OctDDRxztL_F0_IO  (J,1:8,N),&
                    OctDUVxztL_F0_IO  (J,1:8,N),&
                    OctDDUV1xztL_F0_IO(J,1:8,N),&
                    OctDDUV2xztL_F0_IO(J,1:8,N),&
                    OctDVzxztL_F0_IO  (J,1:8,N),&
                    OctDTKxztL_F0_IO  (J,1:8,N)
            

        END DO
        CLOSE(TECFLG_ReyAG(1))
        
        !=======================octant terms====T=================
        
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.OctT.Reynd.'//TRIM(PNTIM)//'.tec'
        
        OPEN (TECFLG_ReyAG(1), FILE=TRIM(ADJUSTL(FLNM)))
        WRITE(TECFLG_ReyAG(1),'(A)') 'TITLE = " Reynolds Averged Flow (10+16*9 variables)" '
        J=0;                          write(TECFLG_ReyAG(1),'(A)',advance="no") 'VARIABLES = '
        
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Tauw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Dw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'DInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'D",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Muw",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'MuInt",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'M",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'F_A",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'YCC",' 
        J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ2//'Ywd",' 
        
        DO N=1, QUADHN
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTHV",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR7",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDR8",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV7",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTUV8",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV11",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV12",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV13",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV14",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV15",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV16",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV17",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV18",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV21",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV22",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV23",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV24",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV25",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV26",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV27",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTDUV28",' 
        
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz7",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTVz8",' 
            
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK1",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK2",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK3",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK4",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK5",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK6",' 
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK7",' 
            IF (N==QUADHN) THEN
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)'             ) '"'//SJ3//'OctTTK8",' 
            ELSE
            J=J+1; WRITE(SJ3,'(1I3.3)') J; write(TECFLG_ReyAG(1),'(A)',advance="no") '"'//SJ3//'OctTTK8",' 
            END IF
        END DO
        
        WRITE(TECFLG_ReyAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
    
        DO J=1,NCL2
            
            WRITE(TECFLG_ReyAG(1),'(10ES20.12)',advance="no") &
                TauwSD(J), &
                DensSD(J), DenAvew, D1xztL_F0_io(J), &
                ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                F_A, YCC(J), YWdiSD(J)
                
            DO N=1, QUADHN-1
                WRITE(TECFLG_ReyAG(1),'(49ES20.12)',advance="no") &
                    QuadHV(N), &
                    OctTDRxztL_F0_IO  (J,1:8,N),&
                    OctTUVxztL_F0_IO  (J,1:8,N),&
                    OctTDUV1xztL_F0_IO(J,1:8,N),&
                    OctTDUV2xztL_F0_IO(J,1:8,N),&
                    OctTVzxztL_F0_IO  (J,1:8,N),&
                    OctTTKxztL_F0_IO  (J,1:8,N)
            END DO
            N = QUADHN
            WRITE(TECFLG_ReyAG(1),'(49ES20.12)') &
                    QuadHV(N), &
                    OctTDRxztL_F0_IO  (J,1:8,N),&
                    OctTUVxztL_F0_IO  (J,1:8,N),&
                    OctTDUV1xztL_F0_IO(J,1:8,N),&
                    OctTDUV2xztL_F0_IO(J,1:8,N),&
                    OctTVzxztL_F0_IO  (J,1:8,N),&
                    OctTTKxztL_F0_IO  (J,1:8,N)
            

        END DO
        CLOSE(TECFLG_ReyAG(1))
        
        RETURN
    END SUBROUTINE


    
!**********************************************************************************************************************************
    SUBROUTINE WRT_HEAT_FA_Profile_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        CHARACTER(128) :: FLNM, FLNN
        CHARACTER(5)   :: STDIM(2)
        CHARACTER(4)   :: STFASTRESS(3)
        CHARACTER(2)   :: SJ2
        REAL(WP)   :: ux, uy, uz, p, &
                      tke, Ruu, Rvv, Rww, Ruv, Ruw, Rvw, &
                      den, tem, enh_ra, enh_fa, drms, trms, &
                      hrms_ra, hrms_fa,thfx_ra, thfy_ra, thfz_ra, &
                      thfx_fa, thfy_fa, thfz_fa
        REAL(WP)   :: CpT(NCL2), Cp_eval
        INTEGER(4) :: I, J, TECFLG_FavAG(2), L, N, Nmax
        REAL(WP)   :: COE2, COE3
        REAL(WP)   :: HTEMP, Ptemp, Dtemp, Ttemp, Mtemp, Ktemp, Cptemp,  Btemp, Prtmp,qflux,du,dv,dw,dh
        REAL(WP)   :: scaling1, scaling2, scaling3, scaling4,scaling5, scaling6, scaling7
        REAL(WP)   :: du_per, dv_per, dw_per, dh_per
        REAL(WP)   :: MutOverPrt(3), Pruv(3)
        REAL(WP)   :: FCT(0:NCL2+1)

        Nmax = 1
        STDIM(1) = 'undim'
        TECFLG_FavAG(1) = 101
        IF(ppdim==1) THEN
            Nmax = 2
            STDIM(2) = 'dimen'
            TECFLG_FavAG(2) = 102
        END IF
!==================profiles================================================
        DO N=1, Nmax
            WRITE(PNTIM,'(1ES15.9)') phyTIME_io
            FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(N))//'.Profile.Heat.Transfer.'//TRIM(PNTIM)//'.tec'
                  
            OPEN (TECFLG_FavAG(N), FILE=TRIM(ADJUSTL(FLNM)))
            WRITE(TECFLG_FavAG(N),'(A)') 'TITLE = " Thermal variables (29 variables)" '
            J=0;                          write(TECFLG_FavAG(N),'(A)',advance="no") 'VARIABLES = '
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Tauw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Dw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'DInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'D",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Muw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'MuInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'M",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'F_A",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'YCC",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Ywd",' 
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Tw",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Hw",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Qw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Cpw",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'T",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'HRA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'HFA",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'D(T)",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'M(T)",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'K(T)",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Cp(T)",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Pr(T)",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Drms",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'Trms",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'HRArms",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'HFArms",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ThfRAx",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ThfRAy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ThfRAz",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ThfFAx",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ThfFAy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'ThfFAz",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'QfluxX",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'QfluxY",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'QfluxZ",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'du_per",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dv_per",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dw_per",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dh_per",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dTdy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dHRAdy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)',advance="no") '"'//SJ2//'dHFAdy",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(N),'(A)'             ) '"'//SJ2//'dDdy",'
            WRITE(TECFLG_FavAG(N),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        END DO
        
        DO J=1,NCL2
            
            drms    = DSQRT(dabs( D2xztL_F0_io(J) - D1xztL_F0_io(J) * D1xztL_F0_io(J) ) ) ! same
            trms    = DSQRT(dabs( T2xztL_F0_io(J) - T1xztL_F0_io(J) * T1xztL_F0_io(J) ) ) ! same
            hrms_ra = DSQRT(dabs( H2xztL_F0_io(J) - H1xztL_F0_io(J) * H1xztL_F0_io(J) ) )
            hrms_fa = DSQRT(dabs( H2xztL_F0_io(J) - DHxztL_F0_io(J)/D1xztL_F0_io(J)* &
                      (2.0_wp*H1xztL_F0_io(J) - DHxztL_F0_io(J)/D1xztL_F0_io(J) ) ) )
                      
            IF(GRAVDIR.eq.0) then
                du_per=0.0_wp
                dv_per=0.0_wp
                dw_per=0.0_wp
                dh_per=0.0_wp
            else
                du_per=( G1xztL_F0_io(J,1) - D1xztL_F0_io(J) * U1xztL_F0_io(J,1) )*F_A
                dv_per=( G1xztL_F0_io(J,2) - D1xztL_F0_io(J) * U1xztL_F0_io(J,2) )*F_A
                dw_per=( G1xztL_F0_io(J,3) - D1xztL_F0_io(J) * U1xztL_F0_io(J,3) )*F_A
                dh_per=( DHxztL_F0_io(J)   - D1xztL_F0_io(J) * H1xztL_F0_io(J)   )*F_A
            end if
            
            TTEMP = T1xztL_F0_io(J)
            IF(thermoStat==search_table) Then
                call NIST_SLEVAL_TD (Ttemp,Dtemp)
                call NIST_SLEVAL_TM (Ttemp,Mtemp)
                call NIST_SLEVAL_TK (Ttemp,Ktemp)
                call NIST_SLEVAL_TCP(Ttemp,CPtemp)
            ELSE IF(thermoStat==idealgas_law) Then
                Ptemp = U1xztL_F0_io(J,4)
                !call idealgas_TD (Ttemp,Dtemp,Ptemp)
                !call idealgas_TM (Ttemp,Mtemp)
                !call idealgas_TK (Ttemp,Ktemp)
                !call idealgas_TCP(Ttemp,CPtemp)
            ELSE
            END IF
            Prtmp=CPtemp*Mtemp/Ktemp       
            CpT(J) =     CPtemp 
            
            
            WRITE(TECFLG_FavAG(1),'(43ES20.12)') &
                    TauwSD(J), &
                    DensSD(J), DenAvew, D1xztL_F0_io(J), &
                    ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                    F_A, YCC(J), YWdiSD(J), &
                    TwSD(J), HwSD(J), &
                    QwSD(J), CpSD(J), &
                    T1xztL_F0_io(J), &
                    H1xztL_F0_io(J), &
                    H_FA(J), &
                    Dtemp, Mtemp, Ktemp, Cptemp, Prtmp, &
                    drms, trms, hrms_ra, hrms_fa, &
                    ufhfd_RA(J,1:3), &
                    uffhffd_FA(J,1:3),&
                    -DTDLKxztL_F0_io(J,1:3)*CTHECD, &
                    du_per,dv_per,dw_per,dh_per, &
                    dTdX(J,2), dHdX_RA(J,2), dHdX_FA(J,2), dDdX(J,2)
!            IF(ppdim==1) THEN        
!                WRITE(TECFLG_FavAG(2),'(29ES20.12)') YCC(J)*L0,(1.0_wp-dabs(YCC(J)))*REN*COE, &
!                        D1xztL_F0_io(J)*D0,T1xztL_F0_io(J)*T0, &
!                        H1xztL_F0_io(J)*CP0*T0+H0, &
!                        H_FA(J)*CP0*T0+H0, &
!                        M1xztL_F0_io(J)*M0, &
!                        drms*D0, trms*T0, &
!                        hrms_ra*CP0*T0+H0, hrms_fa*CP0*T0+H0, &
!                        Dtemp*D0, Mtemp*M0, Ktemp*K0, Cptemp*CP0, Prtmp, &
!                        ufhfd_RA(J,1:3)*U0*D0*CP0*T0+U0*D0*H0, uffhffd_FA(J,1:3)*U0*D0*CP0*T0+U0*D0*H0,&
!                        -DTDLKxztL_F0_io(J,1:3)*CTHECD*T0/L0*K0, &
!                        du_per*D0*U0,dv_per*D0*U0,dw_per*D0*U0,dh_per*D0*U0
!            END IF
        END DO
        DO N=1,Nmax
            CLOSE(TECFLG_FavAG(N))
        END DO
        
!=========================budgets=============================================================
        STFASTRESS(1) = 'THFx'
        STFASTRESS(2) = 'THFy'
        STFASTRESS(3) = 'THFz'
!        scaling1 = CP0*T0*D0*U0*U0/L0
!        scaling2 = H0*D0*U0*U0/L0
!        scaling3 = U0*K0*T0/L0/L0
!        scaling4 = CP0*T0*M0*U0/L0/L0
!        scaling5 = H0*M0*U0/L0/L0
!        scaling6 = CP0*T0*D0/L0
!        scaling7 = H0*D0/L0
        DO L=1, 3
            
            TECFLG_FavAG(1) = TECFLG_FavAG(1) + 1
            
            !DO N=1, Nmax
            WRITE(PNTIM,'(1ES15.9)') phyTIME_io
            FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.BUDGETS.Favre.'//TRIM(STFASTRESS(L))//  &
                  '.'//TRIM(PNTIM)//'.tec'
            OPEN (TECFLG_FavAG(1), FILE=TRIM(ADJUSTL(FLNM)))
            WRITE(TECFLG_FavAG(1),'(A)') 'TITLE = " Favre Averged Flow (26 variables)" '
            
            J=0;                          write(TECFLG_FavAG(1),'(A)',advance="no") 'VARIABLES = '
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Tauw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Dw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'DInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'D",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Muw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'MuInt",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'M",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'F_A",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'YCC",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Ywd",' 
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Tw",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Hw",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Qw",' 
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Cpw",'
            
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'prodc_stres",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'prodc_enthg",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'Turbu_diffu",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'dphdx_diffu",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'pdhdx_stran",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'ConHF_diffu",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'ConHF_dissp",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'viscs_diffu",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'viscs_dissp",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'press_accel",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'ConHF_accel",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'viscs_accel",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)',advance="no") '"'//SJ2//'prodc_dgbdf",'
            J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG_FavAG(1),'(A)'             ) '"'//SJ2//'balance"'
            WRITE(TECFLG_FavAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
            !END DO
        
            DO J=1,NCL2

                WRITE(TECFLG_FavAG(1),'(28ES20.12)') &
                        TauwSD(J), &
                        DensSD(J), DenAvew, D1xztL_F0_io(J), &
                        ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                        F_A, YCC(J), YWdiSD(J), &
                        TwSD(J), HwSD(J), &
                        QwSD(J), CpSD(J), &
                        BUDG_prodc_stres_thf(J,L), BUDG_prodc_enthg_thf(J,L), &
                        BUDG_Turbu_diffu_thf(J,L), BUDG_dphdx_diffu_thf(J,L), & 
                        BUDG_pdhdx_stran_thf(J,L), BUDG_ConHF_diffu_thf(J,L), &
                        BUDG_ConHF_dissp_thf(J,L), BUDG_viscs_diffu_thf(J,L), &
                        BUDG_viscs_dissp_thf(J,L), BUDG_press_accl1_thf(J,L), &
                        BUDG_ConHF_accel_thf(J,L), BUDG_viscs_accl1_thf(J,L), &
                        BUDG_prodc_gvfc2_thf(J),   BUDG_balance1_thf(J,L)
!                IF(ppdim==1 ) THEN      
!                    WRITE(TECFLG_FavAG(2),'(17ES20.12)') YCC(J)*L0, (1.0_wp-dabs(YCC(J)))*REN*COE, COE*U0, &
!                        BUDG_prodc_stres_thf(J,L)*scaling1+scaling2, BUDG_prodc_enthg_thf(J,L)*scaling1+scaling2, &
!                        BUDG_Turbu_diffu_thf(J,L)*scaling1+scaling2, BUDG_press_accl1_thf(J,L)*scaling1+scaling2, &
!                        BUDG_dphdx_diffu_thf(J,L)*scaling1+scaling2, BUDG_pdhdx_stran_thf(J,L)*scaling1+scaling2, &
!                        BUDG_ConHF_accel_thf(J,L)*scaling3, &
!                        BUDG_ConHF_diffu_thf(J,L)*scaling3, &
!                        BUDG_ConHF_dissp_thf(J,L)*scaling3, &
!                        BUDG_viscs_accl1_thf(J,L)*scaling4+scaling5, BUDG_viscs_diffu_thf(J,L)*scaling4+scaling5, &
!                        BUDG_viscs_dissp_thf(J,L)*scaling4+scaling5, &
!                        BUDG_prodc_gvfc2_thf(J)*scaling6+scaling7, BUDG_balance1_thf(J,L)*scaling4+scaling5
!                END IF
            END DO
            CLOSE(TECFLG_FavAG(1))
        END DO
        

        
        RETURN
    END SUBROUTINE
        
!**********************************************************************************************************************************
    SUBROUTINE WRT_HeatTransfer_Table_XZ_IO
    !=================================================================
        !======defination based on Joong Hun Bae, 2005===========
        !===Mdot = mass flow rate
        !===Mbuk = bulk mass flux
        !===Hdot = enthalpy flow rate
        !===Hbuk = bulk enthalpy
        !==
        !===Tbuk = bulk temperature
        !===Dbuk = bulk density
        !===Mbuk = bulk viscousity
        !===Kbuk = bulk thermal conductivity
        !===Cpbk = bulk Cp
        !===
        !===Ubuk = bulk velocity
        !===Rebk = bulk Reynolds no.
        !===Prbk = bulk Prandtl no.
        !===Nubk = bulk Nusselt number
        !===Bobk = Boyancy parameter
        !===Cfbk = local friction coefficient
        
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        CHARACTER(128) :: FLNM
        INTEGER(4)     :: I, J, IP, JP, JM, K, N
        REAL(WP)        :: Hw, Tw, Dw, Cpw, Mw, Kw
        REAL(WP)        :: HTEMP, Dtemp, Ttemp, Mtemp, Ktemp, Cptemp, Btemp, Ptemp, Ttemp1, Ttemp2, Ktemp1, Ktemp2, Ltemp(2)
        REAL(WP)        :: Buoya_1, Buoya_2, hc_Dwve(2),NubkTsdave(2), Nubk_dTw_ave
        INTEGER(4) :: TECFLG = 200
        REAL(WP)      :: DuDyL, DuDyU
        REAL(WP)      :: rtmpmin, TempDiff
        REAL(WP)      :: CMWORK, CMWORK_d ! work done by effective pressure gradient
        REAL(WP)      :: dqw_d, dqw
        REAL(WP)      :: Buoyasd_1(2), Buoyasd_2(2)
        
        
        !===================================================================
        !++++++++++++++++sided bulk values++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL PP_J4TbukTsd
        !WRITE(*,*) '# L4TbkSsd' ,L4TbkSsd(1:2)
        !WRITE(*,*) '# Ldist_io',Ldist_io(1:2)
        !WRITE(*,*) '# DbukSsd  ',DbukSsd(1:2)
        !WRITE(*,*) '# MbukSsd  ',MbukSsd(1:2)
        
        !=========================================================================
        TECFLG = 200
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        OPEN(TECFLG,FILE=TRIM(filepath4)//'Result.IO.Table.WallandBulk.'//TRIM(PNTIM)//'.tec')!, POSITION='APPEND')
        
        WRITE(TECFLG, '(A)         ') '%######## Reference States (dimensionaless)####################'   
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF Re                 =', REN
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF PRT0               =', PRT0
        WRITE(TECFLG, '(A)') '   '
        WRITE(TECFLG, '(A)         ') '%######## Reference States (dimensional)  #####################'   
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF L0(m)              =', L0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF U0(m/s)            =', U0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF Mdot(Kg/m2s)       =', U0*D0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF P0(Pa)             =', P0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF T0(K)              =', T0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF T0(C)              =', T0-TK2C
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF H0(J)              =', H0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF D0(Kg/m3)          =', D0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF K0(W/m-K)          =', K0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF M0(Pa-s)           =', M0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF Cp0(J/Kg/K)        =', CP0
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF B0(1/K)            =', B0
        WRITE(TECFLG, '(A, 1ES17.9)') 'gR/U0^2                =', F_A
        WRITE(TECFLG, '(A)') '   '
        
        !==========wall parameters==================================================
        HWAL_RA_d(:)  = HWAL_RA(:)*T0*CP0+H0
        Dwal_d(:)  = Dwal(:)*D0
        Twal_d(:)  = Twal(:)*T0
        Mwal_d(:)  = Mwal(:)*M0
        Kwal_d(:)  = Kwal(:)*K0
        Cpwal_d(:) = CPwal(:)*CP0
        
        Ttau(1:2)         = dabs(qw(1:2)/(Dwal(1:2)*Cpwal(1:2)*Utaw_io(1:2)))!Ttau_d(1:2)/T0
        Ttau_d(1:2)       = dabs(qw_d(1:2)/(Dwal_d(1:2)*Cpwal_d(1:2)*Utaw_io(1:2)*U0))
        
        WRITE(TECFLG, '(A)         ') '%########Variables (Wall)######################################'
        WRITE(TECFLG, '(A,24X,2(A,6X),2(A,4X))') '%#','dimensional','dimensional','dimensionless','dimensionless'
        WRITE(TECFLG, '(A,24X,4(A,6X))') '%#','wall(y=-1) ','wall(y=+1) ','wall(y=-1) ','wall(y=+1) '
        WRITE(TECFLG, '(A, 4ES17.9)') 'Tw: Temperature (K)     =', Twal_d(1:2), Twal(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)') 'Tw: Temperature (C)     =', Twal_d(1:2)-TK2C
        WRITE(TECFLG, '(A, 4ES17.9)') 'Tt: Friction Ttau (K)   =', Ttau_d(1:2), Ttau(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') 'Dw: Density (Kg/m3)     =', Dwal_d(1:2), Dwal(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') 'Hw: Enthalpy(J/Kg)      =', Hwal_RA_d(1:2), Hwal_RA(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') 'Mw: Viscousity(Pa-s)    =', Mwal_d(1:2), Mwal(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') 'Kw: Conductivity(W/m-K) =', Kwal_d(1:2), Kwal(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') 'Cpw: Specific heat(Cp)  =', Cpwal_d(1:2), Cpwal(1:2)
        WRITE(TECFLG, '(A)') '   '
        
        !============find the location of Tcp=========================
        rtmpmin = 1.0E+10_WP
        IF(DMAX1(Twal(1),Twal(2)).GT.T4CpMax) THEN
            DO J=1, NCL2
                TempDiff=dabs(T1xztL_F0_io(J)-T4CpMax)
                IF(TempDiff.LT.rtmpmin) THEN
                    rtmpmin = TempDiff
                    J4Tpc   = J
                END IF
            END DO
        ELSE
            J4Tpc = 0
        END IF
        
        !===============bulk values=======================================================================
        !==================calculate bulk mass flux, and bulk enthalpy=========
        Mdot = 0.0_WP
        Hdot = 0.0_WP
        DO J=1,NCL2
            DO K=1,NCL3
                Mdot   = Mdot + G1xztL_F0_io(J,1)/RCCI1(J)/DYFI(J)/DZI
                Hdot   = Hdot + GHxztL_F0_io(J,1)/RCCI1(J)/DYFI(J)/DZI
            ENDDO
        ENDDO

        Gbuk =  Mdot /AREA_INLET
        Hbuk =  Hdot /Mdot
        
        HTEMP = Hbuk
                
        IF(thermoStat==search_table) Then
            call NIST_SLEVAL_HT(Htemp,Ttemp)
            call NIST_SLEVAL_TD(Ttemp,Dtemp)
        
            call NIST_SLEVAL_TM(Ttemp,Mtemp)
            call NIST_SLEVAL_TK(Ttemp,Ktemp)
            call NIST_SLEVAL_TCP(Ttemp,CPtemp)
            call NIST_SLEVAL_TB(Ttemp,Btemp)
        ELSE IF(thermoStat==idealgas_law) Then
            Ptemp = 1.0_wp
            !call idealgas_HT(Htemp,Ttemp)
            !call idealgas_TD(Ttemp,Dtemp,Ptemp)
        
            !call idealgas_TM(Ttemp,Mtemp)
            !call idealgas_TK(Ttemp,Ktemp)
            !call idealgas_TCP(Ttemp,CPtemp)
            !call idealgas_TB(Ttemp,Btemp)
        ELSE
        END IF
        
        Dbuk  = Dtemp
        Tbuk  = Ttemp
        Mbuk  = Mtemp
        Kbuk  = Ktemp
        Cpbk  = Cptemp
        Bbuk  = Btemp
        Ubuk  = Gbuk/Dbuk
        
        Gbuk_d = Gbuk * D0 * U0
        Hbuk_d = Hbuk *T0*CP0+H0
        Ubuk_d= Ubuk*U0
        
        Dbuk_d  = Dbuk*D0
        Tbuk_d  = Tbuk*T0
        Mbuk_d  = Mbuk*M0
        Kbuk_d  = Kbuk*K0
        Cpbk_d  = CPbk*CP0
        Bbuk_d  = Bbuk*B0
        
        WRITE(TECFLG, '(A)         ') '%########Case Set up##########################################'
        IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
        WRITE(TECFLG, '(A, A)      ') 'Wall thermal conditions         =', ' Constant Wall temperature(K), (C)'
        WRITE(TECFLG, '(A, 2ES17.9)') '     Wall temperature (K) & (C) at (y=-1) =', Twal_d(1), Twal_d(1)-TK2C
        WRITE(TECFLG, '(A, 2ES17.9)') '     Wall temperature (K) & (C) at (y=+1) =', Twal_d(2), Twal_d(2)-TK2C
        WRITE(TECFLG, '(A, 1ES17.9)') 'Constant Mass Flux (Kg/m2-s)        =', Gbuk_d
        WRITE(TECFLG, '(A)') '   '
        END IF
        
        WRITE(TECFLG, '(A)         ') '%########Variables (Bulk)#####################################'
        WRITE(TECFLG, '(A,24X,4(A,6X))') '%#','dimensional','dimensionless'
        WRITE(TECFLG, '(A, 2ES17.9)') 'Mdot:Mass flux(Kg/m2-s)      =', Gbuk_d, Gbuk
        WRITE(TECFLG, '(A, 2ES17.9)') 'Ub  :Velocity (m/s)          =', Ubuk_d, Ubuk
        WRITE(TECFLG, '(A, 2ES17.9)') 'Hb  :Enthalpy(J/Kg)          =', Hbuk_d, Hbuk
        WRITE(TECFLG, '(A, 2ES17.9)') 'Tb: Temperature (K)          =', Tbuk_d, Tbuk
        WRITE(TECFLG, '(A, 1ES17.9)') 'Tb: Temperature (C)          =', Tbuk_d-TK2C
        WRITE(TECFLG, '(A, 2ES17.9)') 'Db: Density (Kg/m3)          =', Dbuk_d, Dbuk
        WRITE(TECFLG, '(A, 2ES17.9)') 'Mub:Viscousity(Pa-s)         =', Mbuk_d, Mbuk
        WRITE(TECFLG, '(A, 2ES17.9)') 'Kb:Conductivity(W/m-K)       =', Kbuk_d, Kbuk
        WRITE(TECFLG, '(A, 2ES17.9)') 'Cpb:Specific Heat(Cp)        =', Cpbk_d, Cpbk
        WRITE(TECFLG, '(A)') '   '
        
        
        !===============search===table==========
        !================================
        Rebk = Dbuk*Ubuk*2.0_WP*ALX2/Mbuk * REN ! based on diameter or the whole height of channel.
        Prbk = Mbuk*Cpbk/Kbuk      * Prt0
        
        WRITE(TECFLG, '(A)         ') '%########Variables (Bulk), SSzero y=-1 side, y=1 side undimensional######################'
        WRITE(TECFLG, '(A, 3ES17.9)') 'Reynolds No                  =', Rebk, RebkSsd(1:2)
        WRITE(TECFLG, '(A, 3ES17.9)') 'Prandlt No                   =', Prbk, PrbkSsd(1:2)
        !====================================================================

        IF( BCWALLHEAT(itopwall)==isoFluxWall .and. BCWALLHEAT(ibotwall)==isoFluxWall) THEN
            qw_d(1:2)   = WHEAT0_DIM(1:2)
            !====convective heat transfer rate====hc=qw/(Tw-Tb)====
            hc_d(1:2) = qw_d(1:2)/(Twal_d(1:2)-Tbuk_d)
            hc_Dwve(1:2) = 0.5_wp*(dabs(qw_d(1))+dabs(qw_d(2)))/(Twal_d(1:2)-Tbuk_d)
            !====bulk Grashof number===========
            
            Grbk(1:2)   = G_A*Bbuk_d*qw_d(1:2)*(L0**4)*Dbuk_d*Dbuk_d/Mbuk_d/Mbuk_d/Kbuk_d
            !WRITE(*,*) '# Grbk', Grbk(1:2)
        END IF
        IF( BCWALLHEAT(itopwall)==isoThermalWall .and. BCWALLHEAT(ibotwall)==isoThermalWall) THEN
            !====wall heat flux========== here is not corret, as bar over (kdT/dy) /= kdbarT/dy
            qw_d(1)= qw(1)*D0*U0*CP0*T0 !Kwal_d(1)*(Twal_d(1)-T1xztL_F0_io(1)*T0   )/((YCC(1)-YND(1))*L0)
            qw_d(2)= qw(2)*D0*U0*CP0*T0 !Kwal_d(2)*(Twal_d(2)-T1xztL_F0_io(NCL2)*T0)/((YND(NND2)-YCC(NCL2))*L0)
            
            !test
            !!WRITE(*,*) 'qwc', qw_d(1), Kwal_d(1)*(Twal_d(1)-T1xztL_F0_io(1)*T0   )/((YCC(1)-YND(1))*L0)       !Same
            !!WRITE(*,*) 'qwh', qw_d(2), Kwal_d(2)*(Twal_d(2)-T1xztL_F0_io(NCL2)*T0)/((YND(NND2)-YCC(NCL2))*L0) !Same
            
            qw_d_ave = 0.5_wp*(dabs(qw_d(1))+dabs(qw_d(2)))
            !====convective heat transfer rate====hc=qw/(Tw-Tb)====
            hc_dTw_d(1:2)= qw_d(1:2)/dabs(Twal_d(1)-Twal_d(2))
            hc_dTw_d_ave = qw_d_ave /dabs(Twal_d(1)-Twal_d(2))
            
            hc_d(1:2)    = qw_d(1:2)/DABS(Twal_d(1:2)-Tbuk_d)
            hc_Dwve(1:2) = dabs(qw_d_ave /DABS(Twal_d(1:2)-Tbuk_d))
            
            
            
            !====bulk Grashof number===========
            Grbk(1) = G_A*Bbuk_d*dabs(Twal_d(1)-Twal_d(2))*((2.0_wp*L0)**3)*Dbuk_d*Dbuk_d/Mbuk_d/Mbuk_d
            !WRITE(*,*) '# Grbk', Grbk(1)
        END IF
        
        !===============Nasselt number==============
        Nubk_dTw(1:2)     = 2.0_wp*L0*hc_dTw_d(1:2)/Kbuk_d  ! based on q_w(h,c)/(Tc-Th)
        Nubk_dTw_ave      = 2.0_wp*L0*hc_dTw_d_ave/Kbuk_d   ! based on averaged q_w(h+c)/(Tc-Th)
        
        NubkTsd(1:2)         = L4Tbk(1:2)*hc_d(1:2)/KbukTsd(1:2)/K0    ! based on q_w(h,c)/(Tw-Tb)
        NubkTsdave(1:2)      = L4Tbk(1:2)*hc_Dwve(1:2)/KbukTsd(1:2)/K0  ! based on averaged q_w(h+c)/(Tw-Tb)
        
        
        
        !============work done by shear stress=======================
        !CMWORK_d = dabs( (2.0_wp*Tauw_d_ave_io)/(DenAvew*D0)*(U0*D0) )
        CMWORK_d=Tauw_d_ave_io*Ubuk_d
        dqw_d    = dabs( dabs(qw_d(1))-dabs(qw_d(2)))
        
        !CMWORK   = dabs( (2.0_wp*Tauw_ave_io)/DenAvew * Gbuk)
        CMWORK   =Tauw_ave_io*Ubuk
        dqw      = dabs( dabs(qw(1))-dabs(qw(2)))*(T0*Cp0/U0/U0)
        
        WRITE(TECFLG, '(A)         ') '%########Energy Conservation#############################'
        WRITE(TECFLG, '(A,24X,4(A,6X))') '%#','dimensionless','dimensional'
        WRITE(TECFLG, '(A, 2ES19.9)') 'Net Heat Flux =           ', dqw,    dqw_d
        WRITE(TECFLG, '(A, 2ES19.9)') 'Work Done by (dp/dx)_eff =', CMWORK, CMWORK_d
        WRITE(TECFLG, '(A, 2ES19.9)') 'Above two diff =          ', dqw-CMWORK, dqw_d-CMWORK_d
        WRITE(TECFLG,'(A)') '  '
        
        WRITE(TECFLG, '(A)         ') '%########Heat transfer parameters#############################'
        WRITE(TECFLG,'(A,30X,2(A,7X),A,1X,A)') '%#', 'wall(y=-1)','wall(y=+1)','Two-wall average','Error(\%)=(W1-W2)/(W1+W2)*100'
        WRITE(TECFLG, '(A, 4ES19.9)') 'Qw: Wall heat flux(undim)    =', qw(1:2),0.5_wp*(qw(1)+qw(2)), &
                                                                        dabs(dabs(qw(1))-dabs(qw(2)))/(qw(1)+qw(2))*100.0_wp
        WRITE(TECFLG, '(A, 4ES19.9)') 'Qw: Wall heat flux(W/m2)     =', qw_d(1:2),    qw_d_ave, &
                                                                        0.5_wp*dabs(dabs(qw_d(1))-dabs(qw_d(2)))/qw_d_ave*100.0_wp
        WRITE(TECFLG,'(A)') '  '                                                                
        
        !===========integrated density=and Nussult number=========
        D_int= 0.0_wp
        Nu_int= 0.0_WP
        DO J=1,NND2
            IF (J==1) THEN
                TTEMP = T1xztL_F0_io(J)
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_TD(Ttemp,Dtemp)
                    call NIST_SLEVAL_TK(Ttemp,Ktemp)
                else if (thermoStat==idealgas_law) then
                    Ptemp = U1xztL_F0_io(J,4)
                    !call idealgas_TD(Ttemp,Dtemp,Ptemp)
                    !call idealgas_TK(Ttemp,Ktemp,Ptemp)
                else
                end if
                D_int  = D_int  + ( T1xztL_F0_io(J)-Twal(1) )*Dtemp
                Nu_int = Nu_int + YCC(J)/(Ktemp+Kwal(1))*(YCC(J)-0.0_WP)
            ELSE IF(J==NND2) THEN
                Dtemp = Dwal(2)
                TTEMP = T1xztL_F0_io(J-1)
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_TK(Ttemp,Ktemp)
                else if (thermoStat==idealgas_law) then
                    Ptemp = U1xztL_F0_io(J,4)
                    !call idealgas_TK(Ttemp,Ktemp,Ptemp)
                else
                end if
                D_int = D_int + ( Twal(2)-T1xztL_F0_io(J-1) )*Dtemp
                Nu_int = Nu_int + (YCC(J-1)+YND(J))/(Ktemp+Kwal(2))*(YND(J)-YCC(J-1))
            ELSE
                TTEMP1 = T1xztL_F0_io(J-1)
                TTEMP2 = T1xztL_F0_io(J)
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_TD(Ttemp2,Dtemp)
                    call NIST_SLEVAL_TK(Ttemp1,Ktemp1)
                    call NIST_SLEVAL_TK(Ttemp2,Ktemp2)
                else if (thermoStat==idealgas_law) then
                    Ptemp = U1xztL_F0_io(J,4)
                    !call idealgas_TD(Ttemp2,Dtemp,Ptemp)
                    !call idealgas_TK(Ttemp1,Ktemp1,Ptemp)
                    !call idealgas_TK(Ttemp2,Ktemp2,Ptemp)
                else
                end if
                D_int  = D_int + ( T1xztL_F0_io(J)-T1xztL_F0_io(J-1) )*Dtemp
                Nu_int = Nu_int + YND(J)/((YCL2ND_WFF(J)*Ktemp2+YCL2ND_WFB(J)*Ktemp1))*((YCC(J)-YCC(J-1)))
            END IF
        END DO
        Nu_int = Nu_int/(2.0_WP)*L0/K0*qw_d_ave/((Twal(2)-Twal(1))*T0)
        !WRITE(*,*) '# Nu_int', Nu_int
        
        IF(dabs(Twal(2)-Twal(1)).LT.1.0E-12_wp) THEN
            D_int = Dtemp
        ELSE
            D_int = D_int/(Twal(2)-Twal(1))
        END IF
        !WRITE(*,*) '# D_int', D_int
        D_int_d = D_int * D0
        
        Grbk_drho = (Dbuk_d-D_int_d)/Dbuk_d * G_A*(L0**3)*Dbuk_d*Dbuk_d/Mbuk_d/Mbuk_d
        !WRITE(*,*) '# Grbk_drho', Grbk_drho
        
        !====================================================================
        Buoya_1 = Grbk(1) / (Rebk**3.425_wp) * (Prbk**0.8_wp)
        Buoya_2 = Grbk(1) / (Rebk**2.7_wp)
        !WRITE(*,*) '# Buoya_12', Buoya_1, Buoya_2
        
        
        
        Dsd_int= 0.0_wp
        DO J=1,J4TbukSsd(1)
            IF (J==1) THEN
                TTEMP = T1xztL_F0_io(J)
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_TD(Ttemp,Dtemp)
                else if (thermoStat==idealgas_law) then
                    Ptemp = U1xztL_F0_io(J,4)
                    !call idealgas_TD(Ttemp,Dtemp,Ptemp)
                else
                end if
                Dsd_int(1)  = Dsd_int(1)  + ( T1xztL_F0_io(J)-Twal(1) )*Dtemp
            ELSE IF(J==NND2) THEN
                Dtemp = Dwal(2)
                Dsd_int(1) = Dsd_int(1) + ( Twal(2)-T1xztL_F0_io(J-1) )*Dtemp
            ELSE
                TTEMP2 = T1xztL_F0_io(J)
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_TD(Ttemp2,Dtemp)
                else if (thermoStat==idealgas_law) then
                    Ptemp = U1xztL_F0_io(J,4)
                    !call idealgas_TD(Ttemp2,Dtemp,Ptemp)
                else
                end if
                Dsd_int(1)  = Dsd_int(1) + ( T1xztL_F0_io(J)-T1xztL_F0_io(J-1) )*Dtemp
            END IF
        END DO
        
        IF(dabs(Twal(2)-Twal(1)).LT.1.0E-12_wp) THEN
            Dsd_int(1) = Dtemp
        ELSE
            Dsd_int(1) = Dsd_int(1)/(T1xztL_F0_io(J4TbukSsd(1))-Twal(1))
        END IF
        
        
        DO J=J4TbukSsd(2), NND2
            IF (J==1) THEN
                TTEMP = T1xztL_F0_io(J)
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_TD(Ttemp,Dtemp)
                else if (thermoStat==idealgas_law) then
                    Ptemp = U1xztL_F0_io(J,4)
                    !call idealgas_TD(Ttemp,Dtemp,Ptemp)
                else
                end if
                Dsd_int(2)  = Dsd_int(2)  + ( T1xztL_F0_io(J)-Twal(1) )*Dtemp
            ELSE IF(J==NND2) THEN
                Dtemp = Dwal(2)
                Dsd_int(2) = Dsd_int(2) + ( Twal(2)-T1xztL_F0_io(J-1) )*Dtemp
            ELSE
                TTEMP2 = T1xztL_F0_io(J)
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_TD(Ttemp2,Dtemp)
                else if (thermoStat==idealgas_law) then
                    Ptemp = U1xztL_F0_io(J,4)
                    !call idealgas_TD(Ttemp2,Dtemp,Ptemp)
                else
                end if
                Dsd_int(2)  = Dsd_int(2) + ( T1xztL_F0_io(J)-T1xztL_F0_io(J-1) )*Dtemp
            END IF
        END DO
        
        IF(dabs(Twal(2)-Twal(1)).LT.1.0E-12_wp) THEN
            Dsd_int(2) = Dtemp
        ELSE
            Dsd_int(2)= Dsd_int(2)/(Twal(2)-T1xztL_F0_io(J4TbukSsd(2)))
        END IF

        !WRITE(*,*) '# Dsd_int', Dsd_int(1:2)
    
        !Ltemp(1:2) = Ldist_io(1:2)
        Ltemp(1:2) = Ldist_io(1:2)
        
        GrbkSsd(1:2) = G_A*((Ltemp(1:2)*L0)**3)*(DbukSsd(1:2)*D0*DbukSsd(1:2)*D0)/(MbukSsd(1:2)*M0*MbukSsd(1:2)*M0)* &
                      (DbukSsd(1:2)-Dsd_int(1:2))/DbukSsd(1:2)
        !!WRITE(*,*) '# GrbkSsd', GrbkSsd(1:2)
        
        hcsd_d(1:2) = qw_d(1:2)/(Twal_d(1:2)-TbukSsd(1:2)*T0)        
        NubkSsd(1:2) = Ltemp(1:2)*L0*hcsd_d(1:2)/(KbukSsd(1:2)*K0)    ! based on q_w(h,c)/(Tw-Tb)
        !!WRITE(*,*) '# NubkSsd', NubkSsd(1:2)
        
        Buoyasd_1(1:2) = GrbkSsd(1:2) / ((RebkSsd(1:2)/Ldist_io(1:2)*Ltemp(1:2))**3.425_wp) * (PrbkSsd(1:2)**0.8_wp)
        Buoyasd_2(1:2) = GrbkSsd(1:2) / ((RebkSsd(1:2)/Ldist_io(1:2)*Ltemp(1:2))**2.7_wp)
        !WRITE(*,*) '# Buoyasd_1', Buoyasd_1(1:2)
        !WRITE(*,*) '# Buoyasd_2', Buoyasd_2(1:2)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        WRITE(TECFLG, '(A)         ') '%########integral, dimensionless(c,h), dimensinoal(c,h) ######################'
        WRITE(TECFLG, '(A, 2ES17.9)') 'Nu:      Integral                =', Nu_int
        WRITE(TECFLG,'(A)') '  '
        WRITE(TECFLG, '(A)         ') '%########sided bulk average(Zero shear stress balance)#######'
        WRITE(TECFLG, '(A,2X,I3.1,2ES17.9)')'Jindex,Ycc,  =', J4SS0, YCC(J4SS0)
        WRITE(TECFLG,'(A,15X,2(A,7X),A)')   '%#', 'side(y=-1)','side(y=+1)'
        WRITE(TECFLG, '(A, 3ES17.9)')       'Y+max        =', (YCC(J4SS0)+1.0_WP)*Ret_io(1)/Ldist_io(1), &
                                                              (1.0_WP-YCC(J4SS0))*Ret_io(2)/Ldist_io(2)
        WRITE(TECFLG, '(A, 2ES17.9)')       'Reynolds tau  each side:', Ret_io(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)')       'Reynolds bulk each side:', RebkSsd(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)')       'Velocity bulk each side:', UbukSsd(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)')       'Charac L bulk each side:', Ldist_io(1:2)
        WRITE(TECFLG,'(A)') '  '
        
        WRITE(TECFLG, '(A)         ') '%########dimensional (1)on the cold wall, (2)on the hot wall, (3)averaged (cold, hot)###'
        WRITE(TECFLG, '(A, 3ES17.9)') 'Qw: Wall heat flux(W/m2)           =', qw_d(1:2),    qw_d_ave
        WRITE(TECFLG, '(A         )') '1:==based on whole domain L===     '
        WRITE(TECFLG, '(A, 2ES17.9)') '1:Dint:    Integral density(Kg/m3) =', D_int, D_int_d
        WRITE(TECFLG, '(A, 3ES17.9)') '1:hc=qw/(Th-Tc) (W/m2/K)           =', hc_dTw_d(1:2), hc_dTw_d_ave
        WRITE(TECFLG, '(A, 3ES17.9)') '1:Nu=hc*2R/Kb                      =', Nubk_dTw(1:2), Nubk_dTw_ave
        WRITE(TECFLG, '(A, 1ES17.9)') '1:Grb_drho                         =', Grbk_drho
        WRITE(TECFLG, '(A, 1ES17.9)') '1:Grb/Re2.7                        =', Buoya_2
        WRITE(TECFLG, '(A, 1ES17.9)') '1:Grb/Re3.4/Pr0.8                  =', Buoya_1
        WRITE(TECFLG,'(A)') '  '
        WRITE(TECFLG, '(A         )') '2:==based on Tbulk  sided=====    '
        WRITE(TECFLG, '(A, I3.1,1ES17.9)')'2:Jindex,Ycc                       =', J4Tbk, YCC(J4Tbk)
        WRITE(TECFLG, '(A, 2ES17.9)') '2:Ldist:                           =', L4TbkTsd(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') '2:hc=qw/(T_{c/h}-Tb) (W/m2/K)      =', hc_d(1:2),    hc_Dwve(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') '2:Nu=hc*Lt/Kbt(Tbulk sep)          =', NubkTsd(1:2), NubkTsdave(1:2)
        WRITE(TECFLG,'(A)') '  '
        WRITE(TECFLG, '(A         )') '3:==based on SSzero sided=====    '
        WRITE(TECFLG, '(A, I3.1,1ES17.9)')'3:J4SS0,Ycc                        =', J4SS0, YCC(J4SS0)
        WRITE(TECFLG, '(A, 2ES17.9)') '3:Ldist:                           =', Ldist_io(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)') '3:Ldist to Tbulk each side         =', L4TbkSsd(1:2)
        WRITE(TECFLG, '(A, 2I3.1,2ES17.9)')'3:J4TbukSsd,Ycc                    =', J4TbukSsd(1:2), YCC(J4TbukSsd(1:2))
        WRITE(TECFLG, '(A, 2ES17.9)') '3:TbukSsd:                         =', TbukSsd(1:2)
        WRITE(TECFLG, '(A, 4ES17.9)') '3:DintSsd: Integral density(Kg/m3) =', Dsd_int(1:2), Dsd_int(1:2)*D0
        WRITE(TECFLG, '(A, 2ES17.9)') '3:hc=qw/(T_{wc/h}-T{bc/h)  (W/m2/K)=', hcsd_d(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)') '3:Nu=hc*Lt/Kbt(SS0 sep)            =', NubkSsd(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)') '3:Grb_drho(lc,rhobc)               =', GrbkSsd(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)') '3:Grb/Re2.7                        =', Buoyasd_2(1:2)
        WRITE(TECFLG, '(A, 2ES17.9)') '3:Grb/Re3.4/Pr0.8                  =', Buoyasd_1(1:2)
        
        
        WRITE(TECFLG,'(A)') '  '
        
        
        
        RETURN
    END SUBROUTINE
    
!************************************************************************************************************** 
    SUBROUTINE WRT_Cf_Table_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: TECFLG  = 200
        INTEGER(4) :: TECFLG2 = 202
        INTEGER(4)    :: I, J
        REAL(WP)      :: DuDy(2)
        REAL(WP)      :: dxplus(2), dzplus(2), yMINplus(2), yMAXplus(2), dyMIN, dyMAX
        REAL(WP)      :: DENtemp, ddenintg,denmintg
        
        REAL(WP), ALLOCATABLE :: dxplusj(:,:), dzplusj(:,:), dyplusj(:,:)
        REAL(WP)      :: InnerL(2)
        REAL(WP)      :: Tau_diff, Tau_avag, Llocal(2)
         
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        OPEN(TECFLG,FILE=TRIM(filepath4)//'Result.IO.Table.WallandBulk.'//TRIM(PNTIM)//'.tec', POSITION='APPEND')
        
        Tauw_io(1:2) = DABS(Tauw_io(1:2))
        Cf0_io(1:2)  = 2.0_WP*dabs(Tauw_io(1:2))
        IF(thermlflg ==0) THEN
            Cfbk_io(1:2)   = Cf0_io(1:2)
        END IF
        IF(thermlflg ==1) THEN
            Cfbk_io(1:2)   = 2.0_WP*Tauw_io(1:2)/Dbuk/Ubuk/Ubuk
        END IF
        
        Cf0_ave_io = 0.5*(Cf0_io(1)     + Cf0_io(2) )
        
        InnerL(1:2)   = 1.0_WP/Ret_io(1:2)
        
        WRITE(TECFLG,'(A)') '  '
        WRITE(TECFLG, '(A)         ') '%########Wall Coefficient ####################################'
        WRITE(TECFLG,'(A,21X,2(A,7X),A,1X,A)') '%#', 'wall(y=-1)','wall(y=+1)','Two walls average','Error(\%)=(W1-W2)/(W1+W2)*100'
        WRITE(TECFLG, '(A,4ES17.9) ') 'Cf0                 =',Cf0_io(1),  Cf0_io(2),   Cf0_ave_IO, &
                                                              dabs(Cf0_io(1)-Cf0_io(2))/(Cf0_io(1)+Cf0_io(2))*100.0_WP
        WRITE(TECFLG, '(A,3ES17.9) ') 'Cfbk                =',Cfbk_io(1), Cfbk_io(2),  Cfbk_ave_IO
        WRITE(TECFLG, '(A)') '  '
        WRITE(TECFLG, '(A)         ') '%########Wall Shear, undimensional ###########################'
        WRITE(TECFLG,'(A,21X,2(A,7X),A)') '%#', 'wall(y=-1)','wall(y=+1)'
        WRITE(TECFLG, '(A,2ES17.9) ') 'Tau                 =',Tauw_io(1:2)
        WRITE(TECFLG, '(A,2ES17.9) ') 'Utau                =',Utaw_io(1:2)
        WRITE(TECFLG, '(A,2ES17.9) ') 'Re_tau              =',Ret_io(1:2)
        WRITE(TECFLG, '(A,2ES17.9) ') 'delta_nu(inner L)   =',InnerL(1:2)
        WRITE(TECFLG,'(A)') '  '
        
        IF(ICASE==ICHANL) THEN
            dyMIN = 1.0_wp/DYFI(1)
            dyMAX = 1.0_WP/DYFI(NCL2/2)
        ELSE
            dyMIN = 1.0_wp/DYFI(NCL2)
            dyMAX = 1.0_WP/DYFI(1)
        END IF
        
        dxplus(1:2) = DX * Ret_io(1:2)/Ldist_io(1:2)
        dzplus(1:2) = DZ * Ret_io(1:2)/Ldist_io(1:2)
        
        yMINplus(1:2) = dyMIN * Ret_io(1:2)/Ldist_io(1:2)
        yMAXplus(1:2) = dyMAX * Ret_io(1:2)/Ldist_io(1:2)
        
          
        WRITE(TECFLG, '(A)         ') '%########Real mesh resolution ###############################'
        WRITE(TECFLG,'(A,4X,2(A,7X),A,1X,A)') '%#Based on Utau at', ' wall(y=-1)','wall(y=+1), and original'
        WRITE(TECFLG, '(A,2ES17.9) ') 'dx+                 =',dxplus(1:2)
        WRITE(TECFLG, '(A,2ES17.9) ') 'dz+                 =',dzplus(1:2)
        WRITE(TECFLG, '(A,2ES17.9) ') 'dy1+                =',yMINplus(1:2)
        WRITE(TECFLG, '(A,2ES17.9) ') 'dy1max+             =',yMAXplus(1:2)
        
        IF(thermlflg==1) THEN
        
            ALLOCATE (dxplusj(NCL2,3))
            ALLOCATE (dyplusj(NCL2,3))
            ALLOCATE (dzplusj(NCL2,3))
            OPEN(TECFLG2,FILE=TRIM(filepath4)//'Result.IO.SemiLocal.Yplus.'//TRIM(PNTIM)//'.tec')
            WRITE(TECFLG2,'(A)') 'TITLE = " Local Yplus (26 variables)" '
            
            write(TECFLG2,'(A)') 'VARIABLES = "YCC", "dxn+", "dyn+", "dzn+", "dxp+", "dyhp+", "dzp+", "dxa+", "dya+", "dza+" '
            WRITE(TECFLG2,'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
            DO J=1, NCL2
                
                dxplusj(J,1) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_io(1) )/M1xztL_F0_io(J)*DX
                dyplusj(J,1) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_io(1) )/M1xztL_F0_io(J)/DYFI(J)
                dzplusj(J,1) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_io(1) )/M1xztL_F0_io(J)*DZ
                
                dxplusj(J,2) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_io(2) )/M1xztL_F0_io(J)*DX
                dyplusj(J,2) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_io(2) )/M1xztL_F0_io(J)/DYFI(J)
                dzplusj(J,2) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_io(2) )/M1xztL_F0_io(J)*DZ
                
                dxplusj(J,3) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_ave_io )/M1xztL_F0_io(J)*DX
                dyplusj(J,3) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_ave_io )/M1xztL_F0_io(J)/DYFI(J)
                dzplusj(J,3) = REN*DSQRT(D1xztL_F0_io(J)*Tauw_ave_io )/M1xztL_F0_io(J)*DZ
                WRITE(TECFLG2,'(10ES17.9)') YCC(J), dxplusj(J,1), dyplusj(J,1), dzplusj(J,1), &
                                                    dxplusj(J,2), dyplusj(J,2), dzplusj(J,2), &
                                                    dxplusj(J,3), dyplusj(J,3), dzplusj(J,3)
                                                    
            END DO
            CLOSE(TECFLG2)
            
!            dxplus(1:2) = REN*DSQRT(DWAL(1:2)*Tauw_io(1:2) )/MWal(1:2)*DX
!            dzplus(1:2) = REN*DSQRT(DWAL(1:2)*Tauw_io(1:2) )/MWal(1:2)*DZ
        
!            yMINplus(1:2) = dyMIN * REN*DSQRT(DWAL(1:2)*Tauw_io(1:2) )/MWal(1:2)
!            yMAXplus(1:2) = dyMAX * REN*DSQRT(DWAL(1:2)*Tauw_io(1:2) )/MWal(1:2)
            
!            WRITE(TECFLG, '(A)         ') '%########Real mesh resolution ###############################'
!            WRITE(TECFLG,'(A,4X,2(A,7X),A,1X,A)') '%#Based on Utau at', ' wall(y=-1)','wall(y=+1), and original'
!            WRITE(TECFLG, '(A,2ES17.9) ') 'dx+                 =',dxplus(1:2)
!            WRITE(TECFLG, '(A,2ES17.9) ') 'dz+                 =',dzplus(1:2)
!            WRITE(TECFLG, '(A,2ES17.9) ') 'dy1+                =',yMINplus(1:2)
!            WRITE(TECFLG, '(A,2ES17.9) ') 'dy1max+             =',yMAXplus(1:2)
        
            
            WRITE(TECFLG,'(A,4X,2(A,7X),A,1X,A)') '%#Based on Utau_w at', ' wall(y=-1)','wall(y=+1), and semi-local'
            WRITE(TECFLG,'(A)'                  )              'min.    max.     min.  max.'
            WRITE(TECFLG, '(A,4ES17.9) ') 'dx+                 =',MINVAL(dxplusj(:,1)), MAXVAL(dxplusj(:,1)), &
                                                                  MINVAL(dxplusj(:,2)), MAXVAL(dxplusj(:,2))
            WRITE(TECFLG, '(A,4ES17.9) ') 'dz+                 =',MINVAL(dzplusj(:,1)), MAXVAL(dzplusj(:,1)), &
                                                                  MINVAL(dzplusj(:,2)), MAXVAL(dzplusj(:,2))
            WRITE(TECFLG, '(A,4ES17.9) ') 'dy+                 =',MINVAL(dyplusj(:,1)), MAXVAL(dyplusj(:,1)), &
                                                                  MINVAL(dyplusj(:,2)), MAXVAL(dyplusj(:,2))
                                                                  
            WRITE(TECFLG,'(A,4X,2(A,7X),A,1X,A)') '%#Based on Utau_ave at', ' wall(y=-1)','wall(y=+1), and semi-local'
            WRITE(TECFLG,'(A)'                  )              'min.    max.     min.  max.'
            WRITE(TECFLG, '(A,2ES17.9) ') 'dx+                 =',MINVAL(dxplusj(:,3)), MAXVAL(dxplusj(:,3))
            WRITE(TECFLG, '(A,2ES17.9) ') 'dz+                 =',MINVAL(dzplusj(:,3)), MAXVAL(dzplusj(:,3))
            WRITE(TECFLG, '(A,2ES17.9) ') 'dy+                 =',MINVAL(dyplusj(:,3)), MAXVAL(dyplusj(:,3))
        
        END IF

        
        !=============================================================
        ddenintg = 0.0_WP
        denmintg = 0.0_WP
        bdfcintg= 0.0_WP
        densintg= 0.0_WP
        DO J=1, NCL2
            ! second order intgeral 
!            IF(J==1) THEN
!                DENtemp =0.5_WP*( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J) + YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) +             &
!                                  Dwal(1) )
!            ELSE IF(J==NCL2) THEN
!                DENtemp =0.5_WP*( Dwal(2) +             &
!                                  ( YCL2ND_WFF(J)  *D1xztL_F0_io(J) + YCL2ND_WFB(J)  *D1xztL_F0_io(J-1) ) )
!            ELSE
!                DENtemp =0.5_WP*( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J) + YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) +             &
!                                  ( YCL2ND_WFF(J)  *D1xztL_F0_io(J) + YCL2ND_WFB(J)  *D1xztL_F0_io(J-1) ) )
!            END IF
            !first order integral
            DENtemp = D1xztL_F0_io(J)
            
            denmintg = denmintg+DENtemp/DYFI(J)
            densintg(J) = denmintg !; !WRITE(*,*) J, DENtemp, densintg(J)
            
            ddenintg = ddenintg + (DENtemp-DenAvew)/DYFI(J)
            bdfcintg(J) = F_A*ddenintg
        END DO
        CALL CHKHDL('      ==>Calculated bodyforce distribution',myid)
        
        
!        WRITE(TECFLG,'(A)') '  '
!        if(thermlflg==1) THEN
!            WRITE(TECFLG, '(A)         ') '%########Wall Shear, dimensional ########################'
!            WRITE(TECFLG,'(A,21X,2(A,7X),A)') '%#', 'wall(y=-1)','wall(y=+1)'
!            WRITE(TECFLG, '(A,2ES17.9) ') 'Tau(Kg/m/s^2)       =',Tauw_d_io(1), Tauw_d_io(2)
!            WRITE(TECFLG, '(A,2ES17.9) ') 'Utau(m/s)           =',Utaw_d_io(1), Utaw_d_io(2)
!        END IF
        
        CLOSE(TECFLG)
        
        IF(thermlflg==1) THEN
            DEALLOCATE (dxplusj)
            DEALLOCATE (dyplusj)
            DEALLOCATE (dzplusj)
        END IF
        
        
        RETURN
    END SUBROUTINE 
    
    
    SUBROUTINE WRT_Checking_TABLE_XZ_IO
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: TECFLG = 200
        INTEGER(4)    :: I
        
        
        
        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
        OPEN(TECFLG,FILE=TRIM(filepath4)//'Result.IO.Table.WallandBulk.'//TRIM(PNTIM)//'.tec', POSITION='APPEND')
        
        WRITE(TECFLG, '(A)         ') '%########Reference States (dimensionaless)###################'   
        WRITE(TECFLG, '(A, 1ES17.9)') 'REF Re                       =', REN
        WRITE(TECFLG,'(A)') '  '
        WRITE(TECFLG, '(A)         ') '%########Checking force balance (wall unit scaled) ##########'
        WRITE(TECFLG, '(A,2ES17.9) ') 'Total Force in X direction (RA and FA)= ',NSFbalt_RA(1),NSFbalt_FA(1)
        WRITE(TECFLG, '(A,2ES17.9) ') 'Total Force in Y direction (RA and FA)= ',NSFbalt_RA(2),NSFbalt_FA(2)
        WRITE(TECFLG, '(A,2ES17.9) ') 'Total Force in Z direction (RA and FA)= ',NSFbalt_RA(3),NSFbalt_FA(3)
        WRITE(TECFLG,'(A)') '  '
        
        
        if( thermlflg==1 ) THEN
            WRITE(TECFLG, '(A)         ') '%########Checking energy balance in #########################'
            WRITE(TECFLG, '(A,2ES17.9) ') 'Total energy in(FA)                   = ',ENEbalt_FA
            WRITE(TECFLG,'(A)') '  '
            
            WRITE(TECFLG,'(A)') '  '
            WRITE(TECFLG, '(A)         ') '%########Location of Tbk ######################################'
            WRITE(TECFLG, '(A, 1I5.1  )') 'J index  for Tbk       =', J4Tbk
            WRITE(TECFLG, '(A, 1ES17.9)') 'y  for Tpc             =', YCC(J4Tbk)
            WRITE(TECFLG, '(A)') '   '
            WRITE(TECFLG, '(A)         ') '%########Location of Tpc ######################################'
            WRITE(TECFLG, '(A, 1ES17.9)') 'Tpc(K)                 =', T4CpMax*T0
            WRITE(TECFLG, '(A, 1ES17.9)') 'Tpc(C)                 =', T4CpMax*T0-TK2C
            WRITE(TECFLG, '(A, 1ES17.9)') 'Tpc(undim)             =', T4CpMax
            IF(DMAX1(Twal(1),Twal(2)).GT.T4CpMax) THEN
                WRITE(TECFLG, '(A, 1I5.1  )') 'J index  for Tpc       =', J4Tpc
                WRITE(TECFLG, '(A, 1ES17.9)') 'y  for Tpc             =', YCC(J4Tpc)
                IF(YCC(J4Tpc).GT.0.0_wp) THEN
                    WRITE(TECFLG, '(A, 1ES17.9)') 'y+                     =', (1.0_wp-dabs(YCC(J4Tpc)))*Ret_io(2)/Ldist_io(2)
                ELSE
                    WRITE(TECFLG, '(A, 1ES17.9)') 'y+                     =', (1.0_wp-dabs(YCC(J4Tpc)))*Ret_io(1)/Ldist_io(1)
                END IF
            ELSE
                WRITE(TECFLG, '(A)'         ) 'The current Temperature ange does not cover Tpc.'
            END IF
            WRITE(TECFLG,'(A)') '  '
            
            WRITE(TECFLG, '(A)         ') '%########sided bulk average(Tbuk seperation)###################'
            WRITE(TECFLG,'(A)          ') '%#            (based on sided bulk T, K and Qw (Kasagi1996))#'
            WRITE(TECFLG, '(A,2X,I3.1,2ES17.9)')'Jindex,Ycc,  =', J4Tbk, YCC(J4Tbk)
            WRITE(TECFLG,'(A,15X,2(A,7X),A)')   '%#', 'side(y=-1)','side(y=+1)','gloabl'
            WRITE(TECFLG, '(A, 3ES17.9)')       'Y+max        =', (YCC(J4Tbk)+1.0_WP)*Ret_io(1)/Ldist_io(1), &
                                                                  (1.0_WP-YCC(J4Tbk))*Ret_io(2)/Ldist_io(2)
            WRITE(TECFLG,'(A)') '  '
            WRITE(TECFLG, '(A)         ') '%########based on Int{\rho}dy (arithmetic mean Density) #####'
            WRITE(TECFLG, '(A,1ES17.9) ') 'Tau     =0.5*(Tauw1+Tauw2)        =',Tauw_ave_io
            WRITE(TECFLG, '(A,1ES17.9) ') 'GravityF=2.0*F_A*Int{RHO}dy       =',DenAvew*2.0_WP*F_A
            WRITE(TECFLG, '(A,1ES17.9) ') 'Density =1/Ly*Int{RHO}dy          =',DenAvew
            WRITE(TECFLG, '(A,1ES17.9) ') 'Viscity =1/Ly*Int{Mu}dy           =',VisAvew
            WRITE(TECFLG, '(A,1ES17.9) ') 'Utau    =SQRT(Tau_ave/DenAve)     =',Utaw_ave_io
            WRITE(TECFLG, '(A,1ES17.9) ') 'Re_tau  =REN*Utawave*DenAve/MuAve =',Ret_ave_io
        END IF
        
        CLOSE(TECFLG)
        
        RETURN
    END SUBROUTINE 
!===============================================================================    
    SUBROUTINE PP_SSzero_SIDED(FLG)
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        INTEGER(4) :: FLG
        INTEGER(4) :: J
        REAL(WP)   :: SSmax, SSFA, Umax, UFA
        REAL(WP)   :: Mdott(2), Hdott(2), Ltt(2), Gbukt(2), Hbukt(2)
        REAL(WP)   :: Htemp, Ttemp, Dtemp, Ktemp, Ptemp, Mtemp, Cptemp, Btemp
        
        
        !============sided by zero stress balance===================================
        !=========Find the zero stress balance value.==================================
        SSmax = 1.0E14_wp
        DO J=1, NCL2
            IF(FLG==1 .and. thermlflg.eq.1) THEN
                SSFA = Tau_Mean_RA(J,1,2)-uff2d_FA(J,1,2)+bdfcintg(J)
            ELSE IF (FLG==2 .and. thermlflg.ne.1) THEN
                SSFA = Tau_Mean_RA(J,1,2)-uf2d_RA(J,1,2)+bdfcintg(J)
            END IF
            IF (dabs(SSFA).LT.SSmax) THEN
                SSmax = SSFA 
                J4SS0 = J
            END IF
            
        END DO
        !WRITE(*,*) '# zero ss', J4SS0, Ycc(J4SS0)
        
        DO J=1, NCL2
            IF(J.LT.J4SS0) THEN ! near y=-1
                TauwSD(J) = Tauw_io(1)
                DensSD(J) = Dwal(1)
                ViscSD(J) = Mwal(1)
                YWdiSD(J) = YCC(J)+1.0_WP
            else                ! near y=+1
                TauwSD(J) = Tauw_io(2)
                DensSD(J) = Dwal(2)
                ViscSD(J) = Mwal(2)
                YWdiSD(J) = 1.0_WP-YCC(J)
            end if
        END DO  
        
        IF(thermlflg==1) THEN
            DO J=1, NCL2
                IF(J.LT.J4SS0) THEN ! near y=-1
                    CpSD(J)= Cpwal(1)
                    QwSD(J)= Qw(1)
                    TwSD(J)= Twal(1)
                    HwSD(J)= Hwal_RA(1)
                else                ! near y=+1
                    CpSD(J)= Cpwal(2)  
                    QwSD(J)= Qw(2)
                    TwSD(J)= Twal(2)
                    HwSD(J)= Hwal_RA(2)
                end if
            END DO
            
            !=====sided bulk g and h=============
            Mdott  = 0.0_WP
            Hdott  = 0.0_WP
            UbukSsd = 0.0_WP
            MbukSsd = 0.0_WP
            DbukSsd = 0.0_WP
            Ltt    = 0.0_WP
            
            DO J=1,J4SS0
                Mdott(1) = Mdott(1) + G1xztL_F0_io(J,1)/DYFI(J)
                Hdott(1) = Hdott(1) + GHxztL_F0_io(J,1)/DYFI(J)
                UbukSsd(1) = UbukSsd(1) + U1xztL_F0_io(J,1)/DYFI(J)
                MbukSsd(1) = MbukSsd(1) + M1xztL_F0_io(J)/DYFI(J)
                DbukSsd(1) = DbukSsd(1) + D1xztL_F0_io(J)/DYFI(J)
                Ltt(1)   = Ltt(1)   + 1.0_wp/DYFI(J)
            ENDDO
            DO J=J4SS0,NCL2
                Mdott(2) = Mdott(2) + G1xztL_F0_io(J,1)/DYFI(J)
                Hdott(2) = Hdott(2) + GHxztL_F0_io(J,1)/DYFI(J)
                UbukSsd(2) = UbukSsd(2) + U1xztL_F0_io(J,1)/DYFI(J)
                MbukSsd(2) = MbukSsd(2) + M1xztL_F0_io(J)/DYFI(J)
                DbukSsd(2) = DbukSsd(2) + D1xztL_F0_io(J)/DYFI(J)
                Ltt(2)   = Ltt(2)   + 1.0_wp/DYFI(J)
            ENDDO
    
            Gbukt(1:2) =  Mdott(1:2) /Ltt(1:2)
            Hbukt(1:2) =  Hdott(1:2) /Mdott(1:2)
            UbukSsd(1:2) =  UbukSsd(1:2) /Ltt(1:2)
            MbukSsd(1:2) =  MbukSsd(1:2) /Ltt(1:2)
            DbukSsd(1:2) =  DbukSsd(1:2) /Ltt(1:2)

            !write(*,'(A,4ES13.5)') '# Mdot, Hdot, Gbuk, Hbuk', Mdot, Hdot, Gbuk, Hbuk
            !===============search===table==========
            DO J=1, 2
                HTEMP = Hbukt(J)
                        
                IF(thermoStat==search_table) Then
                    call NIST_SLEVAL_HT(Htemp,Ttemp)
                    !call NIST_SLEVAL_TD(Ttemp,Dtemp)
                    call NIST_SLEVAL_TK(Ttemp,Ktemp)
                    !call NIST_SLEVAL_TM(Ttemp,Mtemp)
                    call NIST_SLEVAL_TCP(Ttemp,CPtemp)
                    call NIST_SLEVAL_TB(Ttemp,Btemp)
                ELSE IF(thermoStat==idealgas_law) Then
                    Ptemp = 1.0_wp
                    !call idealgas_HT(Htemp,Ttemp)
                    !call idealgas_TD(Ttemp,Dtemp,Ptemp)
                    !call idealgas_TK(Ttemp,Ktemp)
                    !call idealgas_TB(Ttemp,Btemp)
                ELSE
                END IF
                TbukSsd(J) = Ttemp
                KbukSsd(J) = Ktemp
                !DbukSsd(J) = Dtemp
                !MbukSsd(J) = Mtemp
                CpbkSsd(J) = Cptemp
                BbukSsd(J) = Btemp
                !UbukSsd(J) = Gbukt(J)/Dtemp
            END DO
        
            RebkSsd(1:2) = DbukSsd(1:2)*UbukSsd(1:2)*Ldist_io(1:2)/MbukSsd(1:2) * REN ! based on diameter or the whole height of channel.
            PrbkSsd(1:2) = MbukSsd(1:2)*CpbkSsd(1:2)/KbukSsd(1:2)  * Prt0
            
            !WRITE(*,*) '# RebkSsd',RebkSsd(1:2)
            !WRITE(*,*) '# PrbkSsd',PrbkSsd(1:2)
            
        END IF
            
        RETURN
    END SUBROUTINE
    
    SUBROUTINE PP_J4TbukTsd
        USE VARS_AVERAGED_XZ_IO
        
        IMPLICIT NONE
        INTEGER(4) :: J
        REAL(WP)   :: TTmax, TTdif
        INTEGER(4) :: FLG
        REAL(WP)   :: SSmax, SSFA, Umax, UFA
        REAL(WP)   :: Mdott(2), Hdott(2), Ltt(2), Gbukt(2), Hbukt(2)
        REAL(WP)   :: Htemp, Ttemp, Dtemp, Ktemp, Ptemp, Mtemp, Cptemp, Btemp

        
        !=========Find the Tbulk==================================
        TTmax = 1.0E14_wp
        DO J=1, NCL2
            TTdif = dabs(T1xztL_F0_io(J)-Tbuk)
            IF (TTdif.LT.TTmax) THEN
                TTmax = TTdif 
                J4Tbk = J
            END IF
        END DO
        
        L4Tbk(1) = YCC(J4Tbk)-(-1.0_WP)
        L4Tbk(2) = 1.0_wp - YCC(J4Tbk)
        !WRITE(*,*) '# J4Tbk', J, L4Tbk(1), L4Tbk(2)
        
        Mdott  = 0.0_WP
        Hdott  = 0.0_WP
        UbukTsd = 0.0_WP
        MbukTsd = 0.0_WP
        DbukTsd = 0.0_WP
        Ltt    = 0.0_WP
        
        DO J=1,J4Tbk
            Mdott(1) = Mdott(1) + G1xztL_F0_io(J,1)/DYFI(J)
            Hdott(1) = Hdott(1) + GHxztL_F0_io(J,1)/DYFI(J)
            UbukTsd(1) = UbukTsd(1) + U1xztL_F0_io(J,1)/DYFI(J)
            MbukTsd(1) = MbukTsd(1) + M1xztL_F0_io(J)/DYFI(J)
            DbukTsd(1) = DbukTsd(1) + D1xztL_F0_io(J)/DYFI(J)
            Ltt(1)   = Ltt(1)   + 1.0_wp/DYFI(J)
        ENDDO
        DO J=J4Tbk,NCL2
            Mdott(2) = Mdott(2) + G1xztL_F0_io(J,1)/DYFI(J)
            Hdott(2) = Hdott(2) + GHxztL_F0_io(J,1)/DYFI(J)
            UbukTsd(2) = UbukTsd(2) + U1xztL_F0_io(J,1)/DYFI(J)
            MbukTsd(2) = MbukTsd(2) + M1xztL_F0_io(J)/DYFI(J)
            DbukTsd(2) = DbukTsd(2) + D1xztL_F0_io(J)/DYFI(J)
            Ltt(2)   = Ltt(2)   + 1.0_wp/DYFI(J)
        ENDDO

        Gbukt(1:2) =  Mdott(1:2) /Ltt(1:2)
        Hbukt(1:2) =  Hdott(1:2) /Mdott(1:2)
        UbukTsd(1:2) =  UbukTsd(1:2) /Ltt(1:2)
        MbukTsd(1:2) =  MbukTsd(1:2) /Ltt(1:2)
        DbukTsd(1:2) =  DbukTsd(1:2) /Ltt(1:2)

        !write(*,'(A,4ES13.5)') '# Mdot, Hdot, Gbuk, Hbuk', Mdot, Hdot, Gbuk, Hbuk
        !===============search===table==========
        DO J=1, 2
            HTEMP = Hbukt(J)
                    
            IF(thermoStat==search_table) Then
                call NIST_SLEVAL_HT(Htemp,Ttemp)
                !call NIST_SLEVAL_TD(Ttemp,Dtemp)
                call NIST_SLEVAL_TK(Ttemp,Ktemp)
                !call NIST_SLEVAL_TM(Ttemp,Mtemp)
                call NIST_SLEVAL_TCP(Ttemp,CPtemp)
                call NIST_SLEVAL_TB(Ttemp,Btemp)
            ELSE IF(thermoStat==idealgas_law) Then
                Ptemp = 1.0_wp
                !call idealgas_HT(Htemp,Ttemp)
                !call idealgas_TD(Ttemp,Dtemp,Ptemp)
                !call idealgas_TK(Ttemp,Ktemp)
                !call idealgas_TB(Ttemp,Btemp)
            ELSE
            END IF
            TbukTsd(J) = Ttemp
            KbukTsd(J) = Ktemp
            !DbukTsd(J) = Dtemp
            !MbukTsd(J) = Mtemp
            CpbkTsd(J) = Cptemp
            BbukTsd(J) = Btemp
            !UbukTsd(J) = Gbukt(J)/Dtemp
        END DO
        
        
        !============sided by zero stress balance===================================
        !=========Find the zero stress balance value.==================================
        TTmax = 1.0E14_wp
    ! find where TbukSsd J
        DO J=1, J4SS0
            TTdif = TbukSsd(1)-T1xztL_F0_io(J)
            IF (dabs(TTdif).LT.TTmax) THEN
                TTmax = TTdif
                J4TbukSsd(1) = J
            END IF
        END DO
        TTmax = 1.0E14_wp
        DO J=J4SS0, NCL2
            TTdif = TbukSsd(2)-T1xztL_F0_io(J)
            IF (dabs(TTdif).LT.TTmax) THEN
                TTmax = TTdif
                J4TbukSsd(2) = J
            END IF
        END DO
        !WRITE(*,*) '# J4TbukSsd', J4TbukSsd(1:2), Ycc(J4TbukSsd(1:2)), TbukSsd(1:2)
        
        L4TbkSsd(1) = YCC(J4TbukSsd(1))-(-1.0_WP)
        L4TbkSsd(2) = 1.0_wp - YCC(J4TbukSsd(2))
        
            
    RETURN
    END SUBROUTINE 
        
       

    
    !===============================================================================    
    SUBROUTINE PP_Umax_SIDED
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        INTEGER(4) :: J
        REAL(WP)   :: SSmax, SSFA, Umax, UFA
        REAL(WP)   :: Mdott(2), Hdott(2), Ltt(2), Gbukt(2), Hbukt(2)
        REAL(WP)   :: Htemp, Ttemp, Dtemp, Ktemp, Ptemp
        
!===============sd by max. velocity========================
        !=========Find the max. Velocity==================================
        Umax = 0.0_wp
        DO J=1, NCL2
            UFA = G1xztL_F0_io(J,1)/D1xztL_F0_io(J)
            IF (UFA.GT.Umax) THEN
                Umax = UFA 
                J4MaxU = J
            END IF
        END DO
        

        RETURN
    END SUBROUTINE
    
    
   
    
    SUBROUTINE WRT_FLOW_budgets_Profile_XZ_IO(STR)
        USE VARS_AVERAGED_XZ_IO
        IMPLICIT NONE
        
        CHARACTER(5),INTENT(IN)   :: STR  ! Favre or Reynd
        CHARACTER(128) :: FLNM, FLNN
        CHARACTER(5)   :: STDIM(2)
        CHARACTER(4)   :: STFASTRESS(8)
        CHARACTER(2)   :: SJ2
        REAL(WP)   :: scaling1, scaling2, scaling3, scaling4,scaling5, scaling6, scaling7
        REAL(WP)   :: tke,Ruv_vis
        INTEGER(4)    :: I, J, N, L, TECFLG(2), TECFLG3, Nmax
        REAL(wp)   ::   COE, tem, Mtemp
        REAL(WP)   ::   TKE2, EPPSI, Fmu
        REAL(wp)   ::   Cmu, K2De
        
        REAL(WP)   :: FCT(0:NND2,NDV)
        REAL(WP)   :: COE1, COE2, DENtemp, Dintg, intgbdfc
        
        IF(MYID.EQ.0) CALL CHKRLHDL('     **>WRT_OUT_BUDGET '//STR// ' at ',myid,phyTIME_io)
!============================budgets=============================================
        !============================budgets=============================================
        if(TRIM(STR)=='Favre') THEN
            STFASTRESS(1) = 'FAuu'
            STFASTRESS(2) = 'FAuv'
            STFASTRESS(3) = 'FAuw'
            STFASTRESS(4) = 'FAvv'
            STFASTRESS(5) = 'FAvw'
            STFASTRESS(6) = 'FAww'
            STFASTRESS(7) = 'FATK'
            STFASTRESS(8) = 'FAMK'
        else if (TRIM(STR)=='Reynd') THEN
            STFASTRESS(1) = 'RAuu'
            STFASTRESS(2) = 'RAuv'
            STFASTRESS(3) = 'RAuw'
            STFASTRESS(4) = 'RAvv'
            STFASTRESS(5) = 'RAvw'
            STFASTRESS(6) = 'RAww'
            STFASTRESS(7) = 'RATK'
            STFASTRESS(8) = 'RAMK'
        else
        end if
        
        Nmax = 1
        STDIM(1) = 'undim'
        TECFLG(1) = 101
        IF(ppdim==1) THEN
            Nmax = 2
            STDIM(2) = 'dimen'
            TECFLG(2) = 102
        END IF
        
        DO L=1, 8
        
            DO N=1, Nmax
                WRITE(PNTIM,'(1ES15.9)') phyTIME_io
                FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(N))//'.Profile.BUDGETS.'//STR//'.'//TRIM(STFASTRESS(L))//  &
                      '.'//TRIM(PNTIM)//'.tec'
                OPEN (TECFLG(N), FILE=TRIM(ADJUSTL(FLNM)))
                WRITE(TECFLG(N),'(A)') 'TITLE = "'//STR//' Averged Flow (28 variables)" '
                J=0;                          write(TECFLG(N),'(A)',advance="no") 'VARIABLES = '
                
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'Tauw",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'Dw",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'DInt",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'D",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'Muw",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'MuInt",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'M",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'F_A",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'YCC",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'Ywd",' 
                
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'prodc_stres",' 
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'viscs_dissp",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'pdudx_stran",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'Turbu_diffu",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'dpudx_diffu",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'viscs_diffu",'

                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'press_accel1",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'viscs_accel1",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'prodc_dgbdf1",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'prodc_drvfc1",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'balance1"'

                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'turss_accl2",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'prodc_dgbdf2",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'prodc_drvfc2",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'balance2"'

                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'pressure3",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)',advance="no") '"'//SJ2//'vistress3",'
                J=J+1; WRITE(SJ2,'(1I2.2)') J; write(TECFLG(N),'(A)'             ) '"'//SJ2//'balance3"'
                
                WRITE(TECFLG(N),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
            END DO
        
            DO J=1,NCL2
                IF (L==7) THEN
                    WRITE(TECFLG(1),'(28ES20.12)') &
                    TauwSD(J), &
                    DensSD(J), DenAvew, D1xztL_F0_io(J), &
                    ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                    F_A, YCC(J), YWdiSD(J), &
                    BUDG_prodc_stres_tke(J), BUDG_viscs_dissp_tke(J), BUDG_pdudx_stran_tke(J), &
                    BUDG_Turbu_diffu_tke(J), BUDG_dpudx_diffu_tke(J), BUDG_viscs_diffu_tke(J), &
                    BUDG_press_accl1_tke(J), BUDG_viscs_accl1_tke(J), 0.0_WP, BUDG_prodc_dvfc1_tke(J), BUDG_balance1_tke(J), &
                    BUDG_turss_accl2_tke(J), BUDG_prodc_gvfc2_tke(J),         BUDG_prodc_dvfc1_tke(J), BUDG_balance2_tke(J), &
                    BUDG_pressure3_tke(J),   BUDG_vistress3_tke(J),   BUDG_balance3_tke(J)
                        
                ELSE IF (L==8) THEN
                    WRITE(TECFLG(1),'(28ES20.12)') &
                    TauwSD(J), &
                    DensSD(J), DenAvew, D1xztL_F0_io(J), &
                    ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                    F_A, YCC(J), YWdiSD(J), &
                    BUDG_prodc_stres_mke(J), BUDG_viscs_dissp_mke(J), BUDG_pdudx_stran_mke(J), &
                    BUDG_Turbu_diffu_mke(J), BUDG_dpudx_diffu_mke(J), BUDG_viscs_diffu_mke(J), &
                    BUDG_press_accl1_mke(J), BUDG_viscs_accl1_mke(J), 0.0_WP, BUDG_prodc_dvfc1_mke(J), BUDG_balance1_mke(J), &
                    BUDG_turss_accl2_mke(J), BUDG_prodc_gvfc2_mke(J),         BUDG_prodc_dvfc1_mke(J), BUDG_balance2_mke(J), &
                    BUDG_pressure3_mke(J),   BUDG_vistress3_mke(J),   BUDG_balance3_mke(J)
                ELSE
                    WRITE(TECFLG(1),'(28ES20.12)') &
                    TauwSD(J), &
                    DensSD(J), DenAvew, D1xztL_F0_io(J), &
                    ViscSD(J), VisAvew, M1xztL_F0_io(J), &
                    F_A, YCC(J), YWdiSD(J), &
                    BUDG_prodc_stres_duiuj(J,L), BUDG_viscs_dissp_duiuj(J,L), BUDG_pdudx_stran_duiuj(J,L), &
                    BUDG_Turbu_diffu_duiuj(J,L), BUDG_dpudx_diffu_duiuj(J,L), BUDG_viscs_diffu_duiuj(J,L), &
                    BUDG_press_accl1_duiuj(J,L), BUDG_viscs_accl1_duiuj(J,L), 0.0_WP, &
                    BUDG_prodc_dvfc1_duiuj(J,L), BUDG_balance1_duiuj(J,L), &
                    BUDG_turss_accl2_duiuj(J,L), BUDG_prodc_gvfc2_duiuj(J,L),         &
                    BUDG_prodc_dvfc1_duiuj(J,L), BUDG_balance2_duiuj(J,L), &
                    BUDG_pressure3_duiuj(J,L),   BUDG_vistress3_duiuj(J,L),   BUDG_balance3_duiuj(J,L)
                END IF
!                IF(ppdim==1) THEN        
!                    WRITE(TECFLG(2),'(13ES20.12)') YCC(J)*L0, (1.0_wp-dabs(YCC(J)))*REN*COE, COE*U0, &
!                        BUDG_prodc_stres_duiuj(J,L)*scaling1, BUDG_Turbu_diffu_duiuj(J,L)*scaling1, &
!                        BUDG_dpudx_diffu_duiuj(J,L)*scaling1, BUDG_pdudx_stran_duiuj(J,L)*scaling1, &
!                        BUDG_press_accl1_duiuj(J,L)*scaling1, BUDG_viscs_accl1_duiuj(J,L)*scaling2, &
!                        BUDG_viscs_diffu_duiuj(J,L)*scaling2, BUDG_viscs_dissp_duiuj(J,L)*scaling2, &
!                        BUDG_prodc_gvfc2_duiuj(J,L)*D0*U0, BUDG_balance1_duiuj(J,L)*scaling2
!                END IF
            END DO
            
            DO N=1,Nmax
                CLOSE(TECFLG(N))
            END DO
        END DO
        
        RETURN
    END SUBROUTINE
    
    
    SUBROUTINE WRITE_SPECO_AVE_PROFILE(STR)
    !   Refer to: Page 213, Chapter 9.9e in Orland book.
        use mesh_info
        use init_info
        use SPECO_info
        use postprocess_info
        USE WRT_INFO
        
        IMPLICIT NONE
        CHARACTER(4),INTENT(IN)  :: STR
        
        CHARACTER(15) :: PNTIM
        CHARACTER(15) :: PNLOC
        INTEGER(4)    :: DFLG(MGRID)
        character(256) :: FLNAME
        !REAL(WP)      :: Ret_ave, U_tau_ave
        REAL(WP)      :: AKE, WDS(2)
        
        INTEGER(4)    :: N, JJ, L, KC, IC, M
        
        IF(MYID.NE.0) RETURN
        
        IF(TRIM(STR)=='FLOW') M=1
        IF(TRIM(STR)=='HEAT') M=2
        
        N3MH=NCL3/2+1
        N3MD=NCL3+2
        N1MH=NCL1_io/2+1
        N1MD=NCL1_io+2
        !CALL PP_wall_thermal_shear(flgxzt)
        
        !==========test for a ASCII output==============
        !Ret_ave  = 0.5_wp*(Ret_io(1)+Ret_io(2)) !;!WRITE(*,*) Ret_ave, Ret_io(1), Ret_io(2) !test
        !U_tau_ave= 0.5_wp*(Utaw_io(1)+Utaw_io(2))
        OPEN(100, FILE=TRIM(filepath4)//'CHK_PROBE_for_spectra_instant.dat')
        WRITE(100,'(A)') '## MGRID, JJ, YCC, Yplus1, Yplus2, Utaw_io(1:2), Ret_io(1:2), Ret_ave_io, Ldist_io(1:2)'
        DO N=1, MGRID
            JJ=JGMOV(N)
            DFLG(N) = 100+N
            
            WRITE(PNTIM,'(1ES15.9)') phyTIME_io
            WRITE(PNLOC,'(1I3.3)')   N
            
            !==============correlations======================================
            DO L=1,2
                IF(L==1) THEN
                    FLNAME= TRIM(filepath4)//'Result.IO.Spectral.'//TRIM(STR)//'.2PCorrelation.X.T'  &
                            //TRIM(PNTIM)//'.YLC'//TRIM(PNLOC)//'.tec'
                    OPEN(DFLG(N), FILE=TRIM(ADJUSTL(FLNAME)))
                    WRITE(DFLG(N),'(A)') 'TITLE = " Correlation (Streamwise separation)" '
                    WRITE(DFLG(N),'(A)',advance="no") 'VARIABLES = "1X", "U11", "U22", "U33", "U12", "U13", "U23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"O11", "O22", "O33", "O12", "O13", "O23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"VO11", "VO12", "VO13", '
                    WRITE(DFLG(N),'(A)',advance="no") '"VO21", "VO22", "VO23", '
                    WRITE(DFLG(N), '(A)'            ) '"VO31", "VO32", "VO33"  '
                    
                    WRITE(DFLG(N),'(A,1ES13.5,A,2ES13.5,A,3ES13.5,A)') &
                        'ZONE T= "Rkk(X) at y/delta= ',YCC(JJ), ' Utauw12= ',Utaw_io(1:2), ' Ret12= ',Ret_io(1:2),Ret_ave_io,' " '
                        
                    WDS(1) = YCC(JJ)-(-1.0_WP)
                    WDS(2) = 1.0_wp-YCC(JJ)
                    WRITE(100,'(2I4.2, 10ES13.5)') N, JJ, YCC(JJ), WDS(1:2)*Ret_io(1:2)/Ldist_io(1:2),  &
                                                Utaw_io(1:2), Ret_io(1:2), Ret_ave_io, Ldist_io(1:2) 
                                                
                    DO IC=1,N1MH
                        WRITE(DFLG(N),'(22ES13.5)') 0.5_WP*( XND_io(IC) + XND_io(IC+1) ), &
                                    R11X1_xztLa (JJ,IC,M), R22X1_xztLa (JJ,IC,M), R33X1_xztLa (JJ,IC,M), &
                                    R12X1_xztLa (JJ,IC,M), R13X1_xztLa (JJ,IC,M), R23X1_xztLa (JJ,IC,M), &
                                    V11X1_xztLa (JJ,IC,M), V22X1_xztLa (JJ,IC,M), V33X1_xztLa (JJ,IC,M), &
                                    V12X1_xztLa (JJ,IC,M), V13X1_xztLa (JJ,IC,M), V23X1_xztLa (JJ,IC,M), &
                                    VO11X1_xztLa(JJ,IC,M), VO12X1_xztLa(JJ,IC,M), VO13X1_xztLa(JJ,IC,M), &
                                    VO21X1_xztLa(JJ,IC,M), VO22X1_xztLa(JJ,IC,M), VO23X1_xztLa(JJ,IC,M), &
                                    VO31X1_xztLa(JJ,IC,M), VO32X1_xztLa(JJ,IC,M), VO33X1_xztLa(JJ,IC,M)
                    ENDDO
                    CLOSE(DFLG(N))
                END IF
                
                IF(L==2) THEN
                    FLNAME= TRIM(filepath4)//'Result.IO.Spectral.'//TRIM(STR)//'.2PCorrelation.Z.T'  &
                            //TRIM(PNTIM)//'.YLC'//TRIM(PNLOC)//'.tec'
                    OPEN(DFLG(N), FILE=TRIM(ADJUSTL(FLNAME)))
                    WRITE(DFLG(N),'(A)') 'TITLE = " Correlation (Spanwise separation)" '
                    WRITE(DFLG(N),'(A)',advance="no") 'VARIABLES = "3Z", "U11", "U22", "U33", "U12", "U13", "U23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"O11", "O22", "O33", "O12", "O13", "O23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"VO11", "VO12", "VO13", '
                    WRITE(DFLG(N),'(A)',advance="no") '"VO21", "VO22", "VO23", '
                    WRITE(DFLG(N), '(A)'            ) '"VO31", "VO32", "VO33"  '
                    WRITE(DFLG(N),'(A,1ES13.5,A,2ES13.5,A,3ES13.5,A)') &
                        'ZONE T= "Rkk(Z) at y/delta= ',YCC(JJ), ' Utauw12= ',Utaw_io(1:2), ' Ret12= ',Ret_io(1:2),Ret_ave_io,' " '
                        
                    DO KC=1,N3MH
                        WRITE(DFLG(N),'(22ES13.5)') 0.5_WP*( ZND(KC) + ZND(KC+1) ), &
                                R11X3_xztLa (JJ,KC,M), R22X3_xztLa (JJ,KC,M), R33X3_xztLa (JJ,KC,M), &
                                R12X3_xztLa (JJ,KC,M), R13X3_xztLa (JJ,KC,M), R23X3_xztLa (JJ,KC,M), &
                                V11X3_xztLa (JJ,KC,M), V22X3_xztLa (JJ,KC,M), V33X3_xztLa (JJ,KC,M), &
                                V12X3_xztLa (JJ,KC,M), V13X3_xztLa (JJ,KC,M), V23X3_xztLa (JJ,KC,M), &
                                VO11X3_xztLa(JJ,KC,M), VO12X3_xztLa(JJ,KC,M), VO13X3_xztLa(JJ,KC,M), &
                                VO21X3_xztLa(JJ,KC,M), VO22X3_xztLa(JJ,KC,M), VO23X3_xztLa(JJ,KC,M), &
                                VO31X3_xztLa(JJ,KC,M), VO32X3_xztLa(JJ,KC,M), VO33X3_xztLa(JJ,KC,M)
                    ENDDO
                    CLOSE(DFLG(N))
                
                END IF
            END DO
            
            !==============energy espectra======================================
            DO L=1,2
                IF(L==1) THEN
                    FLNAME= TRIM(filepath4)//'Result.IO.Spectral.'//TRIM(STR)//'.energy.Xwavenumber.T'  &
                            //TRIM(PNTIM)//'.YLC'//TRIM(PNLOC)//'.tec'
                    OPEN(DFLG(N), FILE=TRIM(ADJUSTL(FLNAME)))
                    WRITE(DFLG(N),'(A)') 'TITLE = " energy.spectra (Streamwise)" '
                    WRITE(DFLG(N),'(A)',advance="no") 'VARIABLES = "1IC", "2WaveNumberX"'  
                    WRITE(DFLG(N),'(A)',advance="no") '"E11", "E22", "E33", "E12", "E13", "E23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"H11", "H22", "H33", "H12", "H13", "H23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"EH11", "EH12", "EH13", '
                    WRITE(DFLG(N),'(A)',advance="no") '"EH21", "EH22", "EH23", '
                    WRITE(DFLG(N), '(A)'            ) '"EH31", "EH32", "EH33"  '
                    
                    WRITE(DFLG(N),'(A,1ES13.5,A,2ES13.5,A,3ES13.5,A)') &
                        'ZONE T= "Ekk(K1) at y/delta= ',YCC(JJ), ' Utauw12= ',Utaw_io(1:2), ' Ret12= ',Ret_io(1:2),Ret_ave_io,' " '
                    DO IC=1,N1MH
                        AKE=( DBLE(IC-1)*2.0_wp*PI/HX_io ) 
                        WRITE(DFLG(N),'(1I8.1, 22ES13.5)') IC, AKE, &
                            ENE11T_xztLa(JJ,IC,M), ENE22T_xztLa(JJ,IC,M), ENE33T_xztLa(JJ,IC,M), &
                            ENE12T_xztLa(JJ,IC,M), ENE13T_xztLa(JJ,IC,M), ENE23T_xztLa(JJ,IC,M), &
                            ENV11T_xztLa(JJ,IC,M), ENV22T_xztLa(JJ,IC,M), ENV33T_xztLa(JJ,IC,M), &
                            ENV12T_xztLa(JJ,IC,M), ENV13T_xztLa(JJ,IC,M), ENV23T_xztLa(JJ,IC,M), &
                            EVO11T_xztLa(JJ,IC,M), EVO12T_xztLa(JJ,IC,M), EVO13T_xztLa(JJ,IC,M), &
                            EVO21T_xztLa(JJ,IC,M), EVO22T_xztLa(JJ,IC,M), EVO23T_xztLa(JJ,IC,M), &
                            EVO31T_xztLa(JJ,IC,M), EVO32T_xztLa(JJ,IC,M), EVO33T_xztLa(JJ,IC,M)
                    ENDDO
                    CLOSE(DFLG(N))
                END IF
                
                IF(L==2) THEN
                    FLNAME= TRIM(filepath4)//'Result.IO.Spectral.'//TRIM(STR)//'.energy.Zwavenumber.T'  &
                            //TRIM(PNTIM)//'.YLC'//TRIM(PNLOC)//'.tec'
                    OPEN(DFLG(N), FILE=TRIM(ADJUSTL(FLNAME)))
                    WRITE(DFLG(N),'(A)') 'TITLE = " energy.spectra (spanwise)" '
                    WRITE(DFLG(N),'(A)',advance="no") 'VARIABLES = "3KC", "2WaveNumberZ"'  
                    WRITE(DFLG(N),'(A)',advance="no") '"E11", "E22", "E33", "E12", "E13", "E23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"H11", "H22", "H33", "H12", "H13", "H23", '
                    WRITE(DFLG(N),'(A)',advance="no") '"EH11", "EH12", "EH13", '
                    WRITE(DFLG(N),'(A)',advance="no") '"EH21", "EH22", "EH23", '
                    WRITE(DFLG(N), '(A)'            ) '"EH31", "EH32", "EH33"  '
                    WRITE(DFLG(N),'(A,1ES13.5,A,2ES13.5,A,3ES13.5,A)') &
                        'ZONE T= "Ekk(K3) at y/delta= ',YCC(JJ), ' Utauw12= ',Utaw_io(1:2), ' Ret12= ',Ret_io(1:2),Ret_ave_io,' " '
                    DO KC=1,N3MH
                        AKE=( DBLE(KC-1)*2.0_wp*PI/HZ ) 
                        WRITE(DFLG(N),'(1I8.1, 22ES13.5)') KC, AKE, &
                            ENE11Z_xztLa(JJ,KC,M), ENE22Z_xztLa(JJ,KC,M), ENE33Z_xztLa(JJ,KC,M), &
                            ENE12Z_xztLa(JJ,KC,M), ENE13Z_xztLa(JJ,KC,M), ENE23Z_xztLa(JJ,KC,M), &
                            ENV11Z_xztLa(JJ,KC,M), ENV22Z_xztLa(JJ,KC,M), ENV33Z_xztLa(JJ,KC,M), &
                            ENV12Z_xztLa(JJ,KC,M), ENV13Z_xztLa(JJ,KC,M), ENV23Z_xztLa(JJ,KC,M), &
                            EVO11Z_xztLa(JJ,KC,M), EVO12Z_xztLa(JJ,KC,M), EVO13Z_xztLa(JJ,KC,M), &
                            EVO21Z_xztLa(JJ,KC,M), EVO22Z_xztLa(JJ,KC,M), EVO23Z_xztLa(JJ,KC,M), &
                            EVO31Z_xztLa(JJ,KC,M), EVO32Z_xztLa(JJ,KC,M), EVO33Z_xztLa(JJ,KC,M)
                    ENDDO
                    CLOSE(DFLG(N))
                END IF
            END DO
            
        END DO
            
        
        RETURN
    END SUBROUTINE
    
    SUBROUTINE WRITE_SPECO_AVE_Contour(STR)
    !   Refer to: Page 213, Chapter 9.9e in Orland book.
        use mesh_info
        use init_info
        use SPECO_info
        use postprocess_info
        USE WRT_INFO
        
        IMPLICIT NONE
        CHARACTER(4),INTENT(IN)  :: STR
        
        CHARACTER(15) :: PNTIM
        CHARACTER(15) :: PNLOC
        INTEGER(4)    :: DFLG=101
        character(256) :: FLNAME
        !REAL(WP)      :: Ret_ave, U_tau_ave
        REAL(WP)      :: AKE
        REAL(WP)      :: yplus
         
        INTEGER(4)    :: N, JJ, L, KC, IC, JJM, JJC, M

        
        IF(MYID.NE.0) RETURN
        
        IF(TRIM(STR)=='FLOW') M=1
        IF(TRIM(STR)=='HEAT') M=2
        
        N3MH=NCL3/2+1
        N3MD=NCL3+2
        N1MH=NCL1_io/2+1
        N1MD=NCL1_io+2
        
        !==========test for a ASCII output==============
        !Ret_ave  = 0.5_wp*(Ret_io(1)+Ret_io(2)) !;!WRITE(*,*) Ret_ave, Ret_io(1), Ret_io(2) !test
        !U_tau_ave= 0.5_wp*(Utaw_io(1)+Utaw_io(2))
        
        !===============plane y-z=====================
        FLNAME = TRIM(filepath4)//'Result.IO.Spectral.'//TRIM(STR)//'.Contours.yz.'//TRIM(PNTIM)//'.tec' 
        OPEN (DFLG, FILE=TRIM(ADJUSTL(FLNAME)))

        WRITE(DFLG,'(A)') 'TITLE = "DNS FLOW YZ-plane"'
        WRITE(DFLG,'(A)',advance="no") 'VARIABLES = "X", "Y", "Z", "Y+", "WaveNo3", ' !5
        WRITE(DFLG,'(A)',advance="no") '"U11", "U22", "U33", "U12", "U13", "U23", '   !6
        WRITE(DFLG,'(A)',advance="no") '"O11", "O22", "O33", "O12", "O13", "O23", '   !6
        WRITE(DFLG,'(A)',advance="no") '"VO11", "VO12", "VO13", ' !3
        WRITE(DFLG,'(A)',advance="no") '"VO21", "VO22", "VO23", ' !3
        WRITE(DFLG,'(A)',advance="no") '"VO31", "VO32", "VO33"  ' !3
        WRITE(DFLG,'(A)',advance="no") '"E11", "E22", "E33", "E12", "E13", "E23", ' !6
        WRITE(DFLG,'(A)',advance="no") '"H11", "H22", "H33", "H12", "H13", "H23", ' !6
        WRITE(DFLG,'(A)',advance="no") '"EH11", "EH12", "EH13", ' !3
        WRITE(DFLG,'(A)',advance="no") '"EH21", "EH22", "EH23", ' !3
        WRITE(DFLG, '(A)'            ) '"EH31", "EH32", "EH33"  ' !3
        
        WRITE(DFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T= " ',ITERG, phyTIME, &
            ' ", I=', 1, ', J=',NND2,', K=',N3MH,', F=POINT'
        
        
        DO KC=1,N3MH
            AKE=( DBLE(KC-1)*2.0_wp*PI/HZ ) / Ret_ave_io
            DO JJ=1,NND2
                JJM=JGMV(JJ)
                JJC=JJ
                yplus=(1.0_wp-dabs(YND(JJ)))*REN*Utaw_ave_io
                
                IF(JJ==1)    JJM=1
                IF(JJ==NND2) JJC =NCL2
                WRITE(DFLG,'(47ES15.7)') 0.0_WP,YND(JJ),ZND(KC), yplus, AKE, &
                    0.5*( R11X3_xztLa(JJC,KC,M)+R11X3_xztLa (JJM,KC,M) ), &
                    0.5*( R22X3_xztLa(JJC,KC,M)+R22X3_xztLa(JJM,KC,M) ), &
                    0.5*( R33X3_xztLa(JJC,KC,M)+R33X3_xztLa(JJM,KC,M) ), &
                    0.5*( R12X3_xztLa(JJC,KC,M)+R12X3_xztLa(JJM,KC,M) ), &
                    0.5*( R13X3_xztLa(JJC,KC,M)+R13X3_xztLa(JJM,KC,M) ), &
                    0.5*( R23X3_xztLa(JJC,KC,M)+R23X3_xztLa(JJM,KC,M) ), &
                    0.5*( V11X3_xztLa(JJC,KC,M)+V11X3_xztLa(JJM,KC,M) ), &
                    0.5*( V22X3_xztLa(JJC,KC,M)+V22X3_xztLa(JJM,KC,M) ), &
                    0.5*( V33X3_xztLa(JJC,KC,M)+V33X3_xztLa(JJM,KC,M) ), &
                    0.5*( V12X3_xztLa(JJC,KC,M)+V12X3_xztLa(JJM,KC,M) ), &
                    0.5*( V13X3_xztLa(JJC,KC,M)+V13X3_xztLa(JJM,KC,M) ), &
                    0.5*( V23X3_xztLa(JJC,KC,M)+V23X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO11X3_xztLa(JJC,KC,M)+VO11X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO12X3_xztLa(JJC,KC,M)+VO12X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO13X3_xztLa(JJC,KC,M)+VO13X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO21X3_xztLa(JJC,KC,M)+VO21X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO22X3_xztLa(JJC,KC,M)+VO22X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO23X3_xztLa(JJC,KC,M)+VO23X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO31X3_xztLa(JJC,KC,M)+VO31X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO32X3_xztLa(JJC,KC,M)+VO32X3_xztLa(JJM,KC,M) ), &
                    0.5*( VO33X3_xztLa(JJC,KC,M)+VO33X3_xztLa(JJM,KC,M) ), &
                    0.5*( ENE11Z_xztLa(JJC,KC,M)+ENE11Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENE22Z_xztLa(JJC,KC,M)+ENE22Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENE33Z_xztLa(JJC,KC,M)+ENE33Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENE12Z_xztLa(JJC,KC,M)+ENE12Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENE13Z_xztLa(JJC,KC,M)+ENE13Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENE23Z_xztLa(JJC,KC,M)+ENE23Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENV11Z_xztLa(JJC,KC,M)+ENV11Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENV22Z_xztLa(JJC,KC,M)+ENV22Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENV33Z_xztLa(JJC,KC,M)+ENV33Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENV12Z_xztLa(JJC,KC,M)+ENV12Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENV13Z_xztLa(JJC,KC,M)+ENV13Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( ENV23Z_xztLa(JJC,KC,M)+ENV23Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO11Z_xztLa(JJC,KC,M)+EVO11Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO12Z_xztLa(JJC,KC,M)+EVO12Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO13Z_xztLa(JJC,KC,M)+EVO13Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO21Z_xztLa(JJC,KC,M)+EVO21Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO22Z_xztLa(JJC,KC,M)+EVO22Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO23Z_xztLa(JJC,KC,M)+EVO23Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO31Z_xztLa(JJC,KC,M)+EVO31Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO32Z_xztLa(JJC,KC,M)+EVO32Z_xztLa(JJM,KC,M) )*Ret_ave_io, &
                    0.5*( EVO33Z_xztLa(JJC,KC,M)+EVO33Z_xztLa(JJM,KC,M) )*Ret_ave_io

            END DO
        END DO
        CLOSE(DFLG)
        !===============plane x-y=====================
        FLNAME = TRIM(filepath4)//'Result.IO.Spectral.'//TRIM(STR)//'.Contours.yx.'//TRIM(PNTIM)//'.tec' 
        OPEN (DFLG, FILE=TRIM(ADJUSTL(FLNAME)))
        
        WRITE(DFLG,'(A)') 'TITLE = "DNS FLOW XY-plane"'
        WRITE(DFLG,'(A)',advance="no") 'VARIABLES = "X", "Y", "Z", "Y+", "WaveNo3", '
        WRITE(DFLG,'(A)',advance="no") '"U11", "U22", "U33", "U12", "U13", "U23", '
        WRITE(DFLG,'(A)',advance="no") '"O11", "O22", "O33", "O12", "O13", "O23", '
        WRITE(DFLG,'(A)',advance="no") '"VO11", "VO12", "VO13", '
        WRITE(DFLG,'(A)',advance="no") '"VO21", "VO22", "VO23", '
        WRITE(DFLG,'(A)',advance="no") '"VO31", "VO32", "VO33"  '
        WRITE(DFLG,'(A)',advance="no") '"E11", "E22", "E33", "E12", "E13", "E23", '
        WRITE(DFLG,'(A)',advance="no") '"H11", "H22", "H33", "H12", "H13", "H23", '
        WRITE(DFLG,'(A)',advance="no") '"EH11", "EH12", "EH13", '
        WRITE(DFLG,'(A)',advance="no") '"EH21", "EH22", "EH23", '
        WRITE(DFLG, '(A)'            ) '"EH31", "EH32", "EH33"  '
        
        WRITE(DFLG,'(A,1I11.1,1ES13.5,A,1I4.1,A,1I4.1,A,1I4.1,A)') 'ZONE T= " ',ITERG,phyTIME, &
            ' ", I=', N1MH, ', J=',NND2,', K=',1,', F=POINT'
        DO JJ=1,NND2
            yplus=(1.0_wp-dabs(YND(JJ)))*REN*Utaw_ave_io
            JJM=JGMV(JJ)
            JJC=JJ
            IF(JJ==1)    JJM=1
            IF(JJ==NND2) JJC =NCL2
            DO IC=1,N1MH
                AKE=( DBLE(IC-1)*2.0_wp*PI/HX_io ) / Ret_ave_io
                
                WRITE(DFLG,'(47ES15.7)') XND_io(IC),YND(JJ),0.0_wp, yplus, AKE, &
                    0.5*( R11X1_xztLa(JJC,IC,M)+R11X1_xztLa(JJM,IC,M) ), &
                    0.5*( R22X1_xztLa(JJC,IC,M)+R22X1_xztLa(JJM,IC,M) ), &
                    0.5*( R33X1_xztLa(JJC,IC,M)+R33X1_xztLa(JJM,IC,M) ), &
                    0.5*( R12X1_xztLa(JJC,IC,M)+R12X1_xztLa(JJM,IC,M) ), &
                    0.5*( R13X1_xztLa(JJC,IC,M)+R13X1_xztLa(JJM,IC,M) ), &
                    0.5*( R23X1_xztLa(JJC,IC,M)+R23X1_xztLa(JJM,IC,M) ), &
                    0.5*( V11X1_xztLa(JJC,IC,M)+V11X1_xztLa(JJM,IC,M) ), &
                    0.5*( V22X1_xztLa(JJC,IC,M)+V22X1_xztLa(JJM,IC,M) ), &
                    0.5*( V33X1_xztLa(JJC,IC,M)+V33X1_xztLa(JJM,IC,M) ), &
                    0.5*( V12X1_xztLa(JJC,IC,M)+V12X1_xztLa(JJM,IC,M) ), &
                    0.5*( V13X1_xztLa(JJC,IC,M)+V13X1_xztLa(JJM,IC,M) ), &
                    0.5*( V23X1_xztLa(JJC,IC,M)+V23X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO11X1_xztLa(JJC,IC,M)+VO11X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO12X1_xztLa(JJC,IC,M)+VO12X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO13X1_xztLa(JJC,IC,M)+VO13X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO21X1_xztLa(JJC,IC,M)+VO21X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO22X1_xztLa(JJC,IC,M)+VO22X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO23X1_xztLa(JJC,IC,M)+VO23X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO31X1_xztLa(JJC,IC,M)+VO31X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO32X1_xztLa(JJC,IC,M)+VO32X1_xztLa(JJM,IC,M) ), &
                    0.5*( VO33X1_xztLa(JJC,IC,M)+VO33X1_xztLa(JJM,IC,M) ), &
                    0.5*( ENE11T_xztLa(JJC,IC,M)+ENE11T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENE22T_xztLa(JJC,IC,M)+ENE22T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENE33T_xztLa(JJC,IC,M)+ENE33T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENE12T_xztLa(JJC,IC,M)+ENE12T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENE13T_xztLa(JJC,IC,M)+ENE13T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENE23T_xztLa(JJC,IC,M)+ENE23T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENV11T_xztLa(JJC,IC,M)+ENV11T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENV22T_xztLa(JJC,IC,M)+ENV22T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENV33T_xztLa(JJC,IC,M)+ENV33T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENV12T_xztLa(JJC,IC,M)+ENV12T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENV13T_xztLa(JJC,IC,M)+ENV13T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( ENV23T_xztLa(JJC,IC,M)+ENV23T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO11T_xztLa(JJC,IC,M)+EVO11T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO12T_xztLa(JJC,IC,M)+EVO12T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO13T_xztLa(JJC,IC,M)+EVO13T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO21T_xztLa(JJC,IC,M)+EVO21T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO22T_xztLa(JJC,IC,M)+EVO22T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO23T_xztLa(JJC,IC,M)+EVO23T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO31T_xztLa(JJC,IC,M)+EVO31T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO32T_xztLa(JJC,IC,M)+EVO32T_xztLa(JJM,IC,M) )*Ret_ave_io, &
                    0.5*( EVO33T_xztLa(JJC,IC,M)+EVO33T_xztLa(JJM,IC,M) )*Ret_ave_io
            ENDDO
        END DO
        
        CLOSE(DFLG)
            
        
        RETURN
    END SUBROUTINE
    !======================checking...==========================================================
!        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
!        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.Checking.Favre.'//TRIM(PNTIM)//'.tec'
!        OPEN (TECFLG_FavAG(1), FILE=TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1),'(A)') 'TITLE = " Favre Averged Flow (22 variables)" '
!        WRITE(TECFLG_FavAG(1),'(7A)') &
!            'VARIABLES = "1Y", "2Y+", "3Utau", "4Dwal",', &
!            '"5ViscStress_Tau_Mean_uu", "6ViscStress_Tau_Umea_uu", "7ViscStress_Tau_Uper_uu"', &
!            '"8ViscStress_Tau_Mean_uv",  "9ViscStress_Tau_Umea_uv", "10ViscStress_Tau_Uper_uv"', &
!            '"11ViscStress_Tau_Mean_uw", "12ViscStress_Tau_Umea_uw", "13ViscStress_Tau_Uper_uw',  &
!            '"14ViscStress_Tau_Mean_vv", "15ViscStress_Tau_Umea_vv", "16ViscStress_Tau_Uper_vv"', &
!            '"17ViscStress_Tau_Mean_vw", "18ViscStress_Tau_Umea_vw", "19ViscStress_Tau_Uper_vw"', &
!            '"20ViscStress_Tau_Mean_ww", "21ViscStress_Tau_Umea_ww", "22ViscStress_Tau_Uper_ww"'
!        WRITE(TECFLG_FavAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        
!        DO J=1,NCL2
!            WRITE(TECFLG_FavAG(1),'(22ES20.12)') YCC(J), YpluSD(J), UtauSD(J), DensSD(J) &
!                              Tau_Mean_RA(J,1,1), ViscStress_Tau_Umea(J,1,1), ViscStress_Tau_Uper(J,1,1), &
!                              Tau_Mean_RA(J,1,2), ViscStress_Tau_Umea(J,1,2), ViscStress_Tau_Uper(J,1,2), &
!                              Tau_Mean_RA(J,1,3), ViscStress_Tau_Umea(J,1,3), ViscStress_Tau_Uper(J,1,3), &
!                              Tau_Mean_RA(J,2,2), ViscStress_Tau_Umea(J,2,2), ViscStress_Tau_Uper(J,2,2), &
!                              Tau_Mean_RA(J,2,3), ViscStress_Tau_Umea(J,2,3), ViscStress_Tau_Uper(J,2,3), &
!                              Tau_Mean_RA(J,3,3), ViscStress_Tau_Umea(J,3,3), ViscStress_Tau_Uper(J,3,3)
            
!        END DO
!        CLOSE(TECFLG_FavAG(1))
        
        
!        !======================checking==Momentum Equation is xyz direction======================================================
!        BuoyForceTT = 0.0_WP
!        FCT = 0.0_WP
!        DO J=1, NCL2 
!            FCT(J,1) = -D1xztL_F0_io(J)*U_FA(J,1)*U_FA(J,2) /Utaw_ave_io/Utaw_ave_io/DenAvew + &
!                        Tau_Mean_RA(J,1,2)      *REN/Ret_ave_io/Utaw_ave_io/VisAvew- &
!                        uff2d_FA  (J,1,2)      /Utaw_ave_io/Utaw_ave_io/DenAvew
                        
!            FCT(J,2) = -D1xztL_F0_io(J)*U_FA(J,2)*U_FA(J,2)/Utaw_ave_io/Utaw_ave_io/DenAvew + &
!                        Tau_Mean_RA(J,2,2)     *REN/Ret_ave_io/Utaw_ave_io/VisAvew- &
!                        uff2d_FA  (J,2,2)     /Utaw_ave_io/Utaw_ave_io/DenAvew-&
!                        dPdX_RA(J,2)                      /Utaw_ave_io/Utaw_ave_io/D1xztL_F0_io(J)!
                        
!            FCT(J,3) = -D1xztL_F0_io(J)*U_FA(J,3)*U_FA(J,2)/Utaw_ave_io/Utaw_ave_io/DenAvew + &
!                        Tau_Mean_RA(J,2,3)     *REN/Ret_ave_io/Utaw_ave_io/VisAvew- &
!                        uff2d_FA  (J,2,3)      /Utaw_ave_io/Utaw_ave_io/DenAvew
                        
!            BuoyForceTT = BuoyForceTT + F_A*D1xztL_F0_io(J)/DYFI(J)
!        END DO
!        !===at wall surfaces====
!        FCT(0,1)    = Tauw_io(1)*REN/Ret_ave_io/Utaw_ave_io/VisAvew
!        FCT(0,2)    = 0.0_WP
!        FCT(0,3)    = 0.0_WP
        
!        FCT(NND2,1) = Tauw_io(2)*REN/Ret_ave_io/Utaw_ave_io/VisAvew
!        FCT(NND2,2) = 0.0_WP
!        FCT(NND2,3) = 0.0_WP
        
!        !===gradient at cell centre====
!        NSFbal_FA = 0.0_WP
!        DO J=2, NCL2-1
!            DO N=1,NDV
!                NSFbal_FA(J,N)= ( ( YCL2ND_WFB(J+1)*FCT(J  ,N) + YCL2ND_WFF(J+1)*FCT(J+1,N) ) -             &
!                                  ( YCL2ND_WFF(J)  *FCT(J  ,N) + YCL2ND_WFB(J)  *FCT(J-1,N) ) ) * DYFI(J) + &
!                                F_A/Utaw_ave_io/Utaw_ave_io*D1xztL_F0_io(J)/DenAvew*DBLE(IBuoF(N)) - &
!                                F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
!!                !WRITE(*,*) N, J, ( ( YCL2ND_WFB(J+1)*FCT(J  ,N) + YCL2ND_WFF(J+1)*FCT(J+1,N) ) -             &
!!                                  ( YCL2ND_WFF(J)  *FCT(J  ,N) + YCL2ND_WFB(J)  *FCT(J-1,N) ) ) * DYFI(J), &
!!                                  F_A/Utaw_ave_io/Utaw_ave_io*D1xztL_F0_io(J)/DenAvew*DBLE(IBuoF(N)), &
!!                                  F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
!            END DO
!        END DO
!        DO N=1,NDV
!            NSFbal_FA(1,   N)= ( ( YCL2ND_WFB(1+1)*FCT(1  ,N) +                                           &
!                                   YCL2ND_WFF(1+1)*FCT(1+1,N) ) - FCT(0   ,N) ) * DYFI(1) +               &
!                                 F_A/Utaw_ave_io/Utaw_ave_io*D1xztL_F0_io(1)/DenAvew*DBLE(IBuoF(N)) - &
!                                 F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
                                 
!            NSFbal_FA(NCL2,N)= ( FCT(NND2,N) - ( YCL2ND_WFF(NCL2)  *FCT(NCL2  ,N) +                       &
!                                                 YCL2ND_WFB(NCL2)  *FCT(NCL2-1,N) ) ) * DYFI(NCL2) +      &
!                                F_A/Utaw_ave_io/Utaw_ave_io*D1xztL_F0_io(NCL2)/DenAvew*DBLE(IBuoF(N)) - &
!                                F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
!        END DO
        
!        !===integral along y direction===============================
!        NSFbalt_FA = 0.0_WP
!        DO J=1, NCL2
!            DO N=1, NDV
!                NSFbalt_FA(N)   = NSFbalt_FA(N) + NSFbal_FA(J,N)/DYFI(J)
!            END DO
!            !!WRITE(*,*) J, NSFbalt_FA(1:3)
!        END DO
        
        
!        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
!        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.Checking.NSForceBalance.FA.'//TRIM(PNTIM)//'.tec'
!        OPEN (TECFLG_FavAG(1), FILE=TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1),'(A)') 'TITLE = " Favre Averged Flow (21 variables)" '
!        WRITE(TECFLG_FavAG(1),'(6A)') &
!            'VARIABLES = "1Y", "2Y+", "3Utau", "4FC_x", "5FC_y", "6FC_z","7FCGRD_x", "8FCGRD_y", "9FCGRD_z"'
!        WRITE(TECFLG_FavAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        
!        Dintg = 0.0_WP
!        DO J=1,NCL2
            
!            !if(YND(J).LT.0.0_wp) then
!            IF(J.LT.J4SS0) THEN
!                COE = Utaw_io(1)
!            else
!                COE = Utaw_io(2)
!            end if
!            IF(J==1) THEN
!                DENtemp =0.5_WP*( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J) + YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) +             &
!                                  Dwal(1) )
!            ELSE IF(J==NCL2) THEN
!                DENtemp =0.5_WP*( Dwal(2) +             &
!                                  ( YCL2ND_WFF(J)  *D1xztL_F0_io(J) + YCL2ND_WFB(J)  *D1xztL_F0_io(J-1) ) )
!            ELSE
!                DENtemp =0.5_WP*( ( YCL2ND_WFB(J+1)*D1xztL_F0_io(J) + YCL2ND_WFF(J+1)*D1xztL_F0_io(J+1) ) +             &
!                                  ( YCL2ND_WFF(J)  *D1xztL_F0_io(J) + YCL2ND_WFB(J)  *D1xztL_F0_io(J-1) ) )
!            END IF
!            Dintg = Dintg + (DENtemp/DenAvew-1.0_wp)/DYFI(J)
!            WRITE(TECFLG_FavAG(1),'(13ES20.12)') YCC(J), (1.0_wp-dabs(YCC(J)))*REN*COE, COE, &
!                              FCT(J,1:3), NSFbal_FA(J,1:3), &
!                              Tau_Mean_RA(J,1,2)     *REN/Ret_ave_io/Utaw_ave_io, &
!                             -uff2d_FA  (J,1,2)     /Utaw_ave_io/Utaw_ave_io/DenAvew, &
!                              F_A/Utaw_ave_io/Utaw_ave_io*Dintg,&
!                             -F_A/Utaw_ave_io/Utaw_ave_io/DYFI(J)
!        END DO
!        CLOSE(TECFLG_FavAG(1))



!!======================checking==Momentum Equation is xyz direction======================================================
!        FCT = 0.0_WP
!        DO J=1, NCL2 

!            FCT(J,1) = -D1xztL_F0_io(J)*U1xztL_F0_io(J,1)*U1xztL_F0_io(J,2)/Utaw_ave_io/Utaw_ave_io/DenAvew+ &
!                        Tau_Mean_RA(J,1,2)*REN/Ret_ave_io/Utaw_ave_io/VisAvew-&
!                        uf2d_RA  (J,1,2)/Utaw_ave_io/Utaw_ave_io/DenAvew
                        
!            FCT(J,2) = -D1xztL_F0_io(J)*U1xztL_F0_io(J,2)*U1xztL_F0_io(J,2)/Utaw_ave_io/Utaw_ave_io/DenAvew+ &
!                        Tau_Mean_RA(J,2,2)*REN/Ret_ave_io/Utaw_ave_io/VisAvew-&
!                        uf2d_RA  (J,2,2)/Utaw_ave_io/Utaw_ave_io/DenAvew-&
!                        dPdX_RA(J,2)/Utaw_ave_io/Utaw_ave_io/D1xztL_F0_io(J)!
                        
!            FCT(J,3) = -D1xztL_F0_io(J)*U1xztL_F0_io(J,3)*U1xztL_F0_io(J,2)/Utaw_ave_io/Utaw_ave_io/DenAvew+ &
!                        Tau_Mean_RA(J,2,3)*REN/Ret_ave_io/Utaw_ave_io/VisAvew-&
!                        uf2d_RA  (J,2,3)/Utaw_ave_io/Utaw_ave_io/DenAvew
!        END DO
!        !===at wall surfaces====
!        FCT(0,1)    = Tauw_io(1)*REN/Ret_ave_io/Utaw_ave_io/VisAvew
!        FCT(0,2)    = 0.0_WP
!        FCT(0,3)    = 0.0_WP
        
!        FCT(NND2,1) = Tauw_io(2)*REN/Ret_ave_io/Utaw_ave_io/VisAvew
!        FCT(NND2,2) = 0.0_WP
!        FCT(NND2,3) = 0.0_WP
        
!        !===gradient at cell centre====
!        NSFbal_RA = 0.0_WP
!        DO J=2, NCL2-1
!            DO N=1,NDV
!                NSFbal_RA(J,N)= ( ( YCL2ND_WFB(J+1)*FCT(J  ,N) + YCL2ND_WFF(J+1)*FCT(J+1,N) ) -             &
!                                  ( YCL2ND_WFF(J)  *FCT(J  ,N) + YCL2ND_WFB(J)  *FCT(J-1,N) ) ) * DYFI(J) + &
!                                F_A*D1xztL_F0_io(J)/DenAvew/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N)) - &
!                                F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
!                !!WRITE(*,*) N, J, ( ( YCL2ND_WFB(J+1)*FCT(J  ,N) + YCL2ND_WFF(J+1)*FCT(J+1,N) ) -             &
!                !                  ( YCL2ND_WFF(J)  *FCT(J  ,N) + YCL2ND_WFB(J)  *FCT(J-1,N) ) ) * DYFI(J), &
!                !                  F_A*D1xztL_F0_io(J)/DenAvew/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N)), &
!                !                  F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
!            END DO
!        END DO
!        DO N=1,NDV
!            NSFbal_RA(1,   N)= ( ( YCL2ND_WFB(1+1)*FCT(1  ,N) +                                           &
!                                   YCL2ND_WFF(1+1)*FCT(1+1,N) ) - FCT(0   ,N) ) * DYFI(1) +               &
!                                 F_A*D1xztL_F0_io(1)/DenAvew/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N)) - &
!                                 F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
                                 
!            NSFbal_RA(NCL2,N)= ( FCT(NND2,N) - ( YCL2ND_WFF(NCL2)  *FCT(NCL2  ,N) +                       &
!                                                 YCL2ND_WFB(NCL2)  *FCT(NCL2-1,N) ) ) * DYFI(NCL2) +      &
!                                F_A*D1xztL_F0_io(NCL2)/DenAvew/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N)) - &
!                                F_A/Utaw_ave_io/Utaw_ave_io*DBLE(IBuoF(N))
!        END DO
        
!        !===integral along y direction===============================
!        NSFbalt_RA = 0.0_WP
!        DO J=1, NCL2
!            DO N=1, NDV
!                NSFbalt_RA(N)   = NSFbalt_RA(N) + NSFbal_RA(J,N)/RCCI1(J)/DYFI(J)
!            END DO
!            !!WRITE(*,*) J, NSFbalt_RA(1:3)
!        END DO
!!======================checking==Momentum Equation is xyz direction======================================================

!!====================Checking=================================
!        FCT = 0.0_WP
!        DO J=1, NCL2 
!            FCT(J) = -D1xztL_F0_io(J)*U_FA(J,2)*H_FA(J) + &
!                      DTDLKxztL_F0_io(J,2)*CTHECD - &
!                      uffhffd_FA  (J,2)
!        END DO
!        !===at wall surfaces====
!        FCT(0)    = qw(1)
!        FCT(NND2) = qw(2)
        
!        !===gradient at cell centre====
!        ENEbal_FA = 0.0_WP
!        DO J=2, NCL2-1
            
!            ENEbal_FA(J)= ( ( YCL2ND_WFB(J+1)*FCT(J) + YCL2ND_WFF(J+1)*FCT(J+1) ) - &
!                            ( YCL2ND_WFF(J)  *FCT(J) + YCL2ND_WFB(J)  *FCT(J-1) ) ) * DYFI(J)*(-1.0_WP)

!        END DO
!        ENEbal_FA(1)   = ( ( YCL2ND_WFB(1+1) *FCT(1)    + YCL2ND_WFF(1+1) *FCT(1+1)    )- FCT(0)   ) * DYFI(1)*(-1.0_WP)
!        ENEbal_FA(NCL2)= ( ( YCL2ND_WFF(NCL2)*FCT(NCL2) + YCL2ND_WFB(NCL2)*FCT(NCL2-1) )- FCT(NND2)) * DYFI(NCL2)*((-1.0_WP)**2)
                             
!        !===integral along y direction===============================
!        ENEbalt_FA = 0.0_WP
!        DO J=1, NCL2
!            ENEbalt_FA   = ENEbalt_FA + ENEbal_FA(J)/RCCI1(J)/DYFI(J)
!            !!WRITE(*,*) J, NSFbalt_FA(1:3)
!        END DO
        
!!====================calcuate variables for RANS============================================
!        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
!        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.ForRANS.Heat.Transfer.'//TRIM(PNTIM)//'.tec'
!        OPEN (TECFLG_FavAG(1), FILE=TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1),'(A)') 'TITLE = " Favre Averged Flow (9 variables)" '
!        WRITE(TECFLG_FavAG(1),'(A,A)') &
!            'VARIABLES = "1Y", "2Y+", "3Utau", "4MutOverPrt(x)", "5MutOverPrt(y)", "6MutOverPrt(z)", ', &
!            ' "7Prt(x)", "8Prt(y)", "9Prt(z)"'
!        WRITE(TECFLG_FavAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        
!        DO J=1,NCL2

!            MutOverPrt(1:3) = uffhffd_FA(J,1:3)/(dHdX_FA(J,1:3)+REALMIN)
            
!            Pruv(1:3) = RANS_Mut(J,1,2)/MutOverPrt(1:3)

!            WRITE(TECFLG_FavAG(1),'(9ES20.12)') YCC(J), YpluSD(J), UtauSD(J), &
!                              MutOverPrt(1:3), Pruv(1:3)
!        END DO
!        CLOSE(TECFLG_FavAG(1))
        
        
!!======================checking...==========================================================
!        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
!        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.Checking.Heat.Transfer.'//TRIM(PNTIM)//'.tec'
!        OPEN (TECFLG_FavAG(1), FILE=TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1),'(A)') 'TITLE = " Favre Averged Flow (21 variables)" '
!        WRITE(TECFLG_FavAG(1),'(A)') &
!            'VARIABLES = "1Y", "2Y+", "3Utau", "4Cp(<T>)", "5Cp(dh/dT)"'
!        WRITE(TECFLG_FavAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        
!        DO J=1,NCL2
!            Cp_eval = dHdX_FA(J,2)/dTdX(J,2)

!            WRITE(TECFLG_FavAG(1),'(5ES20.12)') YCC(J), YpluSD(J), UtauSD(J), &
!                              CpT(J), Cp_eval
            
!        END DO
!        CLOSE(TECFLG_FavAG(1))
!======================checking...==========================================================
        
        
        
!        WRITE(PNTIM,'(1ES15.9)') phyTIME_io
!        FLNM= TRIM(filepath4)//'Result.IO.'//TRIM(STDIM(1))//'.Profile.Checking.EnergyBalance.FA.'//TRIM(PNTIM)//'.tec'
!        OPEN (TECFLG_FavAG(1), FILE=TRIM(ADJUSTL(FLNM)))
!        WRITE(TECFLG_FavAG(1),'(A)') 'TITLE = " Favre Averged Flow (4 variables)" '
!        WRITE(TECFLG_FavAG(1),'(6A)') &
!            'VARIABLES = "1Y", "2Y+", "3Utau", "4ENEbal"'
!        WRITE(TECFLG_FavAG(1),'(A)') 'ZONE T="'//TRIM(ADJUSTL(teczonename))//' " '
        
!        DO J=1,NCL2
!            WRITE(TECFLG_FavAG(1),'(4ES20.12)') YCC(J), YpluSD(J), UtauSD(J), ENEbal_FA(J)
!        END DO
!        CLOSE(TECFLG_FavAG(1))
        
