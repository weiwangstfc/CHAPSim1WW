    MODULE SPECO_info
        use cparam 
        
        INTEGER(4) :: N1MH
        INTEGER(4) :: N3MH
        INTEGER(4) :: N1MD
        INTEGER(4) :: N3MD
    
        REAL(WP), ALLOCATABLE :: RHS   (:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENEJTF(:, :)
        REAL(WP), ALLOCATABLE :: ENEJME(:, :)
        
        REAL(WP), ALLOCATABLE :: ENE1MA(:)
        REAL(WP), ALLOCATABLE :: ENE3MA(:)
        
        REAL(WP), ALLOCATABLE :: EN1IK (:, :, :)
        REAL(WP), ALLOCATABLE :: EN3KI (:, :, :)
        
        REAL(WP), ALLOCATABLE :: ENE1(:, :)
        REAL(WP), ALLOCATABLE :: ENE3(:, :)
        
        REAL(WP), ALLOCATABLE :: CORX1(:, :)
        REAL(WP), ALLOCATABLE :: CORX3(:, :)
        
        !================Velocity =====================
        REAL(WP), ALLOCATABLE :: R11X1(:, :)
        REAL(WP), ALLOCATABLE :: R22X1(:, :)
        REAL(WP), ALLOCATABLE :: R33X1(:, :)
        REAL(WP), ALLOCATABLE :: R12X1(:, :)
        REAL(WP), ALLOCATABLE :: R13X1(:, :)
        REAL(WP), ALLOCATABLE :: R23X1(:, :)
        
        REAL(WP), ALLOCATABLE :: R11X3(:, :)
        REAL(WP), ALLOCATABLE :: R22X3(:, :)
        REAL(WP), ALLOCATABLE :: R33X3(:, :)
        REAL(WP), ALLOCATABLE :: R12X3(:, :)
        REAL(WP), ALLOCATABLE :: R13X3(:, :)
        REAL(WP), ALLOCATABLE :: R23X3(:, :)
        
        REAL(WP), ALLOCATABLE :: ENE11T(:, :)
        REAL(WP), ALLOCATABLE :: ENE22T(:, :)
        REAL(WP), ALLOCATABLE :: ENE33T(:, :)
        REAL(WP), ALLOCATABLE :: ENE12T(:, :)
        REAL(WP), ALLOCATABLE :: ENE13T(:, :)
        REAL(WP), ALLOCATABLE :: ENE23T(:, :)
        
        REAL(WP), ALLOCATABLE :: ENE11Z(:, :)
        REAL(WP), ALLOCATABLE :: ENE22Z(:, :)
        REAL(WP), ALLOCATABLE :: ENE33Z(:, :)
        REAL(WP), ALLOCATABLE :: ENE12Z(:, :)
        REAL(WP), ALLOCATABLE :: ENE13Z(:, :)
        REAL(WP), ALLOCATABLE :: ENE23Z(:, :)
        
        !================Voriticity ====================
        REAL(WP), ALLOCATABLE :: V11X1(:, :)
        REAL(WP), ALLOCATABLE :: V22X1(:, :)
        REAL(WP), ALLOCATABLE :: V33X1(:, :)
        REAL(WP), ALLOCATABLE :: V12X1(:, :)
        REAL(WP), ALLOCATABLE :: V13X1(:, :)
        REAL(WP), ALLOCATABLE :: V23X1(:, :)
        
        REAL(WP), ALLOCATABLE :: V11X3(:, :)
        REAL(WP), ALLOCATABLE :: V22X3(:, :)
        REAL(WP), ALLOCATABLE :: V33X3(:, :)
        REAL(WP), ALLOCATABLE :: V12X3(:, :)
        REAL(WP), ALLOCATABLE :: V13X3(:, :)
        REAL(WP), ALLOCATABLE :: V23X3(:, :)
        
        REAL(WP), ALLOCATABLE :: ENV11T(:, :)
        REAL(WP), ALLOCATABLE :: ENV22T(:, :)
        REAL(WP), ALLOCATABLE :: ENV33T(:, :)
        REAL(WP), ALLOCATABLE :: ENV12T(:, :)
        REAL(WP), ALLOCATABLE :: ENV13T(:, :)
        REAL(WP), ALLOCATABLE :: ENV23T(:, :)
        
        REAL(WP), ALLOCATABLE :: ENV11Z(:, :)
        REAL(WP), ALLOCATABLE :: ENV22Z(:, :)
        REAL(WP), ALLOCATABLE :: ENV33Z(:, :)
        REAL(WP), ALLOCATABLE :: ENV12Z(:, :)
        REAL(WP), ALLOCATABLE :: ENV13Z(:, :)
        REAL(WP), ALLOCATABLE :: ENV23Z(:, :)
        
        
        !===============Voriticity & Velocity=========================
        REAL(WP), ALLOCATABLE :: VO11X1(:, :)
        REAL(WP), ALLOCATABLE :: VO12X1(:, :)
        REAL(WP), ALLOCATABLE :: VO13X1(:, :)
        
        REAL(WP), ALLOCATABLE :: VO21X1(:, :)
        REAL(WP), ALLOCATABLE :: VO22X1(:, :)
        REAL(WP), ALLOCATABLE :: VO23X1(:, :)
        
        REAL(WP), ALLOCATABLE :: VO31X1(:, :)
        REAL(WP), ALLOCATABLE :: VO32X1(:, :)
        REAL(WP), ALLOCATABLE :: VO33X1(:, :)
        
        REAL(WP), ALLOCATABLE :: VO11X3(:, :)
        REAL(WP), ALLOCATABLE :: VO12X3(:, :)
        REAL(WP), ALLOCATABLE :: VO13X3(:, :)
        
        REAL(WP), ALLOCATABLE :: VO21X3(:, :)
        REAL(WP), ALLOCATABLE :: VO22X3(:, :)
        REAL(WP), ALLOCATABLE :: VO23X3(:, :)
        
        REAL(WP), ALLOCATABLE :: VO31X3(:, :)
        REAL(WP), ALLOCATABLE :: VO32X3(:, :)
        REAL(WP), ALLOCATABLE :: VO33X3(:, :)
        
        REAL(WP), ALLOCATABLE :: EVO11T(:, :)
        REAL(WP), ALLOCATABLE :: EVO12T(:, :)
        REAL(WP), ALLOCATABLE :: EVO13T(:, :)
        
        REAL(WP), ALLOCATABLE :: EVO21T(:, :)
        REAL(WP), ALLOCATABLE :: EVO22T(:, :)
        REAL(WP), ALLOCATABLE :: EVO23T(:, :)
        
        REAL(WP), ALLOCATABLE :: EVO31T(:, :)
        REAL(WP), ALLOCATABLE :: EVO32T(:, :)
        REAL(WP), ALLOCATABLE :: EVO33T(:, :)
        
        REAL(WP), ALLOCATABLE :: EVO11Z(:, :)
        REAL(WP), ALLOCATABLE :: EVO12Z(:, :)
        REAL(WP), ALLOCATABLE :: EVO13Z(:, :)
        
        REAL(WP), ALLOCATABLE :: EVO21Z(:, :)
        REAL(WP), ALLOCATABLE :: EVO22Z(:, :)
        REAL(WP), ALLOCATABLE :: EVO23Z(:, :)
        
        REAL(WP), ALLOCATABLE :: EVO31Z(:, :)
        REAL(WP), ALLOCATABLE :: EVO32Z(:, :)
        REAL(WP), ALLOCATABLE :: EVO33Z(:, :)
        
        !============================================
        
        REAL(WP), ALLOCATABLE :: EE1_11(:, :)
        REAL(WP), ALLOCATABLE :: EE1_22(:, :)
        REAL(WP), ALLOCATABLE :: EE1_33(:, :)
        REAL(WP), ALLOCATABLE :: EE1_12(:, :)
        REAL(WP), ALLOCATABLE :: EE1_13(:, :)
        REAL(WP), ALLOCATABLE :: EE1_23(:, :)
        
        REAL(WP), ALLOCATABLE :: UU1_11(:, :)
        REAL(WP), ALLOCATABLE :: UU1_22(:, :)
        REAL(WP), ALLOCATABLE :: UU1_33(:, :)
        REAL(WP), ALLOCATABLE :: UU1_12(:, :)
        REAL(WP), ALLOCATABLE :: UU1_13(:, :)
        REAL(WP), ALLOCATABLE :: UU1_23(:, :)

        REAL(WP), ALLOCATABLE :: HH1_11(:, :)
        REAL(WP), ALLOCATABLE :: HH1_22(:, :)
        REAL(WP), ALLOCATABLE :: HH1_33(:, :)
        REAL(WP), ALLOCATABLE :: HH1_12(:, :)
        REAL(WP), ALLOCATABLE :: HH1_13(:, :)
        REAL(WP), ALLOCATABLE :: HH1_23(:, :)

        REAL(WP), ALLOCATABLE :: VV1_11(:, :)
        REAL(WP), ALLOCATABLE :: VV1_22(:, :)
        REAL(WP), ALLOCATABLE :: VV1_33(:, :)
        REAL(WP), ALLOCATABLE :: VV1_12(:, :)
        REAL(WP), ALLOCATABLE :: VV1_13(:, :)
        REAL(WP), ALLOCATABLE :: VV1_23(:, :)

        REAL(WP), ALLOCATABLE :: EH1_11(:, :)
        REAL(WP), ALLOCATABLE :: EH1_12(:, :)
        REAL(WP), ALLOCATABLE :: EH1_13(:, :)
        
        REAL(WP), ALLOCATABLE :: EH1_21(:, :)
        REAL(WP), ALLOCATABLE :: EH1_22(:, :)
        REAL(WP), ALLOCATABLE :: EH1_23(:, :)

        REAL(WP), ALLOCATABLE :: EH1_31(:, :)
        REAL(WP), ALLOCATABLE :: EH1_32(:, :)
        REAL(WP), ALLOCATABLE :: EH1_33(:, :)

        REAL(WP), ALLOCATABLE :: UV1_11(:, :)
        REAL(WP), ALLOCATABLE :: UV1_12(:, :)
        REAL(WP), ALLOCATABLE :: UV1_13(:, :)

        REAL(WP), ALLOCATABLE :: UV1_21(:, :)
        REAL(WP), ALLOCATABLE :: UV1_22(:, :)
        REAL(WP), ALLOCATABLE :: UV1_23(:, :)

        REAL(WP), ALLOCATABLE :: UV1_31(:, :)
        REAL(WP), ALLOCATABLE :: UV1_32(:, :)
        REAL(WP), ALLOCATABLE :: UV1_33(:, :)
        
        REAL(WP), ALLOCATABLE :: WORK_EE1_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE1_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE1_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE1_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE1_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE1_23(:, :)
        
        REAL(WP), ALLOCATABLE :: WORK_UU1_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU1_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU1_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU1_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU1_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU1_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_HH1_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH1_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH1_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH1_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH1_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH1_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_VV1_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV1_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV1_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV1_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV1_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV1_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_EH1_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH1_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH1_13(:, :)
        
        REAL(WP), ALLOCATABLE :: WORK_EH1_21(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH1_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH1_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_EH1_31(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH1_32(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH1_33(:, :)

        REAL(WP), ALLOCATABLE :: WORK_UV1_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV1_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV1_13(:, :)

        REAL(WP), ALLOCATABLE :: WORK_UV1_21(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV1_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV1_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_UV1_31(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV1_32(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV1_33(:, :)


        REAL(WP), ALLOCATABLE :: EE3_11(:, :)
        REAL(WP), ALLOCATABLE :: EE3_22(:, :)
        REAL(WP), ALLOCATABLE :: EE3_33(:, :)
        REAL(WP), ALLOCATABLE :: EE3_12(:, :)
        REAL(WP), ALLOCATABLE :: EE3_13(:, :)
        REAL(WP), ALLOCATABLE :: EE3_23(:, :)
        
        REAL(WP), ALLOCATABLE :: UU3_11(:, :)
        REAL(WP), ALLOCATABLE :: UU3_22(:, :)
        REAL(WP), ALLOCATABLE :: UU3_33(:, :)
        REAL(WP), ALLOCATABLE :: UU3_12(:, :)
        REAL(WP), ALLOCATABLE :: UU3_13(:, :)
        REAL(WP), ALLOCATABLE :: UU3_23(:, :)

        REAL(WP), ALLOCATABLE :: HH3_11(:, :)
        REAL(WP), ALLOCATABLE :: HH3_22(:, :)
        REAL(WP), ALLOCATABLE :: HH3_33(:, :)
        REAL(WP), ALLOCATABLE :: HH3_12(:, :)
        REAL(WP), ALLOCATABLE :: HH3_13(:, :)
        REAL(WP), ALLOCATABLE :: HH3_23(:, :)

        REAL(WP), ALLOCATABLE :: VV3_11(:, :)
        REAL(WP), ALLOCATABLE :: VV3_22(:, :)
        REAL(WP), ALLOCATABLE :: VV3_33(:, :)
        REAL(WP), ALLOCATABLE :: VV3_12(:, :)
        REAL(WP), ALLOCATABLE :: VV3_13(:, :)
        REAL(WP), ALLOCATABLE :: VV3_23(:, :)

        REAL(WP), ALLOCATABLE :: EH3_11(:, :)
        REAL(WP), ALLOCATABLE :: EH3_12(:, :)
        REAL(WP), ALLOCATABLE :: EH3_13(:, :)
        
        REAL(WP), ALLOCATABLE :: EH3_21(:, :)
        REAL(WP), ALLOCATABLE :: EH3_22(:, :)
        REAL(WP), ALLOCATABLE :: EH3_23(:, :)

        REAL(WP), ALLOCATABLE :: EH3_31(:, :)
        REAL(WP), ALLOCATABLE :: EH3_32(:, :)
        REAL(WP), ALLOCATABLE :: EH3_33(:, :)

        REAL(WP), ALLOCATABLE :: UV3_11(:, :)
        REAL(WP), ALLOCATABLE :: UV3_12(:, :)
        REAL(WP), ALLOCATABLE :: UV3_13(:, :)

        REAL(WP), ALLOCATABLE :: UV3_21(:, :)
        REAL(WP), ALLOCATABLE :: UV3_22(:, :)
        REAL(WP), ALLOCATABLE :: UV3_23(:, :)

        REAL(WP), ALLOCATABLE :: UV3_31(:, :)
        REAL(WP), ALLOCATABLE :: UV3_32(:, :)
        REAL(WP), ALLOCATABLE :: UV3_33(:, :)
        
        REAL(WP), ALLOCATABLE :: WORK_EE3_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE3_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE3_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE3_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE3_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EE3_23(:, :)
        
        REAL(WP), ALLOCATABLE :: WORK_UU3_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU3_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU3_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU3_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU3_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UU3_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_HH3_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH3_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH3_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH3_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH3_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_HH3_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_VV3_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV3_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV3_33(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV3_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV3_13(:, :)
        REAL(WP), ALLOCATABLE :: WORK_VV3_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_EH3_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH3_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH3_13(:, :)
        
        REAL(WP), ALLOCATABLE :: WORK_EH3_21(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH3_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH3_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_EH3_31(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH3_32(:, :)
        REAL(WP), ALLOCATABLE :: WORK_EH3_33(:, :)

        REAL(WP), ALLOCATABLE :: WORK_UV3_11(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV3_12(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV3_13(:, :)

        REAL(WP), ALLOCATABLE :: WORK_UV3_21(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV3_22(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV3_23(:, :)

        REAL(WP), ALLOCATABLE :: WORK_UV3_31(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV3_32(:, :)
        REAL(WP), ALLOCATABLE :: WORK_UV3_33(:, :)
        
    END MODULE

    SUBROUTINE SPECO_ALLOCATE
        use mesh_info
        use SPECO_info
        IMPLICIT NONE
    
         !===========u'=u-<u>=======================================
        ALLOCATE ( RHS(NCL1_io,N2DO(0),NCL3) ) ; RHS = 0.0_WP 
        
        !===========energy in each direction ( not used )======================
        ! IN PHYSICAL SPACE
        ALLOCATE ( ENEJTF(6,N2DO(0)) ) ; ENEJTF=0.0_WP  ! <u'u'>_zx
        ALLOCATE ( ENEJME(6,N2DO(0)) ) ; ENEJME=0.0_WP  ! <u'>_zx  ? not zero?
        
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( ENE1MA(6        ) ) ; ENE1MA=0.0_WP  ! 
        ALLOCATE ( ENE3MA(6        ) ) ; ENE3MA=0.0_WP  ! 
        
        !==========================
        ! IN WAVE NUMBER SPACE
        !EN1IK = STREAMWISE INVERSE FFT X1 of u'_ndv
        !EN3KI = SPANWISE   INVERSE FFT X3 of u'_ndv
        ALLOCATE ( EN1IK(6,N1MD,NCL3   )) ; EN1IK = 0.0_WP
        ALLOCATE ( EN3KI(6,N3MD,NCL1_io)) ; EN3KI = 0.0_WP
        
        !=================================================
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( ENE1(NCL2,NCL1_io)) ; ENE1(:,:)=0.0_WP !E(k_x)
        ALLOCATE ( ENE3(NCL2,NCL3   )) ; ENE3(:,:)=0.0_WP !E(k_z)
        
        !=================================================
        ! IN PHYSICAL SPACE
        ALLOCATE ( CORX1(NCL2,NCL1_io)) ; CORX1(:,:)=0.0_WP 
        ALLOCATE ( CORX3(NCL2,NCL3   )) ; CORX3(:,:)=0.0_WP
        
        !========CORRELATION===========================
        ! IN PHYSICAL SPACE
        ! CORRELATIONo R(IJ)X(K) = Velocity_i and Velocity_J along the direction K
        ALLOCATE ( R11X1(NCL2,NCL1_io) ) ; R11X1(:,:) = 0.0_WP 
        ALLOCATE ( R22X1(NCL2,NCL1_io) ) ; R22X1(:,:) = 0.0_WP 
        ALLOCATE ( R33X1(NCL2,NCL1_io) ) ; R33X1(:,:) = 0.0_WP 
        ALLOCATE ( R12X1(NCL2,NCL1_io) ) ; R12X1(:,:) = 0.0_WP 
        ALLOCATE ( R13X1(NCL2,NCL1_io) ) ; R13X1(:,:) = 0.0_WP 
        ALLOCATE ( R23X1(NCL2,NCL1_io) ) ; R23X1(:,:) = 0.0_WP 
        
        ALLOCATE ( R11X3(NCL2,NCL3) ) ; R11X3(:,:) = 0.0_WP 
        ALLOCATE ( R22X3(NCL2,NCL3) ) ; R22X3(:,:) = 0.0_WP 
        ALLOCATE ( R33X3(NCL2,NCL3) ) ; R33X3(:,:) = 0.0_WP 
        ALLOCATE ( R12X3(NCL2,NCL3) ) ; R12X3(:,:) = 0.0_WP 
        ALLOCATE ( R13X3(NCL2,NCL3) ) ; R13X3(:,:) = 0.0_WP 
        ALLOCATE ( R23X3(NCL2,NCL3) ) ; R23X3(:,:) = 0.0_WP 
        
        !===============
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( ENE11T(NCL2,NCL1_io) ) ; ENE11T(:,:)=0.0_WP 
        ALLOCATE ( ENE22T(NCL2,NCL1_io) ) ; ENE22T(:,:)=0.0_WP 
        ALLOCATE ( ENE33T(NCL2,NCL1_io) ) ; ENE33T(:,:)=0.0_WP 
        ALLOCATE ( ENE12T(NCL2,NCL1_io) ) ; ENE12T(:,:)=0.0_WP 
        ALLOCATE ( ENE13T(NCL2,NCL1_io) ) ; ENE13T(:,:)=0.0_WP 
        ALLOCATE ( ENE23T(NCL2,NCL1_io) ) ; ENE23T(:,:)=0.0_WP 
        
        !=============
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( ENE11Z(NCL2,NCL3) ) ; ENE11Z(:,:)=0.0_WP 
        ALLOCATE ( ENE22Z(NCL2,NCL3) ) ; ENE22Z(:,:)=0.0_WP 
        ALLOCATE ( ENE33Z(NCL2,NCL3) ) ; ENE33Z(:,:)=0.0_WP  
        ALLOCATE ( ENE12Z(NCL2,NCL3) ) ; ENE12Z(:,:)=0.0_WP 
        ALLOCATE ( ENE13Z(NCL2,NCL3) ) ; ENE13Z(:,:)=0.0_WP 
        ALLOCATE ( ENE23Z(NCL2,NCL3) ) ; ENE23Z(:,:)=0.0_WP 

        !========CORRELATION===========================
        ! IN PHYSICAL SPACE
        ! CORRELATIONo R(IJ)X(K) = voriticity_i and voriticity_J along the direction K
        ALLOCATE ( V11X1(NCL2,NCL1_io) ) ; V11X1(:,:) = 0.0_WP 
        ALLOCATE ( V22X1(NCL2,NCL1_io) ) ; V22X1(:,:) = 0.0_WP 
        ALLOCATE ( V33X1(NCL2,NCL1_io) ) ; V33X1(:,:) = 0.0_WP 
        ALLOCATE ( V12X1(NCL2,NCL1_io) ) ; V12X1(:,:) = 0.0_WP 
        ALLOCATE ( V13X1(NCL2,NCL1_io) ) ; V13X1(:,:) = 0.0_WP 
        ALLOCATE ( V23X1(NCL2,NCL1_io) ) ; V23X1(:,:) = 0.0_WP 
        
        ALLOCATE ( V11X3(NCL2,NCL3) ) ; V11X3(:,:) = 0.0_WP 
        ALLOCATE ( V22X3(NCL2,NCL3) ) ; V22X3(:,:) = 0.0_WP 
        ALLOCATE ( V33X3(NCL2,NCL3) ) ; V33X3(:,:) = 0.0_WP 
        ALLOCATE ( V12X3(NCL2,NCL3) ) ; V12X3(:,:) = 0.0_WP 
        ALLOCATE ( V13X3(NCL2,NCL3) ) ; V13X3(:,:) = 0.0_WP 
        ALLOCATE ( V23X3(NCL2,NCL3) ) ; V23X3(:,:) = 0.0_WP 
        
        !===============
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( ENV11T(NCL2,NCL1_io) ) ; ENV11T(:,:)=0.0_WP 
        ALLOCATE ( ENV22T(NCL2,NCL1_io) ) ; ENV22T(:,:)=0.0_WP 
        ALLOCATE ( ENV33T(NCL2,NCL1_io) ) ; ENV33T(:,:)=0.0_WP 
        ALLOCATE ( ENV12T(NCL2,NCL1_io) ) ; ENV12T(:,:)=0.0_WP 
        ALLOCATE ( ENV13T(NCL2,NCL1_io) ) ; ENV13T(:,:)=0.0_WP 
        ALLOCATE ( ENV23T(NCL2,NCL1_io) ) ; ENV23T(:,:)=0.0_WP 
        
        !=============
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( ENV11Z(NCL2,NCL3) ) ; ENV11Z(:,:)=0.0_WP 
        ALLOCATE ( ENV22Z(NCL2,NCL3) ) ; ENV22Z(:,:)=0.0_WP 
        ALLOCATE ( ENV33Z(NCL2,NCL3) ) ; ENV33Z(:,:)=0.0_WP  
        ALLOCATE ( ENV12Z(NCL2,NCL3) ) ; ENV12Z(:,:)=0.0_WP 
        ALLOCATE ( ENV13Z(NCL2,NCL3) ) ; ENV13Z(:,:)=0.0_WP 
        ALLOCATE ( ENV23Z(NCL2,NCL3) ) ; ENV23Z(:,:)=0.0_WP 
        
        !=============velocity and voriticity==============
        ALLOCATE ( VO11X1(NCL2,NCL1_io) ) ; VO11X1=0.0_WP
        ALLOCATE ( VO12X1(NCL2,NCL1_io) ) ; VO12X1=0.0_WP
        ALLOCATE ( VO13X1(NCL2,NCL1_io) ) ; VO13X1=0.0_WP
        
        ALLOCATE ( VO21X1(NCL2,NCL1_io) ) ; VO21X1=0.0_WP
        ALLOCATE ( VO22X1(NCL2,NCL1_io) ) ; VO22X1=0.0_WP
        ALLOCATE ( VO23X1(NCL2,NCL1_io) ) ; VO23X1=0.0_WP
        
        ALLOCATE ( VO31X1(NCL2,NCL1_io) ) ; VO31X1=0.0_WP
        ALLOCATE ( VO32X1(NCL2,NCL1_io) ) ; VO32X1=0.0_WP
        ALLOCATE ( VO33X1(NCL2,NCL1_io) ) ; VO33X1=0.0_WP
        
        ALLOCATE ( VO11X3(NCL2,NCL3) ) ; VO11X3=0.0_WP
        ALLOCATE ( VO12X3(NCL2,NCL3) ) ; VO12X3=0.0_WP
        ALLOCATE ( VO13X3(NCL2,NCL3) ) ; VO13X3=0.0_WP
        
        ALLOCATE ( VO21X3(NCL2,NCL3) ) ; VO21X3=0.0_WP
        ALLOCATE ( VO22X3(NCL2,NCL3) ) ; VO22X3=0.0_WP
        ALLOCATE ( VO23X3(NCL2,NCL3) ) ; VO23X3=0.0_WP
        
        ALLOCATE ( VO31X3(NCL2,NCL3) ) ; VO31X3=0.0_WP
        ALLOCATE ( VO32X3(NCL2,NCL3) ) ; VO32X3=0.0_WP
        ALLOCATE ( VO33X3(NCL2,NCL3) ) ; VO33X3=0.0_WP
        
        ALLOCATE ( EVO11T(NCL2,NCL1_io) ) ; EVO11T=0.0_WP
        ALLOCATE ( EVO12T(NCL2,NCL1_io) ) ; EVO12T=0.0_WP
        ALLOCATE ( EVO13T(NCL2,NCL1_io) ) ; EVO13T=0.0_WP
        
        ALLOCATE ( EVO21T(NCL2,NCL1_io) ) ; EVO21T=0.0_WP
        ALLOCATE ( EVO22T(NCL2,NCL1_io) ) ; EVO22T=0.0_WP
        ALLOCATE ( EVO23T(NCL2,NCL1_io) ) ; EVO23T=0.0_WP
        
        ALLOCATE ( EVO31T(NCL2,NCL1_io) ) ; EVO31T=0.0_WP
        ALLOCATE ( EVO32T(NCL2,NCL1_io) ) ; EVO32T=0.0_WP
        ALLOCATE ( EVO33T(NCL2,NCL1_io) ) ; EVO33T=0.0_WP
        
        ALLOCATE ( EVO11Z(NCL2,NCL3) ) ; EVO11Z=0.0_WP
        ALLOCATE ( EVO12Z(NCL2,NCL3) ) ; EVO12Z=0.0_WP
        ALLOCATE ( EVO13Z(NCL2,NCL3) ) ; EVO13Z=0.0_WP
        
        ALLOCATE ( EVO21Z(NCL2,NCL3) ) ; EVO21Z=0.0_WP
        ALLOCATE ( EVO22Z(NCL2,NCL3) ) ; EVO22Z=0.0_WP
        ALLOCATE ( EVO23Z(NCL2,NCL3) ) ; EVO23Z=0.0_WP
        
        ALLOCATE ( EVO31Z(NCL2,NCL3) ) ; EVO31Z=0.0_WP
        ALLOCATE ( EVO32Z(NCL2,NCL3) ) ; EVO32Z=0.0_WP
        ALLOCATE ( EVO33Z(NCL2,NCL3) ) ; EVO33Z=0.0_WP
        

        !=============
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( EE1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( EE1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( EE1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( EE1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( EE1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( EE1_23(NCL2, NCL1_io) ) 
        
        ALLOCATE ( UU1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( UU1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( UU1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( UU1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( UU1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( UU1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( HH1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( HH1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( HH1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( HH1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( HH1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( HH1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( VV1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( VV1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( VV1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( VV1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( VV1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( VV1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( EH1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( EH1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( EH1_13(NCL2, NCL1_io) ) 
        
        ALLOCATE ( EH1_21(NCL2, NCL1_io) ) 
        ALLOCATE ( EH1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( EH1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( EH1_31(NCL2, NCL1_io) ) 
        ALLOCATE ( EH1_32(NCL2, NCL1_io) ) 
        ALLOCATE ( EH1_33(NCL2, NCL1_io) ) 

        ALLOCATE ( UV1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( UV1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( UV1_13(NCL2, NCL1_io) ) 

        ALLOCATE ( UV1_21(NCL2, NCL1_io) ) 
        ALLOCATE ( UV1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( UV1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( UV1_31(NCL2, NCL1_io) ) 
        ALLOCATE ( UV1_32(NCL2, NCL1_io) ) 
        ALLOCATE ( UV1_33(NCL2, NCL1_io) ) 
        
        ALLOCATE ( WORK_EE1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EE1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EE1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EE1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EE1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EE1_23(NCL2, NCL1_io) ) 
        
        ALLOCATE ( WORK_UU1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UU1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UU1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UU1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UU1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UU1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( WORK_HH1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_HH1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_HH1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_HH1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_HH1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_HH1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( WORK_VV1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_VV1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_VV1_33(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_VV1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_VV1_13(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_VV1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( WORK_EH1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EH1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EH1_13(NCL2, NCL1_io) ) 
        
        ALLOCATE ( WORK_EH1_21(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EH1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EH1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( WORK_EH1_31(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EH1_32(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_EH1_33(NCL2, NCL1_io) ) 

        ALLOCATE ( WORK_UV1_11(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UV1_12(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UV1_13(NCL2, NCL1_io) ) 

        ALLOCATE ( WORK_UV1_21(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UV1_22(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UV1_23(NCL2, NCL1_io) ) 

        ALLOCATE ( WORK_UV1_31(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UV1_32(NCL2, NCL1_io) ) 
        ALLOCATE ( WORK_UV1_33(NCL2, NCL1_io) ) 
        
        !=============
        ! IN WAVE NUMBER SPACE
        ALLOCATE ( EE3_11(NCL2, NCL3) ) 
        ALLOCATE ( EE3_22(NCL2, NCL3) ) 
        ALLOCATE ( EE3_33(NCL2, NCL3) ) 
        ALLOCATE ( EE3_12(NCL2, NCL3) ) 
        ALLOCATE ( EE3_13(NCL2, NCL3) ) 
        ALLOCATE ( EE3_23(NCL2, NCL3) ) 
        
        ALLOCATE ( UU3_11(NCL2, NCL3) ) 
        ALLOCATE ( UU3_22(NCL2, NCL3) ) 
        ALLOCATE ( UU3_33(NCL2, NCL3) ) 
        ALLOCATE ( UU3_12(NCL2, NCL3) ) 
        ALLOCATE ( UU3_13(NCL2, NCL3) ) 
        ALLOCATE ( UU3_23(NCL2, NCL3) ) 

        ALLOCATE ( HH3_11(NCL2, NCL3) ) 
        ALLOCATE ( HH3_22(NCL2, NCL3) ) 
        ALLOCATE ( HH3_33(NCL2, NCL3) ) 
        ALLOCATE ( HH3_12(NCL2, NCL3) ) 
        ALLOCATE ( HH3_13(NCL2, NCL3) ) 
        ALLOCATE ( HH3_23(NCL2, NCL3) ) 

        ALLOCATE ( VV3_11(NCL2, NCL3) ) 
        ALLOCATE ( VV3_22(NCL2, NCL3) ) 
        ALLOCATE ( VV3_33(NCL2, NCL3) ) 
        ALLOCATE ( VV3_12(NCL2, NCL3) ) 
        ALLOCATE ( VV3_13(NCL2, NCL3) ) 
        ALLOCATE ( VV3_23(NCL2, NCL3) ) 

        ALLOCATE ( EH3_11(NCL2, NCL3) ) 
        ALLOCATE ( EH3_12(NCL2, NCL3) ) 
        ALLOCATE ( EH3_13(NCL2, NCL3) ) 
        
        ALLOCATE ( EH3_21(NCL2, NCL3) ) 
        ALLOCATE ( EH3_22(NCL2, NCL3) ) 
        ALLOCATE ( EH3_23(NCL2, NCL3) ) 

        ALLOCATE ( EH3_31(NCL2, NCL3) ) 
        ALLOCATE ( EH3_32(NCL2, NCL3) ) 
        ALLOCATE ( EH3_33(NCL2, NCL3) ) 

        ALLOCATE ( UV3_11(NCL2, NCL3) ) 
        ALLOCATE ( UV3_12(NCL2, NCL3) ) 
        ALLOCATE ( UV3_13(NCL2, NCL3) ) 

        ALLOCATE ( UV3_21(NCL2, NCL3) ) 
        ALLOCATE ( UV3_22(NCL2, NCL3) ) 
        ALLOCATE ( UV3_23(NCL2, NCL3) ) 

        ALLOCATE ( UV3_31(NCL2, NCL3) ) 
        ALLOCATE ( UV3_32(NCL2, NCL3) ) 
        ALLOCATE ( UV3_33(NCL2, NCL3) ) 
        
        ALLOCATE ( WORK_EE3_11(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EE3_22(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EE3_33(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EE3_12(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EE3_13(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EE3_23(NCL2, NCL3) ) 
        
        ALLOCATE ( WORK_UU3_11(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UU3_22(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UU3_33(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UU3_12(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UU3_13(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UU3_23(NCL2, NCL3) ) 

        ALLOCATE ( WORK_HH3_11(NCL2, NCL3) ) 
        ALLOCATE ( WORK_HH3_22(NCL2, NCL3) ) 
        ALLOCATE ( WORK_HH3_33(NCL2, NCL3) ) 
        ALLOCATE ( WORK_HH3_12(NCL2, NCL3) ) 
        ALLOCATE ( WORK_HH3_13(NCL2, NCL3) ) 
        ALLOCATE ( WORK_HH3_23(NCL2, NCL3) ) 

        ALLOCATE ( WORK_VV3_11(NCL2, NCL3) ) 
        ALLOCATE ( WORK_VV3_22(NCL2, NCL3) ) 
        ALLOCATE ( WORK_VV3_33(NCL2, NCL3) ) 
        ALLOCATE ( WORK_VV3_12(NCL2, NCL3) ) 
        ALLOCATE ( WORK_VV3_13(NCL2, NCL3) ) 
        ALLOCATE ( WORK_VV3_23(NCL2, NCL3) ) 

        ALLOCATE ( WORK_EH3_11(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EH3_12(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EH3_13(NCL2, NCL3) ) 
        
        ALLOCATE ( WORK_EH3_21(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EH3_22(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EH3_23(NCL2, NCL3) ) 

        ALLOCATE ( WORK_EH3_31(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EH3_32(NCL2, NCL3) ) 
        ALLOCATE ( WORK_EH3_33(NCL2, NCL3) ) 

        ALLOCATE ( WORK_UV3_11(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UV3_12(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UV3_13(NCL2, NCL3) ) 

        ALLOCATE ( WORK_UV3_21(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UV3_22(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UV3_23(NCL2, NCL3) ) 

        ALLOCATE ( WORK_UV3_31(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UV3_32(NCL2, NCL3) ) 
        ALLOCATE ( WORK_UV3_33(NCL2, NCL3) ) 
        
    
        RETURN
    END SUBROUTINE
    
    
    SUBROUTINE SPECO_DEALLOCATE
        use SPECO_info
        IMPLICIT NONE
    
         !===========u'=u-<u>=======================================
        DEALLOCATE ( RHS )
        
        !===========energy in each direction ( not used )======================
        DEALLOCATE ( ENEJTF )
        DEALLOCATE ( ENEJME )
        
        DEALLOCATE ( ENE1MA )
        DEALLOCATE ( ENE3MA )
        
        !==========
        ! IN WAVE NUMBER SPACE
        !EN1IK = STREAMWISE INVERSE FFT X1 of u'_ndv
        !EN3KI = SPANWISE   INVERSE FFT X3 of u'_ndv
        DEALLOCATE ( EN1IK )
        DEALLOCATE ( EN3KI )
        
        !=================================================
        ! IN WAVE NUMBER SPACE
        DEALLOCATE ( ENE1 )
        DEALLOCATE ( ENE3 )
        
        !=================================================
        DEALLOCATE ( CORX1 )
        DEALLOCATE ( CORX3 )
        
        !========CORRELATION===========================
        ! CORRELATIONo R(IJ)X(K) = Velocity_i and Velocity_J along the direction K
        DEALLOCATE ( R11X1 )
        DEALLOCATE ( R22X1 )
        DEALLOCATE ( R33X1 )
        DEALLOCATE ( R12X1 )
        DEALLOCATE ( R13X1 )
        DEALLOCATE ( R23X1 )
        
        DEALLOCATE ( R11X3 )
        DEALLOCATE ( R22X3 ) 
        DEALLOCATE ( R33X3 )
        DEALLOCATE ( R12X3 )
        DEALLOCATE ( R13X3 ) 
        DEALLOCATE ( R23X3 )
        
        !===============
        DEALLOCATE ( ENE11T )
        DEALLOCATE ( ENE22T )
        DEALLOCATE ( ENE33T )
        DEALLOCATE ( ENE12T )
        DEALLOCATE ( ENE13T )
        DEALLOCATE ( ENE23T )
 
        !=============
        DEALLOCATE ( ENE11Z )
        DEALLOCATE ( ENE22Z ) 
        DEALLOCATE ( ENE33Z )
        DEALLOCATE ( ENE12Z )
        DEALLOCATE ( ENE13Z )
        DEALLOCATE ( ENE23Z )
        
        DEALLOCATE ( V11X1 )
        DEALLOCATE ( V22X1 )
        DEALLOCATE ( V33X1 )
        DEALLOCATE ( V12X1 )
        DEALLOCATE ( V13X1 )
        DEALLOCATE ( V23X1 )
        
        DEALLOCATE ( V11X3 )
        DEALLOCATE ( V22X3 )
        DEALLOCATE ( V33X3 )
        DEALLOCATE ( V12X3 )
        DEALLOCATE ( V13X3 )
        DEALLOCATE ( V23X3 )
        
        DEALLOCATE ( ENV11T )
        DEALLOCATE ( ENV22T )
        DEALLOCATE ( ENV33T )
        DEALLOCATE ( ENV12T )
        DEALLOCATE ( ENV13T )
        DEALLOCATE ( ENV23T )
        
        DEALLOCATE ( ENV11Z )
        DEALLOCATE ( ENV22Z )
        DEALLOCATE ( ENV33Z )
        DEALLOCATE ( ENV12Z )
        DEALLOCATE ( ENV13Z )
        DEALLOCATE ( ENV23Z )
        
        
        DEALLOCATE ( VO11X1 )
        DEALLOCATE ( VO12X1 )
        DEALLOCATE ( VO13X1 )
        
        DEALLOCATE ( VO21X1 )
        DEALLOCATE ( VO22X1 )
        DEALLOCATE ( VO23X1 )
        
        DEALLOCATE ( VO31X1 )
        DEALLOCATE ( VO32X1 )
        DEALLOCATE ( VO33X1 )
        
        DEALLOCATE ( VO11X3 )
        DEALLOCATE ( VO12X3 )
        DEALLOCATE ( VO13X3 )
        
        DEALLOCATE ( VO21X3 )
        DEALLOCATE ( VO22X3 )
        DEALLOCATE ( VO23X3 )
        
        DEALLOCATE ( VO31X3 )
        DEALLOCATE ( VO32X3 )
        DEALLOCATE ( VO33X3 )
        
        DEALLOCATE ( EVO11T )
        DEALLOCATE ( EVO12T )
        DEALLOCATE ( EVO13T )
        
        DEALLOCATE ( EVO21T )
        DEALLOCATE ( EVO22T )
        DEALLOCATE ( EVO23T )
        
        DEALLOCATE ( EVO31T )
        DEALLOCATE ( EVO32T )
        DEALLOCATE ( EVO33T )
        
        DEALLOCATE ( EVO11Z )
        DEALLOCATE ( EVO12Z )
        DEALLOCATE ( EVO13Z )
        
        DEALLOCATE ( EVO21Z )
        DEALLOCATE ( EVO22Z )
        DEALLOCATE ( EVO23Z )
        
        DEALLOCATE ( EVO31Z )
        DEALLOCATE ( EVO32Z )
        DEALLOCATE ( EVO33Z )
        
        ! IN WAVE NUMBER SPACE
        DEALLOCATE ( EE1_11 ) 
        DEALLOCATE ( EE1_22 ) 
        DEALLOCATE ( EE1_33 ) 
        DEALLOCATE ( EE1_12 ) 
        DEALLOCATE ( EE1_13 ) 
        DEALLOCATE ( EE1_23 ) 
        
        DEALLOCATE ( UU1_11 ) 
        DEALLOCATE ( UU1_22 ) 
        DEALLOCATE ( UU1_33 ) 
        DEALLOCATE ( UU1_12 ) 
        DEALLOCATE ( UU1_13 ) 
        DEALLOCATE ( UU1_23 ) 

        DEALLOCATE ( HH1_11 ) 
        DEALLOCATE ( HH1_22 ) 
        DEALLOCATE ( HH1_33 ) 
        DEALLOCATE ( HH1_12 ) 
        DEALLOCATE ( HH1_13 ) 
        DEALLOCATE ( HH1_23 ) 

        DEALLOCATE ( VV1_11 ) 
        DEALLOCATE ( VV1_22 ) 
        DEALLOCATE ( VV1_33 ) 
        DEALLOCATE ( VV1_12 ) 
        DEALLOCATE ( VV1_13 ) 
        DEALLOCATE ( VV1_23 ) 

        DEALLOCATE ( EH1_11 ) 
        DEALLOCATE ( EH1_12 ) 
        DEALLOCATE ( EH1_13 ) 
        
        DEALLOCATE ( EH1_21 ) 
        DEALLOCATE ( EH1_22 ) 
        DEALLOCATE ( EH1_23 ) 

        DEALLOCATE ( EH1_31 ) 
        DEALLOCATE ( EH1_32 ) 
        DEALLOCATE ( EH1_33 ) 

        DEALLOCATE ( UV1_11 ) 
        DEALLOCATE ( UV1_12 ) 
        DEALLOCATE ( UV1_13 ) 

        DEALLOCATE ( UV1_21 ) 
        DEALLOCATE ( UV1_22 ) 
        DEALLOCATE ( UV1_23 ) 

        DEALLOCATE ( UV1_31 ) 
        DEALLOCATE ( UV1_32 ) 
        DEALLOCATE ( UV1_33 ) 
        
        DEALLOCATE ( WORK_EE1_11 ) 
        DEALLOCATE ( WORK_EE1_22 ) 
        DEALLOCATE ( WORK_EE1_33 ) 
        DEALLOCATE ( WORK_EE1_12 ) 
        DEALLOCATE ( WORK_EE1_13 ) 
        DEALLOCATE ( WORK_EE1_23 ) 
        
        DEALLOCATE ( WORK_UU1_11 ) 
        DEALLOCATE ( WORK_UU1_22 ) 
        DEALLOCATE ( WORK_UU1_33 ) 
        DEALLOCATE ( WORK_UU1_12 ) 
        DEALLOCATE ( WORK_UU1_13 ) 
        DEALLOCATE ( WORK_UU1_23 ) 

        DEALLOCATE ( WORK_HH1_11 ) 
        DEALLOCATE ( WORK_HH1_22 ) 
        DEALLOCATE ( WORK_HH1_33 ) 
        DEALLOCATE ( WORK_HH1_12 ) 
        DEALLOCATE ( WORK_HH1_13 ) 
        DEALLOCATE ( WORK_HH1_23 ) 

        DEALLOCATE ( WORK_VV1_11 ) 
        DEALLOCATE ( WORK_VV1_22 ) 
        DEALLOCATE ( WORK_VV1_33 ) 
        DEALLOCATE ( WORK_VV1_12 ) 
        DEALLOCATE ( WORK_VV1_13 ) 
        DEALLOCATE ( WORK_VV1_23 ) 

        DEALLOCATE ( WORK_EH1_11 ) 
        DEALLOCATE ( WORK_EH1_12 ) 
        DEALLOCATE ( WORK_EH1_13 ) 
        
        DEALLOCATE ( WORK_EH1_21 ) 
        DEALLOCATE ( WORK_EH1_22 ) 
        DEALLOCATE ( WORK_EH1_23 ) 

        DEALLOCATE ( WORK_EH1_31 ) 
        DEALLOCATE ( WORK_EH1_32 ) 
        DEALLOCATE ( WORK_EH1_33 ) 

        DEALLOCATE ( WORK_UV1_11 ) 
        DEALLOCATE ( WORK_UV1_12 ) 
        DEALLOCATE ( WORK_UV1_13 ) 

        DEALLOCATE ( WORK_UV1_21 ) 
        DEALLOCATE ( WORK_UV1_22 ) 
        DEALLOCATE ( WORK_UV1_23 ) 

        DEALLOCATE ( WORK_UV1_31 ) 
        DEALLOCATE ( WORK_UV1_32 ) 
        DEALLOCATE ( WORK_UV1_33 ) 
        
        !=============
        ! IN WAVE NUMBER SPACE
        DEALLOCATE ( EE3_11 ) 
        DEALLOCATE ( EE3_22 ) 
        DEALLOCATE ( EE3_33 ) 
        DEALLOCATE ( EE3_12 ) 
        DEALLOCATE ( EE3_13 ) 
        DEALLOCATE ( EE3_23 ) 
        
        DEALLOCATE ( UU3_11 ) 
        DEALLOCATE ( UU3_22 ) 
        DEALLOCATE ( UU3_33 ) 
        DEALLOCATE ( UU3_12 ) 
        DEALLOCATE ( UU3_13 ) 
        DEALLOCATE ( UU3_23 ) 

        DEALLOCATE ( HH3_11 ) 
        DEALLOCATE ( HH3_22 ) 
        DEALLOCATE ( HH3_33 ) 
        DEALLOCATE ( HH3_12 ) 
        DEALLOCATE ( HH3_13 ) 
        DEALLOCATE ( HH3_23 ) 

        DEALLOCATE ( VV3_11 ) 
        DEALLOCATE ( VV3_22 ) 
        DEALLOCATE ( VV3_33 ) 
        DEALLOCATE ( VV3_12 ) 
        DEALLOCATE ( VV3_13 ) 
        DEALLOCATE ( VV3_23 ) 

        DEALLOCATE ( EH3_11 ) 
        DEALLOCATE ( EH3_12 ) 
        DEALLOCATE ( EH3_13 ) 
        
        DEALLOCATE ( EH3_21 ) 
        DEALLOCATE ( EH3_22 ) 
        DEALLOCATE ( EH3_23 ) 

        DEALLOCATE ( EH3_31 ) 
        DEALLOCATE ( EH3_32 ) 
        DEALLOCATE ( EH3_33 ) 

        DEALLOCATE ( UV3_11 ) 
        DEALLOCATE ( UV3_12 ) 
        DEALLOCATE ( UV3_13 ) 

        DEALLOCATE ( UV3_21 ) 
        DEALLOCATE ( UV3_22 ) 
        DEALLOCATE ( UV3_23 ) 

        DEALLOCATE ( UV3_31 ) 
        DEALLOCATE ( UV3_32 ) 
        DEALLOCATE ( UV3_33 ) 
        
        DEALLOCATE ( WORK_EE3_11 ) 
        DEALLOCATE ( WORK_EE3_22 ) 
        DEALLOCATE ( WORK_EE3_33 ) 
        DEALLOCATE ( WORK_EE3_12 ) 
        DEALLOCATE ( WORK_EE3_13 ) 
        DEALLOCATE ( WORK_EE3_23 ) 
        
        DEALLOCATE ( WORK_UU3_11 ) 
        DEALLOCATE ( WORK_UU3_22 ) 
        DEALLOCATE ( WORK_UU3_33 ) 
        DEALLOCATE ( WORK_UU3_12 ) 
        DEALLOCATE ( WORK_UU3_13 ) 
        DEALLOCATE ( WORK_UU3_23 ) 

        DEALLOCATE ( WORK_HH3_11 ) 
        DEALLOCATE ( WORK_HH3_22 ) 
        DEALLOCATE ( WORK_HH3_33 ) 
        DEALLOCATE ( WORK_HH3_12 ) 
        DEALLOCATE ( WORK_HH3_13 ) 
        DEALLOCATE ( WORK_HH3_23 ) 

        DEALLOCATE ( WORK_VV3_11 ) 
        DEALLOCATE ( WORK_VV3_22 ) 
        DEALLOCATE ( WORK_VV3_33 ) 
        DEALLOCATE ( WORK_VV3_12 ) 
        DEALLOCATE ( WORK_VV3_13 ) 
        DEALLOCATE ( WORK_VV3_23 ) 

        DEALLOCATE ( WORK_EH3_11 ) 
        DEALLOCATE ( WORK_EH3_12 ) 
        DEALLOCATE ( WORK_EH3_13 ) 
        
        DEALLOCATE ( WORK_EH3_21 ) 
        DEALLOCATE ( WORK_EH3_22 ) 
        DEALLOCATE ( WORK_EH3_23 ) 

        DEALLOCATE ( WORK_EH3_31 ) 
        DEALLOCATE ( WORK_EH3_32 ) 
        DEALLOCATE ( WORK_EH3_33 ) 

        DEALLOCATE ( WORK_UV3_11 ) 
        DEALLOCATE ( WORK_UV3_12 ) 
        DEALLOCATE ( WORK_UV3_13 ) 

        DEALLOCATE ( WORK_UV3_21 ) 
        DEALLOCATE ( WORK_UV3_22 ) 
        DEALLOCATE ( WORK_UV3_23 ) 

        DEALLOCATE ( WORK_UV3_31 ) 
        DEALLOCATE ( WORK_UV3_32 ) 
        DEALLOCATE ( WORK_UV3_33 ) 
    
        RETURN
    END SUBROUTINE
    
!====================main code============================================

    SUBROUTINE PP_SPECOSPEC
    !   Refer to: Page 213, Chapter 9.9e in Orland book.
        use mesh_info
        use init_info
        use flow_info
        USE thermal_info
        use SPECO_info
        use postprocess_info
        USE WRT_INFO
        
        IMPLICIT NONE
        
        REAL(WP) :: VARINST(NCL1S:NCL1E,1:N2DO(0),NCL3,6)
        REAL(WP) :: VARAVRG(NCL2,6)
        REAL(WP) :: DVDL_cct(NDV,NDV)
        
        INTEGER(4) :: I, J, K, IP, JP, KP
        
        IF(ppspectra==0) RETURN
        
        !==========COMMON=======================
        N3MH=NCL3/2+1
        N3MD=NCL3+2
        N1MH=NCL1_io/2+1
        N1MD=NCL1_io+2
        CALL PP_wall_thermal_shear(flgxz)
        
        !========OUTPUT FLOW VELOCITY CORRECTION========================
        IF(MYID.EQ.0) CALL CHKRLHDL('16.IO: postprocessing spectra data (FLOW) at',myid,phyTIME_io)
        DO J=1, N2DO(MYID)
            JP = JLPV(J)
            DO I=1, NCL1_IO
                IP = IPV_IO(I)
                DO K=1, NCL3
                    KP = KPV(K)
                    VARINST(I,J,K,1)=(Q_io(I,J,K,1)+Q_io(IP,J,K,1))*0.5_WP ! Cell centre
                    VARINST(I,J,K,2)=(Q_io(I,J,K,2)+Q_io(I,JP,K,2))*YND2CL
                    VARINST(I,J,K,3)=(Q_io(I,J,K,3)+Q_io(I,J,KP,3))*0.5_WP
                    CALL INST_Qio_GRADIENT_CELLCENTRED(I, J, K, DVDL_cct)
                    VARINST(I,J,K,4)= DVDL_cct(3,2)-DVDL_cct(2,3)
                    VARINST(I,J,K,5)= DVDL_cct(1,3)-DVDL_cct(3,1)
                    VARINST(I,J,K,6)= DVDL_cct(2,1)-DVDL_cct(1,2)
                END DO
            END DO
            VARAVRG(J,1)    =U1xzL_io(J,1)
            VARAVRG(J,2)    =U1xzL_io(J,2)
            VARAVRG(J,3)    =U1xzL_io(J,3)
            VARAVRG(J,4)    =DVDL1xzL_io(J,3,2)-DVDL1xzL_io(J,2,3)
            VARAVRG(J,5)    =DVDL1xzL_io(J,1,3)-DVDL1xzL_io(J,3,1)
            VARAVRG(J,6)    =DVDL1xzL_io(J,2,1)-DVDL1xzL_io(J,1,2)
        END DO
        
        CALL SPECO_ALLOCATE
        CALL SPECO_CORR_ENEG(VARINST,VARAVRG)
        CALL SPECO_FOR_AVERAGE('FLOW')
        CALL SPECO_WRITE_PROFILE('FLOW')
        CALL SPECO_WRITE_Contour('FLOW')
        CALL SPECO_DEALLOCATE
        IF(MYID.EQ.0) CALL CHKRLHDL('       finished postprocessing spectra data (FLOW) at',myid,phyTIME_io)
        
    
        !========OUTPUT THERMAL CORRECTION========================
        IF(thermlflg==1) THEN
            IF(MYID.EQ.0) CALL CHKRLHDL('16.IO: postprocessing spectra data (HEAT) at',myid,phyTIME_io)
            DO J=1, N2DO(MYID)
                JP = JLPV(J)
                DO I=1, NCL1_IO
                    IP = IPV_IO(I)
                    DO K=1, NCL3
                        KP = KPV(K)
                        VARINST(I,J,K,1)=(Q_io(I,J,K,1)+Q_io(IP,J,K,1))*0.5_WP ! Cell centre
                        VARINST(I,J,K,2)=(Q_io(I,J,K,2)+Q_io(I,JP,K,2))*YND2CL
                        VARINST(I,J,K,3)=(Q_io(I,J,K,3)+Q_io(I,J,KP,3))*0.5_WP
                        VARINST(I,J,K,4)= density    (I,J,K)
                        VARINST(I,J,K,5)= ENTHALPY   (I,J,K)
                        VARINST(I,J,K,6)= TEMPERATURE(I,J,K)
                    END DO
                END DO
                VARAVRG(J,1)    = U1xzL_io(J,1)
                VARAVRG(J,2)    = U1xzL_io(J,2)
                VARAVRG(J,3)    = U1xzL_io(J,3)
                VARAVRG(J,4)    = D1XzL_io(J)
                VARAVRG(J,5)    = H1XzL_io(J)
                VARAVRG(J,6)    = T1XzL_io(J)
            END DO
            
            CALL SPECO_ALLOCATE
            CALL SPECO_CORR_ENEG(VARINST,VARAVRG)
            CALL SPECO_FOR_AVERAGE('HEAT')
            CALL SPECO_WRITE_PROFILE('HEAT')
            CALL SPECO_WRITE_Contour('HEAT')
            CALL SPECO_DEALLOCATE
            IF(MYID.EQ.0) CALL CHKRLHDL('       finished postprocessing spectra data (HEAT) at',myid,phyTIME_io)
        END IF
        
        
        RETURN
    END SUBROUTINE

    SUBROUTINE SPECO_FOR_AVERAGE(STR)
        use mesh_info
        use init_info
        use SPECO_info
        use postprocess_info
        USE WRT_INFO
        
        IMPLICIT NONE
        CHARACTER(4),INTENT(IN)  :: STR
        INTEGER(4) :: N
        
        IF(MYID.NE.0) RETURN
        N=100
        
        IF(TRIM(STR)=='FLOW') N=1
        IF(TRIM(STR)=='HEAT') N=2
        
        !================Velocity =====================
        R11X1_xzLa(:, :, N)  = R11X1(:,:)
        R22X1_xzLa(:, :, N)  = R22X1(:,:)
        R33X1_xzLa(:, :, N)  = R33X1(:,:)
        R12X1_xzLa(:, :, N)  = R12X1(:,:)
        R13X1_xzLa(:, :, N)  = R13X1(:,:)
        R23X1_xzLa(:, :, N)  = R23X1(:,:)
        
        R11X3_xzLa(:, :, N)  = R11X3(:,:)
        R22X3_xzLa(:, :, N)  = R22X3(:,:)
        R33X3_xzLa(:, :, N)  = R33X3(:,:)
        R12X3_xzLa(:, :, N)  = R12X3(:,:)
        R13X3_xzLa(:, :, N)  = R13X3(:,:)
        R23X3_xzLa(:, :, N)  = R23X3(:,:)
        
        ENE11T_xzLa(:, :, N)  = ENE11T(:,:)
        ENE22T_xzLa(:, :, N)  = ENE22T(:,:)
        ENE33T_xzLa(:, :, N)  = ENE33T(:,:)
        ENE12T_xzLa(:, :, N)  = ENE12T(:,:)
        ENE13T_xzLa(:, :, N)  = ENE13T(:,:)
        ENE23T_xzLa(:, :, N)  = ENE23T(:,:)
        
        ENE11Z_xzLa(:, :, N)  = ENE11Z(:,:)
        ENE22Z_xzLa(:, :, N)  = ENE22Z(:,:)
        ENE33Z_xzLa(:, :, N)  = ENE33Z(:,:)
        ENE12Z_xzLa(:, :, N)  = ENE12Z(:,:)
        ENE13Z_xzLa(:, :, N)  = ENE13Z(:,:)
        ENE23Z_xzLa(:, :, N)  = ENE23Z(:,:)
        
        !================Voriticity ====================
        V11X1_xzLa(:, :, N)  = V11X1(:,:)
        V22X1_xzLa(:, :, N)  = V22X1(:,:)
        V33X1_xzLa(:, :, N)  = V33X1(:,:)
        V12X1_xzLa(:, :, N)  = V12X1(:,:)
        V13X1_xzLa(:, :, N)  = V13X1(:,:)
        V23X1_xzLa(:, :, N)  = V23X1(:,:)
        
        V11X3_xzLa(:, :, N)  = V11X3(:,:)
        V22X3_xzLa(:, :, N)  = V22X3(:,:)
        V33X3_xzLa(:, :, N)  = V33X3(:,:)
        V12X3_xzLa(:, :, N)  = V12X3(:,:)
        V13X3_xzLa(:, :, N)  = V13X3(:,:)
        V23X3_xzLa(:, :, N)  = V23X3(:,:)
        
        ENV11T_xzLa(:, :, N)  = ENV11T(:,:)
        ENV22T_xzLa(:, :, N)  = ENV22T(:,:)
        ENV33T_xzLa(:, :, N)  = ENV33T(:,:)
        ENV12T_xzLa(:, :, N)  = ENV12T(:,:)
        ENV13T_xzLa(:, :, N)  = ENV13T(:,:)
        ENV23T_xzLa(:, :, N)  = ENV23T(:,:)
        
        ENV11Z_xzLa(:, :, N)  = ENV11Z(:,:)
        ENV22Z_xzLa(:, :, N)  = ENV22Z(:,:)
        ENV33Z_xzLa(:, :, N)  = ENV33Z(:,:)
        ENV12Z_xzLa(:, :, N)  = ENV12Z(:,:)
        ENV13Z_xzLa(:, :, N)  = ENV13Z(:,:)
        ENV23Z_xzLa(:, :, N)  = ENV23Z(:,:)
        
        
        !===============Voriticity & Velocity=========================
        VO11X1_xzLa(:, :, N)  = VO11X1(:,:)
        VO12X1_xzLa(:, :, N)  = VO12X1(:,:)
        VO13X1_xzLa(:, :, N)  = VO13X1(:,:)
        
        VO21X1_xzLa(:, :, N)  = VO21X1(:,:)
        VO22X1_xzLa(:, :, N)  = VO22X1(:,:)
        VO23X1_xzLa(:, :, N)  = VO23X1(:,:)
        
        VO31X1_xzLa(:, :, N)  = VO31X1(:,:)
        VO32X1_xzLa(:, :, N)  = VO32X1(:,:)
        VO33X1_xzLa(:, :, N)  = VO33X1(:,:)
        
        VO11X3_xzLa(:, :, N)  = VO11X3(:,:)
        VO12X3_xzLa(:, :, N)  = VO12X3(:,:)
        VO13X3_xzLa(:, :, N)  = VO13X3(:,:)
        
        VO21X3_xzLa(:, :, N)  = VO21X3(:,:)
        VO22X3_xzLa(:, :, N)  = VO22X3(:,:)
        VO23X3_xzLa(:, :, N)  = VO23X3(:,:)
        
        VO31X3_xzLa(:, :, N)  = VO31X3(:,:)
        VO32X3_xzLa(:, :, N)  = VO32X3(:,:)
        VO33X3_xzLa(:, :, N)  = VO33X3(:,:)
        
        EVO11T_xzLa(:, :, N)  = EVO11T(:,:)
        EVO12T_xzLa(:, :, N)  = EVO12T(:,:)
        EVO13T_xzLa(:, :, N)  = EVO13T(:,:)
        
        EVO21T_xzLa(:, :, N)  = EVO21T(:,:)
        EVO22T_xzLa(:, :, N)  = EVO22T(:,:)
        EVO23T_xzLa(:, :, N)  = EVO23T(:,:)
        
        EVO31T_xzLa(:, :, N)  = EVO31T(:,:)
        EVO32T_xzLa(:, :, N)  = EVO32T(:,:)
        EVO33T_xzLa(:, :, N)  = EVO33T(:,:)
        
        EVO11Z_xzLa(:, :, N)  = EVO11Z(:,:)
        EVO12Z_xzLa(:, :, N)  = EVO12Z(:,:)
        EVO13Z_xzLa(:, :, N)  = EVO13Z(:,:)
        
        EVO21Z_xzLa(:, :, N)  = EVO21Z(:,:)
        EVO22Z_xzLa(:, :, N)  = EVO22Z(:,:)
        EVO23Z_xzLa(:, :, N)  = EVO23Z(:,:)
        
        EVO31Z_xzLa(:, :, N)  = EVO31Z(:,:)
        EVO32Z_xzLa(:, :, N)  = EVO32Z(:,:)
        EVO33Z_xzLa(:, :, N)  = EVO33Z(:,:)
        
        
        
        RETURN
    END SUBROUTINE
!=============================================================================================
    SUBROUTINE SPECO_WRITE_PROFILE(STR)
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
        
        INTEGER(4)    :: N, JJ, L, KC, IC
        
        IF(MYID.NE.0) RETURN
        !==========test for a ASCII output==============
        !Ret_ave  = 0.5_wp*(Ret_io(1)+Ret_io(2)) !;write(*,*) Ret_ave, Ret_io(1), Ret_io(2) !test
        !U_tau_ave= 0.5_wp*(Utaw_io(1)+Utaw_io(2))
        OPEN(100, FILE=TRIM(filepath5)//'CHK_PROBE_for_spectra_instant.dat')
        WRITE(100,'(A)') '## MGRID, JJ, YCC, Yplus1, Yplus2, Utaw_io(1:2), Ret_io(1:2), Ret_ave_io, Ldist_io(1:2)'
        DO N=1, MGRID
            JJ=JGMOV(N)
            DFLG(N) = 100+N
            
            WRITE(PNTIM,'(1ES15.9)') phyTIME_io
            WRITE(PNLOC,'(1I3.3)')   N
            
            !==============correlations======================================
            DO L=1,2
                IF(L==1) THEN
                    FLNAME= TRIM(filepath5)//'Result.IO.Spectral.'//TRIM(STR)//'.2PCorrelation.X.T'  &
                            //TRIM(PNTIM)//'.YLC'//TRIM(PNLOC)//'.instant.tec'
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
                    WRITE(100,'(2I4.2, 10ES13.5)') N, JJ, YCC(JJ), WDS(1:2)*Ret_io(1:2)/Ldist_io(1:2), &
                                                   Utaw_io(1:2), Ret_io(1:2), Ret_ave_io, Ldist_io(1:2) 
                    DO IC=1,N1MH
                        WRITE(DFLG(N),'(22ES13.5)') 0.5_WP*( XND_io(IC) + XND_io(IC+1) ), &
                                    R11X1 (JJ,IC), R22X1 (JJ,IC), R33X1 (JJ,IC), &
                                    R12X1 (JJ,IC), R13X1 (JJ,IC), R23X1 (JJ,IC), &
                                    V11X1 (JJ,IC), V22X1 (JJ,IC), V33X1 (JJ,IC), &
                                    V12X1 (JJ,IC), V13X1 (JJ,IC), V23X1 (JJ,IC), &
                                    VO11X1(JJ,IC), VO12X1(JJ,IC), VO13X1(JJ,IC), &
                                    VO21X1(JJ,IC), VO22X1(JJ,IC), VO23X1(JJ,IC), &
                                    VO31X1(JJ,IC), VO32X1(JJ,IC), VO33X1(JJ,IC)
                    ENDDO
                    CLOSE(DFLG(N))
                END IF
                
                IF(L==2) THEN
                    FLNAME= TRIM(filepath5)//'Result.IO.Spectral.'//TRIM(STR)//'.2PCorrelation.Z.T'  &
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
                                R11X3 (JJ,KC), R22X3 (JJ,KC), R33X3 (JJ,KC), &
                                R12X3 (JJ,KC), R13X3 (JJ,KC), R23X3 (JJ,KC), &
                                V11X3 (JJ,KC), V22X3 (JJ,KC), V33X3 (JJ,KC), &
                                V12X3 (JJ,KC), V13X3 (JJ,KC), V23X3 (JJ,KC), &
                                VO11X3(JJ,KC), VO12X3(JJ,KC), VO13X3(JJ,KC), &
                                VO21X3(JJ,KC), VO22X3(JJ,KC), VO23X3(JJ,KC), &
                                VO31X3(JJ,KC), VO32X3(JJ,KC), VO33X3(JJ,KC)
                    ENDDO
                    CLOSE(DFLG(N))
                
                END IF
            END DO
            
            !==============energy espectra======================================
            DO L=1,2
                IF(L==1) THEN
                    FLNAME= TRIM(filepath5)//'Result.IO.Spectral.'//TRIM(STR)//'.energy.Xwavenumber.T'  &
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
                        'ZONE T= "Ekk(K1) at y/delta= ',YCC(JJ), ' Utauw12= ',Utaw_io(1:2), ' Ret12a= ',Ret_io(1:2),Ret_ave_io,' " '
                    DO IC=1,N1MH
                        AKE=( DBLE(IC-1)*2.0_wp*PI/HX_io )
                        WRITE(DFLG(N),'(1I8.1, 22ES13.5)') IC, AKE, &
                                ENE11T(JJ,IC), ENE22T(JJ,IC), ENE33T(JJ,IC), &
                                ENE12T(JJ,IC), ENE13T(JJ,IC), ENE23T(JJ,IC), &
                                ENV11T(JJ,IC), ENV22T(JJ,IC), ENV33T(JJ,IC), &
                                ENV12T(JJ,IC), ENV13T(JJ,IC), ENV23T(JJ,IC), &
                                EVO11T(JJ,IC), EVO12T(JJ,IC), EVO13T(JJ,IC), &
                                EVO21T(JJ,IC), EVO22T(JJ,IC), EVO23T(JJ,IC), &
                                EVO31T(JJ,IC), EVO32T(JJ,IC), EVO33T(JJ,IC)
                    ENDDO
                    CLOSE(DFLG(N))
                END IF
                
                IF(L==2) THEN
                    FLNAME= TRIM(filepath5)//'Result.IO.Spectral.'//TRIM(STR)//'.energy.Zwavenumber.T'  &
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
                        'ZONE T= "Ekk(K3) at y/delta= ',YCC(JJ), ' Utauw12= ',Utaw_io(1:2), ' Ret12a= ',Ret_io(1:2),Ret_ave_io,' " '
                    DO KC=1,N3MH
                        AKE=( DBLE(KC-1)*2.0_wp*PI/HZ )
                        WRITE(DFLG(N),'(1I8.1, 22ES13.5)') KC, AKE, &
                                ENE11Z(JJ,KC), ENE22Z(JJ,KC), ENE33Z(JJ,KC), &
                                ENE12Z(JJ,KC), ENE13Z(JJ,KC), ENE23Z(JJ,KC), &
                                ENV11Z(JJ,KC), ENV22Z(JJ,KC), ENV33Z(JJ,KC), &
                                ENV12Z(JJ,KC), ENV13Z(JJ,KC), ENV23Z(JJ,KC), &
                                EVO11Z(JJ,KC), EVO12Z(JJ,KC), EVO13Z(JJ,KC), &
                                EVO21Z(JJ,KC), EVO22Z(JJ,KC), EVO23Z(JJ,KC), &
                                EVO31Z(JJ,KC), EVO32Z(JJ,KC), EVO33Z(JJ,KC)
                    ENDDO
                    CLOSE(DFLG(N))
                END IF
            END DO
            
        END DO
        CLOSE(100)
            
        
        RETURN
    END SUBROUTINE
    
    
!=============================================================================================
    SUBROUTINE SPECO_WRITE_Contour(STR)
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
         
        INTEGER(4)    :: N, JJ, L, KC, IC, JJM, JJC
        
        IF(MYID.NE.0) RETURN
        
        !==========test for a ASCII output==============
        !Ret_ave  = 0.5_wp*(Ret_io(1)+Ret_io(2)) !;write(*,*) Ret_ave, Ret_io(1), Ret_io(2) !test
        !U_tau_ave= 0.5_wp*(Utaw_io(1)+Utaw_io(2))
        
        !===============plane y-z=====================
        FLNAME = TRIM(filepath5)//'Result.IO.Spectral.'//TRIM(STR)//'.Contours.yz.'//TRIM(PNTIM)//'.tec' 
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
            AKE=( DBLE(KC-1)*2.0_wp*PI/HZ ) 
            DO JJ=1,NND2
                JJM=JGMV(JJ)
                JJC=JJ
                yplus=(1.0_wp-dabs(YND(JJ)))*REN*Utaw_ave_io
                
                IF(JJ==1)    JJM=1
                IF(JJ==NND2) JJC =NCL2
                WRITE(DFLG,'(47ES15.7)') 0.0_WP,YND(JJ),ZND(KC), yplus, AKE, &
                    0.5*( R11X3 (JJC,KC)+R11X3 (JJM,KC) ), &
                    0.5*( R22X3 (JJC,KC)+R22X3 (JJM,KC) ), &
                    0.5*( R33X3 (JJC,KC)+R33X3 (JJM,KC) ), &
                    0.5*( R12X3 (JJC,KC)+R12X3 (JJM,KC) ), &
                    0.5*( R13X3 (JJC,KC)+R13X3 (JJM,KC) ), &
                    0.5*( R23X3 (JJC,KC)+R23X3 (JJM,KC) ), &
                    0.5*( V11X3 (JJC,KC)+V11X3 (JJM,KC) ), &
                    0.5*( V22X3 (JJC,KC)+V22X3 (JJM,KC) ), &
                    0.5*( V33X3 (JJC,KC)+V33X3 (JJM,KC) ), &
                    0.5*( V12X3 (JJC,KC)+V12X3 (JJM,KC) ), &
                    0.5*( V13X3 (JJC,KC)+V13X3 (JJM,KC) ), &
                    0.5*( V23X3 (JJC,KC)+V23X3 (JJM,KC) ), &
                    0.5*( VO11X3 (JJC,KC)+VO11X3 (JJM,KC) ), &
                    0.5*( VO12X3 (JJC,KC)+VO12X3 (JJM,KC) ), &
                    0.5*( VO13X3 (JJC,KC)+VO13X3 (JJM,KC) ), &
                    0.5*( VO21X3 (JJC,KC)+VO21X3 (JJM,KC) ), &
                    0.5*( VO22X3 (JJC,KC)+VO22X3 (JJM,KC) ), &
                    0.5*( VO23X3 (JJC,KC)+VO23X3 (JJM,KC) ), &
                    0.5*( VO31X3 (JJC,KC)+VO31X3 (JJM,KC) ), &
                    0.5*( VO32X3 (JJC,KC)+VO32X3 (JJM,KC) ), &
                    0.5*( VO33X3 (JJC,KC)+VO33X3 (JJM,KC) ), &
                    0.5*( ENE11Z(JJC,KC)+ENE11Z(JJM,KC) ), &
                    0.5*( ENE22Z(JJC,KC)+ENE22Z(JJM,KC) ), &
                    0.5*( ENE33Z(JJC,KC)+ENE33Z(JJM,KC) ), &
                    0.5*( ENE12Z(JJC,KC)+ENE12Z(JJM,KC) ), &
                    0.5*( ENE13Z(JJC,KC)+ENE13Z(JJM,KC) ), &
                    0.5*( ENE23Z(JJC,KC)+ENE23Z(JJM,KC) ), &
                    0.5*( ENV11Z(JJC,KC)+ENV11Z(JJM,KC) ), &
                    0.5*( ENV22Z(JJC,KC)+ENV22Z(JJM,KC) ), &
                    0.5*( ENV33Z(JJC,KC)+ENV33Z(JJM,KC) ), &
                    0.5*( ENV12Z(JJC,KC)+ENV12Z(JJM,KC) ), &
                    0.5*( ENV13Z(JJC,KC)+ENV13Z(JJM,KC) ), &
                    0.5*( ENV23Z(JJC,KC)+ENV23Z(JJM,KC) ), &
                    0.5*( EVO11Z(JJC,KC)+EVO11Z(JJM,KC) ), &
                    0.5*( EVO12Z(JJC,KC)+EVO12Z(JJM,KC) ), &
                    0.5*( EVO13Z(JJC,KC)+EVO13Z(JJM,KC) ), &
                    0.5*( EVO21Z(JJC,KC)+EVO21Z(JJM,KC) ), &
                    0.5*( EVO22Z(JJC,KC)+EVO22Z(JJM,KC) ), &
                    0.5*( EVO23Z(JJC,KC)+EVO23Z(JJM,KC) ), &
                    0.5*( EVO31Z(JJC,KC)+EVO31Z(JJM,KC) ), &
                    0.5*( EVO32Z(JJC,KC)+EVO32Z(JJM,KC) ), &
                    0.5*( EVO33Z(JJC,KC)+EVO33Z(JJM,KC) )

            END DO
        END DO
        CLOSE(DFLG)
        !===============plane x-y=====================
        FLNAME = TRIM(filepath5)//'Result.IO.Spectral.'//TRIM(STR)//'.Contours.yx.'//TRIM(PNTIM)//'.tec' 
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
                AKE=( DBLE(IC-1)*2.0_wp*PI/HX_io ) 
                
                WRITE(DFLG,'(47ES15.7)') XND_io(IC),YND(JJ),0.0_wp, yplus, AKE, &
                    0.5*( R11X1 (JJC,IC)+R11X1 (JJM,IC) ), &
                    0.5*( R22X1 (JJC,IC)+R22X1 (JJM,IC) ), &
                    0.5*( R33X1 (JJC,IC)+R33X1 (JJM,IC) ), &
                    0.5*( R12X1 (JJC,IC)+R12X1 (JJM,IC) ), &
                    0.5*( R13X1 (JJC,IC)+R13X1 (JJM,IC) ), &
                    0.5*( R23X1 (JJC,IC)+R23X1 (JJM,IC) ), &
                    0.5*( V11X1 (JJC,IC)+V11X1 (JJM,IC) ), &
                    0.5*( V22X1 (JJC,IC)+V22X1 (JJM,IC) ), &
                    0.5*( V33X1 (JJC,IC)+V33X1 (JJM,IC) ), &
                    0.5*( V12X1 (JJC,IC)+V12X1 (JJM,IC) ), &
                    0.5*( V13X1 (JJC,IC)+V13X1 (JJM,IC) ), &
                    0.5*( V23X1 (JJC,IC)+V23X1 (JJM,IC) ), &
                    0.5*( VO11X1 (JJC,IC)+VO11X1 (JJM,IC) ), &
                    0.5*( VO12X1 (JJC,IC)+VO12X1 (JJM,IC) ), &
                    0.5*( VO13X1 (JJC,IC)+VO13X1 (JJM,IC) ), &
                    0.5*( VO21X1 (JJC,IC)+VO21X1 (JJM,IC) ), &
                    0.5*( VO22X1 (JJC,IC)+VO22X1 (JJM,IC) ), &
                    0.5*( VO23X1 (JJC,IC)+VO23X1 (JJM,IC) ), &
                    0.5*( VO31X1 (JJC,IC)+VO31X1 (JJM,IC) ), &
                    0.5*( VO32X1 (JJC,IC)+VO32X1 (JJM,IC) ), &
                    0.5*( VO33X1 (JJC,IC)+VO33X1 (JJM,IC) ), &
                    0.5*( ENE11T(JJC,IC)+ENE11T(JJM,IC) ), &
                    0.5*( ENE22T(JJC,IC)+ENE22T(JJM,IC) ), &
                    0.5*( ENE33T(JJC,IC)+ENE33T(JJM,IC) ), &
                    0.5*( ENE12T(JJC,IC)+ENE12T(JJM,IC) ), &
                    0.5*( ENE13T(JJC,IC)+ENE13T(JJM,IC) ), &
                    0.5*( ENE23T(JJC,IC)+ENE23T(JJM,IC) ), &
                    0.5*( ENV11T(JJC,IC)+ENV11T(JJM,IC) ), &
                    0.5*( ENV22T(JJC,IC)+ENV22T(JJM,IC) ), &
                    0.5*( ENV33T(JJC,IC)+ENV33T(JJM,IC) ), &
                    0.5*( ENV12T(JJC,IC)+ENV12T(JJM,IC) ), &
                    0.5*( ENV13T(JJC,IC)+ENV13T(JJM,IC) ), &
                    0.5*( ENV23T(JJC,IC)+ENV23T(JJM,IC) ), &
                    0.5*( EVO11T(JJC,IC)+EVO11T(JJM,IC) ), &
                    0.5*( EVO12T(JJC,IC)+EVO12T(JJM,IC) ), &
                    0.5*( EVO13T(JJC,IC)+EVO13T(JJM,IC) ), &
                    0.5*( EVO21T(JJC,IC)+EVO21T(JJM,IC) ), &
                    0.5*( EVO22T(JJC,IC)+EVO22T(JJM,IC) ), &
                    0.5*( EVO23T(JJC,IC)+EVO23T(JJM,IC) ), &
                    0.5*( EVO31T(JJC,IC)+EVO31T(JJM,IC) ), &
                    0.5*( EVO32T(JJC,IC)+EVO32T(JJM,IC) ), &
                    0.5*( EVO33T(JJC,IC)+EVO33T(JJM,IC) )
            ENDDO
        END DO
        
        CLOSE(DFLG)
            
        
        RETURN
    END SUBROUTINE
        
        
!==========================================================================
    SUBROUTINE SPECO_CORR_ENEG(VARINST,VARAVRG)
    !   Refer to: Page 213, Chapter 9.9e in Orland book.
        use mesh_info
        use SPECO_info
        
        IMPLICIT NONE
        
        REAL(WP),INTENT(IN) :: VARINST(NCL1S:NCL1E,1:N2DO(0),NCL3,6)
        REAL(WP),INTENT(IN) :: VARAVRG(NCL3,6)
        
        INTEGER(4) :: K, KK, KC, KD, KP
        INTEGER(4) :: I, II, IC, ID, IP
        INTEGER(4) :: JJ, J, L
        INTEGER(4) :: INUMB, KNUMB
        
        REAL(WP)   :: RHS11, RHS22, RHS33, RHS44, RHS55, RHS66
        

        DO J=1, N2DO(MYID)
            JJ = JCL2G(J)
            !=============================
            ENEJTF = 0.0_WP
            ENEJME = 0.0_WP
            EN3KI  = 0.0_WP
            EN1IK  = 0.0_WP
            
            !===========FLUCTUATION STREAMWISE VELOCITY(Q1)=========================
            RHS = 0.0_wp
            DO I=1,NCL1_io
                DO K=1,NCL3
                    RHS11=VARINST(I,J,K,1)-VARAVRG(J,1)   ! u'
                    RHS(I,J,K)=RHS11                      ! u'
                    ENEJTF(1,J)=ENEJTF(1,J)+RHS11**2      ! SIGMA_xz(u'*u')
                    ENEJME(1,J)=ENEJME(1,J)+RHS11         ! SIGMA_xz(u') ! zero??
                ENDDO
            ENDDO
            ENEJTF(1,J)=ENEJTF(1,J)*VL1313_io ! <u'u'>_xz
            ENEJME(1,J)=ENEJME(1,J)*VL1313_io ! <u'>_xz ? Zero??
            
            !-----SPANWISE   INVERSE FFT X3   Q1 VELOCITY
            CALL RUUX3(J,1) ! get EN1IK(1,K,I)
            !-----STREAMWISE INVERSE FFT X1   Q1 VELOCITY
            CALL RUUX1(J,1) ! get EN3KI(1,I,K)
            
            !===========FLUCTUATION NORMAL VELOCITY(Q2)===========================
            RHS = 0.0_wp
            DO K=1,NCL3
                DO I=1,NCL1_io
                    RHS22=VARINST(I,J,K,2)-VARAVRG(J,2) !v'
                    RHS(I,J,K)=RHS22                    !v'
                    ENEJTF(2,J)=ENEJTF(2,J)+RHS22**2    ! SIGMA_xz(v'*v')
                    ENEJME(2,J)=ENEJME(2,J)+RHS22       ! SIGMA_xz(u')
                ENDDO
            ENDDO
            ENEJTF(2,J)=ENEJTF(2,J)*VL1313_io           ! <v'v'>_xz
            ENEJME(2,J)=ENEJME(2,J)*VL1313_io           ! <v'>_xz
            
            !-----SPANWISE   INVERSE FFT X3  Q2 VELOCITY
            CALL RUUX3(J,2) ! get EN1IK(2,K,I)
            !-----STREAMWISE INVERSE FFT X1  Q2 VELOCITY
            CALL RUUX1(J,2) ! get EN3KI(2,I,K)
            
            !===========FLUCTUATION SPANWISE VELOCITY(Q3)===========================
            RHS = 0.0_wp
            DO I=1,NCL1_io
                DO K=1,NCL3
                    RHS33=VARINST(I,J,K,3)-VARAVRG(J,3) !w'
                    RHS(I,J,K)=RHS33
                    ENEJTF(3,J)=ENEJTF(3,J)+RHS33**2
                    ENEJME(3,J)=ENEJME(3,J)+RHS33
                ENDDO
            ENDDO
            ENEJTF(3,J)=ENEJTF(3,J)*VL1313_io
            ENEJME(3,J)=ENEJME(3,J)*VL1313_io
            
            !-----STREAMWISE INVERSE FFT X1  Q3 VELOCITY
            CALL RUUX1(J,3) ! get EN1IK
            !-----SPANWISE   INVERSE FFT X3  Q3 VELOCITY
            CALL RUUX3(J,3) ! get EN3KI
            
            
            !===========FLUCTUATION STREAMWISE VORTICITY(Omega1)=========================
            RHS = 0.0_wp
            DO I=1,NCL1_io
                DO K=1,NCL3
                    RHS44=VARINST(I,J,K,4)-VARAVRG(J,4)   ! omega1'
                    RHS(I,J,K)=RHS44                      ! omega1'
                    ENEJTF(4,J)=ENEJTF(4,J)+RHS44**2      ! SIGMA_xz(omega1'*omega1')
                    ENEJME(4,J)=ENEJME(4,J)+RHS44         ! SIGMA_xz(u') ! zero??
                ENDDO
            ENDDO
            ENEJTF(4,J)=ENEJTF(4,J)*VL1313_io             ! <omega1'omega1'>_xz
            ENEJME(4,J)=ENEJME(4,J)*VL1313_io             ! <omega1'>_xz 
            
            !-----SPANWISE   INVERSE FFT X3   Omega1 VORTICITY
            CALL RUUX3(J,4) ! get EN1IK 
            !-----STREAMWISE INVERSE FFT X1   Omega1 VORTICITY
            CALL RUUX1(J,4) ! get EN3KI
            
            !===========FLUCTUATION STREAMWISE VORTICITY(Omega2)=========================
            RHS = 0.0_wp
            DO I=1,NCL1_io
                DO K=1,NCL3
                    RHS55=VARINST(I,J,K,5)-VARAVRG(J,5)   ! omega2'
                    RHS(I,J,K)=RHS55                      ! omega2'
                    ENEJTF(5,J)=ENEJTF(5,J)+RHS55**2      ! SIGMA_xz(omega2'*omega2')
                    ENEJME(5,J)=ENEJME(5,J)+RHS55         ! SIGMA_xz(omega2') ! zero??
                ENDDO
            ENDDO
            ENEJTF(5,J)=ENEJTF(5,J)*VL1313_io ! <omega2'omega2'>_xz
            ENEJME(5,J)=ENEJME(5,J)*VL1313_io ! <omega2'>_xz ? Zero??
            
            !-----SPANWISE   INVERSE FFT X3   Omega1 VORTICITY
            CALL RUUX3(J,5) ! get EN1IK 
            !-----STREAMWISE INVERSE FFT X1   Omega1 VORTICITY
            CALL RUUX1(J,5) ! get EN3KI
            
            !===========FLUCTUATION STREAMWISE VORTICITY(Omega3)=========================
            RHS = 0.0_wp
            DO I=1,NCL1_io
                DO K=1,NCL3
                    RHS66=VARINST(I,J,K,6)-VARAVRG(J,6)   ! omega3'
                    RHS(I,J,K)=RHS66                      ! omega3'
                    ENEJTF(6,J)=ENEJTF(6,J)+RHS66**2      ! SIGMA_xz(omega3'*omega3')
                    ENEJME(6,J)=ENEJME(6,J)+RHS66         ! SIGMA_xz(omega13') !
                ENDDO
            ENDDO
            ENEJTF(6,J)=ENEJTF(6,J)*VL1313_io ! <u'u'>_xz
            ENEJME(6,J)=ENEJME(6,J)*VL1313_io ! <u'>_xz ? Zero??
            
            !-----SPANWISE   INVERSE FFT X3   Omega1 VORTICITY
            CALL RUUX3(J,6) ! get EN1IK 
            !-----STREAMWISE INVERSE FFT X1   Omega1 VORTICITY
            CALL RUUX1(J,6) ! get EN3KI
            
            !===================get max. ENE3MA (not used)=====================
            DO L=1,6
                ENE3MA(L)=0.0_wp
                ENE1MA(L)=0.0_wp
                DO K=1,N3MH-1
                    KP=2*K
                    KD=2*K-1
                    DO I=1,NCL1_io
                        ENE3MA(L)=DMAX1(DABS(EN3KI(L,KD,I)),ENE3MA(L))
                        ENE3MA(L)=DMAX1(DABS(EN3KI(L,KP,I)),ENE3MA(L))
                    ENDDO
                ENDDO
                
                DO I=1,N1MH-1
                    IP=2*I
                    ID=2*I-1
                    DO K=1,NCL3
                        ENE1MA(L)=MAX(ABS(EN1IK(L,ID,K)),ENE1MA(L))
                        ENE1MA(L)=MAX(ABS(EN1IK(L,IP,K)),ENE1MA(L))
                    ENDDO
                ENDDO
            ENDDO 
            
            ! test=====
            !WRITE(*,'(A,2I4.1,6ES13.5)') '#J, JJ, ENE1MA,ENE3MA', J, JJ, ENE1MA,ENE3MA
            !  
            !===========================================
        
            CALL COSPEC(J,JJ)
        
            !=====================================================
            DO II=1,NCL1_io
            
                EE1_11(JJ,II)=ENE11T(JJ,II)
                EE1_22(JJ,II)=ENE22T(JJ,II)
                EE1_33(JJ,II)=ENE33T(JJ,II)
                EE1_12(JJ,II)=ENE12T(JJ,II)
                EE1_13(JJ,II)=ENE13T(JJ,II)
                EE1_23(JJ,II)=ENE23T(JJ,II)
                
                UU1_11(JJ,II)=R11X1(JJ,II)
                UU1_22(JJ,II)=R22X1(JJ,II)
                UU1_33(JJ,II)=R33X1(JJ,II)
                UU1_12(JJ,II)=R12X1(JJ,II)
                UU1_13(JJ,II)=R13X1(JJ,II)
                UU1_23(JJ,II)=R23X1(JJ,II)
                
                HH1_11(JJ,II)=ENV11T(JJ,II)
                HH1_22(JJ,II)=ENV22T(JJ,II)
                HH1_33(JJ,II)=ENV33T(JJ,II)
                HH1_12(JJ,II)=ENV12T(JJ,II)
                HH1_13(JJ,II)=ENV13T(JJ,II)
                HH1_23(JJ,II)=ENV23T(JJ,II)
                
                VV1_11(JJ,II)=V11X1(JJ,II)
                VV1_22(JJ,II)=V22X1(JJ,II)
                VV1_33(JJ,II)=V33X1(JJ,II)
                VV1_12(JJ,II)=V12X1(JJ,II)
                VV1_13(JJ,II)=V13X1(JJ,II)
                VV1_23(JJ,II)=V23X1(JJ,II)
                
                EH1_11(JJ,II)=EVO11T(JJ,II)
                EH1_12(JJ,II)=EVO12T(JJ,II)
                EH1_13(JJ,II)=EVO13T(JJ,II)
                
                EH1_21(JJ,II)=EVO21T(JJ,II)
                EH1_22(JJ,II)=EVO22T(JJ,II)
                EH1_23(JJ,II)=EVO23T(JJ,II)
                
                EH1_31(JJ,II)=EVO31T(JJ,II)
                EH1_32(JJ,II)=EVO32T(JJ,II)
                EH1_33(JJ,II)=EVO33T(JJ,II)
                
                UV1_11(JJ,II)=VO11X1(JJ,II)
                UV1_12(JJ,II)=VO12X1(JJ,II)
                UV1_13(JJ,II)=VO13X1(JJ,II)
                
                UV1_21(JJ,II)=VO21X1(JJ,II)
                UV1_22(JJ,II)=VO22X1(JJ,II)
                UV1_23(JJ,II)=VO23X1(JJ,II)
                
                UV1_31(JJ,II)=VO31X1(JJ,II)
                UV1_32(JJ,II)=VO32X1(JJ,II)
                UV1_33(JJ,II)=VO33X1(JJ,II)
                
                
            END DO
            
            DO KK=1,NCL3
                EE3_11(JJ,KK)=ENE11Z(JJ,KK)
                EE3_22(JJ,KK)=ENE22Z(JJ,KK)
                EE3_33(JJ,KK)=ENE33Z(JJ,KK)
                EE3_12(JJ,KK)=ENE12Z(JJ,KK)
                EE3_13(JJ,KK)=ENE13Z(JJ,KK)
                EE3_23(JJ,KK)=ENE23Z(JJ,KK)
                
                UU3_11(JJ,KK)=R11X3(JJ,KK)
                UU3_22(JJ,KK)=R22X3(JJ,KK)
                UU3_33(JJ,KK)=R33X3(JJ,KK)
                UU3_12(JJ,KK)=R12X3(JJ,KK)
                UU3_13(JJ,KK)=R13X3(JJ,KK)
                UU3_23(JJ,KK)=R23X3(JJ,KK)
                
                HH3_11(JJ,KK)=ENV11Z(JJ,KK)
                HH3_22(JJ,KK)=ENV22Z(JJ,KK)
                HH3_33(JJ,KK)=ENV33Z(JJ,KK)
                HH3_12(JJ,KK)=ENV12Z(JJ,KK)
                HH3_13(JJ,KK)=ENV13Z(JJ,KK)
                HH3_23(JJ,KK)=ENV23Z(JJ,KK)
                
                VV3_11(JJ,KK)=V11X3(JJ,KK)
                VV3_22(JJ,KK)=V22X3(JJ,KK)
                VV3_33(JJ,KK)=V33X3(JJ,KK)
                VV3_12(JJ,KK)=V12X3(JJ,KK)
                VV3_13(JJ,KK)=V13X3(JJ,KK)
                VV3_23(JJ,KK)=V23X3(JJ,KK)
                
                EH3_11(JJ,KK)=EVO11Z(JJ,KK)
                EH3_12(JJ,KK)=EVO12Z(JJ,KK)
                EH3_13(JJ,KK)=EVO13Z(JJ,KK)
                
                EH3_21(JJ,KK)=EVO21Z(JJ,KK)
                EH3_22(JJ,KK)=EVO22Z(JJ,KK)
                EH3_23(JJ,KK)=EVO23Z(JJ,KK)
                
                EH3_31(JJ,KK)=EVO31Z(JJ,KK)
                EH3_32(JJ,KK)=EVO32Z(JJ,KK)
                EH3_33(JJ,KK)=EVO33Z(JJ,KK)
                
                UV3_11(JJ,KK)=VO11X3(JJ,KK)
                UV3_12(JJ,KK)=VO12X3(JJ,KK)
                UV3_13(JJ,KK)=VO13X3(JJ,KK)
                
                UV3_21(JJ,KK)=VO21X3(JJ,KK)
                UV3_22(JJ,KK)=VO22X3(JJ,KK)
                UV3_23(JJ,KK)=VO23X3(JJ,KK)
                
                UV3_31(JJ,KK)=VO31X3(JJ,KK)
                UV3_32(JJ,KK)=VO32X3(JJ,KK)
                UV3_33(JJ,KK)=VO33X3(JJ,KK)
            END DO
        END DO
        

        
        !===================================================================================================
        CALL MPI_BARRIER(ICOMM,IERROR)
        INUMB=NCL2*NCL1_io
        CALL MPI_ALLREDUCE(EE1_11(1,1), WORK_EE1_11(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE1_22(1,1), WORK_EE1_22(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE1_33(1,1), WORK_EE1_33(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE1_12(1,1), WORK_EE1_12(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE1_13(1,1), WORK_EE1_13(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE1_23(1,1), WORK_EE1_23(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UU1_11(1,1), WORK_UU1_11(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU1_22(1,1), WORK_UU1_22(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU1_33(1,1), WORK_UU1_33(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU1_12(1,1), WORK_UU1_12(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU1_13(1,1), WORK_UU1_13(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU1_23(1,1), WORK_UU1_23(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(HH1_11(1,1), WORK_HH1_11(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH1_22(1,1), WORK_HH1_22(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH1_33(1,1), WORK_HH1_33(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH1_12(1,1), WORK_HH1_12(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH1_13(1,1), WORK_HH1_13(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH1_23(1,1), WORK_HH1_23(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(VV1_11(1,1), WORK_VV1_11(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV1_22(1,1), WORK_VV1_22(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV1_33(1,1), WORK_VV1_33(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV1_12(1,1), WORK_VV1_12(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV1_13(1,1), WORK_VV1_13(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV1_23(1,1), WORK_VV1_23(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(EH1_11(1,1), WORK_EH1_11(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH1_12(1,1), WORK_EH1_12(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH1_13(1,1), WORK_EH1_13(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(EH1_21(1,1), WORK_EH1_21(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH1_22(1,1), WORK_EH1_22(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH1_23(1,1), WORK_EH1_23(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(EH1_31(1,1), WORK_EH1_31(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH1_32(1,1), WORK_EH1_32(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH1_33(1,1), WORK_EH1_33(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UV1_11(1,1), WORK_UV1_11(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV1_12(1,1), WORK_UV1_12(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV1_13(1,1), WORK_UV1_13(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UV1_21(1,1), WORK_UV1_21(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV1_22(1,1), WORK_UV1_22(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV1_23(1,1), WORK_UV1_23(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UV1_31(1,1), WORK_UV1_31(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV1_32(1,1), WORK_UV1_32(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV1_33(1,1), WORK_UV1_33(1,1), INUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        !===================================================================================================
        CALL MPI_BARRIER(ICOMM,IERROR)
        KNUMB=NCL2*NCL3
        CALL MPI_ALLREDUCE(EE3_11(1,1), WORK_EE3_11(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE3_22(1,1), WORK_EE3_22(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE3_33(1,1), WORK_EE3_33(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE3_12(1,1), WORK_EE3_12(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE3_13(1,1), WORK_EE3_13(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EE3_23(1,1), WORK_EE3_23(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UU3_11(1,1), WORK_UU3_11(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU3_22(1,1), WORK_UU3_22(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU3_33(1,1), WORK_UU3_33(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU3_12(1,1), WORK_UU3_12(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU3_13(1,1), WORK_UU3_13(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UU3_23(1,1), WORK_UU3_23(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(HH3_11(1,1), WORK_HH3_11(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH3_22(1,1), WORK_HH3_22(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH3_33(1,1), WORK_HH3_33(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH3_12(1,1), WORK_HH3_12(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH3_13(1,1), WORK_HH3_13(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(HH3_23(1,1), WORK_HH3_23(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(VV3_11(1,1), WORK_VV3_11(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV3_22(1,1), WORK_VV3_22(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV3_33(1,1), WORK_VV3_33(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV3_12(1,1), WORK_VV3_12(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV3_13(1,1), WORK_VV3_13(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(VV3_23(1,1), WORK_VV3_23(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(EH3_11(1,1), WORK_EH3_11(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH3_12(1,1), WORK_EH3_12(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH3_13(1,1), WORK_EH3_13(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(EH3_21(1,1), WORK_EH3_21(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH3_22(1,1), WORK_EH3_22(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH3_23(1,1), WORK_EH3_23(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(EH3_31(1,1), WORK_EH3_31(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH3_32(1,1), WORK_EH3_32(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(EH3_33(1,1), WORK_EH3_33(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UV3_11(1,1), WORK_UV3_11(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV3_12(1,1), WORK_UV3_12(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV3_13(1,1), WORK_UV3_13(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UV3_21(1,1), WORK_UV3_21(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV3_22(1,1), WORK_UV3_22(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV3_23(1,1), WORK_UV3_23(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        CALL MPI_ALLREDUCE(UV3_31(1,1), WORK_UV3_31(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV3_32(1,1), WORK_UV3_32(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        CALL MPI_ALLREDUCE(UV3_33(1,1), WORK_UV3_33(1,1), KNUMB, MPI_DOUBLE_PRECISION, MPI_SUM, ICOMM, IERROR)
        
        IF (MYID.EQ.0) THEN
            DO JJ=1,NCL2
                DO IC=1,NCL1_io
                    ENE11T(JJ,IC) = WORK_EE1_11(JJ,IC)
                    ENE22T(JJ,IC) = WORK_EE1_22(JJ,IC)
                    ENE33T(JJ,IC) = WORK_EE1_33(JJ,IC)
                    ENE12T(JJ,IC) = WORK_EE1_12(JJ,IC)
                    ENE13T(JJ,IC) = WORK_EE1_13(JJ,IC)
                    ENE23T(JJ,IC) = WORK_EE1_23(JJ,IC)
                     
                    R11X1(JJ,IC)  = WORK_UU1_11(JJ,IC)
                    R22X1(JJ,IC)  = WORK_UU1_22(JJ,IC)
                    R33X1(JJ,IC)  = WORK_UU1_33(JJ,IC)
                    R12X1(JJ,IC)  = WORK_UU1_12(JJ,IC)
                    R13X1(JJ,IC)  = WORK_UU1_13(JJ,IC)
                    R23X1(JJ,IC)  = WORK_UU1_23(JJ,IC)
                     
                    ENV11T(JJ,IC) = WORK_HH1_11(JJ,IC)
                    ENV22T(JJ,IC) = WORK_HH1_22(JJ,IC)
                    ENV33T(JJ,IC) = WORK_HH1_33(JJ,IC)
                    ENV12T(JJ,IC) = WORK_HH1_12(JJ,IC)
                    ENV13T(JJ,IC) = WORK_HH1_13(JJ,IC)
                    ENV23T(JJ,IC) = WORK_HH1_23(JJ,IC)
                     
                    V11X1(JJ,IC)  = WORK_VV1_11(JJ,IC)
                    V22X1(JJ,IC)  = WORK_VV1_22(JJ,IC)
                    V33X1(JJ,IC)  = WORK_VV1_33(JJ,IC)
                    V12X1(JJ,IC)  = WORK_VV1_12(JJ,IC)
                    V13X1(JJ,IC)  = WORK_VV1_13(JJ,IC)
                    V23X1(JJ,IC)  = WORK_VV1_23(JJ,IC)
                     
                    EVO11T(JJ,IC) = WORK_EH1_11(JJ,IC)
                    EVO12T(JJ,IC) = WORK_EH1_12(JJ,IC)
                    EVO13T(JJ,IC) = WORK_EH1_13(JJ,IC)
                     
                    EVO21T(JJ,IC) = WORK_EH1_21(JJ,IC)
                    EVO22T(JJ,IC) = WORK_EH1_22(JJ,IC)
                    EVO23T(JJ,IC) = WORK_EH1_23(JJ,IC)
                     
                    EVO31T(JJ,IC) = WORK_EH1_31(JJ,IC)
                    EVO32T(JJ,IC) = WORK_EH1_32(JJ,IC)
                    EVO33T(JJ,IC) = WORK_EH1_33(JJ,IC)
                     
                    VO11X1(JJ,IC) = WORK_UV1_11(JJ,IC)
                    VO12X1(JJ,IC) = WORK_UV1_12(JJ,IC)
                    VO13X1(JJ,IC) = WORK_UV1_13(JJ,IC)
                     
                    VO21X1(JJ,IC) = WORK_UV1_21(JJ,IC)
                    VO22X1(JJ,IC) = WORK_UV1_22(JJ,IC)
                    VO23X1(JJ,IC) = WORK_UV1_23(JJ,IC)
                     
                    VO31X1(JJ,IC) = WORK_UV1_31(JJ,IC)
                    VO32X1(JJ,IC) = WORK_UV1_32(JJ,IC)
                    VO33X1(JJ,IC) = WORK_UV1_33(JJ,IC)
                ENDDO

                DO KC=1,NCL3
                    ENE11Z(JJ,KC) = WORK_EE3_11(JJ,KC)
                    ENE22Z(JJ,KC) = WORK_EE3_22(JJ,KC)
                    ENE33Z(JJ,KC) = WORK_EE3_33(JJ,KC)
                    ENE12Z(JJ,KC) = WORK_EE3_12(JJ,KC)
                    ENE13Z(JJ,KC) = WORK_EE3_13(JJ,KC)
                    ENE23Z(JJ,KC) = WORK_EE3_23(JJ,KC)
                     
                    R11X3(JJ,KC)  = WORK_UU3_11(JJ,KC)
                    R22X3(JJ,KC)  = WORK_UU3_22(JJ,KC)
                    R33X3(JJ,KC)  = WORK_UU3_33(JJ,KC)
                    R12X3(JJ,KC)  = WORK_UU3_12(JJ,KC)
                    R13X3(JJ,KC)  = WORK_UU3_13(JJ,KC)
                    R23X3(JJ,KC)  = WORK_UU3_23(JJ,KC)
                     
                    ENV11Z(JJ,KC) = WORK_HH3_11(JJ,KC)
                    ENV22Z(JJ,KC) = WORK_HH3_22(JJ,KC)
                    ENV33Z(JJ,KC) = WORK_HH3_33(JJ,KC)
                    ENV12Z(JJ,KC) = WORK_HH3_12(JJ,KC)
                    ENV13Z(JJ,KC) = WORK_HH3_13(JJ,KC)
                    ENV23Z(JJ,KC) = WORK_HH3_23(JJ,KC)
                     
                    V11X3(JJ,KC)  = WORK_VV3_11(JJ,KC)
                    V22X3(JJ,KC)  = WORK_VV3_22(JJ,KC)
                    V33X3(JJ,KC)  = WORK_VV3_33(JJ,KC)
                    V12X3(JJ,KC)  = WORK_VV3_12(JJ,KC)
                    V13X3(JJ,KC)  = WORK_VV3_13(JJ,KC)
                    V23X3(JJ,KC)  = WORK_VV3_23(JJ,KC)
                     
                    EVO11Z(JJ,KC) = WORK_EH3_11(JJ,KC)
                    EVO12Z(JJ,KC) = WORK_EH3_12(JJ,KC)
                    EVO13Z(JJ,KC) = WORK_EH3_13(JJ,KC)
                     
                    EVO21Z(JJ,KC) = WORK_EH3_21(JJ,KC)
                    EVO22Z(JJ,KC) = WORK_EH3_22(JJ,KC)
                    EVO23Z(JJ,KC) = WORK_EH3_23(JJ,KC)
                     
                    EVO31Z(JJ,KC) = WORK_EH3_31(JJ,KC)
                    EVO32Z(JJ,KC) = WORK_EH3_32(JJ,KC)
                    EVO33Z(JJ,KC) = WORK_EH3_33(JJ,KC)
                     
                    VO11X3(JJ,KC) = WORK_UV3_11(JJ,KC)
                    VO12X3(JJ,KC) = WORK_UV3_12(JJ,KC)
                    VO13X3(JJ,KC) = WORK_UV3_13(JJ,KC)
                     
                    VO21X3(JJ,KC) = WORK_UV3_21(JJ,KC)
                    VO22X3(JJ,KC) = WORK_UV3_22(JJ,KC)
                    VO23X3(JJ,KC) = WORK_UV3_23(JJ,KC)
                     
                    VO31X3(JJ,KC) = WORK_UV3_31(JJ,KC)
                    VO32X3(JJ,KC) = WORK_UV3_32(JJ,KC)
                    VO33X3(JJ,KC) = WORK_UV3_33(JJ,KC)
                ENDDO
            ENDDO
        END IF
        
        RETURN                                                            
    END SUBROUTINE
      
!======================================================================================   
    ! input :               EN1IK, EN3KI
    ! intermediate output:  ENE1, ENE3, CORX1, CORX3
    ! output:               ENE?ZO, R??Z?
    SUBROUTINE COSPEC(J,JWORK)
        use mesh_info
        use SPECO_info
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: J
        INTEGER(4),INTENT(IN) :: JWORK
        INTEGER(4) :: K
        
        !===========================================================================
        !=========SPECTRA AND CORRELATIONS  VELOCITY 11 COMPONENT-----------------
        CALL COSPX3(J,JWORK,1,1)
        DO K=1,N3MH ! check half or total?
            ENE11Z(JWORK,K)=ENE11Z(JWORK,K)+ENE3 (JWORK,K)
            R11X3 (JWORK,K)=R11X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,1,1)
        DO K=1,N1MH
            ENE11T(JWORK,K)=ENE11T(JWORK,K)+ENE1 (JWORK,K)
            R11X1 (JWORK,K)=R11X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO

        !=========SPECTRA AND CORRELATIONS  VELOCITY 22 COMPONENT-----------------
        CALL COSPX3(J,JWORK,2,2)
        !WRITE(6,*)'IN COSPEC 2 ',J,R22X3(J,1),CORX3(J,1)
        DO K=1,N3MH
            ENE22Z(JWORK,K)=ENE22Z(JWORK,K)+ENE3 (JWORK,K)
            R22X3 (JWORK,K)=R22X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,2,2)
        DO K=1,N1MH
            ENE22T(JWORK,K)=ENE22T(JWORK,K)+ENE1 (JWORK,K)
            R22X1 (JWORK,K)=R22X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO

        !=========SPECTRA AND CORRELATIONS  VELOCITY 33 COMPONENT-----------------
        CALL COSPX3(J,JWORK,3,3)
        DO K=1,N3MH
            ENE33Z(JWORK,K)=ENE33Z(JWORK,K)+ENE3 (JWORK,K)
            R33X3 (JWORK,K)=R33X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,3,3)
        DO K=1,N1MH
            ENE33T(JWORK,K)=ENE33T(JWORK,K)+ENE1 (JWORK,K)
            R33X1 (JWORK,K)=R33X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========COSPECTRA AND CORRELATIONS  12 VELOCITY COMPONENT--------------
        CALL COSPX3(J,JWORK,1,2)
        DO K=1,N3MH
            ENE12Z(JWORK,K)=ENE12Z(JWORK,K)+ENE3 (JWORK,K)
            R12X3 (JWORK,K)=R12X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,1,2)
        DO K=1,N1MH
            ENE12T(JWORK,K)=ENE12T(JWORK,K)+ENE1 (JWORK,K)
            R12X1 (JWORK,K)=R12X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO

        !=========COSPECTRA AND CORRELATIONS  13 VELOCITY COMPONENT--------------
        CALL COSPX3(J,JWORK,1,3)
        DO K=1,N3MH
            ENE13Z(JWORK,K)=ENE13Z(JWORK,K)+ENE3 (JWORK,K)
            R13X3 (JWORK,K)=R13X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,1,3)
        DO K=1,N1MH
            ENE13T(JWORK,K)=ENE13T(JWORK,K)+ENE1 (JWORK,K)
            R13X1 (JWORK,K)=R13X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO

        !=========COSPECTRA AND CORRELATIONS  23 VELOCITY COMPONENT--------------
        CALL COSPX3(J,JWORK,2,3)
        DO K=1,N3MH
            ENE23Z(JWORK,K)=ENE23Z(JWORK,K)+ENE3 (JWORK,K)
            R23X3 (JWORK,K)=R23X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,2,3)
        DO K=1,N1MH
            ENE23T(JWORK,K)=ENE23T(JWORK,K)+ENE1 (JWORK,K)
            R23X1 (JWORK,K)=R23X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        
        !===========================================================================
        !=========SPECTRA AND CORRELATIONS  VORTICITY 11 COMPONENT-----------------
        CALL COSPX3(J,JWORK,4,4)
        DO K=1,N3MH
            ENV11Z(JWORK,K)=ENV11Z(JWORK,K)+ENE3 (JWORK,K)
            V11X3 (JWORK,K)=V11X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,4,4)
        DO K=1,N1MH
            ENV11T(JWORK,K)=ENV11T(JWORK,K)+ENE1 (JWORK,K)
            V11X1 (JWORK,K)=V11X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VORTICITY 22 COMPONENT-----------------
        CALL COSPX3(J,JWORK,5,5)
        DO K=1,N3MH
            ENV22Z(JWORK,K)=ENV22Z(JWORK,K)+ENE3 (JWORK,K)
            V22X3 (JWORK,K)=V22X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,5,5)
        DO K=1,N1MH
            ENV22T(JWORK,K)=ENV22T(JWORK,K)+ENE1 (JWORK,K)
            V22X1 (JWORK,K)=V22X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VORTICITY 33 COMPONENT-----------------
        CALL COSPX3(J,JWORK,6,6)
        DO K=1,N3MH
            ENV33Z(JWORK,K)=ENV33Z(JWORK,K)+ENE3 (JWORK,K)
            V33X3 (JWORK,K)=V33X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,6,6)
        DO K=1,N1MH
            ENV33T(JWORK,K)=ENV33T(JWORK,K)+ENE1 (JWORK,K)
            V33X1 (JWORK,K)=V33X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO


        !=========SPECTRA AND CORRELATIONS  VORTICITY 12 COMPONENT-----------------
        CALL COSPX3(J,JWORK,4,5)
        DO K=1,N3MH
            ENV12Z(JWORK,K)=ENV12Z(JWORK,K)+ENE3 (JWORK,K)
            V12X3 (JWORK,K)=V12X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,4,5)
        DO K=1,N1MH
            ENV12T(JWORK,K)=ENV12T(JWORK,K)+ENE1 (JWORK,K)
            V12X1 (JWORK,K)=V12X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VORTICITY 13 COMPONENT-----------------
        CALL COSPX3(J,JWORK,4,6)
        DO K=1,N3MH
            ENV13Z(JWORK,K)=ENV13Z(JWORK,K)+ENE3 (JWORK,K)
            V13X3 (JWORK,K)=V13X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,4,6)
        DO K=1,N1MH
            ENV13T(JWORK,K)=ENV13T(JWORK,K)+ENE1 (JWORK,K)
            V13X1 (JWORK,K)=V13X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        
        !=========SPECTRA AND CORRELATIONS  VORTICITY 23 COMPONENT-----------------
        CALL COSPX3(J,JWORK,5,6)
        DO K=1,N3MH
            ENV23Z(JWORK,K)=ENV23Z(JWORK,K)+ENE3 (JWORK,K)
            V23X3 (JWORK,K)=V23X3 (JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,5,6)
        DO K=1,N1MH
            ENV23T(JWORK,K)=ENV23T(JWORK,K)+ENE1 (JWORK,K)
            V23X1 (JWORK,K)=V23X1 (JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !===========================================================================
        !=========SPECTRA AND CORRELATIONS  VELOCITY 1 & VORTICITY 1-----------------
        CALL COSPX3(J,JWORK,1,4)
        DO K=1,N3MH
            EVO11Z(JWORK,K)=EVO11Z(JWORK,K)+ENE3 (JWORK,K)
            VO11X3(JWORK,K)=VO11X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,1,4)
        DO K=1,N1MH
            EVO11T(JWORK,K)=EVO11T(JWORK,K)+ENE1 (JWORK,K)
            VO11X1(JWORK,K)=VO11X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 1 & VORTICITY 2-----------------
        CALL COSPX3(J,JWORK,1,5)
        DO K=1,N3MH
            EVO12Z(JWORK,K)=EVO12Z(JWORK,K)+ENE3 (JWORK,K)
            VO12X3(JWORK,K)=VO12X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,1,5)
        DO K=1,N1MH
            EVO12T(JWORK,K)=EVO12T(JWORK,K)+ENE1 (JWORK,K)
            VO12X1(JWORK,K)=VO12X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 1 & VORTICITY 3-----------------
        CALL COSPX3(J,JWORK,1,6)
        DO K=1,N3MH
            EVO13Z(JWORK,K)=EVO13Z(JWORK,K)+ENE3 (JWORK,K)
            VO13X3(JWORK,K)=VO13X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,1,6)
        DO K=1,N1MH
            EVO13T(JWORK,K)=EVO13T(JWORK,K)+ENE1 (JWORK,K)
            VO13X1(JWORK,K)=VO13X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 2 & VORTICITY 1-----------------
        CALL COSPX3(J,JWORK,2,4)
        DO K=1,N3MH
            EVO21Z(JWORK,K)=EVO21Z(JWORK,K)+ENE3 (JWORK,K)
            VO21X3(JWORK,K)=VO21X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,2,4)
        DO K=1,N1MH
            EVO21T(JWORK,K)=EVO21T(JWORK,K)+ENE1 (JWORK,K)
            VO21X1(JWORK,K)=VO21X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 2 & VORTICITY 2-----------------
        CALL COSPX3(J,JWORK,2,5)
        DO K=1,N3MH
            EVO22Z(JWORK,K)=EVO22Z(JWORK,K)+ENE3 (JWORK,K)
            VO22X3(JWORK,K)=VO22X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,2,5)
        DO K=1,N1MH
            EVO22T(JWORK,K)=EVO22T(JWORK,K)+ENE1 (JWORK,K)
            VO22X1(JWORK,K)=VO22X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 2 & VORTICITY 3-----------------
        CALL COSPX3(J,JWORK,2,6)
        DO K=1,N3MH
            EVO23Z(JWORK,K)=EVO23Z(JWORK,K)+ENE3 (JWORK,K)
            VO23X3(JWORK,K)=VO23X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,2,6)
        DO K=1,N1MH
            EVO23T(JWORK,K)=EVO23T(JWORK,K)+ENE1 (JWORK,K)
            VO23X1(JWORK,K)=VO23X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 3 & VORTICITY 1-----------------
        CALL COSPX3(J,JWORK,3,4)
        DO K=1,N3MH
            EVO31Z(JWORK,K)=EVO31Z(JWORK,K)+ENE3 (JWORK,K)
            VO31X3(JWORK,K)=VO31X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,3,4)
        DO K=1,N1MH
            EVO31T(JWORK,K)=EVO31T(JWORK,K)+ENE1 (JWORK,K)
            VO31X1(JWORK,K)=VO31X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 3 & VORTICITY 2-----------------
        CALL COSPX3(J,JWORK,3,5)
        DO K=1,N3MH
            EVO32Z(JWORK,K)=EVO32Z(JWORK,K)+ENE3 (JWORK,K)
            VO32X3(JWORK,K)=VO32X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,3,5)
        DO K=1,N1MH
            EVO32T(JWORK,K)=EVO32T(JWORK,K)+ENE1 (JWORK,K)
            VO32X1(JWORK,K)=VO32X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        !=========SPECTRA AND CORRELATIONS  VELOCITY 3 & VORTICITY 3-----------------
        CALL COSPX3(J,JWORK,3,6)
        DO K=1,N3MH
            EVO33Z(JWORK,K)=EVO33Z(JWORK,K)+ENE3 (JWORK,K)
            VO33X3(JWORK,K)=VO33X3(JWORK,K)+CORX3(JWORK,K)
        ENDDO
        
        CALL COSPX1(J,JWORK,3,6)
        DO K=1,N1MH
            EVO33T(JWORK,K)=EVO33T(JWORK,K)+ENE1 (JWORK,K)
            VO33X1(JWORK,K)=VO33X1(JWORK,K)+CORX1(JWORK,K)
        ENDDO
        
        
        RETURN 
    END SUBROUTINE
    
    
    
!****************************************************************************
! ***************************************************************************
! *     SBR. RUUX1  :INVERSE FFT IN X1 OF A GENERAL QUANTITY RHS IT IS      *
! *                       STORED IN EN1IK THAT IS USED  IN THE ROUTINE COSPX1          *                                                                                          *
! INPUT:  RHS  (1:NCL1, J,      1:NCL3 )
! OUTPUT: EN1IK(L,      1:NCL1, 1:NCL3 ) WAVE NUMBER SPACE
!****************************************************************************
    SUBROUTINE RUUX1(J,L1)
        use mesh_info
        use SPECO_info
        USE FFT99_info1
        IMPLICIT NONE

        INTEGER(4),INTENT(IN) :: J, L1
        
        INTEGER(4) :: K, I, IS, IP, ID
        REAL(WP)   :: XR1  (N1MD,NCL3)
        REAL(WP)   :: WORK1(N1MD,NCL3)
        
        DO K=1,NCL3
            XR1(1,K)=RHS(NCL1_io,J,K)
            DO I=1,NCL1_io
                IS=I+1
                XR1(IS,K)=RHS(I,J,K)
            ENDDO
            XR1(N1MD,K)=RHS(1,J,K)
        ENDDO
        
        !C   2D REAL  FFT APPLIED TO THE RHS BY FFT99 ALONG X1
        !C   FROM PHYSICAL TO WAVE NUMBER SPACE
        CALL FFT99(XR1,WORK1,TRIGXX1,IFXX1,1,N1MD,NCL1_io,NCL3,-1)
        
        DO I=1,N1MH
            IP=2*I
            ID=2*I-1
            DO K=1,NCL3
                EN1IK(L1,ID,K)=XR1(ID,K)
                EN1IK(L1,IP,K)=XR1(IP,K)
            ENDDO
        ENDDO
        
        RETURN
    END SUBROUTINE
    
    
!****************************************************************************
! ***************************************************************************
! *     SBR. RUUX3  :INVERSE FFT IN X3 OF A GENERAL QUANTITY RHS IT IS      *
! *                       STORED IN EN3KI THAT IS USED  IN THE ROUTINE COSPX3          *                                                                                          *
! INPUT:  RHS  (1:NCL1, J,      1:NCL3 )
! OUTPUT: EN3KI(L,      1:NCL3, 1:NCL1 ) !WAVE NUMBER SPACE
!****************************************************************************
    SUBROUTINE RUUX3(J,L1)
        use mesh_info
        use SPECO_info
        USE FFT99_info1
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: J
        INTEGER(4),INTENT(IN) :: L1
        
        INTEGER(4) :: I, K, KS, KP, KD
        REAL(WP)   :: XR1  (N3MD,NCL1_io)
        REAL(WP)   :: WORK1(N3MD,NCL1_io)
        
        !====================put RHS to XR============
        DO I=1,NCL1_io
            XR1(1,I)=RHS(I,J,NCL3)
            DO K=1,NCL3
                KS=K+1
                XR1(KS,I)=RHS(I,J,K)
            ENDDO
            XR1(N3MD,I)=RHS(I,J,1)
        ENDDO
        
        !==============2D REAL  FFT APPLIED TO THE RHS BY FFT99 ALONG X3
        !==============FROM PHYSICAL TO WAVE NUMBER SPACE
        CALL FFT99(XR1,WORK1,TRIGXX3,IFXX3,1,N3MD,NCL3,NCL1_io,-1)
        
        !==============
        DO K=1,N3MH
            KP=2*K
            KD=2*K-1
            DO I=1,NCL1_io
                EN3KI(L1,KD,I)=XR1(KD,I)
                EN3KI(L1,KP,I)=XR1(KP,I)
            ENDDO
        ENDDO
        
        RETURN
    END SUBROUTINE
      
      
      
!C****************************************************************************
!C ***************************************************************************
!C *     SBR. COSPX3  :COSPECTRA  IN THE SPANWISE DIRECTION                  *
!C *                                                                                                          *
!C ***************************************************************************
! INPUT : EN3KI(L,1:NCL3,1:NCL1) !WAVE NUMBER SPACE
! OUTPUT: ENE3 (J, 1:L3MH)       !WAVE NUMBER SPACE
!         CORX3(J, 1:NCL3)       !PHYSICAL SPACE
!C****************************************************************************
    SUBROUTINE COSPX3(J,JWORK,L1,M1)
        use mesh_info
        use SPECO_info
        USE FFT99_info1
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: J
        INTEGER(4),INTENT(IN) :: JWORK
        INTEGER(4),INTENT(IN) :: L1,M1
        
        REAL(WP)  :: CORR3(N3MD,1) !!WAVE NUMBER SPACE
        REAL(WP)  :: WORK1(N3MD,1)
        INTEGER(4) :: KE, K, KP, KD, KS, I
        REAL(WP)   :: ENEL
      
        
        CORR3(:,:)=0.0_WP
        ENE3(JWORK,:) = 0.0_WP
        
        DO K=2,N3MH
            KP=2*K
            KD=2*K-1
            KE=K
            DO I=1,NCL1_io
                ENEL=(EN3KI(L1,KD,I)*EN3KI(M1,KD,I) &
                     +EN3KI(L1,KP,I)*EN3KI(M1,KP,I))/DBLE(NCL1_IO)*2.0_WP
                ENE3(JWORK,KE)=ENE3(JWORK,KE)+ENEL  ! output
                CORR3(KD,1)=CORR3(KD,1)+ENEL   
                CORR3(KP,1)=0.0_WP
            ENDDO
        ENDDO
        
        K=1
        KP=2*K
        KD=2*K-1
        KE=K
        DO I=1,NCL1_IO
            ENEL=(EN3KI(L1,KD,I)*EN3KI(M1,KD,I) &
                 +EN3KI(L1,KP,I)*EN3KI(M1,KP,I))/DBLE(NCL1_IO)
            ENE3(JWORK,KE)=ENE3(JWORK,KE)+ENEL ! output
            CORR3(KD,1)=CORR3(KD,1)+ENEL
            CORR3(KP,1)=0.0_WP
        ENDDO
!C
!C   FFT FROM THE WAVE NUMBER TO PHYSICAL SPACE
!C
        CALL FFT99(CORR3,WORK1,TRIGXX3,IFXX3,1,N3MD,NCL3,1,+1)
        DO K=1,NCL3
            KS=K+1
            CORX3(JWORK,K)=CORR3(KS,1)/CORR3(2,1) ! output
        ENDDO
        
        RETURN
      END
!C****************************************************************************
!C ***************************************************************************
!C *     SBR. COSPX1  :COSPECTRA  IN THE STREAMWISE DIRECTION                *
!C *                                                                                                          *
!C ***************************************************************************
! INPUT : EN1IK(L,1:NCL1,1:NCL3)  !!WAVE NUMBER SPACE
! OUTPUT: ENE1 (J, 1:N1MH)        !!WAVE NUMBER SPACE
!         CORX1(J, 1:NCL1)        !!PHYSICAL SPACE
!C****************************************************************************
    SUBROUTINE COSPX1(J,JWORK,L1,M1)
        use mesh_info
        use SPECO_info
        USE FFT99_info1
        IMPLICIT NONE
        
        INTEGER(4),INTENT(IN) :: J
        INTEGER(4),INTENT(IN) :: JWORK
        INTEGER(4),INTENT(IN) :: L1,M1
        
        REAL(WP)  :: CORR1(N1MD,1)
        REAL(WP)  :: WORK1(N1MD,1)
        INTEGER(4) :: I, KE, IP, ID, IS, K
        REAL(WP)   :: ENEL
        
        

        DO KE=1,N1MD
            CORR1(KE,1)=0.0_wp
        ENDDO
        DO KE=1,NCL1_io
            ENE1(JWORK,KE)=0.0_wp
        ENDDO
        
        DO I=2,N1MH
            IP=2*I
            ID=2*I-1
            KE=I
            DO K=1,NCL3
                ENEL=(EN1IK(L1,ID,K)*EN1IK(M1,ID,K) &
                     +EN1IK(L1,IP,K)*EN1IK(M1,IP,K))/DBLE(NCL3)*2.0_wp
                ENE1(JWORK,KE)=ENE1(JWORK,KE)+ENEL !output
                CORR1(ID,1)=CORR1(ID,1)+ENEL
                CORR1(IP,1)=0.0_wp
            ENDDO
        ENDDO
        
        I=1
        IP=2*I
        ID=2*I-1
        KE=I
        DO K=1,NCL3
            ENEL=(EN1IK(L1,ID,K)*EN1IK(M1,ID,K) &
                 +EN1IK(L1,IP,K)*EN1IK(M1,IP,K))/ DBLE(NCL3)
            ENE1(JWORK,KE)=ENE1(JWORK,KE)+ENEL !output
            CORR1(ID,1)=CORR1(ID,1)+ENEL
            CORR1(IP,1)=0.0_wp
        ENDDO
        
        CALL FFT99(CORR1,WORK1,TRIGXX1,IFXX1,1,N1MD,NCL1_io,1,+1)
        
        DO I=1,NCL1_io
            IS=I+1
            CORX1(JWORK,I)=CORR1(IS,1)/CORR1(2,1) !output
        ENDDO
        RETURN
    END SUBROUTINE

