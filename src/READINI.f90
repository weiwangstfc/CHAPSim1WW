!***********************************************************************
!> @author 
!> Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!> Read in initial configuration file
!
!> @todo
!> Nothing left to do.
!
! REVISION HISTORY:
! 06/12/2013- Initial Version, by Wei Wang
!**********************************************************************
    SUBROUTINE READINI
!>     @note
!>     Read in the flow setting file. *.ini. 
!>     only in master

        use init_info
        use mpi_info
        use mesh_info
        use thermal_info
        IMPLICIT NONE
       
        CHARACTER(128) :: sect          ! section names
        CHARACTER(128) :: skey          ! key names
        INTEGER(4) :: IOS=0
        INTEGER(4) :: INI=13             !file I/O id
        INTEGER(4) :: lens   
        INTEGER(4) :: itmp
        INTEGER(4) :: N
        real(wp)   :: rtmp
        character(256) :: stmp
        
        
        PI=2.0_WP*(DASIN(1.0_WP))
        
        IF(MYID.NE.0) RETURN !   
       
        OPEN(INI,FILE='readdata.ini', status='old', iostat=ios)
        if(ios .NE. 0)  &
        call ERRHDL(' File readdata.ini cannot be found.',MYID)

        !==========Ignore initial introduction==================
        read(ini,*) sect
        do while(sect(1:1).EQ.';' .or. sect(1:1).EQ.'#'.or. sect(1:1).EQ.' ')
            read(ini,*) sect   ! skipping the comments lines.
        end do
       
        !==========Readin Section [geometry]======================
        lens = len_trim(sect)
        call CHKHDL(' Read in '//sect(1:lens),MYID)
        if(sect(1:lens) /= '[geometry]')  &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID)
       
        read(ini,*) skey, icase
            IF(icase==ICHANL) call CHKHDL   ('   icase=                 Plane Channel (Cartesian Cooridnates) ',MYID) 
            IF(icase==IBOX3P) call CHKHDL   ('   icase=                 BOX with 3 periodic b.c (Cartesian Cooridnates) ',MYID) 
            IF(icase==IPIPEC) call CHKHDL   ('   icase=                 Circular Tube (cylindrical Cooridnates)',MYID) 
            IF(icase==IANNUL) call CHKHDL   ('   icase=                 Annular  Tube (cylindrical Cooridnates)',MYID) 
        read(ini,*) skey, HX_tg, Hx_io
            IF(ICASE.EQ.IBOX3P) HX_io=2.0_WP*PI
            call CHKRLHDL  ('   HX_tg=          ',MYID,HX_tg)
            call CHKRLHDL  ('   HX_io=          ',MYID,HX_io)
        read(ini,*) skey, HZ
            IF(icase.EQ.IPIPEC) HZ=2.0_WP*(DASIN(1.0_WP))*2.0_WP
            IF(ICASE.EQ.IBOX3P) HZ=2.0_WP*PI
            call CHKRLHDL  ('   HZ=             ',MYID,HZ)
        read(ini,*) skey, HYB, HYT
            call CHKRLHDL  ('   Coordinate of Y (bottom)=          ',MYID,HYB)
            call CHKRLHDL  ('   Coordinate of Y (top)=             ',MYID,HYT)
            IF(ICASE==IBOX3P) THEN
                HYB = -PI
                HYT =  PI
            END IF
            
        read(ini,*) skey, NFLOW 
            if(NFLOW==1) call CHKHDL ('   NFLOW=                 X Streamwise Flow',MYID)
            if(NFLOW==2) call CHKHDL ('   NFLOW=                 Y Streamwise Flow',MYID)
            if(NFLOW==3) call CHKHDL ('   NFLOW=                 Z Streamwise Flow',MYID)
       
        !==========Readin Section [mesh]======================     
        read(ini,*) sect
        lens = len_trim(sect)
         call CHKHDL(' Read in '//sect(1:lens),MYID)  
        if(sect(1:lens) /= '[mesh]')       &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID)     
       
        read(ini,*) skey, NCL1_tg, NCL1_io
            call CHKINTHDL ('   NCL1= (TurGe)',MYID,NCL1_tg )
            call CHKINTHDL ('   NCL1= (MAIN) ',MYID,NCL1_io )
        
        read(ini,*) skey, NCL2
            call CHKINTHDL ('   NCL2=        ',MYID,NCL2)
        read(ini,*) skey, NCL3
            call CHKINTHDL ('   NCL3=        ',MYID,NCL3)
        
        read(ini,*) skey, STR2
            call CHKRLHDL  ('   STR2=        ',MYID,STR2)
        read(ini,*) skey, ISTR2
            call CHKINTHDL ('   ISTR2=       ',MYID,ISTR2)
       
        !==========Readin Section [boundary]======================     
        read(ini,*) sect
        lens = len_trim(sect)
         call CHKHDL(' Read in '//sect(1:lens),MYID)  
        if(sect(1:lens) /= '[boundary]')       &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID)     
       
        read(ini,*) skey, BCX_tg(1), BCX_tg(2)
            if( BCX_tg(1) == 3 ) BCX_tg(2) = 3
            if( BCX_tg(2) == 3 ) BCX_tg(1) = 3
        
        read(ini,*) skey, BCX_io(1), BCX_io(2) 
            if( BCX_io(1) == 3 ) BCX_io(2) = 3
            if( BCX_io(2) == 3 ) BCX_io(1) = 3
        
            IF(NCL1_io .GT. 2) THEN
                IOFLOWflg =  .TRUE. 
            ELSE
                IOFLOWflg =  .FALSE. 
            END IF
            
            IF( (BCX_io(1) == 3) .and. (NCL1_io .GT. 2)) THEN
                TGFLOWflg =  .FALSE. 
            ELSE
                TGFLOWflg =  .TRUE. 
            END IF
        
            IF (TGFLOWflg) THEN
                IF(BCX_tg(1)==1) THEN
                    call CHKHDL('   BC: TG INLET=           DIRICHELET ',MYID)
                ELSE IF(BCX_tg(1)==2) THEN
                    call CHKHDL('   BC: TG INLET=           NEUMANN ',MYID)
                ELSE IF(BCX_tg(1)==3) THEN
                    call CHKHDL('   BC: TG INLET=           PERIODIC ',MYID)
                ELSE
                    CALL ERRHDL('   No Such b.c. for BCX_tg(1)',MYID)     
                END IF
                IF(BCX_tg(2)==1) THEN
                    call CHKHDL('   BC: TG OUTLET=          DIRICHELET ',MYID)
                ELSE IF(BCX_tg(2)==2) THEN
                    call CHKHDL('   BC: TG OUTLET=          NEUMANN ',MYID)
                ELSE IF(BCX_tg(2)==3) THEN
                    call CHKHDL('   BC: TG OUTLET=          PERIODIC ',MYID)
                ELSE
                    CALL ERRHDL('   No Such b.c. for BCX_tg(2)',MYID)     
                END IF
            END IF
        
            IF (IOFLOWflg ) THEN
                IF(BCX_io(1)==1) THEN
                    call CHKHDL('   BC: MAIN DOMAIN INLET=  DIRICHELET ',MYID)
                ELSE IF(BCX_io(1)==2) THEN
                    call CHKHDL('   BC: MAIN DOMAIN INLET=  NEUMANN ',MYID)
                ELSE IF(BCX_io(1)==3) THEN
                    call CHKHDL('   BC: MAIN DOMAIN INLET=  PERIODIC ',MYID)
                ELSE
                    CALL ERRHDL('   No Such b.c. for BCX_IO(1)',MYID)     
                END IF
                
                IF(BCX_io(2)==1) THEN
                    call CHKHDL('   BC: MAIN DOMAIN OUTLET= DIRICHELET ',MYID)
                ELSE IF(BCX_io(2)==2) THEN
                    call CHKHDL('   BC: MAIN DOMAIN OUTLET= NEUMANN ',MYID)
                ELSE IF(BCX_io(2)==3) THEN
                    call CHKHDL('   BC: MAIN DOMAIN OUTLET= PERIODIC ',MYID)
                ELSE
                    CALL ERRHDL('   No Such b.c. for BCX_IO(2)',MYID)     
                END IF
            END IF
       
       
        read(ini,*) skey, BCZ(1), BCZ(2)
            if( BCZ(1) == 3 ) BCZ(2) = 3
            if( BCZ(2) == 3 ) BCZ(1) = 3
       
            IF(BCZ(1)/=3) THEN 
                call ERRHDL('Z MUST BE PERIODIC',MYID)
            ELSE
                call CHKHDL('   BC: SPANWISE=           PERIODIC ',MYID)
            END IF
        
            IF(TGFLOWflg .AND. IOFLOWflg) THEN
                call CHKHDL('   FLOW DOMAIN:            PERIODIC TG + INLET/OUTLET MAIN DOMAIN',MYID)
            END IF
            IF(TGFLOWflg .AND. (.NOT.IOFLOWflg)) THEN
                call CHKHDL('   FLOW DOMAIN:            PERIODIC FLOW (TG) ONLY',MYID)
            END IF
            IF(IOFLOWflg .AND. (.NOT.TGFLOWflg)) THEN
                call CHKHDL('   FLOW DOMAIN:            PERIODIC FLOW (IO) ONLY',MYID)
            END IF
            
        read(ini,*) skey, BCY(1), BCY(2)
            IF(BCY(1)==3) THEN 
                call CHKHDL('   BC: Y=                  PERIODIC ',MYID)
            ELSE IF(BCY(1)==2) THEN
                call CHKHDL('   BC: Y=                  SYMMETIC ',MYID)
            ELSE IF(BCY(1)==1) THEN
                call CHKHDL('   BC: Y=                  WALL ',MYID)
            ELSE
                call ERRHDL('ERROR! NO SUCH A B.C. IN Y DIRECTION!',MYID)
            END IF
       
        !==========Readin Section [numerical]===================
        read(ini,*) sect
        lens = len_trim(sect)
         call CHKHDL(' Read in '//sect(1:lens),MYID)  
        if(sect(1:lens) /= '[numerical]')       &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID)     
       
        read(ini,*) skey, DT, DTMIN
             call CHKRLHDL  ('   DT=              ',MYID,DT)
             call CHKRLHDL  ('   DTMIN=           ',MYID,DTMIN)
        read(ini,*) skey, CFLGV
             call CHKRLHDL  ('   CFLGV=           ',MYID,CFLGV)
            
        read(ini,*) skey, weightedpressure
            IF(weightedpressure == 1) call CHKHDL ('   weightedpressure=       yes',MYID)
            IF(Weightedpressure == 0) call CHKHDL ('   weightedpressure=       no',MYID)
            
        read(ini,*) skey, visthemflg
            IF(visthemflg == visimplicit) call CHKHDL ('   visthemflg=             implicit',MYID)
            IF(visthemflg == visexplicit) call CHKHDL ('   visthemflg=             explicit',MYID)

        read(ini,*) skey, VPERG
            call CHKRLHDL  ('   VPERG=           ',MYID,VPERG)
        read(ini,*) skey, SVPERG
            call CHKRLHDL  ('   SVPERG=          ',MYID,SVPERG)
        
        read(ini,*) skey, RSTflg_tg, RSTflg_io
        read(ini,*) skey, RSTtim_tg, RSTtim_io
        read(ini,*) skey, RST_type_flg, RST_time_set_flg
            IF(TGFLOWflg) THEN
                IF(RSTflg_tg==0)      call CHKHDL ('   RSTflg_tg=              TG: Start from random flow field',MYID)
                IF(RSTflg_tg==1)      call CHKHDL ('   RSTflg_tg=              TG: Start from interpolation of a coarse mesh',MYID)
                IF(RSTflg_tg==2)      call CHKHDL ('   RSTflg_tg=              TG: Restart from last step',MYID)
            END IF
            
            IF(IOFLOWflg) THEN
                IF(RSTflg_io==0)      call CHKHDL ('   RSTflg_io=              IO: Start from random flow field',MYID)
                IF(RSTflg_io==1)      call CHKHDL ('   RSTflg_io=              IO: Start from interpolation of a coarse mesh',MYID)
                IF(RSTflg_io==2)      call CHKHDL ('   RSTflg_io=              IO: Restart from last step',MYID)
            END IF
            
            IF(RST_type_flg==0)       call CHKHDL ('   RST_type_flg=           IO: Restart both, flow and thermal fields',MYID)
            IF(RST_type_flg==1)       call CHKHDL ('   RST_type_flg=           IO: Restart only flow field, no thermal field',MYID)
            IF(RST_time_set_flg==0)   call CHKHDL ('   RST_time_set_flg=       IO: Restart following previous time.',MYID)
            IF(RST_time_set_flg==1)   call CHKHDL ('   RST_time_set_flg=       IO: Restart with re-setting time to zero',MYID)
        
        !==========Readin Section [fluid]=====================
        read(ini,*) sect
        lens = len_trim(sect)
         call CHKHDL(' Read in '//sect(1:lens),MYID)  
        if(sect(1:lens) /= '[fluid]') &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID) 
       
        read(ini,*) skey, REN 
             call CHKRLHDL  ('   Re No.=          ',MYID,REN)
        read(ini,*) skey, REINI, TLGRE
             call CHKRLHDL  ('   REINI=           ',MYID,REINI)
             call CHKRLHDL  ('   TLGRE=           ',MYID,TLGRE)
        read(ini,*) skey, FLOWDRVTP
            if(FLOWDRVTP==1) call CHKHDL ('   FLOWDRVTP=              Constant mass flux driven',MYID)
            if(FLOWDRVTP==2) call CHKHDL ('   FLOWDRVTP=              Constant pressure gradient driven',MYID)
        read(ini,*) skey, CFGV
            call CHKRLHDL  ('   CFGV=           ',MYID,CFGV)
      

        !==========Readin Section [thermal]=====================
        read(ini,*) sect
        lens = len_trim(sect)
         call CHKHDL(' Read in '//sect(1:lens),MYID)  
        if(sect(1:lens) /= '[thermal]') &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID)

        read(ini,*) skey, thermlflg
            if(thermlflg==0) then
                call CHKHDL    ('   Thermlflg=              Only Velocity+P field',MYID)
            else if(thermlflg==1) then
                call CHKHDL    ('   Thermlflg=              Velocity+Thermal field',MYID)
            end if
            
        read(ini,*) skey, FLOWDIR
            IF (thermlflg==1) THEN
                if(FLOWDIR==0) call CHKHDL ('   FLOWDIR=                Horizontal flow',MYID)
                if(FLOWDIR==1) call CHKHDL ('   FLOWDIR=                Upwards flow',MYID)
                if(FLOWDIR==2) call CHKHDL ('   FLOWDIR=                Downwards flow',MYID)
            END IF
        read(ini,*) skey, GRAVDIR
            IF (thermlflg==1) THEN
                if(GRAVDIR==0) CALL CHKHDL ('   GRAVITY=                Not considered',MYID)
                if(GRAVDIR==1) CALL CHKHDL ('   GRAVITY=                X DIRECTION',MYID)
                if(GRAVDIR==2) CALL CHKHDL ('   GRAVITY=                Y DIRECTION',MYID)
                if(GRAVDIR==3) CALL CHKHDL ('   GRAVITY=                Z DIRECTION',MYID)
            END IF
        
        read(ini,*) skey, BCWALLHEAT(Ibotwall), BCWALLHEAT(Itopwall)   
        read(ini,*) skey, WHEAT0_DIM(Ibotwall), WHEAT0_DIM(Itopwall)
            IF (thermlflg==1) THEN
                IF(ICASE==IPIPEC) THEN
                    BCWALLHEAT(Ibotwall) = 0
                    WHEAT0_DIM(Ibotwall) = 0.0_WP
                END IF
                
                IF(BCWALLHEAT(Ibotwall)==isoFluxWall) THEN
                    CALL CHKHDL    ('   BCWALLHEAT (bottomwall) = Constant Wall Heat Flux (isoflux-wall)',MYID)
                    call CHKRLHDL  ('     Wall Heat flux (W/m2) on the bottom wall= ',MYID,WHEAT0_DIM(Ibotwall))
                ELSE IF (BCWALLHEAT(Ibotwall)==isoThermalWall) THEN
                    CALL CHKHDL    ('   BCWALLHEAT (bottomwall) = Constant Wall Temperature(isothermal-wall)',MYID)
                    call CHKRLHDL  ('     Wall Temperature (K)  on the bottom wall=  ',MYID,WHEAT0_DIM(Ibotwall))
                ELSE
                END IF
                
                IF(BCWALLHEAT(Itopwall)==isoFluxWall) THEN
                    CALL CHKHDL    ('   BCWALLHEAT (topwall)   = Constant Wall Heat Flux (isoflux-wall)',MYID)
                    call CHKRLHDL  ('     Wall Heat flux (W/m2) on the top wall=     ',MYID,WHEAT0_DIM(Itopwall))
                ELSE IF (BCWALLHEAT(Itopwall)==isoThermalWall) THEN
                    CALL CHKHDL    ('   BCWALLHEAT (topmwall)  = Constant Wall Temperature(isothermal-wall)',MYID)
                    call CHKRLHDL  ('     Wall Temperature (K)  on the top wall=     ',MYID,WHEAT0_DIM(Itopwall))
                ELSE
                END IF
            END IF
        
        read(ini,*) skey, L0
        read(ini,*) skey, T0
        read(ini,*) skey, P0
            IF (thermlflg==1) THEN
                call CHKRLHDL  ('   Ref: Length (m) =            ',MYID,L0)
                call CHKRLHDL  ('   Ref: Temperature (K) T0=     ',MYID,T0)
                call CHKRLHDL  ('   Ref: Pressure (Pa) P0=       ',MYID,P0)
            END IF
        
        read(ini,*) skey, thermoStat
            IF (thermlflg==1) THEN
                IF(thermoStat==search_table) CALL CHKHDL    ('   thermoStat=                    Searching Table',MYID)
                IF(thermoStat==idealgas_law) CALL CHKHDL    ('   thermoStat=                    Perfect gas power relation',MYID)
            END IF
            IF(thermoStat==search_table) THEN
                read(ini,*) skey, itmp
                read(ini,*) skey, NISTFLNM
                IF (thermlflg==1) THEN
                    call CHKHDL    ('   PROPREF FILE=                '//TRIM(NISTFLNM),MYID)
                END IF
            ELSE IF(thermoStat==idealgas_law) THEN
                read(ini,*) skey, IdealGasType, fstatetype
                read(ini,*) skey, stmp
                IF (MYID.EQ.0 .and. thermlflg==1) THEN
                    IF(IdealGasType==monatomic_gas)  CALL CHKHDL    ('   IdealGasType=                    monatomic_gas',MYID)
                    IF(IdealGasType==diatomic_gas)   CALL CHKHDL    ('   IdealGasType=                    diatomic_gas (air)',MYID)
                    IF(IdealGasType==trivalence_gas) CALL CHKHDL    ('   IdealGasType=                    trivalence_gas',MYID)
                    IF(fstatetype==sutherlandflg)    CALL CHKHDL    ('   Property Law=                    Sutherland Law',MYID)
                    IF(fstatetype==powerlawflg)      CALL CHKHDL    ('   Property Law=                    Power Law',MYID)
                END IF
            ELSE
            END IF

        !==========Readin Section [statistics]===============
        read(ini,*) sect
        lens = len_trim(sect)
         call CHKHDL(' Read in '//sect(1:lens),MYID)  
        if(sect(1:lens) /= '[statistics]')  &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID) 
       
        read(ini,*) skey, TSTOP 
             call CHKRLHDL  ('   TSTOP           =',MYID,TSTOP)
        read(ini,*) skey, TSTAV1, TSTAV_RESET
             call CHKRLHDL  ('   TSTAV1          =',MYID,TSTAV1)
             call CHKRLHDL  ('   TSTAV_RESET     =',MYID,TSTAV_RESET)
        !IF((TSTAV_RESET+1.0E-2_WP).LT.TSTAV1) CALL ERRHDL('TSTAV_RESET should not be smaller than TSTAV1',myid)
        read(ini,*) skey, DTSTATEC
             call CHKRLHDL  ('   DTSTA_TEC       =',MYID,DTSTATEC)
        read(ini,*) skey, DTTECCK
             call CHKRLHDL  ('   DTTECCK         =',MYID,DTTECCK)
        read(ini,*) skey, DTSAVE1
             call CHKRLHDL  ('   DTSAVE1         =',MYID,DTSAVE1)
        read(ini,*) skey, ITSCN
             call CHKINTHDL ('   ITSCN           =',MYID,ITSCN)
        read(ini,*) skey, MGRID, JINI
             call CHKINTHDL ('   MGRID           =',MYID,MGRID)
             call CHKINTHDL ('   JINI            =',MYID,JINI)
            
            
        !==========Readin Section [postprocess]=====================
        read(ini,*) sect
        lens = len_trim(sect)
         call CHKHDL(' Read in '//sect(1:lens),MYID)  
        if(sect(1:lens) /= '[postprocess]') &
        call ERRHDL(' Reading fails: '//sect(1:lens),MYID) 
       
        read(ini,*) skey, pprocessonly
            IF(pprocessonly==0)   call CHKHDL ('   pprocessonly=      normal restart',MYID)
            IF(pprocessonly==1)   call CHKHDL ('   pprocessonly=      postprocessing data only',MYID)
            IF(pprocessonly==2)   call CHKHDL ('   pprocessonly=      postprocessing data from given instantanous flow field',MYID)

        read(ini,*) skey, ppinst
            IF(ppinst==1)         call CHKHDL ('   ppinst=            postprocessing instantanous isosurface',MYID)
            IF(ppinst==0)         call CHKHDL ('   ppinst=            not postprocessing instantanous isosurface',MYID)
            
        read(ini,*) skey, ppspectra
            IF(ppspectra==1)      call CHKHDL ('   ppspectra=         postprocessing energy spectra and correlations',MYID)
            IF(ppspectra==0)      call CHKHDL ('   ppspectra=         not postprocessing energy spectra and correlations',MYID)

        read(ini,*) skey, ppdim
            IF(ppdim==1)          call CHKHDL ('   ppdim=             output both dim and undim results',MYID)
            IF(ppdim==0)          call CHKHDL ('   ppdim=             output only undim results',MYID)

        read(ini,*) skey, teczonename
            call CHKHDL    ('   teczonename=                '//TRIM(teczonename),MYID)
        
        read(ini,*) skey, pp_instn_sz
            IF(pprocessonly==2) THEN
                call CHKINTHDL ('   Steps of given instantanous flow field =',MYID,pp_instn_sz)
                IF(pprocessonly==2 .and. pp_instn_sz.LT.1) THEN
                    call ERRHDL(' ppinstansz is less than 1 for postprocessing instantanous flow fields',MYID) 
                END IF
            
                ALLOCATE (pp_instn_tim(pp_instn_sz)); pp_instn_tim=0.0_WP
                DO N=1, pp_instn_sz
                    read(ini,*) pp_instn_tim(N)
                    call CHKRLHDL  ('       pp_instn_tim   =',MYID,pp_instn_tim(N))
                END DO
                
            ELSE 
            
                DO N=1, pp_instn_sz
                    read(ini,*) rtmp
                END DO
    
            END IF
        
        close(ini)

        !================end of reading in data================================
        
        IF(MOD(MGRID,2)==0) MGRID =MGRID+1

        NND1_tg=NCL1_tg+1
        NND1_io=NCL1_io+1
        NND2=NCL2+1
        NND3=NCL3+1
        
        IF(IOFLOWflg) THEN
            IF(TGFLOWflg) THEN
                NCL1S = 0
                NCL1E = NCL1_io+1
            ELSE
                NCL1S = 1
                NCL1E = NCL1_io
            END IF
        END IF
        
        
        RETURN
             
    END SUBROUTINE

    SUBROUTINE BCAST_READINI
        use init_info
        use mpi_info
        use mesh_info
        use thermal_info
        IMPLICIT NONE
        
        !================[geometry]===================
        CALL MPI_BCAST( ICASE, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( HX_tg, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( HX_io, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( HZ,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( HYB,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( HYT,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( NFLOW, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
    
        !================[mesh]===================
        CALL MPI_BCAST( NCL1_tg, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NCL1_io, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NCL2,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NCL3,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( STR2,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ISTR2,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
        
        !================[boundary]===================
        CALL MPI_BCAST( BCX_tg,      2, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( BCX_io,      2, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( BCZ,         2, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( BCY,         2, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( IOFLOWflg,   1, MPI_LOGICAL,          0, ICOMM, IERROR )
        CALL MPI_BCAST( TGFLOWflg,   1, MPI_LOGICAL,          0, ICOMM, IERROR )
        
        !================[numerical]===================
        CALL MPI_BCAST( DT,           1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( CFLGV,        1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DTmin,        1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( weightedpressure, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( visthemflg,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( VPERG,        1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( SVPERG,       1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( RSTflg_tg,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( RSTflg_io,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( RSTtim_tg,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( RSTtim_io,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( RST_type_flg, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( RST_time_set_flg,1, MPI_INTEGER4,      0, ICOMM, IERROR )
        
        
        
        !================[fluid]===================
        CALL MPI_BCAST( REN,       1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( REINI,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( TLGRE,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( FLOWDRVTP, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( CFGV,      1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        !================[thermal]===================
        CALL MPI_BCAST( thermlflg,   1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( FLOWDIR,     1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( GRAVDIR,     1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( BCWALLHEAT,  2, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( WHEAT0_DIM,  2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( L0,          1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( T0,          1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( P0,          1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( thermoStat,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( IdealGasType,1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( fstatetype,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NISTFLNM,   64, MPI_CHARACTER,        0, ICOMM, IERROR )

        !================[statistics]===================
        CALL MPI_BCAST( TSTOP,      1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( TSTAV1,     1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( TSTAV_RESET,1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DTSTATEC,   1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DTTECCK,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( DTSAVE1,    1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ITSCN,      1, MPI_INTEGER4, 0, ICOMM, IERROR )
        CALL MPI_BCAST( MGRID,      1, MPI_INTEGER4, 0, ICOMM, IERROR )
        CALL MPI_BCAST( JINI,       1, MPI_INTEGER4, 0, ICOMM, IERROR )
        
        !================[postprocess]===================
        CALL MPI_BCAST( pprocessonly, 1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( ppinst,       1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( ppspectra,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( ppdim,        1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( teczonename, 64, MPI_CHARACTER,        0, ICOMM, IERROR )
        CALL MPI_BCAST( pp_instn_sz,  1, MPI_INTEGER4,         0, ICOMM, IERROR )
        IF(pprocessonly==2) THEN
            if(myid.ne.0) then
                allocate(pp_instn_tim(pp_instn_sz))
                pp_instn_tim=0.0_WP
            end if
            CALL MPI_BCAST( pp_instn_tim, pp_instn_sz, MPI_DOUBLE_PRECISION,         0, ICOMM, IERROR )
        END IF
            
        
        !================all others=================================================
        CALL MPI_BCAST( NND1_tg,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NND1_io,    1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NND2,       1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NND3,       1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NCL1S,      1, MPI_INTEGER4,         0, ICOMM, IERROR )
        CALL MPI_BCAST( NCL1E,      1, MPI_INTEGER4,         0, ICOMM, IERROR )
        
        !IF(myid==1) write(*,*) NISTFLNM !test
        
        RETURN
    END SUBROUTINE
