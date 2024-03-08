#################################################################################################
# Makefile for DNSpipe, by Wei Wang, Dec 2013.                                                                                              #
# Usage:                                                                                        #
#       make all        to make all files with -O2                                              #
#       make cfg=dbg    for debugging                                                           #
#       make cfg=prof   for profiling                                                           #
#                                                                                               #
# For debugging run:                                                                            #
# mpirun -np 4 valgrind --leak-check=full --track-origins=yes \                                 #
#                       --log-file=valgrind_output.txt ./dgdes* < solver_input > solver_output  #
#                                                                                               #
#################################################################################################
FMAKE=make
SHELL  = /bin/sh

VERSION = $(shell date +%d.%m.%Y)
##OUTPUT = DNS-$(VERSION)
OUTPUT= HEFT_DNS

SRCDIR= ../src
PPCDIR= ../src/processing
FFT99DIR= ../libs/fft99
FISHPACK= ../libs/fishpack4.1
EIGENV33= ../libs/eigen33

FC = mpif90 #-vt       

#-------------------------------------------------------------------------------------#
# CAUTION: After this line, NOTHING should be really changed!!!                       #
#-------------------------------------------------------------------------------------#

FOPT = -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch #-intel
#OPT = -r8 -fallow-invalid-boz -fallow-argument-mismatch ###!check gfortran
ifeq ($(cfg), gnu)
#-O0 -DDEBUG -fbacktrace -fbounds-check -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow -finit-real=nan -fsignaling-nans -Warray-temporaries -Wall -Waliasing -Wampersand -Warray-bounds -Warray-temporaries -Wcharacter-truncation -Wconversion -Wconversion-extra -Wextra -Wline-truncation -Wintrinsics-std -Wintrinsic-shadow -Wno-align-commons -Wreal-q-constant -Wunused-parameter -Wsurprising -Wunderflow -pedantic -pedantic-errors
    FOPT12 = \
                	 -Og\
                	 -DDEBUG\
                	 -fbacktrace\
                	 -fbounds-check\
                	 -fcheck=all\
                	 -fdump-core\
                	 -ffpe-trap=invalid,zero,overflow\
                	 -finit-real=nan\
                	 -fsignaling-nans\
                	 -Warray-temporaries\
                	 -Wall\
                	 -Waliasing\
                	 -Wampersand\
                	 -Warray-bounds\
                	 -Warray-temporaries\
                	 -Wcharacter-truncation\
                	 -Wconversion\
                	 -Wconversion-extra\
                	 -Wextra\
                	 -Wline-truncation\
                	 -Wintrinsics-std\
                	 -Wintrinsic-shadow\
                	 -Wno-align-commons\
                	 -Wreal-q-constant\
                	 -Wunused-parameter\
                	 -Wsurprising\
                	 -Wunderflow\
                	 -pedantic\
                	 -pedantic-errors\
                	 -fallow-argument-mismatch
             #-mcmodel medium #
    FOPT3 = -Og -g -fno-range-check -fbacktrace -fbounds-check
    FOPT4 = -Og -g -fno-range-check -fbacktrace -fbounds-check 
    
else ifeq ($(cfg), intel)
    FOPT12 = -Og -g\
                	 -traceback\
                	 -check all\
                	 -check bounds\
                	 -check uninit\
                	 -ftrapuv\
                	 -debug all\
                	 -fpe0\
                	 -fpe3\
                	 -fpe-all=3\
                	 -assume ieee_fpe_flags\
                	 -ftz\
			 -fallow-argument-mismatch
              
    FOPT3 = -Og -g -traceback -check all -check bounds 
    FOPT4 = -Og -g -traceback -check all -check bounds  

else
    FOPT12 = -O2 -fallow-argument-mismatch -fbacktrace
    FOPT3 =  -O2 
    FOPT4 =  -O2   

endif

OBJ1= \
        modules.o\
        MEMALLOCT.o\
        DNS_THERMAL.o\
        BC_COUTLET_ENEG_RK3.o\
        BC_COUTLET_MOM_RK3.o\
        BC_WALL.o\
        BC_INLET.o\
        CFL.o\
        CONSPARA.o\
        CONS_RKCOEF.o\
        CONVECTION_X.o\
        CONVECTION_Y.o\
        CONVECTION_Z.o\
        COORDJR.o\
        DIVGCK.o\
        DIVG.o\
        Flow_State.o\
        FLOWSTART.o\
        functions.o\
        INDEXL2G.o\
        INDEXSET.o\
        INITIAL_FL_FLD_TG.o\
        INITIAL_FL_THEML_FLD_io.o\
        INTFC_VARS.o\
        INTFC_THERMAL_io.o\
        LAMPOISLPROF.o\
        TDMA_COEF.o\
        MASSFLUX_CALC_IO.o\
        mesh_decomp.o\
        possion3d_FFT99.o\
        possion3d_fishpack.o\
        MOMFA_tg.o\
        MOMFA_io.o\
        PRESSURE_CORRECTION.o\
        READINI.o\
        RHS_CvLpGpS_tg.o\
        RHS_ENERGY_EXPLICIT.o\
        RHS_MOM_EXPLICIT_io.o\
        SOLVE.o\
        SOLVERRK3_ENG_IO.o\
        SOLVERRK3_MOM_IO.o\
        SOLVERRK3_MOM_tg.o\
        start_mpi.o\
        TDMAIJI_CYC.o\
	TDMAIJJ_CYC.o\
        TDMAIJI_nonCYC.o\
        TDMAIJJ_nonCYC.o\
        THERM_PROP_UPDATE.o\
        TRASP23.o\
        VELO_CORRECTION.o\
        VISCOUS_X_EXPLICIT_io.o\
        VISCOUS_Y_EXPLICIT_io.o\
        VISCOUS_Z_EXPLICIT_io.o\
        VMAV.o\
        WRTHDL.o\

OBJ2= \
        CHK_conservation.o\
        INITIAL_MEAN.o\
        INTERPOLATION_INITIAL.o\
        READ_INTP_3D_VARS.o\
        POSTPROESS.o\
        POSTPROESS_integral_instants.o\
        PP_SCN_MONITOR.o\
        PP_SPECOSPEC.o\
        PP_MEAN_T.o\
        PP_MEAN_Z.o\
        PP_MEAN_ZX.o\
        PP_TEC360_DATA_CHECK.o\
        PP_QuadrantAnalysis.o\
        WRT_TEC_AVERAGE_XZperiodic_IO.o\
        WRT_TEC_AVERAGE_Zperiodic_IO.o\
        WRT_TEC_AVERAGE_XZperiodic_TG.o\
        WRT_VARS_TG.o\
        WRT_VARS_IO.o\
        RESTART_io.o\
        RESTART_tg.o\
        randomgen.o\
             
OBJ3= \
        CFFT99.o\
        CFTFAX.o\
        CFTRIG.o\
        FACT.o\
        FAX.o\
        FFT99A.o\
        FFT99B.o\
        FFT99.o\
        FFTFAX.o\
        FFTRIG.o\
        VPASSM.o\
        
OBJ4= \
        fftpack.o\
        
OBJ5= \
        dsyev2.o\
        dsyevc3.o\
        dsyevd3.o\
        dsyevh3.o\
        dsyevj3.o\
        dsyevq3.o\
        dsyevv3.o\
        dsytrd3.o\
        slvsec3.o\
        sqrabs.o\
        zheevc3.o\
        zheevd3.o\
        zheevh3.o\
        zheevj3.o\
        zheevq3.o\
        zheevv3.o\
        zhetrd3.o\
        cbabk2.o\
		cbal.o\
		cdiv.o\
		cg.o\
		comqr2.o\
		comqr.o\
		corth.o\
		csroot.o\
		pythag.o\
   
default:
	@cd bin; echo $(shell pwd); make $(OUTPUT) -f ../Makefile

$(OUTPUT): $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5)
	@echo -n "Linking... "
	$(FC) -o $(OUTPUT) *.o $(FOPT12) $(FOPT3) $(FOPT4) $(FOPT)
	@echo -e "Done.\a"
                	
$(OBJ1): %.o: $(SRCDIR)/%.f90
	$(FC) -c   $(FOPT12) $(FOPT)$<

$(OBJ2): %.o: $(PPCDIR)/%.f90
	$(FC) -c  $(FOPT12) $(FOPT)$<
	
$(OBJ3): %.o: $(FFT99DIR)/%.f90
	$(FC) -c  $(FOPT3) $(FOPT12) $(FOPT)$<

$(OBJ4): %.o: $(FISHPACK)/%.f
	$(FC) -c  $(FOPT4) $(FOPT)$<
	
$(OBJ5): %.o: $(EIGENV33)/%.f
	$(FC) -c  $(FOPT4) $(FOPT)$<

clean:
	rm -f bin/*

#clear
#@echo "All binary files and executeable program deleted."

all:
	@make clean -f Makefile

	@cd bin; make $(OUTPUT) -f ../Makefile
