# Universal machine independent makefile for ROMS model
#========== ======= =========== ======== === ==== =====
# Set machine dependent definitions and rules.

include Makedefs.inc


# Build ROMS source list from unprocessed Fortran files
   SRCS_F77 = $(wildcard *$(UPF77_ext))
   SRCS_F90 = $(wildcard *$(UPF90_ext))
   SRCS = $(SRCS_F77) $(SRCS_F90)
# Exclude files that are not part of main build (just f77 for now)
   EXCL = mpi_test srcscheck
   EXCL +=  int_r3d sediment visc3d_GP visc3d_S_FS
   EXCL := $(addsuffix $(UPF77_ext), $(EXCL))
   SRCS_F77 := $(filter-out $(EXCL), $(SRCS_F77))

# Include files that do not exist yet but
# 	are generated as part of the build (just f77 for now):
   INCL = check_srcs check_switches1
   INCL := $(addsuffix $(UPF77_ext), $(INCL))
   SRCS_F77 := $(SRCS_F77) $(INCL)
   SRCS_F77 := $(SRCS_F77)

# Object and pre-processed source files are source list with extension changes
   RCS_F77 = $(SRCS_F77:$(UPF77_ext)=$(PPF77_ext))
   RCS_F90 = $(SRCS_F90:$(UPF90_ext)=$(PPF90_ext))

   OBJS_F77 = $(RCS_F77:$(PPF77_ext)=$(OBJ_ext))
   OBJS_F90 = $(RCS_F90:$(PPF90_ext)=$(OBJ_ext))
   OBJS = $(OBJS_F77) $(OBJS_F90)
# Sort OBJS to remove duplicates
  OBJS := $(sort $(OBJS))
  $(info DEBUG: OBJS = $(OBJS))

  SBIN = roms
 LROMS = libroms.a
 LROMS2 = 2/$(LROMS)
 LROMS3 = 2/$(LROMS) 3/$(LROMS)
 LROMS4 = 2/$(LROMS) 3/$(LROMS) 4/$(LROMS)

# Non-hydrostatic solver library (NHMG): always compiled and linked, and
# switched on at run time by the namelist key nhmg (group NHMG_SETTINGS).
# It is built in the subdirectory nhmg/ with this build's compiler and flags,
# so the source tree is never written to and every build has its own copy.
   NHMG_SRC   = $(abspath $(ROMS_ROOT)/NHMG/src)
   NHMG_BUILD = nhmg
   NHMG_LIB   = $(NHMG_BUILD)/libnhmg.a
   NHMG_NAMES = mg_cst mg_tictoc mg_mpi mg_namelist mg_netcdf_out mg_grids\
                mg_mpi_exchange mg_gather mg_intergrids mg_relax mg_horiz_grids\
                mg_vert_grids mg_projection mg_solvers mg_autotune nhmg_debug nhmg
   NHMG_OBJS  = $(addprefix $(NHMG_BUILD)/, $(addsuffix .o, $(NHMG_NAMES)))
# ROMS sources that use an NHMG module, which must wait for the library
   NHMG_USERS = $(addsuffix .o, $(basename $(shell grep -liE \
                '^[[:space:]]*use[[:space:]]+(nhmg|mg_[a-z_]+)[[:space:]]*(,|$$)' $(SRCS) 2>/dev/null)))
   LCDF += -I$(NHMG_BUILD)

###########
# TARGETS #
###########

# ROMS excutable (This is the first target and hence the default):
$(SBIN): $(OBJS) $(NHMG_LIB)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(NHMG_LIB) $(LCDF) $(LMPI)
	mv a.out $(SBIN)

# NHMG library. A compiler writes its .mod files to the directory it runs in,
# so each file is compiled from inside $(NHMG_BUILD). Every object waits for
# that directory (it is on the include path); the NHMG users wait for the library.
$(OBJS): | $(NHMG_BUILD)
# Changing a build option (make MARBL=yes, PIO=yes ...) recompiles everything
$(OBJS): build_options.txt
ifneq ($(strip $(NHMG_USERS)),)
$(NHMG_USERS): | $(NHMG_LIB)
endif

$(NHMG_BUILD):
	mkdir -p $@

$(NHMG_LIB): $(NHMG_OBJS)
	/bin/rm -f $@
	ar rcs $@ $(NHMG_OBJS)

$(NHMG_BUILD)/%.o: $(NHMG_SRC)/%.f90 | $(NHMG_BUILD)
	cd $(NHMG_BUILD) && $(CFT) -c $(FFLAGS) -o $(@F) $<

# Module dependencies within NHMG (NHMG/src/depend.mk)
NHMG_DEP = $(addprefix $(NHMG_BUILD)/, $(addsuffix .o, $(1)))
$(call NHMG_DEP, mg_mpi mg_namelist mg_netcdf_out): $(call NHMG_DEP, mg_cst)
$(call NHMG_DEP, mg_namelist mg_grids): $(call NHMG_DEP, mg_tictoc)
$(call NHMG_DEP, mg_grids nhmg_debug): $(call NHMG_DEP, mg_mpi)
$(call NHMG_DEP, mg_grids): $(call NHMG_DEP, mg_namelist)
$(call NHMG_DEP, mg_mpi_exchange mg_gather): $(call NHMG_DEP, mg_grids)
$(call NHMG_DEP, mg_intergrids mg_horiz_grids mg_vert_grids mg_projection): $(call NHMG_DEP, mg_gather mg_mpi_exchange)
$(call NHMG_DEP, mg_relax): $(call NHMG_DEP, mg_mpi_exchange)
$(call NHMG_DEP, mg_vert_grids mg_projection mg_solvers nhmg_debug nhmg): $(call NHMG_DEP, mg_netcdf_out)
$(call NHMG_DEP, mg_solvers): $(call NHMG_DEP, mg_relax mg_intergrids)
$(call NHMG_DEP, mg_autotune): $(call NHMG_DEP, mg_solvers mg_projection)
$(call NHMG_DEP, nhmg): $(call NHMG_DEP, mg_autotune mg_horiz_grids mg_vert_grids)


# Multilevel libraries

$(LROMS): $(OBJS)
	/bin/rm -f $(LROMS)
	ar r $(LROMS) $(OBJS)


check_forces: check_forces$(OBJ_ext) $(SBIN)
	cp -pv check_forces$(OBJ_ext) main$(OBJ_ext)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(NHMG_LIB) $(LCDF) $(LMPI)
	mv a.out check_forces

# Everything
.PHONY: all
all: tools depend $(SBIN)

# A program to test MPI halo exchange routines.
# - ------- -- ---- --- ---- -------- ---------
MPI_TEST = mpi_test	mpi_setup	exchange	mpi_exchange4\
							mpi_exchange8WA
MPI_TEST_RCS := $(addsuffix $(PPF77_ext), $(MPI_TEST))
MPI_TEST_OBJ := $(addsuffix $(OBJ_ext)  , $(MPI_TEST))

mpi_test: $(MPI_TEST_OBJ)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpi_test $(MPI_TEST_OBJ) $(LCDF) $(LMPI)

# Auxiliary utility programs
# --------- ------- --------
 TOOLS = srcscheck redefs

tools: $(TOOLS)

    TMP = .

srcscheck: srcscheck$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o srcscheck srcscheck$(OBJ_ext)
srcscheck$(OBJ_ext): strings$(OBJ_ext)

redefs: redefs$(OBJ_ext)
	$(CFT) $(FFLAGS) $(LDFLAGS) -o redefs redefs$(OBJ_ext)



checkdefs: check_switches1$(UPF77_ext)

check_switches1$(UPF77_ext): cppdefs.opt
	@${ROMS_ROOT}/Tools-Roms/cppcheck.py
check_srcs$(UPF77_ext): srcscheck Makefile
	./srcscheck

.PHONY: depend
depend: check_switches1$(UPF77_ext)
        # makedepf90 is preferable but does not work on e.g. osx-arm64. Offer python fallback.
	@${ROMS_ROOT}/Tools-Roms/makedep.py  Make.depend  || \
	(echo "failed to generate dependency list with makedep.py, attempting to use makedepf90"; \
	 ${ROMS_ROOT}/Tools-Roms/makedepf90 $(SRCS) > Make.depend)
	@echo 'Updated Make.depend (dependency list)'

# Target to create tar file.
# ------ -- ------ --- -----
tarfile: clean
	tar cvf roms.tar Make* *$(UPF77_ext) *.h etalon_data.* README.*  *.in* *.mod

# Cleaning targets
# -------- -------
.PHONY: clean
clean:
	/bin/rm -rf *$(PPF77_ext) *$(OBJ_ext) *.a *.s *.mod *.trace *~ $(COMP_FILES) $(NHMG_BUILD)

.PHONY: allclean
allclean: clean
	cd 2; make -f ./Makefile clean; cd ..
	cd 3; make -f ./Makefile clean; cd ..
	cd 4; make -f ./Makefile clean; cd ..

.PHONY: clobber
clobber: clean
	/bin/rm -f check_switches1$(UPF77_ext) check_srcs$(UPF77_ext) build_options.txt
	/bin/rm -f $(SBIN) $(TOOLS) nsub
	/bin/rm -f core core.*
	/bin/rm -f *_genmod.mod *_genmod$(PPF90_ext)

.PHONY: help
help:
	@echo "Universal make procedure for ucla-roms, compatible with most UNIX systems."
	@echo "The following user options are available"
	@echo "COMPILER (fortran compiler):"
	@echo "  gnu :    default, use gfortran as compiler"
	@echo "  intel :  use ifx as compiler"
	@echo "  ifort : use ifort (deprecated) as compiler"
	@echo "MPI_WRAPPER (mpi wrapper command)"
	@echo "   auto : default, automatically detects wrapper by querying"
	@echo "   mpifort  : use generic mpifort wrapper around Fortran compiler"
	@echo "   mpiifx   : use Intel MPI's mpiifx wrapper around Fortran compiler"
	@echo "   mpiifort : use Intel MPI's (deprecated) mpiifort wrapper around Fortran compiler"
	@echo "BUILD_MODE:"
	@echo "   regular : default, uses optimisation level 3"
	@echo "   debug   : compile with debug flags"
	@echo "   strict  : compile with -fp-model strict (intel only)"
	@echo "   vtune   : compile with flags for profiling with vtune (intel only)"
	@echo "   grof    : compile with flags for profiling with grof (intel only)"
	@echo "MARBL (build option):"
	@echo "    no     : default, build without the MARBL biogeochemistry library"
	@echo "    yes    : compile and link MARBL (needs MARBL_ROOT); use_marbl switches it on"
	@echo "PIO (build option):"
	@echo "    no     : default, files with one piece per rank (namelist: parallel_io = .false.)"
	@echo "    yes    : compile and link ParallelIO (needs PIO_ROOT, PNETCDFHOME); joined files"
	@echo "KEEP_PPSRC:"
	@echo "    true   : keep pre-processed source code after processing"
	@echo "    false  : delete pre-processed source code after processing"

# Automatically generated dependency list:
include Make.depend

