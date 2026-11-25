/* This is "global_definitions.h":  A set of predetermined CPP macros
 and switches inserted into other files of the code during compilation
 process. General user is strongly discouraged from attempts to modify
 anything below this line.
 ------------------------------------------------------------------ */

#include "retired_cpp_keys.h"

/* MPI parallelization is always on: the code cannot be built without it. */

#define MPI

/* The model is always three-dimensional, saves two time levels in its
 restart files (exact restart), couples the barotropic mode at the
 corrector stage with extrapolated barotropic fluxes, solves the bottom
 drag implicitly, uses a nonuniform density in the barotropic pressure
 gradient, and writes one file per MPI rank.  These used to be CPP
 switches (SOLVE3D, EXACT_RESTART, CORR_COUPLED_MODE, EXTRAP_BAR_FLUXES,
 IMPLCT_NO_SLIP_BTTM_BC, VAR_RHO_2D, PARALLEL_FILES, ...); the source no
 longer tests them.  SOLVE3D and EXACT_RESTART stay defined so that
 analytical set-up files (ana_*.h) copied from older versions, which
 test them, still compile the same code. */

#define SOLVE3D
#define EXACT_RESTART

/* Switch to recompute Coriolis terms in barotropic momentum equations
 at every barotropic time step (corrector-coupled mode with forward
 extrapolation of barotropic fluxes). Normally kept undefined. */

!#define KEEP_CORIOLIS

/* Normally, whenever nonlinear EOS is selected, use stiffened and
 split form,
                rho(S,T,z) = rho1(S,T) + qp1(S,T)*[z+qp2*z^2]

 where all pressure dependency is explicit via "z", while "rho1","qp1"
 do not depend on pressure (i.e., "z");  "qp2" is just a constant.
 The above form facilitates computation of adiabatic differences in
 pressure-gradient scheme, prsgrd32AC1.F, which is critical for its
 accuracy and robustness.  Overall this follows Shchepetkin and
 McWilliams (2003, 2011), which in its turn is an adaptation of
 Dukowicz (2001) and earlier work of Sun, Bleck, Rooth, Dukowicz,
 Chassignet, and Killworth (1999) with respect to this matter.   */

#ifdef NONLIN_EOS
# define SPLIT_EOS
#endif

/* Switch ON/OFF double precision for real type variables (since this
 is mostly controlled by mpc and/or compiler options, this CPP-switch
 affects only on the correct choice of netCDF functions, see below)
 and the use QUAD precision for global summation variables, which is
 always desirable, but some compilers do not support it.  */

#define DBLEPREC

/* Define standard dimensions for the model arrays (vertical dimensions
 are inserted explicitly in the code, when needed).  Array dimensions
 correspond to a portion of the physical domain, with two ghost zones on
 each side, whatever the periodicity. The data for these two ghost zones
 is exchanged by message passing.  */

#ifdef MPI
# define GLOBAL_2D_ARRAY -1:Lm+2,-1:Mm+2
# define START_2D_ARRAY -1,-1
# define GLOBAL_1DX_ARRAY -1:Lm+2
# define GLOBAL_1DY_ARRAY -1:Mm+2
!# define START_1D -1
!# define END_1DX Lm+2
!# define END_1DY Mm+2
#else
# define GLOBAL_2D_ARRAY -1:Lm+2,-1:Mm+2
# define START_2D_ARRAY -1,-1
#endif

#define PRIVATE_1D_SCRATCH_ARRAY istr-2:iend+2
#define PRIVATE_2D_SCRATCH_ARRAY istr-2:iend+2,jstr-2:jend+2

/* The following macros contain logical expressions which answer
 the question: ''Am I a thread working on subdomain (tile) which is
 adjacent to WESTERN[EASTERN, SOUTHERN, NORTHERN] edge (i.e. physical
 boundary) of the model domain?''  Note that ghost points associated
 with periodicity are NOT considered as physical boundary points by
 these macros: in a periodic direction (ew_periodic, ns_periodic in
 open_boundary_mod.F90) they are .false. */

#ifdef MPI
# define WESTERN_EDGE (inode==0.and..not.ew_periodic)
# define EASTERN_EDGE (inode==np_xi-1.and..not.ew_periodic)
# define SOUTHERN_EDGE (jnode==0.and..not.ns_periodic)
# define NORTHERN_EDGE (jnode==np_eta-1.and..not.ns_periodic)
#else
# define WESTERN_EDGE (.not.ew_periodic)
# define EASTERN_EDGE (.not.ew_periodic)
# define SOUTHERN_EDGE (.not.ns_periodic)
# define NORTHERN_EDGE (.not.ns_periodic)
#endif

/* The following four macros identify position of an MPI-node
 relatively to the edge of the physical grid. They are similar to
 the above, except that: (i) that they apply only to MPI subdomain
 decomposition, hence do not refer to tiling bonds istr,...,jend;
 and, in addition to  that, (ii) their state does not depend on
 periodicity.  */

#ifdef MPI
# define WESTERN_MPI_EDGE iwest+iSW_corn==1
# define EASTERN_MPI_EDGE ieast+iSW_corn==LLm
# define SOUTHERN_MPI_EDGE jsouth+jSW_corn==1
# define NORTHERN_MPI_EDGE jnorth+jSW_corn==MMm
#endif


/* Sometimes an operation needs to be restricted to single MPI process
 (master process). Typically this happens when it is desirable to avoid
 redundant write of the same message by all MPI processes into stdout.
 There are two levels of verbosity avoidance: "mpi_master_only" and
 "mpi_nonexit_warn"; both restrict output to the master process, so
 ranks > 0 print only error messages.  MYID, appended (without a
 preceding comma) to the end of a message, would print the MPI-node
 number; it is disabled.   */

#ifdef MPI
# define mpi_master_only if (mynode==0)
# define mpi_nonexit_warn if (mynode==0)
# define MYID !
#else
# define mpi_master_only
# define mpi_nonexit_warn
# define MYID !
#endif

/* Similarly, if operation needed to be done by one thread only, e.g.,
 copy a redundantly computed private scalar into shared scalar, or to
 write an error message in situation where it is guaranteed that the
 error condition is discovered redundantly by every thread (and is
 the same for all) and only one needs to complain.  ZEROTH_TILE is
 intended to restrict the operation only to thread which is working
 on south-western tile.

 Occasionally a subroutine designed to process a tile may be called
 to process the whole domain. If it is necessary to distinguish
 whether it is being called for the whole domain (SINGLE_TILE_MODE)
 or a tile.

 All these switches are the same for MPI/nonMPI code.  */

#ifdef MPI
# define ZEROTH_TILE (istr==iwest .and. jstr==jsouth)
# define SINGLE_TILE_MODE (iend-istr==ieast-iwest .and. \
 jend-jstr==jnorth-jsouth)
#else
# define ZEROTH_TILE (istr==1 .and. jstr==1)
# define SINGLE_TILE_MODE (iend-istr==Lm-1 .and.+jend-jstr==Mm-1)
#endif


/* The first time step of a run started from a single time record
 uses a two-time-level scheme, "forw_start" (RK2: forward Euler +
 trapezoidal correction).  Restart files keep two consecutive time
 steps, so a restarted run continues exactly as an uninterrupted one. */

#define FIRST_TIME_STEP iic==forw_start
#define FIRST_2D_STEP iif==1

/* Macros to signal that the code is in predictor or corrector stage
 of 3D time stepping, and to signal the stage at which coupling with
 barotropic mode takes place.   The latter typically appears in the
 context of computing 3D --> 2D influences and is done immediately
 before the barotropic stepping (the corrector stage).  */

#define PRED_STAGE nnew==3
#define CORR_STAGE nrhs==3

#define COUPLED_STAGE_ONLY nrhs==3


/* Real kind of the global-sum variables: double precision (kind 8). */

#define QUAD 8
#define QuadZero 0.D0

/* Turn ON/OFF double precision for real type variables, associated
 intrinsic functions and netCDF library functions. It should be noted
 that because ROMS relies on compiler options and "mpc" program (see
 mpc.F) to generate double precision executable from default precision
 source code, this switch actually does NOT affect the size of real
 data and precision of the computation. Its main effect is to select
 the correct netCDF function (nf_xxxx_double/nf_xxxx_float, see below)
 to work properly, so it must be set consistently with mpc settings
 and compiler flags (if any) according to the intended accuracy.
  */

#ifdef DBLEPREC
# define float dble
!#define sqrt dsqrt
!#define SQRT dsqrt
!#define exp dexp
!#define EXP dexp
!#define dtanh dtanh
!#define TANH dtanh
#endif

/* Select the double/single precision types for netCDF input (NF_FTYPE)
 and output (NF_FOUT) functions.  History and averages are written in
 single precision; the restart file keeps the precision of the code. */

#ifdef DBLEPREC
# define NF_FTYPE nf90_double
# define nf_get_var_FTYPE nf_get_var_double
# define nf_get_var1_FTYPE nf_get_var1_double
# define nf_put_var1_FTYPE nf_put_var1_double
# define nf_get_vara_FTYPE nf_get_vara_double
# define nf_put_vara_FTYPE nf_put_vara_double
# define nf_get_att_FTYPE nf_get_att_double
# define nf_put_att_FTYPE nf_put_att_double
# define NF_FOUT nf90_float
# define nf_put_att_FOUT nf_put_att_real
#else
# define NF_FTYPE nf90_float
# define nf_get_var_FTYPE nf_get_var_real
# define nf_get_var1_FTYPE nf_get_var1_real
# define nf_put_var1_FTYPE nf_put_var1_real
# define nf_get_vara_FTYPE nf_get_vara_real
# define nf_put_vara_FTYPE nf_put_vara_real
# define nf_get_att_FTYPE nf_get_att_real
# define nf_put_att_FTYPE nf_put_att_real
# define NF_FOUT nf90_float
#endif
