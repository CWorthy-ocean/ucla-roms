include ../src/Makedefs.inc

# File extensions are set as variables in Makedefs.inc to
# account for case-insensitive filesystems. Default values are:

#UPF77_ext = .F   #unprocessed  Fortran 77 code
#UPF90_ext = .F90 #unprocessed  Fortran 90 code
#PPF77_ext = .f   #preprocessed Fortran 77 code
#PPF90_ext = .f90 #preprocessed Fortran 90 code
#OBJ_ext = .o     #object files


# - Can't use -ip optimization for tools (partit won't compile)
ifeq ($(findstring -ip,$(FFLAGS)),-ip)
    FFLAGS :=$(subst -ip,,$(FFLAGS))
endif
# --------------------------------------------------------------------

# Everything
ALL: depend partit partmap optimal_part ncjoin

all:	$(info MAKING ROMS TOOLS: ) mpc $(ALL)

#install: $(ALL)
#	mv $(ALL) ~/bin

# partit
PARTIT = partit roms_part nc_read_write
PARTIT := $(addsuffix $(OBJ_ext), $(PARTIT))
partit: $(PARTIT)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit $(PARTIT) $(LCDF)

# ncjoin
NCJOIN = ncjoin roms_part nc_read_write
NCJOIN := $(addsuffix $(OBJ_ext), $(NCJOIN))
ncjoin: $(NCJOIN)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin $(NCJOIN) $(LCDF)

# partmap
PARTMAP = partmap roms_part nc_read_write
PARTMAP := $(addsuffix $(OBJ_ext), $(PARTMAP))
partmap: $(PARTMAP)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partmap $(PARTMAP) $(LCDF)

# optimal_part
OPTIMAL_PART = optimal_part roms_part nc_read_write
OPTIMAL_PART := $(addsuffix $(OBJ_ext), $(OPTIMAL_PART))
optimal_part: $(OPTIMAL_PART)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o optimal_part $(OPTIMAL_PART) $(LCDF)


PSRC = partit.F roms_part.F nc_read_write.F
MSRC = partmap.F roms_part.F nc_read_write.F
OSRC = optimal_part.F roms_part.F nc_read_write.F
JSRC = ncjoin.F roms_part.F nc_read_write.F
SRCS = $(PSRC) $(MSRC) $(OSRC) $(JSRC)

  PCS = $(PSRC:.F=.f)
  PBJS = $(PSRC:.F=.o)
  MCS = $(MSRC:.F=.f)
  MBJS = $(MSRC:.F=.o)
  OCS = $(OSRC:.F=.f)
  OBJS = $(OSRC:.F=.o)
  JCS = $(JSRC:.F=.f)
  JBJS = $(JSRC:.F=.o)

.PHONY: depend
depend:  # makedepf90 is preferable but does not work on e.g. osx-arm64. Offer python fallback.
	@${ROMS_ROOT}/Tools-Roms/makedepf90 $(SRCS) > Make.depend  || \
	(echo "failed to generate dependency list with makedepf90, using python fallback"; \
	 ${ROMS_ROOT}/Tools-Roms/makedep.py  Make.depend)
	@echo 'Updated Make.depend (dependency list)'

# clean
clean:
	/bin/rm -f *$(PPF77_ext) *$(OBJ_ext) *.a *.s *.trace *.mod
# clobber
clobber: clean
	/bin/rm -f $(ALL) *.bak

# Note: the following .f: dependencies are needed to make everything
# go through CPP and take custom defined CPPFLAGS rather than compiling
# executables from .F directly.

chkindent: chkindent$(PPF77_ext)

# The following two are built from the same source
# and differ only CPP macro setting:

mreplace$(PPF77_ext): mreplace$(UPF77_ext)
	$(CPP) -P $(CPPFLAGS) -DREPLACE mreplace$(UPF77_ext) > mreplace$(PPF77_ext)

include Make.depend
