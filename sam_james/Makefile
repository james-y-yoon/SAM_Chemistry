# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

SAM = SAM_$(ADV_DIR)_$(SGS_DIR)_$(RAD_DIR)_$(MICRO_DIR)

# Determine platform 
PLATFORM := $(shell uname -s)

#----------------------------------
# AIX (tested only on IBM SP)
#

ifeq ($(PLATFORM),AIX)

#INC_MPI      := /usr/local/include
#LIB_MPI      := /usr/local/lib
INC_NETCDF   := /usr/local/include
LIB_NETCDF   := /usr/local/lib


FF77 = mpxlf90_r -c -qsuffix=f=f -qfixed=132
FF90 = mpxlf90_r -c -qsuffix=f=f90
CC = cc -c -DAIX
FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q -I${INC_NETCDF}
#FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en -I${INC_NETCDF}
LD = mpxlf90_r
LDFLAGS = -bmaxdata:512000000 -bmaxstack:256000000 -L${LIB_NETCDF} -lnetcdf

endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(PLATFORM),IRIX64)

INC_MPI      := /usr/local/include
LIB_MPI      := /usr/local/lib
INC_NETCDF   := /usr/local/include
LIB_NETCDF   := /usr/local/lib

FF77 = f90 -c -fixedform  -extend_source
FF90 = f90 -c -freeform
CC = cc -c -DIRIX64
FFLAGS = -O3 
#FFLAGS = -g -DEBUG:subscript_check=ON:trap_uninitialized=ON 
FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
LD = f90 
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

endif
#----------------------------------------------------------------------
# Linux, Yellowstone, Cheyenne
# #

# ifeq ($(PLATFORM),Linux)
#
# FF77 = mpif90 -c -fixed -extend_source -r8
# FF90 = mpif90 -c -free -r8
# CC = mpicc -c -O3 -DLINUX
#
# FFLAGS = -O2 -fp-model source
# FFLAGS_NOOPT = -O0 -g -ftrapuv -check all -traceback -debug -gen-interfaces -warn interfaces -fp-model source
##FFLAGS = ${FFLAGS_NOOPT}
#
# LD = mpif90  -mcmodel=large
# LDFLAGS = -lnetcdff 
#
#----------------------------------------------------------------------
# Linux, Intel Compiler
#

# ifeq ($(PLATFORM),Linux)
#
# LIB_MPI = /usr/local/pkg/iopenmpi/lib
# INC_MPI = /usr/local/pkg/iopenmpi/include
# INC_NETCDF = /nfs/user08/marat/local/include
# LIB_NETCDF = /nfs/user08/marat/local/lib
#
#
# FF77 = /usr/local/pkg/iopenmpi/bin/mpif90 -c -fixed -extend_source
# FF90 = /usr/local/pkg/iopenmpi/bin/mpif90 -c
# CC = mpicc -c -DLINUX
#
#
# FFLAGS = -O3 
#FFLAGS = -g -ftrapuv -check all
#
# FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
# LD = /usr/local/pkg/iopenmpi/bin/mpif90
# LDFLAGS = -L${LIB_NETCDF} -lnetcdf
#
#
# endif
#----------------------------------------------------------------------
# Linux, XLF compiler, Bluegene at San Diego SC
#
#ifeq ($(PLATFORM),Linux)

#INC_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/include
#LIB_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/lib

#FF77 = mpxlf90  -qarch=440 -qsuffix=f=f -qfixed=132
#FF90 = mpxlf90  -qarch=440 -qsuffix=f=f90
#CC = mpcc -c -DLinux
#FFLAGS = -c -O3 -qtune=440 -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF}
#LD = mpxlf90
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf

#endif
#----------------------------------------------------------------------
# Linux, XLF compiler, Bluegene at BNL
#

#ifeq ($(PLATFORM),Linux)
#
#ifeq ($(HOST),fen)
## Compute nodes:
#INC_NETCDF   := /gpfs/home2/marat/local/netcdf/include
#LIB_NETCDF   := /gpfs/home2/marat/local/netcdf/lib
#FF77 = mpixlf90 -qarch=440 -qsuffix=f=f -qfixed=132
#FF90 = mpixlf90 -qarch=440 -qsuffix=f=f90
#CC = mpicc -c -DLinux
#FFLAGS = -c -O3 -g -qtune=440 -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF}
#LD = mpixlf90  
#LDFLAGS = -L${LIB_NETCDF} -L/bgl/local/lib/ -L/opt/ibmmath/lib -lnetcdf -llapack.rts -lesslbg
#
#else
## front end nodes and vis nodes
#INC_NETCDF   := /gpfs/home2/marat/netcdf-vis-3.6.2/include
#LIB_NETCDF   := /gpfs/home2/marat/netcdf-vis-3.6.2/lib
#LIB_MPI = /usr/local/mpich2-1.0.5p4/lib
#INC_MPI = /usr/local/mpich2-1.0.5p4/include
#FF77 = xlf90 -qsuffix=f=f -qfixed=132
#FF90 = xlf90 -qsuffix=f=f90
#CC = cc -c -DLinux
#FFLAGS = -c -O3  -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF} -I${INC_MPI}
#LD = xlf90 
#LDFLAGS = -L${LIB_NETCDF} -L${LIB_MPI} -L/gpfs/home1/slatest/blas/BLAS  -L/gpfs/home1/slatest/vislapack-2/lapack-3.2 -lmpich -lnetcdf -llapack_RS6K -lblas_LINUX 
#endif
#
#
#endif
#----------------------------------------------------------------------
# Linux, Intel Compiler (Modified for setup at Univ. Washington)

ifeq ($(PLATFORM),Linux)

FF77 = mpif90 -c -fixed -extend_source -r8
FF90 = mpif90 -c -r8
CC = mpicc -c -DLINUX 

# Trim HOSTNAME to strip off .atmos.washington.edu
HOSTTRIM = $(basename $(basename $(basename $(basename $(HOSTNAME)))))

# Setup to choose from multiple machines
#ifneq (,$(filter $(HOSTTRIM),olympus challenger))
# Set up flags depending on machine

ifeq ($(HOSTTRIM),olympus)
  # Compiling on olympus or challenger
  # module load intel/19.0.2 netcdf/4.6.1.i19 openmpi/3.1.3
  FFLAGS = -O2 -axSSE4.2 -fp-model source -debug -traceback 
  FFLAGS_O1 = -O1 -fp-model source 
  NCPATH = ${NETCDF}
else ifeq ($(HOSTTRIM),hermes)
	FFLAGS = -O2 -axSSE4.2 -fp-model source -debug -traceback -lnetcdf
	FFLAGS_O1 = -O1 -fp-model source
	NCPATH = ${NETCDF}
else
  # Compiling on challenger or one of the olympus nodes.  Use older compiler setup.
  FFLAGS = -O2 -xHOST -fp-model source -mcmodel=large 
  FFLAGS_O1 = -O1 -fp-model source #	
  NCPATH = ${NETCDF}
endif

FFLAGS_NOOPT = -O0 -g -ftrapuv -fpe0 -check all -traceback -debug -gen-interfaces -warn interfaces -fp-model source
#FFLAGS = ${FFLAGS_NOOPT}
#FFLAGS_O1 = ${FFLAGS_NOOPT}

LD = mpif90

FFLAGS += -I$(NCPATH)/include
FFLAGS_O1 += -I$(NCPATH)/include
FFLAGS_NOOPT +=  -I$(NCPATH)/include
LDFLAGS = -L$(NCPATH)/lib -Wl,-rpath $(NCPATH)/lib -lnetcdff -lnetcdf

ifeq ($(MICRO_DIR),MICRO_M2005_PA)
  FFLAGS_MICRO = ${FFLAGS_O1}
else
  FFLAGS_MICRO = ${FFLAGS}
endif

#
# Old Lahey compiler setup
    ifeq ($(HOSTTRIM),rex)
      # UNCOMMENT TO USE LAHEY COMPILER -- USEFUL FOR DEBUGGING
      LIB_MPI = /usr/local/mpich-lf/lib
      INC_MPI = /usr/local/mpich-lf/include
      FF77 = lf95
      FF90 = lf95
      CC = gcc -c -DLINUX
      FFLAGS = -g --trap --chk aesu #--trap # --o2 #
      FFLAGS_NOOPT = -g --trap --chk aesu #--trap # --o2 #
      LD = lf95
      LDFLAGS = -L${LIB_MPI} -lmpich
    endif
endif
#--------------------------------------------
# Apple Mac OS X (Darwin) (Absoft Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/absoft/include
#LIB_NETCDF   := /usr/local/absoft/lib
#INC_MPI      := /usr/local/absoft/include
#LIB_MPI       := /usr/local/absoft/lib

#FF77 = f90 -c -f fixed
#FF90 = f90 -c -f free
#CC = cc -c -DMACOSX

#FFLAGS = -O3 -noconsole -nowdir -YEXT_NAMES=LCS -s -YEXT_SFX=_  -z4
#LD = f90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf

#endif

#--------------------------------------------
# Apple Mac OS X (Darwin) (NAG Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/nag/include
#LIB_NETCDF   := /usr/local/nag/lib
#INC_MPI      := /usr/local/nag/include
#LIB_MPI       := /usr/local/nag/lib

#FF77 = f95 -c -fixed -kind=byte
#FF90 = f95 -c -free -kind=byte
#CC = cc -c -DMACOSX

#FFLAGS =      # don't use any optimization -O* option! Will crash!
##FFLAGS = -gline -C=all -C=undefined  # use for debugging
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = f95 
##LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf 
#LDFLAGS =  

#endif


#----------------------------------
# Apple Mac OS X (Darwin) (XLF compiler)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/xlf/include
#LIB_NETCDF   := /usr/local/xlf/lib
#INC_MPI      := /usr/local/xlf/include
#LIB_MPI       := /usr/local/xlf/lib

#FF77 = xlf90 -c -qsuffix=f=f -qfixed=132
#FF90 = xlf90 -c -qsuffix=f=f90
#CC = cc -c -DMACOSX 
#FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = xlf90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

#endif

#----------------------------------

#----------------------------------
# Apple Mac OS X (Darwin) (Intel compiler)
#

ifeq ($(PLATFORM),Darwin)

INC_NETCDF      := /usr/local/netcdf/include
LIB_NETCDF       := /usr/local/netcdf/lib


FF77 = mpif90 -c -fixed -extend_source
FF90 = mpif90 -c 
CC = mpicc -c -DLINUX


FFLAGS = -Os -pad 
#FFLAGS = -g -ftrapuv -check all

FFLAGS += -I${INC_NETCDF}
LD = mpif90 
LDFLAGS = -L${LIB_NETCDF} -lnetcdf 

endif
#
#----------------------------------
# Apple Mac OS X (Darwin) (GNU compiler)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF := /usr/local/include
#LIB_NETCDF := /usr/local/lib
#
#FF77 = gfortran -c -ffixed-form -ffixed-line-length-0
#FF90 = gfortran -c -ffree-form -ffree-line-length-0
#CC = gcc -c -DLINUX
#
#
#FFLAGS = -O3
##FFLAGS = -g -fcheck=all
#
#FFLAGS += -I${INC_NETCDF}
#LD = gfortran
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf
#
#endif


#----------------------------------
#----------------------------------------------
# you dont need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .c .o



all: $(SAM_DIR)/$(SAM)


SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

rrtmg_lw_k_g.o: $(SAM_DIR)/SRC/RAD_RRTM_CFMIP/rrtmg_lw_k_g.f90
	${FF90}  ${FFLAGS_NOOPT} $<

rrtmg_sw_k_g.o: $(SAM_DIR)/SRC/RAD_RRTM_CFMIP/rrtmg_sw_k_g.f90
	${FF90}  ${FFLAGS_NOOPT} $<

microphysics.o: $(SAM_DIR)/SRC/$(MICRO_DIR)/microphysics.f90
	${FF90}  ${FFLAGS_MICRO} $<

.f90.o:
	${FF90}  ${FFLAGS} $<
.f.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<



include Depends



clean: 
	rm ./OBJ/*


