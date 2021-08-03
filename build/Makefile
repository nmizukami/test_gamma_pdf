#========================================================================
# Makefile to compile the routing model
#========================================================================
#
# Need to edit FC, FC_EXE, EXE, isOpenMP (if OMP is used), F_MASTER, NCDF_PATH
# Feel free to modify compiler flag (FLAGS)
# Do not leave space after your edits
#
#========================================================================
# User configure part
#========================================================================
# Define fortran compiler - gnu, intel or pgi
FC  = 

# Define the compiler exe, e.g., gnu=>gfortran, intel=>ifort, pgi=>pgf90
FC_EXE = 

# Define the compiled executable
EXE = 

# Define optional setting
# fast:      Enables optimizations
# debug:     Minimum debug options, still
# profile:   Enables profiling
MODE = debug

# Define core directory below which everything resides
# parent directory of the 'build' directory
# do not put space at the end of path
F_MASTER = 

# Define the NetCDF libraries and path to include files
# currently no netcdf library needed
ifeq "$(FC)" "gnu"
 NCDF_PATH = 
endif
ifeq "$(FC)" "intel"
 NCDF_PATH = 
endif
ifeq "$(FC)" "pgi"
 NCDF_PATH =
endif

LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
            -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
INCNETCDF = -I$(NCDF_PATH)/include

#========================================================================
# Define flags
#========================================================================
LIBOPENMP =

ifeq "$(FC)" "gnu"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -fmax-errors=0 -ffree-line-length-none
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -Wall -fmax-errors=0 -fbacktrace -fcheck=all -ffpe-trap=zero -ffree-line-length-none
  endif
  ifeq "$(MODE)" "profile"
    FLAGS = -p -g -Wall -fmax-errors=0 -fbacktrace -fcheck=all -ffpe-trap=zero -ffree-line-length-none
  endif
endif

ifeq "$(FC)" "intel"
  ifeq "$(MODE)" "fast"
    FLAGS = -O3 -FR -auto -fltconsistency -fpe0
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -debug all -warn all -check all -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency
  endif
  ifeq "$(MODE)" "profile"
    FLAGS = -p -g -debug -warn all -check all -noerror_limit -FR -O0 -auto -WB -fpe0 -traceback -fltconsistency
  endif
endif

ifeq "$(FC)" "pgi"
  ifeq "$(MODE)" "fast"
    FLAGS = -fast -O3 -Mdclchk
  endif
  ifeq "$(MODE)" "debug"
    FLAGS = -g -Mbounds -Mlist -Minfo -Mdclchk
  endif
endif

#========================================================================
# Define directory paths
#========================================================================
# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)build/src/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)build/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)

#========================================================================
# Assemble all of the sub-routines
#========================================================================
# Define subroutines
MODU = \
    nrtype.f90 \
    nr_utility.f90 \
    ascii_util.f90 \
    gamma_func.f90 \
    read_param.f90 \
    comp_gamma.f90

# concatanate model subroutines
TEMP_MODSUB = $(MODU)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# Define driver programs
#========================================================================
# Define the driver routine
DRIVER = $(F_KORE_DIR)main.f90

#========================================================================
# Compile the puppy
#========================================================================
# tasks
all: compile install clean

# compile
compile:
	$(FC_EXE) $(FLAGS) $(MODSUB) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) -o $(EXE)

# Remove object files
clean:
						rm -f *.o
						rm -f *.lst
						rm -f *.mod
						rm -f *__genmod.f90
#
# Copy the executable to the bin directory
install:
	@mkdir -p $(EXE_PATH)
	@mv $(EXE) $(EXE_PATH)
