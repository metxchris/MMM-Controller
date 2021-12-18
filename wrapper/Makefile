#======================================
#    MMM - Makefile for mmm_wrapper
#======================================

lib_mmm    := ../src/libmmm.a
libraries  := $(lib_mmm)
f90module  := -I../src
mmmwrapper := mmm

# All object files are listed here, modules are always built first
objects = mmm_wrapper.o

# Compiler and flags
MMMFC  ?= gfortran
FFLAGS := -O3 -march=native -cpp $(f90module)

all : $(mmmwrapper)

$(lib_mmm) :
	cd ../src; make; cd -

$(mmmwrapper) : $(libraries) $(objects)
	$(MMMFC) $(FFLAGS) -o $(mmmwrapper) $(objects) $(libraries)

# Replace optimization flag with debugging flags
# TODO: How to make this use 'make debug' in src folder instead of just 'make'?
.PHONY : debug
debug : FFLAGS := $(subst -O3,-g -ffpe-trap=invalid -ffpe-trap=zero -ffpe-trap=overflow,$(FFLAGS))
debug : $(mmmwrapper)

# Clean local and src directories
.PHONY : cleanall
cleanall :
	rm -f $(mmmwrapper) $(objects) *.MOD *.mod *.o *.a; cd ../src; make clean; cd -

# Clean local directory only
.PHONY : clean
clean :
	rm -f $(mmmwrapper) $(objects) *.MOD *.mod *.o *.a

%.o : %.f
	$(MMMFC) -c $(FFLAGS) $< -o $@

%.o : %.f90
	$(MMMFC) -c $(FFLAGS) $< -o $@

# Compiler Flag Options:
# -g: produces symbolic debugging info and may disable optimization (only use for make debug)
# -cpp: enables preprocessing (conditional compilation) for any kind of Fortran file
# -O2, -O3: optimization levels
# -fbacktrace: (not needed) default flag that enables backtraces during critical errors
# -fdollar-ok: enables feature of older compilers and should be avoided
# -fno-common: (not needed) default flag concerning the placement of global variables without an initializer
# -fno-range-check: Disables compile errors due to division by 0 and replaces with +-inf
#                   Also causes integer overflow errors to wrap around
# -ffpe-trap=invalid,zero,overflow  # Floating point exception traps, for debugging
