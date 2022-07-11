#! /usr/bin/make

prep   = cencalc_prep
ccmla  = cencalc_ccmla
all : $(prep) $(ccmla)

# DISLIN: scientific data plotting software at https://www.dislin.de/
#
# DISLIN path. Make sure that dislin.mod is precompiled
# with the same compiler version used here. Add DISLIN to LD_LIBRARY_PATH for runntime
#  
DISLIN = /opt/apps/SL7/dislin/

# GNU fortran version 8.3.1 or higher is required
FCOMPL = gfortran
FFLAGC = -fopenmp  -I $(DISLIN)/gf/real64/ 

# For production 
FFLAGC_EXTRA = -O3 -funroll-loops -finline -finline-functions -ftree-vectorize -ffast-math  

# For debugging 
#FFLAGC_EXTRA = -fbacktrace -finit-real=nan -g -pg -fcheck=all -fdump-core -fmax-errors=1 -ffpe-trap=invalid -Wconversion-extra
# 
LDFLAG_PREP = $(opt) -fopenmp  -L $(DISLIN) -ldislin_d
LDFLAG = $(opt) -fopenmp  

# Implicit rules
.SUFFIXES: .vec .f .c .s .o .f90

.f90.o:
	$(FCOMPL) -c $(FFLAGC) $(FFLAGC_EXTRA) -o $@ $<

# Object files
OBJECTS_CCMLA = parameters.o qsort.o cencalc_ccmla.o 
OBJECTS_PREP = parameters.o cencalc_prep.o 

# Clean objects
clean:
	@rm -f core $(OBJECTS_CCMLA) $(OBJECTS_PREP) *.mod $(obj)

# Program targets
$(prep): $(OBJECTS_PREP) 
	$(FCOMPL) $(FFLAGC) -o $(prep) $(OBJECTS_PREP) $(LDFLAG_PREP) 

$(ccmla): $(OBJECTS_CCMLA) 
	$(FCOMPL) $(FFLAGC) -o $(ccmla) $(OBJECTS_CCMLA) $(LDFLAG) 

