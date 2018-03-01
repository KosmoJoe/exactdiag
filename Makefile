 #############################################################################################
 ######## Makefile 
 #############################################################################################
 ######## 26.7.2013 Johannes Schurer
 #############################################################################################

PROG =  threeBody.x

PATH_OBJ=./object
PATH_SRC=./src

SRCS := module.f90 main.f90 swap.f90 grid.f90 saveIn.f90
OBJS := $(addprefix $(PATH_OBJ)/,$(patsubst %.f90, %.o, $(SRCS)))

LIBS =	 -lmkl_lapack95_lp64  -lmkl_intel_lp64  -lmkl_core  -lmkl_sequential -lpthread

F90 = ifort
F90FLAGS = -O3 -module $(PATH_OBJ)
MKLPATHSIS = /afs/physnet.uni-hamburg.de/software/softsrv/intel/mkl/10.0.011/lib/em64t
LDFLAGS =     

OUTERR= > comp.log 2>&1

#-Wl,-rpath,$(MKLPATHSIS) 
#-openmp -check bounds -parallel

all: $(PROG)

$(PROG): $(OBJS)
	@echo -n Linking $@ ...
	@$(F90) $(LDFLAGS) -o  $@ $(OBJS) -L $(LD_LIBRARY_PATH) $(LIBS)
	@echo ' done'
	
clean:
	@echo -n Cleaning ...
	@rm -f $(PROG) $(OBJS) $(PATH_OBJ)/*.mod
	@echo ' done'
	
.SUFFIXES: $(SUFFIXES) .f90

$(PATH_OBJ)/%.o: $(PATH_SRC)/%.f90
	@echo -n Compiling $< ...
	@$(F90) $(F90FLAGS) -c -o $@ $< $(OUTERR)
	@echo ' done'
	

