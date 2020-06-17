# Makefile to compile FDMNES with call to MUMPS, LAPACK and BLAS libraries, for sequential calculations
# MUMPS comes with associated other libraries.
# Works with the gfortran Linux Compiler (for other compiler change the FC command)

FC = gfortran-mp-6
OPTLVL = 3 

EXEC = ../fdmnes

# For intel compiler, it seems that probems at execution are avoided when compiling
# sphere.f90 with O1 option and the other routines with O2 option.

BIBDIR = /opt/local/lib

INCDIR = include

FFLAGS = -c  -O$(OPTLVL) -I$(INCDIR)

OBJ = main.o clemf0.o coabs.o convolution.o diffraction.o dirac.o fdm.o fdmx.o fprime.o fprime_data.o general.o lecture.o mat.o metric.o \
      minim.o optic.o potential.o selec.o scf.o spgroup.o sphere.o tab_data.o tddft.o tensor.o \
      mat_solve_mumps.o RIXS.o

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) -o $@ $^  -L$(BIBDIR) -ldmumps -lzmumps -lmumps_common -lmpiseq  -lesmumps -lmetis -lpord \
                                   -lscotch -lscotcherr -lpthread -llapack -lblas 
#$(EXEC): $(OBJ)
	#$(FC) -o $@ $^ -L$(BIBDIR) -ldmumps -lzmumps -lmumps_common -lesmumps -lmetis -lpord \
                                   -lscotch -lscotcherr -lpthread -llapack -lblas 
#$(EXEC): $(OBJ)
#	$(FC) -o $@ $^ -L$(BIBDIR) -ldmumps -lzmumps -lmumps_common -lpord \
                                   -lscotch -lscotcherr -lpthread -llapack -lblas 

%.o: %.f90
	$(FC) -o  $@ $(FFLAGS) $? 

clean:
	rm -f *.o $(EXEC)
	rm -f *.mod	

