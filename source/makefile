# Makefile to compile FDMNES with call to MUMPS, LAPACK and BLAS libraries, for sequential calculations
# MUMPS comes with associated other libraries.
# Works with the gfortran Linux Compiler (for other compiler change the FC command)

FC = gfortran
OPTLVL = 3 

EXEC = ../fdmnes

# For intel compiler, it seems that probems at execution are avoided when compiling
# sphere.f90 with O1 option and the other routines with O2 option.

BIBDIR = bib

FFLAGS = -c  -O$(OPTLVL)  

OBJ = main.o clemf0.o coabs.o convolution.o dirac.o fdm.o fprime.o general.o lecture.o mat.o metric.o \
      minim.o optic.o potential.o selec.o scf.o spgroup.o sphere.o tab_data.o tddft.o tensor.o \
      mat_solve_mumps.o

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) -o $@ $^ -L$(BIBDIR) -ldmumps -lzmumps -lmumps_common -lmpiseq -lmetis -lpord \
                                   -lesmumps -lscotch -lscotcherr -lpthread -llapack -lblas 

%.o: %.f90
	$(FC) -o $@ $(FFLAGS) $? 

clean:
	rm -f *.o $(EXEC)
	rm -f *.mod	

