Here are the different fdmnes routines:
   main.f90
   clemf0.f90
   coabs.f90
   convolution.f90
   dirac.f90
   fdm.f90
   fprime.f90
   general.f90
   leture.f90
   mat.f90
   mat_solve_mumps.f90
   metric.f90
   minim.f90
   optic.f90
   potential.f90
   scf.f90
   selec.f90
   spgroup.f90
   sphere.f90
   tab_data.f90
   tddft.f90
   tensor.f90

They must be compiled and linked together with the MUMPS, LAPACK and BLAS library.

If no MUMS library are avalable
   mat_solve_mumps.f90 must be replaced by mat_solve_gaussian.f90
   With sequential calculations not_mpi.f90 must also be compiled and linked.

   when no LAPAK and BLAS libraries are avalable, they can be replaced by
   sub_util.f (only without MUMPS)

IMPORTANT !

  The *.h files are files which are automatically included during the compilation.

  These files are for the MUMPS version 5.0.0
  When using MUMPS library already linked on your server, this can give an uncompatibility
  (misalignement...) because MUMPS was compiled with another version of these .h files. 
  In this case ask the informatician the same .h files or download them from the MUMPS site.


Makefile is an example of very simple makefile for sequential code using MUMPS library.
Makefile_gaussian is an example of very simple makefile for sequential code using
Gaussian solver without any library.

For intel compiler, it seems that probems at execution are avoided when compiling
sphere.f90 with O1 option in place of O2 option.
