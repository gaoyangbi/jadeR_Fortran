FC = ifx
OMPFC = mpif90
MPIFC = mpiifort
LIBPATH = -L/data1/share/software/intel/oneapi/mkl/2022.0.2/lib/intel64 
SCALAPACK = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
LIBS =  -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64  -lmkl_intel_thread  -lmkl_core -lpthread -liomp5
#COPTS =   -qopenmp -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000
COPTS =  -fopenmp
#all: shale01 shale02 shale03 shale04 shale05 shale06 shalecg02 shalecg03 shalecg04 shalecg05 shalecg06 shalecg07 shalecg08
all: jadeR_Fortran
#ifeq ($(FC),gfortran)
#	COPTS = -fopenmp
#	LIBS = -llapack -lblas
#endif
jadeR_Fortran: jadeR_Fortran.f90
	$(FC) -o jadeR_Fortran *.f90 -qmkl -static-intel $(LIBPATH) $(LIBS) 

clean:
	rm -f *.o jadeR_Fortran

