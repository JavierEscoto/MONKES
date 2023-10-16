#################################################################################
############################# MONKES Makefile ####################################
##################################################################################

# *** Compiler selection
f90comp=ifort
NETCDF_DIR=/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-17.0.8/netcdf-fortran-4.4.4-7aFYnuhdtreECAnwcH2MYPC7vfd7fUhEpv4r8ohz1zaL
FFLAGS   =-fpp -O3 -g -traceback -xHost -shared-intel -warn all -mcmodel large --with-debugging=yes --with-scalar-type=real --with-precision=double -r8 -start_in_debugger -heap-arrays -qopt-matmul
#
IFLAGS   =-I$(NETCDF_DIR)/include 
LFLAGS   =-L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -L/${MKL_HOME}/lib/intel64 -mkl=parallel -Wl,-rpath,${MKL_HOME}/lib/intel64 
##NETCDF_DIR=/mnt/lustre/ohpc/admin/spack/0.12.1/opt/spack/linux-centos7-x86_64/intel-17.0.8/netcdf-fortran-4.4.4-7aFYnuhdtreECAnwcH2MYPC7vfd7fUhEpv4r8ohz1zaL
##FFLAGS   =-coarray=shared -fpp -DIPPorNIFS -DMPIandPETSc -O3 -g -traceback -xHost -shared-intel -warn all -mcmodel large --with-debugging=yes --with-scalar-type=real --with-precision=double -r8  -start_in_debugger
##IFLAGS   =-I${PETSC_DIR}/lib/petsc/conf/variables -I${PETSC_DIR}/lib/petsc/conf/rules -I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -qopenmp -I$(NETCDF_DIR)/include 
##LFLAGS   =-L${PETSC_DIR}/lib -lpetsc -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -L/${MKL_HOME}/lib/intel64 -mkl -Wl,-rpath,${MKL_HOME}/lib/intel64,-rpath,${PETSC_DIR}/lib

# *** Defining objects for: modules, examples and main in that order
objects = \
Utilities.o \
plots.o \
Sparse_Matrices.o \
Linear_systems.o \
Lagrange_interpolation.o \
Finite_differences.o \
Magnetic_configuration.o \
Barycentric_grids.o \
DKE_BTD_Solution_Legendre.o \
DKE_BTD_Solution_Legendre_Orthonormal.o \
Adjoint_DKE_BTD_Solution_Legendre.o \
API_Example_Barycentric_grids.o \
API_Example_DKE_BTD_Solution_Legendre.o \
API_Example_DKE_BTD_Solution_Legendre_Orthonormal.o \
API_Example_Adjoint_DKE_BTD_Solution_Legendre.o \
API_Example_Coarray.o \
API_Example_DKE_Discretization.o \
main_monkes.o

#gfortran
#-freal-4-real-8 -g -fbacktrace -llapack -lblas  
# *** Executable file
main_monkes.x: $(objects)
	$(f90comp) -o main_monkes.x $(objects) ${COMP_LIB} ${LFLAGS} ${IFLAGS} ${FFLAGS}

# *** Files for libraries
Adjoint_DKE_BTD_Solution_Legendre.o: ../sources/libraries/Adjoint_DKE_BTD_Solution_Legendre.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Adjoint_DKE_BTD_Solution_Legendre.f90

DKE_BTD_Solution_Legendre.o: ../sources/libraries/DKE_BTD_Solution_Legendre.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/DKE_BTD_Solution_Legendre.f90 

DKE_BTD_Solution_Legendre_Orthonormal.o: ../sources/libraries/DKE_BTD_Solution_Legendre_Orthonormal.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/DKE_BTD_Solution_Legendre_Orthonormal.f90 

Barycentric_grids.o: ../sources/libraries/Barycentric_grids.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Barycentric_grids.f90 

Finite_differences.o: ../sources/libraries/Finite_differences.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Finite_differences.f90 

Lagrange_interpolation.o: ../sources/libraries/Lagrange_interpolation.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Lagrange_interpolation.f90 

Linear_systems.o: ../sources/libraries/Linear_systems.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Linear_systems.f90 

Magnetic_configuration.o: ../sources/libraries/Magnetic_configuration.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Magnetic_configuration.f90

plots.o: ../sources/libraries/plots.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/plots.f90

Utilities.o: ../sources/libraries/Utilities.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Utilities.f90  

Sparse_Matrices.o: ../sources/libraries/Sparse_Matrices.f90
	$(f90comp) -c $(FFLAGS) ../sources/libraries/Sparse_Matrices.f90  

# *** Files for examples
API_Example_Coarray.o: ../sources/examples/API_Example_Coarray.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_Coarray.f90

API_Example_Barycentric_grids.o: ../sources/examples/API_Example_Barycentric_grids.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_Barycentric_grids.f90

API_Example_DKE_BTD_Solution_Legendre.o: ../sources/examples/API_Example_DKE_BTD_Solution_Legendre.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_DKE_BTD_Solution_Legendre.f90

API_Example_DKE_BTD_Solution_Legendre_Orthonormal.o: ../sources/examples/API_Example_DKE_BTD_Solution_Legendre_Orthonormal.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_DKE_BTD_Solution_Legendre_Orthonormal.f90

API_Example_Adjoint_DKE_BTD_Solution_Legendre.o: ../sources/examples/API_Example_Adjoint_DKE_BTD_Solution_Legendre.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_Adjoint_DKE_BTD_Solution_Legendre.f90


API_Example_DKE_Discretization.o: ../sources/examples/API_Example_DKE_Discretization.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_DKE_Discretization.f90

# *** File for main
main_monkes.o: ../sources/main_monkes.f90
	$(f90comp) -c $(FFLAGS) ../sources/main_monkes.f90 

# *** Cleaning everything
clean:
	rm main_monkes.x
	rm *.o
	rm *.mod
