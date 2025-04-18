#################################################################################
############################# MONKES Makefile ####################################
##################################################################################

# *** Compiler selection
f90comp=gfortran
NETCDF_DIR=/usr/include/ 

FFLAGS=-I$(NETCDF_DIR) -mcmodel=large -O2 -freal-4-real-8 -ffree-line-length-512 
LFLAGS=-llapack -lnetcdf -lnetcdff -lblas 

# *** Defining objects for: modules, examples and main in that order
objects = \
Utilities.o \
Sparse_Matrices.o \
Linear_systems.o \
Lagrange_interpolation.o \
Finite_differences.o \
Magnetic_configuration.o \
Barycentric_grids.o \
DKE_BTD_Solution_Legendre.o \
API_Example_DKE_BTD_Solution_Legendre.o \
main_monkes.o

#gfortran
#-freal-4-real-8 -g -fbacktrace -llapack -lblas  
# *** Executable file
main_monkes.x: $(objects)
	$(f90comp) -o main_monkes.x $(objects) ${LFLAGS} ${FFLAGS}

# *** Files for libraries 

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

API_Example_Barycentric_grids.o: ../sources/examples/API_Example_Barycentric_grids.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_Barycentric_grids.f90

API_Example_DKE_BTD_Solution_Legendre.o: ../sources/examples/API_Example_DKE_BTD_Solution_Legendre.f90
	$(f90comp) -c $(FFLAGS) ../sources/examples/API_Example_DKE_BTD_Solution_Legendre.f90
   
# *** File for main
main_monkes.o: ../sources/main_monkes.f90
	$(f90comp) -c $(FFLAGS) ../sources/main_monkes.f90 

# *** Cleaning everything
clean:
	rm *.o
	rm *.mod
	rm main_monkes.x
