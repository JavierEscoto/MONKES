#################################################################################
############################# MONKES Makefile ####################################
##################################################################################

# *** Compiler selection
f90comp = mpiifort

##################################################################################
# NetCDF
##################################################################################

NETCDF_INCLUDE = /mpcdf/soft/SLE_15/packages/skylake/netcdf-serial/intel_21.7.1-2021.7.1/4.4.1/include
NETCDF_LIB     = /mpcdf/soft/SLE_15/packages/skylake/netcdf-serial/intel_21.7.1-2021.7.1/4.4.1/lib

IFLAGS = -I$(NETCDF_INCLUDE)
LFLAGS = -L$(NETCDF_LIB) -lnetcdff -lnetcdf -mkl=parallel -Wl,-rpath,$(NETCDF_LIB)

##################################################################################
# Compilation flags
##################################################################################

FFLAGS = -mcmodel=large -fpp -O3 -g -traceback -xHost -shared-intel -warn all \
         -r8 -heap-arrays -qopt-matmul

##################################################################################
# Objects
##################################################################################

objects = \
quadpack.o \
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

##################################################################################
# Executable
##################################################################################

main_monkes.x: $(objects)
	$(f90comp) -o $@ $(objects) $(LFLAGS) $(IFLAGS) $(FFLAGS) -Wl,-z,noexecstack

##################################################################################
# Rules
##################################################################################

%.o: ../sources/libraries/%.f90
	$(f90comp) -c $(FFLAGS) $(IFLAGS) $<

%.o: ../sources/libraries/%.f
	$(f90comp) -c $(FFLAGS) $(IFLAGS) $<

API_Example_%.o: ../sources/examples/API_Example_%.f90
	$(f90comp) -c $(FFLAGS) $(IFLAGS) $<

main_monkes.o: ../sources/main_monkes.f90
	$(f90comp) -c $(FFLAGS) $(IFLAGS) $<

##################################################################################
# Cleaning
##################################################################################

clean:
	rm -f *.o *.mod main_monkes.x

