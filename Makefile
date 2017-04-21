#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
FC= ifort -assume byterecl
# FC= gfortran -I/home/zqj/apps/fftw/3.3.6/include/ -L/home/zqj/apps/fftw/3.3.6/lib/ -lfftw3
# FC= gfortran -I/home/zqj/apps/fftw/3.3.6/include/ -lfftw3
FC= gfortran -lfftw3
FFLAGS= -g -O2
MAKE = make

#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC= prec.f90 lattice.f90 info.f90 wave.f90 gvector.f90 \
	 wavefft.f90 main.f90

SRC_GAM= prec.f90 lattice.f90 info.f90 wave.f90 gvector_gam.f90 \
	 wavefft_gam.f90 main_gam.f90

OBJ = $(SRC:.f90=.o)
OBJ_GAM = $(SRC_GAM:.f90=.o)
EXE = vaspwfc
EXE_gam = vaspwfc_gam

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
wfc:	$(OBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ) 

gam:	$(OBJ_GAM)
	$(FC) $(FFLAGS) -o $(EXE_gam) $(OBJ_GAM) 

clean:
	rm -f *.mod *.a *.o
	rm -f $(OBJ) $(OBJ_GAM)
