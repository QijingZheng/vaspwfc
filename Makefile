#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
FC= ifort -assume byterecl
# FC= gfortran -I/home/zqj/apps/fftw/3.3.6/include/ -L/home/zqj/apps/fftw/3.3.6/lib/ -lfftw3
FC= gfortran -I/home/zqj/apps/fftw/3.3.6/include/ -lfftw3
FFLAGS= -g -O2
MAKE = make

#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC= prec.f90 lattice.f90 info.f90 wave.f90 gvector.f90 \
	 wavefft.f90 main.f90

OBJ = $(SRC:.f90=.o)
OBJ_PARGAMMA = $(SRC_PARGAMMA:.f90=.o)
EXE = vaspwfc

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
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ) $(SPGLIB)  

clean:
	rm -f *.mod *.a
	rm -f $(OBJ) $(EXE) $(OBJ_PARGAMMA)