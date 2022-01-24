#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
FC= ifort -assume byterecl
# FC= gfortran

FFLAGS= -g -O2

#LDFLAGS= -L/home/zqj/apps/fftw/3.3.6/lib/ 
#LIBS= -lfftw3
#FFLAGS+= -I/home/zqj/apps/fftw/3.3.6/include/

# Have problem compiling with intel psxe/onemkl? uncomment following 3 lines
#LDFLAGS=  -L${MKLROOT}/lib/intel64 -L${MKLROOT}/interfaces/fftw3xf/
#LIBS= -lpthread -lm -ldl -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#FFLAGS+= -I${MKLROOT}/include/fftw/

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
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS)

gam:	$(OBJ_GAM)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXE_gam) $(OBJ_GAM) $(LIBS)

clean:
	rm -f *.mod *.a *.o
	rm -f $(OBJ) $(OBJ_GAM)
