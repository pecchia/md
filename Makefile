
FC = gfortran
#FC = ifort
FLAGS = -O3 -g -fmax-errors=3
OMP = -fopenmp
#FC = ifort 
#FLAGS = -g -CB -check all 
LAPACKDIR = ../../../../lapack/lapack-3.12.0/
#LAPACKDIR = /home/paolo/Scaricati/lapack-3.12.0
LAPACK = -L$(LAPACKDIR) -llapack -lrefblas
#MKL_LIBDIR=/home/pecchia/mkl/
#LAPACK=-L$(MKL_LIBDIR) -Wl,-rpath,$(MKL_LIBDIR) -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
LIBS=$(LAPACK)

TARGET1 = md 

SOURCES1 = precision.f90 \
         constants.f90 \
	   parameters.f90 \
	   list.f90 \
	   boxes.f90 \
	   input.f90 \
	   forces.f90 \
	   simulations.f90 \
	   dynamics.f90 \
	   clock.f90 \
	   md.f90 \
	   lyap_n.f90 \
	   vector.f90
OBJS1 = $(SOURCES1:.f90=.o)

%.o: %.f90 
	$(FC) $(OMP) $(FLAGS) -c  $<

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(FC) $(OMP) -o $(TARGET1) $(OBJS1) $(LIBS)


test: list.o clock.o
	$(FC) $(FLAGS) -c test.f90
	$(FC) -o test test.o list.o clock.o



clean:
	rm *.o *.mod $(TARGET1) 

constants.o : precision.o
parameters.o : constants.o
boxes.o : constants.o
lyap_n.o : constants.o parameters.o
vector.o : constants.o parameters.o
simulations.o : constants.o parameters.o forces.o dynamics.o clock.o vector.o
forces.o : constants.o list.o boxes.o 	
dynamics.o : constants.o parameters.o boxes.o lyap_n.o	
md.o : constants.o parameters.o list.o boxes.o input.o simulations.o
