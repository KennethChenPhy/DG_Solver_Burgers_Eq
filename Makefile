.SUFFIXES :
.SUFFIXES : .f90 .o
SYSTEM = 
FC = gfortran 
LAPATH = /home/kennethchenphy/lapack-3.8.0
LIB = $(LAPATH)/liblapack.a $(LAPATH)/librefblas.a
F_SRC = Global_var.o Global_math.o Nodal1D.o Grid1D.o SlopeLimiter.o\
Burgers1D.o DG_BurgersDriver1D.o
FFLAGS =  
.f90.o:
	$(FC) -c $(FFLAGS) $<
Executable : $(F_SRC)
	$(FC) -o Executable $(FFLAGS) $(F_SRC) $(LIB) 
Global_var.o: Global_var.f90
Global_math.o:Global_math.f90
Nodal1D.o: Nodal1D.f90
Grid1D.o: Grid1D.f90
SlopeLimiter.o:SlopeLimiter.f90
Burgers1D.o: Burgers1D.f90
DG_BurgersDriver1D.o: DG_BurgersDriver1D.f90
clean :
	rm -rf *.o
	rm -rf Executable
	rm -rf *.dat
	rm -rf *.mod
