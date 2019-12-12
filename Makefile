.SUFFIXES :
.SUFFIXES : .f90 .o
SYSTEM = `uname`
FC = ifort
FFLAGS = -mkl -check all
F_SRC = Global_var.o Global_math.o Nodal1D.o Grid1D.o SlopeLimiter.o\
	     Advect1D.o Burgers1D.o DG_BurgersDriver1D.o
.f90.o:
	$(FC) -c $(FFLAGS) $<
Executable : $(F_SRC)
	$(FC) -o Executable $(FFLAGS) $(F_SRC)
Global_var.o: Global_var.f90
Global_math.o:Global_math.f90
Nodal1D.o: Nodal1D.f90
Grid1D.o: Grid1D.f90
SlopeLimiter.o:SlopeLimiter.f90
Advect1D.o:Advect1D.f90
Burgers1D.o: Burgers1D.f90
DG_BurgersDriver1D.o: DG_BurgersDriver1D.f90
clean :
	rm -rf *.o
	rm -rf Executable
	rm -rf *.dat
	rm -rf *.mod
