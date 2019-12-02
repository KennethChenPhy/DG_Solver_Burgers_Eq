CP = cp
MV = mv
CD = cd
MAKE = make

# c++ compiler
CXX = g++

INCLUDES = -I./Include
#CXXOPTIONS = -c
#CXXFLAGS = $(CXXOPTIONS) $(OPTFLAGS) $(INCLUDES)
CXXFLAGS = $(INCLUDES)

.SUFFIXES: .cpp

OBJS = \
	Src/Codes1D/GradJacobiP.o    \
	Src/Codes1D/JacobiGL.o        \
	Src/Codes1D/JacobiGQ.o         \
	Src/Codes1D/JacobiP.o           \
	Src/Codes1D/Vandermonde1D.o      \

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $*.o -c $*.cpp
Executable:
	$(CXX) -o Executable $(CXXFLAGS) $(OBJS)

clean:
	rm -rf $(OBJS)
