#---------------------------------------------------------------------------------------
# File:      utilities/makefile
#
# Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
#                      a C++ library for high precision integration of singular 
#                      polynomials of non-integer degree
#
# Authors:   Guido Lombardi, Davide Papapicco
#
# Institute: Politecnico di Torino
#            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
#            Department of Electronics and Telecommunications (DET)
#            Electromagnetic modelling and applications Research Group
#---------------------------------------------------------------------------------------

CXX = g++

BLDDIR = ../build
INCDIR = ../include
SRCDIR = ../src
OBJDIR = $(BLDDIR)/obj

OBJ = $(OBJDIR)/VecOps.o $(OBJDIR)/MonMap.o

CXXOPTIONS = -g -ansi -std=c++17
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR) 

default: $(OBJ)

$(OBJDIR)/VecOps.o: VecOps.cpp $(INCDIR)/*.h
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/MonMap.o:
	cd $(SRCDIR) && $(MAKE)

tools: ErrTools.cpp $(INCDIR)/Utils.h
	$(CXX) -o Tools $< -lmpfr -I $(INCDIR)
	./Tools

clean:
	rm -f Tools && if [ -d "estimate" ]; then rm -r estimate; fi
	rm -f $(OBJDIR)/*.o