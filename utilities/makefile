#------------------------------------------------------------------------------------------------------------------
# File:      utilities/makefile
#
# Library:   QUASIMONT - QUAdrature of SIngular polynomials using MONomial Transformations:
#                        a C++ library for high precision quadrature of singular polynomials of non-integer degree
#
# Authors:   Davide Papapicco, Guido Lombardi, PhD
#
# Institute: Politecnico di Torino
#            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
#            Department of Electronics and Telecommunications (DET)
#            Electromagnetic modelling and applications Research Group
#------------------------------------------------------------------------------------------------------------------

CXX = g++

BLDDIR = ../build
INCDIR = ../include
SRCDIR = ../src
OBJDIR = $(BLDDIR)/obj

OBJ = $(OBJDIR)/VectorManip.o $(OBJDIR)/Plot.o $(OBJDIR)/MonMap.o

CXXOPTIONS = -g
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR) 

default: $(OBJ)

$(OBJDIR)/VectorManip.o: VectorManip.cpp $(INCDIR)/*.h
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/Plot.o: Plot.cpp $(SRCDIR)/MonMap.cpp $(INCDIR)/*.h
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/MonMap.o:
	cd $(SRCDIR) && $(MAKE)

clean:
	rm -f $(OBJDIR)/*.o