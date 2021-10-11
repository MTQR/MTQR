#---------------------------------------------------------------------------------------
# File:      build/makefile
#
# Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
#                      a C++ library for high precision integration of singular 
#                      polynomials of non-integer degree
#
# Authors:   Guido Lombardi, PhD, Davide Papapicco
#
# Institute: Politecnico di Torino
#            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
#            Department of Electronics and Telecommunications (DET)
#            Electromagnetic modelling and applications Research Group
#---------------------------------------------------------------------------------------

CXX = g++

INCDIR = ../include
SRCDIR = ../src
OBJDIR = obj

OBJ = $(OBJDIR)/Quasimont.o $(OBJDIR)/DatIo.o $(OBJDIR)/MonMap.o

CXXOPTIONS = -g
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR) 

default: $(OBJ)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/*.h
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(OBJDIR)/*.o