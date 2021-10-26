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
UTLDIR = ../utilities
OBJDIR = obj

UTL = $(UTLDIR)/Utils.o
OBJ = $(OBJDIR)/Quasimont.o $(OBJDIR)/DatIo.o $(OBJDIR)/MonMap.o

CXXOPTIONS = -g -ansi -std=c++17
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR) 

dafault: $(OBJ)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/*.h
	$(CXX) -c $< -o $@ $(CXXFLAGS)

static: $(OBJ) $(OBJDIR)/Utils.o
	ar rcs libquasimont.a $^

$(OBJDIR)/Utils.o: $(UTLDIR)/Utils.cpp $(INCDIR)/Utils.h
	cd $(UTLDIR) && $(MAKE)

clean:
	rm -f $(OBJDIR)/*.o