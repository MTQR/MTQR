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

OBJ = $(OBJDIR)/Utils.o $(OBJDIR)/MonMap.o

CXXOPTIONS = -g -ansi -std=c++17
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR) 

default: $(OBJ)

$(OBJDIR)/Utils.o: Utils.cpp $(INCDIR)/*.h
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/MonMap.o:
	cd $(SRCDIR) && $(MAKE)

clean:
	rm -f $(OBJDIR)/*.o