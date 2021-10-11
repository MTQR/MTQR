#---------------------------------------------------------------------------------------
# File:      test/makefile
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

BLDDIR = ../build
INCDIR = ../include
SRCDIR = ../src
UTLDIR = ../utilities
OBJDIR = $(BLDDIR)/obj

OBJ = $(OBJDIR)/Main.o $(OBJDIR)/Quasimont.o $(OBJDIR)/Utils.o

CXXOPTIONS = -g #-Wall
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR)

LDLIBS = -lm -lgsl -lgslcblas -lgmp -lquadmath
LDFLAGS = $(LDLIBS)

default: $(OBJ)
	$(CXX) $(OBJDIR)/*.o -o Test $(LDFLAGS)
	./Test

$(OBJDIR)/Main.o: Main.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/Quasimont.o:
	cd $(BLDDIR) && $(MAKE)

$(OBJDIR)/Utils.o:
	cd $(UTLDIR) && $(MAKE) $@

clean:
	rm -f $(OBJDIR)/Main.o Test

allclean:
	cd $(BLDDIR) && $(MAKE) clean