#---------------------------------------------------------------------------------------
# File:      test/makefile
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
UTLDIR = ../utilities
OBJDIR = $(BLDDIR)/obj

OBJ = $(OBJDIR)/Main.o $(OBJDIR)/Quasimont.o $(OBJDIR)/VecOps.o

CXXOPTIONS = -ansi -std=c++17 #-fext-numeric-literals
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR)

LDLIBS = -lm -lquadmath -lgsl -lgslcblas
LDFLAGS = $(LDLIBS)

default: $(OBJ)
	$(CXX) $(OBJDIR)/*.o -o Test $(LDFLAGS)
	./Test

$(OBJDIR)/Main.o: Main.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/Quasimont.o:
	cd $(BLDDIR) && $(MAKE)

$(OBJDIR)/VecOps.o:
	cd $(UTLDIR) && $(MAKE) $@

tools:
	cd $(UTLDIR) && $(MAKE) $@

clean:
	rm -f $(OBJDIR)/Main.o Test
	if [ -d "output" ]; then rm -r output; fi 

allclean:
	$(MAKE) clean
	cd $(UTLDIR) && $(MAKE) clean
	cd $(BLDDIR) && $(MAKE) clean