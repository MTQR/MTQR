#---------------------------------------------------------------------------------------
# File:      build/linux/makefile
#
# Library:   MTQR - Monomial Transformation Quadrature Rule:
#                   a C++ library for high-precision integration of 
#                   generalised polynomials of non-integer degree
#
# Authors:   Guido Lombardi, Davide Papapicco
#
# Institute: Politecnico di Torino
#            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
#            Department of Electronics and Telecommunications (DET)
#            Electromagnetic modelling and applications Research Group
#---------------------------------------------------------------------------------------

CXX = g++

INCDIR = ../../include
SRCDIR = ../../src
OBJDIR = obj

OBJ = $(OBJDIR)/mtqr.o $(OBJDIR)/data_management.o $(OBJDIR)/monomial_transformation.o $(OBJDIR)/vector_operations.o

CXXOPTIONS = -g -ansi -std=c++17 #-fext-numeric-literals
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR) 

dafault: $(OBJ)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/*.h
	@mkdir obj
	$(CXX) -c $< -o $@ $(CXXFLAGS)

static: $(OBJ)
	ar rcs ../static/libmtqr.a $^

clean:
	rm -f $(OBJDIR)