CXX = g++

BLDDIR = ../build
INCDIR = ../include
SRCDIR = ../src
OBJDIR = $(BLDDIR)/obj

OBJ = $(OBJDIR)/main.o $(OBJDIR)/hiperquad.o

CXXOPTIONS = -g -Wall
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR)

LDLIBS = -lm -lgsl -lgslcblas -lgmp -lquadmath
LDFLAGS = $(LDLIBS)

default: $(OBJ)
	$(CXX) $(OBJDIR)/*.o -o HIPERquad $(LDFLAGS)

$(OBJDIR)/main.o: main.cpp
	cd $(BLDDIR) && $(MAKE) clean
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/hiperquad.o:
	cd $(BLDDIR) && $(MAKE)

clean:
	rm -f $(OBJDIR)/main.o HIPERquad