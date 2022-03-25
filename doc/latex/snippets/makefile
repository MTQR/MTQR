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