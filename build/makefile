CXX = g++

INCDIR = ../include
SRCDIR = ../src
OBJDIR = obj

OBJ = $(OBJDIR)/hiperquad.o $(OBJDIR)/datio.o $(OBJDIR)/monmap.o

CXXOPTIONS = -g
CXXFLAGS = $(CXXOPTIONS) -I $(INCDIR) 

LDLIBS = -lm -lgsl -lgslcblas -lgmp -lquadmath
LDFLAGS = $(LDLIBS)

default: $(OBJ)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/*.h
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(OBJDIR)/*.o

