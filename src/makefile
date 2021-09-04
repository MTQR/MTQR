CC = g++
TARGET = main

IDIR =../include
LDIR = ../utilities
ODIR = obj
LIBS = -lm -lgsl -lgslcblas -lgmp -lquadmath

CFLAGS = -Wall --std=gnu++11 -I$(IDIR) -I$(LDIR)

_DEPENDS = Prototypes.h GenerateTabValues.h # insert all other header files here
DEPENDS = $(patsubst %,$(IDIR)/%,$(_DEPENDS))

_OBJ = main.o # insert all other source files here
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPENDS) # This rule applies to all the .o files (%.o) in the obj direcotory ($(ODIR)/). It compiles all the .cpp files in the dependencies list (OBJ and DEPENDS) to generate the obj file (-c flag)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TARGET) $(ODIR)/*.o *~ core $(IDIR)/*~