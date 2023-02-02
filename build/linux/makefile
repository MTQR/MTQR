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

-include macros.mk

.PHONY = all install library clean uninstall 

default: install library clean

all: install library tests clean

install: $(OBJ)

library:
	ar rcs $(LIB) $(OBJ)
	@echo ""

tests:
	@echo "Compiling tests..."
	@$(foreach dir, $(wildcard $(TST_DIR)), $(CC) -o $(notdir $(dir)) $(dir)/main.cpp $(INCLUDE) -L$(LIB_DIR) $(LDFLAGS);)
	@echo "...done!"
	@echo ""

clean:
	@echo "Cleaning object files..."
	@rm -r $(OBJ_DIR)
	@echo "...done!"

uninstall:
	@rm -r $(LIB)

$(OBJ_DIR)/mtqr.o: $(SRC_DIR)/mtqr.cpp
	@mkdir $(OBJ_DIR)
	$(CC) -c $^ -o $@ $(CCFLAGS)

$(OBJ_DIR)/data_management.o: $(SRC_DIR)/data_management.cpp
	$(CC) -c $^ -o $@ $(CCFLAGS)

$(OBJ_DIR)/monomial_transformation.o: $(SRC_DIR)/monomial_transformation.cpp
	$(CC) -c $^ -o $@ $(CCFLAGS)

$(OBJ_DIR)/vector_operations.o: $(SRC_DIR)/vector_operations.cpp
	$(CC) -c $^ -o $@ $(CCFLAGS)