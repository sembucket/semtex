# GNUmakefile

# G95
#FC = g95
#FOPTS = -ffree-form -ftrace=full

# Intel Fortran Compiler
FC = ifort
FOPTS = -FR -O -traceback

MODULES    = \
 constants.o \
 basic_io.o \
 spectral_element_vertex.o \
 boundary_curve.o \
 boundary_condition.o \
 spectral_element_edge.o \
 spectral_element.o \
 boundary_segment.o \
 export_functions.o \
 data_structure.o \
 import_grid.o \
 control_data.o

PROGRAM = gambit2semtex

.PHONY: default clean tidy

default: $(PROGRAM)

$(PROGRAM): $(MODULES) $(PROGRAM).o
	$(FC) -o $@ $@.o $(MODULES)

clean:
	rm -f *.o *.mod *~

tidy:   clean
	rm -f $(PROGRAM)

# rule for objects
%.o : %.f
	$(FC) $(FOPTS) -c $<
