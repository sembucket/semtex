#
# Spectral Element Library (speclib)
#

include $(PRISM)/system/defaults

LIB      = libspec.a
LIBDIR   = $(PRISM)/lib
INCDIR   = $(PRISM)/include
CPPFLAGS = -I$(PRISM)/support

# Module names and headers

ALL      = speclib   utilities  frame   manager  matrix \
           solve     operators  splib   isomesh  session \
	   solve_cg  fieldfile  tree    family   bc \
	   probe     residual   field   edge

HEADER   = speclib.h \
	bc.h config.h edge.h element.h family.h field.h \
	fieldfile.h frame.h isomesh.h matrix.h operators.h \
	probe.h splib.h tree.h 

# Generate rules for building the modules 

LMOD     = $(foreach routine, $(ALL), $(LIB)($(routine).o))
OBJS     = $(foreach routine, $(ALL), $(routine).o)
TARFILE  = $(patsubst %.a, %.tar, $(LIB))

# -----------------------------------------------------------

$(LIB): $(LMOD) 
	$(RANLIB) $(LIB)

install: $(LIB)
	-cp $(LIB) $(LIBDIR)/$(LIB)
	if [ ! -d $(INCDIR)/speclib ]; \
	then mkdir $(INCDIR)/speclib; \
	else rm -f $(INCDIR)/speclib/* ; \
	fi
	-cp $(HEADER) $(INCDIR)/speclib/. 

uninstall:
	-rm -f  $(LIBDIR)/$(LIB)
	-rm -rf $(INCDIR)/speclib

tar: 
	-rm -f $(TARFILE)
	tar cvf $(TARFILE) README Makefile *.c *.y *.h

tidy:
	-rm -f *~ \#*\#
clean: tidy
	-rm -f *.o
empty: clean
	-rm -f $(LIB) $(TARFILE)

# ------------------------------------------------------------------
#                 Special rules for certain files
# ------------------------------------------------------------------

$(OBJS): config.h

