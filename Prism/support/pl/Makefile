#
# Makefile for libpl.a
#
# $Id$

include $(PRISM)/system/defaults

LIB      = libpl.a
LIBDIR   = $(PRISM)/lib
INCDIR   = $(PRISM)/include
CPPFLAGS = -I$(PRISM)/support $(PL_INC)

# Module names and headers

ALL      = pl driver
HEADER   = pl.h devices.h ptype.h ltype.h config.h

# Generate rules for building the modules

LMOD     = $(foreach routine, $(ALL), $(LIB)($(routine).o))
OBJS     = $(foreach routine, $(ALL), $(routine).o)
TARFILE  = $(patsubst %.a, %.tar, $(LIB))

# ------------------------------------------------------------------

$(LIB): $(LMOD) 
	$(RANLIB) $(LIB)

install: $(LIB)
	-mv $(LIB) $(LIBDIR)/$(LIB)
	if [ ! -d $(INCDIR)/pl ]; \
	then mkdir $(INCDIR)/pl; \
	else rm -f $(INCDIR)/pl/*; \
	fi
	-cp $(HEADER) $(INCDIR)/pl/.

uninstall:
	-rm -f  $(LIBDIR)/$(LIB)
	-rm -rf $(INCDIR)/pl

tar: 
	-rm -f $(TARFILE)
	tar cvf $(TARFILE) README Makefile *.c *.h

tidy:
	-rm -f *~ \#*\#
clean: tidy
	-rm -f *.o
empty: clean
	-rm -f $(LIB) $(TARFILE)

# ------------------------------------------------------------------

$(OBJS): pl.h config.h
