#
# A portable communications library
#

include $(PRISM)/system/defaults

LIB      = libcomm.a
LIBDIR   = $(PRISM)/lib
INCDIR   = $(PRISM)/include
CPPFLAGS = -I$(PRISM)/support

# Module names and headers

ALL      = comm
HEADER   = comm.h config.h

# Generate rules for building the modules

LMOD     = $(foreach routine, $(ALL), $(LIB)($(routine).o))
OBJS     = $(foreach routine, $(ALL), $(routine).o)
TARFILE  = $(patsubst %.a, %.tar, $(LIB))

# ------------------------------------------------------------------

$(LIB): $(LMOD) 
	$(RANLIB) $(LIB)

install: $(LIB)
	-mv $(LIB) $(LIBDIR)/$(LIB)
	if [ ! -d $(INCDIR)/comm ]; \
	then mkdir $(INCDIR)/comm; \
	else rm -f $(INCDIR)/comm/*; \
	fi
	cp $(HEADER) $(INCDIR)/comm/.

uninstall:
	-rm -f  $(LIBDIR)/$(LIB)
	-rm -rf $(INCDIR)/comm 

tar: 
	-rm -f $(TARFILE)
	tar cvf $(TARFILE) README Makefile *.c *.y include/*.h

tidy:
	-rm -f *~ \#*\#
clean: tidy
	-rm -f *.o
empty: clean
	-rm -f $(LIB) $(TARFILE)
