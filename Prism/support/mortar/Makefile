#
# Special library functions for mortar-based patching
#

include $(PRISM)/system/defaults

LIB      = libmortar.a
LIBDIR   = $(PRISM)/lib
INCDIR   = $(PRISM)/include
CPPFLAGS = -I$(PRISM)/support

# Module names and headers

ALL     = mortar matrix_q solve_q solve_cg_q residual_q
HEADER  = mortar.h

# Generate rules for building the modules

LMOD    = $(foreach routine, $(ALL), $(LIB)($(routine).o))
OBJS    = $(foreach routine, $(ALL), $(routine).o)
TARFILE = $(patsubst %.a, %.tar, $(LIB))

# ------------------------------------------------------------------

$(LIB): $(LMOD) 
	$(RANLIB) $(LIB)

install: $(LIB)
	-cp $(LIB) $(LIBDIR)/$(LIB)
	if [ ! -d $(INCDIR)/mortar ]; \
	then mkdir $(INCDIR)/mortar; \
	else $(RM) $(INCDIR)/mortar/*; \
	fi
	-cp $(HEADER) $(INCDIR)/mortar/. 

uninstall:
	$(RM)    $(LIBDIR)/$(LIB)
	$(RM) -r $(INCDIR)/mortar

tar:
	$(RM) $(TARFILE)
	tar cvf $(TARFILE) README Makefile *.c include/*.h

tidy:
	$(RM) *~ \#*\#
clean: tidy
	$(RM) *.o
empty: clean
	$(RM) $(LIB) $(TARFILE)

# ------------------------------------------------------------------

$(OBJS): $(PRISM)/prism/config.h
