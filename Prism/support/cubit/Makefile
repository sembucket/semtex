#
# Compile this crap into a library...
#

include $(PRISM)/system/defaults

LIBDIR   = $(PRISM)/lib
INCDIR   = $(PRISM)/include
CPPFLAGS = -I$(PRISM)/support

ALL = cubit utilities manager tree splib session \
	mesh meshio isomesh curve edge vertex element family alias connect bc \
	field fieldfile probe surface \
	error \
	matrix sparse cg residual operators \
	param keyword 

HEADERS = cubit.h \
	alias.h bc.h cubitypes.h curve.h edge.h element.h error.h family.h \
	field.h fieldfile.h isomesh.h keyword.h manager.h matrix.h mesh.h \
	param.h probe.h session.h splib.h tree.h vecerr.h vertex.h

OBJS = $(foreach routine, $(ALL), $(routine).o)
SRC  = $(foreach routine, $(ALL), $(routine).c)
LMOD = $(foreach routine, $(ALL), $(LIB)($(routine).o))
LIB  = libcubit.a

# -------------------------------------------------------------------------

$(LIB): $(LMOD) 

$(OBJS): $(HEADERS)

install: $(LIB)
	-cp $(LIB) $(LIBDIR)/$(LIB)
	if [ ! -d $(INCDIR)/cubit ]; \
	then mkdir $(INCDIR)/cubit; \
	else $(RM) $(INCDIR)/cubit/*; \
	fi
	-cp $(HEADERS) $(INCDIR)/cubit/.

uninstall:
	$(RM)    $(LIBDIR)/$(LIB)
	$(RM) -r $(INCDIR)/cubit
tidy:
	$(RM) *~ \#*\#
clean: tidy
	$(RM) *.o
empty: clean
	$(RM) *.a

