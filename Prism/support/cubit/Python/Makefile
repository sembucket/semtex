#
# Build the python module for cubit
#

include $(PRISM)/system/defaults

CPPFLAGS += $(PY_INC)

OBJS = \
	_cubitmodule.o \
	meshobject.o \
	fieldobject.o 

LIBS = -L$(PRISM)/lib $(CUBIT) $(VDB) $(VECLIB) $(LAPACK) $(BLAS) \
	$(COMM) $(SYSTEM)

_cubitmodule.so: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -shared $(OBJS) $(LIBS) -o $@

tidy:
	$(RM) *~ \#*\#
clean: tidy 
	$(RM) *.o core so_locations
empty: clean
	$(RM) *.so *.pyc

