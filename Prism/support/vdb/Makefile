#
# Makefile for the VDB library
#

include $(PRISM)/system/defaults

LIB      = libvdb.a
LIBDIR   = $(PRISM)/lib
INCDIR   = $(PRISM)/include
CPPFLAGS = -I$(PRISM)/support

# Object modules

MODS     = vdb sync combine voxel alias
HEADER   = vdb.h 

LMODS = $(foreach routine, $(MODS), $(LIB)($(routine).o))
OBJS  = $(foreach routine, $(MODS), $(routine).o)
SRC   = $(foreach routine, $(MODS), $(routine).c)

# ------------------------------------------------------------------

$(LIB): $(LMODS)

install: $(LIB)
	-cp $(LIB) $(LIBDIR)/$(LIB) 
	if [ ! -d $(INCDIR)/vdb ]; \
	then mkdir $(INCDIR)/vdb; \
	else $(RM) $(INCDIR)/vdb/*; \
	fi
	cp $(HEADER) $(INCDIR)/vdb/.

uninstall:
	$(RM)    $(LIBDIR)/$(LIB)
	$(RM) -r $(INCDIR)/vdb

tidy:
	$(RM) *~ \#*\#
clean: tidy
	$(RM) *.o core
empty: clean
	$(RM) $(LIB)

# ------------------------------------------------------------------

vdb.o: vdb.h voxel.h 
sync.o: vdb.h voxel.h
combine.o: vdb.h voxel.h
voxel.o: voxel.h
