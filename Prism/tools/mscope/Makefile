#
# Makefile for Mscope
#

include $(PRISM)/system/defaults

MODULES  = \
	mscope \
	display \
	actions \
	user \
	domain \
	autorefine

CPPFLAGS += -I./src

OBJS   = $(foreach module, $(MODULES), $(module).o)
SRC    = $(foreach module, $(MODULES), $(module).c)
VPATH  = src
LIBS   = -L$(PRISM)/lib $(CUBIT) $(VDB) $(VECLIB) $(LAPACK) $(BLAS) \
	$(COMM) $(PL) $(SYSTEM)

# This fixes a problem on IRIX64 machines when mscope is linked against
# the GNU libplot libraries.

ifeq ($(ARCH),IRIX64)
CFLAGS  += -multigot
LDFLAGS += -multigot
endif

# ----------------------------------------------------------------------

mscope: $(OBJS)
	$(LD) $(LDFLAGS) -o mscope $(OBJS) $(LIBS)

install: mscope
	-cp mscope $(PRISM)/bin/mscope
uninstall:
	$(RM) $(PRISM)/bin/mscope
tidy:
	$(RM) *~ \#*\#
clean: tidy
	$(RM) *.o core
empty: clean
	$(RM) mscope

