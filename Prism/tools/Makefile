#
# Tools
#
# $Id$
# ----------------------------------------------------------------------------

include $(PRISM)/system/defaults

TOOLS    = \
	curl \
	ff \
	fminmax \
	fnrm2 \
	meshgen \
	perturb \
	profile \
        sem2fast \
        sem2tec \
        xfreq \
	xplane \
	xsep \
	forces 

SCRIPTS    = save renumber
PY_SCRIPTS = save.py

CPPFLAGS = -I$(PRISM)/include
LIBS     = -L$(PRISM)/lib $(MORTAR) $(SPECLIB) $(VECLIB) $(LAPACK) $(BLAS) \
		$(COMM) $(SYSTEM)

all: $(TOOLS)

profile: profile.o
	$(LD) $(LDFLAGS) -o profile profile.o $(LIBS)
fnrm2: fnrm2.o
	$(LD) $(LDFLAGS) -o fnrm2 fnrm2.o $(LIBS)
fminmax: fminmax.o
	$(LD) $(LDFLAGS) -o fminmax fminmax.o $(LIBS)
curl: curl.o
	$(LD) $(LDFLAGS) -o curl curl.o $(LIBS)
ff: ff.o
	$(LD) $(LDFLAGS) -o ff ff.o $(LIBS)
sem2fast: sem2fast.o
	$(LD) $(LDFLAGS) -o sem2fast sem2fast.o $(LIBS)
sem2tec: sem2tec.o sem_utils.o mesh.o
	$(LD) $(LDFLAGS) -o sem2tec sem2tec.o sem_utils.o mesh.o $(LIBS)
meshgen: meshgen.o mesh.o
	$(LD) $(LDFLAGS) -o meshgen meshgen.o mesh.o $(LIBS)
xplane: xplane.o
	$(LD) $(LDFLAGS) -o xplane xplane.o $(LIBS)
sinterp: sinterp.o mesh.o sem_utils.o 
	$(LD) $(LDFLAGS) -o sinterp sinterp.o mesh.o sem_utils.o $(LIBS)
perturb: perturb.o
	$(LD) $(LDFLAGS) -o perturb perturb.o $(LIBS)
forces: forces.o
	$(LD) $(LDFLAGS) -o forces forces.o $(LIBS)
xsep: xsep.o
	$(LD) $(LDFLAGS) -o xsep xsep.o $(LIBS)



# ----------------------------------------------------------------------------

install: $(TOOLS)
	for f in $(TOOLS); do \
		$(RM) $(PRISM)/bin/$$f; \
		cp $$f $(PRISM)/bin/$$f; \
	done
	for f in $(SCRIPTS); do \
		$(RM) $(PRISM)/bin/$$f; \
		cp etc/$$f $(PRISM)/bin/$$f; \
	done
	for f in $(PY_SCRIPTS); do \
		$(RM) $(PRISM)/bin/$$f; \
		cp Python/$$f $(PRISM)/bin/$$f; \
	done
	( cd mason  ; $(MAKE) install )
	( cd mscope ; $(MAKE) install )
uninstall:
	for f in $(TOOLS) $(SCRIPTS); do \
		$(RM) $(PRISM)/bin/$$f; \
	done	
	( cd mason  ; $(MAKE) uninstall )
	( cd mscope ; $(MAKE) uninstall )
tidy:
	$(RM) *~ \#*\#
	( cd mason  ; $(MAKE) tidy )
	( cd mscope ; $(MAKE) tidy )
clean: tidy
	$(RM) *.o
	( cd mason  ; $(MAKE) clean )
	( cd mscope ; $(MAKE) clean )
empty: clean
	$(RM) $(TOOLS)
	( cd mason  ; $(MAKE) empty )
	( cd mscope ; $(MAKE) empty )

# ----------------------------------------------------------------------------

mesh.o: mesh.c mesh.h
sem_utils.o: sem_utils.c sem_utils.h
