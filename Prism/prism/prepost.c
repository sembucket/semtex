/* 
 * PREPOST: data processing functions
 *
 * $Id$
 *
 * Author: R. D. Henderson
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "prism/prism.h"
#include "prism/config.h"
#include "veclib/veclib.h"

static char *prog  = "prism";
static char *usage = "[options] session[.rea]";

BSystem *BuildMatrix (Field *, Bedge *, char *, double);

static struct {              /* Default options for Prism         */
  char *name;
  int   value;
  char *descrip;
} prism_ops[] = {
  "checkpt", 0, "-chk      ... checkpoint field dumps",
  "direct" , 2, "-i#       ... set the Helmholtz solver mode (# = [0-2])",
  "logfile", 0, "-log file ... redirect log output to a named file",
  "norder" , 0, "-n #      ... run at a specified polynomial order",
  "core"   , 1, "-swap     ... swap matrix files from disk",
  "verbose", 1, "-v#       ... set the \"verbosity\" level (# = [0-2])",
  "prep"   , 0, "-P        ... run in preprocessor mode, implies -k",
  "matdir" , 0, "-D dir    ... use an alternate directory for matrix files",
   NULL, 0, NULL
};

static struct {
  char *name;
  int   value;
} prism_iparam[] = {
  "Rotational",      Rotational,
  "SkewSymmetric",   SkewSymmetric,
  "Stokes",          Stokes,
  "Convective",      Convective,
  "DIM",             2,
  "TORDER",          2,
   NULL,             0
};

static struct {
  char  *name;
  double value;
} prism_dparam[] = {
  "Re",         1.,       /* Reynolds number           */
  "DENSITY",    1.,       /* Fluid density             */
  "VELOCITY",   1.,       /* Base velocity scale       */
  "LENGTH",     1.,       /* Base length scale         */
  "TIME",       1.,       /* Time interval             */
  "TIME_0",     0.,       /* Initial time              */
  "TIME_N",     0.,       /* Final time                */
  "DT",         0.001,    /* Time step                 */
  "FFZ",        0.,       /* Forcing                   */
  "FFY",        0.,       /*                           */
  "FFX",        0.,       /*                           */
  "FLOWRATE",   0.,       /* Flowrate                  */
  "BETA",       1.,       /* Spanwise wavenumber       */
   NULL,        0.
};

/*
 * Prism's initialization and shutdown routines
 */

void Prism_init (int *argc, char **argv[])
{
  int n;

  /* initialize libraries */

  speclib_init();
#ifdef PARALLEL
  comm_init (argc,argv);
  option_set("procid", comm_rank());
  option_set("nprocs", comm_size());
#endif

  /* set default values for options and parameters */

  for (n = 0; prism_ops[n].name; n++) 
    option_set (prism_ops[n].name, prism_ops[n].value);
  for (n = 0; prism_iparam[n].name; n++)
    iparam_set (prism_iparam[n].name, prism_iparam[n].value);
  for (n = 0; prism_dparam[n].name; n++)
    dparam_set (prism_dparam[n].name, prism_dparam[n].value);
}

/* Calling this function shuts down execution */

void Prism_exit(void) {
  speclib_exit();
#ifdef PARALLEL
  comm_exit();
#endif
  exit(0);
}

/* ------------------------------------------------------------------------- *
 * Prism's error handling and message logs.  All output should be directed   *
 * through these functions.                                                  *
 *                                                                           *
 * Do not print directly to stdout or stderr from within an application!     *
 *                                                                           *
 * PARALLEL: Log output is only allowed from processor 0.  All other         *
 *           processors must generate output as an error message.            *
 * ------------------------------------------------------------------------- */

static FILE *logfp = stdout;

void Prism_logfile (const char *log) 
{
  ROOTONLY {
    FILE *fp = fopen(log,"w");
    assert(fp);
    Prism_logfp(fp);
  }
}

void Prism_logfp (FILE *fp) {
  if (logfp != stdout) 
    fclose (logfp);
  logfp = fp;
}

void Prism_log (int level, const char *fmt, ...) 
{
  va_list ap;
  va_start(ap,fmt);
  ROOTONLY {
    if (level <= option("verbose"))
      vfprintf (logfp, fmt, ap);
  }
  va_end(ap);
}

void Prism_error (const char *fmt, ...) 
{
  va_list ap;

  va_start(ap,fmt);
  fprintf (stderr, "error: ");
  vfprintf(stderr, fmt, ap);
  va_end  (ap);

  Prism_exit();
}

int Prism_options() {
  speclib_options(logfp);
  return 0;
}

int Prism_params() {
  speclib_params(logfp);
  return 0;
}

/*
 * Parse command-line arguments
 */

void parse_args(int argc, char *argv[])
{
  int n;

  if (argc == 1) {
    ROOTONLY
      fprintf (stderr, "usage: %s %s\n", prog, usage);
    Prism_exit();
  }

  for (n = 1; n < argc; n++) {
    if (*argv[n] == '-') {
      char c = argv[n][1];
      switch (c) {
      case 'h':
	ROOTONLY {
	  fprintf (stderr, "usage: %s %s\n\n", prog, usage);
	  fprintf (stderr, "options:\n");
	  for (n = 0; prism_ops[n].name; n++)
	    fprintf (stderr, "%s\n", prism_ops[n].descrip);
	}
	Prism_exit();
	break;
      case 'n':
	option_set("norder", atoi(argv[++n]));
	break;
      case 'c': /* -chk */
	option_set("checkpt", 1);
	break;
      case 's': /* -swap */
	option_set("core", 0);
	break;
      case 'm': /* undocumented */
	option_set("mpack", 0);
	break;
      case 'i':
	if (isdigit(argv[n][2]))
	  option_set("direct", 2-atoi(&argv[n][2]));
	else
	  option_set("direct", option("direct")-1);
	break;
      case 'l':
	Prism_logfile(argv[++n]);
	break;
      case 'v':
	if (isdigit(argv[n][2]))
	  option_set("verbose", atoi(&argv[n][2]));
	else
	  option_set("verbose", option("verbose")+1);
	break;
      case 'P':
	option_set("prep", 1);
	break;
      case 'D':
	fprintf (stderr, "%s: option -D not implemented yet\n", prog);
	break;
      default:
	fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
	break;
      }
    }
  }


  /* check the last argument for a ".rea" suffix */

  if (strstr(argv[n=argc-1],".rea")) {
    char *p = strstr(argv[n],".rea");
    *p = '\0';
  }

  /* Now everything is initialize and we're ready to start.  *
   * Put the opening banner and version info onto the log.   */

  Prism_log(0,"%s%dd Version %s [%s]\n"
	      "Copyright (c) 1994-1998 R. D. Henderson\n\n",
	       prog, DIM, VERSION, ARCH);

#ifdef PARALLEL
  Prism_log(0,"Running on %d processors\n\n", option("nprocs"));
#endif  
}

/* ----------------------------------------------------------------------- *
 * Summary() - Print a summary of the input data                           *
 *                                                                         *
 * This prints a summary of the symbols that have been stored in the       *
 * symbol table manager.                                                   *
 * ----------------------------------------------------------------------- */

void Summary (void)
{
  char buf[BUFSIZ];

  Prism_log(0, "Reynolds number    : %g\n", dparam("Re"));
  Prism_log(0, "Elements           : %d\n", iparam("ELEMENTS"));
  Prism_log(0, "Polynomial order   : %d\n", iparam("NORDER")-1);

#if DIM == 3
  Prism_log(0, "Number of z-planes : %d\n", iparam("NZ"));
  Prism_log(0, "Length in z        : %g, Beta = %g\n", 
	    dparam("LZ"), dparam("BETA"));
#endif

  Prism_log(0, "Integration order  : %d\n", iparam("TORDER"));
  Prism_log(0, "Time step          : %g\n", dparam("DT"));
  Prism_log(0, "Integration time   : %g, or %d steps\n", 
	    dparam("TIME"), iparam("NSTEPS"));

  Prism_log(0, "I/O time for saves : %g, or %d steps",
	    dparam("IO_FLD_DT"), iparam("IO_FLD"));
  if (option("checkpt"))
    Prism_log(0, " [checkpoint]");
  Prism_log(0, "\n");

  Prism_log(0, "Equation type      : ");
  switch (iparam("EQTYPE")) {
  case SkewSymmetric:
    Prism_log(0, "Navier-Stokes (Skew-Symmetric)\n");
    break;
  case Rotational:
    Prism_log(0, "Navier-Stokes (Rotational)\n");
    break;
  case Stokes:
    Prism_log(0, "Stokes\n");
    break;
  default:
    Prism_log(0, "User-defined\n");
    break;
  }

  switch (option("direct")) {
  case 0:
    Prism_log(0, "Pressure tolerance : %.1e\n", dparam("TOLP"));
    /* drop through */
  case 1:
    Prism_log(0, "Velocity tolerance : %.1e\n", dparam("TOLU"));
    break;
  default:
    break;
  }
  
  if (dparam ("XSCALE") != 1. || dparam("YSCALE") != 1. ||
      dparam ("XSHIFT") != 0. || dparam("YSHIFT") != 0. )  
    {
      Prism_log(0, "Mesh adjustments   : x(%g,%g), +(%g,%g)\n",
		dparam ("XSCALE"), dparam ("YSCALE"),
		dparam ("XSHIFT"), dparam ("YSHIFT"));
    }

  if (iparam("NSTEPS") / iparam("IO_FLD") > 10) {
    Prism_log(2,"\nATTENTION: You have more than 10 dumps..."
	      "did you want to use -chk?\n");
  }
}

/*
 *  Post-processing
 */

void PostProcess(Domain *omega)
{
  int nsteps = iparam("NSTEPS");
  int iostep = iparam("IO_FLD");

  /* Check to see if solution should be dumped */
  
  if (nsteps % iostep)
    Domain_save(omega);
}

/* ----------------------------------------------------------------------- *
 * ReadICs() - Read Initial Conditions                                     *
 *                                                                         *
 * This function reads in the initial conditions for a set of fields, as   *
 * determined by one of the following keywords in the line following the   *
 * begining of SEC_ICS:                                                    *
 *                                                                         *
 *     Default            Every field has the value 0                      *
 *     Given              Function specification for each field            *
 *     Restart            Match field values from the data file specified  *
 *                        on the next line.                                *
 *                                                                         *
 * Note: For a restart, the type of each field is matched to the types     *
 *       stored in the file.  If any field is missing, ReadICs generates   *
 *       an error and exits.                                               *
 * ----------------------------------------------------------------------- */

static void Restart (char *name, int nfields, Element *U[]);
static char *SEC_ICS = "INITIAL CONDITIONS";

void ReadICs (FILE *fp, int nfields, ...)
{
  int     n, i, type;
  char    buf[BUFSIZ];
  Field  *U [_MAX_FIELDS];
  va_list ap;

  static  char *ictypes[] = { "PrepOnly", "Restart", "Default", "Given" };

  /* Initialize the Field pointers */

  if (nfields > _MAX_FIELDS)
    Prism_error("Prism: too many fields in ReadICs\n");

  va_start(ap, nfields);
  for (i = 0; i < nfields; i++) U[i] = va_arg(ap, Field *);
  va_end  (ap);

  if (!findSection (SEC_ICS, buf, fp))
    Prism_error("ReadICs -- initial conditions not specified:\n%s", buf);

  if (sscanf(buf, "%d", &n) != 1) {
    fgets (buf, BUFSIZ, fp);
    if (sscanf(buf, "%d", &n) != 1) {
      Prism_error("ReadICs: cannot determine number of initial conditions");
    }
  }

  if (option("prep"))              /* Preprocessor */
    type = 0;
  else if (n == 0)                 /* Default */
    type = 2;
  else {                           /* Find it... */
    n--; fgets (buf, BUFSIZ, fp);
    for (i = 1; i < 4; i++) 
      if (strstr(buf, ictypes[type = i])) break;    
  }
  
  switch (type) {
  case 0: break;                      /* Pre-processor only */
  case 1: {                           /* Restart from a file */
    fgets   (buf, BUFSIZ, fp); n--;
    Restart (buf, nfields, U);
    break;
  }
  case 2: {                           /* Default Initial Conditions */
    const int
    ntotz  = U[0]->nr * U[0]->ns * U[0]->nz * Field_count(U[0]);
    for (i = 0; i < nfields; i++) 
      dzero (ntotz, *U[i]->base, 1);

    Prism_log(0,"Initial condition  : Default\n");
    break;
  }
  case 3: {                          /* Function(x,y,z) Initial conditions */
    double  *z, *field;
    char    *s;
    register int nz, ntot, k, pid;

    ntot = U[0]->nr * U[0]->ns * Field_count(*U);
    nz   = U[0]->nz;
    pid  = option("procid");
    z    = zmesh(iparam("NZ")) + pid*nz;

    Prism_log(0,"Initial condition  : Given\n");
      
    for (i = 0; i < nfields && n--; i++) {

      if(!(s = strchr(fgets(buf,BUFSIZ,fp),'=')))
	Prism_error("ReadICs -- the following function is invalid:\n%s", s);

      while (isspace(*++s));
      field = *U[i]->base;
      vector_def("x y", s);
      for(k = 0; k < nz; k++, field += ntot) {
	scalar_set ("z", z[k]);
	vector_set (ntot, *U[i]->xmesh, *U[i]->ymesh, field);
      }

      Prism_log(0,"\t%c = %s", FIELD_TYPE(U[i]), s);
    }
    free (z - pid*nz);
    break;
  }
  default:
    Prism_error("Prism: invalid initial conditions\n");
    break;
  }
  
  while (n--) fgets(buf, BUFSIZ, fp);    /* read to the end of the section */

  return;
}

/* Initialize variable from a data file */

static void Restart (char *name, int nfields, Field *U[])
{
  FieldFile *f;
  FILE *fp;
  char *p;
  char fname[FILENAME_MAX];
  char buf[BUFSIZ];
  int  nr, ns, nz;
  int  i, dump;

  if (sscanf(name, "%s", fname) != 1)
    Prism_error("ReadICs: couldn't read the restart file name\n");
  if ((fp = fopen(fname,"r")) == (FILE *) NULL)
    Prism_error("ReadICs -- restart file not found: %s", fname);
  
  dump  = 0;
  nr    = U[0]->nr;
  ns    = U[0]->ns;
  nz    = U[0]->nz;
  f     = FieldFile_alloc();

  /* Read to the end of the file */

  while (!feof(fp) && FieldFile_read(f,fp) >= 0) 
    dump++;

  if (!dump) 
    Prism_error("Prism: no dumps read from restart file\n");
  
  /* Project to the current resolution */
  
  Prism_log(0,"Initial condition  : %s at t = %g", fname, f->time);

  if (nz != f->nz) {
    Prism_log(0,", z:[%d -> %d]", FIELDFILE_NZ(f), nz);
    FieldFile_projectz(f, nz);
  }

  if (f->nr != nr || f->ns != ns) {
    Prism_log(0,", xy:[%d x %d -> %d x %d]", f->nr, f->ns, nr, ns);
    FieldFile_project (f, nr, ns);
  }

  Prism_log(0,"\n");

  /* Load each solution from the FieldFile */

  for (i = 0; i < nfields; i++)
    FieldFile_get (f, U[i]);
  
  dparam_set("TIME_0", f->time);
  dparam_set("TIME_N", dparam("TIME") + f->time);
  FieldFile_free (f);
  fclose (fp);
}

/* ----------------------------------------------------------------------- *
 * ReadDF() - Read Drive Force data                                        *
 *                                                                         *
 * This is really one function designed to handle two different types of   *
 * drive force specification.  If the "steady" flag is active, the drive   *
 * force can be an arbitary function.  Otherwise, it must be constant but  *
 * can have any direction.                                                 *
 *                                                                         *
 * If the flowrate option is active, this section is skipped entirely.     *
 *                                                                         * 
 * Example:                                                                *
 *                                                                         *
 *   ***** DRIVE FORCE DATA ***** PRESSURE GRAD, FLOW, Q                   *
 * 4                   Lines of Drive force data follow                    *
 *      FFX = 0.                                                           *
 *      FFY = 0.                                                           *
 *      FFZ = 2. * KINVIS                                                  *
 * C                                                                       *
 *                                                                         *
 * The names given to these lines does not matter, but their order is      *
 * taken as (x,y,z).  The '=' MUST be present, and everything to the right *
 * of it determines the forcing functions.                                 *
 * ----------------------------------------------------------------------- */

#define eq_Helmholtz  0         /* Kludge */
#define eq_Poisson    1
#define eq_Laplace    2

void ReadDF (FILE *fp, int nforces, ...)
{
  int  nlines, i;
  char buf[BUFSIZ], *p;
  static char *ftypes[] = { "FFX", "FFY", "FFZ" };
  static char *SEC_DF   = "DRIVE FORCE";

  if (!findSection (SEC_DF, buf, fp))
    return;
  if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &nlines) != 1)
    Prism_error("Prism: can't read the # of lines in DRIVE FORCE\n");

  /* Don't process this section if running in pre-processor mode  */
  /* Don't process this section if running in parallel & not root */

  if (option("prep") || option("procid")) {
    while (nlines--) fgets(buf, BUFSIZ, fp);
    return;
  }

  /* Check to see if the problem is steady (Helmholtz) or unsteady */

  if (option("steady")) {
    int type;
    va_list ap;
    Field *U[_MAX_FIELDS];
    
    va_start (ap, nforces);
    for (i = 0; i < nforces; i++) U[i] = va_arg(ap, Field *);
    va_end (ap);


    switch (nlines ? iparam("EQTYPE") : eq_Laplace) {
      case eq_Poisson: 
      case eq_Helmholtz:
	Prism_log(0,"Forcing function   : Given\n");
	for (i = 0; i < nforces && nlines--; i++) {
	  if (p = strchr(fgets(buf, BUFSIZ, fp), '=')) {
	    while (isspace(*++p));
	    Prism_log(0,"\t%c = %s\n", FIELD_TYPE(U[i]), p);

	    vector_def ("x y", p);
	    vector_set (U[i]->nr * U[i]->ns * Field_count(U[i]), 
		       *U[i]->xmesh, *U[i]->ymesh, *U[i]->base);
	  } else {
	    Prism_error("ReadDF -- invalid function:\n%s", buf);
	  }
	}
	break;
      case eq_Laplace:
	for  (i = 0; i < nforces; i++)
	  dzero (U[i]->nr * U[i]->ns * Field_count(U[i]), *U[i]->base, 1);
	break;
      }

  } else {    /* ..... Force specification for unsteady problems ..... */

    if (option("procid")) {              /* check the processor number */
      dparam_set ("FFX", 0.);
      dparam_set ("FFY", 0.);
      dparam_set ("FFZ", 0.);
    }
    else if (dparam("FLOWRATE") != 0.)           /* check for flowrate */
      option_set ("flowrate", 1);    
    else if (nlines) {                             /* read the list... */
      for (i = 0; i < nforces && nlines--; i++) {
	for(p = fgets(buf, BUFSIZ, fp); isspace(*p); p++);
	if (toupper(*p) == 'C') {
	  /* do nothing */
	} else {
	  if (p = strchr(buf,'='))
	    dparam_set(ftypes[i], scalar(++p));
	  else
	    Prism_error("ReadDF -- invalid function:\n%s", buf);
	}
      }
    }

    if (option("flowrate")) {
      Prism_log(0,"Flowrate           : %g\n", dparam("FLOWRATE"));
    } else {
      Prism_log(0,"Drive Force        : [%g, %g, %g]\n", 
		dparam("FFX"), dparam("FFY"), dparam("FFZ"));
    }
  }

  while (nlines--) fgets(buf, BUFSIZ, fp);    /* finish the section */
}

/* ----------------------------------------------------------------------- *
 * ReadHisData() - Read History and Integral data                          *
 *                                                                         *
 * This function reads in and sets up structures for the processing of     *
 * history points.  Each "history point" defines a location in the mesh    *
 * where field values will be sampled and written to the history file.     *
 * Each history point's specification has the format:                      *
 *                                                                         *
 *                      fields  type  location  
 *                                                                         *
 * Each token is explained below.                                          *
 *                                                                         *
 * fields       String of characters specifying the variable types to      *
 *              sample, general a subset of "UVWP" (case doesn't matter)   *
 *                                                                         *
 * type         H = standard history point, P = probe (see below)          *
 *                                                                         *
 *   location = i j k m (History Point)                                    *
 *              For a standard history point, this gives the (ijkm)-mesh   *
 *              index of the node:  (ij) are local node indices, (k) is    *
 *              the element number, and (m) is the frame number.           *
 *                                                                         *
 *   location = x y m   (Probe)                                            *
 *              For a probe, you just specify the (x,y)-coordinates and    *
 *              the frame number.  Data is generated by evaluating the     *
 *              local polynomial basis.                                    *
 *                                                                         *
 * Notes                                                                   *
 * -----                                                                   *
 * frame        Frame numbers can refer to either a physical z-plane       *
 *              (0 <= m <  NZ) or a Fourier mode (0 <= m < NZ/2).   To     *
 *              sample a Fourier mode, include the letter 'r' or 'i' after *
 *              the mode number to get the real or imaginary component.    *
 *                                                                         *
 * Examples                                                                *
 * --------                                                                *
 *     UVP   H   4   8    37  1     # node (4,8), frame 1, in element 37   *
 *     UV    P   1.  1.   1         # frame 1 at location x=1.,y=1.        *
 *     UVW   P   5.  0.   8r        # real part of mode 8 at x=5.,y=0.     *
 * ----------------------------------------------------------------------- */

static HisPoint *appendHisPoint (HisPoint *list, HisPoint *hp);

void ReadHisData (FILE *fp, Domain *omega)
{
  char      buf[BUFSIZ], *p, *mspec;
  int       npts, norder, nz, nel, pid;
  HisPoint *hp, *his_list = (HisPoint*)NULL;
  register int n;

  static char *SEC_HIST = "HISTORY";

  if (!findSection (SEC_HIST, buf, fp))
    return;
  if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &npts) != 1)
    Prism_error("Prism: can't read the number of HISTORY points\n");

  GSYNC;

  pid    = option("procid");
  nz     = iparam("NZ");
  norder = iparam("NORDER");
  nel    = iparam("ELEMENTS");

  for (n = 1; n <= npts; n++) {
    hp     = (HisPoint*) calloc(1, sizeof(HisPoint));
    hp->id = n;
    fgets(buf, BUFSIZ, fp);

    /* Break the line into tokens separated by white space */
    
    p = strtok(buf, "\t ");
    while (*p = tolower(*p))
      strncat (hp->fields, p++, 1);
    
    /* Read the location info */
    
    switch (*(p = strtok (NULL, "\t "))) {
    case 'H': {
      int i = atoi(strtok(NULL,"\t "));
      int j = atoi(strtok(NULL,"\t "));
      int k = atoi(strtok(NULL,"\t "));
      
      const double x = omega->U[k].xmesh[j][i];
      const double y = omega->U[k].ymesh[j][i];

      if (i-- > norder || j-- > norder || k-- > nel)
	Prism_error("bad history point");

      hp->locator = Probe_alloc(omega->U, PROBE_INDEX, x, y);
      break;
    }
    case 'P': {
      const double x = atof(strtok(NULL,"\t "));
      const double y = atof(strtok(NULL,"\t "));

      hp->locator = Probe_alloc(omega->U, PROBE_XP, x, y);
      break;
    }
    default:
      Prism_error("bad history point");
      break;
    }
    
# if DIM == 2
    hp->frame = 0;
    hp->mode  = Fourier;
# else
    hp->frame = strtol(strtok(NULL,"\t "), &p, 10);
    
    if (isalpha (*p)) {
      if (hp->frame < 0 || hp->frame > nz/2-1 || !strchr ("ri", *p))
	Prism_error("bad history point");
      hp->mode  = Fourier;
      hp->frame = 2*hp->frame + (*p == 'r' ? 0 : 1);
    } else {
      if (hp->frame < 0 || hp->frame >= nz)
	Prism_error("bad history point");
      hp->mode  = Physical;
    }
# endif
    
    if (!hp->locator)
      Prism_error("bad history point");
    else
      his_list = appendHisPoint (his_list, hp);
  }

  /* Check the output frequency */
  
  if (omega->his_list = his_list) {
    if (! iparam("IO_HIS"))
      iparam_set("IO_HIS", 1);
    
    /* Echo the history points */

    Prism_log(0,"\nHistory points:\n");
    for (hp = his_list; hp; hp = hp->next) {
      Prism_log(0,"%2d: %6s @ %#7.4lf %#7.4lf ", 
		hp->id, hp->fields, hp->locator->x, hp->locator->y);
# if DIM == 3
#	
      if (hp->mode == Physical)
	Prism_log(0,"frame %2d [phys]", hp->frame);
      else
	Prism_log(0,"mode  %2d [%s]", hp->frame >> 1, 
		  (hp->frame & 1 ? "imag" : "real"));
#
# endif
      Prism_log(0,"\n");
    }
  }
}

static HisPoint *appendHisPoint (HisPoint *list, HisPoint *hp)
{
  HisPoint *h = list;

  if (h) {
    while (h->next) h = h->next;
    h->next = hp;
  } else
    list = hp;

  return list;
}

/* ---------------------------------------------------------------------- *
 * File_backup() -- Create a backup copy of a file                        *
 * ---------------------------------------------------------------------- */

extern int unlink (const char *path);
extern int link   (const char *path1, const char *path2);

int File_backup (char *path1)
{
  char path2[FILENAME_MAX];
  int  stat;

  sprintf (path2, "%s.bak", path1);
  unlink  (path2);                    /* unlink path2 regardless    */
  if (!(stat = link(path1, path2)))   /* try to link path1 -> path2 */
    unlink (path1);                   /* unlink path1 only if the   */
  return stat;                        /* link was sucessful         */
}
