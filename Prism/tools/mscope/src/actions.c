/*
 * Mscope command parser
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * Copyright (c) 1994-1998 R. D. Henderson and Caltech
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "vdb/vdb.h"
#include "cubit/cubit.h"
#include "cubit/isomesh.h"
#include "cubit/tree.h"
#include "veclib/veclib.h"
#include "pl/pl.h"

#include "mscope.h"
#include "grammar.h"

/* Shared mesh and field data */

extern Domain Geometry;

/* Prototypes */

static Token token_list[] = {
  INPUT,    "input",    "load a command file",         DoInput,
  DOMAIN,   "domain",   "create a new domain",         DoDomain,
  BOX,      "box",      "plot current axes",           DoBox,
  GRID,     "grid",     "plot element boundaries",     DoGrid,
  CONTOUR,  "contour",  "plot field contours",         DoContour,
  STREAM,   "stream",   "plot streamlines",            DoStream,
  VECTORS,  "vectors",  "plot vector fields",          DoVectors,
  PROFILE,  "profile",  "plot field profiles",         DoProfile,
  SPLOT,    "splot",    "plot elemental spectra",      DoSplot,
  GRIDNUM,  "gridnum",  "plot element numbers",        DoGridNum,
  LOGICAL,  "logical",  "plot logical mesh elements",  DoLogical,
  MESH,     "mesh",     "plot the mesh",               DoMesh,
  BCDRAW,   "bc",       "plot boundary conditions",    DoBCs,
  ID,       "id",       "plot element id numbers",     DoID,
  FAMILY,   "family",   "plot element family number",  DoFamily,
  KEY,      "key",      "plot element key values",     DoKey,
  LEVEL,    "level",    "plot element levels",         DoLevel,
  MINFO,    "minfo",    "print mesh info",             DoMinfo,
  EINFO,    "einfo",    "print element info",          DoEinfo,
  BINFO,    "binfo",    "print boundary info",         DoBinfo,
  NORM,     "norm",     "print field norms",           DoNorms,
  ENORM,    "enorm",    "print error norms",           DoEnorms,
  SHOW,     "show",     "(see online help)",           DoShow,
  LIST,     "list",     "(see online help)",           DoList,
  PARAM,    "param",    "define a parameter",          DoParam,
  OPTION,   "option",   "define an option",            DoOption,
  FORCE,    "force",    "define a drive force",        DoForce,
  SOLUTION, "solution", "define an exact solution",    DoSolution,
  DEFINE,   "define",   "define an mscope constant",   DoDefine,
  BCOND,    "bcond",    "define a boundary condition", DoBcond,
  ICOND,    "icond",    "define initial conditions",   DoIcond,
  FIELDS,   "fields",   "define solution vector",      DoFields,
  USER,     "user",     "define a string for export",  DoUser,
  ELEMENT,  "element",  "define an element",           DoElement,
  NODE,     "node",     "define a node",               DoNode,
  CURVE,    "curve",    "define a curved edge",        DoCurve,
  HISTORY,  "history",  "define a history point",      DoHistory,
  MATRIX,   "matrix",   "build the stiffness matrix",  DoMatrix,
  SOLVE,    "solve",    "solve the Helmholtz equation",DoSolve,
  LOAD,     "load",     "see \"meshin\"",              DoLoad,
  SAVE,     "save",     "see \"meshout\"",             DoSave,
  MESHIN,   "meshin",   "read a mesh file",            DoMeshIn,
  FLOWIN,   "flowin",   "read a field file",           DoFlowin,
  MESHOUT,  "meshout",  "write a mesh file",           DoMeshOut,
  FLOWOUT,  "flowout",  "write a field file",          DoFlowout,
  REFINE,   "refine",   "refine elements",             DoRefine,
  SPLIT,    "split",    "split elements",              DoSplit,
  GRADIENT, "gradient", "refine on gradients",         DoGradient,
  SPECTRUM, "spectrum", "refine on spectra",           DoSpectrum,
  REGRESS,  "regress",  "refine on regression",        DoRegress,
  INSTALL,  "install",  "refine on node keys",         DoInstall,
  RESTRICT, "restrict", "apply refinement rules",      DoRestrict,
  PRUNE,    "prune",    "prune a family",              DoPrune,
  CONNECT,  "connect",  "establish connectivity",      DoConnect,
  LOCK,     "lock",     "lock vertex coordinates",     DoLock,
  KEYS,     "keys",     "install keys from a file",    DoKeys,
  DO,       "do",       "do loop",                     DoLoop,

  FSCAL,    "fscal",    "compute v *= d",              DoFscal,
  FABS,     "fabs",     "compute v = abs(u)",          DoFabs,
  FAXPY,    "faxpy",    "compute v += d*u",            DoFaxpy,
  FCOPY,    "fcopy",    "compute v = u",               DoFcopy,
  FADD,     "fadd",     "compute w = u + v",           DoFadd,
  FMULT,    "fmult",    "compute w = u * v",           DoFmult,
  FDIV,     "fdiv",     "compute w = u / v",           DoFdiv,
  FGRAD,    "fgrad",    "compute v = D[u,x_i]",        DoFgrad,
  FFT,      "fft",      "compute v = F[i]",            DoFFT,
  FINT,     "fint",     "display Integrate[u]",        DoFint,

  DEV,      "dev",      "change the display device",     DoDevice,
  LIMITS,   "limits",   "change display limits",         GetLimits,
  ERASE,    "erase",    "erase the display",             DoErase,
  ZOOM,     "zoom",     "zoom in",                       DoZoom,
  UNZOOM,   "unzoom",   "reset the display limits",      DoUnzoom,
  RELOCATE, "relocate", "relocate the pen",              DoRelocate,
  DRAW,     "draw",     "draw a line to a point",        DoDraw,
  LTYPE,    "ltype",    "set the line type for drawing", DoLtype,
  LABEL,    "label",    "draw a label on the display",   DoLabel,

  WARRANTY, "warranty", "describe no warranty",        DoWarranty,
  COPYRIGHT,"copyright","describe the copyright",      DoCopyright,
  SET,      "set",      "set a solution field",        DoSet,

  EXEC,     "!",        "execute a shell command",     NULL,
  SLEEP,    "sleep",    "sleep for a while",           NULL,
  COMMAND,  "command",  "return to command mode",      NULL,
  HELP,     "help",     "list available commands",     NULL,
  QUIT,     "quit",     "quit the program",            NULL,
  0,         NULL,       NULL,                         NULL
};

static Tree* token_tree;

/* Private Functions */

/* Copy the edge points + offset into the arrays (xb,yb) */

static void get_edge_pts 
   (Element *U, Edge *edge, double offset, double xb[], double yb[])
{
  const int np    = edge->np;
  const int start = edge->start;
  const int skip  = edge->skip;

  double *nx = edge->unx;
  double *ny = edge->uny;
  double *x  = *U->xmesh;
  double *y  = *U->ymesh;

  int i, j;

  for (i = 0, j = start; i < np; i++, j += skip) {
    xb [i] = offset * nx[i] + x[j];
    yb [i] = offset * ny[i] + y[j];
  }
}

/* Grab all of the boundary points for an element */

static void get_boundary (Element *U, double xb[], double yb[])
{
  const int nb = 2 * (U->nr + U->ns - 2);

  dgathr (nb, *U->xmesh, U->emap, xb);
  dgathr (nb, *U->ymesh, U->emap, yb);

  xb [nb] = xb[0];  /* Wrap around the edges */
  yb [nb] = yb[0];
}
  
static int level (int key) {
  int n = 0;
  while (key >>= 2) n++;
  return n;
}

/* Split a token list < ... > into separate strings */

static int token_split (int nmax, char *list, char *tokens[])
{
  char *p = list;
  int   n = 0;

  while (*p != '\0' && n <= nmax) {
    char *q = p;
    char buf[BUFSIZ];

    while (*q != ',' && *q != '\0') {
      buf[q-p] = *q;
      q++;
    }
    buf   [q-p] = '\0';
    tokens[n++] = strdup(buf);
    p = (*q != '\0') ? ++q : q;
  }

  return n;
}

/* Strip leading and trailing whitespace */

static char *strip(char *s)
{
  char *p = s;
  char *t = p + strlen(s)-1;

  while (isspace(*(s=p)))
    *(p++)='\0';
  while (isspace(*t))
    *(t--)='\0';

  return s;
}

/* ------------------------------------------------------------------------- *
 * DoParse                                                                   *
 *                                                                           *
 * Top-level parsing function.  On the first call, DoParse() initializes the *
 * symbol table of commands.  It then looks up a command and branches to an  *
 * action routine.   There are two ways to handle a command:                 *
 *                                                                           *
 *   1.  If the given token has an associated function call, then that func- *
 *       tion is called.  It has argument (void) but can read options using  *
 *       strtok (NULL, separators) from the command line.                    *
 *                                                                           *
 *   2.  If the given token does not have an associated function, then it    *
 *       should be handled in-line by the switch that follows.  All commands *
 *       that require arguments can be handled this way.                     *
 *                                                                           *
 * ------------------------------------------------------------------------- */

/* Print the help string for a token */

void show_token (Token *token) {
  if (token->help)
    printf ("%-10s ----- %s\n", token->name, token->help);
}

/* Clear the command line.  This is necessary so that any functions that    *
 * use strtok() for parsing their own strings can return to the main parser *
 * with a clean state.                                                      */

void clear_command_line () {
  while (strtok(NULL," \t;\n") != NULL)
    continue;
}

/* ------------------------------------------------------------------------- *
 * expand_tokens()                                                           *
 *                                                                           *
 * This function expands tokens on the command line.  Two forms of $expr     *
 * are parsed:                                                               *
 *                                                                           *
 *     $token      Expands a number-valued token from the internal symbol    *
 *                 tables                                                    *
 *                                                                           *
 *     ${token}    Expands an environment variable                           *
 *                                                                           *
 * ------------------------------------------------------------------------- */

void expand_tokens (char *orig, char *expanded)
{
  char *p = orig;
  char *s = expanded;

  char buf[BUFSIZ];

  while (*p) {
    if (*p=='$') {
      char *wp = ++p;
      if (*wp=='{') {
	for (wp = ++p; *wp != '}'; wp++)
	  buf[wp-p] = *wp;
	buf[wp-p] = '\0';

	p  = ++wp;
	wp = getenv(buf);
	strcpy(s,wp);
	s += strlen(wp);
	
      } else {
	char valu[64];
	
	while ((isalpha(*wp) || *wp==':')) {
	  buf[wp-p] = *wp;
	  wp++;
	}

	buf[wp-p] = '\0';
	sprintf(valu,"%.16g",dparam(buf));
	strcpy (s, valu);
	s += strlen(valu);
	p  = wp;
      }
    }
    *s++ = *p++;
  }

  *s = '\0';
}

/* Parse a string into words separated by characters from the input string   *
 * "delim".  This works just like strtok() except that get back a pointer to *
 * the end of the current word so that you can step through the string w/out *
 * stepping on strtok's static memory.                                       */

static char *get_word (char *p, char **endp, const char *delim)
{
  char *ptr = p;
  char *beg = p;

  while (*ptr != '\0' && strchr(delim,*ptr))   /* skip leading spaces */
    ptr++;

  beg = ptr;  /* save pointer to the beginning of this word */

  while (*ptr != '\0' && !strchr(delim,*ptr))  /* search to end of word */
    ptr++;

  *ptr  = '\0';   /* mark word end */
  *endp = ++ptr;  /* save next word */

  return beg;
}


int DoParse (char *command_line)
{
  Node*  np;
  char*  command;
  Token  token;
  int    status = 0;

  char buf[BUFSIZ];

  /* Initialize the tree on the first call */

  if (!token_tree) {
    int i;
    token_tree = create_tree (show_token, (PFV) free);

    for (i = 0; token_list[i].name; i++) {
      (np = create_node (token_list[i].name))->other 
	  = & token_list[i];
      tree_insert (token_tree, np);
    }
  }

  expand_tokens(command_line, buf);

  command = strtok(buf, " \t;\n");   /* Grab the first token */

  while (command) {

    if (*command == '#') {    /* Check for comment line */
      clear_command_line();
      command = NULL;
      continue;
    }

    if ((np = tree_search (token_tree->root, command))) {
      Token *foo = (Token*) np->other;      /* ANSI C workaround */
      token = *foo;
      
      if (token.action) 
	(*token.action) ();
      else 
	switch (token.type) {
	case EXEC:
	  system(strtok(NULL,"\n"));
	  break;
	case SLEEP:
	  sleep (atoi(strtok(NULL, " ;")));
	  break;
	case UNZOOM:
	  SetLimits (&Geometry);
	  break;
	case COMMAND:
	  clear_command_line();
	  status = 1;
	  break;
	case HELP:
	  printf ("Current commands:\n");
	  tree_walk (token_tree);
	  break;
	case QUIT:
	  return EOF;
	default:
	  break;
	}
    } else
      printf ("Mscope: unknown command: %s\n", command);
    
    command = strtok(0, " ;\n");  /* Grab the next command */
  }
  
  return status;
}

/* ------------------------------------------------------------------------- *
// Action Commands                                                           //
//                                                                           //
// Each of the following functions implements one "keyword".                 //
 * ------------------------------------------------------------------------- */

void DoMnodes() { return; }
void DoSnodes() { return; }

/* ------------------------------------------------------------------------- */

void DoBox () {
  pl_box (1, 2, 0, 0);
  pl_gflush();
}

void DoInput ()
{
  char *p = strtok(NULL, " ;\n");

  if (p) {
    char buf[BUFSIZ];
    FILE *fp = fopen(p, "r");
    FILE *save_command_stream;

    extern FILE *mscope_command_stream;

    if (!fp) {
      fprintf (stderr, "mscope: can't open script -- %s\n", p);
      return;
    }

    printf("Reading commands from file %s\n", p);
    
    save_command_stream   = mscope_command_stream; 
    mscope_command_stream = fp;

    while (fgets(buf, BUFSIZ, fp))
      if (DoParse(buf))
	break;
    
    fclose(fp);
    mscope_command_stream = save_command_stream;
  } else
    printf("usage: input command_file\n");
}

void DoDomain()
{
  char *p = strtok(NULL, " ;\n");
  if (p) {
    Domain_init(p);
  } else 
    fprintf (stderr, "usage: domain name\n");
}

static void MeshIn_prism_1 (char *session) 
{
  Mesh *mesh;
  FILE *fp;
  char filename [FILENAME_MAX];
  char buf [BUFSIZ];

  Domain_reset();

  strcpy (filename, session);
  if (!(fp = fopen (filename, "r"))) {
    strcat(filename, ".rea");
    if (!(fp = fopen(filename, "r"))) {
      fprintf (stderr, "Unable to load %s[.rea]\n", session);
      memset (&Geometry, '\0', sizeof(Domain));
      return;
    }
  }

  /* Load the parameters */

  Geometry.name  = strdup(session);
  Geometry.param = param_alloc(32);
  param_import(Geometry.param, fp);
  
  /* Allocate a new mesh and load its definition from a file */
  
  Mesh_import (mesh = Mesh_alloc(), fp);
  
  /* Set up the computational domain */
  
  Geometry.mesh     = mesh;
  Geometry.A        = Matrix_alloc();
  Geometry.force    = NULL;
  Geometry.solution = NULL;
  
  Domain_bbox();
  SetLimits(&Geometry);
  
  /* Try loading the other keywords */
  
  rewind(fp); if (findSection("HISTORY", buf, fp)) {
    int i, npts;
    Geometry.history = keyword_alloc("HISTORY POINTS");

    /* Import the data but modify each line by deleting the 'P' flag*/

    keyword_import (Geometry.history, fp);

    npts = KEYWORD_COUNT(Geometry.history);
    for (i = 0; i < npts; i++) {
      char *s1, *s2;
      for (s1 = KEYWORD_INFO(Geometry.history,i); !isspace(*s1); s1++)
	continue;
      for (s2 = s1+1; !isspace(*s2); *s2++)
	continue;
      memmove(s1, s2, strlen(s2)+1);
    }
  }
  
  rewind(fp); if (findSection("DRIVE FORCE", buf, fp)) {
    Geometry.force = keyword_alloc("DRIVE FORCE");
    keyword_import (Geometry.force, fp);
  }
  
  rewind(fp); if (findSection("SOLUTION", buf, fp)) {
    Geometry.solution = keyword_alloc("SOLUTION");
    keyword_import (Geometry.solution, fp);
  }
  
  rewind(fp); if (findSection("INITIAL CONDITIONS", buf, fp)) {
    Geometry.ic = keyword_alloc("INITIAL CONDITIONS");
    keyword_import (Geometry.ic, fp);
  }
  
  fclose(fp);
}

/* Read a FEML file */

static void foo (Element *elmt)
{
  int nr = ELEMENT_NR(elmt);
  int ns = ELEMENT_NS(elmt);
  int i, j;
  
  printf("[%d]\n", ELEMENT_ID(elmt));
  for (i = 0; i < ns; i++) {
    for (j = 0; j < nr; j++) {
      printf ("%#10.6g ", elmt->xmesh[ns-1-i][j]);
    }
    printf ("\n");
  }
}

#define MAX_GROUPS          8
#define MAX_GROUPS_STRLEN  64
#define MAX_SURFACES       64
#define MAX_BCS             8

static void MeshIn_semtex (char *session) 
{
  Mesh *mesh;
  FILE *fp;
  char buf[BUFSIZ];
  char *p;

  struct {
    double x;
    double y;
  } node[_MAX_NEL*4];

  struct {
    int family;
    int key;
  } etab[_MAX_NEL];

  struct {
    char tag;
    char name[MAX_GROUPS_STRLEN];
  } gtab[MAX_GROUPS];

  BC *btab[MAX_BCS];

  Domain_reset();

  Geometry.name  = strdup(session);
  Geometry.param = param_alloc(32);

  if (!(fp=fopen(session,"r"))) {
    fprintf (stderr, "Unable to load %s\n", session);
    memset(&Geometry, '\0', sizeof(Domain));
    return;
  } else {
    printf ("Loading %s -- please wait\n", session);

    rewind(fp); if (findSection("<TOKENS", buf, fp)) {
      char token[64];
      char expr [64];
      int  done = 0;

      Geometry.param = param_alloc(32);

      do {
	fgets(buf, BUFSIZ, fp);
	
	if (strchr(buf, '=')) {
	  sscanf(buf, "%s%*s%s", token, expr);

	  param_define(Geometry.param, token, expr);

	  if (strcmp(token,"N_POLY")==0)
	    iparam_set("NORDER", iparam("N_POLY"));
	  if (strcmp(token,"N_Z")==0)
	    iparam_set("NZ", iparam("N_Z"));
	}

	if (strstr(buf, "</TOKENS>")) done = 1;

      } while (!done);
    }

    rewind(fp); if (findSection("<USER", buf, fp)) {

      Geometry.user = keyword_alloc("USER");

      fgets(buf, BUFSIZ, fp);
      while (!strstr(buf,"/USER")) {
	keyword_add(Geometry.user, strip(buf));
	fgets(buf, BUFSIZ, fp);
      }
    }

    rewind(fp); if (findSection("<FIELDS", buf, fp)) {
      int i;

      Geometry.fields = keyword_alloc("FIELDS");

      fgets(buf, BUFSIZ, fp);
      while (!strstr(buf, "</FIELDS>")) {
	char *info[16];
	int  n = token_split(16, buf, info);
	for (i = 0; i < n; i++)
	  keyword_add(Geometry.fields, strip(info[i]));
	fgets(buf, BUFSIZ, fp);
      }
    }

    rewind(fp); if (findSection("<HISTORY", buf, fp)) {

      Geometry.history = keyword_alloc("HISTORY");

      fgets(buf, BUFSIZ, fp);
      while (!strstr(buf, "</HISTORY>")) {
	double x, y, z;
	sscanf (buf, "%*d%lf%lf%lf", &x, &y, &z);
	sprintf(buf, "u %.16g %.16g %.16g", x, y, z);
	keyword_add(Geometry.history, buf);
	fgets(buf, BUFSIZ, fp);
      }
    }

    /* ------------------------------------------------------------------------- */

    Geometry.mesh = mesh = Mesh_alloc();
    Geometry.A    = Matrix_alloc();

    rewind(fp); if (findSection("<NODES", buf, fp)) {
      int i, id, n;
      double x, y, z;

      if (!(p=strchr(buf,'='))) {
	fprintf (stderr, "NUMBER attr missing from NODES section\n");
	return;
      }

      n = atoi(++p);

      for (i = 0; i < n; i++) {
	fscanf(fp, "%d%lf%lf%lf", &id, &x, &y, &z);
	node[id].x = x;
	node[id].y = y;
      }
    }

    rewind(fp); if (findSection("<ELEMENTS", buf, fp)) {
      int n;
      int id;
      int inode[4];
      int i, k;

      const int nr = iparam("NORDER");
      const int ns = iparam("NORDER");
      const int nz = iparam("NZ");

      double xc[4], yc[4];

      if (!(p=strchr(buf,'='))) {
	fprintf (stderr, "NUMBER attr missing from ELEMENTS section\n");
	return;
      }

      n = atoi(++p);

      for (k = 0; k < n; k++) {
	Element *elmt;

	fscanf(fp, "%d%*s%d%d%d%d", &id, inode, inode+1, inode+2, inode+3);
	fgets (buf, BUFSIZ, fp);

	elmt = Element_alloc(nr, ns, nz);    /* allocate a new element */
	for (i = 0; i < 4; i++) {            /* extract its nodes */
	  xc[i] = node[inode[i]].x;
	  yc[i] = node[inode[i]].y;
	}

	Element_setGeometry(elmt, xc, yc);   /* initialize the geometry */
	Mesh_add(mesh,elmt);                 /* add it to the mesh */

	/* Since there is no guarantee about the element id's assigned by    *
	 * the mesh, we store the element's (family,key) pair in the table   *
	 * of ids so we can x-reference curves and BCs.                      */

	etab[id].family = ELEMENT_FAMILY(elmt);
	etab[id].key    = ELEMENT_KEY   (elmt);
      }
    }

    rewind(fp); if (findSection("<CURVES", buf, fp)) {
      int  i, id, k, n, side;
      char tag[8];
      Element *elmt;
      Edge    *edge;
      Curve   *c;

      if (!(p=strchr(buf,'='))) {
	fprintf (stderr, "NUMBER attr missing from CURVES section\n");
	return;
      }

      n = atoi(++p);

      for (i = 0; i < n; i++) {
	fscanf(fp, "%*d%d%d%s", &k, &side, tag);
	
	side--;
	id   = Mesh_lookup(mesh, etab[k].family, etab[k].key);
	elmt = MESH_ELEMENT(mesh,id);
	edge = ELEMENT_EDGE(elmt,side);
	c    = edge->curve = (Curve*) calloc(1,sizeof(Curve));

	if (strcmp(tag,"<ARC>")==0) {
	  double radius;
	  fscanf (fp, "%lf", &radius);
	  c->type = Arc;
	  c->info.arc = make_arc (elmt, edge, radius);
	  Element_genxy(elmt);
	} else {
	  fprintf (stderr, "uknown CURVE tag -- %s\n", tag);
	  return;
	}

	fgets(buf, BUFSIZ, fp);
      }
    }

    /* ------------------------------------------------------------------------- */

    rewind(fp); if (findSection("<GROUPS", buf, fp)) {
      int i, n;

      if (!(p=strchr(buf,'='))) {
	fprintf (stderr, "NUMBER attr missing from GROUPS section\n");
	return;
      }

      n = atoi(++p);

      memset(gtab, '\0', sizeof(gtab));

      for (i = 0; i < n; i++) {
	fgets (buf, BUFSIZ, fp);
	sscanf(buf, "%*d%1s%s", &gtab[i].tag, gtab[i].name);
      }
    }

    /* ------------------------------------------------------------------------- */

    rewind(fp); if (findSection("<BCS", buf, fp)) {
      int i, j, n;

      if (!(p=strchr(buf,'='))) {
	fprintf (stderr, "NUMBER attr missing from BCS section\n");
	return;
      }

      n = atoi(++p);

      memset(btab, '\0', sizeof(btab));

      for (i = 0; i < n; i++) {
	char group[2];
	int  nitems;

	fgets (buf, BUFSIZ, fp);
	sscanf(buf, "%*d%1s%d", group, &nitems);

	btab[i] = BC_alloc(*group, nitems);
	
	for (j = 0; j < nitems; j++) {
	  char tagstr[5];
	  char tag;
	  char *p;

	  fscanf(fp, "%s", tagstr);
	  fgets (buf, BUFSIZ, fp);
	  
	  tag = tagstr[1];
	  sprintf(tagstr,"</%c>", tag);
	  if (!(p=strstr(buf,tagstr))) 
	    fprintf (stderr, "no closing tag -- %s\n", buf);
	  
	  *p = '\0'; p = strip(buf);

	  btab[i]->info[j].expr = malloc(strlen(p)+3);
	  sprintf(btab[i]->info[j].expr, "%c:%s", tag, p);
	}

	btab[i]->next = mesh->bc;      /* add to the mesh */
	mesh->bc      = btab[i];
      }
    }

    /* --------------------------------------------------------------------- *
     * We just use the SURFACES section to attach boundary conditions to     *
     * specific elements.  The section gives the element/edge pair and the   *
     * group tag for the boundary condition.  The group tag is looked up in  *
     * btab, a BC is allocated, and finally the BC is attached to the speci- *
     * fied edge.                                                            *
     * --------------------------------------------------------------------- */

    rewind(fp); if (findSection("<SURFACES", buf, fp)) {
      int  i, k, n, id, side;
      char tag[8];
      
      if (!(p=strchr(buf,'='))) {
	fprintf (stderr, "NUMBER attr missing from SURFACES section\n");
	return;
      }

      n = atoi(++p);
      
      for (i = 0; i < n; i++) {
	Element *elmt;
	Edge    *edge;

	fscanf(fp, "%*d%d%d%s", &k, &side, tag);
	
	/* Swap local id for the id assigned by the mesh */

	id    = Mesh_lookup(mesh, etab[k].family, etab[k].key);
	side -= 1;

	if (id < 0) {
	  fprintf (stderr, "unable to local element %d:(%d,%d)\n",
		   k, etab[k].family, etab[k].key);
	  return;
	} else {
	  elmt = mesh->list[id];
	  edge = ELEMENT_EDGE(elmt,side);
	}

	if (strcmp(tag,"<B>")==0) {
	  char group[2];
	  int j, ok = 0;

	  /* find the matching group */

	  fscanf(fp,"%1s", group);
	  for (j = 0; btab[j]; j++) {
	    if (btab[j]->type == *group) {
	      ok = 1;
	      break;
	    }
	  }

	  if (!ok) {
	    fprintf (stderr, "unable to locate group -- %c\n", group);
	    return;
	  }

	  BC_attach(btab[j], elmt, edge);

	} else if (strcmp(tag,"<P>")==0) {
	  struct {
	    int      k;
	    int      side;
	    int      id;
	    Element *elmt;
	    Edge    *edge;
	  } to;

	  BC *bc0 = BC_alloc('P', 2);
	  BC *bc1 = BC_alloc('P', 2);

	  fscanf(fp, "%d%d", &to.k, &to.side);
	  to.side--;

	  if ((to.id = Mesh_lookup(mesh,etab[to.k].family, etab[to.k].key))<0) {
	    fprintf (stderr, "unable to locate element %d:(%d,%d)\n",
		     to.k, etab[to.k].family, etab[to.k].key);
	    return;
	  }

	  to.elmt = MESH_ELEMENT(mesh,to.id);
	  to.edge = ELEMENT_EDGE(to.elmt,to.side);

	  if (fabs(edge->pos[0] - to.edge->pos[0]) < dparam("TOLVDB"))
	    bc0->type = bc1->type = 'Y';
	  else if (fabs(edge->pos[1] - to.edge->pos[1]) < dparam("TOLVDB"))
	    bc0->type = bc1->type = 'X';
	  else {
	    fprintf (stderr, "bad periodic boundary condition\n");
	    return;
	  }

	  bc0->info[0].number = to.id;
	  bc0->info[1].number = to.side;
	  bc0->next           = mesh->bc;
	  mesh->bc            = bc0;
	  BC_attach(bc0, to.elmt, to.edge);

	  bc1->info[0].number = ELEMENT_ID(elmt);
	  bc1->info[1].number = side;
	  bc1->next           = mesh->bc;
	  mesh->bc            = bc1;
	  BC_attach(bc1, elmt, edge);

	} else {
	  fprintf (stderr, "unknown surface tag -- %s\n", tag);
	  return;
	}

	fgets(buf, BUFSIZ, fp);  /* finish the line */
      }
    }

    Domain_bbox ();
    SetLimits   (&Geometry);
    Mesh_connect(mesh);
  }

  fclose(fp);
}

void DoMeshIn() 
{
  char *session;

  if ((session = strtok(NULL, " ;\n"))) {

    if (strstr(session,".rea"))
      MeshIn_prism_1(session);
    else
      MeshIn_semtex(session);

    if (option("autoRedraw")) {
      DoErase();
      DoGrid ();
    }
  } else {
    puts ("usage: meshin session[.rea]");
  }

  clear_command_line();
}

void DoLoad() {
  DoMeshIn();
}

void DoFlowin()
{
  if (Domain_require()) {
    const char *fname = strtok(NULL," ;\n");
    
    if (!fname) {
      fprintf (stderr, "usage: flowin file[.fld]\n");
    } else {
      FieldFile *ff = FieldFile_alloc();

      FILE *fp;
      char buf[BUFSIZ];
      int i;
      
      strcpy(buf,fname);
      
      if (!(fp=fopen(buf,"r"))) {
	strcat(buf,".fld");
	if (!(fp=fopen(buf,"r"))) {
	  fprintf(stderr, "can't open the solution file -- %s or %s\n", 
		  fname, buf);
	}
      }
      
      printf("Reading solution vector from %s\n", buf);
      
      FieldFile_read(ff, fp);

      printf("Fields = [%s]\n", ff->type);
      
      /* check dimensions */
      /* check ... */

      for (i = 0; i < FieldFile_getFieldCount(ff); i++) {
	const char type = ff->type[i];

	if (Domain_chkField(type)) {
	  Field *u = Domain_getField(type);
	  FieldFile_load (ff, u);
	} else {
	  Field *u = Domain_addField(type);
	  FieldFile_load (ff, u);
	}
      }	

      /* save header info using the session as a "namespace" */

      sprintf(buf, "%s:STEP", FIELDFILE_NAME(ff));
      iparam_set(buf,         FIELDFILE_STEP(ff));
      sprintf(buf, "%s:TIME", FIELDFILE_NAME(ff));
      dparam_set(buf,         FIELDFILE_TIME(ff));
      sprintf(buf, "%s:DT",   FIELDFILE_NAME(ff));
      dparam_set(buf,         FIELDFILE_DT  (ff));
      sprintf(buf, "%s:Re",   FIELDFILE_NAME(ff));
      dparam_set(buf,         FIELDFILE_RE  (ff));
      sprintf(buf, "%s:BETA", FIELDFILE_NAME(ff));
      dparam_set(buf,         FIELDFILE_BETA(ff));

      FieldFile_free(ff);
      fclose(fp);
    }
  }
}

void DoFlowout() 
{
  if (Domain_require()) {
    const char *fname = strtok(NULL," ;\n");

    if (!fname) {
      fprintf (stderr, "usage: flowout fname\n");
    } else {
      FieldFile *ff = FieldFile_alloc();
      FILE *fp;
      int i;
      char buf[BUFSIZ];

      if (!(fp=fopen(fname,"w")))
	fprintf (stderr, "can't open the output file -- %s\n", fname);
      
      printf ("Writing solution to %s\n", fname);

      /* restore header info */

      FieldFile_setName (ff, Geometry.name);

      sprintf(buf, "%s:STEP", FIELDFILE_NAME(ff));
      FIELDFILE_STEP(ff) = iparam(buf);
      sprintf(buf, "%s:TIME", FIELDFILE_NAME(ff));
      FIELDFILE_TIME(ff) = dparam(buf);
      sprintf(buf, "%s:DT",   FIELDFILE_NAME(ff));
      FIELDFILE_DT  (ff) = dparam(buf);
      sprintf(buf, "%s:Re",   FIELDFILE_NAME(ff));
      FIELDFILE_RE  (ff) = dparam(buf);
      sprintf(buf, "%s:BETA", FIELDFILE_NAME(ff));
      FIELDFILE_BETA(ff) = dparam(buf);

      for (i = 0; i < Geometry.nfields; i++)
	FieldFile_store (ff,Geometry.solVector[i]);

      FieldFile_write(ff,fp);
      FieldFile_free (ff);
      
      fclose(fp);
    }
  }
}


void DoGrid ()
{
  if (Domain_require()) {
    Element  *uu = Geometry.mesh->head;
    const int nr = Geometry.mesh->nr;
    const int ns = Geometry.mesh->ns;
    const int nb = 2 * (nr + ns - 2);
    
    double xb[_MAX_NB + 1];
    double yb[_MAX_NB + 1];

    pl_graphics();
    while (uu) {
      get_boundary (uu,   xb, yb);
      pl_connect   (nb+1, xb, yb);
      uu = uu->next;
    }
    
    pl_gflush();
    pl_alpha ();
  }
}

void mesh_labels(int);

void DoGridNum () { mesh_labels (GRIDNUM); }
void DoID      () { mesh_labels (ID);      }
void DoFamily  () { mesh_labels (FAMILY);  }
void DoLevel   () { mesh_labels (LEVEL);   }
void DoKey     () { mesh_labels (KEY);     }

void mesh_labels (int labeltype)
{
  if (Domain_require()) {
    const int N = Geometry.mesh->nr;
    char  buf[BUFSIZ];
    float  xp, yp;
    
    Element *uu = Geometry.mesh->head;
    
    pl_graphics();
    while (uu) {
      
      if (N & 1) {
	const int i = (N-1)/2;
	xp = uu -> xmesh[i][i];
	yp = uu -> ymesh[i][i];
      } else {
	const int i = N/2 - 1;
	xp = .25 * (uu->xmesh[ i ][i] + uu->xmesh[ i ][i+1] + 
		    uu->xmesh[i+1][i] + uu->xmesh[i+1][i+1]);
	yp = .25 * (uu->ymesh[ i ][i] + uu->ymesh[ i ][i+1] + 
		    uu->ymesh[i+1][i] + uu->ymesh[i+1][i+1]);
      }      
      
      switch (labeltype) {
      case GRIDNUM: 
	sprintf (buf, "%d", uu->id+1);
	break;
      case ID:
	sprintf (buf, "%d", uu->id);
	break;
      case FAMILY:
	sprintf (buf, "%d", uu->family);
	break;
      case KEY:
	sprintf (buf, "%d", uu->key);
	break;
      case LEVEL:
	sprintf (buf, "%d", level(uu->key));
	break;
      default:
	break;
      }
      
      pl_relocate (xp, yp);
      pl_label    (buf);
      
      uu = uu->next;
    }
    
    pl_gflush();
    pl_alpha ();
  }
}

void DoBCs()
{
  if (Domain_require()) {
    Element *elmt;
    Edge    *edge;
    BC      *bc;
    
    static char *BC_type_Dirichlet = "DTtUuVvWwQq";
    static char *BC_type_Neumann   = "AaFfINOo";
    static char *BC_type_Periodic  = "PXY";
    
    double xb[_MAX_NB], yb[_MAX_NB];
    double offset;
    int i;
    
    offset = 0.01 * MIN(scalar("mscope_xmax-mscope_xmin"),
			scalar("mscope_ymax-mscope_ymin"));
    
    pl_graphics();
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
      for (edge = elmt->edge_list; edge; edge = edge->next) {
	if ((bc = edge->bc)) {
	  const int np    = edge->np;
	  const int start = edge->offset;
	  
	  get_boundary (elmt, xb, yb);
	  
	  /* Dirichlet Boundary */
	  
	  if (strchr (BC_type_Dirichlet, bc->type)) {
	    pl_ltype (0);
	    switch (bc->type) {
	    case 'w':
	    case 'W': 
	    case 'T': {
	      dcopy (np, xb + start, 1, xb, 1);
	      dcopy (np, yb + start, 1, yb, 1);
	      
	      for (i = 0; i < np; i++) {
		xb[np*2 - 1 - i] = offset * (edge->unx)[i] + xb[i];
		yb[np*2 - 1 - i] = offset * (edge->uny)[i] + yb[i];
	      }
	      
	      xb[np*2] = xb[0];
	      yb[np*2] = yb[0];
	      
	      pl_angle(45);
	      pl_shade(np*2+1, xb, yb, 250);
	      break;
	    }
	    default:
	      pl_connect (np, xb + start, yb + start);
	      break;
	    }
	  }
	  
	  /* Neumann boundary */
	  
	  if (strchr (BC_type_Neumann, bc->type)) {
	    pl_ltype  (2);
	    pl_angle  (0);
	    pl_connect(np, xb + start, yb + start);
	  }
	  
	  /* Periodic boundary */
	  
	  if (strchr (BC_type_Periodic, bc->type)) {
	    pl_ltype  (1);
	    pl_angle  (0);
	    pl_connect(np, xb + start, yb + start);
	  }
	}
      }
    
    pl_angle (0);
    pl_ltype (0);
    pl_gflush();
    pl_alpha ();
  }
}

/* Draw the mesh */

void DoMesh ()
{
  if (Domain_require()) {
    const int nr  = Geometry.mesh->nr;
    const int ns  = Geometry.mesh->ns;
    double xp[_MAX_NORDER];
    double yp[_MAX_NORDER];
    int i;

    Element *uu = Geometry.mesh->head;

    pl_graphics();
    while (uu) {
      for (i = 0; i < ns; i++) {
	dcopy(nr, uu->xmesh[i], 1, xp, 1);
	dcopy(nr, uu->ymesh[i], 1, yp, 1);
	pl_connect (nr, xp, yp);
      }
      
      for (i = 0; i < nr; i++) {
	dcopy(ns, uu->xmesh[0]+i, nr, xp, 1);
	dcopy(ns, uu->ymesh[0]+i, nr, yp, 1);
	pl_connect (ns, xp, yp);
      }
      
      uu = uu->next;
    }
    
    pl_gflush();
    pl_alpha ();
  }
}

/* ------------------------------------------------------------------------- *
 * vnode_(): node dictionary                                                 *
 *                                                                           *
 * This set of functions manipulates a table of nodes.  Each node is simply  *
 * a labeled position in space.  A node is added to a table with a the com-  *
 * mand                                                                      *
 *                                                                           *
 *			  node <label, x, y>                                 *
 *                                                                           *
 * where label is a string (label = V0) and x,y are expressions (x=PI,y=5).  *
 * Once created, a node can be used to declare an indexed element.           *
 * ------------------------------------------------------------------------- */

static Tree* nodes;

typedef struct {
  double x;
  double y;
  char   label[32];
} VNode;

static void vnode_show (VNode *p) {
  printf ("<%s, %g, %g>\n", p->label, p->x, p->y);
}

static int vnode_lookup (char *label, double *x, double *y)
{
  int ok = 0;

  if (nodes) {
    Node *np = tree_search(nodes->root, label);
    
    if (np) {
      Point *p = (Point*) np->other;
      *x = p->x;
      *y = p->y;
      ok = 1;
    }
  }

  return ok;
}

static void vnode_add (char *label, double x, double y)
{
  Node *np;

  if (nodes == (Tree*) NULL) {
    nodes = create_tree (vnode_show, (PFV) free);
  }

  if ((np = tree_search(nodes->root, label))) {
    VNode *p = (VNode*) np->other;
    p->x = x;
    p->y = y;
    strcpy(p->label, label);
  } else {
    VNode *p = (VNode*) malloc(sizeof(VNode));
    p->x = x;
    p->y = y;
    strcpy(p->label, label);
    np = create_node(label);
    np->other = (void*) p;
    tree_insert (nodes, np);
  }
}

void DoNode ()
{
  char *p = strtok(NULL, "<>");
    
  if (!p) {
    fputs ("usage: node <label, x, y>\n", stderr);
  } else {
    int count = 1;
    char *info[3];
    
    info[0] = strip(strtok(p, ","));
    while ((p=strtok(NULL,",")) != NULL && count < 3)
      info[count++] = strip(p);
    
    if (count != 3) {
      fprintf (stderr, "node: expected 3 parameters, got %d\n", count);
      return;
    }
    vnode_add (info[0], scalar(info[1]), scalar(info[2]));
  }
}

void DoElement ()
{
  static char *help =
    "options:\n"
    "element <x, y, dx, dy>\n"
    "        create a new rectangular element whose origin is at the point\n"
    "        (x,y) and whose size is (dx,dy)\n"
    "element <+, y, dx, dy>\n"
    "        create a new rectangular element using the maximum x coordinate\n"
    "        of the previous one and a given y coordinate\n"
    "element <x, +, dx, dy>\n"
    "        create a new rectangular element using the maximum y coordinate\n"
    "        of the previous one and a given x coordinate\n"
    "element <v0, v1, v2, v3>\n"
    "        create a new indexed element, where the labels v0 to v3 refer\n"
    "        to positions stored in the node dictionary (see \"node\")\n";

  if (Domain_require()) {
    char *p = strtok(NULL, "<>");

    if (!p) {
      fputs ("usage: element [options]\n", stderr);
      fputs (help, stderr);
    } else {
      int count = 1;
      char *info[4];
      double  xc[4];
      double  yc[4];
      Element *elmt;
      int nr, ns, nz;

      info[0] = strip(strtok(p, ","));
      while ((p=strtok(NULL,",")) != NULL && count < 4)
	info[count++] = strip(p);
      
      if (count != 4) {
	fprintf (stderr, "element: expected 4 parameters, got %d\n", count);
	return;
      }

      /* Get resolution parameters from the Mesh */

      if (Geometry.mesh==NULL) 
	Geometry.mesh = Mesh_alloc();

      nr = MESH_NR(Geometry.mesh);
      ns = MESH_NS(Geometry.mesh);
      nz = MESH_NZ(Geometry.mesh);

      /* The main switch is between indexed and non-indexed specs */

      if (vnode_lookup(info[0], xc, yc)) {
	int i; for (i = 1; i < 4; i++)
	  if(!vnode_lookup(info[i], &xc[i], &yc[i])) {
	    fprintf (stderr, "element: undefined node -- %s\n", info[i]);
	    return;
	  }
      } else {
	if (*info[0] == '+') {
	  Element *elmt = Geometry.mesh->head;
	  while (elmt->next)
	    elmt = elmt->next;
	  xc[0] = elmt->xmesh[ns-1][nr-1];
	} else
	  xc[0] = scalar(info[0]);

	if (*info[1] == '+') {
	  Element *elmt = MESH_HEAD(Geometry.mesh);
	  while (elmt->next)
	    elmt = elmt->next;
	  yc[0] = elmt->ymesh[ns-1][nr-1];
	} else
	  yc[0] = scalar(info[1]);

	xc[1] = xc[0] + scalar(info[2]);
	yc[1] = yc[0];
	xc[2] = xc[1];
	yc[2] = yc[1] + scalar(info[3]);
	xc[3] = xc[0];
	yc[3] = yc[2];
      }

      /* done */

      elmt = Element_alloc(nr, ns, nz);
      Element_setGeometry(elmt,xc, yc);
      Mesh_add (Geometry.mesh, elmt);
      Domain_bbox();

      if (option("autoConnect"))
	Mesh_connect(Geometry.mesh);
      if (option("autoRedraw")) {
	SetLimits (&Geometry);
	DoErase();
	DoGrid ();
      }
    }
  } else
    clear_command_line();
}

void DoCurve()
{
  static char *help = 
    "options:\n"
    "curve <edge, arc, radius>\n"
    "        create an arc through the edge endpoints; edge will be convex\n"
    "        if radius > 0 and concave if radius < 0\n"
    "curve <edge, bezier, x0, y0, x1, y1>\n"
    "        create a Bezier spline using the edge endpoints and two add-\n"
    "        tional control points in the center\n"
    "curve <edge, fit, filename>\n"
    "        create a curved edge by fitting it to a parametric geometry\n"
    "        specified in a file.  The file should contain a list of (x,y)-\n"
    "        coordinates that define a smooth shape like an airfoil\n";

  if (Domain_require()) {

    char *p = strtok(NULL, "<>");

    if (!p) {
      fputs ("usage: curve [options]\n", stderr);
      fputs (help, stderr);
    } else {
      int   id    = atoi(strtok(p,","));
      char *type  = strip(strtok(NULL,","));

      Edge    *edge;
      Curve   *curve;
      Element *elmt;
      
      if (strcmp(type,"arc")==0) {
	if (!(p=strtok(NULL,","))) {
	  fprintf (stderr, "curve: arc: expected 1 parameter\n");
	  return;
	}

	curve = Curve_alloc(Arc, scalar(p));

      } else if (strcmp(type,"bezier")==0) {
	double pos[4];
	int i; for (i = 0; i < 4; i++) {
	  if (!(p=strtok(NULL,","))) {
	    fprintf (stderr, "curve: bezier: expected 4 parameters\n");
	    return;
	  } else 
	    pos[i] = scalar(p);
	}

	curve = Curve_alloc(Spline, pos[0], pos[1], pos[2], pos[3]);

      } else if (strcmp(type,"fit")==0) {
	if (!(p = strtok(NULL,","))) {
	  fprintf (stderr, "curve: fit: expected 1 parameter\n");
	  return;
	} 
	
	curve = Curve_alloc(File, strip(p));
      } else {
	fprintf (stderr, "curve: unknown type -- %s\n", type);
	return;
      }

      elmt = Geometry.mesh->head;
      while (elmt->next)
	elmt = elmt->next;
      edge = &elmt->edge_list[id];

      /* Now attach the curve to an edge and regenerate the element */

      Curve_attach (curve, elmt, edge);
      Element_genxy(elmt);

      if (option("autoScale")) {
	Domain_bbox();
	SetLimits (&Geometry);
      }
      if (option("autoRedraw")) {
	DoErase();
	DoGrid ();
      }
    }
  } else
    clear_command_line();
}

void DoHistory()
{
  static char *help = 
    "history <typelist, x, y, z>\n"
    "<typelist> = list of fields to sample\n"
    "x, y, z    = position of the history point in space\n";
    
  if (Domain_require()) {
    char *p = strtok(NULL,"<>");
    if (!p) {
      fprintf (stderr, "usage: %s\n", help);
      return;
    } else {
      char *info[4];
      char  buf[BUFSIZ];
      int i;

      info[0] = strip(strtok(p,","));
      for (i = 1; i < 3; i++) {
	if (!(p=strtok(NULL,","))) {
	  fprintf (stderr, "usage: %s\n", help);
	  return;
	} else
	  info[i] = strip(p);
      }

      /* The frame specification is optional */

      if (!(p=strtok(NULL,",")))
	info[3] = NULL;
      else
	info[3] = strip(p);

      /* Allocate history points */

      if (!Geometry.history) 
	Geometry.history = keyword_alloc("HISTORY");

      sprintf(buf, "%s %.16g %.16g %s",
	      info[0], 
	      scalar(info[1]), 
	      scalar(info[2]), 
	      info[3] ? info[3] : "0");

      keyword_add(Geometry.history, buf);
    }
  }
}

void DoBcond() 
{
  if (Domain_require()) {
    char *p = strtok(NULL, "<>");
    
    if (!p)
      fputs ("usage: bcond <edge, type, ...>\n", stderr);
    else {
      char *info[16];
      char  type;
      int   i, id, nitems;

      BC      *bc;
      Edge    *edge;
      Element *elmt;

      int n = token_split(16, p, info);

      if (n < 2) {
	fprintf (stderr, "bcond: please specify <edge> and <type>\n");
	return;
      }

      id     = atoi(info[0]);
      type   = *strip(info[1]);
      nitems = n-2;

      bc = BC_alloc(type, nitems);
      if (isupper(type))
	for (i = 0; i < nitems; i++)
	  BC_INFO(bc,i).value = scalar(info[i+2]);
      else
	for (i = 0; i < nitems; i++)
	  BC_INFO(bc,i).expr = strdup(strip(info[i+2]));

      /* Now attach to an edge */

      elmt = Geometry.mesh->head;
      while (elmt->next)
	elmt = elmt->next;
      edge = &elmt->edge_list[id];

      BC_attach(bc,elmt,edge);

      for (i = 0; i < n; i++)
	free(info[i]);
    }
  } else
    clear_command_line();
}

void DoIcond() 
{
  if (Domain_require()) {
    char *p = strtok(NULL, "<>");

    if (!p)
      fputs ("usage: icond <type, ...>\n", stderr);
    else {
      int   i;
      char *info[16];
      char  buf [BUFSIZ];

      int   n    = token_split(16, p, info);
      char *type = info[0];

      if (Geometry.ic) keyword_free(Geometry.ic);

      Geometry.ic = keyword_alloc("INITIAL CONDITIONS");
      keyword_add(Geometry.ic, type);

      switch (toupper(*type)) {

      case 'G': 
	for (i = 1; i < n; i++) {
	  sprintf (buf, "\tv%d = %s", i-1, strip(info[i]));
	  keyword_add(Geometry.ic, buf);
	}
	break;

      case 'R': /* as in Restart */
      case 'D': /* as in Default */
      default:
	for (i = 1; i < n; i++) {
	  sprintf (buf, "\t%s", strip(info[i]));
	  keyword_add(Geometry.ic, buf);
	}
	break;
      }

      for (i = 0; i < n; i++)
	free(info[i]);
    }
  } else
    clear_command_line();
}

void DoKeys ()
{
  if (Domain_require()) {
    char *p = strtok(NULL, " \n");
    FILE *fp;
    int family, key;
    char buf[BUFSIZ];
    
    if (!p) {
      fprintf (stderr, "usage: keys session.rea\n");
      return;
    }
    
    if (!(fp = fopen(p, "r"))) {
      fprintf (stderr, "DoKeys: unable to open keyfile\n");
      return;
    }
    
    if (!findSection("KEYS", buf, fp)) {
      fprintf (stderr, "DoKeys: no keys in the keyfile\n");
      return;
    }
    
    while (fgets(buf, BUFSIZ, fp)) {
      if (sscanf (buf, "%d%d", &family, &key) != 2)
	break;
      Mesh_install (Geometry.mesh, family, key);
    }
  }  
}

void DoInstall ()
{
  if (Domain_require()) {
    int family = -1;
    int key    = -1;
    char *p;

    if ((p = strtok(NULL, " \n")))
      family = atoi(p);
    if ((p = strtok(NULL, " \n")))
      key    = atoi(p);
    
    if (family < 0 || key < 0)
      printf ("usage: install family key\n");
    else 
      Mesh_install (Geometry.mesh, family, key);
  }
}


void DoPrune ()
{
  if (Domain_require()) {
    int family = -1, maxlevel = -1;
    char *p;
    
    if ((p = strtok(NULL, " \n")))
      family   = atoi(p);
    if ((p = strtok(NULL, " \n")))
      maxlevel = atoi(p);
    
    if (family < 0 || maxlevel < 0)
      printf ("usage: prune family maxlevel\n");
    else {
      printf ("Prune: setting maximum level for family %d to %d\n", 
	      family, maxlevel);
      Mesh_prune (Geometry.mesh, family, maxlevel);
    }
  }
}

void DoRestrict () {
  Mesh_restrict(Geometry.mesh);
  return;
}

/* ------------------------------------------------------------------------- */

void DoRefine ()
{
  if (Domain_require()) {
    char *p;

    Mesh_resetFlags(Geometry.mesh);
    
    if ((p = strtok (NULL, " \n"))) {
      do {
	const int id = atoi(p);
	Element *elmt;
	
	if (id < 0) 
	  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
	    Element_setRefineOn(elmt);
	else {
	  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
	    if (elmt->id == id) {
	      Element_setRefineOn(elmt);
	      break;
	    }
	  }
	  if (elmt == (Element*) NULL)
	    printf ("Refine: element %d not found\n", id);
	}
	
      } while ((p = strtok (NULL, " \n")));
      
      Mesh_refine (Geometry.mesh);

      if (option("autoRedraw")) {
	DoErase();
	DoGrid ();
      }

    } else printf ("usage: refine id\n");
  }
}

void DoSplit()
{
  char *p = strtok(NULL,":\n");
  static char *help = "usage: split id axis [sigma=0.5]\n"
    "Split element <id> along a given <axis>=[0,1], using a grading factor \n"
    "<sigma>, where 0 < sigma < 1 gives the location of the split.\n";

  if (!p) 
    fputs(help,stderr);
  else {
    int id;
    int axis;
    double sigma = 0.5;
    char *ptr;
    char *endp = p;
    const char *delim = " \t";

    Mesh_resetFlags(Geometry.mesh);

    if ((ptr = get_word(endp,&endp,delim))) id    = atoi  (ptr);
    if ((ptr = get_word(endp,&endp,delim))) axis  = atoi  (ptr);
    if ((ptr = get_word(endp,&endp,delim))) sigma = scalar(ptr);

    do {
      Element *elmt;
      
      if (id < 0) 
	for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
	  Element_setRefineOn(elmt);
      else {
	for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
	  if (elmt->id == id) {
	    Element_setRefineOn(elmt);
	    break;
	  }
	}
	if (elmt == (Element*) NULL)
	  printf ("split: element %d not found\n", id);
      }
    } while ((p = strtok (NULL, " \n")));
     
    fprintf (stderr, "split <%d,%d,%.16g>\n", id, axis, sigma);

    if (option("autoRedraw")) {
      DoErase();
      DoGrid ();
    }
  }
}

void DoGradient ()
{
  char *p = strtok(NULL, ";\n");

  if (!p) {
    printf ("usage: gradient <type> eps maxDepth maxIterations\n");
    return;
  }

  if (Domain_require()) {
    Mesh    *mesh  = Geometry.mesh;
    error_t *error = Error_alloc (mesh, E_GRADIENT);
    Field   *u     = NULL;

    int  maxit    = 0;
    int  maxdepth = 0;
    double tol    = 1.;
    char type;

    sscanf(p, "%1s%lf%d%d", &type, &tol, &maxdepth, &maxit);
    
    if (maxit < 0)
      maxit = 1000;
    if (maxdepth < 0)
      maxdepth = 15;

    Mesh_resetFlags(Geometry.mesh);
    
    if ((u = Domain_getField(type))) {
      int nref, pass = 0;
      
      do {
	Error_compute (error, u); 
	Error_info(error, tol, stdout);
	nref = Error_adapt (error, tol, maxdepth);
	
	DoErase ();
	DoGrid  ();
      } while (pass++ < maxit && nref > 0);
    }
    
    Error_free (error);
  }
}

void DoSpectrum ()
{
  char *p = strtok(NULL, ";\n");

  if (!p) {
    printf ("usage: spectrum <type> eps maxDepth maxIterations\n");
    return;
  }

  if (Domain_require()) {
    Mesh    *mesh  = Geometry.mesh;
    error_t *error = Error_alloc (mesh, E_SPECTRUM);
    Field   *u     = NULL;

    int  maxit    = 0;
    int  maxdepth = 0;
    double tol    = 1.;
    char type;

    sscanf(p, "%1s%lf%d%d", &type, &tol, &maxdepth, &maxit);
    
    if (maxit < 0)
      maxit = 1000;
    if (maxdepth < 0)
      maxdepth = 15;

    Mesh_resetFlags(Geometry.mesh);
    
    if ((u = Domain_getField(type))) {
      int nref, pass = 0;
      
      do {
	Error_compute (error, u); 
	Error_info(error, tol, stdout);
	nref = Error_adapt (error, tol, maxdepth);

      } while (pass++ < maxit && nref > 0);
    }
    
    Error_free (error);
  }
}

void DoRegress ()
{
  char *p = strtok(NULL, ";\n");

  if (!p) {
    printf ("usage: regress <type> eps maxDepth maxIterations\n");
    return;
  }

  if (Domain_require()) {
    Mesh    *mesh  = Geometry.mesh;
    error_t *error = Error_alloc (mesh, E_REGRESS);
    Field   *u     = NULL;

    int  maxit    = 0;
    int  maxdepth = 0;
    double tol    = 1.;
    char type;

    sscanf(p, "%1s%lf%d%d", &type, &tol, &maxdepth, &maxit);
    
    if (maxit < 0)
      maxit = 1000;
    if (maxdepth < 0)
      maxdepth = 15;

    Mesh_resetFlags(Geometry.mesh);
    
    if ((u = Domain_getField(type))) {
      int nref, pass = 0;
      
      do {
	Error_compute (error, u); 
	Error_info(error, tol, stdout);
	nref = Error_adapt (error, tol, maxdepth);

      } while (pass++ < maxit && nref > 0);
    }
    
    Error_free (error);
  }
}


/* Draw the logical structure (vertices, edges, elements) of the mesh */

void DoLogical ()
{
  if (Domain_require()) {
    Element *U = Geometry.mesh->head;
    const int    N = Geometry.mesh->nr;
    Edge    *edge;
    Vertex  *vert;
    double  x[_MAX_NB];
    double  y[_MAX_NB];
    int i, j, p;
    
    pl_graphics();
    while (U) {
      
      pl_expand (1.);
      pl_angle  (0.);
      
      /* Vertices */
      
      pl_ptype (PL_FILLED_CIRCLE);
      for (vert = U->vert_list, p = 0; vert; vert = vert->next) { 
	if (vert->alias == NULL) {
	  x[p] = vert->pos[0]; 
	  y[p] = vert->pos[1]; 
	  p++;
	}
      }
      pl_points (p, x, y);
      
      /* Edges */
      
      for (edge = U->edge_list; edge; edge = edge->next) {
	if ((edge->alias ? edge->alias->type : 0) == 0) {
	  get_edge_pts (U, edge, 0., x, y);
	  pl_connect   (edge->np-2, x+1, y+1);
	}
      }
      
      /* Element */
      
      p = 0;
      for (i = 1, j = 1; j < N-2; j++, p++)
	{ x[p] = U->xmesh[i][j]; y[p] = U->ymesh[i][j]; }
      for (i = 1, j = N-2; i < N-2; i++, p++)
	{ x[p] = U->xmesh[i][j]; y[p] = U->ymesh[i][j]; }
      for (i = j = N-2; 0 < j; j--, p++)
	{ x[p] = U->xmesh[i][j]; y[p] = U->ymesh[i][j]; }
      for (i = N-2, j = 1; 0 < i; i--, p++)
	{ x[p] = U->xmesh[i][j]; y[p] = U->ymesh[i][j]; }
      
      pl_angle (45.);
      pl_shade (p, x, y, 500);
      
      U = U->next;
    }
    
    pl_angle (0.);
    pl_expand(1.);
    pl_gflush();
    pl_alpha ();
  }
}

void DoConnect() {
  if (Domain_require())
    Mesh_connect (Geometry.mesh);
}

void DoMinfo() {
  if (Domain_require()) 
    Mesh_info (Geometry.mesh);
}

void DoEinfo()
{
  if (Domain_require()) {
    int id;
    Element *elmt;
    Edge    *edge;
    Vertex  *vert;
    
    id   = atoi(strtok(NULL, " \n"));
    
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
      if (elmt->id == id) {
	const int nr = ELEMENT_NR(elmt);
	const int ns = ELEMENT_NS(elmt);
	const int nz = ELEMENT_NZ(elmt);
	
	const double area = dsum(nr*ns, *elmt->mass, 1);

	printf ("Element    : %d\n", id);
	printf ("Resolution : %d %d %d\n", nr, ns, nz);
	printf ("Area       : %g\n", area);

	printf ("Vertices\n");
	for (vert = elmt->vert_list; vert; vert = vert->next) {
	  printf ("\t%d  %3d  %14.7g %14.7g ", 
		  vert->id, vert->key, vert->pos[0], vert->pos[1]);
	  if (vert->alias) printf ("aliased to %d", vert->alias->key.edge);
	  printf ("\n");
	}
	
	printf ("Edges\n");
	for (edge = elmt->edge_list; edge; edge = edge->next) {
	  printf ("\t%d  %3d  %14.7g %14.7g ", 
		  edge->id, edge->key, edge->pos[0], edge->pos[1]);
	  if (edge->alias) printf ("aliased to %d", edge->alias->key.edge);
	  printf ("\n");
	  if (edge->curve)
	    switch (edge->curve->type) {
	    case Arc:
	      printf ("\t\tArc: radius %g, from %g range %g\n",
		      edge->curve->info.arc.radius,
		      edge->curve->info.arc.theta.start,
		      edge->curve->info.arc.theta.range);
	      break;
	    case Spline:
	      printf ("\t\tBezier Spline:\n");
	      break;
	    case Line:
	    case File:
	    default:
	      break;
	    }
	}
	
	break;
      }
  }
}

void DoBinfo() 
{
  if (Domain_require()) {
    Element *elmt = Geometry.mesh->head;
    printf ("Elmt Edge Type ...\n");
    while (elmt) {
      Edge *edge = elmt->edge_list;
      while (edge) {
	if (edge->bc) {
	  const BC *bc = edge->bc;
	  int i;
	  printf ("%4d %3d    %c  ", elmt->id, edge->id, bc->type);
	  switch(bc->type) {
	  case 'X':
	  case 'Y':
	  case 'O':
	  case 'W':
	    putchar('\n');
	    break;
	  default:
	    if (isupper(bc->type)) {
	      for (i = 0; i < BC_ITEMS(bc); i++) 
		printf ("%.16g ", BC_INFO(bc,i).value);
	      printf ("\n");
	    } else {
	      printf ("\n");
	      for (i = 0; i < BC_ITEMS(bc); i++) {
		char *p = strchr(BC_INFO(bc,i).expr,':');
		if (p)
		  printf ("\t%s\n", BC_INFO(bc,i).expr);
		else
		  printf ("\t%c%d = %s\n", bc->type, i, BC_INFO(bc,i).expr);
	      }
	    }
	  }
	}
	edge = edge->next;
      }
      elmt = elmt->next;
    }
  }
}

void DoParam() 
{
  const char *token = strtok(NULL, " =\n");
  const char *value = strtok(NULL, " =\n");
    
  if (Domain_require()) {
    if (!token || !value) {
      printf ("usage: param token = expr\n");
    } else {
      param_define(Geometry.param, token, value);

      /* Duplicate values between Prism and SEMTEX */

      if (strcmp(token,"N_POLY")==0)
	iparam_set("NORDER", iparam("N_POLY"));
      if (strcmp(token,"N_Z")==0)
	iparam_set("NZ", iparam("N_Z"));
    }
  }
}

void DoOption()
{
  char *token  = strtok(NULL, " \n");
  char *status = strtok(NULL, " \n");
    
  if (!token || !status) {
    printf ("usage: option <token> <on|off>\n");
  } else {
    switch(status[1]) {
    case 'n':
    case 'N':
      option_set(token,1);
      printf ("option: %s is ON\n", token);
      break;
    case 'F':
    case 'f':
      option_set(token,0);
      printf ("option: %s is OFF\n", token);
      break;
    default:
      option_set(token,scalar(status));
      break;
    }
  }
}

void DoDefine() 
{
  char *token = strtok(NULL, " =\n");
  char *expr  = strtok(NULL, " =\n");
    
  if (!token || !expr) {
    printf ("usage: define token = expr\n");
  } else
    dparam_set (token, scalar(expr));
}

static int count_curves (Element *list)
{
  int count = 0;
  Element *elmt = list;
  while (elmt) {
    Edge *edge = elmt->edge_list;
    while(edge) {
      if (edge->curve) count++;
      edge = edge->next;
    }
    elmt = elmt->next;
  }
  return count;
}

static void MeshOut_prism_1 (char *name)
{
  int id;
  FILE *fp;
    
  Element *elmt;
  Edge    *edge;
  Vertex  *vert;
  
  struct xref {         /* connectivity info */
    struct {
      int elmt;
      int edge;
    } from;
    struct {
      int elmt;
      int edge;
    } to;
  } xreftab [4*_MAX_NEL];
  
  memset (xreftab, 0, sizeof(xreftab));
  
  /* Renumber the elements */
  
  id = 1;
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
    elmt->workspace.ei[0] = elmt->id;
    elmt->id = id++;
  }
  
  /* Fill out the cross-reference table */
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
    for (edge = elmt->edge_list; edge; edge = edge->next)
      if (edge->bc == NULL && edge->alias == NULL)    /* Internal */
	if (xreftab[edge->key].from.elmt == 0) {
	  xreftab[edge->key].from.elmt = elmt->id;
	  xreftab[edge->key].from.edge = edge->id + 1;
	} else {
	  xreftab[edge->key].to.elmt = elmt->id;
	  xreftab[edge->key].to.edge = edge->id + 1;
	}
  
  /* Periodic BCs are special */
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
    for (edge = elmt->edge_list; edge; edge = edge->next)
      if (edge->bc && strchr("PXY", edge->bc->type)) {
	double pos[2];
	int key;
	
	if (edge->bc->type == 'X') {
	  pos[0] = -100000.;
	  pos[1] = edge->pos[1];
	} else {
	  pos[0] = edge->pos[0];
	  pos[1] = -100000.;
	}
	
	key = vdbdjoin(Geometry.mesh->edges,pos);
	
	if (xreftab[key].from.elmt == 0) {
	  xreftab[key].from.elmt = elmt->id;
	  xreftab[key].from.edge = edge->id;
	} else {
	  xreftab[key].to.elmt = elmt->id;
	  xreftab[key].to.edge = edge->id;
	}
      }
  
  /* Open the file (no parameters) */
  
  fp = fopen(name, "w");
  if (fp) 
    printf ("Saving current session to %s [PRISM-1]\n", name);
  else {
    printf ("Unable to open file %s for writing\n", name);
    return;
  }
  
  param_export(Geometry.param,fp);
  
  fprintf (fp, "***** MESH *****\n"
	   "%d Elements\n", MESH_NELMT(Geometry.mesh));
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
    fprintf (fp, "\tELEMENT %4d [family %d key %d]\n", 
	     elmt->id, elmt->family, elmt->key);
    for (vert = elmt->vert_list; vert; vert = vert->next)
      fprintf (fp, "%.16g ", vert->pos[0]);
    fputc ('\n', fp);
    for (vert = elmt->vert_list; vert; vert = vert->next)
      fprintf (fp, "%.16g ", vert->pos[1]);
    fputc ('\n', fp);
  }
  
  fprintf (fp, "***** CURVES *****\n"
	   "%d Curves\n", count_curves(Geometry.mesh->head));
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
    for (edge = elmt->edge_list; edge; edge = edge->next) {
      if (edge->curve) {
	fprintf (fp, "%d %d ", edge->id+1, elmt->id);
	switch (edge->curve->type) {
	case Arc:
	  fprintf (fp, "%.16g C\n", 
		   edge->curve->info.arc.theta.range > 0 ? 
		   edge->curve->info.arc.radius :
		   -edge->curve->info.arc.radius);
	  break;
	case Spline: {
	  double xp, yp;
	  spline_node(&edge->curve->info.spline, 1, &xp, &yp);
	  fprintf (fp, "%.16g %.16g ", xp, yp);
	  spline_node(&edge->curve->info.spline, 2, &xp, &yp);
	  fprintf (fp, "%.16g %.16g ", xp, yp);
	  fprintf (fp, "S\n");
	  break;
	}
	case File:
	  fprintf (fp, "%s F\n",
		   edge->curve->info.file.name);
	  break;
	case Line:
	  break;
	}
      }
    }
  }
  
  
  fprintf (fp, "***** BOUNDARY CONDITIONS *****\n");
  fprintf (fp, "***** FLUID BOUNDARY CONDITIONS *****\n");
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
    for (edge = elmt->edge_list; edge; edge = edge->next)
	
      if (edge->bc) {
	BC *bc = edge->bc;
	switch (bc->type) {
	case 'X':
	case 'Y': {
	  double pos[2];
	  int key;
	    
	  if (edge->bc->type == 'X') {
	    pos[0] = -100000.;
	    pos[1] = edge->pos[1];
	  } else {
	    pos[0] = edge->pos[0];
	    pos[1] = -100000.;
	  }
	    
	  key = vdbdpeek (Geometry.mesh->edges,pos);
	    
	  fprintf (fp, "P %d %d ", elmt->id, edge->id+1);
	  if (xreftab[key].from.elmt == elmt->id && 
	      xreftab[key].from.edge == edge->id )
	    fprintf (fp, "%d %d\n",
		     xreftab[key].to.elmt, xreftab[key].to.edge+1);
	  else
	    fprintf (fp, "%d %d\n",
		     xreftab[key].from.elmt, xreftab[key].from.edge+1);
	  break;
	}
	
	  
	  /* Symmetry: translated to an axis boundary */
	case 's':
	  fprintf (fp, "A %d %d\n", elmt->id, edge->id+1);
	  break;

	  /* WALL: must be constant velocity with u=v=w=0 */
	case 'w': 
	  fprintf (fp, "W %d %d\n", elmt->id, edge->id+1);
	  break;
	  
	  /* Outflow */
	case 'o':
	  fprintf (fp, "O %d %d\n", elmt->id, edge->id+1);
	  break;

	default:
	  fprintf (fp, "%c %d %d ", bc->type, elmt->id, edge->id+1);
	    
	  if (isupper(edge->bc->type)) {
	    for (id = 0; id < edge->bc->nitems; id++)
	      fprintf (fp, "%.16g ", edge->bc->info[id].value);
	    fputc ('\n', fp);
	  } else {
	    fputc ('\n', fp);
	    for (id = 0; id < edge->bc->nitems; id++) {
	      char *p = strchr(edge->bc->info[id].expr,':');
	      if (!p)
		fprintf (fp, "\t%c%d = %s\n", 
			 edge->bc->type, id, edge->bc->info[id].expr);
	      else
		fprintf (fp, "\t%s\n", ++p);
	    }
	  }
	}
      } else if (edge->alias) {
	switch (edge->alias->type) {
	case 0:
	  fprintf (fp, "M %d %d %d 1\n", elmt->id, edge->id+1, 
		   edge->alias->key.edge+1);
	  break;
	case 1:
	case 2:
	  fprintf (fp, "S %d %d %d %d\n", elmt->id, edge->id+1,
		   edge->alias->key.edge+1, edge->alias->type);
	  break;
	}
      } else {
	fprintf (fp, "E %d %d ", elmt->id, edge->id+1);
	  
	if (xreftab[edge->key].from.elmt == elmt->id)
	  fprintf (fp, "%d %d\n", 
		   xreftab[edge->key].to.elmt, xreftab[edge->key].to.edge);
	else
	  fprintf (fp, "%d %d\n", 
		   xreftab[edge->key].from.elmt, xreftab[edge->key].from.edge);
      }
  }
    
  if (Geometry.ic)
    keyword_export (Geometry.ic, fp);
  if (Geometry.force)
    keyword_export (Geometry.force, fp);
  if (Geometry.solution)
    keyword_export (Geometry.solution, fp);

  if (Geometry.history) {
    keyword_t *his = Geometry.history;
    const int npts = KEYWORD_COUNT(his);
    char buf[BUFSIZ];
    int i;

    fprintf (fp, "***** HISTORY POINTS *****\n");
    fprintf (fp, "%d item%c\n", npts, npts > 1 ? 's' : ' ');
    for (i = 0; i < npts; i++) {
      char *s = buf;
      strcpy (buf, KEYWORD_INFO(his,i));

      while (!isspace(*s)) s++;
      *(s++) = (char) NULL;

      fprintf (fp, "%s P %s\n", buf, s);
    }
  }      
    
  /* Active Keys */
    
  fprintf (fp, "***** KEYS *****\n");
  fprintf (fp, "%d items\n", Geometry.mesh->active);
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
    fprintf (fp, "%d %d\n", elmt->family, elmt->key);
  
  /* Restore ID numbers */
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
    elmt->id = elmt->workspace.ei[0];
  
  fclose(fp);
}

static void MeshOut_semtex (char *name)
{
  FILE *fp;
  Element *elmt;
  Edge *edge;
  int id;
  char buf[BUFSIZ];

  struct xref {         /* connectivity info */
    struct {
      int elmt;
      int edge;
    } from;
    struct {
      int elmt;
      int edge;
    } to;
  } xreftab [4*_MAX_NEL];
  
  memset (xreftab, 0, sizeof(xreftab));

  /* Renumber the elements */
  
  id = 1;
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
    elmt->workspace.ei[0] = elmt->id;
    elmt->id = id++;
  }

  /* Periodic BCs */
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
    for (edge = elmt->edge_list; edge; edge = edge->next) {
      if (edge->bc && strchr("PXY", BC_TYPE(edge->bc))) {
	double pos[2];
	int key;
	
	if (BC_TYPE(edge->bc) == 'X') {
	  pos[0] = -100000.;
	  pos[1] = edge->pos[1];
	} else {
	  pos[0] = edge->pos[0];
	  pos[1] = -100000.;
	}

	key = vdbdjoin(Geometry.mesh->edges,pos);
	elmt->workspace.ei[1+EDGE_ID(edge)] = key;

	if (xreftab[key].from.elmt == 0) {
	  xreftab[key].from.elmt = elmt->id;
	  xreftab[key].from.edge = edge->id;
	} else {
	  xreftab[key].to.elmt = elmt->id;
	  xreftab[key].to.edge = edge->id;
	}
      }
    }
  }

  /* Open the file (no parameters) */
  
  fp = fopen(name, "w");
  if (fp) 
    printf ("Saving current session to %s [SEMTEX]\n", name);
  else {
    printf ("Unable to open file %s for writing\n", name);
    return;
  }


  /* TOKENS */

  if (Geometry.param->count > 0) {
    param_t *p = Geometry.param;
    const int count = PARAM_COUNT(p);
    int i;
    fprintf (fp, "\n<TOKENS>\n");
    for (i = 0; i < count; i++)
      fprintf (fp, "   %s = %s\n", param_name(p,i), param_expr(p,i));
    fprintf (fp, "</TOKENS>\n");
  }

  /* USER */

  if (Geometry.user) {
    keyword_t *u = Geometry.user;
    const int count = KEYWORD_COUNT(u);
    int i;

    fprintf (fp, "\n<USER>\n");
    for (i = 0; i < count; i++)
      fprintf (fp, "   %s\n", KEYWORD_INFO(u,i));
    fprintf (fp, "</USER>\n");
  }

  /* FIELDS */

  if (Geometry.fields) {
    keyword_t *f = Geometry.fields;
    const int count = KEYWORD_COUNT(f);
    int i;

    fprintf (fp, "\n<FIELDS>\n    ");
    for (i = 0; i < count; i++)
      fprintf (fp, " %s", KEYWORD_INFO(f,i));
    fprintf (fp, "\n</FIELDS>\n");
  }

  /* HISTORY POINTS */

  if (Geometry.history) {
    keyword_t *his = Geometry.history;
    const int npts = KEYWORD_COUNT(his);
    int i;

    fprintf (fp, "\n<HISTORY NUMBER=%d>\n", npts);
    for (i = 0; i < npts; i++)
      fprintf (fp, "   %d %s\n", i+1, strchr(KEYWORD_INFO(his,i), ' '));
    fprintf (fp, "</HISTORY>\n");
  }


  /* NODES */

  if (vdbnmember(Geometry.mesh->vertx,NULL)>0) {
    const int nnodes = vdbnmember(Geometry.mesh->vertx, NULL);
    int i;

    struct {
      double x;
      double y;
    } nodetab[4*_MAX_NEL];

    Element *elmt;
    Vertex  *vert;
    
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
      for (vert = elmt->vert_list; vert; vert = vert->next) {
	const int key  = VERTEX_KEY(vert);
	nodetab[key].x = VERTEX_POS(vert,0);
	nodetab[key].y = VERTEX_POS(vert,1);
      }
    }

    fprintf (fp, "\n<NODES NUMBER=%d>\n", nnodes);
    for (i = 0; i < nnodes; i++)
      fprintf (fp, "   %d %.16g %.16g 0.\n", i+1, nodetab[i].x, nodetab[i].y);
    fprintf (fp, "</NODES>\n");
  }

  if (MESH_NELMT(Geometry.mesh)>0) {
    const int nel = MESH_NELMT(Geometry.mesh);

    Element *elmt;
    Vertex  *vert;

    fprintf (fp, "\n<ELEMENTS NUMBER=%d>\n", nel);
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
      fprintf (fp, "   %d <Q>", ELEMENT_ID(elmt));
      for (vert = elmt->vert_list; vert; vert = vert->next) {
	fprintf (fp, " %d", VERTEX_KEY(vert)+1);
      }
      fprintf (fp, " </Q>\n");
    }
    fprintf (fp, "</ELEMENTS>\n");
  }

  /* CURVES */

  if (count_curves(Geometry.mesh->head)) {
    const int ncurves = count_curves(Geometry.mesh->head);

    Element *elmt;
    Edge    *edge;
    int i = 1;

    fprintf (fp, "\n<CURVES NUMBER=%d>\n", ncurves);

    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
      for (edge = elmt->edge_list; edge; edge = edge->next) {
	if (edge->curve) {
	  const Curve *curve = edge->curve;

	  fprintf (fp, "   %d %d %d ", 
		   i++, ELEMENT_ID(elmt), EDGE_ID(edge)+1);

	  switch (curve->type) {

	  case Arc:
	    fprintf (fp, "<ARC> %.16g </ARC>\n", 
		     curve->info.arc.theta.range > 0. ?
		     curve->info.arc.radius : -curve->info.arc.radius);
	    break;

	  case Spline: {
	    double xp, yp;
	    spline_node(&curve->info.spline, 1, &xp, &yp);
	    fprintf (fp, "<SPLINE> %.16g %.16g", xp, yp);
	    spline_node(&curve->info.spline, 2, &xp, &yp);
	    fprintf (fp, " %.16g %.16g </SPLINE>\n", xp, yp);
	    break;
	  }

	  case File:
	    fprintf (fp, "<FILE> %s </FILE>\n", curve->info.file.name);
	    break;

	  default:
	    break;
	  }
	}
      }
    }
    fprintf (fp, "</CURVES>\n");
  }

  /* GROUPS + BCs + SURFACES */

  {
    struct {
      int  count;
      char *info[BC_MAXINFO];
    } bc_xreftab[52];

    int i;
    int ngroups   = 0;
    int nbcs      = 0;
    int nsurfaces = 0;

    /* clear table */
    memset(bc_xreftab, 0, sizeof(bc_xreftab));

    /* set table entries */
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
      for (edge = elmt->edge_list; edge; edge = edge->next) {
	if (edge->bc) {
	  const BC *bc = edge->bc;
	  const int n  = BC_TYPE(bc) - 'a';
	  if (strchr("XY", BC_TYPE(bc))==NULL) {
	    nsurfaces++;
	    if (bc_xreftab[n].count==0) {
	      bc_xreftab[n].count = BC_ITEMS(bc);
	      for (i = 0; i < BC_ITEMS(bc); i++)
		bc_xreftab[n].info[i] = BC_INFO(bc,i).expr;
	    }
	  } else {
	    int key = elmt->workspace.ei[1+EDGE_ID(edge)];
	    if (xreftab[key].to.elmt > ELEMENT_ID(elmt) ||
		(xreftab[key].to.elmt == ELEMENT_ID(elmt) &&
		 xreftab[key].to.edge > EDGE_ID(edge)))
	      nsurfaces++;
	  }
	}
      }
    }

    /* write groups */

    for (i = 0; i < 52; i++) {
      if (bc_xreftab[i].count>0)
	ngroups++;
    }
    nbcs = ngroups;

    fprintf (fp, "\n<GROUPS NUMBER=%d>\n", ngroups);

    ngroups = 1;
    for (i = 0; i < 52; i++) {
      if (bc_xreftab[i].count > 0) {
	const char c = 'a'+i;
	fprintf (fp, "   %d %c ", ngroups++, c);
	switch (c) {
	case 'v':
	  fprintf (fp, "velocity\n");
	  break;
	case 'w':
	  fprintf (fp, "wall\n");
	  break;
	case 'a':
	  fprintf (fp, "axis\n");
	  break;
	case 'o':
	  fprintf (fp, "outflow\n");
	  break;
	case 's':
	  fprintf (fp, "symmetry\n");
	  break;
	default:
	  fprintf (fp, "unknown\n");
	  break;
	}
      }
    }

    fprintf (fp, "</GROUPS>\n");

    /* write BCs */
    
    fprintf (fp, "\n<BCS NUMBER=%d>\n", nbcs);

    nbcs = 1;
    for (i = 0; i < 52; i++) {
      if (bc_xreftab[i].count>0) {
	const char group = 'a'+i;
	int j;
	fprintf (fp, "   %d %c %d\n", nbcs++, group, bc_xreftab[i].count);
	for (j = 0; j < bc_xreftab[i].count; j++) {
	  char *p;
	  char type;

	  strcpy(buf, bc_xreftab[i].info[j]);
	  if (p = strchr(buf,':')) {
	    type = buf[0];
	    fprintf (fp, "      <%c> %s </%c>\n", type, ++p, type);
	  } else {
	    fprintf (fp, "      <D> %c%d = %s </D>\n", group, j, buf);
	  }
	}
      }
    }

    fprintf (fp, "</BCS>\n");

    /* write SURFACES */

    fprintf (fp, "\n<SURFACES NUMBER=%d>\n", nsurfaces);

    nsurfaces = 1;
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
      for (edge = elmt->edge_list; edge; edge = edge->next) {
	if (edge->bc) {
	  const BC *bc = edge->bc;
	  if (strchr("XY", BC_TYPE(bc))==NULL) {
	    fprintf (fp, "   %d %d %d <B> %c </B>\n", 
		     nsurfaces++, ELEMENT_ID(elmt), EDGE_ID(edge)+1, BC_TYPE(bc));
	  } else {
	    int key = elmt->workspace.ei[1+EDGE_ID(edge)];
	    if (xreftab[key].to.elmt > ELEMENT_ID(elmt) ||
		(xreftab[key].to.elmt == ELEMENT_ID(elmt) &&
		 xreftab[key].to.edge > EDGE_ID(edge))) {
	      fprintf (fp, "   %d %d %d <P> %d %d </P>\n", 
		       nsurfaces++, ELEMENT_ID(elmt), EDGE_ID(edge)+1, 
		       xreftab[key].to.elmt, xreftab[key].to.edge+1);
	    }
	  }
	}
      }
    }
    fprintf (fp, "</SURFACES>\n");
  }

  /* Restore ID numbers */
  
  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
    ELEMENT_ID(elmt) = elmt->workspace.ei[0];


  fclose(fp);
}


void DoMeshOut() 
{
  char *name;

  if (Domain_require()) {
    
    if (!(name=strtok(NULL, " ;\n"))) {
      char c;
      printf ("Overwrite %s (y/n) ? ", Geometry.name);
      scanf("%1s", &c);
      getchar();
      if (c != 'y') {
	return;
      } else {
	name = Geometry.name;
      }
    }

    if (strstr(name,".rea"))
      MeshOut_prism_1(name);
    else if (strstr(name,".feml"))
      MeshOut_semtex(name);
    else if (option("format")==0)
      MeshOut_prism_1(name);
    else if (option("format")==1)
      MeshOut_semtex(name);
    else
      fprintf (stderr, "I'm confused -- what format was that?\n");
  }
}

void DoSave() {
  DoMeshOut();
}


void DoLock()
{
  if (Domain_require()) {
    int     n;
    double  *pos;
    Element *elmt;
    Vertex  *vert;
    
    n   = vdbnmember (Geometry.mesh->vertx, NULL);;
    pos = (double*) calloc (2*n, sizeof(double));
    
    /* Write */
    
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
      for (vert = elmt->vert_list; vert; vert = vert->next)
	memcpy (pos + 2*vert->key, vert->pos, 2*sizeof(double));
    
    /* Read */
    
    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next)
      for (vert = elmt->vert_list; vert; vert = vert->next)
	memcpy (vert->pos, pos + 2*vert->key, 2*sizeof(double));
    
    fprintf (stderr, "Vertex coordinates now agree exactly.\n"
	     "Please save and re-load this mesh.\n");
    
    free (pos);
  }
}

void DoMatrix ()
{
  if (Domain_require()) {
    Matrix *A = Geometry.A;

    if (A == NULL) A = Geometry.A = Matrix_alloc();
    
    Matrix_update (A, Geometry.mesh);
  }
}

static double rms (const Field *u)
{
  const double nrm2 = Field_nrm2(u);
  const int    npts = Field_npts(u);
  return sqrt(nrm2*nrm2/npts);
}

void DoSolve ()
{
  char *p = strtok(NULL,";\n");
  char  utype, ftype;

  static char *usage = "usage: solve <utype> <ftype>\n";

  if (!p) {
    fputs(usage,stderr);
    return;
  }

  while (*p && isspace(*p)) p++;
  utype = *p++;
  while (*p && isspace(*p)) p++;
  ftype = *p;

  if (!utype || !ftype) {
    fputs(usage,stderr);
    return;
  }

  if (Domain_require()) {
    Matrix *A    = Geometry.A;
    Field  *F    = Domain_getField(ftype);
    Field  *U    = Domain_chkField(utype);
    double t0, t1;
    
    if (!F) return;
    if (!U) U = Domain_addField(utype);
    if (!A) A = Geometry.A = Matrix_alloc();

    Field_realloc(U);
    Field_realloc(F);
    
    Matrix_update (A, Geometry.mesh);
    
    if (option("timer")) {
      t0 = dclock();
      Solve (U, F, Geometry.mesh, Geometry.A);
      t1 = dclock();
      printf ("solution time = %g msec\n", (t1-t0)*1000.);
    } else {
      Solve (U, F, Geometry.mesh, Geometry.A);
    }

    clear_command_line();
  }
}

void DoSolution ()
{
  if (Domain_require()) {
    char *p = strtok(NULL, ";\n");
    if (p) {
      if (Geometry.solution) keyword_free (Geometry.solution);
      Geometry.solution = keyword_alloc("SOLUTION");
      keyword_add (Geometry.solution, p);
    } else {
      printf ("solution = ");
      if (Geometry.solution) 
	puts (Geometry.solution->info[0]);
      else
	puts ("undefined");
    }
  }
}

void DoForce ()
{
  if (Domain_require()) {
    char *p = strtok(NULL,"<>");

    if (!p)
      fputs ("usage: force <f0, f1, ...>\n", stderr);
    else {
      char *info[16];
      char buf[BUFSIZ];

      int nitems = token_split(16, p, info);
      int i;

      if (Geometry.force) keyword_free(Geometry.force);

      Geometry.force = keyword_alloc("DRIVE FORCE");
      
      for (i = 0; i < nitems; i++) {
	sprintf (buf, "\tf%d = %s", i, strip(info[i]));
	keyword_add(Geometry.force, buf);
      }
    }
  }
}

void DoFields ()
{
  if (Domain_require()) {
    char *p = strtok(NULL,"<>");
    
    if (!p)
      fputs ("usage: fields <char, char, ...>\n", stderr);
    else {
      char *info[16];
      int n = token_split(16, p, info);
      int i;

      if (Geometry.fields) keyword_free(Geometry.fields);

      Geometry.fields = keyword_alloc("FIELDS");

      for (i = 0; i < n; i++) {
	keyword_add(Geometry.fields, strip(info[i]));
      }
    }
  }
}

void DoUser ()
{
  if (Domain_require()) {
    char *p = strtok(NULL,"<>");

    if (!p) {
      fputs ("usage: user < string >\n", stderr);
    } else {
      if (Geometry.user==NULL) 
	Geometry.user = keyword_alloc("USER");
      
      keyword_add(Geometry.user, strip(p));
    }
  }
}

void DoShow() 
{
  char *p = strtok(NULL, ";\n");
  char type[2];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%1s", type) != 1)
    error = 1;
  
  if (error) {
    fprintf (stderr, 
	     "usage: show <field>\n"
	     "Displays all grid-point values in <field>\n");
    return;
  }
  
  if (Domain_require()) {
    Field *u = Domain_getField(*type);
    
    if (!u)
      return;

    Field_show(u);
  }
}

void DoList()
{
  static char *help =
    "usage: list <type>\n"
    "type:\n"  
    "bc        Boundary conditions\n"
    "force     Drive forces\n"
    "history   History points\n"
    "icond     Initial conditions\n"
    "nodes     Node dictionary\n"
    "options   Options\n"
    "param     Parameters\n";

  if (Domain_require()) {
    char *p = strtok(NULL," \n");

    if (!p) {
      fputs (help, stderr);
    } else {
      const char c = *p;
      switch (c) {
      case 'b':
	DoBinfo();
	break;
      case 'f':
	if (Geometry.force)
	  keyword_export (Geometry.force, stdout);
	else 
	  printf ("list: no drive force has been specified\n");
	break;
      case 'i':
	if (Geometry.ic)
	  keyword_export (Geometry.ic, stdout);
	else
	  printf ("list: no initial conditions have been specified\n");
	break;
      case 'n':
	if (nodes) {
	  fputs ("Node dictionary:\n", stdout);
	  tree_walk(nodes);
	}
	break;
      case 'o':
	printf("Note: on=1, off=0");
	show_options(stdout);
	break;
      case 'p':
	if (Geometry.param) {
	  param_t *p = Geometry.param;
	  int i, n = PARAM_COUNT(p);
	  printf ("Parameters:\n");
	  for (i = 0; i < n; i++)
	    printf ("%s = %s\n", param_name(p,i), param_expr(p,i));
	} else
	  printf ("list: no parameters have been specified\n");
	break;
      case 'h':
	if (Geometry.history) 
	  keyword_export (Geometry.history, stdout);
	else
	  printf ("list: no history points have been specified\n");
	break;
      default:
	fprintf (stderr, "list: unknown item -- %s\n", p);
	break;
      }
    }
  }
}

static char warranty[] = "";

void DoWarranty() {
  fputs (warranty, stdout);
  fputs ("See the GNU General Public License for further details.\n", stdout);
}

static char copyright[] = "";

void DoCopyright() {
  fputs (copyright, stdout);
  fputs ("See the GNU General Public License for further details.\n", stdout);
}

/* ------------------------------------------------------------------------- */

#define DSTACKNAME 64
#define DSTACKSIZE  4
#define DSTACKBODY 32

typedef struct {
  double var;
  double end;
  double inc;
  char   name[DSTACKNAME];
  char  *body[DSTACKBODY];
} DSTACK;

static DSTACK *dstack;
static int   n_dstack;
static int     dindex = -1;

static char*  do_name;
static double do_var, do_end, do_inc;

void push_dstack (char *name, double start, double end, double inc)
{
  if (n_dstack==0) {
    n_dstack += DSTACKSIZE;
    if((dstack = (DSTACK*) calloc(n_dstack,sizeof(DSTACK))) == NULL)
      fputs("can't alloc do stack!",stderr);
  } else if (dindex >= 0) {
    dstack[dindex].var = do_var;
  }

  if(++dindex >= n_dstack) {
    n_dstack += DSTACKSIZE;
    if(!(dstack = (DSTACK*)realloc((void*)dstack, n_dstack*sizeof(DSTACK))))
      fputs("can't realloc do stack!", stderr);
  }

  strcpy(dstack[dindex].name,name);
  do_name = dstack[dindex].name;
  dstack[dindex].var = do_var = start;
  dstack[dindex].end = do_end = end;
  dstack[dindex].inc = do_inc = inc;
}  

void pop_dstack()
{
  if (--dindex >= 0) {
    do_name = dstack[dindex].name;
    do_var  = dstack[dindex].var;
    do_end  = dstack[dindex].end;
    do_inc  = dstack[dindex].inc;
  }
}

int next_do()
{
  if (do_var/do_inc > do_end/do_inc + 1.e-4) {
    pop_dstack();
    return -1;
  } else {
    dparam_set(do_name,do_var);
    do_var += do_inc;
    return 0;
  }
}

void DoLoop()
{
  char *p = strtok(NULL,"");
  char *info[3];
  double start, end, inc;

  char name[DSTACKNAME];
  char buf [BUFSIZ];
  int  n = 0;
  int  count = 0;

  extern FILE *mscope_command_stream;

  if (!p) {
    fputs ("usage: do var = start, end [,inc] { ... }\n",
	   stderr);
    clear_command_line();
  }

  strcpy(name, strip(strtok(p, " =")));

  while ((p=strtok(NULL,"= \t,{\n")) != NULL && count < 3)
    info[count++] = strip(p);
  start = atof(info[0]);
  end   = atof(info[1]);
  inc   = (count==3) ? atof(info[2]) : 1.;

  push_dstack(name, start, end, inc);

  /* Read the body of the do loop */

  do {
    fgets(buf, BUFSIZ, mscope_command_stream);
    dstack[dindex].body[n++] = strdup(buf);
  } while (!strchr(buf,'}'));

  free(dstack[dindex].body[--n]);
  dstack[dindex].body[n] = NULL;

  while (next_do() != -1) {
    for (n = 0; dstack[dindex].body[n] != NULL; n++) {
      strcpy (buf, dstack[dindex].body[n]);
      DoParse(buf);
    }
  }
}
  

/* ------------------------------------------------------------------------- *
 *                 P R I V A T E     F U N C T I O N S                       *
 * ------------------------------------------------------------------------- */

int scopy (int n, float *x, int incx, float *y, int incy) {
  int i;
  for (i = 0; i < n; i++, x += incx, y += incy)
    *y = *x;
  return 0;
}

static void set_levels (int n, double min, double max)
{
  double *levels = (double*) malloc(n*sizeof(double));
  double  du     = (max-min)/(n-1.);
  dramp (n, &min, &du, levels, 1);
  pl_levels(n, levels);
  free(levels);
}

void DoContour()
{
  char *p = strtok(NULL,";\n");

  if (!p)
    fputs("usage: contour type [min max nlevs]\n", stderr);
  else {
    char    type;
    int     nlevs;
    double  min, max;
    Element *elmt;
    Field   *u;

    if (sscanf(p, "%1s%lf%lf%d", &type, &min, &max, &nlevs) != 4)
      nlevs = 0;

    if (!(u = Domain_getField(type)))
      return;

    if (nlevs==0) {
      nlevs = 15;
      min   = Field_min(u);
      max   = Field_max(u);
    }

    set_levels (nlevs, min, max);
    fprintf(stderr,"contours of %c from %g to %g\n", type, min, max);

    for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
      int     np = elmt->nr;
      double *xp = elmt->xmesh[0];
      double *yp = elmt->ymesh[0];
      double *up = u->data[elmt->id];
      
      pl_contour (np, np, xp, yp, up);
    }
    pl_gflush();
  }
}

/*
 * Streamlines
 */


#define STREAM_MAX 1024

typedef struct {
  int    npts;
  double *x;
  double *y;
} stream_t;

static struct {
  int nstream;
  stream_t *s[STREAM_MAX];
} stream_list;

static stream_t *stream_alloc (double x, double y, double ds, 
			       int npts, int dir)
{
  stream_t *s = (stream_t*) calloc(1,sizeof(stream_t));

  Field *u = Domain_getField('u');
  Field *v = Domain_getField('v');
  Probe *p = Probe_alloc(u, PROBE_XP, x, y);

  double coef[4], uk[4], vk[4], ddir;
  int i;

  if (u==NULL || v==NULL || p==NULL) {
    fprintf(stderr, "stream: initial probe allocation failed\n");
    free(s);
    return NULL;
  }

  ddir = dir > 0 ? 1. : -1.;

  s->npts = npts;
  s->x    = dvector(0,npts-1);
  s->y    = dvector(0,npts-1);
  s->x[0] = x;
  s->y[0] = y;

  /* Use RK4 to integrate */

  coef[0] = ds/6.;
  coef[1] = ds/3.;
  coef[2] = ds/3.;
  coef[3] = ds/6.;

  pl_relocate(x,y);

  for (i = 1; i < npts; i++) {

    Probe_move(p,x,y);
    uk[0] = Probe_eval(p,u);
    vk[0] = Probe_eval(p,v);

    Probe_move(p, x + 0.5*ds*uk[0], y + 0.5*ds*vk[0]);
    uk[1] = Probe_eval(p,u);
    vk[1] = Probe_eval(p,v);

    Probe_move(p, x + 0.5*ds*uk[1], y + 0.5*ds*vk[1]);
    uk[2] = Probe_eval(p,u);
    vk[2] = Probe_eval(p,v);

    Probe_move(p, x + ds*uk[2], y + ds*vk[2]);
    uk[3] = Probe_eval(p,u);
    vk[3] = Probe_eval(p,v);

    x += ddir*ddot(4, coef, 1, uk, 1);
    y += ddir*ddot(4, coef, 1, vk, 1);

    if (Probe_move(p,x,y) == -1) {
      fprintf (stderr, "stream: terminated early at (%g,%g)\n", x, y);
      break;
    }

    s->x[i] = x;
    s->y[i] = y;

    pl_draw(x,y);
  }

  s->npts = i;
  Probe_free(p);
  return s;
}

static void stream_free (stream_t *s) {
  free(s->x);
  free(s->y);
  free(s);
}
  
static void stream_draw (const stream_t *s) 
{
  const int npts = s->npts;
  const double *x = s->x;
  const double *y = s->y;
  int i;

  pl_relocate(x[0], y[0]);
  for (i = 1; i < npts; i++)
    pl_draw(x[i], y[i]);
  pl_gflush  ();
}

static void stream_add (stream_t *s)
{
  if (stream_list.nstream < STREAM_MAX)
    stream_list.s[stream_list.nstream++] = s;
  else
    fprintf (stderr, "stream: buffer is full -- not added\n");
}

void DoStream()
{
  char *p = strtok(NULL,";\n");
  static char *usage = 
    "usage: stream [draw|list|clear|add x y ds npts dir]\n";

  if (!p)
    fputs(usage, stderr);
  else {
    int i, npts=256, dir=1;
    double x0, y0;
    double ds=0.1;
    char command[16];

    if (sscanf(p, "%s%lf%lf%lf%d%d", command, &x0, &y0, &ds, &npts, &dir)==0){
      fputs(usage,stderr);
    } else {
      switch (tolower(command[0])) {
      case 'a': {
	stream_t *s;

	if (dir != 0) {
	  s = stream_alloc(x0,y0,ds,npts,dir);
	  stream_draw(s);
	  stream_add (s);
	} else {
	  s = stream_alloc(x0,y0,ds,npts,1);
	  stream_draw(s);
	  stream_add (s);
	  s = stream_alloc(x0,y0,ds,npts,-1);
	  stream_draw(s);
	  stream_add (s);
	}

	break;
      }
      case 'c':
	for (i = 0; i < stream_list.nstream; i++)
	  stream_free(stream_list.s[i]);
	stream_list.nstream = 0;
	break;
      case 'd':
	for (i = 0; i < stream_list.nstream; i++)
	  stream_draw(stream_list.s[i]);
	break;
      case 'l':
	for (i = 0; i < stream_list.nstream; i++)
	  printf ("%4d %4d %#10.6g %#10.6g\n", i, 
		  stream_list.s[i]->npts,
		  stream_list.s[i]->x[0], 
		  stream_list.s[i]->y[0]);
	break;
      default:
	fputs(usage,stderr);
	return;
      }
    }
  }
}

/* ------------------------------------------------------------------------- */

void DoVectors()
{
  char *p = strtok(NULL,";\n");

  Field *u = Domain_getField('u');
  Field *v = Domain_getField('v');
  double scal;
  Element *elmt;

  static char *usage = "vectors <scale>";

  if (!u || !v) {
    fprintf (stderr, "usage: %s\n", usage);
    return;
  }
      
  if (p) {
    scal = atof(p);
  } else {
    double umax = Field_amax(u);
    double vmax = Field_amax(v);
    scal = scalar("0.1/sqrt((XMAX-XMIN)^2+(YMAX-YMIN)^2)") / 
      MAX(umax, vmax);
  }

  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    double *xp = elmt->xmesh[0];
    double *yp = elmt->ymesh[0];
    double *up = ELEMENT_DATA(elmt,u);
    double *vp = ELEMENT_DATA(elmt,v);
    double  np = elmt->nr;

    pl_vectors (np, np, xp, yp, scal, up, vp);
  }
  pl_gflush();
}

/* ------------------------------------------------------------------------- */

void DoProfile()
{
  char *p = strtok(NULL,";\n");
  static char *usage = "profile type x0 y0 dx dy [npts]";

  if (!p)
    fprintf(stderr, "usage: %s\n", usage);
  else {
    double x0, y0;
    double dx, dy;
    char type;
    int  npts;
    int  i;

    Field  *u;
    Probe  *probe;
    double *xp;
    double *yp;
    double *sp;
    double *up;

    const char *delim = " \t";
    char *ptr, *endp = p;

    /* set "flag" values */

    dx   = -1.;
    dy   = -1.;
    npts = 0;

    if (ptr = get_word(endp,&endp,delim)) type = *ptr;
    if (ptr = get_word(endp,&endp,delim)) x0   = scalar(ptr);
    if (ptr = get_word(endp,&endp,delim)) y0   = scalar(ptr);
    if (ptr = get_word(endp,&endp,delim)) dx   = scalar(ptr);
    if (ptr = get_word(endp,&endp,delim)) dy   = scalar(ptr);
    if (ptr = get_word(endp,&endp,delim)) npts = atoi  (ptr);

    if (dx < 0. || dy < 0.) {
      fprintf (stderr, "profile: bad input string\n");
      return;
    }
    if (npts == 0) npts = 64;

    if (!(u = Domain_getField(type))) {
      fprintf (stderr, "profile: type %c is undefined\n", type);
      return;
    }

    if (!(probe = Probe_alloc(u, PROBE_XP, x0, y0))) {
      fprintf (stderr, "profile: unable to allocate a probe\n");
      return;
    }

    xp  = (double*) malloc(npts*sizeof(double));
    yp  = (double*) malloc(npts*sizeof(double));
    sp  = (double*) malloc(npts*sizeof(double));
    up  = (double*) malloc(npts*sizeof(double));
    dx /= (npts-1.);
    dy /= (npts-1.);
    for (i = 0; i < npts; i++) {
      xp[i] = x0 + i*dx;
      yp[i] = y0 + i*dy;
      sp[i] = hypot(xp[i]-x0, yp[i]-y0);
      if (Probe_move(probe, xp[i], yp[i])==-1) {
	printf ("warning: point (%g,%g) not found\n", xp[i], yp[i]);
	up[i] = 0.;
      } else
	up[i] = Probe_eval(probe, u);
    }

    /* Compute display limits */

    if (option("autoProfile")) {
      const double smid = (sp[npts-1] + sp[0]) * 0.5;
      const double slen = (sp[npts-1] - sp[0]) * 1.05 / 2.;
      double umin, umax;

      char buf[8];
      sprintf (buf, "%c(s)", type);

      umin = up[idmin(npts,up,1)];
      umax = up[idmax(npts,up,1)];

      pl_erase  ();
      pl_limits (smid-slen, smid+slen, 
		 umin-.05*(umax-umin), umax+.05*(umax-umin));
      pl_box    ();
      pl_xlabel ("s");
      pl_ylabel (buf);
    }

    pl_connect(npts, sp, up);
    
    free(xp);
    free(yp);
    free(up);
    free(sp);
    Probe_free(probe);
  }

  pl_gflush();
}

/* ------------------------------------------------------------------------- */

void DoNorms()
{
  char *p = strtok(NULL, ";\n");

  if (!p) {
    fprintf (stderr, "usage: norm <type>\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*p);
    if (u) {
      printf ("norm(%c) = [%g %g %g]\n",
	      FIELD_TYPE(u),
	      Field_amax(u),
	      Field_L2  (u),
	      Field_H1  (u));
    }
  }
}

void DoEnorms()
{
  char *p = strtok(NULL, ";\n");
  static char *usage = "usage: enorm <type> <expr>\n";

  if (!p) {
    fputs(usage,stderr);
    return;
  }

  if (Domain_require()) {
    char type;
    char expr[BUFSIZ];
    Field *u;

    if (sscanf(p,"%1s%s", &type, expr) != 2) {
      fputs(usage,stderr);
      return;
    }

    if ((u = Domain_getField(type))) {
      Field *v = Field_dup(u);

      Field_set (v,expr);
      Field_axpy(-1., u, v);

      printf ("Errors in field '%c': [%g %g %g]\n", 
	      FIELD_TYPE(u), 
	      Field_amax(v),
	      Field_L2  (v),
	      Field_H1  (v));

      Field_free(v);
    }
    clear_command_line();
  }
}

void DoSet()
{
  char *p = strtok(NULL,";\n");
  char *type;
  char *expr;

  if (!p) {
    fprintf (stderr, "usage: set <type> = <expr>\n");
    return;
  }

  for (type = p; *type && isspace(*type); type++) 
    continue;
  expr = strchr(p,'=')+1;

  if ((!type) || (!expr)) {
    fprintf (stderr, "usage: set <type> = <expr>\n");
    return;
  }

  if (Domain_require()) {
    char type = *p;
    Field *u  = Domain_chkField(type);

    if (!u)
      u = Domain_addField(type);

    Field_set (u, expr);
  }
}

void DoFscal()
{
  char *p = strtok(NULL,";\n");
  char type[2];
  char expr[BUFSIZ];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%s%1s", expr, type) != 2)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fscal <scalar> <field>\n"
	     "computes <field> *= <scalar>\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type);
    double d = scalar(expr);

    if (!u) return;

    Field_scal(d,u);
  }
}

void DoFint()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%1s", type_u) != 1)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fint <field>\n"
	     "computes the integral of <field> over the domain\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);

    if (!u) return;

    printf("Integral[%c] = %g\n", FIELD_TYPE(u), Field_integral(u));
  }
}

void DoFabs()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%1s%1s", type_u, type_v) != 2)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fabs <field-1> <field-2>\n"
	     "computes <field-2> = abs(<field-1>)\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_chkField(*type_v);

    if (!u) 
      return;
    if (!v)
      v = Domain_addField(*type_v);

    Field_abs(u,v);
  }
}

void DoFaxpy()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  char expr[BUFSIZ];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%s%1s%1s", expr, type_u, type_v) != 3)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: faxpy <scalar> <field-1> <field-2>\n"
	     "computes <field-2> += <scalar> * <field-1>\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_getField(*type_v);
    double d = scalar(expr);

    if (!u || !v) return;

    Field_axpy(d,u,v);
  }
}

void DoFcopy()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%1s%1s", type_u, type_v) != 2)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fcopy <field-1> <field-2>\n"
	     "computes <field-2> = <field-1>\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_chkField(*type_v);

    if (!u) 
      return;
    if (!v)
      v = Domain_addField(*type_v);

    Field_copy(u,v);
  }
}

void DoFadd()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  char type_w[2];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%1s%1s%1s", type_u, type_v, type_w) != 3)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fadd <field-1> <field-2> <field-3>\n"
	     "computes <field-3> = <field-1> + <field-2>\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_getField(*type_v);
    Field *w = Domain_chkField(*type_w);

    if (!u || !v) 
      return;
    if (!w) 
      w = Domain_addField(*type_w);

    Field_add(u,v,w);
  }
}

void DoFmult()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  char type_w[2];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%1s%1s%1s", type_u, type_v, type_w) != 3)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fmult <field-1> <field-2> <field-3>\n"
	     "computes <field-3> = <field-1> * <field-2>\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_getField(*type_v);
    Field *w = Domain_chkField(*type_w);

    if (!u || !v) 
      return;
    if (!w) 
      w = Domain_addField(*type_w);

    Field_mult(u,v,w);
  }
}

void DoFdiv()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  char type_w[2];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%1s%1s%1s", type_u, type_v, type_w) != 3)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fdiv <field-1> <field-2> <field-3>\n"
	     "computes <field-3> = <field-1> / <field-2>\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_getField(*type_v);
    Field *w = Domain_chkField(*type_w);

    if (!u || !v) 
      return;
    if (!w) 
      w = Domain_addField(*type_w);

    Field_div(u,v,w);
  }
}

void DoFgrad()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  int  dir;
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%i%1s%1s", &dir, type_u, type_v) != 3)
    error = 1;
  else if (dir < 0 || dir > 2)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fgrad <dir> <field-1> <field-2>\n"
	     "computes <field-2> = D[<field-1>, x_<dir>], dir=[0..2]\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_chkField(*type_v);

    if (!u)
      return;
    if (!v) 
      v = Domain_addField(*type_v);

    Field_gradient(u,v,dir);
  }
}

void DoFFT()
{
  char *p = strtok(NULL,";\n");
  char type_u[2];
  char type_v[2];
  char expr[BUFSIZ];
  int  error = 0;

  if (!p)
    error = 1;
  else if (sscanf(p, "%s%1s%1s", expr, type_u, type_v) != 3)
    error = 1;

  if (error) {
    fprintf (stderr, 
	     "usage: fft <dir> <field-1> <field-2>\n"
	     "computes <field-2> = FourierTransform[<field-1>], where\n"
	     "dir = -1 (Fourier) or dir = +1 (Physical)\n");
    return;
  }

  if (Domain_require()) {
    Field *u = Domain_getField(*type_u);
    Field *v = Domain_chkField(*type_v);
    int  dir = scalar(expr);

    if (!u) 
      return;
    if (!v)
      v = Domain_addField(*type_v);

    Field_FFT(u,v,dir);
  }
}

/* Plot the combined spectra of all elements */

void DoSplot()
{
  char *p = strtok(NULL,";\n");

  if (Domain_require()) {
    char type = *p;
    Field *u  = Domain_getField(type);
    Element *elmt = FIELD_HEAD(u);

    const int np = MESH_NR(Geometry.mesh);
    double **a = dmatrix(0,np-1,0,np-1);
    int i, j;

    while (elmt) {
      legcoef(np,np,ELEMENT_DATA(elmt,u),*a);

	for (j = 0; j < np;j++) {
      for (i = 0; i < np; i++) {
	  double n = sqrt(i*i + j*j);
	  pl_relocate (n, log(fabs(a[i][j])));
	  pl_dot();
	}
      }
      elmt = elmt->next;
    }

    pl_gflush ();
  }
}
  
