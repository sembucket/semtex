/*                                                   
 * Spectral Element Data Model
 *
 * Based (loosely) on NEKTON's view of finite element discretizations
 *                                                                          
 * This file contains functions to translate a spectral element input file 
 * in the NEKTON format into the set of data structures used by Cubit.
 * The file format is fairly specific, expecting the following sections in 
 * order:                                             
 *                                                                          
 *     PARAMETERS          . Symbol table parameters                        
 *     MESH                . Element geometry                               
 *     CURVES              . Element boundary curves, part of the mesh
 *     BOUNDARY CONDITIONS . Boundary conditions for all fields             
 *                                                                          
 * Section names in the file are always compared by ignoring case, so you can
 * use upper, lower, or mixed-case section names in the file.   The functions 
 * in this file process each section to define the following parameters:
 *
 *
 *     ReadParams   ---  Load symbols into the parser's symbol table
 *     ReadMesh     ---  Create a Field (array of elements) with the given
 *                       geometry
 *     ReadBCs      ---  Create an array of BCs (boundary edges) that
 *                       define the boundary conditions for a Field
 *
 * Notes:                                                                   
 *                                                                          
 *   - Parameters may be either integers or reals.  Fortran naming conven-  
 *     tions are used to decide for UPPERCASE parameters, but any lower-   
 *     case parameter is assumed to be real.  There are execeptions (see    
 *     below) for certain standard variable names.                          
 *                                                                          
 *   - Boundary conditions are a bit tricky, so please read the comments    
 *     in ReadBCs if you aren't sure how call this function.           
 *                                                                          
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

#include "veclib/veclib.h"

#include "cubit.h"
#include "element.h"
#include "edge.h"
#include "vertex.h"
#include "bc.h"
#include "curve.h"
#include "isomesh.h"
#include "session.h"

/* Private functions */

static int rescale   (Element *);
static int checklist (char *name, char *list[]);
static int bedgcmp   (const void *, const void *);

/* .....  Section names ..... */
#define  SEC_PARAMS   sections[0]             /* Parameters                 */
#define  SEC_MESH     sections[1]             /* Mesh definition            */
#define  SEC_BCS      sections[2]             /* Boundary conditions        */
#define  SEC_ICS      sections[3]             /* Initial conditions         */
#define  SEC_DF       sections[4]             /* Drive force                */
#define  SEC_HIST     sections[5]             /* History points             */

static char *sections[] = { 
  "PARAMETERS", 
  "MESH", 
  "BOUNDARY CONDITIONS", 
  "INITIAL CONDITIONS", 
  "FORCE", 
  "HISTORY" 
};

void ReadParams (FILE *fp)
{
  int    n;
  char   buf[BUFSIZ], value[25], name[25], c;
  double val;

  static char *dspecial[] = { "LAMBDA", "LZ", 0 };
  static char *ispecial[] = { "EQTYPE", "TORDER", 0 };

  rewind (fp);
  if (!findSection (SEC_PARAMS, buf, fp))
    cubit_err ("no \"parameters\" section was found");
  
  /*
  // SEC_PARAMS begins with the program version, a dimension statement, 
  // and then the number of parameters.  Skip the first two and leave
  // the last one stored in buf.
  */

  for (n = 0; n < 3; n++) fgets(buf, BUFSIZ, fp);

  if (sscanf(buf, "%d", &n) != 1) 
    cubit_err ("can't read the # of parameters");
  
  while (n--) {
    fgets (buf, BUFSIZ, fp);
    if(sscanf(buf, "%25s%25s", value, name) == 2) {
      if (checklist (name, dspecial))
	dparam_set (name, parse(value));
      else if (checklist (name, ispecial))
	iparam_set (name, parse(value));
      else if (isupper(c = *name) && 'I' <= c && c <= 'N')
	iparam_set (name, parse(value));
      else
	dparam_set (name, parse(value));
    }
  }

  /* Running at a different NORDER ? */

  if ((n = iparam ("NORDER-req")) && n != UNSET)
           iparam_set("NORDER", n);
  
  iparam_set ("NR", iparam("NORDER"));
  iparam_set ("NS", iparam("NORDER"));

  return;
}

static int checklist (char *name, char *list[])
{
  do
    if (strcmp(name,*list) == 0)
      return 1;
  while 
    (*++list);

  return 0;
}

static Element **elmt_array_alloc (int nr, int ns, int nz, int nel)
{
  Element **list = (Element**) malloc(nel*sizeof(Element*));
  int k;

  for (k = 0; k < nel; k++)
    list[k] = Element_alloc(nr,ns,nz);

  for (k = 0; k < nel-1; k++)
    list[k]->next = list[k+1];

  return list;
}

/* Read the mesh description */

Element *ReadMesh (FILE *fp)
{
  double v;
  Element **list, *head;
  char buf[BUFSIZ];
  int  nr, ns, nz, nel, i, k, status;
  
  status = findSection(SEC_MESH, buf, fp) != NULL; 
  assert (status);
  status = sscanf(fgets(buf,BUFSIZ,fp), "%d", &nel);
  assert (status == 1);

  iparam_set ("ELEMENTS", nel);
  
  /* Set up a new array of elements */

  nr   = iparam("NORDER");         
  ns   = iparam("NORDER");
  nz   = iparam("NZ");
  list = elmt_array_alloc (nr, ns, nz, nel);

  /* Read in mesh information
   *
   * Initially we assume that every element has straight sides.  This
   * can be overridden when curved-side specifications are loaded after
   * the basic element specs have been read.                            
   */

  for (k = 0; k < nel; k++) {
    Element *elmt = list[k];

    /* Read the header and try to decode [family,key] info */

    fgets(buf, BUFSIZ, fp);  
    if (sscanf(buf,"%*s%*d%*s%d%*s%d", &elmt->family, &elmt->key) != 2)
      elmt->family = -1;

    /* Read the (x,y) coordinates for the this element */

    for (i = 0; i < 4; i++)
      fscanf (fp, "%lf", (*elmt->xmesh) + elmt->edge_list[i].start);
    for (i = 0; i < 4; i++)
      fscanf (fp, "%lf", (*elmt->ymesh) + elmt->edge_list[i].start);
  
    fgets (buf, BUFSIZ, fp);
  }

  /* Read curved side declarations */

  fgets (buf, BUFSIZ, fp);
  fgets (buf, BUFSIZ, fp);
  sscanf(buf, "%d", &k);   /* k = number of curved sides */

  if (k) {
    int     iside, iel;
    char    type, *p;
    Curve   *c;
    Edge    *edge;
    Element *elmt;

    while (k--) {
      fgets  (buf, BUFSIZ, fp);
      status = sscanf(buf, "%1d%d", &iside, &iel);
      assert (status == 2);
      
      iel--; iside--;
      p     = buf + strlen(buf); while (isspace(*--p));
      type  = toupper(*p);
      elmt  = list[iel];
      edge  = elmt->edge_list + iside;
      c = edge->curve = (Curve*) calloc (1, sizeof(Curve));

      switch (type) {

      case 'L':                     /* Line */
	break;
	
      case 'C': {                   /* Circular arc */
	double radius;
	sscanf (buf, "%*1d%*d%lf", &radius);
	c->type = Arc;
	c->info.arc = make_arc (elmt, edge, radius);
	break;
      }

      case 'S':  {                  /* Bezier cubic spline */
	double pos[4];
	sscanf (buf, "%*1d%*d%lf%lf%lf%lf", pos, pos+1, pos+2, pos+3);
	c->type = Spline;
	c->info.spline = spline_alloc (elmt, edge, pos);
	break;
      }

      case 'F': {                   /* Geometry File */
	char name[FILENAME_MAX];
	sscanf (buf, "%*d%*d%s", name);
	c->type = File;
	c->info.file = geofile_alloc(name);
	/* Since this type of curve actually modifies the geometry of the *
	 * element, the boundary points need to generated immediately.    *
	 * Otherwise, straight sides generated nearby can get screwed up. */
	curve(elmt,edge);
	break;
      }

      default: 
	sprintf (cubit_err_msg, "unknown curve type -- %c", type);
	cubit_err (NULL);
	break;
      }
    }
  }

  /* Now the geometric description for every edge in the mesh is stored *
   * in the corresponding Curve, so we just need to generate the mesh   *
   * and mapping factors for each element.                              */

  for (k = 0; k < nel; k++) {
    shape   (list[k]);                   /* Boundaries        */
    blend   (list[k]);                   /* Interior          */
    rescale (list[k]);                   /* Scale and shift   */
    map     (list[k]);                   /* Mapping factors   */
    normals (list[k]);                   /* Outward normals   */
  }
  
  head = list[0];
  free(list);
  return head;
}

static int rescale (Element *U)
{
  int    nrns = U->nr * U->ns;
  double sx   = dparam("XSCALE"),
         sy   = dparam("YSCALE"),
         xo   = dparam("XSHIFT"),
         yo   = dparam("YSHIFT");
  double *x   = U->xmesh[0],
         *y   = U->ymesh[0];
  
  if (sx == 1. && sy == 1. && xo == 0. && yo == 0.) return 0;

  dscal (nrns, sx, x, 1)      ;  dscal (nrns, sy, y, 1);
  dsadd (nrns, xo, x, 1, x, 1);  dsadd (nrns, yo, y, 1, y, 1);

  return 1;
}

/* Restore the original key values or set defaults */

void ReadKeys (FILE *fp, Element *list)
{
  Element *elmt = list;

  char buf[BUFSIZ];
  int  status;

  if(findSection("KEYS", buf, fp)) {
    fgets (buf, BUFSIZ, fp);
    sscanf(buf, "%d", &status);
    while (elmt && fgets(buf, BUFSIZ, fp)) {
      sscanf (buf, "%d%d", &elmt->family, &elmt->key);
      elmt = elmt->next;
    }
  } else
    while (elmt) {
      elmt->family = elmt->id;
      elmt->key    = 1;
      elmt         = elmt->next;
    }
}

/* ------------------------------------------------------------------------ *
 * ReadBCs() - Read Boundary Conditions                                     *
 *                                                                          *
 * This function reads boundary information for a spectral element field    *
 * and sets up the appropriate boundary structures.  The following boundary *
 * conditions are currently recognized:                                     *
 *                                                                          *
 *   Flag      Type         Field #1       Field #2       Field #3          *
 *  -----    --------      ----------     ----------     ----------         *
 *    V      Velocity      U-velocity     V-Velocity     W-Velocity         *
 *    W      Wall          = 0            = 0            = 0                *
 *    v      v             U(x,y,z)       V(x,y,z)       W(x,y,z)           *
 *    A      Axis          dU/dn = 0      dV/dn = 0      dW/dn = 0          *
 *    F      Flux          U' = f1        = 0            = 0                *
 *    f      flux          U'(x,y,z)      V'(x,y,z)      W'(x,y,z)          *
 *    O      Outflow       U'.n = 0       V'.n = 0       W'.n = 0           *
 *    E      Element       w/ ID          w/ EDGE                           *
 *    P      Periodic      w/ ID          w/ EDGE                           *
 *    D      Dirichlet     patch number   segment number                    *
 *    N      Neumann       patch number   segment number                    *
 *    M      Master        patch number   segment number                    *
 *    S      Slave         patch number   segment number                    *
 *                                                                          *
 * NOTES:                                                                   *
 *                                                                          *
 *   - The argument "fnum" specifes which Field Number to read, and it      *
 *     must be less than DIM.  If fnum = 1, then ReadBCs will scan through  *
 *     the file until it finds another section of boundary conditions.      *
 *     Otherwise, it reads back through the current section and selects     *
 *     boundary conditions from the specified field.                        *
 *                                                                          *
 *   - For boundaries with a given function of (x,y,z),  the boundary       *
 *     conditions are read from the lines following rather than the         *
 *     columns.  If ReadBCs finds another boundary condition before it      *
 *     reads "fnum" lines, it exits with an error.                          *
 *                                                                          *
 *     Here is an example (element 1, edge 1):                              *
 *                                                                          *
 *            1  1  v   0.0     0.0     0.0                                 *
 *                 ux = (1 - y^2)*sin(t)                                    *
 *                 uy = (1 - y^2)*cos(t)                                    *
 *                 uz = 0.                                                  *
 *            1  2  E   2.0     1.0     0.0                                 *
 *                                                                          *
 *     The "=" MUST be present, and the function is defined by the string   *
 *     to the right of it.                                                  *
 *                                                                          *
 *   - NO TIME DEPENDENT BOUNDARY CONDITIONS.                               *
 *                                                                          *
 * ------------------------------------------------------------------------ */

BC *ReadBCs (FILE *fp, int group, Element *U)
{
  int nbcs   = 0;
  BC *bclist = (BC*) NULL, *bc;

  char type, buf[BUFSIZ];
  Element *elmt = U;

  findSection ("BOUNDARY CONDITIONS", buf, fp);
  fgets(buf, BUFSIZ, fp);

  while (elmt) {
    Edge *edge;
    for  (edge = elmt->edge_list; edge; edge = edge->next) {
      char *p = fgets (buf, BUFSIZ, fp);
      int iel;
      int iedge;
      
      /* Mark the limits of the string */

      if (p == (char*) NULL) exit(-1);
      buf [strcspn(p,"#\n")] = '\0';

      p     = strtok(buf," \t");
      type  = *p;
      p     = strtok(NULL," \t");
      iel   = atoi(p);
      p     = strtok(NULL," \t");
      iedge = atoi(p);

      if ((iel - elmt->id) != 1 || (iedge - edge->id) != 1) {
	printf ("mismatched element/side: loop = %d %d, file = %d %d\n",
		elmt->id+1, edge->id+1, iel, iedge);
	exit (-1);
      }

      switch (type) {

	/* Dummy side types */
	
      case 'N': case 'M': case 'S': case 'E':
	break;

	
	/* Allocate a new boundary condition */
      
      default:

	bc       = (BC*) calloc (1, sizeof(BC));
	bc->type = type;
	bc->next = bclist;   /* Push it onto the list */
	bclist   = bc;
	edge->bc = bc;       /* Link to the edge */

	if (type == 'P') {   /* Periodic boundaries are special */
	  struct {
	    Element *elmt;
	    Edge    *edge;
	  } to;

	  int jel   = atoi(strtok(NULL," \t"))-1;
	  int jedge = atoi(strtok(NULL," \t"))-1;
	  
	  for (to.elmt = U; to.elmt; to.elmt = to.elmt->next)
	    if (to.elmt->id == jel)
	      break;
	  
	  if (to.elmt) {
	    for (to.edge=to.elmt->edge_list; to.edge; to.edge = to.edge->next)
	      if (to.edge->id == jedge)
		break;
	  }

	  if (to.elmt == NULL || to.edge == NULL) {
	    cubit_err ("bad periodic boundary condition");
	  }

	  if (fabs(edge->pos[0] - to.edge->pos[0]) < 1.e-6)
	    bc->type = 'Y';
	  else if (fabs(edge->pos[1] - to.edge->pos[1]) < 1.e-6)
	    bc->type = 'X';
	  else
	    cubit_err ("bad periodic boundary condition");

	  bc->nitems = 2;
	  bc->info[0].number = jel;
	  bc->info[1].number = jedge;

	} else if (isupper(type))   /* Read data values from the buffer */
	  while ((p = strtok(NULL," \t")) && (bc->nitems < BC_MAXINFO))
	    bc->info[bc->nitems++].value = atof(p);
	else {
	  fpos_t mark;

	  do {
	    fgetpos(fp, &mark);
	    if (p = fgets(buf, BUFSIZ, fp)) {
	      if (p = strchr(buf, '=')) {
		while (isspace(*++p));
		p[strlen(p)-1] = '\0';
		bc->info[bc->nitems++].expr = strdup(p);
	      } else
		p = (char*) NULL;
	    }
	  } while (p && (bc->nitems < BC_MAXINFO));

	  fsetpos (fp, &mark);
	  break;
	}
      }
      while (p = strtok(NULL, " "));
    }
    elmt = elmt->next;
  }
  return bclist;
}

/* Find a section header in the input file */

char *findSection (const char *name, char *buf, FILE *fp)
{
  int i;
  char *p;

  while (p = fgets(buf, BUFSIZ, fp)) {
    for (i = 0; *p && i < BUFSIZ; i++, p++)
      buf[i] = toupper(*p);
    if (strstr (buf, name)) 
      break;
  }
  
  return p;
}
