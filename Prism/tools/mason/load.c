/*
 * Functions for loading the input file
 *
 * RCS Information
 * ------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ------------------------------------------------------------------------ */


#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "mason.h"

typedef struct {      /* ..............  BC Info .............. */
  char  type    ;     /* The character type flag                */
  int   element ;     /* Element number                         */  
  int   face    ;     /* Face number                            */
  union {             /* Boundary data                          */
    float val   ;     /*    floating-point value                */
    char* str   ;     /*    function string                     */
  } data [8]    ;     /*                                        */
} bcInfo;

static bcInfo loadBC      (FILE *fp);
static void   findSection (FILE *fp, char *secname);

/* ------------------------------------------------------------------------ *
 * loadParams() -- Load the session parameters                              *
 *                                                                          *
 * This function's only real purpose is to scan through the parameter list  *
 * and set the N-order (unless set from the command line).                  *
 * ------------------------------------------------------------------------ */

void loadParams (Domain *omega, FILE *fp)
{
  int  n;
  char buf  [BUFSIZ];
  char name [32], value[32];

  for (n = 0; n < 4; n++)
    fgets (buf, BUFSIZ, fp);

  if (sscanf (buf, "%d", &n) != 1)
    error_msg ("can't read the # of parameters");

  /* Scan the parameter list and set the N-order (if not set already) */

  while (n--) {
    fgets (buf, BUFSIZ, fp);
    if (sscanf (buf, "%25s%25s", value, name) == 2)
      if (!omega->norder && 
	  (strstr (name, "NORDER") || strstr (name, "MODES")))
	omega->norder = atoi (value);
  }

  return;
}

/* ------------------------------------------------------------------------ *
 * loadMesh() -- Load the vertex coordinates                                *
 *                                                                          *
 * This function load the coordinates of the mesh.  If the -tri option is   *
 * in effect, it only load the first three points as the corners of the     *
 * triangle.  Otherwise, it loads all four for a rectilinear element.       *
 *                                                                          *
 * Performs memory allocation for the slave vertices.                       *
 * ------------------------------------------------------------------------ */
 
void loadMesh (Domain *omega, FILE *fp)
{
  char  buf    [BUFSIZ];
  Point corner [MAX_VERTX];
  int   nel, i, k;
  Element *U;

  findSection (fp, "MESH");

  if (sscanf (fgets(buf, BUFSIZ, fp), "%d", &nel) != 1)
    error_msg ("unable to read the number of elements");
  
  /* Set up a new element vector */

  omega->elements = nel;
  omega->U =  U   = (Element *) calloc (nel, sizeof(Element));

  /* Read in the element coordinates.  The following section reads *
   * coordiantes one at a time and assigns the element a pointer   *
   * to a Vertex based on the geometry.  In the process, we const- *
   * ruct the list of master vertices.                             */

  for (k = 0; k < nel; k++) {

    fgets (buf, BUFSIZ, fp);  /* Element header */

    U[k].id       =  k + 1;
    U[k].nr       =  omega->norder;
    U[k].ns       =  omega->norder;
    U[k].vertices = (omega->tri) ? 3 : 4;
    U[k].vlist    = (VertexP *) calloc (U[k].vertices, sizeof(VertexP));

    /* Read the coordinates of the corners */

    for (i = 0; i < U[k].vertices; i++)
      fscanf (fp, "%lf", & corner[i].x);
    fgets (buf, BUFSIZ, fp);
    for (i = 0; i < U[k].vertices; i++)
      fscanf (fp, "%lf", & corner[i].y);
    fgets (buf, BUFSIZ, fp);


    /* Do the vertex assignments.  If a match isn't found *
     * then a new vertex is created and inserted at the   *
     * top of the list.                                   */

    for (i = 0; i < U[k].vertices; i++)
      if (!(U[k].vlist[i] = findVertex (omega->master, corner[i])))
	U[k].vlist[i] = omega->master = makeVertex (omega->master, corner[i]);
    
    U[k].next = &U[k] + 1;
  }      
  
  U[--k].next = (Element *) NULL;

  return;
}

/* ------------------------------------------------------------------------ *
 * loadBCs() -- Load boundary conditions                                    *
 *                                                                          *
 * This function loads the boundary conditions and allocates memory for the *
 * edges.  Boundary conditions are translated into one of three possible    *
 * states: Fixed (physical boundaries), Joined (element-element boundaries),*
 * or Patched (master-slave boundaries).                                    *
 *                                                                          *
 * NOTE: Triangular meshes must still define four boundary conditions per   *
 *       element, even though only the first three are used.                *
 * ------------------------------------------------------------------------ */

void loadBCs (Domain *omega, FILE *fp)
{
  Element *U;
  bcInfo   bc;
  char     buf [BUFSIZ];
  int f;

  findSection (fp, "BOUNDARY");
  do
    fgets (buf, BUFSIZ, fp);
  while
    (strstr (buf, "NO"));

  
  for (U = omega->U; U ; U = U->next) {

    U->edges = (omega->tri) ? 3 : 4;
    U->elist = (Edge *) calloc (U->edges, sizeof(Edge));

    for (f = 0; f < (*U).edges; f++) {
      
      /* Initialze the edge */

      U -> elist[f].id  = f + 1;
      U -> elist[f].iel = U->id;
      U -> elist[f].np  = f & 1 ? (*U).nr-2 : (*U).ns-2;


      /* Load the boundary condition for this edge and verify */

      bc = loadBC (fp);

      if (bc.element != U->id || bc.face != (f+1)) {
	sprintf (buf, "Mismatched element/side -- got (%d,%d), "
		      "expected (%d,%d)\n", bc.element, bc.face, U->id, f+1);
	error_msg(buf);
      }


      /* Condense the boundary conditions to { [M|S], D, [E|P] } */

      switch (toupper(bc.type)) {
      case 'M': case 'S': {
	int    patch   = (int) bc.data[0].val;
	int    segment = (int) bc.data[1].val;
	char   type    = toupper(bc.type);
	Patch *P;

	if (!(P = findPatch (patch, omega)))
	  omega->P = P = makePatch (patch, omega->P);
	
	U -> elist[f].type         = type;
	U -> elist[f].bc.p.patch   = patch;
	U -> elist[f].bc.p.segment = segment;
	   
	if (type == 'M')
	  P->masters = 
	    makeSegment(segment, P->masters, U->elist+f);
	else
	  P->slaves  = 
	    makeSegment(segment, P->slaves , U->elist+f);
	break;
      }

      case 'D': case 'N': 
      case 'T': case 'V': case 'W': 
      case 'F': case 'O':
      case 'A':
	U -> elist[f].type         = 'D';
	U -> elist[f].bc.s         =  bc.data[0].val;
	break;

      case 'E': case 'P':
	U -> elist[f].type         =  bc.type;
	U -> elist[f].bc.c.element = (int) bc.data[0].val;
	U -> elist[f].bc.c.face    = (int) bc.data[1].val;
	break;

      default:
	fprintf 
	  (stderr, "mason: unknown boundary condition -- %c\n", bc.type);

	U -> elist[f].type         = 'D';
	U -> elist[f].bc.s         =  0.;
	break;
      }
    }

    if (omega->tri)      /* For triangular elements, load the dummy */
      bc = loadBC (fp);  /* boundary condition.                     */
  }

  return;
}


/* ------------------------------------------------------------------------ *
 *                    P R I V A T E     F U N C T I O N S                   *
 * ------------------------------------------------------------------------ */

/* Search through a file for a section keyword */

static void findSection (FILE *fp, char *secname)
{
  char buf [BUFSIZ];

  while (fgets (buf, BUFSIZ, fp))
    if (strstr (buf, secname))
      break;

  return;
}

/* Load a boundary condition line */

static bcInfo loadBC (FILE *fp)
{
  int    n = 0;
  char   buf [BUFSIZ], *p;
  bcInfo new;
  fpos_t pos;

  p = fgets (buf, BUFSIZ, fp);
  while (isspace(*p))
    p++;

  
  /* If the boundary condition is given,  then use it.  Otherwise, 
  // flag it with a '?' and let the patch numbering section take 
  // care of it.                                                     */

  if (isalpha(*p))
    new.type = *(p++);
  else
    new.type = '?';

  /* Read the rest of the BC info */

  sscanf (p, "%d%d%f%f%f", & new.element, & new.face, 
	  & new.data[0].val, & new.data[1].val, & new.data[2].val);
  
  if (islower(new.type)) {
    do {
      fgetpos (fp, &pos);
      fgets   (buf, BUFSIZ, fp);
      if ((p=strchr(buf,'=')))
	new.data[n].str = strdup (++p);
    } while
      ( p );

    fsetpos (fp, &pos);
  }
  
  return new;
}
  
