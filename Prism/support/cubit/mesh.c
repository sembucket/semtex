/*
 * Dynamic Mesh Generation
 *
 * Copyright (c) 1995 R. D. Henderson and Caltech
 *
 * The following functions maintain an adaptive discretization of a two-
 * dimensional region of space.  The space is broken up into quadrilateral
 * regions called `elements'.  An element is defined by four vertices and
 * four edges.  The edges can be straight, arcs, or Bezier splines.
 *
 * The mesh is stored as a collection of quadtrees (multiple trees are 
 * required to represent complex geometries).  The terminating nodes of
 * the tree represent ACTIVE elements.  Other (itermediate) nodes represent
 * IDLE elements.  Anytime the mesh is updated --- by refinement or 
 * coarsening --- a pointer is set to a particular element at the head of a 
 * linked list through all of the active elements in the mesh.
 *
 * Elements are identified by three numbers: an id number, a family, and a
 * key.  The id number is simply an integer that identifies the element as
 * part of a sequence.  A family identifies elements that are derived from
 * a common root.  A key is a binary number that identifies an element, its
 * siblings, and its ancestors in the mesh.  The root of the tree is called
 * `level zero' and the `level' of successive generations increases.
 *
 * --------------------------------------------------------------------------
 *
 * The library works by using the following operations:
 *
 *      install                           Install a key
 *
 *      prune                             Prune the tree (coarsen)
 *
 *      refine                            Bisect an element
 *
 *      restrict                          Refine coarse neighbors
 *
 * Connectivity in the mesh is maintained by using a Voxel Data Base for 
 * the positions of the vertices and edges.  Every position in the data-
 * base is assigned a unique index that can be used to translate a position
 * to an array index.  Indices can be optimized to group nearest-neighbors
 * together, for example to minimize the bandwidth of a finite element 
 * matrix.
 *
 * Voxel Data Base software based on the original ideas of Roy Williams
 *
 * RCS Information
 * --------------------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include <stdio.h>

#include "vdb/vdb.h"
#include "veclib/veclib.h"

#include "cubit.h"
#include "edge.h"
#include "vertex.h"
#include "alias.h"
#include "element.h"
#include "bc.h"
#include "curve.h"
#include "surface.h"

/* ------------------------------------------------------------------------- */

#define IDLE          0     /* Status flags    */
#define ACTIVE        1

#define REFINE_FLAG(elmt) \
    ((elmt)->workspace.ei[0])
#define STATUS_FLAG(elmt) \
    ((elmt)->workspace.ei[1])

/* User call-back functions */

int (*user_refine)() = NULL;      /* Refine an element    */
int (*user_bc)()     = NULL;      /* Get BC information   */
int (*user_prune)()  = NULL;      /* Prune an element     */
int (*user_perm)()   = NULL;      /* Permute two elements */

/* Private Functions ( i_* = internal recurvsive ) */

static int edge_update   (Element *child[], Edge *edge);

static int keycheck  (int key);
static int level     (int key);

static int lookup    (Mesh *mesh, int family, int key);
static int refine    (Mesh *mesh, int family, int key);
static int unrefine  (Mesh *mesh, int family, int key);

static int i_prune   (Mesh *mesh, int family, int key, int maxlevel);

/* ------------------------------------------------------------------------- *
 * Mesh_alloc, Mesh_free                                                     *
 *                                                                           *
 * These are the functions that initially create and eventually destroy a    *
 * mesh.                                                                     *
 * ------------------------------------------------------------------------- */

Mesh *Mesh_alloc()
{
  Mesh *mesh = (Mesh*) calloc(1,sizeof(Mesh));
  const int bufsize = 1 << 15; /* 32 Kb buffer */

  mesh->elements = 0;
  mesh->active   = 0;
  mesh->nr       = iparam("NORDER");
  mesh->ns       = iparam("NORDER");
  mesh->nz       = iparam("NZ");

  mesh->vertx    = NULL;
  mesh->edges    = NULL;
  mesh->elmts    = vdbcreate(2, 0.1, 0, 0, _MAX_NEL, malloc(bufsize), bufsize);

  mesh->bc       = NULL;
  mesh->head     = NULL;

  return mesh;
}

void Mesh_free (Mesh *mesh)
{
  int k;

  if (mesh->vertx) vdbdestroy (mesh->vertx);
  if (mesh->edges) vdbdestroy (mesh->edges);
  if (mesh->elmts) vdbdestroy (mesh->elmts);

  /* free Elements */
  /* free BCs */

  free (mesh);
}

/* ------------------------------------------------------------------------- *
 * Mesh_add() -- Add a single element to the mesh                            *
 *                                                                           *
 * This is primarily a start-up function to generate the initial list of     *
 * elements.                                                                 *
 * ------------------------------------------------------------------------- */

void Mesh_add (Mesh *mesh, Element *elmt)
{
  float pos[2];

  if (elmt->family == -1) {
    elmt->family = mesh->elements;
    elmt->key    = 1;
  }

  pos[0]   = elmt->family;
  pos[1]   = elmt->key;
  elmt->id = vdbfjoin (mesh->elmts, pos);

  mesh->list[elmt->id] = elmt;
  mesh->active = ++mesh->elements;

  /* Sweep to the end of the active list and append this element */

  if (mesh->head) {
    Element *last = mesh->head;
    while (last->next) 
      last = last->next;
    last->next = elmt;
    elmt->prev = last;
    elmt->next = NULL;
  } else {
    mesh->head = elmt;
    elmt->prev = NULL;
    elmt->next = NULL;
  }
}

/* ------------------------------------------------------------------------- *
 * Mesh_install() -- Install a key                                           *
 *                                                                           *
 * This is used during the start-up phase to generate the fine mesh from a   *
 * coarse one plus a set of "keys" that activate a certain element at the    *
 * deepest level of the tree.  Keys are installed recursively to generate    *
 * all of the intermediate levels in the tree as well.                       *
 * ------------------------------------------------------------------------- */

static i_install (Mesh *mesh, int family, int key)
{
  int status = 0;

  if (key <= 0) {
    printf ("i_install: tried to install a nil-key\n");
    return status;
  }

  if (keycheck(key)) {
    printf ("i_install: %d is an invalid key\n", key);
    return status;
  }

  /* Is it there? */
  
  if ((status = lookup(mesh, family, key)) == -1) {
    if (i_install (mesh, family, key >> 2)) { /* No?  Install my parent */
          refine  (mesh, family, key >> 2);   /* Refine my parent...    */
	  status = 1;                         /*  ...which installs me! */
    }
  } 
  
  return 0 <= status;
}

void Mesh_install (Mesh *mesh, int family, int key) {
  i_install (mesh, family, key);
}

/* ------------------------------------------------------------------------- *
 * Mesh_refine() -- Refine the mesh                                          *
 *                                                                           *
 * This function loops through the mesh and refines all elements that have   *
 * the refinement flag set.  During mesh refinement, any nonconforming       *
 * neighbors that would violate the 2:1 ratio along an edge are flagged for  *
 * refinement in a second pass.  Mesh_refine() performs a loop until no      *
 * flagged elements remain.                                                  *
 * ------------------------------------------------------------------------- */

static int *check = NULL;

void Mesh_refine (Mesh *mesh)
{
  int nref = 1;
  Element *elmt;
  Edge    *edge;

  check = (int*) malloc (4*_MAX_NEL*sizeof(int));  /* buffer */
  
  while (nref > 0) {

    nref = 0;
    elmt = mesh->head;
    memset (check, '\0', 4*_MAX_NEL*sizeof(int));
  
    while (elmt) {
      if (REFINE_FLAG(elmt)) {
	Element *tmp = elmt->next;
	refine (mesh, elmt->family, elmt->key);
	elmt = tmp;
	nref++;
      } else
	elmt = elmt->next;
    }

    printf ("Refine: pass 1 [%d set]\n", nref);

    vdbsync    (mesh->vertx);
    vdbsync    (mesh->edges);
    vdbcombine (mesh->edges, check, vdbimax, sizeof(int), 1, sizeof(int));

    nref = 0;
    for (elmt = mesh->head; elmt; elmt = elmt->next) {
      for (edge = elmt->edge_list; edge; edge = edge->next)
	if (check[edge->key] == 1) {
	  Element_setRefineOn(elmt);
	  nref++;
	}
    }

    if (nref) 
      printf ("Refine: pass 2 [%d flags] -- looping\n", nref);
    else
      printf ("Refine: done\n");
  }

  free(check); check = NULL;
  Mesh_connect (mesh);
  return;
}

/* Set the "refine" flag for an element */

int Element_setRefineOn (Element *elmt) {
  REFINE_FLAG(elmt) = 1;
  return 1;
}

void Mesh_resetFlags (Mesh *mesh)
{
  Element *elmt;

  for (elmt = mesh->head; elmt; elmt = elmt->next)
    memset (&(elmt->workspace), 0, sizeof(elmt->workspace.ef));
}

/* ------------------------------------------------------------------------- *
 * Mesh_restrict() -- Impose a smoothness criteria on the mesh               *
 *                                                                           *
 * This function performs iterative refinement until the mesh satisfies a    *
 * set of restrictions on allowable topology.  The following criteria are    *
 * checked:                                                                  *
 *                                                                           *
 * PERIODIC: Nonconforming edges are not allowed across periodic boundaries. *
 * This is enforced by triggering refinement in the periodic element when-   *
 * ever one of the pair is refined.                                          *
 *                                                                           *
 * LEVELS: Adjancent elements in the mesh must be no more than 1 level apart *
 * in the refinement tree.  This means that nonconforming edges can only     *
 * split two-to-one.                                                         *
 * ------------------------------------------------------------------------- */

static double* getpos (char type, const double pos[], double xpos[])
{
  const int p = type-'X';

  xpos[0] = pos[0];
  xpos[1] = pos[1];
  xpos[p] = FLT_MAX;   /* set the periodic direction to FLT_MAX */

  return xpos;
}

#define MAXENTRY (4*_MAX_NEL)

static int periodic (Mesh *mesh)
{
  VDB nodes;
  int nnodes;
  int nref = 0;

  const double tol = dparam("TOLVDB");
  const int bufsiz = 1 << 15;
  char *buf = (char*) malloc(bufsiz);

  int rflag[MAXENTRY];
  int nmult[MAXENTRY];
  double pos[2];

  Surface *s = Surface_alloc(mesh, "XYZ");

  /* Clear flags */

  memset(rflag, '\0', sizeof(rflag));

  /* Create periodic nodes database */

  nodes = vdbcreate(2, tol, 0, 0, 4*_MAX_NEL, buf, bufsiz);

  Surface_begin(s);
  while (!Surface_end(s)) {
    const Edge* edge = Surface_edge(s);
    const char  type = edge->bc->type;

    vdbdjoin(nodes, getpos(type, edge->pos,    pos));
    vdbdjoin(nodes, getpos(type, edge->a->pos, pos));
    vdbdjoin(nodes, getpos(type, edge->b->pos, pos));

    Surface_next(s);
  }

  nnodes = vdbsync(nodes);
  vdbnmember(nodes,nmult);

  /* Now check for bad periodic edges */

  Surface_begin(s);
  while (!Surface_end(s)) {
    const Edge* edge = Surface_edge(s);
    const char  type = edge->bc->type;
    const int   key  = vdbdpeek(nodes, getpos(type,edge->pos,pos));

    if (nmult[key]==1) {
      rflag[vdbdpeek(nodes,getpos(type,edge->a->pos,pos))]=1;
      rflag[vdbdpeek(nodes,getpos(type,edge->b->pos,pos))]=1;
    }

    Surface_next(s);
  }

  vdbcombine(nodes, rflag, vdbimax, sizeof(int), 1, sizeof(int));

  /* Scan back through and set refinement flags */

  Surface_begin(s);
  while (!Surface_end(s)) {
    const Edge* edge = Surface_edge(s);
    const char  type = edge->bc->type;
    const int   key  = vdbdpeek(nodes, getpos(type,edge->pos,pos));

    if (rflag[key] == 1) {
      Element_setRefineOn(Surface_elmt(s));
      nref++;
    }

    Surface_next(s);
  }

  if (nref) printf ("Restrict: %d periodic edge violation%c\n", 
		    nref, nref > 1 ? 's' : ' ');

  free(buf);
  vdbdestroy(nodes);
  Surface_free(s);

  return nref;
}

static int levels (Mesh *mesh) 
{
  int nvert, key, nref;
  int rflag[MAXENTRY];

  Element *elmt;
  Vertex  *vert;

  nref   = 0;
  nvert  = vdbnmember (mesh->vertx, NULL);
  memset(rflag, '\0', sizeof(rflag));

  for (elmt = mesh->head; elmt != NULL; elmt = elmt->next) {
    const int ilev = level(elmt->key);
    for (vert = elmt->vert_list; vert; vert = vert->next) {
      key        = vert->key;
      rflag[key] = MAX(rflag[key], ilev);
    }
  }

  vdbcombine (mesh->vertx, rflag, vdbimax, sizeof(int), nvert, 1);

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    const int ilev = level(elmt->key);
    for (vert = elmt->vert_list; vert; vert = vert->next) {
      if (abs(ilev - rflag[vert->key]) > 1) { 
	Element_setRefineOn(elmt); 
	nref++; 
      }
    }
  }

  if (nref) printf ("Restrict: %d level violation%c\n", nref, nref>1?'s':' ');

  return nref;
}

#undef MAXENTRY

void Mesh_restrict (Mesh *mesh)
{
  int nref = 0;

  nref += levels  (mesh);
  nref += periodic(mesh);

  /* ----------------- Possible Recursion ----------------- */

  if (nref) {
    printf ("Restrict: marked %d elements\n", nref);
    Mesh_refine   (mesh);
    Mesh_restrict (mesh);
  } else
    printf ("Restrict: all done\n");

  return;
}

/* ------------------------------------------------------------------------- *
 * Mesh_prune() -- Prune a family                                            *
 *                                                                           *
 * This is a high-level function to coarsen the mesh.  Rather than select    *
 * specific elements you select the maximum depth of the tree for a certain  *
 * family.  The refinement tree is then pruned at this level.                *
 * ------------------------------------------------------------------------- */

void Mesh_prune (Mesh *mesh, int family, int maxlevel)
{
  Element *elmt;
  int i;

  /* First mark the active elements */

  for (i = 0; (elmt = mesh->list[i]) != NULL; i++)
    STATUS_FLAG(elmt) = IDLE;
  for (elmt = mesh->head; elmt; elmt = elmt->next)
    STATUS_FLAG(elmt) = ACTIVE;

  i_prune (mesh, family, 1, maxlevel);

  return;
}

/* ------------------------------------------------------------------------- *
 * Mesh_lookup() -- Look up the ID of a [family,key] pair                    *
 * ------------------------------------------------------------------------- */

int Mesh_lookup (Mesh *mesh, int family, int key) 
{ return lookup (mesh, family, key); }

/* ------------------ P R I V A T E    F U N C T I O N S ------------------- */

static i_prune (Mesh *mesh, int family, int key, int maxlevel)
{
  const int me = lookup(mesh, family, key);
  int   status = IDLE;
  int k;

  if (me > -1) {   /* Do I exist? */

    /* Unless I'm ACTIVE, prune my children and get their group status */

    if (STATUS_FLAG(mesh->list[me]) != ACTIVE) {
      status = IDLE;
      for (k = 0; k < 4; k++) 
	status |= i_prune (mesh, family, (key << 2) + k, maxlevel);
    }
    
    /* Now we're guaranteed that the four children of (family,key) *
     * are the only decendents along this branch of the tree.  If  *
     * I'm IDLE and I have ACTIVE children at a higher level, then *
     * unrefine them.                                              */

    if (level(mesh->list[me]->key) >= maxlevel && status == ACTIVE)
      unrefine (mesh, family, key);

    /* I just gave myself the former status of my children.  If I  *
     * don't have any children, return MY status to my parent.     */

    else status = STATUS_FLAG(mesh->list[me]);
  }

  return status;
}

/* Unrefine the four children of (family, key) */

static unrefine (Mesh *mesh, int family, int key)
{
  Element *parent, *child[4];
  Edge    *edge;
  Vertex  *vert;
  int k, loc;

  loc    = lookup (mesh, family, key);
  parent = mesh->list[loc];
  
  /* Some checks */

  if (parent == NULL) {
    printf ("unrefine: no element exists at position %d\n", loc);
    exit   (-1);
  }
  else if (lookup(mesh, family, key << 2) < -1) {
    printf ("unrefine: element (%d,%d) has no children\n", family, key);
    exit   (-1);
  } else
    printf ("unrefine: re-activating element %d\n", loc);

  /* Locate the four children */

  for (k = 0; k < 4; k++) {
    if ((loc = lookup (mesh, family, (key << 2) + k)) < -1) {
      printf ("unrefine: lookup failed for child %d\n", k);
      exit   (-1);
    } else
      child[k] = mesh->list[loc];
  }

  /* Check active vs. idle status */

  for (k = 0; k < 4; k++) {
    if (STATUS_FLAG(child[k]) != ACTIVE) {
      printf ("unrefine: child %d is not active\n", k);
      exit(-1);
    } else
      STATUS_FLAG(child[k]) = IDLE;
  }

  STATUS_FLAG(parent) = ACTIVE;

  
  /* --------------------------------- */
  /* Map children's data to the parent */
  /* --------------------------------- */
  

  /* Unlink children from the mesh, relink the parent */

  for (k = 0; k < 4; k++) {
    for (edge = child[k]->edge_list; edge; edge = edge->next)
      vdbddelete (mesh->edges, edge->pos);
    for (vert = child[k]->vert_list; vert; vert = vert->next)
      vdbddelete (mesh->vertx, vert->pos);
  }

  if (parent->prev = child[0]->prev)
    parent->prev->next = parent;
  if (parent->next = child[3]->next)
    parent->next->prev = parent;


  /* Check global links */

  if (parent->prev == NULL) mesh->head = parent;

  return 0;
}

/* ------------------------------------------------------------------------- *
 * Note that "gamma_p" is the length of the parent edge, and "gamma_c" is    *
 * length of the child edge.  The order of the solution is "m" along the     *
 * parent edge and "n" along the child.                                      *
 *                                                                           *
 *                        | o    o     o     o    o |  child                 *
 *                                                                           *
 *   |-------  s0  ------>|=========================|  integration strip     *
 *       (i.e, s0 > 0)                                                       *
 *   |..........................................................| parent     *
 *                                                                           *
 * ------------------------------------------------------------------------- */


double **project (int m, int n, double s0, double gamma_p, double gamma_c)
{
  const int np = n+1;
  const int mp = m+1;

  double *z1 = dvector(0,n);
  double **p = dmatrix(0,m,0,n);

  double *z, *w, fac;
  double *zp;
  int i, j, k;

  /* Get collocation points & weights for the integration */

  getops (np, &z,  &w,   NULL, NULL);
  getops (mp, &zp, NULL, NULL, NULL);
  dzero  (mp*np, *p, 1);

  for (i = 0; i < np; i++) {
    z1[i] = (1 + z[i]) * gamma_c / gamma_p - 1. + 2. * s0/gamma_p;
  }

  printf ("zp, z, w = \n");
  for (i = 0; i < np; i++)
    printf ("%#10.4g %#10.4g %#10.4g\n", z1[i], z[i], w[i]);

  fac = gamma_c / 2.;
  for (i = 0; i <= m; i++) {
    for (j = 0; j <= n; j++) {
      p[i][j] = fac * w[j] * hgll(i, z1[j], zp, mp);
    }
  }

  printf ("projection matrix:\n");
  for (i = 0; i <= m; i++) {
    for (j = 0; j <= n; j++) {
      printf ("%#10.4g ", p[i][j]);
    }
    printf("\n");
  }

  return p;
}

/* Check to see if a key is valid */

static keycheck (int key)
{
  int status = 0;
  int k = key, p = 0;
    
  if (key < 1) 

    status = 1;  /* Key must be > 1 */

  else if (key > 1) {

    /* Key must not be an odd power of 2 (3, 9, 31, etc) */

    while (k >>= 1) p++;
    if ( (p & 1) && (1 << p) == key )   
      status = 1;
  }

  return status;
}

/* Check to see if <family, key> is in the forest */

static lookup (Mesh *mesh, int family, int key)
{
  float pos[2];

  pos[0] = family;
  pos[1] = key;

  return vdbfpeek(mesh->elmts, pos);
}

/* ------------------------------------------------------------------------- */

static refine (Mesh *mesh, int family, int key)
{
  float    pos[2];
  Element *child[4], *parent;
  Edge    *edge;
  Vertex  *vert;

  int i, k, loc;

  /* Lookup the (family, key) we're supposed to refine */

  parent = mesh->list[loc = lookup(mesh, family, key)];
  if (parent == NULL) {
    printf ("refine: no element exists at position %d\n", loc);
    exit (-1);
  }

  /* Flag any aliased edges */

  if (check) {
    for (edge = parent->edge_list; edge; edge = edge->next)
      if (edge->alias != NULL)
	if (edge->alias->type > 0)
	  check[edge->alias->key.edge] = 1;
  }

  /* Remove the parent's edges and vertices from the database */

  for (edge = parent->edge_list; edge; edge = edge->next)
    vdbddelete (mesh->edges, edge->pos);
  for (vert = parent->vert_list; vert; vert = vert->next)
    vdbddelete (mesh->vertx, vert->pos);


  /* See if the children need to be created or simply reactivated */

  if (lookup (mesh, family, key << 2) > -1) {
    for (k = 0; k < 4; k++) {
      if ((loc  = lookup(mesh, family, (key << 2) + k)) > -1)
	child[k] = mesh->list[loc];
      else {
	printf ("refine: failed on lookup of child %d\n", k);
	exit   (-1);
      }
    }

  } else {
    const int nr = ELEMENT_NR(parent);
    const int ns = ELEMENT_NS(parent);
    const int nz = ELEMENT_NZ(parent);

    for (k = 0; k < 4; k++) {
      pos[0] = family;
      pos[1] = (key << 2) + k;
      loc    = vdbfjoin (mesh->elmts, pos);
  
#if 0
      child[k] = elmt_alloc (nr, ns, loc);
#else
      child[k] = Element_alloc(nr,ns,nz);
#endif

      child[k]->id     =  loc;
      child[k]->emap   =  parent->emap;
      child[k]->family =  parent->family;
      child[k]->key    = (parent->key << 2) + k;
      
      elmt_project (parent, *parent->xmesh, *child[k]->xmesh, k);
      elmt_project (parent, *parent->ymesh, *child[k]->ymesh, k);

      mesh->list[loc] = child[k];
    }

    /* Update BC's and curved sides */

    for (edge = parent->edge_list; edge; edge = edge->next)
      edge_update (child, edge);

    /* Update the geometry and join the VDB database */
    
    for (k = 0; k < 4; k++) {
#if 0
      /* Generate exact boundaries. This turns off the propogation of *
       * curvilinear boundaries to children unless they lie along an  *
       * external boundary with a high-level curve description.  To   *
       * recover that behaviour just comment out this loop.           */
      int i;
      for (i = 0; i < 4; i++)
	curve (child[k], child[k]->edge_list+i);
#endif

      blend      (child[k]);   /* Blend the mesh                 */
      map        (child[k]);   /* Generate isoparametric mapping */
      normals    (child[k]);   /* Compute edge outward normals   */
      
      for (vert = child[k]->vert_list; vert; vert = vert->next)
	vert->key = vdbdjoin (mesh->vertx, vert->pos);
      for (edge = child[k]->edge_list; edge; edge = edge->next)
	edge->key = vdbdjoin (mesh->edges, edge->pos);
    }
  }

  /* Link the children into the mesh, unlink the parent */
  
  if (child[0]->prev = parent->prev)
    parent->prev->next = child[0];
  
  child[0]->next = child[1];
  child[1]->next = child[2];
  child[2]->next = child[3];
  
  child[1]->prev = child[0];      
  child[2]->prev = child[1];      
  child[3]->prev = child[2];
  
  if (child[3]->next = parent->next)
    parent->next->prev = child[3];

  /* Check global links */
  
  if (child[0]->prev == NULL) mesh->head = child[0];
  
  /* --------------------------------- */
  /* Map parent's data to the children */
  /* --------------------------------- */

  if (user_refine != NULL) (*user_refine)(parent, child);

  mesh->elements += 4;
  mesh->active   += 4 - 1;
  return 0;
}

/* Map some data... */

int elmt_project (Element *elmt, double *in, double *out, int quad)
{
  const int nr = elmt->nr;
  const int ns = elmt->ns;
  double *z, *low, *high, *rr, *ss, *hr, hs;
  int i, j, p, q;

  quad &= 0x3;               /* Only 2 bits, so you can pass a key */
  low   = dvector (0, nr);   /* New locations */
  high  = dvector (0, nr);
  hr    = dvector (0, nr);   /* Interpolants  */

  getops (nr, &z, 0, 0, 0);

  for (i = 0; i < nr; i++) {       
    high[i] = (1. + z[i]) * .5;    /* "high" maps (0,1) */
    low [i] = high[i] - 1.;        /* "low" maps (-1,0) */
  }
  
  rr = (quad == 0 || quad == 3) ? low : high;
  ss = (quad == 0 || quad == 1) ? low : high;

  dzero (nr*ns, out, 1);

  for (i = 0; i < ns; i++) {
    for (j = 0; j < nr; j++) {
      for (q = 0; q < nr; q++)               /* sum factorization */
	hr[q] = hgll (q, rr[j], z, nr);      
      for (p = 0; p < ns; p++) {
	hs         = hgll(p, ss[i], z, ns); 
	out[i*nr+j] += hs * ddot (nr, in + p*nr, 1, hr, 1);
      }
    }
  }

  free (low);
  free (high);
  free (hr);

  return 0;
}

/* ------------------------------------------------------------------------- */

void Mesh_info (Mesh *mesh)
{
  int alias  = 0;
  int patch  = 0;
  int depth  = 0;
  int active = 0;

  Element *elmt;
  Edge    *edge;

  for (elmt = mesh->head; elmt; elmt = elmt->next, active++) {
    int ilevel = Element_getLevel(elmt);
    if (ilevel > depth)
      depth = ilevel;
    for (edge = elmt->edge_list; edge; edge = edge->next)
      if (edge->alias) {
	if (edge->alias->type == 0)
	  patch++;
	else
	  alias++;
      }
  }

  printf ("Elements  : %d [active] %d [total]\n", active, mesh->elements);
  printf ("Resolution: %d x %d x %d\n", mesh->nr, mesh->ns, mesh->nz);
  printf ("Levels    : %d\n", depth+1);
  printf ("Patches   : %d\n", patch);
  printf ("Aliases   : %d\n", alias);
  printf ("Edges     : %d\n", vdbnmember(mesh->edges, NULL));
  printf ("Vertices  : %d\n", vdbnmember(mesh->vertx, NULL));

  return;
}

static edge_update (Element *child[], Edge *edge)
{
  Element *c1 = child[ edge->id];
  Element *c2 = child[(edge->id+1)%4];
  Edge    *e1 = c1->edge_list + edge->id;
  Edge    *e2 = c2->edge_list + edge->id;

  if (edge->bc) {
    e1->bc = edge->bc; 
    e2->bc = edge->bc;
  }

  if (edge->curve) {

    memcpy (e1->curve = (Curve*) malloc (sizeof(Curve)),
	    edge->curve, sizeof(Curve));
    memcpy (e2->curve = (Curve*) malloc (sizeof(Curve)),
	    edge->curve, sizeof(Curve));

    switch (edge->curve->type) {

    case Arc: {
      double dt = edge->curve->info.arc.theta.range;
      e1->curve->info.arc.theta.range  = .5 * dt;
      e2->curve->info.arc.theta.start += .5 * dt;
      e2->curve->info.arc.theta.range  = .5 * dt;
      break;
    }

    case Spline: {
      spline_t s1, s2;
      spline_cut (&edge->curve->info.spline, &s1, &s2);
      e1->curve->info.spline = s1;
      e2->curve->info.spline = s2;
      break;
    }

    default:
      break;
    }

    curve (c1, e1);
    curve (c2, e2);
  }

  return 0;
}

static int level (int key)
{
  int n = 0;
  while (key >>= 2) n++;
  return n;
}

/* ------------------------------------------------------------------------- */

int Element_getFamily (const Element *elmt)
{ return elmt->family; }

int Element_getKey (const Element *elmt)
{ return elmt->key; }

int Element_getLevel (const Element *elmt)
{ return level(elmt->key); }

int Element_getNpts (const Element *elmt)
{ return elmt->nr * elmt->ns; }

