/*
 * Build Connectivity Information
 *
 * RCS Information
 * ---------------------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "vdb/vdb.h"
#include "veclib/veclib.h"

#include "cubit.h"
#include "alias.h"
#include "vertex.h"
#include "edge.h"
#include "element.h"

#define VERTEX   0
#define EDGE     1
#define BC       2

#define MAXENTRY    8*_MAX_NEL

#define edge_get(elmt,id)   (elmt->edges  + (id))
#define vertex_get(elmt,id) (elmt->vertex + (id))

/* Private functions */

static void connect (Mesh *mesh);
static void rebuild (Mesh *mesh);
static void resolve (Mesh *mesh, VDB nodes, int nnodes);

/* ------------------------------------------------------------------------- *
 * Mesh_connect() -- Establish Connectivity                                  * 
 *                                                                           *
 * This function establishes the connectivity of the mesh, including both    *
 * conforming and nonconforming edges.  It requires significant communica-   *
 * tion between processors to coordinate the establishment of "aliases".     *
 * ------------------------------------------------------------------------- */

void Mesh_connect (Mesh *mesh)
{
  if (mesh->edges == NULL || mesh->vertx == NULL) 
    rebuild (mesh);

  connect (mesh);
  return;
}

static void connect (Mesh *mesh)
{
  VDB  nodes;
  int  nnodes, nedges;
  int  ntype[MAXENTRY];
  int  nmult[MAXENTRY];
  int  emult[MAXENTRY];

  Element *elmt;
  Edge    *edge;
  Vertex  *vert;

  int i;
  const int bufsiz = 1 << 15;

  double tol = dparam("TOLVDB"), tsec;
  void*  buf = malloc(bufsiz);
  assert(buf != NULL);

  /* Free previous aliases */

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    vert = elmt->vert_list;
    edge = elmt->edge_list;
    for (i = 0; i < 4; i++) {
      vert[i].alias = NULL;
      if (edge[i].alias) {
	free (edge[i].alias);
	edge[i].alias = NULL;
      }
    }
  }
    
  /* Create NODES (edges + vertices) data base */

  tsec  = dclock();
  nodes = vdbcreate (2, tol, 0, 0, MAXENTRY, buf, bufsiz);
  memset (ntype, '\0', MAXENTRY*sizeof(int));

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    for (i = 0, vert = elmt->vert_list; i < 4; i++)
      vdbdjoin (nodes, vert[i].pos);
    for (i = 0, edge = elmt->edge_list; i < 4; i++) {
      const int key = vdbdjoin (nodes, edge[i].pos);
      ntype[key] = edge[i].bc ? BC : EDGE;
    }
  }

  nnodes = vdbsync (nodes);
  vdbcombine (nodes, ntype, vdbimax, sizeof(int), 1, sizeof(int));
  vdbnmember (nodes, nmult);

  /* Compute the edge multiplicity based on ACTIVE edges */

  nedges = vdbnmember(mesh->edges, NULL);
  memset (emult, '\0', nedges*sizeof(int));

  for (elmt = mesh->head; elmt; elmt = elmt->next)
    for (i = 0, edge = elmt->edge_list; i < 4; i++)
      emult[edge[i].key]++;
  
  vdbcombine (mesh->edges, emult, vdbisum, sizeof(int), 1, sizeof(int));
  
  /* Scan for unlinked edges */
  
  for (elmt = mesh->head; elmt; elmt = elmt->next)
    for (i = 0, edge = elmt->edge_list; i < 4; i++) {
      
      if (emult[edge[i].key] == 1 && edge[i].bc == NULL) {
	if (nmult[vdbdpeek(nodes,edge[i].pos)] > 1) {
	  edge[i].alias = Alias_alloc ();
	  Alias_build (edge[i].alias, edge[i].pos, 0, 0);
	} else {
	  Vertex *a = edge[i].a;
	  Vertex *b = edge[i].b;
	  if (ntype[vdbdpeek(nodes,a->pos)] == EDGE) {
	    edge[i].alias = a->alias = Alias_alloc ();
	    Alias_build (edge[i].alias, a->pos, 1, edge[i].np);
	  } else {
	    edge[i].alias = b->alias = Alias_alloc ();
	    Alias_build (edge[i].alias, b->pos, 2, edge[i].np);
	  }
	}
      }
    }

  /* Rebuild the database of edge and vertex positions, clean up */

  rebuild    (mesh);
  resolve    (mesh, nodes, nnodes);
  free       (buf);
  vdbdestroy (nodes);

  tsec = dclock() - tsec;
}

/* ------------------------------------------------------------------------- *
 * rebuild() -- Rebuild edge and vertex database                             *
 * ------------------------------------------------------------------------- */

static void rebuild (Mesh *mesh)
{
  Element *elmt;
  Edge    *edge;
  Vertex  *vert;
  int i;
  
  const int bufsiz = 1 << 15;  /* 32 Kb buffer */
  static void *buf = NULL;

  double tol = dparam("TOLVDB");

  if (tol == 0.) dparam_set("TOLVDB", tol = FLT_EPSILON);

  if (buf == NULL) {       /* Internal buffer */
    buf = malloc (bufsiz);
    assert (buf != NULL);
  }

  if (mesh->edges) vdbdestroy (mesh->edges);
  if (mesh->vertx) vdbdestroy (mesh->vertx);
  
  mesh->vertx = vdbcreate (2, tol, 0, 0, 4*_MAX_NEL, buf, bufsiz);
  mesh->edges = vdbcreate (2, tol, 0, 0, 4*_MAX_NEL, buf, bufsiz);

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    for (i = 0, vert = elmt->vert_list; i < 4; i++)
      vert[i].key = vdbdjoin (mesh->vertx, vert[i].pos);
    for (i = 0, edge = elmt->edge_list; i < 4; i++) {
      edge[i].key = vdbdjoin (mesh->edges, edge[i].pos);
      if (edge[i].alias || edge[i].bc)
	edge[i].dir = 1;
      else 
	edge[i].dir = edge[i].b->key > edge[i].a->key ? 1 : -1;
    }
  }

  vdbsync (mesh->vertx);
  vdbsync (mesh->edges);

  return;
}

/* ------------------------------------------------------------------------- *
 * resolve() -- Resolve aliases                                              *
 *                                                                           *
 * For a distributed mesh, the "position" of an alias is the midpoint of the *
 * edge that defines it.  It is always known locally, but its vertices may   *
 * not be.  To resolve all aliases, the positions of their vertices are      *
 * written to an array (four doubles per node in 2D), the array of positions *
 * is combined, and vertex positions are then read by each aliased edge.     *
 * ------------------------------------------------------------------------- */

static void resolve (Mesh *mesh, VDB nodes, int nnodes)
{
  Element *elmt;
  Edge    *edge;
  Vertex  *vert;
  double  *pos;
  int i;

  const size_t s = 2 * sizeof(double);
  
  pos = (double*) calloc (4*nnodes, sizeof(double));

  /* Write node position data for the parent edges */
  
  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    for (i = 0, edge = elmt->edge_list; i < 4; i++) 
      if (edge[i].alias && edge[i].alias->type == 0) {
	Vertex *a = edge[i].a;
	Vertex *b = edge[i].b;
	double *p = pos + 4 * vdbdpeek(nodes, edge[i].pos);
	
	memcpy (p,     a->pos, s);
	memcpy (p + 2, b->pos, s);

	edge[i].alias->key.a    = a->key;
	edge[i].alias->key.edge = edge[i].key;
	edge[i].alias->key.b    = b->key;
      }
  }

  /* ----- Combine ----- */
  
  vdbcombine (nodes, pos, vdbdsum, sizeof(double), 4, 4*sizeof(double));

  /* Read node position data into children */

  for (elmt = mesh->head; elmt; elmt = elmt->next) {
    for (i = 0, edge = elmt->edge_list; i < 4; i++) 
      if (edge[i].alias && edge[i].alias->type > 0) {
	double *p = pos + 4 * vdbdpeek(nodes, edge[i].alias->pos);
	
	edge[i].alias->key.a    = vdbdjoin (mesh->vertx, p);
	edge[i].alias->key.edge = vdbdjoin (mesh->edges, edge[i].alias->pos);
	edge[i].alias->key.b    = vdbdjoin (mesh->vertx, p + 2);
      }
  }

  vdbsync (mesh->vertx);
  vdbsync (mesh->edges);

  free (pos);
  return;
}
