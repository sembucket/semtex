#ifndef MESH_H
#define MESH_H

/* MESH
 *
 * $Revision$
 *
 * Author:     R. D. Henderson
 *
 * This is the data structure used to represent and manipulate a spectral
 * element mesh.  A mesh consists of a collection of elements at different
 * refinement levels.  The levels are organized as a quadtree.  The tree
 * structure is only used to define the connections between various levels.
 * There is a subset of elements which are considered "active".  These de-
 * fine the leaf nodes of the tree that make up the current discretization 
 * of space.  Other elements are "inactive" and generally represent parents 
 * of active leaf nodes in the tree.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "element.h"

/* Some limits... */

#define _MAX_NEL       16384   /* Elements in a single mesh                */
#define _MAX_NORDER       16   /* Maximum Polynomial order in space        */
#define _MAX_NB          256   /* Points along the edge of an element      */
#define _MAX_FIELDS       24   /* Maximum fields in a solution file        */

/* Macros */

#define MESH_NELMT(mesh) (mesh->active)
#define MESH_NVERT(mesh) (0)
#define MESH_NEDGE(mesh) (0)
#define MESH_NR(mesh)    (mesh->nr)
#define MESH_NS(mesh)    (mesh->ns)
#define MESH_NZ(mesh)    (mesh->nz)
#define MESH_NELMT(mesh) (mesh->active)
#define MESH_NPTS(mesh)  (MESH_NR(mesh)*MESH_NS(mesh)*MESH_NELMT(mesh))

#define MESH_HEAD(mesh)  (mesh->head)
#define MESH_ELEMENT(mesh,k) (mesh->list[k])

typedef struct mesh {               /* ......... MESH Definition ......... */
  int            elements       ;   /* Number of elements (total)          */
  int            active         ;   /* Number of elements (active)         */
  int            nr, ns, nz     ;   /* Number of points in each direction  */

  struct vdb*    vertx          ;   /*    VOXEL Database of Positions      */
  struct vdb*    edges          ;
  struct vdb*    elmts          ;

  struct bc*      bc            ;   /* Boundary conditions                 */
  struct element* head          ;   /* Head of the active element list     */
  struct element* list[_MAX_NEL];   /* List of all elements in the mesh    */
} Mesh;

/* Prototypes */

Mesh* Mesh_alloc      ();
void  Mesh_free       (Mesh *mesh);
void  Mesh_info       (Mesh *mesh);

void  Mesh_add        (Mesh *mesh, Element *elmt);
void  Mesh_connect    (Mesh *mesh);
int   Mesh_import     (Mesh *mesh, FILE *fp);
void  Mesh_install    (Mesh *mesh, int family, int key);
int   Mesh_lookup     (Mesh *mesh, int family, int key);
void  Mesh_prune      (Mesh *mesh, int family, int depth);
void  Mesh_refine     (Mesh *mesh);
void  Mesh_resetFlags (Mesh *mesh);
void  Mesh_restrict   (Mesh *mesh);

void  Mesh_getSurface (Mesh *mesh, char type);

#endif
  
