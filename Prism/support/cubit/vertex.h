#ifndef VERTEX_H
#define VERTEX_H

/* VERTEX
 *
 * $Revision$
 *
 * Author:    R. D. Henderson
 *
 * ------------------------------------------------------------------------- */

#define VERTEX_ID(vert)    (vert->id)
#define VERTEX_KEY(vert)   (vert->key)
#define VERTEX_POS(vert,i) (vert->pos[i])
#define VERTEX_BC(vert)    (vert->bc)
#define VERTEX_ALIAS(vert) (vert->alias)

typedef struct vertex {             /* ........ VERTEX definition ........ */
  int             id         ;      /* ID number                           */
  int             key        ;      /* Key to the global vertex database   */
  double          pos[2]     ;      /* Position                            */
  struct bc      *bc         ;      /* Pointer to a boundary condition     */
  struct alias   *alias      ;      /* Pointer to an alias                 */
  struct vertex  *next       ;      /* Link to the next one                */
} Vertex;

/* Prototypes */

Vertex *Vertex_valloc();
void    Vertex_vfree (Vertex *vertx);

#endif
