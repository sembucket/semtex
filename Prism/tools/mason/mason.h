#ifndef MASON_H
#define MASON_H

/*
 * Type definitions and Prototypes for MASON
 *
 * RCS Information
 * --------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>

#define MAX_EDGES   4
#define MAX_VERTX   4

extern char *error_buf;            /* Character array for error messages  */
extern int   oplevel, opstart,     /* Optimization parameters             */
             rflag;                /* Randomize flag (anti-optimization)  */

typedef struct Point  {            /* ......... A Point in R2 ........... */
  double   x, y;                   /* (x,y)-coordinates                   */
} Point;

typedef int  Node;                 /* .............. A Node ............. */

typedef int                        /* ....... Boundary Conditions ....... */
  Simple;                          /* Simple Dirichlet/Neumann boundary   */
typedef struct {                   /*                                     */
  int   element;                   /*    to element                       */
  int   face   ;                   /*    to face                          */
} Connected;                       /* Element-Element boundary            */
typedef struct {                   /*                                     */
  int   patch  ;                   /*    patch number                     */
  int   segment;                   /*    segment number                   */
} Patched;                         /* Patched boundary                    */
typedef union {                    /* ................................... */
  Simple     s ;                   /*                                     */
  Connected  c ;		   /*                                     */
  Patched    p ;		   /*                                     */
} Binfo;

typedef struct Vertex {            /* ............. Vertex .............. */
  int             id       ;       /* ID number of the vertex             */
  int             mult     ;       /* Multiplicity                        */
  Point           vc       ;       /* Vertex coordinate                   */
  Node*           node     ;       /* Global node number                  */
  struct Vertex*  next     ;       /* Link to another vertex   [optional] */
} Vertex, *VertexP;

typedef struct Edge {              /* .............. Edge ............... */
  int             id       ;       /* ID number of the edge               */
  int             iel      ;       /* Element the edge is attached to     */
  int             np       ;       /* Number of points along the edge     */
  int             dir      ;       /* Direction to move along the edge    */
  char            type     ;       /* Boundary type for the edge          */
  Binfo           bc       ;       /* Boundary information                */
  Node*           nodes    ;       /* List of node numbers                */
  VertexP         right    ;       /* Pointer to the right vertex         */
  VertexP         left     ;       /* Pointer to the left vertex          */
  struct Edge*    next     ;       /* Link to another edge     [optional] */
} Edge, *EdgeP;
  
typedef struct Element {           /* ............ Element .............. */
  int             id       ;       /* ID number of the element            */
  int             nr, ns   ;       /* Number of points in (r,s)           */
  int             vertices ;       /* Number of vertices, 3 or 4          */
  int             edges    ;       /* Number of edges, 3 or 4             */
  VertexP*        vlist    ;       /* List of pointers-to-vertices        */
  Edge*           elist    ;       /* List of edges                       */
  struct Element* next     ;       /* Link to another element  [optional] */
} Element;         

/* ---------------------------------------------------------------------- *
 *                             P A T C H I N G                            *
 * ---------------------------------------------------------------------- */

typedef struct Segment {           /* ............ Segment .............. */  
  int             id       ;       /* ID number for the segment           */
  double          s0       ;       /* Segment origin in patch-coordinates */
  double          length   ;       /* Length of the segment               */
  int             branches ;       /* Number of connected branches        */
  int*            branch_ID;       /* ID numbers of the connections       */
  Edge*           edge     ;       /* Edge connected to this segment      */
  struct Segment* next     ;       /* Link to the next segment            */
} Segment, *SegmentP;

typedef struct Patch {             /* ............. Patch ............... */
  int             id       ;       /* Patch ID number                     */
  Point           origin   ;       /* Patch origin in (x,y)-coordinates   */
  Segment        *masters  ;       /* Pointers to master segments         */
  Segment        *slaves   ;       /* Pointers to slave segments          */
  struct Patch   *next     ;       /* Next patch in the list              */
} Patch, *PatchP;



typedef struct Domain {            /* ............. Domain .............. */
  int             norder   ;       /* Number of points in each direction  */
  int             nodes    ;       /* Number of global boundary nodes     */
  int             dof      ;       /* Number of true boundary nodes       */
  int             elements ;       /* Number of elements in the mesh      */
  int             verbose  ;       /* Verbose flag                        */
  int             tri      ;       /* Triangular element flag             */
  VertexP         master   ;       /* List of master vertices             */
  Element*        U        ;       /* List of elements                    */
  Patch*          P        ;       /* List of patches                     */
} Domain;


/* ------------------------------------------------------------------------ *
 *                          P R O T O T Y P E S                             *
 * ------------------------------------------------------------------------ */

int      main          (int argc, char *argv[]);
Domain  *parse_args    (int argc, char *argv[]);
void     error_msg     (char *msg);

void     loadParams    (Domain *omega, FILE *fp);
void     loadMesh      (Domain *omega, FILE *fp);
void     loadBCs       (Domain *omega, FILE *fp);

Patch*   makePatch     (int id, Patch   *link);
Segment* makeSegment   (int id, Segment *link, Edge *edge);
Patch*   findPatch     (int id, Domain  *omega);
Segment* findSegment   (int id, Segment *slist);
void     showPatch     (Domain *omega);
void     showSegment   (Segment *s);
void     linkPatches   (Domain *omega);

Vertex*  findVertex    (Vertex *list, Point p);
Vertex*  makeVertex    (Vertex *link, Point p);
void     linkVertex    (Domain *omega);
void     showVertex    (Domain *omega);

void     linkEdges     (Domain *omega);
Edge*    findEdge      (Domain *omega, int element, int face);
void     showEdges     (Domain *omega);
void     linkNodes     (Domain *omega);
void     linkNodesE    (Domain *omega, Element *U);

void     optimize      (Domain *omega);
void     output        (Domain *omega, FILE *fp);
int      bandwidth     (Domain *omega);
void     RCMop         (Domain *omega);
void     greedy        (Domain *omega);

#endif
