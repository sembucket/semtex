#ifndef BC_H
#define BC_H

#include "edge.h"
#include "element.h"

#define BC_MAXINFO 8  /* Maximum number of info items */

#define BC_TYPE(bc)   (bc->type)
#define BC_ITEMS(bc)  (bc->nitems)
#define BC_INFO(bc,i) (bc->info[i])

typedef enum {
  BC_ESSENTIAL,                    /* Essential (Dirchlet) BC data */
  BC_NATURAL,                      /* Natural (Neumann) BC data */
  BC_MIXED,                        /* Mixed (D + N) BC data */
} BC_class;

    
typedef struct bc {                /* ....... BOUNDARY Definition ....... */
  char            type       ;     /* Type                                */
  int             nitems     ;     /* How many pieces of information      */
  union {                          /* Boundary condition info...          */
    int           number     ;     /*    - number                         */
    double        value      ;     /*    - value                          */
    char*         expr       ;     /*    - expression                     */
    void*         other      ;     /*    - anything else                  */
  } info [BC_MAXINFO];             /*                                     */
  struct bc       *next      ;     /* Pointer to the next one             */
} BC;


/* Prototypes */

BC*  BC_alloc  (char type, int nitems);
void BC_free   (BC *bc);
void BC_attach (BC *bc, Element *elmt, Edge *edge);

#endif

