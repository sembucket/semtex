#ifndef ALIAS_H
#define ALIAS_H

typedef struct alias {              /* ........ ALIAS definition ......... */
  int             type       ;      /* Node type (0, 1, or 2)              */
  struct {                          /*                                     */
    int a, b, edge;                 /* Key values of the edge and vertices */
  } key;                            /*                                     */
  double          pos[2]     ;      /* Node position                       */
  double**        Z          ;      /* Projection matrix                   */
} Alias;

Alias *Alias_alloc (void);
void   Alias_free  (Alias *alias);
void   Alias_build (Alias *alias, double *pos, int type, int np);

#endif
