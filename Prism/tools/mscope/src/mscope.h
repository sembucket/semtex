#ifndef MSCOPE_H
#define MSCOPE_H

/* Mscope Definitions and Prototype declarations 
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include "cubit/cubit.h"
#include "domain.h"

/* Definitions */

typedef struct Point {
  double x, y;
} Point;

/* Prototypes */

void parse_args (int argc, char *argv[]);

int mscope_refine (Element *parent, Element *child[]);
int mscope_bc     ();

int  DoParse (char *command_line);

void GetLimits (void);
void SetLimits (Domain *gp, ...);
void ZoomIn    (void);

#endif



