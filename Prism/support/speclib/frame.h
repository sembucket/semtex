#ifndef FRAME_H
#define FRAME_H

/* FRAME
 * 
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 *
 * A frame is the data structure that holds the global collection of 
 * elements together.  It contains all of the global information, such as
 * connectivity information, etc.
 * ------------------------------------------------------------------------- */

#include "speclib/field.h"
#include "speclib/matrix.h"

typedef struct frame {
  int    elements;
  int    families;

  int    npts, ndof;
  int    bpts, bdof;
  int    ipts;

  double *mass;
  double *massinv;
  void   *other;
} Frame;

/* Prototypes */

int  Frame_init    (Field *U, BSystem *B, int frame, const char *name);
void Frame_load    (Field *U, BSystem *B, int frame);
void Frame_save    (Field *U, BSystem *B, int frame);
void Frame_set     (int frame, int nfields, ... /* list of fields */);
void Frame_set_one (int frame, Field *U);

#endif
