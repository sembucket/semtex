#ifndef HISTORY_H
#define HISTORY_H

#include "prism/constants.h"
#include "speclib/speclib.h"
#include "speclib/probe.h"

typedef struct hpnt {             /* ......... HISTORY POINT ......... */
  int         id         ;        /* ID number of the history point    */
  int         frame      ;        /* Frame number to sample on         */
  Probe*      locator    ;        /* Location in the mesh              */
  char        fields              /* Fields to sample at this point    */
            [_MAX_FIELDS];        /*                                   */
  ACTION      mode       ;        /* Physical or Fourier space         */   
  struct hpnt *next      ;        /* Pointer to the next point         */
} HisPoint;

#endif
