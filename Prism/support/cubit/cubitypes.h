/*
 * Miscellaneous type declarations
 *
 * RCS Information
 * ----------------------------
 * $Author$
 * $Date$
 * $Revision$
 * ----------------------------------------------------------------------- */

#ifndef CUBITYPES_H           
#define CUBITYPES_H

#include <stdio.h>
#include "element.h"
#include "matrix.h"

typedef enum {               /* --------- List of flags --------- */
  OFF            =  0,       /* General-purpose flags for on/off  */
  ON             =  1,       /*                                   */
  ERROR          =  2,       /* Non-recoverable error             */
  WARNING        =  3,       /* Not too serious                   */
  UNKNOWN        = -1,       /* Unknown boundary condition        */
  DIRICHLET      =  0,       /* Essential   "        "            */
  NEUMANN        =  1        /* Natural     "        "            */
} FLAG;

#endif                       /* END OF SEMTYPES.H DECLARATIONS */

