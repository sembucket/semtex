#ifndef DOMAIN_H
#define DOMAIN_H

/* Domain
 *
 * $Id$
 *
 * This is the data structure that represents the complete computational 
 * domain, including fields for the primary variables, workspace, linear
 * systems, flags, Greens functions and output files.
 *
 * The only way to allocate a domain is to provide the name of a session
 * file.  The function Domain_alloc() reads all information from the session
 * file and constructs a new mesh, builds matrices, sets boundary conditions
 * and computes initial conditions.
 *
 * Copyright (c) 1994-1998 R. D. Henderson
 * ------------------------------------------------------------------------- */

#include <stdio.h>

#include "prism/constants.h"
#include "prism/greens.h"
#include "prism/history.h"

#include "speclib/speclib.h"
#include "speclib/field.h"

typedef struct domain {           /* ............ DOMAIN ............. */
  char     *name                ; /* Session name                      */
  
  Field    *U,  *V,  *W,  *P    ; /* Velocity and Pressure fields      */

  Field    *Us[_MAX_TORDER]     ; /* --------------------------------- */
  Field    *Vs[_MAX_TORDER]     ; /*                                   */
  Field    *Ws[_MAX_TORDER]     ; /*                                   */
  Field    *Uf[_MAX_TORDER]     ; /*        Multi-step storage         */
  Field    *Vf[_MAX_TORDER]     ; /*                                   */
  Field    *Wf[_MAX_TORDER]     ; /*                                   */
                                  /* --------------------------------- */
  Bedge    *Ubc, *Vbc, *Wbc     ; /*        Boundary conditions        */
  Bedge    *Pbc                 ; /* --------------------------------- */
  BSystem  *Velocity            ; /* Boundary matrix system (velocity) */
  BSystem  *Pressure            ; /* Boundary matrix system (pressure) */

  GreensF  *G                   ; /* Green's Function       (optional) */

  Field    *Ux, *Uy, *Uz        ; /* Velocity field (saved - optional) */
  Field    *Qx, *Qy, *Qz        ; /* Vorticity field        (optional) */

  /* The rest of this is output data, including output streams and the *
   * parameters that control how often output is performed.            */

  FILE      *fld_file           ; /* Solution file                     */
  FILE      *his_file           ; /* History point file                */
  HisPoint  *his_list           ; /* History point list                */

  struct measure *mea           ; /* Measurements data structure       */
  struct stat    *stats         ; /* Statistics data structure         */

  struct {                        /* Output frequencies                */
    int     cfl;
    int     history;
    int     measure;
    int     field;
    int     stats;
  } step;

} Domain;

/* Prototypes */

Domain *Domain_alloc   (const char *session);
void    Domain_free    (Domain *domain);

void    Domain_save    (Domain *domain);
void    Domain_checkpt (Domain *domain);

#endif
