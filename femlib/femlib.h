/*****************************************************************************
 *         FUNCTION PROTOTYPES FOR ROUTINES IN LIBRARY FEMBASE.A             *
 *****************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/

/* $Id$ */


#ifndef FEM_H		/* BEGIN fembase.h DECLARATIONS */
#define FEMBASE_H

/* ------------------------------------------------------------------------- *
 * Routines from initial.c                                                   *
 * ------------------------------------------------------------------------- */

void   initialize  (void);
double interpret   (char *);

void   set_option  (char *, int);
int    get_option  (char *);
void   show_option (FILE *);

void   set_iparam  (char *, int);
int    get_iparam  (char *);
void   show_iparam (FILE *);

void   set_dparam  (char *, double);
double get_dparam  (char *);
void   show_dparam (FILE *);

#endif
