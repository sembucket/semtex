#ifndef DOMAIN_H
#define DOMAIN_H

/* A computational domain
 *
 * $Id$
 * ------------------------------------------------------------------------- */

typedef enum {     /* formats for storing a mesh */
  M_NEKTON,
  M_SEMTEX,
  M_NONE
} domain_fmt_t;

typedef struct {
  char      *name;
  param_t   *param;

  keyword_t *force;
  keyword_t *solution;
  keyword_t *ic;
  keyword_t *history;
  keyword_t *fields;
  keyword_t *user;

  Mesh      *mesh;
  BC        *BC;
  Matrix    *A;

  int        nfields;
  Field     *solVector[16];

} Domain;

int Domain_init     (const char *name);
int Domain_reset    ();
int Domain_save     ();
int Domain_load     ();
int Domain_check    ();
int Domain_require  ();
int Domain_bbox     ();

Field* Domain_getField (char type);
Field* Domain_chkField (char type);
Field* Domain_addField (char type);

#endif

