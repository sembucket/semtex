/* Implementation of param_t
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "manager.h"
#include "param.h"

struct token {
  char *name;
  char *expr;
};

param_t *param_alloc (int max)
{
  param_t *param = (param_t*) malloc(sizeof(param_t));

  param->count = 0;
  param->max   = max;
  param->token = (struct token*) calloc(max,sizeof(struct token));

  return param;
}

void param_free (param_t *param)
{
  param_clear(param);

  free (param->token);
  free (param);
}

/* Clear the list */

int param_clear (param_t *param)
{
  int i;
  for (i = 0; i < param->count; i++) {
    free (param->token[i].name);
    free (param->token[i].expr);
  }
  return 0;
}

/* Method to give a parameter's symbolic definition.  Also sets the value of *
 * the parameter in the global symbol table based on the deduced type.       */

int param_define (param_t *param, const char *name, const char *expr)
{
  const int count = param->count;
  int i, n;

  static char *dspecial[] = { "LAMBDA", "LZ", 0 };
  static char *ispecial[] = { "EQTYPE", "TORDER", 0 };

  for (n = 0; n < count; n++) {
    if (strcmp(param->token[n].name,name)==0) {
      free (param->token[n].expr);
      param->token[n].expr = strdup(expr);
      break;
    }
  }

  if (n == count) {
    if (count == param->max)
      cubit_err ("param: too many parameters!");
    param->token[n].name = strdup(name);
    param->token[n].expr = strdup(expr);
    param->count++;
  }

  /* Here we store the param's symbolic definition, but we also need to assign
   * its numerical value in one of the global symbol tables.  This should be
   * either the integer (iparam) or floating-point (dparam) table.           */

  for (i = 0; ispecial[i] != NULL; i++) {
    if (strcmp(name,ispecial[i])==0) {
      iparam_set(name,(int)parse(expr)); 
      return n;
    }
  }

  for (i = 0; dspecial[i] != NULL; i++) {
    if (strcmp(name,dspecial[i])==0) {
      dparam_set(name,parse(expr)); 
      return n;
    }
  }

  if ('I' <= *name && *name <= 'N')
    iparam_set(name, (int) parse(expr));
  else
    dparam_set(name, parse(expr));
  
  return n;
}

/* Return a pointer to a parameter's name */

char *param_name (const param_t *p, int i)
{
  if (0 <= i <= p->count)
    return p->token[i].name;
  else
    return NULL;
}

/* Return a pointer to a parameter's expression */

char *param_expr (const param_t *p, int i)
{
  if (0 <= i <= p->count)
    return p->token[i].expr;
  else
    return NULL;
}

/* Import a list of parameters from an rea file */

int param_import (param_t *param, FILE *fp)
{
  int count = param->count;

  int i, n;
  char buf[BUFSIZ];
  char *fundSection(char*,char*,FILE*);

  rewind(fp); findSection("PARAMETERS", buf, fp);

  if (option("format")==0) {
    for (i = 0; i < 2; i++)      /* Read the two comment lines */
      fgets(buf, BUFSIZ, fp);
  }

  fgets(buf, BUFSIZ, fp);
  if (sscanf(buf, "%d", &n) != 1)
    cubit_err ("param_import: can't read the # of parameters");
  if (count + n > param->max)
    cubit_err ("param_import: too many parameters!");

  count += n;
  for (i = count-n; i < count; i++) {
    char name[32];
    char expr[32];
    fgets (buf, BUFSIZ, fp);
    sscanf(buf, "%25s%25s", expr, name);
    param_define(param, name, expr);
  }

  return param->count;
}

/* Export a list of parameters to an rea file */

int param_export (const param_t *param, FILE *fp)
{
  const int count = param->count;

  if (option("format")==0) {
    fprintf (fp, "***** PARAMETERS *****\n");
    fprintf (fp, "# comment line -- do not remove\n");
    fprintf (fp, "# comment line -- do not remove\n");
    fprintf (fp, "%d item%c\n", count, count != 1 ? 's' : ' ');
  } else {
    fprintf (fp, "***** PARAMETERS *****\n");
    fprintf (fp, "%d item%c\n", count, count != 1 ? 's' : ' ');
  }

  return param_export1(param,fp);
}

/* Export form:  token expr */

int param_export1 (const param_t *param, FILE *fp)
{
  const int count = param->count;
  int i;

  for (i = 0; i < count; i++)
    fprintf (fp, "%s %s\n", param->token[i].expr, param->token[i].name);
  
  return 0;
}

