#ifndef PARAM_H
#define PARAM_H

#ifdef __cpluscplus
extern "C" {
#endif

/* Parameters
 *
 * $Revision$
 *
 * Author:     R. D. Henderson
 *
 * This data structure simply stores parameters from the input file.  For an
 * ISL input file, the list of parameters has the form
 *
 *      param {
 *         token = expr
 *         ...
 *      }
 * 
 * while for a NEKTON input file it looks like
 *
 *      ***** PARAMETERS *****
 *      ----- comment
 *      ----- comment
 *      10 PARAMETERS FOLLOW
 *      expr  token
 *      ...
 *
 * The data structure param_t is allocated to store a maximum number of 
 * parameter values.  Once allocated, you can do the following things to it:
 *
 *      _define = define tokens
 *      _import = import tokens into the list from a FILE
 *      _export = export tokens to a FILE
 *
 * Note that whenever parameters are imported, they are automatically fed
 * into cubit's global symbol table.
 * ------------------------------------------------------------------------- */

#define PARAM_COUNT(p) (p->count)

typedef struct {
  int           max;
  int           count;
  struct token *token;
} param_t;

param_t *param_alloc (int max);
void     param_free  (param_t *p);

int param_define (param_t *p, const char *token, const char *expr);

int param_export (const param_t *p, FILE *fp);
int param_export1(const param_t *p, FILE *fp);
int param_import (      param_t *p, FILE *fp);

char *param_name (const param_t *p, int i);
char *param_expr (const param_t *p, int i);

#ifdef __cpluscplus
}
#endif
#endif

