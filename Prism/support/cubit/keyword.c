/* Implementation of keyword
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "keyword.h"

keyword_t *keyword_alloc(const char *name)
{
  keyword_t *k = (keyword_t*) malloc(sizeof(keyword_t));
  k->count = 0;
  k->name  = strdup(name);
  return k;
}

void keyword_free (keyword_t *k)
{
  int n;
  for (n = 0; n < k->count; n++)
    free (k->info[n]);
  free (k->name);
  free (k);
}

int keyword_add (keyword_t *k, const char *info)
{
  int n = k->count;
  if (n < KEYWORD_MAXINFO) {
    k->count++;
    k->info[n] = strdup(info);
  } else {
    fprintf (stderr, "keyword: too much information!\n");
    n = -1;
  }
  return n;
}

int keyword_export (const keyword_t *k, FILE *fp)
{
  const int count = k->count;
  int n;
  fprintf (fp, "***** %s *****\n", k->name);
  fprintf (fp, "%d item%c\n", count, count==1 ? ' ' : 's');
  for (n = 0; n < count; n++)
    fprintf (fp, "%s\n", k->info[n]);
  return n;
}

int keyword_import (keyword_t *k, FILE *fp)
{
  char buf[BUFSIZ], *p;
  int n;

  fgets(buf, BUFSIZ, fp);
  if (sscanf(buf, "%d", &k->count) != 1) {
    sprintf (buf, "keyword: can't read number of lines for %s\n", k->name);
    cubit_err(buf);
  } else {
    for (n = 0; n < k->count; n++) {
      fgets(buf, BUFSIZ, fp);
      buf[strlen(buf)-1]='\0';
      k->info[n] = strdup(buf);
    }
  }

  return n;
}
