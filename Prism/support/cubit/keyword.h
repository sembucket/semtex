#ifndef KEYWORD_H
#define KEYWORD_H

/* Keyword
 *
 * $Revision$
 *
 * Author: R. D. Henderson
 *
 * This is a generic data structure for storing keyword information from an
 * input file.  A keyword declaration and its assocated data appear in the
 * input file as follows:
 *
 *      ***** KEYWORD *****
 *      <n> lines of information
 *      line number 1
 *      ...
 *      line number N
 *
 * When a keyword is exported, the above block of information is written to
 * the given stream.  When a keyword is imported, the stream should be posi-
 * tioned so that the next item is the information count.  The keyword 
 * identification line is not read.
 *
 * Example:
 *
 *      ***** COLORS *****
 *      3 lines of information
 *      red   = 100
 *      blue  = 200
 *      green = 300
 *
 * ------------------------------------------------------------------------- */

#include <stdio.h>

#define KEYWORD_MAXINFO 32
#define KEYWORD_COUNT(k)   ((k)->count)
#define KEYWORD_INFO(k,i)  ((k)->info[(i)])

typedef struct {
  int   count;
  char *name;
  char *info[KEYWORD_MAXINFO];
} keyword_t;

keyword_t *keyword_alloc(const char *name);
void       keyword_free (keyword_t *keyword);

int keyword_add (keyword_t *keyword, const char *info);

int keyword_export (const keyword_t *keyword, FILE *fp);
int keyword_import (      keyword_t *keyword, FILE *fp);

#endif
