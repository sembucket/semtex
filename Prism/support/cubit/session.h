#ifndef SESSION_H
#define SESSION_H

#include <stdio.h>
#include "mesh.h"
#include "param.h"

typedef void (*KeywordProcessor)(FILE *fp);

/* Parse the input file.  */
void session_parse (FILE *fp, param_t *param, Mesh *mesh);

/* Map a keyword to a function for processing it.  Each time the keyword is  *
 * encountered in the input file, the stream is passed to the given func-    *
 * tion for handling the keyword's data.  The keyword will be skipped over   *
 * if it doesn't map to any known function.                                  */
int session_defKeyword (const char *keyword, KeywordProcessor ptr);

#endif
