/*
 * alias_t implementation
 *
 * $Id$
 *
 * Author: R. D. Henderson
 *
 * Copyright (c) 1998-1999 R. D. Henderson and Caltech.
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include "alias.h"

alias_t *alias_alloc() {
  alias_t *alias = (alias_t*) calloc(1,sizeof(alias_t));
  return alias;
}

void alias_free (alias_t *alias) {
  free(alias);
}

