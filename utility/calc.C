/*****************************************************************************
 * CALC: a basic calculator using the function parser.
 *
 * Usage: calc [-]
 *****************************************************************************/

// $Id$ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Femlib.h>

static char buf[FILENAME_MAX];


int main (int argc, char *argv[])
/* ========================================================================= *
 * Driver.
 * ========================================================================= */
{
  char  *c;

  initialize();

  while (gets (buf)) {
    if ( (c = strstr (buf, "=")) )
      interpret(buf);
    else
       printf("%-.17g\n", interpret(buf));
  }

  return EXIT_SUCCESS;
}
