/*****************************************************************************
 * BOUNDARY.C:  utilities for dealing with boundary conditions.              *
 *                                                                           *
 * This module maintains a list of tagged velocity boundary conditions.      *
 *                                                                           *
 * Boundary conditions are given numeric tags by which the element-wise      *
 * description can access the B-C information.  Tag-number zero is reserved  *
 * to label inter-element boundaries.  It is implicit in the descriptions    *
 * below that on each tagged boundary the same kind of condition is used     *
 * for each variable.                                                        *
 *                                                                           *
 * The tag line is also used to declare the kind of boundary condition.      *
 * Available kinds are (case-insensitive):                                   *
 *   ESSENTIAL                                                               *
 *   VALUE                                                                   *
 *   NATURAL                                                                 *
 *   FLUX                                                                    *
 *   WALL                                                                    *
 *   OUTFLOW                                                                 *
 *                                                                           *
 * For WALL and OUTFLOW BCs, no further information is needed, since they    *
 * implemented internally: a wall is just a special case of a VALUE BC on    *
 * which all velocity components are zero.  OUTFLOW boundaries set du/dn=0   *
 * for all velocity components (this may be ammended later), and so are a    *
 * special case of a FLUX boundary condition.                                *
 *                                                                           *
 * On VALUE boundaries, DIM velocity components need to be set, all on the   *
 * tag line.  Either numerical or interpretable functions values can be      *
 * used.  The order in which they are given must be the same as used for     *
 * space dimensions.  ESSENTIAL is the same as VALUE.                        *
 *                                                                           *
 * On FLUX boundaries, DIM normal gradients of velocities need to be set, in *
 * the same way as for VALUE boundaries.  NATURAL is the same as FLUX.       *
 *                                                                           *
 * Boundary conditions for secondary variables (e.g. pressure) are set       *
 * automatically.                                                            *
 *                                                                           *
 * The list of tagged BCs should be terminated with a blank line.            *
 *                                                                           *
 * Example for 2D (DIM=2):                                                   *
 * NB: value BCs are evaluated once, and may use variables "x" & "y".        *
 *                                                                           *
 * BOUNDARY CONDITIONS                                                       *
 *                                                                           *
 * 1   VALUE      1.0  sin(PI^2*x*y)                                         *
 * 2   WALL                                                                  *
 * 3   OUTFLOW                                                               *
 *                                                                           *
 *****************************************************************************/

static char
RCSid[] = "$Id$";


#include "Fem.h"


#define  readLine(fp)  fgets(buf, STR_MAX, (fp))
#define  echo          if (verb) fputs(buf, stdout)
#define  upperCase(s)  { char *z=(s); while (*z=toupper(*z)) z++; }


typedef struct b_c         {	/* ------------- B-C data type ------------- */
  int              id      ;	/* Identification tag                        */
  char            *name    ;	/* Array of character variable names         */
  BcRecord        *bCond   ;	/* Array of corresponding BCs                */
  struct b_c      *next    ;	/* Link to next one.                         */
} B_C;				/* ----------------------------------------- */

static B_C *BCList = NULL;
static char buf[STR_MAX];





void  readBCs (FILE *fp)
/* ========================================================================= *
 * Create & maintain a list of boundary condition structures.                *
 * ========================================================================= */
{
  char    routine[] = "readBCs";
  char    sep[]     = " \t,:;";
  char   *token, *suffix = NULL;
  B_C    *p;
  int     i, nvar = iparam ("NVAR"), verb = iparam ("VERBOSE");
  double  val;


  readLine (fp);

  while (readLine(fp)) {
    if (*buf == '\n') return;
    echo;

    p = (B_C *) calloc (1, sizeof (B_C));

    p -> name = (char *) calloc (STR_MAX, sizeof (char));
    switch (nvar) {
    case 1: strcpy (p -> name, "u"  );  break;
    case 2: strcpy (p -> name, "uvp" ); break;
    case 3: strcpy (p -> name, "uvwp"); break;
    default: 
      sprintf (buf, "number of variables must be in [1..3]: %1d", nvar);
      message (routine, buf, ERROR);
      break;
    }

    if (!(token = strtok (buf, sep)))
      message (routine, "empty line instead of BC", ERROR);

    if (sscanf (token, "%d", &p->id) != 1)
      message (routine, strcat(buf, "?: can't parse integer tag"), ERROR);
    if (p -> id < 1)
      message (routine, "BC tags must be >= 1", ERROR);
    p -> id--;

    if (!(token = strtok (NULL, sep)))
      message (routine, "can't parse kind of BC", ERROR);
    upperCase (token);

    p -> bCond = (BcRecord *) calloc (nvar, sizeof (BcRecord));

    if (strstr (token, "VALUE") || strstr (token, "ESSENTIAL")) {

      for (i = 0; i < nvar; i++) {
	if (!(token = strtok (NULL, sep)))
	  message (routine, "can't parse enough BC values", ERROR);
	
	val = strtod (token, &suffix);
	if (*suffix && *suffix != '\n') {
	  p -> bCond[i].kind         = ESSENTIAL_FN;
	  p -> bCond[i].value.interp =
	    (char *) calloc (strlen (token) + 1, sizeof (char));
	  strcpy (p -> bCond[i].value.interp, token);
	} else {
	  p -> bCond[i].kind        = ESSENTIAL;
	  p -> bCond[i].value.given = val;
	}
      }

    } else if (strstr (token, "FLUX") || strstr (token, "NATURAL")) {

      for (i = 0; i < nvar; i++) {
	if (!(token = strtok (NULL, sep)))
	  message (routine, "can't parse enough BC values", ERROR);
	
	val = strtod (token, &suffix);
	if (*suffix) {
	  p -> bCond[i].kind         = NATURAL_FN;
	  p -> bCond[i].value.interp =
	    (char *) calloc (strlen (token) + 1, sizeof (char));
	  strcpy (p -> bCond[i].value.interp, token);
	} else {
	  p -> bCond[i].kind        = NATURAL;
	  p -> bCond[i].value.given = val;
	}
      }

    } else if (strstr (token, "WALL")) {

      for (i = 0; i < nvar; i++) {
	p -> bCond[i].kind        = WALL;
	p -> bCond[i].value.given = 0.0;
      }
    
    } else if (strstr (token, "OUTFLOW")) {

      for (i = 0; i < nvar; i++) {
	p -> bCond[i].kind        = OUTFLOW;
	p -> bCond[i].value.given = 0.0;
      }

    } else
      message (routine, strcat (token, "?: unrecognized BC type"), ERROR);

    if (BCList) p -> next = BCList;
    BCList = p;
  }

}





void  printBCs ()
/* ========================================================================= *
 * Utility to print list of tagged boundary conditions.                      *
 * ========================================================================= */
{
  char  routine[] = "printBCs";
  int   i, nvar = iparam ("NVAR");
  B_C  *p;


  if (! BCList) message (routine, "empty list of BCs", WARNING);

  for (p = BCList; p; p = p -> next) {

    printf ("id %1d, names: %s\n", p -> id, p -> name);

    for (i = 0; i<nvar; i++) {

      switch (p -> bCond[i].kind) {

      case ESSENTIAL:
      default:
	printf ("  ESSENTIAL:    %g\n", p -> bCond[i].value.given );
	break;

      case ESSENTIAL_FN:
        printf ("  ESSENTIAL_FN: %s\n", p -> bCond[i].value.interp);
        break;

      case NATURAL:
	printf ("  NATURAL:      %g\n", p -> bCond[i].value.given );
	break;

      case NATURAL_FN:
        printf ("  NATURAL_FN:   %s\n", p -> bCond[i].value.interp);
        break;

      case WALL:
	printf ("  WALL\n");
	break;

      case OUTFLOW:
	printf ("  OUTFLOW\n");
	break;
      }
    }
  }
}





int  countBCs ()
/* ========================================================================= *
 * Count number of tagged BCs.                                               *
 * ========================================================================= */
{
  B_C *p;
  int  num = 0;

  for (p = BCList; p ; p = p -> next) num++;
  return num;
}





BcRecord  *BCSearch (int id, char name)
/* ========================================================================= *
 * Search list for matching tag & return boundary condition record.          *
 * ========================================================================= */
{
  char   routine[] = "BCSearch";
  B_C   *p =  NULL;


  for (p = BCList; p; p = p -> next) 

    if (p -> id == id) {

      if (!(strchr (p -> name, name))) {
	sprintf (buf, "variable \"%c\" not found in \"%s\"", name, p -> name); 
	message (routine, buf, ERROR);
      }

      switch (name) {
      case 'u': return p -> bCond    ;
      case 'v': return p -> bCond + 1;
      case 'w': return p -> bCond + 2;
      case 'p': return 0;
      default :
	message (routine, "blank name character ?", ERROR);
	break;
      }

    }
  
  sprintf (buf, "no boundary condition set for tag %1d", id);
  message (routine, buf, ERROR);
  return  0;
}





int freeBoundary (const Bedge *B)
/* ========================================================================= *
 * Return true if all boundary conditions are OUTFLOW or NATURAL type.       *
 * Helmholtz problem: if all boundary conditions are type OUTFLOW, the       *
 * problem is underdetermined, if also lambda2 = 0.0.                        *
 *                                                                           *
 * The response (carried out elsewhere) will be to arbitrarily assign one of *
 * the nodal values, and reduce the order of the problem by one, thus making *
 * the problem determinate.                                                  *
 *                                                                           *
 * ========================================================================= */
{
  for (; B; B = B -> next)
    if (   B->bc->kind == ESSENTIAL
	|| B->bc->kind == ESSENTIAL_FN
	|| B->bc->kind == WALL        )
      return 0;

  return 1;
}





void  evaluateBC (Bedge *B)
/* ========================================================================= *
 * Load boundary condition storage areas with numeric values.                *
 *                                                                           *
 * This routine is not called for pressure BCs (evaluatePBC instead).        *
 * ========================================================================= */
{
  int    np;


  for (; B; B = B -> next) {
    np = B -> np;

    switch (B -> bc -> kind) {

    case ESSENTIAL:
    case NATURAL:
      Veclib::fill (np, B -> bc -> value.given, B -> value, 1);
      break;

    case ESSENTIAL_FN:
    case NATURAL_FN:
      {
	double *x, *y;
	assert ((x = dvector (0, np-1)) != 0);
	assert ((y = dvector (0, np-1)) != 0);

	Veclib::copy (np, *B -> elmt -> xmesh + B -> estart, B -> eskip, x, 1);
	Veclib::copy (np, *B -> elmt -> ymesh + B -> estart, B -> eskip, y, 1);

	vecInit   ("x y", B -> bc -> value.interp);
	vecInterp (np, x, y, B -> value);

	freeDvector (x, 0);
	freeDvector (y, 0);
      }
      break;

    case OUTFLOW: 
    case WALL:
    default:
      Veclib::fill (np, 0.0, B -> value, 1);
      break;
    }
  }
}
