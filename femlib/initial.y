%{
/*****************************************************************************
 * Synopsis
 * --------
 * initial.y: yacc code for a simple function interpreter, which also allows
 * lookup for named tokens (all are double precision).
 *
 * Modelled on "hoc3" in Chapter 8 of "The UNIX Programming Environment"
 * by Kernighan & Pike, Prentice-Hall 1984.  A hash-table data structure is
 * used in place of the linked list of Symbols which they employed.
 * Hashing code from Chapter 6 of "The C Programming Language", 2nd Edn,
 * by Kernighan & Ritchie, Prentice-Hall 1988.
 *
 * Summary
 * -------
 * void   yy_initialize (void);
 * void   yy_help       (void);
 * void   yy_show       (void);
 * int    yy_dump       (char*, const int);
 *
 * double yy_interpret  (const char*);
 *
 * void   yy_vec_init   (const char*, const char*);
 * void   yy_vec_interp (const int, ...);
 *
 * Notes
 * -----
 * 1. yy_initialize must be called before other routines will work.
 * 2. yy_help prints a summary of available functions on stdout.
 * 3. yy_show prints a summary of installed variables on stdout.
 * 4. yy_dump is similar to yy_show, but prints into a string of given length.
 * 5. yy_interpret is the central routine.  It can be used:
 *    (a), to install a named symbol in internal tables: e.g. "name = value";
 *    (b), to retrieve the value of an installed symbol: e.g. "name";
 *    (c), to evaluate a function string: e.g. "cos(x)*exp(-t)".
 * 6. yy_vec_init is used to set up the interpreter for vector evaluation.
 * 7. yy_vec_interp subsequently used for "vectorized" calls to yy_interpret.
 *
 * Operators
 * ---------
 * Unary:     -
 * Binary:    -, +, *, /, ^ (exponentiation), ~ (atan2), % (fmod)
 * Functions: sin,  cos,  tan,  abs, floor, ceil, int,   heav (Heaviside),
 *            asin, acos, atan, log, log10, exp,  sqrt,  lgamma,
 *            sinh, cosh, tanh, erf, erfc,  j0,   j1,    y0,  y1
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>

#include <femdef.h>
#include <veclib.h>

#if 1
#define HASHSIZE 199
#else
#define HASHSIZE 37
#endif
#define HASHSEED 31
#define VEC_MAX  32

typedef double (*PFD) (double);

typedef struct symbol {
  char* name;
  short type;
  union {
    double val;			/* -- If VAR.   */
    PFD    ptr;			/* -- If BLTIN. */
  } u;
  struct symbol* next;
} Symbol;

static double  Log  (double),  Log10    (double),
               Asin (double),  Acos     (double),
               Sqrt (double),  Pow      (double, double),
               Heavi(double),  errcheck (const double, const char*);

static unsigned hash     (const char*);
static Symbol*  lookup   (const char*);
static Symbol*  install  (const char*, const int, const double);
static void*    emalloc  (const size_t);

       int      yyparse (void);
static int      yylex   (void);
static void     yyerror (char*);

static double  value;
static Symbol* hashtab[HASHSIZE];
static char    func_string[STR_MAX], *cur_string;
static int     nvec = 0;
static Symbol* vs[VEC_MAX];
extern int     errno;

static struct {			    /* -- Built-in functions. */
  char* name;
  PFD   func;
} builtin[] = {
  "sin"   ,  sin     ,
  "cos"   ,  cos     ,
  "tan"   ,  tan     ,
  "atan"  ,  atan    ,
  "exp"   ,  exp     ,
  "sinh"  ,  sinh    ,
  "cosh"  ,  cosh    ,
  "tanh"  ,  tanh    ,
  "erf"   ,  erf     ,
  "erfc"  ,  erfc    ,
  "int"   ,  rint    ,
  "abs"   ,  fabs    ,
  "floor" ,  floor   ,
  "ceil"  ,  ceil    ,
  "asin"  ,  asin    ,
  "acos"  ,  acos    ,
  "log"   ,  log     ,
  "log10" ,  log10   ,
  "sqrt"  ,  sqrt    ,
  "heav"  ,  Heavi   ,
  "j0"    ,  j0      ,
  "j1"    ,  j1      ,
  "y0"    ,  y0      ,
  "y1"    ,  y1      ,
  "lgamma",  lgamma  ,

  NULL    ,  NULL
};

#include "defaults.h"

%}
/* -- Yacc grammar follows: */
%union {			/* -- yacc stack type      */
  double  val;			/* -- actual value         */
  Symbol* sym;			/* -- symbol table pointer */
}
%token <val>   NUMBER
%token <sym>   VAR BLTIN UNDEF
%type  <val>   expr asgn
%right '='
%left  '+' '-'
%left  '*' '/'
%left  UNARYMINUS
%right '^' '~'			/* -- exponentiation, atan2 */
%%
list:     /* nothing */
        | list '\n'
        | list asgn '\n'
        | list expr '\n'     { value = $2; }
        ;
asgn:     VAR '=' expr       { $$=$1->u.val=$3; $1->type = VAR; }
        ;
expr:     NUMBER
        | VAR                { if ($1->type == UNDEF) {
				 message ("yyparse: undefined variable ",
					  $1->name, WARNING);
			       }
			       $$ = $1->u.val;
			     }
        | asgn
        | BLTIN '(' expr ')' { $$ = (*($1->u.ptr))($3); }
        | expr '+' expr      { $$ = $1 + $3; }
        | expr '-' expr      { $$ = $1 - $3; }
        | expr '*' expr      { $$ = $1 * $3; }
        | expr '/' expr      { if ($3 == 0.0) 
				 message ("yyparse",
					  "division by zero", ERROR);
			       $$ = $1 / $3;
			     }
        | expr '^' expr      { $$ = pow   ($1, $3); }
        | expr '~' expr      { $$ = atan2 ($1, $3); }
        | expr '%' expr      { $$ = fmod  ($1, $3); }
        | '(' expr ')'       { $$ = $2; }
        | '-' expr %prec UNARYMINUS { $$ = -$2; }
        ;
%%


void yy_initialize (void)
/* ------------------------------------------------------------------------- *
 * Load lookup tables and symbol table with default values.
 *
 * This routine should be called at start of run-time.
 * ------------------------------------------------------------------------- */
{
  static   int initialized = 0;
  register int i;
  register Symbol* s;

#if 0			/* Enable SGI floating-point traps. */
#include <sigfpe.h>
#include <sgidefs.h>     
#ifndef _MIPS_SIM_ABI64         /* No 64-bit version of libfpe.a yet. */
  sigfpe_[ _OVERFL].trace = 1;
  sigfpe_[ _OVERFL].exit  = 1;

  sigfpe_[_DIVZERO].trace = 1;
  sigfpe_[_DIVZERO].exit  = 1;

  sigfpe_[_INVALID].trace = 1;
  sigfpe_[_INVALID].exit  = 1;

  sigfpe_[_UNDERFL].repls = _ZERO;

  handle_sigfpes (_ON,
		  _EN_OVERFL | _EN_DIVZERO | _EN_INVALID | _EN_UNDERFL,
		  0,
		  _ABORT_ON_ERROR,
		  0);
#endif
#endif

  if (!initialized) {
    for (i = 0; consts[i].name; i++)
      install (consts[i].name, VAR, consts[i].cval);

    for (i = 0; builtin[i].name; i++) {
      s = install (builtin[i].name, BLTIN, 0.0);
      s -> u.ptr = builtin[i].func;
    }

    initialized = 1;
  }
}


double yy_interpret (const char* s)
/* ------------------------------------------------------------------------- *
 * Given a string, interpret it as a function using yacc-generated yyparse.
 * ------------------------------------------------------------------------- */
{
  if (strlen (s) > STR_MAX)
    message ("yy_interpret: too many characters passed:\n", s, ERROR);
  
  strcat (strcpy (func_string, s), "\n");
  
  cur_string = func_string;
  
  yyparse ();
  return value;
}


void yy_vec_init (const char* names,
		  const char* fn   )
/* ------------------------------------------------------------------------- *
 * Set up the vector parser.
 *
 * names contains a list of variable names  e.g. "x y z",
 * fn    contains a function for evaluation e.g. "sin(x)*cos(y)*exp(z)".
 *
 * Valid separator characters in name are space, tab, comma, (semi-)colon.
 * Function string can contain previously-defined symbols (e.g. PI).
 * ------------------------------------------------------------------------- */
{
  char    routine   [] = "vecInit()";
  char    separator [] = " ,:;\t";
  char    tmp       [STR_MAX];
  char*   p;
  Symbol* s;

  if (strlen (fn) == 0)
    message (routine, "empty function string", ERROR);
  else if (strlen (fn) > STR_MAX)
    message (routine, "too many characters in function string", ERROR);

  nvec = 0;
  strcpy (tmp, names);

  p = strtok (tmp, separator);
  do {
    if (nvec++ > VEC_MAX) message (routine, "too many variables", ERROR); 
    vs[nvec-1] = (s = lookup (p)) ? s : install (p, VAR, 0.0);
  } while (p = strtok (NULL, separator));

  
  strcat (strcpy (func_string, fn), "\n");
}


void yy_vec_interp (const int ntot, ...)
/* ------------------------------------------------------------------------- *
 * Vector parser.  Following ntot there should be passed a number of
 * pointers to double (vectors), of which there should be in number the
 * number of variables named previously to vecInit, plus one: the result of
 * continually re-parsing the string "fn" is placed in the last vector, for
 * a total of ntot parsings.
 *
 * To follow on from the previous example, four vectors would be passed,
 * i.e.  vecInterp(ntot, x, y, z, u); the result fn(x,y,z) is placed in u.
 * ------------------------------------------------------------------------- */
{
  char         routine[] = "yy_vec_interp";
  register int i, n;
  double*      x[VEC_MAX];
  double*      fx = NULL;
  va_list      ap;
  
  va_start (ap, ntot);
  for (i = 0; i < nvec; i++) {
    x[i] = NULL;
    if (!(x[i] = va_arg (ap, double*)))
	message (routine, "not enough vectors passed..1", ERROR);
  }
  if (!(fx = va_arg (ap, double*)))
    message (routine, "not enough vectors passed..2", ERROR);
  va_end (ap);

  for (n = 0; n < ntot; n++) {
    cur_string = func_string;
    for (i = 0; i < nvec; i++) vs[i]->u.val = x[i][n];
    yyparse ();
    fx[n] = value;
  }
}


void yy_help (void)
/* ------------------------------------------------------------------------- *
 * Print details of callable functions to stderr.
 * ------------------------------------------------------------------------- */
{
  fprintf 
    (stderr, 
     "Unary    : -\n"
     "Binary   : -, +, *, /, ^ (exponentiation), ~ (atan2), %% (fmod)\n"
     "Functions: sin,  cos,  tan,  abs, floor, ceil, int,  heav (Heaviside),\n"
     "           asin, acos, atan, log, log10, exp,  sqrt, lgamma,\n"
     "           sinh, cosh, tanh, erf, erfc,  j0,   j1,   y0,  y1\n");
}


void yy_show (void)
/* ------------------------------------------------------------------------- *
 * Print details of installed variables to stderr.
 * ------------------------------------------------------------------------- */
{
  register int     i;
  register Symbol* sp;

  for (i = 0; i < HASHSIZE; i++)
    for (sp = hashtab[i]; sp; sp = sp -> next)
      if (sp -> type == VAR) 
	fprintf (stderr, "%-15s = %-.17g\n", sp -> name, sp -> u.val);
}


int yydump (char*     str,
	    const int max)
/* ------------------------------------------------------------------------- *
 * Load description of internal variables into string, to length max.
 * If string overflows, return 0, else 1.
 * ------------------------------------------------------------------------- */
{
  register int     i, n = 0;
  register Symbol* sp;
  char             buf[FILENAME_MAX];

  for (i = 0; i < HASHSIZE; i++)
    for (sp = hashtab[i]; sp; sp = sp -> next)
      if (sp -> type == VAR) {
	sprintf (buf, "%-15s = %-.17g\n", sp -> name, sp -> u.val);
	if ((n += strlen (buf)) > max - 2)
	  return 0;
	else
	  strcat (str, buf);
      }

  return 1;
}


static int yylex (void)
/* ------------------------------------------------------------------------- *
 * Lexical analysis routine called by yyparse, using string loaded by
 * yy_interpret.
 * ------------------------------------------------------------------------- */
{
  register int c;

  while ((c = *cur_string++) == ' ' || c == '\t');

  if (c == EOF) return 0;

  if (c == '.' || isdigit (c)) {
    yylval.val = strtod (--cur_string, &cur_string);
    return NUMBER;
  }

  if (isalpha (c)) {
    register Symbol* s;
    char             sbuf[STR_MAX];
    register char*   p = sbuf;
    do {
      *p++ = c;
    } while ((c = *cur_string++) != EOF && (isalnum (c) || c == '_'));
    cur_string--;
    *p = '\0';
    if ((s = lookup (sbuf)) == NULL) s = install (sbuf, UNDEF, 0.0);
    yylval.sym = s;
    return (s -> type == UNDEF) ? VAR : s -> type;
  }

  return c;
}


static unsigned hash (const char* s)
/* ------------------------------------------------------------------------- *
 * Generate hash table index.
 * ------------------------------------------------------------------------- */
{
  register unsigned hashval;

  for (hashval = 0; *s != '\0'; s++) hashval = *s + HASHSEED * hashval;
  
  return hashval % HASHSIZE;
}


static Symbol* lookup (const char* s)
/* ------------------------------------------------------------------------- *
 * Find s in symbol hashtable.
 * ------------------------------------------------------------------------- */
{
  register Symbol* sp;
  
  for (sp = hashtab[hash (s)]; sp; sp = sp->next)
    if (strcmp (s, sp->name) == 0) return sp;

  return NULL;
}


static Symbol* install (const char*  s,
			const int    t,
			const double d)
/* ------------------------------------------------------------------------- *
 * Install s in symbol hashtable.
 * ------------------------------------------------------------------------- */
{
  register Symbol*  sp;
  register unsigned hashval;

  if (!(sp = lookup (s))) {	/* -- Not found, install in hashtab. */
    sp = (Symbol *) emalloc (sizeof (Symbol));
    if (sp == NULL || (sp -> name = strdup (s)) == NULL) return NULL;
    hashval          = hash (s);
    sp -> next       = hashtab[hashval];
    hashtab[hashval] = sp;
  }

  sp -> type  = t;
  sp -> u.val = d;
  
  return sp;
}


static void *emalloc (const size_t n)
/* ------------------------------------------------------------------------- *
 * Check return from malloc.
 * ------------------------------------------------------------------------- */
{
  void* p;

  if (!(p = (void *) malloc (n))) message ("emalloc", "out of memory", ERROR);
 
  return p;
}


static void yyerror (char *s)
/* ------------------------------------------------------------------------- *
 * Handler for yyparse syntax errors.
 * ------------------------------------------------------------------------- */
{
  message ("yyparse", s, WARNING);
}


static double errcheck (const double d,
			const char*  s)
/* ------------------------------------------------------------------------- *
 * Check result of math library call.
 * ------------------------------------------------------------------------- */
{
  if (errno == EDOM) {
    errno = 0;
    message ("errcheck: argument out of domain in call to", s, ERROR);
  } else if (errno == ERANGE) {
    errno = 0;
    message ("errcheck: result out of range in call to", s, ERROR);
  }

  return d;
}


static double Heavi (double x) { return (x >= 0.0) ? 1.0 : 0.0; }

#if 0
/* ------------------------------------------------------------------------- *
 * Error checking is done for cases where we can have illegal input
 * values (typically negative), otherwise we accept exception returns.
 *
 * After much trouble with underflow errors, I've turned these off for now.
 * HMB -- 1/6/2001.
 * ------------------------------------------------------------------------- */

static double Log (double x)
{ return errcheck (log (x), "log"); }


static double Log10 (double x)
{ return errcheck (log10 (x), "log10"); }


static double Sqrt (double x)
{ return errcheck (sqrt (x), "sqrt"); }


static double Pow (double x, double y)
{ return errcheck (pow (x, y), "exponentiation"); }


static double Acos (double x)
{ return errcheck (acos (x), "acos"); }


static double Asin (double x)
{ return errcheck (asin (x), "asin"); }
#endif
