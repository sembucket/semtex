%{
/*****************************************************************************
 * INITIAL: yacc code for maintaining lists used for parameter and flag
 * lookup, along with a simple function interpreter.  Initialize() routine
 * must be called before the other parts of the code will work.
 *
 * Modelled on hoc3 in Chapter 8 of "The UNIX Programming Environment" by
 * Kernighan & Pike.
 *
 * We maintain 3 externally-accessible lists:
 *   iparam:  named integer parameters;
 *   oparam:  named integer options (flags);
 *   dparam:  named double parameters.
 * The lists get preloaded with default values and useful constants by
 * initial().  Once initialized, the lists are accessible only through calls
 * to routines located in this file.
 *
 * In addition a list of symbols, symlist, is maintained for use by the
 * function interpreter.  This is externally accessible through the routine
 * interpret(), which returns double.  Operations include single-argument
 * functions, unary minus and the binary operators ^ (exponentiation) and
 * ~ (atan2), and the Heaviside function heav.
 *
 * Summary of externally-accessible functions:
 *
 * void   initialize (void);
 * double interpret  (const char *);
 *
 * void   vecInit    (const char *, const char *);
 * void   vecInterp  (int   , ...   );
 *
 * void   setOption  (const char *, int);
 * int    option     (const char *);
 * void   showOption (void);
 *
 * void   setIparam  (const char *, int);
 * int    iparam     (const char *);
 * void   showIparam (void);
 *
 * void   setDparam  (const char *, double);
 * double dparam     (const char *);
 * void   showDparam (void);
 *
 * $Id$
 *****************************************************************************/


#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <errno.h>
extern int errno;

#include <alplib.h>

int yyparse (void);


/* ------------------------------------------------------------------------- *
 * File-scope type definitions.
 * ------------------------------------------------------------------------- */

typedef double (*PFD) (double);	/* Pointer to function returning double */

typedef struct symbol {		/* Symbol table entry */
  char *name;
  int   type;			/* VAR, BLTIN, UNDEF, DPARAM, IPARAM, OPTION */
  union {
    double dval;		/* If VAR, DPARAM    */
    int    ival;		/* If IPARAM, OPTION */
    PFD    ptr;			/* If BLTIN          */
  } u;
  struct symbol *next;
} Symbol;

    
/* ------------------------------------------------------------------------- *
 * File-scope prototypes.
 * ------------------------------------------------------------------------- */

static double  Log  (double),          Log10   (double),
               Exp  (double),          Sqrt    (double),
               Asin (double),          Acos    (double),
               Pow  (double, double),  integer (double),
               Sinh (double),          Cosh    (double),
               Tanh (double),          Heavi   (double);
						     

static Symbol *lookup  (const char *);
static Symbol *install (const char *, int, ...);
static void   *emalloc (size_t);


/* ------------------------------------------------------------------------- *
 * File-scope variables.
 * ------------------------------------------------------------------------- */

static Symbol *symlist = NULL;      /* Internal use for function interpreter */
static Symbol *dlist   = NULL;      /* Storage of double lookup parameters   */
static Symbol *ilist   = NULL;      /* Storage of integer lookup parameters  */
static Symbol *olist   = NULL;      /* Storage for option lookup             */

static char   func_string[STR_MAX], *cur_string;
static double value;

#define VEC_MAX 16
static int     nvec = 0;
static Symbol *vs[VEC_MAX];	    /* Storage for vector parser variables   */

static struct {			    /* Built-in functions                    */
  char *name;
  PFD   func;
} builtin[] = {
  "sin"   ,  sin     ,
  "cos"   ,  cos     ,
  "tan"   ,  tan     ,
  "atan"  ,  atan    ,
  "int"   ,  integer ,
  "abs"   ,  fabs    ,
  "floor" ,  floor   ,
  "ceil"  ,  ceil    ,
  "heav"  ,  Heavi   ,
  "asin"  ,  Asin    ,		    /* Rest do error-checking */
  "acos"  ,  Acos    ,
  "log"   ,  Log     ,
  "log10" ,  Log10   ,
  "exp"   ,  Exp     ,
  "sqrt"  ,  Sqrt    ,
  "sinh"  ,  Sinh    ,
  "cosh"  ,  Cosh    ,
  "tanh"  ,  Tanh    ,
  NULL    ,  NULL
};

#include "defaults.h"

%}
/* ------------------------------------------------------------------------- *
 * Yacc grammar follows.
 * ------------------------------------------------------------------------- */
%union {			/* yacc stack type      */
  double  val;			/* actual value         */
  Symbol *sym;			/* symbol table pointer */
}
%token <val>   NUMBER
%token <sym>   VAR BLTIN UNDEF DPARAM IPARAM OPTION
%type  <val>   expr asgn
%right '='
%left  '+' '-'
%left  '*' '/'
%left  UNARYMINUS
%right '^' '~'			/* exponentiation, atan2 */
%%
list:     /* nothing */
        | list '\n'
        | list asgn '\n'
        | list expr '\n'     { value = $2; }
        ;
asgn:     VAR '=' expr       { $$=$1->u.dval=$3; $1->type = VAR; }
        ;
expr:     NUMBER
        | VAR                { if ($1->type == UNDEF) {
				 message("in yyparse(): undefined variable ",
					 $1->name, WARNING);
			       }
			       $$ = $1->u.dval;
			     }
        | asgn
        | BLTIN '(' expr ')' { $$ = (*($1->u.ptr))($3); }
        | expr '+' expr      { $$ = $1 + $3; }
        | expr '-' expr      { $$ = $1 - $3; }
        | expr '*' expr      { $$ = $1 * $3; }
        | expr '/' expr      { if ($3 == 0.0)
				 message("in yyparse()",
					 " division by zero.", ERROR);
			       $$ = $1 / $3;
			     }
        | expr '^' expr      { $$ = Pow($1, $3); }
        | expr '~' expr      { $$ = atan2($1, $3); }
        | '(' expr ')'       { $$ = $2; }
        | '-' expr %prec UNARYMINUS { $$ = -$2; }
        ;
%%


/*****************************************************************************
 * Externally-visible functions follow.
 *****************************************************************************/


void initialize (void)
/* ========================================================================= *
 * Load lookup tables and symbol table with default values.
 *
 * This routine should be called at start of run-time.
 * ========================================================================= */
{
  int i;

#ifdef __sgi			/* Enable floating-point traps. */
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

  for (i=0; consts[i].name; i++)
    install (consts[i].name, DPARAM, consts[i].cval);

  for (i=0; option_init[i].name; i++) 
    install (option_init[i].name, OPTION, option_init[i].oval);

  for (i=0; iparam_init[i].name; i++)
    install (iparam_init[i].name, IPARAM, iparam_init[i].ival);

  for (i=0; dparam_init[i].name; i++)
    install (dparam_init[i].name, DPARAM, dparam_init[i].dval);

  for (i=0; builtin[i].name; i++)
    install (builtin[i].name, BLTIN, builtin[i].func);
}


double interpret (const char *s)
/* ========================================================================= *
 * Given a string, interpret it as a function using yacc-generated yyparse().
 * ========================================================================= */
{
  if (strlen (s) > STR_MAX)
    message ("in interpret(): too many characters passed\n", s, ERROR);
  
  strcat (strcpy (func_string, s), "\n");
  
  cur_string = func_string;
  
  yyparse ();
  return value;
}


void vecInit (const char *names, const char *fn)
/* ========================================================================= *
 * Set up the vector parser.
 *
 * names contains a list of variable names  e.g. "x y z",
 * fn    contains a function for evaluation e.g. "sin(x)*cos(y)*exp(z)".
 *
 * Valid separator characters in name are space, tab, comma, (semi-)colon.
 * Function string can contain previously-defined symbols (e.g. PI).
 * ========================================================================= */
{
  char    routine   [] = "vecInit()";
  char    separator [] = " ,:;\t";
  char    tmp       [STR_MAX];
  char   *p;
  Symbol *s;

  if (strlen (fn) == 0)
    message (routine, "empty function string", ERROR);
  else if (strlen (fn) > STR_MAX)
    message (routine, "too many characters in function string", ERROR);

  nvec = 0;
  strcpy (tmp, names);

  p = strtok (tmp, separator);
  do {
    if (nvec++ > VEC_MAX) message (routine, "too many variables", ERROR); 
    vs[nvec-1] = (s=lookup(p)) ? s : install (p, VAR, 0.0);
  } while (p = strtok (NULL, separator));

  
  strcat (strcpy(func_string, fn), "\n");
}


void vecInterp (int ntot, ...)
/* ========================================================================= *
 * Vector parser.  Following ntot there should be passed a number of
 * pointers to double (vectors), of which there should be in number the
 * number of variables named previously to vecInit, plus one: the result of
 * continually re-parsing the string "fn" is placed in the last vector, for
 * a total of ntot parsings.
 *
 * To follow on from the previous example, four vectors would be passed,
 * i.e.  vecInterp(ntot, x, y, z, u); the result fn(x,y,z) is placed in u.
 *
 * ========================================================================= */
{
  char       routine[] = "vecInterp()";
  register   i, n;
  double    *x[VEC_MAX];
  double    *fx = NULL;
  va_list    ap;
  
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
    for (i = 0; i < nvec; i++) vs[i]->u.dval = x[i][n];
    yyparse ();
    fx[n] = value;
  }
}


void setOption (const char *s, int v)
/* ========================================================================= *
 * Set option on list true/false, or install it.
 * ========================================================================= */
{
  Symbol *sp;
 
  for (sp = olist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0) {
      sp->u.ival = v;
      break;
    }
  
  if (!sp) sp = install (s, OPTION, v);
}


int option (const char *s)
/* ========================================================================= *
 * Retrieve value from option list.
 * ========================================================================= */
{
  Symbol *sp;
 
  for (sp = olist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0)
      return sp->u.ival;

  message ("option(): name not found", s, WARNING);
  return 0;
}
    

void showOption (void)
/* ========================================================================= *
 * Echo option list to stdout.
 * ========================================================================= */
{
  Symbol *sp;

  for (sp = olist; sp; sp = sp->next)
    printf ("%-12s%-3d\n", sp->name, sp->u.ival);
}


void setIparam (const char *s, int v)
/* ========================================================================= *
 * Set or install iparam on list.
 * ========================================================================= */
{
  Symbol *sp;

  for (sp = ilist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0) {
      sp->u.ival = v;
      break;
    }
  
  if (!sp) sp = install (s, IPARAM, v);
}


int iparam (const char *s)
/* ========================================================================= *
 * Retrieve value from iparam list.
 * ========================================================================= */
{
  Symbol *sp;
 
  for (sp = ilist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0)
      return sp->u.ival;

  message ("iparam(): name not found", s, WARNING);
  return 0;
}


void showIparam (void)
/* ========================================================================= *
 * Echo iparam list to stdout.
 * ========================================================================= */
{
  Symbol *sp;
  
  for (sp = ilist; sp; sp = sp->next)
    printf ("%-12s%-3d\n", sp->name, sp->u.ival);
}


void setDparam (const char *s, double v)
/* ========================================================================= *
 * Set or install dparam on list.
 * ========================================================================= */
{
  Symbol *sp;

  for (sp = dlist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0) {
      sp->u.dval = v;
      break;
    }
  
  for (sp = symlist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0) {
      sp->u.dval = v;
      break;
    }
  
  if (!sp) sp = install (s, DPARAM, v);
}


double dparam (const char *s)
/* ========================================================================= *
 * Retrieve value from dparam list.
 * ========================================================================= */
{
  Symbol *sp;
  
  for (sp = dlist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0)
      return sp->u.dval;

  message ("dparam(): name not found", s, WARNING);
  return 0;
}


void showDparam (void)
/* ========================================================================= *
 * Echo dparam list to stdout.
 * ========================================================================= */
{
  Symbol *sp;
  
  for (sp = dlist; sp; sp = sp->next)
    printf ("%-12s%-.6g\n", sp->name, sp->u.dval);
}


/* ##################### FILE-SCOPE ROUTINES FOLLOW ######################## */


static int yylex (void)
/* ========================================================================= *
 * Lexical analysis routine called by yyparse, using string loaded by
 * interpret().
 * ========================================================================= */
{
  int c;

  while ((c = *cur_string++) == ' ' || c == '\t');

  if (c == EOF)
    return 0;
  if (c == '.' || isdigit (c)) {	/* number */
    yylval.val = strtod (--cur_string, &cur_string);
    return NUMBER;
  }
  if (isalpha (c)) {		/* symbol */
    Symbol *s;
    char    sbuf[STR_MAX], *p = sbuf;
    do {
      *p++ = c;
    } while ((c = *cur_string++) != EOF && (isalnum (c) || c == '_'));
    cur_string--;
    *p = '\0';
    if ((s = lookup (sbuf)) == NULL)
      s = install (sbuf, UNDEF, 0.0);
    yylval.sym = s;
    return (s->type == UNDEF || s->type == DPARAM) ? VAR : s->type;
  }
  return c;
}


static Symbol *lookup(const char *s)
/* ========================================================================= *
 * Find s in symbol table.
 * ========================================================================= */
{
  Symbol *sp;
  
  for (sp = symlist; sp; sp = sp->next)
    if (strcmp (sp->name, s) == 0)
      return sp;

  return NULL;
}


static Symbol *install(const char *s, int t, ...)
/* ========================================================================= *
 * Install value of variable type into appropriate lists.
 * Note that variables of type DPARAM get mounted in both symlist and dlist.
 * ========================================================================= */
{
  Symbol *sp;
  va_list ap;
  
  va_start (ap, t);

  sp       = (Symbol *) emalloc (sizeof(Symbol));
  sp->name = emalloc (strlen(s)+1);
  strcpy (sp->name, s);

  switch (sp->type = t) {
  case OPTION:
    sp->u.ival = va_arg (ap, int);
    sp->next   = olist;
    olist      = sp;
    break;
  case IPARAM:
    sp->u.ival = va_arg (ap, int);
    sp->next   = ilist;
    ilist      = sp;
    break;
  case DPARAM:
    sp->u.dval = va_arg (ap, double);
    sp->next   = dlist;
    dlist      = sp;
    sp         = (Symbol *) emalloc (sizeof(Symbol));
    memcpy (sp, dlist, sizeof (Symbol));
    sp->name = emalloc (strlen (s)+1);
    strcpy (sp->name, s);
    sp->next   = symlist;
    symlist    = sp;
    break;
  case VAR:
  case UNDEF:
    sp->u.dval = va_arg (ap, double);
    sp->next   = symlist;
    symlist    = sp;
    break;
  case BLTIN:
    sp->u.ptr  = va_arg (ap, PFD);
    sp->next   = symlist;
    symlist    = sp;
    break;
  default:
    sp->u.dval = va_arg (ap, double);
    sp->type   = UNDEF;
    sp->next   = symlist;
    symlist    = sp;
    break;
  }

  va_end (ap);

  return sp;
}


static void *emalloc (size_t n)
/* ========================================================================= *
 * Check return from malloc().
 * ========================================================================= */
{
  void *p;

  if (!(p = (void *) malloc (n)))
    message("emalloc()", "out of memory", ERROR);
 
  return p;
}


static void yyerror (char *s)
/* ========================================================================= *
 * Handler for yyparse() syntax errors.
 * ========================================================================= */
{
  message("yyparse()", s, WARNING);
}


static double errcheck (double d, char *s)
/* ========================================================================= *
 * Check result of math library call.
 * ========================================================================= */
{
  if (errno == EDOM) {
    errno = 0;
    message ("errcheck(): argument out of domain in call to", s, ERROR);
  } else if (errno == ERANGE) {
    errno = 0;
    message ("errcheck(): result out of range in call to", s, ERROR);
  }

  return d;
}


/*****************************************************************************
 * Remaining routines do error-checking calls to math library routines.
 *****************************************************************************/


static double Heavi (double x)
{
  return (x >= 0.0) ? 1.0 : 0.0;
}


static double Log (double x)
{
  return errcheck (log(x), "log");
}


static double Log10 (double x)
{
  return errcheck (log10(x), "log10");
}


static double Exp (double x)
{
  return errcheck (exp(x), "exp");
}


static double Sqrt (double x)
{
  return errcheck (sqrt(x), "sqrt");
}


static double Pow (double x, double y)
{
  return errcheck (pow(x, y), "exponentiation");
}


static double Acos (double x)
{
  return errcheck (acos(x), "acos");
}



static double Asin (double x)
{
  return errcheck (asin(x), "asin");
}


static double integer (double x)
{
  return (double) (long)x;
}


static double Sinh (double x)
{
  return errcheck (sinh(x), "sinh");
}


static double Cosh (double x)
{
  return errcheck (cosh(x), "cosh");
}


static double Tanh (double x)
{
  return errcheck (tanh(x), "tanh");
}
