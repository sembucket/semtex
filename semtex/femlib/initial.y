%{
/*****************************************************************************
 * INITIAL: yacc code for maintaining lists used for parameter and flag      *
 * lookup, along with a simple function interpreter.  Initialize() routine   *
 * must be called before the other parts of the code will work.              *
 *                                                                           *
 * Modelled on hoc3 in Chapter 8 of "The UNIX Programming Environment".      *
 * A(nother) good idea of Ron Henderson's...thanks, Ron.                     *
 *                                                                           *
 * We maintain 3 externally-accessible lists:                                *
 *   iparam:  named integer parameters;                                      *
 *   oparam:  named integer options (flags);                                 *
 *   dparam:  named double parameters.                                       *
 * The lists get preloaded with default values and useful constants by       *
 * initial().  Once initialized, the lists are accessible only through calls *
 * to routines located in this file.                                         *
 *                                                                           *
 * In addition a list of symbols, symlist, is maintained for use by the      *
 * function interpreter.  This is externally accessible through the routine  *
 * interpret(), which returns double.  Operations include single-argument    *
 * functions, unary minus and the binary operators ^ (exponentiation) and    *
 * ~ (atan2). Flow control, procedures and functions are left for revisions. *
 *                                                                           *
 * Summary of externally-accessible functions:                               *
 *                                                                           *
 * void   initialize  (void);                                                *
 * double interpret   (char *);                                              *
 *                                                                           *
 * void   set_option  (char *, int);                                         *
 * int    get_option  (char *);                                              *
 * void   show_option (FILE *);                                              *
 *                                                                           *
 * void   set_iparam  (char *, int);                                         *
 * int    get_iparam  (char *);                                              *
 * void   show_iparam (FILE *);                                              *
 *                                                                           *
 * void   set_dparam  (char *, double);                                      *
 * double get_dparam  (char *);                                              *
 * void   show_dparam (FILE *);                                              *
 *                                                                           *
 *****************************************************************************/

/*------------------*   
 * RCS Information: *
 *------------------*/
static char
  rcsid00[] = "$Id$";


#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <errno.h>

#include <linalp.h>

extern int errno;

#ifdef _mips			/* We have a MIPS floating-point accelerator */
  #include <sigfpe.h>
  extern struct sigfpe_template sigfpe_[_N_EXCEPTION_TYPES+1];
#endif


/* ------------------------------------------------------------------------- *
 * File-scope type definitions.                                              *
 * ------------------------------------------------------------------------- */

typedef double (*PFD)(double);	/* Pointer to function returning double */

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
 * File-scope prototypes.                                                    *
 * ------------------------------------------------------------------------- */

static double  Log (double),           Log10   (double),
               Exp (double),           Sqrt    (double),
               Asin (double),          Acos    (double),
               Pow  (double, double),  integer (double);

static Symbol *lookup(char *);
static Symbol *install(char *, int, ...);
static void   *emalloc(size_t);



/* ------------------------------------------------------------------------- *
 * File-scope variables.                                                     *
 * ------------------------------------------------------------------------- */

static Symbol *symlist = NULL;	/* Internal use for function interpreter */
static Symbol *dlist   = NULL;	/* Storage of double lookup parameters   */
static Symbol *ilist   = NULL;	/* Storage of integer lookup parameters  */
static Symbol *olist   = NULL;	/* Storage for option lookup             */

static char   func_string[FILENAME_MAX], *cur_string;
static double value;

static struct {			/* Built-in functions */
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
  "asin"  ,  Asin    ,		/* Rest do error-checking */
  "acos"  ,  Acos    ,
  "log"   ,  Log     ,
  "log10" ,  Log10   ,
  "exp"   ,  Exp     ,
  "sqrt"  ,  Sqrt    ,
  NULL    ,  NULL
};

#include "defaults.h"

%}
/* ------------------------------------------------------------------------- *
 * Yacc grammar follows.                                                     *
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
 * Externally-visible functions follow.                                      *
 *****************************************************************************/


void initialize(void)
/* ========================================================================= *
 * Load lookup tables and symbol table with default values.                  *
 * Floating-point errors cause abortion; we have to take special action to   *
 * make this happen on a machine with current release of MIPS FPU.           *
 * This routine should be called at start of run-time.                       *
 * ========================================================================= */
{
  int i;
 

#ifdef _mips

  sigfpe_[_UNDERFL].repls = _ZERO;
  sigfpe_[ _OVERFL].abort = 1;
  sigfpe_[_DIVZERO].abort = 1;
  sigfpe_[_INVALID].abort = 1;

  handle_sigfpe(_ON, _EN_UNDERFL | _EN_OVERFL | _EN_DIVZERO | _EN_INVALID,
		NULL, _ABORT_ON_ERROR, NULL);
#endif

  for (i=0; consts[i].name; i++)
    install(consts[i].name, DPARAM, consts[i].cval);

  for (i=0; option[i].name; i++) 
    install(option[i].name, OPTION, option[i].oval);

  for (i=0; iparam[i].name; i++)
    install(iparam[i].name, IPARAM, iparam[i].ival);

  for (i=0; dparam[i].name; i++)
    install(dparam[i].name, DPARAM, dparam[i].dval);

  for (i=0; builtin[i].name; i++)
    install(builtin[i].name, BLTIN, builtin[i].func);
}





double interpret(char *s)
/* ========================================================================= *
 * Given a string, interpret it as a function using yacc-generated yyparse().*
 * ========================================================================= */
{
  if (strlen(s) > FILENAME_MAX)
    message("in interpret(): too many characters passed:\n", s, ERROR);
  
  strcpy(func_string, s);
  strcat(func_string, "\n");
  
  cur_string = func_string;
  
  yyparse();
  return value;
}





void set_option(char *s, int v)
/* ========================================================================= *
 * Set option on list true/false, or install it.                             *
 * ========================================================================= */
{
  Symbol *sp;
 

  for (sp = olist; sp; sp = sp->next)
    if (strcmp(sp->name, s) == 0) {
      sp->u.ival = v;
      break;
    }
  
  if (!sp) sp = install(s, OPTION, v);
}





int get_option(char *s)
/* ========================================================================= *
 * Retrieve value from option list.                                          *
 * ========================================================================= */
{
  Symbol *sp;
 
 
  for (sp = olist; sp; sp = sp->next)
    if (strcmp(sp->name, s) == 0)
      return sp->u.ival;

  message("get_option(): name not found: ", s, WARNING);
}
    




void show_option(FILE *fp)
/* ========================================================================= *
 * Echo option list to fp.                                                   *
 * ========================================================================= */
{
  Symbol *sp;

  
  for (sp = olist; sp; sp = sp->next)
    fprintf(fp, "%-12s%-3d\n", sp->name, sp->u.ival);
}





void set_iparam(char *s, int v)
/* ========================================================================= *
 * Set or install iparam on list.                                            *
 * ========================================================================= */
{
  Symbol *sp;

 
  for (sp = ilist; sp; sp = sp->next)
    if (strcmp(sp->name, s) == 0) {
      sp->u.ival = v;
      break;
    }
  
  if (!sp) sp = install(s, IPARAM, v);
}





int get_iparam(char *s)
/* ========================================================================= *
 * Retrieve value from iparam list.                                          *
 * ========================================================================= */
{
  Symbol *sp;
 
 
  for (sp = ilist; sp; sp = sp->next)
    if (strcmp(sp->name, s) == 0)
      return sp->u.ival;

  message("get_iparam(): name not found:", s, WARNING);
}
    




void show_iparam(FILE *fp)
/* ========================================================================= *
 * Echo iparam list to fp.                                                   *
 * ========================================================================= */
{
  Symbol *sp;

  
  for (sp = ilist; sp; sp = sp->next)
    fprintf(fp, "%-12s%-3d\n", sp->name, sp->u.ival);
}





void set_dparam(char *s, double v)
/* ========================================================================= *
 * Set or install dparam on list.                                            *
 * ========================================================================= */
{
  Symbol *sp;
 

  for (sp = dlist; sp; sp = sp->next)
    if (strcmp(sp->name, s) == 0) {
      sp->u.dval = v;
      break;
    }
  
  if (!sp) sp = install(s, DPARAM, v);
}





double get_dparam(char *s)
/* ========================================================================= *
 * Retrieve value from dparam list.                                          *
 * ========================================================================= */
{
  Symbol *sp;
  
  for (sp = dlist; sp; sp = sp->next)
    if (strcmp(sp->name, s) == 0)
      return sp->u.dval;

  message("get_dparam(): name not found:", s, WARNING);
}
    




void show_dparam(FILE *fp)
/* ========================================================================= *
 * Echo dparam list to fp.                                                   *
 * ========================================================================= */
{
  Symbol *sp;

  
  for (sp = dlist; sp; sp = sp->next)
    fprintf(fp, "%-12s%-.6g\n", sp->name, sp->u.dval);
}

/*****************************************************************************
 * Remaining routines are accessible only in this module.                    *
 *****************************************************************************/




static int yylex(void)
/* ========================================================================= *
 * Lexical analysis routine called by yyparse, using string loaded by        *
 * interpret().                                                              *
 * ========================================================================= */
{
  int c;


  while ((c = *cur_string++) == ' ' || c == 't');

  if (c == EOF)
    return 0;
  if (c == '.' || isdigit(c)) {	/* number */
    yylval.val = strtod(--cur_string, &cur_string);
    return NUMBER;
  }
  if (isalpha(c)) {		/* symbol */
    Symbol *s;
    char    sbuf[FILENAME_MAX], *p = sbuf;
    do {
      *p++ = c;
    } while ((c = *cur_string++) != EOF && (isalnum(c) || c == '_'));
    cur_string--;
    *p = '\0';
    if ((s=lookup(sbuf)) == NULL)
      s = install(sbuf, UNDEF, 0.0);
    yylval.sym = s;
    return (s->type == UNDEF || s->type == DPARAM) ? VAR : s->type;
  }
  return c;
}





static Symbol *lookup(char *s)
/* ========================================================================= *
 * Find s in symbol table.                                                   *
 * ========================================================================= */
{
  Symbol *sp;
  
  
  for (sp = symlist; sp; sp = sp->next)
    if (strcmp(sp->name, s) == 0)
      return sp;

  return NULL;
}





static Symbol *install(char *s, int t, ...)
/* ========================================================================= *
 * Install value of variable type into appropriate lists.                    *
 * Note that variables of type DPARAM get mounted in both symlist and dlist. *
 * ========================================================================= */
{
  Symbol *sp;
  va_list ap;
  
  
  va_start(ap, t);

  sp       = (Symbol *) emalloc(sizeof(Symbol));
  sp->name = emalloc(strlen(s)+1);
  strcpy(sp->name, s);

  switch (sp->type = t) {
  case OPTION:
    sp->u.ival = va_arg(ap, int);
    sp->next   = olist;
    olist      = sp;
    break;
  case IPARAM:
    sp->u.ival = va_arg(ap, int);
    sp->next   = ilist;
    ilist      = sp;
    break;
  case DPARAM:
    sp->u.dval = va_arg(ap, double);
    sp->next   = dlist;
    dlist      = sp;
    sp         = (Symbol *) emalloc(sizeof(Symbol));
    memcpy(sp, dlist, sizeof(Symbol));
    sp->name = emalloc(strlen(s)+1);
    strcpy(sp->name, s);
    symlist    = sp;
    break;
  case VAR:
  case UNDEF:
    sp->u.dval = va_arg(ap, double);
    sp->next   = symlist;
    symlist    = sp;
    break;
  case BLTIN:
    sp->u.ptr  = va_arg(ap, PFD);
    sp->next   = symlist;
    symlist    = sp;
    break;
  default:
    sp->u.dval = va_arg(ap, double);
    sp->type   = UNDEF;
    sp->next   = symlist;
    symlist    = sp;
    break;
  }

  va_end(ap);

  return (sp);
}





static void *emalloc(size_t n)
/* ========================================================================= *
 * Check return from malloc().                                               *
 * ========================================================================= */
{
  void *p;


  if (!(p = (void *) malloc(n)))
    message("emalloc()", "out of memory", ERROR);
 
  return p;
}





static void yyerror(char *s)
/* ========================================================================= *
 * Handler for yyparse() syntax errors.                                      *
 * ========================================================================= */
{
  message("yyparse()", s, WARNING);
}





static double errcheck(double d, char *s)
/* ========================================================================= *
 * Check result of math library call.                                        *
 * ========================================================================= */
{
  if (errno == EDOM) {
    errno = 0;
    message("errcheck(): argument out of domain in call to ", s, ERROR);
  } else if (errno == ERANGE) {
    errno = 0;
    message("errcheck(): result out of range in call to ", s, ERROR);
  }

  return d;
}


/*****************************************************************************
 * Remaining routines do error-checking calls to math library routines.      *
 *****************************************************************************/





static double Log(double x)
{
  return errcheck(log(x), " log");
}


static double Log10(double x)
{
  return errcheck(log10(x), " log10");
}


static double Exp(double x)
{
  return errcheck(exp(x), " exp");
}


static double Sqrt(double x)
{
  return errcheck(sqrt(x), " sqrt");
}


static double Pow(double x, double y)
{
  return errcheck(pow(x, y), " exponentiation");
}


static double Acos(double x)
{
  return errcheck(acos(x), " acos");
}



static double Asin(double x)
{
  return errcheck(asin(x), " asin");
}


static double integer(double x)
{
  return (double)(long)x;
}
