%{
/* --------------------------------------------------------------------- *
 * Manager:  Symbol table manager and parser                             *
 *                                                                       *
 * This file contains functions for managing a set of lookup tables and  *
 * several parsers.  The tables maintained are: (1) a global symbol tab- *
 * le for the parser containing mathematical constants and functions,    *
 * (2) a parameter table used by the spectral element solver, and (3) an *
 * option table used maintaining integer-valued options.  The interface  *
 * routines are as follows:                                              *
 *                                                                       *
 *                                                                       *
 * Internal Symbol Table                                                 *
 * ---------------------                                                 * 
 * Symbol *install(char *name, int type, ...)                            *
 * Symbol *lookup (char *name)                                           *
 *                                                                       *
 * Parameter Symbol Table                                                *
 * ----------------------                                                *
 * int     iparam     (char *name)                                       *
 * int     iparam_set (char *name, int value)                            *
 *                                                                       *
 * double  dparam     (char *name)                                       *
 * double  dparam_set (char *name, double value)                         *
 *                                                                       *
 * Options Table                                                         *
 * -------------                                                         *
 * int     option     (char *name)                                       *
 * int     option_set (char *name, int status)                           *
 *                                                                       *
 *                                                                       *
 * Vector/Scalar Parser                                                  *
 * --------------------                                                  *
 * The parsers provide two types of function-string parsing based on the *
 * type of access involved: a vector parser for forcing functions and    *
 * boundary conditions and a scalar parser for miscellaneous applica-    *   
 * tions.  The interfaces for these routines are:                        *
 *                                                                       *
 * void    vector_def (char *vlist, char *function)                      *
 * void    vector_set (int   vsize, v1, v2, ..., f(v))                   *
 *                                                                       *
 * double  scalar     (char *name)                                       *
 * double  parse      (char *expression)                                 *
 * --------------------------------------------------------------------- */
 
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <signal.h>
#include <setjmp.h>
#include <errno.h>

#include "tree.h"

#define SIZE   128       /* Maximum number of function string characters */

typedef double (*PFD)(); /* Pointer to a function returning double */

typedef struct Symbol {  /* Symbol table entry */
        char *name   ;   /* Symbol name        */
        short type   ;   /* VAR, BLTIN, UNDEF, DPARAM, IPARAM, OPTION */
	short status ;   /* (see status definitions) */
        union {
		int    num ;
                double val ;
		PFD    ptr ;
	      } u;
      } Symbol;

/* --------------------------------------------------------------------- *
 *                 Function declarations and prototypes                  *
 * --------------------------------------------------------------------- */

#include "veclib/veclib.h"
#include "manager.h"

/* static FILE *symbol_stream = stdout; Changed by hmb Jan 2002 */

static FILE *symbol_stream;

/* Internal Prototypes */

static Symbol *install(const char*, int, ...);       /* Table management */
static Symbol *lookup (const char*);

static double 
  Sqrt   (double),                           /* Operators (mathematical) */
  Rand   (double), 
  Integer(double), 
  Log    (double), 
  Log10  (double), 
  Exp    (double), 
  Coth   (double),
  Sech   (double),
  Csch   (double),
  Radius (double,double),                    /* ... binary operators ... */
  Angle  (double,double),
  Pow    (double,double),
  Sign   (double,double);


extern void show_symbol (Symbol *s),         /* Print symbol's value     */
            manager_init(void);

extern void warning  (const char *s, const char *t);
extern void yyerror  (char *s);
extern void execerror(char *s, char *t);
extern void execerrnr(char *s, char *t);
extern void fpecatch (void);

extern int yyparse  (void);
extern int yylex    (void);

/* --------------------------------------------------------------------- *
 *                                                                       *
 *          P R O G R A M    D E F A U L T   P A R A M E T E R S         *
 *                                                                       *
 * The following are the default values for options and parameters used  *
 * in the Helmholtz and Navier-Stokes solvers.                           *
 * --------------------------------------------------------------------- */

static struct {
  char   *name;
  int     oval;
} O_default[] = {              /* Options */
	 "binary",     1,
	 "direct",     1,
	 "core",       1,
	 "procid",     0,
	 "nprocs",     1,
	  0, 0
};

static struct {                /* Parameters (integer) */
  char   *name;
  int     pval;
} I_default[] = {
         "DIM",           2,
	 "NZ",            1,
	 "NORDER",        5,
	 "Rotational",    0,     /* ...should be overridden */
	 "SkewSymmetric", 1,
	 "Stokes",        2,
	 "Convective",    3,
	  0, 0
};

static struct {                  /* Parameters (double) */
  char   *name;
  double  pval;
} D_default[] = {
	 "TOL",        1.e-8,    /* General tolerance         */
	 "TOLFAM",     1.e-8,    /* Family matching tolerance */
	 "TOLCURV",    1.e-8,    /* Curve fitting tolerance   */
	 "TOLABS",     1.e-8,    /* Absolute difference tol.  */
	 "TOLREL",     1.e-8,    /* Relative difference tol.  */
	 "XSCALE",     1.,
	 "YSCALE",     1.,
	 
	 "Re",         1.,       /* ..... Navier-Stokes ..... */
	 "DENSITY",    1.,       /* Fluid density             */
	 "VELOCITY",   1.,       /* Base velocity scale       */
	 "LENGTH",     1.,       /* Base length scale         */
	 "TIME",       1.,       /* Time interval             */
	 "TIME_0",     0.,       /* Initial time              */
	 "TIME_N",     0.,       /* Final time                */
         "DT",         0.001,    /* Time step                 */
	 "FFZ",        0.,       /* Forcing                   */
	 "FFY",        0.,       /*                           */
	 "FFX",        0.,       /*                           */
	 "FLOWRATE",   0.,       /* Flowrate                  */
	 "BETA",       1.,       /* Spanwise wavenumber       */
	  0,           0.
};


static struct {                  /* Constants */
  char    *name;
  double   cval;
} consts[] = {
        "PI",     3.14159265358979323846,   /* Pi */
        "TWOPI",  6.28318530717958647688,   /* 2*Pi */
	"E",      2.71828182845904523536,   /* Natural logarithm */
	"GAMMA",  0.57721566490153286060,   /* Euler constant */
	"DEG",   57.29577951308232087680,   /* deg/radian */
	"PHI",    1.61803398874989484820,   /* golden ratio */
	 0,       0
};

static struct {                /* Built-ins */
  char    *name;               /* Function name */
  short    args;               /* # of arguments */
  PFD      func;               /* Pointer to the function */
} builtins[] = {
  "acos",  1,  acos,
  "asin",  1,  asin,
  "atan",  1,  atan,
  "cos",   1,  cos,
  "sin",   1,  sin,
  "tan",   1,  tan,
  "cosh",  1,  cosh,
  "sinh",  1,  sinh,
  "tanh",  1,  tanh,
  "csch",  1,  Csch,
  "sech",  1,  Sech,
  "coth",  1,  Coth,
  "exp",   1,  Exp,         
  "log",   1,  Log,         
  "log10", 1,  Log10,       
  "sqrt",  1,  Sqrt,        
  "abs",   1,  fabs,
  "int",   1,  Integer,     
  "rand",  1,  Rand,        /* random number (input the magnitude) */
  "rad",   2,  Radius,      /* rad = sqrt(x^2 + y^2) */
  "ang",   2,  Angle,       /* ang = atan2(x,y)      */
  "sign",  2,  Sign,        /* transfer of sign      */
  "atan2", 2,  atan2,
  0,       0,  0
};

/* External variables */

Tree*    Symbols  = 0;     /* Symbol table     */
Tree*    Options  = 0;     /* Option table     */
Tree*    Params   = 0;     /* Parameters table */
jmp_buf  begin;            

static char     func_string[SIZE], 
                *cur_string;
static double   stack_value;

extern int errno;

%}
%union {                /* stack type */
	double  val;    /* actual value */
	Symbol *sym;    /* symbol table pointer */
}
%token	<val>	NUMBER
%token	<sym>	VAR BLTIN_UNARY BLTIN_BINARY UNDEF DPARAM IPARAM OPTION
%type	<val>	expr asgn
%right	'='
%left	'+' '-'        	/* left associative, same precedence */
%left	'*' '/'        	/* left associative, higher precedence */
%left	UNARYMINUS
%right	'^'		/* exponentiation */
%%
list:  /* nothing */
	| list '\n'
	| list asgn  '\n'
	| list expr  '\n'     	{ stack_value = $2; }
	| list error '\n'	{ yyerrok; }	
	;
asgn:	  VAR '=' expr { $$=$1->u.val=$3; $1->type = VAR; }
	;
expr:	  NUMBER  { $$ = $1; }
	| VAR     { if ($1->type == UNDEF)
		      execerrnr("undefined variable -- ",$1->name);
		    $$ = $1->u.val; }
        | IPARAM  { $$ = (double) $1->u.num; }
        | DPARAM  { $$ = $1->u.val; }
	| asgn
	| BLTIN_UNARY  '(' expr ')'	
            { $$ = (*($1->u.ptr))($3); }
	| BLTIN_BINARY '(' expr ',' expr ')'	
            { $$ = (*($1->u.ptr))($3,$5); }
	| expr '+' expr { $$ = $1 + $3; }
	| expr '-' expr { $$ = $1 - $3; }
	| expr '*' expr { $$ = $1 * $3; }
	| expr '/' expr {
		if ($3 == 0.0)
			execerror("division by zero","");
		 $$ = $1 / $3; }
	| expr '^' expr	{ $$ = Pow($1,$3); }
	| '(' expr ')'	{ $$ = $2; }
	| '-' expr %prec UNARYMINUS { $$ = -$2; }
	;
%%
	/* end of grammer */

/* --------------------------------------------------------------------- *
 *                                                                       *
 *                              P A R S E R                              *
 *                                                                       *
 * --------------------------------------------------------------------- */

int yylex (void)
{
	int c;

	while((c = *cur_string++) == ' ' || c == '\t');

	if(c == EOF)
		return 0;
	if(c == '.' || isdigit(c)) {                      /* number */
	        char *p;
	        yylval.val = strtod(--cur_string, &p);
		cur_string = p;
		return NUMBER;
	}
	if(isalpha(c)) {                                  /* symbol */
		Symbol *s;
		char sbuf[100], *p = sbuf;

		do
		  *p++ = c;
		while
		  ((c = *cur_string++) != EOF && (isalnum(c) || c == '_')); 

		cur_string--;
		*p = '\0';
		if(!(s=lookup(sbuf))) 
		  s = install(sbuf, UNDEF, 0.);
		yylval.sym = s;
		return (s->type == UNDEF) ? VAR : s->type;
	}

	return c;
}

void yyerror(char *s)      /* called for yacc syntax error */
{
  warning (s, (char *) 0);
}

void warning(const char *s, const char *t)    /* print warning message */
{
  fprintf(stderr,"parser: %s",s);
  if (t)
    fprintf(stderr," %s\n",t);
  else
    fprintf(stderr," in function string %s\n",func_string);
}

void execerror(char *s, char *t)    /* recover from run-time error */
{
  warning (s,t);
  longjmp (begin,0);
}

void execerrnr(char *s, char *t)   /* run-time error, no recovery */
{
  warning(s,t);
  fprintf(stderr,"exiting to system...\n");
  exit(-1);
}

void fpecatch (void)           /* catch floating point exceptions */
{
  fputs ("cubit: floating point exception\n"
	 "exiting to system...\n", stderr);
  exit  (-1);
}

/* --------------------------------------------------------------------- *
 * Vector/Scalar parser                                                  *
 *                                                                       *
 * The scalar and vector parsers are the interfaces to the arithmetic    *
 * routines.  The scalar parser evaluates a single expression using      *
 * variables that have been defined as PARAM's or VAR's.                 *
 *                                                                       *
 * The vector parser is just a faster way to call the scalar parser.     *
 * Using the vector parser involves two steps: a call to vector_def() to *
 * declare the names of the vectors and the function, and a call to      *
 * vector_set() to evaluate it.                                          *
 *                                                                       *
 * Example:    vector_def ("x y z", "sin(x)*cos(y)*exp(z)");             *
 *             vector_set (100, x, y, z, u);                             *
 *                                                                       *
 * In this example, "x y z" is the space-separated list of vector names  *
 * referenced in the function string "sin(x)...".  The number 100 is the *
 * length of the vectors to be processed.  The function is evaluted as:  *
 *                                                                       *
 *             u[i] = sin(x[i])*cos(y[i])*exp(z[i])                      *
 *                                                                       *
 * --------------------------------------------------------------------- */

#define  VMAX     10    /* maximum number of vectors in a single call */
#define  VLEN   SIZE    /* maximum vector name string length          */

double parse (const char *expr)
{
  Symbol *sp;

  if (strlen(expr) > SIZE-1)
    execerrnr ("Too many characters in expression:\n", (char*) expr);
  
  sprintf (cur_string = func_string, "%s\n", expr);
  yyparse ();
  
  return stack_value;
}

double scalar (const char *name)
{
  Symbol *sp;

  if (strlen(name) > SIZE-1)
    execerrnr ("Too many characters in expression:\n", (char*) name);
  
  sprintf (cur_string = func_string, "%s\n", name);
  yyparse ();
  
  return stack_value;
}

double scalar_set (const char *name, double val)
{
  Node   *np;
  Symbol *sp;

  if (np = tree_search (Symbols->root, name)) {
    if ((sp = (Symbol*) np->other)->type == VAR)
      sp->u.val = val;
    else
      warning (name, "has a type other than VAR.  Not set.");
  } else
    install (name, VAR, val);

  return val;
}

static int     nvec;
static Symbol *vs[VMAX];
static double *vv[VMAX];

void vector_def (char *vlist, char *function)
{
  Symbol  *s;
  char    *name, buf[VLEN];

  if (strlen(vlist) > SIZE)
    execerrnr("name string is too long:\n", vlist);
  else
    strcpy(buf, vlist);

  /* install the vector names in the symbol table */

  name = strtok(buf, " ");
  nvec = 0;
  while (name && nvec < VMAX) {
    if (!(s=lookup(name))) 
      s = install (name, VAR, 0.);
    vs[nvec++] = s;
    name  = strtok((char*) NULL, " ");
  }

  if (strlen(function) > SIZE-1)
    execerrnr("too many characters in function:\n", function);

  sprintf (func_string, "%s\n", function);

  return;
}

void vector_set (int n, ...)
{
  va_list  ap;
  double   *fv;
  int i;

  /* initialize the vectors */

  va_start(ap, n);
  for (i = 0; i < nvec; i++) vv[i] = va_arg(ap, double*);
  fv = va_arg(ap, double*);
  va_end(ap);

  /* evaluate the function */

  while (n--) {
    for (i = 0; i < nvec; i++) vs[i]->u.val = *(vv[i]++);    
    cur_string = func_string; 
    yyparse();
    *(fv++)    = stack_value;
  }

  return;
}

#undef VMAX
#undef VLEN

/* --------------------------------------------------------------------- *
 * Parameters and Options                                                *
 *                                                                       *
 * The following functions simply set and lookup values from the tables  *
 * of variables.   If a symbol isn't found, they silently return zero.   *
 * --------------------------------------------------------------------- */

int iparam (const char *name)
{
  Node   *np;
  Symbol *sp;
  int    num = 0;

  if ((np = tree_search (Params->root, name)) &&
      (sp = (Symbol*) np->other)->type == IPARAM)
    num = sp->u.num;

  return num;
}

int iparam_set (const char *name, int num)
{
  Node   *np;
  Symbol *sp;

  if (np = tree_search (Params->root, name)) {
    if ((sp = (Symbol*) np->other)->type == IPARAM)
      sp->u.num = num;
    else
      warning (name, "has a type other than IPARAM.  Not set.");
  } else
    install (name, IPARAM, num);

  return num;
}

double dparam (const char *name)
{
  Node   *np;
  Symbol *sp;
  double val = 0.;

  if ((np = tree_search (Params->root, name)) &&
      (sp = (Symbol*) np->other)->type == DPARAM)
    val = sp->u.val;

  return val;
}

double dparam_set (const char *name, double val)
{
  Node   *np;
  Symbol *sp;

  if (np = tree_search (Params->root, name)) {
    if ((sp = (Symbol*) np->other)->type == DPARAM)
      sp->u.val = val;
    else
      warning (name, "has a type other than DPARAM.  Not set.");
  } else
    install (name, DPARAM, val);

  return val;
}

int option (const char *name)
{
  Node *np;
  int   status = 0;
  
  if (np = tree_search (Options->root, name))
    status = ((Symbol*) np->other)->u.num;

  return status;
}

int option_set (const char *name, int status)
{
  Node *np;

  if (np = tree_search (Options->root, name))
    ((Symbol*) np->other)->u.num = status;
  else
    install (name, OPTION, status);
  
  return status;
}

/* --------------------------------------------------------------------- *
 * manager_init()                                                        *
 *                                                                       *
 * The following function must be called before any other parser func-   *
 * tions to install the symbol tables and builtin functions.             *
 * --------------------------------------------------------------------- */

void manager_init (void)
{
  int i;

  /* initialize the trees */

  Symbols = create_tree (show_symbol, free);
  Options = create_tree (show_symbol, free);
  Params  = create_tree (show_symbol, free);

  symbol_stream = stdout;	/* Changed by hmb 2002 */

  /* initialize the signal manager */

  setjmp(begin);
  signal(SIGFPE, (void(*)()) fpecatch);

  /* options and parameters */

  for(i = 0; O_default[i].name; i++)
     install(O_default[i].name,OPTION,O_default[i].oval);
  for(i = 0; I_default[i].name; i++) 
     install(I_default[i].name,IPARAM,I_default[i].pval);
  for(i = 0; D_default[i].name; i++)
     install(D_default[i].name,DPARAM,D_default[i].pval);

  /* constants and built-ins */

  for(i = 0; consts[i].name; i++)
    install (consts[i].name,VAR,consts[i].cval);
  for(i = 0; builtins[i].name; i++) {
    switch  (builtins[i].args) {
    case 1:
      install (builtins[i].name, BLTIN_UNARY, builtins[i].func);
      break;
    case 2:
      install (builtins[i].name, BLTIN_BINARY, builtins[i].func);
      break;
    default:
      execerrnr ("too many arguments for builtin:", builtins[i].name);
      break;
    }
  }
  
  return;
}

/* Is this symbol defined? */

int defined (char *symbol)
{
  if (lookup(symbol) != (Symbol*) NULL)
    return 1;
  else
    return 0;
}


/* Print parameter, option, and symbol tables */



void show_symbols(FILE *fp) 
   { fputs ("\nSymbol table:\n", symbol_stream = fp); tree_walk (Symbols); }
void show_options(FILE *fp) 
   { fputs ("\nOptions:\n",      symbol_stream = fp); tree_walk (Options); }
void show_params (FILE *fp) 
   { fputs ("\nParameters:\n",   symbol_stream = fp); tree_walk (Params);  }

/* Print a Symbol */

void show_symbol (Symbol *s)
{
  switch (s->type) {
  case OPTION:
  case IPARAM:
    fprintf (symbol_stream, "%-15s -- %d\n", s->name, s->u.num);
    break;
  case DPARAM:
  case VAR:
    fprintf (symbol_stream, "%-15s -- %g\n", s->name, s->u.val);
    break;
  default:
    break;
  }

  return;
}

/* ..........  Symbol Table Functions  .......... */

static Symbol *lookup (const char *key)
{
  Node *np;

  if (np = tree_search (Symbols->root, key))
    return (Symbol*) np->other;

  if (np = tree_search (Params ->root, key))
    return (Symbol*) np->other;

  if (np = tree_search (Options->root, key))
    return (Symbol*) np->other;

  return (Symbol*) NULL;     /* not found */
}                  

/* 
 * install "key" in a symbol table 
 */

static Symbol *install (const char *key, int type, ...)     
{
  Node   *np;
  Symbol *sp;
  Tree   *tp;
  va_list ap;

  va_start (ap, type);
  
  /* Get a node for this key and create a new symbol */

  np       = create_node (key);
  sp       = (Symbol *) malloc(sizeof(Symbol));
  sp->name = np->name;

  switch (sp->type = type) {
  case OPTION:
    tp        = Options;
    sp->u.num = va_arg(ap, int);
    break;
  case IPARAM:
    tp        = Params;
    sp->u.num = va_arg(ap, int);
    break;
  case DPARAM:
    tp        = Params;
    sp->u.val = va_arg(ap, double);
    break;
  case VAR: 
  case UNDEF:  
    tp        = Symbols;
    sp->u.val = va_arg(ap, double);
    break;
  case BLTIN_UNARY:
  case BLTIN_BINARY:
    tp        = Symbols;
    sp->u.ptr = va_arg(ap, PFD);
    break;
  default:
    tp        = Symbols;
    sp->u.val = va_arg(ap, double);
    sp->type  = UNDEF;
    break;
  }

  va_end (ap);

  np->other = (void *) sp;     /* Save the symbol */
  tree_insert (tp, np);        /* Insert the node */

  return sp;
}

/*
 *  Math Functions
 *  --------------  */

static double errcheck (double d, char *s)
{
  if (errno == EDOM) {
    errno = 0                              ;
    execerror(s, "argument out of domain") ;
  }
  else if (errno == ERANGE) {
    errno = 0                           ;
    execerror(s, "result out of range") ;
  }
  return d;
}

static double Log (double x)
{
  return errcheck(log(x), "log") ;
}

static double Log10 (double x) 
{
  return errcheck(log10(x), "log10") ;
}

static double Exp (double x)
{
  return errcheck(exp(x), "exp") ;
}

static double Sqrt (double x)
{
  return errcheck(sqrt(x), "sqrt") ;
}

static double Pow (double x, double y)
{
  const
  double yn = floor(y + .5);
  double px = 1.;

  if (yn >= 0 && yn == y) {     /* Do it inline if y is an integer power */
      int n = yn;
      while (n--) 
         px *= x;
  } else  
      px = errcheck (pow(x,y), "exponentiation");

  return px;
}

static double Integer (double x)
{
  return (double) (long) x;
} 

static double Rand (double x)
{
  return x * drand();
}

static double Sign (double x, double y)
{
  if (y > 0)
    return  fabs(x);
  else if (y < 0)
    return -fabs(x);
  else
    return  0.;
}

static double Radius (double x, double y)
{
  if (x != 0. || y != 0.)
    return sqrt (x*x + y*y);
  else
    return 0.;
}

#ifndef M_PI
#define M_PI  consts[0].cval
#endif

static double Angle (double x, double y)
{
  double theta = 0.;

  if (x != 0.)
    theta =  atan2 (y,x);
  else if (y > 0.)
    theta =  M_PI / 2.;
  else if (y < 0.)
    theta = -M_PI / 2.;
  
  return theta;
}

static double Coth (double z) {
  return cosh(z)/sinh(z);
}

static double Sech (double z) {
  return 1./cosh(z);
}

static double Csch (double z) {
  return 1/sinh(z);
}

#undef M_PI
