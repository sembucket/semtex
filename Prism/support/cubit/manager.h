#ifndef MANAGER_H
#define MANAGER_H

/* Global Symbol Table
 *
 * The library maintains a global symbol table for a number of things:
 * parsing boundary conditions, defining parameter values, etc.  Since the
 * table is global (all symbols are stored internally by the library), there
 * is no public data structure.  Instead there is simply a set of functions
 * that let you look up symbols and parse symbolic expressions. 
 *
 * This should be self-explanatory.
 * ------------------------------------------------------------------------ */
 
double parse      (const char *expr);

int option        (const char *name);
int option_set    (const char *name, int value);

int iparam        (const char *name);
int iparam_set    (const char *name, int value);

double dparam     (const char *name);
double dparam_set (const char *name, double value);

double scalar     (const char *expr);
double scalar_set (const char *name, double value);

void vector_def   (char *vlist, char *expr);
void vector_set   (int   vsize, ... /* v1, v2, ..., f(v) */ );

void show_symbols (FILE *fp);
void show_params  (FILE *fp);
void show_options (FILE *fp);

#endif
