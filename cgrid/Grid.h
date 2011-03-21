/*****************************************************************************
 * Header file for C & O-grid generators.
 * 
 * $Id$
 *****************************************************************************/


typedef struct {
  double  t;
  double  a[5];
} NACAparam;


void   odeint   (double* , int, double, double, double, double, double, int*,
	         int*, int, int*, double**, double,
	         void(*)(double, double*, double*, void*), void*,
	         void(*)(double*, double*, int, double*, double,
		         double, double*, double*, double*,
		         void(*) (double, double*, double*, void*), void*));

void   rkqs     (double*, double*, int, double*, double, double, double*,
	         double*, double*,
	         void(*)(double, double*, double*, void*), void*);

void   NACAdsdx (double, double*, double*, void*);

double NACAfoil (double, void*);

double stretch1 (double, double, double);

