#ifndef PRISM_H
#define PRISM_H

/*
 * Prototypes and macros for PRISM
 *
 * $Id$
 * 
 * Copyright (c) 1994-1998 R. D. Henderson
 *
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdarg.h>

#ifndef DIM
#error  Please define DIM to be 2 or 3!
#endif

#include "speclib/speclib.h"
#include "comm/comm.h"

#include "prism/config.h"
#include "prism/constants.h"
#include "prism/history.h"
#include "prism/measure.h"
#include "prism/greens.h"
#include "prism/domain.h"
#include "prism/hooks.h"
#include "prism/stats.h"

#ifdef MORTAR
#include "mortar/mortar.h"
#endif

#ifdef PARALLEL
#  define ROOTONLY if(!comm_rank())
#  define GSYNC    comm_sync()
#else
#  define ROOTONLY /* nothing */
#  define GSYNC    /* nothing */
#endif

#define MEMCHK { char *t=malloc(1); free (t); puts("Memory is OK"); }

#define Velocity_sys Velocity     /* Backwards compatibility */ */
#define Pressure_sys Pressure

typedef Field* Scalar;            /* ......... Scalar Field .......... */

typedef struct Vector {           /* ......... Vector in R3 .......... */
  Scalar  x, y, z;
} Vector;

/* ------------------------------- Prototypes ------------------------------ */

/* These are the functions that initialize and shutdown the application */

void Prism_init (int *argc, char **argv[]);
void Prism_exit (void);

/* Log files and error handling.  All diagnostic output should be directed   *
 * through either Prism_log() or Prism_error().  These provides a single set *
 * of functions that make it easy to redirect output to a new log file, for  *
 * example.  The usage of both Prism_log() and Prism_error() mimics the      *
 * standard I/O routines (printf).                                           *
 *                                                                           *
 * The function Prism_log() also acts as a filter based on the value of the  *
 * parameter <level> and the option "verbose".  Any message whose level is   *
 * greater than the value of "verbose" will not appear on the log output.    *
 * For example, if verbose = 0 then only level 0 messages are printed, while *
 * level 1 and above are ignored.  The usual values of <level> are:          *
 *                                                                           *
 *          0   warnings                                                     *
 *          1   standard diagnostics [default]                               *
 *          2   additional information that might be useful                  *
 *          9   debugging information for development                        */

void Prism_error  (const char *fmt, ...);
void Prism_log    (int level, const char *fmt, ...);
void Prism_logfile(const char *log);
void Prism_logfp  (FILE *fp);

/* Other stuff */

void parse_args   (int argc, char *argv[]);
void Summary      (void);
void ReadICs      (FILE*, int, ...);
void ReadDF       (FILE*, int, ...);
void ReadHisData  (FILE*, Domain*);
void PostProcess  (Domain *);

/* DRIVE */

void March_dt     (Domain *omega, double t, double dt);
void March_step   (Domain *omega, double dt, int nsteps);
void March_cfl    (Domain *omega, double t, double cfl);

int  Navier_Stokes (Domain *omega, double t0, double t1);
void MakeF         (Domain *omega, ACTION What, int step, double time);
void MultiSolve    (Domain *omega, Field   **U, Field   ***Uf, Bedge **Ubc,
		    BSystem *B, int nfields, int step);
void Startup       (Domain *omega, double *time0, int *step0);

/* ANALYZE */

void   Analyzer   (Domain *omega, double time, int step);
int    History    (Domain *omega, double time);
double Divergence (Domain *omega);
double Momentum   (Domain *omega);
void   Vorticity  (Domain *omega);

/* FLOWRATE */

void   flowrate_init (Domain *omega);
void   flowrate      (Domain *omega);


/* PRESSURE */

Bedge* BuildPBCs   (Field *U, Bedge *Ubc);
void   ComputePBCs (Domain  *omega);
void   SetPBCs     (Domain  *omega);

/* TIMESTEP */

void   Integrate    (Field *U, Field *Us[], Field *Fs[], int Je, double dt);
void   set_Itype    (int type);
void   set_order    (int Je);
void   get_alpha    (double *alpha);
void   get_beta     (double *beta);
double get_gamma    (int Je);
void   estimate_CFL (Domain *omega);

/* Advection Operators */

void SkewSymm (Domain *omega);
void VxOmega  (Domain *omega);
void StokesBC (Domain *omega);
 
/* FOURIER [3D Only] */

void Transform       (Field *U, double *Fu, ACTION direction);
void Transform32     (Field *U, double *Fu, ACTION direction);
void Field_gradz     (Field *U, Field *Uz);
void Field_davg_3D   (Field *U, BSystem *B);
void Field_grad_3D   (Field *U, Field *Ux, Field *Uy, double *wx, double *wy);
void Element_grad_3D (Element *U, Element *Ux, Element *Uy, int ntot);
void Element_gradz   (Element *U, Element *Ux, int ntot);

#endif
