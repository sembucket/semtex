#ifndef MSCOPE_GRAMMAR_H
#define MSCOPE_GRAMMAR_H

/* Mscope's grammar definitions
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include "cubit/tree.h"

#ifdef DOMAIN   /* maybe in <math.h>? */
#undef DOMAIN
#endif

typedef struct {        /* ...... Token Definition ...... */
  int        type;      /* type, from the TokenType list  */
  char*      name;      /* name to be called              */
  char*      help;      /* help string (optional)         */
  PFV        action;    /* function to call for action    */
} Token;
  
typedef enum {          /* Commands understood by Mscope */

  /* I/O */

  LOAD,
  SAVE, 
  MESHIN,
  MESHOUT,
  FLOWIN,
  FLOWOUT,    
  INPUT,
                       
  /* computation */

  NORM,
  ENORM,
  SET,
  FSCAL,
  FABS,
  FAXPY,
  FCOPY,
  FADD,
  FMULT,
  FDIV,
  FGRAD,
  FFT,
  FINT,

  /* Drawing commands */

  ERASE,
  BOX,
  GRID, 
  GRIDNUM, 
  MESH, 
  LOGICAL,
  BCDRAW,
  INFO, 
  CONTOUR,
  VECTORS,
  PROFILE,
  SPLOT,
  LTYPE,
  RELOCATE,
  DRAW,
  LABEL,
  STREAM,
  
  ID, 
  FAMILY, 
  KEY, 
  LEVEL,
  MINFO,
  EINFO,
  BINFO,
  PARAM,
  DEFINE,
  OPTION,
  SHOW,
  LIST,

  DOMAIN,
  INSTALL, 
  RESTRICT, 
  PRUNE, 
  CONNECT, 
  LOCK, 
  KEYS,
  ELEMENT,
  NODE,
  CURVE,
  BCOND,
  ICOND,
  FIELDS,
  USER,

  /* Adaptive mesh generation and mesh refinement */

  REFINE, 
  SPLIT,
  GRADIENT,
  SPECTRUM,
  REGRESS,

  /* Helmholtz solver (manual and solution-adaptive) */

  MATRIX, 
  SOLVE, 
  FORCE, 
  HISTORY,
  SOLUTION,
  AUTO,

  /* Misc. */

  DO,
  DEV,
  LIMITS,  
  ZOOM,      
  UNZOOM,
  SLEEP,
  WARRANTY,
  COPYRIGHT,
  HELP,
  COMMAND,
  EXEC,
  QUIT

} TokenType;

/* Prototypes */

void DoErase    ();
void DoLoad     ();
void DoSave     ();
void DoMeshIn   ();
void DoMeshOut  ();
void DoInput    ();
void DoBox      ();
void DoGrid     ();
void DoContour  ();
void DoVectors  ();
void DoProfile  ();
void DoGridNum  ();
void DoMesh     ();
void DoBCs      ();
void DoRefine   ();
void DoSplit    ();
void DoGradient ();
void DoSpectrum ();
void DoRegress  ();
void DoInstall  ();
void DoKeys     ();
void DoRestrict ();
void DoID       ();
void DoFamily   ();
void DoKey      ();
void DoLevel    ();
void DoLogical  ();
void DoPrune    ();
void DoConnect  ();
void DoMinfo    ();
void DoEinfo    ();
void DoBinfo    ();
void DoParam    ();
void DoDefine   ();
void DoOption   ();
void DoLock     ();
void DoStream   ();

void DoFlowin   ();
void DoFlowout  ();
void DoNorms    ();
void DoEnorms   ();
void DoSet      ();

void DoFscal    ();
void DoFabs     ();
void DoFaxpy    ();
void DoFcopy    ();
void DoFadd     ();
void DoFmult    ();
void DoFdiv     ();
void DoFgrad    ();
void DoFint     ();
void DoFFT      ();

void DoSplot    ();

void DoElement  ();
void DoNode     ();
void DoCurve    ();
void DoBcond    ();
void DoIcond    ();
void DoFields   ();
void DoUser     ();
void DoShow     ();
void DoMatrix   ();
void DoSolve    ();
void DoForce    ();
void DoHistory  ();
void DoSolution ();
void DoWarranty ();
void DoCopyright();
void DoLoop     ();
void DoZoom     ();
void DoUnzoom   ();
void DoDomain   ();
void DoList     ();
void DoCommand  ();
void DoLtype    ();
void DoRelocate ();
void DoDraw     ();
void DoDevice   ();
void DoLabel    ();

#endif



