/*                                                   
 * Spectral Element Data Model
 *
 * Based (loosely) on NEKTON's view of finite element discretizations
 *                                                                          
 * This file contains functions to translate a spectral element input file 
 * in the NEKTON format into the set of data structures used by SpecLib.
 * The file format is fairly specific, expecting the following sections in 
 * order:                                             
 *                                                                          
 *     PARAMETERS          . Symbol table parameters                        
 *     Passive Scalars     . Not used                                       
 *     LOGICAL             . Not used                                       
 *     MESH                . Element geometry                               
 *     CURVED SIDES        . Element boundary curves, part of the mesh
 *     BOUNDARY CONDITIONS . Boundary conditions for all fields             
 *                                                                          
 * Section names in the file are always compared by ignoring case, so you can
 * use upper, lower, or mixed-case section names in the file.   The functions 
 * in this file process each section to define the following parameters:
 *
 *
 *     ReadParams   ---  Load symbols into the parser's symbol table
 *     ReadPscals   ---          
 *     ReadLogics   ---             
 *     ReadMesh     ---  Create a Field (array of elements) with the given
 *                       geometry
 *     ReadBCs      ---  Create an array of Bedges (boundary edges) that
 *                       define the boundary conditions for a Field
 *     ReadVPs      ---
 *
 *
 * Notes:                                                                   
 *                                                                          
 *   - Parameters may be either integers or reals.  Fortran naming conven-  
 *     tions are used to decide for UPPERCASE parameters, but any lower-   
 *     case parameter is assumed to be real.  There are execeptions (see    
 *     below) for certain standard variable names.                          
 *                                                                          
 *   - Boundary conditions are a bit tricky, so please read the comments    
 *     in ReadBCs if you aren't sure how to call this function.           
 *                                                                          
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "speclib/isomesh.h"

/* Private functions */

static int checklist (char *name, char *list[]);
static int bedgcmp   (const void *, const void *);
                                              /* .....  Section names ..... */
#define  SEC_PARAMS   sections[0]             /* Parameters                 */
#define  SEC_PSCALS   sections[1]             /* Passive scalars            */
#define  SEC_LOGICS   sections[2]             /* Logical switches           */
#define  SEC_MESH     sections[3]             /* Mesh definition            */
#define  SEC_CURVE    sections[4]             /* Curved sides               */
#define  SEC_BCS      sections[5]             /* Boundary conditions        */
#define  SEC_ICS      sections[6]             /* Initial conditions         */
#define  SEC_DF       sections[7]             /* Drive force                */
#define  SEC_HIST     sections[8]             /* History points             */

/* The NEKTON crowd chose these... */

static char *sections[] = { 
  "PARAMETERS", "PASSIVE SCALARS", "LOGICAL", "MESH", "CURVED SIDES",
  "BOUNDARY", "INITIAL CONDITIONS", "DRIVE FORCE", "HISTORY" };

void ReadParams (FILE *fp)
{
  int    n;
  char   buf[BUFSIZ], value[25], name[25], c;

  static char *dspecial[] = { "LAMBDA", "LZ", 0 };
  static char *ispecial[] = { "EQTYPE", "TORDER", 0 };

  rewind (fp);
  if (!findSection (SEC_PARAMS, buf, fp))
    speclib_error("session: can't find PARAMETERS in the input file\n");
  
  /*
  // SEC_PARAMS begins with the program version, a dimension statement, 
  // and then the number of parameters.  Skip the first two and leave
  // the last one stored in buf.
  */

  for (n = 0; n < 3; n++) fgets(buf, BUFSIZ, fp);

  if (sscanf(buf, "%d", &n) != 1)
    speclib_error("session [%s]: can't read the # of parameters", SEC_PARAMS);
  
  while (n--) {
    fgets (buf, BUFSIZ, fp);
    if(sscanf(buf, "%25s%25s", value, name) == 2) {
      if (checklist (name, dspecial))
	dparam_set (name, scalar(value));
      else if (checklist (name, ispecial))
	iparam_set (name, scalar(value));
      else if (isupper(c = *name) && 'I' <= c && c <= 'N')
	iparam_set (name, scalar(value));
      else
	dparam_set (name, scalar(value));
    }
  }

  /* Running at a different NORDER ? */

  if ((n = option("norder")))
    iparam_set("NORDER", n);
  
  iparam_set ("NR", iparam("NORDER"));
  iparam_set ("NS", iparam("NORDER"));

  return;
}

static int checklist (char *name, char *list[])
{
  do
    if (strcmp(name,*list) == 0)
      return 1;
  while 
    (*++list);

  return 0;
}

/* Read the mesh description */

Field *ReadMesh (FILE *fp)
{
  const int nr = iparam("NR");
  const int ns = iparam("NS");
  const int nz = iparam("NZ");
  const int np = option("nprocs");

  Field *u;

  /* Automatic Parallel Domain Mapping 
   * ---------------------------------                                    *
   * This is an important point.  The input file parameter NZ is used to  *
   * specify how many GLOBAL frames will be created.  The options pid and *
   * nprocs specify the ID and number of processors this domain is dis-   *
   * tributed over.  Most SpecLib functions don't care about the GLOBAL   *
   * size of the domain (they work on a single frame), but in a few cases *
   * they do.  I've tried to mark them with the above string, so you can  *
   * grep for that in the source files to see where to be careful.        *
   *                                                                      *
   * Summary:      NZ       iparam      Global domain size                *
   *               nprocs   option      Number of processors              *
   *               pid      option      Local processor ID number         *
   * -------------------------------------------------------------------- */

  if (np > 1)
    u = Field_build(fp, nr, ns, nz/np, 'u');
  else
    u = Field_build(fp, nr, ns, nz, 'u');

  iparam_set ("ELEMENTS", Field_count(u));
  
  return u;
}

/* ------------------------------------------------------------------------ *
 * ReadBCs() - Read Boundary Conditions                                     *
 *                                                                          *
 * This function reads boundary information for a spectral element field    *
 * and sets up the appropriate boundary structures.  The following boundary *
 * conditions are currently recognized:                                     *
 *                                                                          *
 *   Flag      Type         Field #0       Field #1       Field #2          *
 *  -----    --------      ----------     ----------     ----------         *
 *    V      Velocity      U-velocity     V-Velocity     W-Velocity         *
 *    W      Wall          = 0            = 0            = 0                *
 *    v      v             U(x,y,z)       V(x,y,z)       W(x,y,z)           *
 *    F      Flux          U' = f1        V' = f2        W' = f3            *
 *    f      flux          U'(x,y,z)      V'(x,y,z)      W'(x,y,z)          *
 *    O      Outflow       U'.n = 0       V'.n = 0       W'.n = 0           *
 *    E      Element       w/ ID          w/ EDGE                           *
 *    P      Periodic      w/ ID          w/ EDGE                           *
 *    D      Dirichlet     patch number   segment number                    *
 *    N      Neumann       patch number   segment number                    *
 *    M      Master        patch number   segment number                    *
 *    S      Slave         patch number   segment number                    *
 *                                                                          *
 * NOTES:                                                                   *
 *                                                                          *
 *   - The argument "fnum" specifes which bc slot to read from.             *
 *     If fnum = 0, then ReadBCs will scan through the file until it finds  *
 *     another section of boundary conditions.  Otherwise, it reads back    *
 *     through the current section and selects boundary conditions from     *
 *     the specified field.                                                 *
 *                                                                          *
 *   - For boundaries with a given function of (x,y,z),  the boundary       *
 *     conditions are read from the lines following rather than the         *
 *     columns.  If ReadBCs finds another boundary condition before it      *
 *     reads "fnum" lines, it exits with an error.                          *
 *                                                                          *
 *     Here is an example (element 1, edge 1):                              *
 *                                                                          *
 *            1  1  v   0.0     0.0     0.0                                 *
 *                 ux = (1 - y^2)*sin(t)                                    *
 *                 uy = (1 - y^2)*cos(t)                                    *
 *                 uz = 0.                                                  *
 *            1  2  E   2.0     1.0     0.0                                 *
 *                                                                          *
 *     The "=" MUST be present, and the function is defined by the string   *
 *     to the right of it.                                                  *
 *                                                                          *
 *   - NO TIME DEPENDENT BOUNDARY CONDITIONS.                               *
 *                                                                          *
 *   - This routine will read any type of boundary condition, FLUID or      *
 *     otherwise, as long as the file is positioned correctly.  For examp-  *
 *     le, to read FLUID and THERMAL boundary conditions, you would call:   *
 *                                                                          *
 *            Ubc = ReadBCs (fp, 1, U);                                     *
 *            Vbc = ReadBCs (fp, 2, V);                                     *
 *            Wbc = ReadBCs (fp, 3, W);   //   Fluid BC's //                *
 *                                                                          *
 *            Tbc = ReadBCs (fp, 1, T);   // Thermal BC's //                *
 *                                                                          *
 *     After the last call, the file position is at the end of the second   *
 *     set of boundary conditions.  Each time ReadBCs() is called with      *
 *     fnum = 1, it scans through the file for another section of BCs.      *
 * ------------------------------------------------------------------------ */

#define  on  1
#define  off 0

static fpos_t bc_start;
static void   prescan (FILE *fp, fpos_t *p, char *buf);

Bedge *ReadBCs (FILE *fp, int fnum, Field *U)
{
  char     bc, buf[BUFSIZ];
  int      iel, iside, i,
           nbcs       = 0,
           dim        = iparam ("DIM");
  Bedge   *bndry_list = (Bedge *) NULL,
          *new_bndry  = (Bedge *) NULL;
  Bedge*   bsort (Bedge *BC, int nBCs);   /* for sorting the BC's      */

  if (fnum == 0)
    prescan (fp, &bc_start, buf);  /* Look for a new section of BC's   */
  else
    fsetpos (fp, &bc_start);       /* Re-read the current section      */


  while (U) {     /* ............ Begin Reading Boundary Conditions ...... */

    Edge  *edg = U->edges;
    double f[_MAX_FIELDS];

    for (edg = U->edges; edg; edg = edg->next) {
      char *p = fgets (buf, BUFSIZ, fp);
      while (isspace(*p)) p++; bc = *p++;
      sscanf(p, "%d%d%lf%lf%lf", &iel, &iside, &f[0], &f[1], &f[2]);
 
      if (iel - U->id != 1 || iside - edg->id != 1)
	speclib_error("Mismatched element/side -- got %d,%d, expected %d,%d",
		      iel, iside, U->id+1, edg->id+1);

      switch (bc) {

	/* Patching boundaries */

      case 'N': case 'M': case 'S':
      case 'D':
	new_bndry = BC_make (bc, U, edg, f[0]);
	new_bndry->bc.value[0] = f[0];
	new_bndry->bc.value[1] = f[1];
	new_bndry->bc.value[2] = f[2];
	break;

	/* Standard Velocity and Temperature boundaries */

      case 'I': case 'O':
      case 'T': case 'W':
	f[fnum] = 0.;
      case 'V':
	new_bndry = BC_make (bc, U, edg, f[fnum]);
	break;
      case 'A':
        f[fnum] = 0.;
      case 'F':
	new_bndry = BC_make (bc, U, edg, f[fnum]);
	break;

	/* Functional boundaries */

      case 'f':
      case 't':
      case 'v': {
	fpos_t pos;
	
	for (i = 0; i <= fnum; i++) fgets(buf, BUFSIZ, fp);
	new_bndry = BC_make (bc, U, edg, buf);
	
	do {                         
	  fgetpos (fp, &pos);         /* Search through the following    */
	  fgets   (buf, BUFSIZ, fp);  /* lines for the first one without */
	} while (strchr (buf, '='));  /* an '=' (function specification) */
	
	fsetpos (fp, &pos);
	break;
      }
	
	/* Element and Periodic boundaries */

      case 'P':
	if (dim == 2) {
	                          /* Stored for 2-D Flowrate */
	  new_bndry = BC_make (bc, U, edg, f[fnum]);
	  break;
	}
      case 'E':  
	new_bndry = (Bedge *) NULL;
	break;
      default:
	speclib_error("session [%s]: unknown bc type -- %c", SEC_BCS, bc);
	break;
      }

      if (new_bndry) {                       /* Insert the new boundary */
	nbcs++;
	new_bndry->next = bndry_list;
	bndry_list      = new_bndry;
      }
    }
    U = U->next;                     /* next element... */
  }

  /* Return the sorted list */

  return bndry_list != NULL ? bsort(bndry_list, nbcs) : NULL;
}

/* Scan to the next set of BC's and record the position */

static void prescan (FILE *fp, fpos_t *p, char *buf)
{
  do
    findSection (SEC_BCS, buf, fp);
  while
    (strstr(buf, "NO"));

  /* OK, we found a SEC_BCS, but is it the right one?  Check to see if *
   * the NEXT line also includes a SEC_BCS.  If so, position the file  *
   * there.  Otherwise, leave it here.                                 */

  fgetpos (fp, p);        
  if (strstr(fgets(buf,BUFSIZ,fp),SEC_BCS))
    fgetpos (fp, p);                         /* Keep this one */
  else
    fsetpos (fp, p);                         /* Move back one */

  return;
}

/* ------------------------------------------------------------------------ *
 * bsort() - Sort Boundary Conditions                                       *
 *                                                                          *
 * The following function sorts boundary conditions so they are set accord- *
 * ing to precedence and not mesh order.  The following is the "hierarchy"  *
 * of boundary conditions :                                                 *
 *                                                                          *
 *    W          : Walls (zero velocity) are the most binding               *
 *    V          : Velocity boundary conditions (constant)                  *
 *    v          : Velocity function boundary conditions                    *
 *    all others                                                            *
 *                                                                          *
 * This function returns the sordid boundary list (ha ha).                  *
 * ------------------------------------------------------------------------ */

Bedge *bsort(Bedge *bndry_list, int nbcs)
{
  Bedge   **tmp  = (Bedge**) malloc( nbcs * sizeof(Bedge*) );
  Bedge   *bedg  = bndry_list;
  register int i;

  /* Copy the linked list into a regular array and sort it */

  for(i = 0 ; i < nbcs; i++, bedg = bedg->next) tmp[i] = bedg;
  qsort(tmp, nbcs, sizeof(Bedge*), bedgcmp);

  /* Create a new version of the ordered list */
  
  bedg  = (Bedge*) malloc (nbcs * sizeof(Bedge));
  for(i = 0; i < nbcs; i++) {
    memcpy (bedg + i, tmp[i], sizeof(Bedge));
    bedg[i].id   = i;
    bedg[i].next = bedg + i + 1;
  }
  bedg[nbcs-1].next = (Bedge*) NULL;

  /* Free up the space occupied by the old list */

  bndry_list = bedg;
  for(i = 0; i < nbcs; i++) free (tmp[i]);
  free (tmp);

  return bndry_list;
}

/*
 * Boundaries are sorted by their ASCII character values except for
 * types 'D' and 'N', which are sorted by their element ID numbers.
 */

int bedgcmp(const void *b1, const void *b2)
{
  Bedge *bedg_1 = *((Bedge**) b1),
        *bedg_2 = *((Bedge**) b2);
  char  btype_1 = bedg_1 -> type,
        btype_2 = bedg_2 -> type;

  /* Convert {D,N} -> {G,H} to get the ASCII precedence right */

  if (btype_1 == 'D') btype_1 = 'G';
  if (btype_2 == 'D') btype_2 = 'G';
  if (btype_1 == 'N') btype_1 = 'H';
  if (btype_2 == 'N') btype_2 = 'H';
  
  if      (btype_1 > btype_2)              /* Check ASCII code */
    return  1;
  else if (btype_1 < btype_2)
    return -1;
  else                                     /* Check element ID */
    {
      int id_1 = bedg_1 -> elmt -> id,
          id_2 = bedg_2 -> elmt -> id;

      if       (id_1 > id_2)
	return  1;
      else if  (id_1 < id_2)
	return -1;
      else                                 /* Check edge ID    */
	{
	  id_1 = bedg_1 -> edge -> id;
	  id_2 = bedg_2 -> edge -> id;

	  if      (id_1 > id_2)
	    return  1;
	  else if (id_1 < id_2)
	    return -1;
	}
    }

  /* If we get this far there is an error in the boundary conditions */
  
  speclib_error("session [%s]: duplicate boundary\n", SEC_BCS);
  return 0;
}

/* Find a section header in the input file */

char *findSection (char *name, char *buf, FILE *fp)
{
  int   i;
  char *p, lbuf[128];

  while ((p = fgets (buf, BUFSIZ, fp))) {
    for (i = 0; *p && i < sizeof(lbuf); i++, p++)
      lbuf[i] = toupper(*p);
    lbuf[i] = '\0';
    if (strstr (lbuf, name)) 
      break;
  }

  return p;
}

