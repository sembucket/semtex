/* 
 * Field implementation
 *
 * Copyright (c) 1994-98 Ronald Dean Henderson
 *
 * $Id$
 * ----------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>

#include "veclib/veclib.h"

#include "speclib/speclib.h"
#include "speclib/field.h"
#include "speclib/isomesh.h"

/* Private functions */

static Element *make_element  (char type, int nr, int ns, int nz, int k, 
			       Element *link);
static Element *sort_elements (Element *list);
static Edge    *make_edge     (int face, int Nr, int Ns, Edge *link);
static Edge    *sort_edges    (Edge *list);
static void     build_emap    (Field *U);
static double **frame_alloc   (int nr, int ns, int nz, int nel);
static void     frame_project (int nr, int ns, int nz, int nel, double *data,
			       int nx, int ny, double *proj);

/* -------------------------------------------------------------------- * 
 * Field_alloc() - Create a new Field                                   *
 *                                                                      *
 * This function allocates memory for a new Field (element array) and   *
 * initializes each element according to the input parameters for the   *
 * type and resolution.  The following storage areas are allocated:     *
 *                                                                      *
 *     - Mesh storage.  A continuous block of memory to hold the (x,y)- *
 *       coordinates of the spectral element mesh.  No initialization.  *
 *       Memory is allocated as                                         *
 *                                                                      *
 *               MESH (NR,NS,K)      |   U[K].mesh[NS][NR]              *
 *                                                                      *
 *       where MESH = ( xmesh, ymesh ).                                 *
 *                                                                      *
 *     - Data storage.  A continuous (3D) block of memory for storage   *
 *       of the data associated with the element.  The memory block is  *
 *       dimensioned (FORTRAN-style) as:                                *
 *                                                                      *
 *               DATA (NR,NS,K,NZ)   |   U[K].base[NS*NZ][NR]           *
 *                                                                      *
 *       where (NR,NS) are the mesh dimensions in the (r,s)-direction,  *
 *       K is the number of elements, and NZ is the number of "frames". *
 *                                                                      *
 *     - Edge pointers.  This includes the actual structure storage     *
 *       plus memory for the unit normals and area arrays.              *
 *                                                                      *
 *                                                                      *
 *     Total Memory:                                                    *
 *                                                                      *
 *               K  sizeof(Element) +     [ Overhead ]                  *
 *                                                                      *
 *                2                                                     *
 *           K M N  sizeof(double)  +   [ Data Storage ]                *
 *                                                                      *
 *                2                                                     *
 *             K N  sizeof(double)  +   [ Mesh Storage ]                *
 *                                                                      *
 *                                                                      *
 *           4 K N (sizeof(Edge) + sizeof(int) + 3 sizeof(double) )     *
 *                                                                      *
 * -------------------------------------------------------------------- */

static int rescale (Element *U)
{
  int    nrns = U->nr * U->ns;
  double sx   = dparam("XSCALE"),
         sy   = dparam("YSCALE"),
         xo   = dparam("XSHIFT"),
         yo   = dparam("YSHIFT");
  double *x   = U->xmesh[0],
         *y   = U->ymesh[0];
  
  if (sx == 1. && sy == 1. && xo == 0. && yo == 0.) return 0;

  dscal (nrns, sx, x, 1)      ;  dscal (nrns, sy, y, 1);
  dsadd (nrns, xo, x, 1, x, 1);  dsadd (nrns, yo, y, 1, y, 1);

  return 1;
}

Field *Field_build (FILE *fp, int nr, int ns, int nz, char type)
{
  int nel;
  int i, k;
  char buf[BUFSIZ];
  Curve BndryShape[_MAX_NEL];
  double **p;

  Field   *u    = (Field*) NULL;
  Element *elmt = (Element *) NULL;

  if (!findSection("MESH", buf, fp))
    speclib_error("session: can't find a MESH in the input file\n");

  if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &nel) != 1)
    speclib_error("session: can't read the # of elements\n");

  /* Check limits */

  if (nr > _MAX_NORDER || ns > _MAX_NORDER)
    speclib_error("session: NORDER must be less than %d", _MAX_NORDER);
  if (nel > _MAX_NEL)
    speclib_error("session: # elements must be less than %d", _MAX_NEL);
  
  /* Create a continuous block of element structures */

  for (k = 0; k < nel; k++)
    elmt = make_element (type, nr, ns, nz, k, elmt);
  elmt = sort_elements (elmt);
  u    = elmt;

  /* * * *    M E S H   S T O R A G E   * * * */

  p = frame_alloc (nr, ns, 2, nel);
  for (k = 0; k < nel; k++, p += ns)
    elmt[k].xmesh = p;
  for (k = 0; k < nel; k++, p += ns)
    elmt[k].ymesh = p;

  /* * * *    D A T A  S T O R A G E   * * * */

  p = frame_alloc (nr, ns, nz, nel);
  for (k = 0; k < nel; k++, p += ns)
    elmt[k].base = elmt[k].field = p;

  build_emap(elmt);  /* Build the B/I scatter array    */
  coef(nr);          /* r-direction spectral operators */
  coef(ns);          /* s-direction spectral operators */

  /* Read in mesh information.  Initially we assume that every element   *
   * has straight sides.  This can be overridden when curved-side spec-  *
   * ifications are loaded after the basic element specs have been read. */
  
  for (k = 0; k < nel; k++) {
    fgets(buf, BUFSIZ, fp);  
    for (i = 0; i < 4; i++)
      fscanf (fp, "%lf", (*u[k].xmesh) + u[k].edges[i].start);
    for (i = 0; i < 4; i++)
      fscanf (fp, "%lf", (*u[k].ymesh) + u[k].edges[i].start);
    for (i = 0; i < 4; i++)
      BndryShape[k].type[i] = T_Strait;    /* Initially, at least... */
    fgets (buf, BUFSIZ, fp);
  }

  /* Read the curved side information */

  i = 2; while (i--) fgets (buf, BUFSIZ, fp);
  sscanf(buf, "%d", &k);   /* k = number of curved sides */

  if (k) {
    int iside;
    int iel;
    char type;
    char *p;
    Curve *curv;

    while (k--) {
      fgets(buf, BUFSIZ, fp);
      if (sscanf(buf, "%1d%d", &iside, &iel) != 2) {
	speclib_error
	  ("ReadMesh -- unable to process the following curve:\n%s", buf);
      }
      
      iel--; iside--; p = buf + strlen(buf); 
      while (isspace(*--p))
	{ /* empty while */ }
      type = toupper(*p);
      curv = &BndryShape[iel];

      switch (type) {
      case 'C':                                       /* Circular arc */
	curv->type[iside] = T_Arc;
	sscanf(buf, "%*1d%*d%lf", &(curv->info[iside].arc.radius));
	break;

      case 'W':                                       /* Sinusoidal wave */
	curv->type[iside] = T_SinWave;
	sscanf(buf, "%*1d%*d%lf%lf",
	       &(curv->info[iside].wave.wavenum), 
	       &(curv->info[iside].wave.magnitude));
	break;

      case 'P': {                                     /* Parametric form */
	char *c;
	curv->type[iside] = T_Parametric;
	if (!(c = strchr(buf, '='))) {
	  speclib_error
	    ("ReadMesh -- invalid parametric form:\n%s", buf);
	}

	strncpy(curv->info[iside].map.function = (char *) malloc (p - c),
		c + 1, p - c - 1);
	break;
      }

      case 'F': {                                     /* Read from a file */
	char fname[FILENAME_MAX];
	sscanf(buf,"%*1d%*d%s", fname);	
	strcpy( curv->info[iside].file.name = 
		(char*) malloc(strlen(fname)+1), fname);
	
	/* Check to see if any shifting parameters have been defined. *
	 * The sscanf() should simply fail if it can't convert two    *
	 * additional double's in the string.                         */
	
	sscanf(buf,"%*1d%*d%*s%lf%lf", 
	       &(curv->info[iside].file.xoffset),
	       &(curv->info[iside].file.yoffset));
	  
	curv->type[iside] = T_File;
	break;
      }
	
      default: 
	speclib_error("session: unknown curve type -- %c", type);
	break;
      }
    }
  }
  
  /* Now the geometric description for every edge in the mesh is stored *
   * in the array BndryShape[], so we just need to generate the mesh    *
   * and mapping factors for each one.                                  */

  for (k = 0; k < nel; k++) {
    genxy   (&u[k], &BndryShape[k]); /* generate collocation points */
    rescale (&u[k]);                 /* rescale and shift origin */
    geomap  (&u[k]);                 /* compute isoparametric map */
    normals (&u[k]);                 /* compute outward normal vectors */
  }

  return u;
}

Field *Field_alloc (char type, int nr, int ns, int nz, int nel) 
{
  Element *elmt = NULL;
  double **p;
  int k;

  /* Check limits */

  if (nr > _MAX_NORDER || ns > _MAX_NORDER)
    speclib_error("session: NORDER must be less than %d", _MAX_NORDER);
  if (nel > _MAX_NEL)
    speclib_error("session: # elements must be less than %d", _MAX_NEL);
  
  /* Create a continuous block of element structures */

  for (k = 0; k < nel; k++)
    elmt = make_element (type, nr, ns, nz, k, elmt);
  elmt = sort_elements (elmt);

  /* * * *    M E S H   S T O R A G E   * * * */

  p = frame_alloc (nr, ns, 2, nel);
  for (k = 0; k < nel; k++, p += ns)
    elmt[k].xmesh = p;
  for (k = 0; k < nel; k++, p += ns)
    elmt[k].ymesh = p;

  /* * * *    D A T A  S T O R A G E   * * * */

  p = frame_alloc (nr, ns, nz, nel);
  for (k = 0; k < nel; k++, p += ns)
    elmt[k].base = elmt[k].field = p;

  build_emap(elmt);  /* Build the B/I scatter array    */
  coef(nr);          /* r-direction spectral operators */
  coef(ns);          /* s-direction spectral operators */

  return elmt;
}

/* free memory for the Field data structure and data storage */

void Field_free (Field *u)
{
  if (u->base) {
    free(*u->base);
    free( u->base);
  }
  free(u);
}

/* ------------------------------------------------------------------------ *
 * BLAS-type operations                                                     *
 * ------------------------------------------------------------------------ */

double Field_dot (const Field *u, const Field *v) { 
  return ddot(Field_size(u), *u->base, 1, *v->base, 1); 
}

void Field_axpy (double d, const Field *u, Field *v) { 
  daxpy(Field_size(u), d, *u->base, 1, *v->base, 1); 
}

void Field_copy (const Field *u, Field *v) { 
  dcopy(Field_size(u), *u->base, 1, *v->base, 1); 
}

void Field_swap (Field *u, Field *v) { 
  dswap(Field_size(u), *u->base, 1, *v->base, 1); 
}

double Field_nrm2 (const Field *u) { 
  return dnrm2(Field_size(u), *u->base, 1); 
}

double Field_asum (const Field *u) { 
  return dasum(Field_size(u), *u->base, 1); 
}

void Field_scal (double d, Field *u) { 
  dscal(Field_size(u), d, *u->base, 1); 
}

double Field_amax (const Field *u) { 
  return fabs((*u->base)[idamax(Field_size(u), *u->base, 1)]);
}

double Field_max (const Field *u) {
  return (*u->base)[idmax(Field_size(u), *u->base, 1)];
}

double Field_min (const Field *u) {
  return (*u->base)[idmin(Field_size(u), *u->base, 1)];
}

void Field_mul (const Field *u, const Field *v, Field *w) {
  dvmul(Field_size(u), *u->base, 1, *v->base, 1, *w->base, 1);
}

void Field_vvtp (const Field *u, const Field *v, Field *w) {
  dvvtvp(Field_size(u), *u->base, 1, *v->base, 1, *w->base, 1, *w->base, 1);
}

/* Function norms */

double Field_L2 (const Field *u) 
{
  const int npts = u->nr * u->ns;
  double sum = 0.;

  Element *elmt;
  int i;

  for (elmt = Field_head(u); elmt; elmt = elmt->next) {
    const double *mass = *elmt->mass;
    const double *valu = *elmt->field;
    for (i = 0; i < npts; i++)
      sum += mass[i] * valu[i] * valu[i];
  }

  return sqrt(sum);
}

double Field_H1 (const Field *u) 
{
  Field *tmp = Field_dup(u);
  double comp[3];

  comp[0] = Field_L2(u);

  Field_dx (u, tmp);
  comp[1] = Field_L2(tmp);

  Field_dy (u, tmp);
  comp[2] = Field_L2(tmp);

  Field_free(tmp);
  return dnrm2(3,comp,1);
}

/* Integral of a field over the domain */

double Field_integral (Field *u)
{
  double sum = 0.;

  const int npts = FIELD_NR(u)*FIELD_NS(u)*FIELD_NELMT(u);
  const int nrns = FIELD_NR(u)*FIELD_NS(u);
  const int nz   = FIELD_NZ(u);

  Element *elmt;
  int k;

  for (elmt = Field_head(u); elmt; elmt = elmt->next) {
    const int id = elmt->id;
    for (k = 0; k < nz; k++)
      sum += ddot(nrns, *u[id].mass, 1, *u[id].base+k*npts, 1);
  }

  if (nz > 1) sum *= dparam("LZ")/nz;

  return sum;
}

/* -----------------  P R I V A T E    F U N C T I O N S  ----------------- */

static Element *make_element
 (char type, int nr, int ns, int nz, int id, Element *link)
{
  Edge*    edge = NULL;
  Element* elmt = (Element *) calloc (1, sizeof(Element));

  /* Build a set of edges for the element */

  int  face;
  for (face = 0; face < 4; face++)
    edge = make_edge (face, nr, ns, edge);

  /* Initialize the data for the element */

  elmt -> id           =  id;
  elmt -> type         =  type;
  elmt -> nr           =  nr;
  elmt -> ns           =  ns;
  elmt -> nz           =  nz;
  elmt -> frame        =  0;
  elmt -> edges        =  sort_edges (edge);
  elmt -> next         =  link;

  /* Geometry Factors to map this element's geometry to the standard  *
   * domain (-1,1) x (-1,1) are computed by genxy(), as well as the   *
   * factors for the unit normals and the surface integral array.     */
  
  return elmt;
}

static Edge *make_edge (int edge_no, int nr, int ns, Edge *link)
{
  const int np    = (edge_no & 1) ? ns : nr;
  Edge*     edge  = (Edge*) malloc(sizeof(Edge));
  double*   local = dvector (0, 3 * np - 1);

  edge->id     =  edge_no;
  edge->np     =  np;
  edge->unx    =  local;
  edge->uny    =  local + np;
  edge->area   =  local + np * 2;
  edge->next   =  link;

  switch (edge_no) {
  case 0:
    edge->bindex = 0;
    edge->start  = 0;
    edge->skip   = 1;
    break;
  case 1:
    edge->bindex = nr - 1;
    edge->start  = nr - 1;
    edge->skip   = nr;
    break;
  case 2:
    edge->bindex =  nr + ns  - 2;
    edge->start  = (nr * ns) - 1;
    edge->skip   = -1;
    break;
  case 3:
    edge->bindex = (nr << 1) + ns - 3;
    edge->start  = -nr * (1 - ns);
    edge->skip   = -nr;
    break;
  default:
    break;
  }

  return edge;
}

static Element *sort_elements (Element *linkedlist)
{
  Element *list = (Element*) NULL;

  if (linkedlist) {
    int i;

    const int nel = Field_count(linkedlist);
    Element *s = linkedlist;

    list = (Element*) calloc(nel, sizeof(Element));
    
    for (i = 0; i < nel; i++) {
      Element *tmp = s->next;
      memcpy (list + i, s, sizeof(Element));
      free   (s);
      s = tmp;
    }

    qsort (list, nel, sizeof(Element), idcomp);
    for (i = 0; i < nel-1; i++) 
      list[i].next = &list[i+1];
    list[i].next = (Element*) NULL;
  }

  return list;
}

static Edge *sort_edges (Edge *list)
{
  if (list) {
    const int nedge = 4;
    Edge* s = list;
    Edge* n = (Edge*) calloc (nedge, sizeof(Edge));
    
    int  i;
    for (i = 0; i < nedge; i++) {
      Edge* t = s->next;
      memcpy (n + i, s, sizeof(Edge));
      free   (s);
      s = t;
    }
    
    qsort (list = n, nedge, sizeof(Edge), idcomp);
    for (i = 0; i < nedge-1; i++) 
      list[i].next = list + i + 1;
    list[i].next = (Edge*) NULL;
  }
  
  return list;
}

/*
 * Generic id-comparison function ... just requires that the first member
 * be an integer which will be the basis for the sort.
 */

int idcomp (const void *p1, const void *p2)
{
  int  id1   = *((int*) p1);
  int  id2   = *((int*) p2);

  if (id1 > id2) return  1;
  if (id1 < id2) return -1;

  return 0;
}

/* Allocate memory for a multi-frame field */

static double **frame_alloc (int nr, int ns, int nz, int nel)
{
  double*  dblock = NULL;
  double** pblock = (double**) malloc (ns * nz * nel * sizeof(double*));

  int ntot = nr*ns*nz*nel;  /* Total storage for all frames */
  int i;

#ifdef PARALLEL
  /* The parallel transpose used by Prism requires each field to    *
   * allocate "p" additional words such that nr*ns*nel+p is a mult- *
   * iple of the number of processors.                              */

  if (nz > 1) ntot += ((nr*ns*nel)%comm_size())*nz;
#endif

  dblock = (double*) calloc(ntot, sizeof(double));

  if (!pblock || !dblock)
    speclib_error ("frame_alloc(): out of memory");
  
  for (i = 0; i < ns*nel*nz; i++, dblock += nr)
    pblock[i] = dblock;

  return pblock;
}

/* ------------------------------------------------------------------------ *
 * build_emap() - Compute the Elemental Boundary/Interior Scatter Vector    *
 *                                                                          *
 * This function computes the index map which takes a sequentially ordered  *
 * matrix into a boundary/interior ordered matrix.  This version assumes    *
 * that all elements have the same scatter properties.                      *
 * ------------------------------------------------------------------------ */

static void build_emap (Field *U)
{
  const int nr  = U->nr;
  const int ns  = U->ns;
  const int nel = Field_count (U);

  int  pos   = 0;
  int* emap  = ivector (0, nr * ns - 1);
  register int i, j, k;

  /* Boundary scatter vector */

  for (j = 0, i = 0; j < nr; j++)
    emap[pos++] = i * nr + j;
  for (j = nr-1, i = 1; i < ns; i++)
    emap[pos++] = i * nr + j;
  for (j = nr-2, i = ns-1; j >= 0; j--)
    emap[pos++] = i * nr + j;
  for(j = 0, i = ns-2;i > 0;i--)
    emap[pos++] = i * nr + j;

  /* Interior scatter vector */

  for (i = 1;i < ns-1; i++)
    for (j = 1;j < nr-1; j++)
      emap[pos++] = i * nr + j;

  /* Now assign this emap to every element in the mesh */

  for (k = 0; k < nel; k++) 
    U[k].emap = emap;
  
  return;
}


/* Project frame data from DATA(nr,ns,nz,K) to PROJ(nx,ny,nz,K) */

static void frame_project
  (int nr, int ns, int nz, int nel, double *data,
   int nx, int ny,                  double *proj)
{

  double *zr, *zs, *zx, *zy;
  double tmp [_MAX_NORDER * _MAX_NORDER];
  register int k, m;

  /* Allocate space for the interpolation matrices */

  const int nrns = nr * ns;
  const int nxny = nx * ny;

  double **imr   = dmatrix (0, nx-1, 0, nr-1);
  double **itmr  = dmatrix (0, nr-1, 0, nx-1);
  double **ims   = dmatrix (0, ny-1, 0, ns-1);
  double **itms  = dmatrix (0, ns-1, 0, ny-1);

  /* Compute the GLL points */

  coef (nr); getops (nr, &zr, 0, 0, 0);
  coef (ns); getops (ns, &zs, 0, 0, 0);
  coef (nx); getops (nx, &zx, 0, 0, 0);
  coef (ny); getops (ny, &zy, 0, 0, 0);

  /* Compute the interpolation matrices */

  igllm (imr, itmr, zr, zx, nr, nx);
  igllm (ims, itms, zs, zy, ns, ny);

  /* ----- Project ----- */

  for (m = 0; m < nz; m++) {
    for (k = 0; k < nel; k++, data += nrns, proj += nxny) {
      dgemm ('N', 'N', nr, ny, ns, 1., data, nr, *ims, ns, 0., tmp,  nr);
      dgemm ('T', 'N', nx, ny, nr, 1., *imr, nr,  tmp, nr, 0., proj, nx);
    }
  }

  free_dmatrix (imr,  0, 0);
  free_dmatrix (itmr, 0, 0);
  free_dmatrix (ims,  0, 0);
  free_dmatrix (itms, 0, 0);
}

/* ------------------------------------------------------------------------- *
 * Field_[size functions]                                                    *
 *                                                                           *
 * The following collection of functions return various size parameters for  *
 * a Field.  There are four basic size parameters: NR, NS (polynomial order  *
 * in each direction), NZ (number of 3D frames), and NEL (number of ele-     *
 * ments).  The functions below return various combinations of these:        *
 *                                                                           *
 *   count       number of elements in the mesh, = NEL                       *
 *                                                                           *
 *   npts        number of grid points in the mesh, = NR*NS*NZ*NEL           *
 *                                                                           *
 *   frameSize   number of grid points in a single frame, = NR*NS*NEL        *
 *                                                                           *
 *   frameCount  number of frames in a 3D Field, = NZ                        *
 * ------------------------------------------------------------------------- */

int Field_count (const Field *u)
{
  static int    count = 0;
  static Field* last  = NULL;

  if (u != last) {
    Element* elmt = Field_head(u);
    for (count = 0; elmt != NULL; elmt = elmt->next)
      count++;
    last = (Field*) u;
  }

  return count;
}

int Field_frameSize (const Field *u)
{ return u->nr * u->ns * Field_count(u); }

int Field_frameCount (const Field *u)
{ return u->nz; }

int Field_npts (const Field *u)
{ return u->nr * u->ns * u->nz * Field_count(u); }

int Field_size (const Field *u)
{ return u->nr * u->ns * u->nz * Field_count(u); }

/* ------------------------------------------------------------------------- */

int Field_setFrame (Field *u, int m) 
{
  const int nel  = Field_count(u);
  const int skip = u->ns * nel * m;
  int k;

  for (k = 0; k < nel; k++) {
    u[k].field = u[k].base + skip;
    u[k].frame = m;
  }

  return m;
}

int Field_setFrameMulti (int m, int nfields, ...) 
{
  va_list  ap;
  va_start(ap, nfields);

  while (nfields--) {
    Field *u = (Field*) va_arg(ap,Field*);
    Field_setFrame(u, m);
  }

  va_end(ap);
  return m;
}

/* ------------------------------------------------------------------------- */

char Field_setType (Field *U, char type)
{ return U->type = type; }

char Field_getType (const Field *U)
{ return U->type; }

double* Field_base (const Field *U)
{ return *U->base; }

Element* Field_head (const Field *U)
{ return (Element*) U; }

/* ------------------------------------------------------------------------ *
 * Field_aux() - Create an auxillary field                                  *
 *                                                                          *
 * Create a new Field from an existing one, using the "parent" to define    *
 * the geometry of the finite element domain.  If the new field is compat-  *
 * ible (has the same polynomial order), then it receives pointers to the   *
 * mesh and geometric factors of the parent.  Otherwise, the mesh is        *
 * "projected" and new geometric factors are computed.                      *
 *                                                                          *
 * If the number of frames in the parent and child fields are the same,     *
 * the parent's solution is "projected" to the child for initialization.    *
 * Otherwise, the new field is initialized to zero.                         *
 * ------------------------------------------------------------------------ */

Field *Field_aux (const Field *parent, int nr, int ns, int nz, char type)
{
  int k;
  Field* child;
  extern void geomap (Element*);

  /* Compatibility Flags */

  const int compatible_rs = ((parent->nr == nr) && (parent->ns == ns));
  const int compatible_z  = ((parent->nz == nz));
  
  /* Allocate space for the new Field */

  const int nel = Field_count(parent);

  if (compatible_rs) {
    double** data  = frame_alloc (nr, ns, nz, nel);

    child = (Field*) 
      memcpy(malloc(nel*sizeof(Element)), parent, nel*sizeof(Element));

    /* Initialize type and pointers to field, base, and "next" */

    for (k = 0; k < nel; k++, data += ns) {
      child[k].type  = type;
      child[k].field = child[k].base = data;
    }
    for (k = 0; k < (nel-1); k++)
      child[k].next  = &child[k+1];

  } else {

    /* Allocate everything (edges, etc.) from scratch and   //
    // interpolate the (x,y)-mesh.                          */

    child = Field_alloc (type, nr, ns, nz, nel);

    frame_project (parent->nr, parent->ns, 1, nel, *parent->xmesh,
		   child ->nr, child ->ns,         *child ->xmesh);
    frame_project (parent->nr, parent->ns, 1, nel, *parent->ymesh,
		   child ->nr, child ->ns,         *child ->ymesh);
    
    for (k = 0; k < nel; k++) {       /* Generate new mapping factors */
      geomap (& child[k]);
      normals(& child[k]);
    }
  }

  /* If the two fields are (3D) compatible, project the "parent" to //
  // the child.  Otherwise, set the child's initial value to zero.  */

  if (compatible_z)

    frame_project (parent->nr, parent->ns, nz, nel, *parent->base,
		   child ->nr, child ->ns,          *child ->base);

  else dzero (nr * ns * nz * nel, *child->base, 1);

  return child;
}

/* Convenient wrapper around Field_aux() */

Field *Field_dup (const Field *u) { 
  return Field_aux(u, u->nr, u->ns, u->nz, u->type); 
}

/* ------------------------------------------------------------------------ *
 * Field_project()                                                          *
 *                                                                          *
 * Projects a field to a higher or lower polynomial space using spectral    *
 * interpolation.  For multi-frame fields, the projection is done one frame *
 * at a time.  If the input (U) and output (V) fields have a different      *
 * number of frames, it exits with an error condition.                      *
 * ------------------------------------------------------------------------ */

int Field_project (const Field *U, Field *V)
{  
  const int nz  = U->nz;
  const int nel = Field_count(U);

  /* Make sure they're compatible */

  if (nz != V->nz || nel != Field_count(V))
    speclib_error("Field_project: U and V have different frame dimensions");

  
  /* Copy if they're FULLY compatible... */

  if ((U->nr == V->nr) && (U->ns == V->ns))
    dcopy (U->nr * U->ns * nz * nel, *U->base, 1, *V->base, 1);

  /* ...otherwise, do the projection */
  
  else 

    frame_project (U->nr, U->ns, nz, nel, *U->base,
		   V->nr, V->ns,          *V->base);

  return 0;
}

/* ------------------------------------------------------------------------ *
 * Field_set() - Evaluate a function over an element field                  *
 *                                                                          *
 * Takes a Field and a character string which it then evaluates at each     *
 * point in the field.                                                      *
 * ------------------------------------------------------------------------ */

void Field_set (Field *u, const char *expr)
{
  const int n  = FIELD_NR(u)*FIELD_NS(u)*FIELD_NELMT(u);
  const int nz = FIELD_NZ(u);

  if (nz==1) {
    vector_def ("x y", expr);
    vector_set (n, *u->xmesh, *u->ymesh, *u->field);
  } else {
    int k;
#ifdef PARALLEL
    /* The parallel version uses local z coordinates.  If we shift the start *
     * of the array "z", then z[0:nz-1] contains the correct local values.   */

    double *z = zmesh(nz*comm_size());
    z += nz*comm_rank();
#else
    double *z = zmesh(nz);
#endif

    vector_def ("x y", expr);
    
    for (k = 0; k < nz; k++) {
      Field_setFrame (u, k);
      scalar_set ("z", z[k]);
      vector_set (n, *u->xmesh, *u->ymesh, *u->field);
    }

#ifdef PARALLEL
    z -= nz*comm_rank();
#endif
    free(z);
  }
}

void Field_set_3D (Field *U, const char *expr) {
  Field_set(U,expr);
}

/* Compute the Fourier transform of a field in place.  The direction of the  *
 * transform is dir = -1 (Fourier) or dir = +1 (Physical).                   */

void Field_FFT (Field *u, int dir)
{
  const int m = u->nr * u->ns * Field_count(u);
  const int n = u->nz;

  int i, j, ierr;

  const int irev  = 1;
  const int ireal = 1;

  static int     ncmplx = 0;
  static int    *bitrev = NULL;
  static double *factor = NULL;
  
  double **tmp = dmatrix (0, m-1, 0, n+1);

  /* Initialize the FFT */

  if (ncmplx != n/2) {
    if (factor) free (factor);
    if (bitrev) free (bitrev);

    ncmplx = n/2;
    factor = (double*) malloc(6*ncmplx*sizeof(double));
    bitrev = (int*) malloc(6*ncmplx*sizeof(int));
    fftdf(*tmp, ncmplx, 1, 1, 1, 0, factor, irev, bitrev, ierr, ireal);
  }

  for (i = 0; i < m; i++)
    for (j = n; j < n+2; j++)
      tmp[i][j] = 0.;

  /* Begin transform */

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      tmp[j][i] = (*u->base)[i*m + j];

  fftdf(*tmp,ncmplx,1,m,ncmplx+1,dir,factor,irev,bitrev,ierr,ireal);

  if (dir == -1)
    dscal (m*(n+2), 1./(4.*ncmplx), *tmp, 1);

  /* Transpose back */

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      (*u->base)[i*m+j] = tmp[j][i];

  free_dmatrix(tmp, 0, 0);
}

/* ------------------------------------------------------------------------- *
 * Field_grad() - Scalar gradient                                            *
 *                                                                           *
 * This function computes the gradient of a scalar field.  Since the deriva- *
 * tives are computed using the chain rule, it's more efficient to compute   *
 * both x and y derivatives at the same time (if you need them).             *
 *                                                                           *
 * If only one of the two components (x,y) are desired, a NULL pointer can   *
 * be passed for the undesired scalar field and it will be ignored.  Since   *
 * the derivative operators are defined with respect to a square (r,s) coor- *
 * dinate sytem, the derivatives in (x,y)-space are given by:                *
 *                                                                           *
 *               dU/dx = dU/dr * dr/dx + dU/ds * ds/dx                       *
 *               dU/dy = dU/dr * dr/dy + dU/ds * ds/dy                       *
 *                                                                           *
 * A single scalar field may be used for both input and output if only one   *
 * component is needed.                                                      *
 * ------------------------------------------------------------------------- */

void Field_grad (const Field *U, Field *X, Field *Y)
{
  const int nr   = U->nr;
  const int ns   = U->ns;
  const int nrns = nr * ns;
  const int nel  = Field_count(U);
  const int ntot = nrns * nel;
  int i, k;

  tempVector (Ur, ntot);
  tempVector (Us, ntot);

  double  **dr, **ds;
  getops (nr, 0, 0, &dr, 0);
  getops (ns, 0, 0, &ds, 0);

  /* Compute (r,s) partials */

  dgemm ('T','N', nr, ns*nel, nr, 1., *dr, nr, *U->field, nr, 0., Ur, nr);
  for (k = 0; k < nel; k++)
    dgemm ('N','N', nr, ns, ns, 1., *U->field+k*nrns, nr, *ds, ns, 
	   0., Us+k*nrns, nr);

  /* Compute -X- derivatives */

  if (X)
    for (k = 0; k < nel; k++)
      if (X[k].sx)
	for (i = 0; i < nrns; i++)
	  (*X[k].field)[i] = (Ur + k*nrns)[i] * (*X[k].rx)[i] +
	                     (Us + k*nrns)[i] * (*X[k].sx)[i];
      else
	for (i = 0; i < nrns; i++)
	  (*X[k].field)[i] = (Ur + k*nrns)[i] * (*X[k].rx)[i];

  
  /* Compute -Y- derivatives */
  
  if (Y)
    for (k = 0; k < nel; k++)
      if (Y[k].ry)
	for (i = 0; i < nrns; i++)
	  (*Y[k].field)[i] = (Ur + k*nrns)[i] * (*Y[k].ry)[i] +
	                     (Us + k*nrns)[i] * (*Y[k].sy)[i];
      else
	for (i = 0; i < nrns; i++)
	  (*Y[k].field)[i] = (Us + k*nrns)[i] * (*Y[k].sy)[i];


  freeVector (Ur); 
  freeVector (Us);
}

void Field_dx (const Field *u, Field *dx) { 
  Field_grad (u, dx, (Field*) NULL); 
}

void Field_dy (const Field *u, Field *dy) { 
  Field_grad (u, (Field*) NULL, dy); 
}

/* Note: this function computes D[u,z] for u in physical space, and returns  *
 * du in physical space.  Use _dzhat if you have u in Fourier space and want *
 * the Fourier modes of its derivative.                                      */

void Field_dz (const Field *u, Field *du) 
{
  const int m = FIELD_NR(u)*FIELD_NS(u)*Field_count(u);
  const int n = FIELD_NZ(u);

  const double beta = dparam("BETA");
  int k;  

  Field_copy (u, du);
  Field_FFT (du, -1);

  /* Zero the first two planes */

  dzero (m+m, *du->base, 1);

  /* For the remaining planes we have:
   *
   *    du/dz = I beta_k (ur + I ui) = beta_k (-ui + I ur)
   */

  for (k = 2; k < n; k += 2) {
    double scal = beta * (k >> 1);
    double *ur  = *du->base + k*m;
    double *ui  = *du->base + k*m + m;

    dswap (m, ur, 1, ui, 1);
    dscal (m, -scal, ur, 1);
    dscal (m,  scal, ui, 1);
  }

  Field_FFT (du, +1);
}

void Field_dzhat (const Field *u, Field *du) 
{
  const int m = FIELD_NR(u)*FIELD_NS(u)*Field_count(u);
  const int n = FIELD_NZ(u);

  const double beta = dparam("BETA");
  int k;

  Field_copy (u, du);

  /* Zero the first two planes */

  dzero (m+m, *du->base, 1);

  /* For the remaining planes we have:
   *
   *    du/dz = I beta_k (ur + I ui) = beta_k (-ui + I ur)
   */

  for (k = 2; k < n; k += 2) {
    double scal = beta * (k >> 1);
    double *ur  = *du->base + k*m;
    double *ui  = *du->base + k*m + m;

    dswap (m, ur, 1, ui, 1);
    dscal (m, -scal, ur, 1);
    dscal (m,  scal, ui, 1);
  }
}

/* ------------------------------------------------------------------------- */

int Field_show (const Field *u)
{
  const int nr  = FIELD_NR(u);
  const int ns  = FIELD_NS(u);
  const int nz  = FIELD_NZ(u);
  const int nel = FIELD_NELMT(u);

  Element *elmt;
  int i, j, k;

  printf ("\nField %c:\n", FIELD_TYPE(u));

  for (elmt = Field_head(u); elmt; elmt = elmt->next) {
    const double *uptr = *elmt->base;

    printf ("Element %d:\n", ELEMENT_ID(elmt));
    for (k = 0; k < nz; k++) {
      for (i = 0; i < ns; i++) {
	for (j = 0; j < nr; j++) {
	  printf (" %#10.6g", uptr[k*(nr*ns*nel) + i*ns + j]);
	}
	printf ("\n");
      }
      printf ("--\n");
    }
  }

  return 0;
}

