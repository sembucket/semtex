///////////////////////////////////////////////////////////////////////////////
// isowire.C: generate and manipulate triangulated wireframe
// approximations to isosurfaces.
//
// Much of this code is reproduced verbatim from sview_2.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sview.h>

const int VERT_MAX = 100000;
const int POLY_MAX = 100000;

// -- File-scope polygon globals.

static int    NUM_VERTICES;
static int    NUM_POLYGONS;
static float* XVERTICES = NULL;
static float* YVERTICES = NULL;
static float* ZVERTICES = NULL;
static int*   POLYGONS  = NULL;
static int*   first_poly;
static int*   last_poly;
static int    goff[2][2][2];
static int    too_many_polys;

#define v_mag(v) (sqrt((v)[0]*(v)[0]+(v)[1]*(v)[1]+(v)[2]*(v)[2]))

// -- Local prototypes.

static int   add_vertex   (int, float*, int);
static int   equal_vertex (float*, int);
static float tri_lin_int  (float [3], int, int, int, float*);
static int   add_polygon  (float* [3], int, int, int, int, int, int,
			   float*, float*, float*);
static float calc_norm    (float*, int*, float*);

// -- NCSA "marching cubes" routine prototypes.

static void calc_index_and_temps (register float*, int, int, int, int, int,
				  int, register float, int*);
static void get_cell_verts       (int, int, int, int, float, float [13][3]);
static void get_cell_polys       (int, float [13][3], int, int, int, int,
				  int, int, float*, float*, float*);


Iso* makeSurf (int     nel  ,	// -- Number of spectral elements.
	       float** xgrid,	// -- Array of x mesh locations for each elmt.
	       float** ygrid,	// -- ...      y ...
	       float** zgrid,	// -- ...      z ...
	       float** data ,	// -- Array of data for each element.
	       int*    idim ,	// -- I dimension of each element.
	       int*    jdim ,	// -- J dimension of each element.
	       int*    kdim ,	// -- K dimension of each element.
	       float   iso  ,	// -- Isosurface value.
	       int     both )	// -- Flag for isosurfacing +/- values.
/* ------------------------------------------------------------------------- *
 * Generate a wire-frame interpolant to data at value iso.
 *
 * In this version of the code we assume that there will not be more
 * than VERT_MAX vertices or POLY_MAX polygons - if there are, don't
 * make wire-frame.
 * ------------------------------------------------------------------------- */
{
  char routine[] = "makeSurf";
  register int   i, j, k, l;
  register float *vect, mag;
  float          crossings[13][3];
  int            index, npoints, nzero;
  Iso*           S = 0;

  too_many_polys = 0;
  NUM_VERTICES   = 0; 
  NUM_POLYGONS   = 0;

  if (!(XVERTICES = (float*) calloc (VERT_MAX, sizeof (float)))) {
    message (routine, "insufficient memory, no surface made", WARNING);
    return 0;
  }
  if (!(YVERTICES = (float*) calloc (VERT_MAX, sizeof (float)))) {
    message (routine, "insufficient memory, no surface made", WARNING);
    free (XVERTICES);
    return 0;
  }
  if (!(ZVERTICES = (float *) calloc (VERT_MAX, sizeof (float)))) {
    message (routine, "insufficient memory, no surface made", WARNING);
    free (YVERTICES);
    free (XVERTICES);
    return 0;
  }
  if (!(POLYGONS = (int*) calloc (3*POLY_MAX, sizeof (int)))) {
    message (routine, "insufficient memory, no surface made", WARNING);
    free (ZVERTICES);
    free (YVERTICES);
    free (XVERTICES);
    return 0;
  }
  
  /* -- Generate polygons - loop through each macro element. */
    
  for (l = 0; l < nel; l++) {
    
    npoints = idim[l] * jdim[l] * kdim[l];
    
    /* -- Set the element offset array for grid interpolation. */
    
    goff[0][0][0] = 0;
    goff[1][0][0] = 1;
    goff[0][1][0] = idim[l];
    goff[1][1][0] = idim[l]+1;
    goff[0][0][1] = idim[l]*jdim[l];
    goff[1][0][1] = idim[l]*jdim[l]+1;
    goff[0][1][1] = idim[l]*jdim[l]+idim[l];
    goff[1][1][1] = idim[l]*jdim[l]+idim[l]+1;
    
    if (!(first_poly = (int*) calloc (npoints, sizeof (int)))) {
      message (routine, "insufficient memory, no surface made", WARNING);
      free (POLYGONS);
      free (ZVERTICES);
      free (YVERTICES);
      free (XVERTICES);
      return 0;
    }
    if (!(last_poly = (int*) calloc (npoints, sizeof(int)))) {
      message (routine, "insufficient memory, no surface made", WARNING);
      free (first_poly);
      free (POLYGONS);
      free (ZVERTICES);
      free (YVERTICES);
      free (XVERTICES);
      return 0;
    }

    for (k = 0; k < kdim[l] - 1; k++)
      for (j = 0; j < jdim[l] - 1; j++)
	for (i = 0 ; i < idim[l] - 1; i++) {
	  last_poly [(k*jdim[l]+j)*idim[l]+i] = 0;
	  first_poly[(k*jdim[l]+j)*idim[l]+i] = NUM_POLYGONS + 1;
	  calc_index_and_temps (data[l], i, j, k, 
				idim[l], jdim[l], kdim[l],
				iso, &index);
	  if (index) {
	    get_cell_verts (index, i, j, k, iso, crossings);
	    get_cell_polys (index, crossings, i, j, k,
			    idim[l],  jdim[l],  kdim[l],
			    xgrid[l], ygrid[l], zgrid[l]);
	  }
     	}

    if (!both) continue;
    iso = -iso;
    
    for (k = 0; k < kdim[l] - 1; k++)
      for (j = 0; j < jdim[l] - 1; j++)
	for (i = 0 ; i < idim[l] - 1; i++) {
	  last_poly [(k*jdim[l]+j)*idim[l]+i] = 0;
	  first_poly[(k*jdim[l]+j)*idim[l]+i] = NUM_POLYGONS + 1;
	  calc_index_and_temps (data[l], i, j, k,
				idim[l], jdim[l], kdim[l],
				iso, &index);
	  if (index) {
	    get_cell_verts (index, i, j, k, iso, crossings);
	    get_cell_polys (index, crossings, i, j, k,
			    idim[l],  jdim[l],  kdim[l],
			    xgrid[l], ygrid[l], zgrid[l]);
	  }
	}

    free (first_poly);
    free (last_poly);
  }

  if (!NUM_POLYGONS)
    message (routine, "no polygons extracted --- bad iso-value", WARNING);
  else if (too_many_polys)
    message (routine, "memory allocation failure, no surface made", WARNING);
  else {
    S = (Iso*) malloc (sizeof (Iso));
    S -> nvert = NUM_VERTICES;
    S -> npoly = NUM_POLYGONS;
    S -> pxyz  = (float*) malloc (3 * NUM_VERTICES * sizeof (float));
    S -> nxyz  = (float*) calloc (3 * NUM_VERTICES,  sizeof (float));
    S -> plist = (int*)   malloc (3 * NUM_POLYGONS * sizeof (int));

    for (i = 0; i < NUM_VERTICES; i++) {
      S -> pxyz[3 * i    ] = XVERTICES[i];
      S -> pxyz[3 * i + 1] = YVERTICES[i];
      S -> pxyz[3 * i + 2] = ZVERTICES[i];
    }

    j = 3 * NUM_POLYGONS;
    for (i = 0; i < j; i++)
      S -> plist[i] = POLYGONS[i] - 1;

    j = NUM_POLYGONS;
    for (nzero = 0, i = 0; i < j; i++)
      if (calc_norm (S -> pxyz, S -> plist + 3 * i, S -> nxyz) == 0.0) nzero++;

    j = NUM_VERTICES;
    for (i = 0; i < j; i++) {
      vect = S -> nxyz + 3 * i;
      mag  = v_mag (vect);
      mag  = (mag == 0.0) ? mag = 1.0 : 1.0 / mag;
      vect[0] *= mag;
      vect[1] *= mag;
      vect[2] *= mag;
    }

    cout << "Isosurface has "
	 << NUM_VERTICES << " vertices, "
	 << NUM_POLYGONS << " triangles, "
	 << nzero        << " zero length normals" << endl;
  }
  
  free (POLYGONS);
  free (ZVERTICES);
  free (YVERTICES);
  free (XVERTICES);

  return S;
}


void flipNorms (Iso* S)
/* ------------------------------------------------------------------------- *
 * Negate directions of all surface normals for S.
 * ------------------------------------------------------------------------- */
{
  register int    i;
  const int       ntot      = 3 * S -> nvert;
  register float* component = S -> nxyz;

  for (i = 0; i < ntot; i++) component[i] *= -1.0F;
}


static int add_polygon (float* p[3],
			int    i   ,
			int    j   ,
			int    k   ,
			int    idim,
			int    jdim,
			int    kdim,
			float* gx  ,
			float* gy  ,
			float* gz  )
/* ------------------------------------------------------------------------- *
 * This subroutine stores a polygon (triangle) in a list of vertices
 * and POLYGONS.
 * ------------------------------------------------------------------------- */
{
  char routine[] = "add_polygon";
  int   l,m,n,cell,ioff,joff,koff,offset;
  int   old_vertex;
  float point[3];

  if (NUM_VERTICES >= (VERT_MAX-3)) {
    message (routine, "too many vertices, aborting", WARNING);
    too_many_polys = 1;
    return 0;
  }

  if (NUM_POLYGONS>= (POLY_MAX-1)) {
    message (routine, "too many polygons, aborting", WARNING);
    too_many_polys = 1;
    return 0;
  }

  /* -- Calculate the grid cell offset for cell (i,j,k). */

  offset = (k*jdim+j)*idim+i;

  /* -- For each of the 3 vertices in the new polygon: */

  for (n = 0; n < 3; n++) {
    
    /* We must first convert the integral-grid vertex values to
     * physical values based on grid coordinates - messy tri-linear
     * interpolation NB: a rather approximate estimate of the physical
     * vertex co-ords.  */
    
    point[0] = tri_lin_int (p[n],i,j,k,gx+offset);
    point[1] = tri_lin_int (p[n],i,j,k,gy+offset);
    point[2] = tri_lin_int (p[n],i,j,k,gz+offset);
    
    /* Next: check 9 cells in the previous (k-1) layer, 3 cells in the
     * previous (j-1) row, and the previous (i-1) cell and this cell
     * (i) to see if the 3 new vertices correspond to previous
     * vertices --- if they do, then we don't want to create new
     * vertices, just reference the old ones, for each of the thirteen
     * possible cells that the vertex could be joined to.  */

    old_vertex = 0;

    for (l = 0; !old_vertex && l >= -13 ; l--)    {
      koff = l / 9;
      joff = l / 3 - 3 * koff;
      ioff = l - 3 * (joff + 3 * koff);
      if (k+koff >= 0 && j+joff >= 0 && i+ioff >= 0) {
        cell = ((k+koff)*jdim + (j+joff))*idim + i+ioff;

	/* -- For each of the polygons in cell: */

        for (m = first_poly[cell]; m <= last_poly[cell]; m++)
          if (old_vertex = add_vertex (n, point, m)) break;
      }
    }

    /* -- If old_vertex == 0, we have a new vertex - add it. */

    if (!old_vertex) {
      XVERTICES[NUM_VERTICES] = point[0];
      YVERTICES[NUM_VERTICES] = point[1];
      ZVERTICES[NUM_VERTICES] = point[2];
      NUM_VERTICES++;
      POLYGONS[NUM_POLYGONS*3+n] = NUM_VERTICES;
    }
  }

  NUM_POLYGONS++;
  last_poly[(k*jdim+j)*idim+i] = NUM_POLYGONS;
        
  return 1;
}


static int add_vertex (int    n     ,
		       float* vertex,
		       int    poly  )
/* ------------------------------------------------------------------------- *
 * poly = the number of the polygon we are checking - the numbers
 * of its three vertices are stored in v  N.B. the polygon numbers
 * start at 1, whereas the storage in POLYGONS starts at zero
 * hence, the indexing below.
 * ------------------------------------------------------------------------- */
{
  int i, v[3];

  v[0] = POLYGONS[3*poly-3];
  v[1] = POLYGONS[3*poly-2];
  v[2] = POLYGONS[3*poly-1];

  /* -- Check against the 3 vertices of polygon poly. */

  for (i = 0; i < 3 ; i++)
    if (equal_vertex (vertex, v[i])) {
      POLYGONS[NUM_POLYGONS*3+n] = v[i];
      return 1;
    }

  return 0;
}


static int equal_vertex (float* new_vertex,
			 int    old_vertex)
/* ------------------------------------------------------------------------- *
 *  This routine checks to see if the vertex referenced by pointer
 *  new_vertex is equal to XVERTICES[old_vertex-1],
 *  YVERTICES[old_vertex-1], Z....  If equal, return 1.
 *
 *  4 decimal accuracy is considered to be sufficient.
 * ------------------------------------------------------------------------- */
{
  if ((int)(10000*new_vertex[0]) == (int)(10000*XVERTICES[old_vertex-1]) &&
      (int)(10000*new_vertex[1]) == (int)(10000*YVERTICES[old_vertex-1]) &&
      (int)(10000*new_vertex[2]) == (int)(10000*ZVERTICES[old_vertex-1])  )
    return 1;
  else
    return 0;
}


static float tri_lin_int (float  p[3],
			  int    i   ,
			  int    j   ,
			  int    k   ,
			  float* tg  )
/* ------------------------------------------------------------------------- *
 * Tri-linear interpolation.  Inline this?
 * ------------------------------------------------------------------------- */
{
  float x, y, z, value;

  x = p[0];
  y = p[1];
  z = p[2];

  value =
    ((float)(i+1)-x)*((float)(j+1)-y)*((float)(k+1)-z) * tg[goff[0][0][0]] +
    (  x-(float)(i))*((float)(j+1)-y)*((float)(k+1)-z) * tg[goff[1][0][0]] +
    ((float)(i+1)-x)*(  y-(float)(j))*((float)(k+1)-z) * tg[goff[0][1][0]] +
    (  x-(float)(i))*(  y-(float)(j))*((float)(k+1)-z) * tg[goff[1][1][0]] +
    ((float)(i+1)-x)*((float)(j+1)-y)*(  z-(float)(k)) * tg[goff[0][0][1]] +
    (  x-(float)(i))*((float)(j+1)-y)*(  z-(float)(k)) * tg[goff[1][0][1]] +
    ((float)(i+1)-x)*(  y-(float)(j))*(  z-(float)(k)) * tg[goff[0][1][1]] +
    (  x-(float)(i))*(  y-(float)(j))*(  z-(float)(k)) * tg[goff[1][1][1]] ;
  
  return value;
}

/*****************************************************************************
 * Remainder of this file contains the marching cubes algorithm used to
 * generate wireframe meshes.
 *
 * The marching cubes surface tiler is described by
 * Lorensen & Cline in the SIGGRAPH 87 Conference Proceedings.
 * 
 * The wire frame generator part this code was developed at the
 * National Center for Supercomputing Applications at the University
 * of Illinois at Urbana-Champaign.
 * 
 * THE UNIVERSITY OF ILLINOIS GIVE NO WARRANTY, EXPRESSED OR IMPLIED
 * FOR THE SOFTWARE AND/OR DOCUMENTATION PROVIDED, INCLUDING, WITHOUT
 * LIMITATION, WARRANTY OF MERCHANTABILITY AND WARRANTY OF FITNESS FOR
 * A PARTICULAR PURPOSE
 *
 * Written by Mike Krogh, NCSA, Feb. 2, 1990.
 *****************************************************************************/

#include "cell_table.h"

float DATA1,DATA2,DATA3,DATA4,DATA5,DATA6,DATA7,DATA8;


static void calc_index_and_temps (register float* data ,
				  int             x1   ,
				  int             y1   ,
				  int             z1   ,
				  int             xdim ,
				  int             ydim ,
				  int             zdim ,
				  register float  diso ,
				  int*            index)
/* ------------------------------------------------------------------------- *
 * This subroutine calculates the index and creates some global
 * temporary variables (for speed).
 * ------------------------------------------------------------------------- */
{
  register float *tmp;

  *index = 0;

  tmp = data + (z1*xdim*ydim) + (y1*xdim) + x1;

  *index += (diso <= (DATA1 = *(tmp)));
  *index += (diso <= (DATA2 = *(tmp + 1))) * 2;

  tmp += xdim;
  *index += (diso <= (DATA3 = *(tmp + 1))) * 4;
  *index += (diso <= (DATA4 = *(tmp))) * 8;

  tmp = tmp - xdim + xdim*ydim;
  *index += (diso <= (DATA5 = *(tmp))) * 16;
  *index += (diso <= (DATA6 = *(tmp + 1))) * 32;

  tmp += xdim;
  *index += (diso <= (DATA7 = *(tmp + 1))) * 64;
  *index += (diso <= (DATA8 = *(tmp))) * 128;
}


static void get_cell_verts (int   index           ,
			    int   x1              ,
			    int   y1              ,
			    int   z1              ,
			    float diso            ,
			    float crossings[13][3])
/* ------------------------------------------------------------------------- *
 * Load values of crossings.
 * ------------------------------------------------------------------------- */
{
  register int i;
  register int x2,y2,z2;
  int          nedges;
  int          crnt_edge;
 
#define linterp(a1,a2,a,b1,b2) \
               ((float)(((a-a1)*(float)(b2-b1)/(a2-a1))+(float)b1))

  x2 = x1+1;
  y2 = y1+1;
  z2 = z1+1;

  nedges = cell_table[index].nedges;
  for (i = 0; i < nedges; i++) {
     crnt_edge = cell_table[index].edges[i];
     switch (crnt_edge) {
     case 1:
       crossings[1][0] = linterp(DATA1,DATA2,diso,x1,x2);
       crossings[1][1] = (float)y1;
       crossings[1][2] = (float)z1;
       break;

     case 2:
       crossings[2][1] = linterp(DATA2,DATA3,diso,y1,y2);
       crossings[2][0] = (float)x2;
       crossings[2][2] = (float)z1;
       break;

     case 3:
       crossings[3][0] = linterp(DATA4,DATA3,diso,x1,x2);
       crossings[3][1] = (float)y2;
       crossings[3][2] = (float)z1;
       break;

     case 4:
       crossings[4][1] = linterp(DATA1,DATA4,diso,y1,y2);
       crossings[4][0] = (float)x1;
       crossings[4][2] = (float)z1;
       break;

     case 5:
       crossings[5][0] = linterp(DATA5,DATA6,diso,x1,x2);
       crossings[5][1] = (float)y1;
       crossings[5][2] = (float)z2;
       break;

     case 6:
       crossings[6][1] = linterp(DATA6,DATA7,diso,y1,y2);
       crossings[6][0] = (float)x2;
       crossings[6][2] = (float)z2;
       break;

     case 7:
       crossings[7][0] = linterp(DATA8,DATA7,diso,x1,x2);
       crossings[7][1] = (float)y2;
       crossings[7][2] = (float)z2;
       break;

     case 8:
       crossings[8][1] = linterp(DATA5,DATA8,diso,y1,y2);
       crossings[8][0] = (float)x1;
       crossings[8][2] = (float)z2;
       break;

     case 9:
       crossings[9][2] = linterp(DATA1,DATA5,diso,z1,z2);
       crossings[9][1] = (float)y1;
       crossings[9][0] = (float)x1;
       break;

     case 10:
       crossings[10][2] = linterp(DATA2,DATA6,diso,z1,z2);
       crossings[10][1] = (float)y1;
       crossings[10][0] = (float)x2;
       break;
       
     case 11:
       crossings[11][2] = linterp(DATA4,DATA8,diso,z1,z2);
       crossings[11][1] = (float)y2;
       crossings[11][0] = (float)x1;
       break;

     case 12:
       crossings[12][2] = linterp(DATA3,DATA7,diso,z1,z2);
       crossings[12][1] = (float)y2;
       crossings[12][0] = (float)x2;
       break;

     }
  }
}


static void get_cell_polys (int    index           ,
			    float  crossings[13][3],
			    int    i               ,
			    int    j               ,
			    int    k               ,
			    int    idim            ,
			    int    jdim            ,
			    int    kdim            ,
			    float* gx              ,
			    float* gy              ,
			    float* gz              )
/* ------------------------------------------------------------------------- *
 * This subroutine will calculate the polygons.
 * ------------------------------------------------------------------------- */
{
  register int num_o_polys;
  register int poly;
  float        *p[3];
  float        n1[3],n2[3],n3[3];
  
  num_o_polys = cell_table[index].npolys;
  
  for (poly = 0; poly < num_o_polys; poly++) {

    p[0] = &crossings[cell_table[index].polys[poly*3]  ][0];
    p[1] = &crossings[cell_table[index].polys[poly*3+1]][0];
    p[2] = &crossings[cell_table[index].polys[poly*3+2]][0];
    
    /* store the polygons */
    
    if (!(add_polygon (p, i, j, k, idim, jdim, kdim, gx, gy, gz))) {
      printf ("\n    *** error from add_polygon\n");
      exit   (1);
    }
  }
}


static float calc_norm (float* poly,
			int*   vert,
			float* norm)
/* ------------------------------------------------------------------------- *
 * For each triangle, calculate a unit surface normal that points away
 * from the origin using any two of the sides - ADD this normal to any
 * unit normals that have been previously calculated for EACH of the 3
 * vertices in the triangle (if the vertex is AT the origin or if its
 * normal is perpendicular to the radius vector, ensure the normal has
 * positive x co-ord (or positive y-co-ord if x-co-ord=0 (or ...))
 * ------------------------------------------------------------------------- */
{
  float r1[3],r2[3],r3[3], mag;
  
  r1[0] = *(poly+3*vert[0])   - *(poly+3*vert[2]);
  r2[0] = *(poly+3*vert[1])   - *(poly+3*vert[2]);
  r1[1] = *(poly+3*vert[0]+1) - *(poly+3*vert[2]+1);
  r2[1] = *(poly+3*vert[1]+1) - *(poly+3*vert[2]+1);
  r1[2] = *(poly+3*vert[0]+2) - *(poly+3*vert[2]+2);
  r2[2] = *(poly+3*vert[1]+2) - *(poly+3*vert[2]+2);
  
  r3[0] = r1[1]*r2[2] - r2[1]*r1[2];
  r3[1] = r2[0]*r1[2] - r1[0]*r2[2];
  r3[2] = r1[0]*r2[1] - r2[0]*r1[1];
  
  if ((mag = v_mag(r3)) == 0.0) return mag;

  /* -- Add the new normal to any existing normals at each triangle vertex. */

  *(norm+3*vert[0])   += r3[0];
  *(norm+3*vert[0]+1) += r3[1];
  *(norm+3*vert[0]+2) += r3[2];
  
  *(norm+3*vert[1])   += r3[0];
  *(norm+3*vert[1]+1) += r3[1];
  *(norm+3*vert[1]+2) += r3[2];
  
  *(norm+3*vert[2])   += r3[0];
  *(norm+3*vert[2]+1) += r3[1];
  *(norm+3*vert[2]+2) += r3[2];
  
  return mag;
}
