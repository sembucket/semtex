/*****************************************************************************
 * MESH.C:  Routines for mesh I/O, geometric factors.                        *
 *****************************************************************************/

static char
RCSid[] = "$Id$";

#include "Fem.h"

#define  readLine(fp)  fgets(buf, STR_MAX, (fp))
#define  echo          if (verb) fputs(buf, stdout)
#define  skipComments  if (*buf != '#' && *buf != '\n') break
#define  upperCase(s)  { char *z=(s); while (*z=toupper(*z)) z++; }





Element *readMesh(FILE  *fp)
/* ========================================================================= *
 * Read mesh description from fp, return a pointer to a list of elements     *
 * with unique storage locations for element data structure allocated (this  *
 * does not include geometric partials and the like, which can be shared).   *
 *                                                                           *
 * Mesh node locations (xmesh & ymesh) are computed.                         *
 *                                                                           *
 * NB: SIDE EFFECTS: iparam "NEL" is set here.                               *
 *                                                                           *
 * MESH FILE DESCRIPTION:                                                    *
 * ---------------------                                                     *
 * The first line of the file contains the number of element corner          *
 * vertices: <num> VERTICES.  Following this is a list of corner vertices    *
 * given in Cartesian form, one vertex per line.                             *
 *                                                                           *
 * The list of vertices is followed by a blank line, then a section which    *
 * describes the elements, starting with the line <num> ELEMENTS.            *
 * This is also followed by a blank line and domain information is then      *
 * given on an element-wise basis.  No global node-numbering is needed; this *
 * is generated automatically (& dealt with elsewhere).                      *
 *                                                                           *
 * Corner vertex numbering proceeds CCW round the element.  Side 1 lies      *
 * between corner vertices 1 & 2 and side numbering is also CCW-sequential   *
 * around the element.  Sides can either lie on the domain boundary (BC) or  *
 * mate with another element (EL).  Internal numbering of elements and sides *
 * starts at zero.                                                           *
 *                                                                           *
 * Information for each element is followed by a blank line, except that the *
 * last element may be followed by EOF.                                      *
 *                                                                           *
 * Element data for a biquadratic element:                                   *
 * ELEMENT number ORDER 2                                                    *
 * node1 node2 node3 node4       // these are indices in the vertex list     *
 * SIDE 1  BC  tag-num                                                       *
 * SIDE 2  EL  mate-el SIDE mate-side                                        *
 * SIDE 3  EL  mate-el SIDE mate-side                                        *
 * SIDE 4  EL  mate-el SIDE mate-side                                        *
 *                                                                           *
 * ========================================================================= */
{
  char       routine[] = "readMesh", s[STR_MAX];
  char       buf[STR_MAX];
  int        i, elmt;

  int        verb = iparam("VERBOSE");
  int        nvert, side, np, nq, nel, ntot, nedg;

  Point     *vertex;
  Element   *E;
  int      **elmtVert;


  /* -- Input vertex information. */

  while (readLine (fp)) {
    echo; upperCase (buf);
    if ( (strstr (buf, "VERTICES")) && (sscanf (buf, "%d", &nvert))) {
      readLine (fp);
      break;
    } else {
      sprintf (s, "couldn't get number of vertices: %s", buf);
      message (routine, s, ERROR);
    }
  }

  vertex = (Point *) malloc (nvert*sizeof(Point));

  i = 0;
  while (readLine (fp) && i < nvert) {
    echo;
    if ( sscanf (buf, "%*s%lf%lf", &vertex[i].x, &vertex[i].y) != 2) {
      sprintf (s, "expected info for vertex %1d, got: %s", i+1, buf);
      message (routine, s, ERROR);
    }
    i++;
  }
  
  /* -- Check that vertex input completed OK. */

  if (i < nvert) {
    sprintf (s, "expected another vertex, got newline at vertex #%1d", i);
    message (routine, s, ERROR);
  } else if (*buf != '\n') {
    sprintf (s, "read all vertices, but next line not blank: %s", buf);
    message (routine, s, ERROR);
  }


  /* -- Input element information. */

  while (readLine(fp)) {
   echo; upperCase(buf);
    if ( (strstr(buf, "ELEMENTS")) && (sscanf(buf, "%d", &nel)))
      break;
    else{
      sprintf (s, "couldn't get number of elements: %s", buf);
      message (routine, s, ERROR);
    }
  }

  elmtVert = imatrix (0, nel-1, 0, 3);
  
  E = (Element *) calloc (nel, sizeof (Element));
  for (i = 0; i < nel - 1; i++) (E + i) -> next = E + i + 1;

  for (i = 0; i < nel; i++) {
    readLine (fp);
    readLine (fp);
    echo; upperCase (buf);

    if (strstr (buf, "ELEMENT")) {
      char *sp;
	    
      sscanf (buf, "%*s %d %*s %d", &elmt, &np);
      if (elmt > nel)
	message (routine, "element id tag exceeds number of elements", ERROR);
      elmt--;
      (E + elmt) -> id      = elmt;
      (E + elmt) -> np      = ++np;
      (E + elmt) -> fldtype = VELOCITY;
      (E + elmt) -> name    = 'u';

      readLine (fp); echo;
      sscanf (sp = strtok (buf, " \t"), "%d", elmtVert[elmt]);
      elmtVert[elmt][0]--;
      for (side = 1; side < 4; side++) {
	sscanf (sp = strtok (NULL, " \t\0"), "%d", elmtVert[elmt] + side);
	elmtVert[elmt][side]--;
      }
	   
      (E + elmt) -> mate = (Link *) calloc (4, sizeof (Link));

      for (side = 0; side < 4; side++) {
	readLine (fp); echo; upperCase (buf);

	if (strstr (buf, "BC")) {

	  sscanf (buf, "%*s %*s %*s %d", &(E + elmt) -> mate[side].facetag);
	  (E + elmt) -> mate[side].elmt = DOMAIN_BOUNDARY;
	  (E + elmt) -> mate[side].facetag--;
	 
	} else if (strstr (buf, "EL")) {

	  sscanf (buf, "%*s %*s %*s %d %*s %d",
		  &(E + elmt) -> mate[side].elmt,
		  &(E + elmt) -> mate[side].facetag);
	  (E + elmt) -> mate[side].elmt--;
	  (E + elmt) -> mate[side].facetag--;


	} else {
	  sprintf (s, "can't parse element %1d edge %1d data: %s",
		   elmt + 1, side + 1, buf);
	  message (routine, s, ERROR);
	}
      }
    } else {
      sprintf (s, "element descriptor for #%1d? : %s", i+1, buf);
      message (routine, s, ERROR);
    }
  }

  /* -- Check element input and enforce equal orders in all elements. */

  np = E -> np;
  nq = (option ("RULE") == GL) ? quadComplete (DIM, np) : np;

  for (elmt = 0; elmt < nel; elmt++) {
    if ( (E + elmt) -> id != elmt ) {
      sprintf (s, "element %1d not input", elmt+1);
      message (routine, s, ERROR);
    }
    (E + elmt) -> np = np;
    (E + elmt) -> nq = nq;
  }

  /* -- Check nominated connectivity. */

  {
    int localElmt, remoteElmt, localSide, remoteSide;


    for (elmt = 0; elmt < nel; elmt++)
      for (side = 0; side < 4; side++) {
	remoteElmt = (E + elmt) -> mate[side].elmt;
	remoteSide = (E + elmt) -> mate[side].facetag;
	  if (remoteElmt != DOMAIN_BOUNDARY) {
	    localElmt = (E + remoteElmt) -> mate[remoteSide].elmt;
	    localSide = (E + remoteElmt) -> mate[remoteSide].facetag;
	    if (localElmt != elmt || localSide != side) {
	      sprintf (buf, "connectivity problem:"
		       " %1d.%1d -> %1d.%1d -> %1d.%1d",
		       elmt       + 1, side       + 1,
		       remoteElmt + 1, remoteSide + 1,
		       localElmt  + 1, localSide  + 1);
	      message (routine, buf, ERROR);
	    }
	  }
      }
  }

  /* -- All the mesh input file has been read.  Allocate storage.   *
   *      ntot: total number of storage locations;                  *
   *      nedg: total number of element-edge storage locations;     *
   *      np  : number of points on an edge (polynomial order + 1). *
    --------------------------------------------------------------- */

  setIparam ("NEL", nel);
  ntot = nel * SQR (np);
  nedg = nel * 4 * (np - 1);
  
  E -> value    = (double **) malloc (np * nel * sizeof (double*));
  E -> xmesh    = (double **) malloc (np * nel * sizeof (double*));
  E -> ymesh    = (double **) malloc (np * nel * sizeof (double*));

  E -> value[0] = dvector (0, ntot-1);
  E -> xmesh[0] = dvector (0, ntot-1);
  E -> ymesh[0] = dvector (0, ntot-1);

  E -> solve    = ivector (0, nedg-1);
  E -> bmap     = ivector (0, nedg-1);
  E -> emap     = getEmap (np);

  if (!(   E -> value[0] && E -> xmesh[0] && E -> ymesh[0]
	&& E -> value    && E -> xmesh    && E -> ymesh
	&& E -> solve    && E -> bmap     && E -> emap    ))
    message (routine, "mesh storage allocation failed", ERROR);

  /* -- Set pointers for internal arrays. This could do with a clean. */

  for (elmt = 1; elmt < nel; elmt++) {
    (E + elmt) -> value    = (E + elmt - 1) -> value    + np;
    (E + elmt) -> xmesh    = (E + elmt - 1) -> xmesh    + np;
    (E + elmt) -> ymesh    = (E + elmt - 1) -> ymesh    + np;
    (E + elmt) -> value[0] = (E + elmt - 1) -> value[0] + SQR (np);
    (E + elmt) -> xmesh[0] = (E + elmt - 1) -> xmesh[0] + SQR (np);
    (E + elmt) -> ymesh[0] = (E + elmt - 1) -> ymesh[0] + SQR (np);
    (E + elmt) -> solve    = (E + elmt - 1) -> solve    + 4*(np - 1);
    (E + elmt) -> bmap     = (E + elmt - 1) -> bmap     + 4*(np - 1);
    (E + elmt) -> emap     = E -> emap;
  }
 
  for (elmt = 0; elmt < nel; elmt++)
    for (i = 1; i < np; i++) {
      (E + elmt) -> value[i] = (E + elmt) -> value[i - 1] + np;
      (E + elmt) -> xmesh[i] = (E + elmt) -> xmesh[i - 1] + np;
      (E + elmt) -> ymesh[i] = (E + elmt) -> ymesh[i - 1] + np;
    }
 
  /* -- Fill mesh internal node locations. */

  for (elmt = 0; elmt < nel; elmt++)
    meshElement (E + elmt, elmtVert[elmt], vertex);

  freeImatrix (elmtVert, 0, 0);
  free (vertex);

  return E;
}





void  printMesh (FILE      *fp,
		 Element   *E )
/* ========================================================================= *
 * Information is written out element-by-element.                            *
 * ========================================================================= */
{
  int      nel  = iparam ("NEL");
  int      np   = E -> np;
  int      ntot = nel * SQR (np);
  double  *x    = *E -> xmesh,
          *y    = *E -> ymesh;

  
  fprintf(fp, "%1d %1d %1d %1d NR NS NZ NEL\n", np, np, 1, nel);

  printDvector(fp, 15, 6, ntot, 1, 2, x, y);
}





int  quadComplete(int dim ,	              /* Number of space dimensions  */
		  int np  )	              /* Number of points on an edge */
/* ========================================================================= *
 * Return the number of Gauss-Legendre quadrature points sufficient to       *
 * achieve the full rate of convergence for tensor-product element bases.    *
 *                                                                           *
 * References: Hughes \S 4.1, Strang & Fix \S 4.3.                           *
 * ========================================================================= */
{
  int  n, ktot;
  

  ktot = (dim + 1)*(np - 1) - 2;
  n = (ktot & 0) ? ktot + 2 : ktot + 1;
  n >>= 1;

  return MAX (n, 2);
}
