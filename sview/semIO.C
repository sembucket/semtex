///////////////////////////////////////////////////////////////////////////////
// Input SEM mesh and field data.
// 
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sview.h>
#include <limits.h>


Sem* loadMesh (ifstream& mfile)
// ---------------------------------------------------------------------------
// Mesh header is of form:
//
// 9 9 32 206 NR NS NZ NEL
//
// then follows (x,y) locations of each element's 2D mesh points in
// row-major order, ASCII format.  In all local implementations,
// NR=NS=np.  Finally the NZ + 1 z locations of the 3D mesh are
// supplied.
//
// NB: the mean values of extent in the x, y, & z directions are subtracted
// from the mesh locations.
// ---------------------------------------------------------------------------
{
  char     routine[] = "loadMesh", buf[StrMax], err[StrMax];
  register int  i, n;
  register double xavg, yavg, zavg;
  Sem*     M = new Sem;

  mfile >> n >> M -> np >> M -> nz >> M -> nel;
  mfile.getline (buf, StrMax);

  if (!strstr (buf, "NR NS NZ NEL")) {
    sprintf (err, "mesh header line should include NR NS NZ NEL: %s", buf);
    message (routine, err, ERROR);
  }
  if (n != M -> np) {
    sprintf (err, "Element NR, NS orders must be equal: %1d: %1d", n, M -> np);
    message (routine, err, ERROR);
  }
  if (M -> nz < 2) {
    sprintf (err, "Mesh mush be 3D (NZ >= 2): %1d", M -> nz);
    message (routine, err, ERROR);
  }
  
  M -> ntot = M -> np * M -> np * M -> nel;

  M -> xmesh = new float [M -> ntot];
  M -> ymesh = new float [M -> ntot];
  M -> zmesh = new float [M -> nz + 1];

  State.xmin = 1.0e6;
  State.xmax = -State.xmin;
  State.ymin =  State.xmin;
  State.ymax =  State.xmax;
  State.zmin =  State.zmin;
  State.zmax =  State.zmax;

  n = M -> ntot;
  for (i = 0; i < n; i++) {
    mfile >> M -> xmesh[i] >> M -> ymesh[i];
    State.xmin = min (State.xmin, M -> xmesh[i]);
    State.ymin = min (State.ymin, M -> ymesh[i]);
    State.xmax = max (State.xmax, M -> xmesh[i]);
    State.ymax = max (State.ymax, M -> ymesh[i]);
  }
  xavg = 0.5 * (State.xmin + State.xmax);
  yavg = 0.5 * (State.ymin + State.ymax);
  for (i = 0; i < n; i++) {
    M -> xmesh[i] -= xavg;
    M -> ymesh[i] -= yavg;
  }

  n = M -> nz + 1;
  for (i = 0; i < n; i++) {
    mfile >> M -> zmesh[i];
    State.zmin = min (State.zmin, M -> zmesh[i]);
    State.zmax = max (State.zmax, M -> zmesh[i]);
  }
  zavg = 0.5 * (State.zmin + State.zmax);
  for (i = 0; i < n; i++)
    M -> zmesh[i] -= zavg;

  State.xmin -= xavg;
  State.xmax -= xavg;
  State.ymax -= yavg;
  State.ymax -= yavg;
  State.zmin -= zavg;
  State.zmax -= zavg;

  State.length = max (hypot (xavg, yavg), hypot (xavg, zavg));
  State.length = max (State.length,       hypot (yavg, zavg));
  State.length *= 2.0;

  return M;
}
