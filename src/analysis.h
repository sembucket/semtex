#ifndef ANALYSIS_H
#define ANALYSIS_H


class Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow
// solver.  This is designed to be overridden at implementation level
// if needed.
// ===========================================================================
{
public:
  Analyser  (Domain*, FEML*);
  ~Analyser () { }

  void analyse (AuxField**);

protected:
  Domain*               src      ; // Source information.
  ofstream              par_strm ; // File for particle tracking.
  ofstream              his_strm ; // File for history points.
  ofstream              mdl_strm ; // File for modal energies.
  vector<HistoryPoint*> history  ; // Locations, etc. of history points.
  list<FluidParticle*>  particle ; // List of fluid particles.
  vector<Point*>        initial  ; // Starting locations of particles.
  Statistics*           stats    ; // Field average statistics.

  void modalEnergy ();
  void divergence  (AuxField**) const;
  void estimateCFL ()           const;
};

#endif
