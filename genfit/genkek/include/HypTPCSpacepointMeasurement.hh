//  Authors: Wooseung Jung

#ifndef HYPTPCSPACEPOINTMEASUREMENT_HH
#define HYPTPCSPACEPOINTMEASUREMENT_HH

//GenFit2
#include <SpacepointMeasurement.h>
#include <TrackCandHit.h>

//k18-analyzer
#include "TPCLTrackHit.hh"

//class TVector3;

class HypTPCSpacepointMeasurement : public genfit::SpacepointMeasurement{

public:
  HypTPCSpacepointMeasurement() : genfit::SpacepointMeasurement() {}
  HypTPCSpacepointMeasurement(const TPCLTrackHit* dethit, const genfit::TrackCandHit* hit); //TPCHit Hitpos & Resolution

  virtual HypTPCSpacepointMeasurement* clone() const { return new HypTPCSpacepointMeasurement(*this); }
  double GetCharge() { return fCharge; }

private:
  double fCharge;

  ClassDef(HypTPCSpacepointMeasurement,1)
};

#endif // HypTPCSpacepointMeasurement_hh
