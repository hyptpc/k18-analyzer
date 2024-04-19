//  Authors: Wooseung Jung

#ifndef HYPTPCSPACEPOINTMEASUREMENT_HH
#define HYPTPCSPACEPOINTMEASUREMENT_HH

//GenFit
#include <SpacepointMeasurement.h>
#include <TrackCandHit.h>

//k18-analyzer
#include "TPCLTrackHit.hh"
#include "HypTPCHit.hh"

namespace genfit {

  class HypTPCSpacepointMeasurement : public SpacepointMeasurement{

  public:
    HypTPCSpacepointMeasurement() : SpacepointMeasurement() {}
    HypTPCSpacepointMeasurement(const HypTPCHit* dethit, const TrackCandHit* hit); //TPCLTrackHit Hitpos & Resolution
    ~HypTPCSpacepointMeasurement(){}

    virtual AbsMeasurement* clone() const { return new HypTPCSpacepointMeasurement(*this); }
    double GetCharge() { return fCharge; }

  private:
    double fCharge; //currently not used

    ClassDef(HypTPCSpacepointMeasurement,1)
  };

} /* End of namespace genfit */
#endif // HypTPCSpacepointMeasurement_hh
