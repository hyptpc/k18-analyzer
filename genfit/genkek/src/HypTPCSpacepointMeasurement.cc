//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCSpacepointMeasurement.hh"

//ROOT
#include <TMatrixDSymfwd.h>                // for TMatrixDSym
#include <TVector3.h>

ClassImp(genfit::HypTPCSpacepointMeasurement)

namespace genfit{
  HypTPCSpacepointMeasurement::HypTPCSpacepointMeasurement(
const HypTPCHit* dethit, const TrackCandHit* hit)
    : SpacepointMeasurement()
  {
    const TPCLTrackHit& tpchit = dethit -> GetHit();
    const TVector3& res_vect = tpchit.GetResolutionVect();

    //GenFit Units : GeV/c, ns, cm, kGauss
    //K1.8Ana Units : GeV/c, ns, mm, T
    int nDim = 3;
    TMatrixDSym hitCov(nDim);
    hitCov.Zero();

    //Transverse & Vertical position resolution
    double resX = 0.1*res_vect.X();
    double resY = 0.1*res_vect.Y();
    double resZ = 0.1*res_vect.Z();
    hitCov(0, 0) = resX*resX;
    hitCov(1, 1) = resY*resY;
    hitCov(2, 2) = resZ*resZ;

    const TVector3& pos = tpchit.GetLocalHitPos();
    rawHitCoords_(0) = pos.X()/10.; //mm -> cm
    rawHitCoords_(1) = pos.Y()/10.;
    rawHitCoords_(2) = pos.Z()/10.;
    rawHitCov_ = hitCov;
    detId_ = hit -> getDetId();
    hitId_ = hit -> getHitId();
    this -> initG();
  }

} /* End of namespace genfit */
