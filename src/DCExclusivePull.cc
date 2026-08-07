// -*- C++ -*-

#include "DCExclusivePull.hh"

#include <vector>

#include <TMath.h>

#include "DCGeomMan.hh"
#include "DCLTrackHit.hh"
#include "DCLocalTrack.hh"
#include "MathTools.hh"

namespace DCExclusivePull
{
namespace
{
const auto& gGeom = DCGeomMan::GetInstance();
constexpr Int_t MinHitsExclusive = 5; // leave-one-out needs >=4 hits remaining
}

//_____________________________________________________________________________
Double_t
Residual(const DCLocalTrack& track, Int_t ihit)
{
  if(!track.IsFitted()) return TMath::QuietNaN();

  const auto& hits = track.GetHitArray();
  const Int_t n = static_cast<Int_t>(hits.size());
  if(ihit < 0 || ihit >= n || n < MinHitsExclusive)
    return TMath::QuietNaN();

  const Int_t n_use = n - 1;
  std::vector<Double_t> z(n_use), w(n_use), s(n_use), ct(n_use), st(n_use);
  Int_t k = 0;
  for(Int_t i = 0; i < n; ++i){
    if(i == ihit) continue;
    const DCLTrackHit* hitp = hits[i];
    if(!hitp) return TMath::QuietNaN();
    const Int_t lnum = hitp->GetLayer();
    const Double_t ww = gGeom.GetResolution(lnum);
    if(ww <= 0.) return TMath::QuietNaN();
    const Double_t aa = hitp->GetTiltAngle() * TMath::DegToRad();
    z[k]  = hitp->GetZ();
    w[k]  = 1. / (ww * ww);
    s[k]  = hitp->GetLocalHitPos();
    ct[k] = TMath::Cos(aa);
    st[k] = TMath::Sin(aa);
    ++k;
  }
  if(k != n_use) return TMath::QuietNaN();

  Double_t x0 = 0., u0 = 0., y0 = 0., v0 = 0.;
  if(!MathTools::SolveGaussJordan(z, w, s, ct, st, x0, u0, y0, v0))
    return TMath::QuietNaN();

  const DCLTrackHit* hit = hits[ihit];
  if(!hit) return TMath::QuietNaN();
  const Double_t aa = hit->GetTiltAngle() * TMath::DegToRad();
  const Double_t zz = hit->GetZ();
  const Double_t scal = (x0 + u0 * zz) * TMath::Cos(aa)
                      + (y0 + v0 * zz) * TMath::Sin(aa);
  const Double_t dsdz = u0 * TMath::Cos(aa) + v0 * TMath::Sin(aa);
  const Double_t coss = hit->IsHoneycomb()
    ? TMath::Cos(TMath::ATan(dsdz)) : 1.;
  const Double_t wp = hit->GetWirePosition();
  const Double_t ss = wp + (hit->GetLocalHitPos() - wp) / coss;
  return (ss - scal) * coss;
}

//_____________________________________________________________________________
Double_t
Pull(const DCLocalTrack& track, Int_t ihit)
{
  const Double_t res = Residual(track, ihit);
  if(!TMath::Finite(res)) return TMath::QuietNaN();

  const auto& hits = track.GetHitArray();
  if(ihit < 0 || ihit >= static_cast<Int_t>(hits.size()) || !hits[ihit])
    return TMath::QuietNaN();

  const Double_t sigma = gGeom.GetResolution(hits[ihit]->GetLayer());
  if(sigma <= 0.) return TMath::QuietNaN();
  return res / sigma;
}

}
