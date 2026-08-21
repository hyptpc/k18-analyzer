// -*- C++ -*-
// Helpers for GenfitKpTopology next-stage analysis (V0 quality, Armenteros,
// dual-vertex DCA, q-frame kinematics). Header-only; no link dependency.
#ifndef KP_TOPOLOGY_ANALYSIS_HH
#define KP_TOPOLOGY_ANALYSIS_HH

#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector3.h>
#include <cmath>
#include <vector>

#include "Kinematics.hh"

namespace kptopo
{
inline const Double_t qnan = TMath::QuietNaN();

// Topology-first V0 quality: lower is better. Mass is NOT used.
// Components: daughter DCA (mm), pointing (1-cos), fiducial penalty, track quality.
inline Double_t V0QualityScore(Double_t dauDcaMm,
			       Double_t cosPointing,
			       Bool_t inFiducial,
			       Double_t chisqrSum,
			       Int_t nclustMin)
{
  Double_t s = 0.;
  if(TMath::IsNaN(dauDcaMm)) s += 1e3;
  else s += dauDcaMm;
  if(TMath::IsNaN(cosPointing)) s += 50.;
  else s += 20. * (1. - cosPointing);
  if(!inFiducial) s += 30.;
  if(!TMath::IsNaN(chisqrSum)) s += 0.5 * chisqrSum;
  if(nclustMin < 8) s += 5. * (8 - nclustMin);
  return s;
}

inline void ArmenterosPodolanski(const TVector3 &pPos, const TVector3 &pNeg,
				 Double_t &alpha, Double_t &qT)
{
  const TVector3 pV0 = pPos + pNeg;
  const Double_t pV0Mag = pV0.Mag();
  if(pV0Mag < 1e-9){ alpha = qnan; qT = qnan; return; }
  const TVector3 u = pV0.Unit();
  const Double_t qlPos = pPos.Dot(u);
  const Double_t qlNeg = pNeg.Dot(u);
  const Double_t den = qlPos + qlNeg;
  alpha = (std::abs(den) < 1e-12) ? qnan : (qlPos - qlNeg) / den;
  qT = (pPos - qlPos * u).Mag();
}

inline Double_t CosPointing(const TVector3 &primary, const TVector3 &decay,
			    const TVector3 &pV0)
{
  const TVector3 flight = decay - primary;
  if(flight.Mag() < 1e-6 || pV0.Mag() < 1e-6) return qnan;
  return flight.Unit().Dot(pV0.Unit());
}

inline Bool_t InTargetFiducial(const TVector3 &v, Double_t zTarget,
			       Double_t rMax=30., Double_t zHalf=60.)
{
  if(TMath::IsNaN(v.x()) || TMath::IsNaN(v.y()) || TMath::IsNaN(v.z())) return false;
  return (TMath::Hypot(v.x(), v.y()) < rMax)
    && (TMath::Abs(v.z() - zTarget) < zHalf);
}

// Project lab momentum onto q = p_K - p_forward.
inline void ProjectOnQ(const TVector3 &pLab, const TVector3 &q,
		       Double_t &pPar, Double_t &pT, Double_t &cosTh)
{
  const Double_t qMag = q.Mag();
  if(qMag < 1e-9 || TMath::IsNaN(pLab.x())){
    pPar = qnan; pT = qnan; cosTh = qnan; return;
  }
  const TVector3 u = q.Unit();
  pPar = pLab.Dot(u);
  pT = (pLab - pPar * u).Mag();
  cosTh = pLab.Unit().Dot(u);
}

// Robust multitrack vertex: iterative outlier rejection by DCA to fitted vertex.
// Arrays must already exclude K18/forward/beam/accidental/V0 daughters.
inline TVector3 RobustMultitrackVertex(Int_t nIn,
				       Double_t *x0, Double_t *y0,
				       Double_t *u0, Double_t *v0,
				       Double_t dcaCutMm,
				       Int_t maxIter,
				       Int_t &nUsed,
				       Bool_t &valid)
{
  valid = false;
  nUsed = 0;
  if(nIn < 1) return TVector3(qnan, qnan, qnan);
  if(nIn == 1){
    // Single track: no unique 3D vertex; mark invalid.
    return TVector3(qnan, qnan, qnan);
  }

  std::vector<Int_t> keep(nIn, 1);
  TVector3 vtx(qnan, qnan, qnan);
  for(Int_t iter=0; iter<maxIter; ++iter){
    Double_t xx[100], yy[100], uu[100], vv[100];
    Int_t n = 0;
    for(Int_t i=0; i<nIn && n<100; ++i){
      if(!keep[i]) continue;
      xx[n]=x0[i]; yy[n]=y0[i]; uu[n]=u0[i]; vv[n]=v0[i];
      ++n;
    }
    if(n < 2) break;
    vtx = Kinematics::MultitrackVertex(n, xx, yy, uu, vv);
    if(TMath::IsNaN(vtx.x())) break;

    // Approximate DCA of each straight-line track to vtx in xy+z using slopes.
    Int_t nOut = 0;
    for(Int_t i=0; i<nIn; ++i){
      if(!keep[i]) continue;
      const Double_t dx = x0[i] - vtx.x();
      const Double_t dy = y0[i] - vtx.y();
      // ignore z slope residual at first order; use transverse impact
      const Double_t dca = TMath::Hypot(dx, dy);
      if(dca > dcaCutMm){ keep[i]=0; ++nOut; }
    }
    nUsed = 0;
    for(Int_t i=0; i<nIn; ++i) if(keep[i]) ++nUsed;
    if(nOut==0 || nUsed<2) break;
  }
  valid = (nUsed >= 2 && !TMath::IsNaN(vtx.x()));
  if(!valid) return TVector3(qnan, qnan, qnan);
  return vtx;
}

struct V0Cand {
  Int_t idPos{-1};
  Int_t idNeg{-1};
  Double_t mass{qnan};
  Double_t dauDca{qnan};
  Double_t decayLen{qnan};
  Double_t cosPoint{qnan};
  Double_t alpha{qnan};
  Double_t qT{qnan};
  Double_t quality{qnan};
  Double_t vtx_x{qnan}, vtx_y{qnan}, vtx_z{qnan};
  Double_t mom{qnan}, mom_x{qnan}, mom_y{qnan}, mom_z{qnan};
  Double_t pPos{qnan}, pNeg{qnan};
  Int_t region{0}; // 0 none, 1 signal, 2 left SB, 3 right SB
  Bool_t fiducial{false};
};

inline Int_t MassRegionLambda(Double_t m, Double_t m0,
			      Double_t sigWin, Double_t sbLow, Double_t sbHigh)
{
  const Double_t dm = m - m0;
  const Double_t adm = TMath::Abs(dm);
  if(adm < sigWin) return 1;
  if(adm >= sbLow && adm < sbHigh) return (dm < 0) ? 2 : 3;
  return 0;
}

inline Int_t MassRegionK0S(Double_t m, Double_t m0,
			   Double_t sigWin, Double_t sbLow, Double_t sbHigh)
{
  return MassRegionLambda(m, m0, sigWin, sbLow, sbHigh);
}

// Sideband scale for flat BG under peak: signal half-width / sideband half-total-width.
// Signal |dm|<sigWin → width 2*sigWin; sideband |dm| in [sbLow,sbHigh) → width 2*(sbHigh-sbLow).
inline Double_t FlatSidebandScale(Double_t sigWin, Double_t sbLow, Double_t sbHigh)
{
  const Double_t wSig = 2. * sigWin;
  const Double_t wSb  = 2. * (sbHigh - sbLow);
  if(wSb <= 0) return qnan;
  return wSig / wSb;
}

} // namespace kptopo

#endif
