// -*- C++ -*-

#include "HSTrack.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

#include <TMath.h>

#include "FieldMan.hh"
#include "RungeKuttaUtilities.hh"
#include "TPCPadHelper.hh"

namespace
{
const Double_t kInvalidVal = -9999.;
const Double_t kHSStepZ =   10.;
const Int_t    kMaxStep = 10000;
const Bool_t   kPrintVPHSField = false;

// RK q = kShsMapChargeSign * charge / p (FieldMan SHS map).
constexpr Double_t kShsMapChargeSign = -1.;

// Legacy TPCHS / default sampling: upstream residuals + 10 mm through TPC.
const std::vector<Double_t>& VPZPlanesImpl()
{
  static const std::vector<Double_t> planes = [] {
    std::vector<Double_t> out = {-243., -235., -225., -214., -206.};
    for (Double_t z = -196.; z <= 244.; z += 10.) out.push_back(z);
    out.push_back(250.);
    return out;
  }();
  return planes;
}

// Kp option: same as standard outside the target band; 1 mm inside
// [Z_TARGET - TARGET_RADIUS, Z_TARGET + TARGET_RADIUS], always include Z_TARGET.
const std::vector<Double_t>& VPZPlanesDenseTarget()
{
  static const std::vector<Double_t> planes = [] {
    const Double_t z_lo = tpc::Z_TARGET - tpc::TARGET_RADIUS;
    const Double_t z_hi = tpc::Z_TARGET + tpc::TARGET_RADIUS;
    std::vector<Double_t> out = {-243., -235., -225., -214., -206.};
    for (Double_t z = -196.; z <= 244. + 1.e-9; z += 10.) {
      if (z >= z_lo - 1.e-9 && z <= z_hi + 1.e-9)
        continue;
      out.push_back(z);
    }
    for (Double_t z = z_lo; z <= z_hi + 1.e-9; z += 1.)
      out.push_back(z);
    out.push_back(tpc::Z_TARGET);
    out.push_back(250.);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(),
                          [](Double_t a, Double_t b) {
                            return TMath::Abs(a - b) < 1.e-6;
                          }),
              out.end());
    return out;
  }();
  return planes;
}

Bool_t IsInvalidVector(const TVector3& v)
{
  return v.X() == kInvalidVal && v.Y() == kInvalidVal && v.Z() == kInvalidVal;
}
}

//_____________________________________________________________________________
HSTrack::HSTrack(Double_t xout, Double_t yout,
                 Double_t uout, Double_t vout, Double_t p)
  : m_xout(xout),
    m_yout(yout),
    m_uout(uout),
    m_vout(vout),
    m_momentum(p),
    m_start_z(DefaultStartZ()),
    m_charge(-1),
    m_status(kInit),
    m_dense_target_vp(false),
    m_vp_position(),
    m_vp_momentum()
{
  FillInvalidResults();
}

//_____________________________________________________________________________
const std::vector<Double_t>&
HSTrack::VPZPlanes()
{
  return VPZPlanesImpl();
}

//_____________________________________________________________________________
const std::vector<Double_t>&
HSTrack::ActiveVPZPlanes() const
{
  return m_dense_target_vp ? VPZPlanesDenseTarget() : VPZPlanesImpl();
}

//_____________________________________________________________________________
void
HSTrack::ClearResults()
{
  m_vp_position.clear();
  m_vp_momentum.clear();
}

//_____________________________________________________________________________
void
HSTrack::FillInvalidResults()
{
  ClearResults();
  const auto& planes = ActiveVPZPlanes();
  m_vp_position.assign(planes.size(), TVector3(kInvalidVal, kInvalidVal, kInvalidVal));
  m_vp_momentum.assign(planes.size(), TVector3(kInvalidVal, kInvalidVal, kInvalidVal));
}

//_____________________________________________________________________________
Bool_t
HSTrack::Propagate()
{
  FillInvalidResults();
  m_status = kInit;

  if(!FieldMan::GetInstance().IsReady()){
    m_status = kNoField;
    return false;
  }
  if(m_momentum <= 0. || !std::isfinite(m_momentum) ||
     !std::isfinite(m_xout) || !std::isfinite(m_yout) ||
     !std::isfinite(m_uout) || !std::isfinite(m_vout) ||
     !std::isfinite(m_start_z) ||
     (m_charge != -1 && m_charge != 1)){
    m_status = kInvalidInput;
    return false;
  }

  const auto& planes = ActiveVPZPlanes();
  // Start at (x,y,u,v) on m_start_z (default VO ≈ -1300 mm).
  const RKCordParameter ini(m_xout, m_yout, m_start_z, m_uout, m_vout,
                            kShsMapChargeSign * static_cast<Double_t>(m_charge)
                              / m_momentum);
  RKTrajectoryPoint prev(ini,
                         1., 0., 0., 0., 0.,
                         0., 1., 0., 0., 0.,
                         0., 0., 1., 0., 0.,
                         0., 0., 0., 1., 0.,
                         0.);

  std::size_t ivp = 0;
  for(Int_t istep=0; istep<kMaxStep && ivp<planes.size(); ++istep){
    const RKTrajectoryPoint next = RK::PropagateOnce(kHSStepZ, prev);
    const ThreeVector pos1 = prev.PositionInGlobal();
    const ThreeVector pos2 = next.PositionInGlobal();

    while(ivp<planes.size() &&
          (pos1.z() - planes[ivp])*(pos2.z() - planes[ivp]) <= 0.){
      const Double_t denom = pos2.z() - pos1.z();
      if(TMath::Abs(denom) < 1.e-9)
        break;
      const Double_t f = (planes[ivp] - pos1.z())/denom;
      const ThreeVector mom1 = prev.MomentumInGlobal();
      const ThreeVector mom2 = next.MomentumInGlobal();
      const ThreeVector pos = pos1 + f*(pos2 - pos1);
      const ThreeVector mom = mom1 + f*(mom2 - mom1);
      m_vp_position[ivp].SetXYZ(pos.x(), pos.y(), planes[ivp]);
      m_vp_momentum[ivp].SetXYZ(mom.x(), mom.y(), mom.z());
      if(kPrintVPHSField){
        const ThreeVector field = FieldMan::GetInstance().GetField(pos);
        std::cout << "[HSTrack] VPHS" << ivp+1
                  << " pos=(" << pos.x() << ", " << pos.y()
                  << ", " << planes[ivp] << ")"
                  << " B=(" << field.x() << ", " << field.y()
                  << ", " << field.z() << ")"
                  << std::endl;
      }
      ++ivp;
    }
    prev = next;
  }

  if (ivp == planes.size()) {
    std::vector<TVector3> valid_position;
    std::vector<TVector3> valid_momentum;
    valid_position.reserve(m_vp_position.size());
    valid_momentum.reserve(m_vp_momentum.size());
    for (std::size_t i = 0; i < m_vp_position.size(); ++i) {
      const TVector3& pos = m_vp_position[i];
      const Bool_t in_target_band =
        TMath::Abs(pos.Z() - tpc::Z_TARGET) <= tpc::TARGET_RADIUS + 1.e-6;
      // Default: require an actual pad (skip inter-layer gaps).
      // Dense-target mode: keep padless VPs inside the target band (beam @ vtx).
      if (!(m_dense_target_vp && in_target_band)
          && tpc::FindPadID(pos.Z(), pos.X()) < 0)
        continue;
      valid_position.push_back(pos);
      valid_momentum.push_back(m_vp_momentum[i]);
    }
    m_vp_position.swap(valid_position);
    m_vp_momentum.swap(valid_momentum);
  }
  m_status = (ivp == planes.size()) ? kPassed : kExceedMaxStep;
  return IsPassed();
}

//_____________________________________________________________________________
std::vector<Double_t>
HSTrack::VPX() const
{
  std::vector<Double_t> out;
  out.reserve(m_vp_position.size());
  for(const auto& pos : m_vp_position) out.push_back(pos.X());
  return out;
}

//_____________________________________________________________________________
std::vector<Double_t>
HSTrack::VPY() const
{
  std::vector<Double_t> out;
  out.reserve(m_vp_position.size());
  for(const auto& pos : m_vp_position) out.push_back(pos.Y());
  return out;
}

//_____________________________________________________________________________
std::vector<Double_t>
HSTrack::VPZ() const
{
  std::vector<Double_t> out;
  out.reserve(m_vp_position.size());
  for(const auto& pos : m_vp_position) out.push_back(pos.Z());
  return out;
}

//_____________________________________________________________________________
std::vector<Double_t>
HSTrack::VPU() const
{
  std::vector<Double_t> out;
  out.reserve(m_vp_momentum.size());
  for(const auto& mom : m_vp_momentum)
    out.push_back(!IsInvalidVector(mom) && mom.Z() != 0. ? mom.X()/mom.Z() : kInvalidVal);
  return out;
}

//_____________________________________________________________________________
std::vector<Double_t>
HSTrack::VPV() const
{
  std::vector<Double_t> out;
  out.reserve(m_vp_momentum.size());
  for(const auto& mom : m_vp_momentum)
    out.push_back(!IsInvalidVector(mom) && mom.Z() != 0. ? mom.Y()/mom.Z() : kInvalidVal);
  return out;
}

//_____________________________________________________________________________
std::vector<Double_t>
HSTrack::VPP() const
{
  std::vector<Double_t> out;
  out.reserve(m_vp_momentum.size());
  for(const auto& mom : m_vp_momentum)
    out.push_back(!IsInvalidVector(mom) ? mom.Mag() : kInvalidVal);
  return out;
}
