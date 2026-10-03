// -*- C++ -*-

#ifndef HSTRACK_HH
#define HSTRACK_HH

#include <vector>

#include <Rtypes.h>
#include <TVector3.h>

class HSTrack
{
public:
  enum Status {
    kInit = 0,
    kPassed,
    kNoField,
    kInvalidInput,
    kExceedMaxStep,
    nStatus
  };

  HSTrack(Double_t xout, Double_t yout,
          Double_t uout, Double_t vout, Double_t p);

  // Default false: legacy 10 mm VP grid + pad-required filter (other users).
  // true: 1 mm VP in target z-band (incl. Z_TARGET) and keep padless VPs there.
  void EnableDenseTargetVP(Bool_t enable = true) { m_dense_target_vp = enable; }
  Bool_t DenseTargetVP() const { return m_dense_target_vp; }

  // RK start z [mm]. Default is DefaultStartZ() (-1300). Constructor (x,y,u,v)
  // are interpreted at this plane.
  static Double_t DefaultStartZ() { return -1300.; }
  void SetStartZ(Double_t z) { m_start_z = z; }
  Double_t StartZ() const { return m_start_z; }

  // Beam charge (±1). Default -1. q = kShsMapChargeSign * charge / p in Propagate.
  void SetCharge(Int_t charge) { m_charge = charge; }
  Int_t Charge() const { return m_charge; }

  Bool_t Propagate();

  Int_t StatusCode() const { return m_status; }
  Bool_t IsPassed() const { return m_status == kPassed; }

  const std::vector<TVector3>& VPPosition() const { return m_vp_position; }
  const std::vector<TVector3>& VPMomentum() const { return m_vp_momentum; }

  std::vector<Double_t> VPX() const;
  std::vector<Double_t> VPY() const;
  std::vector<Double_t> VPZ() const;
  std::vector<Double_t> VPU() const;
  std::vector<Double_t> VPV() const;
  std::vector<Double_t> VPP() const;

  // Legacy / default VP z list (no dense target band).
  static const std::vector<Double_t>& VPZPlanes();

private:
  void ClearResults();
  void FillInvalidResults();
  const std::vector<Double_t>& ActiveVPZPlanes() const;

  Double_t m_xout;
  Double_t m_yout;
  Double_t m_uout;
  Double_t m_vout;
  Double_t m_momentum;
  Double_t m_start_z;
  Int_t m_charge;
  Int_t m_status;
  Bool_t m_dense_target_vp;

  std::vector<TVector3> m_vp_position;
  std::vector<TVector3> m_vp_momentum;
};

#endif
