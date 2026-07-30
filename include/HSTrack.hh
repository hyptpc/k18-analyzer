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

  static const std::vector<Double_t>& VPZPlanes();

private:
  void ClearResults();
  void FillInvalidResults();

  Double_t m_xout;
  Double_t m_yout;
  Double_t m_uout;
  Double_t m_vout;
  Double_t m_momentum;
  Int_t m_status;

  std::vector<TVector3> m_vp_position;
  std::vector<TVector3> m_vp_momentum;
};

#endif
