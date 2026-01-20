// -*- C++ -*-

#ifndef CHERENKOV_HIT_HH
#define CHERENKOV_HIT_HH

#include <TString.h>

#include "DebugCounter.hh"
#include "HodoHit.hh"
#include "HodoRawHit.hh"

//_____________________________________________________________________________
class CherenkovHit : public HodoHit
{
public:
  static const TString& ClassName();
  CherenkovHit(HodoRawHit* rhit, Double_t max_time_diff=10.);
  ~CherenkovHit();

public:
  Bool_t   Calculate() override;
  Double_t Npe(Int_t j=0) const override;  // DeltaE-like. 1ch: that ch. 2ch: U+D. KVC: ch0-3. Offline sum uses Npe; Event/User sum over hits for BAC.
  Double_t NpeSum(Int_t j=0) const;        // Online sum (hardware SUM). KVC: ch4. BAC: seg4 only, else NaN.
  Double_t GetNpe(Int_t ch, Int_t j=0) const;
  const std::vector<Double_t>& GetArrayNpe(Int_t ch) const;
  Double_t DeltaE(Int_t j=0) const override { return Npe(j); }  // Compatibility: returns Npe (semantics differ from HodoHit's dE)
  Double_t GetDeltaEHighGain(Int_t i, Int_t j=0) const override;  // Returns Npe; i=EChannel. KVC: i=2 -> kSUM(ch=4).
  virtual Bool_t ReCalc(Bool_t applyRecursively=false) override { return Calculate(); }

private:
  CherenkovHit(const CherenkovHit&);
  CherenkovHit& operator=(const CherenkovHit&);

protected:
  // [ch][hit]. AdcHigh only (no m_npe_low; AdcLow unused as in Hodo)
  std::vector<std::vector<Double_t>> m_npe_high;
};

//_____________________________________________________________________________
inline const TString&
CherenkovHit::ClassName()
{
  static TString s_name("CherenkovHit");
  return s_name;
}

#endif
