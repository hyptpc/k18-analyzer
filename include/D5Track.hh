// -*- C++ -*-

#ifndef D5_TRACK_HH
#define D5_TRACK_HH

#include <TString.h>
#include <vector>

#include "D5TransMatrix.hh"

class DCLocalTrack;

namespace ROOT {
namespace Minuit2 {
class Minuit2Minimizer;
}
}

//_____________________________________________________________________________
class D5Track
{
public:
  static const TString& ClassName();
  D5Track(const DCLocalTrack* blc1, const DCLocalTrack* blc2);
  ~D5Track();

private:
  D5Track();
  D5Track(const D5Track&);
  D5Track& operator=(const D5Track&);

private:
  const DCLocalTrack* m_trk_blc1;
  const DCLocalTrack* m_trk_blc2;

  Bool_t   m_is_fitted;
  Double_t m_delta;
  Double_t m_momentum;
  Double_t m_d5_chi2;
  Int_t    m_d5_ndf;

  Double_t m_fit_x0;
  Double_t m_fit_u0;
  Double_t m_fit_y0;
  Double_t m_fit_v0;

  Double_t m_p0;

  Double_t m_mtxout_x;
  Double_t m_mtxout_u; // dx/dz (same as BLC), not mrad
  Double_t m_mtxout_y;
  Double_t m_mtxout_v; // dy/dz (same as BLC), not mrad

  std::vector<Double_t> m_par;

public:
  Bool_t   CalcMomentum();
  Bool_t   IsFitted() const { return m_is_fitted; }

  Double_t GetDelta() const { return m_delta; }
  Double_t GetMomentum() const { return m_momentum; }
  Double_t GetD5Chi2() const { return m_d5_chi2; }
  Int_t    GetD5Ndf() const { return m_d5_ndf; }
  Double_t GetD5Chi2Ndf() const
  {
    return m_d5_ndf > 0 ? m_d5_chi2 / static_cast<Double_t>(m_d5_ndf) : -1.0;
  }

  Double_t GetFitX0() const { return m_fit_x0; }
  Double_t GetFitU0() const { return m_fit_u0; }
  Double_t GetFitY0() const { return m_fit_y0; }
  Double_t GetFitV0() const { return m_fit_v0; }

  const DCLocalTrack* GetTrkBlc1() const { return m_trk_blc1; }
  const DCLocalTrack* GetTrkBlc2() const { return m_trk_blc2; }

  Double_t GetMtxoutX() const { return m_mtxout_x; }
  Double_t GetMtxoutU() const { return m_mtxout_u; } // dx/dz (BLC)
  Double_t GetMtxoutY() const { return m_mtxout_y; }
  Double_t GetMtxoutV() const { return m_mtxout_v; } // dy/dz (BLC)

  Double_t GetD5ZIn() const;
  Double_t GetD5ZOut() const;

  Double_t GetResidualX() const;
  Double_t GetResidualY() const;

  Double_t   operator()(const Double_t* par);

private:
  static constexpr Int_t kNumFitPar = 5;

  Int_t      CalcNDF() const;
  void       SetupMinuitVariables(ROOT::Minuit2::Minuit2Minimizer& minimizer) const;
  void       ReadFitParameters(const Double_t* xs);
  Double_t   CalcChi2(const Double_t* par);
};

//_____________________________________________________________________________
inline const TString&
D5Track::ClassName()
{
  static TString s_name("D5Track");
  return s_name;
}

#endif
