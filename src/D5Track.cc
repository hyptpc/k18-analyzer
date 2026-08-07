// -*- C++ -*-

#include "D5Track.hh"

#include <iostream>

#include <TMath.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "D5TransMatrix.hh"

#define DEBUG_D5 0

namespace
{
  auto& gD5Mtx = D5TransMatrix::GetInstance();

  void
  BuildTransportInput(Double_t x0, Double_t u0, Double_t y0, Double_t v0,
                      Double_t delta, Double_t z_ref_in, Double_t* in)
  {
    in[0] = x0 + u0 * z_ref_in;
    in[1] = TMath::ATan(u0) * 1000.0;
    in[2] = y0 + v0 * z_ref_in;
    in[3] = TMath::ATan(v0) * 1000.0;
    in[4] = delta;
  }

  Double_t
  WireResidual(const DCLTrackHit* hit, Double_t x_pred, Double_t y_pred)
  {
    const Double_t a = hit->GetTiltAngle() * TMath::DegToRad();
    const Double_t scal =
      x_pred * TMath::Cos(a) + y_pred * TMath::Sin(a);
    return hit->GetLocalHitPos() - scal;
  }

  Double_t
  WireResidualChi2(const DCLTrackHit* hit, Double_t x_pred, Double_t y_pred)
  {
    const Double_t resi = WireResidual(hit, x_pred, y_pred);
    const Double_t reso = hit->GetResolution();
    return (resi * resi) / (reso * reso);
  }

  Double_t
  AddTrackWireChi2(const DCLocalTrack* track,
                   Double_t x0, Double_t u0, Double_t y0, Double_t v0)
  {
    Double_t chisqr = 0.0;
    const Int_t nhit = track->GetNHit();
    for (Int_t i=0; i<nhit; ++i) {
      const auto* hit = track->GetHit(i);
      const Double_t z_hit = hit->GetZ();
      const Double_t x_pred = x0 + u0 * z_hit;
      const Double_t y_pred = y0 + v0 * z_hit;
      chisqr += WireResidualChi2(hit, x_pred, y_pred);
    }
    return chisqr;
  }
}

//_____________________________________________________________________________
D5Track::D5Track(const DCLocalTrack* blc1, const DCLocalTrack* blc2)
  : m_trk_blc1(blc1),
    m_trk_blc2(blc2),
    m_is_fitted(false),
    m_delta(0.0),
    m_momentum(0.0),
    m_d5_chi2(-1.0),
    m_d5_ndf(0),
    m_fit_x0(0.0),
    m_fit_u0(0.0),
    m_fit_y0(0.0),
    m_fit_v0(0.0),
    m_p0(gD5Mtx.GetCentralMomentum()),
    m_mtxout_x(0.0),
    m_mtxout_u(0.0),
    m_mtxout_y(0.0),
    m_mtxout_v(0.0)
{
}

//_____________________________________________________________________________
D5Track::~D5Track()
{
}

//_____________________________________________________________________________
Int_t
D5Track::CalcNDF() const
{
  const Int_t nhit_blc1 = m_trk_blc1 ? m_trk_blc1->GetNHit() : 0;
  const Int_t nhit_blc2 = m_trk_blc2 ? m_trk_blc2->GetNHit() : 0;
  return nhit_blc1 + nhit_blc2 - kNumFitPar;
}

//_____________________________________________________________________________
void
D5Track::SetupMinuitVariables(ROOT::Minuit2::Minuit2Minimizer& minimizer) const
{
  minimizer.SetVariable(0, "x_blc1", m_trk_blc1->GetX0(), 0.1);
  minimizer.SetVariable(1, "u_blc1", m_trk_blc1->GetU0(), 0.001);
  minimizer.SetVariable(2, "y_blc1", m_trk_blc1->GetY0(), 0.1);
  minimizer.SetVariable(3, "v_blc1", m_trk_blc1->GetV0(), 0.001);
  minimizer.SetVariable(4, "delta", 0.0, 0.01);
}

//_____________________________________________________________________________
void
D5Track::ReadFitParameters(const Double_t* xs)
{
  m_par.assign(xs, xs + kNumFitPar);
  m_fit_x0 = xs[0];
  m_fit_u0 = xs[1];
  m_fit_y0 = xs[2];
  m_fit_v0 = xs[3];
  m_delta = xs[4];
}

//_____________________________________________________________________________
Double_t
D5Track::GetD5ZIn() const
{
  return gD5Mtx.IsRefPlaneReady() ? gD5Mtx.GetD5ZIn() : D5TransMatrix::RefZIn;
}

//_____________________________________________________________________________
Double_t
D5Track::GetD5ZOut() const
{
  return gD5Mtx.IsRefPlaneReady() ? gD5Mtx.GetD5ZOut() : D5TransMatrix::RefZOut;
}

//_____________________________________________________________________________
Bool_t
D5Track::CalcMomentum()
{
  if (!m_trk_blc1 || !m_trk_blc2 || !gD5Mtx.IsReady()) return false;

  m_d5_ndf = CalcNDF();
  if (m_d5_ndf < 1) return false;

  ROOT::Minuit2::Minuit2Minimizer minimizer("Migrad");
  minimizer.SetMaxFunctionCalls(1000000);
  minimizer.SetTolerance(0.001);
  minimizer.SetPrintLevel(0);

  ROOT::Math::Functor f(this, &D5Track::operator(), kNumFitPar);
  minimizer.SetFunction(f);
  SetupMinuitVariables(minimizer);
  minimizer.Minimize();

  const Double_t* xs = minimizer.X();
  ReadFitParameters(xs);

  m_momentum = m_p0 * (1.0 + m_delta/100.0);
  m_d5_chi2 = minimizer.MinValue();

  const Double_t z_ref_in = gD5Mtx.GetD5ZIn();
  Double_t in[5];
  BuildTransportInput(m_fit_x0, m_fit_u0, m_fit_y0, m_fit_v0,
                      m_delta, z_ref_in, in);
  Double_t out[4];
  gD5Mtx.Transport(in, out);
  m_mtxout_x = out[0];
  m_mtxout_u = TMath::Tan(out[1] / 1000.0); // mrad -> dx/dz
  m_mtxout_y = out[2];
  m_mtxout_v = TMath::Tan(out[3] / 1000.0);

  m_is_fitted = true;

#if DEBUG_D5
  const Double_t z_ref_out = gD5Mtx.GetD5ZOut();
  const Int_t nhit_blc1 = m_trk_blc1->GetNHit();
  const Int_t nhit_blc2 = m_trk_blc2->GetNHit();
  std::cout << "[D5 Layer Debug] | D5ZIn: " << z_ref_in
            << " | D5ZOut: " << z_ref_out << std::endl;
  std::cout << "   BLC1 hits: " << nhit_blc1 << std::endl;
  for (Int_t i=0; i<nhit_blc1; ++i) {
    const auto* hit = m_trk_blc1->GetHit(i);
    Double_t z_hit = hit->GetZ();
    Double_t x_pred = m_fit_x0 + m_fit_u0 * z_hit;
    Double_t y_pred = m_fit_y0 + m_fit_v0 * z_hit;
    const Double_t a = hit->GetTiltAngle() * TMath::DegToRad();
    const Double_t scal =
      x_pred * TMath::Cos(a) + y_pred * TMath::Sin(a);
    std::cout << "   BLC1 L#" << hit->GetLayer()
              << " | Z: " << z_hit
              << " | meas: " << hit->GetLocalHitPos()
              << " | pred: " << scal
              << " | resi: " << (hit->GetLocalHitPos() - scal)
              << std::endl;
  }
  std::cout << "   BLC2 hits: " << nhit_blc2 << std::endl;
  for (Int_t i=0; i<nhit_blc2; ++i) {
    const auto* hit = m_trk_blc2->GetHit(i);
    Double_t z_hit = hit->GetZ();
    Double_t dz = z_hit - z_ref_out;
    Double_t x_pred = m_mtxout_x + m_mtxout_u * dz;
    Double_t y_pred = m_mtxout_y + m_mtxout_v * dz;
    const Double_t a = hit->GetTiltAngle() * TMath::DegToRad();
    const Double_t scal =
      x_pred * TMath::Cos(a) + y_pred * TMath::Sin(a);
    std::cout << "   BLC2 L#" << hit->GetLayer()
              << " | Z: " << z_hit
              << " | meas: " << hit->GetLocalHitPos()
              << " | pred: " << scal
              << " | resi: " << (hit->GetLocalHitPos() - scal)
              << std::endl;
  }
#endif

  return true;
}

//_____________________________________________________________________________
Double_t
D5Track::GetResidualX() const
{
  if (!m_is_fitted) return -999.0;
  const Double_t z_ref_out = GetD5ZOut();
  return m_trk_blc2->GetX(z_ref_out) - m_mtxout_x;
}

//_____________________________________________________________________________
Double_t
D5Track::GetResidualY() const
{
  if (!m_is_fitted) return -999.0;
  const Double_t z_ref_out = GetD5ZOut();
  return m_trk_blc2->GetY(z_ref_out) - m_mtxout_y;
}

//_____________________________________________________________________________
Double_t
D5Track::operator()(const Double_t* par)
{
  return CalcChi2(par);
}

//_____________________________________________________________________________
Double_t
D5Track::CalcChi2(const Double_t* par)
{
  const Double_t z_ref_in  = gD5Mtx.GetD5ZIn();
  const Double_t z_ref_out = gD5Mtx.GetD5ZOut();

  const Double_t x0 = par[0];
  const Double_t u0 = par[1];
  const Double_t y0 = par[2];
  const Double_t v0 = par[3];
  const Double_t delta = par[4];

  Double_t chisqr = AddTrackWireChi2(m_trk_blc1, x0, u0, y0, v0);

  Double_t in[5];
  BuildTransportInput(x0, u0, y0, v0, delta, z_ref_in, in);
  Double_t out[4];
  if (!gD5Mtx.Transport(in, out)) return 1.0e10;

  const Double_t mu = TMath::Tan(out[1] / 1000.0);
  const Double_t mv = TMath::Tan(out[3] / 1000.0);

  const Int_t nhit = m_trk_blc2->GetNHit();
  for (Int_t i=0; i<nhit; ++i) {
    const auto* hit = m_trk_blc2->GetHit(i);
    const Double_t z_hit = hit->GetZ();
    const Double_t dz = z_hit - z_ref_out;
    const Double_t x_pred = out[0] + mu * dz;
    const Double_t y_pred = out[2] + mv * dz;
    chisqr += WireResidualChi2(hit, x_pred, y_pred);
  }

  return chisqr;
}
