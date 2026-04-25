// -*- C++ -*-

#include "TPCReconstructor.hh"

#include <iomanip>
#include <iostream>
#include <TPDGCode.h>
#include <TLorentzVector.h>

#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "Kinematics.hh"
#include "UserParamMan.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"

#define TPC_RECO_ENABLE_LAMBDA_REFIT 0

namespace
{
  const Double_t PPI_CLOSE_MM = 10.; // mm, same scale as TPCVertex / TPCTrackSearch
  const Double_t LAMBDA_MASS_WINDOW = 0.1; // GeV/c^2
  const Double_t K0_MASS_WINDOW = 0.05; // GeV/c^2
  const TVector3 VERTEX_RES(0.6, 0.6, 1.0); // just follow e42 code

  // Common preselection for any 2-track decay reconstruction.
  Bool_t PassScatterTrackFlags(const TPCVertex* vertex)
  {
    if (!vertex || vertex->GetNTracks() < 2)
      return false;
    const UInt_t scatter_track_flags = vertex->GetScatterTrackFlags();
    if ((scatter_track_flags & TPCVertex::kScatterTrackFlagsUnknown) != 0)
      return false;
    constexpr UInt_t k_reject_mask =
      TPCVertex::kTrack1IsK18 | TPCVertex::kTrack1IsBeam | TPCVertex::kTrack1IsAccidental |
      TPCVertex::kTrack2IsK18 | TPCVertex::kTrack2IsBeam | TPCVertex::kTrack2IsAccidental;
    if ((scatter_track_flags & k_reject_mask) != 0)
      return false;
    return true;
  }

  Bool_t PassLambdaVertexPreselection(const TPCVertex* vertex)
  {
    if (!vertex || !vertex->IsCalculated() || vertex->GetNTracks() < 2)
      return false;
    if (!PassScatterTrackFlags(vertex))
      return false;
    if (vertex->GetClosestDist() > PPI_CLOSE_MM)
      return false;
    return true;
  }

  // p_idx = proton-like, pi_idx = pi-; writes invariant mass when returning true.
  Bool_t CheckLambda(const TPCVertex* vertex, Int_t p_idx, Int_t pi_idx, Double_t& mass)
  {
    const Int_t p_charge = vertex->GetTrackCharge(p_idx);
    const Int_t pi_charge = vertex->GetTrackCharge(pi_idx);
    const Int_t p_pid = vertex->GetTrackPid(p_idx);
    const Int_t pi_pid = vertex->GetTrackPid(pi_idx);
    if ((p_pid & 0x4) != 0x4 || p_charge != 1)
      return false;
    if ((pi_pid & 0x1) != 0x1 || pi_charge != -1)
      return false;

    const TVector3 p_mom = vertex->GetTrackMom(p_idx);
    const TVector3 pi_mom = vertex->GetTrackMom(pi_idx);
    const Double_t m_proton = pdg::ProtonMass();
    const Double_t m_pion = pdg::PionMass();
    const TLorentzVector lv_p(p_mom.X(), p_mom.Y(), p_mom.Z(),
                              TMath::Sqrt(p_mom.Mag2() + m_proton * m_proton));
    const TLorentzVector lv_pi(pi_mom.X(), pi_mom.Y(), pi_mom.Z(),
                               TMath::Sqrt(pi_mom.Mag2() + m_pion * m_pion));
    const TLorentzVector lv_lambda = lv_p + lv_pi;
    mass = lv_lambda.M();
    return TMath::Abs(mass - pdg::LambdaMass()) <= LAMBDA_MASS_WINDOW;
  }

}

//_____________________________________________________________________________
TPCRecoCandidate::TPCRecoCandidate(Int_t mother_pdg, Double_t mass,
                                   const TVector3& vertex, const TVector3& momentum,
                                   Double_t closest_dist,
                                   Int_t track_id_1, Int_t track_id_2,
                                   Int_t track_1_pid, Int_t track_2_pid,
                                   Int_t track_1_charge, Int_t track_2_charge)
  : m_mother_pdg(mother_pdg),
    m_mass(mass),
    m_vertex(vertex),
    m_momentum(momentum),
    m_closest_dist(closest_dist),
    m_track_id_1(track_id_1),
    m_track_id_2(track_id_2),
    m_track_1_pid(track_1_pid),
    m_track_2_pid(track_2_pid),
    m_track_1_charge(track_1_charge),
    m_track_2_charge(track_2_charge)
{
}

//_____________________________________________________________________________
void
TPCRecoCandidate::Print(const TString& label) const
{
  const std::ios::fmtflags old_flags = std::cout.flags();
  const std::streamsize old_precision = std::cout.precision();
  std::cout << std::fixed << std::setprecision(4);

  std::cout << "#D";
  if (!label.IsNull())
    std::cout << " [" << label << "]";
  std::cout << "\n"
            << "   pdg   = " << m_mother_pdg << "\n"
            << "   mass  = " << m_mass << "\n"
            << "   dist  = " << m_closest_dist << "\n"
            << "   trkId = (" << m_track_id_1
            << (m_track_1_charge > 0 ? "+" : (m_track_1_charge < 0 ? "-" : "0")) << ","
            << m_track_id_2
            << (m_track_2_charge > 0 ? "+" : (m_track_2_charge < 0 ? "-" : "0")) << ")\n"
            << "   pid   = (" << m_track_1_pid << "," << m_track_2_pid << ")\n"
            << "   vtx   = (" << m_vertex.X() << "," << m_vertex.Y() << "," << m_vertex.Z() << ")\n"
            << "   mom   = (" << m_momentum.X() << "," << m_momentum.Y() << "," << m_momentum.Z() << ")\n"
            << "   |p|   = " << m_momentum.Mag()
            << std::endl;

  std::cout.flags(old_flags);
  std::cout.precision(old_precision);
}

//_____________________________________________________________________________
TPCReconstructor::TPCReconstructor()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCReconstructor::~TPCReconstructor()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCReconstructor::Clear()
{
}

//_____________________________________________________________________________
void
TPCReconstructor::ReconstructLambda(TPCVertex* vertex,
                                     const std::vector<TPCLocalTrackHelix*>* helix_tracks)
{
  if (!PassLambdaVertexPreselection(vertex))
    return;

  const Double_t closest_dist_orig = vertex->GetClosestDist();

  vertex->SetHasLambdaCandidate(false);

  for (Int_t i = 0; i < 2; ++i) {
    const Int_t p_idx = i;
    const Int_t pi_idx = 1 - i;
    Double_t mass = TMath::QuietNaN();
    if (!CheckLambda(vertex, p_idx, pi_idx, mass))
      continue;

    // After Lambda is identified, optionally refit helices to the vertex (same order as legacy).
    TVector3 lambda_vtx = vertex->GetVertex();
    Double_t closest_dist = closest_dist_orig;
    TVector3 out_mom1 = vertex->GetTrackMom(0);
    TVector3 out_mom2 = vertex->GetTrackMom(1);
    Double_t out_mass = mass;
#if TPC_RECO_ENABLE_LAMBDA_REFIT
    TVector3 refit_vtx;
    TVector3 refit_mom_1;
    TVector3 refit_mom_2;
    Double_t refit_dist = TMath::QuietNaN();
    if (RefitLambdaWithVertex(vertex, helix_tracks, refit_vtx, refit_mom_1, refit_mom_2, refit_dist)) {
      lambda_vtx = refit_vtx;
      closest_dist = refit_dist;
      out_mom1 = refit_mom_1;
      out_mom2 = refit_mom_2;
      const TVector3 p_mom_f = (p_idx == 0) ? out_mom1 : out_mom2;
      const TVector3 pi_mom_f = (pi_idx == 0) ? out_mom1 : out_mom2;
      const Double_t m_proton = pdg::ProtonMass();
      const Double_t m_pion = pdg::PionMass();
      const TLorentzVector lv_pf(p_mom_f.X(), p_mom_f.Y(), p_mom_f.Z(),
                                 TMath::Sqrt(p_mom_f.Mag2() + m_proton * m_proton));
      const TLorentzVector lv_pif(pi_mom_f.X(), pi_mom_f.Y(), pi_mom_f.Z(),
                                  TMath::Sqrt(pi_mom_f.Mag2() + m_pion * m_pion));
      out_mass = (lv_pf + lv_pif).M();
    }
#else
    (void)helix_tracks;
#endif
    const TVector3 p_mom_out = (p_idx == 0) ? out_mom1 : out_mom2;
    const TVector3 pi_mom_out = (pi_idx == 0) ? out_mom1 : out_mom2;
    const TVector3 mom = p_mom_out + pi_mom_out;
    vertex->SetHasLambdaCandidate(true);
    const TPCRecoCandidate cand(
      kLambda0,
      out_mass,
      lambda_vtx,
      mom,
      closest_dist,
      vertex->GetTrackId(0),
      vertex->GetTrackId(1),
      vertex->GetTrackPid(0),
      vertex->GetTrackPid(1),
      vertex->GetTrackCharge(0),
      vertex->GetTrackCharge(1)
    );
    vertex->AddRecoCandidate(cand);
    cand.Print("Lambda");
  }
}

//_____________________________________________________________________________
Bool_t
TPCReconstructor::RefitLambdaWithVertex(const TPCVertex* vertex,
                                       const std::vector<TPCLocalTrackHelix*>* helix_tracks,
                                       TVector3& out_vertex,
                                       TVector3& out_mom1,
                                       TVector3& out_mom2,
                                       Double_t& out_dist) const
{
  if (!vertex || !helix_tracks || vertex->GetNTracks() < 2)
    return false;

  const Int_t track_id_1 = vertex->GetTrackId(0);
  const Int_t track_id_2 = vertex->GetTrackId(1);
  if (track_id_1 < 0 || track_id_2 < 0 ||
      track_id_1 >= static_cast<Int_t>(helix_tracks->size()) ||
      track_id_2 >= static_cast<Int_t>(helix_tracks->size()))
    return false;

  TPCLocalTrackHelix* base_track_1 = helix_tracks->at(track_id_1);
  TPCLocalTrackHelix* base_track_2 = helix_tracks->at(track_id_2);
  if (!base_track_1 || !base_track_2)
    return false;

  auto* track_1 = new TPCLocalTrackHelix(base_track_1);
  if (!track_1->DoFitTrackwVertex(vertex->GetVertex(), VERTEX_RES)) {
    delete track_1;
    return false;
  }
  auto* track_2 = new TPCLocalTrackHelix(base_track_2);
  if (!track_2->DoFitTrackwVertex(vertex->GetVertex(), VERTEX_RES)) {
    delete track_1;
    delete track_2;
    return false;
  }

  const Double_t vertex_scan_range = UserParamMan::GetInstance().GetParameter("VertexScanRange");
  Double_t par_1[5];
  Double_t par_2[5];
  track_1->GetParam(par_1);
  track_2->GetParam(par_2);
  const Double_t scan_theta_1 = vertex_scan_range/par_1[3];
  const Double_t scan_theta_2 = vertex_scan_range/par_2[3];
  const Double_t range_theta_1[2] = {track_1->GetMint() - scan_theta_1,
                                     track_1->GetMaxt() + scan_theta_1};
  const Double_t range_theta_2[2] = {track_2->GetMint() - scan_theta_2,
                                     track_2->GetMaxt() + scan_theta_2};

  Double_t theta_1 = 0.0;
  Double_t theta_2 = 0.0;
  Double_t dist = 0.0;
  const TVector3 vertex_refit = Kinematics::VertexPointHelix(
    par_1, par_2,
    range_theta_1[0], range_theta_1[1],
    range_theta_2[0], range_theta_2[1],
    theta_1, theta_2, dist
  );
  if (!std::isfinite(dist) || dist > PPI_CLOSE_MM) {
    delete track_1;
    delete track_2;
    return false;
  }

  out_vertex = vertex_refit;
  out_dist = dist;
  out_mom1 = track_1->CalcHelixMom(par_1, theta_1);
  out_mom2 = track_2->CalcHelixMom(par_2, theta_2);
  delete track_1;
  delete track_2;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCReconstructor::HasLambdaCandidate(const TPCVertex* vertex)
{
  if (!PassLambdaVertexPreselection(vertex))
    return false;
  for (Int_t i = 0; i < 2; ++i) {
    Double_t mass = TMath::QuietNaN();
    if (CheckLambda(vertex, i, 1 - i, mass))
      return true;
  }
  return false;
}

//_____________________________________________________________________________
void
TPCReconstructor::ReconstructK0Short(TPCVertex* vertex)
{
  if (!vertex || !vertex->IsCalculated() || vertex->GetNTracks() < 2)
    return;
  if (!PassScatterTrackFlags(vertex))
    return;
  const Double_t dist = vertex->GetClosestDist();
  if (dist > PPI_CLOSE_MM)
    return;

  const Int_t q0 = vertex->GetTrackCharge(0);
  const Int_t q1 = vertex->GetTrackCharge(1);
  if (q0 * q1 >= 0)
    return;

  if ((vertex->GetTrackPid(0) & 1) != 1 || (vertex->GetTrackPid(1) & 1) != 1)
    return;

  const TVector3 p0 = vertex->GetTrackMom(0);
  const TVector3 p1 = vertex->GetTrackMom(1);
  TLorentzVector lv_pi0;
  TLorentzVector lv_pi1;
  lv_pi0.SetXYZM(p0.X(), p0.Y(), p0.Z(), pdg::PionMass());
  lv_pi1.SetXYZM(p1.X(), p1.Y(), p1.Z(), pdg::PionMass());
  const TLorentzVector lv_ks = lv_pi0 + lv_pi1;
  const Double_t m_ks = lv_ks.M();
  const Double_t m0_ks = pdg::Mass(kK0Short);
  if (m0_ks <= 0. || TMath::Abs(m_ks - m0_ks) > K0_MASS_WINDOW)
    return;

  const TVector3 vtx = vertex->GetVertex();
  const TVector3 mom = lv_ks.Vect();
  const TPCRecoCandidate cand(
    kK0Short,
    m_ks,
    vtx,
    mom,
    dist,
    vertex->GetTrackId(0),
    vertex->GetTrackId(1),
    vertex->GetTrackPid(0),
    vertex->GetTrackPid(1),
    vertex->GetTrackCharge(0),
    vertex->GetTrackCharge(1)
  );
  vertex->AddRecoCandidate(cand);
  cand.Print("K0Short");
}

//_____________________________________________________________________________
Bool_t
TPCReconstructor::HasK0ShortCandidate(const TPCVertex* vertex)
{
  // Fast check used for ambiguity/debug checks; no candidate is stored.
  if (!vertex || !vertex->IsCalculated() || vertex->GetNTracks() < 2)
    return false;
  if (!PassScatterTrackFlags(vertex))
    return false;

  const Double_t dist = vertex->GetClosestDist();
  if (dist > PPI_CLOSE_MM)
    return false;

  const Int_t q0 = vertex->GetTrackCharge(0);
  const Int_t q1 = vertex->GetTrackCharge(1);
  if (q0 * q1 >= 0)
    return false;

  if ((vertex->GetTrackPid(0) & 1) != 1 || (vertex->GetTrackPid(1) & 1) != 1)
    return false;

  const TVector3 p0 = vertex->GetTrackMom(0);
  const TVector3 p1 = vertex->GetTrackMom(1);
  TLorentzVector lv_pi0;
  TLorentzVector lv_pi1;
  lv_pi0.SetXYZM(p0.X(), p0.Y(), p0.Z(), pdg::PionMass());
  lv_pi1.SetXYZM(p1.X(), p1.Y(), p1.Z(), pdg::PionMass());
  const TLorentzVector lv_ks = lv_pi0 + lv_pi1;
  const Double_t m_ks = lv_ks.M();
  const Double_t m0_ks = pdg::Mass(kK0Short);
  if (m0_ks <= 0. || TMath::Abs(m_ks - m0_ks) > K0_MASS_WINDOW)
    return false;

  return true;
}

//_____________________________________________________________________________
void
TPCReconstructor::Calculate(const std::vector<TPCVertex*>& vertices,
                            const std::vector<TPCLocalTrackHelix*>* helix_tracks,
                            UInt_t reco_mode)
{
  Clear();
  const Bool_t do_lambda = (reco_mode & kRecoLambda) != 0u;
  const Bool_t do_k0 = (reco_mode & kRecoK0Short) != 0u;
  for (auto* vtx : vertices) {
    if (!vtx || !vtx->IsCalculated())
      continue;
    vtx->ClearRecoCandidates();
    if (do_lambda)
      ReconstructLambda(vtx, helix_tracks);
    if (do_k0)
      ReconstructK0Short(vtx);
  }
}
