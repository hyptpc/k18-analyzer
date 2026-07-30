// -*- C++ -*-

#include "TPCReconstructor.hh"

#include <iomanip>
#include <iostream>

#include <TLorentzVector.h>
#include <TPDGCode.h>

#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "Kinematics.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "UserParamMan.hh"

#define TPC_RECO_ENABLE_LAMBDA_REFIT 0
#define TPC_RECO_PRINT_CANDIDATES 0

namespace
{
  const Double_t MIN_CLOSE_DIST_LAMBDA = 10.; // mm, same scale as TPCVertex / TPCTrackSearch
  const Double_t MIN_CLOSE_DIST_K0 = 10.; // mm, same scale as TPCVertex / TPCTrackSearch
  const Double_t LAMBDA_MASS_WINDOW = 0.3; // GeV/c^2
  const Double_t K0_MASS_WINDOW = 0.3; // GeV/c^2
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
    if (vertex->GetClosestDist() > MIN_CLOSE_DIST_LAMBDA)
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
                                   Int_t track1_id, Int_t track2_id,
                                   Int_t track1_pid, Int_t track2_pid,
                                   Int_t track1_charge, Int_t track2_charge)
  : m_mother_pdg(mother_pdg),
    m_mass(mass),
    m_vertex(vertex),
    m_momentum(momentum),
    m_closest_dist(closest_dist),
    m_track1_id(track1_id),
    m_track2_id(track2_id),
    m_track1_pid(track1_pid),
    m_track2_pid(track2_pid),
    m_track1_charge(track1_charge),
    m_track2_charge(track2_charge)
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
            << "   trkId = (" << m_track1_id
            << (m_track1_charge > 0 ? "+" : (m_track1_charge < 0 ? "-" : "0")) << ","
            << m_track2_id
            << (m_track2_charge > 0 ? "+" : (m_track2_charge < 0 ? "-" : "0")) << ")\n"
            << "   pid   = (" << m_track1_pid << "," << m_track2_pid << ")\n"
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
#if TPC_RECO_PRINT_CANDIDATES
    cand.Print("Lambda");
#endif
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

  const Int_t track1_id = vertex->GetTrackId(0);
  const Int_t track2_id = vertex->GetTrackId(1);
  if (track1_id < 0 || track2_id < 0 ||
      track1_id >= static_cast<Int_t>(helix_tracks->size()) ||
      track2_id >= static_cast<Int_t>(helix_tracks->size()))
    return false;

  TPCLocalTrackHelix* base_track1 = helix_tracks->at(track1_id);
  TPCLocalTrackHelix* base_track2 = helix_tracks->at(track2_id);
  if (!base_track1 || !base_track2)
    return false;

  auto* track1 = new TPCLocalTrackHelix(base_track1);
  if (!track1->DoFitTrackwVertex(vertex->GetVertex(), VERTEX_RES)) {
    delete track1;
    return false;
  }
  auto* track2 = new TPCLocalTrackHelix(base_track2);
  if (!track2->DoFitTrackwVertex(vertex->GetVertex(), VERTEX_RES)) {
    delete track1;
    delete track2;
    return false;
  }

  const Double_t vertex_scan_range = UserParamMan::GetInstance().GetParameter("VertexScanRange");
  Double_t par1[5];
  Double_t par2[5];
  track1->GetParam(par1);
  track2->GetParam(par2);
  const Double_t scan_theta1 = vertex_scan_range/par1[kHelixR];
  const Double_t scan_theta2 = vertex_scan_range/par2[kHelixR];
  const Double_t range_theta1[2] = {track1->GetMint() - scan_theta1,
                                    track1->GetMaxt() + scan_theta1};
  const Double_t range_theta2[2] = {track2->GetMint() - scan_theta2,
                                    track2->GetMaxt() + scan_theta2};

  Double_t theta1 = 0.0;
  Double_t theta2 = 0.0;
  Double_t dist = 0.0;
  const TVector3 vertex_refit = Kinematics::VertexPointHelix(
    par1, par2,
    range_theta1[0], range_theta1[1],
    range_theta2[0], range_theta2[1],
    theta1, theta2, dist
  );
  if (!std::isfinite(dist) || dist > MIN_CLOSE_DIST_LAMBDA) {
    delete track1;
    delete track2;
    return false;
  }

  out_vertex = vertex_refit;
  out_dist = dist;
  out_mom1 = track1->CalcHelixMom(par1, theta1);
  out_mom2 = track2->CalcHelixMom(par2, theta2);
  delete track1;
  delete track2;
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
  if (dist > MIN_CLOSE_DIST_K0)
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
#if TPC_RECO_PRINT_CANDIDATES
  cand.Print("K0Short");
#endif
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
  if (dist > MIN_CLOSE_DIST_K0)
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
