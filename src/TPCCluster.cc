 // -*- C++ -*-

#include "TPCCluster.hh"

#include <iomanip>
#include <iostream>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"
#include "ThreeVector.hh"
#include "TPCPadHelper.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#include <escape_sequence.hh>
#include <spdlog/spdlog.h>
#include <std_ostream.hh>

#define TPC_CLUSTER_WRAP_DEBUG 0

namespace
{
const auto& gTPCPos = TPCPositionCorrector::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
TPCCluster::TPCCluster(Int_t layer, const TPCHitContainer& HitCont)
  : m_is_good(false),
    m_is_onframe(false),
    m_layer(layer),
    m_cluster_de(),
    m_cluster_position(),
    m_hit_array(HitCont), // shallow copy
    m_mean_row(),
    m_mean_theta(),
    m_center_hitid(-1),
    m_mean_hit(new TPCHit(layer, TMath::QuietNaN()))
{
  auto itr = m_hit_array.begin();
  while(itr != m_hit_array.end()){
    if(!*itr || !(*itr)->IsGood()){
      itr = m_hit_array.erase(itr);
    }else{
      itr++;
    }
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCCluster::~TPCCluster()
{
  ClearTPCHits();
  delete m_mean_hit;
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCCluster::ClearTPCHits()
{
  m_hit_array.clear();
}

//_____________________________________________________________________________
void
TPCCluster::AddTPCHit(TPCHit* hit)
{
  if(hit) m_hit_array.push_back(hit);
  m_is_good = false;
  m_is_onframe = false;
}

//_____________________________________________________________________________
void
TPCCluster::CheckClusterOnTheFrame()
{
  m_is_onframe = false;
  if(m_layer<8) return; // no frame

  Int_t low_row = 10000; Int_t high_row = -1; //get edge pad ids of the cluster
  for(Int_t i=0; i<m_hit_array.size(); ++i){
    if(!m_hit_array[i]) continue;
    Int_t row = m_hit_array[i] -> GetRow();
    if(row<low_row) low_row = row;
    if(row>high_row) high_row = row;
  }

  Bool_t status = false;
  // FrameHighEdge and FrameLowEdge are [NumOfLayersTPC][5] arrays
  // The second dimension size is 5 (maximum number of frame edges per layer)
  static const Int_t NumFrameEdges = 5;
  for(Int_t i=0; i<NumFrameEdges; ++i){
    if(TMath::Abs(tpc::frameHighEdge[m_layer][i] - low_row) <= tpc::MAX_ROW_DIF_TPC) status = true;
    if(TMath::Abs(tpc::frameLowEdge[m_layer][i] - high_row) <= tpc::MAX_ROW_DIF_TPC) status = true;
  }

  m_is_onframe = status;
}

//_____________________________________________________________________________
Bool_t
TPCCluster::Calculate()
{
  static const TVector2 target_center(0., tpc::Z_TARGET); // (X, Z)
  Int_t max_row = static_cast<Int_t>(tpc::padParameter[m_layer][tpc::kNumOfPad]);
  m_cluster_de = 0.;
  m_cluster_position.SetXYZ(0., 0., 0.);

  // w/o position correction
  TVector2 xz_vector_hs0(0., 0.);

  // w/ position correction
  m_mean_row = 0.;
  m_mean_theta = 0.;
  Double_t mean_y = 0.;
  TVector2 xz_vector_hs(0., 0.);
  for(const auto& hit: m_hit_array){
    const auto& pos = hit->GetPosition();
    TVector2 xz_vector(pos.X(), pos.Z());
    xz_vector -= target_center;

    const Double_t cde = hit->GetCDe();
    mean_y += pos.Y() * cde;
    m_cluster_de += cde;
    xz_vector_hs += cde * xz_vector;

    Int_t row = hit->GetRow();
    Int_t pad_id = tpc::GetPadId(m_layer, row);
    auto pos0 = tpc::GetPosition(pad_id);
    TVector2 xz_vector0(pos0.X(), pos0.Z());
    xz_vector0 -= target_center;
    xz_vector_hs0 += cde * xz_vector0;
  }

  mean_y /= m_cluster_de;
  xz_vector_hs /= m_cluster_de;
  m_mean_theta = xz_vector_hs.Phi();
  TVector2 xz_vector = xz_vector_hs + target_center;
  m_cluster_position.SetXYZ(xz_vector.X(), mean_y, xz_vector.Y());
  m_mean_row = tpc::GetMrow(m_layer, m_mean_theta*TMath::RadToDeg());

  const Double_t raw_mean_row = m_mean_row;
  const Bool_t is_inner_layer = tpc::IsCircularLayer(m_layer);
  Double_t resolved_mean_row = TMath::QuietNaN();
  if (!tpc::TryResolveMRow(m_layer, raw_mean_row, resolved_mean_row)) {
#if TPC_CLUSTER_WRAP_DEBUG
    const Double_t n_div = tpc::padParameter[m_layer][tpc::kNumOfDivision];
    std::cout << "\n"
              << "[TPCCluster::WrapDebug][reject m_mean_row]\n"
              << "   layer       = " << m_layer << "\n"
              << "   n_pad       = " << max_row << "\n"
              << "   n_div       = " << n_div << "\n"
              << "   phi_deg     = " << m_mean_theta*TMath::RadToDeg() << "\n"
              << "   m_row_raw   = " << raw_mean_row << "\n"
              << "   weighted_xz = (" << xz_vector_hs.X() << ", " << xz_vector_hs.Y()
              << ")\n"
              << "   cluster_de  = " << m_cluster_de << "\n"
              << "   size        = " << m_hit_array.size()
              << std::endl;
    Int_t printed = 0;
    for (const auto& hit : m_hit_array) {
      if (!hit) continue;
      const auto& pos = hit->GetPosition();
      const Int_t hit_padid = tpc::GetPadId(m_layer, hit->GetRow());
      const TVector3 nominal_pos = tpc::GetPosition(hit_padid);
      const Double_t hit_phi = (TVector2(pos.X(), pos.Z()) - target_center).Phi()*TMath::RadToDeg();
      const Double_t nominal_phi = (TVector2(nominal_pos.X(), nominal_pos.Z()) - target_center).Phi()*TMath::RadToDeg();
      std::cout << "   hit[" << printed << "]\n"
                << "      row             = " << hit->GetRow() << "\n"
                << "      cde             = " << hit->GetCDe() << "\n"
                << "      pos             = (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ")\n"
                << "      phi_deg_corr    = " << hit_phi << "\n"
                << "      phi_deg_nominal = " << nominal_phi
                << std::endl;
      ++printed;
    }
#endif
    spdlog::warn(
      "[TPCCluster::Calculate] invalid m_mean_row={} for layer {}. Rejecting cluster.",
      raw_mean_row, m_layer);
    return false;
  }

  Int_t row_id = -1;
  if (!tpc::TryResolveRow(m_layer, raw_mean_row, row_id)) {
    spdlog::warn(
      "[TPCCluster::Calculate] failed to resolve row from m_mean_row={} for layer {}. Rejecting cluster.",
      raw_mean_row, m_layer);
    return false;
  }

  m_mean_row = resolved_mean_row;
  m_mean_hit->SetPad(tpc::GetPadId(m_layer, row_id));
  m_mean_hit->AddHit(0., 0.);
  m_mean_hit->SetMRow(m_mean_row);
  m_mean_hit->SetPadLength(tpc::padParameter[m_layer][tpc::kLength]);
  m_mean_hit->SetPadTheta(tpc::GetTheta(m_layer, m_mean_row)*TMath::DegToRad());
  m_mean_hit->SetDe(m_cluster_de);
  m_mean_hit->SetPosition(m_cluster_position);
  m_mean_hit->SetParentCluster(this);
  m_mean_hit->SetIsGood(true);

  // Nearest cluster hit to mean_row0 (nominal-weighted phi). 
  // If none within MaxCenterRowDiffTPC, m_center_hitid stays -1; GetCenterHit() then uses m_mean_hit.
  Double_t mean_phi0 = xz_vector_hs0.Phi();
  Double_t mean_row0 = tpc::GetMrow(m_layer, mean_phi0*TMath::RadToDeg());
  if (!tpc::TryResolveMRow(m_layer, mean_row0, mean_row0)) {
#if TPC_CLUSTER_WRAP_DEBUG
    std::cout << "\n"
              << "[TPCCluster::WrapDebug][reject center mean_row]\n"
              << "   layer         = " << m_layer << "\n"
              << "   n_pad         = " << max_row << "\n"
              << "   mean_phi0_deg = " << mean_phi0*TMath::RadToDeg() << "\n"
              << "   mean_row0_raw = " << tpc::GetMrow(m_layer, mean_phi0*TMath::RadToDeg()) << "\n"
              << "   size          = " << m_hit_array.size()
              << std::endl;
#endif
    spdlog::warn(
      "[TPCCluster::Calculate] invalid center mean_row={} for layer {}. Rejecting cluster.",
      mean_row0, m_layer);
    return false;
  }

  // MaxCenterRowDiffTPC: max |row - mean_row| to select center hit (default 10 if absent)
  static const Double_t max_center_row_diff =
    gUser.Has("MaxCenterRowDiffTPC")
      ? gUser.GetParameter("MaxCenterRowDiffTPC")
      : 10.;
  Int_t center_hit_id = -1;
  Double_t best_row_diff = max_center_row_diff;
  for (Int_t i = 0; i < static_cast<Int_t>(m_hit_array.size()); ++i) {
    if (!m_hit_array[i]) continue;
    const Double_t row = static_cast<Double_t>(m_hit_array[i]->GetRow());
    Double_t candidate_row_diff = std::abs(mean_row0 - row);
    if (is_inner_layer) {
      candidate_row_diff = std::min(candidate_row_diff, max_row - candidate_row_diff);
    }
    if (candidate_row_diff < best_row_diff) {
      best_row_diff = candidate_row_diff;
      center_hit_id = i;
    }
  }
  m_center_hitid = center_hit_id;
  CheckClusterOnTheFrame(); //check whether the cluster on the frame or not
  m_is_good = true;
  return true;
}

//_____________________________________________________________________________
// True if a real cluster hit was chosen as the nominal center pad; 
// false when m_center_hitid is -1 (too far from mean_row0).
Bool_t
TPCCluster::HasCenterHit() const
{
  return m_center_hitid >= 0
         && m_center_hitid < static_cast<Int_t>(m_hit_array.size())
         && m_hit_array[m_center_hitid];
}

//_____________________________________________________________________________
// Center pad hit for Dst. If none was selected, returns m_mean_hit so callers still get a valid TPCHit*.
TPCHit*
TPCCluster::GetCenterHit() const
{
  if (HasCenterHit())
    return m_hit_array[m_center_hitid];
  return m_mean_hit;
}

//_____________________________________________________________________________
void
TPCCluster::Print(Option_t* opt) const
{
  PrintHelper helper(1, std::ios::fixed);
  const Double_t R = tpc::GetRadius(m_layer);
  hddaq::cout << FUNC_NAME << " " << std::endl
              << "is good = "  << m_is_good  << std::endl
              << "de = "  << m_cluster_de << std::endl
              << "position = " << m_cluster_position << std::endl
              << "Radius = " << R << std::endl
              << "mean row = " << m_mean_row << std::endl
              << "mean phi = " << m_mean_theta*TMath::RadToDeg() << " (in XZ plane)" << std::endl;
  hddaq::cout << "layer" << std::setw(2) << m_layer <<" size="
              << std::setw(3) << m_hit_array.size() << "  ";
  for(const auto& hit: m_hit_array){
    hddaq::cout << "(de=" << hddaq::unpacker::esc::k_purple
                << hit->GetCDe() << hddaq::unpacker::esc::k_default_color
                << ", y=" << hddaq::unpacker::esc::k_cyan
                << hit->GetDriftLength() << hddaq::unpacker::esc::k_default_color
                << ")" << " ";
    const auto& pos = hit->GetPosition();
    TVector2 xz_vector(pos.X(), pos.Z());
    TVector2 target_position(0., tpc::Z_TARGET);
    auto phi = (xz_vector - target_position).Phi()*TMath::RadToDeg();
    auto residual = tpc::ArcLength(m_layer, hit->GetRow(), m_mean_row);
    hddaq::cout << " pos=" << pos << ", phi=" << phi << ", res=" << residual << std::endl;
  }
  hddaq::cout << std::endl;
}
