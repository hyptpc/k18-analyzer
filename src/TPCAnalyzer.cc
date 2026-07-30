// -*- C++ -*-

#include "TPCAnalyzer.hh"

#include <iostream>

#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "TPCCluster.hh"
#include "TPCHit.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCPadHelper.hh"
#include "TPCRawData.hh"
#include "TPCRawHit.hh"
#include "TPCTrackSearch.hh"
#include "TPCVertex.hh"
#include "UserParamMan.hh"

/* TPCTracking */
#define UseTpcCluster 1 // 1 : Common clustering method, 0 : Cluster size=1 no clustering

namespace
{
const auto& gUser   = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
TPCAnalyzer::TPCAnalyzer()
  : m_is_decoded(n_type),
    m_TPCHitCont(NumOfLayersTPC+1),
    m_TPCClCont(NumOfLayersTPC)
{
  for(Int_t i=0; i<n_type; ++i){
    m_is_decoded[i] = false;
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCAnalyzer::~TPCAnalyzer()
{
  ClearTPCHits();
  ClearTPCClusters();
  ClearTPCTracks();
  ClearTPCVertices();
  ClearTPCK18Tracks();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DecodeTPCHits(TPCRawData &TPCrawData,
			   const std::vector<Double_t> clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTPCHits();

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    for(const auto& rhit: TPCrawData.GetTPCCorHitContainer(layer)){
      auto hit = new TPCHit(rhit);
      Int_t row  = hit -> GetRow();
      Int_t cobo_id = tpc::GetCoBoId(layer, row);
      if(cobo_id >= clock.size()){
        delete hit;
        return false; //No cobo input
      }

      if(hit->DoFit() && hit->Calculate(clock[cobo_id])){
        m_TPCHitCont[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::MakeUpTPCClusters(const TPCHitContainer& hit_cont,
			       TPCClusterContainer& cl_cont,
			       Double_t max_dy)
{
  static const Double_t min_cluster_de   = gUser.GetParameter("MinClusterDeTPC");
  static const Int_t    min_cluster_size = gUser.GetParameter("MinClusterSizeTPC");
  static const Double_t min_cluster_ypos = gUser.GetParameter("MinClusterYPosTPC");
  static const Double_t max_cluster_ypos = gUser.GetParameter("MaxClusterYPosTPC");

  const auto nh = hit_cont.size();
  if(nh==0) return false;

  std::vector<Int_t> joined(nh, 0);
  for(Int_t i=0; i<nh; ++i){
    if(joined[i] > 0) continue;
    TPCHitContainer cand_cont;
    TPCHit* hit = hit_cont[i];
    if(!hit || !hit->IsGood()) continue;
    Int_t layer = hit->GetLayer();
    cand_cont.push_back(hit);
    joined[i]++;
    Double_t padlength = hit->GetPadLength();
    TVector3 dist2tgt = hit->GetPosition() - TVector3(0., 0., tpc::Z_TARGET);
    Double_t vertical_pathlength =
      padlength*dist2tgt.y()/TMath::Hypot(dist2tgt.x(), dist2tgt.z());
    max_dy = TMath::Max(max_dy, vertical_pathlength);
#if UseTpcCluster
    for(Int_t j=0; j<nh; ++j){
      if(i==j || joined[j]>0) continue;
      TPCHit* thit = hit_cont[j];
      if(!thit || !thit->IsGood()) continue;
      Int_t row_id = thit->GetRow();
      for(const auto& c_hit: cand_cont){
        Int_t c_row_id = c_hit->GetRow();
        if(tpc::IsClusterable(layer, row_id, c_row_id)
           && TMath::Abs(thit->GetY() - c_hit->GetY()) < max_dy){
          cand_cont.push_back(thit);
          joined[j]++;
          break;
        }
      }
    }
#endif

    TPCCluster* cluster = new TPCCluster(layer, cand_cont);
    if(!cluster) continue;
    if(cluster->Calculate()
       && cluster->GetDe()>=min_cluster_de 
       && cluster->GetClusterSize()>=min_cluster_size
       && cluster->GetY()>=min_cluster_ypos 
       && cluster->GetY()<=max_cluster_ypos)
    {
      cl_cont.push_back(cluster);
    }else{
      delete cluster;
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::ReCalcTPCHits(const Int_t nhits,
			   const std::vector<Int_t>& pad,
			   const std::vector<Double_t>& time,
			   const std::vector<Double_t>& de,
			   const std::vector<Double_t>& clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  // MinCDeTPC: minimum corrected dE for a hit before clustering
  static const Double_t min_cde = gUser.GetParameter("MinCDeTPC");

  ClearTPCHits();
  ClearTPCClusters();

  if(nhits != pad.size() || nhits != time.size() || nhits != de.size()){
    hddaq::cerr << FUNC_NAME << " vector size mismatch" << std::endl;
    return false;
  }

  static const Bool_t NoiseOff = gUser.Has("NoiseOffClusterTPC")
                              ? (gUser.GetParameter("NoiseOffClusterTPC") == 1.)
                              : false;
  
  for(Int_t ih=0; ih<nhits; ih++){
    const Int_t layer  = tpc::getLayerID(pad[ih]);
    const Double_t row = tpc::getRowID(pad[ih]);
    auto hit = new TPCHit(layer, row);
    hit->AddHit(de[ih], time[ih]);
    Int_t cobo_id = tpc::GetCoBoId(layer, row);
    if(cobo_id >= clock.size()){
      delete hit;
      return false; //No cobo input
    }

    // Remove Noise Pad (padAbnormalWaveform_E72)
    Bool_t noise_pad = false;
    if(NoiseOff){
      if(tpc::Noise(layer,row)){
      	noise_pad = true;
      }
    }
        
    if(hit->Calculate(clock[cobo_id]) && hit->GetCDe()>=min_cde && hit->IsGood() && !noise_pad){
      m_TPCHitCont[layer].push_back(hit);
    }else{
      delete hit;
    }
  }
  

#if 1
  // MaxYDifClusterTPC: max |dY| between hits merged into one cluster (per layer)
  static const Double_t max_y_dif = gUser.GetParameter("MaxYDifClusterTPC");
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], max_y_dif);
  }
#endif

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
// Geant4 version of ReCalcTPCHits:
//   - Y position is set directly from ytpc_pad (no clock-based drift calculation)
//   - dE is used as-is (no gain calibration assumed for Geant4)
Bool_t
TPCAnalyzer::ReCalcTPCHitsGeant4(const Int_t nhits,
                                  const std::vector<Int_t>& pad,
                                  const std::vector<Double_t>& de,
                                  const std::vector<Double_t>& ytpc_pad)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  static const Double_t min_cde = gUser.GetParameter("MinCDeTPC");

  ClearTPCHits();
  ClearTPCClusters();

  if(nhits != (Int_t)pad.size() || nhits != (Int_t)de.size() || nhits != (Int_t)ytpc_pad.size()){
    hddaq::cerr << FUNC_NAME << " vector size mismatch" << std::endl;
    return false;
  }

  for(Int_t ih = 0; ih < nhits; ih++){
    const Int_t layer = tpc::getLayerID(pad[ih]);
    const Int_t row   = tpc::getRowID(pad[ih]);
    auto hit = new TPCHit(layer, row);
    hit->AddHit(de[ih], 0.);      // time=0 (unused); appends to m_de/m_cde/m_position
    hit->SetDe(de[ih]);           // explicitly set m_cde[0] = de[ih]
    TVector3 pad_pos = tpc::GetPosition(pad[ih]);  // (x, 0, z) from geometry
    pad_pos.SetY(ytpc_pad[ih]);   // override Y with Geant4 drift position
    hit->SetPosition(pad_pos);

    if(hit->GetCDe() >= min_cde && hit->IsGood()){
      m_TPCHitCont[layer].push_back(hit);
    }else{
      delete hit;
    }
  }

  static const Double_t max_y_dif = gUser.GetParameter("MaxYDifClusterTPC");
  for(Int_t layer = 0; layer < NumOfLayersTPC; ++layer){
    MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], max_y_dif);
  }

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
// Geant4 cluster-level input: keep detector-response cluster positions intact
// and make a single-hit TPCCluster for each entry, so that the common helix
// tracker starts from m_TPCClCont without re-clustering.
Bool_t
TPCAnalyzer::ReCalcTPCHitsGeant4(const std::vector<Int_t>& pad,
                                 const std::vector<Double_t>& de,
                                 const std::vector<Double_t>& x,
                                 const std::vector<Double_t>& y,
                                 const std::vector<Double_t>& z)
{
  if (m_is_decoded[kTPC]) {
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }
  const auto nhits = pad.size();
  if (de.size() != nhits || x.size() != nhits || y.size() != nhits || z.size() != nhits) {
    hddaq::cerr << FUNC_NAME << " vector size mismatch" << std::endl;
    return false;
  }

  ClearTPCHits();
  ClearTPCClusters();
  for (std::size_t ih = 0; ih < nhits; ++ih) {
    // Geant4 pad IDs are zero-based.  Ignore only rare edge entries that are
    // outside the active analyzer geometry before getLayerID can throw.
    if (pad[ih] < 0 || pad[ih] >= NumOfPadTPC) continue;
    const Int_t layer = tpc::getLayerID(pad[ih]);
    const Int_t row = tpc::getRowID(pad[ih]);
    if (layer < 0 || layer >= NumOfLayersTPC) continue;

    auto* hit = new TPCHit(layer, row);
    hit->AddHit(de[ih], 0.);
    hit->SetDe(de[ih]);
    hit->SetPad(pad[ih]);
    hit->SetPosition(TVector3(x[ih], y[ih], z[ih]));
    hit->SetIsGood(true);
    m_TPCHitCont[layer].push_back(hit);

    TPCHitContainer cluster_hits{hit};
    auto* cluster = new TPCCluster(layer, cluster_hits);
    if (cluster->Calculate()) {
      cluster->SetClusterSizeG4(1);
      m_TPCClCont[layer].push_back(cluster);
    } else {
      delete cluster;
    }
  }
  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
//HS-OFF: Track searching
Bool_t
TPCAnalyzer::TrackSearchTPC(Bool_t exclusive)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  // MinLayerTPC: min cluster count for Hough/fit (not number of TPC layers)
  static const Int_t min_num_of_hits = gUser.GetParameter("MinLayerTPC");

  tpc::LocalTrackSearch(m_TPCClCont, m_TPCTC, m_TPCTCFailed, exclusive, min_num_of_hits);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
//HS-On: Track Searching
Bool_t
TPCAnalyzer::TrackSearchTPCHelix(Bool_t exclusive, UInt_t reco_mode)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  // MinLayerTPC: min cluster count for Hough/fit (not number of TPC layers)
  static const Int_t min_num_of_hits = gUser.GetParameter("MinLayerTPC");

  tpc::LocalTrackSearchHelix(m_TPCClCont, m_TPCTCHelix, m_TPCTCHelixInverted, m_TPCTCHelixFailed, m_TPCVC, m_TPCVCClustered, exclusive, min_num_of_hits);

  if (reco_mode != TPCReconstructor::kRecoNone) {
    TPCReconstructor reconstructor;
    reconstructor.Calculate(m_TPCVC, &m_TPCTCHelix, reco_mode);
  }

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
//HS-On: Track Searching with BcOut Track
Bool_t
TPCAnalyzer::TrackSearchTPCHelix(std::vector<std::vector<TVector3>> K18VPs,
				 Bool_t exclusive)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  // MinLayerTPC: min cluster count for Hough/fit (not number of TPC layers)
  static const Int_t min_num_of_hits = gUser.GetParameter("MinLayerTPC");
  
  tpc::LocalTrackSearchHelix(K18VPs, m_TPCClCont, m_TPCTCHelix, m_TPCTCHelixInverted, m_TPCTCVP, m_TPCTCHelixFailed, m_TPCVC, m_TPCVCClustered, exclusive, min_num_of_hits);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
TPCVertex*
TPCAnalyzer::FindVertexTPC(Int_t id1, Int_t id2) const
{
  TPCVertex* found_vertex = nullptr;
  Int_t n_matches = 0;
  for (const auto& vertex : m_TPCVC) {
    if (!vertex || vertex->GetNTracks() < 2)
      continue;
    const Int_t track_id1 = vertex->GetTrackId(0);
    const Int_t track_id2 = vertex->GetTrackId(1);
    if ((track_id1 == id1 && track_id2 == id2) ||
        (track_id1 == id2 && track_id2 == id1)) {
      if (!found_vertex)
        found_vertex = vertex;
      ++n_matches;
    }
  }
  if (n_matches > 1) {
    std::cerr << "#W " << FUNC_NAME
              << " duplicate vertices for pair=(" << id1 << "," << id2 << ")"
              << " n_matches=" << n_matches
              << " (return first match)" << std::endl;
  }
  return found_vertex;
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCHits()
{
  del::ClearContainerAll(m_TPCHitCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCClusters()
{
  del::ClearContainerAll(m_TPCClCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCTracks()
{
  del::ClearContainer(m_TPCTC);
  del::ClearContainer(m_TPCTCFailed);
  del::ClearContainer(m_TPCTCHelix);
  del::ClearContainer(m_TPCTCHelixInverted);
  del::ClearContainer(m_TPCTCHelixFailed);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCVertices()
{
  del::ClearContainer(m_TPCVC);
  del::ClearContainer(m_TPCVCClustered);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCK18Tracks()
{
  del::ClearContainer(m_TPCK18TC);
}
