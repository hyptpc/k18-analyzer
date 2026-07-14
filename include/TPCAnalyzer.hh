// -*- C++ -*-

#ifndef TPC_ANALYZER_HH
#define TPC_ANALYZER_HH

#include <vector>

#include <TString.h>

#include "DetectorID.hh"
#include "TPCReconstructor.hh"

class TPCRawData;
class TPCHit;
class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrackHelix;
class TPCVertex;

typedef std::vector<TPCHit*>             TPCHitContainer;
typedef std::vector<TPCCluster*>         TPCClusterContainer;
typedef std::vector<TPCLocalTrack*>      TPCLocalTrackContainer;
typedef std::vector<TPCLocalTrackHelix*> TPCLocalTrackHelixContainer;
typedef std::vector<TPCVertex*>          TPCVertexContainer;

//_____________________________________________________________________________
class TPCAnalyzer
{
public:
  static const TString& ClassName();
  TPCAnalyzer();
  ~TPCAnalyzer();

private:
  TPCAnalyzer(const TPCAnalyzer&);
  TPCAnalyzer& operator =(const TPCAnalyzer&);

private:
  enum e_type
  { kTPC, kTPCTracking, kTPCK18, n_type };
  std::vector<Bool_t>                m_is_decoded;
  std::vector<TPCHitContainer>       m_TPCHitCont;
  std::vector<TPCClusterContainer>   m_TPCClCont;
  TPCLocalTrackContainer             m_TPCTC;
  TPCLocalTrackContainer             m_TPCTCFailed;
  TPCLocalTrackHelixContainer        m_TPCTCHelix;
  TPCLocalTrackHelixContainer        m_TPCTCHelixInverted;
  TPCLocalTrackHelixContainer        m_TPCTCHelixFailed;
  TPCLocalTrackHelixContainer        m_TPCTCVP;
  TPCVertexContainer                 m_TPCVC; //vertex between two tracks
  TPCVertexContainer                 m_TPCVCClustered; //clusted position of multi-tracks
  TPCLocalTrackContainer             m_TPCK18TC;
  


public:

  //TPC Hit&Cluster
  Bool_t DecodeTPCHits(TPCRawData &TPCrawData, const std::vector<Double_t> clock);
  Bool_t ReCalcTPCHits(const Int_t nhits,
                       const std::vector<Int_t>& pad,
                       const std::vector<Double_t>& time,
                       const std::vector<Double_t>& de,
                       const std::vector<Double_t>& clock);
  // Geant4 version: Y position set directly from ytpc_pad (no clock-based drift)
  Bool_t ReCalcTPCHitsGeant4(const Int_t nhits,
                             const std::vector<Int_t>& pad,
                             const std::vector<Double_t>& de,
                             const std::vector<Double_t>& ytpc_pad);

  //HS-Off
  Bool_t TrackSearchTPC(Bool_t exclusive=false);
  Int_t GetNTracksTPC() const { return m_TPCTC.size(); }
  Int_t GetNTracksTPCFailed() const { return m_TPCTCFailed.size(); }
  TPCLocalTrack* GetTrackTPC(Int_t l) const { return m_TPCTC.at(l); }
  TPCLocalTrack* GetTrackTPCFailed(Int_t l) const { return m_TPCTCFailed.at(l); }

  //HS-On
  Bool_t TrackSearchTPCHelix(Bool_t exclusive=false,
                             UInt_t reco_mode=TPCReconstructor::kRecoAll);
  Bool_t TrackSearchTPCHelix(std::vector<std::vector<TVector3>> K18BRVPs,
			     Bool_t exclusive=false);
  Int_t GetNTracksTPCHelix() const { return m_TPCTCHelix.size(); }
  TPCLocalTrackHelix* GetTrackTPCHelix(Int_t l) const { return m_TPCTCHelix.at(l); }
  
  Int_t GetNVerticesTPC() const { return m_TPCVC.size(); }
  TPCVertex* GetVertexTPC(Int_t i) const { return m_TPCVC.at(i); }
  TPCVertex* FindVertexTPC(Int_t id1, Int_t id2) const;

  const TPCHitContainer& GetTPCHC(Int_t l) const { return m_TPCHitCont.at(l); }
  const TPCClusterContainer& GetTPCClCont(Int_t l) const { return m_TPCClCont.at(l); }

private:

protected:

  void ClearTPCHits();
  void ClearTPCClusters();
  void ClearTPCTracks();
  void ClearTPCVertices();
  void ClearTPCK18Tracks();
  static Bool_t MakeUpTPCClusters(const TPCHitContainer& hit_cont,
                                  TPCClusterContainer& cl_cont,
                                  Double_t max_dy);

};

//_____________________________________________________________________________
inline const TString&
TPCAnalyzer::ClassName()
{
  static TString s_name("TPCAnalyzer");
  return s_name;
}

#endif
