// -*- C++ -*-

#ifndef TPC_ANALYZER_HH
#define TPC_ANALYZER_HH

#include <vector>

#include <TString.h>
#include <TVector3.h>

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
  // Geant4 cluster-level input.  One input element is kept as one cluster;
  // do not run the pad-hit clustering step.
  Bool_t ReCalcTPCHitsGeant4(const std::vector<Int_t>& pad,
                             const std::vector<Double_t>& de,
                             const std::vector<Double_t>& x,
                             const std::vector<Double_t>& y,
                             const std::vector<Double_t>& z);

  //HS-Off
  Bool_t TrackSearchTPC(Bool_t exclusive=false);
  Int_t GetNTracksTPC() const { return m_TPCTC.size(); }
  Int_t GetNTracksTPCFailed() const { return m_TPCTCFailed.size(); }
  TPCLocalTrack* GetTrackTPC(Int_t l) const { return m_TPCTC.at(l); }
  TPCLocalTrack* GetTrackTPCFailed(Int_t l) const { return m_TPCTCFailed.at(l); }

  //HS-On
  Bool_t TrackSearchTPCHelix(Bool_t exclusive=false,
                             UInt_t reco_mode=TPCReconstructor::kRecoAll);
  Bool_t TrackSearchTPCHelix(std::vector<std::vector<TVector3>> K18VPs,
			     Bool_t exclusive=false);
  Int_t GetNTracksTPCHelix() const { return m_TPCTCHelix.size(); }
  TPCLocalTrackHelix* GetTrackTPCHelix(Int_t l) const { return m_TPCTCHelix.at(l); }
  Int_t GetNTracksTPCHelixVP() const { return m_TPCTCVP.size(); }
  TPCLocalTrackHelix* GetTrackTPCHelixVP(Int_t l) const { return m_TPCTCVP.at(l); }
  
  Int_t GetNVerticesTPC() const { return m_TPCVC.size(); }
  TPCVertex* GetVertexTPC(Int_t i) const { return m_TPCVC.at(i); }
  TPCVertex* FindVertexTPC(Int_t id1, Int_t id2) const;

  const TPCHitContainer& GetTPCHC(Int_t l) const { return m_TPCHitCont.at(l); }
  const TPCClusterContainer& GetTPCClCont(Int_t l) const { return m_TPCClCont.at(l); }

  // Extrapolate helix to HTOF / target (mm).
  Bool_t ExtrapolateToTarget(const TPCLocalTrackHelix* track,
                             TVector3& pos, TVector3& mom,
                             Double_t& len, Double_t& dist) const;
  Bool_t ExtrapolateToHTOF(const TPCLocalTrackHelix* track,
                           std::vector<Int_t>& segid,
                           std::vector<TVector3>& pos,
                           std::vector<TVector3>& mom,
                           std::vector<Double_t>& tracklen,
                           std::vector<Int_t>& plane_id,
                           std::vector<Double_t>& horizontal,
                           std::vector<Double_t>& vertical) const;

private:

  // HTOF plane geometry (see dist_htof_mm in .cc).
  TVector3 m_htof_origin[NumOfPlanesHTOF];
  TVector3 m_htof_normal[NumOfPlanesHTOF];

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
