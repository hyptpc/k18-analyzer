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

  // Extrapolate helix to the target (mm).
  Bool_t ExtrapolateToTarget(const TPCLocalTrackHelix* track,
                             TVector3& pos, TVector3& mom,
                             Double_t& len, Double_t& dist) const;
  // Extrapolate helix to HTOF (mm). Candidates returned via out-vectors.
  Bool_t ExtrapolateToHTOF(const TPCLocalTrackHelix* track,
                           std::vector<Int_t>& segid,
                           std::vector<TVector3>& pos,
                           std::vector<TVector3>& mom,
                           std::vector<Double_t>& tracklen,
                           std::vector<Int_t>& plane_id,
                           std::vector<Double_t>& horizontal,
                           std::vector<Double_t>& vertical) const;

  // HTOF cluster/raw seg matching (extrap seg ↔ MeanSeg; |Δ|<=1 for cluster).
  static Int_t MatchHtofCluster(Double_t extrap_seg,
                                const std::vector<Double_t>& cl_seg);
  // Same, but among clusters with the same |Δseg| prefer the one whose time is closest to
  // time0 (never rejects a cluster by time; NaN time0 / cluster time -> seg order only).
  static Int_t MatchHtofCluster(Double_t extrap_seg,
                                const std::vector<Double_t>& cl_seg,
                                const std::vector<Double_t>& cl_time,
                                Double_t time0);
  // Closest seg index; no max-distance cut (ADC/dE attach).
  static Int_t MatchHtofBySeg(Double_t cl_seg,
                              const std::vector<Double_t>& segs);
  // Plane residual for match QA: abs_s = |(pos-origin)·n|, dRho = |ρ_xz - L|.
  Bool_t HtofMatchResidual(Int_t plane_id, const TVector3& pos,
                           Double_t& abs_s, Double_t& drho) const;
  // Center-to-plane distance [mm] of HTOF face plane_id (dist_htof_mm[]); NaN if out of range.
  static Double_t HtofPlaneDistance(Int_t plane_id);
  // Unit normal of HTOF face plane_id (TPC frame); NaN vector if out of range.
  TVector3 HtofPlaneNormal(Int_t plane_id) const;

private:

  // HTOF planes (origin + normal); distances from dist_htof_mm[].
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
