// -*- C++ -*-

#ifndef TPC_ANALYZER_HH
#define TPC_ANALYZER_HH

#include <vector>
#include <TString.h>
#include <TVector3.h>

#include "DetectorID.hh"

class RawData;
class TPCHit;
class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrackHelix;

typedef std::vector<TPCHit*>        TPCHitContainer;
typedef std::vector<TPCCluster*>    TPCClusterContainer;
typedef std::vector<TPCLocalTrack*> TPCLocalTrackContainer;
typedef std::vector<TPCLocalTrackHelix*> TPCLocalTrackHelixContainer;

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
  { kTPC, n_type };
  std::vector<Bool_t>                m_is_decoded;
  std::vector<TPCHitContainer>       m_TPCHitCont;
  std::vector<TPCHitContainer>       m_TempTPCHitCont;
  std::vector<TPCClusterContainer>   m_TPCClCont;
  TPCLocalTrackContainer             m_TPCTC;
  TPCLocalTrackContainer             m_TPCTCFailed;
  TPCLocalTrackHelixContainer        m_TPCTCHelix;
  TPCLocalTrackHelixContainer        m_TPCTCHelixFailed;

public:
  Bool_t DecodeTPCHitsGeant4(const Int_t nhits,
                             const Double_t *x, const Double_t *y,
                             const Double_t *z, const Double_t *de);
  Bool_t DecodeTPCHits(RawData* rawData, Double_t clock=0.);
  const TPCHitContainer& GetTPCHC(Int_t l) const { return m_TPCHitCont.at(l); }
  const TPCClusterContainer& GetTPCClCont(Int_t l) const { return m_TPCClCont.at(l); }

  Bool_t TrackSearchTPC(Bool_t exclusive=false);
  Bool_t TrackSearchTPCHelix(Bool_t exclusive=false);
  Bool_t TrackSearchTPCHelix(std::vector<std::vector<TVector3>> K18VPs,
			     std::vector<std::vector<TVector3>> KuramaVPs,
			     Bool_t exclusive=false);
  Bool_t TestHoughTransform();
  Bool_t TestHoughTransformHelix();

  Int_t GetNTracksTPC() const { return m_TPCTC.size(); }
  Int_t GetNTracksTPCFailed() const { return m_TPCTCFailed.size(); }
  Int_t GetNTracksTPCHelix() const { return m_TPCTCHelix.size(); }
  Int_t GetNTracksTPCHelixFailed() const { return m_TPCTCHelixFailed.size(); }
  TPCLocalTrack* GetTrackTPC(Int_t l) const { return m_TPCTC.at(l); }
  TPCLocalTrack* GetTrackTPCFailed(Int_t l) const { return m_TPCTCFailed.at(l); }
  TPCLocalTrackHelix* GetTrackTPCHelix(Int_t l) const { return m_TPCTCHelix.at(l); }
  TPCLocalTrackHelix* GetTrackTPCHelixFailed(Int_t l) const { return m_TPCTCHelixFailed.at(l); }

  Bool_t ReCalcTPCHits(const Int_t nhits,
                       const std::vector<Int_t>& pad,
                       const std::vector<Double_t>& time,
                       const std::vector<Double_t>& de,
                       Double_t clock=0.);

  void ResetTracksTPC() { ClearTracksTPC(); }

  protected:

  void ClearTPCHits();
  void ClearTPCClusters();
  void ClearTracksTPC();
  static Bool_t MakeUpTPCClusters(const TPCHitContainer& HitCont,
                                  TPCClusterContainer& ClCont,
                                  Double_t maxdy);

};

//_____________________________________________________________________________
inline const TString&
TPCAnalyzer::ClassName()
{
  static TString s_name("TPCAnalyzer");
  return s_name;
}

#endif
