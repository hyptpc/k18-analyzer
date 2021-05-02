// -*- C++ -*-

#ifndef DC_ANALYZER_HH
#define DC_ANALYZER_HH

#include <vector>
#include <TString.h>

#include "DetectorID.hh"
#include "ThreeVector.hh"

class DCHit;
class DCLocalTrack;
class K18TrackU2D;
class K18TrackD2U;
class KuramaTrack;
class RawData;
class MWPCCluster;
class FiberCluster;
class HodoCluster;

class TPCHit;
class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrack_Helix;

class Hodo1Hit;
class Hodo2Hit;
class HodoAnalyzer;

typedef std::vector<DCHit*>        DCHitContainer;
typedef std::vector<MWPCCluster*>  MWPCClusterContainer;
typedef std::vector<DCLocalTrack*> DCLocalTrackContainer;
typedef std::vector<K18TrackU2D*>  K18TrackU2DContainer;
typedef std::vector<K18TrackD2U*>  K18TrackD2UContainer;
typedef std::vector<KuramaTrack*>  KuramaTrackContainer;

typedef std::vector<TPCHit*>        TPCHitContainer;
typedef std::vector<TPCCluster*>    TPCClusterContainer;
typedef std::vector<TPCLocalTrack*> TPCLocalTrackContainer;
typedef std::vector<TPCLocalTrack_Helix*> TPCLocalTrack_HelixContainer;

typedef std::vector<Hodo1Hit*> Hodo1HitContainer;
typedef std::vector<Hodo2Hit*> Hodo2HitContainer;
typedef std::vector<HodoCluster*> HodoClusterContainer;

//_____________________________________________________________________________
class DCAnalyzer
{
public:
  static TString& ClassName();
  DCAnalyzer();
  ~DCAnalyzer();

private:
  DCAnalyzer(const DCAnalyzer&);
  DCAnalyzer& operator =(const DCAnalyzer&);

private:
  enum e_type { k_BcIn,  k_BcOut,
		k_SdcIn, k_SdcOut,
		k_SsdIn, k_SsdOut,
		k_TPC,
		k_TOF, n_type };
  std::vector<bool>     m_is_decoded;
  std::vector<int>      m_much_combi;
  std::vector<MWPCClusterContainer> m_MWPCClCont;
  std::vector<DCHitContainer>       m_TempBcInHC;
  std::vector<DCHitContainer>       m_BcInHC;
  std::vector<DCHitContainer>       m_BcOutHC;
  std::vector<DCHitContainer>       m_SdcInHC;
  std::vector<DCHitContainer>       m_SdcOutHC;
  std::vector<TPCHitContainer>      m_TPCHitCont;
  std::vector<TPCHitContainer>      m_TempTPCHitCont;
  std::vector<TPCHitContainer>      m_TPCClCont;

  DCHitContainer        m_TOFHC;
  DCHitContainer        m_VtxPoint;
  DCLocalTrackContainer m_BcInTC;
  DCLocalTrackContainer m_BcOutTC;
  DCLocalTrackContainer m_SdcInTC;
  DCLocalTrackContainer m_SdcOutTC;

  TPCLocalTrackContainer m_TPCTC;
  TPCLocalTrack_HelixContainer m_TPCTC_Helix;

  K18TrackU2DContainer  m_K18U2DTC;
  K18TrackD2UContainer  m_K18D2UTC;
  KuramaTrackContainer  m_KuramaTC;
  DCLocalTrackContainer m_BcOutSdcInTC;
  DCLocalTrackContainer m_SdcInSdcOutTC;
  // Exclusive Tracks
  std::vector<DCLocalTrackContainer> m_SdcInExTC;
  std::vector<DCLocalTrackContainer> m_SdcOutExTC;

public:
  int  MuchCombinationSdcIn() const { return m_much_combi[k_SdcIn]; }
  bool DecodeRawHits(RawData* rawData);
  // bool DecodeFiberHits(FiberCluster* FiberCl, int layer);
  bool DecodeFiberHits(RawData* rawData);
  bool DecodeBcInHits(RawData* rawData);
  bool DecodeBcOutHits(RawData* rawData);
  bool DecodeTPCHitsGeant4(const int nhits,
     			    const double *x, const double *y, const double *z, const double *de);
  bool DecodeTPCHits(RawData* rawData);
  bool DecodeSdcInHits(RawData* rawData);
  bool DecodeSdcOutHits(RawData* rawData, double ofs_dt=0.);
  bool DecodeTOFHits(const Hodo2HitContainer& HitCont);
  bool DecodeTOFHits(const HodoClusterContainer& ClCont);
  // bool DecodeSimuHits(SimuData *simuData);
  int  ClusterizeMWPCHit(const DCHitContainer& hits,
			  MWPCClusterContainer& clusters);
  bool  ClusterizeTPC(int layerID, const TPCHitContainer& HitCont,
		       TPCClusterContainer& ClCont);

  inline const DCHitContainer& GetTempBcInHC(int layer) const;
  inline const DCHitContainer& GetBcInHC(int layer) const;
  inline const DCHitContainer& GetBcOutHC(int layer) const;
  inline const DCHitContainer& GetSdcInHC(int layer) const;
  inline const DCHitContainer& GetSdcOutHC(int layer) const;
  inline const DCHitContainer& GetTOFHC() const;
  inline const TPCHitContainer& GetTPCHC(int layer) const;
  inline const TPCHitContainer& GetTPCClCont(int layer) const;

  bool TrackSearchBcIn();
  bool TrackSearchBcIn(const std::vector< std::vector<DCHitContainer> >& hc);
  bool TrackSearchBcOut(Int_t T0Seg=-1);
  bool TrackSearchBcOut(const std::vector< std::vector<DCHitContainer> >& hc, int T0Seg);
  bool TrackSearchSdcIn();
  bool TrackSearchSdcInFiber();
  bool TrackSearchSdcOut();
  bool TrackSearchSdcOut(const Hodo2HitContainer& HitCont);
  bool TrackSearchSdcOut(const HodoClusterContainer& ClCont);
  bool TrackSearchTPC();
  bool TrackSearchTPC_Helix();

  int GetNtracksBcIn()   const { return m_BcInTC.size(); }
  int GetNtracksBcOut()  const { return m_BcOutTC.size(); }
  int GetNtracksSdcIn()  const { return m_SdcInTC.size(); }
  int GetNtracksSdcOut() const { return m_SdcOutTC.size(); }
  int GetNTracksTPC() const { return m_TPCTC.size(); }
  int GetNTracksTPC_Helix() const { return m_TPCTC_Helix.size(); }
  // Exclusive Tracks
  int GetNtracksSdcInEx(int layer) const { return m_SdcInExTC[layer].size(); }
  int GetNtracksSdcOutEx(int layer) const { return m_SdcOutExTC[layer].size(); }

  inline DCLocalTrack* GetTrackBcIn(int i) const;
  inline DCLocalTrack* GetTrackBcOut(int i) const;
  inline DCLocalTrack* GetTrackSdcIn(int i) const;
  inline DCLocalTrack* GetTrackSdcOut(int i) const;
  inline TPCLocalTrack* GetTrackTPC(int i) const;
  inline TPCLocalTrack_Helix* GetTrackTPC_Helix(int i) const;
  // Exclusive Tracks
  inline DCLocalTrack* GetTrackSdcInEx(int layer, int i) const;
  inline DCLocalTrack* GetTrackSdcOutEx(int layer, int i) const;

  bool TrackSearchK18U2D();
  bool TrackSearchK18D2U(const std::vector<double>& XinCont);
  bool TrackSearchKurama(double initial_momentum);
  bool TrackSearchKurama();

  void ChiSqrCutBcOut(double chisqr);
  void ChiSqrCutSdcIn(double chisqr);
  void ChiSqrCutSdcOut(double chisqr);

  void TotCutBCOut(double min_tot);
  void TotCutSDC1(double min_tot);
  void TotCutSDC2(double min_tot);
  void TotCutSDC3(double min_tot);
  void TotCutSDC4(double min_tot);

  void DriftTimeCutBC34(double min_dt, double max_dt);
  void DriftTimeCutSDC1(double min_dt, double max_dt);
  void DriftTimeCutSDC2(double min_dt, double max_dt);
  void DriftTimeCutSDC3(double min_dt, double max_dt);
  void DriftTimeCutSDC4(double min_dt, double max_dt);

  int GetNTracksK18U2D() const { return m_K18U2DTC.size(); }
  int GetNTracksK18D2U() const { return m_K18D2UTC.size(); }
  int GetNTracksKurama() const { return m_KuramaTC.size(); }


  inline K18TrackU2D  * GetK18TrackU2D(int i) const;
  inline K18TrackD2U  * GetK18TrackD2U(int i) const;
  inline KuramaTrack  * GetKuramaTrack(int i)    const;

  int GetNClustersMWPC(int layer) const { return m_MWPCClCont[layer].size(); };

  inline const MWPCClusterContainer & GetClusterMWPC(int layer) const;

  void PrintKurama(const std::string& arg="") const;

  bool ReCalcMWPCHits(std::vector<DCHitContainer>& cont,
		       bool applyRecursively=false);
  bool ReCalcDCHits(std::vector<DCHitContainer>& cont,
		     bool applyRecursively=false);
  bool ReCalcDCHits(bool applyRecursively=false);
  bool ReCalcTPCHits(const int nhits,
                      const std::vector<int>& padid,
                      const std::vector<double>& time,
                      const std::vector<double>& de,
                      Bool_t do_clusterize=true);
  void HoughYCut(double min_y, double max_y);
  bool ReCalcTrack(DCLocalTrackContainer& cont, bool applyRecursively=false);
  bool ReCalcTrack(K18TrackD2UContainer& cont, bool applyRecursively=false);
  bool ReCalcTrack(KuramaTrackContainer& cont, bool applyRecursively=false);

  bool ReCalcTrackBcIn(bool applyRecursively=false);
  bool ReCalcTrackBcOut(bool applyRecursively=false);
  bool ReCalcTrackSdcIn(bool applyRecursively=false);
  bool ReCalcTrackSdcOut(bool applyRecursively=false);

  bool ReCalcK18TrackD2U(bool applyRecursively=false);
  // bool ReCalcK18TrackU2D(bool applyRecursively=false);
  bool ReCalcKuramaTrack(bool applyRecursively=false);

  bool ReCalcAll();

  bool TrackSearchBcOutSdcIn();
  bool TrackSearchSdcInSdcOut();
  int GetNtracksBcOutSdcIn() const { return m_BcOutSdcInTC.size(); }
  int GetNtracksSdcInSdcOut() const { return m_SdcInSdcOutTC.size(); }
  inline DCLocalTrack * GetTrackBcOutSdcIn(int i) const;
  inline DCLocalTrack * GetTrackSdcInSdcOut(int i) const;

  bool MakeBH2DCHit(int t0seg);

protected:
  void ClearDCHits();
  void ClearBcInHits();
  void ClearBcOutHits();
  void ClearSdcInHits();
  void ClearSdcOutHits();

  void ClearTOFHits();
  void ClearVtxHits();

  void ClearTPCHits();
  void ClearTPCClusters();

  void ClearTracksBcIn();
  void ClearTracksBcOut();
  void ClearTracksSdcIn();
  void ClearTracksSdcOut();
  void ClearTracksTPC();
  void ClearTracksBcOutSdcIn();
  void ClearTracksSdcInSdcOut();
  void ClearK18TracksU2D();
  void ClearK18TracksD2U();
  void ClearKuramaTracks();
  void ChiSqrCut(DCLocalTrackContainer& cont, double chisqr);
  void TotCut(DCHitContainer& cont, double min_tot, bool adopt_nan);
  void DriftTimeCut(DCHitContainer& cont, double min_dt, double max_dt, bool select_1st);
  static int MakeUpMWPCClusters(const DCHitContainer& HitCont,
  				 MWPCClusterContainer& ClusterCont,
  				 double maxTimeDif);
public:
  void ResetTracksBcIn()        { ClearTracksBcIn();        }
  void ResetTracksBcOut()       { ClearTracksBcOut();       }
  void ResetTracksSdcIn()       { ClearTracksSdcIn();       }
  void ResetTracksSdcOut()      { ClearTracksSdcOut();      }
  void ResetTracksBcOutSdcIn()  { ClearTracksBcOutSdcIn();  }
  void ResetTracksSdcInSdcOut()  { ClearTracksSdcInSdcOut();  }
  void ApplyBh1SegmentCut(const std::vector<double>& validBh1Cluster);
  void ApplyBh2SegmentCut(const double Time0_Cluster);

};

//_____________________________________________________________________________
inline TString&
DCAnalyzer::ClassName()
{
  static TString s_name("DCAnalyzer");
  return s_name;
}

//_____________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetTempBcInHC(int layer) const
{
  if(layer>NumOfLayersBcIn) layer=0;
  return m_TempBcInHC[layer];
}

//_____________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetBcInHC(int layer) const
{
  if(layer>NumOfLayersBcIn) layer=0;
  return m_BcInHC[layer];
}

//_____________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetBcOutHC(int layer) const
{
  if(layer>NumOfLayersBcOut+1) layer=0;
  return m_BcOutHC[layer];
}

//_____________________________________________________________________________
inline const TPCHitContainer&
DCAnalyzer::GetTPCHC(int layer) const
{
  if(layer>NumOfLayersTPC) layer=NumOfLayersTPC;
  return m_TPCHitCont[layer];
}

//_____________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetSdcInHC(int layer) const
{
  if(layer>NumOfLayersSdcIn) layer=0;
  return m_SdcInHC[layer];
}

//_____________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetSdcOutHC(int layer) const
{
  if(layer>NumOfLayersSdcOut) layer=0;
  return m_SdcOutHC[layer];
}

//_____________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetTOFHC() const
{
  return m_TOFHC;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackBcIn(int i) const
{
  if(i<m_BcInTC.size())
    return m_BcInTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackBcOut(int i) const
{
  if(i<m_BcOutTC.size())
    return m_BcOutTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcIn(int i) const
{
  if(i<m_SdcInTC.size())
    return m_SdcInTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcOut(int i) const
{
  if(i<m_SdcOutTC.size())
    return m_SdcOutTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline TPCLocalTrack*
DCAnalyzer::GetTrackTPC(int i) const
{
  if(i<m_TPCTC.size())
    return m_TPCTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline TPCLocalTrack_Helix*
DCAnalyzer::GetTrackTPC_Helix(int i) const
{
  if(i<m_TPCTC_Helix.size())
    return m_TPCTC_Helix[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcInEx(int layer, int i) const
{
  if(i<m_SdcInExTC[layer].size())
    return m_SdcInExTC[layer][i];
  else
    return 0;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcOutEx(int layer, int i) const
{
  if(i<m_SdcOutExTC[layer].size())
    return m_SdcOutExTC[layer][i];
  else
    return 0;
}

//_____________________________________________________________________________
inline K18TrackU2D*
DCAnalyzer::GetK18TrackU2D(int i) const
{
  if(i<m_K18U2DTC.size())
    return m_K18U2DTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline K18TrackD2U*
DCAnalyzer::GetK18TrackD2U(int i) const
{
  if(i<m_K18D2UTC.size())
    return m_K18D2UTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline KuramaTrack*
DCAnalyzer::GetKuramaTrack(int i) const
{
  if(i<m_KuramaTC.size())
    return m_KuramaTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackBcOutSdcIn(int i) const
{
  if(i<m_BcOutSdcInTC.size())
    return m_BcOutSdcInTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcInSdcOut(int i) const
{
  if(i<m_SdcInSdcOutTC.size())
    return m_SdcInSdcOutTC[i];
  else
    return 0;
}

//_____________________________________________________________________________
inline const MWPCClusterContainer&
DCAnalyzer::GetClusterMWPC(int layer) const
{
  if(layer>NumOfLayersBcIn) layer=0;
  return m_MWPCClCont[layer];
}

//_____________________________________________________________________________
inline const TPCHitContainer&
DCAnalyzer::GetTPCClCont(int layer) const
{
  if(layer>NumOfLayersTPC) layer=NumOfLayersTPC;
  return m_TPCClCont[layer];
}

#endif
