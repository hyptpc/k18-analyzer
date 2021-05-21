// -*- C++ -*-

#ifndef DC_TRACK_SEARCH_HH
#define DC_TRACK_SEARCH_HH

#include "DCAnalyzer.hh"

#include <vector>

#include <TString.h>

struct DCPairPlaneInfo;
class DCPairHitCluster;
class DCLocalTrack;
class DCLTrackHit;
class MWPCCluster;
class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrack_Helix;

// typedef std::vector<TPCCluster*>       TPCClusterList;
typedef std::vector<DCPairHitCluster*> ClusterList;
typedef std::vector<int>               IndexList;

namespace track
{
inline const TString& ClassName()
{
  static TString s_name("DCTrackSearch");
  return s_name;
}

std::vector<IndexList> MakeIndex(int ndim, const int *index1, bool& status);
std::vector<IndexList> MakeIndex(int ndim, const IndexList& index1, bool& status);
std::vector<IndexList> MakeIndex_VXU(int ndim,int maximumHit, const int *index1);
std::vector<IndexList> MakeIndex_VXU(int ndim,int maximumHit, const IndexList& index1);
DCLocalTrack*          MakeTrack(const std::vector<ClusterList>& CandCont,
                                 const IndexList& combination);

int LocalTrackSearch(const std::vector<DCHitContainer>& HC,
                     const DCPairPlaneInfo *PpInfo,
                     int npp, std::vector<DCLocalTrack*>& TrackCont,
                     int MinNumOfHits=6, int T0Seg = -1);
int LocalTrackSearch(const std::vector< std::vector<DCHitContainer> >& hcAssemble,
                     const DCPairPlaneInfo *PpInfo,
                     int npp, std::vector<DCLocalTrack*>& TrackCont,
                     int MinNumOfHits=6, int T0Seg = -1);
int LocalTrackSearchSdcInFiber(const std::vector<DCHitContainer>& HC,
                               const DCPairPlaneInfo *PpInfo,
                               int npp, std::vector<DCLocalTrack*>& trackCont,
                               int MinNumOfHits=6);
int LocalTrackSearchVUX(const std::vector<DCHitContainer>& HC,
                        const DCPairPlaneInfo *PpInfo,
                        int npp, std::vector<DCLocalTrack*>& TrackCont,
                        int MinNumOfHits=6);
int LocalTrackSearchSdcOut(const std::vector<DCHitContainer>& SdcOutHC,
                           const DCPairPlaneInfo *PpInfo, int npp,
                           std::vector<DCLocalTrack*>& TrackCont,
                           int MinNumOfHits=6);
int LocalTrackSearchSdcOut(const DCHitContainer& TOFHC,
                           const std::vector<DCHitContainer>& SdcOutHC,
                           const DCPairPlaneInfo *PpInfo,
                           int npp,
                           std::vector<DCLocalTrack*>& TrackCont,
                           int MinNumOfHits=6);
int LocalTrackSearchBcOutSdcIn(const std::vector<DCHitContainer>& BcHC,
                               const DCPairPlaneInfo *BcPpInfo,
                               const std::vector<DCHitContainer>& SdcHC,
                               const DCPairPlaneInfo *SdcPpInfo,
                               int BcNpp, int SdcNpp,
                               std::vector<DCLocalTrack*>& TrackCont,
                               int MinNumOfHits=18);
int LocalTrackSearchSdcInSdcOut(const std::vector<DCHitContainer>& SdcInHC,
                                const DCPairPlaneInfo *SdcInPpInfo,
                                const std::vector<DCHitContainer>& SdcOutHC,
                                const DCPairPlaneInfo *SdcOutPpInfo,
                                int SdcInNpp, int SdcOutNpp,
                                std::vector<DCLocalTrack*>& TrackCont,
                                int MinNumOfHits=12);
int MWPCLocalTrackSearch(const std::vector<DCHitContainer>& HC,
                         std::vector<DCLocalTrack*>& TrackCont);
int MWPCLocalTrackSearch(const std::vector< std::vector<DCHitContainer> >& hcList,
                         std::vector<DCLocalTrack*>& trackCont);
int LocalTrackSearchCFT(const std::vector<DCHitContainer>& HC,
                        const DCPairPlaneInfo *PpInfo,
                        int npp, std::vector<DCLocalTrack*>& trackCont,
                        int MinNumOfHits=3);
int LocalTrackSearchCFTppPhi(const std::vector<DCHitContainer>& HC,
                             const DCPairPlaneInfo *PpInfo,
                             int npp, std::vector<DCLocalTrack*>& trackCont,
                             int MinNumOfHits=3);
int LocalTrackSearchTPC(const std::vector<TPCHitContainer>& TPCHC,
                        std::vector<TPCLocalTrack*>& TrackCont,
                        int MinNumOfHits=8);
int LocalTrackSearchTPC_Helix(const std::vector<TPCHitContainer>& TPCHC,
                              std::vector<TPCLocalTrack_Helix*>& TrackCont,
                              int MinNumOfHits=8);

}

#endif
