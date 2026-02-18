// -*- C++ -*-
#include "DCTrackSearch.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TH2D.h>
#include <TH3D.h>
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "DCPairHitCluster.hh"
#include "DCParameters.hh"
#include "DebugTimer.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "TrackMaker.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>

namespace
{
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zBH2  = gGeom.LocalZ("BH2");
const Double_t MaxChisquare       = 2000.; // Set to be More than 30
const Double_t MaxNumOfCluster = 10.;    // Set to be Less than 30
const Double_t MaxCombi = 1.0e6;    // Set to be Less than 10^6
// BcOut for XUV Tracking routine
const Double_t MaxChisquareVXU = 50.;//
const Double_t ChisquareCutVXU = 50.;//

// BH2 Geometry Parameters
const Double_t Bh2Width  = 14.0;
const Double_t Bh2Height = 100.0;
const Double_t Bh2YMin   = -50.0;
const Double_t Bh2YMax   =  50.0;

const Double_t localPosBh2X_dX = 0.;
inline std::vector<Double_t> initBh2XPos() {
  std::vector<Double_t> pos;
  pos.reserve(NumOfSegBH2);
  Double_t offset = (NumOfSegBH2 - 1) / 2.0;
  for(Int_t i=0; i<NumOfSegBH2; ++i) {
    pos.push_back((i - offset) * Bh2Width + localPosBh2X_dX);
  }
  return pos;
}

inline std::vector<Double_t> initBh2XAcc() {
  return std::vector<Double_t>(NumOfSegBH2, Bh2Width/2.0 + 1.0);
}

const std::vector<Double_t> Bh2XPos = initBh2XPos();
const std::vector<Double_t> Bh2XAcc = initBh2XAcc();


//_____________________________________________________________________________
template <typename T> void
CalcTracks(std::vector<T*>& trackCont)
{
  for(auto& track: trackCont) track->Calculate();
}

//_____________________________________________________________________________
template <typename T> void
ClearFlags(std::vector<T*>& trackCont)
{
  for(const auto& track: trackCont){
    if(!track) continue;
    Int_t nh = track->GetNHit();
    for(Int_t j=0; j<nh; ++j) track->GetHit(j)->QuitTrack();
  }
}

//_____________________________________________________________________________
inline void
DeleteDuplicatedTracks(DCLocalTC& trackCont, Double_t ChisqrCut=0.)
{
  // evaluate container size in every iteration
  for(Int_t i=0; i<(Int_t)trackCont.size(); ++i){
    const auto& tp = trackCont[i];
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();
    for(Int_t i2=(Int_t)trackCont.size()-1; i2>i; --i2){
      const DCLocalTrack* tp2 = trackCont[i2];
      Int_t nh2 = tp2->GetNHit(), flag=0;
      Double_t chisqr = tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag>0 && chisqr>ChisqrCut){
        delete tp2;
        trackCont.erase(trackCont.begin()+i2);
      }
    }
  }
}

//_____________________________________________________________________________
inline void
DeleteDuplicatedTracks(DCLocalTC& trackCont,
                       Int_t first, Int_t second, Double_t ChisqrCut=0.)
{
  std::vector<Int_t> delete_index;
  // evaluate container size in every iteration
  for(Int_t i=first; i<=second; ++i){

    auto itr = std::find(delete_index.begin(), delete_index.end(), i);
    if(itr != delete_index.end())
      continue;

    const DCLocalTrack* const tp = trackCont[i];
    if(!tp) continue;

    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

    for(Int_t i2=second; i2>i; --i2){
      auto itr = std::find(delete_index.begin(), delete_index.end(), i2);
      if(itr != delete_index.end())
        continue;

      const DCLocalTrack* tp2 = trackCont[i2];
      Int_t nh2 = tp2->GetNHit(), flag=0;
      Double_t chisqr = tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag>0 && chisqr>ChisqrCut){
        // delete tp2; // Be careful, if already deleted in previous loop?
        // tp2 = nullptr;
        delete_index.push_back(i2);
      }
    }
  }

  // sort from bigger order
  std::sort(delete_index.begin(), delete_index.end(), std::greater<Int_t>());
  for(Int_t i=0; i<(Int_t)delete_index.size(); i++) {
    delete trackCont[delete_index[i]];
    trackCont.erase(trackCont.begin()+delete_index[i]);
  }

  // reset hit record of DCHit
  for(Int_t i=0; i<(Int_t)trackCont.size(); ++i){
    const DCLocalTrack* const tp = trackCont[i];
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->QuitTrack();
  }
}

//_____________________________________________________________________________
[[maybe_unused]] void
DebugPrint(const IndexList& nCombi,
           const std::vector<ClusterList>& CandCont,
           const TString& arg="")
{
  if(SPDLOG_ACTIVE_LEVEL > spdlog::level::debug) return;
  spdlog::debug("{} #Hits of each group", arg.Data());
  Int_t np = nCombi.size();
  Int_t nn = 1;
  std::stringstream ss;
  for(Int_t i=0; i<np; ++i){
    ss << std::setw(4) << nCombi[i];
    nn *= nCombi[i] + 1;
  }
  ss << " -> " << nn-1 << " Combinations";
  spdlog::debug(ss.str());
  for(Int_t i=0; i<np; ++i){
    Int_t n=CandCont[i].size();
    ss.str("");
    ss << "[" << std::setw(3) << i << "]: "
       << std::setw(3) << n << " ";
    for(Int_t j=0; j<n; ++j){
      const auto pair = CandCont[i][j];
      const auto nh = pair->NumberOfHits();
      ss << "{";
      for(Int_t k=0; k<nh; ++k){
        const auto hit = pair->GetHit(k);
        ss << hit->GetLayer()
           << ":" << hit->GetWire() << " ";
      }
      ss << "} ";
    }
    spdlog::debug(ss.str());
  }
}

//_____________________________________________________________________________
template <class Functor>
inline void
FinalizeTrack(const TString& arg,
              DCLocalTC& trackCont,
              Functor comp,
              std::vector<ClusterList>& candCont,
              Bool_t delete_flag=true)
{
  ClearFlags(trackCont);

#if 0
  DebugPrint(trackCont, arg+" Before Sorting ");
#endif

  std::stable_sort(trackCont.begin(), trackCont.end(), DCLTrackComp_Nhit());

#if 0
  DebugPrint(trackCont, arg+" After Sorting (Nhit) ");
#endif

  typedef std::pair<Int_t, Int_t> index_pair;
  std::vector<index_pair> index_pair_vec;

  std::vector<Int_t> nhit_vec;

  for(Int_t i=0; i<trackCont.size(); i++) {
    Int_t nhit = trackCont[i]->GetNHit();
    nhit_vec.push_back(nhit);
  }

  if(!nhit_vec.empty()) {
    Int_t max_nhit = nhit_vec.front();
    Int_t min_nhit = nhit_vec.back();
    for(Int_t nhit=max_nhit; nhit>=min_nhit; nhit--) {
      auto itr1 = std::find(nhit_vec.begin(), nhit_vec.end(), nhit);
      if(itr1 == nhit_vec.end())
        continue;

      size_t index1 = std::distance(nhit_vec.begin(), itr1);

      auto itr2 = std::find(nhit_vec.rbegin(), nhit_vec.rend(), nhit);
      size_t index2 = nhit_vec.size() - std::distance(nhit_vec.rbegin(), itr2) - 1;

      index_pair_vec.push_back(index_pair(index1, index2));
    }
  }

  for(Int_t i=0; i<index_pair_vec.size(); i++) {
    std::stable_sort(trackCont.begin() + index_pair_vec[i].first,
                     trackCont.begin() + index_pair_vec[i].second + 1, DCLTrackComp_Chisqr());
  }

#if 0
  DebugPrint(trackCont, arg+" After Sorting (chisqr)");
#endif

  if(delete_flag) {
    for(Int_t i = index_pair_vec.size()-1; i>=0; --i) {
      DeleteDuplicatedTracks(trackCont, index_pair_vec[i].first, index_pair_vec[i].second, 0.);
    }
  }

#if 0
  DebugPrint(trackCont, arg+" After Deleting in each hit number");
#endif

  std::stable_sort(trackCont.begin(), trackCont.end(), comp);

#if 0
  DebugPrint(trackCont, arg+" After Sorting with comp func ");
#endif

  if(delete_flag) DeleteDuplicatedTracks(trackCont);

#if 0
  DebugPrint(trackCont, arg+" After Deleting ");
#endif

  CalcTracks(trackCont);
  del::ClearContainerAll(candCont);
}

//_____________________________________________________________________________
// MakeCluster ________________________________________________________________
//_____________________________________________________________________________
Bool_t
MakePairPlaneHitCluster(const DCHC & HC1,
                        const DCHC & HC2,
                        Double_t CellSize,
                        ClusterList& Cont,
                        Bool_t honeycomb=false)
{
  Int_t nh1=HC1.size(), nh2=HC2.size();
  std::vector<Int_t> UsedFlag(nh2,0);
  for(Int_t i1=0; i1<nh1; ++i1){
    DCHit *hit1=HC1[i1];
    Double_t wp1=hit1->GetWirePosition();
    Bool_t flag=false;
    for(Int_t i2=0; i2<nh2; ++i2){
      DCHit *hit2=HC2[i2];
      Double_t wp2=hit2->GetWirePosition();
      if(std::abs(wp1-wp2)<=CellSize){
        Int_t multi1 = hit1->GetEntries();
        Int_t multi2 = hit2->GetEntries();
        for(Int_t m1=0; m1<multi1; ++m1) {
          if(!hit1->IsGood(m1))
            continue;
          for(Int_t m2=0; m2<multi2; ++m2) {
            if(!hit2->IsGood(m2))
              continue;
            Double_t x1,x2;
            if(wp1<wp2){
              x1=wp1+hit1->DriftLength(m1);
              x2=wp2-hit2->DriftLength(m2);
            }
            else {
              x1=wp1-hit1->DriftLength(m1);
              x2=wp2+hit2->DriftLength(m2);
            }
            DCPairHitCluster *cluster =
              new DCPairHitCluster(new DCLTrackHit(hit1,x1,m1),
                                   new DCLTrackHit(hit2,x2,m2));
            cluster->SetHoneycomb(honeycomb);
            Cont.push_back(cluster);
            flag=true; ++UsedFlag[i2];
          }
        }
      }
    }
#if 1
    if(!flag){
      Int_t multi1 = hit1->GetEntries();
      for(Int_t m1=0; m1<multi1; m1++) {
        if(!(hit1->IsGood(m1))) continue;
        Double_t dl=hit1->DriftLength(m1);
        DCPairHitCluster *cluster1 = new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl,m1));
        DCPairHitCluster *cluster2 = new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl,m1));
        cluster1->SetHoneycomb(honeycomb);
        cluster2->SetHoneycomb(honeycomb);
        Cont.push_back(cluster1);
        Cont.push_back(cluster2);
      }
    }
#endif
  }
#if 1
  for(Int_t i2=0; i2<nh2; ++i2){
    if(UsedFlag[i2]==0) {
      DCHit *hit2=HC2[i2];
      Int_t multi2 = hit2->GetEntries();
      for(Int_t m2=0; m2<multi2; m2++) {
        if(!(hit2->IsGood(m2))) continue;
        Double_t wp=hit2->GetWirePosition();
        Double_t dl=hit2->DriftLength(m2);
        DCPairHitCluster *cluster1 = new DCPairHitCluster(new DCLTrackHit(hit2,wp+dl,m2));
        DCPairHitCluster *cluster2 = new DCPairHitCluster(new DCLTrackHit(hit2,wp-dl,m2));
        cluster1->SetHoneycomb(honeycomb);
        cluster2->SetHoneycomb(honeycomb);
        Cont.push_back(cluster1);
        Cont.push_back(cluster2);
      }
    }
  }
#endif
  return true;
}

//_____________________________________________________________________________
Bool_t
MakeUnPairPlaneHitCluster(const DCHC& HC,
                          ClusterList& Cont,
                          Bool_t honeycomb=false)
{
  const std::size_t nh = HC.size();
  for(std::size_t i=0; i<nh; ++i){
    DCHit *hit = HC[i];
    if(!hit) continue;
    std::size_t mh = hit->GetEntries();
    for(std::size_t m=0; m<mh; ++m) {
      if(!hit->IsGood(m)) continue;
      Double_t wp = hit->GetWirePosition();
      Double_t dl = hit->DriftLength(m);
      DCPairHitCluster *cluster1 =
        new DCPairHitCluster(new DCLTrackHit(hit,wp+dl,m));
      DCPairHitCluster *cluster2 =
        new DCPairHitCluster(new DCLTrackHit(hit,wp-dl,m));
      cluster1->SetHoneycomb(honeycomb);
      cluster2->SetHoneycomb(honeycomb);
      Cont.push_back(cluster1);
      Cont.push_back(cluster2);
    }
  }

  return true;
}


//_____________________________________________________________________________
Bool_t
MakePairPlaneHitClusterVUX(const DCHC& HC1,
                           const DCHC& HC2,
                           Double_t CellSize,
                           ClusterList& Cont,
                           Bool_t honeycomb=false)
{
  Int_t nh1=HC1.size(), nh2=HC2.size();
  std::vector<Int_t> UsedFlag(nh2,0);

  for(Int_t i1=0; i1<nh1; ++i1){
    DCHit *hit1=HC1[i1];

    Double_t wp1=hit1->GetWirePosition();
    Bool_t flag=false;
    for(Int_t i2=0; i2<nh2; ++i2){
      DCHit *hit2=HC2[i2];
      Double_t wp2=hit2->GetWirePosition();
      if(std::abs(wp1-wp2)<CellSize){

        Int_t multi1 = hit1->GetEntries();
        Int_t multi2 = hit2->GetEntries();
        for(Int_t m1=0; m1<multi1; m1++) {
          if(!(hit1->IsGood(m1))) continue;
          for(Int_t m2=0; m2<multi2; m2++) {
            if(!(hit2->IsGood(m2))) continue;
            Double_t dl1=hit1->DriftLength(m1);
            Double_t dl2=hit2->DriftLength(m2);

            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl1,m1),
                                                new DCLTrackHit(hit2,wp2+dl2,m2)));
            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl1,m1),
                                                new DCLTrackHit(hit2,wp2-dl2,m2)));
            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl1,m1),
                                                new DCLTrackHit(hit2,wp2+dl2,m2)));
            Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl1,m1),
                                                new DCLTrackHit(hit2,wp2-dl2,m2)));

            flag=true; ++UsedFlag[i2];
          }
        }
      }
    }
    if(!flag){
      Int_t multi1 = hit1->GetEntries();
      for(Int_t m1=0; m1<multi1; m1++) {
        if(!(hit1->IsGood(m1))) continue;
        Double_t dl=hit1->DriftLength(m1);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl,m1)));
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl,m1)));
      }
    }
  }
  for(Int_t i2=0; i2<nh2; ++i2){
    if(UsedFlag[i2]==0) {
      DCHit *hit2=HC2[i2];
      Int_t multi2 = hit2->GetEntries();
      for(Int_t m2=0; m2<multi2; m2++) {
        if(!(hit2->IsGood(m2))) continue;

        Double_t wp=hit2->GetWirePosition();
        Double_t dl=hit2->DriftLength(m2);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit2,wp+dl,m2)));
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit2,wp-dl,m2)));
      }
    }
  }

  return true;
}

} // namespace

//_____________________________________________________________________________
namespace track
{

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex(Int_t ndim, const Int_t *index1, Bool_t& status)
{
  if(ndim==1){
    std::vector<IndexList> index2;
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector<IndexList> index2 = MakeIndex(ndim-1, index1+1, status);

  std::vector<IndexList> index;
  Int_t n2=index2.size();
  for(Int_t j=0; j<n2; ++j){
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem;
      Int_t n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      for(Int_t k=0; k<n3; ++k)
        elem.push_back(index2[j][k]);
      index.push_back(elem);
      Int_t size1=index.size();
      if(size1>MaxCombi){
        status = false;
#if 0
        hddaq::cout
          // << hddaq::unpacker::esc::k_yellow
          << FUNC_NAME << " too much combinations... " << n2
          // << hddaq::unpacker::esc::k_default_color
          << std::endl;
#endif
        return std::vector<IndexList>(0);
      }
    }
  }

  return index;
}

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex(Int_t ndim, const IndexList& index1, Bool_t& status)
{
  return MakeIndex(ndim, &(index1[0]), status);
}

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex_VXU(Int_t ndim,Int_t maximumHit, const Int_t *index1)
{
  if(ndim==1){
    std::vector<IndexList> index2;
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector<IndexList>
    index2=MakeIndex_VXU(ndim-1, maximumHit, index1+1);

  std::vector<IndexList> index;
  Int_t n2=index2.size();
  for(Int_t j=0; j<n2; ++j){
    for(Int_t i=-1; i<index1[0]; ++i){
      IndexList elem;
      Int_t validHitNum=0;
      Int_t n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if(i != -1)
        validHitNum++;
      for(Int_t k=0; k<n3; ++k){
        elem.push_back(index2[j][k]);
        if(index2[j][k] != -1)
          validHitNum++;
      }
      if(validHitNum <= maximumHit)
        index.push_back(elem);
    }
  }

  return index;
}

//_____________________________________________________________________________
std::vector<IndexList>
MakeIndex_VXU(Int_t ndim,Int_t maximumHit, const IndexList& index1)
{
  return MakeIndex_VXU(ndim, maximumHit, &(index1[0]));
}

//_____________________________________________________________________________
DCLocalTrack*
MakeTrack(const std::vector<ClusterList>& CandCont,
          const IndexList& combination)
{
  DCLocalTrack *tp = new DCLocalTrack;
  for(std::size_t i=0, n=CandCont.size(); i<n; ++i){
    Int_t m = combination[i];
    if(m<0) continue;
    DCPairHitCluster *cluster = CandCont[i][m];
    if(!cluster) continue;
    Int_t mm = cluster->NumberOfHits();
    for(Int_t j=0; j<mm; ++j){
      DCLTrackHit *hitp = cluster->GetHit(j);
      if(!hitp) continue;
      tp->AddHit(hitp);
    }
#if 0
    hddaq::cout << FUNC_NAME << ":" << std::setw(3)
                << i << std::setw(3) << m  << " "
                << CandCont[i][m] << " " << mm << std::endl;
#endif
  }
  return tp;
}

//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<DCHC>& HC,
                 const DCPairPlaneInfo * PpInfo,
                 Int_t npp, DCLocalTC& TrackCont,
                 Int_t MinNumOfHits, Int_t T0Seg)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Bool_t fiber     = PpInfo[i].fiber;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    if(ppFlag && !fiber){
      MakePairPlaneHitCluster(HC[layer1], HC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else if(!ppFlag && fiber){
      PpInfo[i].Print(FUNC_NAME+" invalid parameter", hddaq::cerr);
    }else{
      MakeUnPairPlaneHitCluster(HC[layer1], CandCont[i], honeycomb);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);
  if (!status) DebugPrint(nCombi, CandCont, FUNC_NAME);

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits
       && track->DoFit()
       && track->GetChiSquare()<MaxChisquare){
      if(0 <= T0Seg && T0Seg < (Int_t)Bh2XPos.size()) {
        Double_t xbh2 = track->GetX(zBH2);
        Double_t ybh2 = track->GetY(zBH2);
        Double_t difPosBh2 = Bh2XPos[T0Seg] - xbh2;
        if(true
           && std::abs(difPosBh2) < Bh2XAcc[T0Seg]
           && (Bh2YMin < ybh2 && ybh2 < Bh2YMax)
          ){
          TrackCont.push_back(track);
        }else{
          delete track;
        }
      }else{
        TrackCont.push_back(track);
      }
    }
    else{
      delete track;
    }
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status ? TrackCont.size() : -1;
}

// BC3&4 Tracking ___________________________________________
Int_t
LocalTrackSearchVUX(const std::vector<DCHC>& HC,
                    const DCPairPlaneInfo* PpInfo,
                    Int_t npp, DCLocalTC& TrackCont,
                    Int_t MinNumOfHits /*=6*/, Int_t T0Seg /*=-1*/)
{
  DCLocalTC TrackContV;
  DCLocalTC TrackContX;
  DCLocalTC TrackContU;
  std::vector<ClusterList>   CandCont(npp);
  std::vector<ClusterList>   CandContV(npp);
  std::vector<ClusterList>   CandContX(npp);
  std::vector<ClusterList>   CandContU(npp);

  //  Int_t NumOfLayersDC_12 = 12 ;
  Int_t iV=0, iX=0, iU=0;
  Int_t nV=0, nX=0, nU=0;
  for(Int_t i=0; i<npp; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    Int_t  nh1       = HC[layer1].size();
    Int_t  nh2       = HC[layer2].size();
    Double_t TiltAngle = 0.;
    if(nh1>0) TiltAngle = HC[layer1][0]->GetTiltAngle();
    if(ppFlag && nh1==0 && nh2>0)
      TiltAngle = HC[layer2][0]->GetTiltAngle();
    if(ppFlag && nh1>0 && nh2>0){
      if(TiltAngle<0){
        MakePairPlaneHitClusterVUX(HC[layer1], HC[layer2],
                                   PpInfo[i].CellSize, CandContV[iV],
                                   honeycomb);
        ++iV; nV = nV+2;
      }
      if(TiltAngle==0){
        MakePairPlaneHitClusterVUX(HC[layer1], HC[layer2],
                                   PpInfo[i].CellSize, CandContX[iX],
                                   honeycomb);
        ++iX; nX = nX+2;
      }
      if(TiltAngle>0){
        MakePairPlaneHitClusterVUX(HC[layer1], HC[layer2],
                                   PpInfo[i].CellSize, CandContU[iU],
                                   honeycomb);
        ++iU; nU = nU+2;
      }
    }
    if(!ppFlag){
      if(TiltAngle<0){
        MakeUnPairPlaneHitCluster(HC[layer1], CandContV[iV], honeycomb);
        ++nV; ++iV;
      }
      if(TiltAngle==0){
        MakeUnPairPlaneHitCluster(HC[layer1], CandContX[iX], honeycomb);
        ++nX; ++iX;
      }
      if(TiltAngle>0){
        MakeUnPairPlaneHitCluster(HC[layer1], CandContU[iU], honeycomb);
        ++nU; ++iU;
      }
    }
    if(ppFlag && nh1==0 && nh2>0){
      if(TiltAngle<0){
        MakeUnPairPlaneHitCluster(HC[layer2], CandContV[iV], honeycomb);
        ++nV; ++iV;
      }
      if(TiltAngle==0){
        MakeUnPairPlaneHitCluster(HC[layer2], CandContX[iX], honeycomb);
        ++nX; ++iX;
      }
      if(TiltAngle>0){
        MakeUnPairPlaneHitCluster(HC[layer2], CandContU[iU], honeycomb);
        ++nU; ++iU;
      }
    }
  }

  IndexList nCombi(npp);
  IndexList nCombiV(npp);
  IndexList nCombiX(npp);
  IndexList nCombiU(npp);

  for(Int_t i=0; i<npp; ++i){
    nCombiV[i]=(CandContV[i]).size();
    nCombiX[i]=(CandContX[i]).size();
    nCombiU[i]=(CandContU[i]).size();
  }

#if 0
  DebugPrint(nCombiV, CandContV, FUNC_NAME+" V");
  DebugPrint(nCombiX, CandContX, FUNC_NAME+" X");
  DebugPrint(nCombiU, CandContU, FUNC_NAME+" U");
#endif

  Bool_t status[3] = {true, true, true};
  std::vector<IndexList> CombiIndexV = MakeIndex(npp, nCombiV, status[0]);
  Int_t nnCombiV=CombiIndexV.size();
  std::vector<IndexList> CombiIndexX = MakeIndex(npp, nCombiX, status[1]);
  Int_t nnCombiX=CombiIndexX.size();
  std::vector<IndexList> CombiIndexU = MakeIndex(npp, nCombiU, status[2]);
  Int_t nnCombiU=CombiIndexU.size();

  for(Int_t i=0; i<nnCombiV; ++i){
    DCLocalTrack *track = MakeTrack(CandContV, CombiIndexV[i]);
    if(!track) continue;
    if(track->GetNHit()>=3 && track->DoFitVXU() &&
       track->GetChiSquare()<MaxChisquareVXU){
      TrackContV.push_back(track);
    }
    else{
      delete track;
    }
  }

  for(Int_t i=0; i<nnCombiX; ++i){
    DCLocalTrack *track = MakeTrack(CandContX, CombiIndexX[i]);
    if(!track) continue;
    if(track->GetNHit()>=3 && track->DoFitVXU() &&
       track->GetChiSquare()<MaxChisquareVXU){
      TrackContX.push_back(track);
    }
    else{
      delete track;
    }
  }

  for(Int_t i=0; i<nnCombiU; ++i){
    DCLocalTrack *track = MakeTrack(CandContU, CombiIndexU[i]);
    if(!track) continue;
    if(track->GetNHit()>=3 && track->DoFitVXU() &&
       track->GetChiSquare()<MaxChisquareVXU){
      TrackContU.push_back(track);
    }
    else{
      delete track;
    }
  }

  // Clear Flags
  if(nV>3) ClearFlags(TrackContV);
  if(nX>3) ClearFlags(TrackContX);
  if(nU>3) ClearFlags(TrackContU);

  std::stable_sort(TrackContV.begin(), TrackContV.end(), DCLTrackComp1());
  std::stable_sort(TrackContX.begin(), TrackContX.end(), DCLTrackComp1());
  std::stable_sort(TrackContU.begin(), TrackContU.end(), DCLTrackComp1());

#if 0
  DebugPrint(TrackContV, FUNC_NAME+" V After Sorting.");
  DebugPrint(TrackContX, FUNC_NAME+" X After Sorting.");
  DebugPrint(TrackContU, FUNC_NAME+" U After Sorting.");
#endif

  // Delete Duplicated Tracks (cut chisqr>100 & flag)
  Double_t chiV = ChisquareCutVXU;
  Double_t chiX = ChisquareCutVXU;
  Double_t chiU = ChisquareCutVXU;

  DeleteDuplicatedTracks(TrackContV, chiV);
  DeleteDuplicatedTracks(TrackContX, chiX);
  DeleteDuplicatedTracks(TrackContU, chiU);
  CalcTracks(TrackContV);
  CalcTracks(TrackContX);
  CalcTracks(TrackContU);

#if 0
  DebugPrint(TrackContV, FUNC_NAME+" V After Delete.");
  DebugPrint(TrackContX, FUNC_NAME+" X After Delete.");
  DebugPrint(TrackContU, FUNC_NAME+" U After Delete.");
#endif

  Int_t nnV=1, nnX=1, nnU=1;
  Int_t nkV=0, nkX=0, nkU=0;
  Int_t nnVT=1, nnXT=1, nnUT=1;
  Int_t checkV=0, checkX=0, checkU=0;

  Int_t cV=TrackContV.size();
  if(chiV>1.5 || cV<5) checkV++;
  Int_t cX=TrackContX.size();
  if(chiX>1.5 || cX<5) checkX++;
  Int_t cU=TrackContU.size();
  if(chiU>1.5 || cU<5) checkU++;

  std::vector<IndexList> CombiIndexSV;
  std::vector<IndexList> CombiIndexSX;
  std::vector<IndexList> CombiIndexSU;

  {
    if((nV>=3) && (cV)){
      nnV = TrackContV.size();
      nnVT = TrackContV.size();
      ++nkV;
    }
    if((nV>0) && checkV){
      CombiIndexSV = MakeIndex_VXU(npp, 2, nCombiV);
      nnV = nnV +  CombiIndexSV.size();
    }
  }

  {
    if((nX>=3) && (cX)){
      nnX = TrackContX.size();
      nnXT = TrackContX.size();
      ++nkX;
    }
    if((nX>0) && checkX){
      CombiIndexSX = MakeIndex_VXU(npp, 2, nCombiX);
      nnX = nnX + CombiIndexSX.size();
    }
  }

  {
    if((nU>=3) && (cU)){
      nnU = TrackContU.size();
      nnUT = TrackContU.size();
      ++nkU;
    }
    if((nU>0) && checkU){
      CombiIndexSU = MakeIndex_VXU(npp, 2, nCombiU);
      nnU = nnU + CombiIndexSU.size();
    }
  }

  // Double_t DifVXU=0.0;
  // Double_t Av=0.0;
  // Double_t Ax=0.0;
  // Double_t Au=0.0;

  Double_t chiv, chix, chiu;

#if 0
  for(Int_t i=0; i<nnV; ++i){
    for(Int_t j=0; j<nnX; ++j){
      for(Int_t k=0; k<nnU; ++k){

        chiv=-1.0,chix=-1.0,chiu=-1.0;

        DCLocalTrack *track = new DCLocalTrack();

        /* V Plane  */
        if(nkV){
          DCLocalTrack *trackV=TrackContV[i];
          //Av=trackV->GetVXU_A();
          chiv=trackV->GetChiSquare();
          for(Int_t l=0; l<(trackV->GetNHit()); ++l){
            DCLTrackHit *hitpV=trackV->GetHit(l);
            if(hitpV){
              track->AddHit(hitpV) ;
            }
          }
          //delete trackV;
        }
        if((!nkV) && (nV>0)){
          DCLocalTrack *trackV = MakeTrack(CandContV, CombiIndexV[i]);
          for(Int_t l=0; l<(trackV->GetNHit()); ++l){
            DCLTrackHit *hitpV=trackV->GetHit(l);
            if(hitpV){
              track->AddHit(hitpV) ;
            }
          }
          delete trackV;
        }

        /* X Plane  */
        if(nkX){
          DCLocalTrack *trackX=TrackContX[j];
          //Ax=trackX->GetVXU_A();
          chix=trackX->GetChiSquare();
          for(Int_t l=0; l<(trackX->GetNHit()); ++l){
            DCLTrackHit *hitpX=trackX->GetHit(l);
            if(hitpX){
              track->AddHit(hitpX) ;
            }
          }
          //delete trackX;
        }
        if((!nkX) && (nX>0)){
          DCLocalTrack *trackX = MakeTrack(CandContX, CombiIndexX[j]);
          for(Int_t l=0; l<(trackX->GetNHit()); ++l){
            DCLTrackHit *hitpX=trackX->GetHit(l);
            if(hitpX){
              track->AddHit(hitpX) ;
            }
          }
          delete trackX;
        }

        /* U Plane  */
        if(nkU){
          DCLocalTrack *trackU=TrackContU[k];
          //Au=trackU->GetVXU_A();
          chiu=trackU->GetChiSquare();
          for(Int_t l=0; l<(trackU->GetNHit()); ++l){
            DCLTrackHit *hitpU=trackU->GetHit(l);
            if(hitpU){
              track->AddHit(hitpU) ;
            }
          }
          //delete trackU;
        }
        if((!nkU) && (nU>0)){
          DCLocalTrack *trackU = MakeTrack(CandContU, CombiIndexU[k]);
          for(Int_t l=0; l<(trackU->GetNHit()); ++l){
            DCLTrackHit *hitpU=trackU->GetHit(l);
            if(hitpU){
              track->AddHit(hitpU) ;
            }
          }
          delete trackU;
        }
        //track->SetAv(Av);
        //track->SetAx(Ax);
        //track->SetAu(Au);
        //DifVXU = track->GetDifVXU();

        track->SetChiv(chiv);
        track->SetChix(chix);
        track->SetChiu(chiu);

        if(!track) continue;
        if(track->GetNHit()>=MinNumOfHits && track->DoFit() &&
           track->GetChiSquare()<MaxChisquare){
          TrackCont.push_back(track);
        }
        else{
          delete track;
        }
      }
    }
  }
#endif

  for(Int_t i=-1; i<nnV; ++i){
    for(Int_t j=-1; j<nnX; ++j){
      for(Int_t k=-1; k<nnU; ++k){
        if(((i+j)==-2) || ((j+k)==-2) || ((k+i)==-2)) continue;

        chiv=-1.0,chix=-1.0,chiu=-1.0;

        DCLocalTrack *track = new DCLocalTrack();

        /* V Plane  */
        if(i>-1){
          if(nkV && i<nnVT){
            DCLocalTrack *trackV=TrackContV[i];
            // Av=trackV->GetVXU_A();
            chiv=trackV->GetChiSquare();
            for(Int_t l=0; l<(trackV->GetNHit()); ++l){
              DCLTrackHit *hitpV=trackV->GetHit(l);
              if(hitpV){
                track->AddHit(hitpV) ;
              }
            }
          }
          if((i>=nnVT) && (nV>0)){
            DCLocalTrack *trackV = MakeTrack(CandContV, CombiIndexSV[i-nnVT]);
            for(Int_t l=0; l<(trackV->GetNHit()); ++l){
              DCLTrackHit *hitpV=trackV->GetHit(l);
              if(hitpV){
                track->AddHit(hitpV) ;
              }
            }
            delete trackV ;
          }
          // Int_t NHitV = track->GetNHit();
        }

        /* X Plane  */
        if(j>-1){
          if(nkX && j<nnXT){
            DCLocalTrack *trackX=TrackContX[j];
            // Ax=trackX->GetVXU_A();
            chix=trackX->GetChiSquare();
            for(Int_t l=0; l<(trackX->GetNHit()); ++l){
              DCLTrackHit *hitpX=trackX->GetHit(l);
              if(hitpX){
                track->AddHit(hitpX) ;
              }
            }
          }
          if((j>=nnXT) && (nX>0)){
            DCLocalTrack *trackX = MakeTrack(CandContX, CombiIndexSX[j-nnXT]);
            for(Int_t l=0; l<(trackX->GetNHit()); ++l){
              DCLTrackHit *hitpX=trackX->GetHit(l);
              if(hitpX){
                track->AddHit(hitpX) ;
              }
            }
            delete trackX;
          }
          // Int_t NHitX = track->GetNHit();
        }
        /* U Plane  */
        if(k>-1){
          if(nkU && k<nnUT){
            DCLocalTrack *trackU=TrackContU[k];
            // Au=trackU->GetVXU_A();
            chiu=trackU->GetChiSquare();
            for(Int_t l=0; l<(trackU->GetNHit()); ++l){
              DCLTrackHit *hitpU=trackU->GetHit(l);
              if(hitpU){
                track->AddHit(hitpU) ;
              }
            }
          }
          if((k>=nnUT) && (nU>0)){
            DCLocalTrack *trackU = MakeTrack(CandContU, CombiIndexSU[k-nnUT]);
            for(Int_t l=0; l<(trackU->GetNHit()); ++l){
              DCLTrackHit *hitpU=trackU->GetHit(l);
              if(hitpU){
                track->AddHit(hitpU) ;
              }
            }
            delete trackU;
          }
          // Int_t NHitU = track->GetNHit();
        }

        //track->SetAv(Av);
        //track->SetAx(Ax);
        //track->SetAu(Au);
        //DifVXU = track->GetDifVXU();

        track->SetChiv(chiv);
        track->SetChix(chix);
        track->SetChiu(chiu);

        if(!track) continue;
        if(track->GetNHit()>=MinNumOfHits && track->DoFit() &&
           track->GetChiSquare()<MaxChisquare){//MaXChisquare
          TrackCont.push_back(track);
        }
        else{
          delete track;
        }
      }
    }
  }

  {

    for(Int_t i=0; i<Int_t(TrackContV.size()); ++i){
      DCLocalTrack *tp=TrackContV[i];
      delete tp;
      TrackContV.erase(TrackContV.begin()+i);
    }
    for(Int_t i=0; i<Int_t(TrackContX.size()); ++i){
      DCLocalTrack *tp=TrackContX[i];
      delete tp;
      TrackContX.erase(TrackContX.begin()+i);
    }
    for(Int_t i=0; i<Int_t(TrackContU.size()); ++i){
      DCLocalTrack *tp=TrackContU[i];
      delete tp;
      TrackContU.erase(TrackContU.begin()+i);
    }

  }

  ClearFlags(TrackCont);
  std::stable_sort(TrackCont.begin(), TrackCont.end(), DCLTrackComp1());

#if 1
  // Delete Tracks about  (Nhit1 > Nhit2+1) (Nhit1 > Nhit2  && chi1 < chi2)
  for(Int_t i=0; i<Int_t(TrackCont.size()); ++i){
    DCLocalTrack *tp=TrackCont[i];
    Int_t nh=tp->GetNHit();
    Double_t chi=tp->GetChiSquare();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

    for(Int_t i2=TrackCont.size()-1; i2>i; --i2){
      DCLocalTrack *tp2=TrackCont[i2];
      Int_t nh2=tp2->GetNHit(), flag=0;
      Double_t chi2=tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if((flag>=2) && ((nh==nh2) || ((nh>nh2) && (chi<chi2)))){
        //      if((flag) && ((nh>nh2+1) || ((nh==nh2) || (nh>nh2) && (chi<chi2)))){
        //if((nh>nh2) && (chi<chi2)){
        delete tp2;
        TrackCont.erase(TrackCont.begin()+i2);
      }
    }
  }
#endif

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);

  std::stable_sort(TrackCont.begin(), TrackCont.end(), DCLTrackComp());

  // Clear Flags
  {
    Int_t nbefore=TrackCont.size();
    for(Int_t i=0; i<nbefore; ++i){
      DCLocalTrack *tp=TrackCont[i];
      Int_t nh=tp->GetNHit();
      for(Int_t j=0; j<nh; ++j){
        tp->GetHit(j)->QuitTrack();
      }
    }
  }

  // Delete Duplicated Tracks
  for(Int_t i=0; i<Int_t(TrackCont.size()); ++i){
    DCLocalTrack *tp=TrackCont[i];
    Int_t nh=tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

    for(Int_t i2=TrackCont.size()-1; i2>i; --i2){
      DCLocalTrack *tp2=TrackCont[i2];
      Int_t nh2=tp2->GetNHit(), flag=0;
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag){
        delete tp2;
        TrackCont.erase(TrackCont.begin()+i2);
      }
    }
  }

  // Clear Flags
  {
    Int_t nbefore=TrackCont.size();
    for(Int_t i=0; i<nbefore; ++i){
      DCLocalTrack *tp=TrackCont[i];
      Int_t nh=tp->GetNHit();
      for(Int_t j=0; j<nh; ++j){
        tp->GetHit(j)->QuitTrack();
      }
    }
  }

  {
    Int_t nn=TrackCont.size();
    for(Int_t i=0; i<nn; ++i){
      DCLocalTrack *tp=TrackCont[i];
      Int_t nh=tp->GetNHit();
      for(Int_t j=0; j<nh; ++j){
        Int_t lnum = tp->GetHit(j)->GetLayer();
        Double_t zz = gGeom.GetLocalZ(lnum);
        tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  if(TrackCont.size()>0){
    Int_t nn=TrackCont.size();
    hddaq::cout << FUNC_NAME << ": After Deleting. #Tracks = " << nn << std::endl;
    for(Int_t i=0; i<nn; ++i){
      DCLocalTrack *track=TrackCont[i];
      hddaq::cout << std::setw(3) << i << " #Hits="
                  << std::setw(2) << track->GetNHit()
                  << " ChiSqr=" << track->GetChiSquare()
                  << std::endl;
      hddaq::cout << std::endl;
      for(Int_t j=0; j<(track->GetNHit()); ++j){
        DCLTrackHit *hit = track->GetHit(j);
        hddaq::cout << "layer = " << hit->GetLayer()+1 << " Res = " << hit->GetResidual() << std::endl ;
        hddaq::cout << "hitp = " << hit->GetLocalHitPos() << " calp = " << hit->GetLocalCalPos() << std::endl ;
        hddaq::cout << "X = " << hit->GetXcal() << " Y = " << hit->GetYcal() << std::endl ;
        hddaq::cout << std::endl;
      }
      hddaq::cout << "*********************************************" << std::endl;
    }
    hddaq::cout << std::endl;
  }
#endif

  del::ClearContainerAll(CandCont);
  del::ClearContainerAll(CandContV);
  del::ClearContainerAll(CandContX);
  del::ClearContainerAll(CandContU);

  Bool_t status_all = true;
  status_all = status_all && status[0];
  status_all = status_all && status[1];
  status_all = status_all && status[2];

  return status_all? TrackCont.size() : -1;
}

} // namespace track