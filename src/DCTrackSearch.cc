// -*- C++ -*-

#include "DCTrackSearch.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TH3D.h>

#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "TPCLTrackHit.hh"
#include "DCPairHitCluster.hh"
#include "DCParameters.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "MathTools.hh"
#include "MWPCCluster.hh"
#include "TrackMaker.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrack_Helix.hh"
#include "TPCCluster.hh"

#include "RootHelper.hh"


namespace
{
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zTarget    = gGeom.LocalZ("Target");
const auto& zK18Target = gGeom.LocalZ("K18Target");
const auto& zBH2       = gGeom.LocalZ("BH2");
const Double_t MaxChisquare       = 2000.; // Set to be More than 30
const Double_t MaxChisquareSdcIn  = 5000.; // Set to be More than 30
const Double_t MaxChisquareSdcOut = 50.;   // Set to be More than 30
const Double_t MaxNumOfCluster = 20.;    // Set to be Less than 30
const Double_t MaxCombi = 1.0e6;    // Set to be Less than 10^6
// SdcIn & BcOut for XUV Tracking routine
const Double_t MaxChisquareVXU = 50.;//
const Double_t ChisquareCutVXU = 50.;//

// TPC Tracking
const Int_t    MaxNumOfTrackTPC = 100;
const auto& valueHall = ConfMan::Get<Double_t>("HSFLDHALL");

const Double_t Bh2SegX[NumOfSegBH2]      = {35./2., 10./2., 7./2., 7./2., 7./2., 7./2., 10./2., 35./2.};
const Double_t Bh2SegXAcc[NumOfSegBH2]   = {20., 6.5, 5., 5., 5., 5., 6.5, 20.};
const Double_t localPosBh2X_dX           = 0.;
const Double_t localPosBh2X[NumOfSegBH2] = {-41.5 + localPosBh2X_dX,
  -19.0 + localPosBh2X_dX,
  -10.5 + localPosBh2X_dX,
  -3.5  + localPosBh2X_dX,
  3.5   + localPosBh2X_dX,
  10.5  + localPosBh2X_dX,
  19.0  + localPosBh2X_dX,
  41.5  + localPosBh2X_dX};

const Double_t zTgtTPC = -143.;

//_____________________________________________________________________________
// Local Functions

//_____________________________________________________________________________
inline void
CalcTracks(std::vector<DCLocalTrack*>& trackCont)
{
  for(auto& track: trackCont) track->Calculate();
}

//_____________________________________________________________________________
inline void
CalcTracksTPC(std::vector<TPCLocalTrack*>& trackCont)
{
  for(auto& track: trackCont) track->Calculate();
}

//_____________________________________________________________________________
inline void
CalcTracksTPC_Helix(std::vector<TPCLocalTrack_Helix*>& trackCont)
{
  for(auto& track: trackCont) track->Calculate();
}

//_____________________________________________________________________________
inline void
ClearFlags(std::vector<DCLocalTrack*>& trackCont)
{
  for(const auto& track: trackCont){
    if(!track) continue;
    Int_t nh = track->GetNHit();
    for(Int_t j=0; j<nh; ++j) track->GetHit(j)->QuitTrack();
  }
}

//_____________________________________________________________________________
inline void
DeleteDuplicatedTracks(std::vector<DCLocalTrack*>& trackCont, Double_t ChisqrCut=0.)
{
  // evaluate container size in every iteration
  for(std::size_t i=0; i<trackCont.size(); ++i){
    const auto& tp = trackCont[i];
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();
    for(std::size_t i2=trackCont.size()-1; i2>i; --i2){
      const DCLocalTrack* tp2 = trackCont[i2];
      Int_t nh2 = tp2->GetNHit(), flag=0;
      Double_t chisqr = tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag>0 && chisqr>ChisqrCut){
        delete tp2;
        tp2 = 0;
        trackCont.erase(trackCont.begin()+i2);
      }
    }
  }
}

//_____________________________________________________________________________
inline void
DeleteDuplicatedTracks(std::vector<DCLocalTrack*>& trackCont,
                       Int_t first, Int_t second, Double_t ChisqrCut=0.)
{
  std::vector<Int_t> delete_index;
  // evaluate container size in every iteration
  for(std::size_t i=first; i<=second; ++i){

    auto itr = std::find(delete_index.begin(), delete_index.end(), i);
    if(itr != delete_index.end())
      continue;

    const DCLocalTrack* const tp = trackCont[i];
    if(!tp) continue;

    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->JoinTrack();

    for(std::size_t i2=second; i2>i; --i2){
      auto itr = std::find(delete_index.begin(), delete_index.end(), i2);
      if(itr != delete_index.end())
        continue;

      const DCLocalTrack* tp2 = trackCont[i2];
      Int_t nh2 = tp2->GetNHit(), flag=0;
      Double_t chisqr = tp2->GetChiSquare();
      for(Int_t j=0; j<nh2; ++j)
        if(tp2->GetHit(j)->BelongToTrack()) ++flag;
      if(flag>0 && chisqr>ChisqrCut){
        delete tp2;
        tp2 = 0;

        delete_index.push_back(i2);
      }
    }
  }

  // sort from bigger order
  std::sort(delete_index.begin(), delete_index.end(), std::greater<Int_t>());
  for(Int_t i=0; i<delete_index.size(); i++) {
    trackCont.erase(trackCont.begin()+delete_index[i]);
  }

  // reset hit record of DCHit
  for(std::size_t i=0; i<trackCont.size(); ++i){
    const DCLocalTrack* const tp = trackCont[i];
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    for(Int_t j=0; j<nh; ++j) tp->GetHit(j)->QuitTrack();
  }
}

//_____________________________________________________________________________
inline void // for SSD PreTrack
DeleteWideTracks(std::vector<DCLocalTrack*>& TrackContX,
                 std::vector<DCLocalTrack*>& TrackContY,
                 const std::size_t Nth=1)
{
  const Double_t xSize = 400./2*1.; // SDC1 X-Size
  const Double_t ySize = 250./2*1.; // SDC1 Y-Size
  const Double_t zPos  = -658.;  // SDC1 Z-Position "V1"
  // X
  if(TrackContX.size()>Nth){
    for(Int_t i=TrackContX.size()-1; i>=0; --i){
      const DCLocalTrack* tp = TrackContX[i];
      if(!tp) continue;
      Double_t x = tp->GetVXU(zPos);
      if(std::abs(x)>xSize){
        delete tp;
        tp = 0;
        TrackContX.erase(TrackContX.begin()+i);
      }
    }
  }
  // Y
  if(TrackContY.size()>Nth){
    for(Int_t i=TrackContY.size()-1; i>=0; --i){
      const DCLocalTrack* tp = TrackContY[i];
      if(!tp) continue;
      Double_t y = tp->GetVXU(zPos);
      if(std::abs(y)>ySize){
        delete tp;
        tp = 0;
        TrackContY.erase(TrackContY.begin()+i);
      }
    }
  }
}

//_____________________________________________________________________________
inline void // for SSD PreTrack
DeleteWideTracks(std::vector<DCLocalTrack*>& TrackContX,
                 const std::vector<DCHitContainer>& SdcInHC,
                 const std::size_t Nth=1)
{
  const Double_t xSize =  100.;    // MaxDiff
  const Double_t zPos  = -635.467; // SDC1 Z-Position "X1/X2"
  // X
  if(TrackContX.size()>Nth){
    for(Int_t i=TrackContX.size()-1; i>=0; --i){
      const DCLocalTrack* tp = TrackContX[i];
      if(!tp) continue;
      Double_t x = tp->GetVXU(zPos);

      Bool_t accept = false;
      for(std::size_t j=0; j<NumOfLayersSDC1; ++j){
        if(j!=2 && j!=3) continue; // only X1/X2
        const std::size_t nh = SdcInHC[j+1].size();
        for(std::size_t k=0; k<nh; ++k){
          const DCHit *hit = SdcInHC[j+1][k];
          if(!hit) continue;
          if(std::abs(x-hit->GetWirePosition()) < xSize){
            accept = true;
          }
        }
      }
      if(!accept){
        delete tp;
        tp = 0;
        TrackContX.erase(TrackContX.begin()+i);
      }
    }
  }

  return;
}

//_____________________________________________________________________________
[[maybe_unused]] void
DebugPrint(const IndexList& nCombi,
           const TString& func_name="",
           const TString& msg="")
{
  Int_t n  =1;
  Int_t nn =1;
  Int_t sum=0;
  hddaq::cout << func_name << ":" ;
  IndexList::const_iterator itr, end = nCombi.end();
  for(itr=nCombi.begin(); itr!=end; ++itr){
    Int_t val = *itr;
    sum += val;
    nn *= (val+1);
    hddaq::cout << " " << val;
    if(val!=0) n *= val;
  }
  if(sum==0)
    n=0;
  hddaq::cout << ": total = " << n << ", " << nn << ", " << std::endl;
  return;
}

//_____________________________________________________________________________
[[maybe_unused]] void
DebugPrint(const std::vector<DCLocalTrack*>& trackCont,
           const TString& arg="")
{
  const Int_t nn = trackCont.size();
  hddaq::cout << arg << " " << nn << std::endl;
  for(Int_t i=0; i<nn; ++i){
    const DCLocalTrack * const track=trackCont[i];
    if(!track) continue;
    Int_t    nh     = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    hddaq::cout << std::setw(4) << i
                << "  #Hits : " << std::setw(2) << nh
                << "  ChiSqr : " << chisqr << std::endl;
  }
  hddaq::cout << std::endl;
}

//_____________________________________________________________________________
[[maybe_unused]] void
DebugPrint(const IndexList& nCombi,
           const std::vector<ClusterList>& CandCont,
           const TString& arg="")
{
  hddaq::cout << arg << " #Hits of each group" << std::endl;
  Int_t np = nCombi.size();
  Int_t nn = 1;
  for(Int_t i=0; i<np; ++i){
    hddaq::cout << std::setw(4) << nCombi[i];
    nn *= nCombi[i] + 1;
  }
  hddaq::cout << " -> " << nn-1 << " Combinations" << std::endl;
  for(Int_t i=0; i<np; ++i){
    Int_t n=CandCont[i].size();
    hddaq::cout << "[" << std::setw(3) << i << "]: "
                << std::setw(3) << n << " ";
    for(Int_t j=0; j<n; ++j){
      hddaq::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire()
                  << "(" << CandCont[i][j]->NumberOfHits() << ")"
                  << " ";
    }
    hddaq::cout << std::endl;
  }
}

//_____________________________________________________________________________
template <class Functor>
inline void
FinalizeTrack(const TString& arg,
              std::vector<DCLocalTrack*>& trackCont,
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
                     trackCont.begin() +  index_pair_vec[i].second + 1, DCLTrackComp_Chisqr());
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
// MakeCluster _________________________________________________________________

//_____________________________________________________________________________
Bool_t
MakePairPlaneHitCluster(const DCHitContainer & HC1,
                        const DCHitContainer & HC2,
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
        Int_t multi1 = hit1->GetDriftLengthSize();
        Int_t multi2 = hit2->GetDriftLengthSize();
        for(Int_t m1=0; m1<multi1; ++m1) {
          if(!hit1->IsWithinRange(m1))
            continue;
          for(Int_t m2=0; m2<multi2; ++m2) {
            if(!hit2->IsWithinRange(m2))
              continue;
            Double_t x1,x2;
            if(wp1<wp2){
              x1=wp1+hit1->GetDriftLength(m1);
              x2=wp2-hit2->GetDriftLength(m2);
            }
            else {
              x1=wp1-hit1->GetDriftLength(m1);
              x2=wp2+hit2->GetDriftLength(m2);
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
      Int_t multi1 = hit1->GetDriftLengthSize();
      for(Int_t m1=0; m1<multi1; m1++) {
        if(!(hit1->IsWithinRange(m1))) continue;
        Double_t dl=hit1->GetDriftLength(m1);
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
      Int_t multi2 = hit2->GetDriftLengthSize();
      for(Int_t m2=0; m2<multi2; m2++) {
        if(!(hit2->IsWithinRange(m2))) continue;
        Double_t wp=hit2->GetWirePosition();
        Double_t dl=hit2->GetDriftLength(m2);
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
MakeUnPairPlaneHitCluster(const DCHitContainer& HC,
                          ClusterList& Cont,
                          Bool_t honeycomb=false)
{
  const std::size_t nh = HC.size();
  for(std::size_t i=0; i<nh; ++i){
    DCHit *hit = HC[i];
    if(!hit) continue;
    std::size_t mh = hit->GetDriftLengthSize();
    for(std::size_t m=0; m<mh; ++m) {
      if(!hit->IsWithinRange(m)) continue;
      Double_t wp = hit->GetWirePosition();
      Double_t dl = hit->GetDriftLength(m);
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
MakeMWPCPairPlaneHitCluster(const DCHitContainer& HC,
                            ClusterList& Cont)
{
  Int_t nh=HC.size();
  for(Int_t i=0; i<nh; ++i){
    DCHit *hit=HC[i];
    if(hit){
      Int_t multi = hit->GetDriftTimeSize();
      for(Int_t m=0; m<multi; m++) {
        if(!(hit->IsWithinRange(m))) continue;
        Double_t wp=hit->GetWirePosition();
        // Double_t dl=hit->GetDriftLength(m);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit, wp, m)));
      }
    }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
MakeTOFHitCluster(const DCHitContainer& HitCont,
                  ClusterList& Cont,
                  Int_t xy)
{
  Int_t nh = HitCont.size();
  for(Int_t i=0; i<nh; ++i){
    if(i%2!=xy) continue;
    DCHit *hit = HitCont[i];
    if(!hit) continue;
    Double_t wp = hit->GetWirePosition();
    Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit, wp, 0)));
  }
  return true;
}

//  //_____________________________________________________________________________
//  Bool_t
//  MakeTPCHitCluster(const TPCHitContainer& HitCont,
//         TPCClusterList& ClCont,
//        Int_t layerID)
//  {
//    static const Double_t ClusterYCut = gUser.GetParameter("ClusterYCut");
//
//    //del::ClearContainer(ClCont);
//
//    const std::size_t nh = HitCont.size();
//    if(nh==0) return false;
//
//    std::vector<Int_t> flag(nh, 0);
//
//    for(std::size_t hiti=0; hiti < nh; hiti++) {
//      if(flag[hiti] > 0) continue;
//      TPCHitContainer CandCont;
//      TPCHit* hit = HitCont[hiti];
//      if(!hit || !hit->IsGoodHit()) continue;
//      CandCont.push_back(hit);
//      flag[hiti]++;
//
//      for(std::size_t hitj=0; hitj < nh; hitj++) {
//        if(hiti==hitj || flag[hitj]>0) continue;
//        TPCHit* thit = HitCont[hitj];
//        if(!thit || !thit->IsGoodHit()) continue;
//        for(Int_t ci=0; ci < CandCont.size(); ci++) {
//          TPCHit* c_hit = CandCont[ci];
//          Int_t rowID = thit->RowId();
//          Int_t c_rowID = c_hit->RowId();
//          if((abs(rowID - c_rowID) <= 2 ||
//         (layerID<10 && abs(rowID - c_rowID)>=tpc::padParameter[layerID][1]-2))
//       && fabs(thit->Y() - c_hit->Y()) < ClusterYCut)
//          {
//            CandCont.push_back(thit);
//            flag[hitj]++;
//            break;
//          }
//        }
//      }
//      ClCont.push_back(new TPCCluster(layerID, CandCont));
//    }
//
//    return true;
//  }

//_____________________________________________________________________________
Bool_t
MakePairPlaneHitClusterVUX(const DCHitContainer& HC1,
                           const DCHitContainer& HC2,
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

        Int_t multi1 = hit1->GetDriftLengthSize();
        Int_t multi2 = hit2->GetDriftLengthSize();
        for(Int_t m1=0; m1<multi1; m1++) {
          if(!(hit1->IsWithinRange(m1))) continue;
          for(Int_t m2=0; m2<multi2; m2++) {
            if(!(hit2->IsWithinRange(m2))) continue;
            Double_t dl1=hit1->GetDriftLength(m1);
            Double_t dl2=hit2->GetDriftLength(m2);

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
      Int_t multi1 = hit1->GetDriftLengthSize();
      for(Int_t m1=0; m1<multi1; m1++) {
        if(!(hit1->IsWithinRange(m1))) continue;
        Double_t dl=hit1->GetDriftLength(m1);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1+dl,m1)));
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit1,wp1-dl,m1)));
      }
    }
  }
  for(Int_t i2=0; i2<nh2; ++i2){
    if(UsedFlag[i2]==0) {
      DCHit *hit2=HC2[i2];
      Int_t multi2 = hit2->GetDriftLengthSize();
      for(Int_t m2=0; m2<multi2; m2++) {
        if(!(hit2->IsWithinRange(m2))) continue;

        Double_t wp=hit2->GetWirePosition();
        Double_t dl=hit2->GetDriftLength(m2);
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit2,wp+dl,m2)));
        Cont.push_back(new DCPairHitCluster(new DCLTrackHit(hit2,wp-dl,m2)));
      }
    }
  }

  return true;
}

}


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
#if 1
        hddaq::cout << FUNC_NAME << " too much combinations... " << n2 << std::endl;
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
      // Int_t size1=index.size();
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
Int_t /* Local Track Search without BH2Filter */
LocalTrackSearch(const std::vector<DCHitContainer>& HC,
                 const DCPairPlaneInfo * PpInfo,
                 Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
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

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits
       && track->DoFit()
       && track->GetChiSquare()<MaxChisquare){

      if(T0Seg>=0 && T0Seg<NumOfSegBH2) {
        Double_t xbh2=track->GetX(zBH2), ybh2=track->GetY(zBH2);
        Double_t difPosBh2 = localPosBh2X[T0Seg] - xbh2;

        //   Double_t xtgt=track->GetX(zTarget), ytgt=track->GetY(zTarget);
        //   Double_t ytgt=track->GetY(zTarget);

        if(true
           && fabs(difPosBh2)<Bh2SegXAcc[T0Seg]
           && (-10 < ybh2 && ybh2 < 40)
           //       && fabs(ytgt)<21.
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
  return status? TrackCont.size() : -1;
}

//_____________________________________________________________________________
Int_t /* Local Track Search with BH2Filter */
LocalTrackSearch(const std::vector<std::vector<DCHitContainer>> &hcAssemble,
                 const DCPairPlaneInfo * PpInfo,
                 Int_t npp, std::vector<DCLocalTrack*> &trackCont,
                 Int_t MinNumOfHits, Int_t T0Seg)
{
  std::vector<std::vector<DCHitContainer>>::const_iterator
    itr, itr_end = hcAssemble.end();

  Int_t status = 0;
  for(itr=hcAssemble.begin(); itr!=itr_end; ++itr){
    const std::vector<DCHitContainer>& l = *itr;
    std::vector<DCLocalTrack*> tc;
    status = LocalTrackSearch(l, PpInfo, npp, tc, MinNumOfHits, T0Seg);
    trackCont.insert(trackCont.end(), tc.begin(), tc.end());
  }

  ClearFlags(trackCont);
  std::stable_sort(trackCont.begin(), trackCont.end(), DCLTrackComp());

  DeleteDuplicatedTracks(trackCont);
  //    CalcTracks(trackCont);

  return status < 0? status : trackCont.size();
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchSdcOut(const std::vector<DCHitContainer>& SdcOutHC,
                       const DCPairPlaneInfo *PpInfo,
                       Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                       Int_t MinNumOfHits /*=6*/)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    if(ppFlag ){
      MakePairPlaneHitCluster(SdcOutHC[layer1], SdcOutHC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else{
      MakeUnPairPlaneHitCluster(SdcOutHC[layer1], CandCont[i], honeycomb);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i) {
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits     &&
       track->GetNHitY() >= 2             &&
       track->DoFit()                     &&
       track->GetChiSquare()<MaxChisquare)
    {
      TrackCont.push_back(track);
    }
    else
    {
      delete track;
    }
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status? TrackCont.size() : -1;
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchSdcOut(const DCHitContainer& TOFHC,
                       const std::vector<DCHitContainer>& SdcOutHC,
                       const DCPairPlaneInfo *PpInfo,
                       Int_t npp,
                       std::vector<DCLocalTrack*>& TrackCont,
                       Int_t MinNumOfHits /*=6*/)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp-2; ++i){
    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;
    if(ppFlag ){ //DC2, 3
      MakePairPlaneHitCluster(SdcOutHC[layer1], SdcOutHC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else{ //FBT
      MakeMWPCPairPlaneHitCluster(SdcOutHC[layer1], CandCont[i]);
      MakeMWPCPairPlaneHitCluster(SdcOutHC[layer2], CandCont[i]);
    }
  }

  // TOF
  MakeTOFHitCluster(TOFHC, CandCont[npp-2], 0);
  MakeTOFHitCluster(TOFHC, CandCont[npp-1], 1);

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;

    static const Int_t IdTOF_UX = gGeom.GetDetectorId("TOF-UX");
    static const Int_t IdTOF_UY = gGeom.GetDetectorId("TOF-UY");
    static const Int_t IdTOF_DX = gGeom.GetDetectorId("TOF-DX");
    static const Int_t IdTOF_DY = gGeom.GetDetectorId("TOF-DY");

    Bool_t TOFSegXYMatching =
      (track->GetWire(IdTOF_UX)==track->GetWire(IdTOF_UY)) ||
      (track->GetWire(IdTOF_DX)==track->GetWire(IdTOF_DY));

    // Int_t Track[20]={0};
    // Int_t layer;
    // for(Int_t i=0; i<(track->GetNHit()); ++i){
    //  layer=track->GetHit(i)->GetLayer();
    //  Track[layer]=1;
    // }

    // Bool_t FBT =
    //  (Track[80]==1 && Track[82]==1) || (Track[81]==1 && Track[83]==1) ||
    //  (Track[84]==1 && Track[86]==1) || (Track[85]==1 && Track[87]==1) ;
    // Bool_t DC23x_off =
    //  (Track[31]==0 && Track[32]==0 && Track[37]==0 && Track[38]==0);


    if(TOFSegXYMatching &&
       //FBT&&
       track->GetNHit()>=MinNumOfHits+2   &&
       track->GetNHitY() >= 2             &&
       track->DoFit()                     &&
       track->GetChiSquare()<MaxChisquare)
    {
      TrackCont.push_back(track);
    }
    else
    {
      delete track;
    }
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status? TrackCont.size() : -1;
}

//_____________________________________________________________________________
Int_t /* Local Track Search SdcIn w/Fiber */
LocalTrackSearchSdcInFiber(const std::vector<DCHitContainer>& HC,
                           const DCPairPlaneInfo *PpInfo,
                           Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                           Int_t MinNumOfHits /*=6*/)
{
  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<npp; ++i){

    Bool_t ppFlag    = PpInfo[i].pair;
    Bool_t honeycomb = PpInfo[i].honeycomb;
    Int_t  layer1    = PpInfo[i].id1;
    Int_t  layer2    = PpInfo[i].id2;

    if(ppFlag) {
      MakePairPlaneHitCluster(HC[layer1], HC[layer2],
                              PpInfo[i].CellSize, CandCont[i], honeycomb);
    }else if(layer1==layer2){
      MakeMWPCPairPlaneHitCluster(HC[layer1], CandCont[i]);
    }else{
      MakeUnPairPlaneHitCluster(HC[layer1], CandCont[i], honeycomb);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i) {
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);

  for(Int_t i=0, n=CombiIndex.size(); i<n; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(true
       && track->GetNHitSFT() > 1
       && track->GetNHit()>=MinNumOfHits
       && track->DoFit()
       && track->GetChiSquare()<MaxChisquare
      ){
      TrackCont.push_back(track);
    }
    else
      delete track;
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackCompSdcInFiber(), CandCont);
  return status? TrackCont.size() : -1;
}

// BC3&4, SDC1 VUX Tracking ___________________________________________
Int_t
LocalTrackSearchVUX(const std::vector<DCHitContainer>& HC,
                    const DCPairPlaneInfo *PpInfo,
                    Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                    Int_t MinNumOfHits /*=6*/)
{
  std::vector<DCLocalTrack*> TrackContV;
  std::vector<DCLocalTrack*> TrackContX;
  std::vector<DCLocalTrack*> TrackContU;
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

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackCompSdcInFiber(), CandCont);

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

//_____________________________________________________________________________
Int_t
LocalTrackSearchBcOutSdcIn(const std::vector<DCHitContainer>& BcHC,
                           const DCPairPlaneInfo *BcPpInfo,
                           const std::vector<DCHitContainer>& SdcHC,
                           const DCPairPlaneInfo *SdcPpInfo,
                           Int_t BcNpp, Int_t SdcNpp,
                           std::vector<DCLocalTrack*>& TrackCont,
                           Int_t MinNumOfHits)
{
  const Int_t npp = BcNpp + SdcNpp;

  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<BcNpp; ++i){
    Bool_t ppFlag=BcPpInfo[i].pair;
    Int_t layer1=BcPpInfo[i].id1, layer2=BcPpInfo[i].id2;
    if(ppFlag)
      MakePairPlaneHitCluster(BcHC[layer1], BcHC[layer2],
                              BcPpInfo[i].CellSize, CandCont[i]);
    else
      MakeUnPairPlaneHitCluster(BcHC[layer1], CandCont[i]);
  }

  for(Int_t i=0; i<SdcNpp; ++i){
    Bool_t ppFlag=SdcPpInfo[i].pair;
    Int_t layer1=SdcPpInfo[i].id1, layer2=SdcPpInfo[i].id2;
    if(ppFlag)
      MakePairPlaneHitCluster(SdcHC[layer1], SdcHC[layer2],
                              SdcPpInfo[i].CellSize, CandCont[i+BcNpp]);
    else
      MakeUnPairPlaneHitCluster(SdcHC[layer1], CandCont[i+BcNpp]);
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
  Int_t nnCombi = CombiIndex.size();

#if 0
  hddaq::cout << " ===> " << nnCombi << " combinations will be checked.."
              << std::endl;
  for(Int_t i=0; i<nnCombi; ++i){
    for(Int_t j=0;j<npp;j++) {
      hddaq::cout << CombiIndex[i][j] << " ";
    }
    hddaq::cout << std::endl;
  }
#endif

  for(Int_t i=0; i<nnCombi; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits && track->DoFitBcSdc() &&
       track->GetChiSquare()<MaxChisquare)
      TrackCont.push_back(track);
    else
      delete track;
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);

  return TrackCont.size();
}

//_____________________________________________________________________________
inline Bool_t
IsDeletionTarget(const std::vector<std::pair<Int_t,Int_t>>& nh,
                 std::size_t NDelete, Int_t layer)
{
  if(NDelete==0) return false;

  std::vector<std::pair<Int_t,Int_t>>::const_iterator itr, end;
  for(itr=nh.begin(), end=nh.end(); itr!=end; ++itr){
    if(itr->first==0 && itr->second==layer)
      return true;
  }

  if(nh[0].second == layer)
    return true;

  return false;
}

//_____________________________________________________________________________
Int_t
LocalTrackSearchSdcInSdcOut(const std::vector<DCHitContainer>& SdcInHC,
                            const DCPairPlaneInfo* SdcInPpInfo,
                            const std::vector<DCHitContainer>& SdcOutHC,
                            const DCPairPlaneInfo* SdcOutPpInfo,
                            Int_t SdcInNpp, Int_t SdcOutNpp,
                            std::vector<DCLocalTrack*>& TrackCont,
                            Int_t MinNumOfHits)
{
  const Int_t npp = SdcInNpp+SdcOutNpp;

  std::vector<ClusterList> CandCont(npp);

  for(Int_t i=0; i<SdcInNpp; ++i){
    Bool_t ppFlag = SdcInPpInfo[i].pair;
    Int_t  layer1 = SdcInPpInfo[i].id1;
    Int_t  layer2 = SdcInPpInfo[i].id2;
    if(ppFlag){
      MakePairPlaneHitCluster(SdcInHC[layer1], SdcInHC[layer2],
                              SdcInPpInfo[i].CellSize, CandCont[i]);
    } else {
      MakeUnPairPlaneHitCluster(SdcInHC[layer1], CandCont[i]);
    }
  }

  for(Int_t i=0; i<SdcOutNpp; ++i){
    Bool_t ppFlag = SdcOutPpInfo[i].pair;
    Int_t  layer1 = SdcOutPpInfo[i].id1;
    Int_t  layer2 = SdcOutPpInfo[i].id2;
    if(ppFlag){
      MakePairPlaneHitCluster(SdcOutHC[layer1], SdcOutHC[layer2],
                              SdcOutPpInfo[i].CellSize, CandCont[i+SdcInNpp]);
    } else {
      MakeUnPairPlaneHitCluster(SdcOutHC[layer1], CandCont[i+SdcInNpp]);
    }
  }

  IndexList nCombi(npp);
  for(Int_t i=0; i<npp; ++i){
    Int_t n = CandCont[i].size();
    nCombi[i] = n>MaxNumOfCluster ? 0 : n;
  }

  Bool_t status = true;
  std::vector<IndexList> CombiIndex = MakeIndex(npp, nCombi, status);
  Int_t nnCombi = CombiIndex.size();

  for(Int_t i=0; i<nnCombi; ++i){
    DCLocalTrack *track = MakeTrack(CandCont, CombiIndex[i]);
    if(!track) continue;
    if(track->GetNHit()>=MinNumOfHits &&
       track->DoFit() &&
       track->GetChiSquare()<MaxChisquare)
      TrackCont.push_back(track);
    else
      delete track;
  }

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return status? TrackCont.size() : -1;
}


//_____________________________________________________________________________
Int_t
LocalTrackSearchTPC(const std::vector<TPCHitContainer>& TPCHC,
                    std::vector<TPCLocalTrack*>& TrackCont,
                    Int_t MinNumOfHits /*=8*/)
{
  static const Double_t HoughWindowCut = gUser.GetParameter("HoughWindowCut");
  Bool_t status = true;

  //    if(valueHall) { // TODO
  //    }

  // y = p0 + p1 * x
  Double_t p0[MaxNumOfTrackTPC];
  Double_t p1[MaxNumOfTrackTPC];
  const Int_t    Li_theta_ndiv = 200;
  const Double_t Li_theta_min  =   0;
  const Double_t Li_theta_max  = 180;
  const Int_t    Li_r_ndiv =  200;
  const Double_t Li_r_min  = -600;
  const Double_t Li_r_max  =  600;

  static const Int_t ClusterSizeCut = gUser.GetParameter("TPCClusterSizeCut");


  //for TPC linear track
  // r = x * cos(theta) + y * sin(theta)
  TH2D Li_hist("hist_linear",";theta (deg.); r (mm)",
               Li_theta_ndiv, Li_theta_min, Li_theta_max,
               Li_r_ndiv, Li_r_min, Li_r_max);

  std::vector<std::vector<Int_t>> flag;
  flag.resize(NumOfLayersTPC);
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    flag[layer].resize(TPCHC[layer].size(), 0);
  }

  std::vector<Double_t> hough_x;
  std::vector<Double_t> hough_y;

  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    Li_hist.Reset();
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=TPCHC[layer].size(); ci<n; ci++){
        if(flag[layer][ci]>0) continue;
        TPCHit* hit = TPCHC[layer][ci];
        TVector3 pos = hit->GetPos();
        for(Int_t ti=0; ti<Li_theta_ndiv; ti++){
          Double_t theta = Li_theta_min+ti*(Li_theta_max-Li_theta_min)/Li_theta_ndiv;
          Li_hist.Fill(theta, cos(theta*acos(-1)/180.)*pos.Z()
                       +sin(theta*acos(-1)/180.)*pos.X());
          if(fabs(cos(theta*acos(-1)/180.)*pos.Z()+sin(theta*acos(-1)/180.)*pos.X())>Li_r_max)
            std::cout<<"Hough: out of range:"<<cos(theta*acos(-1)/180.)*pos.Z()+sin(theta*acos(-1)/180.)*pos.X()<<std::endl;
        }
      } // cluster
    } // layer
      //      if(Li_hist.GetMaximum() < MinNumOfHits){
    if(Li_hist.GetMaximum() < MinNumOfHits/2){
      //Li_hist.Delete();
      Li_hist.Reset();
      break;
    }




    Int_t maxbin = Li_hist.GetMaximumBin();
    Int_t mx,my,mz;
    Li_hist.GetBinXYZ(maxbin, mx, my, mz);
    Double_t mtheta = Li_hist.GetXaxis()->GetBinCenter(mx)*acos(-1)/180.;
    Double_t mr = Li_hist.GetYaxis()->GetBinCenter(my);

    Bool_t hough_flag = true;
    for(Int_t i=0; i<hough_x.size(); ++i){
      Int_t bindiff = fabs(mx-hough_x[i])+fabs(my-hough_y[i]);
      if(bindiff<=4)
        hough_flag = false;
    }
    hough_x.push_back(mx);
    hough_y.push_back(my);
    if(!hough_flag)
      continue;


    TPCLocalTrack *track = new TPCLocalTrack();
    p0[tracki] = mr/sin(mtheta);
    p1[tracki] = -cos(mtheta)/sin(mtheta);


    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=TPCHC[layer].size(); ci<n; ci++){
        //if(flag[layer][ci]>0) continue;
        TPCHit* hit = TPCHC[layer][ci];
        TVector3 pos = hit->GetPos();
        Double_t dist = fabs(p1[tracki]*pos.Z()-pos.X()+p0[tracki])/sqrt(pow(p1[tracki],2)+1);

        if(dist < HoughWindowCut && hit->GetClusterSize()>=ClusterSizeCut){
          track->AddTPCHit(new TPCLTrackHit(hit));
          //     track_he->AddTPCHit(new TPCLTrackHit(hit));
          flag[layer][ci]++;
        }
      }
    }
    //temporary (x0, y0) are position at Target position
    //Double_t zTgtTPC = -143.;
    track->SetAx(p0[tracki]+p1[tracki]*zTgtTPC);
    track->SetAu(p1[tracki]);

    if(track->DoFit(MinNumOfHits)){
      TrackCont.push_back(track);
    }
    else
      delete track;

    Li_hist.Reset();
    //delete track;

  }//track

  CalcTracksTPC(TrackCont);
  // std::cout<<"event end"<<std::endl;
  //getchar();

  return status? TrackCont.size() : -1;

  return 0;
}


//_____________________________________________________________________________
Int_t
LocalTrackSearchTPC_Helix(const std::vector<TPCHitContainer>& TPCHC,
                          std::vector<TPCLocalTrack_Helix*>& TrackCont,
                          Int_t MinNumOfHits /*=8*/)
{
  static const Double_t HoughWindowCut = gUser.GetParameter("HoughWindowCut");
  static const Double_t MaxLayerCut = gUser.GetParameter("TPCMaxLayerCut");
  static const Double_t DECut_TPCTrack = gUser.GetParameter("DECut_TPCTrack");
  Bool_t status = true;

  //    if(valueHall) { // TODO
  //    }
  const Double_t Const = 0.299792458;
  //(x - (r + rd)*cos(theta))^2 + (y - (r + rd)*sin(theta))^2 = r^2
  // p = r * Const; // 1T
  // Parameters
  Double_t hough_p[MaxNumOfTrackTPC];
  Double_t hough_theta[MaxNumOfTrackTPC];
  Double_t hough_rd[MaxNumOfTrackTPC];

  const Int_t nBin_rdiff = 11;
  const Double_t rdiff_min = -110.;
  const Double_t rdiff_max = 110.;

  const Int_t nBin_theta = 180;
  const Double_t theta_min = -1.*acos(-1);
  const Double_t theta_max = acos(-1);

  const Int_t nBin_p = 300;
  const Double_t pmin = 50.;//MeV/c
  const Double_t pmax = 2050.;//MeV/c

  //for TPC circle track
  // static TH3D Ci_hist("hist_circle",";rd (mm); theta (rad); p(MeV/c)",
  //    nBin_rdiff, rdiff_min,  rdiff_max,
  //    nBin_theta, theta_min, theta_max,
  //    nBin_p,pmin,pmax);

  TH3D *Ci_hist=new TH3D("hist_circle",";rd (mm); theta (rad); p(MeV/c)",
                         nBin_rdiff, rdiff_min,  rdiff_max,
                         nBin_theta, theta_min, theta_max,
                         nBin_p,pmin,pmax);

  // std::vector<std::vector<Int_t>> flag;
  // flag.resize(NumOfLayersTPC);
  // for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
  //   flag[layer].resize(TPCHC[layer].size(), 0);
  // }

  std::vector<Double_t> hough_x;
  std::vector<Double_t> hough_y;

  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    Ci_hist->Reset();
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=TPCHC[layer].size(); ci<n; ci++){
	//        if(flag[layer][ci]>0) continue;
        TPCHit* hit = TPCHC[layer][ci];
	if(hit->GetHoughFlag()>0) continue;

        TVector3 pos = hit->GetPos();
        for(Int_t ird=0; ird<nBin_rdiff; ++ird){
          Double_t rd = Ci_hist->GetXaxis()->GetBinCenter(ird+1);
          for(Int_t ip=0; ip<nBin_p; ++ip){
            Double_t x = -pos.x();
            Double_t y = pos.z()-zTgtTPC;
            Double_t p = Ci_hist->GetZaxis()->GetBinCenter(ip+1);
            Double_t r = p/(Const*1.);//1T

            //a*sin(theta) + b*cos(theta) +c = 0
            Double_t a = 2.*(r+rd)*y;
            Double_t b = 2.*(r+rd)*x;
            Double_t c = -1.*(rd*rd + 2.*r*rd + x*x + y*y);

            Double_t r0 = sqrt(a*a + b*b);
            if(fabs(-1.*c/r0)>1.){
              // std::cout<<"No solution, "
              //   <<"x:"<<x<<", y:"<<y
              //   <<", r:"<<r<<", rd:"<<rd<<std::endl;
              continue;
            }
            Double_t theta1_alpha =  asin(-1.*c/r0);
            Double_t theta2_alpha;
            if(theta1_alpha>0.)
              theta2_alpha = acos(-1.) - theta1_alpha;
            else
              theta2_alpha = -1.*acos(-1.) - theta1_alpha;

            Double_t theta_alpha = atan2(b, a);

            Double_t xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
            Double_t ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
            Double_t r_re1 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

            Double_t xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
            Double_t ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
            Double_t r_re2 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

            Double_t theta1 = atan2(ycenter1, xcenter1);
            Double_t theta2 = atan2(ycenter2, xcenter2);

            if(TMath::IsNaN(theta1)){
              std::cout<<"theta1="<<theta1<<", x="<<x<<", y="<<y
                       <<"rd="<<rd<<", r"<<r<<std::endl;
            }


            if(fabs(r-r_re1)>0.01||fabs(r-r_re2)>0.01){
              std::cout<<"r="<<r<<", r_re1="<<r_re1<<", r_re1="<<r_re2<<std::endl;
              std::cout<<"x:"<<x<<", y:"<<y
                       <<", theta1:"<<theta1<<", theta2:"<<theta2
                       <<", theta_alpha:"<<theta_alpha<<std::endl;
            }
            Ci_hist->Fill(rd, theta1, p);
            Ci_hist->Fill(rd, theta2, p);
            // std::cout<<"rd: "<<rd<<", "
            //         <<"theta1: "<<theta1<<", "
            //         <<"theta2: "<<theta2<<", "
            //         <<"p: "<<p<<std::endl;
          }
        }
      }// cluster
    } // layer

    if(Ci_hist->GetMaximum() < MinNumOfHits/2){
      Ci_hist->Reset();
      break;
    }
    //std::cout<<"Maxbin0: "<<Ci_hist.GetMaximum()<<std::endl;


    Int_t maxbin = Ci_hist->GetMaximumBin();
    Int_t mx,my,mz;
    Ci_hist->GetBinXYZ(maxbin, mx, my, mz);

    Bool_t hough_flag = true;
    for(Int_t i=0; i<hough_x.size(); ++i){
      Int_t bindiff = fabs(mx-hough_x[i])+fabs(my-hough_y[i]);
      if(bindiff<=4)
      //if(bindiff<=8)
        hough_flag = false;
    }
    hough_x.push_back(mx);
    hough_y.push_back(my);
    if(!hough_flag)
      continue;

    TPCLocalTrack_Helix *track = new TPCLocalTrack_Helix();
    hough_rd[tracki] = Ci_hist->GetXaxis()->GetBinCenter(mx);
    hough_theta[tracki] = Ci_hist->GetYaxis()->GetBinCenter(my);
    hough_p[tracki] = Ci_hist->GetZaxis()->GetBinCenter(mz);
    Double_t hough_r = hough_p[tracki]/Const;
    Double_t hough_cx = (hough_r + hough_rd[tracki])*cos(hough_theta[tracki]);
    Double_t hough_cy = (hough_r + hough_rd[tracki])*sin(hough_theta[tracki]);


    //std::cout<<"Hough (x,z) Maxbin: "<<Ci_hist.GetMaximum()<<std::endl;
    //     std::cout<<""<<std::endl;
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=TPCHC[layer].size(); ci<n; ci++){
        //if(flag[layer][ci]>0) continue;
        TPCHit* hit = TPCHC[layer][ci];
        TVector3 pos = hit->GetPos();
        Double_t x = -pos.x();
        Double_t y = pos.z()-zTgtTPC;
	Double_t de = hit->GetCharge();
        // Double_t xcenter1 = (hough_r + hough_rd[tracki])*cos(hough_theta[tracki]);
        // Double_t ycenter1 = (hough_r + hough_rd[tracki])*sin(hough_theta[tracki]);
        Double_t r_cal = sqrt(pow(x-hough_cx,2) + pow(y-hough_cy,2));

        Double_t dist = fabs(r_cal - hough_r);
        if(dist < HoughWindowCut && layer < MaxLayerCut && de>DECut_TPCTrack){
          track->AddTPCHit(new TPCLTrackHit(hit));
          //flag[layer][ci]++;
        }
      }
    }
    // track->SetAdrho(hough_rd[tracki]);
    // track->SetAphi0(hough_theta[tracki]);
    // track->SetArho(1./hough_r);
    track->SetAcx(hough_cx);
    track->SetAcy(hough_cy);
    track->SetAr(hough_r);

    if(track->DoFit(MinNumOfHits)){
      TrackCont.push_back(track);
      Int_t nh = track->GetNHit();
      for( int ih=0; ih<nh; ++ih ){
	TPCHit *hit = track->GetHit( ih )->GetHit();
	if( !hit ) continue;
	hit->SetHoughFlag(1+hit->GetHoughFlag());
      }
    }
    else
      delete track;

    //      Ci_hist->Reset();

    //delete track;
  }//track
  CalcTracksTPC_Helix(TrackCont);
  delete Ci_hist;
  return status? TrackCont.size() : -1;

  return 0;
}





//For MWPC
//_____________________________________________________________________________
Int_t MWPCLocalTrackSearch(const std::vector<DCHitContainer>& HC,
                           std::vector<DCLocalTrack*>& TrackCont)

{
  std::vector<ClusterList> CandCont(NumOfLayersBcIn);

  for(Int_t i=0; i<NumOfLayersBcIn; ++i){
    MakeMWPCPairPlaneHitCluster(HC[i], CandCont[i]);
  }

  //   IndexList nCombi(NumOfLayersBcIn);
  //   for(Int_t i=0; i<NumOfLayersBcIn; ++i){
  //     nCombi[i]=(CandCont[i]).size();

  //     // If #Cluster>MaxNumerOfCluster,  error return

  //     if(nCombi[i]>MaxNumOfCluster){
  //       hddaq::cout << FUNC_NAME << " too many clusters " << FUNC_NAME
  //                   << "  layer = " << i << " : " << nCombi[i] << std::endl;
  //       del::ClearContainerAll(CandCont);
  //       return 0;
  //     }
  //   }

  //   DebugPrint(nCombi, FUNC_NAME);

#if 0
  DebugPrint(nCombi, CandCont, FUNC_NAME);
#endif

  const Int_t MinNumOfHitsBcIn   = 6;
  TrackMaker trackMaker(CandCont, MinNumOfHitsBcIn, MaxCombi, MaxChisquare);
  trackMaker.MakeTracks(TrackCont);

  FinalizeTrack(FUNC_NAME, TrackCont, DCLTrackComp(), CandCont);
  return TrackCont.size();
}

//_____________________________________________________________________________
Int_t
MWPCLocalTrackSearch(const std::vector<std::vector<DCHitContainer>>& hcList,
                     std::vector<DCLocalTrack*>& trackCont)
{
  for(auto itr=hcList.begin(), end=hcList.end(); itr!=end; ++itr){
    const std::vector<DCHitContainer>& l = *itr;
    std::vector<DCLocalTrack*> tc;
    MWPCLocalTrackSearch(l, tc);
    trackCont.insert(trackCont.end(), tc.begin(), tc.end());
    // hddaq::cout << " tc " << tc.size()
    //             << " : " << trackCont.size() << std::endl;
  }

  ClearFlags(trackCont);
  DeleteDuplicatedTracks(trackCont);
  CalcTracks(trackCont);
  return trackCont.size();
}
}
