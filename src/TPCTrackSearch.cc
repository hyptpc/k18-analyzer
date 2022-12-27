// -*- C++ -*-

#include "TPCTrackSearch.hh"

#include <chrono>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TH3D.h>

#include "DCGeomMan.hh"
#include "TPCLTrackHit.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCCluster.hh"
#include "TPCBeamRemover.hh"

#include "RootHelper.hh"

#define DebugDisp      0
#define DebugEvDisp    0
#define DummyHitTest   0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
#endif

namespace
{
const std::string& class_name("TPCTrackSearch");
const double& HS_field_0 = 0.9860;
const double& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const double& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& zTarget    = gGeom.LocalZ("Target");
const auto& zK18Target = gGeom.LocalZ("K18Target");
// const auto& zBH2       = gGeom.LocalZ("BH2");
// const Double_t MaxChisquare       = 2000.; // Set to be More than 30
// const Double_t MaxChisquareSdcIn  = 5000.; // Set to be More than 30
// const Double_t MaxNumOfCluster = 20.;    // Set to be Less than 30
// const Double_t MaxCombi = 1.0e6;    // Set to be Less than 10^6
// // SdcIn & BcOut for XUV Tracking routine
// const Double_t MaxChisquareVXU = 50.;//
// const Double_t ChisquareCutVXU = 50.;//

// TPC Tracking
//const Int_t    MaxNumOfTrackTPC = 100;
const Int_t    MaxNumOfTrackTPC = 5;
const auto& valueHall = ConfMan::Get<Double_t>("HSFLDHALL");

//_____________________________________________________________________________
// Local Functions

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
// [[maybe_unused]] void
// DebugPrint(const IndexList& nCombi,
//            const TString& func_name="",
//            const TString& msg="")
// {
//   Int_t n  =1;
//   Int_t nn =1;
//   Int_t sum=0;
//   hddaq::cout << func_name << ":" ;
//   IndexList::const_iterator itr, end = nCombi.end();
//   for(itr=nCombi.begin(); itr!=end; ++itr){
//     Int_t val = *itr;
//     sum += val;
//     nn *= (val+1);
//     hddaq::cout << " " << val;
//     if(val!=0) n *= val;
//   }
//   if(sum==0)
//     n=0;
//   hddaq::cout << ": total = " << n << ", " << nn << ", " << std::endl;
//   return;
// }

//_____________________________________________________________________________
// [[maybe_unused]] void
// DebugPrint(const std::vector<DCLocalTrack*>& trackCont,
//            const TString& arg="")
// {
//   const Int_t nn = trackCont.size();
//   hddaq::cout << arg << " " << nn << std::endl;
//   for(Int_t i=0; i<nn; ++i){
//     const DCLocalTrack * const track=trackCont[i];
//     if(!track) continue;
//     Int_t    nh     = track->GetNHit();
//     Double_t chisqr = track->GetChiSquare();
//     hddaq::cout << std::setw(4) << i
//                 << "  #Hits : " << std::setw(2) << nh
//                 << "  ChiSqr : " << chisqr << std::endl;
//   }
//   hddaq::cout << std::endl;
// }

//_____________________________________________________________________________
// [[maybe_unused]] void
// DebugPrint(const IndexList& nCombi,
//            const std::vector<ClusterList>& CandCont,
//            const TString& arg="")
// {
//   hddaq::cout << arg << " #Hits of each group" << std::endl;
//   Int_t np = nCombi.size();
//   Int_t nn = 1;
//   for(Int_t i=0; i<np; ++i){
//     hddaq::cout << std::setw(4) << nCombi[i];
//     nn *= nCombi[i] + 1;
//   }
//   hddaq::cout << " -> " << nn-1 << " Combinations" << std::endl;
//   for(Int_t i=0; i<np; ++i){
//     Int_t n=CandCont[i].size();
//     hddaq::cout << "[" << std::setw(3) << i << "]: "
//                 << std::setw(3) << n << " ";
//     for(Int_t j=0; j<n; ++j){
//       hddaq::cout << ((DCLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire()
//                   << "(" << CandCont[i][j]->NumberOfHits() << ")"
//                   << " ";
//     }
//     hddaq::cout << std::endl;
//   }
// }

//_____________________________________________________________________________
// template <class Functor>
// inline void
// FinalizeTrack(const TString& arg,
//               std::vector<DCLocalTrack*>& trackCont,
//               Functor comp,
//               std::vector<ClusterList>& candCont,
//               Bool_t delete_flag=true)
// {
//   ClearFlags(trackCont);

// #if 0
//   DebugPrint(trackCont, arg+" Before Sorting ");
// #endif

//   std::stable_sort(trackCont.begin(), trackCont.end(), DCLTrackComp_Nhit());


// #if 0
//   DebugPrint(trackCont, arg+" After Sorting (Nhit) ");
// #endif

//   typedef std::pair<Int_t, Int_t> index_pair;
//   std::vector<index_pair> index_pair_vec;

//   std::vector<Int_t> nhit_vec;

//   for(Int_t i=0; i<trackCont.size(); i++) {
//     Int_t nhit = trackCont[i]->GetNHit();
//     nhit_vec.push_back(nhit);
//   }

//   if(!nhit_vec.empty()) {
//     Int_t max_nhit = nhit_vec.front();
//     Int_t min_nhit = nhit_vec.back();
//     for(Int_t nhit=max_nhit; nhit>=min_nhit; nhit--) {
//       auto itr1 = std::find(nhit_vec.begin(), nhit_vec.end(), nhit);
//       if(itr1 == nhit_vec.end())
//         continue;

//       size_t index1 = std::distance(nhit_vec.begin(), itr1);

//       auto itr2 = std::find(nhit_vec.rbegin(), nhit_vec.rend(), nhit);
//       size_t index2 = nhit_vec.size() - std::distance(nhit_vec.rbegin(), itr2) - 1;

//       index_pair_vec.push_back(index_pair(index1, index2));
//     }
//   }

//   for(Int_t i=0; i<index_pair_vec.size(); i++) {
//     std::stable_sort(trackCont.begin() + index_pair_vec[i].first,
//                      trackCont.begin() +  index_pair_vec[i].second + 1, DCLTrackComp_Chisqr());
//   }

// #if 0
//   DebugPrint(trackCont, arg+" After Sorting (chisqr)");
// #endif

//   if(delete_flag) {
//     for(Int_t i = index_pair_vec.size()-1; i>=0; --i) {
//       DeleteDuplicatedTracks(trackCont, index_pair_vec[i].first, index_pair_vec[i].second, 0.);
//     }
//   }

// #if 0
//   DebugPrint(trackCont, arg+" After Deleting in each hit number");
// #endif


//   std::stable_sort(trackCont.begin(), trackCont.end(), comp);

// #if 0
//   DebugPrint(trackCont, arg+" After Sorting with comp func ");
// #endif

//   if(delete_flag) DeleteDuplicatedTracks(trackCont);

// #if 0
//   DebugPrint(trackCont, arg+" After Deleting ");
// #endif

//   CalcTracks(trackCont);
//   del::ClearContainerAll(candCont);
// }

}

//_____________________________________________________________________________
namespace tpc
{
//_____________________________________________________________________________
// Int_t
// LocalTrackSearch(const std::vector<TPCHitContainer>& HitCont,
//                  std::vector<TPCLocalTrack*>& TrackCont,
//                  Int_t MinNumOfHits /*=8*/)
// {
//   static const Double_t MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
//   static const Double_t MaxLayerCut = gUser.GetParameter("TPCMaxLayerCut");
//   Bool_t status = true;

//   //    if(valueHall) { // TODO
//   //    }

//   // y = p0 + p1 * x
//   Double_t p0[MaxNumOfTrackTPC];
//   Double_t p1[MaxNumOfTrackTPC];
//   static const Int_t    Li_theta_ndiv = 200;
//   static const Double_t Li_theta_min  =   0;
//   static const Double_t Li_theta_max  = 180;
//   static const Int_t    Li_r_ndiv =  200;
//   static const Double_t Li_r_min  = -600;
//   static const Double_t Li_r_max  =  600;
//   static const Int_t MinClusterSize = gUser.GetParameter("MinClusterSizeTPC");

//   //for TPC linear track
//   // r = x * cos(theta) + y * sin(theta)
//   TH2D h1("hist_linear",";theta (deg.); r (mm)",
//           Li_theta_ndiv, Li_theta_min, Li_theta_max,
//           Li_r_ndiv, Li_r_min, Li_r_max);

//   std::vector<std::vector<Int_t>> flag;
//   flag.resize(NumOfLayersTPC);
//   for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
//     flag[layer].resize(HitCont[layer].size(), 0);
//   }

//   std::vector<Double_t> hough_x;
//   std::vector<Double_t> hough_y;

//   for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
//     h1.Reset();
//     for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
//       for(Int_t ci=0, n=HitCont[layer].size(); ci<n; ci++){
//         if(flag[layer][ci]>0) continue;
//         TPCHit* hit = HitCont[layer][ci];
//         TVector3 pos = hit->GetPosition();
//         for(Int_t ti=0; ti<Li_theta_ndiv; ti++){
//           Double_t theta = Li_theta_min+ti*(Li_theta_max-Li_theta_min)/Li_theta_ndiv;
//           h1.Fill(theta, cos(theta*acos(-1)/180.)*pos.Z()
//                   +sin(theta*acos(-1)/180.)*pos.X());
//           if(fabs(cos(theta*acos(-1)/180.)*pos.Z()+sin(theta*acos(-1)/180.)*pos.X())>Li_r_max)
//             std::cout<<"Hough: out of range:"<<cos(theta*acos(-1)/180.)*pos.Z()+sin(theta*acos(-1)/180.)*pos.X()<<std::endl;
//         }
//       } // cluster
//     } // layer
//       //      if(h1.GetMaximum() < MinNumOfHits){
//     if(h1.GetMaximum() < MinNumOfHits/2){
//       //h1.Delete();
//       h1.Reset();
//       break;
//     }




//     Int_t maxbin = h1.GetMaximumBin();
//     Int_t mx,my,mz;
//     h1.GetBinXYZ(maxbin, mx, my, mz);
//     Double_t mtheta = h1.GetXaxis()->GetBinCenter(mx)*acos(-1)/180.;
//     Double_t mr = h1.GetYaxis()->GetBinCenter(my);

//     Bool_t hough_flag = true;
//     for(Int_t i=0; i<hough_x.size(); ++i){
//       Int_t bindiff = fabs(mx-hough_x[i])+fabs(my-hough_y[i]);
//       if(bindiff<=4)
//         hough_flag = false;
//     }
//     hough_x.push_back(mx);
//     hough_y.push_back(my);
//     if(!hough_flag)
//       continue;

//     TPCLocalTrack *track = new TPCLocalTrack();
//     p0[tracki] = mr/sin(mtheta);
//     p1[tracki] = -cos(mtheta)/sin(mtheta);

//     for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
//       for(Int_t ci=0, n=HitCont[layer].size(); ci<n; ci++){
//         //if(flag[layer][ci]>0) continue;
//         TPCHit* hit = HitCont[layer][ci];
//         TVector3 pos = hit->GetPosition();
//         Double_t dist = fabs(p1[tracki]*pos.Z()-pos.X()+p0[tracki])/sqrt(pow(p1[tracki],2)+1);

//         //if(dist < MaxHoughWindow && hit->GetClusterSize()>=MinClusterSize){
// 	if(dist < MaxHoughWindow && hit->GetClusterSize()>=MinClusterSize&& layer < MaxLayerCut){
//           track->AddTPCHit(new TPCLTrackHit(hit));
//           //     track_he->AddTPCHit(new TPCLTrackHit(hit));
//           flag[layer][ci]++;
//         }
//       }
//     }
//     //temporary (x0, y0) are position at Target position
//     track->SetAx(p0[tracki]+p1[tracki]*tpc::ZTarget);
//     track->SetAu(p1[tracki]);

//     if(track->DoFit(MinNumOfHits)){
//       TrackCont.push_back(track);
//     }
//     else
//       delete track;

//     h1.Reset();
//   }//track

//   CalcTracks(TrackCont);
//   // std::cout<<"event end"<<std::endl;

//   return status? TrackCont.size() : -1;

//   return 0;
// }

//_____________________________________________________________________________
Int_t
LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
                 std::vector<TPCLocalTrack*>& TrackCont,
                 Int_t MinNumOfHits /*=8*/)
{
  static const Double_t MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  static const Double_t MaxLayerCut = gUser.GetParameter("TPCMaxLayerCut");
  Bool_t status = true;

  // y = p0 + p1 * x
  Double_t p0[MaxNumOfTrackTPC];
  Double_t p1[MaxNumOfTrackTPC];
  static const Int_t    Li_theta_ndiv = 180;
  static const Double_t Li_theta_min  =   0;
  static const Double_t Li_theta_max  = 180;
  static const Int_t    Li_rho_ndiv =  180;
  static const Double_t Li_rho_min  = -720;
  static const Double_t Li_rho_max  =  720;

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  //for TPC linear track
  // rho = x * cos(theta) + y * sin(theta)
  TH2D h1("hist_linear",
          Form("%s XZ linear hough; #theta (deg.); #rho (mm)", FUNC_NAME.Data()),
          Li_theta_ndiv, Li_theta_min, Li_theta_max,
          Li_rho_ndiv, Li_rho_min, Li_rho_max);

  std::vector<std::vector<Int_t>> flag;
  flag.resize(NumOfLayersTPC);
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    flag[layer].resize(ClCont[layer].size(), 0);
  }

  std::vector<Double_t> hough_x;
  std::vector<Double_t> hough_y;
  std::vector<Double_t> hough_z;

  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    h1.Reset();
    h1.SetTitle(Form("%dth track", tracki));
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
        if(flag[layer][ci]>0) continue;
        const auto& pos = ClCont[layer][ci]->GetPosition();
        for(Int_t ti=0; ti<Li_theta_ndiv; ti++){
          Double_t theta = Li_theta_min +
            (ti+0.5)*(Li_theta_max-Li_theta_min)/Li_theta_ndiv;
          Double_t rho = TMath::Cos(theta*TMath::DegToRad())*pos.Z()
            + TMath::Sin(theta*TMath::DegToRad())*pos.X();
          h1.Fill(theta, rho);
          if(TMath::Abs(rho) > Li_rho_max){
            hddaq::cerr << FUNC_NAME << " Out of Hough range: "
                        << rho <<std::endl;
          }
        } // theta
      } // cluster
    } // layer

#if DebugEvDisp
    h1.Draw("colz");
    gPad->Modified();
    gPad->Update();
    c1.Print("c1.pdf");
    getchar();
#endif

    //      if(h1.GetMaximum() < MinNumOfHits){
    if(h1.GetMaximum() < MinNumOfHits/2){
      h1.Reset();
      break;
    }

    Int_t maxbin = h1.GetMaximumBin();
    Int_t mx,my,mz;
    h1.GetBinXYZ(maxbin, mx, my, mz);
    Double_t mtheta = h1.GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
    Double_t mr = h1.GetYaxis()->GetBinCenter(my);
    Bool_t hough_flag = true;
    for(Int_t i=0; i<hough_x.size(); ++i){
      Int_t bindiff = TMath::Abs(mx-hough_x[i]) + TMath::Abs(my-hough_y[i]) + TMath::Abs(mz-hough_z[i]);
      if(bindiff<=4)
        hough_flag = false;
    }
    hough_x.push_back(mx);
    hough_y.push_back(my);
    hough_z.push_back(mz);
    if(!hough_flag) continue;
    TPCLocalTrack *track = new TPCLocalTrack;
    p0[tracki] = mr/TMath::Sin(mtheta);
    p1[tracki] = -TMath::Cos(mtheta)/TMath::Sin(mtheta);
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
        const auto cl = ClCont[layer][ci];
        //if(flag[layer][ci]>0) continue;
        TPCHit* hit = cl->GetMeanHit();
        const TVector3& pos = cl->GetPosition();
        Double_t dist = TMath::Abs(p1[tracki]*pos.Z()-pos.X() +
                                   p0[tracki])/TMath::Sqrt(TMath::Sq(p1[tracki])+1);
        if(dist > MaxHoughWindow) continue;
	if(layer < MaxLayerCut){
	  hit->SetHoughDist(dist);
          track->AddTPCHit(new TPCLTrackHit(hit));
          flag[layer][ci]++;
        }
      }
    }
    //temporary (x0, y0) are position at Target position
    track->SetAx(p0[tracki]+p1[tracki]*tpc::ZTarget);
    track->SetAu(p1[tracki]);

    if(track->GetNHit() >= MinNumOfHits && track->DoFit(MinNumOfHits)){
      TrackCont.push_back(track);
    }else{
      delete track;
    }

    h1.Reset();
  }//track

  CalcTracks(TrackCont);

  return status? TrackCont.size() : -1;

  return 0;
}

//_____________________________________________________________________________
// Int_t
// LocalTrackSearchHelix(const std::vector<TPCHitContainer>& HitCont,
//                       std::vector<TPCLocalTrackHelix*>& TrackCont,
//                       Int_t MinNumOfHits /*=8*/)
// {
//   static const Double_t MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
//   static const Double_t MaxLayerCut = gUser.GetParameter("TPCMaxLayerCut");
//   // static const Double_t DECut_TPCTrack = gUser.GetParameter("DECut_TPCTrack");
//   Bool_t status = true;

//   //    if(valueHall) { // TODO
//   //    }
//   const Double_t Const = 0.299792458;

//   //(x - (r + rd)*cos(theta))^2 + (y - (r + rd)*sin(theta))^2 = r^2
//   // p = r * Const; // 1T
//   // Parameters
//   Double_t hough_p[MaxNumOfTrackTPC];
//   Double_t hough_theta[MaxNumOfTrackTPC];
//   Double_t hough_rd[MaxNumOfTrackTPC];

//   const Int_t nBin_rdiff = 11;
//   const Double_t rdiff_min = -110.;
//   const Double_t rdiff_max = 110.;

//   const Int_t nBin_theta = 180;
//   //const Int_t nBin_theta = 90;
//   //const Int_t nBin_theta = 360;
//   const Double_t theta_min = -1.*acos(-1);
//   const Double_t theta_max = acos(-1);

//   //  const Int_t nBin_p = 300;
//   //const Int_t nBin_p = 600;
//   const Int_t nBin_p = 300;
//   const Double_t pmin = 50.;//MeV/c
//   //const Double_t pmax = 2050.;//MeV/c
//   const Double_t pmax = 1550.;//MeV/c

//   //for TPC circle track

//   //start from hougy by using the HoughYcut info
//   int Max_tracki_houghY =0;
//   for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
//     for(Int_t ci=0, n=HitCont[layer].size(); ci<n; ci++){
//       TPCHit* hit = HitCont[layer][ci];
//       int ihoughy_size = hit->GetHoughY_num_size();
//       //if(ihoughy_size>1)
//       for(int ih=0; ih<ihoughy_size; ++ih){
// 	int ihoughy =  hit->GetHoughY_num(ih);
// 	if(ihoughy>Max_tracki_houghY)
// 	  Max_tracki_houghY = ihoughy;
//       }
//     }
//   }
//   ++Max_tracki_houghY;
//   //  std::cout<<"Max_tracki="<<Max_tracki_houghY<<std::endl;
//   TH3D *Ci_hist=new TH3D("hist_circle",";rd (mm); theta (rad); p(MeV/c)",
//                          nBin_rdiff, rdiff_min,  rdiff_max,
//                          nBin_theta, theta_min, theta_max,
//                          nBin_p,pmin,pmax);


//   std::vector<Double_t> hough_x;
//   std::vector<Double_t> hough_y;

//   for(Int_t ity=0; ity<Max_tracki_houghY+1; ity++){
//     //for(Int_t ity=0; ity<Max_tracki_houghY; ity++){
//     for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
//       //    for(Int_t ity=0; ity<Max_tracki_houghY+1; ity++){
//       Ci_hist->Reset();
//       for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
// 	for(Int_t ci=0, n=HitCont[layer].size(); ci<n; ci++){
// 	  //        if(flag[layer][ci]>0) continue;
// 	  TPCHit* hit = HitCont[layer][ci];
// 	  if(ity<Max_tracki_houghY){
// 	    bool status_houghy = false;
// 	    int ihoughy_size = hit->GetHoughY_num_size();
// 	    for(int ih=0; ih<ihoughy_size; ++ih){
// 	      int ihoughy = hit->GetHoughY_num(ih);
// 	      //if(ity==ihoughy-1)
// 	      if(ity==ihoughy)
// 		status_houghy = true;
// 	    }
// 	    if(!status_houghy)
// 	      continue;
// 	  }
// 	  if(hit->GetHoughFlag()>0) continue;

// 	  TVector3 pos = hit->GetPosition();
// 	  for(Int_t ird=0; ird<nBin_rdiff; ++ird){
// 	    Double_t rd = Ci_hist->GetXaxis()->GetBinCenter(ird+1);
// 	    for(Int_t ip=0; ip<nBin_p; ++ip){
// 	      Double_t x = -pos.x();
// 	      Double_t y = pos.z()-tpc::ZTarget;
// 	      Double_t p = Ci_hist->GetZaxis()->GetBinCenter(ip+1);
// 	      Double_t r = p/(Const*1.);//1T

// 	      //a*sin(theta) + b*cos(theta) +c = 0
// 	      Double_t a = 2.*(r+rd)*y;
// 	      Double_t b = 2.*(r+rd)*x;
// 	      Double_t c = -1.*(rd*rd + 2.*r*rd + x*x + y*y);

// 	      Double_t r0 = sqrt(a*a + b*b);
// 	      if(fabs(-1.*c/r0)>1.){
// 		// std::cout<<"No solution, "
// 		//   <<"x:"<<x<<", y:"<<y
// 		//   <<", r:"<<r<<", rd:"<<rd<<std::endl;
// 		continue;
// 	      }
// 	      Double_t theta1_alpha =  asin(-1.*c/r0);
// 	      Double_t theta2_alpha;
// 	      if(theta1_alpha>0.)
// 		theta2_alpha = acos(-1.) - theta1_alpha;
// 	      else
// 		theta2_alpha = -1.*acos(-1.) - theta1_alpha;

// 	      Double_t theta_alpha = atan2(b, a);

// 	      Double_t xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
// 	      Double_t ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
// 	      Double_t r_re1 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

// 	      Double_t xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
// 	      Double_t ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
// 	      Double_t r_re2 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

// 	      Double_t theta1 = atan2(ycenter1, xcenter1);
// 	      Double_t theta2 = atan2(ycenter2, xcenter2);

// 	      if(TMath::IsNaN(theta1)){
// 		std::cout<<"theta1="<<theta1<<", x="<<x<<", y="<<y
// 			 <<"rd="<<rd<<", r"<<r<<std::endl;
// 	      }


// 	      if(fabs(r-r_re1)>0.01||fabs(r-r_re2)>0.01){
// 		std::cout<<"r="<<r<<", r_re1="<<r_re1<<", r_re1="<<r_re2<<std::endl;
// 		std::cout<<"x:"<<x<<", y:"<<y
// 			 <<", theta1:"<<theta1<<", theta2:"<<theta2
// 			 <<", theta_alpha:"<<theta_alpha<<std::endl;
// 	      }
// 	      Ci_hist->Fill(rd, theta1, p);
// 	      Ci_hist->Fill(rd, theta2, p);
// 	      // std::cout<<"rd: "<<rd<<", "
// 	      //         <<"theta1: "<<theta1<<", "
// 	      //         <<"theta2: "<<theta2<<", "
// 	      //         <<"p: "<<p<<std::endl;
// 	    }
// 	  }
// 	}// cluster
//       } // layer

//       if(Ci_hist->GetMaximum() < MinNumOfHits/2){
// 	Ci_hist->Reset();
// 	break;
//       }
//       //std::cout<<"Maxbin0: "<<Ci_hist.GetMaximum()<<std::endl;


//       Int_t maxbin = Ci_hist->GetMaximumBin();
//       Int_t mx,my,mz;
//       Ci_hist->GetBinXYZ(maxbin, mx, my, mz);

//       Bool_t hough_flag = true;
//       for(Int_t i=0; i<hough_x.size(); ++i){
// 	Int_t bindiff = fabs(mx-hough_x[i])+fabs(my-hough_y[i]);
// 	if(bindiff<=4)
// 	  //if(bindiff<=8)
// 	  hough_flag = false;
//       }
//       hough_x.push_back(mx);
//       hough_y.push_back(my);
//       if(!hough_flag)
// 	continue;

//       TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
//       hough_rd[tracki] = Ci_hist->GetXaxis()->GetBinCenter(mx);
//       hough_theta[tracki] = Ci_hist->GetYaxis()->GetBinCenter(my);
//       hough_p[tracki] = Ci_hist->GetZaxis()->GetBinCenter(mz);
//       Double_t hough_r = hough_p[tracki]/Const;
//       Double_t hough_cx = (hough_r + hough_rd[tracki])*cos(hough_theta[tracki]);
//       Double_t hough_cy = (hough_r + hough_rd[tracki])*sin(hough_theta[tracki]);


//       //std::cout<<"Hough (x,z) Maxbin: "<<Ci_hist.GetMaximum()<<std::endl;
//       //     std::cout<<""<<std::endl;
//       for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
// 	for(Int_t ci=0, n=HitCont[layer].size(); ci<n; ci++){
// 	  //if(flag[layer][ci]>0) continue;
// 	  TPCHit* hit = HitCont[layer][ci];
// 	  if(hit->GetHoughFlag()>0) continue;
// 	  TVector3 pos = hit->GetPosition();
// 	  Double_t x = -pos.x();
// 	  Double_t y = pos.z()-tpc::ZTarget;
// 	  // Double_t de = hit->GetCharge();
// 	  // Double_t xcenter1 = (hough_r + hough_rd[tracki])*cos(hough_theta[tracki]);
// 	  // Double_t ycenter1 = (hough_r + hough_rd[tracki])*sin(hough_theta[tracki]);
// 	  Double_t r_cal = sqrt(pow(x-hough_cx,2) + pow(y-hough_cy,2));

// 	  Double_t dist = fabs(r_cal - hough_r);
// 	  //if(dist < MaxHoughWindow && layer < MaxLayerCut && de>DECut_TPCTrack){
// 	  if(dist < MaxHoughWindow && layer < MaxLayerCut){
// 	    if(ity<Max_tracki_houghY){
// 	      bool status_houghy = false;
// 	      int ihoughy_size = hit->GetHoughY_num_size();
// 	      for(int ih=0; ih<ihoughy_size; ++ih){
// 		int ihoughy = hit->GetHoughY_num(ih);
// 		//if(ity==ihoughy-1)
// 		if(ity==ihoughy)
// 		  status_houghy = true;
// 	      }
// 	      if(status_houghy)
// 		track->AddTPCHit(new TPCLTrackHit(hit));
// 	    }
// 	    else
// 	      track->AddTPCHit(new TPCLTrackHit(hit));
// 	    //flag[layer][ci]++;
// 	  }
// 	}
//       }
//       // track->SetAdrho(hough_rd[tracki]);
//       // track->SetAphi0(hough_theta[tracki]);
//       // track->SetArho(1./hough_r);
//       track->SetAcx(hough_cx);
//       track->SetAcy(hough_cy);
//       track->SetAr(hough_r);

//       if(ity==0){
// 	if(track->DoFit(3,1)){
// 	  TrackCont.push_back(track);
// 	  Int_t nh = track->GetNHit();
// 	  for( int ih=0; ih<nh; ++ih ){
// 	    TPCHit *hit = track->GetHit( ih )->GetHit();
// 	    if( !hit ) continue;
// 	    hit->SetHoughFlag(1+hit->GetHoughFlag());
// 	  }
// 	}
// 	else
// 	  delete track;
//       }
//       else{
// 	if(track->DoFit(MinNumOfHits,0)){
// 	  TrackCont.push_back(track);
// 	  Int_t nh = track->GetNHit();
// 	  for( int ih=0; ih<nh; ++ih ){
// 	    TPCHit *hit = track->GetHit( ih )->GetHit();
// 	    if( !hit ) continue;
// 	    hit->SetHoughFlag(1+hit->GetHoughFlag());
// 	  }
// 	}
// 	else
// 	  delete track;
//       }
//       //      Ci_hist->Reset();

//       //delete track;
//     }//hough y cut
//   }//track
//   CalcTracks(TrackCont);
//   delete Ci_hist;
//   return status? TrackCont.size() : -1;

//   return 0;
// }

//_____________________________________________________________________________
Int_t
LocalTrackSearchHelix(const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      Int_t MinNumOfHits /*=8*/)
{
  static const Double_t MaxHoughWindow = gUser.GetParameter("MaxHoughWindow");
  static const Double_t MaxHoughWindowY = gUser.GetParameter("MaxHoughWindowY");
  static const Double_t MaxLayerCut = gUser.GetParameter("TPCMaxLayerCut");

  Bool_t status = true;

  const Double_t Const = 0.299792458;
  const Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  //(x - (r + rd)*cos(theta))^2 + (y - (r + rd)*sin(theta))^2 = r^2
  // p = r * Const * dMagneticField;
  // Parameters
  Double_t hough_p[MaxNumOfTrackTPC];
  Double_t hough_theta[MaxNumOfTrackTPC];
  Double_t hough_rd[MaxNumOfTrackTPC];

  const Int_t nBin_rdiff = 11;
  const Double_t rdiff_min = -110.;
  const Double_t rdiff_max = 110.;

  //const Int_t nBin_theta = 180; //previous value
  const Int_t nBin_theta = 720;
  const Double_t theta_min = -1.*acos(-1);
  const Double_t theta_max = acos(-1);

  //const Int_t nBin_p = 300; //previous value
  const Int_t nBin_p = 900;
  const Double_t pmin = 50.;//MeV/c
  //const Double_t pmax = 2050.;//MeV/c
  const Double_t pmax = 1550.;//MeV/c

  //for hough Y
  //const int    thetaY_ndiv = 360; //previous value
  const int    thetaY_ndiv = 1080;
  const double thetaY_min  =   0.;
  const double thetaY_max  = 180.;
  const int    r_ndiv =  2000;
  const double r_min  = -5000.;
  const double r_max  =  5000.;

  //for TPC circle track

  //start from hougy by using the HoughYcut info
  // int Max_tracki_houghY =0;
  // for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
  //   for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
  //     auto cl = ClCont[layer][ci];
  //     TPCHit* hit = cl->GetMeanHit();
  //     int ihoughy_size = hit->GetHoughY_num_size();
  //     //if(ihoughy_size>1)
  //     for(int ih=0; ih<ihoughy_size; ++ih){
  // 	int ihoughy =  hit->GetHoughY_num(ih);
  // 	if(ihoughy>Max_tracki_houghY)
  // 	  Max_tracki_houghY = ihoughy;
  //     }
  //   }
  // }
  // ++Max_tracki_houghY;
  //  std::cout<<"Max_tracki="<<Max_tracki_houghY<<std::endl;
	bool RemoveBeam = true;
	TPCBeamRemover BeamRemover(ClCont);
	BeamRemover.Enable(RemoveBeam);
	BeamRemover.SearchAccidentalBeam(-50,50);
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(auto cl:ClCont[layer]){
			auto hit = cl->GetMeanHit();
			if(BeamRemover.IsAccidental(hit)) hit->SetHoughFlag(2);//2 -> Accidental Flag;
		}
	}
#if 0
  //Accidental beam kill
  //TrackBeam->SetHoughFlag(4);
#endif
  //Hough-transform
  TH3D *Ci_hist=new TH3D("hist_circle",";rd (mm); theta (rad); p(MeV/c)",
			 nBin_rdiff, rdiff_min,  rdiff_max,
			 nBin_theta, theta_min, theta_max,
			 nBin_p, pmin, pmax);
  TH2D *histY=new TH2D("histY",";theta (deg.); r (mm)",
		       thetaY_ndiv, thetaY_min, thetaY_max, r_ndiv, r_min, r_max);

  std::vector<Double_t> hough_x;
  std::vector<Double_t> hough_y;
  std::vector<Double_t> hough_z;

  //std::cout<<"before tracks : "<<TrackCont.size()<<" "<<TrackContFailed.size()<<std::endl;

  bool prev_add = true;
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    std::chrono::seconds sec;
    auto before_hough = std::chrono::high_resolution_clock::now();
    if(!prev_add) continue;
    prev_add = false;
    Ci_hist->Reset();
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	// if(ity<Max_tracki_houghY){
	//   bool status_houghy = false;
	//   int ihoughy_size = hit->GetHoughY_num_size();
	//   for(int ih=0; ih<ihoughy_size; ++ih){
	//     int ihoughy = hit->GetHoughY_num(ih);
	//     //if(ity==ihoughy-1)
	//     if(ity==ihoughy)
	// 	status_houghy = true;
	//   }
	//   if(!status_houghy)
	//     continue;
	// }
	if(hit->GetHoughFlag()>0)
	  continue;

	TVector3 pos = cl->GetPosition();
	for(Int_t ird=0; ird<nBin_rdiff; ++ird){
	  Double_t rd = Ci_hist->GetXaxis()->GetBinCenter(ird+1);
	  for(Int_t ip=0; ip<nBin_p; ++ip){
	    Double_t x = -pos.x();
	    Double_t y = pos.z()-tpc::ZTarget;
	    Double_t p = Ci_hist->GetZaxis()->GetBinCenter(ip+1);
	    Double_t r = p/(Const*dMagneticField);

	    //a*sin(theta) + b*cos(theta) + c = 0
	    Double_t a = 2.*(r+rd)*y;
	    Double_t b = 2.*(r+rd)*x;
	    Double_t c = -1.*(rd*rd + 2.*r*rd + x*x + y*y);

	    Double_t r0 = sqrt(a*a + b*b);
	    if(fabs(-1.*c/r0)>1.){
	      continue;
	    }
	    Double_t theta1_alpha = asin(-1.*c/r0);
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
	  } // p bins
	} // rdiff bins
      } // cluster
    } // layer

#if DebugDisp
    std::cout<<"xz hough vote : "<<Ci_hist->GetMaximum()<<std::endl;
#endif

    if(Ci_hist->GetMaximum() < MinNumOfHits/2){
      Ci_hist->Reset();
      break;
    }
    Int_t maxbin = Ci_hist->GetMaximumBin();
    Int_t mx,my,mz;
    Ci_hist->GetBinXYZ(maxbin, mx, my, mz);

    Bool_t hough_flag = true;
    for(Int_t i=0; i<hough_x.size(); ++i){
      Int_t bindiff = fabs(mx-hough_x[i])+fabs(my-hough_y[i])+fabs(mz-hough_z[i]);
      if(bindiff<=2){
	hough_flag = false;
#if DebugDisp
	std::cout<<"i:"<<i
		 <<", hough_x:"<<hough_x[i]
		 <<", mx:"<<mx
		 <<", hough_y:"<<hough_y[i]
		 <<", my:"<<my
		 <<", hough_z:"<<hough_z[i]
		 <<", mz:"<<mz<<std::endl;
#endif
      }
    }
    hough_x.push_back(mx);
    hough_y.push_back(my);
    hough_z.push_back(mz);
    if(!hough_flag)
      continue;

    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    hough_rd[tracki] = Ci_hist->GetXaxis()->GetBinCenter(mx);
    hough_theta[tracki] = Ci_hist->GetYaxis()->GetBinCenter(my);
    hough_p[tracki] = Ci_hist->GetZaxis()->GetBinCenter(mz);
#if DebugDisp
    std::cout<<FUNC_NAME+" Circle hough maxbin:"<<Ci_hist->GetBinContent(maxbin)<<std::endl;
    std::cout<<FUNC_NAME+" Hough(rd, theta, p)=("<<hough_rd[tracki]<<", "
	     <<hough_theta[tracki]<<", "<<hough_p[tracki]<<")"<<std::endl;
#endif

    Double_t hough_r = hough_p[tracki]/(Const*dMagneticField);
    Double_t hough_cx = (hough_r + hough_rd[tracki])*cos(hough_theta[tracki]);
    Double_t hough_cy = (hough_r + hough_rd[tracki])*sin(hough_theta[tracki]);

    histY->Reset();
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	Double_t x = -pos.x();
	Double_t y = pos.z()-tpc::ZTarget;
	Double_t r_cal = sqrt(pow(x-hough_cx,2) + pow(y-hough_cy,2));
	Double_t dist = fabs(r_cal - hough_r);
	if(dist < MaxHoughWindow && layer < MaxLayerCut){
	  for(int ti=0; ti<thetaY_ndiv; ti++){
	    double theta = thetaY_min+ti*(thetaY_max-thetaY_min)/thetaY_ndiv;
	    double tmpx = -pos.x();
	    double tmpy = pos.z()-tpc::ZTarget;
	    double tmpz = pos.y();
	    double tmp_t = atan2(tmpy - hough_cy,
				 tmpx - hough_cx);
	    double tmp_xval = hough_r * tmp_t;
	    histY->Fill(theta, cos(theta*acos(-1.)/180.)*tmp_xval
			+sin(theta*acos(-1.)/180.)*tmpz);
	  }
	}
      }
    }

    int maxbinY = histY->GetMaximumBin();
    int mxY,myY,mzY;
    histY->GetBinXYZ(maxbinY, mxY, myY, mzY);

    double mtheta = histY->GetXaxis()->GetBinCenter(mxY)*acos(-1.)/180.;
    double mr = histY->GetYaxis()->GetBinCenter(myY);
    double p0 = mr/sin(mtheta);
    double p1 = -cos(mtheta)/sin(mtheta);

    //Hough space
    //mr = tmp_xval*cos(mtheta) + tmpz*sin(mtheta)
    //tmpz = p0 + p1*tmp_xval
#if DebugDisp
    std::cout<<FUNC_NAME+" (Y,phi) hough maxbin:"<<histY->GetBinContent(maxbinY)<<std::endl;
    std::cout<<FUNC_NAME+" Hough (theta, r)=("<<mtheta<<", "<<mr<<")"<<std::endl;
#endif

    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()>0) continue;
	TVector3 pos = cl->GetPosition();
	Double_t x = -pos.x();
	Double_t y = pos.z() - tpc::ZTarget;
	Double_t r_cal = sqrt(pow(x - hough_cx, 2) + pow(y - hough_cy, 2));
	Double_t dist = fabs(r_cal - hough_r);
	if(dist < MaxHoughWindow && layer < MaxLayerCut){

	  double tmpx = -pos.x();
	  double tmpy = pos.z() - tpc::ZTarget;
	  double tmpz = pos.y();
	  double tmp_t = atan2(tmpy - hough_cy,
			       tmpx - hough_cx);
	  double tmp_xval = hough_r * tmp_t;
	  double distY = fabs(p1*tmp_xval - tmpz + p0)/sqrt(pow(p1, 2) + 1);
	  if(distY < MaxHoughWindowY){
	    hit->SetHoughDist(dist);
	    hit->SetHoughDistY(distY);
	    track->AddTPCHit(new TPCLTrackHit(hit));
	    prev_add = true;
#if DebugDisp
	    track->Print(FUNC_NAME+"AddTPCHit");
#endif
	  }
	}
      } //ci
    } //layer

    track->SetAcx(hough_cx);
    track->SetAcy(hough_cy);
    track->SetAr(hough_r);
    track->SetAz0(p0);
    track->SetAdz(p1);

    auto after_hough = std::chrono::high_resolution_clock::now();
    sec = std::chrono::duration_cast<std::chrono::seconds>(before_hough - after_hough);
    track->SetHoughTime((int) sec.count());
    if(track->DoFit(MinNumOfHits,0)){ //first fit
      track->SetHoughFlag(1);
#if DebugDisp
      track->Print(FUNC_NAME+"SetHoughFlag");
#endif
      //Residual check with other hits
      bool Add_rescheck = false;
      for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
	for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	  auto cl = ClCont[layer][ci];
	  TPCHit* hit = cl->GetMeanHit();
	  if(hit->GetHoughFlag()==1) continue;
	  if(hit->GetHoughFlag()==2) continue;
	  TVector3 pos = cl->GetPosition();
	  TVector3 res =  TVector3(cl->ResolutionX(),
				   cl->ResolutionY(),
				   cl->ResolutionZ());
	  double resi=0.;
	  if(track->Residual_check(pos,res,resi)){
	    track->AddTPCHit(new TPCLTrackHit(hit));
	    Add_rescheck = true;
	  }
	}
      } //for (layer)
      if(!Add_rescheck){
	auto nohitadd = std::chrono::high_resolution_clock::now();
	sec = std::chrono::duration_cast<std::chrono::seconds>(nohitadd - after_hough);
	track->SetFitTime((int) sec.count());
	track->SetFlag(1);
	TrackCont.push_back(track);
      }
      else{
	if(track->DoFit(MinNumOfHits,0)){
	  auto hitadd = std::chrono::high_resolution_clock::now();
	  sec = std::chrono::duration_cast<std::chrono::seconds>(hitadd - after_hough);
	  track->SetFitTime((int) sec.count());
	  track->SetFlag(2);
	  TrackCont.push_back(track);
#if DebugDisp
	  track->Print(FUNC_NAME+" Add (residual check) SetHoughFlag");
#endif
	}
	else{
	  std::cout<<"edge case"<<std::endl;
	  auto hitadd_fail = std::chrono::high_resolution_clock::now();
	  sec = std::chrono::duration_cast<std::chrono::seconds>(hitadd_fail - after_hough);
	  track->SetFitTime((int) sec.count());
	  track->SetHelixUsingHoughParam();
	  track->SetFlag(3);
	  if(track->GetNHit() > 0.5*MinNumOfHits) TrackContFailed.push_back(track);
	  else delete track;
#if DebugDisp
	  track->Print(FUNC_NAME+" Error !! (TPCTrackSearch)");
	  getchar();
#endif
	}
      }
    }
    else{ //first fitting failed
      std::cout<<std::endl;
      auto after_1stfit_fail = std::chrono::high_resolution_clock::now();
      sec = std::chrono::duration_cast<std::chrono::seconds>(after_1stfit_fail - after_hough);
      track->SetFitTime((int) sec.count());
      track->SetHoughFlag(100);
      track->SetHelixUsingHoughParam();
      track->SetFlag(4);
      if(track->GetNHit() > 0.5*MinNumOfHits) TrackContFailed.push_back(track);
      else delete track;
#if DebugDisp
      track->Print(FUNC_NAME+" Delete tarcks");
#endif
    }
  }//tracki

#if DummyHitTest
  for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
    TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
	auto cl = ClCont[layer][ci];
	TPCHit* hit = cl->GetMeanHit();
	if(hit->GetHoughFlag()==0||hit->GetHoughFlag()==100)
	  track->AddTPCHit(new TPCLTrackHit(hit));
      }
    }
    if(track->DoFit(MinNumOfHits,0)){
      track->SetHoughFlag(1);
      TrackCont.push_back(track);
#if DebugDisp
      track->Print(FUNC_NAME+" SetHoughFlag");
#endif
    }
    else{
      delete track;
      break;
    }
  }
#endif
  CalcTracks(TrackCont);
  CalcTracks(TrackContFailed);
  delete Ci_hist;
  delete histY;
  return status? TrackCont.size() : -1;

  return 0;

  // for(Int_t ity=0; ity<Max_tracki_houghY+1; ity++){
  //   //for(Int_t ity=0; ity<Max_tracki_houghY; ity++){
  //   for(Int_t tracki=0; tracki<MaxNumOfTrackTPC; tracki++){
  //     //    for(Int_t ity=0; ity<Max_tracki_houghY+1; ity++){
  //     Ci_hist->Reset();
  //     for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
  // 	for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
  // 	  auto cl = ClCont[layer][ci];
  // 	  //        if(flag[layer][ci]>0) continue;
  // 	  TPCHit* hit = cl->GetMeanHit();
  // 	  if(ity<Max_tracki_houghY){
  // 	    bool status_houghy = false;
  // 	    int ihoughy_size = hit->GetHoughY_num_size();
  // 	    for(int ih=0; ih<ihoughy_size; ++ih){
  // 	      int ihoughy = hit->GetHoughY_num(ih);
  // 	      //if(ity==ihoughy-1)
  // 	      if(ity==ihoughy)
  // 		status_houghy = true;
  // 	    }
  // 	    if(!status_houghy)
  // 	      continue;
  // 	  }
  // 	  if(hit->GetHoughFlag()>0) continue;

  // 	  TVector3 pos = cl->GetPosition();
  // 	  for(Int_t ird=0; ird<nBin_rdiff; ++ird){
  // 	    Double_t rd = Ci_hist->GetXaxis()->GetBinCenter(ird+1);
  // 	    for(Int_t ip=0; ip<nBin_p; ++ip){
  // 	      Double_t x = -pos.x();
  // 	      Double_t y = pos.z()-tpc::ZTarget;
  // 	      Double_t p = Ci_hist->GetZaxis()->GetBinCenter(ip+1);
  // 	      Double_t r = p/(Const*1.);//1T

  // 	      //a*sin(theta) + b*cos(theta) +c = 0
  // 	      Double_t a = 2.*(r+rd)*y;
  // 	      Double_t b = 2.*(r+rd)*x;
  // 	      Double_t c = -1.*(rd*rd + 2.*r*rd + x*x + y*y);

  // 	      Double_t r0 = sqrt(a*a + b*b);
  // 	      if(fabs(-1.*c/r0)>1.){
  // 		// std::cout<<"No solution, "
  // 		//   <<"x:"<<x<<", y:"<<y
  // 		//   <<", r:"<<r<<", rd:"<<rd<<std::endl;
  // 		continue;
  // 	      }
  // 	      Double_t theta1_alpha =  asin(-1.*c/r0);
  // 	      Double_t theta2_alpha;
  // 	      if(theta1_alpha>0.)
  // 		theta2_alpha = acos(-1.) - theta1_alpha;
  // 	      else
  // 		theta2_alpha = -1.*acos(-1.) - theta1_alpha;

  // 	      Double_t theta_alpha = atan2(b, a);

  // 	      Double_t xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
  // 	      Double_t ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
  // 	      Double_t r_re1 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

  // 	      Double_t xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
  // 	      Double_t ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
  // 	      Double_t r_re2 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

  // 	      Double_t theta1 = atan2(ycenter1, xcenter1);
  // 	      Double_t theta2 = atan2(ycenter2, xcenter2);

  // 	      if(TMath::IsNaN(theta1)){
  // 		std::cout<<"theta1="<<theta1<<", x="<<x<<", y="<<y
  // 			 <<"rd="<<rd<<", r"<<r<<std::endl;
  // 	      }


  // 	      if(fabs(r-r_re1)>0.01||fabs(r-r_re2)>0.01){
  // 		std::cout<<"r="<<r<<", r_re1="<<r_re1<<", r_re1="<<r_re2<<std::endl;
  // 		std::cout<<"x:"<<x<<", y:"<<y
  // 			 <<", theta1:"<<theta1<<", theta2:"<<theta2
  // 			 <<", theta_alpha:"<<theta_alpha<<std::endl;
  // 	      }
  // 	      Ci_hist->Fill(rd, theta1, p);
  // 	      Ci_hist->Fill(rd, theta2, p);
  // 	      // std::cout<<"rd: "<<rd<<", "
  // 	      //         <<"theta1: "<<theta1<<", "
  // 	      //         <<"theta2: "<<theta2<<", "
  // 	      //         <<"p: "<<p<<std::endl;
  // 	    }
  // 	  }
  // 	}// cluster
  //     } // layer

  //     if(Ci_hist->GetMaximum() < MinNumOfHits/2){
  // 	Ci_hist->Reset();
  // 	break;
  //     }
  //     //std::cout<<"Maxbin0: "<<Ci_hist.GetMaximum()<<std::endl;


  //     Int_t maxbin = Ci_hist->GetMaximumBin();
  //     Int_t mx,my,mz;
  //     Ci_hist->GetBinXYZ(maxbin, mx, my, mz);

  //     Bool_t hough_flag = true;
  //     for(Int_t i=0; i<hough_x.size(); ++i){
  // 	Int_t bindiff = fabs(mx-hough_x[i])+fabs(my-hough_y[i]);
  // 	if(bindiff<=4)
  // 	  //if(bindiff<=8)
  // 	  hough_flag = false;
  //     }
  //     hough_x.push_back(mx);
  //     hough_y.push_back(my);
  //     if(!hough_flag)
  // 	continue;

  //     TPCLocalTrackHelix *track = new TPCLocalTrackHelix();
  //     hough_rd[tracki] = Ci_hist->GetXaxis()->GetBinCenter(mx);
  //     hough_theta[tracki] = Ci_hist->GetYaxis()->GetBinCenter(my);
  //     hough_p[tracki] = Ci_hist->GetZaxis()->GetBinCenter(mz);
  //     Double_t hough_r = hough_p[tracki]/Const;
  //     Double_t hough_cx = (hough_r + hough_rd[tracki])*cos(hough_theta[tracki]);
  //     Double_t hough_cy = (hough_r + hough_rd[tracki])*sin(hough_theta[tracki]);


  //     //std::cout<<"Hough (x,z) Maxbin: "<<Ci_hist.GetMaximum()<<std::endl;
  //     //     std::cout<<""<<std::endl;
  //     for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
  // 	for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
  // 	  //if(flag[layer][ci]>0) continue;
  // 	  auto cl = ClCont[layer][ci];
  // 	  TPCHit* hit = cl->GetMeanHit();
  // 	  if(hit->GetHoughFlag()>0) continue;
  // 	  TVector3 pos = cl->GetPosition();
  // 	  Double_t x = -pos.x();
  // 	  Double_t y = pos.z()-tpc::ZTarget;
  // 	  // Double_t de = hit->GetCharge();
  // 	  // Double_t xcenter1 = (hough_r + hough_rd[tracki])*cos(hough_theta[tracki]);
  // 	  // Double_t ycenter1 = (hough_r + hough_rd[tracki])*sin(hough_theta[tracki]);
  // 	  Double_t r_cal = sqrt(pow(x-hough_cx,2) + pow(y-hough_cy,2));

  // 	  Double_t dist = fabs(r_cal - hough_r);
  // 	  //if(dist < MaxHoughWindow && layer < MaxLayerCut && de>DECut_TPCTrack){
  // 	  if(dist < MaxHoughWindow && layer < MaxLayerCut){
  // 	    if(ity<Max_tracki_houghY){
  // 	      bool status_houghy = false;
  // 	      int ihoughy_size = hit->GetHoughY_num_size();
  // 	      for(int ih=0; ih<ihoughy_size; ++ih){
  // 		int ihoughy = hit->GetHoughY_num(ih);
  // 		//if(ity==ihoughy-1)
  // 		if(ity==ihoughy)
  // 		  status_houghy = true;
  // 	      }
  // 	      if(status_houghy)
  // 		track->AddTPCHit(new TPCLTrackHit(hit));
  // 	    }
  // 	    else
  // 	      track->AddTPCHit(new TPCLTrackHit(hit));
  // 	    //flag[layer][ci]++;
  // 	  }
  // 	}
  //     }
  //     // track->SetAdrho(hough_rd[tracki]);
  //     // track->SetAphi0(hough_theta[tracki]);
  //     // track->SetArho(1./hough_r);
  //     track->SetAcx(hough_cx);
  //     track->SetAcy(hough_cy);
  //     track->SetAr(hough_r);

  //     if(ity==0){
  // 	if(track->DoFit(3,1)){
  // 	  TrackCont.push_back(track);
  // 	  Int_t nh = track->GetNHit();
  // 	  for( int ih=0; ih<nh; ++ih ){
  // 	    TPCHit *hit = track->GetHit(ih)->GetHit();
  // 	    if(!hit) continue;
  // 	    hit->SetHoughFlag(1+hit->GetHoughFlag());
  // 	  }
  // 	}
  // 	else
  // 	  delete track;
  //     }
  //     else{
  // 	if(track->DoFit(MinNumOfHits,0)){
  // 	  TrackCont.push_back(track);
  // 	  Int_t nh = track->GetNHit();
  // 	  for( int ih=0; ih<nh; ++ih ){
  // 	    TPCHit *hit = track->GetHit(ih)->GetHit();
  // 	    if(!hit) continue;
  // 	    hit->SetHoughFlag(1+hit->GetHoughFlag());
  // 	  }
  // 	}
  // 	else
  // 	  delete track;
  //     }
  //     //      Ci_hist->Reset();

  //     //delete track;
  //   }//hough y cut
  // }//track
  // CalcTracks(TrackCont);
  // delete Ci_hist;
  // return status? TrackCont.size() : -1;

  // return 0;
}

}
