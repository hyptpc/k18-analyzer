// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define TrackSearch 0
#define Gain_center 1
#define DebugEvDisp 0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
TCanvas c1("c1", "c1", 800, 800);
#endif

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
//position cut for gain histogram
//  const double min_ycut = -15.;//mm
//const double max_ycut = 15.;//mm
const double min_ycut = -50.;//mm
const double max_ycut = 50.;//mm
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit,  kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]",  "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  Int_t runnum;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;
  Int_t nhTpc;
  Int_t nclTpc;
  std::vector<Double_t> xTpc;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Double_t> cluster_hitpos_x;
  std::vector<Double_t> cluster_hitpos_y;
  std::vector<Double_t> cluster_hitpos_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_size;
  std::vector<Int_t> cluster_layer;
  std::vector<Int_t> cluster_row;
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Double_t> cluster_hitpos_center_x;
  std::vector<Double_t> cluster_hitpos_center_y;
  std::vector<Double_t> cluster_hitpos_center_z;
  std::vector<Double_t> clkTpc;

  void clear()
    {
      runnum = 0;
      evnum = 0;
      status = 0;
      nhTpc = 0;
      nclTpc = 0;
      xTpc.clear();
      raw_hitpos_y.clear();
      raw_hitpos_z.clear();
      raw_de.clear();
      raw_padid.clear();
      cluster_hitpos_x.clear();
      cluster_hitpos_y.clear();
      cluster_hitpos_z.clear();
      cluster_de.clear();
      cluster_size.clear();
      cluster_layer.clear();
      cluster_row.clear();
      cluster_mrow.clear();
      cluster_de_center.clear();
      cluster_hitpos_center_x.clear();
      cluster_hitpos_center_y.clear();
      cluster_hitpos_center_z.clear();
      trigpat.clear();
      trigflag.clear();
      clkTpc.clear();
    }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;
  TTreeReaderValue<std::vector<Int_t>>* padTpc;
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;
  TTreeReaderValue<std::vector<Double_t>>* deTpc;
  TTreeReaderValue<std::vector<Double_t>>* tTpc;
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc;
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid {
  PadHid    = 100000,
};
}

//_____________________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  int open_file = 0;
  int open_tree = 0;
  for(Int_t i=0; i<nArgc; ++i){
    if(i==kProcess || i==kConfFile || i==kOutFile) continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=open_tree || open_file!=nArgc-3)
    return false;
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(int ievent)
{
  if(ievent%100==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  HC2Poly(10);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.nhTpc = **src.nhTpc;
  event.clkTpc = **src.clkTpc;

  HF1(1, event.status++);

  if(**src.nhTpc == 0)
    return true;

  HF1(1, event.status++);

  Double_t clock = event.clkTpc.at(0);

  DCAnalyzer DCAna;
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);

  Int_t nhTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = DCAna.GetTPCHC(layer);
    Int_t nh = 0;
    for(const auto& hit: hc){
      if(!hit || !hit->IsGood())
        continue;
      const auto& pos = hit->GetPosition();
      Double_t x = pos.X();
      Double_t y = pos.Y();
      Double_t z = pos.Z();
      Double_t de = hit->GetCDe();
      Double_t pad = hit->GetPad();
      event.xTpc.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);
      HF2Poly(10, z, x, de);
      HF1(101, de);
      HF1(1000*(layer)+101, de);
      HF1(102, y);
      HF1(1000*(layer)+102, y);
      ++nh;
    }
    nhTpc += nh;
    HF1(1000*(layer+1)+100, nh);
  }
  event.nhTpc = nhTpc;
  HF1(100, nhTpc);

#if DebugEvDisp
  TGraph g1;
  g1.SetMarkerStyle(5);
  g1.SetMarkerColor(kRed+1);
#endif

  Int_t nclTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const auto cont = DCAna.GetTPCClCont(layer);
    const Double_t R = tpc::GetRadius(layer);
    Int_t ncl = 0;
    for(const auto& cl: cont){
      if(!cl || !cl->IsGood()) continue;
      Double_t clde = cl->GetDe();
      // if(clde < 200) continue;
      // if(clde < 600) continue;
      auto& clpos = cl->GetPosition();
      Int_t cs = cl->GetClusterSize();
      // Double_t x = clpos.X();
      Double_t cly = clpos.Y();
      // Double_t z = clpos.Z();
      HF1(201, clde);
      HF1(1000*(layer)+201, clde);
      HF1(202, cly);
      HF1(1000*(layer)+202, cly);
      HF1(203, cs);
      HF1(1000*(layer)+203, cs);
      HF2(204, cs, clde);
      HF2(1000*(layer)+204, cs, clde);
      if(cs > 2){
        Double_t mtheta = tpc::getTheta(layer, cl->MeanRow());
        for(const auto& hit: cl->GetHitContainer()){
          Double_t theta = tpc::getTheta(layer, hit->GetRow());
          if(theta-mtheta > 180) theta -= 360.;
          Double_t Rt = R*(theta-mtheta)*TMath::DegToRad();
          Double_t de = hit->GetCDe();
          Double_t y = hit->GetY();
          HF2(205, Rt, de/clde);
          HF2(1000*(layer+1)+205, Rt, de/clde);
          HFProf(215, Rt, de/clde);
          HFProf(1000*(layer+1)+215, Rt, de/clde);
          HF2(206, Rt, y-cly);
          HF2(1000*(layer+1)+206, Rt, y-cly);
        }
      }
      ++ncl;
#if DebugEvDisp
      g1.SetPoint(g1.GetN(), clpos.Z(), clpos.X());
#endif
    }
    nclTpc += ncl;
    HF1(1000*(layer+1)+200, ncl);
  }
  HF1(200, nclTpc);

#if DebugEvDisp
  root::h[101]->Draw("colz");
  g1.Draw("P");
  gPad->SetLogz();
  gPad->Update();
  gPad->Modified();
  gSystem->ProcessEvents();
  c1.Print("c1.pdf");
  getchar();
  // ::sleep(1);
  gPad->Clear();
#endif

//   Int_t nh_cl_Tpc = 0;
//   for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
//     auto hc = DCAna.GetTPCClCont(layer);
//     for(const auto& hit : hc){
//       if(!hit || !hit->IsGood())
//         continue;
//       Double_t x = hit->GetX();
//       Double_t y = hit->GetY();
//       Double_t z = hit->GetZ();
//       Double_t de = hit->GetCharge();
//       Int_t cl_size = hit->GetClusterSize();
//       Int_t row = hit->GetRow();
//       Double_t mrow = hit->GetMRow();
//       Double_t de_center = hit->GetCharge_center();
//       TVector3 pos_center = hit->GetPos_center();
//       event.cluster_hitpos_x.push_back(x);
//       event.cluster_hitpos_y.push_back(y);
//       event.cluster_hitpos_z.push_back(z);
//       event.cluster_de.push_back(de);
//       event.cluster_size.push_back(cl_size);
//       event.cluster_layer.push_back(layer);
//       event.cluster_row.push_back(row);
//       event.cluster_mrow.push_back(mrow);
//       event.cluster_de_center.push_back(de_center);
//       event.cluster_hitpos_center_x.push_back(pos_center.X());
//       event.cluster_hitpos_center_y.push_back(pos_center.Y());
//       event.cluster_hitpos_center_z.push_back(pos_center.Z());
// #if Gain_center
//       //	if(69.<time&&time<85.&&nhit==1)
//       if(cl_size>1){
// 	if(min_ycut<y&&y<max_ycut)
// 	  HF1(PadHid + layer*1000 + row, de_center);
//       }
// #endif

//       ++nh_cl_Tpc;
//     }
//   }
//   event.nclTpc = nh_cl_Tpc;

  HF1(1, event.status++);

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for(Int_t i=0; i<n; ++i){
    if(TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 2000.;
  const Int_t    NbinClDe = 1000;
  const Double_t MinClDe  =    0.;
  const Double_t MaxClDe  = 2000.;
  const Int_t    NbinY =  800;
  const Double_t MinY  = -400.;
  const Double_t MaxY  =  400.;
  const Int_t    NbinCs = 40;
  const Double_t MinCs  =  0.;
  const Double_t MaxCs  = 40.;
  const Int_t    NbinRt = 100;
  const Double_t MinRt  = -25.;
  const Double_t MaxRt  =  25.;
  const Int_t    NbinDy = 200;
  const Double_t MinDy  = -10.;
  const Double_t MaxDy  =  10.;
  HB1(1, "Status", 21, 0., 21.);
  HB2Poly(10, "TPC EvDisp");
  HB1(100, "NHits TPC", 200, 0., 200.);
  HB1(101, "DeltaE TPCHit", NbinDe, MinDe, MaxDe);
  HB1(102, "Y TPCHit", NbinY, MinY, MaxY);
  HB1(200, "NClusters TPC", 200, 0., 200.);
  HB1(201, "DeltaE TPCCluster", NbinClDe, MinClDe, MaxClDe);
  HB1(202, "Y TPCCluster", NbinY, MinY, MaxY);
  HB1(203, "ClusterSize TPCCluster", NbinCs, MinCs, MaxCs);
  HB2(204, "DeltaE%ClusterSize TPCCluster",
      NbinCs, MinCs, MaxCs, NbinClDe/10, MinClDe, MaxClDe);
  HB2(205, "DeltaE/DeltaE_{total}%Rd#theta TPCCluster;[mm];",
      NbinRt, MinRt, MaxRt, 100, 0., 1.);
  HB2(206, "Y%Rd#theta TPCCluster;[mm];",
      NbinRt, MinRt, MaxRt, NbinDy, MinDy, MaxDy);
  HBProf(215, "DeltaE/DeltaE_{total}%Rd#theta (Prof) TPCCluster;[mm];",
         NbinRt, MinRt, MaxRt, 0., 1.);
  for(Int_t l=0; l<NumOfLayersTPC; ++l){
    HB1(1000*(l+1)+100, Form("NHits Layer#%d", l), 40, 0., 40.);
    HB1(1000*(l+1)+101, Form("DeltaE TPCHit Layer#%d", l), NbinDe, MinDe, MaxDe);
    HB1(1000*(l+1)+102, Form("Y TPCHit Layer#%d", l), NbinY, MinY, MaxY);
    HB1(1000*(l+1)+200, Form("NClusters Layer#%d", l), 40, 0., 40.);
    HB1(1000*(l+1)+201, Form("DeltaE TPCCluster Layer#%d", l),
        NbinClDe, MinClDe, MaxClDe);
    HB1(1000*(l+1)+202, Form("Y TPCCluster Layer#%d", l), NbinY, MinY, MaxY);
    HB1(1000*(l+1)+203, Form("ClusterSize TPCCluster Layer#%d", l),
        NbinCs, MinCs, MaxCs);
    HB2(1000*(l+1)+204, Form("DeltaE%%ClusterSize TPCCluster Layer#%d", l),
        NbinCs, MinCs, MaxCs, NbinClDe/10, MinClDe, MaxClDe);
    HB2(1000*(l+1)+205, Form("DeltaE/DeltaE_{total}%%Rd#theta TPCCluster Layer#%d;[mm];", l),
        NbinRt, MinRt, MaxRt, 100, 0., 1.);
    HB2(1000*(l+1)+206, Form("Y%%Rd#theta TPCCluster Layer#%d;[mm];", l),
        NbinRt, MinRt, MaxRt, NbinDy, MinDy, MaxDy);
    HBProf(1000*(l+1)+215, Form("DeltaE/DeltaE_{total}%%Rd#theta (Prof) TPCCluster Layer#%d;[mm];", l),
           NbinRt, MinRt, MaxRt, 0., 1.);
  }
  tpc::InitializeHistograms();

#if DebugEvDisp
  gStyle->SetOptStat(0);
  c1.cd();
#endif

  HBTree("tpc", "tree of DstTPCCluster");

  tree->Branch("status", &event.status);
  tree->Branch("runnum", &event.runnum);
  tree->Branch("evnum", &event.evnum);
  tree->Branch("trigpat", &event.trigpat);
  tree->Branch("trigflag", &event.trigflag);
  tree->Branch("nhTpc", &event.nhTpc);
  tree->Branch("nclTpc", &event.nclTpc);
  tree->Branch("xTpc", &event.xTpc);
  tree->Branch("raw_hitpos_y", &event.raw_hitpos_y);
  tree->Branch("raw_hitpos_z", &event.raw_hitpos_z);
  tree->Branch("raw_de", &event.raw_de);
  tree->Branch("raw_padid", &event.raw_padid);
  tree->Branch("cluster_hitpos_x", &event.cluster_hitpos_x);
  tree->Branch("cluster_hitpos_y", &event.cluster_hitpos_y);
  tree->Branch("cluster_hitpos_z", &event.cluster_hitpos_z);
  tree->Branch("cluster_de", &event.cluster_de);
  tree->Branch("cluster_size", &event.cluster_size);
  tree->Branch("cluster_layer", &event.cluster_layer);
  tree->Branch("cluster_row", &event.cluster_row);
  tree->Branch("cluster_mrow", &event.cluster_mrow);
  tree->Branch("cluster_de_center", &event.cluster_de_center);
  tree->Branch("cluster_hitpos_center_x", &event.cluster_hitpos_center_x);
  tree->Branch("cluster_hitpos_center_y", &event.cluster_hitpos_center_y);
  tree->Branch("cluster_hitpos_center_z", &event.cluster_hitpos_center_z);
  tree->Branch("clkTpc", &event.clkTpc);

// #if Gain_center
//   for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
//     const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
//     for(Int_t r=0; r<NumOfRow; ++r){
//       HB1(PadHid + layer*1000 + r , "TPC DeltaE_center", NbinDe, MinDe, MaxDe);
//     }
//   }
// #endif

  TTreeReaderCont[kTpcHit] = new TTreeReader("tpc", TFileCont[kTpcHit]);
  const auto& reader = TTreeReaderCont[kTpcHit];
  src.runnum = new TTreeReaderValue<Int_t>(*reader, "runnum");
  src.evnum = new TTreeReaderValue<Int_t>(*reader, "evnum");
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigpat");
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigflag");
  src.nhTpc = new TTreeReaderValue<Int_t>(*reader, "nhTpc");
  src.layerTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "layerTpc");
  src.rowTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "rowTpc");
  src.padTpc = new TTreeReaderValue<std::vector<Int_t>>(*reader, "padTpc");
  src.pedTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "pedTpc");
  src.rmsTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "rmsTpc");
  src.deTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "deTpc");
  src.tTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "tTpc");
  src.ctTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ctTpc");
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "chisqrTpc");
  src.clkTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "clkTpc");

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
