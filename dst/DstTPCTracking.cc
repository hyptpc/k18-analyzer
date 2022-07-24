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
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define Gain_center 1
#define HoughYcut 0


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
  std::vector<Double_t> clkTpc;
  Int_t nhTpc;
  Int_t nclTpc;
  std::vector<Double_t> raw_hitpos_x;
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

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Double_t> chisqr;
  std::vector<Double_t> x0;
  std::vector<Double_t> y0;
  std::vector<Double_t> u0;
  std::vector<Double_t> v0;
  std::vector<Double_t> theta;
  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;

  void clear()
    {
      runnum = 0;
      evnum = 0;
      status = 0;
      clkTpc.clear();
      nhTpc = 0;
      nclTpc = 0;
      raw_hitpos_x.clear();
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
      ntTpc = 0;
      trigpat.clear();
      trigflag.clear();
      nhtrack.clear();
      chisqr.clear();
      x0.clear();
      y0.clear();
      u0.clear();
      v0.clear();
      theta.clear();
      hitlayer.clear();
      hitpos_x.clear();
      hitpos_y.clear();
      hitpos_z.clear();
      calpos_x.clear();
      calpos_y.clear();
      calpos_z.clear();
      residual.clear();
      residual_x.clear();
      residual_y.clear();
      residual_z.clear();
    }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;
  TTreeReaderValue<Int_t>* npadTpc;   // number of pads
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  // vector (size=nhTpc)
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;     // layer id
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;       // row id
  TTreeReaderValue<std::vector<Int_t>>* padTpc;       // pad id
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;    // pedestal
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;    // rms
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
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

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  HF1(1, event.status++);

  if(**src.nhTpc == 0)
    return true;
  event.clkTpc = **src.clkTpc;

  Double_t clock = event.clkTpc.at(0);

  DCAnalyzer DCAna;
#if 1
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);
#else
  static const Int_t exclusive_layer = gUser.GetParameter("TPCExclusive");
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock, exclusive_layer);
#endif
#if 0 //HoughYcut
  DCAna.HoughYCut(min_ycut, max_ycut);
#endif

  Int_t nhTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = DCAna.GetTPCHC(layer);
    for(const auto& hit : hc){
      if(!hit || !hit->IsGood())
	continue;
      const auto& pos = hit->GetPosition();
      Double_t x = pos.X();
      Double_t y = pos.Y();
      Double_t z = pos.Z();
      Double_t de = hit->GetCDe();
      Double_t pad = hit->GetPad();
      event.raw_hitpos_x.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);
      ++nhTpc;
    }
  }
  event.nhTpc = nhTpc;

  Int_t nclTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = DCAna.GetTPCClCont(layer);
    for(const auto& cl : hc){
      if(!cl || !cl->IsGood())
	continue;
      Double_t x = cl->GetX();
      Double_t y = cl->GetY();
      Double_t z = cl->GetZ();
      Double_t clde = cl->GetDe();
      Int_t cs = cl->GetClusterSize();
      Double_t row = cl->MeanRow();
      Double_t mrow = cl->MeanRow(); // same
      //       Double_t de_center = cl->GetDe_center();
      //       TVector3 pos_center = cl->GetPos_center();
      event.cluster_hitpos_x.push_back(x);
      event.cluster_hitpos_y.push_back(y);
      event.cluster_hitpos_z.push_back(z);
      event.cluster_de.push_back(clde);
      event.cluster_size.push_back(cs);
      event.cluster_layer.push_back(layer);
      event.cluster_row.push_back(row);
      event.cluster_mrow.push_back(mrow);
      //       event.cluster_de_center.push_back(de_center);
      //       event.cluster_hitpos_center_x.push_back(pos_center.X());
      //       event.cluster_hitpos_center_y.push_back(pos_center.Y());
      //       event.cluster_hitpos_center_z.push_back(pos_center.Z());
      // #if Gain_center
      //       //	if(69.<time&&time<85.&&nhit==1)
      //       if(cs>1){
      // 	if(min_ycut<y&&y<max_ycut)
      // 	  HF1(PadHid + layer*1000 + row, de_center);
      //       }
      // #endif
      ++nclTpc;
    }
  }
  event.nclTpc = nclTpc;

  HF1(1, event.status++);

  DCAna.TrackSearchTPC();

  Int_t ntTpc = DCAna.GetNTracksTPC();
  event.ntTpc = ntTpc;
  HF1(10, ntTpc);
  if(event.ntTpc == 0)
     return true;

  HF1(1, event.status++);

  event.nhtrack.resize(ntTpc);
  event.chisqr.resize(ntTpc);
  event.x0.resize(ntTpc);
  event.y0.resize(ntTpc);
  event.u0.resize(ntTpc);
  event.v0.resize(ntTpc);
  event.theta.resize(ntTpc);
  event.hitlayer.resize(ntTpc);
  event.hitpos_x.resize(ntTpc);
  event.hitpos_y.resize(ntTpc);
  event.hitpos_z.resize(ntTpc);
  event.calpos_x.resize(ntTpc);
  event.calpos_y.resize(ntTpc);
  event.calpos_z.resize(ntTpc);
  event.residual.resize(ntTpc);
  event.residual_x.resize(ntTpc);
  event.residual_y.resize(ntTpc);
  event.residual_z.resize(ntTpc);

  for(Int_t it=0; it<ntTpc; ++it){
    auto track = DCAna.GetTrackTPC(it);
    if(!track) continue;
    Int_t nh = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    Double_t x0=track->GetX0(), y0=track->GetY0();
    Double_t u0=track->GetU0(), v0=track->GetV0();
    Double_t theta = track->GetTheta();
    HF1(11, nh);
    HF1(12, chisqr);
    HF1(14, x0);
    HF1(15, y0);
    HF1(16, u0);
    HF1(17, v0);
    HF2(18, x0, u0);
    HF2(19, y0, v0);
    HF2(20, x0, y0);
    event.nhtrack[it] = nh;
    event.chisqr[it] = chisqr;
    event.x0[it] = x0;
    event.y0[it] = y0;
    event.u0[it] = u0;
    event.v0[it] = v0;
    event.theta[it] = theta;
    event.hitlayer[it].resize(nh);
    event.hitpos_x[it].resize(nh);
    event.hitpos_y[it].resize(nh);
    event.hitpos_z[it].resize(nh);
    event.calpos_x[it].resize(nh);
    event.calpos_y[it].resize(nh);
    event.calpos_z[it].resize(nh);
    event.residual[it].resize(nh);
    event.residual_x[it].resize(nh);
    event.residual_y[it].resize(nh);
    event.residual_z[it].resize(nh);

    for(int ih=0; ih<nh; ++ih){
      auto hit = track->GetHit(ih);
      if(!hit) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPos();
      const TVector3& res_vect = hit->GetResidualVect();
      Double_t residual = hit->GetResidual();
      HF1(13, layer);
      HF1(100*(layer+1)+15, residual);
      HF1(100*(layer+1)+31, res_vect.X());
      HF1(100*(layer+1)+32, res_vect.Y());
      HF1(100*(layer+1)+33, res_vect.Z());
      event.hitlayer[it][ih] = layer;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = res_vect.x();
      event.residual_y[it][ih] = res_vect.y();
      event.residual_z[it][ih] = res_vect.z();
    }

  }

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
  HB1(1, "Status", 21, 0., 21.);
  HB1(10, "#Tracks TPC", 40, 0., 40.);
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 500.);
  HB1(13, "LayerId TPC", 35, 0., 35.);
  HB1(14, "X0 TPC", 400, -100., 100.);
  HB1(15, "Y0 TPC", 400, -100., 100.);
  HB1(16, "U0 TPC", 200, -0.20, 0.20);
  HB1(17, "V0 TPC", 200, -0.20, 0.20);
  HB2(18, "U0%X0 TPC", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 TPC", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 TPC", 100, -100., 100., 100, -100, 100);

  for(Int_t i=1; i<=NumOfLayersTPC; ++i){
    // Tracking Histgrams
    TString title11 = Form("HitPat TPC%2d [Track]", i);
    TString title14 = Form("Position TPC%2d", i);
    TString title15 = Form("Residual TPC%2d", i);
    TString title16 = Form("Resid%%Pos TPC%2d", i);
    TString title17 = Form("Y%%Xcal TPC%2d", i);
    TString title31 = Form("ResidualX TPC%2d", i);
    TString title32 = Form("ResidualY TPC%2d", i);
    TString title33 = Form("ResidualZ TPC%2d", i);
    HB1(100*i+11, title11, 400, 0., 400.);
    HB1(100*i+14, title14, 200, -250., 250.);
    HB1(100*i+15, title15, 200, 0.0, 10.0);
    HB2(100*i+16, title16, 250, -250., 250., 100, -1.0, 1.0);
    HB2(100*i+17, title17, 100, -250., 250., 100, -250., 250.);
    HB1(100*i+31, title31, 200, -2.0, 2.0);
    HB1(100*i+32, title32, 200, -2.0, 2.0);
    HB1(100*i+33, title33, 200, -2.0, 2.0);
  }

  HBTree("tpc", "tree of DstTPCTracking");

  tree->Branch("status", &event.status);
  tree->Branch("runnum", &event.runnum);
  tree->Branch("evnum", &event.evnum);
  tree->Branch("trigpat", &event.trigpat);
  tree->Branch("trigflag", &event.trigflag);
  tree->Branch("clkTpc", &event.clkTpc);
  tree->Branch("nhTpc", &event.nhTpc);
  tree->Branch("nclTpc", &event.nclTpc);
  tree->Branch("raw_hitpos_x", &event.raw_hitpos_x);
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

  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("nhtrack", &event.nhtrack);
  tree->Branch("chisqr", &event.chisqr);
  tree->Branch("x0", &event.x0);
  tree->Branch("y0", &event.y0);
  tree->Branch("u0", &event.u0);
  tree->Branch("v0", &event.v0);
  tree->Branch("theta", &event.theta);
  tree->Branch("hitlayer", &event.hitlayer);
  tree->Branch("hitpos_x", &event.hitpos_x);
  tree->Branch("hitpos_y", &event.hitpos_y);
  tree->Branch("hitpos_z", &event.hitpos_z);
  tree->Branch("calpos_x", &event.calpos_x);
  tree->Branch("calpos_y", &event.calpos_y);
  tree->Branch("calpos_z", &event.calpos_z);
  tree->Branch("residual", &event.residual);
  tree->Branch("residual_x", &event.residual_x);
  tree->Branch("residual_y", &event.residual_y);
  tree->Branch("residual_z", &event.residual_z);

#if 0 //Gain_center
  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 2000.;

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for(Int_t r=0; r<NumOfRow; ++r){
      HB1(PadHid + layer*1000 + r , "TPC DeltaE_center", NbinDe, MinDe, MaxDe);
    }
  }
#endif



  TTreeReaderCont[kTpcHit] = new TTreeReader("tpc", TFileCont[kTpcHit]);
  const auto& reader = TTreeReaderCont[kTpcHit];
  src.runnum = new TTreeReaderValue<Int_t>(*reader, "runnum");
  src.evnum = new TTreeReaderValue<Int_t>(*reader, "evnum");
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigpat");
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigflag");
  src.npadTpc = new TTreeReaderValue<Int_t>(*reader, "npadTpc");
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
