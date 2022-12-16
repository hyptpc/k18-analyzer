// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <TSystem.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "FieldMan.hh"
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

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

namespace
{
using namespace root;
using namespace dst;
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gTPCPositionCorrector  = TPCPositionCorrector::GetInstance();
const Double_t& zK18HS = gGeom.LocalZ("K18HS");

//For GenFit Setting
const bool Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");
const Int_t& pdgcode = ConfMan::Get<Int_t>("PDGcode");
const Double_t& PK18 = ConfMan::Get<Double_t>("PK18");
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit,  kK18Tracking, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]", "[K18Tracking]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc","k18track" ,"" };
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
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Double_t> cluster_hitpos_center_x;
  std::vector<Double_t> cluster_hitpos_center_y;
  std::vector<Double_t> cluster_hitpos_center_z;

  //BcOut info
  Int_t ntBcOut;
  std::vector<Double_t> chisqrBcOut;
  std::vector<Double_t> x0BcOut;
  std::vector<Double_t> y0BcOut;
  std::vector<Double_t> u0BcOut;
  std::vector<Double_t> v0BcOut;
  std::vector<Double_t> xBcOut; //HS coordinate
  std::vector<Double_t> yBcOut;
  std::vector<Double_t> uBcOut;
  std::vector<Double_t> vBcOut;

  // K18
  Int_t ntK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> xout;
  std::vector<Double_t> yout;
  std::vector<Double_t> uout;
  std::vector<Double_t> vout;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Double_t> chisqrTpc;
  std::vector<Double_t> x0Tpc;
  std::vector<Double_t> y0Tpc;
  std::vector<Double_t> u0Tpc;
  std::vector<Double_t> v0Tpc;
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
  std::vector<std::vector<Double_t>> residual_wbcout_x;
  std::vector<std::vector<Double_t>> residual_wbcout_y;
  std::vector<std::vector<Double_t>> residual_trackwbcout_x;
  std::vector<std::vector<Double_t>> residual_trackwbcout_y;

  std::vector<Double_t> xCorVec; // correction vector x
  std::vector<Double_t> yCorVec; // correction vector y
  std::vector<Double_t> zCorVec; // correction vector z
  std::vector<Double_t> xCorPos; // position x
  std::vector<Double_t> yCorPos; // position y
  std::vector<Double_t> zCorPos; // position z

  std::vector<Double_t> clkTpc;

  Int_t GFstatus;
  Int_t GFntTpc;
  std::vector<Int_t> GFfitstatus;
  std::vector<Int_t> GFndf;
  std::vector<Int_t> GFnhits;
  std::vector<Double_t> GFchisqr;
  std::vector<Double_t> GFcharge;
  std::vector<Double_t> GFtof;
  std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFpval;
  std::vector<Double_t> GFxCorVec; // Kalman filtered correction vector x
  std::vector<Double_t> GFyCorVec; // Kalman filtered correction vector y
  std::vector<Double_t> GFzCorVec; // Kalman filtered correction vector z
  std::vector<Double_t> GFxCorPos; // Kalman filtered position x
  std::vector<Double_t> GFyCorPos; // Kalman filtered position y
  std::vector<Double_t> GFzCorPos; // Kalman filtered position z
  std::vector<std::vector<Double_t>> GFpos_x;
  std::vector<std::vector<Double_t>> GFpos_y;
  std::vector<std::vector<Double_t>> GFpos_z;
  std::vector<std::vector<Double_t>> GFmom;
  std::vector<std::vector<Double_t>> GFmom_x;
  std::vector<std::vector<Double_t>> GFmom_y;
  std::vector<std::vector<Double_t>> GFmom_z;
  std::vector<std::vector<Double_t>> GFresidual_x;
  std::vector<std::vector<Double_t>> GFresidual_y;
  std::vector<std::vector<Double_t>> GFresidual_z;
  std::vector<std::vector<Double_t>> GFlayer;

  void clear()
  {
    runnum = 0;
    evnum = 0;
    status = 0;
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
    cluster_mrow.clear();
    cluster_de_center.clear();
    cluster_hitpos_center_x.clear();
    cluster_hitpos_center_y.clear();
    cluster_hitpos_center_z.clear();
    ntTpc = 0;
    trigpat.clear();
    trigflag.clear();
    clkTpc.clear();

    ntBcOut = 0;
    chisqrBcOut.clear();
    x0BcOut.clear();
    y0BcOut.clear();
    u0BcOut.clear();
    v0BcOut.clear();
    xBcOut.clear();
    yBcOut.clear();
    uBcOut.clear();
    vBcOut.clear();

    ntK18 = 0;
    chisqrK18.clear();
    pK18.clear();
    xout.clear();
    yout.clear();
    uout.clear();
    vout.clear();

    nhtrack.clear();
    chisqrTpc.clear();
    x0Tpc.clear();
    y0Tpc.clear();
    u0Tpc.clear();
    v0Tpc.clear();
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
    residual_wbcout_x.clear();
    residual_wbcout_y.clear();
    residual_trackwbcout_x.clear();
    residual_trackwbcout_y.clear();
    residual_z.clear();

    xCorVec.clear();
    yCorVec.clear();
    zCorVec.clear();
    xCorPos.clear();
    yCorPos.clear();
    zCorPos.clear();

    GFstatus = 0;
    GFntTpc = 0;
    GFchisqr.clear();
    GFcharge.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFfitstatus.clear();
    GFndf.clear();
    GFnhits.clear();
    GFpos_x.clear();
    GFpos_y.clear();
    GFpos_z.clear();
    GFmom.clear();
    GFmom_x.clear();
    GFmom_y.clear();
    GFmom_z.clear();
    GFresidual_x.clear();
    GFresidual_y.clear();
    GFresidual_z.clear();
    GFlayer.clear();
    GFxCorVec.clear();
    GFyCorVec.clear();
    GFzCorVec.clear();
    GFxCorPos.clear();
    GFyCorPos.clear();
    GFzCorPos.clear();

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
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;      // time

  //BcOut input
  Int_t ntBcOut;
  Double_t chisqrBcOut[MaxHits];
  Double_t x0BcOut[MaxHits];
  Double_t y0BcOut[MaxHits];
  Double_t u0BcOut[MaxHits];
  Double_t v0BcOut[MaxHits];

  // K18
  Int_t ntK18;
  Double_t chisqrK18[MaxHits];
  Double_t pK18[MaxHits];
  Double_t xout[MaxHits];
  Double_t yout[MaxHits];
  Double_t uout[MaxHits];
  Double_t vout[MaxHits];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
enum eDetHid {
  XCorrectionMapHid = 1000000,
  YCorrectionMapHid = 2000000,
  GFXCorrectionMapHid = 3000000,
  GFYCorrectionMapHid = 4000000,
  AfterXCorrectionMapHid = 5000000,
  AfterYCorrectionMapHid = 6000000,
  GFAfterXCorrectionMapHid = 7000000,
  GFAfterYCorrectionMapHid = 8000000,
  TPCDeHid   = 100000,
  TPCClDeHid = 200000,
  genfitHid = 300000
};

const Int_t MinPosMapXZ = -300;
const Int_t MaxPosMapXZ = 300;
const Int_t MinPosMapY = -200;
const Int_t MaxPosMapY = 200;
const Int_t Meshsize = 20;
const Int_t NumOfDivXZ = ((MaxPosMapXZ - MinPosMapXZ)/Meshsize) + 1;
const Int_t NumOfDivY = ((MaxPosMapY - MinPosMapY)/Meshsize) + 1;

//_____________________________________________________________________________
Int_t
XyzToHid(Double_t x, Double_t y, Double_t z)
{
  Int_t ix = TMath::Nint((x - MinPosMapXZ)/Meshsize);
  Int_t iy = TMath::Nint((y - MinPosMapY)/Meshsize);
  Int_t iz = TMath::Nint((z - MinPosMapXZ)/Meshsize);
  return ix*NumOfDivY*NumOfDivXZ + iy*NumOfDivXZ + iz;
}

//_____________________________________________________________________________
Int_t
XyzToHid(const TVector3& pos)
{
  return XyzToHid(pos.X(), pos.Y(), pos.Z());
}

//_____________________________________________________________________________
TVector3
HidToXyz(Int_t hid)
{
  Int_t ix = TMath::FloorNint(hid/NumOfDivY/NumOfDivXZ);
  Int_t iy = TMath::FloorNint((hid - ix*NumOfDivY*NumOfDivXZ)/NumOfDivXZ);
  Int_t iz = TMath::FloorNint((hid - ix*NumOfDivY*NumOfDivXZ - iy*NumOfDivXZ));
  return TVector3(ix*Meshsize + MinPosMapXZ,
                  iy*Meshsize + MinPosMapY,
                  iz*Meshsize + MinPosMapXZ);
}
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

  //Initiallize Geometry, Field, Fitter
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  GFtracks.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;

#if 0
  GFtracks.DebugMode();
#endif

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;
  DstClose();

  delete fitter;
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

  event.clkTpc = **src.clkTpc;

  event.ntBcOut = src.ntBcOut;
  for(int it=0; it<src.ntBcOut; ++it){
    event.chisqrBcOut.push_back(src.chisqrBcOut[it]);
    event.x0BcOut.push_back(src.x0BcOut[it]);
    event.y0BcOut.push_back(src.y0BcOut[it]);
    event.u0BcOut.push_back(src.u0BcOut[it]);
    event.v0BcOut.push_back(src.v0BcOut[it]);
    event.xBcOut.push_back(src.x0BcOut[it]+zK18HS*src.u0BcOut[it]);
    event.yBcOut.push_back(src.y0BcOut[it]+zK18HS*src.v0BcOut[it]);
    event.uBcOut.push_back(src.u0BcOut[it]);
    event.vBcOut.push_back(src.v0BcOut[it]);
  }

  event.ntK18 = src.ntK18;
  for(int it=0; it<src.ntK18; ++it){
    event.chisqrK18.push_back(src.chisqrK18[it]);
    event.pK18.push_back(src.pK18[it]);
    event.xout.push_back(src.xout[it]);
    event.yout.push_back(src.yout[it]);
    event.uout.push_back(src.uout[it]);
    event.vout.push_back(src.vout[it]);
  }

  HF1(1, event.status++);

  if(**src.nhTpc == 0)
    return true;

  HF1(1, event.status++);

  if(event.clkTpc.size() != 1){
    std::cerr << "something is wrong" << std::endl;
    return true;
  }
  Double_t clock = event.clkTpc.at(0);

  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  GFtracks.Init();

  HF1( 2, event.GFstatus++ );

  DCAnalyzer DCAna;
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);

  HF1(1, event.status++);

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

  HF1(1, event.status++);
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
      Double_t mrow = cl->MeanRow(); // same
      event.cluster_hitpos_x.push_back(x);
      event.cluster_hitpos_y.push_back(y);
      event.cluster_hitpos_z.push_back(z);
      event.cluster_de.push_back(clde);
      event.cluster_size.push_back(cs);
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(mrow);
      // event.cluster_de_center.push_back(de_center);

      ///// Compare with BcOut
      if(src.ntBcOut == 1){
	Double_t x0BcOut = src.x0BcOut[0];
	Double_t u0BcOut = src.u0BcOut[0];
	Double_t y0BcOut = src.y0BcOut[0];
	Double_t v0BcOut = src.v0BcOut[0];
	Double_t zTPC = zK18HS + z;
	Double_t xBcOut = x0BcOut + zTPC*u0BcOut;
	Double_t yBcOut = y0BcOut + zTPC*v0BcOut;
        Int_t hid = XyzToHid(x, y, z);
        event.xCorVec.push_back(xBcOut - x);
        event.yCorVec.push_back(yBcOut - y);
        event.zCorVec.push_back(0.);
        HF1(XCorrectionMapHid+hid, xBcOut - x);
        HF1(YCorrectionMapHid+hid, yBcOut - y);
	if(y<=0.5*Meshsize && y>-0.5*Meshsize){
	  HF2(XCorrectionMapHid+100000, z, x);
	  HF2(YCorrectionMapHid+100000, z, y);
	}

	TVector3 pos(x,y,z);
	TVector3 corPos = gTPCPositionCorrector.GetCorrectionVector(pos);
        event.xCorPos.push_back(corPos.x());
        event.yCorPos.push_back(corPos.y());
        event.zCorPos.push_back(corPos.z());
	HF1(AfterXCorrectionMapHid+hid, xBcOut - corPos.x());
        HF1(AfterYCorrectionMapHid+hid, yBcOut - corPos.y());
	if(corPos.y()<=0.5*Meshsize && corPos.y()>-0.5*Meshsize){
	  HF2(AfterXCorrectionMapHid+100000, corPos.z(), corPos.x());
	  HF2(AfterYCorrectionMapHid+100000, corPos.z(), corPos.y());
	}
      }
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
  event.chisqrTpc.resize(ntTpc);
  event.x0Tpc.resize(ntTpc);
  event.y0Tpc.resize(ntTpc);
  event.u0Tpc.resize(ntTpc);
  event.v0Tpc.resize(ntTpc);
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
  event.residual_wbcout_x.resize(ntTpc);
  event.residual_wbcout_y.resize(ntTpc);
  event.residual_trackwbcout_x.resize(ntTpc);
  event.residual_trackwbcout_y.resize(ntTpc);
  event.residual_z.resize(ntTpc);

  for(Int_t it=0; it<ntTpc; ++it){
    TPCLocalTrack *tp = DCAna.GetTrackTPC(it);
    if(!tp) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t x0Tpc=tp->GetX0(), y0Tpc=tp->GetY0();
    Double_t u0Tpc=tp->GetU0(), v0Tpc=tp->GetV0();
    Double_t theta = tp->GetTheta();
    HF1(11, nh);
    HF1(12, chisqr);
    HF1(14, x0Tpc);
    HF1(15, y0Tpc);
    HF1(16, u0Tpc);
    HF1(17, v0Tpc);
    HF2(18, x0Tpc, u0Tpc);
    HF2(19, y0Tpc, v0Tpc);
    HF2(20, x0Tpc, y0Tpc);
    if(src.ntBcOut==1){
      HF1(21, x0Tpc - event.xBcOut[0]);
      HF1(22, y0Tpc - event.yBcOut[0]);
      HF1(23, u0Tpc - event.uBcOut[0]);
      HF1(24, v0Tpc - event.vBcOut[0]);
    }
    event.nhtrack[it] = nh;
    event.chisqrTpc[it] = chisqr;
    event.x0Tpc[it] = x0Tpc;
    event.y0Tpc[it] = y0Tpc;
    event.u0Tpc[it] = u0Tpc;
    event.v0Tpc[it] = v0Tpc;
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
    event.residual_wbcout_x[it].resize(nh);
    event.residual_wbcout_y[it].resize(nh);
    event.residual_trackwbcout_x[it].resize(nh);
    event.residual_trackwbcout_y[it].resize(nh);
    event.residual_z[it].resize(nh);

    for(int ih=0; ih<nh; ++ih){
      TPCLTrackHit *hit = tp->GetHit(ih);
      if(!hit) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPos();
      const TVector3& res_vect = hit->GetResidualVect();
      Double_t residual = hit->GetResidual();
      HF1(13, layer);
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
      if(src.ntBcOut==1){
	Double_t x0BcOut = src.x0BcOut[0];
	Double_t u0BcOut = src.u0BcOut[0];
	Double_t y0BcOut = src.y0BcOut[0];
	Double_t v0BcOut = src.v0BcOut[0];
	Double_t zTPC = zK18HS + hitpos.z();
	Double_t xBcOut = x0BcOut + zTPC*u0BcOut;
	Double_t yBcOut = y0BcOut + zTPC*v0BcOut;
	event.residual_wbcout_x[it][ih] = hitpos.x() - xBcOut;
	event.residual_wbcout_y[it][ih] = hitpos.y() - yBcOut;
	zTPC = zK18HS + calpos.z();
	xBcOut = x0BcOut + zTPC*u0BcOut;
	yBcOut = y0BcOut + zTPC*v0BcOut;
	event.residual_trackwbcout_x[it][ih] = calpos.x() - xBcOut;
	event.residual_trackwbcout_y[it][ih] = calpos.y() - yBcOut;
      }
    }
  }

  if(event.ntK18!=1 || event.ntTpc!=1) return true;
  for(Int_t it=0; it<event.ntTpc; ++it){
    TPCLocalTrack *tp = DCAna.GetTrackTPC(it);
    if(!tp) continue;
    GFtracks.AddLinearTrack(pdgcode, tp, PK18);
  }
  GFtracks.FitTracks();

  int GFntTpc = GFtracks.GetNTrack();
  event.GFntTpc = GFntTpc;
  event.GFchisqr.resize(GFntTpc);
  event.GFcharge.resize(GFntTpc);
  event.GFtof.resize(GFntTpc);
  event.GFtracklen.resize(GFntTpc);
  event.GFpval.resize(GFntTpc);
  event.GFfitstatus.resize(GFntTpc);
  event.GFndf.resize(GFntTpc);
  event.GFnhits.resize(GFntTpc);
  event.GFmom.resize(GFntTpc);
  event.GFmom_x.resize(GFntTpc);
  event.GFmom_y.resize(GFntTpc);
  event.GFmom_z.resize(GFntTpc);
  event.GFmom.resize(GFntTpc);
  event.GFpos_x.resize(GFntTpc);
  event.GFpos_y.resize(GFntTpc);
  event.GFpos_z.resize(GFntTpc);
  event.GFresidual_x.resize(GFntTpc);
  event.GFresidual_y.resize(GFntTpc);
  event.GFresidual_z.resize(GFntTpc);
  event.GFlayer.resize(GFntTpc);
  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    event.GFfitstatus[igf] = (int)GFtracks.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
    if(!GFtracks.TrackCheck(igf)) continue;
    int nh = GFtracks.GetNHits(igf);
    event.GFmom_x[igf].resize(nh);
    event.GFmom_y[igf].resize(nh);
    event.GFmom_z[igf].resize(nh);
    event.GFmom[igf].resize(nh);
    event.GFpos_x[igf].resize(nh);
    event.GFpos_y[igf].resize(nh);
    event.GFpos_z[igf].resize(nh);
    event.GFresidual_x[igf].resize(nh);
    event.GFresidual_y[igf].resize(nh);
    event.GFresidual_z[igf].resize(nh);
    event.GFlayer[igf].resize(nh);
    for( Int_t ihit=0; ihit<nh; ++ihit ){
      TVector3 hit = GFtracks.GetPos(igf, ihit);
      TVector3 mom = GFtracks.GetMom(igf, ihit);
      event.GFmom_x[igf][ihit] = mom.x();
      event.GFmom_y[igf][ihit] = mom.y();
      event.GFmom_z[igf][ihit] = mom.z();
      event.GFmom[igf][ihit] = mom.Mag();
      event.GFpos_x[igf][ihit] = hit.x();
      event.GFpos_y[igf][ihit] = hit.y();
      event.GFpos_z[igf][ihit] = hit.z();
      event.GFresidual_x[igf][ihit] = hit.x() - event.hitpos_x[igf][ihit];
      event.GFresidual_y[igf][ihit] = hit.y() - event.hitpos_y[igf][ihit];
      event.GFresidual_z[igf][ihit] = hit.z() - event.hitpos_z[igf][ihit];
      event.GFlayer[igf][ihit] = event.hitlayer[igf][ihit];
      HF1( genfitHid+1000*(1+event.hitlayer[igf][ihit]), event.GFresidual_x[igf][ihit]);
      HF1( genfitHid+1000*(1+event.hitlayer[igf][ihit]+1), event.GFresidual_y[igf][ihit]);
      HF1( genfitHid+1000*(1+event.hitlayer[igf][ihit]+2), event.GFresidual_z[igf][ihit]);
      if(src.ntBcOut == 1){
	Double_t x0BcOut = src.x0BcOut[0];
	Double_t u0BcOut = src.u0BcOut[0];
	Double_t y0BcOut = src.y0BcOut[0];
	Double_t v0BcOut = src.v0BcOut[0];
	Double_t zTPC = zK18HS + hit.z();
	Double_t xBcOut = x0BcOut + zTPC*u0BcOut;
	Double_t yBcOut = y0BcOut + zTPC*v0BcOut;
        Int_t hid = XyzToHid(hit.x(), hit.y(), hit.z());
	event.GFxCorVec.push_back(xBcOut - hit.x());
	event.GFyCorVec.push_back(yBcOut - hit.y());
	event.GFzCorVec.push_back(0.);
	HF1(GFXCorrectionMapHid+hid, xBcOut - hit.x());
        HF1(GFYCorrectionMapHid+hid, yBcOut - hit.y());
	if(hit.y()<=0.5*Meshsize && hit.y()>-0.5*Meshsize){
	  HF2(GFXCorrectionMapHid+100000, hit.z(), hit.x());
	  HF2(GFYCorrectionMapHid+100000, hit.z(), hit.y());
	}
	TVector3 corPos = gTPCPositionCorrector.GetCorrectionVector(hit);
        event.GFxCorPos.push_back(corPos.x());
        event.GFyCorPos.push_back(corPos.y());
        event.GFzCorPos.push_back(corPos.z());
	HF1(GFAfterXCorrectionMapHid+hid, xBcOut - corPos.x());
        HF1(GFAfterYCorrectionMapHid+hid, yBcOut - corPos.y());
	if(corPos.y()<=0.5*Meshsize && corPos.y()>-0.5*Meshsize){
	  HF2(GFAfterXCorrectionMapHid+100000, corPos.z(), corPos.x());
	  HF2(GFAfterYCorrectionMapHid+100000, corPos.z(), corPos.y());
	}
      }
    } //ihit
    event.GFchisqr[igf] = GFtracks.GetChi2NDF(igf);
    event.GFcharge[igf] = GFtracks.GetCharge(igf);
    event.GFtof[igf] = GFtracks.GetTrackTOF(igf);
    event.GFtracklen[igf] = GFtracks.GetTrackLength(igf);
    event.GFpval[igf] = GFtracks.GetPvalue(igf);
    event.GFndf[igf] = GFtracks.GetNDF(igf);
    event.GFnhits[igf] = GFtracks.GetNHits(igf);
    HF1( genfitHid+1, event.GFchisqr[igf]);
    HF1( genfitHid+2, event.GFpval[igf]);
    HF1( genfitHid+3, event.GFcharge[igf]);
    HF1( genfitHid+4, event.GFnhits[igf]);

  } //igf
  HF1( genfitHid, event.GFntTpc );

  HF1(2, event.GFstatus++);

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
  const Int_t    NbinRes = 400;
  const Double_t MinRes  = -10.;
  const Double_t MaxRes  =  10.;
  const Int_t    NbinPos = 400;
  const Double_t MinPos  = -100.;
  const Double_t MaxPos  = 100.;

  HB1(1, "Status", 21, 0., 21.);
  HB1(2, "GFStatus", 21, 0., 21.);
  HB1(3, "Genfit Fit Status", 2, 0., 2. );
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
  HB2(20, "Y0%X0 TPC", 100, -100., 100., 100, -100, 100);
  HB1(21, "X0 TPC - BcOut", 400, -100., 100.);
  HB1(22, "Y0 TPC - BcOut", 400, -100., 100.);
  HB1(23, "U0 TPC - BcOut", 200, -0.20, 0.20);
  HB1(24, "V0 TPC - BcOut", 200, -0.20, 0.20);

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(100 + layer, "dE TPC", NbinDe, MinDe, MaxDe);
  }
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for(Int_t r=0; r<NumOfRow; ++r){
      HB1(TPCDeHid + layer*1000 + r , "TPC hit dE", NbinDe, MinDe, MaxDe);
      HB1(TPCClDeHid + layer*1000 + r , "TPC dE_center", NbinDe, MinDe, MaxDe);
    }
  }

  ///// Correction Map
  HB2(XCorrectionMapHid+100000, "Before TPC XCorrection at Y=0; Z[mm], X[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ);
  HB2(YCorrectionMapHid+100000, "Before TPC YCorrection at Y=0; Z[mm], Y[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivY, MinPosMapY, MaxPosMapY);
  HB2(AfterXCorrectionMapHid+100000, "After TPC XCorrection at Y=0; Z[mm], X[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ);
  HB2(AfterYCorrectionMapHid+100000, "After TPC YCorrection at Y=0; Z[mm], Y[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivY, MinPosMapY, MaxPosMapY);
  HB2(GFXCorrectionMapHid+100000, "[GenFit] Before TPC XCorrection at Y=0; Z[mm], X[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ);
  HB2(GFYCorrectionMapHid+100000, "[GenFit] Before TPC YCorrection at Y=0; Z[mm], Y[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivY, MinPosMapY, MaxPosMapY);
  HB2(GFAfterXCorrectionMapHid+100000, "[GenFit] After TPC XCorrection at Y=0; Z[mm], X[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ);
  HB2(GFAfterYCorrectionMapHid+100000, "[GenFit] AfterTPC YCorrection at Y=0; Z[mm], Y[mm]",
      NumOfDivXZ, MinPosMapXZ, MaxPosMapXZ,
      NumOfDivY, MinPosMapY, MaxPosMapY);

  Int_t hid_test = 0;
  for(Double_t x=MinPosMapXZ; x<=MaxPosMapXZ; x+=Meshsize){
    for(Double_t y=MinPosMapY; y<=MaxPosMapY; y+=Meshsize){
      for(Double_t z=MinPosMapXZ; z<=MaxPosMapXZ; z+=Meshsize){
        TVector3 pos(x, y, z);
        Int_t hid = XyzToHid(pos);
        if(pos != HidToXyz(hid)){ return false; } // check functions
        std::stringstream ss; ss << pos;
        HB1(XCorrectionMapHid+hid,
            "TPC XCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);
        HB1(YCorrectionMapHid+hid,
            "TPC YCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);
        HB1(GFXCorrectionMapHid+hid,
            "[GenFit] TPC XCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);
        HB1(GFYCorrectionMapHid+hid,
            "[GenFit] TPC YCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);
	HB1(AfterXCorrectionMapHid+hid,
            "After TPC XCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);
        HB1(AfterYCorrectionMapHid+hid,
            "After TPC YCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);
        HB1(GFAfterXCorrectionMapHid+hid,
            "After [GenFit] TPC XCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);
        HB1(GFAfterYCorrectionMapHid+hid,
            "After [GenFit] TPC YCorrection at "+ss.str(), NbinPos, MinPos, MaxPos);

        ++hid_test;
      }
    }
  }
  std::cout << hid_test << " passed" <<  std::endl;

  HB1(genfitHid, "[GenFit] #Track TPC", 10, 0., 10. );
  HB1(genfitHid+1, "[GenFit] Chisqr/ndf;", 200, 0, 10 );
  HB1(genfitHid+2, "[GenFit] p-value;p-value; Number of Tracks", 100, -0.05, 1.05);
  HB1(genfitHid+3, "[GenFit] Charge;", 6, -3, 3 );
  HB1(genfitHid+4, "[GenFit] #Hits of Track TPC", 33, 0., 33. );
  for(Int_t layer=1; layer<NumOfLayersTPC+1; ++layer){
    HB1(genfitHid+1000*layer, Form("[GenFit] TPC X Residual Layer%d (TPCHit);[mm];Counts", layer),NbinRes, MinRes, MaxRes);
    HB1(genfitHid+1000*layer+1, Form("[GenFit] TPC Y Residual Layer%d (TPCHit);[mm];Counts", layer),NbinRes, MinRes, MaxRes);
    HB1(genfitHid+1000*layer+2, Form("[GenFit] TPC Z Residual Layer%d (TPCHit);[mm];Counts", layer),NbinRes, MinRes, MaxRes);
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
  tree->Branch("cluster_mrow", &event.cluster_mrow);
  tree->Branch("cluster_de_center", &event.cluster_de_center);
  tree->Branch("cluster_hitpos_center_x", &event.cluster_hitpos_center_x);
  tree->Branch("cluster_hitpos_center_y", &event.cluster_hitpos_center_y);
  tree->Branch("cluster_hitpos_center_z", &event.cluster_hitpos_center_z);

  tree->Branch("ntBcOut", &event.ntBcOut);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);
  tree->Branch("xBcOut", &event.xBcOut);
  tree->Branch("yBcOut", &event.yBcOut);
  tree->Branch("uBcOut", &event.uBcOut);
  tree->Branch("vBcOut", &event.vBcOut);

  tree->Branch("ntK18",      &event.ntK18);
  tree->Branch("chisqrK18",  &event.chisqrK18);
  tree->Branch("pK18",       &event.pK18);
  tree->Branch("xout",    &event.xout);
  tree->Branch("yout",    &event.yout);
  tree->Branch("uout",    &event.uout);
  tree->Branch("vout",    &event.vout);

  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("nhtrack", &event.nhtrack);
  tree->Branch("chisqrTpc", &event.chisqrTpc);
  tree->Branch("x0Tpc", &event.x0Tpc);
  tree->Branch("y0Tpc", &event.y0Tpc);
  tree->Branch("u0Tpc", &event.u0Tpc);
  tree->Branch("v0Tpc", &event.v0Tpc);
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
  tree->Branch("residual_wbcout_x", &event.residual_wbcout_x);
  tree->Branch("residual_wbcout_y", &event.residual_wbcout_y);
  tree->Branch("residual_trackwbcout_x", &event.residual_trackwbcout_x);
  tree->Branch("residual_trackwbcout_y", &event.residual_trackwbcout_y);
  tree->Branch("residual_z", &event.residual_z);

  tree->Branch("xCorVec", &event.xCorVec);
  tree->Branch("yCorVec", &event.yCorVec);
  tree->Branch("zCorVec", &event.zCorVec);
  tree->Branch("xCorPos", &event.xCorPos);
  tree->Branch("yCorPos", &event.yCorPos);
  tree->Branch("zCorPos", &event.zCorPos);

  tree->Branch("GFstatus",&event.GFstatus);
  tree->Branch("GFntTpc",&event.GFntTpc);
  tree->Branch("GFchisqr",&event.GFchisqr);
  tree->Branch("GFpval",&event.GFpval);
  tree->Branch("GFfitstatus", &event.GFfitstatus);
  tree->Branch("GFmom",&event.GFmom);
  tree->Branch("GFmom_x",&event.GFmom_x);
  tree->Branch("GFmom_y",&event.GFmom_y);
  tree->Branch("GFmom_z",&event.GFmom_z);
  tree->Branch("GFpos_x",&event.GFpos_x);
  tree->Branch("GFpos_y",&event.GFpos_y);
  tree->Branch("GFpos_z",&event.GFpos_z);
  tree->Branch("GFresidual_x",&event.GFresidual_x);
  tree->Branch("GFresidual_y",&event.GFresidual_y);
  tree->Branch("GFresidual_z",&event.GFresidual_z);
  tree->Branch("GFlayer",&event.GFlayer);
  tree->Branch("GFcharge",&event.GFcharge);
  tree->Branch("GFtracklen",&event.GFtracklen);
  tree->Branch("GFtof",&event.GFtof);
  tree->Branch("GFndf",&event.GFndf);
  tree->Branch("GFnhits",&event.GFnhits);
  tree->Branch("GFxCorVec", &event.GFxCorVec);
  tree->Branch("GFyCorVec", &event.GFyCorVec);
  tree->Branch("GFzCorVec", &event.GFzCorVec);
  tree->Branch("GFxCorPos", &event.GFxCorPos);
  tree->Branch("GFyCorPos", &event.GFyCorPos);
  tree->Branch("GFzCorPos", &event.GFzCorPos);

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

  TTreeCont[kK18Tracking]->SetBranchStatus("*", 0);
  TTreeCont[kK18Tracking]->SetBranchStatus("ntBcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("chisqrBcOut", 1);
  TTreeCont[kK18Tracking]->SetBranchStatus("x0BcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("y0BcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("u0BcOut",     1);
  TTreeCont[kK18Tracking]->SetBranchStatus("v0BcOut",     1);

  TTreeCont[kK18Tracking]->SetBranchStatus("ntK18",       1);
  TTreeCont[kK18Tracking]->SetBranchStatus("chisqrK18",   1);
  TTreeCont[kK18Tracking]->SetBranchStatus("p_3rd",       1);
  TTreeCont[kK18Tracking]->SetBranchStatus("xout",        1);
  TTreeCont[kK18Tracking]->SetBranchStatus("yout",        1);
  TTreeCont[kK18Tracking]->SetBranchStatus("uout",        1);
  TTreeCont[kK18Tracking]->SetBranchStatus("vout",        1);

  TTreeCont[kK18Tracking]->SetBranchAddress("ntBcOut",     &src.ntBcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrBcOut", &src.chisqrBcOut);
  TTreeCont[kK18Tracking]->SetBranchAddress("x0BcOut",     &src.x0BcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("y0BcOut",     &src.y0BcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("u0BcOut",     &src.u0BcOut    );
  TTreeCont[kK18Tracking]->SetBranchAddress("v0BcOut",     &src.v0BcOut    );

  TTreeCont[kK18Tracking]->SetBranchAddress("ntK18",     &src.ntK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrK18", &src.chisqrK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("p_3rd",     &src.pK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("xout",   &src.xout);
  TTreeCont[kK18Tracking]->SetBranchAddress("yout",   &src.yout);
  TTreeCont[kK18Tracking]->SetBranchAddress("uout",   &src.uout);
  TTreeCont[kK18Tracking]->SetBranchAddress("vout",   &src.vout);

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
     InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
