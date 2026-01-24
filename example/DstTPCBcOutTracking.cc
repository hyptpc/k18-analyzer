// -*- C++ -*-
// To Do: consider BcOut nhit (not important)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <type_traits>

#include "TVector3.h"
#include "Math/Rotation3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Vector3D.h"

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
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

#define Exclusive 1

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf     = ConfMan::GetInstance();
const auto& gGeom     = DCGeomMan::GetInstance();
const auto& gUser     = UserParamMan::GetInstance();
const auto& gCounter  = debug::ObjectCounter::GetInstance();
const Double_t truncatedMean = 0.8; //80%
const Double_t MaxChisqrTpc = 50.0;
const Double_t MaxChisqrBcOut = 5.0;
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcTracking, kBcOut, kHodo, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCTracking]", "[BcOut]", "[Hodo]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "bcout", "hodo", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
std::vector<UInt_t> evnumPerFile;
Bool_t SetupReader();
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  UInt_t runnum;
  UInt_t evnum;
  std::vector<Double_t> trigpat;
  std::vector<std::vector<Double_t>> trigflag;
  Int_t beamflag;
  std::vector<Double_t> clkTpc;

  // TPC Tracking Info
  Int_t ntTpc;
  std::vector<Int_t> nhtrack;
  std::vector<Double_t> chisqrTpc, x0Tpc, y0Tpc, u0Tpc, v0Tpc, theta;
  std::vector<Double_t> dE, dEdx;
    
  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x, hitpos_y, hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x, calpos_y, calpos_z;
  std::vector<std::vector<Double_t>> residual, residual_x, residual_y, residual_z;
  std::vector<std::vector<Double_t>> residual_horizontal, residual_vertical;
  std::vector<std::vector<Double_t>> resolution_x, resolution_y, resolution_z;
  std::vector<std::vector<Double_t>> resolution_horizontal, resolution_vertical;
  std::vector<std::vector<Double_t>> pathhit, theta_diff;
  std::vector<std::vector<Double_t>> track_cluster_de, track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  std::vector<std::vector<Double_t>> exresidual, exresidual_x, exresidual_y, exresidual_z;
  std::vector<std::vector<Double_t>> exresidual_horizontal, exresidual_vertical;

  Int_t ntTpc_inside;
  Double_t prodvtx_x, prodvtx_y, prodvtx_z;

  // BcOut Info
  Int_t ntBcOut;
  // std::vector<Int_t> nhBcOut;
  std::vector<Double_t> chisqrBcOut;
  std::vector<Double_t> x0BcOut, y0BcOut, u0BcOut, v0BcOut;
  std::vector<Double_t> xtgtBcOut, ytgtBcOut, utgtBcOut, vtgtBcOut;

  // Hodo Info
  Double_t btof;
  Double_t ftof;

  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc);
  }
  
  void clearTPCTracks() {
    ntTpc = 0;
    dst::clear_all(
      nhtrack, chisqrTpc, 
      x0Tpc, y0Tpc, u0Tpc, v0Tpc, 
      theta, dE, dEdx,
      
      hitlayer, 
      hitpos_x, hitpos_y, hitpos_z, 
      calpos_x, calpos_y, calpos_z,
      
      residual, 
      residual_x, residual_y, residual_z, 
      residual_horizontal, residual_vertical,
      
      resolution_x, resolution_y, resolution_z, 
      resolution_horizontal, resolution_vertical,
      
      pathhit, theta_diff, track_cluster_de, track_cluster_size,
      track_cluster_mrow, track_cluster_de_center,
      track_cluster_x_center, track_cluster_y_center, track_cluster_z_center,
      track_cluster_row_center,

      exresidual, 
      exresidual_x, exresidual_y, exresidual_z, 
      exresidual_horizontal, exresidual_vertical
    );
  }

  void clearTPCVertex() {
    ntTpc_inside = 0;
    prodvtx_x = TMath::QuietNaN(); 
    prodvtx_y = TMath::QuietNaN(); 
    prodvtx_z = TMath::QuietNaN();
  }

  void clearBcOut() {
    ntBcOut = 0;
    dst::clear_all(
      // nhBcOut, 
      chisqrBcOut, 
      x0BcOut, y0BcOut, u0BcOut, v0BcOut, 
      xtgtBcOut, ytgtBcOut, utgtBcOut, vtgtBcOut
    );
  }

  void clear() {
    clearBasicInfo();
    clearTPCTracks();
    clearTPCVertex();
    clearBcOut();
    btof = TMath::QuietNaN();
    ftof = TMath::QuietNaN();
  }

};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<UInt_t>* runnum;
  TTreeReaderValue<UInt_t>* evnum_tpc;
  TTreeReaderValue<std::vector<Double_t>>* trigpat;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* trigflag;
  TTreeReaderValue<Int_t>* beamflag;
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;
  TTreeReaderValue<Int_t>* ntTpc;
  TTreeReaderValue<std::vector<Int_t>>* nhtrack;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc;
  TTreeReaderValue<std::vector<Double_t>>* x0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* y0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* u0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* v0Tpc;
  TTreeReaderValue<std::vector<Double_t>>* theta;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx;

  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitlayer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* hitpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* calpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_vertical;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_vertical;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* theta_diff;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_horizontal;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* exresidual_vertical;

  TTreeReaderValue<Int_t>* ntTpc_inside;
  TTreeReaderValue<Double_t>* prodvtx_x;
  TTreeReaderValue<Double_t>* prodvtx_y;
  TTreeReaderValue<Double_t>* prodvtx_z;

  TTreeReaderValue<Double_t>* btof;
  TTreeReaderValue<Double_t>* ftof;

  TTreeReaderValue<UInt_t>* evnum_bcout;
  TTreeReaderValue<UInt_t>* evnum_hodo;

  TTreeReaderValue<Int_t>* ntBcOut;
  // TTreeReaderValue<std::vector<Int_t>>* nhBcOut;
  TTreeReaderValue<std::vector<Double_t>>* chisqrBcOut;
  TTreeReaderValue<std::vector<Double_t>>* x0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* y0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* u0BcOut;
  TTreeReaderValue<std::vector<Double_t>>* v0BcOut;
};

//_____________________________________________________________________________
namespace root
{
Event event;
Src   src;
TTree *tree;

void CopyTPCTrackingData(Event& event, const Src& src)
{
  event.ntTpc     = **src.ntTpc;
  event.nhtrack   = **src.nhtrack;
  event.chisqrTpc = **src.chisqrTpc;
  event.x0Tpc     = **src.x0Tpc;
  event.y0Tpc     = **src.y0Tpc;
  event.u0Tpc     = **src.u0Tpc;
  event.v0Tpc     = **src.v0Tpc;
  event.theta     = **src.theta;
  event.dE        = **src.dE;
  event.dEdx      = **src.dEdx;

  event.hitlayer   = **src.hitlayer;
  event.hitpos_x   = **src.hitpos_x;
  event.hitpos_y   = **src.hitpos_y;
  event.hitpos_z   = **src.hitpos_z;
  event.calpos_x   = **src.calpos_x;
  event.calpos_y   = **src.calpos_y;
  event.calpos_z   = **src.calpos_z;
  event.residual   = **src.residual;
  event.residual_x = **src.residual_x;
  event.residual_y = **src.residual_y;
  event.residual_z = **src.residual_z;
  event.residual_horizontal = **src.residual_horizontal;
  event.residual_vertical   = **src.residual_vertical;
  event.resolution_x        = **src.resolution_x;
  event.resolution_y        = **src.resolution_y;
  event.resolution_z        = **src.resolution_z;
  event.resolution_horizontal = **src.resolution_horizontal;
  event.resolution_vertical   = **src.resolution_vertical;
  event.pathhit             = **src.pathhit;
  event.theta_diff          = **src.theta_diff;
  event.track_cluster_de    = **src.track_cluster_de;
  event.track_cluster_size  = **src.track_cluster_size;
  event.track_cluster_mrow       = **src.track_cluster_mrow;
  event.track_cluster_de_center  = **src.track_cluster_de_center;
  event.track_cluster_x_center   = **src.track_cluster_x_center;
  event.track_cluster_y_center   = **src.track_cluster_y_center;
  event.track_cluster_z_center   = **src.track_cluster_z_center;
  event.track_cluster_row_center = **src.track_cluster_row_center;

#if Exclusive
  event.exresidual   = **src.exresidual;
  event.exresidual_x = **src.exresidual_x;
  event.exresidual_y = **src.exresidual_y;
  event.exresidual_z = **src.exresidual_z;
  event.exresidual_horizontal = **src.exresidual_horizontal;
  event.exresidual_vertical   = **src.exresidual_vertical;
#endif

  event.ntTpc_inside = **src.ntTpc_inside;
  event.prodvtx_x    = **src.prodvtx_x;
  event.prodvtx_y    = **src.prodvtx_y;
  event.prodvtx_z    = **src.prodvtx_z;
}

void CopyBcOutData(Event& event, const Src& src)
{
  static const Double_t ztgtGlobal = gGeom.GetGlobalPosition("Target").z();

  // BcOut Reading
  const Int_t ntBcOut = **src.ntBcOut; 
  event.ntBcOut = ntBcOut;
  // event.nhBcOut = **src.nhBcOut;
  event.chisqrBcOut = **src.chisqrBcOut;
  event.x0BcOut     = **src.x0BcOut; 
  event.y0BcOut     = **src.y0BcOut;
  event.u0BcOut     = **src.u0BcOut; 
  event.v0BcOut     = **src.v0BcOut;
  event.utgtBcOut   = **src.u0BcOut; 
  event.vtgtBcOut   = **src.v0BcOut;

  event.xtgtBcOut.reserve(ntBcOut);
  event.ytgtBcOut.reserve(ntBcOut);
  for(Int_t it=0; it<ntBcOut; ++it){
    event.xtgtBcOut.push_back((**src.u0BcOut)[it]*ztgtGlobal + (**src.x0BcOut)[it]);
    event.ytgtBcOut.push_back((**src.v0BcOut)[it]*ztgtGlobal + (**src.y0BcOut)[it]);
  }
}
} // namespace root

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
  if(!gConf.InitializeHistograms())
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;
  if(!dst::SetupReader())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if(max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();
  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
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
  Int_t n_input_files = 0;
  for(const auto& name : TreeName) if(name != "") n_input_files++;

  Int_t open_file = 0, open_tree = 0;
  for(Int_t i=0; i<nArgc; ++i){
    if(TreeName[i] == "") continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=n_input_files || open_tree!=n_input_files){
    spdlog::error("DstOpen Failed: opened files/trees mismatch based on TreeName definitions."
                  " expected: {}, open_file: {}, open_tree: {}",
                  n_input_files, open_file, open_tree);
    return false;
  }
  if(!CheckEntries(TTreeCont)) return false;
  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");
  return true;
}


//_____________________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if(ievent%1000==0) {
    std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  evnumPerFile[0] = **src.evnum_tpc;
  evnumPerFile[1] = **src.evnum_bcout;
  evnumPerFile[2] = **src.evnum_hodo;
  if(!dst::CheckEventNumbers(evnumPerFile, ievent, 
                              {TreeName[kTpcTracking], TreeName[kBcOut], TreeName[kHodo]})){
    return false;
  }

  event.runnum   = **src.runnum;
  event.evnum    = **src.evnum_tpc;
  event.trigpat  = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.clkTpc   = **src.clkTpc;
  HF1("Status", event.status++);

  if(**src.ntTpc == 0) return true;
  HF1("Status", event.status++);

  CopyTPCTrackingData(event, src);
  HF1("Status", event.status++);

  CopyBcOutData(event, src);
  HF1("Status", event.status++);

  event.btof = **src.btof;
  event.ftof = **src.ftof;
  HF1("Status", event.status++);

  // Select best tracks (minimum chisqr) for correlation analysis
  Int_t bestTpcIdx = 0;
  Int_t bestBcOutIdx = 0;
  if(event.ntTpc > 1){
    Double_t minChisqrTpc = event.chisqrTpc[0];
    for(Int_t it=1; it<event.ntTpc; ++it){
      if(event.chisqrTpc[it] < minChisqrTpc){
        minChisqrTpc = event.chisqrTpc[it];
        bestTpcIdx = it;
      }
    }
  }
  if(event.ntBcOut > 1){
    Double_t minChisqrBcOut = event.chisqrBcOut[0];
    for(Int_t it=1; it<event.ntBcOut; ++it){
      if(event.chisqrBcOut[it] < minChisqrBcOut){
        minChisqrBcOut = event.chisqrBcOut[it];
        bestBcOutIdx = it;
      }
    }
  }

  // Require at least one track in each detector
  if(event.ntTpc == 0 || event.ntBcOut == 0) return true;
  HF1("Status", event.status++);

  // Check quality of selected tracks
  if(event.chisqrTpc[bestTpcIdx] > MaxChisqrTpc || event.chisqrBcOut[bestBcOutIdx] > MaxChisqrBcOut) return true;
  HF1("Status", event.status++);

  // TPCTracking
  HF1("Num_Track_TPC", event.ntTpc);
  for(Int_t it=0; it<event.ntTpc; ++it){
    HF1("Num_Track_TPC_Hits", event.nhtrack[it]);
    HF1("Chisqr_TPC", event.chisqrTpc[it]);
    HF1("X0_TPC", event.x0Tpc[it]);
    HF1("Y0_TPC", event.y0Tpc[it]);
    HF1("U0_TPC", event.u0Tpc[it]);
    HF1("V0_TPC", event.v0Tpc[it]);
    HF2("X0_vs_U0_TPC", event.x0Tpc[it], event.u0Tpc[it]);
    HF2("Y0_vs_V0_TPC", event.y0Tpc[it], event.v0Tpc[it]);
    HF2("X0_vs_Y0_TPC", event.x0Tpc[it], event.y0Tpc[it]);
    HF2("X0_vs_atanU0_TPC", event.x0Tpc[it], std::atan(event.u0Tpc[it])*TMath::RadToDeg());
    HF2("Y0_vs_atanV0_TPC", event.y0Tpc[it], std::atan(event.v0Tpc[it])*TMath::RadToDeg());
  }
  
  // BcOut
  const Int_t ntBcOut = **src.ntBcOut;
  HF1("Num_Track_BcOut", ntBcOut);
  for(Int_t it=0; it<ntBcOut; ++it){
    // HF1("Num_Track_BcOut_Hits", event.nhBcOut[it]);
    HF1("Chisqr_BcOut", event.chisqrBcOut[it]);
    HF1("X0_BcOut", event.x0BcOut[it]);
    HF1("Y0_BcOut", event.y0BcOut[it]);
    HF1("U0_BcOut", event.u0BcOut[it]);
    HF1("V0_BcOut", event.v0BcOut[it]);
    HF1("Xtgt_BcOut", event.xtgtBcOut[it]);
    HF1("Ytgt_BcOut", event.ytgtBcOut[it]);
    HF1("Utgt_BcOut", event.utgtBcOut[it]);
    HF1("Vtgt_BcOut", event.vtgtBcOut[it]);
    HF2("Xtgt_vs_Utgt_BcOut", event.xtgtBcOut[it], event.utgtBcOut[it]);
    HF2("Ytgt_vs_Vtgt_BcOut", event.ytgtBcOut[it], event.vtgtBcOut[it]);
    HF2("Xtgt_vs_Ytgt_BcOut", event.xtgtBcOut[it], event.ytgtBcOut[it]);
  }

  // Correlations
  // BcOut conversion
  std::vector<Double_t> xyuvBcOut_TPCcoor(4);
  {
    static const Double_t ztgtGlobal = gGeom.GetGlobalPosition("Target").z();
    static const Double_t RA1        = gGeom.GetRotAngle1("HypTPC");
    static const Double_t RA2        = gGeom.GetRotAngle2("HypTPC");
    static const Double_t tanRA1     = std::tan(RA1*TMath::DegToRad());
    static const Double_t tanRA2     = std::tan(RA2*TMath::DegToRad());

    ThreeVector globalpos(
      (**src.u0BcOut)[bestBcOutIdx]*ztgtGlobal + (**src.x0BcOut)[bestBcOutIdx],
      (**src.v0BcOut)[bestBcOutIdx]*ztgtGlobal + (**src.y0BcOut)[bestBcOutIdx],
      ztgtGlobal
    );

    // On the TPC coordinate
    ThreeVector localpos = gGeom.Global2LocalPos("HypTPC", globalpos);
    xyuvBcOut_TPCcoor[0] = localpos.x();
    xyuvBcOut_TPCcoor[1] = localpos.y();
    xyuvBcOut_TPCcoor[2] = ((**src.u0BcOut)[bestBcOutIdx] - tanRA2) / (1. + (**src.u0BcOut)[bestBcOutIdx]*tanRA2);
    xyuvBcOut_TPCcoor[3] = ((**src.v0BcOut)[bestBcOutIdx] + tanRA1) / (1. - (**src.v0BcOut)[bestBcOutIdx]*tanRA1);
  }
  const Double_t XtgtDiff = xyuvBcOut_TPCcoor[0]-event.x0Tpc[bestTpcIdx];
  const Double_t YtgtDiff = xyuvBcOut_TPCcoor[1]-event.y0Tpc[bestTpcIdx];
  const Double_t UtgtDiff = xyuvBcOut_TPCcoor[2]-event.u0Tpc[bestTpcIdx];
  const Double_t VtgtDiff = xyuvBcOut_TPCcoor[3]-event.v0Tpc[bestTpcIdx];

  HF2("Xtgt_BcOut_vs_Tpc", event.x0Tpc[bestTpcIdx], xyuvBcOut_TPCcoor[0]);
  HF2("Ytgt_BcOut_vs_Tpc", event.y0Tpc[bestTpcIdx], xyuvBcOut_TPCcoor[1]);
  HF2("Utgt_BcOut_vs_Tpc", event.u0Tpc[bestTpcIdx], xyuvBcOut_TPCcoor[2]);
  HF2("Vtgt_BcOut_vs_Tpc", event.v0Tpc[bestTpcIdx], xyuvBcOut_TPCcoor[3]);
  HF1("Xtgt_Diff", XtgtDiff);
  HF1("Ytgt_Diff", YtgtDiff);
  HF1("Utgt_Diff", UtgtDiff);
  HF1("Vtgt_Diff", VtgtDiff);
  HF2("Xtgt_Diff_vs_Xtgt_BcOut", xyuvBcOut_TPCcoor[0], XtgtDiff);
  HF2("Ytgt_Diff_vs_Ytgt_BcOut", xyuvBcOut_TPCcoor[1], YtgtDiff);
  HF2("Utgt_Diff_vs_Utgt_BcOut", xyuvBcOut_TPCcoor[2], UtgtDiff);
  HF2("Vtgt_Diff_vs_Vtgt_BcOut", xyuvBcOut_TPCcoor[3], VtgtDiff);

  // Residuals
  for(Int_t it=0; it<event.ntTpc; ++it){
    for(Int_t ih=0; ih<event.nhtrack[it]; ++ih){
      Int_t layer = static_cast<Int_t>(event.hitlayer[it][ih]);
      HF2("Layer_vs_ResY", event.hitlayer[it][ih], event.residual_y[it][ih]);
      if (event.resolution_x[it][ih] > 0.9e5 
          && event.resolution_y[it][ih] > 0.9e5 
          && event.resolution_z[it][ih] > 0.9e5) continue; // exclude dummy hits

      // TPC residual
      HF1(Form("TPC_Layer%02d_X_Residual", layer), event.residual_x[it][ih]);
      HF1(Form("TPC_Layer%02d_Y_Residual", layer), event.residual_y[it][ih]);
      HF1(Form("TPC_Layer%02d_Z_Residual", layer), event.residual_z[it][ih]);
      HF1(Form("TPC_Layer%02d_Local_X_Residual", layer), event.residual_horizontal[it][ih]);
      HF1(Form("TPC_Layer%02d_Local_Y_Residual", layer), event.residual_vertical[it][ih]);
      HF1(Form("TPC_Layer%02d_XZ_Residual", layer), 
          std::hypot(event.residual_x[it][ih], event.residual_z[it][ih]));

      // TPC pull
      HF1(Form("TPC_Layer%02d_X_Pull", layer),
          event.residual_x[it][ih]/event.resolution_x[it][ih]);
      HF1(Form("TPC_Layer%02d_Y_Pull", layer),
          event.residual_y[it][ih]/event.resolution_y[it][ih]);
      HF1(Form("TPC_Layer%02d_Z_Pull", layer),
          event.residual_z[it][ih]/event.resolution_z[it][ih]);
      HF1(Form("TPC_Layer%02d_Local_X_Pull", layer), 
          event.residual_horizontal[it][ih]/event.resolution_horizontal[it][ih]);
      HF1(Form("TPC_Layer%02d_Local_Y_Pull", layer),
          event.residual_vertical[it][ih]/event.resolution_vertical[it][ih]);

      // TPC - BcOut residual
      ThreeVector localposTpc(
        event.hitpos_x[it][ih], event.hitpos_y[it][ih], event.hitpos_z[it][ih]
      );
      ThreeVector globalposTpc = gGeom.Local2GlobalPos("HypTPC", localposTpc);
      Double_t resX_BcOut = globalposTpc.x() - ((**src.u0BcOut)[it]*globalposTpc.z() + (**src.x0BcOut)[it]);
      Double_t resY_BcOut = globalposTpc.y() - ((**src.v0BcOut)[it]*globalposTpc.z() + (**src.y0BcOut)[it]);
      HF1(Form("TPC_Layer%02d_BcOut_X_Residual", layer), resX_BcOut);
      HF1(Form("TPC_Layer%02d_BcOut_Y_Residual", layer), resY_BcOut);
      // Parameter tuning: Position correction (BcOut reference)
      HF2(Form("TPC_ResidualX_vs_X_BcOut_Layer%02d", layer), globalposTpc.x(), resX_BcOut);
      HF2(Form("TPC_ResidualY_vs_Y_BcOut_Layer%02d", layer), globalposTpc.y(), resY_BcOut);
    }
  }
  HF1("Status", event.status++);

  return true;
}


//_____________________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  spdlog::info(" Close : {}", TFileCont[kOutFile]->GetName());
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
dst::SetupReader()
{
  // -------------------------------------------------------
  // TPCTracking
  // -------------------------------------------------------
  if (!dst::SetupReader(kTpcTracking, "kTpcTracking")) return false;

  // Explicitly using TTreeReaderCont[kTpcTracking]
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "run_number",   src.runnum);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "event_number", src.evnum_tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "trig_pat",     src.trigpat);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "trig_flag",    src.trigflag);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "beam_flag",    src.beamflag);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "clkTpc",       src.clkTpc);

  dst::SetBranch(TTreeReaderCont[kTpcTracking], "ntTpc",        src.ntTpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "nhtrack",      src.nhtrack);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "chisqrTpc",    src.chisqrTpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "x0Tpc",        src.x0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "y0Tpc",        src.y0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "u0Tpc",        src.u0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "v0Tpc",        src.v0Tpc);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "theta",        src.theta);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "dE",           src.dE);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "dEdx",         src.dEdx);

  // Hit & Cal Pos
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitlayer",     src.hitlayer);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitpos_x",     src.hitpos_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitpos_y",     src.hitpos_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "hitpos_z",     src.hitpos_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "calpos_x",     src.calpos_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "calpos_y",     src.calpos_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "calpos_z",     src.calpos_z);

  // Residuals
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual",            src.residual);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_x",          src.residual_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_y",          src.residual_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_z",          src.residual_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_horizontal", src.residual_horizontal);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "residual_vertical",   src.residual_vertical);

  // Resolution
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_x",          src.resolution_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_y",          src.resolution_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_z",          src.resolution_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_horizontal", src.resolution_horizontal);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "resolution_vertical",   src.resolution_vertical);

  // Track Clusters
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "pathhit",                  src.pathhit);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "theta_diff",               src.theta_diff);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_de",         src.track_cluster_de);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_size",       src.track_cluster_size);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_mrow",       src.track_cluster_mrow);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_de_center",  src.track_cluster_de_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_x_center",   src.track_cluster_x_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_y_center",   src.track_cluster_y_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_z_center",   src.track_cluster_z_center);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "track_cluster_row_center", src.track_cluster_row_center);

  // ExResiduals
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual",            src.exresidual);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_x",          src.exresidual_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_y",          src.exresidual_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_z",          src.exresidual_z);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_horizontal", src.exresidual_horizontal);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "exresidual_vertical",   src.exresidual_vertical);

  dst::SetBranch(TTreeReaderCont[kTpcTracking], "ntTpc_target", src.ntTpc_inside);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "prodvtx_x",    src.prodvtx_x);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "prodvtx_y",    src.prodvtx_y);
  dst::SetBranch(TTreeReaderCont[kTpcTracking], "prodvtx_z",    src.prodvtx_z);

  // -------------------------------------------------------
  // Hodoscope
  // -------------------------------------------------------
  if (!dst::SetupReader(kHodo, "kHodo")) return false;

  dst::SetBranch(TTreeReaderCont[kHodo], "event_number", src.evnum_hodo);
  dst::SetBranch(TTreeReaderCont[kHodo], "btof0", src.btof);
  dst::SetBranch(TTreeReaderCont[kHodo], "ftof0", src.ftof);

  // -------------------------------------------------------
  // BcOut
  // -------------------------------------------------------
  if (!dst::SetupReader(kBcOut, "kBcOut")) return false;

  dst::SetBranch(TTreeReaderCont[kBcOut], "event_number", src.evnum_bcout);
  dst::SetBranch(TTreeReaderCont[kBcOut], "ntrack", src.ntBcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "chisqr", src.chisqrBcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "x0",     src.x0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "y0",     src.y0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "u0",     src.u0BcOut);
  dst::SetBranch(TTreeReaderCont[kBcOut], "v0",     src.v0BcOut);

  evnumPerFile.resize(3);  // [0]=TPC, [1]=BcOut, [2]=Hodo
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  hist::BuildStatus();
  hist::BuildTPCBasic();
  hist::BuildTPCBcOutTracking();

  tree = new TTree("tpc", "tree of DstTPCBcOutTracking");
  tree->Branch("status",       &event.status);
  tree->Branch("run_number",   &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat",     &event.trigpat);
  tree->Branch("trig_flag",    &event.trigflag);
  tree->Branch("beam_flag",    &event.beamflag);
  tree->Branch("clkTpc",       &event.clkTpc);

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
  tree->Branch("residual_z", &event.residual_z);
  tree->Branch("residual_horizontal", &event.residual_horizontal);
  tree->Branch("residual_vertical", &event.residual_vertical);
  tree->Branch("resolution_x", &event.resolution_x);
  tree->Branch("resolution_y", &event.resolution_y);
  tree->Branch("resolution_z", &event.resolution_z);
  tree->Branch("resolution_horizontal", &event.resolution_horizontal);
  tree->Branch("resolution_vertical", &event.resolution_vertical);
  tree->Branch("dE", &event.dE );
  tree->Branch("dEdx", &event.dEdx );
  tree->Branch("pathhit", &event.pathhit);
  tree->Branch("theta_diff", &event.theta_diff);
  tree->Branch("track_cluster_de", &event.track_cluster_de);
  tree->Branch("track_cluster_size", &event.track_cluster_size);
  tree->Branch("track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch("track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch("track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch("track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch("track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch("track_cluster_row_center", &event.track_cluster_row_center);

#if Exclusive
  tree->Branch("exresidual", &event.exresidual);
  tree->Branch("exresidual_x", &event.exresidual_x);
  tree->Branch("exresidual_y", &event.exresidual_y);
  tree->Branch("exresidual_z", &event.exresidual_z);
  tree->Branch("exresidual_horizontal", &event.exresidual_horizontal);
  tree->Branch("exresidual_vertical", &event.exresidual_vertical);
#endif

  tree->Branch("btof", &event.btof);
  tree->Branch("ftof", &event.ftof);

  tree->Branch("ntBcOut", &event.ntBcOut);
  // tree->Branch("nhBcOut", &event.nhBcOut);
  tree->Branch("chisqrBcOut", &event.chisqrBcOut);
  tree->Branch("x0BcOut", &event.x0BcOut);
  tree->Branch("y0BcOut", &event.y0BcOut);
  tree->Branch("u0BcOut", &event.u0BcOut);
  tree->Branch("v0BcOut", &event.v0BcOut);
  tree->Branch("xtgtBcOut", &event.xtgtBcOut);
  tree->Branch("ytgtBcOut", &event.ytgtBcOut);
  tree->Branch("utgtBcOut", &event.utgtBcOut);
  tree->Branch("vtgtBcOut", &event.vtgtBcOut);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<TPCParamMan>("TPCPRM") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
