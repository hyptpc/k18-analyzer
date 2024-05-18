// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "FieldMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertexHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#define RawHit 0
#define TrackCluster 1

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& zHSCenter = gGeom.LocalZ("HS");
const Double_t truncatedMean = 0.8; //80%

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const Double_t vtx_scan_range = 150.;
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit, kHodoscope, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[TPCHit]", "[Hodoscope]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "hodo", "" };
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
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;

  Int_t nclTpc;
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_size;
  std::vector<Int_t> cluster_layer;
  std::vector<Double_t> cluster_mrow;
  std::vector<Double_t> cluster_de_center;
  std::vector<Double_t> cluster_x_center;
  std::vector<Double_t> cluster_y_center;
  std::vector<Double_t> cluster_z_center;
  std::vector<Int_t> cluster_row_center;
  std::vector<Int_t> cluster_houghflag;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> isBeam; // isBeam: 1 = Beam, 0 = Scat
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> fittime;  //usec
  std::vector<Int_t> searchtime; //usec
  std::vector<Int_t> niteration; //usec
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dz_factor;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx;
  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0
  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path

  std::vector<Int_t> pid;
  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
  std::vector<std::vector<Double_t>> calpos_x;
  std::vector<std::vector<Double_t>> calpos_y;
  std::vector<std::vector<Double_t>> calpos_z;
  std::vector<std::vector<Double_t>> mom_x;
  std::vector<std::vector<Double_t>> mom_y;
  std::vector<std::vector<Double_t>> mom_z;
  std::vector<std::vector<Double_t>> residual;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> alpha;
  std::vector<std::vector<Double_t>> houghflag;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  Int_t nvtxTpc;
  std::vector<Double_t> vtx_x;
  std::vector<Double_t> vtx_y;
  std::vector<Double_t> vtx_z;
  std::vector<Double_t> vtx_dist;
  std::vector<Double_t> vtx_angle;
  std::vector<std::vector<Double_t>> vtxid;
  std::vector<std::vector<Double_t>> vtxmom_theta;
  std::vector<std::vector<Double_t>> vtxpos_x;
  std::vector<std::vector<Double_t>> vtxpos_y;
  std::vector<std::vector<Double_t>> vtxpos_z;
  std::vector<std::vector<Double_t>> vtxmom_x;
  std::vector<std::vector<Double_t>> vtxmom_y;
  std::vector<std::vector<Double_t>> vtxmom_z;

  Int_t nhHtof;
  std::vector<Int_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;

  Int_t GFstatus;
  Int_t GFntTpc;
  std::vector<Int_t> GFfitstatus;
  std::vector<Int_t> GFpdgcode;
  std::vector<Int_t> GFnhtrack;
  std::vector<Double_t> GFcharge;
  std::vector<Double_t> GFchisqr;
  std::vector<Double_t> GFtof;
  std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFpval;
  std::vector<std::vector<Double_t>> GFlayer;
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
  std::vector<std::vector<Double_t>> GFresidual_p;
  std::vector<std::vector<Double_t>> GFresidual_px;
  std::vector<std::vector<Double_t>> GFresidual_py;
  std::vector<std::vector<Double_t>> GFresidual_pz;

  std::vector<Int_t> GFinside;
  std::vector<Double_t> GFxTgt;
  std::vector<Double_t> GFyTgt;
  std::vector<Double_t> GFzTgt;
  std::vector<Double_t> GFmomTgt;
  std::vector<Double_t> GFmomxTgt;
  std::vector<Double_t> GFmomyTgt;
  std::vector<Double_t> GFmomzTgt;

  std::vector<Double_t> GFtracklenTgt; //extrapolate to the target
  std::vector<Double_t> GFtofTgt;

  Int_t GFnvtxTpc;
  std::vector<Double_t> GFvtx_x;
  std::vector<Double_t> GFvtx_y;
  std::vector<Double_t> GFvtx_z;
  std::vector<Double_t> GFvtx_dist;
  std::vector<Double_t> GFvtx_angle;
  std::vector<std::vector<Double_t>> GFvtxid;
  std::vector<std::vector<Double_t>> GFvtxmom_x;
  std::vector<std::vector<Double_t>> GFvtxmom_y;
  std::vector<std::vector<Double_t>> GFvtxmom_z;

  std::vector<Int_t> GFnhHtof; //not real hit, extrapolation candidate.
  std::vector<std::vector<Double_t>> GFsegHtof;
  std::vector<std::vector<Double_t>> GFxHtof;
  std::vector<std::vector<Double_t>> GFyHtof;
  std::vector<std::vector<Double_t>> GFzHtof;
  std::vector<std::vector<Double_t>> GFtracklenHtof; //extrapolate to the htof
  std::vector<std::vector<Double_t>> GFtofHtof;

  std::vector<Int_t> GFfitstatus_p;
  std::vector<Double_t> GFmom_p;
  std::vector<Double_t> GFtof_p;
  std::vector<Double_t> GFtracklen_p;

  std::vector<Int_t> GFinside_p;
  std::vector<Double_t> GFxTgt_p;
  std::vector<Double_t> GFyTgt_p;
  std::vector<Double_t> GFzTgt_p;
  std::vector<Double_t> GFmomTgt_p;
  std::vector<Double_t> GFmomxTgt_p;
  std::vector<Double_t> GFmomyTgt_p;
  std::vector<Double_t> GFmomzTgt_p;
  std::vector<Double_t> GFtracklenTgt_p; //extrapolate to the target
  std::vector<Double_t> GFtofTgt_p;

  std::vector<Int_t> GFnhHtof_p;
  std::vector<std::vector<Double_t>> GFsegHtof_p;
  std::vector<std::vector<Double_t>> GFtracklenHtof_p;
  std::vector<std::vector<Double_t>> GFtofHtof_p;

  std::vector<Int_t> GFfitstatus_pi;
  std::vector<Double_t> GFmom_pi;
  std::vector<Double_t> GFtof_pi;
  std::vector<Double_t> GFtracklen_pi;

  std::vector<Int_t> GFinside_pi;
  std::vector<Double_t> GFxTgt_pi;
  std::vector<Double_t> GFyTgt_pi;
  std::vector<Double_t> GFzTgt_pi;
  std::vector<Double_t> GFmomTgt_pi;
  std::vector<Double_t> GFmomxTgt_pi;
  std::vector<Double_t> GFmomyTgt_pi;
  std::vector<Double_t> GFmomzTgt_pi;
  std::vector<Double_t> GFtracklenTgt_pi; //extrapolate to the target
  std::vector<Double_t> GFtofTgt_pi;

  std::vector<Int_t> GFnhHtof_pi;
  std::vector<std::vector<Double_t>> GFsegHtof_pi;
  std::vector<std::vector<Double_t>> GFtracklenHtof_pi;
  std::vector<std::vector<Double_t>> GFtofHtof_pi;

  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();
    clkTpc.clear();

    nhTpc = 0;
    raw_hitpos_x.clear();
    raw_hitpos_y.clear();
    raw_hitpos_z.clear();
    raw_de.clear();
    raw_padid.clear();
    raw_layer.clear();
    raw_row.clear();

    nclTpc = 0;
    cluster_x.clear();
    cluster_y.clear();
    cluster_z.clear();
    cluster_de.clear();
    cluster_size.clear();
    cluster_layer.clear();
    cluster_mrow.clear();
    cluster_de_center.clear();
    cluster_x_center.clear();
    cluster_y_center.clear();
    cluster_z_center.clear();
    cluster_row_center.clear();
    cluster_houghflag.clear();

    ntTpc = 0;
    nhtrack.clear();
    isBeam.clear();
    isAccidental.clear();
    fittime.clear();
    searchtime.clear();
    niteration.clear();
    chisqr.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();

    dz_factor.clear();
    mom0_x.clear();
    mom0_y.clear();
    mom0_z.clear();
    mom0.clear();
    charge.clear();
    path.clear();
    pid.clear();

    hitlayer.clear();
    hitpos_x.clear();
    hitpos_y.clear();
    hitpos_z.clear();
    calpos_x.clear();
    calpos_y.clear();
    calpos_z.clear();
    mom_x.clear();
    mom_y.clear();
    mom_z.clear();
    residual.clear();
    residual_x.clear();
    residual_y.clear();
    residual_z.clear();
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    helix_t.clear();
    pathhit.clear();
    alpha.clear();
    houghflag.clear();

    track_cluster_de.clear();
    track_cluster_size.clear();
    track_cluster_mrow.clear();
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();

    nvtxTpc = 0;
    vtx_x.clear();
    vtx_y.clear();
    vtx_z.clear();
    vtx_dist.clear();
    vtx_angle.clear();
    vtxid.clear();
    vtxmom_theta.clear();
    vtxpos_x.clear();
    vtxpos_y.clear();
    vtxpos_z.clear();
    vtxmom_x.clear();
    vtxmom_y.clear();
    vtxmom_z.clear();

    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    dtHtof.clear();
    deHtof.clear();
    posHtof.clear();

    GFstatus = 0;
    GFntTpc = 0;
    GFcharge.clear();
    GFchisqr.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFfitstatus.clear();
    GFpdgcode.clear();
    GFnhtrack.clear();
    GFlayer.clear();
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
    GFresidual_p.clear();
    GFresidual_px.clear();
    GFresidual_py.clear();
    GFresidual_pz.clear();

    GFinside.clear();
    GFxTgt.clear();
    GFyTgt.clear();
    GFzTgt.clear();
    GFmomxTgt.clear();
    GFmomxTgt.clear();
    GFmomyTgt.clear();
    GFmomzTgt.clear();
    GFtracklenTgt.clear();
    GFtofTgt.clear();

    GFnvtxTpc = -1;
    GFvtx_x.clear();
    GFvtx_y.clear();
    GFvtx_z.clear();
    GFvtx_dist.clear();
    GFvtx_angle.clear();
    GFvtxid.clear();
    GFvtxmom_x.clear();
    GFvtxmom_y.clear();
    GFvtxmom_z.clear();

    GFnhHtof.clear();
    GFsegHtof.clear();
    GFxHtof.clear();
    GFyHtof.clear();
    GFzHtof.clear();
    GFtracklenHtof.clear();
    GFtofHtof.clear();

    GFfitstatus_p.clear();
    GFmom_p.clear();
    GFtof_p.clear();
    GFtracklen_p.clear();

    GFinside_p.clear();
    GFxTgt_p.clear();
    GFyTgt_p.clear();
    GFzTgt_p.clear();
    GFmomTgt_p.clear();
    GFmomxTgt_p.clear();
    GFmomyTgt_p.clear();
    GFmomzTgt_p.clear();
    GFtracklenTgt_p.clear();
    GFtofTgt_p.clear();

    GFnhHtof_p.clear();
    GFsegHtof_p.clear();
    GFtofHtof_p.clear();
    GFtracklenHtof_p.clear();

    GFfitstatus_pi.clear();
    GFmom_pi.clear();
    GFtof_pi.clear();
    GFtracklen_pi.clear();

    GFinside_pi.clear();
    GFxTgt_pi.clear();
    GFyTgt_pi.clear();
    GFzTgt_pi.clear();
    GFmomTgt_pi.clear();
    GFmomxTgt_pi.clear();
    GFmomyTgt_pi.clear();
    GFmomzTgt_pi.clear();
    GFtracklenTgt_pi.clear();
    GFtofTgt_pi.clear();

    GFnhHtof_pi.clear();
    GFsegHtof_pi.clear();
    GFtofHtof_pi.clear();
    GFtracklenHtof_pi.clear();
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
  TTreeReaderValue<std::vector<Double_t>>* cdeTpc;    // cdE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;     // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;    // clock time

  Int_t    nhHtof;
  Int_t    csHtof[NumOfSegHTOF*MaxDepth];
  Double_t HtofSeg[NumOfSegHTOF*MaxDepth];
  Double_t tHtof[NumOfSegHTOF*MaxDepth];
  Double_t dtHtof[NumOfSegHTOF*MaxDepth];
  Double_t deHtof[NumOfSegHTOF*MaxDepth];
  Double_t posHtof[NumOfSegHTOF*MaxDepth];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    PadHid    = 100000,
    genfitHid = 200000
  };
}

//_____________________________________________________________________________
int
main( int argc, char **argv )
{
  std::vector<std::string> arg( argv, argv+argc );

  if( !CheckArg( arg ) )
    return EXIT_FAILURE;
  if( !DstOpen( arg ) )
    return EXIT_FAILURE;
  if( !gConf.Initialize( arg[kConfFile] ) )
    return EXIT_FAILURE;
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
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
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  delete fitter;
  return EXIT_SUCCESS;
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent( void )
{
  event.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen( std::vector<std::string> arg )
{
  int open_file = 0;
  int open_tree = 0;
  for( Int_t i=0; i<nArgc; ++i ){
    if( i==kProcess || i==kConfFile || i==kOutFile ) continue;
    open_file += OpenFile( TFileCont[i], arg[i] );
    open_tree += OpenTree( TFileCont[i], TTreeCont[i], TreeName[i] );
  }

  if( open_file!=open_tree || open_file!=nArgc-3 )
    return false;
  if( !CheckEntries( TTreeCont ) )
    return false;

  TFileCont[kOutFile] = new TFile( arg[kOutFile].c_str(), "recreate" );

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead( int ievent )
{

  Double_t vtx_scan_range = gUser.GetParameter("VertexScanRange");

  //if( ievent%1000==0 ){
  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.clkTpc = **src.clkTpc;
  event.nhHtof = src.nhHtof;
  for(Int_t it=0; it<event.nhHtof; it++){
    event.HtofSeg.push_back(src.HtofSeg[it]);
    event.tHtof.push_back(src.tHtof[it]);
    event.dtHtof.push_back(src.dtHtof[it]);
    event.deHtof.push_back(src.deHtof[it]);
    event.posHtof.push_back(src.posHtof[it]);
  }
  HF1( 1, event.status++ );

  if( **src.nhTpc < 5 ) return true;

  if(event.clkTpc.size() != 1){
    std::cerr << "something is wrong" << std::endl;
    return true;
  }

  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  //GFtracks.Init();
  HF1( 2, event.GFstatus++ );
  Double_t clock = event.clkTpc.at(0);
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);

  HF1( 1, event.status++ );

#if RawHit
  Int_t nh_Tpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = TPCAna.GetTPCHC( layer );
    for( const auto& hit : hc ){
      if( !hit || !hit->IsGood() )
        continue;
      Double_t x = hit->GetX();
      Double_t y = hit->GetY();
      Double_t z = hit->GetZ();
      Double_t de = hit->GetCDe();
      Double_t pad = hit->GetPad();
      Int_t row = hit->GetRow();
      event.raw_hitpos_x.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);
      event.raw_layer.push_back(layer);
      event.raw_row.push_back(row);
      ++nh_Tpc;
    }
  }
  event.nhTpc = nh_Tpc;
  HF1( 1, event.status++ );
#endif

  Int_t nclTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = TPCAna.GetTPCClCont( layer );
    for( const auto& cl : hc ){
      if( !cl || !cl->IsGood() )
        continue;
      Double_t x = cl->GetX();
      Double_t y = cl->GetY();
      Double_t z = cl->GetZ();
      Double_t de = cl->GetDe();
      Int_t cl_size = cl->GetClusterSize();
      Double_t mrow = cl->MeanRow();
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      TPCHit* meanHit = cl->GetMeanHit();
      Int_t houghflag = meanHit->GetHoughFlag();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      event.cluster_x.push_back(x);
      event.cluster_y.push_back(y);
      event.cluster_z.push_back(z);
      event.cluster_de.push_back(de);
      event.cluster_size.push_back(cl_size);
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(mrow);
      event.cluster_de_center.push_back(centerDe);
      event.cluster_x_center.push_back(centerPos.X());
      event.cluster_y_center.push_back(centerPos.Y());
      event.cluster_z_center.push_back(centerPos.Z());
      event.cluster_row_center.push_back(centerRow);
      event.cluster_houghflag.push_back(houghflag);
      ++nclTpc;
    }
  }
  event.nclTpc = nclTpc;
  HF1( 1, event.status++ );

  TPCAna.TrackSearchTPCHelix();

  Int_t ntTpc = TPCAna.GetNTracksTPCHelix();
  event.ntTpc = ntTpc;
  HF1( 10, ntTpc );
  if( event.ntTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  event.nhtrack.resize( ntTpc );
  event.isBeam.resize( ntTpc );
  event.isAccidental.resize( ntTpc );
  event.fittime.resize( ntTpc );
  event.searchtime.resize( ntTpc );
  event.niteration.resize( ntTpc );
  event.chisqr.resize( ntTpc );

  event.helix_cx.resize( ntTpc );
  event.helix_cy.resize( ntTpc );
  event.helix_z0.resize( ntTpc );
  event.helix_r.resize( ntTpc );
  event.helix_dz.resize( ntTpc );
  event.mom0_x.resize( ntTpc );
  event.mom0_y.resize( ntTpc );
  event.mom0_z.resize( ntTpc );
  event.mom0.resize( ntTpc );
  event.dE.resize( ntTpc );
  event.dEdx.resize( ntTpc );
  event.dz_factor.resize( ntTpc );
  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );
  event.pid.resize( ntTpc );

  event.hitlayer.resize( ntTpc );
  event.hitpos_x.resize( ntTpc );
  event.hitpos_y.resize( ntTpc );
  event.hitpos_z.resize( ntTpc );
  event.calpos_x.resize( ntTpc );
  event.calpos_y.resize( ntTpc );
  event.calpos_z.resize( ntTpc );
  event.mom_x.resize( ntTpc );
  event.mom_y.resize( ntTpc );
  event.mom_z.resize( ntTpc );
  event.residual.resize( ntTpc );
  event.residual_x.resize( ntTpc );
  event.residual_y.resize( ntTpc );
  event.residual_z.resize( ntTpc );
  event.resolution_x.resize( ntTpc );
  event.resolution_y.resize( ntTpc );
  event.resolution_z.resize( ntTpc );
  event.helix_t.resize( ntTpc );
  event.pathhit.resize( ntTpc );
  event.alpha.resize(ntTpc);
  event.houghflag.resize(ntTpc);
  event.track_cluster_de.resize(ntTpc);
  event.track_cluster_size.resize(ntTpc);
  event.track_cluster_mrow.resize(ntTpc);
  event.track_cluster_de_center.resize(ntTpc);
  event.track_cluster_x_center.resize(ntTpc);
  event.track_cluster_y_center.resize(ntTpc);
  event.track_cluster_z_center.resize(ntTpc);
  event.track_cluster_row_center.resize(ntTpc);

  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    Double_t helix_z0=tp->Getz0(), helix_r=tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
    Int_t isbeam = tp->GetIsBeam();
    Int_t isaccidental = tp->GetIsAccidental();
    double fittime = tp->GetFitTime();
    double searchtime = tp->GetSearchTime();
    Int_t charge = tp->GetCharge();
    Double_t pathlen = tp->GetPath();
    Int_t iteration = tp->GetNIteration();

    event.nhtrack[it] = nh;
    event.isBeam[it] = isbeam;
    event.isAccidental[it] = isaccidental;
    event.fittime[it] = fittime;
    event.charge[it] = charge;
    event.path[it] = pathlen;
    event.chisqr[it] = chisqr;
    event.niteration[it] = iteration;
    event.searchtime[it] = searchtime;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r ;
    event.helix_dz[it] = helix_dz;
    event.mom0_x[it] = Mom0.x();;
    event.mom0_y[it] = Mom0.y();;
    event.mom0_z[it] = Mom0.z();;
    event.mom0[it] = Mom0.Mag();;
    event.dE[it] = tp->GetTrackdE();
    event.dEdx[it] = tp->GetdEdx(truncatedMean);
    event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));

    event.hitlayer[it].resize( nh );
    event.hitpos_x[it].resize( nh );
    event.hitpos_y[it].resize( nh );
    event.hitpos_z[it].resize( nh );
    event.calpos_x[it].resize( nh );
    event.calpos_y[it].resize( nh );
    event.calpos_z[it].resize( nh );
    event.mom_x[it].resize( nh );
    event.mom_y[it].resize( nh );
    event.mom_z[it].resize( nh );
    event.residual[it].resize( nh );
    event.residual_x[it].resize( nh );
    event.residual_y[it].resize( nh );
    event.residual_z[it].resize( nh );
    event.resolution_x[it].resize( nh );
    event.resolution_y[it].resize( nh );
    event.resolution_z[it].resize( nh );
    event.helix_t[it].resize( nh );
    event.pathhit[it].resize(nh);
    event.alpha[it].resize(nh);
    event.houghflag[it].resize(nh);

    event.track_cluster_de[it].resize(nh);
    event.track_cluster_size[it].resize(nh);
    event.track_cluster_mrow[it].resize(nh);
    event.track_cluster_de_center[it].resize(nh);
    event.track_cluster_x_center[it].resize(nh);
    event.track_cluster_y_center[it].resize(nh);
    event.track_cluster_z_center[it].resize(nh);
    event.track_cluster_row_center[it].resize(nh);

    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHitInOrder( ih );
      if( !hit ) continue;
      Int_t layer = hit->GetLayer();
      Int_t houghflag = hit->GetHoughFlag();
      Double_t residual = hit->GetResidual();
      const TVector3& resi_vect = hit->GetResidualVect();
      const TVector3& res_vect = hit->GetResolutionVect();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& mom = hit->GetMomentumHelix(charge);

      Double_t clde = hit->GetDe();
      Double_t mrow = hit->GetMRow(); // same
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_mrow[it][ih] = mrow;
      event.alpha[it][ih] = tp->GetAlpha(ih);

#if TrackCluster
      TPCHit *clhit = hit -> GetHit();
      TPCCluster *cl = clhit -> GetParentCluster();
      Int_t clsize = cl->GetClusterSize();
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;
#endif
      event.pathhit[it][ih] = hit->GetPathHelix();

      event.hitlayer[it][ih] = (double)layer;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.mom_x[it][ih] = mom.x();
      event.mom_y[it][ih] = mom.y();
      event.mom_z[it][ih] = mom.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = resi_vect.x();
      event.residual_y[it][ih] = resi_vect.y();
      event.residual_z[it][ih] = resi_vect.z();
      event.resolution_x[it][ih] = res_vect.x();
      event.resolution_y[it][ih] = res_vect.y();
      event.resolution_z[it][ih] = res_vect.z();
      event.houghflag[it][ih] = houghflag;
      event.helix_t[it][ih] = hit->GetTheta();
    }
    event.pid[it]=tp->GetPid();
    std::vector<Int_t> pdgcode;
    if(event.charge[it]>0){
      pdgcode.push_back(211);
      pdgcode.push_back(2212);
    }
    else{
      pdgcode.push_back(-211);
    }
    GFtracks.AddHelixTrack(pdgcode, tp);
  }
  HF1( 2, event.GFstatus++ );

  Int_t nvtxTpc = TPCAna.GetNVerticesTPCHelix();
  event.nvtxTpc = nvtxTpc;
  event.vtx_x.resize(nvtxTpc);
  event.vtx_y.resize(nvtxTpc);
  event.vtx_z.resize(nvtxTpc);
  event.vtx_dist.resize(nvtxTpc);
  event.vtx_angle.resize(nvtxTpc);
  event.vtxid.resize(nvtxTpc);
  event.vtxmom_theta.resize(nvtxTpc);
  event.vtxpos_x.resize(nvtxTpc);
  event.vtxpos_y.resize(nvtxTpc);
  event.vtxpos_z.resize(nvtxTpc);
  event.vtxmom_x.resize(nvtxTpc);
  event.vtxmom_y.resize(nvtxTpc);
  event.vtxmom_z.resize(nvtxTpc);

  for( Int_t it=0; it<nvtxTpc; ++it ){
    TPCVertexHelix *vp = TPCAna.GetTPCVertexHelix( it );
    if( !vp ) continue;
    event.vtx_x[it] = vp -> GetVertex().x();
    event.vtx_y[it] = vp -> GetVertex().y();
    event.vtx_z[it] = vp -> GetVertex().z();
    event.vtx_dist[it] = vp -> GetClosestDist();
    event.vtx_angle[it] = vp -> GetOpeningAngle();

    event.vtxid[it].resize(2);
    event.vtxmom_theta[it].resize(2);
    event.vtxpos_x[it].resize(2);
    event.vtxpos_y[it].resize(2);
    event.vtxpos_z[it].resize(2);
    event.vtxmom_x[it].resize(2);
    event.vtxmom_y[it].resize(2);
    event.vtxmom_z[it].resize(2);

    Int_t ivtx1 = vp -> GetTrackId(0);
    Int_t ivtx2 = vp -> GetTrackId(1);
    event.vtxid[it][0] = ivtx1;
    event.vtxmom_theta[it][0] = vp -> GetTrackTheta(0);
    event.vtxpos_x[it][0] = vp -> GetTrackPos(0).x();
    event.vtxpos_y[it][0] = vp -> GetTrackPos(0).y();
    event.vtxpos_z[it][0] = vp -> GetTrackPos(0).z();
    event.vtxmom_x[it][0] = vp -> GetTrackMom(0).x();
    event.vtxmom_y[it][0] = vp -> GetTrackMom(0).y();
    event.vtxmom_z[it][0] = vp -> GetTrackMom(0).z();

    event.vtxid[it][1] = ivtx2;
    event.vtxmom_theta[it][1] = vp -> GetTrackTheta(1);
    event.vtxpos_x[it][1] = vp -> GetTrackPos(1).x();
    event.vtxpos_y[it][1] = vp -> GetTrackPos(1).y();
    event.vtxpos_z[it][1] = vp -> GetTrackPos(1).z();
    event.vtxmom_x[it][1] = vp -> GetTrackMom(1).x();
    event.vtxmom_y[it][1] = vp -> GetTrackMom(1).y();
    event.vtxmom_z[it][1] = vp -> GetTrackMom(1).z();

    if(TMath::Abs(event.vtx_dist[it]) > 3) continue;
    if(event.isBeam[ivtx1]==1) continue;
    if(event.isAccidental[ivtx1]==1) continue;
    if(event.isBeam[ivtx2]==1) continue;
    if(event.isAccidental[ivtx2]==1) continue;

    std::vector<Double_t> GFextrapolation_decays(2);
    std::vector<TVector3> GFmom_decays(2);
    Double_t GFdist;
    TVector3 GFvertex;
    Bool_t vtxcut = GFtracks.FindVertex(ivtx1, ivtx2, -1, -1,
					GFextrapolation_decays[0],
					GFextrapolation_decays[1],
					GFmom_decays[0],
					GFmom_decays[1],
					GFdist, GFvertex,
					vtx_scan_range);
    if(!vtxcut) continue;
    if(GFdist > 10.) continue;
    std::vector<Double_t> vtxid{(Double_t) ivtx1, (Double_t) ivtx2};
    std::vector<Double_t> vtxmom_x{GFmom_decays[0].x(), GFmom_decays[1].x()};
    std::vector<Double_t> vtxmom_y{GFmom_decays[0].y(), GFmom_decays[1].y()};
    std::vector<Double_t> vtxmom_z{GFmom_decays[0].z(), GFmom_decays[1].z()};

    event.GFvtx_x.push_back(GFvertex.x());
    event.GFvtx_y.push_back(GFvertex.y());
    event.GFvtx_z.push_back(GFvertex.z());
    event.GFvtx_dist.push_back(GFdist);
    event.GFvtx_angle.push_back(GFmom_decays[0].Angle(GFmom_decays[1]));
    event.GFvtxid.push_back(vtxid);
    event.GFvtxmom_x.push_back(vtxmom_x);
    event.GFvtxmom_y.push_back(vtxmom_y);
    event.GFvtxmom_z.push_back(vtxmom_z);
  }
  event.GFnvtxTpc = event.GFvtx_dist.size();

  GFtracks.FitTracks();
  int GFntTpc = GFtracks.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }
  HF1( 2, event.GFstatus++ );

  HF1( genfitHid, GFntTpc);
  event.GFntTpc = GFntTpc;
  event.GFcharge.resize(GFntTpc);
  event.GFchisqr.resize(GFntTpc);
  event.GFtof.resize(GFntTpc);
  event.GFtracklen.resize(GFntTpc);
  event.GFpval.resize(GFntTpc);
  event.GFfitstatus.resize(GFntTpc);
  event.GFpdgcode.resize(GFntTpc);
  event.GFnhtrack.resize(GFntTpc);
  event.GFlayer.resize(GFntTpc);
  event.GFpos_x.resize(GFntTpc);
  event.GFpos_y.resize(GFntTpc);
  event.GFpos_z.resize(GFntTpc);
  event.GFmom.resize(GFntTpc);
  event.GFmom_x.resize(GFntTpc);
  event.GFmom_y.resize(GFntTpc);
  event.GFmom_z.resize(GFntTpc);
  event.GFresidual_x.resize(GFntTpc);
  event.GFresidual_y.resize(GFntTpc);
  event.GFresidual_z.resize(GFntTpc);
  event.GFresidual_p.resize(GFntTpc);
  event.GFresidual_px.resize(GFntTpc);
  event.GFresidual_py.resize(GFntTpc);
  event.GFresidual_pz.resize(GFntTpc);
  event.GFinside.resize(GFntTpc);
  event.GFxTgt.resize(GFntTpc);
  event.GFyTgt.resize(GFntTpc);
  event.GFzTgt.resize(GFntTpc);
  event.GFmomTgt.resize(GFntTpc);
  event.GFmomxTgt.resize(GFntTpc);
  event.GFmomyTgt.resize(GFntTpc);
  event.GFmomzTgt.resize(GFntTpc);
  event.GFtracklenTgt.resize(GFntTpc);
  event.GFtofTgt.resize(GFntTpc);

  event.GFnhHtof.resize(GFntTpc);
  event.GFsegHtof.resize(GFntTpc);
  event.GFxHtof.resize(GFntTpc);
  event.GFyHtof.resize(GFntTpc);
  event.GFzHtof.resize(GFntTpc);
  event.GFtracklenHtof.resize(GFntTpc);
  event.GFtofHtof.resize(GFntTpc);

  event.GFfitstatus_p.resize(GFntTpc);
  event.GFmom_p.resize(GFntTpc);
  event.GFtof_p.resize(GFntTpc);
  event.GFtracklen_p.resize(GFntTpc);

  event.GFinside_p.resize(GFntTpc);
  event.GFxTgt_p.resize(GFntTpc);
  event.GFyTgt_p.resize(GFntTpc);
  event.GFzTgt_p.resize(GFntTpc);
  event.GFmomTgt_p.resize(GFntTpc);
  event.GFmomxTgt_p.resize(GFntTpc);
  event.GFmomyTgt_p.resize(GFntTpc);
  event.GFmomzTgt_p.resize(GFntTpc);
  event.GFtracklenTgt_p.resize(GFntTpc);
  event.GFtofTgt_p.resize(GFntTpc);

  event.GFnhHtof_p.resize(GFntTpc);
  event.GFtofHtof_p.resize(GFntTpc);
  event.GFtracklenHtof_p.resize(GFntTpc);
  event.GFsegHtof_p.resize(GFntTpc);

  event.GFfitstatus_pi.resize(GFntTpc);
  event.GFmom_pi.resize(GFntTpc);
  event.GFtof_pi.resize(GFntTpc);
  event.GFtracklen_pi.resize(GFntTpc);

  event.GFinside_pi.resize(GFntTpc);
  event.GFxTgt_pi.resize(GFntTpc);
  event.GFyTgt_pi.resize(GFntTpc);
  event.GFzTgt_pi.resize(GFntTpc);
  event.GFmomTgt_pi.resize(GFntTpc);
  event.GFmomxTgt_pi.resize(GFntTpc);
  event.GFmomyTgt_pi.resize(GFntTpc);
  event.GFmomzTgt_pi.resize(GFntTpc);
  event.GFtracklenTgt_pi.resize(GFntTpc);
  event.GFtofTgt_pi.resize(GFntTpc);

  event.GFnhHtof_pi.resize(GFntTpc);
  event.GFtofHtof_pi.resize(GFntTpc);
  event.GFtracklenHtof_pi.resize(GFntTpc);
  event.GFsegHtof_pi.resize(GFntTpc);

  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    event.GFfitstatus[igf] = (int)GFtracks.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
    if(!GFtracks.TrackCheck(igf)) continue;
    int nh = GFtracks.GetNHits(igf);
    event.GFlayer[igf].resize(nh);
    event.GFpos_x[igf].resize(nh);
    event.GFpos_y[igf].resize(nh);
    event.GFpos_z[igf].resize(nh);
    event.GFmom[igf].resize(nh);
    event.GFmom_x[igf].resize(nh);
    event.GFmom_y[igf].resize(nh);
    event.GFmom_z[igf].resize(nh);
    event.GFresidual_x[igf].resize(nh);
    event.GFresidual_y[igf].resize(nh);
    event.GFresidual_z[igf].resize(nh);
    event.GFresidual_p[igf].resize(nh);
    event.GFresidual_px[igf].resize(nh);
    event.GFresidual_py[igf].resize(nh);
    event.GFresidual_pz[igf].resize(nh);

    event.GFchisqr[igf] = GFtracks.GetChi2NDF(igf);
    event.GFcharge[igf] = GFtracks.GetCharge(igf);
    event.GFtof[igf] = GFtracks.GetTrackTOF(igf, 0, -1);
    event.GFtracklen[igf] = GFtracks.GetTrackLength(igf, 0, -1);
    event.GFpval[igf] = GFtracks.GetPvalue(igf);
    event.GFnhtrack[igf] = GFtracks.GetNHits(igf);
    event.GFpdgcode[igf] = GFtracks.GetPDGcode(igf);

    HF1( genfitHid+1, event.GFchisqr[igf]);
    HF1( genfitHid+2, event.GFpval[igf]);
    HF1( genfitHid+3, event.GFcharge[igf]);
    HF1( genfitHid+4, event.GFnhtrack[igf]);
    HF1( genfitHid+5, event.GFtracklen[igf]);
    HF1( genfitHid+6, event.GFtof[igf]);
    for( Int_t ihit=0; ihit<nh; ++ihit ){
      TVector3 hit = GFtracks.GetPos(igf, ihit);
      TVector3 mom = GFtracks.GetMom(igf, ihit);
      Int_t layer = (int)event.hitlayer[igf][ihit];
      event.GFlayer[igf][ihit] = layer;
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

      double chargetest = event.GFcharge[igf]*event.charge[igf];
      event.GFresidual_p[igf][ihit] = mom.Mag() - event.mom0[igf];
      event.GFresidual_px[igf][ihit] = mom.x() - chargetest*event.mom_x[igf][ihit];
      event.GFresidual_py[igf][ihit] = mom.y() - chargetest*event.mom_y[igf][ihit];
      event.GFresidual_pz[igf][ihit] = mom.z() - chargetest*event.mom_z[igf][ihit];
      if(ihit==0) HF1( genfitHid+7, event.GFmom[igf][0]);
      HF1( genfitHid+8, event.GFlayer[igf][ihit]);
      HF1( genfitHid+10, event.GFresidual_x[igf][ihit]);
      HF1( genfitHid+11, event.GFresidual_y[igf][ihit]);
      HF1( genfitHid+12, event.GFresidual_z[igf][ihit]);
      HF1( genfitHid+13, event.GFresidual_p[igf][ihit]);
      HF1( genfitHid+14, event.GFresidual_px[igf][ihit]);
      HF1( genfitHid+15, event.GFresidual_py[igf][ihit]);
      HF1( genfitHid+16, event.GFresidual_pz[igf][ihit]);
      HF1( genfitHid+1000*(layer+1), event.GFresidual_x[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+1, event.GFresidual_y[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+2, event.GFresidual_z[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+3, event.GFresidual_p[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+4, event.GFresidual_px[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+5, event.GFresidual_py[igf][ihit]);
      HF1( genfitHid+1000*(layer+1)+6, event.GFresidual_pz[igf][ihit]);
    } //ihit

    //Extrapolation
    if(GFtracks.IsInsideTarget(igf)){
      event.GFinside[igf] = 1;
      TVector3 posv; TVector3 momv; double len; double tof;
      GFtracks.ExtrapolateToTarget(igf, posv, momv, len, tof);
      HF1( genfitHid+20, posv.x());
      HF1( genfitHid+21, posv.y());
      HF1( genfitHid+22, posv.z());
      event.GFxTgt[igf] = posv.x();
      event.GFyTgt[igf] = posv.y();
      event.GFzTgt[igf] = posv.z();
      event.GFmomTgt[igf] = momv.Mag();
      event.GFmomxTgt[igf] = momv.x();
      event.GFmomyTgt[igf] = momv.y();
      event.GFmomzTgt[igf] = momv.z();
      event.GFtracklenTgt[igf] = len;
      event.GFtofTgt[igf] = tof;
    }
    else event.GFinside[igf] = 0;

    int candidates = 0; int htofid[8]; TVector3 htofpos[8]; TVector3 htofmom[8]; double htoflen[8]; double htoftof[8];
    if(GFtracks.ExtrapolateToHTOF(igf, candidates, htofid, htofpos, htofmom, htoflen, htoftof)){
      event.GFnhHtof[igf]=candidates;
      event.GFsegHtof[igf].resize(candidates);
      event.GFxHtof[igf].resize(candidates);
      event.GFyHtof[igf].resize(candidates);
      event.GFzHtof[igf].resize(candidates);
      event.GFtracklenHtof[igf].resize(candidates);
      event.GFtofHtof[igf].resize(candidates);
      for( Int_t ihit=0; ihit<candidates; ++ihit ){
	event.GFsegHtof[igf][ihit]=htofid[ihit];
	event.GFxHtof[igf][ihit]=htofpos[ihit].x();
	event.GFyHtof[igf][ihit]=htofpos[ihit].y();
	event.GFzHtof[igf][ihit]=htofpos[ihit].z();
	event.GFtracklenHtof[igf][ihit]=htoflen[ihit];
	event.GFtofHtof[igf][ihit]=htoftof[ihit];

	HF1( genfitHid+24, htoflen[ihit]);
	HF1( genfitHid+25, htoftof[ihit]);
	HF1( genfitHid+26, htofpos[ihit].x());
	HF1( genfitHid+27, htofpos[ihit].y());
	HF1( genfitHid+28, htofpos[ihit].z());
	HF1( genfitHid+29, htofid[ihit]);
      }
    }
    HF1( genfitHid+23, event.GFnhHtof[igf]);
  } //igf



  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    Int_t repid = 0; //pion
    TVector3 posv; TVector3 momv; double len; double tof;
    if(GFtracks.TrackCheck(igf, repid) && GFtracks.IsInsideTarget(igf, repid) &&
       GFtracks.ExtrapolateToTarget(igf, posv, momv, len, tof, repid)){

      event.GFxTgt_pi[igf] = posv.x();
      event.GFyTgt_pi[igf] = posv.y();
      event.GFzTgt_pi[igf] = posv.z();
      event.GFmomTgt_pi[igf] = momv.Mag();
      event.GFmomxTgt_pi[igf] = momv.x();
      event.GFmomyTgt_pi[igf] = momv.y();
      event.GFmomzTgt_pi[igf] = momv.z();
      event.GFtracklenTgt_pi[igf] = len;
      event.GFtofTgt_pi[igf] = tof;

      event.GFtof_pi[igf] = GFtracks.GetTrackTOF(igf, 0, -1, repid);
      event.GFmom_pi[igf] = GFtracks.GetMom(igf, 0, repid).Mag();
      event.GFtracklen_pi[igf] = GFtracks.GetTrackLength(igf, 0, -1, repid);

      int candidates = 0; int htofid[8]; TVector3 htofpos[8]; TVector3 htofmom[8]; double htoflen[8]; double htoftof[8];
      if(GFtracks.ExtrapolateToHTOF(igf, candidates, htofid, htofpos, htofmom, htoflen, htoftof, repid)){
	event.GFnhHtof_pi[igf]=candidates;
	event.GFsegHtof_pi[igf].resize(candidates);
	event.GFtracklenHtof_pi[igf].resize(candidates);
	event.GFtofHtof_pi[igf].resize(candidates);
	for( Int_t ihit=0; ihit<candidates; ++ihit ){
	  event.GFsegHtof_pi[igf][ihit]=htofid[ihit];
	  event.GFtracklenHtof_pi[igf][ihit]=htoflen[ihit];
	  event.GFtofHtof_pi[igf][ihit]=htoftof[ihit];
	}
      }
    } //pion

    if(event.charge[igf]==1){ //proton
      repid = 1;
      if(GFtracks.TrackCheck(igf, repid) && GFtracks.IsInsideTarget(igf, repid) &&
	 GFtracks.ExtrapolateToTarget(igf, posv, momv, len, tof, repid)){

	event.GFxTgt_p[igf] = posv.x();
	event.GFyTgt_p[igf] = posv.y();
	event.GFzTgt_p[igf] = posv.z();
	event.GFmomTgt_p[igf] = momv.Mag();
	event.GFmomxTgt_p[igf] = momv.x();
	event.GFmomyTgt_p[igf] = momv.y();
	event.GFmomzTgt_p[igf] = momv.z();
	event.GFtracklenTgt_p[igf] = len;
	event.GFtofTgt_p[igf] = tof;

	event.GFtof_p[igf] = GFtracks.GetTrackTOF(igf, 0, -1, repid);
	event.GFmom_p[igf] = GFtracks.GetMom(igf, 0, repid).Mag();
	event.GFtracklen_p[igf] = GFtracks.GetTrackLength(igf, 0, -1, repid);

	int candidates = 0; int htofid[8]; TVector3 htofpos[8]; TVector3 htofmom[8]; double htoflen[8]; double htoftof[8];
	if(GFtracks.ExtrapolateToHTOF(igf, candidates, htofid, htofpos, htofmom, htoflen, htoftof, repid)){
	  event.GFnhHtof_p[igf]=candidates;
	  event.GFsegHtof_p[igf].resize(candidates);
	  event.GFtracklenHtof_p[igf].resize(candidates);
	  event.GFtofHtof_p[igf].resize(candidates);
	  for( Int_t ihit=0; ihit<candidates; ++ihit ){
	    event.GFsegHtof_p[igf][ihit]=htofid[ihit];
	    event.GFtracklenHtof_p[igf][ihit]=htoflen[ihit];
	    event.GFtofHtof_p[igf][ihit]=htoftof[ihit];
	  }
	}
      }
    } //proton
  } //igf

  HF1( 2, event.GFstatus++);
  HF1( 1, event.status++ );

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for( Int_t i=0; i<n; ++i ){
    if( TTreeReaderCont[i] ) delete TTreeReaderCont[i];
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{
  /*
    HB1( 1, "Status", 21, 0., 21. );
    HB1( 2, "Genfit Status", 20, 0., 20. );
    HB1( 3, "Genfit Fit Status", 2, 0., 2. );
    HB1( 10, "NTrack TPC", 20, 0., 20. );

    HB1(genfitHid, "[GenFit] #Track TPC; #Track; Counts", 20, 0., 20. );
    HB1(genfitHid+1, "[GenFit] Chisqr/ndf; ; Counts", 200, 0, 10 );
    HB1(genfitHid+2, "[GenFit] p-value; p-value; Counts", 100, -0.05, 1.05);
    HB1(genfitHid+3, "[GenFit] Charge;", 6, -3, 3 );
    HB1(genfitHid+4, "[GenFit] #Hits of Track", 50, 0., 50. );
    HB1(genfitHid+5, "[GenFit] Track Length; Length [mm]; Counts [/1 mm]", 500, 0, 500 );
    HB1(genfitHid+6, "[GenFit] Tof; Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
    HB1(genfitHid+7, "[GenFit] Reconstructed P; P [GeV/c]; Counts [/0.001 GeV/c]", 1500, 0., 1.5 );
    HB1(genfitHid+8, "[GenFit] LayerID", 33, 0., 33. );

    //Residuals
    HB1(genfitHid+10, "[GenFit] Residual x; Residual x [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+11, "[GenFit] Residual y; Residual y [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+12, "[GenFit] Residual z; Residual z [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+13, "[GenFit] Residual P; Residual P [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
    HB1(genfitHid+14, "[GenFit] Residual Px; Residual Px [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(genfitHid+15, "[GenFit] Residual Py; Residual Py [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(genfitHid+16, "[GenFit] Residual Pz; Residual Pz [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    HB1(genfitHid+(layer+1)*1000, Form("[GenFit] Residual x Layer%d; Residual x [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+1, Form("[GenFit] Residual y Layer%d; Residual y [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+2, Form("[GenFit] Residual z Layer%d; Residual z [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+3, Form("[GenFit] Residual P Layer%d; Residual P [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+4, Form("[GenFit] Residual Px Layer%d; Residual Px [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+5, Form("[GenFit] Residual Py Layer%d; Residual Py [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+6, Form("[GenFit] Residual Pz Layer%d; Residual Pz [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
  }

  //Extrapolation
  HB1(genfitHid+20, "[GenFit] Vertex X (extrapolated to the Target); X [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+21, "[GenFit] Vertex Y (extrapolated to the Target); Vertex Y [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+22, "[GenFit] Vertex Z (extrapolated to the Target); Vertex Z [mm]; Counts [/0.1 mm]", 500, -143 -25, -143 + 25 );
  HB1(genfitHid+23, "[GenFit] #Hits (extrapolated to the HTOF); #Hits; Counts", 10, 0, 10 );
  HB1(genfitHid+24, "[GenFit] Track Length (extrapolated to the HTOF); Length [mm]; Counts [/1 mm]", 1000, -500, 500 );
  HB1(genfitHid+25, "[GenFit] Tof (extrapolated to the HTOF); Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
  HB1(genfitHid+26, "[GenFit] X (extrapolated to the HTOF); X [mm]; Counts [/0.1 mm]", 10000, -500, 500 );
  HB1(genfitHid+27, "[GenFit] Y (extrapolated to the HTOF); Y [mm]; Counts [/0.1 mm]", 8000, -400, 400 );
  HB1(genfitHid+28, "[GenFit] Z (extrapolated to the HTOF); Z [mm]; Counts [/0.1 mm]", 10000, -500, 500 );
  HB1(genfitHid+29, "[GenFit] HTOF ID; #ID ; Counts", 36, 0, 36 );
*/
  HBTree( "tpc", "tree of DstTPCTracking" );

  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
  tree->Branch( "clkTpc", &event.clkTpc);

  tree->Branch( "nhHtof", &event.nhHtof );
  tree->Branch( "HtofSeg", &event.HtofSeg );
  tree->Branch( "tHtof", &event.tHtof );
  tree->Branch( "dtHtof", &event.dtHtof );
  tree->Branch( "deHtof", &event.deHtof );
  tree->Branch( "posHtof", &event.posHtof );
#if RawHit
  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
#endif

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "cluster_x", &event.cluster_x );
  tree->Branch( "cluster_y", &event.cluster_y );
  tree->Branch( "cluster_z", &event.cluster_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_size", &event.cluster_size );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_row_center", &event.cluster_row_center );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_de_center", &event.cluster_de_center );
  tree->Branch( "cluster_x_center", &event.cluster_x_center );
  tree->Branch( "cluster_y_center", &event.cluster_y_center );
  tree->Branch( "cluster_z_center", &event.cluster_z_center );
  tree->Branch( "cluster_houghflag", &event.cluster_houghflag );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "fittime", &event.fittime );
  tree->Branch( "searchtime", &event.searchtime );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "niteration", &event.niteration );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "mom0_x", &event.mom0_x );
  tree->Branch( "mom0_y", &event.mom0_y );
  tree->Branch( "mom0_z", &event.mom0_z );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );

  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );
  tree->Branch( "pid", &event.pid );

  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  tree->Branch( "mom_x", &event.mom_x );
  tree->Branch( "mom_y", &event.mom_y );
  tree->Branch( "mom_z", &event.mom_z );
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "resolution_x", &event.resolution_x );
  tree->Branch( "resolution_y", &event.resolution_y );
  tree->Branch( "resolution_z", &event.resolution_z );
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "houghflag", &event.houghflag );
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
#if TrackCluster
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
#endif

  tree->Branch( "nvtxTpc", &event.nvtxTpc );
  tree->Branch( "vtx_x", &event.vtx_x );
  tree->Branch( "vtx_y", &event.vtx_y );
  tree->Branch( "vtx_z", &event.vtx_z );
  tree->Branch( "vtx_dist", &event.vtx_dist );
  tree->Branch( "vtx_angle", &event.vtx_angle );
  tree->Branch( "vtxid", &event.vtxid );
  tree->Branch( "vtxmom_theta", &event.vtxmom_theta );
  tree->Branch( "vtxpos_x", &event.vtxpos_x );
  tree->Branch( "vtxpos_y", &event.vtxpos_y );
  tree->Branch( "vtxpos_z", &event.vtxpos_z );
  tree->Branch( "vtxmom_x", &event.vtxmom_x );
  tree->Branch( "vtxmom_y", &event.vtxmom_y );
  tree->Branch( "vtxmom_z", &event.vtxmom_z );

  tree->Branch("GFnhHtof", &event.GFnhHtof);
  tree->Branch("GFsegHtof", &event.GFsegHtof);
  tree->Branch("GFxHtof", &event.GFxHtof);
  tree->Branch("GFyHtof", &event.GFyHtof);
  tree->Branch("GFzHtof", &event.GFzHtof);
  tree->Branch("GFtracklenHtof", &event.GFtracklenHtof);
  tree->Branch("GFtofHtof", &event.GFtofHtof);

  //track fitting results
  tree->Branch("GFstatus", &event.GFstatus);
  tree->Branch("GFntTpc", &event.GFntTpc);
  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFtof", &event.GFtof);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFfitstatus", &event.GFfitstatus);
  tree->Branch("GFpdgcode", &event.GFpdgcode);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFlayer", &event.GFlayer);
  tree->Branch("GFpos_x", &event.GFpos_x);
  tree->Branch("GFpos_y", &event.GFpos_y);
  tree->Branch("GFpos_z", &event.GFpos_z);
  tree->Branch("GFmom", &event.GFmom);
  tree->Branch("GFmom_x", &event.GFmom_x);
  tree->Branch("GFmom_y", &event.GFmom_y);
  tree->Branch("GFmom_z", &event.GFmom_z);
  tree->Branch("GFresidual_x", &event.GFresidual_x);
  tree->Branch("GFresidual_y", &event.GFresidual_y);
  tree->Branch("GFresidual_z", &event.GFresidual_z);
  tree->Branch("GFresidual_p", &event.GFresidual_p);
  tree->Branch("GFresidual_px", &event.GFresidual_px);
  tree->Branch("GFresidual_py", &event.GFresidual_py);
  tree->Branch("GFresidual_pz", &event.GFresidual_pz);

  //extrapolation
  tree->Branch("GFinside", &event.GFinside);
  tree->Branch("GFxTgt", &event.GFxTgt);
  tree->Branch("GFyTgt", &event.GFyTgt);
  tree->Branch("GFzTgt", &event.GFzTgt);
  tree->Branch("GFmomTgt", &event.GFmomTgt);
  tree->Branch("GFmomxTgt", &event.GFmomxTgt);
  tree->Branch("GFmomyTgt", &event.GFmomyTgt);
  tree->Branch("GFmomzTgt", &event.GFmomzTgt);
  tree->Branch("GFtracklenTgt", &event.GFtracklenTgt);
  tree->Branch("GFtofTgt", &event.GFtofTgt);

  tree->Branch("GFnvtxTpc", &event.GFnvtxTpc);
  tree->Branch("GFvtx_x", &event.GFvtx_x);
  tree->Branch("GFvtx_y", &event.GFvtx_y);
  tree->Branch("GFvtx_z", &event.GFvtx_z);
  tree->Branch("GFvtx_dist", &event.GFvtx_dist);
  tree->Branch("GFvtx_angle", &event.GFvtx_angle);
  tree->Branch("GFvtxid", &event.GFvtxid);
  tree->Branch("GFvtxmom_x", &event.GFvtxmom_x);
  tree->Branch("GFvtxmom_y", &event.GFvtxmom_y);
  tree->Branch("GFvtxmom_z", &event.GFvtxmom_z);

  tree->Branch("GFfitstatus_p", &event.GFfitstatus_p);
  tree->Branch("GFtof_p", &event.GFtof_p);
  tree->Branch("GFmom_p", &event.GFmom_p);
  tree->Branch("GFtracklen_p", &event.GFtracklen_p);

  tree->Branch("GFinside_p", &event.GFinside_p);
  tree->Branch("GFxTgt_p", &event.GFxTgt_p);
  tree->Branch("GFyTgt_p", &event.GFyTgt_p);
  tree->Branch("GFzTgt_p", &event.GFzTgt_p);
  tree->Branch("GFmomTgt_p", &event.GFmomTgt_p);
  tree->Branch("GFmomxTgt_p", &event.GFmomxTgt_p);
  tree->Branch("GFmomyTgt_p", &event.GFmomyTgt_p);
  tree->Branch("GFmomzTgt_p", &event.GFmomzTgt_p);
  tree->Branch("GFtracklenTgt_p", &event.GFtracklenTgt_p);
  tree->Branch("GFtofTgt_p", &event.GFtofTgt_p);

  tree->Branch("GFnhHtof_p", &event.GFnhHtof_p);
  tree->Branch("GFsegHtof_p", &event.GFsegHtof_p);
  tree->Branch("GFtofHtof_p", &event.GFtofHtof_p);
  tree->Branch("GFtracklenHtof_p", &event.GFtracklenHtof_p);

  tree->Branch("GFfitstatus_pi", &event.GFfitstatus_pi);
  tree->Branch("GFtof_pi", &event.GFtof_pi);
  tree->Branch("GFmom_pi", &event.GFmom_pi);
  tree->Branch("GFtracklen_pi", &event.GFtracklen_pi);

  tree->Branch("GFinside_pi", &event.GFinside_pi);
  tree->Branch("GFxTgt_pi", &event.GFxTgt_pi);
  tree->Branch("GFyTgt_pi", &event.GFyTgt_pi);
  tree->Branch("GFzTgt_pi", &event.GFzTgt_pi);
  tree->Branch("GFmomTgt_pi", &event.GFmomTgt_pi);
  tree->Branch("GFmomxTgt_pi", &event.GFmomxTgt_pi);
  tree->Branch("GFmomyTgt_pi", &event.GFmomyTgt_pi);
  tree->Branch("GFmomzTgt_pi", &event.GFmomzTgt_pi);
  tree->Branch("GFtracklenTgt_pi", &event.GFtracklenTgt_pi);
  tree->Branch("GFtofTgt_pi", &event.GFtofTgt_pi);

  tree->Branch("GFnhHtof_pi", &event.GFnhHtof_pi);
  tree->Branch("GFsegHtof_pi", &event.GFsegHtof_pi);
  tree->Branch("GFtofHtof_pi", &event.GFtofHtof_pi);
  tree->Branch("GFtracklenHtof_pi", &event.GFtracklenHtof_pi);

  TTreeReaderCont[kTpcHit] = new TTreeReader( "tpc", TFileCont[kTpcHit] );
  const auto& reader = TTreeReaderCont[kTpcHit];
  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );
  src.npadTpc = new TTreeReaderValue<Int_t>( *reader, "npadTpc" );
  src.nhTpc = new TTreeReaderValue<Int_t>( *reader, "nhTpc" );
  src.layerTpc = new TTreeReaderValue<std::vector<Int_t>>( *reader, "layerTpc" );
  src.rowTpc = new TTreeReaderValue<std::vector<Int_t>>( *reader, "rowTpc" );
  src.padTpc = new TTreeReaderValue<std::vector<Int_t>>( *reader, "padTpc" );
  src.pedTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pedTpc" );
  src.rmsTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "rmsTpc" );
  src.deTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "deTpc" );
  src.cdeTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cdeTpc" );
  src.tTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tTpc" );
  src.ctTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ctTpc" );
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTpc" );
  src.clkTpc = new TTreeReaderValue<std::vector<Double_t>>(*reader, "clkTpc");

  ////////// Bring Address From Dst
  TTreeCont[kHodoscope]->SetBranchStatus("*", 0);

  TTreeCont[kHodoscope]->SetBranchStatus("nhHtof",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("csHtof",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("HtofSeg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("tHtof",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtHtof",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("deHtof",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("posHtof",   1);

  TTreeCont[kHodoscope]->SetBranchAddress("nhHtof",   &src.nhHtof);
  TTreeCont[kHodoscope]->SetBranchAddress("HtofSeg",   src.HtofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tHtof",     src.tHtof);
  TTreeCont[kHodoscope]->SetBranchAddress("dtHtof",    src.dtHtof);
  TTreeCont[kHodoscope]->SetBranchAddress("deHtof",    src.deHtof);
  TTreeCont[kHodoscope]->SetBranchAddress("posHtof",    src.posHtof);
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
