// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <Math/ProbFunc.h>

#include <TGeoPhysicalConstants.h>

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
#include "DCDriftParamMan.hh"
#include "DCParameters.hh"
#include "DCTdcCalibMan.hh"
#include "DCLTrackHit.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "NuclearMass.hh"
#include "RootHelper.hh"
#include "LorentzVector.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRKTrack.hh"
#include "UserParamMan.hh"

#define SaveHistograms 0
#define SaveRawHit 0
#define SaveCluster 1
#define TrackClusterHist 0
#define TruncatedMean 0
#define TrackSearchFailed 0

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gGeom  = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
//const auto& gPHC  = HodoPHCMan::GetInstance();
const Double_t truncatedMean = 0.8; //E42 reference 80%
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[TPCHit]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc","" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{
  UInt_t run_number;
  UInt_t event_number;
  std::vector<std::vector<Double_t>> trig_flag;
  std::vector<Double_t> trig_pat;
  std::vector<Double_t> clkTpc;
  std::vector<Double_t> cobo_id;

  Int_t nhTpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;

  Int_t nclTpc;
  Int_t remain_nclTpc;
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

  Int_t ntTpc; // Number of tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> nhtrackEff; // Number of Hits actually used in tracking.
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> isMultiloop;
  std::vector<Int_t> flag;
  std::vector<Int_t> fittime;  //usec
  std::vector<Int_t> searchtime; //usec
  std::vector<Int_t> niteration; //usec
  std::vector<Double_t> chisqr;
  std::vector<Double_t> pval;
  std::vector<Double_t> distTgt;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx

  std::vector<Double_t> dEdx_0;
  std::vector<Double_t> dEdx_10;
  std::vector<Double_t> dEdx_20;
  std::vector<Double_t> dEdx_30;
  std::vector<Double_t> dEdx_40;
  std::vector<Double_t> dEdx_50;
  std::vector<Double_t> dEdx_60;

  std::vector<Double_t> dz_factor;
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
  std::vector<std::vector<Double_t>> residual_t;
  std::vector<std::vector<Double_t>> residual_x;
  std::vector<std::vector<Double_t>> residual_y;
  std::vector<std::vector<Double_t>> residual_z;
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> pull;
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
  //Exclusive tracking
  std::vector<std::vector<Double_t>> exresidual_t;
  std::vector<std::vector<Double_t>> exresidual_x;
  std::vector<std::vector<Double_t>> exresidual_y;
  std::vector<std::vector<Double_t>> exresidual_z;
  //Geometric mean of inclusive & exclusive residual
  std::vector<std::vector<Double_t>> intrinsic_residual_t;
  std::vector<std::vector<Double_t>> intrinsic_residual_x;
  std::vector<std::vector<Double_t>> intrinsic_residual_y;
  std::vector<std::vector<Double_t>> intrinsic_residual_z;

  //Inverted charge track
  std::vector<Int_t> chargeIndistinguishable;
  std::vector<Double_t> chisqr_inverted;
  std::vector<Double_t> pval_inverted;
  std::vector<Double_t> helix_cx_inverted;
  std::vector<Double_t> helix_cy_inverted;
  std::vector<Double_t> helix_z0_inverted;
  std::vector<Double_t> helix_r_inverted;
  std::vector<Double_t> helix_dz_inverted;
  std::vector<Double_t> mom0_x_inverted;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y_inverted;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z_inverted;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_inverted;//Helix momentum at Y = 0
  std::vector<Int_t> pid_inverted;

  Int_t failed_ntTpc; // Number of Tracks
  std::vector<Int_t> failed_nhtrack;
  std::vector<Int_t> failed_isBeam; // isBeam: 1 = Beam, 0 = Scat
  std::vector<Int_t> failed_nclbeforetgt;
  std::vector<Int_t> failed_isAccidental;
  std::vector<Int_t> failed_flag;
  std::vector<Int_t> failed_fittime; //usec
  std::vector<Int_t> failed_searchtime; //usec
  std::vector<Int_t> failed_niteration; //usec
  std::vector<Double_t> failed_helix_cx;
  std::vector<Double_t> failed_helix_cy;
  std::vector<Double_t> failed_helix_z0;
  std::vector<Double_t> failed_helix_r;
  std::vector<Double_t> failed_helix_dz;
  std::vector<Double_t> failed_mom0;//Helix momentum at Y = 0
  std::vector<Int_t> failed_charge;//Helix charge

  std::vector<std::vector<Double_t>> failed_hitlayer;
  std::vector<std::vector<Double_t>> failed_hitpos_x;
  std::vector<std::vector<Double_t>> failed_hitpos_y;
  std::vector<std::vector<Double_t>> failed_hitpos_z;
  std::vector<std::vector<Double_t>> failed_calpos_x;
  std::vector<std::vector<Double_t>> failed_calpos_y;
  std::vector<std::vector<Double_t>> failed_calpos_z;
  std::vector<std::vector<Double_t>> failed_helix_t;
  std::vector<std::vector<Double_t>> failed_residual;
  std::vector<std::vector<Double_t>> failed_residual_x;
  std::vector<std::vector<Double_t>> failed_residual_y;
  std::vector<std::vector<Double_t>> failed_residual_z;
  std::vector<std::vector<Double_t>> failed_track_cluster_de;
  std::vector<std::vector<Double_t>> failed_track_cluster_size;
  std::vector<std::vector<Double_t>> failed_track_cluster_mrow;

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

  std::vector<Int_t> isLambda;
  std::vector<Int_t> ncombiLambda;
  std::vector<Double_t> distLambda;
  std::vector<Double_t> angleLambda;
  std::vector<Double_t> bestmassLambda;
  std::vector<std::vector<Double_t>> massLambda;
  std::vector<std::vector<Double_t>> vtxLambda_x;
  std::vector<std::vector<Double_t>> vtxLambda_y;
  std::vector<std::vector<Double_t>> vtxLambda_z;
  std::vector<std::vector<Double_t>> momLambda;
  std::vector<std::vector<Double_t>> momLambda_x;
  std::vector<std::vector<Double_t>> momLambda_y;
  std::vector<std::vector<Double_t>> momLambda_z;
  std::vector<std::vector<Double_t>> decaysidLambda;
  std::vector<std::vector<Double_t>> decaysmomLambda;
  std::vector<std::vector<Double_t>> decaysmomLambda_x;
  std::vector<std::vector<Double_t>> decaysmomLambda_y;
  std::vector<std::vector<Double_t>> decaysmomLambda_z;

  Int_t nvtxTpcClustered;
  std::vector<Double_t> Clusteredvtx_x;
  std::vector<Double_t> Clusteredvtx_y;
  std::vector<Double_t> Clusteredvtx_z;
  std::vector<std::vector<Double_t>> Clusteredvtxid;

  void clear( void )
  {
    run_number = 0;
    event_number = 0;
    trig_pat.clear();
    trig_flag.clear();
    clkTpc.clear();
    cobo_id.clear();

    nhTpc = 0;
    raw_hitpos_x.clear();
    raw_hitpos_y.clear();
    raw_hitpos_z.clear();
    raw_de.clear();
    raw_padid.clear();
    raw_layer.clear();
    raw_row.clear();

    nclTpc = 0;
    remain_nclTpc = 0;
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
    nhtrackEff.clear();
    isBeam.clear();
    isAccidental.clear();
    isMultiloop.clear();
    flag.clear();
    fittime.clear();
    searchtime.clear();
    niteration.clear();
    chisqr.clear();
    pval.clear();
    distTgt.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();
#if TruncatedMean
    dEdx_0.clear();
    dEdx_10.clear();
    dEdx_20.clear();
    dEdx_30.clear();
    dEdx_40.clear();
    dEdx_50.clear();
    dEdx_60.clear();
#endif
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
    residual_t.clear();
    residual_x.clear();
    residual_y.clear();
    residual_z.clear();
    resolution_x.clear();
    resolution_y.clear();
    resolution_z.clear();
    pull.clear();
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

    exresidual_t.clear();
    exresidual_x.clear();
    exresidual_y.clear();
    exresidual_z.clear();
    intrinsic_residual_t.clear();
    intrinsic_residual_x.clear();
    intrinsic_residual_y.clear();
    intrinsic_residual_z.clear();

    chargeIndistinguishable.clear();
    chisqr_inverted.clear();
    pval_inverted.clear();
    helix_cx_inverted.clear();
    helix_cy_inverted.clear();
    helix_z0_inverted.clear();
    helix_r_inverted.clear();
    helix_dz_inverted.clear();
    mom0_x_inverted.clear();
    mom0_y_inverted.clear();
    mom0_z_inverted.clear();
    mom0_inverted.clear();
    pid_inverted.clear();

    failed_ntTpc = 0;
    failed_nhtrack.clear();
    failed_isBeam.clear();
    failed_nclbeforetgt.clear();
    failed_isAccidental.clear();
    failed_flag.clear();
    failed_fittime.clear();
    failed_searchtime.clear();
    failed_niteration.clear();

    failed_helix_cx.clear();
    failed_helix_cy.clear();
    failed_helix_z0.clear();
    failed_helix_r.clear();
    failed_helix_dz.clear();
    failed_mom0.clear();
    failed_charge.clear();

    failed_hitlayer.clear();
    failed_hitpos_x.clear();
    failed_hitpos_y.clear();
    failed_hitpos_z.clear();
    failed_calpos_x.clear();
    failed_calpos_y.clear();
    failed_calpos_z.clear();
    failed_helix_t.clear();
    failed_residual.clear();
    failed_residual_x.clear();
    failed_residual_y.clear();
    failed_residual_z.clear();
    failed_track_cluster_de.clear();
    failed_track_cluster_size.clear();
    failed_track_cluster_mrow.clear();

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

    nvtxTpcClustered = 0;
    Clusteredvtx_x.clear();
    Clusteredvtx_y.clear();
    Clusteredvtx_z.clear();
    Clusteredvtxid.clear();

    isLambda.clear();
    ncombiLambda.clear();
    distLambda.clear();
    angleLambda.clear();
    bestmassLambda.clear();
    massLambda.clear();
    vtxLambda_x.clear();
    vtxLambda_y.clear();
    vtxLambda_z.clear();
    momLambda.clear();
    momLambda_x.clear();
    momLambda_y.clear();
    momLambda_z.clear();
    decaysidLambda.clear();
    decaysmomLambda.clear();
    decaysmomLambda_x.clear();
    decaysmomLambda_y.clear();
    decaysmomLambda_z.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<UInt_t>* run_number;
  TTreeReaderValue<UInt_t>* event_number;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* trig_flag;
  TTreeReaderValue<std::vector<Double_t>>* trig_pat;
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
  TTreeReaderValue<std::vector<Double_t>>* cobo_id;   // cobo id
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    TPCClHid = 500000,
    TPCInclusiveHid = 600000,
    TPCExclusiveHid = 700000,
    TPCIntrinsicHid = 800000
  };

Double_t
TranseverseDistance(Double_t x_center, Double_t z_center, Double_t x, Double_t z)
{
  Double_t dummy = TMath::Sqrt((x-x_center)*(x-x_center) + (z-z_center)*(z-z_center));
  Double_t dist;
  if(x_center-x<0) dist=-1.*dummy;
  else dist=dummy;
  return dist;
}

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
  if( !gConf.InitializeHistograms() )
    return EXIT_FAILURE;
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if(max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();
  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

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
  auto start = std::chrono::high_resolution_clock::now();
  static const auto KaonMass    = pdg::KaonMass();
  static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();
  static const Bool_t ExclusiveTracking = gUser.GetParameter("ExclusiveTracking");
  GetEntry(ievent);
  event.run_number = **src.run_number;
  event.event_number = **src.event_number;
  event.trig_pat = **src.trig_pat;
  event.trig_flag = **src.trig_flag;
  event.clkTpc = **src.clkTpc;
  event.cobo_id = **src.cobo_id;

  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, **src.clkTpc);
  TPCAna.TrackSearchTPCHelix(ExclusiveTracking);
  Int_t ntTpc = TPCAna.GetNTracksTPCHelix();
  event.ntTpc = ntTpc;

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
      Int_t pad = hit->GetPad();
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

#if RawCluster
  Int_t nclTpc = 0;
  Int_t remain_nclTpc = 0;
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
      TPCHit* meanHit = cl->GetMeanHit();
      Int_t houghflag = meanHit->GetHoughFlag();
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
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

      if(houghflag!=100&&houghflag!=200) ++remain_nclTpc; //Clusters without track
    }
  }
  event.nclTpc = nclTpc;
  event.remain_nclTpc = remain_nclTpc;
#endif

  //HF1( 10, ntTpc );
  //if( event.ntTpc == 0 ) return true;

  event.nhtrack.resize( ntTpc );
  event.nhtrackEff.resize( ntTpc );
  event.flag.resize( ntTpc );
  event.isBeam.resize( ntTpc );
  event.isAccidental.resize( ntTpc );
  event.isMultiloop.resize( ntTpc );
  event.fittime.resize( ntTpc );
  event.searchtime.resize( ntTpc );
  event.niteration.resize( ntTpc );
  event.chisqr.resize( ntTpc );
  event.pval.resize( ntTpc );
  event.distTgt.resize( ntTpc );

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
#if TruncatedMean
  event.dEdx_0.resize( ntTpc );
  event.dEdx_10.resize( ntTpc );
  event.dEdx_20.resize( ntTpc );
  event.dEdx_30.resize( ntTpc );
  event.dEdx_40.resize( ntTpc );
  event.dEdx_50.resize( ntTpc );
  event.dEdx_60.resize( ntTpc );
#endif
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
  event.residual_t.resize( ntTpc );
  event.residual_x.resize( ntTpc );
  event.residual_y.resize( ntTpc );
  event.residual_z.resize( ntTpc );
  event.resolution_x.resize( ntTpc );
  event.resolution_y.resize( ntTpc );
  event.resolution_z.resize( ntTpc );
  event.pull.resize( ntTpc );
  event.helix_t.resize( ntTpc );
  event.pathhit.resize(ntTpc);
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

  event.exresidual_t.resize( ntTpc );
  event.exresidual_x.resize( ntTpc );
  event.exresidual_y.resize( ntTpc );
  event.exresidual_z.resize( ntTpc );
  event.intrinsic_residual_t.resize( ntTpc );
  event.intrinsic_residual_x.resize( ntTpc );
  event.intrinsic_residual_y.resize( ntTpc );
  event.intrinsic_residual_z.resize( ntTpc );

  event.chargeIndistinguishable.resize( ntTpc );
  event.chisqr_inverted.resize( ntTpc );
  event.pval_inverted.resize( ntTpc );
  event.helix_cx_inverted.resize( ntTpc );
  event.helix_cy_inverted.resize( ntTpc );
  event.helix_z0_inverted.resize( ntTpc );
  event.helix_r_inverted.resize( ntTpc );
  event.helix_dz_inverted.resize( ntTpc );
  event.mom0_x_inverted.resize( ntTpc );
  event.mom0_y_inverted.resize( ntTpc );
  event.mom0_z_inverted.resize( ntTpc );
  event.mom0_inverted.resize( ntTpc );
  event.pid_inverted.resize( ntTpc );

  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Int_t nhEff = tp->GetNHitsEffective();
    Double_t chisqr = tp->GetChiSquare();
    Double_t pval = 1-ROOT::Math::chisquared_cdf(chisqr*(2*nhEff-5), 2*nhEff-5);
    Double_t helix_cx = tp->Getcx(), helix_cy = tp->Getcy();
    Double_t helix_z0 = tp->Getz0(), helix_r = tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 mom0 = tp->GetMom0();
    Int_t flag = tp->GetFitFlag();
    Int_t isbeam = tp->GetIsBeam();
    Int_t isaccidental = tp->GetIsAccidental();
    Int_t ismultiloop = tp->GetIsMultiloop();
    Int_t charge = tp->GetCharge();
    Int_t pid = tp->GetPid();
    Int_t iteration = tp->GetNIteration();
    Double_t fittime = tp->GetFitTime();
    Double_t searchtime = tp->GetSearchTime();
    Double_t pathlen = tp->GetPath();
    Double_t distTgt = tp->GetClosestDist();

    //HF1(11, nh);
    //HF1(12, chisqr);

#if TruncatedMean
    event.dEdx_0[it]=tp->GetdEdx(1.0);
    event.dEdx_10[it]=tp->GetdEdx(0.9);
    event.dEdx_20[it]=tp->GetdEdx(0.8);
    event.dEdx_30[it]=tp->GetdEdx(0.7);
    event.dEdx_40[it]=tp->GetdEdx(0.6);
    event.dEdx_50[it]=tp->GetdEdx(0.5);
    event.dEdx_60[it]=tp->GetdEdx(0.4);
#endif
    event.nhtrack[it] = nh;
    event.nhtrackEff[it] = nhEff;
    event.flag[it] = flag;
    event.isBeam[it] = isbeam;
    event.isAccidental[it] = isaccidental;
    event.isMultiloop[it] = ismultiloop;
    event.fittime[it] = fittime;
    event.charge[it] = charge;
    event.path[it] = pathlen;
    event.distTgt[it] = distTgt;
    event.chisqr[it] = chisqr;
    event.pval[it] = pval;
    event.niteration[it] = iteration;
    event.searchtime[it] = searchtime;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r ;
    event.helix_dz[it] = helix_dz;
    event.mom0_x[it] = mom0.x();
    event.mom0_y[it] = mom0.y();
    event.mom0_z[it] = mom0.z();
    event.mom0[it] = mom0.Mag();
    event.pid[it] = pid;
    event.dE[it] = tp->GetTrackdE();
    event.dEdx[it] = tp->GetdEdx(truncatedMean);
    event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));

    //HF2(20, event.mom0[it]*event.charge[it], event.dEdx[it]);
    if(event.charge[it]<0){
      //HF2(21, -event.mom0[it]*event.charge[it], event.dEdx[it]);
    } else {
      //HF2(22, event.mom0[it]*event.charge[it], event.dEdx[it]);
    }

    //HF1(15, event.mom0[it]);

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
    event.residual_t[it].resize( nh );
    event.residual_x[it].resize( nh );
    event.residual_y[it].resize( nh );
    event.residual_z[it].resize( nh );
    event.resolution_x[it].resize( nh );
    event.resolution_y[it].resize( nh );
    event.resolution_z[it].resize( nh );
    event.pull[it].resize( nh );
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
    event.exresidual_t[it].resize( nh );
    event.exresidual_x[it].resize( nh );
    event.exresidual_y[it].resize( nh );
    event.exresidual_z[it].resize( nh );
    event.intrinsic_residual_t[it].resize( nh );
    event.intrinsic_residual_x[it].resize( nh );
    event.intrinsic_residual_y[it].resize( nh );
    event.intrinsic_residual_z[it].resize( nh );

    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHitInOrder( ih );
      if( !hit ) continue;

      //HF1( 2, hit->GetHoughDist()); HF1( 3, hit->GetHoughDistY());
      Int_t layer = hit->GetLayer();
      Int_t houghflag = hit->GetHoughFlag();
      Double_t residual = hit->GetResidual();
      const TVector3& resi_vect = hit->GetResidualVect();
      const TVector3& res_vect = hit->GetResolutionVect();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& mom = hit->GetMomentumHelix(charge);

      //HF1(13, layer);

      Double_t clde = hit->GetDe();
      Double_t mrow = hit->GetMRow();
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_mrow[it][ih] = mrow;
      event.alpha[it][ih] = tp->GetAlpha(ih);

      TPCHit *clhit = hit->GetHit();
      TPCCluster *cl = clhit->GetParentCluster();
      Int_t clsize = cl->GetClusterSize();
      //Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();
      /*
      HF1(TPCClHid, clsize);
      HF1(TPCClHid+(layer+1)*1000, clsize);
      HF1(TPCClHid+1, clde);
      HF1(TPCClHid+(layer+1)*1000+1, clde);
      */
      const TPCHitContainer& hc = cl -> GetHitContainer();
      for(const auto& hits : hc){
	if(!hits || !hits->IsGood()) continue;
	const TVector3& pos = hits->GetPosition();
	Double_t de = hits->GetCDe();
	Double_t transDist = TranseverseDistance(hitpos.x(), hitpos.z(), pos.x(), pos.z());
	Double_t ratio = de/clde;
	//HF2(TPCClHid+2, transDist, ratio); HF2(TPCClHid+(layer+1)*1000+2, transDist, ratio);
      }
      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;
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

      Double_t resi_theta = TMath::ATan2(resi_vect.z(), -resi_vect.x());
      TVector3 dir_hit(TMath::Cos(hit->GetTheta()), TMath::Sin(hit->GetTheta()), 0);
      TVector3 dir_resi(TMath::Cos(resi_theta), TMath::Sin(resi_theta), 0);
      Double_t sign = 1.;
      if(dir_hit*dir_resi<0) sign = -1.;
      Double_t resi_t = sign*TMath::Hypot(resi_vect.x(), resi_vect.z());
      event.residual_t[it][ih] = resi_t;
      Double_t res_t = TMath::Hypot(res_vect.x(), res_vect.z());
      event.pull[it][ih] = hypot(resi_t/res_t, resi_vect.y()/resi_vect.y());
      /*
      HF1(TPCInclusiveHid+layer,resi_t);
      HF1(TPCInclusiveHid+100+layer,resi_vect.x());
      HF1(TPCInclusiveHid+200+layer,resi_vect.y());
      HF1(TPCInclusiveHid+300+layer,resi_vect.z());
      HF1(TPCInclusiveHid+1000+layer,resi_t/res_t);
      HF1(TPCInclusiveHid+1100+layer,resi_vect.x()/res_vect.x());
      HF1(TPCInclusiveHid+1200+layer,resi_vect.y()/res_vect.y());
      HF1(TPCInclusiveHid+1300+layer,resi_vect.z()/res_vect.z());

      HF1(TPCInclusiveHid+32,resi_t);
      HF1(TPCInclusiveHid+100+32,resi_vect.x());
      HF1(TPCInclusiveHid+200+32,resi_vect.y());
      HF1(TPCInclusiveHid+300+32,resi_vect.z());
      HF1(TPCInclusiveHid+1000+32,resi_t/res_t);
      HF1(TPCInclusiveHid+1100+32,resi_vect.x()/res_vect.x());
      HF1(TPCInclusiveHid+1200+32,resi_vect.y()/res_vect.y());
      HF1(TPCInclusiveHid+1300+32,resi_vect.z()/res_vect.z());

      HF2(TPCInclusiveHid+10000,layer,resi_t);
      HF2(TPCInclusiveHid+100+10000,layer,resi_vect.x());
      HF2(TPCInclusiveHid+200+10000,layer,resi_vect.y());
      HF2(TPCInclusiveHid+300+10000,layer,resi_vect.z());
      HF2(TPCInclusiveHid+1000+10000,layer,resi_t/res_t);
      HF2(TPCInclusiveHid+1100+10000,layer,resi_vect.x()/res_vect.x());
      HF2(TPCInclusiveHid+1200+10000,layer,resi_vect.y()/res_vect.y());
      HF2(TPCInclusiveHid+1300+10000,layer,resi_vect.z()/res_vect.z());
      */
      if(ExclusiveTracking){
	const TVector3& exresi_vect = hit->GetResidualVectExclusive();
	Double_t exresi_x = exresi_vect.x();
	Double_t exresi_y = exresi_vect.y();
	Double_t exresi_z = exresi_vect.z();
	Double_t exresi_t = sign*hypot(exresi_x,exresi_z);
	event.exresidual_t[it][ih] = exresi_t;
	event.exresidual_x[it][ih] = exresi_x;
	event.exresidual_y[it][ih] = exresi_y;
	event.exresidual_z[it][ih] = exresi_z;
	/*
	HF1(TPCExclusiveHid+layer,exresi_t);
	HF1(TPCExclusiveHid+100+layer,exresi_x);
	HF1(TPCExclusiveHid+200+layer,exresi_y);
	HF1(TPCExclusiveHid+300+layer,exresi_z);
	HF1(TPCExclusiveHid+1000+layer,exresi_t/res_t);
	HF1(TPCExclusiveHid+1100+layer,exresi_x/res_vect.x());
	HF1(TPCExclusiveHid+1200+layer,exresi_y/res_vect.y());
	HF1(TPCExclusiveHid+1300+layer,exresi_z/res_vect.z());

	HF1(TPCExclusiveHid+32,exresi_t);
	HF1(TPCExclusiveHid+100+32,exresi_x);
	HF1(TPCExclusiveHid+200+32,exresi_y);
	HF1(TPCExclusiveHid+300+32,exresi_z);
	HF1(TPCExclusiveHid+1000+32,exresi_t/res_t);
	HF1(TPCExclusiveHid+1100+32,exresi_x/res_vect.x());
	HF1(TPCExclusiveHid+1200+32,exresi_y/res_vect.y());
	HF1(TPCExclusiveHid+1300+32,exresi_z/res_vect.z());

	HF2(TPCExclusiveHid+10000,layer,exresi_t);
	HF2(TPCExclusiveHid+100+10000,layer,exresi_x);
	HF2(TPCExclusiveHid+200+10000,layer,exresi_y);
	HF2(TPCExclusiveHid+300+10000,layer,exresi_z);
	HF2(TPCExclusiveHid+1000+10000,layer,exresi_t/res_t);
	HF2(TPCExclusiveHid+1100+10000,layer,exresi_x/res_vect.x());
	HF2(TPCExclusiveHid+1200+10000,layer,exresi_y/res_vect.y());
	HF2(TPCExclusiveHid+1300+10000,layer,exresi_z/res_vect.z());
	*/
	Double_t intrinsic_resi_t = sqrt(abs(resi_t*exresi_t));
	Double_t intrinsic_resi_x = sqrt(abs(resi_vect.x()*exresi_x));
	Double_t intrinsic_resi_y = sqrt(abs(resi_vect.y()*exresi_y));
	Double_t intrinsic_resi_z = sqrt(abs(resi_vect.z()*exresi_z));
	if(resi_t<0) intrinsic_resi_t*=-1;
	if(resi_vect.x()<0) intrinsic_resi_x*=-1;
	if(resi_vect.y()<0) intrinsic_resi_y*=-1;
	if(resi_vect.z()<0) intrinsic_resi_z*=-1;
	event.intrinsic_residual_t[it][ih] = intrinsic_resi_t;
	event.intrinsic_residual_x[it][ih] = intrinsic_resi_x;
	event.intrinsic_residual_y[it][ih] = intrinsic_resi_y;
	event.intrinsic_residual_z[it][ih] = intrinsic_resi_z;
	/*
	HF1(TPCIntrinsicHid+layer,intrinsic_resi_t);
	HF1(TPCIntrinsicHid+100+layer,intrinsic_resi_x);
	HF1(TPCIntrinsicHid+200+layer,intrinsic_resi_y);
	HF1(TPCIntrinsicHid+300+layer,intrinsic_resi_z);
	HF1(TPCIntrinsicHid+1000+layer,intrinsic_resi_t/res_t);
	HF1(TPCIntrinsicHid+1100+layer,intrinsic_resi_x/res_vect.x());
	HF1(TPCIntrinsicHid+1200+layer,intrinsic_resi_y/res_vect.y());
	HF1(TPCIntrinsicHid+1300+layer,intrinsic_resi_z/res_vect.z());

	HF1(TPCIntrinsicHid+32,intrinsic_resi_t);
	HF1(TPCIntrinsicHid+100+32,intrinsic_resi_x);
	HF1(TPCIntrinsicHid+200+32,intrinsic_resi_y);
	HF1(TPCIntrinsicHid+300+32,intrinsic_resi_z);
	HF1(TPCIntrinsicHid+1000+32,intrinsic_resi_t/res_t);
	HF1(TPCIntrinsicHid+1100+32,intrinsic_resi_x/res_vect.x());
	HF1(TPCIntrinsicHid+1200+32,intrinsic_resi_y/res_vect.y());
	HF1(TPCIntrinsicHid+1300+32,intrinsic_resi_z/res_vect.z());

	HF2(TPCIntrinsicHid+10000,layer,intrinsic_resi_t);
	HF2(TPCIntrinsicHid+100+10000,layer,intrinsic_resi_x);
	HF2(TPCIntrinsicHid+200+10000,layer,intrinsic_resi_y);
	HF2(TPCIntrinsicHid+300+10000,layer,intrinsic_resi_z);
	HF2(TPCIntrinsicHid+1000+10000,layer,intrinsic_resi_t/res_t);
	HF2(TPCIntrinsicHid+1100+10000,layer,intrinsic_resi_x/res_vect.x());
	HF2(TPCIntrinsicHid+1200+10000,layer,intrinsic_resi_y/res_vect.y());
	HF2(TPCIntrinsicHid+1300+10000,layer,intrinsic_resi_z/res_vect.z());
	*/
      }
    }

    //Inverted charge tracks
    TPCLocalTrackHelix *tp_inverted = TPCAna.GetTrackTPCHelixChargeInverted( it );
    if( !tp_inverted ) event.chargeIndistinguishable[it] = 0;
    else{
      Double_t chisqr = tp_inverted->GetChiSquare();
      Double_t pval = 1-ROOT::Math::chisquared_cdf(chisqr*(2*nhEff-5), 2*nhEff-5);
      Double_t helix_cx = tp_inverted->Getcx(), helix_cy = tp_inverted->Getcy();
      Double_t helix_z0 = tp_inverted->Getz0(), helix_r = tp_inverted->Getr();
      Double_t helix_dz = tp_inverted->Getdz();
      TVector3 mom0 = tp_inverted->GetMom0();
      Int_t charge = tp_inverted->GetCharge();
      Int_t pid = tp_inverted->GetPid();

      event.chargeIndistinguishable[it] = 1;
      event.chisqr_inverted[it] = chisqr;
      event.pval_inverted[it] = pval;
      event.helix_cx_inverted[it] = helix_cx;
      event.helix_cy_inverted[it] = helix_cy;
      event.helix_z0_inverted[it] = helix_z0;
      event.helix_r_inverted[it] = helix_r ;
      event.helix_dz_inverted[it] = helix_dz;
      event.mom0_x_inverted[it] = mom0.x();
      event.mom0_y_inverted[it] = mom0.y();
      event.mom0_z_inverted[it] = mom0.z();
      event.mom0_inverted[it] = mom0.Mag();
      event.pid_inverted[it] = pid;
      continue;
    }
  }

  Int_t nvtxTpc = TPCAna.GetNVerticesTPC();
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

  event.isLambda.resize(nvtxTpc);
  event.ncombiLambda.resize(nvtxTpc);
  event.distLambda.resize(nvtxTpc);
  event.angleLambda.resize(nvtxTpc);
  event.bestmassLambda.resize(nvtxTpc);
  event.massLambda.resize(nvtxTpc);
  event.vtxLambda_x.resize(nvtxTpc);
  event.vtxLambda_y.resize(nvtxTpc);
  event.vtxLambda_z.resize(nvtxTpc);
  event.momLambda.resize(nvtxTpc);
  event.momLambda_x.resize(nvtxTpc);
  event.momLambda_y.resize(nvtxTpc);
  event.momLambda_z.resize(nvtxTpc);
  event.decaysidLambda.resize(nvtxTpc);
  event.decaysmomLambda.resize(nvtxTpc);
  event.decaysmomLambda_x.resize(nvtxTpc);
  event.decaysmomLambda_y.resize(nvtxTpc);
  event.decaysmomLambda_z.resize(nvtxTpc);
  for( Int_t it=0; it<nvtxTpc; ++it ){
    TPCVertex *vp = TPCAna.GetTPCVertex( it );
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

    event.vtxid[it][0] = vp -> GetTrackId(0);
    event.vtxmom_theta[it][0] = vp -> GetTrackTheta(0);
    event.vtxpos_x[it][0] = vp -> GetTrackPos(0).x();
    event.vtxpos_y[it][0] = vp -> GetTrackPos(0).y();
    event.vtxpos_z[it][0] = vp -> GetTrackPos(0).z();
    event.vtxmom_x[it][0] = vp -> GetTrackMom(0).x();
    event.vtxmom_y[it][0] = vp -> GetTrackMom(0).y();
    event.vtxmom_z[it][0] = vp -> GetTrackMom(0).z();

    event.vtxid[it][1] = vp -> GetTrackId(1);
    event.vtxmom_theta[it][1] = vp -> GetTrackTheta(1);
    event.vtxpos_x[it][1] = vp -> GetTrackPos(1).x();
    event.vtxpos_y[it][1] = vp -> GetTrackPos(1).y();
    event.vtxpos_z[it][1] = vp -> GetTrackPos(1).z();
    event.vtxmom_x[it][1] = vp -> GetTrackMom(1).x();
    event.vtxmom_y[it][1] = vp -> GetTrackMom(1).y();
    event.vtxmom_z[it][1] = vp -> GetTrackMom(1).z();

    event.isLambda[it] = vp -> GetIsLambda();
    event.distLambda[it] = vp -> GetClosestDistLambda();
    event.angleLambda[it] = vp -> GetOpeningAngleLambda();
    if(event.isLambda[it]){
      Int_t ncombi = event.ncombiLambda[it] = vp -> GetNcombiLambda();
      Double_t best_lmass = 9999;
      for( Int_t combi=0; combi<ncombi; ++combi ){

	Double_t lmass = vp -> GetMassLambda(combi);
	event.massLambda[it].push_back(lmass);
	Double_t diff = TMath::Abs(lmass - LambdaMass);
	Double_t best_diff = TMath::Abs(best_lmass - LambdaMass);
	if(diff < best_diff) best_lmass = lmass;

	TVector3 vtx = vp -> GetVertexLambda(combi);
	event.vtxLambda_x[it].push_back(vtx.x());
	event.vtxLambda_y[it].push_back(vtx.y());
	event.vtxLambda_z[it].push_back(vtx.z());

	TVector3 lmom = vp -> GetMomLambda(combi);
	event.momLambda[it].push_back(lmom.Mag());
	event.momLambda_x[it].push_back(lmom.x());
	event.momLambda_y[it].push_back(lmom.y());
	event.momLambda_z[it].push_back(lmom.z());

	Int_t pid = vp -> GetProtonIdLambda(combi);
	TVector3 pmom = vp -> GetProtonMomLambda(combi);
	event.decaysidLambda[it].push_back(pid);
	event.decaysmomLambda[it].push_back(pmom.Mag());
	event.decaysmomLambda_x[it].push_back(pmom.x());
	event.decaysmomLambda_y[it].push_back(pmom.y());
	event.decaysmomLambda_z[it].push_back(pmom.z());

	Int_t piid = vp -> GetPionIdLambda(combi);
	TVector3 pimom = vp -> GetPionMomLambda(combi);
	event.decaysidLambda[it].push_back(piid);
	event.decaysmomLambda[it].push_back(pimom.Mag());
	event.decaysmomLambda_x[it].push_back(pimom.x());
	event.decaysmomLambda_y[it].push_back(pimom.y());
	event.decaysmomLambda_z[it].push_back(pimom.z());
      }
      event.bestmassLambda[it] = best_lmass;
    }
  }

  Int_t nvtxTpcClustered = TPCAna.GetNVerticesTPCClustered();
  event.nvtxTpcClustered = nvtxTpcClustered;
  event.Clusteredvtx_x.resize(nvtxTpcClustered);
  event.Clusteredvtx_y.resize(nvtxTpcClustered);
  event.Clusteredvtx_z.resize(nvtxTpcClustered);
  event.Clusteredvtxid.resize(nvtxTpcClustered);
  for( Int_t ivtx=0; ivtx<nvtxTpcClustered; ++ivtx ){
    TPCVertex *vp = TPCAna.GetTPCVertexClustered( ivtx );
    if( !vp ) continue;
    event.Clusteredvtx_x[ivtx] = vp -> GetVertex().x();
    event.Clusteredvtx_y[ivtx] = vp -> GetVertex().y();
    event.Clusteredvtx_z[ivtx] = vp -> GetVertex().z();

    Int_t ntracks = vp -> GetNTracks();
    event.Clusteredvtxid[ivtx].resize(ntracks);
    for( Int_t it=0; it<ntracks; ++it ){
      event.Clusteredvtxid[ivtx][it] = vp -> GetTrackId(it);
    }
  }

#if TrackSearchFailed
  Int_t failed_ntTpc = TPCAna.GetNTracksTPCHelixFailed();
  event.failed_ntTpc = failed_ntTpc;
  event.failed_nhtrack.resize( failed_ntTpc );
  event.failed_flag.resize( failed_ntTpc );
  event.failed_isBeam.resize( failed_ntTpc );
  event.failed_isKurama.resize( failed_ntTpc );
  event.failed_nclbeforetgt.resize( failed_ntTpc );
  event.failed_isAccidental.resize( failed_ntTpc );
  event.failed_fittime.resize( failed_ntTpc );
  event.failed_searchtime.resize( failed_ntTpc );
  event.failed_niteration.resize( failed_ntTpc );

  event.failed_helix_cx.resize( failed_ntTpc );
  event.failed_helix_cy.resize( failed_ntTpc );
  event.failed_helix_z0.resize( failed_ntTpc );
  event.failed_helix_r.resize( failed_ntTpc );
  event.failed_helix_dz.resize( failed_ntTpc );
  event.failed_mom0.resize( failed_ntTpc );
  event.failed_charge.resize( failed_ntTpc );

  event.failed_hitlayer.resize( failed_ntTpc );
  event.failed_hitpos_x.resize( failed_ntTpc );
  event.failed_hitpos_y.resize( failed_ntTpc );
  event.failed_hitpos_z.resize( failed_ntTpc );
  event.failed_calpos_x.resize( failed_ntTpc );
  event.failed_calpos_y.resize( failed_ntTpc );
  event.failed_calpos_z.resize( failed_ntTpc );
  event.failed_helix_t.resize( failed_ntTpc );
  event.failed_residual.resize( failed_ntTpc );
  event.failed_residual_x.resize( failed_ntTpc );
  event.failed_residual_y.resize( failed_ntTpc );
  event.failed_residual_z.resize( failed_ntTpc );
  event.failed_track_cluster_de.resize( failed_ntTpc );
  event.failed_track_cluster_size.resize( failed_ntTpc );
  event.failed_track_cluster_mrow.resize( failed_ntTpc );

  for( Int_t it=0; it<failed_ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelixFailed( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    Double_t helix_z0=tp->Getz0(), helix_r=tp->Getr();
    Double_t helix_dz=tp->Getdz();
    TVector3 mom0 = tp->GetMom0();
    Int_t flag = tp->GetFitFlag();
    Int_t isbeam = tp->GetIsBeam();
    Int_t nclbeforetgt = tp->GetNclBeforeTgt();
    Int_t isaccidental = tp->GetIsAccidental();
    Int_t fittime = tp->GetFitTime();
    Int_t charge = tp->GetCharge();
    Int_t iteration = tp->GetNIteration();
    event.failed_nhtrack[it] = nh;
    event.failed_flag[it] = flag;
    event.failed_isBeam[it] = isbeam;
    event.failed_nclbeforetgt[it] = nclbeforetgt;
    event.failed_isAccidental[it] = isaccidental;
    event.failed_fittime[it] = fittime;
    event.failed_searchtime[it] = fittime;
    event.failed_niteration[it] = iteration;

    event.failed_helix_cx[it] = helix_cx;
    event.failed_helix_cy[it] = helix_cy;
    event.failed_helix_z0[it] = helix_z0;
    event.failed_helix_r[it] = helix_r ;
    event.failed_helix_dz[it] = helix_dz;
    event.failed_mom0[it] = mom0.Mag();
    event.failed_charge[it] = charge;

    event.failed_hitlayer[it].resize( nh );
    event.failed_hitpos_x[it].resize( nh );
    event.failed_hitpos_y[it].resize( nh );
    event.failed_hitpos_z[it].resize( nh );
    event.failed_calpos_x[it].resize( nh );
    event.failed_calpos_y[it].resize( nh );
    event.failed_calpos_z[it].resize( nh );
    event.failed_helix_t[it].resize( nh );
    event.failed_residual[it].resize( nh );
    event.failed_residual_x[it].resize( nh );
    event.failed_residual_y[it].resize( nh );
    event.failed_residual_z[it].resize( nh );
    event.failed_track_cluster_de[it].resize( nh );
    event.failed_track_cluster_size[it].resize( nh );
    event.failed_track_cluster_mrow[it].resize( nh );

    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& res_vect = hit->GetResidualVect();
      Int_t layer = hit->GetLayer();
      Double_t residual = hit->GetResidual();
      Double_t clde = hit->GetDe();
      Double_t mrow = hit->GetMRow();
      TPCCluster *cl = hit->GetHit()->GetParentCluster();
      Int_t clsize = cl->GetClusterSize();

      event.failed_hitlayer[it][ih] = (double)layer;
      event.failed_hitpos_x[it][ih] = hitpos.x();
      event.failed_hitpos_y[it][ih] = hitpos.y();
      event.failed_hitpos_z[it][ih] = hitpos.z();
      event.failed_calpos_x[it][ih] = calpos.x();
      event.failed_calpos_y[it][ih] = calpos.y();
      event.failed_calpos_z[it][ih] = calpos.z();

      event.failed_residual[it][ih] = residual;
      event.failed_residual_x[it][ih] = res_vect.x();
      event.failed_residual_y[it][ih] = res_vect.y();
      event.failed_residual_z[it][ih] = res_vect.z();
      event.failed_track_cluster_de[it][ih] = clde;
      event.failed_track_cluster_size[it][ih] = clsize;
      event.failed_track_cluster_mrow[it][ih] = mrow;
    }
  }
#endif
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);

  //HF1( 4, sec.count() );
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

#if SaveHistograms
  HB1(1, "Status", 21, 0., 21. );
  HB1(2, "Hough Dist [mm]", 500, 0., 50 );
  HB1(3, "Hough DistY [mm]", 1000, 0., 100 );
  HB1(4, "Process Time [sec]", 1000, 0., 100 );

  HB1(10, "NTrack TPC", 40, 0., 40. );
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 500.);
  HB1(13, "LayerId TPC", 35, 0., 35.);
  HB1(14, "pHS", 1000, 0., 2.5);
  HB1(15, "mom0", 1000, 0., 2.5);
  HB1(16, "mom0-pHS", 1000, -1.25, 1.25);

  const Int_t nbinpoq = 1000;
  const Double_t minpoq = -2.0;
  const Double_t maxpoq = 2.0;
  const Int_t nbindedx = 1000;
  const Double_t mindedx = 0.;
  const Double_t maxdedx = 350.;

  HB2(20, "<dE/dx>;p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2(21, "<dE/dx>;-p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HB2(22, "<dE/dx>;+p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq/2, 0.0, maxpoq, nbindedx, mindedx, maxdedx);
  HBProf(100, ";#beta#gamma;<-dE/dx> [MeVg^{-1}cm^{2}]", 1000000, 0.1, 10000, 1., 10. );
  HBProf(101, "#pi, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(102, "K, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(103, "p, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(104, "#Xi^{-}, CH_{2}; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(105, "#pi, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(106, "K, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(107, "p, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  HBProf(108, "#Xi^{-}, Carbon; Momentum [GeV/#font[12]{c}]; Stopping Power [MeVcm^{-1}]", 100000, 0.1, maxpoq, 1., 200. );
  const Double_t npoints = 1000000;
  for(Int_t i=0; i<npoints; ++i){
    Double_t x = (Double_t) i/npoints;
    Double_t x1 = TMath::Power(10., 5.*x - 1.);
    Double_t beta = TMath::Sqrt(x1*x1/(x1*x1+1.));
    Double_t Carbon = 3.223;
    HFProf(100, x1, Kinematics::HypTPCdEdx(2, 1000.*pdg::KaonMass(), beta)/Carbon);

    Double_t x2 = TMath::Power(10., (TMath::Log10(maxpoq) - TMath::Log10(0.1))*x + TMath::Log10(0.1));
    HFProf(101, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::PionMass(), x2/TMath::Sqrt(x2*x2 + pdg::PionMass()*pdg::PionMass())));
    HFProf(102, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::KaonMass(), x2/TMath::Sqrt(x2*x2 + pdg::KaonMass()*pdg::KaonMass())));
    HFProf(103, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::ProtonMass(), x2/TMath::Sqrt(x2*x2 + pdg::ProtonMass()*pdg::ProtonMass())));
    HFProf(104, x2, Kinematics::HypTPCdEdx(1, 1000.*pdg::XiMinusMass(), x2/TMath::Sqrt(x2*x2 + pdg::XiMinusMass()*pdg::XiMinusMass())));
    HFProf(105, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::PionMass(), x2/TMath::Sqrt(x2*x2 + pdg::PionMass()*pdg::PionMass())));
    HFProf(106, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::KaonMass(), x2/TMath::Sqrt(x2*x2 + pdg::KaonMass()*pdg::KaonMass())));
    HFProf(107, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::ProtonMass(), x2/TMath::Sqrt(x2*x2 + pdg::ProtonMass()*pdg::ProtonMass())));
    HFProf(108, x2, Kinematics::HypTPCdEdx(2, 1000.*pdg::XiMinusMass(), x2/TMath::Sqrt(x2*x2 + pdg::XiMinusMass()*pdg::XiMinusMass())));
  }


#if TrackClusterHist
  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 2000.;

  const Int_t NbinClSize = 25;
  const Double_t MinClSize = 0;
  const Double_t MaxClSize = 25;
  const Int_t NbinDist = 60;
  const Double_t MinDist = -15.;
  const Double_t MaxDist = 15.;
  const Int_t NbinRatio = 100;
  const Double_t MinRatio = 0.;
  const Double_t MaxRatio = 1.;

  HB1(TPCClHid, "Cluster size;Cluster size;Counts", NbinClSize, MinClSize, MaxClSize);
  HB1(TPCClHid+1, "Cluster dE;Cluster dE;Counts", NbinDe, MinDe, MaxDe);
  HB2(TPCClHid+2, "Transverse diffusion;X_{cluster_center}-X_{pad};A/A_{sum}", NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCClHid+(layer+1)*1000, Form("Cluster size layer%d;Cluster size;Counts",layer), NbinClSize, MinClSize, MaxClSize);
    HB1(TPCClHid+(layer+1)*1000+1, Form("Cluster dE layer%d;Cluster dE;Counts",layer), NbinDe, MinDe, MaxDe);
    HB2(TPCClHid+(layer+1)*1000+2, Form("Transverse diffusion Layer%d;X_{cluster_center}-X_{pad};A/A_{sum}",layer), NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  }
#endif

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(TPCExclusiveHid+layer,	Form("TPC ExclusiveResidual T[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+100+layer,	Form("TPC ExclusiveResidual X[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+200+layer,	Form("TPC ExclusiveResidual Y[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+300+layer,	Form("TPC ExclusiveResidual Z[mm];layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+layer,	Form("TPC ExclusivePull T;layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+100+layer,	Form("TPC ExclusivePull X;layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+200+layer,	Form("TPC ExclusivePull Y;layer%d",layer),1000,-10,10);
    HB1(TPCExclusiveHid+1000+300+layer,	Form("TPC ExclusivePull Z;layer%d",layer),1000,-10,10);

    HB1(TPCIntrinsicHid+layer,	Form("TPC IntrinsicResidual T[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+100+layer,	Form("TPC IntrinsicResidual X[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+200+layer,	Form("TPC IntrinsicResidual Y[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+300+layer,	Form("TPC IntrinsicResidual Z[mm];layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+layer,	Form("TPC IntrinsicPull T;layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+100+layer,	Form("TPC IntrinsicPull X;layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+200+layer,	Form("TPC IntrinsicPull Y;layer%d",layer),1000,-10,10);
    HB1(TPCIntrinsicHid+1000+300+layer,	Form("TPC IntrinsicPull Z;layer%d",layer),1000,-10,10);
  }

  HB1(TPCInclusiveHid+32,	Form("TPC InclusiveResidual T[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+100+32,	Form("TPC InclusiveResidual X[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+200+32,	Form("TPC InclusiveResidual Y[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+300+32,	Form("TPC InclusiveResidual Z[mm];AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+32,	Form("TPC InclusivePull T;AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+100+32,	Form("TPC InclusivePull X;AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+200+32,	Form("TPC InclusivePull Y;AllLayer"),1000,-10,10);
  HB1(TPCInclusiveHid+1000+300+32,	Form("TPC InclusivePull Z;AllLayer"),1000,-10,10);

  HB1(TPCExclusiveHid+32,	Form("TPC ExclusiveResidual T[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+100+32,	Form("TPC ExclusiveResidual X[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+200+32,	Form("TPC ExclusiveResidual Y[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+300+32,	Form("TPC ExclusiveResidual Z[mm];AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+32,	Form("TPC ExclusivePull T;AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+100+32,	Form("TPC ExclusivePull X;AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+200+32,	Form("TPC ExclusivePull Y;AllLayer"),1000,-10,10);
  HB1(TPCExclusiveHid+1000+300+32,	Form("TPC ExclusivePull Z;AllLayer"),1000,-10,10);

  HB2(TPCInclusiveHid+10000,	Form("TPC InclusiveResidual T[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+100+10000,	Form("TPC InclusiveResidual X[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+200+10000,	Form("TPC InclusiveResidual Y[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+300+10000,	Form("TPC InclusiveResidual Z[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+10000,	Form("TPC InclusivePull T:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+100+10000,	Form("TPC InclusivePull X:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+200+10000,	Form("TPC InclusivePull Y:Layer"),31,0,31,1000,-10,10);
  HB2(TPCInclusiveHid+1000+300+10000,	Form("TPC InclusivePull Z:Layer"),31,0,31,1000,-10,10);

  HB2(TPCExclusiveHid+10000,	Form("TPC ExclusiveResidual T[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+100+10000,	Form("TPC ExclusiveResidual X[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+200+10000,	Form("TPC ExclusiveResidual Y[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+300+10000,	Form("TPC ExclusiveResidual Z[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+10000,	Form("TPC ExclusivePull T:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+100+10000,	Form("TPC ExclusivePull X:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+200+10000,	Form("TPC ExclusivePull Y:Layer"),31,0,31,1000,-10,10);
  HB2(TPCExclusiveHid+1000+300+10000,	Form("TPC ExclusivePull Z:Layer"),31,0,31,1000,-10,10);

  HB2(TPCIntrinsicHid+10000,	Form("TPC IntrinsicResidual T[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+100+10000,	Form("TPC IntrinsicResidual X[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+200+10000,	Form("TPC IntrinsicResidual Y[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+300+10000,	Form("TPC IntrinsicResidual Z[mm]:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+10000,	Form("TPC IntrinsicPull T:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+100+10000,	Form("TPC IntrinsicPull X:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+200+10000,	Form("TPC IntrinsicPull Y:Layer"),31,0,31,1000,-10,10);
  HB2(TPCIntrinsicHid+1000+300+10000,	Form("TPC IntrinsicPull Z:Layer"),31,0,31,1000,-10,10);
#endif

  //HBTree( "tpc", "tree of DstTPCHelixTracking" );
  tree = new TTree("tpc", "tree of DstTPCHelixTracking");
  tree->Branch( "run_number", &event.run_number );
  tree->Branch( "event_number", &event.event_number );
  tree->Branch( "trig_pat", &event.trig_pat );
  tree->Branch( "trig_flag", &event.trig_flag );
  tree->Branch( "clkTpc", &event.clkTpc);
  tree->Branch( "cobo_id", &event.cobo_id);

  tree->Branch( "nhTpc", &event.nhTpc );
#if SaveRawHit
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
#endif

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "remain_nclTpc", &event.remain_nclTpc );
#if SaveCluster
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
#endif

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "nhtrackEff", &event.nhtrackEff );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
  tree->Branch( "flag", &event.flag );
  tree->Branch( "fittime", &event.fittime );
  tree->Branch( "searchtime", &event.searchtime );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  tree->Branch( "distTgt", &event.distTgt );
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
#if TruncatedMean
  tree->Branch( "dEdx_0", &event.dEdx_0 );
  tree->Branch( "dEdx_10", &event.dEdx_10 );
  tree->Branch( "dEdx_20", &event.dEdx_20 );
  tree->Branch( "dEdx_30", &event.dEdx_30 );
  tree->Branch( "dEdx_40", &event.dEdx_40 );
  tree->Branch( "dEdx_50", &event.dEdx_50 );
  tree->Branch( "dEdx_60", &event.dEdx_60 );
#endif
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
  tree->Branch( "residual_t", &event.residual_t );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);
  tree->Branch( "pull", &event.pull);
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "houghflag", &event.houghflag );
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
  tree->Branch( "exresidual_t", &event.exresidual_t );
  tree->Branch( "exresidual_x", &event.exresidual_x );
  tree->Branch( "exresidual_y", &event.exresidual_y );
  tree->Branch( "exresidual_z", &event.exresidual_z );
  tree->Branch( "intrinsic_residual_t", &event.intrinsic_residual_t );
  tree->Branch( "intrinsic_residual_x", &event.intrinsic_residual_x );
  tree->Branch( "intrinsic_residual_y", &event.intrinsic_residual_y );
  tree->Branch( "intrinsic_residual_z", &event.intrinsic_residual_z );

  tree->Branch( "chargeIndistinguishable", &event.chargeIndistinguishable );
  tree->Branch( "chisqr_inverted", &event.chisqr_inverted );
  tree->Branch( "pval_inverted", &event.pval_inverted );
  tree->Branch( "helix_cx_inverted", &event.helix_cx_inverted );
  tree->Branch( "helix_cy_inverted", &event.helix_cy_inverted );
  tree->Branch( "helix_z0_inverted", &event.helix_z0_inverted );
  tree->Branch( "helix_r_inverted", &event.helix_r_inverted );
  tree->Branch( "helix_dz_inverted", &event.helix_dz_inverted );
  tree->Branch( "mom0_x_inverted", &event.mom0_x_inverted );
  tree->Branch( "mom0_y_inverted", &event.mom0_y_inverted );
  tree->Branch( "mom0_z_inverted", &event.mom0_z_inverted );
  tree->Branch( "mom0_inverted", &event.mom0_inverted );
  tree->Branch( "pid_inverted", &event.pid_inverted );

#if TrackSearchFailed
  tree->Branch( "failed_ntTpc", &event.failed_ntTpc );
  tree->Branch( "failed_nhtrack", &event.failed_nhtrack );
  tree->Branch( "failed_flag", &event.failed_flag );
  tree->Branch( "failed_isBeam", &event.failed_isBeam );
  tree->Branch( "failed_nclbeforetgt", &event.failed_nclbeforetgt );
  tree->Branch( "failed_isAccidental", &event.failed_isAccidental );
  tree->Branch( "failed_fittime", &event.failed_fittime );
  tree->Branch( "failed_searchtime", &event.failed_searchtime );
  tree->Branch( "failed_niteration", &event.failed_niteration );

  tree->Branch( "failed_helix_cx", &event.failed_helix_cx );
  tree->Branch( "failed_helix_cy", &event.failed_helix_cy );
  tree->Branch( "failed_helix_z0", &event.failed_helix_z0 );
  tree->Branch( "failed_helix_r", &event.failed_helix_r );
  tree->Branch( "failed_helix_dz", &event.failed_helix_dz );
  tree->Branch( "failed_mom0", &event.failed_mom0 );
  tree->Branch( "failed_charge", &event.failed_charge );

  tree->Branch( "failed_hitlayer", &event.failed_hitlayer );
  tree->Branch( "failed_hitpos_x", &event.failed_hitpos_x );
  tree->Branch( "failed_hitpos_y", &event.failed_hitpos_y );
  tree->Branch( "failed_hitpos_z", &event.failed_hitpos_z );
  tree->Branch( "failed_calpos_x", &event.failed_calpos_x );
  tree->Branch( "failed_calpos_y", &event.failed_calpos_y );
  tree->Branch( "failed_calpos_z", &event.failed_calpos_z );
  tree->Branch( "failed_helix_t", &event.failed_helix_t );
  tree->Branch( "failed_residual", &event.failed_residual );
  tree->Branch( "failed_residual_x", &event.failed_residual_x );
  tree->Branch( "failed_residual_y", &event.failed_residual_y );
  tree->Branch( "failed_residual_z", &event.failed_residual_z );
  tree->Branch( "failed_track_cluster_de", &event.failed_track_cluster_de);
  tree->Branch( "failed_track_cluster_size", &event.failed_track_cluster_size);
  tree->Branch( "failed_track_cluster_mrow", &event.failed_track_cluster_mrow);
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

  tree->Branch( "isLambda", &event.isLambda );
  tree->Branch( "ncombiLambda", &event.ncombiLambda );
  tree->Branch( "distLambda", &event.distLambda );
  tree->Branch( "angleLambda", &event.angleLambda );
  tree->Branch( "bestmassLambda", &event.bestmassLambda );
  tree->Branch( "massLambda", &event.massLambda );
  tree->Branch( "vtxLambda_x", &event.vtxLambda_x );
  tree->Branch( "vtxLambda_y", &event.vtxLambda_y );
  tree->Branch( "vtxLambda_z", &event.vtxLambda_z );
  tree->Branch( "momLambda", &event.momLambda );
  tree->Branch( "momLambda_x", &event.momLambda_x );
  tree->Branch( "momLambda_y", &event.momLambda_y );
  tree->Branch( "momLambda_z", &event.momLambda_z );
  tree->Branch( "decaysidLambda", &event.decaysidLambda );
  tree->Branch( "decaysmomLambda", &event.decaysmomLambda );
  tree->Branch( "decaysmomLambda_x", &event.decaysmomLambda_x );
  tree->Branch( "decaysmomLambda_y", &event.decaysmomLambda_y );
  tree->Branch( "decaysmomLambda_z", &event.decaysmomLambda_z );

  tree->Branch( "nvtxTpcClustered", &event.nvtxTpcClustered );
  tree->Branch( "clusteredVtx_x", &event.Clusteredvtx_x );
  tree->Branch( "clusteredVtx_y", &event.Clusteredvtx_y );
  tree->Branch( "clusteredVtx_z", &event.Clusteredvtx_z );
  tree->Branch( "clusteredVtxid", &event.Clusteredvtxid );

  std::cout<<"initialize!!!!!!!!!!!!!!!!!"<<std::endl;
  TTreeReaderCont[kTpcHit] = new TTreeReader( "tpc", TFileCont[kTpcHit] );
  const auto& reader = TTreeReaderCont[kTpcHit];
  src.run_number = new TTreeReaderValue<UInt_t>( *reader, "run_number" );
  src.event_number = new TTreeReaderValue<UInt_t>( *reader, "event_number" );
  src.trig_flag = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "trig_flag" );
  src.trig_pat = new TTreeReaderValue<std::vector<Double_t>>( *reader, "trig_pat" );
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
  src.cobo_id = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cobo_id");

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
      InitializeParameter<FieldMan>("FLDMAP") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
