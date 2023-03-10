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
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define TrigA 0 //if 1, TrigA is required
#define TrigB 0
#define TrigC 0
#define TrigD 0

#define TrackCluster 1
#define TruncatedMean 0
#define TrackSearchFailed 1

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
const double truncatedMean = 0.8; //80%
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit, kK18Tracking, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]", "[K18Tracking]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "k18track","" };
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

  // K18
  Int_t ntK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> delta_3rd;
  std::vector<Double_t> xout;
  std::vector<Double_t> yout;
  std::vector<Double_t> uout;
  std::vector<Double_t> vout;

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
  std::vector<Int_t> isKurama; // isKurama: 1 = Beam, 0 = Scat
  std::vector<Int_t> flag;
  std::vector<Int_t> fittime;  //sec
  std::vector<Double_t> chisqr;
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

  std::vector<std::vector<Double_t>> combi_id; //track number of combi
  std::vector<std::vector<Double_t>> closeDistTpc;
  std::vector<std::vector<Double_t>> vtxTpc;
  std::vector<std::vector<Double_t>> vtyTpc;
  std::vector<std::vector<Double_t>> vtzTpc;
  std::vector<std::vector<Double_t>> mom_vtx;//Helix momentum at vtx
  std::vector<std::vector<Double_t>> mom_vty;//Helix momentum at vtx
  std::vector<std::vector<Double_t>> mom_vtz;//Helix momentum at vtx

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
  std::vector<std::vector<Double_t>> helix_t;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> theta_diff;
  std::vector<std::vector<Double_t>> houghflag;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  Int_t failed_ntTpc; // Number of Tracks
  std::vector<Int_t> failed_nhtrack;
  std::vector<Int_t> failed_flag;
  std::vector<Int_t> failed_isBeam;
  std::vector<Int_t> failed_isKurama;
  std::vector<Int_t> failed_fittime;  //sec
  std::vector<Double_t> failed_chisqr;
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
  std::vector<std::vector<Double_t>> failed_residual;
  std::vector<std::vector<Double_t>> failed_residual_x;
  std::vector<std::vector<Double_t>> failed_residual_y;
  std::vector<std::vector<Double_t>> failed_residual_z;
  std::vector<std::vector<Double_t>> failed_track_cluster_de;
  std::vector<std::vector<Double_t>> failed_track_cluster_size;
  std::vector<std::vector<Double_t>> failed_track_cluster_mrow;

  void clear( void )
  {
    trigpat.clear();
    trigflag.clear();
    clkTpc.clear();

    runnum = 0;
    evnum = 0;
    status = 0;

    ntK18 = 0;
    chisqrK18.clear();
    pK18.clear();
    delta_3rd.clear();
    xout.clear();
    yout.clear();
    uout.clear();
    vout.clear();

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
    isKurama.clear();
    flag.clear();
    fittime.clear();
    chisqr.clear();
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

    combi_id.clear();
    mom_vtx.clear();
    mom_vty.clear();
    mom_vtz.clear();

    pid.clear();
    vtxTpc.clear();
    vtyTpc.clear();
    vtzTpc.clear();
    closeDistTpc.clear();

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
    helix_t.clear();
    pathhit.clear();
    theta_diff.clear();
    houghflag.clear();

    track_cluster_de.clear();
    track_cluster_size.clear();
    track_cluster_mrow.clear();
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();

    failed_ntTpc = 0;
    failed_nhtrack.clear();
    failed_flag.clear();
    failed_isBeam.clear();
    failed_isKurama.clear();
    failed_fittime.clear();
    failed_chisqr.clear();
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
    failed_residual.clear();
    failed_residual_x.clear();
    failed_residual_y.clear();
    failed_residual_z.clear();
    failed_track_cluster_de.clear();
    failed_track_cluster_size.clear();
    failed_track_cluster_mrow.clear();
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
  TTreeReaderValue<std::vector<Double_t>>* cdeTpc;     // cdE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;      //clock time

  // K18
  Int_t ntK18;
  Double_t chisqrK18[MaxHits];
  Double_t pK18[MaxHits];
  Double_t delta_3rd[MaxHits];
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
    TPCGainHid    = 100000,
    TPCClHid = 200000,
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
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

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

  if( ievent%1000==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.clkTpc = **src.clkTpc;

#if TrigA
  if(event.trigflag[20]<0) return true;
#endif
#if TrigB
  if(event.trigflag[21]<0) return true;
#endif
#if TrigC
  if(event.trigflag[22]<0) return true;
#endif
#if TrigD
  if(event.trigflag[23]<0) return true;
#endif

  event.ntK18 = src.ntK18;
  for(int it=0; it<src.ntK18; ++it){
    event.chisqrK18.push_back(src.chisqrK18[it]);
    event.pK18.push_back(src.pK18[it]);
    event.delta_3rd.push_back(src.delta_3rd[it]);
    event.xout.push_back(src.xout[it]);
    event.yout.push_back(src.yout[it]);
    event.uout.push_back(src.uout[it]);
    event.vout.push_back(src.vout[it]);
    HF1(14, src.pK18[it]);
  }

  HF1( 1, event.status++ );

  if( **src.nhTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  if(event.clkTpc.size() != 1){
    std::cerr << "something is wrong" << std::endl;
    return true;
  }
  Double_t clock = event.clkTpc.at(0);
  DCAnalyzer DCAna;
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, clock);

  DCAna.TrackSearchTPCHelix();

  HF1( 1, event.status++ );
  Int_t nh_Tpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = DCAna.GetTPCHC( layer );
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
  HF1(1, event.status++);

  Int_t nclTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = DCAna.GetTPCClCont( layer );
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
    }
  }
  event.nclTpc = nclTpc;
  HF1( 1, event.status++ );

  Int_t ntTpc = DCAna.GetNTracksTPCHelix();
  event.ntTpc = ntTpc;
  HF1( 10, ntTpc );
  if( event.ntTpc == 0 )
    return true;

  HF1( 1, event.status++ );
  event.nhtrack.resize( ntTpc );
  event.isBeam.resize( ntTpc );
  event.isKurama.resize( ntTpc );
  event.flag.resize( ntTpc );
  event.fittime.resize( ntTpc );
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
  event.combi_id.resize( ntTpc );
  event.closeDistTpc.resize( ntTpc );
  event.vtxTpc.resize( ntTpc );
  event.vtyTpc.resize( ntTpc );
  event.vtzTpc.resize( ntTpc );
  event.mom_vtx.resize( ntTpc );
  event.mom_vty.resize( ntTpc );
  event.mom_vtz.resize( ntTpc );

  event.hitlayer.resize( ntTpc );
  event.hitpos_x.resize( ntTpc );
  event.hitpos_y.resize( ntTpc );
  event.hitpos_z.resize( ntTpc );
  event.calpos_x.resize( ntTpc );
  event.calpos_y.resize( ntTpc );
  event.calpos_z.resize( ntTpc );
  event.residual.resize( ntTpc );
  event.residual_x.resize( ntTpc );
  event.residual_y.resize( ntTpc );
  event.residual_z.resize( ntTpc );
  event.helix_t.resize( ntTpc );
  event.pathhit.resize(ntTpc);
  event.theta_diff.resize(ntTpc);
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
    TPCLocalTrackHelix *tp = DCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t helix_cx = tp->Getcx(), helix_cy = tp->Getcy();
    Double_t helix_z0 = tp->Getz0(), helix_r = tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
    Int_t flag = tp->GetFitFlag();
    Int_t isbeam = tp->GetIsBeam();
    Int_t iskurama = tp->GetIsKurama();
    Int_t fittime = tp->GetFitTime();
    Int_t charge = tp->GetCharge();
    Double_t pathlen = tp->GetPath();

    HF1(11, nh);
    HF1(12, chisqr);

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
    event.isBeam[it] = isbeam;
    event.isKurama[it] = iskurama;
    event.flag[it] = flag;
    event.fittime[it] = fittime;
    event.charge[it] = charge;
    event.path[it] = pathlen;
    event.chisqr[it] = chisqr;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r ;
    event.helix_dz[it] = helix_dz;
    event.mom0_x[it] = Mom0.x();
    event.mom0_y[it] = Mom0.y();
    event.mom0_z[it] = Mom0.z();
    event.mom0[it] = Mom0.Mag();
    event.dE[it] = tp->GetTrackdE();
    event.dEdx[it] = tp->GetdEdx(truncatedMean);
    event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));

    HF2(20, event.mom0[it]*event.charge[it], event.dEdx[it]);

    HF1(15, event.mom0[it]);
    if(src.ntK18==1) HF1(16, event.mom0[it]-src.pK18[0]);
    int particleID = Kinematics::HypTPCdEdxPID_temp(event.dEdx[it], event.mom0[it]*event.charge[it]);
    event.pid[it]=particleID;

    event.combi_id[it].resize( ntTpc );
    event.closeDistTpc[it].resize( ntTpc );
    event.vtxTpc[it].resize( ntTpc );
    event.vtyTpc[it].resize( ntTpc );
    event.vtzTpc[it].resize( ntTpc );
    event.mom_vtx[it].resize( ntTpc );
    event.mom_vty[it].resize( ntTpc );
    event.mom_vtz[it].resize( ntTpc );
    event.hitlayer[it].resize( nh );
    event.hitpos_x[it].resize( nh );
    event.hitpos_y[it].resize( nh );
    event.hitpos_z[it].resize( nh );
    event.calpos_x[it].resize( nh );
    event.calpos_y[it].resize( nh );
    event.calpos_z[it].resize( nh );
    event.residual[it].resize( nh );
    event.residual_x[it].resize( nh );
    event.residual_y[it].resize( nh );
    event.residual_z[it].resize( nh );
    event.helix_t[it].resize( nh );
    event.pathhit[it].resize(nh);
    event.theta_diff[it].resize(nh);
    event.houghflag[it].resize(nh);
    event.track_cluster_de[it].resize(nh);
    event.track_cluster_size[it].resize(nh);
    event.track_cluster_mrow[it].resize(nh);
    event.track_cluster_de_center[it].resize(nh);
    event.track_cluster_x_center[it].resize(nh);
    event.track_cluster_y_center[it].resize(nh);
    event.track_cluster_z_center[it].resize(nh);
    event.track_cluster_row_center[it].resize(nh);

    double par1[5]={helix_cx, helix_cy, helix_z0,
		    helix_r, helix_dz};
    for( Int_t it2=0; it2<ntTpc; ++it2 ){
      if(it2==it) continue;
      TPCLocalTrackHelix *tp2 = DCAna.GetTrackTPCHelix( it2 );
      if( !tp2 ) continue;
      Double_t helix_cx2 = tp2->Getcx(), helix_cy2 = tp2->Getcy();
      Double_t helix_z02 = tp2->Getz0(), helix_r2 = tp2->Getr();
      Double_t helix_dz2 = tp2->Getdz();

      double par2[5]={helix_cx2, helix_cy2, helix_z02,
		      helix_r2, helix_dz2};
      double closeDistTpc, t1, t2;
      TVector3 vert = Kinematics::VertexPointHelix(par1, par2,
						   closeDistTpc, t1, t2);
      event.combi_id[it][it2] = (double)it2;
      event.closeDistTpc[it][it2] = closeDistTpc;
      event.vtxTpc[it][it2] = vert.x();
      event.vtyTpc[it][it2] = vert.y();
      event.vtzTpc[it][it2] = vert.z();

      TVector3 mom_vtx = tp->CalcHelixMom(par1, vert.y());
      event.mom_vtx[it][it2] = mom_vtx.x();
      event.mom_vty[it][it2] = mom_vtx.y();
      event.mom_vtz[it][it2] = mom_vtx.z();
    }

    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      HF1( 2, hit->GetHoughDist());
      HF1( 3, hit->GetHoughDistY());
      Int_t layer = hit->GetLayer();
      Int_t houghflag = hit->GetHoughFlag();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& res_vect = hit->GetResidualVect();
      HF1(13, layer);

#if TrackCluster
      TPCHit *clhit = hit->GetHit();
      Double_t clde = hit->GetDe();
      TPCCluster *cl = clhit->GetParentCluster();
      Int_t clsize = cl->GetClusterSize();
      Double_t mrow = cl->MeanRow(); // same
      TPCHit* centerHit = cl->GetCenterHit();
      const TVector3& centerPos = centerHit->GetPosition();
      Double_t centerDe = centerHit->GetCDe();
      Int_t centerRow = centerHit->GetRow();

      HF1(TPCClHid, clsize);
      HF1(TPCClHid+(layer+1)*1000, clsize);
      HF1(TPCClHid+1, clde);
      HF1(TPCClHid+(layer+1)*1000+1, clde);
      const TPCHitContainer& hc = cl -> GetHitContainer();
      for(const auto& hits : hc){
	if(!hits || !hits->IsGood()) continue;
	const TVector3& pos = hits->GetPosition();
	Double_t de = hits->GetCDe();
	Double_t transDist = TranseverseDistance(hitpos.x(), hitpos.z(), pos.x(), pos.z());
	Double_t ratio = de/clde;
	HF2(TPCClHid+2, transDist, ratio);
	HF2(TPCClHid+(layer+1)*1000+2, transDist, ratio);
      }
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_size[it][ih] = clsize;
      event.track_cluster_mrow[it][ih] = mrow;
      event.track_cluster_de_center[it][ih] = centerDe;
      event.track_cluster_x_center[it][ih] = centerPos.X();
      event.track_cluster_y_center[it][ih] = centerPos.Y();
      event.track_cluster_z_center[it][ih] = centerPos.Z();
      event.track_cluster_row_center[it][ih] = centerRow;
#endif
      event.theta_diff[it][ih] = hit->GetPadTrackAngleHelix();
      event.pathhit[it][ih] = hit->GetPathHelix();

      Double_t residual = hit->GetResidual();
      event.hitlayer[it][ih] = (double)layer;
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
      event.helix_t[it][ih] = hit->GetTcal();
      event.houghflag[it][ih] = houghflag;
    }
  }

#if TrackSearchFailed
  Int_t failed_ntTpc = DCAna.GetNTracksTPCHelixFailed();
  event.failed_ntTpc = failed_ntTpc;
  event.failed_nhtrack.resize( failed_ntTpc );
  event.failed_flag.resize( failed_ntTpc );
  event.failed_isBeam.resize( failed_ntTpc );
  event.failed_isKurama.resize( failed_ntTpc );
  event.failed_fittime.resize( failed_ntTpc );
  event.failed_chisqr.resize( failed_ntTpc );
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
  event.failed_residual.resize( failed_ntTpc );
  event.failed_residual_x.resize( failed_ntTpc );
  event.failed_residual_y.resize( failed_ntTpc );
  event.failed_residual_z.resize( failed_ntTpc );
  event.failed_track_cluster_de.resize( failed_ntTpc );
  event.failed_track_cluster_size.resize( failed_ntTpc );
  event.failed_track_cluster_mrow.resize( failed_ntTpc );

  for( Int_t it=0; it<failed_ntTpc; ++it ){
    TPCLocalTrackHelix *tp = DCAna.GetTrackTPCHelixFailed( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    Double_t helix_z0=tp->Getz0(), helix_r=tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
    Int_t flag = tp->GetFitFlag();
    Int_t isbeam = tp->GetIsBeam();
    Int_t iskurama = tp->GetIsKurama();
    Int_t fittime = tp->GetFitTime();
    Int_t charge = tp->GetCharge();
    event.failed_nhtrack[it] = nh;
    event.failed_flag[it] = flag;
    event.failed_isBeam[it] = isbeam;
    event.failed_isKurama[it] = iskurama;
    event.failed_fittime[it] = fittime;
    event.failed_chisqr[it] = chisqr;
    event.failed_helix_cx[it] = helix_cx;
    event.failed_helix_cy[it] = helix_cy;
    event.failed_helix_z0[it] = helix_z0;
    event.failed_helix_r[it] = helix_r ;
    event.failed_helix_dz[it] = helix_dz;
    event.failed_mom0[it] = Mom0.Mag();
    event.failed_charge[it] = charge;

    event.failed_hitlayer[it].resize( nh );
    event.failed_hitpos_x[it].resize( nh );
    event.failed_hitpos_y[it].resize( nh );
    event.failed_hitpos_z[it].resize( nh );
    event.failed_calpos_x[it].resize( nh );
    event.failed_calpos_y[it].resize( nh );
    event.failed_calpos_z[it].resize( nh );
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
  HF1( 1, event.status++ );

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);

  HF1( 4, sec.count() );

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

  HB1(1, "Status", 21, 0., 21. );
  HB1(2, "Hough Dist [mm]", 500, 0., 50 );
  HB1(3, "Hough DistY [mm]", 1000, 0., 100 );
  HB1(4, "Process Time [sec]", 4000, 0., 4000 );

  HB1(10, "NTrack TPC", 40, 0., 40. );
  HB1(11, "#Hits of Track TPC", 50, 0., 50.);
  HB1(12, "Chisqr TPC", 500, 0., 500.);
  HB1(13, "LayerId TPC", 35, 0., 35.);
  HB1(14, "pK18", 1000, 0., 2.5);
  HB1(15, "mom0", 1000, 0., 2.5);
  HB1(16, "mom0-pK18", 1000, -1.25, 1.25);

  const Int_t nbinpoq = 1000;
  const Double_t minpoq = -2.0;
  const Double_t maxpoq = 2.0;
  const Int_t nbindedx = 1000;
  const Double_t mindedx = 0.;
  const Double_t maxdedx = 350.;

  HB2(20, "<dE/dx>;p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);

#if TrackCluster
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

  HBTree( "tpc", "tree of DstTPCHelixK18Tracking" );

  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
  tree->Branch( "clkTpc", &event.clkTpc);

  tree->Branch( "ntK18",      &event.ntK18);
  tree->Branch( "chisqrK18",  &event.chisqrK18);
  tree->Branch( "pK18",       &event.pK18);
  tree->Branch( "delta_3rd",  &event.delta_3rd);
  tree->Branch( "xout",    &event.xout);
  tree->Branch( "yout",    &event.yout);
  tree->Branch( "uout",    &event.uout);
  tree->Branch( "vout",    &event.vout);

  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
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
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "flag", &event.flag );
  tree->Branch( "fittime", &event.fittime );
  tree->Branch( "chisqr", &event.chisqr );
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
  tree->Branch( "combi_id", &event.combi_id );
  tree->Branch( "closeDistTpc", &event.closeDistTpc );
  tree->Branch( "vtxTpc", &event.vtxTpc );
  tree->Branch( "vtyTpc", &event.vtyTpc );
  tree->Branch( "vtzTpc", &event.vtzTpc );
  tree->Branch( "mom_vtx", &event.mom_vtx );
  tree->Branch( "mom_vty", &event.mom_vty );
  tree->Branch( "mom_vtz", &event.mom_vtz );

  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "theta_diff", &event.theta_diff);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "houghflag", &event.houghflag );
#if TrackCluster
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
#endif
#if TrackSearchFailed
  tree->Branch( "failed_ntTpc", &event.failed_ntTpc );
  tree->Branch( "failed_nhtrack", &event.failed_nhtrack );
  tree->Branch( "failed_isBeam", &event.failed_isBeam );
  tree->Branch( "failed_isKurama", &event.failed_isKurama );
  tree->Branch( "failed_flag", &event.failed_flag );
  tree->Branch( "failed_fittime", &event.failed_fittime );
  tree->Branch( "failed_chisqr", &event.failed_chisqr );
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
  tree->Branch( "failed_residual", &event.failed_residual );
  tree->Branch( "failed_residual_x", &event.failed_residual_x );
  tree->Branch( "failed_residual_y", &event.failed_residual_y );
  tree->Branch( "failed_residual_z", &event.failed_residual_z );
  tree->Branch( "failed_track_cluster_de", &event.failed_track_cluster_de);
  tree->Branch( "failed_track_cluster_size", &event.failed_track_cluster_size);
  tree->Branch( "failed_track_cluster_mrow", &event.failed_track_cluster_mrow);
#endif

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

  TTreeCont[kK18Tracking]->SetBranchStatus("*", 0);
  TTreeCont[kK18Tracking]->SetBranchStatus("ntK18",       1);
  TTreeCont[kK18Tracking]->SetBranchStatus("chisqrK18",   1);
  TTreeCont[kK18Tracking]->SetBranchStatus("p_3rd",       1);
  TTreeCont[kK18Tracking]->SetBranchStatus("delta_3rd",   1);
  TTreeCont[kK18Tracking]->SetBranchStatus("xout",        1);
  TTreeCont[kK18Tracking]->SetBranchStatus("yout",        1);
  TTreeCont[kK18Tracking]->SetBranchStatus("uout",        1);
  TTreeCont[kK18Tracking]->SetBranchStatus("vout",        1);

  TTreeCont[kK18Tracking]->SetBranchAddress("ntK18",     &src.ntK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrK18", &src.chisqrK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("p_3rd",     &src.pK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("delta_3rd", &src.delta_3rd);
  TTreeCont[kK18Tracking]->SetBranchAddress("xout",   &src.xout);
  TTreeCont[kK18Tracking]->SetBranchAddress("yout",   &src.yout);
  TTreeCont[kK18Tracking]->SetBranchAddress("uout",   &src.uout);
  TTreeCont[kK18Tracking]->SetBranchAddress("vout",   &src.vout);

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
