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

#define TrackSearch 1
#define Gain_center 0
#define HoughYcut 1
#define DE_padparam 1


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
const double chisqrCut = 10.;//mm
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
  Int_t nh_cluster_Tpc;
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;
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
  std::vector<Int_t> nhtrack_clmulti; // Number of Hits (clsize>1)
  std::vector<Int_t> isBeam; // isBeam: 1 = Beam, 0 = Scat
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx;
  std::vector<Double_t> dEdx2;
  std::vector<Double_t> dEdx_20;
  std::vector<Double_t> dEdx_10;
  std::vector<Double_t> dEdx_rms;
  std::vector<Double_t> dEdx_rms20;
  std::vector<Double_t> dEdx_rms10;

  std::vector<Double_t> dEdx_geo;
  std::vector<Double_t> dEdx_geo_20;
  std::vector<Double_t> dEdx_geo_10;
  std::vector<Double_t> dEdx_harm;
  std::vector<Double_t> dEdx_harm_20;
  std::vector<Double_t> dEdx_harm_10;

  std::vector<Double_t> dE_clmulti;
  std::vector<Double_t> dEdx_clmulti;
  std::vector<Double_t> dEdx_clmulti_20;
  std::vector<Double_t> dEdx_clmulti_10;
  std::vector<Double_t> dE_cor;
  std::vector<Double_t> dEdx_cor;
  std::vector<Double_t> dEdx2_cor;
  std::vector<Double_t> dEdx_cor_20;
  std::vector<Double_t> dEdx_cor_10;


  std::vector<Double_t> dz_factor;
  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0

  std::vector<Double_t> mom0_cor_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_cor_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_cor_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_cor;//Helix momentum at Y = 0


  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path

  // under dev
  std::vector<Int_t> combi_id; //track number of combi
  std::vector<Double_t> mom_vtx;//Helix momentum at vtx
  std::vector<Double_t> mom_vty;//Helix momentum at vtx
  std::vector<Double_t> mom_vtz;//Helix momentum at vtx
  std::vector<Double_t> mom_cor_vtx;//Helix momentum at vtx
  std::vector<Double_t> mom_cor_vty;//Helix momentum at vtx
  std::vector<Double_t> mom_cor_vtz;//Helix momentum at vtx

  std::vector<Int_t> pid;
  std::vector<Double_t> vtx;
  std::vector<Double_t> vty;
  std::vector<Double_t> vtz;
  std::vector<Int_t> chisqr_flag;
  std::vector<Double_t> closeDist;


  std::vector<Double_t> M_Lambda;
  std::vector<Double_t> Lvtx;
  std::vector<Double_t> Lvty;
  std::vector<Double_t> Lvtz;
  std::vector<Double_t> LcloseDist;
  std::vector<Double_t> Mom_Lambda;
  std::vector<Double_t> Mom_Lambdax;
  std::vector<Double_t> Mom_Lambday;
  std::vector<Double_t> Mom_Lambdaz;

  std::vector<Double_t> M_Ks;
  std::vector<Double_t> Ksvtx;
  std::vector<Double_t> Ksvty;
  std::vector<Double_t> Ksvtz;
  std::vector<Double_t> KscloseDist;
  std::vector<Double_t> Mom_Ks;
  std::vector<Double_t> Mom_Ksx;
  std::vector<Double_t> Mom_Ksy;
  std::vector<Double_t> Mom_Ksz;


  std::vector<std::vector<Double_t>> hitlayer;
  std::vector<std::vector<Double_t>> hitmrow;
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
  std::vector<std::vector<Double_t>> hitde;
  std::vector<std::vector<Double_t>> hitClsize;

  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    nhTpc = 0;
    nh_cluster_Tpc = 0;
    raw_hitpos_x.clear();
    raw_hitpos_y.clear();
    raw_hitpos_z.clear();
    raw_de.clear();
    raw_padid.clear();
    raw_layer.clear();
    raw_row.clear();
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
    nhtrack_clmulti.clear();
    isBeam.clear();
    chisqr.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();
    dEdx2.clear();
    dEdx_20.clear();
    dEdx_10.clear();
    dEdx_rms.clear();
    dEdx_rms20.clear();
    dEdx_rms10.clear();

    dEdx_geo.clear();
    dEdx_geo_20.clear();
    dEdx_geo_10.clear();
    dEdx_harm.clear();
    dEdx_harm_20.clear();
    dEdx_harm_10.clear();
    dE_clmulti.clear();
    dEdx_clmulti.clear();
    dEdx_clmulti_20.clear();
    dEdx_clmulti_10.clear();
    dE_cor.clear();
    dEdx_cor.clear();
    dEdx2_cor.clear();
    dEdx_cor_20.clear();
    dEdx_cor_10.clear();

    dz_factor.clear();

    mom0_x.clear();
    mom0_y.clear();
    mom0_z.clear();
    mom0.clear();

    mom0_cor_x.clear();
    mom0_cor_y.clear();
    mom0_cor_z.clear();
    mom0_cor.clear();

    charge.clear();
    path.clear();

    combi_id.clear();
    mom_vtx.clear();
    mom_vty.clear();
    mom_vtz.clear();
    mom_cor_vtx.clear();
    mom_cor_vty.clear();
    mom_cor_vtz.clear();

    pid.clear();
    vtx.clear();
    vty.clear();
    vtz.clear();
    closeDist.clear();
    chisqr_flag.clear();

    M_Lambda.clear();
    Lvtx.clear();
    Lvty.clear();
    Lvtz.clear();
    LcloseDist.clear();
    Mom_Lambda.clear();
    Mom_Lambdax.clear();
    Mom_Lambday.clear();
    Mom_Lambdaz.clear();

    M_Ks.clear();
    Ksvtx.clear();
    Ksvty.clear();
    Ksvtz.clear();
    KscloseDist.clear();
    Mom_Ks.clear();
    Mom_Ksx.clear();
    Mom_Ksy.clear();
    Mom_Ksz.clear();

    hitlayer.clear();
    hitmrow.clear();
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
    hitde.clear();
    hitClsize.clear();
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
  if( ievent%100==0 ){
    //if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  // static const auto KaonMass    = pdg::KaonMass();
  static const auto PionMass    = pdg::PionMass();
  static const auto ProtonMass  = pdg::ProtonMass();


  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  //  event.nhTpc = **src.nhTpc;

  HF1( 1, event.status++ );

  if( **src.nhTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  DCAnalyzer DCAna;
  //  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc);
  //DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.ctTpc, **src.deTpc);
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.ctTpc, **src.cdeTpc);
#if HoughYcut
  DCAna.HoughYCut(min_ycut, max_ycut);
#endif


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


  Int_t nh_cl_Tpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = DCAna.GetTPCClCont( layer );
    for( const auto& hit : hc ){
      if( !hit || !hit->IsGood() )
        continue;
      Double_t x = hit->GetX();
      Double_t y = hit->GetY();
      Double_t z = hit->GetZ();
      Double_t de = hit->GetDe();
      Int_t cl_size = hit->GetClusterSize();
      Int_t row = hit->MeanRow();
      Double_t mrow = hit->MeanRow();
      // Double_t de_center = hit->GetDe_center();
      // TVector3 pos_center = hit->GetPos_center();
      event.cluster_hitpos_x.push_back(x);
      event.cluster_hitpos_y.push_back(y);
      event.cluster_hitpos_z.push_back(z);
      event.cluster_de.push_back(de);
      event.cluster_size.push_back(cl_size);
      event.cluster_layer.push_back(layer);
      event.cluster_row.push_back(row);
      event.cluster_mrow.push_back(mrow);
      // event.cluster_de_center.push_back(de_center);
      // event.cluster_hitpos_center_x.push_back(pos_center.X());
      // event.cluster_hitpos_center_y.push_back(pos_center.Y());
      // event.cluster_hitpos_center_z.push_back(pos_center.Z());
#if Gain_center
      //	if(69.<time&&time<85.&&nhit==1)
      if(cl_size>1){
	if(min_ycut<y&&y<max_ycut)
	  HF1(PadHid + layer*1000 + row, de_center);
      }
#endif

      ++nh_cl_Tpc;
    }
  }
  event.nh_cluster_Tpc = nh_cl_Tpc;
#if TrackSearch
  DCAna.TrackSearchTPCHelix();
#endif

  Int_t ntTpc = DCAna.GetNTracksTPCHelix();

  event.ntTpc = ntTpc;

  //std::cout<<"ntTpc:"<<ntTpc<<std::endl;
  HF1( 10, ntTpc );
  if( event.ntTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  event.nhtrack.resize( ntTpc );
  event.nhtrack_clmulti.resize( ntTpc );
  event.isBeam.resize( ntTpc );
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
  event.mom0_cor_x.resize( ntTpc );
  event.mom0_cor_y.resize( ntTpc );
  event.mom0_cor_z.resize( ntTpc );
  event.mom0_cor.resize( ntTpc );


  event.dE.resize( ntTpc );
  event.dEdx.resize( ntTpc );
  event.dEdx2.resize( ntTpc );

  event.dEdx_rms.resize( ntTpc );
  event.dEdx_rms20.resize( ntTpc );
  event.dEdx_rms10.resize( ntTpc );

  event.dEdx_20.resize( ntTpc );
  event.dEdx_10.resize( ntTpc );
  event.dEdx_geo.resize( ntTpc );
  event.dEdx_geo_20.resize( ntTpc );
  event.dEdx_geo_10.resize( ntTpc );
  event.dEdx_harm.resize( ntTpc );
  event.dEdx_harm_20.resize( ntTpc );
  event.dEdx_harm_10.resize( ntTpc );
  event.dE_clmulti.resize( ntTpc );
  event.dEdx_clmulti.resize( ntTpc );
  event.dEdx_clmulti_20.resize( ntTpc );
  event.dEdx_clmulti_10.resize( ntTpc );
  event.dE_cor.resize( ntTpc );
  event.dEdx_cor.resize( ntTpc );
  event.dEdx2_cor.resize( ntTpc );
  event.dEdx_cor_20.resize( ntTpc );
  event.dEdx_cor_10.resize( ntTpc );
  event.dz_factor.resize( ntTpc );
  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );

  // under dev
  event.combi_id.resize( ntTpc );
  event.mom_vtx.resize( ntTpc );
  event.mom_vty.resize( ntTpc );
  event.mom_vtz.resize( ntTpc );
  event.mom_cor_vtx.resize( ntTpc );
  event.mom_cor_vty.resize( ntTpc );
  event.mom_cor_vtz.resize( ntTpc );
  event.pid.resize( ntTpc );
  event.vtx.resize( ntTpc );
  event.vty.resize( ntTpc );
  event.vtz.resize( ntTpc );
  event.closeDist.resize( ntTpc );
  event.chisqr_flag.resize( ntTpc );

  event.hitlayer.resize( ntTpc );
  event.hitmrow.resize( ntTpc );
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
  event.hitde.resize( ntTpc );
  event.hitClsize.resize( ntTpc );


  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = DCAna.GetTrackTPCHelix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    Double_t helix_z0=tp->Getz0(), helix_r=tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
    Int_t isbeam = tp->GetIsBeam();

    event.nhtrack[it] = nh;
    event.isBeam[it] = isbeam;
    event.chisqr[it] = chisqr;
    event.helix_cx[it] = helix_cx;
    event.helix_cy[it] = helix_cy;
    event.helix_z0[it] = helix_z0;
    event.helix_r[it] = helix_r ;
    event.helix_dz[it] = helix_dz;
    event.mom0_x[it] = Mom0.x();;
    event.mom0_y[it] = Mom0.y();;
    event.mom0_z[it] = Mom0.z();;
    event.mom0[it] = Mom0.Mag();;

    event.hitlayer[it].resize( nh );
    event.hitmrow[it].resize( nh );
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
    event.hitde[it].resize( nh );
    event.hitClsize[it].resize( nh );

    double min_closeDist = 1000000.;
    double par1[5]={helix_cx, helix_cy, helix_z0,
		    helix_r, helix_dz};

    for( Int_t it2=0; it2<ntTpc; ++it2 ){
      if(it==it2)
	continue;
      TPCLocalTrackHelix *tp2 = DCAna.GetTrackTPCHelix( it2 );
      if( !tp2 ) continue;
      Double_t helix_cx2=tp2->Getcx(), helix_cy2=tp2->Getcy();
      Double_t helix_z02=tp2->Getz0(), helix_r2=tp2->Getr();
      Double_t helix_dz2 = tp2->Getdz();
      Double_t chisqr2 = tp2->GetChiSquare();

      double par2[5]={helix_cx2, helix_cy2, helix_z02,
		      helix_r2, helix_dz2};
      double closeDist, t1, t2;
      TVector3 vert = Kinematics::VertexPointHelix(par1, par2,
						    closeDist, t1, t2);
      if(closeDist<min_closeDist){
	event.combi_id[it] = it2;
	event.vtx[it] = vert.x();
	event.vty[it] = vert.y();
	event.vtz[it] = vert.z();
	event.closeDist[it] = closeDist;
	TVector3 mom_vtx = tp->CalcHelixMom(par1, vert.y());
	//TVector3 mom_vtx = tp->CalcHelixMom_t(par1, t1);
	event.mom_vtx[it] = mom_vtx.x();
	event.mom_vty[it] = mom_vtx.y();
	event.mom_vtz[it] = mom_vtx.z();
	event.closeDist[it] = closeDist;
	if(chisqr<chisqrCut&&chisqr2<chisqrCut)
	  event.chisqr_flag[it] = 1;
	else
	  event.chisqr_flag[it] = 0;
      }
    }

    Int_t nh_clmulti =0;
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      int clsize_hit = hit->GetHit()->GetClusterSize();
      if(clsize_hit>1)
	++nh_clmulti;
    }

    double min_t = 10000.;
    double max_t = -10000.;
    Int_t min_layer=-1, max_layer=33;
    double min_layer_t=0., max_layer_t=0.;
    // double max_layer_y=0.;
    double de=0., path_dEdx=0.;
    double de_clmulti=0., path_clmulti_dEdx=0.;
    double path_dEdx_cor=0.;

    std::vector<Double_t> dEdx_vect;
    std::vector<Double_t> dEdx_cor_vect;
    std::vector<Double_t> dEdx_clmulti_vect;


    int ih_clmulti=0;
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& res_vect = hit->GetResidualVect();
      double de_hit = hit->GetHit()->GetCDe();
      int clsize_hit = hit->GetHit()->GetClusterSize();
      double path_hit = tpc::padParameter[layer][5];

      double mrow = hit->GetHit()->GetMRow();
      double pad_theta = tpc::getTheta(layer, mrow)*acos(-1)/180.;
      double t_calc = hit->GetTcal();
      double theta_diff = t_calc - pad_theta;
      //Approximation of path_hit correction
      double path_hit_cor = (path_hit/fabs(cos(theta_diff)))*sqrt(1.+(pow(helix_dz,2)));
      // std::cout<<"cor cos(thete_diff)="<<cos(theta_diff)<<", "
      // 	       <<"pad_theta="<<pad_theta*180./acos(-1)<<", "
      // 	       <<"t_calc="<<t_calc*180./acos(-1)<<", "
      // 	       <<"dz_fact="<<sqrt(1.+(pow(helix_dz,2)))<<std::endl;
      // std::cout<<"path_hit="<<path_hit<<", path_hitcor="<<path_hit_cor<<std::endl;

      //double path_hit_cor=0.;

#if DE_padparam
      path_hit = 1.;
#endif
      de += de_hit;
      path_dEdx += path_hit;
      path_dEdx_cor += path_hit_cor;

      double dEdx = de_hit/path_hit;
      double dEdx_cor = de_hit/path_hit_cor;
      dEdx_vect.push_back(dEdx);
      dEdx_cor_vect.push_back(dEdx_cor);

      if(clsize_hit>1){
       	de_clmulti += de_hit;
	path_clmulti_dEdx += path_hit;
	dEdx_clmulti_vect.push_back(de_hit/path_hit);
	++ih_clmulti;
      }

      Double_t t_cal = hit->GetTcal();
      if(min_t>t_cal)
	min_t = t_cal;
      if(max_t<t_cal)
	max_t = t_cal;
      if(layer<min_layer){
	min_layer = layer;
	min_layer_t = t_cal;
      }
      if(layer>max_layer){
	max_layer = layer;
	max_layer_t = t_cal;
      }

      Double_t residual = hit->GetResidual();

      // for(int ii=0; ii<ih; ++ii){
      // 	if((int)event.hitlayer[it][ii]==layer){
      // 	  std::cout<<"layer="<<layer<<std::endl;
      // 	  std::cout<<"mrow1="<<event.hitmrow[it][ii]<<
      // 	    ", "<<mrow<<std::endl;
      // 	  std::cout<<"hitpos1=("<<event.hitpos_x[it][ii]<<", "
      // 		   <<event.hitpos_y[it][ii]<<", "
      // 		   <<event.hitpos_z[it][ii]<<"), hitpos2="
      // 		   <<hitpos<<std::endl;
      // 	}
      // }
      event.hitlayer[it][ih] = (double)layer;
      event.hitmrow[it][ih] = mrow;
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
      event.helix_t[it][ih] = t_cal;
      event.hitde[it][ih] = de_hit;
      event.hitClsize[it][ih] = (double)clsize_hit;
    }
    if(min_layer_t<max_layer_t)
      event.charge[it] = 1;
    else
      event.charge[it] = -1;

    TVector3 Mom0_cor, mom_cor_vtx;
    if(event.charge[it] == 1){
      Mom0_cor = tp->GetMom0_corP();
      mom_cor_vtx = tp->CalcHelixMom_corP(par1, event.vty[it]);
    }
    if(event.charge[it] == -1){
      Mom0_cor = tp->GetMom0_corN();
      mom_cor_vtx = tp->CalcHelixMom_corN(par1, event.vty[it]);
    }

    event.mom0_cor_x[it]=Mom0_cor.x();
    event.mom0_cor_y[it]=Mom0_cor.y();
    event.mom0_cor_z[it]=Mom0_cor.z();
    event.mom0_cor[it]=Mom0_cor.Mag();

    event.mom_cor_vtx[it]=mom_cor_vtx.x();
    event.mom_cor_vty[it]=mom_cor_vtx.y();
    event.mom_cor_vtz[it]=mom_cor_vtx.z();


    Double_t pathlen = (max_t - min_t)*sqrt(helix_r*helix_r*(1.+helix_dz*helix_dz));
    // Int_t htofseg = tp->GetHTOFSeg(min_layer_t, max_layer_t, max_layer_y);
    // std::cout<<"min_t="<<min_t<<", max_t="<<max_t<<", helix_r="<<helix_r
    //          <<", path="<<pathlen<<std::endl;
    event.path[it] = pathlen;
    event.dE[it] = de;
    event.dEdx[it] = de/path_dEdx;
    event.dE_cor[it] = 0.;
    event.dEdx_cor[it] = de/path_dEdx_cor;

    std::sort(dEdx_vect.begin(), dEdx_vect.end());
    std::sort(dEdx_cor_vect.begin(), dEdx_cor_vect.end());
    std::sort(dEdx_clmulti_vect.begin(), dEdx_clmulti_vect.end());

    event.dEdx2[it] = TMath::Mean(dEdx_vect.size(),dEdx_vect.data());
    event.dEdx2_cor[it] = TMath::Mean(dEdx_cor_vect.size(),dEdx_cor_vect.data());
    event.dE_clmulti[it] = de_clmulti;
    event.dEdx_clmulti[it] = de_clmulti/path_clmulti_dEdx;

    event.dEdx_geo[it]=1.;
    event.dEdx_geo_20[it]=1.;
    event.dEdx_geo_10[it]=1.;
    event.dEdx_harm[it]=0.;
    event.dEdx_harm_20[it]=0.;
    event.dEdx_harm_10[it]=0.;

    event.dEdx_20[it]=0.;
    event.dEdx_cor_20[it]=0.;
    event.dEdx_10[it]=0.;
    event.dEdx_cor_10[it]=0.;

    double dEdxRMS= TMath::RMS(dEdx_vect.size(),dEdx_vect.data());
    event.dEdx_rms[it]=dEdxRMS;
    int n_rms10=0;
    int n_rms20=0;
    event.dEdx_rms10[it]=0.;
    event.dEdx_rms20[it]=0.;

    int n_20 = (int)(dEdx_vect.size()*0.8);
    int n_10 = (int)(dEdx_vect.size()*0.9);
    for( int ih=0; ih<dEdx_vect.size(); ++ih ){
      event.dEdx_geo[it] *= dEdx_vect[it];
      event.dEdx_harm[it] += 1./dEdx_vect[it];
      if(ih<n_20){
	event.dEdx_20[it]+=dEdx_vect[ih];
	event.dEdx_cor_20[it] +=dEdx_cor_vect[ih];
	event.dEdx_geo_20[it] *=dEdx_vect[ih];
	event.dEdx_harm_20[it] += 1./dEdx_vect[it];
      }
      if(ih<n_10){
	event.dEdx_10[it]+=dEdx_vect[ih];
	event.dEdx_cor_10[it]+=dEdx_cor_vect[ih];
	event.dEdx_geo_10[it] *=dEdx_vect[ih];
	event.dEdx_harm_10[it] += 1./dEdx_vect[it];
      }
      //80% in gaus
      // if(event.dEdx2[it]-1.28*dEdxRMS<dEdx_vect[ih]
      // 	 &&dEdx_vect[ih]<event.dEdx2[it]+1.28*dEdxRMS){
      if(dEdx_vect[ih]<event.dEdx2[it]+1.28*dEdxRMS){
	event.dEdx_rms20[it]+=dEdx_vect[ih];
	++n_rms20;
      }
      //90% in gaus
      // if(event.dEdx2[it]-1.64*dEdxRMS<dEdx_vect[ih]
      // 	 &&dEdx_vect[ih]<event.dEdx2[it]+1.64*dEdxRMS){
      if(dEdx_vect[ih]<event.dEdx2[it]+1.64*dEdxRMS){
	event.dEdx_rms10[it]+=dEdx_vect[ih];
	++n_rms10;
      }

      // std::cout<<"inside for, ih="<<ih<<", geo="
      // 	       <<event.dEdx_geo[it]
      // 	       <<", harm="<<event.dEdx_harm[it]
      // 	       <<", harm_20="<<event.dEdx_harm_20[it]<<std::endl;
    }

    event.dEdx_20[it]/=(double)n_20;
    event.dEdx_cor_20[it]/=(double)n_20;
    event.dEdx_10[it]/=(double)n_10;
    event.dEdx_cor_10[it]/=(double)n_10;

    event.dEdx_geo[it] = pow(event.dEdx_geo[it], 1./(double)nh);
    event.dEdx_geo_20[it] = pow(event.dEdx_geo_20[it], 1./(double)n_20);
    event.dEdx_geo_10[it] = pow(event.dEdx_geo_10[it], 1./(double)n_10);
    event.dEdx_harm[it] = ((double)nh)/event.dEdx_harm[it];
    event.dEdx_harm_20[it] = ((double)n_20)/event.dEdx_harm_20[it];
    event.dEdx_harm_10[it] = ((double)n_10)/event.dEdx_harm_10[it];

    event.dEdx_rms20[it]/=(double)n_rms20;
    event.dEdx_rms10[it]/=(double)n_rms10;
    // std::cout<<"geo="<<event.dEdx_geo[it]
    // 	     <<", harm="<<event.dEdx_harm[it]
    // 	     <<", harm_20="<<event.dEdx_harm_20[it]<<std::endl;

    // std::cout<<"dEdx="<<event.dEdx[it]<<", "
    // 	     <<"dEdx2="<<event.dEdx2[it]<<", "
    // 	     <<"dEdx_cor="<<event.dEdx_cor[it]<<", "
    // 	     <<"dEdx_20="<<event.dEdx_20[it]<<", "
    // 	     <<"dEdx_cor_20="<<event.dEdx_cor_20[it]<<std::endl;

    event.dEdx_clmulti_20[it]=0.;
    double n_clmulti_20 = dEdx_clmulti_vect.size()*0.8;
    double n_clmulti_10 = dEdx_clmulti_vect.size()*0.9;
    if(n_clmulti_20>1.){
      for( int ih=0; ih<(int)n_clmulti_20; ++ih ){
     	event.dEdx_clmulti_20[it]+=dEdx_clmulti_vect[ih];
      }
      event.dEdx_clmulti_20[it]/=(double)n_clmulti_20;
    }
    event.dEdx_clmulti_10[it]=0.;
    if(n_clmulti_10>1.){
      for( int ih=0; ih<(int)n_clmulti_10; ++ih ){
     	event.dEdx_clmulti_10[it]+=dEdx_clmulti_vect[ih];
      }
      event.dEdx_clmulti_10[it]/=(double)n_clmulti_10;
    }

    //double dz_factor = sqrt(1.+(pow(helix_dz,2)));
    event.dz_factor[it] = sqrt(1.+(pow(helix_dz,2)));
    // event.dE_cor[it]=event.dE[it]/dz_factor;
    // event.dEdx_cor[it]=event.dEdx[it]/dz_factor;
    // event.dEdx_20_cor[it]=event.dEdx_20[it]/dz_factor;
    // event.dEdx_10_cor[it]=event.dEdx_10[it]/dz_factor;
    int particle= Kinematics::PID_HypTPC_dEdx(event.dEdx_rms20[it],
					      event.mom0_cor[it],
					      event.charge[it]);
    event.pid[it]=particle;
  }


  // rough estimation for Lambda and Ks event
  if(ntTpc>1){
    for( Int_t it=0; it<ntTpc-1; ++it ){
      for( Int_t it2=it+1; it2<ntTpc; ++it2 ){
	if(event.isBeam[it]==0&&event.isBeam[it2]==0
	   &&event.combi_id[it]==it2&&event.combi_id[it2]==it
	   &&event.charge[it]==-1*event.charge[it2]
	   &&event.chisqr_flag[it]==1&&event.chisqr_flag[it2]==1){
	  if(event.charge[it]==1){
	    TVector3 mom_pos(event.mom_cor_vtx[it], event.mom_cor_vty[it], event.mom_cor_vtz[it]);
	    TLorentzVector Lp(mom_pos, std::sqrt(ProtonMass*ProtonMass+mom_pos.Mag2()));
	    TLorentzVector Lpip(mom_pos, std::sqrt(PionMass*PionMass+mom_pos.Mag2()));
	    TVector3 mom_neg(event.mom_cor_vtx[it2], event.mom_cor_vty[it2], event.mom_cor_vtz[it2]);
	    TLorentzVector Lpi(mom_neg, std::sqrt(PionMass*PionMass+mom_neg.Mag2()));
	    TLorentzVector LLambda = Lp + Lpi;
	    TLorentzVector LKs = Lpip + Lpi;
	    //	    if(event.dEdx_20_cor[it]>event.dEdx_20_cor[it2]*1.5){
	    if(event.pid[it]==2&&event.pid[it2]==-1){
	      event.M_Lambda.push_back(LLambda.M());
	      event.Lvtx.push_back(event.vtx[it]);
	      event.Lvty.push_back(event.vty[it]);
	      event.Lvtz.push_back(event.vtz[it]);
	      event.LcloseDist.push_back(event.closeDist[it]);
	      event.Mom_Lambda.push_back(LLambda.P());
	      event.Mom_Lambdax.push_back(LLambda.Px());
	      event.Mom_Lambday.push_back(LLambda.Py());
	      event.Mom_Lambdaz.push_back(LLambda.Pz());
	    }
	    if(event.pid[it]==1&&event.pid[it2]==-1){
	      event.M_Ks.push_back(LKs.M());
	      event.Ksvtx.push_back(event.vtx[it]);
	      event.Ksvty.push_back(event.vty[it]);
	      event.Ksvtz.push_back(event.vtz[it]);
	      event.KscloseDist.push_back(event.closeDist[it]);
	      event.Mom_Ks.push_back(LKs.P());
	      event.Mom_Ksx.push_back(LKs.Px());
	      event.Mom_Ksx.push_back(LKs.Py());
	      event.Mom_Ksx.push_back(LKs.Pz());
	    }
	  }
	  else{
	    TVector3 mom_pos(event.mom_cor_vtx[it2], event.mom_cor_vty[it2], event.mom_cor_vtz[it2]);
	    TLorentzVector Lp(mom_pos, std::sqrt(ProtonMass*ProtonMass+mom_pos.Mag2()));
	    TLorentzVector Lpip(mom_pos, std::sqrt(PionMass*PionMass+mom_pos.Mag2()));
	    TVector3 mom_neg(event.mom_cor_vtx[it], event.mom_cor_vty[it], event.mom_cor_vtz[it]);
	    TLorentzVector Lpi(mom_neg, std::sqrt(PionMass*PionMass+mom_neg.Mag2()));
	    TLorentzVector LLambda = Lp + Lpi;
	    TLorentzVector LKs = Lpip + Lpi;
	    //	    if(event.dEdx_20_cor[it2]>event.dEdx_20_cor[it]*1.5){
	    if(event.pid[it2]==2&&event.pid[it]==-1){
	      event.M_Lambda.push_back(LLambda.M());
	      event.Lvtx.push_back(event.mom_vtx[it]);
	      event.Lvty.push_back(event.mom_vty[it]);
	      event.Lvtz.push_back(event.mom_vtz[it]);
	      event.LcloseDist.push_back(event.closeDist[it]);
	      event.Mom_Lambda.push_back(LLambda.P());
	      event.Mom_Lambdax.push_back(LLambda.Px());
	      event.Mom_Lambday.push_back(LLambda.Py());
	      event.Mom_Lambdaz.push_back(LLambda.Pz());
	    }
	    if(event.pid[it2]==1&&event.pid[it]==-1){
	      event.M_Ks.push_back(LKs.M());
	      event.Ksvtx.push_back(event.mom_vtx[it]);
	      event.Ksvty.push_back(event.mom_vty[it]);
	      event.Ksvtz.push_back(event.mom_vtz[it]);
	      event.KscloseDist.push_back(event.closeDist[it]);
	      event.Mom_Ks.push_back(LKs.P());
	      event.Mom_Ksx.push_back(LKs.Px());
	      event.Mom_Ksy.push_back(LKs.Py());
	      event.Mom_Ksz.push_back(LKs.Pz());
	    }
	  }
	}
      }
    }
  }

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
  // const Int_t    NbinDe = 1000;
  // const Double_t MinDe  =    0.;
  // const Double_t MaxDe  = 2000.;
  HB1( 1, "Status", 21, 0., 21. );
  HB1( 10, "NTrack TPC", 40, 0., 40. );

  HBTree( "tpc", "tree of DstTPCTracking" );

  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "nh_cluster_Tpc", &event.nh_cluster_Tpc );
  tree->Branch( "raw_hitpos_x", &event.raw_hitpos_x );
  tree->Branch( "raw_hitpos_y", &event.raw_hitpos_y );
  tree->Branch( "raw_hitpos_z", &event.raw_hitpos_z );
  tree->Branch( "raw_de", &event.raw_de );
  tree->Branch( "raw_padid", &event.raw_padid );
  tree->Branch( "raw_layer", &event.raw_layer );
  tree->Branch( "raw_row", &event.raw_row );
  tree->Branch( "cluster_hitpos_x", &event.cluster_hitpos_x );
  tree->Branch( "cluster_hitpos_y", &event.cluster_hitpos_y );
  tree->Branch( "cluster_hitpos_z", &event.cluster_hitpos_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_size", &event.cluster_size );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_row", &event.cluster_row );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_de_center", &event.cluster_de_center );
  tree->Branch( "cluster_hitpos_center_x", &event.cluster_hitpos_center_x );
  tree->Branch( "cluster_hitpos_center_y", &event.cluster_hitpos_center_y );
  tree->Branch( "cluster_hitpos_center_z", &event.cluster_hitpos_center_z );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "isBeam", &event.isBeam );
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
  tree->Branch( "mom0_cor_x", &event.mom0_cor_x );
  tree->Branch( "mom0_cor_y", &event.mom0_cor_y );
  tree->Branch( "mom0_cor_z", &event.mom0_cor_z );
  tree->Branch( "mom0_cor", &event.mom0_cor );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "dEdx2", &event.dEdx2 );
  tree->Branch( "dEdx2_cor", &event.dEdx2_cor );
  tree->Branch( "dEdx_20", &event.dEdx_20 );
  tree->Branch( "dEdx_10", &event.dEdx_10 );

  tree->Branch( "dEdx_geo", &event.dEdx_geo );
  tree->Branch( "dEdx_geo_20", &event.dEdx_geo_20 );
  tree->Branch( "dEdx_geo_10", &event.dEdx_geo_10 );
  tree->Branch( "dEdx_harm", &event.dEdx_harm );
  tree->Branch( "dEdx_harm_20", &event.dEdx_harm_20 );
  tree->Branch( "dEdx_harm_10", &event.dEdx_harm_10 );

  tree->Branch( "dEdx_rms", &event.dEdx_rms );
  tree->Branch( "dEdx_rms20", &event.dEdx_rms20 );
  tree->Branch( "dEdx_rms10", &event.dEdx_rms10 );

  tree->Branch( "dE_clmulti", &event.dE_clmulti );
  tree->Branch( "dEdx_clmulti", &event.dEdx_clmulti );
  tree->Branch( "dEdx_clmulti_20", &event.dEdx_clmulti_20 );
  tree->Branch( "dEdx_clmulti_10", &event.dEdx_clmulti_10 );
  tree->Branch( "dE_cor", &event.dE_cor );
  tree->Branch( "dEdx_cor", &event.dEdx_cor );
  tree->Branch( "dEdx_cor_20", &event.dEdx_cor_20 );
  tree->Branch( "dEdx_cor_10", &event.dEdx_cor_10 );
  tree->Branch( "dz_factor", &event.dz_factor );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );

  tree->Branch( "pid", &event.pid );
  tree->Branch( "combi_id", &event.combi_id );
  tree->Branch( "chisqr_flag", &event.chisqr_flag );
  tree->Branch( "mom_vtx", &event.mom_vtx );
  tree->Branch( "mom_vty", &event.mom_vty );
  tree->Branch( "mom_vtz", &event.mom_vtz );
  tree->Branch( "mom_cor_vtx", &event.mom_cor_vtx );
  tree->Branch( "mom_cor_vty", &event.mom_cor_vty );
  tree->Branch( "mom_cor_vtz", &event.mom_cor_vtz );
  tree->Branch( "vtx", &event.vtx );
  tree->Branch( "vty", &event.vty );
  tree->Branch( "vtz", &event.vtz );
  tree->Branch( "closeDist", &event.closeDist );
  tree->Branch( "M_Lambda", &event.M_Lambda );
  tree->Branch( "Lvtx", &event.Lvtx );
  tree->Branch( "Lvty", &event.Lvty );
  tree->Branch( "Lvtz", &event.Lvtz );
  tree->Branch( "LcloseDist", &event.LcloseDist );
  tree->Branch( "Mom_Lambda", &event.Mom_Lambda );
  tree->Branch( "Mom_Lambdax", &event.Mom_Lambdax );
  tree->Branch( "Mom_Lambday", &event.Mom_Lambday );
  tree->Branch( "Mom_Lambdaz", &event.Mom_Lambdaz );
  tree->Branch( "M_Ks", &event.M_Ks );
  tree->Branch( "Ksvtx", &event.Ksvtx );
  tree->Branch( "Ksvty", &event.Ksvty );
  tree->Branch( "Ksvtz", &event.Ksvtz );
  tree->Branch( "KscloseDist", &event.KscloseDist );
  tree->Branch( "Mom_Ks", &event.Mom_Ks );
  tree->Branch( "Mom_Ksx", &event.Mom_Ksx );
  tree->Branch( "Mom_Ksy", &event.Mom_Ksy );
  tree->Branch( "Mom_Ksz", &event.Mom_Ksz );

  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitmrow", &event.hitmrow );
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
  tree->Branch( "hitde", &event.hitde );
  tree->Branch( "hitClsize", &event.hitClsize );

#if Gain_center
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for( Int_t r=0; r<NumOfRow; ++r ){
      HB1(PadHid + layer*1000 + r , "TPC DeltaE_center", NbinDe, MinDe, MaxDe );
    }
  }
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
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
