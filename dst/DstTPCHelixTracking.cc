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
#include "TPCLocalTrack_Helix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define TrackSearch 1
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
  Int_t nhTpc;
  Int_t nh_cluster_Tpc;
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

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0

  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path

  // not yet prepared
  std::vector<Double_t> mom_vtx;//Helix momentum at vtx
  std::vector<Double_t> mom_vty;//Helix momentum at vtx
  std::vector<Double_t> mom_vtz;//Helix momentum at vtx
  std::vector<Double_t> vtx;
  std::vector<Double_t> vty;
  std::vector<Double_t> vtz;
  std::vector<Double_t> closeDist;



  
  std::vector<std::vector<Int_t>> hitlayer;
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
    cluster_hitpos_x.clear();
    cluster_hitpos_y.clear();
    cluster_hitpos_z.clear();
    cluster_de.clear();
    cluster_size.clear();
    cluster_layer.clear();
    cluster_row.clear();
    cluster_mrow.clear();
    cluster_de_center.clear();
    ntTpc = 0;
    trigpat.clear();
    trigflag.clear();
    nhtrack.clear();
    chisqr.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    mom0_x.clear();
    mom0_y.clear();
    mom0_z.clear();
    mom0.clear();

    charge.clear();
    path.clear();

    mom_vtx.clear();
    mom_vty.clear();
    mom_vtz.clear();
    vtx.clear();
    vty.clear();
    vtz.clear();
    closeDist.clear();

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
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

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
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc);
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
      Double_t de = hit->GetCharge();
      Double_t pad = hit->GetPad();
      event.raw_hitpos_x.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);
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
      Double_t de = hit->GetCharge();
      Int_t cl_size = hit->GetClusterSize();
      Int_t row = hit->GetRow();
      Double_t mrow = hit->GetMRow();
      Double_t de_center = hit->GetCharge_center();
      event.cluster_hitpos_x.push_back(x);
      event.cluster_hitpos_y.push_back(y);
      event.cluster_hitpos_z.push_back(z);
      event.cluster_de.push_back(de);
      event.cluster_size.push_back(cl_size);
      event.cluster_layer.push_back(layer);
      event.cluster_row.push_back(row);
      event.cluster_mrow.push_back(mrow);
      event.cluster_de_center.push_back(de_center);
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
  DCAna.TrackSearchTPC_Helix();
#endif 

  Int_t ntTpc = DCAna.GetNTracksTPC_Helix();
  event.ntTpc = ntTpc;
  HF1( 10, ntTpc );
  if( event.ntTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  event.nhtrack.resize( ntTpc );
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
  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );
  
  // not yet prepared
  event.mom_vtx.resize( ntTpc );
  event.mom_vty.resize( ntTpc );
  event.mom_vtz.resize( ntTpc );
  event.vtx.resize( ntTpc );
  event.vty.resize( ntTpc );
  event.vtz.resize( ntTpc );
  event.closeDist.resize( ntTpc );

  
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
  
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrack_Helix *tp = DCAna.GetTrackTPC_Helix( it );
    if( !tp ) continue;
    Int_t nh = tp->GetNHit();
    Double_t chisqr = tp->GetChiSquare();
    Double_t helix_cx=tp->Getcx(), helix_cy=tp->Getcy();
    Double_t helix_z0=tp->Getz0(), helix_r=tp->Getr();
    Double_t helix_dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();
    event.nhtrack[it] = nh;
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

    double min_t = 10000.;
    double max_t = -10000.;
    double min_layer_t=0., max_layer_t=0.;
    
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPos_Helix();
      const TVector3& res_vect = hit->GetResidualVect();
      Double_t t_cal = hit->GetTcal();
      if(min_t>t_cal)
	min_t = t_cal;
      if(max_t<t_cal)
	max_t = t_cal;
      if(ih==0)
	min_layer_t = t_cal;
      if(ih==nh-1)
	max_layer_t = t_cal;
      Double_t residual = hit->GetResidual();
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
      event.helix_t[it][ih] = t_cal;
    }
    if(min_layer_t<max_layer_t)
      event.charge[it] = 1;
    else
      event.charge[it] = -1;
    double pathlen = (max_t - min_t)*sqrt(helix_r*helix_r + helix_dz*helix_r);
    //std::cout<<"min_t="<<min_t<<", max_t="<<max_t<<", helix_r="<<helix_r<<", path="<<pathlen<<std::endl;
    event.path[it] = pathlen;
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
  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 2000.;
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
  tree->Branch( "cluster_hitpos_x", &event.cluster_hitpos_x );
  tree->Branch( "cluster_hitpos_y", &event.cluster_hitpos_y );
  tree->Branch( "cluster_hitpos_z", &event.cluster_hitpos_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_size", &event.cluster_size );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_row", &event.cluster_row );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_de_center", &event.cluster_de_center );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
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
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );
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
  src.tTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tTpc" );
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
