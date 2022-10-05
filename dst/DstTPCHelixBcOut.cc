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
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#define TrackSearch 1
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
const Int_t MinPosMapXZ = -270;
const Int_t MaxPosMapXZ = 270;
const Int_t MinPosMapY = -20;
const Int_t MaxPosMapY = 20;
const Int_t Meshsize = 10;

const Double_t& zK18HS = gGeom.LocalZ("K18HS");
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit,  kBcOut, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]", "[BcOut]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc","bcout" ,"" };
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
  std::vector<Double_t> cluster_hitpos_center_x;
  std::vector<Double_t> cluster_hitpos_center_y;
  std::vector<Double_t> cluster_hitpos_center_z;


  //BcOut info
  Int_t ntrack_bcout;
  std::vector<Double_t> chisqr_bcout;
  std::vector<Double_t> x0_bcout;
  std::vector<Double_t> y0_bcout;
  std::vector<Double_t> u0_bcout;
  std::vector<Double_t> v0_bcout;

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

  std::vector<Double_t> mom0_x;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_y;//Helix momentum at Y = 0
  std::vector<Double_t> mom0_z;//Helix momentum at Y = 0
  std::vector<Double_t> mom0;//Helix momentum at Y = 0

  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path

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

    ntrack_bcout = 0;
    chisqr_bcout.clear();
    x0_bcout.clear();
    y0_bcout.clear();
    u0_bcout.clear();
    v0_bcout.clear();

    nhtrack.clear();
    nhtrack_clmulti.clear();
    isBeam.clear();
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
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* ctTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting

  //BcOut input
  Int_t ntrack;
  Double_t chisqr[MaxHits];
  Double_t x0[MaxHits];
  Double_t y0[MaxHits];
  Double_t u0[MaxHits];
  Double_t v0[MaxHits];
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    PosXHid    = 1000000,
    PosYHid    = 2000000,
    PosYPadHid    = 3000000,
    HitdEHid    = 100000,
    ClCenterdEHid    = 200000,
  };
}

Int_t GetHistNum(Double_t x, Double_t y, Double_t z){
  int ix = (int)((x - (MinPosMapXZ - Meshsize/2.)))/Meshsize;
  int iy = (int)((y - (MinPosMapY - Meshsize/2.)))/Meshsize;
  int iz = (int)((z - (MinPosMapXZ - Meshsize/2.)))/Meshsize;

  Int_t NumOfDivXZ = ((MaxPosMapXZ - MinPosMapXZ)/Meshsize) + 1;
  Int_t NumOfDivY = ((MaxPosMapY - MinPosMapY)/Meshsize) + 1;
  int histnum = ix*NumOfDivY*NumOfDivXZ + iy*NumOfDivXZ + iz;
  return histnum;
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

  event.ntrack_bcout = src.ntrack;
  //  event.nhTpc = **src.nhTpc;

  // event.chisqr_bcout.resize( event.ntrack_bcout );
  // event.x0_bcout.resize( event.ntrack_bcout );
  // event.y0_bcout.resize( event.ntrack_bcout );
  // event.u0_bcout.resize( event.ntrack_bcout );
  // event.v0_bcout.resize( event.ntrack_bcout );

  for(int it=0; it<src.ntrack; ++it){
    event.chisqr_bcout.push_back(src.chisqr[it]);
    event.x0_bcout.push_back(src.x0[it]);
    event.y0_bcout.push_back(src.y0[it]);
    event.u0_bcout.push_back(src.u0[it]);
    event.v0_bcout.push_back(src.v0[it]);
  }


  HF1( 1, event.status++ );

  if( **src.nhTpc == 0 )
    return true;

  HF1( 1, event.status++ );

  DCAnalyzer DCAna;
  //  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc);
  DCAna.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.ctTpc, **src.deTpc);
#if HoughYcut
  DCAna.HoughYCut(min_ycut, max_ycut);
#endif

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
      Double_t pad = hit->GetPad();
      Int_t row = hit->GetRow();
      event.raw_hitpos_x.push_back(x);
      event.raw_hitpos_y.push_back(y);
      event.raw_hitpos_z.push_back(z);
      event.raw_de.push_back(de);
      event.raw_padid.push_back(pad);

      for(int it=0; it<src.ntrack; ++it){
	Double_t x0_BC = src.x0[it];
     	Double_t u0_BC = src.u0[it];
     	Double_t y0_BC = src.y0[it];
     	Double_t v0_BC = src.v0[it];

     	Double_t zTPC = zK18HS + z;
     	Double_t x_BC = x0_BC + zTPC*u0_BC;
     	Double_t y_BC = y0_BC + zTPC*v0_BC;
     	if(fabs(x_BC - x)<100.&&de>100.){
     	  HF1(PosYPadHid + layer*1000+ row,  y - y_BC);
     	}
     	if(fabs(x_BC - x)<100.&& fabs(y_BC - y)< 60){
     	  HF1(HitdEHid + layer*1000+ row,  de);
     	}
     }

      ++nh_Tpc;
    }
  }
  event.nhTpc = nh_Tpc;

  HF1( 1, event.status++ );
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

      //Compared with BcOut
      for(int it=0; it<src.ntrack; ++it){
	Double_t x0_BC = src.x0[it];
	Double_t u0_BC = src.u0[it];
	Double_t y0_BC = src.y0[it];
	Double_t v0_BC = src.v0[it];

	Double_t zTPC = zK18HS + z;
	Double_t x_BC = x0_BC + zTPC*u0_BC;
	Double_t y_BC = y0_BC + zTPC*v0_BC;
	if(MinPosMapXZ<x && x<MaxPosMapXZ&&
	   MinPosMapY<y && y<MaxPosMapY&&
	   MinPosMapXZ<z && z<MaxPosMapXZ&&
	   cl_size>=2){
	  int histnum = GetHistNum(x, y, z);
	  if(fabs(y_BC - y)<60.)
	    HF1(PosXHid + histnum, x_BC - x);
	  if(fabs(x_BC - x)<100.){
	    HF1(PosYHid + histnum, y_BC - y);
	    //	    HF1(PosYPadHid + layer*1000+ row,  pos_center.Y() - y_BC);
	  }
	  if(fabs(y_BC - y)<60.&&fabs(x_BC - x)<100.){
	    HF1(100 + layer, de);
	    // HF1(ClCenterdEHid + layer*1000+ row,  de_center);
	  }
	}
      }

      ++nh_cl_Tpc;
    }
  }
  event.nh_cluster_Tpc = nh_cl_Tpc;

  HF1( 1, event.status++ );





#if TrackSearch
  DCAna.TrackSearchTPCHelix();
#endif

  Int_t ntTpc = DCAna.GetNTracksTPCHelix();
  event.ntTpc = ntTpc;
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

  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );

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

    Int_t nh_clmulti =0;
    double min_t = 10000.;
    double max_t = -10000.;
    double min_layer_t=0., max_layer_t=0.;
    Int_t min_layer=-1, max_layer=33;
    // double max_layer_y=0.;

    // int ih_clmulti=0;
    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit = tp->GetHit( ih );
      if( !hit ) continue;
      Int_t layer = hit->GetLayer();
      const TVector3& hitpos = hit->GetLocalHitPos();
      const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& res_vect = hit->GetResidualVect();
      double de_hit = hit->GetHit()->GetDe();
      int clsize_hit = hit->GetHit()->GetClusterSize();
      if(clsize_hit>1)
	++nh_clmulti;
      // double path_hit = tpc::padParameter[layer][5];

      double mrow = hit->GetHit()->GetMRow();
      // double pad_theta = tpc::getTheta(layer, mrow)*acos(-1)/180.;
      // double theta_diff = t_calc - pad_theta;

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

    Double_t pathlen = (max_t - min_t)*sqrt(helix_r*helix_r*(1.+helix_dz*helix_dz));
    event.path[it] = pathlen;
  }
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
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    HB1( 100 + layer, "dE TPC", NbinDe, MinDe, MaxDe );
  }
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for( Int_t r=0; r<NumOfRow; ++r ){
      HB1(HitdEHid + layer*1000 + r , "TPC hit dE", NbinDe, MinDe, MaxDe );
      HB1(ClCenterdEHid + layer*1000 + r , "TPC dE_center", NbinDe, MinDe, MaxDe );
    }
  }
  const Int_t    NbinPos = 1600;
  const Double_t MinPos  = -40.;
  const Double_t MaxPos  = 40.;

  Int_t NumOfDivXZ = ((MaxPosMapXZ - MinPosMapXZ)/Meshsize) + 1;
  Int_t NumOfDivY = ((MaxPosMapY - MinPosMapY)/Meshsize) + 1;

  int histnum =0;
  for(Int_t ix=0; ix<NumOfDivXZ; ++ix){
    for(Int_t iy=0; iy<NumOfDivY; ++iy){
      for(Int_t iz=0; iz<NumOfDivXZ; ++iz){
	HB1(PosXHid + histnum, "TPC Pos XCor", NbinPos, MinPos, MaxPos );
	HB1(PosYHid + histnum, "TPC Pos YCor", NbinPos, MinPos, MaxPos );
	++histnum;
      }
    }
  }
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for( Int_t r=0; r<NumOfRow; ++r ){
      HB1(PosYPadHid + layer*1000 + r , "TPC Pad Y Cor", NbinPos, MinPos, MaxPos );
    }
  }


  // for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
  //   const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
  //   for( Int_t r=0; r<NumOfRow; ++r ){
  //     HB1(PadHid + layer*1000 + r , "TPC DeltaE_center", NbinDe, MinDe, MaxDe );
  //   }
  // }


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
  tree->Branch( "cluster_hitpos_center_x", &event.cluster_hitpos_center_x );
  tree->Branch( "cluster_hitpos_center_y", &event.cluster_hitpos_center_y );
  tree->Branch( "cluster_hitpos_center_z", &event.cluster_hitpos_center_z );

  tree->Branch( "ntrack_bcout", &event.ntrack_bcout );
  tree->Branch( "chisqr_bcout", &event.chisqr_bcout );
  tree->Branch( "x0_bcout", &event.x0_bcout );
  tree->Branch( "y0_bcout", &event.y0_bcout );
  tree->Branch( "u0_bcout", &event.u0_bcout );
  tree->Branch( "v0_bcout", &event.v0_bcout );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "nhtrack_clmulti", &event.nhtrack_clmulti );
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

  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );

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

  // for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
  //   const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
  //   for( Int_t r=0; r<NumOfRow; ++r ){
  //     HB1(PadHid + layer*1000 + r , "TPC DeltaE_center", NbinDe, MinDe, MaxDe );
  //   }
  // }


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
  src.ctTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ctTpc" );
  src.chisqrTpc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTpc" );

  TTreeCont[kBcOut]->SetBranchStatus("*", 0);
  TTreeCont[kBcOut]->SetBranchStatus("ntrack",  1);
  TTreeCont[kBcOut]->SetBranchStatus("chisqr",  1);
  TTreeCont[kBcOut]->SetBranchStatus("x0",  1);
  TTreeCont[kBcOut]->SetBranchStatus("y0",  1);
  TTreeCont[kBcOut]->SetBranchStatus("u0",  1);
  TTreeCont[kBcOut]->SetBranchStatus("v0",  1);

  TTreeCont[kBcOut]->SetBranchAddress("ntrack",  &src.ntrack);
  TTreeCont[kBcOut]->SetBranchAddress("chisqr",  src.chisqr);
  TTreeCont[kBcOut]->SetBranchAddress("x0",  src.x0);
  TTreeCont[kBcOut]->SetBranchAddress("y0",  src.y0);
  TTreeCont[kBcOut]->SetBranchAddress("u0",  src.u0);
  TTreeCont[kBcOut]->SetBranchAddress("v0",  src.v0);


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
