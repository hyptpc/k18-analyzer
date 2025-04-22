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
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"
#include "FourVectorFitter.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#define SaveRawData 1
#define DebugDisp 0
#define DoKF 1
#define SaveHistograms 1

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gTPC  = TPCParamMan::GetInstance();
const Int_t MaxTPCHits = 10000;

//For GenFit Setting
const bool Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//0~3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 150.;
const Double_t vtx_scan_rangeInsidePi = 50.;
const Double_t target_tarnsverse_dist_cut = 45.; //ref

const Double_t phi_masscut = 0.1; //ref
const Double_t k_vtx_distcut = 300;
const Double_t kk_distcut = 10.; //ref

const Double_t GFkk_distcut = 10.;



const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

int TPCToG4TrackID(std::vector<TVector3>TPCHit, int nhG4,int* tidG4, double* xG4,double* yG4,double* zG4 ,int& nhits){
  std::vector<TVector3> G4Hits;
  for(int ih=0;ih<nhG4;++ih){
    TVector3 G4Hit(xG4[ih],yG4[ih],zG4[ih]);
    G4Hits.push_back(G4Hit);
  }
  int MaxTracks = 1000;
  TH1I Counter("counter","counter",MaxTracks,0,MaxTracks);
  for(auto hit:TPCHit){
    double dl = 5000;
    int G4ID = -1;
    for(int ih=0;ih<nhG4;++ih){
      auto G4Hit = G4Hits.at(ih);
      double dist = (G4Hit - hit).Mag();
      if(dist < dl){
        dl = dist;
        G4ID = tidG4[ih];
      }
    }
    Counter.Fill(G4ID);
  }
  nhits = Counter.GetMaximum();
  int G4id = Counter.GetMaximumBin()-1;
  return G4id;
}
TVector3 GetG4Mom(TVector3 TPCHit, vector<TVector3> G4Hits,vector<TVector3>G4Moms){
  int nh = G4Hits.size();
  double dl = 5000;
  TVector3 mom;
  for(int ih=0;ih<nh;++ih){
    auto G4Hit = G4Hits.at(ih);
    double dist = (G4Hit - TPCHit).Mag();
    if(dist < dl){
      dl = dist;
      mom = G4Moms.at(ih);
    }
  }
  return mom;
}
int CountHits(int id, double* posx,double* posz,int* tidG4, int nh){
  int count = 0;
  for(int ih=0;ih<nh;++ih){
    double x = posx[ih],z=posz[ih];
    int pad = tpc::findPadID(z,x);
    int layer = tpc::getLayerID(pad);
    int row = tpc::getRowID(pad);
    double val = 0;
    gTPC.GetCDe(layer,row,1,val);
    if(val==0) continue;
    if(tidG4[ih] == id) count++;
  }
  return count;
}
TLorentzVector ToHelix(TLorentzVector GlobalLV){
  double E = GlobalLV.E();
  double X = -GlobalLV.X();
  double Y = GlobalLV.Z();
  double Z = GlobalLV.Y();
  return TLorentzVector(X,Y,Z,E);
}
TLorentzVector ToGlobal(TLorentzVector HelixLV){
  return ToHelix(HelixLV);
}
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kE42, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[E42]", "[OutFile]" };
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

  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
  std::vector<Double_t> tHtof;
  std::vector<Double_t> dtHtof;
  std::vector<Double_t> deHtof;
  std::vector<Double_t> posHtof;

  Int_t nhFtof;
  std::vector<Double_t> FtofSeg;
  std::vector<Double_t> tFtof;
  std::vector<Double_t> deFtof;
  std::vector<Double_t> posFtof;
  
  Int_t ntK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> xtgtK18;
  std::vector<Double_t> ytgtK18;
  std::vector<Double_t> utgtK18;
  std::vector<Double_t> vtgtK18;

  Int_t NumberOfTracks;
  Int_t PIDOfTrack[1000];
  Int_t ParentIDOfTrack[1000];
  Double_t VertexOfTrack_x[1000];
  Double_t VertexOfTrack_y[1000];
  Double_t VertexOfTrack_z[1000];
  Double_t MomentumOfTrack[1000];
  Double_t MomentumOfTrack_x[1000];
  Double_t MomentumOfTrack_y[1000];
  Double_t MomentumOfTrack_z[1000];
  
  std::vector<double> g4phimass;
  std::vector<double> g4phimom;
  std::vector<double> g4phimom_x;
  std::vector<double> g4phimom_y;
  std::vector<double> g4phimom_z;
  std::vector<double> g4phivtx_x;
  std::vector<double> g4phivtx_y;
  std::vector<double> g4phivtx_z;
  std::vector<double> g4phidecay_mom;
  std::vector<double> g4phidecay_g4tid;
  std::vector<double> g4phidecay_mom_x;
  std::vector<double> g4phidecay_mom_y;
  std::vector<double> g4phidecay_mom_z;

  int nPhi;
  std::vector<Bool_t> phiflag; 
  std::vector<Bool_t> good_phi; 
  std::vector<Bool_t> swap_phi; 
  std::vector<double> phimass;
  std::vector<double> phimom;
  std::vector<double> phimom_x;
  std::vector<double> phimom_y;
  std::vector<double> phimom_z;
  std::vector<double> phivtx_x;
  std::vector<double> phivtx_y;
  std::vector<double> phivtx_z;
  std::vector<double> kpkmdist;
  std::vector<double> phiphidist;
  
  std::vector<int> phidecay_id;
  std::vector<int> phidecay_purity;
  std::vector<int> phidecay_efficiency;
  std::vector<int> phidecay_g4tid;//kp1,km1,kp2,km2
  std::vector<int> phidecay_charge;//For clarification
  std::vector<double> phidecay_mom;
  std::vector<double> phidecay_mom_x;
  std::vector<double> phidecay_mom_y;
  std::vector<double> phidecay_mom_z;
  
  double MissMass;
  double MissMom;
  double MissMom_x;
  double MissMom_y;
  double MissMom_z;
  double MissMom_t;

  double delM;
  double delE;
  double delM4Pi;

  
  
  Int_t nclTpc;
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_layer;
  std::vector<Double_t> cluster_mrow;
  std::vector<Int_t> cluster_houghflag;
  std::vector<Int_t> cluster_G4tid;
  std::vector<Int_t> cluster_G4pid;
  
  Int_t remain_nclTpc; // Number of remain clusters not occupied in the tracks
  std::vector<Double_t> remain_cluster_x;
  std::vector<Double_t> remain_cluster_y;
  std::vector<Double_t> remain_cluster_z;
  std::vector<Double_t> remain_cluster_de;
  std::vector<Int_t> remain_cluster_layer;
  std::vector<Int_t> remain_cluster_houghflag;
  std::vector<Int_t> remain_cluster_G4tid;
  std::vector<Int_t> remain_cluster_G4pid;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> trackid; //for Kurama K1.8 tracks
  std::vector<Int_t> isXi;
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> isMultiloop;
  std::vector<Int_t> charge; //Helix charge
  std::vector<Int_t> pid;
  std::vector<Double_t> chisqr;
  std::vector<Double_t> pval;
  std::vector<Double_t> purity;
  std::vector<Double_t> efficiency;
  std::vector<Int_t> G4tid;//Geant4 Track ID of helix track
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx
  std::vector<Double_t> mom0; //Helix momentum at Y = 0
  std::vector<Double_t> path; //Helix path

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
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_mrow;




  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();


    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    dtHtof.clear();
    deHtof.clear();
    posHtof.clear();

    ntK18 = 0;
    pK18.clear();
    chisqrK18.clear();
    xtgtK18.clear();
    ytgtK18.clear();
    utgtK18.clear();
    vtgtK18.clear();
   
    g4phimass.clear();
    g4phimom.clear();
    g4phimom_x.clear();
    g4phimom_y.clear();
    g4phimom_z.clear();
    g4phivtx_x.clear();
    g4phivtx_y.clear();
    g4phivtx_z.clear();
    g4phidecay_mom.clear();
    g4phidecay_g4tid.clear();
    g4phidecay_mom_x.clear();
    g4phidecay_mom_y.clear();
    g4phidecay_mom_z.clear();



    nPhi = 0;
    phiflag.clear();
    swap_phi.clear();
    good_phi.clear();
    phimass.clear();
    phivtx_x.clear();
    phivtx_y.clear();
    phivtx_z.clear();
    phidecay_id.clear();
    phidecay_purity.clear();
    phidecay_efficiency.clear();
    phidecay_g4tid.clear();
    phidecay_charge.clear();
    phidecay_mom.clear();
    phidecay_mom_x.clear();
    phidecay_mom_y.clear();
    phidecay_mom_z.clear();


   
    
    
    nclTpc = 0;
    cluster_x.clear();
    cluster_y.clear();
    cluster_z.clear();
    cluster_de.clear();
    cluster_layer.clear();
    cluster_mrow.clear();
    cluster_houghflag.clear();
    cluster_G4tid.clear();
    cluster_G4pid.clear();
    
    remain_nclTpc=0;
    remain_cluster_x.clear();
    remain_cluster_y.clear();
    remain_cluster_z.clear();
    remain_cluster_de.clear();
    remain_cluster_layer.clear();
    remain_cluster_houghflag.clear(); 
    remain_cluster_G4tid.clear();
    remain_cluster_G4pid.clear();

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isXi.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    charge.clear();
    pid.clear();

    chisqr.clear();
    pval.clear();
    purity.clear();
    efficiency.clear();
    G4tid.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();
    mom0.clear();
    path.clear();



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
    track_cluster_size.clear();
    track_cluster_de.clear();
    track_cluster_mrow.clear();

    MissMass = qnan;
    MissMom = qnan;
    MissMom_x = qnan;
    MissMom_y = qnan;
    MissMom_z = qnan;
    MissMom_t = qnan;
    delM = qnan;
    delE = qnan;
    delM4Pi = qnan;
    for(int i=0;i<1000;++i){
      PIDOfTrack[i]=qnan;
      ParentIDOfTrack[i]=qnan;
      VertexOfTrack_x[i]=qnan;
      VertexOfTrack_y[i]=qnan;
      VertexOfTrack_z[i]=qnan;
      MomentumOfTrack[i]=qnan;
      MomentumOfTrack_x[i]=qnan;
      MomentumOfTrack_y[i]=qnan;
      MomentumOfTrack_z[i]=qnan;
    }

  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* evnum;
  
  TTreeReaderValue<Int_t>* nhHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;
  
	int nhFtof;
  int didFtof[500];
  double tFtof[500];
  double deFtof[500];
  double yFtof[500];//->posHtof

  int nhittpc;
  Int_t ititpc[MaxTPCHits];
  Double_t xtpc[MaxTPCHits];//with resolution
  Double_t ytpc[MaxTPCHits];//with resolution
  Double_t ztpc[MaxTPCHits];//with resolution
  Double_t pxtpc[MaxTPCHits];//with resolution
  Double_t pytpc[MaxTPCHits];//with resolution
  Double_t pztpc[MaxTPCHits];//with resolution

	int inc;

  Int_t NumberOfTracks;
  Int_t PIDOfTrack[1000];
  Int_t ParentIDOfTrack[1000];
  Double_t VertexOfTrack_x[1000];
  Double_t VertexOfTrack_y[1000];
  Double_t VertexOfTrack_z[1000];
  Double_t MomentumOfTrack[1000];
  Double_t MomentumOfTrack_x[1000];//In MeV scale
  Double_t MomentumOfTrack_y[1000];
  Double_t MomentumOfTrack_z[1000];




  TTreeReaderValue<Int_t>* nclTpc;
  TTreeReaderValue<std::vector<Double_t>>* cluster_x;  
  TTreeReaderValue<std::vector<Double_t>>* cluster_y;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de;
  TTreeReaderValue<std::vector<Int_t>>* cluster_layer;
  TTreeReaderValue<std::vector<Double_t>>* cluster_mrow;
  TTreeReaderValue<std::vector<Int_t>>* cluster_houghflag;
  TTreeReaderValue<std::vector<Int_t>>* cluster_G4tid;

  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isXi;
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;

  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<Double_t>>* purity;
  TTreeReaderValue<std::vector<Double_t>>* efficiency;
  TTreeReaderValue<std::vector<Int_t>>* G4tid;
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
  TTreeReaderValue<std::vector<Double_t>>* pval;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx; //reference dedx
  TTreeReaderValue<std::vector<Double_t>>* mom0;//Helix momentum at Y = 0
  TTreeReaderValue<std::vector<Double_t>>* path;//Helix path
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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* helix_t;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* alpha;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
	TTreeReaderValue<vector<double>>* xtgtHS=nullptr;
	TTreeReaderValue<vector<double>>* ytgtHS=nullptr;
	TTreeReaderValue<vector<double>>* xtgtKurama=nullptr;
	TTreeReaderValue<vector<double>>* ytgtKurama=nullptr;

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    TPCHid    = 100000,
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
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(), Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  GFTrackCont.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFTrackCont.GetVerbosity()<<std::endl;

#if 0
  GFTrackCont.DebugMode();
#endif

  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
    GFTrackCont.Clear();
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
  Int_t open_file = 0;
  Int_t open_tree = 0;
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
dst::DstRead( Int_t ievent )
{
  static const auto PhiMass = pdg::Mass(333);
  static const auto PionMass    = pdg::PionMass();
  static const auto KaonMass    = pdg::KaonMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto me = 0.001*0.5109989461;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};
  TVector3 tgtpos(0, 0, tpc::ZTarget);
  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);
  event.evnum = **src.evnum;
  HF1( 1, event.status++ );
  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;
  event.nhFtof = src.nhFtof;
  for(int ih=0;ih<event.nhFtof;++ih){
    event.FtofSeg.push_back(src.didFtof[ih]);
    event.tFtof.push_back(src.tFtof[ih]);
    event.deFtof.push_back(src.deFtof[ih]);
    event.posFtof.push_back(src.yFtof[ih]);
  }
  event.nclTpc = **src.nclTpc;
  event.cluster_x = **src.cluster_x;
  event.cluster_y = **src.cluster_y;
  event.cluster_z = **src.cluster_z;
  event.cluster_de = **src.cluster_de;
  event.cluster_layer = **src.cluster_layer;
  event.cluster_mrow = **src.cluster_mrow;
  event.cluster_houghflag = **src.cluster_houghflag;
  event.cluster_G4tid = **src.cluster_G4tid;
  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;
  Int_t icl_remain = 0;
  for( Int_t icl=0; icl<event.nclTpc; ++icl ){
    if(event.cluster_houghflag[icl]!=0) continue;
    event.remain_cluster_de.push_back(event.cluster_de[icl]);
    event.remain_cluster_x.push_back(event.cluster_x[icl]);
    event.remain_cluster_y.push_back(event.cluster_y[icl]);
    event.remain_cluster_z.push_back(event.cluster_z[icl]);
    event.remain_cluster_layer.push_back(event.cluster_layer[icl]);
    event.remain_cluster_houghflag.push_back(event.cluster_houghflag[icl]);
    event.remain_cluster_G4tid.push_back(event.cluster_G4tid[icl]); 
    /*
    event.remain_cluster_x[icl_remain] = event.cluster_x[icl];
    event.remain_cluster_y[icl_remain] = event.cluster_y[icl];
    event.remain_cluster_z[icl_remain] = event.cluster_z[icl];
    event.remain_cluster_de[icl_remain] = event.cluster_de[icl];
    event.remain_cluster_layer[icl_remain] = event.cluster_layer[icl];
    event.remain_cluster_houghflag[icl_remain] = event.cluster_houghflag[icl];
    */
    icl_remain++;
  }
  event.remain_nclTpc = icl_remain;
  HF1( 1, event.status++ );
  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isXi = **src.isXi;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.charge = **src.charge;
  event.pid = **src.pid;
  event.purity = **src.purity;
  event.G4tid = **src.G4tid;
  event.efficiency = **src.efficiency;
  event.chisqr = **src.chisqr;
  event.pval = **src.pval;
  event.helix_cx = **src.helix_cx;
  event.helix_cy = **src.helix_cy;
  event.helix_z0 = **src.helix_z0;
  event.helix_r = **src.helix_r;
  event.helix_dz = **src.helix_dz;
  event.dE = **src.dE;
  event.dEdx = **src.dEdx;
  event.mom0 = **src.mom0;
  event.path = **src.path;
  event.hitlayer = **src.hitlayer;
  event.hitpos_x = **src.hitpos_x;
  event.hitpos_y = **src.hitpos_y;
  event.hitpos_z = **src.hitpos_z;
  event.calpos_x = **src.calpos_x;
  event.calpos_y = **src.calpos_y;
  event.calpos_z = **src.calpos_z;
  event.residual = **src.residual;
  event.residual_x = **src.residual_x;
  event.residual_y = **src.residual_y;
  event.residual_z = **src.residual_z;
  event.resolution_x = **src.resolution_x;
  event.resolution_y = **src.resolution_y;
  event.resolution_z = **src.resolution_z;
  event.helix_t = **src.helix_t;
  event.alpha = **src.alpha;
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_mrow = **src.track_cluster_mrow;
  vector<TVector3> G4Hits;
  vector<TVector3> G4Moms;
  for(int ih=0;ih<src.nhittpc;++ih){
    TVector3 G4Hit(src.xtpc[ih],src.ytpc[ih],src.ztpc[ih]);
    TVector3 G4Mom(src.pxtpc[ih],src.pytpc[ih],src.pztpc[ih]);
    G4Hits.push_back(G4Hit);
    G4Moms.push_back(G4Mom);
  }
  for(int it=1;it<=src.NumberOfTracks;++it){
    int parent = src.ParentIDOfTrack[it];
    if(parent<0) continue;
    int pid = src.PIDOfTrack[it];
    int parent_pid = src.PIDOfTrack[parent];
    double mom = src.MomentumOfTrack[it]/1000;//MeV to GeV
    double mom_x = src.MomentumOfTrack_x[it]/1000;
    double mom_y = src.MomentumOfTrack_y[it]/1000;
    double mom_z = src.MomentumOfTrack_z[it]/1000;
    double vert_x = src.VertexOfTrack_x[it];
    double vert_y = src.VertexOfTrack_y[it];
    double vert_z = src.VertexOfTrack_z[it];

    if(it == 1){
      event.ntK18 = 1;
      event.pK18.push_back(mom);
      event.xtgtK18.push_back(vert_x);
      event.ytgtK18.push_back(vert_y);
      event.utgtK18.push_back(mom_x/mom);
      event.vtgtK18.push_back(mom_y/mom);
      event.chisqrK18.push_back(1);
    }
    if(pid == 333){
      if(it == 2){
        event.g4phimom.push_back(mom);
        event.g4phimom_x.push_back(mom_x);
        event.g4phimom_y.push_back(mom_y);
        event.g4phimom_z.push_back(mom_z);
        event.g4phivtx_x.push_back(vert_x);
        event.g4phivtx_y.push_back(vert_y);
        event.g4phivtx_z.push_back(vert_z);
      }
      if(it == 3){
        event.g4phimom.push_back(mom);
        event.g4phimom_x.push_back(mom_x);
        event.g4phimom_y.push_back(mom_y);
        event.g4phimom_z.push_back(mom_z);
        event.g4phivtx_x.push_back(vert_x);
        event.g4phivtx_y.push_back(vert_y);
        event.g4phivtx_z.push_back(vert_z);
      }
    }
    if(parent_pid == 333 and parent == 2 and pid == 321){
      event.g4phidecay_mom.push_back(mom);
      event.g4phidecay_mom_x.push_back(mom_x);
      event.g4phidecay_mom_y.push_back(mom_y);
      event.g4phidecay_mom_z.push_back(mom_z);
      event.g4phidecay_g4tid.push_back(it);
    }
    if(parent_pid == 333 and parent == 2 and pid == -321){
      event.g4phidecay_mom.push_back(mom);
      event.g4phidecay_mom_x.push_back(mom_x);
      event.g4phidecay_mom_y.push_back(mom_y);
      event.g4phidecay_mom_z.push_back(mom_z);
      event.g4phidecay_g4tid.push_back(it);
    }
    if(parent_pid == 333 and parent == 3 and pid == 321){
      event.g4phidecay_mom.push_back(mom);
      event.g4phidecay_mom_x.push_back(mom_x);
      event.g4phidecay_mom_y.push_back(mom_y);
      event.g4phidecay_mom_z.push_back(mom_z);
      event.g4phidecay_g4tid.push_back(it);
    }
    if(parent_pid == 333 and parent == 3 and pid == -321){
      event.g4phidecay_mom.push_back(mom);
      event.g4phidecay_mom_x.push_back(mom_x);
      event.g4phidecay_mom_y.push_back(mom_y);
      event.g4phidecay_mom_z.push_back(mom_z);
      event.g4phidecay_g4tid.push_back(it);
    }
  }
  if(event.g4phidecay_mom.size()==4){
    TVector3 g4kp1_mom(event.g4phidecay_mom_x[0],
        event.g4phidecay_mom_y[0],
        event.g4phidecay_mom_z[0]);
    TVector3 g4km1_mom(event.g4phidecay_mom_x[1],
        event.g4phidecay_mom_y[1],
        event.g4phidecay_mom_z[1]);
    TVector3 g4kp2_mom(event.g4phidecay_mom_x[2],
        event.g4phidecay_mom_y[2],
        event.g4phidecay_mom_z[2]);
    TVector3 g4km2_mom(event.g4phidecay_mom_x[3],
        event.g4phidecay_mom_y[3],
        event.g4phidecay_mom_z[3]);
    TLorentzVector g4kp1(g4kp1_mom, sqrt(g4kp1_mom.Mag2()+KaonMass*KaonMass));
    TLorentzVector g4km1(g4km1_mom, sqrt(g4km1_mom.Mag2()+KaonMass*KaonMass));
    TLorentzVector g4kp2(g4kp2_mom, sqrt(g4kp2_mom.Mag2()+KaonMass*KaonMass));
    TLorentzVector g4km2(g4km2_mom, sqrt(g4km2_mom.Mag2()+KaonMass*KaonMass));
    TLorentzVector g4phi1 = g4kp1 + g4km1;
    TLorentzVector g4phi2 = g4kp2 + g4km2;
    event.g4phimass.push_back(g4phi1.M());
    event.g4phimass.push_back(g4phi2.M());
  }
  vector<int> G4TrackID;
  vector<int> PureHits;
  HF1( 1, event.status++ );
  if( ntTpc == 0 ) return true;
  HF1( 1, event.status++ );
  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracksGeant4(**src.ntTpc,
      **src.charge, **src.nhtrack, **src.helix_cx,
      **src.helix_cy, **src.helix_z0, **src.helix_r,
      **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
      **src.helix_t, **src.track_cluster_de, **src.resolution_x,
      **src.resolution_y, **src.resolution_z, **src.hitpos_x,
      **src.hitpos_y, **src.hitpos_z);
  HF1( 1, event.status++ );
  HypTPCTask& GFTrackCont = HypTPCTask::GetInstance();
  for( Int_t it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
    GFTrackCont.AddHelixTrack(pdgcode, tp);
    vector<TVector3> TPCHit;
    for(int ih=0;ih<event.hitpos_x.at(it).size();++ih){
      TPCHit.push_back(TVector3(
            event.hitpos_x.at(it).at(ih),
            event.hitpos_y.at(it).at(ih),
            event.hitpos_z.at(it).at(ih)));
    }
    int nPureHit;
    int G4tid = TPCToG4TrackID(TPCHit,src.nhittpc,src.ititpc,src.xtpc,src.ytpc,src.ztpc,nPureHit);
    G4TrackID.push_back(G4tid);
    PureHits.push_back(nPureHit);
  }
  GFTrackCont.FitTracks();
  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }
  HF1( 1, event.status++ );



  //Recon Process//


  //phi_*_container are for single phi events.
  vector<int> phi_kp_id_container;
  vector<int> phi_kp_repid_container;
  vector<int> phi_kp_g4tid_container;
  vector<double> phi_kp_purity_container;
  vector<double> phi_kp_efficiency_container;
  vector<TVector3> phi_kp_mom_container;

  vector<int> phi_km_id_container;
  vector<int> phi_km_repid_container;
  vector<int> phi_km_g4tid_container;
  vector<double> phi_km_purity_container;
  vector<double> phi_km_efficiency_container;
  vector<TVector3> phi_km_mom_container;
  
  vector<double> phi_mass_container;
  vector<TVector3> phi_mom_container;
  vector<TVector3> phi_vertex_container;
  vector<double> phi_kkdist_container;

  //d(ouble)phi_*_container are for double phi events.
  vector<int> dphi_kp_id_container;
  vector<int> dphi_kp_repid_container;
  vector<int> dphi_kp_g4tid_container;
  vector<double> dphi_kp_purity_container;
  vector<double> dphi_kp_efficiency_container;
  vector<TVector3> dphi_kp_mom_container;

  vector<int> dphi_km_id_container;
  vector<int> dphi_km_repid_container;
  vector<int> dphi_km_g4tid_container;
  vector<double> dphi_km_purity_container;
  vector<double> dphi_km_efficiency_container;
  vector<TVector3> dphi_km_mom_container;
  
  vector<double> dphi_mass_container;
  vector<TVector3> dphi_mom_container;
  vector<TVector3> dphi_vertex_container;
  vector<double> dphi_kkdist_container;
  {

    for(Int_t it1=0;it1<ntTpc;it1++){
    #if DebugDisp
      cout<<Form("Processing track %d / %d",it1,ntTpc)<<endl;
    #endif
      if((event.pid[it1]&2)!=2) continue; //select k+like
      if(event.charge[it1]!=1) continue;

      Int_t repid_kp1 = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it1];
	if(temp==flag) repid_kp1 += 1;
	flag*=2;
      }
      cout<<Form("Checking Track %d, repid = %d",it1,repid_kp1)<<endl;
      if(!GFTrackCont.TrackCheck(it1, repid_kp1)) continue;
      cout<<Form("Track %d, repid = %d checked",it1,repid_kp1)<<endl;

      Double_t kp1_par[5];
      kp1_par[0] = event.helix_cx[it1];
      kp1_par[1] = event.helix_cy[it1];
      kp1_par[2] = event.helix_z0[it1];
      kp1_par[3] = event.helix_r[it1];
      kp1_par[4] = event.helix_dz[it1];
      Int_t kp1_nh = event.helix_t[it1].size();
      Double_t kp1_theta_min = event.helix_t[it1][0] - vtx_scan_range/kp1_par[3];
      Double_t kp1_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_range/kp1_par[3], event.helix_t[it1][kp1_nh-1]);
      TVector3 kp1_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 kp1_end = TVector3(event.calpos_x[it1][kp1_nh-1], event.calpos_y[it1][kp1_nh-1], event.calpos_z[it1][kp1_nh-1]);

      for(Int_t it2=0;it2<ntTpc;it2++){
    #if DebugDisp
      cout<<Form("Processing track (%d,%d) / %d",it1,it2,ntTpc)<<endl;
    #endif
	if(it1==it2) continue;
	if((event.pid[it2]&2)!=2) continue; //select k- like
	if(event.charge[it2]!=-1) continue;

	Int_t repid_km1 = 0;
  flag = 1;
  for(Int_t i=0;i<2;i++){
    Int_t temp = flag&event.pid[it2];
    if(temp==flag) repid_km1 += 1;
    flag*=2;
  }
  cout<<Form("Checking Track %d, repid = %d",it2,repid_km1)<<endl;
	if(!GFTrackCont.TrackCheck(it2, repid_km1)) continue;
  cout<<Form("Track %d, repid = %d checked",it2,repid_km1)<<endl;

	Double_t km1_par[5];
	km1_par[0] = event.helix_cx[it2];
	km1_par[1] = event.helix_cy[it2];
	km1_par[2] = event.helix_z0[it2];
	km1_par[3] = event.helix_r[it2];
	km1_par[4] = event.helix_dz[it2];

	Int_t km1_nh = event.helix_t[it2].size();
	Double_t km1_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/km1_par[3], event.helix_t[it2][km1_nh-1]);
	Double_t km1_theta_max = event.helix_t[it2][0] + vtx_scan_range/km1_par[3];
	TVector3 km1_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 km1_end = TVector3(event.calpos_x[it2][km1_nh-1], event.calpos_y[it2][km1_nh-1], event.calpos_z[it2][km1_nh-1]);

	Double_t kp1km1_dist = 10000.;
	TVector3 kp1_mom; TVector3 km1_mom; TVector3 phi1_mom;
	TVector3 phi1_vert = Kinematics::LambdaVertex(dMagneticField, kp1_par, km1_par,
  kp1_theta_min, kp1_theta_max, km1_theta_min, km1_theta_max,
  kp1_mom, km1_mom, phi1_mom, kp1km1_dist);
	if(TMath::IsNaN(kp1km1_dist)) continue;
	phi1_mom = km1_mom + kp1_mom;

	TLorentzVector Lkp1(kp1_mom, TMath::Hypot(kp1_mom.Mag(), KaonMass));
	TLorentzVector Lkm1(km1_mom, TMath::Hypot(km1_mom.Mag(), KaonMass));
	TLorentzVector Lphi1 = Lkp1+Lkm1;

	if(TMath::Abs(phi1_vert.x()) > 250. || TMath::Abs(phi1_vert.z()) > 250. || TMath::Abs(phi1_vert.y()) > 250.) continue; //Vertex cut
	if(kp1km1_dist > kk_distcut) continue;
	Double_t km1_vertex_dist; Double_t kp1_vertex_dist;
	if(!Kinematics::HelixDirection(phi1_vert, kp1_start, kp1_end, kp1_vertex_dist) ||
	   !Kinematics::HelixDirection(phi1_vert, km1_start, km1_end, km1_vertex_dist)) continue;
	
	if(km1_vertex_dist > k_vtx_distcut) continue;
 	if(kp1_vertex_dist > k_vtx_distcut) continue;
	if(TMath::Abs(Lphi1.M() - PhiMass) > phi_masscut) continue;

	Double_t phi1target_dist;
	TVector3 phi1target_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
							       phi1_vert,
							       phi1_mom,
							       phi1target_dist);
//  if(hypot(phi1target_vtx.x(), phi1target_vtx.z()) > target_tarnsverse_dist_cut) continue;
  phi_kp_id_container.push_back(it1);
  phi_kp_repid_container.push_back(repid_kp1);
  phi_kp_purity_container.push_back(event.purity[it1]);
  phi_kp_efficiency_container.push_back(event.efficiency[it1]);
  phi_kp_g4tid_container.push_back(G4TrackID[it1]);
  phi_kp_mom_container.push_back(kp1_mom);
  
  phi_km_id_container.push_back(it2);
  phi_km_repid_container.push_back(repid_km1);
  phi_km_purity_container.push_back(event.purity[it2]);
  phi_km_efficiency_container.push_back(event.efficiency[it2]);
  phi_km_g4tid_container.push_back(G4TrackID[it2]);
  phi_km_mom_container.push_back(km1_mom);
  
  phi_mass_container.push_back(Lphi1.M());
  phi_mom_container.push_back(phi1_mom);
  phi_kkdist_container.push_back(kp1km1_dist);
  phi_vertex_container.push_back(phi1_vert);                 

	for(int it3=0;it3<ntTpc;it3++){
    #if DebugDisp
      cout<<Form("Processing track (%d,%d,%d) / %d",it1,it2,it3,ntTpc)<<endl;
    #endif
	  if(it3==it2 or it3==it1) continue;
	  if((event.pid[it3]&2)!=2) continue; //select p like
	  if(event.charge[it3]!=1) continue;
      
    Int_t repid_kp2 = 0;
    flag = 1;
    for(Int_t i=0;i<2;i++){
      Int_t temp = flag&event.pid[it3];
      if(temp==flag) repid_kp2 += 1;
      flag*=2;
    }
    cout<<Form("Checking Track %d, repid = %d",it3,repid_kp2)<<endl;
    if(!GFTrackCont.TrackCheck(it3, repid_kp2)) continue;
    cout<<Form("Track %d, repid = %d checked",it3,repid_kp2)<<endl;

    Double_t kp2_par[5];
    kp2_par[0] = event.helix_cx[it3];
    kp2_par[1] = event.helix_cy[it3];
    kp2_par[2] = event.helix_z0[it3];
    kp2_par[3] = event.helix_r[it3];
    kp2_par[4] = event.helix_dz[it3];
    Int_t kp2_nh = event.helix_t[it3].size();
    Double_t kp2_theta_min = event.helix_t[it3][0] - vtx_scan_range/kp2_par[3];
    Double_t kp2_theta_max = TMath::Min(event.helix_t[it3][0] + vtx_scan_range/kp2_par[3], event.helix_t[it3][kp2_nh-1]);
    TVector3 kp2_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
    TVector3 kp2_end = TVector3(event.calpos_x[it3][kp2_nh-1], event.calpos_y[it3][kp2_nh-1], event.calpos_z[it3][kp2_nh-1]);
    for(int it4=0;it4<ntTpc;++it4){
      #if DebugDisp
      cout<<Form("Processing track (%d,%d,%d,%d) / %d",it1,it2,it3,it4,ntTpc)<<endl;
      #endif
      if(it4==it3 or it4==it2 or it4==it1) continue;
      if((event.pid[it4]&2)!=2) continue; //select K-like
      if(event.charge[it4]!=-1) continue;

      Int_t repid_km2 = 0;
      flag = 1;
      for(Int_t i=0;i<2;i++){
        Int_t temp = flag&event.pid[it4];
        if(temp==flag) repid_km2 += 1;
        flag*=2;
      }
      cout<<Form("Checking Track %d, repid = %d",it4,repid_km2)<<endl;
      if(!GFTrackCont.TrackCheck(it4, repid_km2)) continue;
      cout<<Form("Track %d, repid = %d checked",it4,repid_km2)<<endl;

      Double_t km2_par[5];
      km2_par[0] = event.helix_cx[it4];
      km2_par[1] = event.helix_cy[it4];
      km2_par[2] = event.helix_z0[it4];
      km2_par[3] = event.helix_r[it4];
      km2_par[4] = event.helix_dz[it4];

      Int_t km2_nh = event.helix_t[it4].size();
      Double_t km2_theta_min = TMath::Max(event.helix_t[it4][0] - vtx_scan_rangeInsideL/km2_par[3], event.helix_t[it4][km2_nh-1]);
      Double_t km2_theta_max = event.helix_t[it4][0] + vtx_scan_range/km2_par[3];
      TVector3 km2_start = TVector3(event.calpos_x[it4][0], event.calpos_y[it4][0], event.calpos_z[it4][0]);
      TVector3 km2_end = TVector3(event.calpos_x[it4][km2_nh-1], event.calpos_y[it4][km2_nh-1], event.calpos_z[it4][km2_nh-1]);
    
      Double_t kp2km2_dist = 10000.;
      TVector3 kp2_mom; TVector3 km2_mom; TVector3 phi2_mom;
      TVector3 phi2_vert = Kinematics::LambdaVertex(dMagneticField, kp2_par, km2_par,
      kp2_theta_min, kp2_theta_max, km2_theta_min, km2_theta_max,
      kp2_mom, km2_mom, phi2_mom, kp2km2_dist);
      if(TMath::IsNaN(kp2km2_dist)) continue;
      phi2_mom = km2_mom + kp2_mom;

      TLorentzVector Lkp2(kp2_mom, TMath::Hypot(kp2_mom.Mag(), KaonMass));
      TLorentzVector Lkm2(km2_mom, TMath::Hypot(km2_mom.Mag(), KaonMass));
      TLorentzVector Lphi2 = Lkp2+Lkm2;

      if(TMath::Abs(phi2_vert.x()) > 250. || TMath::Abs(phi2_vert.z()) > 250. || TMath::Abs(phi2_vert.y()) > 250.) continue; //Vertex cut
      if(kp2km2_dist > kk_distcut) continue;
      Double_t km2_vertex_dist; Double_t kp2_vertex_dist;
      if(!Kinematics::HelixDirection(phi2_vert, kp2_start, kp2_end, kp2_vertex_dist) ||
        !Kinematics::HelixDirection(phi2_vert, km2_start, km2_end, km2_vertex_dist)) continue;
      
      if(km2_vertex_dist > k_vtx_distcut) continue;
      if(kp2_vertex_dist > k_vtx_distcut) continue;
      if(TMath::Abs(Lphi2.M() - PhiMass) > phi_masscut) continue;

      Double_t phi2target_dist;
      TVector3 phi2target_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
                        phi2_vert,
                        phi2_mom,
                        phi2target_dist);
//      if(hypot(phi2target_vtx.x(), phi2target_vtx.z()) > target_tarnsverse_dist_cut) continue;
      dphi_kp_id_container.push_back(it1);
      dphi_kp_repid_container.push_back(repid_kp1);
      dphi_kp_purity_container.push_back(event.purity[it1]);
      dphi_kp_efficiency_container.push_back(event.efficiency[it1]);
      dphi_kp_g4tid_container.push_back(G4TrackID[it1]);
      dphi_kp_mom_container.push_back(kp1_mom);
      dphi_kp_id_container.push_back(it3);
      dphi_kp_repid_container.push_back(repid_kp2);
      dphi_kp_purity_container.push_back(event.purity[it3]);
      dphi_kp_efficiency_container.push_back(event.efficiency[it3]);
      dphi_kp_g4tid_container.push_back(G4TrackID[it3]);
      dphi_kp_mom_container.push_back(kp2_mom);
      
      dphi_km_id_container.push_back(it2);
      dphi_km_repid_container.push_back(repid_km1);
      dphi_km_purity_container.push_back(event.purity[it2]);
      dphi_km_efficiency_container.push_back(event.efficiency[it2]);
      dphi_km_g4tid_container.push_back(G4TrackID[it2]);
      dphi_km_mom_container.push_back(km1_mom);
      dphi_km_id_container.push_back(it4);
      dphi_km_repid_container.push_back(repid_km2);
      dphi_km_purity_container.push_back(event.purity[it4]);
      dphi_km_efficiency_container.push_back(event.efficiency[it4]);
      dphi_km_g4tid_container.push_back(G4TrackID[it4]);
      dphi_km_mom_container.push_back(km2_mom);
      
      dphi_mass_container.push_back(Lphi1.M());
      dphi_mom_container.push_back(phi1_mom);
      dphi_kp_mom_container.push_back(kp1_mom);
      dphi_km_mom_container.push_back(km1_mom);
      dphi_kkdist_container.push_back(kp1km1_dist);
      dphi_vertex_container.push_back(phi1_vert);                 
      
      dphi_mass_container.push_back(Lphi2.M());
      dphi_mom_container.push_back(phi2_mom);
      dphi_kp_mom_container.push_back(kp2_mom);
      dphi_km_mom_container.push_back(km2_mom);
      dphi_kkdist_container.push_back(kp2km2_dist);
      dphi_vertex_container.push_back(phi2_vert);                 
    
    }//it4
	} //it3
      } //it2
    } //it1
  }
#if DebugDisp
  std::cout<<Form("Ks, KsP candidates searching finished: nks, nksP = (%d,%d) ",ks_candidates,ksp_candidates)<<std::endl;
#endif
  //Candidate Sorting w Helix Fit


  //GF for single Phi
  int phi_candidates = phi_kp_id_container.size();
  std::vector<Int_t> GFphi_kp_id_container(phi_candidates, -1);
  std::vector<Int_t> GFphi_kp_g4tid_container(phi_candidates, -1);
  std::vector<Int_t> GFphi_kp_repid_container(phi_candidates, qnan);
  std::vector<Double_t> GFphi_kp_purity_container(phi_candidates, qnan);
  std::vector<Double_t> GFphi_kp_efficiency_container(phi_candidates, qnan);
  std::vector<TVector3> GFphi_kp_mom_container(phi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFphi_kp_extrapolation_container(phi_candidates, qnan);

  std::vector<Int_t> GFphi_km_id_container(phi_candidates, -1);
  std::vector<Int_t> GFphi_km_g4tid_container(phi_candidates, -1);
  std::vector<Int_t> GFphi_km_repid_container(phi_candidates, qnan);
  std::vector<Double_t> GFphi_km_purity_container(phi_candidates, qnan);
  std::vector<Double_t> GFphi_km_efficiency_container(phi_candidates, qnan);
  std::vector<TVector3> GFphi_km_mom_container(phi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFphi_km_extrapolation_container(phi_candidates, qnan);
  
  std::vector<Double_t> GFphi_mass_container(phi_candidates, qnan);
  std::vector<Double_t> GFphi_kkdist_container(phi_candidates, qnan);
  std::vector<TVector3> GFphi_mom_container(phi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFphi_vtx_container(phi_candidates, TVector3(qnan, qnan, qnan));

  if(phi_candidates>0){
#if DebugDisp
  cout<<"Sorting KS Candidates"<<endl;
#endif
    //Reconstructed Single phi
    for(int idp=0;idp<phi_candidates;idp++){
      #if DebugDisp
      cout<<Form("Processing Ks candidate %d / %d",idp,ks_candidates)<<endl;
      #endif
      Int_t kp1_id = phi_kp_id_container[idp];
      Int_t kp1_repid = phi_kp_repid_container[idp];
      Int_t kp1_g4tid = phi_kp_g4tid_container[idp];
      Double_t kp1_purity = phi_kp_purity_container[idp];
      Double_t kp1_efficiency = phi_kp_efficiency_container[idp];
      
      Int_t km1_id = phi_km_id_container[idp];
      Int_t km1_repid = phi_km_repid_container[idp];
      Int_t km1_g4tid = phi_km_g4tid_container[idp];
      Double_t km1_purity = phi_km_purity_container[idp];
      Double_t km1_efficiency = phi_km_efficiency_container[idp];
     
      
      Double_t kp1_extrapolation; Double_t km1_extrapolation;
      TVector3 kp1_mom; TVector3 km1_mom;
      Double_t kp1km1_dist; TVector3 phi_vertex;
      Bool_t vtxcut = (GFTrackCont.FindVertex(kp1_id, km1_id,
					      kp1_repid, km1_repid,
					      kp1_extrapolation, km1_extrapolation,
					      kp1_mom, km1_mom,
					      kp1km1_dist, phi_vertex,
					      vtx_scan_range)
		       and kp1km1_dist < GFkk_distcut and !TMath::IsNaN(kp1km1_dist) );
      if(!vtxcut) continue;

      TVector3 phi_mom = kp1_mom + km1_mom;
      TLorentzVector GFLkp1(kp1_mom, TMath::Hypot(kp1_mom.Mag(), KaonMass));
      TLorentzVector GFLkm1(km1_mom, TMath::Hypot(km1_mom.Mag(), KaonMass));
      TLorentzVector GFLks = GFLkp1 + GFLkm1;

      GFphi_kp_id_container[idp] = kp1_id;
      GFphi_kp_repid_container[idp] = kp1_repid;
      GFphi_kp_g4tid_container[idp] = kp1_g4tid;
      GFphi_kp_purity_container[idp] = kp1_purity;
      GFphi_kp_efficiency_container[idp] = kp1_efficiency;
      GFphi_kp_mom_container[idp] = kp1_mom;
      GFphi_kp_extrapolation_container[idp] = kp1_extrapolation;
      
      GFphi_km_id_container[idp] = km1_id;
      GFphi_km_repid_container[idp] = km1_repid;
      GFphi_km_g4tid_container[idp] = km1_g4tid;
      GFphi_km_purity_container[idp] = km1_purity;
      GFphi_km_efficiency_container[idp] = km1_efficiency;
      GFphi_km_mom_container[idp] = km1_mom;
      GFphi_km_extrapolation_container[idp] = km1_extrapolation;
      
      GFphi_mass_container[idp] = GFLks.M();
      GFphi_mom_container[idp] = phi_mom;
      GFphi_kkdist_container[idp] = kp1km1_dist;
      GFphi_vtx_container[idp] = phi_vertex;
    }
  } //phi_candidates

  int dphi_candidates = dphi_kp_id_container.size()/2;
  std::vector<Int_t> GFdphi_kp_id_container(2*dphi_candidates, -1);
  std::vector<Int_t> GFdphi_kp_g4tid_container(2*dphi_candidates, -1);
  std::vector<Int_t> GFdphi_kp_repid_container(2*dphi_candidates, -1);
  std::vector<Double_t> GFdphi_kp_purity_container(2*dphi_candidates, qnan);
  std::vector<Double_t> GFdphi_kp_efficiency_container(2*dphi_candidates, qnan);
  std::vector<TVector3> GFdphi_kp_mom_container(2*dphi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFdphi_kp_extrapolation_container(2*dphi_candidates, qnan);

  std::vector<Int_t> GFdphi_km_id_container(2*dphi_candidates, -1);
  std::vector<Int_t> GFdphi_km_g4tid_container(2*dphi_candidates, -1);
  std::vector<Int_t> GFdphi_km_repid_container(2*dphi_candidates, -1);
  std::vector<Double_t> GFdphi_km_purity_container(2*dphi_candidates, qnan);
  std::vector<Double_t> GFdphi_km_efficiency_container(2*dphi_candidates, qnan);
  std::vector<TVector3> GFdphi_km_mom_container(2*dphi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFdphi_km_extrapolation_container(2*dphi_candidates, qnan);
  
  std::vector<Double_t> GFdphi_mass_container(2*dphi_candidates, qnan);
  std::vector<TVector3> GFdphi_mom_container(2*dphi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFdphi_vtx_container(2*dphi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFdphi_kkdist_container(2*dphi_candidates, qnan);
  std::vector<Double_t> GFdphi_phiphidist_container(dphi_candidates, qnan);


  for(Int_t candi=0;candi<dphi_candidates;candi++){
    #if DebugDisp
    cout<<Form("Processing KsP candidate %d / %d",candi,ksp_candidates)<<endl;
    #endif
    Int_t kp1_id = dphi_kp_id_container[2*candi];
    Int_t kp1_repid = dphi_kp_repid_container[2*candi];
    Int_t kp1_g4tid = dphi_kp_g4tid_container[2*candi];
    Double_t kp1_purity = dphi_kp_purity_container[2*candi];
    Double_t kp1_efficiency = dphi_kp_efficiency_container[2*candi];

    Int_t km1_id = dphi_km_id_container[2*candi];
    Int_t km1_repid = dphi_km_repid_container[2*candi];
    Int_t km1_g4tid = dphi_km_g4tid_container[2*candi];
    Double_t km1_purity = dphi_km_purity_container[2*candi];
    Double_t km1_efficiency = dphi_km_efficiency_container[2*candi];

    Int_t kp2_id = dphi_kp_id_container[2*candi+1];
    Int_t kp2_repid = dphi_kp_repid_container[2*candi+1];
    Int_t kp2_g4tid = dphi_kp_g4tid_container[2*candi+1];
    Double_t kp2_purity = dphi_kp_purity_container[2*candi+1];
    Double_t kp2_efficiency = dphi_kp_efficiency_container[2*candi+1];

    Int_t km2_id = dphi_km_id_container[2*candi+1];
    Int_t km2_repid = dphi_km_repid_container[2*candi+1];
    Int_t km2_g4tid = dphi_km_g4tid_container[2*candi+1];
    Double_t km2_purity = dphi_km_purity_container[2*candi+1];
    Double_t km2_efficiency = dphi_km_efficiency_container[2*candi+1];
      
    Double_t kp1_extrapolation; Double_t km1_extrapolation;
    TVector3 kp1_mom; TVector3 km1_mom;
    Double_t kp1km1_dist; TVector3 phi1_vertex;
    Bool_t vtxcut1 = (GFTrackCont.FindVertex(kp1_id, km1_id,
              kp1_repid, km1_repid,
              kp1_extrapolation, km1_extrapolation,
              kp1_mom, km1_mom,
              kp1km1_dist, phi1_vertex,
              vtx_scan_range)
          and kp1km1_dist < GFkk_distcut and !TMath::IsNaN(kp1km1_dist) );
    if(!vtxcut1) continue;
    
    Double_t kp2_extrapolation; Double_t km2_extrapolation;
    TVector3 kp2_mom; TVector3 km2_mom;
    Double_t kp2km2_dist; TVector3 phi2_vertex;
    Bool_t vtxcut2 = (GFTrackCont.FindVertex(kp2_id, km2_id,
              kp2_repid, km2_repid,
              kp2_extrapolation, km2_extrapolation,
              kp2_mom, km2_mom,
              kp2km2_dist, phi2_vertex,
              vtx_scan_range)
          and kp2km2_dist < GFkk_distcut and !TMath::IsNaN(kp2km2_dist) );
    if(!vtxcut2) continue;
    TVector3 phi1_mom = kp1_mom + km1_mom;
    TLorentzVector GFLkp1(kp1_mom, TMath::Hypot(kp1_mom.Mag(), KaonMass));
    TLorentzVector GFLkm1(km1_mom, TMath::Hypot(km1_mom.Mag(), KaonMass));
    TLorentzVector GFLphi1 = GFLkp1 + GFLkm1;
    
    TVector3 phi2_mom = kp2_mom + km2_mom;
    TLorentzVector GFLkp2(kp2_mom, TMath::Hypot(kp2_mom.Mag(), KaonMass));
    TLorentzVector GFLkm2(km2_mom, TMath::Hypot(km2_mom.Mag(), KaonMass));
    TLorentzVector GFLphi2 = GFLkp2 + GFLkm2;

    GFdphi_kp_id_container[2*candi] = kp1_id;
    GFdphi_kp_repid_container[2*candi] = kp1_repid;
    GFdphi_kp_g4tid_container[2*candi] = kp1_g4tid;
    GFdphi_kp_purity_container[2*candi] = kp1_purity;
    GFdphi_kp_efficiency_container[2*candi] = kp1_efficiency;
    GFdphi_kp_mom_container[2*candi] = kp1_mom;
    GFdphi_kp_extrapolation_container[2*candi] = kp1_extrapolation;

    GFdphi_km_id_container[2*candi] = km1_id;
    GFdphi_km_repid_container[2*candi] = km1_repid;
    GFdphi_km_g4tid_container[2*candi] = km1_g4tid;
    GFdphi_km_purity_container[2*candi] = km1_purity;
    GFdphi_km_efficiency_container[2*candi] = km1_efficiency;
    GFdphi_km_mom_container[2*candi] = km1_mom;
    GFdphi_km_extrapolation_container[2*candi] = km1_extrapolation;

    GFdphi_kp_id_container[2*candi+1] = kp2_id;
    GFdphi_kp_repid_container[2*candi+1] = kp2_repid;
    GFdphi_kp_g4tid_container[2*candi+1] = kp2_g4tid;
    GFdphi_kp_purity_container[2*candi+1] = kp2_purity;
    GFdphi_kp_efficiency_container[2*candi+1] = kp2_efficiency;
    GFdphi_kp_mom_container[2*candi+1] = kp2_mom;
    GFdphi_kp_extrapolation_container[2*candi+1] = kp2_extrapolation;

    GFdphi_km_id_container[2*candi+1] = km2_id;
    GFdphi_km_repid_container[2*candi+1] = km2_repid;
    GFdphi_km_g4tid_container[2*candi+1] = km2_g4tid;
    GFdphi_km_purity_container[2*candi+1] = km2_purity;
    GFdphi_km_efficiency_container[2*candi+1] = km2_efficiency;
    GFdphi_km_mom_container[2*candi+1] = km2_mom;
    GFdphi_km_extrapolation_container[2*candi+1] = km2_extrapolation;

    GFdphi_mass_container[2*candi] = GFLphi1.M();
    GFdphi_mass_container[2*candi+1] = GFLphi2.M();
    GFdphi_mom_container[2*candi] = phi1_mom;
    GFdphi_mom_container[2*candi+1] = phi2_mom;
    GFdphi_kkdist_container[2*candi] = kp1km1_dist;
    GFdphi_kkdist_container[2*candi+1] = kp2km2_dist;
    GFdphi_vtx_container[2*candi] = phi1_vertex;
    GFdphi_vtx_container[2*candi+1] = phi2_vertex;
    GFdphi_phiphidist_container[candi] = (phi1_vertex - phi2_vertex).Mag();

  } //dphicandidates
  double dm = 9999;
  int best_phi = -1;
  for(int ip=0;ip<phi_candidates;++ip){
    if(TMath::Abs(GFphi_mass_container[ip] - PhiMass) < dm){
      dm = TMath::Abs(GFphi_mass_container[ip] - PhiMass);
      best_phi = ip;
    }
  }
  double dm2 = 9999;
  int best_dphi = -1;
  for(int ip=0;ip<dphi_candidates;++ip){
    if(hypot(TMath::Abs(GFdphi_mass_container[2*ip] - PhiMass),TMath::Abs(GFdphi_mass_container[2*ip+1]- PhiMass) )< dm2){
      dm2 = hypot(TMath::Abs(GFdphi_mass_container[2*ip] - PhiMass),TMath::Abs(GFdphi_mass_container[2*ip+1]- PhiMass) );  
      best_dphi = ip;
    }
  }
  if(best_phi>=0 and best_dphi == -1){
    event.nPhi = 1;
    event.phiflag.push_back(1);
    TVector3 phi_mom = GFphi_mom_container[best_phi];
    TVector3 phi_vtx = GFphi_vtx_container[best_phi];
    event.phimass.push_back(GFphi_mass_container[best_phi]);
    event.phimom.push_back(phi_mom.Mag());
    event.phimom_x.push_back(phi_mom.x());
    event.phimom_y.push_back(phi_mom.y());
    event.phimom_z.push_back(phi_mom.z());
    event.phivtx_x.push_back(phi_vtx.x());
    event.phivtx_y.push_back(phi_vtx.y());
    event.phivtx_z.push_back(phi_vtx.z());
    event.kpkmdist.push_back(GFphi_kkdist_container[best_phi]);

    event.phidecay_id.push_back(GFphi_kp_id_container[best_phi]);
    event.phidecay_purity.push_back(GFphi_kp_purity_container[best_phi]);
    event.phidecay_efficiency.push_back(GFphi_kp_efficiency_container[best_phi]);
    event.phidecay_g4tid.push_back(GFphi_kp_g4tid_container[best_phi]);
    event.phidecay_charge.push_back(1);
    TVector3 kp1_mom = GFphi_kp_mom_container[best_phi];
    event.phidecay_mom.push_back(kp1_mom.Mag());
    event.phidecay_mom_x.push_back(kp1_mom.x());
    event.phidecay_mom_y.push_back(kp1_mom.y());
    event.phidecay_mom_z.push_back(kp1_mom.z());

    event.phidecay_id.push_back(GFphi_km_id_container[best_phi]);
    event.phidecay_purity.push_back(GFphi_km_purity_container[best_phi]);
    event.phidecay_efficiency.push_back(GFphi_km_efficiency_container[best_phi]);
    event.phidecay_g4tid.push_back(GFphi_km_g4tid_container[best_phi]);
    event.phidecay_charge.push_back(-1);
    TVector3 km1_mom = GFphi_km_mom_container[best_phi];
    event.phidecay_mom.push_back(km1_mom.Mag());
    event.phidecay_mom_x.push_back(km1_mom.x());
    event.phidecay_mom_y.push_back(km1_mom.y());
    event.phidecay_mom_z.push_back(km1_mom.z());
   

  }
  else if(best_dphi>=0){
    event.nPhi = 2;
    event.phiflag.push_back(1);
    event.phiflag.push_back(1);
    int g4kp1_g4tid = event.g4phidecay_g4tid[0];
    int g4kp2_g4tid = event.g4phidecay_g4tid[2];
    int kp1_g4tid = GFdphi_kp_g4tid_container[2*best_dphi];
    int kp2_g4tid = GFdphi_kp_g4tid_container[2*best_dphi+1];
    int phi1_id,phi2_id;
    if(g4kp1_g4tid == kp2_g4tid or g4kp2_g4tid == kp1_g4tid){
      phi2_id = 2*best_dphi+1;
      phi1_id = 2*best_dphi;
    }
    else{
      phi1_id = 2*best_dphi;
      phi2_id = 2*best_dphi+1;
    }

    TVector3 phi1_mom = GFdphi_mom_container[phi1_id];
    TVector3 phi1_vtx = GFdphi_vtx_container[phi1_id];
    event.phimass.push_back(GFdphi_mass_container[phi1_id]);
    event.phimom.push_back(phi1_mom.Mag());
    event.phimom_x.push_back(phi1_mom.x());
    event.phimom_y.push_back(phi1_mom.y());
    event.phimom_z.push_back(phi1_mom.z());
    event.phivtx_x.push_back(phi1_vtx.x());
    event.phivtx_y.push_back(phi1_vtx.y());
    event.phivtx_z.push_back(phi1_vtx.z());
    event.kpkmdist.push_back(GFdphi_kkdist_container[phi1_id]);

    TVector3 phi2_mom = GFdphi_mom_container[phi2_id];
    TVector3 phi2_vtx = GFdphi_vtx_container[phi2_id];
    event.phimass.push_back(GFdphi_mass_container[phi2_id]);
    event.phimom.push_back(phi2_mom.Mag());
    event.phimom_x.push_back(phi2_mom.x());
    event.phimom_y.push_back(phi2_mom.y());
    event.phimom_z.push_back(phi2_mom.z());
    event.phivtx_x.push_back(phi2_vtx.x());
    event.phivtx_y.push_back(phi2_vtx.y());
    event.phivtx_z.push_back(phi2_vtx.z());
    event.kpkmdist.push_back(GFdphi_kkdist_container[phi2_id]);
    event.phiphidist.push_back(GFdphi_phiphidist_container[best_dphi]);

    TVector3 kp1_mom = GFdphi_kp_mom_container[phi1_id];
    event.phidecay_id.push_back(GFdphi_kp_id_container[phi1_id]);
    event.phidecay_purity.push_back(GFdphi_kp_purity_container[phi1_id]);
    event.phidecay_efficiency.push_back(GFdphi_kp_efficiency_container[phi1_id]);
    event.phidecay_g4tid.push_back(GFdphi_kp_g4tid_container[phi1_id]);
    event.phidecay_charge.push_back(1);
    event.phidecay_mom.push_back(kp1_mom.Mag());
    event.phidecay_mom_x.push_back(kp1_mom.x());
    event.phidecay_mom_y.push_back(kp1_mom.y());
    event.phidecay_mom_z.push_back(kp1_mom.z());

    TVector3 km1_mom = GFdphi_km_mom_container[phi1_id];
    event.phidecay_id.push_back(GFdphi_km_id_container[phi1_id]);
    event.phidecay_purity.push_back(GFdphi_km_purity_container[phi1_id]);
    event.phidecay_efficiency.push_back(GFdphi_km_efficiency_container[phi1_id]);
    event.phidecay_g4tid.push_back(GFdphi_km_g4tid_container[phi1_id]);
    event.phidecay_charge.push_back(-1);
    event.phidecay_mom.push_back(km1_mom.Mag());
    event.phidecay_mom_x.push_back(km1_mom.x());
    event.phidecay_mom_y.push_back(km1_mom.y());
    event.phidecay_mom_z.push_back(km1_mom.z());

    TVector3 kp2_mom = GFdphi_kp_mom_container[phi2_id];
    event.phidecay_id.push_back(GFdphi_kp_id_container[phi2_id]);
    event.phidecay_purity.push_back(GFdphi_kp_purity_container[phi2_id]);
    event.phidecay_efficiency.push_back(GFdphi_kp_efficiency_container[phi2_id]);
    event.phidecay_g4tid.push_back(GFdphi_kp_g4tid_container[phi2_id]);
    event.phidecay_charge.push_back(1);
    event.phidecay_mom.push_back(kp2_mom.Mag());
    event.phidecay_mom_x.push_back(kp2_mom.x());
    event.phidecay_mom_y.push_back(kp2_mom.y());
    event.phidecay_mom_z.push_back(kp2_mom.z());

    TVector3 km2_mom = GFdphi_km_mom_container[phi2_id];
    event.phidecay_id.push_back(GFdphi_km_id_container[phi2_id]);
    event.phidecay_purity.push_back(GFdphi_km_purity_container[phi2_id]);
    event.phidecay_efficiency.push_back(GFdphi_km_efficiency_container[phi2_id]);
    event.phidecay_g4tid.push_back(GFdphi_km_g4tid_container[phi2_id]);
    event.phidecay_charge.push_back(-1);
    event.phidecay_mom.push_back(km2_mom.Mag());
    event.phidecay_mom_x.push_back(km2_mom.x());
    event.phidecay_mom_y.push_back(km2_mom.y());
    event.phidecay_mom_z.push_back(km2_mom.z());
  }
  if(event.nPhi == 1){
    if( event.phidecay_g4tid[0] == event.g4phidecay_g4tid[0] 
    and event.phidecay_g4tid[1] == event.g4phidecay_g4tid[1] ){
      event.good_phi.push_back(1);
    }
    else if( event.phidecay_g4tid[0] == event.g4phidecay_g4tid[2]
    and event.phidecay_g4tid[1] == event.g4phidecay_g4tid[3] ){
      event.good_phi.push_back(1);
    }
    else if( event.phidecay_g4tid[0] == event.g4phidecay_g4tid[0]
    and event.phidecay_g4tid[1] == event.g4phidecay_g4tid[3] ){
      event.swap_phi.push_back(1);
    }
    else if( event.phidecay_g4tid[0] == event.g4phidecay_g4tid[2]
    and event.phidecay_g4tid[1] == event.g4phidecay_g4tid[1] ){
      event.swap_phi.push_back(1);
    }
  }
  if(event.nPhi == 2){//double phi analysis. g4-recon matching is based on kp g4tid.
    event.good_phi.resize(2,0);
    event.swap_phi.resize(2,0);
    if( event.phidecay_g4tid[0] == event.g4phidecay_g4tid[0]
    and event.phidecay_g4tid[1] == event.g4phidecay_g4tid[1] ){
      event.good_phi[0] = 1;
    }
    else if( event.phidecay_g4tid[0] == event.g4phidecay_g4tid[0]
    and event.phidecay_g4tid[1] == event.g4phidecay_g4tid[3] ){
      event.swap_phi[0] = 1;
    }

    if( event.phidecay_g4tid[2] == event.g4phidecay_g4tid[2]
    and event.phidecay_g4tid[3] == event.g4phidecay_g4tid[3] ){
      event.good_phi[1] = 1;
    }
    else if( event.phidecay_g4tid[2] == event.g4phidecay_g4tid[2]
    and event.phidecay_g4tid[3] == event.g4phidecay_g4tid[1] ){
      event.good_phi[1] = 1;
    }


    TVector3 kp1_mom(event.phidecay_mom_x[0], event.phidecay_mom_y[0], event.phidecay_mom_z[0]);
    TVector3 km1_mom(event.phidecay_mom_x[0], event.phidecay_mom_y[0], event.phidecay_mom_z[0]);
    TVector3 kp2_mom(event.phidecay_mom_x[2], event.phidecay_mom_y[2], event.phidecay_mom_z[2]);
    TVector3 km2_mom(event.phidecay_mom_x[2], event.phidecay_mom_y[2], event.phidecay_mom_z[2]);

    double pbeam = event.pK18[0],ubeam = event.utgtK18[0],vbeam = event.vtgtK18[0];
    double nbeam = 1./sqrt(1+ubeam*ubeam+vbeam*vbeam);
    TVector3 beam_mom = TVector3(ubeam*pbeam,vbeam*pbeam,pbeam*nbeam);
    TLorentzVector LVbeam = TLorentzVector(beam_mom, sqrt(beam_mom.Mag2() + ProtonMass*ProtonMass));
    TLorentzVector LVtarget = TLorentzVector(0,0,0,ProtonMass);
    TLorentzVector LVevent = LVbeam + LVtarget;
    
    TLorentzVector LVkp1 = TLorentzVector(kp1_mom, TMath::Hypot(kp1_mom.Mag(), KaonMass));
    TLorentzVector LVkm1 = TLorentzVector(km1_mom, TMath::Hypot(km1_mom.Mag(), KaonMass));
    TLorentzVector LVkp2 = TLorentzVector(kp2_mom, TMath::Hypot(kp2_mom.Mag(), KaonMass));
    TLorentzVector LVkm2 = TLorentzVector(km2_mom, TMath::Hypot(km2_mom.Mag(), KaonMass));

    TLorentzVector LVphi1 = LVkp1 + LVkm1;
    TLorentzVector LVphi2 = LVkp2 + LVkm2;

    event.delM = LVevent.Mag() - (LVphi1 + LVphi2).Mag();
    event.delE = LVevent.E() - (LVphi1 + LVphi2).E();
    event.MissMass = (LVevent - LVphi1 - LVphi2).Mag();
    event.MissMom = LVevent.Vect().Mag() - (LVphi1 + LVphi2).Vect().Mag();
    event.MissMom_x = LVevent.Vect().x() - (LVphi1 + LVphi2).Vect().x();
    event.MissMom_y = LVevent.Vect().y() - (LVphi1 + LVphi2).Vect().y();
    event.MissMom_z = LVevent.Vect().z() - (LVphi1 + LVphi2).Vect().z();
    event.MissMom_t = hypot(event.MissMom_x, event.MissMom_y);

    TLorentzVector LVpip1 = TLorentzVector(kp1_mom, TMath::Hypot(kp1_mom.Mag(), PionMass));
    TLorentzVector LVpim1 = TLorentzVector(km1_mom, TMath::Hypot(km1_mom.Mag(), PionMass));
    TLorentzVector LVpip2 = TLorentzVector(kp2_mom, TMath::Hypot(kp2_mom.Mag(), PionMass));
    TLorentzVector LVpim2 = TLorentzVector(km2_mom, TMath::Hypot(km2_mom.Mag(), PionMass));

    TLorentzVector LVpipi1 = LVpip1 + LVpim1;
    TLorentzVector LVpipi2 = LVpip2 + LVpim2;

    event.delM4Pi = LVpipi1.Mag() + LVpipi2.Mag() - (LVphi1 + LVphi2).Mag();

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
#if SaveHistograms
  /*
  HB1( 1, "Status", 21, 0., 21. );
  HB1( 3, "Genfit Fit Status", 2, 0., 2. );

  HB1( 10, "Missing Mass [KURAMA]; Missing mass (GeV/c^{2}); Counts (/10 MeV/#font[12]{c}^{2})", 320, 1., 1.8);
  HB1( 11," #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/0.002 GeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 12," #Xi Invariant Mass; M_{#Lambda#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/0.002 GeV/#font[12]{c}^{2})", 125, 1.2 ,1.45);
  HB1( 13, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 14, "#pi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 15, "#pi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 16, "#Lambda Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 17, "#Xi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 18, "Closest Dist. p#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 19, "Closest Dist. #Lambda#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);

  HB1( 21, "[GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 22, "[GenFit] #Xi Invariant Mass; M_{#Lambda#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 125, 1.2 ,1.45);
  HB1( 23, "[GenFit] p_{#Lambda} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 24, "[GenFit] #pi_{#Lambda} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 25, "[GenFit] #pi_{#Xi} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 26, "[GenFit] #Lambda Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 27, "[GenFit] #Xi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 28, "[GenFit] Closest Dist. p#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 29, "[GenFit] Closest Dist. #Lambda#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 30, "[GenFit] p_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 31, "[GenFit] pi_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);
  HB1( 32, "[GenFit] pi_{#Xi} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0., 0.5);
  HB1( 33, "[GenFit] p_{#Lambda} TOF_{HTOF mt.} - TOF_{pi calc.}; TOF_{HTOF mt.} - TOF_{pi calc.} (ns);Counts [0.01 ns]", 1000, -5, 5);
  HB1( 34, "[GenFit] pi_{#Lambda} TOF_{HTOF mt.} - TOF_{pi calc.}; TOF_{HTOF mt.} - TOF_{pi calc.} (ns);Counts [0.01 ns]", 1000, -5, 5);
  HB1( 35, "[GenFit] pi_{#Xi} TOF_{HTOF mt.} - TOF_{pi calc.}; TOF_{HTOF mt.} - TOF_{pi calc.} (ns);Counts [0.01 ns]", 1000, -5, 5);
  HB2( 36, "[GenFit] p_{#Lambda} Mass/Charge;Mass/Charge (GeV/#font[12]{c}^{2}/q);Momentum (GeV/#font[12]{c})",1000, -1.0, 1.5, 1000, 0, 2.0);
  HB2( 37, "[GenFit] pi_{#Lambda} Mass/Charge;Mass/Charge (GeV/#font[12]{c}^{2}/q);Momentum (GeV/#font[12]{c})",1000, -1.0, 1.5, 1000, 0, 2.0);
  HB2( 38, "[GenFit] pi_{#Xi} Mass/Charge;Mass/Charge (GeV/#font[12]{c}^{2}/q);Momentum (GeV/#font[12]{c})",1000, -1.0, 1.5, 1000, 0, 2.0);

  HB1( 100, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 101, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 102, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 103, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 110, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 111, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 112, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 113, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 120, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 121, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 122, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 123, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 130, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 131, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 132, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 133, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 134, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 135, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 136, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 137, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 140, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 141, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 142, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);
  HB1( 143, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/10 MeV/#font[12]{c}^{2})", 140, -400., 1000.);

  HB1( 200, ";^{12}C(K^{-}, K^{+}#Xi^{-})X MM - M(^{11}B) (GeV/#font[12]{c}^{2});Counts (/5 MeV/#font[12]{c}^{2})", 80, -0.1, 0.3);

  HB1( 500, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/5 MeV/#font[12]{c}^{2})", 280, -400., 1000.);
  HB2( 501, ";T_{#piL} (MeV/#font[12]{c}); T_{#piH} (MeV/#font[12]{c})", 60, 0, 60, 60, 0, 60);
  HB2( 502, ";M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); ;M_{p#pi^{-}} (GeV/#font[12]{c}^{2})", 80, 1.04, 1.2, 80, 1.04, 1.2);
  HB2( 503, "[GenFit] L Vs L Invariant Mass;M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); ;M_{p#pi^{-}} (GeV/#font[12]{c}^{2})", 80, 1.04, 1.2, 80, 1.04, 1.2);
  HB1( 504, "#Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 505, "[GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);

  HB1( 506, "p_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 507, "pi_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);
  HB1( 508, "[GenFit] p_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 509, "[GenFit] pi_{#Lambda} MassSquare; M^{2} (GeV/#font[12]{c}^{2})^{2}; Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);

  HB1( 1000, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1001, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1002, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1100, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 1101, "#pi^{+} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 1102, "#pi^{-} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);

  HB1( 1010, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1011, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1012, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1110, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 1111, "#pi^{+} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 1112, "#pi^{-} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);

  HB1( 1020, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1021, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1022, "; Multiplicity;Counts ", 5,0,5);

  HB1( 1030, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1031, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1032, "; Multiplicity;Counts ", 5,0,5);
  */
  
  
  
  HB1(10000,"KF#{Lambda} pvalue",100,0,1);
  HB1(10001,"KF#{Lambda} chisqr",1000,0,15);
  HB1(10002,"KF#{Lambda} mass",1000,pdg::LambdaMass()-0.2,pdg::LambdaMass()+0.2);
  HB1(10003,"KF#{Lambda} dist",1000,0,20);
  HB1(10010,"KF#{Lambda} pull p_{p}",100,-5,5);
  HB1(10011,"KF#{Lambda} pull p_{#theta}",100,-5,5);
  HB1(10012,"KF#{Lambda} pull p_{#phi}",100,-5,5);
  HB1(10013,"KF#{Lambda} pull #pi_{p}",100,-5,5);
  HB1(10014,"KF#{Lambda} pull #pi_{#theta}",100,-5,5);
  HB1(10015,"KF#{Lambda} pull #pi_{#phi}",100,-5,5);
  HB1(10016,"KF#{Lambda} pull #Lambda_{p}",100,-5,5);
  HB1(10017,"KF#{Lambda} pull #Lambda_{#theta}",100,-5,5);
  HB1(10018,"KF#{Lambda} pull #Lambda_{#phi}",100,-5,5);

  HB1(10020,"KF#{Lambda} residual p_{p}",1000,-1,1);
  HB1(10021,"KF#{Lambda} residual p_{#theta}",1000,-0.1,0.1);
  HB1(10022,"KF#{Lambda} residual p_{#phi}",1000,-0.1,0.1);
  HB1(10023,"KF#{Lambda} residual #pi_{p}",1000,-0.3,0.3);
  HB1(10024,"KF#{Lambda} residual #pi_{#theta}",1000,-0.1,0.1);
  HB1(10025,"KF#{Lambda} residual #pi_{#phi}",1000,-0.1,0.1);
  HB1(10026,"KF#{Lambda} residual #Lambda_{p}",1000,-0.3,0.3);  
  HB1(10027,"KF#{Lambda} residual #Lambda_{#theta}",1000,-0.1,0.1);
  HB1(10028,"KF#{Lambda} residual #Lambda_{#phi}",1000,-0.3,0.3);
  
  HB1(10030,"#{Lambda} residual p_{p}",1000,-1,1);
  HB1(10031,"#{Lambda} residual p_{#theta}",1000,-0.1,0.1);
  HB1(10032,"#{Lambda} residual p_{#phi}",1000,-0.1,0.1);
  HB1(10033,"#{Lambda} residual #pi_{p}",1000,-0.3,0.3);
  HB1(10034,"#{Lambda} residual #pi_{#theta}",1000,-0.1,0.1);
  HB1(10035,"#{Lambda} residual #pi_{#phi}",1000,-0.1,0.1);
  HB1(10036,"#{Lambda} residual #Lambda_{p}",1000,-0.3,0.3);  
  HB1(10037,"#{Lambda} residual #Lambda_{#theta}",1000,-0.1,0.1);
  HB1(10038,"#{Lambda} residual #Lambda_{#phi}",1000,-0.3,0.3);

  HB1(20000,"KF#{Xi}^{-} pvalue",100,0,1);
  HB1(20001,"KF#{Xi}^{-} chisqr",1000,0,15);
  HB1(20002,"KF#{Xi} mass",1000,pdg::XiMinusMass()-0.2,pdg::XiMinusMass()+0.2);
  HB1(20010,"KF#{Xi}^{-} pull #Lambda_{p}",100,-5,5);
  HB1(20011,"KF#{Xi}^{-} pull #Lambda_{#theta}",100,-5,5);
  HB1(20012,"KF#{Xi}^{-} pull #Lambda_{#phi}",100,-5,5);
  HB1(20013,"KF#{Xi}^{-} pull #pi_{p}",100,-5,5);
  HB1(20014,"KF#{Xi}^{-} pull #pi_{#theta}",100,-5,5);
  HB1(20015,"KF#{Xi}^{-} pull #pi_{#phi}",100,-5,5);

  HB1(20020,"KF#{Xi}^{-} residual #Lambda_{p}",1000,-1,1);
  HB1(20021,"KF#{Xi}^{-} residual #Lambda_{#theta}",1000,-0.5,0.5);
  HB1(20022,"KF#{Xi}^{-} residual #Lambda_{#phi}",1000,-0.5,0.5);
  HB1(20023,"KF#{Xi}^{-} residual #pi_{p}",1000,-0.3,0.3);
  HB1(20024,"KF#{Xi}^{-} residual #pi_{#theta}",1000,-1,1);
  HB1(20025,"KF#{Xi}^{-} residual #pi_{#phi}",1000,-1,1);
  HB1(20026,"KF#{Xi}^{-} residual #Xi_{p}",1000,-0.3,0.3);
  HB1(20027,"KF#{Xi}^{-} residual #Xi_{#theta}",1000,-0.1,0.1);
  HB1(20028,"KF#{Xi}^{-} residual #Xi_{#phi}",1000,-0.3,0.3);
  
  HB1(20030,"#{Xi}^{-} residual #Lambda_{p}",1000,-1,1);
  HB1(20031,"#{Xi}^{-} residual #Lambda_{#theta}",1000,-0.5,0.5);
  HB1(20032,"#{Xi}^{-} residual #Lambda_{#phi}",1000,-0.5,0.5);
  HB1(20033,"#{Xi}^{-} residual #pi_{p}",1000,-0.3,0.3);
  HB1(20034,"#{Xi}^{-} residual #pi_{#theta}",1000,-1,1);
  HB1(20035,"#{Xi}^{-} residual #pi_{#phi}",1000,-1,1);
  HB1(20036,"#{Xi}^{-} residual #Xi_{p}",1000,-0.3,0.3);
  HB1(20037,"#{Xi}^{-} residual #Xi_{#theta}",1000,-0.1,0.1);
  HB1(20038,"#{Xi}^{-} residual #Xi_{#phi}",1000,-0.3,0.3);

  HB1(30000,"p_{p} pull",1000,-5,5);
  HB1(30001,"p_{#theta} pull",1000,-5,5);
  HB1(30002,"p_{#phi} pull",1000,-5,5);
  HB1(30003,"#pi_{1,p} pull",1000,-5,5);
  HB1(30004,"#pi_{1,#theta} pull",1000,-5,5);
  HB1(30005,"#pi_{1,#phi} pull",1000,-5,5);
  HB1(30006,"#pi_{2,p} pull",1000,-5,5);
  HB1(30007,"#pi_{2,#theta} pull",1000,-5,5);
  HB1(30008,"#pi_{2,#phi} pull",1000,-5,5);


#endif
  HBTree( "tpc", "tree of GenfitCarbon" );
  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );

  tree->Branch( "nhHtof", &event.nhHtof );
  tree->Branch( "HtofSeg", &event.HtofSeg );
  tree->Branch( "tHtof", &event.tHtof );
  tree->Branch( "dtHtof", &event.dtHtof );
  tree->Branch( "deHtof", &event.deHtof );
  tree->Branch( "posHtof", &event.posHtof );

  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "cluster_x", &event.cluster_x ); 
  tree->Branch( "cluster_y", &event.cluster_y );
  tree->Branch( "cluster_z", &event.cluster_z );
  tree->Branch( "cluster_de", &event.cluster_de );
  tree->Branch( "cluster_layer", &event.cluster_layer );
  tree->Branch( "cluster_mrow", &event.cluster_mrow );
  tree->Branch( "cluster_houghflag", &event.cluster_houghflag );
  tree->Branch( "cluster_G4tid", &event.cluster_G4tid );
  
  tree->Branch( "remain_nclTpc", &event.remain_nclTpc );
  tree->Branch( "remain_cluster_x", &event.remain_cluster_x ); 
  tree->Branch( "remain_cluster_y", &event.remain_cluster_y );
  tree->Branch( "remain_cluster_z", &event.remain_cluster_z );
  tree->Branch( "remain_cluster_de", &event.remain_cluster_de );
  tree->Branch( "remain_cluster_layer", &event.remain_cluster_layer );
  tree->Branch( "remain_cluster_houghflag", &event.remain_cluster_houghflag );
  tree->Branch( "remain_cluster_G4tid", &event.remain_cluster_G4tid );

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isXi", &event.isXi );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "pid", &event.pid );
  tree->Branch( "purity", &event.purity );
  tree->Branch( "efficiency", &event.efficiency );
  tree->Branch(  "G4tid", &event.G4tid );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "path", &event.path );
#if SaveRawData
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
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
#endif
  tree->Branch( "helix_t", &event.helix_t );



  TTreeReaderCont[kE42] = new TTreeReader( "tpc", TFileCont[kE42] );
  const auto& reader = TTreeReaderCont[kE42];
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );

  TTreeCont[kE42]->SetBranchAddress("nhFtof",&src.nhFtof);
  TTreeCont[kE42]->SetBranchAddress("didFtof",src.didFtof);
  TTreeCont[kE42]->SetBranchAddress("tFtof",src.tFtof);
  TTreeCont[kE42]->SetBranchAddress("deFtof",src.deFtof);
  TTreeCont[kE42]->SetBranchAddress("yFtof",src.yFtof);
  TTreeCont[kE42]->SetBranchAddress("nhittpc",&src.nhittpc);
  TTreeCont[kE42]->SetBranchAddress("ititpc",src.ititpc);
  TTreeCont[kE42]->SetBranchAddress("xtpc",src.xtpc);
  TTreeCont[kE42]->SetBranchAddress("ytpc",src.ytpc);
  TTreeCont[kE42]->SetBranchAddress("ztpc",src.ztpc);
  TTreeCont[kE42]->SetBranchAddress("pxtpc",src.pxtpc);
  TTreeCont[kE42]->SetBranchAddress("pytpc",src.pytpc);
  TTreeCont[kE42]->SetBranchAddress("pztpc",src.pztpc);

  TTreeCont[kE42]->SetBranchAddress("NumberOfTracks",&src.NumberOfTracks);
  TTreeCont[kE42]->SetBranchAddress("PIDOfTrack",src.PIDOfTrack);
  TTreeCont[kE42]->SetBranchAddress("ParentIDOfTrack",src.ParentIDOfTrack);
  TTreeCont[kE42]->SetBranchAddress("VertexOfTrack_x",src.VertexOfTrack_x);
  TTreeCont[kE42]->SetBranchAddress("VertexOfTrack_y",src.VertexOfTrack_y);
  TTreeCont[kE42]->SetBranchAddress("VertexOfTrack_z",src.VertexOfTrack_z);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack",src.MomentumOfTrack);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack_x",src.MomentumOfTrack_x);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack_y",src.MomentumOfTrack_y);
  TTreeCont[kE42]->SetBranchAddress("MomentumOfTrack_z",src.MomentumOfTrack_z);
  
  src.nhHtof = new TTreeReaderValue<Int_t>( *reader, "nhHtof" );
  src.HtofSeg = new TTreeReaderValue<std::vector<Double_t>>( *reader, "HtofSeg" );
  src.tHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tHtof" );
  src.dtHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dtHtof" );
  src.deHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "deHtof" );
  src.posHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "posHtof" );


  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" ); 
  src.cluster_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x" );  
  src.cluster_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y" );
  src.cluster_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z" );
  src.cluster_de = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de" );
  src.cluster_layer = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_layer" );
  src.cluster_mrow = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_mrow" );  
  src.cluster_houghflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_houghflag" );  
  src.cluster_G4tid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_G4tid" );  

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isXi = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isXi" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );
  src.purity = new TTreeReaderValue<std::vector<Double_t>>( *reader, "purity" );
  src.efficiency = new TTreeReaderValue<std::vector<Double_t>>( *reader, "efficiency" );
  src.G4tid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "G4tid" );
  src.chisqr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr" );
  src.pval = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pval" );
  src.helix_cx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx" );
  src.helix_cy = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy" );
  src.helix_z0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0" );
  src.helix_r = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r" );
  src.helix_dz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz" );
  src.dE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dE" );
  src.dEdx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dEdx" );
  src.mom0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0" );
  src.path = new TTreeReaderValue<std::vector<Double_t>>( *reader, "path" );
  src.hitlayer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitlayer" );
  src.hitpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_x" );
  src.hitpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_y" );
  src.hitpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "hitpos_z" );
  src.calpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_x" );
  src.calpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_y" );
  src.calpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "calpos_z" );
  src.residual = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual" );
  src.residual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_x" );
  src.residual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_y" );
  src.residual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_z" );
  src.resolution_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_x" );
  src.resolution_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_y" );
  src.resolution_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_z" );
  src.helix_t = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "helix_t" );
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "alpha" );
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_size" );
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );
  src.xtgtHS = new TTreeReaderValue<vector<Double_t>>( *reader, "xtgtHS" );
  src.ytgtHS = new TTreeReaderValue<vector<Double_t>>( *reader, "ytgtHS" );
  src.xtgtKurama = new TTreeReaderValue<vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<vector<Double_t>>( *reader, "ytgtKurama" );

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
