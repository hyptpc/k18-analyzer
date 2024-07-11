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
#define DebugDisp 1
#define DoKF 1

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
//const Int_t verbosity = 1;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;

//const Double_t xi_masscut = 0.15; const Double_t lambda_masscut = 0.1; //ref
const Double_t xi_masscut = 0.25; const Double_t lambda_masscut = 0.1; //ref
//const Double_t xi_masscut = 0.03; const Double_t lambda_masscut = 0.03;
const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t pi2_vtx_distcut = 300;
const Double_t ppi_distcut = 10.; //ref
const Double_t lpi_distcut = 15.; //ref
const Double_t xitarget_distcut = 50.; //ref
const Double_t ltarget_distcut = 50.;

const Double_t GFppi_distcut = 10.;
const Double_t GFlpi_distcut = 15.;
const Double_t GFxitarget_distcut = 50.;
const Double_t GFltarget_distcut = 50.;
const Double_t GFxitarget_ycut = 20.;
const Double_t GFltarget_ycut = 20.;

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
{ "[Process]", "[ConfFile]", "[DstTPCTrackingHelixGeant4]", "[OutFile]" };
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

  int G4kmid;
  int G4kmtid;
  double G4kmvtx_x;// Production vertex, identical to xi vert.
  double G4kmvtx_y;
  double G4kmvtx_z;
  double G4kmmom;
  double G4kmmom_x;
  double G4kmmom_y;
  double G4kmmom_z;

  int G4kpid;
  int G4kptid;
  double G4kpvtx_x;// Production vertex, identical to xi vert.
  double G4kpvtx_y;
  double G4kpvtx_z;
  double G4kpmom;
  double G4kpmom_x;
  double G4kpmom_y;
  double G4kpmom_z;

  int G4xiid;
  double G4xivtx_x;//Production vertex
  double G4xivtx_y;
  double G4xivtx_z;
  double G4ximom;
  double G4ximom_x;//Momentum at production vtx
  double G4ximom_y;
  double G4ximom_z;

  double xivtx_x;
  double xivtx_y;
  double xivtx_z;
  double ximom;

  int G4lid;
  double G4lvtx_x;//Production vertex, identical to xi decay vert
  double G4lvtx_y;
  double G4lvtx_z;
  double G4lmom;
  double G4lmom_x;
  double G4lmom_y;
  double G4lmom_z;
  double lvtx_x;
  double lvtx_y;
  double lvtx_z;
  double lmom;

  int G4pid;//Id of Geant4 Track
  int G4ptid;//G4 Id of TPC Track
  int G4pnh;//Number of G4 Hits
  int G4ptnh;//Number of correct G4 Hits
  double G4pvtx_x;//Production vertex, identical to l decay vert
  double G4pvtx_y;
  double G4pvtx_z;
  double G4pmom;
  double G4pmom_x;
  double G4pmom_y;
  double G4pmom_z;

  int ptid;
  int pnh;
  double pvtx_x;
  double pvtx_y;
  double pvtx_z;
  double pmom;
  double pmom_x;
  double pmom_y;
  double pmom_z;
  double GFpmom;
  double GFpmom_x;
  double GFpmom_y;
  double GFpmom_z;


  int G4pi1id;
  int G4pi1tid;
  int G4pi1nh;
  int G4pi1tnh;
  double G4pi1vtx_x;
  double G4pi1vtx_y;
  double G4pi1vtx_z;
  double G4pi1mom;
  double G4pi1mom_x;
  double G4pi1mom_y;
  double G4pi1mom_z;

  int pi1tid;
  int pi1nh;
  double pi1vtx_x;
  double pi1vtx_y;
  double pi1vtx_z;
  double pi1mom;
  double pi1mom_x;
  double pi1mom_y;
  double pi1mom_z;
  double GFpi1mom;
  double GFpi1mom_x;
  double GFpi1mom_y;
  double GFpi1mom_z;


  int G4pi2id;
  int G4pi2tid;
  int G4pi2nh;
  int G4pi2tnh;
  double G4pi2vtx_x;
  double G4pi2vtx_y;
  double G4pi2vtx_z;
  double G4pi2mom;
  double G4pi2mom_x;
  double G4pi2mom_y;
  double G4pi2mom_z;

  int pi2tid;
  int pi2nh;
  double pi2vtx_x;
  double pi2vtx_y;
  double pi2vtx_z;
  double pi2mom;
  double pi2mom_x;
  double pi2mom_y;
  double pi2mom_z;
  double GFpi2mom;
  double GFpi2mom_x;
  double GFpi2mom_y;
  double GFpi2mom_z;


  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> trackid; //for Kurama K1.8 tracks
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
  std::vector<Int_t> G4tid;
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
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_mrow;

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

  //Multi-track production vertex
  Double_t GFprodvtx_x;
  Double_t GFprodvtx_y;
  Double_t GFprodvtx_z;

  Bool_t lpiflag;
  Bool_t lflag;
//
  Bool_t llflag;
  Double_t ltarget_dist1;
  Double_t ltargetvtx_x1;
  Double_t ltargetvtx_y1;
  Double_t ltargetvtx_z1;
  Double_t lmass1;
  Double_t ldecayvtx_x1;
  Double_t ldecayvtx_y1;
  Double_t ldecayvtx_z1;
  Double_t lmom1;
  Double_t lmom_x1;
  Double_t lmom_y1;
  Double_t lmom_z1;
  Double_t ppi_dist1;
  Double_t ltarget_dist2;
  Double_t ltargetvtx_x2;
  Double_t ltargetvtx_y2;
  Double_t ltargetvtx_z2;
  Double_t lmass2;
  Double_t ldecayvtx_x2;
  Double_t ldecayvtx_y2;
  Double_t ldecayvtx_z2;
  Double_t lmom2;
  Double_t lmom_x2;
  Double_t lmom_y2;
  Double_t lmom_z2;
  Double_t ppi_dist2;

  Double_t GFlmass1;
  Double_t GFldecayvtx_x1;
  Double_t GFldecayvtx_y1;
  Double_t GFldecayvtx_z1;
  Double_t GFlmom1;
  Double_t GFlmom_x1;
  Double_t GFlmom_y1;
  Double_t GFlmom_z1;
  Double_t GFppi_dist1;
  Double_t GFltarget_dist1;
  Double_t GFltargetvtx_x1;
  Double_t GFltargetvtx_y1;
  Double_t GFltargetvtx_z1;
  Double_t GFltargetcenter_x1;
  Double_t GFltargetcenter_y1;
  Double_t GFltargetcenter_z1;
  Double_t GFltargetcenter_dist1;
  Double_t GFlprodvtx_x1;
  Double_t GFlprodvtx_y1;
  Double_t GFlprodvtx_z1;
  Double_t GFlprodvtx_dist1;
  Double_t GFltracklen1;
  Double_t GFltof1;

  Double_t GFlmass2;
  Double_t GFldecayvtx_x2;
  Double_t GFldecayvtx_y2;
  Double_t GFldecayvtx_z2;
  Double_t GFlmom2;
  Double_t GFlmom_x2;
  Double_t GFlmom_y2;
  Double_t GFlmom_z2;
  Double_t GFppi_dist2;
  Double_t GFltarget_dist2;
  Double_t GFltargetvtx_x2;
  Double_t GFltargetvtx_y2;
  Double_t GFltargetvtx_z2;
  Double_t GFltargetcenter_x2;
  Double_t GFltargetcenter_y2;
  Double_t GFltargetcenter_z2;
  Double_t GFltargetcenter_dist2;
  Double_t GFlprodvtx_x2;
  Double_t GFlprodvtx_y2;
  Double_t GFlprodvtx_z2;
  Double_t GFlprodvtx_dist2;
  Double_t GFltracklen2;
  Double_t GFltof2;

  Double_t llvtx_x;
  Double_t llvtx_y;
  Double_t llvtx_z;
  Double_t lldist;
  Double_t GFllvtx_x;
  Double_t GFllvtx_y;
  Double_t GFllvtx_z;
  Double_t GFlldist;
//
  int best;
  Bool_t xiflag;
  Double_t ximass;
  Double_t xidecayvtx_x;
  Double_t xidecayvtx_y;
  Double_t xidecayvtx_z;
  Double_t ximom;
  Double_t ximom_x;
  Double_t ximom_y;
  Double_t ximom_z;
  Double_t lpi_dist;
  Double_t xitargetvtx_x;
  Double_t xitargetvtx_y;
  Double_t xitargetvtx_z;
  Double_t xitargetmom;
  Double_t xitargetmom_x;
  Double_t xitargetmom_y;
  Double_t xitargetmom_z;
  Double_t xitarget_dist;

  Double_t lmass;
  Double_t ldecayvtx_x;
  Double_t ldecayvtx_y;
  Double_t ldecayvtx_z;
  Double_t lmom;
  Double_t lmom_x;
  Double_t lmom_y;
  Double_t lmom_z;
  Double_t ppi_dist;
  Double_t ltarget_dist;
  Double_t ltargetvtx_x;
  Double_t ltargetvtx_y;
  Double_t ltargetvtx_z;
//???
  Double_t lmass_vtx;
  Double_t ldecayvtx_x_vtx;
  Double_t ldecayvtx_y_vtx;
  Double_t ldecayvtx_z_vtx;
  Double_t lmom_vtx;
  Double_t lmom_x_vtx;
  Double_t lmom_y_vtx;
  Double_t lmom_z_vtx;
  Double_t ppi_dist_vtx;
//
  
  Bool_t GFxiflag;
  Double_t GFximass;
  Double_t GFxidecayvtx_x;
  Double_t GFxidecayvtx_y;
  Double_t GFxidecayvtx_z;
  Double_t GFximom;
  Double_t GFximom_x;
  Double_t GFximom_y;
  Double_t GFximom_z;
  

  Double_t GFxikkvtx_x;
  Double_t GFxikkvtx_y;
  Double_t GFxikkvtx_z;
  Double_t GFxikkmom;
  Double_t GFxikkmom_x;
  Double_t GFxikkmom_y;
  Double_t GFxikkmom_z;
  Double_t GFxikkvtx_dist;

  Double_t GFxiprodvtx_x;
  Double_t GFxiprodvtx_y;
  Double_t GFxiprodvtx_z;
  Double_t GFxiprodmom;
  Double_t GFxiprodmom_x;
  Double_t GFxiprodmom_y;
  Double_t GFxiprodmom_z;
  Double_t GFxiprodvtx_dist;
  Double_t GFxitracklen;
  Double_t GFxitof;
  Double_t GFlpi_dist;
  Double_t GFximomloss;
  Double_t GFxiexcitation;

  Double_t GFxitargetvtx_x;
  Double_t GFxitargetvtx_y;
  Double_t GFxitargetvtx_z;
  Double_t GFxitargetmom;
  Double_t GFxitargetmom_x;
  Double_t GFxitargetmom_y;
  Double_t GFxitargetmom_z;
  Double_t GFxitarget_dist;

  Double_t GFxitargetcenter_x;
  Double_t GFxitargetcenter_y;
  Double_t GFxitargetcenter_z;
  Double_t GFxitargetcentermom;
  Double_t GFxitargetcentermom_x;
  Double_t GFxitargetcentermom_y;
  Double_t GFxitargetcentermom_z;
  Double_t GFxitargetcenter_dist;

  Double_t GFlmass;
  Double_t GFldecayvtx_x;
  Double_t GFldecayvtx_y;
  Double_t GFldecayvtx_z;
  Double_t GFlmom;
  Double_t GFlmom_x;
  Double_t GFlmom_y;
  Double_t GFlmom_z;
  Double_t GFppi_dist;
  Double_t GFltarget_dist;
  Double_t GFltargetvtx_x;
  Double_t GFltargetvtx_y;
  Double_t GFltargetvtx_z;
  Double_t GFltargetcenter_x;
  Double_t GFltargetcenter_y;
  Double_t GFltargetcenter_z;
  Double_t GFltargetcenter_dist;
  Double_t GFlprodvtx_x;
  Double_t GFlprodvtx_y;
  Double_t GFlprodvtx_z;
  Double_t GFlprodvtx_dist;
  Double_t GFltracklen;
  Double_t GFltof;

  Bool_t lphiflag;
  Double_t phimass;
  Double_t phidecayvtx_x;
  Double_t phidecayvtx_y;
  Double_t phidecayvtx_z;
  Double_t phimom;
  Double_t phimom_x;
  Double_t phimom_y;
  Double_t phimom_z;
  Double_t kk_dist;
  Double_t phi_km_mass2;

  Double_t GFphimass;
  Double_t GFphidecayvtx_x;
  Double_t GFphidecayvtx_y;
  Double_t GFphidecayvtx_z;
  Double_t GFphimom;
  Double_t GFphimom_x;
  Double_t GFphimom_y;
  Double_t GFphimom_z;
  Double_t GFkk_dist;
  Double_t GFphiprodvtx_dist;
  std::vector<Int_t> phidecays_id;
  std::vector<Double_t> phidecays_mom;
  std::vector<Double_t> phidecays_mom_x;
  std::vector<Double_t> phidecays_mom_y;
  std::vector<Double_t> phidecays_mom_z;
  std::vector<Double_t> GFphidecays_mom;
  std::vector<Double_t> GFphidecays_mom_x;
  std::vector<Double_t> GFphidecays_mom_y;
  std::vector<Double_t> GFphidecays_mom_z;

  Int_t GFntdecays;
  std::vector<Int_t> GFdecays_htofid;
  std::vector<Double_t> GFdecays_tracklen;
  std::vector<Double_t> GFdecays_tof;
  std::vector<Double_t> GFdecays_mass2;
  std::vector<Double_t> GFdecays_mom;
  std::vector<Double_t> GFdecays_mom_x;
  std::vector<Double_t> GFdecays_mom_y;
  std::vector<Double_t> GFdecays_mom_z;
  std::vector<Double_t> GFmomloss;
  std::vector<Double_t> GFeloss;

  Int_t GFntTpc;
  std::vector<Int_t> GFpdgcode;
  std::vector<Int_t> GFnhtrack;
  std::vector<Double_t> GFcharge;
  std::vector<Double_t> GFchisqr;
  std::vector<Double_t> GFtof;
  std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFpval;
  
  std::vector<Int_t> decays_id;
  std::vector<Double_t> decays_purity;
  std::vector<Double_t> decays_efficiency;
  std::vector<Int_t> decays_G4tid;
  std::vector<Double_t> decays_mom;
  std::vector<Double_t> decays_mom_x;
  std::vector<Double_t> decays_mom_y;
  std::vector<Double_t> decays_mom_z;


  Bool_t pipiflag;

  Int_t residual_multi;
  Int_t pim_multi;
  Int_t pip_multi;
  Int_t p_multi;
  Int_t ppip_multi;
  Int_t accdient_multi;

  std::vector<Int_t> residual_id;
  std::vector<Double_t> residual_mass2;
  std::vector<Double_t> residual_mom;
  std::vector<Double_t> residual_mom_x;
  std::vector<Double_t> residual_mom_y;
  std::vector<Double_t> residual_mom_z;
  std::vector<Double_t> residual_charge;


  Double_t KFlchisqr;
  Double_t KFlpval;
  std::vector<Double_t> KFlpull;
  Double_t KFlmass;//MassBefore;
  Double_t KFlmom;
  Double_t KFlmom_x;
  Double_t KFlmom_y;
  Double_t KFlmom_z;
  Double_t KFlpi_dist;
  
  Double_t KFxichisqr;
  Double_t KFxipval;
  std::vector<Double_t> KFxipull;
  Double_t KFximass;
  Double_t KFximom_x;
  Double_t KFximom_y;
  Double_t KFximom_z;
  
  vector<Double_t> KFdecays_mom;
  vector<Double_t> KFdecays_mom_x;
  vector<Double_t> KFdecays_mom_y;
  vector<Double_t> KFdecays_mom_z;

  void clear( void )
  {
    status = 0;
    runnum = 0;
    evnum = 0;
    trigpat.clear();
    trigflag.clear();
    
    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    deHtof.clear();
    posHtof.clear();

    nhFtof = 0;
    FtofSeg.clear();
    tFtof.clear();
    deFtof.clear();
    posFtof.clear();
    
    G4kmid = -1;
    G4kmtid = -1;
    G4kmvtx_x = qnan;
    G4kmvtx_y = qnan;
    G4kmvtx_z = qnan;
    G4kmmom = qnan;
    G4kmmom_x = qnan;
    G4kmmom_y = qnan;
    G4kmmom_z = qnan;
    
    G4kpid = -1;
    G4kptid = -1;
    G4kpvtx_x = qnan;
    G4kpvtx_y = qnan;
    G4kpvtx_z = qnan;
    G4kpmom = qnan;
    G4kpmom_x = qnan;
    G4kpmom_y = qnan;
    G4kpmom_z = qnan;
    
    
    G4xiid = -1;
    G4xivtx_x = qnan;
    G4xivtx_y = qnan;
    G4xivtx_z = qnan;
    G4ximom = qnan;
    G4ximom_x = qnan;
    G4ximom_y = qnan;
    G4ximom_z = qnan;
    xivtx_x = qnan;
    xivtx_y = qnan;
    xivtx_z = qnan;
    ximom = qnan;

    G4lid = -1;
    G4lvtx_x = qnan;
    G4lvtx_y = qnan;
    G4lvtx_z = qnan;
    G4lmom = qnan;
    G4lmom_x = qnan;
    G4lmom_y = qnan;
    G4lmom_z = qnan;
    lvtx_x = qnan;
    lvtx_y = qnan;
    lvtx_z = qnan;
    lmom = qnan;
    
    G4pid = -1;
    G4ptid = -1;
    G4pnh = -1;
    G4ptnh = -1;
    G4pvtx_x = qnan;
    G4pvtx_y = qnan;
    G4pvtx_z = qnan;
    G4pmom = qnan;
    G4pmom_x = qnan;
    G4pmom_y = qnan;
    G4pmom_z = qnan;
    ptid = -1;
    pnh = -1;
    pvtx_x = qnan;
    pvtx_y = qnan;
    pvtx_z = qnan;
    pmom = qnan;
    pmom_x = qnan;
    pmom_y = qnan;
    pmom_z = qnan;
    GFpmom = qnan;
    GFpmom_x = qnan;
    GFpmom_y = qnan;
    GFpmom_z = qnan;

    G4pi1id = -1;
    G4pi1tid = -1;
    G4pi1nh = -1;
    G4pi1tnh = -1;
    G4pi1vtx_x = qnan;
    G4pi1vtx_y = qnan;
    G4pi1vtx_z = qnan;
    G4pi1mom = qnan;
    G4pi1mom_x = qnan;
    G4pi1mom_y = qnan;
    G4pi1mom_z = qnan;
    pi1tid = -1;
    pi1nh = -1;
    pi1vtx_x = qnan;
    pi1vtx_y = qnan;
    pi1vtx_z = qnan;
    pi1mom = qnan;
    pi1mom_x = qnan;
    pi1mom_y = qnan;
    pi1mom_z = qnan;
    pi1res_x = qnan;
    pi1res_y = qnan;
    pi1res_z = qnan;
    pi1cov_xy = qnan;
    pi1cov_yz = qnan;
    pi1cov_zx = qnan;
    GFpi1mom = qnan;
    GFpi1mom_x = qnan;
    GFpi1mom_y = qnan;
    GFpi1mom_z = qnan;

    G4pi2id = -1;
    G4pi2tid = -1;
    G4pi2nh = -1;
    G4pi2tnh = -1;
    G4pi2vtx_x = qnan;
    G4pi2vtx_y = qnan;
    G4pi2vtx_z = qnan;
    G4pi2mom = qnan;
    G4pi2mom_x = qnan;
    G4pi2mom_y = qnan;
    G4pi2mom_z = qnan;
    pi2tid = -1;
    pi2nh = -1;
    pi2vtx_x = qnan;
    pi2vtx_y = qnan;
    pi2vtx_z = qnan;
    pi2mom = qnan;
    pi2mom_x = qnan;
    pi2mom_y = qnan;
    pi2mom_z = qnan;
    pi2res_x = qnan;
    pi2res_y = qnan;
    pi2res_z = qnan;
    pi2cov_xy = qnan;
    pi2cov_yz = qnan;
    pi2cov_zx = qnan;
    GFpi2mom = qnan;
    GFpi2mom_x = qnan;
    GFpi2mom_y = qnan;
    GFpi2mom_z = qnan;

    nhHtof = 0;
    HtofSeg.clear();
    tHtof.clear();
    dtHtof.clear();
    deHtof.clear();
    posHtof.clear();

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    charge.clear();
    pid.clear();

    chisqr.clear();
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
    track_cluster_de.clear();
    track_cluster_mrow.clear();

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

    GFprodvtx_x = qnan;
    GFprodvtx_y = qnan;
    GFprodvtx_z = qnan;

    GFntTpc = 0;
    GFcharge.clear();
    GFchisqr.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFpdgcode.clear();
    GFnhtrack.clear();

    GFntdecays = 0;
    GFdecays_htofid.clear();
    GFdecays_tracklen.clear();
    GFdecays_tof.clear();
    GFdecays_mass2.clear();
    GFdecays_mom.clear();
    GFdecays_mom_x.clear();
    GFdecays_mom_y.clear();
    GFdecays_mom_z.clear();
    GFmomloss.clear();
    GFeloss.clear();

    xiflag = false;
    ximass = qnan;
    xidecayvtx_x = qnan;
    xidecayvtx_y = qnan;
    xidecayvtx_z = qnan;
    ximom = qnan;
    ximom_x = qnan;
    ximom_y = qnan;
    ximom_z = qnan;
    lpi_dist = qnan;

    lmass = qnan;
    ldecayvtx_x = qnan;
    ldecayvtx_y = qnan;
    ldecayvtx_z = qnan;
    lmom = qnan;
    lmom_x = qnan;
    lmom_y = qnan;
    lmom_z = qnan;
    ppi_dist = qnan;
    ltarget_dist = qnan;
    ltargetvtx_x = qnan;
    ltargetvtx_y = qnan;
    ltargetvtx_z = qnan;

    lmass_vtx = qnan;
    ldecayvtx_x_vtx = qnan;
    ldecayvtx_y_vtx = qnan;
    ldecayvtx_z_vtx = qnan;
    lmom_vtx = qnan;
    lmom_x_vtx = qnan;
    lmom_y_vtx = qnan;
    lmom_z_vtx = qnan;
    ppi_dist_vtx = qnan;

    xitargetvtx_x = qnan;
    xitargetvtx_y = qnan;
    xitargetvtx_z = qnan;
    xitargetmom = qnan;
    xitargetmom_x = qnan;
    xitargetmom_y = qnan;
    xitargetmom_z = qnan;
    xitarget_dist = qnan;

    GFximass = qnan;
    GFxidecayvtx_x = qnan;
    GFxidecayvtx_y = qnan;
    GFxidecayvtx_z = qnan;
    GFximom = qnan;
    GFximom_x = qnan;
    GFximom_y = qnan;
    GFximom_z = qnan;
    GFlpi_dist = qnan;

    GFlmass = qnan;
    GFldecayvtx_x = qnan;
    GFldecayvtx_y = qnan;
    GFldecayvtx_z = qnan;
    GFlmom = qnan;
    GFlmom_x = qnan;
    GFlmom_y = qnan;
    GFlmom_z = qnan;
    GFltracklen = qnan;
    GFltof = qnan;
    GFppi_dist = qnan;
    GFltarget_dist = qnan;
    GFltargetvtx_x = qnan;
    GFltargetvtx_y = qnan;
    GFltargetvtx_z = qnan;
    GFltargetcenter_dist = qnan;
    GFltargetcenter_x = qnan;
    GFltargetcenter_y = qnan;
    GFltargetcenter_z = qnan;
    GFlprodvtx_x = qnan;
    GFlprodvtx_y = qnan;
    GFlprodvtx_z = qnan;
    GFlprodvtx_dist = qnan;
    GFltracklen = qnan;
    GFltof = qnan;

    GFxitargetvtx_x = qnan;
    GFxitargetvtx_y = qnan;
    GFxitargetvtx_z = qnan;
    GFxitargetmom = qnan;
    GFxitargetmom_x = qnan;
    GFxitargetmom_y = qnan;
    GFxitargetmom_z = qnan;
    GFxitarget_dist = qnan;

    GFxitargetcenter_x = qnan;
    GFxitargetcenter_y = qnan;
    GFxitargetcenter_z = qnan;
    GFxitargetcentermom = qnan;
    GFxitargetcentermom_x = qnan;
    GFxitargetcentermom_y = qnan;
    GFxitargetcentermom_z = qnan;
    GFxitargetcenter_dist = qnan;

    GFxikkvtx_x = qnan;
    GFxikkvtx_y = qnan;
    GFxikkvtx_z = qnan;
    GFxikkmom = qnan;
    GFxikkmom_x = qnan;
    GFxikkmom_y = qnan;
    GFxikkmom_z = qnan;
    GFxikkvtx_dist = qnan;

    GFxiprodvtx_x = qnan;
    GFxiprodvtx_y = qnan;
    GFxiprodvtx_z = qnan;
    GFxiprodmom = qnan;
    GFxiprodmom_x = qnan;
    GFxiprodmom_y = qnan;
    GFxiprodmom_z = qnan;
    GFxiprodvtx_dist = qnan;
    GFxitracklen = qnan;
    GFxitof = qnan;
    GFximomloss = qnan;
    GFxiexcitation = qnan;

    lphiflag = false;
    phimass  = qnan;
    phidecayvtx_x  = qnan;
    phidecayvtx_y = qnan;
    phidecayvtx_z = qnan;
    phimom = qnan;
    phimom_x = qnan;
    phimom_y = qnan;
    phimom_z = qnan;
    kk_dist = qnan;
    phi_km_mass2 = qnan;
    GFphimass = qnan;
    GFphidecayvtx_x = qnan;
    GFphidecayvtx_y = qnan;
    GFphidecayvtx_z = qnan;
    GFphimom = qnan;
    GFphimom_x = qnan;
    GFphimom_y = qnan;
    GFphimom_z = qnan;
    GFkk_dist = qnan;
    GFphiprodvtx_dist = qnan;

    phidecays_id.clear();;
    phidecays_mom.clear();;
    phidecays_mom_x.clear();;
    phidecays_mom_y.clear();;
    phidecays_mom_z.clear();;
    GFphidecays_mom.clear();;
    GFphidecays_mom_x.clear();;
    GFphidecays_mom_y.clear();;
    GFphidecays_mom_z.clear();;

    lpiflag = false;
    lflag = false;
    llflag = false;
    pipiflag = false;
    lmass1 = qnan;
    ltarget_dist1 = qnan;
    ltargetvtx_x1 = qnan;
    ltargetvtx_y1 = qnan;
    ltargetvtx_z1 = qnan;
    ldecayvtx_x1 = qnan;
    ldecayvtx_y1 = qnan;
    ldecayvtx_z1 = qnan;
    lmom1 = qnan;
    lmom_x1 = qnan;
    lmom_y1 = qnan;
    lmom_z1 = qnan;
    ppi_dist1 = qnan;
    lmass2 = qnan;
    ltarget_dist2 = qnan;
    ltargetvtx_x2 = qnan;
    ltargetvtx_y2 = qnan;
    ltargetvtx_z2 = qnan;

    ldecayvtx_x2 = qnan;
    ldecayvtx_y2 = qnan;
    ldecayvtx_z2 = qnan;
    lmom2 = qnan;
    lmom_x2 = qnan;
    lmom_y2 = qnan;
    lmom_z2 = qnan;
    ppi_dist2 = qnan;

    GFlmass1 = qnan;
    GFldecayvtx_x1 = qnan;
    GFldecayvtx_y1 = qnan;
    GFldecayvtx_z1 = qnan;
    GFlmom1 = qnan;
    GFlmom_x1 = qnan;
    GFlmom_y1 = qnan;
    GFlmom_z1 = qnan;
    GFppi_dist1 = qnan;
    GFltarget_dist1 = qnan;
    GFltargetvtx_x1 = qnan;
    GFltargetvtx_y1 = qnan;
    GFltargetvtx_z1 = qnan;
    GFltargetcenter_dist1 = qnan;
    GFltargetcenter_x1 = qnan;
    GFltargetcenter_y1 = qnan;
    GFltargetcenter_z1 = qnan;
    GFlprodvtx_x1 = qnan;
    GFlprodvtx_y1 = qnan;
    GFlprodvtx_z1 = qnan;
    GFlprodvtx_dist1 = qnan;
    GFltracklen1 = qnan;
    GFltof1 = qnan;

    GFlmass2 = qnan;
    GFldecayvtx_x2 = qnan;
    GFldecayvtx_y2 = qnan;
    GFldecayvtx_z2 = qnan;
    GFlmom2 = qnan;
    GFlmom_x2 = qnan;
    GFlmom_y2 = qnan;
    GFlmom_z2 = qnan;
    GFppi_dist2 = qnan;
    GFltarget_dist2 = qnan;
    GFltargetvtx_x2 = qnan;
    GFltargetvtx_y2 = qnan;
    GFltargetvtx_z2 = qnan;
    GFltargetcenter_dist2 = qnan;
    GFltargetcenter_x2 = qnan;
    GFltargetcenter_y2 = qnan;
    GFltargetcenter_z2 = qnan;
    GFlprodvtx_x2 = qnan;
    GFlprodvtx_y2 = qnan;
    GFlprodvtx_z2 = qnan;
    GFlprodvtx_dist2 = qnan;
    GFltracklen2 = qnan;
    GFltof2 = qnan;

    llvtx_x = qnan;
    llvtx_y = qnan;
    llvtx_z = qnan;
    lldist = qnan;
    GFllvtx_x = qnan;
    GFllvtx_y = qnan;
    GFllvtx_z = qnan;
    GFlldist = qnan;

    decays_id.clear();
    decays_mom.clear();
    decays_mom_x.clear();
    decays_mom_y.clear();
    decays_mom_z.clear();

    residual_multi = 0;
    pim_multi = 0;
    pip_multi = 0;
    p_multi = 0;
    ppip_multi = 0;
    accdient_multi = 0;

    residual_id.clear();
    residual_mass2.clear();
    residual_mom.clear();
    residual_mom_x.clear();
    residual_mom_y.clear();
    residual_mom_z.clear();
    residual_charge.clear();

    KFlchisqr = qnan;
    KFlpval = qnan;
    KFxichisqr = qnan;
    KFxipval = qnan;

    KFximass = qnan;
    KFlpi_dist = qnan;
    KFlpull.clear();
    KFxipull.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* evnum;
  
  int nhHtof;
  int didHtof[500];
  double tHtof[500];
  double deHtof[500];
  double yHtof[500];//->posHtof

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
  
  

  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Int_t>>* isMultiloop;
  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* mom_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* residual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* resolution_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* helix_t;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* alpha;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;


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
  //skip = 456640; max_loop = 1;
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

  static const auto PionMass    = pdg::PionMass();
  static const auto KaonMass    = pdg::KaonMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  static const auto m12C = 11.174864;
  static const auto m11B = 10.252548;
  static const auto me = 0.001*0.5109989461;
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};

  TVector3 tgtpos(0, 0, tpc::ZTarget);

  //if( ievent%100000==0 ){
  //if( ievent%1000==0 ){
  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  event.evnum = **src.evnum;

  event.nhHtof = src.nhHtof;
  for(int ih=0;ih<event.nhHtof;++ih){
    event.HtofSeg.push_back(src.didHtof[ih]);
    event.tHtof.push_back(src.tHtof[ih]);
    event.deHtof.push_back(src.deHtof[ih]);
    event.posHtof.push_back(src.yHtof[ih]);
  }

  event.nhFtof = src.nhFtof;
  for(int ih=0;ih<event.nhFtof;++ih){
    event.FtofSeg.push_back(src.didFtof[ih]);
    event.tFtof.push_back(src.tFtof[ih]);
    event.deFtof.push_back(src.deFtof[ih]);
    event.posFtof.push_back(src.yFtof[ih]);
  }


  Int_t ntTpc = **src.ntTpc;
  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
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
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_mrow = **src.track_cluster_mrow;
	event.xtgtHS = **src.xtgtHS;
	event.ytgtHS = **src.ytgtHS;
	event.xtgtKurama = **src.xtgtKurama;
	event.ytgtKurama = **src.ytgtKurama;

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
    double mom = src.MomentumOfTrack[it]/1000;//MeV to GeV
    double mom_x = src.MomentumOfTrack_x[it]/1000;
    double mom_y = src.MomentumOfTrack_y[it]/1000;
    double mom_z = src.MomentumOfTrack_z[it]/1000;
    double vert_x = src.VertexOfTrack_x[it];
    double vert_y = src.VertexOfTrack_y[it];
    double vert_z = src.VertexOfTrack_z[it];
    if(abs(pid) == 321 and parent ==0){
      if(mom_z < 0){
        event.G4kmid = it;
        event.G4kmvtx_x = vert_x;
        event.G4kmvtx_y = vert_y;
        event.G4kmvtx_z = vert_z;
        event.G4kmmom = mom;
        event.G4kmmom_x = -mom_x;
        event.G4kmmom_y = -mom_y;
        event.G4kmmom_z = -mom_z;
      }
      else {
        event.G4kpid = it;
        event.G4kpvtx_x = vert_x;
        event.G4kpvtx_y = vert_y;
        event.G4kpvtx_z = vert_z;
        event.G4kpmom = mom;
        event.G4kpmom_x = mom_x;
        event.G4kpmom_y = mom_y;
        event.G4kpmom_z = mom_z;
      }
    }
    if(pid==3312){//Xi
        event.G4xiid = it;
        event.G4xivtx_x = vert_x;
        event.G4xivtx_y = vert_y;
        event.G4xivtx_z = vert_z;
        event.G4ximom = mom;
        event.G4ximom_x = mom_x;
        event.G4ximom_y = mom_y;
        event.G4ximom_z = mom_z;
    }
    if(pid==3122 and src.PIDOfTrack[parent] == 3312){//Lambda
        event.G4lid = it;
        event.G4lvtx_x = vert_x;
        event.G4lvtx_y = vert_y;
        event.G4lvtx_z = vert_z;
        event.G4lmom = mom;
        event.G4lmom_x = mom_x;
        event.G4lmom_y = mom_y;
        event.G4lmom_z = mom_z;
    }
    if(pid == 2212 and src.PIDOfTrack[parent] == 3122){//Proton
        event.G4pid = it;
        event.G4pvtx_x = vert_x;
        event.G4pvtx_y = vert_y;
        event.G4pvtx_z = vert_z;
        event.G4pmom = mom;
        event.G4pmom_x = mom_x;
        event.G4pmom_y = mom_y;
        event.G4pmom_z = mom_z;
    }
    if(abs(pid) == 211){
      if(src.PIDOfTrack[parent] == 3312){//Pi form Xi
        event.G4pi2id = it;
        event.G4pi2vtx_x = vert_x;
        event.G4pi2vtx_y = vert_y;
        event.G4pi2vtx_z = vert_z;
        event.G4pi2mom = mom;
        event.G4pi2mom_x = mom_x;
        event.G4pi2mom_y = mom_y;
        event.G4pi2mom_z = mom_z;
      }
      else if(src.PIDOfTrack[parent] == 3122){//Pi from L
        event.G4pi1id = it;
        event.G4pi1vtx_x = vert_x;
        event.G4pi1vtx_y = vert_y;
        event.G4pi1vtx_z = vert_z;
        event.G4pi1mom = mom;
        event.G4pi1mom_x = mom_x;
        event.G4pi1mom_y = mom_y;
        event.G4pi1mom_z = mom_z;
      }
    }
  }



  HF1( 1, event.status++ );

  if( ntTpc == 0 ) return true;
    std::cout << "#D Event Number: "
        << std::setw(6) << ievent << std::endl;
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
  }
  GFTrackCont.FitTracks();
  HF1( 1, event.status++ );

  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }

#if DebugDisp
  std::cout<<"0. Before particle reconstruction, checking each track is originated from the target or not."<<std::endl;
#endif

  std::vector<Int_t> target_accidental_id_container;
  std::vector<Int_t> target_k_id_container;
  std::vector<Double_t> target_k_mass2_container;
  std::vector<TVector3> target_k_mom_container;
  std::vector<Int_t> target_p_id_container, target_pip_id_container, target_pim_id_container, target_ppip_id_container;
  std::vector<Double_t> target_p_mass2_container, target_pip_mass2_container, target_pim_mass2_container, target_ppip_mass2_container;
  std::vector<TVector3> target_p_mom_container, target_pip_mom_container, target_pim_mom_container, target_ppip_mom_container;
  for(Int_t it=0;it<ntTpc;it++){ //1st searching
    if(event.isK18[it]==1) continue;
    if(event.isKurama[it]==1) continue;
    if(event.isAccidental[it]==1) continue;
    if(!GFTrackCont.IsInsideTarget(it)) continue;
    if(event.isBeam[it]==1){
      target_accidental_id_container.push_back(it); //Accidental beam on the target
      continue;
    }
    if((event.pid[it]&2)==2 && event.charge[it]==-1){ //kaon

      Int_t repid = 1;
      if(GFTrackCont.TrackCheck(it, repid)){
	Int_t htofhitid_k; Double_t tracklen_k;
	Double_t tof; TVector3 pos; Double_t track2tgt_dist;
	Bool_t htofextrapolation_k = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								      event.HtofSeg, event.posHtof,
								      htofhitid_k, tof, tracklen_k,
								      pos, track2tgt_dist);
	Double_t mass2 = qnan;
	TVector3 mom = GFTrackCont.GetMom(it, 0, repid);
	if(htofextrapolation_k) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_k, event.tHtof[htofhitid_k]);
	target_k_id_container.push_back(it);
	target_k_mass2_container.push_back(mass2);
	target_k_mom_container.push_back(mom);
      }
    }

    if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton
      Int_t repid = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it];
	if(temp==flag) repid += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it, repid)) continue;

      Int_t htofhitid_p; Double_t tracklen_p;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								    event.HtofSeg, event.posHtof,
								    htofhitid_p, tof, tracklen_p,
								    pos, track2tgt_dist);
      Double_t mass2 = qnan;
      TVector3 mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_p) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_p, event.tHtof[htofhitid_p]);
      target_p_id_container.push_back(it);
      target_p_mass2_container.push_back(mass2);
      target_p_mom_container.push_back(mom);
    }
    else if((event.pid[it]&1)==1 && event.charge[it]==-1){ //pi-
      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(it, repid)) continue;

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      Double_t mass2 = qnan;
      TVector3 mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      //if(htofextrapolation_pi && mass2 > 0.25) continue;
      target_pim_id_container.push_back(it);
      target_pim_mass2_container.push_back(mass2);
      target_pim_mom_container.push_back(mom);
    }
    else if((event.pid[it]&4)!=4 && (event.pid[it]&1)==1 && event.charge[it]==1){ //pi+
      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(it, repid)) continue;

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								     event.HtofSeg, event.posHtof,
								     htofhitid_pi, tof, tracklen_pi,
								     pos, track2tgt_dist);
      Double_t mass2 = qnan;
      TVector3 mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation_pi) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      //if(htofextrapolation_pi && mass2 > 0.25) continue;
      target_pip_id_container.push_back(it);
      target_pip_mass2_container.push_back(mass2);
      target_pip_mom_container.push_back(mom);
    }
    else{ //p or pi
      Int_t repid = 0;
      if(!GFTrackCont.TrackCheck(it, repid)) continue;

      Int_t htofhitid_ppi; Double_t tracklen_ppi;
      Double_t tof; TVector3 pos; Double_t track2tgt_dist;
      Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(it, repid, tgtpos,
								  event.HtofSeg, event.posHtof,
								  htofhitid_ppi, tof, tracklen_ppi,
								  pos, track2tgt_dist);
      Double_t mass2 = qnan;
      TVector3 mom = GFTrackCont.GetMom(it, 0, repid);
      if(htofextrapolation) mass2 = Kinematics::MassSquare(mom.Mag(), tracklen_ppi, event.tHtof[htofhitid_ppi]);
      target_ppip_id_container.push_back(it);
      target_ppip_mass2_container.push_back(mass2);
      target_ppip_mom_container.push_back(mom);
    }
  }

  //for L, LL
  std::vector<Int_t> L_p_id_container, L_pi_id_container;
  std::vector<Int_t> L_p_repid_container, L_pi_repid_container;
  std::vector<TVector3> L_p_mom_container, L_pi_mom_container;
  std::vector<Double_t> L_mass_container;
  std::vector<Double_t> L_ppidist_container;
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_mom_container, L_vtx_container, L_targetvtx_container;

  //Xi/L candidates searching
#if DebugDisp
  std::cout<<"1. Xi/L candidate searching"<<std::endl;
#endif
  Int_t l_candidates = 0;
  Int_t xi_candidates = 0;
  std::vector<Int_t> xi_l_container;
  std::vector<Int_t> xi_p_container, xi_pi_container, xi_pi2_container;
  std::vector<Int_t> p_repid_container, pi_repid_container, pi2_repid_container;
  std::vector<TVector3> xi_mom_container, xi_decayvertex_container;
  std::vector<TVector3> l_mom_container, l_vert_container;
  std::vector<Double_t> xi_mass_container, lambda_mass_container;
  std::vector<TVector3> xi_p_mom_container, xi_pi_mom_container, xi_pi2_mom_container;
  std::vector<Double_t> ppi_closedist; std::vector<Double_t> lpi_closedist;
  std::vector<Double_t> xi_targetdist_container;
  std::vector<TVector3> xi_targetvtx_container, xi_targetmom_container;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&4)!=4) continue; //select p-like
      if(event.charge[it1]!=1) continue;

      Int_t repid_p = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it1];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it1, repid_p)) continue;

      Double_t p_par[5];
      p_par[0] = event.helix_cx[it1];
      p_par[1] = event.helix_cy[it1];
      p_par[2] = event.helix_z0[it1];
      p_par[3] = event.helix_r[it1];
      p_par[4] = event.helix_dz[it1];
      Int_t p_nh = event.helix_t[it1].size();
      Double_t p_theta_min = event.helix_t[it1][0] - vtx_scan_range/p_par[3];
      Double_t p_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_rangeInsideL/p_par[3], event.helix_t[it1][p_nh-1]);
      TVector3 p_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 p_end = TVector3(event.calpos_x[it1][p_nh-1], event.calpos_y[it1][p_nh-1], event.calpos_z[it1][p_nh-1]);

      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if((event.pid[it2]&1)!=1) continue; //select pi-like
	//if((event.pid[it2]&4)==4) continue; //veto p-like
	if(event.charge[it2]!=-1) continue;

	Int_t repid_pi = 0;
	if(!GFTrackCont.TrackCheck(it2, repid_pi)) continue;

	Double_t pi_par[5];
	pi_par[0] = event.helix_cx[it2];
	pi_par[1] = event.helix_cy[it2];
	pi_par[2] = event.helix_z0[it2];
	pi_par[3] = event.helix_r[it2];
	pi_par[4] = event.helix_dz[it2];

	Int_t pi_nh = event.helix_t[it2].size();
	Double_t pi_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/pi_par[3], event.helix_t[it2][pi_nh-1]);
	Double_t pi_theta_max = event.helix_t[it2][0] + vtx_scan_range/pi_par[3];
	TVector3 pi_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 pi_end = TVector3(event.calpos_x[it2][pi_nh-1], event.calpos_y[it2][pi_nh-1], event.calpos_z[it2][pi_nh-1]);

	Double_t ppi_dist = 10000.;
	TVector3 p_mom; TVector3 pi_mom; TVector3 lambda_mom;
	TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, p_mom, pi_mom, lambda_mom, ppi_dist);
	if(TMath::IsNaN(ppi_dist)) continue;
	lambda_mom = pi_mom + p_mom;

	TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
	TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
	TLorentzVector Llambda = Lp + Lpi;

	if(TMath::Abs(lambda_vert.x()) > 250. || TMath::Abs(lambda_vert.z()) > 250. || TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut
	if(ppi_dist > ppi_distcut) continue;
	Double_t pi_vertex_dist; Double_t p_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

	if(pi_vertex_dist > pi_vtx_distcut) continue;
 	if(p_vertex_dist > p_vtx_distcut) continue;
	if(TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;

	Double_t ltarget_dist;
	TVector3 ltarget_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
							       lambda_vert,
							       lambda_mom,
							       ltarget_dist);

	L_p_id_container.push_back(it1);
	L_pi_id_container.push_back(it2);
	L_mass_container.push_back(Llambda.M());
	L_mom_container.push_back(lambda_mom);
	L_p_mom_container.push_back(p_mom);
	L_pi_mom_container.push_back(pi_mom);
	L_ppidist_container.push_back(ppi_dist);
	L_vtx_container.push_back(lambda_vert);
	L_p_repid_container.push_back(repid_p);
	L_pi_repid_container.push_back(repid_pi);
	L_targetdist_container.push_back(ltarget_dist);
	L_targetvtx_container.push_back(ltarget_vtx);
	l_candidates++;

	for(int it3=0;it3<ntTpc;it3++){
	  if(it3==it2 || it3==it1) continue;
	  if(event.isK18[it3]==1) continue;
	  if(event.isKurama[it3]==1) continue;
	  if(event.isBeam[it3]==1) continue;
	  if(event.isAccidental[it3]==1) continue;
	  if((event.pid[it3]&1)!=1) continue; //select pi like
	  //if((event.pid[it3]&4)==4) continue; //veto p-like
	  if(event.charge[it3]!=-1) continue;

	  Int_t repid_pi2 = 0;
	  if(!GFTrackCont.TrackCheck(it3, repid_pi2)) continue;

	  Double_t pi2_par[5];
	  pi2_par[0] = event.helix_cx[it3];
	  pi2_par[1] = event.helix_cy[it3];
	  pi2_par[2] = event.helix_z0[it3];
	  pi2_par[3] = event.helix_r[it3];
	  pi2_par[4] = event.helix_dz[it3];

	  Int_t pi2_nh = event.helix_t[it3].size();
	  Double_t pi2_theta_min = TMath::Max(event.helix_t[it3][0] - vtx_scan_rangeInsidePi/pi2_par[3], event.helix_t[it3][pi2_nh-1]);
	  Double_t pi2_theta_max = event.helix_t[it3][0] + vtx_scan_range/pi2_par[3];

	  TVector3 pi2_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
	  TVector3 pi2_end = TVector3(event.calpos_x[it3][pi2_nh-1], event.calpos_y[it3][pi2_nh-1], event.calpos_z[it3][pi2_nh-1]);

	  TVector3 pi2_mom; Double_t lpi_dist;

	  TVector3 xi_vert = Kinematics::XiVertex(dMagneticField,
						  pi2_par,
						  pi2_theta_min,
						  pi2_theta_max,
						  lambda_vert,
						  lambda_mom,
						  pi2_mom,
						  lpi_dist);
	  if(TMath::IsNaN(lpi_dist)) continue;

	  TLorentzVector Lpi2(pi2_mom, TMath::Sqrt(pi2_mom.Mag()*pi2_mom.Mag() + PionMass*PionMass));
	  TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Sqrt(lambda_mom.Mag()*lambda_mom.Mag() + LambdaMass*LambdaMass));
	  TLorentzVector Lxi = Llambda_fixedmass + Lpi2;
	  TVector3 xi_mom = TVector3(Lxi.Px(), Lxi.Py(), Lxi.Pz());
	  Double_t pi2_vertex_dist;
	  if(lpi_dist > lpi_distcut) continue;
	  if(!Kinematics::HelixDirection(xi_vert, pi2_start, pi2_end, pi2_vertex_dist)) continue;
	  if(pi2_vertex_dist > pi2_vtx_distcut) continue;

	  if(TMath::Abs(Lxi.M() - XiMinusMass) > xi_masscut) continue; //Check reconstructed mass cut
	  Double_t xitarget_dist; TVector3 xi_mom_target;
	  TVector3 xi_vert_target = Kinematics::CalcCloseDistXi(tgtpos,
								dMagneticField,
								xi_vert,
								xi_mom,
								xi_mom_target,
								xitarget_dist);
	  if(xitarget_dist > xitarget_distcut) continue; //Closest distance between the xi and the target cut

	  xi_targetdist_container.push_back(xitarget_dist);
	  xi_targetvtx_container.push_back(xi_vert_target);
	  xi_targetmom_container.push_back(xi_mom_target);

	  ppi_closedist.push_back(ppi_dist);
	  lpi_closedist.push_back(lpi_dist);

	  xi_l_container.push_back(l_candidates - 1);
	  xi_p_container.push_back(it1);
	  xi_pi_container.push_back(it2);
	  xi_pi2_container.push_back(it3);

	  p_repid_container.push_back(repid_p);
	  pi_repid_container.push_back(repid_pi);
	  pi2_repid_container.push_back(repid_pi2);

	  xi_mass_container.push_back(Lxi.M());
	  lambda_mass_container.push_back(Llambda.M());

	  xi_mom_container.push_back(xi_mom);
	  l_mom_container.push_back(lambda_mom);
	  xi_p_mom_container.push_back(p_mom);
	  xi_pi_mom_container.push_back(pi_mom);
	  xi_pi2_mom_container.push_back(pi2_mom);

	  xi_decayvertex_container.push_back(xi_vert);
	  l_vert_container.push_back(lambda_vert);

	  xi_candidates++;
	} //it3
      } //it2
    } //it1
  }
#if DebugDisp
  std::cout<<"Xi/L candidates searching ends"<<std::endl;
#endif

  std::vector<Double_t> GFL_mass_container(l_candidates, qnan);
  std::vector<Double_t> GFL_ppidist_container(l_candidates, qnan);
  std::vector<Double_t> GFL_targetdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetvtx_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFL_targetcenterdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetcentervtx_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFL_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFL_vtx_container(l_candidates, TVector3(qnan, qnan, qnan));

  std::vector<Int_t> GFL_p_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_p_repid_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_repid_container(l_candidates, qnan);
  std::vector<TVector3> GFL_p_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFL_pi_mom_container(l_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFL_p_extrapolation_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_extrapolation_container(l_candidates, qnan);
  std::vector<Int_t> GFL_p_htofid_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_htofid_container(l_candidates, qnan);
  std::vector<Double_t> GFL_p_tracklen_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_tracklen_container(l_candidates, qnan);
  std::vector<Double_t> GFL_p_tof_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_tof_container(l_candidates, qnan);
  std::vector<Double_t> GFL_p_mass2_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_mass2_container(l_candidates, qnan);
  if(l_candidates>0){
#if DebugDisp
    std::cout<<"2. Single L, LL candidates searching starts"<<std::endl;
#endif
    //Reconstructed real Lambdas
    for(int idp=0;idp<l_candidates;idp++){
      if(L_targetdist_container[idp] > ltarget_distcut) continue;
      Int_t p_id = L_p_id_container[idp];
      Int_t pi_id = L_pi_id_container[idp];
      Int_t p_repid = L_p_repid_container[idp];
      Int_t pi_repid = L_pi_repid_container[idp];
      Double_t p_extrapolation; Double_t pi_extrapolation;
      TVector3 p_mom; TVector3 pi_mom;
      Double_t ppi_dist; TVector3 l_vertex;

      Bool_t vtxcut = (GFTrackCont.FindVertex(p_id, pi_id,
					      p_repid, pi_repid,
					      p_extrapolation, pi_extrapolation,
					      p_mom, pi_mom,
					      ppi_dist, l_vertex,
					      vtx_scan_range)
		       && ppi_dist < GFppi_distcut);
      if(!vtxcut) continue;

      TVector3 l_mom = p_mom + pi_mom;
      Double_t l_target_dist;
      TVector3 l_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, l_vertex, l_mom, l_target_dist);
      TVector3 l_fight = l_vertex - l_pos_tgt;
      Double_t l_tof = Kinematics::CalcTimeOfFlight(l_mom.Mag(), l_fight.Mag(), pdg::LambdaMass());
      if(l_target_dist > GFltarget_distcut) continue;


      Double_t l_targetcenter_dist;
      TVector3 l_pos_tgtcenter = Kinematics::LambdaTargetCenter(l_vertex, l_mom, l_targetcenter_dist);
      if(TMath::Abs(l_pos_tgtcenter.y()) > GFltarget_ycut) continue;

      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t dummy;
      Bool_t htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(p_id, p_repid, l_vertex,
								  event.HtofSeg, event.posHtof,
								  hitid_htof, tof_htof,
								  tracklen_htof, pos_htof, dummy);
      if(htofextrapolation){
	GFL_p_htofid_container[idp] = hitid_htof;
	GFL_p_tracklen_container[idp] = tracklen_htof;
	GFL_p_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_p_mass2_container[idp] = Kinematics::MassSquare(p_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	//if(GFL_p_mass2_container[idp] < 0.25) continue;
      }

      htofextrapolation = GFTrackCont.TPCHTOFTrackMatching(pi_id, pi_repid, l_vertex,
							   event.HtofSeg, event.posHtof,
							   hitid_htof, tof_htof,
							   tracklen_htof, pos_htof, dummy);
      if(htofextrapolation){
	GFL_pi_htofid_container[idp] = hitid_htof;
	GFL_pi_tracklen_container[idp] = tracklen_htof;
	GFL_pi_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_pi_mass2_container[idp] = Kinematics::MassSquare(pi_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	//if(GFL_pi_mass2_container[idp] > 0.25) continue;
      }

      TLorentzVector GFLp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
      TLorentzVector GFLpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
      TLorentzVector GFLlambda = GFLp + GFLpi;

      GFL_p_id_container[idp] = p_id;
      GFL_pi_id_container[idp] = pi_id;
      GFL_p_repid_container[idp] = p_repid;
      GFL_pi_repid_container[idp] = pi_repid;
      GFL_p_mom_container[idp] = p_mom;
      GFL_pi_mom_container[idp] = pi_mom;
      GFL_p_extrapolation_container[idp] = p_extrapolation;
      GFL_pi_extrapolation_container[idp] = pi_extrapolation;
      GFL_mass_container[idp] = GFLlambda.M();
      GFL_mom_container[idp] = l_mom;
      GFL_ppidist_container[idp] = ppi_dist;
      GFL_vtx_container[idp] = l_vertex;
      GFL_targetdist_container[idp] = l_target_dist;
      GFL_targetvtx_container[idp] = l_pos_tgt;
      GFL_targetcenterdist_container[idp] = l_targetcenter_dist;
      GFL_targetcentervtx_container[idp] = l_pos_tgtcenter;
    }

    //order : L1, L2
    std::vector<std::vector<Int_t>> GFLLid_container;
    std::vector<std::vector<TVector3>> GFLLvtx_container;
    std::vector<std::vector<TVector3>> GFLLmom_container;
    std::vector<std::vector<Double_t>> GFLLmass_container;
    std::vector<std::vector<Double_t>> GFLLppidist_container;
    std::vector<std::vector<Double_t>> GFLLtargetdist_container;
    std::vector<std::vector<TVector3>> GFLLtargetvtx_container;
    std::vector<std::vector<Double_t>> GFLLtargetcenterdist_container;
    std::vector<std::vector<TVector3>> GFLLtargetcentervtx_container;
    //order : p, pi(L1), p, pi(L2)
    std::vector<std::vector<Int_t>> GFLLdecays_trackid_container;
    std::vector<std::vector<Int_t>> GFLLdecays_repid_container;
    std::vector<std::vector<Int_t>> GFLLdecays_htofid_container;
    std::vector<std::vector<Double_t>> GFLLdecays_tof_container;
    std::vector<std::vector<TVector3>> GFLLdecays_mom_container;
    std::vector<std::vector<Double_t>> GFLLdecays_tracklen_container;
    std::vector<std::vector<Double_t>> GFLLdecays_mass2_container;

    Int_t LL_best = -1; Int_t LLcount = 0;
    Double_t prev_LLmassdiff = 9999.;
    for(int idp1=0;idp1<l_candidates;idp1++){
      if(TMath::IsNaN(GFL_mass_container[idp1])) continue;

      //For LL searching
      for(int idp2=idp1+1;idp2<l_candidates;idp2++){
	if(TMath::IsNaN(GFL_mass_container[idp2])) continue;
	if(GFL_p_id_container[idp1]==GFL_p_id_container[idp2] ||
	   GFL_p_id_container[idp1]==GFL_pi_id_container[idp2]) continue;
	if(GFL_pi_id_container[idp1]==GFL_p_id_container[idp2] ||
	   GFL_pi_id_container[idp1]==GFL_pi_id_container[idp2]) continue;

	std::vector<Int_t> GFlambda_containerid(2);
	std::vector<TVector3> GFlambda_vert(2);
	std::vector<TVector3> GFlambda_mom(2);
	std::vector<Double_t> GFlambda_mass(2);
	std::vector<Double_t> GFppi_dist(2);
	std::vector<Double_t> GFlambdatgt_dist(2);
	std::vector<TVector3> GFlambdatgt_vtx(2);
	std::vector<Double_t> GFlambdatgtcenter_dist(2);
	std::vector<TVector3> GFlambdatgtcenter_vtx(2);

	std::vector<Int_t> GFtrackid_decays(4);
	std::vector<Int_t> GFrepid_decays(4);
	std::vector<Int_t> GFhtofid_decays(4);
	std::vector<Double_t> GFtof_decays(4);
	std::vector<TVector3> GFmom_decays(4);
	std::vector<Double_t> GFtracklen_decays(4);
	std::vector<Double_t> GFmass2_decays(4);

	GFlambda_containerid[0] = idp1;
	GFlambda_containerid[1] = idp2;
	GFlambda_vert[0] = GFL_vtx_container[idp1];
	GFlambda_vert[1] = GFL_vtx_container[idp2];
	GFlambda_mom[0] = GFL_mom_container[idp1];
	GFlambda_mom[1] = GFL_mom_container[idp2];
	GFlambda_mass[0] = GFL_mass_container[idp1];
	GFlambda_mass[1] = GFL_mass_container[idp2];
	GFppi_dist[0] = GFL_ppidist_container[idp1];
	GFppi_dist[1] = GFL_ppidist_container[idp2];
	GFlambdatgt_dist[0] = GFL_targetdist_container[idp1];
	GFlambdatgt_dist[1] = GFL_targetdist_container[idp2];
	GFlambdatgt_vtx[0] = GFL_targetvtx_container[idp1];
	GFlambdatgt_vtx[1] = GFL_targetvtx_container[idp2];
	GFlambdatgtcenter_dist[0] = GFL_targetcenterdist_container[idp1];
	GFlambdatgtcenter_dist[1] = GFL_targetcenterdist_container[idp2];
	GFlambdatgtcenter_vtx[0] = GFL_targetcentervtx_container[idp1];
	GFlambdatgtcenter_vtx[1] = GFL_targetcentervtx_container[idp2];

	GFtrackid_decays[0] = GFL_p_id_container[idp1];
	GFtrackid_decays[1] = GFL_pi_id_container[idp1];
	GFtrackid_decays[2] = GFL_p_id_container[idp2];
	GFtrackid_decays[3] = GFL_pi_id_container[idp2];
	GFrepid_decays[0] = GFL_p_repid_container[idp1];
	GFrepid_decays[1] = GFL_pi_repid_container[idp1];
	GFrepid_decays[2] = GFL_p_repid_container[idp2];
	GFrepid_decays[3] = GFL_pi_repid_container[idp2];
	GFhtofid_decays[0] = GFL_p_htofid_container[idp1];
	GFhtofid_decays[1] = GFL_pi_htofid_container[idp1];
	GFhtofid_decays[2] = GFL_p_htofid_container[idp2];
	GFhtofid_decays[3] = GFL_pi_htofid_container[idp2];
	GFtof_decays[0] = GFL_p_tof_container[idp1];
	GFtof_decays[1] = GFL_pi_tof_container[idp1];
	GFtof_decays[2] = GFL_p_tof_container[idp2];
	GFtof_decays[3] = GFL_pi_tof_container[idp2];
	GFmom_decays[0] = GFL_p_mom_container[idp1];
	GFmom_decays[1] = GFL_pi_mom_container[idp1];
	GFmom_decays[2] = GFL_p_mom_container[idp2];
	GFmom_decays[3] = GFL_pi_mom_container[idp2];
	GFtracklen_decays[0] = GFL_p_tracklen_container[idp1];
	GFtracklen_decays[1] = GFL_pi_tracklen_container[idp1];
	GFtracklen_decays[2] = GFL_p_tracklen_container[idp2];
	GFtracklen_decays[3] = GFL_pi_tracklen_container[idp2];
	GFmass2_decays[0] = GFL_p_mass2_container[idp1];
	GFmass2_decays[1] = GFL_pi_mass2_container[idp1];
	GFmass2_decays[2] = GFL_p_mass2_container[idp2];
	GFmass2_decays[3] = GFL_pi_mass2_container[idp2];

	//order : L1, L2
	GFLLid_container.push_back(GFlambda_containerid);
	GFLLvtx_container.push_back(GFlambda_vert);
	GFLLmom_container.push_back(GFlambda_mom);
	GFLLmass_container.push_back(GFlambda_mass);
	GFLLppidist_container.push_back(GFppi_dist);
	GFLLtargetdist_container.push_back(GFlambdatgt_dist);
	GFLLtargetvtx_container.push_back(GFlambdatgt_vtx);
	GFLLtargetcenterdist_container.push_back(GFlambdatgtcenter_dist);
	GFLLtargetcentervtx_container.push_back(GFlambdatgtcenter_vtx);

	//order : p, pi(L1), p, pi(L2)
	GFLLdecays_trackid_container.push_back(GFtrackid_decays);
	GFLLdecays_repid_container.push_back(GFrepid_decays);
	GFLLdecays_htofid_container.push_back(GFhtofid_decays);
	GFLLdecays_tof_container.push_back(GFtof_decays);
	GFLLdecays_mom_container.push_back(GFmom_decays);
	GFLLdecays_tracklen_container.push_back(GFtracklen_decays);
	GFLLdecays_mass2_container.push_back(GFmass2_decays);

	event.llflag = true;
	Double_t diff = TMath::Hypot(GFlambda_mass[0] - LambdaMass, GFlambda_mass[1] - LambdaMass);
	//Double_t diff = TMath::Hypot(L_mass_container[idp1] - LambdaMass, L_mass_container[idp2] - LambdaMass);

	if(prev_LLmassdiff > diff){
	  prev_LLmassdiff = diff;
	  LL_best = LLcount;
	}
	LLcount++;
      } //L combi1
    } //L combi2

    //for the LL event
    if(event.llflag){
#if DebugDisp
    std::cout<<"3. Calculate and save the best LL pair"<<std::endl;
#endif
      Int_t L1_id = GFLLid_container[LL_best][0];
      Int_t L2_id = GFLLid_container[LL_best][1];

      event.lmass1 = L_mass_container.at(L1_id);
      event.ldecayvtx_x1 = L_vtx_container.at(L1_id).x();
      event.ldecayvtx_y1 = L_vtx_container.at(L1_id).y();
      event.ldecayvtx_z1 = L_vtx_container.at(L1_id).z();
      event.lmom1 = L_mom_container.at(L1_id).Mag();
      event.lmom_x1 = L_mom_container.at(L1_id).x();
      event.lmom_y1 = L_mom_container.at(L1_id).y();
      event.lmom_z1 = L_mom_container.at(L1_id).z();
      event.ppi_dist1 = L_ppidist_container.at(L1_id);
      event.ltarget_dist1 = L_targetdist_container.at(L1_id);
      event.ltargetvtx_x1 = L_targetvtx_container.at(L1_id).x();
      event.ltargetvtx_y1 = L_targetvtx_container.at(L1_id).y();
      event.ltargetvtx_z1 = L_targetvtx_container.at(L1_id).z();

      event.lmass2 = L_mass_container.at(L2_id);
      event.ldecayvtx_x2 = L_vtx_container.at(L2_id).x();
      event.ldecayvtx_y2 = L_vtx_container.at(L2_id).y();
      event.ldecayvtx_z2 = L_vtx_container.at(L2_id).z();
      event.lmom2 = L_mom_container.at(L2_id).Mag();
      event.lmom_x2 = L_mom_container.at(L2_id).x();
      event.lmom_y2 = L_mom_container.at(L2_id).y();
      event.lmom_z2 = L_mom_container.at(L2_id).z();
      event.ppi_dist2 = L_ppidist_container.at(L2_id);
      event.ltarget_dist2 = L_targetdist_container.at(L2_id);
      event.ltargetvtx_x2 = L_targetvtx_container.at(L2_id).x();
      event.ltargetvtx_y2 = L_targetvtx_container.at(L2_id).y();
      event.ltargetvtx_z2 = L_targetvtx_container.at(L2_id).z();

      event.GFlmass1 = GFLLmass_container[LL_best][0];
      event.GFldecayvtx_x1 = GFLLvtx_container[LL_best][0].x();
      event.GFldecayvtx_y1 = GFLLvtx_container[LL_best][0].y();
      event.GFldecayvtx_z1 = GFLLvtx_container[LL_best][0].z();
      event.GFlmom1 = GFLLmom_container[LL_best][0].Mag();
      event.GFlmom_x1 = GFLLmom_container[LL_best][0].x();
      event.GFlmom_y1 = GFLLmom_container[LL_best][0].y();
      event.GFlmom_z1 = GFLLmom_container[LL_best][0].z();
      event.GFppi_dist1 = GFLLppidist_container[LL_best][0];
      event.GFltarget_dist1 = GFLLtargetdist_container[LL_best][0];
      event.GFltargetvtx_x1 = GFLLtargetvtx_container[LL_best][0].x();
      event.GFltargetvtx_y1 = GFLLtargetvtx_container[LL_best][0].y();
      event.GFltargetvtx_z1 = GFLLtargetvtx_container[LL_best][0].z();
      event.GFltargetcenter_dist1 = GFLLtargetcenterdist_container[LL_best][0];
      event.GFltargetcenter_x1 = GFLLtargetcentervtx_container[LL_best][0].x();
      event.GFltargetcenter_y1 = GFLLtargetcentervtx_container[LL_best][0].y();
      event.GFltargetcenter_z1 = GFLLtargetcentervtx_container[LL_best][0].z();

      event.GFlmass2 = GFLLmass_container[LL_best][1];
      event.GFldecayvtx_x2 = GFLLvtx_container[LL_best][1].x();
      event.GFldecayvtx_y2 = GFLLvtx_container[LL_best][1].y();
      event.GFldecayvtx_z2 = GFLLvtx_container[LL_best][1].z();
      event.GFlmom2 = GFLLmom_container[LL_best][1].Mag();
      event.GFlmom_x2 = GFLLmom_container[LL_best][1].x();
      event.GFlmom_y2 = GFLLmom_container[LL_best][1].y();
      event.GFlmom_z2 = GFLLmom_container[LL_best][1].z();
      event.GFppi_dist2 = GFLLppidist_container[LL_best][1];
      event.GFltarget_dist2 = GFLLtargetdist_container[LL_best][1];
      event.GFltargetvtx_x2 = GFLLtargetvtx_container[LL_best][1].x();
      event.GFltargetvtx_y2 = GFLLtargetvtx_container[LL_best][1].y();
      event.GFltargetvtx_z2 = GFLLtargetvtx_container[LL_best][1].z();
      event.GFltargetcenter_dist2 = GFLLtargetcenterdist_container[LL_best][1];
      event.GFltargetcenter_x2 = GFLLtargetcentervtx_container[LL_best][1].x();
      event.GFltargetcenter_y2 = GFLLtargetcentervtx_container[LL_best][1].y();
      event.GFltargetcenter_z2 = GFLLtargetcentervtx_container[LL_best][1].z();

	Double_t ll_dist;
	TVector3 ll_vtx1, ll_vtx2;
	TVector3 ll_vtx
	  = Kinematics::LambdaLambdaVertex(L_vtx_container.at(L1_id),
					   L_mom_container.at(L1_id),
					   L_vtx_container.at(L2_id),
					   L_mom_container.at(L2_id),
					   ll_vtx1, ll_vtx2, ll_dist);
	event.llvtx_x = ll_vtx.x();
	event.llvtx_y = ll_vtx.y();
	event.llvtx_z = ll_vtx.z();
	event.lldist = ll_dist;

	Double_t GFll_dist;
	TVector3 GFll_vtx1, GFll_vtx2;
	TVector3 GFll_vtx
	  = Kinematics::LambdaLambdaVertex(GFLLvtx_container[LL_best][0],
					   GFLLmom_container[LL_best][0],
					   GFLLvtx_container[LL_best][1],
					   GFLLmom_container[LL_best][1],
					   GFll_vtx1, GFll_vtx2, GFll_dist);
	event.GFllvtx_x = GFll_vtx.x();
	event.GFllvtx_y = GFll_vtx.y();
	event.GFllvtx_z = GFll_vtx.z();
	event.GFlldist = GFll_dist;

	const Int_t ntrack = 4;
	Double_t x0[ntrack] = {event.xtgtK18[0], event.xtgtTPCKurama[0],
			       event.GFltargetcenter_x1, event.GFltargetcenter_x2};
	Double_t y0[ntrack] = {event.ytgtK18[0], event.ytgtTPCKurama[0],
			       event.GFltargetcenter_y1, event.GFltargetcenter_y2};
	Double_t u0[ntrack] = {event.utgtK18[0], event.utgtTPCKurama[0],
			       event.GFlmom_x1/event.GFlmom_z1,
			       event.GFlmom_x2/event.GFlmom_z2};
	Double_t v0[ntrack] = {event.vtgtK18[0], event.vtgtTPCKurama[0],
			       event.GFlmom_y1/event.GFlmom_z1,
			       event.GFlmom_y2/event.GFlmom_z2};
	TVector3 kk_ll_vertex = Kinematics::MultitrackVertex(ntrack, x0, y0, u0, v0);

	event.GFprodvtx_x = kk_ll_vertex.x();
	event.GFprodvtx_y = kk_ll_vertex.y();
	event.GFprodvtx_z = kk_ll_vertex.z();
	Double_t prodvtx_closedist1 = qnan;
	TVector3 prodvtx_closest1 = Kinematics::CalcCloseDistLambda(kk_ll_vertex,
								    GFLLvtx_container[LL_best][0],
								    GFLLmom_container[LL_best][0],
								    prodvtx_closedist1);

	Double_t prodvtx_closedist2 = qnan;
	TVector3 prodvtx_closest2 = Kinematics::CalcCloseDistLambda(kk_ll_vertex,
								    GFLLvtx_container[LL_best][1],
								    GFLLmom_container[LL_best][1],
								    prodvtx_closedist2);

#if DebugDisp
	std::cout<<"K-K+LL vertex "<<kk_ll_vertex
		 <<" Lambda1's closest point to the vertex"<<prodvtx_closest1
		 <<" Lambda2's closest point to the vertex"<<prodvtx_closest2
		 <<std::endl;
#endif
	event.GFlprodvtx_x1 = prodvtx_closest1.x();
	event.GFlprodvtx_y1 = prodvtx_closest1.y();
	event.GFlprodvtx_z1 = prodvtx_closest1.z();
	event.GFlprodvtx_dist1 = prodvtx_closedist1;

	event.GFlprodvtx_x2 = prodvtx_closest2.x();
	event.GFlprodvtx_y2 = prodvtx_closest2.y();
	event.GFlprodvtx_z2 = prodvtx_closest2.z();
	event.GFlprodvtx_dist2 = prodvtx_closedist2;

	TVector3 lambda_tracklen = kk_ll_vertex - GFLLvtx_container[LL_best][0];
	event.GFltracklen1 = lambda_tracklen.Mag();
	event.GFltof1 = Kinematics::CalcTimeOfFlight(event.GFlmom1, lambda_tracklen.Mag(), pdg::LambdaMass());

	lambda_tracklen = kk_ll_vertex - GFLLvtx_container[LL_best][1];
	event.GFltracklen2 = lambda_tracklen.Mag();
	event.GFltof2 = Kinematics::CalcTimeOfFlight(event.GFlmom2, lambda_tracklen.Mag(), pdg::LambdaMass());

      //p, pi(L1), p, pi(L2) tracks
      event.GFntdecays = 4;
      event.GFnhtrack.resize(event.GFntdecays);
      event.GFchisqr.resize(event.GFntdecays);
      event.GFcharge.resize(event.GFntdecays);
      event.GFtof.resize(event.GFntdecays);
      event.GFtracklen.resize(event.GFntdecays);
      event.GFpval.resize(event.GFntdecays);
      event.GFpdgcode.resize(event.GFntdecays);

      event.GFdecays_htofid.resize(event.GFntdecays);
      event.GFdecays_tracklen.resize(event.GFntdecays);
      event.GFdecays_tof.resize(event.GFntdecays);
      event.GFdecays_mass2.resize(event.GFntdecays);
      event.GFdecays_mom.resize(event.GFntdecays);
      event.GFdecays_mom_x.resize(event.GFntdecays);
      event.GFdecays_mom_y.resize(event.GFntdecays);
      event.GFdecays_mom_z.resize(event.GFntdecays);
      event.GFmomloss.resize(event.GFntdecays);
      event.GFeloss.resize(event.GFntdecays);

      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][0]);
      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][1]);
      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][2]);
      event.decays_id.push_back(GFLLdecays_trackid_container[LL_best][3]);

      event.decays_mom.push_back(L_p_mom_container[L1_id].Mag());
      event.decays_mom.push_back(L_pi_mom_container[L1_id].Mag());
      event.decays_mom.push_back(L_p_mom_container[L2_id].Mag());
      event.decays_mom.push_back(L_pi_mom_container[L2_id].Mag());
      event.decays_mom_x.push_back(L_p_mom_container[L1_id].x());
      event.decays_mom_x.push_back(L_pi_mom_container[L1_id].x());
      event.decays_mom_x.push_back(L_p_mom_container[L2_id].x());
      event.decays_mom_x.push_back(L_pi_mom_container[L2_id].x());
      event.decays_mom_y.push_back(L_p_mom_container[L1_id].y());
      event.decays_mom_y.push_back(L_pi_mom_container[L1_id].y());
      event.decays_mom_y.push_back(L_p_mom_container[L2_id].y());
      event.decays_mom_y.push_back(L_pi_mom_container[L2_id].y());
      event.decays_mom_z.push_back(L_p_mom_container[L1_id].z());
      event.decays_mom_z.push_back(L_pi_mom_container[L1_id].z());
      event.decays_mom_z.push_back(L_p_mom_container[L2_id].z());
      event.decays_mom_z.push_back(L_pi_mom_container[L2_id].z());

      for(int j=0;j<event.GFntdecays;j++){
	Int_t igf = GFLLdecays_trackid_container[LL_best][j];
	Int_t repid = GFLLdecays_repid_container[LL_best][j];
	event.GFnhtrack[j] = GFTrackCont.GetNHits(igf);
	event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
	event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
	event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
	event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
	event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
	event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);

	event.GFdecays_htofid[j] = GFLLdecays_htofid_container[LL_best][j];
	event.GFdecays_tracklen[j] = GFLLdecays_tracklen_container[LL_best][j];
	event.GFdecays_tof[j] = GFLLdecays_tof_container[LL_best][j];
	event.GFdecays_mass2[j] = GFLLdecays_mass2_container[LL_best][j];
	event.GFdecays_mom[j] = GFLLdecays_mom_container[LL_best][j].Mag();
	event.GFdecays_mom_x[j] = GFLLdecays_mom_container[LL_best][j].x();
	event.GFdecays_mom_y[j] = GFLLdecays_mom_container[LL_best][j].y();
	event.GFdecays_mom_z[j] = GFLLdecays_mom_container[LL_best][j].z();
	event.GFmomloss[j] = qnan;
	event.GFeloss[j] = qnan;
      } //GFntdecays

      HF2( 502, event.lmass1, event.lmass2);
      HF2( 503, event.GFlmass1, event.GFlmass2);
      HF1( 504, event.lmass1);
      HF1( 504, event.lmass2);
      HF1( 505, event.GFlmass1);
      HF1( 505, event.GFlmass2);

      //Remaining tracks
#if DebugDisp
  std::cout<<"Save remaining tracks"<<std::endl;
#endif
      event.pip_multi = target_pip_id_container.size();
      event.p_multi = target_p_id_container.size();
      event.pim_multi = target_pim_id_container.size();
      event.ppip_multi = target_ppip_id_container.size();
      event.accdient_multi = target_accidental_id_container.size();
      for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
	Int_t id_p = target_p_id_container[itp];
	if(id_p==event.decays_id[0] || id_p==event.decays_id[2]){
	  event.p_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_p);
	event.residual_mass2.push_back(target_p_mass2_container[itp]);
	event.residual_mom.push_back(target_p_mom_container[itp].Mag());
	event.residual_mom_x.push_back(target_p_mom_container[itp].x());
	event.residual_mom_y.push_back(target_p_mom_container[itp].y());
	event.residual_mom_z.push_back(target_p_mom_container[itp].z());
	event.residual_charge.push_back(1);
	HF1( 1110, target_p_mom_container[itp].Mag());
      }
      for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
	Int_t id_pip = target_pip_id_container[itpip];
	if(id_pip==event.decays_id[0] || id_pip==event.decays_id[2]){
	  event.pip_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_pip);
	event.residual_mass2.push_back(target_pip_mass2_container[itpip]);
	event.residual_mom.push_back(target_pip_mom_container[itpip].Mag());
	event.residual_mom_x.push_back(target_pip_mom_container[itpip].x());
	event.residual_mom_y.push_back(target_pip_mom_container[itpip].y());
	event.residual_mom_z.push_back(target_pip_mom_container[itpip].z());
	event.residual_charge.push_back(1);
	HF1( 1111, target_pip_mom_container[itpip].Mag());
      }
      for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
	Int_t id_pim = target_pim_id_container[itpim];
	if(id_pim==event.decays_id[1] || id_pim==event.decays_id[3]){
	  event.pim_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_pim);
	event.residual_mass2.push_back(target_pim_mass2_container[itpim]);
	event.residual_mom.push_back(target_pim_mom_container[itpim].Mag());
	event.residual_mom_x.push_back(target_pim_mom_container[itpim].x());
	event.residual_mom_y.push_back(target_pim_mom_container[itpim].y());
	event.residual_mom_z.push_back(target_pim_mom_container[itpim].z());
	event.residual_charge.push_back(-1);
	HF1( 1112, target_pim_mom_container[itpim].Mag());
      }
      for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
	Int_t id_ppip = target_ppip_id_container[itppip];
	if(id_ppip==event.decays_id[0] || id_ppip==event.decays_id[2]){
	  event.ppip_multi -= 1;
	  continue;
	}
	event.residual_id.push_back(id_ppip);
	event.residual_mass2.push_back(target_ppip_mass2_container[itppip]);
	event.residual_mom.push_back(target_ppip_mom_container[itppip].Mag());
	event.residual_mom_x.push_back(target_ppip_mom_container[itppip].x());
	event.residual_mom_y.push_back(target_ppip_mom_container[itppip].y());
	event.residual_mom_z.push_back(target_ppip_mom_container[itppip].z());
	event.residual_charge.push_back(event.charge[id_ppip]);
      }
      event.residual_multi = event.pip_multi + event.p_multi + event.pim_multi + event.ppip_multi + event.accdient_multi;

      HF1( 1010, event.pim_multi);
      HF1( 1011, event.p_multi);
      HF1( 1012, event.pip_multi);

      HF1( 110, -event.BE[0]);
      if(event.isgoodTPC[0] == 1) HF1( 112, -event.BETPC[0]);
      return true;
    } //LL flag
  } //l_candidates

  std::vector<Int_t> GFxi_p_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_id_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_p_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_pi2_rep_container(xi_candidates, -1);
  std::vector<Int_t> GFxi_id_container(xi_candidates, -1);
  std::vector<TVector3> GFxi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFxi_decayvertex_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFxi_mass_container(xi_candidates, qnan);
  std::vector<TVector3> GFxi_pos_targetcenter_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFxi_mom_targetcenter_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFxi_dist_targetcenter_container(xi_candidates, qnan);
  std::vector<TVector3> GFl_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFl_vert_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFl_mass_container(xi_candidates, qnan);
  std::vector<Double_t> GFl_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFl_tof_container(xi_candidates, qnan);
  std::vector<TVector3> GFp_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpi2_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFp_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi2_mass2_container(xi_candidates, qnan);
  std::vector<Double_t> GFp_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi2_tof_container(xi_candidates, qnan);
  std::vector<Double_t> GFp_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi_tracklen_container(xi_candidates, qnan);
  std::vector<Double_t> GFpi2_tracklen_container(xi_candidates, qnan);
  std::vector<Int_t> GFp_htofid_container(xi_candidates, qnan);
  std::vector<Int_t> GFpi_htofid_container(xi_candidates, qnan);
  std::vector<Int_t> GFpi2_htofid_container(xi_candidates, qnan);
  std::vector<Double_t> GFppi_closedist_container(xi_candidates, qnan);
  std::vector<Double_t> GFlpi_closedist_container(xi_candidates, qnan);

  std::vector<Double_t> KFlchisqr_container(xi_candidates, qnan); 
  std::vector<Double_t> KFlpval_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_chisqr_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_pval_container(xi_candidates, qnan);
  std::vector<Double_t> KFxi_mass_container(xi_candidates, qnan);
  std::vector<Double_t> KFlpi_closedist_container(xi_candidates, qnan);
  std::vector<std::vector<Double_t>>KFlpull_container(xi_candidates, std::vector<Double_t>(6,qnan));
  std::vector<std::vector<Double_t>>KFxipull_container(xi_candidates, std::vector<Double_t>(6,qnan));
#if DebugDisp
  if(xi_candidates>0) std::cout<<"4. Xi candidates searching starts"<<std::endl;
#endif
  Int_t best_xi = -1; Double_t prev_Lmassdiff = 9999.;
  for(Int_t candi=0;candi<xi_candidates;candi++){
    Int_t trackid_p = xi_p_container[candi];
    Int_t trackid_pi = xi_pi_container[candi];
    Int_t trackid_pi2 = xi_pi2_container[candi];
    Int_t repid_p = p_repid_container[candi];
    Int_t repid_pi = pi_repid_container[candi];
    Int_t repid_pi2 = pi2_repid_container[candi];
    Int_t l_id = xi_l_container[candi];
    if(TMath::IsNaN(GFL_ppidist_container[l_id])) continue; //Genfit's fitting was succeeded.

    double GFppi_dist = GFL_ppidist_container[l_id];
    TVector3 GFlambda_vert = GFL_vtx_container[l_id];
    Double_t GFextrapolation_decays[3];
    GFextrapolation_decays[0] = GFL_p_extrapolation_container[l_id];
    GFextrapolation_decays[1] = GFL_pi_extrapolation_container[l_id];
    TVector3 GFmom_decays[3];
    GFmom_decays[0] = GFL_p_mom_container[l_id];
    GFmom_decays[1] = GFL_pi_mom_container[l_id];
    TLorentzVector GFLp(GFmom_decays[0], TMath::Hypot(GFmom_decays[0].Mag(), ProtonMass));
    TLorentzVector GFLpi(GFmom_decays[1], TMath::Hypot(GFmom_decays[1].Mag(), PionMass));
    TLorentzVector GFLlambda = GFLp + GFLpi;
    TVector3 GFlambda_mom = GFmom_decays[0] + GFmom_decays[1];
#if DoKF
    Double_t KFchisqrxi=-1;
    Double_t KFpvalxi=-1;
    Double_t KFchisqrl=-1;
    Double_t KFpvall=-1;
    TPCLocalTrackHelix *track_p = TPCAna.GetTrackTPCHelix(trackid_p);
    auto Vp = track_p->GetCovarianceMatrix();
    TPCLocalTrackHelix *track_pi = TPCAna.GetTrackTPCHelix(trackid_pi);
    auto Vpi1 = track_pi->GetCovarianceMatrix();
    double Diag_ppi1[6]={
      Vp(0,0),Vp(1,1),Vp(2,2),Vpi1(0,0),Vpi1(1,1),Vpi1(2,2)
    };
    auto Offdiag_ppi1 = MathTools::MergeOffdiagonals(Vp,Vpi1); 
    //Y and Z coordinates should be swapped in KF. 
    TVector3 HTVP(GFLp.X(),GFLp.Z(),GFLp.Y());
    TVector3 HTVPi1(GFLpi.X(),GFLpi.Z(),GFLpi.Y());
    TVector3 HTVLd = HTVP+HTVPi1;
    TLorentzVector HLVP(HTVP,hypot(ProtonMass,HTVP.Mag()));
    TLorentzVector HLVPi1(HTVPi1,hypot(PionMass,HTVPi1.Mag()));
    TLorentzVector HLVLd(HTVLd,hypot(LambdaMass,HTVLd.Mag()));

    cout<<"KFLd"<<endl;
    FourVectorFitter KFLd(HLVP,HLVPi1,HLVLd);
    KFLd.SetInvMass(LambdaMass);
    KFLd.SetMaximumStep(5);
    KFLd.SetVariance(Diag_ppi1);
    KFLd.AddOffdiagonals(Offdiag_ppi1);
    KFchisqrl = KFLd.DoKinematicFit();
    KFpvall = KFLd.GetPValue();
    auto HcontLd = KFLd.GetFittedLV();
    auto PullLd = KFLd.GetPull();
    auto KFHLVP = HcontLd.at(0);
    auto KFHLVPi1 = HcontLd.at(1);
    auto KFHLVLd = HcontLd.at(2);
    auto KFlambda_mom = TVector3(KFHLVLd.X(),KFHLVLd.Z(),KFHLVLd.Y());
    TLorentzVector KFLlambda_fixed(KFlambda_mom,hypot(KFlambda_mom.Mag(),LambdaMass));  
    auto VLd = KFLd.GetUnmeasuredCovariance();
    TVector3 KFxi_vert; Double_t KFlpi_dist = qnan; Double_t KFlambda_tracklen;
    double l_res_x,l_res_y,l_phi;
    MathTools::DecomposeResolution(VLd,l_res_x,l_res_y,l_phi); 
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, KFlambda_mom,
    KFlambda_tracklen, GFextrapolation_decays[2], GFmom_decays[2], KFlpi_dist, KFxi_vert, vtx_scan_range,l_res_x,l_res_y,l_phi) 
      || KFlpi_dist > GFlpi_distcut) continue;
    TLorentzVector KFLpi2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + PionMass*PionMass));
    TLorentzVector KFLxi = KFLlambda_fixed + KFLpi2;
    TVector3 KFxi_mom = KFlambda_mom + GFmom_decays[2];

    Double_t KFlambda_tof = Kinematics::CalcTimeOfFlight(KFlambda_mom.Mag(), KFlambda_tracklen, pdg::LambdaMass());
    GFlambda_mom = KFlambda_mom;
    KFlpull_container[candi] = PullLd;
#endif
    TLorentzVector GFLlambda_fixed(GFlambda_mom, TMath::Sqrt(GFlambda_mom.Mag()*GFlambda_mom.Mag() + LambdaMass*LambdaMass));

    TVector3 GFxi_vert; Double_t GFlpi_dist = qnan; Double_t GFlambda_tracklen;
    cout<<Form("Resolution: %f %f %f",l_res_x,l_res_y,l_phi)<<endl;
    if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, GFlambda_mom, GFlambda_tracklen,
     GFextrapolation_decays[2], GFmom_decays[2], GFlpi_dist, GFxi_vert, vtx_scan_range)
      || GFlpi_dist > GFlpi_distcut) continue;

    TLorentzVector GFLpi2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + PionMass*PionMass));
    TLorentzVector GFLxi = GFLlambda_fixed + GFLpi2;
    TVector3 GFxi_mom = GFlambda_mom + GFmom_decays[2];

    Double_t GFlambda_tof = Kinematics::CalcTimeOfFlight(GFlambda_mom.Mag(), GFlambda_tracklen, pdg::LambdaMass());

    Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof;
    Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(trackid_p, repid_p,
                  event.HtofSeg, event.posHtof,
                  hitid_htof, tof_htof,
                  tracklen_htof, pos_htof);
    if(htofextrapolation_p){
      GFp_htofid_container[candi] = hitid_htof;
      //GFp_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
      GFp_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof;
      GFp_tracklen_container[candi] = tracklen_htof - GFextrapolation_decays[0];
      GFp_mass2_container[candi] = Kinematics::MassSquare(GFmom_decays[0].Mag(), GFp_tracklen_container[candi], GFp_tof_container[candi]);
      //if(GFL_p_mass2_container[candi] < 0.25) continue;
    }

    Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(trackid_pi, repid_pi,
                    event.HtofSeg, event.posHtof,
                    hitid_htof, tof_htof,
                    tracklen_htof, pos_htof);
    if(htofextrapolation_pi){
      GFpi_htofid_container[candi] = hitid_htof;
      //GFpi_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
      GFpi_tof_container[candi] = event.tHtof[hitid_htof] - GFlambda_tof;
      GFpi_tracklen_container[candi] = tracklen_htof - GFextrapolation_decays[1];
      GFpi_mass2_container[candi] = Kinematics::MassSquare(GFmom_decays[1].Mag(), GFpi_tracklen_container[candi], GFpi_tof_container[candi]);
      //if(GFL_pi_mass2_container[candi] > 0.25) continue;
    }

    Bool_t htofextrapolation_pi2 = GFTrackCont.TPCHTOFTrackMatching(trackid_pi2, repid_pi2,
                    event.HtofSeg, event.posHtof,
                    hitid_htof, tof_htof,
                    tracklen_htof, pos_htof);
    if(htofextrapolation_pi2){
      GFpi2_htofid_container[candi] = hitid_htof;
      //GFpi2_tof_container[candi] = event.tHtof[hitid_htof] - GFxi_tof;
      GFpi2_tof_container[candi] = event.tHtof[hitid_htof];
      GFpi2_tracklen_container[candi] = tracklen_htof - GFextrapolation_decays[2];
      GFpi2_mass2_container[candi] = Kinematics::MassSquare(GFmom_decays[2].Mag(), GFpi2_tracklen_container[candi], GFpi2_tof_container[candi]);
      //if(GFL_pi2_mass2_container[candi] > 0.25) continue;
    }
#if DoKF

    auto* track_pi2 = TPCAna.GetTrackTPCHelix(trackid_pi2);
    auto VPi2 = track_pi2->GetCovarianceMatrix();
    double Diag_lpi2[6] =
    {VLd(0,0),VLd(1,1),VLd(2,2),VPi2(0,0),VPi2(1,1),VPi2(2,2)};
    auto Offdiag_lpi2 = MathTools::MergeOffdiagonals(VLd,VPi2);

    TVector3 HTVPi2(GFLpi2.X(),GFLpi2.Z(),GFLpi2.Y());
    TLorentzVector HLVPi2(HTVPi2,hypot(PionMass,HTVPi2.Mag())); 
    auto HLVXi = KFHLVLd + HLVPi2;
    FourVectorFitter KFXi(KFHLVLd,HLVPi2,HLVXi);
    KFXi.SetInvMass(XiMinusMass);
    cout<<"KFXi"<<endl;
    KFXi.SetMaximumStep(5);
    KFXi.SetVariance(Diag_lpi2);
    KFXi.AddOffdiagonals(Offdiag_lpi2);
    KFchisqrxi = KFXi.DoKinematicFit();
    KFpvalxi = KFXi.GetPValue(); 
    auto HcontXi = KFXi.GetFittedLV();
    auto PullXi = KFXi.GetPull();
    auto KFKFHLVLd = HcontXi.at(0);
    auto KFHLVPi2 = HcontXi.at(1);
    auto KFHLVXi = HcontXi.at(2);

    auto KFLVP = TLorentzVector(KFHLVP.X(),KFHLVP.Z(),KFHLVP.Y(),KFHLVP.E());
    auto KFLVPi1 = TLorentzVector(KFHLVPi1.X(),KFHLVPi1.Z(),KFHLVPi1.Y(),KFHLVPi1.E());
    auto KFLVPi2 = TLorentzVector(KFHLVPi2.X(),KFHLVPi2.Z(),KFHLVPi2.Y(),KFHLVPi2.E());
    auto KFLVLd = TLorentzVector(KFKFHLVLd.X(),KFKFHLVLd.Z(),KFKFHLVLd.Y(),KFKFHLVLd.E());
    auto KFLVXi = TLorentzVector(KFHLVXi.X(),KFHLVXi.Z(),KFHLVXi.Y(),KFHLVXi.E());

    auto KFTVP = KFLVP.Vect();
    auto KFTVPi1 = KFLVPi1.Vect();
    auto KFTVPi2 = KFLVPi2.Vect();
    auto KFTVLd = KFLVLd.Vect();
    auto KFTVXi = KFLVXi.Vect();
    GFxi_mom = KFTVXi;
    GFlambda_mom = KFTVLd;
    KFxipull_container[candi] = PullXi;
#endif
    GFTrackCont.AddReconstructedTrack(XiMinusPdgCode, KFxi_vert, GFxi_mom);
    Int_t xi_trackid = GFTrackCont.GetNTrack() - 1;
    GFTrackCont.FitTrack(xi_trackid);

    TVector3 xipos_tgtcenter; TVector3 ximom_tgtcenter;
    Double_t xitracklen_tgtcenter; Double_t xitof_tgtcenter;
    Bool_t xi_extrapolation = GFTrackCont.ExtrapolateToTargetCenter(xi_trackid,
								    xipos_tgtcenter,
								    ximom_tgtcenter,
								    xitracklen_tgtcenter,
								    xitof_tgtcenter);
    if(!xi_extrapolation) continue;
    if(TMath::Abs(xipos_tgtcenter.y()) > GFxitarget_ycut) continue;

    event.xiflag = true;
    TVector3 dist = tgtpos - xipos_tgtcenter;
    GFxi_dist_targetcenter_container[candi] = dist.Mag();
    GFxi_pos_targetcenter_container[candi] = xipos_tgtcenter;
    GFxi_mom_targetcenter_container[candi] = ximom_tgtcenter;
    GFxi_id_container[candi] = xi_trackid;

    GFxi_p_id_container[candi] = trackid_p;
    GFxi_pi_id_container[candi] = trackid_pi;
    GFxi_pi2_id_container[candi] = trackid_pi2;
    GFxi_p_rep_container[candi] = repid_p;
    GFxi_pi_rep_container[candi] = repid_pi;
    GFxi_pi2_rep_container[candi] = repid_pi2;
    GFxi_mom_container[candi] = GFxi_mom;
    GFxi_decayvertex_container[candi] = GFxi_vert;
    GFxi_mass_container[candi] = GFLxi.M();

    GFl_mom_container[candi] = GFlambda_mom;
    GFl_vert_container[candi] = GFlambda_vert;
    GFl_mass_container[candi] = GFLlambda.M();
    GFl_tracklen_container[candi] = GFlambda_tracklen;
    GFl_tof_container[candi] = GFlambda_tof;
    GFp_mom_container[candi] = GFmom_decays[0];
    GFpi_mom_container[candi] = GFmom_decays[1];
    GFpi2_mom_container[candi] = GFmom_decays[2];
    GFppi_closedist_container[candi] = GFppi_dist;
    GFlpi_closedist_container[candi] = GFlpi_dist;
#if DoKF
    KFlchisqr_container[candi] = KFchisqrl;
    KFlpval_container[candi] = KFpvall;
    KFxi_chisqr_container[candi] = KFchisqrxi;
    KFxi_pval_container[candi] = KFpvalxi;
    KFxi_mass_container[candi] = KFLxi.M();
    KFlpi_closedist_container[candi] = KFlpi_dist;
#endif
    Double_t diff = TMath::Abs(GFLlambda.M() - LambdaMass);
    //Double_t diff = TMath::Abs(lambda_mass_container[candi] - LambdaMass);
    if(prev_Lmassdiff > diff){
      prev_Lmassdiff = diff;
      best_xi = candi;
    }
  } //candi

  if(event.xiflag){
#if DebugDisp
    std::cout<<"5. Calculate and save the best Xi combination"<<std::endl;
#endif
    event.ximass = xi_mass_container[best_xi];
    event.xidecayvtx_x = xi_decayvertex_container[best_xi].x();
    event.xidecayvtx_y = xi_decayvertex_container[best_xi].y();
    event.xidecayvtx_z = xi_decayvertex_container[best_xi].z();
    event.ximom = xi_mom_container[best_xi].Mag();
    event.ximom_x = xi_mom_container[best_xi].x();
    event.ximom_y = xi_mom_container[best_xi].y();
    event.ximom_z = xi_mom_container[best_xi].z();
    event.lpi_dist = lpi_closedist[best_xi];
    event.xitarget_dist = xi_targetdist_container[best_xi];
    event.xitargetvtx_x = xi_targetvtx_container[best_xi].x();
    event.xitargetvtx_y = xi_targetvtx_container[best_xi].y();
    event.xitargetvtx_z = xi_targetvtx_container[best_xi].z();
    event.xitargetmom = xi_targetmom_container[best_xi].Mag();
    event.xitargetmom_x = xi_targetmom_container[best_xi].x();
    event.xitargetmom_y = xi_targetmom_container[best_xi].y();
    event.xitargetmom_z = xi_targetmom_container[best_xi].z();
    event.lmass = lambda_mass_container[best_xi];
    event.ldecayvtx_x = l_vert_container[best_xi].x();
    event.ldecayvtx_y = l_vert_container[best_xi].y();
    event.ldecayvtx_z = l_vert_container[best_xi].z();
    event.lmom = l_mom_container[best_xi].Mag();
    event.lmom_x = l_mom_container[best_xi].x();
    event.lmom_y = l_mom_container[best_xi].y();
    event.lmom_z = l_mom_container[best_xi].z();
    event.ppi_dist = ppi_closedist[best_xi];


    event.decays_mom.push_back(xi_p_mom_container[best_xi].Mag());
    event.decays_mom.push_back(xi_pi_mom_container[best_xi].Mag());
    event.decays_mom.push_back(xi_pi2_mom_container[best_xi].Mag());
    event.decays_mom_x.push_back(xi_p_mom_container[best_xi].x());
    event.decays_mom_x.push_back(xi_pi_mom_container[best_xi].x());
    event.decays_mom_x.push_back(xi_pi2_mom_container[best_xi].x());
    event.decays_mom_y.push_back(xi_p_mom_container[best_xi].y());
    event.decays_mom_y.push_back(xi_pi_mom_container[best_xi].y());
    event.decays_mom_y.push_back(xi_pi2_mom_container[best_xi].y());
    event.decays_mom_z.push_back(xi_p_mom_container[best_xi].z());
    event.decays_mom_z.push_back(xi_pi_mom_container[best_xi].z());
    event.decays_mom_z.push_back(xi_pi2_mom_container[best_xi].z());
    event.decays_id.push_back(xi_p_container[best_xi]);
    event.decays_id.push_back(xi_pi_container[best_xi]);
    event.decays_id.push_back(xi_pi2_container[best_xi]);

    event.GFximass = GFxi_mass_container[best_xi];
    event.GFxidecayvtx_x = GFxi_decayvertex_container[best_xi].x();
    event.GFxidecayvtx_y = GFxi_decayvertex_container[best_xi].y();
    event.GFxidecayvtx_z = GFxi_decayvertex_container[best_xi].z();
    event.GFximom = GFxi_mom_container[best_xi].Mag();
    event.GFximom_x = GFxi_mom_container[best_xi].x();
    event.GFximom_y = GFxi_mom_container[best_xi].y();
    event.GFximom_z = GFxi_mom_container[best_xi].z();
    event.GFlpi_dist = GFlpi_closedist_container[best_xi];
    event.GFlmass = GFl_mass_container[best_xi];
    event.GFldecayvtx_x = GFl_vert_container[best_xi].x();
    event.GFldecayvtx_y = GFl_vert_container[best_xi].y();
    event.GFldecayvtx_z = GFl_vert_container[best_xi].z();
    event.GFlmom = GFl_mom_container[best_xi].Mag();
    event.GFlmom_x = GFl_mom_container[best_xi].x();
    event.GFlmom_y = GFl_mom_container[best_xi].y();
    event.GFlmom_z = GFl_mom_container[best_xi].z();
    event.GFltracklen = GFl_tracklen_container[best_xi];
    event.GFltof = GFl_tof_container[best_xi];
    event.GFppi_dist = GFppi_closedist_container[best_xi];
    
    event.KFlchisqr = KFlchisqr_container[best_xi];
    event.KFlpval = KFlpval_container[best_xi];
    event.KFxichisqr = KFxi_chisqr_container[best_xi];
    event.KFxipval = KFxi_pval_container[best_xi];
    event.KFximass = KFxi_mass_container[best_xi];
    event.KFlpi_dist = KFlpi_closedist_container[best_xi];
    event.KFlpull = KFlpull_container[best_xi];
    event.KFxipull = KFxipull_container[best_xi];

    for(int i=0;i<event.KFlpull.size();++i){
      int hn = 10000 + 10 + i;
      auto pull = event.KFlpull.at(i);
      HF1(hn,pull);
    }
    for(int i=0;i<event.KFxipull.size();++i){
      int hn = 20000 + 10 + i;
      auto pull = event.KFxipull.at(i);
      HF1(hn,pull);
    }

    HF1(10000,event.KFlpval);
    HF1(10001,event.KFlchisqr);
    HF1(20002,event.KFximass);
    HF1(10003,event.KFlpi_dist);





    event.GFntdecays = 3;
    event.GFnhtrack.resize(event.GFntdecays);
    event.GFchisqr.resize(event.GFntdecays);
    event.GFcharge.resize(event.GFntdecays);
    event.GFtof.resize(event.GFntdecays);
    event.GFtracklen.resize(event.GFntdecays);
    event.GFpval.resize(event.GFntdecays);
    event.GFpdgcode.resize(event.GFntdecays);

    event.GFdecays_htofid.resize(event.GFntdecays);
    event.GFdecays_tracklen.resize(event.GFntdecays);
    event.GFdecays_tof.resize(event.GFntdecays);
    event.GFdecays_mass2.resize(event.GFntdecays);
    event.GFdecays_mom.resize(event.GFntdecays);
    event.GFdecays_mom_x.resize(event.GFntdecays);
    event.GFdecays_mom_y.resize(event.GFntdecays);
    event.GFdecays_mom_z.resize(event.GFntdecays);
    event.GFmomloss.resize(event.GFntdecays);
    event.GFeloss.resize(event.GFntdecays);

    for(Int_t j=0; j<event.GFntdecays; ++j){
      Int_t igf = GFxi_p_id_container[best_xi];
      if(j==1) igf = GFxi_pi_id_container[best_xi];
      if(j==2) igf = GFxi_pi2_id_container[best_xi];

      Int_t repid = GFxi_p_rep_container[best_xi];
      if(j==1) repid = GFxi_pi_rep_container[best_xi];
      if(j==2) repid = GFxi_pi2_rep_container[best_xi];

      Int_t GFhtofid_decays = GFp_htofid_container[best_xi];
      if(j==1) GFhtofid_decays = GFpi_htofid_container[best_xi];
      if(j==2) GFhtofid_decays = GFpi2_htofid_container[best_xi];

      Double_t GFtof_decays = GFp_tof_container[best_xi];
      if(j==1) GFtof_decays = GFpi_tof_container[best_xi];
      if(j==2) GFtof_decays = GFpi2_tof_container[best_xi];

      Double_t GFtracklen_decays = GFp_tracklen_container[best_xi];
      if(j==1) GFtracklen_decays = GFpi_tracklen_container[best_xi];
      if(j==2) GFtracklen_decays = GFpi2_tracklen_container[best_xi];

      TVector3 GFmom_decays = GFp_mom_container[best_xi];
      if(j==1) GFmom_decays = GFpi_mom_container[best_xi];
      if(j==2) GFmom_decays = GFpi2_mom_container[best_xi];

      Double_t GFmass2_decays = GFp_mass2_container[best_xi];
      if(j==1) GFmass2_decays = GFpi_mass2_container[best_xi];
      if(j==2) GFmass2_decays = GFpi2_mass2_container[best_xi];

      event.GFdecays_htofid[j] = GFhtofid_decays;
      event.GFdecays_tracklen[j] = GFtracklen_decays;
      event.GFdecays_tof[j] = GFtof_decays;
      event.GFdecays_mass2[j] = GFmass2_decays;
      event.GFdecays_mom[j] = GFmom_decays.Mag();
      event.GFdecays_mom_x[j] = GFmom_decays.x();
      event.GFdecays_mom_y[j] = GFmom_decays.y();
      event.GFdecays_mom_z[j] = GFmom_decays.z();
      event.GFmomloss[j] = GFmom_decays.Mag() - GFTrackCont.GetMom(igf, 0, repid).Mag();
      event.GFeloss[j] = TMath::Hypot(GFmom_decays.Mag(), pdgmass[j]) - TMath::Hypot(GFTrackCont.GetMom(igf, 0, repid).Mag(), pdgmass[j]);

      event.GFnhtrack[j] = GFTrackCont.GetNHits(igf);
      event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
      event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
      event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
      event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
      event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
      event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
    } //j : decays

    if(event.isgoodTPC[0] == 1){
      TVector3 xipos_tgtcenter = GFxi_pos_targetcenter_container[best_xi];
      TVector3 ximom_tgtcenter = GFxi_mom_targetcenter_container[best_xi];
      event.GFxitargetcenter_x = xipos_tgtcenter.x();
      event.GFxitargetcenter_y = xipos_tgtcenter.y();
      event.GFxitargetcenter_z = xipos_tgtcenter.z();
      event.GFxitargetcentermom = ximom_tgtcenter.Mag();
      event.GFxitargetcentermom_x = ximom_tgtcenter.x();
      event.GFxitargetcentermom_y = ximom_tgtcenter.y();
      event.GFxitargetcentermom_z = ximom_tgtcenter.z();
      event.GFxitargetcenter_dist = GFxi_dist_targetcenter_container[best_xi];

      const Int_t ntrack = 3;
      Double_t x0[ntrack] = {event.xtgtTPCKurama[0], event.xtgtK18[0],
			     xipos_tgtcenter.x()};
      Double_t y0[ntrack] = {event.ytgtTPCKurama[0], event.ytgtK18[0],
			     xipos_tgtcenter.y()};
      Double_t u0[ntrack] = {event.utgtTPCKurama[0], event.utgtK18[0],
			     ximom_tgtcenter.x()/ximom_tgtcenter.z()};
      Double_t v0[ntrack] = {event.vtgtTPCKurama[0], event.vtgtK18[0],
			     ximom_tgtcenter.y()/ximom_tgtcenter.z()};
      TVector3 kk_xi_vertex = Kinematics::MultitrackVertex(ntrack, x0, y0, u0, v0);
      event.GFprodvtx_x = kk_xi_vertex.x();
      event.GFprodvtx_y = kk_xi_vertex.y();
      event.GFprodvtx_z = kk_xi_vertex.z();

      TVector3 xipos_prodvtx; TVector3 ximom_prodvtx;
      Double_t xitracklen_prodvtx; Double_t xitof_prodvtx;
      Bool_t xi_extrapolation_prodvtx = GFTrackCont.XiDecayToProdVertex(GFxi_id_container[best_xi],
									kk_xi_vertex, xipos_prodvtx,
									ximom_prodvtx,
									xitracklen_prodvtx,
									xitof_prodvtx);

      if(xi_extrapolation_prodvtx){
	TVector3 dist = xipos_prodvtx - kk_xi_vertex;
	event.GFxiprodvtx_dist = dist.Mag();
	event.GFxiprodvtx_x = xipos_prodvtx.x();
	event.GFxiprodvtx_y = xipos_prodvtx.y();
	event.GFxiprodvtx_z = xipos_prodvtx.z();
	event.GFxiprodmom = ximom_prodvtx.Mag();
	event.GFxiprodmom_x = ximom_prodvtx.x();
	event.GFxiprodmom_y = ximom_prodvtx.y();
	event.GFxiprodmom_z = ximom_prodvtx.z();
	event.GFxitracklen = xitracklen_prodvtx;
	event.GFxitof = xitof_prodvtx;
	event.GFximomloss = ximom_prodvtx.Mag() - event.GFximom;

	TLorentzVector LvXi_fixedmass(ximom_prodvtx, TMath::Hypot(ximom_prodvtx.Mag(), XiMinusMass));
	TLorentzVector LvRcXi_fixedmass = LvRcTPC - LvXi_fixedmass;
	event.GFxiexcitation = LvRcXi_fixedmass.M() - m11B - me;

	HF1( 200, LvRcXi_fixedmass.M() - m11B - me);

#if DebugDisp
	std::cout<<"K-K+Xi vertex "<<kk_xi_vertex
		 <<" Xi extrapolate to the vertex "<<xipos_prodvtx<<std::endl;
#endif
      }

      TVector3 xipos_kk; TVector3 ximom_kk;
      Double_t xitracklen_kk; Double_t xitof_kk;
      Bool_t xi_extrapolation = GFTrackCont.XiDecayToProdVertex(GFxi_id_container[best_xi],
								kkvtxTPC,
								xipos_kk,
								ximom_kk,
								xitracklen_kk,
								xitof_kk);
      if(xi_extrapolation){
	TVector3 dist = xipos_kk - kkvtxTPC;
	event.GFxikkvtx_x = xipos_kk.x();
	event.GFxikkvtx_y = xipos_kk.y();
	event.GFxikkvtx_z = xipos_kk.z();
	event.GFxikkmom = ximom_kk.Mag();
	event.GFxikkmom_x = ximom_kk.x();
	event.GFxikkmom_y = ximom_kk.y();
	event.GFxikkmom_z = ximom_kk.z();
	event.GFxikkvtx_dist = dist.Mag();
      }

      TVector3 xipos_tgt; TVector3 ximom_tgt;
      Double_t xitracklen_tgt; Double_t xitof_tgt;
      xi_extrapolation = GFTrackCont.XiDecayToProdVertex(GFxi_id_container[best_xi],
							 tgtpos,
							 xipos_tgt,
							 ximom_tgt,
							 xitracklen_tgt,
							 xitof_tgt);
      if(xi_extrapolation){
	TVector3 dist = tgtpos - xipos_tgt;
	event.GFxitargetvtx_x = xipos_tgt.x();
	event.GFxitargetvtx_y = xipos_tgt.y();
	event.GFxitargetvtx_z = xipos_tgt.z();
	event.GFxitargetmom = ximom_tgt.Mag();
	event.GFxitargetmom_x = ximom_tgt.x();
	event.GFxitargetmom_y = ximom_tgt.y();
	event.GFxitargetmom_z = ximom_tgt.z();
	event.GFxitarget_dist = dist.Mag();
      }
    } //if(event.isgoodTPC[0] == 1)

    HF1( 10, event.MissMass[0] );
    HF1( 11, event.lmass );
    HF1( 12, event.ximass );
    HF1( 13, event.decays_mom[0] );
    HF1( 14, event.decays_mom[1] );
    HF1( 15, event.decays_mom[2] );
    HF1( 16, event.lmom );
    HF1( 17, event.ximom );
    HF1( 18, event.ppi_dist);
    HF1( 19, event.lpi_dist);

    HF1( 21, event.GFlmass);
    HF1( 22, event.GFximass);
    HF1( 23, event.GFdecays_mom[0]);
    HF1( 24, event.GFdecays_mom[1]);
    HF1( 25, event.GFdecays_mom[2]);
    HF1( 26, event.GFlmom);
    HF1( 27, event.GFximom);
    HF1( 28, event.GFppi_dist);
    HF1( 29, event.GFlpi_dist);
    HF1( 30, event.GFdecays_mass2[0]);
    HF1( 31, event.GFdecays_mass2[1]);
    HF1( 32, event.GFdecays_mass2[2]);
    HF1( 33, event.GFdecays_tof[0] - Kinematics::CalcTimeOfFlight(event.GFdecays_mom[0],
								  event.GFdecays_tracklen[0],
								  pdg::ProtonMass()));
    HF1( 34, event.GFdecays_tof[1] - Kinematics::CalcTimeOfFlight(event.GFdecays_mom[1],
								  event.GFdecays_tracklen[1],
								  pdg::PionMass()));
    HF1( 35, event.GFdecays_tof[2] - Kinematics::CalcTimeOfFlight(event.GFdecays_mom[2],
								  event.GFdecays_tracklen[2],
								  pdg::PionMass()));
    HF2( 36, event.GFdecays_mass2[0], event.GFdecays_mom[0]);
    HF2( 37, -event.GFdecays_mass2[1], event.GFdecays_mom[1]);
    HF2( 38, -event.GFdecays_mass2[2], event.GFdecays_mom[2]);

  //Remaining tracks
#if DebugDisp
  std::cout<<"Save remaining tracks"<<std::endl;
#endif
    event.pip_multi = target_pip_id_container.size();
    event.p_multi = target_p_id_container.size();
    event.pim_multi = target_pim_id_container.size();
    event.ppip_multi = target_ppip_id_container.size();
    event.accdient_multi = target_accidental_id_container.size();
    for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
      Int_t id_p = target_p_id_container[itp];
      if(id_p==event.decays_id[0]){
	event.p_multi -= 1;
	continue;
      }
      event.residual_id.push_back(target_p_id_container[itp]);
      event.residual_mass2.push_back(target_p_mass2_container[itp]);
      event.residual_mom.push_back(target_p_mom_container[itp].Mag());
      event.residual_mom_x.push_back(target_p_mom_container[itp].x());
      event.residual_mom_y.push_back(target_p_mom_container[itp].y());
      event.residual_mom_z.push_back(target_p_mom_container[itp].z());
      event.residual_charge.push_back(1);
      HF1( 1120, target_p_mom_container[itp].Mag());
    }
    for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
      Int_t id_pip = target_pip_id_container[itpip];
      if(id_pip==event.decays_id[0]){
	event.pip_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_pip);
      event.residual_mass2.push_back(target_pip_mass2_container[itpip]);
      event.residual_mom.push_back(target_pip_mom_container[itpip].Mag());
      event.residual_mom_x.push_back(target_pip_mom_container[itpip].x());
      event.residual_mom_y.push_back(target_pip_mom_container[itpip].y());
      event.residual_mom_z.push_back(target_pip_mom_container[itpip].z());
      event.residual_charge.push_back(1);
      HF1( 1121, target_pip_mom_container[itpip].Mag());
    }
    for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
      Int_t id_pim = target_pim_id_container[itpim];
      if(id_pim==event.decays_id[1] || id_pim==event.decays_id[2]){
	event.pim_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_pim);
      event.residual_mass2.push_back(target_pim_mass2_container[itpim]);
      event.residual_mom.push_back(target_pim_mom_container[itpim].Mag());
      event.residual_mom_x.push_back(target_pim_mom_container[itpim].x());
      event.residual_mom_y.push_back(target_pim_mom_container[itpim].y());
      event.residual_mom_z.push_back(target_pim_mom_container[itpim].z());
      event.residual_charge.push_back(-1);
      HF1( 1122, target_pim_mom_container[itpim].Mag());
    }
    for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
      Int_t id_ppip = target_ppip_id_container[itppip];
      if(id_ppip==event.decays_id[0]){
	event.ppip_multi -= 1;
	continue;
      }
      event.residual_id.push_back(id_ppip);
      event.residual_mass2.push_back(target_ppip_mass2_container[itppip]);
      event.residual_mom.push_back(target_ppip_mom_container[itppip].Mag());
      event.residual_mom_x.push_back(target_ppip_mom_container[itppip].x());
      event.residual_mom_y.push_back(target_ppip_mom_container[itppip].y());
      event.residual_mom_z.push_back(target_ppip_mom_container[itppip].z());
      event.residual_charge.push_back(event.charge[id_ppip]);
    }
    event.residual_multi = event.pip_multi + event.p_multi + event.pim_multi + event.ppip_multi + event.accdient_multi;

    HF1( 1020, event.pim_multi);
    HF1( 1021, event.p_multi);
    HF1( 1022, event.pip_multi);

    //find Xi-p event
    if(event.residual_multi==1 && event.p_multi==1){
      HF1( 500, -event.BE[0]);
    }

    HF1( 120, -event.BE[0]);
    if(event.isgoodTPC[0] == 1) HF1( 122, -event.BETPC[0]);
    return true;
  } //if(event.xiflag)

#if DebugDisp
  std::cout<<"6. Single L searching starts (No Xi-, LL)"<<std::endl;
#endif
  Int_t GFtrackid_decays[2] = {-1, -1}; Int_t GFrepid_decays[2];

  for(Int_t id=0; id<l_candidates; ++id){
    if(TMath::IsNaN(GFL_mass_container[id])) continue;
    if(TMath::IsNaN(GFL_ppidist_container[id])) continue; //Genfit's fitting was succeeded.
    if(L_targetdist_container[id] > ltarget_distcut) continue; //Select Lambda from the traget
    if(GFL_targetdist_container[id] > GFltarget_distcut) continue;

    event.lflag = true;
    event.lmass = L_mass_container[id];
    event.ldecayvtx_x = L_vtx_container[id].x();
    event.ldecayvtx_y = L_vtx_container[id].y();
    event.ldecayvtx_z = L_vtx_container[id].z();
    event.lmom = L_mom_container[id].Mag();
    event.lmom_x = L_mom_container[id].x();
    event.lmom_y = L_mom_container[id].y();
    event.lmom_z = L_mom_container[id].z();
    event.ppi_dist = L_ppidist_container[id];
    event.ltarget_dist = L_targetdist_container[id];
    event.ltargetvtx_x = L_targetvtx_container[id].x();
    event.ltargetvtx_y = L_targetvtx_container[id].y();
    event.ltargetvtx_z = L_targetvtx_container[id].z();

    event.decays_mom.push_back(L_p_mom_container[id].Mag());
    event.decays_mom.push_back(L_pi_mom_container[id].Mag());
    event.decays_mom_x.push_back(L_p_mom_container[id].x());
    event.decays_mom_x.push_back(L_pi_mom_container[id].x());
    event.decays_mom_y.push_back(L_p_mom_container[id].y());
    event.decays_mom_y.push_back(L_pi_mom_container[id].y());
    event.decays_mom_z.push_back(L_p_mom_container[id].z());
    event.decays_mom_z.push_back(L_pi_mom_container[id].z());
    event.decays_id.push_back(L_p_id_container[id]);
    event.decays_id.push_back(L_pi_id_container[id]);

    event.GFlmass = GFL_mass_container[id];
    event.GFldecayvtx_x = GFL_vtx_container[id].x();
    event.GFldecayvtx_y = GFL_vtx_container[id].y();
    event.GFldecayvtx_z = GFL_vtx_container[id].z();
    event.GFlmom = GFL_mom_container[id].Mag();
    event.GFlmom_x = GFL_mom_container[id].x();
    event.GFlmom_y = GFL_mom_container[id].y();
    event.GFlmom_z = GFL_mom_container[id].z();
    event.GFppi_dist = GFL_ppidist_container[id];
    event.GFltarget_dist = GFL_targetdist_container[id];
    event.GFltargetvtx_x = GFL_targetvtx_container[id].x();
    event.GFltargetvtx_y = GFL_targetvtx_container[id].y();
    event.GFltargetvtx_z = GFL_targetvtx_container[id].z();
    event.GFltargetcenter_dist = GFL_targetcenterdist_container[id];
    event.GFltargetcenter_x = GFL_targetcentervtx_container[id].x();
    event.GFltargetcenter_y = GFL_targetcentervtx_container[id].y();
    event.GFltargetcenter_z = GFL_targetcentervtx_container[id].z();

    GFtrackid_decays[0] = GFL_p_id_container[id];
    GFtrackid_decays[1] = GFL_pi_id_container[id];
    GFrepid_decays[0] = GFL_p_repid_container[id];
    GFrepid_decays[1] = GFL_pi_repid_container[id];

    if(event.isgoodTPC[0] == 1){
      const Int_t ntrack = 3;
      Double_t x0[ntrack] = {event.xtgtK18[0], event.xtgtTPCKurama[0], event.GFltargetcenter_x};
      Double_t y0[ntrack] = {event.ytgtK18[0], event.ytgtTPCKurama[0], event.GFltargetcenter_y};
      Double_t u0[ntrack] = {event.utgtK18[0], event.utgtTPCKurama[0], event.GFlmom_x/event.GFlmom_z};
      Double_t v0[ntrack] = {event.vtgtK18[0], event.vtgtTPCKurama[0], event.GFlmom_y/event.GFlmom_z};
      TVector3 kk_l_vertex = Kinematics::MultitrackVertex(ntrack, x0, y0, u0, v0);

      event.GFprodvtx_x = kk_l_vertex.x();
      event.GFprodvtx_y = kk_l_vertex.y();
      event.GFprodvtx_z = kk_l_vertex.z();
      Double_t prodvtx_closedist = qnan;
      TVector3 prodvtx_closest = Kinematics::CalcCloseDistLambda(kk_l_vertex,
								 GFL_vtx_container[id],
								 GFL_mom_container[id],
								 prodvtx_closedist);

      event.GFlprodvtx_x = prodvtx_closest.x();
      event.GFlprodvtx_y = prodvtx_closest.y();
      event.GFlprodvtx_z = prodvtx_closest.z();
      event.GFlprodvtx_dist = prodvtx_closedist;

      TVector3 lambda_tracklen = kk_l_vertex - GFL_vtx_container[id];
      event.GFltracklen = lambda_tracklen.Mag();
      event.GFltof = Kinematics::CalcTimeOfFlight(event.GFlmom, lambda_tracklen.Mag(), pdg::LambdaMass());
    }

    event.GFdecays_htofid.push_back(GFL_p_htofid_container[id]);
    event.GFdecays_htofid.push_back(GFL_pi_htofid_container[id]);
    event.GFdecays_tracklen.push_back(GFL_p_tracklen_container[id]);
    event.GFdecays_tracklen.push_back(GFL_pi_tracklen_container[id]);
    event.GFdecays_tof.push_back(GFL_p_tof_container[id]);
    event.GFdecays_tof.push_back(GFL_pi_tof_container[id]);
    event.GFdecays_mass2.push_back(GFL_p_mass2_container[id]);
    event.GFdecays_mass2.push_back(GFL_pi_mass2_container[id]);

    event.GFdecays_mom.push_back(GFL_p_mom_container[id].Mag());
    event.GFdecays_mom.push_back(GFL_pi_mom_container[id].Mag());
    event.GFdecays_mom_x.push_back(GFL_p_mom_container[id].x());
    event.GFdecays_mom_x.push_back(GFL_pi_mom_container[id].x());
    event.GFdecays_mom_y.push_back(GFL_p_mom_container[id].y());
    event.GFdecays_mom_y.push_back(GFL_pi_mom_container[id].y());
    event.GFdecays_mom_z.push_back(GFL_p_mom_container[id].z());
    event.GFdecays_mom_z.push_back(GFL_pi_mom_container[id].z());
    event.GFmomloss.push_back(qnan);
    event.GFmomloss.push_back(qnan);
    event.GFeloss.push_back(qnan);
    event.GFeloss.push_back(qnan);

    //order : p, pi
    event.GFntdecays = 2;
    event.GFnhtrack.resize(event.GFntdecays);
    event.GFchisqr.resize(event.GFntdecays);
    event.GFcharge.resize(event.GFntdecays);
    event.GFtof.resize(event.GFntdecays);
    event.GFtracklen.resize(event.GFntdecays);
    event.GFpval.resize(event.GFntdecays);
    event.GFpdgcode.resize(event.GFntdecays);

    for(int j=0;j<event.GFntdecays;j++){
      Int_t igf = GFtrackid_decays[j];
      Int_t repid = GFrepid_decays[j];
      event.GFnhtrack[j] = GFTrackCont.GetNHits(igf);
      event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
      event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
      event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
      event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
      event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
      event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
    } //GFntdecays
    HF1( 130, -event.BE[0]);
    if(event.isgoodTPC[0] == 1) HF1( 132, -event.BETPC[0]);

    //Lphi searching
    if(target_k_id_container.size()==1){
#if DebugDisp
    std::cout<<"optional: L phi searching "<<std::endl;
#endif
      Int_t id_kp = -1;
      for(Int_t k=0;k<ntTpc;k++){
	if(event.isKurama[k]==1) id_kp = k;
      }

      Int_t id_km = target_k_id_container[0];
      Double_t km_mass2 = target_k_mass2_container[0];
      //Double_t km_mom_helix = target_k_mom_container[0];

      Double_t km_par[5];
      km_par[0] = event.helix_cx[id_km];
      km_par[1] = event.helix_cy[id_km];
      km_par[2] = event.helix_z0[id_km];
      km_par[3] = event.helix_r[id_km];
      km_par[4] = event.helix_dz[id_km];

      Int_t km_nh = event.helix_t[id_km].size();
      Double_t km_theta_min = TMath::Max(event.helix_t[id_km][0] - vtx_scan_rangeInsideL/km_par[3], event.helix_t[id_km][km_nh-1]);
      Double_t km_theta_max = event.helix_t[id_km][0] + vtx_scan_range/km_par[3];
      TVector3 km_start = TVector3(event.calpos_x[id_km][0], event.calpos_y[id_km][0], event.calpos_z[id_km][0]);
      TVector3 km_end = TVector3(event.calpos_x[id_km][km_nh-1], event.calpos_y[id_km][km_nh-1], event.calpos_z[id_km][km_nh-1]);

      if(id_kp!=-1){
	Double_t kp_par[5];
	kp_par[0] = event.helix_cx[id_kp];
	kp_par[1] = event.helix_cy[id_kp];
	kp_par[2] = event.helix_z0[id_kp];
	kp_par[3] = event.helix_r[id_kp];
	kp_par[4] = event.helix_dz[id_kp];
	Int_t kp_nh = event.helix_t[id_kp].size();
	Double_t kp_theta_min = event.helix_t[id_kp][0] - vtx_scan_range/kp_par[3];
	Double_t kp_theta_max = TMath::Min(event.helix_t[id_kp][0] + vtx_scan_rangeInsideL/kp_par[3], event.helix_t[id_kp][kp_nh-1]);
	TVector3 kp_start = TVector3(event.calpos_x[id_kp][0], event.calpos_y[id_kp][0], event.calpos_z[id_kp][0]);
	TVector3 kp_end = TVector3(event.calpos_x[id_kp][kp_nh-1], event.calpos_y[id_kp][kp_nh-1], event.calpos_z[id_kp][kp_nh-1]);

	Double_t kk_dist = 10000.;
	TVector3 kp_mom; TVector3 km_mom;
	TVector3 phi_mom;
	TVector3 phi_vert = Kinematics::LambdaVertex(dMagneticField, kp_par, km_par,
						     kp_theta_min, kp_theta_max,
						     km_theta_min, km_theta_max,
						     kp_mom, km_mom,
						     phi_mom, kk_dist);
	phi_mom = km_mom + kp_mom;

	TLorentzVector Lkp(kp_mom, TMath::Hypot(kp_mom.Mag(), KaonMass));
	TLorentzVector Lkm(km_mom, TMath::Hypot(km_mom.Mag(), KaonMass));
	TLorentzVector Lphi = Lkp + Lkm;

	Int_t flag = 1;
	Int_t temp = flag&event.pid[id_kp];
	if(temp == flag && !TMath::IsNaN(kk_dist) && TMath::Abs(phi_vert.x()) < 30. && TMath::Abs(phi_vert.z() - tpc::ZTarget) < 30. && TMath::Abs(phi_vert.y()) < 30.){
	  if(GFTrackCont.TrackCheck(id_km, 1) && GFTrackCont.TrackCheck(id_kp, 1)){
	    Double_t GFkk_dist = 10000.;
	    Double_t GFextrapolation_kk[2];
	    TVector3 GFmom_kk[2];
	    TVector3 GFphi_vert;
	    GFTrackCont.FindVertex(id_kp, id_km,
				   1, 1,
				   GFextrapolation_kk[0], GFextrapolation_kk[1],
				   GFmom_kk[0], GFmom_kk[1],
				   GFkk_dist, GFphi_vert,
				   vtx_scan_range);

	    TLorentzVector GFkm(GFmom_kk[0],
				TMath::Hypot(GFmom_kk[0].Mag(), KaonMass));
	    TLorentzVector GFkp(GFmom_kk[1],
				TMath::Hypot(GFmom_kk[1].Mag(), KaonMass));
	    TLorentzVector GFphi = GFkm + GFkp;
	    TVector3 GFmom_phi = GFmom_kk[0] + GFmom_kk[1];
	    TVector3 phivtx_dist(event.GFprodvtx_x - GFphi_vert.x(),
				 event.GFprodvtx_y - GFphi_vert.y(),
				 event.GFprodvtx_z - GFphi_vert.z());

	    event.lphiflag = true;
	    event.GFphimass = GFphi.M();
	    event.GFphidecayvtx_x = GFphi_vert.x();
	    event.GFphidecayvtx_y = GFphi_vert.y();
	    event.GFphidecayvtx_z = GFphi_vert.z();
	    event.GFphimom = GFmom_phi.Mag();
	    event.GFphimom_x = GFmom_phi.x();
	    event.GFphimom_y = GFmom_phi.y();
	    event.GFphimom_z = GFmom_phi.z();
	    event.GFkk_dist = GFkk_dist;
	    if(event.isgoodTPC[0]==1) event.GFphiprodvtx_dist = phivtx_dist.Mag();

	    event.phimass = Lphi.M();
	    event.phidecayvtx_x = phi_vert.x();
	    event.phidecayvtx_y = phi_vert.y();
	    event.phidecayvtx_z = phi_vert.z();
	    event.phimom = phi_mom.Mag();
	    event.phimom_x = phi_mom.x();
	    event.phimom_y = phi_mom.y();
	    event.phimom_z = phi_mom.z();
	    event.kk_dist = kk_dist;
	    event.phi_km_mass2 = km_mass2;

	    event.phidecays_id.push_back(id_km);
	    event.phidecays_id.push_back(id_km);
	    event.phidecays_mom.push_back(kp_mom.Mag());
	    event.phidecays_mom.push_back(km_mom.Mag());
	    event.phidecays_mom_x.push_back(kp_mom.x());
	    event.phidecays_mom_x.push_back(km_mom.x());
	    event.phidecays_mom_y.push_back(kp_mom.y());
	    event.phidecays_mom_y.push_back(km_mom.y());
	    event.phidecays_mom_z.push_back(kp_mom.z());
	    event.phidecays_mom_z.push_back(km_mom.z());

	    event.GFphidecays_mom.push_back(GFmom_kk[1].Mag());
	    event.GFphidecays_mom.push_back(GFmom_kk[0].Mag());
	    event.GFphidecays_mom_x.push_back(GFmom_kk[1].x());
	    event.GFphidecays_mom_x.push_back(GFmom_kk[0].x());
	    event.GFphidecays_mom_y.push_back(GFmom_kk[1].y());
	    event.GFphidecays_mom_y.push_back(GFmom_kk[0].y());
	    event.GFphidecays_mom_z.push_back(GFmom_kk[1].z());
	    event.GFphidecays_mom_z.push_back(GFmom_kk[0].z());
	  }
	}
      }
    } //lphiflag
  } //lflag

  //Remaining tracks
#if DebugDisp
  std::cout<<"Save remaining tracks"<<std::endl;
#endif
  std::vector<Double_t> twopi;
  event.pip_multi = target_pip_id_container.size();
  event.p_multi = target_p_id_container.size();
  event.pim_multi = target_pim_id_container.size();
  event.ppip_multi = target_ppip_id_container.size();
  event.accdient_multi = target_accidental_id_container.size();
  for(Int_t itp=0; itp<target_p_id_container.size(); ++itp){
    Int_t id_p = target_p_id_container[itp];
    if(id_p==GFtrackid_decays[0]){
      event.p_multi -= 1;
      continue;
    }
    event.residual_id.push_back(id_p);
    event.residual_mass2.push_back(target_p_mass2_container[itp]);
    event.residual_mom.push_back(target_p_mom_container[itp].Mag());
    event.residual_mom_x.push_back(target_p_mom_container[itp].x());
    event.residual_mom_y.push_back(target_p_mom_container[itp].y());
    event.residual_mom_z.push_back(target_p_mom_container[itp].z());
    event.residual_charge.push_back(1);
    if(event.lflag) HF1( 1130, target_p_mom_container[itp].Mag());
    else HF1( 1100, target_p_mom_container[itp].Mag());
  }
  for(Int_t itpip=0; itpip<target_pip_id_container.size(); ++itpip){
    Int_t id_pip = target_pip_id_container[itpip];
    if(id_pip==GFtrackid_decays[0]){
      event.pip_multi -= 1;
      continue;
    }
    event.residual_id.push_back(id_pip);
    event.residual_mass2.push_back(target_pip_mass2_container[itpip]);
    event.residual_mom.push_back(target_pip_mom_container[itpip].Mag());
    event.residual_mom_x.push_back(target_pip_mom_container[itpip].x());
    event.residual_mom_y.push_back(target_pip_mom_container[itpip].y());
    event.residual_mom_z.push_back(target_pip_mom_container[itpip].z());
    event.residual_charge.push_back(1);
    if(event.lflag) HF1( 1131, target_pip_mom_container[itpip].Mag());
    else HF1( 1101, target_pip_mom_container[itpip].Mag());
  }
  for(Int_t itpim=0; itpim<target_pim_id_container.size(); ++itpim){
    Int_t id_pim = target_pim_id_container[itpim];
    if(id_pim==GFtrackid_decays[1]){
      event.pim_multi -= 1;
      continue;
    }
    double KE = 1000.*(TMath::Hypot(target_pim_mom_container[itpim].Mag(), PionMass) - PionMass); //MeV/c2
    if(KE < 60.) twopi.push_back(KE);
    event.residual_id.push_back(id_pim);
    event.residual_mass2.push_back(target_pim_mass2_container[itpim]);
    event.residual_mom.push_back(target_pim_mom_container[itpim].Mag());
    event.residual_mom_x.push_back(target_pim_mom_container[itpim].x());
    event.residual_mom_y.push_back(target_pim_mom_container[itpim].y());
    event.residual_mom_z.push_back(target_pim_mom_container[itpim].z());
    event.residual_charge.push_back(-1);
    if(event.lflag) HF1( 1132, target_pim_mom_container[itpim].Mag());
    else HF1( 1102, target_pim_mom_container[itpim].Mag());
  }
  for(Int_t itppip=0; itppip<target_ppip_id_container.size(); ++itppip){
    Int_t id_ppip = target_ppip_id_container[itppip];
    if(id_ppip==GFtrackid_decays[0]){
      event.ppip_multi -= 1;
      continue;
    }
    event.residual_id.push_back(id_ppip);
    event.residual_mass2.push_back(target_ppip_mass2_container[itppip]);
    event.residual_mom.push_back(target_ppip_mom_container[itppip].Mag());
    event.residual_mom_x.push_back(target_ppip_mom_container[itppip].x());
    event.residual_mom_y.push_back(target_ppip_mom_container[itppip].y());
    event.residual_mom_z.push_back(target_ppip_mom_container[itppip].z());
    event.residual_charge.push_back(event.charge[id_ppip]);
  }
  event.residual_multi = event.pip_multi + event.p_multi + event.pim_multi + event.ppip_multi + event.accdient_multi;

  if(event.lflag){
    HF1( 1030, event.pim_multi);
    HF1( 1031, event.p_multi);
    HF1( 1032, event.pip_multi);

    HF1( 130, -event.BE[0]);
    if(event.isgoodTPC[0] == 1) HF1( 132, -event.BETPC[0]);

    if(event.pim_multi>0){
      HF1( 134, -event.BE[0]);
      if(event.isgoodTPC[0] == 1) HF1( 136, -event.BETPC[0]);
    }
    if(event.pim_multi==1) event.lpiflag = true;
  }
  else if(twopi.size()>=2){
    for(int pi1 = 0;pi1<twopi.size();pi1++){
      Double_t prevKE = twopi[pi1];
      for(int pi2 = pi1+1;pi1<twopi.size();pi1++){
	Double_t currentKE = twopi[pi2];
	HF2( 501, TMath::Min(prevKE, currentKE), TMath::Max(prevKE, currentKE) );
	event.pipiflag = true;
      }
    }
    HF1( 140, -event.BE[0]);
    if(event.isgoodTPC[0] == 1) HF1( 142, -event.BETPC[0]);
  }
  else{
    HF1( 1000, event.pim_multi);
    HF1( 1001, event.p_multi);
    HF1( 1002, event.pip_multi);
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
  tree->Branch( "remain_nclTpc", &event.nclTpc );
#if SaveRawData
  tree->Branch( "remain_cluster_x", &event.remain_cluster_x );
  tree->Branch( "remain_cluster_y", &event.remain_cluster_y );
  tree->Branch( "remain_cluster_z", &event.remain_cluster_z );
  tree->Branch( "remain_cluster_de", &event.remain_cluster_de );
  tree->Branch( "remain_cluster_size", &event.remain_cluster_size );
  tree->Branch( "remain_cluster_layer", &event.remain_cluster_layer );
  tree->Branch( "remain_cluster_row_center", &event.remain_cluster_row_center );
  tree->Branch( "remain_cluster_mrow", &event.remain_cluster_mrow );
  tree->Branch( "remain_cluster_de_center", &event.remain_cluster_de_center );
  tree->Branch( "remain_cluster_x_center", &event.remain_cluster_x_center );
  tree->Branch( "remain_cluster_y_center", &event.remain_cluster_y_center );
  tree->Branch( "remain_cluster_z_center", &event.remain_cluster_z_center );
  tree->Branch( "remain_cluster_houghflag", &event.remain_cluster_houghflag );
#endif

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "pid", &event.pid );
  tree->Branch( "chisqr", &event.chisqr );
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
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
#endif
  tree->Branch( "helix_t", &event.helix_t );

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

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "utgtK18", &event.utgtK18);
  tree->Branch( "vtgtK18", &event.vtgtK18);

  tree->Branch( "isgoodTPCK18", &event.isgoodTPCK18);
  tree->Branch( "chisqrTPCK18", &event.chisqrTPCK18);
  tree->Branch( "qTPCK18", &event.qTPCK18);
  tree->Branch( "pTPCK18", &event.pTPCK18);
  tree->Branch( "xtgtTPCK18", &event.xtgtTPCK18);
  tree->Branch( "ytgtTPCK18", &event.ytgtTPCK18);
  tree->Branch( "utgtTPCK18", &event.utgtTPCK18);
  tree->Branch( "vtgtTPCK18", &event.vtgtTPCK18);
  tree->Branch( "thetaTPCK18", &event.thetaTPCK18);
  tree->Branch( "lhtofTPCK18", &event.lhtofTPCK18);
  tree->Branch( "xhtofTPCK18", &event.xhtofTPCK18);
  tree->Branch( "yhtofTPCK18", &event.yhtofTPCK18);
  tree->Branch( "lvpTPCK18", &event.lvpTPCK18);
  tree->Branch( "xvpTPCK18", &event.xvpTPCK18);
  tree->Branch( "yvpTPCK18", &event.yvpTPCK18);

  tree->Branch( "ntKurama", &event.ntKurama);
  tree->Branch( "pKurama", &event.pKurama);
  tree->Branch( "qKurama", &event.qKurama);
  tree->Branch( "chisqrKurama", &event.chisqrKurama);
  tree->Branch( "m2Kurama", &event.m2);
  tree->Branch( "xtgtKurama", &event.xtgtKurama);
  tree->Branch( "ytgtKurama", &event.ytgtKurama);
  tree->Branch( "utgtKurama", &event.utgtKurama);
  tree->Branch( "vtgtKurama", &event.vtgtKurama);
  tree->Branch( "thetaKurama", &event.thetaKurama);
  tree->Branch( "pathwcKurama", &event.pathwcKurama);
  tree->Branch( "xin", &event.xin);
  tree->Branch( "yin", &event.yin);
  tree->Branch( "zin", &event.zin);
  tree->Branch( "pxin", &event.pxin);
  tree->Branch( "pyin", &event.pyin);
  tree->Branch( "pzin", &event.pzin);
  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "kflagTPCKurama", &event.kflagTPCKurama);
  tree->Branch( "chisqrTPCKurama", &event.chisqrTPCKurama);
  tree->Branch( "pTPCKurama", &event.pTPCKurama);
  tree->Branch( "qTPCKurama", &event.qTPCKurama);
  tree->Branch( "m2TPCKurama", &event.m2TPCKurama);
  tree->Branch( "xtgtTPCKurama", &event.xtgtTPCKurama);
  tree->Branch( "ytgtTPCKurama", &event.ytgtTPCKurama);
  tree->Branch( "utgtTPCKurama", &event.utgtTPCKurama);
  tree->Branch( "vtgtTPCKurama", &event.vtgtTPCKurama);
  tree->Branch( "thetaTPCKurama", &event.thetaTPCKurama);
  tree->Branch( "pathTPCKurama", &event.pathTPCKurama);
  tree->Branch( "lhtofTPCKurama", &event.lhtofTPCKurama);
  tree->Branch( "xhtofTPCKurama", &event.xhtofTPCKurama);
  tree->Branch( "yhtofTPCKurama", &event.yhtofTPCKurama);
  tree->Branch( "lvpTPCKurama", &event.lvpTPCKurama);
  tree->Branch( "xvpTPCKurama", &event.xvpTPCKurama);
  tree->Branch( "yvpTPCKurama", &event.yvpTPCKurama);

  tree->Branch( "isgoodTPC", &event.isgoodTPC);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "vtxTPC", &event.vtxTPC);
  tree->Branch( "vtyTPC", &event.vtyTPC);
  tree->Branch( "vtzTPC", &event.vtzTPC);
  tree->Branch( "closeDistTPC", &event.closeDistTPC);
  tree->Branch( "MissMassTPC", &event.MissMassTPC);
  tree->Branch( "MissMassCorrTPC", &event.MissMassCorrTPC);
  tree->Branch( "MissMassCorrDETPC", &event.MissMassCorrDETPC);
  tree->Branch( "pOrgTPC", &event.pOrgTPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "pCorrTPC", &event.pCorrTPC);
  tree->Branch( "pCorrDETPC", &event.pCorrDETPC);
  tree->Branch("xbTPC", &event.xbTPC);
  tree->Branch("ybTPC", &event.ybTPC);
  tree->Branch("ubTPC", &event.ubTPC);
  tree->Branch("vbTPC", &event.vbTPC);
  tree->Branch("xsTPC", &event.xsTPC);
  tree->Branch("ysTPC", &event.ysTPC);
  tree->Branch("usTPC", &event.usTPC);
  tree->Branch("vsTPC", &event.vsTPC);

  tree->Branch("nKK", &event.nKK);
  tree->Branch("Kflag", &event.Kflag);
  tree->Branch("MissMass", &event.MissMass);
  tree->Branch("MissMassCorr", &event.MissMassCorr);
  tree->Branch("MissMassCorrDE", &event.MissMassCorrDE);
  tree->Branch("vtx", &event.vtx);
  tree->Branch("vty", &event.vty);
  tree->Branch("vtz", &event.vtz);
  tree->Branch("pOrg", &event.pOrg);
  tree->Branch("pCalc", &event.pCalc);
  tree->Branch("pCorr", &event.pCorr);
  tree->Branch("pCorrDE", &event.pCorrDE);
  tree->Branch("xb", &event.xb);
  tree->Branch("yb", &event.yb);
  tree->Branch("ub", &event.ub);
  tree->Branch("vb", &event.vb);
  tree->Branch("xs", &event.xs);
  tree->Branch("ys", &event.ys);
  tree->Branch("us", &event.us);
  tree->Branch("vs", &event.vs);

  tree->Branch("KmMom_x", &event.km_mom_x);
  tree->Branch("KmMom_y", &event.km_mom_y);
  tree->Branch("KmMom_z", &event.km_mom_z);
  tree->Branch("KpMom_x", &event.kp_mom_x);
  tree->Branch("KpMom_y", &event.kp_mom_y);
  tree->Branch("KpMom_z", &event.kp_mom_z);

  tree->Branch("BE", &event.BE);
  tree->Branch("BEcorrDE", &event.BEcorrDE);
  tree->Branch("BETPC", &event.BETPC);
  tree->Branch("BEcorrDETPC", &event.BEcorrDETPC);

  //Multi-track vertex
  tree->Branch("GFProductionVtx_x", &event.GFprodvtx_x);
  tree->Branch("GFProductionVtx_y", &event.GFprodvtx_y);
  tree->Branch("GFProductionVtx_z", &event.GFprodvtx_z);

  //track fitting results
  //for all tracks
  //tree->Branch("GFntTpc", &event.GFntTpc);
  //for decay particles
  tree->Branch("GFntDecays", &event.GFntdecays);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFtof", &event.GFtof);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFpdgcode", &event.GFpdgcode);

  tree->Branch("GFDecaysHtofId", &event.GFdecays_htofid);
  tree->Branch("GFDecaysTrackLen", &event.GFdecays_tracklen);
  tree->Branch("GFDecaysTof", &event.GFdecays_tof);
  tree->Branch("GFDecaysMassSquare", &event.GFdecays_mass2);
  tree->Branch("GFDecaysMom", &event.GFdecays_mom);
  tree->Branch("GFDecaysMom_x", &event.GFdecays_mom_x);
  tree->Branch("GFDecaysMom_y", &event.GFdecays_mom_y);
  tree->Branch("GFDecaysMom_z", &event.GFdecays_mom_z);
  tree->Branch("GFMomLoss", &event.GFmomloss);
  tree->Branch("GFELoss", &event.GFeloss);

  tree->Branch("Xiflag", &event.xiflag);
  tree->Branch("XiMass", &event.ximass);
  tree->Branch("XiDecayVtx_x", &event.xidecayvtx_x);
  tree->Branch("XiDecayVtx_y", &event.xidecayvtx_y);
  tree->Branch("XiDecayVtx_z", &event.xidecayvtx_z);
  tree->Branch("XiMom", &event.ximom);
  tree->Branch("XiMom_x", &event.ximom_x);
  tree->Branch("XiMom_y", &event.ximom_y);
  tree->Branch("XiMom_z", &event.ximom_z);
  tree->Branch("XiVtxCloseDist", &event.lpi_dist);
  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom", &event.lmom);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  tree->Branch("XiTarget_x", &event.xitargetvtx_x);
  tree->Branch("XiTarget_y", &event.xitargetvtx_y);
  tree->Branch("XiTarget_z", &event.xitargetvtx_z);
  tree->Branch("XiTargetMom", &event.xitargetmom);
  tree->Branch("XiTargetMom_x", &event.xitargetmom_x);
  tree->Branch("XiTargetMom_y", &event.xitargetmom_y);
  tree->Branch("XiTargetMom_z", &event.xitargetmom_z);
  tree->Branch("XiTargetCloseDist", &event.xitarget_dist);

  tree->Branch("GFXiMass", &event.GFximass);
  tree->Branch("GFXiDecayVtx_x", &event.GFxidecayvtx_x);
  tree->Branch("GFXiDecayVtx_y", &event.GFxidecayvtx_y);
  tree->Branch("GFXiDecayVtx_z", &event.GFxidecayvtx_z);
  tree->Branch("GFXiMom", &event.GFximom);
  tree->Branch("GFXiMom_x", &event.GFximom_x);
  tree->Branch("GFXiMom_y", &event.GFximom_y);
  tree->Branch("GFXiMom_z", &event.GFximom_z);
  tree->Branch("GFXiVtxCloseDist", &event.GFlpi_dist);

  tree->Branch("KFXiMass", &event.KFximass);//Test branch; Xi InvMass with KF Lambda + GF pi
  tree->Branch("KFXiVtxCloseDist", &event.KFlpi_dist);
  //extrapolation
  //to the K, K vertex
  tree->Branch("GFXiKKVtx_x", &event.GFxikkvtx_x);
  tree->Branch("GFXiKKVtx_y", &event.GFxikkvtx_y);
  tree->Branch("GFXiKKVtx_z", &event.GFxikkvtx_z);
  tree->Branch("GFXiKKVtxMom", &event.GFxikkmom);
  tree->Branch("GFXiKKVtxMom_x", &event.GFxikkmom_x);
  tree->Branch("GFXiKKVtxMom_y", &event.GFxikkmom_y);
  tree->Branch("GFXiKKVtxMom_z", &event.GFxikkmom_z);
  tree->Branch("GFXiKKVtxCloseDist", &event.GFxikkvtx_dist);
  //to the production vertex
  tree->Branch("GFXiProductionVtx_x", &event.GFxiprodvtx_x);
  tree->Branch("GFXiProductionVtx_y", &event.GFxiprodvtx_y);
  tree->Branch("GFXiProductionVtx_z", &event.GFxiprodvtx_z);
  tree->Branch("GFXiProductionVtxMom", &event.GFxiprodmom);
  tree->Branch("GFXiProductionVtxMom_x", &event.GFxiprodmom_x);
  tree->Branch("GFXiProductionVtxMom_y", &event.GFxiprodmom_y);
  tree->Branch("GFXiProductionVtxMom_z", &event.GFxiprodmom_z);
  tree->Branch("GFXiProductionVtxCloseDist", &event.GFxiprodvtx_dist);
  tree->Branch("GFXiTrackLen", &event.GFxitracklen);
  tree->Branch("GFXiTof", &event.GFxitof);
  tree->Branch("GFXiMomLoss", &event.GFximomloss);
  tree->Branch("GFXiExcitation", &event.GFxiexcitation);
  //closest point to the target
  tree->Branch("GFXiTarget_x", &event.GFxitargetvtx_x);
  tree->Branch("GFXiTarget_y", &event.GFxitargetvtx_y);
  tree->Branch("GFXiTarget_z", &event.GFxitargetvtx_z);
  tree->Branch("GFXiTargetMom", &event.GFxitargetmom);
  tree->Branch("GFXiTargetMom_x", &event.GFxitargetmom_x);
  tree->Branch("GFXiTargetMom_y", &event.GFxitargetmom_y);
  tree->Branch("GFXiTargetMom_z", &event.GFxitargetmom_z);
  tree->Branch("GFXiTargetCloseDist", &event.GFxitarget_dist);
  //at z=z_target
  tree->Branch("GFXiTargetCenter_x", &event.GFxitargetcenter_x);
  tree->Branch("GFXiTargetCenter_y", &event.GFxitargetcenter_y);
  tree->Branch("GFXiTargetCenter_z", &event.GFxitargetcenter_z);
  tree->Branch("GFXiTargetCenterMom", &event.GFxitargetcentermom);
  tree->Branch("GFXiTargetCenterMom_x", &event.GFxitargetcentermom_x);
  tree->Branch("GFXiTargetCenterMom_y", &event.GFxitargetcentermom_y);
  tree->Branch("GFXiTargetCenterMom_z", &event.GFxitargetcentermom_z);
  tree->Branch("GFXiTargetCenterCloseDist", &event.GFxitargetcenter_dist);

  tree->Branch("LLflag", &event.llflag);
  tree->Branch("LambdaLambdaVtx_x", &event.llvtx_x);
  tree->Branch("LambdaLambdaVtx_y", &event.llvtx_y);
  tree->Branch("LambdaLambdaVtx_z", &event.llvtx_z);
  tree->Branch("LambdaLambdaCloseDist", &event.lldist);

  tree->Branch("LambdaMass1", &event.lmass1);
  tree->Branch("LambdaDecayVtx_x1", &event.ldecayvtx_x1);
  tree->Branch("LambdaDecayVtx_y1", &event.ldecayvtx_y1);
  tree->Branch("LambdaDecayVtx_z1", &event.ldecayvtx_z1);
  tree->Branch("LambdaMom1", &event.lmom1);
  tree->Branch("LambdaMom_x1", &event.lmom_x1);
  tree->Branch("LambdaMom_y1", &event.lmom_y1);
  tree->Branch("LambdaMom_z1", &event.lmom_z1);
  tree->Branch("LambdaVtxCloseDist1", &event.ppi_dist1);
  tree->Branch("LambdaTargetCloseDist1", &event.ltarget_dist1);
  tree->Branch("LambdaTargetCloseVtx_x1", &event.ltargetvtx_x1);
  tree->Branch("LambdaTargetCloseVtx_y1", &event.ltargetvtx_y1);
  tree->Branch("LambdaTargetCloseVtx_z1", &event.ltargetvtx_z1);

  tree->Branch("LambdaMass2", &event.lmass2);
  tree->Branch("LambdaDecayVtx_x2", &event.ldecayvtx_x2);
  tree->Branch("LambdaDecayVtx_y2", &event.ldecayvtx_y2);
  tree->Branch("LambdaDecayVtx_z2", &event.ldecayvtx_z2);
  tree->Branch("LambdaMom2", &event.lmom2);
  tree->Branch("LambdaMom_x2", &event.lmom_x2);
  tree->Branch("LambdaMom_y2", &event.lmom_y2);
  tree->Branch("LambdaMom_z2", &event.lmom_z2);
  tree->Branch("LambdaVtxCloseDist2", &event.ppi_dist2);
  tree->Branch("LambdaTargetCloseDist2", &event.ltarget_dist2);
  tree->Branch("LambdaTargetCloseVtx_x2", &event.ltargetvtx_x2);
  tree->Branch("LambdaTargetCloseVtx_y2", &event.ltargetvtx_y2);
  tree->Branch("LambdaTargetCloseVtx_z2", &event.ltargetvtx_z2);

  tree->Branch("GFLambdaMass1", &event.GFlmass1);
  tree->Branch("GFLambdaDecayVtx_x1", &event.GFldecayvtx_x1);
  tree->Branch("GFLambdaDecayVtx_y1", &event.GFldecayvtx_y1);
  tree->Branch("GFLambdaDecayVtx_z1", &event.GFldecayvtx_z1);
  tree->Branch("GFLambdaMom1", &event.GFlmom1);
  tree->Branch("GFLambdaMom_x1", &event.GFlmom_x1);
  tree->Branch("GFLambdaMom_y1", &event.GFlmom_y1);
  tree->Branch("GFLambdaMom_z1", &event.GFlmom_z1);
  tree->Branch("GFLambdaVtxCloseDist1", &event.GFppi_dist1);

  tree->Branch("GFLambdaMass2", &event.GFlmass2);
  tree->Branch("GFLambdaDecayVtx_x2", &event.GFldecayvtx_x2);
  tree->Branch("GFLambdaDecayVtx_y2", &event.GFldecayvtx_y2);
  tree->Branch("GFLambdaDecayVtx_z2", &event.GFldecayvtx_z2);
  tree->Branch("GFLambdaMom2", &event.GFlmom2);
  tree->Branch("GFLambdaMom_x2", &event.GFlmom_x2);
  tree->Branch("GFLambdaMom_y2", &event.GFlmom_y2);
  tree->Branch("GFLambdaMom_z2", &event.GFlmom_z2);
  tree->Branch("GFLambdaVtxCloseDist2", &event.GFppi_dist2);

  //L, L vertex
  tree->Branch("GFLambdaLambdaCloseDist", &event.GFlldist);
  tree->Branch("GFLambdaLambdaVtx_x", &event.GFllvtx_x);
  tree->Branch("GFLambdaLambdaVtx_y", &event.GFllvtx_y);
  tree->Branch("GFLambdaLambdaVtx_z", &event.GFllvtx_z);

  //extrapolation
  //to closest point to the target
  tree->Branch("GFLambdaTarget_x1", &event.GFltargetvtx_x1);
  tree->Branch("GFLambdaTarget_y1", &event.GFltargetvtx_y1);
  tree->Branch("GFLambdaTarget_z1", &event.GFltargetvtx_z1);
  tree->Branch("GFLambdaTargetCloseDist1", &event.GFltarget_dist1);
  tree->Branch("GFLambdaTarget_x2", &event.GFltargetvtx_x2);
  tree->Branch("GFLambdaTarget_y2", &event.GFltargetvtx_y2);
  tree->Branch("GFLambdaTarget_z2", &event.GFltargetvtx_z2);
  tree->Branch("GFLambdaTargetCloseDist2", &event.GFltarget_dist2);
  //at z=z_target
  tree->Branch("GFLambdaTargetCenter_x1", &event.GFltargetcenter_x1);
  tree->Branch("GFLambdaTargetCenter_y1", &event.GFltargetcenter_y1);
  tree->Branch("GFLambdaTargetCenter_z1", &event.GFltargetcenter_z1);
  tree->Branch("GFLambdaTargetCenterCloseDist1", &event.GFltargetcenter_dist1);
  tree->Branch("GFLambdaTargetCenter_x2", &event.GFltargetcenter_x2);
  tree->Branch("GFLambdaTargetCenter_y2", &event.GFltargetcenter_y2);
  tree->Branch("GFLambdaTargetCenter_z2", &event.GFltargetcenter_z2);
  tree->Branch("GFLambdaTargetCenterCloseDist2", &event.GFltargetcenter_dist2);
  //to the production vertex
  tree->Branch("GFLambdaProductionVtx_x1", &event.GFlprodvtx_x1);
  tree->Branch("GFLambdaProductionVtx_y1", &event.GFlprodvtx_y1);
  tree->Branch("GFLambdaProductionVtx_z1", &event.GFlprodvtx_z1);
  tree->Branch("GFLambdaProductionVtxCloseDist1", &event.GFlprodvtx_dist1);
  tree->Branch("GFLambdaTrackLen1", &event.GFltracklen1);
  tree->Branch("GFLambdaTof1", &event.GFltof1);
  tree->Branch("GFLambdaProductionVtx_x2", &event.GFlprodvtx_x2);
  tree->Branch("GFLambdaProductionVtx_y2", &event.GFlprodvtx_y2);
  tree->Branch("GFLambdaProductionVtx_z2", &event.GFlprodvtx_z2);
  tree->Branch("GFLambdaProductionVtxCloseDist2", &event.GFlprodvtx_dist2);
  tree->Branch("GFLambdaTrackLen2", &event.GFltracklen2);
  tree->Branch("GFLambdaTof2", &event.GFltof2);

  tree->Branch("Lflag", &event.lflag);
  tree->Branch("LambdaTargetCloseDist", &event.ltarget_dist);
  tree->Branch("LambdaTargetCloseVtx_x", &event.ltargetvtx_x);
  tree->Branch("LambdaTargetCloseVtx_y", &event.ltargetvtx_y);
  tree->Branch("LambdaTargetCloseVtx_z", &event.ltargetvtx_z);

  tree->Branch("GFLambdaMass", &event.GFlmass);
  tree->Branch("GFLambdaDecayVtx_x", &event.GFldecayvtx_x);
  tree->Branch("GFLambdaDecayVtx_y", &event.GFldecayvtx_y);
  tree->Branch("GFLambdaDecayVtx_z", &event.GFldecayvtx_z);
  tree->Branch("GFLambdaMom", &event.GFlmom);
  tree->Branch("GFLambdaMom_x", &event.GFlmom_x);
  tree->Branch("GFLambdaMom_y", &event.GFlmom_y);
  tree->Branch("GFLambdaMom_z", &event.GFlmom_z);
  tree->Branch("GFLambdaVtxCloseDist", &event.GFppi_dist);
  tree->Branch("GFLambdaTargetCloseDist", &event.GFltarget_dist);
  tree->Branch("GFLambdaTarget_x", &event.GFltargetvtx_x);
  tree->Branch("GFLambdaTarget_y", &event.GFltargetvtx_y);
  tree->Branch("GFLambdaTarget_z", &event.GFltargetvtx_z);
  tree->Branch("GFLambdaTargetCenter_x", &event.GFltargetcenter_x);
  tree->Branch("GFLambdaTargetCenter_y", &event.GFltargetcenter_y);
  tree->Branch("GFLambdaTargetCenter_z", &event.GFltargetcenter_z);
  tree->Branch("GFLambdaTargetCenterCloseDist", &event.GFltargetcenter_dist);
  tree->Branch("GFLambdaProductionVtx_x", &event.GFlprodvtx_x);
  tree->Branch("GFLambdaProductionVtx_y", &event.GFlprodvtx_y);
  tree->Branch("GFLambdaProductionVtx_z", &event.GFlprodvtx_z);
  tree->Branch("GFLambdaProductionVtxCloseDist", &event.GFlprodvtx_dist);
  tree->Branch("GFLambdaTrackLen", &event.GFltracklen);
  tree->Branch("GFLambdaTof", &event.GFltof);

  tree->Branch("LPhiflag", &event.lphiflag);
  tree->Branch("PhiMass", &event.phimass);
  tree->Branch("PhiDecayVtx_x", &event.phidecayvtx_x);
  tree->Branch("PhiDecayVtx_y", &event.phidecayvtx_y);
  tree->Branch("PhiDecayVtx_z", &event.phidecayvtx_z);
  tree->Branch("PhiMom", &event.phimom);
  tree->Branch("PhiMom_x", &event.phimom_x);
  tree->Branch("PhiMom_y", &event.phimom_y);
  tree->Branch("PhiMom_z", &event.phimom_z);
  tree->Branch("PhiVtxCloseDist", &event.kk_dist);
  tree->Branch("PhiDecaysTrackId", &event.phidecays_id);
  tree->Branch("PhiDecaysMom", &event.phidecays_mom);
  tree->Branch("PhiDecaysMom_x", &event.phidecays_mom_x);
  tree->Branch("PhiDecaysMom_y", &event.phidecays_mom_y);
  tree->Branch("PhiDecaysMom_z", &event.phidecays_mom_z);

  tree->Branch("GFPhiMass", &event.GFphimass);
  tree->Branch("GFPhiDecayVtx_x", &event.GFphidecayvtx_x);
  tree->Branch("GFPhiDecayVtx_y", &event.GFphidecayvtx_y);
  tree->Branch("GFPhiDecayVtx_z", &event.GFphidecayvtx_z);
  tree->Branch("GFPhiMom", &event.GFphimom);
  tree->Branch("GFPhiMom_x", &event.GFphimom_x);
  tree->Branch("GFPhiMom_y", &event.GFphimom_y);
  tree->Branch("GFPhiMom_z", &event.GFphimom_z);
  tree->Branch("GFPhiVtxCloseDist", &event.GFkk_dist);
  tree->Branch("GFPhiProductionVtxCloseDist", &event.GFphiprodvtx_dist);
  tree->Branch("KmMass2", &event.phi_km_mass2);
  tree->Branch("GFPhiDecaysMom", &event.GFphidecays_mom);
  tree->Branch("GFPhiDecaysMom_x", &event.GFphidecays_mom_x);
  tree->Branch("GFPhiDecaysMom_y", &event.GFphidecays_mom_y);
  tree->Branch("GFPhiDecaysMom_z", &event.GFphidecays_mom_z);

  tree->Branch("LPiflag", &event.lpiflag);

  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);

  tree->Branch("PiPiflag", &event.pipiflag);

  //Remaining p, pi after Xi, L searching
  //Multiplicity means tracks comes from the target
  tree->Branch("ResidualsMultiplicity", &event.residual_multi);
  tree->Branch("pipMultiplicity", &event.pip_multi);
  tree->Branch("pMultiplicity", &event.p_multi);
  tree->Branch("pimMultiplicity", &event.pim_multi);
  tree->Branch("ppipMultiplicity", &event.ppip_multi);
  tree->Branch("AccidentalMultiplicity", &event.accdient_multi);
  tree->Branch("ResidualsTrackId", &event.residual_id);
  tree->Branch("ResidualsMassSquare", &event.residual_mass2);
  tree->Branch("ResidualsMom", &event.residual_mom);
  tree->Branch("ResidualsMom_x", &event.residual_mom_x);
  tree->Branch("ResidualsMom_y", &event.residual_mom_y);
  tree->Branch("ResidualsMom_z", &event.residual_mom_z);
  tree->Branch("ResidualsCharge", &event.residual_charge);

#if DoKF
  tree->Branch("KFLambdaChisqr", &event.KFlchisqr);
  tree->Branch("KFLambdaPval", &event.KFlpval);
  tree->Branch("KFXiChisqr", &event.KFxichisqr);
  tree->Branch("KFXiPval", &event.KFxipval);
  tree->Branch("KFXiMass",&event.KFximass);

  tree->Branch("KFLambdaPull",&event.KFlpull);
  tree->Branch("KFXiPull",&event.KFxipull);
#endif

  TTreeReaderCont[kE42] = new TTreeReader( "tpc", TFileCont[kE42] );
  const auto& reader = TTreeReaderCont[kE42];
  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );

  src.nhHtof = new TTreeReaderValue<Int_t>( *reader, "nhHtof" );
  src.HtofSeg = new TTreeReaderValue<std::vector<Double_t>>( *reader, "HtofSeg" );
  src.tHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "tHtof" );
  src.dtHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dtHtof" );
  src.deHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "deHtof" );
  src.posHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "posHtof" );

  src.ntK18 = new TTreeReaderValue<Int_t>( *reader, "ntK18" );
  src.pK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pK18" );
  src.chisqrK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrK18" );
  src.xtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtK18" );
  src.ytgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtK18" );
  src.utgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtK18" );
  src.vtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtK18" );

  src.ntKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.chisqrKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrKurama" );
  src.pKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pKurama" );
  src.qKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qKurama" );
  src.m2 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2" );
  src.thetaKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaKurama" );
  src.xtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );
  src.pathwcKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pathwcKurama" );
  src.xin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xin" );
  src.yin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yin" );
  src.zin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "zin" );
  src.pxin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pxin" );
  src.pyin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pyin" );
  src.pzin = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pzin" );
  src.xout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xout" );
  src.yout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yout" );
  src.zout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "zout" );
  src.pxout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pxout" );
  src.pyout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pyout" );
  src.pzout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pzout" );

  src.nKm = new TTreeReaderValue<Int_t>( *reader, "nKm" );
  src.nKp = new TTreeReaderValue<Int_t>( *reader, "nKp" );
  src.nKK = new TTreeReaderValue<Int_t>( *reader, "nKK" );
  src.inside = new TTreeReaderValue<std::vector<Int_t>>( *reader, "inside" );
  src.vtx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx" );
  src.vty = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vty" );
  src.vtz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtz" );
  src.closeDist = new TTreeReaderValue<std::vector<Double_t>>( *reader, "closeDist" );
  src.MissMass = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMass" );
  src.MissMassCorr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorr" );
  src.MissMassCorrDE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrDE" );

  src.pOrg = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pOrg" );
  src.pCalc = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCalc" );
  src.pCorr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorr" );
  src.pCorrDE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrDE" );
  src.xb = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xb" );
  src.yb = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yb" );
  src.ub = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ub" );
  src.vb = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vb" );
  src.xs = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xs" );
  src.ys = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ys" );
  src.us = new TTreeReaderValue<std::vector<Double_t>>( *reader, "us" );
  src.vs = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vs" );
  src.Kflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "Kflag" );

  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" );
  src.remain_nclTpc = new TTreeReaderValue<Int_t>( *reader, "remain_nclTpc" );
#if SaveRawData
  src.cluster_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x" );
  src.cluster_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y" );
  src.cluster_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z" );
  src.cluster_de = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de" );
  src.cluster_size = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_size" );
  src.cluster_layer = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_layer" );
  src.cluster_mrow = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_mrow" );
  src.cluster_de_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_de_center" );
  src.cluster_x_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_x_center" );
  src.cluster_y_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_y_center" );
  src.cluster_z_center = new TTreeReaderValue<std::vector<Double_t>>( *reader, "cluster_z_center" );
  src.cluster_row_center = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_row_center" );
  src.cluster_houghflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "cluster_houghflag" );
#endif

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.isMultiloop = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isMultiloop" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );
  src.chisqr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr" );
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
  src.mom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_x" );
  src.mom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_y" );
  src.mom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "mom_z" );
  src.residual = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual" );
  src.residual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_x" );
  src.residual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_y" );
  src.residual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "residual_z" );
  src.resolution_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_x" );
  src.resolution_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_y" );
  src.resolution_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "resolution_z" );
  src.helix_t = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "helix_t" );
  src.pathhit = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "pathhit" );
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "alpha" );
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );

  src.nvtxTpc = new TTreeReaderValue<Int_t>(*reader,"nvtxTpc");
  src.vtx_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_x" );
  src.vtx_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_y" );
  src.vtx_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_z" );
  src.vtx_dist = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_dist" );
  src.vtx_angle = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx_angle" );
  src.vtxid = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxid" );
  src.vtxmom_theta = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_theta" );
  src.vtxpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_x" );
  src.vtxpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_y" );
  src.vtxpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxpos_z" );
  src.vtxmom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_x" );
  src.vtxmom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_y" );
  src.vtxmom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxmom_z" );

  src.isLambda = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isLambda" );
  src.ncombiLambda = new TTreeReaderValue<std::vector<Int_t>>( *reader, "ncombiLambda" );
  src.distLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "distLambda" );
  src.angleLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "angleLambda" );
  src.bestmassLambda = new TTreeReaderValue<std::vector<Double_t>>( *reader, "bestmassLambda" );
  src.massLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "massLambda" );
  src.vtxLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_x" );
  src.vtxLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_y" );
  src.vtxLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "vtxLambda_z" );
  src.momLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda" );
  src.momLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_x" );
  src.momLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_y" );
  src.momLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "momLambda_z" );
  src.decaysidLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysidLambda" );
  src.decaysmomLambda = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda" );
  src.decaysmomLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_x" );
  src.decaysmomLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_y" );
  src.decaysmomLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "decaysmomLambda_z" );

#if SaveTPCK18
  src.isgoodTPCK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCK18" );
  src.chisqrTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTPCK18" );
  src.qTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCK18");
  src.pTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCK18");
  src.xtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtTPCK18" );
  src.ytgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtTPCK18" );
  src.utgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtTPCK18" );
  src.vtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtTPCK18" );
  src.thetaTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPCK18" );
  src.lhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "lhtofTPCK18" );
  src.xhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofTPCK18" );
  src.yhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofTPCK18" );
  src.lvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "lvpTPCK18" );
  src.xvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpTPCK18" );
  src.yvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpTPCK18" );
  src.xhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofHS" );
  src.yhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofHS" );
#endif

  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
  src.kflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "kflagTPCKurama" );
  src.chisqrTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTPCKurama" );
  src.pTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCKurama" );
  src.qTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCKurama" );
  src.m2TPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2TPCKurama" );
  src.xtgtTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtTPCKurama" );
  src.ytgtTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtTPCKurama" );
  src.utgtTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtTPCKurama" );
  src.vtgtTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtTPCKurama" );
  src.thetaTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPCKurama" );
  src.pathTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pathTPCKurama" );
  src.lhtofTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "lhtofTPCKurama" );
  src.xhtofTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofTPCKurama" );
  src.yhtofTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofTPCKurama" );
  src.lvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "lvpTPCKurama" );
  src.xvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpTPCKurama" );
  src.yvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpTPCKurama" );

  src.isgoodTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPC" );
  src.insideTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "insideTPC" );
  src.vtxTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtxTPC" );
  src.vtyTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtyTPC" );
  src.vtzTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtzTPC" );
  src.closeDistTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "closeDistTPC" );
  src.MissMassTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassTPC" );
  src.MissMassCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrTPC" );
  src.MissMassCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrDETPC" );
  src.pOrgTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pOrgTPC" );
  src.pCalcTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCalcTPC" );
  src.pCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrTPC" );
  src.pCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrDETPC" );
  src.thetaTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPC" );
  src.thetaCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaCMTPC" );
  src.costCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "costCMTPC" );
  src.xbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xbTPC" );
  src.ybTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ybTPC" );
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );
  src.xsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xsTPC" );
  src.ysTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ysTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );

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
