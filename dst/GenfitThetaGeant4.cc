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
//const Int_t verbosity = 1;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 150.;
const Double_t vtx_scan_rangeInsidePi = 50.;

const Double_t ks_masscut = 0.2; //ref
//const Double_t xi_masscut = 0.25; const Double_t lambda_masscut = 0.1; //ref
//const Double_t xi_masscut = 0.03; const Double_t lambda_masscut = 0.03;
const Double_t pip_vtx_distcut = 300;
const Double_t pim_vtx_distcut = 300;
const Double_t p_vtx_distcut = 300;
const Double_t pipi_distcut = 10.; //ref
//const Double_t lpi_distcut = 10.; //ref
const Double_t ksp_distcut = 15.; //ref

const Double_t GFpipi_distcut = 10.;
const Double_t GFksp_distcut = 15.;
//const Double_t GFlpi_distcut = 15.;



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

	Int_t inc;


  int G4thetaid;
  double G4thetavtx_x;//Production vertex
  double G4thetavtx_y;
  double G4thetavtx_z;
  double G4thetamom;
  double G4thetamom_x;//Momentum at production vtx
  double G4thetamom_y;
  double G4thetamom_z;
  double G4ksp1mass;
  
  double G4ksp2mass;

  int G4ksid;
  double G4ksvtx_x;//Production vertex
  double G4ksvtx_y;
  double G4ksvtx_z;
  double G4ksmom;
  double G4ksmom_x;
  double G4ksmom_y;
  double G4ksmom_z;
  double G4ksmass;

  int G4p1id;//Id of Geant4 Track
  int G4p1tid;//G4 Id of TPC Track
  int G4p1nh;//Number of G4 Hits
  int G4p1tnh;//Number of correct G4 Hits
  double G4p1vtx_x;//Production vertex
  double G4p1vtx_y;
  double G4p1vtx_z;
  double G4p1mom;
  double G4p1mom_x;
  double G4p1mom_y;
  double G4p1mom_z;

  int p1tid;
  int p1nh;
  double p1vtx_x;
  double p1vtx_y;
  double p1vtx_z;
  double p1mom;
  double p1mom_x;
  double p1mom_y;
  double p1mom_z;
  double GFp1mom;
  double GFp1mom_x;
  double GFp1mom_y;
  double GFp1mom_z;
  
  int G4p2id;//Id of Geant4 Track
  int G4p2tid;//G4 Id of TPC Track
  int G4p2nh;//Number of G4 Hits
  int G4p2tnh;//Number of correct G4 Hits
  double G4p2vtx_x;//Production vertex
  double G4p2vtx_y;
  double G4p2vtx_z;
  double G4p2mom;
  double G4p2mom_x;
  double G4p2mom_y;
  double G4p2mom_z;

  int p2tid;
  int p2nh;
  double p2vtx_x;
  double p2vtx_y;
  double p2vtx_z;
  double p2mom;
  double p2mom_x;
  double p2mom_y;
  double p2mom_z;
  double GFp2mom;
  double GFp2mom_x;
  double GFp2mom_y;
  double GFp2mom_z;


  int G4pimid;
  int G4pimtid;
  int G4pimnh;
  int G4pimtnh;
  double G4pimvtx_x;//Production vertex, identical to Ks decay vertex
  double G4pimvtx_y;
  double G4pimvtx_z;
  double G4pimmom;
  double G4pimmom_x;
  double G4pimmom_y;
  double G4pimmom_z;

  int pimtid;
  int pimnh;
  double pimvtx_x;
  double pimvtx_y;
  double pimvtx_z;
  double pimmom;
  double pimmom_x;
  double pimmom_y;
  double pimmom_z;
  double GFpimmom;
  double GFpimmom_x;
  double GFpimmom_y;
  double GFpimmom_z;
  double KFpimmom;
  double KFpimmom_x;
  double KFpimmom_y;
  double KFpimmom_z;


  int G4pipid;
  int G4piptid;
  int G4pipnh;
  int G4piptnh;
  double G4pipvtx_x;
  double G4pipvtx_y;
  double G4pipvtx_z;
  double G4pipmom;
  double G4pipmom_x;
  double G4pipmom_y;
  double G4pipmom_z;

  int piptid;
  int pipnh;
  double pipvtx_x;
  double pipvtx_y;
  double pipvtx_z;
  double pipmom;
  double pipmom_x;
  double pipmom_y;
  double pipmom_z;
  double GFpipmom;
  double GFpipmom_x;
  double GFpipmom_y;
  double GFpipmom_z;
  double KFpipmom;
  double KFpipmom_x;
  double KFpipmom_y;
  double KFpipmom_z;


  bool Ksflag;
  int Kspipid;
  int Kspimid;
  double Ksmass;
  double Ksdecayvtx_x;
  double Ksdecayvtx_y;
  double Ksdecayvtx_z;
  double Ksmom;
  double Ksmom_x;
  double Ksmom_y;
  double Ksmom_z;
  double Ksdist;
  double KFpvalKs;
  std::vector<double> KFpullKs;
  double KFKsmom;
  double KFKsmom_x;
  double KFKsmom_y;
  double KFKsmom_z;

  int KsPcnt;
  bool KsP1flag;
  int KsP1pid;
  double KsP1mass;
  double KsP1decayvtx_x;
  double KsP1decayvtx_y;
  double KsP1decayvtx_z;
  double KsP1mom;
  double KsP1mom_x;
  double KsP1mom_y;
  double KsP1mom_z;
  double KsP1dist;

  bool KsP2flag;
  int KsP2pid;
  double KsP2mass;
  double KsP2decayvtx_x;
  double KsP2decayvtx_y;
  double KsP2decayvtx_z;
  double KsP2mom;
  double KsP2mom_x;
  double KsP2mom_y;
  double KsP2mom_z;
  double KsP2dist;

  bool ForcedKsflag;
  bool GFForcedKsflag;
  double ForcedKFpvalKs;
  std::vector<double> ForcedKFpullKs;
  double ForcedKsmass;
  double ForcedKsdecayvtx_x;
  double ForcedKsdecayvtx_y;
  double ForcedKsdecayvtx_z;
  double ForcedKsmom;
  double ForcedKsmom_x;
  double ForcedKsmom_y;
  double ForcedKsmom_z;
  double ForcedKsdist;
  double Forcedpipvertexdist;
  double Forcedpimvertexdist;
  double Forcedpipmom;
  double Forcedpipmom_x;
  double Forcedpipmom_y;
  double Forcedpipmom_z;
  double Forcedpimmom;
  double Forcedpimmom_x;
  double Forcedpimmom_y;
  double Forcedpimmom_z;
  
  double ForcedKFKsmom;
  double ForcedKFKsmom_x;
  double ForcedKFKsmom_y;
  double ForcedKFKsmom_z;
  double ForcedKFpipmom;
  double ForcedKFpipmom_x;
  double ForcedKFpipmom_y;
  double ForcedKFpipmom_z;
  double ForcedKFpimmom;
  double ForcedKFpimmom_x;
  double ForcedKFpimmom_y;
  double ForcedKFpimmom_z;


  bool ForcedKsP1flag;
  bool GFForcedKsP1flag;
  double ForcedKsP1mass;
  double ForcedKsP1decayvtx_x;
  double ForcedKsP1decayvtx_y;
  double ForcedKsP1decayvtx_z;
  double ForcedKsP1mom;
  double ForcedKsP1mom_x;
  double ForcedKsP1mom_y;
  double ForcedKsP1mom_z;
  double ForcedKsP1dist;
  double ForcedKsP1vertexdist;


  bool ForcedKsP2flag;
  bool GFForcedKsP2flag;
  double ForcedKsP2mass;
  double ForcedKsP2decayvtx_x;
  double ForcedKsP2decayvtx_y;
  double ForcedKsP2decayvtx_z;
  double ForcedKsP2mom;
  double ForcedKsP2mom_x;
  double ForcedKsP2mom_y;
  double ForcedKsP2mom_z;
  double ForcedKsP2dist;
  double ForcedKsP2vertexdist;
  
  
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
		inc = -1;



    G4thetaid = -1;
    G4thetavtx_x = qnan;
    G4thetavtx_y = qnan;
    G4thetavtx_z = qnan;
    G4thetamom = qnan;
    G4thetamom_x = qnan;
    G4thetamom_y = qnan;
    G4thetamom_z = qnan;
    G4ksp1mass = qnan;
    G4ksp2mass = qnan;

    G4ksid = -1;
    G4ksvtx_x = qnan;
    G4ksvtx_y = qnan;
    G4ksvtx_z = qnan;
    G4ksmom = qnan;
    G4ksmom_x = qnan;
    G4ksmom_y = qnan;
    G4ksmom_z = qnan;
    G4ksmass = qnan;

    G4p1id = -1;
    G4p1tid = -1;
    G4p1nh = 0;
    G4p1tnh = 0;
    G4p1vtx_x = qnan;
    G4p1vtx_y = qnan;
    G4p1vtx_z = qnan;
    G4p1mom = qnan;
    G4p1mom_x = qnan;
    G4p1mom_y = qnan;
    G4p1mom_z = qnan;
    
    p1tid = -1;
    p1nh = 0;
    p1vtx_x = qnan;
    p1vtx_y = qnan;
    p1vtx_z = qnan;
    p1mom = qnan;
    p1mom_x = qnan;
    p1mom_y = qnan;
    p1mom_z = qnan;
    GFp1mom = qnan;
    GFp1mom_x = qnan;
    GFp1mom_y = qnan;
    GFp1mom_z = qnan;

    G4p2id = -1;
    G4p2tid = -1;
    G4p2nh = 0;
    G4p2tnh = 0;
    G4p2vtx_x = qnan;
    G4p2vtx_y = qnan;
    G4p2vtx_z = qnan;
    G4p2mom = qnan;
    G4p2mom_x = qnan;
    G4p2mom_y = qnan;
    G4p2mom_z = qnan;
    
    p2tid = -1;
    p2nh = 0;
    p2vtx_x = qnan;
    p2vtx_y = qnan;
    p2vtx_z = qnan;
    p2mom = qnan;
    p2mom_x = qnan;
    p2mom_y = qnan;
    p2mom_z = qnan;
    GFp2mom = qnan;
    GFp2mom_x = qnan;
    GFp2mom_y = qnan;
    GFp2mom_z = qnan;

    G4pimid = -1;
    G4pimtid = -1;
    G4pimnh = 0;
    G4pimtnh = 0;
    G4pimvtx_x = qnan;
    G4pimvtx_y = qnan;
    G4pimvtx_z = qnan;
    G4pimmom = qnan;
    G4pimmom_x = qnan;
    G4pimmom_y = qnan;
    G4pimmom_z = qnan;
    
    pimtid = -1;
    pimnh = 0;
    pimvtx_x = qnan;
    pimvtx_y = qnan;
    pimvtx_z = qnan;
    pimmom = qnan;
    pimmom_x = qnan;
    pimmom_y = qnan;
    pimmom_z = qnan;
    GFpimmom = qnan;
    GFpimmom_x = qnan;
    GFpimmom_y = qnan;
    GFpimmom_z = qnan;
    KFpimmom = qnan;
    KFpimmom_x = qnan;
    KFpimmom_y = qnan;
    KFpimmom_z = qnan;

    G4pipid = -1;
    G4piptid = -1;
    G4pipnh = 0;
    G4piptnh = 0;
    G4pipvtx_x = qnan;
    G4pipvtx_y = qnan;
    G4pipvtx_z = qnan;
    G4pipmom = qnan;
    G4pipmom_x = qnan;
    G4pipmom_y = qnan;
    G4pipmom_z = qnan;

    piptid = -1;
    pipnh = 0;
    pipvtx_x = qnan;
    pipvtx_y = qnan;
    pipvtx_z = qnan;
    pipmom = qnan;
    pipmom_x = qnan;
    pipmom_y = qnan;
    pipmom_z = qnan;
    GFpipmom = qnan;
    GFpipmom_x = qnan;
    GFpipmom_y = qnan;
    GFpipmom_z = qnan;
    KFpipmom = qnan;
    KFpipmom_x = qnan;
    KFpipmom_y = qnan;
    KFpipmom_z = qnan;

    Ksflag = false;
    Kspipid = -1;
    Ksmass = qnan;
    Ksdecayvtx_x = qnan;
    Ksdecayvtx_y = qnan;
    Ksdecayvtx_z = qnan;
    Ksmom = qnan;
    Ksmom_x = qnan;
    Ksmom_y = qnan;
    Ksmom_z = qnan;
    Ksdist = qnan;

    KsPcnt = 0;

    KsP1flag = false;
    KsP1pid = -1;
    KsP1mass = qnan;
    KsP1decayvtx_x = qnan;
    KsP1decayvtx_y = qnan;
    KsP1decayvtx_z = qnan;
    KsP1mom = qnan;
    KsP1mom_x = qnan;
    KsP1mom_y = qnan;
    KsP1mom_z = qnan;
    KsP1dist = qnan;

    KsP2flag = false;
    KsP2pid = -1;
    KsP2mass = qnan;
    KsP2decayvtx_x = qnan;
    KsP2decayvtx_y = qnan;
    KsP2decayvtx_z = qnan;
    KsP2mom = qnan;
    KsP2mom_x = qnan;
    KsP2mom_y = qnan;
    KsP2mom_z = qnan;
    KsP2dist = qnan;
    
    ForcedKsflag = false;
    GFForcedKsflag = false;
    ForcedKsmass = qnan;
    ForcedKsdecayvtx_x = qnan;
    ForcedKsdecayvtx_y = qnan;
    ForcedKsdecayvtx_z = qnan;
    ForcedKsmom = qnan;
    ForcedKsmom_x = qnan;
    ForcedKsmom_y = qnan;
    ForcedKsmom_z = qnan;
    ForcedKsdist = qnan;
    Forcedpipvertexdist = qnan;
    Forcedpimvertexdist = qnan;
    Forcedpimmom = qnan;
    Forcedpimmom_x = qnan;
    Forcedpimmom_y = qnan;  
    Forcedpimmom_z = qnan;
    Forcedpipmom = qnan;
    Forcedpipmom_x = qnan;
    Forcedpipmom_y = qnan;  
    Forcedpipmom_z = qnan;
    
    ForcedKFKsmom = qnan;
    ForcedKFKsmom_x = qnan;
    ForcedKFKsmom_y = qnan;  
    ForcedKFKsmom_z = qnan;
    ForcedKFpimmom = qnan;
    ForcedKFpimmom_x = qnan;
    ForcedKFpimmom_y = qnan;  
    ForcedKFpimmom_z = qnan;
    ForcedKFpipmom = qnan;
    ForcedKFpipmom_x = qnan;
    ForcedKFpipmom_y = qnan;  
    ForcedKFpipmom_z = qnan;
    

    ForcedKsP1flag = false;
    GFForcedKsP1flag = false;
    ForcedKsP1mass = qnan;
    ForcedKsP1decayvtx_x = qnan;
    ForcedKsP1decayvtx_y = qnan;
    ForcedKsP1decayvtx_z = qnan;
    ForcedKsP1mom = qnan;
    ForcedKsP1mom_x = qnan;
    ForcedKsP1mom_y = qnan;
    ForcedKsP1mom_z = qnan;
    ForcedKsP1dist = qnan;
    ForcedKsP1vertexdist = qnan;
    
    ForcedKsP2flag = false;
    GFForcedKsP2flag = false;
    ForcedKsP2mass = qnan;
    ForcedKsP2decayvtx_x = qnan;
    ForcedKsP2decayvtx_y = qnan;
    ForcedKsP2decayvtx_z = qnan;
    ForcedKsP2mom = qnan;
    ForcedKsP2mom_x = qnan;
    ForcedKsP2mom_y = qnan;
    ForcedKsP2mom_z = qnan;
    ForcedKsP2dist = qnan;
    ForcedKsP2vertexdist = qnan;



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
  Double_t MomentumOfTrack_x[1000];
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

  static const auto KSMass    = 0.497662;
  static const auto ThetaMass = 1.524;
  static const auto PionMass    = pdg::PionMass();
  static const auto KaonMass    = pdg::KaonMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto me = 0.001*0.5109989461;
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

  HF1( 1, event.status++ );

  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;
  event.inc = src.inc;
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
    double mom = src.MomentumOfTrack[it]/1000;//MeV to GeV
    double mom_x = src.MomentumOfTrack_x[it]/1000;
    double mom_y = src.MomentumOfTrack_y[it]/1000;
    double mom_z = src.MomentumOfTrack_z[it]/1000;
    double vert_x = src.VertexOfTrack_x[it];
    double vert_y = src.VertexOfTrack_y[it];
    double vert_z = src.VertexOfTrack_z[it];
    if(pid == -211 and parent ==1){
      event.G4pimid = it;
      event.G4pimvtx_x = vert_x;
      event.G4pimvtx_y = vert_y;
      event.G4pimvtx_z = vert_z;
      event.G4pimmom = mom;
      event.G4pimmom_x = mom_x;
      event.G4pimmom_y = mom_y;
      event.G4pimmom_z = mom_z;

    }
    if(pid == 310){
      event.G4ksid = it;
      event.G4ksvtx_x = vert_x;
      event.G4ksvtx_y = vert_y;
      event.G4ksvtx_z = vert_z;
      event.G4ksmom = mom;
      event.G4ksmom_x = mom_x;
      event.G4ksmom_y = mom_y;
      event.G4ksmom_z = mom_z;
    }
    if(pid == 211 and parent ==1){
      event.G4pipid = it;
      event.G4pipvtx_x = vert_x;
      event.G4pipvtx_y = vert_y;
      event.G4pipvtx_z = vert_z;
      event.G4pipmom = mom;
      event.G4pipmom_x = mom_x;
      event.G4pipmom_y = mom_y;
      event.G4pipmom_z = mom_z;
    }
    if(pid==2212 and parent == 0 and event.G4p1id==-1){
      event.G4p1id = it;
      event.G4p1vtx_x = vert_x;
      event.G4p1vtx_y = vert_y;
      event.G4p1vtx_z = vert_z;
      event.G4p1mom = mom;
      event.G4p1mom_x = mom_x;
      event.G4p1mom_y = mom_y;
      event.G4p1mom_z = mom_z;
    }
    if(pid == 2212 and parent == 0 and event.G4p1id!=-1 and it != event.G4p1id){
      event.G4p2id = it;
      event.G4p2vtx_x = vert_x;
      event.G4p2vtx_y = vert_y;
      event.G4p2vtx_z = vert_z;
      event.G4p2mom = mom;
      event.G4p2mom_x = mom_x;
      event.G4p2mom_y = mom_y;
      event.G4p2mom_z = mom_z;
    }
  }
  vector<int> G4TrackID;
  vector<int> PureHits;

  event.G4pimnh = CountHits(event.G4pimid,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4pipnh = CountHits(event.G4pipid,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4p1nh = CountHits(event.G4p1id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);
  event.G4p2nh = CountHits(event.G4p2id,src.xtpc,src.ztpc,src.ititpc,src.nhittpc);

  TVector3 G4TVPim(event.G4pimmom_x,event.G4pimmom_y,event.G4pimmom_z);
  TVector3 G4TVPip(event.G4pipmom_x,event.G4pipmom_y,event.G4pipmom_z);
  TVector3 G4TVP1(event.G4p1mom_x,event.G4p1mom_y,event.G4p1mom_z);
  TVector3 G4TVP2(event.G4p2mom_x,event.G4p2mom_y,event.G4p2mom_z);
  TLorentzVector G4LVPim(G4TVPim,hypot(event.G4pimmom,PionMass));
  TLorentzVector G4LVPip(G4TVPip,hypot(event.G4pipmom,PionMass));
  TLorentzVector G4LVP1(G4TVP1,hypot(event.G4p1mom,ProtonMass));
  TLorentzVector G4LVP2(G4TVP2,hypot(event.G4p2mom,ProtonMass));

  event.G4ksmass = (G4LVPim+G4LVPip).M();
  event.G4ksp1mass = (G4LVPim+G4LVPip+G4LVP1).M();
  event.G4ksp2mass = (G4LVPim+G4LVPip+G4LVP2).M();
  
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
    if(event.G4p1id == G4tid){
      if(event.G4p1tid != -1){
        event.G4p1tid = it;
      }
      else if(event.efficiency[it]>event.efficiency[event.G4p1tid]){
        event.G4p1tid = it;
      }
    }
    if(event.G4p2id == G4tid){
      if(event.G4p2tid != -1){
        event.G4p2tid = it;
      }
      else if(event.efficiency[it]>event.efficiency[event.G4p2tid]){
        event.G4p2tid = it;
      }
    }
    if(event.G4pipid == G4tid){
      if(event.G4piptid != -1){
        event.G4piptid = it;
      }
      else if(event.efficiency[it]>event.efficiency[event.G4piptid]){
        event.G4piptid = it;
      }
    }
    if(event.G4pimid == G4tid){
      if(event.G4pimtid != -1){
        event.G4pimtid = it;
      }
      else if(event.efficiency[it]>event.efficiency[event.G4pimtid]){
        event.G4pimtid = it;
      }
    }
  }
  GFTrackCont.FitTracks();
  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }
  HF1( 1, event.status++ );


  //for Ks -> pi+ pi-
  std::vector<Int_t> Ks_pip_id_container, Ks_pim_id_container;
  std::vector<Int_t> Ks_pip_repid_container, Ks_pim_repid_container;
  std::vector<Int_t> Ks_pip_G4tid_container, Ks_pim_G4tid_container;
  std::vector<Double_t> Ks_pip_purity_container, Ks_pim_purity_container;
  std::vector<Double_t> Ks_pip_efficiency_container, Ks_pim_efficiency_container;
  std::vector<TVector3> Ks_pip_mom_container, Ks_pim_mom_container;
  std::vector<TVector3> Ks_KFpip_mom_container, Ks_KFpim_mom_container;
  std::vector<Double_t> Ks_mass_container;
  std::vector<Double_t> Ks_pval_container;
  std::vector<std::vector<Double_t>> Ks_pull_container;
  std::vector<Double_t> Ks_pipidist_container;
  std::vector<Double_t> Ks_targetdist_container;

  std::vector<TVector3> Ks_mom_container,Ks_KFmom_container, Ks_vtx_container, Ks_targetvtx_container;


  //ks/ksp candidates searching
  Int_t ks_candidates = 0;
  Int_t ksp_candidates = 0;
  std::vector<Int_t> ksp_ks_container;
  std::vector<Int_t> ksp_pim_container,ksp_pip_container,ksp_p_container;
  std::vector<Int_t> ksp_pip_repid_container,ksp_pim_repid_container,ksp_p_repid_container;
  std::vector<Int_t> ksp_pip_G4tid_container,ksp_pim_G4tid_container,ksp_p_G4tid_container;
  std::vector<Double_t> ksp_pip_purity_container,ksp_pim_purity_container,ksp_p_purity_container;
  std::vector<Double_t> ksp_pip_efficiency_container,ksp_pim_efficiency_container,ksp_p_efficiency_container;
  std::vector<TVector3> ksp_pip_mom_container,ksp_pim_mom_container,ksp_p_mom_container;
  std::vector<Double_t> ksp_mass_container;
  std::vector<Double_t> ksp_pipidist_container;
  std::vector<Double_t> ksp_targetdist_container;
  std::vector<TVector3> ksp_mom_container;
  std::vector<TVector3> ksp_vtx_container;
  std::vector<TVector3> ksp_targetvtx_container;
  std::vector<Double_t> ksp_pipidist;
  std::vector<Double_t> ksp_targetdist;
  std::vector<TVector3> ksp_targetvtx;
  std::vector<TVector3> ksp_vtx;
  std::vector<TVector3> ksp_mom;
  std::vector<Double_t> ksp_mass;
  std::vector<TVector3> ksp_pip_mom,ksp_pim_mom,ksp_p_mom;
  std::vector<Double_t> ksp_pip_purity,ksp_pim_purity,ksp_p_purity;
  std::vector<Double_t> ksp_p_efficiency,ksp_pip_efficiency,ksp_pim_efficiency;
  std::vector<Int_t> ksp_pip_id,ksp_pim_id,ksp_p_id;
  std::vector<Int_t> ksp_pip_G4tid,ksp_pim_G4tid,ksp_p_G4tid;
  std::vector<Int_t> ksp_pip_repid,ksp_pim_repid,ksp_p_repid;



  Int_t best_ks_id = -1;
  Double_t mass_diff = 9999.;
  {

    for(Int_t it1=0;it1<ntTpc;it1++){
    #if DebugDisp
      cout<<Form("Processing track %d / %d",it1,ntTpc)<<endl;
    #endif
      if((event.pid[it1]&1)!=1) continue; //select pi+like
      if(event.charge[it1]!=1) continue;

      Int_t repid_pip = 0;
      cout<<Form("Checking Track %d, repid = %d",it1,repid_pip)<<endl;
      if(!GFTrackCont.TrackCheck(it1, repid_pip)) continue;
      cout<<Form("Track %d, repid = %d checked",it1,repid_pip)<<endl;

      Double_t pip_par[5];
      pip_par[0] = event.helix_cx[it1];
      pip_par[1] = event.helix_cy[it1];
      pip_par[2] = event.helix_z0[it1];
      pip_par[3] = event.helix_r[it1];
      pip_par[4] = event.helix_dz[it1];
      Int_t pip_nh = event.helix_t[it1].size();
      Double_t pip_theta_min = event.helix_t[it1][0] - vtx_scan_range/pip_par[3];
      Double_t pip_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_range/pip_par[3], event.helix_t[it1][pip_nh-1]);
      TVector3 pip_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 pip_end = TVector3(event.calpos_x[it1][pip_nh-1], event.calpos_y[it1][pip_nh-1], event.calpos_z[it1][pip_nh-1]);

      for(Int_t it2=0;it2<ntTpc;it2++){
    #if DebugDisp
      cout<<Form("Processing track (%d,%d) / %d",it1,it2,ntTpc)<<endl;
    #endif
	if(it1==it2) continue;
	if((event.pid[it2]&1)!=1) continue; //select pi-like
	if(event.charge[it2]!=-1) continue;

	Int_t repid_pim = 0;
  
  cout<<Form("Checking Track %d, repid = %d",it2,repid_pim)<<endl;
	if(!GFTrackCont.TrackCheck(it2, repid_pim)) continue;
  cout<<Form("Track %d, repid = %d checked",it2,repid_pim)<<endl;

	Double_t pim_par[5];
	pim_par[0] = event.helix_cx[it2];
	pim_par[1] = event.helix_cy[it2];
	pim_par[2] = event.helix_z0[it2];
	pim_par[3] = event.helix_r[it2];
	pim_par[4] = event.helix_dz[it2];

	Int_t pim_nh = event.helix_t[it2].size();
	Double_t pim_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/pim_par[3], event.helix_t[it2][pim_nh-1]);
	Double_t pim_theta_max = event.helix_t[it2][0] + vtx_scan_range/pim_par[3];
	TVector3 pim_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 pim_end = TVector3(event.calpos_x[it2][pim_nh-1], event.calpos_y[it2][pim_nh-1], event.calpos_z[it2][pim_nh-1]);

	Double_t pipi_dist = 10000.;
	TVector3 pip_mom; TVector3 pim_mom; TVector3 ks_mom;
	TVector3 ks_vert = Kinematics::LambdaVertex(dMagneticField, pip_par, pim_par, pip_theta_min, pip_theta_max, pim_theta_min, pim_theta_max, pip_mom, pim_mom, ks_mom, pipi_dist);
	if(TMath::IsNaN(pipi_dist)) continue;
	ks_mom = pim_mom + pip_mom;

	TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
	TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
	TLorentzVector Lks = Lpim + Lpip;

	if(TMath::Abs(ks_vert.x()) > 250. || TMath::Abs(ks_vert.z()) > 250. || TMath::Abs(ks_vert.y()) > 250.) continue; //Vertex cut
	if(pipi_dist > pipi_distcut) continue;
	Double_t pim_vertex_dist; Double_t pip_vertex_dist;
	if(!Kinematics::HelixDirection(ks_vert, pip_start, pip_end, pip_vertex_dist) ||
	   !Kinematics::HelixDirection(ks_vert, pim_start, pim_end, pim_vertex_dist)) continue;
	
	if(pim_vertex_dist > pim_vtx_distcut) continue;
 	if(pip_vertex_dist > pip_vtx_distcut) continue;
	if(TMath::Abs(Lks.M() - KSMass) > ks_masscut) continue;

	Double_t kstarget_dist;
	TVector3 kstarget_vtx = Kinematics::CalcCloseDistLambda(tgtpos,
							       ks_vert,
							       ks_mom,
							       kstarget_dist);


	Ks_pip_id_container.push_back(it1);
	Ks_pim_id_container.push_back(it2);
	Ks_pip_purity_container.push_back(event.purity[it1]);
  Ks_pim_purity_container.push_back(event.purity[it2]);
  Ks_pip_efficiency_container.push_back(event.efficiency[it1]);
  Ks_pim_efficiency_container.push_back(event.efficiency[it2]);
  Ks_pip_G4tid_container.push_back(G4TrackID[it1]);
  Ks_pim_G4tid_container.push_back(G4TrackID[it2]);
  Ks_mass_container.push_back(Lks.M());
	Ks_mom_container.push_back(ks_mom);
	Ks_pip_mom_container.push_back(pip_mom);
	Ks_pim_mom_container.push_back(pim_mom);
	Ks_pipidist_container.push_back(pipi_dist);
	Ks_vtx_container.push_back(ks_vert);
	Ks_pip_repid_container.push_back(repid_pip);
	Ks_pim_repid_container.push_back(repid_pim);
	Ks_targetdist_container.push_back(kstarget_dist);
	Ks_targetvtx_container.push_back(kstarget_vtx);

	ks_candidates++;

  if(TMath::Abs(Lks.M() - KSMass) < mass_diff){
    best_ks_id = ks_candidates - 1;
    mass_diff = TMath::Abs(Lks.M() - KSMass);
  }
	for(int it3=0;it3<ntTpc;it3++){
    #if DebugDisp
      cout<<Form("Processing track (%d,%d,%d) / %d",it1,it2,it3,ntTpc)<<endl;
    #endif
	  if(it3==it2 or it3==it1) continue;
	  if((event.pid[it3]&4)!=4) continue; //select p like
	  if(event.charge[it3]!=1) continue;

	  Int_t repid_p = 0;
    Int_t flag = 1;
    for(Int_t i=0;i<2;i++){
      Int_t temp = flag&event.pid[it3];
      if(temp==flag) repid_p += 1;
      flag*=2;
    }
    cout<<Form("Checking Track %d, repid = %d",it3,repid_p)<<endl;
	  if(!GFTrackCont.TrackCheck(it3, repid_p)) continue;
    cout<<Form("Track %d, repid = %d checked",it3,repid_p)<<endl;

	  Double_t p_par[5];
	  p_par[0] = event.helix_cx[it3];
	  p_par[1] = event.helix_cy[it3];
	  p_par[2] = event.helix_z0[it3];
	  p_par[3] = event.helix_r[it3];
	  p_par[4] = event.helix_dz[it3];

	  Int_t p_nh = event.helix_t[it3].size();
      
    Double_t p_theta_min = event.helix_t[it3][0] - vtx_scan_range/p_par[3];
    Double_t p_theta_max = TMath::Min(event.helix_t[it3][0] + vtx_scan_range/p_par[3], event.helix_t[it3][p_nh-1]);

	  TVector3 p_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
	  TVector3 p_end = TVector3(event.calpos_x[it3][p_nh-1], event.calpos_y[it3][p_nh-1], event.calpos_z[it3][p_nh-1]);

	  TVector3 p_mom; Double_t ksp_dist;
    cout<<"Searching KsP Vertex"<<endl; 
	  TVector3 ksp_vert = Kinematics::XiVertex(dMagneticField,
						  p_par,
						  p_theta_min,
						  p_theta_max,
						  ks_vert,
              ks_mom,
              p_mom,
						  ksp_dist);
	  if(TMath::IsNaN(ksp_dist)) continue;
    p_mom = -p_mom;//p is positive charge, while XiVertex assumes Lpi-.
	  TLorentzVector Lp(p_mom, TMath::Sqrt(p_mom.Mag()*p_mom.Mag() + ProtonMass*ProtonMass));
	  TLorentzVector Lks_fixedmass(ks_mom, TMath::Sqrt(ks_mom.Mag()*ks_mom.Mag() + KSMass*KSMass));
    TLorentzVector Lksp = Lks_fixedmass + Lp;
	  TVector3 ksp_mom = TVector3(Lksp.Px(), Lksp.Py(), Lksp.Pz());
	  Double_t p_vertex_dist;
	  if(ksp_dist > ksp_distcut) continue;
	  if(!Kinematics::HelixDirection(ksp_vert, p_start, p_end, p_vertex_dist)) continue;
	  if(p_vertex_dist > p_vtx_distcut) continue;




	  ksp_ks_container.push_back(ks_candidates - 1);
	  ksp_pip_container.push_back(it1);
	  ksp_pim_container.push_back(it2);
	  ksp_p_container.push_back(it3);

    ksp_pip_purity_container.push_back(event.purity[it1]);
    ksp_pim_purity_container.push_back(event.purity[it2]);
    ksp_p_purity_container.push_back(event.purity[it3]);
    ksp_pip_efficiency_container.push_back(event.efficiency[it1]);
    ksp_pim_efficiency_container.push_back(event.efficiency[it2]);
    ksp_p_efficiency_container.push_back(event.efficiency[it3]);
    ksp_pip_G4tid_container.push_back(G4TrackID[it1]);
    ksp_pim_G4tid_container.push_back(G4TrackID[it2]);
    ksp_p_G4tid_container.push_back(G4TrackID[it3]);


	  ksp_pip_repid_container.push_back(repid_pip);
	  ksp_pim_repid_container.push_back(repid_pim);
	  ksp_p_repid_container.push_back(repid_p);

	  ksp_mass_container.push_back(Lksp.M());
	  Ks_mass_container.push_back(Lks.M());
	  ksp_mom_container.push_back(ksp_mom);
	  Ks_mom_container.push_back(ks_mom);
	  ksp_pip_mom_container.push_back(pip_mom);
	  ksp_pim_mom_container.push_back(pim_mom);
	  ksp_p_mom_container.push_back(p_mom);
//	  ksp_decayvertex_container.push_back(ksp_vert);
//	  Ks_vert_container.push_back(ks_vert);

	  ksp_candidates++;
	} //it3
      } //it2
    } //it1
  }
#if DebugDisp
  std::cout<<Form("Ks, KsP candidates searching finished: nks, nksP = (%d,%d) ",ks_candidates,ksp_candidates)<<std::endl;
#endif

  std::vector<Double_t> GFKs_mass_container(ks_candidates, qnan);
  std::vector<Double_t> GFKs_pipidist_container(ks_candidates, qnan);
  std::vector<Double_t> GFKs_targetdist_container(ks_candidates, qnan);
  std::vector<TVector3> GFKs_targetvtx_container(ks_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFKs_targetcenterdist_container(ks_candidates, qnan);
  std::vector<TVector3> GFKs_targetcentervtx_container(ks_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFKs_mom_container(ks_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFKs_vtx_container(ks_candidates, TVector3(qnan, qnan, qnan));

  std::vector<Int_t> GFKs_pip_id_container(ks_candidates, -1);
  std::vector<Int_t> GFKs_pim_id_container(ks_candidates, -1);
  std::vector<Int_t> GFKs_pip_repid_container(ks_candidates, qnan);
  std::vector<Int_t> GFKs_pim_repid_container(ks_candidates, qnan);
  std::vector<TVector3> GFKs_pip_mom_container(ks_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFKs_pim_mom_container(ks_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFKs_pip_extrapolation_container(ks_candidates, qnan);
  std::vector<Double_t> GFKs_pim_extrapolation_container(ks_candidates, qnan);
  std::vector<Double_t> GFKs_pip_tracklen_container(ks_candidates, qnan);
  std::vector<Double_t> GFKs_pim_tracklen_container(ks_candidates, qnan);
  if(ks_candidates>0){
#if DebugDisp
  cout<<"Sorting KS Candidates"<<endl;
#endif
    //Reconstructed Ks
    for(int idp=0;idp<ks_candidates;idp++){

      #if DebugDisp
      cout<<Form("Processing Ks candidate %d / %d",idp,ks_candidates)<<endl;
      #endif
      Int_t pip_id = Ks_pip_id_container[idp];
      Int_t pim_id = Ks_pim_id_container[idp];
      Int_t pip_repid = Ks_pip_repid_container[idp];
      Int_t pim_repid = Ks_pim_repid_container[idp];
      Double_t pip_extrapolation; Double_t pim_extrapolation;
      TVector3 pip_mom; TVector3 pim_mom;
      Double_t pipi_dist; TVector3 ks_vertex;

      Bool_t vtxcut = (GFTrackCont.FindVertex(pip_id, pim_id,
					      pip_repid, pim_repid,
					      pip_extrapolation, pim_extrapolation,
					      pip_mom, pim_mom,
					      pipi_dist, ks_vertex,
					      vtx_scan_range)
		       and pipi_dist < GFpipi_distcut and !TMath::IsNaN(pipi_dist) );
      if(!vtxcut) continue;

      TVector3 ks_mom = pip_mom + pim_mom;



      TLorentzVector GFLpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
      TLorentzVector GFLpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
      TLorentzVector GFLks = GFLpip + GFLpim;

      GFKs_pip_id_container[idp] = pip_id;
      GFKs_pim_id_container[idp] = pim_id;
      GFKs_pip_repid_container[idp] = pip_repid;
      GFKs_pim_repid_container[idp] = pim_repid;
      GFKs_pip_mom_container[idp] = pip_mom;
      GFKs_pim_mom_container[idp] = pim_mom;
      GFKs_pip_extrapolation_container[idp] = pip_extrapolation;
      GFKs_pim_extrapolation_container[idp] = pim_extrapolation;
      GFKs_mass_container[idp] = GFLks.M();
      GFKs_mom_container[idp] = ks_mom;
      GFKs_pipidist_container[idp] = pipi_dist;
      GFKs_vtx_container[idp] = ks_vertex;
    }
  } //ks_candidates
  int pip_id=-1,pim_id=-1;
  if(ks_candidates>0 ){
    pip_id = GFKs_pip_id_container[best_ks_id];
    pim_id = GFKs_pim_id_container[best_ks_id];
  }
  if(ks_candidates>0 and pip_id>-1 and pim_id>-1){
    #if DebugDisp
    cout<<"Selecting Best Ks"<<endl;
    #endif
#if DoKF
#if DebugDisp
    cout<<"Doing Kinematic Fit"<<endl;
#endif
      auto trpip = TPCAna.GetTrackTPCHelix(pip_id);
      auto trpim = TPCAna.GetTrackTPCHelix(pim_id);
      #if DebugDisp
      cout<<"Helix loaded"<<endl;
      #endif
      auto Vpip = trpip->GetCovarianceMatrix(0);
      auto Vpim = trpim->GetCovarianceMatrix(0);
      #if DebugDisp
      cout<<"Covariance loaded"<<endl;
      #endif
      double Diag_pipi[6]= {
        Vpip(0,0),Vpip(1,1),Vpip(2,2),Vpim(0,0),Vpim(1,1),Vpim(2,2)
      };
      auto Offdiag_pipi = MathTools::MergeOffdiagonals(Vpip,Vpim);
      auto pip_mom = GFKs_pip_mom_container[best_ks_id];
      auto pim_mom = GFKs_pim_mom_container[best_ks_id];
      auto ks_mom = GFKs_mom_container[best_ks_id];
      TVector3 HTVPip(pip_mom.x(),pip_mom.z(),pip_mom.y());
      TVector3 HTVPim(pim_mom.x(),pim_mom.z(),pim_mom.y());
      TVector3 HTVKs(ks_mom.x(),ks_mom.z(),ks_mom.y());
      TLorentzVector HLVPip(HTVPip,TMath::Hypot(pip_mom.Mag(),PionMass));
      TLorentzVector HLVPim(HTVPim,TMath::Hypot(pim_mom.Mag(),PionMass));
      TLorentzVector HLVKs(HTVKs,TMath::Hypot(ks_mom.Mag(),KSMass));
      FourVectorFitter KsFitter(HLVPip,HLVPim,HLVKs);
      KsFitter.SetInvMass(KSMass);
      KsFitter.SetVariance(Diag_pipi);
      KsFitter.AddOffdiagonals(Offdiag_pipi);
      double KFchi2Ks = KsFitter.DoKinematicFit();
      double KFpvalKs = KsFitter.GetPValue();
      auto KFpullKs = KsFitter.GetPull();
      #if DebugDisp
      cout<<"KFchi2Ks = "<<KFchi2Ks<<endl;
      cout<<"KFpvalKs = "<<KFpvalKs<<endl;
      cout<<"KFPull : ";
      for(auto p:KFpullKs) cout<<p<<" ";
      #endif 
      auto HLVcontKs = KsFitter.GetFittedLV();
      auto KFHLVPip = HLVcontKs[0];
      auto KFHLVPim = HLVcontKs[1];
      auto KFHLVKs = HLVcontKs[2];
      auto KFLVPip = TLorentzVector(KFHLVPip.Px(),KFHLVPip.Pz(),KFHLVPip.Py(),KFHLVPip.E());
      auto KFLVPim = TLorentzVector(KFHLVPim.Px(),KFHLVPim.Pz(),KFHLVPim.Py(),KFHLVPim.E());
      auto KFLVKs = TLorentzVector(KFHLVKs.Px(),KFHLVKs.Pz(),KFHLVKs.Py(),KFHLVKs.E());  
      auto KFpipmom = TVector3(KFLVPip.Px(),KFLVPip.Py(),KFLVPip.Pz());
      auto KFpimmom = TVector3(KFLVPim.Px(),KFLVPim.Py(),KFLVPim.Pz());
      auto KFKsmom = TVector3(KFLVKs.Px(),KFLVKs.Py(),KFLVKs.Pz());

#endif

      event.Ksflag = 1;
      event.Ksmass = GFKs_mass_container[best_ks_id];
      event.Ksdecayvtx_x = GFKs_vtx_container[best_ks_id].X();
      event.Ksdecayvtx_y = GFKs_vtx_container[best_ks_id].Y();
      event.Ksdecayvtx_z = GFKs_vtx_container[best_ks_id].Z();
      event.Ksmom = GFKs_mom_container[best_ks_id].Mag();
      event.Ksmom_x = GFKs_mom_container[best_ks_id].X();
      event.Ksmom_y = GFKs_mom_container[best_ks_id].Y();
      event.Ksmom_z = GFKs_mom_container[best_ks_id].Z();
      event.Ksdist = GFKs_pipidist_container[best_ks_id];

      event.Kspipid = GFKs_pip_id_container[best_ks_id];
      event.GFpimmom = GFKs_pim_mom_container[best_ks_id].Mag();
      event.GFpimmom_x = GFKs_pim_mom_container[best_ks_id].X();
      event.GFpimmom_y = GFKs_pim_mom_container[best_ks_id].Y();
      event.GFpimmom_z = GFKs_pim_mom_container[best_ks_id].Z();
      event.pimmom = Ks_pim_mom_container[best_ks_id].Mag();
      event.pimmom_x = Ks_pim_mom_container[best_ks_id].X();
      event.pimmom_y = Ks_pim_mom_container[best_ks_id].Y();
      event.pimmom_z = Ks_pim_mom_container[best_ks_id].Z();

      event.Kspimid = GFKs_pim_id_container[best_ks_id];
      event.GFpipmom = GFKs_pip_mom_container[best_ks_id].Mag();
      event.GFpipmom_x = GFKs_pip_mom_container[best_ks_id].X();
      event.GFpipmom_y = GFKs_pip_mom_container[best_ks_id].Y();
      event.GFpipmom_z = GFKs_pip_mom_container[best_ks_id].Z();
      event.pipmom = Ks_pip_mom_container[best_ks_id].Mag();
      event.pipmom_x = Ks_pip_mom_container[best_ks_id].X();
      event.pipmom_y = Ks_pip_mom_container[best_ks_id].Y();
      event.pipmom_z = Ks_pip_mom_container[best_ks_id].Z();
#if DoKF
      event.KFpvalKs = KFpvalKs;
      event.KFpullKs = KFpullKs;
      event.KFKsmom = KFKsmom.Mag();
      event.KFKsmom_x = KFKsmom.X();
      event.KFKsmom_y = KFKsmom.Y();
      event.KFKsmom_z = KFKsmom.Z();
      event.KFpimmom = KFpimmom.Mag();
      event.KFpimmom_x = KFpimmom.X();
      event.KFpimmom_y = KFpimmom.Y();
      event.KFpimmom_z = KFpimmom.Z();
      event.KFpipmom = KFpipmom.Mag();
      event.KFpipmom_x = KFpipmom.X();
      event.KFpipmom_y = KFpipmom.Y();
      event.KFpipmom_z = KFpipmom.Z();
#endif
  }


  std::vector<Int_t> GFksp_pip_id_container(ksp_candidates, -1);
  std::vector<Int_t> GFksp_pim_id_container(ksp_candidates, -1);
  std::vector<Int_t> GFksp_p_id_container(ksp_candidates, -1);
  std::vector<Int_t> GFksp_pip_rep_container(ksp_candidates, -1);
  std::vector<Int_t> GFksp_pim_rep_container(ksp_candidates, -1);
  std::vector<Int_t> GFksp_p_rep_container(ksp_candidates, -1);
  std::vector<Int_t> GFksp_id_container(ksp_candidates, -1);
  std::vector<TVector3> GFksp_mom_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFksp_decayvertex_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFksp_mass_container(ksp_candidates, qnan);
  std::vector<TVector3> GFksp_pos_targetcenter_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFksp_mom_targetcenter_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFksp_dist_targetcenter_container(ksp_candidates, qnan);
  std::vector<TVector3> GFks_mom_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFks_vert_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFks_mass_container(ksp_candidates, qnan);
  std::vector<Double_t> GFks_tracklen_container(ksp_candidates, qnan);
  std::vector<Double_t> GFks_tof_container(ksp_candidates, qnan);
  std::vector<TVector3> GFpip_mom_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFpim_mom_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<TVector3> GFp_mom_container(ksp_candidates, TVector3(qnan, qnan, qnan));
  std::vector<Double_t> GFpipi_closedist_container(ksp_candidates, qnan);
  std::vector<Double_t> GFksp_closedist_container(ksp_candidates, qnan);

  for(Int_t candi=0;candi<ksp_candidates;candi++){
    #if DebugDisp
    cout<<Form("Processing KsP candidate %d / %d",candi,ksp_candidates)<<endl;
    #endif
    Int_t trackid_pip = ksp_pip_container[candi];
    Int_t trackid_pim = ksp_pim_container[candi];
    Int_t trackid_p = ksp_p_container[candi];
    Int_t repid_pip = ksp_pip_repid_container[candi];
    Int_t repid_pim = ksp_pim_repid_container[candi];
    Int_t repid_p = ksp_p_repid_container[candi];
    Int_t ks_id = ksp_ks_container[candi];
    if(TMath::IsNaN(GFKs_pipidist_container[ks_id])) continue; //Genfit's fitting was succeeded.
    double GFpipi_dist = GFKs_pipidist_container[ks_id];
    TVector3 GFks_vert = GFKs_vtx_container[ks_id];
    Double_t GFextrapolation_decays[3];
    GFextrapolation_decays[0] = GFKs_pip_extrapolation_container[ks_id];
    GFextrapolation_decays[1] = GFKs_pim_extrapolation_container[ks_id];
    TVector3 GFmom_decays[3];
    GFmom_decays[0] = GFKs_pip_mom_container[ks_id];
    GFmom_decays[1] = GFKs_pim_mom_container[ks_id];
    TLorentzVector GFKSpip(GFmom_decays[0], TMath::Hypot(GFmom_decays[0].Mag(), PionMass));
    TLorentzVector GFKSpim(GFmom_decays[1], TMath::Hypot(GFmom_decays[1].Mag(), PionMass));
    TLorentzVector GFKSks = GFKSpip + GFKSpim;
    TVector3 GFks_mom = GFmom_decays[0] + GFmom_decays[1];
#if DoKF
    GFks_mom = TVector3(event.KFKsmom_x,event.KFKsmom_y,event.KFKsmom_z);
#endif
    TLorentzVector GFKSks_fixed(GFks_mom, hypot(GFks_mom.Mag(), KSMass));
    TVector3 GFksp_vert; Double_t GFksp_dist = qnan; Double_t GFks_tracklen;
    if(!GFTrackCont.FindVertexXi(trackid_p, repid_p, GFks_vert, GFks_mom, GFks_tracklen,
				 GFextrapolation_decays[2], GFmom_decays[2], GFksp_dist, GFksp_vert, vtx_scan_range)
       || GFksp_dist > GFksp_distcut) continue;

    TLorentzVector GFKSp(GFmom_decays[2],hypot(GFmom_decays[2].Mag(),ProtonMass));
    TLorentzVector GFKSksp = GFKSks_fixed + GFKSp;
    TVector3 GFksp_mom = GFks_mom + GFmom_decays[2];



    GFks_mom = GFmom_decays[0] + GFmom_decays[1];
    GFksp_pip_id_container[candi] = trackid_pip;
    GFksp_pim_id_container[candi] = trackid_pim;
    GFksp_p_id_container[candi] = trackid_p;
    GFksp_pip_rep_container[candi] = repid_pip;
    GFksp_pim_rep_container[candi] = repid_pim;
    GFksp_p_rep_container[candi] = repid_p;
    GFksp_mom_container[candi] = GFksp_mom;
    GFksp_decayvertex_container[candi] = GFksp_vert;
    GFksp_mass_container[candi] = GFKSksp.M();

    GFks_mom_container[candi] = GFks_mom;
    GFks_vert_container[candi] = GFks_vert;
    GFks_mass_container[candi] = GFKSks.M();
    GFks_tracklen_container[candi] = GFks_tracklen;
    GFpip_mom_container[candi] = GFmom_decays[0];
    GFpim_mom_container[candi] = GFmom_decays[1];
    GFp_mom_container[candi] = GFmom_decays[2];
    GFpipi_closedist_container[candi] = GFpipi_dist;
    GFksp_closedist_container[candi] = GFksp_dist;
    //Double_t diff = TMath::Abs(lambda_mass_container[candi] - LambdaMass);
  } //kspcandi
  int ksp_cnt=0,ksp1_id=-1, ksp2_id=-1;
  for(int ic=0;ic<ksp_candidates;ic++){
    #if DebugDisp
    cout<<Form("Processing KsP candidate %d / %d",ic,ksp_candidates)<<endl;
    #endif
    if(ksp_ks_container[ic] != best_ks_id) continue;
    ksp_cnt++;
    if(ksp_cnt==1){
      event.KsP1flag = 1;
      event.KsP1mass = GFksp_mass_container[ic];
      event.KsP1pid = ksp_p_container[ic];
      event.KsP1decayvtx_x = GFksp_decayvertex_container[ic].x();
      event.KsP1decayvtx_y = GFksp_decayvertex_container[ic].y();
      event.KsP1decayvtx_z = GFksp_decayvertex_container[ic].z();
      event.KsP1mom = GFksp_mom_container[ic].Mag();
      event.KsP1mom_x = GFksp_mom_container[ic].x();
      event.KsP1mom_y = GFksp_mom_container[ic].y();
      event.KsP1mom_z = GFksp_mom_container[ic].z();
      event.KsP1dist = GFksp_closedist_container[ic];
      event.GFp1mom = GFp_mom_container[ic].Mag();
      event.GFp1mom_x = GFp_mom_container[ic].x();
      event.GFp1mom_y = GFp_mom_container[ic].y();
      event.GFp1mom_z = GFp_mom_container[ic].z();
      event.p1mom = ksp_p_mom_container[ic].Mag();
      event.p1mom_x = ksp_p_mom_container[ic].x();
      event.p1mom_y = ksp_p_mom_container[ic].y();
      event.p1mom_z = ksp_p_mom_container[ic].z();
    }
    if(ksp_cnt==2){
      event.KsP2flag = 1;
      event.KsP2mass = GFksp_mass_container[ic];
      event.KsP2pid = ksp_p_container[ic];
      event.KsP2decayvtx_x = GFksp_decayvertex_container[ic].x();
      event.KsP2decayvtx_y = GFksp_decayvertex_container[ic].y();
      event.KsP2decayvtx_z = GFksp_decayvertex_container[ic].z();
      event.KsP2mom = GFksp_mom_container[ic].Mag();
      event.KsP2mom_x = GFksp_mom_container[ic].x();
      event.KsP2mom_y = GFksp_mom_container[ic].y();
      event.KsP2mom_z = GFksp_mom_container[ic].z();
      event.KsP2dist = GFksp_closedist_container[ic];
      event.GFp2mom = GFp_mom_container[ic].Mag();
      event.GFp2mom_x = GFp_mom_container[ic].x();
      event.GFp2mom_y = GFp_mom_container[ic].y();
      event.GFp2mom_z = GFp_mom_container[ic].z();
      event.p2mom = ksp_p_mom_container[ic].Mag();
      event.p2mom_x = ksp_p_mom_container[ic].x();
      event.p2mom_y = ksp_p_mom_container[ic].y();
      event.p2mom_z = ksp_p_mom_container[ic].z();
    }
  }
  event.KsPcnt = ksp_cnt;
  Int_t GFtrackid_decays[2] = {-1, -1}; Int_t GFrepid_decays[2];

#if DebugDisp
  cout<<"Searching Forced Ks"<<endl;
#endif
  //Force Recon//
  if(event.G4piptid != -1 and event.G4pimtid != -1){
      bool ForcedKSflag = true;
      int it1 = event.G4piptid;
      int it2 = event.G4pimtid;
      #if DebugDisp
      cout<<Form("Processing Forced Ks: %d, %d",it1,it2)<<endl;
      #endif
      if((event.pid[it1]&1)!=1) ForcedKSflag = false; //select pi+like
      if(event.charge[it1]!=1) ForcedKSflag = false;
      if((event.pid[it2]&1)!=1) ForcedKSflag = false; //select pi-like
      if(event.charge[it2]!=-1) ForcedKSflag = false;
      if(!GFTrackCont.TrackCheck(it1,0)) ForcedKSflag = false;
      if(!GFTrackCont.TrackCheck(it2,0)) ForcedKSflag = false;
      if(ForcedKSflag){
        #if DebugDisp
        cout<<Form("Forced Ks: %d, %d",it1,it2)<<endl;
        #endif
        Double_t pip_par[5];
        pip_par[0] = event.helix_cx[it1];
        pip_par[1] = event.helix_cy[it1];
        pip_par[2] = event.helix_z0[it1];
        pip_par[3] = event.helix_r[it1];
        pip_par[4] = event.helix_dz[it1];
        Int_t pip_nh = event.helix_t[it1].size();
        Double_t pip_theta_min = event.helix_t[it1][0] - vtx_scan_range/pip_par[3];
        Double_t pip_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_range/pip_par[3], event.helix_t[it1][pip_nh-1]);
        TVector3 pip_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
        TVector3 pip_end = TVector3(event.calpos_x[it1][pip_nh-1], event.calpos_y[it1][pip_nh-1], event.calpos_z[it1][pip_nh-1]);

        Double_t pim_par[5];
        pim_par[0] = event.helix_cx[it2];
        pim_par[1] = event.helix_cy[it2];
        pim_par[2] = event.helix_z0[it2];
        pim_par[3] = event.helix_r[it2];
        pim_par[4] = event.helix_dz[it2];
        Int_t pim_nh = event.helix_t[it2].size();
        Double_t pim_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/pim_par[3], event.helix_t[it2][pim_nh-1]);
        Double_t pim_theta_max = event.helix_t[it2][0] + vtx_scan_range/pim_par[3];
        TVector3 pim_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
        TVector3 pim_end = TVector3(event.calpos_x[it2][pim_nh-1], event.calpos_y[it2][pim_nh-1], event.calpos_z[it2][pim_nh-1]);
    
    
        Double_t pipi_dist = 10000.;
        TVector3 pip_mom; TVector3 pim_mom; TVector3 ks_mom;
        TVector3 ks_vert = Kinematics::LambdaVertex(dMagneticField, pip_par, pim_par, pip_theta_min, pip_theta_max, pim_theta_min, pim_theta_max, pip_mom, pim_mom, ks_mom, pipi_dist);
        if(TMath::IsNaN(pipi_dist)) ForcedKSflag = false;
        ks_mom = pim_mom + pip_mom;

        Double_t pim_vertex_dist; Double_t pip_vertex_dist;
        Kinematics::HelixDirection(ks_vert, pip_start, pip_end, pip_vertex_dist);
        Kinematics::HelixDirection(ks_vert, pim_start, pim_end, pim_vertex_dist);
        TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
        TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
        TLorentzVector Lks = Lpim + Lpip;

        int pip_repid = 0, pim_repid = 0;
        Double_t pip_extrapolation; Double_t pim_extrapolation;
        TVector3 ks_vertex;
        Bool_t vtxcut = (GFTrackCont.FindVertex(it1, it2,
                  pip_repid, pim_repid,
                  pip_extrapolation, pim_extrapolation,
                  pip_mom, pim_mom,
                  pipi_dist, ks_vertex,
                  vtx_scan_range)
            and pipi_dist < GFpipi_distcut and !TMath::IsNaN(pipi_dist) );
        if(!vtxcut){
           event.GFForcedKsflag = false;
        }
        else{
          event.GFForcedKsflag = true;
        }
      if(event.GFForcedKsflag){
#if DoKF
        auto Vpip = TPCAna.GetTrackTPCHelix(it1)->GetCovarianceMatrix();
        auto Vpim = TPCAna.GetTrackTPCHelix(it2)->GetCovarianceMatrix();
        double Diag_pipi[6]= {
          Vpip(0,0),Vpip(1,1),Vpip(2,2),Vpim(0,0),Vpim(1,1),Vpim(2,2)
        };
        auto Offdiag_pipi = MathTools::MergeOffdiagonals(Vpip,Vpim);
        auto ks_mom = pip_mom + pim_mom;
        TVector3 HTVPip(pip_mom.x(),pip_mom.z(),pip_mom.y());
        TVector3 HTVPim(pim_mom.x(),pim_mom.z(),pim_mom.y());
        TVector3 HTVKs(ks_mom.x(),ks_mom.z(),ks_mom.y());
        TLorentzVector HLVPip(HTVPip,TMath::Hypot(pip_mom.Mag(),PionMass));
        TLorentzVector HLVPim(HTVPim,TMath::Hypot(pim_mom.Mag(),PionMass));
        TLorentzVector HLVKs(HTVKs,TMath::Hypot(ks_mom.Mag(),KSMass));
        FourVectorFitter KsFitter(HLVPip,HLVPim,HLVKs);
        KsFitter.SetInvMass(KSMass);
        KsFitter.SetVariance(Diag_pipi);
        KsFitter.AddOffdiagonals(Offdiag_pipi);
        double KFchi2Ks = KsFitter.DoKinematicFit();
        double KFpvalKs = KsFitter.GetPValue();
        auto KFpullKs = KsFitter.GetPull();
        auto HLVcontKs = KsFitter.GetFittedLV();
        auto KFHLVPip = HLVcontKs[0];
        auto KFHLVPim = HLVcontKs[1];
        auto KFHLVKs = HLVcontKs[2];
        auto KFLVPip = TLorentzVector(KFHLVPip.Px(),KFHLVPip.Pz(),KFHLVPip.Py(),KFHLVPip.E());
        auto KFLVPim = TLorentzVector(KFHLVPim.Px(),KFHLVPim.Pz(),KFHLVPim.Py(),KFHLVPim.E());
        auto KFLVKs = TLorentzVector(KFHLVKs.Px(),KFHLVKs.Pz(),KFHLVKs.Py(),KFHLVKs.E());  
        auto KFpipmom = TVector3(KFLVPip.Px(),KFLVPip.Py(),KFLVPip.Pz());
        auto KFpimmom = TVector3(KFLVPim.Px(),KFLVPim.Py(),KFLVPim.Pz());
        auto KFKsmom = TVector3(KFLVKs.Px(),KFLVKs.Py(),KFLVKs.Pz());
        event.ForcedKFpvalKs = KFpvalKs;
        event.ForcedKFpullKs = KFpullKs;
        event.ForcedKFKsmom = KFKsmom.Mag();
        event.ForcedKFKsmom_x = KFKsmom.X();
        event.ForcedKFKsmom_y = KFKsmom.Y();
        event.ForcedKFKsmom_z = KFKsmom.Z();
        event.ForcedKFpipmom = KFpipmom.Mag();
        event.ForcedKFpipmom_x = KFpipmom.X();
        event.ForcedKFpipmom_y = KFpipmom.Y();
        event.ForcedKFpipmom_z = KFpipmom.Z();
        event.ForcedKFpimmom = KFpimmom.Mag();
        event.ForcedKFpimmom_x = KFpimmom.X();
        event.ForcedKFpimmom_y = KFpimmom.Y();
        event.ForcedKFpimmom_z = KFpimmom.Z();
#endif
      }
        event.ForcedKsflag = ForcedKSflag;
        event.ForcedKsmass = Lks.M();
        event.ForcedKsdecayvtx_x = ks_vert.x();
        event.ForcedKsdecayvtx_y = ks_vert.y();
        event.ForcedKsdecayvtx_z = ks_vert.z();
        event.ForcedKsmom = ks_mom.Mag();
        event.ForcedKsmom_x = ks_mom.x();
        event.ForcedKsmom_y = ks_mom.y();
        event.ForcedKsmom_z = ks_mom.z();
        event.ForcedKsdist = pipi_dist;
        event.Forcedpipvertexdist = pip_vertex_dist;
        event.Forcedpimvertexdist = pim_vertex_dist;
        event.Forcedpipmom = pip_mom.Mag();
        event.Forcedpipmom_x = pip_mom.x();
        event.Forcedpipmom_y = pip_mom.y();
        event.Forcedpipmom_z = pip_mom.z();
        event.Forcedpimmom = pim_mom.Mag();
        event.Forcedpimmom_x = pim_mom.x();
        event.Forcedpimmom_y = pim_mom.y();
        event.Forcedpimmom_z = pim_mom.z();
    

        if(event.G4p1tid != -1){
          int it3 = event.G4p1tid;
          bool ForcedKsP1flag = true;
          #if DebugDisp
          cout<<Form("Processing Forced KsP1: %d",it3)<<endl;
          #endif
          if((event.pid[it3]&4)!=4) ForcedKsP1flag = false; //select p like
          if(event.charge[it3]!=1) ForcedKsP1flag = false;
          if(ForcedKsP1flag){
            Int_t repid_p = 0;
            Int_t flag = 1;
            for(Int_t i=0;i<2;i++){
              Int_t temp = flag&event.pid[it3];
              if(temp==flag) repid_p += 1;
              flag*=2;
            }
            #if DebugDisp
            cout<<Form("Checking Forced Track %d, repid = %d",it3,repid_p)<<endl;
            #endif
            if(!GFTrackCont.TrackCheck(it3, repid_p)) ForcedKsP1flag = false;
            #if DebugDisp
            cout<<Form("Track %d, repid = %d checked",it3,repid_p)<<endl;
            #endif
          }
          if(ForcedKsP1flag){
            #if DebugDisp
            cout<<Form("Forced KsP1: %d",it3)<<endl;
            #endif
            Double_t p_par[5];
            p_par[0] = event.helix_cx[it3];
            p_par[1] = event.helix_cy[it3];
            p_par[2] = event.helix_z0[it3];
            p_par[3] = event.helix_r[it3];
            p_par[4] = event.helix_dz[it3];
            Int_t p_nh = event.helix_t[it3].size();
              
            Double_t p_theta_min = event.helix_t[it3][0] - vtx_scan_range/p_par[3];
            Double_t p_theta_max = TMath::Min(event.helix_t[it3][0] + vtx_scan_range/p_par[3], event.helix_t[it3][p_nh-1]);

            TVector3 p_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
            TVector3 p_end = TVector3(event.calpos_x[it3][p_nh-1], event.calpos_y[it3][p_nh-1], event.calpos_z[it3][p_nh-1]);

            TVector3 p_mom; Double_t ksp_dist;
            #if DoKF
              ks_mom = TVector3(event.ForcedKFKsmom_x,event.ForcedKFKsmom_y,event.ForcedKFKsmom_z);
            #endif
            TVector3 ksp_vert = Kinematics::XiVertex(dMagneticField,
                      p_par,
                      p_theta_min,
                      p_theta_max,
                      ks_vert,
                      ks_mom,
                      p_mom,
                      ksp_dist);
            if(TMath::IsNaN(ksp_dist))ForcedKsP1flag = false;
            p_mom = -p_mom;//p is positive charge, while XiVertex assumes Lpi-.
            TLorentzVector Lp(p_mom, TMath::Sqrt(p_mom.Mag()*p_mom.Mag() + ProtonMass*ProtonMass));
            TLorentzVector Lks_fixedmass(ks_mom, TMath::Sqrt(ks_mom.Mag()*ks_mom.Mag() + KSMass*KSMass));
            TLorentzVector Lksp = Lks_fixedmass + Lp;
            TVector3 ksp_mom = TVector3(Lksp.Px(), Lksp.Py(), Lksp.Pz());
            Double_t p_vertex_dist;
            Kinematics::HelixDirection(ksp_vert, p_start, p_end, p_vertex_dist);
            event.ForcedKsP1flag = ForcedKsP1flag;
            event.ForcedKsP1mass = Lksp.M();
            event.ForcedKsP1decayvtx_x = ksp_vert.x();
            event.ForcedKsP1decayvtx_y = ksp_vert.y();
            event.ForcedKsP1decayvtx_z = ksp_vert.z();
            event.ForcedKsP1mom = ksp_mom.Mag();
            event.ForcedKsP1mom_x = ksp_mom.x();
            event.ForcedKsP1mom_y = ksp_mom.y();
            event.ForcedKsP1mom_z = ksp_mom.z();
            event.ForcedKsP1dist = ksp_dist;
            event.ForcedKsP1vertexdist = p_vertex_dist;
          }
        }
        if(event.G4p2tid != -1){
          int it3 = event.G4p2tid;
          bool ForcedKsP2flag = true;
          #if DebugDisp
          cout<<Form("Processing Forced KsP2: %d",it3)<<endl;
          #endif
          if((event.pid[it3]&4)!=4) ForcedKsP2flag = false; //select p like
          if(event.charge[it3]!=1) ForcedKsP2flag = false;
          if(ForcedKsP2flag){
            Int_t repid_p = 0;
            Int_t flag = 1;
            for(Int_t i=0;i<2;i++){
              Int_t temp = flag&event.pid[it3];
              if(temp==flag) repid_p += 1;
              flag*=2;
            }
            #if DebugDisp
            cout<<Form("Checking Forced Track %d, repid = %d",it3,repid_p)<<endl;
            #endif
            if(!GFTrackCont.TrackCheck(it3, repid_p)) ForcedKsP2flag = false;
            #if DebugDisp
            cout<<Form("Track %d, repid = %d checked",it3,repid_p)<<endl;
            #endif
          }
          if(ForcedKsP2flag){
            #if DebugDisp
            cout<<Form("Forced KsP2: %d",it3)<<endl;
            #endif
            Double_t p_par[5];
            p_par[0] = event.helix_cx[it3];
            p_par[1] = event.helix_cy[it3];
            p_par[2] = event.helix_z0[it3];
            p_par[3] = event.helix_r[it3];
            p_par[4] = event.helix_dz[it3];
            Int_t p_nh = event.helix_t[it3].size();
              
            Double_t p_theta_min = event.helix_t[it3][0] - vtx_scan_range/p_par[3];
            Double_t p_theta_max = TMath::Min(event.helix_t[it3][0] + vtx_scan_range/p_par[3], event.helix_t[it3][p_nh-1]);

            TVector3 p_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
            TVector3 p_end = TVector3(event.calpos_x[it3][p_nh-1], event.calpos_y[it3][p_nh-1], event.calpos_z[it3][p_nh-1]);

            TVector3 p_mom; Double_t ksp_dist;
            #if DoKF
              ks_mom = TVector3(event.ForcedKFKsmom_x,event.ForcedKFKsmom_y,event.ForcedKFKsmom_z);
            #endif
            TVector3 ksp_vert = Kinematics::XiVertex(dMagneticField,
                      p_par,
                      p_theta_min,
                      p_theta_max,
                      ks_vert,
                      ks_mom,
                      p_mom,
                      ksp_dist);
            if(TMath::IsNaN(ksp_dist))ForcedKsP2flag = false;
            p_mom = -p_mom;//p is positive charge, while XiVertex assumes Lpi-.
            TLorentzVector Lp(p_mom, TMath::Sqrt(p_mom.Mag()*p_mom.Mag() + ProtonMass*ProtonMass));
            TLorentzVector Lks_fixedmass(ks_mom, TMath::Sqrt(ks_mom.Mag()*ks_mom.Mag() + KSMass*KSMass));
            TLorentzVector Lksp = Lks_fixedmass + Lp;
            TVector3 ksp_mom = TVector3(Lksp.Px(), Lksp.Py(), Lksp.Pz());
            Double_t p_vertex_dist;
            Kinematics::HelixDirection(ksp_vert, p_start, p_end, p_vertex_dist);
            event.ForcedKsP2flag = ForcedKsP2flag;
            event.ForcedKsP2mass = Lksp.M();
            event.ForcedKsP2decayvtx_x = ksp_vert.x();
            event.ForcedKsP2decayvtx_y = ksp_vert.y();
            event.ForcedKsP2decayvtx_z = ksp_vert.z();
            event.ForcedKsP2mom = ksp_mom.Mag();
            event.ForcedKsP2mom_x = ksp_mom.x();
            event.ForcedKsP2mom_y = ksp_mom.y();
            event.ForcedKsP2mom_z = ksp_mom.z();
            event.ForcedKsP2dist = ksp_dist;
            event.ForcedKsP2vertexdist = p_vertex_dist;
          }
        }
      }
  }

  //Remaining tracks
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
	tree->Branch( "inc" , &event.inc );
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

  tree->Branch( "G4thetaid", &event.G4thetaid );
  tree->Branch( "G4thetamom", &event.G4thetamom );
  tree->Branch( "G4thetamom_x", &event.G4thetamom_x );
  tree->Branch( "G4thetamom_y", &event.G4thetamom_y );
  tree->Branch( "G4thetamom_z", &event.G4thetamom_z );
  tree->Branch( "G4thetavtx_x", &event.G4thetavtx_x );
  tree->Branch( "G4thetavtx_y", &event.G4thetavtx_y );
  tree->Branch( "G4thetavtx_z", &event.G4thetavtx_z );
  
  tree->Branch( "G4ksp1mass", &event.G4ksp1mass );
  tree->Branch( "G4ksp2mass", &event.G4ksp2mass );


  tree->Branch( "G4ksid", &event.G4ksid );
  tree->Branch( "G4ksmom", &event.G4ksmom );
  tree->Branch( "G4ksmom_x", &event.G4ksmom_x );
  tree->Branch( "G4ksmom_y", &event.G4ksmom_y );
  tree->Branch( "G4ksmom_z", &event.G4ksmom_z );
  tree->Branch( "G4ksvtx_x", &event.G4ksvtx_x );
  tree->Branch( "G4ksvtx_y", &event.G4ksvtx_y );
  tree->Branch( "G4ksvtx_z", &event.G4ksvtx_z );
  tree->Branch( "G4ksmass", &event.G4ksmass);

  tree->Branch( "G4pimid", &event.G4pimid );
  tree->Branch( "G4pimtid", &event.G4pimtid );
  tree->Branch( "G4pimnh", &event.G4pimnh );
  tree->Branch( "G4pimtnh", &event.G4pimtnh );
  tree->Branch( "G4pimvtx_x", &event.G4pimvtx_x );
  tree->Branch( "G4pimvtx_y", &event.G4pimvtx_y );
  tree->Branch( "G4pimvtx_z", &event.G4pimvtx_z );
  tree->Branch( "G4pimmom", &event.G4pimmom );
  tree->Branch( "G4pimmom_x", &event.G4pimmom_x );
  tree->Branch( "G4pimmom_y", &event.G4pimmom_y );
  tree->Branch( "G4pimmom_z", &event.G4pimmom_z );
  tree->Branch( "pimtid", &event.pimtid );
  tree->Branch( "pimnh", &event.pimnh );
  tree->Branch( "pimvtx_x", &event.pimvtx_x );
  tree->Branch( "pimvtx_y", &event.pimvtx_y );
  tree->Branch( "pimvtx_z", &event.pimvtx_z );
  tree->Branch( "pimmom", &event.pimmom );
  tree->Branch( "pimmom_x", &event.pimmom_x );
  tree->Branch( "pimmom_y", &event.pimmom_y );
  tree->Branch( "pimmom_z", &event.pimmom_z );
  tree->Branch( "GFpimmom", &event.GFpimmom );
  tree->Branch( "GFpimmom_x", &event.GFpimmom_x );
  tree->Branch( "GFpimmom_y", &event.GFpimmom_y );
  tree->Branch( "GFpimmom_z", &event.GFpimmom_z );
  tree->Branch( "KFpimmom", &event.KFpimmom );
  tree->Branch( "KFpimmom_x", &event.KFpimmom_x );
  tree->Branch( "KFpimmom_y", &event.KFpimmom_y );
  tree->Branch( "KFpimmom_z", &event.KFpimmom_z );

  tree->Branch( "G4pipid", &event.G4pipid );
  tree->Branch( "G4piptid", &event.G4piptid );
  tree->Branch( "G4pipnh", &event.G4pipnh );
  tree->Branch( "G4piptnh", &event.G4piptnh );
  tree->Branch( "G4pipvtx_x", &event.G4pipvtx_x );
  tree->Branch( "G4pipvtx_y", &event.G4pipvtx_y );
  tree->Branch( "G4pipvtx_z", &event.G4pipvtx_z );
  tree->Branch( "G4pipmom", &event.G4pipmom );
  tree->Branch( "G4pipmom_x", &event.G4pipmom_x );
  tree->Branch( "G4pipmom_y", &event.G4pipmom_y );
  tree->Branch( "G4pipmom_z", &event.G4pipmom_z );
  tree->Branch( "piptid", &event.piptid );
  tree->Branch( "pipmom", &event.pipmom );
  tree->Branch( "pipmom_x", &event.pipmom_x );
  tree->Branch( "pipmom_y", &event.pipmom_y );
  tree->Branch( "pipmom_z", &event.pipmom_z );
  tree->Branch( "GFpipmom", &event.GFpipmom );
  tree->Branch( "GFpipmom_x", &event.GFpipmom_x );
  tree->Branch( "GFpipmom_y", &event.GFpipmom_y );
  tree->Branch( "GFpipmom_z", &event.GFpipmom_z );
  tree->Branch( "KFpipmom", &event.KFpipmom );
  tree->Branch( "KFpipmom_x", &event.KFpipmom_x );
  tree->Branch( "KFpipmom_y", &event.KFpipmom_y );
  tree->Branch( "KFpipmom_z", &event.KFpipmom_z );

  tree->Branch("G4p1id", &event.G4p1id);
  tree->Branch("G4p1tid", &event.G4p1tid);
  tree->Branch("G4p1nh", &event.G4p1nh);
  tree->Branch("G4p1tnh", &event.G4p1tnh);
  tree->Branch("G4p1vtx_x", &event.G4p1vtx_x);
  tree->Branch("G4p1vtx_y", &event.G4p1vtx_y);
  tree->Branch("G4p1vtx_z", &event.G4p1vtx_z);
  tree->Branch("G4p1mom", &event.G4p1mom);
  tree->Branch("G4p1mom_x", &event.G4p1mom_x);
  tree->Branch("G4p1mom_y", &event.G4p1mom_y);
  tree->Branch("G4p1mom_z", &event.G4p1mom_z);
  tree->Branch("p1tid", &event.p1tid);
  tree->Branch("p1mom", &event.p1mom);
  tree->Branch("p1mom_x", &event.p1mom_x);
  tree->Branch("p1mom_y", &event.p1mom_y);
  tree->Branch("p1mom_z", &event.p1mom_z);
  tree->Branch("GFp1mom", &event.GFp1mom);
  tree->Branch("GFp1mom_x", &event.GFp1mom_x);
  tree->Branch("GFp1mom_y", &event.GFp1mom_y);
  tree->Branch("GFp1mom_z", &event.GFp1mom_z);

  tree->Branch("G4p2id", &event.G4p2id);
  tree->Branch("G4p2tid", &event.G4p2tid);
  tree->Branch("G4p2nh", &event.G4p2nh);
  tree->Branch("G4p2tnh", &event.G4p2tnh);
  tree->Branch("G4p2vtx_x", &event.G4p2vtx_x);
  tree->Branch("G4p2vtx_y", &event.G4p2vtx_y);
  tree->Branch("G4p2vtx_z", &event.G4p2vtx_z);
  tree->Branch("G4p2mom", &event.G4p2mom);
  tree->Branch("G4p2mom_x", &event.G4p2mom_x);
  tree->Branch("G4p2mom_y", &event.G4p2mom_y);
  tree->Branch("G4p2mom_z", &event.G4p2mom_z);
  tree->Branch("p2tid", &event.p2tid);
  tree->Branch("p2mom", &event.p2mom);
  tree->Branch("p2mom_x", &event.p2mom_x);
  tree->Branch("p2mom_y", &event.p2mom_y);
  tree->Branch("p2mom_z", &event.p2mom_z);
  tree->Branch("GFp2mom", &event.GFp2mom);
  tree->Branch("GFp2mom_x", &event.GFp2mom_x);
  tree->Branch("GFp2mom_y", &event.GFp2mom_y);
  tree->Branch("GFp2mom_z", &event.GFp2mom_z);

  tree->Branch("Ksflag", &event.Ksflag);
  tree->Branch("Kspipid", &event.Kspipid);
  tree->Branch("Kspimid", &event.Kspimid);
  tree->Branch("Ksmass", &event.Ksmass);
  tree->Branch("Ksdecayvtx_x", &event.Ksdecayvtx_x);
  tree->Branch("Ksdecayvtx_y", &event.Ksdecayvtx_y);
  tree->Branch("Ksdecayvtx_z", &event.Ksdecayvtx_z);
  tree->Branch("Ksmom", &event.Ksmom);
  tree->Branch("Ksmom_x", &event.Ksmom_x);
  tree->Branch("Ksmom_y", &event.Ksmom_y);
  tree->Branch("Ksmom_z", &event.Ksmom_z);
  tree->Branch("KsDist", &event.Ksdist);
  tree->Branch("KFpvalKs", &event.KFpvalKs);
  tree->Branch("KFpullKs", &event.KFpullKs);
  tree->Branch("KFKsmom", &event.KFKsmom);
  tree->Branch("KFKsmom_x", &event.KFKsmom_x);
  tree->Branch("KFKsmom_y", &event.KFKsmom_y);
  tree->Branch("KFKsmom_z", &event.KFKsmom_z);

  tree->Branch( "KsPcnt", &event.KsPcnt );
  tree->Branch( "KsP1flag", &event.KsP1flag );
  tree->Branch( "KsP1pid", &event.KsP1pid );
  tree->Branch( "KsP1mom", &event.KsP1mom );
  tree->Branch( "KsP1mom_x", &event.KsP1mom_x );
  tree->Branch( "KsP1mom_y", &event.KsP1mom_y );
  tree->Branch( "KsP1mom_z", &event.KsP1mom_z );
  tree->Branch( "KsP1mass", &event.KsP1mass );
  tree->Branch( "KsP1decayvtx_x", &event.KsP1decayvtx_x );
  tree->Branch( "KsP1decayvtx_y", &event.KsP1decayvtx_y );
  tree->Branch( "KsP1decayvtx_z", &event.KsP1decayvtx_z );
  tree->Branch( "KsP1dist", &event.KsP1dist );

  tree->Branch( "KsP2flag", &event.KsP2flag );
  tree->Branch( "KsP2pid", &event.KsP2pid );
  tree->Branch( "KsP2mom", &event.KsP2mom );
  tree->Branch( "KsP2mom_x", &event.KsP2mom_x );
  tree->Branch( "KsP2mom_y", &event.KsP2mom_y );
  tree->Branch( "KsP2mom_z", &event.KsP2mom_z );
  tree->Branch( "KsP2mass", &event.KsP2mass );
  tree->Branch( "KsP2decayvtx_x", &event.KsP2decayvtx_x );
  tree->Branch( "KsP2decayvtx_y", &event.KsP2decayvtx_y );
  tree->Branch( "KsP2decayvtx_z", &event.KsP2decayvtx_z );  
  tree->Branch( "KsP2dist", &event.KsP2dist );
  

  tree->Branch("ForcedKsflag", &event.ForcedKsflag);
  tree->Branch("GFForcedKsflag", &event.GFForcedKsflag);
  tree->Branch("ForcedKsmass", &event.ForcedKsmass);
  tree->Branch("ForcedKsdecayvtx_x", &event.ForcedKsdecayvtx_x);
  tree->Branch("ForcedKsdecayvtx_y", &event.ForcedKsdecayvtx_y);
  tree->Branch("ForcedKsdecayvtx_z", &event.ForcedKsdecayvtx_z);
  tree->Branch("ForcedKsmom", &event.ForcedKsmom);
  tree->Branch("ForcedKsmom_x", &event.ForcedKsmom_x);
  tree->Branch("ForcedKsmom_y", &event.ForcedKsmom_y);
  tree->Branch("ForcedKsmom_z", &event.ForcedKsmom_z);
  tree->Branch("ForcedKsDist", &event.ForcedKsdist);
  tree->Branch("Forcedpipvertexdist", &event.Forcedpipvertexdist);
  tree->Branch("Forcedpimvertexdist", &event.Forcedpimvertexdist);
  tree->Branch("Forcedpipmom", &event.Forcedpipmom);
  tree->Branch("Forcedpipmom_x", &event.Forcedpipmom_x);
  tree->Branch("Forcedpipmom_y", &event.Forcedpipmom_y);
  tree->Branch("Forcedpipmom_z", &event.Forcedpipmom_z);
  tree->Branch("Forcedpimmom", &event.Forcedpimmom);
  tree->Branch("Forcedpimmom_x", &event.Forcedpimmom_x);
  tree->Branch("Forcedpimmom_y", &event.Forcedpimmom_y);
  tree->Branch("Forcedpimmom_z", &event.Forcedpimmom_z);
  tree->Branch("ForcedKFpvalKs", &event.ForcedKFpvalKs);
  tree->Branch("ForcedKFpullKs", &event.ForcedKFpullKs);
  tree->Branch("ForcedKFKsmom", &event.ForcedKFKsmom);
  tree->Branch("ForcedKFKsmom_x", &event.ForcedKFKsmom_x);
  tree->Branch("ForcedKFKsmom_y", &event.ForcedKFKsmom_y);
  tree->Branch("ForcedKFKsmom_z", &event.ForcedKFKsmom_z);
  tree->Branch("ForcedKFpipmom", &event.ForcedKFpipmom);
  tree->Branch("ForcedKFpipmom_x", &event.ForcedKFpipmom_x);
  tree->Branch("ForcedKFpipmom_y", &event.ForcedKFpipmom_y);
  tree->Branch("ForcedKFpipmom_z", &event.ForcedKFpipmom_z);
  tree->Branch("ForcedKFpimmom", &event.ForcedKFpimmom);
  tree->Branch("ForcedKFpimmom_x", &event.ForcedKFpimmom_x);
  tree->Branch("ForcedKFpimmom_y", &event.ForcedKFpimmom_y);
  tree->Branch("ForcedKFpimmom_z", &event.ForcedKFpimmom_z);



  tree->Branch("ForcedKsP1flag", &event.ForcedKsP1flag);
  tree->Branch("ForcedKsP1mass", &event.ForcedKsP1mass);
  tree->Branch("ForcedKsP1decayvtx_x", &event.ForcedKsP1decayvtx_x);
  tree->Branch("ForcedKsP1decayvtx_y", &event.ForcedKsP1decayvtx_y);
  tree->Branch("ForcedKsP1decayvtx_z", &event.ForcedKsP1decayvtx_z);
  tree->Branch("ForcedKsP1mom", &event.ForcedKsP1mom);
  tree->Branch("ForcedKsP1mom_x", &event.ForcedKsP1mom_x);
  tree->Branch("ForcedKsP1mom_y", &event.ForcedKsP1mom_y);
  tree->Branch("ForcedKsP1mom_z", &event.ForcedKsP1mom_z);
  tree->Branch("ForcedKsP1dist", &event.ForcedKsP1dist);
  tree->Branch("ForcedKsP1vertexdist", &event.ForcedKsP1vertexdist);
  
  tree->Branch("ForcedKsP2flag", &event.ForcedKsP2flag);
  tree->Branch("ForcedKsP2mass", &event.ForcedKsP2mass);
  tree->Branch("ForcedKsP2decayvtx_x", &event.ForcedKsP2decayvtx_x);
  tree->Branch("ForcedKsP2decayvtx_y", &event.ForcedKsP2decayvtx_y);
  tree->Branch("ForcedKsP2decayvtx_z", &event.ForcedKsP2decayvtx_z);
  tree->Branch("ForcedKsP2mom", &event.ForcedKsP2mom);
  tree->Branch("ForcedKsP2mom_x", &event.ForcedKsP2mom_x);
  tree->Branch("ForcedKsP2mom_y", &event.ForcedKsP2mom_y);
  tree->Branch("ForcedKsP2mom_z", &event.ForcedKsP2mom_z);
  tree->Branch("ForcedKsP2dist", &event.ForcedKsP2dist);
  tree->Branch("ForcedKsP2vertexdist", &event.ForcedKsP2vertexdist);




  TTreeReaderCont[kE42] = new TTreeReader( "tpc", TFileCont[kE42] );
  const auto& reader = TTreeReaderCont[kE42];
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );

  TTreeCont[kE42]->SetBranchAddress("nhFtof",&src.nhFtof);
  TTreeCont[kE42]->SetBranchAddress("didFtof",src.didFtof);
  TTreeCont[kE42]->SetBranchAddress("tFtof",src.tFtof);
  TTreeCont[kE42]->SetBranchAddress("deFtof",src.deFtof);
  TTreeCont[kE42]->SetBranchAddress("yFtof",src.yFtof);
  TTreeCont[kE42]->SetBranchAddress("nhittpc",&src.nhittpc);
  TTreeCont[kE42]->SetBranchAddress("inc",&src.inc);
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
