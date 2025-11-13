// -*- C++ -*-
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <TLorentzVector.h>

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
#include "TPCVertex.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#include "TF1.h"
#include "PidCommon.hh"

#define MakePidFig 1

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
  //const auto& gPidLike = PidLikelihoodMan::GetInstance();
const auto& zHSCenter = gGeom.LocalZ("HS");
const Double_t truncatedMean = 0.8; //80%

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");
  
const int nbininvbeta = 100;
const double mininvbeta = 0.;
const double maxinvbeta = 5.;
const int nbinbe = 400;
const double minbe = -0.5; //GeV
const double maxbe = 0.5;
const int pdfHid = 1000000;
// const double momstep = 0.05; //GeV/c
// const int momId = (Int_t)(maxpoq/momstep); // 0.05 GeV/c step?

const int nbinpoq = pidlikeli::nbinpoq;
const double minpoq = pidlikeli::minpoq;
const double maxpoq = pidlikeli::maxpoq; //GeV/c
const int nbindedx = pidlikeli::nbindedx;
const double mindedx = pidlikeli::mindedx;
const double maxdedx = pidlikeli::maxdedx;
const int nbinm2 = pidlikeli::nbinm2;
const double minm2 = pidlikeli::minm2;
const double maxm2 = pidlikeli::maxm2;
  
const double bestep = 0.20; //GeV // should be chaged to USER parameter 
const int BEId = (Int_t)((maxbe-minbe)/bestep); // should be less than 

const int kPmin = pidlikeli::kPmin;
const int kPmax = pidlikeli::kPmax;
const double kMomstep = pidlikeli::kDP;
const int kNmom = pidlikeli::kNmom;
const int kNpid = pidlikeli::kNpid;
const int piHid = pidlikeli::kPion;
const int kHid = pidlikeli::kKaon;
const int pHid = pidlikeli::kProton;  
const int dHid = pidlikeli::kDeutron;
const int eHid = pidlikeli::kElectron;
const int allHid = pidlikeli::kAllParticles;      
const int kNchg = pidlikeli::kNchg;
const int plusHid = pidlikeli::kPlus;
const int minusHid = pidlikeli::kMinus;
const int kNtype = pidlikeli::kNtype;
const int kNbe = pidlikeli::kNbe;
const auto& plist = pidlikeli::plist;
const auto& clist = pidlikeli::clist;
const auto& type = pidlikeli::type;
const int fac_t = pidlikeli::fac_t;
const int fac_p = pidlikeli::fac_p;
const int fac_c = pidlikeli::fac_c;
const int fac_b = pidlikeli::fac_b;
const int fac_m = pidlikeli::fac_m;
const int typeGenHid = pidlikeli::kTypeGen;  
const int typeLHid = pidlikeli::kTypeLmd;
const int typeK0Hid = pidlikeli::kTypeK0;
const int typeKmHid = pidlikeli::kTypeKm;

const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;

const Double_t lambda_masscut = 0.1;
const Double_t lambda_masscut_final = 0.02; //final
const Double_t xi_masscut = 0.1;  
const Double_t k0_masscut = 0.1; //final
//const Double_t xi_masscut = 0.15; const Double_t lambda_masscut = 0.1; //ref
const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t e_vtx_distcut = 100;
const Double_t ppi_distcut = 10.; //ref
//const Double_t lpi_distcut = 10.;
const Double_t ltarget_distcut = 25.;
const Double_t ltarget_ycut = 20.;  
const Double_t pip_vtx_distcut = 300;
const Double_t pim_vtx_distcut = 300;
const Double_t pipi_distcut = 10.; //ref  
const Double_t k0target_distcut = 25.;

const Double_t lpi_distcut = 15.;
const Double_t pi2_vtx_distcut = 300;
const Double_t xitarget_distcut = 50.; //ref
  
const Double_t GFppi_distcut = 10.;  
const Double_t GFlpi_distcut = 10.;
//const Double_t GFlpi_distcut = 15.;
const Double_t GFltarget_distcut = 25.;
const Double_t GFltarget_ycut = 20.;
const Double_t GFpipi_distcut = 10.;
const Double_t GFk0target_distcut = 25.;  
const Double_t GFk0target_ycut = 20.;
const Double_t GFxitarget_ycut = 20.;

const Double_t residual_track_distcut = 25.;
const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

const Double_t cutm2proton = 0.3;
const Double_t tTPCcut = 10.0;

const Bool_t kp=false;
const Bool_t ch2kk = true;
const Bool_t forGeant4 = true;

const double minMM = 0.38;
const double maxMM = 0.58;
// const double minM2Km = -0.20;
// const double maxM2Km =  0.50;

const double minThetaKP = 3.5;
const double maxThetaKP = 4.5;

  //const Int_t MaxHits = 2000;

constexpr double cTau0L = 7.89; // cTauLambda [cm]
constexpr double c_cm_per_ns = 29.9792458; // c[cm/ns]

const auto& psTrigA = ConfMan::Get<Int_t>("PSTRGA");
const auto& psTrigB = ConfMan::Get<Int_t>("PSTRGB");

  const double maxbek = 0.300;
  const double minbek = -0.040;  
  const double onebinbek = 0.005*4.;
  const double numbinbek = (maxbek-minbek)/onebinbek;

  int GetBinIndexBEk(double x){
    double xmin=-0.3;
    double xmax=0.3;
    double nbins=120.;
    
    if(x<xmin||x>xmax) return -1;
    double binwidth = (xmax-xmin)/nbins;    
    int bin = static_cast<int>(std::floor((x-xmin)/binwidth))+1;
    if(bin>nbins) bin = nbins;
    return bin;
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
  { "[Process]", "[ConfFile]", "[DstE42]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct DstG4
{
  // for Geant4
  Int_t ich;
  Double_t bpx;
  Double_t bpy;
  Double_t bpz;
  Double_t kppx;
  Double_t kppy;
  Double_t kppz;  
  Int_t np;
  Int_t pidpdg[MaxHits];
  Double_t px[MaxHits];
  Double_t py[MaxHits];
  Double_t pz[MaxHits];
  // std::vector<Double_t> px;
  // std::vector<Double_t> py;
  // std::vector<Double_t> pz;

  void clear( void )
  {
    
    ich = -1;
    bpx = qnan;
    bpy = qnan;
    bpz = qnan;
    kppx = qnan;
    kppy = qnan;
    kppz = qnan;    
    np = -1;
    for(Int_t ip=0; ip<MaxHits; ++ip){
      pidpdg[ip] = -1; 
      px[ip] = qnan;
      py[ip] = qnan;
      pz[ip] = qnan;
    }
  }  
};

struct Event
{
  Int_t status;
  Int_t runnum;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;

  Int_t ntK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> thetaK18;  
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> xtgtK18;
  std::vector<Double_t> ytgtK18;
  std::vector<Double_t> utgtK18;
  std::vector<Double_t> vtgtK18;

  Int_t ntKurama;
  std::vector<Double_t> chisqrKurama;
  std::vector<Double_t> pKurama;
  std::vector<Double_t> qKurama;
  std::vector<Double_t> m2;
  std::vector<Double_t> xtgtKurama;
  std::vector<Double_t> ytgtKurama;
  std::vector<Double_t> utgtKurama;
  std::vector<Double_t> vtgtKurama;

  Int_t nKm;
  Int_t nKp;
  Int_t nKK;
  std::vector<Int_t> inside;
  std::vector<Double_t> vtx;
  std::vector<Double_t> vty;
  std::vector<Double_t> vtz;
  std::vector<Double_t> closeDist;
  std::vector<Double_t> MissMass;
  std::vector<Double_t> MissMassCorr;
  std::vector<Double_t> MissMassCorrDE;
  std::vector<Int_t> Kflag;
  std::vector<Int_t> Pflag;
  std::vector<Int_t> Heavyflag;  
  
  Int_t ntTPCK18;
  std::vector<Int_t> isgoodTPCK18;
  std::vector<Int_t> tpcidTPCK18;
  std::vector<Double_t> chisqrTPCK18;
  std::vector<Double_t> pTPCK18;
  std::vector<Double_t> qTPCK18;
  std::vector<Double_t> thetaTPCK18;  
  std::vector<Double_t> xtgtTPCK18;
  std::vector<Double_t> ytgtTPCK18;
  std::vector<Double_t> utgtTPCK18;
  std::vector<Double_t> vtgtTPCK18;
  std::vector<Double_t> lhtofTPCK18;
  std::vector<Double_t> xhtofTPCK18;
  std::vector<Double_t> yhtofTPCK18;
  std::vector<std::vector<Double_t>> lvpTPCK18;
  std::vector<std::vector<Double_t>> xvpTPCK18;
  std::vector<std::vector<Double_t>> yvpTPCK18;

  std::vector<Int_t> isgoodTPCKurama;
  std::vector<Int_t> tpcidTPCKurama;
  std::vector<Int_t> kflagTPCKurama;
  std::vector<Int_t> pflagTPCKurama;
  std::vector<Double_t> chisqrTPCKurama;
  std::vector<Double_t> pTPCKurama;
  std::vector<Double_t> qTPCKurama;
  std::vector<Double_t> m2TPCKurama;
  std::vector<Double_t> xtgtTPCKurama;
  std::vector<Double_t> ytgtTPCKurama;
  std::vector<Double_t> utgtTPCKurama;
  std::vector<Double_t> vtgtTPCKurama;
  std::vector<Double_t> thetaTPCKurama;
  std::vector<Double_t> pathTPCKurama;
  std::vector<Double_t> lhtofTPCKurama;
  std::vector<Double_t> xhtofTPCKurama;
  std::vector<Double_t> yhtofTPCKurama;
  std::vector<std::vector<Double_t>> lvpTPCKurama;
  std::vector<std::vector<Double_t>> xvpTPCKurama;
  std::vector<std::vector<Double_t>> yvpTPCKurama;

  std::vector<Int_t> isgoodTPC;
  std::vector<Int_t> insideTPC;
  std::vector<Double_t> vtxTPC;
  std::vector<Double_t> vtyTPC;
  std::vector<Double_t> vtzTPC;
  std::vector<Double_t> closeDistTPC;
  std::vector<Double_t> MissMassTPC;
  std::vector<Double_t> MissMassCorrTPC;
  std::vector<Double_t> MissMassCorrDETPC;
  std::vector<Double_t> MissMassNuclTPC;
  std::vector<Double_t> MissMassNuclCorrTPC;
  std::vector<Double_t> MissMassNuclCorrDETPC;
  std::vector<Double_t> BEkaonTPC;

  std::vector<Double_t> pOrg;
  std::vector<Double_t> pCalc;
  std::vector<Double_t> pCorr;
  std::vector<Double_t> pCorrDE;
  std::vector<Double_t> xb;
  std::vector<Double_t> yb;
  std::vector<Double_t> ub;
  std::vector<Double_t> vb;
  std::vector<Double_t> xs;
  std::vector<Double_t> ys;
  std::vector<Double_t> us;
  std::vector<Double_t> vs;

  std::vector<Double_t> pOrgTPC;
  std::vector<Double_t> pCalcTPC;
  std::vector<Double_t> pCorrTPC;
  std::vector<Double_t> pCorrDETPC;
  std::vector<Double_t> thetaTPC;
  std::vector<Double_t> thetaCMTPC;
  std::vector<Double_t> costCMTPC;
  std::vector<Double_t> xbTPC;
  std::vector<Double_t> ybTPC;
  std::vector<Double_t> ubTPC;
  std::vector<Double_t> vbTPC;
  std::vector<Double_t> xsTPC;
  std::vector<Double_t> ysTPC;
  std::vector<Double_t> usTPC;
  std::vector<Double_t> vsTPC;

  std::vector<Double_t> BE;
  std::vector<Double_t> BETPC;
  std::vector<Double_t> BE_LL;
  std::vector<Double_t> BETPC_LL;
  std::vector<Double_t> km_mom_x;
  std::vector<Double_t> km_mom_y;
  std::vector<Double_t> km_mom_z;
  std::vector<Double_t> kp_mom_x;
  std::vector<Double_t> kp_mom_y;
  std::vector<Double_t> kp_mom_z;  
  
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
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Double_t> chisqr;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx;
  std::vector<Double_t> mom0;//Helix momentum at Y = 0
  std::vector<Int_t> charge;//Helix charge
  std::vector<Double_t> path;//Helix path
  std::vector<Int_t> pid;
  std::vector<Int_t> isElectron;
  std::vector<Double_t> nsigma_triton;
  std::vector<Double_t> nsigma_deutron;
  std::vector<Double_t> nsigma_proton;
  std::vector<Double_t> nsigma_kaon;
  std::vector<Double_t> nsigma_pion;
  std::vector<Double_t> nsigma_electron;

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

  Int_t lflag;
  Int_t kuramalflag;

  Int_t k0flag;
  Double_t k0mass;
  Double_t k0decayvtx_x;
  Double_t k0decayvtx_y;
  Double_t k0decayvtx_z;
  Double_t k0mom_x;
  Double_t k0mom_y;
  Double_t k0mom_z;
  Double_t pipi_dist;
  Double_t pipiangle;    
  std::vector<Int_t> k0decays_id;
  std::vector<Double_t> k0decays_mom;
  std::vector<Double_t> k0decays_mom_x;
  std::vector<Double_t> k0decays_mom_y;
  std::vector<Double_t> k0decays_mom_z;    
  std::vector<Double_t> k0decays_mass2;
  std::vector<Double_t> k0decays_tracklen;
  std::vector<Double_t> k0decays_htofpos_x;
  std::vector<Double_t> k0decays_htofpos_y;
  std::vector<Double_t> k0decays_htofpos_z;

  Int_t nhHtof;
  std::vector<Double_t> HtofSeg;
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

  Int_t kmflag;
  Int_t kminc;    
  Int_t    GFkmid;  
  Double_t GFkmmom;
  Double_t GFkmmom_x;
  Double_t GFkmmom_y;
  Double_t GFkmmom_z;
  Double_t GFkmtheta;
  Double_t GFkmphi;    
  Double_t GFkmtarget_dist;
  Double_t GFkmtargetvtx_x;
  Double_t GFkmtargetvtx_y;
  Double_t GFkmtargetvtx_z;
  Double_t GFkmtargetcenter_x;
  Double_t GFkmtargetcenter_y;
  Double_t GFkmtargetcenter_z;
  Double_t GFkmtargetcenter_dist;
  Double_t GFkmprodvtx_x;
  Double_t GFkmprodvtx_y;
  Double_t GFkmprodvtx_z;
  Double_t GFkmprodvtx_dist;
  Int_t    GFkmhtofid;
  Int_t    GFkmhtofseg;    
  Double_t GFkmtracklen;
  Double_t GFkmtof;
  Double_t GFkmmass2;
  Double_t GFkminvbeta;
  std::vector<Double_t> GFkmposHtof;
    
  Double_t GFk0mass;
  Double_t GFk0decayvtx_x;
  Double_t GFk0decayvtx_y;
  Double_t GFk0decayvtx_z;
  Double_t GFk0mom;
  Double_t GFk0mom_x;
  Double_t GFk0mom_y;
  Double_t GFk0mom_z;
  Double_t GFk0pipi_dist;
  Double_t GFk0target_dist;
  Double_t GFk0targetvtx_x;
  Double_t GFk0targetvtx_y;
  Double_t GFk0targetvtx_z;
  Double_t GFk0targetcenter_x;
  Double_t GFk0targetcenter_y;
  Double_t GFk0targetcenter_z;
  Double_t GFk0targetcenter_dist;
  Double_t GFk0prodvtx_x;
  Double_t GFk0prodvtx_y;
  Double_t GFk0prodvtx_z;
  Double_t GFk0prodvtx_dist;
  Double_t GFk0tracklen;
  Double_t GFk0tof;

  std::vector<Double_t> GFk0decays_id;
  std::vector<Double_t> GFk0decays_mass2;
  std::vector<Double_t> GFk0decays_invbeta;
  std::vector<Double_t> GFk0decays_mom;
  
  Int_t GFinsideKurama;
  Int_t GFKuramaVtxOutTgt;
  Int_t GFFwdYDecayTrack;
  std::vector<Int_t> GFinside;
  std::vector<Int_t> GFfromVtx;  

  Int_t GFntTpc_inside;
  Double_t GFprodvtx_x;
  Double_t GFprodvtx_y;
  Double_t GFprodvtx_z;

  std::vector<Int_t> GFextrapolationHtof;  
  std::vector<Double_t> GFtracklen;
  std::vector<Double_t> GFtrack2vtxdist;
  std::vector<Double_t> GFcalctof;
  std::vector<Double_t> GFsegHtof;
  std::vector<Double_t> GFtofHtof;
  std::vector<Double_t> GFtdiffHtof;
  std::vector<Double_t> GFposHtof;
  std::vector<Double_t> GFposx;
  std::vector<Double_t> GFposy;
  std::vector<Double_t> GFposz;
  std::vector<Double_t> GFinvbeta;
  std::vector<Double_t> GFm2;
  std::vector<Double_t> nsigma_tritonHtof;
  std::vector<Double_t> nsigma_deutronHtof;
  std::vector<Double_t> nsigma_protonHtof;
  std::vector<Double_t> nsigma_kaonHtof;
  std::vector<Double_t> nsigma_pionHtof;
  std::vector<Double_t> nsigma_electronHtof;  

  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();

    ntK18 = 0;
    pK18.clear();
    thetaK18.clear();  
    chisqrK18.clear();
    xtgtK18.clear();
    ytgtK18.clear();
    utgtK18.clear();
    vtgtK18.clear();

    ntKurama = 0;
    chisqrKurama.clear();    
    pKurama.clear();
    qKurama.clear();    
    m2.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();
    
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
    isK18.clear();
    isAccidental.clear();
    chisqr.clear();
    helix_cx.clear();
    helix_cy.clear();
    helix_z0.clear();
    helix_r.clear();
    helix_dz.clear();
    dE.clear();
    dEdx.clear();
    mom0.clear();
    charge.clear();
    path.clear();
    pid.clear();
    isElectron.clear();
    nsigma_triton.clear();
    nsigma_deutron.clear();
    nsigma_proton.clear();
    nsigma_kaon.clear();
    nsigma_pion.clear();
    nsigma_electron.clear();

    nKp = 0;        
    nKm = 0;    
    nKK = 0;
    inside.clear();
    vtx.clear();
    vty.clear();
    vtz.clear();
    closeDist.clear();
    MissMass.clear();
    MissMassCorr.clear();
    MissMassCorrDE.clear();
    Kflag.clear();
    Pflag.clear();
    Heavyflag.clear();    
    BE.clear();

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
    track_cluster_size.clear();
    track_cluster_mrow.clear();
    track_cluster_de_center.clear();
    track_cluster_x_center.clear();
    track_cluster_y_center.clear();
    track_cluster_z_center.clear();
    track_cluster_row_center.clear();

    isgoodTPCKurama.clear();
    isgoodTPCKurama.clear();
    tpcidTPCKurama.clear();
    kflagTPCKurama.clear();
    pflagTPCKurama.clear();
    chisqrTPCKurama.clear();
    
    pTPCKurama.clear();
    qTPCKurama.clear();
    m2TPCKurama.clear();    

    isgoodTPC.clear();
    insideTPC.clear();
    vtxTPC.clear();
    vtyTPC.clear();
    vtzTPC.clear();
    closeDistTPC.clear();
    MissMassTPC.clear();
    MissMassCorrTPC.clear();
    MissMassCorrDETPC.clear();
    MissMassNuclTPC.clear();
    MissMassNuclCorrTPC.clear();
    MissMassNuclCorrDETPC.clear();
    BEkaonTPC.clear();

        pOrg.clear();
    pCalc.clear();
    pCorr.clear();
    pCorrDE.clear();
    xb.clear();
    yb.clear();
    ub.clear();
    vb.clear();
    xs.clear();
    ys.clear();
    us.clear();
    vs.clear();

    pOrgTPC.clear();
    pCalcTPC.clear();
    pCorrTPC.clear();
    pCorrDETPC.clear();
    thetaTPC.clear();
    thetaCMTPC.clear();
    costCMTPC.clear();
    xbTPC.clear();
    ybTPC.clear();
    ubTPC.clear();
    vbTPC.clear();
    xsTPC.clear();
    ysTPC.clear();
    usTPC.clear();
    vsTPC.clear();    

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

    lflag = false;
    kuramalflag = false;
    
    k0flag = false;
    k0mass = qnan;
    k0decayvtx_x = qnan;
    k0decayvtx_y = qnan;
    k0decayvtx_z = qnan;
    k0mom_x = qnan;
    k0mom_y = qnan;
    k0mom_z = qnan;
    pipi_dist = qnan;
    pipiangle = qnan;
    k0decays_id.clear();
    k0decays_mom.clear();
    k0decays_mom_x.clear();
    k0decays_mom_y.clear();
    k0decays_mom_z.clear();
    k0decays_mass2.clear();
    k0decays_tracklen.clear();
    k0decays_htofpos_x.clear();
    k0decays_htofpos_y.clear();
    k0decays_htofpos_z.clear();

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

    kmflag = false;
    kminc = false;    
    GFkmid = -1;    
    GFkmmom = qnan;
    GFkmmom_x = qnan;
    GFkmmom_y = qnan;
    GFkmmom_z = qnan;
    GFkmtheta = qnan;
    GFkmphi = qnan;        
    GFkmtarget_dist = qnan;
    GFkmtargetvtx_x = qnan;
    GFkmtargetvtx_y = qnan;
    GFkmtargetvtx_z = qnan;
    GFkmtargetcenter_dist = qnan;
    GFkmtargetcenter_x = qnan;
    GFkmtargetcenter_y = qnan;
    GFkmtargetcenter_z = qnan;
    GFkmprodvtx_x = qnan;
    GFkmprodvtx_y = qnan;
    GFkmprodvtx_z = qnan;
    GFkmprodvtx_dist = qnan;
    GFkmhtofid = -1;
    GFkmhtofseg = -1;            
    GFkmtracklen = qnan;    
    GFkmtof = qnan;
    GFkmmass2 = qnan;
    GFkminvbeta = qnan;
    GFkmposHtof.clear();

    GFk0mass = qnan;
    GFk0decayvtx_x = qnan;
    GFk0decayvtx_y = qnan;
    GFk0decayvtx_z = qnan;
    GFk0mom = qnan;
    GFk0mom_x = qnan;
    GFk0mom_y = qnan;
    GFk0mom_z = qnan;
    GFk0tracklen = qnan;
    GFk0tof = qnan;
    GFk0pipi_dist = qnan;
    GFk0target_dist = qnan;
    GFk0targetvtx_x = qnan;
    GFk0targetvtx_y = qnan;
    GFk0targetvtx_z = qnan;
    GFk0targetcenter_dist = qnan;
    GFk0targetcenter_x = qnan;
    GFk0targetcenter_y = qnan;
    GFk0targetcenter_z = qnan;
    GFk0prodvtx_x = qnan;
    GFk0prodvtx_y = qnan;
    GFk0prodvtx_z = qnan;
    GFk0prodvtx_dist = qnan;
    GFk0tracklen = qnan;
    GFk0tof = qnan;

    GFk0decays_id.clear();    
    GFk0decays_mass2.clear();
    GFk0decays_invbeta.clear();
    GFk0decays_mom.clear();
    
    GFinsideKurama = 0;
    GFKuramaVtxOutTgt = 0;
    GFFwdYDecayTrack = -1;    
    GFinside.clear();
    GFfromVtx.clear();    

    GFntTpc_inside = 0;
    GFprodvtx_x = qnan;
    GFprodvtx_y = qnan;
    GFprodvtx_z = qnan;

    GFextrapolationHtof.clear();    
    GFtracklen.clear();
    GFtrack2vtxdist.clear();
    GFcalctof.clear();
    GFsegHtof.clear();
    GFtofHtof.clear();
    GFtdiffHtof.clear();
    GFposHtof.clear();
    GFposx.clear();
    GFposy.clear();
    GFposz.clear();
    GFinvbeta.clear();
    GFm2.clear();
    nsigma_tritonHtof.clear();
    nsigma_deutronHtof.clear();
    nsigma_protonHtof.clear();
    nsigma_kaonHtof.clear();
    nsigma_pionHtof.clear();
    nsigma_electronHtof.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

  TTreeReaderValue<Int_t>* ntK18;
  TTreeReaderValue<std::vector<Double_t>>* pK18;
  TTreeReaderValue<std::vector<Double_t>>* thetaK18;  
  TTreeReaderValue<std::vector<Double_t>>* chisqrK18;
  TTreeReaderValue<std::vector<Double_t>>* xtgtK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtK18;

  TTreeReaderValue<Int_t>* ntKurama;
  TTreeReaderValue<std::vector<Double_t>>* chisqrKurama;
  TTreeReaderValue<std::vector<Double_t>>* pKurama;
  TTreeReaderValue<std::vector<Double_t>>* qKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2;
  TTreeReaderValue<std::vector<Double_t>>* xtgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* ytgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* utgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* vtgtKurama;

  TTreeReaderValue<Int_t>* nKm;
  TTreeReaderValue<Int_t>* nKp;  
  TTreeReaderValue<Int_t>* nKK;
  TTreeReaderValue<std::vector<Int_t>>* inside;
  TTreeReaderValue<std::vector<Double_t>>* vtx;
  TTreeReaderValue<std::vector<Double_t>>* vty;
  TTreeReaderValue<std::vector<Double_t>>* vtz;
  TTreeReaderValue<std::vector<Double_t>>* closeDist;
  TTreeReaderValue<std::vector<Double_t>>* MissMass;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorr;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorrDE;
  TTreeReaderValue<std::vector<Int_t>>* Kflag;
  TTreeReaderValue<std::vector<Int_t>>* Pflag;
  TTreeReaderValue<std::vector<Int_t>>* Heavyflag;

  TTreeReaderValue<Int_t>* ntTPCK18; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* tpcidTPCK18;
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* qTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* pTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPCK18;  
  TTreeReaderValue<std::vector<Double_t>>* xtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* lhtofTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xhtofTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* yhtofTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* lvpTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpTPCK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xhtofK18;
  TTreeReaderValue<std::vector<Double_t>>* yhtofK18;

  TTreeReaderValue<Int_t>* ntTPCKurama; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* tpcidTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* kflagTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* pflagTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* pTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* qTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2TPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xsTPC;
  TTreeReaderValue<std::vector<Double_t>>* ysTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;
  TTreeReaderValue<std::vector<Double_t>>* xbTPC;
  TTreeReaderValue<std::vector<Double_t>>* ybTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;  
  TTreeReaderValue<std::vector<Double_t>>* xtgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* ytgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* utgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* vtgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* pathTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* lhtofTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xhtofTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* yhtofTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* lvpTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpTPCKurama;

  TTreeReaderValue<std::vector<Int_t>>* isgoodTPC;
  TTreeReaderValue<std::vector<Int_t>>* insideTPC;
  TTreeReaderValue<std::vector<Double_t>>* vtxTPC;
  TTreeReaderValue<std::vector<Double_t>>* vtyTPC;
  TTreeReaderValue<std::vector<Double_t>>* vtzTPC;
  TTreeReaderValue<std::vector<Double_t>>* closeDistTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassNuclTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassNuclCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* MissMassNuclCorrDETPC;

  TTreeReaderValue<std::vector<Double_t>>* pOrg;
  TTreeReaderValue<std::vector<Double_t>>* pCalc;
  TTreeReaderValue<std::vector<Double_t>>* pCorr;
  TTreeReaderValue<std::vector<Double_t>>* pCorrDE;
  TTreeReaderValue<std::vector<Double_t>>* xb;
  TTreeReaderValue<std::vector<Double_t>>* yb;
  TTreeReaderValue<std::vector<Double_t>>* ub;
  TTreeReaderValue<std::vector<Double_t>>* vb;
  TTreeReaderValue<std::vector<Double_t>>* xs;
  TTreeReaderValue<std::vector<Double_t>>* ys;
  TTreeReaderValue<std::vector<Double_t>>* us;
  TTreeReaderValue<std::vector<Double_t>>* vs;
  
  TTreeReaderValue<std::vector<Double_t>>* pOrgTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCalcTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* costCMTPC;

  TTreeReaderValue<Int_t>* nclTpc; // Number of clusters
  TTreeReaderValue<std::vector<Double_t>>* cluster_x;
  TTreeReaderValue<std::vector<Double_t>>* cluster_y;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de;
  TTreeReaderValue<std::vector<Int_t>>* cluster_size;
  TTreeReaderValue<std::vector<Int_t>>* cluster_layer;
  TTreeReaderValue<std::vector<Double_t>>* cluster_mrow;
  TTreeReaderValue<std::vector<Double_t>>* cluster_de_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_x_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_y_center;
  TTreeReaderValue<std::vector<Double_t>>* cluster_z_center;
  TTreeReaderValue<std::vector<Int_t>>* cluster_row_center;
  TTreeReaderValue<std::vector<Int_t>>* cluster_houghflag;

  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge  
  TTreeReaderValue<std::vector<Int_t>>* pid; 
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
  TTreeReaderValue<std::vector<Int_t>>* isElectron;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_triton;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_deutron;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_proton;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_kaon;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_pion;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_electron;  

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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* alpha;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* pathhit;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<Int_t>* nvtxTpc;
  TTreeReaderValue<std::vector<Double_t>>* vtx_x;
  TTreeReaderValue<std::vector<Double_t>>* vtx_y;
  TTreeReaderValue<std::vector<Double_t>>* vtx_z;
  TTreeReaderValue<std::vector<Double_t>>* vtx_dist;
  TTreeReaderValue<std::vector<Double_t>>* vtx_angle;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxid;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_theta;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxmom_z;

  TTreeReaderValue<Int_t>* nhHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;

  TTreeReaderValue<std::vector<Int_t>>* isLambda;
  TTreeReaderValue<std::vector<Int_t>>* ncombiLambda;
  TTreeReaderValue<std::vector<Double_t>>* distLambda;
  TTreeReaderValue<std::vector<Double_t>>* angleLambda;
  TTreeReaderValue<std::vector<Double_t>>* bestmassLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* massLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* momLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysidLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* decaysmomLambda_z;

  TTreeReaderValue<Int_t>*   lflag;  
};

namespace root
{
Event  event;
DstG4  dstg4;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
TTree *treetpc;
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
    if( DstRead( ievent ) ){
      tree->Fill();
      treetpc->Fill();
    }
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
  dstg4.clear();  

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

  static const auto ElectronMass = pdg::ElectronMass();
  static const auto PionMass = pdg::PionMass();
  static const auto KaonMass = pdg::KaonMass();
  static const auto K0Mass = pdg::K0Mass();  
  static const auto PhiMass = pdg::Mass(333);
  static const auto ProtonMass = pdg::ProtonMass();
  static const auto LambdaMass = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  static const auto XiMass = XiMinusMass;
  static const Double_t Carbon12Mass = 12.*TGeoUnit::amu_c2 - 6.*ElectronMass;
  static const Double_t Boron11Mass  = 11.009305167*TGeoUnit::amu_c2 - 5.*ElectronMass;
  //std::cout << " MassK: " << KaonMass << " Mass11B: " << Boron11Mass << std::endl;
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};
  TVector3 tgtpos(0, 0, tpc::ZTarget);
  TVector3 qnan_vec = TVector3(qnan, qnan, qnan);

  static const auto KKEvent = gUser.GetParameter("KKEvent");
  static const auto KPEvent = gUser.GetParameter("KPEvent");
  static const auto KHeavyEvent = gUser.GetParameter("KHeavyEvent");

  static const auto mindEdxSigKaon = gUser.GetParameter("MindEdXSigKaon");
  static const auto maxdEdxSigKaon = gUser.GetParameter("MaxdEdXSigKaon");
  static const auto minM2Km = gUser.GetParameter("MinM2Kaon");
  static const auto maxM2Km = gUser.GetParameter("MaxM2Kaon");    
  //const Double_t& mindEdxSigKaon = ConfMan::Get<Double_t>("MindEdXSigKaon");
  //const Double_t& maxdEdxSigKaon = ConfMan::Get<Double_t>("MaxdEdXSigKaon");    
  
  Double_t vtx_scan_range = gUser.GetParameter("VertexScanRange"); 
  
  if( ievent%1000==0 ){
  //  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  bool debug = true;
  //
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  event.nclTpc = **src.nclTpc;
  event.cluster_x = **src.cluster_x;
  event.cluster_y = **src.cluster_y;
  event.cluster_z = **src.cluster_z;
  event.cluster_de = **src.cluster_de;
  event.cluster_size = **src.cluster_size;
  event.cluster_layer = **src.cluster_layer;
  event.cluster_mrow = **src.cluster_mrow;
  event.cluster_de_center = **src.cluster_de_center;
  event.cluster_x_center = **src.cluster_x_center;
  event.cluster_y_center = **src.cluster_y_center;
  event.cluster_z_center = **src.cluster_z_center;
  event.cluster_row_center = **src.cluster_row_center;
  event.cluster_houghflag = **src.cluster_houghflag;

  event.nKm = **src.nKm;
  event.nKp = **src.nKp;
  event.nKK = **src.nKK;
  event.inside = **src.inside;
  event.vtx = **src.vtx;
  event.vty = **src.vty;
  event.vtz = **src.vtz;
  event.closeDist = **src.closeDist;  
  //  event.MissMass = **src.MissMass;
  event.MissMassCorr = **src.MissMassCorr;
  event.MissMassCorrDE = **src.MissMassCorrDE;
  event.Pflag = **src.Pflag;
  event.Kflag = **src.Kflag;
  event.Heavyflag = **src.Heavyflag;

  event.vtx = **src.vtx;
  event.vty = **src.vty;
  event.vtz = **src.vtz;
  event.closeDist = **src.closeDist;

  event.ntK18 = **src.ntK18;
  event.ntTPCK18 = **src.ntTPCK18;
  event.chisqrK18 = **src.chisqrK18;
  event.pK18 = **src.pK18;
  event.thetaK18 = **src.thetaK18;  
  event.xtgtK18 = **src.xtgtK18;
  event.ytgtK18 = **src.ytgtK18;
  event.utgtK18 = **src.utgtK18;
  event.vtgtK18 = **src.vtgtK18;

  event.hitlayer = **src.hitlayer;
  event.hitpos_x = **src.hitpos_x;
  event.hitpos_y = **src.hitpos_y;
  event.hitpos_z = **src.hitpos_z;
  event.calpos_x = **src.calpos_x;
  event.calpos_y = **src.calpos_y;
  event.calpos_z = **src.calpos_z;
  event.mom_x = **src.mom_x;
  event.mom_y = **src.mom_y;
  event.mom_z = **src.mom_z;
  event.residual = **src.residual;
  event.residual_x = **src.residual_x;
  event.residual_y = **src.residual_y;
  event.residual_z = **src.residual_z;
  event.resolution_x = **src.resolution_x;
  event.resolution_y = **src.resolution_y;
  event.resolution_z = **src.resolution_z;
  event.helix_t = **src.helix_t;
  event.alpha = **src.alpha;
  event.pathhit = **src.pathhit;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;
  event.track_cluster_de_center = **src.track_cluster_de_center;
  event.track_cluster_x_center = **src.track_cluster_x_center;
  event.track_cluster_y_center = **src.track_cluster_y_center;
  event.track_cluster_z_center = **src.track_cluster_z_center;

  event.track_cluster_row_center = **src.track_cluster_row_center;

  event.isgoodTPCK18 = **src.isgoodTPCK18;
  event.chisqrTPCK18 = **src.chisqrTPCK18;
  event.pTPCK18 = **src.pTPCK18;
  event.qTPCK18 = **src.qTPCK18;
  event.thetaTPCK18 = **src.thetaTPCK18;  
  event.xtgtTPCK18 = **src.xtgtTPCK18;
  event.ytgtTPCK18 = **src.ytgtTPCK18;
  event.utgtTPCK18 = **src.utgtTPCK18;
  event.vtgtTPCK18 = **src.vtgtTPCK18;
  event.lhtofTPCK18 = **src.lhtofTPCK18;
  event.xhtofTPCK18 = **src.xhtofTPCK18;
  event.yhtofTPCK18 = **src.yhtofTPCK18;
  event.lvpTPCK18 = **src.lvpTPCK18;
  event.xvpTPCK18 = **src.xvpTPCK18;
  event.yvpTPCK18 = **src.yvpTPCK18;

  event.ntKurama = **src.ntKurama;
  event.chisqrKurama = **src.chisqrKurama;
  event.pKurama = **src.pKurama;
  event.qKurama = **src.qKurama;
  event.xtgtKurama = **src.xtgtKurama;
  event.ytgtKurama = **src.ytgtKurama;
  event.utgtKurama = **src.utgtKurama;
  event.vtgtKurama = **src.vtgtKurama;
  
  event.tpcidTPCKurama = **src.tpcidTPCKurama;
  event.isgoodTPCKurama = **src.isgoodTPCKurama;
  event.kflagTPCKurama = **src.kflagTPCKurama;
  event.pflagTPCKurama = **src.pflagTPCKurama;
  event.chisqrTPCKurama = **src.chisqrTPCKurama;
  event.pTPCKurama = **src.pTPCKurama;
  event.qTPCKurama = **src.qTPCKurama;
  event.m2TPCKurama = **src.m2TPCKurama;
  event.xtgtTPCKurama = **src.xtgtTPCKurama;
  event.ytgtTPCKurama = **src.ytgtTPCKurama;
  event.utgtTPCKurama = **src.utgtTPCKurama;
  event.vtgtTPCKurama = **src.vtgtTPCKurama;
  event.thetaTPCKurama = **src.thetaTPCKurama;

  event.xsTPC = **src.xsTPC;
  event.ysTPC = **src.ysTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;
  event.xbTPC = **src.xbTPC;
  event.ybTPC = **src.ybTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;

  event.nvtxTpc = **src.nvtxTpc;
  event.vtx_x = **src.vtx_x;
  event.vtx_y = **src.vtx_y;
  event.vtx_z = **src.vtx_z;
  event.vtx_dist = **src.vtx_dist;
  event.vtx_angle = **src.vtx_angle;
  event.vtxid = **src.vtxid;
  event.vtxmom_theta = **src.vtxmom_theta;
  event.vtxpos_x = **src.vtxpos_x;
  event.vtxpos_y = **src.vtxpos_y;
  event.vtxpos_z = **src.vtxpos_z;
  event.vtxmom_x = **src.vtxmom_x;
  event.vtxmom_y = **src.vtxmom_y;
  event.vtxmom_z = **src.vtxmom_z;
  
  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;

  event.pathTPCKurama = **src.pathTPCKurama;
  event.lhtofTPCKurama = **src.lhtofTPCKurama;
  event.xhtofTPCKurama = **src.xhtofTPCKurama;
  event.yhtofTPCKurama = **src.yhtofTPCKurama;
  event.lvpTPCKurama = **src.lvpTPCKurama;
  event.xvpTPCKurama = **src.xvpTPCKurama;
  event.yvpTPCKurama = **src.yvpTPCKurama;

  event.isgoodTPC = **src.isgoodTPC;
  event.insideTPC = **src.insideTPC;
  event.vtxTPC = **src.vtxTPC;
  event.vtyTPC = **src.vtyTPC;
  event.vtzTPC = **src.vtzTPC;
  event.closeDistTPC = **src.closeDistTPC;
  event.MissMassTPC = **src.MissMassTPC;
  event.MissMassCorrTPC = **src.MissMassCorrTPC;
  event.MissMassCorrDETPC = **src.MissMassCorrDETPC;
  event.MissMassNuclTPC = **src.MissMassNuclTPC;
  event.MissMassNuclCorrTPC = **src.MissMassNuclCorrTPC;
  event.MissMassNuclCorrDETPC = **src.MissMassNuclCorrDETPC;

  event.pOrg = **src.pOrg;
  event.pCalc = **src.pCalc;
  event.pCorr = **src.pCorr;
  event.pCorrDE = **src.pCorrDE;
  event.ub = **src.ub;
  event.vb = **src.vb;
  event.us = **src.us;
  event.vs = **src.vs;

  event.pOrgTPC = **src.pOrgTPC;
  event.pCalcTPC = **src.pCalcTPC;
  event.pCorrTPC = **src.pCorrTPC;
  event.pCorrDETPC = **src.pCorrDETPC;
  event.thetaTPC = **src.thetaTPC;
  event.thetaCMTPC = **src.thetaCMTPC;
  event.costCMTPC = **src.costCMTPC;

  event.BE.resize(event.nKK);
  event.BETPC.resize(event.nKK);
  event.BE_LL.resize(event.nKK);
  event.BETPC_LL.resize(event.nKK);

  Int_t psfac = 0;  
  bool trigA = (event.trigflag[20]>0);
  bool trigB = (event.trigflag[21]>0);  
  if(!trigA&&!trigB) return false;
  if(trigA) psfac = psTrigA;
  if(trigB) psfac = psTrigB;

  for(int ips=0; ips<psfac; ips++){
    HF1( 1, event.status ); // debug 0    
  }
  event.status++;
    
  if(event.nKK != 1) return true;
  double BE = 0.;
  double thetaTPC = 0.;  
  for(Int_t iKK=0; iKK<event.nKK; iKK++){
    //BE = event.MissMassNuclCorrDETPC[iKK] - KaonMass - Boron11Mass - 0.075;
    BE = event.MissMassNuclCorrDETPC[iKK] - KaonMass - Boron11Mass;
    thetaTPC = event.thetaTPC[0];
  }
  if( !(event.runnum >= 5641 && event.runnum <= 5666) ){// not CH2
    if(!(thetaTPC>minThetaKP && thetaTPC<maxThetaKP)) return true;
  }
  if( event.isgoodTPCKurama.size()!=1 ) return false;
  if( event.isgoodTPCKurama[0]!=1 ) return false;
  if( event.insideTPC[0] != 1) return false;
  
  //  if(src.chisqrKurama[0] > MaxChisqrKurama || src.chisqrK18[0] > MaxChisqrBcOut) return true;
  // if(KKEvent && event.Kflag[0] != 1){
  //   if(event.kflagTPCKurama[0]!=1) return false;
  //   return true; //precut with Kurama tracking
  // }
  if(KPEvent && event.Pflag[0] != 1){
    if(event.pflagTPCKurama[0]!=1) return false;    
    return true; //precut with Kurama tracking
  }
  // if(KHeavyEvent && event.Heavyflag[0] != 1){
  //   return true; //precut with Kurama tracking
  // }

  // TVector3 kkvtxTPC(event.vtxTPC[0], event.vtyTPC[0], event.vtzTPC[0] + tpc::ZTarget);
  // event.BE[0] = 1000.*binding_energy; //MeV/c2
  // Double_t binding_energy_LL = m10Be + 2.*LambdaMass - mm_12C; //GeV/c2
  // event.BE_LL[0] = 1000.*binding_energy_LL; //MeV/c2

  TLorentzVector LvRcTPC;
  TVector3 km_unit = TVector3(event.utgtK18[0], event.vtgtK18[0], 1.).Unit();
  TVector3 km_momTPC = km_unit*event.pK18[0];

  TVector3 kp_unit = TVector3(event.usTPC[0], event.vsTPC[0], 1.).Unit();
  TVector3 kp_momTPC = kp_unit*event.pCorrDETPC[0];

  TVector3 miss_momTPC = km_momTPC - kp_momTPC;
  
  //Double_t thetaTPC = event.thetaTPC[0];

  TLorentzVector LvKmTPC(km_momTPC, TMath::Hypot(km_momTPC.Mag(), KaonMass));
  TLorentzVector LvScatPTPC(kp_momTPC, TMath::Hypot(kp_momTPC.Mag(), ProtonMass));
  TLorentzVector LvCTPC(0., 0., 0., Carbon12Mass);
  TLorentzVector LvPTPC(0., 0., 0., ProtonMass);
  //std::cout << "m12C:" << m12C << " ProtonMass:" << ProtonMass << std::endl;
  TLorentzVector LvScatKmTPC = LvKmTPC + LvPTPC - LvScatPTPC;
  LvRcTPC = LvKmTPC + LvCTPC - LvScatPTPC;

  double mm_12CTPC = LvRcTPC.M();
  //double binding_energyTPC = Boron11Mass + KaonMass - (mm_12CTPC - 0.120); //GeV/c2
  double binding_energyTPC = Boron11Mass + KaonMass - mm_12CTPC; //GeV/c2
  event.BETPC[0] = binding_energyTPC; //MeV/c2
  
  //event.MissMassCorrDETPC[0] = event.MissMassCorrDETPC[0] - 0.120 ;
  event.MissMassCorrDETPC[0] = event.MissMassCorrDETPC[0];  
  double missmass = event.MissMassCorrDETPC[0];
  
  if(missmass>minMM&&missmass<maxMM){
    for(int ips=0; ips<psfac; ips++) HF1( 1, event.status ); // debug 1
  }
  event.status++;
  
  for(int ips=0; ips<psfac; ips++){
    HF1(3900,event.MissMassCorrDETPC[0]);
    HF1(13900,-event.BETPC[0]);
  }

  // for Geant4
  Double_t pKp = event.pCorrDETPC[0];
  Double_t uKp = event.utgtTPCKurama[0];  
  Double_t vKp = event.vtgtTPCKurama[0];
  Double_t ptKp = pKp/std::sqrt(1.+uKp*uKp+vKp*vKp);
  TVector3 smom(ptKp*uKp, ptKp*vKp, ptKp);

  Double_t pKm = event.pK18[0];
  Double_t uKm = event.utgtK18[0];
  Double_t vKm = event.vtgtK18[0];
  Double_t ptKm = pKm/std::sqrt(1.+uKm*uKm+vKm*vKm);
  TVector3 bmom(ptKm*uKm, ptKm*vKm, ptKm);  
  dstg4.ich = 0;
  dstg4.np = 1; // only K-
  dstg4.pidpdg[0] = -321;
  dstg4.bpx = bmom[0];
  dstg4.bpy = bmom[1];
  dstg4.bpz = bmom[2];
  dstg4.kppx = smom[0];
  dstg4.kppy = smom[1];
  dstg4.kppz = smom[2];
  TVector3 missmom = bmom - smom;
  dstg4.px[0] = missmom.X();
  dstg4.py[0] = missmom.Y();
  dstg4.pz[0] = missmom.Z();
  // std::cout << " debug: " << __FILE__ << " " << __LINE__ << " "
  // 	    << " bmom(" << bmom[0] << "," << bmom[1] << "," << bmom[2] << ") "
  // 	    << " smom(" << smom[0] << "," << smom[1] << "," << smom[2] << ") "
  // 	    << " px[0],py[0],pz[0]:" << dstg4.px[0] << " " << dstg4.py[0] << " " << dstg4.pz[0] << std::endl;
  // std::cout << " debug: " << __FILE__ << " " << __LINE__ << " "
  // 	    << " LvScatKmTPC.M():" << LvScatKmTPC.M() << " LvScatKmTPC.P():" << LvScatKmTPC.P()
  // 	    << " LvRcTPC.M():" << LvRcTPC.M() << " LvRcTPC.P():" << LvRcTPC.P() << std::endl;
      

  int ntTpc = **src.ntTpc;  
  if( ntTpc == 0 )
    return true;
  if(missmass>minMM&&missmass<maxMM){
    for(int ips=0; ips<psfac; ips++) HF1( 1, event.status ); // debug 2
  }
  event.status++;
  
  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.chisqr = **src.chisqr;
  event.helix_cx = **src.helix_cx;
  event.helix_cy = **src.helix_cy;
  event.helix_z0 = **src.helix_z0;
  event.helix_r = **src.helix_r;
  event.helix_dz = **src.helix_dz;
  event.dE = **src.dE;
  event.dEdx = **src.dEdx;
  event.mom0 = **src.mom0;
  event.charge = **src.charge;
  event.path = **src.path;
  event.pid = **src.pid;
  event.isElectron.resize(ntTpc);
  event.nsigma_triton.resize(ntTpc);
  event.nsigma_deutron.resize(ntTpc);
  event.nsigma_proton.resize(ntTpc);
  event.nsigma_kaon.resize(ntTpc);
  event.nsigma_pion.resize(ntTpc);
  event.nsigma_electron.resize(ntTpc);
  for(int it=0; it<ntTpc; ++it){
    event.isElectron[it] = Kinematics::HypTPCdEdxElectron(event.dEdx[it], event.mom0[it]);
    event.nsigma_triton[it] = Kinematics::HypTPCdEdxNsigmaTriton(event.dEdx[it], event.mom0[it]);
    event.nsigma_deutron[it] = Kinematics::HypTPCdEdxNsigmaDeutron(event.dEdx[it], event.mom0[it]);
    event.nsigma_proton[it] = Kinematics::HypTPCdEdxNsigmaProton(event.dEdx[it], event.mom0[it]);
    event.nsigma_kaon[it]  = Kinematics::HypTPCdEdxNsigmaKaon(event.dEdx[it], event.mom0[it]);
    event.nsigma_pion[it] = Kinematics::HypTPCdEdxNsigmaPion(event.dEdx[it], event.mom0[it]);
    event.nsigma_electron[it] = Kinematics::HypTPCdEdxNsigmaElectron(event.dEdx[it], event.mom0[it]);
  };

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);  
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
			 **src.charge, **src.nhtrack, **src.helix_cx,
			 **src.helix_cy, **src.helix_z0, **src.helix_r,
			 **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
			 **src.helix_t, **src.track_cluster_de, **src.resolution_x,
			 **src.resolution_y, **src.resolution_z, **src.hitpos_x,
			 **src.hitpos_y, **src.hitpos_z);

  
  HypTPCTask& GFtrackCont = HypTPCTask::GetInstance();
  
  for(int it=0; it<event.ntTpc; ++it){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    if(KPEvent&&event.isKurama[it]==1) GFtrackCont.AddHelixTrack(2212, tp);
    else if(KKEvent&&event.isKurama[it]==1) GFtrackCont.AddHelixTrack(321, tp);
    else if(event.isElectron[it]==1) GFtrackCont.AddHelixTrack(event.charge[it]*(-11), tp);
    else{    
      event.pid[it] = tp ->GetPid();
      std::vector<Int_t> pdgcode;
      Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
      if((event.pid[it]&4)!=4 && TMath::Abs(event.nsigma_electron[it]) < 3)
	pdgcode.push_back(event.charge[it]*(-11));	
      GFtrackCont.AddHelixTrack(pdgcode, tp);
    }
  }
  GFtrackCont.FitTracks();
  if(missmass>minMM&&missmass<maxMM){
    for(int ips=0; ips<psfac; ips++) HF1( 1, event.status ); // debug 3
  }
  event.status++;
  
  HF1( 2, event.GFstatus++ );
  
  int GFntTpc = GFtrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    std::cout << "ntTpc:" << ntTpc << " GFntTpc:" << GFntTpc << std::endl;
    return true;
  }
  if(missmass>minMM&&missmass<maxMM){
    for(int ips=0; ips<psfac; ips++) HF1( 1, event.status ); // debug 4
  }
  event.status++;
  
  HF1( 2, event.GFstatus++ );  
  HF1( genfitHid, GFntTpc);
  event.GFntTpc = GFntTpc;
  event.GFfitstatus.resize(GFntTpc);    
  event.GFcharge.resize(GFntTpc);
  event.GFchisqr.resize(GFntTpc);
  event.GFtof.resize(GFntTpc);
  event.GFpval.resize(GFntTpc);
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
  event.GFfromVtx.resize(GFntTpc);  

  event.GFextrapolationHtof.resize(GFntTpc);  
  event.GFtracklen.resize(GFntTpc);
  event.GFtrack2vtxdist.resize(GFntTpc);
  event.GFcalctof.resize(GFntTpc);
  event.GFsegHtof.resize(GFntTpc);
  event.GFtofHtof.resize(GFntTpc);
  event.GFtdiffHtof.resize(GFntTpc);
  event.GFposHtof.resize(GFntTpc);
  event.GFposx.resize(GFntTpc);
  event.GFposy.resize(GFntTpc);
  event.GFposz.resize(GFntTpc);
  event.GFinvbeta.resize(GFntTpc);
  event.GFm2.resize(GFntTpc);
  event.nsigma_tritonHtof.resize(ntTpc);
  event.nsigma_deutronHtof.resize(ntTpc);
  event.nsigma_protonHtof.resize(ntTpc);
  event.nsigma_kaonHtof.resize(ntTpc);
  event.nsigma_pionHtof.resize(ntTpc);
  event.nsigma_electronHtof.resize(ntTpc);
  
  Int_t ntrack_intarget = 0;
  Double_t x0[100] = {0};
  Double_t y0[100] = {0};
  Double_t u0[100] = {0};
  Double_t v0[100] = {0};
  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    event.GFfitstatus[igf] = (int)GFtrackCont.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
    if(GFtrackCont.TrackCheck(igf)) {
      int nh = GFtrackCont.GetNHits(igf);
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

      event.GFchisqr[igf] = GFtrackCont.GetChi2NDF(igf);
      event.GFcharge[igf] = GFtrackCont.GetCharge(igf);
      event.GFtof[igf] = GFtrackCont.GetTrackTOF(igf, 0, -1);
      event.GFpval[igf] = GFtrackCont.GetPvalue(igf);
      event.GFnhtrack[igf] = GFtrackCont.GetNHits(igf);
      event.GFpdgcode[igf] = GFtrackCont.GetPDGcode(igf);

      HF1( genfitHid+1, event.GFchisqr[igf]);
      HF1( genfitHid+2, event.GFpval[igf]);
      HF1( genfitHid+3, event.GFcharge[igf]);
      HF1( genfitHid+4, event.GFnhtrack[igf]);
      HF1( genfitHid+5, event.GFtracklen[igf]);
      HF1( genfitHid+6, event.GFtof[igf]);
      for( Int_t ihit=0; ihit<nh; ++ihit ){
	TVector3 hit = GFtrackCont.GetPos(igf, ihit);
	TVector3 mom = GFtrackCont.GetMom(igf, ihit);
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
      if( event.isBeam[igf]==1 || event.isK18[igf]==1 || event.isAccidental[igf]==1 ) continue;      
      if(event.isKurama[igf]==1){
	for(int jgf=0; jgf<GFntTpc; jgf++){
	  if(igf==jgf) continue;
	  if( event.charge[jgf]==1 ) continue;
	  if( event.isBeam[jgf]==1 || event.isK18[jgf]==1 || event.isAccidental[jgf]==1 ) continue;
	  int repidKurama = 0;	  
	  Int_t repid = 0;
	  if((event.pid[jgf]&1)==1) repid = 0;
	  else if((event.pid[jgf]&2)==2) {
	    if((1&event.pid[jgf])==1) repid += 1;
	  }
	  else if((event.pid[jgf]&4)==4){
	    Int_t flag = 1;
	    for(Int_t i=0;i<2;i++){
	      Int_t temp = flag&event.pid[jgf];
	      if(temp==flag) repid += 1;
	      flag*=2;
	    }
	  }	  
	  double thetafac = 0.3;
	  Double_t kurama_par[5];
	  kurama_par[0] = event.helix_cx[igf];
	  kurama_par[1] = event.helix_cy[igf];
	  kurama_par[2] = event.helix_z0[igf];
	  kurama_par[3] = event.helix_r[igf];
	  kurama_par[4] = event.helix_dz[igf];
	  Int_t kurama_nh = event.helix_t[igf].size();
	  Double_t kurama_theta_range = event.helix_t[igf][kurama_nh-1] - event.helix_t[igf][0];
	  Double_t kurama_theta_min = event.helix_t[igf][0] - thetafac*kurama_theta_range;
	  Double_t kurama_theta_max = event.helix_t[igf][kurama_nh-1] + thetafac*kurama_theta_range;
	  
	  Double_t other_par[5];
	  other_par[0] = event.helix_cx[jgf];
	  other_par[1] = event.helix_cy[jgf];
	  other_par[2] = event.helix_z0[jgf];
	  other_par[3] = event.helix_r[jgf];
	  other_par[4] = event.helix_dz[jgf];
	  Int_t other_nh = event.helix_t[jgf].size();
	  Double_t other_theta_range = event.helix_t[jgf][0] - event.helix_t[jgf][other_nh-1];
	  Double_t other_theta_min = event.helix_t[jgf][other_nh-1] - thetafac*other_theta_range;
	  Double_t other_theta_max = event.helix_t[jgf][0] + thetafac*other_theta_range;
	  double thetaKurama,thetaOther;
	  double dist = qnan;
	  TVector3 vertex;	  
	  vertex = Kinematics::VertexPointHelix(kurama_par,other_par,
				    kurama_theta_min,kurama_theta_max,
				    other_theta_min,other_theta_max,
				    thetaKurama,thetaOther,
				    dist);				    
	  // if( !(GFtrackCont.FindVertex(igf, jgf,
	  // 			       repidKurama, repid,
	  // 			       extrapKurama,extrapOther,
	  // 			       momKurama,momOther,
	  // 			       dist, vertex,
	  // 			       vtx_scan_range) ) )  continue;
	  Bool_t vtxouttgt
	    = dist < ppi_distcut
	    && ( TMath::Abs(vertex.x()) > 25. 
		 || TMath::Abs(vertex.y()) > 25. 
		 || TMath::Abs(vertex.z() - tpc::ZTarget) > 50.) 
	    && ( vertex.z() - tpc::ZTarget>0 );
	  //std::cout << " VertexZ: " << vertex.z() << " VertexZ-tgtZ: " << vertex.z() - tpc::ZTarget << std::endl;
	  event.GFKuramaVtxOutTgt = vtxouttgt;
	  if(vtxouttgt){
	    event.GFFwdYDecayTrack = jgf;	    
	    //std::cout << " ForwardYDecayLikeTrack: " << jgf << std::endl;
	  }
	  if(vtxouttgt) break;
	}
      }
      if(GFtrackCont.IsInsideTarget(igf)){
	event.GFinside[igf] = 1;
	if(event.isKurama[igf]){
	  event.GFinsideKurama = 1;
	}
	TVector3 post; TVector3 momt; double lentgt; double toftgt;
	if(GFtrackCont.ExtrapolateToTargetCenter(igf, post, momt, lentgt, toftgt)){
	  x0[ntrack_intarget] = post.x();
	  y0[ntrack_intarget] = post.y();
	  u0[ntrack_intarget] = momt.x()/momt.z();
	  v0[ntrack_intarget] = momt.y()/momt.z();
	  ntrack_intarget++;
	}
	TVector3 vertex = Kinematics::MultitrackVertex(ntrack_intarget, x0, y0, u0, v0);
	event.GFntTpc_inside = ntrack_intarget;	
	event.GFprodvtx_x = vertex.x();
	event.GFprodvtx_y = vertex.y();
	event.GFprodvtx_z = vertex.z();    
	//TVector3 vertex(event.vtxTPC[igf],event.vtyTPC[igf],event.vtzTPC[igf]+tpc::ZTarget);
	Double_t mom = event.GFmom[igf][0];
	Int_t repid=-1;
	Int_t hitid_htof; Double_t tof; Double_t len;
	TVector3 pos_htof; Double_t track2tgt_dist;
	Bool_t htofextrapvtx =
	  GFtrackCont.TPCHTOFTrackMatching(igf, repid, vertex,
					   event.HtofSeg, event.posHtof,
					   hitid_htof, tof,
					   len, pos_htof, track2tgt_dist);
	if(htofextrapvtx){
	  event.GFfromVtx[igf] = 1;
	  event.GFtracklen[igf] = len;
	  //event.GFtrack2vtxdist[igf] = track2tgt_dist;
	  event.GFcalctof[igf] = tof;
	  event.GFposx[igf] = pos_htof.x();
	  event.GFposy[igf] = pos_htof.y();
	  event.GFposz[igf] = pos_htof.z();
	  event.GFsegHtof[igf] = event.HtofSeg[hitid_htof];
	  event.GFtofHtof[igf] = event.tHtof[hitid_htof];
	  event.GFposHtof[igf] = event.posHtof[hitid_htof];
	  Double_t beta = len/event.tHtof[hitid_htof]/MathTools::C();	  
	  event.GFinvbeta[igf] = 1./beta;
	  Double_t mass2 = Kinematics::MassSquare(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	  event.GFm2[igf] = mass2;
	  event.nsigma_tritonHtof[igf] = Kinematics::HypTPCHTOFNsigmaTriton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	  event.nsigma_deutronHtof[igf] = Kinematics::HypTPCHTOFNsigmaDeutron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);                                                                          
	  event.nsigma_protonHtof[igf] = Kinematics::HypTPCHTOFNsigmaProton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);                                                                            
	  event.nsigma_kaonHtof[igf] = Kinematics::HypTPCHTOFNsigmaKaon(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	  event.nsigma_pionHtof[igf] = Kinematics::HypTPCHTOFNsigmaPion(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
	  event.nsigma_electronHtof[igf] = Kinematics::HypTPCHTOFNsigmaElectron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]); 	  
	}	
      } else {
	event.GFinside[igf] = 0;
      }
    } else {
      event.GFnhtrack[igf] = 0;
      event.GFpdgcode[igf] = -9999;
        
      event.GFchisqr[igf] = TMath::QuietNaN();
      event.GFcharge[igf] = TMath::QuietNaN();
      event.GFtof[igf]    = TMath::QuietNaN();
      event.GFpval[igf]   = TMath::QuietNaN();

      event.GFlayer[igf].clear();
      event.GFpos_x[igf].clear();
      event.GFpos_y[igf].clear();
      event.GFpos_z[igf].clear();
      event.GFmom[igf].clear();
      event.GFmom_x[igf].clear();
      event.GFmom_y[igf].clear();
      event.GFmom_z[igf].clear();
      event.GFresidual_x[igf].clear();
      event.GFresidual_y[igf].clear();
      event.GFresidual_z[igf].clear();
      event.GFresidual_p[igf].clear();
      event.GFresidual_px[igf].clear();
      event.GFresidual_py[igf].clear();
      event.GFresidual_pz[igf].clear();

      event.GFinside[igf] = -1;
    }   
  }

  
  std::vector<Int_t> k_id_container;
  std::vector<Int_t> k_repid_container;
  std::vector<TVector3> k_mom_container, k_vert_container;  
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_targetvtx_container;
  std::vector<Int_t>    k_htofextrap_container;
  std::vector<Int_t>    k_htofhitid_container;
  std::vector<Int_t>    k_htofseg_container;
  std::vector<Double_t> k_mass2_container;
  std::vector<Double_t> k_posY_container;
  std::vector<TVector3> k_htofpos_container;
  std::vector<Double_t> k_diffY_container;
  std::vector<Double_t> k_tracklen_container;

  Int_t numK18=0; Int_t numKurama=0; Int_t numBeam=0; Int_t numAcc=0;
  Int_t numKm=0; Int_t numPim=0; Int_t numPip=0; Int_t numP=0; Int_t numPPip=0;
  Int_t numEp=0; Int_t numEm=0;
  Int_t idKm=-1;
  bool kmflag_dedxpid = false;    
  {
    for(Int_t it=0;it<ntTpc;it++){ // kaon
      if( event.isK18[it]==1 ){
	numK18++;
      } else if( event.isKurama[it]==1 ) {
	numKurama++;
      } else if( event.isBeam[it]==1 ) {
	numBeam++;
      } else if( event.isAccidental[it]==1 ) {
	numAcc++;
      }
      if(event.isElectron[it]==1) continue;
      if(event.isK18[it]==1) continue;
      if(event.isKurama[it]==1) continue;
      if(event.isBeam[it]==1) continue;
      if(event.isAccidental[it]==1) continue;
      if(event.GFinside[it]!=1) continue;
      
      //if((event.pid[it]&2)==2 && event.charge[it]==-1){ //k-
      if(event.charge[it]==-1){ //k-
	if( event.nsigma_kaon[it]>mindEdxSigKaon && event.nsigma_kaon[it]<maxdEdxSigKaon ){ 
	  numKm++; 
	  idKm = it; 
	  continue; 
	}
      }
      if(event.isElectron[it]==1){ //e+, e-
	if(event.charge[it]==1){
	  numEp++;
	  continue;
	} else {
	  numEm++;
	  continue;
	}
      }    
      else if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton
	numP++;
	continue;      
      }
      else if((event.pid[it]&1)==1 && event.charge[it]==-1){ //pi-
	Double_t slope = event.helix_dz[it];
	Double_t helixmom = event.mom0[it];
	Int_t pi_nh = event.helix_t[it].size();      
	TVector3 pi_start = TVector3(event.calpos_x[it][0],
				     event.calpos_y[it][0],
				     event.calpos_z[it][0]);
	TVector3 pi_end = TVector3(event.calpos_x[it][pi_nh-1],
				   event.calpos_y[it][pi_nh-1],
				   event.calpos_z[it][pi_nh-1]);
	double pi_vertex_dist=-999.;
	if(TMath::Abs(slope)<0.05 && TMath::Abs(helixmom)>0.5 &&
	   !(Kinematics::HelixDirection(tgtpos,pi_start,pi_end,pi_vertex_dist))&&
	   (pi_start.x()-pi_end.x())>-10 && (pi_start.x()-pi_end.x())<50.){
	  event.isAccidental[it] = 1;
	  numAcc++;
	  //target_accidental_id_container.push_back(it); //Accidental beam on the target
	  continue; //Accidental K-
	}      
	numPim++;
	continue;
      }
      else if((event.pid[it]&4)!=4 && (event.pid[it]&1)==1 && event.charge[it]==1){ //pi+
	numPip++;
	continue;
      }
      else if(((event.pid[it]&4)==4 || (event.pid[it]&1)==1) && event.charge[it]==1){ //p or pi+ with high-mom
	numPPip++;
	continue;
      }
    }
  }

  int l_candidates=1;
  std::vector<Double_t> GFk_targetdist_container(l_candidates, qnan);
  std::vector<Double_t> GFk_tof_container(l_candidates, qnan);
  std::vector<Double_t> GFk_ctau_container(l_candidates, qnan);
  std::vector<TVector3> GFk_targetvtx_container(l_candidates, qnan_vec);
  std::vector<Double_t> GFk_targetcenterdist_container(l_candidates, qnan);
  std::vector<TVector3> GFk_targetcentervtx_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFk_mom_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFk_vert_container(l_candidates, qnan_vec);
  std::vector<Int_t> GFk_id_container(l_candidates, qnan);
  std::vector<Int_t> GFk_repid_container(l_candidates, qnan);
  std::vector<Int_t> GFk_extrap_container(l_candidates, qnan);
  std::vector<Int_t> GFk_htofhitid_container(l_candidates, qnan);
  std::vector<Int_t> GFk_htofseg_container(l_candidates, qnan);
  std::vector<Double_t> GFk_tracklen_container(l_candidates, qnan);
  std::vector<TVector3> GFk_poshtof_container(l_candidates, qnan_vec);
  std::vector<Double_t> GFk_mass2_container(l_candidates, qnan);
  std::vector<Double_t> GFk_invbeta_container(l_candidates, qnan);  

  if(numKm>0){
    // std::cout << " debug " << __FILE__ << " " << __LINE__
    // 	      << " numKm:" << numKm << " numPPip:" << numPPip << " numPip:" << numPip << " numPim:" << numPim << " numP:" << numP
    // 	      << " numEm:" << numEm << " numEp:" << numEp << std::endl;
    event.kminc = true;
  }
  if( numPPip==0&&numPip==0&&numPim==0&&numP==0&&numEm==0&&numEp==0 ){
    if(missmass>minMM&&missmass<maxMM){
      for(int ips=0; ips<psfac; ips++){
	HF1( 1, event.status ); // debug 5
	if(ips==psfac-1) event.status++;	
      }
    }
    for(int it=0; it<ntTpc; it++){
      if(idKm!=it) continue;
      kmflag_dedxpid=true;            
      if ( !event.GFfitstatus[it] ) continue;       
      if ( event.isElectron[it]==1 ) continue; 
      if ( event.isK18[it]==1 ) continue; 
      if ( event.isKurama[it]==1 ) continue;
      if ( event.isBeam[it]==1 ) continue;
      if ( event.isAccidental[it]==1 ) continue;      
      Int_t repid_km = -1;
      //if(!GFtrackCont.TrackCheck(it, repid_km)) continue;
      if(event.GFinside[it]!=1) continue;
      Int_t km_nh = event.helix_t[it].size();
      TVector3 km_start = TVector3(event.calpos_x[it][0], event.calpos_y[it][0], event.calpos_z[it][0]);
      TVector3 km_end = TVector3(event.calpos_x[it][km_nh-1], event.calpos_y[it][km_nh-1], event.calpos_z[it][km_nh-1]);
      double vertex_dist = 0.;
      if(!Kinematics::HelixDirection(tgtpos,km_start,km_end,vertex_dist)) continue;
      TVector3 post; TVector3 km_mom; double lentgt; double toftgt;
      Int_t extrap=-1;
      if(GFtrackCont.ExtrapolateToTargetCenter(it, post, km_mom, lentgt, toftgt)){
	
      }

      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist; Int_t htofseg;
      Bool_t km_htofextrap =
	GFtrackCont.TPCHTOFTrackMatching(it, repid_km, tgtpos,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(km_htofextrap){
	GFk_htofhitid_container[0] = hitid_htof;
	GFk_htofseg_container[0] = event.HtofSeg[hitid_htof];
	GFk_tracklen_container[0] = tracklen_htof;
	GFk_poshtof_container[0] = pos_htof;		
	GFk_tof_container[0] = event.tHtof[hitid_htof];
	GFk_mass2_container[0] =
	  Kinematics::MassSquare(km_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof])+0.05;
	GFk_invbeta_container[0] =
	  MathTools::C()*event.tHtof[hitid_htof]/tracklen_htof;
      }
      GFk_id_container[0]=idKm;
      GFk_repid_container[0]=repid_km;
      GFk_mom_container[0]=km_mom;
      GFk_targetvtx_container[0]=post;
      GFk_targetcentervtx_container[0]=post-tgtpos;
      TVector3 dist = post-tgtpos;
      GFk_targetcenterdist_container[0]=dist.Mag();
      event.kmflag = true;
    }
  }

  if(kmflag_dedxpid){
    for(int ips=0; ips<psfac; ips++){
      HF2(10, event.charge[idKm]*event.mom0[idKm], event.dEdx[idKm]);    
      HF1(3950,event.MissMassCorrDETPC[0]);
      HF1(13950,-event.BETPC[0]);      
    }
  }
  
  if(event.kmflag){
    int id=0;
    event.GFkmid = GFk_id_container[id];  
    event.GFkmmom = GFk_mom_container[id].Mag();
    event.GFkmmom_x = GFk_mom_container[id].x();
    event.GFkmmom_y = GFk_mom_container[id].y();
    event.GFkmmom_z = GFk_mom_container[id].z();

    TVector3 scat_km_mom(event.GFkmmom_x, event.GFkmmom_y, event.GFkmmom_z);
    TVector3 diff_kmmom = scat_km_mom - miss_momTPC;
    event.GFkmtheta = scat_km_mom.Theta(); // rad
    event.GFkmphi   = scat_km_mom.Phi();   // rad
    double GFkmCosTheta = scat_km_mom.CosTheta();
    
    event.GFkmtarget_dist = GFk_targetdist_container[id];
    event.GFkmtargetvtx_x = GFk_targetvtx_container[id].x();
    event.GFkmtargetvtx_y = GFk_targetvtx_container[id].y();
    event.GFkmtargetvtx_z = GFk_targetvtx_container[id].z();
    event.GFkmtargetcenter_dist = GFk_targetcenterdist_container[id];
    event.GFkmtargetcenter_x = GFk_targetcentervtx_container[id].x();
    event.GFkmtargetcenter_y = GFk_targetcentervtx_container[id].y();
    event.GFkmtargetcenter_z = GFk_targetcentervtx_container[id].z();
  
    event.GFkmhtofid = GFk_htofhitid_container[id];
    event.GFkmhtofseg = GFk_htofseg_container[id];
    event.GFkmposHtof.push_back(GFk_poshtof_container[id].x());
    event.GFkmposHtof.push_back(GFk_poshtof_container[id].y());
    event.GFkmposHtof.push_back(GFk_poshtof_container[id].z());
    event.GFkmmass2 = GFk_mass2_container[id];
    event.GFkmtracklen = GFk_tracklen_container[id];
    event.GFkminvbeta = GFk_invbeta_container[id];
    event.GFkmtof = GFk_tof_container[id];
    
    if(event.GFkmmom>0.01){
      for(int ips=0; ips<psfac; ips++){
	HF1(20, event.GFkmmass2);
	for(int i=0; i<5; i++){
	  if(double(i)*0.2<event.GFkmmom && double(i+1)*0.2>event.GFkmmom) HF1(121+i, event.GFkmmass2);
	}
	if(0.2<event.GFkmmom && event.GFkmmom<0.7) HF1(126, event.GFkmmass2);
	if(0.7<event.GFkmmom) HF1(127, event.GFkmmass2);
	
	HF1(21, event.GFkmtracklen);
	HF1(22, event.GFkmtof);
	HF1(23, event.dEdx[event.GFkmid]);      
	HF1(3050, event.GFkmmom);
	HF1(3100, event.GFkmtheta*TMath::RadToDeg());
	HF1(3110, GFkmCosTheta);
	HF1(3120, event.GFkmphi);      
	HF1(3121, event.GFkmphi*TMath::RadToDeg());
	HF2(3500, event.GFkmtheta*TMath::RadToDeg(), event.GFkmmom);
	HF2(3510, event.GFkmphi*TMath::RadToDeg(), event.GFkmmom);
	HF2(3520, event.GFkmphi*TMath::RadToDeg(), event.GFkmtheta*TMath::RadToDeg());
	for(int i=0; i<numbinbek; i++){
	  if(-event.BETPC[0]>minbek+i*onebinbek && -event.BETPC[0]<minbek+(i+1)*onebinbek) HF1(201+i,event.GFkmmass2);
	}
	if(event.GFkmmom>0.01&&event.GFkmmass2>minM2Km&&event.GFkmmass2<maxM2Km){
	  if(event.GFkmmom<1.0){
	    HF1(3951, event.MissMassCorrDETPC[0]);
	    HF1(13951, -event.BETPC[0]);
	    // hist
	    HF1(4050, event.GFkmmom);
	    HF1(4060, diff_kmmom.Mag());
	    HF1(4070, miss_momTPC.Mag());
	    double bek = -event.BETPC[0];
	    for(int i=0; i<8; i++){
	      if( double(i)*0.050-0.100<bek && double(i+1)*0.050-0.100>bek ){
		std::cout << "-bek(8rgn):" << bek <<  " binbek(8rgn):" << i << " diffMissMom:" << diff_kmmom.Mag() << std::endl;
		HF1(4062+i, diff_kmmom.Mag());
	      }
	    }
	    int binbek = GetBinIndexBEk(bek);	
	    HF1(14000+binbek, diff_kmmom.Mag());
	    std::cout << "-bek:" << bek <<  " binbek:" << binbek << " diffMissMom:" << diff_kmmom.Mag() << std::endl;
	    if(missmass>minMM&&missmass<maxMM){
	      HF1(4051, event.GFkmmom);
	      HF1(4061, diff_kmmom.Mag());
	      HF1(4071, miss_momTPC.Mag());	      
	    }
	    HF1(4100, event.GFkmtheta*TMath::RadToDeg());
	    HF1(4110, GFkmCosTheta);
	    HF1(4120, event.GFkmphi);      
	    HF1(4121, event.GFkmphi*TMath::RadToDeg());
	    HF2(4500, event.GFkmtheta*TMath::RadToDeg(), event.GFkmmom);
	    HF2(4510, event.GFkmphi*TMath::RadToDeg(), event.GFkmmom);
	    HF2(4520, event.GFkmphi*TMath::RadToDeg(), event.GFkmtheta*TMath::RadToDeg());
	  }
	}
      }
    }    
    double pmom = kp_momTPC.Mag();
    double ptheta = kp_momTPC.Theta(); //rad
    double pcost = kp_momTPC.CosTheta();
    double pphi = kp_momTPC.Phi(); //rad
    for(int ips=0; ips<psfac; ips++){    
      HF1(1050, pmom);
      HF1(1100, ptheta*TMath::RadToDeg());
      HF1(1110, pcost);
      HF1(1120, pphi);      
      HF1(1121, pphi*TMath::RadToDeg());
      HF2(1500, ptheta*TMath::RadToDeg(), pmom);
      HF2(1510, pphi*TMath::RadToDeg(), pmom);
      HF2(1520, pphi*TMath::RadToDeg(), ptheta*TMath::RadToDeg());
      if(event.GFkmmom>0.01&&event.GFkmmom<1.0&&event.GFkmmass2>minM2Km&&event.GFkmmass2<maxM2Km){
	HF1(2050, pmom);
	HF1(2100, ptheta*TMath::RadToDeg());
	HF1(2110, pcost);
	HF1(2120, pphi);      
	HF1(2121, pphi*TMath::RadToDeg());
	HF2(2500, ptheta*TMath::RadToDeg(), pmom);
	HF2(2510, pphi*TMath::RadToDeg(), pmom);
	HF2(2520, pphi*TMath::RadToDeg(), ptheta*TMath::RadToDeg());
      }
    }    
  }
  HF1( 2, event.GFstatus++);
  //HF1( 1, event.status++ ); // debug 9
  GFtrackCont.Clear();
  
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
  Int_t nbinpoq = 1000;
  Int_t minpoq = -1.5;
  Int_t maxpoq = 1.5;
  Int_t nbindedx = 1000;
  Int_t mindedx = 0;
  Int_t maxdedx = 350;
  Int_t nbinmass2 = 100;
  Int_t minmass2 = -1;
  Int_t maxmass2 = 1;
  Int_t nbininvbeta = 100;
  Int_t mininvbeta = 0;
  Int_t maxinvbeta = 5;
  
  static const auto KPEvent = gUser.GetParameter("KPEvent");
  static const auto KKEvent = gUser.GetParameter("KKEvent");
  HB1(1, "Status", 21, 0., 21. );
  HB2(10, "AnalysisKm <dE/dx>;p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",nbinpoq,minpoq,maxpoq,nbindedx,mindedx,maxdedx);
  HB1(20, "Kaon M2; MassSquare [GeV]; counts ", nbinmass2, minmass2, maxmass2);
  HB1(21, "Kaon tracklen; tracklen [mm]; counts ", 1000, 0, 1000);
  HB1(22, "Kaon tof; tof [nsec]; counts ", 500, 0, 10);
  HB1(23, "Kaon dEdx; dEdx [arb.unit]; counts ", 1000, 0, 350);
  HB1(50, "Mom of beam K- [total] ; momenutm [GeV/c]; counts", 500, 1.50, 2.0);

  for(int i=0; i<5; i++){
    HB1(121+i, Form("Kaon M2 (%f<mom<%f [GeV/c]); MassSquare [GeV]; counts ",double(i)*0.2, double(i+1)*0.2), nbinmass2, minmass2, maxmass2);
  }
  HB1(126, "Kaon M2 (0.2<mom<0.7 [GeV/c]); MassSquare [GeV]; counts ", nbinmass2, minmass2, maxmass2);
  HB1(127, "Kaon M2 (0.7<mom<1.0 [GeV/c]); MassSquare [GeV]; counts ", nbinmass2, minmass2, maxmass2);

  for(int i=0; i<numbinbek; i++){
    HB1(201+i, Form("M2 (%.2f<-BEk<%.2f [GeV]); MassSquare [GeV]; counts ",double(i)*onebinbek+minbek, double(i+1)*onebinbek+minbek), nbinmass2, minmass2, maxmass2);
  }  
  
  HB1(1050, " Mom of scatP ; momentum [GeV/c]; coutns", 500, 1.50, 2.0);
  HB1(1100, " Theta of scatP [deg]; #theta [deg]; counts", 300, 0, 30); 
  HB1(1110, " CosTheta of scatP ; Cos(#theta); counts", 200, -1, 1);
  HB1(1120, " Phi of scatP [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(1121, " Phi of scatP [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(1500, " Mom vs Theta scatP; #theta [deg]; momentum [GeV/c]", 300, 0, 30, 500, 1.5, 2.0);
  HB2(1510, " Mom vs Phi scatP; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 500, 1.5, 2.0);
  HB2(1520, " Theta vs Phi scatP; #phi [deg]; #theta [deg]", 3600, -180, 180, 300, 0, 30);
  
  HB1(2050, "[M2] Mom of scatP ; momentum [GeV/c]; coutns", 500, 1.5, 2.0);
  HB1(2100, "[M2] Theta of scatP [deg]; #theta [deg]; counts", 300, 0, 30); 
  HB1(2110, "[M2] CosTheta of scatP ; Cos(#theta); counts", 200, -1, 1);
  HB1(2120, "[M2] Phi of scatP [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(2121, "[M2] Phi of scatP [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(2500, "[M2] Mom vs Theta scatP; #theta [deg]; momentum [GeV/c]", 300, 0, 30, 500, 1.5, 2.0);
  HB2(2510, "[M2] Mom vs Phi scatP; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 500, 1.5, 2.0);
  HB2(2520, "[M2] Theta vs Phi scatP; #phi [deg]; #theta [deg]", 3600, -180, 180, 300, 0, 30);  

  HB1(3050, " Mom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  HB1(3100, " Theta of scatK- [deg]; #theta [deg]; counts", 1800, 0, 180); 
  HB1(3110, " CosTheta of scatK- ; Cos(#theta); counts", 200, -1, 1);
  HB1(3120, " Phi of scatK- [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(3121, " Phi of scatK- [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(3500, " Mom vs Theta scatK-; #theta [deg]; momentum [GeV/c]", 1800, 0, 180, 500, 0, 1.0);
  HB2(3510, " Mom vs Phi scatK-; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 500, 0, 1.0);
  HB2(3520, " Theta vs Phi scatK-; #phi [deg]; #theta [deg]", 3600, -180, 180, 1800, 0, 180);
  HB1(3900, " MissingMass KP inclusive; MissingMass [GeV]; Counts", 280, 0., 1.4);
  HB1(3950, " MissingMass KP exclusive [dEdxPID]; MissingMass [GeV]; Counts", 280, 0., 1.4);
  HB1(3951, " MissingMass KP exclusive [M2PID]; MissingMass [GeV]; Counts", 280, 0., 1.4);

  HB1(4050, "[M2] Mom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0); 
  HB1(4051, "[M2][missmass] Mom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0); 
  HB1(4060, "[M2] (MissMom - Mom) of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0); 
  HB1(4061, "[M2][missmass] (MissMom - Mom) of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);
  for(int i=0; i<8; i++){
    HB1(4062+i, Form("[M2][bek] (MissMom - Mom) of scatK- (%.2f<-BE_K<%.2f); momentum [GeV/c]; coutns",double(i)*0.050-0.100, double(i+1)*0.050-0.100), 500, 0., 1.0);
  }
  HB1(4070, "[M2] MissMom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0); 
  HB1(4071, "[M2][missmass] MissMom of scatK- ; momentum [GeV/c]; coutns", 500, 0., 1.0);

  HB1(4100, "[M2] Theta of scatK- [deg]; #theta [deg]; counts", 1800, 0, 180); 
  HB1(4110, "[M2] CosTheta of scatK- ; Cos(#theta); counts", 200, -1, 1);
  HB1(4120, "[M2] Phi of scatK- [rad]; #phi [rad]; counts", 500, -TMath::Pi(), TMath::Pi());
  HB1(4121, "[M2] Phi of scatK- [deg]; #phi [deg]; counts", 3600, -180, 180);
  HB2(4500, "[M2] Mom vs Theta scatK-; #theta [deg]; momentum [GeV/c]", 1800, 0, 180, 500, 0, 1.0);
  HB2(4510, "[M2] Mom vs Phi scatK-; #phi [deg]; momentum [GeV/c]", 3600, -180, 180, 500, 0, 1.0);
  HB2(4520, "[M2] Theta vs Phi scatK-; #phi [deg]; #theta [deg]", 3600, -180, 180, 1800, 0, 180);

  HB1(13900, " BE KP inclusive; MissingMass [GeV]; Counts", 120, -0.3, 0.3);
  HB1(13950, " BE KP exclusive [dEdxPID]; MissingMass [GeV]; Counts", 120, -0.3, 0.3);
  HB1(13951, " BE KP exclusive [M2PID]; MissingMass [GeV]; Counts", 120, -0.3, 0.3);

  for(int i=0; i<120; i++){
    HB1(14000+i, Form("[M2][bek] (MissMom - Mom) of scatK- (%.3f<-BE_K<%.3f); momentum [GeV/c]; coutns", double(i)*0.005-0.300, double(i+1)*0.005-0.300), 500, 0., 1.0);
  }  
  
  HBTree( "tree", "tree of GenfitQFKaon" );
  tree->Branch( "ich", &dstg4.ich );
  tree->Branch( "bpx", &dstg4.bpx );
  tree->Branch( "bpy", &dstg4.bpy );
  tree->Branch( "bpz", &dstg4.bpz );
  tree->Branch( "kppx", &dstg4.kppx );
  tree->Branch( "kppy", &dstg4.kppy );
  tree->Branch( "kppz", &dstg4.kppz );
  tree->Branch( "np", &dstg4.np );
  tree->Branch( "pid", &dstg4.pidpdg, Form("pid[%d]/I",MaxHits) );
  tree->Branch( "px", &dstg4.px, Form("px[%d]/D",MaxHits) );
  tree->Branch( "py", &dstg4.py, Form("py[%d]/D",MaxHits) );      
  tree->Branch( "pz", &dstg4.pz, Form("pz[%d]/D",MaxHits) );  

  treetpc = new TTree("tpc","Data Summary Table of GenfitQFKaon");    
  treetpc->Branch( "status", &event.status );
  treetpc->Branch( "runnum", &event.runnum );
  treetpc->Branch( "evnum", &event.evnum );
  treetpc->Branch( "trigpat", &event.trigpat );
  treetpc->Branch( "trigflag", &event.trigflag );
  treetpc->Branch( "nhHtof", &event.nhHtof );
  treetpc->Branch( "HtofSeg", &event.HtofSeg );
  treetpc->Branch( "tHtof", &event.tHtof );
  treetpc->Branch( "dtHtof", &event.dtHtof );
  treetpc->Branch( "deHtof", &event.deHtof );
  treetpc->Branch( "posHtof", &event.posHtof );
  treetpc->Branch( "nclTpc", &event.nclTpc );
  treetpc->Branch( "cluster_x", &event.cluster_x );
  treetpc->Branch( "cluster_y", &event.cluster_y );
  treetpc->Branch( "cluster_z", &event.cluster_z );
  treetpc->Branch( "cluster_de", &event.cluster_de );
  treetpc->Branch( "cluster_size", &event.cluster_size );
  treetpc->Branch( "cluster_layer", &event.cluster_layer );
  treetpc->Branch( "cluster_row_center", &event.cluster_row_center );
  treetpc->Branch( "cluster_mrow", &event.cluster_mrow );
  treetpc->Branch( "cluster_de_center", &event.cluster_de_center );
  treetpc->Branch( "cluster_x_center", &event.cluster_x_center );
  treetpc->Branch( "cluster_y_center", &event.cluster_y_center );
  treetpc->Branch( "cluster_z_center", &event.cluster_z_center );
  treetpc->Branch( "cluster_houghflag", &event.cluster_houghflag );
  treetpc->Branch( "ntTpc", &event.ntTpc );
  treetpc->Branch( "nhtrack", &event.nhtrack );
  treetpc->Branch( "isBeam", &event.isBeam );
  treetpc->Branch( "isK18", &event.isK18 );
  treetpc->Branch( "isKurama", &event.isKurama );
  treetpc->Branch( "isAccidental", &event.isAccidental );
  treetpc->Branch( "charge", &event.charge );
  treetpc->Branch( "pid", &event.pid );  
  treetpc->Branch( "chisqr", &event.chisqr );
  treetpc->Branch( "helix_cx", &event.helix_cx );
  treetpc->Branch( "helix_cy", &event.helix_cy );
  treetpc->Branch( "helix_z0", &event.helix_z0 );
  treetpc->Branch( "helix_r", &event.helix_r );
  treetpc->Branch( "helix_dz", &event.helix_dz );
  treetpc->Branch( "mom0", &event.mom0 );
  treetpc->Branch( "path", &event.path );  
  treetpc->Branch( "dE", &event.dE );
  treetpc->Branch( "dEdx", &event.dEdx );  
  treetpc->Branch( "isElectron", &event.isElectron );
  treetpc->Branch( "nsigma_triton", &event.nsigma_triton );
  treetpc->Branch( "nsigma_deutron", &event.nsigma_deutron );
  treetpc->Branch( "nsigma_proton", &event.nsigma_proton );
  treetpc->Branch( "nsigma_kaon", &event.nsigma_kaon );
  treetpc->Branch( "nsigma_pion", &event.nsigma_pion );
  treetpc->Branch( "nsigma_electron", &event.nsigma_electron );
  treetpc->Branch( "hitlayer", &event.hitlayer );
  treetpc->Branch( "hitpos_x", &event.hitpos_x );
  treetpc->Branch( "hitpos_y", &event.hitpos_y );
  treetpc->Branch( "hitpos_z", &event.hitpos_z );
  treetpc->Branch( "calpos_x", &event.calpos_x );
  treetpc->Branch( "calpos_y", &event.calpos_y );
  treetpc->Branch( "calpos_z", &event.calpos_z );
  treetpc->Branch( "mom_x", &event.mom_x );
  treetpc->Branch( "mom_y", &event.mom_y );
  treetpc->Branch( "mom_z", &event.mom_z );
  treetpc->Branch( "residual", &event.residual );
  treetpc->Branch( "residual_x", &event.residual_x );
  treetpc->Branch( "residual_y", &event.residual_y );
  treetpc->Branch( "residual_z", &event.residual_z );
  treetpc->Branch( "resolution_x", &event.resolution_x );
  treetpc->Branch( "resolution_y", &event.resolution_y );
  treetpc->Branch( "resolution_z", &event.resolution_z );
  treetpc->Branch( "helix_t", &event.helix_t );
  treetpc->Branch( "alpha", &event.alpha);
  treetpc->Branch( "pathhit", &event.pathhit);
  treetpc->Branch( "track_cluster_de", &event.track_cluster_de);
  treetpc->Branch( "track_cluster_size", &event.track_cluster_size);
  treetpc->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  treetpc->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  treetpc->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  treetpc->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  treetpc->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  treetpc->Branch( "track_cluster_row_center", &event.track_cluster_row_center);
  treetpc->Branch( "ntK18", &event.ntK18);
  treetpc->Branch( "chisqrK18", &event.chisqrK18);
  treetpc->Branch( "pK18", &event.pK18);
  treetpc->Branch( "xtgtK18", &event.xtgtK18);
  treetpc->Branch( "ytgtK18", &event.ytgtK18);
  treetpc->Branch( "utgtK18", &event.utgtK18);
  treetpc->Branch( "vtgtK18", &event.vtgtK18);
  treetpc->Branch( "thetaK18", &event.thetaK18);
  treetpc->Branch( "ntKurama",     &event.ntKurama);
  treetpc->Branch( "chisqrKurama", &event.chisqrKurama);
  treetpc->Branch( "pKurama",      &event.pKurama);
  treetpc->Branch( "qKurama",      &event.qKurama);
  treetpc->Branch( "m2",           &event.m2);
  //dsttpc->Branch( "m2Org",        &event.m2Org);
  treetpc->Branch( "xtgtKurama",   &event.xtgtKurama);
  treetpc->Branch( "ytgtKurama",   &event.ytgtKurama);
  treetpc->Branch( "utgtKurama",   &event.utgtKurama);
  treetpc->Branch( "vtgtKurama",   &event.vtgtKurama);  
  treetpc->Branch( "tpcidTPCKurama", &event.tpcidTPCKurama);  
  treetpc->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  treetpc->Branch( "kflagTPCKurama", &event.kflagTPCKurama);
  treetpc->Branch( "pflagTPCKurama", &event.pflagTPCKurama);
  treetpc->Branch( "chisqrTPCKurama", &event.chisqrTPCKurama);
  treetpc->Branch( "pTPCKurama", &event.pTPCKurama);
  treetpc->Branch( "qTPCKurama", &event.qTPCKurama);
  treetpc->Branch( "m2TPCKurama", &event.m2TPCKurama);
  treetpc->Branch( "xsTPC", &event.xsTPC);
  treetpc->Branch( "ysTPC", &event.ysTPC);
  treetpc->Branch( "usTPC", &event.usTPC);
  treetpc->Branch( "vsTPC", &event.vsTPC);
  treetpc->Branch( "xbTPC", &event.xbTPC);
  treetpc->Branch( "ybTPC", &event.ybTPC);
  treetpc->Branch( "ubTPC", &event.ubTPC);
  treetpc->Branch( "vbTPC", &event.vbTPC);
  treetpc->Branch( "xtgtTPCKurama",   &event.xtgtTPCKurama);
  treetpc->Branch( "ytgtTPCKurama",   &event.ytgtTPCKurama);
  treetpc->Branch( "utgtTPCKurama",   &event.utgtTPCKurama);
  treetpc->Branch( "vtgtTPCKurama",   &event.vtgtTPCKurama);  
  treetpc->Branch("nKm",           &event.nKm);
  treetpc->Branch("nKp",           &event.nKp);
  treetpc->Branch("nKK",           &event.nKK);
  treetpc->Branch("vtx",           &event.vtx);
  treetpc->Branch("vty",           &event.vty);
  treetpc->Branch("vtz",           &event.vtz);
  treetpc->Branch("closeDist",     &event.closeDist);
  treetpc->Branch("inside",        &event.inside);
  treetpc->Branch("MissMass",      &event.MissMass);
  treetpc->Branch("MissMassCorr",  &event.MissMassCorr);
  treetpc->Branch("MissMassCorrDE", &event.MissMassCorrDE);
  treetpc->Branch("pOrg",       &event.pOrg);
  treetpc->Branch("pCalc",      &event.pCalc);
  treetpc->Branch("pCorr",      &event.pCorr);
  treetpc->Branch("pCorrDE",    &event.pCorrDE);
  treetpc->Branch("BE", &event.BE);
  treetpc->Branch("BETPC", &event.BETPC);
  treetpc->Branch("BE_LL", &event.BE_LL);
  treetpc->Branch("BETPC_LL", &event.BETPC_LL);  
  treetpc->Branch("Kflag",      &event.Kflag);
  treetpc->Branch("Pflag",      &event.Pflag);
  treetpc->Branch("Heavyflag",  &event.Heavyflag);  
  treetpc->Branch( "nvtxTpc", &event.nvtxTpc );
  treetpc->Branch( "vtx_x", &event.vtx_x );
  treetpc->Branch( "vtx_y", &event.vtx_y );
  treetpc->Branch( "vtx_z", &event.vtx_z );
  treetpc->Branch( "vtx_dist", &event.vtx_dist );
  treetpc->Branch( "vtx_angle", &event.vtx_angle );
  treetpc->Branch( "vtxid", &event.vtxid );
  treetpc->Branch( "vtxmom_theta", &event.vtxmom_theta );
  treetpc->Branch( "vtxpos_x", &event.vtxpos_x );
  treetpc->Branch( "vtxpos_y", &event.vtxpos_y );
  treetpc->Branch( "vtxpos_z", &event.vtxpos_z );
  treetpc->Branch( "vtxmom_x", &event.vtxmom_x );
  treetpc->Branch( "vtxmom_y", &event.vtxmom_y );
  treetpc->Branch( "vtxmom_z", &event.vtxmom_z );
  treetpc->Branch( "isgoodTPC", &event.isgoodTPC);
  treetpc->Branch( "insideTPC", &event.insideTPC);
  treetpc->Branch( "vtxTPC", &event.vtxTPC);
  treetpc->Branch( "vtyTPC", &event.vtyTPC);
  treetpc->Branch( "vtzTPC", &event.vtzTPC);
  treetpc->Branch( "closeDistTPC", &event.closeDistTPC);
  treetpc->Branch( "MissMassTPC", &event.MissMassTPC);
  treetpc->Branch( "MissMassCorrTPC", &event.MissMassCorrTPC);
  treetpc->Branch( "MissMassCorrDETPC", &event.MissMassCorrDETPC);
  treetpc->Branch( "MissMassNuclTPC", &event.MissMassNuclTPC);
  treetpc->Branch( "MissMassNuclCorrTPC", &event.MissMassNuclCorrTPC);
  treetpc->Branch( "MissMassNuclCorrDETPC", &event.MissMassNuclCorrDETPC);
  treetpc->Branch( "pOrgTPC", &event.pOrgTPC);
  treetpc->Branch( "pCorrTPC", &event.pCorrTPC);
  treetpc->Branch( "pCorrDETPC", &event.pCorrDETPC);
  treetpc->Branch( "pCalcTPC", &event.pCalcTPC);
  treetpc->Branch( "thetaCMTPC", &event.thetaCMTPC);
  treetpc->Branch( "costCMTPC", &event.costCMTPC);
  treetpc->Branch( "thetaTPC", &event.thetaTPC);
  treetpc->Branch( "xbTPC", &event.xbTPC);
  treetpc->Branch( "ybTPC", &event.ybTPC);
  treetpc->Branch( "ubTPC", &event.ubTPC);
  treetpc->Branch( "vbTPC", &event.vbTPC);
  treetpc->Branch( "xsTPC", &event.xsTPC);
  treetpc->Branch( "ysTPC", &event.ysTPC);
  treetpc->Branch( "usTPC", &event.usTPC);
  treetpc->Branch( "vsTPC", &event.vsTPC);
  treetpc->Branch("Lflag", &event.lflag);
  treetpc->Branch("LToKuramaPflag", &event.kuramalflag);  
  treetpc->Branch("K0flag", &event.k0flag);
  treetpc->Branch("K0Mass", &event.k0mass);
  treetpc->Branch("K0DecayVtx_x", &event.k0decayvtx_x);
  treetpc->Branch("K0DecayVtx_y", &event.k0decayvtx_y);
  treetpc->Branch("K0DecayVtx_z", &event.k0decayvtx_z);
  treetpc->Branch("K0Mom_x", &event.k0mom_x);
  treetpc->Branch("K0Mom_y", &event.k0mom_y);
  treetpc->Branch("K0Mom_z", &event.k0mom_z);
  treetpc->Branch("K0VtxCloseDist", &event.GFk0pipi_dist);
  treetpc->Branch("K0DecaysTrackId", &event.k0decays_id);
  treetpc->Branch("K0DecaysMom", &event.k0decays_mom);
  treetpc->Branch("K0DecaysMom_x", &event.k0decays_mom_x);
  treetpc->Branch("K0DecaysMom_y", &event.k0decays_mom_y);
  treetpc->Branch("K0DecaysMom_z", &event.k0decays_mom_z);  
  treetpc->Branch("K0DecaysHtofPos_x", &event.k0decays_htofpos_x);
  treetpc->Branch("K0DecaysHtofPos_y", &event.k0decays_htofpos_y);
  treetpc->Branch("K0DecaysHtofPos_z", &event.k0decays_htofpos_z);
  treetpc->Branch("K0DecaysMass2", &event.k0decays_mass2);
  treetpc->Branch("K0DecaysTrackLen", &event.k0decays_tracklen);
  treetpc->Branch("GFstatus", &event.GFstatus);
  treetpc->Branch("GFntTpc", &event.GFntTpc);
  treetpc->Branch("GFcharge", &event.GFcharge);
  treetpc->Branch("GFchisqr", &event.GFchisqr);
  treetpc->Branch("GFtof", &event.GFtof);
  treetpc->Branch("GFpval", &event.GFpval);
  treetpc->Branch("GFfitstatus", &event.GFfitstatus);
  treetpc->Branch("GFpdgcode", &event.GFpdgcode);
  treetpc->Branch("GFnhtrack", &event.GFnhtrack);
  treetpc->Branch("GFlayer", &event.GFlayer);
  treetpc->Branch("GFpos_x", &event.GFpos_x);
  treetpc->Branch("GFpos_y", &event.GFpos_y);
  treetpc->Branch("GFpos_z", &event.GFpos_z);
  treetpc->Branch("GFmom", &event.GFmom);
  treetpc->Branch("GFmom_x", &event.GFmom_x);
  treetpc->Branch("GFmom_y", &event.GFmom_y);
  treetpc->Branch("GFmom_z", &event.GFmom_z);
  treetpc->Branch("GFresidual_x", &event.GFresidual_x);
  treetpc->Branch("GFresidual_y", &event.GFresidual_y);
  treetpc->Branch("GFresidual_z", &event.GFresidual_z);
  treetpc->Branch("GFresidual_p", &event.GFresidual_p);
  treetpc->Branch("GFresidual_px", &event.GFresidual_px);
  treetpc->Branch("GFresidual_py", &event.GFresidual_py);
  treetpc->Branch("GFresidual_pz", &event.GFresidual_pz);
  treetpc->Branch("KmFlag", &event.kmflag);
  treetpc->Branch("KmIncFlag", &event.kminc);    
  treetpc->Branch("GFKmTrackId", &event.GFkmid);  
  treetpc->Branch("GFKmMom", &event.GFkmmom);
  treetpc->Branch("GFKmMom_x", &event.GFkmmom_x);
  treetpc->Branch("GFKmMom_y", &event.GFkmmom_y);
  treetpc->Branch("GFKmMom_z", &event.GFkmmom_z);
  treetpc->Branch("GFKmTheta", &event.GFkmtheta);
  treetpc->Branch("GFKmPhi", &event.GFkmphi);    
  treetpc->Branch("GFKmTargetCloseDist", &event.GFkmtarget_dist);
  treetpc->Branch("GFKmTarget_x", &event.GFkmtargetvtx_x);
  treetpc->Branch("GFKmTarget_y", &event.GFkmtargetvtx_y);
  treetpc->Branch("GFKmTarget_z", &event.GFkmtargetvtx_z);
  treetpc->Branch("GFKmTargetCenter_x", &event.GFkmtargetcenter_x);
  treetpc->Branch("GFKmTargetCenter_y", &event.GFkmtargetcenter_y);
  treetpc->Branch("GFKmTargetCenter_z", &event.GFkmtargetcenter_z);
  treetpc->Branch("GFKmTargetCenterCloseDist", &event.GFkmtargetcenter_dist);
  treetpc->Branch("GFKmHtofId", &event.GFkmhtofid);
  treetpc->Branch("GFKmHtofSeg", &event.GFkmhtofseg);
  treetpc->Branch("GFKmHtofPos", &event.GFkmposHtof);  
  treetpc->Branch("GFKmMassSquare", &event.GFkmmass2);
  treetpc->Branch("GFKmInvBeta", &event.GFkminvbeta);         
  treetpc->Branch("GFKmTrackLen", &event.GFkmtracklen);
  treetpc->Branch("GFKmTof", &event.GFkmtof);
  treetpc->Branch("GFK0Mass", &event.GFk0mass);
  treetpc->Branch("GFK0DecayVtx_x", &event.GFk0decayvtx_x);
  treetpc->Branch("GFK0DecayVtx_y", &event.GFk0decayvtx_y);
  treetpc->Branch("GFK0DecayVtx_z", &event.GFk0decayvtx_z);
  treetpc->Branch("GFK0Mom", &event.GFk0mom);
  treetpc->Branch("GFK0Mom_x", &event.GFk0mom_x);
  treetpc->Branch("GFK0Mom_y", &event.GFk0mom_y);
  treetpc->Branch("GFK0Mom_z", &event.GFk0mom_z);
  treetpc->Branch("GFK0VtxCloseDist", &event.GFk0pipi_dist);
  treetpc->Branch("GFK0TargetCloseDist", &event.GFk0target_dist);
  treetpc->Branch("GFK0Target_x", &event.GFk0targetvtx_x);
  treetpc->Branch("GFK0Target_y", &event.GFk0targetvtx_y);
  treetpc->Branch("GFK0Target_z", &event.GFk0targetvtx_z);
  treetpc->Branch("GFK0TargetCenter_x", &event.GFk0targetcenter_x);
  treetpc->Branch("GFK0TargetCenter_y", &event.GFk0targetcenter_y);
  treetpc->Branch("GFK0TargetCenter_z", &event.GFk0targetcenter_z);
  treetpc->Branch("GFK0TargetCenterCloseDist", &event.GFk0targetcenter_dist);
  treetpc->Branch("GFK0ProductionVtx_x", &event.GFk0prodvtx_x);
  treetpc->Branch("GFK0ProductionVtx_y", &event.GFk0prodvtx_y);
  treetpc->Branch("GFK0ProductionVtx_z", &event.GFk0prodvtx_z);
  treetpc->Branch("GFK0ProductionVtxCloseDist", &event.GFk0prodvtx_dist);
  treetpc->Branch("GFK0TrackLen", &event.GFk0tracklen);
  treetpc->Branch("GFK0Tof", &event.GFk0tof);  
  treetpc->Branch("GFntTpc_inside", &event.GFntTpc_inside);
  treetpc->Branch("GFprodvtx_x", &event.GFprodvtx_x);
  treetpc->Branch("GFprodvtx_y", &event.GFprodvtx_y);
  treetpc->Branch("GFprodvtx_z", &event.GFprodvtx_z);
  treetpc->Branch("GFinside", &event.GFinside);
  treetpc->Branch("GFinsideTgtToKurama", &event.GFinsideKurama);
  treetpc->Branch("GFKuramaHasVtxOutOfTgt", &event.GFKuramaVtxOutTgt);
  treetpc->Branch("GFForwardYDecayLikeTrack", &event.GFFwdYDecayTrack);  
  treetpc->Branch("GFfromVtx", &event.GFfromVtx);  
  treetpc->Branch("GFextrapolationHtof", &event.GFextrapolationHtof); 
  treetpc->Branch("GFtracklen", &event.GFtracklen);
  treetpc->Branch("GFtrack2vtxdist", &event.GFtrack2vtxdist);
  treetpc->Branch("GFcalctof", &event.GFcalctof);
  treetpc->Branch("GFsegHtof", &event.GFsegHtof);
  treetpc->Branch("GFtofHtof", &event.GFtofHtof);
  treetpc->Branch("GFtdiffHtof", &event.GFtdiffHtof);
  treetpc->Branch("GFposHtof", &event.GFposHtof);
  treetpc->Branch("GFposx", &event.GFposx);
  treetpc->Branch("GFposy", &event.GFposy);
  treetpc->Branch("GFposz", &event.GFposz);
  treetpc->Branch("GFinvbeta", &event.GFinvbeta);
  treetpc->Branch("GFm2", &event.GFm2);
  treetpc->Branch("nsigma_tritonHtof", &event.nsigma_tritonHtof);
  treetpc->Branch("nsigma_deutronHtof", &event.nsigma_deutronHtof);
  treetpc->Branch("nsigma_protonHtof", &event.nsigma_protonHtof);
  treetpc->Branch("nsigma_kaonHtof", &event.nsigma_kaonHtof);
  treetpc->Branch("nsigma_pionHtof", &event.nsigma_pionHtof);
  treetpc->Branch("nsigma_electronHtof", &event.nsigma_electronHtof);

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
  src.thetaK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaK18" );  
  src.chisqrK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrK18" );
  src.xtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtK18" );
  src.ytgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtK18" );
  src.utgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtK18" );
  src.vtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtK18" );

  src.ntKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.chisqrKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrKurama" );
  src.pKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pKurama" );
  src.qKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qKurama" );
  src.xtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );
  
  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" );
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

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.chisqr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr" );
  src.helix_cx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx" );
  src.helix_cy = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy" );
  src.helix_z0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0" );
  src.helix_r = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r" );
  src.helix_dz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz" );
  src.dE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dE" );
  src.dEdx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "dEdx" );
  src.mom0 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.path = new TTreeReaderValue<std::vector<Double_t>>( *reader, "path" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );
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
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "alpha" );
  src.pathhit = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "pathhit" );
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de" );
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_size" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );
  src.track_cluster_de_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de_center" );
  src.track_cluster_x_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_x_center" );
  src.track_cluster_y_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_y_center" );
  src.track_cluster_z_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_z_center" );
  src.track_cluster_row_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_row_center" );

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
  src.Kflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "Kflag" );
  src.Pflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "Pflag" );
  src.Heavyflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "Heavyflag" );

  src.ntTPCK18 = new TTreeReaderValue<Int_t>( *reader, "ntK18" );
  src.chisqrK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrK18" );
  src.xtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtK18" );
  src.ytgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtK18" );
  src.utgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtK18" );
  src.vtgtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtK18" );
  src.isgoodTPCK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCK18" );
  src.chisqrTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrTPCK18" );
  src.qTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCK18");
  src.pTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCK18");
  src.thetaTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPCK18");  
  src.xtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtTPCK18" );
  src.ytgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtTPCK18" );
  src.utgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtTPCK18" );
  src.vtgtTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtTPCK18" );
  src.lhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "lhtofTPCK18" );
  src.xhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofTPCK18" );
  src.yhtofTPCK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofTPCK18" );
  src.lvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "lvpTPCK18" );
  src.xvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpTPCK18" );
  src.yvpTPCK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpTPCK18" );
  
  src.ntTPCKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.chisqrKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqrKurama" );
  src.pKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pKurama" );
  src.qKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qKurama" );
  //  src.m2  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2" );
  src.xtgtKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );
  src.tpcidTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "tpcidTPCKurama" );
  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
  src.kflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "kflagTPCKurama" );
  src.pflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pflagTPCKurama" );
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
  src.xsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xsTPC" );
  src.ysTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ysTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );

  src.isgoodTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPC" );
  src.insideTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "insideTPC" );
  src.vtxTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtxTPC" );
  src.vtyTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtyTPC" );
  src.vtzTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtzTPC" );
  src.closeDistTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "closeDistTPC" );
  src.MissMassTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassTPC" );
  src.MissMassCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrTPC" );
  src.MissMassCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrDETPC" );
  src.MissMassNuclTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassNuclTPC" );
  src.MissMassNuclCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassNuclCorrTPC" );
  src.MissMassNuclCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassNuclCorrDETPC" );

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

  //src.lflag		  = new TTreeReaderValue<Int_t>(*reader,"Lflag");			
  // src.lmass		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMass");			
  // src.ldecayvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_x");		
  // src.ldecayvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_y");		
  // src.ldecayvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_z");		
  // src.lmom		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom");			
  // src.lmom_x		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_x");			
  // src.lmom_y		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_y");			
  // src.lmom_z		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_z");			
  // src.ppi_dist          = new TTreeReaderValue<Double_t>(*reader,"LambdaVtxCloseDist");		
  // // src.ltarget_dist	  = new TTreeReaderValue<Double_t>(*reader,"LambdaTargetCloseDist");	
  // // src.ltargetvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"LambdaTarget_x");		
  // // src.ltargetvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"LambdaTarget_y");		
  // // src.ltargetvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"LambdaTarget_z");		
  // // src.ltargetcenter_x	  = new TTreeReaderValue<Double_t>(*reader,"LambdaTargetCenter_x");	
  // // src.ltargetcenter_y	  = new TTreeReaderValue<Double_t>(*reader,"LambdaTargetCenter_y");	
  // // src.ltargetcenter_z	  = new TTreeReaderValue<Double_t>(*reader,"LambdaTargetCenter_z");	
  // // src.ltargetcenter_dist= new TTreeReaderValue<Double_t>(*reader,"LambdaTargetCenterCloseDist");	
  // // src.lprodvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"LambdaProductionVtx_x");	
  // // src.lprodvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"LambdaProductionVtx_y");	
  // // src.lprodvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"LambdaProductionVtx_z");	
  // // src.lprodvtx_dist	  = new TTreeReaderValue<Double_t>(*reader,"LambdaProductionVtxCloseDist");
  // // src.ltracklen         = new TTreeReaderValue<Double_t>(*reader,"LambdaTrackLen");
  // // src.ltof              = new TTreeReaderValue<Double_t>(*reader,"LambdaTof");                     
  // src.ldecays_id        = new TTreeReaderValue<std::vector<Int_t>>(*reader,"LDecaysTrackId");
  // src.ldecays_mom       = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysMom");
  // src.ldecays_mom_x     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysMom_x");
  // src.ldecays_mom_y     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysMom_y");
  // src.ldecays_mom_z     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysMom_z");    
  // src.ldecays_htofhitid = new TTreeReaderValue<std::vector<Int_t>>(*reader,"LDecaysHtofHitId");
  // src.ldecays_htofseg   = new TTreeReaderValue<std::vector<Int_t>>(*reader,"LDecaysHtofSeg");    
  // src.ldecays_tracklen  = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysTrackLen");
  // src.ldecays_invbeta   = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysInvBeta");
  // src.ldecays_mass2     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysMass2");
  // src.ldecays_htofpos_x = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysHtofPos_x");
  // src.ldecays_htofpos_y = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysHtofPos_y");
  // src.ldecays_htofpos_z = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysHtofPos_z");  
  
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO") &&
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
