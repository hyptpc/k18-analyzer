// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
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
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCVertex.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"
#include "PidLikelihoodMan.hh"

#define SaveHistograms 1
#define RawCluster 1
#define UsePidLH 1

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPidLike = PidLikelihoodMan::GetInstance();  
bool usePidLH = false;

const int nbinpoq = pidlikeli::nbinpoq;
const double minpoq = pidlikeli::minpoq;
const double maxpoq = pidlikeli::maxpoq; //GeV/c
const int nbindedx = pidlikeli::nbindedx;
const double mindedx = pidlikeli::mindedx;
const double maxdedx = pidlikeli::maxdedx;
const int nbinm2 = pidlikeli::nbinm2;
const double minm2 = pidlikeli::minm2;
const double maxm2 = pidlikeli::maxm2;
  
const int kPmin = pidlikeli::kPmin;
const int kPmax = pidlikeli::kPmax;
const double kMomstep = pidlikeli::kDP;
const int kNmom = pidlikeli::kNmom;  
const int fac_t = pidlikeli::fac_t;
const int fac_p = pidlikeli::fac_p;
const int fac_c = pidlikeli::fac_c;
const int fac_b = pidlikeli::fac_b;
const int fac_m = pidlikeli::fac_m;
const auto& type = pidlikeli::type;  
const int typeGenHid = pidlikeli::kTypeGen;  
const int typeLHid = pidlikeli::kTypeLmd;
const int typeK0Hid = pidlikeli::kTypeK0;
const int typeKmHid = pidlikeli::kTypeKm;
const auto& plist = pidlikeli::plist;
const int kNpid = static_cast<int>(pidlikeli::Pid::COUNT);  
const int kPidPi = static_cast<int>(pidlikeli::Pid::Pi);
const int kPidK  = static_cast<int>(pidlikeli::Pid::K);
const int kPidP  = static_cast<int>(pidlikeli::Pid::P);
const int kPidD  = static_cast<int>(pidlikeli::Pid::D);
const int kPidE  = static_cast<int>(pidlikeli::Pid::E);
const int kPidAll = pidlikeli::kAllParticles;  
const auto& clist = pidlikeli::clist;
const int kNbe = static_cast<int>(pidlikeli::BE::COUNT);//pidlikeli::kNbe;
const int kNtype = pidlikeli::kNtype;
const int kNchg = pidlikeli::kNchg;
const int plusHid = pidlikeli::kPlus;
const int minusHid = pidlikeli::kMinus;
  
  
const double minbe = -0.5; //GeV
const double maxbe = 0.5;
const double bestep = 0.20; //GeV // should be chaged to USER parameter
  
const Double_t cutm2proton = 0.3;

const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;
  
const Double_t vtx_cut_twotrack = 50.;

const Double_t lambda_masscut = 0.1;
const Double_t lambda_masscut_final = 0.02; //final  
const Double_t k0_masscut = 0.1; //final
// const Double_t p_vtx_distcut = 300;
// const Double_t pi_vtx_distcut = 300;
// const Double_t p_vtx_distcut = 300;
// const Double_t pi_vtx_distcut = 300;
const Double_t p_vtx_distcut = 100;
const Double_t pi_vtx_distcut = 100;
const Double_t e_vtx_distcut = 300;
  //const Double_t ppi_distcut = 100; //ref
const Double_t ppi_distcut = 10; //ref
  
const Double_t ltarget_distcut = 25.;
// const Double_t pip_vtx_distcut = 300;
// const Double_t pim_vtx_distcut = 300;
const Double_t pip_vtx_distcut = 50;
const Double_t pim_vtx_distcut = 50;  
const Double_t pipi_distcut = 10.; //ref  
const Double_t k0target_distcut = 25.;

const Double_t residual_track_distcut = 25.;
const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");  
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpc, kKScat, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[TPCKuramaK18Tracking]", "[KScat]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "kk","" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
Double_t nsigma_m2 = 3.; //M2 cut for TPCKurama RK
static TString eq_sigmaM2 = "4*TMath::Power([3], 2)*(1.+[3]/(x*x))*[0]+4*TMath::Power([3], 2)*x*x*[1]+4*x*x*(x*x+[3])*[2]";
static TString eq_M2 = "x*x*[2]+x*[1]+[0]";
}

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  Int_t runnum;
  Int_t evnum;
  std::vector<Int_t> trigpat;
  std::vector<Int_t> trigflag;

  Int_t nhTpc; // Number of clusters
  Int_t nclTpc; // Number of clusters
  Int_t remain_nclTpc; // Number of clusters without tracks
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
  Int_t ntKuramaCandidate; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)
  std::vector<Int_t> isKuramaCandidate;
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> trackid; //for Kurama K1.8 tracks
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isXi;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> isMultiloop;
  std::vector<Int_t> charge; //Helix charge
  std::vector<Int_t> pid;
  std::vector<Double_t> chisqr;
  std::vector<Double_t> pval;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx
  std::vector<Double_t> mom0; //Helix momentum at Y = 0
  std::vector<Double_t> path; //Helix path
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

  std::vector<Int_t> chargeIndistinguishable;
  std::vector<Double_t> chisqr_inverted;
  std::vector<Double_t> pval_inverted;
  std::vector<Double_t> helix_cx_inverted;
  std::vector<Double_t> helix_cy_inverted;
  std::vector<Double_t> helix_z0_inverted;
  std::vector<Double_t> helix_r_inverted;
  std::vector<Double_t> helix_dz_inverted;
  std::vector<Double_t> mom0_inverted;//Helix momentum at Y = 0
  std::vector<Int_t> pid_inverted;

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
  std::vector<Double_t> clusteredVtx_x;
  std::vector<Double_t> clusteredVtx_y;
  std::vector<Double_t> clusteredVtx_z;
  std::vector<std::vector<Double_t>> clusteredVtxid;

  Int_t ntK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> p_3rd;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> xoutK18;
  std::vector<Double_t> youtK18;
  std::vector<Double_t> uoutK18;
  std::vector<Double_t> voutK18;
  std::vector<Double_t> xtgtK18;
  std::vector<Double_t> ytgtK18;
  std::vector<Double_t> utgtK18;
  std::vector<Double_t> vtgtK18;
  std::vector<Double_t> thetaK18;
  std::vector<Double_t> xhtofK18;
  std::vector<Double_t> yhtofK18;
  std::vector<std::vector<Double_t>> xvpHS;
  std::vector<std::vector<Double_t>> yvpHS;
  std::vector<std::vector<Double_t>> zvpHS;
  std::vector<Double_t> xtgtHS;
  std::vector<Double_t> ytgtHS;
  std::vector<Double_t> ztgtHS;
  std::vector<std::vector<Double_t>> layerK18;
  std::vector<std::vector<Double_t>> wireK18;
  std::vector<std::vector<Double_t>> localhitposK18;
  std::vector<std::vector<Double_t>> wposK18;

  Int_t ntKurama;
  std::vector<Double_t> chisqrKurama;
  std::vector<Double_t> pKurama;
  std::vector<Double_t> qKurama;
  std::vector<Double_t> m2;
  std::vector<Double_t> m2Org;
  std::vector<Double_t> xtgtKurama;
  std::vector<Double_t> ytgtKurama;
  std::vector<Double_t> utgtKurama;
  std::vector<Double_t> vtgtKurama;
  std::vector<Double_t> thetaKurama;
  std::vector<Double_t> pathKurama;
  std::vector<Double_t> xhtofKurama;
  std::vector<Double_t> yhtofKurama;
  std::vector<Double_t> cstof;
  std::vector<Double_t> tofsegKurama;
  std::vector<Double_t> pathwcKurama;
  std::vector<std::vector<Double_t>> xvpKurama;
  std::vector<std::vector<Double_t>> yvpKurama;
  std::vector<std::vector<Double_t>> zvpKurama;
  std::vector<Double_t> xin;
  std::vector<Double_t> yin;
  std::vector<Double_t> zin;
  std::vector<Double_t> pxin;
  std::vector<Double_t> pyin;
  std::vector<Double_t> pzin;
  std::vector<Double_t> xout;
  std::vector<Double_t> yout;
  std::vector<Double_t> zout;
  std::vector<Double_t> pxout;
  std::vector<Double_t> pyout;
  std::vector<Double_t> pzout;
  std::vector<std::vector<Double_t>> layer;
  std::vector<std::vector<Double_t>> wire;
  std::vector<std::vector<Double_t>> localhitpos;
  std::vector<std::vector<Double_t>> wpos;

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
  std::vector<Double_t> pOrg;
  std::vector<Double_t> pCalc;
  std::vector<Double_t> pCorr;
  std::vector<Double_t> pCorrDE;

  std::vector<Double_t> xkm;
  std::vector<Double_t> ykm;
  std::vector<Double_t> ukm;
  std::vector<Double_t> vkm;
  std::vector<Double_t> xkp;
  std::vector<Double_t> ykp;
  std::vector<Double_t> ukp;
  std::vector<Double_t> vkp;
  std::vector<Int_t> Kflag;
  std::vector<Int_t> Pflag;
  std::vector<Int_t> Heavyflag;

  //TPC RK
  std::vector<Int_t> isgoodTPCK18;
  std::vector<Int_t> tpcidTPCK18;
  std::vector<Double_t> chisqrTPCK18;
  std::vector<Double_t> pTPCK18;
  std::vector<Double_t> qTPCK18;
  std::vector<Double_t> xtgtTPCK18;
  std::vector<Double_t> ytgtTPCK18;
  std::vector<Double_t> utgtTPCK18;
  std::vector<Double_t> vtgtTPCK18;
  std::vector<Double_t> thetaTPCK18;
  std::vector<Double_t> lhtofTPCK18;
  std::vector<Double_t> xhtofTPCK18;
  std::vector<Double_t> yhtofTPCK18;
  std::vector<std::vector<Double_t>> lvpTPCK18;
  std::vector<std::vector<Double_t>> xvpTPCK18;
  std::vector<std::vector<Double_t>> yvpTPCK18;

  std::vector<Int_t> tpcidTPCKurama;
  std::vector<Int_t> isgoodTPCKurama;
  //std::vector<Int_t> piflagTPCKurama;
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
  std::vector<Double_t> lgasvesselTPCKurama;
  std::vector<Double_t> xgasvesselTPCKurama;
  std::vector<Double_t> ygasvesselTPCKurama;
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
  std::vector<Double_t> pOrgTPC;
  std::vector<Double_t> pCorrTPC;
  std::vector<Double_t> pCorrDETPC;
  std::vector<Double_t> pCalcTPC;
  std::vector<Double_t> thetaCMTPC;
  std::vector<Double_t> costCMTPC;
  std::vector<Double_t> pCalcDETPC;
  std::vector<Double_t> thetaCMDETPC;
  std::vector<Double_t> costCMDETPC;
  std::vector<Double_t> xistarpCalcDETPC;
  std::vector<Double_t> xistarthetaCMDETPC;
  std::vector<Double_t> xistarcostCMDETPC;
  std::vector<Double_t> kpscatpCalcTPC;
  std::vector<Double_t> kpscatthetaCMTPC;
  std::vector<Double_t> kpscatcostCMTPC;
  std::vector<Double_t> kpscatpCalcDETPC;
  std::vector<Double_t> kpscatthetaCMDETPC;
  std::vector<Double_t> kpscatcostCMDETPC;
  std::vector<Double_t> thetaTPC;
  std::vector<Double_t> xbTPC;
  std::vector<Double_t> ybTPC;
  std::vector<Double_t> ubTPC;
  std::vector<Double_t> vbTPC;
  std::vector<Double_t> xsTPC;
  std::vector<Double_t> ysTPC;
  std::vector<Double_t> usTPC;
  std::vector<Double_t> vsTPC;

  std::vector<Bool_t>   insideTgt;
  std::vector<Bool_t>   extrapTgt;
  std::vector<Bool_t>   extrapVtx;
  Int_t ntTpc_inside;
  Double_t prodvtx_x;
  Double_t prodvtx_y;
  Double_t prodvtx_z;  
  std::vector<Double_t> m2HtofVtx;
  std::vector<Int_t> hitidHtof;
  std::vector<Int_t> segHtof;
  std::vector<Double_t> tofHtof;
  std::vector<Double_t> invbetaHtof;    
  std::vector<Double_t> tracklenHtof;
  std::vector<Double_t> posHtof_x;
  std::vector<Double_t> posHtof_y;
  std::vector<Double_t> posHtof_z;
  std::vector<Double_t> distVtx;

  Int_t lflag;
  Double_t lmass;
  Double_t ldecayvtx_x;
  Double_t ldecayvtx_y;
  Double_t ldecayvtx_z;
  Double_t lmom;  
  Double_t lmom_x;
  Double_t lmom_y;
  Double_t lmom_z;
  Double_t ppi_dist;
  Double_t ppiangle;    
  std::vector<Int_t> ldecays_id;
  std::vector<Double_t> ldecays_mom;
  std::vector<Double_t> ldecays_mom_x;
  std::vector<Double_t> ldecays_mom_y;
  std::vector<Double_t> ldecays_mom_z;
  std::vector<Bool_t>   ldecays_htofextrap;  
  std::vector<Int_t>    ldecays_htofhitid;
  std::vector<Int_t>    ldecays_htofseg;    
  std::vector<Double_t> ldecays_mass2;
  std::vector<Double_t> ldecays_invbeta;  
  std::vector<Double_t> ldecays_tracklen;
  std::vector<Double_t> ldecays_htofpos_x;
  std::vector<Double_t> ldecays_htofpos_y;
  std::vector<Double_t> ldecays_htofpos_z;   

  Bool_t k0flag;
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
  std::vector<Bool_t>   k0decays_htofextrap;  
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
  std::vector<Double_t> m2Htof;


  std::vector<Double_t> nsigma_tritonHtof;
  std::vector<Double_t> nsigma_deutronHtof;
  std::vector<Double_t> nsigma_protonHtof;
  std::vector<Double_t> nsigma_kaonHtof;
  std::vector<Double_t> nsigma_pionHtof;
  std::vector<Double_t> nsigma_electronHtof;    

  void clear( void )
  {
    trigpat.clear();
    trigflag.clear();

    runnum = 0;
    evnum = 0;
    status = 0;

    nhTpc = 0;
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
    ntKuramaCandidate = 0; //Numer of tracks which are kurama track candidates(before TPCKurama tracking) 
    isKuramaCandidate.clear();
    nhtrack.clear();
    trackid.clear();
    isBeam.clear();
    isXi.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    charge.clear();
    pid.clear();
    
    chisqr.clear();
    pval.clear();
    helix_cx.clear(); 
    helix_cy.clear(); 
    helix_z0.clear(); 
    helix_r.clear();
    helix_dz.clear(); 
    dE.clear(); 
    dEdx.clear();
    mom0.clear();
    path.clear();
    isElectron.clear();
    nsigma_triton.clear();
    nsigma_deutron.clear();
    nsigma_proton.clear();
    nsigma_kaon.clear();
    nsigma_pion.clear();
    nsigma_electron.clear();
    
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

    chargeIndistinguishable.clear();
    chisqr_inverted.clear();
    pval_inverted.clear();
    helix_cx_inverted.clear();
    helix_cy_inverted.clear();
    helix_z0_inverted.clear();
    helix_r_inverted.clear();
    helix_dz_inverted.clear();
    mom0_inverted.clear();
    pid_inverted.clear();

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

    nvtxTpcClustered = 0;
    clusteredVtx_x.clear();
    clusteredVtx_y.clear();
    clusteredVtx_z.clear();
    clusteredVtxid.clear();

    ntK18 = 0;
    pK18.clear();
    chisqrK18.clear();
    p_3rd.clear();
    xtgtK18.clear();
    ytgtK18.clear();
    utgtK18.clear();
    vtgtK18.clear();
    thetaK18.clear();
    xhtofK18.clear();
    yhtofK18.clear();
    xoutK18.clear();
    youtK18.clear();
    uoutK18.clear();
    voutK18.clear();
    xtgtHS.clear();
    ytgtHS.clear();
    ztgtHS.clear();
    xvpHS.clear();
    yvpHS.clear();
    zvpHS.clear();
    layerK18.clear();
    wireK18.clear();
    localhitposK18.clear();

    ntKurama = 0;
    chisqrKurama.clear();
    pKurama.clear();
    qKurama.clear();
    m2.clear();
    m2Org.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();
    thetaKurama.clear();
    pathKurama.clear();
    xhtofKurama.clear();
    yhtofKurama.clear();
    cstof.clear();
    tofsegKurama.clear();
    pathwcKurama.clear();
    xvpKurama.clear();
    yvpKurama.clear();
    zvpKurama.clear();
    layer.clear();
    wire.clear();
    localhitpos.clear();
    wpos.clear();
    xin.clear();
    yin.clear();
    zin.clear();
    pxin.clear();
    pyin.clear();
    pzin.clear();
    xout.clear();
    yout.clear();
    zout.clear();
    pxout.clear();
    pyout.clear();
    pzout.clear();

    nKm = 0;
    nKp = 0;
    nKK = 0;
    inside.clear();
    vtx.clear();
    vty.clear();
    vtz.clear();
    closeDist.clear();
    MissMass.clear();
    MissMassCorr.clear();
    MissMassCorrDE.clear();
    pOrg.clear();
    pCalc.clear();
    pCorr.clear();
    pCorrDE.clear();

    xkm.clear();
    ykm.clear();
    ukm.clear();
    vkm.clear();
    xkp.clear();
    ykp.clear();
    ukp.clear();
    vkp.clear();
    Kflag.clear();
    Pflag.clear();
    Heavyflag.clear();

    tpcidTPCK18.clear();
    isgoodTPCK18.clear();
    chisqrTPCK18.clear();
    qTPCK18.clear();
    pTPCK18.clear();
    xtgtTPCK18.clear();
    ytgtTPCK18.clear();
    utgtTPCK18.clear();
    vtgtTPCK18.clear();
    thetaTPCK18.clear();
    lhtofTPCK18.clear();
    xhtofTPCK18.clear();
    yhtofTPCK18.clear();
    lvpTPCK18.clear();
    xvpTPCK18.clear();
    yvpTPCK18.clear();

    tpcidTPCKurama.clear();
    isgoodTPCKurama.clear();
    //    piflagTPCKurama.clear();
    kflagTPCKurama.clear();
    pflagTPCKurama.clear();
    chisqrTPCKurama.clear();
    pTPCKurama.clear();
    qTPCKurama.clear();
    m2TPCKurama.clear();
    xtgtTPCKurama.clear();
    ytgtTPCKurama.clear();
    utgtTPCKurama.clear();
    vtgtTPCKurama.clear();
    thetaTPCKurama.clear();
    pathTPCKurama.clear();
    lhtofTPCKurama.clear();
    xhtofTPCKurama.clear();
    yhtofTPCKurama.clear();
    lgasvesselTPCKurama.clear();
    xgasvesselTPCKurama.clear();
    ygasvesselTPCKurama.clear();
    lvpTPCKurama.clear();
    xvpTPCKurama.clear();
    yvpTPCKurama.clear();

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
    pOrgTPC.clear();
    pCorrTPC.clear();
    pCorrDETPC.clear();
    pCalcTPC.clear();
    thetaCMTPC.clear();
    costCMTPC.clear();
    pCalcDETPC.clear();
    thetaCMDETPC.clear();
    costCMDETPC.clear();
    xistarpCalcDETPC.clear();
    xistarthetaCMDETPC.clear();
    xistarcostCMDETPC.clear();
    kpscatpCalcTPC.clear();
    kpscatthetaCMTPC.clear();
    kpscatcostCMTPC.clear();
    kpscatpCalcDETPC.clear();
    kpscatthetaCMDETPC.clear();
    kpscatcostCMDETPC.clear();
    thetaTPC.clear();
    xbTPC.clear();
    ybTPC.clear();
    ubTPC.clear();
    vbTPC.clear();
    xsTPC.clear();
    ysTPC.clear();
    usTPC.clear();
    vsTPC.clear();

    ntTpc_inside = qnan;
    insideTgt.clear();
    extrapTgt.clear();
    extrapVtx.clear();
    prodvtx_x = qnan;
    prodvtx_y = qnan;    
    prodvtx_z = qnan;      
    m2HtofVtx.clear();
    hitidHtof.clear();
    segHtof.clear();
    tofHtof.clear();
    invbetaHtof.clear(); 
    posHtof_x.clear();
    posHtof_y.clear();            
    posHtof_z.clear();              
    tracklenHtof.clear();
    distVtx.clear();

    lflag = false;
    lmass  = qnan;
    ldecayvtx_x = qnan;
    ldecayvtx_y = qnan;
    ldecayvtx_z = qnan;
    lmom   = qnan;    
    lmom_x = qnan;
    lmom_y = qnan;
    lmom_z = qnan;
    ppi_dist = qnan;
    ppiangle = qnan;    
    ldecays_id.clear();
    ldecays_mom.clear();
    ldecays_mom_x.clear();
    ldecays_mom_y.clear();
    ldecays_mom_z.clear();
    ldecays_htofextrap.clear();            
    ldecays_htofhitid.clear();
    ldecays_htofseg.clear();        
    ldecays_mass2.clear();    
    ldecays_tracklen.clear();
    ldecays_htofpos_x.clear();
    ldecays_htofpos_y.clear();
    ldecays_htofpos_z.clear();                
    
    k0flag = false;
    k0mass  = qnan;
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
    k0decays_htofextrap.clear();    
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
    m2Htof.clear();    

    nsigma_tritonHtof.clear();
    nsigma_deutronHtof.clear();
    nsigma_protonHtof.clear();
    nsigma_kaonHtof.clear();
    nsigma_pionHtof.clear();        
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

  TTreeReaderValue<Int_t>* nhTpc; // Number of clusters
  TTreeReaderValue<Int_t>* nclTpc; // Number of clusters
  TTreeReaderValue<Int_t>* remain_nclTpc; // Number of clusters without tracks
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

  TTreeReaderValue<Int_t>* ntTpc; // Number of tracks
  TTreeReaderValue<Int_t>* ntKuramaCandidate; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)
  TTreeReaderValue<std::vector<Int_t>>* isKuramaCandidate;
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isXi;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Int_t>>* isMultiloop;
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
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_size;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_mrow;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_de_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_x_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_y_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_z_center;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* track_cluster_row_center;

  TTreeReaderValue<std::vector<Int_t>>* chargeIndistinguishable;
  TTreeReaderValue<std::vector<Int_t>>* pid_inverted;
  TTreeReaderValue<std::vector<Double_t>>* chisqr_inverted;
  TTreeReaderValue<std::vector<Double_t>>* pval_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_r_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz_inverted;
  TTreeReaderValue<std::vector<Double_t>>* mom0_inverted;//Helix momentum at Y = 0

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

  TTreeReaderValue<Int_t>* nvtxTpcClustered;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_x;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_y;
  TTreeReaderValue<std::vector<Double_t>>* clusteredVtx_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* clusteredVtxid;

  TTreeReaderValue<Int_t>* ntTPCK18; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* tpcidTPCK18;
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* qTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* pTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* xtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPCK18;
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
  //TTreeReaderValue<std::vector<Int_t>>* piflagTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* kflagTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* pflagTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* pTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* qTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2TPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xtgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* ytgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* utgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* vtgtTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* pathTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* lhtofTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xhtofTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* yhtofTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* lgasvesselTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xgasvesselTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* ygasvesselTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* lvpTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpTPCKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xhtofKurama;
  TTreeReaderValue<std::vector<Double_t>>* yhtofKurama;

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
  TTreeReaderValue<std::vector<Double_t>>* pOrgTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* pCalcTPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* costCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCalcDETPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaCMDETPC;
  TTreeReaderValue<std::vector<Double_t>>* costCMDETPC;
  TTreeReaderValue<std::vector<Double_t>>* xistarpCalcDETPC;
  TTreeReaderValue<std::vector<Double_t>>* xistarthetaCMDETPC;
  TTreeReaderValue<std::vector<Double_t>>* xistarcostCMDETPC;
  TTreeReaderValue<std::vector<Double_t>>* kpscatpCalcTPC;
  TTreeReaderValue<std::vector<Double_t>>* kpscatthetaCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* kpscatcostCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* kpscatpCalcDETPC;
  TTreeReaderValue<std::vector<Double_t>>* kpscatthetaCMDETPC;
  TTreeReaderValue<std::vector<Double_t>>* kpscatcostCMDETPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;

  Int_t    ntK18;
  Double_t chisqrK18[MaxHits];
  Double_t pK18[MaxHits];
  Double_t xtgtK18[MaxHits];
  Double_t ytgtK18[MaxHits];
  Double_t utgtK18[MaxHits];
  Double_t vtgtK18[MaxHits];
  Double_t thetaK18[MaxHits];
  //From TPCKuramaK18Tracking
  TTreeReaderValue<std::vector<Double_t>>* p_3rd;
  TTreeReaderValue<std::vector<Double_t>>* xoutK18;
  TTreeReaderValue<std::vector<Double_t>>* youtK18;
  TTreeReaderValue<std::vector<Double_t>>* uoutK18;
  TTreeReaderValue<std::vector<Double_t>>* voutK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* zvpHS;
  TTreeReaderValue<std::vector<Double_t>>* xtgtHS;
  TTreeReaderValue<std::vector<Double_t>>* ytgtHS;
  TTreeReaderValue<std::vector<Double_t>>* ztgtHS;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* layerK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* wireK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* localhitposK18;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* wposK18;

  Int_t    ntKurama;
  Double_t chisqrKurama[MaxHits];
  Double_t pKurama[MaxHits];
  Double_t qKurama[MaxHits];
  Double_t m2[MaxHits];
  Double_t m2Org[MaxHits];
  Double_t xtgtKurama[MaxHits];
  Double_t ytgtKurama[MaxHits];
  Double_t utgtKurama[MaxHits];
  Double_t vtgtKurama[MaxHits];
  Double_t thetaKurama[MaxHits];
  Double_t pathKurama[MaxHits];
  Double_t cstof[MaxHits];
  Double_t tofsegKurama[MaxHits];
  Double_t pathwcKurama[MaxHits];
  Double_t xin[MaxHits];
  Double_t yin[MaxHits];
  Double_t zin[MaxHits];
  Double_t pxin[MaxHits];
  Double_t pyin[MaxHits];
  Double_t pzin[MaxHits];
  Double_t xout[MaxHits];
  Double_t yout[MaxHits];
  Double_t zout[MaxHits];
  Double_t pxout[MaxHits];
  Double_t pyout[MaxHits];
  Double_t pzout[MaxHits];

  //From TPCKuramaTracking
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* xvpKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* yvpKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* zvpKurama;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* layer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* wire;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* localhitpos;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* wpos;

  Int_t    nhHtof;
  Double_t HtofSeg[MaxHits];
  Double_t tHtof[MaxHits];
  Double_t dtHtof[MaxHits];
  Double_t deHtof[MaxHits];
  Double_t posHtof[MaxHits];

  //Reaction
  Int_t    nKm;
  Int_t    nKp;
  Int_t    nKK;
  Int_t    inside[MaxHits];
  Double_t vtx[MaxHits];
  Double_t vty[MaxHits];
  Double_t vtz[MaxHits];
  Double_t closeDist[MaxHits];
  Double_t MissMass[MaxHits];
  Double_t MissMassCorr[MaxHits];
  Double_t MissMassCorrDE[MaxHits];
  Double_t pOrg[MaxHits];
  Double_t pCalc[MaxHits];
  Double_t pCorr[MaxHits];
  Double_t pCorrDE[MaxHits];
  Double_t xkm[MaxHits];
  Double_t ykm[MaxHits];
  Double_t ukm[MaxHits];
  Double_t vkm[MaxHits];
  Double_t xkp[MaxHits];
  Double_t ykp[MaxHits];
  Double_t ukp[MaxHits];
  Double_t vkp[MaxHits];
  Int_t Kflag[MaxHits];
  Int_t Pflag[MaxHits];
  Int_t Heavyflag[MaxHits];

};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    TPCHid    = 100000,
  };
  Double_t tofaddjustment[24] = {0};

  Bool_t PionSelection(Double_t mass2, Double_t mom, Double_t nsigma)
  {
    if(TMath::IsNaN(mass2) || mass2<0) return false;

    Double_t pdgmass2 = pdg::PionMass()*pdg::PionMass(); //(GeV/c2)^2

    //Measured values(sigma of M2 spectrum)
    TF1 *f_sigmaM2 = new TF1("f_sigmaM2", eq_sigmaM2.Data(), 0., 5.);
    f_sigmaM2 -> FixParameter(0, 0.00427278);
    f_sigmaM2 -> FixParameter(1, -0.00916854);
    f_sigmaM2 -> FixParameter(2, 0.000130298);
    f_sigmaM2 -> FixParameter(3, pdgmass2);

    TF1 *f_M2 = new TF1("f_M2", eq_M2.Data(), 0., 1.5);
    f_M2 -> FixParameter(0, 0.0233145);
    f_M2 -> FixParameter(1, -0.0100331);
    f_M2 -> FixParameter(2, 0.013772);

    Double_t m2cut = nsigma*TMath::Sqrt(f_sigmaM2 -> Eval(mom)); //nsigma cut for M^2
    Double_t measured_m2 = f_M2 -> Eval(mom); //Measured M^2
    return (TMath::Abs(mass2 - measured_m2) < m2cut);
  }

  Bool_t KaonSelection(Double_t mass2, Double_t mom, Double_t nsigma)
  {
    if(TMath::IsNaN(mass2) || mass2<0) return false;
    if(mom > 1.4) return false;
    Double_t pdgmass2 = pdg::KaonMass()*pdg::KaonMass(); //(GeV/c2)^2

    //Measured values(sigma of M2 spectrum)
    TF1 *f_sigmaM2 = new TF1("f_sigmaM2", eq_sigmaM2.Data(), 0., 5.);
    f_sigmaM2 -> FixParameter(0, 0.000216333);
    f_sigmaM2 -> FixParameter(1, -0.00036971);
    f_sigmaM2 -> FixParameter(2, 0.000144059);
    f_sigmaM2 -> FixParameter(3, pdgmass2);

    TF1 *f_M2 = new TF1("f_M2", eq_M2.Data(), 0., 1.5);
    f_M2 -> FixParameter(0, 0.246792);
    f_M2 -> FixParameter(1, -0.0192738);
    f_M2 -> FixParameter(2, 0.0206932);

    Double_t m2cut = nsigma*TMath::Sqrt(f_sigmaM2 -> Eval(mom)); //nsigma cut for M^2
    Double_t measured_m2 = f_M2 -> Eval(mom); //Measured M^2
    return (TMath::Abs(mass2 - measured_m2) < m2cut);
  }

  Bool_t ProtonSelection(Double_t mass2, Double_t mom, Double_t nsigma)
  {
    if(TMath::IsNaN(mass2) || mass2<0) return false;

    Double_t pdgmass2 = pdg::ProtonMass()*pdg::ProtonMass(); //(GeV/c2)^2

    //Measured values(sigma of M2 spectrum)
    TF1 *f_sigmaM2 = new TF1("f_sigmaM2", eq_sigmaM2.Data(), 0., 5.);
    f_sigmaM2 -> FixParameter(0, 0.000184038);
    f_sigmaM2 -> FixParameter(1, 9.99371e-05);
    f_sigmaM2 -> FixParameter(2, 6.91217e-05);
    f_sigmaM2 -> FixParameter(3, pdgmass2);

    TF1 *f_M2 = new TF1("f_M2", eq_M2.Data(), 0., 1.5);
    f_M2 -> FixParameter(0, 0.988934);
    f_M2 -> FixParameter(1, -0.226358);
    f_M2 -> FixParameter(2, 0.113411);

    Double_t m2cut = nsigma*TMath::Sqrt(f_sigmaM2 -> Eval(mom)); //nsigma cut for M^2
    Double_t measured_m2 = f_M2 -> Eval(mom); //Measured M^2
    return (TMath::Abs(mass2 - measured_m2) < m2cut);
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
  static const auto MaxChisqrBcOut = gUser.GetParameter("MaxChisqrBcOut");
  static const auto MaxChisqrKurama = gUser.GetParameter("MaxChisqrKurama");
  static const auto KKEvent = gUser.GetParameter("KKEvent");
  static const auto KPEvent = gUser.GetParameter("KPEvent");
  static const auto KHeavyEvent = gUser.GetParameter("KHeavyEvent");
  static const Double_t ScatMomCut = gUser.GetParameter("ScatMomCut");

  static const auto ElectronMass = pdg::ElectronMass();
  static const auto PionMass = pdg::PionMass();
  static const auto KaonMass = pdg::KaonMass();
  static const auto K0Mass = pdg::K0Mass();  
  static const auto PhiMass = pdg::Mass(333);
  static const auto ProtonMass = pdg::ProtonMass();
  static const auto LambdaMass = pdg::LambdaMass();
  static const auto XiMinusMass = pdg::XiMinusMass();
  static const Double_t Carbon12Mass = 12.*TGeoUnit::amu_c2 - 6.*ElectronMass;
  static const Double_t Boron11Mass  = 11.009305167*TGeoUnit::amu_c2 - 5.*ElectronMass;
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};
  TVector3 tgtpos(0, 0, tpc::ZTarget);
  TVector3 qnan_vec = TVector3(qnan, qnan, qnan);  

  if( ievent%1000==0 ){
    //if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.nhHtof = src.nhHtof;
  for(Int_t it=0; it<event.nhHtof; it++){
    event.HtofSeg.push_back(src.HtofSeg[it]);
    event.tHtof.push_back(src.tHtof[it]);
    event.dtHtof.push_back(src.dtHtof[it]);
    event.deHtof.push_back(src.deHtof[it]);
    event.posHtof.push_back(src.posHtof[it]);
  }

  if(src.nKK != 1) return true;
  if(src.chisqrKurama[0] > MaxChisqrKurama || src.chisqrK18[0] > MaxChisqrBcOut) return true;
  if(KKEvent && src.Kflag[0] != 1){
    return true; //precut with Kurama tracking
  }
  if(KPEvent && src.Pflag[0] != 1){
    return true; //precut with Kurama tracking
  }
  if(KHeavyEvent && src.Heavyflag[0] != 1){
    return true; //precut with Kurama tracking
  }
  if(ScatMomCut>0 && src.pKurama[0]<ScatMomCut) return true;

  if(src.ntKurama != **src.ntTPCKurama)
    std::cerr << "Kurama Event Missmatching : DstTPCKuramaK18Tracking(" << **src.ntTPCKurama
	      << ") <-> DstKScat(" << src.ntKurama << ") " << std::endl;
  if(src.ntK18 != **src.ntTPCK18)
    std::cerr << "K18 Event Missmatching : DstTPCKuramaK18Tracking(" << **src.ntTPCK18
	      <<  ") <-> DstKScat(" << src.ntK18 << ") " << std::endl;


  // <-> DstKScat" << std::endl;

#if UsePidLH
  usePidLH = true;
#endif  
  HF1( 1, event.status++ );
  event.ntK18 = src.ntK18;
  for(int it=0; it<src.ntK18; ++it){
    event.pK18.push_back(src.pK18[it]);
    event.chisqrK18.push_back(src.chisqrK18[it]);
    event.xtgtK18.push_back(src.xtgtK18[it]);
    event.ytgtK18.push_back(src.ytgtK18[it]);
    event.utgtK18.push_back(src.utgtK18[it]);
    event.vtgtK18.push_back(src.vtgtK18[it]);
    event.thetaK18.push_back(src.thetaK18[it]);
  }
  event.xhtofK18 = **src.xhtofK18;
  event.yhtofK18 = **src.yhtofK18;
  event.tpcidTPCK18 = **src.tpcidTPCK18;
  event.isgoodTPCK18 = **src.isgoodTPCK18;
  event.chisqrTPCK18 = **src.chisqrTPCK18;
  event.pTPCK18 = **src.pTPCK18;
  event.qTPCK18 = **src.qTPCK18;
  event.xtgtTPCK18 = **src.xtgtTPCK18;
  event.ytgtTPCK18 = **src.ytgtTPCK18;
  event.utgtTPCK18 = **src.utgtTPCK18;
  event.vtgtTPCK18 = **src.vtgtTPCK18;
  event.thetaTPCK18 = **src.thetaTPCK18;
  event.lhtofTPCK18 = **src.lhtofTPCK18;
  event.xhtofTPCK18 = **src.xhtofTPCK18;
  event.yhtofTPCK18 = **src.yhtofTPCK18;
  event.lvpTPCK18 = **src.lvpTPCK18;
  event.xvpTPCK18 = **src.xvpTPCK18;
  event.yvpTPCK18 = **src.yvpTPCK18;
  event.p_3rd = **src.p_3rd;
  event.xoutK18 = **src.xoutK18;
  event.youtK18 = **src.youtK18;
  event.uoutK18 = **src.uoutK18;
  event.voutK18 = **src.voutK18;
  event.xtgtHS = **src.xtgtHS;
  event.ytgtHS = **src.ytgtHS;
  event.ztgtHS = **src.ztgtHS;
  event.xvpHS = **src.xvpHS;
  event.yvpHS = **src.yvpHS;
  event.zvpHS = **src.zvpHS;
  event.layerK18 = **src.layerK18;
  event.wireK18 = **src.wireK18;
  event.localhitposK18 = **src.localhitposK18;
  event.wposK18 = **src.wposK18;

  event.ntKurama = src.ntKurama;
  event.tpcidTPCKurama = **src.tpcidTPCKurama;
  event.isgoodTPCKurama = **src.isgoodTPCKurama;
  event.chisqrTPCKurama = **src.chisqrTPCKurama;
  event.pTPCKurama = **src.pTPCKurama;
  event.qTPCKurama = **src.qTPCKurama;
  event.xtgtTPCKurama = **src.xtgtTPCKurama;
  event.ytgtTPCKurama = **src.ytgtTPCKurama;
  event.utgtTPCKurama = **src.utgtTPCKurama;
  event.vtgtTPCKurama = **src.vtgtTPCKurama;
  event.thetaTPCKurama = **src.thetaTPCKurama;
  event.pathTPCKurama = **src.pathTPCKurama;
  event.lhtofTPCKurama = **src.lhtofTPCKurama;
  event.xhtofTPCKurama = **src.xhtofTPCKurama;
  event.yhtofTPCKurama = **src.yhtofTPCKurama;
  event.lgasvesselTPCKurama = **src.lgasvesselTPCKurama;
  event.xgasvesselTPCKurama = **src.xgasvesselTPCKurama;
  event.ygasvesselTPCKurama = **src.ygasvesselTPCKurama;
  event.lvpTPCKurama = **src.lvpTPCKurama;
  event.xvpTPCKurama = **src.xvpTPCKurama;
  event.yvpTPCKurama = **src.yvpTPCKurama;
  event.xhtofKurama = **src.xhtofKurama;
  event.yhtofKurama = **src.yhtofKurama;
  event.xvpKurama = **src.xvpKurama;
  event.yvpKurama = **src.yvpKurama;
  event.zvpKurama = **src.zvpKurama;
  event.layer = **src.layer;
  event.wire = **src.wire;
  event.localhitpos = **src.localhitpos;
  event.wpos = **src.wpos;

  event.insideTPC = **src.insideTPC;
  event.m2TPCKurama = **src.m2TPCKurama;
  
  //event.piflagTPCKurama.resize(src.ntKurama);
  event.kflagTPCKurama.resize(src.ntKurama);
  event.pflagTPCKurama.resize(src.ntKurama);
  const int nKurama = src.ntKurama;
  event.chisqrKurama.resize(nKurama);
  event.pKurama.resize(nKurama);
  event.qKurama.resize(nKurama);
  event.m2.resize(nKurama);
  event.m2Org.resize(nKurama);
  event.xtgtKurama.resize(nKurama);
  event.ytgtKurama.resize(nKurama);
  event.utgtKurama.resize(nKurama);
  event.vtgtKurama.resize(nKurama);
  event.thetaKurama.resize(nKurama);
  event.pathKurama.resize(nKurama);
  event.tofsegKurama.resize(nKurama);
  event.cstof.resize(nKurama);
  event.pathwcKurama.resize(nKurama);
  event.xin.resize(nKurama);
  event.yin.resize(nKurama);
  event.zin.resize(nKurama);
  event.pxin.resize(nKurama);
  event.pyin.resize(nKurama);
  event.pzin.resize(nKurama);
  event.xout.resize(nKurama);
  event.yout.resize(nKurama);
  event.zout.resize(nKurama);
  event.pxout.resize(nKurama);
  event.pyout.resize(nKurama);
  event.pzout.resize(nKurama);
  for(Int_t it=0; it<src.ntKurama; ++it){
    event.chisqrKurama[it] = src.chisqrKurama[it];
    event.pKurama[it] = src.pKurama[it];
    event.qKurama[it] = src.qKurama[it];
    event.m2[it] = src.m2[it];
    event.m2Org[it] = src.m2Org[it];
    event.xtgtKurama[it] = src.xtgtKurama[it];
    event.ytgtKurama[it] = src.ytgtKurama[it];
    event.utgtKurama[it] = src.utgtKurama[it];
    event.vtgtKurama[it] = src.vtgtKurama[it];
    event.thetaKurama[it] = src.thetaKurama[it];
    event.pathKurama[it] = src.pathKurama[it];
    event.tofsegKurama[it] = src.tofsegKurama[it];
    event.cstof[it] = src.cstof[it];
    std::cout << __LINE__ << " " << event.cstof[it] << std::endl;
    event.pathwcKurama[it] = src.pathwcKurama[it];
    event.xin[it] = src.xin[it];
    event.yin[it] = src.yin[it];
    event.zin[it] = src.zin[it];
    event.pxin[it] = src.pxin[it];
    event.pyin[it] = src.pyin[it];
    event.pzin[it] = src.pzin[it];
    event.xout[it] = src.xout[it];
    event.yout[it] = src.yout[it];
    event.zout[it] = src.zout[it];
    event.pxout[it] = src.pxout[it];
    event.pyout[it] = src.pyout[it];
    event.pzout[it] = src.pzout[it];
    int seg = src.tofsegKurama[it] - 1;
    if(event.cstof[it] + tofaddjustment[seg] > 0.)
      event.m2TPCKurama[it]
	= Kinematics::MassSquare(event.pTPCKurama[it],event.pathTPCKurama[it],event.cstof[it] + tofaddjustment[seg]);
    else event.m2TPCKurama[it] = TMath::QuietNaN();
    // event.piflagTPCKurama[it]
    //   = PionSelection(event.m2TPCKurama[it], event.pTPCKurama[it], nsigma_m2);
    event.kflagTPCKurama[it]
      = KaonSelection(event.m2TPCKurama[it], event.pTPCKurama[it], nsigma_m2);
    std::cout << __FILE__ << " " << __LINE__
	      << " m2TPCKurama:" << event.m2TPCKurama[it] << " event.pTPCKurama[it]:" << event.pTPCKurama[it]  << " kflag:" <<  event.kflagTPCKurama[it] << std::endl;

    event.pflagTPCKurama[it]
      = ProtonSelection(event.m2TPCKurama[it], event.pTPCKurama[it], nsigma_m2);
  }
  HF1( 1, event.status++ );

  event.vtxTPC = **src.vtxTPC;
  event.vtyTPC = **src.vtyTPC;
  event.vtzTPC = **src.vtzTPC;
  event.closeDistTPC = **src.closeDistTPC;

  event.nKm = src.nKm;
  event.nKp = src.nKp;
  event.nKK = src.nKK;
  for(Int_t it=0; it<src.nKK; ++it){
    event.vtx.push_back(src.vtx[it]);
    event.vty.push_back(src.vty[it]);
    event.vtz.push_back(src.vtz[it]);
    event.closeDist.push_back(src.closeDist[it]);
    event.inside.push_back(src.inside[it]);
    event.MissMass.push_back(src.MissMass[it]);
    event.MissMassCorr.push_back(src.MissMassCorr[it]);
    event.MissMassCorrDE.push_back(src.MissMassCorrDE[it]);
    event.pOrg.push_back(src.pOrg[it]);
    event.pCalc.push_back(src.pCalc[it]);
    event.pCorr.push_back(src.pCorr[it]);
    event.pCorrDE.push_back(src.pCorrDE[it]);
    event.xkm.push_back(src.xkm[it]);
    event.ykm.push_back(src.ykm[it]);
    event.ukm.push_back(src.ukm[it]);
    event.vkm.push_back(src.vkm[it]);
    event.xkp.push_back(src.xkp[it]);
    event.ykp.push_back(src.ykp[it]);
    event.ukp.push_back(src.ukp[it]);
    event.vkp.push_back(src.vkp[it]);
    event.Kflag.push_back(src.Kflag[it]);
    event.Pflag.push_back(src.Pflag[it]);
    event.Heavyflag.push_back(src.Heavyflag[it]);
  }
  event.isgoodTPC = **src.isgoodTPC;

  event.MissMassTPC = **src.MissMassTPC;
  event.MissMassCorrTPC = **src.MissMassCorrTPC;
  event.MissMassCorrDETPC = **src.MissMassCorrDETPC;
  event.MissMassNuclTPC = **src.MissMassNuclTPC;
  event.MissMassNuclCorrTPC = **src.MissMassNuclCorrTPC;
  event.MissMassNuclCorrDETPC = **src.MissMassNuclCorrDETPC;
  event.pOrgTPC = **src.pOrgTPC;
  event.pCorrTPC = **src.pCorrTPC;
  event.pCorrDETPC = **src.pCorrDETPC;
  event.pCalcTPC = **src.pCalcTPC;
  event.thetaCMTPC = **src.thetaCMTPC;
  event.costCMTPC = **src.costCMTPC;
  event.pCalcDETPC = **src.pCalcDETPC;
  event.thetaCMDETPC = **src.thetaCMDETPC;
  event.costCMDETPC = **src.costCMDETPC;
  event.xistarpCalcDETPC = **src.xistarpCalcDETPC;
  event.xistarthetaCMDETPC = **src.xistarthetaCMDETPC;
  event.xistarcostCMDETPC = **src.xistarcostCMDETPC;
  event.kpscatpCalcTPC = **src.kpscatpCalcTPC;
  event.kpscatthetaCMTPC = **src.kpscatthetaCMTPC;
  event.kpscatcostCMTPC = **src.kpscatcostCMTPC;
  event.kpscatpCalcDETPC = **src.kpscatpCalcDETPC;
  event.kpscatthetaCMDETPC = **src.kpscatthetaCMDETPC;
  event.kpscatcostCMDETPC = **src.kpscatcostCMDETPC;
  event.thetaTPC = **src.thetaTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

  double BE = 0.;
  for(Int_t iKK=0; iKK<event.nKK; iKK++){
    if(KPEvent){
      if(event.MissMassNuclCorrDETPC[0]<0.1) return true;
      //BE = event.MissMassNuclCorrDETPC[iKK] - KaonMass - Boron11Mass - 0.075;
      BE = event.MissMassNuclCorrDETPC[iKK] - KaonMass - Boron11Mass;
    }
    else if (KKEvent){
      BE = 0.; // choose as you like
    }
  }  
  
  event.xbTPC.resize(src.nKK);
  event.ybTPC.resize(src.nKK);
  event.ubTPC.resize(src.nKK);
  event.vbTPC.resize(src.nKK);
  event.xsTPC.resize(src.nKK);
  event.ysTPC.resize(src.nKK);
  event.usTPC.resize(src.nKK);
  event.vsTPC.resize(src.nKK);

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
			 **src.charge, **src.nhtrack, **src.helix_cx,
			 **src.helix_cy, **src.helix_z0, **src.helix_r,
			 **src.helix_dz, **src.hitlayer, **src.track_cluster_mrow,
			 **src.helix_t, **src.track_cluster_de, **src.resolution_x,
			 **src.resolution_y, **src.resolution_z, **src.hitpos_x,
			 **src.hitpos_y, **src.hitpos_z);
  
  for(Int_t idScat=0; idScat<event.ntKurama; ++idScat){
    for(Int_t idKm=0; idKm<event.ntK18; ++idKm){
      Int_t id = idScat*src.ntK18 + idKm;
      Int_t inside = event.inside[id];
      Double_t us = event.ukp[id];
      Double_t vs = event.vkp[id];
      Double_t closeDist = event.closeDist[id];
      Double_t kkvertx = event.vtx[id];
      Double_t kkverty = event.vty[id];
      Double_t kkvertz = event.vtz[id];

      Double_t KaonMom = event.pCalc[id];
      Double_t pScat = event.pOrg[id];
      Double_t pScatCorr = event.pCorr[id];
      Double_t pScatCorrDE = event.pCorrDE[id];
      Double_t MissMass = event.MissMass[id];
      Double_t MissMassCorr = event.MissMassCorr[id];
      Double_t MissMassCorrDE = event.MissMassCorrDE[id];
      Double_t MissMassNuclTPC = event.MissMassNuclTPC[id];
      Double_t MissMassNuclCorrTPC = event.MissMassNuclCorrTPC[id];
      Double_t MissMassNuclCorrDETPC = event.MissMassNuclCorrDETPC[id];
      Double_t thetaTPC = event.thetaTPC[id];

      if(event.chisqrK18[idKm] < MaxChisqrBcOut && event.chisqrKurama[idScat] < MaxChisqrKurama){
	HF1(12, event.isgoodTPCK18[idKm]);
	HF1(13, event.isgoodTPCKurama[idScat]);
	HF1(14, event.isgoodTPC[id]);
	if(inside==1){
	  HF1(22, event.isgoodTPCK18[idKm]);
	  HF1(23, event.isgoodTPCKurama[idScat]);
	  HF1(24, event.isgoodTPC[id]);
	  if(event.qKurama[idScat] > 0 &&
	     event.m2[idScat] > 0.12 && event.m2[idScat] < 0.3){
	    HF1(32, event.isgoodTPCK18[idKm]);
	    HF1(33, event.isgoodTPCKurama[idScat]);
	    HF1(34, event.isgoodTPC[id]);
	    if(event.pKurama[idScat] < 1.4 && event.pKurama[idScat] > 1.1){
	      HF1(42, event.isgoodTPCK18[idKm]);
	      HF1(43, event.isgoodTPCKurama[idScat]);
	      HF1(44, event.isgoodTPC[id]);
	    }
	  }
	}

	HF1(1001, event.pK18[idKm]);
	HF1(1002, event.pKurama[idScat]);
	HF1(1003, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]));
	HF2(1004, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]), event.pKurama[idScat]);
	if(event.qKurama[idScat] > 0){
	  HF1(1005, TMath::Sqrt(event.m2Org[idScat]));
	  HF1(1006, TMath::Sqrt(event.m2[idScat]));
	}
	HF2(1010, kkvertz + tpc::ZTarget, kkvertx);
	HF1(1011, closeDist);
	HF1(1012, kkvertx);
	HF1(1013, kkverty);
	HF1(1014, kkvertz);
	if(event.insideTPC[id]==1){
	  HF1(2001, event.pK18[idKm]);
	  HF1(2002, event.pKurama[idScat]);
	  HF1(2003, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]));
	  HF2(2004, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]), event.pKurama[idScat]);
	  HF2(2010, kkvertz + tpc::ZTarget, kkvertx);
	  HF1(2011, closeDist);
	  HF1(2012, kkvertx);
	  HF1(2013, kkverty);
	  HF1(2014, kkvertz);

	  //K+
	  if(event.qTPCKurama[idScat] > 0 && event.pTPCKurama[idScat] < 1.4 &&
	     event.m2TPCKurama[idScat] > 0.15 && event.m2TPCKurama[idScat] < 0.35){
	    HF1(3001, event.pK18[idKm]);
	    HF1(3002, event.pKurama[idScat]);
	    HF1(3003, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]));
	    HF2(3004, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]), event.pKurama[idScat]);
	    HF1(3010, closeDist);
	    HF1(3011, kkvertx);
	    HF1(3012, kkverty);
	    HF1(3013, kkvertz);
	    HF1(3014, MissMassCorr);
	    HF1(3015, MissMassCorrDE);
	    HF2(3016, us, pScat - KaonMom);
	    HF2(3017, us, pScatCorr - KaonMom);
	    HF2(3018, us, pScatCorrDE - KaonMom);
	    HF2(3019, us, MissMass);
	    HF2(3020, us, MissMassCorr);
	    HF2(3021, us, MissMassCorrDE);
	    HF2(3022, vs, pScat - KaonMom);
	    HF2(3023, vs, pScatCorr - KaonMom);
	    HF2(3024, vs, pScatCorrDE - KaonMom);
	    HF2(3025, vs, MissMass);
	    HF2(3026, vs, MissMassCorr);
	    HF2(3027, vs, MissMassCorrDE);
	    HF2(3028, KaonMom, pScat - KaonMom);
	    HF2(3029, KaonMom, pScatCorr - KaonMom);
	    HF2(3030, KaonMom, pScatCorrDE - KaonMom);

	    HFProf(3116, us, pScat - KaonMom);
	    HFProf(3117, us, pScatCorr - KaonMom);
	    HFProf(3118, us, pScatCorrDE - KaonMom);
	    HFProf(3119, us, MissMass);
	    HFProf(3120, us, MissMassCorr);
	    HFProf(3121, us, MissMassCorrDE);
	    HFProf(3122, vs, pScat - KaonMom);
	    HFProf(3123, vs, pScatCorr - KaonMom);
	    HFProf(3124, vs, pScatCorrDE - KaonMom);
	    HFProf(3125, vs, MissMass);
	    HFProf(3126, vs, MissMassCorr);
	    HFProf(3127, vs, MissMassCorrDE);
	    HFProf(3128, KaonMom, pScat - KaonMom);
	    HFProf(3129, KaonMom, pScatCorr - KaonMom);
	    HFProf(3130, KaonMom, pScatCorrDE - KaonMom);

	    if(event.pTPCKurama[idScat] > 1.1){
	      HF1(4001, event.pK18[idKm]);
	      HF1(4002, event.pKurama[idScat]);
	      HF1(4003, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]));
	      HF2(4004, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]), event.pKurama[idScat]);
	      HF1(4010, closeDist);
	      HF1(4011, kkvertx);
	      HF1(4012, kkverty);
	      HF1(4013, kkvertz);
	      HF1(4014, MissMassCorr);
	      HF1(4015, MissMassCorrDE);
	      HF2(4016, us, pScat - KaonMom);
	      HF2(4017, us, pScatCorr - KaonMom);
	      HF2(4018, us, pScatCorrDE - KaonMom);
	      HF2(4019, us, MissMass);
	      HF2(4020, us, MissMassCorr);
	      HF2(4021, us, MissMassCorrDE);
	      HF2(4022, vs, pScat - KaonMom);
	      HF2(4023, vs, pScatCorr - KaonMom);
	      HF2(4024, vs, pScatCorrDE - KaonMom);
	      HF2(4025, vs, MissMass);
	      HF2(4026, vs, MissMassCorr);
	      HF2(4027, vs, MissMassCorrDE);
	      HF2(4028, KaonMom, pScat - KaonMom);
	      HF2(4029, KaonMom, pScatCorr - KaonMom);
	      HF2(4030, KaonMom, pScatCorrDE - KaonMom);

	      HFProf(4116, us, pScat - KaonMom);
	      HFProf(4117, us, pScatCorr - KaonMom);
	      HFProf(4118, us, pScatCorrDE - KaonMom);
	      HFProf(4119, us, MissMass);
	      HFProf(4120, us, MissMassCorr);
	      HFProf(4121, us, MissMassCorrDE);
	      HFProf(4122, vs, pScat - KaonMom);
	      HFProf(4123, vs, pScatCorr - KaonMom);
	      HFProf(4124, vs, pScatCorrDE - KaonMom);
	      HFProf(4125, vs, MissMass);
	      HFProf(4126, vs, MissMassCorr);
	      HFProf(4127, vs, MissMassCorrDE);
	      HFProf(4128, KaonMom, pScat - KaonMom);
	      HFProf(4129, KaonMom, pScatCorr - KaonMom);
	      HFProf(4130, KaonMom, pScatCorrDE - KaonMom);
	    }
	  }

	  //Proton
	  if(event.qTPCKurama[idScat] > 0 &&
	     event.m2TPCKurama[idScat] > 0.5 && event.m2TPCKurama[idScat] < 1.5){

	    HF1(4201, event.pK18[idKm]);
	    HF1(4202, event.pKurama[idScat]);
	    HF1(4203, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]));
	    HF2(4204, event.qKurama[idScat]*TMath::Sqrt(event.m2[idScat]), event.pKurama[idScat]);
	    HF1(4210, closeDist);
	    HF1(4211, kkvertx);
	    HF1(4212, kkverty);
	    HF1(4213, kkvertz);
	    HF1(4214, MissMassCorr);
	    HF1(4215, MissMassCorrDE);
	    //HF2(4216, us, pScat - ProtonMom);
	    //HF2(4217, us, pScatCorr - ProtonMom);
	    //HF2(4218, us, pScatCorrDE - ProtonMom);
	    HF2(4219, us, MissMass);
	    HF2(4220, us, MissMassCorr);
	    HF2(4221, us, MissMassCorrDE);
	    //HF2(4222, vs, pScat - ProtonMom);
	    //HF2(4223, vs, pScatCorr - ProtonMom);
	    //HF2(4224, vs, pScatCorrDE - ProtonMom);
	    HF2(4225, vs, MissMass);
	    HF2(4226, vs, MissMassCorr);
	    HF2(4227, vs, MissMassCorrDE);

	    //HFProf(4316, us, pScat - ProtonMom);
	    //HFProf(4317, us, pScatCorr - ProtonMom);
	    //HFProf(4318, us, pScatCorrDE - ProtonMom);
	    HFProf(4319, us, MissMass);
	    HFProf(4320, us, MissMassCorr);
	    HFProf(4321, us, MissMassCorrDE);
	    //HFProf(4322, vs, pScat - ProtonMom);
	    //HFProf(4323, vs, pScatCorr - ProtonMom);
	    //HFProf(4324, vs, pScatCorrDE - ProtonMom);
	    HFProf(4325, vs, MissMass);
	    HFProf(4326, vs, MissMassCorr);
	    HFProf(4327, vs, MissMassCorrDE);
	    // HF1( 101, MissMassNucl );
	    // HF2( 102, MissMassNuclCorr );
	    // HF2( 103, MissMassNuclCorrDE );
	    HF1( 104, MissMassNuclTPC );
	    HF1( 105, MissMassNuclCorrTPC );
	    HF1( 106, MissMassNuclCorrDETPC );
	    if(thetaTPC>3.5 && thetaTPC<4.5){
	      HF1(106, MissMassNuclTPC);
	      // std::cout << "MissMassNuclTPC: " << MissMassNuclTPC << std::endl;
	      HF1(107, MissMassNuclCorrTPC);
	      // std::cout << "MissMassNuclCorrTPC: " << MissMassNuclCorrTPC << std::endl;
	      HF1(108, MissMassNuclCorrDETPC);
	      // std::cout << "MissMassNuclCorrDETPC: " << MissMassNuclCorrDETPC << std::endl;
	    }
	  }
	}
      }
    }
  }

  for(Int_t idScat=0; idScat<event.ntKurama; ++idScat){
    if(event.isgoodTPCKurama[idScat]!=1) continue;
    for(Int_t idKm=0; idKm<event.ntK18; ++idKm){
      Int_t id = idScat*src.ntK18 + idKm;
      
      /*
	if(event.isgoodTPC[idScat]!=1) continue;
	event.xbTPC[id] = event.xtgtTPCKurama[idKm];
	event.ybTPC[id] = event.ytgtTPCKurama[idKm];
	event.ubTPC[id] = event.utgtTPCKurama[idKm];
	event.vbTPC[id] = event.vtgtTPCKurama[idKm];
      */

      //Temporary
      event.xbTPC[id] = src.xkm[idKm];
      event.ybTPC[id] = src.ykm[idKm];
      event.ubTPC[id] = src.ukm[idKm];
      event.vbTPC[id] = src.vkm[idKm];

      event.xsTPC[id] = event.xtgtTPCKurama[idScat];
      event.ysTPC[id] = event.ytgtTPCKurama[idScat];
      event.usTPC[id] = event.utgtTPCKurama[idScat];
      event.vsTPC[id] = event.vtgtTPCKurama[idScat];

      Int_t inside = event.insideTPC[id];
      Double_t us = event.usTPC[id];
      Double_t vs = event.vsTPC[id];
      Double_t closeDist = event.closeDistTPC[id];
      Double_t kkvertx = event.vtxTPC[id];
      Double_t kkverty = event.vtyTPC[id];
      Double_t kkvertz = event.vtzTPC[id];
      Double_t KaonMom = event.pCalcDETPC[id];
      Double_t KaonMomCorrDE = event.pCalcDETPC[id];
      Double_t pScat = event.pOrgTPC[id];
      Double_t pScatCorr = event.pCorrTPC[id];
      Double_t pScatCorrDE = event.pCorrDETPC[id];
      Double_t MissMass = event.MissMassTPC[id];
      Double_t MissMassCorr = event.MissMassCorrTPC[id];
      Double_t MissMassCorrDE = event.MissMassCorrDETPC[id];
      Double_t MissMassNucl = event.MissMassNuclTPC[id];
      Double_t MissMassNuclCorr = event.MissMassNuclCorrTPC[id];
      Double_t MissMassNuclCorrDE = event.MissMassNuclCorrDETPC[id];
      Double_t XiStarKaonMomCorrDE = event.xistarpCalcDETPC[id];
      Double_t ProtonMom = event.kpscatpCalcDETPC[id];
      Double_t ProtonMomCorrDE = event.kpscatpCalcDETPC[id];

      if(event.chisqrK18[idKm] < MaxChisqrBcOut && event.chisqrKurama[idScat] < MaxChisqrKurama){
	HF1(5001, event.pTPCK18[idKm]);
	HF1(5002, event.pTPCKurama[idScat]);
	HF1(5003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	HF2(5004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	HF1(5007, event.closeDist[id]);
	HF2(5010, kkvertz + tpc::ZTarget, kkvertx);
	HF1(5011, closeDist);
	HF1(5012, kkvertx);
	HF1(5013, kkverty);
	HF1(5014, kkvertz);
	if(inside==1){
	  HF1(6001, event.pTPCK18[idKm]);
	  HF1(6002, event.pTPCKurama[idScat]);
	  HF1(6003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	  HF2(6004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	  HF2(6010, kkvertz + tpc::ZTarget, kkvertx);
	  HF1(6011, closeDist);
	  HF1(6012, kkvertx);
	  HF1(6013, kkverty);
	  HF1(6014, kkvertz);
	  HF2(6041, src.xkm[idKm], event.xbTPC[id]);
	  HF2(6042, src.ykm[idKm], event.ybTPC[id]);
	  HF2(6043, src.ukm[idKm], event.ubTPC[id]);
	  HF2(6044, src.vkm[idKm], event.vbTPC[id]);
	  HF2(6045, src.xkp[idKm], event.xsTPC[id]);
	  HF2(6046, src.ykp[idKm], event.ysTPC[id]);
	  HF2(6047, src.ukp[idKm], event.usTPC[id]);
	  HF2(6048, src.vkp[idKm], event.vsTPC[id]);
	  //K+
	  if(event.qTPCKurama[idScat] > 0 && event.pTPCKurama[idScat] < 1.4 &&
	     event.m2TPCKurama[idScat] > 0.15 && event.m2TPCKurama[idScat] < 0.35){
	    HF1(7001, event.pTPCK18[idKm]);
	    HF1(7002, event.pTPCKurama[idScat]);
	    HF1(7003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	    HF2(7004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	    HF1(7010, closeDist);
	    HF1(7011, kkvertx);
	    HF1(7012, kkverty);
	    HF1(7013, kkvertz);
	    HF1(7014, MissMassCorr);
	    HF1(7015, MissMassCorrDE);
	    HF2(7016, us, pScat - KaonMom);
	    HF2(7017, us, pScatCorr - KaonMom);
	    HF2(7018, us, pScatCorrDE - KaonMomCorrDE);
	    HF2(7019, us, MissMass);
	    HF2(7020, us, MissMassCorr);
	    HF2(7021, us, MissMassCorrDE);
	    HF2(7022, vs, pScat - KaonMom);
	    HF2(7023, vs, pScatCorr - KaonMom);
	    HF2(7024, vs, pScatCorrDE - KaonMomCorrDE);
	    HF2(7025, vs, MissMass);
	    HF2(7026, vs, MissMassCorr);
	    HF2(7027, vs, MissMassCorrDE);
	    HF2(7028, KaonMom, pScat - KaonMom);
	    HF2(7029, KaonMom, pScatCorr - KaonMom);
	    HF2(7030, KaonMomCorrDE, pScatCorrDE - KaonMomCorrDE);
	    HF2(7031, KaonMom, MissMass);
	    HF2(7032, KaonMom, MissMassCorr);
	    HF2(7033, KaonMomCorrDE, MissMassCorrDE);
	    HF2(7036, XiStarKaonMomCorrDE, pScatCorrDE - XiStarKaonMomCorrDE);
	    HF2(7039, XiStarKaonMomCorrDE, MissMassCorrDE);

	    HF2(7041, src.xkm[idKm], event.xbTPC[id]);
	    HF2(7042, src.ykm[idKm], event.ybTPC[id]);
	    HF2(7043, src.ukm[idKm], event.ubTPC[id]);
	    HF2(7044, src.vkm[idKm], event.vbTPC[id]);
	    HF2(7045, src.xkp[idKm], event.xsTPC[id]);
	    HF2(7046, src.ykp[idKm], event.ysTPC[id]);
	    HF2(7047, src.ukp[idKm], event.usTPC[id]);
	    HF2(7048, src.vkp[idKm], event.vsTPC[id]);

	    HFProf(7116, us, pScat - KaonMom);
	    HFProf(7117, us, pScatCorr - KaonMom);
	    HFProf(7118, us, pScatCorrDE - KaonMomCorrDE);
	    HFProf(7119, us, MissMass);
	    HFProf(7120, us, MissMassCorr);
	    HFProf(7121, us, MissMassCorrDE);
	    HFProf(7122, vs, pScat - KaonMom);
	    HFProf(7123, vs, pScatCorr - KaonMom);
	    HFProf(7124, vs, pScatCorrDE - KaonMomCorrDE);
	    HFProf(7125, vs, MissMass);
	    HFProf(7126, vs, MissMassCorr);
	    HFProf(7127, vs, MissMassCorrDE);
	    HFProf(7128, KaonMom, pScat - KaonMom);
	    HFProf(7129, KaonMom, pScatCorr - KaonMom);
	    HFProf(7130, KaonMomCorrDE, pScatCorrDE - KaonMomCorrDE);
	    HFProf(7131, KaonMom, MissMass);
	    HFProf(7132, KaonMom, MissMassCorr);
	    HFProf(7133, KaonMomCorrDE, MissMassCorrDE);
	    HFProf(7136, XiStarKaonMomCorrDE, pScatCorrDE - XiStarKaonMomCorrDE);
	    HFProf(7139, XiStarKaonMomCorrDE, MissMassCorrDE);

	    if(event.pTPCKurama[idScat] > 1.1){
	      HF1(8001, event.pTPCK18[idKm]);
	      HF1(8002, event.pTPCKurama[idScat]);
	      HF1(8003, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	      HF2(8004, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	      HF1(8010, closeDist);
	      HF1(8011, kkvertx);
	      HF1(8012, kkverty);
	      HF1(8013, kkvertz);
	      HF1(8014, MissMassCorr);
	      HF1(8015, MissMassCorrDE);
	      HF2(8016, us, pScat - KaonMom);
	      HF2(8017, us, pScatCorr - KaonMom);
	      HF2(8018, us, pScatCorrDE - KaonMomCorrDE);
	      HF2(8019, us, MissMass);
	      HF2(8020, us, MissMassCorr);
	      HF2(8021, us, MissMassCorrDE);
	      HF2(8022, vs, pScat - KaonMom);
	      HF2(8023, vs, pScatCorr - KaonMom);
	      HF2(8024, vs, pScatCorrDE - KaonMomCorrDE);
	      HF2(8025, vs, MissMass);
	      HF2(8026, vs, MissMassCorr);
	      HF2(8027, vs, MissMassCorrDE);
	      HF2(8028, KaonMom, pScat - KaonMom);
	      HF2(8029, KaonMom, pScatCorr - KaonMom);
	      HF2(8030, KaonMomCorrDE, pScatCorrDE - KaonMomCorrDE);
	      HF2(8031, KaonMom, MissMass);
	      HF2(8032, KaonMom, MissMassCorr);
	      HF2(8033, KaonMomCorrDE, MissMassCorrDE);

	      HFProf(8116, us, pScat - KaonMom);
	      HFProf(8117, us, pScatCorr - KaonMom);
	      HFProf(8118, us, pScatCorrDE - KaonMomCorrDE);
	      HFProf(8119, us, MissMass);
	      HFProf(8120, us, MissMassCorr);
	      HFProf(8121, us, MissMassCorrDE);
	      HFProf(8122, vs, pScat - KaonMom);
	      HFProf(8123, vs, pScatCorr - KaonMom);
	      HFProf(8124, vs, pScatCorrDE - KaonMomCorrDE);
	      HFProf(8125, vs, MissMass);
	      HFProf(8126, vs, MissMassCorr);
	      HFProf(8127, vs, MissMassCorrDE);
	      HFProf(8128, KaonMom, pScat - KaonMom);
	      HFProf(8129, KaonMom, pScatCorr - KaonMom);
	      HFProf(8130, KaonMomCorrDE, pScatCorrDE - KaonMomCorrDE);
	      HFProf(8131, KaonMom, MissMass);
	      HFProf(8132, KaonMom, MissMassCorr);
	      HFProf(8133, KaonMomCorrDE, MissMassCorrDE);

	      HF2(8041, src.xkm[idKm], event.xbTPC[id]);
	      HF2(8042, src.ykm[idKm], event.ybTPC[id]);
	      HF2(8043, src.ukm[idKm], event.ubTPC[id]);
	      HF2(8044, src.vkm[idKm], event.vbTPC[id]);
	      HF2(8045, src.xkp[idKm], event.xsTPC[id]);
	      HF2(8046, src.ykp[idKm], event.ysTPC[id]);
	      HF2(8047, src.ukp[idKm], event.usTPC[id]);
	      HF2(8048, src.vkp[idKm], event.vsTPC[id]);
	    }
	    else{
	      HF2(8036, XiStarKaonMomCorrDE, pScatCorrDE - XiStarKaonMomCorrDE);
	      HF2(8039, XiStarKaonMomCorrDE, MissMassCorrDE);

	      HFProf(8136, XiStarKaonMomCorrDE, pScatCorrDE - XiStarKaonMomCorrDE);
	      HFProf(8139, XiStarKaonMomCorrDE, MissMassCorrDE);
	    }
	  }

	  //Proton
	  if(event.qTPCKurama[idScat] > 0 &&
	     event.m2TPCKurama[idScat] > 0.5 && event.m2TPCKurama[idScat] < 1.5){
	    HF1(8201, event.pTPCK18[idKm]);
	    HF1(8202, event.pTPCKurama[idScat]);
	    HF1(8203, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]));
	    HF2(8204, event.qTPCKurama[idScat]*TMath::Sqrt(event.m2TPCKurama[idScat]), event.pTPCKurama[idScat]);
	    HF1(8210, closeDist);
	    HF1(8211, kkvertx);
	    HF1(8212, kkverty);
	    HF1(8213, kkvertz);
	    HF1(8214, MissMassCorr);
	    HF1(8215, MissMassCorrDE);
	    HF2(8216, us, pScat - ProtonMom);
	    HF2(8217, us, pScatCorr - ProtonMom);
	    HF2(8218, us, pScatCorrDE - ProtonMomCorrDE);
	    HF2(8219, us, MissMass);
	    HF2(8220, us, MissMassCorr);
	    HF2(8221, us, MissMassCorrDE);
	    HF2(8222, vs, pScat - ProtonMom);
	    HF2(8223, vs, pScatCorr - ProtonMom);
	    HF2(8224, vs, pScatCorrDE - ProtonMomCorrDE);
	    HF2(8225, vs, MissMass);
	    HF2(8226, vs, MissMassCorr);
	    HF2(8227, vs, MissMassCorrDE);
	    HF2(8228, ProtonMom, pScat - ProtonMom);
	    HF2(8229, ProtonMom, pScatCorr - ProtonMom);
	    HF2(8230, ProtonMomCorrDE, pScatCorrDE - ProtonMomCorrDE);
	    HF2(8231, ProtonMom, MissMass);
	    HF2(8232, ProtonMom, MissMassCorr);
	    HF2(8233, ProtonMomCorrDE, MissMassCorrDE);

	    HFProf(8316, us, pScat - ProtonMom);
	    HFProf(8317, us, pScatCorr - ProtonMom);
	    HFProf(8318, us, pScatCorrDE - ProtonMomCorrDE);
	    HFProf(8319, us, MissMass);
	    HFProf(8320, us, MissMassCorr);
	    HFProf(8321, us, MissMassCorrDE);
	    HFProf(8322, vs, pScat - ProtonMom);
	    HFProf(8323, vs, pScatCorr - ProtonMom);
	    HFProf(8324, vs, pScatCorrDE - ProtonMomCorrDE);
	    HFProf(8325, vs, MissMass);
	    HFProf(8326, vs, MissMassCorr);
	    HFProf(8327, vs, MissMassCorrDE);
	    HFProf(8328, ProtonMom, pScat - ProtonMom);
	    HFProf(8329, ProtonMom, pScatCorr - ProtonMom);
	    HFProf(8330, ProtonMomCorrDE, pScatCorrDE - ProtonMomCorrDE);
	    HFProf(8331, ProtonMom, MissMass);
	    HFProf(8332, ProtonMom, MissMassCorr);
	    HFProf(8333, ProtonMomCorrDE, MissMassCorrDE);

	    HF2(8241, src.xkm[idKm], event.xbTPC[id]);
	    HF2(8242, src.ykm[idKm], event.ybTPC[id]);
	    HF2(8243, src.ukm[idKm], event.ubTPC[id]);
	    HF2(8244, src.vkm[idKm], event.vbTPC[id]);
	    HF2(8245, src.xkp[idKm], event.xsTPC[id]);
	    HF2(8246, src.ykp[idKm], event.ysTPC[id]);
	    HF2(8247, src.ukp[idKm], event.usTPC[id]);
	    HF2(8248, src.vkp[idKm], event.vsTPC[id]);
	  }
	}
      }
    }
  }

  HF1( 1, event.status++ );

  for(Int_t it=0; it<src.nKK; ++it){
    Int_t inside = 0;
    if(event.inside[it]==1 && event.isgoodTPCKurama[it]==1){
      if(TMath::Abs(event.vtzTPC[it]) < 150.) inside = 1;
    }
    event.insideTPC[it] = inside;
  }

  event.nhTpc = **src.nhTpc;
  //if( **src.nhTpc == 0 ) return true;

  event.nclTpc = **src.nclTpc;
  event.remain_nclTpc = **src.remain_nclTpc;
#if RawCluster
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
#endif
  if( **src.ntTpc == 0 )
    return true;
  Int_t ntTpc = **src.ntTpc;
  event.ntTpc = ntTpc;
  event.ntKuramaCandidate = **src.ntKuramaCandidate;
  event.isKuramaCandidate = **src.isKuramaCandidate;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isBeam = **src.isBeam;
  event.isXi = **src.isXi;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.isMultiloop = **src.isMultiloop;
  event.charge = **src.charge;
  event.pid = **src.pid;
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

    for(int ih=0; ih<event.nhHtof; ++ih){
      HF2(20, event.mom0[it]*event.charge[it], event.path[it]/event.tHtof[ih]/MathTools::C());
    }
  }
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
  event.pathhit = **src.pathhit;
  event.alpha = **src.alpha;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_size = **src.track_cluster_size;
  event.track_cluster_mrow = **src.track_cluster_mrow;
  event.track_cluster_de_center = **src.track_cluster_de_center;
  event.track_cluster_x_center = **src.track_cluster_x_center;
  event.track_cluster_y_center = **src.track_cluster_y_center;
  event.track_cluster_z_center = **src.track_cluster_z_center;
  event.track_cluster_row_center = **src.track_cluster_row_center;

  event.chargeIndistinguishable = **src.chargeIndistinguishable;
  event.pid_inverted = **src.pid_inverted;
  event.chisqr_inverted = **src.chisqr_inverted;
  event.pval_inverted = **src.pval_inverted;
  event.helix_cx_inverted = **src.helix_cx_inverted;
  event.helix_cy_inverted = **src.helix_cy_inverted;
  event.helix_z0_inverted = **src.helix_z0_inverted;
  event.helix_r_inverted = **src.helix_r_inverted;
  event.helix_dz_inverted = **src.helix_dz_inverted;
  event.mom0_inverted = **src.mom0_inverted;

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

  event.isLambda = **src.isLambda;
  event.ncombiLambda = **src.ncombiLambda;
  event.distLambda = **src.distLambda;
  event.angleLambda = **src.angleLambda;
  event.bestmassLambda = **src.bestmassLambda;
  event.massLambda = **src.massLambda;
  event.vtxLambda_x = **src.vtxLambda_x;
  event.vtxLambda_y = **src.vtxLambda_y;
  event.vtxLambda_z = **src.vtxLambda_z;
  event.momLambda = **src.momLambda;
  event.momLambda_x = **src.momLambda_x;
  event.momLambda_y = **src.momLambda_y;
  event.momLambda_z = **src.momLambda_z;
  event.decaysidLambda = **src.decaysidLambda;
  event.decaysmomLambda = **src.decaysmomLambda;
  event.decaysmomLambda_x = **src.decaysmomLambda_x;
  event.decaysmomLambda_y = **src.decaysmomLambda_y;
  event.decaysmomLambda_z = **src.decaysmomLambda_z;

  event.nvtxTpcClustered = **src.nvtxTpcClustered;
  event.clusteredVtx_x = **src.clusteredVtx_x;
  event.clusteredVtx_y = **src.clusteredVtx_y;
  event.clusteredVtx_z = **src.clusteredVtx_z;
  event.clusteredVtxid = **src.clusteredVtxid;

  // Make PidPdf for LikelihoodPid
  //all
  event.hitidHtof.resize(event.ntTpc);
  event.tracklenHtof.resize(event.ntTpc);
  event.m2HtofVtx.resize(event.ntTpc);
  event.m2Htof.resize(event.ntTpc);
  event.extrapTgt.resize(event.ntTpc);
  event.extrapVtx.resize(event.ntTpc);
  event.tofHtof.resize(event.ntTpc);
  event.segHtof.resize(event.ntTpc);  
  event.invbetaHtof.resize(event.ntTpc);
  event.tracklenHtof.resize(event.ntTpc);
  event.posHtof_x.resize(event.ntTpc);
  event.posHtof_y.resize(event.ntTpc);
  event.posHtof_z.resize(event.ntTpc);
  event.insideTgt.resize(event.ntTpc);
  event.nsigma_tritonHtof.resize(event.ntTpc);
  event.nsigma_deutronHtof.resize(event.ntTpc);
  event.nsigma_protonHtof.resize(event.ntTpc);
  event.nsigma_kaonHtof.resize(event.ntTpc);
  event.nsigma_pionHtof.resize(event.ntTpc);
  event.nsigma_electronHtof.resize(event.ntTpc);      
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  Int_t ntrack_intarget = 0;
  Double_t x0[100] = {0};
  Double_t y0[100] = {0};
  Double_t u0[100] = {0};
  Double_t v0[100] = {0};
  for( Int_t itTpc=0; itTpc<event.ntTpc; ++itTpc ){
    if(event.isBeam[itTpc] || event.isK18[itTpc] || event.isAccidental[itTpc]) continue;
    double m2Htof    = qnan;
    double nsigma_t  = event.nsigma_triton[itTpc]   ;
    double nsigma_d  = event.nsigma_deutron[itTpc]  ;
    double nsigma_p  = event.nsigma_proton[itTpc]   ;
    double nsigma_k  = event.nsigma_kaon[itTpc]     ;
    double nsigma_pi = event.nsigma_pion[itTpc]     ;
    double nsigma_e  = event.nsigma_electron[itTpc] ;
    
    // double nsigmaHtof_t  = event.nsigma_tritonHtof[itTpc]   ;
    // double nsigmaHtof_d  = event.nsigma_deutronHtof[itTpc]  ;
    // double nsigmaHtof_p  = event.nsigma_protonHtof[itTpc]   ;
    // double nsigmaHtof_k  = event.nsigma_kaonHtof[itTpc]     ;
    // double nsigmaHtof_pi = event.nsigma_pionHtof[itTpc]     ;
    // double nsigmaHtof_e  = event.nsigma_electronHtof[itTpc] ;
    TPCLocalTrackHelix* track = TPCAna.GetTrackTPCHelix(itTpc);
    if (!track) continue;
    int htof_seg= -1;
    double tracklen_htof = -1.0;
    TVector3 pos_htof;
    bool isInsideTarget = false;
    // if(event.isKurama[itTpc]==1){
    //   for(int jgf=0; jgf<GFntTpc; jgf++){
    // 	if(itTpc==jgf) continue;
    // 	if( event.charge[jgf]==1 ) continue;
    // 	if( event.isBeam[jgf]==1
    // 	    || event.isK18[jgf]==1
    // 	    || event.isAccidental[jgf]==1 ) continue;
    // 	double extrapKurama = qnan;
    // 	double extrapOther = qnan;
    // 	double dist = qnan;
    // 	TVector3 momKurama; TVector3 momOther;	  
    // 	TVector3 vertex;

    // 	double thetafac = 0.3;
    // 	Double_t kurama_par[5];
    // 	kurama_par[0] = event.helix_cx[itTpc];
    // 	kurama_par[1] = event.helix_cy[itTpc];
    // 	kurama_par[2] = event.helix_z0[itTpc];
    // 	kurama_par[3] = event.helix_r[itTpc];
    // 	kurama_par[4] = event.helix_dz[itTpc];
    // 	Int_t kurama_nh = event.helix_t[itTpc].size();
    // 	Double_t kurama_theta_range = event.helix_t[itTpc][kurama_nh-1] - event.helix_t[itTpc][0];
    // 	Double_t kurama_theta_min = event.helix_t[itTpc][0] - thetafac*kurama_theta_range;
    // 	Double_t kurama_theta_max = event.helix_t[itTpc][kurama_nh-1] + thetafac*kurama_theta_range;
	  
    // 	Double_t other_par[5];
    // 	other_par[0] = event.helix_cx[jgf];
    // 	other_par[1] = event.helix_cy[jgf];
    // 	other_par[2] = event.helix_z0[jgf];
    // 	other_par[3] = event.helix_r[jgf];
    // 	other_par[4] = event.helix_dz[jgf];
    // 	Int_t other_nh = event.helix_t[jgf].size();
    // 	Double_t other_theta_range = event.helix_t[jgf][0] - event.helix_t[jgf][other_nh-1];
    // 	Double_t other_theta_min = event.helix_t[jgf][other_nh-1] - thetafac*other_theta_range;
    // 	Double_t other_theta_max = event.helix_t[jgf][0] + thetafac*other_theta_range;
    // 	double thetaKurama,thetaOther;
    // 	vertex = Kinematics::VertexPointHelix(kurama_par,other_par,
    // 					      kurama_theta_min,kurama_theta_max,
    // 					      other_theta_min,other_theta_max,
    // 					      thetaKurama,thetaOther,
    // 					      dist);			    
    // 	Bool_t vtxouttgt
    // 	  = dist < ppi_distcut 
    // 	  && ( TMath::Abs(vertex.x()) > 25. 
    // 	       || TMath::Abs(vertex.y()) > 25. 
    // 	       || TMath::Abs(vertex.z() - tpc::ZTarget) > 50.) 
    // 	  && ( vertex.z() - tpc::ZTarget>0 );
    // 	event.GFKuramaVtxOutTgt = vtxouttgt;
    // 	if(vtxouttgt) break;
    //   }
    // }
    if ( TPCAna.IsInsideTarget(track) ) {      
      event.insideTgt[itTpc] = true;
      // if(event.isKurama[itTpc]){
      // 	event.GFKuramaFromTgt = 1;
      // }            
      //TVector3 reaction_vertex(event.vtxTPC[itTpc],event.vtyTPC[itTpc],event.vtzTPC[itTpc]+tpc::ZTarget);
      TVector3 post; TVector3 momt; double lentgt; double toftgt;      
      if ( TPCAna.ExtrapolateToTargetCenter(track, post, momt, lentgt) ){
          x0[ntrack_intarget] = post.x();
          y0[ntrack_intarget] = post.y();
          u0[ntrack_intarget] = momt.x()/momt.z();
          v0[ntrack_intarget] = momt.y()/momt.z();
          ntrack_intarget++;
      }
    } else event.insideTgt[itTpc] = 0;

    {
      TVector3 vertex = Kinematics::MultitrackVertex(ntrack_intarget, x0, y0, u0, v0);
      event.ntTpc_inside = ntrack_intarget;
      event.prodvtx_x = vertex.x();
      event.prodvtx_y = vertex.y();
      event.prodvtx_z = vertex.z();
      //TVector3 vertex(event.vtxTPC[igf],event.vtyTPC[igf],event.vtzTPC[igf]+tpc::ZTarget);
      double mom = track->GetMom0().Mag();
      Int_t htof_hitid; Double_t len_htof; TVector3 pos_htof;
      Bool_t htofextrapvtx =
	TPCAna.TPCHTOFTrackMatching(itTpc, vertex,
				    event.HtofSeg, event.posHtof,
				    htof_hitid, len_htof, pos_htof);      
      //      if ( TPCAna.TPCHTOFTrackMatching(itTpc,tgt_center,event.HtofSeg,event.posHtof,
      //				       htof_seg, tracklen_htof, pos_htof) ) {
      if(htofextrapvtx){
	std::cout << " htof_hitid: " <<  htof_hitid << " len_htof: "  << len_htof << " pos_htof: " << pos_htof << std::endl;
	std::cout << " nhHtof: " << event.nhHtof << std::endl;
      }
      for ( int ih=0; ih<event.nhHtof; ++ih ) {	
	if ( ih == htof_hitid && htofextrapvtx ) {
	  int hitid_htof = ih;
	  isInsideTarget = true;
	  //event.m2Htof[itTpc] = m2Htof;
	  event.hitidHtof[itTpc] = ih;
	  event.tracklenHtof[itTpc] = len_htof;
	  //event.GFtrack2vtxdist[igf] = track2tgt_dist;
	  event.posHtof_x[itTpc] = pos_htof.x();
	  event.posHtof_y[itTpc] = pos_htof.y();
	  event.posHtof_z[itTpc] = pos_htof.z();
	  event.segHtof[itTpc] = event.HtofSeg[ih];
	  event.tofHtof[itTpc] = event.tHtof[ih];
	  Double_t beta = len_htof/event.tHtof[ih]/MathTools::C();
	  Double_t tof = event.tofHtof[itTpc];
	  event.invbetaHtof[itTpc] = 1./beta;
	  event.m2Htof[itTpc] = Kinematics::MassSquare(mom, len_htof, tof);
	  //Double_t mass2 = Kinematics::MassSquare(event.GFmom[itTpc][0], len_htof, event.tHtof[hitid_htof]);
	  //          event.m2Htof[itTpc] = mass2;
	  event.nsigma_tritonHtof[itTpc] = Kinematics::HypTPCHTOFNsigmaTriton(mom, len_htof, event.tHtof[hitid_htof]);
	  event.nsigma_deutronHtof[itTpc] = Kinematics::HypTPCHTOFNsigmaDeutron(mom, len_htof, event.tHtof[hitid_htof]);
	  event.nsigma_protonHtof[itTpc] = Kinematics::HypTPCHTOFNsigmaProton(mom, len_htof, event.tHtof[hitid_htof]);
	  event.nsigma_kaonHtof[itTpc] = Kinematics::HypTPCHTOFNsigmaKaon(mom, len_htof, event.tHtof[hitid_htof]);
	  event.nsigma_pionHtof[itTpc] = Kinematics::HypTPCHTOFNsigmaPion(mom, len_htof, event.tHtof[hitid_htof]);
	  event.nsigma_electronHtof[itTpc] = Kinematics::HypTPCHTOFNsigmaElectron(mom, len_htof, event.tHtof[hitid_htof]);
	  break;
	}
      }
    }
    //    }
    //event.tgt2htofflag[itTpc] = isInsideTarget;
    if(!isInsideTarget) continue;
    int seghtof = event.HtofSeg[itTpc];
    Int_t charge = event.charge[itTpc];
    double dEdxtpc = event.dEdx[itTpc];
    double be = 0.;
    double cut_tofpid_min=3.0;   //should be better written in USER param
    double cut_tofpid_max=3.0;   //should be better written in USER param  
    double cut_dedxpid_min=3.0;   //should be better written in USER param
    double cut_dedxpid_max=3.0;   //should be better written in USER param
    double cut_dedxpid_veto=2.0;   //should be better written in USER param    
    Bool_t flagp[kNpid+1] = {};
    if( nsigma_pi > -cut_dedxpid_min && nsigma_pi < cut_dedxpid_max
	&& !event.isElectron[itTpc] ) {
      flagp[kPidPi] = true;
    }
    if( nsigma_k > -cut_dedxpid_min && nsigma_k < cut_dedxpid_max
	//&& nsigma_pi > cut_dedxpid_veto
	&& !event.isElectron[itTpc] ){	
      flagp[kPidK] = true;  
    }
    if( nsigma_p > -cut_dedxpid_min && nsigma_p < cut_dedxpid_max
	&& nsigma_pi > cut_dedxpid_veto	
	&& !event.isElectron[itTpc] ){     
      flagp[kPidP] = true;            
    }
    if( nsigma_d > -cut_dedxpid_min && nsigma_d < cut_dedxpid_max
	&& !event.isElectron[itTpc] ){
      flagp[kPidD] = true;            
    }    
    if( event.isElectron[itTpc] ){
      flagp[kPidE] = true;                  
    }
    if( !event.isElectron[itTpc] ){
      flagp[kPidAll] = true;
    }
   
    Int_t type=typeGenHid;
    int chargeid = (charge>0) ? 0 : (charge<0) ? 1 : -2;
    if(chargeid==-2) continue;    
    int beid = 0;	    
    double mom = event.mom0[itTpc];
    int momid = pidlikeli::MomToBin(mom);    
    if (momid < 0) continue;
    if (momid >= nbinpoq) momid = nbinpoq - 1;
    
    HF1(11300+chargeid*10, m2Htof*charge);
    
    for(int ip=0; ip<kNpid+1; ip++){
      if(flagp[ip]){
	HF1(11301+chargeid*10+ip+1, m2Htof*charge);
	HF2((type+1)*fac_t+ip*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2Htof*charge, dEdxtpc);
      }
    }
    Int_t pid=-1;
    if(KPEvent){
      for(int ibe=0; ibe<kNbe; ibe++){ // for selecting quasi-free Kaon and K star short
	std::cout << __FILE__ << " " << __LINE__ << " ibe:" << ibe << " BE:" << BE << " m2*charge:" << m2Htof*charge <<  " dEdx:" << dEdxtpc << std::endl;
	if(ibe==0) continue;
	if(flagp[kPidK]&&BE>pidlikeli::cutbemin[ibe]&&BE<pidlikeli::cutbemax[ibe]){
	  type=typeKmHid;
	  int pid=kPidK;
	  std::cout << __FILE__ << " " << __LINE__ << " BE:" << ibe << " m2*charge:" << m2Htof*charge <<  " dEdx:" << dEdxtpc << std::endl;
	  HF2((type+1)*fac_t+pid*fac_p+chargeid*fac_c+ibe*fac_b+momid*fac_m, m2Htof*charge, dEdxtpc);
	}
      }
    }
  }
  std::cout << __FILE__ << " " << __LINE__ << std::endl;  
  //Lambda reconstruction
  std::vector<Int_t> L_p_id_container, L_pi_id_container;
  std::vector<TVector3> L_p_mom_container, L_pi_mom_container;
  std::vector<TVector3> L_mom_container, L_vert_container;  
  std::vector<Double_t> L_mass_container;
  std::vector<Double_t> L_ppidist_container;
  std::vector<Double_t> L_ppiangle_container;  
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_targetvtx_container;
  std::vector<Bool_t> L_p_htofextrap_container,L_pi_htofextrap_container;  
  std::vector<Double_t> L_p_mass2_container,L_pi_mass2_container;
  std::vector<Double_t> L_p_posY_container,L_pi_posY_container;  
  std::vector<TVector3> L_p_htofpos_container,L_pi_htofpos_container;
  std::vector<Int_t> L_p_htofhitid_container,L_pi_htofhitid_container;
  std::vector<Int_t> L_p_htofseg_container,L_pi_htofseg_container;    
  std::vector<Double_t> L_p_diffY_container,L_pi_diffY_container;
  std::vector<Double_t> L_p_tracklen_container,L_pi_tracklen_container;
  //L candidates searching
  Int_t l_candidates = 0;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){ // proton
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if(event.charge[it1]!=1) continue;
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
	if(event.isElectron[it2]==1) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if(event.charge[it2]!=-1) continue;
	Double_t pi_par[5];
	pi_par[0] = event.helix_cx[it2];
	pi_par[1] = event.helix_cy[it2];
	pi_par[2] = event.helix_z0[it2];
	pi_par[3] = event.helix_r[it2];
	pi_par[4] = event.helix_dz[it2];
	Int_t pi_nh = event.helix_t[it2].size();
	Double_t pi_theta_min = TMath::Max(event.helix_t[it2][pi_nh-1] - vtx_scan_rangeInsideL/pi_par[3], event.helix_t[it2][0]);
	Double_t pi_theta_max = event.helix_t[it2][0] + vtx_scan_range/pi_par[3];
	TVector3 pi_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 pi_end = TVector3(event.calpos_x[it2][pi_nh-1], event.calpos_y[it2][pi_nh-1], event.calpos_z[it2][pi_nh-1]);	
	Double_t ppi_dist = 10000.;
	TVector3 p_mom; TVector3 pi_mom; TVector3 l_mom;
	TVector3 l_vertex = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, p_mom, pi_mom, l_mom, ppi_dist);
	if(TMath::IsNaN(ppi_dist)) continue;
	Int_t hitidHtof1,hitidHtof2;
	Double_t tracklenHtof1,tracklenHtof2;
	TVector3 posHtof1,posHtof2;
	bool match1 = TPCAna.TPCHTOFTrackMatching(it1,l_vertex,event.HtofSeg,event.posHtof,hitidHtof1,tracklenHtof1,posHtof1);
	if(match1&&tracklenHtof1<10){
	  std::cout << __FILE__ << " " << __LINE__ << " "
		    << " hitidHtof1:" <<  hitidHtof1 << " tracklenHtof1:" << tracklenHtof1
		    << " extrapolated posY1:" <<  posHtof1.y()  << " posHtof1[Y]:" << event.posHtof[hitidHtof1]
		    << std::endl;
	}
	bool match2 = TPCAna.TPCHTOFTrackMatching(it2,l_vertex,event.HtofSeg,event.posHtof,hitidHtof2,tracklenHtof2,posHtof2);
	//if ( !match2 ) continue;
	// std::cout << __FILE__ << " " << __LINE__ << " "
	// 	  << " hitidHtof1:" <<  hitidHtof2 << " tracklenHtof1:" << tracklenHtof2
	// 	  << " extrapolated posY2:" <<  posHtof2.y()  << " posHtof2[Y]:" << event.posHtof[hitidHtof2]
	// 	  << std::endl;	
	// if((event.pid[it1]&4)!=4) continue;
	// if((event.pid[it2]&1)!=1) continue; //select pi-like
	bool pidcondP = event.nsigma_proton[it1]<3&&event.nsigma_proton[it1]>-3&&event.nsigma_pion[it1]>2;
	bool pidcondPi = event.nsigma_pion[it2]<3&&event.nsigma_pion[it2]>-3;
	if(!pidcondP) continue;
	if(!pidcondPi) continue;
	//std::cout << __FILE__ << " " << __LINE__ << " PID was done using dEdx info" << std::endl;
	l_mom = pi_mom + p_mom;
	Double_t l_target_dist;
	TVector3 l_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, l_vertex, l_mom, l_target_dist);
	TVector3 l_flight = l_vertex - l_pos_tgt;
	Double_t l_tof = Kinematics::CalcTimeOfFlight(l_mom.Mag(), l_flight.Mag(), pdg::LambdaMass());
	double mass2p = -999.; double mass2pi = -999.;
	if(match1&&match2){
	  mass2p  = Kinematics::MassSquare(p_mom.Mag(), tracklenHtof1, event.tHtof[hitidHtof1] - l_tof);
	  mass2pi = Kinematics::MassSquare(pi_mom.Mag(), tracklenHtof2, event.tHtof[hitidHtof2] - l_tof);
	}
	TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
	TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));	
	TLorentzVector Llambda = Lp + Lpi;
	if(TMath::Abs(l_vertex.x()) > 250. ||
	   TMath::Abs(l_vertex.z()) > 250. ||
	   TMath::Abs(l_vertex.y()) > 250.) continue; //Vertex cut	
	Double_t pi_vertex_dist; Double_t p_vertex_dist;
	if(!Kinematics::HelixDirection(l_vertex, p_start, p_end, p_vertex_dist) ||
	   !Kinematics::HelixDirection(l_vertex, pi_start, pi_end, pi_vertex_dist)) continue;
	// std::cout << __FILE__ << " " << __LINE__ << " "
	// 	  << " pi_vertex_dist: " << pi_vertex_dist
	// 	  << " p_vertex_dist: " << p_vertex_dist
	// 	  << " ppi_dist: " << ppi_dist	  
	// 	  << std::endl;;
	if(pi_vertex_dist > pi_vtx_distcut) continue;
	if(p_vertex_dist > p_vtx_distcut) continue;
	if(ppi_dist > ppi_distcut) continue;
	if(TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut)continue;
	event.lflag = true;
	Double_t ltarget_dist;
	TVector3 ltarget_vtx = 
	  Kinematics::CalcCloseDistLambda(tgtpos,
					  l_vertex,
					  l_mom,
					  ltarget_dist);
	L_p_id_container.push_back(it1);
	L_pi_id_container.push_back(it2);
	L_mass_container.push_back(Llambda.M());
	L_mom_container.push_back(l_mom);
	L_p_mom_container.push_back(p_mom);
	L_pi_mom_container.push_back(pi_mom);
	L_ppidist_container.push_back(ppi_dist);
	L_vert_container.push_back(l_vertex);
	L_targetdist_container.push_back(ltarget_dist);
	L_targetvtx_container.push_back(ltarget_vtx);
	L_p_htofextrap_container.push_back(match1);	
	L_pi_htofextrap_container.push_back(match2);		
	if(match1&&match2){
	  double diffY1 = posHtof1.y() - event.posHtof[hitidHtof1];
	  double diffY2 = posHtof2.y() - event.posHtof[hitidHtof2];
	  double posY1 = posHtof1.y();
	  double posY2 = posHtof2.y();
	  Int_t hitsegHtof1 = event.HtofSeg[hitidHtof1];
	  Int_t hitsegHtof2 = event.HtofSeg[hitidHtof2];
	  L_p_mass2_container.push_back(mass2p);
	  L_pi_mass2_container.push_back(mass2pi);
	  L_p_tracklen_container.push_back(tracklenHtof1);
	  if(tracklenHtof1<10){
	    std::cout << __FILE__ << " " << __LINE__ << " "
		      << " hitidHtof1:" <<  hitidHtof1
		      << " tracklenHtof1:" << tracklenHtof1
		      << " extrapolated posY1:" <<  posHtof1.y()
		      << " posHtof1[Y]:" << event.posHtof[hitidHtof1]
		      << std::endl;
	  }
	  L_pi_tracklen_container.push_back(tracklenHtof2);
	  L_p_htofpos_container.push_back(posHtof1);	
	  L_pi_htofpos_container.push_back(posHtof2);
	  L_p_htofhitid_container.push_back(hitidHtof1);
	  L_pi_htofhitid_container.push_back(hitidHtof2);
	  L_p_htofseg_container.push_back(hitsegHtof1);
	  L_pi_htofseg_container.push_back(hitsegHtof2);		
	  L_p_posY_container.push_back(posY1);
	  L_pi_posY_container.push_back(posY2);
	  L_p_diffY_container.push_back(diffY1);
	  L_pi_diffY_container.push_back(diffY2);	  
	}
	l_candidates++;	
	//if(TMath::Abs(Llambda.M()-LambdaMass)<lambda_masscut_final)
	{
	  //hist	
	  TVector3 lvtxdist = l_vertex - tgtpos;
	  HF1(20010, ppi_dist);
	  HF1(20011, p_vertex_dist);
	  HF1(20012, pi_vertex_dist);			
	  HF1(20013, lvtxdist.Mag());
	}
      } //it2
    } //it1
  }
  std::cout << __FILE__ << " " << __LINE__ << std::endl;  
  Int_t best_l = -1; Double_t prev_massdiff_l = 9999.;
  for(Int_t candi=0;candi<l_candidates;candi++){
    Double_t diff = TMath::Abs(L_mass_container[candi] - LambdaMass);
    if(prev_massdiff_l > diff){
      prev_massdiff_l = diff;
      //if(diff<lambda_masscut_final){
      best_l = candi;
      std::cout << "best lambda: " << best_l << std::endl;
      //}
    }
  }
  if(best_l!=-1){
    event.lmass = L_mass_container[best_l];
    event.ldecayvtx_x = L_vert_container[best_l].x();
    event.ldecayvtx_y = L_vert_container[best_l].y();
    event.ldecayvtx_z = L_vert_container[best_l].z();
    event.lmom   = L_mom_container[best_l].Mag();    
    event.lmom_x = L_mom_container[best_l].x();
    event.lmom_y = L_mom_container[best_l].y();
    event.lmom_z = L_mom_container[best_l].z();
    event.ppi_dist = L_ppidist_container[best_l];
    //event.ppiangle = L_ppiangle_container[best_l];
    double momp  = L_p_mom_container[best_l].Mag();
    double mompi = L_pi_mom_container[best_l].Mag();
    event.ldecays_id.push_back(L_p_id_container[best_l]);
    event.ldecays_mom.push_back(momp);
    event.ldecays_mom_x.push_back(L_p_mom_container[best_l].x());
    event.ldecays_mom_y.push_back(L_p_mom_container[best_l].y());
    event.ldecays_mom_z.push_back(L_p_mom_container[best_l].z());
    //  event.ldecays_theta.push_back(l_p_mom_container[best].Theta()*TMath::RadToDeg());
    event.ldecays_id.push_back(L_pi_id_container[best_l]);
    event.ldecays_mom.push_back(mompi);
    event.ldecays_mom_x.push_back(L_pi_mom_container[best_l].x());
    event.ldecays_mom_y.push_back(L_pi_mom_container[best_l].y());
    event.ldecays_mom_z.push_back(L_pi_mom_container[best_l].z());
    //htof ana
    bool match1 = L_p_htofextrap_container[best_l];
    bool match2 = L_pi_htofextrap_container[best_l]; 
    event.ldecays_htofextrap.push_back(match1);
    event.ldecays_htofextrap.push_back(match2);
    if(match1&&match2){
      double mass2p  = L_p_mass2_container[best_l];
      double mass2pi = L_pi_mass2_container[best_l];
      event.ldecays_mass2.push_back(mass2p);
      event.ldecays_mass2.push_back(mass2pi);    
      double len1 = L_p_tracklen_container[best_l];
      double len2 = L_pi_tracklen_container[best_l];
      TVector3 poshtofp  = L_p_htofpos_container[best_l];
      TVector3 poshtofpi = L_pi_htofpos_container[best_l];    
      double posY1 = L_p_posY_container[best_l];
      double posY2 = L_pi_posY_container[best_l];
      double diffY1 = L_p_diffY_container[best_l];
      double diffY2 = L_pi_diffY_container[best_l];
      int hitid1 = L_p_htofhitid_container[best_l];
      int hitid2 = L_pi_htofhitid_container[best_l];
      int hitseg1 = L_p_htofseg_container[best_l];
      int hitseg2 = L_pi_htofseg_container[best_l];                
      event.ldecays_tracklen.push_back(len1);
      event.ldecays_tracklen.push_back(len2);
      event.ldecays_htofpos_x.push_back(poshtofp.x());
      event.ldecays_htofpos_y.push_back(poshtofp.y());
      event.ldecays_htofpos_z.push_back(poshtofp.z());
      event.ldecays_htofpos_x.push_back(poshtofpi.x());
      event.ldecays_htofpos_y.push_back(poshtofpi.y());
      event.ldecays_htofpos_z.push_back(poshtofpi.z());
      event.ldecays_htofhitid.push_back(hitid1);
      event.ldecays_htofhitid.push_back(hitid2);
      event.ldecays_htofseg.push_back(hitseg1);
      event.ldecays_htofseg.push_back(hitseg2);                            
      if(TMath::Abs(event.lmass-LambdaMass)<lambda_masscut_final){
	HF1(20301, mass2p);
	HF1(20302, mass2pi);    
	HF1(10101, len1);
	HF1(10102, posY1); 
	HF1(10103, diffY1);
	HF1(10201, len2);
	HF1(10202, posY2);
	HF1(10203, diffY2);
      }    
    }
    {
      bool match1 = L_p_htofextrap_container[best_l];
      bool match2 = L_pi_htofextrap_container[best_l];       
      int p_id = L_p_id_container[best_l];
      double dedx = event.dEdx[p_id];
      int charge = event.charge[p_id];
      int type = typeLHid; 
      int PID = kPidP;
      int chargeid = plusHid;
      int momid = pidlikeli::MomToBin(momp);
      if(match1&&match2){
	double mass2p  = L_p_mass2_container[best_l];	
	HF2((type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2p*charge, dedx);
      }
    }
    {
      int pi_id = L_pi_id_container[best_l];
      double dedx = event.dEdx[pi_id];
      int charge = event.charge[pi_id];
      int type = typeLHid;
      int PID = kPidPi;
      int chargeid = minusHid;
      int momid = pidlikeli::MomToBin(mompi);
      if(match1&&match2){
	double mass2pi = L_pi_mass2_container[best_l];	
	// std::cout << "Fill hist " << (type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m
	// 	  << " with " << mass2pi*charge
	// 	  << ", " << dedx
	// 	  << std::endl;
	HF2((type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2pi*charge,dedx);
      }
    }
  }
  std::cout << __FILE__ << " " << __LINE__ << std::endl;  
  //K0 reconstruction   
  std::vector<Int_t> K0_pip_id_container, K0_pim_id_container;
  std::vector<TVector3> K0_pip_mom_container, K0_pim_mom_container;
  std::vector<Double_t> K0_mass_container;
  std::vector<Double_t> K0_pip_mass2_container, K0_pim_mass2_container;
  std::vector<Double_t> K0_pipidist_container;
  std::vector<Double_t> K0_pipiangle_container;  
  std::vector<Double_t> K0_targetdist_container;
  std::vector<Bool_t>   K0_pip_htofextrap_container, K0_pim_htofextrap_container;  
  std::vector<Double_t> K0_pip_tracklen_container, K0_pim_tracklen_container;
  std::vector<TVector3> K0_pip_htofpos_container, K0_pim_htofpos_container;  
  std::vector<TVector3> K0_mom_container, K0_vert_container, K0_targetvtx_container;
  //L candidates searching  
  Int_t k0_candidates = 0;
  {
    for(Int_t it1=0;it1<ntTpc;it1++){ // pi+
      if(!KPEvent) continue;
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if(event.charge[it1]!=1) continue;
      Double_t pip_par[5];
      pip_par[0] = event.helix_cx[it1];
      pip_par[1] = event.helix_cy[it1];
      pip_par[2] = event.helix_z0[it1];
      pip_par[3] = event.helix_r[it1];
      pip_par[4] = event.helix_dz[it1];
      Int_t pip_nh = event.helix_t[it1].size();
      Double_t pip_theta_min = event.helix_t[it1][0] - vtx_scan_range/pip_par[3];
      Double_t pip_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_rangeInsideL/pip_par[3], event.helix_t[it1][pip_nh-1]);
      TVector3 pip_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 pip_end = TVector3(event.calpos_x[it1][pip_nh-1], event.calpos_y[it1][pip_nh-1], event.calpos_z[it1][pip_nh-1]);
      for(Int_t it2=0;it2<ntTpc;it2++){ // pi-
	if(it1==it2) continue;
	if(event.isElectron[it2]==1) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if(event.charge[it2]!=-1) continue;	
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
	TVector3 pip_mom; TVector3 pim_mom; TVector3 k0_mom;
	TVector3 k0_vertex = Kinematics::LambdaVertex(dMagneticField, pip_par, pim_par, pip_theta_min, pip_theta_max, pim_theta_min, pim_theta_max, pip_mom, pim_mom, k0_mom, pipi_dist);
	if(TMath::IsNaN(pipi_dist)) continue;
	Int_t hitidHtof1,hitidHtof2;
	Double_t tracklenHtof1,tracklenHtof2;
	TVector3 posHtof1,posHtof2;
	bool match1 = TPCAna.TPCHTOFTrackMatching(it1,k0_vertex,event.HtofSeg,event.posHtof,hitidHtof1,tracklenHtof1,posHtof1);
	//if ( !match1 ) continue;
	// std::cout << __FILE__ << " " << __LINE__ << " "
	// 	  << " hitidHtof1:" <<  hitidHtof1 << " tracklenHtof1:" << tracklenHtof1
	// 	  << " extrapolated posY1:" <<  posHtof1.y()  << " posHtof1[Y]:" << event.posHtof[hitidHtof1]
	// 	  << std::endl;
	// std::cout << __FILE__ << " " << __LINE__ << " HTOFTrackMatching2:"<< std::endl;
	bool match2 = TPCAna.TPCHTOFTrackMatching(it2,k0_vertex,event.HtofSeg,event.posHtof,hitidHtof2,tracklenHtof2,posHtof2);
	//if ( !match2 ) continue;
	// std::cout << __FILE__ << " " << __LINE__ << " "
	// 	  << " hitidHtof1:" <<  hitidHtof2 << " tracklenHtof1:" << tracklenHtof2
	// 	  << " extrapolated posY2:" <<  posHtof2.y()  << " posHtof2[Y]:" << event.posHtof[hitidHtof2]
	// 	  << std::endl;	
	bool pidcondPip = event.nsigma_pion[it1]<3&&event.nsigma_pion[it1]>-3;
	bool pidcondPim = event.nsigma_pion[it2]<3&&event.nsigma_pion[it2]>-3;
	if(!pidcondPip) continue;
	if(!pidcondPim) continue;	
	// if((event.pid[it1]&1)!=1) continue;
	// if((event.pid[it2]&1)!=1) continue;
	k0_mom = pip_mom + pim_mom;
	Double_t k0_target_dist;
	TVector3 k0_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, k0_vertex, k0_mom, k0_target_dist);
	TVector3 k0_flight = k0_vertex - k0_pos_tgt;
	Double_t k0_tof = Kinematics::CalcTimeOfFlight(k0_mom.Mag(), k0_flight.Mag(), pdg::K0Mass());
	double mass2pip = -999.0; double mass2pim = -999.0;
	if(match1&&match2){
	  mass2pip = Kinematics::MassSquare(pip_mom.Mag(), tracklenHtof1, event.tHtof[hitidHtof1] - k0_tof);
	  mass2pim = Kinematics::MassSquare(pim_mom.Mag(), tracklenHtof2, event.tHtof[hitidHtof2] - k0_tof);
	}
	TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
	TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));	
	TLorentzVector Lkaon0 = Lpip + Lpim;
	if(TMath::Abs(k0_vertex.x()) > 250. ||
	   TMath::Abs(k0_vertex.z()) > 250. ||
	   TMath::Abs(k0_vertex.y()) > 250.) continue; //Vertex cut

	Double_t pip_vertex_dist; Double_t pim_vertex_dist;
	if(!Kinematics::HelixDirection(k0_vertex, pip_start, pip_end, pip_vertex_dist) ||
	   !Kinematics::HelixDirection(k0_vertex, pim_start, pim_end, pim_vertex_dist)) continue;

	if(pip_vertex_dist > pip_vtx_distcut) continue;
	if(pim_vertex_dist > pim_vtx_distcut) continue;
	if(pipi_dist > pipi_distcut || TMath::Abs(Lkaon0.M() - K0Mass) > k0_masscut) continue;

	Double_t k0target_dist;
	TVector3 k0target_vtx =
	  Kinematics::CalcCloseDistLambda(tgtpos,
					  k0_vertex,
					  k0_mom,
					  k0target_dist);
	
	K0_pip_id_container.push_back(it1);
	K0_pim_id_container.push_back(it2);
	K0_mass_container.push_back(Lkaon0.M());
	K0_mom_container.push_back(k0_mom);
	K0_pip_mom_container.push_back(pip_mom);
	K0_pim_mom_container.push_back(pim_mom);
	K0_pipidist_container.push_back(pipi_dist);
	K0_pip_htofextrap_container.push_back(match1);
	K0_pim_htofextrap_container.push_back(match2);
	if(match1&&match2){
	  K0_pip_mass2_container.push_back(mass2pip);
	  K0_pim_mass2_container.push_back(mass2pim);	
	  K0_pip_htofpos_container.push_back(posHtof1);
	  K0_pim_htofpos_container.push_back(posHtof2);
	  K0_pip_tracklen_container.push_back(tracklenHtof1);
	  K0_pim_tracklen_container.push_back(tracklenHtof2);
	}
	//	K0_pipiangle_container.push_back(pipi_anlge);	
	K0_vert_container.push_back(k0_vertex);
	K0_targetdist_container.push_back(k0target_dist);
	K0_targetvtx_container.push_back(k0target_vtx);
	k0_candidates++;
	HF1( 30010, pipi_dist );
	HF1( 30011, pip_vertex_dist );
	HF1( 30012, pim_vertex_dist );
	HF1( 30013, k0_vertex.Mag() );			
      } //it2
    } //it1
  }
  Int_t best_k0 = -1; Double_t prev_massdiff_k0 = 9999.;    
  for(Int_t candi=0;candi<k0_candidates;candi++){
    Double_t diff = TMath::Abs(K0_mass_container[candi] - K0Mass);
    if(prev_massdiff_k0 > diff){
      prev_massdiff_k0 = diff;
      best_k0 = candi;
      std::cout << "best Kaon0: " << best_k0 << std::endl;
    }
  }
  
  if(best_k0!=-1){
    event.k0mass = K0_mass_container[best_k0];
    event.k0decayvtx_x = K0_vert_container[best_k0].x();
    event.k0decayvtx_y = K0_vert_container[best_k0].y();
    event.k0decayvtx_z = K0_vert_container[best_k0].z();
    event.k0mom_x = K0_mom_container[best_k0].x();
    event.k0mom_y = K0_mom_container[best_k0].y();
    event.k0mom_z = K0_mom_container[best_k0].z();
    event.pipi_dist = K0_pipidist_container[best_k0];
    // event.pipiangle = K0_pipiangle_container[best_k0];
    double k0mass = K0_mass_container[best_k0];       
    event.k0decays_id.push_back(K0_pip_id_container[best_k0]);
    event.k0decays_mom.push_back(K0_pip_mom_container[best_k0].Mag());
    event.k0decays_mom_x.push_back(K0_pip_mom_container[best_k0].x());
    event.k0decays_mom_y.push_back(K0_pip_mom_container[best_k0].y());
    event.k0decays_mom_z.push_back(K0_pip_mom_container[best_k0].z());
    //  event.ldecays_theta.push_back(l_p_mom_container[best].Theta()*TMath::RadToDeg());  
    event.k0decays_id.push_back(K0_pim_id_container[best_k0]);
    event.k0decays_mom.push_back(K0_pim_mom_container[best_k0].Mag());
    event.k0decays_mom_x.push_back(K0_pim_mom_container[best_k0].x());
    event.k0decays_mom_y.push_back(K0_pim_mom_container[best_k0].y());
    event.k0decays_mom_z.push_back(K0_pim_mom_container[best_k0].z());
    bool match1 = K0_pip_htofextrap_container[best_k0];
    bool match2 = K0_pim_htofextrap_container[best_k0];
    event.k0decays_htofextrap.push_back(match1);
    event.k0decays_htofextrap.push_back(match2);        
    if(match1&&match2){
      double mass2pip = K0_pip_mass2_container[best_k0];
      double mass2pim = K0_pim_mass2_container[best_k0];      
      event.k0decays_mass2.push_back(K0_pip_mass2_container[best_k0]);
      event.k0decays_mass2.push_back(K0_pim_mass2_container[best_k0]);
      event.k0decays_tracklen.push_back(K0_pip_tracklen_container[best_k0]);
      event.k0decays_tracklen.push_back(K0_pim_tracklen_container[best_k0]);
      event.k0decays_htofpos_x.push_back(K0_pip_htofpos_container[best_k0].x());
      event.k0decays_htofpos_y.push_back(K0_pip_htofpos_container[best_k0].y());
      event.k0decays_htofpos_z.push_back(K0_pip_htofpos_container[best_k0].z());
      HF1( 30301, mass2pip ); 
      HF1( 30302, mass2pim );          
    }
    HF1( 30201, k0mass );
    {      
      if(match1&&match2){
	double mass2 = K0_pip_mass2_container[best_k0];
	int pip_id = K0_pip_id_container[best_k0];
	double dedx = event.dEdx[pip_id];
	double mom = K0_pip_mom_container[best_k0].Mag();
	int charge = event.charge[pip_id];	
	int type = typeK0Hid;
	int PID = kPidPi;
	int chargeid = plusHid;
	int momid = pidlikeli::MomToBin(mom);
	// std::cout << "Fill hist " << (type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m
	// 		<< " with " << mass2pip*charge
	// 		<< ", " << dedx
	// 		<< std::endl;	
	HF2( (type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2*charge,dedx);
      }
    }
    {
      if(match1&&match2){
	double mass2 = K0_pim_mass2_container[best_k0];
	int pim_id = K0_pim_id_container[best_k0];
	double dedx = event.dEdx[pim_id];
	double mom = K0_pim_mom_container[best_k0].Mag();
	int charge = event.charge[pim_id];	
	int type = typeK0Hid;
	int PID = kPidPi;
	int chargeid = minusHid;
	int momid = pidlikeli::MomToBin(mom);
	HF2( (type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2*charge,dedx);
      }
    }                
  }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
  
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

  static const auto KKEvent = gUser.GetParameter("KKEvent");
  static const auto KPEvent = gUser.GetParameter("KPEvent");
  static const auto KHeavyEvent = gUser.GetParameter("KHeavyEvent");

  const Double_t nbinIML = 80;  
  const Double_t minIML = 1.04;
  const Double_t maxIML = 1.20;  
  const Double_t enestep = 0.002;
  const Double_t minIMK0 = 0.3;
  const Double_t maxIMK0 = 0.7;
  const Double_t nbinIMK0 = (maxIMK0-minIMK0)/enestep;
  const Double_t nbinM2P = 630;  
  const Double_t minM2P = 0.;
  const Double_t maxM2P = 1.26;  
  const Double_t nbinM2Pi = 630;  
  const Double_t minM2Pi = 0.;
  const Double_t maxM2Pi = 1.26;      
  
#if SaveHistograms
  HB1(1, "Status", 21, 0., 21. );
  HB1(12, "K18 TPC tagging", 2, 0., 2. );
  HB1(13, "Kurama TPC tagging", 2, 0., 2. );
  HB1(14, "KK TPC tagging", 2, 0., 2. );
  HB1(22, "K18 TPC tagging", 2, 0., 2. );
  HB1(23, "Kurama TPC tagging", 2, 0., 2. );
  HB1(24, "KK TPC tagging", 2, 0., 2. );
  HB1(32, "K18 TPC tagging", 2, 0., 2. );
  HB1(33, "Kurama TPC tagging", 2, 0., 2. );
  HB1(34, "KK TPC tagging", 2, 0., 2. );
  HB1(42, "K18 TPC tagging", 2, 0., 2. );
  HB1(43, "Kurama TPC tagging", 2, 0., 2. );
  HB1(44, "KK TPC tagging", 2, 0., 2. );

  HB2(20, "1/#beta;p/q [GeV/#font[12]{c}];1/#beta", 1000, -2.0, 2.0, 1000, 0.0, 5.0);

  // missing mass with scat angle, vtx,
  HB1(100, "MissMass", 3600, -1.0, 17. );
  HB1(101, "MissMassCorr", 3600, -1.0, 17. );
  HB1(102, "MissMassCorrDE", 3600, -1.0, 17. );
  HB1(103, "MissMassTPC", 3600, -1.0, 17. );
  HB1(104, "MissMassCorrTPC", 3600, -1.0, 17. );
  HB1(105, "MissMassCorrDETPC", 3600, -1.0, 17. );
  HB1(106, "MissMassNuclTPC (3.5<thetaTPC<4.5)", 3600, -1.0, 17. );
  HB1(107, "MissMassNuclCorrTPC (3.5<thetaTPC<4.5)", 3600, -1.0, 17. );
  HB1(108, "MissMassNuclCorrDETPC (3.5<thetaTPC<4.5)", 3600, -1.0, 17. );

  HB1(1001, "P K18", 800, 1.4, 2.2);
  HB1(1002, "P Kurama", 600, 0, 3);
  HB1(1003, "Charge*Mass", 600, -1., 2.5);
  HB2(1004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(1005, "Sqaured Mass", 200, -0.2, 1.2);
  HB1(1006, "Sqaured Mass [P Corrected]", 200, -0.2, 1.2);
  HB2(1010, "Vertex", 400, -500, 500, 400, -150, 150);
  HB1(1011, "Closest distance", 50, 0., 50.);
  HB1(1012, "Vertex X", 200, -100, 100);
  HB1(1013, "Vertex Y", 200, -100, 100);
  HB1(1014, "Vertex Z", 200, -200, 200);
  HB1(2001, "P K18", 800, 1.4, 2.2);
  HB1(2002, "P Kurama", 600, 0, 3);
  HB1(2003, "Charge*Mass", 600, -1., 2.5);
  HB2(2004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB2(2010, "Vertex", 400, -500, 500, 400, -150, 150);
  HB1(2011, "Closest distance", 50, 0., 50.);
  HB1(2012, "Vertex X", 200, -100, 100);
  HB1(2013, "Vertex Y", 200, -100, 100);
  HB1(2014, "Vertex Z", 200, -200, 200);

  HB1(3001, "P K18", 800, 1.4, 2.2);
  HB1(3002, "P Kurama", 600, 0, 3);
  HB1(3003, "Charge*Mass", 600, -1., 2.);
  HB2(3004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(3010, "Closest distance", 50, 0., 50.);
  HB1(3011, "Vertex X", 200, -100, 100);
  HB1(3012, "Vertex Y", 200, -100, 100);
  HB1(3013, "Vertex Z", 200, -200, 200);
  HB1(3014, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB1(3015, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB2(3016, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(3017, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(3018, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(3019, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(3020, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(3021, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(3022, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(3023, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(3024, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(3025, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(3026, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(3027, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(3028, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(3029, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(3030, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(3031, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(3032, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(3033, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);

  HBProf(3116, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(3117, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(3118, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(3119, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(3120, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(3121, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(3122, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(3123, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(3124, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(3125, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(3126, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(3127, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(3128, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(3129, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(3130, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(3131, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(3132, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(3133, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);

  HB1(4001, "P K18", 800, 1.4, 2.2);
  HB1(4002, "P Kurama", 600, 0, 3);
  HB1(4003, "Charge*Mass", 600, -1., 2.);
  HB2(4004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(4010, "Closest distance", 50, 0., 50.);
  HB1(4011, "Vertex X", 200, -100, 100);
  HB1(4012, "Vertex Y", 200, -100, 100);
  HB1(4013, "Vertex Z", 200, -200, 200);
  HB1(4014, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB1(4015, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB2(4016, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4017, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4018, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4019, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(4020, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(4021, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(4022, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4023, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4024, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4025, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(4026, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(4027, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(4028, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(4029, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(4030, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(4031, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(4032, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(4033, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);

  HBProf(4116, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(4117, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(4118, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(4119, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(4120, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(4121, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(4122, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(4123, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(4124, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(4125, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(4126, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(4127, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(4128, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(4129, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(4130, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(4131, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(4132, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(4133, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);

  HB1(4201, "P K18", 800, 1.4, 2.2);
  HB1(4202, "P Kurama", 600, 0, 3);
  HB1(4203, "Charge*Mass", 700, -1., 2.5);
  HB2(4204, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(4210, "Closest distance", 50, 0., 50.);
  HB1(4211, "Vertex X", 200, -100, 100);
  HB1(4212, "Vertex Y", 200, -100, 100);
  HB1(4213, "Vertex Z", 200, -200, 200);
  HB1(4214, "MissingMass P", 200, 0.3, 1.3);
  HB1(4215, "MissingMass P", 200, 0.3, 1.3);
  HB2(4216, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4217, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4218, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4219, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(4220, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(4221, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(4222, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4223, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4224, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(4225, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(4226, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(4227, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);

  HBProf(4316, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(4317, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(4318, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(4319, "MissingMass%U Prof ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(4320, "MissingMass%U Prof ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(4321, "MissingMass%U Prof ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(4322, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(4323, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(4324, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(4325, "MissingMass%V Prof ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(4326, "MissingMass%V Prof ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(4327, "MissingMass%V Prof ", 160, -0.4, 0.4, 0.3, 1.3);

  HB1(5001, "P K18", 800, 1.4, 2.2);
  HB1(5002, "P Kurama", 600, 0, 3);
  HB1(5003, "Charge*Mass", 600, -1., 2.);
  HB2(5004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB2(5010, "Vertex", 400, -500, 500, 400, -150, 150);
  HB1(5011, "Closest distance", 200, 0., 200.);
  HB1(5012, "Vertex X", 400, -200, 200);
  HB1(5013, "Vertex Y", 400, -200, 200);
  HB1(5014, "Vertex Z", 500, -500, 500);

  HB1(6001, "P K18", 800, 1.4, 2.2);
  HB1(6002, "P Kurama", 600, 0, 3);
  HB1(6003, "Charge*Mass", 600, -1., 2.);
  HB2(6004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB2(6010, "Vertex", 400, -500, 500, 400, -150, 150);
  HB1(6011, "Closest distance", 50, 0., 50.);
  HB1(6012, "Vertex X", 200, -100, 100);
  HB1(6013, "Vertex Y", 200, -100, 100);
  HB1(6014, "Vertex Z", 200, -200, 200);

  HB2(6041, "K-, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(6042, "K-, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(6043, "K-, Utgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);
  HB2(6044, "K-, Vtgt w/ vs w/o TPC", 200, -0.1, 0.1, 200, -0.1, 0.1);
  HB2(6045, "K+, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(6046, "K+, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(6047, "K+, Utgt w/ vs w/o TPC", 200, -0.5, 0.5, 200, -0.5, 0.5);
  HB2(6048, "K+, Vtgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);

  HB1(7001, "P K18", 800, 1.4, 2.2);
  HB1(7002, "P Kurama", 600, 0, 3);
  HB1(7003, "Charge*Mass", 600, -1., 2.);
  HB2(7004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(7010, "Closest distance", 50, 0., 50.);
  HB1(7011, "Vertex X", 200, -100, 100);
  HB1(7012, "Vertex Y", 200, -100, 100);
  HB1(7013, "Vertex Z", 200, -200, 200);
  HB1(7014, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB1(7015, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB2(7016, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7017, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7018, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7019, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(7020, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(7021, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(7022, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7023, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7024, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(7025, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(7026, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(7027, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(7028, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(7029, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(7030, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(7031, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(7032, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(7033, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(7036, "#DeltaP%P_{calc.} ", 160, 0.3, 1.1, 200, -1., 1.);
  HB2(7039, "MissingMass%P_{calc.} ", 160, 0.3, 1.1, 140, 1.1, 1.8);

  HB2(7041, "K-, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(7042, "K-, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(7043, "K-, Utgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);
  HB2(7044, "K-, Vtgt w/ vs w/o TPC", 200, -0.1, 0.1, 200, -0.1, 0.1);
  HB2(7045, "K+, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(7046, "K+, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(7047, "K+, Utgt w/ vs w/o TPC", 200, -0.5, 0.5, 200, -0.5, 0.5);
  HB2(7048, "K+, Vtgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);

  HBProf(7116, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(7117, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(7118, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(7119, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(7120, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(7121, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(7122, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(7123, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(7124, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(7125, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(7126, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(7127, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(7128, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(7129, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(7130, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(7131, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(7132, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(7133, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(7136, "#DeltaP%P_{calc.} Prof ", 160, 0.3, 1.1, -1., 1.);
  HBProf(7139, "MissingMass%P_{calc.} Prof ", 160, 0.3, 1.1, 1.1, 1.8);

  HB1(8001, "P K18", 800, 1.4, 2.2);
  HB1(8002, "P Kurama", 600, 0, 3);
  HB1(8003, "Charge*Mass", 600, -1., 2.);
  HB2(8004, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(8010, "Closest distance", 50, 0., 50.);
  HB1(8011, "Vertex X", 200, -100, 100);
  HB1(8012, "Vertex Y", 200, -100, 100);
  HB1(8013, "Vertex Z", 200, -200, 200);
  HB1(8014, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB1(8015, "MissingMass K^{+}", 140, 1.1, 1.8);
  HB2(8016, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8017, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8018, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8019, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(8020, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(8021, "MissingMass%U ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(8022, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8023, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8024, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8025, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(8026, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(8027, "MissingMass%V ", 160, -0.4, 0.4, 200, 1., 1.8);
  HB2(8028, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(8029, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(8030, "#DeltaP%P_{calc.} ", 80, 1.1, 1.5, 200, -1., 1.);
  HB2(8031, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(8032, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(8033, "MissingMass%P_{calc.} ", 80, 1.1, 1.5, 200, 1., 1.8);
  HB2(8036, "#DeltaP%P_{calc.} ", 160, 0.3, 1.1, 200, -1., 1.);
  HB2(8039, "MissingMass%P_{calc.} ", 160, 0.3, 1.1, 140, 1.1, 1.8);

  HB2(8041, "K-, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8042, "K-, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8043, "K-, Utgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);
  HB2(8044, "K-, Vtgt w/ vs w/o TPC", 200, -0.1, 0.1, 200, -0.1, 0.1);
  HB2(8045, "K+, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8046, "K+, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8047, "K+, Utgt w/ vs w/o TPC", 200, -0.5, 0.5, 200, -0.5, 0.5);
  HB2(8048, "K+, Vtgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);

  HBProf(8116, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8117, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8118, "#DeltaP%U Prof ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8119, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8120, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8121, "MissingMass%U Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8122, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(8123, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(8124, "#DeltaP%V Prof ", 160, -0.5, 0.5, -1., 1.);
  HBProf(8125, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8126, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8127, "MissingMass%V Prof ", 160, -0.4, 0.4, 1., 1.8);
  HBProf(8128, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(8129, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(8130, "#DeltaP%P_{calc.} Prof ", 160, 1.1, 1.5, -1., 1.);
  HBProf(8131, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(8132, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(8133, "MissingMass%P_{calc.} Prof ", 160, 1.1, 1.5, 1., 1.8);
  HBProf(8136, "#DeltaP%P_{calc.} Prof ", 160, 0.3, 1.1, -1., 1.);
  HBProf(8139, "MissingMass%P_{calc.} Prof ", 160, 0.3, 1.1, 1.1, 1.8);

  HB1(8201, "P K18", 800, 1.4, 2.2);
  HB1(8202, "P Kurama", 600, 0, 3);
  HB1(8203, "Charge*Mass", 600, -1., 2.);
  HB2(8204, "P Kurama%Charge*Mass", 400, -0.8, 2.2, 400, 0, 2.5);
  HB1(8210, "Closest distance", 50, 0., 50.);
  HB1(8211, "Vertex X", 200, -100, 100);
  HB1(8212, "Vertex Y", 200, -100, 100);
  HB1(8213, "Vertex Z", 200, -200, 200);
  HB1(8214, "MissingMass p", 200, 0.3, 1.3);
  HB1(8215, "MissingMass p", 200, 0.3, 1.3);
  HB2(8216, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8217, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8218, "#DeltaP%U ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8219, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8220, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8221, "MissingMass%U ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8222, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8223, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8224, "#DeltaP%V ", 160, -0.4, 0.4, 200, -1., 1.);
  HB2(8225, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8226, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8227, "MissingMass%V ", 160, -0.4, 0.4, 200, 0.3, 1.3);
  HB2(8228, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, 200, -1., 1.);
  HB2(8229, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, 200, -1., 1.);
  HB2(8230, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, 200, -1., 1.);
  HB2(8231, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 200, 0.3, 1.3);
  HB2(8232, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 200, 0.3, 1.3);
  HB2(8233, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 200, 0.3, 1.3);

  HBProf(8316, "#DeltaP%U ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8317, "#DeltaP%U ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8318, "#DeltaP%U ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8319, "MissingMass%U ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8320, "MissingMass%U ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8321, "MissingMass%U ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8322, "#DeltaP%V ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8323, "#DeltaP%V ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8324, "#DeltaP%V ", 160, -0.4, 0.4, -1., 1.);
  HBProf(8325, "MissingMass%V ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8326, "MissingMass%V ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8327, "MissingMass%V ", 160, -0.4, 0.4, 0.3, 1.3);
  HBProf(8328, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, -1., 1.);
  HBProf(8329, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, -1., 1.);
  HBProf(8330, "#DeltaP%P_{calc.} ", 440, 0.3, 2.5, -1., 1.);
  HBProf(8331, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 0.3, 1.3);
  HBProf(8332, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 0.3, 1.3);
  HBProf(8333, "MissingMass%P_{calc.} ", 440, 0.3, 2.5, 0.3, 1.3);

  HB2(8241, "K-, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8242, "K-, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8243, "K-, Utgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);
  HB2(8244, "K-, Vtgt w/ vs w/o TPC", 200, -0.1, 0.1, 200, -0.1, 0.1);
  HB2(8245, "K+, Xtgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8246, "K+, Ytgt w/ vs w/o TPC", 200, -100., 100., 200, -100., 100.);
  HB2(8247, "K+, Utgt w/ vs w/o TPC", 200, -0.5, 0.5, 200, -0.5, 0.5);
  HB2(8248, "K+, Vtgt w/ vs w/o TPC", 200, -0.4, 0.4, 200, -0.4, 0.4);

  HB1(9011, "HTOF tracklen", 140, 1.1, 1.8);
  HB1(10001, "HTOF tracklen [InTgt]", 1000, 0.0, 1000);
  HB1(10002, "HTOF posY     [InTgt]", 1000, -500, 500);
  HB1(10003, "HTOF diffPosY [InTgt]", 1000, -500, 500);    
  HB1(10101, "HTOF tracklen [p_{#Lambda}]", 1000, 0.0, 1000);
  HB1(10102, "HTOF posY     [p_{#Lambda}]", 1000, -500, 500);
  HB1(10103, "HTOF diffPosY [p_{#Lambda}]", 1000, -500, 500);    
  HB1(10201, "HTOF tracklen [#pi_{#Lambda}]", 1000, 0.0, 1000);
  HB1(10202, "HTOF posY     [#pi_{#Lambda}]", 1000, -500, 500);
  HB1(10203, "HTOF diffPosY [#pi_{#Lambda}]", 1000, -500, 500);      
  HB1(10301, "HTOF tracklen [p_{K0}]", 1000, 0.0, 1000);
  HB1(10302, "HTOF posY     [p_{K0}]", 1000, -500, 500);
  HB1(10303, "HTOF diffPosY [p_{K0}]", 1000, -500, 500);        
  HB1(10401, "HTOF tracklen [#pi_{K0}]", 1000, 0.0, 1000);
  HB1(10402, "HTOF posY     [#pi_{K0}]", 1000, -500, 500);
  HB1(10403, "HTOF diffPosY [#pi_{K0}]", 1000, -500, 500);

  //HB1(20010, "p#pi_vtxdist", 600, -300, 300);
  HB1(11011, "vtx_dist [plus]", 600, -300, 300);
  HB1(11012, "vtx_dist [minus]", 600, -300, 300);
  //HB1(20013, "#Lambda_vtxdist", 600, -300, 300);
  //HB1(20201, "Lambda Mass Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", nbinIML,minIML,maxIML);
  HB1(11300, "Mass2 [plus]", 2*nbinM2P, -maxM2P, maxM2P);
  HB1(11301, "Mass2 [#pi^{+}]", 2*nbinM2P, -maxM2P, maxM2P);
  HB1(11302, "Mass2 [K^{+}]", 2*nbinM2P, -maxM2P, maxM2P);
  HB1(11303, "Mass2 [p]", 2*nbinM2P, -maxM2P, maxM2P);
  HB1(11310, "Mass2 [minus]", 2*nbinM2Pi, -maxM2P, maxM2P);
  HB1(11311, "Mass2 [#pi^{-}]", 2*nbinM2P, -maxM2P, maxM2P);
  HB1(11312, "Mass2 [K^{-}]", 2*nbinM2P, -maxM2P, maxM2P);
  HB1(11313, "Mass2 [p^{-}]", 2*nbinM2P, -maxM2P, maxM2P);    
  
  HB1(20010, "p#pi_vtxdist", 600, -300, 300);
  HB1(20011, "vtx_dist [p_{#Lambda}]", 600, -300, 300);
  HB1(20012, "vtx_dist [#pi_{#Lambda}]", 600, -300, 300);
  HB1(20013, "#Lambda_vtxdist", 600, -300, 300);
  HB1(20201, "Lambda Mass Invariant Mass; M_{p#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", nbinIML,minIML,maxIML);
  HB1(20301, "Mass2 [p_{#Lambda}]", 2*nbinM2P, -maxM2P, maxM2P);
  HB1(20302, "Mass2 [#pi_{#Lambda}]", 2*nbinM2Pi, -maxM2Pi, maxM2Pi);

  HB1(30010, "#pi#pi_vtxdist", 600, -300, 300);    
  HB1(30011, "vtx_dist [#pi^{+}_{K0}]", 600, -300, 300);
  HB1(30012, "vtx_dist [#pi^{+}_{K0}]", 600, -300, 300);
  HB1(30013, "K0_vtxdist", 600, -300, 300);  
  HB1(30201, "K0 Mass Invariant Mass; M_{#pi^{+}#pi^{-}} [GeV/#font[12]{c}^{2}]; Counts [/0.002 GeV/#font[12]{c}^{2}]", nbinIMK0,minIMK0,maxIMK0);
  HB1(30301, "Mass2 [#pi^+_{K^0}]",2*nbinM2Pi,-maxM2Pi,maxM2Pi); 
  HB1(30302, "Mass2 [#pi^-_{K^0}]",2*nbinM2Pi,-maxM2Pi,maxM2Pi);
  
#endif
  
  for(int itype=0; itype<kNtype; itype++){//0:general, 1:Lambda reconstruct, 2:K0 reconstruct, 3: K- 
    if(!KPEvent) continue;
    if(itype>pidlikeli::kTypeKm) continue;
    for(int ipid=0; ipid<kNpid+1; ipid++){ // pi,K,p,d,e,all
      for(int icharge=0; icharge<kNchg; icharge++){ // posi,nega
	//if( !( (ipid==1&&icharge==1)||(ipid==2&&icharge==0) ) ) continue; 
        for(int ibe=0; ibe<kNbe; ibe++){ // default: beid=0	  
	  if( !( (itype==pidlikeli::kTypeKm&&ibe!=0) || (itype!=pidlikeli::kTypeKm&&ibe==0) ) ) continue;
	  if( !( ( itype==pidlikeli::kTypeGen&&ipid==pidlikeli::kAllParticles ) // 1
		 || ( itype==pidlikeli::kTypeLmd&&((ipid==pidlikeli::kPion&&icharge==pidlikeli::kMinus)||(ipid==pidlikeli::kProton&&icharge==pidlikeli::kPlus)) ) // 2
		 || ( itype==pidlikeli::kTypeK0&&ipid==pidlikeli::kPion ) // 2
		 || ( itype==pidlikeli::kTypeKm&&ipid==pidlikeli::kKaon&&icharge==pidlikeli::kMinus) //1
		 ) 
	      ) continue; // total=7
	  
          for(int imom=0; imom<kNmom; imom++){
	    if(ipid<kNpid){
	      HB2( (itype+1)*fac_t + ipid*fac_p + icharge*fac_c + ibe*fac_b + imom*fac_m,
		   Form("PDF %s %s %s BE:%s mom=%.4fGeV/c; mass2 ; dEdx",
			type[itype].Data(), plist[ipid].Data(),clist[icharge].Data(),
			pidlikeli::typebe[ibe].Data(),kMomstep*(Double_t(imom))),
		   nbinm2, minm2, maxm2, nbindedx, mindedx, maxdedx);
	    } else {
	      HB2( (itype+1)*fac_t + ipid*fac_p + icharge*fac_c + ibe*fac_b + imom*fac_m,
		   Form("PDF %s %s %s BE:%s mom=%.4fGeV/c; mass2 ; dEdx",
			type[itype].Data(), "all", clist[icharge].Data(),
			pidlikeli::typebe[ibe].Data(),kMomstep*(Double_t(imom))),			
		   nbinm2, minm2, maxm2, nbindedx, mindedx, maxdedx);
	    }
          }
        }
      }
    }
  }

  for(int itype=0; itype<kNtype; itype++){//0:general, 1:Lambda reconstruct, 2:K0 reconstruct, 3: K- 
    if(!KKEvent) continue;
    for(int ipid=0; ipid<kNpid+1; ipid++){ // pi,K,p,d,e,all
      for(int icharge=0; icharge<kNchg; icharge++){ // posi,nega
        for(int ibe=0; ibe<kNbe; ibe++){ // default: beid=0
	  if( !(itype!=pidlikeli::kTypeGen
		&& ( ipid==pidlikeli::kAllParticles||ipid==pidlikeli::kPion||ipid==pidlikeli::kKaon||ipid==pidlikeli::kProton)
		&& (ibe==0)
		) ) continue;
          for(int imom=0; imom<kNmom; imom++){
	    if(ipid<kNpid){
	      HB2( (itype+1)*fac_t + ipid*fac_p + icharge*fac_c + ibe*fac_b + imom*fac_m,
		   Form("PDF %s %s %s BE:%s mom=%.4fGeV/c; mass2 ; dEdx",
			type[itype].Data(), plist[ipid].Data(),clist[icharge].Data(),
			pidlikeli::typebe[ibe].Data(),kMomstep*(Double_t(imom))),
		   nbinm2, minm2, maxm2, nbindedx, mindedx, maxdedx);
	    } else {
	      HB2( (itype+1)*fac_t + ipid*fac_p + icharge*fac_c + ibe*fac_b + imom*fac_m,
		   Form("PDF %s %s %s BE:%s mom=%.4fGeV/c; mass2 ; dEdx",
			type[itype].Data(), "all", clist[icharge].Data(),
			pidlikeli::typebe[ibe].Data(),kMomstep*(Double_t(imom))),			
		   nbinm2, minm2, maxm2, nbindedx, mindedx, maxdedx);
	    }
          }
        }
      }
    }
  }
  

  for(int itype=0; itype<kNtype; itype++){//0:general, 1:Lambda reconstruct, 2:K0 reconstruct, 3: K- 
    if(KPEvent||KKEvent) continue;
    if(itype>0) continue;
    for(int ipid=0; ipid<kNpid+1; ipid++){ // pi,K,p,d,e,all
      for(int icharge=0; icharge<kNchg; icharge++){ // posi,nega
        for(int ibe=0; ibe<kNbe; ibe++){ // default: beid=0
	  if(ibe!=0) continue;
          for(int imom=0; imom<kNmom; imom++){
	    if(ipid<kNpid){
	      HB2( (itype+1)*fac_t + ipid*fac_p + icharge*fac_c + ibe*fac_b + imom*fac_m,
		   Form("PDF %s %s %s BE=%.3fGeV mom=%.4fGeV/c; mass2 ; dEdx",
			type[itype].Data(), plist[ipid].Data(),clist[icharge].Data(),
			minbe+bestep*(double)(ibe),kMomstep*(Double_t(imom))),
		   nbinm2, minm2, maxm2, nbindedx, mindedx, maxdedx);
	    } else {
	      HB2( (itype+1)*fac_t + ipid*fac_p + icharge*fac_c + ibe*fac_b + imom*fac_m,
		   Form("PDF %s %s %s BE=%.3fGeV mom=%.4fGeV/c; mass2 ; dEdx",
			type[itype].Data(), "all", clist[icharge].Data(),
			minbe+bestep*(double)(ibe),kMomstep*(Double_t(imom))),
		   nbinm2, minm2, maxm2, nbindedx, mindedx, maxdedx);
	    }
          }
        }
      }
    }
  }  

  HBTree( "tpc", "tree of E42" );
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
  tree->Branch( "m2Htof", &event.m2Htof );  

  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "nclTpc", &event.nclTpc );
  tree->Branch( "remain_nclTpc", &event.remain_nclTpc );
#if RawCluster
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
  tree->Branch( "ntKuramaCandidate", &event.ntKuramaCandidate );
  tree->Branch( "isKuramaCandidate", &event.isKuramaCandidate );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isXi", &event.isXi );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "isMultiloop", &event.isMultiloop );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "pid", &event.pid );
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
  tree->Branch( "isElectron", &event.isElectron );
  tree->Branch( "nsigma_triton", &event.nsigma_triton );
  tree->Branch( "nsigma_deutron", &event.nsigma_deutron );
  tree->Branch( "nsigma_proton", &event.nsigma_proton );
  tree->Branch( "nsigma_kaon", &event.nsigma_kaon );
  tree->Branch( "nsigma_pion", &event.nsigma_pion );
  tree->Branch( "nsigma_electron", &event.nsigma_electron );

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
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);

  tree->Branch( "chargeIndistinguishable", &event.chargeIndistinguishable );
  tree->Branch( "chisqr_inverted", &event.chisqr_inverted );
  tree->Branch( "pval_inverted", &event.pval_inverted );
  tree->Branch( "helix_cx_inverted", &event.helix_cx_inverted );
  tree->Branch( "helix_cy_inverted", &event.helix_cy_inverted );
  tree->Branch( "helix_z0_inverted", &event.helix_z0_inverted );
  tree->Branch( "helix_r_inverted", &event.helix_r_inverted );
  tree->Branch( "helix_dz_inverted", &event.helix_dz_inverted );
  tree->Branch( "mom0_inverted", &event.mom0_inverted );
  tree->Branch( "pid_inverted", &event.pid_inverted );

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
  tree->Branch( "clusteredVtx_x", &event.clusteredVtx_x );
  tree->Branch( "clusteredVtx_y", &event.clusteredVtx_y );
  tree->Branch( "clusteredVtx_z", &event.clusteredVtx_z );
  tree->Branch( "clusteredVtxid", &event.clusteredVtxid );

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "p_3rd" , &event.p_3rd);
  tree->Branch( "xoutK18", &event.xoutK18);
  tree->Branch( "youtK18", &event.youtK18);
  tree->Branch( "uoutK18", &event.uoutK18);
  tree->Branch( "voutK18", &event.voutK18);
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "utgtK18", &event.utgtK18);
  tree->Branch( "vtgtK18", &event.vtgtK18);
  tree->Branch( "thetaK18", &event.thetaK18);
  tree->Branch( "xhtofK18", &event.xhtofK18);
  tree->Branch( "yhtofK18", &event.yhtofK18);
  tree->Branch( "xvpHS", &event.xvpHS);
  tree->Branch( "yvpHS", &event.yvpHS);
  tree->Branch( "zvpHS", &event.zvpHS);
  tree->Branch( "xtgtHS", &event.xtgtHS);
  tree->Branch( "ytgtHS", &event.ytgtHS);
  tree->Branch( "ztgtHS", &event.ztgtHS);
  tree->Branch( "layerK18", &event.layerK18);
  tree->Branch( "wireK18", &event.wireK18);
  tree->Branch( "localhitposK18", &event.localhitposK18);
  tree->Branch( "wposK18", &event.wposK18);

  tree->Branch( "tpcidTPCK18", &event.tpcidTPCK18);
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

  tree->Branch( "ntKurama",     &event.ntKurama);
  tree->Branch( "chisqrKurama", &event.chisqrKurama);
  tree->Branch( "pKurama",      &event.pKurama);
  tree->Branch( "qKurama",      &event.qKurama);
  tree->Branch( "m2",           &event.m2);
  tree->Branch( "m2Org",        &event.m2Org);
  tree->Branch( "xtgtKurama",   &event.xtgtKurama);
  tree->Branch( "ytgtKurama",   &event.ytgtKurama);
  tree->Branch( "utgtKurama",   &event.utgtKurama);
  tree->Branch( "vtgtKurama",   &event.vtgtKurama);
  tree->Branch( "thetaKurama",  &event.thetaKurama);
  tree->Branch( "pathKurama",  &event.pathKurama);
  tree->Branch( "cstof",  &event.cstof);
  tree->Branch( "pathwcKurama",  &event.pathwcKurama);
  tree->Branch( "xvpKurama",    &event.xvpKurama);
  tree->Branch( "yvpKurama",    &event.yvpKurama);
  tree->Branch( "zvpKurama",    &event.zvpKurama);
  tree->Branch( "xin",  &event.xin);
  tree->Branch( "yin",  &event.yin);
  tree->Branch( "zin",  &event.zin);
  tree->Branch( "pxin",  &event.pxin);
  tree->Branch( "pyin",  &event.pyin);
  tree->Branch( "pzin",  &event.pzin);
  tree->Branch( "xout",  &event.xout);
  tree->Branch( "yout",  &event.yout);
  tree->Branch( "zout",  &event.zout);
  tree->Branch( "pxout",  &event.pxout);
  tree->Branch( "pyout",  &event.pyout);
  tree->Branch( "pzout",  &event.pzout);
  tree->Branch( "layer",  &event.layer);
  tree->Branch( "wire",  &event.wire);
  tree->Branch( "localhitpos",  &event.localhitpos);
  tree->Branch( "wpos",  &event.wpos);

  tree->Branch( "xhtofKurama", &event.xhtofKurama);
  tree->Branch( "yhtofKurama", &event.yhtofKurama);
  tree->Branch( "tpcidTPCKurama", &event.tpcidTPCKurama);
  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  //tree->Branch( "piflagTPCKurama", &event.piflagTPCKurama);
  tree->Branch( "kflagTPCKurama", &event.kflagTPCKurama);
  tree->Branch( "pflagTPCKurama", &event.pflagTPCKurama);
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
  tree->Branch( "lgasvesselTPCKurama", &event.lgasvesselTPCKurama);
  tree->Branch( "xgasvesselTPCKurama", &event.xgasvesselTPCKurama);
  tree->Branch( "ygasvesselTPCKurama", &event.ygasvesselTPCKurama);
  tree->Branch( "lvpTPCKurama", &event.lvpTPCKurama);
  tree->Branch( "xvpTPCKurama", &event.xvpTPCKurama);
  tree->Branch( "yvpTPCKurama", &event.yvpTPCKurama);

  tree->Branch("nKm",           &event.nKm);
  tree->Branch("nKp",           &event.nKp);
  tree->Branch("nKK",           &event.nKK);
  tree->Branch("vtx",           &event.vtx);
  tree->Branch("vty",           &event.vty);
  tree->Branch("vtz",           &event.vtz);
  tree->Branch("closeDist",     &event.closeDist);
  tree->Branch("inside",        &event.inside);
  tree->Branch("MissMass",      &event.MissMass);
  tree->Branch("MissMassCorr",  &event.MissMassCorr);
  tree->Branch("MissMassCorrDE", &event.MissMassCorrDE);
  tree->Branch("pOrg",       &event.pOrg);
  tree->Branch("pCalc",      &event.pCalc);
  tree->Branch("pCorr",      &event.pCorr);
  tree->Branch("pCorrDE",    &event.pCorrDE);
  tree->Branch("xb",         &event.xkm);
  tree->Branch("yb",         &event.ykm);
  tree->Branch("ub",         &event.ukm);
  tree->Branch("vb",         &event.vkm);
  tree->Branch("xs",         &event.xkp);
  tree->Branch("ys",         &event.ykp);
  tree->Branch("us",         &event.ukp);
  tree->Branch("vs",         &event.vkp);
  tree->Branch("Kflag",      &event.Kflag);
  tree->Branch("Pflag",      &event.Pflag);
  tree->Branch("Heavyflag",  &event.Heavyflag);

  tree->Branch( "isgoodTPC", &event.isgoodTPC);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "vtxTPC", &event.vtxTPC);
  tree->Branch( "vtyTPC", &event.vtyTPC);
  tree->Branch( "vtzTPC", &event.vtzTPC);
  tree->Branch( "closeDistTPC", &event.closeDistTPC);
  tree->Branch( "MissMassTPC", &event.MissMassTPC);
  tree->Branch( "MissMassCorrTPC", &event.MissMassCorrTPC);
  tree->Branch( "MissMassCorrDETPC", &event.MissMassCorrDETPC);
  tree->Branch( "MissMassNuclTPC", &event.MissMassNuclTPC);
  tree->Branch( "MissMassNuclCorrTPC", &event.MissMassNuclCorrTPC);
  tree->Branch( "MissMassNuclCorrDETPC", &event.MissMassNuclCorrDETPC);
  tree->Branch( "pOrgTPC", &event.pOrgTPC);
  tree->Branch( "pCorrTPC", &event.pCorrTPC);
  tree->Branch( "pCorrDETPC", &event.pCorrDETPC);
  tree->Branch( "pCalcTPC", &event.pCalcTPC);
  tree->Branch( "thetaCMTPC", &event.thetaCMTPC);
  tree->Branch( "costCMTPC", &event.costCMTPC);
  tree->Branch( "pCalcDETPC", &event.pCalcDETPC);
  tree->Branch( "thetaCMDETPC", &event.thetaCMDETPC);
  tree->Branch( "costCMDETPC", &event.costCMDETPC);
  tree->Branch( "xistarpCalcDETPC", &event.xistarpCalcDETPC);
  tree->Branch( "xistarthetaCMDETPC", &event.xistarthetaCMDETPC);
  tree->Branch( "xistarcostCMDETPC", &event.xistarcostCMDETPC);
  tree->Branch( "kpscatpCalcTPC", &event.kpscatpCalcTPC);
  tree->Branch( "kpscatthetaCMTPC", &event.kpscatthetaCMTPC);
  tree->Branch( "kpscatcostCMTPC", &event.kpscatcostCMTPC);
  tree->Branch( "kpscatpCalcDETPC", &event.kpscatpCalcDETPC);
  tree->Branch( "kpscatthetaCMDETPC", &event.kpscatthetaCMDETPC);
  tree->Branch( "kpscatcostCMDETPC", &event.kpscatcostCMDETPC);
  tree->Branch( "thetaTPC", &event.thetaTPC);
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);

  tree->Branch("insideTgt", &event.insideTgt);
  tree->Branch("extrapolateTgt", &event.extrapTgt);
  tree->Branch("extrapolateVtx", &event.extrapVtx);
  tree->Branch("ntTpcInside", &event.ntTpc_inside);  
  tree->Branch("m2HtofVtx", &event.m2HtofVtx);
  tree->Branch("HitIdHtof", &event.hitidHtof);
  tree->Branch("SegHtof", &event.segHtof);
  tree->Branch("TofHtof", &event.tofHtof);
  tree->Branch("InvBetaHtof", &event.invbetaHtof);  
  tree->Branch("ProdVtx_x", &event.prodvtx_x);  
  tree->Branch("ProdVtx_y", &event.prodvtx_y);
  tree->Branch("ProdVtx_z", &event.prodvtx_z);
  tree->Branch("posHtof_x", &event.posHtof_x);    
  tree->Branch("posHtof_y", &event.posHtof_y);
  tree->Branch("posHtof_z", &event.posHtof_z);
  tree->Branch("TrackLenHtofTgt", &event.tracklenHtof);
  tree->Branch("DistFromVtxTgt", &event.distVtx);
  
  tree->Branch("Lflag", &event.lflag);
  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom", &event.lmom);  
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  tree->Branch("LambdaPPiAngle", &event.ppiangle);
  //tree->Branch("LambdaTrackLen", &event.ldecays_tracklen);              
  tree->Branch("LDecaysTrackId", &event.ldecays_id);
  tree->Branch("LDecaysMom", &event.ldecays_mom);
  tree->Branch("LDecaysMom_x", &event.ldecays_mom_x);
  tree->Branch("LDecaysMom_y", &event.ldecays_mom_y);
  tree->Branch("LDecaysMom_z", &event.ldecays_mom_z);
  tree->Branch("LDecaysHtofExtrapolate", &event.ldecays_htofextrap);  
  tree->Branch("LDecaysHtofHitId", &event.ldecays_htofhitid);
  tree->Branch("LDecaysHtofSeg", &event.ldecays_htofseg);
  tree->Branch("LDecaysTrackLen", &event.ldecays_tracklen);
  tree->Branch("LDecaysInvBeta", &event.ldecays_invbeta);
  tree->Branch("LDecaysMass2", &event.ldecays_mass2);  
  tree->Branch("LDecaysHtofPos_x", &event.ldecays_htofpos_x);
  tree->Branch("LDecaysHtofPos_y", &event.ldecays_htofpos_y);
  tree->Branch("LDecaysHtofPos_z", &event.ldecays_htofpos_z);
  
  tree->Branch("K0flag", &event.k0flag);
  tree->Branch("K0Mass", &event.k0mass);
  tree->Branch("K0DecayVtx_x", &event.k0decayvtx_x);
  tree->Branch("K0DecayVtx_y", &event.k0decayvtx_y);
  tree->Branch("K0DecayVtx_z", &event.k0decayvtx_z);
  tree->Branch("K0Mom_x", &event.k0mom_x);
  tree->Branch("K0Mom_y", &event.k0mom_y);
  tree->Branch("K0Mom_z", &event.k0mom_z);
  tree->Branch("K0VtxCloseDist", &event.pipi_dist);
  tree->Branch("K0PPiAngle", &event.pipiangle);  
  tree->Branch("K0DecaysTrackId", &event.k0decays_id);
  tree->Branch("K0DecaysMom", &event.k0decays_mom);
  tree->Branch("K0DecaysMom_x", &event.k0decays_mom_x);
  tree->Branch("K0DecaysMom_y", &event.k0decays_mom_y);
  tree->Branch("K0DecaysMom_z", &event.k0decays_mom_z);
  tree->Branch("K0DecaysHtofExtrapolate", &event.k0decays_htofextrap);  
  tree->Branch("K0DecaysHtofPos_x", &event.k0decays_htofpos_x);
  tree->Branch("K0DecaysHtofPos_y", &event.k0decays_htofpos_y);
  tree->Branch("K0DecaysHtofPos_z", &event.k0decays_htofpos_z);
  tree->Branch("K0DecaysMass2", &event.k0decays_mass2);
  tree->Branch("K0DecaysTrackLen", &event.k0decays_tracklen);          
  
  // tree->Branch("K0VtxCloseDist", &event.GFk0pipi_dist);
  // tree->Branch("K0PiPiAngle", &event.pipiangle);  
  tree->Branch("K0DecaysTrackId", &event.k0decays_id);
  tree->Branch("K0DecaysMom", &event.k0decays_mom);
  tree->Branch("K0DecaysMom_x", &event.k0decays_mom_x);
  tree->Branch("K0DecaysMom_y", &event.k0decays_mom_y);
  tree->Branch("K0DecaysMom_z", &event.k0decays_mom_z);
  tree->Branch("K0DecaysHtofPos_x", &event.k0decays_htofpos_x);
  tree->Branch("K0DecaysHtofPos_y", &event.k0decays_htofpos_y);
  tree->Branch("K0DecaysHtofPos_z", &event.k0decays_htofpos_z);      
  tree->Branch("K0DecaysMass2", &event.k0decays_mass2);
  tree->Branch("K0DecaysTrackLen", &event.k0decays_tracklen);        

  TTreeReaderCont[kTpc] = new TTreeReader( "tpc", TFileCont[kTpc] );
  const auto& reader = TTreeReaderCont[kTpc];

  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );
  src.nhTpc = new TTreeReaderValue<Int_t>( *reader, "nhTpc" );
  src.nclTpc = new TTreeReaderValue<Int_t>( *reader, "nclTpc" );
  src.remain_nclTpc = new TTreeReaderValue<Int_t>( *reader, "remain_nclTpc" );
#if RawCluster
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
  src.ntKuramaCandidate = new TTreeReaderValue<Int_t>( *reader, "ntKuramaCandidate" );
  src.isKuramaCandidate = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKuramaCandidate" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isXi = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isXi" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
  src.isMultiloop = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isMultiloop" );

  src.charge = new TTreeReaderValue<std::vector<Int_t>>( *reader, "charge" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid" );
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
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_size" );
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_mrow" );
  src.track_cluster_de_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_de_center" );
  src.track_cluster_x_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_x_center" );
  src.track_cluster_y_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_y_center" );
  src.track_cluster_z_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_z_center" );
  src.track_cluster_row_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "track_cluster_row_center" );

  src.chargeIndistinguishable = new TTreeReaderValue<std::vector<Int_t>>( *reader, "chargeIndistinguishable" );
  src.pid_inverted = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pid_inverted" );
  src.chisqr_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "chisqr_inverted" );
  src.pval_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pval_inverted" );
  src.helix_cx_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cx_inverted" );
  src.helix_cy_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_cy_inverted" );
  src.helix_z0_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_z0_inverted" );
  src.helix_r_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_r_inverted" );
  src.helix_dz_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "helix_dz_inverted" );
  src.mom0_inverted = new TTreeReaderValue<std::vector<Double_t>>( *reader, "mom0_inverted" );

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

  src.nvtxTpcClustered = new TTreeReaderValue<Int_t>(*reader,"nvtxTpcClustered");
  src.clusteredVtx_x = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_x" );
  src.clusteredVtx_y = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_y" );
  src.clusteredVtx_z = new TTreeReaderValue<std::vector<Double_t>>( *reader, "clusteredVtx_z" );
  src.clusteredVtxid = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "clusteredVtxid" );

  src.ntTPCK18 = new TTreeReaderValue<Int_t>( *reader, "ntK18" );
  src.tpcidTPCK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "tpcidTPCK18" );
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

  src.ntTPCKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
  src.tpcidTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "tpcidTPCKurama" );
  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
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
  src.lgasvesselTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "lgasvesselTPCKurama" );
  src.xgasvesselTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xgasvesselTPCKurama" );
  src.ygasvesselTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ygasvesselTPCKurama" );
  src.lvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "lvpTPCKurama" );
  src.xvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpTPCKurama" );
  src.yvpTPCKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpTPCKurama" );
  src.xhtofKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofKurama" );
  src.yhtofKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofKurama" );

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
  src.kflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "kflagTPCKurama" );
  src.pflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pflagTPCKurama" );    
  src.pOrgTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pOrgTPC" );
  src.pCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrTPC" );
  src.pCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrDETPC" );
  src.pCalcTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCalcTPC" );
  src.thetaCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaCMTPC" );
  src.costCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "costCMTPC" );
  src.pCalcDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCalcDETPC" );
  src.thetaCMDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaCMDETPC" );
  src.costCMDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "costCMDETPC" );
  src.xistarpCalcDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xistarpCalcDETPC" );
  src.xistarthetaCMDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xistarthetaCMDETPC" );
  src.xistarcostCMDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xistarcostCMDETPC" );
  src.kpscatpCalcTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "kpscatpCalcTPC" );
  src.kpscatthetaCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "kpscatthetaCMTPC" );
  src.kpscatcostCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "kpscatcostCMTPC" );
  src.kpscatpCalcDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "kpscatpCalcDETPC" );
  src.kpscatthetaCMDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "kpscatthetaCMDETPC" );
  src.kpscatcostCMDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "kpscatcostCMDETPC" );
  src.thetaTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPC" );
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );

  //For Kurama + K18 part in G4input
  src.p_3rd = new TTreeReaderValue<std::vector<Double_t>>( *reader, "p_3rd" );
  src.xoutK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xoutK18" );
  src.youtK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "youtK18" );
  src.uoutK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "uoutK18" );
  src.voutK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "voutK18" );
  src.xvpHS  =  new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpHS" );
  src.yvpHS  =  new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpHS" );
  src.zvpHS  =  new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "zvpHS" );
  src.xtgtHS = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtHS" );
  src.ytgtHS = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtHS" );
  src.ztgtHS = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ztgtHS" );
  src.layerK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "layerK18" );
  src.wireK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "wireK18" );
  src.localhitposK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "localhitposK18" );
  src.wposK18 = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "wposK18" );

  src.xvpKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "xvpKurama" );
  src.yvpKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "yvpKurama" );
  src.zvpKurama = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "zvpKurama" );
  src.layer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "layer" );
  src.wire = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "wire" );
  src.localhitpos = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "localhitpos" );
  src.wpos = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "wpos" );

  TTreeCont[kKScat]->SetBranchStatus("*", 0);
  TTreeCont[kKScat]->SetBranchStatus("ntK18",          1);
  TTreeCont[kKScat]->SetBranchStatus("chisqrK18",      1);
  TTreeCont[kKScat]->SetBranchStatus("pK18",           1);
  TTreeCont[kKScat]->SetBranchStatus("xtgtK18",        1);
  TTreeCont[kKScat]->SetBranchStatus("ytgtK18",        1);
  TTreeCont[kKScat]->SetBranchStatus("utgtK18",        1);
  TTreeCont[kKScat]->SetBranchStatus("vtgtK18",        1);

  TTreeCont[kKScat]->SetBranchStatus("ntKurama",       1);
  TTreeCont[kKScat]->SetBranchStatus("chisqrKurama",   1);
  TTreeCont[kKScat]->SetBranchStatus("pKurama",        1);
  TTreeCont[kKScat]->SetBranchStatus("qKurama",        1);
  TTreeCont[kKScat]->SetBranchStatus("m2",             1);
  TTreeCont[kKScat]->SetBranchStatus("m2Org",          1);
  TTreeCont[kKScat]->SetBranchStatus("xtgtKurama",     1);
  TTreeCont[kKScat]->SetBranchStatus("ytgtKurama",     1);
  TTreeCont[kKScat]->SetBranchStatus("utgtKurama",     1);
  TTreeCont[kKScat]->SetBranchStatus("vtgtKurama",     1);
  TTreeCont[kKScat]->SetBranchStatus("thetaKurama",    1);
  TTreeCont[kKScat]->SetBranchStatus("tofsegKurama",   1);
  TTreeCont[kKScat]->SetBranchStatus("path",           1);
  TTreeCont[kKScat]->SetBranchStatus("pathwcKurama",   1);
  TTreeCont[kKScat]->SetBranchStatus("xout",        1);
  TTreeCont[kKScat]->SetBranchStatus("yout",        1);
  TTreeCont[kKScat]->SetBranchStatus("zout",        1);
  TTreeCont[kKScat]->SetBranchStatus("pxout",       1);
  TTreeCont[kKScat]->SetBranchStatus("pyout",       1);
  TTreeCont[kKScat]->SetBranchStatus("pzout",       1);
  TTreeCont[kKScat]->SetBranchStatus("xin",        1);
  TTreeCont[kKScat]->SetBranchStatus("yin",        1);
  TTreeCont[kKScat]->SetBranchStatus("zin",        1);
  TTreeCont[kKScat]->SetBranchStatus("pxin",       1);
  TTreeCont[kKScat]->SetBranchStatus("pyin",       1);
  TTreeCont[kKScat]->SetBranchStatus("pzin",       1);

  TTreeCont[kKScat]->SetBranchStatus("cstof",          1);
  TTreeCont[kKScat]->SetBranchStatus("nhHtof",         1);
  TTreeCont[kKScat]->SetBranchStatus("HtofSeg",        1);
  TTreeCont[kKScat]->SetBranchStatus("tHtof",          1);
  TTreeCont[kKScat]->SetBranchStatus("dtHtof",         1);
  TTreeCont[kKScat]->SetBranchStatus("deHtof",         1);
  TTreeCont[kKScat]->SetBranchStatus("posHtof",        1);

  TTreeCont[kKScat]->SetBranchStatus("nKm",            1);
  TTreeCont[kKScat]->SetBranchStatus("nKp",            1);
  TTreeCont[kKScat]->SetBranchStatus("nKK",            1);
  TTreeCont[kKScat]->SetBranchStatus("inside",         1);
  TTreeCont[kKScat]->SetBranchStatus("vtx",            1);
  TTreeCont[kKScat]->SetBranchStatus("vty",            1);
  TTreeCont[kKScat]->SetBranchStatus("vtz",            1);
  TTreeCont[kKScat]->SetBranchStatus("closeDist",      1);
  TTreeCont[kKScat]->SetBranchStatus("MissMass",       1);
  TTreeCont[kKScat]->SetBranchStatus("MissMassCorr",   1);
  TTreeCont[kKScat]->SetBranchStatus("MissMassCorrDE", 1);
  TTreeCont[kKScat]->SetBranchStatus("pOrg",           1);
  TTreeCont[kKScat]->SetBranchStatus("pCalc",          1);
  TTreeCont[kKScat]->SetBranchStatus("pCorr",          1);
  TTreeCont[kKScat]->SetBranchStatus("pCorrDE",        1);
  TTreeCont[kKScat]->SetBranchStatus("xkm",            1);
  TTreeCont[kKScat]->SetBranchStatus("ykm",            1);
  TTreeCont[kKScat]->SetBranchStatus("ukm",            1);
  TTreeCont[kKScat]->SetBranchStatus("vkm",            1);
  TTreeCont[kKScat]->SetBranchStatus("xkp",            1);
  TTreeCont[kKScat]->SetBranchStatus("ykp",            1);
  TTreeCont[kKScat]->SetBranchStatus("ukp",            1);
  TTreeCont[kKScat]->SetBranchStatus("vkp",            1);
  TTreeCont[kKScat]->SetBranchStatus("Kflag",          1);
  TTreeCont[kKScat]->SetBranchStatus("Pflag",          1);
  TTreeCont[kKScat]->SetBranchStatus("Heavyflag",      1);

  TTreeCont[kKScat]->SetBranchAddress("ntK18",    &src.ntK18);
  TTreeCont[kKScat]->SetBranchAddress("chisqrK18", src.chisqrK18);
  TTreeCont[kKScat]->SetBranchAddress("pK18",      src.pK18);
  TTreeCont[kKScat]->SetBranchAddress("xtgtK18",   src.xtgtK18);
  TTreeCont[kKScat]->SetBranchAddress("ytgtK18",   src.ytgtK18);
  TTreeCont[kKScat]->SetBranchAddress("utgtK18",   src.utgtK18);
  TTreeCont[kKScat]->SetBranchAddress("vtgtK18",   src.vtgtK18);
  TTreeCont[kKScat]->SetBranchAddress("thetaK18",  src.thetaK18);

  TTreeCont[kKScat]->SetBranchAddress("ntKurama",    &src.ntKurama);
  TTreeCont[kKScat]->SetBranchAddress("chisqrKurama", src.chisqrKurama);
  TTreeCont[kKScat]->SetBranchAddress("pKurama",      src.pKurama);
  TTreeCont[kKScat]->SetBranchAddress("qKurama",      src.qKurama);
  TTreeCont[kKScat]->SetBranchAddress("m2",           src.m2);
  TTreeCont[kKScat]->SetBranchAddress("m2Org",        src.m2Org);
  TTreeCont[kKScat]->SetBranchAddress("xtgtKurama",   src.xtgtKurama);
  TTreeCont[kKScat]->SetBranchAddress("ytgtKurama",   src.ytgtKurama);
  TTreeCont[kKScat]->SetBranchAddress("utgtKurama",   src.utgtKurama);
  TTreeCont[kKScat]->SetBranchAddress("vtgtKurama",   src.vtgtKurama);
  TTreeCont[kKScat]->SetBranchAddress("thetaKurama",  src.thetaKurama);
  TTreeCont[kKScat]->SetBranchAddress("tofsegKurama", src.tofsegKurama);
  TTreeCont[kKScat]->SetBranchAddress("path",         src.pathKurama);
  TTreeCont[kKScat]->SetBranchAddress("pathwcKurama", src.pathwcKurama);
  TTreeCont[kKScat]->SetBranchAddress("cstof",        src.cstof);
  TTreeCont[kKScat]->SetBranchAddress("xout",  src.xout);
  TTreeCont[kKScat]->SetBranchAddress("yout",  src.yout);
  TTreeCont[kKScat]->SetBranchAddress("zout",  src.zout);
  TTreeCont[kKScat]->SetBranchAddress("pxout", src.pxout);
  TTreeCont[kKScat]->SetBranchAddress("pyout", src.pyout);
  TTreeCont[kKScat]->SetBranchAddress("pzout", src.pzout);
  TTreeCont[kKScat]->SetBranchAddress("xin",  src.xin);
  TTreeCont[kKScat]->SetBranchAddress("yin",  src.yin);
  TTreeCont[kKScat]->SetBranchAddress("zin",  src.zin);
  TTreeCont[kKScat]->SetBranchAddress("pxin", src.pxin);
  TTreeCont[kKScat]->SetBranchAddress("pyin", src.pyin);
  TTreeCont[kKScat]->SetBranchAddress("pzin", src.pzin);

  TTreeCont[kKScat]->SetBranchAddress("nhHtof", &src.nhHtof);
  TTreeCont[kKScat]->SetBranchAddress("HtofSeg", src.HtofSeg);
  TTreeCont[kKScat]->SetBranchAddress("tHtof", src.tHtof);
  TTreeCont[kKScat]->SetBranchAddress("dtHtof", src.dtHtof);
  TTreeCont[kKScat]->SetBranchAddress("deHtof", src.deHtof);
  TTreeCont[kKScat]->SetBranchAddress("posHtof", src.posHtof);

  TTreeCont[kKScat]->SetBranchAddress("nKm",       &src.nKm);
  TTreeCont[kKScat]->SetBranchAddress("nKp",       &src.nKp);
  TTreeCont[kKScat]->SetBranchAddress("nKK",       &src.nKK);
  TTreeCont[kKScat]->SetBranchAddress("inside",    src.inside);
  TTreeCont[kKScat]->SetBranchAddress("vtx",       src.vtx);
  TTreeCont[kKScat]->SetBranchAddress("vty",       src.vty);
  TTreeCont[kKScat]->SetBranchAddress("vtz",       src.vtz);
  TTreeCont[kKScat]->SetBranchAddress("closeDist", src.closeDist);
  TTreeCont[kKScat]->SetBranchAddress("MissMass",  src.MissMass);
  TTreeCont[kKScat]->SetBranchAddress("MissMassCorr", src.MissMassCorr);
  TTreeCont[kKScat]->SetBranchAddress("MissMassCorrDE", src.MissMassCorrDE);
  TTreeCont[kKScat]->SetBranchAddress("pOrg",      src.pOrg);
  TTreeCont[kKScat]->SetBranchAddress("pCalc",     src.pCalc);
  TTreeCont[kKScat]->SetBranchAddress("pCorr",     src.pCorr);
  TTreeCont[kKScat]->SetBranchAddress("pCorrDE",   src.pCorrDE);
  TTreeCont[kKScat]->SetBranchAddress("xkm",       src.xkm);
  TTreeCont[kKScat]->SetBranchAddress("ykm",       src.ykm);
  TTreeCont[kKScat]->SetBranchAddress("ukm",       src.ukm);
  TTreeCont[kKScat]->SetBranchAddress("vkm",       src.vkm);
  TTreeCont[kKScat]->SetBranchAddress("xkp",       src.xkp);
  TTreeCont[kKScat]->SetBranchAddress("ykp",       src.ykp);
  TTreeCont[kKScat]->SetBranchAddress("ukp",       src.ukp);
  TTreeCont[kKScat]->SetBranchAddress("vkp",       src.vkp);
  TTreeCont[kKScat]->SetBranchAddress("Kflag",     src.Kflag);
  TTreeCont[kKScat]->SetBranchAddress("Pflag",     src.Pflag);
  TTreeCont[kKScat]->SetBranchAddress("Heavyflag", src.Heavyflag);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&      
      InitializeParameter<UserParamMan>("USER") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
