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
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "PidLikelihoodMan.hh"
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
auto&       gPidLike = PidLikelihoodMan::GetInstance();  
const auto& zHSCenter = gGeom.LocalZ("HS");
const Double_t truncatedMean = 0.8; //80%
const Double_t dist_cut = 5;
const Double_t chisqrcut = 100;
const Double_t decutppi = 2.0;
const Double_t gfposycut = 100.0;  

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");
const auto& psTrigA = ConfMan::Get<Int_t>("PSTRGA");
const auto& psTrigB = ConfMan::Get<Int_t>("PSTRGB");    

const int kPmin = pidlikeli::kPmin;
const int kPmax = pidlikeli::kPmax;
const int kMomstep = pidlikeli::kDP;
const int kNmom = pidlikeli::kNmom;
const int kNpid = pidlikeli::kNpid; 
const int kNchg = pidlikeli::kNchg;
const int kNtype = pidlikeli::kNtype;
const int kNbe = pidlikeli::kNbe;
const auto& plist = pidlikeli::plist;
const auto& clist = pidlikeli::clist;
const auto& type = pidlikeli::type;
const int kPidPi = static_cast<int>(pidlikeli::Pid::Pi);
const int kPidK  = static_cast<int>(pidlikeli::Pid::K);
const int kPidP  = static_cast<int>(pidlikeli::Pid::P);    
bool debugflag=true;

const Double_t lambda_masscut = 0.1;
const Double_t lambda_masscut_final = 0.02; //final  
  
const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;

const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t e_vtx_distcut = 300;
const Double_t ppi_distcut = 10.; //ref

const Double_t residual_track_distcut = 25.;
const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
  
bool use_pidlikeli = false;
const Double_t min_mass2_p = 0.20;
const Double_t max_mass2_p = 1.25;
const Double_t min_mass2_pi = -0.05;
const Double_t max_mass2_pi =  0.05;
const Double_t min_mass2_k  = 0.10;
const Double_t max_mass2_k  = 0.40 ;
  
const Double_t lmass_cut_min  = 1.105;
const Double_t lmass_cut_max  = 1.125;
const Double_t lctau_cut_min  = 3.0; //cm

  const int hid_inc = 90000;   //Lmd:0,K0:0,gamma:0  
  const int hid_exc = 100000;   //Lmd:0,K0:0,gamma:0
  const int hid_exc_Lmd = 200000;   //Lmd:1,K0:0,gamma:0
  const int hid_exc_gamma = 500000;   //Lmd:0,K0:0,gamma:1
  const int hid_semiexc = 800000;   //Lmd:0,K0:0,gamma:0  
}

//enum class Species : int { PPIP, EP, EM, P, KP, KM, PIP, PIM, N_SPECIES };
enum class Species : int { PIM, PIP, KM, KP, P, EM, EP, PPIP, N_SPECIES };
using Counts = std::array<int, static_cast<int>(Species::N_SPECIES)>;

static constexpr uint32_t BASE_TRACK = 5;

static constexpr std::array<uint32_t, static_cast<int>(Species::N_SPECIES)> W_TRACK = []{
  std::array<uint32_t, static_cast<int>(Species::N_SPECIES)> w{};
  uint32_t f = 1;
  for (int i=0;i<(int)Species::N_SPECIES;++i){ w[i]=f; f*=BASE_TRACK; }
  return w;
}();

inline uint32_t encode_pid_code(const Counts& n) {
    uint32_t code = 0;
    for (int i=0;i<(int)Species::N_SPECIES;++i){
        int nk = n[i];
        if (nk < 0) nk = 0;
        if (nk >= (int)BASE_TRACK) nk = BASE_TRACK-1;
        code += (uint32_t)nk * W_TRACK[i];
    }
    return code;
}

inline void decode_pid_code(uint32_t code, Counts& n_out){
  for (int i=0;i<(int)Species::N_SPECIES;++i){ n_out[i] = code % BASE_TRACK; code /= BASE_TRACK; }
}

enum class Reco : int { LMD, K0, GAMMA, N_RECO };
using RecoCounts = std::array<int, static_cast<int>(Reco::N_RECO)>;

static constexpr std::array<int, (int)Reco::N_RECO> BASES_RECO = { 2, 2, 3 }; // Λ×K0×γ → 2×2×3 = 12 通り

static constexpr std::array<int, (int)Reco::N_RECO> W_RECO = []{
    std::array<int, (int)Reco::N_RECO> w{};
    int f = 1;
    for (int i=0;i<(int)Reco::N_RECO;++i){ w[i]=f; f*=BASES_RECO[i]; }
    return w;
}();

inline std::string decimal_to_base5(uint64_t value) {
    if (value == 0) return "0";
    std::string out;
    while (value > 0) {
        int digit = value % 5;
        out.insert(out.begin(), '0' + digit);
        value /= 5;
    }
    return out;
}

inline uint64_t base5_to_decimal(const std::string& base5str) {
    uint64_t value = 0;
    for (char c : base5str) {
        if (c < '0' || c > '4')
            throw std::runtime_error("invalid digit in base5 string");
        int digit = c - '0';
        value = value * 5 + digit;
    }
    return value;
}

inline int encode_reco_code(const RecoCounts& r){
    int code = 0;
    for (int i=0;i<(int)Reco::N_RECO;++i){
        int ni = std::clamp(r[i], 0, BASES_RECO[i]-1);
        code += ni * W_RECO[i];
    }
    return code;
}

inline void decode_reco_code(int code, RecoCounts& r_out){
    for (int i=0;i<(int)Reco::N_RECO;++i){
        r_out[i] = (code / W_RECO[i]) % BASES_RECO[i];
    }
}

inline int encode_reco_id_1based(const RecoCounts& r){ return encode_reco_code(r) + 1; }

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kE42, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[GenfitE42]", "[OutFile]" };
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

  Int_t ntK18;
  std::vector<Double_t> pK18;
  std::vector<Double_t> chisqrK18;
  std::vector<Double_t> xtgtK18;
  std::vector<Double_t> ytgtK18;
  std::vector<Double_t> utgtK18;
  std::vector<Double_t> vtgtK18;

  Int_t ntKurama;
  std::vector<Double_t> chisqrKurama;
  std::vector<Double_t> pKurama;
  std::vector<Double_t> qKurama;
  std::vector<Double_t> m2Kurama;
  std::vector<Double_t> thetaKurama;
  std::vector<Double_t> xtgtKurama;
  std::vector<Double_t> ytgtKurama;
  std::vector<Double_t> utgtKurama;
  std::vector<Double_t> vtgtKurama;
  std::vector<Double_t> pathwcKurama;
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
  std::vector<Double_t> xb;
  std::vector<Double_t> yb;
  std::vector<Double_t> ub;
  std::vector<Double_t> vb;
  std::vector<Double_t> xs;
  std::vector<Double_t> ys;
  std::vector<Double_t> us;
  std::vector<Double_t> vs;
  std::vector<Int_t> Kflag;
  std::vector<Int_t> Pflag;
  std::vector<Int_t> Heavyflag;

  //TPC RK
  std::vector<Int_t> isgoodTPCK18;
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

  std::vector<Int_t> isgoodTPCKurama;
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

  Int_t remain_nclTpc;
  std::vector<Double_t> remain_cluster_x;
  std::vector<Double_t> remain_cluster_y;
  std::vector<Double_t> remain_cluster_z;
  std::vector<Double_t> remain_cluster_de;
  std::vector<Int_t> remain_cluster_size;
  std::vector<Int_t> remain_cluster_layer;
  std::vector<Double_t> remain_cluster_mrow;
  std::vector<Double_t> remain_cluster_de_center;
  std::vector<Double_t> remain_cluster_x_center;
  std::vector<Double_t> remain_cluster_y_center;
  std::vector<Double_t> remain_cluster_z_center;
  std::vector<Int_t> remain_cluster_row_center;
  std::vector<Int_t> remain_cluster_houghflag;

  Int_t ntTpc;
  std::vector<Int_t> nhtrack;
  std::vector<Int_t> trackid;
  std::vector<Int_t> isXi;
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> isMultiloop;
  std::vector<Int_t> isInTarget;
  std::vector<Int_t> charge;
  std::vector<Int_t> pid;
  std::vector<Double_t> chisqr;
  std::vector<Double_t> pval;
  std::vector<Double_t> helix_cx;
  std::vector<Double_t> helix_cy;
  std::vector<Double_t> helix_z0;
  std::vector<Double_t> helix_r;
  std::vector<Double_t> helix_dz;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx;
  std::vector<Double_t> mom0;
  std::vector<Double_t> path;
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
  std::vector<Double_t> mom0_inverted;
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

  Int_t ncombiLreconfailed;
  std::vector<Int_t> pidLreconfailed;
  std::vector<Int_t> piidLreconfailed;
  std::vector<Double_t> LdecayvtxLreconfailed_x;
  std::vector<Double_t> LdecayvtxLreconfailed_y;
  std::vector<Double_t> LdecayvtxLreconfailed_z;
  std::vector<Double_t> LmassLreconfailed;
  std::vector<Double_t> LmomLreconfailed;
  std::vector<Double_t> LmomLreconfailed_x;
  std::vector<Double_t> LmomLreconfailed_y;
  std::vector<Double_t> LmomLreconfailed_z;
  std::vector<Double_t> pmomLreconfailed;
  std::vector<Double_t> pmomLreconfailed_x;
  std::vector<Double_t> pmomLreconfailed_y;
  std::vector<Double_t> pmomLreconfailed_z;
  std::vector<Double_t> pimomLreconfailed;
  std::vector<Double_t> pimomLreconfailed_x;
  std::vector<Double_t> pimomLreconfailed_y;
  std::vector<Double_t> pimomLreconfailed_z;
  std::vector<Double_t> ppidistLreconfailed;
  
  Int_t nEscapeKm;
  std::vector<Int_t> kmid;
  std::vector<Double_t> kmmom;
  std::vector<Double_t> kmmom_x;
  std::vector<Double_t> kmmom_y;
  std::vector<Double_t> kmmom_z;
  
  std::vector<Double_t> GFkmdecayvtx_x;
  std::vector<Double_t> GFkmdecayvtx_y;
  std::vector<Double_t> GFkmdecayvtx_z;
  std::vector<Double_t> GFkmmom;
  std::vector<Double_t> GFkmmom_x;
  std::vector<Double_t> GFkmmom_y;
  std::vector<Double_t> GFkmmom_z;
  std::vector<Double_t> GFkmtarget_dist;
  std::vector<Double_t> GFkmtargetvtx_x;
  std::vector<Double_t> GFkmtargetvtx_y;
  std::vector<Double_t> GFkmtargetvtx_z;
  std::vector<Double_t> GFkmtargetcenter_x;
  std::vector<Double_t> GFkmtargetcenter_y;
  std::vector<Double_t> GFkmtargetcenter_z;
  std::vector<Double_t> GFkmtargetcenter_dist;
  std::vector<Double_t> GFkmm2;  
  std::vector<Double_t> GFkmtracklen;
  std::vector<Double_t> GFkmtof;  

  Int_t ncombiPipair;
  std::vector<Int_t> pipidPipair;
  std::vector<Int_t> pimidPipair;
  std::vector<Double_t> pipmomPipair;
  std::vector<Double_t> pipmomPipair_x;
  std::vector<Double_t> pipmomPipair_y;
  std::vector<Double_t> pipmomPipair_z;
  std::vector<Double_t> pimmomPipair;
  std::vector<Double_t> pimmomPipair_x;
  std::vector<Double_t> pimmomPipair_y;
  std::vector<Double_t> pimmomPipair_z;
  std::vector<Double_t> momPipair;
  std::vector<Double_t> momPipair_x;
  std::vector<Double_t> momPipair_y;
  std::vector<Double_t> momPipair_z;
  std::vector<Double_t> reconLmassPipair;
  std::vector<Double_t> reconmassPipair;
  std::vector<Double_t> pipidistPipair;

  // Genfit Track's information
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

  std::vector<Int_t> GFinside;
  Int_t GFKuramaFromTgt;
  Int_t GFKuramaVtxOutTgt;  

  Int_t GFntTpc_target;
  Double_t GFprodvtx_x;
  Double_t GFprodvtx_y;
  Double_t GFprodvtx_z;

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

  //Multi-track production vertex
  Double_t GFprodvtx_x_ll;
  Double_t GFprodvtx_y_ll;
  Double_t GFprodvtx_z_ll;
  Double_t GFprodvtx_x_l1;
  Double_t GFprodvtx_y_l1;
  Double_t GFprodvtx_z_l1;
  Double_t GFprodvtx_x_l2;
  Double_t GFprodvtx_y_l2;
  Double_t GFprodvtx_z_l2;
  Double_t GFprodvtx_x_l;
  Double_t GFprodvtx_y_l;
  Double_t GFprodvtx_z_l;

  Bool_t emptyflag;
  Bool_t pimflag;
  Bool_t lpiflag;
  Bool_t lpflag;
  Bool_t lflag;
  Bool_t kuramalflag;

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

  Double_t GFllexcitation;
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
  Double_t GFltargetcenter_dist1;
  Double_t GFltargetcenter_x1;
  Double_t GFltargetcenter_y1;
  Double_t GFltargetcenter_z1;

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
  Double_t GFltargetcenter_dist2;
  Double_t GFltargetcenter_x2;
  Double_t GFltargetcenter_y2;
  Double_t GFltargetcenter_z2;

  Double_t GFlprodvtx_x2;
  Double_t GFlprodvtx_y2;
  Double_t GFlprodvtx_z2;
  Double_t GFlprodvtx_dist2;
  Double_t GFltracklen2;
  Double_t GFltof2;

  Double_t GFlmass_alter1;
  Double_t GFldecayvtx_x_alter1;
  Double_t GFldecayvtx_y_alter1;
  Double_t GFldecayvtx_z_alter1;
  Double_t GFlmom_alter1;
  Double_t GFlmom_x_alter1;
  Double_t GFlmom_y_alter1;
  Double_t GFlmom_z_alter1;
  Double_t GFppi_dist_alter1;
  Double_t GFltarget_dist_alter1;
  Double_t GFltargetvtx_x_alter1;
  Double_t GFltargetvtx_y_alter1;
  Double_t GFltargetvtx_z_alter1;
  Double_t GFltargetcenter_dist_alter1;
  Double_t GFltargetcenter_x_alter1;
  Double_t GFltargetcenter_y_alter1;
  Double_t GFltargetcenter_z_alter1;

  Double_t GFlmass_alter2;
  Double_t GFldecayvtx_x_alter2;
  Double_t GFldecayvtx_y_alter2;
  Double_t GFldecayvtx_z_alter2;
  Double_t GFlmom_alter2;
  Double_t GFlmom_x_alter2;
  Double_t GFlmom_y_alter2;
  Double_t GFlmom_z_alter2;
  Double_t GFppi_dist_alter2;
  Double_t GFltarget_dist_alter2;
  Double_t GFltargetvtx_x_alter2;
  Double_t GFltargetvtx_y_alter2;
  Double_t GFltargetvtx_z_alter2;
  Double_t GFltargetcenter_dist_alter2;
  Double_t GFltargetcenter_x_alter2;
  Double_t GFltargetcenter_y_alter2;
  Double_t GFltargetcenter_z_alter2;

  Double_t llvtx_x;
  Double_t llvtx_y;
  Double_t llvtx_z;
  Double_t lldist;
  Double_t GFllvtx_x;
  Double_t GFllvtx_y;
  Double_t GFllvtx_z;
  Double_t GFlldist;

  Double_t KFllexcitation;
  Double_t KFlpval1;
  std::vector<Double_t> KFlpull1;
  Double_t KFlchisqr1;
  Double_t KFlmom1;
  Double_t KFlmom_x1;
  Double_t KFlmom_y1;
  Double_t KFlmom_z1;
  Double_t KFlpval2;
  std::vector<Double_t> KFlpull2;
  Double_t KFlchisqr2;
  Double_t KFlmom2;
  Double_t KFlmom_x2;
  Double_t KFlmom_y2;
  Double_t KFlmom_z2;

  Double_t KFprodvtx_chisqr_ll;
  Double_t KFprodvtx_x_ll;
  Double_t KFprodvtx_y_ll;
  Double_t KFprodvtx_z_ll;
  Double_t KFprodvtx_x_l1;
  Double_t KFprodvtx_y_l1;
  Double_t KFprodvtx_z_l1;
  Double_t KFprodvtx_x_l2;
  Double_t KFprodvtx_y_l2;
  Double_t KFprodvtx_z_l2;
  Double_t KFprodvtx_x_l;
  Double_t KFprodvtx_y_l;
  Double_t KFprodvtx_z_l;

  Double_t KFllvtx_x;
  Double_t KFllvtx_y;
  Double_t KFllvtx_z;
  Double_t KFlldist;

  Double_t KFlprodvtx_x1;
  Double_t KFlprodvtx_y1;
  Double_t KFlprodvtx_z1;
  Double_t KFlprodvtx_dist1;
  Double_t KFltracklen1;
  Double_t KFltof1;

  Double_t KFlprodvtx_x2;
  Double_t KFlprodvtx_y2;
  Double_t KFlprodvtx_z2;
  Double_t KFlprodvtx_dist2;
  Double_t KFltracklen2;
  Double_t KFltof2;

  Double_t KFlprodvtx_x;
  Double_t KFlprodvtx_y;
  Double_t KFlprodvtx_z;
  Double_t KFlprodvtx_dist;
  Double_t KFltracklen;
  Double_t KFltof;

  Bool_t xiflag;
  Bool_t xipflag;
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

  Double_t lmass_vtx;
  Double_t ldecayvtx_x_vtx;
  Double_t ldecayvtx_y_vtx;
  Double_t ldecayvtx_z_vtx;
  Double_t lmom_vtx;
  Double_t lmom_x_vtx;
  Double_t lmom_y_vtx;
  Double_t lmom_z_vtx;
  Double_t ppi_dist_vtx;

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
  Double_t GFXiexcitation;

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
  Double_t GFlctau;

  //Multi-track production vertex
  Double_t GFprodvtx_x_kkxi;
  Double_t GFprodvtx_y_kkxi;
  Double_t GFprodvtx_z_kkxi;

  Bool_t lphiflag;
  Double_t phimass;
  Double_t phidecayvtx_x;
  Double_t phidecayvtx_y;
  Double_t phidecayvtx_z;
  Double_t phicosKK;
  Double_t phimom;
  Double_t phimom_x;
  Double_t phimom_y;
  Double_t phimom_z;
  Double_t kk_dist;
  Double_t GFphi_km_mass2;
  Double_t GFphi_km_invbeta;
  Double_t GFphi_kp_mass2;
  Double_t GFphi_kp_invbeta;

  Double_t GFphimass;
  Double_t GFphidecayvtx_x;
  Double_t GFphidecayvtx_y;
  Double_t GFphidecayvtx_z;
  Double_t GFphicosKK;
  Double_t GFphimom;
  Double_t GFphimom_x;
  Double_t GFphimom_y;
  Double_t GFphimom_z;
  Double_t GFkk_dist;
  Double_t GFphiprodvtx_dist;

  Double_t GFphimass_wKurama;
  Double_t GFphimom_wKurama;
  Double_t GFphimom_x_wKurama;
  Double_t GFphimom_y_wKurama;
  Double_t GFphimom_z_wKurama;
  std::vector<Int_t> phidecays_id;
  std::vector<Double_t> phidecays_mom;
  std::vector<Double_t> phidecays_mom_x;
  std::vector<Double_t> phidecays_mom_y;
  std::vector<Double_t> phidecays_mom_z;
  std::vector<Double_t> GFphidecays_mom;
  std::vector<Double_t> GFphidecays_mom_x;
  std::vector<Double_t> GFphidecays_mom_y;
  std::vector<Double_t> GFphidecays_mom_z;

  std::vector<Int_t>     GFlldecays_pdgcode;
  std::vector<Int_t>     GFlldecays_nhtrack;
  std::vector<Double_t>  GFlldecays_charge;
  std::vector<Double_t>  GFlldecays_chisqr;
  std::vector<Double_t>  GFlldecays_pval;
  std::vector<Int_t>     GFlldecays_htofid;
  std::vector<Double_t>  GFlldecays_tracklen;
  std::vector<Double_t>  GFlldecays_tof;
  std::vector<Double_t>  GFlldecays_mass2;
  std::vector<Double_t>  GFlldecays_invbeta;
  std::vector<Double_t>  GFlldecays_mom;
  std::vector<Double_t>  GFlldecays_mom_x;
  std::vector<Double_t>  GFlldecays_mom_y;
  std::vector<Double_t>  GFlldecays_mom_z;
  std::vector<Double_t>  GFlldecays_CMmom;
  std::vector<Double_t>  GFlldecays_CMmom_x;
  std::vector<Double_t>  GFlldecays_CMmom_y;
  std::vector<Double_t>  GFlldecays_CMmom_z;
  std::vector<Double_t>  GFlldecays_momloss;
  std::vector<Double_t>  GFlldecays_eloss;

  std::vector<Int_t> lldecays_id;
  std::vector<Double_t> lldecays_mom;
  std::vector<Double_t> lldecays_mom_x;
  std::vector<Double_t> lldecays_mom_y;
  std::vector<Double_t> lldecays_mom_z;
  std::vector<Double_t> lldecays_CMmom;
  std::vector<Double_t> lldecays_CMmom_x;
  std::vector<Double_t> lldecays_CMmom_y;
  std::vector<Double_t> lldecays_CMmom_z;

  std::vector<Int_t> GFxidecays_pdgcode;
  std::vector<Int_t> GFxidecays_nhtrack;
  std::vector<Double_t> GFxidecays_charge;
  std::vector<Double_t> GFxidecays_chisqr;
  std::vector<Double_t> GFxidecays_tracktof;
  std::vector<Double_t> GFxidecays_pval;
  std::vector<Int_t> GFxidecays_htofid;
  std::vector<Double_t> GFxidecays_tracklen;
  std::vector<Double_t> GFxidecays_tof;
  std::vector<Double_t> GFxidecays_mass2;
  std::vector<Double_t> GFxidecays_invbeta;
  std::vector<Double_t> GFxidecays_mom;
  std::vector<Double_t> GFxidecays_mom_x;
  std::vector<Double_t> GFxidecays_mom_y;
  std::vector<Double_t> GFxidecays_mom_z;
  std::vector<Double_t> GFxidecays_CMmom;
  std::vector<Double_t> GFxidecays_CMmom_x;
  std::vector<Double_t> GFxidecays_CMmom_y;
  std::vector<Double_t> GFxidecays_CMmom_z;
  std::vector<Double_t> GFxidecays_momloss;
  std::vector<Double_t> GFxidecays_eloss;

  std::vector<Int_t> xidecays_id;
  std::vector<Double_t> xidecays_mom;
  std::vector<Double_t> xidecays_mom_x;
  std::vector<Double_t> xidecays_mom_y;
  std::vector<Double_t> xidecays_mom_z;
  std::vector<Double_t> xidecays_CMmom;
  std::vector<Double_t> xidecays_CMmom_x;
  std::vector<Double_t> xidecays_CMmom_y;
  std::vector<Double_t> xidecays_CMmom_z;

  std::vector<Int_t>    GFdecays_pdgcode; 
  std::vector<Int_t>    GFdecays_nhtrack; 
  std::vector<Double_t> GFdecays_charge;  
  std::vector<Double_t> GFdecays_chisqr;  
  std::vector<Double_t> GFdecays_pval;	  
  std::vector<Int_t>    GFdecays_htofid;  
  std::vector<Double_t> GFdecays_tracklen;
  std::vector<Double_t> GFdecays_tof;	  
  std::vector<Double_t> GFdecays_mass2;	  
  std::vector<Double_t> GFdecays_invbeta; 
  std::vector<Double_t> GFdecays_mom;	  
  std::vector<Double_t> GFdecays_mom_x;	  
  std::vector<Double_t> GFdecays_mom_y;	  
  std::vector<Double_t> GFdecays_mom_z;	  
  std::vector<Double_t> GFdecays_CMmom;	  
  std::vector<Double_t> GFdecays_CMmom_x; 
  std::vector<Double_t> GFdecays_CMmom_y; 
  std::vector<Double_t> GFdecays_CMmom_z; 
  std::vector<Double_t> GFdecays_momloss; 
  std::vector<Double_t> GFdecays_eloss;   

  std::vector<Int_t>    decays_id;     
  std::vector<Double_t> decays_mom;    
  std::vector<Double_t> decays_mom_x;  
  std::vector<Double_t> decays_mom_y;  
  std::vector<Double_t> decays_mom_z;  
  std::vector<Double_t> decays_CMmom;  
  std::vector<Double_t> decays_CMmom_x;
  std::vector<Double_t> decays_CMmom_y;
  std::vector<Double_t> decays_CMmom_z;

  Bool_t pipiflag;

  Int_t accident_multi;
  std::vector<Int_t> accident_id;

  Int_t xiresidual_multi;
  Int_t xipim_multi;
  Int_t xipip_multi;
  Int_t xiem_multi;
  Int_t xiep_multi;
  Int_t xip_multi;
  Int_t xippip_multi;

  std::vector<Int_t> xiresidual_id;
  std::vector<Double_t> xiresidual_dist2tgt;
  std::vector<Double_t> xiresidual_KFdist2prodvtx;
  std::vector<Double_t> xiresidual_GFdist2prodvtx;
  std::vector<Double_t> xiresidual_mass2;
  std::vector<Double_t> xiresidual_invbeta;
  std::vector<Double_t> xiresidual_mom;
  std::vector<Double_t> xiresidual_mom_x;
  std::vector<Double_t> xiresidual_mom_y;
  std::vector<Double_t> xiresidual_mom_z;
  std::vector<Double_t> xiresidual_charge;

  Int_t llresidual_multi;
  Int_t llpim_multi;
  Int_t llpip_multi;
  Int_t llem_multi;
  Int_t llep_multi;
  Int_t llp_multi;
  Int_t llppip_multi;

  std::vector<Int_t> llresidual_id;
  std::vector<Double_t> llresidual_dist2tgt;
  std::vector<Double_t> llresidual_GFdist2prodvtx;
  std::vector<Double_t> llresidual_KFdist2prodvtx;
  std::vector<Double_t> llresidual_mass2;
  std::vector<Double_t> llresidual_invbeta;
  std::vector<Double_t> llresidual_mom;
  std::vector<Double_t> llresidual_mom_x;
  std::vector<Double_t> llresidual_mom_y;
  std::vector<Double_t> llresidual_mom_z;
  std::vector<Double_t> llresidual_charge;

  Int_t residual_multi;
  Int_t pim_multi;
  Int_t pip_multi;
  Int_t em_multi;
  Int_t ep_multi;
  Int_t p_multi;
  Int_t ppip_multi;
  std::vector<Int_t> residual_id;
  std::vector<Double_t> residual_dist2tgt;
  std::vector<Double_t> residual_GFdist2prodvtx;
  std::vector<Double_t> residual_KFdist2prodvtx;
  std::vector<Double_t> residual_mass2;
  std::vector<Double_t> residual_invbeta;
  std::vector<Double_t> residual_mom;
  std::vector<Double_t> residual_mom_x;
  std::vector<Double_t> residual_mom_y;
  std::vector<Double_t> residual_mom_z;
  std::vector<Double_t> residual_charge;

  //Kinematic fitting
  Double_t KFlmom0;
  Double_t KFlmom_x0;
  Double_t KFlmom_y0;
  Double_t KFlmom_z0;
  Double_t KFlmom;
  Double_t KFlmom_x;
  Double_t KFlmom_y;
  Double_t KFlmom_z;
  Double_t KFlchisqr;
  Double_t KFlpval;
  std::vector<std::vector<Double_t>> KFlCovMatrix;
  Double_t KFlpi_dist;
  Double_t KFximom;
  Double_t KFximom_x;
  Double_t KFximom_y;
  Double_t KFximom_z;
  Double_t KFxichisqr;
  Double_t KFxipval;
  std::vector<std::vector<Double_t>> KFxiCovMatrix;
  Double_t KFximass;
  Double_t KFxidecayvtx_x;
  Double_t KFxidecayvtx_y;
  Double_t KFxidecayvtx_z;
  std::vector<Double_t> KFlpull;
  std::vector<Double_t> KFxipull;

  //Multi-track vertex
  Double_t KFprodvtx_chisqr_kkxi;
  Double_t KFprodvtx_x_kkxi;
  Double_t KFprodvtx_y_kkxi;
  Double_t KFprodvtx_z_kkxi;
  Double_t KFprodvtx_x_kpxi;
  Double_t KFprodvtx_y_kpxi;
  Double_t KFprodvtx_z_kpxi;

  Double_t KFxiprodvtx_x;
  Double_t KFxiprodvtx_y;
  Double_t KFxiprodvtx_z;
  Double_t KFxiprodmom;
  Double_t KFxiprodmom_x;
  Double_t KFxiprodmom_y;
  Double_t KFxiprodmom_z;
  Double_t KFxiprodvtx_dist;
  Double_t KFxitracklen;
  Double_t KFxitof;
  Double_t KFximomloss;
  Double_t KFXiexcitation;

  Double_t KFxi_kkvtx_x;
  Double_t KFxi_kkvtx_y;
  Double_t KFxi_kkvtx_z;
  Double_t KFxi_kkvtx_mom;
  Double_t KFxi_kkvtx_mom_x;
  Double_t KFxi_kkvtx_mom_y;
  Double_t KFxi_kkvtx_mom_z;
  Double_t KFxi_kkvtx_dist;

  Double_t KFxi_kpxiprodvtx_x;
  Double_t KFxi_kpxiprodvtx_y;
  Double_t KFxi_kpxiprodvtx_z;
  Double_t KFxi_kpxiprodmom;
  Double_t KFxi_kpxiprodmom_x;
  Double_t KFxi_kpxiprodmom_y;
  Double_t KFxi_kpxiprodmom_z;
  Double_t KFxi_kpxiprodvtx_dist;

  Double_t KFxitargetvtx_x;
  Double_t KFxitargetvtx_y;
  Double_t KFxitargetvtx_z;
  Double_t KFxitargetmom;
  Double_t KFxitargetmom_x;
  Double_t KFxitargetmom_y;
  Double_t KFxitargetmom_z;
  Double_t KFxitarget_dist;

  Double_t KFxitargetcenter_x;
  Double_t KFxitargetcenter_y;
  Double_t KFxitargetcenter_z;
  Double_t KFxitargetcentermom;
  Double_t KFxitargetcentermom_x;
  Double_t KFxitargetcentermom_y;
  Double_t KFxitargetcentermom_z;
  Double_t KFxitargetcenter_dist;

  std::vector<Double_t> KFlldecays_mom;
  std::vector<Double_t> KFlldecays_mom_x;
  std::vector<Double_t> KFlldecays_mom_y;
  std::vector<Double_t> KFlldecays_mom_z;
  std::vector<Double_t> KFlldecays_CMmom;
  std::vector<Double_t> KFlldecays_CMmom_x;
  std::vector<Double_t> KFlldecays_CMmom_y;
  std::vector<Double_t> KFlldecays_CMmom_z;

  std::vector<Double_t> KFxidecays_mom;
  std::vector<Double_t> KFxidecays_mom_x;
  std::vector<Double_t> KFxidecays_mom_y;
  std::vector<Double_t> KFxidecays_mom_z;
  std::vector<Double_t> KFxidecays_CMmom;
  std::vector<Double_t> KFxidecays_CMmom_x;
  std::vector<Double_t> KFxidecays_CMmom_y;
  std::vector<Double_t> KFxidecays_CMmom_z;

  std::vector<Double_t> KFdecays_mom;
  std::vector<Double_t> KFdecays_mom_x;
  std::vector<Double_t> KFdecays_mom_y;
  std::vector<Double_t> KFdecays_mom_z;
  std::vector<Double_t> KFdecays_CMmom;
  std::vector<Double_t> KFdecays_CMmom_x;
  std::vector<Double_t> KFdecays_CMmom_y;
  std::vector<Double_t> KFdecays_CMmom_z;

  //For gamma searching
  Int_t llg_multi;
  std::vector<Int_t> llepidgamma;
  std::vector<Int_t> llemidgamma;
  std::vector<Double_t> llepmomgamma;
  std::vector<Double_t> llepmomgamma_x;
  std::vector<Double_t> llepmomgamma_y;
  std::vector<Double_t> llepmomgamma_z;
  std::vector<Double_t> llemmomgamma;
  std::vector<Double_t> llemmomgamma_x;
  std::vector<Double_t> llemmomgamma_y;
  std::vector<Double_t> llemmomgamma_z;
  std::vector<Double_t> llmomgamma;
  std::vector<Double_t> llmomgamma_x;
  std::vector<Double_t> llmomgamma_y;
  std::vector<Double_t> llmomgamma_z;
  std::vector<Double_t> llepidistgamma;
  std::vector<Double_t> llvtxgamma_x;
  std::vector<Double_t> llvtxgamma_y;
  std::vector<Double_t> llvtxgamma_z;

  Int_t xig_multi;
  std::vector<Int_t> xiepidgamma;
  std::vector<Int_t> xiemidgamma;
  std::vector<Double_t> xiepmomgamma;
  std::vector<Double_t> xiepmomgamma_x;
  std::vector<Double_t> xiepmomgamma_y;
  std::vector<Double_t> xiepmomgamma_z;
  std::vector<Double_t> xiemmomgamma;
  std::vector<Double_t> xiemmomgamma_x;
  std::vector<Double_t> xiemmomgamma_y;
  std::vector<Double_t> xiemmomgamma_z;
  std::vector<Double_t> ximomgamma;
  std::vector<Double_t> ximomgamma_x;
  std::vector<Double_t> ximomgamma_y;
  std::vector<Double_t> ximomgamma_z;
  std::vector<Double_t> xiepidistgamma;
  std::vector<Double_t> xivtxgamma_x;
  std::vector<Double_t> xivtxgamma_y;
  std::vector<Double_t> xivtxgamma_z;

  Int_t g_multi;
  std::vector<Int_t> epidgamma;
  std::vector<Int_t> emidgamma;
  std::vector<Double_t> epmomgamma;
  std::vector<Double_t> epmomgamma_x;
  std::vector<Double_t> epmomgamma_y;
  std::vector<Double_t> epmomgamma_z;
  std::vector<Double_t> emmomgamma;
  std::vector<Double_t> emmomgamma_x;
  std::vector<Double_t> emmomgamma_y;
  std::vector<Double_t> emmomgamma_z;
  std::vector<Double_t> momgamma;
  std::vector<Double_t> momgamma_x;
  std::vector<Double_t> momgamma_y;
  std::vector<Double_t> momgamma_z;
  std::vector<Double_t> epidistgamma;
  std::vector<Double_t> vtxgamma_x;
  std::vector<Double_t> vtxgamma_y;
  std::vector<Double_t> vtxgamma_z;

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

    ntKurama = 0;
    chisqrKurama.clear();
    pKurama.clear();
    qKurama.clear();
    m2Kurama.clear();
    thetaKurama.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();
    pathwcKurama.clear();
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
    xb.clear();
    yb.clear();
    ub.clear();
    vb.clear();
    xs.clear();
    ys.clear();
    us.clear();
    vs.clear();
    Kflag.clear();
    Pflag.clear();
    Heavyflag.clear();    

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

    isgoodTPCKurama.clear();
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

    BE.clear();
    BETPC.clear();
    BE_LL.clear();
    BETPC_LL.clear();
    km_mom_x.clear();
    km_mom_y.clear();
    km_mom_z.clear();
    kp_mom_x.clear();
    kp_mom_y.clear();
    kp_mom_z.clear();

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

    remain_nclTpc = 0;
    remain_cluster_x.clear();
    remain_cluster_y.clear();
    remain_cluster_z.clear();
    remain_cluster_de.clear();
    remain_cluster_size.clear();
    remain_cluster_layer.clear();
    remain_cluster_mrow.clear();
    remain_cluster_de_center.clear();
    remain_cluster_x_center.clear();
    remain_cluster_y_center.clear();
    remain_cluster_z_center.clear();
    remain_cluster_row_center.clear();
    remain_cluster_houghflag.clear();

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isXi.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
    isMultiloop.clear();
    isInTarget.clear();
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

    nvtxTpcClustered = 0;
    clusteredVtx_x.clear();
    clusteredVtx_y.clear();
    clusteredVtx_z.clear();
    clusteredVtxid.clear();

    ncombiLreconfailed = 0;
    pidLreconfailed.clear();
    piidLreconfailed.clear();
    LdecayvtxLreconfailed_x.clear();
    LdecayvtxLreconfailed_y.clear();
    LdecayvtxLreconfailed_z.clear();
    LmassLreconfailed.clear();
    LmomLreconfailed.clear();
    LmomLreconfailed_x.clear();
    LmomLreconfailed_y.clear();
    LmomLreconfailed_z.clear();
    pmomLreconfailed.clear();
    pmomLreconfailed_x.clear();
    pmomLreconfailed_y.clear();
    pmomLreconfailed_z.clear();
    pimomLreconfailed.clear();
    pimomLreconfailed_x.clear();
    pimomLreconfailed_y.clear();
    pimomLreconfailed_z.clear();
    ppidistLreconfailed.clear();

    nEscapeKm = 0;
    kmid.clear();
    kmmom.clear();
    kmmom_x.clear();
    kmmom_y.clear();
    kmmom_z.clear();

    GFkmdecayvtx_x.clear();
    GFkmdecayvtx_y.clear();
    GFkmdecayvtx_z.clear();
    GFkmmom.clear();
    GFkmmom_x.clear();
    GFkmmom_y.clear();
    GFkmmom_z.clear();
    //GFkmpi_dist;
    GFkmtarget_dist.clear();
    GFkmtargetvtx_x.clear();
    GFkmtargetvtx_y.clear();
    GFkmtargetvtx_z.clear();
    GFkmtargetcenter_x.clear();
    GFkmtargetcenter_y.clear();
    GFkmtargetcenter_z.clear();
    GFkmtargetcenter_dist.clear();

    GFkmm2.clear();
    GFkmtracklen.clear();		       
    GFkmtof.clear();

    ncombiPipair = 0;
    pipidPipair.clear();
    pimidPipair.clear();
    pipmomPipair.clear();
    pipmomPipair_x.clear();
    pipmomPipair_y.clear();
    pipmomPipair_z.clear();
    pimmomPipair.clear();
    pimmomPipair_x.clear();
    pimmomPipair_y.clear();
    pimmomPipair_z.clear();
    momPipair.clear();
    momPipair_x.clear();
    momPipair_y.clear();
    momPipair_z.clear();
    reconLmassPipair.clear();
    reconmassPipair.clear();
    pipidistPipair.clear();

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

    GFinside.clear();
    GFKuramaFromTgt = 0;    
    GFKuramaVtxOutTgt = 0;  
    
    GFntTpc_target = 0;
    GFprodvtx_x = qnan;
    GFprodvtx_y = qnan;
    GFprodvtx_z = qnan;

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

    GFprodvtx_x_ll = qnan;
    GFprodvtx_y_ll = qnan;
    GFprodvtx_z_ll = qnan;
    GFprodvtx_x_l1 = qnan;
    GFprodvtx_y_l1 = qnan;
    GFprodvtx_z_l1 = qnan;
    GFprodvtx_x_l2 = qnan;
    GFprodvtx_y_l2 = qnan;
    GFprodvtx_z_l2 = qnan;
    GFprodvtx_x_l = qnan;
    GFprodvtx_y_l = qnan;
    GFprodvtx_z_l = qnan;

    GFdecays_pdgcode.clear();
    GFdecays_nhtrack.clear();
    GFdecays_charge.clear();
    GFdecays_chisqr.clear();
    GFdecays_pval.clear();
    GFdecays_htofid.clear();
    GFdecays_tracklen.clear();
    GFdecays_tof.clear();
    GFdecays_mass2.clear();
    GFdecays_invbeta.clear();
    GFdecays_mom.clear();
    GFdecays_mom_x.clear();
    GFdecays_mom_y.clear();
    GFdecays_mom_z.clear();
    GFdecays_CMmom.clear();
    GFdecays_CMmom_x.clear();
    GFdecays_CMmom_y.clear();
    GFdecays_CMmom_z.clear();
    GFdecays_momloss.clear();
    GFdecays_eloss.clear();

    decays_id.clear();
    decays_mom.clear();
    decays_mom_x.clear();
    decays_mom_y.clear();
    decays_mom_z.clear();
    decays_CMmom.clear();
    decays_CMmom_x.clear();
    decays_CMmom_y.clear();
    decays_CMmom_z.clear();

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

    GFlmass = qnan;
    GFldecayvtx_x = qnan;
    GFldecayvtx_y = qnan;
    GFldecayvtx_z = qnan;
    GFlmom = qnan;
    GFlmom_x = qnan;
    GFlmom_y = qnan;
    GFlmom_z = qnan;
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
    GFlctau = qnan;    

    GFprodvtx_x_kkxi = qnan;
    GFprodvtx_y_kkxi = qnan;
    GFprodvtx_z_kkxi = qnan;

    emptyflag = false;
    pimflag = false;
    lpiflag = false;
    lpflag = false;    
    lflag = false;
    kuramalflag = false;    
    llflag = false;
    pipiflag = false;

    accident_multi = 0;
    accident_id.clear();

    residual_multi = 0;
    pim_multi = 0;
    pip_multi = 0;
    em_multi = 0;
    ep_multi = 0;
    p_multi = 0;
    ppip_multi = 0;

    residual_id.clear();
    residual_dist2tgt.clear();
    residual_GFdist2prodvtx.clear();
    residual_KFdist2prodvtx.clear();
    residual_mass2.clear();
    residual_invbeta.clear();
    residual_mom.clear();
    residual_mom_x.clear();
    residual_mom_y.clear();
    residual_mom_z.clear();
    residual_charge.clear();

    g_multi = 0;
    epidgamma.clear();
    emidgamma.clear();
    epmomgamma.clear();
    epmomgamma_x.clear();
    epmomgamma_y.clear();
    epmomgamma_z.clear();
    emmomgamma.clear();
    emmomgamma_x.clear();
    emmomgamma_y.clear();
    emmomgamma_z.clear();
    momgamma.clear();
    momgamma_x.clear();
    momgamma_y.clear();
    momgamma_z.clear();
    epidistgamma.clear();
    vtxgamma_x.clear();
    vtxgamma_y.clear();
    vtxgamma_z.clear();      
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* status;
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

  TTreeReaderValue<Int_t>* nhHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;

  TTreeReaderValue<Int_t>* ntK18;
  TTreeReaderValue<std::vector<Double_t>>* pK18;
  TTreeReaderValue<std::vector<Double_t>>* chisqrK18;
  TTreeReaderValue<std::vector<Double_t>>* xtgtK18;
  TTreeReaderValue<std::vector<Double_t>>* ytgtK18;
  TTreeReaderValue<std::vector<Double_t>>* utgtK18;
  TTreeReaderValue<std::vector<Double_t>>* vtgtK18;

  TTreeReaderValue<Int_t>* ntKurama;
  TTreeReaderValue<std::vector<Double_t>>* chisqrKurama;
  TTreeReaderValue<std::vector<Double_t>>* pKurama;
  TTreeReaderValue<std::vector<Double_t>>* qKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2Kurama;
  TTreeReaderValue<std::vector<Double_t>>* thetaKurama;
  TTreeReaderValue<std::vector<Double_t>>* xtgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* ytgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* utgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* vtgtKurama;
  TTreeReaderValue<std::vector<Double_t>>* pathwcKurama;
  TTreeReaderValue<std::vector<Double_t>>* xin;
  TTreeReaderValue<std::vector<Double_t>>* yin;
  TTreeReaderValue<std::vector<Double_t>>* zin;
  TTreeReaderValue<std::vector<Double_t>>* pxin;
  TTreeReaderValue<std::vector<Double_t>>* pyin;
  TTreeReaderValue<std::vector<Double_t>>* pzin;
  TTreeReaderValue<std::vector<Double_t>>* xout;
  TTreeReaderValue<std::vector<Double_t>>* yout;
  TTreeReaderValue<std::vector<Double_t>>* zout;
  TTreeReaderValue<std::vector<Double_t>>* pxout;
  TTreeReaderValue<std::vector<Double_t>>* pyout;
  TTreeReaderValue<std::vector<Double_t>>* pzout;

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
  TTreeReaderValue<std::vector<Int_t>>* Kflag;
  TTreeReaderValue<std::vector<Int_t>>* Pflag;
  TTreeReaderValue<std::vector<Int_t>>* Heavyflag;

  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* chisqrTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* pTPCK18;
  TTreeReaderValue<std::vector<Double_t>>* qTPCK18;
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

  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCKurama;
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
  TTreeReaderValue<std::vector<Double_t>>* pOrgTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCalcTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* costCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* xbTPC;
  TTreeReaderValue<std::vector<Double_t>>* ybTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;
  TTreeReaderValue<std::vector<Double_t>>* xsTPC;
  TTreeReaderValue<std::vector<Double_t>>* ysTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;

  TTreeReaderValue<std::vector<Double_t>>* BE;
  TTreeReaderValue<std::vector<Double_t>>* BETPC;
  TTreeReaderValue<std::vector<Double_t>>* BE_LL;
  TTreeReaderValue<std::vector<Double_t>>* BETPC_LL;
  TTreeReaderValue<std::vector<Double_t>>* km_mom_x;
  TTreeReaderValue<std::vector<Double_t>>* km_mom_y;
  TTreeReaderValue<std::vector<Double_t>>* km_mom_z;
  TTreeReaderValue<std::vector<Double_t>>* kp_mom_x;
  TTreeReaderValue<std::vector<Double_t>>* kp_mom_y;
  TTreeReaderValue<std::vector<Double_t>>* kp_mom_z;

  TTreeReaderValue<Int_t>* nclTpc;
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

  TTreeReaderValue<Int_t>* remain_nclTpc;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_x;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_y;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_z;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_de;
  TTreeReaderValue<std::vector<Int_t>>* remain_cluster_size;
  TTreeReaderValue<std::vector<Int_t>>* remain_cluster_layer;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_mrow;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_de_center;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_x_center;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_y_center;
  TTreeReaderValue<std::vector<Double_t>>* remain_cluster_z_center;
  TTreeReaderValue<std::vector<Int_t>>* remain_cluster_row_center;
  TTreeReaderValue<std::vector<Int_t>>* remain_cluster_houghflag;

  TTreeReaderValue<Int_t>* ntTpc;
  TTreeReaderValue<std::vector<Int_t>>* nhtrack;
  TTreeReaderValue<std::vector<Int_t>>* trackid;
  TTreeReaderValue<std::vector<Int_t>>* isXi;
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
  TTreeReaderValue<std::vector<Int_t>>* isMultiloop;
  TTreeReaderValue<std::vector<Int_t>>* isInTarget;
  TTreeReaderValue<std::vector<Int_t>>* charge;
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
  TTreeReaderValue<std::vector<Double_t>>* pval;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx;
  TTreeReaderValue<std::vector<Double_t>>* mom0;
  TTreeReaderValue<std::vector<Double_t>>* path;
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
  TTreeReaderValue<std::vector<Double_t>>* chisqr_inverted;
  TTreeReaderValue<std::vector<Double_t>>* pval_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_r_inverted;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz_inverted;
  TTreeReaderValue<std::vector<Double_t>>* mom0_inverted;
  TTreeReaderValue<std::vector<Int_t>>* pid_inverted;

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

  TTreeReaderValue<Int_t>* ncombiLreconfailed;
  TTreeReaderValue<std::vector<Int_t>>* pidLreconfailed;
  TTreeReaderValue<std::vector<Int_t>>* piidLreconfailed;
  TTreeReaderValue<std::vector<Double_t>>* LdecayvtxLreconfailed_x;
  TTreeReaderValue<std::vector<Double_t>>* LdecayvtxLreconfailed_y;
  TTreeReaderValue<std::vector<Double_t>>* LdecayvtxLreconfailed_z;
  TTreeReaderValue<std::vector<Double_t>>* LmassLreconfailed;
  TTreeReaderValue<std::vector<Double_t>>* LmomLreconfailed;
  TTreeReaderValue<std::vector<Double_t>>* LmomLreconfailed_x;
  TTreeReaderValue<std::vector<Double_t>>* LmomLreconfailed_y;
  TTreeReaderValue<std::vector<Double_t>>* LmomLreconfailed_z;
  TTreeReaderValue<std::vector<Double_t>>* pmomLreconfailed;
  TTreeReaderValue<std::vector<Double_t>>* pmomLreconfailed_x;
  TTreeReaderValue<std::vector<Double_t>>* pmomLreconfailed_y;
  TTreeReaderValue<std::vector<Double_t>>* pmomLreconfailed_z;
  TTreeReaderValue<std::vector<Double_t>>* pimomLreconfailed;
  TTreeReaderValue<std::vector<Double_t>>* pimomLreconfailed_x;
  TTreeReaderValue<std::vector<Double_t>>* pimomLreconfailed_y;
  TTreeReaderValue<std::vector<Double_t>>* pimomLreconfailed_z;
  TTreeReaderValue<std::vector<Double_t>>* ppidistLreconfailed;
  
  TTreeReaderValue<Int_t>* nEscapeKm;
  TTreeReaderValue<std::vector<Int_t>>* kmid;
  TTreeReaderValue<std::vector<Double_t>>* kmmom;
  TTreeReaderValue<std::vector<Double_t>>* kmmom_x;
  TTreeReaderValue<std::vector<Double_t>>* kmmom_y;
  TTreeReaderValue<std::vector<Double_t>>* kmmom_z;
  
  TTreeReaderValue<std::vector<Double_t>>* GFkmdecayvtx_x;
  TTreeReaderValue<std::vector<Double_t>>* GFkmdecayvtx_y;
  TTreeReaderValue<std::vector<Double_t>>* GFkmdecayvtx_z;
  TTreeReaderValue<std::vector<Double_t>>* GFkmmom;
  TTreeReaderValue<std::vector<Double_t>>* GFkmmom_x;
  TTreeReaderValue<std::vector<Double_t>>* GFkmmom_y;
  TTreeReaderValue<std::vector<Double_t>>* GFkmmom_z;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtarget_dist;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtargetvtx_x;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtargetvtx_y;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtargetvtx_z;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtargetcenter_x;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtargetcenter_y;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtargetcenter_z;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtargetcenter_dist;
  TTreeReaderValue<std::vector<Double_t>>* GFkmm2;  
  TTreeReaderValue<std::vector<Double_t>>* GFkmtracklen;
  TTreeReaderValue<std::vector<Double_t>>* GFkmtof;  

  TTreeReaderValue<Int_t>* ncombiPipair;
  TTreeReaderValue<std::vector<Int_t>>* pipidPipair;
  TTreeReaderValue<std::vector<Int_t>>* pimidPipair;
  TTreeReaderValue<std::vector<Double_t>>* pipmomPipair;
  TTreeReaderValue<std::vector<Double_t>>* pipmomPipair_x;
  TTreeReaderValue<std::vector<Double_t>>* pipmomPipair_y;
  TTreeReaderValue<std::vector<Double_t>>* pipmomPipair_z;
  TTreeReaderValue<std::vector<Double_t>>* pimmomPipair;
  TTreeReaderValue<std::vector<Double_t>>* pimmomPipair_x;
  TTreeReaderValue<std::vector<Double_t>>* pimmomPipair_y;
  TTreeReaderValue<std::vector<Double_t>>* pimmomPipair_z;
  TTreeReaderValue<std::vector<Double_t>>* momPipair;
  TTreeReaderValue<std::vector<Double_t>>* momPipair_x;
  TTreeReaderValue<std::vector<Double_t>>* momPipair_y;
  TTreeReaderValue<std::vector<Double_t>>* momPipair_z;
  TTreeReaderValue<std::vector<Double_t>>* reconLmassPipair;
  TTreeReaderValue<std::vector<Double_t>>* reconmassPipair;
  TTreeReaderValue<std::vector<Double_t>>* pipidistPipair;
  
  TTreeReaderValue<std::vector<Int_t>>* GFfitstatus;
  TTreeReaderValue<std::vector<Int_t>>* GFpdgcode;
  TTreeReaderValue<std::vector<Int_t>>* GFnhtrack;
  TTreeReaderValue<std::vector<Double_t>>* GFcharge;
  TTreeReaderValue<std::vector<Double_t>>* GFchisqr;
  TTreeReaderValue<std::vector<Double_t>>* GFtof;
  TTreeReaderValue<std::vector<Double_t>>* GFpval;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFlayer;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFpos_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFpos_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFpos_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFmom;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFmom_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFmom_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFmom_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFresidual_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFresidual_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFresidual_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFresidual_p;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFresidual_px;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFresidual_py;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* GFresidual_pz;

  TTreeReaderValue<std::vector<Int_t>>* GFinside;
  TTreeReaderValue<Int_t>* GFKuramaFromTgt;
  TTreeReaderValue<Int_t>* GFKuramaVtxOutTgt;  

  TTreeReaderValue<Int_t>* GFntTpc_target;
  TTreeReaderValue<Double_t>* GFprodvtx_x;
  TTreeReaderValue<Double_t>* GFprodvtx_y;
  TTreeReaderValue<Double_t>* GFprodvtx_z;

  TTreeReaderValue<std::vector<Double_t>>* GFtracklen;
  TTreeReaderValue<std::vector<Double_t>>* GFtrack2vtxdist;
  TTreeReaderValue<std::vector<Double_t>>* GFcalctof;
  TTreeReaderValue<std::vector<Double_t>>* GFsegHtof;
  TTreeReaderValue<std::vector<Double_t>>* GFtofHtof;
  TTreeReaderValue<std::vector<Double_t>>* GFtdiffHtof;
  TTreeReaderValue<std::vector<Double_t>>* GFposHtof;
  TTreeReaderValue<std::vector<Double_t>>* GFposx;
  TTreeReaderValue<std::vector<Double_t>>* GFposy;
  TTreeReaderValue<std::vector<Double_t>>* GFposz;
  TTreeReaderValue<std::vector<Double_t>>* GFinvbeta;
  TTreeReaderValue<std::vector<Double_t>>* GFm2;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_tritonHtof;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_deutronHtof;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_protonHtof;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_kaonHtof;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_pionHtof;
  TTreeReaderValue<std::vector<Double_t>>* nsigma_electronHtof;

  TTreeReaderValue<Double_t>* GFprodvtx_x_ll;
  TTreeReaderValue<Double_t>* GFprodvtx_y_ll;
  TTreeReaderValue<Double_t>* GFprodvtx_z_ll;
  TTreeReaderValue<Double_t>* GFprodvtx_x_l1;
  TTreeReaderValue<Double_t>* GFprodvtx_y_l1;
  TTreeReaderValue<Double_t>* GFprodvtx_z_l1;
  TTreeReaderValue<Double_t>* GFprodvtx_x_l2;
  TTreeReaderValue<Double_t>* GFprodvtx_y_l2;
  TTreeReaderValue<Double_t>* GFprodvtx_z_l2;
  TTreeReaderValue<Double_t>* GFprodvtx_x_l;
  TTreeReaderValue<Double_t>* GFprodvtx_y_l;
  TTreeReaderValue<Double_t>* GFprodvtx_z_l;

  TTreeReaderValue<Bool_t>* emptyflag;
  TTreeReaderValue<Bool_t>* pimflag;
  TTreeReaderValue<Bool_t>* lpiflag;
  TTreeReaderValue<Bool_t>* lpflag;
  TTreeReaderValue<Bool_t>* lflag;
  TTreeReaderValue<Bool_t>* kuramalflag;

  TTreeReaderValue<Double_t>* lmass;
  TTreeReaderValue<Double_t>* ldecayvtx_x;
  TTreeReaderValue<Double_t>* ldecayvtx_y;
  TTreeReaderValue<Double_t>* ldecayvtx_z;
  TTreeReaderValue<Double_t>* lmom;
  TTreeReaderValue<Double_t>* lmom_x;
  TTreeReaderValue<Double_t>* lmom_y;
  TTreeReaderValue<Double_t>* lmom_z;
  TTreeReaderValue<Double_t>* ppi_dist;
  TTreeReaderValue<Double_t>* ltarget_dist;
  TTreeReaderValue<Double_t>* ltargetvtx_x;
  TTreeReaderValue<Double_t>* ltargetvtx_y;
  TTreeReaderValue<Double_t>* ltargetvtx_z;

  TTreeReaderValue<Double_t>* lmass_vtx;
  TTreeReaderValue<Double_t>* ldecayvtx_x_vtx;
  TTreeReaderValue<Double_t>* ldecayvtx_y_vtx;
  TTreeReaderValue<Double_t>* ldecayvtx_z_vtx;
  TTreeReaderValue<Double_t>* lmom_vtx;
  TTreeReaderValue<Double_t>* lmom_x_vtx;
  TTreeReaderValue<Double_t>* lmom_y_vtx;
  TTreeReaderValue<Double_t>* lmom_z_vtx;
  TTreeReaderValue<Double_t>* ppi_dist_vtx;

  TTreeReaderValue<Double_t>* GFlmass;
  TTreeReaderValue<Double_t>* GFldecayvtx_x;
  TTreeReaderValue<Double_t>* GFldecayvtx_y;
  TTreeReaderValue<Double_t>* GFldecayvtx_z;
  TTreeReaderValue<Double_t>* GFlmom;
  TTreeReaderValue<Double_t>* GFlmom_x;
  TTreeReaderValue<Double_t>* GFlmom_y;
  TTreeReaderValue<Double_t>* GFlmom_z;
  TTreeReaderValue<Double_t>* GFppi_dist;
  TTreeReaderValue<Double_t>* GFltarget_dist;
  TTreeReaderValue<Double_t>* GFltargetvtx_x;
  TTreeReaderValue<Double_t>* GFltargetvtx_y;
  TTreeReaderValue<Double_t>* GFltargetvtx_z;
  TTreeReaderValue<Double_t>* GFltargetcenter_x;
  TTreeReaderValue<Double_t>* GFltargetcenter_y;
  TTreeReaderValue<Double_t>* GFltargetcenter_z;
  TTreeReaderValue<Double_t>* GFltargetcenter_dist;
  TTreeReaderValue<Double_t>* GFlprodvtx_x;
  TTreeReaderValue<Double_t>* GFlprodvtx_y;
  TTreeReaderValue<Double_t>* GFlprodvtx_z;
  TTreeReaderValue<Double_t>* GFlprodvtx_dist;
  TTreeReaderValue<Double_t>* GFltracklen;
  TTreeReaderValue<Double_t>* GFltof;
  TTreeReaderValue<Double_t>* GFlctau;

  TTreeReaderValue<Double_t>* GFprodvtx_x_kkxi;
  TTreeReaderValue<Double_t>* GFprodvtx_y_kkxi;
  TTreeReaderValue<Double_t>* GFprodvtx_z_kkxi;

  TTreeReaderValue<std::vector<Int_t>>* GFdecays_pdgcode; 
  TTreeReaderValue<std::vector<Int_t>>* GFdecays_nhtrack; 
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_charge;  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_chisqr;  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_pval;	  
  TTreeReaderValue<std::vector<Int_t>>* GFdecays_htofid;  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_tracklen;
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_tof;	  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_mass2;	  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_invbeta; 
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_mom;	  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_mom_x;	  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_mom_y;	  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_mom_z;	  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_CMmom;	  
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_CMmom_x; 
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_CMmom_y; 
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_CMmom_z; 
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_momloss; 
  TTreeReaderValue<std::vector<Double_t>>* GFdecays_eloss;   

  TTreeReaderValue<std::vector<Int_t>>* decays_id;     
  TTreeReaderValue<std::vector<Double_t>>* decays_mom;    
  TTreeReaderValue<std::vector<Double_t>>* decays_mom_x;  
  TTreeReaderValue<std::vector<Double_t>>* decays_mom_y;  
  TTreeReaderValue<std::vector<Double_t>>* decays_mom_z;  
  TTreeReaderValue<std::vector<Double_t>>* decays_CMmom;  
  TTreeReaderValue<std::vector<Double_t>>* decays_CMmom_x;
  TTreeReaderValue<std::vector<Double_t>>* decays_CMmom_y;
  TTreeReaderValue<std::vector<Double_t>>* decays_CMmom_z;

  TTreeReaderValue<Bool_t>* pipiflag;

  TTreeReaderValue<Int_t>* accident_multi;
  TTreeReaderValue<std::vector<Int_t>>* accident_id;

  TTreeReaderValue<Int_t>* residual_multi;
  TTreeReaderValue<Int_t>* pim_multi;
  TTreeReaderValue<Int_t>* pip_multi;
  TTreeReaderValue<Int_t>* em_multi;
  TTreeReaderValue<Int_t>* ep_multi;
  TTreeReaderValue<Int_t>* p_multi;
  TTreeReaderValue<Int_t>* ppip_multi;
  TTreeReaderValue<std::vector<Int_t>>* residual_id;
  TTreeReaderValue<std::vector<Double_t>>* residual_dist2tgt;
  TTreeReaderValue<std::vector<Double_t>>* residual_GFdist2prodvtx;
  TTreeReaderValue<std::vector<Double_t>>* residual_KFdist2prodvtx;
  TTreeReaderValue<std::vector<Double_t>>* residual_mass2;
  TTreeReaderValue<std::vector<Double_t>>* residual_invbeta;
  TTreeReaderValue<std::vector<Double_t>>* residual_mom;
  TTreeReaderValue<std::vector<Double_t>>* residual_mom_x;
  TTreeReaderValue<std::vector<Double_t>>* residual_mom_y;
  TTreeReaderValue<std::vector<Double_t>>* residual_mom_z;
  TTreeReaderValue<std::vector<Double_t>>* residual_charge;

  TTreeReaderValue<Int_t>* g_multi;
  TTreeReaderValue<std::vector<Int_t>>* epidgamma;
  TTreeReaderValue<std::vector<Int_t>>* emidgamma;
  TTreeReaderValue<std::vector<Double_t>>* epmomgamma;
  TTreeReaderValue<std::vector<Double_t>>* epmomgamma_x;
  TTreeReaderValue<std::vector<Double_t>>* epmomgamma_y;
  TTreeReaderValue<std::vector<Double_t>>* epmomgamma_z;
  TTreeReaderValue<std::vector<Double_t>>* emmomgamma;
  TTreeReaderValue<std::vector<Double_t>>* emmomgamma_x;
  TTreeReaderValue<std::vector<Double_t>>* emmomgamma_y;
  TTreeReaderValue<std::vector<Double_t>>* emmomgamma_z;
  TTreeReaderValue<std::vector<Double_t>>* momgamma;
  TTreeReaderValue<std::vector<Double_t>>* momgamma_x;
  TTreeReaderValue<std::vector<Double_t>>* momgamma_y;
  TTreeReaderValue<std::vector<Double_t>>* momgamma_z;
  TTreeReaderValue<std::vector<Double_t>>* epidistgamma;
  TTreeReaderValue<std::vector<Double_t>>* vtxgamma_x;
  TTreeReaderValue<std::vector<Double_t>>* vtxgamma_y;
  TTreeReaderValue<std::vector<Double_t>>* vtxgamma_z;
  
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
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
  
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;
  
  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  //Initiallize Geometry, Field, Fitter
  //HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  //Initiallize the genfit track container
  //HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  //GFtracks.SetVerbosity(verbosity);
  //std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  //std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;

#if 0
  //GFtracks.DebugMode();
#endif
  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
    //GFtracks.Clear();    
  }
  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;    

  DstClose();
  
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;    
  //delete fitter;
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
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;
    
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
  static const Double_t Carbon12Mass = 12.*TGeoUnit::amu_c2 - 6.*ElectronMass;
  static const Double_t Boron11Mass  = 11.009305167*TGeoUnit::amu_c2 - 5.*ElectronMass;  
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};
  TVector3 tgtpos(0, 0, tpc::ZTarget);
  TVector3 qnan_vec = TVector3(qnan, qnan, qnan);
  
  static const auto KKEvent = gUser.GetParameter("KKEvent");
  static const auto KPEvent = gUser.GetParameter("KPEvent");
  static const auto KHeavyEvent = gUser.GetParameter("KHeavyEvent");
  
  
  Double_t vtx_scan_range = gUser.GetParameter("VertexScanRange");
  use_pidlikeli = gUser.GetParameter("UsePidLikeli");  

  //if( ievent%1000==0 ){
  if( ievent%100==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  //event.status = **src.status;
  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;

  event.ntK18 = **src.ntK18;
  event.pK18 = **src.pK18;
  event.chisqrK18 = **src.chisqrK18;
  event.xtgtK18 = **src.xtgtK18;
  event.ytgtK18 = **src.ytgtK18;
  event.utgtK18 = **src.utgtK18;
  event.vtgtK18 = **src.vtgtK18;

  event.ntKurama = **src.ntKurama;
  event.chisqrKurama = **src.chisqrKurama;
  event.pKurama = **src.pKurama;
  event.qKurama = **src.qKurama;
  //event.m2Kurama = **src.m2Kurama;
  event.thetaKurama = **src.thetaKurama;
  event.xtgtKurama = **src.xtgtKurama;
  event.ytgtKurama = **src.ytgtKurama;
  event.utgtKurama = **src.utgtKurama;
  event.vtgtKurama = **src.vtgtKurama;
  event.pathwcKurama = **src.pathwcKurama;
  event.xin = **src.xin;
  event.yin = **src.yin;
  event.zin = **src.zin;
  event.pxin = **src.pxin;
  event.pyin = **src.pyin;
  event.pzin = **src.pzin;
  // event.xout = **src.xout;
  // event.yout = **src.yout;
  // event.zout = **src.zout;
  // event.pxout = **src.pxout;
  // event.pyout = **src.pyout;
  // event.pzout = **src.pzout;

  // event.nKm = **src.nKm;
  // event.nKp = **src.nKp;
  event.nKK = **src.nKK;
  // event.inside = **src.inside;
  event.vtx = **src.vtx;
  event.vty = **src.vty;
  event.vtz = **src.vtz;
  // event.closeDist = **src.closeDist;
  event.MissMass = **src.MissMass;
  event.MissMassCorr = **src.MissMassCorr;
  event.MissMassCorrDE = **src.MissMassCorrDE;
  event.pOrg = **src.pOrg;
  event.pCalc = **src.pCalc;
  event.pCorr = **src.pCorr;
  event.pCorrDE = **src.pCorrDE;
  event.xb = **src.xb;
  event.yb = **src.yb;
  event.ub = **src.ub;
  event.vb = **src.vb;
  event.xs = **src.xs;
  event.ys = **src.ys;
  event.us = **src.us;
  event.vs = **src.vs;
  event.Kflag = **src.Kflag;
  //event.Pflag = **src.Pflag;
  //event.Heavyflag = **src.Heavyflag;

  // DstE42からの情報 (TPC-RK)
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

  event.isgoodTPCKurama = **src.isgoodTPCKurama;
  event.kflagTPCKurama = **src.kflagTPCKurama;
  //event.pflagTPCKurama = **src.pflagTPCKurama;
  event.chisqrTPCKurama = **src.chisqrTPCKurama;
  event.pTPCKurama = **src.pTPCKurama;
  event.qTPCKurama = **src.qTPCKurama;
  event.m2TPCKurama = **src.m2TPCKurama;
  event.xtgtTPCKurama = **src.xtgtTPCKurama;
  event.ytgtTPCKurama = **src.ytgtTPCKurama;
  event.utgtTPCKurama = **src.utgtTPCKurama;
  event.vtgtTPCKurama = **src.vtgtTPCKurama;
  event.thetaTPCKurama = **src.thetaTPCKurama;
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
  // event.MissMassNuclTPC = **src.MissMassNuclTPC;
  // event.MissMassNuclCorrTPC = **src.MissMassNuclCorrTPC;
  // event.MissMassNuclCorrDETPC = **src.MissMassNuclCorrDETPC;
  event.pOrgTPC = **src.pOrgTPC;
  event.pCalcTPC = **src.pCalcTPC;
  event.pCorrTPC = **src.pCorrTPC;
  event.pCorrDETPC = **src.pCorrDETPC;
  event.thetaTPC = **src.thetaTPC;
  //event.thetaCMTPC = **src.thetaCMTPC;
  //event.costCMTPC = **src.costCMTPC;
  event.xbTPC = **src.xbTPC;
  event.ybTPC = **src.ybTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;
  event.xsTPC = **src.xsTPC;
  event.ysTPC = **src.ysTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

  event.BE = **src.BE;
  event.BETPC = **src.BETPC;
  event.BE_LL = **src.BE_LL;
  event.BETPC_LL = **src.BETPC_LL;
  event.km_mom_x = **src.km_mom_x;
  event.km_mom_y = **src.km_mom_y;
  event.km_mom_z = **src.km_mom_z;
  event.kp_mom_x = **src.kp_mom_x;
  event.kp_mom_y = **src.kp_mom_y;
  event.kp_mom_z = **src.kp_mom_z;

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

  event.remain_nclTpc = **src.remain_nclTpc;
  event.remain_cluster_x = **src.remain_cluster_x;
  event.remain_cluster_y = **src.remain_cluster_y;
  event.remain_cluster_z = **src.remain_cluster_z;
  event.remain_cluster_de = **src.remain_cluster_de;
  event.remain_cluster_size = **src.remain_cluster_size;
  event.remain_cluster_layer = **src.remain_cluster_layer;
  event.remain_cluster_mrow = **src.remain_cluster_mrow;
  event.remain_cluster_de_center = **src.remain_cluster_de_center;
  event.remain_cluster_x_center = **src.remain_cluster_x_center;
  event.remain_cluster_y_center = **src.remain_cluster_y_center;
  event.remain_cluster_z_center = **src.remain_cluster_z_center;
  event.remain_cluster_row_center = **src.remain_cluster_row_center;
  event.remain_cluster_houghflag = **src.remain_cluster_houghflag;

  event.ntTpc = **src.ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isXi = **src.isXi;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.isMultiloop = **src.isMultiloop;
  event.isInTarget = **src.isInTarget;
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
  event.isElectron = **src.isElectron;
  event.nsigma_triton = **src.nsigma_triton;
  event.nsigma_deutron = **src.nsigma_deutron;
  event.nsigma_proton = **src.nsigma_proton;
  event.nsigma_kaon = **src.nsigma_kaon;
  event.nsigma_pion = **src.nsigma_pion;
  event.nsigma_electron = **src.nsigma_electron;

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
  event.chisqr_inverted = **src.chisqr_inverted;
  event.pval_inverted = **src.pval_inverted;
  event.helix_cx_inverted = **src.helix_cx_inverted;
  event.helix_cy_inverted = **src.helix_cy_inverted;
  event.helix_z0_inverted = **src.helix_z0_inverted;
  event.helix_r_inverted = **src.helix_r_inverted;
  event.helix_dz_inverted = **src.helix_dz_inverted;
  event.mom0_inverted = **src.mom0_inverted;
  event.pid_inverted = **src.pid_inverted;

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

  event.ncombiLreconfailed = **src.ncombiLreconfailed;
  event.pidLreconfailed = **src.pidLreconfailed;
  event.piidLreconfailed = **src.piidLreconfailed;
  event.LdecayvtxLreconfailed_x = **src.LdecayvtxLreconfailed_x;
  event.LdecayvtxLreconfailed_y = **src.LdecayvtxLreconfailed_y;
  event.LdecayvtxLreconfailed_z = **src.LdecayvtxLreconfailed_z;
  event.LmassLreconfailed = **src.LmassLreconfailed;
  event.LmomLreconfailed = **src.LmomLreconfailed;
  event.LmomLreconfailed_x = **src.LmomLreconfailed_x;
  event.LmomLreconfailed_y = **src.LmomLreconfailed_y;
  event.LmomLreconfailed_z = **src.LmomLreconfailed_z;
  event.pmomLreconfailed = **src.pmomLreconfailed;
  event.pmomLreconfailed_x = **src.pmomLreconfailed_x;
  event.pmomLreconfailed_y = **src.pmomLreconfailed_y;
  event.pmomLreconfailed_z = **src.pmomLreconfailed_z;
  event.pimomLreconfailed = **src.pimomLreconfailed;
  event.pimomLreconfailed_x = **src.pimomLreconfailed_x;
  event.pimomLreconfailed_y = **src.pimomLreconfailed_y;
  event.pimomLreconfailed_z = **src.pimomLreconfailed_z;
  event.ppidistLreconfailed = **src.ppidistLreconfailed;
  
  event.nEscapeKm = **src.nEscapeKm;
  event.kmid = **src.kmid;
  event.kmmom = **src.kmmom;
  event.kmmom_x = **src.kmmom_x;
  event.kmmom_y = **src.kmmom_y;
  event.kmmom_z = **src.kmmom_z;
  
  event.GFkmdecayvtx_x = **src.GFkmdecayvtx_x;
  event.GFkmdecayvtx_y = **src.GFkmdecayvtx_y;
  event.GFkmdecayvtx_z = **src.GFkmdecayvtx_z;
  event.GFkmmom = **src.GFkmmom;
  event.GFkmmom_x = **src.GFkmmom_x;
  event.GFkmmom_y = **src.GFkmmom_y;
  event.GFkmmom_z = **src.GFkmmom_z;
  event.GFkmtarget_dist = **src.GFkmtarget_dist;
  event.GFkmtargetvtx_x = **src.GFkmtargetvtx_x;
  event.GFkmtargetvtx_y = **src.GFkmtargetvtx_y;
  event.GFkmtargetvtx_z = **src.GFkmtargetvtx_z;
  event.GFkmtargetcenter_x = **src.GFkmtargetcenter_x;
  event.GFkmtargetcenter_y = **src.GFkmtargetcenter_y;
  event.GFkmtargetcenter_z = **src.GFkmtargetcenter_z;
  event.GFkmtargetcenter_dist = **src.GFkmtargetcenter_dist;
  event.GFkmm2 = **src.GFkmm2;  
  event.GFkmtracklen = **src.GFkmtracklen;
  event.GFkmtof = **src.GFkmtof;  

  event.ncombiPipair = **src.ncombiPipair;
  event.pipidPipair = **src.pipidPipair;
  event.pimidPipair = **src.pimidPipair;
  event.pipmomPipair = **src.pipmomPipair;
  event.pipmomPipair_x = **src.pipmomPipair_x;
  event.pipmomPipair_y = **src.pipmomPipair_y;
  event.pipmomPipair_z = **src.pipmomPipair_z;
  event.pimmomPipair = **src.pimmomPipair;
  event.pimmomPipair_x = **src.pimmomPipair_x;
  event.pimmomPipair_y = **src.pimmomPipair_y;
  event.pimmomPipair_z = **src.pimmomPipair_z;
  event.momPipair = **src.momPipair;
  event.momPipair_x = **src.momPipair_x;
  event.momPipair_y = **src.momPipair_y;
  event.momPipair_z = **src.momPipair_z;
  event.reconLmassPipair = **src.reconLmassPipair;
  event.reconmassPipair = **src.reconmassPipair;
  event.pipidistPipair = **src.pipidistPipair;

  // Genfit Track's information
  event.GFfitstatus = **src.GFfitstatus;
  event.GFpdgcode = **src.GFpdgcode;
  event.GFnhtrack = **src.GFnhtrack;
  event.GFcharge = **src.GFcharge;
  event.GFchisqr = **src.GFchisqr;
  event.GFtof = **src.GFtof;
  event.GFpval = **src.GFpval;
  event.GFlayer = **src.GFlayer;
  event.GFpos_x = **src.GFpos_x;
  event.GFpos_y = **src.GFpos_y;
  event.GFpos_z = **src.GFpos_z;
  event.GFmom = **src.GFmom;
  event.GFmom_x = **src.GFmom_x;
  event.GFmom_y = **src.GFmom_y;
  event.GFmom_z = **src.GFmom_z;
  event.GFresidual_x = **src.GFresidual_x;
  event.GFresidual_y = **src.GFresidual_y;
  event.GFresidual_z = **src.GFresidual_z;
  event.GFresidual_p = **src.GFresidual_p;
  event.GFresidual_px = **src.GFresidual_px;
  event.GFresidual_py = **src.GFresidual_py;
  event.GFresidual_pz = **src.GFresidual_pz;

  event.GFinside = **src.GFinside;
  event.GFKuramaFromTgt = **src.GFKuramaFromTgt;
  event.GFKuramaVtxOutTgt = **src.GFKuramaVtxOutTgt;  

  event.GFntTpc_target = **src.GFntTpc_target;    
  event.GFprodvtx_x = **src.GFprodvtx_x;
  event.GFprodvtx_y = **src.GFprodvtx_y;
  event.GFprodvtx_z = **src.GFprodvtx_z;

  event.GFtracklen = **src.GFtracklen;
  event.GFtrack2vtxdist = **src.GFtrack2vtxdist;
  event.GFcalctof = **src.GFcalctof;
  event.GFsegHtof = **src.GFsegHtof;
  event.GFtofHtof = **src.GFtofHtof;
  event.GFtdiffHtof = **src.GFtdiffHtof;
  event.GFposHtof = **src.GFposHtof;
  event.GFposx = **src.GFposx;
  event.GFposy = **src.GFposy;
  event.GFposz = **src.GFposz;
  event.GFinvbeta = **src.GFinvbeta;
  event.GFm2 = **src.GFm2;
  event.nsigma_tritonHtof = **src.nsigma_tritonHtof;
  event.nsigma_deutronHtof = **src.nsigma_deutronHtof;
  event.nsigma_protonHtof = **src.nsigma_protonHtof;
  event.nsigma_kaonHtof = **src.nsigma_kaonHtof;
  event.nsigma_pionHtof = **src.nsigma_pionHtof;
  event.nsigma_electronHtof = **src.nsigma_electronHtof;

  // event.GFprodvtx_x_ll = **src.GFprodvtx_x_ll;
  // event.GFprodvtx_y_ll = **src.GFprodvtx_y_ll;
  // event.GFprodvtx_z_ll = **src.GFprodvtx_z_ll;
  // event.GFprodvtx_x_l1 = **src.GFprodvtx_x_l1;
  // event.GFprodvtx_y_l1 = **src.GFprodvtx_y_l1;
  // event.GFprodvtx_z_l1 = **src.GFprodvtx_z_l1;
  // event.GFprodvtx_x_l2 = **src.GFprodvtx_x_l2;
  // event.GFprodvtx_y_l2 = **src.GFprodvtx_y_l2;
  // event.GFprodvtx_z_l2 = **src.GFprodvtx_z_l2;
  // event.GFprodvtx_x_l = **src.GFprodvtx_x_l;
  // event.GFprodvtx_y_l = **src.GFprodvtx_y_l;
  // event.GFprodvtx_z_l = **src.GFprodvtx_z_l;

  event.emptyflag = **src.emptyflag;
  event.pimflag = **src.pimflag;
  event.lpiflag = **src.lpiflag;
  event.lpflag = **src.lpflag;
  event.lflag = **src.lflag;
  event.kuramalflag = **src.kuramalflag;

  event.lmass = **src.lmass;
  event.ldecayvtx_x = **src.ldecayvtx_x;
  event.ldecayvtx_y = **src.ldecayvtx_y;
  event.ldecayvtx_z = **src.ldecayvtx_z;
  event.lmom = **src.lmom;
  event.lmom_x = **src.lmom_x;
  event.lmom_y = **src.lmom_y;
  event.lmom_z = **src.lmom_z;
  event.ppi_dist = **src.ppi_dist;
  event.ltarget_dist = **src.ltarget_dist;
  event.ltargetvtx_x = **src.ltargetvtx_x;
  event.ltargetvtx_y = **src.ltargetvtx_y;
  event.ltargetvtx_z = **src.ltargetvtx_z;

  // event.lmass_vtx = **src.lmass_vtx;
  // event.ldecayvtx_x_vtx = **src.ldecayvtx_x_vtx;
  // event.ldecayvtx_y_vtx = **src.ldecayvtx_y_vtx;
  // event.ldecayvtx_z_vtx = **src.ldecayvtx_z_vtx;
  // event.lmom_vtx = **src.lmom_vtx;
  // event.lmom_x_vtx = **src.lmom_x_vtx;
  // event.lmom_y_vtx = **src.lmom_y_vtx;
  // event.lmom_z_vtx = **src.lmom_z_vtx;
  // event.ppi_dist_vtx = **src.ppi_dist_vtx;

  event.GFlmass = **src.GFlmass;
  event.GFldecayvtx_x = **src.GFldecayvtx_x;
  event.GFldecayvtx_y = **src.GFldecayvtx_y;
  event.GFldecayvtx_z = **src.GFldecayvtx_z;
  event.GFlmom = **src.GFlmom;
  event.GFlmom_x = **src.GFlmom_x;
  event.GFlmom_y = **src.GFlmom_y;
  event.GFlmom_z = **src.GFlmom_z;
  event.GFppi_dist = **src.GFppi_dist;
  event.GFltarget_dist = **src.GFltarget_dist;
  event.GFltargetvtx_x = **src.GFltargetvtx_x;
  event.GFltargetvtx_y = **src.GFltargetvtx_y;
  event.GFltargetvtx_z = **src.GFltargetvtx_z;
  event.GFltargetcenter_x = **src.GFltargetcenter_x;
  event.GFltargetcenter_y = **src.GFltargetcenter_y;
  event.GFltargetcenter_z = **src.GFltargetcenter_z;
  event.GFltargetcenter_dist = **src.GFltargetcenter_dist;
  event.GFlprodvtx_x = **src.GFlprodvtx_x;
  event.GFlprodvtx_y = **src.GFlprodvtx_y;
  event.GFlprodvtx_z = **src.GFlprodvtx_z;
  event.GFlprodvtx_dist = **src.GFlprodvtx_dist;
  event.GFltracklen = **src.GFltracklen;
  event.GFltof = **src.GFltof;
  event.GFlctau = **src.GFlctau;

  // event.GFprodvtx_x_kkxi = **src.GFprodvtx_x_kkxi;
  // event.GFprodvtx_y_kkxi = **src.GFprodvtx_y_kkxi;
  // event.GFprodvtx_z_kkxi = **src.GFprodvtx_z_kkxi;

  event.GFdecays_pdgcode = **src.GFdecays_pdgcode;
  event.GFdecays_nhtrack = **src.GFdecays_nhtrack;
  event.GFdecays_charge = **src.GFdecays_charge;
  event.GFdecays_chisqr = **src.GFdecays_chisqr;
  event.GFdecays_pval = **src.GFdecays_pval;
  event.GFdecays_htofid = **src.GFdecays_htofid;
  event.GFdecays_tracklen = **src.GFdecays_tracklen;
  event.GFdecays_tof = **src.GFdecays_tof;
  event.GFdecays_mass2 = **src.GFdecays_mass2;
  event.GFdecays_invbeta = **src.GFdecays_invbeta;
  event.GFdecays_mom = **src.GFdecays_mom;
  event.GFdecays_mom_x = **src.GFdecays_mom_x;
  event.GFdecays_mom_y = **src.GFdecays_mom_y;
  event.GFdecays_mom_z = **src.GFdecays_mom_z;
  event.GFdecays_CMmom = **src.GFdecays_CMmom;
  event.GFdecays_CMmom_x = **src.GFdecays_CMmom_x;
  event.GFdecays_CMmom_y = **src.GFdecays_CMmom_y;
  event.GFdecays_CMmom_z = **src.GFdecays_CMmom_z;
  event.GFdecays_momloss = **src.GFdecays_momloss;
  event.GFdecays_eloss = **src.GFdecays_eloss;

  event.decays_id = **src.decays_id;
  event.decays_mom = **src.decays_mom;
  event.decays_mom_x = **src.decays_mom_x;
  event.decays_mom_y = **src.decays_mom_y;
  event.decays_mom_z = **src.decays_mom_z;
  event.decays_CMmom = **src.decays_CMmom;
  event.decays_CMmom_x = **src.decays_CMmom_x;
  event.decays_CMmom_y = **src.decays_CMmom_y;
  event.decays_CMmom_z = **src.decays_CMmom_z;

  event.pipiflag = **src.pipiflag;

  event.accident_multi = **src.accident_multi;
  event.accident_id = **src.accident_id;

  event.residual_multi = **src.residual_multi;
  event.pim_multi = **src.pim_multi;
  event.pip_multi = **src.pip_multi;
  event.em_multi = **src.em_multi;
  event.ep_multi = **src.ep_multi;
  event.p_multi = **src.p_multi;
  event.ppip_multi = **src.ppip_multi;
  event.residual_id = **src.residual_id;
  event.residual_dist2tgt = **src.residual_dist2tgt;
  event.residual_GFdist2prodvtx = **src.residual_GFdist2prodvtx;
  event.residual_KFdist2prodvtx = **src.residual_KFdist2prodvtx;
  event.residual_mass2 = **src.residual_mass2;
  event.residual_invbeta = **src.residual_invbeta;
  event.residual_mom = **src.residual_mom;
  event.residual_mom_x = **src.residual_mom_x;
  event.residual_mom_y = **src.residual_mom_y;
  event.residual_mom_z = **src.residual_mom_z;
  event.residual_charge = **src.residual_charge;

  event.g_multi = **src.g_multi;
  event.epidgamma = **src.epidgamma;
  event.emidgamma = **src.emidgamma;
  event.epmomgamma = **src.epmomgamma;
  event.epmomgamma_x = **src.epmomgamma_x;
  event.epmomgamma_y = **src.epmomgamma_y;
  event.epmomgamma_z = **src.epmomgamma_z;
  event.emmomgamma = **src.emmomgamma;
  event.emmomgamma_x = **src.emmomgamma_x;
  event.emmomgamma_y = **src.emmomgamma_y;
  event.emmomgamma_z = **src.emmomgamma_z;
  event.momgamma = **src.momgamma;
  event.momgamma_x = **src.momgamma_x;
  event.momgamma_y = **src.momgamma_y;
  event.momgamma_z = **src.momgamma_z;
  event.epidistgamma = **src.epidistgamma;
  event.vtxgamma_x = **src.vtxgamma_x;
  event.vtxgamma_y = **src.vtxgamma_y;
  event.vtxgamma_z = **src.vtxgamma_z;
  
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;
  Int_t ntTpc = event.ntTpc;  
  HF1( 1, event.status++ );  
  HF1( 10, ntTpc );  
  if( event.ntTpc == 0 ) return true;
  HF1( 1, event.status++ );  


  HF1( 1, event.status++ );  
  if( event.nKK != 1 ) return true;
  HF1( 1, event.status++ );
  
  double BEkaon = 0.;
  double thetaTPC = 0.;
  double MissMass = 0.;
  for(Int_t iKK=0; iKK<event.nKK; iKK++){
    if(KKEvent&&event.kflagTPCKurama[iKK]!=1) return true;
    //BEkaon = event.MissMassNuclCorrDETPC[iKK] - KaonMass - Boron11Mass - 0.075;
    thetaTPC = event.thetaTPC[0];
    MissMass = event.MissMassCorrDETPC[iKK];
  }

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Int_t psfac = 0;
  bool trigA = (event.trigflag[20]>0);
  bool trigB = (event.trigflag[21]>0);  
  if(trigA) psfac = psTrigA;
  else psfac = 0;
  
  std::cout << " debug " << __FILE__ << " " << __LINE__
	    << " ntTpc: " << ntTpc << std::endl;
  // spectrum decomposition
  double be = -event.BETPC[0]-75;
  int BErgn = -1;
  int BE1rgn = -1;  
  if(be>-200&&be<0){
    BErgn = 0;
    if(be>-200&&be<-140){
      BE1rgn = 3;
    } else if(be<-60){
      BE1rgn = 4;
    } else {
      BE1rgn = 5;
    }
  } else if(be<180){ 
    BErgn = 1;
  } else if(be<400){
    BErgn = 2;
  }
  // additional particles
  int numpim=0; int numpip=0; int numkm=0; int numkp=0; int nump=0; int numem=0; int numep=0; int numppip=0;
  int numk18=0; int numkurama=0; int numbeam=0; int numacc=0;
  for(int it=0; it<ntTpc; it++){
    if( event.isK18[it]==1 ){
      numk18++;
    } else if( event.isKurama[it]==1 ) {
      numkurama++;
    } else if( event.isBeam[it]==1 ) {
      numbeam++;
    } else if( event.isAccidental[it]==1 ) {
      numacc++;
    }
    
    if ( event.isK18[it]==1 || event.isKurama[it]==1 || event.isBeam[it]==1 || event.isAccidental[it]==1 ) continue;
    if ( event.lflag ){
      if (it==event.decays_id[0] || it==event.decays_id[1] ) continue;
    }
    bool eflag = false;
    for( int ie=0; ie<event.epidgamma.size(); ie++){
      if ( it==event.epidgamma[ie] || it==event.emidgamma[ie] ) eflag=true;
    }
    if(eflag) continue;
    Int_t nhit = event.nhtrack[it];
    TVector3 start(event.hitpos_x[it][0], event.hitpos_y[it][0], event.hitpos_z[it][0]);
    TVector3 end(event.hitpos_x[it][nhit-1], event.hitpos_y[it][nhit-1], event.hitpos_z[it][nhit-1]); 
    
    bool intarget=false;
    for(int iresi=0; iresi<event.residual_id.size(); iresi++){
      if(it != event.residual_id[iresi]) continue;
      else{
	intarget = true;
	break;
      }
    }
    if(intarget) continue;
    
    // follow pid criteria in GenfitE42
    if((event.pid[it]&2)==2 && event.charge[it]==-1){ //k-
      numkm++;
      continue;      
    }
    if(event.isElectron[it]==1){ //e+, e-
      if(event.charge[it]==1){
	numep++;
	continue;
      } else {
	numem++;
	continue;
      }
    }    
    else if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton
      nump++;
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
	//if(TMath::Abs(slope)<0.1 && TMath::Abs(helixmom)>0.5 &&
	(start.x()-end.x())>-10 && (start.x()-end.x())<50.){
	event.isAccidental[it] = 1;
	//target_accidental_id_container.push_back(it); //Accidental beam on the target
	continue; //Accidental K-
      }      
      numpim++;
      continue;
    }
    else if((event.pid[it]&4)!=4 && (event.pid[it]&1)==1 && event.charge[it]==1){ //pi+
      numpip++;
      continue;
    }
    else if(((event.pid[it]&4)==4 || (event.pid[it]&1)==1) && event.charge[it]==1){ //p or pi+ with high-mom
      numppip++;
      continue;
    }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;    
  }
  
  Counts n = {}; RecoCounts nReco = {};
  n[(int)Species::PIP] = event.pip_multi+numpip;
  n[(int)Species::KM] = event.nEscapeKm;
  n[(int)Species::KP] = 0+numkp;
  n[(int)Species::EM] = event.em_multi+numem;
  n[(int)Species::EP] = event.ep_multi+numep;
  n[(int)Species::PPIP] = event.ppip_multi+numppip;
  if( event.lflag==1
      &&event.GFlmass>lmass_cut_min&&event.GFlmass<lmass_cut_max
      &&event.GFlctau>lctau_cut_min ){
    nReco[(int)Reco::LMD] = 1;
    n[(int)Species::P] = event.p_multi+nump;
    n[(int)Species::PIM] = event.pim_multi+numpim;
  } else if (event.lflag==1){
    n[(int)Species::P] = event.p_multi+1+nump;
    n[(int)Species::PIM] = event.pim_multi+1+numpim;
  } else {
    n[(int)Species::P] = event.p_multi+nump;
    n[(int)Species::PIM] = event.pim_multi+numpim;
  }
  if( event.g_multi>0 ){
    n[(int)Species::EM] = event.em_multi-event.g_multi+numem;        
    n[(int)Species::EP] = event.ep_multi-event.g_multi+numep;    
  }
  nReco[(int)Reco::K0] = 0;
  nReco[(int)Reco::GAMMA] = event.g_multi;

  int reco_code = encode_reco_code(nReco);
  int reco_id = encode_reco_id_1based(nReco);
  uint32_t pidcombi = encode_pid_code(n);
  //std::cout << "accidental: " << event.accident_multi << std::endl;
  std::cout << "K18: " << numk18 << " Kurama:" << numkurama
	    << " beam:" << numbeam << " accidental:" << numacc << std::endl;  
  std::string s = decimal_to_base5(pidcombi);
  std::cout << pidcombi << " (decimal) = " << s << " (base5)" << std::endl;
  std::cout << "recoid: " << reco_code << std::endl;
  for ( int ips=0; ips<psfac; ips++ ) {
    HF1(BErgn*100+reco_id+100, pidcombi);
    if(BE1rgn>0){
      HF1(BE1rgn*100+reco_id+100, pidcombi);
    }  
    if( thetaTPC<10 ){
      HF1(hid_inc, be);
      if( nReco[(int)Reco::LMD]==0&&event.g_multi==0 ){
	if( pidcombi==base5_to_decimal("0") ){ // empty      
	  HF1(hid_exc+pidcombi,be);
	} else if( pidcombi==base5_to_decimal("1") ){
	  HF1(hid_exc+pidcombi,be);    
	} else if( pidcombi==base5_to_decimal("2") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("10") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("11") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("100") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("101") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("102") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("10000") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("10001") ){
	  HF1(hid_exc+pidcombi,be);        
	} else if( pidcombi==base5_to_decimal("10002") ){
	  HF1(hid_exc+pidcombi,be);
	} else {
	  HF1(hid_exc+10000, be); 
	}
      } else if( nReco[(int)Reco::LMD]==1&&event.g_multi==0&&pidcombi==base5_to_decimal("0") ){
	HF1(hid_exc_Lmd+pidcombi,be);
      } else if( nReco[(int)Reco::LMD]==0&&event.g_multi==1&&pidcombi==base5_to_decimal("11") ){
	HF1(hid_exc_gamma+pidcombi,be);
      } else {
	HF1(hid_exc+10000, be); 
      }
      // semi-exclusive
      if(n[(int)Species::PIM]>=1){
	int id = base5_to_decimal("1");
	HF1(hid_semiexc+id,be);
      }
      if(n[(int)Species::PIP]>=1){
	int id = base5_to_decimal("10");
	HF1(hid_semiexc+id,be);	
      }
      if(n[(int)Species::KM]>=1){
	int id = base5_to_decimal("100");
	HF1(hid_semiexc+id,be);	
      }
      if(n[(int)Species::P]>=1){
	int id = base5_to_decimal("10000");
	HF1(hid_semiexc+id,be);	
      }
      if(nReco[(int)Reco::LMD]>=1){
	HF1(hid_semiexc+20000,be);	
      }
      if(n[(int)Species::PIM]>=1||n[(int)Species::KM]>=1){ // nega
	HF1(hid_semiexc+10000,be);	
      }            
      if(n[(int)Species::PIP]>=1||n[(int)Species::P]>=1){ // posi
	HF1(hid_semiexc+10001,be);	
      }      
    }
  }

  if(event.emptyflag){
    std::cout << " emptyflag:" << event.emptyflag << " pidcombi:" << pidcombi << std::endl;
  }
  if(pidcombi==0){
    std::cout << " emptyflag:" << event.emptyflag << " pidcombi:" << pidcombi << std::endl;
    event.emptyflag=true;
  } else {
    event.emptyflag=false;    
  }
  
  double mint = 3.5;
  double maxt = 4.5;
  int lcandi = 1;
  {
    HF1(11201,BEkaon);
    if(thetaTPC>mint&&thetaTPC<maxt) HF1(11202,BEkaon); 
    if(thetaTPC<10) HF1(11203,BEkaon);            
    HF1(12001,event.GFlmass);                    
    if(BEkaon<-0.1) HF1(12002,event.GFlmass);     
    else if(BEkaon<0.) HF1(12003,event.GFlmass);  
    else if(BEkaon<0.1) HF1(12004,event.GFlmass);         
    else if(BEkaon<0.2) HF1(12005,event.GFlmass); 
    else if(BEkaon<0.3) HF1(12006,event.GFlmass);              
  }
  for(int icandi=0; icandi<lcandi; icandi++){
    if(!event.lflag) continue;
    HF1(11204,BEkaon);
    if(thetaTPC>mint&&thetaTPC<maxt) HF1(11205,BEkaon);
    if(thetaTPC<10) HF1(11206,BEkaon);
    for(int i=0;i<2;i++){
      double id = event.decays_id[i];
      double mom = event.GFdecays_mom[i];      
      double ch = event.charge[id];
      double dedx = event.dEdx[id];
      if(TMath::Abs(event.GFlmass-LambdaMass)<0.01){
	HF2(3100,ch*mom,dedx); HF2(3101+i,ch*mom,dedx);
      } else {
	HF2(3110,ch*mom,dedx); HF2(3111+i,ch*mom,dedx);
      }
    }
    // if(event.GFldecays_mass2[0]<min_mass2_p||event.GFldecays_mass2[0]>max_mass2_p) continue; 
    // if(event.GFldecays_mass2[1]<min_mass2_pi||event.GFldecays_mass2[1]>max_mass2_pi) continue; 
    for(int it=0; it<ntTpc; it++){ // proton 
      if ( !event.GFfitstatus[it] ) continue; 
      if ( event.isElectron[it]==1 ) continue; 
      if ( event.isK18[it]==1 ) continue; 
      if ( event.isKurama[it]==1 ) continue;
      if ( event.isBeam[it]==1 ) continue;
      if ( event.isAccidental[it]==1 ) continue;
      if ( event.charge[it]!=1 ) continue;
      if ( it==event.decays_id[0] ) continue;
      if ( (event.pid[it]&4)!=4 ) continue;
      HF1(11207,BEkaon);
      if(thetaTPC>mint&&thetaTPC<maxt) HF1(11208,BEkaon);
      if(thetaTPC<10) HF1(11209,BEkaon);
      HF1(12011,event.GFlmass);
      if(BEkaon<-0.1) HF1(12012,event.GFlmass);
      else if(BEkaon<0. )  HF1(12013,event.GFlmass);
      else if(BEkaon<0.1) HF1(12014,event.GFlmass);        
      else if(BEkaon<0.2) HF1(12015,event.GFlmass);
      else if(BEkaon<0.3) HF1(12016,event.GFlmass);
      // TVector3 lmom(event.GFldecays_mom_x[0][0]+event.GFldecays_mom_x[1][0],
      // 		    event.GFldecays_mom_y[0][0]+event.GFldecays_mom_y[1][0],
      // 		    event.GFldecays_mom_z[0][0]+event.GFldecays_mom_z[1][0]);
      TVector3 lmom(event.decays_mom_x[0]+event.decays_mom_x[1],
		    event.decays_mom_y[0]+event.decays_mom_y[1],
		    event.decays_mom_z[0]+event.decays_mom_z[1]);      
      TLorentzVector GFLlmd(lmom, LambdaMass);
      TVector3 p2mom(event.GFmom_x[it][0],event.GFmom_y[it][0],event.GFmom_z[it][0]);
      TLorentzVector GFLp2(p2mom, ProtonMass);
      TLorentzVector GFLlp2 = GFLp2 + GFLlmd;
      double lp2mass = GFLlp2.M();
      if(BEkaon<-0.1) HF1(12022,lp2mass);
      else if(BEkaon<0. ) HF1(12023,lp2mass);
      else if(BEkaon<0.1) HF1(12024,lp2mass);        
      else if(BEkaon<0.2) HF1(12025,lp2mass);
      else if(BEkaon<0.3) HF1(12026,lp2mass);      
    }
  }
  // with m2 cut 
  for(int icandi=0; icandi<lcandi; icandi++){
    if(!event.lflag) continue;
    // if(event.GFldecays_mass2[0]<min_mass2_p||event.GFldecays_mass2[0]>max_mass2_p) continue;
    // if(event.GFldecays_mass2[1]<min_mass2_pi||event.GFldecays_mass2[1]>max_mass2_pi) continue;    
    HF1(11304,BEkaon);
    if(thetaTPC>mint&&thetaTPC<maxt) HF1(11305,BEkaon);
    if(thetaTPC<10) HF1(11306,BEkaon);
    for(int it=0; it<ntTpc; it++){ // proton
      if(!event.GFfitstatus[it]) continue;
      if(event.isElectron[it]==1) continue; 
      if(event.isK18[it]==1) continue; 
      if(event.isKurama[it]==1) continue;
      if(event.isBeam[it]==1) continue;
      if(event.isAccidental[it]==1) continue;
      if(event.charge[it]!=1) continue;
      if(it==event.decays_id[0]) continue;
      // if((event.pid[it]&4)!=4) continue;
      double nsigmap = event.nsigma_proton[it];
      if(nsigmap<-2.0||nsigmap>2.0) continue;
      //if(!event.GFfromVtx[it]) continue;
      double mass2p = event.GFm2[it];
      if(mass2p<min_mass2_p||mass2p>max_mass2_p) continue; // proton found      
      HF1(11307,BEkaon);
      if(thetaTPC>mint&&thetaTPC<maxt) HF1(11308,BEkaon);
      if(thetaTPC<10) HF1(11309,BEkaon);
      HF1(12311,event.GFlmass);
      if(BEkaon<-0.1) HF1(12312,event.GFlmass);
      else if(BEkaon<0. ) HF1(12313,event.GFlmass);
      else if(BEkaon<0.1) HF1(12314,event.GFlmass);
      else if(BEkaon<0.2) HF1(12315,event.GFlmass);
      else if(BEkaon<0.3) HF1(12316,event.GFlmass);
      // TVector3 lmom(event.GFldecays_mom_x[0]+event.GFldecays_mom_x[1],
      // 		    event.GFldecays_mom_y[0]+event.GFldecays_mom_y[1],
      // 		    event.GFldecays_mom_z[0]+event.GFldecays_mom_z[1]);
      TVector3 lmom(event.decays_mom_x[0]+event.decays_mom_x[1],
		    event.decays_mom_y[0]+event.decays_mom_y[1],
		    event.decays_mom_z[0]+event.decays_mom_z[1]);
      TLorentzVector GFLlmd(lmom, LambdaMass);
      TVector3 p2mom(event.GFmom_x[it][0],event.GFmom_y[it][0],event.GFmom_z[it][0]);
      TLorentzVector GFLp2(p2mom, ProtonMass);
      TLorentzVector GFLlp2 = GFLp2 + GFLlmd;
      double lp2mass = GFLlp2.M();
      if(BEkaon<-0.1) HF1(12322,lp2mass);
      else if(BEkaon<0. ) HF1(12323,lp2mass);
      else if(BEkaon<0.1) HF1(12324,lp2mass);        
      else if(BEkaon<0.2) HF1(12325,lp2mass);
      else if(BEkaon<0.3) HF1(12326,lp2mass);            
    }
  }
    
  if(use_pidlikeli){ // Helix,dEdx Pid
    for (int it1 = 0; it1 < ntTpc; ++it1) { // proton loop
      if(!event.GFfitstatus[it1]) continue;
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if(event.charge[it1] != 1) continue;
      if((event.pid[it1] & 4) != 4) continue; // check proton-like
      for (int it2 = 0; it2 < ntTpc; ++it2) { // pion loop
	if(!event.GFfitstatus[it2]) continue;
	if (it1 == it2) continue;
	std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;
	if (event.charge[it2] != -1) continue;
	std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	
	if ((event.pid[it2] & 1) != 1) continue; // check pion-like
	if( !(std::abs(event.nsigma_proton[it1]) < 2.0 &&
	      std::abs(event.nsigma_pion[it2]) < 2.0) ) continue;
	if(event.lflag){
	  Double_t mom1 = event.GFmom[it1][0];
	  Double_t mom2 = event.GFmom[it2][0];	  
	  HF1(41100, event.lmass);
	  HF2(1100, mom1*event.charge[it1], event.dEdx[it1]);
	  HF2(1100, mom2*event.charge[it2], event.dEdx[it2]);
	  HF2(1101, mom1*event.charge[it1], event.dEdx[it1]);
	  HF2(1102, mom2*event.charge[it2], event.dEdx[it2]);	  	  
	}
	//	HF1(40001,event.lmass);
      }
      // lambda mass
    }
  }
  if(use_pidlikeli){ 
    // Helix,LikelihoodPid
    std::vector<Double_t> L_mass_container_LH;
    std::vector<Int_t> l_p_container_LH, l_pi_container_LH;    
    Int_t l_candidates_LH = 0;
    const double pid_threshold = 0.80;
    for (int it1 = 0; it1 < ntTpc; ++it1) { // proton candidate
      if(!event.GFfitstatus[it1]) continue;
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if(event.charge[it1] != 1) continue;
      //auto prob1 = gPidLike.CalculatePosterior(event.charge[it1], event.mom0[it1], event.GFm2[it1], event.dEdx[it1]);
      //if (event.GFmom[it1].empty()) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__
		<< " GFm2: " << event.GFm2[it1] << std::endl;      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	            
      //if ( event.GFextrapolationHtof[it1]!=1 ) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;
      auto prob1 = gPidLike.CalculatePosterior(event.charge[it1], event.GFmom[it1][0], event.GFm2[it1], event.dEdx[it1]);
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	            
      if(prob1.size()>0){
	std::cout << prob1[kPidP] << std::endl;
      } else { std::cout << "NaN" << std::endl;}
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	   
      if (prob1.empty() || prob1[kPidP] < pid_threshold) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;
      Double_t p_par[5];
      //if (event.helix_t[it1].empty()) continue;
      p_par[0] = event.helix_cx[it1];
      p_par[1] = event.helix_cy[it1];
      p_par[2] = event.helix_z0[it1];
      p_par[3] = event.helix_r[it1];
      p_par[4] = event.helix_dz[it1];
      Int_t p_nh = event.helix_t[it1].size();
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      Double_t p_theta_min = event.helix_t[it1][0] - vtx_scan_range/p_par[3];
      Double_t p_theta_max = TMath::Min(event.helix_t[it1][0] + vtx_scan_rangeInsideL/p_par[3], event.helix_t[it1][p_nh-1]);
      if (event.calpos_x[it1].empty() || event.calpos_y[it1].empty() || event.calpos_z[it1].empty()) continue;
      TVector3 p_start = TVector3(event.calpos_x[it1][0],
				  event.calpos_y[it1][0],
				  event.calpos_z[it1][0]);
      TVector3 p_end = TVector3(event.calpos_x[it1][p_nh-1],
				event.calpos_y[it1][p_nh-1],
				event.calpos_z[it1][p_nh-1]);
      for (int it2 = 0; it2 < ntTpc; ++it2) { // pion candidate
	std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	
	if(!event.GFfitstatus[it2]) continue;      	
	if(event.isElectron[it2]==1) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;	
	if (event.charge[it2] != -1) continue;
	//if ( event.GFextrapolationHtof[it2]!=1 ) continue;	
	if (it1 == it2) continue;	
	//	auto prob2 = gPidLike.CalculatePosterior(event.charge[it2], event.mom0[it2], event.GFm2[it1], event.dEdx[it2]);
	//if (event.GFmom[it2].empty()) continue;
	std::cout << " debug " << __FILE__ << " " << __LINE__
		  << " GFm2: " << event.GFm2[it2] << std::endl;
	auto prob2 = gPidLike.CalculatePosterior(event.charge[it2], event.GFmom[it2][0], event.GFm2[it1], event.dEdx[it2]);
	if (prob2.empty() || prob2[kPidPi] < pid_threshold) continue;
	std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;      		
	Double_t pi_par[5];
	if (event.helix_t[it2].empty()) continue;	
	pi_par[0] = event.helix_cx[it2];
	pi_par[1] = event.helix_cy[it2];
	pi_par[2] = event.helix_z0[it2];
	pi_par[3] = event.helix_r[it2];
	pi_par[4] = event.helix_dz[it2];
	Int_t pi_nh = event.helix_t[it2].size();
	Double_t pi_theta_min = TMath::Max(event.helix_t[it2][0] - vtx_scan_rangeInsideL/pi_par[3], event.helix_t[it2][pi_nh-1]);
	Double_t pi_theta_max = event.helix_t[it2][0] + vtx_scan_range/pi_par[3];
	TVector3 pi_start = TVector3(event.calpos_x[it2][0],
				     event.calpos_y[it2][0],
				     event.calpos_z[it2][0]);
	TVector3 pi_end = TVector3(event.calpos_x[it2][pi_nh-1],
				   event.calpos_y[it2][pi_nh-1],
				   event.calpos_z[it2][pi_nh-1]);
	Double_t ppi_dist = 10000.;
	TVector3 p_mom; TVector3 pi_mom; TVector3 lambda_mom;
	TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField,
							p_par,pi_par,
							p_theta_min,p_theta_max,
							pi_theta_min,pi_theta_max,
							p_mom,pi_mom,lambda_mom,
							ppi_dist);
	if(TMath::IsNaN(ppi_dist)) continue;
	lambda_mom = pi_mom + p_mom;
	TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
	TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));	
	TLorentzVector Llambda = Lp + Lpi;
	if(TMath::Abs(lambda_vert.x()) > 250. ||
	   TMath::Abs(lambda_vert.z()) > 250. ||
	   TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut
	Double_t pi_vertex_dist; Double_t p_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

	if(pi_vertex_dist > pi_vtx_distcut) continue;
	if(p_vertex_dist > p_vtx_distcut) continue;
	if(ppi_dist > ppi_distcut || TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;
	event.lflag = true;
	Double_t ltarget_dist;
	TVector3 ltarget_vtx =
	  Kinematics::CalcCloseDistLambda(tgtpos,
					  lambda_vert,
					  lambda_mom,
					  ltarget_dist);
	L_mass_container_LH.push_back(Llambda.M());
	l_p_container_LH.push_back(it1);
	l_pi_container_LH.push_back(it2);
	l_candidates_LH++;	
	std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;
      }
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;            
    }
    std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;          
    Int_t best_l = -1; Double_t prev_massdiff_l = 9999.;
    for(Int_t candi=0;candi<l_candidates_LH;candi++){
      Double_t diff = TMath::Abs(L_mass_container_LH[candi] - LambdaMass);
      if(prev_massdiff_l > diff){
	prev_massdiff_l = diff;
	best_l = candi;
	std::cout << "best lambda: " << best_l << std::endl;
      }
    }
    std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;
    if(best_l!=-1){
      double bestmass = L_mass_container_LH[best_l];
      int it1 = l_p_container_LH[best_l];
      int it2 = l_pi_container_LH[best_l];
      Double_t mom1 = event.GFmom[it1][0];
      Double_t mom2 = event.GFmom[it2][0];
      HF1(41300, bestmass);
      HF2(1300, mom1*event.charge[it1], event.dEdx[it1]);
      HF2(1300, mom2*event.charge[it2], event.dEdx[it1]);
      HF2(1301, mom1*event.charge[it1], event.dEdx[it1]);
      HF2(1302, mom2*event.charge[it2], event.dEdx[it2]);	  	        
    }
  }
  if(use_pidlikeli){ //Genfit, dEdxPid
    const auto& allDecays = **src.decaysidLambda;  // std::vector<std::vector<Double_t>>&
    if(allDecays.size() >= 1 && allDecays[0].size() >= 2) {
      const auto& decays = allDecays[0];          // std::vector<Double_t>&
      Int_t pId  = static_cast<Int_t>(decays[0]);
      Int_t piId = static_cast<Int_t>(decays[1]);
      if(event.GFfitstatus[pId]&&event.GFfitstatus[piId]){   		      
	if(std::abs(event.nsigma_proton[pId]) < 2.0 &&
	   std::abs(event.nsigma_pion [piId]) < 2.0) {
	  Double_t lamM = **src.GFlmass;
	  if (!event.GFmom[pId].empty() && 
	      !event.GFmom[piId].empty() ){
	    Double_t mom1 = event.GFmom[pId][0];
	    Double_t mom2 = event.GFmom[piId][0];	
	    HF1(42100, lamM);
	    HF2(2100, mom1*event.charge[pId], event.dEdx[pId]);
	    HF2(2100, mom2*event.charge[piId], event.dEdx[piId]);
	    HF2(2101, mom1*event.charge[pId], event.dEdx[pId]);
	    HF2(2102, mom2*event.charge[piId], event.dEdx[piId]);
	  }
	}
      }
    }
  }
  if(use_pidlikeli){ //Genfit, LikelihoodPid
    const double pid_thresh = 0.80;
    const auto& allDecays = **src.decaysidLambda;   // std::vector<std::vector<Double_t>>&
    if(allDecays.size() >= 1 && allDecays[0].size() >= 2) {
      const auto& decays = allDecays[0];            
      Int_t pId  = static_cast<Int_t>(decays[0]);
      Int_t piId = static_cast<Int_t>(decays[1]);
      if ( !event.GFmom[pId].empty() && 
	   !event.GFmom[piId].empty() ) {
	Double_t mom1 = event.GFmom[pId][0];
	Double_t mom2 = event.GFmom[piId][0];	      
	auto prob_p  = gPidLike.CalculatePosterior(
						   event.charge[pId],
						   mom1,
						   event.GFm2[  pId],
						   event.dEdx[  pId]
						   );
	auto prob_pi = gPidLike.CalculatePosterior(
						   event.charge[piId],
						   mom2,
						   event.GFm2[  piId],
						   event.dEdx[  piId]
						   );	
	if(!prob_p.empty() && !prob_pi.empty() &&
	   prob_p[ kPidP ] >= pid_thresh &&
	   prob_pi[kPidPi] >= pid_thresh) {
	  Double_t lamM = **src.GFlmass; 
	  HF1(42300, lamM); 
	  HF2(2300, mom1*event.charge[pId], event.dEdx[pId]);   
	  HF2(2300, mom2*event.charge[piId], event.dEdx[piId]); 
	  HF2(2301, mom1*event.charge[pId], event.dEdx[pId]);   
	  HF2(2302, mom2*event.charge[piId], event.dEdx[piId]);  
	}
      }
    }
  }
  std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;

  if(KKEvent){
    HF1(50001, MissMass);
    if(event.xiflag) HF1(40001, event.GFximass);
  }
  //std::vector<Double_t> GFL_mass_container(l_candidates_LH, qnan);
    
  HF1( 1, event.status++ );

  //  GFTrackCont.Clear();
  
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;  
  TFileCont[kOutFile]->Close();
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;  
  const Int_t n = TFileCont.size();
  for( Int_t i=0; i<n; ++i ){
    if( TTreeReaderCont[i] ) delete TTreeReaderCont[i];
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;      
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{

  static const auto KKEvent = gUser.GetParameter("KKEvent");
  static const auto KPEvent = gUser.GetParameter("KPEvent");
  static const auto KHeavyEvent = gUser.GetParameter("KHeavyEvent");
  
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << __func__ << std::endl;;
  Int_t nbinpoq = 1000;
  Int_t minpoq = -1.5;
  Int_t maxpoq = 1.5;
  Int_t nbindedx = 1000;
  Int_t mindedx = 0;
  Int_t maxdedx = 350;
  Int_t nbinmass2 = 500;
  Int_t minmass2 = -5;
  Int_t maxmass2 = 5;
  Int_t nbininvbeta = 100;
  Int_t mininvbeta = 0;
  Int_t maxinvbeta = 5;
  double th95 = 0.95;
  double th90 = 0.90;
  double th85 = 0.85;
  
  HB1( 1, "Status", 21, 0., 21. );
  HB1( 2, "Genfit Status", 20, 0., 20. );
  HB1( 3, "Genfit Fit Status", 2, 0., 2. );
  //HB1( 10, "NTrack TPC", 20, 0., 20. );

  int nLmd=2; int nK0=2; int nGamma=3; int nBE = 3;
  int nbin = static_cast<int>(pow(BASE_TRACK, (int)Species::N_SPECIES));
  for(int iBE=0; iBE<nBE; iBE++){  
    for(int iGamma=0; iGamma<nGamma; iGamma++){
      for(int iK0=0; iK0<nK0; iK0++){	
	for(int iLmd=0; iLmd<nLmd; iLmd++){
	  HB1( (iLmd+iK0*nLmd+iGamma*nLmd*nK0+1)+iBE*100+100,
	       Form("PID decomposition BE%d, #Lambda:%d,K0:%d,#gamma:%d",iBE+1,iLmd,iK0,iGamma),
	       nbin, 0, nbin );
	}
      }
    }
  }
  for(int iBE=3; iBE<=5; iBE++){ // BE region 1A/1B/1C
    for(int iGamma=0; iGamma<nGamma; iGamma++){
      for(int iK0=0; iK0<nK0; iK0++){	
	for(int iLmd=0; iLmd<nLmd; iLmd++){
	  if(iBE==3){
	    HB1( (iLmd+iK0*nLmd+iGamma*nLmd*nK0+1)+iBE*100+100,
		 Form("PID decomposition BE1a, #Lambda:%d,K0:%d,#gamma:%d",iLmd,iK0,iGamma),
		 nbin, 0, nbin );
	  } else if(iBE==4){
	    HB1( (iLmd+iK0*nLmd+iGamma*nLmd*nK0+1)+iBE*100+100,
		 Form("PID decomposition BE1b, #Lambda:%d,K0:%d,#gamma:%d",iLmd,iK0,iGamma),
		 nbin, 0, nbin );	    
	  } else {
	    HB1( (iLmd+iK0*nLmd+iGamma*nLmd*nK0+1)+iBE*100+100,
		 Form("PID decomposition BE1c, #Lambda:%d,K0:%d,#gamma:%d",iLmd,iK0,iGamma),
		 nbin, 0, nbin );	    
	  }
	}
      }
    }
  }
  //HB1( 100, "Decomposition TPC track", 3000000, 0., 3000000 );  
  //HB1( 100, "Decomposition TPC track", 3000000, 0., 3000000 );
  //HB1( 100, "Decomposition TPC track", 3*W_LMD, 0., 3*W_LMD );  

  // no cut
  HB2( 1001, "<dE/dx> ; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  // dEdx Helix,dEdxPid
  HB2( 1100, Form("<dE/dx> #Lambda [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1101, Form("<dE/dx> p+_{#Lambda} [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1102, Form("<dE/dx> #pi-_{#Lambda} [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1110, Form("<dE/dx> #Lambda [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1111, Form("<dE/dx> p+_{#Lambda} [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1112, Form("<dE/dx> #pi-_{#Lambda} [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1120, Form("<dE/dx> #Lambda [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1121, Form("<dE/dx> p+_{#Lambda} [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1122, Form("<dE/dx> #pi-_{#Lambda} [Helix,dEdxPid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  // dEdx Helix,M2Pid
  HB2( 1200, Form("<dE/dx> #Lambda [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1201, Form("<dE/dx> p+_{#Lambda} [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1202, Form("<dE/dx> #pi-_{#Lambda} [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1210, Form("<dE/dx> #Lambda [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1211, Form("<dE/dx> p+_{#Lambda} [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1212, Form("<dE/dx> #pi-_{#Lambda} [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1220, Form("<dE/dx> #Lambda [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1221, Form("<dE/dx> p+_{#Lambda} [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1222, Form("<dE/dx> #pi-_{#Lambda} [Helix,M2Pid][%.2f]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  // dEdx Helix,LikelihoodPid
  HB2( 1300, Form("<dE/dx> #Lambda [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1301, Form("<dE/dx> p+_{#Lambda} [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1302, Form("<dE/dx> #pi-_{#Lambda} [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1310, Form("<dE/dx> #Lambda [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1311, Form("<dE/dx> p+_{#Lambda} [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1312, Form("<dE/dx> #pi-_{#Lambda} [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1320, Form("<dE/dx> #Lambda [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 1321, Form("<dE/dx> p+_{#Lambda} [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 1322, Form("<dE/dx> #pi-_{#Lambda} [Helix,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  // dEdx Genfit,dEdxPid
  HB2( 2100, Form("<dE/dx> #Lambda [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2101, Form("<dE/dx> p+_{#Lambda} [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2102, Form("<dE/dx> #pi-_{#Lambda} [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2110, Form("<dE/dx> #Lambda [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2111, Form("<dE/dx> p+_{#Lambda} [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2112, Form("<dE/dx> #pi-_{#Lambda} [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2120, Form("<dE/dx> #Lambda [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2121, Form("<dE/dx> p+_{#Lambda} [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2122, Form("<dE/dx> #pi-_{#Lambda} [Genfit,dEdxPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);    
  // dEdx Genfit,M2Pid
  HB2( 2200, Form("<dE/dx> #Lambda [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95),nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2201, Form("<dE/dx> p+_{#Lambda} [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2202, Form("<dE/dx> #pi-_{#Lambda} [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2210, Form("<dE/dx> #Lambda [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2211, Form("<dE/dx> p+_{#Lambda} [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2212, Form("<dE/dx> #pi-_{#Lambda} [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2220, Form("<dE/dx> #Lambda [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2221, Form("<dE/dx> p+_{#Lambda} [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2222, Form("<dE/dx> #pi-_{#Lambda} [Genfit,M2Pid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);      
  // dEdx Helix,LikelihoodPid
  HB2( 2300, Form("<dE/dx> #Lambda [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2301, Form("<dE/dx> p+_{#Lambda} [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2302, Form("<dE/dx> #pi-_{#Lambda} [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th95), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2310, Form("<dE/dx> #Lambda [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2311, Form("<dE/dx> p+_{#Lambda} [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2312, Form("<dE/dx> #pi-_{#Lambda} [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th90), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2320, Form("<dE/dx> #Lambda [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 2321, Form("<dE/dx> p+_{#Lambda} [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 2322, Form("<dE/dx> #pi-_{#Lambda} [Genfit,LikeliPid]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]",th85), nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);

  HB2( 3100, "<dE/dx> {#delta M_{#Lambda}<0.1} [p#pi]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 3101, "<dE/dx> {#delta M_{#Lambda}<0.1} [p]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);
  HB2( 3102, "<dE/dx> {#delta M_{#Lambda}<0.1} [#pi]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 3110, "<dE/dx> {#delta M_{#Lambda}>0.1} [p#pi]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 3111, "<dE/dx> {#delta M_{#Lambda}>0.1} [p]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  
  HB2( 3112, "<dE/dx> {#delta M_{#Lambda}>0.1} [#pi]; p/q [GeV/#font[12]{c}];<dE/dx> [arb.]", nbinpoq, minpoq, maxpoq, nbindedx, mindedx, maxdedx);  

  HB1( 11201, "[Genfit] Binding Energy of Kaon ; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11202, "[Genfit] Binding Energy of Kaon [3.5#circ<thetaTPC<4.5#circ]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11203, "[Genfit] Binding Energy of Kaon [thetaTPC<10#circ]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11204, "[Genfit] Binding Energy of Kaon [w/ #Lambda]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11205, "[Genfit] Binding Energy of Kaon [3.5#circ<thetaTPC<4.5#circ && w/ #Lambda]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11206, "[Genfit] Binding Energy of Kaon [thetaTPC<10#circ && w/ #Lambda]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11207, "[Genfit] Binding Energy of Kaon [w/ #Lambdap]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11208, "[Genfit] Binding Energy of Kaon [3.5#circ<thetaTPC<4.5#circ && w/ #Lambdap]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11209, "[Genfit] Binding Energy of Kaon [thetaTPC<10#circ && w/ #Lambdap]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);

  HB1( 11301, "[Genfit] Binding Energy of Kaon [m2]; -B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11302, "[Genfit] Binding Energy of Kaon [m2][3.5#circ<thetaTPC<4.5#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11303, "[Genfit] Binding Energy of Kaon [m2][thetaTPC<10#circ]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11304, "[Genfit] Binding Energy of Kaon [m2][w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11305, "[Genfit] Binding Energy of Kaon [m2][3.5#circ<thetaTPC<4.5#circ && w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11306, "[Genfit] Binding Energy of Kaon [m2][thetaTPC<10#circ && w/ #Lambda]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11307, "[Genfit] Binding Energy of Kaon [m2][w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11308, "[Genfit] Binding Energy of Kaon [m2][3.5#circ<thetaTPC<4.5#circ && w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  HB1( 11309, "[Genfit] Binding Energy of Kaon [m2][thetaTPC<10#circ && w/ #Lambdap]; B_{K} [GeV]; Counts [/5 MeV]", 200, -0.5, 0.5);
  
  HB1( 12001, "[Genfit] #Lambda Invariant Mass ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12002, "[Genfit] #Lambda Invariant Mass [-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12003, "[Genfit] #Lambda Invariant Mass [-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12004, "[Genfit] #Lambda Invariant Mass [0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12005, "[Genfit] #Lambda Invariant Mass [0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12006, "[Genfit] #Lambda Invariant Mass [0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12011, "[Genfit] #Lambda Invariant Mass [one more p]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12012, "[Genfit] #Lambda Invariant Mass [one more p][-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12013, "[Genfit] #Lambda Invariant Mass [one more p][-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12014, "[Genfit] #Lambda Invariant Mass [one more p][0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12015, "[Genfit] #Lambda Invariant Mass [one more p][0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12016, "[Genfit] #Lambda Invariant Mass [one more p][0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);    
  HB1( 12021, "[Genfit] #LambdaP Invariant Mass ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12022, "[Genfit] #LambdaP Invariant Mass [-BEk<-0.1(GeV)]   ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12023, "[Genfit] #LambdaP Invariant Mass [-0.1<-BEk<0(GeV)] ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12024, "[Genfit] #LambdaP Invariant Mass [0<-BEk<0.1(GeV)]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12025, "[Genfit] #LambdaP Invariant Mass [0.1<-BEk<0.2(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12026, "[Genfit] #LambdaP Invariant Mass [0.2<-BEk<0.3(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);

  HB1( 12301, "[Genfit] #Lambda Invariant Mass [m2]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12302, "[Genfit] #Lambda Invariant Mass [m2][-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12303, "[Genfit] #Lambda Invariant Mass [m2][-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12304, "[Genfit] #Lambda Invariant Mass [m2][0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12305, "[Genfit] #Lambda Invariant Mass [m2][0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12306, "[Genfit] #Lambda Invariant Mass [m2][0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12311, "[Genfit] #Lambda Invariant Mass [m2][one more p]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12312, "[Genfit] #Lambda Invariant Mass [m2][one more p][-BEk<-0.1(GeV)]   ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12313, "[Genfit] #Lambda Invariant Mass [m2][one more p][-0.1<-BEk<0(GeV)] ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12314, "[Genfit] #Lambda Invariant Mass [m2][one more p][0<-BEk<0.1(GeV)]  ; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12315, "[Genfit] #Lambda Invariant Mass [m2][one more p][0.1<-BEk<0.2(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);
  HB1( 12316, "[Genfit] #Lambda Invariant Mass [m2][one more p][0.2<-BEk<0.3(GeV)]; #Lambda IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 180, 1.04, 1.4);    
  HB1( 12321, "[Genfit] #LambdaP Invariant Mass [m2]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12322, "[Genfit] #LambdaP Invariant Mass [m2][-BEk<-0.1(GeV)]   ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12323, "[Genfit] #LambdaP Invariant Mass [m2][-0.1<-BEk<0(GeV)] ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12324, "[Genfit] #LambdaP Invariant Mass [m2][0<-BEk<0.1(GeV)]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12325, "[Genfit] #LambdaP Invariant Mass [m2][0.1<-BEk<0.2(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 12326, "[Genfit] #LambdaP Invariant Mass [m2][0.2<-BEk<0.3(GeV)]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  
  HB1( 13240, "[Genfit] #LambdaP Invariant Mass ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13241, "[Genfit] #LambdaP Invariant Mass [0#circ<LPOpAngle<20#circ]   ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13242, "[Genfit] #LambdaP Invariant Mass [20#circ<LPOpAngle<40#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13243, "[Genfit] #LambdaP Invariant Mass [40#circ<LPOpAngle<60#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13244, "[Genfit] #LambdaP Invariant Mass [60#circ<LPOpAngle<80#circ]  ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13245, "[Genfit] #LambdaP Invariant Mass [80#circ<LPOpAngle<100#circ] ; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13246, "[Genfit] #LambdaP Invariant Mass [100#circ<LPOpAngle<120#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13247, "[Genfit] #LambdaP Invariant Mass [120#circ<LPOpAngle<140#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13248, "[Genfit] #LambdaP Invariant Mass [140#circ<LPOpAngle<160#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  HB1( 13249, "[Genfit] #LambdaP Invariant Mass [160#circ<LPOpAngle<180#circ]; #LambdaP IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 1750, 2.0, 5.5);
  
  //  HB1(10000, )
  HB1( 40001, "[Genfit] #Xi Invariant Mass ; #Xi IM [GeV]; Counts [/0.002 GeV/#font[12]{c}^{2}]", 380, 1.04, 1.8);  
  HB1( 50001, "MissMass KK on H(CH2); MissMass [GeV]; Counts [/5 MeV]", 600, 0., 3.0);

  // int hid_exc = 100000; //Lmd:0,K0:0,gamma:0
  // int hid_exc_Lmd = 200000;  //Lmd:1,K0:0,gamma:0
  // int hid_exc_gamma = 500000; //Lmd:0,K0:0,gamma:1
  HB1( hid_inc, "-BEk inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 90000
  HB1( hid_exc+base5_to_decimal("0"), "-BEk empty; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000
  HB1( hid_exc+base5_to_decimal("1"), "-BEk 1#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+1
  HB1( hid_exc+base5_to_decimal("2"), "-BEk 2#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+2
  HB1( hid_exc+base5_to_decimal("10"), "-BEk 1#pi+; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+5
  HB1( hid_exc+base5_to_decimal("11"), "-BEk 1#pi+1#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+6
  HB1( hid_exc+base5_to_decimal("100"), "-BEk 1K-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+25
  HB1( hid_exc+base5_to_decimal("101"), "-BEk 1K-1#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+26
  HB1( hid_exc+base5_to_decimal("102"), "-BEk 1K-2#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+27
  HB1( hid_exc+base5_to_decimal("10000"), "-BEk 1p; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 1000000+125
  HB1( hid_exc+base5_to_decimal("10001"), "-BEk 1p1#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 1000000+126
  HB1( hid_exc+base5_to_decimal("10002"), "-BEk 1p2#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000+127      
  HB1( hid_exc_Lmd+base5_to_decimal("0"), "-BEk #Lambda; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300); // 2000000+0
  HB1( hid_exc_gamma+base5_to_decimal("11"), "-BEk #gamma1#pi+1#pi-; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 500000+6
  HB1( hid_exc+10000, "-BEk othres; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 100000 + 10000 
  HB1( hid_semiexc+base5_to_decimal("1"), "-BEk 1#pi- inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 800000+1 
  HB1( hid_semiexc+base5_to_decimal("10"), "-BEk 1#pi+ inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 800000+5 
  HB1( hid_semiexc+base5_to_decimal("100"), "-BEk 1K- inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 800000+25 
  HB1( hid_semiexc+base5_to_decimal("10000"), "-BEk 1p inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 800000+25 
  HB1( hid_semiexc+20000  , " -BEk #Lambda inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 800000+20000+1 
  HB1( hid_semiexc+10000  , " -BEk 1nega inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 800000+10000+0 
  HB1( hid_semiexc+10000+1, " -BEk 1posi inclusive; -BE[MeV]; Counts [/0.002 GeV]", 120,-300,300);  // 800000+10000+0 
  
  HBTree( "tpc", "tree of GenfitKNucleusAna" );
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
  tree->Branch( "remain_nclTpc", &event.remain_nclTpc );
  
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

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "isInTarget", &event.isInTarget );
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
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "pval", &event.pval );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "mom_x", &event.mom_x );
  tree->Branch( "mom_y", &event.mom_y );
  tree->Branch( "mom_z", &event.mom_z );
  
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
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);

  tree->Branch( "helix_t", &event.helix_t );
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

  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFtof", &event.GFtof);
  //tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFfitstatus", &event.GFfitstatus);
  tree->Branch("GFpdgcode", &event.GFpdgcode);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFlayer", &event.GFlayer);
  tree->Branch("GFpos_x", &event.GFpos_x);
  tree->Branch("GFpos_y", &event.GFpos_y);
  tree->Branch("GFpos_z", &event.GFpos_z);
  tree->Branch("GFmom", &event.GFmom);
  tree->Branch("GFmom_x", &event.GFmom_x);
  tree->Branch("GFmom_y", &event.GFmom_y);
  tree->Branch("GFmom_z", &event.GFmom_z);
  tree->Branch("GFresidual_x", &event.GFresidual_x);
  tree->Branch("GFresidual_y", &event.GFresidual_y);
  tree->Branch("GFresidual_z", &event.GFresidual_z);
  tree->Branch("GFresidual_p", &event.GFresidual_p);
  tree->Branch("GFresidual_px", &event.GFresidual_px);
  tree->Branch("GFresidual_py", &event.GFresidual_py);
  tree->Branch("GFresidual_pz", &event.GFresidual_pz);

  tree->Branch("GFntTpc_target", &event.GFntTpc_target);
  tree->Branch("GFprodvtx_x", &event.GFprodvtx_x);
  tree->Branch("GFprodvtx_y", &event.GFprodvtx_y);
  tree->Branch("GFprodvtx_z", &event.GFprodvtx_z);

  //extrapolation
  tree->Branch("GFinside", &event.GFinside);
  tree->Branch("GFinsideTgtToKurama", &event.GFKuramaFromTgt);
  tree->Branch("GFKuramaHasVtxOutOfTgt", &event.GFKuramaVtxOutTgt);  
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFtrack2vtxdist", &event.GFtrack2vtxdist);
  tree->Branch("GFcalctof", &event.GFcalctof);
  tree->Branch("GFsegHtof", &event.GFsegHtof);
  tree->Branch("GFtofHtof", &event.GFtofHtof);
  tree->Branch("GFtdiffHtof", &event.GFtdiffHtof);
  tree->Branch("GFposHtof", &event.GFposHtof);
  tree->Branch("GFposx", &event.GFposx);
  tree->Branch("GFposy", &event.GFposy);
  tree->Branch("GFposz", &event.GFposz);
  tree->Branch("GFinvbeta", &event.GFinvbeta);
  tree->Branch("GFm2", &event.GFm2);
  tree->Branch("nsigma_tritonHtof", &event.nsigma_tritonHtof);
  tree->Branch("nsigma_deutronHtof", &event.nsigma_deutronHtof);
  tree->Branch("nsigma_protonHtof", &event.nsigma_protonHtof);
  tree->Branch("nsigma_kaonHtof", &event.nsigma_kaonHtof);
  tree->Branch("nsigma_pionHtof", &event.nsigma_pionHtof);
  tree->Branch("nsigma_electronHtof", &event.nsigma_electronHtof);

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

  tree->Branch( "ncombiReconFailedLambda", &event.ncombiLreconfailed  );
  tree->Branch( "ReconFailedLambdaPId", &event.pidLreconfailed);
  tree->Branch( "ReconFailedLambdaPiId", &event.piidLreconfailed);
  tree->Branch( "ReconFailedLambdaMass", &event.LmassLreconfailed);
  tree->Branch( "ReconFailedLambdaDecayVtx_x", &event.LdecayvtxLreconfailed_x);
  tree->Branch( "ReconFailedLambdaDecayVtx_y", &event.LdecayvtxLreconfailed_y);
  tree->Branch( "ReconFailedLambdaDecayVtx_z", &event.LdecayvtxLreconfailed_z);
  tree->Branch( "ReconFailedLambdaMom", &event.LmomLreconfailed );
  tree->Branch( "ReconFailedLambdaMom_x", &event.LmomLreconfailed_x );
  tree->Branch( "ReconFailedLambdaMom_y", &event.LmomLreconfailed_y );
  tree->Branch( "ReconFailedLambdaMom_z", &event.LmomLreconfailed_z );
  tree->Branch( "ReconFailedLambdaPMom", &event.pmomLreconfailed );
  tree->Branch( "ReconFailedLambdaPMom_x", &event.pmomLreconfailed_x );
  tree->Branch( "ReconFailedLambdaPMom_y", &event.pmomLreconfailed_y );
  tree->Branch( "ReconFailedLambdaPMom_z", &event.pmomLreconfailed_z );
  tree->Branch( "ReconFailedLambdaPiMom", &event.pimomLreconfailed );
  tree->Branch( "ReconFailedLambdaPiMom_x", &event.pimomLreconfailed_x );
  tree->Branch( "ReconFailedLambdaPiMom_y", &event.pimomLreconfailed_y );
  tree->Branch( "ReconFailedLambdaPiMom_z", &event.pimomLreconfailed_z );
  tree->Branch( "ReconFailedLambdavtxCloseDist", & event.ppidistLreconfailed);

  tree->Branch( "ncombiPiPair", &event.ncombiPipair);
  tree->Branch( "PiPairPipId", &event.pipidPipair);
  tree->Branch( "PiPairPimId", &event.pimidPipair);
  tree->Branch( "PiPairPipMom", &event.pipmomPipair);
  tree->Branch( "PiPairPipMom_x", &event.pipmomPipair_x);
  tree->Branch( "PiPairPipMom_y", &event.pipmomPipair_y);
  tree->Branch( "PiPairPipMom_z", &event.pipmomPipair_z);
  tree->Branch( "PiPairPimMom", &event.pimmomPipair);
  tree->Branch( "PiPairPimMom_x", &event.pimmomPipair_x);
  tree->Branch( "PiPairPimMom_y", &event.pimmomPipair_y);
  tree->Branch( "PiPairPimMom_z", &event.pimmomPipair_z);
  tree->Branch( "PiPairMom", & event.momPipair);
  tree->Branch( "PiPairMom_x", & event.momPipair_x);
  tree->Branch( "PiPairMom_y", & event.momPipair_y);
  tree->Branch( "PiPairMom_z", & event.momPipair_z);
  tree->Branch( "PiPairReconLambdaMass", & event.reconLmassPipair);
  tree->Branch( "PiPairReconMass", & event.reconmassPipair);
  tree->Branch( "PiPairCloseDist", & event.pipidistPipair);

  tree->Branch( "nEscapeKm", &event.nEscapeKm);
  tree->Branch( "KmId", &event.kmid);
  tree->Branch( "KmMom", &event.kmmom);
  tree->Branch( "KmMom_x", &event.kmmom_x);
  tree->Branch( "KmMom_y", &event.kmmom_y);
  tree->Branch( "KmMom_z", &event.kmmom_z);

  tree->Branch( "GFKmDecayVtx_x", &event.GFkmdecayvtx_x);
  tree->Branch( "GFKmDecayVtx_y", &event.GFkmdecayvtx_y);
  tree->Branch( "GFKmDecayVtx_z", &event.GFkmdecayvtx_z);
  tree->Branch( "GFKmMom", &event.GFkmmom);
  tree->Branch( "GFKmMom_x", &event.GFkmmom_x);
  tree->Branch( "GFKmMom_y", &event.GFkmmom_y);
  tree->Branch( "GFKmMom_z", &event.GFkmmom_z);
  //tree->Branch( "GFKmVtxCloseDist", &event.GFppi_dist);
  tree->Branch( "GFKmTargetCloseDist", &event.GFkmtarget_dist);
  tree->Branch( "GFKmTarget_x", &event.GFkmtargetvtx_x);
  tree->Branch( "GFKmTarget_y", &event.GFkmtargetvtx_y);
  tree->Branch( "GFKmTarget_z", &event.GFkmtargetvtx_z);
  tree->Branch( "GFKmTargetCenter_x", &event.GFkmtargetcenter_x);
  tree->Branch( "GFKmTargetCenter_y", &event.GFkmtargetcenter_y);
  tree->Branch( "GFKmTargetCenter_z", &event.GFkmtargetcenter_z);
  tree->Branch( "GFKmTargetCenterCloseDist", &event.GFkmtargetcenter_dist);
  tree->Branch( "GFKmTrackLen", &event.GFkmtracklen);
  tree->Branch( "GFKmMassSquare", &event.GFkmtracklen);  
  tree->Branch( "GFKmTof", &event.GFkmtof);  
  
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
  tree->Branch( "m2Kurama", &event.m2Kurama);
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
  tree->Branch( "thetaTPC", &event.thetaTPC);  
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);

  tree->Branch("nKK", &event.nKK);
  tree->Branch("Kflag", &event.Kflag);
  tree->Branch("Pflag", &event.Pflag);  
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
  tree->Branch("BETPC", &event.BETPC);
  tree->Branch("BE_LL", &event.BE_LL);
  tree->Branch("BETPC_LL", &event.BETPC_LL);

  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom", &event.lmom);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  
  tree->Branch("Lflag", &event.lflag);
  tree->Branch("LambdaTargetCloseDist", &event.ltarget_dist);
  tree->Branch("LambdaTargetCloseVtx_x", &event.ltargetvtx_x);
  tree->Branch("LambdaTargetCloseVtx_y", &event.ltargetvtx_y);
  tree->Branch("LambdaTargetCloseVtx_z", &event.ltargetvtx_z);

  tree->Branch("GFLambdaMass", &event.GFlmass);
  // tree->Branch("GFLambdaDecayVtx_x", &event.GFdecayvtx_x);
  // tree->Branch("GFLambdaDecayVtx_y", &event.GFdecayvtx_y);
  // tree->Branch("GFLambdaDecayVtx_z", &event.GFdecayvtx_z);
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
  tree->Branch("GFLambdaCTau", &event.GFlctau);    

  tree->Branch("LToKuramaPflag", &event.kuramalflag);  
  tree->Branch("LPiflag", &event.lpiflag);
  tree->Branch("LPflag", &event.lpflag);  
  tree->Branch("PiPiflag", &event.pipiflag);
  tree->Branch("Pimflag", &event.pimflag);
  tree->Branch("Emptyflag", &event.emptyflag);

  //for decay particles
  tree->Branch("GFDecaysNhit", &event.GFdecays_nhtrack);
  tree->Branch("GFDecaysChisqr", &event.GFdecays_chisqr);
  tree->Branch("GFDecaysCharge", &event.GFdecays_charge);
  tree->Branch("GFDecaysTof", &event.GFdecays_tof);
  tree->Branch("GFDecaysPval", &event.GFdecays_pval);
  tree->Branch("GFDecaysPdgcode", &event.GFdecays_pdgcode);
  tree->Branch("GFDecaysHtofId", &event.GFdecays_htofid);
  tree->Branch("GFDecaysTrackLen", &event.GFdecays_tracklen);
  tree->Branch("GFDecaysTrackTof", &event.GFdecays_tof);
  tree->Branch("GFDecaysMassSquare", &event.GFdecays_mass2);
  tree->Branch("GFDecaysInvbeta", &event.GFdecays_invbeta);
  tree->Branch("GFDecaysMom", &event.GFdecays_mom);
  tree->Branch("GFDecaysMom_x", &event.GFdecays_mom_x);
  tree->Branch("GFDecaysMom_y", &event.GFdecays_mom_y);
  tree->Branch("GFDecaysMom_z", &event.GFdecays_mom_z);
  tree->Branch("GFDecaysMomCM", &event.GFdecays_CMmom);
  tree->Branch("GFDecaysMomCM_x", &event.GFdecays_CMmom_x);
  tree->Branch("GFDecaysMomCM_y", &event.GFdecays_CMmom_y);
  tree->Branch("GFDecaysMomCM_z", &event.GFdecays_CMmom_z);
  tree->Branch("GFDecaysMomLoss", &event.GFdecays_momloss);
  tree->Branch("GFDecaysELoss", &event.GFdecays_eloss);

  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);
  tree->Branch("DecaysMomCM", &event.decays_CMmom);
  tree->Branch("DecaysMomCM_x", &event.decays_CMmom_x);
  tree->Branch("DecaysMomCM_y", &event.decays_CMmom_y);
  tree->Branch("DecaysMomCM_z", &event.decays_CMmom_z);

  //Remaining p, pi after Xi, L searching
  //Multiplicity means tracks comes from the target
  tree->Branch("AccidentalMultiplicity", &event.accident_multi);
  tree->Branch("AccidentalTrackId", &event.accident_id);

  tree->Branch("ResidualsMultiplicity", &event.residual_multi);
  tree->Branch("pMultiplicity", &event.p_multi);
  tree->Branch("pipMultiplicity", &event.pip_multi);
  tree->Branch("pimMultiplicity", &event.pim_multi);
  tree->Branch("epMultiplicity", &event.ep_multi);
  tree->Branch("emMultiplicity", &event.em_multi);
  tree->Branch("ppipMultiplicity", &event.ppip_multi);
  tree->Branch("ResidualsTrackId", &event.residual_id);
  tree->Branch("ResidualsMassSquare", &event.residual_mass2);
  tree->Branch("ResidualsInvbeta", &event.residual_invbeta);
  tree->Branch("ResidualsCloseDistTgt", &event.residual_dist2tgt);
  tree->Branch("ResidualsGFProductionVtxCloseDist", &event.residual_GFdist2prodvtx);
  tree->Branch("ResidualsKFProductionVtxCloseDist", &event.residual_KFdist2prodvtx);
  tree->Branch("ResidualsMom", &event.residual_mom);
  tree->Branch("ResidualsMom_x", &event.residual_mom_x);
  tree->Branch("ResidualsMom_y", &event.residual_mom_y);
  tree->Branch("ResidualsMom_z", &event.residual_mom_z);
  tree->Branch("ResidualsCharge", &event.residual_charge);

  tree->Branch("nGamma", &event.g_multi);
  tree->Branch("GammaEpTrackId", &event.epidgamma);
  tree->Branch("GammaEmTrackId", &event.emidgamma);
  tree->Branch("GammaMomId", &event.epmomgamma);
  tree->Branch("GammaEpMom_x", &event.epmomgamma_x);
  tree->Branch("GammaEpMom_y", &event.epmomgamma_y);
  tree->Branch("GammaEpMom_z", &event.epmomgamma_z);
  tree->Branch("GammaEmMom", &event.emmomgamma);
  tree->Branch("GammaEmMom_x", &event.emmomgamma_x);
  tree->Branch("GammaEmMom_y", &event.emmomgamma_y);
  tree->Branch("GammaEmMom_z", &event.emmomgamma_z);
  tree->Branch("GammaMom", &event.momgamma);
  tree->Branch("GammaMom_x", &event.momgamma_x);
  tree->Branch("GammaMom_y", &event.momgamma_y);
  tree->Branch("GammaMom_z", &event.momgamma_z);
  tree->Branch("GammaVtxCloseDist", &event.epidistgamma);
  tree->Branch("GammaDecayVtx_x", &event.vtxgamma_x);
  tree->Branch("GammaDecayVtx_y", &event.vtxgamma_y);
  tree->Branch("GammaDecayVtx_z", &event.vtxgamma_z);

  
  TTreeReaderCont[kE42] = new TTreeReader( "tpc", TFileCont[kE42] );
  const auto& reader = TTreeReaderCont[kE42];

  //src.status = new TTreeReaderValue<Int_t>(*reader, "status" );
  src.runnum = new TTreeReaderValue<Int_t>(*reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>(*reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trigflag" );

  src.nhHtof = new TTreeReaderValue<Int_t>(*reader, "nhHtof" );
  src.HtofSeg = new TTreeReaderValue<std::vector<Double_t>>(*reader, "HtofSeg" );
  src.tHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "tHtof" );
  src.dtHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "dtHtof" );
  src.deHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "deHtof" );
  src.posHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "posHtof" );

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
  //src.m2Org = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2Org" );
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
  // src.xout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xout" );
  // src.yout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yout" );
  // src.zout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "zout" );
  // src.pxout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pxout" );
  // src.pyout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pyout" );
  // src.pzout = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pzout" );

  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
  src.kflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "kflagTPCKurama" );
  //src.pflagTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "pflagTPCKurama" );
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

  // src.nKm = new TTreeReaderValue<Int_t>( *reader, "nKm" );
  // src.nKp = new TTreeReaderValue<Int_t>( *reader, "nKp" );
  src.nKK = new TTreeReaderValue<Int_t>( *reader, "nKK" );
  //src.inside = new TTreeReaderValue<std::vector<Int_t>>( *reader, "inside" );
  src.vtx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtx" );
  src.vty = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vty" );
  src.vtz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtz" );
  //src.closeDist = new TTreeReaderValue<std::vector<Double_t>>( *reader, "closeDist" );
  src.MissMass = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMass" );
  src.MissMassCorr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorr" );
  src.MissMassCorrDE = new TTreeReaderValue<std::vector<Double_t>>( *reader, "MissMassCorrDE" );
  
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
  //src.thetaCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaCMTPC" );
  //src.costCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "costCMTPC" );
  src.xbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xbTPC" );
  src.ybTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ybTPC" );
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );
  src.xsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xsTPC" );
  src.ysTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ysTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );  

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

  src.nclTpc = new TTreeReaderValue<Int_t>(*reader, "nclTpc" );
  src.remain_nclTpc = new TTreeReaderValue<Int_t>(*reader, "remain_nclTpc" );
  src.cluster_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_x" );
  src.cluster_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_y" );
  src.cluster_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_z" );
  src.cluster_de = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_de" );
  src.cluster_size = new TTreeReaderValue<std::vector<Int_t>>(*reader, "cluster_size" );
  src.cluster_layer = new TTreeReaderValue<std::vector<Int_t>>(*reader, "cluster_layer" );
  src.cluster_row_center = new TTreeReaderValue<std::vector<Int_t>>(*reader, "cluster_row_center" );
  src.cluster_mrow = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_mrow" );
  src.cluster_de_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_de_center" );
  src.cluster_x_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_x_center" );
  src.cluster_y_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_y_center" );
  src.cluster_z_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "cluster_z_center" );
  src.cluster_houghflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "cluster_houghflag" );
  src.remain_cluster_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_x" );
  src.remain_cluster_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_y" );
  src.remain_cluster_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_z" );
  src.remain_cluster_de = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_de" );
  src.remain_cluster_size = new TTreeReaderValue<std::vector<Int_t>>(*reader, "remain_cluster_size" );
  src.remain_cluster_layer = new TTreeReaderValue<std::vector<Int_t>>(*reader, "remain_cluster_layer" );
  src.remain_cluster_row_center = new TTreeReaderValue<std::vector<Int_t>>(*reader, "remain_cluster_row_center" );
  src.remain_cluster_mrow = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_mrow" );
  src.remain_cluster_de_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_de_center" );
  src.remain_cluster_x_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_x_center" );
  src.remain_cluster_y_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_y_center" );
  src.remain_cluster_z_center = new TTreeReaderValue<std::vector<Double_t>>(*reader, "remain_cluster_z_center" );
  src.remain_cluster_houghflag = new TTreeReaderValue<std::vector<Int_t>>(*reader, "remain_cluster_houghflag" );

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
  // src.xhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xhtofHS" );
  // src.yhtofK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "yhtofHS" );  

  src.ntTpc = new TTreeReaderValue<Int_t>(*reader, "ntTpc" );
  src.isInTarget = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isInTarget" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>(*reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>(*reader, "trackid" );
  src.isXi = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isXi" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isAccidental" );
  src.isMultiloop = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isMultiloop" );
  src.charge = new TTreeReaderValue<std::vector<Int_t>>(*reader, "charge" );
  src.pid = new TTreeReaderValue<std::vector<Int_t>>(*reader, "pid" );
  src.chisqr = new TTreeReaderValue<std::vector<Double_t>>(*reader, "chisqr" );
  src.pval = new TTreeReaderValue<std::vector<Double_t>>(*reader, "pval" );
  src.helix_cx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_cx" );
  src.helix_cy = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_cy" );
  src.helix_z0 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_z0" );
  src.helix_r = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_r" );
  src.helix_dz = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_dz" );
  src.mom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "mom_x" );
  src.mom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "mom_y" );
  src.mom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "mom_z" );
  src.dE = new TTreeReaderValue<std::vector<Double_t>>(*reader, "dE" );
  src.dEdx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "dEdx" );
  src.mom0 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "mom0" );
  src.path = new TTreeReaderValue<std::vector<Double_t>>(*reader, "path" );
  src.isElectron = new TTreeReaderValue<std::vector<Int_t>>(*reader, "isElectron" );
  src.nsigma_triton = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_triton" );
  src.nsigma_deutron = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_deutron" );
  src.nsigma_proton = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_proton" );
  src.nsigma_kaon = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_kaon" );
  src.nsigma_pion = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_pion" );
  src.nsigma_electron = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_electron" );
  src.hitlayer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitlayer" );
  src.hitpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitpos_x" );
  src.hitpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitpos_y" );
  src.hitpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "hitpos_z" );
  src.calpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "calpos_x" );
  src.calpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "calpos_y" );
  src.calpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "calpos_z" );
  src.residual = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual" );
  src.residual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_x" );
  src.residual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_y" );
  src.residual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "residual_z" );
  src.resolution_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_x");
  src.resolution_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_y");
  src.resolution_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "resolution_z");
  src.pathhit = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "pathhit");
  src.alpha = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "alpha");
  src.track_cluster_de = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_de");
  src.track_cluster_size = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_size");
  src.track_cluster_mrow = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_mrow");
  src.track_cluster_de_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_de_center");
  src.track_cluster_x_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_x_center");
  src.track_cluster_y_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_y_center");
  src.track_cluster_z_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_z_center");
  src.track_cluster_row_center = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "track_cluster_row_center");
  src.helix_t = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "helix_t" );
  src.chargeIndistinguishable = new TTreeReaderValue<std::vector<Int_t>>(*reader, "chargeIndistinguishable" );
  src.chisqr_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "chisqr_inverted" );
  src.pval_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "pval_inverted" );
  src.helix_cx_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_cx_inverted" );
  src.helix_cy_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_cy_inverted" );
  src.helix_z0_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_z0_inverted" );
  src.helix_r_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_r_inverted" );
  src.helix_dz_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "helix_dz_inverted" );
  src.mom0_inverted = new TTreeReaderValue<std::vector<Double_t>>(*reader, "mom0_inverted" );
  src.pid_inverted = new TTreeReaderValue<std::vector<Int_t>>(*reader, "pid_inverted" );

  src.nvtxTpc = new TTreeReaderValue<Int_t>(*reader, "nvtxTpc" );
  src.vtx_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "vtx_x" );
  src.vtx_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "vtx_y" );
  src.vtx_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "vtx_z" );
  src.vtx_dist = new TTreeReaderValue<std::vector<Double_t>>(*reader, "vtx_dist" );
  src.vtx_angle = new TTreeReaderValue<std::vector<Double_t>>(*reader, "vtx_angle" );
  src.vtxid = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxid" );
  src.vtxmom_theta = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxmom_theta" );
  src.vtxpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxpos_x" );
  src.vtxpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxpos_y" );
  src.vtxpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxpos_z" );
  src.vtxmom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxmom_x" );
  src.vtxmom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxmom_y" );
  src.vtxmom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "vtxmom_z" );

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

  src.ncombiLreconfailed = new TTreeReaderValue<Int_t>(*reader, "ncombiReconFailedLambda" );
  src.pidLreconfailed = new TTreeReaderValue<std::vector<Int_t>>(*reader, "ReconFailedLambdaPId");
  src.piidLreconfailed = new TTreeReaderValue<std::vector<Int_t>>(*reader, "ReconFailedLambdaPiId");
  src.LmassLreconfailed = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaMass");
  src.LdecayvtxLreconfailed_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaDecayVtx_x");
  src.LdecayvtxLreconfailed_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaDecayVtx_y");
  src.LdecayvtxLreconfailed_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaDecayVtx_z");
  src.LmomLreconfailed = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaMom" );
  src.LmomLreconfailed_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaMom_x" );
  src.LmomLreconfailed_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaMom_y" );
  src.LmomLreconfailed_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaMom_z" );
  src.pmomLreconfailed = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPMom" );
  src.pmomLreconfailed_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPMom_x" );
  src.pmomLreconfailed_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPMom_y" );
  src.pmomLreconfailed_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPMom_z" );
  src.pimomLreconfailed = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPiMom" );
  src.pimomLreconfailed_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPiMom_x" );
  src.pimomLreconfailed_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPiMom_y" );
  src.pimomLreconfailed_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdaPiMom_z" );
  src.ppidistLreconfailed = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ReconFailedLambdavtxCloseDist");  

  src.ncombiPipair = new TTreeReaderValue<Int_t>(*reader, "ncombiPiPair");
  src.pipidPipair = new TTreeReaderValue<std::vector<Int_t>>(*reader, "PiPairPipId");
  src.pimidPipair = new TTreeReaderValue<std::vector<Int_t>>(*reader, "PiPairPimId");
  src.pipmomPipair = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPipMom");
  src.pipmomPipair_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPipMom_x");
  src.pipmomPipair_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPipMom_y");
  src.pipmomPipair_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPipMom_z");
  src.pimmomPipair = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPimMom");
  src.pimmomPipair_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPimMom_x");
  src.pimmomPipair_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPimMom_y");
  src.pimmomPipair_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairPimMom_z");
  src.momPipair = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairMom");
  src.momPipair_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairMom_x");
  src.momPipair_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairMom_y");
  src.momPipair_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairMom_z");
  src.reconLmassPipair = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairReconLambdaMass");
  src.reconmassPipair = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairReconMass");
  src.pipidistPipair = new TTreeReaderValue<std::vector<Double_t>>(*reader, "PiPairCloseDist");

  src.nEscapeKm = new TTreeReaderValue<Int_t>(*reader, "nEscapeKm");
  src.kmid = new TTreeReaderValue<std::vector<Int_t>>(*reader, "KmId");
  src.kmmom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KmMom");
  src.kmmom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KmMom_x");
  src.kmmom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KmMom_y");
  src.kmmom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KmMom_z");

  src.GFkmdecayvtx_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmDecayVtx_x");
  src.GFkmdecayvtx_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmDecayVtx_y");
  src.GFkmdecayvtx_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmDecayVtx_z");
  src.GFkmmom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmMom");
  src.GFkmmom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmMom_x");
  src.GFkmmom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmMom_y");
  src.GFkmmom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmMom_z");
  src.GFkmtarget_dist = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTargetCloseDist");
  src.GFkmtargetvtx_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTarget_x");
  src.GFkmtargetvtx_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTarget_y");
  src.GFkmtargetvtx_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTarget_z");
  src.GFkmtargetcenter_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTargetCenter_x");
  src.GFkmtargetcenter_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTargetCenter_y");
  src.GFkmtargetcenter_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTargetCenter_z");
  src.GFkmtargetcenter_dist = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTargetCenterCloseDist");
  src.GFkmtracklen = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTrackLen");
  src.GFkmm2 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmMassSquare");
  src.GFkmtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFKmTof");  
  
  src.GFcharge = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFcharge");
  src.GFchisqr = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFchisqr");
  src.GFtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFtof");
  src.GFpval = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFpval");
  src.GFfitstatus = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GFfitstatus");
  src.GFpdgcode = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GFpdgcode");
  src.GFnhtrack = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GFnhtrack");
  src.GFlayer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFlayer");
  src.GFpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFpos_x");
  src.GFpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFpos_y");
  src.GFpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFpos_z");
  src.GFmom = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFmom");
  src.GFmom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFmom_x");
  src.GFmom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFmom_y");
  src.GFmom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFmom_z");
  src.GFresidual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFresidual_x");
  src.GFresidual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFresidual_y");
  src.GFresidual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFresidual_z");
  src.GFresidual_p = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFresidual_p");
  src.GFresidual_px = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFresidual_px");
  src.GFresidual_py = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFresidual_py");
  src.GFresidual_pz = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "GFresidual_pz");

  src.GFntTpc_target = new TTreeReaderValue<Int_t>(*reader, "GFntTpc_target" );
  src.GFprodvtx_x = new TTreeReaderValue<Double_t>(*reader, "GFprodvtx_x" );
  src.GFprodvtx_y = new TTreeReaderValue<Double_t>(*reader, "GFprodvtx_y" );
  src.GFprodvtx_z = new TTreeReaderValue<Double_t>(*reader, "GFprodvtx_z" );    

  src.GFinside = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GFinside" );
  src.GFKuramaFromTgt = new TTreeReaderValue<Int_t>(*reader, "GFinsideTgtToKurama" );
  src.GFKuramaVtxOutTgt = new TTreeReaderValue<Int_t>(*reader, "GFKuramaHasVtxOutOfTgt" );
  
  src.GFtracklen = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFtracklen" );
  src.GFtrack2vtxdist = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFtrack2vtxdist" );
  src.GFcalctof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFcalctof" );
  src.GFsegHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFsegHtof" );
  src.GFtofHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFtofHtof" );
  src.GFtdiffHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFtdiffHtof" );
  src.GFposHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFposHtof" );
  src.GFposx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFposx" );
  src.GFposy = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFposy" );
  src.GFposz = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFposz" );
  src.GFinvbeta = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFinvbeta" );
  src.GFm2 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFm2" );
  src.nsigma_tritonHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_tritonHtof" );
  src.nsigma_deutronHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_deutronHtof" );
  src.nsigma_protonHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_protonHtof" );
  src.nsigma_kaonHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_kaonHtof" );
  src.nsigma_pionHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_pionHtof" );
  src.nsigma_electronHtof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "nsigma_electronHtof" );

  // src.GFprodvtx_x_l = new TTreeReaderValue<Double_t>(*reader, "GFprodvtx_x_l" );
  // src.GFprodvtx_y_l = new TTreeReaderValue<Double_t>(*reader, "GFprodvtx_y_l" );
  // src.GFprodvtx_z_l = new TTreeReaderValue<Double_t>(*reader, "GFprodvtx_z_l" );      

  src.km_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KmMom_x" );
  src.km_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KmMom_y" );
  src.km_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KmMom_z" );    
  src.kp_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KpMom_x" );
  src.kp_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KpMom_y" );
  src.kp_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KpMom_z" );    
  
  src.BE = new TTreeReaderValue<std::vector<Double_t>>(*reader, "BE" );
  src.BE_LL = new TTreeReaderValue<std::vector<Double_t>>(*reader, "BE_LL" );  
  src.BETPC = new TTreeReaderValue<std::vector<Double_t>>(*reader, "BETPC" );
  src.BETPC_LL = new TTreeReaderValue<std::vector<Double_t>>(*reader, "BETPC_LL" );    
  
  src.lmass = new TTreeReaderValue<Double_t>(*reader, "LambdaMass" );
  src.ldecayvtx_x = new TTreeReaderValue<Double_t>(*reader, "LambdaDecayVtx_x" );
  src.ldecayvtx_y = new TTreeReaderValue<Double_t>(*reader, "LambdaDecayVtx_y" );
  src.ldecayvtx_z = new TTreeReaderValue<Double_t>(*reader, "LambdaDecayVtx_z" );
  src.lmom = new TTreeReaderValue<Double_t>(*reader, "LambdaMom" );
  src.lmom_x = new TTreeReaderValue<Double_t>(*reader, "LambdaMom_x" );
  src.lmom_y = new TTreeReaderValue<Double_t>(*reader, "LambdaMom_y" );
  src.lmom_z = new TTreeReaderValue<Double_t>(*reader, "LambdaMom_z" );
  src.ppi_dist = new TTreeReaderValue<Double_t>(*reader, "LambdaVtxCloseDist" );

  src.lflag = new TTreeReaderValue<Bool_t>(*reader, "Lflag" );
  src.ltarget_dist = new TTreeReaderValue<Double_t>(*reader, "LambdaTargetCloseDist" );
  src.ltargetvtx_x = new TTreeReaderValue<Double_t>(*reader, "LambdaTargetCloseVtx_x" );
  src.ltargetvtx_y = new TTreeReaderValue<Double_t>(*reader, "LambdaTargetCloseVtx_y" );
  src.ltargetvtx_z = new TTreeReaderValue<Double_t>(*reader, "LambdaTargetCloseVtx_z" );

  

  src.GFlmass = new TTreeReaderValue<Double_t>(*reader, "GFLambdaMass" );
  src.GFldecayvtx_x = new TTreeReaderValue<Double_t>(*reader, "GFLambdaDecayVtx_x" );
  src.GFldecayvtx_y = new TTreeReaderValue<Double_t>(*reader, "GFLambdaDecayVtx_y" );
  src.GFldecayvtx_z = new TTreeReaderValue<Double_t>(*reader, "GFLambdaDecayVtx_z" );
  src.GFlmom = new TTreeReaderValue<Double_t>(*reader, "GFLambdaMom" );
  src.GFlmom_x = new TTreeReaderValue<Double_t>(*reader, "GFLambdaMom_x" );
  src.GFlmom_y = new TTreeReaderValue<Double_t>(*reader, "GFLambdaMom_y" );
  src.GFlmom_z = new TTreeReaderValue<Double_t>(*reader, "GFLambdaMom_z" );
  src.GFppi_dist = new TTreeReaderValue<Double_t>(*reader, "GFLambdaVtxCloseDist" );
  src.GFltarget_dist = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTargetCloseDist" );
  src.GFltargetvtx_x = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTarget_x" );
  src.GFltargetvtx_y = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTarget_y" );
  src.GFltargetvtx_z = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTarget_z" );
  src.GFltargetcenter_x = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTargetCenter_x" );
  src.GFltargetcenter_y = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTargetCenter_y" );
  src.GFltargetcenter_z = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTargetCenter_z" );
  src.GFltargetcenter_dist = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTargetCenterCloseDist" );
  src.GFlprodvtx_x = new TTreeReaderValue<Double_t>(*reader, "GFLambdaProductionVtx_x" );
  src.GFlprodvtx_y = new TTreeReaderValue<Double_t>(*reader, "GFLambdaProductionVtx_y" );
  src.GFlprodvtx_z = new TTreeReaderValue<Double_t>(*reader, "GFLambdaProductionVtx_z" );
  src.GFlprodvtx_dist = new TTreeReaderValue<Double_t>(*reader, "GFLambdaProductionVtxCloseDist" );
  src.GFltracklen = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTrackLen" );
  src.GFltof = new TTreeReaderValue<Double_t>(*reader, "GFLambdaTof" );
  src.GFlctau = new TTreeReaderValue<Double_t>(*reader, "GFLambdaCTau" );

  src.kuramalflag = new TTreeReaderValue<Bool_t>(*reader, "LToKuramaPflag" );
  src.lpiflag = new TTreeReaderValue<Bool_t>(*reader, "LPiflag" );
  src.lpflag = new TTreeReaderValue<Bool_t>(*reader, "LPflag" );
  src.pipiflag = new TTreeReaderValue<Bool_t>(*reader, "PiPiflag" );
  src.pimflag = new TTreeReaderValue<Bool_t>(*reader, "Pimflag" );
  src.emptyflag = new TTreeReaderValue<Bool_t>(*reader, "Emptyflag" );

  src.GFdecays_pdgcode = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GFDecaysPdgcode" );
  src.GFdecays_nhtrack = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GFDecaysNhit" );
  src.GFdecays_charge = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysCharge" );
  src.GFdecays_chisqr = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysChisqr" );
  src.GFdecays_pval = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysPval" );
  src.GFdecays_htofid = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GFDecaysHtofId" );
  src.GFdecays_tracklen = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysTrackLen" );
  src.GFdecays_tof = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysTrackTof" );
  src.GFdecays_mass2 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMassSquare" );
  src.GFdecays_invbeta = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysInvbeta" );
  src.GFdecays_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMom" );
  src.GFdecays_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMom_x" );
  src.GFdecays_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMom_y" );
  src.GFdecays_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMom_z" );
  src.GFdecays_CMmom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMomCM" );
  src.GFdecays_CMmom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMomCM_x" );
  src.GFdecays_CMmom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMomCM_y" );
  src.GFdecays_CMmom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMomCM_z" );
  src.GFdecays_momloss = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysMomLoss" );
  src.GFdecays_eloss = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFDecaysELoss" );

  src.decays_id = new TTreeReaderValue<std::vector<Int_t>>(*reader, "DecaysTrackId" );
  src.decays_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMom" );
  src.decays_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMom_x" );
  src.decays_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMom_y" );
  src.decays_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMom_z" );
  src.decays_CMmom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMomCM" );
  src.decays_CMmom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMomCM_x" );
  src.decays_CMmom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMomCM_y" );
  src.decays_CMmom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "DecaysMomCM_z" );

  src.accident_multi = new TTreeReaderValue<Int_t>(*reader, "AccidentalMultiplicity" );
  src.accident_id = new TTreeReaderValue<std::vector<Int_t>>(*reader, "AccidentalTrackId" );

  // src.xiresidual_multi = new TTreeReaderValue<Int_t>(*reader, "XiResidualsMultiplicity" );
  // src.xipim_multi = new TTreeReaderValue<Int_t>(*reader, "XipimMultiplicity" );
  // src.xipip_multi = new TTreeReaderValue<Int_t>(*reader, "XipipMultiplicity" );
  // src.xiem_multi = new TTreeReaderValue<Int_t>(*reader, "XiemMultiplicity" );
  // src.xiep_multi = new TTreeReaderValue<Int_t>(*reader, "XiepMultiplicity" );
  // src.xip_multi = new TTreeReaderValue<Int_t>(*reader, "XipMultiplicity" );
  // src.xippip_multi = new TTreeReaderValue<Int_t>(*reader, "XippipMultiplicity" );
  // src.xiresidual_id = new TTreeReaderValue<std::vector<Int_t>>(*reader, "XiResidualsTrackId" );
  // src.xiresidual_dist2tgt = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsCloseDistTgt" );
  // src.xiresidual_KFdist2prodvtx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsKFProductionVtxCloseDist" );
  // src.xiresidual_GFdist2prodvtx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsGFProductionVtxCloseDist" );
  // src.xiresidual_mass2 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsMassSquare" );
  // src.xiresidual_invbeta = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsInvbeta" );
  // src.xiresidual_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsMom" );
  // src.xiresidual_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsMom_x" );
  // src.xiresidual_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsMom_y" );
  // src.xiresidual_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsMom_z" );
  // src.xiresidual_charge = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiResidualsCharge" );

  // src.llresidual_multi = new TTreeReaderValue<Int_t>(*reader, "LLResidualsMultiplicity" );
  // src.llpim_multi = new TTreeReaderValue<Int_t>(*reader, "LLpimMultiplicity" );
  // src.llpip_multi = new TTreeReaderValue<Int_t>(*reader, "LLpipMultiplicity" );
  // src.llem_multi = new TTreeReaderValue<Int_t>(*reader, "LLemMultiplicity" );
  // src.llep_multi = new TTreeReaderValue<Int_t>(*reader, "LLepMultiplicity" );
  // src.llp_multi = new TTreeReaderValue<Int_t>(*reader, "LLpMultiplicity" );
  // src.llppip_multi = new TTreeReaderValue<Int_t>(*reader, "LLppipMultiplicity" );
  // src.llresidual_id = new TTreeReaderValue<std::vector<Int_t>>(*reader, "LLResidualsTrackId" );
  // src.llresidual_dist2tgt = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsCloseDistTgt" );
  // src.llresidual_GFdist2prodvtx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsGFProductionVtxCloseDist" );
  // src.llresidual_KFdist2prodvtx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsKFProductionVtxCloseDist" );
  // src.llresidual_mass2 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsMassSquare" );
  // src.llresidual_invbeta = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsInvbeta" );
  // src.llresidual_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsMom" );
  // src.llresidual_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsMom_x" );
  // src.llresidual_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsMom_y" );
  // src.llresidual_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsMom_z" );
  // src.llresidual_charge = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLResidualsCharge" );

  src.residual_multi = new TTreeReaderValue<Int_t>(*reader, "ResidualsMultiplicity" );
  src.pim_multi = new TTreeReaderValue<Int_t>(*reader, "pimMultiplicity" );
  src.pip_multi = new TTreeReaderValue<Int_t>(*reader, "pipMultiplicity" );
  src.em_multi = new TTreeReaderValue<Int_t>(*reader, "emMultiplicity" );
  src.ep_multi = new TTreeReaderValue<Int_t>(*reader, "epMultiplicity" );
  src.p_multi = new TTreeReaderValue<Int_t>(*reader, "pMultiplicity" );
  src.ppip_multi = new TTreeReaderValue<Int_t>(*reader, "ppipMultiplicity" );
  src.residual_id = new TTreeReaderValue<std::vector<Int_t>>(*reader, "ResidualsTrackId" );
  src.residual_dist2tgt = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsCloseDistTgt" );
  src.residual_GFdist2prodvtx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsGFProductionVtxCloseDist" );
  src.residual_KFdist2prodvtx = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsKFProductionVtxCloseDist" );
  src.residual_mass2 = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsMassSquare" );
  src.residual_invbeta = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsInvbeta" );
  src.residual_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsMom" );
  src.residual_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsMom_x" );
  src.residual_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsMom_y" );
  src.residual_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsMom_z" );
  src.residual_charge = new TTreeReaderValue<std::vector<Double_t>>(*reader, "ResidualsCharge" );

  // src.KFlmom0 = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMomPpi" );
  // src.KFlmom_x0 = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMomPpi_x" );
  // src.KFlmom_y0 = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMomPpi_y0" );
  // src.KFlmom_z0 = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMomPpi_z0" );
  // src.KFlmom = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMom" );
  // src.KFlmom_x = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMom_x" );
  // src.KFlmom_y = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMom_y" );
  // src.KFlmom_z = new TTreeReaderValue<Double_t>(*reader, "KFLambdaMom_z" );
  // src.KFlchisqr = new TTreeReaderValue<Double_t>(*reader, "KFLambdaChisqr" );
  // src.KFlpval = new TTreeReaderValue<Double_t>(*reader, "KFLambdaPval" );
  // src.KFlCovMatrix = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "KFLambdaCovMatrix" );
  // src.KFlpi_dist = new TTreeReaderValue<Double_t>(*reader, "KFXiVtxCloseDist" );
  // src.KFximom = new TTreeReaderValue<Double_t>(*reader, "KFXiMom" );
  // src.KFximom_x = new TTreeReaderValue<Double_t>(*reader, "KFXiMom_x" );
  // src.KFximom_y = new TTreeReaderValue<Double_t>(*reader, "KFXiMom_y" );
  // src.KFximom_z = new TTreeReaderValue<Double_t>(*reader, "KFXiMom_z" );
  // src.KFxichisqr = new TTreeReaderValue<Double_t>(*reader, "KFXiChisqr" );
  // src.KFxipval = new TTreeReaderValue<Double_t>(*reader, "KFXiPval" );
  // src.KFxiCovMatrix = new TTreeReaderValue<std::vector<std::vector<Double_t>>>(*reader, "KFXiCovMatrix" );
  // src.KFximass = new TTreeReaderValue<Double_t>(*reader, "KFXiMass" );
  // src.KFxidecayvtx_x = new TTreeReaderValue<Double_t>(*reader, "KFXiDecayVtx_x" );
  // src.KFxidecayvtx_y = new TTreeReaderValue<Double_t>(*reader, "KFXiDecayVtx_y" );
  // src.KFxidecayvtx_z = new TTreeReaderValue<Double_t>(*reader, "KFXiDecayVtx_z" );
  // src.KFlpull = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLambdaPull" );
  // src.KFxipull = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiPull" );

  // src.KFprodvtx_chisqr_kkxi = new TTreeReaderValue<Double_t>(*reader, "KFKKXiProductionVtxChisqr" );
  // src.KFprodvtx_x_kkxi = new TTreeReaderValue<Double_t>(*reader, "KFKKXiProductionVtx_x" );
  // src.KFprodvtx_y_kkxi = new TTreeReaderValue<Double_t>(*reader, "KFKKXiProductionVtx_y" );
  // src.KFprodvtx_z_kkxi = new TTreeReaderValue<Double_t>(*reader, "KFKKXiProductionVtx_z" );
  // src.KFprodvtx_x_kpxi = new TTreeReaderValue<Double_t>(*reader, "KFKpXiProductionVtx_x" );
  // src.KFprodvtx_y_kpxi = new TTreeReaderValue<Double_t>(*reader, "KFKpXiProductionVtx_y" );
  // src.KFprodvtx_z_kpxi = new TTreeReaderValue<Double_t>(*reader, "KFKpXiProductionVtx_z" );

  // src.KFxiprodvtx_x = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_x" );
  // src.KFxiprodvtx_y = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_y" );
  // src.KFxiprodvtx_z = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_z" );
  // src.KFxiprodmom = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom" );
  // src.KFxiprodmom_x = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_x" );
  // src.KFxiprodmom_y = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_y" );
  // src.KFxiprodmom_z = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_z" );
  // src.KFxiprodvtx_dist = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxCloseDist" );
  // src.KFxitracklen = new TTreeReaderValue<Double_t>(*reader, "KFXiTrackLen" );
  // src.KFxitof = new TTreeReaderValue<Double_t>(*reader, "KFXiTof" );
  // src.KFximomloss = new TTreeReaderValue<Double_t>(*reader, "KFXiMomLoss" );
  // src.KFXiexcitation = new TTreeReaderValue<Double_t>(*reader, "KFXiExcitation" );

  // src.KFxi_kkvtx_x = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_x_KK" );
  // src.KFxi_kkvtx_y = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_y_KK" );
  // src.KFxi_kkvtx_z = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_z_KK" );
  // src.KFxi_kkvtx_mom = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_KK" );
  // src.KFxi_kkvtx_mom_x = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_x_KK" );
  // src.KFxi_kkvtx_mom_y = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_y_KK" );
  // src.KFxi_kkvtx_mom_z = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_z_KK" );
  // src.KFxi_kkvtx_dist = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxCloseDist_KK" );

  // src.KFxi_kpxiprodvtx_x = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_x_KpXi" );
  // src.KFxi_kpxiprodvtx_y = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_y_KpXi" );
  // src.KFxi_kpxiprodvtx_z = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtx_z_KpXi" );
  // src.KFxi_kpxiprodmom = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_KpXi" );
  // src.KFxi_kpxiprodmom_x = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_x_KpXi" );
  // src.KFxi_kpxiprodmom_y = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_y_KpXi" );
  // src.KFxi_kpxiprodmom_z = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxMom_z_KpXi" );
  // src.KFxi_kpxiprodvtx_dist = new TTreeReaderValue<Double_t>(*reader, "KFXiProductionVtxCloseDist_KpXi" );

  // src.KFxitargetvtx_x = new TTreeReaderValue<Double_t>(*reader, "KFXiTarget_x" );
  // src.KFxitargetvtx_y = new TTreeReaderValue<Double_t>(*reader, "KFXiTarget_y" );
  // src.KFxitargetvtx_z = new TTreeReaderValue<Double_t>(*reader, "KFXiTarget_z" );
  // src.KFxitargetmom = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetMom" );
  // src.KFxitargetmom_x = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetMom_x" );
  // src.KFxitargetmom_y = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetMom_y" );
  // src.KFxitargetmom_z = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetMom_z" );
  // src.KFxitarget_dist = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCloseDist" );

  // src.KFxitargetcenter_x = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenter_x" );
  // src.KFxitargetcenter_y = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenter_y" );
  // src.KFxitargetcenter_z = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenter_z" );
  // src.KFxitargetcentermom = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenterMom" );
  // src.KFxitargetcentermom_x = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenterMom_x" );
  // src.KFxitargetcentermom_y = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenterMom_y" );
  // src.KFxitargetcentermom_z = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenterMom_z" );
  // src.KFxitargetcenter_dist = new TTreeReaderValue<Double_t>(*reader, "KFXiTargetCenterCloseDist" );

  // src.KFlldecays_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMom" );
  // src.KFlldecays_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMom_x" );
  // src.KFlldecays_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMom_y" );
  // src.KFlldecays_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMom_z" );
  // src.KFlldecays_CMmom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMomCM" );
  // src.KFlldecays_CMmom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMomCM_x" );
  // src.KFlldecays_CMmom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMomCM_y" );
  // src.KFlldecays_CMmom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFLLDecaysMomCM_z" );

  // src.KFxidecays_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMom" );
  // src.KFxidecays_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMom_x" );
  // src.KFxidecays_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMom_y" );
  // src.KFxidecays_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMom_z" );
  // src.KFxidecays_CMmom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMomCM" );
  // src.KFxidecays_CMmom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMomCM_x" );
  // src.KFxidecays_CMmom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMomCM_y" );
  // src.KFxidecays_CMmom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFXiDecaysMomCM_z" );

  // src.KFdecays_mom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMom" );
  // src.KFdecays_mom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMom_x" );
  // src.KFdecays_mom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMom_y" );
  // src.KFdecays_mom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMom_z" );
  // src.KFdecays_CMmom = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMomCM" );
  // src.KFdecays_CMmom_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMomCM_x" );
  // src.KFdecays_CMmom_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMomCM_y" );
  // src.KFdecays_CMmom_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "KFDecaysMomCM_z" );

  // src.llg_multi = new TTreeReaderValue<Int_t>(*reader, "LLnGamma" );
  // src.llepidgamma = new TTreeReaderValue<std::vector<Int_t>>(*reader, "LLGammaEpTrackId" );
  // src.llemidgamma = new TTreeReaderValue<std::vector<Int_t>>(*reader, "LLGammaEmTrackId" );
  // src.llepmomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaMomId" );
  // src.llepmomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaEpMom_x" );
  // src.llepmomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaEpMom_y" );
  // src.llepmomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaEpMom_z" );
  // src.llemmomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaEmMom" );
  // src.llemmomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaEmMom_x" );
  // src.llemmomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaEmMom_y" );
  // src.llemmomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaEmMom_z" );
  // src.llmomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaMom" );
  // src.llmomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaMom_x" );
  // src.llmomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaMom_y" );
  // src.llmomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaMom_z" );
  // src.llepidistgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaVtxCloseDist" );
  // src.llvtxgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaDecayVtx_x" );
  // src.llvtxgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaDecayVtx_y" );
  // src.llvtxgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "LLGammaDecayVtx_z" );

  // src.xig_multi = new TTreeReaderValue<Int_t>(*reader, "XinGamma" );
  // src.xiepidgamma = new TTreeReaderValue<std::vector<Int_t>>(*reader, "XiGammaEpTrackId" );
  // src.xiemidgamma = new TTreeReaderValue<std::vector<Int_t>>(*reader, "XiGammaEmTrackId" );
  // src.xiepmomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaMomId" );
  // src.xiepmomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaEpMom_x" );
  // src.xiepmomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaEpMom_y" );
  // src.xiepmomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaEpMom_z" );
  // src.xiemmomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaEmMom" );
  // src.xiemmomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaEmMom_x" );
  // src.xiemmomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaEmMom_y" );
  // src.xiemmomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaEmMom_z" );
  // src.ximomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaMom" );
  // src.ximomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaMom_x" );
  // src.ximomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaMom_y" );
  // src.ximomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaMom_z" );
  // src.xiepidistgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaVtxCloseDist" );
  // src.xivtxgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaDecayVtx_x" );
  // src.xivtxgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaDecayVtx_y" );
  // src.xivtxgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "XiGammaDecayVtx_z" );

  src.g_multi = new TTreeReaderValue<Int_t>(*reader, "nGamma" );
  src.epidgamma = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GammaEpTrackId" );
  src.emidgamma = new TTreeReaderValue<std::vector<Int_t>>(*reader, "GammaEmTrackId" );
  src.epmomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaMomId" );
  src.epmomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaEpMom_x" );
  src.epmomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaEpMom_y" );
  src.epmomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaEpMom_z" );
  src.emmomgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaEmMom" );
  src.emmomgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaEmMom_x" );
  src.emmomgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaEmMom_y" );
  src.emmomgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaEmMom_z" );
  src.momgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaMom" );
  src.momgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaMom_x" );
  src.momgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaMom_y" );
  src.momgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaMom_z" );
  src.epidistgamma = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaVtxCloseDist" );
  src.vtxgamma_x = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaDecayVtx_x" );
  src.vtxgamma_y = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaDecayVtx_y" );
  src.vtxgamma_z = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GammaDecayVtx_z" );
  
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
      InitializeParameter<PidLikelihoodMan>("PIDLIKE") &&
      InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
