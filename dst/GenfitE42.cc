// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
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
//#include "PidLikelihoodMan.hh"
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
const int typeLHid = 1;
const int typeK0Hid = 2;

const Double_t vtx_scan_range = 150.; //ref
const Double_t vtx_scan_rangeInsideL = 50.;
const Double_t vtx_scan_rangeInsidePi = 50.;

const Double_t lambda_masscut = 0.1;
const Double_t lambda_masscut_final = 0.01; //final  
const Double_t k0_masscut = 0.1; //final
//const Double_t xi_masscut = 0.15; const Double_t lambda_masscut = 0.1; //ref
const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t e_vtx_distcut = 300;
const Double_t ppi_distcut = 10.; //ref
//const Double_t lpi_distcut = 10.;
const Double_t ltarget_distcut = 25.;
//
const Double_t pip_vtx_distcut = 300;
const Double_t pim_vtx_distcut = 300;
const Double_t pipi_distcut = 10.; //ref  
const Double_t k0target_distcut = 25.;

const Double_t GFppi_distcut = 10.;  
const Double_t GFlpi_distcut = 10.;
//const Double_t GFlpi_distcut = 15.;
const Double_t GFltarget_distcut = 25.;
const Double_t GFltarget_ycut = 20.;
const Double_t GFpipi_distcut = 10.;
const Double_t GFk0target_distcut = 25.;  
const Double_t GFk0target_ycut = 20.;

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
  std::vector<Double_t> BE;
  
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

  Bool_t lflag;
  Double_t lmass;
  Double_t ldecayvtx_x;
  Double_t ldecayvtx_y;
  Double_t ldecayvtx_z;
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

  std::vector<Double_t> GFldecays_id;
  std::vector<Double_t> GFldecays_mass2;
  std::vector<Double_t> GFldecays_invbeta;
  std::vector<Double_t> GFldecays_mom;
  

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
  

  std::vector<Int_t> GFinside;

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
    insideTPC.clear();
    pTPCKurama.clear();
    qTPCKurama.clear();
    m2TPCKurama.clear();
    xsTPC.clear();
    ysTPC.clear();
    usTPC.clear();
    vsTPC.clear();

    pK18.clear();
    xbTPC.clear();
    ybTPC.clear();
    ubTPC.clear();
    vbTPC.clear();

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

    lflag = false;
    lmass = qnan;
    ldecayvtx_x = qnan;
    ldecayvtx_y = qnan;
    ldecayvtx_z = qnan;
    lmom_x = qnan;
    lmom_y = qnan;
    lmom_z = qnan;
    ppi_dist = qnan;
    pipiangle = qnan;
    ldecays_id.clear();
    ldecays_mom.clear();  
    ldecays_mom_x.clear();
    ldecays_mom_y.clear();
    ldecays_mom_z.clear();

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

    GFldecays_id.clear();    
    GFldecays_mass2.clear();
    GFldecays_invbeta.clear();
    GFldecays_mom.clear();        

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

    GFinside.clear();

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
    if( DstRead( ievent ) ) tree->Fill();
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
  static const auto m12C = 11.174864;
  static const auto m11B = 10.252548;
  static const auto m10Be = 9.325504;
  static const auto me = 0.001*0.5109989461;
  static const int XiMinusPdgCode = 3312;
  Double_t pdgmass[3] = {ProtonMass, KaonMass, PionMass};
  TVector3 tgtpos(0, 0, tpc::ZTarget);
  TVector3 qnan_vec = TVector3(qnan, qnan, qnan);  

  static const auto KKEvent = gUser.GetParameter("KKEvent");
  static const auto KPEvent = gUser.GetParameter("KPEvent");
  static const auto KHeavyEvent = gUser.GetParameter("KHeavyEvent");
  
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

  event.insideTPC = **src.insideTPC;
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
  event.pOrgTPC = **src.pOrgTPC;
  event.pCalcTPC = **src.pCalcTPC;
  event.pCorrTPC = **src.pCorrTPC;
  event.pCorrDETPC = **src.pCorrDETPC;
  event.thetaTPC = **src.thetaTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;

  if(event.nKK != 1) return true;
  //  if(src.chisqrKurama[0] > MaxChisqrKurama || src.chisqrK18[0] > MaxChisqrBcOut) return true;
  if(KKEvent && event.Kflag[0] != 1){
    return true; //precut with Kurama tracking
  }
  if(KPEvent && event.Pflag[0] != 1){
    return true; //precut with Kurama tracking
  }
  if(KHeavyEvent && event.Heavyflag[0] != 1){
    return true; //precut with Kurama tracking
  }
  HF1( 1, event.status++ );
  if( event.nKK != 1 ) return true;  
  int ntTpc = **src.ntTpc;  
  if( ntTpc == 0 )
    return true;

  
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
  
  HF1( 1, event.status++ );
  for(int it=0; it<event.ntTpc; ++it){
    TPCLocalTrackHelix *tp = TPCAna.GetTrackTPCHelix(it);
    if( !tp ) continue;
    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], event.pid[it], pdgcode);
    GFtrackCont.AddHelixTrack(pdgcode, tp);
  }
  GFtrackCont.FitTracks();
  HF1( 2, event.GFstatus++ );
  
  int GFntTpc = GFtrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }
  HF1( 2, event.GFstatus++ );  
  HF1( genfitHid, GFntTpc);
  event.GFntTpc = GFntTpc;
  event.GFcharge.resize(GFntTpc);
  event.GFchisqr.resize(GFntTpc);
  event.GFtof.resize(GFntTpc);
  event.GFpval.resize(GFntTpc);
  event.GFfitstatus.resize(GFntTpc);
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

  //for Lambda
  std::vector<Int_t> L_p_id_container, L_pi_id_container;
  std::vector<Int_t> L_p_repid_container, L_pi_repid_container;
  std::vector<TVector3> L_p_mom_container, L_pi_mom_container;
  std::vector<TVector3> L_mom_container, L_vert_container;  
  std::vector<Double_t> L_mass_container;
  std::vector<Double_t> L_ppidist_container;
  std::vector<Double_t> L_ppiangle_container;  
  std::vector<Double_t> L_targetdist_container;
  std::vector<TVector3> L_targetvtx_container;
  //L candidates searching
  Int_t l_candidates = 0;  
  std::vector<Int_t> p_repid_container, pi_repid_container, pi2_repid_container;

  {
    for(Int_t it1=0;it1<ntTpc;it1++){ // proton
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&4)!=4) continue;
      if(event.charge[it1]!=1) continue;
      if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                  
      Int_t repid_p = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it1];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                    
      if(!GFtrackCont.TrackCheck(it1, repid_p)) continue;
      if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;              
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
      if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                    
      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isElectron[it2]==1) continue;
	if(event.isK18[it2]==1) continue;
	if(event.isKurama[it2]==1) continue;
	if(event.isBeam[it2]==1) continue;
	if(event.isAccidental[it2]==1) continue;
	if((event.pid[it2]&1)!=1) continue; //select pi-like
	//if((event.pid[it2]&4)==4) continue; //veto p-like 
	if(event.charge[it2]!=-1) continue;
	Int_t repid_pi = 0;
	if(!GFtrackCont.TrackCheck(it2, repid_pi)) continue;
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

	L_p_id_container.push_back(it1);
	L_pi_id_container.push_back(it2);
	L_mass_container.push_back(Llambda.M());
	L_mom_container.push_back(lambda_mom);
	L_p_mom_container.push_back(p_mom);
	L_pi_mom_container.push_back(pi_mom);
	L_ppidist_container.push_back(ppi_dist);
	L_vert_container.push_back(lambda_vert);
	L_p_repid_container.push_back(repid_p);
	L_pi_repid_container.push_back(repid_pi);
	L_targetdist_container.push_back(ltarget_dist);
	L_targetvtx_container.push_back(ltarget_vtx);
	l_candidates++;
      } //it2
    } //it1
  }
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;
  Int_t best_l = -1; Double_t prev_massdiff_l = 9999.;
  for(Int_t candi=0;candi<l_candidates;candi++){
    Double_t diff = TMath::Abs(L_mass_container[candi] - LambdaMass);
    if(prev_massdiff_l > diff){
      prev_massdiff_l = diff;
      best_l = candi;
      std::cout << "best lambda: " << best_l << std::endl;
    }
  }
  if(best_l!=-1){
    event.lmass = L_mass_container[best_l];
    event.ldecayvtx_x = L_vert_container[best_l].x();
    event.ldecayvtx_y = L_vert_container[best_l].y();
    event.ldecayvtx_z = L_vert_container[best_l].z();
    event.lmom_x = L_mom_container[best_l].x();
    event.lmom_y = L_mom_container[best_l].y();
    event.lmom_z = L_mom_container[best_l].z();
    event.ppi_dist = L_ppidist_container[best_l];
    //event.ppiangle = L_ppiangle_container[best_l];

    event.ldecays_id.push_back(L_p_id_container[best_l]);
    event.ldecays_mom.push_back(L_p_mom_container[best_l].Mag());
    event.ldecays_mom_x.push_back(L_p_mom_container[best_l].x());
    event.ldecays_mom_y.push_back(L_p_mom_container[best_l].y());
    event.ldecays_mom_z.push_back(L_p_mom_container[best_l].z());
    //  event.ldecays_theta.push_back(l_p_mom_container[best].Theta()*TMath::RadToDeg());

    event.ldecays_id.push_back(L_pi_id_container[best_l]);
    event.ldecays_mom.push_back(L_pi_mom_container[best_l].Mag());
    event.ldecays_mom_x.push_back(L_pi_mom_container[best_l].x());
    event.ldecays_mom_y.push_back(L_pi_mom_container[best_l].y());
    event.ldecays_mom_z.push_back(L_pi_mom_container[best_l].z());
  }
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;              
  //for K0short
  std::vector<Int_t> K0_pip_id_container, K0_pim_id_container;
  std::vector<Int_t> K0_pip_repid_container, K0_pim_repid_container;
  std::vector<TVector3> K0_pip_mom_container, K0_pim_mom_container;
  std::vector<Double_t> K0_mass_container;
  std::vector<Double_t> K0_pipidist_container;
  std::vector<Double_t> K0_pipiangle_container;  
  std::vector<Double_t> K0_targetdist_container;
  std::vector<TVector3> K0_mom_container, K0_vert_container, K0_targetvtx_container;
  //L candidates searching  
  Int_t k0_candidates = 0;
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                  
  {
    for(Int_t it1=0;it1<ntTpc;it1++){ // pi+
      if(event.isElectron[it1]==1) continue;
      if(event.isK18[it1]==1) continue;
      if(event.isKurama[it1]==1) continue;
      if(event.isBeam[it1]==1) continue;
      if(event.isAccidental[it1]==1) continue;
      if((event.pid[it1]&1)!=1) continue;
      if(event.charge[it1]!=1) continue;
      if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                      
      Int_t repid_pip = 0;
      if(!GFtrackCont.TrackCheck(it1, repid_pip)) continue;
      if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                            
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
	if((event.pid[it2]&1)!=1) continue; //select pi-like
	if(event.charge[it2]!=-1) continue;
	Int_t repid_pim = 0;
	if(!GFtrackCont.TrackCheck(it2, repid_pim)) continue;
	if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                      	
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
	TVector3 k0_vert = Kinematics::LambdaVertex(dMagneticField, pip_par, pim_par, pip_theta_min, pip_theta_max, pim_theta_min, pim_theta_max, pip_mom, pim_mom, k0_mom, pipi_dist);
	if(TMath::IsNaN(pipi_dist)) continue;
	k0_mom = pip_mom + pim_mom;
	TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
	TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));	
	TLorentzVector Lkaon0 = Lpip + Lpim;
	if(TMath::Abs(k0_vert.x()) > 250. ||
	   TMath::Abs(k0_vert.z()) > 250. ||
	   TMath::Abs(k0_vert.y()) > 250.) continue; //Vertex cut

	Double_t pip_vertex_dist; Double_t pim_vertex_dist;
	if(!Kinematics::HelixDirection(k0_vert, pip_start, pip_end, pip_vertex_dist) ||
	   !Kinematics::HelixDirection(k0_vert, pim_start, pim_end, pim_vertex_dist)) continue;

	if(pip_vertex_dist > pip_vtx_distcut) continue;
 	if(pim_vertex_dist > pim_vtx_distcut) continue;
	if(pipi_dist > pipi_distcut || TMath::Abs(Lkaon0.M() - K0Mass) > k0_masscut) continue;

	Double_t k0target_dist;
	TVector3 k0target_vtx =
	  Kinematics::CalcCloseDistLambda(tgtpos,
					  k0_vert,
					  k0_mom,
					  k0target_dist);
	
	K0_pip_id_container.push_back(it1);
	K0_pim_id_container.push_back(it2);
	K0_mass_container.push_back(Lkaon0.M());
	K0_mom_container.push_back(k0_mom);
	K0_pip_mom_container.push_back(pip_mom);
	K0_pim_mom_container.push_back(pim_mom);
	K0_pipidist_container.push_back(pipi_dist);
	//	K0_pipiangle_container.push_back(pipi_anlge);	
	K0_vert_container.push_back(k0_vert);
	K0_pip_repid_container.push_back(repid_pip);
	K0_pim_repid_container.push_back(repid_pim);
	K0_targetdist_container.push_back(k0target_dist);
	K0_targetvtx_container.push_back(k0target_vtx);
	k0_candidates++;
      } //it2
    } //it1
  }
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                      	  
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
  }
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                      	  
  Int_t ntrack_intarget = 0;
  Double_t x0[100] = {0};
  Double_t y0[100] = {0};
  Double_t u0[100] = {0};
  Double_t v0[100] = {0};
  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    event.GFfitstatus[igf] = (int)GFtrackCont.TrackCheck(igf);
    HF1( 3, event.GFfitstatus[igf]);
    if(!GFtrackCont.TrackCheck(igf)) continue;
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
    if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                      	    
    //Extrapolation
    if(event.isBeam[igf]==1) continue;
    if(event.isK18[igf]==1) continue;
    if(event.isAccidental[igf]==1) continue;
    if(GFtrackCont.IsInsideTarget(igf)){
      event.GFinside[igf] = 1;
      TVector3 posv; TVector3 momv; double len; double tof;
      if(GFtrackCont.ExtrapolateToTargetCenter(igf, posv, momv, len, tof)){
	x0[ntrack_intarget] = posv.x();
	y0[ntrack_intarget] = posv.y();
	u0[ntrack_intarget] = momv.x()/momv.z();
	v0[ntrack_intarget] = momv.y()/momv.z();
	ntrack_intarget++;
      }
    } else {
      event.GFinside[igf] = 0;
    }
  }
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;                      	
  std::vector<Double_t> GFL_mass_container(l_candidates, qnan);
  std::vector<Double_t> GFL_ppidist_container(l_candidates, qnan);
  std::vector<Double_t> GFL_targetdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetvtx_container(l_candidates, qnan_vec);
  std::vector<Double_t> GFL_targetcenterdist_container(l_candidates, qnan);
  std::vector<TVector3> GFL_targetcentervtx_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFL_mom_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFL_vert_container(l_candidates, qnan_vec);

  std::vector<Int_t> GFL_p_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_id_container(l_candidates, qnan);
  std::vector<Int_t> GFL_p_repid_container(l_candidates, qnan);
  std::vector<Int_t> GFL_pi_repid_container(l_candidates, qnan);
  std::vector<TVector3> GFL_p_mom_container(l_candidates, qnan_vec);
  std::vector<TVector3> GFL_pi_mom_container(l_candidates, qnan_vec);
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
  std::vector<Double_t> GFL_p_invbeta_container(l_candidates, qnan);
  std::vector<Double_t> GFL_pi_invbeta_container(l_candidates, qnan);
  if(l_candidates>0){
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

      Bool_t vtxcut =
	(GFtrackCont.FindVertex(p_id, pi_id,
				p_repid, pi_repid,
				p_extrapolation,
				pi_extrapolation,
				p_mom, pi_mom,
				ppi_dist, l_vertex,
				vtx_scan_range)
	 && ppi_dist < GFppi_distcut);
      if(!vtxcut) continue;

      TVector3 l_mom = p_mom + pi_mom;
      Double_t l_target_dist;
      TVector3 l_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, l_vertex, l_mom, l_target_dist);
      TVector3 l_flight = l_vertex - l_pos_tgt;
      Double_t l_tof = Kinematics::CalcTimeOfFlight(l_mom.Mag(), l_flight.Mag(), pdg::LambdaMass());
      if(l_target_dist > GFltarget_distcut) continue;

      Double_t l_targetcenter_dist;
      TVector3 l_pos_tgtcenter =
	Kinematics::LambdaTargetCenter(l_vertex, l_mom, l_targetcenter_dist);
      if(TMath::Abs(l_pos_tgtcenter.y()) > GFltarget_ycut) continue;

      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFtrackCont.TPCHTOFTrackMatching(p_id, p_repid, l_vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFL_p_htofid_container[idp] = hitid_htof;
	GFL_p_tracklen_container[idp] = tracklen_htof;
	GFL_p_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_p_mass2_container[idp] =
	  Kinematics::MassSquare(p_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	GFL_p_invbeta_container[idp] =
	  MathTools::C()*(event.tHtof[hitid_htof] - l_tof)/tracklen_htof;
      }

      htofextrapolation =
	GFtrackCont.TPCHTOFTrackMatching(pi_id, pi_repid, l_vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFL_pi_htofid_container[idp] = hitid_htof;
	GFL_pi_tracklen_container[idp] = tracklen_htof;
	GFL_pi_tof_container[idp] = event.tHtof[hitid_htof] - l_tof;
	GFL_pi_mass2_container[idp] =
	  Kinematics::MassSquare(pi_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - l_tof);
	GFL_pi_invbeta_container[idp] =
	  MathTools::C()*(event.tHtof[hitid_htof] - l_tof)/tracklen_htof;
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
      GFL_vert_container[idp] = l_vertex;
      GFL_targetdist_container[idp] = l_target_dist;
      GFL_targetvtx_container[idp] = l_pos_tgt;
      GFL_targetcenterdist_container[idp] = l_targetcenter_dist;
      GFL_targetcentervtx_container[idp] = l_pos_tgtcenter;
    }
  }
  {
    Int_t gfbest_l = -1; Double_t gfprev_massdiff = 9999.;
    Int_t gfbest_massdiff = 9999.;
    for(Int_t id=0; id<l_candidates; ++id){
      if(TMath::IsNaN(GFL_mass_container[id])) continue;
      if(TMath::IsNaN(GFL_ppidist_container[id])) continue; //Genfit's fitting was succeeded.
      if(L_targetdist_container[id] > ltarget_distcut) continue; //Select Lambda from the traget
      if(GFL_targetdist_container[id] > GFltarget_distcut) continue;
      event.lflag = true;
      Double_t diff = TMath::Abs(GFL_mass_container[id] - LambdaMass);
      if(gfprev_massdiff > diff){
	gfprev_massdiff = diff;
	gfbest_massdiff = diff;
	gfbest_l = id;
      }
    }
    if(gfbest_massdiff>lambda_masscut_final) gfbest_l=-1;
    if(gfbest_l!=-1){
      Int_t id = gfbest_l;
      event.GFlmass = GFL_mass_container[id];
      event.GFldecayvtx_x = GFL_vert_container[id].x();
      event.GFldecayvtx_y = GFL_vert_container[id].y();
      event.GFldecayvtx_z = GFL_vert_container[id].z();
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
      event.GFldecays_id.push_back(GFL_p_id_container[id]);
      event.GFldecays_id.push_back(GFL_pi_id_container[id]);         
      event.GFldecays_mass2.push_back(GFL_p_mass2_container[id]);
      event.GFldecays_mass2.push_back(GFL_pi_mass2_container[id]);
      event.GFldecays_mom.push_back(GFL_p_mom_container[id].Mag());
      event.GFldecays_mom.push_back(GFL_pi_mom_container[id].Mag());
      event.GFldecays_invbeta.push_back(GFL_p_invbeta_container[id]);
      event.GFldecays_invbeta.push_back(GFL_pi_invbeta_container[id]);
      double lmom = GFL_mom_container[id].Mag();
      double lim  = GFL_mass_container[id];
      if(lim-LambdaMass)
      HF1( 20001, lmom );
      HF1( 20101, lim );			      
      {
	double mass2 = GFL_p_mass2_container[id];
	int p_id = GFL_p_id_container[id];
	double dedx = event.dEdx[p_id];
	double GFmom = event.GFmom[p_id][0];
	int charge = event.charge[p_id];
	int type = typeLHid;
	int PID = pHid;
	int chargeid = plusHid;
	if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;
	int momid = pidlikeli::MomToBin(GFmom);
	if(!std::isnan(mass2)){
	  HF2((type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2*charge,dedx);
	}
	HF1( 22001, GFmom );
	HF1( 22011, mass2 );			
      }
      {
	double mass2 = GFL_pi_mass2_container[id];
	int pi_id = GFL_pi_id_container[id];
	double dedx = event.dEdx[pi_id];
	double GFmom = event.GFmom[pi_id][0];
	int charge = event.charge[pi_id];
	int type = typeLHid;
	int PID = piHid;
	int chargeid = minusHid;
	if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;
	int momid = pidlikeli::MomToBin(GFmom);
	if(debug){
	  std::cout << "debug " << __FILE__ << " " << __LINE__ << " "
		    << (type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m
		    << " mass2,dedx: " << mass2*charge << "," << dedx
		    << std::endl;
	}
	if(!std::isnan(mass2)){
	  HF2((type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2*charge,dedx);
	}
	HF1( 23001, GFmom );
	HF1( 23011, mass2 );			
      }      
    }
  }
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;
  std::vector<Double_t> GFK0_mass_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pipidist_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_targetdist_container(k0_candidates, qnan);
  std::vector<TVector3> GFK0_targetvtx_container(k0_candidates, qnan_vec);
  std::vector<Double_t> GFK0_targetcenterdist_container(k0_candidates, qnan);
  std::vector<TVector3> GFK0_targetcentervtx_container(k0_candidates, qnan_vec);
  std::vector<TVector3> GFK0_mom_container(k0_candidates, qnan_vec);
  std::vector<TVector3> GFK0_vert_container(k0_candidates, qnan_vec);

  std::vector<Int_t> GFK0_pip_id_container(k0_candidates, qnan);
  std::vector<Int_t> GFK0_pim_id_container(k0_candidates, qnan);
  std::vector<Int_t> GFK0_pip_repid_container(k0_candidates, qnan);
  std::vector<Int_t> GFK0_pim_repid_container(k0_candidates, qnan);
  std::vector<TVector3> GFK0_pip_mom_container(k0_candidates, qnan_vec);
  std::vector<TVector3> GFK0_pim_mom_container(k0_candidates, qnan_vec);
  std::vector<Double_t> GFK0_pip_extrapolation_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pim_extrapolation_container(k0_candidates, qnan);
  std::vector<Int_t> GFK0_pip_htofid_container(k0_candidates, qnan);
  std::vector<Int_t> GFK0_pim_htofid_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pip_tracklen_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pim_tracklen_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pip_tof_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pim_tof_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pip_mass2_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pim_mass2_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pip_invbeta_container(k0_candidates, qnan);
  std::vector<Double_t> GFK0_pim_invbeta_container(k0_candidates, qnan);
  if(k0_candidates>0){
    //Reconstructed real Lambdas
    for(int idp=0;idp<k0_candidates;idp++){
      if(K0_targetdist_container[idp] > k0target_distcut) continue;
      Int_t pip_id = K0_pip_id_container[idp];
      Int_t pim_id = K0_pim_id_container[idp];
      Int_t pip_repid = K0_pip_repid_container[idp];
      Int_t pim_repid = K0_pim_repid_container[idp];
      Double_t pip_extrapolation; Double_t pim_extrapolation;
      TVector3 pip_mom; TVector3 pim_mom;
      Double_t pipi_dist; TVector3 k0_vertex;

      Bool_t vtxcut =
	(GFtrackCont.FindVertex(pip_id, pim_id,
				pip_repid, pim_repid,
				pip_extrapolation,
				pim_extrapolation,
				pip_mom, pim_mom,
				pipi_dist, k0_vertex,
				vtx_scan_range)
	 && pipi_dist < GFpipi_distcut);
      if(!vtxcut) continue;

      TVector3 k0_mom = pip_mom + pim_mom;
      Double_t k0_target_dist;
      TVector3 k0_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, k0_vertex, k0_mom, k0_target_dist);
      TVector3 k0_flight = k0_vertex - k0_pos_tgt;
      Double_t k0_tof = Kinematics::CalcTimeOfFlight(k0_mom.Mag(), k0_flight.Mag(), K0Mass);
      if(k0_target_dist > GFk0target_distcut) continue;

      Double_t k0_targetcenter_dist;
      TVector3 k0_pos_tgtcenter =
	Kinematics::LambdaTargetCenter(k0_vertex, k0_mom, k0_targetcenter_dist);
      if(TMath::Abs(k0_pos_tgtcenter.y()) > GFk0target_ycut) continue;

      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFtrackCont.TPCHTOFTrackMatching(pip_id, pim_repid, k0_vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFK0_pip_htofid_container[idp] = hitid_htof;
	GFK0_pip_tracklen_container[idp] = tracklen_htof;
	GFK0_pip_tof_container[idp] = event.tHtof[hitid_htof] - k0_tof;
	GFK0_pip_mass2_container[idp] =
	  Kinematics::MassSquare(pip_mom.Mag(),tracklen_htof,event.tHtof[hitid_htof] - k0_tof);
	GFK0_pip_invbeta_container[idp] =
	  MathTools::C()*(event.tHtof[hitid_htof] - k0_tof)/tracklen_htof;
      }

      htofextrapolation =
	GFtrackCont.TPCHTOFTrackMatching(pim_id, pim_repid, k0_vertex,
					 event.HtofSeg, event.posHtof,
					 hitid_htof, tof_htof,
					 tracklen_htof, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	GFK0_pim_htofid_container[idp] = hitid_htof;
	GFK0_pim_tracklen_container[idp] = tracklen_htof;
	GFK0_pim_tof_container[idp] = event.tHtof[hitid_htof] - k0_tof;
	GFK0_pim_mass2_container[idp] =
	 Kinematics::MassSquare(pim_mom.Mag(), tracklen_htof, event.tHtof[hitid_htof] - k0_tof);
	GFK0_pim_invbeta_container[idp] =
	  MathTools::C()*(event.tHtof[hitid_htof] - k0_tof)/tracklen_htof;
      }

      TLorentzVector GFLpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
      TLorentzVector GFLpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
      TLorentzVector GFLkaon0 = GFLpip + GFLpim;

      GFK0_pip_id_container[idp] = pip_id;
      GFK0_pim_id_container[idp] = pim_id;
      GFK0_pip_repid_container[idp] = pip_repid;
      GFK0_pim_repid_container[idp] = pim_repid;
      GFK0_pip_mom_container[idp] = pip_mom;
      GFK0_pim_mom_container[idp] = pim_mom;
      GFK0_pip_extrapolation_container[idp] = pip_extrapolation;
      GFK0_pim_extrapolation_container[idp] = pim_extrapolation;
      GFK0_mass_container[idp] = GFLkaon0.M();
      GFK0_mom_container[idp] = k0_mom;
      GFK0_pipidist_container[idp] = pipi_dist;
      GFK0_vert_container[idp] = k0_vertex;
      GFK0_targetdist_container[idp] = k0_target_dist;
      GFK0_targetvtx_container[idp] = k0_pos_tgt;
      GFK0_targetcenterdist_container[idp] = k0_targetcenter_dist;
      GFK0_targetcentervtx_container[idp] = k0_pos_tgtcenter;
    }
  }
  if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;  
  {
    Int_t gfbest_k0 = -1; Double_t gfprev_massdiff_k0 = 9999.;
    for(Int_t id=0; id<k0_candidates; ++id){
      if(TMath::IsNaN(GFK0_mass_container[id])) continue;
      if(TMath::IsNaN(GFK0_pipidist_container[id])) continue; //Genfit's fitting was succeeded.
      if(K0_targetdist_container[id] > k0target_distcut) continue; //Select Lambda from the traget
      if(GFK0_targetdist_container[id] > GFk0target_distcut) continue;
      event.k0flag = true;
      if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;        
      Double_t diff = TMath::Abs(GFK0_mass_container[id] - K0Mass);
      if(gfprev_massdiff_k0 > diff){
	gfprev_massdiff_k0 = diff;
	gfbest_k0 = id;
      }
    }
    if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;      
    if(gfbest_k0!=-1){
      Int_t id = gfbest_k0;
      event.GFk0mass = GFK0_mass_container[id];
      event.GFk0decayvtx_x = GFK0_vert_container[id].x();
      event.GFk0decayvtx_y = GFK0_vert_container[id].y();
      event.GFk0decayvtx_z = GFK0_vert_container[id].z();
      event.GFk0mom = GFK0_mom_container[id].Mag();
      event.GFk0mom_x = GFK0_mom_container[id].x();
      event.GFk0mom_y = GFK0_mom_container[id].y();
      event.GFk0mom_z = GFK0_mom_container[id].z();
      event.GFk0pipi_dist = GFK0_pipidist_container[id];
      event.GFk0target_dist = GFK0_targetdist_container[id];
      event.GFk0targetvtx_x = GFK0_targetvtx_container[id].x();
      event.GFk0targetvtx_y = GFK0_targetvtx_container[id].y();
      event.GFk0targetvtx_z = GFK0_targetvtx_container[id].z();
      event.GFk0targetcenter_dist = GFK0_targetcenterdist_container[id];
      event.GFk0targetcenter_x = GFK0_targetcentervtx_container[id].x();
      event.GFk0targetcenter_y = GFK0_targetcentervtx_container[id].y();
      event.GFk0targetcenter_z = GFK0_targetcentervtx_container[id].z();
      event.GFk0decays_id.push_back(GFK0_pip_id_container[id]);
      event.GFk0decays_id.push_back(GFK0_pim_id_container[id]);         
      event.GFk0decays_mass2.push_back(GFK0_pip_mass2_container[id]);
      event.GFk0decays_mass2.push_back(GFK0_pim_mass2_container[id]);
      event.GFk0decays_mom.push_back(GFK0_pip_mom_container[id].Mag());
      event.GFk0decays_mom.push_back(GFK0_pim_mom_container[id].Mag());
      event.GFk0decays_invbeta.push_back(GFK0_pip_invbeta_container[id]);
      event.GFk0decays_invbeta.push_back(GFK0_pim_invbeta_container[id]);
      
      double k0mom = GFK0_mom_container[id].Mag();
      double k0mass = GFK0_mass_container[id];      
      HF1( 30001, k0mom );
      HF1( 30101, k0mass );		      
      {
	double mass2 = GFK0_pip_mass2_container[id];
	int pip_id = GFK0_pip_id_container[id];
	double dedx = event.dEdx[pip_id];
	double GFmom = event.GFmom[pip_id][0];
	int charge = event.charge[pip_id];	
	int type = typeK0Hid;
	int PID = piHid;
	int chargeid = plusHid;
	if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;
	int momid = pidlikeli::MomToBin(GFmom);
	HF2( (type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2*charge,dedx);
	HF1( 32001, GFmom );
	HF1( 32011, mass2 );		
      }
      {
	double mass2 = GFK0_pim_mass2_container[id];
	int pim_id = GFK0_pim_id_container[id];
	double dedx = event.dEdx[pim_id];
	double GFmom = event.GFmom[pim_id][0];
	int charge = event.charge[pim_id];	
	int type = typeK0Hid;
	int PID = piHid;
	int chargeid = minusHid;
	if(debug) std::cout << "debug " << __FILE__ << " " << __LINE__ << std::endl;
	int momid = pidlikeli::MomToBin(GFmom);
	HF2( (type+1)*fac_t + PID*fac_p + chargeid*fac_c + 0*fac_b + momid*fac_m, mass2*charge,dedx); 
	HF1( 33001, GFmom );
	HF1( 33011, mass2 );			
      }            
    }
  }
  
  TVector3 vertex = Kinematics::MultitrackVertex(ntrack_intarget,x0,y0,u0,v0);
  event.GFntTpc_inside = ntrack_intarget;
  event.GFprodvtx_x = vertex.x();
  event.GFprodvtx_y = vertex.y();
  event.GFprodvtx_z = vertex.z();
  
  for( Int_t igf=0; igf<GFntTpc; ++igf ){  
    // htof ana
    {
      Int_t repid = -1;
      Int_t hitid_htof; Double_t tof; Double_t len;
      TVector3 pos_htof; Double_t track2tgt_dist;
      Bool_t htofextrapolation =
	GFtrackCont.TPCHTOFTrackMatching(igf, repid, vertex,
				      event.HtofSeg, event.posHtof,
				      hitid_htof, tof,
				      len, pos_htof, track2tgt_dist);
      if(htofextrapolation){
	event.GFextrapolationHtof[igf] = 1;
	event.GFtracklen[igf] = len;
	event.GFtrack2vtxdist[igf] = track2tgt_dist;
	event.GFcalctof[igf] = tof;
	event.GFposx[igf] = pos_htof.x();
	event.GFposy[igf] = pos_htof.y();
	event.GFposz[igf] = pos_htof.z();
	event.GFsegHtof[igf] = event.HtofSeg[hitid_htof];
	event.GFtofHtof[igf] = event.tHtof[hitid_htof];
	event.GFtdiffHtof[igf] = event.dtHtof[hitid_htof];
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
	event.nsigma_electronHtof[igf]
	  = Kinematics::HypTPCHTOFNsigmaElectron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
      }
    } //common
  } //igf
    
  // make 2D hist for PDF
  for( Int_t itTpc=0; itTpc<event.GFntTpc; ++itTpc ){
    if(!GFtrackCont.TrackCheck(itTpc)) continue;
    double invbeta = event.GFinvbeta[itTpc];
    if(invbeta<0.1) continue;
    double m2 = event.GFm2[itTpc];
    //    double m2 = event.GFm2[itTpc];
    double nsigma_t  = event.nsigma_triton[itTpc]   ;
    double nsigma_d  = event.nsigma_deutron[itTpc]  ;
    double nsigma_p  = event.nsigma_proton[itTpc]   ;
    double nsigma_k  = event.nsigma_kaon[itTpc]     ;
    double nsigma_pi = event.nsigma_pion[itTpc]     ;
    double nsigma_e  = event.nsigma_electron[itTpc] ;
    
    double nsigmaHtof_t  = event.nsigma_tritonHtof[itTpc]   ;
    double nsigmaHtof_d  = event.nsigma_deutronHtof[itTpc]  ;
    double nsigmaHtof_p  = event.nsigma_protonHtof[itTpc]   ;
    double nsigmaHtof_k  = event.nsigma_kaonHtof[itTpc]     ;
    double nsigmaHtof_pi = event.nsigma_pionHtof[itTpc]     ;
    double nsigmaHtof_e  = event.nsigma_electronHtof[itTpc] ;

    int seghtof = event.GFsegHtof[itTpc];
    // tpc
    // Int_t pid = event.pid[itTpc];
    Int_t charge = event.charge[itTpc];
    double dEdxtpc = event.dEdx[itTpc];
    //      double gfposy = event.GFposy[itTpc];
    double be = 0.;
    double cut_tofpid_min=5.0;   //should be better written in USER param
    double cut_tofpid_max=5.0;   //should be better written in USER param  
    double cut_dedxpid_min=5.0;   //should be better written in USER param
    double cut_dedxpid_max=5.0;   //should be better written in USER param   
    Bool_t flag_pi=false;        //should be better written in USER param
    Bool_t flag_k=false;         //should be better written in USER param
    Bool_t flag_p=false;         //should be better written in USER param
    Bool_t flag_d=false;         //should be better written in USER param    
    Bool_t flag_e=false;         //should be better written in USER param 
    if( nsigmaHtof_pi > -cut_tofpid_min && nsigmaHtof_pi < cut_tofpid_max
	&& nsigma_pi > -cut_dedxpid_min && nsigma_pi < cut_dedxpid_max
	&& !event.isElectron[itTpc] ){	
      flag_pi = true;
    }
    if( nsigmaHtof_k > -cut_tofpid_min && nsigmaHtof_k < cut_tofpid_max
	&& nsigma_k > -cut_dedxpid_min && nsigma_k < cut_dedxpid_max
	&& !event.isElectron[itTpc] ){	
      flag_k = true;
    }
    if( nsigmaHtof_p > -cut_tofpid_min && nsigmaHtof_p < cut_tofpid_max
	&& nsigma_p > -cut_dedxpid_min && nsigma_p < cut_dedxpid_max
	&& !event.isElectron[itTpc] ){
      flag_p = true;
    }
    if( nsigmaHtof_d > -cut_tofpid_min && nsigmaHtof_d < cut_tofpid_max
	&& nsigma_d > -cut_dedxpid_min && nsigma_d < cut_dedxpid_max
	&& !event.isElectron[itTpc] ){
      flag_d = true;
    }    
    if( event.isElectron[itTpc] ){
      flag_e = true;
    }
    double GFmom = event.GFmom[itTpc][0];
    // if(dist > dist_cut) continue;
    // if(chisqr>chisqrcut) continue;
    if(event.isBeam[itTpc] || event.isK18[itTpc] || event.isAccidental[itTpc]) continue;
    // if(!event.GFinside[itTpc]) continue;
    // select type, particle, charge, BE, momentum
    constexpr double eps = std::numeric_limits<double>::epsilon();
    // if (be < minbe || be >= maxbe ) continue;
    //BE
    int beid = -1;
    if(1){ // temp
      beid=0;
    } else {
      beid = static_cast<int>( std::floor( (be - minbe + eps) / bestep ) );
      if (beid < 0)       beid = 0;
      if (beid >= nbinbe) beid = nbinbe - 1;
    }      
    //mom    
    int momid = pidlikeli::MomToBin(GFmom);
    //    std::cout << "GFmom:" << GFmom <<", momid:" << momid << std::endl;
    //static_cast<int>( std::floor( (GFmom + eps) / momstep ) );
    if (momid < 0) continue;
    if (momid >= nbinpoq) momid = nbinpoq - 1;
    Int_t type=1;
    Int_t typetpcxp = 4;
    Int_t typetpcxm = 5;    
    //      if(KPEvent) type=2;
    Int_t pid=-1;
    int chargeid = (charge>0) ? 0 : (charge<0) ? 1 : -2;
    if(chargeid==-2) continue;
    if(flag_pi){
      pid=0;
      // std::cout << type*pdfHid+pid*100000+chargeid*10000+beid*100+momid
      // 		<< " " << invbeta << " " <<  dEdxtpc << std::endl;
      HF2(type*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2, dEdxtpc);
    }
    if(flag_k){
      pid=1;
      HF2(type*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
      if( (seghtof>=12&&seghtof<=19) ){
	HF2(typetpcxm*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
      } else if( (seghtof>=20&&seghtof<=27) ){
	// for TPC +x/-x study	
	//HF2(typetpcxp*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
      }      
    }
    if(flag_p){
      pid=2;
      HF2(type*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
      if( (seghtof>=12&&seghtof<=19) ){
	HF2(typetpcxm*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
      } else if( (seghtof>=20&&seghtof<=27) ){
	// for TPC +x/-x study		
	//HF2(typetpcxp*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
      }
    }
    if(flag_d){
      pid=3;
      HF2(type*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
    }
    if(flag_e){
      pid=4;
      HF2(type*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
    }
    if(!flag_e){
      pid=5; // all but e
      HF2(type*fac_t+pid*fac_p+chargeid*fac_c+beid*fac_b+momid*fac_m, m2*charge, dEdxtpc);
    }
  }  
  HF1( 2, event.GFstatus++);
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
    HB1( 2, "Genfit Status", 20, 0., 20. );
    HB1( 3, "Genfit Fit Status", 2, 0., 2. );
    HB1( 10, "NTrack TPC", 20, 0., 20. );

    HB1(genfitHid, "[GenFit] #Track TPC; #Track; Counts", 20, 0., 20. );
    HB1(genfitHid+1, "[GenFit] Chisqr/ndf; ; Counts", 200, 0, 10 );
    HB1(genfitHid+2, "[GenFit] p-value; p-value; Counts", 100, -0.05, 1.05);
    HB1(genfitHid+3, "[GenFit] Charge;", 6, -3, 3 );
    HB1(genfitHid+4, "[GenFit] #Hits of Track", 50, 0., 50. );
    HB1(genfitHid+5, "[GenFit] Track Length; Length [mm]; Counts [/1 mm]", 500, 0, 500 );
    HB1(genfitHid+6, "[GenFit] Tof; Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
    HB1(genfitHid+7, "[GenFit] Reconstructed P; P [GeV/c]; Counts [/0.001 GeV/c]", 1500, 0., 1.5 );
    HB1(genfitHid+8, "[GenFit] LayerID", 33, 0., 33. );

    //Residuals
    HB1(genfitHid+10, "[GenFit] Residual x; Residual x [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+11, "[GenFit] Residual y; Residual y [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+12, "[GenFit] Residual z; Residual z [mm]; Counts [/0.1 mm]", 400, -20., 20.);
    HB1(genfitHid+13, "[GenFit] Residual P; Residual P [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
    HB1(genfitHid+14, "[GenFit] Residual Px; Residual Px [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(genfitHid+15, "[GenFit] Residual Py; Residual Py [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  HB1(genfitHid+16, "[GenFit] Residual Pz; Residual Pz [GeV/c]; Counts [/0.001 GeV/c]", 400, -0.2, 0.2 );
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    HB1(genfitHid+(layer+1)*1000, Form("[GenFit] Residual x Layer%d; Residual x [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+1, Form("[GenFit] Residual y Layer%d; Residual y [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+2, Form("[GenFit] Residual z Layer%d; Residual z [mm]; Counts [/0.1 mm]", layer), 400, -20., 20.);
    HB1(genfitHid+(layer+1)*1000+3, Form("[GenFit] Residual P Layer%d; Residual P [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+4, Form("[GenFit] Residual Px Layer%d; Residual Px [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+5, Form("[GenFit] Residual Py Layer%d; Residual Py [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
    HB1(genfitHid+(layer+1)*1000+6, Form("[GenFit] Residual Pz Layer%d; Residual Pz [GeV/c]; Counts [/0.001 GeV/c]", layer), 400, -0.2, 0.2 );
  }

  //Extrapolation
  HB1(genfitHid+20, "[GenFit] Vertex X (extrapolated to the Target); X [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+21, "[GenFit] Vertex Y (extrapolated to the Target); Vertex Y [mm]; Counts [/0.1 mm]", 500, -25, 25 );
  HB1(genfitHid+22, "[GenFit] Vertex Z (extrapolated to the Target); Vertex Z [mm]; Counts [/0.1 mm]", 500, -143 -25, -143 + 25 );
  HB1(genfitHid+23, "[GenFit] #Hits (extrapolated to the HTOF); #Hits; Counts", 10, 0, 10 );
  HB1(genfitHid+24, "[GenFit] Track Length (extrapolated to the HTOF); Length [mm]; Counts [/1 mm]", 1000, -500, 500 );
  HB1(genfitHid+25, "[GenFit] Tof (extrapolated to the HTOF); Tof [ns]; Counts [/0.01 ns]", 500, 0, 5 );
  HB1(genfitHid+26, "[GenFit] X (extrapolated to the HTOF); X [mm]; Counts [/0.1 mm]", 10000, -500, 500 );
  HB1(genfitHid+27, "[GenFit] Y (extrapolated to the HTOF); Y [mm]; Counts [/0.1 mm]", 8000, -400, 400 );
  HB1(genfitHid+28, "[GenFit] Z (extrapolated to the HTOF); Z [mm]; Counts [/0.1 mm]", 10000, -500, 500 );
  HB1(genfitHid+29, "[GenFit] HTOF ID; #ID ; Counts", 36, 0, 36 );
*/
  // PDF for pi,K,p for each momentum region
  // int nbinpoq = 1000;
  // double minpoq = -1.5;
  // double maxpoq = 1.5; //GeV/c
  // int nbininvbeta = 100;
  // double mininvbeta = 0.;
  // double maxinvbeta = 5.;
  // int nbindedx = 1000;
  // double mindedx = 0.;
  // double maxdedx = 350.;
  // int nbinbe = 400;
  // double minbe = -0.5; //GeV
  // double maxbe = 0.5;
  // int pdfHid = 1000000;
  // double momstep = 0.05; //GeV/c
  // int momId = (Int_t)(maxpoq/momstep); // 0.05 GeV/c step?
  // double bestep = 0.20; //GeV // should be chaged to USER parameter
  // int BEId = (Int_t)((maxbe-minbe)/bestep); // should be less than
  
  //  TString pid[kNpid+1] = {"pi", "K", "p", "e", "d", "all"};
  //  TString pid[kNpid+1] = {"pi", "K", "p","all"};  
  //  TString charge[kNchg] = {"+", "-"};
  //  TString type[kNtype] = {"general", "KP"};
  
  for(int itype=0; itype<kNtype; itype++){//1:general, 2:Lambda reconstruct,3:K0 reconstruct,
    if(itype>3) continue;
    for(int ipid=0; ipid<kNpid+1; ipid++){ // pi,K,p,d,e,all
      for(int icharge=0; icharge<kNchg; icharge++){ // posi,nega
	if( itype==3||itype==4 ){
	  if( !( (ipid==1&&icharge==1)||(ipid==2&&icharge==0) ) ) continue;
	}	
        for(int ibe=0; ibe<kNbe; ibe++){ // default: beid=0
          for(int imom=0; imom<kNmom; imom++){
	    if(ipid<kNpid){
	      if(itype==typeLHid&&(ipid==kHid||ipid==dHid||ipid==eHid)) continue;
	      if(itype==typeK0Hid&&(ipid==pHid||ipid==kHid||ipid==dHid||ipid==eHid)) continue;	      
	      HB2( (itype+1)*fac_t + ipid*fac_p + icharge*fac_c + ibe*fac_b + imom*fac_m,
		   Form("PDF %s %s %s BE=%.3fGeV mom=%.4fGeV/c; mass2 ; dEdx",
			type[itype].Data(), plist[ipid].Data(),clist[icharge].Data(),
			minbe+bestep*(double)(ibe),kMomstep*(Double_t(imom))),
		   nbinm2, minm2, maxm2, nbindedx, mindedx, maxdedx);
	    } else {
	      if(itype==typeLHid&&(ipid==kHid||ipid==dHid||ipid==eHid)) continue;
	      if(itype==typeK0Hid&&(ipid==pHid||ipid==kHid||ipid==dHid||ipid==eHid)) continue;
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

  HB1(20001, "GF#Lambda mass",1000,pdg::LambdaMass()-0.2,pdg::LambdaMass()+0.2);
  HB1(20002, "GF#Lambda mass selected",1000,pdg::LambdaMass()-0.2,pdg::LambdaMass()+0.2);  
  HB1(20011, "GF#Lambda vtx",100,-5.,5.);
  HB1(20012, "GF#Lambda vty",100,-5.,5.);
  HB1(20013, "GF#Lambda vtz",100,-5.,5.);
  HB1(20101, "GF#Lambda mom",1500,0.,1.5);
  HB1(22001, "GF#Lambda p mom; Mom[GeV/c]",nbinpoq,0.,maxpoq);
  HB1(22011, "GF#Lambda p m2; Mom[GeV/c]",nbinm2,minm2,maxm2);
  HB1(23001, "GF#Lambda #pi mom; Mom[GeV/c]",nbinpoq,0.,maxpoq);
  HB1(23011, "GF#Lambda #pi m2; Mom[GeV/c]",nbinm2,minm2,maxm2);
  
  HB1(30001, "GFK0 mass",1000,pdg::K0Mass()-0.2,pdg::K0Mass()+0.2);
  HB1(30011, "GFK0 vtx",100,-5.,5.);
  HB1(30012, "GFK0 vty",100,-5.,5.);
  HB1(30013, "GFK0 vtz",100,-5.,5.);
  HB1(30101, "GFK0 mom",1500,0.,1.5);
  HB1(32001, "GFK0 #pi+ mom; Mom[GeV/c]",nbinpoq,0.,maxpoq);
  HB1(32011, "GFK0 #pi+ m2; Mom[GeV/c]",nbinm2,minm2,maxm2);
  HB1(33001, "GFK0 #pi- mom; Mom[GeV/c]",nbinpoq,0.,maxpoq);
  HB1(33011, "GFK0 #pi- m2; Mom[GeV/c]",nbinm2,minm2,maxm2);
  
  HBTree( "tpc", "tree of GenfitE42" );
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
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isAccidental", &event.isAccidental );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "pid", &event.pid );  
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "path", &event.path );  
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );  
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
  tree->Branch( "resolution_x", &event.resolution_x );
  tree->Branch( "resolution_y", &event.resolution_y );
  tree->Branch( "resolution_z", &event.resolution_z );
  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "pathhit", &event.pathhit);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_size", &event.track_cluster_size);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch( "track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch( "track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch( "track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch( "track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch( "track_cluster_row_center", &event.track_cluster_row_center);

  tree->Branch( "ntK18", &event.ntK18);
  tree->Branch( "chisqrK18", &event.chisqrK18);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "xtgtK18", &event.xtgtK18);
  tree->Branch( "ytgtK18", &event.ytgtK18);
  tree->Branch( "utgtK18", &event.utgtK18);
  tree->Branch( "vtgtK18", &event.vtgtK18);
  tree->Branch( "thetaK18", &event.thetaK18);

  tree->Branch( "ntKurama",     &event.ntKurama);
  tree->Branch( "chisqrKurama", &event.chisqrKurama);
  tree->Branch( "pKurama",      &event.pKurama);
  tree->Branch( "qKurama",      &event.qKurama);
  tree->Branch( "m2",           &event.m2);
  //  tree->Branch( "m2Org",        &event.m2Org);
  tree->Branch( "xtgtKurama",   &event.xtgtKurama);
  tree->Branch( "ytgtKurama",   &event.ytgtKurama);
  tree->Branch( "utgtKurama",   &event.utgtKurama);
  tree->Branch( "vtgtKurama",   &event.vtgtKurama);  

  tree->Branch( "tpcidTPCKurama", &event.tpcidTPCKurama);  
  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "insideTPC", &event.insideTPC);
  tree->Branch( "pTPCKurama", &event.pTPCKurama);
  tree->Branch( "qTPCKurama", &event.qTPCKurama);
  tree->Branch( "m2TPCKurama", &event.m2TPCKurama);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "xtgtTPCKurama",   &event.xtgtTPCKurama);
  tree->Branch( "ytgtTPCKurama",   &event.ytgtTPCKurama);
  tree->Branch( "utgtTPCKurama",   &event.utgtTPCKurama);
  tree->Branch( "vtgtTPCKurama",   &event.vtgtTPCKurama);  
  
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
  //  tree->Branch("pOrg",       &event.pOrg);
  //  tree->Branch("pCalc",      &event.pCalc);
  //  tree->Branch("pCorr",      &event.pCorr);
  //  tree->Branch("pCorrDE",    &event.pCorrDE);
  tree->Branch("Kflag",      &event.Kflag);
  tree->Branch("Pflag",      &event.Pflag);
  tree->Branch("Heavyflag",  &event.Heavyflag);  

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
  tree->Branch( "thetaTPC", &event.thetaTPC);
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);

  tree->Branch("Lflag", &event.lflag);
  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);
  tree->Branch("LambdaPPiAngle", &event.ppiangle);  
  tree->Branch("LDecaysTrackId", &event.ldecays_id);
  tree->Branch("LDecaysMom", &event.ldecays_mom);
  tree->Branch("LDecaysMom_x", &event.ldecays_mom_x);
  tree->Branch("LDecaysMom_y", &event.ldecays_mom_y);
  tree->Branch("LDecaysMom_z", &event.ldecays_mom_z);  
  
  tree->Branch("K0flag", &event.k0flag);
  tree->Branch("K0Mass", &event.k0mass);
  tree->Branch("K0DecayVtx_x", &event.k0decayvtx_x);
  tree->Branch("K0DecayVtx_y", &event.k0decayvtx_y);
  tree->Branch("K0DecayVtx_z", &event.k0decayvtx_z);
  tree->Branch("K0Mom_x", &event.k0mom_x);
  tree->Branch("K0Mom_y", &event.k0mom_y);
  tree->Branch("K0Mom_z", &event.k0mom_z);
  tree->Branch("K0VtxCloseDist", &event.GFk0pipi_dist);
  //  tree->Branch("K0PiPiAngle", &event.pipiangle);  
  tree->Branch("K0DecaysTrackId", &event.k0decays_id);
  tree->Branch("K0DecaysMom", &event.k0decays_mom);
  tree->Branch("K0DecaysMom_x", &event.k0decays_mom_x);
  tree->Branch("K0DecaysMom_y", &event.k0decays_mom_y);
  tree->Branch("K0DecaysMom_z", &event.k0decays_mom_z);  
  
  //track fitting results
  tree->Branch("GFstatus", &event.GFstatus);
  tree->Branch("GFntTpc", &event.GFntTpc);
  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFtof", &event.GFtof);
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

  tree->Branch("GFK0Mass", &event.GFk0mass);
  tree->Branch("GFK0DecayVtx_x", &event.GFk0decayvtx_x);
  tree->Branch("GFK0DecayVtx_y", &event.GFk0decayvtx_y);
  tree->Branch("GFK0DecayVtx_z", &event.GFk0decayvtx_z);
  tree->Branch("GFK0Mom", &event.GFk0mom);
  tree->Branch("GFK0Mom_x", &event.GFk0mom_x);
  tree->Branch("GFK0Mom_y", &event.GFk0mom_y);
  tree->Branch("GFK0Mom_z", &event.GFk0mom_z);
  tree->Branch("GFK0VtxCloseDist", &event.GFk0pipi_dist);
  tree->Branch("GFK0TargetCloseDist", &event.GFk0target_dist);
  tree->Branch("GFK0Target_x", &event.GFk0targetvtx_x);
  tree->Branch("GFK0Target_y", &event.GFk0targetvtx_y);
  tree->Branch("GFK0Target_z", &event.GFk0targetvtx_z);
  tree->Branch("GFK0TargetCenter_x", &event.GFk0targetcenter_x);
  tree->Branch("GFK0TargetCenter_y", &event.GFk0targetcenter_y);
  tree->Branch("GFK0TargetCenter_z", &event.GFk0targetcenter_z);
  tree->Branch("GFK0TargetCenterCloseDist", &event.GFk0targetcenter_dist);
  tree->Branch("GFK0ProductionVtx_x", &event.GFk0prodvtx_x);
  tree->Branch("GFK0ProductionVtx_y", &event.GFk0prodvtx_y);
  tree->Branch("GFK0ProductionVtx_z", &event.GFk0prodvtx_z);
  tree->Branch("GFK0ProductionVtxCloseDist", &event.GFk0prodvtx_dist);
  tree->Branch("GFK0TrackLen", &event.GFk0tracklen);
  tree->Branch("GFK0Tof", &event.GFk0tof);  
  
  tree->Branch("GFntTpc_inside", &event.GFntTpc_inside);
  tree->Branch("GFprodvtx_x", &event.GFprodvtx_x);
  tree->Branch("GFprodvtx_y", &event.GFprodvtx_y);
  tree->Branch("GFprodvtx_z", &event.GFprodvtx_z);

  //extrapolation
  tree->Branch("GFinside", &event.GFinside);
  tree->Branch("GFextrapolationHtof", &event.GFextrapolationHtof); 
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
  src.insideTPC = new TTreeReaderValue<std::vector<Int_t>>( *reader, "insideTPC" );
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
  src.pOrgTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pOrgTPC" );
  src.pCalcTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCalcTPC" );
  src.pCorrTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrTPC" );
  src.pCorrDETPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pCorrDETPC" );
  src.thetaTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaTPC" );
 
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
