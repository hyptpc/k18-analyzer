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
}


namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kGenfitE42, kOutFile, nArgc
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

  std::vector<Int_t> isgoodTPCKurama;
  std::vector<Double_t> pTPCKurama;
  std::vector<Double_t> qTPCKurama;
  std::vector<Double_t> m2TPCKurama;
  std::vector<Double_t> xsTPC;
  std::vector<Double_t> ysTPC;
  std::vector<Double_t> usTPC;
  std::vector<Double_t> vsTPC;

  std::vector<Double_t> pK18;
  std::vector<Double_t> xbTPC;
  std::vector<Double_t> ybTPC;
  std::vector<Double_t> ubTPC;
  std::vector<Double_t> vbTPC;

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
  //std::vector<Double_t> GFtracklen;
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
  std::vector<Int_t> GFfromVtx;  
  std::vector<Int_t> GFextrapolationHtof;

  Int_t GFntTpc_inside;
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
  std::vector<Int_t>    ldecays_htofhitid;
  std::vector<Int_t>    ldecays_htofseg;    
  std::vector<Double_t> ldecays_mass2;
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

  Bool_t GFlflag;  
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
  
  std::vector<Int_t> GFldecays_id;
  std::vector<Double_t> GFldecays_mom;
  std::vector<Double_t> GFldecays_mom_x;
  std::vector<Double_t> GFldecays_mom_y;
  std::vector<Double_t> GFldecays_mom_z;          
  std::vector<Double_t> GFldecays_mass2;
  std::vector<Double_t> GFldecays_invbeta;
  std::vector<Int_t>    GFldecays_htofextrap;  
  std::vector<Int_t>    GFldecays_htofhitid;
  std::vector<Int_t>    GFldecays_htofseg;      
  std::vector<Double_t> GFldecays_tracklen;
  std::vector<Double_t> GFldecays_htofpos_x;
  std::vector<Double_t> GFldecays_htofpos_y;
  std::vector<Double_t> GFldecays_htofpos_z;     
  
  Double_t GFk0mass;
  Double_t GFk0decayvtx_x;
  Double_t GFk0decayvtx_y;
  Double_t GFk0decayvtx_z;
  Double_t GFk0mom;
  std::vector<Double_t> GFk0decays_id;
  std::vector<Double_t> GFk0decays_mass2;
  std::vector<Double_t> GFk0decays_invbeta;
  std::vector<Double_t> GFk0decays_mom;

  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
    trigpat.clear();
    trigflag.clear();
    
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
    pOrgTPC.clear();
    pCalcTPC.clear();
    pCorrTPC.clear();
    pCorrDETPC.clear();
    thetaTPC.clear();
    thetaCMTPC.clear();
    costCMTPC.clear();

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

    GFinside.clear();
    GFfromVtx.clear();    
    GFextrapolationHtof.clear();  
    GFntTpc_inside = 0;
    GFprodvtx_x = qnan;
    GFprodvtx_y = qnan;
    GFprodvtx_z = qnan;

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
    ppiangle = qnan;
    ldecays_id.clear();
    ldecays_mom.clear();
    ldecays_mom_x.clear();
    ldecays_mom_y.clear();
    ldecays_mom_z.clear();
    ldecays_htofhitid.clear();
    ldecays_htofseg.clear();            
    ldecays_mass2.clear();
    ldecays_tracklen.clear();
    ldecays_htofpos_x.clear();
    ldecays_htofpos_y.clear();
    ldecays_htofpos_z.clear();                        

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

    GFlflag = false;    
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
    GFldecays_htofextrap.clear();    
    GFldecays_htofhitid.clear();
    GFldecays_htofseg.clear();            
    GFldecays_mass2.clear();
    GFldecays_invbeta.clear();
    GFldecays_mom.clear();
    GFldecays_mom_x.clear();
    GFldecays_mom_y.clear();
    GFldecays_mom_z.clear();                    
    GFldecays_tracklen.clear();
    GFldecays_htofpos_x.clear();
    GFldecays_htofpos_y.clear();
    GFldecays_htofpos_z.clear();
    
    GFk0mass = qnan;
    GFk0decayvtx_x = qnan;
    GFk0decayvtx_y = qnan;
    GFk0decayvtx_z = qnan;        
    GFk0decays_id.clear();
    GFk0decays_mass2.clear();
    GFk0decays_invbeta.clear();
    GFk0decays_mom.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* runnum;
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<std::vector<Int_t>>* trigpat;
  TTreeReaderValue<std::vector<Int_t>>* trigflag;

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
  TTreeReaderValue<std::vector<Double_t>>* chisqr;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* dE;
  TTreeReaderValue<std::vector<Double_t>>* dEdx; //reference dedx
  TTreeReaderValue<std::vector<Double_t>>* mom0;//Helix momentum at Y = 0
  TTreeReaderValue<std::vector<Int_t>>* charge;//Helix charge
  TTreeReaderValue<std::vector<Double_t>>* path;//Helix path
  TTreeReaderValue<std::vector<Int_t>>* pid;

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

  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* pTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* qTPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* m2TPCKurama;
  TTreeReaderValue<std::vector<Double_t>>* xsTPC;
  TTreeReaderValue<std::vector<Double_t>>* ysTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;

  TTreeReaderValue<std::vector<Double_t>>* pK18;
  TTreeReaderValue<std::vector<Double_t>>* xbTPC;
  TTreeReaderValue<std::vector<Double_t>>* ybTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;

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

  Double_t Time0;
  Double_t CTime0;

  TTreeReaderValue<Int_t>* nhHtof;
  //  TTreeReaderValue<std::vector<Int_t>* csHtof;
  TTreeReaderValue<std::vector<Double_t>>* HtofSeg;
  TTreeReaderValue<std::vector<Double_t>>* tHtof;
  TTreeReaderValue<std::vector<Double_t>>* dtHtof;
  TTreeReaderValue<std::vector<Double_t>>* deHtof;
  TTreeReaderValue<std::vector<Double_t>>* posHtof;

  TTreeReaderValue<std::vector<Int_t>>            *isLambda;
  TTreeReaderValue<Bool_t>                        *Lflag;
  TTreeReaderValue<Double_t>                      *LambdaMass;
  //  TTreeReaderValue<std::vector<Double_t>>         *LDecaysMom;
  TTreeReaderValue<std::vector<Int_t>>            *ncombiLambda;
  TTreeReaderValue<std::vector<Double_t>>         *distLambda;
  TTreeReaderValue<std::vector<Double_t>>         *angleLambda;
  TTreeReaderValue<std::vector<Double_t>>         *bestmassLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *massLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *vtxLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *vtxLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *vtxLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *momLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *momLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *momLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *momLambda_z;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *decaysidLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *decaysmomLambda;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *decaysmomLambda_x;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *decaysmomLambda_y;
  TTreeReaderValue<std::vector<std::vector<Double_t>>> *decaysmomLambda_z;

  // TTreeReaderValue<Bool_t>* K0flag;
  // TTreeReaderValue<Double_t>* K0Mass; // Branch: "K0Mass"
  //TTreeReaderValue<Double_t>* K0DecayVtx_x; // Branch: "K0DecayVtx_x"
  // TTreeReaderValue<std::vector<Int_t>>* K0DecaysTrackId;
  // TTreeReaderValue<std::vector<Double_t>>* K0DecaysMom;

  TTreeReaderValue<Bool_t>*   lflag;
  TTreeReaderValue<Double_t>* lmass;
  TTreeReaderValue<Double_t>* ldecayvtx_x;
  TTreeReaderValue<Double_t>* ldecayvtx_y;
  TTreeReaderValue<Double_t>* ldecayvtx_z;
  TTreeReaderValue<Double_t>* lmom;
  TTreeReaderValue<Double_t>* lmom_x;
  TTreeReaderValue<Double_t>* lmom_y;
  TTreeReaderValue<Double_t>* lmom_z;
  TTreeReaderValue<Double_t>* ppi_dist;
  // TTreeReaderValue<Double_t>* ltarget_dist;
  // TTreeReaderValue<Double_t>* ltargetvtx_x;
  // TTreeReaderValue<Double_t>* ltargetvtx_y;
  // TTreeReaderValue<Double_t>* ltargetvtx_z;
  // TTreeReaderValue<Double_t>* ltargetcenter_x;
  // TTreeReaderValue<Double_t>* ltargetcenter_y;
  // TTreeReaderValue<Double_t>* ltargetcenter_z;
  // TTreeReaderValue<Double_t>* ltargetcenter_dist;
  // TTreeReaderValue<Double_t>* lprodvtx_x;
  // TTreeReaderValue<Double_t>* lprodvtx_y;
  // TTreeReaderValue<Double_t>* lprodvtx_z;
  // TTreeReaderValue<Double_t>* lprodvtx_dist;
  TTreeReaderValue<Double_t>* ltracklen;
  //TTreeReaderValue<Double_t>* ltof;
  TTreeReaderValue<std::vector<Int_t>>* ldecays_id;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_mom;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_mom_x;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_mom_y;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_mom_z;      
  TTreeReaderValue<std::vector<Int_t>>* ldecays_htofhitid;
  TTreeReaderValue<std::vector<Int_t>>* ldecays_htofseg;    
  TTreeReaderValue<std::vector<Double_t>>* ldecays_mass2;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_invbeta;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_tracklen;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_htofpos_x;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_htofpos_y;
  TTreeReaderValue<std::vector<Double_t>>* ldecays_htofpos_z;                   
  
  TTreeReaderValue<Bool_t>*         k0flag;
  TTreeReaderValue<Double_t>*       k0mass;
  TTreeReaderValue<Double_t>*       k0decayvtx_x;
  TTreeReaderValue<Double_t>*       k0decayvtx_y;
  TTreeReaderValue<Double_t>*       k0decayvtx_z;
  TTreeReaderValue<Double_t>*       k0mom_x;
  TTreeReaderValue<Double_t>*       k0mom_y;
  TTreeReaderValue<Double_t>*       k0mom_z;
  TTreeReaderValue<Double_t>*       pipi_dist;
  TTreeReaderValue<Double_t>*       pipiangle;
  TTreeReaderValue<std::vector<Int_t>>*    k0decays_id;
  TTreeReaderValue<std::vector<Double_t>>* k0decays_mom;
  TTreeReaderValue<std::vector<Double_t>>* k0decays_mom_x;
  TTreeReaderValue<std::vector<Double_t>>* k0decays_mom_y;
  TTreeReaderValue<std::vector<Double_t>>* k0decays_mom_z;

  TTreeReaderValue<Double_t>* GFK0Mass; // Branch: "GFK0Mass"
  TTreeReaderValue<Double_t>* GFK0DecayVtx_x;
  TTreeReaderValue<Double_t>* GFK0DecayVtx_y;
  TTreeReaderValue<Double_t>* GFK0DecayVtx_z;  
  TTreeReaderValue<std::vector<Double_t>>* GFk0decays_id;
  TTreeReaderValue<std::vector<Double_t>>* GFk0decays_mass2;
  TTreeReaderValue<std::vector<Double_t>>* GFk0decays_invbeta;
  TTreeReaderValue<std::vector<Double_t>>* GFk0decays_mom;
  // TTreeReaderValue<Double_t>*       GFK0decayvtx_x;
  // TTreeReaderValue<Double_t>*       GFK0decayvtx_y;
  // TTreeReaderValue<Double_t>*       GFK0decayvtx_z;
  TTreeReaderValue<Double_t>*       GFK0mom;
  TTreeReaderValue<Double_t>*       GFK0mom_x;
  TTreeReaderValue<Double_t>*       GFK0mom_y;
  TTreeReaderValue<Double_t>*       GFK0mom_z;
  TTreeReaderValue<Double_t>*       GFK0pipi_dist;
  TTreeReaderValue<Double_t>*       GFK0target_dist;
  TTreeReaderValue<Double_t>*       GFK0targetvtx_x;
  TTreeReaderValue<Double_t>*       GFK0targetvtx_y;
  TTreeReaderValue<Double_t>*       GFK0targetvtx_z;
  TTreeReaderValue<Double_t>*       GFK0targetcenter_dist;
  TTreeReaderValue<Double_t>*       GFK0targetcenter_x;
  TTreeReaderValue<Double_t>*       GFK0targetcenter_y;
  TTreeReaderValue<Double_t>*       GFK0targetcenter_z;
  TTreeReaderValue<std::vector<Int_t>>*    GFK0decays_id;
  TTreeReaderValue<std::vector<Double_t>>* GFK0decays_mass2;
  TTreeReaderValue<std::vector<Double_t>>* GFK0decays_mom;
  TTreeReaderValue<std::vector<Double_t>>* GFK0decays_invbeta;

  TTreeReaderValue<Int_t>* GFstatus;
  TTreeReaderValue<Int_t>* GFntTpc;
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
  TTreeReaderValue<std::vector<Int_t>>* GFfromVtx;  
  TTreeReaderValue<std::vector<Int_t>>* GFextrapolationHtof;

  TTreeReaderValue<Int_t>* GFntTpc_inside;
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

  TTreeReaderValue<Bool_t>*   GFlflag;
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
  TTreeReaderValue<std::vector<Int_t>>* GFldecays_id;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_mom;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_mom_x;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_mom_y;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_mom_z;      
  TTreeReaderValue<std::vector<Int_t>>* GFldecays_htofhitid;
  TTreeReaderValue<std::vector<Int_t>>* GFldecays_htofextrap;  
  TTreeReaderValue<std::vector<Int_t>>* GFldecays_htofseg;    
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_mass2;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_invbeta;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_tracklen;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_htofpos_x;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_htofpos_y;
  TTreeReaderValue<std::vector<Double_t>>* GFldecays_htofpos_z;                 
  
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
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;    

  DstClose();
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;    
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
  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;; 
  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.nKK = **src.nKK;

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


  int ntTpc = **src.ntTpc;
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

  event.isgoodTPCKurama = **src.isgoodTPCKurama;
  event.pTPCKurama = **src.pTPCKurama;
  event.qTPCKurama = **src.qTPCKurama;
  event.m2TPCKurama = **src.m2TPCKurama;
  event.xsTPC = **src.xsTPC;
  event.ysTPC = **src.ysTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

  event.pK18 = **src.pK18;
  event.xbTPC = **src.xbTPC;
  event.ybTPC = **src.ybTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;

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
  event.thetaCMTPC = **src.thetaCMTPC;
  event.costCMTPC = **src.costCMTPC;  

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
  }

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

  event.HtofSeg = **src.HtofSeg;  
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;

  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;    
  event.GFstatus = **src.GFstatus;
  event.GFntTpc = **src.GFntTpc;
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
  event.GFmom   = **src.GFmom;	
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
  event.GFfromVtx= **src.GFfromVtx;  
  event.GFextrapolationHtof = **src.GFextrapolationHtof;

  event.GFntTpc_inside = **src.GFntTpc_inside;
  event.GFprodvtx_x = **src.GFprodvtx_x;
  event.GFprodvtx_y = **src.GFprodvtx_y;
  event.GFprodvtx_z = **src.GFprodvtx_z;
  event.GFtracklen          = **src.GFtracklen;	    
  event.GFtrack2vtxdist	    = **src.GFtrack2vtxdist;    
  event.GFcalctof	    = **src.GFcalctof;	    
  event.GFsegHtof	    = **src.GFsegHtof;	    
  event.GFtofHtof	    = **src.GFtofHtof;	    
  event.GFtdiffHtof	    = **src.GFtdiffHtof;	    
  event.GFposHtof	    = **src.GFposHtof;	    
  event.GFposx	    	    = **src.GFposx;		    
  event.GFposy		    = **src.GFposy;		    
  event.GFposz		    = **src.GFposz;		    
  event.GFinvbeta	    = **src.GFinvbeta;	    
  event.GFm2	    	    = **src.GFm2;		    
  event.nsigma_tritonHtof   = **src.nsigma_tritonHtof;  
  event.nsigma_deutronHtof  = **src.nsigma_deutronHtof; 
  event.nsigma_protonHtof   = **src.nsigma_protonHtof;  
  event.nsigma_kaonHtof	    = **src.nsigma_kaonHtof;    
  event.nsigma_pionHtof     = **src.nsigma_pionHtof;    
  event.nsigma_electronHtof = **src.nsigma_electronHtof;
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
  
  event.lflag = **src.lflag;
  event.lmass = **src.LambdaMass;
  // event.ldecayvtx_x = **src.LambdaDecayVtx_x;
  // event.ldecayvtx_y = **src.LambdaDecayVtx_y;
  // event.ldecayvtx_z = **src.LambdaDecayVtx_z;
  // event.lmom_x = **src.LambdaMom_x;
  // event.lmom_y = **src.LambdaMom_y;
  // event.lmom_z = **src.LambdaMom_z;
  //event.ppi_dist = **src.LambdaVtxCloseDist;
  //  event.ppiangle = **src.LambdaPPiAngle;
  //event.ldecays_mom = **src.LDecaysMom;
  event.ldecays_mom_x = **src.ldecays_mom_x;
  event.ldecays_mom_y = **src.ldecays_mom_y;
  event.ldecays_mom_z = **src.ldecays_mom_z;
  
  event.k0flag = **src.k0flag;
  event.k0mass = **src.k0mass;
  event.k0decayvtx_x = **src.k0decayvtx_x;
  event.k0decayvtx_y = **src.k0decayvtx_y;
  event.k0decayvtx_z = **src.k0decayvtx_z;
  event.k0mom_x = **src.k0mom_x;
  event.k0mom_y = **src.k0mom_y;
  event.k0mom_z = **src.k0mom_z;
  //event.pipi_dist = **src.K0VtxCloseDist;
  event.k0decays_id = **src.k0decays_id;
  //event.k0decays_mom = **src.K0DecaysMom;
  event.k0decays_mom_x = **src.k0decays_mom_x;
  event.k0decays_mom_y = **src.k0decays_mom_y;
  event.k0decays_mom_z = **src.k0decays_mom_z;

  event.GFlflag             = **src.GFlflag             ;  
  event.GFlmass             = **src.GFlmass             ;
  event.GFldecayvtx_x       = **src.GFldecayvtx_x       ;
  event.GFldecayvtx_y       = **src.GFldecayvtx_y       ;
  event.GFldecayvtx_z       = **src.GFldecayvtx_z       ;
  event.GFlmom              = **src.GFlmom              ;
  event.GFlmom_x            = **src.GFlmom_x            ;
  event.GFlmom_y            = **src.GFlmom_y            ;
  event.GFlmom_z            = **src.GFlmom_z            ;
  event.GFppi_dist          = **src.GFppi_dist          ;
  event.GFltarget_dist      = **src.GFltarget_dist      ;
  event.GFltargetvtx_x      = **src.GFltargetvtx_x      ;
  event.GFltargetvtx_y      = **src.GFltargetvtx_y      ;
  event.GFltargetvtx_z      = **src.GFltargetvtx_z      ;
  event.GFltargetcenter_x   = **src.GFltargetcenter_x   ;
  event.GFltargetcenter_y   = **src.GFltargetcenter_y   ;
  event.GFltargetcenter_z   = **src.GFltargetcenter_z   ;
  event.GFltargetcenter_dist= **src.GFltargetcenter_dist; 
  event.GFlprodvtx_x        = **src.GFlprodvtx_x        ;
  event.GFlprodvtx_y        = **src.GFlprodvtx_y        ;
  event.GFlprodvtx_z        = **src.GFlprodvtx_z        ;
  event.GFlprodvtx_dist     = **src.GFlprodvtx_dist     ;
  event.GFltracklen         = **src.GFltracklen         ;
  event.GFltof              = **src.GFltof              ;
  event.GFldecays_id        = **src.GFldecays_id        ;
  event.GFldecays_htofhitid = **src.GFldecays_htofhitid ;
  event.GFldecays_htofextrap= **src.GFldecays_htofextrap;
  event.GFldecays_htofseg   = **src.GFldecays_htofseg   ;  
  event.GFldecays_mass2     = **src.GFldecays_mass2     ;
  event.GFldecays_invbeta   = **src.GFldecays_invbeta   ;
  event.GFldecays_tracklen  = **src.GFldecays_tracklen  ;
  event.GFldecays_invbeta   = **src.GFldecays_invbeta   ;
  event.GFldecays_mom       = **src.GFldecays_mom       ;
  event.GFldecays_mom_x     = **src.GFldecays_mom_x     ;
  event.GFldecays_mom_y     = **src.GFldecays_mom_y     ;
  event.GFldecays_mom_z     = **src.GFldecays_mom_z     ;      
  event.GFldecays_htofpos_x = **src.GFldecays_htofpos_x ;
  event.GFldecays_htofpos_y = **src.GFldecays_htofpos_y ;
  event.GFldecays_htofpos_z = **src.GFldecays_htofpos_z ;
  
  event.GFk0mass = **src.GFK0Mass;
  event.GFk0decayvtx_x = **src.GFK0DecayVtx_x;
  event.GFk0decayvtx_y = **src.GFK0DecayVtx_y;
  event.GFk0decayvtx_z = **src.GFK0DecayVtx_z;
  //event.GFk0decays_id = **src.GFk0decays_id;
  // event.GFk0decays_mass2 = **src.GFk0decays_mass2;
  // event.GFk0decays_invbeta = **src.GFk0decays_invbeta;
  // event.GFk0decays_mom = **src.GFk0decays_mom;
  HF1( 1, event.status++ );  
  HF1( 10, ntTpc );  
  Int_t GFntTpc = event.GFntTpc;
  if( event.ntTpc == 0 ) return true;
  HF1( 1, event.status++ );    
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;    
  }
  HF1( 1, event.status++ );  
  if( event.nKK != 1 ) return true;
  HF1( 1, event.status++ );    
  double BEkaon = 0.;
  double thetaTPC = 0.;  
  for(Int_t iKK=0; iKK<event.nKK; iKK++){
    BEkaon = event.MissMassNuclCorrDETPC[iKK] - KaonMass - Boron11Mass - 0.075;
    thetaTPC = event.thetaTPC[0];
  }

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);  
  std::cout << " debug " << __FILE__ << " " << __LINE__
	    << " ntTpc: " << ntTpc << std::endl;
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
    std::cout << __FILE__ << " " << __LINE__ << " GFlflag:" << event.GFlflag << std::endl;
    if(!event.GFlflag) continue;
    HF1(11204,BEkaon);
    if(thetaTPC>mint&&thetaTPC<maxt) HF1(11205,BEkaon);
    if(thetaTPC<10) HF1(11206,BEkaon);
    // if(event.GFldecays_mass2[0]<min_mass2_p||event.GFldecays_mass2[0]>max_mass2_p) continue;
    // if(event.GFldecays_mass2[1]<min_mass2_pi||event.GFldecays_mass2[1]>max_mass2_pi) continue;
    for(int it=0; it<ntTpc; it++){ // proton
      if(!event.GFfitstatus[it]) continue; 
      if(event.isElectron[it]==1) continue; 
      if(event.isK18[it]==1) continue; 
      if(event.isKurama[it]==1) continue;
      if(event.isBeam[it]==1) continue;
      if(event.isAccidental[it]==1) continue;
      if(event.charge[it]!=1) continue;
      if(it==event.GFldecays_id[0]) continue;
      if((event.pid[it]&4)!=4) continue;
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
      TVector3 lmom(event.ldecays_mom_x[0]+event.ldecays_mom_x[1],
		    event.ldecays_mom_y[0]+event.ldecays_mom_y[1],
		    event.ldecays_mom_z[0]+event.ldecays_mom_z[1]);      
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
    std::cout << __FILE__ << " " << __LINE__ << " GFlflag:" << event.GFlflag << std::endl;
    if(!event.GFlflag) continue;
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
      if(it==event.GFldecays_id[0]) continue;
      // if((event.pid[it]&4)!=4) continue;
      double nsigmap = event.nsigma_proton[it];
      if(nsigmap<-2.0||nsigmap>2.0) continue;
      if(!event.GFfromVtx[it]) continue;
      double mass2p = event.GFm2[it];
      if(mass2p<min_mass2_p||mass2p>max_mass2_p) continue; // proton found      
      HF1(11307,BEkaon);
      if(thetaTPC>mint&&thetaTPC<maxt) HF1(11308,BEkaon);
      if(thetaTPC<10) HF1(11309,BEkaon);
      std::cout << __FILE__ << " " << __LINE__ << " GFlmass:" << event.GFlmass << std::endl;
      HF1(12311,event.GFlmass);
      if(BEkaon<-0.1) HF1(12312,event.GFlmass);
      else if(BEkaon<0. ) HF1(12313,event.GFlmass);
      else if(BEkaon<0.1) HF1(12314,event.GFlmass);        
      else if(BEkaon<0.2) HF1(12315,event.GFlmass);
      else if(BEkaon<0.3) HF1(12316,event.GFlmass);
      // TVector3 lmom(event.GFldecays_mom_x[0]+event.GFldecays_mom_x[1],
      // 		    event.GFldecays_mom_y[0]+event.GFldecays_mom_y[1],
      // 		    event.GFldecays_mom_z[0]+event.GFldecays_mom_z[1]);
      TVector3 lmom(event.ldecays_mom_x[0]+event.ldecays_mom_x[1],
		    event.ldecays_mom_y[0]+event.ldecays_mom_y[1],
		    event.ldecays_mom_z[0]+event.ldecays_mom_z[1]);            
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
      if (event.charge[it1] != 1) continue;
      if ((event.pid[it1] & 4) != 4) continue; // check proton-like
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
    const double pid_threshold = 0.80; //
    for (int it1 = 0; it1 < ntTpc; ++it1) { // proton candidate
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      if(!event.GFfitstatus[it1]) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      if(event.isElectron[it1]==1) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      if(event.isK18[it1]==1) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      if(event.isKurama[it1]==1) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      if(event.isBeam[it1]==1) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      if(event.isAccidental[it1]==1) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	      
      if (event.charge[it1] != 1) continue;
      //auto prob1 = gPidLike.CalculatePosterior(event.charge[it1], event.mom0[it1], event.GFm2[it1], event.dEdx[it1]);
      //if (event.GFmom[it1].empty()) continue;
      std::cout << " debug " << __FILE__ << " " << __LINE__
		<< " GFm2: " << event.GFm2[it1] << std::endl;      std::cout << " debug " << __FILE__ << " " << __LINE__ << std::endl;                   	            
      if ( event.GFextrapolationHtof[it1]!=1 ) continue;
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
	if ( event.GFextrapolationHtof[it2]!=1 ) continue;	
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
  
  //std::vector<Double_t> GFL_mass_container(l_candidates_LH, qnan);
    
  HF1( genfitHid, GFntTpc);
  HF1( 2, event.GFstatus++);
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
  HB1( 10, "NTrack TPC", 20, 0., 20. );

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

  if(debugflag) std::cout << __FILE__ << " " << __LINE__ << std::endl;;
  
  HBTree( "tpc", "tree of GenfitPidLikelihood" );
  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
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
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  tree->Branch( "isElectron", &event.isElectron );
  tree->Branch( "nsigma_triton", &event.nsigma_triton );
  tree->Branch( "nsigma_deutron", &event.nsigma_deutron );
  tree->Branch( "nsigma_proton", &event.nsigma_proton );
  tree->Branch( "nsigma_kaon", &event.nsigma_kaon );
  tree->Branch( "nsigma_pion", &event.nsigma_pion );
  tree->Branch( "nsigma_electron", &event.nsigma_electron );
  tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );
  tree->Branch( "pid", &event.pid );

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

  tree->Branch( "isgoodTPCKurama", &event.isgoodTPCKurama);
  tree->Branch( "pTPCKurama", &event.pTPCKurama);
  tree->Branch( "qTPCKurama", &event.qTPCKurama);
  tree->Branch( "m2TPCKurama", &event.m2TPCKurama);
  tree->Branch( "xsTPC", &event.xsTPC);
  tree->Branch( "ysTPC", &event.ysTPC);
  tree->Branch( "usTPC", &event.usTPC);
  tree->Branch( "vsTPC", &event.vsTPC);
  tree->Branch( "pK18", &event.pK18);
  tree->Branch( "xbTPC", &event.xbTPC);
  tree->Branch( "ybTPC", &event.ybTPC);
  tree->Branch( "ubTPC", &event.ubTPC);
  tree->Branch( "vbTPC", &event.vbTPC);

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
  tree->Branch( "thetaTPC", &event.thetaTPC);
  tree->Branch( "thetaCMTPC", &event.thetaCMTPC);
  tree->Branch( "costCMTPC", &event.costCMTPC);


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

  //track fitting results
  tree->Branch("GFstatus", &event.GFstatus);
  tree->Branch("GFntTpc", &event.GFntTpc);
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

  tree->Branch("GFntTpc_inside", &event.GFntTpc_inside);
  tree->Branch("GFprodvtx_x", &event.GFprodvtx_x);
  tree->Branch("GFprodvtx_y", &event.GFprodvtx_y);
  tree->Branch("GFprodvtx_z", &event.GFprodvtx_z);
  
  //extrapolation
  tree->Branch("GFinside", &event.GFinside);
  tree->Branch("GFfromVtx", &event.GFfromVtx);  
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

  tree->Branch("GFLflag", &event.GFlflag);  
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
  tree->Branch("GFLambdaDecaysTrackId", &event.GFldecays_id);
  tree->Branch("GFLambdaDecaysMom", &event.GFldecays_mom);
  tree->Branch("GFLambdaDecaysMom_x", &event.GFldecays_mom_x);
  tree->Branch("GFLambdaDecaysMom_y", &event.GFldecays_mom_y);
  tree->Branch("GFLambdaDecaysMom_z", &event.GFldecays_mom_z);      
  tree->Branch("GFLambdaDecaysHtofHitId", &event.GFldecays_htofhitid);
  tree->Branch("GFLambdaDecaysHtofExtrapolate", &event.GFldecays_htofextrap);
  tree->Branch("GFLambdaDecaysHtofSeg", &event.GFldecays_htofseg);    
  tree->Branch("GFLambdaDecaysTrackLen", &event.GFldecays_tracklen);
  tree->Branch("GFLambdaDecaysInvBeta", &event.GFldecays_invbeta);
  tree->Branch("GFLambdaDecaysMass2", &event.GFldecays_mass2);     
  tree->Branch("GFLambdaDecaysHtofPos_x", &event.GFldecays_htofpos_x);
  tree->Branch("GFLambdaDecaysHtofPos_y", &event.GFldecays_htofpos_y);
  tree->Branch("GFLambdaDecaysHtofPos_z", &event.GFldecays_htofpos_z);

  
  TTreeReaderCont[kGenfitE42] = new TTreeReader( "tpc", TFileCont[kGenfitE42] );
  const auto& reader = TTreeReaderCont[kGenfitE42];
  src.runnum = new TTreeReaderValue<Int_t>( *reader, "runnum" );
  src.evnum = new TTreeReaderValue<Int_t>( *reader, "evnum" );
  src.trigpat = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigpat" );
  src.trigflag = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trigflag" );

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

  src.isgoodTPCKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isgoodTPCKurama" );
  src.pTPCKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pTPCKurama" );
  src.qTPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "qTPCKurama" );
  src.m2TPCKurama  = new TTreeReaderValue<std::vector<Double_t>>( *reader, "m2TPCKurama" );
  src.xsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xsTPC" );
  src.ysTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ysTPC" );
  src.usTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "usTPC" );
  src.vsTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vsTPC" );

  src.pK18 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "pK18" );
  src.xbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xbTPC" );
  src.ybTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ybTPC" );
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );

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
  src.thetaCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "thetaCMTPC" );
  src.costCMTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "costCMTPC" );

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

  src.isLambda          = new TTreeReaderValue<std::vector<Int_t>>            (*reader, "isLambda");
  //src.Lflag             = new TTreeReaderValue<bool>                         (*reader, "Lflag");  
  src.LambdaMass        = new TTreeReaderValue<Double_t>                      (*reader,"LambdaMass");  
  src.ncombiLambda      = new TTreeReaderValue<std::vector<Int_t>>            (*reader, "ncombiLambda");
  src.distLambda        = new TTreeReaderValue<std::vector<Double_t>>         (*reader, "distLambda");
  src.angleLambda       = new TTreeReaderValue<std::vector<Double_t>>         (*reader, "angleLambda");
  src.bestmassLambda    = new TTreeReaderValue<std::vector<Double_t>>         (*reader, "bestmassLambda");
  src.massLambda        = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "massLambda");
  src.vtxLambda_x       = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "vtxLambda_x");
  src.vtxLambda_y       = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "vtxLambda_y");
  src.vtxLambda_z       = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "vtxLambda_z");
  src.momLambda         = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "momLambda");
  src.momLambda_x       = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "momLambda_x");
  src.momLambda_y       = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "momLambda_y");
  src.momLambda_z       = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "momLambda_z");
  src.decaysidLambda    = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "decaysidLambda");
  src.decaysmomLambda   = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "decaysmomLambda");
  src.decaysmomLambda_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "decaysmomLambda_x");
  src.decaysmomLambda_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "decaysmomLambda_y");
  src.decaysmomLambda_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>> (*reader, "decaysmomLambda_z");

  src.lflag		  = new TTreeReaderValue<Bool_t>(*reader,"Lflag");			
  src.lmass		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMass");			
  src.ldecayvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_x");		
  src.ldecayvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_y");		
  src.ldecayvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_z");		
  src.lmom		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom");			
  src.lmom_x		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_x");			
  src.lmom_y		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_y");			
  src.lmom_z		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_z");			
  src.ppi_dist          = new TTreeReaderValue<Double_t>(*reader,"LambdaVtxCloseDist");		  
  src.ldecays_id        = new TTreeReaderValue<std::vector<Int_t>>(*reader,"LDecaysTrackId");
  src.ldecays_mom       = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysMom");
  src.ldecays_mom_x     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LambdaDecaysMom_x");
  src.ldecays_mom_y     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LambdaDecaysMom_y");
  src.ldecays_mom_z     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LambdaDecaysMom_z");    
  src.ldecays_htofhitid = new TTreeReaderValue<std::vector<Int_t>>(*reader,"LDecaysHtofHitId");
  src.ldecays_htofseg   = new TTreeReaderValue<std::vector<Int_t>>(*reader,"LDecaysHtofSeg");    
  src.ldecays_tracklen  = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysTrackLen");
  // src.ldecays_invbeta   = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysInvBeta");
  src.ldecays_mass2     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysMass2");
  src.ldecays_htofpos_x = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysHtofPos_x");
  src.ldecays_htofpos_y = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysHtofPos_y");
  src.ldecays_htofpos_z = new TTreeReaderValue<std::vector<Double_t>>(*reader,"LDecaysHtofPos_z");  
  
  src.k0flag        = new TTreeReaderValue<Bool_t>       (*reader, "K0flag");
  src.k0mass        = new TTreeReaderValue<Double_t>     (*reader, "K0Mass");
  src.k0decayvtx_x  = new TTreeReaderValue<Double_t>     (*reader, "K0DecayVtx_x");
  src.k0decayvtx_y  = new TTreeReaderValue<Double_t>     (*reader, "K0DecayVtx_y");
  src.k0decayvtx_z  = new TTreeReaderValue<Double_t>     (*reader, "K0DecayVtx_z");
  src.k0mom_x       = new TTreeReaderValue<Double_t>     (*reader, "K0Mom_x");
  src.k0mom_y       = new TTreeReaderValue<Double_t>     (*reader, "K0Mom_y");
  src.k0mom_z       = new TTreeReaderValue<Double_t>     (*reader, "K0Mom_z");
  src.pipi_dist     = new TTreeReaderValue<Double_t>     (*reader, "K0VtxCloseDist");
  //src.pipiangle     = new TTreeReaderValue<Double_t>     (*reader, "K0PiPiAngle");  // optional
  src.k0decays_id   = new TTreeReaderValue<std::vector<Int_t>>   (*reader, "K0DecaysTrackId");
  src.k0decays_mom  = new TTreeReaderValue<std::vector<Double_t>>(*reader, "K0DecaysMom");
  src.k0decays_mom_x= new TTreeReaderValue<std::vector<Double_t>>(*reader, "K0DecaysMom_x");
  src.k0decays_mom_y= new TTreeReaderValue<std::vector<Double_t>>(*reader, "K0DecaysMom_y");
  src.k0decays_mom_z= new TTreeReaderValue<std::vector<Double_t>>(*reader, "K0DecaysMom_z");

  src.GFK0Mass             = new TTreeReaderValue<Double_t>     (*reader, "GFK0Mass");
  src.GFK0DecayVtx_x       = new TTreeReaderValue<Double_t>     (*reader, "GFK0DecayVtx_x");
  src.GFK0DecayVtx_y       = new TTreeReaderValue<Double_t>     (*reader, "GFK0DecayVtx_y");
  src.GFK0DecayVtx_z       = new TTreeReaderValue<Double_t>     (*reader, "GFK0DecayVtx_z");
  src.GFK0mom              = new TTreeReaderValue<Double_t>     (*reader, "GFK0Mom");
  src.GFK0mom_x            = new TTreeReaderValue<Double_t>     (*reader, "GFK0Mom_x");
  src.GFK0mom_y            = new TTreeReaderValue<Double_t>     (*reader, "GFK0Mom_y");
  src.GFK0mom_z            = new TTreeReaderValue<Double_t>     (*reader, "GFK0Mom_z");
  //  src.GFK0pipi_dist        = new TTreeReaderValue<Double_t>     (*reader, "GFK0pipi_dist");
  src.GFK0target_dist      = new TTreeReaderValue<Double_t>     (*reader, "GFK0TargetCloseDist");
  src.GFK0targetvtx_x      = new TTreeReaderValue<Double_t>     (*reader, "GFK0Target_x");
  src.GFK0targetvtx_y      = new TTreeReaderValue<Double_t>     (*reader, "GFK0Target_y");
  src.GFK0targetvtx_z      = new TTreeReaderValue<Double_t>     (*reader, "GFK0Target_z");
  src.GFK0targetcenter_dist= new TTreeReaderValue<Double_t>     (*reader, "GFK0TargetCenterCloseDist");
  src.GFK0targetcenter_x   = new TTreeReaderValue<Double_t>     (*reader, "GFK0TargetCenter_x");
  src.GFK0targetcenter_y   = new TTreeReaderValue<Double_t>     (*reader, "GFK0TargetCenter_y");
  src.GFK0targetcenter_z   = new TTreeReaderValue<Double_t>     (*reader, "GFK0TargetCenter_z");
  // src.GFK0decays_id        = new TTreeReaderValue<std::vector<Int_t>>   (*reader, "GFK0decays_id");
  // src.GFK0decays_mass2     = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFK0decays_mass2");
  // src.GFK0decays_mom       = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFK0decays_mom");
  // src.GFK0decays_invbeta   = new TTreeReaderValue<std::vector<Double_t>>(*reader, "GFK0decays_invbeta");

  src.GFstatus = new TTreeReaderValue<Int_t>( *reader, "GFstatus" );
  src.GFntTpc = new TTreeReaderValue<Int_t>( *reader, "GFntTpc" );
  src.GFfitstatus = new TTreeReaderValue<std::vector<Int_t>>( *reader, "GFfitstatus" );
  src.GFpdgcode = new TTreeReaderValue<std::vector<Int_t>>( *reader, "GFpdgcode" );
  src.GFnhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "GFnhtrack" );
  src.GFcharge = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFcharge" );
  src.GFchisqr = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFchisqr" );
  src.GFtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFtof" );
  src.GFpval = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFpval" );
  src.GFlayer = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFlayer" );
  src.GFpos_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFpos_x" );
  src.GFpos_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFpos_y" );
  src.GFpos_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFpos_z" );
  src.GFmom = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFmom" );
  src.GFmom_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFmom_x" );
  src.GFmom_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFmom_y" );
  src.GFmom_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFmom_z" );  
  src.GFresidual_x = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFresidual_x" );
  src.GFresidual_y = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFresidual_y" );
  src.GFresidual_z = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFresidual_z" );
  src.GFresidual_p = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFresidual_p" );
  src.GFresidual_px = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFresidual_px" );
  src.GFresidual_py = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFresidual_py" );
  src.GFresidual_pz = new TTreeReaderValue<std::vector<std::vector<Double_t>>>( *reader, "GFresidual_pz" );
  
  src.GFinside = new TTreeReaderValue<std::vector<Int_t>>( *reader, "GFinside" );
  src.GFfromVtx= new TTreeReaderValue<std::vector<Int_t>>( *reader, "GFinside" );  
  src.GFextrapolationHtof = new TTreeReaderValue<std::vector<Int_t>>( *reader, "GFextrapolationHtof" );
  src.GFntTpc_inside = new TTreeReaderValue<Int_t>( *reader, "GFntTpc_inside" );
  src.GFprodvtx_x = new TTreeReaderValue<Double_t>( *reader, "GFprodvtx_x" );
  src.GFprodvtx_y = new TTreeReaderValue<Double_t>( *reader, "GFprodvtx_y" );
  src.GFprodvtx_z = new TTreeReaderValue<Double_t>( *reader, "GFprodvtx_z" );
  src.GFtracklen = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFtracklen" );
  src.GFtrack2vtxdist= new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFtrack2vtxdist" );
  src.GFcalctof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFcalctof" );
  src.GFsegHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFsegHtof" );
  src.GFtofHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFtofHtof" );
  src.GFtdiffHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFtdiffHtof" );
  src.GFposHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFposHtof" );
  src.GFposx = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFposx" );
  src.GFposy = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFposy" );
  src.GFposz = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFposz" );
  src.GFinvbeta = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFinvbeta" );
  src.GFm2 = new TTreeReaderValue<std::vector<Double_t>>( *reader, "GFm2" );
  src.nsigma_tritonHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_tritonHtof" );
  src.nsigma_deutronHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_deutronHtof" );
  src.nsigma_protonHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_protonHtof" );
  src.nsigma_kaonHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_kaonHtof" );
  src.nsigma_pionHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_pionHtof" );
  src.nsigma_electronHtof = new TTreeReaderValue<std::vector<Double_t>>( *reader, "nsigma_electronHtof" );
  
  src.GFlflag		  = new TTreeReaderValue<Bool_t>(*reader,"GFLflag");			
  src.GFlmass		  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaMass");			
  src.GFldecayvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaDecayVtx_x");		
  src.GFldecayvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaDecayVtx_y");		
  src.GFldecayvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaDecayVtx_z");		
  src.GFlmom		  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaMom");			
  src.GFlmom_x		  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaMom_x");			
  src.GFlmom_y		  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaMom_y");			
  src.GFlmom_z		  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaMom_z");			
  src.GFppi_dist          = new TTreeReaderValue<Double_t>(*reader,"GFLambdaVtxCloseDist");		
  src.GFltarget_dist	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTargetCloseDist");	
  src.GFltargetvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTarget_x");		
  src.GFltargetvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTarget_y");		
  src.GFltargetvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTarget_z");		
  src.GFltargetcenter_x	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTargetCenter_x");	
  src.GFltargetcenter_y	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTargetCenter_y");	
  src.GFltargetcenter_z	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTargetCenter_z");	
  src.GFltargetcenter_dist= new TTreeReaderValue<Double_t>(*reader,"GFLambdaTargetCenterCloseDist");	
  src.GFlprodvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaProductionVtx_x");	
  src.GFlprodvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaProductionVtx_y");	
  src.GFlprodvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaProductionVtx_z");	
  src.GFlprodvtx_dist	  = new TTreeReaderValue<Double_t>(*reader,"GFLambdaProductionVtxCloseDist");
  src.GFltracklen         = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTrackLen");		
  src.GFltof              = new TTreeReaderValue<Double_t>(*reader,"GFLambdaTof");                     
  src.GFldecays_id        = new TTreeReaderValue<std::vector<Int_t>>(*reader,"GFLambdaDecaysTrackId");
  src.GFldecays_mom       = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysMom");
  // src.GFldecays_mom_x     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysMom_x");
  // src.GFldecays_mom_y     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysMom_y");
  // src.GFldecays_mom_z     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysMom_z");  
  src.GFldecays_htofhitid = new TTreeReaderValue<std::vector<Int_t>>(*reader,"GFLambdaDecaysHtofHitId");
  src.GFldecays_htofextrap= new TTreeReaderValue<std::vector<Int_t>>(*reader,"GFLambdaDecaysHtofExtrapolation");
  src.GFldecays_htofseg   = new TTreeReaderValue<std::vector<Int_t>>(*reader,"GFLambdaDecaysHtofSeg");
  src.GFldecays_tracklen  = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysTrackLen");
  src.GFldecays_invbeta   = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysInvBeta");
  src.GFldecays_mass2     = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysMass2");
  src.GFldecays_htofpos_x = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysHtofPos_x");
  src.GFldecays_htofpos_y = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysHtofPos_y");
  src.GFldecays_htofpos_z = new TTreeReaderValue<std::vector<Double_t>>(*reader,"GFLambdaDecaysHtofPos_z");
  
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
