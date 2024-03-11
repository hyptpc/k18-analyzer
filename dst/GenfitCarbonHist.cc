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

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"

#define SaveRawData 0

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

//For GenFit Setting
const bool Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//0~3;
  //const Int_t verbosity = 1;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const double XiMassWindow = 0.4;
const double LambdaMassWindow = 0.1;

const Double_t p_vtx_distcut = 300;
const Double_t pi_vtx_distcut = 300;
const Double_t pi2_vtx_distcut = 300;
const Double_t xi_vtx_distcut = 100;
const Double_t vtx_scan_range = 150.;
//const Double_t vtx_scan_range = 50.;

const Double_t ppi_distcut = 10.;
const Double_t lpi_distcut = 10.;
const Double_t GFppi_distcut = 10.;
const Double_t GFlpi_distcut = 10.;

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
  std::vector<Double_t> BEcorrDE;
  std::vector<Double_t> BEcorrDETPC;
  std::vector<Double_t> km_mom_x;
  std::vector<Double_t> km_mom_y;
  std::vector<Double_t> km_mom_z;
  std::vector<Double_t> kp_mom_x;
  std::vector<Double_t> kp_mom_y;
  std::vector<Double_t> kp_mom_z;

  Int_t ntTpc; // Number of Tracks
  std::vector<Int_t> nhtrack; // Number of Hits (in 1 tracks)
  std::vector<Int_t> trackid; //for Kurama K1.8 tracks
  std::vector<Int_t> isBeam;
  std::vector<Int_t> isKurama;
  std::vector<Int_t> isK18;
  std::vector<Int_t> isAccidental;
  std::vector<Int_t> charge; //Helix charge
  std::vector<Int_t> pid;
  std::vector<Double_t> chisqr;
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

  Bool_t lpiflag;
  Bool_t lflag;

  Bool_t llflag;
  Double_t lmass1;
  Double_t ldecayvtx_x1;
  Double_t ldecayvtx_y1;
  Double_t ldecayvtx_z1;
  Double_t lmom_x1;
  Double_t lmom_y1;
  Double_t lmom_z1;
  Double_t ppi_dist1;
  Double_t lmass2;
  Double_t ldecayvtx_x2;
  Double_t ldecayvtx_y2;
  Double_t ldecayvtx_z2;
  Double_t lmom_x2;
  Double_t lmom_y2;
  Double_t lmom_z2;
  Double_t ppi_dist2;

  Double_t GFlmass1;
  Double_t GFldecayvtx_x1;
  Double_t GFldecayvtx_y1;
  Double_t GFldecayvtx_z1;
  Double_t GFlmom_x1;
  Double_t GFlmom_y1;
  Double_t GFlmom_z1;
  Double_t GFppi_dist1;
  Double_t GFlmass2;
  Double_t GFldecayvtx_x2;
  Double_t GFldecayvtx_y2;
  Double_t GFldecayvtx_z2;
  Double_t GFlmom_x2;
  Double_t GFlmom_y2;
  Double_t GFlmom_z2;
  Double_t GFppi_dist2;

  Bool_t xiflag;
  Double_t ximass;
  Double_t xidecayvtx_x;
  Double_t xidecayvtx_y;
  Double_t xidecayvtx_z;
  Double_t ximom_x;
  Double_t ximom_y;
  Double_t ximom_z;
  Double_t lpi_dist;

  Double_t lmass;
  Double_t ldecayvtx_x;
  Double_t ldecayvtx_y;
  Double_t ldecayvtx_z;
  Double_t lmom_x;
  Double_t lmom_y;
  Double_t lmom_z;
  Double_t ppi_dist;
  std::vector<Int_t> decays_id;
  std::vector<Double_t> decays_mom;
  std::vector<Double_t> decays_mom_x;
  std::vector<Double_t> decays_mom_y;
  std::vector<Double_t> decays_mom_z;

  Double_t GFximass;
  Double_t GFxidecayvtx_x;
  Double_t GFxidecayvtx_y;
  Double_t GFxidecayvtx_z;
  Double_t GFximom;
  Double_t GFximom_x;
  Double_t GFximom_y;
  Double_t GFximom_z;
  Double_t GFxiprodvtx_x;
  Double_t GFxiprodvtx_y;
  Double_t GFxiprodvtx_z;
  Double_t GFxiprodmom;
  Double_t GFxiprodmom_x;
  Double_t GFxiprodmom_y;
  Double_t GFxiprodmom_z;
  Double_t GFxitracklen;
  Double_t GFxitof;
  Double_t GFlpi_dist;

  Double_t GFlmass;
  Double_t GFldecayvtx_x;
  Double_t GFldecayvtx_y;
  Double_t GFldecayvtx_z;
  Double_t GFlmom;
  Double_t GFlmom_x;
  Double_t GFlmom_y;
  Double_t GFlmom_z;
  Double_t GFltracklen;
  Double_t GFltof;
  Double_t GFppi_dist;

  Int_t GFntdecays;
  std::vector<Double_t> GFdecays_mass;
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

  Bool_t pipiflag;

  Int_t ppi_multi;
  Int_t p_multi;
  Int_t pi_multi;
  std::vector<Double_t> ppi_mass;
  std::vector<Double_t> ppi_mom;
  std::vector<Double_t> ppi_mom_x;
  std::vector<Double_t> ppi_mom_y;
  std::vector<Double_t> ppi_mom_z;
  std::vector<Double_t> ppi_charge;

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
    m2.clear();
    xtgtKurama.clear();
    ytgtKurama.clear();
    utgtKurama.clear();
    vtgtKurama.clear();

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
    Kflag.clear();

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
    BEcorrDE.clear();
    BEcorrDETPC.clear();
    km_mom_x.clear();
    km_mom_y.clear();
    km_mom_z.clear();
    kp_mom_x.clear();
    kp_mom_y.clear();
    kp_mom_z.clear();

    ntTpc = 0;
    nhtrack.clear();
    trackid.clear();
    isBeam.clear();
    isKurama.clear();
    isK18.clear();
    isAccidental.clear();
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

    GFntTpc = 0;
    GFcharge.clear();
    GFchisqr.clear();
    GFtof.clear();
    GFtracklen.clear();
    GFpval.clear();
    GFpdgcode.clear();
    GFnhtrack.clear();

    GFntdecays = 0;
    GFdecays_mass.clear();
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
    ximom_x = qnan;
    ximom_y = qnan;
    ximom_z = qnan;
    lpi_dist = qnan;
    lmass = qnan;
    ldecayvtx_x = qnan;
    ldecayvtx_y = qnan;
    ldecayvtx_z = qnan;
    lmom_x = qnan;
    lmom_y = qnan;
    lmom_z = qnan;
    ppi_dist = qnan;

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
    GFppi_dist = qnan;

    GFxiprodvtx_x = qnan;
    GFxiprodvtx_y = qnan;
    GFxiprodvtx_z = qnan;
    GFxiprodmom = qnan;
    GFxiprodmom_x = qnan;
    GFxiprodmom_y = qnan;
    GFxiprodmom_z = qnan;
    GFxitracklen = qnan;
    GFxitof = qnan;
    GFltracklen = qnan;
    GFltof = qnan;

    lpiflag = false;
    lflag = false;
    llflag = false;
    pipiflag = false;
    lmass1 = qnan;
    ldecayvtx_x1 = qnan;
    ldecayvtx_y1 = qnan;
    ldecayvtx_z1 = qnan;
    lmom_x1 = qnan;
    lmom_y1 = qnan;
    lmom_z1 = qnan;
    ppi_dist1 = qnan;
    lmass2 = qnan;
    ldecayvtx_x2 = qnan;
    ldecayvtx_y2 = qnan;
    ldecayvtx_z2 = qnan;
    lmom_x2 = qnan;
    lmom_y2 = qnan;
    lmom_z2 = qnan;
    ppi_dist2 = qnan;

    GFlmass1 = qnan;
    GFldecayvtx_x1 = qnan;
    GFldecayvtx_y1 = qnan;
    GFldecayvtx_z1 = qnan;
    GFlmom_x1 = qnan;
    GFlmom_y1 = qnan;
    GFlmom_z1 = qnan;
    GFppi_dist1 = qnan;
    GFlmass2 = qnan;
    GFldecayvtx_x2 = qnan;
    GFldecayvtx_y2 = qnan;
    GFldecayvtx_z2 = qnan;
    GFlmom_x2 = qnan;
    GFlmom_y2 = qnan;
    GFlmom_z2 = qnan;
    GFppi_dist2 = qnan;

    decays_id.clear();
    decays_mom.clear();
    decays_mom_x.clear();
    decays_mom_y.clear();
    decays_mom_z.clear();

    ppi_multi = 0;
    p_multi = 0;
    pi_multi = 0;
    ppi_mass.clear();
    ppi_mom.clear();
    ppi_mom_x.clear();
    ppi_mom_y.clear();
    ppi_mom_z.clear();
    ppi_charge.clear();
  }
};

//_____________________________________________________________________________
struct Src
{
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

  TTreeReaderValue<Int_t>* ntTPCK18; // Number of Tracks
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
  TTreeReaderValue<std::vector<Int_t>>* isgoodTPCKurama;
  TTreeReaderValue<std::vector<Int_t>>* kflagTPCKurama;
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
  TTreeReaderValue<std::vector<Double_t>>* pOrgTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCalcTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrTPC;
  TTreeReaderValue<std::vector<Double_t>>* pCorrDETPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaTPC;
  TTreeReaderValue<std::vector<Double_t>>* thetaCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* costCMTPC;
  TTreeReaderValue<std::vector<Double_t>>* ubTPC;
  TTreeReaderValue<std::vector<Double_t>>* vbTPC;
  TTreeReaderValue<std::vector<Double_t>>* usTPC;
  TTreeReaderValue<std::vector<Double_t>>* vsTPC;

  TTreeReaderValue<Int_t>* ntTpc; // Number of Tracks
  TTreeReaderValue<std::vector<Int_t>>* nhtrack; // Number of Hits (in 1 tracks)
  TTreeReaderValue<std::vector<Int_t>>* trackid; //for Kurama K1.8 tracks
  TTreeReaderValue<std::vector<Int_t>>* isBeam;
  TTreeReaderValue<std::vector<Int_t>>* isKurama;
  TTreeReaderValue<std::vector<Int_t>>* isK18;
  TTreeReaderValue<std::vector<Int_t>>* isAccidental;
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
  skip = 224223; max_loop = 1;
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

  //if( ievent%1000==0 ){
  if( ievent%1==0 ){
  //if( ievent%100000==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  event.runnum = **src.runnum;
  event.evnum = **src.evnum;
  event.trigpat = **src.trigpat;
  event.trigflag = **src.trigflag;

  event.nKK = **src.nKK;
  event.Kflag = **src.Kflag;
  event.MissMass = **src.MissMass;
  event.MissMassCorr = **src.MissMassCorr;
  event.MissMassCorrDE = **src.MissMassCorrDE;

  event.vtx = **src.vtx;
  event.vty = **src.vty;
  event.vtz = **src.vtz;
  event.closeDist = **src.closeDist;
  event.inside = **src.inside;

  event.pKurama = **src.pKurama;
  event.xtgtKurama = **src.xtgtKurama;
  event.ytgtKurama = **src.ytgtKurama;
  event.utgtKurama = **src.utgtKurama;
  event.vtgtKurama = **src.vtgtKurama;
  event.pK18 = **src.pK18;
  event.xtgtK18 = **src.xtgtK18;
  event.ytgtK18 = **src.ytgtK18;
  event.utgtK18 = **src.utgtK18;
  event.vtgtK18 = **src.vtgtK18;

#if SaveTPCK18
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
#endif

  event.ntKurama = **src.ntKurama;
  event.isgoodTPCKurama = **src.isgoodTPCKurama;
  event.kflagTPCKurama = **src.kflagTPCKurama;
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
  event.pOrgTPC = **src.pOrgTPC;
  event.pCalcTPC = **src.pCalcTPC;
  event.pCorrTPC = **src.pCorrTPC;
  event.pCorrDETPC = **src.pCorrDETPC;
  event.thetaTPC = **src.thetaTPC;
  event.thetaCMTPC = **src.thetaCMTPC;
  event.costCMTPC = **src.costCMTPC;
  event.ubTPC = **src.ubTPC;
  event.vbTPC = **src.vbTPC;
  event.usTPC = **src.usTPC;
  event.vsTPC = **src.vsTPC;

  event.BE.resize(event.nKK);
  event.BEcorrDE.resize(event.nKK);
  event.BETPC.resize(event.nKK);
  event.BEcorrDETPC.resize(event.nKK);

  event.km_mom_x.resize(event.nKK);
  event.km_mom_y.resize(event.nKK);
  event.km_mom_z.resize(event.nKK);
  event.kp_mom_x.resize(event.nKK);
  event.kp_mom_y.resize(event.nKK);
  event.kp_mom_z.resize(event.nKK);

  if( event.nKK != 1 || event.Kflag[0]!=1 ) return false;
  if( event.inside[0] != 1) std::cout<<"not inside!!!!!!!!!"<<std::endl;

  Double_t pKp = event.pKurama[0];
  Double_t xKp = event.xtgtKurama[0];
  Double_t yKp = event.ytgtKurama[0];
  Double_t uKp = event.utgtKurama[0];
  Double_t vKp = event.vtgtKurama[0];
  Double_t ptKp = pKp/std::sqrt(1.+uKp*uKp+vKp*vKp);
  TVector3 kp_mom(ptKp*uKp, ptKp*vKp, ptKp);

  Double_t pKm = event.pK18[0];
  Double_t xKm = event.xtgtK18[0];
  Double_t yKm = event.ytgtK18[0];
  Double_t uKm = event.utgtK18[0];
  Double_t vKm = event.vtgtK18[0];
  Double_t ptKm = pKm/std::sqrt(1.+uKm*uKm+vKm*vKm);
  TVector3 km_mom(ptKm*uKm, ptKm*vKm, ptKm);
  Double_t cost = km_mom*kp_mom/(pKm*pKp);
  Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
  TLorentzVector LvKm(km_mom, TMath::Hypot(km_mom.Mag(), KaonMass));
  TLorentzVector LvKp(kp_mom, TMath::Hypot(kp_mom.Mag(), KaonMass));

  TLorentzVector LvC(0., 0., 0., m12C);
  TLorentzVector LvP(0., 0., 0., ProtonMass);
  TLorentzVector LvRproton = LvKm + LvP - LvKp;
  TLorentzVector LvRc = LvKm + LvC - LvKp;
  Double_t mm_12C = LvRc.Mag();
  Double_t binding_energy = m11B + XiMinusMass - mm_12C; //GeV/c2

  event.km_mom_x[0] = km_mom.x();
  event.km_mom_y[0] = km_mom.y();
  event.km_mom_z[0] = km_mom.z();
  event.kp_mom_x[0] = kp_mom.x();
  event.kp_mom_y[0] = kp_mom.y();
  event.kp_mom_z[0] = kp_mom.z();

  event.BE[0] = 1000.*binding_energy; //MeV/c2
  HF1( 100, -event.BE[0]);
  //event.BEcorrDE = ;
  //HB1( 101, event.BEcorrDE[0]);

  if(event.isgoodTPC[0] == 1){
    TVector3 km_unit = TVector3(event.ubTPC[0], event.vbTPC[0], 1.).Unit();
    TVector3 km_momTPC = km_unit*event.pK18[0];

    TVector3 kp_unit = TVector3(event.usTPC[0], event.vsTPC[0], 1.).Unit();
    TVector3 kp_momTPC = kp_unit*event.pCorrTPC[0];
    Double_t thetaTPC = event.thetaTPC[0];

    TLorentzVector LvKmTPC(km_momTPC, TMath::Hypot(km_momTPC.Mag(), KaonMass));
    TLorentzVector LvKpTPC(kp_momTPC, TMath::Hypot(kp_momTPC.Mag(), KaonMass));
    TLorentzVector LvCTPC(0., 0., 0., m12C);
    TLorentzVector LvPTPC(0., 0., 0., ProtonMass);
    TLorentzVector LvRprotonTPC = LvKmTPC + LvPTPC - LvKpTPC;
    TLorentzVector LvRcTPC = LvKmTPC + LvCTPC - LvKpTPC;

    double mm_12CTPC = LvRcTPC.M();
    double binding_energyTPC = m11B + XiMinusMass - mm_12CTPC; //GeV/c2
    event.BETPC[0] = 1000.*binding_energyTPC; //MeV/c2
    HF1( 102, -event.BETPC[0]);
    //event.BEcorrDETPC[0] = ;
    //HB1( 103, event.BEcorrDETPC[0]);
  }

  HF1( 1, event.status++ );

  event.nhHtof = **src.nhHtof;
  event.HtofSeg = **src.HtofSeg;
  event.tHtof = **src.tHtof;
  event.dtHtof = **src.dtHtof;
  event.deHtof = **src.deHtof;
  event.posHtof = **src.posHtof;
  Int_t ntTpc = **src.ntTpc;
  if( ntTpc == 0 ) return true;
  HF1( 1, event.status++ );

  event.ntTpc = ntTpc;
  event.nhtrack = **src.nhtrack;
  event.trackid = **src.trackid;
  event.isBeam = **src.isBeam;
  event.isKurama = **src.isKurama;
  event.isK18 = **src.isK18;
  event.isAccidental = **src.isAccidental;
  event.charge = **src.charge;
  event.pid = **src.pid;
  event.chisqr = **src.chisqr;
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
  event.pathhit = **src.pathhit;
  event.alpha = **src.alpha;
  event.track_cluster_de = **src.track_cluster_de;
  event.track_cluster_mrow = **src.track_cluster_mrow;

  Double_t dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  TPCAnalyzer TPCAna;
  TPCAna.ReCalcTPCTracks(**src.ntTpc, **src.isK18, **src.isKurama,
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

    int dedxPID = Kinematics::HypTPCdEdxPID_temp(event.dEdx[it], event.mom0[it]*event.charge[it]);
    event.pid[it]=dedxPID;

    std::vector<Int_t> pdgcode;
    Kinematics::HypTPCPID_PDGCode(event.charge[it], dedxPID, pdgcode);

    GFTrackCont.AddHelixTrack(pdgcode, tp);
  }
  GFTrackCont.FitTracks();

  Int_t GFntTpc = GFTrackCont.GetNTrack();
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    return true;
  }

  std::vector<Int_t> target_p_container, target_pip_container, target_pim_container;
  std::vector<Double_t> target_p_mass_container, target_pip_mass_container, target_pim_mass_container;
  std::vector<TVector3> target_p_mom_container, target_pip_mom_container, target_pim_mom_container;

  for(Int_t it=0;it<ntTpc;it++){
    if(event.isK18[it]==1 || event.isKurama[it]==1) continue;
    if((event.pid[it]&4)==4 && event.charge[it]==1){
      Int_t repid_p = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it, repid_p)) continue;
      if(!GFTrackCont.IsInsideTarget(it)) continue;

      TVector3 mom = GFTrackCont.GetMom(it, 0, repid_p);
      target_p_container.push_back(it);
      target_p_mass_container.push_back(qnan);
      target_p_mom_container.push_back(mom);
    }
    else if((event.pid[it]&1)==1){
      Int_t repid_pi = 0;
      if(!GFTrackCont.TrackCheck(it, repid_pi)) continue;
      if(!GFTrackCont.IsInsideTarget(it)) continue;

      Int_t htofhitid_pi; Double_t tracklen_pi;
      Double_t tof; TVector3 pos; TVector3 vtx;
      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it, repid_pi, event.HtofSeg,
								     event.posHtof, htofhitid_pi,
								     tof, tracklen_pi, pos, vtx);
      Double_t mass = qnan;
      TVector3 mom = GFTrackCont.GetMom(it, 0, repid_pi);
      if(htofextrapolation_pi) mass = Kinematics::MassSquare(mom.Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
      if(htofextrapolation_pi && mass > 0.25) continue;
      if(event.charge[it]==1){
	target_pip_container.push_back(it);
	target_pip_mass_container.push_back(mass);
	target_pip_mom_container.push_back(mom);
      }
      else{
	target_pim_container.push_back(it);
	target_pim_mass_container.push_back(mass);
	target_pim_mom_container.push_back(mom);
      }
    }
  }

  //for L, LL
  std::vector<Int_t> L_p_container, L_pi_container;
  std::vector<Int_t> L_p_repid_container, L_pi_repid_container;
  std::vector<Double_t> L_mass_container;
  std::vector<Double_t> L_ppidist_container;
  std::vector<Double_t> L_p_mass_container, L_pi_mass_container;
  std::vector<TVector3> L_mom_container, L_vtx_container;
  std::vector<TVector3> L_p_mom_container, L_pi_mom_container;

  //L + pi
  std::vector<Int_t> Lpi_p_container, Lpi_pi_container, Lpi_pi2_container;

  //Xi candidates searching
  std::cout<<"Xi/L Searching "<<std::endl;
  Int_t xi_candidates = 0;
  {
    std::vector<Int_t> xi_p_container, xi_pi_container, xi_pi2_container;
    std::vector<Int_t> p_repid_container, pi_repid_container, pi2_repid_container;
    std::vector<TVector3> xi_mom_container, xi_vert_container;
    std::vector<TVector3> l_mom_container, l_vert_container;
    std::vector<Double_t> xi_mass_container, lambda_mass_container;
    std::vector<TVector3> p_mom_container, pi_mom_container, pi2_mom_container;
    std::vector<Double_t> p_mass_container, pi_mass_container, pi2_mass_container;
    std::vector<Double_t> ppi_closedist; std::vector<Double_t> lpi_closedist;

    for(Int_t it1=0;it1<ntTpc;it1++){
      if(event.isK18[it1]==1 || event.isKurama[it1]==1) continue;
      //select proton like
      if((event.pid[it1]&4)!=4 || event.charge[it1]!=1) continue;
      Int_t repid_p = 0;
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	Int_t temp = flag&event.pid[it1];
	if(temp==flag) repid_p += 1;
	flag*=2;
      }
      if(!GFTrackCont.TrackCheck(it1, repid_p)) continue;

      Double_t GFmass_decays[3];
      Int_t htofhitid_p; Double_t tracklen_p;
      Double_t tof; TVector3 pos; TVector3 vtx;
      Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(it1, repid_p,
								    event.HtofSeg, event.posHtof,
								    htofhitid_p, tof, tracklen_p,
								    pos, vtx);
      //if(htofextrapolation_p){
      if(htofextrapolation_p&&false){
	GFmass_decays[0] = Kinematics::MassSquare(GFTrackCont.GetMom(it1, 0, repid_p).Mag(), tracklen_p, event.tHtof[htofhitid_p]);
	if(GFmass_decays[0] < 0.25) continue;
      }
      else GFmass_decays[0] = qnan;

      Double_t p_par[5];
      p_par[0] = event.helix_cx[it1];
      p_par[1] = event.helix_cy[it1];
      p_par[2] = event.helix_z0[it1];
      p_par[3] = event.helix_r[it1];
      p_par[4] = event.helix_dz[it1];
      Int_t p_nh = event.helix_t[it1].size();
      Double_t p_theta_min = event.helix_t[it1][0] - vtx_scan_range/p_par[3];
      Double_t p_theta_max = event.helix_t[it1][0];
      TVector3 p_start = TVector3(event.calpos_x[it1][0], event.calpos_y[it1][0], event.calpos_z[it1][0]);
      TVector3 p_end = TVector3(event.calpos_x[it1][p_nh-1], event.calpos_y[it1][p_nh-1], event.calpos_z[it1][p_nh-1]);

      for(Int_t it2=0;it2<ntTpc;it2++){
	if(it1==it2) continue;
	if(event.isK18[it2]==1 || event.isKurama[it2]==1) continue;
	//select pi- like
	if((event.pid[it2]&1)!=1) continue;
	if(event.charge[it2]!=-1) continue;
	Int_t repid_pi = 0;
	if(!GFTrackCont.TrackCheck(it2, repid_pi)) continue;
	Int_t htofhitid_pi; Double_t tracklen_pi;
	Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(it2, repid_pi,
								    event.HtofSeg, event.posHtof, htofhitid_pi,
								    tof, tracklen_pi, pos, vtx);
	if(htofextrapolation_pi){
	  GFmass_decays[1] = Kinematics::MassSquare(GFTrackCont.GetMom(it2, 0, repid_pi).Mag(), tracklen_pi, event.tHtof[htofhitid_pi]);
	  if(GFmass_decays[1] > 0.25) continue;
	}
	else GFmass_decays[1] = qnan;

	Double_t pi_par[5];
	pi_par[0] = event.helix_cx[it2];
	pi_par[1] = event.helix_cy[it2];
	pi_par[2] = event.helix_z0[it2];
	pi_par[3] = event.helix_r[it2];
	pi_par[4] = event.helix_dz[it2];

	Int_t pi_nh = event.helix_t[it2].size();
	Double_t pi_theta_min = event.helix_t[it2][0];
	Double_t pi_theta_max = event.helix_t[it2][0] + vtx_scan_range/pi_par[3];
	TVector3 pi_start = TVector3(event.calpos_x[it2][0], event.calpos_y[it2][0], event.calpos_z[it2][0]);
	TVector3 pi_end = TVector3(event.calpos_x[it2][pi_nh-1], event.calpos_y[it2][pi_nh-1], event.calpos_z[it2][pi_nh-1]);

	Double_t ppi_dist = 10000.;
	TVector3 p_mom; TVector3 pi_mom; TVector3 lambda_mom;
	TVector3 lambda_vert = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par, p_theta_min, p_theta_max, pi_theta_min, pi_theta_max, p_mom, pi_mom, lambda_mom, ppi_dist);
	TLorentzVector Lpi(pi_mom, TMath::Sqrt(pi_mom.Mag()*pi_mom.Mag() + PionMass*PionMass));
	TLorentzVector Lp(p_mom, TMath::Sqrt(p_mom.Mag()*p_mom.Mag() + ProtonMass*ProtonMass));
	TLorentzVector Llambda = Lp + Lpi;

	if(TMath::Abs(lambda_vert.x()) > 250. || TMath::Abs(lambda_vert.z()) > 250. || TMath::Abs(lambda_vert.y()) > 250.) continue; //Vertex cut
	if(ppi_dist > ppi_distcut) continue;
	Double_t pi_vertex_dist; Double_t p_vertex_dist;
	if(!Kinematics::HelixDirection(lambda_vert, p_start, p_end, p_vertex_dist) ||
	   !Kinematics::HelixDirection(lambda_vert, pi_start, pi_end, pi_vertex_dist)) continue;

	if(pi_vertex_dist > pi_vtx_distcut) continue;
	if(p_vertex_dist > p_vtx_distcut) continue;
	if(TMath::Abs(Llambda.M() - LambdaMass) > LambdaMassWindow*0.5) continue;

	L_p_container.push_back(it1);
	L_pi_container.push_back(it2);
	L_mass_container.push_back(Llambda.M());
	L_mom_container.push_back(lambda_mom);
	L_p_mom_container.push_back(p_mom);
	L_pi_mom_container.push_back(pi_mom);
	L_p_mass_container.push_back(GFmass_decays[0]);
	L_pi_mass_container.push_back(GFmass_decays[1]);
	L_ppidist_container.push_back(ppi_dist);
	L_vtx_container.push_back(lambda_vert);
	L_p_repid_container.push_back(repid_p);
	L_pi_repid_container.push_back(repid_pi);
	for(int it3=0;it3<ntTpc;it3++){
	  if(it3==it2 || it3==it1) continue;
	  if(event.isK18[it3]==1 || event.isKurama[it3]==1) continue;
	  //select pi like
	  if((event.pid[it3]&1)!=1) continue;
	  if(event.charge[it3]!=-1) continue;
	  Int_t repid_pi2 = 0;
	  if(!GFTrackCont.TrackCheck(it3, repid_pi2)) continue;
	  Int_t htofhitid_pi2; Double_t tracklen_pi2;
	  Bool_t htofextrapolation_pi2 = GFTrackCont.TPCHTOFTrackMatching(it3, repid_pi2,
								       event.HtofSeg, event.posHtof, htofhitid_pi2,
								       tof, tracklen_pi2, pos, vtx);
	  if(htofextrapolation_pi2){
	    GFmass_decays[2] = Kinematics::MassSquare(GFTrackCont.GetMom(it3, 0, repid_pi2).Mag(), tracklen_pi2, event.tHtof[htofhitid_pi2]);
	    if(GFmass_decays[2] > 0.25) continue;

	  }
	  else GFmass_decays[2] = qnan;

	  Double_t pi2_par[5];
	  pi2_par[0] = event.helix_cx[it3];
	  pi2_par[1] = event.helix_cy[it3];
	  pi2_par[2] = event.helix_z0[it3];
	  pi2_par[3] = event.helix_r[it3];
	  pi2_par[4] = event.helix_dz[it3];

	  Int_t pi2_nh = event.helix_t[it3].size();
	  Double_t pi2_theta_min = event.helix_t[it3][0];
	  Double_t pi2_theta_max = event.helix_t[it3][0] + vtx_scan_range/pi2_par[3];
	  TVector3 pi2_start = TVector3(event.calpos_x[it3][0], event.calpos_y[it3][0], event.calpos_z[it3][0]);
	  TVector3 pi2_end = TVector3(event.calpos_x[it3][pi2_nh-1], event.calpos_y[it3][pi2_nh-1], event.calpos_z[it3][pi2_nh-1]);

	  TVector3 pi2_mom; Double_t lpi_dist;
	  TVector3 xi_vert = Kinematics::XiVertex(dMagneticField, pi2_par, pi2_theta_min, pi2_theta_max, lambda_vert, lambda_mom, pi2_mom, lpi_dist);

	  TLorentzVector Lpi2(pi2_mom, TMath::Sqrt(pi2_mom.Mag()*pi2_mom.Mag() + PionMass*PionMass));
	  TLorentzVector Llambda_fixedmass(lambda_mom, TMath::Sqrt(lambda_mom.Mag()*lambda_mom.Mag() + LambdaMass*LambdaMass));
	  TLorentzVector Lxi = Llambda_fixedmass + Lpi2;
	  TVector3 xi_mom = TVector3(Lxi.Px(), Lxi.Py(), Lxi.Pz());
	  Double_t xi_vertex_dist; Double_t pi2_vertex_dist;

	  Lpi_p_container.push_back(it1);
	  Lpi_pi_container.push_back(it2);
	  Lpi_pi2_container.push_back(it3);

	  if(lpi_dist > lpi_distcut) continue;
	  if(!Kinematics::HelixDirection(xi_vert, p_start, p_end, xi_vertex_dist) ||
	     !Kinematics::HelixDirection(xi_vert, pi2_start, pi2_end, pi2_vertex_dist)) continue;
	  if(pi2_vertex_dist > pi2_vtx_distcut) continue;
	  if(xi_vertex_dist > xi_vtx_distcut) continue;
	  if(TMath::Abs(Lxi.M() - XiMinusMass) > 0.5*XiMassWindow) continue;

	  ppi_closedist.push_back(ppi_dist);
	  lpi_closedist.push_back(lpi_dist);

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
	  p_mom_container.push_back(p_mom);
	  pi_mom_container.push_back(pi_mom);
	  pi2_mom_container.push_back(pi2_mom);

	  p_mass_container.push_back(GFmass_decays[0]);
	  pi_mass_container.push_back(GFmass_decays[1]);
	  pi2_mass_container.push_back(GFmass_decays[2]);

	  xi_vert_container.push_back(xi_vert);
	  l_vert_container.push_back(lambda_vert);

	  xi_candidates++;
	} //it3
      } //it2
    } //it1
    std::cout<<"GF Xi Searching ends"<<std::endl;

    std::cout<<"LL Searching starts"<<std::endl;
    if(L_p_container.size()>0){
      bool double_L = false;
      int bestid_p[2] = {-1, -1}; int bestid_pi[2] = {-1, -1}; int bestid[2] = {-1, -1};
      Double_t GFextrapdist_decays[4];
      TVector3 GFmom_decays[4];
      TVector3 GFlambda_vert1; double GFppi_dist1 = qnan;
      TVector3 GFlambda_vert2; double GFppi_dist2 = qnan;

      for(int idp1=0;idp1<L_p_container.size();idp1++){
	int count = 0;
	bestid_p[0] = L_p_container[idp1];
	bestid_pi[0] = L_pi_container[idp1];

	Bool_t vtxcut1 = (GFTrackCont.FindVertex(bestid_p[0], bestid_pi[0],
						 L_p_repid_container[idp1], L_pi_repid_container[idp1],
						 GFextrapdist_decays[0], GFextrapdist_decays[1],
						 GFmom_decays[0], GFmom_decays[1],
						 GFppi_dist1, GFlambda_vert1,
						 vtx_scan_range) && GFppi_dist1 < GFppi_distcut);

	if(!vtxcut1) continue;
	for(int idp2=idp1+1;idp2<L_p_container.size();idp2++){
	  if(L_p_container[idp1]==L_p_container[idp2] || L_p_container[idp1]==L_pi_container[idp2]) continue;
	  if(L_pi_container[idp1]==L_p_container[idp2] || L_pi_container[idp1]==L_pi_container[idp2]) continue;

	  bestid_p[1] = L_p_container[idp2];
	  bestid_pi[1] = L_pi_container[idp2];
	  Bool_t vtxcut2 = (GFTrackCont.FindVertex(bestid_p[1], bestid_pi[1],
						   L_p_repid_container[idp2], L_pi_repid_container[idp2],
						   GFextrapdist_decays[2], GFextrapdist_decays[3],
						   GFmom_decays[2], GFmom_decays[3],
						   GFppi_dist2, GFlambda_vert2,
						   vtx_scan_range) && GFppi_dist2 < GFppi_distcut);
	  if(!vtxcut2) continue;
	  double_L = true;
	  bestid[0] = idp1;
	  bestid[1] = idp2;

	  TLorentzVector GFLp1(GFmom_decays[0],
			       TMath::Sqrt(GFmom_decays[0].Mag()*GFmom_decays[0].Mag() + ProtonMass*ProtonMass));
	  TLorentzVector GFLpi1(GFmom_decays[1],
				TMath::Sqrt(GFmom_decays[1].Mag()*GFmom_decays[1].Mag() + PionMass*PionMass));
	  TLorentzVector GFLp2(GFmom_decays[2],
			       TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + ProtonMass*ProtonMass));
	  TLorentzVector GFLpi2(GFmom_decays[3],
				TMath::Sqrt(GFmom_decays[3].Mag()*GFmom_decays[3].Mag() + PionMass*PionMass));
	  TLorentzVector GFLlambda1 = GFLp1 + GFLpi1;
	  TLorentzVector GFLlambda2 = GFLp2 + GFLpi2;
	  TVector3 GFlambda_mom1 = GFmom_decays[0] + GFmom_decays[1];
	  TVector3 GFlambda_mom2 = GFmom_decays[2] + GFmom_decays[3];

	  event.llflag = true;
	  event.lmass1 = L_mass_container.at(bestid[0]);
	  event.ldecayvtx_x1 = L_vtx_container.at(bestid[0]).x();
	  event.ldecayvtx_y1 = L_vtx_container.at(bestid[0]).y();
	  event.ldecayvtx_z1 = L_vtx_container.at(bestid[0]).z();
	  event.lmom_x1 = L_mom_container.at(bestid[0]).x();
	  event.lmom_y1 = L_mom_container.at(bestid[0]).y();
	  event.lmom_z1 = L_mom_container.at(bestid[0]).z();
	  event.ppi_dist1 = L_ppidist_container.at(bestid[0]);
	  event.lmass2 = L_mass_container.at(bestid[1]);
	  event.ldecayvtx_x2 = L_vtx_container.at(bestid[1]).x();
	  event.ldecayvtx_y2 = L_vtx_container.at(bestid[1]).y();
	  event.ldecayvtx_z2 = L_vtx_container.at(bestid[1]).z();
	  event.lmom_x2 = L_mom_container.at(bestid[1]).x();
	  event.lmom_y2 = L_mom_container.at(bestid[1]).y();
	  event.lmom_z2 = L_mom_container.at(bestid[1]).z();
	  event.ppi_dist2 = L_ppidist_container.at(bestid[1]);

	  event.GFlmass1 = GFLlambda1.M();
	  event.GFldecayvtx_x1 = GFlambda_vert1.x();
	  event.GFldecayvtx_y1 = GFlambda_vert1.y();
	  event.GFldecayvtx_z1 = GFlambda_vert1.z();
	  event.GFlmom_x1 = GFlambda_mom1.x();
	  event.GFlmom_y1 = GFlambda_mom1.y();
	  event.GFlmom_z1 = GFlambda_mom1.z();
	  event.GFppi_dist1 = GFppi_dist1;
	  event.GFlmass2 = GFLlambda2.M();
	  event.GFldecayvtx_x2 = GFlambda_vert2.x();
	  event.GFldecayvtx_y2 = GFlambda_vert2.y();
	  event.GFldecayvtx_z2 = GFlambda_vert2.z();
	  event.GFlmom_x2 = GFlambda_mom2.x();
	  event.GFlmom_y2 = GFlambda_mom2.y();
	  event.GFlmom_z2 = GFlambda_mom2.z();
	  event.GFppi_dist2 = GFppi_dist2;

	  event.GFntdecays = 4;
	  event.decays_mom.push_back(L_p_mom_container[bestid[0]].Mag());
	  event.decays_mom.push_back(L_pi_mom_container[bestid[0]].Mag());
	  event.decays_mom.push_back(L_p_mom_container[bestid[1]].Mag());
	  event.decays_mom.push_back(L_pi_mom_container[bestid[1]].Mag());
	  event.decays_mom_x.push_back(L_p_mom_container[bestid[0]].x());
	  event.decays_mom_x.push_back(L_pi_mom_container[bestid[0]].x());
	  event.decays_mom_x.push_back(L_p_mom_container[bestid[1]].x());
	  event.decays_mom_x.push_back(L_pi_mom_container[bestid[1]].x());
	  event.decays_mom_y.push_back(L_p_mom_container[bestid[0]].y());
	  event.decays_mom_y.push_back(L_pi_mom_container[bestid[0]].y());
	  event.decays_mom_y.push_back(L_p_mom_container[bestid[1]].y());
	  event.decays_mom_y.push_back(L_pi_mom_container[bestid[1]].y());
	  event.decays_mom_z.push_back(L_p_mom_container[bestid[0]].z());
	  event.decays_mom_z.push_back(L_pi_mom_container[bestid[0]].z());
	  event.decays_mom_z.push_back(L_p_mom_container[bestid[1]].z());
	  event.decays_mom_z.push_back(L_pi_mom_container[bestid[1]].z());
	  event.decays_id.push_back(bestid_p[0]);
	  event.decays_id.push_back(bestid_pi[0]);
	  event.decays_id.push_back(bestid_p[1]);
	  event.decays_id.push_back(bestid_pi[1]);

	  event.GFdecays_mass.push_back(L_p_mass_container[bestid[0]]);
	  event.GFdecays_mass.push_back(L_pi_mass_container[bestid[0]]);
	  event.GFdecays_mass.push_back(L_p_mass_container[bestid[1]]);
	  event.GFdecays_mass.push_back(L_pi_mass_container[bestid[1]]);
	  event.GFdecays_mom.push_back(GFmom_decays[0].Mag());
	  event.GFdecays_mom.push_back(GFmom_decays[1].Mag());
	  event.GFdecays_mom.push_back(GFmom_decays[2].Mag());
	  event.GFdecays_mom.push_back(GFmom_decays[3].Mag());
	  event.GFdecays_mom_x.push_back(GFmom_decays[0].x());
	  event.GFdecays_mom_x.push_back(GFmom_decays[1].x());
	  event.GFdecays_mom_x.push_back(GFmom_decays[2].x());
	  event.GFdecays_mom_x.push_back(GFmom_decays[3].x());
	  event.GFdecays_mom_y.push_back(GFmom_decays[0].y());
	  event.GFdecays_mom_y.push_back(GFmom_decays[1].y());
	  event.GFdecays_mom_y.push_back(GFmom_decays[2].y());
	  event.GFdecays_mom_y.push_back(GFmom_decays[3].y());
	  event.GFdecays_mom_z.push_back(GFmom_decays[0].z());
	  event.GFdecays_mom_z.push_back(GFmom_decays[1].z());
	  event.GFdecays_mom_z.push_back(GFmom_decays[2].z());
	  event.GFdecays_mom_z.push_back(GFmom_decays[3].z());

	  HF2( 502, event.lmass1, event.lmass2);
	  HF2( 503, event.GFlmass1, event.GFlmass2);
	  HF1( 504, event.lmass1);
	  HF1( 504, event.lmass2);
	  HF1( 505, event.GFlmass1);
	  HF1( 505, event.GFlmass2);

	  if(double_L) break;
	}
	if(double_L) break;
      }
      if(double_L){
	event.p_multi = target_p_container.size();
	event.pi_multi = target_pip_container.size() + target_pim_container.size();
	for(Int_t itp=0; itp<target_p_container.size(); ++itp){
	  if(itp==bestid_p[0] || itp==bestid_p[1]){
	    event.p_multi -= 1;
	    continue;
	  }
	  event.ppi_mass.push_back(target_p_mass_container[itp]);
	  event.ppi_mom.push_back(target_p_mom_container[itp].Mag());
	  event.ppi_mom_x.push_back(target_p_mom_container[itp].x());
	  event.ppi_mom_y.push_back(target_p_mom_container[itp].y());
	  event.ppi_mom_z.push_back(target_p_mom_container[itp].z());
	  event.ppi_charge.push_back(1);
	  HF1( 1110, target_p_mom_container[itp].Mag());
	}
	for(Int_t itpip=0; itpip<target_pip_container.size(); ++itpip){
	  event.ppi_mass.push_back(target_pip_mass_container[itpip]);
	  event.ppi_mom.push_back(target_pip_mom_container[itpip].Mag());
	  event.ppi_mom_x.push_back(target_pip_mom_container[itpip].x());
	  event.ppi_mom_y.push_back(target_pip_mom_container[itpip].y());
	  event.ppi_mom_z.push_back(target_pip_mom_container[itpip].z());
	  event.ppi_charge.push_back(1);
	  HF1( 1111, target_pip_mom_container[itpip].Mag());
	}
	for(Int_t itpim=0; itpim<target_pim_container.size(); ++itpim){
	  if(itpim==bestid_pi[0] || itpim==bestid_pi[1]){
	    event.pi_multi -= 1;
	    continue;
	  }
	  event.ppi_mass.push_back(target_pim_mass_container[itpim]);
	  event.ppi_mom.push_back(target_pim_mom_container[itpim].Mag());
	  event.ppi_mom_x.push_back(target_pim_mom_container[itpim].x());
	  event.ppi_mom_y.push_back(target_pim_mom_container[itpim].y());
	  event.ppi_mom_z.push_back(target_pim_mom_container[itpim].z());
	  event.ppi_charge.push_back(-1);
	  HF1( 1112, target_pim_mom_container[itpim].Mag());
	}
	event.ppi_multi = event.p_multi + event.pi_multi;

	HF1( 1010, event.ppi_multi);
	HF1( 1011, event.p_multi);
	HF1( 1012, event.pi_multi);

	HF1( 110, -event.BE[0]);
	if(event.isgoodTPC[0] == 1) HF1( 112, -event.BETPC[0]);
	return true;
      }
    } //L_p_container

    std::vector<Int_t> GFxi_p_id_container(xi_candidates, -1);
    std::vector<Int_t> GFxi_pi_id_container(xi_candidates, -1);
    std::vector<Int_t> GFxi_pi2_id_container(xi_candidates, -1);
    std::vector<Int_t> GFxi_p_rep_container(xi_candidates, -1);
    std::vector<Int_t> GFxi_pi_rep_container(xi_candidates, -1);
    std::vector<Int_t> GFxi_pi2_rep_container(xi_candidates, -1);
    std::vector<TVector3> GFxi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFl_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFp_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFpi_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFpi2_mom_container(xi_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFxi_vert_container(xi_candidates, TVector3(qnan, qnan, qnan));
    std::vector<TVector3> GFl_vert_container(xi_candidates, TVector3(qnan, qnan, qnan));
    std::vector<Double_t> GFxi_mass_container(xi_candidates, qnan);
    std::vector<Double_t> GFl_mass_container(xi_candidates, qnan);
    std::vector<Double_t> GFxi_tracklen_container(xi_candidates, qnan);
    std::vector<Double_t> GFl_tracklen_container(xi_candidates, qnan);
    std::vector<Double_t> GFxi_tof_container(xi_candidates, qnan);
    std::vector<Double_t> GFl_tof_container(xi_candidates, qnan);
    std::vector<Double_t> GFp_mass_container(xi_candidates, qnan);
    std::vector<Double_t> GFpi_mass_container(xi_candidates, qnan);
    std::vector<Double_t> GFpi2_mass_container(xi_candidates, qnan);
    std::vector<Double_t> GFppi_closedist_container(xi_candidates, qnan);
    std::vector<Double_t> GFlpi_closedist_container(xi_candidates, qnan);

    Int_t best = -1; Double_t prev_massdiff = 9999.;
    for(Int_t candi=0;candi<xi_candidates;candi++){

      Int_t trackid_p = xi_p_container[candi];
      Int_t trackid_pi = xi_pi_container[candi];
      Int_t trackid_pi2 = xi_pi2_container[candi];
      Int_t repid_p = p_repid_container[candi];
      Int_t repid_pi = pi_repid_container[candi];
      Int_t repid_pi2 = pi2_repid_container[candi];

      Double_t GFextrapdist_decays[3]; Double_t GFmass_decays[3];
      TVector3 GFmom_decays[3]; TVector3 GFlambda_vert; double GFppi_dist = qnan;
      if(!GFTrackCont.FindVertex(trackid_p, trackid_pi,
				 repid_p, repid_pi,
				 GFextrapdist_decays[0], GFextrapdist_decays[1],
				 GFmom_decays[0], GFmom_decays[1],
				 GFppi_dist, GFlambda_vert,
				 vtx_scan_range)
	 || GFppi_dist > GFppi_distcut) continue;

      TLorentzVector GFLp(GFmom_decays[0], TMath::Sqrt(GFmom_decays[0].Mag()*GFmom_decays[0].Mag() + ProtonMass*ProtonMass));
      TLorentzVector GFLpi(GFmom_decays[1], TMath::Sqrt(GFmom_decays[1].Mag()*GFmom_decays[1].Mag() + PionMass*PionMass));
      TLorentzVector GFLlambda = GFLp + GFLpi;
      TVector3 GFlambda_mom = GFmom_decays[0] + GFmom_decays[1];
      TLorentzVector GFLlambda_fixed(GFlambda_mom, TMath::Sqrt(GFlambda_mom.Mag()*GFlambda_mom.Mag() + LambdaMass*LambdaMass));

      TVector3 GFxi_vert; Double_t GFlpi_dist = qnan; Double_t GFlambda_Ltracklen;
      if(!GFTrackCont.FindVertexXi(trackid_pi2, repid_pi2, GFlambda_vert, GFlambda_mom, GFlambda_Ltracklen, GFextrapdist_decays[2], GFmom_decays[2], GFlpi_dist, GFxi_vert, vtx_scan_range) || GFlpi_dist > GFlpi_distcut) continue;

      TLorentzVector GFLpi2(GFmom_decays[2], TMath::Sqrt(GFmom_decays[2].Mag()*GFmom_decays[2].Mag() + PionMass*PionMass));
      TLorentzVector GFLxi = GFLlambda_fixed + GFLpi2;
      TVector3 GFxi_mom = GFlambda_mom + GFmom_decays[2];

      ///////////////////
      std::cout<<"Lambda extrapolation"<<std::endl;
      //TVector3 GFlambdaVect = GFlambda_vert - GFxi_vert;
      Double_t GFlambda_tracklen = TMath::Abs(GFlambda_Ltracklen);
      Double_t GFlambda_tof = Kinematics::CalcTimeOfFlight(GFlambda_mom.Mag(), GFlambda_tracklen, pdg::LambdaMass());

      Int_t hitid_htof; Double_t tof_htof; Double_t tracklen_htof; TVector3 pos_htof;
      Bool_t htofextrapolation_p = GFTrackCont.TPCHTOFTrackMatching(trackid_p, repid_p,
								    event.HtofSeg, event.posHtof,
								    hitid_htof, tof_htof,
								    tracklen_htof, pos_htof);
      if(htofextrapolation_p){
	//Double_t GFp_tof = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
	Double_t GFp_tof = event.tHtof[hitid_htof] - GFlambda_tof;
	Double_t GFp_tracklen = tracklen_htof - GFextrapdist_decays[0];
	GFmass_decays[0] = Kinematics::MassSquare(GFmom_decays[0].Mag(), GFp_tracklen, GFp_tof);
      }
      else GFmass_decays[0] = qnan;

      Bool_t htofextrapolation_pi = GFTrackCont.TPCHTOFTrackMatching(trackid_pi, repid_pi,
								    event.HtofSeg, event.posHtof,
								    hitid_htof, tof_htof,
								    tracklen_htof, pos_htof);
      if(htofextrapolation_pi){
	//Double_t GFpi_tof = event.tHtof[hitid_htof] - GFlambda_tof - GFxi_tof;
	Double_t GFpi_tof = event.tHtof[hitid_htof] - GFlambda_tof;
	Double_t GFpi_tracklen = tracklen_htof - GFextrapdist_decays[1];
	GFmass_decays[1] = Kinematics::MassSquare(GFmom_decays[1].Mag(), GFpi_tracklen, GFpi_tof);
      }
      else GFmass_decays[1] = qnan;

      Bool_t htofextrapolation_pi2 = GFTrackCont.TPCHTOFTrackMatching(trackid_pi2, repid_pi2,
								      event.HtofSeg, event.posHtof,
								      hitid_htof, tof_htof,
								      tracklen_htof, pos_htof);
      if(htofextrapolation_pi2){
	//Double_t GFpi2_tof = event.tHtof[hitid_htof] - GFxi_tof;
	Double_t GFpi2_tof = event.tHtof[hitid_htof];
	Double_t GFpi2_tracklen = tracklen_htof - GFextrapdist_decays[2];
	GFmass_decays[2] = Kinematics::MassSquare(GFmom_decays[2].Mag(), GFpi2_tracklen, GFpi2_tof);
      }
      else GFmass_decays[2] = qnan;

      event.xiflag = true;
      GFxi_p_id_container[candi] = trackid_p;
      GFxi_pi_id_container[candi] = trackid_pi;
      GFxi_pi2_id_container[candi] = trackid_pi2;
      GFxi_p_rep_container[candi] = repid_p;
      GFxi_pi_rep_container[candi] = repid_pi;
      GFxi_pi2_rep_container[candi] = repid_pi2;
      GFxi_mom_container[candi] = GFxi_mom;
      GFxi_vert_container[candi] = GFxi_vert;
      GFxi_mass_container[candi] = GFLxi.M();

      GFl_mom_container[candi] = GFlambda_mom;
      GFl_vert_container[candi] = GFlambda_vert;
      GFl_mass_container[candi] = GFLlambda.M();
      GFl_tracklen_container[candi] = GFlambda_tracklen;
      GFl_tof_container[candi] = GFlambda_tof;

      GFp_mom_container[candi] = GFmom_decays[0];
      GFpi_mom_container[candi] = GFmom_decays[1];
      GFpi2_mom_container[candi] = GFmom_decays[2];
      GFp_mass_container[candi] = GFmass_decays[0];
      GFpi_mass_container[candi] = GFmass_decays[1];
      GFpi2_mass_container[candi] = GFmass_decays[2];

      GFppi_closedist_container[candi] = GFppi_dist;
      GFlpi_closedist_container[candi] = GFlpi_dist;

      Double_t diff = TMath::Abs(GFLlambda.M() - LambdaMass);
      if(prev_massdiff > diff){
	prev_massdiff = diff;
	best = candi;
      }
    } //candi

    std::cout<<"xi best combination "<<best<<std::endl;
    if(event.xiflag){
      event.ximass = xi_mass_container[best];
      event.xidecayvtx_x = xi_vert_container[best].x();
      event.xidecayvtx_y = xi_vert_container[best].y();
      event.xidecayvtx_z = xi_vert_container[best].z();
      event.ximom_x = xi_mom_container[best].x();
      event.ximom_y = xi_mom_container[best].y();
      event.ximom_z = xi_mom_container[best].z();
      event.lpi_dist = lpi_closedist[best];

      event.lmass = lambda_mass_container[best];
      event.ldecayvtx_x = l_vert_container[best].x();
      event.ldecayvtx_y = l_vert_container[best].y();
      event.ldecayvtx_z = l_vert_container[best].z();
      event.lmom_x = l_mom_container[best].x();
      event.lmom_y = l_mom_container[best].y();
      event.lmom_z = l_mom_container[best].z();
      event.ppi_dist = ppi_closedist[best];

      event.decays_mom.push_back(p_mom_container[best].Mag());
      event.decays_mom.push_back(pi_mom_container[best].Mag());
      event.decays_mom.push_back(pi2_mom_container[best].Mag());
      event.decays_mom_x.push_back(p_mom_container[best].x());
      event.decays_mom_x.push_back(pi_mom_container[best].x());
      event.decays_mom_x.push_back(pi2_mom_container[best].x());
      event.decays_mom_y.push_back(p_mom_container[best].y());
      event.decays_mom_y.push_back(pi_mom_container[best].y());
      event.decays_mom_y.push_back(pi2_mom_container[best].y());
      event.decays_mom_z.push_back(p_mom_container[best].z());
      event.decays_mom_z.push_back(pi_mom_container[best].z());
      event.decays_mom_z.push_back(pi2_mom_container[best].z());

      event.decays_id.push_back(xi_p_container[best]);
      event.decays_id.push_back(xi_pi_container[best]);
      event.decays_id.push_back(xi_pi2_container[best]);

      event.GFntTpc = GFntTpc;
      event.GFximass = GFxi_mass_container[best];
      event.GFxidecayvtx_x = GFxi_vert_container[best].x();
      event.GFxidecayvtx_y = GFxi_vert_container[best].y();
      event.GFxidecayvtx_z = GFxi_vert_container[best].z();
      event.GFximom = GFxi_mom_container[best].Mag();
      event.GFximom_x = GFxi_mom_container[best].x();
      event.GFximom_y = GFxi_mom_container[best].y();
      event.GFximom_z = GFxi_mom_container[best].z();
      event.GFlpi_dist = GFlpi_closedist_container[best];
      event.GFlmass = GFl_mass_container[best];
      event.GFldecayvtx_x = GFl_vert_container[best].x();
      event.GFldecayvtx_y = GFl_vert_container[best].y();
      event.GFldecayvtx_z = GFl_vert_container[best].z();
      event.GFlmom = GFl_mom_container[best].Mag();
      event.GFlmom_x = GFl_mom_container[best].x();
      event.GFlmom_y = GFl_mom_container[best].y();
      event.GFlmom_z = GFl_mom_container[best].z();
      event.GFltracklen = GFl_tracklen_container[best];
      event.GFltof = GFl_tof_container[best];
      event.GFppi_dist = GFppi_closedist_container[best];

      event.GFnhtrack.resize(GFntTpc);
      event.GFchisqr.resize(GFntTpc);
      event.GFcharge.resize(GFntTpc);
      event.GFtof.resize(GFntTpc);
      event.GFtracklen.resize(GFntTpc);
      event.GFpval.resize(GFntTpc);
      event.GFpdgcode.resize(GFntTpc);

      event.GFntdecays = 3;
      event.GFdecays_mass.resize(event.GFntdecays);
      event.GFdecays_mom.resize(event.GFntdecays);
      event.GFdecays_mom_x.resize(event.GFntdecays);
      event.GFdecays_mom_y.resize(event.GFntdecays);
      event.GFdecays_mom_z.resize(event.GFntdecays);
      event.GFmomloss.resize(event.GFntdecays);
      event.GFeloss.resize(event.GFntdecays);

      for( Int_t j=0; j<event.GFntdecays; ++j ){
	Int_t igf = GFxi_p_id_container[best];
	if(j==1) igf = GFxi_pi_id_container[best];
	if(j==2) igf = GFxi_pi2_id_container[best];

	Int_t repid = GFxi_p_rep_container[best];
	if(j==1) repid = GFxi_pi_rep_container[best];
	if(j==2) repid = GFxi_pi2_rep_container[best];

	TVector3 GFmom_decays = GFp_mom_container[best];
	if(j==1) GFmom_decays = GFpi_mom_container[best];
	if(j==2) GFmom_decays = GFpi2_mom_container[best];

	Double_t GFmass_decays = GFp_mass_container[best];
	if(j==1) GFmass_decays = GFpi_mass_container[best];
	if(j==2) GFmass_decays = GFpi2_mass_container[best];

	Int_t nh = GFTrackCont.GetNHits(igf);
	event.GFnhtrack[j] = nh;
	event.GFchisqr[j] = GFTrackCont.GetChi2NDF(igf, repid);
	event.GFcharge[j] = GFTrackCont.GetCharge(igf, repid);
	event.GFtof[j] = GFTrackCont.GetTrackTOF(igf, 0, -1, repid);
	event.GFtracklen[j] = GFTrackCont.GetTrackLength(igf, 0, -1, repid);
	event.GFpval[j] = GFTrackCont.GetPvalue(igf, repid);
	event.GFpdgcode[j] = GFTrackCont.GetPDGcode(igf, repid);
	event.GFdecays_mass[j] = GFmass_decays;
	event.GFdecays_mom[j] = GFmom_decays.Mag();
	event.GFdecays_mom_x[j] = GFmom_decays.x();
	event.GFdecays_mom_y[j] = GFmom_decays.y();
	event.GFdecays_mom_z[j] = GFmom_decays.z();
	event.GFmomloss[j] = GFTrackCont.GetMom(igf, 0, repid).Mag() - GFmom_decays.Mag();
	event.GFeloss[j] = TMath::Sqrt(GFTrackCont.GetMom(igf, 0, repid).Mag()*GFTrackCont.GetMom(igf, 0, repid).Mag() + pdgmass[j]*pdgmass[j]) - TMath::Sqrt(GFmom_decays.Mag()*GFmom_decays.Mag() + pdgmass[j]*pdgmass[j]);
      } //j : decays

      std::cout<<"Xi extrapolation"<<std::endl;
      GFTrackCont.AddReconstructedTrack(XiMinusPdgCode, GFxi_vert_container[best],
					GFxi_mom_container[best]);
      GFTrackCont.FitTrack(GFTrackCont.GetNTrack() - 1);

      TVector3 kkvtx(event.vtx[0], event.vty[0], event.vtz[0] + tpc::ZTarget);
      TVector3 xiprodvtx; TVector3 ximom_prodvtx; double xitracklen_kk; double xitof_kk;
      Bool_t xi_flight = GFTrackCont.ExtrapolateToPoint(GFTrackCont.GetNTrack() - 1, kkvtx,
							xiprodvtx, ximom_prodvtx,
							xitracklen_kk, xitof_kk);
      event.GFxiprodvtx_x = xiprodvtx.x();
      event.GFxiprodvtx_y = xiprodvtx.y();
      event.GFxiprodvtx_z = xiprodvtx.z();
      event.GFxiprodmom = ximom_prodvtx.Mag();
      event.GFxiprodmom_x = ximom_prodvtx.x();
      event.GFxiprodmom_y = ximom_prodvtx.y();
      event.GFxiprodmom_z = ximom_prodvtx.z();
      event.GFxitracklen = xitracklen_kk;
      event.GFxitof = xitof_kk;
      ///////////////////////////////////

      HF1( 10, event.MissMass[0] );
      HF1( 11, event.lmass );
      HF1( 12, event.ximass );
      HF1( 13, event.decays_mom[0] );
      HF1( 14, event.decays_mom[1] );
      HF1( 15, event.decays_mom[2] );
      HF1( 16, TMath::Sqrt( event.lmom_x*event.lmom_x +
			    event.lmom_y*event.lmom_y +
			    event.lmom_z*event.lmom_z) );
      HF1( 17, TMath::Sqrt( event.ximom_x*event.ximom_x +
			    event.ximom_y*event.ximom_y +
			    event.ximom_z*event.ximom_z) );
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
      HF1( 30, event.GFdecays_mass[0]);
      HF1( 31, event.GFdecays_mass[1]);
      HF1( 32, event.GFdecays_mass[2]);

      TLorentzVector LvXi_fixedmass(TVector3(event.GFxiprodmom_x, event.GFxiprodmom_y, event.GFxiprodmom_z),
				    TMath::Hypot(TVector3(event.GFxiprodmom_x, event.GFxiprodmom_y, event.GFxiprodmom_z).Mag(), XiMinusMass));
      TLorentzVector LvRcXi_fixedmass = LvRc - LvXi_fixedmass;

      HF1( 201, LvRcXi_fixedmass.M() - m11B - me);

      event.p_multi = target_p_container.size();
      event.pi_multi = target_pip_container.size() + target_pim_container.size();
      for(Int_t itp=0; itp<target_p_container.size(); ++itp){
	if(itp==event.decays_id[0]){
	  event.p_multi -= 1;
	  continue;
	}
	event.ppi_mass.push_back(target_p_mass_container[itp]);
	event.ppi_mom.push_back(target_p_mom_container[itp].Mag());
	event.ppi_mom_x.push_back(target_p_mom_container[itp].x());
	event.ppi_mom_y.push_back(target_p_mom_container[itp].y());
	event.ppi_mom_z.push_back(target_p_mom_container[itp].z());
	event.ppi_charge.push_back(1);
	HF1( 1120, target_p_mom_container[itp].Mag());
      }
      for(Int_t itpip=0; itpip<target_pip_container.size(); ++itpip){
	event.ppi_mass.push_back(target_pip_mass_container[itpip]);
	event.ppi_mom.push_back(target_pip_mom_container[itpip].Mag());
	event.ppi_mom_x.push_back(target_pip_mom_container[itpip].x());
	event.ppi_mom_y.push_back(target_pip_mom_container[itpip].y());
	event.ppi_mom_z.push_back(target_pip_mom_container[itpip].z());
	event.ppi_charge.push_back(1);
	HF1( 1121, target_pip_mom_container[itpip].Mag());
      }
      for(Int_t itpim=0; itpim<target_pim_container.size(); ++itpim){
	if(itpim==event.decays_id[1] || itpim==event.decays_id[2]){
	  event.pi_multi -= 1;
	  continue;
	}
	event.ppi_mass.push_back(target_pim_mass_container[itpim]);
	event.ppi_mom.push_back(target_pim_mom_container[itpim].Mag());
	event.ppi_mom_x.push_back(target_pim_mom_container[itpim].x());
	event.ppi_mom_y.push_back(target_pim_mom_container[itpim].y());
	event.ppi_mom_z.push_back(target_pim_mom_container[itpim].z());
	event.ppi_charge.push_back(-1);
	HF1( 1122, target_pim_mom_container[itpim].Mag());
      }
      event.ppi_multi = event.p_multi + event.pi_multi;

      //find Xi-p event
      if(event.pi_multi==0 && event.p_multi==1){
	HF1( 500, -event.BE[0]);
      }

      HF1( 1020, event.ppi_multi);
      HF1( 1021, event.p_multi);
      HF1( 1022, event.pi_multi);

      HF1( 120, -event.BE[0]);
      if(event.isgoodTPC[0] == 1) HF1( 122,- event.BETPC[0]);
      return true;
    } //if(event.xiflag)
  } //xi searching

  std::cout<<"L Searching starts"<<std::endl;
  if(L_p_container.size()>0){
    int best = -1; double prev_mass = 9999.;
    for(int id=0;id<L_p_container.size();id++){
      if(TMath::Abs(L_mass_container[id] - LambdaMass)<prev_mass){
	best = id; prev_mass = L_mass_container[id];
      }
    }
    if(Lpi_p_container.size()>0){
      event.lpiflag = true;
      HF1( 134, -event.BE[0]);
      if(event.isgoodTPC[0] == 1) HF1( 136, -event.BETPC[0]);
    }
    else event.lflag = true;

    event.lmass = L_mass_container[best];
    event.ldecayvtx_x = L_vtx_container[best].x();
    event.ldecayvtx_y = L_vtx_container[best].y();
    event.ldecayvtx_z = L_vtx_container[best].z();
    event.lmom_x = L_mom_container[best].x();
    event.lmom_y = L_mom_container[best].y();
    event.lmom_z = L_mom_container[best].z();
    event.ppi_dist = L_ppidist_container[best];

    event.p_multi = target_p_container.size();
    event.pi_multi = target_pip_container.size() + target_pim_container.size();
    for(Int_t itp=0; itp<target_p_container.size(); ++itp){
      if(itp==L_p_container[best]){
	event.p_multi -= 1;
	continue;
      }
      event.ppi_mass.push_back(target_p_mass_container[itp]);
      event.ppi_mom.push_back(target_p_mom_container[itp].Mag());
      event.ppi_mom_x.push_back(target_p_mom_container[itp].x());
      event.ppi_mom_y.push_back(target_p_mom_container[itp].y());
      event.ppi_mom_z.push_back(target_p_mom_container[itp].z());
      event.ppi_charge.push_back(1);
      HF1( 1130, target_p_mom_container[itp].Mag());
    }
    for(Int_t itpip=0; itpip<target_pip_container.size(); ++itpip){
      event.ppi_mass.push_back(target_pip_mass_container[itpip]);
      event.ppi_mom.push_back(target_pip_mom_container[itpip].Mag());
      event.ppi_mom_x.push_back(target_pip_mom_container[itpip].x());
      event.ppi_mom_y.push_back(target_pip_mom_container[itpip].y());
      event.ppi_mom_z.push_back(target_pip_mom_container[itpip].z());
      event.ppi_charge.push_back(1);
      HF1( 1131, target_pip_mom_container[itpip].Mag());
    }
    for(Int_t itpim=0; itpim<target_pim_container.size(); ++itpim){
      if(itpim==L_pi_container[best]){
	event.pi_multi -= 1;
	continue;
      }
      event.ppi_mass.push_back(target_pim_mass_container[itpim]);
      event.ppi_mom_x.push_back(target_pim_mom_container[itpim].x());
      event.ppi_mom_y.push_back(target_pim_mom_container[itpim].y());
      event.ppi_mom_z.push_back(target_pim_mom_container[itpim].z());
      event.ppi_charge.push_back(-1);
      HF1( 1132, target_pim_mom_container[itpim].Mag());
    }
    event.ppi_multi = event.p_multi + event.pi_multi;

    HF1( 1030, event.ppi_multi);
    HF1( 1031, event.p_multi);
    HF1( 1032, event.pi_multi);

    HF1( 130, -event.BE[0]);
    if(event.isgoodTPC[0] == 1) HF1( 132, -event.BETPC[0]);

    return true;
  }

  std::cout<<"No Xi-,L"<<std::endl;
  std::vector<Double_t> twopi;
  event.p_multi = target_p_container.size();
  event.pi_multi = target_pip_container.size() + target_pim_container.size();
  for(Int_t itp=0; itp<target_p_container.size(); ++itp){
    event.ppi_mass.push_back(target_p_mass_container[itp]);
    event.ppi_mom.push_back(target_p_mom_container[itp].Mag());
    event.ppi_mom_x.push_back(target_p_mom_container[itp].x());
    event.ppi_mom_y.push_back(target_p_mom_container[itp].y());
    event.ppi_mom_z.push_back(target_p_mom_container[itp].z());
    event.ppi_charge.push_back(1);
    HF1( 1100, target_p_mom_container[itp].Mag());
  }
  for(Int_t itpip=0; itpip<target_pip_container.size(); ++itpip){
    event.ppi_mass.push_back(target_pip_mass_container[itpip]);
    event.ppi_mom.push_back(target_pip_mom_container[itpip].Mag());
    event.ppi_mom_x.push_back(target_pip_mom_container[itpip].x());
    event.ppi_mom_y.push_back(target_pip_mom_container[itpip].y());
    event.ppi_mom_z.push_back(target_pip_mom_container[itpip].z());
    event.ppi_charge.push_back(1);
    HF1( 1101, target_pip_mom_container[itpip].Mag());
  }
  for(Int_t itpim=0; itpim<target_pim_container.size(); ++itpim){
    event.ppi_mass.push_back(target_pim_mass_container[itpim]);
    event.ppi_mom.push_back(target_pim_mom_container[itpim].Mag());
    event.ppi_mom_x.push_back(target_pim_mom_container[itpim].x());
    event.ppi_mom_y.push_back(target_pim_mom_container[itpim].y());
    event.ppi_mom_z.push_back(target_pim_mom_container[itpim].z());
    event.ppi_charge.push_back(-1);
    HF1( 1102, target_pim_mom_container[itpim].Mag());

    double KE = 1000.*(TMath::Hypot(target_pim_mom_container[itpim].Mag(), PionMass) - PionMass); //MeV/c2
    if(KE < 60.) twopi.push_back(KE);
  }
  event.ppi_multi = event.p_multi + event.pi_multi;

  if(twopi.size()>=2){
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

  HF1( 1000, event.ppi_multi);
  HF1( 1001, event.p_multi);
  HF1( 1002, event.pi_multi);

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

  HB1( 21, "(GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 22, "(GenFit] #Xi Invariant Mass; M_{#Lambda#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 125, 1.2 ,1.45);
  HB1( 23, "(GenFit] p_{#Lambda} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 24, "(GenFit] #pi_{#Lambda} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 25, "(GenFit] #pi_{#Xi} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 26, "(GenFit] #Lambda Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 27, "(GenFit] #Xi Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 28, "(GenFit] Closest Dist. p#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 29, "(GenFit] Closest Dist. #Lambda#pi;Distance [mm]; Counts (/0.1 mm)", 1000, 0., 100.);
  HB1( 30, "(GenFit] p_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 31, "(GenFit] pi_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);
  HB1( 32, "(GenFit] pi_{#Xi} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0., 0.5);

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

  HB1( 201, ";^{12}C(K^{-}, K^{+}#Xi^{-})X MM - M(^{11}B) (GeV/#font[12]{c}^{2});Counts (/5 MeV/#font[12]{c}^{2})", 80, -0.1, 0.3);

  HB1( 500, "; -B_{#Xi^{-}} (MeV/#font[12]{c}^{2});Counts (/5 MeV/#font[12]{c}^{2})", 280, -400., 1000.);
  HB2( 501, ";T_{#piL} (MeV/#font[12]{c}); T_{#piH} (MeV/#font[12]{c})", 60, 0, 60, 60, 0, 60);
  HB2( 502, ";M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); ;M_{p#pi^{-}} (GeV/#font[12]{c}^{2})", 80, 1.04, 1.2, 80, 1.04, 1.2);
  HB2( 503, "(Genfit] L Vs L Invariant Mass;M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); ;M_{p#pi^{-}} (GeV/#font[12]{c}^{2})", 80, 1.04, 1.2, 80, 1.04, 1.2);
  HB1( 504, "#Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);
  HB1( 505, "(GenFit] #Lambda Invariant Mass; M_{p#pi^{-}} (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 80, 1.04, 1.2);

  HB1( 506, "p_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 507, "pi_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);
  HB1( 508, "(GenFit] p_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 350, 0.5, 1.2);
  HB1( 509, "(GenFit] pi_{#Lambda} Mass; M (GeV/#font[12]{c}^{2}); Counts (/2 MeV/#font[12]{c}^{2})", 250, 0.0, 0.5);

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
  HB1( 1120, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 1121, "#pi^{+} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 1122, "#pi^{-} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);

  HB1( 1030, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1031, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1032, "; Multiplicity;Counts ", 5,0,5);
  HB1( 1130, "p Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 100, 0., 2.0);
  HB1( 1131, "#pi^{+} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  HB1( 1132, "#pi^{-} Momentum; Momentum (GeV/#font[12]{c}); Counts (/0.02 GeV/#font[12]{c})", 50, 0., 1.0);
  */
  HBTree( "tpc", "tree of DstTPCTracking" );
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

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "trackid", &event.trackid );
  tree->Branch( "isBeam", &event.isBeam );
  tree->Branch( "isKurama", &event.isKurama );
  tree->Branch( "isK18", &event.isK18 );
  tree->Branch( "isAccidental", &event.isAccidental );
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
  tree->Branch( "helix_t", &event.helix_t );

#if SaveRawData
  tree->Branch( "ntK18",      &event.ntK18);
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

  tree->Branch(" ntKurama", &event.ntKurama);
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
#endif
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

  tree->Branch("nKK", &event.nKK);
  tree->Branch("Kflag", &event.Kflag);
  tree->Branch("MissMass", &event.MissMass);
  tree->Branch("MissMassCorr", &event.MissMassCorr);
  tree->Branch("MissMassCorrDE", &event.MissMassCorrDE);
  tree->Branch("vtx", &event.vtx);
  tree->Branch("vty", &event.vty);
  tree->Branch("vtz", &event.vtz);
  tree->Branch("closeDist", &event.closeDist);

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

  //track fitting results
  tree->Branch("GFntTpc", &event.GFntTpc);
  tree->Branch("GFnhtrack", &event.GFnhtrack);
  tree->Branch("GFchisqr", &event.GFchisqr);
  tree->Branch("GFcharge", &event.GFcharge);
  tree->Branch("GFtof", &event.GFtof);
  tree->Branch("GFtracklen", &event.GFtracklen);
  tree->Branch("GFpval", &event.GFpval);
  tree->Branch("GFpdgcode", &event.GFpdgcode);
  tree->Branch("GFntDecays", &event.GFntdecays);
  tree->Branch("GFDecaysMass", &event.GFdecays_mass);
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
  tree->Branch("XiMom_x", &event.ximom_x);
  tree->Branch("XiMom_y", &event.ximom_y);
  tree->Branch("XiMom_z", &event.ximom_z);
  tree->Branch("XiVtxCloseDist", &event.lpi_dist);
  tree->Branch("LambdaMass", &event.lmass);
  tree->Branch("LambdaDecayVtx_x", &event.ldecayvtx_x);
  tree->Branch("LambdaDecayVtx_y", &event.ldecayvtx_y);
  tree->Branch("LambdaDecayVtx_z", &event.ldecayvtx_z);
  tree->Branch("LambdaMom_x", &event.lmom_x);
  tree->Branch("LambdaMom_y", &event.lmom_y);
  tree->Branch("LambdaMom_z", &event.lmom_z);
  tree->Branch("LambdaVtxCloseDist", &event.ppi_dist);

  tree->Branch("GFXiMass", &event.GFximass);
  tree->Branch("GFXiDecayVtx_x", &event.GFxidecayvtx_x);
  tree->Branch("GFXiDecayVtx_y", &event.GFxidecayvtx_y);
  tree->Branch("GFXiDecayVtx_z", &event.GFxidecayvtx_z);
  tree->Branch("GFXiMom", &event.GFximom);
  tree->Branch("GFXiMom_x", &event.GFximom_x);
  tree->Branch("GFXiMom_y", &event.GFximom_y);
  tree->Branch("GFXiMom_z", &event.GFximom_z);
  tree->Branch("GFXiProductionVtx_x", &event.GFxiprodvtx_x);
  tree->Branch("GFXiProductionVtx_y", &event.GFxiprodvtx_y);
  tree->Branch("GFXiProductionVtx_z", &event.GFxiprodvtx_z);
  tree->Branch("GFXiProductionVtxMom", &event.GFxiprodmom);
  tree->Branch("GFXiProductionVtxMom_x", &event.GFxiprodmom_x);
  tree->Branch("GFXiProductionVtxMom_y", &event.GFxiprodmom_y);
  tree->Branch("GFXiProductionVtxMom_z", &event.GFxiprodmom_z);
  tree->Branch("GFXiTrackLen", &event.GFxitracklen);
  tree->Branch("GFXiTof", &event.GFxitof);
  tree->Branch("GFXiVtxCloseDist", &event.GFlpi_dist);

  tree->Branch("GFLambdaMass", &event.GFlmass);
  tree->Branch("GFLambdaDecayVtx_x", &event.GFldecayvtx_x);
  tree->Branch("GFLambdaDecayVtx_y", &event.GFldecayvtx_y);
  tree->Branch("GFLambdaDecayVtx_z", &event.GFldecayvtx_z);
  tree->Branch("GFLambdaMom_x", &event.GFlmom_x);
  tree->Branch("GFLambdaMom_y", &event.GFlmom_y);
  tree->Branch("GFLambdaMom_z", &event.GFlmom_z);
  tree->Branch("GFLambdaTrackLen", &event.GFltracklen);
  tree->Branch("GFLambdaTof", &event.GFltof);
  tree->Branch("GFLambdaVtxCloseDist", &event.GFppi_dist);

  tree->Branch("Lflag", &event.lflag);
  tree->Branch("LPiflag", &event.lpiflag);
  tree->Branch("LLflag", &event.llflag);
  tree->Branch("LambdaMass1", &event.lmass1);
  tree->Branch("LambdaDecayVtx_x1", &event.ldecayvtx_x1);
  tree->Branch("LambdaDecayVtx_y1", &event.ldecayvtx_y1);
  tree->Branch("LambdaDecayVtx_z1", &event.ldecayvtx_z1);
  tree->Branch("LambdaMom_x1", &event.lmom_x1);
  tree->Branch("LambdaMom_y1", &event.lmom_y1);
  tree->Branch("LambdaMom_z1", &event.lmom_z1);
  tree->Branch("LambdaVtxCloseDist1", &event.ppi_dist1);
  tree->Branch("LambdaMass2", &event.lmass2);
  tree->Branch("LambdaDecayVtx_x2", &event.ldecayvtx_x2);
  tree->Branch("LambdaDecayVtx_y2", &event.ldecayvtx_y2);
  tree->Branch("LambdaDecayVtx_z2", &event.ldecayvtx_z2);
  tree->Branch("LambdaMom_x2", &event.lmom_x2);
  tree->Branch("LambdaMom_y2", &event.lmom_y2);
  tree->Branch("LambdaMom_z2", &event.lmom_z2);
  tree->Branch("LambdaVtxCloseDist2", &event.ppi_dist2);

  tree->Branch("GFLambdaMass1", &event.GFlmass1);
  tree->Branch("GFLambdaDecayVtx_x1", &event.GFldecayvtx_x1);
  tree->Branch("GFLambdaDecayVtx_y1", &event.GFldecayvtx_y1);
  tree->Branch("GFLambdaDecayVtx_z1", &event.GFldecayvtx_z1);
  tree->Branch("GFLambdaMom_x1", &event.GFlmom_x1);
  tree->Branch("GFLambdaMom_y1", &event.GFlmom_y1);
  tree->Branch("GFLambdaMom_z1", &event.GFlmom_z1);
  tree->Branch("GFLambdaVtxCloseDist1", &event.GFppi_dist1);
  tree->Branch("GFLambdaMass2", &event.GFlmass2);
  tree->Branch("GFLambdaDecayVtx_x2", &event.GFldecayvtx_x2);
  tree->Branch("GFLambdaDecayVtx_y2", &event.GFldecayvtx_y2);
  tree->Branch("GFLambdaDecayVtx_z2", &event.GFldecayvtx_z2);
  tree->Branch("GFLambdaMom_x2", &event.GFlmom_x2);
  tree->Branch("GFLambdaMom_y2", &event.GFlmom_y2);
  tree->Branch("GFLambdaMom_z2", &event.GFlmom_z2);
  tree->Branch("GFLambdaVtxCloseDist2", &event.GFppi_dist2);

  tree->Branch("PiPiflag", &event.pipiflag);
  tree->Branch("ppiMultiplicity", &event.ppi_multi);
  tree->Branch("pMultiplicity", &event.p_multi);
  tree->Branch("piMultiplicity", &event.pi_multi);
  tree->Branch("ppiMass", &event.ppi_mass);
  tree->Branch("ppiMom", &event.ppi_mom);
  tree->Branch("ppiMom_x", &event.ppi_mom_x);
  tree->Branch("ppiMom_y", &event.ppi_mom_y);
  tree->Branch("ppiMom_z", &event.ppi_mom_z);
  tree->Branch("ppiCharge", &event.ppi_charge);

  tree->Branch("DecaysTrackId", &event.decays_id);
  tree->Branch("DecaysMom", &event.decays_mom);
  tree->Branch("DecaysMom_x", &event.decays_mom_x);
  tree->Branch("DecaysMom_y", &event.decays_mom_y);
  tree->Branch("DecaysMom_z", &event.decays_mom_z);

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
  src.xtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "xtgtKurama" );
  src.ytgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ytgtKurama" );
  src.utgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "utgtKurama" );
  src.vtgtKurama = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vtgtKurama" );

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

  src.ntTpc = new TTreeReaderValue<Int_t>( *reader, "ntTpc" );
  src.nhtrack = new TTreeReaderValue<std::vector<Int_t>>( *reader, "nhtrack" );
  src.trackid = new TTreeReaderValue<std::vector<Int_t>>( *reader, "trackid" );
  src.isBeam = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isBeam" );
  src.isKurama = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isKurama" );
  src.isK18 = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isK18" );
  src.isAccidental = new TTreeReaderValue<std::vector<Int_t>>( *reader, "isAccidental" );
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

#if SaveTPCK18
  src.ntTPCK18 = new TTreeReaderValue<Int_t>( *reader, "ntK18" );
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

  src.ntTPCKurama = new TTreeReaderValue<Int_t>( *reader, "ntKurama" );
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
  src.ubTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "ubTPC" );
  src.vbTPC = new TTreeReaderValue<std::vector<Double_t>>( *reader, "vbTPC" );
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
