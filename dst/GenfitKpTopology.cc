// -*- C++ -*-
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <vector>
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
#include "TH1.h"
#include "TH2.h"
#include "PidCommon.hh"
#include "KpTopologyAnalysis.hh"

#define MakePidFig 1
// Set to 1 to also write the full per-hit DST tree of GenfitQFKaon.cc. Off by
// default: the topology tree is meant to be reprocessed many times with
// different cuts, which is only practical if it stays small.
#define WriteFullDst 0

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
  
const double minThetaKPCH2 = 1.5;
const double maxThetaKPCH2 = 10.5;

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
  //static TH2* gAcceptance = nullptr;
  static TH1* gAcceptance = nullptr;

//_____________________________________________________________________________
// Topology analysis definitions
//
// The sign convention used throughout this file is B_K > 0 on the deeply bound
// side, i.e. BKaon = M(11B) + M(K-) - MM(12C). GenfitQFKaon.cc plots -BETPC and
// calls it "bek"; that variable is deliberately not reused here.

namespace topodef
{
// Role of a TPC track with respect to the primary vertex.
// Existing numeric values must not change; new roles are appended before kNTrackRole.
enum ETrackRole {
  kRoleUnknown = 0,
  kRoleBeam,        // beam-like track, not the tagged K18 track
  kRoleK18,         // incident K- matched to K18
  kRoleForward,     // forward proton matched to KURAMA
  kRoleAccidental,
  kRoleV0Daughter,  // Lambda SIGNAL daughter (nominal exclusion)
  kRolePrompt,      // compatible with the primary vertex
  kRoleDisplaced,   // inside target but not pointing to the primary vertex
  kRoleOutside,     // does not cross the target
  kRoleK0Daughter,  // K0S candidate daughter (appended; do not reorder above)
  kNTrackRole
};

// Particle classes. The dE/dx pid word is a set of allowed hypotheses rather
// than an exclusive identification, so the ambiguous combinations are kept as
// their own classes instead of being forced into one species.
enum EPidClass {
  kPidUnknown = 0,
  kPidProton,       // proton, not compatible with pion
  kPidPiPlus,
  kPidPiMinus,
  kPidProtonOrPiPlus, // positive, compatible with both (= Ambiguous)
  kPidKMinus,
  kPidElectron,
  // Appended (do not reorder above). Not packed into exclusiveVisibleKey v1.
  kPidDeuteron,
  kPidTriton,
  kPidHeavyAmbiguous, // d/t overlap or heavy-fragment ambiguous
  kNPidClass
};

// Coarse exclusive class kept for QA / backward compatibility.
enum ETopoClass {
  kTopoEmpty = 0,
  kTopo1Proton,
  kTopo1PiPlus,
  kTopo1PiMinus,
  kTopo1Other,
  kTopo2Prompt,
  kTopo3PlusPrompt,
  kNTopoClass
};

// Semi-exclusive tags: non-mutually-exclusive bitset (else-if forbidden).
// Naming: LambdaPiMinus = 1NA-enriched; LambdaProton = 2NA-enriched;
// LambdaPiPlus = control tag (not direct 1NA); TwoPiPlus != Sigma+.
enum ESemiTagBit {
  kTagLambda = 0,                 // lIsSignal
  kTagLambdaPiMinus = 1,          // 1NA-enriched
  kTagLambdaPiPlus = 2,           // control tag
  kTagLambdaProton = 3,           // 2NA-enriched
  kTagLambdaProtonNoPion = 4,
  kTagTwoPiPlus = 5,              // residual pi+ >= 2 (not named Sigma+)
  kTagSigmaPlusCharged = 6,       // prompt pi- + displaced pi+  (Sigma+ enriched)
  kTagKminusEscape = 7,           // escFlag
  kTagK0S = 8,
  kTagK0SProton = 9,
  kTagLambdaSideband = 10,
  kTagLambdaSidebandPiMinus = 11,
  kTagLambdaSidebandPiPlus = 12,
  kTagLambdaSidebandProton = 13,
  kTagLambdaSidebandProtonNoPion = 14,
  kTagSigmaMinusCharged = 15,     // prompt pi+ + displaced pi-  (Sigma- enriched)
  // Appended tags (semiTagVersion >= 3). Existing 0–15 bits unchanged.
  kTagLambdaDeuteron = 16,
  kTagLambdaTriton = 17,
  kTagV0Ambiguous = 18,           // overlapping Lambda/K0S assignment
  kNSemiTagBit
};

inline ULong64_t SemiTagMask(ESemiTagBit b){ return (ULong64_t)1 << (Int_t)b; }

// Working points for the multiplicity counting. Index 0 is nominal; the others
// move one axis at a time so that the threshold scan required before any
// physics interpretation can be done on the output tree alone.
const Int_t kNThr = 9;
const Double_t kThrPMin[kNThr]  = { 0.05, 0.00, 0.08, 0.10, 0.15, 0.05, 0.05, 0.05, 0.05 }; // GeV/c
const Int_t    kThrNClust[kNThr]= {    8,    8,    8,    8,    8,    6,   10,    8,    8 };
const Double_t kThrDca[kNThr]   = { 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 10.0, 30.0 }; // mm

// Separation between prompt and displaced tracks is scanned through kThrDca;
// this value only defines the nominal role label written per track.
const Double_t kNominalDca = 20.0; // mm
// Narrow mass window for nominal Lambda physics tags (not the wide input lFlag).
const Double_t kLambdaMassWindow = 0.005; // GeV/c2, |dM| < 5 MeV signal
const Double_t kLambdaSidebandLow = 0.015; // GeV/c2, |dM| in [low, high)
const Double_t kLambdaSidebandHigh = 0.035;
const Double_t kLambdaFwdMassWindow = 0.005; // GeV/c2, M(p_fwd, pi-)
const Double_t kK0SMassWindow = 0.020; // GeV/c2, |M(pi+pi-)-m_K0| signal
const Double_t kK0SSidebandLow = 0.030;
const Double_t kK0SSidebandHigh = 0.050;
// Flat sideband scale factors (signal width / sideband width) for scaled-SB subtraction.
const Double_t kLambdaSidebandScale =
  kptopo::FlatSidebandScale(kLambdaMassWindow, kLambdaSidebandLow, kLambdaSidebandHigh);
const Double_t kK0SSidebandScale =
  kptopo::FlatSidebandScale(kK0SMassWindow, kK0SSidebandLow, kK0SSidebandHigh);

// Classification schema versions stored in the tree.
// roleClassVersion 2: nominal prompt DCA uses K18×forward (kk) vertex; dual vertices.
// exclusiveVisibleVersion 1: packing unchanged (d/t not in key).
// semiTagVersion 3: + LambdaDeuteron/Triton/V0Ambiguous; Lambda residuals exclude K0 daughters.
// v0CandidateVersion / kinematicsVersion: new analysis blocks.
const Int_t kRoleClassVersion = 2;
const Int_t kExclusiveVisibleVersion = 1;
const Int_t kSemiTagVersion = 3;
const Int_t kV0CandidateVersion = 1;
const Int_t kKinematicsVersion = 1;

// Coarse B_K edges for momentum spectroscopy (shared with hist booking).
const Int_t kNBKMom = 5;
const Double_t kBKMomEdge[kNBKMom+1] = {-0.30,-0.15,-0.05,0.00,0.05,0.30};

inline Int_t BKMomBin(Double_t bk)
{
  for(Int_t i=0; i<kNBKMom; ++i)
    if(bk >= kBKMomEdge[i] && bk < kBKMomEdge[i+1]) return i;
  return -1;
}

// Visible-exclusive key packing (ULong64_t):
//   16 fields x 3 bits (cap at 7) = 48 bits.
// Field order:
//   0 nLambdaSignal, 1 nK0S,
//   2..8  prompt: p, pi+, pi-, K-, e, ambiguous, unknown
//   9..15 displaced: same order
inline Int_t CapMult7(Int_t n)
{
  if(n < 0) return 0;
  if(n > 7) return 7;
  return n;
}

inline ULong64_t EncodeExclusiveVisible(
  Int_t nLam, Int_t nK0,
  Int_t pP, Int_t pPip, Int_t pPim, Int_t pKm, Int_t pE, Int_t pA, Int_t pU,
  Int_t dP, Int_t dPip, Int_t dPim, Int_t dKm, Int_t dE, Int_t dA, Int_t dU)
{
  ULong64_t key = 0;
  Int_t sh = 0;
  auto pack = [&](Int_t v){
    key |= (ULong64_t)CapMult7(v) << sh;
    sh += 3;
  };
  pack(nLam); pack(nK0);
  pack(pP); pack(pPip); pack(pPim); pack(pKm); pack(pE); pack(pA); pack(pU);
  pack(dP); pack(dPip); pack(dPim); pack(dKm); pack(dE); pack(dA); pack(dU);
  return key;
}

inline Int_t DecodeExclusiveField(ULong64_t key, Int_t fieldIndex)
{
  return (Int_t)((key >> (3*fieldIndex)) & 0x7ULL);
}

// Named exclusive spectra for QA (not the full ULong64_t key space).
enum EExclNamed {
  kExclNamedEmpty = 0,
  kExclNamed1Proton,
  kExclNamed1PiPlus,
  kExclNamed1PiMinus,
  kExclNamed1Other,
  kExclNamedLambdaOnly,     // Lambda signal, no residual charged
  kExclNamedLambdaPiMinus,  // Lambda + residual prompt pi-
  kExclNamedLambdaProton,   // Lambda + residual prompt p
  kExclNamedLambdaOther,
  kExclNamedK0Only,         // K0S, no residual charged
  kExclNamedK0Proton,       // K0S + residual prompt p
  kExclNamedK0Other,
  kExclNamedMultiOther,
  kNExclNamed
};

inline Int_t ClassifyExclusiveNamed(
  Int_t nLam, Int_t nK0,
  Int_t pP, Int_t pPip, Int_t pPim, Int_t pKm, Int_t pE, Int_t pA, Int_t pU,
  Int_t dP, Int_t dPip, Int_t dPim, Int_t dKm, Int_t dE, Int_t dA, Int_t dU)
{
  const Int_t nP = pP + pPip + pPim + pKm + pE + pA + pU;
  const Int_t nD = dP + dPip + dPim + dKm + dE + dA + dU;
  const Int_t nRes = nP + nD;

  if(nLam==0 && nK0==0){
    if(nRes==0) return kExclNamedEmpty;
    if(nD==0 && nP==1){
      if(pP==1)   return kExclNamed1Proton;
      if(pPip==1) return kExclNamed1PiPlus;
      if(pPim==1) return kExclNamed1PiMinus;
      return kExclNamed1Other;
    }
    return kExclNamedMultiOther;
  }
  if(nLam>=1 && nK0==0){
    if(nRes==0) return kExclNamedLambdaOnly;
    if(nD==0 && pPim>=1 && (nP - pPim)==0) return kExclNamedLambdaPiMinus;
    if(nD==0 && pP>=1 && (nP - pP)==0) return kExclNamedLambdaProton;
    return kExclNamedLambdaOther;
  }
  if(nK0>=1 && nLam==0){
    if(nRes==0) return kExclNamedK0Only;
    if(nD==0 && pP>=1 && (nP - pP)==0) return kExclNamedK0Proton;
    return kExclNamedK0Other;
  }
  return kExclNamedMultiOther;
}
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

//_____________________________________________________________________________
// Slim output of the topology analysis. Everything needed for the B_K-resolved
// topology, particle-momentum spectroscopy and Lambda-tag studies from the
// quasi-free region through the deeply bound side is written here; the heavy
// per-hit containers of Event stay internal to this program.
//
// "Exclusive" here means visible-exclusive (TPC-reconstructed charged tracks),
// not a fully exclusive reaction including neutrals / residual nucleus.
struct Topo
{
  Int_t runnum;
  Int_t evnum;
  Int_t evStatus;
  Int_t trigA;
  Int_t trigB;
  Int_t PSfacTrigA;
  Int_t PSfacTrigB;

  // classification schema versions
  Int_t roleClassVersion;
  Int_t exclusiveVisibleVersion;
  Int_t semiTagVersion;
  Int_t v0CandidateVersion;
  Int_t kinematicsVersion;

  // production coordinate
  Double_t BKaon;          // GeV, positive on the deeply bound side
  Double_t MissMassNucl;   // GeV
  Double_t thetaKP;        // deg
  Double_t qTransfer;      // GeV/c, |p_beam - p_forward|
  Double_t pBeam, pBeam_x, pBeam_y, pBeam_z;
  Double_t pFwd, pFwd_x, pFwd_y, pFwd_z;
  Double_t thetaFwd, phiFwd;
  Double_t PX_x, PX_y, PX_z, PX_E; // missing system four-momentum

  // vertices
  Double_t prodVtx_x, prodVtx_y, prodVtx_z;
  Int_t    prodVtxNTrack;
  Double_t kkVtx_x, kkVtx_y, kkVtx_z; // K18 x forward-proton vertex from input
  Int_t    kkVtxValid;
  Double_t multiVtx_x, multiVtx_y, multiVtx_z;
  Int_t    multiVtxNTrack;
  Int_t    multiVtxValid;
  Double_t vtxDelta_x, vtxDelta_y, vtxDelta_z, vtxDeltaMag;
  Int_t    vtxStatus; // 0=both invalid, 1=kk only, 2=multi only, 3=both
  Int_t    dcaReferenceFallback; // 1 if role DCA uses multi because kk invalid

  Int_t ntTpc;
  Int_t GFntTpc;

  // per-track record, indexed as the TPC tracks
  std::vector<Int_t>    trk_role;
  std::vector<Int_t>    trk_pidClass;
  std::vector<Int_t>    trk_pidbits;
  std::vector<Int_t>    trk_charge;
  std::vector<Int_t>    trk_nclust;
  std::vector<Int_t>    trk_gffit;
  std::vector<Int_t>    trk_inside;
  std::vector<Int_t>    trk_htofReached;
  std::vector<Int_t>    trk_isLambdaCandidateDaughter;
  std::vector<Int_t>    trk_isLambdaSignalDaughter;
  std::vector<Int_t>    trk_isLambdaSidebandDaughter;
  std::vector<Int_t>    trk_isK0Daughter;
  std::vector<Double_t> trk_chisqr;
  std::vector<Double_t> trk_dEdx;
  std::vector<Double_t> trk_mom;    // at the primary vertex
  std::vector<Double_t> trk_mom_x;
  std::vector<Double_t> trk_mom_y;
  std::vector<Double_t> trk_mom_z;
  std::vector<Double_t> trk_theta;  // deg, lab
  std::vector<Double_t> trk_phi;    // deg, lab
  std::vector<Double_t> trk_pstar;  // GeV/c, missing-system rest frame
  std::vector<Double_t> trk_costStar;
  std::vector<Double_t> trk_dcaPrimary; // mm; nominal = trk_dcaKK (or multi fallback)
  std::vector<Double_t> trk_dcaKK;      // mm, distance to K18×forward vertex
  std::vector<Double_t> trk_dcaMulti;  // mm, distance to robust multitrack vertex
  std::vector<Double_t> trk_pParallelQ;
  std::vector<Double_t> trk_pTransverseQ;
  std::vector<Double_t> trk_cosThetaQ;
  std::vector<Double_t> trk_pstarProtonHyp; // ambiguous PID: proton mass hypothesis
  std::vector<Double_t> trk_pstarPionHyp;   // ambiguous PID: pion mass hypothesis
  std::vector<Double_t> trk_m2;
  std::vector<Double_t> trk_invbeta;
  std::vector<Double_t> trk_nsigma_proton;
  std::vector<Double_t> trk_nsigma_kaon;
  std::vector<Double_t> trk_nsigma_pion;
  std::vector<Double_t> trk_nsigma_electron;
  std::vector<Double_t> trk_nsigma_deutron;
  std::vector<Double_t> trk_nsigma_triton;
  std::vector<Double_t> trk_nsigmaHtof_deutron;
  std::vector<Double_t> trk_nsigmaHtof_triton;
  std::vector<Double_t> trk_nsigmaHtof_proton;
  std::vector<Double_t> trk_nsigmaHtof_kaon;
  std::vector<Double_t> trk_nsigmaHtof_pion;

  // multiplicity ladder, one entry per working point of topodef::kNThr
  // Legacy coarse counts (kept for backward compatibility)
  std::vector<Int_t> nchRaw;
  std::vector<Int_t> nchPrompt;
  std::vector<Int_t> nProton;
  std::vector<Int_t> nPiPlus;
  std::vector<Int_t> nPiMinus;
  std::vector<Int_t> nAmbiguous;
  std::vector<Int_t> nDisplaced;
  std::vector<Int_t> topoClass;

  // Full prompt/displaced PID multiplicities (visible-exclusive basis)
  std::vector<Int_t> nPromptProton;
  std::vector<Int_t> nPromptPiPlus;
  std::vector<Int_t> nPromptPiMinus;
  std::vector<Int_t> nPromptKMinus;
  std::vector<Int_t> nPromptElectron;
  std::vector<Int_t> nPromptAmbiguous;
  std::vector<Int_t> nPromptUnknown;
  std::vector<Int_t> nDisplacedProton;
  std::vector<Int_t> nDisplacedPiPlus;
  std::vector<Int_t> nDisplacedPiMinus;
  std::vector<Int_t> nDisplacedKMinus;
  std::vector<Int_t> nDisplacedElectron;
  std::vector<Int_t> nDisplacedAmbiguous;
  std::vector<Int_t> nDisplacedUnknown;
  std::vector<Int_t> nPromptDeuteron;
  std::vector<Int_t> nPromptTriton;
  std::vector<Int_t> nPromptHeavyAmbiguous;
  std::vector<Int_t> nDisplacedDeuteron;
  std::vector<Int_t> nDisplacedTriton;
  std::vector<Int_t> nDisplacedHeavyAmbiguous;
  std::vector<Int_t> nOutsideUsed;
  std::vector<Int_t> nBelowThreshold;
  std::vector<Int_t> nNoMomentum;
  std::vector<Int_t> nPromptClosureDiff;
  std::vector<Int_t> nDisplacedClosureDiff;

  // Visible-exclusive signature and semi-exclusive tags (per WP)
  std::vector<ULong64_t> exclusiveVisibleKey;
  std::vector<ULong64_t> exclusiveVisibleKeyLambdaHypothesis;
  std::vector<ULong64_t> semiTagBits;

  // Lambda, best candidate of the input tree
  Int_t    lFlag;
  Int_t    lIsSignal;
  Int_t    lIsSideband;
  Int_t    lambdaRegion; // 0 none/outside, 1 signal, 2 sideband
  Int_t    nLambdaSignal;
  Int_t    nLambdaSideband;
  Double_t lMass;
  Double_t lMom, lMom_x, lMom_y, lMom_z;
  Double_t lDecayVtx_x, lDecayVtx_y, lDecayVtx_z;
  Double_t lDecayLen;
  Double_t lPPiDist;
  Double_t lPPiAngle;
  Double_t lPStar;
  std::vector<Int_t> lDaughterId;
  Double_t pLambdaDaughters; // |p_p + p_pi-| of Lambda daughters

  // All in-producer V0 candidates (mass-blind quality ranking)
  Int_t nLamCand;
  std::vector<Double_t> lamCand_mass;
  std::vector<Int_t>    lamCand_idPos;
  std::vector<Int_t>    lamCand_idNeg;
  std::vector<Double_t> lamCand_vtx_x;
  std::vector<Double_t> lamCand_vtx_y;
  std::vector<Double_t> lamCand_vtx_z;
  std::vector<Double_t> lamCand_dauDca;
  std::vector<Double_t> lamCand_decayLen;
  std::vector<Double_t> lamCand_cosPoint;
  std::vector<Double_t> lamCand_mom;
  std::vector<Double_t> lamCand_mom_x;
  std::vector<Double_t> lamCand_mom_y;
  std::vector<Double_t> lamCand_mom_z;
  std::vector<Double_t> lamCand_alpha;
  std::vector<Double_t> lamCand_qT;
  std::vector<Double_t> lamCand_quality;
  std::vector<Int_t>    lamCand_region; // 0 none, 1 signal, 2 left SB, 3 right SB
  std::vector<Int_t>    lamCand_fiducial;
  Int_t lamBestId; // index of selected candidate (-1 if none)

  // K0S reconstructed in this producer
  Int_t    k0Flag;
  Int_t    nK0S;
  Double_t k0Mass;
  Double_t k0Mom, k0Mom_x, k0Mom_y, k0Mom_z;
  std::vector<Int_t> k0DaughterId;
  Double_t pK0Daughters; // |p_pi+ + p_pi-| of K0S daughters

  Int_t nK0Cand;
  std::vector<Double_t> k0Cand_mass;
  std::vector<Int_t>    k0Cand_idPos;
  std::vector<Int_t>    k0Cand_idNeg;
  std::vector<Double_t> k0Cand_vtx_x;
  std::vector<Double_t> k0Cand_vtx_y;
  std::vector<Double_t> k0Cand_vtx_z;
  std::vector<Double_t> k0Cand_dauDca;
  std::vector<Double_t> k0Cand_decayLen;
  std::vector<Double_t> k0Cand_cosPoint;
  std::vector<Double_t> k0Cand_mom;
  std::vector<Double_t> k0Cand_mom_x;
  std::vector<Double_t> k0Cand_mom_y;
  std::vector<Double_t> k0Cand_mom_z;
  std::vector<Double_t> k0Cand_alpha;
  std::vector<Double_t> k0Cand_qT;
  std::vector<Double_t> k0Cand_quality;
  std::vector<Int_t>    k0Cand_region;
  std::vector<Int_t>    k0Cand_fiducial;
  Int_t k0BestId;

  Int_t v0OverlapFlag;
  std::vector<Int_t> v0OverlapTrackIds;
  Double_t sidebandScaleLambda; // FlatSidebandScale constant written per event
  Double_t sidebandScaleK0S;
  Int_t lPPiAngleValid; // 0 if lPPiAngle unset and not computed from daughters

  // Visible energy flow (exclude selected V0 daughters)
  Int_t visibleNetCharge;
  Double_t visibleScalarPSum;
  Double_t visibleVectorPSum_x;
  Double_t visibleVectorPSum_y;
  Double_t visibleVectorPSum_z;
  Double_t visibleVectorPSum;
  Double_t visiblePParallelQ;
  Double_t visiblePTransverseQ;
  Double_t visibleEnergyPionHyp;
  Double_t visibleEnergyProtonHyp;

  // Recoil observables: M_recoil = |P_X - sum P_detected| (not nuclear mass ID)
  Double_t mRecoilLambda;
  Double_t mRecoilLambdaP;
  Double_t mRecoilLambdaPim;
  Double_t mRecoilK0S;
  Double_t mRecoilPromptPim;
  Double_t mRecoilVisible;
  Double_t pRecoilLambda;
  Double_t pRecoilLambdaP;
  Double_t pRecoilLambdaPim;
  Double_t pRecoilK0S;
  Double_t pRecoilVisible;
  Double_t eRecoilLambda;
  Double_t eRecoilVisible;

  // Lambda + additional hadron: all pairs (vectors); scalars = best-by-max-p
  std::vector<Int_t>    addProtonIdAll;
  std::vector<Double_t> mLambdaPAll;
  std::vector<Double_t> pLambdaPAll;
  std::vector<Double_t> cosThetaLambdaPAll;
  std::vector<Double_t> mRecoilLambdaPAll;
  std::vector<Double_t> cosThetaLambdaPQAll;
  std::vector<Int_t>    addPimIdAll;
  std::vector<Double_t> mLambdaPimAll;
  std::vector<Double_t> pLambdaPimAll;
  std::vector<Double_t> cosThetaLambdaPimAll;
  std::vector<Double_t> mRecoilLambdaPimAll;
  std::vector<Double_t> cosThetaLambdaPimQAll;
  std::vector<Int_t>    addPipIdAll;
  std::vector<Double_t> mLambdaPipAll;
  std::vector<Double_t> pLambdaPipAll;
  std::vector<Double_t> cosThetaLambdaPipAll;

  // Lambda + additional proton / pi- (scalar best-by-max-p retained)
  Int_t    nAddProton;
  Int_t    addProtonId;
  Double_t addProtonMom;
  Double_t mLambdaP;
  Double_t cosThetaLambdaP;
  Double_t pLambdaP; // |p_Lambda + p_add|
  Int_t    nAddPim;
  Int_t    addPimId;
  Double_t addPimMom;
  Double_t mLambdaPim;
  Double_t cosThetaLambdaPim;
  Double_t pLambdaPim; // |p_Lambda + p_pi-|

  // forward-proton parentage test: M(p_forward, pi-)
  std::vector<Double_t> mFwdPPim;
  Int_t    lFwdFlag;      // any pi- pair evaluated
  Int_t    lFwdIsSignal;  // best M in narrow Lambda window
  Int_t    lFwdPionId;
  Double_t lFwdMass;
  Double_t lFwdPimMom;
  Double_t fwdLamDca;     // DCA of forward-p + pi- vertex (NaN if not found)
  Double_t fwdLamVtx_x;
  Double_t fwdLamVtx_y;
  Double_t fwdLamVtx_z;

  // reproduction of the GenfitQFKaon K- escape selection
  Int_t    kmIncFlag;
  Int_t    escFlag;
  Double_t escKmMom;
  Double_t escKmM2;

  void clear( void )
  {
    runnum = 0; evnum = 0; evStatus = 0;
    trigA = 0; trigB = 0; PSfacTrigA = 0; PSfacTrigB = 0;

    roleClassVersion = topodef::kRoleClassVersion;
    exclusiveVisibleVersion = topodef::kExclusiveVisibleVersion;
    semiTagVersion = topodef::kSemiTagVersion;
    v0CandidateVersion = topodef::kV0CandidateVersion;
    kinematicsVersion = topodef::kKinematicsVersion;

    BKaon = qnan; MissMassNucl = qnan; thetaKP = qnan; qTransfer = qnan;
    pBeam = qnan; pBeam_x = qnan; pBeam_y = qnan; pBeam_z = qnan;
    pFwd = qnan; pFwd_x = qnan; pFwd_y = qnan; pFwd_z = qnan;
    thetaFwd = qnan; phiFwd = qnan;
    PX_x = qnan; PX_y = qnan; PX_z = qnan; PX_E = qnan;

    prodVtx_x = qnan; prodVtx_y = qnan; prodVtx_z = qnan; prodVtxNTrack = 0;
    kkVtx_x = qnan; kkVtx_y = qnan; kkVtx_z = qnan; kkVtxValid = 0;
    multiVtx_x = qnan; multiVtx_y = qnan; multiVtx_z = qnan;
    multiVtxNTrack = 0; multiVtxValid = 0;
    vtxDelta_x = qnan; vtxDelta_y = qnan; vtxDelta_z = qnan; vtxDeltaMag = qnan;
    vtxStatus = 0; dcaReferenceFallback = 0;

    ntTpc = 0; GFntTpc = 0;

    trk_role.clear(); trk_pidClass.clear(); trk_pidbits.clear();
    trk_charge.clear(); trk_nclust.clear(); trk_gffit.clear();
    trk_inside.clear(); trk_htofReached.clear();
    trk_isLambdaCandidateDaughter.clear();
    trk_isLambdaSignalDaughter.clear();
    trk_isLambdaSidebandDaughter.clear();
    trk_isK0Daughter.clear();
    trk_chisqr.clear(); trk_dEdx.clear();
    trk_mom.clear(); trk_mom_x.clear(); trk_mom_y.clear(); trk_mom_z.clear();
    trk_theta.clear(); trk_phi.clear();
    trk_pstar.clear(); trk_costStar.clear();
    trk_dcaPrimary.clear(); trk_dcaKK.clear(); trk_dcaMulti.clear();
    trk_pParallelQ.clear(); trk_pTransverseQ.clear(); trk_cosThetaQ.clear();
    trk_pstarProtonHyp.clear(); trk_pstarPionHyp.clear();
    trk_m2.clear(); trk_invbeta.clear();
    trk_nsigma_proton.clear(); trk_nsigma_kaon.clear(); trk_nsigma_pion.clear();
    trk_nsigma_electron.clear(); trk_nsigma_deutron.clear(); trk_nsigma_triton.clear();
    trk_nsigmaHtof_proton.clear(); trk_nsigmaHtof_kaon.clear();
    trk_nsigmaHtof_pion.clear(); trk_nsigmaHtof_deutron.clear();
    trk_nsigmaHtof_triton.clear();

    nchRaw.assign(topodef::kNThr, 0);
    nchPrompt.assign(topodef::kNThr, 0);
    nProton.assign(topodef::kNThr, 0);
    nPiPlus.assign(topodef::kNThr, 0);
    nPiMinus.assign(topodef::kNThr, 0);
    nAmbiguous.assign(topodef::kNThr, 0);
    nDisplaced.assign(topodef::kNThr, 0);
    topoClass.assign(topodef::kNThr, -1);

    nPromptProton.assign(topodef::kNThr, 0);
    nPromptPiPlus.assign(topodef::kNThr, 0);
    nPromptPiMinus.assign(topodef::kNThr, 0);
    nPromptKMinus.assign(topodef::kNThr, 0);
    nPromptElectron.assign(topodef::kNThr, 0);
    nPromptAmbiguous.assign(topodef::kNThr, 0);
    nPromptUnknown.assign(topodef::kNThr, 0);
    nDisplacedProton.assign(topodef::kNThr, 0);
    nDisplacedPiPlus.assign(topodef::kNThr, 0);
    nDisplacedPiMinus.assign(topodef::kNThr, 0);
    nDisplacedKMinus.assign(topodef::kNThr, 0);
    nDisplacedElectron.assign(topodef::kNThr, 0);
    nDisplacedAmbiguous.assign(topodef::kNThr, 0);
    nDisplacedUnknown.assign(topodef::kNThr, 0);
    nPromptDeuteron.assign(topodef::kNThr, 0);
    nPromptTriton.assign(topodef::kNThr, 0);
    nPromptHeavyAmbiguous.assign(topodef::kNThr, 0);
    nDisplacedDeuteron.assign(topodef::kNThr, 0);
    nDisplacedTriton.assign(topodef::kNThr, 0);
    nDisplacedHeavyAmbiguous.assign(topodef::kNThr, 0);
    nOutsideUsed.assign(topodef::kNThr, 0);
    nBelowThreshold.assign(topodef::kNThr, 0);
    nNoMomentum.assign(topodef::kNThr, 0);
    nPromptClosureDiff.assign(topodef::kNThr, 0);
    nDisplacedClosureDiff.assign(topodef::kNThr, 0);

    exclusiveVisibleKey.assign(topodef::kNThr, 0);
    exclusiveVisibleKeyLambdaHypothesis.assign(topodef::kNThr, 0);
    semiTagBits.assign(topodef::kNThr, 0);

    lFlag = 0; lIsSignal = 0; lIsSideband = 0;
    lambdaRegion = 0; nLambdaSignal = 0; nLambdaSideband = 0;
    lMass = qnan;
    lMom = qnan; lMom_x = qnan; lMom_y = qnan; lMom_z = qnan;
    lDecayVtx_x = qnan; lDecayVtx_y = qnan; lDecayVtx_z = qnan;
    lDecayLen = qnan; lPPiDist = qnan; lPPiAngle = qnan; lPStar = qnan;
    lDaughterId.clear();
    pLambdaDaughters = qnan;

    nLamCand = 0; lamBestId = -1;
    lamCand_mass.clear(); lamCand_idPos.clear(); lamCand_idNeg.clear();
    lamCand_vtx_x.clear(); lamCand_vtx_y.clear(); lamCand_vtx_z.clear();
    lamCand_dauDca.clear(); lamCand_decayLen.clear(); lamCand_cosPoint.clear();
    lamCand_mom.clear(); lamCand_mom_x.clear(); lamCand_mom_y.clear(); lamCand_mom_z.clear();
    lamCand_alpha.clear(); lamCand_qT.clear(); lamCand_quality.clear();
    lamCand_region.clear(); lamCand_fiducial.clear();

    k0Flag = 0; nK0S = 0;
    k0Mass = qnan; k0Mom = qnan; k0Mom_x = qnan; k0Mom_y = qnan; k0Mom_z = qnan;
    k0DaughterId.clear();
    pK0Daughters = qnan;

    nK0Cand = 0; k0BestId = -1;
    k0Cand_mass.clear(); k0Cand_idPos.clear(); k0Cand_idNeg.clear();
    k0Cand_vtx_x.clear(); k0Cand_vtx_y.clear(); k0Cand_vtx_z.clear();
    k0Cand_dauDca.clear(); k0Cand_decayLen.clear(); k0Cand_cosPoint.clear();
    k0Cand_mom.clear(); k0Cand_mom_x.clear(); k0Cand_mom_y.clear(); k0Cand_mom_z.clear();
    k0Cand_alpha.clear(); k0Cand_qT.clear(); k0Cand_quality.clear();
    k0Cand_region.clear(); k0Cand_fiducial.clear();

    v0OverlapFlag = 0; v0OverlapTrackIds.clear();
    sidebandScaleLambda = topodef::kLambdaSidebandScale;
    sidebandScaleK0S = topodef::kK0SSidebandScale;
    lPPiAngleValid = 0;

    visibleNetCharge = 0; visibleScalarPSum = qnan;
    visibleVectorPSum_x = qnan; visibleVectorPSum_y = qnan;
    visibleVectorPSum_z = qnan; visibleVectorPSum = qnan;
    visiblePParallelQ = qnan; visiblePTransverseQ = qnan;
    visibleEnergyPionHyp = qnan; visibleEnergyProtonHyp = qnan;

    mRecoilLambda = qnan; mRecoilLambdaP = qnan; mRecoilLambdaPim = qnan;
    mRecoilK0S = qnan; mRecoilPromptPim = qnan; mRecoilVisible = qnan;
    pRecoilLambda = qnan; pRecoilLambdaP = qnan; pRecoilLambdaPim = qnan;
    pRecoilK0S = qnan; pRecoilVisible = qnan;
    eRecoilLambda = qnan; eRecoilVisible = qnan;

    addProtonIdAll.clear(); mLambdaPAll.clear(); pLambdaPAll.clear();
    cosThetaLambdaPAll.clear(); mRecoilLambdaPAll.clear(); cosThetaLambdaPQAll.clear();
    addPimIdAll.clear(); mLambdaPimAll.clear(); pLambdaPimAll.clear();
    cosThetaLambdaPimAll.clear(); mRecoilLambdaPimAll.clear(); cosThetaLambdaPimQAll.clear();
    addPipIdAll.clear(); mLambdaPipAll.clear(); pLambdaPipAll.clear();
    cosThetaLambdaPipAll.clear();

    nAddProton = 0; addProtonId = -1; addProtonMom = qnan;
    mLambdaP = qnan; cosThetaLambdaP = qnan; pLambdaP = qnan;
    nAddPim = 0; addPimId = -1; addPimMom = qnan;
    mLambdaPim = qnan; cosThetaLambdaPim = qnan; pLambdaPim = qnan;

    mFwdPPim.clear();
    lFwdFlag = 0; lFwdIsSignal = 0; lFwdPionId = -1;
    lFwdMass = qnan; lFwdPimMom = qnan;
    fwdLamDca = qnan; fwdLamVtx_x = qnan; fwdLamVtx_y = qnan; fwdLamVtx_z = qnan;

    kmIncFlag = 0; escFlag = 0; escKmMom = qnan; escKmM2 = qnan;
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
  Double_t kp_phi;
  Double_t kp_theta;      
  
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
  Int_t ntKuramaCandidate; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)  
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

  Int_t IncFlag;
  Int_t EscFlag;
  Int_t PSfacTrigA;
  Int_t PSfacTrigB;  

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
    ntKuramaCandidate = 0; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)     
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

    IncFlag = false;
    EscFlag = false;
    PSfacTrigA = 0;
    PSfacTrigB = 0;

    kp_phi = qnan;
    kp_theta = qnan;    
    
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
  TTreeReaderValue<Int_t>* ntKuramaCandidate; //Numer of tracks which are kurama track candidates(before TPCKurama tracking)    
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

  // Lambda of the DstE42 schema: one best candidate per event
  TTreeReaderValue<Double_t>* lmass;
  TTreeReaderValue<Double_t>* ldecayvtx_x;
  TTreeReaderValue<Double_t>* ldecayvtx_y;
  TTreeReaderValue<Double_t>* ldecayvtx_z;
  TTreeReaderValue<Double_t>* lmom;
  TTreeReaderValue<Double_t>* lmom_x;
  TTreeReaderValue<Double_t>* lmom_y;
  TTreeReaderValue<Double_t>* lmom_z;
  TTreeReaderValue<Double_t>* ppi_dist;
  TTreeReaderValue<Double_t>* ppiangle;
  TTreeReaderValue<std::vector<Int_t>>* ldecays_id;
};

namespace root
{
Event  event;
DstG4  dstg4;
Topo   topo;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
TTree *treetpc;
std::map<ULong64_t, Long64_t> gExclKeyCount; // nominal WP exclusive key frequency
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
  topo.clear();

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

  if (!TFileCont[kOutFile] || !TFileCont[kOutFile]->IsOpen()) {
      std::cerr << "!!! DstOpen: Failed to create output file: " << arg[kOutFile] << std::endl;
      return false;
  }
  TFileCont[kOutFile]->cd();

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
  static const auto DeuteronMass = pdg::DeutronMass();
  static const Double_t TritonMass = 3.0160492 * TGeoUnit::amu_c2;
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
  event.PSfacTrigA = psTrigA;
  event.PSfacTrigB = psTrigB;
  
  topo.runnum = event.runnum;
  topo.evnum = event.evnum;
  topo.PSfacTrigA = psTrigA;
  topo.PSfacTrigB = psTrigB;
  topo.trigA = trigA ? 1 : 0;
  topo.trigB = trigB ? 1 : 0;

  // Only events reaching the end of the selection below are written out, so
  // that the tree itself defines the denominator of every topology fraction.
  if(!trigA&&!trigB) return false;
  topo.evStatus |= (1<<0);

  for(int ips=0; ips<psfac; ips++){
    HF1( 1, event.status ); // debug 0    
  }
  event.status++;
  
  if(event.nKK != 1) return false;
  topo.evStatus |= (1<<1);
  double BE = 0.;
  double thetaTPC = 0.;  
  for(Int_t iKK=0; iKK<event.nKK; iKK++){
    BE = event.MissMassNuclCorrDETPC[iKK] - KaonMass - Boron11Mass;
    thetaTPC = event.thetaTPC[0];
  }
  HF1(500,thetaTPC);
  event.ntKuramaCandidate = **src.ntKuramaCandidate;
  if(event.Pflag[0] != 1) return false;
  topo.evStatus |= (1<<2);
  {
    std::cout << __FILE__ << " " << __LINE__ << " ntKuramaCandidate:" << event.ntKuramaCandidate << std::endl;
    double vtx = event.vtx[0];
    double vty = event.vty[0];
    double vtz = event.vtz[0]; // = vertexZ - targetZ
    HF1(600,vtx); HF1(601,vty); HF1(602,vtz); 
    HF2(700,vtx,vty); HF2(701,vty,vtz); HF2(702,vtz,vtx); 
    for(int i=0; i<12; i++){
      if(double(i)*2.0<thetaTPC&&double(i+1)*2.0>thetaTPC) HF1(610+i,vtz);
      {
	TVector3 km(event.ubTPC[0], event.vbTPC[0], 1.);
	TVector3 kp(event.us[0], event.vs[0], 1.);
	double cos_kp = km*kp;
	cos_kp /= (km.Mag()*kp.Mag());
	double theta_kp = TMath::ACos(cos_kp)*TMath::RadToDeg();
	if(double(i)*2.0<theta_kp&&double(i+1)*2.0>theta_kp){
	  if(event.ntKuramaCandidate==0) HF1(800+i,vtz);
	  if(event.ntKuramaCandidate==1) HF1(820+i,vtz); 
	  if(event.ntKuramaCandidate==2) HF1(840+i,vtz);
	  if(event.ntKuramaCandidate>0)  HF1(860+i,vtz);	  
	}
      }
    }
    if(3.5<thetaTPC&&4.5>thetaTPC){
      HF1(630,vtx);
      HF1(631,vty);
      HF1(632,vtz);            
    }
  }
  if( event.isgoodTPCKurama.size()!=1 ) return false; 
  if( event.isgoodTPCKurama[0]!=1 ) return false; 
  topo.evStatus |= (1<<3);
  {
    double vtxtpc = event.vtxTPC[0];
    double vtytpc = event.vtyTPC[0];
    double vtztpc = event.vtzTPC[0];  // = vertexZ - targetZ
    HF1(650,vtxtpc); HF1(651,vtytpc); HF1(652,vtztpc); 
    HF2(750,vtxtpc,vtytpc); HF2(751,vtytpc,vtztpc); HF2(752,vtztpc,vtxtpc); 
    for(int i=0; i<12; i++){ 
      if(double(i)*2.0<thetaTPC&&double(i+1)*2.0>thetaTPC) HF1(660+i,vtztpc);
      // if(double(i)*2.0<thetaTPC&&double(i+1)*2.0>thetaTPC){
      // 	if(event.ntKuramaCandidate==0) HF1(800+i,vtz); 
      // 	if(event.ntKuramaCandidate==1) HF1(820+i,vtz); 
      // 	if(event.ntKuramaCandidate==2) HF1(840+i,vtz);
      // 	if(event.ntKuramaCandidate>0)  HF1(860+i,vtz); 	
      // }
    }
    if(3.5<thetaTPC&&4.5>thetaTPC){
      HF1(680,vtxtpc);
      HF1(681,vtytpc);
      HF1(682,vtztpc);            
    }    
  }
  if( event.insideTPC[0] != 1) return false;  
  topo.evStatus |= (1<<4);
  if(KPEvent && event.Pflag[0] != 1){
    if(event.pflagTPCKurama[0]!=1) return false;    
    return false; //precut with Kurama tracking
  }  
  if( !(event.runnum >= 5641 && event.runnum <= 5666) ){// not CH2
    if(!(thetaTPC>minThetaKP && thetaTPC<maxThetaKP)) return false;
  } else {
    if(!(thetaTPC>minThetaKPCH2 && thetaTPC<maxThetaKPCH2)) return false;
  }  
  topo.evStatus |= (1<<5);
  event.IncFlag = true;
  // beam
  TLorentzVector LvRcTPC;
  TVector3 km_unit = TVector3(event.utgtK18[0], event.vtgtK18[0], 1.).Unit();
  TVector3 km_momTPC = km_unit*event.pK18[0];
  // scat
  TVector3 kp_unit = TVector3(event.usTPC[0], event.vsTPC[0], 1.).Unit();
  TVector3 kp_momTPC = kp_unit*event.pCorrDETPC[0];
  event.kp_phi = kp_momTPC.Phi();
  event.kp_theta = kp_momTPC.Theta();
  
  // missing
  TVector3 miss_momTPC = km_momTPC - kp_momTPC;  

  TLorentzVector LvKmTPC(km_momTPC, TMath::Hypot(km_momTPC.Mag(), KaonMass));
  TLorentzVector LvScatPTPC(kp_momTPC, TMath::Hypot(kp_momTPC.Mag(), ProtonMass));
  TLorentzVector LvCTPC(0., 0., 0., Carbon12Mass);
  TLorentzVector LvPTPC(0., 0., 0., ProtonMass);
  TLorentzVector LvScatKmTPC = LvKmTPC + LvPTPC - LvScatPTPC;
  LvRcTPC = LvKmTPC + LvCTPC - LvScatPTPC;

  double mm_12CTPC = LvRcTPC.M();
  //double binding_energyTPC = Boron11Mass + KaonMass - (mm_12CTPC - 0.120); //GeV/c2
  double binding_energyTPC = Boron11Mass + KaonMass - mm_12CTPC; //GeV/c2
  event.BETPC[0] = binding_energyTPC; //MeV/c2
  
  //event.MissMassCorrDETPC[0] = event.MissMassCorrDETPC[0] - 0.120 ;
  event.MissMassCorrDETPC[0] = event.MissMassCorrDETPC[0];  
  double missmass = event.MissMassCorrDETPC[0];

  // LvRcTPC is the four-momentum of the missing system X; the rest frame of X
  // is what Analysis 8 of the proposal asks the particle spectra to be shown in.
  const TLorentzVector LvMissX = LvRcTPC;
  const TVector3 boostToX = -LvMissX.BoostVector();

  topo.BKaon = binding_energyTPC;
  topo.MissMassNucl = mm_12CTPC;
  topo.thetaKP = thetaTPC;
  topo.qTransfer = miss_momTPC.Mag();
  topo.pBeam = km_momTPC.Mag();
  topo.pBeam_x = km_momTPC.x();
  topo.pBeam_y = km_momTPC.y();
  topo.pBeam_z = km_momTPC.z();
  topo.pFwd = kp_momTPC.Mag();
  topo.pFwd_x = kp_momTPC.x();
  topo.pFwd_y = kp_momTPC.y();
  topo.pFwd_z = kp_momTPC.z();
  topo.thetaFwd = kp_momTPC.Theta()*TMath::RadToDeg();
  topo.phiFwd = kp_momTPC.Phi()*TMath::RadToDeg();
  topo.PX_x = LvMissX.Px();
  topo.PX_y = LvMissX.Py();
  topo.PX_z = LvMissX.Pz();
  topo.PX_E = LvMissX.E();
  topo.kkVtx_x = event.vtxTPC[0];
  topo.kkVtx_y = event.vtyTPC[0];
  topo.kkVtx_z = event.vtzTPC[0];
  
  if(missmass>minMM&&missmass<maxMM){
    for(int ips=0; ips<psfac; ips++) HF1( 1, event.status ); // debug 1
  }
  event.status++;

  if(trigA){    
    for(int ips=0; ips<event.PSfacTrigA; ips++){
      HF1(3900,event.MissMassCorrDETPC[0]);
      HF1(13900,-event.BETPC[0]);
    }
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
  
  int ntTpc = **src.ntTpc;  
  topo.ntTpc = ntTpc;
  // ntTpc==0 is written out rather than dropped: it is a legitimate member of
  // the denominator, flagged by evStatus bit 6 being unset.
  if( ntTpc == 0 )
    return true;
  topo.evStatus |= (1<<6);
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
  topo.GFntTpc = GFntTpc;
  if(GFntTpc!=event.ntTpc){
    std::cout<<"# of Tracks in Genfit Track Container != # of TPC Tracks"<<std::endl;
    std::cout << "ntTpc:" << ntTpc << " GFntTpc:" << GFntTpc << std::endl;
    GFtrackCont.Clear();
    return true;
  }
  topo.evStatus |= (1<<7);
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
  Int_t ntrack_multi = 0;
  Double_t x0[100] = {0};
  Double_t y0[100] = {0};
  Double_t u0[100] = {0};
  Double_t v0[100] = {0};
  Double_t mx0[100] = {0};
  Double_t my0[100] = {0};
  Double_t mu0[100] = {0};
  Double_t mv0[100] = {0};
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
	// HF1( genfitHid+8, event.GFlayer[igf][ihit]); 
	// HF1( genfitHid+10, event.GFresidual_x[igf][ihit]);
	// HF1( genfitHid+11, event.GFresidual_y[igf][ihit]);
	// HF1( genfitHid+12, event.GFresidual_z[igf][ihit]);
	// HF1( genfitHid+13, event.GFresidual_p[igf][ihit]);
	// HF1( genfitHid+14, event.GFresidual_px[igf][ihit]);
	// HF1( genfitHid+15, event.GFresidual_py[igf][ihit]);
	// HF1( genfitHid+16, event.GFresidual_pz[igf][ihit]);
	// HF1( genfitHid+1000*(layer+1), event.GFresidual_x[igf][ihit]);
	// HF1( genfitHid+1000*(layer+1)+1, event.GFresidual_y[igf][ihit]);
	// HF1( genfitHid+1000*(layer+1)+2, event.GFresidual_z[igf][ihit]);
	// HF1( genfitHid+1000*(layer+1)+3, event.GFresidual_p[igf][ihit]);
	// HF1( genfitHid+1000*(layer+1)+4, event.GFresidual_px[igf][ihit]);
	// HF1( genfitHid+1000*(layer+1)+5, event.GFresidual_py[igf][ihit]);
	// HF1( genfitHid+1000*(layer+1)+6, event.GFresidual_pz[igf][ihit]);
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
	if(GFtrackCont.ExtrapolateToTargetCenter(igf, post, momt, lentgt, toftgt)
	   && ntrack_intarget < 100){
	  x0[ntrack_intarget] = post.x();
	  y0[ntrack_intarget] = post.y();
	  u0[ntrack_intarget] = momt.x()/momt.z();
	  v0[ntrack_intarget] = momt.y()/momt.z();
	  ntrack_intarget++;
	}
	// Robust multitrack vertex: exclude beam-like tracks only (V0 daughters unknown here).
	if(!event.isK18[igf] && !event.isKurama[igf]
	   && !event.isBeam[igf] && !event.isAccidental[igf]
	   && ntrack_multi < 100
	   && GFtrackCont.ExtrapolateToTargetCenter(igf, post, momt, lentgt, toftgt)){
	  mx0[ntrack_multi] = post.x();
	  my0[ntrack_multi] = post.y();
	  mu0[ntrack_multi] = momt.x()/momt.z();
	  mv0[ntrack_multi] = momt.y()/momt.z();
	  ntrack_multi++;
	}
      } else {
	event.GFinside[igf] = 0;
      }
    } else {
      event.GFnhtrack[igf] = 0;
      event.GFpdgcode[igf] = -9999;
      event.GFm2[igf] = TMath::QuietNaN();
        
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

  // Second pass: dual-vertex construction.
  // kkVertex = K18×forward from input (nominal DCA reference for role WP).
  // multiVtx = robust fit excluding beam-like tracks (fallback + QA).
  TVector3 prodVertex = Kinematics::MultitrackVertex(ntrack_intarget, x0, y0, u0, v0);
  TVector3 kkVertex(topo.kkVtx_x, topo.kkVtx_y, topo.kkVtx_z);
  topo.kkVtxValid = (!TMath::IsNaN(kkVertex.x())) ? 1 : 0;
  Int_t multiNUsed = 0;
  Bool_t multiValidFlag = false;
  TVector3 multiVertex = kptopo::RobustMultitrackVertex(
    ntrack_multi, mx0, my0, mu0, mv0, 15., 3, multiNUsed, multiValidFlag);
  topo.multiVtx_x = multiVertex.x();
  topo.multiVtx_y = multiVertex.y();
  topo.multiVtx_z = multiVertex.z();
  topo.multiVtxNTrack = multiNUsed;
  topo.multiVtxValid = multiValidFlag ? 1 : 0;
  if(topo.kkVtxValid && topo.multiVtxValid){
    topo.vtxDelta_x = kkVertex.x() - multiVertex.x();
    topo.vtxDelta_y = kkVertex.y() - multiVertex.y();
    topo.vtxDelta_z = kkVertex.z() - multiVertex.z();
    topo.vtxDeltaMag = TMath::Sqrt(topo.vtxDelta_x*topo.vtxDelta_x
				   + topo.vtxDelta_y*topo.vtxDelta_y
				   + topo.vtxDelta_z*topo.vtxDelta_z);
    topo.vtxStatus = 3;
  } else if(topo.kkVtxValid){
    topo.vtxStatus = 1;
  } else if(topo.multiVtxValid){
    topo.vtxStatus = 2;
  } else {
    topo.vtxStatus = 0;
  }
  TVector3 htofVertex = kkVertex;
  if(!topo.kkVtxValid){
    if(topo.multiVtxValid) htofVertex = multiVertex;
    else htofVertex = TVector3(qnan, qnan, qnan);
  }
  event.GFntTpc_inside = ntrack_intarget;
  event.GFprodvtx_x = prodVertex.x();
  event.GFprodvtx_y = prodVertex.y();
  event.GFprodvtx_z = prodVertex.z();
  topo.prodVtx_x = prodVertex.x();
  topo.prodVtx_y = prodVertex.y();
  topo.prodVtx_z = prodVertex.z();
  topo.prodVtxNTrack = ntrack_intarget;

  for( Int_t igf=0; igf<GFntTpc; ++igf ){
    if( !event.GFfitstatus[igf] ) continue;
    if( event.GFinside[igf]!=1 ) continue;
    if( event.isBeam[igf]==1 || event.isK18[igf]==1 || event.isAccidental[igf]==1 ) continue;
    if(!topo.kkVtxValid && !topo.multiVtxValid) continue;
    Int_t repid=-1;
    Int_t hitid_htof; Double_t tof; Double_t len;
    TVector3 pos_htof; Double_t track2tgt_dist;
    Bool_t htofextrapvtx =
      GFtrackCont.TPCHTOFTrackMatching(igf, repid, htofVertex,
				       event.HtofSeg, event.posHtof,
				       hitid_htof, tof,
				       len, pos_htof, track2tgt_dist);
    if(!htofextrapvtx) continue;
    event.GFfromVtx[igf] = 1;
    event.GFtracklen[igf] = len;
    event.GFcalctof[igf] = tof;
    event.GFposx[igf] = pos_htof.x();
    event.GFposy[igf] = pos_htof.y();
    event.GFposz[igf] = pos_htof.z();
    event.GFsegHtof[igf] = event.HtofSeg[hitid_htof];
    event.GFtofHtof[igf] = event.tHtof[hitid_htof];
    event.GFposHtof[igf] = event.posHtof[hitid_htof];
    Double_t beta = len/event.tHtof[hitid_htof]/MathTools::C();
    event.GFinvbeta[igf] = 1./beta;
    event.GFm2[igf] = Kinematics::MassSquare(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
    event.nsigma_tritonHtof[igf] = Kinematics::HypTPCHTOFNsigmaTriton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
    event.nsigma_deutronHtof[igf] = Kinematics::HypTPCHTOFNsigmaDeutron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
    event.nsigma_protonHtof[igf] = Kinematics::HypTPCHTOFNsigmaProton(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
    event.nsigma_kaonHtof[igf] = Kinematics::HypTPCHTOFNsigmaKaon(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
    event.nsigma_pionHtof[igf] = Kinematics::HypTPCHTOFNsigmaPion(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
    event.nsigma_electronHtof[igf] = Kinematics::HypTPCHTOFNsigmaElectron(event.GFmom[igf][0], len, event.tHtof[hitid_htof]);
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
  Int_t numOutside=0;  
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
      if(event.isElectron[it]==1){
	if(event.charge[it]==1){
	  numEp++;
	  continue;
	} else {
	  numEm++;
	  continue;
	}
	continue;
      }
      if(event.isK18[it]==1) continue;
      if(event.isKurama[it]==1) continue;
      if(event.isBeam[it]==1) continue;
      if(event.isAccidental[it]==1) continue;
      if(event.GFinside[it]!=1){
	numOutside++;
	continue;
      }
      if(event.charge[it]==-1&&event.mom0[it]<0.5) HF1(15, -1.*event.GFm2[it]);
      if( event.charge[it]==-1 && event.mom0[it]<0.5 && event.GFinvbeta[it]>0.
	  && event.nsigma_kaon[it]>-3 && event.nsigma_kaon[it]<3 ) HF1(18,-1.*event.GFm2[it]);
      //if((event.pid[it]&2)==2 && event.charge[it]==-1){ //k-
      if(event.charge[it]==-1){ //k-
	if( event.nsigma_kaon[it]>mindEdxSigKaon && event.nsigma_kaon[it]<maxdEdxSigKaon ){ 
	  numKm++; 
	  idKm = it; 
	  continue; 
	}
      }
      if((event.pid[it]&4)==4 && (event.pid[it]&1)!=1 && event.charge[it]==1){ //proton 
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
	  // target_accidental_id_container.push_back(it); //Accidental beam on the target 
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
    std::cout << " debug " << __FILE__ << " " << __LINE__
	      << " TotalTrack:" << event.ntTpc << " ntK18:" << numK18 << " ntBeam:" << numBeam << " nuKurama:" << numKurama << " ntAcc:" << numAcc << " ntOutTarget:" << numOutside
	      << " numKm:" << numKm << " numPPip:" << numPPip << " numPip:" << numPip << " numPim:" << numPim << " numP:" << numP
	      << " numEm:" << numEm << " numEp:" << numEp << std::endl;
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

  if(kmflag_dedxpid&&trigB){
    for(int ips=0; ips<event.PSfacTrigB; ips++){
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
    
    if(trigB&&event.GFkmmom>0.01){      
      for(int ips=0; ips<event.PSfacTrigB; ips++){	
	HF1(20, event.GFkmmass2);
	if(event.GFkmmom<0.5) HF1(19, -1.*event.GFkmmass2);
	for(int i=0; i<5; i++){
	  if(double(i)*0.2<event.GFkmmom && double(i+1)*0.2>event.GFkmmom) HF1(121+i, event.GFkmmass2);
	}
	if(0.2<event.GFkmmom && event.GFkmmom<0.7) HF1(126, event.GFkmmass2); 
	if(0.7<event.GFkmmom) HF1(127, event.GFkmmass2); 
	for(int i=0; i<10; i++){ 
	  if(double(i)*0.1<event.GFkmmom && double(i+1)*0.1>event.GFkmmom && event.charge[id]<0) HF1(130+i, event.charge[id]*event.GFkmmass2);
	}		
	
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
	  {
	    double thetaTPC = event.thetaTPC[0];
	    double vtx = event.vtx[0];
	    double vty = event.vty[0];
	    double vtz = event.vtz[0];  // = vertexZ - targetZ
	    HF1(10600,vtx); HF1(10601,vty); HF1(10602,vtz); 
	    HF2(10700,vtx,vty); HF2(10701,vty,vtz); HF2(10702,vtz,vtx); 
	    for(int i=0; i<12; i++){
	      if(double(i)*2.0<thetaTPC&&double(i+1)*2.0>thetaTPC) HF1(10610+i,vtz);
	    }
	    double vtxtpc = event.vtxTPC[0];
	    double vtytpc = event.vtyTPC[0];
	    double vtztpc = event.vtzTPC[0]; // = vertexZ - targetZ
	    HF1(10650,vtxtpc); HF1(10651,vtytpc); HF1(10652,vtztpc); 
	    HF2(10750,vtxtpc,vtytpc); HF2(10751,vtytpc,vtztpc); HF2(10752,vtztpc,vtxtpc); 
	    for(int i=0; i<12; i++){
	      if(double(i)*2.0<thetaTPC&&double(i+1)*2.0>thetaTPC) HF1(10660+i,vtztpc);
	    }	    
	  }

	  event.EscFlag = true;
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
	    HF2(5100, bek, event.GFkmmom);
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

  //___________________________________________________________________________
  // Topology analysis
  //
  // Runs after the escape selection above so that isAccidental carries the same
  // reclassification, and so that escFlag reproduces GenfitQFKaon.cc exactly.

  topo.kmIncFlag = event.kminc ? 1 : 0;
  topo.escFlag = event.EscFlag ? 1 : 0;
  if(event.kmflag){
    topo.escKmMom = event.GFkmmom;
    topo.escKmM2 = event.GFkmmass2;
  }

  // Lambda from input DstE42 (legacy). In-producer FindVertex candidates are
  // preferred for nominal physics when available (see V0 block below).
  Int_t inputLFlag = 0;
  Double_t inputLMass = qnan, inputLMom = qnan;
  Double_t inputLMom_x = qnan, inputLMom_y = qnan, inputLMom_z = qnan;
  Double_t inputLDecayVtx_x = qnan, inputLDecayVtx_y = qnan, inputLDecayVtx_z = qnan;
  Double_t inputLPPiDist = qnan, inputLPPiAngle = qnan;
  std::vector<Int_t> lDaughter;
  if(**src.lflag == 1){
    inputLFlag = 1;
    inputLMass = **src.lmass;
    inputLMom = **src.lmom;
    inputLMom_x = **src.lmom_x;
    inputLMom_y = **src.lmom_y;
    inputLMom_z = **src.lmom_z;
    inputLDecayVtx_x = **src.ldecayvtx_x;
    inputLDecayVtx_y = **src.ldecayvtx_y;
    inputLDecayVtx_z = **src.ldecayvtx_z;
    inputLPPiDist = **src.ppi_dist;
    inputLPPiAngle = **src.ppiangle;
    lDaughter = **src.ldecays_id;
  }

  const TVector3 primaryRef = topo.kkVtxValid
    ? TVector3(topo.kkVtx_x, topo.kkVtx_y, topo.kkVtx_z)
    : (topo.multiVtxValid
       ? TVector3(topo.multiVtx_x, topo.multiVtx_y, topo.multiVtx_z)
       : prodVertex);
  topo.dcaReferenceFallback = (topo.kkVtxValid ? 0 : (topo.multiVtxValid ? 1 : 0));
  const TVector3 kkRef(topo.kkVtx_x, topo.kkVtx_y, topo.kkVtx_z);
  const TVector3 multiRef(topo.multiVtx_x, topo.multiVtx_y, topo.multiVtx_z);
  const TVector3 qLab = miss_momTPC;

  auto isInIdList = [](const std::vector<Int_t>& ids, Int_t it){
    for(size_t i=0; i<ids.size(); ++i) if(ids[i]==it) return true;
    return false;
  };
  auto repidFromPid = [](Int_t pidbit, Bool_t isKurama)->Int_t{
    if(isKurama) return 0;
    Int_t repid = 0;
    if((pidbit&1)==1) repid = 0;
    else if((pidbit&2)==2){
      if((1&pidbit)==1) repid += 1;
    } else if((pidbit&4)==4){
      Int_t flag = 1;
      for(Int_t i=0;i<2;i++){
	if((flag&pidbit)==flag) repid += 1;
	flag *= 2;
      }
    }
    return repid;
  };
  auto isProtonLike = [](Int_t cls){
    return cls==topodef::kPidProton || cls==topodef::kPidProtonOrPiPlus;
  };

  // Per-track record. Every TPC track is written; the role label says what it
  // is, so no track is silently removed from the sample.
  topo.trk_role.assign(ntTpc, topodef::kRoleUnknown);
  topo.trk_pidClass.assign(ntTpc, topodef::kPidUnknown);
  topo.trk_pidbits.assign(ntTpc, 0);
  topo.trk_charge.assign(ntTpc, 0);
  topo.trk_nclust.assign(ntTpc, 0);
  topo.trk_gffit.assign(ntTpc, 0);
  topo.trk_inside.assign(ntTpc, 0);
  topo.trk_htofReached.assign(ntTpc, 0);
  topo.trk_isLambdaCandidateDaughter.assign(ntTpc, 0);
  topo.trk_isLambdaSignalDaughter.assign(ntTpc, 0);
  topo.trk_isLambdaSidebandDaughter.assign(ntTpc, 0);
  topo.trk_isK0Daughter.assign(ntTpc, 0);
  topo.trk_chisqr.assign(ntTpc, qnan);
  topo.trk_dEdx.assign(ntTpc, qnan);
  topo.trk_mom.assign(ntTpc, qnan);
  topo.trk_mom_x.assign(ntTpc, qnan);
  topo.trk_mom_y.assign(ntTpc, qnan);
  topo.trk_mom_z.assign(ntTpc, qnan);
  topo.trk_theta.assign(ntTpc, qnan);
  topo.trk_phi.assign(ntTpc, qnan);
  topo.trk_pstar.assign(ntTpc, qnan);
  topo.trk_costStar.assign(ntTpc, qnan);
  topo.trk_dcaPrimary.assign(ntTpc, qnan);
  topo.trk_dcaKK.assign(ntTpc, qnan);
  topo.trk_dcaMulti.assign(ntTpc, qnan);
  topo.trk_pParallelQ.assign(ntTpc, qnan);
  topo.trk_pTransverseQ.assign(ntTpc, qnan);
  topo.trk_cosThetaQ.assign(ntTpc, qnan);
  topo.trk_pstarProtonHyp.assign(ntTpc, qnan);
  topo.trk_pstarPionHyp.assign(ntTpc, qnan);
  topo.trk_m2.assign(ntTpc, qnan);
  topo.trk_invbeta.assign(ntTpc, qnan);
  topo.trk_nsigma_proton.assign(ntTpc, qnan);
  topo.trk_nsigma_kaon.assign(ntTpc, qnan);
  topo.trk_nsigma_pion.assign(ntTpc, qnan);
  topo.trk_nsigma_electron.assign(ntTpc, qnan);
  topo.trk_nsigma_deutron.assign(ntTpc, qnan);
  topo.trk_nsigma_triton.assign(ntTpc, qnan);
  topo.trk_nsigmaHtof_proton.assign(ntTpc, qnan);
  topo.trk_nsigmaHtof_kaon.assign(ntTpc, qnan);
  topo.trk_nsigmaHtof_pion.assign(ntTpc, qnan);
  topo.trk_nsigmaHtof_deutron.assign(ntTpc, qnan);
  topo.trk_nsigmaHtof_triton.assign(ntTpc, qnan);

  for(Int_t it=0; it<ntTpc; ++it){
    topo.trk_charge[it] = event.charge[it];
    topo.trk_pidbits[it] = event.pid[it];
    topo.trk_nclust[it] = event.nhtrack[it];
    topo.trk_chisqr[it] = event.chisqr[it];
    topo.trk_dEdx[it] = event.dEdx[it];
    topo.trk_gffit[it] = event.GFfitstatus[it];
    topo.trk_inside[it] = event.GFinside[it];
    topo.trk_nsigma_proton[it] = event.nsigma_proton[it];
    topo.trk_nsigma_kaon[it] = event.nsigma_kaon[it];
    topo.trk_nsigma_pion[it] = event.nsigma_pion[it];
    topo.trk_nsigma_electron[it] = event.nsigma_electron[it];
    topo.trk_nsigma_deutron[it] = event.nsigma_deutron[it];
    topo.trk_nsigma_triton[it] = event.nsigma_triton[it];

    if(event.GFfromVtx[it]==1){
      topo.trk_htofReached[it] = 1;
      topo.trk_m2[it] = event.GFm2[it];
      topo.trk_invbeta[it] = event.GFinvbeta[it];
      topo.trk_nsigmaHtof_proton[it] = event.nsigma_protonHtof[it];
      topo.trk_nsigmaHtof_kaon[it] = event.nsigma_kaonHtof[it];
      topo.trk_nsigmaHtof_pion[it] = event.nsigma_pionHtof[it];
      topo.trk_nsigmaHtof_deutron[it] = event.nsigma_deutronHtof[it];
      topo.trk_nsigmaHtof_triton[it] = event.nsigma_tritonHtof[it];
    }

    // Momentum at primary reference vertex; dual DCA to kk and multi vertices.
    TVector3 tmom(qnan, qnan, qnan);
    if(event.GFfitstatus[it]){
      TVector3 tpos; Double_t tlen; Double_t ttof;
      if(topo.kkVtxValid
	 && GFtrackCont.ExtrapolateToPoint(it, kkRef, tpos, tmom, tlen, ttof)){
	topo.trk_dcaKK[it] = (tpos - kkRef).Mag();
      }
      if(topo.multiVtxValid){
	TVector3 tposM; TVector3 tmomM;
	if(GFtrackCont.ExtrapolateToPoint(it, multiRef, tposM, tmomM, tlen, ttof)){
	  topo.trk_dcaMulti[it] = (tposM - multiRef).Mag();
	}
      }
      const TVector3 dcaRef = topo.kkVtxValid ? kkRef
	: (topo.multiVtxValid ? multiRef : prodVertex);
      tmom = TVector3(qnan, qnan, qnan);
      if(GFtrackCont.ExtrapolateToPoint(it, dcaRef, tpos, tmom, tlen, ttof)){
	topo.trk_dcaPrimary[it] = topo.kkVtxValid ? topo.trk_dcaKK[it] : topo.trk_dcaMulti[it];
      } else if(!event.GFmom[it].empty()){
	tmom = TVector3(event.GFmom_x[it][0], event.GFmom_y[it][0], event.GFmom_z[it][0]);
      }
    }
    if(TMath::IsNaN(tmom.x())) continue;

    topo.trk_mom[it] = tmom.Mag();
    topo.trk_mom_x[it] = tmom.x();
    topo.trk_mom_y[it] = tmom.y();
    topo.trk_mom_z[it] = tmom.z();
    topo.trk_theta[it] = tmom.Theta()*TMath::RadToDeg();
    topo.trk_phi[it] = tmom.Phi()*TMath::RadToDeg();
  }

  // Particle class before V0 / role assignment (needed for pairing).
  for(Int_t it=0; it<ntTpc; ++it){
    const Int_t pidbit = event.pid[it];
    const Int_t q = event.charge[it];
    Int_t cls = topodef::kPidUnknown;
    if(event.isElectron[it]==1){
      cls = topodef::kPidElectron;
    } else if(q==-1 && event.nsigma_kaon[it]>mindEdxSigKaon
	      && event.nsigma_kaon[it]<maxdEdxSigKaon){
      cls = topodef::kPidKMinus;
    } else if((pidbit&4)==4 && (pidbit&1)!=1 && q==1){
      cls = topodef::kPidProton;
    } else if((pidbit&4)!=4 && (pidbit&1)==1 && q==1){
      cls = topodef::kPidPiPlus;
    } else if(((pidbit&4)==4 || (pidbit&1)==1) && q==1){
      cls = topodef::kPidProtonOrPiPlus;
    } else if((pidbit&1)==1 && q==-1){
      cls = topodef::kPidPiMinus;
    }
    // d/t appended after base classification; not in exclusiveVisibleKey v1.
    const Bool_t clearProtonOnly = ((pidbit&4)==4 && (pidbit&1)!=1 && q==1);
    if(q > 0 && !clearProtonOnly){
      if(!TMath::IsNaN(event.nsigma_deutron[it])
	 && TMath::Abs(event.nsigma_deutron[it]) < 3.){
	cls = topodef::kPidDeuteron;
      } else if(!TMath::IsNaN(event.nsigma_triton[it])
		&& TMath::Abs(event.nsigma_triton[it]) < 3.){
	cls = topodef::kPidTriton;
      } else if(!TMath::IsNaN(event.nsigma_deutron[it])
		&& !TMath::IsNaN(event.nsigma_triton[it])
		&& TMath::Abs(event.nsigma_deutron[it]) < 5.
		&& TMath::Abs(event.nsigma_triton[it]) < 5.){
	cls = topodef::kPidHeavyAmbiguous;
      }
    }
    topo.trk_pidClass[it] = cls;
  }

  // --- In-producer V0 reconstruction (mass-blind quality ranking) ---
  std::vector<kptopo::V0Cand> lamCands, k0Cands;
  {
    for(Int_t ip=0; ip<ntTpc; ++ip){
      if(event.isK18[ip] || event.isKurama[ip] || event.isBeam[ip] || event.isAccidental[ip]) continue;
      if(!isProtonLike(topo.trk_pidClass[ip]) || event.charge[ip]!=1) continue;
      if(!event.GFfitstatus[ip]) continue;
      Int_t repid_p = repidFromPid(event.pid[ip], event.isKurama[ip]);
      if(!GFtrackCont.TrackCheck(ip, repid_p)) continue;
      for(Int_t im=0; im<ntTpc; ++im){
	if(ip==im) continue;
	if(event.isK18[im] || event.isKurama[im] || event.isBeam[im] || event.isAccidental[im]) continue;
	if(topo.trk_pidClass[im]!=topodef::kPidPiMinus || event.charge[im]!=-1) continue;
	if(!event.GFfitstatus[im]) continue;
	Int_t repid_pi = 0;
	if(!GFtrackCont.TrackCheck(im, repid_pi)) continue;

	// Helix pre-filter (mass-blind): cheap pair screening before FindVertex.
	Double_t p_par[5] = {event.helix_cx[ip], event.helix_cy[ip], event.helix_z0[ip],
			     event.helix_r[ip], event.helix_dz[ip]};
	Double_t pi_par[5] = {event.helix_cx[im], event.helix_cy[im], event.helix_z0[im],
			      event.helix_r[im], event.helix_dz[im]};
	Int_t p_nh = event.helix_t[ip].size();
	Int_t pi_nh = event.helix_t[im].size();
	if(p_nh<=0 || pi_nh<=0) continue;
	Double_t p_tmin = event.helix_t[ip][0] - vtx_scan_range/p_par[3];
	Double_t p_tmax = TMath::Min(event.helix_t[ip][0]+vtx_scan_rangeInsideL/p_par[3],
				     event.helix_t[ip][p_nh-1]);
	Double_t pi_tmin = TMath::Max(event.helix_t[im][0]-vtx_scan_rangeInsideL/pi_par[3],
				    event.helix_t[im][pi_nh-1]);
	Double_t pi_tmax = event.helix_t[im][0] + vtx_scan_range/pi_par[3];
	TVector3 h_p_mom, h_pi_mom, h_l_mom;
	Double_t h_ppi_dist = qnan;
	TVector3 h_vertex = Kinematics::LambdaVertex(dMagneticField, p_par, pi_par,
						     p_tmin, p_tmax, pi_tmin, pi_tmax,
						     h_p_mom, h_pi_mom, h_l_mom, h_ppi_dist);
	if(TMath::IsNaN(h_ppi_dist) || h_ppi_dist > ppi_distcut) continue;
	if(TMath::Abs(h_vertex.x())>250. || TMath::Abs(h_vertex.y())>250.
	   || TMath::Abs(h_vertex.z())>250.) continue;
	TVector3 p_start(event.calpos_x[ip][0], event.calpos_y[ip][0], event.calpos_z[ip][0]);
	TVector3 p_end(event.calpos_x[ip][p_nh-1], event.calpos_y[ip][p_nh-1], event.calpos_z[ip][p_nh-1]);
	TVector3 pi_start(event.calpos_x[im][0], event.calpos_y[im][0], event.calpos_z[im][0]);
	TVector3 pi_end(event.calpos_x[im][pi_nh-1], event.calpos_y[im][pi_nh-1], event.calpos_z[im][pi_nh-1]);
	Double_t p_vd = qnan, pi_vd = qnan;
	if(!Kinematics::HelixDirection(h_vertex, p_start, p_end, p_vd)) continue;
	if(!Kinematics::HelixDirection(h_vertex, pi_start, pi_end, pi_vd)) continue;
	if(p_vd > p_vtx_distcut || pi_vd > pi_vtx_distcut) continue;

	TVector3 p_mom = h_p_mom, pi_mom = h_pi_mom;
	Double_t ppi_dist = h_ppi_dist;
	TVector3 l_vertex = h_vertex;
	Double_t p_extrap = qnan, pi_extrap = qnan;
	if(GFtrackCont.FindVertex(ip, im, repid_p, repid_pi,
				  p_extrap, pi_extrap,
				  p_mom, pi_mom, ppi_dist, l_vertex, vtx_scan_range)
	   && ppi_dist < GFppi_distcut){
	  // use GenFit vertex
	} // else keep helix fallback (already assigned)
	TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
	TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
	const TLorentzVector Llam = Lp + Lpi;
	Double_t ltarget_dist = qnan;
	TVector3 l_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, l_vertex, (p_mom+pi_mom), ltarget_dist);
	const Double_t decayLen = (l_vertex - l_pos_tgt).Mag();
	const Double_t cosPoint = kptopo::CosPointing(primaryRef, l_vertex, p_mom+pi_mom);
	const Bool_t inFid = kptopo::InTargetFiducial(l_vertex, tpc::ZTarget);
	Double_t alpha = qnan, qT = qnan;
	kptopo::ArmenterosPodolanski(p_mom, pi_mom, alpha, qT);
	const Double_t chis = topo.trk_chisqr[ip] + topo.trk_chisqr[im];
	const Int_t nclMin = TMath::Min(topo.trk_nclust[ip], topo.trk_nclust[im]);
	kptopo::V0Cand c;
	c.idPos = ip; c.idNeg = im;
	c.mass = Llam.M();
	c.dauDca = ppi_dist;
	c.decayLen = decayLen;
	c.cosPoint = cosPoint;
	c.alpha = alpha; c.qT = qT;
	c.vtx_x = l_vertex.x(); c.vtx_y = l_vertex.y(); c.vtx_z = l_vertex.z();
	c.mom = (p_mom+pi_mom).Mag();
	c.mom_x = (p_mom+pi_mom).x(); c.mom_y = (p_mom+pi_mom).y(); c.mom_z = (p_mom+pi_mom).z();
	c.pPos = p_mom.Mag(); c.pNeg = pi_mom.Mag();
	c.fiducial = inFid;
	c.quality = kptopo::V0QualityScore(ppi_dist, cosPoint, inFid, chis, nclMin);
	c.region = 0;
	lamCands.push_back(c);
	HF1(23020, c.mass);
      }
    }
    for(Int_t ip=0; ip<ntTpc; ++ip){
      if(event.isK18[ip] || event.isKurama[ip] || event.isBeam[ip] || event.isAccidental[ip]) continue;
      if(topo.trk_pidClass[ip]!=topodef::kPidPiPlus || event.charge[ip]!=1) continue;
      if(!event.GFfitstatus[ip]) continue;
      Int_t repid_pip = 0;
      if(!GFtrackCont.TrackCheck(ip, repid_pip)) continue;
      for(Int_t im=0; im<ntTpc; ++im){
	if(ip==im) continue;
	if(event.isK18[im] || event.isKurama[im] || event.isBeam[im] || event.isAccidental[im]) continue;
	if(topo.trk_pidClass[im]!=topodef::kPidPiMinus || event.charge[im]!=-1) continue;
	if(!event.GFfitstatus[im]) continue;
	Int_t repid_pim = 0;
	if(!GFtrackCont.TrackCheck(im, repid_pim)) continue;

	Double_t pip_par[5] = {event.helix_cx[ip], event.helix_cy[ip], event.helix_z0[ip],
			       event.helix_r[ip], event.helix_dz[ip]};
	Double_t pim_par[5] = {event.helix_cx[im], event.helix_cy[im], event.helix_z0[im],
			       event.helix_r[im], event.helix_dz[im]};
	Int_t pip_nh = event.helix_t[ip].size();
	Int_t pim_nh = event.helix_t[im].size();
	if(pip_nh<=0 || pim_nh<=0) continue;
	Double_t pip_tmin = event.helix_t[ip][0] - vtx_scan_range/pip_par[3];
	Double_t pip_tmax = event.helix_t[ip][0] + vtx_scan_range/pip_par[3];
	Double_t pim_tmin = TMath::Max(event.helix_t[im][0]-vtx_scan_rangeInsidePi/pim_par[3],
				       event.helix_t[im][pim_nh-1]);
	Double_t pim_tmax = event.helix_t[im][0] + vtx_scan_range/pim_par[3];
	TVector3 h_pip_mom, h_pim_mom, h_k0_mom;
	Double_t h_pipi_dist = qnan;
	TVector3 h_k0_vertex = Kinematics::LambdaVertex(dMagneticField, pip_par, pim_par,
							pip_tmin, pip_tmax, pim_tmin, pim_tmax,
							h_pip_mom, h_pim_mom, h_k0_mom, h_pipi_dist);
	if(TMath::IsNaN(h_pipi_dist) || h_pipi_dist > pipi_distcut) continue;
	if(TMath::Abs(h_k0_vertex.x())>250. || TMath::Abs(h_k0_vertex.y())>250.
	   || TMath::Abs(h_k0_vertex.z())>250.) continue;

	TVector3 pip_mom = h_pip_mom, pim_mom = h_pim_mom;
	Double_t pipi_dist = h_pipi_dist;
	TVector3 k0_vertex = h_k0_vertex;
	Double_t pip_extrap = qnan, pim_extrap = qnan;
	if(GFtrackCont.FindVertex(ip, im, repid_pip, repid_pim,
				  pip_extrap, pim_extrap,
				  pip_mom, pim_mom, pipi_dist, k0_vertex, vtx_scan_range)
	   && pipi_dist < GFpipi_distcut){
	  // use GenFit vertex
	}
	TLorentzVector Lpip(pip_mom, TMath::Hypot(pip_mom.Mag(), PionMass));
	TLorentzVector Lpim(pim_mom, TMath::Hypot(pim_mom.Mag(), PionMass));
	const TLorentzVector Lk0 = Lpip + Lpim;
	Double_t k0target_dist = qnan;
	TVector3 k0_pos_tgt = Kinematics::CalcCloseDistLambda(tgtpos, k0_vertex, (pip_mom+pim_mom), k0target_dist);
	const Double_t decayLen = (k0_vertex - k0_pos_tgt).Mag();
	const Double_t cosPoint = kptopo::CosPointing(primaryRef, k0_vertex, pip_mom+pim_mom);
	const Bool_t inFid = kptopo::InTargetFiducial(k0_vertex, tpc::ZTarget);
	Double_t alpha = qnan, qT = qnan;
	kptopo::ArmenterosPodolanski(pip_mom, pim_mom, alpha, qT);
	const Double_t chis = topo.trk_chisqr[ip] + topo.trk_chisqr[im];
	const Int_t nclMin = TMath::Min(topo.trk_nclust[ip], topo.trk_nclust[im]);
	kptopo::V0Cand c;
	c.idPos = ip; c.idNeg = im;
	c.mass = Lk0.M();
	c.dauDca = pipi_dist;
	c.decayLen = decayLen;
	c.cosPoint = cosPoint;
	c.alpha = alpha; c.qT = qT;
	c.vtx_x = k0_vertex.x(); c.vtx_y = k0_vertex.y(); c.vtx_z = k0_vertex.z();
	c.mom = (pip_mom+pim_mom).Mag();
	c.mom_x = (pip_mom+pim_mom).x(); c.mom_y = (pip_mom+pim_mom).y(); c.mom_z = (pip_mom+pim_mom).z();
	c.pPos = pip_mom.Mag(); c.pNeg = pim_mom.Mag();
	c.fiducial = inFid;
	c.quality = kptopo::V0QualityScore(pipi_dist, cosPoint, inFid, chis, nclMin);
	c.region = 0;
	k0Cands.push_back(c);
	HF1(21003, c.mass);
      }
    }
    std::sort(lamCands.begin(), lamCands.end(),
	      [](const kptopo::V0Cand& a, const kptopo::V0Cand& b){
		return a.quality < b.quality; });
    std::sort(k0Cands.begin(), k0Cands.end(),
	      [](const kptopo::V0Cand& a, const kptopo::V0Cand& b){
		return a.quality < b.quality; });
    for(size_t i=0; i<lamCands.size(); ++i){
      lamCands[i].region = kptopo::MassRegionLambda(
	lamCands[i].mass, LambdaMass,
	topodef::kLambdaMassWindow, topodef::kLambdaSidebandLow, topodef::kLambdaSidebandHigh);
    }
    for(size_t i=0; i<k0Cands.size(); ++i){
      k0Cands[i].region = kptopo::MassRegionK0S(
	k0Cands[i].mass, K0Mass,
	topodef::kK0SMassWindow, topodef::kK0SSidebandLow, topodef::kK0SSidebandHigh);
    }
    auto fillV0Vec = [](const std::vector<kptopo::V0Cand>& src, Topo& t,
			Int_t& nCand,
			std::vector<Double_t>& mass,
			std::vector<Int_t>& idPos, std::vector<Int_t>& idNeg,
			std::vector<Double_t>& vx, std::vector<Double_t>& vy, std::vector<Double_t>& vz,
			std::vector<Double_t>& dauDca, std::vector<Double_t>& decayLen,
			std::vector<Double_t>& cosPoint,
			std::vector<Double_t>& mom, std::vector<Double_t>& mx,
			std::vector<Double_t>& my, std::vector<Double_t>& mz,
			std::vector<Double_t>& alpha, std::vector<Double_t>& qT,
			std::vector<Double_t>& quality, std::vector<Int_t>& region,
			std::vector<Int_t>& fiducial){
      nCand = (Int_t)src.size();
      for(size_t i=0; i<src.size(); ++i){
	mass.push_back(src[i].mass);
	idPos.push_back(src[i].idPos); idNeg.push_back(src[i].idNeg);
	vx.push_back(src[i].vtx_x); vy.push_back(src[i].vtx_y); vz.push_back(src[i].vtx_z);
	dauDca.push_back(src[i].dauDca); decayLen.push_back(src[i].decayLen);
	cosPoint.push_back(src[i].cosPoint);
	mom.push_back(src[i].mom); mx.push_back(src[i].mom_x);
	my.push_back(src[i].mom_y); mz.push_back(src[i].mom_z);
	alpha.push_back(src[i].alpha); qT.push_back(src[i].qT);
	quality.push_back(src[i].quality); region.push_back(src[i].region);
	fiducial.push_back(src[i].fiducial ? 1 : 0);
      }
    };
    fillV0Vec(lamCands, topo, topo.nLamCand,
	      topo.lamCand_mass, topo.lamCand_idPos, topo.lamCand_idNeg,
	      topo.lamCand_vtx_x, topo.lamCand_vtx_y, topo.lamCand_vtx_z,
	      topo.lamCand_dauDca, topo.lamCand_decayLen, topo.lamCand_cosPoint,
	      topo.lamCand_mom, topo.lamCand_mom_x, topo.lamCand_mom_y, topo.lamCand_mom_z,
	      topo.lamCand_alpha, topo.lamCand_qT, topo.lamCand_quality,
	      topo.lamCand_region, topo.lamCand_fiducial);
    fillV0Vec(k0Cands, topo, topo.nK0Cand,
	      topo.k0Cand_mass, topo.k0Cand_idPos, topo.k0Cand_idNeg,
	      topo.k0Cand_vtx_x, topo.k0Cand_vtx_y, topo.k0Cand_vtx_z,
	      topo.k0Cand_dauDca, topo.k0Cand_decayLen, topo.k0Cand_cosPoint,
	      topo.k0Cand_mom, topo.k0Cand_mom_x, topo.k0Cand_mom_y, topo.k0Cand_mom_z,
	      topo.k0Cand_alpha, topo.k0Cand_qT, topo.k0Cand_quality,
	      topo.k0Cand_region, topo.k0Cand_fiducial);

    Int_t lamSel = -1, k0Sel = -1;
    for(size_t i=0; i<lamCands.size(); ++i) if(lamCands[i].region==1){ lamSel=(Int_t)i; break; }
    if(lamSel<0 && !lamCands.empty()) lamSel = 0;
    for(size_t i=0; i<k0Cands.size(); ++i) if(k0Cands[i].region==1){ k0Sel=(Int_t)i; break; }
    if(k0Sel<0 && !k0Cands.empty()) k0Sel = 0;
    topo.lamBestId = lamSel;
    topo.k0BestId = k0Sel;

    // Overlap resolution: prefer higher-quality (lower score) V0.
    // Do not adopt both overlapping candidates in the exclusive classification.
    Bool_t lamDemotedByOverlap = false;
    if(lamSel>=0 && k0Sel>=0){
      const Int_t lt[2] = {lamCands[lamSel].idPos, lamCands[lamSel].idNeg};
      const Int_t kt[2] = {k0Cands[k0Sel].idPos, k0Cands[k0Sel].idNeg};
      for(Int_t a=0;a<2;++a) for(Int_t b=0;b<2;++b){
	  if(lt[a]==kt[b] && lt[a]>=0){
	    topo.v0OverlapFlag = 1;
	    topo.v0OverlapTrackIds.push_back(lt[a]);
	  }
	}
      if(topo.v0OverlapFlag){
	if(lamCands[lamSel].quality <= k0Cands[k0Sel].quality){
	  k0Sel = -1; topo.k0BestId = -1;
	} else {
	  lamSel = -1; topo.lamBestId = -1;
	  lamDemotedByOverlap = true;
	}
      }
    }

    // Nominal Lambda scalar branches: in-producer preferred over input.
    // If overlap demoted the in-producer Lambda, do NOT fall back to DstE42
    // input (it may share the same contested tracks).
    if(lamSel >= 0){
      const kptopo::V0Cand& c = lamCands[lamSel];
      topo.lFlag = 1;
      topo.lMass = c.mass;
      topo.lMom = c.mom;
      topo.lMom_x = c.mom_x; topo.lMom_y = c.mom_y; topo.lMom_z = c.mom_z;
      topo.lDecayVtx_x = c.vtx_x; topo.lDecayVtx_y = c.vtx_y; topo.lDecayVtx_z = c.vtx_z;
      topo.lPPiDist = c.dauDca;
      topo.lDecayLen = c.decayLen;
      topo.lDaughterId = {c.idPos, c.idNeg};
      lDaughter = topo.lDaughterId;
      topo.lIsSignal = (c.region==1) ? 1 : 0;
      topo.lIsSideband = (c.region==2 || c.region==3) ? 1 : 0;
      topo.lambdaRegion = c.region;
      topo.nLambdaSignal = topo.lIsSignal;
      topo.nLambdaSideband = topo.lIsSideband;
      TLorentzVector LvL(TVector3(c.mom_x,c.mom_y,c.mom_z), TMath::Hypot(c.mom, LambdaMass));
      LvL.Boost(boostToX);
      topo.lPStar = LvL.Vect().Mag();
    } else if(inputLFlag && !lamDemotedByOverlap){
      topo.lFlag = 1;
      topo.lMass = inputLMass; topo.lMom = inputLMom;
      topo.lMom_x = inputLMom_x; topo.lMom_y = inputLMom_y; topo.lMom_z = inputLMom_z;
      topo.lDecayVtx_x = inputLDecayVtx_x; topo.lDecayVtx_y = inputLDecayVtx_y;
      topo.lDecayVtx_z = inputLDecayVtx_z;
      topo.lPPiDist = inputLPPiDist; topo.lPPiAngle = inputLPPiAngle;
      topo.lDaughterId = lDaughter;
      TVector3 ldecay(topo.lDecayVtx_x, topo.lDecayVtx_y, topo.lDecayVtx_z);
      topo.lDecayLen = (ldecay - primaryRef).Mag();
      Double_t dm = TMath::Abs(topo.lMass - LambdaMass);
      topo.lIsSignal = (dm < topodef::kLambdaMassWindow) ? 1 : 0;
      topo.lIsSideband = (dm >= topodef::kLambdaSidebandLow
			  && dm < topodef::kLambdaSidebandHigh) ? 1 : 0;
      topo.lambdaRegion = topo.lIsSignal ? 1 : (topo.lIsSideband ? 2 : 0);
      topo.nLambdaSignal = topo.lIsSignal;
      topo.nLambdaSideband = topo.lIsSideband;
      TLorentzVector LvL(TVector3(topo.lMom_x, topo.lMom_y, topo.lMom_z),
			 TMath::Hypot(topo.lMom, LambdaMass));
      LvL.Boost(boostToX);
      topo.lPStar = LvL.Vect().Mag();
    }

    if(k0Sel >= 0){
      const kptopo::V0Cand& c = k0Cands[k0Sel];
      topo.k0Flag = 1;
      // Exclusive nK0S / daughter roles only for signal-region candidate.
      // Sideband masses remain in k0Mass for SB subtraction / mass fits.
      topo.nK0S = (c.region==1) ? 1 : 0;
      topo.k0Mass = c.mass;
      topo.k0Mom = c.mom;
      topo.k0Mom_x = c.mom_x; topo.k0Mom_y = c.mom_y; topo.k0Mom_z = c.mom_z;
      topo.k0DaughterId = {c.idPos, c.idNeg};
    }
  }

  // LambdaPPiAngle: compute from daughters if input unset.
  if(topo.lFlag && TMath::IsNaN(topo.lPPiAngle) && topo.lDaughterId.size()>=2){
    const Int_t ip = topo.lDaughterId[0];
    const Int_t im = topo.lDaughterId[1];
    if(ip>=0 && ip<ntTpc && im>=0 && im<ntTpc
       && !TMath::IsNaN(topo.trk_mom[ip]) && !TMath::IsNaN(topo.trk_mom[im])){
      TVector3 pp(topo.trk_mom_x[ip], topo.trk_mom_y[ip], topo.trk_mom_z[ip]);
      TVector3 pm(topo.trk_mom_x[im], topo.trk_mom_y[im], topo.trk_mom_z[im]);
      if(pp.Mag()>0. && pm.Mag()>0.)
	topo.lPPiAngle = pp.Unit().Angle(pm.Unit());
      topo.lPPiAngleValid = 1;
    }
  } else if(topo.lFlag && !TMath::IsNaN(topo.lPPiAngle)){
    topo.lPPiAngleValid = 1;
  }

  // Lambda daughter flags (candidate / signal / sideband kept separate).
  for(Int_t it=0; it<ntTpc; ++it){
    if(!isInIdList(lDaughter, it)) continue;
    topo.trk_isLambdaCandidateDaughter[it] = 1;
    if(topo.lIsSignal) topo.trk_isLambdaSignalDaughter[it] = 1;
    if(topo.lIsSideband) topo.trk_isLambdaSidebandDaughter[it] = 1;
  }

  // K0S daughter flags: only signal-region K0S enters exclusive roles.
  if(topo.nK0S >= 1){
    for(size_t i=0; i<topo.k0DaughterId.size(); ++i){
      const Int_t id = topo.k0DaughterId[i];
      if(id>=0 && id<ntTpc) topo.trk_isK0Daughter[id] = 1;
    }
  }

  // Role assignment. Priority:
  // K18/Forward/Beam/Accidental > LambdaSignalDaughter > K0Daughter
  // > Outside > Prompt/Displaced > Unknown
  for(Int_t it=0; it<ntTpc; ++it){
    if(event.isK18[it]==1)             topo.trk_role[it] = topodef::kRoleK18;
    else if(event.isKurama[it]==1)     topo.trk_role[it] = topodef::kRoleForward;
    else if(event.isBeam[it]==1)       topo.trk_role[it] = topodef::kRoleBeam;
    else if(event.isAccidental[it]==1) topo.trk_role[it] = topodef::kRoleAccidental;
    else if(topo.trk_isLambdaSignalDaughter[it])
      topo.trk_role[it] = topodef::kRoleV0Daughter;
    else if(topo.trk_isK0Daughter[it]) topo.trk_role[it] = topodef::kRoleK0Daughter;
    else if(event.GFinside[it]!=1)     topo.trk_role[it] = topodef::kRoleOutside;
    else if(!TMath::IsNaN(topo.trk_dcaPrimary[it])
	    && topo.trk_dcaPrimary[it] < topodef::kNominalDca)
      topo.trk_role[it] = topodef::kRolePrompt;
    else                               topo.trk_role[it] = topodef::kRoleDisplaced;
  }

  // Rest-frame momentum of the missing system + q-frame projections.
  for(Int_t it=0; it<ntTpc; ++it){
    if(TMath::IsNaN(topo.trk_mom[it])) continue;
    const TVector3 pLab(topo.trk_mom_x[it], topo.trk_mom_y[it], topo.trk_mom_z[it]);
    kptopo::ProjectOnQ(pLab, qLab, topo.trk_pParallelQ[it],
		       topo.trk_pTransverseQ[it], topo.trk_cosThetaQ[it]);
    Double_t mass = PionMass;
    switch(topo.trk_pidClass[it]){
    case topodef::kPidProton:  mass = ProtonMass;   break;
    case topodef::kPidKMinus:  mass = KaonMass;     break;
    case topodef::kPidElectron:mass = ElectronMass; break;
    case topodef::kPidDeuteron:mass = DeuteronMass; break;
    case topodef::kPidTriton:  mass = TritonMass;   break;
    default: mass = PionMass; break;
    }
    TLorentzVector Lv(pLab, TMath::Hypot(topo.trk_mom[it], mass));
    Lv.Boost(boostToX);
    if(topo.trk_pidClass[it]==topodef::kPidProtonOrPiPlus
       || topo.trk_pidClass[it]==topodef::kPidHeavyAmbiguous){
      TLorentzVector LvP(pLab, TMath::Hypot(topo.trk_mom[it], ProtonMass));
      TLorentzVector LvPi(pLab, TMath::Hypot(topo.trk_mom[it], PionMass));
      LvP.Boost(boostToX); LvPi.Boost(boostToX);
      topo.trk_pstarProtonHyp[it] = LvP.Vect().Mag();
      topo.trk_pstarPionHyp[it] = LvPi.Vect().Mag();
      topo.trk_pstar[it] = topo.trk_pstarPionHyp[it]; // legacy branch; use hyps for physics
    } else {
      topo.trk_pstar[it] = Lv.Vect().Mag();
      topo.trk_costStar[it] = Lv.Vect().CosTheta();
    }
  }

  auto bumpPid = [](Int_t pid, Int_t& nP, Int_t& nPip, Int_t& nPim,
		    Int_t& nKm, Int_t& nE, Int_t& nA, Int_t& nU){
    switch(pid){
    case topodef::kPidProton:         ++nP;   break;
    case topodef::kPidPiPlus:         ++nPip; break;
    case topodef::kPidPiMinus:        ++nPim; break;
    case topodef::kPidKMinus:         ++nKm;  break;
    case topodef::kPidElectron:       ++nE;   break;
    case topodef::kPidProtonOrPiPlus: ++nA;   break;
    // d/t/heavy-amb map to unknown for exclusiveVisibleKey v1 packing.
    default:                          ++nU;   break;
    }
  };
  auto bumpDt = [](Int_t pid, Int_t& nD, Int_t& nT, Int_t& nH){
    switch(pid){
    case topodef::kPidDeuteron:       ++nD; break;
    case topodef::kPidTriton:         ++nT; break;
    case topodef::kPidHeavyAmbiguous: ++nH; break;
    default: break;
    }
  };

  // Multiplicity ladder + visible-exclusive key + semi-exclusive tags.
  for(Int_t ith=0; ith<topodef::kNThr; ++ith){
    Int_t nraw=0, nprompt=0, np=0, npip=0, npim=0, namb=0, ndisp=0;
    Int_t nOut=0, nBelow=0, nNoMom=0;

    // Nominal exclusive basis: exclude Lambda SIGNAL daughters and K0S daughters.
    Int_t pP=0,pPip=0,pPim=0,pKm=0,pE=0,pA=0,pU=0;
    Int_t dP=0,dPip=0,dPim=0,dKm=0,dE=0,dA=0,dU=0;
    Int_t nExclPrompt=0, nExclDisp=0;
    // Sideband hypothesis: exclude Lambda SIDEBAND daughters (+ K0S daughters).
    Int_t hpP=0,hpPip=0,hpPim=0,hpKm=0,hpE=0,hpA=0,hpU=0;
    Int_t hdP=0,hdPip=0,hdPim=0,hdKm=0,hdE=0,hdA=0,hdU=0;
    // Residual for semi-tags with Lambda SIGNAL daughters removed (no K0 exclusion
    // unless the tag itself is K0-based; residual def = drop K18/Fwd/Lam daughters).
    Int_t rpP=0,rpPip=0,rpPim=0,rpKm=0,rpE=0,rpA=0,rpU=0;
    Int_t rdP=0,rdPip=0,rdPim=0,rdKm=0,rdE=0,rdA=0,rdU=0;
    // Residual with Lambda SIDEBAND daughters removed (symmetric sideband tags).
    Int_t spP=0,spPip=0,spPim=0,spKm=0,spE=0,spA=0,spU=0;
    Int_t sdP=0,sdPip=0,sdPim=0,sdKm=0,sdE=0,sdA=0,sdU=0;
    // Residual with K0S daughters removed (for K0SProton).
    Int_t kpP=0,kpPip=0,kpPim=0,kpKm=0,kpE=0,kpA=0,kpU=0;
    Int_t nPD=0,nPT=0,nPH=0, nDD=0,nDT=0,nDH=0;

    for(Int_t it=0; it<ntTpc; ++it){
      const Int_t role = topo.trk_role[it];
      if(role==topodef::kRoleK18 || role==topodef::kRoleForward
	 || role==topodef::kRoleBeam || role==topodef::kRoleAccidental) continue;

      if(TMath::IsNaN(topo.trk_mom[it])){ ++nNoMom; continue; }
      if(topo.trk_mom[it] < topodef::kThrPMin[ith]
	 || topo.trk_nclust[it] < topodef::kThrNClust[ith]){
	++nBelow; continue;
      }

      nraw++;
      if(role==topodef::kRoleOutside){ ++nOut; continue; }

      const Bool_t promptLike =
	!TMath::IsNaN(topo.trk_dcaPrimary[it])
	&& topo.trk_dcaPrimary[it] < topodef::kThrDca[ith];
      const Int_t pid = topo.trk_pidClass[it];
      const Bool_t isLamSig = topo.trk_isLambdaSignalDaughter[it];
      const Bool_t isLamSb  = topo.trk_isLambdaSidebandDaughter[it];
      const Bool_t isK0     = topo.trk_isK0Daughter[it];

      // Legacy coarse counts: exclude V0/K0 daughter roles from prompt.
      if(role==topodef::kRoleV0Daughter || role==topodef::kRoleK0Daughter){
	// still contribute to raw only
      } else if(!promptLike){
	ndisp++;
      } else {
	nprompt++;
	switch(pid){
	case topodef::kPidProton:         np++;   break;
	case topodef::kPidPiPlus:         npip++; break;
	case topodef::kPidPiMinus:        npim++; break;
	case topodef::kPidProtonOrPiPlus: namb++; break;
	default: break;
	}
      }

      // Visible-exclusive PID (exclude Lam-signal + K0 daughters).
      if(!(isLamSig || isK0)){
	if(promptLike){ ++nExclPrompt; bumpPid(pid,pP,pPip,pPim,pKm,pE,pA,pU); bumpDt(pid,nPD,nPT,nPH); }
	else           { ++nExclDisp;   bumpPid(pid,dP,dPip,dPim,dKm,dE,dA,dU); bumpDt(pid,nDD,nDT,nDH); }
      }
      // Lambda-hypothesis key (exclude Lam-sideband + K0 when sideband region).
      const Bool_t exclHyp = topo.lIsSideband ? (isLamSb || isK0)
	: (topo.lIsSignal ? (isLamSig || isK0) : isK0);
      if(!exclHyp){
	if(promptLike) bumpPid(pid,hpP,hpPip,hpPim,hpKm,hpE,hpA,hpU);
	else           bumpPid(pid,hdP,hdPip,hdPim,hdKm,hdE,hdA,hdU);
      }
      // Residual for Lambda signal tags: drop Lam-signal and K0 daughters.
      if(!isLamSig && !isK0){
	if(promptLike) bumpPid(pid,rpP,rpPip,rpPim,rpKm,rpE,rpA,rpU);
	else           bumpPid(pid,rdP,rdPip,rdPim,rdKm,rdE,rdA,rdU);
      }
      // Residual for Lambda sideband tags: drop Lam-sideband daughters only.
      if(!isLamSb){
	if(promptLike) bumpPid(pid,spP,spPip,spPim,spKm,spE,spA,spU);
	else           bumpPid(pid,sdP,sdPip,sdPim,sdKm,sdE,sdA,sdU);
      }
      // Residual for K0SProton: drop K0 daughters (and not Lam-signal by role).
      if(!isK0 && !isLamSig){
	if(promptLike) bumpPid(pid,kpP,kpPip,kpPim,kpKm,kpE,kpA,kpU);
      }
    }

    topo.nchRaw[ith] = nraw;
    topo.nchPrompt[ith] = nprompt;
    topo.nProton[ith] = np;
    topo.nPiPlus[ith] = npip;
    topo.nPiMinus[ith] = npim;
    topo.nAmbiguous[ith] = namb;
    topo.nDisplaced[ith] = ndisp;
    topo.nOutsideUsed[ith] = nOut;
    topo.nBelowThreshold[ith] = nBelow;
    topo.nNoMomentum[ith] = nNoMom;

    topo.nPromptProton[ith] = pP;
    topo.nPromptPiPlus[ith] = pPip;
    topo.nPromptPiMinus[ith] = pPim;
    topo.nPromptKMinus[ith] = pKm;
    topo.nPromptElectron[ith] = pE;
    topo.nPromptAmbiguous[ith] = pA;
    topo.nPromptUnknown[ith] = pU;
    topo.nDisplacedProton[ith] = dP;
    topo.nDisplacedPiPlus[ith] = dPip;
    topo.nDisplacedPiMinus[ith] = dPim;
    topo.nDisplacedKMinus[ith] = dKm;
    topo.nDisplacedElectron[ith] = dE;
    topo.nDisplacedAmbiguous[ith] = dA;
    topo.nDisplacedUnknown[ith] = dU;
    topo.nPromptDeuteron[ith] = nPD;
    topo.nPromptTriton[ith] = nPT;
    topo.nPromptHeavyAmbiguous[ith] = nPH;
    topo.nDisplacedDeuteron[ith] = nDD;
    topo.nDisplacedTriton[ith] = nDT;
    topo.nDisplacedHeavyAmbiguous[ith] = nDH;

    // Event-by-event PID closure: sum of 7 PID bins must equal exclusive track count.
    topo.nPromptClosureDiff[ith] =
      (pP+pPip+pPim+pKm+pE+pA+pU) - nExclPrompt;
    topo.nDisplacedClosureDiff[ith] =
      (dP+dPip+dPim+dKm+dE+dA+dU) - nExclDisp;

    Int_t cls = topodef::kTopo3PlusPrompt;
    if(nprompt==0)      cls = topodef::kTopoEmpty;
    else if(nprompt==1){
      if(np==1)         cls = topodef::kTopo1Proton;
      else if(npip==1)  cls = topodef::kTopo1PiPlus;
      else if(npim==1)  cls = topodef::kTopo1PiMinus;
      else              cls = topodef::kTopo1Other;
    }
    else if(nprompt==2) cls = topodef::kTopo2Prompt;
    topo.topoClass[ith] = cls;

    const Int_t nLamKey = topo.nLambdaSignal;
    const Int_t nK0Key  = topo.nK0S;
    topo.exclusiveVisibleKey[ith] = topodef::EncodeExclusiveVisible(
      nLamKey, nK0Key,
      pP,pPip,pPim,pKm,pE,pA,pU,
      dP,dPip,dPim,dKm,dE,dA,dU);

    const Int_t nLamHyp = topo.lIsSideband ? 1 : (topo.lIsSignal ? 1 : 0);
    topo.exclusiveVisibleKeyLambdaHypothesis[ith] = topodef::EncodeExclusiveVisible(
      nLamHyp, nK0Key,
      hpP,hpPip,hpPim,hpKm,hpE,hpA,hpU,
      hdP,hdPip,hdPim,hdKm,hdE,hdA,hdU);

    // Semi-exclusive tags: independent ifs (never else-if).
    ULong64_t bits = 0;
    if(topo.lIsSignal){
      bits |= topodef::SemiTagMask(topodef::kTagLambda);
      if(rpPim >= 1) bits |= topodef::SemiTagMask(topodef::kTagLambdaPiMinus); // 1NA-enriched
      if(rpPip >= 1) bits |= topodef::SemiTagMask(topodef::kTagLambdaPiPlus);  // control tag
      if(rpP >= 1)   bits |= topodef::SemiTagMask(topodef::kTagLambdaProton);  // 2NA-enriched
      if(rpP >= 1 && rpPip==0 && rpPim==0 && rpA==0)
	bits |= topodef::SemiTagMask(topodef::kTagLambdaProtonNoPion);
    }
    if(rpPip >= 2) bits |= topodef::SemiTagMask(topodef::kTagTwoPiPlus);
    if(rpPim >= 1 && rdPip >= 1)
      bits |= topodef::SemiTagMask(topodef::kTagSigmaPlusCharged);
    if(rpPip >= 1 && rdPim >= 1)
      bits |= topodef::SemiTagMask(topodef::kTagSigmaMinusCharged);
    if(topo.escFlag)
      bits |= topodef::SemiTagMask(topodef::kTagKminusEscape);
    if(topo.nK0S >= 1){
      bits |= topodef::SemiTagMask(topodef::kTagK0S);
      if(kpP >= 1) bits |= topodef::SemiTagMask(topodef::kTagK0SProton);
    }
    if(topo.lIsSideband){
      bits |= topodef::SemiTagMask(topodef::kTagLambdaSideband);
      if(spPim >= 1) bits |= topodef::SemiTagMask(topodef::kTagLambdaSidebandPiMinus);
      if(spPip >= 1) bits |= topodef::SemiTagMask(topodef::kTagLambdaSidebandPiPlus);
      if(spP >= 1)   bits |= topodef::SemiTagMask(topodef::kTagLambdaSidebandProton);
      if(spP >= 1 && spPip==0 && spPim==0 && spA==0)
	bits |= topodef::SemiTagMask(topodef::kTagLambdaSidebandProtonNoPion);
    }
    if(topo.lIsSignal){
      Int_t nResD = 0, nResT = 0;
      for(Int_t it=0; it<ntTpc; ++it){
	if(topo.trk_isLambdaSignalDaughter[it] || topo.trk_isK0Daughter[it]) continue;
	if(topo.trk_role[it]!=topodef::kRolePrompt && topo.trk_role[it]!=topodef::kRoleDisplaced) continue;
	if(topo.trk_pidClass[it]==topodef::kPidDeuteron) ++nResD;
	if(topo.trk_pidClass[it]==topodef::kPidTriton) ++nResT;
      }
      if(nResD >= 1) bits |= topodef::SemiTagMask(topodef::kTagLambdaDeuteron);
      if(nResT >= 1) bits |= topodef::SemiTagMask(topodef::kTagLambdaTriton);
    }
    if(topo.v0OverlapFlag)
      bits |= topodef::SemiTagMask(topodef::kTagV0Ambiguous);
    topo.semiTagBits[ith] = bits;
  }

  // Lambda + additional hadrons: all prompt pairs; scalars = best-by-max-p.
  if(topo.lFlag==1){
    Double_t bestmom = -1., bestpimmom = -1., bestpipmom = -1.;
    TVector3 lmom3(topo.lMom_x, topo.lMom_y, topo.lMom_z);
    TLorentzVector LvL(lmom3, TMath::Hypot(topo.lMom, LambdaMass));
    auto recoilM = [&](const TLorentzVector& Ladd)->Double_t{
      return (LvMissX - (LvL + Ladd)).M();
    };
    for(Int_t it=0; it<ntTpc; ++it){
      if(topo.trk_role[it]!=topodef::kRolePrompt) continue;
      if(topo.trk_isLambdaCandidateDaughter[it]) continue;
      if(topo.trk_isK0Daughter[it]) continue;
      if(TMath::IsNaN(topo.trk_mom[it])) continue;
      if(topo.trk_mom[it] < topodef::kThrPMin[0]) continue;
      TVector3 hmom(topo.trk_mom_x[it], topo.trk_mom_y[it], topo.trk_mom_z[it]);
      if(topo.trk_pidClass[it]==topodef::kPidProton){
	topo.nAddProton++;
	TLorentzVector LvP(hmom, TMath::Hypot(hmom.Mag(), ProtonMass));
	TLorentzVector LvLP = LvL + LvP;
	topo.addProtonIdAll.push_back(it);
	topo.mLambdaPAll.push_back(LvLP.M());
	topo.pLambdaPAll.push_back((lmom3+hmom).Mag());
	topo.cosThetaLambdaPAll.push_back((lmom3.Mag()>0.&&hmom.Mag()>0.)
					  ? lmom3.Unit().Dot(hmom.Unit()) : qnan);
	topo.mRecoilLambdaPAll.push_back(recoilM(LvP));
	Double_t cThQ = qnan;
	kptopo::ProjectOnQ(hmom, qLab, topo.trk_pParallelQ[it], topo.trk_pTransverseQ[it], cThQ);
	topo.cosThetaLambdaPQAll.push_back(cThQ);
	if(topo.trk_mom[it] > bestmom){
	  bestmom = topo.trk_mom[it]; topo.addProtonId = it;
	}
      }
      if(topo.trk_pidClass[it]==topodef::kPidPiMinus){
	topo.nAddPim++;
	TLorentzVector LvPim(hmom, TMath::Hypot(hmom.Mag(), PionMass));
	TLorentzVector LvLPim = LvL + LvPim;
	topo.addPimIdAll.push_back(it);
	topo.mLambdaPimAll.push_back(LvLPim.M());
	topo.pLambdaPimAll.push_back((lmom3+hmom).Mag());
	topo.cosThetaLambdaPimAll.push_back((lmom3.Mag()>0.&&hmom.Mag()>0.)
					    ? lmom3.Unit().Dot(hmom.Unit()) : qnan);
	topo.mRecoilLambdaPimAll.push_back(recoilM(LvPim));
	Double_t cThQ = qnan;
	kptopo::ProjectOnQ(hmom, qLab, topo.trk_pParallelQ[it], topo.trk_pTransverseQ[it], cThQ);
	topo.cosThetaLambdaPimQAll.push_back(cThQ);
	if(topo.trk_mom[it] > bestpimmom){
	  bestpimmom = topo.trk_mom[it]; topo.addPimId = it;
	}
      }
      if(topo.trk_pidClass[it]==topodef::kPidPiPlus){
	TLorentzVector LvPip(hmom, TMath::Hypot(hmom.Mag(), PionMass));
	TLorentzVector LvLPip = LvL + LvPip;
	topo.addPipIdAll.push_back(it);
	topo.mLambdaPipAll.push_back(LvLPip.M());
	topo.pLambdaPipAll.push_back((lmom3+hmom).Mag());
	topo.cosThetaLambdaPipAll.push_back((lmom3.Mag()>0.&&hmom.Mag()>0.)
					     ? lmom3.Unit().Dot(hmom.Unit()) : qnan);
	if(topo.trk_mom[it] > bestpipmom) bestpipmom = topo.trk_mom[it];
      }
    }
    if(topo.addProtonId >= 0){
      const Int_t ip = topo.addProtonId;
      TVector3 pmom(topo.trk_mom_x[ip], topo.trk_mom_y[ip], topo.trk_mom_z[ip]);
      TLorentzVector LvP(pmom, TMath::Hypot(pmom.Mag(), ProtonMass));
      TLorentzVector LvLP = LvL + LvP;
      topo.addProtonMom = pmom.Mag();
      topo.mLambdaP = LvLP.M();
      topo.pLambdaP = (lmom3 + pmom).Mag();
      topo.cosThetaLambdaP = (lmom3.Mag()>0. && pmom.Mag()>0.)
	? lmom3.Unit().Dot(pmom.Unit()) : qnan;
      topo.mRecoilLambdaP = recoilM(LvP);
    }
    if(topo.addPimId >= 0){
      const Int_t ip = topo.addPimId;
      TVector3 pimom(topo.trk_mom_x[ip], topo.trk_mom_y[ip], topo.trk_mom_z[ip]);
      TLorentzVector LvPim(pimom, TMath::Hypot(pimom.Mag(), PionMass));
      TLorentzVector LvLPim = LvL + LvPim;
      topo.addPimMom = pimom.Mag();
      topo.mLambdaPim = LvLPim.M();
      topo.pLambdaPim = (lmom3 + pimom).Mag();
      topo.cosThetaLambdaPim = (lmom3.Mag()>0. && pimom.Mag()>0.)
	? lmom3.Unit().Dot(pimom.Unit()) : qnan;
      topo.mRecoilLambdaPim = recoilM(LvPim);
    }
    topo.mRecoilLambda = (LvMissX - LvL).M();
    topo.pRecoilLambda = (LvMissX - LvL).P();
    topo.eRecoilLambda = (LvMissX - LvL).E();
  }

  // Visible energy flow: sum detected tracks excluding selected V0 daughters.
  {
    TVector3 pSum(0,0,0);
    Double_t ePi = 0., eP = 0.;
    Double_t cThQVis = qnan;
    Int_t qSum = 0;
    TLorentzVector LvDet(0,0,0,0);
    for(Int_t it=0; it<ntTpc; ++it){
      if(topo.trk_isLambdaSignalDaughter[it] || topo.trk_isK0Daughter[it]) continue;
      if(topo.trk_role[it]==topodef::kRoleK18 || topo.trk_role[it]==topodef::kRoleForward
	 || topo.trk_role[it]==topodef::kRoleBeam || topo.trk_role[it]==topodef::kRoleAccidental) continue;
      if(TMath::IsNaN(topo.trk_mom[it])) continue;
      TVector3 p(topo.trk_mom_x[it], topo.trk_mom_y[it], topo.trk_mom_z[it]);
      pSum += p;
      qSum += topo.trk_charge[it];
      TLorentzVector LvPi(p, TMath::Hypot(p.Mag(), PionMass));
      TLorentzVector LvPr(p, TMath::Hypot(p.Mag(), ProtonMass));
      ePi += LvPi.E();
      eP += LvPr.E();
      Double_t m = PionMass;
      if(topo.trk_pidClass[it]==topodef::kPidProton) m = ProtonMass;
      else if(topo.trk_pidClass[it]==topodef::kPidKMinus) m = KaonMass;
      else if(topo.trk_pidClass[it]==topodef::kPidElectron) m = ElectronMass;
      LvDet += TLorentzVector(p, TMath::Hypot(p.Mag(), m));
    }
    topo.visibleNetCharge = qSum;
    topo.visibleScalarPSum = pSum.Mag();
    topo.visibleVectorPSum_x = pSum.x();
    topo.visibleVectorPSum_y = pSum.y();
    topo.visibleVectorPSum_z = pSum.z();
    topo.visibleVectorPSum = pSum.Mag();
    kptopo::ProjectOnQ(pSum, qLab, topo.visiblePParallelQ, topo.visiblePTransverseQ, cThQVis);
    topo.visibleEnergyPionHyp = ePi;
    topo.visibleEnergyProtonHyp = eP;
    const TLorentzVector LvRecVis = LvMissX - LvDet;
    topo.mRecoilVisible = LvRecVis.M();
    topo.pRecoilVisible = LvRecVis.P();
    topo.eRecoilVisible = LvRecVis.E();
  }
  if(topo.k0Flag){
    TLorentzVector LvK0(TVector3(topo.k0Mom_x, topo.k0Mom_y, topo.k0Mom_z),
			TMath::Hypot(topo.k0Mom, K0Mass));
    topo.mRecoilK0S = (LvMissX - LvK0).M();
    topo.pRecoilK0S = (LvMissX - LvK0).P();
  }
  // mRecoilPromptPim: missing minus all prompt pi- (excluding V0 daughters).
  {
    TLorentzVector LvDet(0,0,0,0);
    for(Int_t it=0; it<ntTpc; ++it){
      if(topo.trk_isLambdaSignalDaughter[it] || topo.trk_isK0Daughter[it]) continue;
      if(topo.trk_role[it]!=topodef::kRolePrompt) continue;
      if(topo.trk_pidClass[it]!=topodef::kPidPiMinus) continue;
      if(TMath::IsNaN(topo.trk_mom[it])) continue;
      TVector3 p(topo.trk_mom_x[it], topo.trk_mom_y[it], topo.trk_mom_z[it]);
      LvDet += TLorentzVector(p, TMath::Hypot(p.Mag(), PionMass));
    }
    topo.mRecoilPromptPim = (LvMissX - LvDet).M();
  }

  // Daughter momentum-vector sums |sum p_i|.
  if(topo.lFlag==1 && lDaughter.size()>=2){
    TVector3 sum(0,0,0);
    Bool_t ok = true;
    for(size_t i=0; i<lDaughter.size(); ++i){
      const Int_t id = lDaughter[i];
      if(id<0 || id>=ntTpc || TMath::IsNaN(topo.trk_mom[id])){ ok=false; break; }
      sum += TVector3(topo.trk_mom_x[id], topo.trk_mom_y[id], topo.trk_mom_z[id]);
    }
    if(ok) topo.pLambdaDaughters = sum.Mag();
  }
  if(topo.k0Flag==1 && topo.k0DaughterId.size()>=2){
    TVector3 sum(0,0,0);
    Bool_t ok = true;
    for(size_t i=0; i<topo.k0DaughterId.size(); ++i){
      const Int_t id = topo.k0DaughterId[i];
      if(id<0 || id>=ntTpc || TMath::IsNaN(topo.trk_mom[id])){ ok=false; break; }
      sum += TVector3(topo.trk_mom_x[id], topo.trk_mom_y[id], topo.trk_mom_z[id]);
    }
    if(ok) topo.pK0Daughters = sum.Mag();
  }

  // Forward-proton parentage: all M(p_fwd, pi-) pairs (fit mass spectrum; no yield cut).
  {
    TLorentzVector LvFwdP(kp_momTPC, TMath::Hypot(kp_momTPC.Mag(), ProtonMass));
    Double_t best_dm = 1e9;
    Int_t kuramaGf = -1;
    for(Int_t igf=0; igf<GFntTpc; ++igf) if(event.isKurama[igf]){ kuramaGf=igf; break; }
    for(Int_t it=0; it<ntTpc; ++it){
      if(topo.trk_pidClass[it]!=topodef::kPidPiMinus) continue;
      if(topo.trk_role[it]!=topodef::kRolePrompt
	 && topo.trk_role[it]!=topodef::kRoleDisplaced) continue;
      if(TMath::IsNaN(topo.trk_mom[it])) continue;
      TVector3 pimom(topo.trk_mom_x[it], topo.trk_mom_y[it], topo.trk_mom_z[it]);
      TLorentzVector LvPim(pimom, TMath::Hypot(pimom.Mag(), PionMass));
      const Double_t m = (LvFwdP + LvPim).M();
      topo.mFwdPPim.push_back(m);
      HF1(23040, m);
      topo.lFwdFlag = 1;
      const Double_t dm = TMath::Abs(m - LambdaMass);
      if(dm < best_dm){
	best_dm = dm;
	topo.lFwdMass = m;
	topo.lFwdPionId = it;
	topo.lFwdPimMom = pimom.Mag();
      }
      if(kuramaGf>=0 && event.GFfitstatus[it]){
	TVector3 p_mom, pi_mom; Double_t ppi_dist = qnan; TVector3 vtx;
	Double_t p_extrap = qnan, pi_extrap = qnan;
	if(GFtrackCont.FindVertex(kuramaGf, it, 0, 0,
				  p_extrap, pi_extrap,
				  p_mom, pi_mom, ppi_dist, vtx, vtx_scan_range)
	   && ppi_dist < GFppi_distcut){
	  topo.fwdLamDca = ppi_dist;
	  topo.fwdLamVtx_x = vtx.x(); topo.fwdLamVtx_y = vtx.y(); topo.fwdLamVtx_z = vtx.z();
	}
      }
    }
    if(topo.lFwdPionId>=0
       && TMath::Abs(topo.lFwdMass - LambdaMass) < topodef::kLambdaFwdMassWindow){
      topo.lFwdIsSignal = 1;
    }
  }

  // ===================================================================
  //  Fill topology histograms (ID 20000–29999)
  // ===================================================================
  {
    const Double_t bk = topo.BKaon;
    const Int_t iNom = 0;

    // --- Closure ---
    HF1(20000, bk);
    HF1(20002, topo.MissMassNucl);
    if(topo.escFlag){
      HF1(20001, bk);
      HF1(20003, topo.MissMassNucl);
      HF1(20010, topo.escKmMom);
      HF1(20011, topo.escKmM2);
    }
    HF1(20012, topo.pFwd);
    HF1(20013, topo.thetaFwd);

    // --- Track classification ---
    for(Int_t it=0; it<ntTpc; ++it){
      HF1(20100, topo.trk_role[it]);
      HF1(20101, topo.trk_pidClass[it]);
      const Double_t poq = (topo.trk_charge[it]!=0)
	? topo.trk_mom[it]/topo.trk_charge[it] : topo.trk_mom[it];
      if(!TMath::IsNaN(topo.trk_dEdx[it])){
	HF2(20102, poq, topo.trk_dEdx[it]);
	if(topo.trk_role[it]==topodef::kRolePrompt)
	  HF2(20103, poq, topo.trk_dEdx[it]);
      }
      if(!TMath::IsNaN(topo.trk_dcaPrimary[it]))
	HF1(20110, topo.trk_dcaPrimary[it]);
      HF1(20111, topo.trk_nclust[it]);
      if(topo.trk_htofReached[it] && !TMath::IsNaN(topo.trk_m2[it]))
	HF2(20120, topo.trk_mom[it], topo.trk_m2[it]);
    }

    // --- Primary vertex ---
    HF1(20200, topo.prodVtx_x);
    HF1(20201, topo.prodVtx_y);
    HF1(20202, topo.prodVtx_z);

    // --- Multiplicity (nominal = index 0) ---
    HF1(20300, topo.nchRaw[iNom]);
    HF1(20301, topo.nchPrompt[iNom]);
    HF1(20302, topo.nProton[iNom]);
    HF1(20303, topo.nPiPlus[iNom]);
    HF1(20304, topo.nPiMinus[iNom]);
    HF1(20305, topo.nDisplaced[iNom]);
    HF1(20306, topo.topoClass[iNom]);
    HF2(20310, topo.nchPrompt[iNom], bk);
    HF2(20311, topo.topoClass[iNom], bk);

    for(int ith=0; ith<topodef::kNThr; ++ith){
      HF1(20320+ith, topo.nchPrompt[ith]);
      HF2(20340+ith, topo.nchPrompt[ith], bk);
    }

    // --- Physics: B_K-resolved multiplicity ---
    for(int ibk=0; ibk<5; ++ibk){
      const Double_t kBK5Edge[6] = {-0.30,-0.15,-0.05,0.00,0.05,0.30};
      if(bk >= kBK5Edge[ibk] && bk < kBK5Edge[ibk+1]){
	HF1(20500+ibk, topo.nchPrompt[iNom]);
	HF1(20510+ibk, topo.topoClass[iNom]);
	break;
      }
    }

    // --- Lambda ---
    if(topo.lFlag){
      HF1(20600, topo.lMass);
      HF1(20601, topo.lDecayLen);
      if(topo.lIsSignal) HF1(20610, bk);
      if(topo.lIsSideband) HF1(20611, bk);
    }
    HF2(20602, topo.lFlag, bk);

    // --- M(p_fwd pi-) ---
    for(size_t i=0; i<topo.mFwdPPim.size(); ++i){
      HF1(20700, topo.mFwdPPim[i]);
      HF2(20701, topo.mFwdPPim[i], bk);
    }

    // --- M(Lambda p), cos(theta) ---
    if(topo.lFlag && topo.addProtonId>=0){
      HF1(20710, topo.mLambdaP);
      HF1(20711, topo.cosThetaLambdaP);
      HF2(20712, topo.mLambdaP, bk);
    }

    // --- 21000+: exclusive / semi-tag QA ---
    HF1(21000, topo.lambdaRegion);
    HF1(21001, topo.nK0S);
    if(topo.k0Flag){
      HF1(21002, topo.k0Mass);
      HF2(21004, topo.k0Mass, bk);
    }

    HF1(21010, topo.nPromptProton[iNom]);
    HF1(21011, topo.nPromptPiPlus[iNom]);
    HF1(21012, topo.nPromptPiMinus[iNom]);
    HF1(21013, topo.nPromptKMinus[iNom]);
    HF1(21014, topo.nPromptElectron[iNom]);
    HF1(21015, topo.nPromptAmbiguous[iNom]);
    HF1(21016, topo.nPromptUnknown[iNom]);
    HF1(21020, topo.nDisplacedProton[iNom]);
    HF1(21021, topo.nDisplacedPiPlus[iNom]);
    HF1(21022, topo.nDisplacedPiMinus[iNom]);
    HF1(21023, topo.nDisplacedKMinus[iNom]);
    HF1(21024, topo.nDisplacedElectron[iNom]);
    HF1(21025, topo.nDisplacedAmbiguous[iNom]);
    HF1(21026, topo.nDisplacedUnknown[iNom]);

    {
      const Int_t sumP = topo.nPromptProton[iNom]+topo.nPromptPiPlus[iNom]
	+topo.nPromptPiMinus[iNom]+topo.nPromptKMinus[iNom]
	+topo.nPromptElectron[iNom]+topo.nPromptAmbiguous[iNom]
	+topo.nPromptUnknown[iNom];
      const Int_t sumD = topo.nDisplacedProton[iNom]+topo.nDisplacedPiPlus[iNom]
	+topo.nDisplacedPiMinus[iNom]+topo.nDisplacedKMinus[iNom]
	+topo.nDisplacedElectron[iNom]+topo.nDisplacedAmbiguous[iNom]
	+topo.nDisplacedUnknown[iNom];
      HF1(21030, sumP - topo.nchPrompt[iNom]);
      HF1(21031, sumD - topo.nDisplaced[iNom]);
      HF1(21032, topo.nPromptClosureDiff[iNom]);
      HF1(21033, topo.nDisplacedClosureDiff[iNom]);
    }

    gExclKeyCount[topo.exclusiveVisibleKey[iNom]]++;
    for(Int_t f=0; f<16; ++f){
      HF1(23000+f, (Double_t)topodef::DecodeExclusiveField(topo.exclusiveVisibleKey[iNom], f));
    }
    Int_t exclNamed = topodef::ClassifyExclusiveNamed(
      topo.nLambdaSignal, topo.nK0S,
      topo.nPromptProton[iNom], topo.nPromptPiPlus[iNom], topo.nPromptPiMinus[iNom],
      topo.nPromptKMinus[iNom], topo.nPromptElectron[iNom],
      topo.nPromptAmbiguous[iNom], topo.nPromptUnknown[iNom],
      topo.nDisplacedProton[iNom], topo.nDisplacedPiPlus[iNom], topo.nDisplacedPiMinus[iNom],
      topo.nDisplacedKMinus[iNom], topo.nDisplacedElectron[iNom],
      topo.nDisplacedAmbiguous[iNom], topo.nDisplacedUnknown[iNom]);
    HF1(21450, exclNamed);
    HF2(21451, exclNamed, bk);
    HF1(21500 + exclNamed, bk);

    const ULong64_t bits = topo.semiTagBits[iNom];
    Int_t nBits = 0;
    for(Int_t b=0; b<topodef::kNSemiTagBit; ++b){
      if(bits & topodef::SemiTagMask((topodef::ESemiTagBit)b)){
	++nBits;
	HF1(21100+b, bk);
	HF1(21200, b);
      }
    }
    HF1(21050, nBits);
    for(Int_t b1=0; b1<topodef::kNSemiTagBit; ++b1){
      if(!(bits & topodef::SemiTagMask((topodef::ESemiTagBit)b1))) continue;
      for(Int_t b2=0; b2<topodef::kNSemiTagBit; ++b2){
	if(!(bits & topodef::SemiTagMask((topodef::ESemiTagBit)b2))) continue;
	HF2(21060, b1, b2);
      }
    }

    if(topo.lIsSignal) HF1(21070, bk);
    if(topo.lIsSideband) HF1(21071, bk);
    if(bits & topodef::SemiTagMask(topodef::kTagLambdaPiMinus)) HF1(21072, bk);
    if(bits & topodef::SemiTagMask(topodef::kTagLambdaSidebandPiMinus)) HF1(21073, bk);
    if(bits & topodef::SemiTagMask(topodef::kTagLambdaProton)) HF1(21074, bk);
    if(bits & topodef::SemiTagMask(topodef::kTagLambdaSidebandProton)) HF1(21075, bk);

    for(int ith=0; ith<topodef::kNThr; ++ith){
      Int_t nOn = 0;
      for(Int_t b=0; b<topodef::kNSemiTagBit; ++b)
	if(topo.semiTagBits[ith] & topodef::SemiTagMask((topodef::ESemiTagBit)b)) ++nOn;
      HF1(21300+ith, nOn);
      HF1(21400+ith, topo.topoClass[ith]);
    }

    // --- 22000+: kinematics / B_K momentum / forward-Lambda ---
    if(!TMath::IsNaN(topo.pLambdaDaughters)) HF1(22000, topo.pLambdaDaughters);
    if(!TMath::IsNaN(topo.pK0Daughters)) HF1(22001, topo.pK0Daughters);
    if(topo.lIsSignal && !TMath::IsNaN(topo.lMom)) HF1(22002, topo.lMom);
    if(topo.k0Flag && !TMath::IsNaN(topo.k0Mom)) HF1(22003, topo.k0Mom);

    if(topo.addProtonId>=0 && !TMath::IsNaN(topo.pLambdaP)){
      if(topo.lIsSignal){ HF1(22010, topo.pLambdaP); HF1(22020, topo.cosThetaLambdaP); HF1(22030, topo.mLambdaP); }
      else if(topo.lIsSideband){ HF1(22011, topo.pLambdaP); HF1(22021, topo.cosThetaLambdaP); HF1(22031, topo.mLambdaP); }
      else if(topo.lFlag){ HF1(22012, topo.pLambdaP); HF1(22022, topo.cosThetaLambdaP); HF1(22032, topo.mLambdaP); }
    }
    if(topo.addPimId>=0 && !TMath::IsNaN(topo.pLambdaPim)){
      if(topo.lIsSignal){ HF1(22040, topo.pLambdaPim); HF1(22050, topo.cosThetaLambdaPim); HF1(22060, topo.mLambdaPim); }
      else if(topo.lIsSideband){ HF1(22041, topo.pLambdaPim); HF1(22051, topo.cosThetaLambdaPim); HF1(22061, topo.mLambdaPim); }
      else if(topo.lFlag){ HF1(22042, topo.pLambdaPim); HF1(22052, topo.cosThetaLambdaPim); HF1(22062, topo.mLambdaPim); }
    }

    {
      const Int_t ibk = topodef::BKMomBin(bk);
      if(topo.lIsSignal && !TMath::IsNaN(topo.lMom)){
	HF2(22200, topo.lMom, bk);
	if(ibk>=0) HF1(22100+ibk, topo.lMom);
      }
      if(topo.k0Flag && !TMath::IsNaN(topo.k0Mom)){
	HF2(22201, topo.k0Mom, bk);
	if(ibk>=0) HF1(22110+ibk, topo.k0Mom);
      }
      for(Int_t it=0; it<ntTpc; ++it){
	if(topo.trk_isLambdaSignalDaughter[it] || topo.trk_isK0Daughter[it]) continue;
	if(topo.trk_role[it]!=topodef::kRolePrompt
	   && topo.trk_role[it]!=topodef::kRoleDisplaced) continue;
	if(TMath::IsNaN(topo.trk_mom[it])) continue;
	if(topo.trk_mom[it] < topodef::kThrPMin[0]) continue;
	const Int_t pid = topo.trk_pidClass[it];
	if(pid==topodef::kPidProton){
	  HF2(22202, topo.trk_mom[it], bk);
	  if(ibk>=0) HF1(22120+ibk, topo.trk_mom[it]);
	} else if(pid==topodef::kPidPiPlus){
	  HF2(22203, topo.trk_mom[it], bk);
	  if(ibk>=0) HF1(22130+ibk, topo.trk_mom[it]);
	} else if(pid==topodef::kPidPiMinus){
	  HF2(22204, topo.trk_mom[it], bk);
	  if(ibk>=0) HF1(22140+ibk, topo.trk_mom[it]);
	}
      }
    }

    if(!TMath::IsNaN(topo.lFwdMass)) HF1(22300, topo.lFwdMass);
    if(topo.lFwdIsSignal){
      HF1(22301, bk);
      HF1(22310, topo.nchPrompt[iNom]);
      HF1(22311, topo.nPromptPiMinus[iNom]);
      HF1(22312, topo.nPromptProton[iNom]);
      HF1(22313, topo.nPromptPiPlus[iNom]);
      HF1(22314, topo.nDisplaced[iNom]);
      if(exclNamed>=0) HF1(22320, exclNamed);
      for(Int_t b=0; b<topodef::kNSemiTagBit; ++b)
	if(bits & topodef::SemiTagMask((topodef::ESemiTagBit)b)) HF1(22350+b, bk);
      if(topo.nPromptPiMinus[iNom]==1){
	HF1(22340, bk);
	if(exclNamed>=0) HF1(22341, exclNamed);
      }
    } else {
      HF1(22302, bk);
      if(topo.nPromptPiMinus[iNom]==1) HF1(22342, bk);
    }

    // --- 23000+: phase-3 topology / V0 / kinematics QA [raw] ---
    if(topo.vtxStatus==3 && !TMath::IsNaN(topo.vtxDeltaMag)){
      HF1(23100, topo.vtxDeltaMag);
      HF2(23101, topo.vtxDeltaMag, bk);
      HF2(23102, topo.vtxDeltaMag, topo.prodVtx_z);
    }
    HF1(23110, topo.vtxStatus);
    HF2(23111, topo.topoClass[iNom], topo.vtxDeltaMag);
    for(Int_t ic=0; ic<topo.nLamCand; ++ic){
      HF1(23200, topo.lamCand_mass[ic]);
      HF1(23201, topo.lamCand_dauDca[ic]);
      HF1(23202, topo.lamCand_decayLen[ic]);
      HF1(23203, topo.lamCand_cosPoint[ic]);
      if(topo.lamCand_region[ic]==1) HF1(23210, topo.lamCand_mass[ic]);
      else if(topo.lamCand_region[ic]>=2) HF1(23211, topo.lamCand_mass[ic]);
    }
    for(Int_t ic=0; ic<topo.nK0Cand; ++ic){
      HF1(23220, topo.k0Cand_mass[ic]);
      if(topo.k0Cand_region[ic]==1) HF1(23221, topo.k0Cand_mass[ic]);
      else if(topo.k0Cand_region[ic]>=2) HF1(23222, topo.k0Cand_mass[ic]);
    }
    for(Int_t it=0; it<ntTpc; ++it){
      if(topo.trk_role[it]!=topodef::kRolePrompt && topo.trk_role[it]!=topodef::kRoleDisplaced) continue;
      if(!TMath::IsNaN(topo.trk_pParallelQ[it])){
	HF2(23300, topo.trk_pParallelQ[it], topo.trk_pTransverseQ[it]);
	HF2(23301, bk, topo.trk_pstar[it]);
      }
      if(topo.trk_pidClass[it]==topodef::kPidDeuteron || topo.trk_pidClass[it]==topodef::kPidTriton){
	HF2(23310, bk, topo.trk_pstar[it]);
	if(!TMath::IsNaN(topo.trk_nsigma_deutron[it])) HF1(23311, topo.trk_nsigma_deutron[it]);
	if(!TMath::IsNaN(topo.trk_nsigma_triton[it])) HF1(23312, topo.trk_nsigma_triton[it]);
      }
    }
    if(!TMath::IsNaN(topo.mRecoilVisible)) HF1(23400, topo.mRecoilVisible);
    if(!TMath::IsNaN(topo.mRecoilLambda)) HF1(23401, topo.mRecoilLambda);
    if(!TMath::IsNaN(topo.mRecoilLambdaP)) HF2(23402, topo.mRecoilLambdaP, bk);
    HF1(23410, topo.visibleNetCharge);
    HF2(23411, topo.visiblePParallelQ, topo.visiblePTransverseQ);
    HF2(23500, (Double_t)topo.runnum, bk);
    HF1(23501, topo.phiFwd);
    HF2(23502, topo.prodVtx_z, bk);

    // --- π charge asymmetry [raw, detector-level; no PID efficiency] ---
    {
      const Double_t npimP = topo.nPromptPiMinus[iNom];
      const Double_t npipP = topo.nPromptPiPlus[iNom];
      const Double_t npimD = topo.nDisplacedPiMinus[iNom];
      const Double_t npipD = topo.nDisplacedPiPlus[iNom];
      HF2(23600, bk, npimP);
      HF2(23601, bk, npipP);
      HF2(23602, bk, npimD);
      HF2(23603, bk, npipD);
      if(npimP+npipP > 0){
	HF2(23611, bk, (npimP-npipP)/(npimP+npipP));
	if(npipP > 0) HF2(23612, bk, npimP/npipP);
      }
      if(npimD+npipD > 0){
	HF2(23613, bk, (npimD-npipD)/(npimD+npipD));
	if(npipD > 0) HF2(23614, bk, npimD/npipD);
      }
    }

    // --- escape / non-escape exclusive named [raw] ---
    if(exclNamed >= 0){
      if(topo.escFlag){ HF1(23620, exclNamed); HF2(23621, exclNamed, bk); }
      else            { HF1(23622, exclNamed); HF2(23623, exclNamed, bk); }
    }

    // --- WP deep/QF nchPrompt for double-ratio [raw] ---
    for(Int_t ith=0; ith<topodef::kNThr; ++ith){
      if(bk > 0.030){
	HF1(23630+ith, topo.nchPrompt[ith]);
	if(topo.nchPrompt[ith]==0) HF1(23650+ith, 1);
      } else if(bk < 0.0){
	HF1(23640+ith, topo.nchPrompt[ith]);
	if(topo.nchPrompt[ith]==0) HF1(23660+ith, 1);
      }
    }

    // --- trigger-split B_K [raw; weight via PSfac in analysis] ---
    if(topo.trigA) HF1(23670, bk);
    if(topo.trigB) HF1(23671, bk);

    // --- Λp all-pair projections [raw] ---
    for(size_t ip=0; ip<topo.mLambdaPAll.size(); ++ip){
      if(topo.lIsSignal){
	HF2(23680, bk, topo.mLambdaPAll[ip]);
	HF2(23681, topo.mLambdaPAll[ip], topo.pLambdaPAll[ip]);
	HF2(23682, topo.mLambdaPAll[ip], topo.cosThetaLambdaPAll[ip]);
	if(ip < topo.mRecoilLambdaPAll.size())
	  HF2(23683, bk, topo.mRecoilLambdaPAll[ip]);
      } else if(topo.lIsSideband){
	HF2(23684, bk, topo.mLambdaPAll[ip]);
      }
    }
  }

  GFtrackCont.Clear();
  
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose( void )
{
  // Write exclusive-key frequency summary before the main Write().
  if(TFileCont[kOutFile]){
    TFileCont[kOutFile]->cd();
    TTree *tmap = new TTree("exclKeyMap",
			    "visible-exclusive key frequency (nominal WP)");
    ULong64_t key = 0;
    Long64_t count = 0;
    Int_t rank = 0;
    tmap->Branch("key", &key);
    tmap->Branch("count", &count);
    tmap->Branch("rank", &rank);
    // Sort by descending count for stable rank.
    std::vector<std::pair<Long64_t, ULong64_t>> sorted;
    sorted.reserve(gExclKeyCount.size());
    for(const auto& kv : gExclKeyCount)
      sorted.emplace_back(kv.second, kv.first);
    std::sort(sorted.begin(), sorted.end(),
	      [](const auto& a, const auto& b){
		if(a.first != b.first) return a.first > b.first;
		return a.second < b.second;
	      });
    for(size_t i=0; i<sorted.size(); ++i){
      rank = (Int_t)i;
      count = sorted[i].first;
      key = sorted[i].second;
      tmap->Fill();
    }
  }

  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  std::cout << "#D exclKeyMap entries : " << gExclKeyCount.size() << std::endl;
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

  if (TFileCont[kOutFile]) {
      TFileCont[kOutFile]->cd();
  } else {
      std::cerr << "!!! ConfMan::InitializeHistograms: Output file (TFileCont[kOutFile]) is not open!" << std::endl;
  }
  
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
  HB1(15, "M2 [<1.0GeV/c];  #it{M^{2}} [GeV]; counts ", nbinmass2, minmass2, maxmass2);
  HB1(18, "Kaon M2 [<1.0GeV/c][3.0#sigma dEdx cut];  #it{M^{2}} [GeV]; counts ", nbinmass2, minmass2, maxmass2);    
  HB1(19, "Kaon M2 [<1.0GeV/c][1.7#sigma dEdx cut];  #it{M^{2}} [GeV]; counts ", nbinmass2, minmass2, maxmass2);  
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
  for(int i=0; i<10; i++){
    HB1(130+i, Form("Kaon M2*charge (%f<mom<%f [GeV/c]); #it{M^{2}} [GeV/#it{c^{2}}]; counts ",double(i)*0.1, double(i+1)*0.1), nbinmass2, minmass2, maxmass2);
  }

  for(int i=0; i<numbinbek; i++){
    HB1(201+i, Form("M2 (%.2f<-BEk<%.2f [GeV]); MassSquare [GeV]; counts ",double(i)*onebinbek+minbek, double(i+1)*onebinbek+minbek), nbinmass2, minmass2, maxmass2);
  }

  HB1(500, " Theta of KP [deg]; #theta [deg]; counts", 300, 0, 30);

  HB1(600, "[Inc] Vtx [mm]; vertex X [mm]; counts", 300, -150, 150);
  HB1(601, "[Inc] Vty [mm]; vertex Y [mm]; counts", 300, -150, 150);
  HB1(602, "[Inc] Vtz [mm]; vertex Z [mm]; counts", 800, -300, 500);
  for(int i=0; i<12; i++){
    HB1(610+i, Form("[Inc] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);
  }
  HB1(630, "[Inc] Vtx [mm] (3.5< #theta_{#it{Kp}} <4.5); vertex X [mm]; counts", 600, -300, 300);
  HB1(631, "[Inc] Vty [mm] (3.5< #theta_{#it{Kp}} <4.5); vertex Y [mm]; counts", 600, -300, 300);    
  HB1(632, "[Inc] Vtz [mm] (3.5< #theta_{#it{Kp}} <4.5); vertex Z [mm]; counts", 800, -300, 500);
  HB1(650, "[Inc][TPC] Vtx [mm]; vertex X [mm]; counts", 300, -150, 150);
  HB1(651, "[Inc][TPC] Vty [mm]; vertex Y [mm]; counts", 300, -150, 150);
  HB1(652, "[Inc][TPC] Vtz [mm]; vertex Z [mm]; counts", 800, -300, 500);
  for(int i=0; i<12; i++){
    HB1(660+i, Form("[Inc][TPC] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);
  }
  HB1(680, "[Inc][TPC] Vtx [mm] (3.5< #theta_{#it{Kp}} <4.5); vertex X [mm]; counts", 600, -300, 300);
  HB1(681, "[Inc][TPC] Vty [mm] (3.5< #theta_{#it{Kp}} <4.5); vertex Y [mm]; counts", 600, -300, 300);    
  HB1(682, "[Inc][TPC] Vtz [mm] (3.5< #theta_{#it{Kp}} <4.5); vertex Z [mm]; counts", 800, -300, 500);  

  HB2(700, "[Inc] Vtx [mm] vs Vty [mm]; vertex X [mm]; vertex Y [mm]", 300, -150, 150, 300, -150, 150);
  HB2(701, "[Inc] Vty [mm] vs Vtz [mm]; vertex Y [mm]; vertex Z [mm]", 300, -150, 150, 800, -300, 500);
  HB2(702, "[Inc] Vtz [mm] vs Vtx [mm]; vertex Z [mm]; vertex X [mm]", 800, -300, 500, 300, -150, 150);  
  HB2(750, "[Inc][TPC] Vtx [mm] vs Vty [mm]; vertex X [mm]; vertex Y [mm]", 300, -150, 150, 300, -150, 150);
  HB2(751, "[Inc][TPC] Vty [mm] vs Vtz [mm]; vertex Y [mm]; vertex Z [mm]", 300, -150, 150, 800, -300, 500);
  HB2(752, "[Inc][TPC] Vtz [mm] vs Vtx [mm]; vertex Z [mm]; vertex X [mm]", 800, -300, 500, 300, -150, 150);

  for(int i=0; i<12; i++){
    HB1(800+i, Form("[Inc][TPC][N=0] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);
    HB1(820+i, Form("[Inc][TPC][N=1] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);
    HB1(840+i, Form("[Inc][TPC][N=2] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);
    HB1(860+i, Form("[Inc][TPC][N>0] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);    
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

  HB2(5100, "[M2] Mom vs -B_{K} scatK-; -#it{B_{K}} [GeV]; #it{P_{Scat K^{-}}} [GeV/c]", 120, -0.3, 0.3, 500, 0, 1.0);

  HB1(10600, "[Exc] Vtx [mm]; vertex X [mm]; counts", 300, -150, 150);
  HB1(10601, "[Exc] Vty [mm]; vertex Y [mm]; counts", 300, -150, 150);
  HB1(10602, "[Exc] Vtz [mm]; vertex Z [mm]; counts", 800, -300, 500);
  for(int i=0; i<12; i++){
    HB1(10610+i, Form("[Exc] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);
  }
  
  HB2(10700, "[Exc] Vtx [mm] vs Vty [mm]; vertex X [mm]; vertex Y [mm]", 300, -150, 150, 300, -150, 150);
  HB2(10701, "[Exc] Vty [mm] vs Vtz [mm]; vertex Y [mm]; vertex Z [mm]", 300, -150, 150, 800, -300, 500);
  HB2(10702, "[Exc] Vtz [mm] vs Vtx [mm]; vertex Z [mm]; vertex X [mm]", 800, -300, 500, 300, -150, 150);    

  HB1(10650, "[Exc] Vtx [mm]; vertex X [mm]; counts", 300, -150, 150);
  HB1(10651, "[Exc] Vty [mm]; vertex Y [mm]; counts", 300, -150, 150);
  HB1(10652, "[Exc] Vtz [mm]; vertex Z [mm]; counts", 800, -300, 500);
  for(int i=0; i<12; i++){
    HB1(10660+i, Form("[Exc] Vtz [mm] (%.1f<#theta_{#it{Kp}}<%.1f); vertex Z [mm]; counts",double(i)*2.0,double(i+1)*2.0), 800, -300, 500);
  }  
  HB2(10750, "[Exc] Vtx [mm] vs Vty [mm]; vertex X [mm]; vertex Y [mm]", 300, -150, 150, 300, -150, 150);
  HB2(10751, "[Exc] Vty [mm] vs Vtz [mm]; vertex Y [mm]; vertex Z [mm]", 300, -150, 150, 800, -300, 500);
  HB2(10752, "[Exc] Vtz [mm] vs Vtx [mm]; vertex Z [mm]; vertex X [mm]", 800, -300, 500, 300, -150, 150);    
  
  HB1(13900, " BE KP inclusive; #minusB_{K} [GeV]; Counts", 120, -0.3, 0.3);
  HB1(13910, " BE KP inclusive (all angle); #minusB_{K} [GeV]; Counts", 120, -0.3, 0.3);
  HB1(13950, " BE KP exclusive [dEdxPID]; #minusB_{K} [GeV]; Counts", 120, -0.3, 0.3);
  HB1(13951, " BE KP exclusive [M2PID]; #minusB_{K} [GeV]; Counts", 120, -0.3, 0.3);

  for(int i=0; i<120; i++){
    HB1(14000+i, Form("[M2][bek] (MissMom - Mom) of scatK- (%.3f<-BE_K<%.3f); momentum [GeV/c]; coutns", double(i)*0.005-0.300, double(i+1)*0.005-0.300), 500, 0., 1.0);
  }  

  // ===================================================================
  //  Topology histograms  (ID 20000 – 29999)
  //  Filled in the event loop so that full-run hadd is sufficient.
  // ===================================================================
  const Int_t    kBKBins = 120;
  const Double_t kBKLo   = -0.30;
  const Double_t kBKHi   =  0.30;
  const Int_t    kNBK5   = 5; // coarse bins for 2D
  const Double_t kBK5Edge[kNBK5+1] = {-0.30,-0.15,-0.05,0.00,0.05,0.30};

  // --- 20000–20099: Closure ---
  HB1(20000, "[topo] B_{K} inclusive; B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(20001, "[topo] B_{K} K-escape; B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(20002, "[topo] MM inclusive; MM [GeV/#it{c}^{2}]; Counts", 300, 10.0, 12.0);
  HB1(20003, "[topo] MM K-escape; MM [GeV/#it{c}^{2}]; Counts", 300, 10.0, 12.0);
  HB1(20010, "[topo] K- escape mom; p [GeV/#it{c}]; Counts", 500, 0., 1.0);
  HB1(20011, "[topo] K- escape M^{2}; M^{2} [GeV^{2}]; Counts", 200, -0.5, 0.5);
  HB1(20012, "[topo] forward p mom; p [GeV/#it{c}]; Counts", 500, 0., 2.0);
  HB1(20013, "[topo] #theta_{Kp}; #theta [deg]; Counts", 300, 0., 30.);

  // --- 20100–20299: Track classification sanity ---
  HB1(20100, "[topo] track role; role; Counts", topodef::kNTrackRole, 0, topodef::kNTrackRole);
  HB1(20101, "[topo] PID class; pid; Counts", topodef::kNPidClass, 0, topodef::kNPidClass);
  HB2(20102, "[topo] dE/dx vs p/q (all); p/q [GeV/#it{c}]; dE/dx [arb.]",
       200, -1.0, 1.0, 200, 0., 350.);
  HB2(20103, "[topo] dE/dx vs p/q (prompt); p/q [GeV/#it{c}]; dE/dx [arb.]",
       200, -1.0, 1.0, 200, 0., 350.);
  HB1(20110, "[topo] DCA to prim vtx; DCA [mm]; Counts", 200, 0., 100.);
  HB1(20111, "[topo] nClust; nClust; Counts", 100, 0, 100);
  HB2(20120, "[topo] M^{2} vs p (HTOF-matched); p [GeV/#it{c}]; M^{2} [GeV^{2}]",
       200, 0., 1.0, 200, -0.5, 1.5);

  // --- 20200–20299: Primary vertex ---
  HB1(20200, "[topo] vtx_x; x [mm]; Counts", 200, -100, 100);
  HB1(20201, "[topo] vtx_y; y [mm]; Counts", 200, -100, 100);
  HB1(20202, "[topo] vtx_z; z [mm]; Counts", 400, -200, 200);

  // --- 20300–20499: Multiplicity (nominal working point = index 0) ---
  HB1(20300, "[topo] nch raw (nom); N_{ch}^{raw}; Counts", 15, 0, 15);
  HB1(20301, "[topo] nch prompt (nom); N_{ch}^{prompt}; Counts", 15, 0, 15);
  HB1(20302, "[topo] n proton (nom); N_{p}; Counts", 10, 0, 10);
  HB1(20303, "[topo] n #pi^{+} (nom); N_{#pi^{+}}; Counts", 10, 0, 10);
  HB1(20304, "[topo] n #pi^{-} (nom); N_{#pi^{-}}; Counts", 10, 0, 10);
  HB1(20305, "[topo] n displaced (nom); N_{disp}; Counts", 10, 0, 10);
  HB1(20306, "[topo] topo class (nom); class; Counts",
       topodef::kNTopoClass, 0, topodef::kNTopoClass);
  HB2(20310, "[topo] B_{K} vs nch prompt (nom); N_{ch}^{prompt}; B_{K} [GeV]",
       15, 0, 15, kBKBins, kBKLo, kBKHi);
  HB2(20311, "[topo] B_{K} vs topo class (nom); class; B_{K} [GeV]",
       topodef::kNTopoClass, 0, topodef::kNTopoClass, kBKBins, kBKLo, kBKHi);

  for(int ith=0; ith<topodef::kNThr; ++ith){
    HB1(20320+ith, Form("[topo] nch prompt (thr %d); N_{ch}^{prompt}; Counts", ith),
	 15, 0, 15);
    HB2(20340+ith, Form("[topo] B_{K} vs nch prompt (thr %d); N_{ch}^{prompt}; B_{K} [GeV]", ith),
	 15, 0, 15, kBKBins, kBKLo, kBKHi);
  }

  // --- 20500–20699: Physics-ready ---
  // B_K in coarse B_K bins x multiplicity
  for(int ibk=0; ibk<kNBK5; ++ibk){
    HB1(20500+ibk, Form("[topo] nch prompt (%.2f<B_{K}<%.2f); N_{ch}^{prompt}; Counts",
			 kBK5Edge[ibk], kBK5Edge[ibk+1]), 15, 0, 15);
    HB1(20510+ibk, Form("[topo] topo class (%.2f<B_{K}<%.2f); class; Counts",
			 kBK5Edge[ibk], kBK5Edge[ibk+1]),
	 topodef::kNTopoClass, 0, topodef::kNTopoClass);
  }

  // Lambda
  HB1(20600, "[topo] #Lambda mass; M(p#pi^{-}) [GeV/#it{c}^{2}]; Counts", 200, 1.05, 1.25);
  HB1(20601, "[topo] #Lambda decay length; L [mm]; Counts", 200, 0., 200.);
  HB2(20602, "[topo] B_{K} vs #Lambda flag; #Lambda flag; B_{K} [GeV]",
       3, 0, 3, kBKBins, kBKLo, kBKHi);
  HB1(20610, "[topo] B_{K} (#Lambda signal); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(20611, "[topo] B_{K} (#Lambda sideband); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);

  // M(p_fwd pi-)
  HB1(20700, "[topo] M(p_{fwd}#pi^{-}); M [GeV/#it{c}^{2}]; Counts", 200, 1.0, 1.5);
  HB2(20701, "[topo] B_{K} vs M(p_{fwd}#pi^{-}); M [GeV/#it{c}^{2}]; B_{K} [GeV]",
       200, 1.0, 1.5, kBKBins, kBKLo, kBKHi);

  // M(Lambda p), cos(theta)
  HB1(20710, "[topo] M(#Lambda p); M [GeV/#it{c}^{2}]; Counts", 200, 2.0, 3.0);
  HB1(20711, "[topo] cos#theta_{#Lambda p}; cos#theta; Counts", 100, -1., 1.);
  HB2(20712, "[topo] B_{K} vs M(#Lambda p); M [GeV/#it{c}^{2}]; B_{K} [GeV]",
       200, 2.0, 3.0, kBKBins, kBKLo, kBKHi);

  // --- 21000–21499: visible-exclusive / semi-tag QA ---
  HB1(21000, "[topo] lambdaRegion; region; Counts", 3, 0, 3);
  HB1(21001, "[topo] nK0S; nK0S; Counts", 5, 0, 5);
  HB1(21002, "[topo] K0S mass (selected); M(#pi^{+}#pi^{-}) [GeV/#it{c}^{2}]; Counts", 200, 0.4, 0.6);
  HB1(21003, "[topo] M(#pi^{+}#pi^{-}) all pairs; M [GeV/#it{c}^{2}]; Counts", 200, 0.3, 0.7);
  HB2(21004, "[topo] B_{K} vs K0S mass (selected); M [GeV/#it{c}^{2}]; B_{K} [GeV]",
       100, 0.4, 0.6, kBKBins, kBKLo, kBKHi);
  HB1(21005, "[topo] M(#pi^{+}#pi^{-}) best any; M [GeV/#it{c}^{2}]; Counts", 200, 0.3, 0.7);

  HB1(21010, "[topo] nPromptProton (nom); N; Counts", 10, 0, 10);
  HB1(21011, "[topo] nPromptPiPlus (nom); N; Counts", 10, 0, 10);
  HB1(21012, "[topo] nPromptPiMinus (nom); N; Counts", 10, 0, 10);
  HB1(21013, "[topo] nPromptKMinus (nom); N; Counts", 10, 0, 10);
  HB1(21014, "[topo] nPromptElectron (nom); N; Counts", 10, 0, 10);
  HB1(21015, "[topo] nPromptAmbiguous (nom); N; Counts", 10, 0, 10);
  HB1(21016, "[topo] nPromptUnknown (nom); N; Counts", 10, 0, 10);
  HB1(21020, "[topo] nDisplacedProton (nom); N; Counts", 10, 0, 10);
  HB1(21021, "[topo] nDisplacedPiPlus (nom); N; Counts", 10, 0, 10);
  HB1(21022, "[topo] nDisplacedPiMinus (nom); N; Counts", 10, 0, 10);
  HB1(21023, "[topo] nDisplacedKMinus (nom); N; Counts", 10, 0, 10);
  HB1(21024, "[topo] nDisplacedElectron (nom); N; Counts", 10, 0, 10);
  HB1(21025, "[topo] nDisplacedAmbiguous (nom); N; Counts", 10, 0, 10);
  HB1(21026, "[topo] nDisplacedUnknown (nom); N; Counts", 10, 0, 10);

  HB1(21030, "[topo] exclPromptSum - nchPrompt; diff; Counts", 21, -10, 11);
  HB1(21031, "[topo] exclDispSum - nDisplaced; diff; Counts", 21, -10, 11);
  HB1(21032, "[topo] prompt PID closure diff; diff; Counts", 11, -5, 6);
  HB1(21033, "[topo] displaced PID closure diff; diff; Counts", 11, -5, 6);

  for(int f=0; f<16; ++f){
    HB1(23000+f, Form("[topo] excl field %d (DecodeExclusiveField); mult; Counts", f),
	 8, 0, 8);
  }
  HB1(21050, "[topo] n semiTag bits (nom); nBits; Counts", 20, 0, 20);
  HB2(21060, "[topo] semiTag overlap; bit1; bit2",
       topodef::kNSemiTagBit, 0, topodef::kNSemiTagBit,
       topodef::kNSemiTagBit, 0, topodef::kNSemiTagBit);

  HB1(21070, "[topo] B_{K} (#Lambda signal); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(21071, "[topo] B_{K} (#Lambda sideband); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(21072, "[topo] B_{K} (#Lambda#pi^{-} 1NA-enriched); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(21073, "[topo] B_{K} (#Lambda_{sb}#pi^{-}); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(21074, "[topo] B_{K} (#Lambda p 2NA-enriched); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(21075, "[topo] B_{K} (#Lambda_{sb}p); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);

  for(int b=0; b<topodef::kNSemiTagBit; ++b){
    HB1(21100+b, Form("[topo] B_{K} (semiTag bit %d); B_{K} [GeV]; Counts", b),
	 kBKBins, kBKLo, kBKHi);
  }
  HB1(21200, "[topo] semiTag bit fired; bit; Counts",
       topodef::kNSemiTagBit, 0, topodef::kNSemiTagBit);

  for(int ith=0; ith<topodef::kNThr; ++ith){
    HB1(21300+ith, Form("[topo] n semiTag bits (thr %d); nBits; Counts", ith), 16, 0, 16);
    HB1(21400+ith, Form("[topo] topoClass (thr %d); class; Counts", ith),
	 topodef::kNTopoClass, 0, topodef::kNTopoClass);
  }

  // Named exclusive class frequency and B_K spectra
  HB1(21450, "[topo] exclusive named class; class; Counts",
       topodef::kNExclNamed, 0, topodef::kNExclNamed);
  HB2(21451, "[topo] B_{K} vs exclusive named class; class; B_{K} [GeV]",
       topodef::kNExclNamed, 0, topodef::kNExclNamed, kBKBins, kBKLo, kBKHi);
  const char* exclNamedTitle[topodef::kNExclNamed] = {
    "empty", "1p", "1#pi^{+}", "1#pi^{-}", "1other",
    "#Lambda only", "#Lambda#pi^{-}", "#Lambda p", "#Lambda other",
    "K0 only", "K0 p", "K0 other", "multi/other"
  };
  for(int ic=0; ic<topodef::kNExclNamed; ++ic){
    HB1(21500+ic, Form("[topo] B_{K} (excl %s); B_{K} [GeV]; Counts", exclNamedTitle[ic]),
	 kBKBins, kBKLo, kBKHi);
  }

  // --- 22000+: daughter sums, Lambda+hadron, B_K momenta, forward-Lambda ---
  HB1(22000, "[topo] |p_{#Lambda daughters}|; p [GeV/#it{c}]; Counts", 200, 0., 2.0);
  HB1(22001, "[topo] |p_{K0 daughters}|; p [GeV/#it{c}]; Counts", 200, 0., 2.0);
  HB1(22002, "[topo] p_{#Lambda} (signal); p [GeV/#it{c}]; Counts", 200, 0., 2.0);
  HB1(22003, "[topo] p_{K0S}; p [GeV/#it{c}]; Counts", 200, 0., 2.0);

  HB1(22010, "[topo] |p_{#Lambda}+p_{p}| (signal); p [GeV/#it{c}]; Counts", 200, 0., 2.5);
  HB1(22011, "[topo] |p_{#Lambda}+p_{p}| (sideband); p [GeV/#it{c}]; Counts", 200, 0., 2.5);
  HB1(22012, "[topo] |p_{#Lambda}+p_{p}| (other lFlag); p [GeV/#it{c}]; Counts", 200, 0., 2.5);
  HB1(22020, "[topo] cos#theta_{#Lambda p} (signal); cos#theta; Counts", 100, -1., 1.);
  HB1(22021, "[topo] cos#theta_{#Lambda p} (sideband); cos#theta; Counts", 100, -1., 1.);
  HB1(22022, "[topo] cos#theta_{#Lambda p} (other lFlag); cos#theta; Counts", 100, -1., 1.);
  HB1(22030, "[topo] M(#Lambda p) (signal); M [GeV/#it{c}^{2}]; Counts", 200, 2.0, 3.0);
  HB1(22031, "[topo] M(#Lambda p) (sideband); M [GeV/#it{c}^{2}]; Counts", 200, 2.0, 3.0);
  HB1(22032, "[topo] M(#Lambda p) (other lFlag); M [GeV/#it{c}^{2}]; Counts", 200, 2.0, 3.0);

  HB1(22040, "[topo] |p_{#Lambda}+p_{#pi^{-}}| (signal); p [GeV/#it{c}]; Counts", 200, 0., 2.5);
  HB1(22041, "[topo] |p_{#Lambda}+p_{#pi^{-}}| (sideband); p [GeV/#it{c}]; Counts", 200, 0., 2.5);
  HB1(22042, "[topo] |p_{#Lambda}+p_{#pi^{-}}| (other lFlag); p [GeV/#it{c}]; Counts", 200, 0., 2.5);
  HB1(22050, "[topo] cos#theta_{#Lambda #pi^{-}} (signal); cos#theta; Counts", 100, -1., 1.);
  HB1(22051, "[topo] cos#theta_{#Lambda #pi^{-}} (sideband); cos#theta; Counts", 100, -1., 1.);
  HB1(22052, "[topo] cos#theta_{#Lambda #pi^{-}} (other lFlag); cos#theta; Counts", 100, -1., 1.);
  HB1(22060, "[topo] M(#Lambda #pi^{-}) (signal); M [GeV/#it{c}^{2}]; Counts", 200, 1.2, 2.5);
  HB1(22061, "[topo] M(#Lambda #pi^{-}) (sideband); M [GeV/#it{c}^{2}]; Counts", 200, 1.2, 2.5);
  HB1(22062, "[topo] M(#Lambda #pi^{-}) (other lFlag); M [GeV/#it{c}^{2}]; Counts", 200, 1.2, 2.5);

  for(int ibk=0; ibk<topodef::kNBKMom; ++ibk){
    HB1(22100+ibk, Form("[topo] p_{#Lambda} (%.2f<B_{K}<%.2f); p [GeV/#it{c}]; Counts",
			 topodef::kBKMomEdge[ibk], topodef::kBKMomEdge[ibk+1]), 100, 0., 1.5);
    HB1(22110+ibk, Form("[topo] p_{K0S} (%.2f<B_{K}<%.2f); p [GeV/#it{c}]; Counts",
			 topodef::kBKMomEdge[ibk], topodef::kBKMomEdge[ibk+1]), 100, 0., 1.5);
    HB1(22120+ibk, Form("[topo] p_{p non-dau} (%.2f<B_{K}<%.2f); p [GeV/#it{c}]; Counts",
			 topodef::kBKMomEdge[ibk], topodef::kBKMomEdge[ibk+1]), 100, 0., 1.5);
    HB1(22130+ibk, Form("[topo] p_{#pi^{+} non-dau} (%.2f<B_{K}<%.2f); p [GeV/#it{c}]; Counts",
			 topodef::kBKMomEdge[ibk], topodef::kBKMomEdge[ibk+1]), 100, 0., 1.5);
    HB1(22140+ibk, Form("[topo] p_{#pi^{-} non-dau} (%.2f<B_{K}<%.2f); p [GeV/#it{c}]; Counts",
			 topodef::kBKMomEdge[ibk], topodef::kBKMomEdge[ibk+1]), 100, 0., 1.5);
  }
  HB2(22200, "[topo] B_{K} vs p_{#Lambda}; p [GeV/#it{c}]; B_{K} [GeV]",
       100, 0., 1.5, kBKBins, kBKLo, kBKHi);
  HB2(22201, "[topo] B_{K} vs p_{K0S}; p [GeV/#it{c}]; B_{K} [GeV]",
       100, 0., 1.5, kBKBins, kBKLo, kBKHi);
  HB2(22202, "[topo] B_{K} vs p_{p non-dau}; p [GeV/#it{c}]; B_{K} [GeV]",
       100, 0., 1.5, kBKBins, kBKLo, kBKHi);
  HB2(22203, "[topo] B_{K} vs p_{#pi^{+} non-dau}; p [GeV/#it{c}]; B_{K} [GeV]",
       100, 0., 1.5, kBKBins, kBKLo, kBKHi);
  HB2(22204, "[topo] B_{K} vs p_{#pi^{-} non-dau}; p [GeV/#it{c}]; B_{K} [GeV]",
       100, 0., 1.5, kBKBins, kBKLo, kBKHi);

  HB1(22300, "[topo] M(p_{fwd}#pi^{-}) best; M [GeV/#it{c}^{2}]; Counts", 200, 1.0, 1.5);
  HB1(22301, "[topo] B_{K} (fwd-#Lambda signal); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(22302, "[topo] B_{K} (no fwd-#Lambda); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(22310, "[topo] nchPrompt (fwd-#Lambda); N; Counts", 15, 0, 15);
  HB1(22311, "[topo] nPrompt#pi^{-} (fwd-#Lambda); N; Counts", 10, 0, 10);
  HB1(22312, "[topo] nPrompt p (fwd-#Lambda); N; Counts", 10, 0, 10);
  HB1(22313, "[topo] nPrompt#pi^{+} (fwd-#Lambda); N; Counts", 10, 0, 10);
  HB1(22314, "[topo] nDisplaced (fwd-#Lambda); N; Counts", 10, 0, 10);
  HB1(22320, "[topo] excl named (fwd-#Lambda); class; Counts",
       topodef::kNExclNamed, 0, topodef::kNExclNamed);
  HB1(22340, "[topo] B_{K} (fwd-#Lambda && 1#pi^{-}); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(22341, "[topo] excl named (fwd-#Lambda && 1#pi^{-}); class; Counts",
       topodef::kNExclNamed, 0, topodef::kNExclNamed);
  HB1(22342, "[topo] B_{K} (no fwd-#Lambda && 1#pi^{-}); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  for(int b=0; b<topodef::kNSemiTagBit; ++b){
    HB1(22350+b, Form("[topo] B_{K} (fwd-#Lambda && semiTag %d); B_{K} [GeV]; Counts", b),
	 kBKBins, kBKLo, kBKHi);
  }

  // --- 23000–24599: phase-3 V0 / vertex / kinematics [raw] ---
  HB1(23100, "[raw] |#Delta vtx|; |#Delta| [mm]; Counts", 100, 0., 50.);
  HB2(23101, "[raw] B_{K} vs |#Delta vtx|; |#Delta| [mm]; B_{K} [GeV]",
       50, 0., 50., kBKBins, kBKLo, kBKHi);
  HB2(23102, "[raw] prodVtx_z vs |#Delta vtx|; z [mm]; |#Delta| [mm]",
       100, -200., 200., 50, 0., 50.);
  HB1(23110, "[raw] vtxStatus; status; Counts", 4, 0, 4);
  HB2(23111, "[raw] topoClass vs |#Delta vtx|; class; |#Delta| [mm]",
       topodef::kNTopoClass, 0, topodef::kNTopoClass, 50, 0., 50.);
  HB1(23200, "[raw] #Lambda cand mass (all); M [GeV/#it{c}^{2}]; Counts", 200, 1.05, 1.25);
  HB1(23201, "[raw] #Lambda dau DCA; DCA [mm]; Counts", 100, 0., 20.);
  HB1(23202, "[raw] #Lambda decay length; L [mm]; Counts", 200, 0., 200.);
  HB1(23203, "[raw] #Lambda cos pointing; cos; Counts", 100, 0.9, 1.01);
  HB1(23210, "[raw] #Lambda mass (signal region); M [GeV/#it{c}^{2}]; Counts", 200, 1.05, 1.25);
  HB1(23211, "[raw] #Lambda mass (sideband); M [GeV/#it{c}^{2}]; Counts", 200, 1.05, 1.25);
  HB1(23220, "[raw] K0S cand mass (all); M [GeV/#it{c}^{2}]; Counts", 200, 0.3, 0.7);
  HB1(23221, "[raw] K0S mass (signal); M [GeV/#it{c}^{2}]; Counts", 200, 0.4, 0.6);
  HB1(23222, "[raw] K0S mass (sideband); M [GeV/#it{c}^{2}]; Counts", 200, 0.3, 0.7);
  HB2(23300, "[raw] p_{#parallel q} vs p_{T q}; p [GeV/#it{c}]; p [GeV/#it{c}]",
       100, -1.5, 1.5, 100, 0., 1.5);
  HB2(23301, "[raw] B_{K} vs p^{*}; B_{K} [GeV]; p^{*} [GeV/#it{c}]",
       kBKBins, kBKLo, kBKHi, 100, 0., 1.5);
  HB2(23310, "[raw] B_{K} vs p^{*} (d/t); B_{K} [GeV]; p^{*} [GeV/#it{c}]",
       kBKBins, kBKLo, kBKHi, 100, 0., 1.5);
  HB1(23311, "[raw] n#sigma deuteron (d/t tracks); n#sigma; Counts", 100, -5., 5.);
  HB1(23312, "[raw] n#sigma triton (d/t tracks); n#sigma; Counts", 100, -5., 5.);
  HB1(23400, "[raw] M_{recoil}^{visible}; M [GeV/#it{c}^{2}]; Counts", 200, 0.5, 2.5);
  HB1(23401, "[raw] M_{recoil}^{#Lambda}; M [GeV/#it{c}^{2}]; Counts", 200, 0.5, 2.5);
  HB2(23402, "[raw] B_{K} vs M_{recoil}(#Lambda p); B_{K} [GeV]; M [GeV/#it{c}^{2}]",
       kBKBins, kBKLo, kBKHi, 200, 2.0, 3.0);
  HB1(23410, "[raw] visible net charge; Q; Counts", 11, -5, 6);
  HB2(23411, "[raw] visible p_{#parallel q} vs p_{T q}; p [GeV/#it{c}]; p [GeV/#it{c}]",
       100, -1.5, 1.5, 100, 0., 1.5);
  HB1(23020, "[raw] #Lambda cand mass (FindVertex); M [GeV/#it{c}^{2}]; Counts", 200, 1.05, 1.25);
  HB1(23040, "[raw] M(p_{fwd}#pi^{-}) all pairs; M [GeV/#it{c}^{2}]; Counts", 200, 1.0, 1.5);
  HB2(23500, "[raw] B_{K} vs run; run; B_{K} [GeV]", 200, 5500, 5700, kBKBins, kBKLo, kBKHi);
  HB1(23501, "[raw] forward #phi; #phi [deg]; Counts", 72, -180., 180.);
  HB2(23502, "[raw] prodVtx_z vs B_{K}; z [mm]; B_{K} [GeV]", 100, -200., 200., kBKBins, kBKLo, kBKHi);

  // Phase-3 extensions: charge asymmetry, escape split, WP double-ratio, Λp
  HB2(23600, "[raw] B_{K} vs N(#pi^{-} prompt); B_{K}; N", kBKBins, kBKLo, kBKHi, 8, 0, 8);
  HB2(23601, "[raw] B_{K} vs N(#pi^{+} prompt); B_{K}; N", kBKBins, kBKLo, kBKHi, 8, 0, 8);
  HB2(23602, "[raw] B_{K} vs N(#pi^{-} disp); B_{K}; N", kBKBins, kBKLo, kBKHi, 8, 0, 8);
  HB2(23603, "[raw] B_{K} vs N(#pi^{+} disp); B_{K}; N", kBKBins, kBKLo, kBKHi, 8, 0, 8);
  HB2(23611, "[raw] B_{K} vs A_#pi prompt; B_{K}; A", kBKBins, kBKLo, kBKHi, 40, -1., 1.);
  HB2(23612, "[raw] B_{K} vs R(#pi^{-}/#pi^{+}) prompt; B_{K}; R", kBKBins, kBKLo, kBKHi, 40, 0., 10.);
  HB2(23613, "[raw] B_{K} vs A_#pi disp; B_{K}; A", kBKBins, kBKLo, kBKHi, 40, -1., 1.);
  HB2(23614, "[raw] B_{K} vs R(#pi^{-}/#pi^{+}) disp; B_{K}; R", kBKBins, kBKLo, kBKHi, 40, 0., 10.);
  HB1(23620, "[raw] excl named | escape; class; Counts", topodef::kNExclNamed, 0, topodef::kNExclNamed);
  HB2(23621, "[raw] B_{K} vs excl | escape; class; B_{K}", topodef::kNExclNamed, 0, topodef::kNExclNamed, kBKBins, kBKLo, kBKHi);
  HB1(23622, "[raw] excl named | non-escape; class; Counts", topodef::kNExclNamed, 0, topodef::kNExclNamed);
  HB2(23623, "[raw] B_{K} vs excl | non-escape; class; B_{K}", topodef::kNExclNamed, 0, topodef::kNExclNamed, kBKBins, kBKLo, kBKHi);
  for(int ith=0; ith<topodef::kNThr; ++ith){
    HB1(23630+ith, Form("[raw] nchPrompt deep WP%d; N; Counts", ith), 15, 0, 15);
    HB1(23640+ith, Form("[raw] nchPrompt QF WP%d; N; Counts", ith), 15, 0, 15);
    HB1(23650+ith, Form("[raw] empty deep WP%d; dummy; Counts", ith), 2, 0, 2);
    HB1(23660+ith, Form("[raw] empty QF WP%d; dummy; Counts", ith), 2, 0, 2);
  }
  HB1(23670, "[raw] B_{K} (trigA); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB1(23671, "[raw] B_{K} (trigB); B_{K} [GeV]; Counts", kBKBins, kBKLo, kBKHi);
  HB2(23680, "[raw] B_{K} vs M(#Lambda p) all pairs signal; B_{K}; M",
       kBKBins, kBKLo, kBKHi, 100, 2.0, 3.0);
  HB2(23681, "[raw] M(#Lambda p) vs |P|; M; P", 100, 2.0, 3.0, 80, 0., 2.5);
  HB2(23682, "[raw] M(#Lambda p) vs cos#theta; M; cos", 100, 2.0, 3.0, 50, -1., 1.);
  HB2(23683, "[raw] B_{K} vs Mrecoil(#Lambda p); B_{K}; M", kBKBins, kBKLo, kBKHi, 100, 0.5, 2.5);
  HB2(23684, "[raw] B_{K} vs M(#Lambda p) sideband; B_{K}; M",
       kBKBins, kBKLo, kBKHi, 100, 2.0, 3.0);

  HBTree( "topo", "topology summary of GenfitKpTopology" );

  tree->Branch("runnum", &topo.runnum);
  tree->Branch("evnum", &topo.evnum);
  tree->Branch("evStatus", &topo.evStatus);
  tree->Branch("trigA", &topo.trigA);
  tree->Branch("trigB", &topo.trigB);
  tree->Branch("PSfacTrigA", &topo.PSfacTrigA);
  tree->Branch("PSfacTrigB", &topo.PSfacTrigB);

  tree->Branch("roleClassVersion", &topo.roleClassVersion);
  tree->Branch("exclusiveVisibleVersion", &topo.exclusiveVisibleVersion);
  tree->Branch("semiTagVersion", &topo.semiTagVersion);
  tree->Branch("v0CandidateVersion", &topo.v0CandidateVersion);
  tree->Branch("kinematicsVersion", &topo.kinematicsVersion);

  tree->Branch("BKaon", &topo.BKaon);
  tree->Branch("MissMassNucl", &topo.MissMassNucl);
  tree->Branch("thetaKP", &topo.thetaKP);
  tree->Branch("qTransfer", &topo.qTransfer);
  tree->Branch("pBeam", &topo.pBeam);
  tree->Branch("pBeam_x", &topo.pBeam_x);
  tree->Branch("pBeam_y", &topo.pBeam_y);
  tree->Branch("pBeam_z", &topo.pBeam_z);
  tree->Branch("pFwd", &topo.pFwd);
  tree->Branch("pFwd_x", &topo.pFwd_x);
  tree->Branch("pFwd_y", &topo.pFwd_y);
  tree->Branch("pFwd_z", &topo.pFwd_z);
  tree->Branch("thetaFwd", &topo.thetaFwd);
  tree->Branch("phiFwd", &topo.phiFwd);
  tree->Branch("PX_x", &topo.PX_x);
  tree->Branch("PX_y", &topo.PX_y);
  tree->Branch("PX_z", &topo.PX_z);
  tree->Branch("PX_E", &topo.PX_E);

  tree->Branch("prodVtx_x", &topo.prodVtx_x);
  tree->Branch("prodVtx_y", &topo.prodVtx_y);
  tree->Branch("prodVtx_z", &topo.prodVtx_z);
  tree->Branch("prodVtxNTrack", &topo.prodVtxNTrack);
  tree->Branch("kkVtx_x", &topo.kkVtx_x);
  tree->Branch("kkVtx_y", &topo.kkVtx_y);
  tree->Branch("kkVtx_z", &topo.kkVtx_z);
  tree->Branch("kkVtxValid", &topo.kkVtxValid);
  tree->Branch("multiVtx_x", &topo.multiVtx_x);
  tree->Branch("multiVtx_y", &topo.multiVtx_y);
  tree->Branch("multiVtx_z", &topo.multiVtx_z);
  tree->Branch("multiVtxNTrack", &topo.multiVtxNTrack);
  tree->Branch("multiVtxValid", &topo.multiVtxValid);
  tree->Branch("vtxDelta_x", &topo.vtxDelta_x);
  tree->Branch("vtxDelta_y", &topo.vtxDelta_y);
  tree->Branch("vtxDelta_z", &topo.vtxDelta_z);
  tree->Branch("vtxDeltaMag", &topo.vtxDeltaMag);
  tree->Branch("vtxStatus", &topo.vtxStatus);
  tree->Branch("dcaReferenceFallback", &topo.dcaReferenceFallback);

  tree->Branch("ntTpc", &topo.ntTpc);
  tree->Branch("GFntTpc", &topo.GFntTpc);

  tree->Branch("trk_role", &topo.trk_role);
  tree->Branch("trk_pidClass", &topo.trk_pidClass);
  tree->Branch("trk_pidbits", &topo.trk_pidbits);
  tree->Branch("trk_charge", &topo.trk_charge);
  tree->Branch("trk_nclust", &topo.trk_nclust);
  tree->Branch("trk_gffit", &topo.trk_gffit);
  tree->Branch("trk_inside", &topo.trk_inside);
  tree->Branch("trk_htofReached", &topo.trk_htofReached);
  tree->Branch("trk_isLambdaCandidateDaughter", &topo.trk_isLambdaCandidateDaughter);
  tree->Branch("trk_isLambdaSignalDaughter", &topo.trk_isLambdaSignalDaughter);
  tree->Branch("trk_isLambdaSidebandDaughter", &topo.trk_isLambdaSidebandDaughter);
  tree->Branch("trk_isK0Daughter", &topo.trk_isK0Daughter);
  tree->Branch("trk_chisqr", &topo.trk_chisqr);
  tree->Branch("trk_dEdx", &topo.trk_dEdx);
  tree->Branch("trk_mom", &topo.trk_mom);
  tree->Branch("trk_mom_x", &topo.trk_mom_x);
  tree->Branch("trk_mom_y", &topo.trk_mom_y);
  tree->Branch("trk_mom_z", &topo.trk_mom_z);
  tree->Branch("trk_theta", &topo.trk_theta);
  tree->Branch("trk_phi", &topo.trk_phi);
  tree->Branch("trk_pstar", &topo.trk_pstar);
  tree->Branch("trk_costStar", &topo.trk_costStar);
  tree->Branch("trk_dcaPrimary", &topo.trk_dcaPrimary);
  tree->Branch("trk_dcaKK", &topo.trk_dcaKK);
  tree->Branch("trk_dcaMulti", &topo.trk_dcaMulti);
  tree->Branch("trk_pParallelQ", &topo.trk_pParallelQ);
  tree->Branch("trk_pTransverseQ", &topo.trk_pTransverseQ);
  tree->Branch("trk_cosThetaQ", &topo.trk_cosThetaQ);
  tree->Branch("trk_pstarProtonHyp", &topo.trk_pstarProtonHyp);
  tree->Branch("trk_pstarPionHyp", &topo.trk_pstarPionHyp);
  tree->Branch("trk_m2", &topo.trk_m2);
  tree->Branch("trk_invbeta", &topo.trk_invbeta);
  tree->Branch("trk_nsigma_proton", &topo.trk_nsigma_proton);
  tree->Branch("trk_nsigma_kaon", &topo.trk_nsigma_kaon);
  tree->Branch("trk_nsigma_pion", &topo.trk_nsigma_pion);
  tree->Branch("trk_nsigma_electron", &topo.trk_nsigma_electron);
  tree->Branch("trk_nsigma_deutron", &topo.trk_nsigma_deutron);
  tree->Branch("trk_nsigma_triton", &topo.trk_nsigma_triton);
  tree->Branch("trk_nsigmaHtof_proton", &topo.trk_nsigmaHtof_proton);
  tree->Branch("trk_nsigmaHtof_kaon", &topo.trk_nsigmaHtof_kaon);
  tree->Branch("trk_nsigmaHtof_pion", &topo.trk_nsigmaHtof_pion);
  tree->Branch("trk_nsigmaHtof_deutron", &topo.trk_nsigmaHtof_deutron);
  tree->Branch("trk_nsigmaHtof_triton", &topo.trk_nsigmaHtof_triton);

  tree->Branch("nchRaw", &topo.nchRaw);
  tree->Branch("nchPrompt", &topo.nchPrompt);
  tree->Branch("nProton", &topo.nProton);
  tree->Branch("nPiPlus", &topo.nPiPlus);
  tree->Branch("nPiMinus", &topo.nPiMinus);
  tree->Branch("nAmbiguous", &topo.nAmbiguous);
  tree->Branch("nDisplaced", &topo.nDisplaced);
  tree->Branch("topoClass", &topo.topoClass);

  tree->Branch("nPromptProton", &topo.nPromptProton);
  tree->Branch("nPromptPiPlus", &topo.nPromptPiPlus);
  tree->Branch("nPromptPiMinus", &topo.nPromptPiMinus);
  tree->Branch("nPromptKMinus", &topo.nPromptKMinus);
  tree->Branch("nPromptElectron", &topo.nPromptElectron);
  tree->Branch("nPromptAmbiguous", &topo.nPromptAmbiguous);
  tree->Branch("nPromptUnknown", &topo.nPromptUnknown);
  tree->Branch("nDisplacedProton", &topo.nDisplacedProton);
  tree->Branch("nDisplacedPiPlus", &topo.nDisplacedPiPlus);
  tree->Branch("nDisplacedPiMinus", &topo.nDisplacedPiMinus);
  tree->Branch("nDisplacedKMinus", &topo.nDisplacedKMinus);
  tree->Branch("nDisplacedElectron", &topo.nDisplacedElectron);
  tree->Branch("nDisplacedAmbiguous", &topo.nDisplacedAmbiguous);
  tree->Branch("nDisplacedUnknown", &topo.nDisplacedUnknown);
  tree->Branch("nPromptDeuteron", &topo.nPromptDeuteron);
  tree->Branch("nPromptTriton", &topo.nPromptTriton);
  tree->Branch("nPromptHeavyAmbiguous", &topo.nPromptHeavyAmbiguous);
  tree->Branch("nDisplacedDeuteron", &topo.nDisplacedDeuteron);
  tree->Branch("nDisplacedTriton", &topo.nDisplacedTriton);
  tree->Branch("nDisplacedHeavyAmbiguous", &topo.nDisplacedHeavyAmbiguous);
  tree->Branch("nOutsideUsed", &topo.nOutsideUsed);
  tree->Branch("nBelowThreshold", &topo.nBelowThreshold);
  tree->Branch("nNoMomentum", &topo.nNoMomentum);
  tree->Branch("nPromptClosureDiff", &topo.nPromptClosureDiff);
  tree->Branch("nDisplacedClosureDiff", &topo.nDisplacedClosureDiff);

  tree->Branch("exclusiveVisibleKey", &topo.exclusiveVisibleKey);
  tree->Branch("exclusiveVisibleKeyLambdaHypothesis", &topo.exclusiveVisibleKeyLambdaHypothesis);
  tree->Branch("semiTagBits", &topo.semiTagBits);

  tree->Branch("lFlag", &topo.lFlag);
  tree->Branch("lIsSignal", &topo.lIsSignal);
  tree->Branch("lIsSideband", &topo.lIsSideband);
  tree->Branch("lambdaRegion", &topo.lambdaRegion);
  tree->Branch("nLambdaSignal", &topo.nLambdaSignal);
  tree->Branch("nLambdaSideband", &topo.nLambdaSideband);
  tree->Branch("lMass", &topo.lMass);
  tree->Branch("lMom", &topo.lMom);
  tree->Branch("lMom_x", &topo.lMom_x);
  tree->Branch("lMom_y", &topo.lMom_y);
  tree->Branch("lMom_z", &topo.lMom_z);
  tree->Branch("lDecayVtx_x", &topo.lDecayVtx_x);
  tree->Branch("lDecayVtx_y", &topo.lDecayVtx_y);
  tree->Branch("lDecayVtx_z", &topo.lDecayVtx_z);
  tree->Branch("lDecayLen", &topo.lDecayLen);
  tree->Branch("lPPiDist", &topo.lPPiDist);
  tree->Branch("lPPiAngle", &topo.lPPiAngle);
  tree->Branch("lPStar", &topo.lPStar);
  tree->Branch("lDaughterId", &topo.lDaughterId);
  tree->Branch("pLambdaDaughters", &topo.pLambdaDaughters);

  tree->Branch("nLamCand", &topo.nLamCand);
  tree->Branch("lamCand_mass", &topo.lamCand_mass);
  tree->Branch("lamCand_idPos", &topo.lamCand_idPos);
  tree->Branch("lamCand_idNeg", &topo.lamCand_idNeg);
  tree->Branch("lamCand_vtx_x", &topo.lamCand_vtx_x);
  tree->Branch("lamCand_vtx_y", &topo.lamCand_vtx_y);
  tree->Branch("lamCand_vtx_z", &topo.lamCand_vtx_z);
  tree->Branch("lamCand_dauDca", &topo.lamCand_dauDca);
  tree->Branch("lamCand_decayLen", &topo.lamCand_decayLen);
  tree->Branch("lamCand_cosPoint", &topo.lamCand_cosPoint);
  tree->Branch("lamCand_mom", &topo.lamCand_mom);
  tree->Branch("lamCand_mom_x", &topo.lamCand_mom_x);
  tree->Branch("lamCand_mom_y", &topo.lamCand_mom_y);
  tree->Branch("lamCand_mom_z", &topo.lamCand_mom_z);
  tree->Branch("lamCand_alpha", &topo.lamCand_alpha);
  tree->Branch("lamCand_qT", &topo.lamCand_qT);
  tree->Branch("lamCand_quality", &topo.lamCand_quality);
  tree->Branch("lamCand_region", &topo.lamCand_region);
  tree->Branch("lamCand_fiducial", &topo.lamCand_fiducial);
  tree->Branch("lamBestId", &topo.lamBestId);

  tree->Branch("k0Flag", &topo.k0Flag);
  tree->Branch("nK0S", &topo.nK0S);
  tree->Branch("k0Mass", &topo.k0Mass);
  tree->Branch("k0Mom", &topo.k0Mom);
  tree->Branch("k0Mom_x", &topo.k0Mom_x);
  tree->Branch("k0Mom_y", &topo.k0Mom_y);
  tree->Branch("k0Mom_z", &topo.k0Mom_z);
  tree->Branch("k0DaughterId", &topo.k0DaughterId);
  tree->Branch("pK0Daughters", &topo.pK0Daughters);

  tree->Branch("nK0Cand", &topo.nK0Cand);
  tree->Branch("k0Cand_mass", &topo.k0Cand_mass);
  tree->Branch("k0Cand_idPos", &topo.k0Cand_idPos);
  tree->Branch("k0Cand_idNeg", &topo.k0Cand_idNeg);
  tree->Branch("k0Cand_vtx_x", &topo.k0Cand_vtx_x);
  tree->Branch("k0Cand_vtx_y", &topo.k0Cand_vtx_y);
  tree->Branch("k0Cand_vtx_z", &topo.k0Cand_vtx_z);
  tree->Branch("k0Cand_dauDca", &topo.k0Cand_dauDca);
  tree->Branch("k0Cand_decayLen", &topo.k0Cand_decayLen);
  tree->Branch("k0Cand_cosPoint", &topo.k0Cand_cosPoint);
  tree->Branch("k0Cand_mom", &topo.k0Cand_mom);
  tree->Branch("k0Cand_mom_x", &topo.k0Cand_mom_x);
  tree->Branch("k0Cand_mom_y", &topo.k0Cand_mom_y);
  tree->Branch("k0Cand_mom_z", &topo.k0Cand_mom_z);
  tree->Branch("k0Cand_alpha", &topo.k0Cand_alpha);
  tree->Branch("k0Cand_qT", &topo.k0Cand_qT);
  tree->Branch("k0Cand_quality", &topo.k0Cand_quality);
  tree->Branch("k0Cand_region", &topo.k0Cand_region);
  tree->Branch("k0Cand_fiducial", &topo.k0Cand_fiducial);
  tree->Branch("k0BestId", &topo.k0BestId);
  tree->Branch("v0OverlapFlag", &topo.v0OverlapFlag);
  tree->Branch("v0OverlapTrackIds", &topo.v0OverlapTrackIds);
  tree->Branch("sidebandScaleLambda", &topo.sidebandScaleLambda);
  tree->Branch("sidebandScaleK0S", &topo.sidebandScaleK0S);
  tree->Branch("lPPiAngleValid", &topo.lPPiAngleValid);

  tree->Branch("visibleNetCharge", &topo.visibleNetCharge);
  tree->Branch("visibleScalarPSum", &topo.visibleScalarPSum);
  tree->Branch("visibleVectorPSum_x", &topo.visibleVectorPSum_x);
  tree->Branch("visibleVectorPSum_y", &topo.visibleVectorPSum_y);
  tree->Branch("visibleVectorPSum_z", &topo.visibleVectorPSum_z);
  tree->Branch("visibleVectorPSum", &topo.visibleVectorPSum);
  tree->Branch("visiblePParallelQ", &topo.visiblePParallelQ);
  tree->Branch("visiblePTransverseQ", &topo.visiblePTransverseQ);
  tree->Branch("visibleEnergyPionHyp", &topo.visibleEnergyPionHyp);
  tree->Branch("visibleEnergyProtonHyp", &topo.visibleEnergyProtonHyp);

  tree->Branch("mRecoilLambda", &topo.mRecoilLambda);
  tree->Branch("mRecoilLambdaP", &topo.mRecoilLambdaP);
  tree->Branch("mRecoilLambdaPim", &topo.mRecoilLambdaPim);
  tree->Branch("mRecoilK0S", &topo.mRecoilK0S);
  tree->Branch("mRecoilPromptPim", &topo.mRecoilPromptPim);
  tree->Branch("mRecoilVisible", &topo.mRecoilVisible);
  tree->Branch("pRecoilLambda", &topo.pRecoilLambda);
  tree->Branch("pRecoilLambdaP", &topo.pRecoilLambdaP);
  tree->Branch("pRecoilLambdaPim", &topo.pRecoilLambdaPim);
  tree->Branch("pRecoilK0S", &topo.pRecoilK0S);
  tree->Branch("pRecoilVisible", &topo.pRecoilVisible);
  tree->Branch("eRecoilLambda", &topo.eRecoilLambda);
  tree->Branch("eRecoilVisible", &topo.eRecoilVisible);

  tree->Branch("addProtonIdAll", &topo.addProtonIdAll);
  tree->Branch("mLambdaPAll", &topo.mLambdaPAll);
  tree->Branch("pLambdaPAll", &topo.pLambdaPAll);
  tree->Branch("cosThetaLambdaPAll", &topo.cosThetaLambdaPAll);
  tree->Branch("mRecoilLambdaPAll", &topo.mRecoilLambdaPAll);
  tree->Branch("cosThetaLambdaPQAll", &topo.cosThetaLambdaPQAll);
  tree->Branch("addPimIdAll", &topo.addPimIdAll);
  tree->Branch("mLambdaPimAll", &topo.mLambdaPimAll);
  tree->Branch("pLambdaPimAll", &topo.pLambdaPimAll);
  tree->Branch("cosThetaLambdaPimAll", &topo.cosThetaLambdaPimAll);
  tree->Branch("mRecoilLambdaPimAll", &topo.mRecoilLambdaPimAll);
  tree->Branch("cosThetaLambdaPimQAll", &topo.cosThetaLambdaPimQAll);
  tree->Branch("addPipIdAll", &topo.addPipIdAll);
  tree->Branch("mLambdaPipAll", &topo.mLambdaPipAll);
  tree->Branch("pLambdaPipAll", &topo.pLambdaPipAll);
  tree->Branch("cosThetaLambdaPipAll", &topo.cosThetaLambdaPipAll);

  tree->Branch("nAddProton", &topo.nAddProton);
  tree->Branch("addProtonId", &topo.addProtonId);
  tree->Branch("addProtonMom", &topo.addProtonMom);
  tree->Branch("mLambdaP", &topo.mLambdaP);
  tree->Branch("cosThetaLambdaP", &topo.cosThetaLambdaP);
  tree->Branch("pLambdaP", &topo.pLambdaP);
  tree->Branch("nAddPim", &topo.nAddPim);
  tree->Branch("addPimId", &topo.addPimId);
  tree->Branch("addPimMom", &topo.addPimMom);
  tree->Branch("mLambdaPim", &topo.mLambdaPim);
  tree->Branch("cosThetaLambdaPim", &topo.cosThetaLambdaPim);
  tree->Branch("pLambdaPim", &topo.pLambdaPim);

  tree->Branch("mFwdPPim", &topo.mFwdPPim);
  tree->Branch("lFwdFlag", &topo.lFwdFlag);
  tree->Branch("lFwdIsSignal", &topo.lFwdIsSignal);
  tree->Branch("lFwdPionId", &topo.lFwdPionId);
  tree->Branch("lFwdMass", &topo.lFwdMass);
  tree->Branch("lFwdPimMom", &topo.lFwdPimMom);
  tree->Branch("fwdLamDca", &topo.fwdLamDca);
  tree->Branch("fwdLamVtx_x", &topo.fwdLamVtx_x);
  tree->Branch("fwdLamVtx_y", &topo.fwdLamVtx_y);
  tree->Branch("fwdLamVtx_z", &topo.fwdLamVtx_z);

  tree->Branch("kmIncFlag", &topo.kmIncFlag);
  tree->Branch("escFlag", &topo.escFlag);
  tree->Branch("escKmMom", &topo.escKmMom);
  tree->Branch("escKmM2", &topo.escKmM2);

#if WriteFullDst
  treetpc = new TTree("tpc","Data Summary Table of GenfitKpTopology");    
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
  treetpc->Branch( "ntKuramaCandidate", &event.ntKuramaCandidate );  
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
  treetpc->Branch("PSTrigA", &event.PSfacTrigA);  
  treetpc->Branch("PSTrigB", &event.PSfacTrigB);
  treetpc->Branch("ScatPPhi", &event.kp_phi);
  treetpc->Branch("ScatPTheta", &event.kp_theta);    
  treetpc->Branch("InclusiveFlag", &event.IncFlag);      
  treetpc->Branch("EscapeFlag", &event.EscFlag);    
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
#endif // WriteFullDst

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
  src.ntKuramaCandidate = new TTreeReaderValue<Int_t>( *reader, "ntKuramaCandidate" );  
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

  src.lflag		  = new TTreeReaderValue<Int_t>(*reader,"Lflag");
  src.lmass		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMass");
  src.ldecayvtx_x	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_x");
  src.ldecayvtx_y	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_y");
  src.ldecayvtx_z	  = new TTreeReaderValue<Double_t>(*reader,"LambdaDecayVtx_z");
  src.lmom		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom");
  src.lmom_x		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_x");
  src.lmom_y		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_y");
  src.lmom_z		  = new TTreeReaderValue<Double_t>(*reader,"LambdaMom_z");
  src.ppi_dist            = new TTreeReaderValue<Double_t>(*reader,"LambdaVtxCloseDist");
  src.ppiangle            = new TTreeReaderValue<Double_t>(*reader,"LambdaPPiAngle");
  src.ldecays_id          = new TTreeReaderValue<std::vector<Int_t>>(*reader,"LDecaysTrackId");
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
