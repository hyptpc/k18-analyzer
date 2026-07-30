// -*- C++ -*-

/*
//Comment by Ichikawa
//TPCLocalTrackHelix.cc is for Helix fit
//Pre circle fit (TMinuit -> reduced chi2 method)

//Comment by Wooseung
Please see the discription in the TPCTrackSearch.cc
The track coordinate origin is the target center, ***NOT TPC center***

//Equation of helix (HelixParIndex in TPCLocalTrackHelix.hh)
x = -X, y = Z - tpc::Z_TARGET, z = Y;
x = p[kHelixCx] + p[kHelixR]*cos(theta);
y = p[kHelixCy] + p[kHelixR]*sin(theta);
z = p[kHelixZ0] + p[kHelixDz]*p[kHelixR]*(theta);

//Containers
1. m_hit_array : cluster(hit) container
2. m_hit_t : theta value of the cluster
3. m_hit_order : order of clusters along the track (ascending order in theta)

//FCN functions for chisqr2/ndf minimization by Minuit
1. fcn_circle
2. fcn_line
3. fcn_helix

//for beam tracks are fitted with the fixed momentum parameter.

//Fitting process
1. Preliminary fitting (Get helix initial parameters)
1-1. Circle fiting on the horizontal plane

+Optional step Line searching Hough-Transform on the vertical plane (>= 10 ms per track)
1-2. Straight-line fitting on the vertical plane

2. Helix fitting (track chisqr minimization in the 3-D)
*/

#include "TPCLocalTrackHelix.hh"

#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TMinuit.h>

#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "FuncName.hh"
#include "HoughTransform.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "TPCPadHelper.hh"
#include "UserParamMan.hh"

#include <std_ostream.hh>

#define DebugDisp 0  // 0=off, 1=on (#if DebugDisp)
#define IterativeResolution 1
#define CircCross 0

namespace
{
  const auto& gUser = UserParamMan::GetInstance();

  // B-field
  const Double_t& HS_field_0 = ConfMan::Get<Double_t>("HSFLDCALIB");
  const Double_t& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const Double_t& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

  // For Minuit / fit workspace (static globals for fcn_* callbacks)
  static Int_t gNumOfHits;
  static std::vector<TVector3> gHitPos;
  static std::vector<TVector3> gRes;
  static std::vector<Double_t> gHelixTheta;
  static std::vector<Int_t> gLayer;
  static std::vector<Double_t> gPadTheta;
  static std::vector<std::vector<Double_t>> gResParam;
  static TVector3 gVertex;
  static TVector3 gVertexRes;
  static Double_t gPar[5] = {0};
  static Double_t gChisqr = 1.e+10;
  static Bool_t gMomConstraint = false; //w or w/o momentum constraint
  static Bool_t gMultiLoop = false;
  static Int_t gBadHits;
  static Int_t gMinuitStatus;

  constexpr Double_t InvalidChi2 = 1.e+10; // CalcChi2 sentinel when ndf < 6

  // For helix fit parameter limits (cx, cy, z0, r, dz)
  const Double_t FitStep[5] = {0.1, 0.1, 0.1, 0.1, 0.0001};
  const Double_t LowLimitBeam[5] = { -50000., -50000.,  -1000.,   2000, -15. }; // about 0.6 GeV/c
  const Double_t LowLimit[5]     = { -50000., -50000., -15000.,     0., -15. };
  const Double_t UpLimit[5]      = {  50000.,  50000.,  15000., 10000.,  15. }; // 3.0 GeV/c

  // For hit selection window
  //Add hits into the track within the window.
  //(ResidualWindowPull > residual/resolution), (ResidualWindowXZ > residual_xz)
  //const Double_t ResidualWindowPullXZ = 10.; //bad
  const Double_t ResidualWindowPullXZ = 6.; //ref
  const Double_t ResidualWindowPullY = 6.; //ref

  const Double_t ResidualWindowUnderTgtXZ = 5; //[mm]
  const Double_t ResidualWindowInXZ = 10; //[mm]
  const Double_t ResidualWindowOutXZ = 10; //[mm]

  // For good hit selection, hypot(pull_t,pull_y) < PullWindow
  const Double_t PullWindow = 3;

  const Double_t ArcLengthWindow = 10; //[mm] helix arc length for local theta search window
  const Double_t ThetaNSigma = 5;

  // For CalcHelixTheta()
  const Double_t ThetaShiftK = 0.5; // sqrt(dvar/n) > ThetaShiftK * max(|2*pi*r*dz|, eps)
  const Double_t ThetaShiftPitchEps = 1.e-3; //[mm]

  // Hit array capacity (32 layers * ~10 hits/layer)
  const Int_t ReservedNumOfHits = 32 * 10;

  // For SeparateClustersWithGap()
  const Double_t MaxGapBtwClusters = 100.; //[mm]

  // For DoFit() iteration limit
  const Int_t MaxIteration = 100;
  //const Int_t MaxIteration = 500; //ref

  // For dE/dx (truncated mean in Calculate): keep the lowest 80% of per-hit dE/dx
  const Double_t TruncatedMeanRatio = 0.8;

  // For fit quality cut (DoFit / FinalizeTrack)
  const Double_t MaxChisqr = 400.;
  //const Double_t MaxChisqr = 1000.;

  // For DoFit(): MinHits override when IsBackward() && GetIsBeam()==1
  constexpr Int_t MinHitsBackwardBeam = 3;

  // HelixFit(): accept fit when Minuit fails but reduced chi2 is below this
  const Double_t GoodChisqr = 1.5;
  const Double_t GoodChisqrBeam = 4.5; // beam-track fallback (IsBeam==1)

  // For Minuit retry (HelixFit / HelixFitwVertex)
  const Int_t MaxTryMinuit = 3;
  //const Int_t MaxTryMinuit = 0;

  // For SeparateTracksAtTarget() beam/scatter boundary (mirror TPCLocalTrack)
  // Gate: (m_dz > MinSlope OR pT < MaxPt) && gap > MinGap && beam-like flips.
  // pT from helix radius at fixed B: pT[GeV/c] = 0.3 * B[T] * R[m],  R = |m_r|*1e-3.
  const Double_t MinGapForBeamScatterSep = 30.; // [mm], same as TPCLocalTrack
  const Double_t MinSlopeForBeamScatterSep = 0.01; // helix dz, cf. m_v0 in TPCLocalTrack

  // TODO: try 0.6 GeV/c on Kaon-beam (~0.65) runs
  const Double_t MaxPtForBeamScatterSep = 1.0; // [GeV/c], ~ m_r < 3.3 m at 1 T

  // Beam-like hit box upstream of target (keep in sync with TPCLocalTrack.cc)
  constexpr Double_t BeamLikeXMin = -30.; // [mm]
  constexpr Double_t BeamLikeXMax =  40.;
  constexpr Double_t BeamLikeYMin = -45.;
  constexpr Double_t BeamLikeYMax =  45.;

  // IsBackward(): |lab X| at extrapolation plane Z = -250 mm
  constexpr Double_t BackwardMaxAbsX = 75.; // [mm]

  // SetIsBeam(): helix pitch |dz| below this is treated as straight beam
  constexpr Double_t BeamLikeMaxAbsDz = 0.05;

  // CheckIsAccidental(): nominal beam momentum [GeV/c].
  // Hardcoded per build for now (e72_735 → 0.735).
  // TODO (future): gUser.GetParameter("BeamMom") and/or TPCLocalTrackHelix::SetBeamMom()
  //   from Dst BeginRun (run-dependent momentum).
  constexpr Double_t BeamMom = 0.735;
  // constexpr Double_t BeamMom = 0.933;

  // p_t band lower edge: tag accidental only if p_t >= BeamMom - BeamMomOffset [GeV/c]
  constexpr Double_t BeamMomOffset = 0.2;

  // min helix start/end z difference [mm]; TPC z span ~500 mm so just half of it
  constexpr Double_t AccidentalMinDiffZ = 250.;

  //Horizontal resolution function
  //x : alpha(track-pad angle), y : y pos of cluster (y+300 : Drift length)
  //[0] : Intrinsic XZ resolution, [1] : Attenuation term, [2] : Diffusion coefficient, [3] : Effective # of signal electrons, [4] : Pad length, [5] : Effective # of electron clusters
  static TString eq_horizontal =
    "TMath::Sqrt("
    "TMath::Power([0],2.) + " // Intrinsic XZ resolution
    "TMath::Power([2],2.)*(y+300.)/([3]*TMath::Exp(-[1]*(y+300.))) + " // Attenuation term
    "TMath::Power([4]*TMath::Tan(x),2.)/(12.*[5])" // angular term
    ")";
  static TF2 *f_horizontal = new TF2("f_horizontal", eq_horizontal.Data(), -4., 4., -300., 300.);

  //Vertical resolution function
  //x : x pos of cluster (x+300 : Drift length)
  //[0] : Intrinsic Y resolution, [1] : Attenuation term, [2] : Diffusion coefficient, [3] : Effective # of signal electrons
  static TString eq_vertical =
    "TMath::Sqrt("
    "TMath::Power([0],2.) + " // Intrinsic Y resolution
    "TMath::Power([2],2.)*(x+300.)/([3]*TMath::Exp(-[1]*(x+300.)))" // Attenuation term
    ")";
  static TF1 *f_drift = new TF1("f_drift", eq_vertical.Data(), -300., 300.);

  // --- TF1: scan helix parameter theta (x) to minimize distance to a point ---
  // Default theta range (overridden by SetRange in EvalTheta*)
  const Double_t HELIX_THETA_HALF_W = 10. * TMath::Pi();

  // 3D squared distance (target - point on helix).
  // TF1 [0..4] = helix cx,cy,z0,r,dz; [5][6][7] = target x,y,z
  static TF1 fint(
      "fint",
      "TMath::Power([5] - ([0] + [3]*TMath::Cos(x)), 2.) + "
      "TMath::Power([6] - ([1] + [3]*TMath::Sin(x)), 2.) + "
      "TMath::Power([7] - ([2] + [3]*[4]*x), 2.)",
      -HELIX_THETA_HALF_W, HELIX_THETA_HALF_W);

  // XY projection: squared circle–point distance. 
  // [0] cx, [1] cy, [2] target x, [3] R, [4] target y
  static TF1 fintXZ(
      "fintXZ",
      "TMath::Power([2] - ([0] + [3]*TMath::Cos(x)), 2.) + "
      "TMath::Power([4] - ([1] + [3]*TMath::Sin(x)), 2.)",
      -HELIX_THETA_HALF_W, HELIX_THETA_HALF_W);

  // Squared concentric-circle crossing term (theta refinement). 
  // [0] cx, [1] cy, [2] R, [3] target radius hypot(x,y)
  // Inner: cx^2 + cy^2 + R^2 + 2*R*(cx*cos(x)+cy*sin(x)) - rho^2
  static TF1 fcir_cross(
      "fcir_cross",
      "TMath::Power("
      "[0]*[0] + [1]*[1] + [2]*[2] + "
      "2.0*[2]*([0]*TMath::Cos(x) + [1]*TMath::Sin(x)) - [3]*[3], 2.)",
      -HELIX_THETA_HALF_W, HELIX_THETA_HALF_W);
  
  //______________________________________________________________________________
  // SetNpx scales with |Δθ|; HELIX_THETA_SCAN_STEP is the θ step in radians.
  static constexpr Double_t HELIX_THETA_SCAN_STEP = 1. * TMath::DegToRad();
  static constexpr Int_t HELIX_THETA_NPX_MIN = 64;
  static constexpr Int_t HELIX_THETA_NPX_MAX = 4096;
  static inline Int_t HelixThetaScanNpx(Double_t window_low, Double_t window_up)
  {
    const Double_t ref_span = 2. * HELIX_THETA_HALF_W;
    Double_t span = window_up - window_low;
    if (span <= 0. || !TMath::Finite(span)) span = ref_span;
    Int_t n = TMath::Nint(span / HELIX_THETA_SCAN_STEP);
    if (n > HELIX_THETA_NPX_MAX) n = HELIX_THETA_NPX_MAX;
    if (n < HELIX_THETA_NPX_MIN) n = HELIX_THETA_NPX_MIN;
    return n;
  }

  //______________________________________________________________________________
  static inline void NormalizeThetaWindow(Double_t& window_low, Double_t& window_up)
  {
    if (!TMath::Finite(window_low) || !TMath::Finite(window_up)) {
      window_low = -HELIX_THETA_HALF_W;
      window_up = +HELIX_THETA_HALF_W;
      return;
    }
    if (window_low > window_up) std::swap(window_low, window_up);
    if (MathTools::Equal(window_low, window_up)) {
      const Double_t center = 0.5 * (window_low + window_up);
      window_low = center - 1.e-3;
      window_up  = center + 1.e-3;
    }
  }

  //______________________________________________________________________________
  // local theta search window [rad] = ArcLengthWindow / r; cap r at ~1.5 GeV/c, 1 T assumed
  static inline Double_t HelixThetaSearchWindow(Double_t r)
  {
    const Double_t min_pt = 1.5; // [GeV/c]
    return ArcLengthWindow / TMath::Min(r, min_pt / 0.3 * 1000.);
  }

  //______________________________________________________________________________
  Double_t HelixRToPt(Double_t helix_r) // helix_r [mm] -> p_t [GeV/c]
  {
    const Double_t B = HS_field_0 * (HS_field_Hall / HS_field_Hall_calc); // [T]
    return TMath::Abs(helix_r) * (tpc::C_LIGHT * B) * 0.001; // [GeV/c]
  }

  //______________________________________________________________________________
  void LogFatalPrecondition(const char* caller, const char* message)
  {
    std::cout << caller << " Fatal error : " << message << std::endl;
  }

  //______________________________________________________________________________
  Bool_t ValidateCalcHelixTheta(const TPCLocalTrackHelix& trk, const char* caller)
  {
    if(trk.IsThetaCalculated()) return true;
    LogFatalPrecondition(caller,
      "No helix theta information!!! CalcHelixTheta() should be run in front of this");
    return false;
  }

  //______________________________________________________________________________
  void WarnValidateCalcHelixTheta(const TPCLocalTrackHelix& trk, const char* caller)
  {
    if(!trk.IsThetaCalculated())
      LogFatalPrecondition(caller,
        "No helix theta information!!! CalcHelixTheta() should be run in front of this");
  }
  
  //______________________________________________________________________________
  Bool_t ValidateTrackHasHits(const TPCLocalTrackHelix& trk, const char* caller)
  {
    if(trk.GetNHit() > 0) return true;
    LogFatalPrecondition(caller, "Empty track!");
    return false;
  }

  //______________________________________________________________________________
  // One field: actual must equal expected (logs and returns false on mismatch).
  Bool_t ValidateEqualSize(const char* caller, std::size_t expected,
                           const char* name, std::size_t actual)
  {
    if(expected == actual) return true;
    hddaq::cerr << caller << " size mismatch: expected=" << expected
                << " " << name << "=" << actual << std::endl;
    return false;
  }

  // Recursion stop: no more "name", size pairs after caller and expected.
  Bool_t ValidateEqualSizes(const char* caller, std::size_t expected)
  {
    return true;
  }

  // caller, expected are fixed; then "name", size, "name", size, ... (remaining...).
  template<typename... Remaining>
  Bool_t ValidateEqualSizes(const char* caller, std::size_t expected,
                            const char* name, std::size_t actual,
                            Remaining... remaining)
  {
    if(!ValidateEqualSize(caller, expected, name, actual)) return false;
    return ValidateEqualSizes(caller, expected, remaining...);
  }

  //______________________________________________________________________________
  // m_hit_order[k] = m_hit_array index of the k-th hit along track (permutation)
  void RemapHitOrderAfterErase(std::vector<Int_t>& hit_order, Int_t delete_hit)
  {
    std::vector<Int_t> new_order;
    new_order.reserve(hit_order.size());
    for(std::size_t k = 0; k < hit_order.size(); ++k){
      const Int_t idx = hit_order[k];
      if(idx == delete_hit) continue;
      new_order.push_back(idx > delete_hit ? idx - 1 : idx);
    }
    hit_order.swap(new_order);
  }

  //______________________________________________________________________________
  static void PrintHelixFitAbnormalWarning(
    const TString& label,
    Double_t chi2, Int_t nhit,
    Int_t itry = -1)
  {
    std::cout << "#W [" << label << "] abnormal chi2" << std::endl
              << "   chi2 = " << chi2 << std::endl;
    if(itry >= 0)
      std::cout << "   itry = " << itry << std::endl;
    std::cout << "   nhit = " << nhit << std::endl;
  }

#if DebugDisp
  inline void DebugHelixPar(const TString& comment,
                            Double_t cx, Double_t cy, Double_t z0,
                            Double_t r, Double_t dz)
  {
    std::cout << comment << std::endl
              << "    cx: " << cx << std::endl
              << "    cy: " << cy << std::endl
              << "    z0: " << z0 << std::endl
              << "     r: " << r << std::endl
              << "    dz: " << dz << std::endl;
  }

  //______________________________________________________________________________
  inline void DebugHelixPar(const TString& comment, const Double_t par[5])
  {
    DebugHelixPar(comment, par[kHelixCx], par[kHelixCy], par[kHelixZ0], par[kHelixR], par[kHelixDz]);
  }

  //______________________________________________________________________________
  // Circle fit updates cx, cy, r (= par[kHelixR]) only; debug print
  inline void DebugHelixPar(const TString& comment,
                            Double_t cx, Double_t cy, Double_t r)
  {
    std::cout << comment << std::endl
              << "    cx: " << cx << std::endl
              << "    cy: " << cy << std::endl
              << "     r: " << r << std::endl;
  }
#endif
}

//______________________________________________________________________________
static inline Bool_t CompareY(const Int_t a, const Int_t b){
  return gHitPos[a].y() < gHitPos[b].y();
}

//______________________________________________________________________________
static inline Bool_t CompareTheta(const Int_t a, const Int_t b){
  return gHelixTheta[a] < gHelixTheta[b];
}

//______________________________________________________________________________
static inline TVector3 GlobalToLocal(TVector3 pos){
  return TVector3(-pos.X(), pos.Z() - tpc::Z_TARGET, pos.Y());
}

//______________________________________________________________________________
static inline TVector3 LocalToGlobal(TVector3 pos){
  return TVector3(-pos.X(), pos.Z(), pos.Y() + tpc::Z_TARGET);
}

//______________________________________________________________________________
static inline TVector3 LocalPosition(const Double_t par[5], Double_t t){

  //TPC local coordinate
  //This is the eqation of Helix
  Double_t x = par[kHelixCx] + par[kHelixR]*cos(t);
  Double_t y = par[kHelixCy] + par[kHelixR]*sin(t);
  Double_t z = par[kHelixZ0] + (par[kHelixDz]*par[kHelixR]*t);
  return TVector3(x, y, z);
}

//______________________________________________________________________________
static inline TVector3 GlobalPosition(const Double_t par[5], Double_t t){

  TVector3 pos = LocalPosition(par, t);
  return LocalToGlobal(pos);
}

//______________________________________________________________________________
static inline Double_t EvalTheta(Double_t par[5], TVector3 pos, Double_t window_low, Double_t window_up)
{
  NormalizeThetaWindow(window_low, window_up);

  fint.SetRange(window_low, window_up);
  Double_t fpar[8];
  TVector3 localpos = GlobalToLocal(pos);
  // fint TF1 [0..4] = helix cx,cy,z0,r,dz (same order as HelixParIndex)
  fpar[0] = par[kHelixCx];
  fpar[1] = par[kHelixCy];
  fpar[2] = par[kHelixZ0];
  fpar[3] = par[kHelixR];
  fpar[4] = par[kHelixDz];
  fpar[5] = localpos.X();
  fpar[6] = localpos.Y();
  fpar[7] = localpos.Z();

  Int_t steps = HelixThetaScanNpx(window_low, window_up);
  fint.SetParameters(fpar);
  fint.SetNpx(steps);
  const Double_t theta_fint = fint.GetMinimumX();

#if CircCross
  fcir_cross.SetRange(window_low,window_up);
  Double_t cpar[4];
  cpar[0] = par[kHelixCx];
  cpar[1] = par[kHelixCy];
  cpar[2] = par[kHelixR];
  cpar[3] = TMath::Hypot(localpos.X(), localpos.Y());
  fcir_cross.SetParameters(cpar);
  fcir_cross.SetNpx(steps);
  const Double_t theta_fcir = fcir_cross.GetMinimumX();
  if (MathTools::Equal(theta_fcir, window_low) ||
      MathTools::Equal(theta_fcir, window_up)) {
    return theta_fint;
  }
  return theta_fcir;
#else
  return theta_fint;
#endif
}

//______________________________________________________________________________
static inline Double_t EvalThetaXZ(Double_t par[5], TVector3 pos, Double_t window_low, Double_t window_up){
  NormalizeThetaWindow(window_low, window_up);

  Double_t fpar[5];
  TVector3 localpos = GlobalToLocal(pos);
  // fintXZ TF1: [0,1,3]=helix cx,cy,r; [2,4]=hit x,y (not HelixParIndex)
  fpar[0] = par[kHelixCx];
  fpar[1] = par[kHelixCy];
  fpar[2] = localpos.X();
  fpar[3] = par[kHelixR];
  fpar[4] = localpos.Y();

  Int_t steps = HelixThetaScanNpx(window_low, window_up);
  fintXZ.SetRange(window_low, window_up);
  fintXZ.SetParameters(fpar);
  fintXZ.SetNpx(steps);
  Double_t min_t = fintXZ.GetMinimumX();
  return min_t;
}
//______________________________________________________________________________
static inline TVector3 ResidualVect(Double_t par[5], TVector3 pos, Double_t theta){ //Closest distance on

  // for Helix tracking (par[kHelixCx..kHelixDz])
  Double_t window = HelixThetaSearchWindow(par[kHelixR]);
  Double_t evaltheta = EvalTheta(par, pos, theta - 0.5*window, theta + 0.5*window);
  TVector3 calpos = GlobalPosition(par, evaltheta);
  TVector3 d = pos - calpos;
  return d;
}

//______________________________________________________________________________
static inline TVector3 ResidualVect(Double_t par[5], TVector3 pos, Double_t theta_min, Double_t theta_max){ //Closest distance on

  // for Helix tracking (par[kHelixCx..kHelixDz])
  Double_t theta = EvalTheta(par, pos, theta_min, theta_max);
  TVector3 calpos = GlobalPosition(par, theta);
  TVector3 d = pos - calpos;
  return d;
}

//______________________________________________________________________________
static inline TVector3 ResidualVectXZ(Double_t par[5], TVector3 pos){ //Closest distance on the y=pos.y() plane

  TVector3 localpos = GlobalToLocal(pos);
  Double_t x = localpos.x(); Double_t y = localpos.y();
  Double_t magnitude = TMath::Hypot(x - par[kHelixCx], y - par[kHelixCy]) - par[kHelixR];
  TVector3 direction(x - par[kHelixCx], y - par[kHelixCy], 0.);
  TVector3 resi_local = magnitude*direction.Unit();
  TVector3 resi_global(-resi_local.x(), 0., resi_local.y());
  return resi_global;
}

//______________________________________________________________________________
static inline TVector3 CalcResolution(Double_t par[5], Int_t layer, TVector3 pos, Double_t padTheta, Double_t theta, std::vector<Double_t> resparam, Bool_t vetoBadClusters){

  Double_t cos_pad    = TMath::Cos(padTheta);
  Double_t sin_pad    = TMath::Sin(padTheta);
  Double_t tan_pad    = TMath::Tan(padTheta);
  Double_t pad_length = tpc::padParameter[layer][tpc::kLength];
  Double_t pad_radius = tpc::padParameter[layer][tpc::kRadius];
  TVector3 closest_dist_to_track_xz = ResidualVectXZ(par, TVector3(0., 0., tpc::Z_TARGET));

  // check whether the track is crossing the layer or not
  if (vetoBadClusters &&
      TMath::Abs(closest_dist_to_track_xz.Mag() - pad_radius) < 0.5 * pad_length)
    return TVector3(1.e+10, 1.e+10, 1.e+10);

  //alpha : pad - track angle
  TVector3 localpos = GlobalToLocal(pos);
  Double_t tan_track = (localpos.y()-par[kHelixCy])/(localpos.x()-par[kHelixCx]);
  Double_t tan_diff  = (tan_pad-tan_track)/(1.+tan_pad*tan_track);
  Double_t alpha     = TMath::ATan(tan_diff);

  //Calculate resolution
  //horizontal resolution
  Double_t param_horizontal[6] = {resparam[0], resparam[1], resparam[2], resparam[3], resparam[4], resparam[5]};
  f_horizontal->SetParameters(param_horizontal);
  Double_t res_horizontal = f_horizontal->Eval(alpha, pos.y());

  //vertical resolution
  Double_t param_y[4] = {resparam[6], resparam[1], resparam[7], resparam[8]};
  f_drift->SetParameters(param_y);
  Double_t res_drift = f_drift->Eval(pos.y());

  TVector3 res(res_horizontal*TMath::Abs(cos_pad), res_drift, res_horizontal*TMath::Abs(sin_pad));

  //Residual/resolution < window
  if(vetoBadClusters){
    TVector3 resi = ResidualVect(par, pos, theta);
    Double_t residual_vertical   = resi.y();
    Double_t residual_horizontal = TMath::Hypot(resi.x(), resi.z());
    Double_t pull_t = residual_horizontal / res_horizontal;
    Double_t pull_y = residual_vertical / res_drift;
    if (TMath::Hypot(pull_t, pull_y) > PullWindow) return TVector3(2.e+10, 2.e+10, 2.e+10);
  }

  return res;
}

//______________________________________________________________________________
static inline void fcn_helix(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t chisqr=0.; Int_t dof = 0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(par, gHitPos[i], gHelixTheta[i]);
    if (tpc::IsDummyResolutionVec(gRes[i])) continue; // exclude dummy hits in calculation
    chisqr += TMath::Sq( TMath::Hypot(d.x(), d.z()) / TMath::Hypot(gRes[i].x(), gRes[i].z()) );
    dof++;
    chisqr += TMath::Sq( d.y() / gRes[i].y() );
    dof++;
  }
  if(gMomConstraint) dof += 1; //if there is a momentum constraint
  f = chisqr/(Double_t)(dof - 5);
}

//______________________________________________________________________________
static inline void fcn_helixwVertex(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t chisqr=0.; Int_t dof = 0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(par, gHitPos[i], gHelixTheta[i]);
    if (tpc::IsDummyResolutionVec(gRes[i])) continue; // exclude dummy hits in calculation
    chisqr += TMath::Sq( TMath::Hypot(d.x(), d.z()) / TMath::Hypot(gRes[i].x(), gRes[i].z()) );
    dof++;
    chisqr += TMath::Sq( d.y() / gRes[i].y() );
    dof++;
  }
  TVector3 d = ResidualVect(par, gVertex, -2.*TMath::Pi(), 2.*TMath::Pi());
  TVector3 res = gVertexRes;
  chisqr += TMath::Sq( TMath::Hypot(d.x(), d.z()) / TMath::Hypot(res.x(), res.z()) );
  dof++;
  chisqr += TMath::Sq( d.y() / res.y() );
  dof++;

  f = chisqr/(Double_t)(dof - 5);
}

//______________________________________________________________________________
static inline void fcn_line(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // par[0]=z0, par[1]=dz (2-param line fit; not HelixParIndex)
  Double_t chisqr = 0.;
  Int_t dof = 0;

  Bool_t flipcheck = false;
  Double_t prev_theta = 0.; Double_t theta0 = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 pos = GlobalToLocal(gHitPos[i]);
    //Check ATan2 function's theta flip (-pi ~ pi)
    Double_t tmp_theta = TMath::ATan2(pos.Y() - gPar[kHelixCy], pos.X() - gPar[kHelixCx]);
    if(i==0) theta0 = tmp_theta;
    if(TMath::Abs(prev_theta - tmp_theta) > TMath::Pi()) flipcheck = true;
    if(flipcheck){
      if(theta0>0. && tmp_theta<0.) tmp_theta += 2.*TMath::Pi();
      if(theta0<0. && tmp_theta>0.) tmp_theta -= 2.*TMath::Pi();
    }
    prev_theta = tmp_theta;
    const Double_t diff_z = pos.Z() - (par[0] + par[1]*gPar[kHelixR]*tmp_theta);
    chisqr += TMath::Sq(diff_z / gRes[i].y());
    dof++;
  }
  if(gMomConstraint) dof += 1; //if there is a momentum constraint
  f = chisqr/(Double_t)(dof - 2);
}

//______________________________________________________________________________
static inline Double_t CalcChi2(Double_t *HelixPar, Int_t &ndf, Bool_t vetoBadClusters)
{

  if(!ValidateEqualSizes("TPCLocalTrackHelix::CalcChi2",
                        static_cast<std::size_t>(gNumOfHits),
                        "gHitPos", gHitPos.size(),
                        "gHelixTheta", gHelixTheta.size(),
                        "gLayer", gLayer.size(),
                        "gPadTheta", gPadTheta.size(),
                        "gResParam", gResParam.size()))
    return TMath::QuietNaN();

  ndf = 0; Double_t chisqr = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(HelixPar, gHitPos[i], gHelixTheta[i]);
    TVector3 res = CalcResolution(HelixPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    if (tpc::IsDummyResolutionVec(res)) continue; // exclude bad clusters
    chisqr += TMath::Sq( TMath::Hypot(d.x(), d.z()) / TMath::Hypot(res.x(), res.z()) );
    ndf++;
    chisqr += TMath::Sq( d.y()/res.y() );
    ndf++;
  }
  if(gMomConstraint) ndf += 1; //if there is a momentum constraint
  if(ndf < 6) return 1.e+10;
  return chisqr/(Double_t)(ndf-5);
}

// gBadHits: ndf = 2 * (good hits) + ndfExtra (mom +1, vertex +2)
static inline Int_t NBadFromNdf(Int_t nHits, Int_t ndf, Int_t ndfExtra)
{
  return nHits - (ndf - ndfExtra) / 2;
}

//______________________________________________________________________________
// Maybe Taubin algebraic circle fit (centroid frame); returns sum of squared radial residuals.
static inline Double_t CircleFit(const Double_t *mX, const Double_t *mY, const Int_t npoints,
                                 Double_t* mXCenter, Double_t* mYCenter, Double_t* mRadius)
{
  if(npoints < 4){
#if DebugDisp
    hddaq::cerr << "#W CircleFit: npoints=" << npoints
                << " (algebraic fit needs >= 4 points)" << std::endl;
#endif
    return -1.;
  }
  if(npoints > ReservedNumOfHits){
#if DebugDisp
    hddaq::cerr << "#W CircleFit: npoints=" << npoints
                << " > ReservedNumOfHits=" << ReservedNumOfHits << std::endl;
#endif
    return -1.;
  }

  const Double_t n = static_cast<Double_t>(npoints);

  Double_t x_gravity = 0., y_gravity = 0.;
  for(Int_t i = 0; i < npoints; ++i){
    x_gravity += mX[i];
    y_gravity += mY[i];
  }
  x_gravity /= n;
  y_gravity /= n;

  Double_t sum_xx = 0., sum_yy = 0., sum_xy = 0.;
  Double_t sum_xz = 0., sum_yz = 0., sum_zz = 0.;
  for(Int_t i = 0; i < npoints; ++i){
    const Double_t x = mX[i] - x_gravity;
    const Double_t y = mY[i] - y_gravity;
    const Double_t z = x*x + y*y;
    sum_xx += x*x;
    sum_yy += y*y;
    sum_xy += x*y;
    sum_xz += x*z;
    sum_yz += y*z;
    sum_zz += z*z;
  }

  if(TMath::Abs(sum_xx) < 0.0001 || TMath::Abs(sum_yy) < 0.0001 || TMath::Abs(sum_xy) < 0.0001){
    hddaq::cerr << "#W CircleFit: x2=" << sum_xx << " y2=" << sum_yy << " xy=" << sum_xy
                << " grav=(" << x_gravity << ", " << y_gravity << ")" << std::endl;
    return -1.;
  }

  // Centroid moments.
  const Double_t f = (3.*sum_xx + sum_yy) / n;
  const Double_t g = (sum_xx + 3.*sum_yy) / n;
  const Double_t h = 2.*sum_xy / n;
  const Double_t p = sum_xz / n;
  const Double_t q = sum_yz / n;
  const Double_t t = sum_zz / n;
  const Double_t g0 = (sum_xx + sum_yy) / n;
  const Double_t g0_sq = g0*g0;
  const Double_t g0_qu = g0_sq*g0_sq;

  // Taubin: solve P(xroot)=0 with xroot = 1 + 2*lambda/g0 (lambda = algebraic eigenvalue; xroot ~ 1).
  const Double_t a = -4.;
  const Double_t b = (f*g - t - h*h) / g0_sq;
  const Double_t c = (t*(f + g) - 2.*(p*p + q*q)) / (g0_sq*g0);
  const Double_t d = (t*(h*h - f*g) + 2.*(p*p*g + q*q*f) - 4.*p*q*h) / g0_qu;

  Double_t xroot = 1.;
  for(Int_t iter = 0; iter < 5; ++iter){
    const Double_t P = (((xroot + a)*xroot + b)*xroot + c)*xroot + d;
    const Double_t dP = ((4.*xroot + 3.*a)*xroot + 2.*b)*xroot + c;
    xroot -= P/dP;
  }

  const Double_t g1 = xroot*g0;
  const Double_t xnom1 = (g - g1)*(f - g1) - h*h;
  if(TMath::Abs(xnom1) < 0.0001 || TMath::IsNaN(xnom1)){
    hddaq::cerr << "#W CircleFit: xnom1=" << xnom1 << std::endl;
    return -1.;
  }

  const Double_t yd = (q*(f - g1) - h*p) / xnom1;
  const Double_t xnom2 = f - g1;
  if(TMath::Abs(xnom2) < 0.0001 || TMath::IsNaN(xnom2)){
    hddaq::cerr << "#W CircleFit: xnom2=" << xnom2 << std::endl;
    return -1.;
  }

  const Double_t xd = (p - h*yd) / xnom2;
  const Double_t radius2 = xd*xd + yd*yd + g1;
  *mXCenter = xd + x_gravity;
  *mYCenter = yd + y_gravity;
  *mRadius = TMath::Sqrt(radius2);

  Double_t geom_chi2 = 0.;
  for(Int_t i = 0; i < npoints; ++i){
    const Double_t dx = mX[i] - (*mXCenter);
    const Double_t dy = mY[i] - (*mYCenter);
    const Double_t rho2 = dx*dx + dy*dy;
    const Double_t dr = TMath::Sqrt(rho2) - (*mRadius);
    geom_chi2 += dr*dr;
  }

  if(geom_chi2 < 0.){
    hddaq::cerr << "#W CircleFit: variance=" << geom_chi2 << std::endl;
    return -1.;
  }

#if DebugDisp
  std::cout << "CircleFit fitting :"
            << " variance: " << geom_chi2
            << ", radius: " << *mRadius
            << std::endl;
#endif

  return geom_chi2;
}

//______________________________________________________________________________
static inline Bool_t StraightLineFit()
{

  if(!ValidateEqualSizes("TPCLocalTrackHelix::StraightLineFit",
                        static_cast<std::size_t>(gNumOfHits),
                        "gHitPos", gHitPos.size(),
                        "gHelixTheta", gHelixTheta.size(),
                        "gLayer", gLayer.size(),
                        "gPadTheta", gPadTheta.size(),
                        "gResParam", gResParam.size()))
    return false;

  //pre t-y fit
  // par_li[2]: Minuit line fit (z0, dz) = helix par[kHelixZ0, kHelixDz]
  Double_t par_li[2] = {gPar[kHelixZ0], gPar[kHelixDz]};
  Double_t err_li[2] = {-999., -999.};
  TMinuit *minuit = new TMinuit(2);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_line);

  Int_t ierflg_li = 0;
  Double_t arglist_li[10];
  arglist_li[0] = 2.3;
  minuit->mnexcm("SET ERR", arglist_li,1,ierflg_li); // No warnings
  arglist_li[0] = 1;
  minuit->mnexcm("SET NOW", arglist_li,1,ierflg_li); // No warnings

  TString name_li[2] = {"z0", "dz"};
  minuit->mnparm(0, name_li[0], par_li[0], FitStep[kHelixZ0], LowLimit[kHelixZ0], UpLimit[kHelixZ0], ierflg_li);
  minuit->mnparm(1, name_li[1], par_li[1], FitStep[kHelixDz], LowLimit[kHelixDz], UpLimit[kHelixDz], ierflg_li);

  minuit->Command("SET STRategy 0");
  arglist_li[0] = 1000*5*5*5;
  arglist_li[1] = 0.1/(10.*10.*10.);
  minuit->mnexcm("MIGRAD", arglist_li, 2, ierflg_li);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  Int_t Err_li;
  Double_t bnd1_li, bnd2_li;
  for(Int_t i=0; i<2; i++){
    minuit->mnpout(i, name_li[i], par_li[i], err_li[i], bnd1_li, bnd2_li, Err_li);
  }

  //Double_t grad[5]; Double_t Chisqr;
  //minuit->Eval(5, grad, Chisqr, par_li, 0);
  delete minuit;
#if DebugDisp
  if(icstat==0) std::cout<<"StraightLineFit() icstat=="<<icstat<<std::endl;
#endif

  gPar[kHelixZ0] = par_li[0];
  gPar[kHelixDz] = par_li[1];
  gMinuitStatus = icstat;
  return true;
}

//______________________________________________________________________________
static inline Bool_t HelixFit(Int_t IsBeam, Bool_t vetoBadClusters, Bool_t ExclusiveFlag = false){
  // Loose fail-safe: reject only pathological fit states.
  constexpr Double_t max_allowed_chi2_helix = 1.e8;

  if(!ValidateEqualSizes("TPCLocalTrackHelix::HelixFit",
                        static_cast<std::size_t>(gNumOfHits),
                        "gHitPos", gHitPos.size(),
                        "gHelixTheta", gHelixTheta.size(),
                        "gLayer", gLayer.size(),
                        "gPadTheta", gPadTheta.size(),
                        "gResParam", gResParam.size()))
    return false;

  Double_t par[5] = {gPar[kHelixCx], gPar[kHelixCy], gPar[kHelixZ0], gPar[kHelixR], gPar[kHelixDz]};
  Double_t err[5] = {-999., -999., -999., -999., -999.};

  gRes.clear();
  for(Int_t i=0; i<gNumOfHits; i++){
    TVector3 res = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    gRes.push_back(res);
  }

  // NDF value becomes different. Initialize gChisqr and fit again.
  Int_t ndf;
  if(vetoBadClusters){
    gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);
    gBadHits = NBadFromNdf(gNumOfHits, ndf, gMomConstraint ? 1 : 0);
    if(!TMath::Finite(gChisqr) || gChisqr >= InvalidChi2) return false;
  }

  //for exclusive tracking, previous fitting result should not affect current fitting
  if (ExclusiveFlag) gChisqr = 1.e+10;

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_helix);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist, 1, ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(Int_t i=0; i<5; i++){
    if(IsBeam==1) minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimitBeam[i], UpLimit[i], ierflg);
    else minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  if(gMomConstraint) minuit -> FixParameter(3);

  minuit->Command("SET STRategy 0");
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;

  Bool_t status = false;
  Bool_t hit_abnormal = false;
  Int_t itry=0; gMinuitStatus = 0;
  while(!status){
    if(itry>MaxTryMinuit) break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

    //minuit->mnprin(4, amin);
    for(Int_t i=0; i<5; i++){
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    }

    //Double_t grad[5];
    //minuit -> Eval(5, grad, Chisqr, par, 0);
    Double_t Chisqr = CalcChi2(par, ndf, vetoBadClusters);
    if(!TMath::Finite(Chisqr) || Chisqr > max_allowed_chi2_helix){
      PrintHelixFitAbnormalWarning("HelixFit", Chisqr, gNumOfHits, itry);
      hit_abnormal = true;
      status = false;
      break;
    }
    if(gChisqr>=Chisqr || TMath::Abs(gChisqr-Chisqr) < 0.01){
      gChisqr = Chisqr;
      gPar[kHelixCx] = par[kHelixCx];
      gPar[kHelixCy] = par[kHelixCy];
      gPar[kHelixZ0] = par[kHelixZ0];
      gPar[kHelixR] = par[kHelixR];
      gPar[kHelixDz] = par[kHelixDz];
      gMinuitStatus = icstat;
      gBadHits = NBadFromNdf(gNumOfHits, ndf, gMomConstraint ? 1 : 0);
      status = true;
    }
    arglist[0] = arglist[0]*5;
    arglist[1] = arglist[1]*0.1;
    ++itry;
  }
  delete minuit;
  if(hit_abnormal) return false;
  if(IsBeam==1 && !status && gChisqr < GoodChisqrBeam) status = true;
  if(gChisqr < GoodChisqr) status = true;

#if DebugDisp
  std::cout<<"HelixFit() status="<<status<<" gChisqr "<<gChisqr<<" isBeam "<<IsBeam<<std::endl;
  if(gMinuitStatus==0) std::cout<<"HelixFit() icstat==0"<<std::endl;
#endif
  return status;
}

//______________________________________________________________________________
static inline Bool_t HelixFitInvertCharge(){
  // Fail fast on pathological fit states to avoid MIGRAD hanging.
  constexpr Double_t max_allowed_chi2 = 1.e6;

  if(!ValidateEqualSizes("TPCLocalTrackHelix::HelixFitInvertCharge",
                        static_cast<std::size_t>(gNumOfHits),
                        "gHitPos", gHitPos.size(),
                        "gHelixTheta", gHelixTheta.size(),
                        "gLayer", gLayer.size(),
                        "gPadTheta", gPadTheta.size(),
                        "gResParam", gResParam.size(),
                        "gRes", gRes.size()))
    return false;

  Double_t par[5] = {gPar[kHelixCx], gPar[kHelixCy], gPar[kHelixZ0], gPar[kHelixR], gPar[kHelixDz]};
  Double_t err[5] = {-999., -999., -999., -999., -999.};

  Double_t lowLimit[5];
  lowLimit[0] = gPar[kHelixCx] - 0.5*gPar[kHelixR];
  lowLimit[1] = gPar[kHelixCy] - 0.5*gPar[kHelixR];
  lowLimit[2] = gPar[kHelixZ0] - 5000.;
  lowLimit[3] = 0.;
  lowLimit[4] = gPar[kHelixDz] - 5.;
  Double_t upLimit[5];
  upLimit[0] = gPar[kHelixCx] + 0.5*gPar[kHelixR];
  upLimit[1] = gPar[kHelixCy] + 0.5*gPar[kHelixR];
  upLimit[2] = gPar[kHelixZ0] + 5000.;
  upLimit[3] = 10000.; // 3 GeV/c at 1 T
  upLimit[4] = gPar[kHelixDz] + 5.;

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_helix);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist, 1, ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], lowLimit[i], upLimit[i], ierflg);
  }

  minuit->Command("SET STRategy 0");
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;
  Int_t ndf;
  gMinuitStatus = 0;

  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  minuit->mnimpr();
  //minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  //minuit->mnprin(4, amin);
  for(Int_t i=0; i<5; i++){
    minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
  }

  gPar[kHelixCx] = par[kHelixCx];
  gPar[kHelixCy] = par[kHelixCy];
  gPar[kHelixZ0] = par[kHelixZ0];
  gPar[kHelixR]  = par[kHelixR];
  gPar[kHelixDz] = par[kHelixDz];
  gChisqr = CalcChi2(gPar, ndf, false);
  if(!TMath::Finite(gChisqr) || gChisqr > max_allowed_chi2){
    PrintHelixFitAbnormalWarning("HelixFitInvertCharge", gChisqr, gNumOfHits);
    delete minuit;
    return false;
  }

  Int_t itry=0;
  Bool_t status = false;
  while(!status){
    if(itry>MaxTryMinuit) break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

    //minuit->mnprin(4, amin);
    for(Int_t i=0; i<5; i++){
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    }

    //Double_t grad[5];
    //minuit -> Eval(5, grad, Chisqr, par, 0);
    Double_t Chisqr = CalcChi2(par, ndf, false);
    if(!TMath::Finite(Chisqr) || Chisqr > max_allowed_chi2){
      PrintHelixFitAbnormalWarning("HelixFitInvertCharge", Chisqr, gNumOfHits, itry);
      delete minuit;
      return false;
    }
    if(gChisqr>=Chisqr || TMath::Abs(gChisqr-Chisqr) < 0.01){
      gChisqr = Chisqr;
      gPar[kHelixCx] = par[kHelixCx];
      gPar[kHelixCy] = par[kHelixCy];
      gPar[kHelixZ0] = par[kHelixZ0];
      gPar[kHelixR]  = par[kHelixR];
      gPar[kHelixDz] = par[kHelixDz];
      gMinuitStatus = icstat;
      gBadHits = NBadFromNdf(gNumOfHits, ndf, gMomConstraint ? 1 : 0);
      status = true;
    }
    arglist[0] = arglist[0]*5.;
    arglist[1] = arglist[1]*0.1;
    ++itry;
  }
  delete minuit;

#if DebugDisp
  std::cout<<"HelixFitInvertCharge() status="<<status<<" gChisqr "<<gChisqr<<std::endl;
  if(gMinuitStatus==0) std::cout<<"HelixFit() icstat==0"<<std::endl;
#endif
  return status;
}

//______________________________________________________________________________
static inline Double_t CalcChi2wVertex(Double_t *HelixPar, Int_t &ndf)
{

  if(!ValidateEqualSizes("TPCLocalTrackHelix::CalcChi2wVertex",
                        static_cast<std::size_t>(gNumOfHits),
                        "gHitPos", gHitPos.size(),
                        "gHelixTheta", gHelixTheta.size(),
                        "gLayer", gLayer.size(),
                        "gPadTheta", gPadTheta.size(),
                        "gResParam", gResParam.size()))
    return TMath::QuietNaN();

  ndf = 0; Double_t chisqr = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 d = ResidualVect(HelixPar, gHitPos[i], gHelixTheta[i]);
    TVector3 res = gRes[i];
    if (tpc::IsDummyResolutionVec(res)) continue; // exclude bad clusters
    chisqr += TMath::Sq( TMath::Hypot(d.x(), d.z()) / TMath::Hypot(res.x(), res.z()) );
    ndf++;
    chisqr += TMath::Sq( d.y() / res.y() );
    ndf++;
  }
  TVector3 d = ResidualVect(HelixPar, gVertex, -2.*TMath::Pi(), 2.*TMath::Pi());
  TVector3 res = gVertexRes;
  chisqr += TMath::Sq( TMath::Hypot(d.x(), d.z()) / TMath::Hypot(res.x(), res.z()) );
  ndf++;
  chisqr += TMath::Sq( d.y() / res.y() );
  ndf++;

  if(ndf < 6) return 1.e+10;
  return chisqr/(Double_t)(ndf-5);
}

//______________________________________________________________________________
static inline Bool_t HelixFitwVertex(){

  if(!ValidateEqualSizes("TPCLocalTrackHelix::HelixFitwVertex",
                        static_cast<std::size_t>(gNumOfHits),
                        "gHitPos", gHitPos.size(),
                        "gHelixTheta", gHelixTheta.size(),
                        "gLayer", gLayer.size(),
                        "gPadTheta", gPadTheta.size(),
                        "gResParam", gResParam.size()))
    return false;

  Double_t par[5] = {gPar[kHelixCx], gPar[kHelixCy], gPar[kHelixZ0], gPar[kHelixR], gPar[kHelixDz]};
  Double_t err[5] = {-999., -999., -999., -999., -999.};

  gChisqr = 1.e+10;

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_helixwVertex);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist, 1, ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  minuit->Command("SET STRategy 0");
  //arglist[0] = 5000.;
  //arglist[1] = 0.01;
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  Int_t Err;
  Double_t bnd1, bnd2;

  Bool_t status = false;
  Int_t ndf=0; Int_t itry=0; gMinuitStatus = 0;
  while(true){
    if(itry>MaxTryMinuit) break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

    //minuit->mnprin(4, amin);
    for(Int_t i=0; i<5; i++){
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    }

    //Double_t grad[5];
    //minuit -> Eval(5, grad, Chisqr, par, 0);
    Double_t Chisqr = CalcChi2wVertex(par, ndf);
    if(gChisqr>=Chisqr || TMath::Abs(gChisqr-Chisqr) < 0.01){
      gChisqr = Chisqr;
      gPar[kHelixCx] = par[kHelixCx];
      gPar[kHelixCy] = par[kHelixCy];
      gPar[kHelixZ0] = par[kHelixZ0];
      gPar[kHelixR] = par[kHelixR];
      gPar[kHelixDz] = par[kHelixDz];
      gMinuitStatus = icstat;
      gBadHits = NBadFromNdf(gNumOfHits, ndf, 2);
      status = true;
    }
    arglist[0] = arglist[0]*5.;
    arglist[1] = arglist[1]*0.1;
    ++itry;
  }
  delete minuit;

#if DebugDisp
  std::cout<<"HelixFitwVertex() status="<<status<<" gChisqr "<<gChisqr<<std::endl;
  if(gMinuitStatus==0) std::cout<<"HelixFit() icstat==0"<<std::endl;
#endif

  return status;
}

//Only for momentum constraint fitting
//______________________________________________________________________________
static inline void fcn_circle(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // par[0]=cx, par[1]=cy, par[2]=r (3-param circle fit; maps to helix [kHelixCx, kHelixCy, kHelixR])
  f = 0.;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 pos = GlobalToLocal(gHitPos[i]);

    TVector3 localpos(pos.x(), pos.y(), 0.);
    TVector3 center(par[0], par[1], 0.);
    TVector3 from_center = localpos - center; // radial vector: center -> hit
    f += TMath::Sq(from_center.Mag() - par[2]);
  }
}

//______________________________________________________________________________
static inline Double_t CircleFit(Int_t IsBeam)
{

  if(!ValidateEqualSizes("TPCLocalTrackHelix::CircleFit",
                        static_cast<std::size_t>(gNumOfHits),
                        "gHitPos", gHitPos.size()))
    return -1;

  // par_circ[3]: Minuit circle fit (cx, cy, r) = helix par[kHelixCx, kHelixCy, kHelixR]
  Double_t par_circ[3] = {gPar[kHelixCx], gPar[kHelixCy], gPar[kHelixR]};
  Double_t err_circ[3] = {-999., -999., -999.};

#if DebugDisp
  DebugHelixPar(Form("TPCLocalTrackHelix::CircleFit before MIGRAD (IsBeam=%d)", IsBeam),
                par_circ[0], par_circ[1], par_circ[2]);
#endif

  //cx, cy, r
  Double_t fit_step[3] = {1.0e-4, 1.0e-4, 1.0e-4};
  Double_t low_limit_beam[3] = {gPar[kHelixCx] - 500., gPar[kHelixCy] - 500., 2000.};// about 0.6 GeV/c
  Double_t up_limit_beam[3]  = {gPar[kHelixCx] + 500., gPar[kHelixCy] + 500., 4167.};// about 1.25 GeV/c
  Double_t low_limit[3] = {gPar[kHelixCx] - 500., gPar[kHelixCy] - 500., gPar[kHelixR] - 3000.};
  Double_t up_limit[3]  = {gPar[kHelixCx] + 500., gPar[kHelixCy] + 500., gPar[kHelixR] + 3000.};

  TMinuit *minuit = new TMinuit(3);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_circle);

  Int_t ierflg_circ = 0;
  Double_t arglist_circ[10];
  arglist_circ[0] = 3.52; //for 3 parameters
  minuit->mnexcm("SET ERR", arglist_circ, 1, ierflg_circ);
  arglist_circ[0] = 1;
  minuit->mnexcm("SET NOW", arglist_circ, 1, ierflg_circ); // No warnings

  TString name[3] = {"cx", "cy", "r"};
  for(Int_t i=0; i<3; i++){
    if(IsBeam==1) minuit->mnparm(i, name[i], par_circ[i], fit_step[i], low_limit_beam[i], up_limit_beam[i], ierflg_circ);
    else minuit->mnparm(i, name[i], par_circ[i], fit_step[i], low_limit[i], up_limit[i], ierflg_circ);
  }
  if(gMomConstraint) minuit->FixParameter(2);

  minuit->Command("SET STRategy 0");

  arglist_circ[0] = 1000*5*5*5;
  arglist_circ[1] = 0.1/(10.*10.*10.);
  Int_t err_circ_status;
  Double_t bnd1_circ, bnd2_circ;

  minuit->mnexcm("MIGRAD", arglist_circ, 2, ierflg_circ);
  for(Int_t i=0; i<3; i++){
    minuit->mnpout(i, name[i], par_circ[i], err_circ[i], bnd1_circ, bnd2_circ, err_circ_status);
  }
  Double_t grad[3];
  Double_t chisqr = 0.;
  minuit->Eval(3, grad, chisqr, par_circ, 0);
  delete minuit;

  gPar[kHelixCx] = par_circ[0];
  gPar[kHelixCy] = par_circ[1];
  gPar[kHelixR] = par_circ[2];

#if DebugDisp
  DebugHelixPar(Form("TPCLocalTrackHelix::CircleFit after MIGRAD (IsBeam=%d)\n    chisqr: %g\n    err status (last par): %d",
                     IsBeam, chisqr, err_circ_status),
                par_circ[0], par_circ[1], par_circ[2]);
#endif

  return chisqr;
}

//______________________________________________________________________________
TPCLocalTrackHelix::TPCLocalTrackHelix()
  : m_is_fitted(false),
    m_is_calculated(false),
    m_is_theta_calculated(false),
    m_is_fitted_exclusive(false),
    m_is_multiloop(false),
    m_hit_order(), m_hit_t(),
    m_cx(0.), m_cy(0.), m_z0(0.), m_r(0.), m_dz(0.),
    m_pid(0),
    m_closedist(1.e+10, 1.e+10, 1.e+10),
    m_closedistXZ(1.e+10, 1.e+10, 1.e+10),
    m_chisqr(1.e+10),
    m_minuit(0),
    m_n_iteration(0),
    m_mom0(0.,0.,0.),
    m_edgepoint(0.,0.,-143.),
    m_min_t(0.), m_max_t(0.),
    m_path(0.), m_transverse_path(0.),
    m_charge(0), m_fitflag(0), m_vtxflag(0),
    m_isBeam(0), m_isK18(0), m_isAccidental(0), m_allow_target_crossing_merge(false),
    m_trackid(-1),
    m_ncl_beforetgt(-1),
    m_searchtime(0), m_fittime(0),
    m_MomResScale(-1),
    m_dZResScale(-1),
    m_PhResScale(-1),
    m_cx_exclusive(), m_cy_exclusive(), m_z0_exclusive(),
    m_r_exclusive(), m_dz_exclusive(),
    m_chisqr_exclusive(),
    m_t_exclusive(),
    m_vp()
{

  static const Double_t MomResScale = gUser.GetParameter("MomResScale");
  static const Double_t dZResScale = gUser.GetParameter("dZResScale");
  static const Double_t PhiResScale = gUser.GetParameter("PhiResScale");
  m_MomResScale = MomResScale;
  m_dZResScale = dZResScale;
  m_PhResScale = PhiResScale;

  m_hit_array.reserve(ReservedNumOfHits);
  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
TPCLocalTrackHelix::~TPCLocalTrackHelix()
{
  debug::ObjectCounter::decrease(ClassName());
}

//______________________________________________________________________________
TPCLocalTrackHelix::TPCLocalTrackHelix(TPCLocalTrackHelix *init){

  this -> m_is_fitted = false;
  this -> m_is_calculated = false;
  this -> m_is_fitted_exclusive = false;

  for(Int_t i=0;i<init -> m_hit_array.size();i++){
    this -> m_hit_array.push_back(new TPCLTrackHit(init -> m_hit_array[i] -> GetHit()));
  }

  this -> m_cx = init -> m_cx ;
  this -> m_cy = init -> m_cy ;
  this -> m_z0 = init -> m_z0 ;
  this -> m_r = init -> m_r ;
  this -> m_dz = init -> m_dz ;
  this -> m_pid = init -> m_pid ;
  this -> m_closedist = init -> m_closedist ;
  this -> m_closedistXZ = init -> m_closedistXZ ;
  this -> m_chisqr = init -> m_chisqr ;
  this -> m_minuit = init -> m_minuit ;
  this -> m_n_iteration = init -> m_n_iteration ;
  this -> m_mom0 = init -> m_mom0 ;
  this -> m_edgepoint = init -> m_edgepoint ;
  this -> m_min_t = init -> m_min_t ;
  this -> m_max_t = init -> m_max_t ;
  this -> m_path = init -> m_path ;
  this -> m_transverse_path = init -> m_transverse_path ;
  this -> m_charge = init -> m_charge ;
  this -> m_fitflag = init -> m_fitflag ;
  this -> m_vtxflag = init -> m_vtxflag ;
  this -> m_isBeam = init -> m_isBeam ;
  this -> m_isK18 = init -> m_isK18 ;
  this -> m_isAccidental = init -> m_isAccidental ;
  this -> m_allow_target_crossing_merge = init -> m_allow_target_crossing_merge ;
  this -> m_trackid = init -> m_trackid ;
  this -> m_ncl_beforetgt = init -> m_ncl_beforetgt ;
  this -> m_searchtime = init -> m_searchtime ; //millisec
  this -> m_fittime = init -> m_fittime ; //millisec
  this -> m_is_multiloop = init -> m_is_multiloop ;

  for(Int_t i=0;i<init -> m_hit_order.size();i++){
    this -> m_hit_order.push_back(init -> m_hit_order[i]);
  }

  this -> m_is_theta_calculated = init -> m_is_theta_calculated;

  for(Int_t i=0;i<init -> m_hit_t.size();i++){
    this -> m_hit_t.push_back(init -> m_hit_t[i]);
  }

  for(Int_t i=0;i<init -> m_cx_exclusive.size();i++){
    this -> m_cx_exclusive.push_back(init -> m_cx_exclusive[i]);
    this -> m_cy_exclusive.push_back(init -> m_cy_exclusive[i]);
    this -> m_z0_exclusive.push_back(init -> m_z0_exclusive[i]);
    this -> m_r_exclusive.push_back(init -> m_r_exclusive[i]);
    this -> m_dz_exclusive.push_back(init -> m_dz_exclusive[i]);
    this -> m_chisqr_exclusive.push_back(init -> m_chisqr_exclusive[i]);
    this -> m_t_exclusive.push_back(init -> m_t_exclusive[i]);

    this -> m_hit_array[i] -> SetCalHelixExclusive(init -> m_cx_exclusive[i],
						   init -> m_cy_exclusive[i],
						   init -> m_z0_exclusive[i],
						   init -> m_r_exclusive[i],
						   init -> m_dz_exclusive[i]);
    this -> m_hit_array[i] -> SetThetaExclusive(init -> m_t_exclusive[i]);
    this -> m_hit_array[i] -> SetCalPositionExclusive(init -> m_hit_array[i] -> GetLocalCalPosHelixExclusive());
  }

  for(Int_t i=0;i<init -> m_vp.size();i++){
    this -> m_vp.push_back(init -> m_vp[i]);
  }

  for(Int_t i=0;i<this -> m_hit_array.size();i++){
    this -> m_hit_array[i] -> SetCalHelix(init -> m_cx, init -> m_cy, init -> m_z0, init -> m_r, init -> m_dz);
    this -> m_hit_array[i] -> SetTheta(init -> m_hit_array[i] -> GetTheta());
    this -> m_hit_array[i] -> SetCalPosition(init -> m_hit_array[i] -> GetLocalCalPosHelix());
    this -> m_hit_array[i] -> SetResolution(init -> m_hit_array[i] -> GetResolutionVect());
  }

  static const Double_t MomResScale = gUser.GetParameter("MomResScale");
  static const Double_t dZResScale = gUser.GetParameter("dZResScale");
  static const Double_t PhiResScale = gUser.GetParameter("PhiResScale");
  m_MomResScale = MomResScale;
  m_dZResScale = dZResScale;
  m_PhResScale = PhiResScale;

  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::ClearHits()
{
  m_hit_array.clear();
  m_hit_order.clear();
  m_hit_t.clear();
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetNHitsEffective() const
{

  Int_t n = 0;
  for(auto h:m_hit_array){
    if(!h->IsGoodForTracking()) continue;
    ++n;
  }
  return n;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::AddTPCHit(TPCLTrackHit *hit)
{
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};

  if(hit->IsGood()){
    m_hit_order.push_back(m_hit_array.size());
    m_hit_array.push_back(hit);
    if(m_is_theta_calculated){
      //Calculate Atan2(y,x) and turns of helix. And convert them into the theta of helix
      TVector3 pos = hit->GetLocalHitPos();
      TVector3 localpos = GlobalToLocal(pos);
      Double_t theta = TMath::ATan2(localpos.y() - m_cy, localpos.x() - m_cx);
      Double_t pitch = 2.*TMath::Pi()*m_r*m_dz;
      Double_t ref_ypos = GetPosition(par, theta).y();
      Double_t turns = TMath::Nint((pos.y() - ref_ypos)/pitch);
      theta += 2.*TMath::Pi()*turns;
      m_hit_t.push_back(theta);
    }
    else m_hit_t.push_back(hit->GetTheta());
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::Calculate()
{
  if(IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME << " "
                << "already called" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalHelix(m_cx, m_cy, m_z0, m_r, m_dz);
    hitp->SetTheta(m_hit_t[i]);
    hitp->SetCalPosition(hitp->GetLocalCalPosHelix());
    if(m_is_fitted){
      hitp->SetResolution(GetResolutionVect(i, true));
    }
  }

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  m_mom0 = CalcHelixMom(par, 0.);
  m_pid = Kinematics::HypTPCdEdxPID(
    GetdEdx(TruncatedMeanRatio),
    static_cast<Double_t>(m_charge)*m_mom0.Mag()
  );
  m_is_calculated = true;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::CalculateExclusive()
{

  if(!IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME
		<< "No inclusive calculation" << std::endl;
    return;
  }

  if(!m_is_fitted_exclusive){
    hddaq::cerr << "#W " << FUNC_NAME
		<< "No exclusive fitting" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalHelixExclusive(
      m_cx_exclusive[i], m_cy_exclusive[i],
      m_z0_exclusive[i], m_r_exclusive[i],
      m_dz_exclusive[i]
    );
    hitp->SetThetaExclusive(m_t_exclusive[i]);
    hitp->SetCalPositionExclusive(hitp->GetLocalCalPosHelixExclusive());
  }
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetNPad() const
{
  // #hits < MinHits
  const std::size_t n = m_hit_array.size();
  std::set<Int_t> pads;
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    pads.insert(hit->GetHit()->GetPad());
  }
  return static_cast<Int_t>(pads.size());
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetNDF() const
{

  Int_t nhit = GetNHit();
  Int_t ndf = 2*nhit - 5;
  if(gMomConstraint) ndf += 1;
  return ndf;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetPosition(const Double_t par[5], Double_t t) const
{
  return GlobalPosition(par, t);
}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::CalcClosestDistTgt()
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return;

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 tgt(0., 0., tpc::Z_TARGET);
  Double_t theta_margin = 0.5*TMath::Pi();
  Double_t theta_tgt = EvalTheta(par, tgt, m_min_t - theta_margin, m_max_t + theta_margin);
  TVector3 closest_point = GlobalPosition(par, theta_tgt);
  TVector3 dist = closest_point - tgt;
  m_closedist = dist;

  TVector3 residual_xz = ResidualVectXZ(par, tgt);
  TVector3 horizontal_closest_point = residual_xz + tgt;
  Double_t theta_tgtXZ = EvalThetaXZ(par, horizontal_closest_point, m_min_t - theta_margin, m_max_t + theta_margin);
  TVector3 closest_pointXZ = GlobalPosition(par, theta_tgtXZ);
  TVector3 distXZ = closest_pointXZ - tgt;
  m_closedistXZ = distXZ;
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetResolutionVect(TPCLTrackHit* hit, Bool_t vetoBadClusters){

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME))
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());

  if(vetoBadClusters && hit->GetHit()->GetParentCluster()->IsOnTheFrame()) return TVector3(3.e+10, 3.e+10, 3.e+10);

  TVector3 pos       = hit->GetLocalHitPos();
  Int_t layer        = hit->GetLayer();
  Double_t pad_theta = hit->GetPadTheta();
  Double_t theta     = hit->GetTheta();
  std::vector<Double_t> res_param = hit->GetResolutionParams();
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};

  //Convert resolution along the row direction into closets point's x, y, z resolutions
  TVector3 res = CalcResolution(par, layer, pos, pad_theta, theta, res_param, vetoBadClusters);
  return res;
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetResolutionVect(Int_t i, Bool_t vetoBadClusters){

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME))
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());

  return GetResolutionVect(m_hit_array[i], vetoBadClusters);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetResolutionY(TPCHit* hit){
  //vertical resolution
  std::vector<Double_t> resparam = hit->GetResolutionParams();
  Double_t param_y[4] = {resparam[6], resparam[1], resparam[7], resparam[8]};
  f_drift -> SetParameters(param_y);
  TVector3 pos = hit->GetPosition();
  Double_t res_drift = f_drift -> Eval(pos.y());
  return res_drift;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetResolutionY(Int_t i){
  TPCHit *hit = m_hit_array[i]->GetHit();
  return GetResolutionY(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetAlpha(Int_t i) const //for a point on the track
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return -1;

  TPCHit *hit = m_hit_array[i]->GetHit();
  return GetAlpha(hit);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetAlpha(TPCHit* hit) const
{

  //find expected position on the track
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 pos = hit->GetPosition();
  TVector3 localpos = GlobalToLocal(pos);
  Double_t tanTrack = (localpos.y()-par[kHelixCy])/(localpos.x()-par[kHelixCx]);
  Double_t padTheta = hit->GetPadTheta();
  Double_t tanPad = TMath::Tan(padTheta);

  //alpha : pad - track angle
  Double_t tanDiff = (tanPad-tanTrack)/(1.+tanPad*tanTrack);
  Double_t alpha = TMath::ATan(tanDiff);
  return alpha;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::CalcResidual(TVector3 pos)
{

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  return ResidualVect(par, pos, -2.*TMath::Pi(), 2.*TMath::Pi());
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::CalcHelixMomCenter(const Double_t par[5]) const
{

  Double_t p_t = HelixRToPt(par[kHelixR]); // [GeV/c]
  // z0 + r*dz*theta = 0  (track z = lab Y)
  Double_t theta = -par[kHelixZ0]/(par[kHelixR]*par[kHelixDz]); // y = 0

  // helix-fit plane [GeV/c]: x = r*cos(theta), y = r*sin(theta)
  Double_t px_local = -p_t * TMath::Sin(theta);
  Double_t py_local =  p_t * TMath::Cos(theta);
  Double_t pz_local =  p_t * par[kHelixDz];

  // lab/track frame [GeV/c]: x=-X, y=Z-z_tgt, z=Y
  Double_t px_global = -px_local;
  Double_t py_global =  pz_local;
  Double_t pz_global =  py_local;

  TVector3 p = TVector3(px_global, py_global, pz_global);
  if(m_charge > 0) p *= -1.;
  return p;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::CalcHelixMom(const Double_t par[5], Double_t theta) const
{

  Double_t p_t = HelixRToPt(par[kHelixR]); // [GeV/c]

  // helix-fit plane [GeV/c]: x = r*cos(theta), y = r*sin(theta)
  Double_t px_local = -p_t * TMath::Sin(theta);
  Double_t py_local =  p_t * TMath::Cos(theta);
  Double_t pz_local =  p_t * par[kHelixDz];

  // lab/track frame [GeV/c]: x=-X, y=Z-z_tgt, z=Y
  Double_t px_global = -px_local;
  Double_t py_global =  pz_local;
  Double_t pz_global =  py_local;

  TVector3 p = TVector3(px_global, py_global, pz_global);
  if(m_charge > 0) p *= -1.;
  return p;
}

//______________________________________________________________________________
TPCLTrackHit*
TPCLocalTrackHelix::GetHit(std::size_t nth) const
{
  if(nth<m_hit_array.size() && nth>=0)
    return m_hit_array[nth];
  else
    return 0;
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::GetOrder(Int_t i) const
{
  Int_t size = m_hit_order.size();
  if(i>=size || i<0) return -1;
  Int_t id = i;
  // E72 configuration.
  if(m_charge > 0) id = size - i - 1;
  Int_t order = m_hit_order[id];
  return order;
}

//______________________________________________________________________________
TPCLTrackHit*
TPCLocalTrackHelix::GetHitInOrder(std::size_t nth) const
{
  Int_t size = m_hit_order.size();
  Int_t id = nth;
  // E72 configuration.
  if(m_charge > 0) id = size - nth - 1;
  Int_t order = m_hit_order[id];
  if(nth>=m_hit_array.size() || nth<0) return 0;
  return m_hit_array[order];
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::DeleteNullHit()
{
  const std::size_t before = m_hit_array.size();
  if(before == 0) return;
  if(!ValidateEqualSizes(FUNC_NAME, before,
                        "m_hit_t", m_hit_t.size(),
                        "m_hit_order", m_hit_order.size()))
    return;

  std::size_t n_removed = 0;
  for(Int_t i = static_cast<Int_t>(before) - 1; i >= 0; --i){
    if(m_hit_array[i]->IsGood()) continue;
    m_hit_array.erase(m_hit_array.begin() + i);
    m_hit_t.erase(m_hit_t.begin() + i);
    RemapHitOrderAfterErase(m_hit_order, i);
    ++n_removed;
  }
  if(n_removed == 0) return;

  m_is_theta_calculated = false;
  hddaq::cout << FUNC_NAME << " "
              << n_removed << " null hit(s) deleted"
              << std::endl;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SetClustersHoughFlag(Int_t hough_flag)
{
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TPCHit *hit = hitp->GetHit();
    if( !hit ) continue;
    hit->SetHoughFlag(hough_flag);
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SetFlag(Int_t flag)
{

  if(flag==0){ m_isBeam=0; m_isK18=0; m_isAccidental=0; } //reset
  if((flag&1)==1) m_isBeam=1;
  if((flag&2)==2) m_isK18=1;
  if((flag&8)==8) m_isAccidental=1;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Int_t MinHits)
{

  gMomConstraint = false; //No momentum constraint
  if(!IsGoodForTracking()) return false;

  Bool_t status = DoHelixTrackFit(); //track chisqr minimization
  m_is_fitted = status;

  if(!status || m_chisqr > MaxChisqr) return false;

  SeparateTracksAtTarget();

  //If track speration is performed, m_is_fitted flag is set to false
  if(!m_is_fitted) return DoFit(MinHits); //Do chisqr minimization again after separation

  //Minimum # of clusters
  Int_t nhit = GetNHit();
#if 1
  if(IsBackward() && GetIsBeam()==1) MinHits = MinHitsBackwardBeam;
  if(nhit<MinHits) return false;
#else
  Int_t nbadhit = gBadHits;
  if(nhit-nbadhit<MinHits) return false;
#endif

  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoHelixTrackFit(Double_t RKpar[5]) //for beam track fitting
{

  if(!gMomConstraint){
    std::cout<<FUNC_NAME+" Fatal error : No momentum constraint"<<std::endl;
    return false;
  }

  Bool_t vetoBadClusters = false;

  DeleteNullHit();
  const std::size_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gHelixTheta.clear();
  gBadHits = 0;
  gMinuitStatus = 0;
  gChisqr = 1.e+10;
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
    Int_t layer = hitp->GetLayer();
    gLayer.push_back(layer);
    Double_t padTheta = hitp->GetPadTheta();
    gPadTheta.push_back(padTheta);
    std::vector<Double_t> resparam = hitp->GetResolutionParams();
    gResParam.push_back(resparam);
  }

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " Inital helix params before fitting",
                m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  //1. calculate chisqr by using RK parameters
  SetParam(RKpar);
  CalcHelixTheta();
  Int_t ndf;
  Double_t RKChisqr = CalcChi2(RKpar, ndf, vetoBadClusters);
  if(RKChisqr>10.){
    //2. calculate chisqr with pre-fitting parameters (gChisqr)
    if(!DoPreFit(RKpar)) return false;

#if DebugDisp
    DebugHelixPar(FUNC_NAME + Form(" Helix params after pre-fitting (w/ constraints from RK)\n    chisqr: %g", gChisqr),
                  gPar);
#endif
  }

  //Compare 1 & 2 and choose the better one
  if(gChisqr > RKChisqr){
    gChisqr = RKChisqr;
    gPar[kHelixCx] = RKpar[kHelixCx];
    gPar[kHelixCy] = RKpar[kHelixCy];
    gPar[kHelixZ0] = RKpar[kHelixZ0];
    gPar[kHelixR] = RKpar[kHelixR];
    gPar[kHelixDz] = RKpar[kHelixDz];
  }
  SetParam(gPar);
  CalcHelixTheta();

  //Helix fitting
  if(!DoHelixFit(gPar, vetoBadClusters)) return false;

#if DebugDisp
  DebugHelixPar(FUNC_NAME + Form(" After helix fitting\n    n_iteration: %d\n    chisqr: %g",
                                 m_n_iteration, gChisqr), gPar);
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
  if(false_layer < 0) return false;
  Bool_t goodtrack = true;
  if(false_layer > 0 || m_chisqr > MaxChisqr){
    goodtrack = false;

#if DebugDisp
    TPCLTrackHit *hitp = m_hit_array[delete_hit];
    TVector3 pos = hitp->GetLocalHitPos();
    std::cout<<"delete hits ["<<delete_hit<<"]=("
    	     <<pos.x()<<", "
      	     <<pos.y()<<", "
      	     <<pos.z()<<")"<<std::endl;
#endif

    EraseHit(delete_hit);
    gHelixTheta.erase(gHelixTheta.begin()+delete_hit);
  }

  m_n_iteration++;
  if(m_n_iteration > MaxIteration) return false;
  if(goodtrack){ //Tracking is over.
#if IterativeResolution
    vetoBadClusters = true;

    //Now excluding bad clusters and fitting again.
    DoHelixFit(gPar, vetoBadClusters);
#if DebugDisp
    DebugHelixPar(FUNC_NAME + Form(" After excluding bad hits\n    n_iteration: %d\n    chisqr: %g",
                                   m_n_iteration, gChisqr), gPar);
#endif
#endif
    return true;
  }
  else return DoHelixTrackFit(RKpar);
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Double_t RKpar[5], Int_t MinHits)
{

  gMomConstraint = true; //fit with the momentum constraint
  if(!IsGoodForTracking()) return false;

  Bool_t status = DoHelixTrackFit(RKpar); //track chisqr minimization
  m_is_fitted = status;
  if(!status || m_chisqr > MaxChisqr) return false;

  //SeparateTracksAtTarget();

  //Minimum # of clusters
  Int_t nhit = GetNHit();
#if 1
  if(IsBackward() && GetIsBeam()==1) MinHits = MinHitsBackwardBeam;
  if(nhit<MinHits) return false;
#else
  Int_t nbadhit = gBadHits;
  if(nhit-nbadhit<MinHits) return false;
#endif
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Double_t RKCharge, Double_t RKpar[5], Int_t MinHits)
{

  gMomConstraint = true; //fit with a momentum constraint
  Bool_t status = DoFit(RKpar, MinHits);
  Int_t charge = (int) RKCharge;
  if(m_charge != charge) return false; //fit with a charge constraint
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFit(Double_t RKCharge, Int_t MinHits)
{

  Bool_t status = DoFit(MinHits);
  Int_t charge = (int) RKCharge;
  if(m_charge != charge) return false; //fit with a charge constraint
  return status;
}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::EraseHits(std::vector<Int_t> delete_hits)
{
  std::sort(delete_hits.begin(), delete_hits.end());
  for(Int_t i=0; i<delete_hits.size(); ++i){
    //Reset houghflag
    EraseHit(delete_hits[i]-i);
  }
  m_is_theta_calculated = false;
}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::EraseHit(Int_t delete_hit)
{

  //Reset houghflag
  TPCLTrackHit *hitp = m_hit_array[delete_hit];
  TPCHit *hit = hitp->GetHit();
  hit->SetHoughFlag(0);
  m_hit_array.erase(m_hit_array.begin()+delete_hit);
  m_hit_t.erase(m_hit_t.begin()+delete_hit);
  RemapHitOrderAfterErase(m_hit_order, delete_hit);

  m_is_theta_calculated = false;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsGoodForTracking()
{

  const std::size_t n = m_hit_array.size();
  if(GetNDF()<1){
#if DebugDisp
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "Min layer should be > NDF" << std::endl;
#endif
    return false;
  }

  if(n>ReservedNumOfHits){
#if DebugDisp
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "n > ReservedNumOfHits" << std::endl;
#endif
    return false;
  }

  Int_t npad = GetNPad();
  if(npad<=1) return false; //frequently many noise hits appear in a single pad

  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoCircleFit(Double_t *par)
{

  DeleteNullHit();
  if(!IsGoodForTracking()) return false;
  const std::size_t n = m_hit_array.size();
  gNumOfHits = n;

  Double_t xp[n]; Double_t yp[n];
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 localpos = GlobalToLocal(pos);
    xp[i] = localpos.X();
    yp[i] = localpos.Y();
  }

  Double_t par_circ[3]={0.};
  Int_t npad = GetNPad();
  if(npad==3) return true;
  if(npad<3 || CircleFit(xp, yp, n, &par_circ[0], &par_circ[1], &par_circ[2]) < 0.) return false;

  par[kHelixCx] = par_circ[0];
  par[kHelixCy] = par_circ[1];
  par[kHelixR]  = par_circ[2];

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " Circlefit results", par[kHelixCx], par[kHelixCy], par[kHelixR]);
#endif

  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoCircleFitwMomConstraint(Double_t *par)
{

  if(!gMomConstraint){
    std::cout<<FUNC_NAME+" Fatal error : No momentum constraint"<<std::endl;
    return false;
  }

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " Helix params before circle fitting w/ the momentum constraint",
                par[kHelixCx], par[kHelixCy], par[kHelixZ0], par[kHelixR], par[kHelixDz]);
#endif

  DeleteNullHit();
  if(!IsGoodForTracking()) return false;
  const std::size_t n = m_hit_array.size();
  gNumOfHits = n;
  gHitPos.clear();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
  }

  gPar[kHelixCx] = par[kHelixCx];
  gPar[kHelixCy] = par[kHelixCy];
  gPar[kHelixR] = par[kHelixR];
  CircleFit(m_isBeam);

  par[kHelixCx] = gPar[kHelixCx];
  par[kHelixCy] = gPar[kHelixCy];
  par[kHelixR] = gPar[kHelixR];

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " Circlefit results", par[kHelixCx], par[kHelixCy], par[kHelixR]);
#endif

  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoStraightLineFit(Double_t *par)
{
  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

  DeleteNullHit();
  const std::size_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
    Int_t layer = hitp->GetLayer();
    gLayer.push_back(layer);
    Double_t padTheta = hitp->GetPadTheta();
    gPadTheta.push_back(padTheta);
    std::vector<Double_t> resparam = hitp->GetResolutionParams();
    gResParam.push_back(resparam);
  }

  Bool_t vetoBadClusters = false;
  gRes.clear();
  for(Int_t i=0; i<gNumOfHits; i++){
    Double_t dummy_theta = 0.;
    TVector3 res = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], dummy_theta, gResParam[i], vetoBadClusters);
    gRes.push_back(res);
  }

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " Helix params before straight-line fitting", gPar);
#endif

  gPar[kHelixCx] = par[kHelixCx];
  gPar[kHelixCy] = par[kHelixCy];
  gPar[kHelixZ0] = par[kHelixZ0];
  gPar[kHelixR] = par[kHelixR];
  gPar[kHelixDz] = par[kHelixDz];

#if 1 //Spiral-like track with low p_T
  if(gMultiLoop) return true;
#endif
  return StraightLineFit();
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoPreFit(Double_t par[5])
{

  Bool_t vetoBadClusters = false;

  gPar[kHelixCx] = par[kHelixCx];
  gPar[kHelixCy] = par[kHelixCy];
  gPar[kHelixZ0] = par[kHelixZ0];
  gPar[kHelixR] = par[kHelixR];
  gPar[kHelixDz] = par[kHelixDz];

  //pre circle fit
  Bool_t pass = DoCircleFit(gPar);
  if(!pass && gMomConstraint) pass = DoCircleFitwMomConstraint(gPar);
  if(!pass) return false;

#if 1
  const Int_t n = m_hit_array.size();
  gNumOfHits = n;
  gHitPos.clear();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
  }

  Int_t MaxBin[3];
  if(pass && (gPar[kHelixCx] < LowLimit[kHelixCx] ||
              gPar[kHelixCx] > UpLimit[kHelixCx]  ||
              gPar[kHelixCy] < LowLimit[kHelixCy] ||
              gPar[kHelixCy] > UpLimit[kHelixCy]  ||
              gPar[kHelixR]  < LowLimit[kHelixR] ||
              gPar[kHelixR]  > UpLimit[kHelixR]))
    pass = tpc::HoughTransformCircleXZ(gHitPos, MaxBin, gPar, 3);
  if(!pass) return false; //For very high momentum tracks
#endif
  SetParam(gPar);
  CalcHelixTheta();

#if 1 //Optional
  Int_t MaxBinY[3];
  if(pass && (gPar[kHelixZ0] < LowLimit[kHelixZ0] ||
              gPar[kHelixZ0] > UpLimit[kHelixZ0]  ||
              gPar[kHelixDz] < LowLimit[kHelixDz] ||
              gPar[kHelixDz] > UpLimit[kHelixDz]))
    tpc::HoughTransformLineYTheta(gHitPos, MaxBinY, gPar, 1000.);
#endif
  if(!DoStraightLineFit(gPar)) return false;
  SetParam(gPar);

  Int_t ndf;
  gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoHelixFit(Double_t *par, Bool_t vetoBadClusters)
{
  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

  gPar[kHelixCx] = par[kHelixCx];
  gPar[kHelixCy] = par[kHelixCy];
  gPar[kHelixZ0] = par[kHelixZ0];
  gPar[kHelixR] = par[kHelixR];
  gPar[kHelixDz] = par[kHelixDz];

  IsBackward(); // may SetIsBeam() for HelixFit limits; return value unused

  Bool_t status = false;
  if(HelixFit(m_isBeam, vetoBadClusters)){
    status = true;
    SetParam(gPar);
    CalcHelixTheta();
    m_chisqr = gChisqr;
    m_minuit = gMinuitStatus;
    Double_t window = HelixThetaSearchWindow(gPar[kHelixR]);
    for(Int_t i=0; i<gNumOfHits; ++i){
      Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
      m_hit_t[i] = theta;
      gHelixTheta[i] = theta;
    }
  }
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoHelixTrackFit()
{

  if(gMomConstraint){
    std::cout<<FUNC_NAME+" Fatal error : Momentum constraint is applied"<<std::endl;
    return false;
  }

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " Helix params (Hough-Transform)", m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  Bool_t vetoBadClusters = false;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gBadHits = 0;
  gMinuitStatus = 0;

  //Initialization params with Hough-transform result or previous tracking result
  gPar[kHelixCx] = m_cx;
  gPar[kHelixCy] = m_cy;
  gPar[kHelixZ0] = m_z0;
  gPar[kHelixR] = m_r;
  gPar[kHelixDz] = m_dz;

  if(!DoPreFit(gPar)) return false;
#if DebugDisp
  DebugHelixPar(FUNC_NAME + Form(" After pre fitting\n    n_iteration: %d\n    chisqr: %g",
                                 m_n_iteration, gChisqr), gPar);
#endif

  // Helix fitting
  if(!DoHelixFit(gPar, vetoBadClusters)) return false;
#if DebugDisp
  DebugHelixPar(FUNC_NAME + Form(" After helix fitting\n    n_iteration: %d\n    chisqr: %g",
                                 m_n_iteration, gChisqr), gPar);
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
  if(false_layer < 0) return false;

  Bool_t goodtrack = true;
  if(false_layer > 0 || m_chisqr > MaxChisqr){
    goodtrack = false;
#if DebugDisp
    TPCLTrackHit *hitp = m_hit_array[delete_hit];
    TVector3 pos = hitp->GetLocalHitPos();
    std::cout<<"delete hit #"<<delete_hit<<"=("
	     <<pos.x()<<", "
	     <<pos.y()<<", "
	     <<pos.z()<<")"<<std::endl;
#endif
    EraseHit(delete_hit);
    gHelixTheta.erase(gHelixTheta.begin()+delete_hit);
  }

  m_n_iteration++;
  if(m_n_iteration > MaxIteration) return false;
  if(goodtrack){ //Tracking is over.
#if IterativeResolution
    vetoBadClusters = true;

    //Now excluding bad clusters and fitting again.
    DoHelixFit(gPar, vetoBadClusters);
#if DebugDisp
    DebugHelixPar(FUNC_NAME + Form(" with precise resolution calculation\n    n_iteration: %d\n    chisqr: %g",
                                   m_n_iteration, gChisqr), gPar);
#endif

#endif
    return true;
  }
  else return DoHelixTrackFit();
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ResidualCheck(Int_t i, Double_t &residual)
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

  TPCHit *hit = m_hit_array[i] -> GetHit();
  Int_t layer = hit->GetLayer();
  Double_t pad_theta = hit->GetPadTheta();
  std::vector<Double_t> res_param = hit->GetResolutionParams();

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 position = m_hit_array[i]->GetLocalHitPos();
  TVector3 res  = CalcResolution(par, layer, position, pad_theta, gHelixTheta[i], res_param, false);
  TVector3 resi = ResidualVect(par, position, gHelixTheta[i]); // Closest distance
  residual = resi.Mag();

  //XZ residual < window
  TVector3 residual_xz = ResidualVectXZ(par, position);
  if(m_is_multiloop){
    if(residual_xz.Mag() > ResidualWindowOutXZ) return false;
  }
  else if(layer <= tpc::LAST_TGT_LAYER && TMath::Abs(position.y()) < tpc::TARGET_HALF_Y){
    if(residual_xz.Mag() > ResidualWindowUnderTgtXZ) return false;
  }
  else if(layer < 10 && residual_xz.Mag() > ResidualWindowInXZ) return false;
  else if(layer >=10 && residual_xz.Mag() > ResidualWindowOutXZ) return false;

  //Residual/resolution < window
  Double_t residual_vertical = resi.y();
  Double_t residual_horizontal = TMath::Hypot(resi.x(), resi.z());
  Double_t resolution_vertical = res.y();
  Double_t resolution_horizontal = TMath::Hypot(res.x(), res.z());
  if(TMath::Abs(residual_vertical) > resolution_vertical*ResidualWindowPullY) return false;
  if(TMath::Abs(residual_horizontal) > resolution_horizontal*ResidualWindowPullXZ) return false;
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsGoodHitToAdd(TPCHit *hit, Double_t &residual, Bool_t nolimitation)
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

  Int_t layer = hit->GetLayer();
  Double_t padTheta = hit->GetPadTheta();
  std::vector<Double_t> resparam = hit->GetResolutionParams();

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 position = hit->GetPosition();
  Double_t dummy_theta = 0;
  TVector3 res = CalcResolution(par, layer, position, padTheta, dummy_theta, resparam, false);

  // to do: it might be better to optimize for E72 target geometry
  Int_t upstream_tgt = -1;
  if(IsBeamLikeHit(position)) upstream_tgt = 0; // beam section
  else if(position.z() < tpc::Z_TARGET) upstream_tgt = 1; // upstream, outside beam box

  Int_t section = tpc::GetSection(position.x(), position.z());

  Double_t factor = 1.;
  Int_t nhit_upstream_tgt = 0;
  std::vector<Int_t> gSection;
  const std::size_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    const Int_t track_section = hitp->GetSection();
    gSection.push_back(track_section);

    TVector3 pos = hitp -> GetLocalHitPos();
    if(IsBeamLikeHit(pos)) nhit_upstream_tgt++;
  }

  Double_t max_scanrange = 40.;
  if(find(gSection.begin(), gSection.end(), section) == gSection.end()){ //Crossing to the other sections
    if(TMath::Abs(position.y()) > 50) factor = 2.5;
    else factor = 2.;

    // track_lacks_gem[i]: no hit on this track in GEM section (i + 1), i = 0..3.
    Bool_t track_lacks_gem[4];
    for (Int_t i = 0; i < 4; ++i)
      track_lacks_gem[i] = (find(gSection.begin(), gSection.end(), i + 1) == gSection.end());

    // When crossing sectors, widen theta scan less if the "across" pair is already missing.
    if (section == 1 || section == 3) {
      if (track_lacks_gem[1] || track_lacks_gem[3]) max_scanrange += 60.; // sections 2 and 4
      else max_scanrange += 120.;
    }
    else if (section == 2 || section == 4) {
      if (track_lacks_gem[0] || track_lacks_gem[2]) max_scanrange += 60.; // sections 1 and 3
      else max_scanrange += 120.;
    }
  }
  else if(nhit_upstream_tgt==0 && upstream_tgt==1 &&
         (TMath::Hypot(m_closedist.x(), m_closedist.z())<25. &&
          TMath::Abs(m_closedist.y())>10.)){ //Under or over the target (todo: modify for E72 target geometry)
    factor = 1.;
    max_scanrange += 40.;
  }
  else if(nhit_upstream_tgt==0 && upstream_tgt==0){ //Beam section
    factor = 1.;
    max_scanrange = 20.;
  }

  if(nolimitation){
    max_scanrange = 300.;
    factor = 2.;
  }
  max_scanrange /= m_r;

  Double_t theta_range[2];
  Double_t max_scanrange_theta = m_is_multiloop ? TMath::Pi() : 0.5*TMath::Pi();
  theta_range[0] = m_min_t - TMath::Min(max_scanrange_theta, max_scanrange);
  theta_range[1] = m_max_t + TMath::Min(max_scanrange_theta, max_scanrange);

  //XZ residual < window
  TVector3 residual_xz = ResidualVectXZ(par, position);
#if 1
  if(m_r < 250.){
    if(residual_xz.Mag() > 3.*ResidualWindowOutXZ) return false;
  }
#else
  if(m_is_multiloop){
    if(residual_xz.Mag() > factor*ResidualWindowOutXZ) return false;
  }
#endif
  else if(layer <= tpc::LAST_TGT_LAYER && TMath::Abs(position.y()) < tpc::TARGET_HALF_Y){
    if(residual_xz.Mag() > factor*ResidualWindowUnderTgtXZ) return false;
  }
  else if(layer <  10 && residual_xz.Mag() > factor*ResidualWindowInXZ) return false;
  else if(layer >= 10 && residual_xz.Mag() > factor*ResidualWindowOutXZ) return false;

  TVector3 resi = ResidualVect(par, position, theta_range[0], theta_range[1]);
  residual = resi.Mag();

  //Residual/resolution < window
  Double_t residual_vertical = resi.y();
  Double_t residual_horizontal = TMath::Hypot(resi.x(), resi.z());
  Double_t resolution_vertical = res.y();
  Double_t resolution_horizontal = TMath::Hypot(res.x(), res.z());

#if 1
  if(nolimitation){
    if(TMath::Abs(residual_horizontal) > 5.*resolution_horizontal*ResidualWindowPullXZ) return false;
    if(TMath::Abs(residual_vertical) > 3.*resolution_vertical*ResidualWindowPullY) return false;
  }
  else if(m_r < 250.){ //wider cut condition for low-momentum tracks (currently not supported)
    if(TMath::Abs(residual_horizontal) > factor*resolution_horizontal*ResidualWindowPullXZ) return false;
    if(TMath::Abs(residual_vertical) > factor*resolution_vertical*ResidualWindowPullY) return false;
  }
  else{
    if(TMath::Abs(residual_horizontal) > factor*resolution_horizontal*ResidualWindowPullXZ) return false;
    if(TMath::Abs(residual_vertical) > factor*resolution_vertical*ResidualWindowPullY) return false;
  }
#else
  if(TMath::Abs(residual_horizontal) > factor*resolution_horizontal*ResidualWindowPullXZ) return false;
  if(TMath::Abs(residual_vertical) > factor*resolution_vertical*ResidualWindowPullY) return false;
#endif
  return true;
}

//_____________________________________________________________________________
Int_t
TPCLocalTrackHelix::Side(TVector3 hitpos)
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return 0;

  if(!ValidateTrackHasHits(*this, FUNC_NAME)) return 0;

  //TPC local coordinate
  Int_t flag = -1;
  TVector3 Vect1(-hitpos.X() - m_cx, hitpos.Z() - tpc::Z_TARGET - m_cy, 0.); //Vect1(Hit - Helix center)
  TVector3 Vect2(-m_cx, -m_cy, 0.); //Vec2(Tgt - Helix center)

  TVector3 norm = Vect1.Cross(Vect2); //Vect1 X Vect2
  if(norm.Z()>0.) flag = 1;
  return flag;
}

//______________________________________________________________________________
Double_t
TPCLocalTrackHelix::CalcThetaExclusive(Int_t ith) const
{

  if(!m_is_fitted_exclusive){
    std::cout<<FUNC_NAME+" Fatal error : Do exclusive fitting first!"<<std::endl;
    return TMath::QuietNaN();
  }

  TPCLTrackHit *hitp = m_hit_array[ith];
  TVector3 pos = hitp -> GetLocalHitPos();
  TVector3 localpos = GlobalToLocal(pos);
  Double_t par[5] = {m_cx_exclusive[ith], m_cy_exclusive[ith], m_z0_exclusive[ith], m_r_exclusive[ith], m_dz_exclusive[ith]};
  Double_t window = HelixThetaSearchWindow(par[kHelixR]);
  Double_t theta = EvalTheta(par, pos, m_hit_t[ith] - 0.5*window, m_hit_t[ith] + 0.5*window);

  return theta;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::CalcHelixTheta()
{

  if(!ValidateTrackHasHits(*this, FUNC_NAME)) return;

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
#if DebugDisp
  DebugHelixPar(FUNC_NAME, m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  gHelixTheta.clear();
  std::vector<Double_t> helix_theta;

  // Compare a variance of y positions to determine theta with Atan2 or Atan2 +/- pi.
  Double_t var_base = 0.;
  Double_t var_shifted = 0.;

  Bool_t theta_flip = false;
  Double_t prev_theta = 0.;
  Double_t theta0 = 0.;
  m_min_t = 9999.;
  m_max_t = -9999.;
  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();

    // theta_base: Atan2 output with range (-pi, pi)
    TVector3 localpos = GlobalToLocal(pos);
    Double_t theta_base = TMath::ATan2(localpos.y() - m_cy, localpos.x() - m_cx);
    if(m_is_multiloop){
      // Calculate Atan2(y,x) and turns of helix. And convert them into the theta of helix.
      Double_t pitch = 2.*TMath::Pi()*m_r*m_dz;
      Double_t ref_ypos = GetPosition(par, theta_base).y();
      Double_t turns = TMath::Nint((pos.y() - ref_ypos)/pitch);
      theta_base += 2.*TMath::Pi()*turns;
    }
    else{
      // Check ATan2 function's theta flip (-pi ~ pi) within a loop.
      if(i==0) theta0 = theta_base;
      else if(TMath::Abs(prev_theta - theta_base) > TMath::Pi()) theta_flip = true;

      if(theta_flip){
        if(theta0>0 && theta_base<0) theta_base += 2.*TMath::Pi();
        if(theta0<0 && theta_base>0) theta_base -= 2.*TMath::Pi();
      }
      // resi_*: y residual; theta_shifted: base theta shifted by +/- 2*pi.
      Double_t resi_base = GetPosition(par, theta_base).y() - pos.y();
      Double_t theta_shifted = (theta0 > 0.)
        ? theta_base - 2.*TMath::Pi()
        : theta_base + 2.*TMath::Pi();
      Double_t resi_shifted = GetPosition(par, theta_shifted).y() - pos.y();
      var_base    += TMath::Sq(resi_base);
      var_shifted += TMath::Sq(resi_shifted);
    }
    prev_theta = theta_base;
    helix_theta.push_back(theta_base);
  }

  Bool_t apply_2pi_shift = false;
  if(!m_is_multiloop && n > 0){
    Double_t pitch = TMath::Abs(2.*TMath::Pi()*m_r*m_dz);
    // |2*pi*r*dz| -> 0 for r~0 or dz~0; without a floor,
    // k*pitch ~ 0 and apply_2pi_shift fires too easily.
    pitch = TMath::Max(pitch, ThetaShiftPitchEps);
    apply_2pi_shift = var_base > var_shifted
      && TMath::Sqrt((var_base - var_shifted) / static_cast<Double_t>(n)) > ThetaShiftK * pitch;
  }

  for(std::size_t i=0; i<n; ++i){
    Double_t theta = helix_theta[i];
    if(apply_2pi_shift){
      // Apply the +/- 2*pi shift chosen from var_base vs var_shifted.
      theta += (helix_theta[0] > 0.) ? -2.*TMath::Pi() : 2.*TMath::Pi();
    }
    if(theta < m_min_t) m_min_t = theta;
    if(theta > m_max_t) m_max_t = theta;
    gHelixTheta.push_back(theta);
  }

#if DebugDisp
  std::cout<<"helix theta calculated: theta min "<<m_min_t<<" max "<<m_max_t<<std::endl;
#endif
  m_is_theta_calculated = true;

}

//______________________________________________________________________________
void
TPCLocalTrackHelix::DoFitExclusive()
{
  const std::size_t n = m_hit_array.size();

  if(m_hit_t.size() <= 1 || !m_is_fitted || !m_is_calculated) std::cout<<FUNC_NAME+" Fatal error : Wrong track. Please check"<<std::endl;

  gMultiLoop = m_is_multiloop;

  m_is_fitted_exclusive = true;
  m_cx_exclusive.resize(n);
  m_cy_exclusive.resize(n);
  m_z0_exclusive.resize(n);
  m_r_exclusive.resize(n);
  m_dz_exclusive.resize(n);
  m_t_exclusive.resize(n);
  m_chisqr_exclusive.resize(n);

  gNumOfHits = n-1;
  for(Int_t ihit=0; ihit<n; ++ihit){ //exclude the ith hit
    gHitPos.clear();
    gLayer.clear();
    gPadTheta.clear();
    gResParam.clear();
    gHelixTheta.clear();
    gHitPos.resize(n-1);
    gLayer.resize(n-1);
    gPadTheta.resize(n-1);
    gResParam.resize(n-1);
    gHelixTheta.resize(n-1);

    gChisqr = 1.e+10;
    gPar[kHelixCx] = m_cx;
    gPar[kHelixCy] = m_cy;
    gPar[kHelixZ0] = m_z0;
    gPar[kHelixR] = m_r;
    gPar[kHelixDz] = m_dz;

    Int_t flag=0;
    for(Int_t j=0; j<n; ++j){
      if(j==ihit) continue; //exclude the ith hit
      TPCLTrackHit *hitp = m_hit_array[j];
      TVector3 pos = hitp->GetLocalHitPos();
      Int_t layer = hitp->GetLayer();
      Double_t padTheta = hitp->GetPadTheta();
      std::vector<Double_t> resparam = hitp->GetResolutionParams();
      gHitPos[flag] = pos;
      gLayer[flag] = layer;
      gPadTheta[flag] = padTheta;
      gHelixTheta[flag] = hitp->GetTheta();
      gResParam[flag] = resparam;
      flag++;
    } //j

    //Helix fitting
    Bool_t vetoBadClusters = true;
    Bool_t exclusive = true;
    HelixFit(m_isBeam, vetoBadClusters, exclusive);
    m_chisqr_exclusive[ihit] = gChisqr;
    m_cx_exclusive[ihit] = gPar[kHelixCx];
    m_cy_exclusive[ihit] = gPar[kHelixCy];
    m_z0_exclusive[ihit] = gPar[kHelixZ0];
    m_r_exclusive[ihit]  = gPar[kHelixR];
    m_dz_exclusive[ihit] = gPar[kHelixDz];
    m_t_exclusive[ihit] = CalcThetaExclusive(ihit);

#if DebugDisp
    std::cout<<"exclusive "<<ihit<<" th chisqr : "<<m_chisqr_exclusive[ihit]<<
      " par : "<<m_cx_exclusive[ihit]<<
      " "<<m_cy_exclusive[ihit]<<
      " "<<m_z0_exclusive[ihit]<<
      " "<<m_r_exclusive[ihit]<<
      " "<<m_dz_exclusive[ihit]<<
      " "<<m_t_exclusive[ihit]<<std::endl;
#endif
  } //ihit
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::Print(const TString& arg, Bool_t print_allhits) const
{
  TString tracksize = Form(" #clusters = %d", (Int_t)m_hit_array.size());
  TString trackpar = Form(", params cx:%f cy:%f z0:%f r:%f dz:%f", m_cx, m_cy, m_z0, m_r, m_dz);
  std::cout<<arg.Data()<<std::endl;
  std::cout<<"Track info : "<<tracksize.Data()<<trackpar.Data()<<std::endl;
  if(print_allhits){
    for(std::size_t i=0; i<m_hit_array.size(); ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      if( !hitp ) continue;
      hitp->Print();
      std::cout<<"theta : "<<m_hit_t[i]<<std::endl;
    }
  }
}

//______________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetdEdx(Double_t truncatedMeanRatio)
{

  std::vector<Double_t> dEdx_vect;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    Double_t clde = hitp->GetDe();
    Double_t pathHit = hitp->GetPathHelix();
    Double_t dEdx_cor = clde/pathHit;
    dEdx_vect.push_back(dEdx_cor);
  }

  Double_t dEdx = 0.;
  std::sort(dEdx_vect.begin(), dEdx_vect.end());
  Int_t n_truncated = static_cast<Int_t>(dEdx_vect.size()*truncatedMeanRatio);
  if(n_truncated<=0 || n_truncated>static_cast<Int_t>(dEdx_vect.size())){
    hddaq::cerr << "#W " << FUNC_NAME << " "
                << "invalid n_truncated (n_truncated=" << n_truncated
                << ", nhit=" << dEdx_vect.size()
                << ", truncatedMeanRatio=" << truncatedMeanRatio << ")" << std::endl;
    return dEdx;
  }
  for( Int_t ih=0; ih<n_truncated; ++ih ){
    dEdx += dEdx_vect[ih];
  }
  dEdx /= static_cast<Double_t>(n_truncated);

  return dEdx;
}

//______________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTrackdE()
{

  Double_t dE = 0;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    dE += hitp->GetDe();
  }
  return dE;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::AddVPHit(TVector3 vp)
{
  m_hit_order.push_back(m_vp.size());
  m_vp.push_back(vp);
  m_hit_t.push_back(TMath::QuietNaN());

}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoVPFit()
{

  const std::size_t n = m_vp.size();
  gMultiLoop = false;
  gNumOfHits = n;
  gHitPos.clear();
  gHitPos.resize(n);
  gRes.clear();
  gRes.resize(n);
  gLayer.clear();
  gLayer.resize(n);
  gPadTheta.clear();
  gPadTheta.resize(n);
  gResParam.clear();
  gResParam.resize(n);
  gHelixTheta.clear();
  gHelixTheta.resize(n);
  gBadHits = 0;
  gMinuitStatus = 0;

  //for pre circle fitting
  Double_t xp[n]; Double_t yp[n];
  for(std::size_t i=0; i<n; ++i){
    gHitPos[i] = m_vp[i];
    xp[i] = -m_vp[i].X();
    yp[i] = m_vp[i].Z() - tpc::Z_TARGET;
    gRes[i] = TVector3(1., 1., 1.); //dummy values
  }

  Double_t par_circ[3]={0};
  //1st circle fitting to get initial parameters
  if(CircleFit(xp, yp, n, &par_circ[0], &par_circ[1], &par_circ[2])<0) return false;
  gPar[kHelixCx] = par_circ[0];
  gPar[kHelixCy] = par_circ[1];
  gPar[kHelixZ0] = 0.;
  gPar[kHelixR]  = par_circ[2];
  gPar[kHelixDz] = 0.;
  SetParam(gPar);

  Bool_t theta_flip = false;
  Double_t prev_theta = 0.; Double_t theta0 = 0;
  m_min_t = 9999; m_max_t = -9999;
  for(std::size_t i=0; i<n; ++i){
    Double_t theta = TMath::ATan2(m_vp[i].Z() - tpc::Z_TARGET - gPar[kHelixCy], -m_vp[i].X() - gPar[kHelixCx]);
    //Check ATan2 function's theta flip (-pi ~ pi)
    if(i==0) theta0 = theta;
    if(TMath::Abs(prev_theta - theta) > TMath::Pi()) theta_flip = true;
    if(theta_flip){
      if(theta0>0 && theta<0) theta += 2.*TMath::Pi();
      if(theta0<0 && theta>0) theta -= 2.*TMath::Pi();
    }
    prev_theta = theta;
    if(theta < m_min_t) m_min_t = theta;
    if(theta > m_max_t) m_max_t = theta;
    gHelixTheta[i] = theta;
    m_hit_t[i] = theta;
  }
  m_is_theta_calculated = true;

#if 0 //Optional
  Int_t MaxBinY[3];
  tpc::HoughTransformLineYTheta(gHitPos, MaxBinY, gPar, 1000.);
#endif
  if(!StraightLineFit()) return false;
  SetParam(gPar);
  m_minuit = gMinuitStatus;

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " RK helix params", m_cx, m_cy, m_z0, m_r, m_dz);
  for(std::size_t i=0; i<n; ++i){
    TVector3 tmp = GlobalPosition(gPar, m_hit_t[i]);
    TVector3 diff = m_vp[i] - tmp;
    std::cout<<FUNC_NAME+" VP - reconstructed VP : "<<diff.Mag()<<" mm"<<std::endl;
  }
#endif

  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ResidualCheck(TVector3 pos, Double_t xzwindow, Double_t ywindow, Double_t &resi)
{

  Bool_t status = false;
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  Double_t theta = EvalTheta(par, pos, m_min_t - 0.5*TMath::Pi(), m_max_t + 0.5*TMath::Pi());
  TVector3 fittmp = GlobalPosition(par, theta);
  TVector3 d = pos - fittmp;
  resi = d.Mag();
  Double_t xz_resi = TMath::Sqrt(d.x()*d.x()+d.z()*d.z());
  Double_t y_resi = TMath::Sqrt(d.y()*d.y());

  if(xz_resi<xzwindow && y_resi<ywindow) status = true;
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ResidualCheck(TVector3 pos, Double_t xzwindow, Double_t ywindow)
{
  Double_t resi;
  return ResidualCheck(pos, xzwindow, ywindow, resi);
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SortHitOrder()
{

  WarnValidateCalcHelixTheta(*this, FUNC_NAME);

  m_hit_order.clear();
  for(std::size_t i=0; i<m_hit_t.size(); ++i){
    m_hit_order.push_back(i);
  }
  for(std::size_t i=0; i<m_hit_order.size(); ++i){
    std::sort(m_hit_order.begin(), m_hit_order.end(), CompareTheta);
  }
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DetermineCharge()
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

  Double_t minlayer_t = 0., maxlayer_t = 0.;
  Int_t minlayer = 33, maxlayer = -1;
  if(m_hit_array.size()==0) return false;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    Int_t layer = hitp->GetLayer();
    if(layer<minlayer){
      minlayer = layer;
      minlayer_t = m_hit_t[i];
    }
    if(layer>maxlayer){
      maxlayer = layer;
      maxlayer_t = m_hit_t[i];
    }
  }

  // For E72 conditions: the magnetic field direction is along the Z-axis (u = (0, 0, 1))
  if(minlayer_t<maxlayer_t) m_charge = -1;
  else m_charge = 1;

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  m_edgepoint = GlobalPosition(par, maxlayer_t);

  if(m_isK18) m_charge = -1; //for K1.8 tracking.
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::VertexAtTarget()
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

  Bool_t status = false;
  if(m_closedist.Mag() < tpc::TARGET_RADIUS) status = true;
  return status;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsBackward()
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

  // track is starting from the target
  if(TMath::Abs(TMath::Hypot(m_cx, m_cy) - m_r) > tpc::TARGET_RADIUS) return false;

  // Track exist before the target position
  // if(m_edgepoint.z() > tpc::Z_TARGET) return false;
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 start_point  = GlobalPosition(par, m_min_t);
  TVector3 end_point    = GlobalPosition(par, m_max_t);
  if(start_point.z() > tpc::Z_TARGET || end_point.z() > tpc::Z_TARGET) return false;

  // upstream end of the track is within window
  // Solve the circle equation at Z = -250; this is the X discriminant (must be >= 0).
  const Double_t y_local = -250. - tpc::Z_TARGET;
  const Double_t x_disc  = m_r*m_r - TMath::Sq(y_local - m_cy);
  if(x_disc < 0) return false;

  // Z=-250 has two circle branches; pick the one continuous with the measured arc.
  Double_t theta_ref = (GlobalPosition(par, m_max_t).z() < GlobalPosition(par, m_min_t).z())
    ? m_max_t : m_min_t;
  const std::size_t n_hit = m_hit_array.size();
  // Anchor to the most upstream hit theta when hits are available.
  if(n_hit > 0 && gHelixTheta.size() == n_hit){
    Double_t z_upstream = m_hit_array[0]->GetLocalHitPos().z();
    std::size_t i_upstream = 0;
    for(std::size_t i=1; i<n_hit; ++i){
      const Double_t z = m_hit_array[i]->GetLocalHitPos().z();
      if(z < z_upstream){
        z_upstream = z;
        i_upstream = i;
      }
    }
    theta_ref = gHelixTheta[i_upstream];
  }

  const Double_t s = TMath::Sqrt(x_disc);
  const Double_t theta_plus = TMath::ATan2(y_local - m_cy, s);
  const Double_t theta_minus = TMath::ATan2(y_local - m_cy, -s);
  // Shortest angular separation from theta_ref in (-pi, pi].
  Double_t d_theta_plus = theta_plus - theta_ref;
  while(d_theta_plus > TMath::Pi()) d_theta_plus -= 2.*TMath::Pi();
  while(d_theta_plus < -TMath::Pi()) d_theta_plus += 2.*TMath::Pi();
  Double_t d_theta_minus = theta_minus - theta_ref;
  while(d_theta_minus > TMath::Pi()) d_theta_minus -= 2.*TMath::Pi();
  while(d_theta_minus < -TMath::Pi()) d_theta_minus += 2.*TMath::Pi();
  const Double_t theta_cut = (TMath::Abs(d_theta_minus) < TMath::Abs(d_theta_plus))
    ? theta_minus : theta_plus;

  const TVector3 pos_cut = GlobalPosition(par, theta_cut);
  const Double_t extrap_point_x = TMath::Abs(pos_cut.x()); // Abs(X) at Z=-250 on measured branch
  if(extrap_point_x > BackwardMaxAbsX) return false;

  // Tag straight accidental-beam backward tracks (eases helix fit).
  if(TMath::Abs(m_dz) < BeamLikeMaxAbsDz &&
     IsBeamLikeHit(GlobalToLocal(start_point)) &&
     IsBeamLikeHit(GlobalToLocal(end_point)))
    SetIsBeam();

  return true;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::IsMultiLoop()
{
  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return;
  //High pT spiral-like track
  //Pitch > NSigma * y_resolution
  //Radius < 250.mm (TPC radius)
  //loop is larger then half circle
  Double_t pitch = TMath::Abs(2.*TMath::Pi()*m_r*m_dz);
  if(!m_is_multiloop){
    if(pitch > ThetaNSigma*GetResolutionY(0) &&
       m_r < 250. &&
       (m_max_t - m_min_t) > TMath::Pi())
      m_is_multiloop = true;
#if DebugDisp
    if(m_is_multiloop) std::cout<< " Multi-loop track!!"<<std::endl;
#endif
  }
  else{
    if(pitch < ThetaNSigma*GetResolutionY(0) || m_r > 250.) m_is_multiloop = false;
#if DebugDisp
    if(!m_is_multiloop) std::cout<< " Not Multi-loop track!!"<<std::endl;
#endif
  }
  gMultiLoop = m_is_multiloop;
}

//______________________________________________________________________________
Int_t
TPCLocalTrackHelix::FinalizeTrack(Int_t &delete_hit)
{
  // Returns: 0 = OK, >0 = bad hit count, -1 = fatal

  if(m_hit_array.size()==0) return -1;

  if(!ValidateEqualSizes(FUNC_NAME, m_hit_array.size(),
                         "m_hit_order", m_hit_order.size(),
                         "m_hit_t", m_hit_t.size()))
    return -1;

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return -1;

  Int_t false_layer = 0;
  Double_t max_residual = -100.;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    Double_t resi = 0.;
    if(!ResidualCheck(i, resi)) ++false_layer;
    if(max_residual<resi){
      max_residual = resi;
      delete_hit = i;
    }
  }

  SortHitOrder();

  // Global coordinates:
  // x = x_0 + r * cos(theta)
  // y = y_0 + r * sin(theta)
  // z = z_0 + m_dz * r * theta
  // dl = sqrt( dx^2 + dy^2 + dz^2 ) = sqrt( r^2 + r^2*m_dz^2 ) * d(theta)
  const Double_t dtheta = m_max_t - m_min_t;
  m_path = dtheta*TMath::Hypot(m_r, m_r*m_dz);
  m_transverse_path = dtheta*m_r;
  if(false_layer!=0 || m_chisqr > MaxChisqr) return false_layer;

  if(!DetermineCharge()) return -1;
  m_mom0 = CalcHelixMom(gPar, 0.);
  IsMultiLoop();

#if DebugDisp
  DebugHelixPar(FUNC_NAME + Form(" FinalizeTrack\n    chisqr: %g\n    charge: %d\n    path: %g\n    |m_min_t-m_max_t|: %g",
                                 m_chisqr, m_charge, m_path, TMath::Abs(m_min_t - m_max_t)),
                m_cx, m_cy, m_z0, m_r, m_dz);
  if(m_path>550.) std::cout<<FUNC_NAME+" too long track!!! : m_path="<<m_path<<std::endl;
#endif

  return false_layer;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::ConvertParam(Double_t *linear_par)
{

  gPar[kHelixCx] = 0;
  gPar[kHelixCy] = 0;
  gPar[kHelixZ0] = linear_par[1];
  gPar[kHelixR] = 0;
  gPar[kHelixDz] = linear_par[3];

  gMultiLoop = false;
  gMomConstraint = false; //No momentum constraint

  if(!DoPreFit(gPar)) return false;

  //Vtx in the target, need to check whether two tracks are merged or not
  if(!SeparateTracksAtTarget()) return false;

#if DebugDisp
  DebugHelixPar(FUNC_NAME + Form(" Converted track\n    chisqr: %g\n    theta min: %g\n    theta max: %g\n    hit_t size: %zu",
                                 gChisqr, m_min_t, m_max_t, m_hit_t.size()), gPar);
#endif

  return true;
}

//______________________________________________________________________________
//If two tracks are merged at the target, separate them and recalculate params.
Bool_t
TPCLocalTrackHelix::SeparateTracksAtTarget()
{
  // If the target is at the middle of the inintial track, 
  // that track is separated with different side flag. (Two tracks are reconized as a single track)
  // Please see the Side() function.

  Bool_t status = true;
  m_vtxflag = 0;

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1.);
  if(BeamThroughTPC || m_isAccidental==1){
    return status;
  }

  WarnValidateCalcHelixTheta(*this, FUNC_NAME);

  //High pT spiral-like track
  if(m_is_multiloop && (m_max_t - m_min_t) > 2.*TMath::Pi()) return status;

  CalcClosestDistTgt(); //Distance between the target & the track
  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  if(VertexAtTarget()){ //If track is crossing the target.
    //Exclude the beam hit from the scattered track of the other scattered track's hit.
    const std::size_t n = m_hit_array.size();
    for(std::size_t i=0; i<n; ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      TVector3 pos = hitp->GetLocalHitPos();
      if(Side(pos)==1) side1_hits.push_back(i);
      else if(Side(pos)==-1) side2_hits.push_back(i);
    }
  }
  else{ //If track is not crossing the target.
    SortHitOrder();

    const Double_t p_t = HelixRToPt(m_r); // [GeV/c]

    //Exclude the beam hit from the scattered track.
    Bool_t flag = false;
    TVector3 prev_pos;
    Bool_t prev_is_beam_hit = false;
    Bool_t is_beam_hit = false;
    const std::size_t n = m_hit_array.size();
    for(Int_t i=0; i<n; ++i){
      Int_t id = m_hit_order[i];
      TPCLTrackHit *hitp = m_hit_array[id];
      TVector3 pos = hitp->GetLocalHitPos();
      is_beam_hit = IsBeamLikeHit(pos);
      TVector3 gap = pos - prev_pos;
#if DebugDisp
      std::cout << FUNC_NAME
                << " i=" << i << " gap=" << gap.Mag()
                << " pos=(" << pos.x() << "," << pos.y() << "," << pos.z() << ")"
                << std::endl;
#endif
      // Skip separation for nearly straight + high-pT (low curvature) tracks.
      if((m_dz > MinSlopeForBeamScatterSep || p_t < MaxPtForBeamScatterSep)
         && i!=0 && gap.Mag() > MinGapForBeamScatterSep
         && prev_is_beam_hit != is_beam_hit) flag = true;
      if(!flag) side1_hits.push_back(id);
      else side2_hits.push_back(id);
      prev_pos = pos;
      prev_is_beam_hit = is_beam_hit;
    }
  }

  if(side1_hits.size()==0 || side2_hits.size()==0) return status;
  else{ //exclude the shorter side
    m_is_fitted = false; //Need to do minimization again
    m_n_iteration = 0;

    if(side1_hits.size() >= side2_hits.size()) EraseHits(side2_hits);
    else EraseHits(side1_hits);

    Double_t par[5]={0.};
    status = DoPreFit(par);
    if(status){
      CalcClosestDistTgt();
      TVector3 pos = m_hit_array[0]->GetLocalHitPos();
      if(VertexAtTarget()) m_vtxflag = Side(pos);
      else m_vtxflag = 0;
    }
  }

#if DebugDisp
  if(!status) std::cout<<FUNC_NAME+" Separated track is not good for tracking"<<std::endl;
  std::cout<<FUNC_NAME+" Close distance from the target to the track : "<<m_closedist.Mag()<<std::endl;
#endif

  return status;
}

//______________________________________________________________________________
//If two tracks are merged at the target, separate them and recalculate params.
Bool_t
TPCLocalTrackHelix::SeparateClustersWithGap()
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;
  SortHitOrder();

  //Check a gap between clusters.
  Bool_t flag = false;
  TVector3 prev_pos;
  std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
  const std::size_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    Int_t id = m_hit_order[i];
    TPCLTrackHit *hitp = m_hit_array[id];
    TVector3 pos = hitp -> GetLocalHitPos();
    TVector3 gap = pos - prev_pos;
#if DebugDisp
    std::cout << FUNC_NAME
              << " i=" << i << " gap=" << gap.Mag()
              << " pos=(" << pos.x() << "," << pos.y() << "," << pos.z() << ")"
              << std::endl;
#endif
    if(i!=0 && gap.Mag() > MaxGapBtwClusters) flag = true;
    if(!flag) side1_hits.push_back(id);
    else side2_hits.push_back(id);
    prev_pos = pos;
  }

  if(side1_hits.size()!=0 && side2_hits.size()!=0){
    /*
    std::cout<<"size "<< side1_hits.size()<<" "<<side2_hits.size()<<std::endl;
    for(Int_t i=0; i<side1_hits.size(); ++i){
      std::cout<<i<<"/"<<side1_hits.size()<<" "<<m_hit_array[side1_hits[i]] -> GetLocalHitPos()<<std::endl;
    }
    for(Int_t i=0; i<side2_hits.size(); ++i){
      std::cout<<i<<"/"<<side2_hits.size()<<" "<<m_hit_array[side2_hits[i]] -> GetLocalHitPos()<<std::endl;
    }
    */
    if(side1_hits.size() >= side2_hits.size()) EraseHits(side2_hits);
    else EraseHits(side1_hits);
  }

  gHitPos.clear();
  return flag;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::TestMergedTrack()
{

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " longer track's Helix params", m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  //Initialization of helix params
  gPar[kHelixCx] = m_cx;
  gPar[kHelixCy] = m_cy;
  gPar[kHelixZ0] = m_z0;
  gPar[kHelixR]  = m_r;
  gPar[kHelixDz] = m_dz;

  Bool_t vetoBadClusters = false;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gBadHits = 0;
  gMinuitStatus = 0;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gRes.clear();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
    Int_t layer = hitp->GetLayer();
    gLayer.push_back(layer);
    Double_t padTheta = hitp->GetPadTheta();
    gPadTheta.push_back(padTheta);
    std::vector<Double_t> resparam = hitp->GetResolutionParams();
    gResParam.push_back(resparam);
    Double_t dummy_theta = 0;
    TVector3 res = CalcResolution(gPar, layer, pos, padTheta, dummy_theta, resparam, vetoBadClusters);
    gRes.push_back(res);
  }
  CalcHelixTheta();

  Int_t ndf;
  gChisqr = CalcChi2(gPar, ndf, vetoBadClusters);

  //pre fitting
  if(!DoPreFit(gPar)) return false;

  //Helix fitting
  if(!DoHelixFit(gPar, vetoBadClusters)) return false;
#if DebugDisp
  DebugHelixPar(FUNC_NAME + Form(" Helix fitting of the merged track\n    n_iteration: %d\n    chisqr: %g\n    # of clusters: %d",
                                 m_n_iteration, gChisqr, GetNHit()), gPar);
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
  if(false_layer < 0) return false;
#if DebugDisp
  std::cout << FUNC_NAME << " # of bad clusters : " << false_layer << std::endl;
#endif

  // Empirical merge pre-check: bad merges fail clearly, good ones pass cleanly.
  // Not rigorously optimized, but OK in practice.
  if(false_layer > 3 && (Double_t) false_layer/GetNHit() > 0.2) return false;

  //Check whether the track passing the target or not
  CalcClosestDistTgt(); //Distance between the target & the track
  CheckIsAccidental(); //check whether it is accidental beam

  if(VertexAtTarget()){
    //case1. The merged track is not crossing the target. (no problem)
    //case2. The merged track is an accidental beam crossing the target. (no problem)
    //case3. The merged track is wrongly merged two individual scattered tracks near the target <- veto

    std::vector<Int_t> side1_hits; std::vector<Int_t> side2_hits;
    for(Int_t i=0; i<n; ++i){
      Double_t resi;
      if(!ResidualCheck(i, resi)) continue;
      TPCLTrackHit *hitp = m_hit_array[i];
      TVector3 pos = hitp->GetLocalHitPos();
      if(Side(pos)==1) side1_hits.push_back(i);
      else if(Side(pos)==-1) side2_hits.push_back(i);
    }

    if(side1_hits.size() > 0 && side2_hits.size() > 0){
      // The caller validates a geometry-selected fragment pair with a full
      // refit. During that test only, target crossing is allowed before its
      // final beam classification is assigned.
      if(m_isBeam!=1 && !m_allow_target_crossing_merge) return false;
      /*
      if(TMath::Abs(m_closedist.x())<15. &&
         TMath::Abs(m_closedist.y())<10. &&
         TMath::Abs(m_closedist.z())<10.){
        if(m_isAccidental!=1 || m_isBeam!=1) return false; //(case3) if it is not accidental beam
      }
      */
    }
  }

  //Set the vtxflag 0 to make SeparateTracksAtTarget() off.
  m_vtxflag=0;

#if IterativeResolution
  vetoBadClusters = true;

  //Now excluding bad clusters and fitting again.
  DoHelixFit(gPar, vetoBadClusters);
#endif
  return true;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::DoFitTrackwVertex(TVector3 vertex_pos, TVector3 vertex_res)
{

  Bool_t status = false;
  if(!IsGoodForTracking()) return status;

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " Before vertex constraint fitting, track's Helix params",
                m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  //Initialization of helix params
  gPar[kHelixCx] = m_cx;
  gPar[kHelixCy] = m_cy;
  gPar[kHelixZ0] = m_z0;
  gPar[kHelixR] = m_r;
  gPar[kHelixDz] = m_dz;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();

  //Vertex information
  gVertex = vertex_pos;
  gVertexRes = vertex_res;

  gMultiLoop = m_is_multiloop;
  gNumOfHits = 0;
  gBadHits = 0;
  gMinuitStatus = 0;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  gRes.clear();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    gNumOfHits++;
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
    Int_t layer = hitp->GetLayer();
    gLayer.push_back(layer);
    Double_t padTheta = hitp->GetPadTheta();
    gPadTheta.push_back(padTheta);
    std::vector<Double_t> resparam = hitp->GetResolutionParams();
    gResParam.push_back(resparam);
    const TVector3& res = hitp->GetResolutionVect();
    gRes.push_back(res);
  }

  CalcHelixTheta();
  if(HelixFitwVertex()){
    status = true;
    SetParam(gPar);
    CalcHelixTheta();
    Bool_t vetoBadClusters = true;
    Int_t ndf;
    m_chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
    m_minuit = gMinuitStatus;
    Double_t window = HelixThetaSearchWindow(gPar[kHelixR]);
    for(Int_t i=0; i<gNumOfHits; ++i){
      Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
      m_hit_t[i] = theta;
      gHelixTheta[i] = theta;
    }

    Int_t delete_hit = -1;
    if(FinalizeTrack(delete_hit) < 0) status = false;
  }

  return status;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::RecalcTrack()
{
  if (m_hit_array.empty() || m_hit_t.empty()) {
    return;
  }
  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  m_is_theta_calculated = true;
  m_is_calculated = false;
  m_mom0 = CalcHelixMom(par, 0.);
  Calculate();
  m_is_fitted = true;
  m_min_t = m_hit_t[0];
  m_max_t = m_hit_t[m_hit_t.size() - 1];
  CalcClosestDistTgt();
  const Double_t dtheta = m_max_t - m_min_t;
  m_path = dtheta * TMath::Hypot(m_r, m_r * m_dz);
  m_transverse_path = dtheta * m_r;
  IsMultiLoop();
}

//_____________________________________________________________________________
void
TPCLocalTrackHelix::CheckIsAccidental()
{

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return;
  //if(!m_is_fitted) return;
  if(m_is_multiloop) return; //A multiloop track is obviously not accidental beam

  static const Bool_t BeamThroughTPC = (gUser.GetParameter("BeamThroughTPC") == 1.);
  if(BeamThroughTPC) return;
  m_isAccidental = 0;

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 start_point = GlobalPosition(par, m_min_t);
  TVector3 end_point   = GlobalPosition(par, m_max_t);

  const Double_t diff_z = TMath::Abs(start_point.z() - end_point.z());
  if(diff_z < AccidentalMinDiffZ) return;

  const Double_t p_t = HelixRToPt(m_r); // [GeV/c]
  if(p_t < BeamMom - BeamMomOffset) return;

  Int_t nhit_upstream_tgt = 0;
  Int_t nhit_downstream_tgt = 0;
  const std::size_t n = m_hit_array.size();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    if(pos.z() < tpc::Z_TARGET) nhit_upstream_tgt++;
    if(pos.z() > tpc::Z_TARGET) nhit_downstream_tgt++;
  }

  if(nhit_upstream_tgt >= 1 && nhit_downstream_tgt >= 5){
    m_isAccidental = 1;
    if(TMath::Abs(m_dz) < BeamLikeMaxAbsDz) m_isBeam = 1;
  }
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrackHelix::IsBeamLikeHit(const TVector3& pos) const
{
  const Bool_t in_x = (pos.x() > BeamLikeXMin) && (pos.x() < BeamLikeXMax);
  const Bool_t in_y = (pos.y() > BeamLikeYMin) && (pos.y() < BeamLikeYMax);
  const Bool_t upstream_of_target = (pos.z() < tpc::Z_TARGET);
  return in_x && in_y && upstream_of_target;
}

//______________________________________________________________________________
Bool_t
TPCLocalTrackHelix::TestInvertCharge()
{

  if(m_is_multiloop || m_isBeam==1 || m_isK18==1 || m_isAccidental==1) return false;

  if(!ValidateCalcHelixTheta(*this, FUNC_NAME)) return false;

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " track for testing", m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  //Initialization of helix params
  gPar[kHelixCx] = m_cx;
  gPar[kHelixCy] = m_cy;
  gPar[kHelixZ0] = m_z0;
  gPar[kHelixR]  = m_r;
  gPar[kHelixDz] = m_dz;

  Bool_t vetoBadClusters = false;
  Int_t prev_charge = m_charge;

  DeleteNullHit();
  const Int_t n = m_hit_array.size();
  if(!IsGoodForTracking()) return false;
  gMultiLoop = m_is_multiloop;
  gNumOfHits = n;
  gBadHits = 0;
  gMinuitStatus = 0;
  gHitPos.clear();
  gLayer.clear();
  gPadTheta.clear();
  gResParam.clear();
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    gHitPos.push_back(pos);
    Int_t layer = hitp->GetLayer();
    gLayer.push_back(layer);
    Double_t padTheta = hitp->GetPadTheta();
    gPadTheta.push_back(padTheta);
    std::vector<Double_t> resparam = hitp->GetResolutionParams();
    gResParam.push_back(resparam);
  }

  TVector3 localpos_start = GlobalToLocal(gHitPos[0]);
  TVector3 localpos_end   = GlobalToLocal(gHitPos[n-1]);
  TVector3 localpos_avg   = 0.5*(localpos_start + localpos_end);
  Double_t x_avg = localpos_avg.x();
  Double_t y_avg = localpos_avg.y();
  Double_t z_avg = localpos_avg.z();
  // Mirror the helix center across pos_avg(A):
  // C = (m_x, m_y), A = (C + C')/2 -> C' = 2*A - C
  gPar[kHelixCx] = 2.*x_avg - m_cx;
  gPar[kHelixCy] = 2.*y_avg - m_cy;
  Double_t tmp_theta = TMath::ATan2(y_avg - gPar[kHelixCy], x_avg - gPar[kHelixCx]);
  gPar[kHelixZ0] = z_avg - tmp_theta*m_r*(-m_dz);
  gPar[kHelixR]  = m_r;
  gPar[kHelixDz] = -m_dz;
  SetParam(gPar);
  CalcHelixTheta();

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " before helixfit for inverting", m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  gRes.clear();
  for(Int_t i=0; i<n; ++i){
    TVector3 res = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    gRes.push_back(res);
  }

  Bool_t status = false;
  if(HelixFitInvertCharge()){
    status = true;
    SetParam(gPar);
    CalcHelixTheta();

    Double_t window = HelixThetaSearchWindow(gPar[kHelixR]);
    for(Int_t i=0; i<gNumOfHits; ++i){
      Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
      m_hit_t[i] = theta;
      gHelixTheta[i] = theta;
      gRes[i] = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
    }

    Int_t ndf;
    m_chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
    m_minuit = gMinuitStatus;
  }
  else return false;

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " after helixfit for inverting", m_cx, m_cy, m_z0, m_r, m_dz);
#endif

  Int_t delete_hit = -1;
  Int_t false_layer = FinalizeTrack(delete_hit);
  if(false_layer < 0) return false;
  if(!DetermineCharge() || prev_charge == m_charge) return false; //charge is not converted

#if DebugDisp
  std::cout << FUNC_NAME << " # of bad clusters : " << false_layer << std::endl;
#endif

  if(false_layer > 0){
    Int_t maxloop = 0;
    while(false_layer > 0 && maxloop < 30){
      if(!HelixFitInvertCharge()) return false;
      else{
        Double_t window = HelixThetaSearchWindow(gPar[kHelixR]);
        for(Int_t i=0; i<gNumOfHits; ++i){
          Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
          m_hit_t[i] = theta;
          gHelixTheta[i] = theta;
          gRes[i] = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);
        }

        Int_t ndf;
        Double_t chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
        if(chisqr >= m_chisqr) break;
        SetParam(gPar);
        CalcHelixTheta();
        m_chisqr = CalcChi2(gPar, ndf, vetoBadClusters);
        m_minuit = gMinuitStatus;
        false_layer = FinalizeTrack(delete_hit);
      }
      maxloop++;
    } //while

    if(false_layer < 0 || false_layer > 2) return false;
    if(!DetermineCharge() || prev_charge == m_charge) return false;
  }

  //Check whether the track passing the target or not
  CalcClosestDistTgt(); //Distance between the target & the track
  CheckIsAccidental(); //check whether it is accidental beam
  VertexAtTarget();

#if IterativeResolution
  vetoBadClusters = true;

  for(Int_t i=0; i<gNumOfHits; ++i) gRes[i] = CalcResolution(gPar, gLayer[i], gHitPos[i], gPadTheta[i], gHelixTheta[i], gResParam[i], vetoBadClusters);

  //Now excluding bad clusters and fitting again.
  HelixFitInvertCharge();
  SetParam(gPar);
  CalcHelixTheta();
  m_chisqr = gChisqr;
  m_minuit = gMinuitStatus;
  Double_t window = HelixThetaSearchWindow(gPar[kHelixR]);
  for(Int_t i=0; i<gNumOfHits; ++i){
    Double_t theta = EvalTheta(gPar, gHitPos[i], gHelixTheta[i] - 0.5*window, gHelixTheta[i] + 0.5*window);
    m_hit_t[i] = theta;
    gHelixTheta[i] = theta;
  }

#if DebugDisp
  DebugHelixPar(FUNC_NAME + " after refitting helixfit for inverting", m_cx, m_cy, m_z0, m_r, m_dz);
  std::cout<<"prev current charge "<<prev_charge<<" "<<m_charge<<std::endl;
#endif
#endif

  return status;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetClosestPositionTgt()
{
  TVector3 tgt(0., 0., tpc::Z_TARGET);
  TVector3 pos = m_closedist + tgt;
  return pos;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetClosestPositionTgtXZ()
{
  TVector3 tgt(0., 0., tpc::Z_TARGET);
  TVector3 pos = m_closedistXZ + tgt;
  return pos;
}

//Functions For Kinematic Fit

//Functions for momentum resolution
//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumResolutionVectT(Double_t t, Double_t MomScale, Double_t PhiScale, Double_t dZScale){

  /*
    For a Multivaraible function F(M), the Covariance is given as V(F) = JV(M)J^T, where J = dFdM is a Jacobian matrix.
    We will calculate the momentum resolution from pt, theta, and dZ resolution.
  */

  Double_t p_t = HelixRToPt(m_r); // [GeV/c]
  //Double_t pz = p_t*(cos(t));
  //Double_t py = p_t*m_dz;
  //Double_t px = p_t*(sin(t));

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 calcmom = CalcHelixMom(par, t);
  Double_t ddZ = dZScale * GetdZResolution();
  Double_t dt = PhiScale * GetTransverseAngularResolution(t);
  Double_t dp_t = MomScale * GetTransverseMomentumResolution(); // [GeV/c]
  Double_t Vp_t = dp_t*dp_t;
  Double_t Vt = dt*dt;
  Double_t VdZ = ddZ*ddZ;
  /*
    p_x = p_t*sin(t);
    p_y = p_t*dZ;
    p_z = p_t*cos(t);

    Cov(px,pz,py) = J V(p_t,t,dZ)J^T
    **Note! Cov(pt,t) = res_pt*res_t, because res_t ~ res_pt!

             f
		 J = x| dx/df|

		    p_t    t     dZ
    px| sin  pt*cos   0 |;
  J=py|  dZ    0     pt |;
    pz| cos -pt*sin   0 |;
 */

  Double_t cov_pt_t = MomScale * PhiScale * GetTransverseMomentumAngularCovariance(t);
  Double_t V[3][3] =
    { Vp_t,   cov_pt_t,    0,
      cov_pt_t, Vt,        0,
      0,         0,       VdZ};
  Double_t J[3][3] =
    { sin(t),   p_t*cos(t), 0,
      m_dz,     0,          p_t,
      cos(t), -p_t*sin(t), 0};

  Double_t JT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      JT[row][col] = J[col][row];
    }
  }
  Double_t VJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	VJT[row][col]+=V[row][itr]*JT[itr][col];
      }
    }
  }
  Double_t JVJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	JVJT[row][col]+=J[row][itr]*VJT[itr][col];
      }
    }
  }
  Double_t Vpx = JVJT[0][0];
  Double_t Vpy = JVJT[1][1];
  Double_t Vpz = JVJT[2][2];
  if(Vpx<0 or Vpy<0 or Vpz < 0 or std::isnan(Vpx) or std::isnan(Vpy) or std::isnan(Vpz)){
    std::cout<<Form("MomVar = (%g,%g,%g)",Vpx,Vpy,Vpz)<<std::endl;
    std::cout<<Form("dPt,dt,ddZ = (%g,%g,%g)",dp_t,dt,ddZ)<<std::endl;
  }
  return TVector3(sqrt(Vpx), sqrt(Vpy), sqrt(Vpz));
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumResolutionVect(Int_t i, Double_t MomScale, Double_t PhiScale, Double_t dZScale){
  Double_t t = GetHitInOrder(i) -> GetTheta();
  return GetMomentumResolutionVectT(t, MomScale, PhiScale, dZScale);
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumResolutionVect(){

  Double_t t0 = GetHitInOrder(0) -> GetTheta();
  return GetMomentumResolutionVectT(t0, m_MomResScale, m_PhResScale, m_dZResScale);
}

//_____________________________________________________________________________
double
TPCLocalTrackHelix::GetTransverseMomentumAngularCovariance(Double_t t){


  if(t == -9999){
    t = GetHitInOrder(0)->GetTheta();
  }

  Double_t sign = 1;
  if(m_charge>0)sign = -1;
  Double_t dp_t =  GetTransverseMomentumResolution();
  Double_t dt = GetTransverseAngularResolution(t,0);
  return sign * dp_t * dt;
}

//_____________________________________________________________________________
double
TPCLocalTrackHelix::GetMomentumPitchAngleCovariance(){
  Double_t p_t = HelixRToPt(m_r); // [GeV/c]
  Double_t pitch = atan2(1,m_dz);//dYdZ angle, dZ= 0 -> should return pi/2
  Double_t res_pitch = GetThetaResolution();
  return p_t*cos(pitch)/sin(pitch)/sin(pitch)*res_pitch*res_pitch;

}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumCovarianceVectT(Double_t t, Double_t MomScale, Double_t PhiScale, Double_t dZScale){

  /*
    For a Multivaraible function F(M), the Covariance is given as V(F) = JV(M)J^T, where J = dFdM is a Jacobian matrix.
    We will calculate the momentum resolution from pt, theta, and dZ resolution.
  */
  Double_t p_t = HelixRToPt(m_r); // [GeV/c]
  //Double_t pz = p_t*(cos(t));
  //Double_t py = p_t*m_dz;
  //Double_t px = p_t*(sin(t));

  Double_t par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 calcmom = CalcHelixMom(par, t);
  Double_t ddZ = dZScale * GetdZResolution();
  Double_t dt = PhiScale * GetTransverseAngularResolution(t);
  Double_t dp_t = MomScale * GetTransverseMomentumResolution(); // [GeV/c]
  Double_t Vp_t = dp_t*dp_t;
  Double_t Vt = dt*dt;
  Double_t VdZ = ddZ*ddZ;
  /*
    p_x = p_t*sin(t);
    p_y = p_t*dZ;
    p_z = p_t*cos(t);

    Cov(px,pz,py) = J V(p_t,t,dZ)J^T
    **Note! Cov(pt,t) = res_pt*res_t, because res_t ~ res_pt!
  */
  Double_t cov_pt_t = MomScale*PhiScale*GetTransverseMomentumAngularCovariance(t);
  Double_t V[3][3] =
    { Vp_t,     cov_pt_t,     0,
      cov_pt_t, Vt,           0,
      0,        0,            VdZ};
  Double_t J[3][3] =
    { sin(t),   p_t*cos(t),   0,
      m_dz,     0,            p_t,
      cos(t),   -p_t*sin(t),  0};

  Double_t JT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      JT[row][col] = J[col][row];
    }
  }
  Double_t VJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	VJT[row][col]+=V[row][itr]*JT[itr][col];
      }
    }
  }
  Double_t JVJT[3][3]={0};
  for(Int_t row=0;row<3;++row){
    for(Int_t col=0;col<3;++col){
      for(Int_t itr=0;itr<3;++itr){
	JVJT[row][col]+=J[row][itr]*VJT[itr][col];
      }
    }
  }
  Double_t Cxy = JVJT[0][1];
  Double_t Cyz = JVJT[1][2];
  Double_t Czx = JVJT[2][0];
  return TVector3((Cxy), (Cyz), (Czx));
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumCovarianceVect(Int_t i, Double_t MomScale, Double_t PhiScale, Double_t dZScale){

  Double_t t = GetHitInOrder(i) -> GetTheta();
  return GetMomentumCovarianceVectT(t, MomScale, PhiScale, dZScale);
}

//_____________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetMomentumCovarianceVect(){

  Double_t t0 = GetHitInOrder(0) -> GetTheta();
  return GetMomentumCovarianceVectT(t0, m_MomResScale, m_PhResScale, m_dZResScale);

}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTransverseMomentumResolution(){
  Double_t B = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
  Double_t dt = abs(m_max_t - m_min_t);
  if(dt > 2*acos(-1) )dt = 2*acos(-1);

  Double_t L = 2 * sin(0.5*dt)*m_r;//String length, not Arc length
  L*=0.001;//mm->m;
  Double_t res = 0;
  Int_t nh = m_hit_array.size();
  for(Int_t ih=0;ih<m_hit_array.size();++ih){
    Int_t id = m_hit_order[ih];
    TPCLTrackHit *hitp = m_hit_array[id];
    auto ResV = hitp -> GetResolutionVect();
    Double_t res_T = hypot(ResV.X(),ResV.Z());
    if(!hitp->IsGoodForTracking()){
      nh--;
      continue;
    }
    res+=res_T*res_T;
  }
  Double_t p_t = HelixRToPt(m_r); // [GeV/c]
  if(nh<4) return p_t*0.1;
  res = sqrt(3./2) * sqrt(res / nh)* 0.001;//mm-> m
  Double_t dPOverP = p_t / (0.3*L*L*B)*sqrt(720./(nh+4))*res;
  if(std::isnan(dPOverP) || dPOverP < 0){
    std::cout<<Form("dPt error! nh = %d, res = %g",nh,res)<<std::endl;
  }
  return m_MomResScale*p_t*dPOverP;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTransverseAngularResolution(Double_t t, Double_t sig0){
  //Transverse Angle Definition: atan2(pz,px);
  Double_t p_t = HelixRToPt(m_r); // [GeV/c]
  Double_t dp = GetTransverseMomentumResolution(); // [GeV/c]
  Double_t dr = m_r * dp/p_t;
  Double_t t_avg = 0.5*(m_max_t + m_min_t);
  Double_t dt = (t-t_avg);
  if(dt >acos(-1))dt = acos(-1);
  Double_t path = m_r * abs(dt);

  return hypot(m_PhResScale*path * dr / m_r / m_r,sig0);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetTransverseAngularResolution(){

  Double_t t0 = GetHitInOrder(0) -> GetTheta();
  Double_t sig0 = 0.001;
  return GetTransverseAngularResolution(t0, sig0);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetdZResolution(){
  Double_t res2 = 0;
  Int_t nh = m_hit_array.size();
  for(Int_t ih=0;ih<m_hit_array.size();++ih){
    Int_t id = m_hit_order[ih];
    TPCLTrackHit *hitp = m_hit_array[id];
    auto ResV = hitp -> GetResolutionVect();
    Double_t res_Y = ResV.Y();
    if(!hitp->IsGoodForTracking()){
      nh--;
      continue;
    }
    res2 += res_Y*res_Y;
  }
  Double_t dt = abs(m_max_t - m_min_t);
  //  if(dt > 2*acos(-1)) dt = 2*acos(-1);
  Double_t path = m_r * dt;
  Double_t path_dev = path*path*nh/12;
  if(nh < 3) return 0.01;
  Double_t d_slope = 1./(nh-2)*res2/path_dev;
  if(std::isnan(d_slope) || std::isinf(d_slope)){
    std::cout<<Form("Nan || inf dZ resol! nh = %d, dt = %g, res = %g", nh,dt,res2)<<std::endl;
  }
  return m_dZResScale*sqrt(d_slope);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetThetaResolution(){
  Double_t d_slope = GetdZResolution();
  return d_slope / (1+m_dz*m_dz);
}

//_____________________________________________________________________________
Double_t
TPCLocalTrackHelix::GetMomentumResolution(){
  Double_t d_slope = GetdZResolution();
  return GetTransverseMomentumResolution()*hypot(1,d_slope);
}

//_____________________________________________________________________________
TMatrixD
TPCLocalTrackHelix::GetCovarianceMatrix(){
  double Elements[3*3]={0};
  double cov_mom_th = GetMomentumPitchAngleCovariance();
  double cov_mom_ph = GetTransverseMomentumAngularCovariance();
  double res_mom = GetMomentumResolution();
  double res_th = GetThetaResolution();
  double res_ph = GetTransverseAngularResolution();
  Elements[0+3*0] = res_mom*res_mom;
  Elements[1+3*1] = res_th*res_th;
  Elements[2+3*2] = res_ph*res_ph;
  Elements[0*3+1] = cov_mom_th;
  Elements[1*3+0] = cov_mom_th;
  Elements[0*3+2] = cov_mom_ph;
  Elements[2*3+0] = cov_mom_ph;
  TMatrixD CovMat(3,3,Elements);
  return CovMat;
}
