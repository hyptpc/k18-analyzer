// -*- C++ -*-
//
// HTOF time of flight for TPC helix tracks.
// Inputs (entry-aligned): DstTPCHelixTracking ("tpc") and UserHodoscope ("hodo").
// Per track: vertex -> helix extrapolated to HTOF -> HTOF cluster match -> TOF and m2.
// Helix parameters are taken from the input tree (no re-tracking).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <TMath.h>
#include <TVector3.h>

#include "BeamRunCondition.hh"
#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "HodoParamMan.hh"
#include "Kinematics.hh"
#include "Lh2ElossMan.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCPadHelper.hh"
#include "ThreeVector.hh"
#include "UserParamMan.hh"

#include <UnpackerManager.hh>
#include <spdlog/spdlog.h>

// 1: also fill the LH2 TOF branches for beam tracks (kVtxBeamHelix) into separate *_beam branches.
#define HTOF_LH2_INCLUDE_BEAM 0

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
auto&       gLh2Eloss = Lh2ElossMan::GetInstance();

const Double_t qnan = TMath::QuietNaN();
const std::vector<TString> kUserParamKeys = {};

// Where the vertex of a track comes from.
enum EVtxSrc : Int_t { kVtxNone = 0, kVtxPair = 1, kVtxBeamHelix = 2 };

// How far a track got (HTOF_Stage histogram); one entry per track that reaches the HTOF step.
enum EStage : Int_t {
  kStageNoCandidate = 0, // helix does not reach an HTOF face beyond the vertex / track start
  kStageNoCluster   = 1, // no HTOF cluster in time on the same or an adjacent segment
  kStageNoVertex    = 2, // matched, but no vertex: match and ADC only
  kStageBadTime     = 3, // time0 or cluster time not finite
  kStageNoBeam      = 4, // beam particle / momentum of the run unresolved, or t_beam not finite
  kStageBadTof      = 5, // t_sec <= 0 or p <= 0
  kStageTof         = 6, // TOF and m2 filled
  kNStage
};

// Cluster match window relative to time0 [ns].
const Double_t kHtofMatchTmin = -3.;
const Double_t kHtofMatchTmax = 30.;

// BH2 time0_seg correction: HDPRM rows with AorT = 1, UorD = 10 (time offset per segment).
constexpr Int_t kTime0SegCorrUorD = 10;
constexpr Int_t kBh2PlId = 0;
Double_t gTime0SegCorr[NumOfSegBH2] = {};
Bool_t   gTime0SegCorrEnabled = false;

// BH2 segment geometry (as in src/DCTrackSearch.cc): 0-origin, -X to +X, centre = seg 7.
const Double_t kBh2SegPitch = 14.0 + 0.5; // width + gap [mm]
}

namespace dst
{
enum kArgc {
  kProcess, kConfFile,
  kHelix, kHodo, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[HelixTracking]", "[Hodoscope]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "hodo", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
std::vector<UInt_t> evnumPerFile;
Bool_t SetupReaders();
}

//_____________________________________________________________________________
struct Event
{
  // --- run / event ---
  Int_t status;
  UInt_t runnum;
  UInt_t evnum;
  Int_t beamflag;

  // --- trigger (from helix) ---
  std::vector<Double_t> trigpat;
  std::vector<std::vector<Double_t>> trigflag;

  // --- BH2 time0 (from hodo) ---
  Double_t time0;
  Double_t time0_seg;
  Double_t time0_seg_corr; // subtracted from ctof_htof

  // --- per-track: copied from helix ---
  Int_t ntTpc;
  std::vector<Int_t> charge;
  std::vector<Int_t> is_beam;
  std::vector<Int_t> is_accidental;
  std::vector<Int_t> pid; // dE/dx PID bits: bit0 = pi, bit1 = K, bit2 = p
  std::vector<Double_t> dEdx;
  std::vector<Double_t> mom0;

  // --- per-track: vertex ---
  std::vector<Int_t> vertex_source; // EVtxSrc
  std::vector<Double_t> vtx_x;
  std::vector<Double_t> vtx_y;
  std::vector<Double_t> vtx_z;
  std::vector<Double_t> p_vtx; // = mom0 when a vertex exists

  // --- per-track: HTOF match ---
  std::vector<Int_t> seg_match;
  std::vector<Double_t> extrap_seg;
  std::vector<Int_t> extrap_plane;
  std::vector<Double_t> extrap_x;
  std::vector<Double_t> extrap_y;
  std::vector<Double_t> extrap_z;
  std::vector<Double_t> extrap_horizontal;    // on-face tangential [mm], signed
  std::vector<Double_t> extrap_vertical;      // on-face Y - y_offset [mm]
  std::vector<Double_t> extrap_cos_incidence; // |cos| to the face normal
  std::vector<Double_t> match_cl_seg;
  std::vector<Double_t> match_raw_seg; // segment of htof_adc_*
  std::vector<Double_t> match_hit_seg; // segment of htof_de_*
  std::vector<Double_t> htof_adc_u;
  std::vector<Double_t> htof_adc_d;
  std::vector<Double_t> htof_adc_s;
  std::vector<Double_t> htof_de_u;
  std::vector<Double_t> htof_de_d;
  std::vector<Double_t> htof_de_s;

  // --- per-track: TOF ---
  std::vector<Int_t> tof_valid;
  std::vector<Double_t> ctof_htof; // HTOF - time0 - time0_seg_corr [ns]
  std::vector<Double_t> L_beam;    // BH2 -> vertex [mm]
  std::vector<Double_t> L_sec;     // vertex -> HTOF along the helix [mm]
  std::vector<Double_t> t_beam;    // BH2 -> vertex flight of the beam [ns]
  std::vector<Double_t> t_sec;     // ctof_htof - t_beam [ns]
  std::vector<Double_t> m2;
  std::vector<Double_t> dt_pi;     // t_sec - expected time if pi (not a species assignment)
  std::vector<Double_t> dt_k;
  std::vector<Double_t> dt_p;
  std::vector<Double_t> nsigma_pi;
  std::vector<Double_t> nsigma_k;
  std::vector<Double_t> nsigma_p;

  // --- per-track: TOF with the LH2 segment (conf LH2ELOSS, pair vertex only) ---
  std::vector<Double_t> L_lh2;
  std::vector<Int_t> lh2_exit;
  std::vector<Double_t> dt_pi_lh2;
  std::vector<Double_t> dt_k_lh2;
  std::vector<Double_t> dt_p_lh2;
#if HTOF_LH2_INCLUDE_BEAM
  std::vector<Double_t> L_lh2_beam;
  std::vector<Int_t> lh2_exit_beam;
  std::vector<Double_t> dt_pi_lh2_beam;
  std::vector<Double_t> dt_k_lh2_beam;
  std::vector<Double_t> dt_p_lh2_beam;
#endif

  void clear()
  {
    status = 0;
    runnum = 0;
    evnum = 0;
    beamflag = beam::kUnknown;
    time0 = qnan;
    time0_seg = qnan;
    time0_seg_corr = 0.;
    ntTpc = 0;
    dst::clear_all(trigpat, trigflag,
                   charge, is_beam, is_accidental, pid, dEdx, mom0,
                   vertex_source, vtx_x, vtx_y, vtx_z, p_vtx,
                   seg_match, extrap_seg, extrap_plane, extrap_x, extrap_y, extrap_z,
                   extrap_horizontal, extrap_vertical, extrap_cos_incidence,
                   match_cl_seg, match_raw_seg, match_hit_seg,
                   htof_adc_u, htof_adc_d, htof_adc_s, htof_de_u, htof_de_d, htof_de_s,
                   tof_valid, ctof_htof, L_beam, L_sec, t_beam, t_sec, m2,
                   dt_pi, dt_k, dt_p, nsigma_pi, nsigma_k, nsigma_p,
                   L_lh2, lh2_exit, dt_pi_lh2, dt_k_lh2, dt_p_lh2);
#if HTOF_LH2_INCLUDE_BEAM
    dst::clear_all(L_lh2_beam, lh2_exit_beam, dt_pi_lh2_beam, dt_k_lh2_beam, dt_p_lh2_beam);
#endif
  }

  void resizeTracks(Int_t n)
  {
    ntTpc = n;
    dst::resize_all(n,
                    charge, is_beam, is_accidental, pid, dEdx, mom0,
                    vertex_source, vtx_x, vtx_y, vtx_z, p_vtx,
                    seg_match, extrap_seg, extrap_plane, extrap_x, extrap_y, extrap_z,
                    extrap_horizontal, extrap_vertical, extrap_cos_incidence,
                    match_cl_seg, match_raw_seg, match_hit_seg,
                    htof_adc_u, htof_adc_d, htof_adc_s, htof_de_u, htof_de_d, htof_de_s,
                    tof_valid, ctof_htof, L_beam, L_sec, t_beam, t_sec, m2,
                    dt_pi, dt_k, dt_p, nsigma_pi, nsigma_k, nsigma_p,
                    L_lh2, lh2_exit, dt_pi_lh2, dt_k_lh2, dt_p_lh2);
#if HTOF_LH2_INCLUDE_BEAM
    dst::resize_all(n, L_lh2_beam, lh2_exit_beam, dt_pi_lh2_beam, dt_k_lh2_beam, dt_p_lh2_beam);
#endif
    for (Int_t i = 0; i < n; ++i) {
      charge[i] = is_beam[i] = is_accidental[i] = pid[i] = 0;
      dEdx[i] = mom0[i] = qnan;
      vertex_source[i] = kVtxNone;
      vtx_x[i] = vtx_y[i] = vtx_z[i] = p_vtx[i] = qnan;
      seg_match[i] = 0;
      extrap_plane[i] = -1;
      extrap_seg[i] = extrap_x[i] = extrap_y[i] = extrap_z[i] = qnan;
      extrap_horizontal[i] = extrap_vertical[i] = extrap_cos_incidence[i] = qnan;
      match_cl_seg[i] = match_raw_seg[i] = match_hit_seg[i] = qnan;
      htof_adc_u[i] = htof_adc_d[i] = htof_adc_s[i] = qnan;
      htof_de_u[i] = htof_de_d[i] = htof_de_s[i] = qnan;
      tof_valid[i] = 0;
      ctof_htof[i] = L_beam[i] = L_sec[i] = t_beam[i] = t_sec[i] = m2[i] = qnan;
      dt_pi[i] = dt_k[i] = dt_p[i] = qnan;
      nsigma_pi[i] = nsigma_k[i] = nsigma_p[i] = qnan;
      L_lh2[i] = dt_pi_lh2[i] = dt_k_lh2[i] = dt_p_lh2[i] = qnan;
      lh2_exit[i] = -1;
#if HTOF_LH2_INCLUDE_BEAM
      L_lh2_beam[i] = dt_pi_lh2_beam[i] = dt_k_lh2_beam[i] = dt_p_lh2_beam[i] = qnan;
      lh2_exit_beam[i] = -1;
#endif
    }
  }
};

//_____________________________________________________________________________
struct Src
{
  // hodo
  TTreeReaderValue<UInt_t>* runnum_hodo;
  TTreeReaderValue<UInt_t>* evnum_hodo;
  TTreeReaderValue<Int_t>* beamflag;
  TTreeReaderValue<Double_t>* time0;
  TTreeReaderValue<Double_t>* time0_seg;
  TTreeReaderValue<std::vector<Double_t>>* htof_cl_seg;
  TTreeReaderValue<std::vector<Double_t>>* htof_cl_time;
  TTreeReaderValue<std::vector<Double_t>>* htof_hit_seg;
  TTreeReaderValue<std::vector<Double_t>>* htof_raw_seg;
  TTreeReaderValue<std::vector<Double_t>>* htof_adc_u;
  TTreeReaderValue<std::vector<Double_t>>* htof_adc_d;
  TTreeReaderValue<std::vector<Double_t>>* htof_adc_s;
  TTreeReaderValue<std::vector<Double_t>>* htof_de_u;
  TTreeReaderValue<std::vector<Double_t>>* htof_de_d;
  TTreeReaderValue<std::vector<Double_t>>* htof_de_s;

  // helix (tpc tree)
  TTreeReaderValue<UInt_t>* runnum_tpc;
  TTreeReaderValue<UInt_t>* evnum_tpc;
  TTreeReaderValue<std::vector<Double_t>>* trigpat;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* trigflag;
  TTreeReaderValue<Int_t>* ntTpc;
  TTreeReaderValue<std::vector<Double_t>>* helix_cx;
  TTreeReaderValue<std::vector<Double_t>>* helix_cy;
  TTreeReaderValue<std::vector<Double_t>>* helix_z0;
  TTreeReaderValue<std::vector<Double_t>>* helix_r;
  TTreeReaderValue<std::vector<Double_t>>* helix_dz;
  TTreeReaderValue<std::vector<Double_t>>* helix_theta_min;
  TTreeReaderValue<std::vector<Double_t>>* helix_theta_max;
  TTreeReaderValue<std::vector<Int_t>>* charge;
  TTreeReaderValue<std::vector<Int_t>>* is_beam;
  TTreeReaderValue<std::vector<Int_t>>* is_accidental;
  TTreeReaderValue<std::vector<Int_t>>* pid;
  TTreeReaderValue<std::vector<Double_t>>* dEdx;
  TTreeReaderValue<std::vector<Double_t>>* mom0;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* closeDistTpc;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtxTpc;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtyTpc;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* vtzTpc;
};

namespace root
{
Event  event;
Src    src;
TTree* tree;
}

namespace
{
using namespace root;

// Beam particle and momentum of the run (for the BH2 -> vertex flight time).
struct BeamInfo
{
  Double_t mass = qnan; // [GeV/c2]
  Double_t mom  = qnan; // [GeV/c]
  Bool_t Valid() const { return TMath::Finite(mass) && TMath::Finite(mom); }
};

// HTOF candidates of one track from TPCAnalyzer::ExtrapolateToHTOF.
struct HtofCands
{
  std::vector<Int_t> seg;
  std::vector<TVector3> pos;
  std::vector<TVector3> mom;
  std::vector<Double_t> len;
  std::vector<Int_t> plane;
  std::vector<Double_t> h;
  std::vector<Double_t> v;
};

// Chosen (candidate, cluster) pair; dseg = 0 same segment, 1 adjacent.
struct HtofMatch
{
  Int_t cand = -1;
  Int_t cl = -1;
  Int_t dseg = 2;
  Double_t dt = std::numeric_limits<Double_t>::max(); // |cluster time - time0|
  Bool_t Found() const { return cand >= 0; }
};

//_____________________________________________________________________________
Bool_t
LoadTime0SegCorr()
{
  const auto& gHodo = HodoParamMan::GetInstance();
  Int_t n_loaded = 0;
  TString msg;
  for (Int_t seg = 0; seg < NumOfSegBH2; ++seg) {
    Double_t t = 0.;
    gTime0SegCorr[seg] = 0.;
    if (gHodo.GetTime(DetIdBH2, kBh2PlId, seg, kTime0SegCorrUorD, 0., t)) {
      gTime0SegCorr[seg] = -t;
      ++n_loaded;
    }
    msg += Form(" %.4f", gTime0SegCorr[seg]);
  }
  gTime0SegCorrEnabled = (n_loaded > 0);
  if (gTime0SegCorrEnabled)
    spdlog::info("Time0SegCorr enabled ({}/{} segs):{}", n_loaded, NumOfSegBH2, msg.Data());
  return gTime0SegCorrEnabled;
}

//_____________________________________________________________________________
// Correction at a (possibly fractional) time0_seg, linear between segments.
Double_t
Time0SegCorrAt(Double_t time0_seg)
{
  if (!gTime0SegCorrEnabled || !TMath::Finite(time0_seg)) return 0.;
  if (time0_seg < 0. || time0_seg > NumOfSegBH2 - 1) return 0.;
  const Int_t i0 = static_cast<Int_t>(TMath::Floor(time0_seg));
  if (i0 >= NumOfSegBH2 - 1) return gTime0SegCorr[NumOfSegBH2 - 1];
  const Double_t f = time0_seg - i0;
  return (1. - f) * gTime0SegCorr[i0] + f * gTime0SegCorr[i0 + 1];
}

//_____________________________________________________________________________
// BH2 segment position in the TPC frame: BH2 local -> global -> TPC ("HypTPC").
TVector3
Bh2SegPositionTpc(Double_t time0_seg)
{
  const Double_t seg_center = 0.5 * (NumOfSegBH2 - 1);
  const Double_t x_local = TMath::IsNaN(time0_seg) ? 0. : kBh2SegPitch * (time0_seg - seg_center);
  const ThreeVector global = gGeom.Local2GlobalPos("BH2", ThreeVector(x_local, 0., 0.));
  const ThreeVector tpc = gGeom.Global2LocalPos("HypTPC", global);
  return TVector3(tpc.x(), tpc.y(), tpc.z());
}

//_____________________________________________________________________________
BeamInfo
GetBeamInfo(UInt_t runnum)
{
  BeamInfo info;
  const ResolvedBeamCondition beam = ResolveBeamCondition(runnum);
  if (!beam.ValidParticle() || !beam.ValidMomentum()) return info;
  const auto mass = BeamMassGeV(beam.particle);
  if (!mass) return info;
  info.mass = *mass;
  info.mom = *beam.momentum_gev_c;
  return info;
}

//_____________________________________________________________________________
Bool_t
IsBeamOrAccidental(Int_t it)
{
  return event.is_beam[it] == 1 || event.is_accidental[it] == 1;
}

//_____________________________________________________________________________
// HTOF segments sharing an edge. Plane 0: 0 | 1 (above the beam window), 2 (below) |
// 3 (above), 4 (below) | 5. Seg 33 and seg 0 meet at the corner of planes 7 and 0.
Bool_t
HtofSegAdjacent(Int_t a, Int_t b)
{
  if (a == b) return true;
  if (a > b) std::swap(a, b);
  if (a <= 5 && b <= 5) {
    static const Int_t pairs[][2] = { {0, 1}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 5} };
    for (const auto& q : pairs) if (a == q[0] && b == q[1]) return true;
    return false;
  }
  if (a == 0 && b == NumOfSegHTOF - 1) return true;
  return (b - a == 1);
}

//_____________________________________________________________________________
// Track copied from the helix DST into a helix object; false if the theta range is missing.
Bool_t
MakeHelix(Int_t it, TPCLocalTrackHelix& helix)
{
  const Double_t tmin = (**src.helix_theta_min)[it];
  const Double_t tmax = (**src.helix_theta_max)[it];
  if (TMath::IsNaN(tmin) || TMath::IsNaN(tmax)) return false;
  Double_t par[5] = {
    (**src.helix_cx)[it], (**src.helix_cy)[it], (**src.helix_z0)[it],
    (**src.helix_r)[it], (**src.helix_dz)[it]
  };
  helix.SetParam(par);
  helix.SetCharge(event.charge[it]);
  helix.SetMint(tmin);
  helix.SetMaxt(tmax);
  helix.SetIsThetaCalculated(true);
  return true;
}

//_____________________________________________________________________________
// Pair vertex of secondary track it: the non-beam, non-accidental partner with the smallest
// closest distance in the helix DST.
Bool_t
FindPairVertex(Int_t it, TVector3& vertex)
{
  if (it >= static_cast<Int_t>((**src.closeDistTpc).size())) return false;
  const auto& dist = (**src.closeDistTpc)[it];
  Double_t best = std::numeric_limits<Double_t>::max();
  for (Int_t j = 0; j < event.ntTpc && j < static_cast<Int_t>(dist.size()); ++j) {
    if (j == it || IsBeamOrAccidental(j) || !(dist[j] < best)) continue;
    const TVector3 v((**src.vtxTpc)[it][j], (**src.vtyTpc)[it][j], (**src.vtzTpc)[it][j]);
    if (TMath::IsNaN(v.X()) || TMath::IsNaN(v.Y()) || TMath::IsNaN(v.Z())) continue;
    best = dist[j];
    vertex = v;
  }
  return best < std::numeric_limits<Double_t>::max();
}

//_____________________________________________________________________________
// Beam / accidental tracks: helix extrapolated to the target. Secondaries: pair vertex.
// Either way the vertex must be inside the target cylinder.
Int_t
FindVertex(Int_t it, const TPCLocalTrackHelix& helix, TVector3& vertex)
{
  if (IsBeamOrAccidental(it)) {
    TVector3 mom;
    Double_t len = 0., dist = 0.;
    if (helix.ExtrapolateToTarget(vertex, mom, len, dist) && tpc::IsInsideTargetCylinder(vertex))
      return kVtxBeamHelix;
    return kVtxNone;
  }
  if (FindPairVertex(it, vertex) && tpc::IsInsideTargetCylinder(vertex))
    return kVtxPair;
  return kVtxNone;
}

//_____________________________________________________________________________
// Best (candidate, cluster) pair: candidates beyond len_min along the helix; clusters in the
// time window on the same or an adjacent segment; same segment first, then closest to time0.
HtofMatch
MatchHtof(const HtofCands& c, Double_t len_min)
{
  const auto& cl_seg = **src.htof_cl_seg;
  const auto& cl_time = **src.htof_cl_time;
  HtofMatch best;
  for (std::size_t ic = 0; ic < c.seg.size(); ++ic) {
    if (c.len[ic] <= len_min) continue;
    for (std::size_t icl = 0; icl < cl_seg.size() && icl < cl_time.size(); ++icl) {
      const Double_t ctof = cl_time[icl] - event.time0;
      if (TMath::IsNaN(cl_seg[icl]) || !TMath::Finite(ctof)) continue;
      if (ctof < kHtofMatchTmin || ctof > kHtofMatchTmax) continue;
      const Int_t s = TMath::Nint(cl_seg[icl]);
      if (!HtofSegAdjacent(c.seg[ic], s)) continue;
      const Int_t d = (s == c.seg[ic]) ? 0 : 1;
      if (d > best.dseg || (d == best.dseg && TMath::Abs(ctof) >= best.dt)) continue;
      best = { static_cast<Int_t>(ic), static_cast<Int_t>(icl), d, TMath::Abs(ctof) };
    }
  }
  return best;
}

//_____________________________________________________________________________
void
FillCandHist(TPCEventAnalyzer& event_ana, const TPCAnalyzer& tpc_ana, const HtofCands& c)
{
  const std::size_t n = c.seg.size();
  std::vector<Double_t> abs_s(n, qnan), drho(n, qnan);
  for (std::size_t i = 0; i < n; ++i)
    tpc_ana.HtofMatchResidual(c.plane[i], c.pos[i], abs_s[i], drho[i]);
  event_ana.FillHelixHtofCandHist(c.seg, c.len, abs_s, drho, c.h, c.v);
}

//_____________________________________________________________________________
void
FillMatch(Int_t it, const TPCAnalyzer& tpc_ana, const HtofCands& c, const HtofMatch& m)
{
  const Int_t i = m.cand;
  event.seg_match[it] = 1;
  event.extrap_seg[it] = c.seg[i];
  event.extrap_plane[it] = c.plane[i];
  event.extrap_x[it] = c.pos[i].X();
  event.extrap_y[it] = c.pos[i].Y();
  event.extrap_z[it] = c.pos[i].Z();
  event.extrap_horizontal[it] = c.h[i];
  event.extrap_vertical[it] = c.v[i];
  if (c.mom[i].Mag() > 0.)
    event.extrap_cos_incidence[it]
      = TMath::Abs(c.mom[i].Unit().Dot(tpc_ana.HtofPlaneNormal(c.plane[i])));
  event.match_cl_seg[it] = (**src.htof_cl_seg)[m.cl];
}

//_____________________________________________________________________________
// ADC (raw-indexed) and dE (hit-indexed) of the segment closest to the matched cluster;
// match_raw_seg / match_hit_seg tell which segment they belong to.
void
FillAdcDe(Int_t it)
{
  const Double_t cl_seg = event.match_cl_seg[it];
  const Int_t ir = TPCAnalyzer::MatchHtofBySeg(cl_seg, **src.htof_raw_seg);
  if (ir >= 0 && ir < static_cast<Int_t>((**src.htof_adc_s).size())) {
    event.match_raw_seg[it] = (**src.htof_raw_seg)[ir];
    event.htof_adc_u[it] = (**src.htof_adc_u)[ir];
    event.htof_adc_d[it] = (**src.htof_adc_d)[ir];
    event.htof_adc_s[it] = (**src.htof_adc_s)[ir];
  }
  const Int_t ih = TPCAnalyzer::MatchHtofBySeg(cl_seg, **src.htof_hit_seg);
  if (ih >= 0 && ih < static_cast<Int_t>((**src.htof_de_s).size())) {
    event.match_hit_seg[it] = (**src.htof_hit_seg)[ih];
    event.htof_de_u[it] = (**src.htof_de_u)[ih];
    event.htof_de_d[it] = (**src.htof_de_d)[ih];
    event.htof_de_s[it] = (**src.htof_de_s)[ih];
  }
}

//_____________________________________________________________________________
void
FillMatchHist(TPCEventAnalyzer& event_ana, const TPCAnalyzer& tpc_ana,
              const HtofCands& c, const HtofMatch& m, Double_t L_sec)
{
  Double_t abs_s = qnan, drho = qnan;
  tpc_ana.HtofMatchResidual(c.plane[m.cand], c.pos[m.cand], abs_s, drho);
  event_ana.FillHelixHtofMatchHist(abs_s, drho, c.h[m.cand], c.v[m.cand], m.dseg, L_sec);
}

//_____________________________________________________________________________
// t_sec minus the expected time of particle X when the LH2 part of the path is flown at the
// momentum before the LH2 loss (trapezoid between vertex and LH2 exit) and the rest at the
// TPC momentum p. Returns dt_plain when there is no LH2 path.
Double_t
DtWithLh2(Double_t t_sec, Double_t p, Double_t L_sec, Double_t dt_plain,
          Double_t mass, Int_t particle, const Lh2ElossMan::Track& trk)
{
  if (!gLh2Eloss.HasParticle(particle)) return qnan;
  if (!(trk.L > 0.)) return dt_plain;
  if (!(trk.L < L_sec)) return qnan;

  Lh2ElossMan::Track trk_exit = trk;
  trk_exit.L = 0.;
  const Double_t dp = gLh2Eloss.Loss(trk, p * 1.e3, particle);           // [MeV/c]
  const Double_t dp_exit = gLh2Eloss.Loss(trk_exit, p * 1.e3, particle); // wall only
  if (!TMath::Finite(dp) || !TMath::Finite(dp_exit)) return qnan;

  const Double_t p_vtx = p + dp * 1.e-3;
  const Double_t p_exit = p + dp_exit * 1.e-3;
  if (!(p_vtx > 0.) || !(p_exit > 0.)) return qnan;
  const Double_t inv_beta_vtx = TMath::Sqrt(mass * mass + p_vtx * p_vtx) / p_vtx;
  const Double_t inv_beta_exit = TMath::Sqrt(mass * mass + p_exit * p_exit) / p_exit;
  const Double_t t_lh2 = 0.5 * trk.L * (inv_beta_vtx + inv_beta_exit) / MathTools::C();
  const Double_t t_rest = Kinematics::CalcTimeOfFlight(p, L_sec - trk.L, mass);
  if (!TMath::Finite(t_lh2) || !TMath::Finite(t_rest)) return qnan;
  return t_sec - (t_lh2 + t_rest);
}

//_____________________________________________________________________________
// LH2 TOF branches of track it (pair vertex; also beam tracks with HTOF_LH2_INCLUDE_BEAM).
void
FillLh2Tof(Int_t it, const TVector3& vertex, const TVector3& mom_at_vtx)
{
  if (!gLh2Eloss.IsReady()) return;
  const Bool_t is_pair = (event.vertex_source[it] == kVtxPair);
#if HTOF_LH2_INCLUDE_BEAM
  const Bool_t is_beam = (event.vertex_source[it] == kVtxBeamHelix);
#else
  const Bool_t is_beam = false;
#endif
  if (!is_pair && !is_beam) return;

  const auto trk = gLh2Eloss.Classify(vertex, mom_at_vtx);
  const Double_t p = event.p_vtx[it];
  const Double_t L = event.L_sec[it];
  const Double_t t = event.t_sec[it];
  const Double_t dt_pi = DtWithLh2(t, p, L, event.dt_pi[it], pdg::PionMass(), Lh2ElossMan::kPion, trk);
  const Double_t dt_k = DtWithLh2(t, p, L, event.dt_k[it], pdg::KaonMass(), Lh2ElossMan::kKaon, trk);
  const Double_t dt_p = DtWithLh2(t, p, L, event.dt_p[it], pdg::ProtonMass(), Lh2ElossMan::kProton, trk);
  if (is_pair) {
    event.L_lh2[it] = trk.L;
    event.lh2_exit[it] = trk.exit;
    event.dt_pi_lh2[it] = dt_pi;
    event.dt_k_lh2[it] = dt_k;
    event.dt_p_lh2[it] = dt_p;
  }
#if HTOF_LH2_INCLUDE_BEAM
  else {
    event.L_lh2_beam[it] = trk.L;
    event.lh2_exit_beam[it] = trk.exit;
    event.dt_pi_lh2_beam[it] = dt_pi;
    event.dt_k_lh2_beam[it] = dt_k;
    event.dt_p_lh2_beam[it] = dt_p;
  }
#endif
}

//_____________________________________________________________________________
// TOF of matched track it with a vertex. Secondary flight time = (HTOF - time0 - BH2 correction)
// minus the beam flight BH2 -> vertex. dt_X are "if it were X" residuals, not an assignment.
Int_t
FillTof(Int_t it, const TVector3& vertex, const TVector3& mom_at_vtx,
        const BeamInfo& beam, Double_t cl_time)
{
  if (!TMath::Finite(event.time0) || !TMath::Finite(cl_time)) return kStageBadTime;
  event.ctof_htof[it] = cl_time - event.time0 - event.time0_seg_corr;
  event.L_beam[it] = (vertex - Bh2SegPositionTpc(event.time0_seg)).Mag();
  if (beam.Valid())
    event.t_beam[it] = Kinematics::CalcTimeOfFlight(beam.mom, event.L_beam[it], beam.mass);
  if (!TMath::Finite(event.t_beam[it])) return kStageNoBeam;
  event.t_sec[it] = event.ctof_htof[it] - event.t_beam[it];

  const Double_t p = event.p_vtx[it], L = event.L_sec[it], t = event.t_sec[it];
  if (!(t > 0.) || !(p > 0.)) return kStageBadTof;
  const Double_t poq = p * event.charge[it];
  event.m2[it] = Kinematics::MassSquare(p, L, t);
  event.dt_pi[it] = t - Kinematics::CalcTimeOfFlight(p, L, pdg::PionMass());
  event.dt_k[it]  = t - Kinematics::CalcTimeOfFlight(p, L, pdg::KaonMass());
  event.dt_p[it]  = t - Kinematics::CalcTimeOfFlight(p, L, pdg::ProtonMass());
  event.nsigma_pi[it] = Kinematics::HypTPCHTOFNsigmaPion(poq, L, t);
  event.nsigma_k[it]  = Kinematics::HypTPCHTOFNsigmaKaon(poq, L, t);
  event.nsigma_p[it]  = Kinematics::HypTPCHTOFNsigmaProton(poq, L, t);
  FillLh2Tof(it, vertex, mom_at_vtx);
  event.tof_valid[it] = 1;
  return kStageTof;
}

//_____________________________________________________________________________
void
FillTofHist(TPCEventAnalyzer& event_ana, Int_t it)
{
  event_ana.FillHelixHtofTofHist(event.ctof_htof[it], event.L_sec[it], event.L_beam[it],
                                 event.m2[it], event.vertex_source[it],
                                 event.p_vtx[it], event.charge[it], event.pid[it],
                                 event.t_sec[it],
                                 event.dt_pi[it], event.dt_k[it], event.dt_p[it],
                                 event.extrap_seg[it]);
  if (event.vertex_source[it] == kVtxPair)
    event_ana.FillHelixHtofLh2Hist(event.L_lh2[it], event.dt_pi_lh2[it],
                                   event.dt_p_lh2[it], event.dt_p[it], event.p_vtx[it]);
}

//_____________________________________________________________________________
// Track part after the vertex: HTOF candidates, cluster match, ADC/dE and TOF.
// Returns the stage reached (EStage).
Int_t
ProcessHtof(Int_t it, const TPCLocalTrackHelix& helix,
            TPCEventAnalyzer& event_ana, const TPCAnalyzer& tpc_ana,
            const TVector3& vertex, Double_t len_vtx, const TVector3& mom_at_vtx,
            const BeamInfo& beam)
{
  HtofCands c;
  if (!tpc_ana.ExtrapolateToHTOF(&helix, c.seg, c.pos, c.mom, c.len, c.plane, c.h, c.v))
    c = HtofCands{};
  FillCandHist(event_ana, tpc_ana, c);
  const Bool_t have_vtx = (event.vertex_source[it] != kVtxNone);
  // With a vertex: only candidates beyond it; without: only ahead of the track start.
  const HtofMatch m = MatchHtof(c, have_vtx ? len_vtx : 0.);
  if (!m.Found())
    return c.seg.empty() ? kStageNoCandidate : kStageNoCluster;
  FillMatch(it, tpc_ana, c, m);
  FillAdcDe(it);
  if (have_vtx) event.L_sec[it] = c.len[m.cand] - len_vtx;
  FillMatchHist(event_ana, tpc_ana, c, m, event.L_sec[it]);
  if (!have_vtx) return kStageNoVertex;
  return FillTof(it, vertex, mom_at_vtx, beam, (**src.htof_cl_time)[m.cl]);
}

//_____________________________________________________________________________
void
ProcessTrack(Int_t it, TPCEventAnalyzer& event_ana, const TPCAnalyzer& tpc_ana,
             const BeamInfo& beam)
{
  TPCLocalTrackHelix helix;
  if (!MakeHelix(it, helix)) return;

  TVector3 vertex, mom_at_vtx;
  event.vertex_source[it] = FindVertex(it, helix, vertex);
  Double_t len_vtx = 0.;
  if (event.vertex_source[it] != kVtxNone) {
    event.vtx_x[it] = vertex.X();
    event.vtx_y[it] = vertex.Y();
    event.vtx_z[it] = vertex.Z();
    TVector3 pos;
    Double_t dist = 0.;
    if (!helix.ExtrapolateToPoint(vertex, pos, mom_at_vtx, len_vtx, dist)) return;
    // |p| is constant along the helix: use the helix DST value (mom0).
    event.p_vtx[it] = event.mom0[it];
  }

  const Int_t stage = ProcessHtof(it, helix, event_ana, tpc_ana, vertex, len_vtx, mom_at_vtx, beam);
  event_ana.FillHelixHtofStage(stage);
  if (event.tof_valid[it]) FillTofHist(event_ana, it);
}
}

//_____________________________________________________________________________
Int_t
main(Int_t argc, Char_t** argv)
{
  std::vector<std::string> arg(argv, argv + argc);

  if (!CheckArg(arg))
    return EXIT_FAILURE;
  if (!DstOpen(arg))
    return EXIT_FAILURE;
  if (!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if (!dst::ValidateUserParams(gUser, kUserParamKeys))
    return EXIT_FAILURE;
  if (!gConf.InitializeHistograms())
    return EXIT_FAILURE;
  if (!gConf.InitializeUnpacker())
    return EXIT_FAILURE;
  if (!dst::SetupReaders())
    return EXIT_FAILURE;

  Long64_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  const Long64_t max_loop = gUnpacker.get_max_loop();
  const Long64_t nentries = GetEntries(TTreeCont);
  Long64_t nevent = nentries;
  if (max_loop > 0) nevent = std::min(nentries, skip + max_loop);

  CatchSignal::Set();

  Int_t exit_code = EXIT_SUCCESS;
  Long64_t ievent = skip;
  for (; ievent < nevent && !CatchSignal::Stop(); ++ievent) {
    gCounter.check();
    InitializeEvent();
    const Bool_t ok = DstRead(ievent);
    if (dst::LastGetEntryFailed() || dst::LastSyncFailed()) {
      exit_code = EXIT_FAILURE;
      break;
    }
    if (ok) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  std::_Exit(exit_code);
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t n_input_files = 0;
  for (const auto& name : TreeName) if (name != "") n_input_files++;

  Int_t open_file = 0;
  Int_t open_tree = 0;
  for (Int_t i = 0; i < nArgc; ++i) {
    if (TreeName[i] == "") continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if (open_file != n_input_files || open_tree != n_input_files) {
    spdlog::error("DstOpen Failed: expected {}, open_file {}, open_tree {}",
                  n_input_files, open_file, open_tree);
    return false;
  }
  if (!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");
  dst::ResetRunNumberCheck();
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(Long64_t ievent)
{
  if (ievent % 1000 == 0) {
    std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  }
  if (!GetEntry(ievent))
    return false;

  if (!dst::EnsureRunNumbersMatch(
          { **src.runnum_tpc, **src.runnum_hodo }, ievent,
          { TreeName[kHelix], TreeName[kHodo] }))
    return false;

  evnumPerFile = { **src.evnum_tpc, **src.evnum_hodo };
  if (!dst::CheckEventNumbers(evnumPerFile, ievent,
                              { TreeName[kHelix], TreeName[kHodo] })) {
    dst::LastSyncFailedFlag() = true;
    return false;
  }

  event.runnum    = **src.runnum_tpc;
  event.evnum     = **src.evnum_tpc;
  event.trigpat   = **src.trigpat;
  event.trigflag  = **src.trigflag;
  event.beamflag  = **src.beamflag;
  event.time0     = **src.time0;
  event.time0_seg = **src.time0_seg;
  event.time0_seg_corr = Time0SegCorrAt(event.time0_seg);
  HF1("Status", event.status++);

  const Int_t nt = **src.ntTpc;
  if (nt <= 0) return true;
  event.resizeTracks(nt);
  // All tracks are copied first: the pair vertex of a track looks at the tags of the others.
  for (Int_t it = 0; it < nt; ++it) {
    event.charge[it]        = (**src.charge)[it];
    event.is_beam[it]       = (**src.is_beam)[it];
    event.is_accidental[it] = (**src.is_accidental)[it];
    event.pid[it]           = (**src.pid)[it];
    event.dEdx[it]          = (**src.dEdx)[it];
    event.mom0[it]          = (**src.mom0)[it];
  }
  HF1("Status", event.status++);

  static TPCEventAnalyzer event_ana;
  const TPCAnalyzer tpc_ana;
  const BeamInfo beam = GetBeamInfo(event.runnum);
  for (Int_t it = 0; it < nt; ++it)
    ProcessTrack(it, event_ana, tpc_ana, beam);
  HF1("Status", event.status++);

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  spdlog::info(" Close : {}", TFileCont[kOutFile]->GetName());
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for (Int_t i = 0; i < n; ++i) {
    if (TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if (TTreeCont[i]) delete TTreeCont[i];
    if (TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::SetupReaders()
{
  if (!dst::SetupReader(kHelix, "kHelix")) return false;
  dst::SetBranch(TTreeReaderCont[kHelix], "run_number", src.runnum_tpc);
  dst::SetBranch(TTreeReaderCont[kHelix], "event_number", src.evnum_tpc);
  dst::SetBranch(TTreeReaderCont[kHelix], "trig_pat", src.trigpat);
  dst::SetBranch(TTreeReaderCont[kHelix], "trig_flag", src.trigflag);
  dst::SetBranch(TTreeReaderCont[kHelix], "ntTpc", src.ntTpc);
  dst::SetBranch(TTreeReaderCont[kHelix], "helix_cx", src.helix_cx);
  dst::SetBranch(TTreeReaderCont[kHelix], "helix_cy", src.helix_cy);
  dst::SetBranch(TTreeReaderCont[kHelix], "helix_z0", src.helix_z0);
  dst::SetBranch(TTreeReaderCont[kHelix], "helix_r", src.helix_r);
  dst::SetBranch(TTreeReaderCont[kHelix], "helix_dz", src.helix_dz);
  dst::SetBranch(TTreeReaderCont[kHelix], "helix_theta_min", src.helix_theta_min);
  dst::SetBranch(TTreeReaderCont[kHelix], "helix_theta_max", src.helix_theta_max);
  dst::SetBranch(TTreeReaderCont[kHelix], "charge", src.charge);
  dst::SetBranch(TTreeReaderCont[kHelix], "is_beam", src.is_beam);
  dst::SetBranch(TTreeReaderCont[kHelix], "is_accidental", src.is_accidental);
  dst::SetBranch(TTreeReaderCont[kHelix], "pid", src.pid);
  dst::SetBranch(TTreeReaderCont[kHelix], "dEdx", src.dEdx);
  dst::SetBranch(TTreeReaderCont[kHelix], "mom0", src.mom0);
  dst::SetBranch(TTreeReaderCont[kHelix], "closeDistTpc", src.closeDistTpc);
  dst::SetBranch(TTreeReaderCont[kHelix], "vtxTpc", src.vtxTpc);
  dst::SetBranch(TTreeReaderCont[kHelix], "vtyTpc", src.vtyTpc);
  dst::SetBranch(TTreeReaderCont[kHelix], "vtzTpc", src.vtzTpc);

  if (!dst::SetupReader(kHodo, "kHodo")) return false;
  dst::SetBranch(TTreeReaderCont[kHodo], "run_number", src.runnum_hodo);
  dst::SetBranch(TTreeReaderCont[kHodo], "event_number", src.evnum_hodo);
  dst::SetBranch(TTreeReaderCont[kHodo], "beam_flag", src.beamflag);
  dst::SetBranch(TTreeReaderCont[kHodo], "time0", src.time0);
  dst::SetBranch(TTreeReaderCont[kHodo], "time0_seg", src.time0_seg);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_cl_seg", src.htof_cl_seg);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_cl_time", src.htof_cl_time);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_hit_seg", src.htof_hit_seg);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_raw_seg", src.htof_raw_seg);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_adc_u", src.htof_adc_u);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_adc_d", src.htof_adc_d);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_adc_s", src.htof_adc_s);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_de_u", src.htof_de_u);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_de_d", src.htof_de_d);
  dst::SetBranch(TTreeReaderCont[kHodo], "htof_de_s", src.htof_de_s);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  hist::BuildStatus();
  hist::BuildTPCHelixHTOF();

  tree = new TTree("htof", "tree of DstTPCHelixHTOF");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat", &event.trigpat);
  tree->Branch("trig_flag", &event.trigflag);
  tree->Branch("beam_flag", &event.beamflag);
  tree->Branch("time0", &event.time0);
  tree->Branch("time0_seg", &event.time0_seg);
  tree->Branch("time0_seg_corr", &event.time0_seg_corr);

  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("charge", &event.charge);
  tree->Branch("is_beam", &event.is_beam);
  tree->Branch("is_accidental", &event.is_accidental);
  tree->Branch("pid", &event.pid);
  tree->Branch("dEdx", &event.dEdx);
  tree->Branch("mom0", &event.mom0);

  tree->Branch("vertex_source", &event.vertex_source);
  tree->Branch("vtx_x", &event.vtx_x);
  tree->Branch("vtx_y", &event.vtx_y);
  tree->Branch("vtx_z", &event.vtx_z);
  tree->Branch("p_vtx", &event.p_vtx);

  tree->Branch("seg_match", &event.seg_match);
  tree->Branch("extrap_seg", &event.extrap_seg);
  tree->Branch("extrap_plane", &event.extrap_plane);
  tree->Branch("extrap_x", &event.extrap_x);
  tree->Branch("extrap_y", &event.extrap_y);
  tree->Branch("extrap_z", &event.extrap_z);
  tree->Branch("extrap_horizontal", &event.extrap_horizontal);
  tree->Branch("extrap_vertical", &event.extrap_vertical);
  tree->Branch("extrap_cos_incidence", &event.extrap_cos_incidence);
  tree->Branch("match_cl_seg", &event.match_cl_seg);
  tree->Branch("match_raw_seg", &event.match_raw_seg);
  tree->Branch("match_hit_seg", &event.match_hit_seg);
  tree->Branch("htof_adc_u", &event.htof_adc_u);
  tree->Branch("htof_adc_d", &event.htof_adc_d);
  tree->Branch("htof_adc_s", &event.htof_adc_s);
  tree->Branch("htof_de_u", &event.htof_de_u);
  tree->Branch("htof_de_d", &event.htof_de_d);
  tree->Branch("htof_de_s", &event.htof_de_s);

  tree->Branch("tof_valid", &event.tof_valid);
  tree->Branch("ctof_htof", &event.ctof_htof);
  tree->Branch("L_beam", &event.L_beam);
  tree->Branch("L_sec", &event.L_sec);
  tree->Branch("t_beam", &event.t_beam);
  tree->Branch("t_sec", &event.t_sec);
  tree->Branch("m2", &event.m2);
  tree->Branch("dt_pi", &event.dt_pi);
  tree->Branch("dt_k", &event.dt_k);
  tree->Branch("dt_p", &event.dt_p);
  tree->Branch("nsigma_pi", &event.nsigma_pi);
  tree->Branch("nsigma_k", &event.nsigma_k);
  tree->Branch("nsigma_p", &event.nsigma_p);

  tree->Branch("L_lh2", &event.L_lh2);
  tree->Branch("lh2_exit", &event.lh2_exit);
  tree->Branch("dt_pi_lh2", &event.dt_pi_lh2);
  tree->Branch("dt_k_lh2", &event.dt_k_lh2);
  tree->Branch("dt_p_lh2", &event.dt_p_lh2);
#if HTOF_LH2_INCLUDE_BEAM
  tree->Branch("L_lh2_beam", &event.L_lh2_beam);
  tree->Branch("lh2_exit_beam", &event.lh2_exit_beam);
  tree->Branch("dt_pi_lh2_beam", &event.dt_pi_lh2_beam);
  tree->Branch("dt_k_lh2_beam", &event.dt_k_lh2_beam);
  tree->Branch("dt_p_lh2_beam", &event.dt_p_lh2_beam);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  if (!(InitializeParameter<DCGeomMan>("DCGEO") &&
        InitializeParameter<UserParamMan>("USER")))
    return false;
  if (HasKey("HDPRM")) {
    if (!InitializeParameter<HodoParamMan>("HDPRM"))
      return false;
    LoadTime0SegCorr();
  }
  if (HasKey("LH2ELOSS")) {
    if (!InitializeParameter<Lh2ElossMan>("LH2ELOSS"))
      return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
