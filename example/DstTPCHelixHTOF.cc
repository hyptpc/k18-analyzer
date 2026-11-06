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
#include <optional>
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
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCPadHelper.hh"
#include "ThreeVector.hh"
#include "UserParamMan.hh"

#include <UnpackerManager.hh>
#include <spdlog/spdlog.h>

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

// Per-track output branches: X(type, name, default value).
#define HTOF_TRACK_BRANCHES(X)            \
  X(Int_t,    charge,               0)    \
  X(Int_t,    is_beam,              0)    \
  X(Int_t,    is_accidental,        0)    \
  X(Int_t,    pid,                  0)    \
  X(Double_t, dEdx,                 qnan) \
  X(Double_t, mom0,                 qnan) \
  X(Int_t,    vertex_source,        kVtxNone) \
  X(Double_t, vtx_x,                qnan) \
  X(Double_t, vtx_y,                qnan) \
  X(Double_t, vtx_z,                qnan) \
  X(Double_t, p_vtx,                qnan) \
  X(Int_t,    seg_match,            0)    \
  X(Double_t, extrap_seg,           qnan) \
  X(Int_t,    extrap_plane,         -1)   \
  X(Double_t, extrap_x,             qnan) \
  X(Double_t, extrap_y,             qnan) \
  X(Double_t, extrap_z,             qnan) \
  X(Double_t, extrap_horizontal,    qnan) \
  X(Double_t, extrap_vertical,      qnan) \
  X(Double_t, extrap_cos_incidence, qnan) \
  X(Double_t, match_cl_seg,         qnan) \
  X(Double_t, match_raw_seg,        qnan) \
  X(Double_t, match_hit_seg,        qnan) \
  X(Double_t, htof_adc_u,           qnan) \
  X(Double_t, htof_adc_d,           qnan) \
  X(Double_t, htof_adc_s,           qnan) \
  X(Double_t, htof_de_u,            qnan) \
  X(Double_t, htof_de_d,            qnan) \
  X(Double_t, htof_de_s,            qnan) \
  X(Int_t,    tof_valid,            0)    \
  X(Double_t, ctof_htof,            qnan) \
  X(Double_t, L_beam,               qnan) \
  X(Double_t, L_sec,                qnan) \
  X(Double_t, t_beam,               qnan) \
  X(Double_t, t_sec,                qnan) \
  X(Double_t, m2,                   qnan) \
  X(Double_t, dt_pi,                qnan) \
  X(Double_t, dt_k,                 qnan) \
  X(Double_t, dt_p,                 qnan) \
  X(Double_t, nsigma_pi,            qnan) \
  X(Double_t, nsigma_k,             qnan) \
  X(Double_t, nsigma_p,             qnan)

// One track of the output.
struct TrackRow
{
#define X(type, name, init) type name = init;
  HTOF_TRACK_BRANCHES(X)
#undef X
};
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
  Int_t status;
  UInt_t runnum;
  UInt_t evnum;
  std::vector<Double_t> trigpat;
  std::vector<std::vector<Double_t>> trigflag;
  Int_t beamflag;
  Double_t time0;
  Double_t time0_seg;
  Double_t time0_seg_corr;
  Int_t ntTpc;
#define X(type, name, init) std::vector<type> name;
  HTOF_TRACK_BRANCHES(X)
#undef X

  void clear()
  {
    status = 0;
    runnum = 0;
    evnum = 0;
    trigpat.clear();
    trigflag.clear();
    beamflag = beam::kUnknown;
    time0 = qnan;
    time0_seg = qnan;
    time0_seg_corr = 0.;
    ntTpc = 0;
#define X(type, name, init) name.clear();
    HTOF_TRACK_BRANCHES(X)
#undef X
  }

  void push(const TrackRow& row)
  {
    ++ntTpc;
#define X(type, name, init) name.push_back(row.name);
    HTOF_TRACK_BRANCHES(X)
#undef X
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
IsBeamOrAccidental(Int_t i)
{
  return (**src.is_beam)[i] == 1 || (**src.is_accidental)[i] == 1;
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
  Double_t par[5] = { (**src.helix_cx)[it], (**src.helix_cy)[it], (**src.helix_z0)[it],
                      (**src.helix_r)[it], (**src.helix_dz)[it] };
  helix.SetParam(par);
  helix.SetCharge((**src.charge)[it]);
  helix.SetMint(tmin);
  helix.SetMaxt(tmax);
  helix.SetIsThetaCalculated(true);
  return true;
}

//_____________________________________________________________________________
// Pair vertex of secondary track it: the non-beam, non-accidental partner with the smallest
// closest distance in the helix DST.
Bool_t
FindPairVertex(Int_t it, Int_t nt, TVector3& vertex)
{
  if (it >= static_cast<Int_t>((**src.closeDistTpc).size())) return false;
  const auto& dist = (**src.closeDistTpc)[it];
  Double_t best = std::numeric_limits<Double_t>::max();
  for (Int_t j = 0; j < nt && j < static_cast<Int_t>(dist.size()); ++j) {
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
FindVertex(Int_t it, Int_t nt, const TPCLocalTrackHelix& helix, TVector3& vertex)
{
  if (IsBeamOrAccidental(it)) {
    TVector3 mom;
    Double_t len = 0., dist = 0.;
    if (helix.ExtrapolateToTarget(vertex, mom, len, dist) && tpc::IsInsideTargetCylinder(vertex))
      return kVtxBeamHelix;
    return kVtxNone;
  }
  if (FindPairVertex(it, nt, vertex) && tpc::IsInsideTargetCylinder(vertex))
    return kVtxPair;
  return kVtxNone;
}

//_____________________________________________________________________________
void
FillCandHist(const TPCAnalyzer& tpc_ana, const HtofCands& c)
{
  HF1("HTOF_Cand_N", c.seg.size());
  for (std::size_t i = 0; i < c.seg.size(); ++i) {
    Double_t abs_s = qnan, drho = qnan;
    tpc_ana.HtofMatchResidual(c.plane[i], c.pos[i], abs_s, drho);
    HF1("HTOF_Cand_Seg", c.seg[i]);
    HF1("HTOF_Cand_Len", c.len[i]);
    HF1("HTOF_Cand_AbsS", abs_s);
    HF1("HTOF_Cand_dRho", drho);
    HF1("HTOF_Cand_H", c.h[i]);
    HF1("HTOF_Cand_V", c.v[i]);
  }
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
FillMatch(TrackRow& row, const TPCAnalyzer& tpc_ana, const HtofCands& c, const HtofMatch& m)
{
  const Int_t i = m.cand;
  row.seg_match = 1;
  row.extrap_seg = c.seg[i];
  row.extrap_plane = c.plane[i];
  row.extrap_x = c.pos[i].X();
  row.extrap_y = c.pos[i].Y();
  row.extrap_z = c.pos[i].Z();
  row.extrap_horizontal = c.h[i];
  row.extrap_vertical = c.v[i];
  if (c.mom[i].Mag() > 0.)
    row.extrap_cos_incidence = TMath::Abs(c.mom[i].Unit().Dot(tpc_ana.HtofPlaneNormal(c.plane[i])));
  row.match_cl_seg = (**src.htof_cl_seg)[m.cl];
}

//_____________________________________________________________________________
// ADC (raw-indexed) and dE (hit-indexed) of the segment closest to the matched cluster;
// match_raw_seg / match_hit_seg tell which segment they belong to.
void
FillAdcDe(TrackRow& row)
{
  const Int_t ir = TPCAnalyzer::MatchHtofBySeg(row.match_cl_seg, **src.htof_raw_seg);
  if (ir >= 0 && ir < static_cast<Int_t>((**src.htof_adc_s).size())) {
    row.match_raw_seg = (**src.htof_raw_seg)[ir];
    row.htof_adc_u = (**src.htof_adc_u)[ir];
    row.htof_adc_d = (**src.htof_adc_d)[ir];
    row.htof_adc_s = (**src.htof_adc_s)[ir];
  }
  const Int_t ih = TPCAnalyzer::MatchHtofBySeg(row.match_cl_seg, **src.htof_hit_seg);
  if (ih >= 0 && ih < static_cast<Int_t>((**src.htof_de_s).size())) {
    row.match_hit_seg = (**src.htof_hit_seg)[ih];
    row.htof_de_u = (**src.htof_de_u)[ih];
    row.htof_de_d = (**src.htof_de_d)[ih];
    row.htof_de_s = (**src.htof_de_s)[ih];
  }
}

//_____________________________________________________________________________
void
FillMatchHist(const TPCAnalyzer& tpc_ana, const HtofCands& c, const HtofMatch& m, Double_t L_sec)
{
  Double_t abs_s = qnan, drho = qnan;
  tpc_ana.HtofMatchResidual(c.plane[m.cand], c.pos[m.cand], abs_s, drho);
  HF1("HTOF_Match_AbsS", abs_s);
  HF1("HTOF_Match_dRho", drho);
  HF1("HTOF_Match_H", c.h[m.cand]);
  HF1("HTOF_Match_V", c.v[m.cand]);
  HF1("HTOF_Match_dSeg", m.dseg);
  if (TMath::Finite(L_sec)) HF1("HTOF_Match_Lsec", L_sec);
}

//_____________________________________________________________________________
// TOF of a matched track with a vertex. Secondary flight time = (HTOF - time0 - BH2 correction)
// minus the beam flight BH2 -> vertex. dt_X are "if it were X" residuals, not an assignment.
Int_t
FillTof(TrackRow& row, const TVector3& vertex, const BeamInfo& beam, Double_t cl_time)
{
  if (!TMath::Finite(event.time0) || !TMath::Finite(cl_time)) return kStageBadTime;
  row.ctof_htof = cl_time - event.time0 - event.time0_seg_corr;
  row.L_beam = (vertex - Bh2SegPositionTpc(event.time0_seg)).Mag();
  if (beam.Valid())
    row.t_beam = Kinematics::CalcTimeOfFlight(beam.mom, row.L_beam, beam.mass);
  if (!TMath::Finite(row.t_beam)) return kStageNoBeam;
  row.t_sec = row.ctof_htof - row.t_beam;
  const Double_t p = row.p_vtx, L = row.L_sec, t = row.t_sec;
  if (!(t > 0.) || !(p > 0.)) return kStageBadTof;
  const Double_t poq = p * row.charge;
  row.m2 = Kinematics::MassSquare(p, L, t);
  row.dt_pi = t - Kinematics::CalcTimeOfFlight(p, L, pdg::PionMass());
  row.dt_k  = t - Kinematics::CalcTimeOfFlight(p, L, pdg::KaonMass());
  row.dt_p  = t - Kinematics::CalcTimeOfFlight(p, L, pdg::ProtonMass());
  row.nsigma_pi = Kinematics::HypTPCHTOFNsigmaPion(poq, L, t);
  row.nsigma_k  = Kinematics::HypTPCHTOFNsigmaKaon(poq, L, t);
  row.nsigma_p  = Kinematics::HypTPCHTOFNsigmaProton(poq, L, t);
  row.tof_valid = 1;
  return kStageTof;
}

//_____________________________________________________________________________
void
FillTofHist(const TrackRow& row)
{
  HF1("HTOF_Ctof", row.ctof_htof);
  HF1("HTOF_Lsec", row.L_sec);
  HF1("HTOF_Lbeam", row.L_beam);
  HF1("HTOF_M2", row.m2);
  HF1(Form("HTOF_M2_VtxSrc%d", row.vertex_source), row.m2);
  HF1("HTOF_dT_Pi", row.dt_pi);
  HF1("HTOF_dT_K", row.dt_k);
  HF1("HTOF_dT_P", row.dt_p);
  HF2("HTOF_dT_Pi_vs_Seg", row.extrap_seg, row.dt_pi);
  HF2("HTOF_dT_Pi_vs_Lsec", row.L_sec, row.dt_pi);
  const Double_t poq = row.p_vtx * row.charge;
  const Double_t inv_beta = MathTools::C() * row.t_sec / row.L_sec;
  HF2("HTOF_InvBeta_vs_PoQ", poq, inv_beta);
  if (row.pid & 0x1) HF2("HTOF_InvBeta_vs_PoQ_PiPid", poq, inv_beta);
  if (row.pid & 0x2) HF2("HTOF_InvBeta_vs_PoQ_KPid", poq, inv_beta);
  if (row.pid & 0x4) HF2("HTOF_InvBeta_vs_PoQ_PPid", poq, inv_beta);
}

//_____________________________________________________________________________
// Track part after the vertex: HTOF candidates, cluster match, ADC/dE and TOF.
// Returns the stage reached (EStage).
Int_t
ProcessHtof(TrackRow& row, const TPCLocalTrackHelix& helix, const TPCAnalyzer& tpc_ana,
            const TVector3& vertex, Double_t len_vtx, const BeamInfo& beam)
{
  HtofCands c;
  if (!tpc_ana.ExtrapolateToHTOF(&helix, c.seg, c.pos, c.mom, c.len, c.plane, c.h, c.v))
    c = HtofCands{};
  FillCandHist(tpc_ana, c);
  const Bool_t have_vtx = (row.vertex_source != kVtxNone);
  // With a vertex: only candidates beyond it; without: only ahead of the track start.
  const HtofMatch m = MatchHtof(c, have_vtx ? len_vtx : 0.);
  if (!m.Found())
    return c.seg.empty() ? kStageNoCandidate : kStageNoCluster;
  FillMatch(row, tpc_ana, c, m);
  FillAdcDe(row);
  if (have_vtx) row.L_sec = c.len[m.cand] - len_vtx;
  FillMatchHist(tpc_ana, c, m, row.L_sec);
  if (!have_vtx) return kStageNoVertex;
  return FillTof(row, vertex, beam, (**src.htof_cl_time)[m.cl]);
}

//_____________________________________________________________________________
TrackRow
ProcessTrack(Int_t it, Int_t nt, const TPCAnalyzer& tpc_ana, const BeamInfo& beam)
{
  TrackRow row;
  row.charge        = (**src.charge)[it];
  row.is_beam       = (**src.is_beam)[it];
  row.is_accidental = (**src.is_accidental)[it];
  row.pid           = (**src.pid)[it];
  row.dEdx          = (**src.dEdx)[it];
  row.mom0          = (**src.mom0)[it];

  TPCLocalTrackHelix helix;
  if (!MakeHelix(it, helix)) return row;

  TVector3 vertex;
  row.vertex_source = FindVertex(it, nt, helix, vertex);
  Double_t len_vtx = 0.;
  if (row.vertex_source != kVtxNone) {
    row.vtx_x = vertex.X();
    row.vtx_y = vertex.Y();
    row.vtx_z = vertex.Z();
    TVector3 pos, mom;
    Double_t dist = 0.;
    if (!helix.ExtrapolateToPoint(vertex, pos, mom, len_vtx, dist)) return row;
    // |p| is constant along the helix: use the helix DST value (mom0).
    row.p_vtx = row.mom0;
  }

  const Int_t stage = ProcessHtof(row, helix, tpc_ana, vertex, len_vtx, beam);
  HF1("HTOF_Stage", stage);
  if (row.tof_valid) FillTofHist(row);
  return row;
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
  HF1("Status", event.status++);

  const TPCAnalyzer tpc_ana;
  const BeamInfo beam = GetBeamInfo(event.runnum);
  for (Int_t it = 0; it < nt; ++it)
    event.push(ProcessTrack(it, nt, tpc_ana, beam));
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
  auto* rh = TTreeReaderCont[kHelix];
  dst::SetBranch(rh, "run_number", src.runnum_tpc);
  dst::SetBranch(rh, "event_number", src.evnum_tpc);
  dst::SetBranch(rh, "trig_pat", src.trigpat);
  dst::SetBranch(rh, "trig_flag", src.trigflag);
  dst::SetBranch(rh, "ntTpc", src.ntTpc);
  dst::SetBranch(rh, "helix_cx", src.helix_cx);
  dst::SetBranch(rh, "helix_cy", src.helix_cy);
  dst::SetBranch(rh, "helix_z0", src.helix_z0);
  dst::SetBranch(rh, "helix_r", src.helix_r);
  dst::SetBranch(rh, "helix_dz", src.helix_dz);
  dst::SetBranch(rh, "helix_theta_min", src.helix_theta_min);
  dst::SetBranch(rh, "helix_theta_max", src.helix_theta_max);
  dst::SetBranch(rh, "charge", src.charge);
  dst::SetBranch(rh, "is_beam", src.is_beam);
  dst::SetBranch(rh, "is_accidental", src.is_accidental);
  dst::SetBranch(rh, "pid", src.pid);
  dst::SetBranch(rh, "dEdx", src.dEdx);
  dst::SetBranch(rh, "mom0", src.mom0);
  dst::SetBranch(rh, "closeDistTpc", src.closeDistTpc);
  dst::SetBranch(rh, "vtxTpc", src.vtxTpc);
  dst::SetBranch(rh, "vtyTpc", src.vtyTpc);
  dst::SetBranch(rh, "vtzTpc", src.vtzTpc);

  if (!dst::SetupReader(kHodo, "kHodo")) return false;
  auto* rd = TTreeReaderCont[kHodo];
  dst::SetBranch(rd, "run_number", src.runnum_hodo);
  dst::SetBranch(rd, "event_number", src.evnum_hodo);
  dst::SetBranch(rd, "beam_flag", src.beamflag);
  dst::SetBranch(rd, "time0", src.time0);
  dst::SetBranch(rd, "time0_seg", src.time0_seg);
  dst::SetBranch(rd, "htof_cl_seg", src.htof_cl_seg);
  dst::SetBranch(rd, "htof_cl_time", src.htof_cl_time);
  dst::SetBranch(rd, "htof_hit_seg", src.htof_hit_seg);
  dst::SetBranch(rd, "htof_raw_seg", src.htof_raw_seg);
  dst::SetBranch(rd, "htof_adc_u", src.htof_adc_u);
  dst::SetBranch(rd, "htof_adc_d", src.htof_adc_d);
  dst::SetBranch(rd, "htof_adc_s", src.htof_adc_s);
  dst::SetBranch(rd, "htof_de_u", src.htof_de_u);
  dst::SetBranch(rd, "htof_de_d", src.htof_de_d);
  dst::SetBranch(rd, "htof_de_s", src.htof_de_s);

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
#define X(type, name, init) tree->Branch(#name, &event.name);
  HTOF_TRACK_BRANCHES(X)
#undef X

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
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
