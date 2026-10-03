// -*- C++ -*-
//
// Combine DstTPCHelixTracking (tpc) + UserHodoscope (hodo) for HTOF TOF PID.
// Reuses helix parameters from the input tree (no re-tracking).

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
#include "Kinematics.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCPadHelper.hh"
#include "ThreeVector.hh"
#include "UserParamMan.hh"

#include <UnpackerManager.hh>
#include <spdlog/spdlog.h>

// Compile-time switches (same style as DstTPCHelixTracking).
// ExtrapAll: 0 = only tracks with a vertex are extrapolated and matched to HTOF;
//            1 = all tracks (no-vertex tracks get the match and ADC only, no TOF).
#define ExtrapAll 1

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

// vertex_source codes
enum EVtxSrc : Int_t { kVtxNone = 0, kVtxPair = 1, kVtxBeamHelix = 2 };

const std::vector<TString> kUserParamKeys = {};
// BeamMom / BeamParticle: optional via Has(); resolved by
// ResolveBeamCondition(run_num). BeamParticle codes: 0=Pion, 1=Kaon, 2=Proton.

// BH2 segment position in the TPC frame (same frame as the TPC vertex / helix).
// BH2 local (segment x) -> global (DCGEO "BH2") -> TPC frame (DCGEO "HypTPC").
TVector3
Bh2SegPositionTpc(Double_t time0_seg)
{
  // BH2: NumOfSegBH2 (= 15) segs, 0-origin (time0_seg = 0..14), -X→+X; pitch = width + gap;
  // center = (NumOfSegBH2 - 1)/2 = seg 7, as in src/DCTrackSearch.cc.
  const Double_t seg_width = 14.0;
  const Double_t seg_gap = 0.5;
  const Double_t seg_pitch = seg_width + seg_gap; // 14.5
  const Double_t seg_center = 0.5 * (NumOfSegBH2 - 1); // = 7 (0-origin)

  Double_t x_local = 0.;
  if (!TMath::IsNaN(time0_seg))
    x_local = +seg_pitch * (time0_seg - seg_center);

  const ThreeVector bh2_local(x_local, 0., 0.);
  const ThreeVector bh2_global = gGeom.Local2GlobalPos("BH2", bh2_local);
  const ThreeVector bh2_tpc = gGeom.Global2LocalPos("HypTPC", bh2_global);
  return TVector3(bh2_tpc.x(), bh2_tpc.y(), bh2_tpc.z());
}
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

  Int_t ntTpc;
  std::vector<Int_t> charge;
  std::vector<Int_t> is_beam;
  std::vector<Int_t> is_accidental;
  std::vector<Int_t> seg_match;
  std::vector<Int_t> tof_valid;
  std::vector<Int_t> vertex_source;
  // HypTPC dE/dx PID bit pattern (from HelixTracking): bit0=π, bit1=K, bit2=p
  std::vector<Int_t> pid;
  std::vector<Double_t> dEdx;

  std::vector<Double_t> mom0;
  std::vector<Double_t> p_vtx;
  std::vector<Double_t> p_htof;

  std::vector<Double_t> ctof_htof;
  std::vector<Double_t> L_beam;
  std::vector<Double_t> L_sec;
  std::vector<Double_t> t_beam;
  std::vector<Double_t> t_sec;
  // if-X hypotheses (not a species assignment)
  std::vector<Double_t> tof_calc_pi;
  std::vector<Double_t> tof_calc_k;
  std::vector<Double_t> tof_calc_p;
  std::vector<Double_t> dt_pi;
  std::vector<Double_t> dt_k;
  std::vector<Double_t> dt_p;

  std::vector<Double_t> m2;
  std::vector<Double_t> nsigma_p;
  std::vector<Double_t> nsigma_k;
  std::vector<Double_t> nsigma_pi;

  std::vector<Double_t> extrap_seg;
  std::vector<Double_t> match_cl_seg;
  std::vector<Double_t> match_hit_seg;
  std::vector<Double_t> match_raw_seg;
  std::vector<Double_t> htof_adc_u;
  std::vector<Double_t> htof_adc_d;
  std::vector<Double_t> htof_adc_s;
  std::vector<Double_t> htof_de_u;
  std::vector<Double_t> htof_de_d;
  std::vector<Double_t> htof_de_s;
  std::vector<Double_t> extrap_x;
  std::vector<Double_t> extrap_y;
  std::vector<Double_t> extrap_z;
  std::vector<Int_t> extrap_plane;
  std::vector<Double_t> extrap_horizontal; // on-face tangential [mm], signed
  std::vector<Double_t> extrap_vertical;   // on-face Y - y_offset [mm]
  std::vector<Double_t> extrap_px;         // unit momentum direction at the HTOF point (TPC frame)
  std::vector<Double_t> extrap_py;
  std::vector<Double_t> extrap_pz;
  std::vector<Double_t> extrap_cos_incidence; // |cos| of the angle to the face normal (path ~ 1/cos)
  std::vector<Double_t> vtx_x;
  std::vector<Double_t> vtx_y;
  std::vector<Double_t> vtx_z;

  void clear()
  {
    status = 0;
    runnum = 0;
    evnum = 0;
    beamflag = beam::kUnknown;
    time0 = qnan;
    time0_seg = qnan;
    ntTpc = 0;
    dst::clear_all(trigpat, trigflag,
                   charge, is_beam, is_accidental, seg_match, tof_valid, vertex_source,
                   pid, dEdx,
                   mom0, p_vtx, p_htof,
                   ctof_htof, L_beam, L_sec, t_beam, t_sec,
                   tof_calc_pi, tof_calc_k, tof_calc_p,
                   dt_pi, dt_k, dt_p,
                   m2, nsigma_p, nsigma_k, nsigma_pi,
                   extrap_seg, match_cl_seg, match_hit_seg,
                   match_raw_seg,
                   htof_adc_u, htof_adc_d, htof_adc_s,
                   htof_de_u, htof_de_d, htof_de_s,
                   extrap_x, extrap_y, extrap_z, extrap_plane,
                   extrap_horizontal, extrap_vertical,
                   extrap_px, extrap_py, extrap_pz, extrap_cos_incidence,
                   vtx_x, vtx_y, vtx_z);
  }

  void resizeTracks(Int_t n)
  {
    ntTpc = n;
    dst::resize_all(n, charge, is_beam, is_accidental, seg_match, tof_valid, vertex_source,
                    pid, dEdx,
                    mom0, p_vtx, p_htof,
                    ctof_htof, L_beam, L_sec, t_beam, t_sec,
                    tof_calc_pi, tof_calc_k, tof_calc_p,
                    dt_pi, dt_k, dt_p,
                    m2, nsigma_p, nsigma_k, nsigma_pi,
                    extrap_seg, match_cl_seg, match_hit_seg,
                    match_raw_seg,
                    htof_adc_u, htof_adc_d, htof_adc_s,
                    htof_de_u, htof_de_d, htof_de_s,
                    extrap_x, extrap_y, extrap_z, extrap_plane,
                    extrap_horizontal, extrap_vertical,
                    extrap_px, extrap_py, extrap_pz, extrap_cos_incidence,
                    vtx_x, vtx_y, vtx_z);
    for (Int_t i = 0; i < n; ++i) {
      charge[i] = 0;
      is_beam[i] = 0;
      is_accidental[i] = 0;
      seg_match[i] = 0;
      tof_valid[i] = 0;
      vertex_source[i] = kVtxNone;
      pid[i] = 0;
      dEdx[i] = qnan;
      mom0[i] = p_vtx[i] = p_htof[i] = qnan;
      ctof_htof[i] = L_beam[i] = L_sec[i] = t_beam[i] = t_sec[i] = qnan;
      tof_calc_pi[i] = tof_calc_k[i] = tof_calc_p[i] = qnan;
      dt_pi[i] = dt_k[i] = dt_p[i] = qnan;
      m2[i] = nsigma_p[i] = nsigma_k[i] = nsigma_pi[i] = qnan;
      extrap_seg[i] = match_cl_seg[i] = match_hit_seg[i] = qnan;
      match_raw_seg[i] = qnan;
      htof_adc_u[i] = htof_adc_d[i] = htof_adc_s[i] = qnan;
      htof_de_u[i] = htof_de_d[i] = htof_de_s[i] = qnan;
      extrap_x[i] = extrap_y[i] = extrap_z[i] = qnan;
      extrap_plane[i] = -1;
      extrap_horizontal[i] = extrap_vertical[i] = qnan;
      extrap_px[i] = extrap_py[i] = extrap_pz[i] = extrap_cos_incidence[i] = qnan;
      vtx_x[i] = vtx_y[i] = vtx_z[i] = qnan;
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

Bool_t
IsBeamOrAccidental(Int_t i,
                   const std::vector<Int_t>& is_beam,
                   const std::vector<Int_t>& is_accidental)
{
  if (i < 0 || i >= static_cast<Int_t>(is_beam.size())) return false;
  if (is_beam[i] == 1) return true;
  if (i < static_cast<Int_t>(is_accidental.size()) && is_accidental[i] == 1)
    return true;
  return false;
}

// True if the secondary path length, flight time and momentum needed for the TOF PID are all
// finite and positive (a non-positive t_sec means the matched HTOF cluster is too early).
Bool_t
HasValidTofPidInputs(Double_t L_sec, Double_t t_sec, Double_t p_vtx)
{
  return TMath::Finite(L_sec) && L_sec > 0.
      && TMath::Finite(t_sec) && t_sec > 0.
      && TMath::Finite(p_vtx) && p_vtx > 0.;
}

// Vertex of secondary track `it` from the track-pair vertices of the helix DST.
// close_dist[it][j] and (vx, vy, vz)[it][j] are the closest distance and vertex point of the
// pair (it, j) (NaN if not calculated). Among the non-beam, non-accidental partners j, the pair
// with the smallest closest distance is taken. The target-cylinder check is done by the caller.
Bool_t
FindPairVertex(Int_t it, Int_t nt,
               const std::vector<std::vector<Double_t>>& close_dist,
               const std::vector<std::vector<Double_t>>& vx,
               const std::vector<std::vector<Double_t>>& vy,
               const std::vector<std::vector<Double_t>>& vz,
               const std::vector<Int_t>& is_beam,
               const std::vector<Int_t>& is_accidental,
               TVector3& vertex)
{
  Double_t best = std::numeric_limits<Double_t>::max();
  Bool_t found = false;
  if (it < 0 || it >= static_cast<Int_t>(close_dist.size())) return false;
  // Caller must be a secondary; skip beam/accidental partners.
  for (Int_t j = 0; j < nt; ++j) {
    if (j == it) continue;
    if (IsBeamOrAccidental(j, is_beam, is_accidental)) continue;
    if (j >= static_cast<Int_t>(close_dist[it].size())) continue;
    const Double_t d = close_dist[it][j];
    if (TMath::IsNaN(d)) continue;  // no vertex for this pair
    if (d >= best) continue;        // keep the closest pair only
    const Double_t x = vx[it][j], y = vy[it][j], z = vz[it][j];
    if (TMath::IsNaN(x) || TMath::IsNaN(y) || TMath::IsNaN(z)) continue;
    best = d;
    vertex.SetXYZ(x, y, z);
    found = true;
  }
  return found;
}

void
ProcessTrack(Int_t it, TPCAnalyzer& tpc_ana, TPCEventAnalyzer& event_ana)
{
  const Double_t tmin = (**src.helix_theta_min)[it];
  const Double_t tmax = (**src.helix_theta_max)[it];
  if (TMath::IsNaN(tmin) || TMath::IsNaN(tmax)) return;

  TPCLocalTrackHelix helix;
  Double_t par[5] = {
    (**src.helix_cx)[it], (**src.helix_cy)[it], (**src.helix_z0)[it],
    (**src.helix_r)[it], (**src.helix_dz)[it]
  };
  helix.SetParam(par);
  helix.SetCharge((**src.charge)[it]);
  helix.SetMint(tmin);
  helix.SetMaxt(tmax);
  helix.SetIsThetaCalculated(true);

  event.charge[it]        = (**src.charge)[it];
  event.is_beam[it]       = (**src.is_beam)[it];
  event.is_accidental[it] = (**src.is_accidental)[it];
  event.pid[it]           = (**src.pid)[it];
  event.dEdx[it]          = (**src.dEdx)[it];
  event.mom0[it]          = (**src.mom0)[it];

  // ---- vertex ----
  // Beam/accidental: ExtrapolateToTarget only (not used in pair vertexing).
  // Secondary: pair with non-beam/non-accidental partners only; no BeamHelix fallback.
  TVector3 vertex;
  Int_t vtx_src = kVtxNone;
  if (IsBeamOrAccidental(it, event.is_beam, event.is_accidental)) {
    TVector3 pos_t, mom_t;
    Double_t len_t = 0., dist_t = 0.;
    if (helix.ExtrapolateToTarget(pos_t, mom_t, len_t, dist_t)
        && tpc::IsInsideTargetCylinder(pos_t)) {
      vertex = pos_t;
      vtx_src = kVtxBeamHelix;
    }
  } else if (FindPairVertex(it, event.ntTpc,
                            **src.closeDistTpc, **src.vtxTpc, **src.vtyTpc, **src.vtzTpc,
                            event.is_beam, event.is_accidental, vertex)
             && tpc::IsInsideTargetCylinder(vertex)) {
    vtx_src = kVtxPair;
  }
  event.vertex_source[it] = vtx_src;
  const Bool_t have_vtx = (vtx_src != kVtxNone);
#if !ExtrapAll
  if (!have_vtx) return;
#endif

  Double_t len_vtx = 0.;
  if (have_vtx) {
    event.vtx_x[it] = vertex.X();
    event.vtx_y[it] = vertex.Y();
    event.vtx_z[it] = vertex.Z();

    TVector3 pos_at_vtx, mom_at_vtx;
    Double_t dist_vtx = 0.;
    if (!helix.ExtrapolateToPoint(vertex, pos_at_vtx, mom_at_vtx, len_vtx, dist_vtx))
      return;
    event.p_vtx[it] = mom_at_vtx.Mag();
    if (TMath::IsNaN(event.p_vtx[it]) || event.p_vtx[it] <= 0.)
      event.p_vtx[it] = event.mom0[it];
  }

  // ---- HTOF extrapolation ----
  std::vector<Int_t> cand_seg;
  std::vector<TVector3> cand_pos, cand_mom;
  std::vector<Double_t> cand_len;
  std::vector<Int_t> cand_plane;
  std::vector<Double_t> cand_h, cand_v;
  if (!tpc_ana.ExtrapolateToHTOF(&helix, cand_seg, cand_pos, cand_mom, cand_len,
                                 cand_plane, cand_h, cand_v)
      || cand_seg.empty()) {
    event_ana.FillHelixHtofExtrapHist(0, cand_seg, cand_pos, cand_len);
    event_ana.FillHelixHtofPathStage(0); // no HTOF candidate
    return;
  }
  event_ana.FillHelixHtofExtrapHist(static_cast<Int_t>(cand_seg.size()),
                                    cand_seg, cand_pos, cand_len,
                                    cand_plane, cand_h, cand_v);

  // Take the candidate whose seg is closest to any HTOF cluster (over all candidates).
  Int_t best_candi       = -1;
  Int_t best_cl          = -1;
  Double_t best_seg_diff = 1.e9;
  Double_t best_dt       = 1.e9;
  const auto& cl_seg     = **src.htof_cl_seg;
  const auto& cl_time    = **src.htof_cl_time;
  for (std::size_t ic = 0; ic < cand_seg.size(); ++ic) {
    // Vertex tracks: only candidates beyond the vertex along the helix (L_sec > 0); a backward
    // candidate matched to another particle's cluster would otherwise win and lose the TOF.
    if (have_vtx && cand_len[ic] <= len_vtx) continue;
    const Int_t icl = TPCAnalyzer::MatchHtofCluster(static_cast<Double_t>(cand_seg[ic]), cl_seg,
                                                    cl_time, event.time0);
    if (icl < 0) continue;
    const Double_t d = TMath::Abs(static_cast<Double_t>(cand_seg[ic]) - cl_seg[icl]);
    // same |Δseg| across candidates (planes): prefer the cluster closer in time to time0
    const Double_t dt = (TMath::Finite(cl_time[icl]) && TMath::Finite(event.time0))
                        ? TMath::Abs(cl_time[icl] - event.time0) : 1.e9;
    if (d < best_seg_diff - 1e-9 || (TMath::Abs(d - best_seg_diff) <= 1e-9 && dt < best_dt)) {
      best_seg_diff = d;
      best_dt = dt;
      best_candi = static_cast<Int_t>(ic);
      best_cl = icl;
    }
  }
  if (best_candi < 0 || best_cl < 0) {
    event_ana.FillHelixHtofMatchHist(false);
    event_ana.FillHelixHtofPathStage(1); // cand but no cluster match
    return;
  }
  event_ana.FillHelixHtofMatchHist(true);

  event.seg_match[it] = 1;
  event.extrap_seg[it] = cand_seg[best_candi];
  event.match_cl_seg[it] = cl_seg[best_cl];
  event.extrap_x[it] = cand_pos[best_candi].X();
  event.extrap_y[it] = cand_pos[best_candi].Y();
  event.extrap_z[it] = cand_pos[best_candi].Z();
  event.extrap_plane[it] = cand_plane[best_candi];
  event.extrap_horizontal[it] = cand_h[best_candi];
  event.extrap_vertical[it] = cand_v[best_candi];
  event.p_htof[it]   = cand_mom[best_candi].Mag();
  if (cand_mom[best_candi].Mag() > 0.) {
    const TVector3 dir = cand_mom[best_candi].Unit();
    event.extrap_px[it] = dir.X();
    event.extrap_py[it] = dir.Y();
    event.extrap_pz[it] = dir.Z();
    event.extrap_cos_incidence[it] =
      TMath::Abs(dir.Dot(tpc_ana.HtofPlaneNormal(cand_plane[best_candi])));
  }
  if (TMath::IsNaN(event.p_htof[it]) || event.p_htof[it] <= 0.)
    event.p_htof[it] = event.mom0[it];

  // ADC is raw-indexed (htof_raw_seg); dE is hit-indexed (htof_hit_seg).
  const Double_t clseg = cl_seg[best_cl];
  const Int_t iraw = TPCAnalyzer::MatchHtofBySeg(clseg, **src.htof_raw_seg);
  if (iraw >= 0) {
    const auto& raw_seg = **src.htof_raw_seg;
    const auto& adc_u = **src.htof_adc_u;
    const auto& adc_d = **src.htof_adc_d;
    const auto& adc_s = **src.htof_adc_s;
    const std::size_t ir = static_cast<std::size_t>(iraw);
    event.match_raw_seg[it] = raw_seg[ir];
    if (ir < adc_u.size()) event.htof_adc_u[it] = adc_u[ir];
    if (ir < adc_d.size()) event.htof_adc_d[it] = adc_d[ir];
    if (ir < adc_s.size()) event.htof_adc_s[it] = adc_s[ir];
  }
  const Int_t ihit = TPCAnalyzer::MatchHtofBySeg(clseg, **src.htof_hit_seg);
  if (ihit >= 0) {
    const auto& hit_seg = **src.htof_hit_seg;
    const auto& de_u = **src.htof_de_u;
    const auto& de_d = **src.htof_de_d;
    const auto& de_s = **src.htof_de_s;
    const std::size_t ih = static_cast<std::size_t>(ihit);
    event.match_hit_seg[it] = hit_seg[ih];
    if (ih < de_u.size()) event.htof_de_u[it] = de_u[ih];
    if (ih < de_d.size()) event.htof_de_d[it] = de_d[ih];
    if (ih < de_s.size()) event.htof_de_s[it] = de_s[ih];
  }

  // ExtrapAll==1 and no vertex: match/ADC only. No L_sec / TOF PID.
  if (!have_vtx) {
    Double_t abs_s = qnan, drho = qnan;
    if (tpc_ana.HtofMatchResidual(cand_plane[best_candi], cand_pos[best_candi],
                                  abs_s, drho)) {
      event_ana.FillHelixHtofMatchQuality(abs_s, drho,
                                          cand_h[best_candi], cand_v[best_candi],
                                          best_seg_diff, qnan);
    }
    event_ana.FillHelixHtofPathStage(2); // match but no vertex: match/ADC only
    return;
  }

  // ---- TOF components ----
  // ctof_htof = HTOF - Time0 (BTOF-style; BH2→HTOF). Secondary flight is vtx→HTOF,
  // so subtract BH2→vtx (t_beam, beam mass/mom). Measured t_sec has no secondary
  // mass. tof_calc_*/dt_* are if-X hypotheses, not a species assignment.
  const Double_t time0 = event.time0;
  if (!TMath::Finite(time0) || !TMath::Finite(cl_time[best_cl])) {
    event_ana.FillHelixHtofPathStage(4); // match but time0 / cluster time invalid
    return;
  }

  event.ctof_htof[it] = cl_time[best_cl] - time0;
  event.L_sec[it] = cand_len[best_candi] - len_vtx;

  {
    Double_t abs_s = qnan, drho = qnan;
    if (tpc_ana.HtofMatchResidual(cand_plane[best_candi], cand_pos[best_candi],
                                  abs_s, drho)) {
      event_ana.FillHelixHtofMatchQuality(abs_s, drho,
                                          cand_h[best_candi], cand_v[best_candi],
                                          best_seg_diff, event.L_sec[it]);
    }
  }

  if (event.L_sec[it] <= 0.) {
    event_ana.FillHelixHtofPathStage(5); // match but L_sec <= 0
    return;
  }

  const TVector3 bh2_tpc = Bh2SegPositionTpc(event.time0_seg); // vertex is in the TPC frame
  event.L_beam[it] = (vertex - bh2_tpc).Mag();
  const ResolvedBeamCondition beam = ResolveBeamCondition(event.runnum);
  if (!beam.ValidParticle() || !beam.ValidMomentum()) {
    event_ana.FillHelixHtofPathStage(6); // beam condition unresolved
    return;
  }
  const auto m_beam = BeamMassGeV(beam.particle);
  if (!m_beam) {
    event_ana.FillHelixHtofPathStage(6);
    return;
  }
  event.t_beam[it] = Kinematics::CalcTimeOfFlight(
      *beam.momentum_gev_c, event.L_beam[it], *m_beam);
  if (!TMath::Finite(event.t_beam[it])) {
    event_ana.FillHelixHtofPathStage(6); // t_beam not finite
    return;
  }

  // Secondary TOF: BH2→HTOF minus BH2→vtx.
  event.t_sec[it] = event.ctof_htof[it] - event.t_beam[it];

  if (!HasValidTofPidInputs(event.L_sec[it], event.t_sec[it], event.p_vtx[it])) {
    event_ana.FillHelixHtofPathStage(7); // L_sec / t_sec / p_vtx not positive
    return;
  }

  const Double_t p_pid = event.p_vtx[it];
  const Double_t poq = p_pid * static_cast<Double_t>(event.charge[it]);
  event.m2[it] = Kinematics::MassSquare(p_pid, event.L_sec[it], event.t_sec[it]);
  event.nsigma_p[it] = Kinematics::HypTPCHTOFNsigmaProton(poq, event.L_sec[it], event.t_sec[it]);
  event.nsigma_k[it] = Kinematics::HypTPCHTOFNsigmaKaon(poq, event.L_sec[it], event.t_sec[it]);
  event.nsigma_pi[it] = Kinematics::HypTPCHTOFNsigmaPion(poq, event.L_sec[it], event.t_sec[it]);

  event.tof_calc_pi[it] = Kinematics::CalcTimeOfFlight(p_pid, event.L_sec[it], pdg::PionMass());
  event.tof_calc_k[it]  = Kinematics::CalcTimeOfFlight(p_pid, event.L_sec[it], pdg::KaonMass());
  event.tof_calc_p[it]  = Kinematics::CalcTimeOfFlight(p_pid, event.L_sec[it], pdg::ProtonMass());
  event.dt_pi[it] = event.t_sec[it] - event.tof_calc_pi[it];
  event.dt_k[it]  = event.t_sec[it] - event.tof_calc_k[it];
  event.dt_p[it]  = event.t_sec[it] - event.tof_calc_p[it];

  event_ana.FillHelixHtofPidHist(event.ctof_htof[it], event.L_sec[it], event.L_beam[it],
                                 event.m2[it], vtx_src,
                                 event.p_vtx[it], event.charge[it], event.pid[it],
                                 event.t_sec[it],
                                 event.dt_pi[it], event.dt_k[it], event.dt_p[it],
                                 event.extrap_seg[it]);
  event.tof_valid[it] = 1;
  event_ana.FillHelixHtofPathStage(3); // TOF / PID filled
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
    const Bool_t ok = DstRead(static_cast<Int_t>(ievent));
    if (dst::LastGetEntryFailed()) {
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
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if (ievent % 1000 == 0) {
    std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  }
  if (!GetEntry(ievent))
    return false;

  evnumPerFile = { **src.evnum_tpc, **src.evnum_hodo };
  if (!dst::CheckEventNumbers(evnumPerFile, ievent,
                              { TreeName[kHelix], TreeName[kHodo] })) {
    return false;
  }

  event.runnum    = **src.runnum_tpc;
  event.evnum     = **src.evnum_tpc;
  event.trigpat   = **src.trigpat;
  event.trigflag  = **src.trigflag;
  event.beamflag  = **src.beamflag;
  event.time0     = **src.time0;
  event.time0_seg = **src.time0_seg;
  HF1("Status", event.status++);

  const Int_t nt = **src.ntTpc;
  if (nt <= 0) return true;
  event.resizeTracks(nt);
  // Copy beam/accidental tags for all tracks before ProcessTrack so FindPairVertex
  // can exclude partners with j > it (and tracks that early-return on NaN theta).
  for (Int_t it = 0; it < nt; ++it) {
    event.is_beam[it]       = (**src.is_beam)[it];
    event.is_accidental[it] = (**src.is_accidental)[it];
  }
  HF1("Status", event.status++);

  TPCAnalyzer tpc_ana;
  static TPCEventAnalyzer event_ana;
  for (Int_t it = 0; it < nt; ++it)
    ProcessTrack(it, tpc_ana, event_ana);
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
  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("charge", &event.charge);
  tree->Branch("is_beam", &event.is_beam);
  tree->Branch("is_accidental", &event.is_accidental);
  tree->Branch("seg_match", &event.seg_match);
  tree->Branch("tof_valid", &event.tof_valid);
  tree->Branch("vertex_source", &event.vertex_source);
  tree->Branch("pid", &event.pid);
  tree->Branch("dEdx", &event.dEdx);
  tree->Branch("mom0", &event.mom0);
  tree->Branch("p_vtx", &event.p_vtx);
  tree->Branch("p_htof", &event.p_htof);
  tree->Branch("ctof_htof", &event.ctof_htof);
  tree->Branch("L_beam", &event.L_beam);
  tree->Branch("L_sec", &event.L_sec);
  tree->Branch("t_beam", &event.t_beam);
  tree->Branch("t_sec", &event.t_sec);
  tree->Branch("tof_calc_pi", &event.tof_calc_pi);
  tree->Branch("tof_calc_k", &event.tof_calc_k);
  tree->Branch("tof_calc_p", &event.tof_calc_p);
  tree->Branch("dt_pi", &event.dt_pi);
  tree->Branch("dt_k", &event.dt_k);
  tree->Branch("dt_p", &event.dt_p);
  tree->Branch("m2", &event.m2);
  tree->Branch("nsigma_p", &event.nsigma_p);
  tree->Branch("nsigma_k", &event.nsigma_k);
  tree->Branch("nsigma_pi", &event.nsigma_pi);
  tree->Branch("extrap_seg", &event.extrap_seg);
  tree->Branch("match_cl_seg", &event.match_cl_seg);
  tree->Branch("match_hit_seg", &event.match_hit_seg);
  tree->Branch("match_raw_seg", &event.match_raw_seg);
  tree->Branch("htof_adc_u", &event.htof_adc_u);
  tree->Branch("htof_adc_d", &event.htof_adc_d);
  tree->Branch("htof_adc_s", &event.htof_adc_s);
  tree->Branch("htof_de_u", &event.htof_de_u);
  tree->Branch("htof_de_d", &event.htof_de_d);
  tree->Branch("htof_de_s", &event.htof_de_s);
  tree->Branch("extrap_x", &event.extrap_x);
  tree->Branch("extrap_y", &event.extrap_y);
  tree->Branch("extrap_z", &event.extrap_z);
  tree->Branch("extrap_plane", &event.extrap_plane);
  tree->Branch("extrap_horizontal", &event.extrap_horizontal);
  tree->Branch("extrap_vertical", &event.extrap_vertical);
  tree->Branch("extrap_px", &event.extrap_px);
  tree->Branch("extrap_py", &event.extrap_py);
  tree->Branch("extrap_pz", &event.extrap_pz);
  tree->Branch("extrap_cos_incidence", &event.extrap_cos_incidence);
  tree->Branch("vtx_x", &event.vtx_x);
  tree->Branch("vtx_y", &event.vtx_y);
  tree->Branch("vtx_z", &event.vtx_z);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
