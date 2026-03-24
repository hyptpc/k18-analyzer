// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DstHelper.hh"
#include "HistTools.hh"
#include "Kinematics.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCEventAnalyzer.hh"
#include "TPCLocalTrack.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

#include <UnpackerManager.hh>

#define Exclusive 1

#define RawHit 0
#define RawCluster 1
#define TrackSearchFailed 1

//_____________________________________________________________________________
struct Event
{
  Int_t status;
  UInt_t runnum;
  UInt_t evnum;
  std::vector<Double_t> trigpat;
  std::vector<std::vector<Double_t>> trigflag;
  Int_t beamflag;
  std::vector<Double_t> clkTpc;
  std::vector<Double_t> cobo_id;

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
  std::vector<Double_t> chisqrTpc;
  std::vector<Double_t> x0Tpc;
  std::vector<Double_t> y0Tpc;
  std::vector<Double_t> u0Tpc;
  std::vector<Double_t> v0Tpc;
  std::vector<Double_t> theta;
  std::vector<Double_t> dE;
  std::vector<Double_t> dEdx; //reference dedx

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
  std::vector<std::vector<Double_t>> residual_horizontal;
  std::vector<std::vector<Double_t>> residual_vertical;
  std::vector<std::vector<Double_t>> resolution_x;
  std::vector<std::vector<Double_t>> resolution_y;
  std::vector<std::vector<Double_t>> resolution_z;
  std::vector<std::vector<Double_t>> resolution_horizontal;
  std::vector<std::vector<Double_t>> resolution_vertical;
  std::vector<std::vector<Double_t>> pathhit;
  std::vector<std::vector<Double_t>> theta_diff;
  std::vector<std::vector<Double_t>> track_cluster_de;
  std::vector<std::vector<Double_t>> track_cluster_size;
  std::vector<std::vector<Double_t>> track_cluster_mrow;
  std::vector<std::vector<Double_t>> track_cluster_de_center;
  std::vector<std::vector<Double_t>> track_cluster_x_center;
  std::vector<std::vector<Double_t>> track_cluster_y_center;
  std::vector<std::vector<Double_t>> track_cluster_z_center;
  std::vector<std::vector<Double_t>> track_cluster_row_center;

  //exclusive
  std::vector<std::vector<Double_t>> exresidual;
  std::vector<std::vector<Double_t>> exresidual_x;
  std::vector<std::vector<Double_t>> exresidual_y;
  std::vector<std::vector<Double_t>> exresidual_z;
  std::vector<std::vector<Double_t>> exresidual_horizontal;
  std::vector<std::vector<Double_t>> exresidual_vertical;

  Int_t ntTpc_inside;
  Double_t prodvtx_x;
  Double_t prodvtx_y;
  Double_t prodvtx_z;

  //failed track container
  Int_t failed_ntTpc;
  std::vector<Int_t> failed_nhtrack;
  std::vector<Double_t> failed_x0Tpc;
  std::vector<Double_t> failed_y0Tpc;
  std::vector<Double_t> failed_u0Tpc;
  std::vector<Double_t> failed_v0Tpc;

  std::vector<std::vector<Double_t>> failed_hitlayer;
  std::vector<std::vector<Double_t>> failed_hitpos_x;
  std::vector<std::vector<Double_t>> failed_hitpos_y;
  std::vector<std::vector<Double_t>> failed_hitpos_z;
  std::vector<std::vector<Double_t>> failed_calpos_x;
  std::vector<std::vector<Double_t>> failed_calpos_y;
  std::vector<std::vector<Double_t>> failed_calpos_z;

  void clearBasicInfo() {
    runnum   = 0;
    evnum    = 0;
    status   = 0;
    beamflag = beam::kUnknown;
    dst::clear_all(trigpat, trigflag, clkTpc, cobo_id);
  }

  void clearRawHits() {
    nhTpc = 0;
    dst::clear_all(raw_hitpos_x, raw_hitpos_y, raw_hitpos_z,
                       raw_de, raw_padid, raw_layer, raw_row);
  }

  void clearClusters() {
    nclTpc = 0;
    dst::clear_all(cluster_x, cluster_y, cluster_z, cluster_de,
                       cluster_size, cluster_layer, cluster_mrow,
                       cluster_de_center, cluster_x_center, cluster_y_center,
                       cluster_z_center, cluster_row_center, cluster_houghflag);
  }

  void clearTracks() {
    ntTpc = 0;
    dst::clear_all(
      nhtrack, chisqrTpc, x0Tpc, y0Tpc, u0Tpc, v0Tpc, theta,

      hitlayer, hitpos_x, hitpos_y, hitpos_z,
      calpos_x, calpos_y, calpos_z,

      residual, residual_x, residual_y, residual_z,
      residual_horizontal, residual_vertical,

      resolution_x, resolution_y, resolution_z,
      resolution_horizontal, resolution_vertical,

      dE, dEdx, pathhit, theta_diff,

      track_cluster_de, track_cluster_size, track_cluster_mrow,
      track_cluster_de_center, track_cluster_x_center, track_cluster_y_center,
      track_cluster_z_center, track_cluster_row_center
    );
  }

  void clearExclusiveTracks() {
    dst::clear_all(exresidual, exresidual_x, exresidual_y, exresidual_z,
                       exresidual_horizontal, exresidual_vertical);
  }

  void clearVertex() {
    ntTpc_inside = 0;
    prodvtx_x    = TMath::QuietNaN();
    prodvtx_y    = TMath::QuietNaN();
    prodvtx_z    = TMath::QuietNaN();
  }

  void clearFailedTracks() {
    failed_ntTpc = 0;
    dst::clear_all(
      failed_nhtrack, failed_x0Tpc, failed_y0Tpc, failed_u0Tpc, failed_v0Tpc,
      failed_hitlayer, failed_hitpos_x, failed_hitpos_y, failed_hitpos_z,
      failed_calpos_x, failed_calpos_y, failed_calpos_z
    );
  }

  void clear()
  {
    clearBasicInfo();
    clearRawHits();
    clearClusters();
    clearTracks();
    clearExclusiveTracks();
    clearVertex();
    clearFailedTracks();
  }

  void resizeTracks(Int_t nTracks) {
    dst::resize_all(nTracks,
      nhtrack, chisqrTpc, x0Tpc, y0Tpc, u0Tpc, v0Tpc, theta,

      hitlayer, hitpos_x, hitpos_y, hitpos_z,
      calpos_x, calpos_y, calpos_z,

      residual, residual_x, residual_y, residual_z,
      residual_horizontal, residual_vertical,

      resolution_x, resolution_y, resolution_z,
      resolution_horizontal, resolution_vertical,

      exresidual, exresidual_x, exresidual_y, exresidual_z,
      exresidual_horizontal, exresidual_vertical,

      dE, dEdx, pathhit, theta_diff,

      track_cluster_de, track_cluster_size, track_cluster_mrow,
      track_cluster_de_center, track_cluster_x_center, track_cluster_y_center,
      track_cluster_z_center, track_cluster_row_center
    );
  }

  void resizeTrackHits(Int_t it, Int_t nh) {
    dst::resize_all(nh,
      hitlayer[it],

      hitpos_x[it], hitpos_y[it], hitpos_z[it],
      calpos_x[it], calpos_y[it], calpos_z[it],

      residual[it],
      residual_x[it], residual_y[it], residual_z[it],
      residual_horizontal[it], residual_vertical[it],

      resolution_x[it], resolution_y[it], resolution_z[it],
      resolution_horizontal[it], resolution_vertical[it],

      exresidual[it], exresidual_x[it], exresidual_y[it], exresidual_z[it],
      exresidual_horizontal[it], exresidual_vertical[it],

      pathhit[it], theta_diff[it],

      track_cluster_de[it], track_cluster_size[it], track_cluster_mrow[it],
      track_cluster_de_center[it],
      track_cluster_x_center[it], track_cluster_y_center[it],
      track_cluster_z_center[it], track_cluster_row_center[it]
    );
  }

  void resizeFailedTracks(Int_t failed_ntTpc) {
    dst::resize_all(failed_ntTpc,
      failed_nhtrack,
      failed_x0Tpc, failed_y0Tpc, failed_u0Tpc, failed_v0Tpc,
      failed_hitlayer,
      failed_hitpos_x, failed_hitpos_y, failed_hitpos_z,
      failed_calpos_x, failed_calpos_y, failed_calpos_z
    );
  }

  void resizeFailedTrackHits(Int_t it, Int_t nh)
  {
    dst::resize_all(nh,
      failed_hitlayer[it],
      failed_hitpos_x[it], failed_hitpos_y[it], failed_hitpos_z[it],
      failed_calpos_x[it], failed_calpos_y[it], failed_calpos_z[it]
    );
  }

};

//_____________________________________________________________________________
struct Src
{
  TTreeReaderValue<UInt_t>* runnum;
  TTreeReaderValue<UInt_t>* evnum;
  TTreeReaderValue<std::vector<Double_t>>* trigpat;
  TTreeReaderValue<std::vector<std::vector<Double_t>>>* trigflag;
  TTreeReaderValue<Int_t>* beamflag;
  TTreeReaderValue<Int_t>* npadTpc;   // number of pads
  TTreeReaderValue<Int_t>* nhTpc;     // number of hits
  // vector (size=nhTpc)
  TTreeReaderValue<std::vector<Int_t>>* layerTpc;     // layer id
  TTreeReaderValue<std::vector<Int_t>>* rowTpc;       // row id
  TTreeReaderValue<std::vector<Int_t>>* padTpc;       // pad id
  TTreeReaderValue<std::vector<Double_t>>* pedTpc;    // pedestal
  TTreeReaderValue<std::vector<Double_t>>* rmsTpc;    // rms
  TTreeReaderValue<std::vector<Double_t>>* deTpc;     // dE
  TTreeReaderValue<std::vector<Double_t>>* tTpc;      // time
  TTreeReaderValue<std::vector<Double_t>>* chisqrTpc; // chi^2 of signal fitting
  TTreeReaderValue<std::vector<Double_t>>* clkTpc;    // clock time
  TTreeReaderValue<std::vector<Double_t>>* cobo_id;   // CoBo ID
};

namespace root
{
  Event  event;
  Src    src;
  TTree *tree;

  Double_t
  TranseverseDistance(Double_t x_center, Double_t z_center, Double_t x, Double_t z)
  {
    Double_t dummy = std::hypot(x-x_center, z-z_center);
    Double_t dist;
    if(x_center-x<0) dist=-1.*dummy;
    else dist=dummy;
    return dist;
  }
}

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gTpcParam = TPCParamMan::GetInstance();
const auto& gCounter  = debug::ObjectCounter::GetInstance();
const double TRUNCATED_MEAN = 0.8;  // 80%

// Fill event with run/ev/trigger/CoBo clock; return false to skip (DstRead then returns true).
Bool_t FillEventBasicInfo(Int_t ievent)
{
  // Copy run/event/trigger/CoBo from source to event; increment Status.
  event.runnum   = **src.runnum;
  event.evnum    = **src.evnum;
  event.trigpat  = **src.trigpat;
  event.trigflag = **src.trigflag;
  event.beamflag = **src.beamflag;
  event.clkTpc   = **src.clkTpc;
  event.cobo_id  = **src.cobo_id;
  HF1("Status", event.status++);

  if (**src.nhTpc == 0)
    return false;

  HF1("Status", event.status++);

  // Reject if CoBo clock vector size is wrong.
  if (event.clkTpc.size() != NumOfSegCOBO) {
    spdlog::warn("something is wrong: event.clkTpc.size() != {}", NumOfSegCOBO);
    return false;
  }
  // Reject if any CoBo clock is NaN/Inf.
  std::vector<Int_t> bad_cobo;
  for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
    if (!std::isfinite(event.clkTpc[cobo])) bad_cobo.push_back(cobo);
  }
  if (!bad_cobo.empty()) {
    std::ostringstream oss;
    for (size_t i = 0; i < bad_cobo.size(); ++i) { oss << (i ? "," : "") << bad_cobo[i]; }
    spdlog::warn("CoBo clock(s) missing (NaN/Inf): cobo={}, skip event", oss.str());
    return false;
  }

  HF1("Status", event.status++);
  return true;
}

#if RawHit
// Fill event with raw TPC hits (position, dE, pad, layer, row).
void FillRawHits(TPCAnalyzer& TPCAna)
{
  Int_t nh_tpc = 0;
  // Loop layers and good hits; append position, dE, pad, layer, row to event.
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    auto hit_cont = TPCAna.GetTPCHC(layer);
    for (const auto& hit : hit_cont) {
      if (!hit || !hit->IsGood()) continue;
      const auto& pos = hit->GetPosition();
      event.raw_hitpos_x.push_back(pos.X());
      event.raw_hitpos_y.push_back(pos.Y());
      event.raw_hitpos_z.push_back(pos.Z());
      event.raw_de.push_back(hit->GetCDe());
      event.raw_padid.push_back(hit->GetPad());
      event.raw_layer.push_back(layer);
      event.raw_row.push_back(hit->GetRow());
      ++nh_tpc;
    }
  }
  event.nhTpc = nh_tpc;
  HF1("Status", event.status++);
}
#endif

#if RawCluster
// Fill event with TPC cluster info (position, dE, size, center hit, etc.).
void FillClusters(TPCAnalyzer& TPCAna)
{
  Int_t ncl_tpc = 0;
  // Loop layers and good clusters; append position, dE, size, center, Hough flag to event.
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    auto cl_cont = TPCAna.GetTPCClCont(layer);
    for (const auto& cl : cl_cont) {
      if (!cl || !cl->IsGood()) continue;
      TPCHit* center_hit = cl->GetCenterHit();
      const TVector3& center_pos = center_hit->GetPosition();
      event.cluster_x.push_back(cl->GetX());
      event.cluster_y.push_back(cl->GetY());
      event.cluster_z.push_back(cl->GetZ());
      event.cluster_de.push_back(cl->GetDe());
      event.cluster_size.push_back(cl->GetClusterSize());
      event.cluster_layer.push_back(layer);
      event.cluster_mrow.push_back(cl->MeanRow());
      event.cluster_houghflag.push_back(cl->GetHoughFlag());
      event.cluster_de_center.push_back(center_hit->GetCDe());
      event.cluster_x_center.push_back(center_pos.X());
      event.cluster_y_center.push_back(center_pos.Y());
      event.cluster_z_center.push_back(center_pos.Z());
      event.cluster_row_center.push_back(center_hit->GetRow());
      ++ncl_tpc;
    }
  }
  event.nclTpc = ncl_tpc;
  HF1("Status", event.status++);
}
#endif

// Fill event and histograms for one track hit (residuals, cluster, CoBo time).
void ProcessOneTrackHit(Int_t it, Int_t ih, TPCLTrackHit* hit, TPCLocalTrack* track,
                        TPCEventAnalyzer& event_ana)
{
  // Hough distance histos.
  HF1("TPCTrk_Hough_Dist",  hit->GetHoughDist());
  HF1("TPCTrk_Hough_DistY", hit->GetHoughDistY());

  // Extract layer, positions, residuals, cluster/center, path length.
  Int_t layer = hit->GetLayer();
  const TVector3& hit_pos   = hit->GetLocalHitPos();
  const TVector3& cal_pos   = hit->GetLocalCalPos();
  const TVector3& resi_vect = hit->GetResidualVect();
  const TVector3& res_vect  = hit->GetResolutionVect();

  TPCHit* cl_hit     = hit->GetHit();
  TPCCluster* cl     = cl_hit->GetParentCluster();
  Int_t cl_size      = cl->GetClusterSize();
  Double_t clde      = cl->GetDe();
  Double_t mrow      = cl->MeanRow();
  TPCHit* center_hit = cl->GetCenterHit();
  const TVector3& center_pos = center_hit->GetPosition();
  Double_t center_de = center_hit->GetCDe();
  Int_t center_row   = center_hit->GetRow();
  Double_t residual  = hit->GetResidual();
  Double_t hit_length = track->GetHitLength(ih);

  // Fill event: cluster, hit/cal positions, residuals, resolution, path.
  event.track_cluster_de[it][ih]   = clde;
  event.track_cluster_size[it][ih] = cl_size;
  event.track_cluster_mrow[it][ih] = mrow;
  event.track_cluster_de_center[it][ih]  = center_de;
  event.track_cluster_x_center[it][ih]   = center_pos.X();
  event.track_cluster_y_center[it][ih]   = center_pos.Y();
  event.track_cluster_z_center[it][ih]   = center_pos.Z();
  event.track_cluster_row_center[it][ih] = center_row;
  event.hitlayer[it][ih] = layer;
  event.hitpos_x[it][ih] = hit_pos.x();
  event.hitpos_y[it][ih] = hit_pos.y();
  event.hitpos_z[it][ih] = hit_pos.z();
  event.calpos_x[it][ih] = cal_pos.x();
  event.calpos_y[it][ih] = cal_pos.y();
  event.calpos_z[it][ih] = cal_pos.z();
  event.residual[it][ih]   = residual;
  event.residual_x[it][ih] = resi_vect.x();
  event.residual_y[it][ih] = resi_vect.y();
  event.residual_z[it][ih] = resi_vect.z();
  event.residual_horizontal[it][ih] = track->GetHorizontalResidual(ih);
  event.residual_vertical[it][ih]   = track->GetVerticalResidual(ih);
  event.resolution_x[it][ih] = res_vect.x();
  event.resolution_y[it][ih] = res_vect.y();
  event.resolution_z[it][ih] = res_vect.z();
  event.resolution_horizontal[it][ih] = track->GetHorizontalResolution(ih);
  event.resolution_vertical[it][ih]   = track->GetVerticalResolution(ih);
  event.theta_diff[it][ih] = track->GetAlpha(ih);
  event.pathhit[it][ih]    = hit_length;

  // Track/hit pattern and residual histos (layer, row, pad, res, res vs x/y).
  HF1("TPCTrk_Layer", layer);
  HF2("TPCTrk_Row_vs_Layer", layer, center_row);
  HF1(Form("TPCHit_HitPat_Layer%02d", layer), center_row);
  Int_t pad_id = tpc::GetPadId(layer, center_row);
  if (pad_id >= 0) {
    Double_t bin_cont = HG2Poly("TPCTrk_HitPat", pad_id + 1);
    HF2Poly("TPCTrk_HitPat", pad_id + 1, bin_cont + 1.);
  }
  HF1(Form("TPCHit_Xhit_Layer%02d", layer), hit_pos.x());
  HF1(Form("TPCTrk_Res_Layer%02d", layer), residual);
  HF2(Form("TPCTrk_Res_vs_Xhit_Layer%02d", layer), hit_pos.x(), residual);
  HF2(Form("TPCHit_Yhit_vs_Xtrk_Layer%02d", layer), cal_pos.x(), hit_pos.y());
  HF1(Form("TPCTrk_ResX_Layer%02d", layer), resi_vect.X());
  HF1(Form("TPCTrk_ResY_Layer%02d", layer), resi_vect.Y());
  HF1(Form("TPCTrk_ResZ_Layer%02d", layer), resi_vect.Z());
  HF1(Form("TPCTrk_ResY_Layer%02d_Row%03d", layer, center_row), resi_vect.Y());
  HF2(Form("TPCTrk_ResY_vs_Y_Layer%02d_Row%03d", layer, center_row), cal_pos.y(), resi_vect.Y());
  HF2("TPCTrk_ResY_vs_Layer_Trk", layer, resi_vect.Y());

  // Cluster size/dE histos; dE ratio vs transverse distance within cluster.
  HF1("TPCCl_Size", cl_size);
  HF1(Form("TPCCl_Size_Layer%02d", layer), cl_size);
  HF1("TPCCl_dE", clde);
  HF1(Form("TPCCl_dE_Layer%02d", layer), clde);
  const TPCHitContainer& hit_cont = cl->GetHitContainer();
  for (const auto& hits : hit_cont) {
    if (!hits || !hits->IsGood()) continue;
    const TVector3& pos = hits->GetPosition();
    Double_t de = hits->GetCDe();
    Double_t trans_dist = TranseverseDistance(hit_pos.x(), hit_pos.z(), pos.x(), pos.z());
    Double_t ratio = de / clde;
    HF2("TPCCl_Ratio_vs_Dist_Diff", trans_dist, ratio);
    HF2(Form("TPCCl_Ratio_vs_Dist_Diff_Layer%02d", layer), trans_dist, ratio);
  }

  // CoBo clock time vs position (event_ana).
  if (center_hit && center_hit->GetCTimeSize() > 0) {
    TVector3 local_cal = hit->GetLocalCalPos();
    ThreeVector local_trk(local_cal.X(), local_cal.Y(), local_cal.Z());
    Double_t y_trk_global = gGeom.Local2GlobalPos("HypTPC", local_trk).y();
    event_ana.FillCoBoClockTime("TPCTrk", layer, center_row,
                                center_hit->GetCTime(0), center_hit->GetPosition(), y_trk_global);
  }

#if Exclusive
  // Exclusive residuals: fill event.exresidual*.
  const TVector3& exres_vect = hit->GetResidualVectExclusive();
  event.exresidual[it][ih]   = hit->GetResidualExclusive();
  event.exresidual_x[it][ih] = exres_vect.x();
  event.exresidual_y[it][ih] = exres_vect.y();
  event.exresidual_z[it][ih] = exres_vect.z();
  event.exresidual_horizontal[it][ih] = track->GetHorizontalResidualExclusive(ih);
  event.exresidual_vertical[it][ih]   = track->GetVerticalResidualExclusive(ih);
#endif
}

// Fill event and histograms for one track (params, vertex candidates, hits, dE/dx).
void ProcessOneTrack(Int_t it, TPCLocalTrack* track,
                     Int_t& ntrack_intarget,
                     std::vector<Double_t>& x0_vtx, std::vector<Double_t>& y0_vtx,
                     std::vector<Double_t>& u0_vtx, std::vector<Double_t>& v0_vtx,
                     TPCEventAnalyzer& event_ana)
{
  // Track parameters and vertex candidate (in-target x0,y0,u0,v0).
  Int_t nhits = track->GetNHit();
  Double_t chisqr = track->GetChiSquare();
  Double_t x0 = track->GetX0(), y0 = track->GetY0();
  Double_t u0 = track->GetU0(), v0 = track->GetV0();
  Double_t theta = track->GetTheta();

  if (TMath::Abs(x0) < 50. && TMath::Abs(y0) < 50.) {
    x0_vtx.push_back(x0);
    y0_vtx.push_back(y0);
    u0_vtx.push_back(u0);
    v0_vtx.push_back(v0);
    ++ntrack_intarget;
  }

  event.nhtrack[it]   = nhits;
  event.chisqrTpc[it] = chisqr;
  event.x0Tpc[it]     = x0;
  event.y0Tpc[it]     = y0;
  event.u0Tpc[it]     = u0;
  event.v0Tpc[it]     = v0;
  event.theta[it]     = theta;
  event.resizeTrackHits(it, nhits);

  // Track-parameter histos.
  HF1("TPCTrk_Num_TrackHits", nhits);
  HF1("TPCTrk_Chisqr", chisqr);
  HF1("TPCTrk_X0", x0);
  HF1("TPCTrk_Y0", y0);
  HF1("TPCTrk_U0", u0);
  HF1("TPCTrk_V0", v0);
  HF2("TPCTrk_U0_vs_X0", x0, u0);
  HF2("TPCTrk_V0_vs_Y0", y0, v0);
  HF2("TPCTrk_Y0_vs_X0", x0, y0);

  // Fit/search metadata histos (iterations, flags, times, Minuit status).
  Int_t niter          = track->GetNIteration();
  Int_t fit_flag       = track->GetFitFlag();
  Int_t search_time    = track->GetSearchTime();
  Int_t fit_time       = track->GetFitTime();
  Int_t minuit_status  = track->GetMinuitStatus();
  HF1("TPCTrk_Num_Iter",     niter);
  HF1("TPCTrk_Fitting_Flag", fit_flag);
  HF1("TPCTrk_Searching_Time", search_time);
  HF1("TPCTrk_Fitting_Time",   fit_time);
  HF1("TPCTrk_Minuit_Status",  minuit_status);

  // Per-hit processing, total dE, and dE/dx per hit for truncated mean.
  Double_t total_clde = 0.;
  std::vector<Double_t> dedx_vect;
  for (Int_t ih = 0; ih < nhits; ++ih) {
    auto trk_hit = track->GetHit(ih);
    if (!trk_hit) continue;
    ProcessOneTrackHit(it, ih, trk_hit, track, event_ana);
    TPCCluster* cl = trk_hit->GetHit()->GetParentCluster();
    total_clde += cl->GetDe();
    dedx_vect.push_back(cl->GetDe() / track->GetHitLength(ih));
  }

  // Truncated-mean dE/dx -> event.dEdx.
  event.dE[it] = total_clde;
  std::sort(dedx_vect.begin(), dedx_vect.end());
  Int_t n_truncated = static_cast<Int_t>(std::floor(dedx_vect.size() * TRUNCATED_MEAN));
  if (n_truncated < 1) {
    spdlog::warn("something is wrong: n_truncated = {}", n_truncated);
  }
  Double_t sum_truncated_de = 0.;
  for (Int_t i = 0; i < n_truncated; ++i) {
    sum_truncated_de += dedx_vect[i];
  }
  event.dEdx[it] = sum_truncated_de / n_truncated;
}

// Compute multi-track vertex from in-target tracks and fill event.
void FillVertex(Int_t ntrack_intarget,
                const std::vector<Double_t>& x0_vtx, const std::vector<Double_t>& y0_vtx,
                const std::vector<Double_t>& u0_vtx, const std::vector<Double_t>& v0_vtx)
{
  // Multitrack vertex; fill ntTpc_inside, prodvtx_*, Status.
  TVector3 vertex = Kinematics::MultitrackVertex(
    ntrack_intarget, x0_vtx, y0_vtx, u0_vtx, v0_vtx);
  event.ntTpc_inside = ntrack_intarget;
  event.prodvtx_x = vertex.x();
  event.prodvtx_y = vertex.y();
  event.prodvtx_z = vertex.z();
  HF1("Status", event.status++);
}

#if TrackSearchFailed
// Fill event with failed track-search candidates (count and params).
void FillFailedTracks(TPCAnalyzer& TPCAna)
{
  Int_t failed_nt_tpc = TPCAna.GetNTracksTPCFailed();
  event.failed_ntTpc = failed_nt_tpc;
  event.resizeFailedTracks(failed_nt_tpc);

  // Per failed track: params and hit positions.
  for (Int_t it = 0; it < failed_nt_tpc; ++it) {
    auto track = TPCAna.GetTrackTPCFailed(it);
    if (!track) continue;
    Int_t nhits = track->GetNHit();
    event.failed_nhtrack[it] = nhits;
    event.failed_x0Tpc[it] = track->GetX0();
    event.failed_y0Tpc[it] = track->GetY0();
    event.failed_u0Tpc[it] = track->GetU0();
    event.failed_v0Tpc[it] = track->GetV0();
    event.resizeFailedTrackHits(it, nhits);

    for (Int_t ih = 0; ih < nhits; ++ih) {
      auto hit = track->GetHit(ih);
      if (!hit) continue;
      event.failed_hitlayer[it][ih] = hit->GetLayer();
      const TVector3& hit_pos = hit->GetLocalHitPos();
      const TVector3& cal_pos = hit->GetLocalCalPos();
      event.failed_hitpos_x[it][ih] = hit_pos.x();
      event.failed_hitpos_y[it][ih] = hit_pos.y();
      event.failed_hitpos_z[it][ih] = hit_pos.z();
      event.failed_calpos_x[it][ih] = cal_pos.x();
      event.failed_calpos_y[it][ih] = cal_pos.y();
      event.failed_calpos_z[it][ih] = cal_pos.z();
    }
  }
  HF1("Status", event.status++);
}
#endif
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kTpcHit,  kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[TPCHit]",  "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
Bool_t SetupReader();
}

//_____________________________________________________________________________
int
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeHistograms())
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;
  if(!dst::SetupReader())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if(skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if(max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
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
  for(const auto& name : TreeName) if(name != "") n_input_files++;

  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(Int_t i=0; i<nArgc; ++i){
    if(TreeName[i] == "") continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=n_input_files || open_tree!=n_input_files){
    spdlog::error("DstOpen Failed: opened files/trees mismatch based on TreeName definitions."
                  " expected: {}, open_file: {}, open_tree: {}",
                  n_input_files, open_file, open_tree);
    return false;
  }
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  if (ievent % 100 == 0) {
    std::cout << "#D Event Number: " << std::setw(6) << ievent << std::endl;
  }
  GetEntry(ievent);

  // Fill basic info and validate CoBo clocks; skip further processing when validation fails.
  if (!FillEventBasicInfo(ievent))
    return true;

  // Recalc TPC hits and set up event analyzer clock.
  TPCAnalyzer tpc_ana;
  tpc_ana.ReCalcTPCHits(**src.nhTpc, **src.padTpc, **src.tTpc, **src.deTpc, **src.clkTpc);
  HF1("Status", event.status++);

  static TPCEventAnalyzer event_ana;
  event_ana.SetClock(**src.clkTpc);

  // Track search (exclusive or standard).
#if Exclusive
  tpc_ana.TrackSearchTPC(true);
#else
  tpc_ana.TrackSearchTPC();
#endif
  HF1("Status", event.status++);

  // Optional: raw hits and/or clusters into event.
#if RawHit
  FillRawHits(tpc_ana);
#endif

#if RawCluster
  FillClusters(tpc_ana);
#endif

  // Track count; skip if no tracks.
  Int_t nt_tpc = tpc_ana.GetNTracksTPC();
  event.ntTpc = nt_tpc;
  HF1("TPCTrk_Num_Track", nt_tpc);
  if (event.ntTpc == 0)
    return true;

  HF1("Status", event.status++);
  event.resizeTracks(nt_tpc);

  // Per-track fill (params, vertex candidates, hits, dE/dx).
  Int_t ntrack_intarget = 0;
  std::vector<Double_t> x0_vtx, y0_vtx, u0_vtx, v0_vtx;
  x0_vtx.reserve(100);
  y0_vtx.reserve(100);
  u0_vtx.reserve(100);
  v0_vtx.reserve(100);

  for (Int_t it = 0; it < nt_tpc; ++it) {
    auto track = tpc_ana.GetTrackTPC(it);
    if (!track) continue;
    ProcessOneTrack(it, track, ntrack_intarget, x0_vtx, y0_vtx, u0_vtx, v0_vtx, event_ana);
  }
  HF1("Status", event.status++);

  // Multi-track vertex; optional failed-track info.
  FillVertex(ntrack_intarget, x0_vtx, y0_vtx, u0_vtx, v0_vtx);

#if TrackSearchFailed
  FillFailedTracks(tpc_ana);
#endif

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
  for(Int_t i=0; i<n; ++i){
    if(TTreeReaderCont[i]) delete TTreeReaderCont[i];
    if(TTreeCont[i])       delete TTreeCont[i];
    if(TFileCont[i])       delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
dst::SetupReader()
{
  if (!dst::SetupReader(kTpcHit, "kTpcHit")) return false;

  const auto& reader = TTreeReaderCont[kTpcHit];

  dst::SetBranch(reader, "run_number",   src.runnum);
  dst::SetBranch(reader, "event_number", src.evnum);
  dst::SetBranch(reader, "trig_pat",     src.trigpat);
  dst::SetBranch(reader, "trig_flag",    src.trigflag);
  dst::SetBranch(reader, "beam_flag",    src.beamflag);
  dst::SetBranch(reader, "npadTpc",      src.npadTpc);
  dst::SetBranch(reader, "nhTpc",        src.nhTpc);
  dst::SetBranch(reader, "layerTpc",     src.layerTpc);
  dst::SetBranch(reader, "rowTpc",       src.rowTpc);
  dst::SetBranch(reader, "padTpc",       src.padTpc);
  dst::SetBranch(reader, "pedTpc",       src.pedTpc);
  dst::SetBranch(reader, "rmsTpc",       src.rmsTpc);
  dst::SetBranch(reader, "deTpc",        src.deTpc);
  dst::SetBranch(reader, "tTpc",         src.tTpc);
  dst::SetBranch(reader, "chisqrTpc",    src.chisqrTpc);
  dst::SetBranch(reader, "clkTpc",       src.clkTpc);
  dst::SetBranch(reader, "cobo_id",      src.cobo_id);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  hist::BuildStatus();
  hist::BuildTPCBasic();
  hist::BuildTPCTracking();

  tree = new TTree("tpc", "tree of DstTPCTracking");
  tree->Branch("status", &event.status);
  tree->Branch("run_number", &event.runnum);
  tree->Branch("event_number", &event.evnum);
  tree->Branch("trig_pat", &event.trigpat);
  tree->Branch("trig_flag", &event.trigflag);
  tree->Branch("beam_flag", &event.beamflag);
  tree->Branch("clkTpc", &event.clkTpc);
  tree->Branch("cobo_id", &event.cobo_id);

#if RawHit
  tree->Branch("nhTpc", &event.nhTpc);
  tree->Branch("raw_hitpos_x", &event.raw_hitpos_x);
  tree->Branch("raw_hitpos_y", &event.raw_hitpos_y);
  tree->Branch("raw_hitpos_z", &event.raw_hitpos_z);
  tree->Branch("raw_de", &event.raw_de);
  tree->Branch("raw_padid", &event.raw_padid);
  tree->Branch("raw_layer", &event.raw_layer);
  tree->Branch("raw_row", &event.raw_row);
#endif

#if RawCluster
  tree->Branch("nclTpc", &event.nclTpc);
  tree->Branch("cluster_x", &event.cluster_x);
  tree->Branch("cluster_y", &event.cluster_y);
  tree->Branch("cluster_z", &event.cluster_z);
  tree->Branch("cluster_de", &event.cluster_de);
  tree->Branch("cluster_size", &event.cluster_size);
  tree->Branch("cluster_layer", &event.cluster_layer);
  tree->Branch("cluster_row_center", &event.cluster_row_center);
  tree->Branch("cluster_mrow", &event.cluster_mrow);
  tree->Branch("cluster_houghflag", &event.cluster_houghflag);
  tree->Branch("cluster_de_center", &event.cluster_de_center);
  tree->Branch("cluster_x_center", &event.cluster_x_center);
  tree->Branch("cluster_y_center", &event.cluster_y_center);
  tree->Branch("cluster_z_center", &event.cluster_z_center);
#endif

  tree->Branch("ntTpc", &event.ntTpc);
  tree->Branch("nhtrack", &event.nhtrack);
  tree->Branch("chisqrTpc", &event.chisqrTpc);
  tree->Branch("x0Tpc", &event.x0Tpc);
  tree->Branch("y0Tpc", &event.y0Tpc);
  tree->Branch("u0Tpc", &event.u0Tpc);
  tree->Branch("v0Tpc", &event.v0Tpc);
  tree->Branch("theta", &event.theta);
  tree->Branch("hitlayer", &event.hitlayer);
  tree->Branch("hitpos_x", &event.hitpos_x);
  tree->Branch("hitpos_y", &event.hitpos_y);
  tree->Branch("hitpos_z", &event.hitpos_z);
  tree->Branch("calpos_x", &event.calpos_x);
  tree->Branch("calpos_y", &event.calpos_y);
  tree->Branch("calpos_z", &event.calpos_z);
  tree->Branch("residual", &event.residual);
  tree->Branch("residual_x", &event.residual_x);
  tree->Branch("residual_y", &event.residual_y);
  tree->Branch("residual_z", &event.residual_z);
  tree->Branch("residual_horizontal", &event.residual_horizontal);
  tree->Branch("residual_vertical", &event.residual_vertical);
  tree->Branch("resolution_x", &event.resolution_x);
  tree->Branch("resolution_y", &event.resolution_y);
  tree->Branch("resolution_z", &event.resolution_z);
  tree->Branch("resolution_horizontal", &event.resolution_horizontal);
  tree->Branch("resolution_vertical", &event.resolution_vertical);
  tree->Branch("dE", &event.dE);
  tree->Branch("dEdx", &event.dEdx);
  tree->Branch("pathhit", &event.pathhit);
  tree->Branch("theta_diff", &event.theta_diff);
  tree->Branch("track_cluster_de", &event.track_cluster_de);
  tree->Branch("track_cluster_size", &event.track_cluster_size);
  tree->Branch("track_cluster_mrow", &event.track_cluster_mrow);
  tree->Branch("track_cluster_de_center", &event.track_cluster_de_center);
  tree->Branch("track_cluster_x_center", &event.track_cluster_x_center);
  tree->Branch("track_cluster_y_center", &event.track_cluster_y_center);
  tree->Branch("track_cluster_z_center", &event.track_cluster_z_center);
  tree->Branch("track_cluster_row_center", &event.track_cluster_row_center);
  tree->Branch("ntTpc_target", &event.ntTpc_inside);
  tree->Branch("prodvtx_x", &event.prodvtx_x);
  tree->Branch("prodvtx_y", &event.prodvtx_y);
  tree->Branch("prodvtx_z", &event.prodvtx_z);

#if Exclusive
  tree->Branch("exresidual", &event.exresidual);
  tree->Branch("exresidual_x", &event.exresidual_x);
  tree->Branch("exresidual_y", &event.exresidual_y);
  tree->Branch("exresidual_z", &event.exresidual_z);
  tree->Branch("exresidual_horizontal", &event.exresidual_horizontal);
  tree->Branch("exresidual_vertical", &event.exresidual_vertical);
#endif

#if TrackSearchFailed
  tree->Branch("failed_ntTpc", &event.failed_ntTpc);
  tree->Branch("failed_nhtrack", &event.failed_nhtrack);
  tree->Branch("failed_x0Tpc", &event.failed_x0Tpc);
  tree->Branch("failed_y0Tpc", &event.failed_y0Tpc);
  tree->Branch("failed_u0Tpc", &event.failed_u0Tpc);
  tree->Branch("failed_v0Tpc", &event.failed_v0Tpc);
  tree->Branch("failed_hitlayer", &event.failed_hitlayer);
  tree->Branch("failed_hitpos_x", &event.failed_hitpos_x);
  tree->Branch("failed_hitpos_y", &event.failed_hitpos_y);
  tree->Branch("failed_hitpos_z", &event.failed_hitpos_z);
  tree->Branch("failed_calpos_x", &event.failed_calpos_x);
  tree->Branch("failed_calpos_y", &event.failed_calpos_y);
  tree->Branch("failed_calpos_z", &event.failed_calpos_z);
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<TPCParamMan>("TPCPRM", "TPCPHASE") &&
     InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
