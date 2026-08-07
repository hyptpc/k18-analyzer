// -*- C++ -*-

#include "VEvent.hh"

#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

#include <TTree.h>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCTdcCalibMan.hh"
#include "D5Track.hh"
#include "D5TransMatrix.hh"
#include "EventAnalyzer.hh"
#include "HistTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

#define DEBUG_D5 0

namespace
{
using hddaq::unpacker::GUnpacker;
auto& gD5Mtx = D5TransMatrix::GetInstance();
auto& gUnpacker = GUnpacker::get_instance();

TTree* tree;

struct Event
{
  UInt_t run_number;
  UInt_t event_number;
  beam::EBeamFlag beam_flag;
  Int_t status;

  Int_t ntrack_blc1;
  Int_t ntrack_blc2;
  Int_t ntrack_d5;

  Double_t d5_z_in;
  Double_t d5_z_out;

  std::vector<Double_t> blc1_x;
  std::vector<Double_t> blc1_y;
  std::vector<Double_t> blc1_u;
  std::vector<Double_t> blc1_v;
  std::vector<Double_t> blc1_chi2;

  std::vector<Double_t> blc2_x;
  std::vector<Double_t> blc2_y;
  std::vector<Double_t> blc2_u;
  std::vector<Double_t> blc2_v;
  std::vector<Double_t> blc2_chi2;

  std::vector<Int_t>    d5_blc1_trk;
  std::vector<Int_t>    d5_blc2_trk;
  std::vector<Double_t> d5_delta;
  std::vector<Double_t> d5_momentum;
  std::vector<Double_t> d5_fit_chi2;
  std::vector<Int_t>    d5_fit_ndf;
  std::vector<Double_t> d5_fit_chi2_ndf;
  std::vector<Int_t>    d5_pair_rank;
  std::vector<Double_t> d5_residual_x;
  std::vector<Double_t> d5_residual_y;
  std::vector<Double_t> d5_blc2out_x;
  std::vector<Double_t> d5_blc2out_y;
  std::vector<Double_t> d5_mtxout_x;
  std::vector<Double_t> d5_mtxout_u;
  std::vector<Double_t> d5_mtxout_y;
  std::vector<Double_t> d5_mtxout_v;
  std::vector<Double_t> d5_fit_x0;
  std::vector<Double_t> d5_fit_u0;
  std::vector<Double_t> d5_fit_y0;
  std::vector<Double_t> d5_fit_v0;
};
Event event;

constexpr Double_t MaxBlcChi2Ndf = 100.0;

Double_t
BlcChi2Ndf(const DCLocalTrack* tr)
{
  if (!tr) return std::numeric_limits<Double_t>::infinity();
  const Int_t ndf = tr->GetNDF();
  if (ndf <= 0) return std::numeric_limits<Double_t>::infinity();
  return tr->GetChiSquare() / static_cast<Double_t>(ndf);
}

void
AssignD5PairRanks()
{
  const Int_t n = static_cast<Int_t>(event.d5_fit_chi2_ndf.size());
  event.d5_pair_rank.assign(n, 0);
  if (n <= 0) return;

  std::vector<Int_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [](Int_t a, Int_t b) {
              const Double_t sa = event.d5_fit_chi2_ndf[a];
              const Double_t sb = event.d5_fit_chi2_ndf[b];
              if (sa != sb) return sa < sb;
              return a < b;
            });
  for (Int_t r=0; r<n; ++r)
    event.d5_pair_rank[order[r]] = r;
}

void
PushD5Result(const D5Track& d5tr, Int_t i, Int_t j, const DCLocalTrack* tr2)
{
  event.d5_blc1_trk.push_back(i);
  event.d5_blc2_trk.push_back(j);
  event.d5_delta.push_back(d5tr.GetDelta());
  event.d5_momentum.push_back(d5tr.GetMomentum());
  event.d5_fit_chi2.push_back(d5tr.GetD5Chi2());
  event.d5_fit_ndf.push_back(d5tr.GetD5Ndf());
  event.d5_fit_chi2_ndf.push_back(d5tr.GetD5Chi2Ndf());
  event.d5_residual_x.push_back(d5tr.GetResidualX());
  event.d5_residual_y.push_back(d5tr.GetResidualY());
  event.d5_blc2out_x.push_back(tr2->GetX(event.d5_z_out));
  event.d5_blc2out_y.push_back(tr2->GetY(event.d5_z_out));
  event.d5_mtxout_x.push_back(d5tr.GetMtxoutX());
  event.d5_mtxout_u.push_back(d5tr.GetMtxoutU());
  event.d5_mtxout_y.push_back(d5tr.GetMtxoutY());
  event.d5_mtxout_v.push_back(d5tr.GetMtxoutV());
  event.d5_fit_x0.push_back(d5tr.GetFitX0());
  event.d5_fit_u0.push_back(d5tr.GetFitU0());
  event.d5_fit_y0.push_back(d5tr.GetFitY0());
  event.d5_fit_v0.push_back(d5tr.GetFitV0());
}
}

//_____________________________________________________________________________
Bool_t
ProcessBegin()
{
  event.run_number = gUnpacker.get_run_number();
  event.event_number = gUnpacker.get_event_number();
  event.beam_flag = beam::kUnknown;
  event.status = 0;

  event.ntrack_blc1 = 0;
  event.ntrack_blc2 = 0;
  event.ntrack_d5 = 0;

  event.blc1_x.clear(); event.blc1_y.clear(); event.blc1_u.clear();
  event.blc1_v.clear(); event.blc1_chi2.clear();
  event.blc2_x.clear(); event.blc2_y.clear(); event.blc2_u.clear();
  event.blc2_v.clear(); event.blc2_chi2.clear();
  event.d5_blc1_trk.clear(); event.d5_blc2_trk.clear();
  event.d5_delta.clear(); event.d5_momentum.clear();
  event.d5_fit_chi2.clear(); event.d5_fit_ndf.clear();
  event.d5_fit_chi2_ndf.clear();
  event.d5_pair_rank.clear();
  event.d5_residual_x.clear(); event.d5_residual_y.clear();
  event.d5_blc2out_x.clear(); event.d5_blc2out_y.clear();
  event.d5_mtxout_x.clear(); event.d5_mtxout_u.clear();
  event.d5_mtxout_y.clear(); event.d5_mtxout_v.clear();
  event.d5_fit_x0.clear(); event.d5_fit_u0.clear();
  event.d5_fit_y0.clear(); event.d5_fit_v0.clear();

  event.d5_z_in  = gD5Mtx.IsRefPlaneReady() ? gD5Mtx.GetD5ZIn()  : D5TransMatrix::RefZIn;
  event.d5_z_out = gD5Mtx.IsRefPlaneReady() ? gD5Mtx.GetD5ZOut() : D5TransMatrix::RefZOut;

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessNormal()
{
  using root::HF1;

  RawData rawData;
  for (const auto& name : DCNameList.at("BcIn"))  rawData.DecodeHits(name);
  for (const auto& name : DCNameList.at("BcOut")) rawData.DecodeHits(name);

  EventAnalyzer evAna;
  HF1("Status", event.status++); // BLC raw decoded

  rawData.DecodeHits("TriggerFlag");
  evAna.TriggerFlag(rawData);
  HF1("Status", event.status++); // TriggerFlag

  rawData.DecodeHits("BAC");
  rawData.DecodeHits("BHT");
  event.beam_flag = evAna.BeamFlag(rawData);
  HF1("Status", event.status++); // BTOF / beam flag

  DCAnalyzer dcAna(rawData);
  dcAna.DecodeBcInHits();
  dcAna.DecodeBcOutHits();
  HF1("Status", event.status++); // DC hits decoded

  dcAna.TotCut("BLC1a"); dcAna.TotCut("BLC1b");
  dcAna.TotCut("BLC2a"); dcAna.TotCut("BLC2b");
  dcAna.DriftTimeCut("BLC1a"); dcAna.DriftTimeCut("BLC1b");
  dcAna.DriftTimeCut("BLC2a"); dcAna.DriftTimeCut("BLC2b");
  HF1("Status", event.status++); // DC hit cuts

  dcAna.TrackSearchBcIn();
  dcAna.TrackSearchBcOut();
  evAna.BcInTracking(dcAna, event.beam_flag);
  evAna.BcOutTracking(dcAna, event.beam_flag);
  HF1("Status", event.status++); // BLC track search

  Int_t nt1 = dcAna.GetNtracksBcIn();
  Int_t nt2 = dcAna.GetNtracksBcOut();

  event.ntrack_blc1 = nt1;
  event.ntrack_blc2 = nt2;
  event.ntrack_d5 = 0;

  for (Int_t i=0; i<nt1; ++i) {
    const auto* tr = dcAna.GetTrackBcIn(i);
    event.blc1_x.push_back(tr->GetX0());
    event.blc1_y.push_back(tr->GetY0());
    event.blc1_u.push_back(tr->GetU0());
    event.blc1_v.push_back(tr->GetV0());
    event.blc1_chi2.push_back(tr->GetChiSquare());
  }
  for (Int_t i=0; i<nt2; ++i) {
    const auto* tr = dcAna.GetTrackBcOut(i);
    event.blc2_x.push_back(tr->GetX0());
    event.blc2_y.push_back(tr->GetY0());
    event.blc2_u.push_back(tr->GetU0());
    event.blc2_v.push_back(tr->GetV0());
    event.blc2_chi2.push_back(tr->GetChiSquare());
  }

  if (nt1 == 0 || nt2 == 0) {
    tree->Fill();
    return true;
  }
  HF1("Status", event.status++); // BLC1 & BLC2 both tracked

  if (!gD5Mtx.IsReady()) {
    tree->Fill();
    return true;
  }

  std::vector<Int_t> good_blc1;
  std::vector<Int_t> good_blc2;
  for (Int_t i=0; i<nt1; ++i) {
    if (BlcChi2Ndf(dcAna.GetTrackBcIn(i)) < MaxBlcChi2Ndf)
      good_blc1.push_back(i);
  }
  for (Int_t j=0; j<nt2; ++j) {
    if (BlcChi2Ndf(dcAna.GetTrackBcOut(j)) < MaxBlcChi2Ndf)
      good_blc2.push_back(j);
  }

  std::vector<std::pair<Int_t, Int_t>> d5_pairs;

  if (!good_blc1.empty() && !good_blc2.empty()) {
    for (Int_t i : good_blc1) {
      for (Int_t j : good_blc2) {
        const auto* tr1 = dcAna.GetTrackBcIn(i);
        const auto* tr2 = dcAna.GetTrackBcOut(j);

        D5Track d5tr(tr1, tr2);
        if (!d5tr.CalcMomentum()) continue;

        PushD5Result(d5tr, i, j, tr2);
        d5_pairs.emplace_back(i, j);

#if DEBUG_D5
        std::cout << "[D5 Debug] Event: " << event.event_number
                  << " | BLC1 Trk: " << i << " -> BLC2 Trk: " << j << "\n"
                  << "   - D5ZIn: " << event.d5_z_in << " | D5ZOut: " << event.d5_z_out << "\n"
                  << "   - BLC1 (x, u) at D5ZIn: ("
                  << tr1->GetX(event.d5_z_in) << ", " << tr1->GetU0()*1000.0 << " mrad)\n"
                  << "   - BLC2 (x, u) at D5ZOut: ("
                  << tr2->GetX(event.d5_z_out) << ", " << tr2->GetU0()*1000.0 << " mrad)\n"
                  << "   - mtxout (x, u): ("
                  << d5tr.GetMtxoutX() << ", " << d5tr.GetMtxoutU() << ")\n"
                  << "   - Reconstructed Delta: " << d5tr.GetDelta() << " %\n"
                  << "   - Momentum: " << d5tr.GetMomentum() << " GeV/c\n"
                  << "   - Chi2/ndf: " << d5tr.GetD5Chi2() << " / "
                  << d5tr.GetD5Ndf() << std::endl;
#endif

        evAna.D5Tracking(d5tr);
      }
    }
    AssignD5PairRanks();

    for (size_t k = 0; k < d5_pairs.size(); ++k) {
      if (event.d5_pair_rank[k] != 0) continue;
      const Int_t i = d5_pairs[k].first;
      const Int_t j = d5_pairs[k].second;
      D5Track d5tr(dcAna.GetTrackBcIn(i), dcAna.GetTrackBcOut(j));
      if (!d5tr.CalcMomentum()) continue;
      evAna.D5WireResiduals(d5tr, event.d5_z_out);
    }

    HF1("Status", event.status++); // D5 fit loop done
  }
  event.ntrack_d5 = static_cast<Int_t>(event.d5_fit_chi2_ndf.size());

  tree->Fill();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessEnd()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  hist::BuildStatus();
  hist::BuildTriggerFlag();
  hist::BuildDCRaw("BcIn", true);
  hist::BuildDCHit("BcIn", true);
  hist::BuildDCTrack("BcIn", true);
  hist::BuildDCRaw("BcOut", true);
  hist::BuildDCHit("BcOut", true);
  hist::BuildDCTrack("BcOut", true);
  hist::BuildD5Tracking(false);

  tree = new TTree("d5", "D5 tracking event data");
  tree->Branch("run_number", &event.run_number);
  tree->Branch("event_number", &event.event_number);
  tree->Branch("beam_flag", &event.beam_flag, "beam_flag/I");
  tree->Branch("status", &event.status);

  tree->Branch("ntrack_blc1", &event.ntrack_blc1);
  tree->Branch("ntrack_blc2", &event.ntrack_blc2);
  tree->Branch("ntrack_d5", &event.ntrack_d5);
  tree->Branch("d5_z_in", &event.d5_z_in);
  tree->Branch("d5_z_out", &event.d5_z_out);

  tree->Branch("blc1_x", &event.blc1_x);
  tree->Branch("blc1_y", &event.blc1_y);
  tree->Branch("blc1_u", &event.blc1_u);
  tree->Branch("blc1_v", &event.blc1_v);
  tree->Branch("blc1_chi2", &event.blc1_chi2);

  tree->Branch("blc2_x", &event.blc2_x);
  tree->Branch("blc2_y", &event.blc2_y);
  tree->Branch("blc2_u", &event.blc2_u);
  tree->Branch("blc2_v", &event.blc2_v);
  tree->Branch("blc2_chi2", &event.blc2_chi2);

  tree->Branch("d5_blc1_trk", &event.d5_blc1_trk);
  tree->Branch("d5_blc2_trk", &event.d5_blc2_trk);
  tree->Branch("d5_delta", &event.d5_delta);
  tree->Branch("d5_momentum", &event.d5_momentum);
  tree->Branch("d5_fit_chi2", &event.d5_fit_chi2);
  tree->Branch("d5_fit_ndf", &event.d5_fit_ndf);
  tree->Branch("d5_fit_chi2_ndf", &event.d5_fit_chi2_ndf);
  tree->Branch("d5_pair_rank", &event.d5_pair_rank);
  tree->Branch("d5_residual_x", &event.d5_residual_x);
  tree->Branch("d5_residual_y", &event.d5_residual_y);
  tree->Branch("d5_blc2out_x", &event.d5_blc2out_x);
  tree->Branch("d5_blc2out_y", &event.d5_blc2out_y);
  tree->Branch("d5_mtxout_x", &event.d5_mtxout_x);
  tree->Branch("d5_mtxout_u", &event.d5_mtxout_u);
  tree->Branch("d5_mtxout_y", &event.d5_mtxout_y);
  tree->Branch("d5_mtxout_v", &event.d5_mtxout_v);
  tree->Branch("d5_fit_x0", &event.d5_fit_x0);
  tree->Branch("d5_fit_u0", &event.d5_fit_u0);
  tree->Branch("d5_fit_y0", &event.d5_fit_y0);
  tree->Branch("d5_fit_v0", &event.d5_fit_v0);

  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCTdcCalibMan>("DCTDC")) &&
    (InitializeParameter<DCDriftParamMan>("DCDRFT")) &&
    (InitializeParameter<DCGeomMan>("DCGEO")) &&
    (InitializeParameter<D5TransMatrix>("D5MTX")) &&
    (InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
