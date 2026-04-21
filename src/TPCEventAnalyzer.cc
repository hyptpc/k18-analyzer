// -*- C++ -*-

#include "TPCEventAnalyzer.hh"

#include <cmath>
#include <sstream>

#include "DCGeomMan.hh"
#include "RootHelper.hh"
#include "ThreeVector.hh"
#include "TPCCluster.hh"
#include "TPCHit.hh"
#include "TPCLTrackHit.hh"
#include "TPCLocalTrack.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCRawData.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gTpcParam = TPCParamMan::GetInstance();
using root::HF1;
using root::HF2;
using root::HF2Poly;
using root::HG2Poly;
}

//_____________________________________________________________________________
TPCEventAnalyzer::TPCEventAnalyzer()
{
}

//_____________________________________________________________________________
TPCEventAnalyzer::~TPCEventAnalyzer()
{
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::TPCRawHit(const TPCRawData& TPCrawData)
{
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");
  static const Int_t MinTimeBucket = gUser.GetParameter("TimeBucketTPC", 0);
  static const Int_t MaxTimeBucket = gUser.GetParameter("TimeBucketTPC", 1);
  Int_t npadTpc_raw = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = TPCrawData.GetTPCRawHits(layer);
    const auto nhit = hc.size();
    npadTpc_raw += nhit;
    for(const auto& rhit : hc){
      auto mean    = rhit->Mean(0, NumOfTimeBucket);
      auto max_adc = rhit->MaxAdc(0, NumOfTimeBucket);
      auto min_adc = rhit->MinAdc(0, NumOfTimeBucket);
      auto rms     = rhit->RMS(0, NumOfTimeBucket);
      auto loc_max = rhit->LocMax(0, NumOfTimeBucket);
      auto row     = rhit->RowId();
      auto padid   = tpc::GetPadId(layer,row);

      HF1("TPC_FADC_Mean", mean);
      HF1("TPC_FADC_Max", max_adc);
      HF1("TPC_FADC_RMS", rms);
      HF1("TPC_FADC_LocMax", loc_max);
      HF1("TPC_FADC_Min", min_adc);
      
      auto gate_open_max_adc = rhit->MaxAdc(0,50);
      auto gate_close_max_adc = rhit->MaxAdc(50,140);
      auto middle_rms = rhit->RMS(50, 140);
      auto open_rms = rhit->RMS(0,50);
      auto nmax_adc = rhit->MaxAdc(50, 140);
      auto nmin_adc = rhit->MinAdc(50, 140);

      Bool_t IsNoise = false;
      
      if(gate_open_max_adc > 600 && gate_open_max_adc < 800 && middle_rms <30 && middle_rms >10 && open_rms > 35 && open_rms < 60 && gate_open_max_adc > nmax_adc && gate_close_max_adc < 800){
        IsNoise = true;
        HF1("TPC_FADC_Noise_Max", gate_open_max_adc);
        HF1("TPC_FADC_Noise_RMSfront", open_rms);
        HF1("TPC_FADC_Noise_RMSmiddle", middle_rms);
        HF1("TPC_FADC_Noise_Adcdiff", nmax_adc - nmin_adc);
      }
      auto fadc = rhit->Fadc();
      for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
        HF2("TPC_FADC_Before", tb, fadc.at(tb));
        if(IsNoise){
          HF2("TPC_FADC_Noise",tb,fadc.at(tb));
        }
	if(tpc::Noise(padid))HF2("TPC_FADC_Frame",tb, fadc.at(tb));
      }
      if(IsNoise){
        Double_t bincont = HG2Poly("TPC_HitPat_Noise",padid+1);
        HF2Poly("TPC_HitPat_Noise",padid+1,bincont+1.);
      }
    } 
  }
  HF1("TPC_Multiplicity_Raw", npadTpc_raw);
}

//_____________________________________________________________________________
TPCBaselineInfo
TPCEventAnalyzer::TPCBaselineHit(const TPCRawData &TPCrawData){
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");
  TPCBaselineInfo baseout;
  
  auto baseline = TPCrawData.GetBaselineTPC();
  if(baseline){
    auto brow     = baseline->RowId();
    auto blayer   = baseline->LayerId();
    auto bmean    = baseline->Mean(0, NumOfTimeBucket);
    auto brms     = baseline->RMS(0, NumOfTimeBucket);
    auto bmax_adc = baseline->MaxAdc(0, NumOfTimeBucket);
    auto bmin_adc = baseline->MinAdc(0, NumOfTimeBucket);
    auto bloc_max = baseline->LocMax(0, NumOfTimeBucket);
    
    baseout.row = brow;
    baseout.layer = blayer;
    baseout.rms = brms;
    baseout.valid = kTRUE;

    HF1("TPC_FADC_Baseline_Mean", bmean);
    HF1("TPC_FADC_Baseline_Max", bmax_adc);
    HF1("TPC_FADC_Baseline_RMS", brms);
    HF1("TPC_FADC_Baseline_LocMax", bloc_max);
    HF1("TPC_FADC_Baseline_Min", bmin_adc);
    
    auto fadc = baseline->Fadc();
    for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
      HF2("TPC_FADC_Baseline", tb, fadc.at(tb));
    }
    Int_t bpadid = tpc::GetPadId(baseline->LayerId(), baseline->RowId());
    Double_t bincont = HG2Poly("TPC_HitPat_Baseline",bpadid+1);
    HF2Poly("TPC_HitPat_Baseline",bpadid+1,bincont+1.);
  }
  
  return baseout;
}


//_____________________________________________________________________________
Int_t
TPCEventAnalyzer::TPCCorHit(const TPCRawData &TPCrawData){
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");
  Int_t npadTpc = 0;
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    auto hc = TPCrawData.GetTPCCorHits(layer);
    const auto nhit = hc.size();
    npadTpc += nhit;

    for(const auto& rhit : hc){
      auto mean    = rhit->Mean(0, NumOfTimeBucket);
      auto max_adc = rhit->MaxAdc(0, NumOfTimeBucket);
      auto min_adc = rhit->MinAdc(0, NumOfTimeBucket);
      auto rms     = rhit->RMS(0, NumOfTimeBucket);
      auto loc_max = rhit->LocMax(0, NumOfTimeBucket);
      auto pars    = rhit->GetParameters();

      HF1("TPC_FADC_Cor_Mean", mean);
      HF1("TPC_FADC_Cor_Max", max_adc);
      HF1("TPC_FADC_Cor_RMS", rms);
      HF1("TPC_FADC_Cor_LocMax", loc_max);
      HF1("TPC_FADC_Cor_Min", min_adc);
      HF1("TPC_FADC_Baseline_p0", pars.at(0));
      HF1("TPC_FADC_Baseline_p1", pars.at(1));
      HF1("TPC_FADC_Baseline_p2", pars.at(2));

      // 2D FADC waveform after correction
      auto fadc = rhit->Fadc();
      for(Int_t tb = 0, ntb = fadc.size(); tb < ntb; ++tb){
        HF2("TPC_FADC_After", tb, fadc.at(tb));
      }
    }
  }

  HF1("TPC_Multiplicity_Cor", npadTpc);
  return npadTpc;
}

//_____________________________________________________________________________
// NOTE: Disabled; TPCHit method name collides with TPCHit type name and harms readability.
// void
// TPCEventAnalyzer::TPCHit(const TPCAnalyzer& TPCAna)
// {
// }

//_____________________________________________________________________________
Bool_t
TPCEventAnalyzer::ValidateCoboClocks(const std::vector<Double_t>& clk_tpc)
{
  if (static_cast<Int_t>(clk_tpc.size()) != NumOfSegCOBO) {
    spdlog::warn("something is wrong: clkTpc.size() != {}", NumOfSegCOBO);
    return false;
  }
  std::vector<Int_t> bad_cobo;
  for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
    if (!std::isfinite(clk_tpc.at(static_cast<size_t>(cobo))))
      bad_cobo.push_back(cobo);
  }
  if (!bad_cobo.empty()) {
    std::ostringstream oss;
    for (size_t i = 0; i < bad_cobo.size(); ++i) {
      oss << (i ? "," : "") << bad_cobo[i];
    }
    spdlog::warn("CoBo clock(s) missing (NaN/Inf): cobo={}, skip event", oss.str());
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Double_t
TPCEventAnalyzer::CalcTruncatedMean(const std::vector<Double_t>& cumulative_vec, Double_t fraction)
{
  if (cumulative_vec.empty()) {
    spdlog::warn("CalcTruncatedMean: empty cumulative");
    return 0.;
  }
  const std::size_t size = cumulative_vec.size() - 1;
  const Int_t n_trunc = static_cast<Int_t>(static_cast<Double_t>(size) * fraction);
  if (n_trunc <= 0 || n_trunc > static_cast<Int_t>(size)) {
    spdlog::warn("CalcTruncatedMean: invalid fraction={}, n_trunc={}, size={}", fraction, n_trunc,
                 size);
    return 0.;
  }
  return cumulative_vec[static_cast<std::size_t>(n_trunc)] / static_cast<Double_t>(n_trunc);
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillTrkHist(const TPCLocalTrack* track)
{
  if (!track)
    return;

  const Int_t nhits = track->GetNHit();
  const Double_t chisqr = track->GetChiSquare();
  const Double_t x0 = track->GetX0();
  const Double_t y0 = track->GetY0();
  const Double_t u0 = track->GetU0();
  const Double_t v0 = track->GetV0();

  HF1("TPCTrk_Num_TrackHits", nhits);
  HF1("TPCTrk_Chisqr", chisqr);
  HF1("TPCTrk_X0", x0);
  HF1("TPCTrk_Y0", y0);
  HF1("TPCTrk_U0", u0);
  HF1("TPCTrk_V0", v0);
  HF2("TPCTrk_U0_vs_X0", x0, u0);
  HF2("TPCTrk_V0_vs_Y0", y0, v0);
  HF2("TPCTrk_Y0_vs_X0", x0, y0);

  HF1("TPCTrk_Num_Iter", track->GetNIteration());
  HF1("TPCTrk_Fitting_Flag", track->GetFitFlag());
  HF1("TPCTrk_Searching_Time", track->GetSearchTime());
  HF1("TPCTrk_Fitting_Time", track->GetFitTime());
  HF1("TPCTrk_Minuit_Status", track->GetMinuitStatus());
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillTrkHitHist(TPCLTrackHit* hit, const TPCLocalTrack* track)
{
  if (!hit || !track)
    return;

  HF1("TPCTrk_Hough_Dist", hit->GetHoughDist());
  HF1("TPCTrk_Hough_DistY", hit->GetHoughDistY());

  const Int_t layer = hit->GetLayer();
  const TVector3& hit_pos = hit->GetLocalHitPos();
  const TVector3& cal_pos = hit->GetLocalCalPos();
  const TVector3& resi_vect = hit->GetResidualVect();

  TPCHit* cl_hit = hit->GetHit();
  if (!cl_hit)
    return;
  TPCCluster* cl = cl_hit->GetParentCluster();
  if (!cl)
    return;

  const Int_t cl_size = cl->GetClusterSize();
  const Double_t clde = cl->GetDe();
  TPCHit* center_hit = cl->GetCenterHit();
  if (!center_hit)
    return;

  const Int_t center_row = center_hit->GetRow();
  const Double_t residual = hit->GetResidual();

  HF1("TPCTrk_Layer", layer);
  HF2("TPCTrk_Row_vs_Layer", layer, center_row);
  HF1(Form("TPCHit_HitPat_Layer%02d", layer), center_row);
  const Int_t pad_id = tpc::GetPadId(layer, center_row);
  if (pad_id >= 0) {
    const Double_t bin_cont = HG2Poly("TPCTrk_HitPat", pad_id + 1);
    HF2Poly("TPCTrk_HitPat", pad_id + 1, bin_cont + 1.);
  }
  HF1(Form("TPCHit_Xhit_Layer%02d", layer), hit_pos.x());
  HF1(Form("TPCTrk_Res_Layer%02d", layer), residual);
  HF2(Form("TPCTrk_Res_vs_Xhit_Layer%02d", layer), hit_pos.x(), residual);
  HF2(Form("TPCHit_Yhit_vs_Xtrk_Layer%02d", layer), cal_pos.x(), hit_pos.y());
  HF1(Form("TPCTrk_ResX_Layer%02d", layer), resi_vect.X());
  HF1(Form("TPCTrk_ResY_Layer%02d", layer), resi_vect.Y());
  HF1(Form("TPCTrk_ResZ_Layer%02d", layer), resi_vect.Z());
  HF2(Form("TPCTrk_ResY_vs_Y_Layer%02d", layer), cal_pos.y(), resi_vect.Y());
  if (GetDstCalibFlag()) {
    HF1(Form("TPCTrk_ResY_Layer%02d_Row%03d", layer, center_row), resi_vect.Y());
    HF2(Form("TPCTrk_ResY_vs_Y_Layer%02d_Row%03d", layer, center_row), cal_pos.y(), resi_vect.Y());
    HF1(Form("TPCCl_dE_Layer%02d_Row%03d", layer, center_row), clde);
  }
  HF2("TPCTrk_ResX_vs_Layer_Trk", layer, resi_vect.X());
  HF2("TPCTrk_ResY_vs_Layer_Trk", layer, resi_vect.Y());
  HF2("TPCTrk_ResZ_vs_Layer_Trk", layer, resi_vect.Z());
  HF2("TPCCl_dE_vs_Layer", layer, clde);

  HF1("TPCCl_Size", cl_size);
  HF1(Form("TPCCl_Size_Layer%02d", layer), cl_size);
  HF1("TPCCl_dE", clde);
  HF1(Form("TPCCl_dE_Layer%02d", layer), clde);
  const TPCHitContainer& hit_cont = cl->GetHitContainer();
  for (const auto& hits : hit_cont) {
    if (!hits || !hits->IsGood())
      continue;
    const TVector3& pos = hits->GetPosition();
    const Double_t de = hits->GetCDe();
    const Double_t dummy = std::hypot(pos.x() - hit_pos.x(), pos.z() - hit_pos.z());
    const Double_t trans_dist = (hit_pos.x() - pos.x() < 0.) ? -dummy : dummy;
    const Double_t ratio = de / clde;
    HF2("TPCCl_Ratio_vs_Dist_Diff", trans_dist, ratio);
    HF2(Form("TPCCl_Ratio_vs_Dist_Diff_Layer%02d", layer), trans_dist, ratio);
  }

  if (center_hit->GetCTimeSize() > 0) {
    const TVector3 local_cal = hit->GetLocalCalPos();
    const ThreeVector local_trk(local_cal.X(), local_cal.Y(), local_cal.Z());
    const Double_t y_trk_global = gGeom.Local2GlobalPos("HypTPC", local_trk).y();
    FillCoBoClockTime("TPCTrk", layer, center_row, center_hit->GetCTime(0), center_hit->GetPosition(),
                      y_trk_global);
  }
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillCoBoClockTime(const TString& prefix,
                                    Int_t layer, Int_t row,
                                    Double_t ctime,
                                    const TVector3& localPos,
                                    Double_t referenceY)
{
  if (m_clkTpc.size() != static_cast<size_t>(NumOfSegCOBO)) {
    spdlog::warn("TPCEventAnalyzer::FillCoBoClockTime[{}]: clkTpc size ({}) != NumOfSegCOBO ({})",
                 prefix.Data(), m_clkTpc.size(), NumOfSegCOBO);
    return;
  }

  const Int_t cobo = tpc::GetCoBoId(layer, row);
  const Int_t asad = tpc::GetASADId(layer, row);
  const Bool_t cobo_valid = (0 <= cobo && cobo < NumOfSegCOBO);

  const ThreeVector globalPos = gGeom.Local2GlobalPos("HypTPC", localPos);
  const Double_t resY = globalPos.y() - referenceY;

  Double_t resY_raw   = TMath::QuietNaN();
  Double_t resY_noclk = TMath::QuietNaN();

  if (!cobo_valid) {
    spdlog::warn("TPCEventAnalyzer::FillCoBoClockTime[{}]: invalid CoBo id (cobo={}) for layer={} row={}",
                 prefix.Data(), cobo, layer, row);
  } else if (!std::isfinite(m_clkTpc.at(cobo))) {
    spdlog::warn("TPCEventAnalyzer::FillCoBoClockTime[{}]: non-finite clkTpc[{}]={} for layer={} row={}",
                 prefix.Data(), cobo, m_clkTpc.at(cobo), layer, row);
  } else {
    const Double_t clk = m_clkTpc.at(cobo);
    Double_t cclk = 0.0;
    if (!gTpcParam.GetCClock(layer, row, clk, cclk)) {
      spdlog::warn("TPCEventAnalyzer::FillCoBoClockTime[{}]: GetCClock failed for layer={} row={} clk={}",
                   prefix.Data(), layer, row, clk);
    } else {
      const Double_t ctime_noclk = ctime - cclk;
      Double_t y_noclk = 0.0;
      Double_t y_raw   = 0.0;
      const Bool_t ok_noclk = gTpcParam.GetDriftLength(layer, row, ctime_noclk, y_noclk);
      const Bool_t ok_raw   = gTpcParam.GetDriftLength(layer, row, ctime_noclk + clk, y_raw);
      if (!ok_noclk || !ok_raw) {
        spdlog::warn("TPCEventAnalyzer::FillCoBoClockTime[{}]: GetDriftLength failed (noclk={}, raw={}) for layer={} row={}",
                     prefix.Data(), ok_noclk, ok_raw, layer, row);
      } else {
        const ThreeVector local_noclk(localPos.x(), y_noclk, localPos.z());
        const ThreeVector local_raw  (localPos.x(), y_raw,   localPos.z());
        const ThreeVector global_noclk = gGeom.Local2GlobalPos("HypTPC", local_noclk);
        const ThreeVector global_raw   = gGeom.Local2GlobalPos("HypTPC", local_raw);
        resY_noclk = global_noclk.y() - referenceY;
        resY_raw   = global_raw.y()   - referenceY;
      }
    }

    HF2(Form("%s_ResY_vs_ClockTime_CoBo%d", prefix.Data(), cobo),
        m_clkTpc.at(cobo), resY);
    if (std::isfinite(resY_raw)) {
      HF2(Form("%s_ResY_vs_ClockTime_CoBo%d_RawClock", prefix.Data(), cobo),
          m_clkTpc.at(cobo), resY_raw);
    }
#ifdef DEBUG_COBO_CLOCK
    if (std::isfinite(resY_noclk)) {
      HF2(Form("%s_ResY_vs_ClockTime_CoBo%d_NoClock", prefix.Data(), cobo),
          m_clkTpc.at(cobo), resY_noclk);
    }
#endif
  }

  const Bool_t asad_valid = (0 <= asad && asad < NumOfAsadTPC);
  if (asad_valid && cobo_valid && std::isfinite(m_clkTpc.at(cobo))) {
    const Double_t clk = m_clkTpc.at(cobo);
    HF2(Form("%s_ResY_vs_ClockTime_Asad%02d", prefix.Data(), asad),
        clk, resY);
    if (std::isfinite(resY_raw)) {
      HF2(Form("%s_ResY_vs_ClockTime_Asad%02d_RawClock", prefix.Data(), asad),
          clk, resY_raw);
    }
#ifdef DEBUG_COBO_CLOCK
    if (std::isfinite(resY_noclk)) {
      HF2(Form("%s_ResY_vs_ClockTime_Asad%02d_NoClock", prefix.Data(), asad),
          clk, resY_noclk);
    }
#endif
  }
}
