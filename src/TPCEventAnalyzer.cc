// -*- C++ -*-

#include "TPCEventAnalyzer.hh"

#include <sstream>

#include <TPDGCode.h>

#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "ThreeVector.hh"
#include "TPCAnalyzer.hh"
#include "TPCCluster.hh"
#include "TPCHit.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCLTrackHit.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCRawData.hh"
#include "TPCRawHit.hh"
#include "TPCReconstructor.hh"
#include "TPCVertex.hh"
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

} // namespace

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
  (void)MinTimeBucket;
  (void)MaxTimeBucket;
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

      //Check Threshold
      double min_max = HG2Poly("TPC_Raw_Max_Poly",padid+1);
      if(min_max > max_adc && max_adc < 2000){
        HF2Poly("TPC_Raw_Max_Poly",padid+1,max_adc);
        std::cout<<max_adc<<std::endl;
      }
      double fmin_adc = HG2Poly("TPC_Raw_ADC_Poly",padid+1);
      if(fmin_adc > (max_adc - mean) && (max_adc - mean) < 2000)
        HF2Poly("TPC_Raw_ADC_Poly",padid+1,max_adc - mean);
      double min_mean = HG2Poly("TPC_Raw_Mean_Poly",padid+1);
      if(min_mean > mean && mean < 2000)
      	HF2Poly("TPC_Raw_Mean_Poly",padid+1,mean);
      
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
      	if(tpc::Noise(padid)) HF2("TPC_FADC_Frame",tb, fadc.at(tb));
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
      auto row     = rhit->RowId();
      auto padid   = tpc::GetPadId(layer,row);

      HF1("TPC_FADC_Cor_Mean", mean);
      HF1("TPC_FADC_Cor_Max", max_adc);
      HF1("TPC_FADC_Cor_RMS", rms);
      HF1("TPC_FADC_Cor_LocMax", loc_max);
      HF1("TPC_FADC_Cor_Min", min_adc);
      HF1("TPC_FADC_Baseline_p0", pars.at(0));
      HF1("TPC_FADC_Baseline_p1", pars.at(1));
      HF1("TPC_FADC_Baseline_p2", pars.at(2));

      //Check Threshold
      double min_max = HG2Poly("TPC_Cor_Max_Poly",padid+1); 
      if(min_max > max_adc && max_adc < 2000)
	HF2Poly("TPC_Cor_Max_Poly",padid+1,max_adc);
      double fmin_adc = HG2Poly("TPC_Cor_ADC_Poly",padid+1);
      if(fmin_adc > (max_adc - mean) && (max_adc - mean) <2000)
	HF2Poly("TPC_Cor_ADC_Poly",padid+1,max_adc - mean);
      double min_mean = HG2Poly("TPC_Cor_Mean_Poly",padid+1);
      if(min_mean > mean && mean < 2000)
	HF2Poly("TPC_Cor_Mean_Poly",padid+1,mean);

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
std::string
TPCEventAnalyzer::DecodePidCandidates(Int_t pid_code)
{
  if (pid_code == 0)
    return "e";

  std::string names;
  if (pid_code & 0x1) names += "pi";
  if (pid_code & 0x2) {
    if (!names.empty()) names += "+";
    names += "K";
  }
  if (pid_code & 0x4) {
    if (!names.empty()) names += "+";
    names += "p";
  }
  if (names.empty())
    names = "none";
  return names;
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

  FillTrkParamHist(nhits, chisqr, x0, y0, u0, v0);

  HF1("TPCTrk_Num_Iter", track->GetNIteration());
  HF1("TPCTrk_Fitting_Flag", track->GetFitFlag());
  HF1("TPCTrk_Searching_Time", track->GetSearchTime());
  HF1("TPCTrk_Fitting_Time", track->GetFitTime());
  HF1("TPCTrk_Minuit_Status", track->GetMinuitStatus());
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillTrkParamHist(
  Int_t nhits, Double_t chisqr,
  Double_t x0, Double_t y0, Double_t u0, Double_t v0)
{
  HF1("TPCTrk_Num_TrackHits", nhits);
  HF1("TPCTrk_Chisqr", chisqr);
  HF1("TPCTrk_X0", x0);
  HF1("TPCTrk_Y0", y0);
  HF1("TPCTrk_U0", u0);
  HF1("TPCTrk_V0", v0);
  HF2("TPCTrk_U0_vs_X0", x0, u0);
  HF2("TPCTrk_V0_vs_Y0", y0, v0);
  HF2("TPCTrk_Y0_vs_X0", x0, y0);
  HF2("TPCTrk_atanU0_vs_X0", x0, TMath::ATan(u0) * TMath::RadToDeg());
  HF2("TPCTrk_atanV0_vs_Y0", y0, TMath::ATan(v0) * TMath::RadToDeg());
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
    HF2("Transverse_Diffusion", trans_dist, ratio);
    HF2(Form("Transverse_Diffusion_Layer%02d", layer), trans_dist, ratio);
  }

  if (center_hit->GetCTimeSize() > 0) {
    const TVector3 local_cal = hit->GetLocalCalPos();
    const ThreeVector local_trk(local_cal.X(), local_cal.Y(), local_cal.Z());
    const Double_t y_trk_global = gGeom.Local2GlobalPos("HypTPC", local_trk).y();
    FillCoBoClockTime("TPCTrk", layer, center_row, center_hit->GetCTime(0), center_hit->GetPosition(),
                      y_trk_global);
  }
}

void
TPCEventAnalyzer::FillBcOutTrackHist(
  Double_t chisqr_bcout,
  Double_t x0_bcout, Double_t y0_bcout,
  Double_t u0_bcout, Double_t v0_bcout,
  Double_t xtgt_bcout, Double_t ytgt_bcout,
  Double_t utgt_bcout, Double_t vtgt_bcout)
{
  HF1("BcOut_Chisqr", chisqr_bcout);
  HF1("BcOut_X0", x0_bcout);
  HF1("BcOut_Y0", y0_bcout);
  HF1("BcOut_U0", u0_bcout);
  HF1("BcOut_V0", v0_bcout);
  HF1("BcOut_XTgt", xtgt_bcout);
  HF1("BcOut_YTgt", ytgt_bcout);
  HF1("BcOut_UTgt", utgt_bcout);
  HF1("BcOut_VTgt", vtgt_bcout);
  HF2("BcOut_UTgt_vs_XTgt", xtgt_bcout, utgt_bcout);
  HF2("BcOut_VTgt_vs_YTgt", ytgt_bcout, vtgt_bcout);
  HF2("BcOut_YTgt_vs_XTgt", xtgt_bcout, ytgt_bcout);
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillTPCBcOutTgtResidualHist(
  Double_t tpc_xtgt, Double_t tpc_ytgt, Double_t tpc_utgt, Double_t tpc_vtgt,
  Double_t bcout_xtgt, Double_t bcout_ytgt, Double_t bcout_utgt, Double_t bcout_vtgt)
{
  const Double_t xtgt_diff = bcout_xtgt - tpc_xtgt;
  const Double_t ytgt_diff = bcout_ytgt - tpc_ytgt;
  const Double_t utgt_diff = bcout_utgt - tpc_utgt;
  const Double_t vtgt_diff = bcout_vtgt - tpc_vtgt;

  HF2("BcOut_vs_TPC_XTgt", tpc_xtgt, bcout_xtgt);
  HF2("BcOut_vs_TPC_YTgt", tpc_ytgt, bcout_ytgt);
  HF2("BcOut_vs_TPC_UTgt", tpc_utgt, bcout_utgt);
  HF2("BcOut_vs_TPC_VTgt", tpc_vtgt, bcout_vtgt);
  HF1("TPCTrk_ResX_Tgt", xtgt_diff);
  HF1("TPCTrk_ResY_Tgt", ytgt_diff);
  HF1("TPCTrk_ResU_Tgt", utgt_diff);
  HF1("TPCTrk_ResV_Tgt", vtgt_diff);
  HF2("TPCTrk_ResX_Tgt_vs_XTgt", bcout_xtgt, xtgt_diff);
  HF2("TPCTrk_ResY_Tgt_vs_YTgt", bcout_ytgt, ytgt_diff);
  HF2("TPCTrk_ResU_Tgt_vs_UTgt", bcout_utgt, utgt_diff);
  HF2("TPCTrk_ResV_Tgt_vs_VTgt", bcout_vtgt, vtgt_diff);
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillResidualHist(const TString& prefix,
                                   Int_t layer, Int_t row, Int_t pad,
                                   Double_t res_x, Double_t res_y,
                                   Double_t ref_x, Double_t ref_y,
                                   Bool_t in_window,
                                   Double_t ctime, const TVector3& local_pos)
{
  HF1(Form("%s_ResX", prefix.Data()), res_x);
  HF1(Form("%s_ResX_Layer%02d", prefix.Data(), layer), res_x);
  HF2(Form("%s_ResX_vs_Layer", prefix.Data()), layer, res_x);
  HF2(Form("%s_ResX_vs_X_Layer%02d", prefix.Data(), layer), ref_x, res_x);

  HF1(Form("%s_ResY", prefix.Data()), res_y);
  HF1(Form("%s_ResY_Layer%02d", prefix.Data(), layer), res_y);
  HF2(Form("%s_ResY_vs_Layer", prefix.Data()), layer, res_y);
  HF2(Form("%s_ResY_vs_Y_Layer%02d", prefix.Data(), layer), ref_y, res_y);
  if (GetDstCalibFlag()) {
    HF1(Form("%s_ResY_Layer%02d_Row%03d", prefix.Data(), layer, row), res_y);
    HF2(Form("%s_ResY_vs_Y_Layer%02d_Row%03d", prefix.Data(), layer, row), ref_y, res_y);
  }

  if (in_window) {
    const Double_t c = HG2Poly(Form("%s_HitPat", prefix.Data()), pad + 1);
    HF2Poly(Form("%s_HitPat", prefix.Data()), pad + 1, c + 1.);
    HF2(Form("%s_Row_vs_Layer", prefix.Data()), layer, row);
  }

  FillCoBoClockTime(prefix, layer, row, ctime, local_pos, ref_y);
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillTPCBcOutTrackingResidualPullHist(
  Int_t layer, Int_t center_row, Bool_t valid_tpc_resolution,
  const TVector3& trk_res_global,
  const TVector3& trk_res_local,
  const TVector3& trk_pull_global,
  const TVector3& trk_pull_local,
  Double_t cl_ref_x, Double_t cl_ref_y_tpc, Double_t cl_ref_y_bcout,
  const TVector3& cl_res)
{
  HF2("TPCTrk_ResY_vs_Layer", layer, trk_res_global.y());
  if (valid_tpc_resolution) {
    HF1(Form("TPCTrk_ResX_Layer%02d", layer), trk_res_global.x());
    HF1(Form("TPCTrk_ResY_Layer%02d", layer), trk_res_global.y());
    HF1(Form("TPCTrk_ResZ_Layer%02d", layer), trk_res_global.z());
    HF1(Form("TPCTrk_ResLocalX_Layer%02d", layer), trk_res_local.x());
    HF1(Form("TPCTrk_ResLocalY_Layer%02d", layer), trk_res_local.y());
    HF1(Form("TPCTrk_ResXY_Layer%02d", layer), std::hypot(trk_res_global.x(), trk_res_global.z()));

    HF1(Form("TPC_PullX_Layer%02d", layer), trk_pull_global.x());
    HF1(Form("TPC_PullY_Layer%02d", layer), trk_pull_global.y());
    HF1(Form("TPC_PullZ_Layer%02d", layer), trk_pull_global.z());
    HF1(Form("TPC_PullLocalX_Layer%02d", layer), trk_pull_local.x());
    HF1(Form("TPC_PullLocalY_Layer%02d", layer), trk_pull_local.y());
  }

  HF1(Form("TPCCl_ResX_Layer%02d", layer), cl_res.x());
  HF1(Form("TPCCl_ResY_Layer%02d", layer), cl_res.y());
  HF2(Form("TPCCl_ResX_vs_X_Layer%02d", layer), cl_ref_x, cl_res.x());
  HF2(Form("TPCCl_ResY_vs_Y_TPC_Layer%02d", layer), cl_ref_y_tpc, cl_res.y());
  HF2(Form("TPCCl_ResY_vs_Y_BcOut_Layer%02d", layer), cl_ref_y_bcout, cl_res.y());
  HF2("TPCCl_ResX_vs_Layer", layer, cl_res.x());
  HF2("TPCCl_ResY_vs_Layer", layer, cl_res.y());
  if (GetDstCalibFlag() && center_row >= 0) {
    HF2(Form("TPCCl_ResY_vs_Y_TPC_Layer%02d_Row%03d", layer, center_row), cl_ref_y_tpc, cl_res.y());
    HF2(Form("TPCCl_ResY_vs_Y_BcOut_Layer%02d_Row%03d", layer, center_row), cl_ref_y_bcout, cl_res.y());
  }
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillHelixHitHist(TPCLTrackHit* hit, Bool_t fill_cluster_detail, Int_t track_pid)
{
  if (!hit)
    return;

  HF1("HoughDist", hit->GetHoughDist());
  HF1("HoughDistY", hit->GetHoughDistY());

  const Int_t layer = hit->GetLayer();
  HF1("TPCTrk_Layer", layer);

  if (!fill_cluster_detail)
    return;

  TPCHit* cl_hit = hit->GetHit();
  if (!cl_hit)
    return;
  
  TPCCluster* cl = cl_hit->GetParentCluster();
  if (!cl)
    return;

  const Int_t cl_size = cl->GetClusterSize();
  const Double_t cl_de = cl->GetDe();
  const TVector3& hit_pos = hit->GetLocalHitPos();
  const TVector3& cal_pos = hit->GetLocalCalPosHelix();
  const TVector3& res_vec = hit->GetResidualVect();
  TPCHit* center_hit = cl->GetCenterHit();
  const Int_t center_row = center_hit ? center_hit->GetRow() : -1;

  if (center_row >= 0) {
    HF2("TPCTrk_Row_vs_Layer", layer, center_row);
    HF1(Form("TPCHit_HitPat_Layer%02d", layer), center_row);
    const Int_t pad_id = tpc::GetPadId(layer, center_row);
    if (pad_id >= 0) {
      const Double_t bin_cont = HG2Poly("TPCTrk_HitPat", pad_id + 1);
      HF2Poly("TPCTrk_HitPat", pad_id + 1, bin_cont + 1.);
    }
  }

  HF1("TPCCl_Size", cl_size);
  HF1(Form("TPCCl_Size_Layer%02d", layer), cl_size);
  HF1("TPCCl_dE", cl_de);
  HF1(Form("TPCCl_dE_Layer%02d", layer), cl_de);
  HF2("TPCCl_dE_vs_Layer", layer, cl_de);
  if (track_pid & 0x1) {
    HF1("TPCCl_dE_Pion", cl_de);
    HF1(Form("TPCCl_dE_Pion_Layer%02d", layer), cl_de);
  }
  HF2(Form("TPCTrk_ResY_vs_Y_Layer%02d", layer), cal_pos.y(), res_vec.y());
  HF2("TPCTrk_ResY_vs_Layer_Trk", layer, res_vec.y());

  if (GetDstCalibFlag() && center_row >= 0) {
    HF1(Form("TPCTrk_ResY_Layer%02d_Row%03d", layer, center_row), res_vec.y());
    HF2(Form("TPCTrk_ResY_vs_Y_Layer%02d_Row%03d", layer, center_row), cal_pos.y(), res_vec.y());
    HF1(Form("TPCCl_dE_Layer%02d_Row%03d", layer, center_row), cl_de);
    if (track_pid & 0x1) {
      HF1(Form("TPCCl_dE_Pion_Layer%02d_Row%03d", layer, center_row), cl_de);
    }
  }

  const TPCHitContainer& hit_cont = cl->GetHitContainer();
  for (const auto& hits : hit_cont) {
    if (!hits || !hits->IsGood() || !hit->IsGoodForTracking())
      continue;
    const TVector3& pos = hits->GetPosition();
    const Double_t pad_de = hits->GetCDe();
    Double_t trans_dist = std::hypot(hit_pos.x() - pos.x(), hit_pos.z() - pos.z());
    if (hit_pos.x() < pos.x()) trans_dist = -1.*trans_dist;
    const Double_t ratio = pad_de/cl_de;
    HF2("Transverse_Diffusion", trans_dist, ratio);
    HF2(Form("Transverse_Diffusion_Layer%02d", layer), trans_dist, ratio);
  }

  if (GetDstCalibFlag() && center_hit) {
    const ThreeVector global_ref = gGeom.Local2GlobalPos("HypTPC", cal_pos);
    FillCoBoClockTime("TPCTrk", layer, center_row, center_hit->GetCTime(0),
                      center_hit->GetPosition(), global_ref.y());
  }
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillHelixPidHist(const TPCLocalTrackHelix* track,
                                   Int_t pid_code,
                                   Double_t dedx)
{
  if(!track)
    return;

  const Double_t mom0     = track->GetMom0().Mag();
  const Int_t    charge   = track->GetCharge();
  const Double_t signed_p = static_cast<Double_t>(charge) * mom0;

  const Int_t is_beam       = track->GetIsBeam();
  const Int_t is_accidental = track->GetIsAccidental();
  const Int_t is_k18        = track->GetIsK18();
  const Bool_t is_scatter  = (is_beam == 0 && is_accidental == 0 && is_k18 == 0);

  if(is_scatter){
    HF1("Mom0", mom0);
    HF1("dEdx_PID", pid_code);
    HF2("PID_dEdx_vs_Mom", mom0, dedx);
    HF2("PID_dEdx_vs_SignedMom", signed_p, dedx);
    if(charge > 0) HF2("PID_dEdx_vs_Mom_pos", mom0, dedx);
    else HF2("PID_dEdx_vs_Mom_neg", mom0, dedx);
    if(pid_code & 0x1) HF2("PID_dEdx_vs_Mom_Pi", mom0, dedx);
    if(pid_code & 0x2) HF2("PID_dEdx_vs_Mom_K", mom0, dedx);
    if(pid_code & 0x4) HF2("PID_dEdx_vs_Mom_Proton", mom0, dedx);
  }
  if(is_beam == 1){
    HF1("Mom0_Beam", mom0);
    HF1("dEdx_PID_Beam", pid_code);
    HF2("PID_dEdx_vs_SignedMom_Beam", signed_p, dedx);
  }
  if(is_accidental == 1){
    HF1("Mom0_Accidental", mom0);
    HF1("dEdx_PID_Accidental", pid_code);
    HF2("PID_dEdx_vs_SignedMom_Accidental", signed_p, dedx);
  }
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillHelixLambdaMassHist(const TPCVertex* vertex)
{
  if (!vertex)
    return;
  const Int_t n_cand = vertex->GetNRecoCandidates();
  for (Int_t ic = 0; ic < n_cand; ++ic) {
    const TPCRecoCandidate& cand = vertex->GetRecoCandidate(ic);
    if (cand.GetMotherPdg() != kLambda0)
      continue;
    const Double_t lambda_mass_value = cand.GetMass();
    const Double_t closest_dist = cand.GetClosestDist();
    const ThreeVector lambda_vtx(cand.GetVertex().X(), cand.GetVertex().Y(), cand.GetVertex().Z());
    const ThreeVector lambda_mom(cand.GetMomentum().X(), cand.GetMomentum().Y(), cand.GetMomentum().Z());
    const ThreeVector target_to_vtx = lambda_vtx - ThreeVector(0., 0., tpc::Z_TARGET);
    const Double_t target_to_vtx_dot_mom =
      (target_to_vtx.Mag() > 0.0 && lambda_mom.Mag() > 0.0)
        ? target_to_vtx.Dot(lambda_mom)/(target_to_vtx.Mag()*lambda_mom.Mag())
        : TMath::QuietNaN();
    HF1("Lambda_Mass", lambda_mass_value);

    HF1("Lambda_CloseDist", closest_dist);
    HF1("Lambda_VtxX", lambda_vtx.x());
    HF1("Lambda_VtxY", lambda_vtx.y());
    HF1("Lambda_VtxZ", lambda_vtx.z());
    HF1("Lambda_MomX", lambda_mom.x());
    HF1("Lambda_MomY", lambda_mom.y());
    HF1("Lambda_MomZ", lambda_mom.z());
    HF1("Lambda_TargetToVtxX", target_to_vtx.x());
    HF1("Lambda_TargetToVtxY", target_to_vtx.y());
    HF1("Lambda_TargetToVtxZ", target_to_vtx.z());
    HF1("Lambda_TargetToVtxDotMom", target_to_vtx_dot_mom);
  }
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillHelixK0ShortMassHist(const TPCVertex* vertex)
{
  if (!vertex)
    return;
  const Int_t n_cand = vertex->GetNRecoCandidates();
  for (Int_t ic = 0; ic < n_cand; ++ic) {
    const TPCRecoCandidate& cand = vertex->GetRecoCandidate(ic);
    if (cand.GetMotherPdg() != kK0Short)
      continue;
    const Double_t k0_mass_value = cand.GetMass();
    const Double_t closest_dist = cand.GetClosestDist();
    const ThreeVector k0_vtx(cand.GetVertex().X(), cand.GetVertex().Y(), cand.GetVertex().Z());
    const ThreeVector k0_mom(cand.GetMomentum().X(), cand.GetMomentum().Y(), cand.GetMomentum().Z());
    const ThreeVector target_to_vtx = k0_vtx - ThreeVector(0., 0., tpc::Z_TARGET);
    const Double_t target_to_vtx_dot_mom =
      (target_to_vtx.Mag() > 0.0 && k0_mom.Mag() > 0.0)
        ? target_to_vtx.Dot(k0_mom)/(target_to_vtx.Mag()*k0_mom.Mag())
        : TMath::QuietNaN();
    HF1("K0_Mass", k0_mass_value);

    HF1("K0_CloseDist", closest_dist);
    HF1("K0_VtxX", k0_vtx.x());
    HF1("K0_VtxY", k0_vtx.y());
    HF1("K0_VtxZ", k0_vtx.z());
    HF1("K0_MomX", k0_mom.x());
    HF1("K0_MomY", k0_mom.y());
    HF1("K0_MomZ", k0_mom.z());
    HF1("K0_TargetToVtxX", target_to_vtx.x());
    HF1("K0_TargetToVtxY", target_to_vtx.y());
    HF1("K0_TargetToVtxZ", target_to_vtx.z());
    HF1("K0_TargetToVtxDotMom", target_to_vtx_dot_mom);
  }
}

//_____________________________________________________________________________
void
TPCEventAnalyzer::FillCoBoClockTime(const TString& prefix,
                                    Int_t layer, Int_t row,
                                    Double_t ctime,
                                    const TVector3& local_pos,
                                    Double_t ref_y)
{
  if (m_clkTpc.size() != static_cast<size_t>(NumOfSegCOBO)) {
    spdlog::warn("TPCEventAnalyzer::FillCoBoClockTime[{}]: clkTpc size ({}) != NumOfSegCOBO ({})",
                 prefix.Data(), m_clkTpc.size(), NumOfSegCOBO);
    return;
  }

  const Int_t cobo = tpc::GetCoBoId(layer, row);
  const Int_t asad = tpc::GetASADId(layer, row);
  const Bool_t cobo_valid = (0 <= cobo && cobo < NumOfSegCOBO);

  const ThreeVector global_pos = gGeom.Local2GlobalPos("HypTPC", local_pos);
  const Double_t res_y = global_pos.y() - ref_y;

  Double_t res_y_raw   = TMath::QuietNaN();
#ifdef DEBUG_COBO_CLOCK
  Double_t res_y_noclk = TMath::QuietNaN();
#endif

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
        const ThreeVector local_raw  (local_pos.x(), y_raw,   local_pos.z());
        const ThreeVector global_raw = gGeom.Local2GlobalPos("HypTPC", local_raw);
        res_y_raw = global_raw.y() - ref_y;
#ifdef DEBUG_COBO_CLOCK
        const ThreeVector local_noclk(local_pos.x(), y_noclk, local_pos.z());
        const ThreeVector global_noclk = gGeom.Local2GlobalPos("HypTPC", local_noclk);
        res_y_noclk = global_noclk.y() - ref_y;
#endif
      }
    }

    HF2(Form("%s_ResY_vs_ClockTime_CoBo%d", prefix.Data(), cobo),
        m_clkTpc.at(cobo), res_y);
    if (std::isfinite(res_y_raw)) {
      HF2(Form("%s_ResY_vs_ClockTime_CoBo%d_RawClock", prefix.Data(), cobo),
          m_clkTpc.at(cobo), res_y_raw);
    }
#ifdef DEBUG_COBO_CLOCK
    if (std::isfinite(res_y_noclk)) {
      HF2(Form("%s_ResY_vs_ClockTime_CoBo%d_NoClock", prefix.Data(), cobo),
          m_clkTpc.at(cobo), res_y_noclk);
    }
#endif
  }

  const Bool_t asad_valid = (0 <= asad && asad < NumOfAsadTPC);
  if (asad_valid && cobo_valid && std::isfinite(m_clkTpc.at(cobo))) {
    const Double_t clk = m_clkTpc.at(cobo);
    HF2(Form("%s_ResY_vs_ClockTime_Asad%02d", prefix.Data(), asad),
        clk, res_y);
    if (std::isfinite(res_y_raw)) {
      HF2(Form("%s_ResY_vs_ClockTime_Asad%02d_RawClock", prefix.Data(), asad),
          clk, res_y_raw);
    }
#ifdef DEBUG_COBO_CLOCK
    if (std::isfinite(res_y_noclk)) {
      HF2(Form("%s_ResY_vs_ClockTime_Asad%02d_NoClock", prefix.Data(), asad),
          clk, res_y_noclk);
    }
#endif
  }
}
