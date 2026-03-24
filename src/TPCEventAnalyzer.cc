// -*- C++ -*-

#include "TPCEventAnalyzer.hh"

#include <DAQNode.hh>
#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#include "DetectorID.hh"
#include "RootHelper.hh"
#include "TPCAnalyzer.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include "DCGeomMan.hh"
#include "TPCParamMan.hh"
#include "TPCLTrackHit.hh"
#include "ThreeVector.hh"
#include "TPCRawData.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
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
void
TPCEventAnalyzer::TPCHit(const TPCAnalyzer& TPCAnalyzer){
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
