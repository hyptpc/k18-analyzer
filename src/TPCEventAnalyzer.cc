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
// #include "TPCParamMan.hh"
#include "TPCRawData.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
const auto& gUser = UserParamMan::GetInstance();
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

      bool IsNoise = false;
      
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
      }
      if(IsNoise){
	double bincont = HG2Poly("TPC_HitPat_Noise",padid+1);
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
    int bpadid = tpc::GetPadId(baseline->RowId(), baseline->LayerId());
    double bincont = HG2Poly("TPC_HitPat_Baseline",bpadid+1);
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
