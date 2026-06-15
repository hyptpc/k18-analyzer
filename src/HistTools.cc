// -*- C++ -*-

#include "HistTools.hh"

#include "DetectorID.hh"
#include "RootHelper.hh"
#include "TPCPadHelper.hh"

#include <DAQNode.hh>
#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
using namespace root; 
}

namespace hist
{
// Raw TDC: {n_bin, x_min, x_max}
const Double_t hr_tdc_bins[3]       = {20000.,  600000.,  800000.};  // default
const Double_t hr_tdc_bins_ftof[3]  = {20000.,  200000.,  600000.};
const Double_t hr_tdc_bins_bht[3]   = {20000., 1200000., 1600000.};
const Double_t hr_tdc_bins_cobo[3]  = {10000., 1500000., 1650000.};
const Double_t hr_tdc_bins_trig[3]  = {10000.,       0., 2000000.};  // TriggerFlag
const Double_t hr_tot_bins[3]       = {5000., 0., 50000.};
const Double_t adc_bins[3]          = {4096., -0.5, 4095.5};
const Double_t mh_tdc_bins[3]       = {2000., 0., 2000.};
const Double_t mh_tot_bins[3]       = {1000., 0., 1000.};
const Double_t mh_tdc_bins_sch[3]   = {1000., 0., 1000.}; //SCH MHTDC
// HodoHit
const Double_t hr_time_bins[3]      = {5000., -50., 50.};
const Double_t mh_time_bins[3]      = {500., -50., 50.};
const Double_t hr_tot_time_bins[3]  = {1000., 0., 200.};
const Double_t de_bins[3]           = {1000., 0., 10.};
const Double_t npe_bins[3]          = {700., -50., 300.};  // for Cherenkov (BAC, KVC, SAC3)


// HB2Poly bounds: {xmin, xmax, ymin, ymax}
const Double_t TPC_EVENT_DISPLAY_BINS[4] = { -300., 300., -300., 300. };

// --- Track / vertex common (axes shared with TPC–external reference plots) ---
const Double_t TPC_BINS_LAYER[3] = { static_cast<Double_t>(NumOfLayersTPC),
                                     -0.5,
                                     static_cast<Double_t>(NumOfLayersTPC) - 0.5 };
const Double_t TPC_BINS_ROW_VS_LAYER[6] = {
  static_cast<Double_t>(NumOfLayersTPC), -0.5, static_cast<Double_t>(NumOfLayersTPC) - 0.5,
  244., -0.5, 243.5
};
const Double_t TPC_BINS_MULT[3]       = {40.,    0.,   40.};
const Double_t TPC_BINS_HITS[3]       = {50.,    0.,   50.};
const Double_t TPC_BINS_CHISQR[3]     = {500.,   0.,  100.};
const Double_t TPC_BINS_POS_X[3]      = {500., -250., 250.};
const Double_t TPC_BINS_POS_Y[3]      = {400., -100., 100.};
const Double_t TPC_BINS_SLOPE[3]      = {200., -0.20, 0.20};
const Double_t TPC_BINS_ATAN_DEG[3]   = {100., -20.0,  20.0};

// --- Residuals & pulls (TPC fit / hit / cluster) ---
const Double_t TPC_BINS_RES_TRACK[3]  = {200.,  -2.0,  2.0};
const Double_t TPC_BINS_RES_MED[3]    = {200.,  -8.0,  8.0};
const Double_t TPC_BINS_RES_HIT[3]    = {100., -20.0, 20.0};
const Double_t TPC_BINS_RES_ABS[3]    = {200.,   0.0, 10.0};
const Double_t TPC_BINS_PULL[3]       = {200.,  -5.0,  5.0};
const Double_t TPC_BINS_RES_TGT[3]    = {640., -16.0,  16.0};
const Double_t TPC_BINS_RES_TGT_SLOPE[3]  = {200., -0.05, 0.05};

// --- External reference track (BcOut) axes in TPC analysis context ---
const Double_t TPC_BINS_EXT_MULT[3]        = {20.,    0.,   20.};
const Double_t TPC_BINS_EXT_CHISQR[3]      = {500.,   0.,  500.};
const Double_t TPC_BINS_EXT_TGT_POS[3]     = {400., -200., 200.};
const Double_t TPC_BINS_EXT_TGT_SLOPE_U[3] = {400., -0.15, 0.15};
const Double_t TPC_BINS_EXT_TGT_SLOPE_V[3] = {400., -0.05, 0.05};

// --- Clock time (CoBo / AsAd) ---
const Double_t TPC_BINS_CLOCK[3]       = {2400., -60.0, 60.0};

// --- BuildTPCTracking: cluster & fitter diagnostics (promoted from function-local) ---
const Double_t TPC_BINS_DE[3]          = {400.,   0.,  800.};
const Double_t TPC_BINS_CL_SIZE[3]     = {25.,    0.,   25.};
const Double_t TPC_BINS_DIFF_PAD[3]    = {60.,  -15.0,  15.0};
const Double_t TPC_BINS_HOUGH[3]       = {500.,   0.,   50.};
const Double_t TPC_BINS_TIME_MS[3]     = {100.,   0.,  100.};
const Double_t TPC_BINS_ITER[3]        = {100.,   0.,  100.};
const Double_t TPC_BINS_FLAG[3]        = {10.,    0.,   10.};
const Double_t TPC_BINS_MINUIT[3]      = {5.,     0.,    5.};
const Double_t TPC_BINS_CLUSTER_RATIO[3] = {100., 0., 1.};
const Double_t TPC_BINS_MOM0[3]        = {1000.,  0.,  1.5};
const Double_t TPC_BINS_PID_CODE[3]    = {16.,   -0.5, 15.5};
const Double_t TPC_BINS_SIGNED_P[3]    = {600.,  -3.0,  3.0};

// --- BuildTPCBcOutTracking: 2D residual vs external target (promoted from function-local) ---
const Double_t TPC_BINS_RES_TGT_2D_POS[3] = {400., -150., 150.};
const Double_t TPC_BINS_RES_TGT_2D_X[3]   = {400., -10., 10.};
const Double_t TPC_BINS_RES_TGT_2D_Y[3]   = {400., -15., 15.};
const Double_t TPC_BINS_RES_TGT_2D_UV[3]  = {400., -0.03, 0.03};
const Double_t TPC_BINS_RES_LOCAL_X[3]    = {400., -16.,  16.};
const Double_t TPC_BINS_RES_XY_NORM[3]    = {100.,   0.,   8.};

//_____________________________________________________________________________
void
BuildStatus()
{
  HB1("Status", 21, -0.5, 20.5);
}

//_____________________________________________________________________________
void
BuildTriggerFlag()
{
  const Char_t* name = "TriggerFlag";
  Double_t patbins[3] = {NumOfSegTrigFlag, -0.5, NumOfSegTrigFlag-0.5};
  for(Int_t i=0; i<NumOfSegTrigFlag; ++i){
    HB1(Form("%s_TDC_seg%d", name, i), hr_tdc_bins_trig);
  }
  HB1(Form("%s_HitPat; Segment; Counts", name), patbins);
  auto h1 = HB1("BeamFlag", beam::kBeamFlag, -0.5, beam::kBeamFlag - 0.5);
  for(Int_t i=0, n=beam::BeamFlagList.size(); i<n; ++i){
    auto label = beam::BeamFlagList.at(i);
    if(label.IsNull()) label = "All";
    else label.ReplaceAll("_", "");
    h1->GetXaxis()->SetBinLabel(i+1, label);
  }
  h1->GetXaxis()->SetBinLabel(beam::kBeamFlag, "Unknown");
}

//_____________________________________________________________________________
void
BuildHodoRaw(Bool_t flag_beam_particle)
{
  for(const auto& beam: beam::BeamFlagList){
    const Char_t* b = beam.Data();
    { // BHT
      const Char_t* name = "BHT";
      Double_t nseg = static_cast<Double_t>(NumOfSegBHT);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"} ){
          const Char_t* ud = uord.Data();
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b),      hr_tdc_bins_bht);
          HB1(Form("%s_Trailing_seg%d%s%s; channel; count", name, i, ud, b), hr_tdc_bins_bht);
          HB1(Form("%s_TOT_seg%d%s%s; channel; count", name, i, ud, b),      hr_tot_bins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s; segment; count", name, ud, b),     seg_bins);
        HB1(Form("%s_Multi_%s%s; multiplicity; count", name, ud, b), mul_bins);
      }
    }

    { // BAC
      const Char_t* name = "BAC";
      Double_t nseg = static_cast<Double_t>(NumOfSegBAC);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_ADC_seg%d%s; channel; count", name, i, b),  adc_bins);
        HB1(Form("%s_AwT_seg%d%s; channel; count", name, i, b),  adc_bins);
        HB1(Form("%s_AwoT_seg%d%s; channel; count", name, i, b), adc_bins);
        HB1(Form("%s_TDC_seg%d%s; channel; count", name, i, b),  hr_tdc_bins);
      }
      HB1(Form("%s_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%s_Multi%s; multiplicity; count", name, b), mul_bins);
    }
    
    // BHT-BAC
    {
      HB2(Form("BAC_ADC_vs_BHT_TDC%s", b),
          200, 720000., 750000., 200, 0., 2000.);
    }

    // Hodoscope (U/D or 1ch). 
    // KVC: a,b,c,d,S only. COBO: TDC(U) only. SCH: MHTDC only → dedicated blocks below.
    for(Int_t ihodo=kBH2; ihodo<kNumHodo; ++ihodo){
      if(ihodo == kKVC || ihodo == kCOBO || ihodo == kSCH) continue;
      auto name = NameHodo[ihodo].Data();
      const Double_t* tdc_bins;
      if (HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Ftof))
        tdc_bins = hr_tdc_bins_ftof;
      else
        tdc_bins = hr_tdc_bins;
      Double_t nseg = static_cast<Double_t>(NumOfSegHodo[ihodo]);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      Bool_t is_one_side = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::OneSideReadout);
      for (const auto& uord : std::vector<TString>{"U", "D"}) {
        if (is_one_side && uord == "D") continue;
        auto ud = uord.Data();
        for (Int_t i = 0; i < nseg; ++i) {
          HB1(Form("%s_ADC_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
          HB1(Form("%s_AwT_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
          HB1(Form("%s_AwoT_seg%d%s%s; channel; count", name, i, ud, b), adc_bins);
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b),  tdc_bins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        if (is_one_side && uord == "AND") continue;
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s; segment; count", name, ud, b),     seg_bins);
        HB1(Form("%s_Multi_%s%s; multiplicity; count", name, ud, b), mul_bins);
      }
    }

    { // HTOF Sum (S): U/D is in the loop above
      auto name = "HTOF";
      Double_t nseg = static_cast<Double_t>(NumOfSegHTOF);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      const Char_t* ud = "S";
      for (Int_t i = 0; i < nseg; ++i) {
        HB1(Form("%s_ADC_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
        HB1(Form("%s_AwT_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
        HB1(Form("%s_AwoT_seg%d%s%s; channel; count", name, i, ud, b), adc_bins);
        HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b),  hr_tdc_bins);
      }
      HB1(Form("%s_HitPat_HT%s; segment; count", name, b),     seg_bins);
      HB1(Form("%s_Multi_HT%s; multiplicity; count", name, b), mul_bins);
    }

    { // KVC: a,b,c,d,S (no U/D). HitPat/Multi from OR/AND.
      auto name = "KVC";
      Double_t nseg = static_cast<Double_t>(NumOfSegKVC);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      for (const auto& uord : std::vector<TString>{"a", "b", "c", "d", "S"}) {
        auto ud = uord.Data();
        for (Int_t i = 0; i < nseg; ++i) {
          HB1(Form("%s_ADC_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
          HB1(Form("%s_AwT_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
          HB1(Form("%s_AwoT_seg%d%s%s; channel; count", name, i, ud, b), adc_bins);
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b),  hr_tdc_bins);
        }
      }
      for (const auto& uord : std::vector<TString>{"OR", "AND"}) {
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s; segment; count", name, ud, b),     seg_bins);
        HB1(Form("%s_Multi_%s%s; multiplicity; count", name, ud, b), mul_bins);
      }
    }

    { // COBO: TDC(U) only
      auto name = "COBO";
      for(Int_t i=0; i<NumOfSegHodo[kCOBO]; ++i)
        HB1(Form("%s_TDC_seg%dU%s; channel; count", name, i, b), hr_tdc_bins_cobo);
    }

    { // SCH: MHTDC only, no ADC, one-side readout
      auto name = "SCH";
      Double_t nseg = static_cast<Double_t>(NumOfSegSCH);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      for(Int_t i=0; i<NumOfSegSCH; ++i){
        HB1(Form("%s_TDC_seg%dU%s; channel; count", name, i, b), mh_tdc_bins_sch);
        HB1(Form("%s_Trailing_seg%dU%s; channel; count", name, i, b), mh_tdc_bins_sch);
        HB1(Form("%s_TOT_seg%dU%s; channel; count", name, i, b), mh_tot_bins);
      }
      HB1(Form("%s_HitPat_OR%s; segment; count", name,  b), seg_bins);
      HB1(Form("%s_Multi_OR%s; multiplicity; count", name, b), mul_bins);
    }
	 
    if(!flag_beam_particle) break;
  }
}

//_____________________________________________________________________________
void
BuildHodoHit(Bool_t flag_beam_particle)
{
  for(const auto& beam: beam::BeamFlagList){
    const Char_t* b = beam.Data();
    { // BHT
      const Char_t* name = "BHT";
      Double_t nseg = static_cast<Double_t>(NumOfSegBHT);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB1(Form("%s_Hit_Time_seg%d%s%s; ns; count", name, i, ud, b),    hr_time_bins);
          HB1(Form("%s_Hit_CTime_seg%d%s%s; ns; count", name, i, ud, b),   hr_time_bins);
          HB1(Form("%s_Hit_TOT_seg%d%s%s; ns; count", name, i, ud, b),     hr_tot_time_bins);
          HB1(Form("%s_Hit_DeltaE_seg%d%s%s; mip; count", name, i, ud, b), de_bins);
        }
        HB1(Form("%s_Hit_MeanTime_seg%d%s; ns; count", name, i, b),  hr_time_bins);
        HB1(Form("%s_Hit_CMeanTime_seg%d%s; ns; count", name, i, b), hr_time_bins);
        HB1(Form("%s_Hit_MeanTOT_seg%d%s; ns; count", name, i, b),   hr_tot_time_bins);
        HB1(Form("%s_Hit_DeltaE_seg%d%s; mip; count", name, i, b),   de_bins);
      }
      HB1(Form("%s_Hit_MeanTime%s; ns; count", name, b),  hr_time_bins);
      HB1(Form("%s_Hit_CMeanTime%s; ns; count", name, b), hr_time_bins);
      HB1(Form("%s_Hit_MeanTOT%s; ns; count", name, b),   hr_tot_time_bins);
      HB1(Form("%s_Hit_DeltaE%s; mip; count", name, b),   de_bins);
      const Double_t hr_time_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        hr_time_bins[0]/10, hr_time_bins[1], hr_time_bins[2] };
      const Double_t hr_tot_time_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        hr_tot_time_bins[0]/5, hr_tot_time_bins[1], hr_tot_time_bins[2] };
      const Double_t de_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        de_bins[0]/10, de_bins[1], de_bins[2] };
      HB2(Form("%s_Hit_MeanTime_vs_HitPat%s; segment; ns", name, b),  hr_time_bins_2d);
      HB2(Form("%s_Hit_CMeanTime_vs_HitPat%s; segment; ns", name, b), hr_time_bins_2d);
      HB2(Form("%s_Hit_MeanTOT_vs_HitPat%s; segment; ns", name, b),   hr_tot_time_bins_2d);
      HB2(Form("%s_Hit_DeltaE_vs_HitPat%s; segment; mip", name, b),   de_bins_2d);
      HB1(Form("%s_Hit_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%s_Hit_Multi%s; multiplicity; count", name, b), mul_bins);
    }

    // Hodoscope (U/D). BAC and KVC use only Sum (S) → "BAC, HTOF, KVC Sum" block below.
    // Cherenkov (SAC3): Hit_Npe / npe_bins. Others: Hit_DeltaE / de_bins.
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      if(ihodo == kBAC || ihodo == kKVC) continue;
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = static_cast<Double_t>(NumOfSegHodo[ihodo]);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      Bool_t is_cherenkov = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Cherenkov);
      Bool_t is_one_side = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::OneSideReadout);
      const Char_t*   dex      = is_cherenkov ? "Npe"    : "DeltaE";
      const Double_t* dex_bins = is_cherenkov ? npe_bins : de_bins;
      const Char_t*   dex_axis = is_cherenkov ? "Npe"    : "mip";
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord: std::vector<TString>{"U", "D"} ){
          if (is_one_side && uord == "D") continue;
          auto ud = uord.Data();
          HB1(Form("%s_Hit_%s_seg%d%s%s; %s; count", name, dex, i, ud, b, dex_axis), dex_bins);
          HB1(Form("%s_Hit_Time_seg%d%s%s; ns; count", name, i, ud, b),  hr_time_bins);
          HB1(Form("%s_Hit_CTime_seg%d%s%s; ns; count", name, i, ud, b), hr_time_bins);
        }
        HB1(Form("%s_Hit_%s_seg%d%s; %s; count", name, dex, i, b, dex_axis), dex_bins);
        HB1(Form("%s_Hit_MeanTime_seg%d%s; ns; count", name, i, b),  hr_time_bins);
        HB1(Form("%s_Hit_CMeanTime_seg%d%s; ns; count", name, i, b), hr_time_bins);
      }
      HB1(Form("%s_Hit_MeanTime%s; ns; count", name, b), hr_time_bins);
      HB1(Form("%s_Hit_CMeanTime%s; ns; count", name, b), hr_time_bins);
      HB1(Form("%s_Hit_%s%s; %s; count", name, dex, b, dex_axis), dex_bins);
      const Double_t hr_time_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        hr_time_bins[0]/10, hr_time_bins[1], hr_time_bins[2] };
      const Double_t de_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        de_bins[0]/10, de_bins[1], de_bins[2] };
      const Double_t npe_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        npe_bins[0]/10, npe_bins[1], npe_bins[2] };
      HB2(Form("%s_Hit_MeanTime_vs_HitPat%s; segment; ns", name, b),  hr_time_bins_2d);
      HB2(Form("%s_Hit_CMeanTime_vs_HitPat%s; segment; ns", name, b), hr_time_bins_2d);
      HB2(Form("%s_Hit_%s_vs_HitPat%s; segment; %s", name, dex, b, dex_axis), is_cherenkov ? npe_bins_2d : de_bins_2d);
      HB1(Form("%s_Hit_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%s_Hit_Multi%s; multiplicity; count", name, b), mul_bins);
    }

    // BAC Sum: S_online (seg4), offline_sum (event). Time/CTime_seg S. No a,b,c,d.
    {
      auto name = "BAC";
      Double_t nseg = static_cast<Double_t>(NumOfSegBAC);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      HB1(Form("%s_Hit_Npe_offline_sum%s; Npe; count", name, b), npe_bins);
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_Hit_Npe_seg%dS_online%s; Npe; count", name, i, b), npe_bins);
        HB1(Form("%s_Hit_Time_seg%dS%s; ns; count", name, i, b),  hr_time_bins);
        HB1(Form("%s_Hit_CTime_seg%dS%s; ns; count", name, i, b), hr_time_bins);
      }
      HB1(Form("%sSum_Hit_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%sSum_Hit_Multi%s; multiplicity; count", name, b), mul_bins);
    }

    // HTOF Sum (S, kExtra): Hit_DeltaE, Time/CTime_seg S, Sum_Hit_HitPat/Multi
    {
      auto name = "HTOF";
      Double_t nseg = static_cast<Double_t>(NumOfSegHTOF);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_Hit_DeltaE_seg%dS%s; mip; count", name, i, b), de_bins);
        HB1(Form("%s_Hit_Time_seg%dS%s; ns; count", name, i, b),    hr_time_bins);
        HB1(Form("%s_Hit_CTime_seg%dS%s; ns; count", name, i, b),   hr_time_bins);
      }
      HB1(Form("%sSum_Hit_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%sSum_Hit_Multi%s; multiplicity; count", name, b), mul_bins);
    }

    // KVC Sum (kExtra): 
    // Hit_Npe_seg S_offline (raw-based per seg), S_online+a,b,c,d. 
    // Time/CTime_seg S, Sum_Hit_HitPat/Multi
    {
      auto name = "KVC";
      Double_t nseg = static_cast<Double_t>(NumOfSegKVC);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      const Char_t* abcd[4] = {"a", "b", "c", "d"};
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_Hit_Npe_seg%dS_offline%s; Npe; count", name, i, b), npe_bins);
        HB1(Form("%s_Hit_Npe_seg%dS_online%s; Npe; count", name, i, b),  npe_bins);
        HB1(Form("%s_Hit_Npe_seg%dS_offline_T1%s; Npe; count", name, i, b), npe_bins);
        HB1(Form("%s_Hit_Npe_seg%dS_online_T1%s; Npe; count", name, i, b),  npe_bins);
        for(Int_t c=0; c<4; ++c)
          HB1(Form("%s_Hit_Npe_seg%d%s%s; Npe; count", name, i, abcd[c], b), npe_bins);
        HB1(Form("%s_Hit_Time_seg%dS%s; ns; count", name, i, b),  hr_time_bins);
        HB1(Form("%s_Hit_CTime_seg%dS%s; ns; count", name, i, b), hr_time_bins);
      }
      HB1(Form("%sSum_Hit_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%sSum_Hit_Multi%s; multiplicity; count", name, b), mul_bins);
    }

    // TOF / BTOF / FTOF
    const Double_t phcbins2d[6] = {100, -0.5, 4.5, 100, -10., 10.};
    { // TOF-HTOF
      for(Int_t i=0; i<NumOfSegHTOF; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("HTOF_seg%d%s_TOF_vs_DeltaE%s; mip; ns", i, ud, b),  phcbins2d);
          HB2(Form("HTOF_seg%d%s_CTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("HTOF_TOF_vs_DeltaE%s; mip; ns", b),  phcbins2d);
      HB2(Form("HTOF_CTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
    }

    { // TOF-T1
      for(Int_t i=0; i<NumOfSegT1; ++i){
        for(const auto& uord : std::vector<TString>{"U"}){
          const Char_t* ud = uord.Data();
          HB2(Form("T1_seg%d%s_TOF_vs_DeltaE%s; mip; ns", i, ud, b),  phcbins2d);
          HB2(Form("T1_seg%d%s_CTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("T1_TOF_vs_DeltaE%s; mip; ns", b),  phcbins2d);
      HB2(Form("T1_CTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
    }

    { // BTOF
      for(Int_t i=0; i<NumOfSegHodo[kBH2]; ++i){
        HB1(Form("T0_seg%d_TimeOffset%s; ns; count", i, b), 2000, -10, 10);
      }
      for(Int_t i=0; i<NumOfSegBHT; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("BHT_seg%d%s_BTOF_vs_DeltaE%s; mip; ns", i, ud, b),  phcbins2d);
          HB2(Form("BHT_seg%d%s_CBTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("BHT_BTOF_vs_DeltaE%s; mip; ns", b),  phcbins2d);
      HB2(Form("BHT_CBTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
      for(Int_t i=0; i<NumOfSegHodo[kBH2]; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("T0_seg%d%s_BTOF_vs_DeltaE%s; mip; ns", i, ud, b),  phcbins2d);
          HB2(Form("T0_seg%d%s_CBTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("T0_BTOF_vs_DeltaE%s; mip; ns", b),  phcbins2d);
      HB2(Form("T0_CBTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
    }

    { // FTOF
      for(const auto& id: std::vector<Int_t>{kCVC}){
        const Char_t* n = NameHodo[id];
        for(Int_t i=0; i<NumOfSegHodo[id]; ++i){
          for(const auto& uord : std::vector<TString>{"U", "D"}){
            const Char_t* ud = uord.Data();
            HB2(Form("%s_seg%d%s_FTOF_vs_DeltaE%s; mip; ns", n, i, ud, b),  phcbins2d);
            HB2(Form("%s_seg%d%s_CFTOF_vs_DeltaE%s; mip; ns", n, i, ud, b), phcbins2d);
          }
        }
        HB2(Form("%s_FTOF_vs_DeltaE%s; mip; ns", n, b),  phcbins2d);
        HB2(Form("%s_CFTOF_vs_DeltaE%s; mip; ns", n, b), phcbins2d);
      }
      for(Int_t i=0; i<NumOfSegHodo[kBH2]; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("T0_seg%d%s_FTOF_vs_DeltaE%s; mip; ns", i, ud, b),  phcbins2d);
          HB2(Form("T0_seg%d%s_CFTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("T0_FTOF_vs_DeltaE%s; mip; ns", b),  phcbins2d);
      HB2(Form("T0_CFTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
    }
    if(!flag_beam_particle) break;
  }
}

//_____________________________________________________________________________
void
BuildHodoCluster(Bool_t flag_beam_particle)
{
  for(const auto& beam: beam::BeamFlagList){
    const Char_t* b = beam.Data();
    const Double_t size_bins[3] = {10+1, -0.5, 10+0.5};
    { // BHT
      const Char_t* name = "BHT";
      Double_t nseg = static_cast<Double_t>(NumOfSegBHT);
      const Double_t hr_time_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        hr_time_bins[0]/10, hr_time_bins[1], hr_time_bins[2] };
      const Double_t de_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        de_bins[0]/10, de_bins[1], de_bins[2] };
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      HB2(Form("%s_Cl_MeanTime_vs_HitPat%s; segment; ns", name, b),  hr_time_bins_2d);
      HB2(Form("%s_Cl_CMeanTime_vs_HitPat%s; segment; ns", name, b), hr_time_bins_2d);
      HB2(Form("%s_Cl_TimeDiff_vs_HitPat%s; segment; ns", name, b),  hr_time_bins_2d);
      HB2(Form("%s_Cl_DeltaE_vs_HitPat%s; segment; mip", name, b),   de_bins_2d);
      HB1(Form("%s_Cl_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%s_Cl_Multi%s; multiplicity; count", name, b), mul_bins);
      HB1(Form("%s_Cl_Size%s; size; count", name, b),          size_bins);
    }
    
    // Hodoscope (exclude NoCluster: BAC, T1, SAC3, SFV, COBO). KVC: Cl_Npe_vs_HitPat.
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      if(HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::NoCluster)) continue;
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = static_cast<Double_t>(NumOfSegHodo[ihodo]);
      Bool_t is_cherenkov = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Cherenkov);
      const Double_t hr_time_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        hr_time_bins[0]/10, hr_time_bins[1], hr_time_bins[2] };
      const Double_t de_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        de_bins[0]/10, de_bins[1], de_bins[2] };
      const Double_t npe_bins_2d[6] = { nseg, -0.5, nseg-0.5,
        npe_bins[0]/10, npe_bins[1], npe_bins[2] };
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      HB2(Form("%s_Cl_MeanTime_vs_HitPat%s; segment; ns", name, b), hr_time_bins_2d);
      HB2(Form("%s_Cl_CMeanTime_vs_HitPat%s; segment; ns", name, b), hr_time_bins_2d);
      HB2(Form("%s_Cl_TimeDiff_vs_HitPat%s; segment; ns", name, b), hr_time_bins_2d);
      HB2(Form("%s_Cl_%s_vs_HitPat%s; segment; %s", name, is_cherenkov?"Npe":"DeltaE", b, is_cherenkov?"Npe":"mip"), is_cherenkov ? npe_bins_2d : de_bins_2d);
      HB1(Form("%s_Cl_HitPat%s; segment; count", name, b),     seg_bins);
      HB1(Form("%s_Cl_Multi%s; multiplicity; count", name, b), mul_bins);
      HB1(Form("%s_Cl_Size%s; size; count", name, b),          size_bins);
    }

    // BTOF
    const Double_t ctime0_bins[3] = {400, -4, 4};
    const Double_t cbtof0_bins[3] = {600, -20, 10};
    const Double_t tof_vs_de_bins[6] = {200, 0, 4, 200, -4, 4};
    HB1(Form("CTime0%s; ns; count", b), ctime0_bins);
    HB1(Form("CBtof0%s; ns; count", b), cbtof0_bins);
    HB2(Form("CBtof0_vs_deT0Seg%s; mip; ns", b), tof_vs_de_bins);
    HB2(Form("CBtof0_vs_deBtof0Seg%s; mip; ns", b), tof_vs_de_bins);
    
    // FTOF
    const Double_t cftof0_bins[3] = {600, -10, 30};
    HB1(Form("CFtof0%s; ns; count", b), cftof0_bins);
    HB2(Form("CFtof0_vs_deT0Seg%s; mip; ns", b), tof_vs_de_bins);
    HB2(Form("CFtof0_vs_deFtof0Seg%s; mip; ns", b), tof_vs_de_bins);
    if(!flag_beam_particle) break;
  }
}

//_____________________________________________________________________________
void
BuildDCRaw(const TString& dcname, Bool_t flag_beam_particle)
{
  const auto& digit_info = gUConf.get_digit_info();
  // const auto& plane_names = digit_info.get_name_list(m_detector_id);
  // m_plane_name = plane_names.at(plane_id);
  // m_dcgeom_layer = gGeom.GetLayerId(m_detector_name+"-"+m_plane_name);

  for(const auto& beam: beam::BeamFlagList){
    const Char_t* b = beam.Data();
    // for(Int_t idc=0; idc<=kBPC2; ++idc){
    //   const Char_t* name = NameDC[idc].Data();
    for (const auto& name_str : DCNameList.at(dcname)) {
      const auto name = name_str.Data();
      const auto detector_id = digit_info.get_device_id(name);
      const Int_t n_plane = digit_info.get_n_plane(detector_id);
      const Double_t n_wire = digit_info.get_n_ch(detector_id);
      const Double_t pat_bins[3]    = {n_wire,   -0.5, n_wire-0.5};
      const Double_t mul_bins[3]    = {n_wire+1, -0.5, n_wire+0.5};
      const Double_t tdc_bins_2d[6] = {n_wire,   -0.5, n_wire-0.5,
        mh_tdc_bins[0], mh_tdc_bins[1], mh_tdc_bins[2]};
      const Double_t tot_bins_2d[6] = {n_wire,   -0.5, n_wire-0.5,
        mh_tot_bins[0], mh_tot_bins[1], mh_tot_bins[2]};
      const Double_t tdc_tot_bins_2d[6] = {
        mh_tdc_bins[0], mh_tdc_bins[1], mh_tdc_bins[2],
        mh_tot_bins[0], mh_tot_bins[1], mh_tot_bins[2]};
      for (Int_t plane = 0; plane < n_plane; ++plane) {
        for (const auto& tot_cut : std::vector<TString>{"", "C"}) {
          const Char_t* suffix = tot_cut.Data();
          HB1(Form("%s_%sTDC_plane%d%s; channel; count", name, suffix, plane, b), mh_tdc_bins);
          HB1(Form("%s_%sTDC1st_plane%d%s; channel; count", name, suffix, plane, b), mh_tdc_bins);
          HB1(Form("%s_%sTrailing_plane%d%s; channel; count", name, suffix, plane, b), mh_tdc_bins);
          HB1(Form("%s_%sTrailing1st_plane%d%s; channel; count", name, suffix, plane, b), mh_tdc_bins);
          HB1(Form("%s_%sTOT_plane%d%s; channel; count", name, suffix, plane, b), mh_tot_bins);
          HB1(Form("%s_%sTOT1st_plane%d%s; channel; count", name, suffix, plane, b), mh_tot_bins);
          HB1(Form("%s_%sHitPat_plane%d%s; wire; count", name, suffix, plane, b), pat_bins);
          HB1(Form("%s_%sMulti_plane%d%s; multiplicity; count", name, suffix, plane, b), mul_bins);
          HB2(Form("%s_%sTOT_vs_TDC_plane%d%s; segment; channel", name, suffix, plane, b), tdc_tot_bins_2d);
          HB2(Form("%s_%sTDC_vs_HitPat_plane%d%s; segment; channel", name, suffix, plane, b), tdc_bins_2d);
          HB2(Form("%s_%sTDC1st_vs_HitPat_plane%d%s; segment; channel", name, suffix, plane, b), tdc_bins_2d);
          HB2(Form("%s_%sTrailing_vs_HitPat_plane%d%s; segment; channel", name, suffix, plane, b), tdc_bins_2d);
          HB2(Form("%s_%sTrailing1st_vs_HitPat_plane%d%s; segment; channel", name, suffix, plane, b), tdc_bins_2d);
          HB2(Form("%s_%sTOT_vs_HitPat_plane%d%s; segment; channel", name, suffix, plane, b), tot_bins_2d);
          HB2(Form("%s_%sTOT1st_vs_HitPat_plane%d%s; segment; channel", name, suffix, plane, b), tot_bins_2d);
        }
      }
    }
    if(!flag_beam_particle) break;
  }
}

//_____________________________________________________________________________
void
BuildDCHit(const TString& dcname, Bool_t flag_beam_particle)
{
  const auto& digit_info = gUConf.get_digit_info();
  for(const auto& beam: beam::BeamFlagList){
    const Char_t* b = beam.Data();
    for (const auto& name_str : DCNameList.at(dcname)) {
      const auto name = name_str.Data();
      const auto detector_id = digit_info.get_device_id(name);
      const Int_t n_plane = digit_info.get_n_plane(detector_id);
      const Double_t n_wire = digit_info.get_n_ch(detector_id);
      const Double_t pat_bins[3] = {n_wire, -0.5, n_wire-0.5};
      const Double_t pat_bins_2d[6] = {n_wire, -0.5, n_wire-0.5, n_wire, -0.5, n_wire-0.5};
      const Double_t d_pat_bins[3] = {n_wire * 2, -n_wire-0.5, n_wire-0.5};
      const Double_t mul_bins[3] = {n_wire+1, -0.5, n_wire+0.5};
      const Double_t dt_bins[3] = {600, -100., 400.};
      const Double_t dl_bins[3] = {60, -1.0, 5.0};
      const Double_t dt_bins_2d[6] = {n_wire, -0.5, n_wire-0.5, dt_bins[0], dt_bins[1], dt_bins[2]};
      const Double_t dl_bins_2d[6] = {n_wire, -0.5, n_wire-0.5, dl_bins[0], dl_bins[1], dl_bins[2]};
      const Double_t dt_tot_bins_2d[6] = {dt_bins[0], dt_bins[1], dt_bins[2], mh_tot_bins[0], mh_tot_bins[1], mh_tot_bins[2]};
      for (Int_t plane = 0; plane < n_plane; ++plane) {
        HB1(Form("%s_Hit_DriftTime_plane%d%s; ns; count", name, plane, b), dt_bins);
        HB1(Form("%s_Hit_DriftLength_plane%d%s; mm; count", name, plane, b), dl_bins);
        HB2(Form("%s_Hit_TOT_vs_DriftTime_plane%d%s; segment; ns", name, plane, b), dt_tot_bins_2d);
        HB2(Form("%s_Hit_DriftTime_vs_HitPat_plane%d%s; segment; ns", name, plane, b), dt_bins_2d);
        HB2(Form("%s_Hit_DriftLength_vs_HitPat_plane%d%s; segment; mm", name, plane, b), dl_bins_2d);
        HB1(Form("%s_Hit_HitPat_plane%d%s; wire; count", name, plane, b), pat_bins);
        if (plane % 2 == 0) {
	  HB2(Form("%s_Raw_HitPat_Pairplane%d%d%s; wire [plane %d]; wire [plane %d]", name, plane, plane+1, b, plane, plane+1), pat_bins_2d);
          HB2(Form("%s_Hit_HitPat_Pairplane%d%d%s; wire [plane %d]; wire [plane %d]", name, plane, plane+1, b, plane, plane+1), pat_bins_2d);
          HB1(Form("%s_Hit_HitPat_PP_Sub%d%d%s; wire of plane %d-wire of plane %d;count", name, plane, plane+1, b, plane, plane+1), d_pat_bins);
        }
        HB1(Form("%s_Hit_Multi_plane%d%s; multiplicity; count", name, plane, b), mul_bins);
      }
    }
    if(!flag_beam_particle) break;
  }
}

//_____________________________________________________________________________
void
BuildDCTrack(const TString& dcname, Bool_t flag_beam_particle)
{
  const auto& digit_info = gUConf.get_digit_info();
  for(const auto& beam: beam::BeamFlagList){
    const Char_t* b = beam.Data();
    const Double_t nhit_bins[3]  = {20, -0.5, 19.5};
    const Double_t chisq_bins[3] = {200, 0.0, 40.0};
    const Double_t xy0_bins[3]   = {200, -500.0, 500.0};
    const Double_t uv0_bins[3]   = {200, -0.5, 0.5};
    HB1(Form("%sTrack_NHit%s; ; count", dcname.Data(), b), nhit_bins);
    HB1(Form("%sTrack_ChiSquare%s; ; count", dcname.Data(), b), chisq_bins);
    HB1(Form("%sTrack_X0%s; ; count", dcname.Data(), b), xy0_bins);
    HB1(Form("%sTrack_Y0%s; ; count", dcname.Data(), b), xy0_bins);
    HB1(Form("%sTrack_U0%s; ; count", dcname.Data(), b), uv0_bins);
    HB1(Form("%sTrack_V0%s; ; count", dcname.Data(), b), uv0_bins);
    for (const auto& name_str : DCNameList.at(dcname)) {
      const auto name = name_str.Data();
      const auto detector_id = digit_info.get_device_id(name);
      const Int_t n_plane = digit_info.get_n_plane(detector_id);
      const Double_t n_wire = digit_info.get_n_ch(detector_id);
      const Double_t pat_bins[3] = {n_wire, -0.5, n_wire-0.5};
      const Double_t dt_bins[3] = {600, -100., 400.};
      const Double_t dl_bins[3] = {60, -1.0, 5.0};
      const Double_t dt_bins_2d[6] = {n_wire, -0.5, n_wire-0.5, dt_bins[0], dt_bins[1], dt_bins[2]};
      const Double_t dl_bins_2d[6] = {n_wire, -0.5, n_wire-0.5, dl_bins[0], dl_bins[1], dl_bins[2]};
      const Double_t res_bins[3] = {400, -2.0, 2.0};
      const Double_t res_dl_bins_2d[6] = {200, -3., 3., 200, -2.0, 2.0};
      for (Int_t plane = 0; plane < n_plane; ++plane) {
        HB1(Form("%s_Track_DriftTime_plane%d%s; ns; count", name, plane, b), dt_bins);
        HB1(Form("%s_Track_DriftLength_plane%d%s; mm; count", name, plane, b), dl_bins);
        HB2(Form("%s_Track_DriftTime_vs_HitPat_plane%d%s; segment; ns", name, plane, b), dt_bins_2d);
        HB2(Form("%s_Track_DriftLength_vs_HitPat_plane%d%s; segment; mm", name, plane, b), dl_bins_2d);
        HB1(Form("%s_Track_HitPat_plane%d%s; wire; count", name, plane, b), pat_bins);
        HB1(Form("%s_Track_Residual_plane%d%s; mm; count", name, plane, b), res_bins);
        HB2(Form("%s_Track_Residual_vs_DriftLength_plane%d%s; mm; count", name, plane, b), res_dl_bins_2d);
      }
    }
    if(!flag_beam_particle) break;
  }
}

//_____________________________________________________________________________
void
BuildDAQ()
{
  std::vector<Int_t> vme_fe_id;
  std::vector<Int_t> hul_fe_id;
  std::vector<Int_t> vea0c_fe_id;
  HB1("EB_DataSize; words", 2000, 0, 20000);
  for(auto&& c : gUnpacker.get_root()->get_child_list()){
    if (!c.second) continue;
    TString name = c.second->get_name();
    auto node_id = c.second->get_id();
    if(name.Contains("vme"))
      vme_fe_id.push_back(node_id);
    if(name.Contains("hul"))
      hul_fe_id.push_back(node_id);
    if(name.Contains("veasiroc"))
      vea0c_fe_id.push_back(node_id);
  }
  HB2("FE_VME_DataSize; ; words",
            vme_fe_id.size(), 0, vme_fe_id.size(), 100, 0, 1000);
  for(Int_t i=0, n=vme_fe_id.size(); i<n; ++i){
    auto h1 = gDirectory->Get<TH2>("FE_VME_DataSize");
    h1->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(vme_fe_id[i], 16));
  }
  HB2("FE_HUL_DataSize; ; words",
            hul_fe_id.size(), 0, hul_fe_id.size(), 100, 0, 2000);
  for(Int_t i=0, n=hul_fe_id.size(); i<n; ++i){
    auto h1 = gDirectory->Get<TH2>("FE_HUL_DataSize");
    h1->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(hul_fe_id[i], 16));
  }
  HB2("FE_VEASIROC_DataSize; ; words",
            vea0c_fe_id.size(), 0, vea0c_fe_id.size(), 100, 0, 1000);
  for(Int_t i=0, n=vea0c_fe_id.size(); i<n; ++i){
    auto h1 = gDirectory->Get<TH2>("FE_VEASIROC_DataSize");
    h1->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(vea0c_fe_id[i], 16));
  }
}

//_____________________________________________________________________________
void
BuildTPCHit()
{
  const Int_t    n_bin_adc     = 4096;
  const Double_t min_adc       = 0.;
  const Double_t max_adc       = 4096.;
  const Int_t    n_bin_rms     = 1000;
  const Double_t min_rms       = 0.;
  const Double_t max_rms       = 1000.;
  const Int_t    n_bin_de      = 1000;
  const Double_t min_de        = 0.;
  const Double_t max_de        = 1000.;
  const Int_t    n_bin_chisqr  = 1000;
  const Double_t min_chisqr    = 0.;
  const Double_t max_chisqr    = 1000.;
  const Int_t    n_bin_time    = 1000;
  const Double_t min_time      = -8000.;
  const Double_t max_time      = 8000.;
  const Int_t    n_bin_dl      = 800;
  const Double_t min_dl        = -400.;
  const Double_t max_dl        = 400.;
  const Int_t    n_bin_sigma   = 500;
  const Double_t min_sigma     = 0.;
  const Double_t max_sigma     = 50.;
  const Int_t    n_time_bucket = 170;

  HB1("TPC_Multiplicity_Raw",     NumOfPadTPC+1,   0.,      NumOfPadTPC+1.);
  HB1("TPC_Multiplicity_Cor",     NumOfPadTPC+1,   0.,      NumOfPadTPC+1.);
  HB1("TPC_FADC_Mean",            n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Max",             n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_RMS",             n_bin_rms,       min_rms, max_rms);
  HB1("TPC_FADC_LocMax",          n_time_bucket+1, 0,       n_time_bucket+1);
  HB1("TPC_FADC_Min",             n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Cor_Mean",        n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Cor_Max",         n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Cor_RMS",         n_bin_rms,       min_rms, max_rms);
  HB1("TPC_FADC_Cor_LocMax",      n_time_bucket+1, 0.,      n_time_bucket+1.);
  HB1("TPC_FADC_Cor_Min",         n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Baseline_p0",     n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Baseline_p1",     120,             -6.,     6.);
  HB1("TPC_FADC_Baseline_p2",     120,             -12.,    12.);
  HB1("TPC_FADC_Baseline_Mean",   n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Baseline_Max",    n_bin_adc,       min_adc, max_adc);
  HB1("TPC_FADC_Baseline_RMS",    n_bin_rms,       min_rms, max_rms);
  HB1("TPC_FADC_Baseline_LocMax", n_time_bucket+1, 0.,      n_time_bucket+1.);
  HB1("TPC_FADC_Baseline_Min",    n_bin_adc,       min_adc, max_adc);
  HB2("TPC_FADC_Baseline",        n_time_bucket+1, 0,       n_time_bucket+1, 
                                  n_bin_adc,       min_adc, max_adc);

  HB1("TPC_Multiplicity_TPCHit",   NumOfPadTPC+1,        0.,         NumOfPadTPC+1.);
  HB1("TPC_Pedestal",              n_bin_adc,            min_adc,    max_adc);
  HB1("TPC_DeltaE",                n_bin_de,             min_de,     max_de);
  HB1("TPC_RMS",                   n_bin_rms,            min_rms,    max_rms);
  HB1("TPC_Time",                  (n_time_bucket+1)*30, 0.,         n_time_bucket+1.);
  HB1("TPC_Chisqr",                n_bin_chisqr,         min_chisqr, max_chisqr);
  HB1("TPC_CDeltaE",               n_bin_de,             min_de,     max_de);
  HB1("TPC_CTime",                 n_bin_time,           min_time,   max_time);
  HB1("TPC_DriftLength",           n_bin_dl,             min_dl,     max_dl);
  HB1("TPC_sigma",                 n_bin_sigma,          min_sigma,  max_sigma);
  HB2("TPC_sigma%%de",             n_bin_de,             min_de,     max_de, 
                                   n_bin_sigma,          min_sigma,  max_sigma);
  HB2("TPC_time%%de",              n_bin_de,             min_de,     max_de,
                                   n_bin_time,           min_time,   max_time);

  HB2("TPC_FADC_Before", n_time_bucket+1, 0.,           n_time_bucket+1.,
                         n_bin_adc,       min_adc,      max_adc);
  HB2("TPC_FADC_After",  n_time_bucket+1, 0.,           n_time_bucket+1.,
                         n_bin_adc,       min_adc-500., max_adc-500.);
  HB2("TPC_FADC_Good",   n_time_bucket+1, 0.,           n_time_bucket+1., 
                         n_bin_adc,       min_adc-500.,  max_adc-500.);
  HB2("TPC_FADC_Noise",  n_time_bucket+1, 0.,           n_time_bucket+1., 
                         n_bin_adc,       min_adc,      max_adc);
  HB2("TPC_FADC_Frame",  n_time_bucket+1, 0.,           n_time_bucket+1., 
                         n_bin_adc,       min_adc,      max_adc);

  HB1("TPC_FADC_Noise_Max",        n_bin_adc, min_adc, max_adc);
  HB1("TPC_FADC_Noise_RMSfront",   n_bin_rms, min_rms, max_rms);
  HB1("TPC_FADC_Noise_RMSmiddle",  n_bin_rms, min_rms, max_rms);
  HB1("TPC_FADC_Noise_Adcdiff",    1000,      -100.,   900.);

  HB1("TPC_Clock_TDC",  200000, 0.,    2000000.);
  HB1("TPC_Clock_Time", 20000,  -100., 100.);

  HB2Poly("TPC_Raw_Max_Poly", TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPC_Raw_ADC_Poly", TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPC_Raw_Mean_Poly",TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPC_Cor_Max_Poly", TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPC_Cor_ADC_Poly", TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPC_Cor_Mean_Poly",TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPC_HitPat_Noise",    TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPC_HitPat_Baseline", TPC_EVENT_DISPLAY_BINS);
  tpc::InitializeHistograms("TPC_Raw_Max_Poly");
  tpc::InitializeHistograms("TPC_Raw_ADC_Poly");
  tpc::InitializeHistograms("TPC_Raw_Mean_Poly");
  tpc::InitializeHistograms("TPC_Cor_Max_Poly");
  tpc::InitializeHistograms("TPC_Cor_ADC_Poly");
  tpc::InitializeHistograms("TPC_Cor_Mean_Poly");
  tpc::InitializeHistograms("TPC_HitPat_Noise");
  tpc::InitializeHistograms("TPC_HitPat_Baseline");
}

//_____________________________________________________________________________
void
BuildTPCTrackingCommon()
{
  HB1("TPCTrk_Num_Track;N_{track} (TPC Track);Counts", TPC_BINS_MULT);
  HB1("TPCTrk_Num_TrackHits;N_{hits} (TPC Track);Counts", TPC_BINS_HITS);
  HB1("TPCTrk_Chisqr;#chi^{2} (TPC Track);Counts", TPC_BINS_CHISQR);
}

//_____________________________________________________________________________
void
BuildTPCLineTrackParam()
{
  HB1("TPCTrk_X0;X_{0} (TPC Track) [mm];Counts", TPC_BINS_POS_X);
  HB1("TPCTrk_Y0;Y_{0} (TPC Track) [mm];Counts", TPC_BINS_POS_Y);
  HB1("TPCTrk_U0;U_{0} (TPC Track) (dY/dX);Counts", TPC_BINS_SLOPE);
  HB1("TPCTrk_V0;V_{0} (TPC Track) (dZ/dX);Counts", TPC_BINS_SLOPE);

  HB2("TPCTrk_U0_vs_X0;X_{0} [mm];U_{0} (dY/dX)", TPC_BINS_POS_X, TPC_BINS_SLOPE);
  HB2("TPCTrk_V0_vs_Y0;Y_{0} [mm];V_{0} (dZ/dX)", TPC_BINS_POS_Y, TPC_BINS_SLOPE);
  HB2("TPCTrk_Y0_vs_X0;X_{0} [mm];Y_{0} [mm]", TPC_BINS_POS_X, TPC_BINS_POS_Y);
  HB2("TPCTrk_atanU0_vs_X0;X_{0} [mm];atan(U_{0}) [deg]", TPC_BINS_POS_X, TPC_BINS_ATAN_DEG);
  HB2("TPCTrk_atanV0_vs_Y0;Y_{0} [mm];atan(V_{0}) [deg]", TPC_BINS_POS_Y, TPC_BINS_ATAN_DEG);
}

//_____________________________________________________________________________
void
BuildTPCTracking(Bool_t calib_flag)
{
  HB1("TPCTrk_Hough_Dist;Hough distance [mm];Counts", TPC_BINS_HOUGH);
  HB1("TPCTrk_Hough_DistY;Hough distance Y [mm];Counts", TPC_BINS_HOUGH);
  HB1("TPCTrk_Num_Iter;N_{iter};Counts", TPC_BINS_ITER);
  HB1("TPCTrk_Fitting_Flag;Flag;Counts", TPC_BINS_FLAG);
  HB1("TPCTrk_Searching_Time;Time [ms];Counts", TPC_BINS_TIME_MS);
  HB1("TPCTrk_Fitting_Time;Time [ms];Counts", TPC_BINS_TIME_MS);
  HB1("TPCTrk_Minuit_Status;Status;Counts", TPC_BINS_MINUIT);

  HB1("TPCTrk_Layer;Layer;Counts", TPC_BINS_LAYER);
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    HB1(Form("TPCHit_HitPat_Layer%02d;Row;Counts", layer), n_pad, -0.5, n_pad - 0.5);
    HB1(Form("TPCHit_Xhit_Layer%02d;X (TPC Hit) [mm];Counts", layer), TPC_BINS_POS_X);
    HB1(Form("TPCTrk_Res_Layer%02d;Residual (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_ABS);
    HB2(Form("TPCTrk_Res_vs_Xhit_Layer%02d;X (TPC Hit) [mm];Residual (TPC Hit - TPC Track) [mm]", layer), TPC_BINS_POS_X, TPC_BINS_RES_TRACK);
    HB2(Form("TPCHit_Yhit_vs_Xtrk_Layer%02d;X (TPC Track) [mm];Y (TPC Hit) [mm]", layer), TPC_BINS_POS_X, TPC_BINS_POS_Y);
    HB1(Form("TPCTrk_ResX_Layer%02d;Residual X (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_TRACK);
    HB1(Form("TPCTrk_ResY_Layer%02d;Residual Y (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_TRACK);
    HB1(Form("TPCTrk_ResZ_Layer%02d;Residual Z (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_TRACK);
  }

  HB1("TPCCl_Size;Cluster size;Counts", TPC_BINS_CL_SIZE);
  HB1("TPCCl_dE;Cluster dE;Counts", TPC_BINS_DE);
  HB2("Transverse_Diffusion;X_{cluster}-X_{pad} [mm];A/A_{sum}", TPC_BINS_DIFF_PAD, TPC_BINS_CLUSTER_RATIO);

  HB2Poly("TPCTrk_HitPat", TPC_EVENT_DISPLAY_BINS);
  tpc::InitializeHistograms("TPCTrk_HitPat");

  HB2("TPCTrk_Row_vs_Layer;Layer;Row", TPC_BINS_ROW_VS_LAYER);

  HB2("TPCTrk_ResX_vs_Layer_Trk;Layer;Residual X (TPC Cluster - TPC Track) [mm]", TPC_BINS_LAYER, TPC_BINS_RES_HIT);
  HB2("TPCTrk_ResY_vs_Layer_Trk;Layer;Residual Y (TPC Cluster - TPC Track) [mm]", TPC_BINS_LAYER, TPC_BINS_RES_HIT);
  HB2("TPCTrk_ResZ_vs_Layer_Trk;Layer;Residual Z (TPC Cluster - TPC Track) [mm]", TPC_BINS_LAYER, TPC_BINS_RES_HIT);
  HB2("TPCCl_dE_vs_Layer;Layer;Cluster dE", TPC_BINS_LAYER, TPC_BINS_DE);
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    HB1(Form("TPCCl_Size_Layer%02d;Cluster size;Counts", layer), TPC_BINS_CL_SIZE);
    HB1(Form("TPCCl_dE_Layer%02d;Cluster dE;Counts", layer), TPC_BINS_DE);
    HB2(Form("Transverse_Diffusion_Layer%02d;X_{cluster}-X_{pad} [mm];A/A_{sum}", layer), TPC_BINS_DIFF_PAD, TPC_BINS_CLUSTER_RATIO);
    HB2(Form("TPCTrk_ResY_vs_Y_Layer%02d;Y_{TPC Track} (local) [mm];Residual Y (TPC Cluster - TPC Track) [mm]", layer), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
  }
  if (calib_flag) {
    for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
      const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
      for (Int_t row = 0; row < n_pad; ++row) {
        HB1(Form("TPCTrk_ResY_Layer%02d_Row%03d;Residual Y (TPC Cluster - TPC Track) [mm];Counts", layer, row), TPC_BINS_RES_HIT);
        HB2(Form("TPCTrk_ResY_vs_Y_Layer%02d_Row%03d;Y_{TPC Track} (local) [mm];Residual Y (TPC Cluster - TPC Track) [mm]", layer, row), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
        HB1(Form("TPCCl_dE_Layer%02d_Row%03d;Cluster dE;Counts", layer, row), TPC_BINS_DE);
      }
    }
  }

  BuildCoBoClockTime(kCoBoClockTime_Track);
}

//_____________________________________________________________________________
void
BuildTPCBcOutTracking(Bool_t calib_flag)
{
  HB1("BcOut_Num_Track;N_{track};Counts", TPC_BINS_EXT_MULT);
  HB1("BcOut_Num_TrackHits;N_{hits};Counts", TPC_BINS_HITS);
  HB1("BcOut_Chisqr;#chi^{2};Counts", TPC_BINS_EXT_CHISQR);
  HB1("BcOut_X0;X_{0} [mm];Counts", TPC_BINS_POS_X);
  HB1("BcOut_Y0;Y_{0} [mm];Counts", TPC_BINS_POS_Y);
  HB1("BcOut_U0;U_{0} (dY/dX);Counts", TPC_BINS_SLOPE);
  HB1("BcOut_V0;V_{0} (dZ/dX);Counts", TPC_BINS_SLOPE);
  HB1("BcOut_XTgt;X_{tgt} [mm];Counts", TPC_BINS_POS_X);
  HB1("BcOut_YTgt;Y_{tgt} [mm];Counts", TPC_BINS_POS_Y);
  HB1("BcOut_UTgt;U_{tgt} (dY/dX);Counts", TPC_BINS_SLOPE);
  HB1("BcOut_VTgt;V_{tgt} (dZ/dX);Counts", TPC_BINS_SLOPE);
  HB2("BcOut_UTgt_vs_XTgt;X_{tgt} [mm];U_{tgt} (dY/dX)", TPC_BINS_POS_X, TPC_BINS_SLOPE);
  HB2("BcOut_VTgt_vs_YTgt;Y_{tgt} [mm];V_{tgt} (dZ/dX)", TPC_BINS_POS_Y, TPC_BINS_SLOPE);
  HB2("BcOut_YTgt_vs_XTgt;X_{tgt} [mm];Y_{tgt} [mm]", TPC_BINS_POS_X, TPC_BINS_POS_Y);

  HB2("BcOut_vs_TPC_XTgt;TPC X_{tgt} [mm];BcOut X_{tgt} [mm]", TPC_BINS_EXT_TGT_POS, TPC_BINS_EXT_TGT_POS);
  HB2("BcOut_vs_TPC_YTgt;TPC Y_{tgt} [mm];BcOut Y_{tgt} [mm]", TPC_BINS_POS_Y, TPC_BINS_POS_Y);
  HB2("BcOut_vs_TPC_UTgt;TPC U_{tgt} (dY/dX);BcOut U_{tgt} (dY/dX)", TPC_BINS_EXT_TGT_SLOPE_U, TPC_BINS_EXT_TGT_SLOPE_U);
  HB2("BcOut_vs_TPC_VTgt;TPC V_{tgt} (dZ/dX);BcOut V_{tgt} (dZ/dX)", TPC_BINS_EXT_TGT_SLOPE_V, TPC_BINS_EXT_TGT_SLOPE_V);
  HB1("TPCTrk_ResX_Tgt;Residual X (TPC Track - BcOut) [mm];Counts", TPC_BINS_RES_TGT);
  HB1("TPCTrk_ResY_Tgt;Residual Y (TPC Track - BcOut) [mm];Counts", TPC_BINS_RES_TGT);
  HB1("TPCTrk_ResU_Tgt;Residual U (TPC Track - BcOut);Counts", TPC_BINS_RES_TGT_SLOPE);
  HB1("TPCTrk_ResV_Tgt;Residual V (TPC Track - BcOut);Counts", TPC_BINS_RES_TGT_SLOPE);

  HB2("TPCTrk_ResX_Tgt_vs_XTgt;BcOut X_{tgt} [mm];TPC Track - BcOut [mm]", TPC_BINS_RES_TGT_2D_POS, TPC_BINS_RES_TGT_2D_X);
  HB2("TPCTrk_ResY_Tgt_vs_YTgt;BcOut Y_{tgt} [mm];TPC Track - BcOut [mm]", TPC_BINS_RES_TGT_2D_POS, TPC_BINS_RES_TGT_2D_Y);
  HB2("TPCTrk_ResU_Tgt_vs_UTgt;BcOut U_{tgt} (dY/dX);TPC Track - BcOut", TPC_BINS_EXT_TGT_SLOPE_U, TPC_BINS_RES_TGT_2D_UV);
  HB2("TPCTrk_ResV_Tgt_vs_VTgt;BcOut V_{tgt} (dZ/dX);TPC Track - BcOut", TPC_BINS_EXT_TGT_SLOPE_V, TPC_BINS_RES_TGT_2D_UV);

  HB2("TPCTrk_ResY_vs_Layer;Layer;Residual Y (TPC Hit - TPC Track) [mm]", TPC_BINS_LAYER, TPC_BINS_RES_MED);
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    HB1(Form("TPC_PullX_Layer%02d;Pull;Counts", layer), TPC_BINS_PULL);
    HB1(Form("TPC_PullY_Layer%02d;Pull;Counts", layer), TPC_BINS_PULL);
    HB1(Form("TPC_PullZ_Layer%02d;Pull;Counts", layer), TPC_BINS_PULL);
    HB1(Form("TPC_PullLocalX_Layer%02d;Pull;Counts", layer), TPC_BINS_PULL);
    HB1(Form("TPC_PullLocalY_Layer%02d;Pull;Counts", layer), TPC_BINS_PULL);

    HB1(Form("TPCTrk_ResX_Layer%02d;Residual X (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_MED);
    HB1(Form("TPCTrk_ResY_Layer%02d;Residual Y (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_MED);
    HB1(Form("TPCTrk_ResZ_Layer%02d;Residual Z (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_MED);
    HB1(Form("TPCTrk_ResLocalX_Layer%02d;Residual local X (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_LOCAL_X);
    HB1(Form("TPCTrk_ResLocalY_Layer%02d;Residual local Y (TPC Hit - TPC Track) [mm];Counts", layer), TPC_BINS_RES_MED);
    HB1(Form("TPCTrk_ResXY_Layer%02d;|TPC Hit - TPC Track (XY)| [mm];Counts", layer), TPC_BINS_RES_XY_NORM);

    HB1(Form("TPCCl_ResX_Layer%02d;Residual X (TPC Cluster - BcOut) [mm];Counts", layer), TPC_BINS_RES_MED);
    HB1(Form("TPCCl_ResY_Layer%02d;Residual Y (TPC Cluster - BcOut) [mm];Counts", layer), TPC_BINS_RES_MED);
    HB2(Form("TPCCl_ResX_vs_X_Layer%02d;X_{global} [mm];Residual X (TPC Cluster - BcOut) [mm]", layer), TPC_BINS_POS_X, TPC_BINS_RES_HIT);
    HB2(Form("TPCCl_ResY_vs_Y_TPC_Layer%02d;Y (TPC Track) [mm];Residual Y (TPC Cluster - BcOut) [mm]", layer), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
    HB2(Form("TPCCl_ResY_vs_Y_BcOut_Layer%02d;Y (BcOut Track) [mm];Residual Y (TPC Cluster - BcOut) [mm]", layer), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
    if (calib_flag) {
      const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
      for (Int_t row = 0; row < n_pad; ++row) {
        HB2(Form("TPCCl_ResY_vs_Y_TPC_Layer%02d_Row%03d;Y (TPC Track) [mm];Residual Y (TPC Cluster - BcOut) [mm]", layer, row), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
        HB2(Form("TPCCl_ResY_vs_Y_BcOut_Layer%02d_Row%03d;Y (BcOut Track) [mm];Residual Y (TPC Cluster - BcOut) [mm]", layer, row), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
      }
    }
  }
  HB2("TPCCl_ResX_vs_Layer;Layer;Residual X (TPC Cluster - BcOut) [mm]", TPC_BINS_LAYER, TPC_BINS_RES_MED);
  HB2("TPCCl_ResY_vs_Layer;Layer;Residual Y (TPC Cluster - BcOut) [mm]", TPC_BINS_LAYER, TPC_BINS_RES_MED);

  BuildCoBoClockTime(kCoBoClockTime_Track);
}

void
BuildTPCHitBcOutTracking(Bool_t calib_flag)
{
  HB1("TPCHit_ResX;Residual X (TPC Hit - BcOut) [mm];Counts", TPC_BINS_RES_HIT);
  HB1("TPCHit_ResY;Residual Y (TPC Hit - BcOut) [mm];Counts", TPC_BINS_RES_HIT);
  HB1("TPCCl_ResX;Residual X (TPC Cluster - BcOut) [mm];Counts", TPC_BINS_RES_HIT);
  HB1("TPCCl_ResY;Residual Y (TPC Cluster - BcOut) [mm];Counts", TPC_BINS_RES_HIT);

  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    HB1(Form("TPCHit_ResX_Layer%02d;Residual X (TPC Hit - BcOut) [mm];Counts", layer), TPC_BINS_RES_HIT);
    HB1(Form("TPCHit_ResY_Layer%02d;Residual Y (TPC Hit - BcOut) [mm];Counts", layer), TPC_BINS_RES_HIT);

    HB2(Form("TPCHit_ResX_vs_X_Layer%02d;X_{BcOut} (global) [mm];Residual X (TPC Hit - BcOut) [mm]", layer), TPC_BINS_POS_X, TPC_BINS_RES_HIT);
    HB2(Form("TPCHit_ResY_vs_Y_Layer%02d;Y_{BcOut} (global) [mm];Residual Y (TPC Hit - BcOut) [mm]", layer), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);

    HB1(Form("TPCCl_ResX_Layer%02d;Residual X (TPC Cluster - BcOut) [mm];Counts", layer), TPC_BINS_RES_HIT);
    HB1(Form("TPCCl_ResY_Layer%02d;Residual Y (TPC Cluster - BcOut) [mm];Counts", layer), TPC_BINS_RES_HIT);

    HB2(Form("TPCCl_ResX_vs_X_Layer%02d;X_{BcOut} (global) [mm];Residual X (TPC Cluster - BcOut) [mm]", layer), TPC_BINS_POS_X, TPC_BINS_RES_HIT);
    HB2(Form("TPCCl_ResY_vs_Y_Layer%02d;Y_{BcOut} (global) [mm];Residual Y (TPC Cluster - BcOut) [mm]", layer), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);

    if (calib_flag) {
      const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
      for (Int_t row = 0; row < n_pad; ++row) {
        HB1(Form("TPCHit_ResY_Layer%02d_Row%03d;Residual Y (TPC Hit - BcOut) [mm];Counts", layer, row), TPC_BINS_RES_HIT);
        HB2(Form("TPCHit_ResY_vs_Y_Layer%02d_Row%03d;Y_{BcOut} (global) [mm];Residual Y (TPC Hit - BcOut) [mm]", layer, row), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);

        HB1(Form("TPCCl_ResY_Layer%02d_Row%03d;Residual Y (TPC Cluster - BcOut) [mm];Counts", layer, row), TPC_BINS_RES_HIT);
        HB2(Form("TPCCl_ResY_vs_Y_Layer%02d_Row%03d;Y_{BcOut} (global) [mm];Residual Y (TPC Cluster - BcOut) [mm]", layer, row), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
      }
    }
  }

  HB2Poly("TPCHit_HitPat", TPC_EVENT_DISPLAY_BINS);
  HB2Poly("TPCCl_HitPat", TPC_EVENT_DISPLAY_BINS);
  tpc::InitializeHistograms("TPCHit_HitPat");
  tpc::InitializeHistograms("TPCCl_HitPat");

  // Layer-dependent
  HB2("TPCHit_ResX_vs_Layer;Layer;Residual X (TPC Hit - BcOut) [mm];", TPC_BINS_LAYER, TPC_BINS_RES_HIT);
  HB2("TPCHit_ResY_vs_Layer;Layer;Residual Y (TPC Hit - BcOut) [mm];", TPC_BINS_LAYER, TPC_BINS_RES_HIT);
  HB2("TPCCl_ResX_vs_Layer;Layer;Residual X (TPC Cluster - BcOut) [mm];", TPC_BINS_LAYER, TPC_BINS_RES_HIT);
  HB2("TPCCl_ResY_vs_Layer;Layer;Residual Y (TPC Cluster - BcOut) [mm];", TPC_BINS_LAYER, TPC_BINS_RES_HIT);

  HB2("TPCHit_Row_vs_Layer;Layer;Row;", TPC_BINS_ROW_VS_LAYER);
  HB2("TPCCl_Row_vs_Layer;Layer;Row;", TPC_BINS_ROW_VS_LAYER);

  // Clock-time dependent
  BuildCoBoClockTime(kCoBoClockTime_Hit | kCoBoClockTime_Cluster);
}

//_____________________________________________________________________________
void
BuildCoBoClockTime(UInt_t flags)
{
  for(Int_t cobo=0; cobo<NumOfSegCOBO; ++cobo){
    if (flags & kCoBoClockTime_Hit) {
      HB2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
      HB2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#ifdef DEBUG_COBO_CLOCK
      HB2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d_NoClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#endif
    }
    if (flags & kCoBoClockTime_Cluster) {
      HB2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
      HB2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d_RawClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#ifdef DEBUG_COBO_CLOCK
      HB2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d_NoClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#endif
    }
    if (flags & kCoBoClockTime_Track) {
      HB2(Form("TPCTrk_ResY_vs_ClockTime_CoBo%d;Clock Time [ns];Residual Y (TPCCl - TPC Track) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
      HB2(Form("TPCTrk_ResY_vs_ClockTime_CoBo%d_RawClock;Clock Time [ns];Residual Y (TPCCl - TPC Track) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#ifdef DEBUG_COBO_CLOCK
      HB2(Form("TPCTrk_ResY_vs_ClockTime_CoBo%d_NoClock;Clock Time [ns];Residual Y (TPCCl - TPC Track) [mm]", cobo), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#endif
    }
  }

  for(Int_t asad=0; asad<NumOfAsadTPC; ++asad){
    if (flags & kCoBoClockTime_Hit) {
      HB2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
      HB2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d_RawClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#ifdef DEBUG_COBO_CLOCK
      HB2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d_NoClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#endif
    }
    if (flags & kCoBoClockTime_Cluster) {
      HB2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
      HB2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d_RawClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#ifdef DEBUG_COBO_CLOCK
      HB2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d_NoClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#endif
    }
    if (flags & kCoBoClockTime_Track) {
      HB2(Form("TPCTrk_ResY_vs_ClockTime_Asad%02d;Clock Time [ns];Residual Y (TPCCl - TPC Track) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
      HB2(Form("TPCTrk_ResY_vs_ClockTime_Asad%02d_RawClock;Clock Time [ns];Residual Y (TPCCl - TPC Track) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#ifdef DEBUG_COBO_CLOCK
      HB2(Form("TPCTrk_ResY_vs_ClockTime_Asad%02d_NoClock;Clock Time [ns];Residual Y (TPCCl - TPC Track) [mm]", asad), TPC_BINS_CLOCK, TPC_BINS_RES_HIT);
#endif
    }
  }
}

//_____________________________________________________________________________
void
BuildTPCHelixTracking(Bool_t calib_flag)
{
  HB1("HoughDist;Hough distance [mm];Counts", TPC_BINS_HOUGH);
  HB1("HoughDistY;Hough distance Y [mm];Counts", TPC_BINS_HOUGH);
  HB1("NTracks_TPC;N_{track};Counts", TPC_BINS_MULT);
  HB1("TPCTrk_Num_TrackHits;N_{hits} (TPC Track);Counts", TPC_BINS_HITS);
  HB1("TPCTrk_Chisqr;#chi^{2} (TPC Track);Counts", TPC_BINS_CHISQR);
  HB1("TPCTrk_Layer;Layer;Counts", TPC_BINS_LAYER);
  HB1("Mom0;p [GeV/c];Counts", TPC_BINS_MOM0);
  HB1("Mom0_Beam;p [GeV/c];Counts", TPC_BINS_MOM0);
  HB1("Mom0_Accidental;p [GeV/c];Counts", TPC_BINS_MOM0);
  HB1("dEdx_PID;dE/dx PID code;Counts", TPC_BINS_PID_CODE);
  HB1("dEdx_PID_Beam;dE/dx PID code;Counts", TPC_BINS_PID_CODE);
  HB1("dEdx_PID_Accidental;dE/dx PID code;Counts", TPC_BINS_PID_CODE);
  HB2("PID_dEdx_vs_Mom;p [GeV/c];dE/dx (a.u.)", TPC_BINS_MOM0, TPC_BINS_DE);
  HB2("PID_dEdx_vs_SignedMom;q#timesp [GeV/c];dE/dx (a.u.)", TPC_BINS_SIGNED_P, TPC_BINS_DE);
  HB2("PID_dEdx_vs_SignedMom_Beam;q#timesp [GeV/c];dE/dx (a.u.)", TPC_BINS_SIGNED_P, TPC_BINS_DE);
  HB2("PID_dEdx_vs_SignedMom_Accidental;q#timesp [GeV/c];dE/dx (a.u.)", TPC_BINS_SIGNED_P, TPC_BINS_DE);
  HB2("PID_dEdx_vs_Mom_pos;p [GeV/c];dE/dx (a.u.)", TPC_BINS_MOM0, TPC_BINS_DE);
  HB2("PID_dEdx_vs_Mom_neg;p [GeV/c];dE/dx (a.u.)", TPC_BINS_MOM0, TPC_BINS_DE);
  HB2("PID_dEdx_vs_Mom_Pi;p [GeV/c];dE/dx (a.u.)", TPC_BINS_MOM0, TPC_BINS_DE);
  HB2("PID_dEdx_vs_Mom_K;p [GeV/c];dE/dx (a.u.)", TPC_BINS_MOM0, TPC_BINS_DE);
  HB2("PID_dEdx_vs_Mom_Proton;p [GeV/c];dE/dx (a.u.)", TPC_BINS_MOM0, TPC_BINS_DE);

  HB2Poly("TPCTrk_HitPat", TPC_EVENT_DISPLAY_BINS);
  tpc::InitializeHistograms("TPCTrk_HitPat");
  HB2("TPCTrk_Row_vs_Layer;Layer;Row", TPC_BINS_ROW_VS_LAYER);
  HB2("TPCTrk_ResY_vs_Layer_Trk;Layer;Residual Y (TPC Cluster - TPC Track) [mm]", TPC_BINS_LAYER, TPC_BINS_RES_HIT);

  HB1("TPCCl_Size;Cluster size;Counts", TPC_BINS_CL_SIZE);
  HB1("TPCCl_dE;Cluster dE;Counts", TPC_BINS_DE);
  HB1("TPCCl_dE_Pion;Cluster dE (pion PID);Counts", TPC_BINS_DE);
  HB2("TPCCl_dE_vs_Layer;Layer;Cluster dE", TPC_BINS_LAYER, TPC_BINS_DE);
  HB2("Transverse_Diffusion;X_{cluster_center}-X_{pad};A/A_{sum}", TPC_BINS_DIFF_PAD, TPC_BINS_CLUSTER_RATIO);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    HB1(Form("TPCHit_HitPat_Layer%02d;Row;Counts", layer), n_pad, -0.5, n_pad - 0.5);
    HB1(Form("TPCCl_Size_Layer%02d;Cluster size;Counts", layer), TPC_BINS_CL_SIZE);
    HB1(Form("TPCCl_dE_Layer%02d;Cluster dE;Counts", layer), TPC_BINS_DE);
    HB1(Form("TPCCl_dE_Pion_Layer%02d;Cluster dE (pion PID);Counts", layer), TPC_BINS_DE);
    HB2(Form("TPCTrk_ResY_vs_Y_Layer%02d;Y_{TPC Track} (local) [mm];Residual Y (TPC Cluster - TPC Track) [mm]", layer), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
    HB2(Form("Transverse_Diffusion_Layer%02d;X_{cluster_center}-X_{pad};A/A_{sum}", layer), TPC_BINS_DIFF_PAD, TPC_BINS_CLUSTER_RATIO);
    if (calib_flag) {
      for (Int_t row = 0; row < n_pad; ++row) {
        HB1(Form("TPCTrk_ResY_Layer%02d_Row%03d;Residual Y (TPC Cluster - TPC Track) [mm];Counts", layer, row), TPC_BINS_RES_HIT);
        HB2(Form("TPCTrk_ResY_vs_Y_Layer%02d_Row%03d;Y_{TPC Track} (local) [mm];Residual Y (TPC Cluster - TPC Track) [mm]", layer, row), TPC_BINS_POS_Y, TPC_BINS_RES_HIT);
        HB1(Form("TPCCl_dE_Layer%02d_Row%03d;Cluster dE;Counts", layer, row), TPC_BINS_DE);
        HB1(Form("TPCCl_dE_Pion_Layer%02d_Row%03d;Cluster dE (pion PID);Counts", layer, row), TPC_BINS_DE);
      }
    }
  }

  if (calib_flag) {
    BuildCoBoClockTime(kCoBoClockTime_Track);
  }
}

//_____________________________________________________________________________
void
BuildTPCHelixLambda()
{
  HB1("Lambda_Mass;M(p#pi^{-}) [GeV/c^{2}];Counts", 500, 1.05, 1.25);
  HB1("Lambda_CloseDist;Closest distance [mm];Counts", 500, 0.0, 50.0);
  HB1("Lambda_VtxX;Vertex X [mm];Counts", TPC_BINS_POS_X);
  HB1("Lambda_VtxY;Vertex Y [mm];Counts", TPC_BINS_POS_Y);
  HB1("Lambda_VtxZ;Vertex Z [mm];Counts", TPC_BINS_POS_X);
  HB1("Lambda_MomX;Momentum X [GeV/c];Counts", TPC_BINS_SIGNED_P);
  HB1("Lambda_MomY;Momentum Y [GeV/c];Counts", TPC_BINS_SIGNED_P);
  HB1("Lambda_MomZ;Momentum Z [GeV/c];Counts", TPC_BINS_SIGNED_P);
  HB1("Lambda_TargetToVtxX;TargetCenter#rightarrowVertex X [mm];Counts", TPC_BINS_POS_X);
  HB1("Lambda_TargetToVtxY;TargetCenter#rightarrowVertex Y [mm];Counts", TPC_BINS_POS_Y);
  HB1("Lambda_TargetToVtxZ;TargetCenter#rightarrowVertex Z [mm];Counts", TPC_BINS_POS_X);
  HB1("Lambda_TargetToVtxDotMom;cos#theta((TargetCenter#rightarrowVertex),P_{#Lambda});Counts", 200, -1.0, 1.0);
}

//_____________________________________________________________________________
void
BuildTPCHelixK0Short()
{
  HB1("K0_Mass;M(#pi^{+}#pi^{-}) [GeV/c^{2}];Counts", 700, 0.15, 0.85);
  HB1("K0_CloseDist;Closest distance [mm];Counts", 500, 0.0, 50.0);
  HB1("K0_VtxX;Vertex X [mm];Counts", TPC_BINS_POS_X);
  HB1("K0_VtxY;Vertex Y [mm];Counts", TPC_BINS_POS_Y);
  HB1("K0_VtxZ;Vertex Z [mm];Counts", TPC_BINS_POS_X);
  HB1("K0_MomX;Momentum X [GeV/c];Counts", TPC_BINS_SIGNED_P);
  HB1("K0_MomY;Momentum Y [GeV/c];Counts", TPC_BINS_SIGNED_P);
  HB1("K0_MomZ;Momentum Z [GeV/c];Counts", TPC_BINS_SIGNED_P);
  HB1("K0_TargetToVtxX;TargetCenter#rightarrowVertex X [mm];Counts", TPC_BINS_POS_X);
  HB1("K0_TargetToVtxY;TargetCenter#rightarrowVertex Y [mm];Counts", TPC_BINS_POS_Y);
  HB1("K0_TargetToVtxZ;TargetCenter#rightarrowVertex Z [mm];Counts", TPC_BINS_POS_X);
  HB1("K0_TargetToVtxDotMom;cos#theta((TargetCenter#rightarrowVertex),P_{K0});Counts", 200, -1.0, 1.0);
}

}
