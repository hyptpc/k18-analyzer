// -*- C++ -*-

#include "HistTools.hh"

#include <TString.h>

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
// HodoHit
const Double_t hr_time_bins[3]      = {5000., -50., 50.};
const Double_t mh_time_bins[3]      = {500., -50., 50.};
const Double_t hr_tot_time_bins[3]  = {1000., 0., 200.};
const Double_t de_bins[3]           = {1000., 0., 10.};
const Double_t npe_bins[3]          = {700., -50., 300.};  // for Cherenkov (BAC, KVC, SAC3)
// TPC
const Double_t tpc_event_display_bins[4] = { -300., 300., -300., 300. };

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
    // KVC: a,b,c,d,S only. COBO: TDC(U) only → dedicated blocks below.
    for(Int_t ihodo=kBH2; ihodo<kNumHodo; ++ihodo){
      if(ihodo == kKVC || ihodo == kCOBO) continue;
      auto name = NameHodo[ihodo].Data();
      const Double_t* tdc_bins;
      if (HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Ftof))
        tdc_bins = hr_tdc_bins_ftof;
      else
        tdc_bins = hr_tdc_bins;
      Double_t nseg = static_cast<Double_t>(NumOfSegHodo[ihodo]);
      const Double_t seg_bins[3] = {nseg, -0.5, nseg-0.5};
      const Double_t mul_bins[3] = {nseg+1, -0.5, nseg+0.5};
      for (const auto& uord : std::vector<TString>{"U", "D"}) {
        auto ud = uord.Data();
        for (Int_t i = 0; i < nseg; ++i) {
          HB1(Form("%s_ADC_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
          HB1(Form("%s_AwT_seg%d%s%s; channel; count", name, i, ud, b),  adc_bins);
          HB1(Form("%s_AwoT_seg%d%s%s; channel; count", name, i, ud, b), adc_bins);
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b),  tdc_bins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
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
      const Char_t*   dex      = is_cherenkov ? "Npe"    : "DeltaE";
      const Double_t* dex_bins = is_cherenkov ? npe_bins : de_bins;
      const Char_t*   dex_axis = is_cherenkov ? "Npe"    : "mip";
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord: std::vector<TString>{"U", "D"} ){
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

  HB2Poly("TPC_HitPat_Noise",    tpc_event_display_bins);
  HB2Poly("TPC_HitPat_Baseline", tpc_event_display_bins);
  tpc::InitializeHistograms("TPC_HitPat_Noise");
  tpc::InitializeHistograms("TPC_HitPat_Baseline");
}

//_____________________________________________________________________________
void
BuildTPCBasic()
{
  // {n_bin, x_min, x_max}
  const Double_t tpc_bins_mult[3]   = {40.,    0.,   40.};   // N_track, N_hits
  const Double_t tpc_bins_hits[3]   = {50.,    0.,   50.};
  const Double_t tpc_bins_chisqr[3] = {500.,   0.,  100.};
  const Double_t tpc_bins_pos[3]    = {400., -100., 100.};   // mm (X0, Y0)
  const Double_t tpc_bins_slope[3]  = {200., -0.20, 0.20};   // dy/dx, dz/dx

  HB1("Num_Track_TPC;N_{track};Counts", tpc_bins_mult);
  HB1("Num_Track_TPC_Hits;N_{hits};Counts", tpc_bins_hits);
  HB1("Chisqr_TPC;#chi^{2};Counts", tpc_bins_chisqr);
  HB1("X0_TPC;X_{0} [mm];Counts", tpc_bins_pos);
  HB1("Y0_TPC;Y_{0} [mm];Counts", tpc_bins_pos);
  HB1("U0_TPC;dY/dX;Counts", tpc_bins_slope);
  HB1("V0_TPC;dZ/dX;Counts", tpc_bins_slope);
}

//_____________________________________________________________________________
void
BuildTPCTracking()
{
  // 1D: {n_bin, x_min, x_max}
  const Double_t trk_bins_hough[3]      = {500.,    0.,   50.};   // mm
  const Double_t trk_bins_time[3]       = {100.,    0.,  100.};   // ms
  const Double_t trk_bins_iter[3]       = {100.,    0.,  100.};
  const Double_t trk_bins_flag[3]       = {10.,     0.,   10.};
  const Double_t trk_bins_minuit[3]     = {5.,      0.,    5.};
  const Double_t trk_bins_pos_wide[3]   = {200., -250.,  250.};   // mm
  const Double_t trk_bins_res_abs[3]    = {200.,    0.,  10.};    // mm |Residual|
  const Double_t trk_bins_res[3]        = {200.,  -2.0,  2.0};    // mm (X,Y,Z)
  const Double_t trk_bins_cl_size[3]    = {25.,     0.,  25.};
  const Double_t trk_bins_de[3]         = {1000.,   0., 2000.};   // keV or ADC
  // 2D: {n_bin_x, x_min, x_max, n_bin_y, y_min, y_max}
  const Double_t trk_bins_2d_pos_slope[6]  = {100., -100., 100., 100., -0.20, 0.20};  // mm, dy/dx
  const Double_t trk_bins_2d_pos_pos[6]    = {100., -100., 100., 100., -100., 100.};  // mm
  const Double_t trk_bins_2d_res_vs_pos[6] = {250., -250., 250., 100.,  -1.0,  1.0};  // mm
  const Double_t trk_bins_2d_pos_wide[6]   = {100., -250., 250., 100., -250., 250.};  // mm
  const Double_t trk_bins_2d_ratio[6]      = {60.,  -15.,  15., 100.,    0.,  1.};   // mm, A/A_sum


  HB1("Hough_Dist;Hough distance [mm];Counts", trk_bins_hough);
  HB1("Hough_Dist_Y;Hough distance Y [mm];Counts", trk_bins_hough);
  HB1("Num_Tracking_Iterations;N_{iter};Counts", trk_bins_iter);
  HB1("Fitting_Flag;Flag;Counts", trk_bins_flag);
  HB1("Track_Searching_Time;Time [ms];Counts", trk_bins_time);
  HB1("Track_Fitting_Time;Time [ms];Counts", trk_bins_time);
  HB1("Minuit_Output_Status;Status;Counts", trk_bins_minuit);

  HB1("Layer_Id_TPC;Layer;Counts", NumOfLayersTPC, -0.5, NumOfLayersTPC - 0.5);
  HB2("U0_vs_X0_TPC;X_{0} [mm];dY/dX", trk_bins_2d_pos_slope);
  HB2("V0_vs_Y0_TPC;Y_{0} [mm];dZ/dX", trk_bins_2d_pos_slope);
  HB2("X0_vs_Y0_TPC;Y_{0} [mm];X_{0} [mm]", trk_bins_2d_pos_pos);

  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    HB1(Form("HitPat_TPC_Layer%02d;Row;Counts", layer), n_pad, -0.5, n_pad - 0.5);
    HB1(Form("Position_TPC_Layer%02d;Position [mm];Counts", layer), trk_bins_pos_wide);
    HB1(Form("Residual_TPC_Layer%02d;|Residual| [mm];Counts", layer), trk_bins_res_abs);
    HB2(Form("Residual_vs_Position_TPC_Layer%02d;Position [mm];Residual [mm]", layer), trk_bins_2d_res_vs_pos);
    HB2(Form("Yhit_vs_Xcal_TPC_Layer%02d;X_{cal} [mm];Y_{hit} [mm]", layer), trk_bins_2d_pos_wide);
    HB1(Form("ResidualX_TPC_Layer%02d;Residual X [mm];Counts", layer), trk_bins_res);
    HB1(Form("ResidualY_TPC_Layer%02d;Residual Y [mm];Counts", layer), trk_bins_res);
    HB1(Form("ResidualZ_TPC_Layer%02d;Residual Z [mm];Counts", layer), trk_bins_res);
  }

  HB1("Cluster_size;Cluster size;Counts", trk_bins_cl_size);
  HB1("Cluster_dE;Cluster dE [keV];Counts", trk_bins_de);
  HB2("Ratio_vs_Dist_Transverse_diffusion;X_{cluster}-X_{pad} [mm];A/A_{sum}", trk_bins_2d_ratio);
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    HB1(Form("Cluster_size_layer%02d;Cluster size;Counts", layer), trk_bins_cl_size);
    HB1(Form("Cluster_dE_layer%02d;Cluster dE [keV];Counts", layer), trk_bins_de);
    HB2(Form("Ratio_vs_Dist_Transverse_diffusion_Layer%02d;X_{cluster}-X_{pad} [mm];A/A_{sum}", layer), trk_bins_2d_ratio);
  }

  // TH2Poly: track-hit pad occupancy for event-display
  HB2Poly("TPC_TrackHitPat", tpc_event_display_bins);
  tpc::InitializeHistograms("TPC_TrackHitPat");
}

//_____________________________________________________________________________
void
BuildTPCBcOutTracking()
{
  // 1D: {n_bin, x_min, x_max}
  const Double_t bc_bins_mult[3]       = {40.,    0.,   40.};
  const Double_t bc_bins_hits[3]       = {50.,    0.,   50.};
  const Double_t bc_bins_chisqr[3]     = {500.,   0.,  500.};
  const Double_t bc_bins_pos_1d[3]     = {400., -100., 100.};   // mm (X0,Y0,Xtgt,Ytgt)
  const Double_t bc_bins_slope_1d[3]   = {200., -0.20, 0.20};   // dY/dX, dZ/dX
  const Double_t bc_bins_diff_xy[3]    = {640.,  -16.,  16.};   // BcOut-Tpc [mm]
  const Double_t bc_bins_diff_slope[3] = {200., -0.05, 0.05};   // BcOut-Tpc
  const Double_t bc_bins_pull[3]       = {200.,  -5.,   5.};
  const Double_t bc_bins_res[3]        = {200.,  -8.,   8.};    // mm
  const Double_t bc_bins_res_local_x[3] = {400., -16.,  16.};   // mm
  const Double_t bc_bins_xz[3]         = {100.,   0.,   8.};    // mm |Residual XZ|

  // 2D: {n_bin_x, x_min, x_max, n_bin_y, y_min, y_max}
  const Double_t bc_bins_2d_pos_slope[6]  = {100., -100., 100., 100., -0.20, 0.20};  // mm, dY/dX
  const Double_t bc_bins_2d_pos_pos[6]    = {100., -100., 100., 100., -100., 100.};  // mm
  const Double_t bc_bins_2d_atan[6]       = {300., -300., 300., 100.,  -20.,  20.};  // mm, atan [deg]
  const Double_t bc_bins_2d_xtgt[6]       = {400., -200., 200., 400., -200., 200.};  // Xtgt BcOut vs Tpc
  const Double_t bc_bins_2d_ytgt[6]       = {400., -100., 100., 400., -100., 100.};  // Ytgt BcOut vs Tpc
  const Double_t bc_bins_2d_utgt[6]       = {400., -0.15, 0.15, 400., -0.15, 0.15};
  const Double_t bc_bins_2d_vtgt[6]       = {400., -0.05, 0.05, 400., -0.05, 0.05};
  const Double_t bc_bins_2d_diff_xtgt[6]  = {400., -150., 150., 400.,  -10.,  10.};  // BcOut-Tpc [mm]
  const Double_t bc_bins_2d_diff_ytgt[6]  = {400., -150., 150., 400.,  -15.,  15.};
  const Double_t bc_bins_2d_diff_utgt[6]  = {400., -0.15, 0.15, 400., -0.03, 0.03};
  const Double_t bc_bins_2d_diff_vtgt[6]  = {400., -0.05, 0.05, 400., -0.03, 0.03};
  // Parameter tuning: Position correction (BcOut reference)
  const Double_t bc_bins_2d_resx_vs_x[6]   = {200., -100., 100., 100., -20.0, 20.0};  // X (global) [mm], Residual X [mm]
  const Double_t bc_bins_2d_resy_vs_y[6]   = {200., -100., 100., 100., -20.0, 20.0};  // Y (TPC/BcOut Tracking) [mm], Residual Y (BcOut) [mm]
  const Double_t bc_bins_2d_layer_res_y[6] = {
    static_cast<Double_t>(NumOfLayersTPC), -0.5, static_cast<Double_t>(NumOfLayersTPC) - 0.5,
    400., -5., 5.
  };  // Layer, Residual Y [mm] (Layer_vs_ResY, TPC_Residual*_BcOut_vs_Layer)
  const Double_t bc_bins_2d_clocktime_resy[6] = { 200., -60., 50., 200., -20., 20. }; // Clock Time [ns], Residual Y (BcOut) [mm]

  HB2("X0_vs_U0_TPC;X_{0} [mm];dY/dX", bc_bins_2d_pos_slope);
  HB2("Y0_vs_V0_TPC;Y_{0} [mm];dZ/dX", bc_bins_2d_pos_slope);
  HB2("X0_vs_Y0_TPC;X_{0} [mm];Y_{0} [mm]", bc_bins_2d_pos_pos);

  HB2("X0_vs_atanU0_TPC;X_{0} [mm];atan(dY/dX) [deg]", bc_bins_2d_atan);
  HB2("Y0_vs_atanV0_TPC;Y_{0} [mm];atan(dZ/dX) [deg]", bc_bins_2d_atan);

  HB1("Num_Track_BcOut;N_{track};Counts", bc_bins_mult);
  HB1("Num_Track_BcOut_Hits;N_{hits};Counts", bc_bins_hits);
  HB1("Chisqr_BcOut;#chi^{2};Counts", bc_bins_chisqr);
  HB1("X0_BcOut;X_{0} [mm];Counts", bc_bins_pos_1d);
  HB1("Y0_BcOut;Y_{0} [mm];Counts", bc_bins_pos_1d);
  HB1("U0_BcOut;dY/dX;Counts", bc_bins_slope_1d);
  HB1("V0_BcOut;dZ/dX;Counts", bc_bins_slope_1d);
  HB1("Xtgt_BcOut;X_{tgt} [mm];Counts", bc_bins_pos_1d);
  HB1("Ytgt_BcOut;Y_{tgt} [mm];Counts", bc_bins_pos_1d);
  HB1("Utgt_BcOut;dY/dX;Counts", bc_bins_slope_1d);
  HB1("Vtgt_BcOut;dZ/dX;Counts", bc_bins_slope_1d);
  HB2("Xtgt_vs_Utgt_BcOut;X_{tgt} [mm];dY/dX", bc_bins_2d_pos_slope);
  HB2("Ytgt_vs_Vtgt_BcOut;Y_{tgt} [mm];dZ/dX", bc_bins_2d_pos_slope);
  HB2("Xtgt_vs_Ytgt_BcOut;X_{tgt} [mm];Y_{tgt} [mm]", bc_bins_2d_pos_pos);

  HB2("Xtgt_BcOut_vs_Tpc;Tpc X_{tgt} [mm];BcOut X_{tgt} [mm]", bc_bins_2d_xtgt);
  HB2("Ytgt_BcOut_vs_Tpc;Tpc Y_{tgt} [mm];BcOut Y_{tgt} [mm]", bc_bins_2d_ytgt);
  HB2("Utgt_BcOut_vs_Tpc;Tpc dY/dX;BcOut dY/dX", bc_bins_2d_utgt);
  HB2("Vtgt_BcOut_vs_Tpc;Tpc dZ/dX;BcOut dZ/dX", bc_bins_2d_vtgt);
  HB1("Xtgt_Diff;BcOut-Tpc [mm];Counts", bc_bins_diff_xy);
  HB1("Ytgt_Diff;BcOut-Tpc [mm];Counts", bc_bins_diff_xy);
  HB1("Utgt_Diff;BcOut-Tpc;Counts", bc_bins_diff_slope);
  HB1("Vtgt_Diff;BcOut-Tpc;Counts", bc_bins_diff_slope);

  HB2("Xtgt_Diff_vs_Xtgt_BcOut;BcOut X_{tgt} [mm];BcOut-Tpc [mm]", bc_bins_2d_diff_xtgt);
  HB2("Ytgt_Diff_vs_Ytgt_BcOut;BcOut Y_{tgt} [mm];BcOut-Tpc [mm]", bc_bins_2d_diff_ytgt);
  HB2("Utgt_Diff_vs_Utgt_BcOut;BcOut dY/dX;BcOut-Tpc", bc_bins_2d_diff_utgt);
  HB2("Vtgt_Diff_vs_Vtgt_BcOut;BcOut dZ/dX;BcOut-Tpc", bc_bins_2d_diff_vtgt);

  HB2("Layer_vs_ResY;Layer;Y Residual [mm]", bc_bins_2d_layer_res_y);
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    HB1(Form("TPC_Layer%02d_X_Pull;Pull;Counts", layer), bc_bins_pull);
    HB1(Form("TPC_Layer%02d_Y_Pull;Pull;Counts", layer), bc_bins_pull);
    HB1(Form("TPC_Layer%02d_Z_Pull;Pull;Counts", layer), bc_bins_pull);
    HB1(Form("TPC_Layer%02d_Local_X_Pull;Pull;Counts", layer), bc_bins_pull);
    HB1(Form("TPC_Layer%02d_Local_Y_Pull;Pull;Counts", layer), bc_bins_pull);

    HB1(Form("TPC_Layer%02d_X_Residual;Residual X [mm];Counts", layer), bc_bins_res);
    HB1(Form("TPC_Layer%02d_Y_Residual;Residual Y [mm];Counts", layer), bc_bins_res);
    HB1(Form("TPC_Layer%02d_Z_Residual;Residual Z [mm];Counts", layer), bc_bins_res);
    HB1(Form("TPC_Layer%02d_Local_X_Residual;Residual local X [mm];Counts", layer), bc_bins_res_local_x);
    HB1(Form("TPC_Layer%02d_Local_Y_Residual;Residual local Y [mm];Counts", layer), bc_bins_res);
    HB1(Form("TPC_Layer%02d_XZ_Residual;|Residual XZ| [mm];Counts", layer), bc_bins_xz);

    HB1(Form("TPC_Layer%02d_BcOut_X_Residual;Residual X [mm];Counts", layer), bc_bins_res);
    HB1(Form("TPC_Layer%02d_BcOut_Y_Residual;Residual Y [mm];Counts", layer), bc_bins_res);
    // Parameter tuning: Position correction (BcOut reference)
    HB2(Form("TPC_ResidualX_vs_X_BcOut_Layer%02d;X (global) [mm];Residual X (BcOut) [mm]", layer), bc_bins_2d_resx_vs_x);
    HB2(Form("TPC_ResidualY_vs_Y_TPC_Layer%02d;Y (TPC Tracking) [mm];Residual Y (BcOut) [mm]", layer), bc_bins_2d_resy_vs_y);
    HB2(Form("TPC_ResidualY_vs_Y_BcOut_Layer%02d;Y (BcOut Tracking) [mm];Residual Y (BcOut) [mm]", layer), bc_bins_2d_resy_vs_y);
    // Parameter tuning: Row-dependent residual (BcOut reference)
    const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < n_pad; ++row) {
      HB2(Form("TPC_ResidualY_vs_Y_TPC_Layer%02d_Row%03d;Y (TPC Tracking) [mm];Residual Y (BcOut) [mm]", layer, row), bc_bins_2d_resy_vs_y);
      HB2(Form("TPC_ResidualY_vs_Y_BcOut_Layer%02d_Row%03d;Y (BcOut Tracking) [mm];Residual Y (BcOut) [mm]", layer, row), bc_bins_2d_resy_vs_y);
    }
  }
  
  // Parameter tuning: Layer-dependent residual distribution (BcOut reference)
  HB2("TPC_ResidualX_BcOut_vs_Layer;Layer;Residual X (BcOut) [mm]", bc_bins_2d_layer_res_y);
  HB2("TPC_ResidualY_BcOut_vs_Layer;Layer;Residual Y (BcOut) [mm]", bc_bins_2d_layer_res_y);

  // Parameter tuning: Clock time vs BcOut Residual Y (CoBo and Asad)
  for(Int_t c=0; c<NumOfSegCOBO; ++c){
    HB2(Form("TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d;Clock Time [ns];Residual Y (BcOut) [mm]", c), bc_bins_2d_clocktime_resy);
    HB2(Form("TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock;Clock Time [ns];Residual Y (BcOut) [mm]", c), bc_bins_2d_clocktime_resy);
#ifdef DEBUG_COBO_CLOCK
    HB2(Form("TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_NoClock;Clock Time [ns];Residual Y (BcOut) [mm]", c), bc_bins_2d_clocktime_resy);
#endif
  }
  for(Int_t a=0; a<NumOfAsadTPC; ++a){
    HB2(Form("TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d;Clock Time [ns];Residual Y (BcOut) [mm]", a), bc_bins_2d_clocktime_resy);
    HB2(Form("TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d_RawClock;Clock Time [ns];Residual Y (BcOut) [mm]", a), bc_bins_2d_clocktime_resy);
#ifdef DEBUG_COBO_CLOCK
    HB2(Form("TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d_NoClock;Clock Time [ns];Residual Y (BcOut) [mm]", a), bc_bins_2d_clocktime_resy);
#endif
  }
}

void
BuildTPCHitBcOutTracking()
{
  const Double_t residual_bins[3]       = {200.,  -20.,  20.};    // mm
  const Double_t residual_bins_2d[6]    = {200., -100., 100., 200., -20.0, 20.0};
  const Double_t residual_layer_bins[6] = {
    static_cast<Double_t>(NumOfLayersTPC), -0.5, static_cast<Double_t>(NumOfLayersTPC) - 0.5,
    200., -20., 20.
  };
  const Double_t residual_clocktime_bins[6] = { 200., -50., 50., 200., -20., 20. };
  const Double_t row_layer_bins[6] = {
    static_cast<Double_t>(NumOfLayersTPC), -0.5, static_cast<Double_t>(NumOfLayersTPC) - 0.5,
    244., -0.5, 243.5
  };

  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {
    // TPC Hit Residuals
    HB1(Form("TPCHit_ResX_Layer%02d;Residual X (TPC Hit - BcOut) [mm];", layer), residual_bins);
    HB1(Form("TPCHit_ResY_Layer%02d;Residual Y (TPC Hit - BcOut) [mm];", layer), residual_bins);

    HB2(Form("TPCHit_ResX_vs_X_Layer%02d;X_{BcOut} (global) [mm];Residual X (TPC Hit - BcOut) [mm];", layer), residual_bins_2d);
    HB2(Form("TPCHit_ResY_vs_Y_Layer%02d;Y_{BcOut} (global) [mm];Residual Y (TPC Hit - BcOut) [mm];", layer), residual_bins_2d);

    // TPC Cluster Residuals
    HB1(Form("TPCCl_ResX_Layer%02d;Residual X (TPC Cluster - BcOut) [mm];", layer), residual_bins);
    HB1(Form("TPCCl_ResY_Layer%02d;Residual Y (TPC Cluster - BcOut) [mm];", layer), residual_bins);

    HB2(Form("TPCCl_ResX_vs_X_Layer%02d;X_{BcOut} (global) [mm];Residual X (TPC Cluster - BcOut) [mm];", layer), residual_bins_2d);
    HB2(Form("TPCCl_ResY_vs_Y_Layer%02d;Y_{BcOut} (global) [mm];Residual Y (TPC Cluster - BcOut) [mm];", layer), residual_bins_2d);

    const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < n_pad; ++row) {
      HB2(Form("TPCHit_ResY_vs_Y_Layer%02d_Row%03d;Y_{BcOut} (global) [mm];Residual Y (TPC Hit - BcOut) [mm];", layer, row), residual_bins_2d);
      HB2(Form("TPCCl_ResY_vs_Y_Layer%02d_Row%03d;Y_{BcOut} (global) [mm];Residual Y (TPC Cluster - BcOut) [mm];", layer, row), residual_bins_2d);
    }
  }

  HB2Poly("TPC_HitPat",         tpc_event_display_bins);
  HB2Poly("TPC_Cluster_HitPat", tpc_event_display_bins);
  tpc::InitializeHistograms("TPC_HitPat");
  tpc::InitializeHistograms("TPC_Cluster_HitPat");

  // Layer-dependent
  HB2("TPCHit_ResX_vs_Layer;Layer;Residual X (TPC Hit - BcOut) [mm];", residual_layer_bins);
  HB2("TPCHit_ResY_vs_Layer;Layer;Residual Y (TPC Hit - BcOut) [mm];", residual_layer_bins);
  HB2("TPCCl_ResX_vs_Layer;Layer;Residual X (TPC Cluster - BcOut) [mm];", residual_layer_bins);
  HB2("TPCCl_ResY_vs_Layer;Layer;Residual Y (TPC Cluster - BcOut) [mm];", residual_layer_bins);

  HB2("TPCHit_Row_vs_Layer;Layer;Row;", row_layer_bins);
  HB2("TPCCl_Row_vs_Layer;Layer;Row;", row_layer_bins);

  // Clock-time dependent
  for(Int_t cobo=0; cobo<NumOfSegCOBO; ++cobo){
    HB2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm];", cobo), residual_clocktime_bins);
    HB2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm];", cobo), residual_clocktime_bins);

    HB2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm];", cobo), residual_clocktime_bins);
    HB2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d_RawClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm];", cobo), residual_clocktime_bins);

#ifdef DEBUG_COBO_CLOCK
    HB2(Form("TPCHit_ResY_vs_ClockTime_CoBo%d_NoClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm];", cobo), residual_clocktime_bins);
    HB2(Form("TPCCl_ResY_vs_ClockTime_CoBo%d_NoClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm];", cobo), residual_clocktime_bins);
#endif
  }
  for(Int_t asad=0; asad<NumOfAsadTPC; ++asad){
    HB2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm];", asad), residual_clocktime_bins);
    HB2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d_RawClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm];", asad), residual_clocktime_bins);

    HB2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm];", asad), residual_clocktime_bins);
    HB2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d_RawClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm];", asad), residual_clocktime_bins);

#ifdef DEBUG_COBO_CLOCK
    HB2(Form("TPCHit_ResY_vs_ClockTime_Asad%02d_NoClock;Clock Time [ns];Residual Y (TPC Hit - BcOut) [mm];", asad), residual_clocktime_bins);
    HB2(Form("TPCCl_ResY_vs_ClockTime_Asad%02d_NoClock;Clock Time [ns];Residual Y (TPC Cluster - BcOut) [mm];", asad), residual_clocktime_bins);
#endif
  }
}

void
BuildTPCHelixTracking()
{
  HB1("HoughDist", 500, 0., 50.);
  HB1("HoughDistY", 500, 0., 50.);
  HB1("NTracks_TPC", 40, 0., 40. );
  HB1("NHits_Track_TPC", 50, 0., 50.);
  HB1("Chisqr_TPC", 500, 0., 500.);
  HB1("LayerId_TPC", 35, 0., 35.);
  HB1("mom0", 1000, 0., 2.5);

  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 2000.;

  const Int_t NbinClSize = 25;
  const Double_t MinClSize = 0;
  const Double_t MaxClSize = 25;
  const Int_t NbinDist = 60;
  const Double_t MinDist = -15.;
  const Double_t MaxDist = 15.;
  const Int_t NbinRatio = 100;
  const Double_t MinRatio = 0.;
  const Double_t MaxRatio = 1.;

  HB1("Cluster_size;Cluster size;Counts", NbinClSize, MinClSize, MaxClSize);
  HB1("Cluster_dE;Cluster dE;Counts", NbinDe, MinDe, MaxDe);
  HB2("Transverse_diffusion;X_{cluster_center}-X_{pad};A/A_{sum}", NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(Form("Cluster_size_layer%d;Cluster size;Counts",layer), NbinClSize, MinClSize, MaxClSize);
    HB1(Form("Cluster_dE_layer%d;Cluster dE;Counts",layer), NbinDe, MinDe, MaxDe);
    HB2(Form("Transverse_diffusion_layer%d;X_{cluster_center}-X_{pad};A/A_{sum}",layer), NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  }
}

}
