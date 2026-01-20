// -*- C++ -*-

#include "HistTools.hh"

#include <TString.h>

#include <DAQNode.hh>
#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#include "DetectorID.hh"
#include "RootHelper.hh"
#include "TPCPadHelper.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
//using root::HB1;
//using root::HB2;
using namespace root; 
}

namespace hist
{
// Raw
const Double_t hrtdcbins1[3] = {20000,  600000,  800000};
const Double_t hrtdcbins2[3] = {20000,  200000,  600000}; // for FTOF
const Double_t hrtdcbins3[3] = {20000, 1200000, 1600000}; // for BHT
const Double_t hrtdcbins4[3] = {10000, 1500000, 1650000}; // for COBO
const Double_t hrtdcbins5[3] = {10000,       0, 2000000}; // TriggerFlag
const Double_t hrtotbins[3]  = {5000, 0, 50000};
const Double_t adcbins[3]    = {4096, -0.5, 4095.5};
const Double_t mhtdcbins[3]  = {2000, 0, 2000};
const Double_t mhtotbins[3]  = {1000, 0, 1000};
// HodoHit
const Double_t hrtimebins[3]    = {5000, -50, 50};
const Double_t mhtimebins[3]    = {500, -50, 50};
const Double_t hrtottimebins[3] = {1000, 0, 200};
const Double_t debins[3]        = {1000, 0, 10};
const Double_t npebins[3]       = {700, -50.0, 300.};  // for Cherenkov (BAC, KVC, SAC3)

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
    HB1(Form("%s_TDC_seg%d", name, i), hrtdcbins5);
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
      Int_t nseg = NumOfSegBHT;
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"} ){
          const Char_t* ud = uord.Data();
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b), hrtdcbins3);
          HB1(Form("%s_Trailing_seg%d%s%s; channel; count", name, i, ud, b), hrtdcbins3);
          HB1(Form("%s_TOT_seg%d%s%s; channel; count", name, i, ud, b), hrtotbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s; segment; count", name, ud, b), nseg, -0.5, nseg - 0.5);
        HB1(Form("%s_Multi_%s%s; multiplicity; count", name, ud, b), nseg + 1, -0.5, nseg + 0.5);
      }
    }

    { // BAC
      const Char_t* name = "BAC";
      Int_t nseg = NumOfSegBAC;
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_ADC_seg%d%s; channel; count", name, i, b), adcbins);
        HB1(Form("%s_AwT_seg%d%s; channel; count", name, i, b), adcbins);
        HB1(Form("%s_AwoT_seg%d%s; channel; count", name, i, b), adcbins);
        HB1(Form("%s_TDC_seg%d%s; channel; count", name, i, b), hrtdcbins1);
      }
      HB1(Form("%s_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }
    
    ///// BHT-BAC
    {
      HB2(Form("BAC_ADC_vs_BHT_TDC%s", b),
          200, 720000., 750000., 200, 0., 2000.);
    }

    // Hodoscope (U/D or 1ch). KVC: a,b,c,d,S only. COBO: TDC(U) only → dedicated blocks below.
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      if(ihodo == kKVC || ihodo == kCOBO) continue;
      auto name = NameHodo[ihodo].Data();
      const Double_t* hrtdcbins = (ihodo == kCVC || ihodo == kSFV || ihodo == kSAC3) ? hrtdcbins2 : hrtdcbins1;
      Int_t nseg = NumOfSegHodo[ihodo];
      for(const auto& uord: std::vector<TString>{"U", "D"}){
        auto ud = uord.Data();
        for(Int_t i=0; i<nseg; ++i){
          HB1(Form("%s_ADC_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_AwT_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_AwoT_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b), hrtdcbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s; segment; count", name, ud, b), nseg, -0.5, nseg - 0.5);
        HB1(Form("%s_Multi_%s%s; multiplicity; count", name, ud, b), nseg + 1, -0.5, nseg + 0.5);
      }
    }

    { // HTOF Sum (S): U/D is in the loop above
      auto name = "HTOF";
      const Double_t* hrtdcbins = hrtdcbins1;
      Int_t nseg = NumOfSegHTOF;
      for(const auto& uord: std::vector<TString>{"S"}){
        auto ud = uord.Data();
        for(Int_t i=0; i<nseg; ++i){
          HB1(Form("%s_ADC_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_AwT_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_AwoT_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b), hrtdcbins);
        }
      }
      HB1(Form("%s_HitPat_HT%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Multi_HT%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }

    { // KVC: a,b,c,d,S (no U/D). HitPat/Multi from OR/AND.
      auto name = "KVC";
      const Double_t* hrtdcbins = hrtdcbins1;
      Int_t nseg = NumOfSegKVC;
      for(const auto& uord: std::vector<TString>{"a", "b", "c", "d", "S"}){
        auto ud = uord.Data();
        for(Int_t i=0; i<nseg; ++i){
          HB1(Form("%s_ADC_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_AwT_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_AwoT_seg%d%s%s; channel; count", name, i, ud, b), adcbins);
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b), hrtdcbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"}){
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s; segment; count", name, ud, b), nseg, -0.5, nseg - 0.5);
        HB1(Form("%s_Multi_%s%s; multiplicity; count", name, ud, b), nseg + 1, -0.5, nseg + 0.5);
      }
    }

    { // COBO: TDC(U) only
      auto name = "COBO";
      Int_t nseg = NumOfSegHodo[kCOBO];
      for(Int_t i=0; i<nseg; ++i)
        HB1(Form("%s_TDC_seg%dU%s; channel; count", name, i, b), hrtdcbins4);
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
      Double_t nseg = NumOfSegBHT;
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB1(Form("%s_Hit_Time_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
          HB1(Form("%s_Hit_CTime_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
          HB1(Form("%s_Hit_TOT_seg%d%s%s; ns; count", name, i, ud, b), hrtottimebins);
          HB1(Form("%s_Hit_DeltaE_seg%d%s%s; mip; count", name, i, ud, b), debins);
        }
        HB1(Form("%s_Hit_MeanTime_seg%d%s; ns; count", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CMeanTime_seg%d%s; ns; count", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_MeanTOT_seg%d%s; ns; count", name, i, b), hrtottimebins);
        HB1(Form("%s_Hit_DeltaE_seg%d%s; mip; count", name, i, b), debins);
      }
      HB1(Form("%s_Hit_MeanTime%s; ns; count", name, b), hrtimebins);
      HB1(Form("%s_Hit_CMeanTime%s; ns; count", name, b), hrtimebins);
      HB1(Form("%s_Hit_MeanTOT%s; ns; count", name, b), hrtottimebins);
      HB1(Form("%s_Hit_DeltaE%s; mip; count", name, b), debins);
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t hrtottimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtottimebins[0]/5, hrtottimebins[1], hrtottimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      HB2(Form("%s_Hit_MeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_CMeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_MeanTOT_vs_HitPat%s; segment; ns", name, b), hrtottimebins2d);
      HB2(Form("%s_Hit_DeltaE_vs_HitPat%s; segment; mip", name, b), debins2d);
      HB1(Form("%s_Hit_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Hit_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }

    // Hodoscope (U/D). BAC and KVC use only Sum (S) → "BAC, HTOF, KVC Sum" block below.
    // Cherenkov (SAC3): Hit_Npe / npebins. Others: Hit_DeltaE / debins.
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      if(ihodo == kBAC || ihodo == kKVC) continue;
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = NumOfSegHodo[ihodo];
      Bool_t is_cherenkov = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Cherenkov);
      const Char_t* dex = is_cherenkov ? "Npe" : "DeltaE";
      const Double_t* dex_bins = is_cherenkov ? npebins : debins;
      const Char_t* dex_axis = is_cherenkov ? "Npe" : "mip";
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord: std::vector<TString>{"U", "D"} ){
          auto ud = uord.Data();
          HB1(Form("%s_Hit_%s_seg%d%s%s; %s; count", name, dex, i, ud, b, dex_axis), dex_bins);
          HB1(Form("%s_Hit_Time_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
          HB1(Form("%s_Hit_CTime_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
        }
        HB1(Form("%s_Hit_%s_seg%d%s; %s; count", name, dex, i, b, dex_axis), dex_bins);
        HB1(Form("%s_Hit_MeanTime_seg%d%s; ns; count", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CMeanTime_seg%d%s; ns; count", name, i, b), hrtimebins);
      }
      HB1(Form("%s_Hit_MeanTime%s; ns; count", name, b), hrtimebins);
      HB1(Form("%s_Hit_CMeanTime%s; ns; count", name, b), hrtimebins);
      HB1(Form("%s_Hit_%s%s; %s; count", name, dex, b, dex_axis), is_cherenkov ? npebins : hrtottimebins);
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      const Double_t npebins2d[6] = { nseg, -0.5, nseg - 0.5,
        npebins[0]/10, npebins[1], npebins[2] };
      HB2(Form("%s_Hit_MeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_CMeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_%s_vs_HitPat%s; segment; %s", name, dex, b, dex_axis), is_cherenkov ? npebins2d : debins2d);
      HB1(Form("%s_Hit_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Hit_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }

    // BAC Sum: S_online (seg4), offline_sum (event). Time/CTime_seg S. No a,b,c,d.
    {
      auto name = "BAC";
      Double_t nseg = NumOfSegBAC;
      HB1(Form("%s_Hit_Npe_offline_sum%s; Npe; count", name, b), npebins);
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_Hit_Npe_seg%dS_online%s; Npe; count", name, i, b), npebins);
        HB1(Form("%s_Hit_Time_seg%dS%s; ns; count", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CTime_seg%dS%s; ns; count", name, i, b), hrtimebins);
      }
      HB1(Form("%sSum_Hit_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%sSum_Hit_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }

    // HTOF Sum (S, kExtra): Hit_DeltaE, Time/CTime_seg S, Sum_Hit_HitPat/Multi
    {
      auto name = "HTOF";
      Double_t nseg = NumOfSegHTOF;
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_Hit_DeltaE_seg%dS%s; mip; count", name, i, b), debins);
        HB1(Form("%s_Hit_Time_seg%dS%s; ns; count", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CTime_seg%dS%s; ns; count", name, i, b), hrtimebins);
      }
      HB1(Form("%sSum_Hit_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%sSum_Hit_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }

    // KVC Sum (kExtra): Hit_Npe_seg S_offline (raw-based per seg), S_online + a,b,c,d. Time/CTime_seg S, Sum_Hit_HitPat/Multi
    {
      auto name = "KVC";
      Double_t nseg = NumOfSegKVC;
      const Char_t* abcd[4] = {"a", "b", "c", "d"};
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_Hit_Npe_seg%dS_offline%s; Npe; count", name, i, b), npebins);
        HB1(Form("%s_Hit_Npe_seg%dS_online%s; Npe; count", name, i, b), npebins);
        for(Int_t c=0; c<4; ++c)
          HB1(Form("%s_Hit_Npe_seg%d%s%s; Npe; count", name, i, abcd[c], b), npebins);
        HB1(Form("%s_Hit_Time_seg%dS%s; ns; count", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CTime_seg%dS%s; ns; count", name, i, b), hrtimebins);
      }
      HB1(Form("%sSum_Hit_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%sSum_Hit_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }

    // TOF
    {
      const Double_t phcbins2d[6] = { 100, -0.5, 4.5, 100, -10., 10. };
      for(Int_t i=0; i<NumOfSegHTOF; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("HTOF_seg%d%s_TOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
          HB2(Form("HTOF_seg%d%s_CTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("HTOF_TOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
      HB2(Form("HTOF_CTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
    }

    // BTOF
    {
      for(Int_t i=0; i<NumOfSegHodo[kBH2]; ++i){
        HB1(Form("T0_seg%d_TimeOffset%s; ns; count", i, b), 2000, -10, 10);
      }
      const Double_t phcbins2d[6] = { 100, -0.5, 4.5, 100, -10., 10. };
      for(Int_t i=0; i<NumOfSegBHT; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("BHT_seg%d%s_BTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
          HB2(Form("BHT_seg%d%s_CBTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("BHT_BTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
      HB2(Form("BHT_CBTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
      for(Int_t i=0; i<NumOfSegHodo[kBH2]; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("T0_seg%d%s_BTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
          HB2(Form("T0_seg%d%s_CBTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("T0_BTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
      HB2(Form("T0_CBTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
    }

    // FTOF
    {
      const Double_t phcbins2d[6] = { 100, -0.5, 4.5, 100, -10., 10. };
      for(const auto& id: std::vector<Int_t>{kCVC}){
        const Char_t* n = NameHodo[id];
        for(Int_t i=0; i<NumOfSegHodo[id]; ++i){
          for(const auto& uord : std::vector<TString>{"U", "D"}){
            const Char_t* ud = uord.Data();
            HB2(Form("%s_seg%d%s_FTOF_vs_DeltaE%s; mip; ns", n, i, ud, b), phcbins2d);
            HB2(Form("%s_seg%d%s_CFTOF_vs_DeltaE%s; mip; ns", n, i, ud, b), phcbins2d);
          }
        }
        HB2(Form("%s_FTOF_vs_DeltaE%s; mip; ns", n, b), phcbins2d);
        HB2(Form("%s_CFTOF_vs_DeltaE%s; mip; ns", n, b), phcbins2d);
      }
      for(Int_t i=0; i<NumOfSegHodo[kBH2]; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("T0_seg%d%s_FTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
          HB2(Form("T0_seg%d%s_CFTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("T0_FTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
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
    { // BHT
      const Char_t* name = "BHT";
      Double_t nseg = NumOfSegBHT;
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      HB2(Form("%s_Cl_MeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_CMeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_TimeDiff_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_DeltaE_vs_HitPat%s; segment; mip", name, b), debins2d);
      HB1(Form("%s_Cl_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Cl_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
      HB1(Form("%s_Cl_Size%s; size; count", name, b), 10 + 1, -0.5, 10 + 0.5);
    }
    // Hodoscope (exclude NoCluster: BAC, T1, SAC3, SFV, COBO). KVC: Cl_Npe_vs_HitPat.
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      if(HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::NoCluster)) continue;
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = NumOfSegHodo[ihodo];
      Bool_t is_cherenkov = HasHodoGroup(HodoGroupMask[ihodo], HodoGroup::Cherenkov);
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      const Double_t npebins2d[6] = { nseg, -0.5, nseg - 0.5,
        npebins[0]/10, npebins[1], npebins[2] };
      HB2(Form("%s_Cl_MeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_CMeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_TimeDiff_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_%s_vs_HitPat%s; segment; %s", name, is_cherenkov?"Npe":"DeltaE", b, is_cherenkov?"Npe":"mip"), is_cherenkov ? npebins2d : debins2d);
      HB1(Form("%s_Cl_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Cl_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
      HB1(Form("%s_Cl_Size%s; size; count", name, b), 10 + 1, -0.5, 10 + 0.5);
    }
    // BTOF
    HB1(Form("CTime0%s; ns; count", b), 400, -4, 4);
    HB1(Form("CBtof0%s; ns; count", b), 600, -20, 10);
    HB2(Form("CBtof0_vs_deT0Seg%s; mip; ns", b), 200, 0, 4, 200, -4, 4);
    HB2(Form("CBtof0_vs_deBtof0Seg%s; mip; ns", b), 200, 0, 4, 200, -4, 4);
    // FTOF
    HB1(Form("CFtof0%s; ns; count", b), 600, -10, 30);
    HB2(Form("CFtof0_vs_deT0Seg%s; mip; ns", b), 200, 0, 4, 200, -4, 4);
    HB2(Form("CFtof0_vs_deFtof0Seg%s; mip; ns", b), 200, 0, 4, 200, -4, 4);
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
      auto detector_id = digit_info.get_device_id(name);
      Int_t nplane = digit_info.get_n_plane(detector_id);
      Double_t nwire = digit_info.get_n_ch(detector_id);
      const Double_t patbins[3] = {nwire, -0.5, nwire - 0.5};
      const Double_t mulbins[3] = {nwire + 1, -0.5, nwire + 0.5};
      const Double_t tdcbins2d[6] = {nwire, -0.5, nwire - 0.5,
        mhtdcbins[0], mhtdcbins[1], mhtdcbins[2] };
      const Double_t totbins2d[6] = {nwire, -0.5, nwire - 0.5,
        mhtotbins[0], mhtotbins[1], mhtotbins[2] };
      const Double_t tdctotbins2d[6] = {
        mhtdcbins[0], mhtdcbins[1], mhtdcbins[2],
        mhtotbins[0], mhtotbins[1], mhtotbins[2] };
      for(Int_t plane=0; plane<nplane; ++plane){
        for(const auto& totcut: std::vector<TString>{"", "C"}){
          auto c = totcut.Data();
          HB1(Form("%s_%sTDC_plane%d%s; channel; count", name, c, plane, b), mhtdcbins);
          HB1(Form("%s_%sTDC1st_plane%d%s; channel; count", name, c, plane, b), mhtdcbins);
          HB1(Form("%s_%sTrailing_plane%d%s; channel; count", name, c, plane, b), mhtdcbins);
          HB1(Form("%s_%sTrailing1st_plane%d%s; channel; count", name, c, plane, b), mhtdcbins);
          HB1(Form("%s_%sTOT_plane%d%s; channel; count", name, c, plane, b), mhtotbins);
          HB1(Form("%s_%sTOT1st_plane%d%s; channel; count", name, c, plane, b), mhtotbins);
          HB1(Form("%s_%sHitPat_plane%d%s; wire; count", name, c, plane, b), patbins);
          HB1(Form("%s_%sMulti_plane%d%s; multiplicity; count", name, c, plane, b), mulbins);
          HB2(Form("%s_%sTOT_vs_TDC_plane%d%s; segment; channel", name, c, plane, b), tdctotbins2d);
          HB2(Form("%s_%sTDC_vs_HitPat_plane%d%s; segment; channel", name, c, plane, b), tdcbins2d);
          HB2(Form("%s_%sTDC1st_vs_HitPat_plane%d%s; segment; channel", name, c, plane, b), tdcbins2d);
          HB2(Form("%s_%sTrailing_vs_HitPat_plane%d%s; segment; channel", name, c, plane, b), tdcbins2d);
          HB2(Form("%s_%sTrailing1st_vs_HitPat_plane%d%s; segment; channel", name, c, plane, b), tdcbins2d);
          HB2(Form("%s_%sTOT_vs_HitPat_plane%d%s; segment; channel", name, c, plane, b), totbins2d);
          HB2(Form("%s_%sTOT1st_vs_HitPat_plane%d%s; segment; channel", name, c, plane, b), totbins2d);
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
      auto detector_id = digit_info.get_device_id(name);
      Int_t nplane = digit_info.get_n_plane(detector_id);
      Double_t nwire = digit_info.get_n_ch(detector_id);
      const Double_t patbins[3] = {nwire, -0.5, nwire - 0.5};
      const Double_t patbins2d[6] = {nwire, -0.5, nwire - 0.5, nwire, -0.5, nwire - 0.5};
      const Double_t dpatbins[3] = {nwire*2,-1*nwire-0.5,nwire-0.5};
      const Double_t mulbins[3] = {nwire + 1, -0.5, nwire + 0.5};
      const Double_t dtbins[3] = {600, -100., 400.};
      const Double_t dlbins[3] = {120/2, -1.0, 5.0};
      const Double_t dtbins2d[6] = {nwire, -0.5, nwire - 0.5,
        dtbins[0], dtbins[1], dtbins[2] };
      const Double_t dlbins2d[6] = {nwire, -0.5, nwire - 0.5,
        dlbins[0], dlbins[1], dlbins[2] };
      const Double_t dttotbins2d[6] = {
        dtbins[0], dtbins[1], dtbins[2],
        mhtotbins[0], mhtotbins[1], mhtotbins[2] };
      for(Int_t plane=0; plane<nplane; ++plane){
        HB1(Form("%s_Hit_DriftTime_plane%d%s; ns; count", name, plane, b), dtbins);
        HB1(Form("%s_Hit_DriftLength_plane%d%s; mm; count", name, plane, b), dlbins);
        HB2(Form("%s_Hit_TOT_vs_DriftTime_plane%d%s; segment; ns", name, plane, b), dttotbins2d);
        HB2(Form("%s_Hit_DriftTime_vs_HitPat_plane%d%s; segment; ns", name, plane, b), dtbins2d);
        HB2(Form("%s_Hit_DriftLength_vs_HitPat_plane%d%s; segment; mm", name, plane, b), dlbins2d);
        HB1(Form("%s_Hit_HitPat_plane%d%s; wire; count", name, plane, b), patbins);
	if(plane%2 == 0){
	  HB2(Form("%s_Hit_HitPat_Pairplane%d%d%s; wire [plane %d]; wire [plane %d]", name, plane, plane + 1, b,plane,plane+1), patbins2d);
	  HB1(Form("%s_Hit_HitPat_PP_Sub%d%d%s; wire of plane %d-wire of plane %d;count",name, plane, plane+1, b, plane, plane+1),dpatbins);
	}
        HB1(Form("%s_Hit_Multi_plane%d%s; multiplicity; count", name, plane, b), mulbins);
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
    HB1(Form("%sTrack_NHit%s; ; count", dcname.Data(), b), 20, -0.5, 19.5);
    HB1(Form("%sTrack_ChiSquare%s; ; count", dcname.Data(), b), 200, 0, 40);
    HB1(Form("%sTrack_X0%s; ; count", dcname.Data(), b), 200, -500, 500);
    HB1(Form("%sTrack_Y0%s; ; count", dcname.Data(), b), 200, -500, 500);
    HB1(Form("%sTrack_U0%s; ; count", dcname.Data(), b), 200, -.5, .5);
    HB1(Form("%sTrack_V0%s; ; count", dcname.Data(), b), 200, -.5, .5);
    for (const auto& name_str : DCNameList.at(dcname)) {
      const auto name = name_str.Data();
      auto detector_id = digit_info.get_device_id(name);
      Int_t nplane = digit_info.get_n_plane(detector_id);
      Double_t nwire = digit_info.get_n_ch(detector_id);
      const Double_t patbins[3] = {nwire, -0.5, nwire - 0.5};
      const Double_t dtbins[3] = {600, -100., 400.};
      const Double_t dlbins[3] = {120/2, -1.0, 5.0};
      const Double_t dtbins2d[6] = {nwire, -0.5, nwire - 0.5,
        dtbins[0], dtbins[1], dtbins[2] };
      const Double_t dlbins2d[6] = {nwire, -0.5, nwire - 0.5,
        dlbins[0], dlbins[1], dlbins[2] };
      const Double_t resbins[3] = {400, -2.0, 2.0};
      const Double_t resdlbins2d[6] = {200, -3., 3., 200, -2.0, 2.0};
      for(Int_t plane=0; plane<nplane; ++plane){
        HB1(Form("%s_Track_DriftTime_plane%d%s; ns; count", name, plane, b), dtbins);
        HB1(Form("%s_Track_DriftLength_plane%d%s; mm; count", name, plane, b), dlbins);
        HB2(Form("%s_Track_DriftTime_vs_HitPat_plane%d%s; segment; ns", name, plane, b), dtbins2d);
        HB2(Form("%s_Track_DriftLength_vs_HitPat_plane%d%s; segment; mm", name, plane, b), dlbins2d);
        HB1(Form("%s_Track_HitPat_plane%d%s; wire; count", name, plane, b), patbins);
        HB1(Form("%s_Track_Residual_plane%d%s; mm; count", name, plane, b), resbins);
        HB2(Form("%s_Track_Residual_vs_DriftLength_plane%d%s; mm; count", name, plane, b), resdlbins2d);
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
  const Int_t    NbinAdc     = 4096;
  const Double_t MinAdc      =    0.;
  const Double_t MaxAdc      = 4096.;
  const Int_t    NbinRms     = 1000;
  const Double_t MinRms      =    0.;
  const Double_t MaxRms      = 1000.;
  const Int_t    NbinDe      = 1000;
  const Double_t MinDe       =    0.;
  const Double_t MaxDe       = 1000.;
  const Int_t    NbinChisqr  = 1000;
  const Double_t MinChisqr   =    0.;
  const Double_t MaxChisqr   = 1000.;
  const Int_t    NbinTime    = 1000;
  const Double_t MinTime     = -8000.;
  const Double_t MaxTime     =  8000.;
  const Int_t    NbinDL      = 800;
  const Double_t MinDL       = -400.;
  const Double_t MaxDL       =  400.;
  const Int_t    NbinSigma   = 500;
  const Double_t MinSigma    =    0.;
  const Double_t MaxSigma    =   50.;
  const Int_t    NTimeBucket = 170;

  // 1D histograms
  HB1("TPC_Multiplicity_Raw",   NumOfPadTPC+1, 0,  NumOfPadTPC+1);
  HB1("TPC_Multiplicity_Cor",   NumOfPadTPC+1, 0,  NumOfPadTPC+1);
  HB1("TPC_FADC_Mean",          NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Max",           NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_RMS",           NbinRms,       MinRms,          MaxRms);
  HB1("TPC_FADC_LocMax",        NTimeBucket+1, 0,   NTimeBucket+1);
  HB1("TPC_FADC_Min",           NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Cor_Mean",      NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Cor_Max",       NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Cor_RMS",       NbinRms,       MinRms,          MaxRms);
  HB1("TPC_FADC_Cor_LocMax",    NTimeBucket+1, 0,   NTimeBucket+1);
  HB1("TPC_FADC_Cor_Min",       NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Baseline_p0",   NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Baseline_p1",   120,           -6,              6);
  HB1("TPC_FADC_Baseline_p2",   120,           -12,             12);
  HB1("TPC_FADC_Baseline_Mean", NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Baseline_Max",  NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Baseline_RMS",  NbinRms,       MinRms,          MaxRms);
  HB1("TPC_FADC_Baseline_LocMax",NTimeBucket+1, 0,   NTimeBucket+1);
  HB1("TPC_FADC_Baseline_Min",  NbinAdc,       MinAdc,          MaxAdc);

  // 2D
  HB2("TPC_FADC_Baseline",
      NTimeBucket+1, 0, NTimeBucket+1,
      NbinAdc, MinAdc, 1000);

  // TPCHit
  HB1("TPC_Multiplicity_TPCHit",   NumOfPadTPC+1, 0,  NumOfPadTPC+1);
  HB1("TPC_Pedestal",              NbinAdc,      MinAdc,  MaxAdc);
  HB1("TPC_DeltaE",                NbinDe,       MinDe,   MaxDe);
  HB1("TPC_RMS",                   NbinRms,      MinRms,  MaxRms);
  HB1("TPC_Time",                 (NTimeBucket+1)*30, 0, NTimeBucket+1);
  HB1("TPC_Chisqr",                NbinChisqr,   MinChisqr, MaxChisqr);
  HB1("TPC_CDeltaE",               NbinDe,       MinDe,   MaxDe);
  HB1("TPC_CTime",                 NbinTime,     MinTime, MaxTime);
  HB1("TPC_DriftLength",           NbinDL,       MinDL,   MaxDL);
  HB1("TPC_sigma",                 NbinSigma,    MinSigma, MaxSigma);

  HB2("TPC_sigma%%de",
      NbinDe, MinDe, MaxDe,
      NbinSigma, MinSigma, MaxSigma);

  HB2("TPC_time%%de",
      NbinDe, MinDe, MaxDe,
      NbinTime, MinTime, MaxTime);

  // FADC waveforms
  HB2("TPC_FADC_Before",
      NTimeBucket+1, 0, NTimeBucket+1,
      NbinAdc, MinAdc, MaxAdc);

  HB2("TPC_FADC_After",
      NTimeBucket+1, 0, NTimeBucket+1,
      NbinAdc, MinAdc-500., MaxAdc-500.);

  HB2("TPC_FADC_Good",
      NTimeBucket+1, 0, NTimeBucket+1,
      NbinAdc, MinAdc, MaxAdc);

  HB2("TPC_FADC_Noise",
      NTimeBucket+1, 0, NTimeBucket+1,
      NbinAdc, MinAdc, 1000);

  HB1("TPC_FADC_Noise_Max",           NbinAdc,       MinAdc,          MaxAdc);
  HB1("TPC_FADC_Noise_RMSfront",      NbinRms,       MinRms,          MaxRms);
  HB1("TPC_FADC_Noise_RMSmiddle",     NbinRms,       MinRms,          MaxRms);
  HB1("TPC_FADC_Noise_Adcdiff",       1000,-100,900);
  
  // Clock
  HB1("TPC_Clock_TDC",   100000, 0.,    1000000.);
  HB1("TPC_Clock_Time",  20000, -100.,  100.);

  HB2Poly("TPC_HitPat_Noise",-300.,300.,-300.,300.);
  HB2Poly("TPC_HitPat_Baseline",-300.,300.,-300.,300.);
  tpc::InitializeHistograms("TPC_HitPat_Noise");
  tpc::InitializeHistograms("TPC_HitPat_Baseline");
  
}

//_____________________________________________________________________________
void
BuildTPCBasic()
{
  HB1("Num_Track_TPC", 40, 0., 40.);
  HB1("Num_Track_TPC_Hits", 50, 0., 50.);
  HB1("Chisqr_TPC", 500, 0., 100.);
  HB1("X0_TPC", 400, -100., 100.);
  HB1("Y0_TPC", 400, -100., 100.);
  HB1("U0_TPC", 200, -0.20, 0.20);
  HB1("V0_TPC", 200, -0.20, 0.20);
}

//_____________________________________________________________________________
void
BuildTPCTracking()
{

  HB1("Hough_Dist", "Hough_Dist", 500, 0., 50.);
  HB1("Hough_Dist_Y", "Hough_Dist_Y", 500, 0., 50.);
  HB1("Num_Tracking_Iterations", 100, 0., 100.);
  HB1("Fitting_Flag", 10, 0., 10.);
  HB1("Track_Searching_Time", "Track_Searching_Time;Time [ms];", 100, 0., 100.);
  HB1("Track_Fitting_Time", "Track_Fitting_Time;Time [ms];", 100, 0., 100.);
  HB1("Minuit_Output_Status", 5, 0., 5.);

  HB1("Layer_Id_TPC", 35, 0., 35.);
  HB2("U0_vs_X0_TPC", "U0_vs_X0_TPC;X0;U0", 100, -100., 100., 100, -0.20, 0.20);
  HB2("V0_vs_Y0_TPC", "V0_vs_Y0_TPC;Y0;V0", 100, -100., 100., 100, -0.20, 0.20);
  HB2("X0_vs_Y0_TPC", "X0_vs_Y0_TPC;Y0;X0", 100, -100., 100., 100, -100., 100.);

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    // Tracking Histgrams
    HB1(Form("HitPat_TPC_Layer%02d", layer), Form("HitPat_TPC_Layer%02d;[Track];", layer), 400, 0., 400.);
    HB1(Form("Position_TPC_Layer%02d", layer), Form("Position_TPC_Layer%02d", layer), 200, -250., 250.);
    HB1(Form("Residual_TPC_Layer%02d", layer), Form("Residual_TPC_Layer%02d", layer), 200, 0.0, 10.0);
    HB2(Form("Residual_vs_Position_TPC_Layer%02d", layer), Form("Residual_vs_Position_TPC_Layer%02d", layer), 250, -250., 250., 100, -1.0, 1.0);
    HB2(Form("Yhit_vs_Xcal_TPC_Layer%02d", layer), Form("Yhit_vs_Xcal_TPC_Layer%02d", layer), 100, -250., 250., 100, -250., 250.);
    HB1(Form("ResidualX_TPC_Layer%02d", layer), Form("ResidualX_TPC_Layer%02d", layer), 200, -2.0, 2.0);
    HB1(Form("ResidualY_TPC_Layer%02d", layer), Form("ResidualY_TPC_Layer%02d", layer), 200, -2.0, 2.0);
    HB1(Form("ResidualZ_TPC_Layer%02d", layer), Form("ResidualZ_TPC_Layer%02d", layer), 200, -2.0, 2.0);
  }

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

  HB1("Cluster_size", "Cluster_size;Cluster size;Counts", NbinClSize, MinClSize, MaxClSize);
  HB1("Cluster_dE", "Cluster_dE;Cluster dE;Counts", NbinDe, MinDe, MaxDe);
  HB2("Ratio_vs_Dist_Transverse_diffusion", "Ratio_vs_Dist_Transverse_diffusion;X_{cluster_center}-X_{pad};A/A_{sum}", 
      NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(Form("Cluster_size_layer%2d",layer), Form("Cluster_size_layer%2d;Cluster size;Counts",layer), 
        NbinClSize, MinClSize, MaxClSize);
    HB1(Form("Cluster_dE_layer%2d",layer), Form("Cluster_dE_layer%2d;Cluster dE;Counts",layer), 
        NbinDe, MinDe, MaxDe);
    HB2(Form("Ratio_vs_Dist_Transverse_diffusion_Layer%02d",layer), 
        Form("Ratio_vs_Dist_Transverse_diffusion_Layer%02d;X_{cluster_center}-X_{pad};A/A_{sum}",layer),
        NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  }
}

//_____________________________________________________________________________
void
BuildTPCBcOutTracking()
{
 
  HB2("X0_vs_U0_TPC", "X0_vs_U0_TPC;X0;U0", 100, -100., 100., 100, -0.20, 0.20);
  HB2("Y0_vs_V0_TPC", "Y0_vs_V0_TPC;Y0;V0", 100, -100., 100., 100, -0.20, 0.20);
  HB2("X0_vs_Y0_TPC", "X0_vs_Y0_TPC;X0;Y0", 100, -100., 100., 100, -100., 100.);

  HB2("X0_vs_atanU0_TPC", "X0_vs_atan(U0)_TPC;X0;atan(U0)", 300, -300., 300., 100, -20, 20);
  HB2("Y0_vs_atanV0_TPC", "Y0_vs_atan(V0)_TPC;Y0;atan(V0)", 300, -300., 300., 100, -20, 20);

  // BcOut Basic
  HB1("Num_Track_BcOut", 40, 0., 40.);
  HB1("Num_Track_BcOut_Hits", 50, 0., 50.);
  HB1("Chisqr_BcOut", 500, 0., 500.);
  HB1("X0_BcOut", 400, -100., 100.);
  HB1("Y0_BcOut", 400, -100., 100.);
  HB1("U0_BcOut", 200, -0.20, 0.20);
  HB1("V0_BcOut", 200, -0.20, 0.20);
  HB1("Xtgt_BcOut", 400, -100., 100.);
  HB1("Ytgt_BcOut", 400, -100., 100.);
  HB1("Utgt_BcOut", 200, -0.20, 0.20);
  HB1("Vtgt_BcOut", 200, -0.20, 0.20);
  HB2("Xtgt_vs_Utgt_BcOut", "Xtgt_vs_Utgt_BcOut;Xtgt;Utgt", 100, -100., 100., 100, -0.20, 0.20);
  HB2("Ytgt_vs_Vtgt_BcOut", "Ytgt_vs_Vtgt_BcOut;Ytgt;Vtgt", 100, -100., 100., 100, -0.20, 0.20);
  HB2("Xtgt_vs_Ytgt_BcOut", "Xtgt_vs_Ytgt_BcOut;Xtgt;Ytgt", 100, -100., 100., 100, -100, 100);

  // Correlations
  HB2("Xtgt_BcOut_vs_Tpc", "Xtgt_BcOut_vs_Tpc;Tpc Xtgt;BcOut Xtgt", 400, -200., 200., 400, -200., 200.);
  HB2("Ytgt_BcOut_vs_Tpc", "Ytgt_BcOut_vs_Tpc;Tpc Ytgt;BcOut Ytgt", 400, -100., 100., 400, -100., 100.);
  HB2("Utgt_BcOut_vs_Tpc", "Utgt_BcOut_vs_Tpc;Tpc Utgt;BcOut Utgt", 400, -0.15, 0.15, 400, -0.15, 0.15);
  HB2("Vtgt_BcOut_vs_Tpc", "Vtgt_BcOut_vs_Tpc;Tpc Vtgt;BcOut Vtgt", 400, -0.05, 0.05, 400, -0.05, 0.05);
  HB1("Xtgt_Diff", "Xtgt_Diff;BcOut-Tpc", 640, -16., 16.);
  HB1("Ytgt_Diff", "Ytgt_Diff;BcOut-Tpc", 640, -16., 16.);
  HB1("Utgt_Diff", "Utgt_Diff;BcOut-Tpc", 200, -0.05, 0.05);
  HB1("Vtgt_Diff", "Vtgt_Diff;BcOut-Tpc", 200, -0.05, 0.05);

  HB2("Xtgt_Diff_vs_Xtgt_BcOut", "Xtgt_Diff_vs_Xtgt_BcOut;BcOut Xtgt;BcOut-Tpc",
       400, -150., 150., 400, -10., 10.);
  HB2("Ytgt_Diff_vs_Ytgt_BcOut", "Ytgt_Diff_vs_Ytgt_BcOut;BcOut Ytgt;BcOut-Tpc",
       400, -150., 150., 400, -15., 15.);
  HB2("Utgt_Diff_vs_Utgt_BcOut", "Utgt_Diff_vs_Utgt_BcOut;BcOut Utgt;BcOut-Tpc",
       400, -150., 150., 400, -0.03, 0.03);
  HB2("Vtgt_Diff_vs_Vtgt_BcOut", "Vtgt_Diff_vs_Vtgt_BcOut;BcOut Vtgt;BcOut-Tpc",
       400, -100., 100., 400, -0.03, 0.03);
  
  // Residuals
  HB2("Layer_vs_ResY", "Layer_vs_ResY;Layer;Y Residual", 32, 0, 32, 400, -5, 5);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(Form("TPC_Layer%02d_X_Pull", layer), 200, -5., 5.);
    HB1(Form("TPC_Layer%02d_Y_Pull", layer), 200, -5., 5.);
    HB1(Form("TPC_Layer%02d_Z_Pull", layer), 200, -5., 5.);
    HB1(Form("TPC_Layer%02d_Local_X_Pull", layer), 200, -5., 5.);
    HB1(Form("TPC_Layer%02d_Local_Y_Pull", layer), 200, -5., 5.);

    HB1(Form("TPC_Layer%02d_X_Residual", layer), 200, -8., 8.);
    HB1(Form("TPC_Layer%02d_Y_Residual", layer), 200, -8., 8.);
    HB1(Form("TPC_Layer%02d_Z_Residual", layer), 200, -8., 8.);
    HB1(Form("TPC_Layer%02d_Local_X_Residual", layer), 400, -16., 16.);
    HB1(Form("TPC_Layer%02d_Local_Y_Residual", layer), 200, -8., 8.);
    HB1(Form("TPC_Layer%02d_XZ_Residual", layer), 100, 0, 8.);

    HB1(Form("TPC_Layer%02d_BcOut_X_Residual", layer), 200, -8., 8.);
    HB1(Form("TPC_Layer%02d_BcOut_Y_Residual", layer), 200, -8., 8.);
  }

}

}