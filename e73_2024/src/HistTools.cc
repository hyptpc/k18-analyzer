// -*- C++ -*-

#include "HistTools.hh"

#include <TString.h>

#include <Unpacker.hh>
#include <UnpackerManager.hh>
#include <DAQNode.hh>

#include "DetectorID.hh"
#include "RootHelper.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
using root::HB1;
using root::HB2;
}

namespace hist
{
// Raw
const Double_t hrtdcbins1[3] = {45000, 950000, 1400000};
const Double_t hrtdcbins2[3] = {50000, 0, 1000000}; // for CVC, NC
const Double_t hrtotbins[3] = {5000, 0, 50000};
const Double_t adcbins[3] = {4096, -0.5, 4095.5};
const Double_t mhtdcbins[3] = {2000, 0, 2000};
const Double_t mhtotbins[3] = {1000, 0, 1000};
// HodoHit
const Double_t hrtimebins[3] = {5000, -50, 50};
const Double_t mhtimebins[3] = {500, -50, 50};
const Double_t hrtottimebins[3] = {1000, 0, 200};
const Double_t debins[3] = {200, 0, 10};

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
    HB1(Form("%s_TDC_seg%d", name, i), mhtdcbins);
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
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB1(Form("%s_TDC_seg%d%s%s", name, i, ud, b), hrtdcbins1);
          HB1(Form("%s_Trailing_seg%d%s%s", name, i, ud, b), hrtdcbins1);
          HB1(Form("%s_TOT_seg%d%s%s", name, i, ud, b), hrtotbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s", name, ud, b), nseg, -0.5, nseg - 0.5);
        HB1(Form("%s_Multi_%s%s", name, ud, b), nseg + 1, -0.5, nseg + 0.5);
      }
    }
    { // AC
      const Char_t* name = "AC";
      Int_t nseg = NumOfSegAC;
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_ADC_seg%d%s", name, i, b), adcbins);
        HB1(Form("%s_AwT_seg%d%s", name, i, b), adcbins);
        HB1(Form("%s_AwoT_seg%d%s", name, i, b), adcbins);
        HB1(Form("%s_TDC_seg%d%s", name, i, b), mhtdcbins);
      }
      HB1(Form("%s_HitPat%s", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Multi%s", name, b), nseg + 1, -0.5, nseg + 0.5);
    }
    // Hodoscope
    for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      const Double_t* hrtdcbins;
      if(NameHodo[ihodo].Contains("CVC") ||
         NameHodo[ihodo].Contains("NC")){
        hrtdcbins = hrtdcbins2;
      }else{
        hrtdcbins = hrtdcbins1;
      }
      Int_t nseg = NumOfSegHodo[ihodo];
      for(const auto& uord: std::vector<TString>{"U", "D"} ){
        auto ud = uord.Data();
        for(Int_t i=0; i<nseg; ++i){
          HB1(Form("%s_ADC_seg%d%s%s", name, i, ud, b), adcbins);
          HB1(Form("%s_AwT_seg%d%s%s", name, i, ud, b), adcbins);
          HB1(Form("%s_AwoT_seg%d%s%s", name, i, ud, b), adcbins);
          HB1(Form("%s_TDC_seg%d%s%s", name, i, ud, b), hrtdcbins);
        }
      }
      for(const auto& uord: std::vector<TString>{"OR", "AND"} ){
        auto ud = uord.Data();
        HB1(Form("%s_HitPat_%s%s", name, ud, b), nseg, -0.5, nseg - 0.5);
        HB1(Form("%s_Multi_%s%s", name, ud, b), nseg + 1, -0.5, nseg + 0.5);
      }
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
          HB1(Form("%s_Hit_Time_seg%d%s%s", name, i, ud, b), hrtimebins);
          HB1(Form("%s_Hit_CTime_seg%d%s%s", name, i, ud, b), hrtimebins);
          HB1(Form("%s_Hit_TOT_seg%d%s%s", name, i, ud, b), hrtottimebins);
          HB1(Form("%s_Hit_DeltaE_seg%d%s%s", name, i, ud, b), debins);
        }
        HB1(Form("%s_Hit_MeanTime_seg%d%s", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CMeanTime_seg%d%s", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_MeanTOT_seg%d%s", name, i, b), hrtottimebins);
        HB1(Form("%s_Hit_DeltaE_seg%d%s", name, i, b), debins);
      }
      HB1(Form("%s_Hit_MeanTime%s", name, b), hrtimebins);
      HB1(Form("%s_Hit_CMeanTime%s", name, b), hrtimebins);
      HB1(Form("%s_Hit_MeanTOT%s", name, b), hrtottimebins);
      HB1(Form("%s_Hit_DeltaE%s", name, b), debins);
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t hrtottimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtottimebins[0]/5, hrtottimebins[1], hrtottimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      HB2(Form("%s_Hit_MeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_CMeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_MeanTOT_vs_HitPat%s", name, b), hrtottimebins2d);
      HB2(Form("%s_Hit_DeltaE_vs_HitPat%s", name, b), debins2d);
      HB1(Form("%s_Hit_HitPat%s", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Hit_Multi%s", name, b), nseg + 1, -0.5, nseg + 0.5);
    }
    { // AC
      const Char_t* name = "AC";
      Int_t nseg = NumOfSegAC;
      for(Int_t i=0; i<nseg; ++i){
        HB1(Form("%s_Hit_DeltaE_seg%d%s", name, i, b), debins);
        HB1(Form("%s_Hit_Time_seg%d%s", name, i, b), mhtimebins);
        HB1(Form("%s_Hit_CTime_seg%d%s", name, i, b), mhtimebins);
      }
      HB1(Form("%s_Hit_HitPat%s", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Hit_Multi%s", name, b), nseg + 1, -0.5, nseg + 0.5);
    }
    // Hodoscope
    for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = NumOfSegHodo[ihodo];
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord: std::vector<TString>{"U", "D"} ){
          auto ud = uord.Data();
          HB1(Form("%s_Hit_DeltaE_seg%d%s%s", name, i, ud, b), debins);
          HB1(Form("%s_Hit_Time_seg%d%s%s", name, i, ud, b), hrtimebins);
          HB1(Form("%s_Hit_CTime_seg%d%s%s", name, i, ud, b), hrtimebins);
        }
        HB1(Form("%s_Hit_DeltaE_seg%d%s", name, i, b), debins);
        HB1(Form("%s_Hit_MeanTime_seg%d%s", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CMeanTime_seg%d%s", name, i, b), hrtimebins);
      }
      HB1(Form("%s_Hit_MeanTime%s", name, b), hrtimebins);
      HB1(Form("%s_Hit_CMeanTime%s", name, b), hrtimebins);
      HB1(Form("%s_Hit_DeltaE%s", name, b), hrtottimebins);
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      HB2(Form("%s_Hit_MeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_CMeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_DeltaE_vs_HitPat%s", name, b), debins2d);
      HB1(Form("%s_Hit_HitPat%s", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Hit_Multi%s", name, b), nseg + 1, -0.5, nseg + 0.5);
    }
    // BTOF
    {
      for(Int_t i=0; i<NumOfSegHodo[kT0]; ++i){
        HB1(Form("T0_seg%d_TimeOffset%s", i, b), 2000, -10, 10);
      }
      const Double_t phcbins2d[6] = { 100, -0.5, 4.5, 100, -10., 10. };
      for(Int_t i=0; i<NumOfSegBHT; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("BHT_seg%d%s_BTOF_vs_DeltaE%s", i, ud, b), phcbins2d);
          HB2(Form("BHT_seg%d%s_CBTOF_vs_DeltaE%s", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("BHT_BTOF_vs_DeltaE%s", b), phcbins2d);
      HB2(Form("BHT_CBTOF_vs_DeltaE%s", b), phcbins2d);
      for(Int_t i=0; i<NumOfSegHodo[kT0]; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("T0_seg%d%s_BTOF_vs_DeltaE%s", i, ud, b), phcbins2d);
          HB2(Form("T0_seg%d%s_CBTOF_vs_DeltaE%s", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("T0_BTOF_vs_DeltaE%s", b), phcbins2d);
      HB2(Form("T0_CBTOF_vs_DeltaE%s", b), phcbins2d);
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
      HB2(Form("%s_Cl_MeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_CMeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_DeltaE_vs_HitPat%s", name, b), debins2d);
      HB1(Form("%s_Cl_HitPat%s", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Cl_Multi%s", name, b), nseg + 1, -0.5, nseg + 0.5);
      HB1(Form("%s_Cl_Size%s", name, b), 10 + 1, -0.5, 10 + 0.5);
    }
    { // AC
      const Char_t* name = "AC";
      Int_t nseg = NumOfSegAC;
      HB1(Form("%s_Cl_DeltaE%s", name, b), debins);
      HB1(Form("%s_Cl_Time%s", name, b), mhtimebins);
      HB1(Form("%s_Cl_CTime%s", name, b), mhtimebins);
      HB1(Form("%s_Cl_HitPat%s", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Cl_Multi%s", name, b), nseg + 1, -0.5, nseg + 0.5);
      HB1(Form("%s_Cl_Size%s", name, b), 10 + 1, -0.5, 10 + 0.5);
    }
    // Hodoscope
    for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = NumOfSegHodo[ihodo];
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      HB2(Form("%s_Cl_MeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_CMeanTime_vs_HitPat%s", name, b), hrtimebins2d);
      HB2(Form("%s_Cl_DeltaE_vs_HitPat%s", name, b), debins2d);
      HB1(Form("%s_Cl_HitPat%s", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Cl_Multi%s", name, b), nseg + 1, -0.5, nseg + 0.5);
      HB1(Form("%s_Cl_Size%s", name, b), 10 + 1, -0.5, 10 + 0.5);
    }
    // BTOF
    HB1(Form("CTime0%s", b), 400, -4, 4);
    HB1(Form("CBtof0%s", b), 1000, -10, 10);
    HB2(Form("CBtof0_vs_deT0Seg%s", b), 200, 0, 4, 200, -4, 4);
    HB2(Form("CBtof0_vs_deBtof0Seg%s", b), 200, 0, 4, 200, -4, 4);
    if(!flag_beam_particle) break;
  }
}

//_____________________________________________________________________________
void
BuildDCRaw(Bool_t flag_beam_particle)
{
  for(const auto& beam: beam::BeamFlagList){
    const Char_t* b = beam.Data();
    for(Int_t idc=0; idc<=kVFT; ++idc){
      const Char_t* name = NameDC[idc].Data();
      Int_t nlayer = NumOfLayerDC[idc];
      Double_t nwire = NumOfWireDC[idc];
      const Double_t patbins[3] = {nwire, -0.5, nwire - 0.5};
      const Double_t mulbins[3] = {nwire + 1, -0.5, nwire + 0.5};
      for(Int_t layer=0; layer<nlayer; ++layer){
        for(const auto& totcut: std::vector<TString>{"", "C"}){
          auto c = totcut.Data();
          HB1(Form("%s_%sTDC_layer%d%s", name, c, layer, b), mhtdcbins);
          HB1(Form("%s_%sTDC1st_layer%d%s", name, c, layer, b), mhtdcbins);
          HB1(Form("%s_%sTrailing_layer%d%s", name, c, layer, b), mhtdcbins);
          HB1(Form("%s_%sTrailing1st_layer%d%s", name, c, layer, b), mhtdcbins);
          HB1(Form("%s_%sTOT_layer%d%s", name, c, layer, b), mhtotbins);
          HB1(Form("%s_%sTOT1st_layer%d%s", name, c, layer, b), mhtotbins);
          HB1(Form("%s_%sHitPat_layer%d%s", name, c, layer, b), patbins);
          HB1(Form("%s_%sMulti_layer%d%s", name, c, layer, b), mulbins);
        }
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
}
