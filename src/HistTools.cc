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
    // Hodoscope
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      const Double_t* hrtdcbins;
      if ( ihodo == kCOBO ) {
        hrtdcbins = hrtdcbins4;
      } else if ( ihodo == kCVC || ihodo == kSFV || ihodo == kSAC3 ) {
        hrtdcbins = hrtdcbins2;
      } else {
        hrtdcbins = hrtdcbins1;
      }
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

    { ///// HTOF
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

    { ///// KVC
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
    // Hodoscope
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = NumOfSegHodo[ihodo];
      for(Int_t i=0; i<nseg; ++i){
        for(const auto& uord: std::vector<TString>{"U", "D"} ){
          auto ud = uord.Data();
          HB1(Form("%s_Hit_DeltaE_seg%d%s%s; mip; count", name, i, ud, b), debins);
          HB1(Form("%s_Hit_Time_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
          HB1(Form("%s_Hit_CTime_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
        }
        HB1(Form("%s_Hit_DeltaE_seg%d%s; mip; count", name, i, b), debins);
        HB1(Form("%s_Hit_MeanTime_seg%d%s; ns; count", name, i, b), hrtimebins);
        HB1(Form("%s_Hit_CMeanTime_seg%d%s; ns; count", name, i, b), hrtimebins);
      }
      HB1(Form("%s_Hit_MeanTime%s; ns; count", name, b), hrtimebins);
      HB1(Form("%s_Hit_CMeanTime%s; ns; count", name, b), hrtimebins);
      HB1(Form("%s_Hit_DeltaE%s; mip; count", name, b), hrtottimebins);
      const Double_t hrtimebins2d[6] = { nseg, -0.5, nseg - 0.5,
        hrtimebins[0]/10, hrtimebins[1], hrtimebins[2] };
      const Double_t debins2d[6] = { nseg, -0.5, nseg - 0.5,
        debins[0]/10, debins[1], debins[2] };
      HB2(Form("%s_Hit_MeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_CMeanTime_vs_HitPat%s; segment; ns", name, b), hrtimebins2d);
      HB2(Form("%s_Hit_DeltaE_vs_HitPat%s; segment; mip", name, b), debins2d);
      HB1(Form("%s_Hit_HitPat%s; segment; count", name, b), nseg, -0.5, nseg - 0.5);
      HB1(Form("%s_Hit_Multi%s; multiplicity; count", name, b), nseg + 1, -0.5, nseg + 0.5);
    }

    // HTOF, KVC Sum
    for (const auto& ihodo: std::vector<Int_t>{kHTOF, kKVC}) {
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = NumOfSegHodo[ihodo];
      for(Int_t i=0; i<nseg; ++i){
        const Char_t* ud = "S";
        HB1(Form("%s_Hit_DeltaE_seg%d%s%s; mip; count", name, i, ud, b), debins);
        HB1(Form("%s_Hit_Time_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
        HB1(Form("%s_Hit_CTime_seg%d%s%s; ns; count", name, i, ud, b), hrtimebins);
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
    // Hodoscope
    for(Int_t ihodo=kBH2; ihodo<kNumHodo;++ihodo){
      if (ihodo == kBAC || ihodo == kT1 || ihodo == kSAC3 || ihodo == kSFV || ihodo == kCOBO) continue;
      auto name = NameHodo[ihodo].Data();
      Double_t nseg = NumOfSegHodo[ihodo];
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
// TODO: clean up this code.
void
BuildTPCTracking()
{

  // HB1(const TString& name, const TString& title,
  //     Int_t nbinx, Double_t xlow, Double_t xhigh)

  // HB2(const TString& name, const TString& title,
  //   Int_t nbinx, Double_t xlow, Double_t xhigh,
  //   Int_t nbiny, Double_t ylow, Double_t yhigh)

  HB1("Hough_Dist", "Hough_Dist", 500, 0., 50.);
  HB1("Hough_Dist_Y", "Hough_Dist_Y", 500, 0., 50.);
  HB1("Num_Tracking_Iterations", 100, 0., 100.);
  HB1("Fitting_Flag", 10, 0., 10.);
  HB1("Track_Searching_Time", "Track_Searching_Time;Time [ms];", 100, 0., 100.);
  HB1("Track_Fitting_Time", "Track_Fitting_Time;Time [ms];", 100, 0., 100.);
  HB1("Minuit_Output_Status", 5, 0., 5.);

  HB1("Num_Track_TPC", 40, 0., 40.);
  HB1("Num_Track_TPC_Hits", 50, 0., 50.);
  HB1("Chisqr_TPC", 500, 0., 100.);
  HB1("Layer_Id_TPC", 35, 0., 35.);
  HB1("X0_TPC", 400, -100., 100.);
  HB1("Y0_TPC", 400, -100., 100.);
  HB1("U0_TPC", 200, -0.20, 0.20);
  HB1("V0_TPC", 200, -0.20, 0.20);
  HB2("U0_X0_TPC", "U0_X0_TPC;X0;U0", 100, -100., 100., 100, -0.20, 0.20);
  HB2("V0_Y0_TPC", "V0_Y0_TPC;Y0;V0", 100, -100., 100., 100, -0.20, 0.20);
  HB2("X0_Y0_TPC", "X0_Y0_TPC;Y0;X0", 100, -100., 100., 100, -100., 100.);

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    // Tracking Histgrams
    HB1(Form("HitPat_TPC_Layer%02d", layer), Form("HitPat_TPC_Layer%02d;[Track];", layer), 400, 0., 400.);
    HB1(Form("Position_TPC_Layer%02d", layer), Form("Position_TPC_Layer%02d", layer), 200, -250., 250.);
    HB1(Form("Residual_TPC_Layer%02d", layer), Form("Residual_TPC_Layer%02d", layer), 200, 0.0, 10.0);
    HB2(Form("Resid_vs_Pos_TPC_Layer%02d", layer), Form("Resid_vs_Pos_TPC_Layer%02d", layer), 250, -250., 250., 100, -1.0, 1.0);
    HB2(Form("Y_vs_Xcal_TPC_Layer%02d", layer), Form("Y_vs_Xcal_TPC_Layer%02d", layer), 100, -250., 250., 100, -250., 250.);
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
  HB2("Transverse_diffusion", "Transverse_diffusion;X_{cluster_center}-X_{pad};A/A_{sum}", 
      NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    HB1(Form("Cluster_size_layer%2d",layer), Form("Cluster_size_layer%2d;Cluster size;Counts",layer), 
        NbinClSize, MinClSize, MaxClSize);
    HB1(Form("Cluster_dE_layer%2d",layer), Form("Cluster_dE_layer%2d;Cluster dE;Counts",layer), 
        NbinDe, MinDe, MaxDe);
    HB2(Form("Transverse_diffusion_Layer%02d",layer), 
        Form("Transverse_diffusion_Layer%02d;X_{cluster_center}-X_{pad};A/A_{sum}",layer),
        NbinDist, MinDist, MaxDist, NbinRatio, MinRatio, MaxRatio);
  }
}

}