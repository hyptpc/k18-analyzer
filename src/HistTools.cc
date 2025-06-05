// -*- C++ -*-

#if ! defined E73_2024

#include "HistTools.hh"

#include <TString.h>

#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>
#include <DAQNode.hh>

#include "DetectorID.hh"
#include "RootHelper.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
using root::HB1;
using root::HB2;
}

namespace hist
{
// Raw
const Double_t hrtdcbins1[3] = {20000, 600000, 800000};
const Double_t hrtdcbins2[3] = {50000, 0, 1000000}; // for CVC, NC
const Double_t hrtotbins[3] = {5000, 0, 50000};
const Double_t adcbins[3] = {4096, -0.5, 4095.5};
const Double_t mhtdcbins[3] = {2000, 0, 2000};
const Double_t mhtotbins[3] = {1000, 0, 1000};
// HodoHit
const Double_t hrtimebins[3] = {5000, -50, 50};
const Double_t mhtimebins[3] = {500, -50, 50};
const Double_t hrtottimebins[3] = {1000, 0, 200};
const Double_t debins[3] = {1000, 0, 10};

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
        for(const auto& uord : std::vector<TString>{"U", "D"} ){
          const Char_t* ud = uord.Data();
          HB1(Form("%s_TDC_seg%d%s%s; channel; count", name, i, ud, b), hrtdcbins1);
          HB1(Form("%s_Trailing_seg%d%s%s; channel; count", name, i, ud, b), hrtdcbins1);
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
    for(Int_t ihodo=kT0; ihodo<kNumHodo;++ihodo){
      auto name = NameHodo[ihodo].Data();
      const Double_t* hrtdcbins;
      // if(NameHodo[ihodo].Contains("CVC") ||
      //    NameHodo[ihodo].Contains("NC")){
      //   hrtdcbins = hrtdcbins2;
      // }else{
        hrtdcbins = hrtdcbins1;
      // }
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

    { ///// KVC2
      auto name = "KVC2";
      const Double_t* hrtdcbins = hrtdcbins1;
      Int_t nseg = NumOfSegKVC2;
      for(const auto& uord: std::vector<TString>{"a", "b", "c", "d"}){
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
    for(Int_t ihodo=kT0; ihodo<kNumHodo;++ihodo){
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
    // TOF
    {
      const Double_t phcbins2d[6] = { 100, -0.5, 4.5, 100, -10., 10. };
      for(Int_t i=0; i<NumOfSegBH2; ++i){
        for(const auto& uord : std::vector<TString>{"U", "D"}){
          const Char_t* ud = uord.Data();
          HB2(Form("BH2_seg%d%s_TOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
          HB2(Form("BH2_seg%d%s_CTOF_vs_DeltaE%s; mip; ns", i, ud, b), phcbins2d);
        }
      }
      HB2(Form("BH2_TOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
      HB2(Form("BH2_CTOF_vs_DeltaE%s; mip; ns", b), phcbins2d);
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
      for(Int_t i=0; i<NumOfSegHodo[kT0]; ++i){
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
      for(Int_t i=0; i<NumOfSegHodo[kT0]; ++i){
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
      for(const auto& id: std::vector<Int_t>{kT0// , kCVC
        }){
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
    for(Int_t ihodo=kT0; ihodo<kNumHodo;++ihodo){
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
      const Double_t mulbins[3] = {nwire + 1, -0.5, nwire + 0.5};
      const Double_t dtbins[3] = {600, -100., 400.};
      const Double_t dlbins[3] = {120/2, -1.0, 5.0};
      const Double_t dtbins2d[6] = {nwire, -0.5, nwire - 0.5,
        dtbins[0], dtbins[1], dtbins[2] };
      const Double_t dlbins2d[6] = {nwire, -0.5, nwire - 0.5,
        dlbins[0], dlbins[1], dlbins[2] };
      for(Int_t plane=0; plane<nplane; ++plane){
        HB1(Form("%s_Hit_DriftTime_plane%d%s; ns; count", name, plane, b), dtbins);
        HB1(Form("%s_Hit_DriftLength_plane%d%s; mm; count", name, plane, b), dlbins);
        HB2(Form("%s_Hit_DriftTime_vs_HitPat_plane%d%s; segment; ns", name, plane, b), dtbins2d);
        HB2(Form("%s_Hit_DriftLength_vs_HitPat_plane%d%s; segment; mm", name, plane, b), dlbins2d);
        HB1(Form("%s_Hit_HitPat_plane%d%s; wire; count", name, plane, b), patbins);
        HB1(Form("%s_Hit_Multi_plane%d%s; multiplicity; count", name, plane, b), mulbins);
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

#endif
