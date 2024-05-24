// -*- C++ -*-

#include "EventAnalyzer.hh"

#include <Unpacker.hh>
#include <UnpackerManager.hh>
#include <DAQNode.hh>

#include "BH2Hit.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "FiberHit.hh"
#include "FiberCluster.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
using root::HF1;
using root::HF2;
}

//_____________________________________________________________________________
EventAnalyzer::EventAnalyzer()
{
}

//_____________________________________________________________________________
EventAnalyzer::~EventAnalyzer()
{
}

//_____________________________________________________________________________
beam::EBeamFlag
EventAnalyzer::BeamFlag(const RawData& rawData)
{
  Bool_t ac_hit = false;
  {
    static const Char_t* name = "AC";
    const auto& cont = rawData.GetHodoRawHC(DetIdAC);
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      auto raw = cont[i];
      if(!raw) continue;
      for(Int_t j=0, m=raw->GetSizeTdcUp(); j<m; ++j){
	auto t = raw->GetTdcUp(j);
        if(gUser.IsInRange(Form("%s_TDC", name), t)){
          ac_hit = true; break;
        }
      }
    }
  }
  Bool_t btof_pi = false;
  Bool_t btof_k = false;
  // BHT
  {
    static const Char_t* name = "BHT";
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      for(Int_t ud=0; ud<2; ++ud){
        for(Int_t i=0, n=hit->GetSizeTdcLeading(ud); i<n; ++i){
          auto t = hit->GetTdc(ud, i);
          if(gUser.IsInRange(Form("%s_TDC_Pi", name), t))
            btof_pi = true;
          if(gUser.IsInRange(Form("%s_TDC_K", name), t))
            btof_k = true;
        }
      }
    }
  }
  beam::EBeamFlag flag;
  if(ac_hit && btof_pi)
    flag = beam::kPion;
  else if(!ac_hit && btof_k)
    flag = beam::kKaon;
  else
    flag = beam::kUnknown;
  HF1("BeamFlag", beam::kAll);
  HF1("BeamFlag", flag);
  return flag;
}

//_____________________________________________________________________________
void
EventAnalyzer::HodoRawHit(const RawData& rawData, beam::EBeamFlag beam_flag)
{
  if(beam_flag == beam::kUnknown) return;
  const Char_t* b = beam::BeamFlagList.at(beam_flag).Data();
  // BHT
  {
    static const Char_t* name = "BHT";
    Int_t multi_or = 0;
    Int_t multi_and = 0;
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      auto seg = hit->SegmentId();
      Int_t ud_good = 0;
      for(Int_t ud=0; ud<2; ++ud){
        for(Int_t i=0, n=hit->GetSizeTdcLeading(ud); i<n; ++i){
          auto t = hit->GetTdc(ud, i);
          if(gUser.IsInRange(Form("%s_TDC", name), t))
            ++ud_good;
          HF1(Form("%s_TDC_seg%d%s%s", name, seg, UorD[ud], b), t);
        }
        for(Int_t i=0, n=hit->GetSizeTdcTrailing(ud); i<n; ++i){
          auto t = hit->GetTdcTrailing(ud, i);
          HF1(Form("%s_Trailing_seg%d%s%s", name, seg, UorD[ud], b), t);
          if(n == hit->GetSizeTdcLeading(ud)){
            auto l = hit->GetTdc(ud, i);
            if(gUser.IsInRange(Form("%s_TDC", name), l)){
              HF1(Form("%s_TOT_seg%d%s%s", name, seg, UorD[ud], b), l - t);
            }
          }
        }
      }
      if(ud_good >= 1){
        HF1(Form("%s_HitPat_OR%s", name, b), seg);
        ++multi_or;
      }
      if(ud_good == 2){
        HF1(Form("%s_HitPat_AND%s", name, b), seg);
        ++multi_and;
      }
    }
    HF1(Form("%s_Multi_OR%s", name, b), multi_or);
    HF1(Form("%s_Multi_AND%s", name, b), multi_and);
  }

  // AC
  {
    static const Char_t* name = "AC";
    Int_t multi = 0;
    Bool_t is_good = false;
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      auto seg = hit->SegmentId();
      for(const auto& t: hit->GetArrayTdcLeading()){
        if(gUser.IsInRange(Form("%s_TDC", name), t))
          is_good = true;
	HF1(Form("%s_TDC_seg%d%s", name, seg, b), t);
      }
      if(is_good && seg == 0){
        HF1(Form("%s_HitPat%s", name, b), seg);
        ++multi;
      }
      auto a = hit->GetAdc();
      HF1(Form("%s_ADC_seg%d%s", name, seg, b), a);
      if(is_good)
        HF1(Form("%s_AwT_seg%d%s", name, seg, b), a);
      else
        HF1(Form("%s_AwoT_seg%d%s", name, seg, b), a);
    }
    HF1(Form("%s_Multi%s", name, b), multi);
  }

  // Hodoscope
  for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
    const Char_t* name = NameHodo[ihodo];
    Int_t multi_or = 0;
    Int_t multi_and = 0;
    for(const auto& hit: rawData.GetHodoRawHC(name)){
      auto seg = hit->SegmentId();
      Int_t ud_good = 0;
      for(Int_t ud=0; ud<2; ++ud){
        Bool_t is_good = false;
        for(const auto& t: hit->GetArrayTdc(ud)){
          if(gUser.IsInRange(Form("%s_TDC", name), t))
            is_good = true;
          HF1(Form("%s_TDC_seg%d%s%s", name, seg, UorD[ud], b), t);
        }
        for(const auto& a: hit->GetArrayAdc(ud)){
          HF1(Form("%s_ADC_seg%d%s%s", name, seg, UorD[ud], b), a);
          if(is_good)
            HF1(Form("%s_AwT_seg%d%s%s", name, seg, UorD[ud], b), a);
          else
            HF1(Form("%s_AwoT_seg%d%s%s", name, seg, UorD[ud], b), a);
        }
        ud_good += is_good;
      }
      if(ud_good >= 1){
        HF1(Form("%s_HitPat_OR%s", name, b), seg);
        ++multi_or;
      }
      if(ud_good == 2){
        HF1(Form("%s_HitPat_AND%s", name, b), seg);
        ++multi_and;
      }
    }
    HF1(Form("%s_Multi_OR%s", name, b), multi_or);
    HF1(Form("%s_Multi_AND%s", name, b), multi_and);
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::HodoHit(const HodoAnalyzer& hodoAna, beam::EBeamFlag beam_flag)
{
  if(beam_flag == beam::kUnknown) return;
  const Char_t* b = beam::BeamFlagList.at(beam_flag).Data();
  // BHT
  {
    static const Char_t* name = "BHT";
    Int_t multi = 0;
    for(Int_t i=0, n=hodoAna.GetNHits(name); i<n; ++i){
      const auto& hit = hodoAna.GetHit<FiberHit>(name, i);
      auto seg = hit->SegmentId();
      Bool_t is_good = false;
      for(Int_t j=0, m=hit->GetEntries(); j<m; ++j){
        for(Int_t ud=0; ud<2; ++ud){
          auto time = hit->GetTimeLeading(ud, j);
          auto ctime = hit->GetTimeLeading(ud, j);
          // auto t = hit->GetTimeTrailing(ud, j);
          auto tot = hit->TOT(ud, j);
          auto de = hit->GetDeltaEHighGain(ud, j);
          HF1(Form("%s_Hit_Time_seg%d%s%s", name, seg, UorD[ud], b), time);
          HF1(Form("%s_Hit_CTime_seg%d%s%s", name, seg, UorD[ud], b), ctime);
          HF1(Form("%s_Hit_TOT_seg%d%s%s", name, seg, UorD[ud], b), tot);
          HF1(Form("%s_Hit_DeltaE_seg%d%s%s", name, seg, UorD[ud], b), de);
        }
        auto mt = hit->MeanTime();
        auto cmt = hit->CMeanTime();
        auto mtot = hit->MeanTOT();
        auto de = hit->DeltaE();
        HF1(Form("%s_Hit_MeanTime_seg%d%s", name, seg, b), mt);
        HF1(Form("%s_Hit_CMeanTime_seg%d%s", name, seg, b), cmt);
        HF1(Form("%s_Hit_MeanTOT_seg%d%s", name, seg, b), mtot);
        HF1(Form("%s_Hit_DeltaE_seg%d%s", name, seg, b), de);
        HF1(Form("%s_Hit_MeanTime%s", name, b), mt);
        HF1(Form("%s_Hit_CMeanTime%s", name, b), cmt);
        HF1(Form("%s_Hit_MeanTOT%s", name, b), mtot);
        HF1(Form("%s_Hit_DeltaE%s", name, b), de);
        HF2(Form("%s_Hit_MeanTime_vs_HitPat%s", name, b), seg, mt);
        HF2(Form("%s_Hit_CMeanTime_vs_HitPat%s", name, b), seg, cmt);
        HF2(Form("%s_Hit_MeanTOT_vs_HitPat%s", name, b), seg, mtot);
        HF2(Form("%s_Hit_DeltaE_vs_HitPat%s", name, b), seg, de);
        is_good = true;
      }
      if(is_good){
        HF1(Form("%s_Hit_HitPat%s", name, b), seg);
        ++multi;
      }
    }
    HF1(Form("%s_Hit_Multi%s", name, b), multi);
  }
  // Hodoscope
  for(Int_t ihodo=kT1; ihodo<kNumHodo;++ihodo){
    const Char_t* name = NameHodo[ihodo];
    Int_t multi = 0;
    for(Int_t i=0, n=hodoAna.GetNHits(name); i<n; ++i){
      const auto& hit = hodoAna.GetHit(name, i);
      auto seg = hit->SegmentId();
      auto de = hit->DeltaE();
      auto ude  = hit->UDeltaE();
      auto dde  = hit->DDeltaE();
      HF1(Form("%s_Hit_DeltaE_seg%dU%s", name, seg, b), ude);
      HF1(Form("%s_Hit_DeltaE_seg%dD%s", name, seg, b), dde);
      HF1(Form("%s_Hit_DeltaE_seg%d%s", name, seg, b), de);
      HF2(Form("%s_Hit_DeltaE_vs_HitPat%s", name, b), seg, de);
      Bool_t is_good = false;
      for(Int_t j=0, m=hit->GetEntries(); j<m; ++j){
        auto tu  = hit->GetTUp(j),   td = hit->GetTDown(j);
        auto ctu = hit->GetCTUp(j), ctd = hit->GetCTDown(j);
        auto mt  = hit->MeanTime(j),cmt = hit->CMeanTime(j);
        HF1(Form("%s_Hit_Time_seg%dU%s", name, seg, b), tu);
        HF1(Form("%s_Hit_Time_seg%dD%s", name, seg, b), td);
        HF1(Form("%s_Hit_CTime_seg%dU%s", name, seg, b), ctu);
        HF1(Form("%s_Hit_CTime_seg%dD%s", name, seg, b), ctd);
        HF1(Form("%s_Hit_MeanTime_seg%d%s", name, seg, b), mt);
        HF1(Form("%s_Hit_CMeanTime_seg%d%s", name, seg, b), cmt);
        HF1(Form("%s_Hit_MeanTime%s", name, b), mt);
        HF1(Form("%s_Hit_CMeanTime%s", name, b), cmt);
        HF2(Form("%s_Hit_MeanTime_vs_HitPat%s", name, b), seg, mt);
        HF2(Form("%s_Hit_CMeanTime_vs_HitPat%s", name, b), seg, cmt);
        is_good = true;
      }
      if(is_good){
        HF1(Form("%s_Hit_HitPat%s", name, b), seg);
        ++multi;
      }
    }
    HF1(Form("%s_Hit_Multi%s", name, b), multi);
  }

  // BTOF / FTOF
  {
    for(Int_t i2=0, n2=hodoAna.GetNHits("T0"); i2<n2; ++i2){
      const auto& hit2 = hodoAna.GetHit<BH2Hit>("T0", i2);
      auto seg2 = hit2->SegmentId();
      auto au2 = hit2->GetAUp(), ad2 = hit2->GetADown(), a2 = hit2->DeltaE();
      for(Int_t j2=0, m2=hit2->GetEntries(); j2<m2; ++j2){
        auto tu2 = hit2->GetTUp(j2), td2 = hit2->GetTDown(j2);
        // auto ctu2 = hit2->GetCTUp(j2), ctd2 = hit2->GetCTDown(j2);
        auto mt2 = hit2->MeanTime(j2); //, cmt2 = hit2->CMeanTime(j2);
        auto t0 = hit2->Time0(j2), ct0 = hit2->CTime0(j2);
        auto tofs = hit2->TimeOffset();
        for(const auto& name: std::vector<TString>{"BHT", "CVC"}){
          const Char_t* n = name.Data();
          const Char_t* key = (name == "BHT") ? "BTOF" : "FTOF";
          for(Int_t i1=0, n1=hodoAna.GetNHits(name); i1<n1; ++i1){
            const auto& hit1 = hodoAna.GetHit(name, i1);
            auto seg1 = hit1->SegmentId();
            auto au1 = hit1->GetAUp(), ad1 = hit1->GetADown(), a1 = hit1->DeltaE();
            for(Int_t j1=0, m1=hit1->GetEntries(); j1<m1; ++j1){
              auto tu1 = hit1->GetTUp(j1), td1 = hit1->GetTDown(j1);
              // auto ctu1 = hit1->GetCTUp(j1), ctd1 = hit1->GetCTDown(j1);
              auto mt1 = hit1->MeanTime(j1), cmt1 = hit1->CMeanTime(j1);
              auto cbtof = ct0 - cmt1;
              if(name == "BHT" && TMath::Abs(seg1 - NumOfSegBHT/2) < 2) // center seg
                HF1(Form("T0_seg%d_TimeOffset%s", seg2, b), mt2 - mt1);
              if(name == "CVC" && !gUser.IsInRange("CVC_TimeDiff", hit1->TimeDiff(j1)))
                 continue;
              HF2(Form("%s_seg%dU_%s_vs_DeltaE%s", n, seg1, key, b), au1, ct0-tu1);
              HF2(Form("%s_seg%dD_%s_vs_DeltaE%s", n, seg1, key, b), ad1, ct0-td1);
              HF2(Form("%s_seg%dU_C%s_vs_DeltaE%s", n, seg1, key, b), au1, cbtof);
              HF2(Form("%s_seg%dD_C%s_vs_DeltaE%s", n, seg1, key, b), ad1, cbtof);
              HF2(Form("%s_%s_vs_DeltaE%s", n, key, b), a1, ct0-mt1);
              HF2(Form("%s_C%s_vs_DeltaE%s", n, key, b), a1, cbtof);
              HF2(Form("T0_seg%dU_%s_vs_DeltaE%s", seg2, key, b), au2, cmt1-tofs-tu2);
              HF2(Form("T0_seg%dD_%s_vs_DeltaE%s", seg2, key, b), ad2, cmt1-tofs-td2);
              HF2(Form("T0_seg%dU_C%s_vs_DeltaE%s", seg2, key, b), au2, cbtof);
              HF2(Form("T0_seg%dD_C%s_vs_DeltaE%s", seg2, key, b), ad2, cbtof);
              HF2(Form("T0_%s_vs_DeltaE%s", key, b), a2, cmt1-t0);
              HF2(Form("T0_C%s_vs_DeltaE%s", key, b), a2, cbtof);
              // HF2(Form("CTime0%s", b), ct0);
              // HF2(Form("CBtof0%s", b), cbtof);
            }
          }
        }
      }
    }
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::HodoCluster(const HodoAnalyzer& hodoAna,
                           beam::EBeamFlag beam_flag)
{
  if(beam_flag == beam::kUnknown) return;
  const Char_t* b = beam::BeamFlagList.at(beam_flag).Data();
  // Hodoscope
  for(Int_t ihodo=kBHT; ihodo<kNumHodo;++ihodo){
    const Char_t* name = NameHodo[ihodo];
    Int_t multi = 0;
    for(Int_t i=0, n=hodoAna.GetNClusters(name); i<n; ++i){
      const auto& cl = hodoAna.GetCluster(name, i);
      auto seg = cl->MeanSeg();
      // auto pos = cl->MeanPosition();
      HF2(Form("%s_Cl_MeanTime_vs_HitPat%s", name, b), seg, cl->MeanTime());
      HF2(Form("%s_Cl_CMeanTime_vs_HitPat%s", name, b), seg, cl->CMeanTime());
      HF2(Form("%s_Cl_TimeDiff_vs_HitPat%s", name, b), seg, cl->TimeDiff());
      HF2(Form("%s_Cl_DeltaE_vs_HitPat%s", name, b), seg, cl->DeltaE());
      HF1(Form("%s_Cl_HitPat%s", name, b), seg);
      HF1(Form("%s_Cl_Size%s", name, b), cl->ClusterSize());
      ++multi;
    }
    HF1(Form("%s_Cl_Multi%s", name, b), multi);
  }

  // BTOF / FTOF
  {
    auto time0 = hodoAna.Time0();
    auto btof0 = hodoAna.Btof0();
    auto ftof0 = hodoAna.Ftof0();
    HF1(Form("CTime0%s", b), time0);
    HF1(Form("CBtof0%s", b), btof0);
    HF1(Form("CFtof0%s", b), ftof0);
    const auto& cl_time0 = hodoAna.GetTime0Cluster();
    const auto& cl_btof0 = hodoAna.GetBtof0Cluster();
    const auto& cl_ftof0 = hodoAna.GetFtof0Cluster();
    if(cl_time0){
      HF2(Form("CBtof0_vs_deT0Seg%s", b), cl_time0->DeltaE(), btof0);
      HF2(Form("CFtof0_vs_deT0Seg%s", b), cl_time0->DeltaE(), ftof0);
    }
    if(cl_btof0)
      HF2(Form("CBtof0_vs_deBtof0Seg%s", b), cl_btof0->DeltaE(), btof0);
    if(cl_ftof0)
      HF2(Form("CFtof0_vs_deFtof0Seg%s", b), cl_ftof0->DeltaE(), ftof0);
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::DCRawHit(const RawData& rawData, beam::EBeamFlag beam_flag)
{
  if(beam_flag == beam::kUnknown) return;
  const Char_t* b = beam::BeamFlagList.at(beam_flag).Data();
  // DC
  for(Int_t idc=0; idc<=kVFT; ++idc){
    const Char_t* name = NameDC[idc].Data();
    auto nlayer = NumOfLayerDC[idc];
    for(Int_t layer=0; layer<nlayer; ++layer){
      const auto& cont = rawData.GetDCRawHC(DetIdDC[idc], layer);
      // hist::H1(Form("%s_Mul_layer%d",name,layer),nh,mulbins);
      Int_t multi = 0;
      Int_t cmulti = 0;
      for(Int_t i=0, n=cont.size(); i<n; ++i){
        auto hit = cont[i];
        if(!hit) continue;
        auto wire = hit->WireId();
        Bool_t is_good = false;
        for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
          if(m != hit->GetTrailingSize()) break;
          auto l = hit->GetTdc(j);
          auto t = hit->GetTrailing(j);
          auto tot = l - t;
          if(gUser.IsInRange(Form("%s_TOT", name), tot) &&
             gUser.IsInRange(Form("%s_TDC", name), l))
            is_good = true;
          for(const auto& totcut: std::vector<TString>{"", "C"}){
            if(totcut == "C" && !gUser.IsInRange(Form("%s_TOT", name), tot))
              continue;
            auto c = totcut.Data();
            HF1(Form("%s_%sTDC_layer%d%s", name, c, layer, b), l);
            HF1(Form("%s_%sTrailing_layer%d%s", name, c, layer, b), t);
            HF1(Form("%s_%sTOT_layer%d%s", name, c, layer, b), tot);
            if(j == 0){
              HF1(Form("%s_%sTDC1st_layer%d%s", name, c, layer, b), l);
              HF1(Form("%s_%sTrailing1st_layer%d%s", name, c, layer, b), t);
              HF1(Form("%s_%sTOT1st_layer%d%s", name, c, layer, b), tot);
            }
          }
        }
        if(hit->GetTdcSize() > 0){
          HF1(Form("%s_HitPat_layer%d%s", name, layer, b), wire);
          ++multi;
        }
        if(is_good){
          HF1(Form("%s_CHitPat_layer%d%s", name, layer, b), wire);
          ++cmulti;
        }
      }
      HF1(Form("%s_Multi_layer%d%s", name, layer, b), multi);
      HF1(Form("%s_CMulti_layer%d%s", name, layer, b), cmulti);
    }
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::TriggerFlag(const RawData& rawData)
{
  static const Char_t* name = "TriggerFlag";
  for(const auto& hit: rawData.GetHodoRawHC(name)){
    auto seg = hit->SegmentId();
    if(seg > NumOfSegTrigFlag)
      continue;
    Bool_t is_hit = false;
    for(const auto& tdc: hit->GetArrayTdcUp()){
      HF1(Form("%s_TDC_seg%d", name, seg), tdc);
      is_hit = true;
    }
    if(is_hit){
      HF1(Form("%s_HitPat", name), seg);
    }
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::DAQ(const RawData& rawData)
{
  static const auto k_data_size = hddaq::unpacker::DAQNode::k_data_size;
  Int_t vme_index = 0;
  Int_t hul_index = 0;
  Int_t vea0c_index = 0;
  for(auto&& c : gUnpacker.get_root()->get_child_list()){
    if (!c.second) continue;
    TString name = c.second->get_name();
    auto node_id = c.second->get_id();
    auto data_size = gUnpacker.get_node_header(node_id, k_data_size);
    if(name.Contains("vme")){
      root::HF2("FE_VME_DataSize", vme_index++, data_size);
    }
    if(name.Contains("hul")){
      root::HF2("FE_HUL_DataSize", hul_index++, data_size);
    }
    if(name.Contains("veasiroc")){
      root::HF2("FE_VEASIROC_DataSize", vea0c_index++, data_size);
    }
  }
  auto node_id = gUnpacker.get_fe_id("k18breb");
  auto data_size = gUnpacker.get_node_header(node_id, k_data_size);
  HF1("EB_DataSize", data_size);
}
