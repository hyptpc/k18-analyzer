// -*- C++ -*-

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include <TF1.h>

#include <std_ostream.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "HodoRawHit.hh"
#include "MathTools.hh"
#include "UserParamMan.hh"

#define OscillationCut 0
#define DebugEvDisp    0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
#endif

namespace
{
using namespace hddaq::unpacker;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gConf     = GConfig::get_instance();
const auto& gUser     = UserParamMan::GetInstance();
enum EUorD { kOneSide=1, kBothSide=2 };
enum EHodoDataType { kHodoAdc, kHodoLeading, kHodoTrailing,
  kHodoCstop, kHodoOverflow, kHodoNDataType };
#if OscillationCut
const Int_t  MaxMultiHitDC  = 16;
#endif

}

//_____________________________________________________________________________
RawData::RawData()
  : m_is_decoded(),
    m_hodo_raw_hit_collection(),
    m_dc_raw_hit_collection(),
    m_ScalerRawHC()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
RawData::~RawData()
{
  for(auto& elem: m_hodo_raw_hit_collection)
    del::ClearContainer(elem.second);
  for(auto& elem: m_dc_raw_hit_collection)
    del::ClearContainer(elem.second);

  m_hodo_raw_hit_collection.clear();
  m_dc_raw_hit_collection.clear();
  ClearAll();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
RawData::Clear(const TString& name)
{
  del::ClearContainer(m_hodo_raw_hit_collection[name]);
  del::ClearContainer(m_dc_raw_hit_collection[name]);
}

//_____________________________________________________________________________
void
RawData::ClearAll()
{
  del::ClearContainer(m_ScalerRawHC);
}

//_____________________________________________________________________________
Bool_t
RawData::DecodeHits(const TString& name)
{
  static const auto& digit_info = gConf.get_digit_info();

  if(m_is_decoded[name]){
    hddaq::cerr << FUNC_NAME << " " << name << " is already decoded"
                << std::endl;
    return false;
  }

  if(name.IsNull()){
    Bool_t ret = true;
    for(const auto& n: digit_info.get_name_list()){
      if(!n.empty())
        ret &= DecodeHits(n);
    }
    return ret;
  }

  if(false
     || name.Contains("null", TString::kIgnoreCase)
     || name.Contains("Scaler", TString::kIgnoreCase)
     || name.Contains("RM", TString::kIgnoreCase)){
    return false;
  }

  Clear(name);

  auto id = digit_info.get_device_id(name.Data());

  for(Int_t plane=0, n_plane=gUnpacker.get_n_plane(id);
      plane<n_plane; ++plane){
    for(Int_t seg=0, n_seg=gUnpacker.get_n_segment(id, plane);
        seg<n_seg; ++seg){
      for(Int_t ch=0, n_ch=gUnpacker.get_n_ch(id, plane, seg);
          ch<n_ch; ++ch){
        for(Int_t data=0, n_data=gUnpacker.get_n_data(id, plane, seg, ch);
            data<n_data; ++data){
          for(Int_t i=0, n=gUnpacker.get_entries(id, plane, seg, ch, data);
              i<n; ++i){
            UInt_t val = gUnpacker.get(id, plane, seg, ch, data, i);
            if(name.Contains("BC") || name.Contains("SDC")){
              AddDCRawHit(m_dc_raw_hit_collection[name], name,
                          plane, seg, ch, data, val);
            }else if(digit_info.get_n_ch(id) <= HodoRawHit::kNChannel){
              AddHodoRawHit(m_hodo_raw_hit_collection[name], name,
                            plane, seg, ch, data, val);
            }
          }
        }
      }
    }
  }

  // Print();

  m_is_decoded[name] = true;
  return true;
}

#if 0

//_____________________________________________________________________________
Bool_t
RawData::DecodeHits()
{
  static const Double_t MinTdcBC3  = gUser.GetParameter("TdcBC3", 0);
  static const Double_t MaxTdcBC3  = gUser.GetParameter("TdcBC3", 1);
  static const Double_t MinTdcBC4  = gUser.GetParameter("TdcBC4", 0);
  static const Double_t MaxTdcBC4  = gUser.GetParameter("TdcBC4", 1);
  static const Double_t MinTdcSDC1 = gUser.GetParameter("TdcSDC1", 0);
  static const Double_t MaxTdcSDC1 = gUser.GetParameter("TdcSDC1", 1);
  static const Double_t MinTdcSDC2 = gUser.GetParameter("TdcSDC2", 0);
  static const Double_t MaxTdcSDC2 = gUser.GetParameter("TdcSDC2", 1);
  static const Double_t MinTdcSDC3 = gUser.GetParameter("TdcSDC3", 0);
  static const Double_t MaxTdcSDC3 = gUser.GetParameter("TdcSDC3", 1);
  static const Double_t MinTdcSDC4 = gUser.GetParameter("TdcSDC4", 0);
  static const Double_t MaxTdcSDC4 = gUser.GetParameter("TdcSDC4", 1);
  static const Double_t MinTdcSDC5 = gUser.GetParameter("TdcSDC5", 0);
  static const Double_t MaxTdcSDC5 = gUser.GetParameter("TdcSDC5", 1);
  static const Double_t MinTrailingSDC1 = gUser.GetParameter("TrailingSDC1", 0);
  static const Double_t MaxTrailingSDC1 = gUser.GetParameter("TrailingSDC1", 1);
  static const Double_t MinTrailingSDC2 = gUser.GetParameter("TrailingSDC2", 0);
  static const Double_t MaxTrailingSDC2 = gUser.GetParameter("TrailingSDC2", 1);
  static const Double_t MinTrailingSDC3 = gUser.GetParameter("TrailingSDC3", 0);
  static const Double_t MaxTrailingSDC3 = gUser.GetParameter("TrailingSDC3", 1);
  static const Double_t MinTrailingSDC4 = gUser.GetParameter("TrailingSDC4", 0);
  static const Double_t MaxTrailingSDC4 = gUser.GetParameter("TrailingSDC4", 1);
  static const Double_t MinTrailingSDC5 = gUser.GetParameter("TrailingSDC5", 0);
  static const Double_t MaxTrailingSDC5 = gUser.GetParameter("TrailingSDC5", 1);

  // if(m_is_decoded[kOthers]){
  //   hddaq::cout << "#D " << FUNC_NAME << " "
  //       	<< "already decoded!" << std::endl;
  //   return false;
  // }

  ClearAll();

  static const auto& digit_info = gConf.get_digit_info();
  static std::vector<std::string> name_list;
  if(name_list.empty()){
    for(const auto& name: digit_info.get_name_list()){
      if(name.empty())
        continue;
      TString upper(name);
      upper.ToUpper();
      if(upper.Contains("NULL") || upper.Contains("RM"))
        continue;
      name_list.push_back(name);
    }
  }

  for(auto& name: name_list){
    auto id = digit_info.get_device_id(name);
    TString n(name);
    if(n.Contains("BC") || n.Contains("SDC")){
      ;
    }else if(digit_info.get_n_ch(id) <= HodoRawHit::kNChannel){
      DecodeHodo(name, m_hodo_raw_hit_collection[name]);
    }
  }
  return true;

  // BC3&BC4 MWDC
  for(Int_t plane=0; plane<NumOfLayersBcOut; ++plane){
    if(plane<NumOfLayersBc){
      for(Int_t wire=0; wire<MaxWireBC3; ++wire){
	for(Int_t lt=0; lt<2; ++lt){
	  UInt_t nhit = gUnpacker.get_entries(DetIdBC3, plane, 0, wire, lt);
#if OscillationCut
	  if(nhit>MaxMultiHitDC) continue;
#endif
	  for(Int_t i=0; i<nhit; i++){
	    UInt_t data = gUnpacker.get(DetIdBC3, plane, 0, wire, lt, i);
	    if(lt == 0 && (data<MinTdcBC3 || MaxTdcBC3<data)) continue;
	    if(lt == 1 && data<MinTdcBC3) continue;
	    AddDCRawHit(m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1,
			data, lt);
	  }
	}
      }
    } else {
#if 0
      // BC3 U-U' is dead E42
      if(gUnpacker.get_run_number() >= 5361
	 && (plane == NumOfLayersBc || plane == NumOfLayersBc+1)){
	continue;
      }
#endif
      for(Int_t wire=0; wire<MaxWireBC4; ++wire){
	for(Int_t lt=0; lt<2; ++lt){
	  UInt_t nhit = gUnpacker.get_entries(DetIdBC4, plane-NumOfLayersBc,
					      0, wire, lt);
#if OscillationCut
	  if(nhit>MaxMultiHitDC) continue;
#endif
	  for(Int_t i=0; i<nhit; i++){
	    UInt_t data =  gUnpacker.get(DetIdBC4, plane-NumOfLayersBc, 0,
					 wire, lt, i);
	    if(lt == 0 && (data<MinTdcBC4 || MaxTdcBC4<data)) continue;
	    if(lt == 1 && data<MinTdcBC4) continue;
	    AddDCRawHit(m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1,
			data, lt);
	  }
	}
      }
    }
  }

  // SDC1
  for(Int_t plane=0; plane<NumOfLayersSDC1; ++plane){
    for(Int_t wire=0; wire<MaxWireSDC1; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
	UInt_t nhit = gUnpacker.get_entries(DetIdSDC1, plane, 0, wire, lt);
#if OscillationCut
	if(nhit > MaxMultiHitDC) continue;
#endif
	for(Int_t i=0; i<nhit; ++i){
	  UInt_t data = gUnpacker.get(DetIdSDC1, plane, 0, wire, lt, i);
	  if(lt == 0 && (data < MinTdcSDC1 || MaxTdcSDC1 < data)) continue;
	  if(lt == 1 && (data < MinTrailingSDC1 || MaxTrailingSDC1 < data)) continue;
	  AddDCRawHit(m_SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data , lt);
	}
      }
    }
  }
  // SDC2
  for(Int_t plane=0; plane<NumOfLayersSDC2; ++plane){
    const Int_t MaxWireSDC2 = (plane < 2) ? MaxWireSDC2X : MaxWireSDC2Y;
    for(Int_t wire=0; wire<MaxWireSDC2; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
	UInt_t nhit = gUnpacker.get_entries(DetIdSDC2, plane, 0, wire, lt);
#if OscillationCut
	if(nhit > MaxMultiHitDC) continue;
#endif
	for(Int_t i=0; i<nhit; ++i){
	  UInt_t data = gUnpacker.get(DetIdSDC2, plane, 0, wire, lt, i);
	  if(lt == 0 && (data < MinTdcSDC2 || MaxTdcSDC2 < data)) continue;
	  if(lt == 1 && (data < MinTrailingSDC2 || MaxTrailingSDC2 < data)) continue;
	  AddDCRawHit(m_SdcInRawHC[plane+NumOfLayersSDC1+1],
		      plane+NumOfLayersSDC1+PlMinSdcIn, wire+1, data , lt);
	}
      }
    }
  }

  // SDC3
  for(Int_t plane=0; plane<NumOfLayersSDC3; ++plane){
    for(Int_t wire=0; wire<MaxWireSDC3; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
	auto nhit = gUnpacker.get_entries(DetIdSDC3, plane, 0, wire, lt);
#if OscillationCut
	if(nhit > MaxMultiHitDC) continue;
#endif
	for(Int_t i=0; i<nhit; ++i){
	  auto data = gUnpacker.get(DetIdSDC3, plane, 0, wire, lt, i);
	  if(lt == 0 && (data < MinTdcSDC3 || MaxTdcSDC3 < data)) continue;
	  if(lt == 1 && (data < MinTrailingSDC3 || MaxTrailingSDC3 < data)) continue;
	  AddDCRawHit(m_SdcOutRawHC[plane+1], plane+PlMinSdcOut, wire+1,
		      data , lt);
	}
      }
    }
  }

  // SDC4
  for(Int_t plane=0; plane<NumOfLayersSDC4; ++plane){
    const Int_t MaxWireSDC4 = (plane < 2) ? MaxWireSDC4Y : MaxWireSDC4X;
    for(Int_t wire=0; wire<MaxWireSDC4; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
	auto nhit = gUnpacker.get_entries(DetIdSDC4, plane, 0, wire, lt);
#if OscillationCut
	if(nhit > MaxMultiHitDC) continue;
#endif
	for(Int_t i=0; i<nhit; ++i){
	  auto data = gUnpacker.get(DetIdSDC4, plane, 0, wire, lt, i);
	  if(lt == 0 && (data < MinTdcSDC4 || MaxTdcSDC4 < data)) continue;
	  if(lt == 1 && (data < MinTrailingSDC4 || MaxTrailingSDC4 < data)) continue;
	  AddDCRawHit(m_SdcOutRawHC[plane+NumOfLayersSDC3+1],
		      plane+NumOfLayersSDC3+PlMinSdcOut, wire+1, data , lt);
	}
      }
    }
  }

  // SDC5
  for(Int_t plane=0; plane<NumOfLayersSDC5; ++plane){
    //const Int_t MaxWireSDC5 = (plane < 2) ? MaxWireSDC5X : MaxWireSDC5Y;
    const Int_t MaxWireSDC5 = (plane < 2) ? MaxWireSDC5Y : MaxWireSDC5X;
    for(Int_t wire=0; wire<MaxWireSDC5; ++wire){
      for(Int_t lt=0; lt<2; ++lt){
	auto nhit = gUnpacker.get_entries(DetIdSDC5, plane, 0, wire, lt);
#if OscillationCut
	if(nhit > MaxMultiHitDC) continue;
#endif
	for(Int_t i=0; i<nhit; ++i){
	  auto data = gUnpacker.get(DetIdSDC5, plane, 0, wire, lt, i);
	  if(lt == 0 && (data < MinTdcSDC5 || MaxTdcSDC5 < data)) continue;
	  if(lt == 1 && (data < MinTrailingSDC5 || MaxTrailingSDC5 < data)) continue;
	  AddDCRawHit(m_SdcOutRawHC[plane+NumOfLayersSDC3+NumOfLayersSDC4+1],
		      plane+NumOfLayersSDC3+NumOfLayersSDC4+PlMinSdcOut, wire+1, data , lt);
	}
      }
    }
  }

  // Scaler
  // for(Int_t l = 0; l<NumOfScaler; ++l){
  //   for(Int_t seg=0; seg<NumOfSegScaler; ++seg){
  //     UInt_t nhit = gUnpacker.get_entries(DetIdScaler, l, 0, seg, 0);
  //     if(nhit == 0) continue;
  //     UInt_t data = gUnpacker.get(DetIdScaler, l, 0, seg, 0);
  //     AddHodoRawHit(m_ScalerRawHC, DetIdScaler, l, seg, 0, 0, data);
  //   }
  // }

  // Trigger Flag
  // DecodeHodo(DetIdTrig, NumOfSegTrig, kOneSide, m_TrigRawHC);

  // m_is_decoded[kOthers] = true;
  return true;
}

#endif

//_____________________________________________________________________________
Bool_t
RawData::AddHodoRawHit(HodoRawHitContainer& cont,
		       const TString& name, Int_t plane, Int_t seg,
		       Int_t ch, Int_t data, Double_t val)
{
  HodoRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    HodoRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId() == plane
       && q->SegmentId() == seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit(name, plane, seg);
    cont.push_back(p);
  }

  // Fiber
  if(name.Contains("BFT") || name.Contains("AFT")){
    if(data == gUnpacker.get_data_id(name, "leading")){
      p->SetTdcLeading(ch, val);
    }else if(data == gUnpacker.get_data_id(name, "trailing")){
      p->SetTdcTrailing(ch, val);
    }else if(data == gUnpacker.get_data_id(name, "highgain")){
      p->SetAdcHigh(ch, val);
    }else if(data == gUnpacker.get_data_id(name, "lowgain")){
      p->SetAdcLow(ch, val);
    }else{
      ;
    }
  }
  // Hodoscope
  else{
    if(data == gUnpacker.get_data_id(name, "adc")){
      p->SetAdc(ch, val);
    }else if(data == gUnpacker.get_data_id(name, "tdc")){
      p->SetTdcLeading(ch, val);
    }else if(data == gUnpacker.get_data_id(name, "trailing")){
      p->SetTdcTrailing(ch, val);
    // }else if(data == gUnpacker.get_data_id(name, "cstop")){
    //   ;
    }else if(data == gUnpacker.get_data_id(name, "overflow")){
      p->SetTdcOverflow(ch, val);
    }
    // else{
    //   hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
    //       	<< "DetectorId = " << id    << std::endl
    //       	<< "Plane      = " << plane << std::endl
    //       	<< "Segment    = " << seg   << std::endl
    //       	<< "Channel    = " << ch    << std::endl
    //       	<< "Data       = " << data  << std::endl;
    //   return false;
    // }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddDCRawHit(DCRawHitContainer& cont,
                     const TString& name, Int_t plane, Int_t seg,
                     Int_t ch, Int_t data, Double_t val)
{
  Int_t wire = ch;
  DCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    DCRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId()==plane
       && q->WireId()==wire){
      p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit(name, plane, wire);
    cont.push_back(p);
  }

  if(data == gUnpacker.get_data_id(name, "leading")){
    p->SetTdc(val);
  }else if(data == gUnpacker.get_data_id(name, "trailing")){
    p->SetTrailing(val);
  }else if(data == gUnpacker.get_data_id(name, "overflow")){
    p->SetTdcOverflow(val);
  }else{
    hddaq::cerr << FUNC_NAME << " unknown data type " << std::endl
		<< "PlaneId    = " << plane << std::endl
		<< "WireId     = " << wire  << std::endl
		<< "DataType   = " << data  << std::endl
		<< "Value      = " << val   << std::endl;
  }
  return true;
}

//_____________________________________________________________________________
const HodoRawHitContainer&
RawData::GetHodoRawHitContainer(const TString& name) const
{
  auto itr = m_hodo_raw_hit_collection.find(name);
  if(itr == m_hodo_raw_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoRawHitContainer null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
const DCRawHitContainer&
RawData::GetDCRawHitContainer(const TString& name) const
{
  auto itr = m_dc_raw_hit_collection.find(name);
  if(itr == m_dc_raw_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static DCRawHitContainer null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
void
RawData::Print(Option_t*) const
{
  for(const auto& p: m_hodo_raw_hit_collection){
    for(const auto& hit: p.second){
      hit->Print();
    }
  }
  for(const auto& p: m_dc_raw_hit_collection){
    for(const auto& hit: p.second){
      hit->Print();
    }
  }
}
