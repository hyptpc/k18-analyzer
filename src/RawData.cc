// -*- C++ -*-

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include <std_ostream.hh>
#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "HodoRawHit.hh"
#include "TPCPadHelper.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

#define OscillationCut 0

namespace
{
using namespace hddaq::unpacker;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser     = UserParamMan::GetInstance();
enum EUorD { kOneSide=1, kBothSide=2 };
enum EHodoDataType { kHodoAdc, kHodoLeading, kHodoTrailing,
  kHodoOverflow, kHodoNDataType };
#if OscillationCut
const Int_t  MaxMultiHitDC  = 16;
#endif
}

//_____________________________________________________________________________
RawData::RawData()
  : m_is_decoded(false),
    m_BH1RawHC(),
    m_BH2RawHC(),
    m_BACRawHC(),
    m_HTOFRawHC(),
    m_SCHRawHC(),
    m_BVHRawHC(),
    m_TOFRawHC(),
    m_LACRawHC(),
    m_WCRawHC(),
    m_WCSUMRawHC(),
    m_BFTRawHC(NumOfPlaneBFT),
    m_BcInRawHC(NumOfLayersBcIn+1),
    m_BcOutRawHC(NumOfLayersBcOut+1),
    m_TPCRawHC(NumOfLayersTPC+1),
    m_TPCCorHC(NumOfLayersTPC+1),
    m_SdcInRawHC(NumOfLayersSdcIn+1),
    m_SdcOutRawHC(NumOfLayersSdcOut+1),
    m_ScalerRawHC(),
    m_TrigRawHC(),
    m_VmeCalibRawHC()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
RawData::~RawData()
{
  ClearAll();
  ClearTPC();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
RawData::ClearAll()
{
  del::ClearContainer(m_BH1RawHC);
  del::ClearContainer(m_BACRawHC);
  del::ClearContainer(m_BH2RawHC);
  del::ClearContainer(m_HTOFRawHC);
  del::ClearContainer(m_SCHRawHC);
  del::ClearContainer(m_BVHRawHC);
  del::ClearContainer(m_TOFRawHC);
  del::ClearContainer(m_LACRawHC);
  del::ClearContainer(m_WCRawHC);
  del::ClearContainer(m_WCSUMRawHC);
  del::ClearContainerAll(m_BFTRawHC);
  del::ClearContainerAll(m_BcInRawHC);
  del::ClearContainerAll(m_BcOutRawHC);
  del::ClearContainerAll(m_SdcInRawHC);
  del::ClearContainerAll(m_SdcOutRawHC);
  del::ClearContainer(m_ScalerRawHC);
  del::ClearContainer(m_TrigRawHC);
  del::ClearContainer(m_VmeCalibRawHC);
}
//_____________________________________________________________________________
void
RawData::ClearTPC()
{
  del::ClearContainerAll(m_TPCRawHC);
  del::ClearContainerAll(m_TPCCorHC);
}

//_____________________________________________________________________________
bool
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
  static const Double_t MinTrailingSDC1 = gUser.GetParameter("TrailingSDC1", 0);
  static const Double_t MaxTrailingSDC1 = gUser.GetParameter("TrailingSDC1", 1);
  static const Double_t MinTrailingSDC2 = gUser.GetParameter("TrailingSDC2", 0);
  static const Double_t MaxTrailingSDC2 = gUser.GetParameter("TrailingSDC2", 1);
  static const Double_t MinTrailingSDC3 = gUser.GetParameter("TrailingSDC3", 0);
  static const Double_t MaxTrailingSDC3 = gUser.GetParameter("TrailingSDC3", 1);
  static const Double_t MinTrailingSDC4 = gUser.GetParameter("TrailingSDC4", 0);
  static const Double_t MaxTrailingSDC4 = gUser.GetParameter("TrailingSDC4", 1);

  if(m_is_decoded){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded!" << std::endl;
    return false;
  }

  ClearAll();

  // BH1
  DecodeHodo(DetIdBH1, NumOfSegBH1, kBothSide, m_BH1RawHC);
  // BH2
  DecodeHodo(DetIdBH2, NumOfSegBH2, kBothSide, m_BH2RawHC);
  // BAC
  DecodeHodo(DetIdBAC, NumOfSegBAC, kOneSide,  m_BACRawHC);
  // HTOF
  DecodeHodo(DetIdHTOF, NumOfSegHTOF, kBothSide, m_HTOFRawHC);
  // SCH
  for(Int_t seg=0; seg<NumOfSegSCH; ++seg){
    UInt_t nhit = gUnpacker.get_entries(DetIdSCH, 0, seg, 0, 0);
    for(Int_t i=0; i<nhit; ++i){
      UInt_t leading  = gUnpacker.get(DetIdSCH, 0, seg, 0, 0, i);
      UInt_t trailing = gUnpacker.get(DetIdSCH, 0, seg, 0, 1, i);
      AddHodoRawHit(m_SCHRawHC, DetIdSCH, 0, seg, 0,
                    kHodoLeading,  leading);
      AddHodoRawHit(m_SCHRawHC, DetIdSCH, 0, seg, 0,
                    kHodoTrailing, trailing);
    }
  }
  // BVH
  DecodeHodo(DetIdBVH, NumOfSegBVH, kOneSide,  m_BVHRawHC);
  // TOF
  DecodeHodo(DetIdTOF, NumOfSegTOF, kBothSide, m_TOFRawHC);
  // LAC
  DecodeHodo(DetIdLAC, NumOfSegLAC, kOneSide, m_LACRawHC);
  // WC
  DecodeHodo(DetIdWC, NumOfSegWC, kBothSide, m_WCRawHC);
  // WC SUM
  for(Int_t seg=0; seg<NumOfSegWC; ++seg){
    for(Int_t AorT=0; AorT<2; ++AorT){
      for(Int_t m=0, nhit=gUnpacker.get_entries(DetIdWC, 0, seg, 2, AorT);
          m<nhit; ++m){
	UInt_t data = gUnpacker.get(DetIdWC, 0, seg, 2, AorT, m);
	AddHodoRawHit(m_WCSUMRawHC, DetIdWCSUM, 0, seg, 0, AorT, data);
      }
    }
  }

  //BFT
  for(Int_t plane=0; plane<NumOfPlaneBFT; ++plane){
    for(Int_t seg=0; seg<NumOfSegBFT; ++seg){
      UInt_t nhit = gUnpacker.get_entries(DetIdBFT, plane, 0, seg, 0);
      for(Int_t i=0; i<nhit; ++i){
	UInt_t leading  = gUnpacker.get(DetIdBFT, plane, 0, seg, 0, i);
	UInt_t trailing = gUnpacker.get(DetIdBFT, plane, 0, seg, 1, i);
	AddHodoRawHit(m_BFTRawHC[plane], DetIdBFT, plane, seg, 0,
                      kHodoLeading,  leading);
	AddHodoRawHit(m_BFTRawHC[plane], DetIdBFT, plane, seg, 0,
                      kHodoTrailing, trailing);
      }
    }
  }

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
	    if (lt == 0 && (data<MinTdcBC3 || MaxTdcBC3<data)) continue;
	    if (lt == 1 && data<MinTdcBC3) continue;
	    AddDCRawHit(m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1,
                        data, lt);
	  }
	}
      }
    } else {
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
	    if (lt == 0 && (data<MinTdcBC4 || MaxTdcBC4<data)) continue;
	    if (lt == 1 && data<MinTdcBC4) continue;
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

  // Scaler
  for(Int_t l = 0; l<NumOfScaler; ++l){
    for(Int_t seg=0; seg<NumOfSegScaler; ++seg){
      UInt_t nhit = gUnpacker.get_entries(DetIdScaler, l, 0, seg, 0);
      if(nhit == 0) continue;
      UInt_t data = gUnpacker.get(DetIdScaler, l, 0, seg, 0);
      AddHodoRawHit(m_ScalerRawHC, DetIdScaler, l, seg, 0, 0, data);
    }
  }

  // Trigger Flag
  DecodeHodo(DetIdTrig, NumOfSegTrig, kOneSide, m_TrigRawHC);

  m_is_decoded = true;
  return true;
}

//_____________________________________________________________________________
bool
RawData::DecodeTPCHits()
{
  static const auto k_tpc = gUnpacker.get_device_id("TPC");
  static const auto k_adc = gUnpacker.get_data_id("TPC", "adc");
  // static const Int_t MinTimeBucket = gUser.GetParameter("TimeBucketTPC", 0);
  // static const Int_t MaxTimeBucket = gUser.GetParameter("TimeBucketTPC", 1);

  if(m_is_decoded){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "already decoded!" << std::endl;
    return false;
  }

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for(Int_t r=0; r<NumOfRow; ++r){
      const auto nhit = gUnpacker.get_entries(k_tpc, layer, 0, r, k_adc);
      for(Int_t i=0; i<nhit; ++i){
        auto adc = gUnpacker.get(k_tpc, layer, 0, r, k_adc, i);
	//temporary
	//if(MinTimeBucket<i&&i<MaxTimeBucket)
	AddTPCRawHit(m_TPCRawHC[layer], layer, r, adc);
      }
    }
  }

  m_is_decoded = true;
  return true;
}

//_____________________________________________________________________________
bool
RawData::RecalcTPCHits()
{
  if(!m_is_decoded){
    hddaq::cout << "#D " << FUNC_NAME << " "
		<< "Rawdata has not been decoded!" << std::endl;
    return false;
  }
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    // Int_t Min_rms_row = -1;
    // Int_t Min_rms_hitnum = -1;
    Int_t Min_rms = 100000;
    // Int_t Min_maxadc_row = -1;
    Int_t Min_maxadc_hitnum = -1;
    Int_t Min_maxadc = 100000;
    const std::size_t nh = m_TPCRawHC[layer].size();
    if(nh==0)
      continue;
    for(std::size_t hiti =0; hiti< nh; hiti++){
      auto hit = m_TPCRawHC[layer][hiti];
      if(hit->RMS() < Min_rms){
	Min_rms = hit->RMS();
	// Min_rms_row = hit->RowId();
	// Min_rms_hitnum = hiti;
      }
      if(hit->MaxAdc() < Min_maxadc){
	Min_maxadc = hit->MaxAdc();
	// Min_maxadc_row = hit->RowId();
	Min_maxadc_hitnum = hiti;
      }
    }

    // std::cout<<"Min_maxadc_hitnum:"<<Min_maxadc_hitnum<<std::endl;
    // std::cout<<"Layer:"<<layer<<", min_row:"<<Min_maxadc_row<<std::endl;
    // auto hit_min = m_TPCRawHC[layer][Min_rms_hitnum];
    auto hit_min = m_TPCRawHC[layer][Min_maxadc_hitnum];
    Int_t n0 = hit_min->Fadc().size();
    Int_t mean0 = (Int_t)hit_min->Mean();
    for(std::size_t hiti =0; hiti< nh; hiti++){
      auto hit = m_TPCRawHC[layer][hiti];
      Int_t row = hit->RowId();
      //std::cout<<"row:"<<row<<std::endl;
      Int_t n = hit->Fadc().size();
      if(n0<n)
	n=n0;
      for(Int_t i=0; i<n; ++i){
	auto adc = hit->Fadc().at(i);
	auto adc0 = hit_min->Fadc().at(i);
	auto cor_adc = adc + (mean0 - adc0);
	if(nh==1)
	  AddTPCRawHit (m_TPCCorHC[layer], layer, row, adc);
	else
	  AddTPCRawHit (m_TPCCorHC[layer], layer, row, cor_adc);
      }
    }
  }
  return true;
}
//_____________________________________________________________________________
bool
RawData::DecodeCalibHits()
{
  del::ClearContainer(m_VmeCalibRawHC);

  for(Int_t plane=0; plane<NumOfPlaneVmeCalib; ++plane){
    DecodeHodo(DetIdVmeCalib, plane, NumOfSegVmeCalib,
               kOneSide, m_VmeCalibRawHC);
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddHodoRawHit(HodoRHitContainer& cont,
                       Int_t id, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t type, Int_t data)
{
  HodoRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    HodoRawHit* q = cont[i];
    if(q->DetectorId() == id &&
       q->PlaneId() == plane &&
       q->SegmentId() == seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit(id, plane, seg);
    cont.push_back(p);
  }

  switch(type){
  case kHodoAdc:
    if(UorD==0) p->SetAdcUp(data);
    else        p->SetAdcDown(data);
    break;
  case kHodoLeading:
    if(UorD==0) p->SetTdcUp(data);
    else        p->SetTdcDown(data);
    break;
  case kHodoTrailing:
    if(UorD==0) p->SetTdcTUp(data);
    else        p->SetTdcTDown(data);
    break;
  case kHodoOverflow:
    p->SetTdcOverflow(data);
    break;
  default:
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
		<< "DetectorId = " << id    << std::endl
		<< "PlaneId    = " << plane << std::endl
		<< "SegId      = " << seg   << std::endl
		<< "AorT       = " << type  << std::endl
		<< "UorD       = " << UorD  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddDCRawHit(DCRHitContainer& cont,
                     Int_t plane, Int_t wire, Int_t data, Int_t type)
{
  DCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    DCRawHit* q = cont[i];
    if(q->PlaneId()==plane &&
       q->WireId()==wire){
      p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit(plane, wire);
    cont.push_back(p);
  }

  switch(type){
  case kDcLeading:
    p->SetTdc(data);
    break;
  case kDcTrailing:
    p->SetTrailing(data);
    break;
  case kDcOverflow:
    p->SetTdcOverflow(data);
    break;
  default:
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
		<< "PlaneId    = " << plane << std::endl
		<< "WireId     = " << wire  << std::endl
		<< "DataType   = " << type  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddTPCRawHit(TPCRHitContainer& cont,
                      Int_t layer, Int_t row, Int_t adc)
{
  TPCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    TPCRawHit* q = cont[i];
    if(q->LayerId() == layer && q->RowId() == row){
      p=q; break;
    }
  }
  if(!p){
    p = new TPCRawHit(layer, row);
    cont.push_back(p);
  }
  p->AddFadc(adc);
  return true;
}

//_____________________________________________________________________________
void
RawData::DecodeHodo(Int_t id, Int_t plane, Int_t nseg, Int_t nch,
                    HodoRHitContainer& cont)
{
  for(Int_t seg=0; seg<nseg; ++seg){
    for(Int_t UorD=0; UorD<nch; ++UorD){
      for(Int_t AorT=0; AorT<2; ++AorT){
	for(Int_t m=0, nhit=gUnpacker.get_entries(id, plane, seg, UorD, AorT);
            m<nhit; ++m){
	  UInt_t data = gUnpacker.get(id, plane, seg, UorD, AorT, m);
	  AddHodoRawHit(cont, id, plane, seg, UorD, AorT, data);
	}
      }
    }
  }
}

//_____________________________________________________________________________
void
RawData::DecodeHodo(Int_t id, Int_t nseg, Int_t nch, HodoRHitContainer& cont)
{
  DecodeHodo(id, 0, nseg, nch, cont);
}
