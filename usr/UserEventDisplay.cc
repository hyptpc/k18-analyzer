// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "EventDisplay.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "FLHit.hh"
#include "FuncName.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
//#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "BH2Filter.hh"

namespace
{
const auto& gGeom   = DCGeomMan::GetInstance();
auto&       gEvDisp = EventDisplay::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();
auto&       gFilter = BH2Filter::GetInstance();
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const Double_t KaonMass   = pdg::KaonMass();
const Double_t ProtonMass = pdg::ProtonMass();
}

//_____________________________________________________________________________
class UserEventDisplay : public VEvent
{
private:
  RawData*      rawData;
  DCAnalyzer*   DCAna;
  HodoAnalyzer* hodoAna;
public:
  UserEventDisplay();
  ~UserEventDisplay();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserEventDisplay::ClassName()
{
  static TString s_name("UserEventDisplay");
  return s_name;
}

//_____________________________________________________________________________
UserEventDisplay::UserEventDisplay()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
UserEventDisplay::~UserEventDisplay()
{
  if(DCAna) delete DCAna;
  if(hodoAna) delete hodoAna;
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplay::ProcessingBegin()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplay::ProcessingNormal()
{
#if 1
  // static const auto MaxMultiHitBcOut  = gUser.GetParameter("MaxMultiHitBcOut");
  static const auto MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
  static const auto MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  static const auto MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const auto MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  static const auto MinTdcSCH  = gUser.GetParameter("TdcSCH", 0);
  static const auto MaxTdcSCH  = gUser.GetParameter("TdcSCH", 1);
  static const auto MinTdcTOF  = gUser.GetParameter("TdcTOF", 0);
  static const auto MaxTdcTOF  = gUser.GetParameter("TdcTOF", 1);

  static const auto OffsetToF  = gUser.GetParameter("OffsetToF");
  // static const auto StopTimeDiffSdcOut = gUser.GetParameter("StopTimeDiffSdcOut");
  // static const auto MinStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 0);
  // static const auto MaxStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 1);
  static const auto MinTotSDC3 = gUser.GetParameter("MinTotSDC3");
  static const auto MinTotSDC4 = gUser.GetParameter("MinTotSDC4");

  // static const Int_t IdBH2 = gGeom.GetDetectorId("BH2");
  static const Int_t IdSCH = gGeom.GetDetectorId("SCH");
  static const Int_t IdTOF = gGeom.GetDetectorId("TOF");
  static const Int_t IdSDC1 = gGeom.DetectorId("SDC1-X1");
  // static const Int_t IdSDC2 = gGeom.DetectorId("SDC2-X1");
  // static const Int_t IdSDC3 = gGeom.DetectorId("SDC3-X1");
  // static const Int_t IdSDC4 = gGeom.DetectorId("SDC4-X1");

  rawData->DecodeHits();

  gEvDisp.DrawText(0.1, 0.3, Form("Run# %5d%4sEvent# %6d",
                                  gUnpacker.get_run_number(), "",
                                  gUnpacker.get_event_number()));

  ///// Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  for(auto& hit: rawData->GetTrigRawHC()){
    Int_t seg = hit->SegmentId();
    trigger_flag.set(seg);
  }

  if(trigger_flag[trigger::kSpillEnd]) return true;

  // BH2
  {
    const HodoRHitContainer &cont = rawData->GetBH2RawHC();
    Int_t nh=cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      if(!hit) continue;
      Int_t seg=hit->SegmentId();
      Int_t mh1  = hit->GetSizeTdcUp();
      Int_t mh2  = hit->GetSizeTdcDown();
      Int_t mh = 0;
      if (mh1 <= mh2)
	mh= mh1;
      else
	mh = mh2;
      for (Int_t j=0; j<mh; j++) {
	Int_t Tu=hit->GetTdcUp(j), Td=hit->GetTdcDown(j);
	std::cout << "BH2 : seg = " << seg << ", Tu = " << Tu << ", Td = " << Td << std::endl;
	if(Tu>0 && Td>0) {
	  gEvDisp.DrawBH2(seg, Td);
	} else if (Tu >0) {
	  gEvDisp.DrawBH2(seg, Tu);
	} else if (Td > 0) {
	  gEvDisp.DrawBH2(seg, Td);
	}
      }
    }
  }

  hodoAna->DecodeBH2Hits(rawData);
  Int_t nhBh2 = hodoAna->GetNHitsBH2();

  if(nhBh2==0) {
    std::cout << "Warning : nhBh2 is 0 !" << std::endl;
    //gEvDisp.GetCommand();
    return true;
  }

  Double_t time0 = -999.;
  Double_t time0_seg = -1;

  Double_t min_time = -999.;
  for(Int_t i=0; i<nhBh2; ++i){
    BH2Hit *hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    Double_t seg = hit->SegmentId()+1;
#if HodoCut
    Double_t de  = hit->DeltaE();
    if(de<MinDeBH2 || MaxDeBH2<de) continue;
#endif

    Int_t multi = hit->GetNumOfHit();
    for (Int_t m=0; m<multi; m++) {
      Double_t mt  = hit->MeanTime(m);
      // Double_t cmt = hit->CMeanTime(m);
      Double_t ct0 = hit->CTime0(m);
      if(std::abs(mt)<std::abs(min_time)){
	min_time = mt;
	time0    = ct0;
	time0_seg = seg;
      }
    }
  }

  ///// TOF
  for(const auto& hit: rawData->GetTOFRawHC()){
    if(!hit) continue;
    Int_t seg = hit->SegmentId();
    Bool_t is_hit_u = false;
    Bool_t is_hit_d = false;
    for(Int_t m=0, n=hit->GetSizeTdcUp(); m<n; ++m){
      Int_t Tu = hit->GetTdcUp(m);
      if(MinTdcTOF < Tu && Tu < MaxTdcTOF){
        is_hit_u = true;
      }
    }
    for(Int_t m=0, n=hit->GetSizeTdcDown(); m<n; ++m){
      Int_t Td = hit->GetTdcDown(m);
      if(MinTdcTOF < Td && Td < MaxTdcTOF){
        is_hit_d = true;
      }
    }
    if(is_hit_u || is_hit_d){
      gEvDisp.DrawHitHodoscope(IdTOF, seg, is_hit_u, is_hit_d);
    }
    // gEvDisp.FillTOF(seg, Tu);
    // std::cout << "TOF : seg " << seg << ", " << Tu << std::endl;
  }

  hodoAna->DecodeTOFHits(rawData);
  const auto& TOFCont = hodoAna->GetHitsTOF();
  if(TOFCont.empty()){
    std::cout << FUNC_NAME << "TOF HitContainer is empty!" << std::endl;
    //gEvDisp.GetCommand();
    return true;
  }

  ///// SCH
  for(const auto& hit: rawData->GetSCHRawHC()){
    if(!hit) continue;
    Int_t seg = hit->SegmentId();
    Bool_t is_hit = false;
    for(Int_t m=0, n=hit->GetSizeTdcUp(); m<n; ++m){
      Int_t Tu = hit->GetTdcUp(m);
      if(Tu>0) gEvDisp.FillSCH(seg, Tu);
      if(MinTdcSCH < Tu && Tu < MaxTdcSCH){
        is_hit = true;
      }
    }
    if(is_hit){
      gEvDisp.DrawHitHodoscope(IdSCH, seg);
    }
  }

  // BC Out
  for (Int_t layer=1; layer<=NumOfLayersBcOut; ++layer) {
    const DCRHitContainer &cont = rawData->GetBcOutRawHC(layer);
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      DCRawHit *hit = cont[i];
      if(!hit) continue;
      Int_t mh  = hit->GetTdcSize();
      Int_t wire = hit->WireId();
      for (Int_t j=0; j<mh; j++) {
	Int_t tdc = hit->GetTdc(j);
	if(tdc>0) gEvDisp.DrawBcOutHit(layer, wire, tdc);
      }
    }
  }

  // SDC1
  {
    const DCRHitContainer &cont = rawData->GetSdcInRawHC(IdSDC1);
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      DCRawHit *hit = cont[i];
      if(!hit) continue;
      Int_t mh  = hit->GetTdcSize();
      Int_t wire = hit->WireId();
      for (Int_t j=0; j<mh; j++) {
	Int_t tdc = hit->GetTdc(j);
	if(tdc>0) gEvDisp.DrawSDC1(wire, tdc);
      }
    }
  }
  // SDC1Xp
  {
    const DCRHitContainer &cont = rawData->GetSdcInRawHC(IdSDC1+1);
    Int_t nh = cont.size();
    for(Int_t i=0; i<nh; ++i){
      DCRawHit *hit = cont[i];
      if(!hit) continue;
      Int_t mh  = hit->GetTdcSize();
      Int_t wire = hit->WireId();
      for (Int_t j=0; j<mh; j++) {
	Int_t tdc = hit->GetTdc(j);
	if (tdc>0) gEvDisp.DrawSDC1p(wire, tdc);
      }
    }
  }

  // SDC2
  // {
  //   const DCRHitContainer &cont = rawData->GetSdcOutRawHC(IdSDC2-30);
  //   Int_t nh = cont.size();
  //   for(Int_t i=0; i<nh; ++i){
  //     DCRawHit *hit = cont[i];
  //     if(!hit) continue;
  //     Int_t mh  = hit->GetTdcSize();
  //     Int_t wire = hit->WireId();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTdc(j);
  //       if(tdc>0) gEvDisp.DrawSDC2_Leading(wire, tdc);
  //     }
  //     mh  = hit->GetTrailingSize();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTrailing(j);
  //       if(tdc>0) gEvDisp.DrawSDC2_Trailing(wire, tdc);
  //     }
  //   }
  // }

  // SDC2Xp
  // {
  //   const DCRHitContainer &cont = rawData->GetSdcOutRawHC(IdSDC2+1-30);
  //   Int_t nh = cont.size();
  //   for(Int_t i=0; i<nh; ++i){
  //     DCRawHit *hit = cont[i];
  //     if(!hit) continue;
  //     Int_t mh  = hit->GetTdcSize();
  //     Int_t wire = hit->WireId();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTdc(j);
  //       if(tdc>0) gEvDisp.DrawSDC2p_Leading(wire, tdc);
  //     }
  //     mh  = hit->GetTrailingSize();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTrailing(j);
  //       if(tdc>0) gEvDisp.DrawSDC2p_Trailing(wire, tdc);
  //     }
  //   }
  // }

  // SDC3
  // {
  //   const DCRHitContainer &cont = rawData->GetSdcOutRawHC(IdSDC3-30);
  //   Int_t nh = cont.size();
  //   for(Int_t i=0; i<nh; ++i){
  //     DCRawHit *hit = cont[i];
  //     if(!hit) continue;
  //     Int_t mh  = hit->GetTdcSize();
  //     Int_t wire = hit->WireId();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTdc(j);
  //       if(tdc>0) gEvDisp.DrawSDC3_Leading(wire, tdc);
  //     }
  //     mh  = hit->GetTrailingSize();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTrailing(j);
  //       if(tdc>0) gEvDisp.DrawSDC3_Trailing(wire, tdc);
  //     }
  //   }
  // }

  // SDC3Xp
  // {
  //   const DCRHitContainer &cont = rawData->GetSdcOutRawHC(IdSDC3+1-30);
  //   Int_t nh = cont.size();
  //   for(Int_t i=0; i<nh; ++i){
  //     DCRawHit *hit = cont[i];
  //     if(!hit) continue;
  //     Int_t mh  = hit->GetTdcSize();
  //     Int_t wire = hit->WireId();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTdc(j);
  //       if(tdc>0) gEvDisp.DrawSDC3p_Leading(wire, tdc);
  //     }
  //     mh  = hit->GetTrailingSize();
  //     for (Int_t j=0; j<mh; j++) {
  //       Int_t tdc = hit->GetTrailing(j);
  //       if(tdc>0) gEvDisp.DrawSDC3p_Trailing(wire, tdc);
  //     }
  //   }
  // }

  /*
  if (1)
    gEvDisp.GetCommand();
  return true;
  */

  // BH1
  {
    const HodoRHitContainer &cont = rawData->GetBH1RawHC();
    Int_t nh=cont.size();
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit *hit = cont[i];
      if(!hit) continue;
      Int_t seg=hit->SegmentId();
      Int_t mh1  = hit->GetSizeTdcUp();
      Int_t mh2  = hit->GetSizeTdcDown();
      Int_t mh = 0;
      if (mh1 <= mh2)
	mh= mh1;
      else
	mh = mh2;
      for (Int_t j=0; j<mh; j++) {
	Int_t Tu=hit->GetTdcUp(j), Td=hit->GetTdcDown(j);
	if(Tu>0 && Td>0) {
	  gEvDisp.DrawBH1(seg, Td);
	} else if (Tu >0) {
	  gEvDisp.DrawBH1(seg, Tu);
	} else if (Td > 0) {
	  gEvDisp.DrawBH1(seg, Td);
	}
      }
    }
  }

  // BFT raw data
  {
    for (Int_t layer=0; layer<NumOfPlaneBFT; layer++) {
      const HodoRHitContainer &cont = rawData->GetBFTRawHC(layer);
      Int_t nh = cont.size();
      for(Int_t i=0; i<nh; ++i){
	HodoRawHit *hit = cont[i];
	if(!hit) continue;
	Int_t mh  = hit->GetSizeTdcUp();

	Int_t seg = hit->SegmentId();
	for (Int_t j=0; j<mh; j++) {
	  Int_t Tu = hit->GetTdcUp(j);
	  //std::cout << "BFT-X : seg " << seg << ", " << Tu << std::endl;
	  if (Tu>0)
	    gEvDisp.DrawBFT(layer, seg, Tu);
	}
      }
    }
  }

  DCAna->DecodeRawHits(rawData);
  //DCAna->DriftTimeCutBC34(-10, 50);
  // DCAna->DriftTimeCutBC34(-100, 150);
  // BcOut
  Double_t multi_BcOut = 0.;
  {
    for(Int_t layer=1; layer<=NumOfLayersBcOut; ++layer){
      const DCHitContainer &contIn =DCAna->GetBcOutHC(layer);
      Int_t nhIn=contIn.size();
      //std::cout << "layer : " << layer << std::endl;
      for(Int_t i=0; i<nhIn; ++i){
	 // DCHit  *hit  = contIn[i];
	 // Double_t  wire = hit->GetWire();
	 // Int_t     mhit = hit->GetTdcSize();
	 //std::cout << "wire " << wire << " : ";

	 //for (Int_t j=0; j<mhit; j++) {
	 //std::cout << hit->GetDriftTime(j) << ", ";
	 //}
	 //std::cout << std::endl;

	// Bool_t    goodFlag = false;
	++multi_BcOut;
	// for (Int_t j=0; j<mhit; j++) {
	//   if (hit->IsWithinRange(j)) {
	//     goodFlag = true;
	//     break;
	//   }
	// }
	// if(goodFlag)
	//   gEvDisp.DrawHitWire(layer+112, Int_t(wire));
      }
    }
  }
  multi_BcOut /= (Double_t)NumOfLayersBcOut;
  /*
  if(multi_BcOut > MaxMultiHitBcOut) {
    gEvDisp.GetCommand();
    return true;
  }
  */

  ///// SdcIn
  Double_t multi_SdcIn = 0.;
  for(Int_t layer=1; layer<=NumOfLayersSdcIn; ++layer){
    const auto& cont = DCAna->GetSdcInHC(layer);
    Int_t n = cont.size();
    multi_SdcIn += n;
    if(n > MaxMultiHitSdcIn) continue;
    for(const auto& hit: cont){
      Double_t  wire = hit->GetWire();
      Int_t     mhit = hit->GetTdcSize();
      Bool_t    goodFlag = false;
      for(Int_t j=0; j<mhit; j++){
        if(hit->IsWithinRange(j)){
          goodFlag = true;
          break;
        }
      }
      if(goodFlag)
        gEvDisp.DrawHitWire(layer, Int_t(wire));
      else
        gEvDisp.DrawHitWire(layer, Int_t(wire), false, false);
    }
  }
  multi_SdcIn /= (Double_t)NumOfLayersSdcIn;
  if(multi_SdcIn > MaxMultiHitSdcIn) {
    std::cout << "multi_SdcIn > " << MaxMultiHitSdcIn << std::endl;
    //return true;
  }

  ///// SdcOut
  // Double_t offset = flag_tof_stop ? 0 : StopTimeDiffSdcOut;
  // DCAna->DecodeSdcOutHits(rawData);
  DCAna->TotCutSDC3(MinTotSDC3);
  DCAna->TotCutSDC4(MinTotSDC4);
  Double_t multi_SdcOut = 0.;
  for(Int_t layer=1; layer<=NumOfLayersSdcOut; ++layer){
    const auto& cont = DCAna->GetSdcOutHC(layer);
    Int_t n = cont.size();
    multi_SdcOut += n;
    if(n > MaxMultiHitSdcOut) continue;
    for(const auto& hit: cont){
      Int_t  wire = (Int_t)(hit->GetWire());
      Int_t  mhit = hit->GetDriftTimeSize();
      Bool_t is_good = false;
      for(Int_t j=0; j<mhit && !is_good; ++j){
        is_good = hit->IsWithinRange(j);
      }
      if(is_good)
        gEvDisp.DrawHitWire(layer+30, wire);
      else
        gEvDisp.DrawHitWire(layer+30, wire, false, false);
    }
  }
  multi_SdcOut /= (Double_t)NumOfLayersSdcOut;
  if(multi_SdcOut > MaxMultiHitSdcOut) {
    std::cout << "multi_SdcOut > " << MaxMultiHitSdcOut << std::endl;
    //gEvDisp.GetCommand();
    //return true;
  }

  Int_t ntBcOut = 0;
  //if(multi_BcOut<MaxMultiHitBcOut){
  if(1){
    BH2Filter::FilterList cands;
    gFilter.Apply((Int_t)time0_seg-1, *DCAna, cands);
    DCAna->TrackSearchBcOut(cands, time0_seg-1);
    //DCAna->TrackSearchBcOut(-1);
    //DCAna->TrackSearchBcOut(time0_seg-1);
    ntBcOut = DCAna->GetNtracksBcOut();
    std::cout << "NtBcOut : " << ntBcOut << std::endl;
    for(Int_t it=0; it<ntBcOut; ++it){
      DCLocalTrack *tp = DCAna->GetTrackBcOut(it);
      if(tp) gEvDisp.DrawBcOutLocalTrack(tp);
    }
  } else {
    std::cout << "Multi_BcOut over limit : " << multi_BcOut << std::endl;
  }

  Int_t ntSdcIn = 0;
  if(multi_SdcIn<MaxMultiHitSdcIn){
    std::cout << "TrackSearchSdcIn()" << std::endl;
    DCAna->TrackSearchSdcIn();
    ntSdcIn = DCAna->GetNtracksSdcIn();
    for(Int_t it=0; it<ntSdcIn; ++it){
      DCLocalTrack *tp = DCAna->GetTrackSdcIn(it);

      Int_t nh=tp->GetNHit();
      Double_t chisqr=tp->GetChiSquare();
      std::cout << "SdcIn " << it << "-th track, chi2 = " << chisqr << std::endl;
      for(Int_t ih=0; ih<nh; ++ih){
	DCLTrackHit *hit=tp->GetHit(ih);
	if(!hit) continue;
	// Int_t layerId = hit->GetLayer();
	// Double_t wire=hit->GetWire();
	// Double_t res=hit->GetResidual();
	//std::cout << "layer = " << layerId << ", wire = " << wire << ", res = " << res << std::endl;
      }
      if(tp) gEvDisp.DrawSdcInLocalTrack(tp);
    }
  }
  /*
  if(ntSdcIn==0) {
    gEvDisp.GetCommand();
    return true;
  }
  */
  Int_t ntSdcOut = 0;
  if(multi_SdcOut<MaxMultiHitSdcOut){
    std::cout << "TrackSearchSdcOut()" << std::endl;
    DCAna->TrackSearchSdcOut(TOFCont);
    ntSdcOut = DCAna->GetNtracksSdcOut();
    for(Int_t it=0; it<ntSdcOut; ++it){
      DCLocalTrack *tp = DCAna->GetTrackSdcOut(it);

      Int_t nh=tp->GetNHit();
      // Double_t chisqr=tp->GetChiSquare();
      //std::cout << "SdcOut " << it << "-th track, chi2 = " << chisqr << std::endl;
      for(Int_t ih=0; ih<nh; ++ih){
	DCLTrackHit *hit=tp->GetHit(ih);
	if(!hit) continue;
	// Int_t layerId = hit->GetLayer();
	// Double_t wire=hit->GetWire();
	// Double_t res=hit->GetResidual();
	//std::cout << "layer = " << layerId << ", wire = " << wire << ", res = " << res << std::endl;
      }
    }

    for(Int_t it=0; it<ntSdcOut; ++it){
      DCLocalTrack *tp = DCAna->GetTrackSdcOut(it);
      if(tp) gEvDisp.DrawSdcOutLocalTrack(tp);
    }
  }
  /*
  if(ntSdcOut==0) {
    gEvDisp.GetCommand();
    return true;
  }
  */
  //if (flagFBT)
  /*
  if (1)
    gEvDisp.GetCommand();
  return true;
  */

  std::vector<ThreeVector> KnPCont, KnXCont;
  std::vector<ThreeVector> KpPCont, KpXCont;

  Int_t ntKurama = 0;
  static Int_t ntKurama_all = 0;
  if(ntSdcIn>0 && ntSdcOut>0){
    //if(ntSdcIn==1 && ntSdcOut==1){

    Bool_t through_target = false;
    DCAna->TrackSearchKurama();
    ntKurama = DCAna->GetNTracksKurama();
    ntKurama_all++;
    for(Int_t it=0; it<ntKurama; ++it){
      KuramaTrack *tp = DCAna->GetKuramaTrack(it);
      if(!tp) continue;
      //tp->PrInt_t("in "+func_name);
      const ThreeVector& postgt = tp->PrimaryPosition();
      const ThreeVector& momtgt = tp->PrimaryMomentum();
      Double_t path = tp->PathLengthToTOF();
      Double_t p    = momtgt.Mag();
      gEvDisp.DrawMomentum(p);
      if(std::abs(postgt.x())<50. &&
	  std::abs(postgt.y())<30.){
	through_target = true;
      }
      // MassSquare
      Double_t tofseg = tp->TofSeg();
      for(Int_t j=0, n=TOFCont.size(); j<n; ++j){
      	Hodo2Hit *hit = TOFCont[j];
      	if(!hit) continue;
      	Double_t seg = hit->SegmentId()+1;
	if(tofseg != seg) continue;
      	Double_t stof = hit->CMeanTime()-time0+OffsetToF;
      	if(stof<=0) continue;
      	Double_t m2 = Kinematics::MassSquare(p, path, stof);
      	gEvDisp.DrawMassSquare(m2);
	KpPCont.push_back(momtgt);
	KpXCont.push_back(postgt);
      }
    }
    if(through_target) gEvDisp.DrawTarget();
    static Double_t KuramaOk = 0.;
    KuramaOk += (ntKurama>0);
  }

  //if (ntBcOut == 0)
  //  if (1)
  //gEvDisp.GetCommand();
  //return true;

  std::vector<Double_t> BftXCont;
  ////////// BFT
  {
    hodoAna->DecodeBFTHits(rawData);
    // Fiber Cluster
    hodoAna->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
    Int_t ncl = hodoAna->GetNClustersBFT();
    for(Int_t i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if(!cl) continue;
      Double_t pos    = cl->MeanPosition();
      BftXCont.push_back(pos);
    }
  }

  // K18TrackingD2U
  DCAna->TrackSearchK18D2U(BftXCont);
  Int_t ntK18=DCAna->GetNTracksK18D2U();
  if(ntK18==0) return true;
  for(Int_t i=0; i<ntK18; ++i){
    K18TrackD2U *tp=DCAna->GetK18TrackD2U(i);
    if(!tp) continue;
    Double_t x = tp->Xtgt(), y = tp->Ytgt();
    Double_t u = tp->Utgt(), v = tp->Vtgt();
    Double_t p = tp->P3rd();
    Double_t pt = p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos(x, y, 0.);
    ThreeVector Mom(pt*u, pt*v, pt);
    KnPCont.push_back(Mom);
    KnXCont.push_back(Pos);
  }

# if 1
  if(KnPCont.size()==1 && KpPCont.size()==1){
    ThreeVector pkp = KpPCont[0];
    ThreeVector pkn = KnPCont[0];
    ThreeVector xkp = KpXCont[0];
    ThreeVector xkn = KnXCont[0];
    ThreeVector vertex = Kinematics::VertexPoint(xkn, xkp, pkn, pkp);
    LorentzVector LvKn(KnPCont[0], std::sqrt(KaonMass*KaonMass+pkn.Mag2()));
    LorentzVector LvKp(KpPCont[0], std::sqrt(KaonMass*KaonMass+pkp.Mag2()));
    LorentzVector LvP(0., 0., 0., ProtonMass);
    LorentzVector LvRp = LvKn+LvP-LvKp;
    ThreeVector MissMom = LvRp.Vect();
    gEvDisp.DrawVertex(vertex);
    gEvDisp.DrawMissingMomentum(MissMom, vertex);
  }
# endif

  gEvDisp.UpdateHist();

  // if (1)
  //   gEvDisp.GetCommand();

#endif
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplay::ProcessingEnd()
{
  // gEvDisp.GetCommand();
  gEvDisp.EndOfEvent();
  // if(utility::UserStop()) gEvDisp.Run();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserEventDisplay;
}

//_____________________________________________________________________________
Bool_t
ConfMan:: InitializeHistograms()
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<DCDriftParamMan>("DCDRFT") &&
      InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>("HDPHC")       &&
      InitializeParameter<FieldMan>("FLDMAP")        &&
      InitializeParameter<K18TransMatrix>("K18TM")   &&
      InitializeParameter<BH2Filter>("BH2FLT")       &&
      InitializeParameter<UserParamMan>("USER")      &&
      InitializeParameter<EventDisplay>()           );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
