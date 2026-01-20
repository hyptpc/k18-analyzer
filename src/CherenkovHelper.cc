// -*- C++ -*-

#include "CherenkovHelper.hh"

#include "DetectorID.hh"
#include "HodoParamMan.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"

namespace CherenkovHelper
{
static const auto& gHodo = HodoParamMan::GetInstance();

//_____________________________________________________________________________
std::vector<Double_t>
Compute(const RawData& rawData, const TString& name)
{
  if(!gHodo.IsReady()) return {};

  if(name == "BAC"){
    Double_t s = 0.;
    const auto& cont = rawData.GetHodoRawHitContainer(name);
    for(const auto* rhit : cont){
      if(!rhit) continue;
      Int_t seg = rhit->SegmentId();
      if(seg < 0 || seg > 3) continue;  // seg4=SUM excluded
      Int_t id = rhit->DetectorId(), pl = rhit->PlaneId();
      for(const auto& adc : rhit->GetArrayAdcHigh(0)){
        Double_t npe;
        if(gHodo.GetNpe(id, pl, seg, 0, adc, npe)) s += npe;
      }
    }
    return {s};
  }

  if(name == "KVC"){
    std::vector<Double_t> out(NumOfSegKVC, 0.);
    const auto& cont = rawData.GetHodoRawHitContainer(name);
    for(const auto* rhit : cont){
      if(!rhit) continue;
      Int_t id = rhit->DetectorId(), pl = rhit->PlaneId(), seg = rhit->SegmentId();
      if(seg < 0 || seg >= NumOfSegKVC) continue;
      for(Int_t ch = 0; ch < 4; ++ch){  // ch4=kSUM excluded
        for(const auto& adc : rhit->GetArrayAdcHigh(ch)){
          Double_t npe;
          if(gHodo.GetNpe(id, pl, seg, ch, adc, npe)) out[seg] += npe;
        }
      }
    }
    return out;
  }

  return {};
}
}
