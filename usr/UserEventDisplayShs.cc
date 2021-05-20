// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>

#include <TMath.h>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "EventDisplayShs.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"
#include "DetectorID.hh"

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
auto& gEvDisp = EventDisplayShs::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
class UserEventDisplayShs : public VEvent
{
private:
  RawData*    rawData;
  DCAnalyzer* DCAna;

public:
  UserEventDisplayShs();
  ~UserEventDisplayShs();
  virtual const TString& ClassName();
  virtual Bool_t         ProcessingBegin();
  virtual Bool_t         ProcessingEnd();
  virtual Bool_t         ProcessingNormal();
};

//_____________________________________________________________________________
inline const TString&
UserEventDisplayShs::ClassName()
{
  static TString s_name("UserEventDisplayShs");
  return s_name;
}


//_____________________________________________________________________________
UserEventDisplayShs::UserEventDisplayShs()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer)
{
}

//_____________________________________________________________________________
UserEventDisplayShs::~UserEventDisplayShs()
{
  if(rawData) delete rawData;
  if(DCAna) delete DCAna;
}

//_____________________________________________________________________________
namespace root
{
TH1* h[MaxHist];
enum eDetHid { PadHid = 100000 };
}

//_____________________________________________________________________________
Bool_t
UserEventDisplayShs::ProcessingBegin()
{
  gEvDisp.Reset();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplayShs::ProcessingNormal()
{
  // const Int_t run_number = gUnpacker.get_root()->get_run_number();
  // const Int_t event_number = gUnpacker.get_event_number();

  rawData->DecodeTPCHits();

  //________________________________________________________
  //___ TPCRawHit
  Int_t npadTpc = 0;
  for (Int_t layer=0; layer<NumOfLayersTPC; ++layer) {
    const auto hc = rawData->GetTPCRawHC(layer);
    const auto nhit = hc.size();
    npadTpc += nhit;
    for (const auto& rhit : hc) {
      Int_t layer = rhit->LayerId();
      Int_t row = rhit->RowId();
      auto mean = rhit->Mean();
      auto max_adc = rhit->MaxAdc();
      // auto rms = rhit->RMS();
      auto loc_max = rhit->LocMax();
      if(loc_max < 25 || 175 <loc_max)
        continue;
      TVector3 pos = tpc::getPosition(layer, row);
      pos.SetY((loc_max - 76.75)*80.0*0.05);
      gEvDisp.FillTPCADC(layer, row, max_adc - mean);
      gEvDisp.FillTPCTDC(layer, row, loc_max);
      gEvDisp.SetTPCMarker(pos);
    }
  }

  gEvDisp.Update();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEventDisplayShs::ProcessingEnd()
{
  // gEvDisp.GetCommand();
  gEvDisp.EndOfEvent();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new UserEventDisplayShs;
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
    (InitializeParameter<DCGeomMan>("DCGEO") &&
      InitializeParameter<EventDisplayShs>() &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
