// -*- C++ -*-

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
#include "VEvent.hh"
#include "DetectorID.hh"

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
auto& gEvDisp = EventDisplayShs::GetInstance();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
VEvent::VEvent( void )
{
}

//_____________________________________________________________________________
VEvent::~VEvent( void )
{
}

//_____________________________________________________________________________
class UserEvent : public VEvent
{
public:
  UserEvent( void );
  ~UserEvent( void );
  Bool_t ProcessingBegin( void );
  Bool_t ProcessingEnd( void );
  Bool_t ProcessingNormal( void );
  Bool_t InitializeHistograms( void );
  void   InitializeEvent( void );

private:
  RawData*    rawData;
  DCAnalyzer* DCAna;
};

//_____________________________________________________________________________
UserEvent::UserEvent( void )
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer)
{
}

//_____________________________________________________________________________
UserEvent::~UserEvent( void )
{
  if( rawData ) delete rawData;
  if( DCAna ) delete DCAna;
}

//_____________________________________________________________________________
namespace root
{
  TH1* h[MaxHist];
  enum eDetHid
  {
    PadHid = 100000
  };
}

//_____________________________________________________________________________
Bool_t
UserEvent::ProcessingBegin( void )
{
  gEvDisp.Reset();
  InitializeEvent();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEvent::ProcessingNormal( void )
{
  const Int_t run_number = gUnpacker.get_root()->get_run_number();
  const Int_t event_number = gUnpacker.get_event_number();

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
      gEvDisp.FillTPCADC(layer, row, max_adc - mean);
      gEvDisp.FillTPCTDC(layer, row, loc_max);
    }
  }

  gEvDisp.Update();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEvent::ProcessingEnd( void )
{
  gEvDisp.GetCommand();
  return true;
}

//_____________________________________________________________________________
void
UserEvent::InitializeEvent( void )
{
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new UserEvent;
}

//_____________________________________________________________________________
Bool_t
ConfMan:: InitializeHistograms( void )
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO") &&
      InitializeParameter<EventDisplayShs>() &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<UserParamMan>("USER") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
