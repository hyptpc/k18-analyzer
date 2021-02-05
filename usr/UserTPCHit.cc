// -*- C++ -*-

#include <iostream>

#include <TMath.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "TPCRawHit.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"
#include "VEvent.hh"

namespace
{
using namespace root;
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser     = UserParamMan::GetInstance();
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
private:
  RawData *rawData;

public:
  UserEvent( void );
  ~UserEvent( void );
  Bool_t ProcessingBegin( void );
  Bool_t ProcessingEnd( void );
  Bool_t ProcessingNormal( void );
  Bool_t InitializeHistograms( void );
  void   InitializeEvent( void );
};

//_____________________________________________________________________________
UserEvent::UserEvent( void )
  : VEvent(),
    rawData( nullptr )
{
}

//_____________________________________________________________________________
UserEvent::~UserEvent( void )
{
  delete rawData;
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
  Int_t evnum;
  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];
};

//_____________________________________________________________________________
namespace root
{
Event event;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________________
Bool_t
UserEvent::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//_____________________________________________________________________________
Bool_t
UserEvent::ProcessingNormal( void )
{
  const Int_t run_number   = gUnpacker.get_root()->get_run_number();
  const Int_t event_number = gUnpacker.get_event_number();

  rawData = new RawData;
  rawData->DecodeTPCHits();

  HF1( 1, 0 );

  // Int_t evnum = gRM.EventNumber();
  event.runnum = run_number;
  event.evnum  = event_number;

  Int_t multiplicity = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = rawData->GetTPCRawHC( layer );
    const auto nhit = hc.size();
    multiplicity += nhit;
    for( const auto& rhit : hc ){
      // auto row = rhit->RowId();
      auto fadc = rhit->Fadc();
      const auto mhit = fadc.size();
      // for( const auto& adc : fadc ){
      // }
      auto max_adc = TMath::MaxElement( mhit, fadc.data() );
      auto rms     = TMath::RMS( mhit, fadc.data() );
      auto loc_max = TMath::LocMax( mhit, fadc.data() );
      HF1( 11, max_adc );
      HF1( 12, rms );
      HF1( 13, loc_max );
    }
  }

  HF1( 10, multiplicity );

  HF1( 1, 19 );

  return true;
}

//_____________________________________________________________________________
Bool_t
UserEvent::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
void
UserEvent::InitializeEvent( void )
{
  event.runnum = 0;
  event.evnum  = 0;

  for( Int_t it=0; it<NumOfSegTrig; ++it ){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }
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
  const Int_t    NbinAdc = 4096;
  const Double_t MinAdc  =    0.;
  const Double_t MaxAdc  = 4096.;
  const Int_t    NbinRms = 1000;
  const Double_t MinRms  =    0.;
  const Double_t MaxRms  = 1000.;

  HB1( 1, "Status", 20, 0., 20. );
  HB1( 10, "TPC Multiplicity", NumOfPadTPC+1, 0, NumOfPadTPC+1 );
  HB1( 11, "TPC Max ADC", NbinAdc, MinAdc, MaxAdc );
  HB1( 12, "TPC FADC RMS", NbinRms, MinRms, MaxRms );
  HB1( 13, "TPC LocMax ADC", NumOfTimeBucket+1, 0, NumOfTimeBucket+1 );

  // Tree
  HBTree( "tpc", "tree of TPCHit" );
  tree->Branch("runnum",   &event.runnum,   "runnum/I");
  tree->Branch("evnum",    &event.evnum,    "evnum/I");
  tree->Branch("trigpat",   event.trigpat,  Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",  event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")    &&
      InitializeParameter<UserParamMan>("USER")  );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
