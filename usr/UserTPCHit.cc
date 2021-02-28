// -*- C++ -*-

#include <iostream>

#include <TMath.h>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
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

//#define GateCalib 1
#define GateCalib 1
#define GainCalib 1
//#define GainCalib 0



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
    rawData( new RawData ),
    DCAna( new DCAnalyzer )
{
}

//_____________________________________________________________________________
UserEvent::~UserEvent( void )
{
  if( rawData ) delete rawData;
  if( DCAna ) delete DCAna;
}

//_____________________________________________________________________________
struct Event
{
  Int_t                 runnum;
  Int_t                 evnum;
  std::vector<Int_t>    trigpat;
  std::vector<Int_t>    trigflag;
  Int_t                 npadTpc;   // number of pads
  Int_t                 nhTpc;     // number of hits
  // vector (size=nhTpc)
  std::vector<Int_t>    layerTpc;  // layer id
  std::vector<Int_t>    rowTpc;    // row id
  std::vector<Int_t>    padTpc;    // pad id
  std::vector<Double_t> pedTpc;    // pedestal
  std::vector<Double_t> rmsTpc;    // rms
  std::vector<Double_t> deTpc;     // dE
  std::vector<Double_t> sigmaTpc;     // sigma
  std::vector<Double_t> tTpc;      // time
  std::vector<Double_t> chisqrTpc; // chi^2 of signal fitting
  std::vector<Double_t> cdeTpc;    // dE
  std::vector<Double_t> ctTpc;     // time
  std::vector<Double_t> dlTpc;     // time

  void clear( void )
  {
    runnum  = 0;
    evnum   = 0;
    npadTpc = 0;
    nhTpc   = 0;
    trigpat.clear();
    trigflag.clear();
    layerTpc.clear();
    rowTpc.clear();
    padTpc.clear();
    pedTpc.clear();
    rmsTpc.clear();
    deTpc.clear();
    sigmaTpc.clear();
    tTpc.clear();
    chisqrTpc.clear();
    cdeTpc.clear();
    ctTpc.clear();
    dlTpc.clear();
  }
};

//_____________________________________________________________________________
namespace root
{
  Event event;
  TH1   *h[MaxHist];
  TTree *tree;
  enum eDetHid {
    PadHid    = 100000,
  };
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

  rawData->DecodeTPCHits();
  //rawData->RecalcTPCHits();

  HF1( 1, 0 );

  event.runnum = run_number;
  event.evnum  = event_number;

  //________________________________________________________
  //___ TPCRawHit
  Int_t npadTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = rawData->GetTPCRawHC( layer );
    const auto nhit = hc.size();
    npadTpc += nhit;
    for( const auto& rhit : hc ){
      auto mean    = rhit->Mean();
      auto max_adc = rhit->MaxAdc();
      auto rms     = rhit->RMS();
      auto loc_max = rhit->LocMax();
      HF1( 11, mean );
      HF1( 12, max_adc );
      HF1( 13, rms );
      HF1( 14, loc_max );
    }
  }

  event.npadTpc = npadTpc;
  HF1( 10, npadTpc );

  HF1( 1, 19 );

  //________________________________________________________
  //___ TPCHit
  DCAna->DecodeTPCHits( rawData );
  Int_t nhTpc = 0;
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = DCAna->GetTPCHC( layer );
    for( const auto& hit : hc ){
      if( !hit || !hit->IsGood() )
        continue;
      // hit->Print();
      //Int_t layer = hit->GetLayer();
      Int_t row = hit->GetWire();
      Int_t pad = tpc::GetPadId( layer, row );
      Double_t ped = hit->GetPedestal();
      Double_t rms = hit->GetRMS();
      HF1( 21, ped );
      HF1( 23, rms );
      Int_t nhit = hit->GetNHits();
      Bool_t good_for_analysis = false;
      for( Int_t i=0; i<nhit; ++i ){
        Double_t de = hit->GetDe( i );
        Double_t time = hit->GetTime( i );
        Double_t chisqr = hit->GetChisqr( i );
        Double_t cde = hit->GetCDe( i );
        Double_t ctime = hit->GetCTime( i );
        Double_t dl = hit->GetDriftLength( i );
        Double_t sigma = hit->GetSigma( i );
        event.layerTpc.push_back( layer );
        event.rowTpc.push_back( row );
        event.padTpc.push_back( pad );
        event.pedTpc.push_back( ped );
        event.rmsTpc.push_back( rms );
        event.deTpc.push_back( de );
        event.tTpc.push_back( time );
        event.chisqrTpc.push_back( chisqr );
        event.cdeTpc.push_back( cde );
        event.ctTpc.push_back( ctime );
        event.dlTpc.push_back( dl );
	event.sigmaTpc.push_back( sigma );
	if(69.<time&&time<85.&&nhit==1)
	  HF1( 22, de );
        HF1( 24, time );
        HF1( 25, chisqr );
        HF1( 26, cde );
        HF1( 27, ctime );
        HF1( 28, dl );
        HF1( 29, sigma);

#if GateCalib
	HF1(PadHid + layer*1000 + row, time);
#endif

#if GainCalib
	if(69.<time&&time<85.&&nhit==1)
	HF1(2*PadHid + layer*1000 + row, de);
#endif

        good_for_analysis = true;
        ++nhTpc;
      }
      if( good_for_analysis ){
        // auto fadc = hit->GetRawHit()->Fadc();
        // for( Int_t tb=0, ntb=fadc.size(); tb<ntb; ++tb ){
        //   HF2( 101, tb, fadc.at( tb ) );
        // }
      }
    }
  }

  event.nhTpc = nhTpc;
  HF1( 20, nhTpc );

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
  event.clear();
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
  const Int_t    NbinDe = 1000;
  const Double_t MinDe  =    0.;
  const Double_t MaxDe  = 1000.;
  const Int_t    NbinChisqr = 1000;
  const Double_t MinChisqr  =    0.;
  const Double_t MaxChisqr  = 1000.;
  const Int_t    NbinTime = 1000;
  const Double_t MinTime  = -8000.;
  const Double_t MaxTime  =  8000.;
  const Int_t    NbinDL = 800;
  const Double_t MinDL  = -400.;
  const Double_t MaxDL  =  400.;
  const Int_t    NbinSigma = 500;
  const Double_t MinSigma  = 0.;
  const Double_t MaxSigma  =  50.;

  HB1( 1, "Status", 20, 0., 20. );
  HB1( 10, "TPC Multiplicity (Raw)", NumOfPadTPC+1, 0, NumOfPadTPC+1 );
  HB1( 11, "TPC FADC Mean", NbinAdc, MinAdc, MaxAdc );
  HB1( 12, "TPC FADC Max", NbinAdc, MinAdc, MaxAdc );
  HB1( 13, "TPC FADC RMS", NbinRms, MinRms, MaxRms );
  HB1( 14, "TPC FADC LocMax", NumOfTimeBucket+1, 0, NumOfTimeBucket+1 );
  HB1( 20, "TPC Multiplicity (TPCHit)", NumOfPadTPC+1, 0, NumOfPadTPC+1 );
  HB1( 21, "TPC Pedestal", NbinAdc, MinAdc, MaxAdc );
  HB1( 22, "TPC DeltaE", NbinDe, MinDe, MaxDe );
  HB1( 23, "TPC RMS", NbinRms, MinRms, MaxRms );
  //  HB1( 24, "TPC Time", NumOfTimeBucket+1, 0, NumOfTimeBucket+1 );
  HB1( 24, "TPC Time", (NumOfTimeBucket+1)*30, 0, NumOfTimeBucket+1 );
  HB1( 25, "TPC Chisqr", NbinChisqr, MinChisqr, MaxChisqr );
  HB1( 26, "TPC CDeltaE", NbinDe, MinDe, MaxDe );
  HB1( 27, "TPC CTime", NbinTime, MinTime, MaxTime );
  HB1( 28, "TPC DriftLength", NbinDL, MinDL, MaxDL );
  HB1( 29, "TPC sigma", NbinSigma, MinSigma, MaxSigma );
  // HB2( 101, "TPC Waveform (good)",
  //      NumOfTimeBucket+1, 0, NumOfTimeBucket+1, NbinAdc, MinAdc, MaxAdc );

#if GateCalib
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for( Int_t r=0; r<NumOfRow; ++r ){
      HB1(PadHid + layer*1000 + r , "TPC Time", (NumOfTimeBucket+1)*30, 0, NumOfTimeBucket+1 );
    }
    }
#endif

#if GainCalib
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
    for( Int_t r=0; r<NumOfRow; ++r ){
      HB1(2*PadHid + layer*1000 + r , "TPC DeltaE", NbinDe, MinDe, MaxDe );
    }
    }
#endif

  // Tree
  HBTree( "tpc", "tree of TPCHit" );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );
  tree->Branch( "trigpat", &event.trigpat );
  tree->Branch( "trigflag", &event.trigflag );
  tree->Branch( "npadTpc", &event.npadTpc );
  tree->Branch( "nhTpc", &event.nhTpc );
  tree->Branch( "layerTpc", &event.layerTpc );
  tree->Branch( "rowTpc", &event.rowTpc );
  tree->Branch( "padTpc", &event.padTpc );
  tree->Branch( "pedTpc", &event.pedTpc );
  tree->Branch( "rmsTpc", &event.rmsTpc );
  tree->Branch( "deTpc", &event.deTpc );
  tree->Branch( "tTpc", &event.tTpc );
  tree->Branch( "chisqrTpc", &event.chisqrTpc );
  tree->Branch( "cdeTpc", &event.cdeTpc );
  tree->Branch( "ctTpc", &event.ctTpc );
  tree->Branch( "dlTpc", &event.dlTpc );
  tree->Branch( "sigmaTpc", &event.sigmaTpc );

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO") &&
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
