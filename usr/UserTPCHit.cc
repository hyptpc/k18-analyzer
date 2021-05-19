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
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"


//#define GateCalib 1
#define GateCalib 0
#define GainCalib 0
#define Srdata 0
//#define GainCalib 0



namespace
{
using namespace root;
const auto qnan = TMath::QuietNaN();
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
  HodoAnalyzer* hodoAna;
  DCAnalyzer* DCAna;
};

//_____________________________________________________________________________
UserEvent::UserEvent( void )
  : VEvent(),
    rawData( new RawData ),
    hodoAna(new HodoAnalyzer),
    DCAna( new DCAnalyzer )
{
}

//_____________________________________________________________________________
UserEvent::~UserEvent( void )
{
  if( rawData ) delete rawData;
  if(hodoAna) delete hodoAna;
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


  Int_t htofnhits;
  Int_t htofhitpat[MaxHits];
  Double_t htofua[NumOfSegHTOF];
  Double_t htofda[NumOfSegHTOF];
  Double_t htofut[NumOfSegHTOF][MaxDepth];
  Double_t htofdt[NumOfSegHTOF][MaxDepth];

  Double_t htofmt[NumOfSegHTOF][MaxDepth];
  Double_t htofde[NumOfSegHTOF];

  Int_t    nhHtof;
  Int_t    csHtof[NumOfSegHTOF*MaxDepth];
  Double_t HtofSeg[NumOfSegHTOF*MaxDepth];
  Double_t tHtof[NumOfSegHTOF*MaxDepth];
  Double_t dtHtof[NumOfSegHTOF*MaxDepth];
  Double_t deHtof[NumOfSegHTOF*MaxDepth];

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
    
    htofnhits =0;
    nhHtof =0;
    for(Int_t it=0; it<MaxHits; ++it){
      htofhitpat[it]  = -1;
    }
    for(Int_t it=0; it<NumOfSegHTOF; it++){
      htofua[it] = qnan;
      htofda[it] = qnan;
      htofde[it] = qnan;
      for(Int_t m=0; m<MaxDepth; ++m){
	htofut[it][m] = qnan;
	htofdt[it][m] = qnan;
	htofmt[it][m] = qnan;
	
	csHtof[MaxDepth*it + m]  = 0;
	HtofSeg[MaxDepth*it + m] = qnan;
	tHtof[MaxDepth*it + m]   = qnan;
	dtHtof[MaxDepth*it + m]  = qnan;
	deHtof[MaxDepth*it + m]  = qnan;
      }
    }

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

  static const auto MinTdcHTOF = gUser.GetParameter("TdcHTOF", 0);
  static const auto MaxTdcHTOF = gUser.GetParameter("TdcHTOF", 1);

  rawData->DecodeHits();
  
  ///// HTOF Raw data
  {
    Int_t htof_nhits = 0;
    const HodoRHitContainer& cont = rawData->GetHTOFRawHC();
    Int_t nh = cont.size();
    //    HF1(HTOFHid, nh);
    Int_t nh1 = 0, nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      HodoRawHit* hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      //HF1(HTOFHid+1, seg-0.5);
      // Up
      Int_t Au = hit->GetAdcUp();
      //HF1(HTOFHid+100*seg+1, Au);
      event.htofua[seg-1] = Au;
      Bool_t is_hit_u = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcUp(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcUp(m);
        //HF1(HTOFHid +100*seg +3, T);
        event.htofut[seg-1][m] = T;
        if(MinTdcHTOF < T && T < MaxTdcHTOF) is_hit_u = true;
      }
      // if(is_hit_u) HF1(HTOFHid+100*seg+5, Au);
      // else         HF1(HTOFHid+100*seg+7, Au);
      // Down
      Int_t Ad = hit->GetAdcDown();
      //      HF1(HTOFHid+100*seg+2, Ad);
      event.htofda[seg-1] = Ad;
      Bool_t is_hit_d = false;
      for(Int_t m=0, n_mhit=hit->GetSizeTdcDown(); m<n_mhit; ++m){
        Int_t T = hit->GetTdcDown(m);
        //HF1(HTOFHid +100*seg +4, T);
        event.htofdt[seg-1][m] = T;
        if(MinTdcHTOF < T && T < MaxTdcHTOF) is_hit_d = true;
      }
      // if(is_hit_d) HF1(HTOFHid+100*seg+6, Ad);
      // else         HF1(HTOFHid+100*seg+8, Ad);
      // HitPat
      if(is_hit_u || is_hit_d){
        ++nh1; //HF1(HTOFHid+3, seg-0.5);
      }
      if(is_hit_u && is_hit_d){
        event.htofhitpat[htof_nhits++] = seg;
        ++nh2; //HF1(HTOFHid+5, seg-0.5);
      }
    }
    //HF1(HTOFHid+2, nh1); HF1(HTOFHid+4, nh2);
    event.htofnhits = htof_nhits;
  }

  ///// HTOF Normalized data
  hodoAna->DecodeHTOFHits(rawData);
  {
    Int_t nh = hodoAna->GetNHitsHTOF();
    //HF1(HTOFHid+10, Double_t(nh));
    Int_t nh2 = 0;
    for(Int_t i=0; i<nh; ++i){
      Hodo2Hit *hit = hodoAna->GetHitHTOF(i);
      if(!hit) continue;
      Int_t seg = hit->SegmentId()+1;
      Int_t n_mhit = hit->GetNumOfHit();
      for(Int_t m=0; m<n_mhit; ++m){
	//HF1(HTOFHid+11, seg-0.5);
	Double_t au = hit->GetAUp(), ad = hit->GetADown();
	Double_t tu = hit->GetTUp(), td = hit->GetTDown();
	Double_t ctu = hit->GetCTUp(), ctd = hit->GetCTDown();
	Double_t mt = hit->MeanTime(), cmt = hit->CMeanTime();
	Double_t de = hit->DeltaE();
	event.htofmt[seg-1][m] = mt;
	event.htofde[seg-1] = de;
	// HF1(HTOFHid+100*seg+11, tu); HF1(HTOFHid+100*seg+12, td);
	// HF1(HTOFHid+100*seg+13, mt);
	// HF1(HTOFHid+100*seg+17, ctu); HF1(HTOFHid+100*seg+18, ctd);
	// HF1(HTOFHid+100*seg+19, cmt); HF1(HTOFHid+100*seg+20, ctu-ctd);
	// HF2(HTOFHid+100*seg+21, tu, au); HF2(HTOFHid+100*seg+22, td, ad);
	// HF2(HTOFHid+100*seg+23, ctu, au); HF2(HTOFHid+100*seg+24, ctd, ad);
	// HF1(HTOFHid+12, cmt);
	// if(m == 0){
	//   HF1(HTOFHid+100*seg+14, au); HF1(HTOFHid+100*seg+15, ad);
	//   HF1(HTOFHid+100*seg+16, de); HF1(HTOFHid+13, de);
	// }
	if(de > 0.5){
	  //HF1(HTOFHid+15, seg-0.5);
	  ++nh2;
	}
      }
    }
    Int_t nc = hodoAna->GetNClustersHTOF();
    //HF1(HTOFHid+30, Double_t(nc));
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cluster = hodoAna->GetClusterHTOF(i);
      if(!cluster) continue;
      Int_t cs = cluster->ClusterSize();
      Double_t ms = cluster->MeanSeg()+1;
      Double_t cmt = cluster->CMeanTime();
      Double_t de = cluster->DeltaE();
      // HF1(HTOFHid+31, Double_t(cs));
      // HF1(HTOFHid+32, ms-0.5);
      // HF1(HTOFHid+33, cmt); HF1(HTOFHid+34, de);
    }
  }

  {
    Int_t nc = hodoAna->GetNClustersHTOF();
    event.nhHtof = nc;
    for(Int_t i=0; i<nc; ++i){
      HodoCluster *cl = hodoAna->GetClusterHTOF(i);
      if(!cl) continue;
      event.csHtof[i] = cl->ClusterSize();
      event.HtofSeg[i] = cl->MeanSeg()+1;
      event.tHtof[i] = cl->CMeanTime();
      event.dtHtof[i] = cl->TimeDif();
      event.deHtof[i] = cl->DeltaE();
    }
  }



  rawData->DecodeTPCHits();
  rawData->RecalcTPCHits();
  
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
      auto min_adc = rhit->MinAdc();
      auto rms     = rhit->RMS();
      auto loc_max = rhit->LocMax();
      HF1( 11, mean );
      HF1( 12, max_adc );
      HF1( 13, rms );
      HF1( 14, loc_max );
      HF1( 15, min_adc );
    }
  }

  //________________________________________________________
  //___ TPCRawHit after subtraction
  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
    auto hc = rawData->GetTPCCorHC( layer );
    for( const auto& rhit : hc ){
      auto mean    = rhit->Mean();
      auto max_adc = rhit->MaxAdc();
      auto min_adc = rhit->MinAdc();
      auto rms     = rhit->RMS();
      auto loc_max = rhit->LocMax();
      HF1( 31, mean );
      HF1( 32, max_adc );
      HF1( 33, rms );
      HF1( 34, loc_max );
      HF1( 35, min_adc );
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
	//	if(69.<time&&time<85.&&nhit==1)
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
	//	if(69.<time&&time<85.&&nhit==1)
	HF1(2*PadHid + layer*1000 + row, de);
#endif

#if Srdata
	if(69.<time&&time<85.&&nhit==1&&layer==31)
	  HF1(3*PadHid + layer*1000 + row, de);
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
  HB1( 15, "TPC FADC Min", NbinAdc, MinAdc, MaxAdc );
  HB1( 31, "TPC FADC Mean Cor", NbinAdc, MinAdc, MaxAdc );
  HB1( 32, "TPC FADC Max Cor", NbinAdc, MinAdc, MaxAdc );
  HB1( 33, "TPC FADC RMS Cor", NbinRms, MinRms, MaxRms );
  HB1( 34, "TPC FADC LocMax Cor", NumOfTimeBucket+1, 0, NumOfTimeBucket+1 );
  HB1( 35, "TPC FADC Min Cor", NbinAdc, MinAdc, MaxAdc );

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

#if Srdata
  //  for( Int_t layer=0; layer<NumOfLayersTPC; ++layer ){
  Int_t layer= 31;
  const Int_t NumOfRow = tpc::padParameter[layer][tpc::kNumOfPad];
  for( Int_t r=0; r<NumOfRow; ++r ){
    HB1(3*PadHid + layer*1000 + r , "TPC DeltaE", NbinDe, MinDe, MaxDe );
  }
  //}
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
  
  //htof
  tree->Branch("htofnhits",   &event.htofnhits,   "htofnhits/I");
  tree->Branch("htofhitpat",   event.htofhitpat,  Form("htofhitpat[%d]/I", NumOfSegHTOF));
  tree->Branch("htofua",       event.htofua,      Form("htofua[%d]/D", NumOfSegHTOF));
  tree->Branch("htofda",       event.htofda,      Form("htofda[%d]/D", NumOfSegHTOF));
  tree->Branch("htofut",       event.htofut,      Form("htofut[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofdt",       event.htofdt,      Form("htofdt[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofmt",   event.htofmt,   Form("htofmt[%d][%d]/D", NumOfSegHTOF, MaxDepth));
  tree->Branch("htofde",   event.htofde,   Form("htofde[%d]/D", NumOfSegHTOF));
  tree->Branch("nhHtof",  &event.nhHtof,     "nhHtof/I");
  tree->Branch("csHtof", event.csHtof, "csHtof[nhHtof]/I");
  tree->Branch("HtofSeg",event.HtofSeg, "HtofSeg[nhHtof]/D");
  tree->Branch("tHtof", event.tHtof,  "tHtof[nhHtof]/D");
  tree->Branch("dtHtof", event.dtHtof, "dtHtof[nhHtof]/D");
  tree->Branch("deHtof", event.deHtof, "deHtof[nhHtof]/D");



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
