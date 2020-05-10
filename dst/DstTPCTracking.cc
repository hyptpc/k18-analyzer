/**
 *  file: DstTPCTracking.cc
 *  date: 2020.04.11
 *
 */

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <filesystem_util.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "Kinematics.hh"
#include "DCGeomMan.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "HodoPHCMan.hh" 
#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "RawData.hh"

#include "DstHelper.hh"

namespace
{
  using namespace root;
  using namespace dst;
  const std::string& class_name("DstKuramaHodoscope");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance(); 
  const int MaxTPCHits = 10000;
  const int NumOfCobos = 8;
}

namespace dst
{
  enum kArgc
  {
    kProcess, kConfFile,
    kTpcHit_0, kTpcHit_1, kTpcHit_2, kTpcHit_3,
    kTpcHit_4, kTpcHit_5, kTpcHit_6, kTpcHit_7,
    kOutFile, nArgc
  };
  std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", 
    "[TPCHit_0]", "[TPCHit_1]", "[TPCHit_2]", "[TPCHit_3]",
    "[TPCHit_4]", "[TPCHit_5]", "[TPCHit_6]", "[TPCHit_7]",
    "[OutFile]" };
  std::vector<TString> TreeName =
  { "", "", 
    "tpchit0", "tpchit1", "tpchit2", "tpchit3", 
    "tpchit4", "tpchit5", "tpchit6", "tpchit7", 
    "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
}

//_____________________________________________________________________
struct Event
{
  int evnum;
  int status;
  int nhittpc; 
};

//_____________________________________________________________________
struct Src
{
  int evnum;
  int nhittpc;           

  int padidtpc[MaxTPCHits];
  double ytpc[MaxTPCHits];
  double chargetpc[MaxTPCHits];
};

namespace root
{
  Event  event;
  Src    src[NumOfCobos];
  TH1   *h[MaxHist];
  TTree *tree;
}

//_____________________________________________________________________
int
main( int argc, char **argv )
{
  std::vector<std::string> arg( argv, argv+argc );

  if( !CheckArg( arg ) )
    return EXIT_FAILURE;
  if( !DstOpen( arg ) )
    return EXIT_FAILURE;
  if( !gConf.Initialize( arg[kConfFile] ) )
    return EXIT_FAILURE;

  int nevent = GetEntries( TTreeCont );

  CatchSignal::Set();

  int ievent = 0;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }

  std::cout << "#D Event Number: " << std::setw(6)
    << ievent << std::endl;

  DstClose();

  return EXIT_SUCCESS;
}

//_____________________________________________________________________
bool
dst::InitializeEvent( void )
{
  event.status   = 0;
  event.evnum = 0;
  event.nhittpc = 0; 


  // for( int i=0; i<MaxTPCHits; ++i ){
  // }
  return true;
}

//_____________________________________________________________________
bool
dst::DstOpen( std::vector<std::string> arg )
{
  int open_file = 0;
  int open_tree = 0;
  for( std::size_t i=0; i<nArgc; ++i ){
    if( i==kProcess || i==kConfFile || i==kOutFile ) continue;
    open_file += OpenFile( TFileCont[i], arg[i] );
    open_tree += OpenTree( TFileCont[i], TTreeCont[i], TreeName[i] );
  }

  if( open_file!=open_tree || open_file!=nArgc-3 )
    return false;
  if( !CheckEntries( TTreeCont ) )
    return false;

  TFileCont[kOutFile] = new TFile( arg[kOutFile].c_str(), "recreate" );

  return true;
}

//_____________________________________________________________________
bool
dst::DstRead( int ievent )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");


  if( ievent%10000==0 ){
    std::cout << "#D Event Number: "
      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1( 1, event.status++ );

  RawData *rawData = new RawData();
  DCAnalyzer *DCAna = new DCAnalyzer();
  for( int coboi=0; coboi < NumOfCobos; coboi++) 
  {
    if( coboi==0 ) event.evnum = src[coboi].evnum;
    else {
      if( src[coboi].evnum != event.evnum ) {
	std::cerr << "#E: TpcHit_" << coboi << " event number does not match! " << std::endl;
	return false;
      }
    }
    event.nhittpc += src[coboi].nhittpc;

    for(int hiti=0; hiti<src[coboi].nhittpc; hiti++) 
    {
      rawData->DecodeTPCHits( src[coboi].padidtpc[hiti], 
	  		   src[coboi].ytpc[hiti], 
			   src[coboi].chargetpc[hiti] );
    }
  }

  DCAna->DecodeTPCHits( rawData );
  DCAna->TrackSearchTPC();


  // for(int i=0;i<src.nhittpc;++i){
  //   event.ititpc[i] = src.ititpc[i];
  // }

  // if( event.nhBh1<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.nhBh2<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.nhTof<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.ntKurama<=0 ) return true;
  // if( event.ntKurama>MaxTPCHits )
  //   event.ntKurama = MaxTPCHits;

  //HF1( 1, event.status++ );

  delete rawData;
  delete DCAna;

  return true;
}

//_____________________________________________________________________
bool
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const std::size_t n = TFileCont.size();
  for( std::size_t i=0; i<n; ++i ){
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 21, 0., 21. );


  HBTree( "ktpc_g", "tree of DstTPC_g" );

  tree->Branch("status", &event.status, "status/I" );  
  tree->Branch("evnum", &event.evnum, "evnum/I" );  
  tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");


  ////////// Bring Address From Dst
  for( int coboi=kTpcHit_0; coboi<kTpcHit_0 + NumOfCobos; coboi++ ) 
  {
    TTreeCont[coboi]->SetBranchStatus("*", 0);
    TTreeCont[coboi]->SetBranchStatus("evnum",  1);
    TTreeCont[coboi]->SetBranchStatus("nhittpc",  1);

    TTreeCont[coboi]->SetBranchStatus("padidtpc", 1);
    TTreeCont[coboi]->SetBranchStatus("ytpc", 1);
    TTreeCont[coboi]->SetBranchStatus("chargetpc", 1);


    TTreeCont[coboi]->SetBranchAddress("evnum", &src[coboi].evnum);
    TTreeCont[coboi]->SetBranchAddress("nhittpc", &src[coboi].nhittpc);

    TTreeCont[coboi]->SetBranchAddress("padidtpc", src[coboi].padidtpc);
    TTreeCont[coboi]->SetBranchAddress("ytpc", src[coboi].ytpc);
    TTreeCont[coboi]->SetBranchAddress("chargetpc", src[coboi].chargetpc);
  }


  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<HodoPHCMan>("HDPHC") ); 
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
