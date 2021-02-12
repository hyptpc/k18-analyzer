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


#include "DstHelper.hh"
#include "DebugCounter.hh"

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
  debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();
  //  const int NumOfCobos = 8;
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kTpcHit,  kOutFile, nArgc
    };
  std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]",
    "[TPCHit]",  "[OutFile]" };
  std::vector<TString> TreeName =
  { "", "",
    "tpc", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
}

//_____________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  int status;
  int nhittpc;
  Int_t nttpc;                   // Number of Tracks
  std::vector<Int_t>    nhit_track; // Number of Hits (in 1 tracks)
  std::vector<Double_t> chisqr;
  std::vector<Double_t> x0;
  std::vector<Double_t> y0;
  std::vector<Double_t> u0;
  std::vector<Double_t> v0;
  std::vector<Double_t> theta;
  std::vector<Int_t>    hitlayer;
  std::vector<Double_t> hitpos_x;
  std::vector<Double_t> hitpos_y;
  std::vector<Double_t> hitpos_z;
  std::vector<Double_t> calpos_x;
  std::vector<Double_t> calpos_y;
  std::vector<Double_t> calpos_z;
  std::vector<Double_t> residual;
  std::vector<Double_t> residual_x;
  std::vector<Double_t> residual_y;
  std::vector<Double_t> residual_z;

};

//_____________________________________________________________________
struct Src
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
  std::vector<Double_t> tTpc;      // time
  std::vector<Double_t> chisqrTpc; // chi^2 of signal fitting
};

namespace root
{
  Event  event;
  Src    src;
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
    gCounter.check();
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
  runnum = 0;
  evnum = 0;
  status = 0;
  nhittpc = 0;
  nttpc = 0;
  nhit_track.clear();
  chisqr.clear();
  x0.clear();
  y0.clear();
  u0.clear();
  v0.clear();
  theta.clear();   
  hitlayer.clear();
  hitpos_x.clear();
  hitpos_y.clear();
  hitpos_z.clear();
  calpos_x.clear();
  calpos_y.clear();
  calpos_z.clear();
  residual.clear();
  residual_x.clear();
  residual_y.clear();
  residual_z.clear();
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


  //  if( ievent%10000==0 ){
  if( ievent%100==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  event.runnum = src.runnum;
  event.evnum = src.evnum;
  event.nhittpc = src.nhTpc;
  event.trigpat = src.trigpat;
  event.trigflag = src.trigflag;
  
  GetEntry(ievent);

  HF1( 1, event.status++ );
  DCAnalyzer *DCAna = new DCAnalyzer();

  DCAna->RecalcTPCHits( src.nhTpc, 
			src.padTpc, src.tTpc, src.deTpc);
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
