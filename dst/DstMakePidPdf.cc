// -*- C++ -*-

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <memory> 

#include <filesystem_util.hh>
#include <UnpackerManager.hh>
#include "TF1.h"

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "FieldMan.hh"
#include "DatabasePDG.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "TPCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DstHelper.hh"
#include "HodoPHCMan.hh"
#include "PidPdfMan.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertex.hh"
#include "TPCLTrackHit.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "UserParamMan.hh"

//#include "HypTPCFitter.hh"
//#include "HypTPCTask.hh"

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
const auto& gUser = UserParamMan::GetInstance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
auto&       gPidPdf = PidPdfMan::GetInstance();  
const auto& zHSCenter = gGeom.LocalZ("HS");
const Double_t truncatedMean = 0.8; //80%
const Double_t dist_cut = 5;
const Double_t chisqrcut = 100;
const Double_t decutppi = 2.0;
const Double_t gfposycut = 100.0;  

//For GenFit Setting
const Bool_t Const_field = false; //Must be false for linear tracking
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

const int kPmin = pidlikeli::kPmin;
const int kPmax = pidlikeli::kPmax;
const int kMomstep = pidlikeli::kDP;
const int kNmom = pidlikeli::kNmom;
const int kNpid = pidlikeli::kNpid; 
const int kNchg = pidlikeli::kNchg;
const int kNtype = pidlikeli::kNtype;
const int kNbe = pidlikeli::kNbe;
const auto& plist = pidlikeli::plist;
const auto& clist = pidlikeli::clist;
const auto& type = pidlikeli::type;  
  
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kGenfitE42, kOutFile, nArgc
};
std::vector<TString> ArgName =
  { "[Process]", "[ConfFile]", "[GenfitE42]", "[OutFile]" };
std::vector<TString> TreeName = { "", "", "tpc", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________________
struct Event
{

  Int_t status;
  Int_t runnum;
  Int_t evnum;
  void clear( void )
  {
    runnum = 0;
    evnum = 0;
    status = 0;
  }
};

//_____________________________________________________________________________
struct Src
{
  std::unique_ptr<TTreeReaderValue<Int_t>> runnum;
  std::unique_ptr<TTreeReaderValue<Int_t>> evnum;
};

namespace root
{
Event  event;
Src    src;
TH1   *h[MaxHist];
TTree *tree;
  enum eDetHid {
    PadHid    = 100000,
    genfitHid = 200000
  };
  
}

//_____________________________________________________________________________
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
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();
  
  //Initiallize Geometry, Field, Fitter
  // HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  // //Initiallize the genfit track container
  // HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  // GFtracks.SetVerbosity(verbosity);
  // std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  // std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;

#if 0
  GFtracks.DebugMode();
#endif
  Int_t ievent = skip;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    gCounter.check();
    InitializeEvent();
    if( DstRead( ievent ) ) tree->Fill();
  }
  std::cout << "#D Event Number: " << std::setw(6)
            << ievent << std::endl;
  
  gPidPdf.WriteToRootfile(TFileCont[kOutFile]);
  DstClose();

  //  delete fitter;
  return EXIT_SUCCESS;  
}

//_____________________________________________________________________________
Bool_t
dst::InitializeEvent( void )
{
  event.clear();

  return true;
}

//_____________________________________________________________________________
Bool_t
dst::DstOpen( std::vector<std::string> arg )
{
  int open_file = 0;
  int open_tree = 0;
  for( Int_t i=0; i<nArgc; ++i ){
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
//_____________________________________________________________________________
Bool_t
dst::DstRead( int ievent )
{
  return true;
  Double_t vtx_scan_range = gUser.GetParameter("VertexScanRange");
  //if( ievent%1000==0 ){
  if( ievent%100==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }
  return true;
  
}

//_____________________________________________________________________________
Bool_t
dst::DstClose( void )
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const Int_t n = TFileCont.size();
  for( Int_t i=0; i<n; ++i ){
    if( TTreeReaderCont[i] ) delete TTreeReaderCont[i];
    if( TTreeCont[i] ) delete TTreeCont[i];
    if( TFileCont[i] ) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{
  HBTree( "tpc", "tree of DstMakePdf" );
  tree->Branch( "status", &event.status );
  tree->Branch( "runnum", &event.runnum );
  tree->Branch( "evnum", &event.evnum );

  TTreeReaderCont[kGenfitE42] = new TTreeReader( "tpc", TFileCont[kGenfitE42] );
  const auto& reader = TTreeReaderCont[kGenfitE42];
  src.runnum = std::make_unique<TTreeReaderValue<Int_t>>( *reader, "runnum" );
  src.evnum = std::make_unique<TTreeReaderValue<Int_t>>( *reader, "evnum" );
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
      InitializeParameter<TPCPositionCorrector>("TPCPOS") &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<PidPdfMan>("PIDPDF") &&
      InitializeParameter<FieldMan>("FLDMAP", "HSFLDMAP") &&
      InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess( void )
{
  return true;
}
