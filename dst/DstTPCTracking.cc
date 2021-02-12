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
  const std::string& class_name("DstTPCTracking");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
  debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();
  const Int_t MaxTPCHits = 10000;
  const Int_t MaxTPCTracks = 100;
  const Int_t MaxTPCnHits = 50;
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
  TBranch *tmp_branch =0;
}

//_____________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  std::vector<Int_t>    trigpat;
  std::vector<Int_t>    trigflag;
  int status;
  int nhittpc;
  Int_t nttpc;                   // Number of Tracks
  Int_t nhit_track[MaxTPCTracks]; // Number of Hits (in 1 tracks)
  Double_t chisqr[MaxTPCTracks];
  Double_t x0[MaxTPCTracks];
  Double_t y0[MaxTPCTracks];
  Double_t u0[MaxTPCTracks];
  Double_t v0[MaxTPCTracks];
  Double_t theta[MaxTPCTracks];
  Int_t hitlayer[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t residual[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_x[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_y[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_z[MaxTPCTracks][MaxTPCnHits];
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
    tTpc.clear();
    deTpc.clear();
    chisqrTpc.clear();
  }
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
  event.runnum = 0;
  event.evnum = 0;
  event.status = 0;
  event.nhittpc = 0;
  event.nttpc = 0;
  event.trigpat.clear();
  event.trigflag.clear();
  for(int i=0; i<MaxTPCTracks; ++i){
    event.nhit_track[i] =0;
    event.chisqr[i] =-9999.;
    event.x0[i] =-9999.;
    event.y0[i] =-9999.;
    event.u0[i] =-9999.;
    event.v0[i] =-9999.;
    event.theta[i] =-9999.;
    for(int j=0; j<MaxTPCnHits; ++j){
      event.hitlayer[i][j] =-999;
      event.hitpos_x[i][j] =-9999.;
      event.hitpos_y[i][j] =-9999.;
      event.hitpos_z[i][j] =-9999.;
      event.calpos_x[i][j] =-9999.;
      event.calpos_y[i][j] =-9999.;
      event.calpos_z[i][j] =-9999.;
      event.residual[i][j] =-9999.;
      event.residual_x[i][j] =-9999.;
      event.residual_y[i][j] =-9999.;
      event.residual_z[i][j] =-9999.;
    }
  }
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
  GetEntry(ievent);
  // Long64_t tentry = TTreeCont[kTpcHit]->LoadTree(ievent);
  // tmp_branch->GetEntry(tentry);
  
  event.runnum = src.runnum;
  event.evnum = src.evnum;
  event.nhittpc = src.nhTpc;
  //event.trigpat = src.trigpat;
  //event.trigflag = src.trigflag;
  
 

  HF1( 1, event.status++ );
  DCAnalyzer *DCAna = new DCAnalyzer();
  
  std::cout<<src.nhTpc<<", "<<src.padTpc.size()<<std::endl;
    //DCAna->RecalcTPCHits(src.nhTpc, src.padTpc, src.tTpc, src.deTpc);
  //DCAna->TrackSearchTPC();



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


  HBTree( "tpc", "tree of DstTPC" );

  tree->Branch("status", &event.status, "status/I" );
  tree->Branch("evnum", &event.evnum, "evnum/I" );
  //tree->Branch("trigpat", &event.trigpat);
  //tree->Branch("trigflag", &event.trigflag);
  tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");
  tree->Branch("nhit_track",event.nhit_track,"nhit_track[nttpc]/I");
  tree->Branch("chisqr",event.chisqr,"chisqr[nttpc]/D");
  tree->Branch("x0",event.x0,"x0[nttpc]/D");
  tree->Branch("y0",event.y0,"y0[nttpc]/D");
  tree->Branch("u0",event.u0,"u0[nttpc]/D");
  tree->Branch("v0",event.v0,"v0[nttpc]/D");
  tree->Branch("theta",event.theta,"theta[nttpc]/D");

  tree->Branch("hitlayer",event.hitlayer,"hitlayer[nttpc][64]/I");
  tree->Branch("hitpos_x",event.hitpos_x,"hitpos_x[nttpc][64]/D");
  tree->Branch("hitpos_y",event.hitpos_y,"hitpos_y[nttpc][64]/D");
  tree->Branch("hitpos_z",event.hitpos_z,"hitpos_z[nttpc][64]/D");
  tree->Branch("calpos_x",event.calpos_x,"calpos_x[nttpc][64]/D");
  tree->Branch("calpos_y",event.calpos_y,"calpos_y[nttpc][64]/D");
  tree->Branch("calpos_z",event.calpos_z,"calpos_z[nttpc][64]/D");
  tree->Branch("residual",event.residual,"residual[nttpc][64]/D");
  tree->Branch("residual_x",event.residual_x,"residual_x[nttpc][64]/D");
  tree->Branch("residual_y",event.residual_y,"residual_y[nttpc][64]/D");
  tree->Branch("residual_z",event.residual_z,"residual_z[nttpc][64]/D");

  TTreeCont[kTpcHit]->SetBranchStatus("*", 0);
  TTreeCont[kTpcHit]->SetBranchStatus("runnum",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("evnum",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("trigpat",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("trigflag",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("npadTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("nhTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("layerTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("rowTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("padTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("pedTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("rmsTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("deTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("tTpc",  1);
  TTreeCont[kTpcHit]->SetBranchStatus("chisqrTpc",  1);
  

  TTreeCont[kTpcHit]->SetBranchAddress("runnum",  &src.runnum);
  TTreeCont[kTpcHit]->SetBranchAddress("evnum",  &src.evnum);
  TTreeCont[kTpcHit]->SetBranchAddress("trigpat",  &src.trigpat);
  TTreeCont[kTpcHit]->SetBranchAddress("trigflag", &src.trigflag);
  TTreeCont[kTpcHit]->SetBranchAddress("npadTpc",  &src.npadTpc, &tmp_branch);
  TTreeCont[kTpcHit]->SetBranchAddress("nhTpc",  &src.nhTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("layerTpc", &src.layerTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("rowTpc",  &src.rowTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("padTpc",  &src.padTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("pedTpc",  &src.pedTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("rmsTpc",  &src.rmsTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("deTpc",  &src.deTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("tTpc",  &src.tTpc);
  TTreeCont[kTpcHit]->SetBranchAddress("chisqrTpc",  &src.chisqrTpc);


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
