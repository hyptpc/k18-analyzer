/**
 *  file: DstTPCTracking_geant.cc
 *  date: 2020.04.02
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

namespace
{
  using namespace root;
  using namespace dst;
  const std::string& class_name("DstKuramaHodoscope");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance(); 
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kTPCGeant, kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[TPCGeant]",
      "[OutFile]" };
  std::vector<TString> TreeName =
    { "", "", "TPC_g", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
}

//_____________________________________________________________________
struct Event
{
  
  int evnum;
  int status;
  int nttpc;                 // Number of Hit in Pads
};

//_____________________________________________________________________
struct Src
{
  int evnum;
  int nttpc;                 // Number of Hit in Pads

  int ititpc[MaxHits];
  int idtpc[MaxHits];
  double xtpc[MaxHits];//with resolution
  double ytpc[MaxHits];//with resolution
  double ztpc[MaxHits];//with resolution
  double x0tpc[MaxHits];//w/o resolution
  double y0tpc[MaxHits];//w/o resolution
  double z0tpc[MaxHits];//w/o resolution
  double resoX[MaxHits];
  double pxtpc[MaxHits];
  double pytpc[MaxHits];
  double pztpc[MaxHits];
  double pptpc[MaxHits];   // total mometum
  double masstpc[MaxHits];   // mass TPC
  double betatpc[MaxHits];
  double edeptpc[MaxHits];
  double dedxtpc[MaxHits];
  double slengthtpc[MaxHits];
  int laytpc[MaxHits];
  int rowtpc[MaxHits];
  int parentID[MaxHits];
  int iPadtpc[MaxHits];//Pad number (0 origin)

  double xtpc_pad[MaxHits];//pad center
  double ytpc_pad[MaxHits];//pad center(dummy)
  double ztpc_pad[MaxHits];//pad center

  double dxtpc_pad[MaxHits];//x0tpc - xtpc
  double dytpc_pad[MaxHits];//y0tpc - ytpc = 0 (dummy)
  double dztpc_pad[MaxHits];//z0tpc - ztpc

 
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
  event.nttpc = 0; 

  
  // for( int i=0; i<MaxHits; ++i ){
  //   event.trpptpc[i]=-9999.;
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
  
  event.evnum = src.evnum;
  event.nttpc = src.nttpc;
  

  

  
  


  
  

  DCAnalyzer DCAna;
  DCAna.DecodeTPCHits(src.nttpc, src.iPadtpc, src.dxtpc_pad, src.dztpc_pad);



#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[nttpc]: "<<std::setw(2)<<src.nttpc<<" "<<std::endl;
#endif

  // for(int i=0;i<src.nttpc;++i){
  //   event.ititpc[i] = src.ititpc[i];
  // }

  // if( event.nhBh1<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.nhBh2<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.nhTof<=0 ) return true;
  // HF1( 1, event.status++ );

  // if( event.ntKurama<=0 ) return true;
  // if( event.ntKurama>MaxHits )
  //   event.ntKurama = MaxHits;

  //HF1( 1, event.status++ );

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
  tree->Branch("nttpc",&event.nttpc,"nttpc/I");


  ////////// Bring Address From Dst
  TTreeCont[kTPCGeant]->SetBranchStatus("*", 0);
  TTreeCont[kTPCGeant]->SetBranchStatus("evnum",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("nttpc",  1);

  TTreeCont[kTPCGeant]->SetBranchStatus("ititpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("idtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ztpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("x0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("y0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("z0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("resoX", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pztpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pptpc", 1);   // total mometum
  TTreeCont[kTPCGeant]->SetBranchStatus("masstpc", 1);   // mass TPC
  TTreeCont[kTPCGeant]->SetBranchStatus("betatpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("edeptpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("dedxtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("slengthtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("laytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("rowtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("parentID", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("iPadtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xtpc_pad", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ytpc_pad", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ztpc_pad", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("dxtpc_pad", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("dytpc_pad", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("dztpc_pad", 1);


  TTreeCont[kTPCGeant]->SetBranchAddress("evnum", &src.evnum);
  TTreeCont[kTPCGeant]->SetBranchAddress("nttpc", &src.nttpc);

  TTreeCont[kTPCGeant]->SetBranchAddress("ititpc", src.ititpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("idtpc", src.idtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("xtpc", src.xtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ytpc", src.ytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ztpc", src.ztpc); 
  TTreeCont[kTPCGeant]->SetBranchAddress("x0tpc", src.x0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("y0tpc", src.y0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("z0tpc", src.z0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("resoX", src.resoX);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxtpc", src.pxtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pytpc", src.pytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pztpc", src.pztpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pptpc", src.pptpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("masstpc", src.masstpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("betatpc", src.betatpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("edeptpc", src.edeptpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("dedxtpc", src.dedxtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("slengthtpc", src.slengthtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("laytpc", src.laytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("rowtpc", src.rowtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("parentID", src.parentID);
  TTreeCont[kTPCGeant]->SetBranchAddress("iPadtpc", src.iPadtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("xtpc_pad", src.xtpc_pad);
  TTreeCont[kTPCGeant]->SetBranchAddress("ytpc_pad", src.ytpc_pad);
  TTreeCont[kTPCGeant]->SetBranchAddress("ztpc_pad", src.ztpc_pad); 
  TTreeCont[kTPCGeant]->SetBranchAddress("dxtpc_pad", src.dxtpc_pad);
  TTreeCont[kTPCGeant]->SetBranchAddress("dytpc_pad", src.dytpc_pad);
  TTreeCont[kTPCGeant]->SetBranchAddress("dztpc_pad", src.dztpc_pad); 


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
