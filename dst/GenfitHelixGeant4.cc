// -*- C++ -*-

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

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
#include "TPCLocalTrackHelix.hh"
#include "TPCHit.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"
//#include "HypTPCTrack.hh"
//#include "HypTPCFitProcess.hh"

#include "DstHelper.hh"
#include <TRandom.h>
#include "DebugCounter.hh"

namespace
{
using namespace root;
using namespace dst;
const std::string& class_name("GenfitHelixGeant4");
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
ConfMan&            gConf = ConfMan::GetInstance();
const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
const UserParamMan& gUser = UserParamMan::GetInstance();
const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

const Int_t MaxTPCHits = 10000;
const Int_t MaxTPCTracks = 100;
const Int_t MaxTPCnHits = 50;

//const bool IsWithRes = false;
const bool IsWithRes = true;

//For GenFit Setting
const bool Const_field = true;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");
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
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{

  Int_t evnum;
  Int_t status;

  //GenFit outputs
  Int_t GFstatus;

  //Geant4 Generated Hits & Tracks
  Int_t g4nhTpc;               // Number of Hits
  Int_t g4ntTpc;               // Number of Tracks
  Int_t g4tidTpc[MaxTPCTracks];  // TrackId
  Int_t g4pidTpc[MaxTPCTracks];  // TrackPid
  Int_t g4nhtrack[MaxTPCHits]; // Number of Hits (in 1 tracks)

  //TrackSearchTPCHelix()
  Int_t ntTpc;                 // Number of Tracks
  Int_t nhtrack[MaxTPCTracks]; // Number of Hits (in 1 tracks)
  Double_t chisqr[MaxTPCTracks];
  Double_t helix_cx[MaxTPCTracks];
  Double_t helix_cy[MaxTPCTracks];
  Double_t helix_z0[MaxTPCTracks];
  Double_t helix_r[MaxTPCTracks];
  Double_t helix_dz[MaxTPCTracks];
  Double_t mom0_x[MaxTPCTracks];//Helix momentum at Y = 0
  Double_t mom0_y[MaxTPCTracks];
  Double_t mom0_z[MaxTPCTracks];
  Int_t hitlayer[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t hitpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_x[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_y[MaxTPCTracks][MaxTPCnHits];
  Double_t calpos_z[MaxTPCTracks][MaxTPCnHits];
  Double_t mom_x[MaxTPCTracks][MaxTPCnHits];
  Double_t mom_y[MaxTPCTracks][MaxTPCnHits];
  Double_t mom_z[MaxTPCTracks][MaxTPCnHits];

  Double_t g4mom_x[MaxTPCTracks][MaxTPCnHits];
  Double_t g4mom_y[MaxTPCTracks][MaxTPCnHits];
  Double_t g4mom_z[MaxTPCTracks][MaxTPCnHits];
  Int_t g4tid[MaxTPCTracks][MaxTPCnHits];
  Int_t g4pid[MaxTPCTracks][MaxTPCnHits];

  Double_t residual[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_x[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_y[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_z[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_px[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_py[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_pz[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_p[MaxTPCTracks][MaxTPCnHits];
};

//_____________________________________________________________________
struct Src
{
  Int_t evnum;
  Int_t nhittpc;                 // Number of Hits

  Int_t nhPrm;
  Double_t pidPrm[MaxTPCTracks];
  Double_t xPrm[MaxTPCTracks];
  Double_t yPrm[MaxTPCTracks];
  Double_t zPrm[MaxTPCTracks];
  Double_t pxPrm[MaxTPCTracks];
  Double_t pyPrm[MaxTPCTracks];
  Double_t pzPrm[MaxTPCTracks];

  Int_t ititpc[MaxTPCTracks];
  Int_t idtpc[MaxTPCTracks];
  Double_t xtpc[MaxTPCHits];//with resolution
  Double_t ytpc[MaxTPCHits];//with resolution
  Double_t ztpc[MaxTPCHits];//with resolution
  Double_t x0tpc[MaxTPCHits];//w/o resolution
  Double_t y0tpc[MaxTPCHits];//w/o resolution
  Double_t z0tpc[MaxTPCHits];//w/o resolution
  //  Double_t resoX[MaxTPCHits];
  Double_t pxtpc[MaxTPCHits];
  Double_t pytpc[MaxTPCHits];
  Double_t pztpc[MaxTPCHits];
  // Double_t pptpc[MaxTPCHits];   // total mometum
  // Double_t masstpc[MaxTPCHits];   // mass TPC
  // Double_t betatpc[MaxTPCHits];
  Double_t edeptpc[MaxTPCHits];
  // Double_t dedxtpc[MaxTPCHits];
  // Double_t slengthtpc[MaxTPCHits];
  // Int_t laytpc[MaxTPCHits];
  // Int_t rowtpc[MaxTPCHits];
  // Int_t parentID[MaxTPCHits];
  // Int_t iPadtpc[MaxTPCHits];//Pad number (0 origin)

  // Double_t xtpc_pad[MaxTPCHits];//pad center
  // Double_t ytpc_pad[MaxTPCHits];//pad center(dummy)
  // Double_t ztpc_pad[MaxTPCHits];//pad center

  // Double_t dxtpc_pad[MaxTPCHits];//x0tpc - xtpc
  // Double_t dytpc_pad[MaxTPCHits];//y0tpc - ytpc = 0 (dummy)
  // Double_t dztpc_pad[MaxTPCHits];//z0tpc - ztpc


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
  if( !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;
  // int nevent = GetEntries( TTreeCont );

  // CatchSignal::Set();

  // int ievent = 0;

  // // for(int ii=0; ii<100; ++ii){
  // //   std::cout<<"ii="<<ii<<std::endl;
  // //   ievent = 0;
  // for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
  //   gCounter.check();
  //   InitializeEvent();
  //   if( DstRead( ievent ) ) tree->Fill();
  // }
  // //  }


  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
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
  event.status = 0;
  event.evnum = 0;
  event.g4nhTpc = 0;
  event.g4ntTpc = 0;
  event.ntTpc = 0;

  event.GFstatus = 0;
  for(int i=0; i<MaxTPCTracks; ++i){
    event.g4tidTpc[i] =0;
    event.g4nhtrack[i] =0;
    event.nhtrack[i] =0;
    event.chisqr[i] =-9999.;
    event.helix_cx[i] =-9999.;
    event.helix_cy[i] =-9999.;
    event.helix_z0[i] =-9999.;
    event.helix_r[i] =-9999.;
    event.helix_dz[i] =-9999.;
    event.mom0_x[i] =-9999.;
    event.mom0_y[i] =-9999.;
    event.mom0_z[i] =-9999.;
    for(int j=0; j<MaxTPCnHits; ++j){
      event.hitlayer[i][j] =-999;
      event.mom_x[i][j] =-9999.;
      event.mom_y[i][j] =-9999.;
      event.mom_z[i][j] =-9999.;
      event.g4mom_x[i][j] =-9999.;
      event.g4mom_y[i][j] =-9999.;
      event.g4mom_z[i][j] =-9999.;
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
      event.residual_px[i][j] =-9999.;
      event.residual_py[i][j] =-9999.;
      event.residual_pz[i][j] =-9999.;
      event.residual_p[i][j] =-9999.;
      event.g4tid[i][j] =0;
      event.g4pid[i][j] =0;
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


  if( ievent%10000==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  event.evnum = src.evnum;
  event.g4nhTpc = src.nhittpc;
  for(int ihit=0; ihit<event.g4nhTpc; ++ihit){
    event.g4tidTpc[ihit] = src.ititpc[ihit];
    event.g4pidTpc[ihit] = src.idtpc[ihit];
    if(event.g4ntTpc<src.ititpc[ihit])
      event.g4ntTpc = src.ititpc[ihit];
    ++event.g4nhtrack[src.ititpc[ihit]-1];
  }

  // double u = src.pxPrm[0]/src.pzPrm[0];
  // double v = src.pyPrm[0]/src.pzPrm[0];
  // double cost = 1./std::sqrt(1.+u*u+v*v);
  // double theta=std::acos(cost)*math::Rad2Deg();
  // if(theta>20.)
  //   return true;

  HF1( 1, event.status++ );

  if(src.nhittpc<5)
    return true;

  HF1( 1, event.status++ );

  //Initiallize Geometry, Field, Fitter
  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  HF1( 2, event.GFstatus++ );

  //Initiallize the genfit track container
  HypTPCTask* GFtracks = HypTPCTask::GetInstance();
  GFtracks -> Init();
  //-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)
  GFtracks -> SetVerbosity(0);
  std::cout<<"verbosity=GFtracks -> GetVerbosity()"<<std::endl;
  HF1( 2, event.GFstatus++ );

  DCAnalyzer *DCAna = new DCAnalyzer();
  if(IsWithRes) DCAna->DecodeTPCHitsGeant4(src.nhittpc, src.xtpc, src.ytpc, src.ztpc, src.edeptpc);
  else DCAna->DecodeTPCHitsGeant4(src.nhittpc, src.x0tpc, src.y0tpc, src.z0tpc, src.edeptpc);

  DCAna->TrackSearchTPCHelix();

  HF1( 1, event.status++ );
  HF1( 2, event.GFstatus++ );

  int ntTpc = DCAna->GetNTracksTPCHelix();
  if( MaxHits<ntTpc ){
    std::cout << "#W " << func_name << " "
      	      << "too many ntTpc " << ntTpc << "/" << MaxHits << std::endl;
    ntTpc = MaxHits;
  }
  HF1( 10, ntTpc );

  event.ntTpc = ntTpc;
  for( int it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp= DCAna->GetTrackTPCHelix(it);
    if(!tp) continue;
    int nh=tp->GetNHit();
    double chisqr=tp->GetChiSquare();
    double cx=tp->Getcx(), cy=tp->Getcy();
    double z0=tp->Getz0(), r=tp->Getr();
    double dz=tp->Getdz();
    TVector3 mom0 = tp->GetMom0();

    event.nhtrack[it] = nh;
    event.chisqr[it] = chisqr;
    event.helix_cx[it] = cx;
    event.helix_cy[it] = cy;
    event.helix_z0[it] = z0;
    event.helix_r[it]  = r;
    event.helix_dz[it] = dz;
    event.mom0_x[it] = mom0.X();
    event.mom0_y[it] = mom0.Y();
    event.mom0_z[it] = mom0.Z();

    for( int ih=0; ih<nh; ++ih ){
      TPCLTrackHit *hit=tp->GetHit(ih);
      if(!hit) continue;
      int layerId = 0;
      layerId = hit->GetLayer();
      TVector3 hitpos = hit->GetLocalHitPos();
      TVector3 calpos = hit->GetLocalCalPosHelix();
      double residual = hit->GetResidual();
      TVector3 res_vect = hit->GetResidualVect();
      TVector3 mom = hit->GetMomentumHelix();

      for( int ih2=0; ih2<src.nhittpc; ++ih2 ){
	TVector3 setpos(src.x0tpc[ih2], src.y0tpc[ih2], src.z0tpc[ih2]);
	TVector3 d = setpos - hitpos;
	if(fabs(d.Mag()<0.1)){
	  event.g4mom_x[it][ih] = src.pxtpc[ih2]*1000.;
	  event.g4mom_y[it][ih] = src.pytpc[ih2]*1000.;
	  event.g4mom_z[it][ih] = src.pztpc[ih2]*1000.;
	  double g4mom_mag = sqrt(src.pxtpc[ih2]*src.pxtpc[ih2]
				 +src.pytpc[ih2]*src.pytpc[ih2]
				 +src.pztpc[ih2]*src.pztpc[ih2])*1000.;
	  event.residual_px[it][ih] = mom.x() - src.pxtpc[ih2]*1000.;//MeV/c
	  event.residual_py[it][ih] = mom.y() - src.pytpc[ih2]*1000.;//MeV/c
	  event.residual_pz[it][ih] = mom.z() - src.pztpc[ih2]*1000.;//MeV/c
	  event.residual_p[it][ih] = mom.Mag() - g4mom_mag;//MeV/c
	  event.g4tid[it][ih] = src.ititpc[ih2];
	  event.g4pid[it][ih] = src.idtpc[ih2];
	  break;
	}
      } //ih2
      event.hitlayer[it][ih] = layerId;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.mom_x[it][ih] = mom.x();
      event.mom_y[it][ih] = mom.y();
      event.mom_z[it][ih] = mom.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = res_vect.x();
      event.residual_y[it][ih] = res_vect.y();
      event.residual_z[it][ih] = res_vect.z();
    } //ih
    //Add tracks into the GenFit TrackCand
    GFtracks -> AddHelixTrack(event.g4tid[it][0], tp);
  } //it
  HF1( 2, event.GFstatus++ );
  GFtracks -> FitTracks();
  HF1( 2, event.GFstatus++ );

#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[g4nhTpc]: "<<std::setw(2)<<src.nhittpc<<" "<<std::endl;
#endif

  HF1( 1, event.status++ );

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
  HB1( 2, "GenFit Status", 21, 0., 21. );
  HB1( 10, "NTrack TPC", 20, 0., 20. );

  HBTree( "tpc", "tree of tpc" );

  tree->Branch("status", &event.status, "status/I" );
  tree->Branch("evnum", &event.evnum, "evnum/I" );

  //Geant4 Generated Hits & Tracks
  tree->Branch("g4nhTpc",&event.g4nhTpc,"g4nhTpc/I");
  tree->Branch("g4ntTpc",&event.g4ntTpc,"g4ntTpc/I");
  tree->Branch("g4tidTpc",event.g4tidTpc,"g4tidTpc[g4nhTpc]/I");
  tree->Branch("g4pidTpc",event.g4pidTpc,"g4pidTpc[g4nhTpc]/I");
  tree->Branch("g4nhtrack",event.g4nhtrack,"g4nhtrack[g4ntTpc]/I");

  //TrackSearchTPCHelix()
  tree->Branch("ntTpc",&event.ntTpc,"ntTpc/I");
  tree->Branch("nhtrack",event.nhtrack,"nhtrack[ntTpc]/I");
  tree->Branch("chisqr",event.chisqr,"chisqr[ntTpc]/D");
  tree->Branch("helix_cx",event.helix_cx,"helix_cx[ntTpc]/D");
  tree->Branch("helix_cy",event.helix_cy,"helix_cy[ntTpc]/D");
  tree->Branch("helix_z0",event.helix_z0,"helix_z0[ntTpc]/D");
  tree->Branch("helix_r" ,event.helix_r, "helix_r[ntTpc]/D");
  tree->Branch("helix_dz",event.helix_dz,"helix_dz[ntTpc]/D");
  tree->Branch("mom0_x",event.mom0_x,"mom0_x[ntTpc]/D"); // Momentum at Y = 0
  tree->Branch("mom0_y",event.mom0_y,"mom0_y[ntTpc]/D"); // Momentum at Y = 0
  tree->Branch("mom0_z",event.mom0_z,"mom0_z[ntTpc]/D"); // Momentum at Y = 0
  tree->Branch("hitlayer",event.hitlayer,"hitlayer[ntTpc][MaxTPCnHits]/I");
  tree->Branch("hitpos_x",event.hitpos_x,"hitpos_x[ntTpc][MaxTPCnHits]/D");
  tree->Branch("hitpos_y",event.hitpos_y,"hitpos_y[ntTpc][MaxTPCnHits]/D");
  tree->Branch("hitpos_z",event.hitpos_z,"hitpos_z[ntTpc][MaxTPCnHits]/D");
  tree->Branch("calpos_x",event.calpos_x,"calpos_x[ntTpc][MaxTPCnHits]/D");
  tree->Branch("calpos_y",event.calpos_y,"calpos_y[ntTpc][MaxTPCnHits]/D");
  tree->Branch("calpos_z",event.calpos_z,"calpos_z[ntTpc][MaxTPCnHits]/D");
  tree->Branch("mom_x",event.mom_x,"mom_x[ntTpc][MaxTPCnHits]/D");
  tree->Branch("mom_y",event.mom_y,"mom_y[ntTpc][MaxTPCnHits]/D");
  tree->Branch("mom_z",event.mom_z,"mom_z[ntTpc][MaxTPCnHits]/D");

  //Geant4 Generated input
  tree->Branch("g4mom_x",event.g4mom_x,"g4mom_x[ntTpc][MaxTPCnHits]/D");
  tree->Branch("g4mom_y",event.g4mom_y,"g4mom_y[ntTpc][MaxTPCnHits]/D");
  tree->Branch("g4mom_z",event.g4mom_z,"g4mom_z[ntTpc][MaxTPCnHits]/D");
  tree->Branch("g4tid",event.g4tid,"g4tid[ntTpc][MaxTPCnHits]/I");
  tree->Branch("g4pid",event.g4pid,"g4pid[ntTpc][MaxTPCnHits]/I");

  tree->Branch("residual",event.residual,"residual[ntTpc][MaxTPCnHits]/D");
  tree->Branch("residual_x",event.residual_x,"residual_x[ntTpc][MaxTPCnHits]/D");
  tree->Branch("residual_y",event.residual_y,"residual_y[ntTpc][MaxTPCnHits]/D");
  tree->Branch("residual_z",event.residual_z,"residual_z[ntTpc][MaxTPCnHits]/D");
  tree->Branch("residual_px",event.residual_px,"residual_px[ntTpc][MaxTPCnHits]/D");
  tree->Branch("residual_py",event.residual_py,"residual_py[ntTpc][MaxTPCnHits]/D");
  tree->Branch("residual_pz",event.residual_pz,"residual_pz[ntTpc][MaxTPCnHits]/D");
  tree->Branch("residual_p",event.residual_p,"residual_p[ntTpc][MaxTPCnHits]/D");

  /*
  tree->Branch("nPrm",&src.nPrm,"nPrm/I");
  tree->Branch("pidPrm",src.pidPrm,"pidPrm[nPrm]/I");
  tree->Branch("xPrm",src.xPrm,"xPrm[nPrm]/D");
  tree->Branch("yPrm",src.yPrm,"yPrm[nPrm]/D");
  tree->Branch("zPrm",src.zPrm,"zPrm[nPrm]/D");
  tree->Branch("pxPrm",src.pxPrm,"pxPrm[nPrm]/D");
  tree->Branch("pyPrm",src.pyPrm,"pyPrm[nPrm]/D");
  tree->Branch("pzPrm",src.pzPrm,"pzPrm[nPrm]/D");
  */
  ////////// Bring Address From Dst
  TTreeCont[kTPCGeant]->SetBranchStatus("*", 0);
  TTreeCont[kTPCGeant]->SetBranchStatus("evnum",  1);
  /*
  TTreeCont[kTPCGeant]->SetBranchStatus("nhPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pyPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzPrm",  1);
  */
  TTreeCont[kTPCGeant]->SetBranchStatus("nhittpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ititpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("idtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ztpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("x0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("y0tpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("z0tpc", 1);
  //  TTreeCont[kTPCGeant]->SetBranchStatus("resoX", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pztpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("pptpc", 1);   // total mometum
  // TTreeCont[kTPCGeant]->SetBranchStatus("masstpc", 1);   // mass TPC
  // TTreeCont[kTPCGeant]->SetBranchStatus("betatpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("edeptpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dedxtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("slengthtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("laytpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("rowtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("parentID", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("iPadtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("xtpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("ytpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("ztpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dxtpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dytpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dztpc_pad", 1);

  TTreeCont[kTPCGeant]->SetBranchAddress("evnum", &src.evnum);
  /*
  TTreeCont[kTPCGeant]->SetBranchAddress("nhPrm", &src.nhPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidPrm", src.pidPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("xPrm", src.xPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("yPrm", src.yPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("zPrm", src.zPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxPrm", src.pxPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyPrm", src.pyPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzPrm", src.pzPrm);
  */
  TTreeCont[kTPCGeant]->SetBranchAddress("nhittpc", &src.nhittpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ititpc", src.ititpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("idtpc", src.idtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("xtpc", src.xtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ytpc", src.ytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ztpc", src.ztpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("x0tpc", src.x0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("y0tpc", src.y0tpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("z0tpc", src.z0tpc);
  //  TTreeCont[kTPCGeant]->SetBranchAddress("resoX", src.resoX);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxtpc", src.pxtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pytpc", src.pytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pztpc", src.pztpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("pptpc", src.pptpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("masstpc", src.masstpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("betatpc", src.betatpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("edeptpc", src.edeptpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dedxtpc", src.dedxtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("slengthtpc", src.slengthtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("laytpc", src.laytpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("rowtpc", src.rowtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("parentID", src.parentID);
  // TTreeCont[kTPCGeant]->SetBranchAddress("iPadtpc", src.iPadtpc);
  // TTreeCont[kTPCGeant]->SetBranchAddress("xtpc_pad", src.xtpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("ytpc_pad", src.ytpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("ztpc_pad", src.ztpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dxtpc_pad", src.dxtpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dytpc_pad", src.dytpc_pad);
  // TTreeCont[kTPCGeant]->SetBranchAddress("dztpc_pad", src.dztpc_pad);

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
