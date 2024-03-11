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
#include "TPCParamMan.hh"
#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCHit.hh"

#include "DstHelper.hh"
#include <TRandom.h>
#include "DebugCounter.hh"

namespace
{
using namespace root;
using namespace dst;
using namespace std;
const std::string& class_name("DstTPCTrackingHelixGeant4");
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
ConfMan&            gConf = ConfMan::GetInstance();
const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
const UserParamMan& gUser = UserParamMan::GetInstance();
const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

const Int_t MaxTPCHits = 10000;
const Int_t MaxTPCTracks = 100;
const Int_t MaxTPCnHits = 50;
const Int_t MaxG4Hits = 500;
const double truncatedMean = 0.8; //80%
  //const bool IsWithRes = false;
const bool IsWithRes = true;
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

  Int_t status;
  Int_t evnum;
  Int_t nhittpc;                 // Number of Hits
  Int_t ntTpc;                   // Number of Tracks
  Int_t max_ititpc;
  Int_t ititpc[MaxTPCHits];
  Int_t nhittpc_iti[MaxTPCHits];
  vector<Int_t> nhtrack;
  vector<Double_t> chisqr;
  vector<Double_t> helix_cx;
  vector<Double_t> helix_cy;
  vector<Double_t> helix_z0;
  vector<Double_t> helix_r;
  vector<Double_t> helix_dz;
  vector<Double_t> dE;
  vector<Double_t> dEdx; //reference dedx
  vector<Double_t> dz_factor;
  vector<Int_t> charge;
  vector<Double_t> path;
  vector<Int_t> pid;
  vector<Double_t> mom0;//Helix momentum at Y = 0
  vector<Double_t> mom0_x;//Helix momentum at Y = 0
  vector<Double_t> mom0_y;
  vector<Double_t> mom0_z;
  vector<vector<Double_t>> hitlayer;
  vector<vector<Double_t>> hitpos_x;
  vector<vector<Double_t>> hitpos_y;
  vector<vector<Double_t>> hitpos_z;
  vector<vector<Double_t>> helix_t;
  vector<vector<Double_t>> calpos_x;
  vector<vector<Double_t>> calpos_y;
  vector<vector<Double_t>> calpos_z;
  vector<vector<Double_t>> mom_x;
  vector<vector<Double_t>> mom_y;
  vector<vector<Double_t>> mom_z;

  vector<vector<Double_t>> residual;
  vector<vector<Double_t>> residual_x;
  vector<vector<Double_t>> residual_y;
  vector<vector<Double_t>> residual_z;
  vector<vector<Double_t>> resolution;
  vector<vector<Double_t>> resolution_x;
  vector<vector<Double_t>> resolution_y;
  vector<vector<Double_t>> resolution_z;
  vector<vector<Double_t>> pathhit;
  vector<vector<Double_t>> alpha;
  vector<vector<Double_t>> track_cluster_de;
  vector<vector<Double_t>> track_cluster_mrow;

	//Geant4
  Int_t iti_g[MaxTPCTracks][MaxTPCnHits];
  Int_t idtpc[MaxTPCHits];
  Int_t ID[MaxTPCHits];
  Int_t PID[MaxTPCHits];
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
	Double_t momg_x[MaxTPCTracks][MaxTPCnHits];
	Double_t momg_y[MaxTPCTracks][MaxTPCnHits];
	Double_t momg_z[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_p[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_px[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_py[MaxTPCTracks][MaxTPCnHits];
  Double_t residual_pz[MaxTPCTracks][MaxTPCnHits];
 	
	Int_t NumberOfTracks;
	Int_t PIDOfTrack[1000];
	Int_t ParentIDOfTrack[1000];
	Double_t VertexOfTrack_x[1000];
	Double_t VertexOfTrack_y[1000];
	Double_t VertexOfTrack_z[1000];
	Double_t MomentumOfTrack[1000];
	Double_t MomentumOfTrack_x[1000];
	Double_t MomentumOfTrack_y[1000];
	Double_t MomentumOfTrack_z[1000];
 
	Double_t MomXi_x,MomXi_y,MomXi_z;//At decay vtx
	Double_t SpinXi_x,SpinXi_y,SpinXi_z;
	Double_t ThXi_CM;
	Double_t MomLd_x,MomLd_y,MomLd_z;
	Double_t SpinLd_x,SpinLd_y,SpinLd_z;
	Double_t ThLd_CM;

	//KpXi Kinematics, Momentum at production vtx
	double PKm,PKm_x,PKm_y,PKm_z;
	double PKp,PKp_x,PKp_y,PKp_z;
	double PXi,PXi_x,PXi_y,PXi_z;
	double PPi2,PPi2_x,PPi2_y,PPi2_z;//Pi from Xi
	
	double PLd,PLd_x,PLd_y,PLd_z;
	double PPi1,PPi1_x,PPi1_y,PPi1_z;//Pi from Ld
	double PP,PP_x,PP_y,PP_z;
	
	//HToF
  Int_t nhHtof;
  Int_t tidHtof[MaxG4Hits];
  Int_t pidHtof[MaxG4Hits];
  Int_t didHtof[MaxG4Hits];
  Int_t prtHtof[MaxG4Hits];
  Int_t qHtof[MaxG4Hits];
  Double_t xHtof[MaxG4Hits];
  Double_t yHtof[MaxG4Hits];
  Double_t zHtof[MaxG4Hits];
  Double_t pxHtof[MaxG4Hits];
  Double_t pyHtof[MaxG4Hits];
  Double_t pzHtof[MaxG4Hits];
  Double_t ppHtof[MaxG4Hits];
  Double_t deHtof[MaxG4Hits];
  Double_t tHtof[MaxG4Hits];
	//



	Int_t nhSch;
  Int_t tidSch[MaxG4Hits];
  Int_t pidSch[MaxG4Hits];
  Int_t didSch[MaxG4Hits];
  Int_t prtSch[MaxG4Hits];
  Int_t qSch[MaxG4Hits];
  Double_t xSch[MaxG4Hits];
  Double_t ySch[MaxG4Hits];
  Double_t zSch[MaxG4Hits];
  Double_t pxSch[MaxG4Hits];
  Double_t pySch[MaxG4Hits];
  Double_t pzSch[MaxG4Hits];
  Double_t ppSch[MaxG4Hits];
  Double_t deSch[MaxG4Hits];
  Double_t tSch[MaxG4Hits];
	
	//Lac	
  Int_t nhLac;
  Int_t tidLac[MaxG4Hits];
  Int_t pidLac[MaxG4Hits];
  Int_t didLac[MaxG4Hits];
  Int_t prtLac[MaxG4Hits];
  Int_t qLac[MaxG4Hits];
  Double_t xLac[MaxG4Hits];
  Double_t yLac[MaxG4Hits];
  Double_t zLac[MaxG4Hits];
  Double_t pxLac[MaxG4Hits];
  Double_t pyLac[MaxG4Hits];
  Double_t pzLac[MaxG4Hits];
  Double_t ppLac[MaxG4Hits];
  Double_t deLac[MaxG4Hits];
  Double_t tLac[MaxG4Hits];
	
	//FToF
  Int_t nhFtof;
  Int_t tidFtof[MaxG4Hits];
  Int_t pidFtof[MaxG4Hits];
  Int_t didFtof[MaxG4Hits];
  Int_t prtFtof[MaxG4Hits];
  Int_t qFtof[MaxG4Hits];
  Double_t xFtof[MaxG4Hits];
  Double_t yFtof[MaxG4Hits];
  Double_t zFtof[MaxG4Hits];
  Double_t pxFtof[MaxG4Hits];
  Double_t pyFtof[MaxG4Hits];
  Double_t pzFtof[MaxG4Hits];
  Double_t ppFtof[MaxG4Hits];
  Double_t deFtof[MaxG4Hits];
  Double_t tFtof[MaxG4Hits];
 
	//SDC
	Int_t nhSdc;
  Int_t tidSdc[MaxG4Hits];
  Int_t pidSdc[MaxG4Hits];
  Int_t didSdc[MaxG4Hits];
  Int_t prtSdc[MaxG4Hits];
  Int_t qSdc[MaxG4Hits];
  Double_t xSdc[MaxG4Hits];
  Double_t ySdc[MaxG4Hits];
  Double_t zSdc[MaxG4Hits];
  Double_t pxSdc[MaxG4Hits];
  Double_t pySdc[MaxG4Hits];
  Double_t pzSdc[MaxG4Hits];
  Double_t ppSdc[MaxG4Hits];
  Double_t deSdc[MaxG4Hits];
  Double_t tSdc[MaxG4Hits];

	//WC
  Int_t nhWc;
  Int_t tidWc[MaxG4Hits];
  Int_t pidWc[MaxG4Hits];
  Int_t didWc[MaxG4Hits];
  Int_t prtWc[MaxG4Hits];
  Int_t qWc[MaxG4Hits];
  Double_t xWc[MaxG4Hits];
  Double_t yWc[MaxG4Hits];
  Double_t zWc[MaxG4Hits];
  Double_t pxWc[MaxG4Hits];
  Double_t pyWc[MaxG4Hits];
  Double_t pzWc[MaxG4Hits];
  Double_t ppWc[MaxG4Hits];
  Double_t deWc[MaxG4Hits];
  Double_t tWc[MaxG4Hits];

	//BVH
  Int_t nhBvh;
  Int_t tidBvh[MaxG4Hits];
  Int_t pidBvh[MaxG4Hits];
  Int_t didBvh[MaxG4Hits];
  Int_t prtBvh[MaxG4Hits];
  Int_t qBvh[MaxG4Hits];
  Double_t xBvh[MaxG4Hits];
  Double_t yBvh[MaxG4Hits];
  Double_t zBvh[MaxG4Hits];
  Double_t pxBvh[MaxG4Hits];
  Double_t pyBvh[MaxG4Hits];
  Double_t pzBvh[MaxG4Hits];
  Double_t ppBvh[MaxG4Hits];
  Double_t deBvh[MaxG4Hits];
  Double_t tBvh[MaxG4Hits];
	void Clear(){
		PKm = qnan;PKm_x=qnan;PKm_y=qnan;PKm_z=qnan;
		PKp = qnan;PKp_x=qnan;PKp_y=qnan;PKp_z=qnan;
		PXi = qnan;PXi_x=qnan;PXi_y=qnan;PXi_z=qnan;
		PPi2 = qnan;PPi2_x=qnan;PPi2_y=qnan;PPi2_z=qnan;
		PLd = qnan;PLd_x=qnan;PLd_y=qnan;PLd_z=qnan;
		PPi1 = qnan;PPi1_x=qnan;PPi1_y=qnan;PPi1_z=qnan;
		PP = qnan;PP_x=qnan;PP_y=qnan;PP_z=qnan;
		nhtrack.clear();
		chisqr.clear();
		helix_cx.clear();
		helix_cy.clear();
		helix_z0.clear();
		helix_r.clear();
		helix_dz.clear();
		dz_factor.clear();
    dE.clear();
    dEdx.clear();
		mom0_x.clear();
		mom0_y.clear();
		mom0_z.clear();
		mom0.clear();
		charge.clear();
		path.clear();
		pid.clear();

		hitlayer.clear();
		hitpos_x.clear();
		hitpos_y.clear();
		hitpos_z.clear();
		calpos_x.clear();
		calpos_y.clear();
		calpos_z.clear();
		resolution.clear();
		resolution_x.clear();
		resolution_y.clear();
		resolution_z.clear();
		residual.clear();
		residual_x.clear();
		residual_y.clear();
		residual_z.clear();
		helix_t.clear();
    
		track_cluster_de.clear();
    track_cluster_mrow.clear();
	};
};

//_____________________________________________________________________
struct Src
{
  Int_t evnum;
  Int_t nhittpc;                 // Number of Hits

  Int_t nhPrm;
  Double_t xPrm[MaxTPCTracks];
  Double_t yPrm[MaxTPCTracks];
  Double_t zPrm[MaxTPCTracks];
  Double_t pxPrm[MaxTPCTracks];
  Double_t pyPrm[MaxTPCTracks];
  Double_t pzPrm[MaxTPCTracks];
  Double_t ppPrm[MaxTPCTracks];

 	Int_t NumberOfTracks;
	Int_t PIDOfTrack[1000];
	Int_t ParentIDOfTrack[1000];
	Double_t VertexOfTrack_x[1000];
	Double_t VertexOfTrack_y[1000];
	Double_t VertexOfTrack_z[1000];
	Double_t MomentumOfTrack[1000];
	Double_t MomentumOfTrack_x[1000];
	Double_t MomentumOfTrack_y[1000];
	Double_t MomentumOfTrack_z[1000];
	Double_t MomXi_x,MomXi_y,MomXi_z;
	Double_t SpinXi_x,SpinXi_y,SpinXi_z;
	Double_t ThXi_CM;
	Double_t MomLd_x,MomLd_y,MomLd_z;
	Double_t SpinLd_x,SpinLd_y,SpinLd_z;
	Double_t ThLd_CM;


  Int_t ititpc[MaxTPCHits];
  Int_t idtpc[MaxTPCHits];
  Int_t parentID[MaxTPCHits];
	
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
  Double_t pptpc[MaxTPCHits];   // total mometum
  
	// Double_t masstpc[MaxTPCHits];   // mass TPC
  // Double_t betatpc[MaxTPCHits];
  Double_t edeptpc[MaxTPCHits];
  // Int_t laytpc[MaxTPCHits];
  // Int_t rowtpc[MaxTPCHits];
  // Int_t iPadtpc[MaxTPCHits];//Pad number (0 origin)
	
	//HToF
  Int_t nhHtof;
  Int_t tidHtof[MaxG4Hits];
  Int_t pidHtof[MaxG4Hits];
  Int_t didHtof[MaxG4Hits];
  Int_t prtHtof[MaxG4Hits];
  Int_t qHtof[MaxG4Hits];
  Double_t xHtof[MaxG4Hits];
  Double_t yHtof[MaxG4Hits];
  Double_t zHtof[MaxG4Hits];
  Double_t pxHtof[MaxG4Hits];
  Double_t pyHtof[MaxG4Hits];
  Double_t pzHtof[MaxG4Hits];
  Double_t ppHtof[MaxG4Hits];
  Double_t deHtof[MaxG4Hits];
  Double_t tHtof[MaxG4Hits];

	//SCH
  Int_t nhSch;
  Int_t tidSch[MaxG4Hits];
  Int_t pidSch[MaxG4Hits];
  Int_t didSch[MaxG4Hits];
  Int_t prtSch[MaxG4Hits];
  Int_t qSch[MaxG4Hits];
  Double_t xSch[MaxG4Hits];
  Double_t ySch[MaxG4Hits];
  Double_t zSch[MaxG4Hits];
  Double_t pxSch[MaxG4Hits];
  Double_t pySch[MaxG4Hits];
  Double_t pzSch[MaxG4Hits];
  Double_t ppSch[MaxG4Hits];
  Double_t deSch[MaxG4Hits];
  Double_t tSch[MaxG4Hits];
	
	//Lac	
  Int_t nhLac;
  Int_t tidLac[MaxG4Hits];
  Int_t pidLac[MaxG4Hits];
  Int_t didLac[MaxG4Hits];
  Int_t prtLac[MaxG4Hits];
  Int_t qLac[MaxG4Hits];
  Double_t xLac[MaxG4Hits];
  Double_t yLac[MaxG4Hits];
  Double_t zLac[MaxG4Hits];
  Double_t pxLac[MaxG4Hits];
  Double_t pyLac[MaxG4Hits];
  Double_t pzLac[MaxG4Hits];
  Double_t ppLac[MaxG4Hits];
  Double_t deLac[MaxG4Hits];
  Double_t tLac[MaxG4Hits];
	
	//FToF
  Int_t nhFtof;
  Int_t tidFtof[MaxG4Hits];
  Int_t pidFtof[MaxG4Hits];
  Int_t didFtof[MaxG4Hits];
  Int_t prtFtof[MaxG4Hits];
  Int_t qFtof[MaxG4Hits];
  Double_t xFtof[MaxG4Hits];
  Double_t yFtof[MaxG4Hits];
  Double_t zFtof[MaxG4Hits];
  Double_t pxFtof[MaxG4Hits];
  Double_t pyFtof[MaxG4Hits];
  Double_t pzFtof[MaxG4Hits];
  Double_t ppFtof[MaxG4Hits];
  Double_t deFtof[MaxG4Hits];
  Double_t tFtof[MaxG4Hits];
  
	//SDC
	Int_t nhSdc;
  Int_t tidSdc[MaxG4Hits];
  Int_t pidSdc[MaxG4Hits];
  Int_t didSdc[MaxG4Hits];
  Int_t prtSdc[MaxG4Hits];
  Int_t qSdc[MaxG4Hits];
  Double_t xSdc[MaxG4Hits];
  Double_t ySdc[MaxG4Hits];
  Double_t zSdc[MaxG4Hits];
  Double_t pxSdc[MaxG4Hits];
  Double_t pySdc[MaxG4Hits];
  Double_t pzSdc[MaxG4Hits];
  Double_t ppSdc[MaxG4Hits];
  Double_t deSdc[MaxG4Hits];
  Double_t tSdc[MaxG4Hits];

	//WC
  Int_t nhWc;
  Int_t tidWc[MaxG4Hits];
  Int_t pidWc[MaxG4Hits];
  Int_t didWc[MaxG4Hits];
  Int_t prtWc[MaxG4Hits];
  Int_t qWc[MaxG4Hits];
  Double_t xWc[MaxG4Hits];
  Double_t yWc[MaxG4Hits];
  Double_t zWc[MaxG4Hits];
  Double_t pxWc[MaxG4Hits];
  Double_t pyWc[MaxG4Hits];
  Double_t pzWc[MaxG4Hits];
  Double_t ppWc[MaxG4Hits];
  Double_t deWc[MaxG4Hits];
  Double_t tWc[MaxG4Hits];
	
	//BVH
  Int_t nhBvh;
  Int_t tidBvh[MaxG4Hits];
  Int_t pidBvh[MaxG4Hits];
  Int_t didBvh[MaxG4Hits];
  Int_t prtBvh[MaxG4Hits];
  Int_t qBvh[MaxG4Hits];
  Double_t xBvh[MaxG4Hits];
  Double_t yBvh[MaxG4Hits];
  Double_t zBvh[MaxG4Hits];
  Double_t pxBvh[MaxG4Hits];
  Double_t pyBvh[MaxG4Hits];
  Double_t pzBvh[MaxG4Hits];
  Double_t ppBvh[MaxG4Hits];
  Double_t deBvh[MaxG4Hits];
  Double_t tBvh[MaxG4Hits];

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
	event.Clear();
  event.status   = 0;
  event.evnum = 0;
  event.nhittpc = 0;
  event.ntTpc = 0;
  event.max_ititpc = 0;

  for(int i=0; i<MaxTPCTracks; ++i){
    event.ititpc[i] =0;
    event.nhittpc_iti[i] =0;
    for(int j=0; j<MaxTPCnHits; ++j){
      event.momg_x[i][j] =-9999.;
      event.momg_y[i][j] =-9999.;
      event.momg_z[i][j] =-9999.;
      event.residual_px[i][j] =-9999.;
      event.residual_py[i][j] =-9999.;
      event.residual_pz[i][j] =-9999.;
      event.residual_p[i][j] =-9999.;
      event.iti_g[i][j] =0;
		}
  }
	event.NumberOfTracks = -1;
	for(int i=0;i<1000;++i){
		event.PIDOfTrack[i]=qnan;
		event.ParentIDOfTrack[i]=qnan;
		event.VertexOfTrack_x[i]=qnan;
		event.VertexOfTrack_y[i]=qnan;
		event.VertexOfTrack_z[i]=qnan;
		event.MomentumOfTrack[i]=qnan;
		event.MomentumOfTrack_x[i]=qnan;
		event.MomentumOfTrack_y[i]=qnan;
		event.MomentumOfTrack_z[i]=qnan;
	}

	event.MomXi_x = 0;
	event.MomXi_y = 0;
	event.MomXi_z = 0;
	event.SpinXi_x = 0;
	event.SpinXi_y = 0;
	event.SpinXi_z = 0;
	event.ThXi_CM = 0;
	
	event.MomLd_x = 0;
	event.MomLd_y = 0;
	event.MomLd_z = 0;
	event.SpinLd_x = 0;
	event.SpinLd_y = 0;
	event.SpinLd_z = 0;
	event.ThLd_CM = 0;

	event.nhHtof=-1;
	event.nhSch = -1;
	event.nhLac = -1;
	event.nhFtof=-1;
	event.nhSdc=-1;
	event.nhWc =-1;
	event.nhBvh=-1;
  for(int i=0;i<500;++i){  
		event.ititpc[i]=qnan;
    event.idtpc[i]=qnan;
    event.ID[i]=qnan;
    event.PID[i]=qnan;
    event.xtpc[i]=qnan;
    event.ytpc[i]=qnan;
    event.ztpc[i]=qnan;
    event.x0tpc[i]=qnan;
    event.y0tpc[i]=qnan;
    event.z0tpc[i]=qnan;
    event.pxtpc[i]=qnan;
    event.pytpc[i]=qnan;
    event.pztpc[i]=qnan;
		
		event.tidHtof[i]=qnan;
		event.pidHtof[i]=qnan;
		event.didHtof[i]=qnan;
		event.prtHtof[i]=qnan;
		event.qHtof[i]=qnan;
		event.xHtof[i]=qnan;
		event.yHtof[i]=qnan;
		event.zHtof[i]=qnan;
		event.pxHtof[i]=qnan;
		event.pyHtof[i]=qnan;
		event.pzHtof[i]=qnan;
		event.ppHtof[i]=qnan;
		event.deHtof[i]=qnan;
		event.tHtof[i]=qnan;
		
		event.tidSch[i]=qnan;
		event.pidSch[i]=qnan;
		event.didSch[i]=qnan;
		event.prtSch[i]=qnan;
		event.qSch[i]=qnan;
		event.xSch[i]=qnan;
		event.ySch[i]=qnan;
		event.zSch[i]=qnan;
		event.pxSch[i]=qnan;
		event.pySch[i]=qnan;
		event.pzSch[i]=qnan;
		event.ppSch[i]=qnan;
		event.deSch[i]=qnan;
		event.tSch[i]=qnan;
		
		event.tidLac[i]=qnan;
		event.pidLac[i]=qnan;
		event.didLac[i]=qnan;
		event.prtLac[i]=qnan;
		event.qLac[i]=qnan;
		event.xLac[i]=qnan;
		event.yLac[i]=qnan;
		event.zLac[i]=qnan;
		event.pxLac[i]=qnan;
		event.pyLac[i]=qnan;
		event.pzLac[i]=qnan;
		event.ppLac[i]=qnan;
		event.deLac[i]=qnan;
		event.tLac[i]=qnan;
		
		event.tidFtof[i]=qnan;
		event.pidFtof[i]=qnan;
		event.didFtof[i]=qnan;
		event.prtFtof[i]=qnan;
		event.qFtof[i]=qnan;
		event.xFtof[i]=qnan;
		event.yFtof[i]=qnan;
		event.zFtof[i]=qnan;
		event.pxFtof[i]=qnan;
		event.pyFtof[i]=qnan;
		event.pzFtof[i]=qnan;
		event.ppFtof[i]=qnan;
		event.deFtof[i]=qnan;
		event.tFtof[i]=qnan;
		
		event.tidSdc[i]=qnan;
		event.pidSdc[i]=qnan;
		event.didSdc[i]=qnan;
		event.prtSdc[i]=qnan;
		event.qSdc[i]=qnan;
		event.xSdc[i]=qnan;
		event.ySdc[i]=qnan;
		event.zSdc[i]=qnan;
		event.pxSdc[i]=qnan;
		event.pySdc[i]=qnan;
		event.pzSdc[i]=qnan;
		event.ppSdc[i]=qnan;
		event.deSdc[i]=qnan;
		event.tSdc[i]=qnan;
		
		event.tidWc[i]=qnan;
		event.pidWc[i]=qnan;
		event.didWc[i]=qnan;
		event.prtWc[i]=qnan;
		event.qWc[i]=qnan;
		event.xWc[i]=qnan;
		event.yWc[i]=qnan;
		event.zWc[i]=qnan;
		event.pxWc[i]=qnan;
		event.pyWc[i]=qnan;
		event.pzWc[i]=qnan;
		event.ppWc[i]=qnan;
		event.deWc[i]=qnan;
		event.tWc[i]=qnan;
		
		event.tidBvh[i]=qnan;
		event.pidBvh[i]=qnan;
		event.didBvh[i]=qnan;
		event.prtBvh[i]=qnan;
		event.qBvh[i]=qnan;
		event.xBvh[i]=qnan;
		event.yBvh[i]=qnan;
		event.zBvh[i]=qnan;
		event.pxBvh[i]=qnan;
		event.pyBvh[i]=qnan;
		event.pzBvh[i]=qnan;
		event.ppBvh[i]=qnan;
		event.deBvh[i]=qnan;
		event.tBvh[i]=qnan;
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


  if( ievent%10==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1( 1, event.status++ );

  event.evnum = src.evnum;
  event.nhittpc = src.nhittpc;
  //debug  std::cout<<"DstTPCTracking Helix Geant4, nhit="<<event.nhittpc<<std::endl;
 	event.NumberOfTracks = src.NumberOfTracks; 
	for(int it=0;it<1000;++it){
		event.PIDOfTrack[it]=src.PIDOfTrack[it];
		event.ParentIDOfTrack[it]=src.ParentIDOfTrack[it];
		event.VertexOfTrack_x[it]=src.VertexOfTrack_x[it];
		event.VertexOfTrack_y[it]=src.VertexOfTrack_y[it];
		event.VertexOfTrack_z[it]=src.VertexOfTrack_z[it];
		event.MomentumOfTrack[it]=src.MomentumOfTrack[it];
		event.MomentumOfTrack_x[it]=src.MomentumOfTrack_x[it];
		event.MomentumOfTrack_y[it]=src.MomentumOfTrack_y[it];
		event.MomentumOfTrack_z[it]=src.MomentumOfTrack_z[it];
	}
	for(int it=0;it<event.NumberOfTracks;++it){
		int tid = it + 1;
		int pid = event.PIDOfTrack[tid];
		int parent = event.ParentIDOfTrack[tid];
		int pid_parent = 0;
		if(parent!= -1) pid_parent = event.PIDOfTrack[parent];
		if(pid == 11) continue;
		double px = event.MomentumOfTrack_x[tid];
		double py = event.MomentumOfTrack_y[tid];
		double pz = event.MomentumOfTrack_z[tid];
		double p = hypot(pz,hypot(px,py));
		if(pid == 321){
			if(pz < 0){
				event.PKm = p;
				event.PKm_x = px;
				event.PKm_y = py;
				event.PKm_z = pz;
			}
			else{
				event.PKp = p;
				event.PKp_x = px;
				event.PKp_y = py;
				event.PKp_z = pz;
			}
		}
		if(pid == 3312){
				event.PXi = p;
				event.PXi_x = px;
				event.PXi_y = py;
				event.PXi_z = pz;
		}
		if(pid == 3122 and pid_parent == 3312){
				event.PLd = p;
				event.PLd_x = px;
				event.PLd_y = py;
				event.PLd_z = pz;
		}
	}
	event.MomXi_x = src.MomXi_x;
	event.MomXi_y = src.MomXi_y;
	event.MomXi_z = src.MomXi_z;
	event.SpinXi_x = src.SpinXi_x;
	event.SpinXi_y = src.SpinXi_y;
	event.SpinXi_z = src.SpinXi_z;
	event.ThXi_CM = src.ThXi_CM;
	event.MomLd_x = src.MomLd_x;
	event.MomLd_y = src.MomLd_y;
	event.MomLd_z = src.MomLd_z;
	event.SpinLd_x = src.SpinLd_x;
	event.SpinLd_y = src.SpinLd_y;
	event.SpinLd_z = src.SpinLd_z;
	event.ThLd_CM = src.ThLd_CM;
	event.nhHtof = src.nhHtof;	
	for(int ihit=0;ihit<event.nhHtof;++ihit){
		event.tidHtof[ihit]=src.tidHtof[ihit];
		event.pidHtof[ihit]=src.pidHtof[ihit];
		event.didHtof[ihit]=src.didHtof[ihit];
		event.prtHtof[ihit]=src.prtHtof[ihit];
		event.qHtof[ihit]=src.qHtof[ihit];
		event.xHtof[ihit]=src.xHtof[ihit];
		event.yHtof[ihit]=src.yHtof[ihit];
		event.zHtof[ihit]=src.zHtof[ihit];
		event.pxHtof[ihit]=src.pxHtof[ihit];
		event.pyHtof[ihit]=src.pyHtof[ihit];
		event.pzHtof[ihit]=src.pzHtof[ihit];
		event.ppHtof[ihit]=src.ppHtof[ihit];
		event.deHtof[ihit]=src.deHtof[ihit];
		event.tHtof[ihit]=src.tHtof[ihit];
	}
	
	event.nhSch = src.nhSch;	
	for(int ihit=0;ihit<event.nhSch;++ihit){
		event.tidSch[ihit]=src.tidSch[ihit];
		event.pidSch[ihit]=src.pidSch[ihit];
		event.didSch[ihit]=src.didSch[ihit];
		event.prtSch[ihit]=src.prtSch[ihit];
		event.qSch[ihit]=src.qSch[ihit];
		event.xSch[ihit]=src.xSch[ihit];
		event.ySch[ihit]=src.ySch[ihit];
		event.zSch[ihit]=src.zSch[ihit];
		event.pxSch[ihit]=src.pxSch[ihit];
		event.pySch[ihit]=src.pySch[ihit];
		event.pzSch[ihit]=src.pzSch[ihit];
		event.ppSch[ihit]=src.ppSch[ihit];
		event.deSch[ihit]=src.deSch[ihit];
		event.tSch[ihit]=src.tSch[ihit];
	}
	event.nhLac = src.nhLac;	
	for(int ihit=0;ihit<event.nhLac;++ihit){
		event.tidLac[ihit]=src.tidLac[ihit];
		event.pidLac[ihit]=src.pidLac[ihit];
		event.didLac[ihit]=src.didLac[ihit];
		event.prtLac[ihit]=src.prtLac[ihit];
		event.qLac[ihit]=src.qLac[ihit];
		event.xLac[ihit]=src.xLac[ihit];
		event.yLac[ihit]=src.yLac[ihit];
		event.zLac[ihit]=src.zLac[ihit];
		event.pxLac[ihit]=src.pxLac[ihit];
		event.pyLac[ihit]=src.pyLac[ihit];
		event.pzLac[ihit]=src.pzLac[ihit];
		event.ppLac[ihit]=src.ppLac[ihit];
		event.deLac[ihit]=src.deLac[ihit];
		event.tLac[ihit]=src.tLac[ihit];
	}
	event.nhFtof = src.nhFtof;	
	for(int ihit=0;ihit<event.nhFtof;++ihit){
		event.tidFtof[ihit]=src.tidFtof[ihit];
		event.pidFtof[ihit]=src.pidFtof[ihit];
		event.didFtof[ihit]=src.didFtof[ihit];
		event.prtFtof[ihit]=src.prtFtof[ihit];
		event.qFtof[ihit]=src.qFtof[ihit];
		event.xFtof[ihit]=src.xFtof[ihit];
		event.yFtof[ihit]=src.yFtof[ihit];
		event.zFtof[ihit]=src.zFtof[ihit];
		event.pxFtof[ihit]=src.pxFtof[ihit];
		event.pyFtof[ihit]=src.pyFtof[ihit];
		event.pzFtof[ihit]=src.pzFtof[ihit];
		event.ppFtof[ihit]=src.ppFtof[ihit];
		event.deFtof[ihit]=src.deFtof[ihit];
		event.tFtof[ihit]=src.tFtof[ihit];
	}
	event.nhSdc = src.nhSdc;	
	for(int ihit=0;ihit<event.nhSdc;++ihit){
		event.tidSdc[ihit]=src.tidSdc[ihit];
		event.pidSdc[ihit]=src.pidSdc[ihit];
		event.didSdc[ihit]=src.didSdc[ihit];
		event.prtSdc[ihit]=src.prtSdc[ihit];
		event.qSdc[ihit]=src.qSdc[ihit];
		event.xSdc[ihit]=src.xSdc[ihit];
		event.ySdc[ihit]=src.ySdc[ihit];
		event.zSdc[ihit]=src.zSdc[ihit];
		event.pxSdc[ihit]=src.pxSdc[ihit];
		event.pySdc[ihit]=src.pySdc[ihit];
		event.pzSdc[ihit]=src.pzSdc[ihit];
		event.ppSdc[ihit]=src.ppSdc[ihit];
		event.deSdc[ihit]=src.deSdc[ihit];
		event.tSdc[ihit]=src.tSdc[ihit];
	}
	
	event.nhBvh = src.nhBvh;	
	for(int ihit=0;ihit<event.nhBvh;++ihit){
		event.tidBvh[ihit]=src.tidBvh[ihit];
		event.pidBvh[ihit]=src.pidBvh[ihit];
		event.didBvh[ihit]=src.didBvh[ihit];
		event.prtBvh[ihit]=src.prtBvh[ihit];
		event.qBvh[ihit]=src.qBvh[ihit];
		event.xBvh[ihit]=src.xBvh[ihit];
		event.yBvh[ihit]=src.yBvh[ihit];
		event.zBvh[ihit]=src.zBvh[ihit];
		event.pxBvh[ihit]=src.pxBvh[ihit];
		event.pyBvh[ihit]=src.pyBvh[ihit];
		event.pzBvh[ihit]=src.pzBvh[ihit];
		event.ppBvh[ihit]=src.ppBvh[ihit];
		event.deBvh[ihit]=src.deBvh[ihit];
		event.tBvh[ihit]=src.tBvh[ihit];
	}

	for(int ihit=0; ihit<event.nhittpc; ++ihit){
    event.ititpc[ihit] = src.ititpc[ihit];

    //for debug
    //debug    std::cout<<"iti:"<<src.ititpc[ihit]<<", pos:"
    //debug     <<"("<<src.xtpc[ihit]<<", "
    //debug     <<src.ytpc[ihit]<<", "
    //debug     <<src.ztpc[ihit]<<")"<<std::endl;

    if(event.max_ititpc<src.ititpc[ihit])
      event.max_ititpc = src.ititpc[ihit];
    ++event.nhittpc_iti[src.ititpc[ihit]-1];
  }
  for(int iti=0; iti<event.max_ititpc; ++iti){
    // debug     std::cout<<"nhit_iti("<<iti<<"): "<<event.nhittpc_iti[iti]<<std::endl;
  }
  // debug     getchar();

  // double u = src.pxPrm[0]/src.pzPrm[0];
  // double v = src.pyPrm[0]/src.pzPrm[0];
  // double cost = 1./std::sqrt(1.+u*u+v*v);
  // double theta=std::acos(cost)*math::Rad2Deg();
  // if(theta>20.)
  //   return true;

  if(src.nhittpc<5)
    return true;


  DCAnalyzer *DCAna = new DCAnalyzer();

  //with stable resolution
  // for(int it=0; it<src.nhittpc; ++it){
  //   src.x0tpc[it] += gRandom->Gaus(0,0.3);
  //   src.y0tpc[it] += gRandom->Gaus(0,0.3);
  //   src.z0tpc[it] += gRandom->Gaus(0,0.3);
  // }
  //for test
  // for(int it=0; it<src.nhittpc; ++it){
  //    // src.xtpc[it] += gRandom->Gaus(0,0.1);
  //    // src.ztpc[it] += gRandom->Gaus(0,0.1);
  //   double xtpc_re = src.x0tpc[it] + (src.xtpc[it] - src.x0tpc[it])*4.;
  //   double ztpc_re = src.z0tpc[it] + (src.ztpc[it] - src.z0tpc[it])*4.;
  //   src.xtpc[it]=xtpc_re;
  //   src.ztpc[it]=ztpc_re;
  // }

  //  std::cout<<"nhittpc"<<src.nhittpc<<std::endl;

  if(IsWithRes){

    DCAna->DecodeTPCHitsGeant4(src.nhittpc,
			       src.xtpc, src.ytpc, src.ztpc, src.edeptpc, src.idtpc);
  }
  else{
    DCAna->DecodeTPCHitsGeant4(src.nhittpc,
      			       src.x0tpc, src.y0tpc, src.z0tpc, src.edeptpc, src.idtpc);
  }
  DCAna->TrackSearchTPCHelix();

  int ntTpc = DCAna->GetNTracksTPCHelix();
  if( 1000<ntTpc ){
    std::cout << "#W " << func_name << " "
      	      << "too many ntTpc " << ntTpc << "/" << 1000 << std::endl;
    ntTpc = 1000;
  }
  event.ntTpc = ntTpc;
  event.nhtrack.resize( ntTpc );
  event.chisqr.resize( ntTpc );
  event.helix_cx.resize( ntTpc );
  event.helix_cy.resize( ntTpc );
  event.helix_z0.resize( ntTpc );
  event.helix_r.resize( ntTpc );
  event.helix_dz.resize( ntTpc );
  event.dE.resize( ntTpc );
  event.dEdx.resize( ntTpc );
  event.dz_factor.resize( ntTpc );
  event.mom0_x.resize( ntTpc );
  event.mom0_y.resize( ntTpc );
  event.mom0_z.resize( ntTpc );
  event.mom0.resize( ntTpc );
  event.mom_x.resize( ntTpc );
  event.mom_y.resize( ntTpc );
  event.mom_z.resize( ntTpc );
  event.charge.resize( ntTpc );
  event.path.resize( ntTpc );
  event.pid.resize( ntTpc );
  
  event.hitlayer.resize( ntTpc );
  event.hitpos_x.resize( ntTpc );
  event.hitpos_y.resize( ntTpc );
  event.hitpos_z.resize( ntTpc );
  event.calpos_x.resize( ntTpc );
  event.calpos_y.resize( ntTpc );
  event.calpos_z.resize( ntTpc );
  event.resolution.resize( ntTpc );
  event.resolution_x.resize( ntTpc );
  event.resolution_y.resize( ntTpc );
  event.resolution_z.resize( ntTpc );
  event.residual.resize( ntTpc );
  event.residual_x.resize( ntTpc );
  event.residual_y.resize( ntTpc );
  event.residual_z.resize( ntTpc );
  event.helix_t.resize( ntTpc );
  event.pathhit.resize(ntTpc);
  event.alpha.resize(ntTpc);
  event.track_cluster_de.resize(ntTpc);
  event.track_cluster_mrow.resize(ntTpc);
  for( int it=0; it<ntTpc; ++it ){
    TPCLocalTrackHelix *tp= DCAna->GetTrackTPCHelix(it);
    if(!tp) continue;
    int nh=tp->GetNHit();
    double chisqr    = tp->GetChiSquare();

    double cx=tp->Getcx(), cy=tp->Getcy();
    double z0=tp->Getz0(), r=tp->Getr();
    double dz = tp->Getdz();
    TVector3 Mom0 = tp->GetMom0();

    event.nhtrack[it] = nh;
    event.chisqr[it] = chisqr;
    event.helix_cx[it] = cx;
    event.helix_cy[it] = cy;
    event.helix_z0[it] = z0;
    event.helix_r[it]  = r;
    event.helix_dz[it] = dz;
    event.dz_factor[it] = sqrt(1.+(pow(dz,2)));
    event.mom0_x[it] = Mom0.X();
    event.mom0_y[it] = Mom0.Y();
    event.mom0_z[it] = Mom0.Z();
    event.mom0[it] = Mom0.Mag();
    event.dE[it] = tp->GetTrackdE();
    event.dEdx[it] = tp->GetdEdx(truncatedMean);
    event.charge[it]= tp->GetCharge();
//		event.pid[it]=particleID;
    event.hitlayer[it].resize( nh );
    event.hitpos_x[it].resize( nh );
    event.hitpos_y[it].resize( nh );
    event.hitpos_z[it].resize( nh );
    event.calpos_x[it].resize( nh );
    event.calpos_y[it].resize( nh );
    event.calpos_z[it].resize( nh );
    event.mom_x[it].resize( nh );
    event.mom_y[it].resize( nh );
    event.mom_z[it].resize( nh );
    event.residual[it].resize( nh );
    event.residual_x[it].resize( nh );
    event.residual_y[it].resize( nh );
    event.residual_z[it].resize( nh );
    event.resolution_x[it].resize( nh );
    event.resolution_y[it].resize( nh );
    event.resolution_z[it].resize( nh );
    event.helix_t[it].resize( nh );
    event.pathhit[it].resize(nh);
    event.alpha[it].resize(nh);
    
		event.track_cluster_de[it].resize(nh);
    event.track_cluster_mrow[it].resize(nh);


    //debug     std::cout<<"nh:"<<nh<<std::endl;
    double min_t = 9999,max_t = -9999;
		for( int ih=0; ih<nh; ++ih ){
			TPCLTrackHit *hit = tp -> GetHitInOrder(ih);
      if(!hit) continue;
      int layerId = 0;
      layerId = hit->GetLayer();
      const TVector3& resi_vect = hit->GetResidualVect();
      const TVector3& res_vect = hit->GetResolutionVect();
      const TVector3& hitpos = hit->GetLocalHitPos();
      Double_t clde = hit->GetDe();
      Double_t mrow = hit->GetMRow();
			const TVector3& calpos = hit->GetLocalCalPosHelix();
      const TVector3& mom = hit->GetMomentumHelix(tp->GetCharge());
			double residual = hit->GetResidual();
			

      for( int ih2=0; ih2<src.nhittpc; ++ih2 ){
	TVector3 setpos;
	if(IsWithRes)
	  setpos = TVector3(src.xtpc[ih2], src.ytpc[ih2], src.ztpc[ih2]);
	else
	  setpos = TVector3(src.x0tpc[ih2], src.y0tpc[ih2], src.z0tpc[ih2]);

	TVector3 d = setpos - hitpos;
	if(fabs(d.Mag()<0.1)){
	  event.momg_x[it][ih] = src.pxtpc[ih2]*1000.;
	  event.momg_y[it][ih] = src.pytpc[ih2]*1000.;
	  event.momg_z[it][ih] = src.pztpc[ih2]*1000.;
	  double momg_mag = sqrt(src.pxtpc[ih2]*src.pxtpc[ih2]
				 +src.pytpc[ih2]*src.pytpc[ih2]
				 +src.pztpc[ih2]*src.pztpc[ih2])*1000.;
	  event.residual_px[it][ih] = mom.x() - src.pxtpc[ih2]*1000.;//MeV/c
	  event.residual_py[it][ih] = mom.y() - src.pytpc[ih2]*1000.;//MeV/c
	  event.residual_pz[it][ih] = mom.z() - src.pztpc[ih2]*1000.;//MeV/c
	  event.residual_p[it][ih] = mom.Mag() - momg_mag;//MeV/c

	  event.iti_g[it][ih] = src.ititpc[ih2];
	  break;
	}
      }

      event.hitlayer[it][ih] = (double)layerId;
      event.hitpos_x[it][ih] = hitpos.x();
      event.hitpos_y[it][ih] = hitpos.y();
      event.hitpos_z[it][ih] = hitpos.z();
      event.helix_t[it][ih] = hit->GetTheta();
 //    	if(event.helix_t[it][ih]<min_t) min_t = event.helix_t[it][ih];
//     	if(event.helix_t[it][ih]>max_t) max_t = event.helix_t[it][ih];
			event.calpos_x[it][ih] = calpos.x();
      event.calpos_y[it][ih] = calpos.y();
      event.calpos_z[it][ih] = calpos.z();
      event.mom_x[it][ih] = mom.x();
      event.mom_y[it][ih] = mom.y();
      event.mom_z[it][ih] = mom.z();
      event.residual[it][ih] = residual;
      event.residual_x[it][ih] = resi_vect.x();
      event.residual_y[it][ih] = resi_vect.y();
      event.residual_z[it][ih] = resi_vect.z();
      event.track_cluster_de[it][ih] = clde;
      event.track_cluster_mrow[it][ih] = mrow;
      event.alpha[it][ih] = tp->GetAlpha(ih);
      event.resolution_x[it][ih] = res_vect.x();
      event.resolution_y[it][ih] = res_vect.y();
      event.resolution_z[it][ih] = res_vect.z();
		}
    event.path[it] = r*(max_t - min_t);

  }
  //debug   std::cout<<"end events"<<std::endl;
  //debug   getchar();


#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[nhittpc]: "<<std::setw(2)<<src.nhittpc<<" "<<std::endl;
#endif

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


  HBTree( "tpc", "tree of DstTPC_g" );

  tree->Branch("status", &event.status, "status/I" );
  tree->Branch("evnum", &event.evnum, "evnum/I" );
  tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");

  tree->Branch("max_ititpc",&event.max_ititpc,"max_ititpc/I");
  tree->Branch("ititpc",event.ititpc,"ititpc[nhittpc]/I");
  tree->Branch("nhittpc_iti",event.nhittpc_iti,"nhittpc_iti[max_ititpc]/I");

  tree->Branch( "ntTpc", &event.ntTpc );
  tree->Branch( "nhtrack", &event.nhtrack );
  tree->Branch( "chisqr", &event.chisqr );
  tree->Branch( "helix_cx", &event.helix_cx );
  tree->Branch( "helix_cy", &event.helix_cy );
  tree->Branch( "helix_z0", &event.helix_z0 );
  tree->Branch( "helix_r", &event.helix_r );
  tree->Branch( "helix_dz", &event.helix_dz );
  tree->Branch( "dz_factor", &event.dz_factor );
  
	// Momentum at Y = 0
  tree->Branch( "mom0_x", &event.mom0_x );
  tree->Branch( "mom0_y", &event.mom0_y );
  tree->Branch( "mom0_z", &event.mom0_z );
  tree->Branch( "mom0", &event.mom0 );
  tree->Branch( "dE", &event.dE );
  tree->Branch( "dEdx", &event.dEdx );
  
	tree->Branch( "charge", &event.charge );
  tree->Branch( "path", &event.path );

  tree->Branch( "pid", &event.pid );
  tree->Branch( "hitlayer", &event.hitlayer );
  tree->Branch( "hitpos_x", &event.hitpos_x );
  tree->Branch( "hitpos_y", &event.hitpos_y );
  tree->Branch( "hitpos_z", &event.hitpos_z );
  tree->Branch( "calpos_x", &event.calpos_x );
  tree->Branch( "calpos_y", &event.calpos_y );
  tree->Branch( "calpos_z", &event.calpos_z );
  tree->Branch( "alpha", &event.alpha);
  tree->Branch( "track_cluster_de", &event.track_cluster_de);
  tree->Branch( "track_cluster_mrow", &event.track_cluster_mrow);

  
	tree->Branch("momg_x",event.momg_x,"momg_x[nttpc][64]/D");
  tree->Branch("momg_y",event.momg_y,"momg_y[nttpc][64]/D");
  tree->Branch("momg_z",event.momg_z,"momg_z[nttpc][64]/D");
  tree->Branch("iti_g",event.iti_g,"iti_g[nttpc][64]/I");

  tree->Branch( "helix_t", &event.helix_t );
  tree->Branch( "residual", &event.residual );
  tree->Branch( "residual_x", &event.residual_x );
  tree->Branch( "residual_y", &event.residual_y );
  tree->Branch( "residual_z", &event.residual_z );
  tree->Branch( "resolution_x", &event.resolution_x);
  tree->Branch( "resolution_y", &event.resolution_y);
  tree->Branch( "resolution_z", &event.resolution_z);
  tree->Branch("residual_p",event.residual_p,"residual_p[nttpc][64]/D");


  tree->Branch("nPrm",&src.nhPrm,"nPrm/I");
  tree->Branch("xPrm",src.xPrm,"xPrm[nPrm]/D");
  tree->Branch("yPrm",src.yPrm,"yPrm[nPrm]/D");
  tree->Branch("zPrm",src.zPrm,"zPrm[nPrm]/D");
  tree->Branch("pxPrm",src.pxPrm,"pxPrm[nPrm]/D");
  tree->Branch("pyPrm",src.pyPrm,"pyPrm[nPrm]/D");
  tree->Branch("pzPrm",src.pzPrm,"pzPrm[nPrm]/D");
	tree->Branch("ppPrm",src.ppPrm,"ppPrm[nPrm]/D");


	tree->Branch("xtpc",src.xtpc,"xtpc[nhittpc]/D");
  tree->Branch("ytpc",src.ytpc,"ytpc[nhittpc]/D");
  tree->Branch("ztpc",src.ztpc,"ztpc[nhittpc]/D");
 
	tree->Branch("x0tpc",src.x0tpc,"x0tpc[nhittpc]/D");
  tree->Branch("y0tpc",src.y0tpc,"y0tpc[nhittpc]/D");
  tree->Branch("z0tpc",src.z0tpc,"z0tpc[nhittpc]/D");
 
	tree->Branch("pxtpc",src.pxtpc,"pxtpc[nhittpc]/D");
  tree->Branch("pytpc",src.pytpc,"pytpc[nhittpc]/D");
  tree->Branch("pztpc",src.pztpc,"pztpc[nhittpc]/D");
  tree->Branch("pptpc",src.pptpc,"pptpc[nhittpc]/D");
	tree->Branch("idtpc",src.idtpc,"idtpc[nhittpc]/I");
	tree->Branch("parentid",src.parentID,"parentID[nhittpc]/I");
  
	tree->Branch("NumberOfTracks",&event.NumberOfTracks,"NumberOfTracks/I");
	tree->Branch("PIDOfTrack",event.PIDOfTrack,"PIDOfTrack[1000]/I");
	tree->Branch("ParentIDOfTrack",event.ParentIDOfTrack,"ParentIDOfTrack[1000]/I");
	tree->Branch("VertexOfTrack_x",event.VertexOfTrack_x,"VertexOfTrack_x[1000]/D");
	tree->Branch("VertexOfTrack_y",event.VertexOfTrack_y,"VertexOfTrack_y[1000]/D");
	tree->Branch("VertexOfTrack_z",event.VertexOfTrack_z,"VertexOfTrack_z[1000]/D");
	tree->Branch("MomentumOfTrack",event.MomentumOfTrack,"MomentumOfTrack[1000]/D");
	tree->Branch("MomentumOfTrack_x",event.MomentumOfTrack_x,"MomentumOfTrack_x[1000]/D");
	tree->Branch("MomentumOfTrack_y",event.MomentumOfTrack_y,"MomentumOfTrack_y[1000]/D");
	tree->Branch("MomentumOfTrack_z",event.MomentumOfTrack_z,"MomentumOfTrack_z[1000]/D");
	
	
	tree->Branch("MomXi_x",&event.MomXi_x,"MomXi_x/D");
	tree->Branch("MomXi_y",&event.MomXi_y,"MomXi_y/D");
	tree->Branch("MomXi_z",&event.MomXi_z,"MomXi_z/D");
	tree->Branch("SpinXi_x",&event.SpinXi_x,"SpinXi_x/D");
	tree->Branch("SpinXi_y",&event.SpinXi_y,"SpinXi_y/D");
	tree->Branch("SpinXi_z",&event.SpinXi_z,"SpinXi_z/D");
	tree->Branch("ThXi_CM",&event.ThXi_CM,"ThXi_CM/D");
	
	tree->Branch("MomLd_x",&event.MomLd_x,"MomLd_x/D");
	tree->Branch("MomLd_y",&event.MomLd_y,"MomLd_y/D");
	tree->Branch("MomLd_z",&event.MomLd_z,"MomLd_z/D");
	tree->Branch("SpinLd_x",&event.SpinLd_x,"SpinLd_x/D");
	tree->Branch("SpinLd_y",&event.SpinLd_y,"SpinLd_y/D");
	tree->Branch("SpinLd_z",&event.SpinLd_z,"SpinLd_z/D");
	tree->Branch("ThLd_CM",&event.ThLd_CM,"ThLd_CM/D");

	tree->Branch("nhHtof", &event.nhHtof,"nhHtof/I");
  tree->Branch("tidHtof", event.tidHtof,"tidHtof[500]/I");
  tree->Branch("pidHtof", event.pidHtof,"pidHtof[500]/I");
  tree->Branch("didHtof", event.didHtof,"didHtof[500]/I");
  tree->Branch("prtHtof", event.prtHtof,"prtHtof[500]/I");
  tree->Branch("qHtof", event.qHtof,"qHtof[500]/I");
  tree->Branch("xHtof", event.xHtof,"xHtof[500]/D");
  tree->Branch("yHtof", event.yHtof,"yHtof[500]/D");
  tree->Branch("zHtof", event.zHtof,"zHtof[500]/D");
  tree->Branch("pxHtof", event.pxHtof,"pxHtof[500]/D");
  tree->Branch("pyHtof", event.pyHtof,"pyHtof[500]/D");
  tree->Branch("pzHtof", event.pzHtof,"pzHtof[500]/D");
  tree->Branch("ppHtof", event.ppHtof,"ppHtof[500]/D");
  tree->Branch("deHtof", event.deHtof,"deHtof[500]/D");
  tree->Branch("tHtof", event.tHtof,"tHtof[500]/D");

  tree->Branch("nhSch", &event.nhSch,"nhSch/I");
  tree->Branch("tidSch", event.tidSch,"tidSch[500]/I");
  tree->Branch("pidSch", event.pidSch,"pidSch[500]/I");
  tree->Branch("didSch", event.didSch,"didSch[500]/I");
  tree->Branch("prtSch", event.prtSch,"prtSch[500]/I");
  tree->Branch("qSch", event.qSch,"qSch[500]/I");
  tree->Branch("xSch", event.xSch,"xSch[500]/D");
  tree->Branch("ySch", event.ySch,"ySch[500]/D");
  tree->Branch("zSch", event.zSch,"zSch[500]/D");
  tree->Branch("pxSch", event.pxSch,"pxSch[500]/D");
  tree->Branch("pySch", event.pySch,"pySch[500]/D");
  tree->Branch("pzSch", event.pzSch,"pzSch[500]/D");
  tree->Branch("ppSch", event.ppSch,"ppSch[500]/D");
  tree->Branch("deSch", event.deSch,"deSch[500]/D");
  tree->Branch("tSch", event.tSch,"tSch[500]/D");
 
  tree->Branch("nhLac", &event.nhLac,"nhLac/I");
  tree->Branch("tidLac", event.tidLac,"tidLac[500]/I");
  tree->Branch("pidLac", event.pidLac,"pidLac[500]/I");
  tree->Branch("didLac", event.didLac,"didLac[500]/I");
  tree->Branch("prtLac", event.prtLac,"prtLac[500]/I");
  tree->Branch("qLac", event.qLac,"qLac[500]/I");
  tree->Branch("xLac", event.xLac,"xLac[500]/D");
  tree->Branch("yLac", event.yLac,"yLac[500]/D");
  tree->Branch("zLac", event.zLac,"zLac[500]/D");
  tree->Branch("pxLac", event.pxLac,"pxLac[500]/D");
  tree->Branch("pyLac", event.pyLac,"pyLac[500]/D");
  tree->Branch("pzLac", event.pzLac,"pzLac[500]/D");
  tree->Branch("ppLac", event.ppLac,"ppLac[500]/D");
  tree->Branch("deLac", event.deLac,"deLac[500]/D");
  tree->Branch("tLac", event.tLac,"tLac[500]/D");

	tree->Branch("nhFtof", &event.nhFtof,"nhFtof/I");
  tree->Branch("tidFtof", event.tidFtof,"tidFtof[500]/I");
  tree->Branch("pidFtof", event.pidFtof,"pidFtof[500]/I");
  tree->Branch("didFtof", event.didFtof,"didFtof[500]/I");
  tree->Branch("prtFtof", event.prtFtof,"prtFtof[500]/I");
  tree->Branch("qFtof", event.qFtof,"qFtof[500]/I");
  tree->Branch("xFtof", event.xFtof,"xFtof[500]/D");
  tree->Branch("yFtof", event.yFtof,"yFtof[500]/D");
  tree->Branch("zFtof", event.zFtof,"zFtof[500]/D");
  tree->Branch("pxFtof", event.pxFtof,"pxFtof[500]/D");
  tree->Branch("pyFtof", event.pyFtof,"pyFtof[500]/D");
  tree->Branch("pzFtof", event.pzFtof,"pzFtof[500]/D");
  tree->Branch("ppFtof", event.ppFtof,"ppFtof[500]/D");
  tree->Branch("deFtof", event.deFtof,"deFtof[500]/D");
  tree->Branch("tFtof", event.tFtof,"tFtof[500]/D");
	
	tree->Branch("nhSdc", &event.nhSdc,"nhSdc/I");
  tree->Branch("tidSdc", event.tidSdc,"tidSdc[500]/I");
  tree->Branch("pidSdc", event.pidSdc,"pidSdc[500]/I");
  tree->Branch("didSdc", event.didSdc,"didSdc[500]/I");
  tree->Branch("prtSdc", event.prtSdc,"prtSdc[500]/I");
  tree->Branch("qSdc", event.qSdc,"qSdc[500]/I");
  tree->Branch("xSdc", event.xSdc,"xSdc[500]/D");
  tree->Branch("ySdc", event.ySdc,"ySdc[500]/D");
  tree->Branch("zSdc", event.zSdc,"zSdc[500]/D");
  tree->Branch("pxSdc", event.pxSdc,"pxSdc[500]/D");
  tree->Branch("pySdc", event.pySdc,"pySdc[500]/D");
  tree->Branch("pzSdc", event.pzSdc,"pzSdc[500]/D");
  tree->Branch("ppSdc", event.ppSdc,"ppSdc[500]/D");
  tree->Branch("deSdc", event.deSdc,"deSdc[500]/D");
  tree->Branch("tSdc", event.tSdc,"tSdc[500]/D");


	tree->Branch("nhBvh", &event.nhBvh,"nhBvh/I");
  tree->Branch("tidBvh", event.tidBvh,"tidBvh[500]/I");
  tree->Branch("pidBvh", event.pidBvh,"pidBvh[500]/I");
  tree->Branch("didBvh", event.didBvh,"didBvh[500]/I");
  tree->Branch("prtBvh", event.prtBvh,"prtBvh[500]/I");
  tree->Branch("qBvh", event.qBvh,"qBvh[500]/I");
  tree->Branch("xBvh", event.xBvh,"xBvh[500]/D");
  tree->Branch("yBvh", event.yBvh,"yBvh[500]/D");
  tree->Branch("zBvh", event.zBvh,"zBvh[500]/D");
  tree->Branch("pxBvh", event.pxBvh,"pxBvh[500]/D");
  tree->Branch("pyBvh", event.pyBvh,"pyBvh[500]/D");
  tree->Branch("pzBvh", event.pzBvh,"pzBvh[500]/D");
  tree->Branch("ppBvh", event.ppBvh,"ppBvh[500]/D");
  tree->Branch("deBvh", event.deBvh,"deBvh[500]/D");
  tree->Branch("tBvh", event.tBvh,"tBvh[500]/D");


  ////////// Bring Address From Dst
  TTreeCont[kTPCGeant]->SetBranchStatus("*", 0);
  TTreeCont[kTPCGeant]->SetBranchStatus("evnum",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("nhittpc",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("nhPrm",  1);

  TTreeCont[kTPCGeant]->SetBranchStatus("xPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pyPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ppPrm",  1);

  TTreeCont[kTPCGeant]->SetBranchStatus("ititpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("idtpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("parentID", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("NumberOfTracks", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("PIDOfTrack", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ParentIDOfTrack", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("VertexOfTrack_x", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("VertexOfTrack_y", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("VertexOfTrack_z", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomentumOfTrack", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomentumOfTrack_x", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomentumOfTrack_y", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomentumOfTrack_z", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("SpinXi_x", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("SpinXi_y", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("SpinXi_z", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ThXi_CM", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomXi_x", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomXi_y", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomXi_z", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("SpinLd_x", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("SpinLd_y", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("SpinLd_z", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ThLd_CM", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomLd_x", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomLd_y", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("MomLd_z", 1);
  
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
  TTreeCont[kTPCGeant]->SetBranchStatus("pptpc", 1);   // total mometum
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

  TTreeCont[kTPCGeant]->SetBranchStatus("nhSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tidSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("didSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("prtSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("qSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ySch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pySch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ppSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("deSch", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tSch", 1);


  TTreeCont[kTPCGeant]->SetBranchStatus("nhLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tidLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("didLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("prtLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("qLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pyLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("deLac", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tLac", 1);
  
	TTreeCont[kTPCGeant]->SetBranchStatus("nhHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tidHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("didHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("prtHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("qHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pyHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ppHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("deHtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tHtof", 1);


  TTreeCont[kTPCGeant]->SetBranchStatus("nhFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tidFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("didFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("prtFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("qFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pyFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ppFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("deFtof", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tFtof", 1);
  
	TTreeCont[kTPCGeant]->SetBranchStatus("nhSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tidSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("didSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("prtSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("qSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ySdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pySdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ppSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("deSdc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tSdc", 1);


  TTreeCont[kTPCGeant]->SetBranchStatus("nhBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tidBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pidBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("didBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("prtBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("qBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("xBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("yBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("zBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pxBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pyBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("pzBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("ppBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("deBvh", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("tBvh", 1);
  
	
	TTreeCont[kTPCGeant]->SetBranchAddress("evnum", &src.evnum);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhittpc", &src.nhittpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhPrm", &src.nhPrm);

	TTreeCont[kTPCGeant]->SetBranchAddress("NumberOfTracks",&src.NumberOfTracks);
	TTreeCont[kTPCGeant]->SetBranchAddress("PIDOfTrack",src.PIDOfTrack);
	TTreeCont[kTPCGeant]->SetBranchAddress("ParentIDOfTrack",src.ParentIDOfTrack);
	TTreeCont[kTPCGeant]->SetBranchAddress("VertexOfTrack_x",src.VertexOfTrack_x);
	TTreeCont[kTPCGeant]->SetBranchAddress("VertexOfTrack_y",src.VertexOfTrack_y);
	TTreeCont[kTPCGeant]->SetBranchAddress("VertexOfTrack_z",src.VertexOfTrack_z);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack",src.MomentumOfTrack);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack_x",src.MomentumOfTrack_x);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack_y",src.MomentumOfTrack_y);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomentumOfTrack_z",src.MomentumOfTrack_z);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomXi_x",&src.MomXi_x);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomXi_y",&src.MomXi_y);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomXi_z",&src.MomXi_z);
	TTreeCont[kTPCGeant]->SetBranchAddress("SpinXi_x",&src.SpinXi_x);
	TTreeCont[kTPCGeant]->SetBranchAddress("SpinXi_y",&src.SpinXi_y);
	TTreeCont[kTPCGeant]->SetBranchAddress("SpinXi_z",&src.SpinXi_z);
	TTreeCont[kTPCGeant]->SetBranchAddress("ThXi_CM",&src.ThXi_CM);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomLd_x",&src.MomLd_x);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomLd_y",&src.MomLd_y);
	TTreeCont[kTPCGeant]->SetBranchAddress("MomLd_z",&src.MomLd_z);
	TTreeCont[kTPCGeant]->SetBranchAddress("SpinLd_x",&src.SpinLd_x);
	TTreeCont[kTPCGeant]->SetBranchAddress("SpinLd_y",&src.SpinLd_y);
	TTreeCont[kTPCGeant]->SetBranchAddress("SpinLd_z",&src.SpinLd_z);
	TTreeCont[kTPCGeant]->SetBranchAddress("ThLd_CM",&src.ThLd_CM);


  TTreeCont[kTPCGeant]->SetBranchAddress("xPrm", src.xPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("yPrm", src.yPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("zPrm", src.zPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxPrm", src.pxPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyPrm", src.pyPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzPrm", src.pzPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppPrm", src.ppPrm);

  TTreeCont[kTPCGeant]->SetBranchAddress("ititpc", src.ititpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("idtpc", src.idtpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("parentID", src.parentID);
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
  TTreeCont[kTPCGeant]->SetBranchAddress("pptpc", src.pptpc);
  //TTreeCont[kTPCGeant]->SetBranchAddress("masstpc", src.masstpc);
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
  
	TTreeCont[kTPCGeant]->SetBranchAddress("nhHtof", &src.nhHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidHtof", src.tidHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidHtof", src.pidHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("didHtof", src.didHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtHtof", src.prtHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("qHtof", src.qHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("xHtof", src.xHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("yHtof", src.yHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("zHtof", src.zHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxHtof", src.pxHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyHtof", src.pyHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzHtof", src.pzHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppHtof", src.ppHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("deHtof", src.deHtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tHtof", src.tHtof);



  TTreeCont[kTPCGeant]->SetBranchAddress("nhSch", &src.nhSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidSch", src.tidSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidSch", src.pidSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("didSch", src.didSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtSch", src.prtSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("qSch", src.qSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("xSch", src.xSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("ySch", src.ySch);
  TTreeCont[kTPCGeant]->SetBranchAddress("zSch", src.zSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxSch", src.pxSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pySch", src.pySch);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzSch", src.pzSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppSch", src.ppSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("deSch", src.deSch);
  TTreeCont[kTPCGeant]->SetBranchAddress("tSch", src.tSch);

  TTreeCont[kTPCGeant]->SetBranchAddress("nhLac", &src.nhLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidLac", src.tidLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidLac", src.pidLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("didLac", src.didLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtLac", src.prtLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("qLac", src.qLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("xLac", src.xLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("yLac", src.yLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("zLac", src.zLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxLac", src.pxLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyLac", src.pyLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzLac", src.pzLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppLac", src.ppLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("deLac", src.deLac);
  TTreeCont[kTPCGeant]->SetBranchAddress("tLac", src.tLac);
	
  TTreeCont[kTPCGeant]->SetBranchAddress("nhFtof", &src.nhFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidFtof", src.tidFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidFtof", src.pidFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("didFtof", src.didFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtFtof", src.prtFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("qFtof", src.qFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("xFtof", src.xFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("yFtof", src.yFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("zFtof", src.zFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxFtof", src.pxFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyFtof", src.pyFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzFtof", src.pzFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppFtof", src.ppFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("deFtof", src.deFtof);
  TTreeCont[kTPCGeant]->SetBranchAddress("tFtof", src.tFtof);
  
	TTreeCont[kTPCGeant]->SetBranchAddress("nhSdc", &src.nhSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidSdc", src.tidSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidSdc", src.pidSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("didSdc", src.didSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtSdc", src.prtSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("qSdc", src.qSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("xSdc", src.xSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ySdc", src.ySdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("zSdc", src.zSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxSdc", src.pxSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pySdc", src.pySdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzSdc", src.pzSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppSdc", src.ppSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("deSdc", src.deSdc);
  TTreeCont[kTPCGeant]->SetBranchAddress("tSdc", src.tSdc);

  TTreeCont[kTPCGeant]->SetBranchAddress("nhBvh", &src.nhBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("tidBvh", src.tidBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pidBvh", src.pidBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("didBvh", src.didBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("prtBvh", src.prtBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("qBvh", src.qBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("xBvh", src.xBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("yBvh", src.yBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("zBvh", src.zBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxBvh", src.pxBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyBvh", src.pyBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzBvh", src.pzBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppBvh", src.ppBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("deBvh", src.deBvh);
  TTreeCont[kTPCGeant]->SetBranchAddress("tBvh", src.tBvh);

	return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<TPCParamMan>("TPCPRM") &&
	  InitializeParameter<HodoPHCMan>("HDPHC") );
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
