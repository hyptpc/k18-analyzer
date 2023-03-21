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
#include "TPCPID.hh"

#include "DstHelper.hh"
#include <TRandom.h>
#include "DebugCounter.hh"

#include "TPCRiemannTrack.hh"
#include "TPCPatternRecognition.hh"
#include "TPCTrack.hh"
#include "TPCVertex.hh"

#include "HypTPCFitter.hh"
#include "HypTPCTask.hh"
#include "HypTPCVertex.hh"

namespace
{
using namespace root;
using namespace dst;
const std::string& class_name("DstTPCTrackingRiemannGeant4");
using hddaq::unpacker::GUnpacker;
const auto& gUnpacker = GUnpacker::get_instance();
ConfMan&            gConf = ConfMan::GetInstance();
const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
const UserParamMan& gUser = UserParamMan::GetInstance();
const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

const Int_t MaxTPCHits = 1000;
const Int_t MaxTPCTracks = 100;
const Int_t MaxTPCnHits = 50;

const auto& valueHall = ConfMan::Get<Double_t>("HSFLDHALL");
const auto& valueCalc = ConfMan::Get<Double_t>("HSFLDCALC");
  
  //const bool IsWithRes = false;
const bool IsWithRes = true;

//For GenFit Setting
//const bool Const_field = true;
const bool Const_field = false;
const Int_t verbosity = 0;//3;
const auto& tpcGeo = ConfMan::Get<TString>("TPCGDML");

  const Double_t truncation = 0.7;
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
  Int_t nhittpc;                 // Number of Hits
  Int_t max_ititpc;
  Int_t ititpc[MaxTPCHits];
  Int_t nhittpc_iti[MaxTPCHits];
  
  Int_t nhPrm;
  Double_t xPrm[MaxTPCTracks];
  Double_t yPrm[MaxTPCTracks];
  Double_t zPrm[MaxTPCTracks];
  Double_t pxPrm[MaxTPCTracks];
  Double_t pyPrm[MaxTPCTracks];
  Double_t pzPrm[MaxTPCTracks];
  Double_t ppPrm[MaxTPCTracks];
  Double_t thetaPrm[MaxTPCTracks];
  Double_t phiPrm[MaxTPCTracks];

  Double_t xtpc[MaxTPCHits];
  Double_t ytpc[MaxTPCHits];
  Double_t ztpc[MaxTPCHits];
  Double_t x0tpc[MaxTPCHits];
  Double_t y0tpc[MaxTPCHits];
  Double_t z0tpc[MaxTPCHits];

  Double_t pxtpc[MaxTPCHits];
  Double_t pytpc[MaxTPCHits];
  Double_t pztpc[MaxTPCHits];
  Double_t pptpc[MaxTPCHits];

  Double_t edeptpc[MaxTPCHits];
  Int_t laytpc[MaxTPCHits];
  Int_t rowtpc[MaxTPCHits];

  Int_t ntTpc; // Number of Tracks
  Int_t nhit[MaxTPCTracks]; // Number of Hits (in 1 tracks)
  Int_t isBeam[MaxTPCTracks]; // isBeam: 1 = Beam, 0 = Scat
  Double_t radius[MaxTPCTracks];
  Double_t xcenter[MaxTPCTracks];
  Double_t zcenter[MaxTPCTracks];
  Double_t yoffset[MaxTPCTracks];
  Double_t slope[MaxTPCTracks];
  Double_t dipangle[MaxTPCTracks];
  Double_t alpha[MaxTPCTracks];
  Double_t alphaHead[MaxTPCTracks];
  Double_t alphaTail[MaxTPCTracks];
  Double_t helixmom[MaxTPCTracks];
  Double_t tracklength[MaxTPCTracks];
  Double_t rmsH[MaxTPCTracks];
  Double_t rmsW[MaxTPCTracks];
  Double_t charge[MaxTPCTracks];
  Double_t xPOCA[MaxTPCTracks];
  Double_t yPOCA[MaxTPCTracks];
  Double_t zPOCA[MaxTPCTracks];
  Double_t closestDist[MaxTPCTracks];
  Double_t xtgtDir[MaxTPCTracks];
  Double_t ytgtDir[MaxTPCTracks];
  Double_t ztgtDir[MaxTPCTracks];

  Int_t ncluster[MaxTPCTracks]; // Number of Clusterss (in 1 tracks)

  std::vector<std::vector<Double_t>>    xhitpos;
  std::vector<std::vector<Double_t>>    yhitpos;
  std::vector<std::vector<Double_t>>    zhitpos;
  std::vector<std::vector<Double_t>>    hittime;
  std::vector<std::vector<Double_t>>    layerhit;
  std::vector<std::vector<Double_t>>    rowhit;

  std::vector<std::vector<Double_t>>    isStable;
  std::vector<std::vector<Double_t>>    xclusterpos;
  std::vector<std::vector<Double_t>>    yclusterpos;
  std::vector<std::vector<Double_t>>    zclusterpos;
  std::vector<std::vector<Double_t>>    layercluster;
  std::vector<std::vector<Double_t>>    nhitcluster;
  std::vector<std::vector<Double_t>>    clusterxdir;
  std::vector<std::vector<Double_t>>    clusterydir;
  std::vector<std::vector<Double_t>>    clusterzdir;
  
  Double_t dEdx[MaxTPCTracks];
  Double_t cdEdx[MaxTPCTracks];
  Double_t dEdx_truncated[MaxTPCTracks];
  Double_t cdEdx_truncated[MaxTPCTracks];

  //GenFit outputs
  Int_t nGFTracks;
  std::vector<Int_t> GFStatus;
  std::vector<Int_t> HelixID;
  std::vector<Int_t> GenFitID;
  std::vector<Int_t> VertexID;
  std::vector<Int_t> GFCharge;

  std::vector<Double_t> GFChisqr;
  std::vector<Double_t> GFpValue;

  std::vector<Double_t> vtxPx;
  std::vector<Double_t> vtxPy;
  std::vector<Double_t> vtxPz;
  
  std::vector<Double_t> vtxPosx;
  std::vector<Double_t> vtxPosy;
  std::vector<Double_t> vtxPosz;
  
  std::vector<Double_t> htofPx;
  std::vector<Double_t> htofPy;
  std::vector<Double_t> htofPz;
  
  std::vector<Double_t> htofPosx;
  std::vector<Double_t> htofPosy;
  std::vector<Double_t> htofPosz;

  std::vector<Double_t> GFTrackLength;
  std::vector<Double_t> GFTrackLengthHTOF;

  std::vector<Double_t> GFTrackTOF;
  std::vector<Double_t> GFTrackTOFHTOF;
  
  //GenFit outputs
  Int_t numVertices;
  std::vector<Int_t> RaveVtxID;
  std::vector<Int_t> VtxNumTrack;
  std::vector<Double_t> RaveVtxPosx;
  std::vector<Double_t> RaveVtxPosy;
  std::vector<Double_t> RaveVtxPosz;
  std::vector<Double_t> RaveChisqr;
  std::vector<Double_t> RaveVtxCovx;
  std::vector<Double_t> RaveVtxCovy;
  std::vector<Double_t> RaveVtxCovz;
  
 
  void clear() {
    status   = 0;
    evnum = 0;
    nhittpc = 0;
    max_ititpc = 0;
    nhPrm = 0.;
    
    for(int i=0; i<MaxTPCHits; ++i){
      xtpc[i] = -9999;
      ytpc[i] = -9999;
      ztpc[i] = -9999;
      x0tpc[i] = -9999;
      y0tpc[i] = -9999;
      z0tpc[i] = -9999;
      pxtpc[i] = -9999;
      pytpc[i] = -9999;
      pztpc[i] = -9999;
      
      edeptpc[i] = -9999;

      laytpc[i] = -9999;
      rowtpc[i] = -9999;
    }
    
    ntTpc = 0;

    for (int i = 0 ; i < MaxTPCTracks; i++) {
      xPrm[i] = -9999;
      yPrm[i] = -9999;
      zPrm[i] = -9999;

      pxPrm[i] = -9999;
      pyPrm[i] = -9999;
      pzPrm[i] = -9999;
      ppPrm[i] = -9999;

      thetaPrm[i] = -9999;
      phiPrm[i] = -9999;

      nhit[i] = -9999;
      isBeam[i] -= -9999;
      radius[i] = -9999;
      xcenter[i] = -9999;
      zcenter[i] = -9999;
      yoffset[i] = -9999;
      slope[i] = -9999;
      dipangle[i] = -9999;
      alpha[i] = -9999;
      alphaHead[i] = -9999;
      alphaTail[i] = -9999;
      helixmom[i] = -9999;

      tracklength[i] = -9999;
      rmsH[i] = -9999;
      rmsW[i] = -9999;
      charge[i] = -9999;

      xPOCA[i] = -9999;
      yPOCA[i] = -9999;
      zPOCA[i] = -9999;

      closestDist[i] = -9999;

      xtgtDir[i] = -9999;
      ytgtDir[i] = -9999;
      ztgtDir[i] = -9999;

      dEdx[i] = -9999;
      cdEdx[i] = -9999;
      dEdx_truncated[i] = -9999;
      cdEdx_truncated[i] = -9999;

      ncluster[i] = -9999;

    }
    
    xhitpos.clear();
    yhitpos.clear();
    zhitpos.clear();
    hittime.clear();
    layerhit.clear();
    rowhit.clear();

    isStable.clear();
    xclusterpos.clear();
    yclusterpos.clear();
    zclusterpos.clear();
    layercluster.clear();
    nhitcluster.clear();
    clusterxdir.clear();
    clusterydir.clear();
    clusterzdir.clear();


    nGFTracks = 0;
    GFStatus.clear();
    HelixID.clear();
    GenFitID.clear();
    VertexID.clear();
    GFCharge.clear();

    GFChisqr.clear();
    GFpValue.clear();

    vtxPx.clear();
    vtxPy.clear();
    vtxPz.clear();
    
    vtxPosx.clear();
    vtxPosy.clear();
    vtxPosz.clear();
    
    htofPx.clear();
    htofPy.clear();
    htofPz.clear();
    
    htofPosx.clear();
    htofPosy.clear();
    htofPosz.clear();

    GFTrackLength.clear();
    GFTrackLengthHTOF.clear();
    GFTrackTOF.clear();
    GFTrackTOFHTOF.clear();

    numVertices = 0;
    RaveVtxID.clear();
    VtxNumTrack.clear();
    RaveChisqr.clear();
    RaveVtxPosx.clear();
    RaveVtxPosy.clear();
    RaveVtxPosz.clear();
    RaveVtxCovx.clear();
    RaveVtxCovy.clear();
    RaveVtxCovz.clear();
    
    
  }
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
  Double_t thetaPrm[MaxTPCTracks];
  Double_t phiPrm[MaxTPCTracks];

  Int_t ititpc[MaxTPCHits];
  Int_t idtpc[MaxTPCHits];
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
  // Double_t dedxtpc[MaxTPCHits];
  // Double_t slengthtpc[MaxTPCHits];
  Int_t laytpc[MaxTPCHits];
  Int_t rowtpc[MaxTPCHits];
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

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries( TTreeCont );
  if (max_loop > 0) nevent = skip + max_loop;
  
  CatchSignal::Set();

  HypTPCFitter* fitter = new HypTPCFitter(tpcGeo.Data(),Const_field);
  //Initiallize the genfit track container
  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  GFtracks.SetVerbosity(verbosity);
  std::cout<<"GenFit verbosity = "<<"-1: Silent, 0: Minimum, 1: Errors only, 2: Errors and Warnings, 3: Verbose mode, long term debugging(default)"<<std::endl;
  std::cout<<"Current verbosity = "<<GFtracks.GetVerbosity()<<std::endl;
  
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
Bool_t
dst::InitializeEvent()
{
  event.clear();

  return true;
}

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

  if( ievent%1==0 ){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1( 1, event.status++ );
  
  event.evnum = src.evnum;
  event.nhittpc = src.nhittpc;
  event.nhPrm = src.nhPrm;

  for (int nt = 0 ; nt < event.nhPrm ; nt++) {

    event.xPrm[nt] = src.xPrm[nt];
    event.yPrm[nt] = src.yPrm[nt];
    event.zPrm[nt] = src.zPrm[nt];

    event.pxPrm[nt] = src.pxPrm[nt];
    event.pyPrm[nt] = src.pyPrm[nt];
    event.pzPrm[nt] = src.pzPrm[nt];
    event.ppPrm[nt] = src.ppPrm[nt];
    //    std::cout << event.ppPrm[nt] << std::endl;
    event.thetaPrm[nt] = src.thetaPrm[nt];
    event.phiPrm[nt] = src.phiPrm[nt];
  }
  

  for(int ihit=0; ihit<event.nhittpc; ++ihit){
    event.ititpc[ihit] = src.ititpc[ihit];
    if(event.max_ititpc<src.ititpc[ihit])
      event.max_ititpc = src.ititpc[ihit];
    ++event.nhittpc_iti[src.ititpc[ihit]-1];

    event.xtpc[ihit] = src.xtpc[ihit];
    event.ytpc[ihit] = src.ytpc[ihit];
    event.ztpc[ihit] = src.ztpc[ihit];

    event.x0tpc[ihit] = src.x0tpc[ihit];
    event.y0tpc[ihit] = src.y0tpc[ihit];
    event.z0tpc[ihit] = src.z0tpc[ihit];

    event.pxtpc[ihit] = src.pxtpc[ihit];
    event.pytpc[ihit] = src.pytpc[ihit];
    event.pztpc[ihit] = src.pztpc[ihit];
    event.pptpc[ihit] = src.pptpc[ihit];

    event.edeptpc[ihit] = src.edeptpc[ihit];

    event.laytpc[ihit] = src.laytpc[ihit];
    event.rowtpc[ihit] = src.rowtpc[ihit];
    
  }

  if(src.nhittpc<3)
    return true;

  //  DCAnalyzer *DCAna = new DCAnalyzer();
  DCAnalyzer DCAna;
  HypTPCTask& GFtracks = HypTPCTask::GetInstance();
  GFtracks.Init();
  
  DCAna.DecodeTPCHitsGeant4(src.nhittpc,
  			    src.xtpc, src.ytpc, src.ztpc, src.edeptpc);
  //  DCAna.DecodeTPCHitsGeant4(src.nhittpc,
  //  src.x0tpc, src.y0tpc, src.z0tpc, src.edeptpc);

  DCAna.TrackSearchTPCRiemann();
  
  Int_t nTrack = DCAna.GetNTracksTPCRiemann();
  event.ntTpc = nTrack;
  
  if (nTrack == 0 || nTrack > MaxTPCTracks) return true;

  event.xhitpos.resize(nTrack);
  event.yhitpos.resize(nTrack);
  event.zhitpos.resize(nTrack);
  event.hittime.resize(nTrack);
  event.layerhit.resize(nTrack);
  event.rowhit.resize(nTrack);

  event.isStable.resize(nTrack);
  event.xclusterpos.resize(nTrack);
  event.yclusterpos.resize(nTrack);
  event.zclusterpos.resize(nTrack);
  event.layercluster.resize(nTrack);
  event.nhitcluster.resize(nTrack);
  event.clusterxdir.resize(nTrack);
  event.clusterydir.resize(nTrack);
  event.clusterzdir.resize(nTrack);

  for (int nt = 0 ; nt < nTrack ; nt++) {
    auto track = DCAna.GetTrackTPCRiemann(nt);

    TPCPID::PID bestPID = TPCPID::kNotDetermined;
    Double_t bestPIDprob = 0.;

    for (Int_t iPID = 0 ; iPID < NUMTPCPID; ++iPID) {
      TPCPID::PID pid = static_cast<TPCPID::PID>(iPID);

      Int_t pdg;
      if (pid == TPCPID::kPion)
	pdg = 211;
    }
    
    //    GFtracks.AddHelixTrack(src.idtpc[0], track);
    
    Double_t Radius = track->GetRadius();
    TVector3 Center = track->GetCenter();
    Double_t Slope = track->GetSlope();
    Double_t Offset = track->GetOffset();
    Double_t Dip = track->GetDip();
    
    Double_t Alpha = track->GetAlpha();
    Double_t AlphaHead = track->GetAlphaHead();
    Double_t AlphaTail = track->GetAlphaTail();

    Double_t helixMom = track->GetMomentum(valueHall);
    Double_t tracklength = track->TrackLength();
    
    Double_t rmsW = track->GetRMSW();
    Double_t rmsH = track->GetRMSH();

    Double_t charge = track->Charge();
    Double_t closestDist = track->GetClosestDist();
    TVector3 POCA = track->GetPOCA();
    
    event.radius[nt] = Radius;
    event.xcenter[nt] = Center.X();
    event.zcenter[nt] = Center.Z();
    event.slope[nt] = Slope;
    event.yoffset[nt] = Offset;
    event.dipangle[nt] = Dip;

    event.alpha[nt] = Alpha;
    event.alphaHead[nt] = AlphaHead;
    event.alphaTail[nt] = AlphaTail;

    event.helixmom[nt] = helixMom;
    event.tracklength[nt] = tracklength;

    event.rmsW[nt] = rmsW;
    event.rmsH[nt] = rmsH;

    event.charge[nt] = charge;
    event.xPOCA[nt] = POCA.X();
    event.yPOCA[nt] = POCA.Y();
    event.zPOCA[nt] = POCA.Z();

    event.closestDist[nt] = closestDist;

    TVector3 dirTgt = track->GetDirectionAtTgt();

    event.xtgtDir[nt] = dirTgt.X();
    event.ytgtDir[nt] = dirTgt.Y();
    event.ztgtDir[nt] = dirTgt.Z();
    
    auto hits = track->GetHits();  
    event.nhit[nt] = hits->size();

    if (hits->size() == 0) return true;
    
    auto clusters = track->GetClusters();

    event.xhitpos[nt].resize(hits->size());
    event.yhitpos[nt].resize(hits->size());
    event.zhitpos[nt].resize(hits->size());
    event.hittime[nt].resize(hits->size());
    event.layerhit[nt].resize(hits->size());
    event.rowhit[nt].resize(hits->size());

    for (int nhit = 0 ; nhit < hits->size() ; nhit++) {

      auto hitpos = hits->at(nhit)->GetPosition();
      auto hittime = hits->at(nhit)->GetCTime();
      auto hitlayer = hits->at(nhit)->GetLayer();
      auto hitrow = hits->at(nhit)->GetRow();
      
      event.xhitpos[nt][nhit] = hitpos.X();
      event.yhitpos[nt][nhit] = hitpos.Y();
      event.zhitpos[nt][nhit] = hitpos.Z();
      event.hittime[nt][nhit] = hittime;
      event.layerhit[nt][nhit] = hitlayer;
      event.rowhit[nt][nhit] = hitrow;
    }

    event.isStable[nt].resize(clusters->size());
    event.xclusterpos[nt].resize(clusters->size());
    event.yclusterpos[nt].resize(clusters->size());
    event.zclusterpos[nt].resize(clusters->size());
    event.layercluster[nt].resize(clusters->size());
    event.nhitcluster[nt].resize(clusters->size());
    event.clusterxdir[nt].resize(clusters->size());
    event.clusterydir[nt].resize(clusters->size());
    event.clusterzdir[nt].resize(clusters->size());

    Int_t numCluster = 0;

    for (int ncluster = 0 ; ncluster < clusters->size() ; ncluster++) {

      auto stable = clusters->at(ncluster)->IsStable();
      auto clusterpos = clusters->at(ncluster)->GetPosition();
      auto layercluster = clusters->at(ncluster)->GetLayer();
      auto nhitcluster = clusters->at(ncluster)->GetNumHits();
      auto clusterdir = clusters->at(ncluster)->GetDirection();
      auto xdir = clusterdir.X();
      auto ydir = clusterdir.Y();
      auto zdir = clusterdir.Z();
	
      event.isStable[nt][ncluster] = stable;
      event.xclusterpos[nt][ncluster] = clusterpos.X();
      event.yclusterpos[nt][ncluster] = clusterpos.Y();
      event.zclusterpos[nt][ncluster] = clusterpos.Z();
      event.layercluster[nt][ncluster] = layercluster;
      event.nhitcluster[nt][ncluster] = nhitcluster;
      event.clusterxdir[nt][ncluster] = xdir;
      event.clusterydir[nt][ncluster] = ydir;
      event.clusterzdir[nt][ncluster] = zdir;

      numCluster++;

    }      

    event.ncluster[nt] = numCluster;

    auto dEdxArray = track->GetdEdxArray();
    auto cdEdxArray = track->GetcdEdxArray();
    std::sort(dEdxArray->begin(), dEdxArray->end());
    std::sort(cdEdxArray->begin(), cdEdxArray->end());

    auto numdEdxArray = dEdxArray->size();
    auto numcdEdxArray = cdEdxArray->size();

    if (numdEdxArray == 0) continue;
    if (numcdEdxArray != numdEdxArray) continue;
    
    auto numdEdxArray_truncated = truncation * numdEdxArray;
    
    auto dEdx = 0.;
    auto cdEdx = 0.;
    auto dEdx_truncated = 0.;
    auto cdEdx_truncated = 0.;

    Int_t numPoint = 0;
    Int_t numPoint_truncated = 0;
    
    for (int idEdx = 0 ; idEdx < numdEdxArray ; idEdx++) {
      dEdx += dEdxArray->at(idEdx);
      cdEdx += cdEdxArray->at(idEdx);
      if (idEdx < numdEdxArray_truncated) {
	dEdx_truncated += dEdxArray->at(idEdx);
	cdEdx_truncated += cdEdxArray->at(idEdx);
	numPoint_truncated++;
      }
      numPoint++;
    }
    
    if (numPoint != 0) {
      event.dEdx[nt] = dEdx/numPoint;
      event.cdEdx[nt] = cdEdx/numPoint;
    }
    if (numPoint_truncated != 0) {
      event.dEdx_truncated[nt] = dEdx_truncated/numPoint_truncated;
      event.cdEdx_truncated[nt] = cdEdx_truncated/numPoint_truncated;
    }    
  }  
  
  DCAna.TrackSearchTPCGenFit(GFtracks,true);
  
  Int_t nGFTracks = DCAna.GetNTracksTPCGenFit();

  event.GFStatus.resize(nGFTracks);
  event.HelixID.resize(nGFTracks);
  event.GenFitID.resize(nGFTracks);
  event.VertexID.resize(nGFTracks);
  event.GFCharge.resize(nGFTracks);

  event.GFChisqr.resize(nGFTracks);
  event.GFpValue.resize(nGFTracks);

  event.vtxPx.resize(nGFTracks);
  event.vtxPy.resize(nGFTracks);
  event.vtxPz.resize(nGFTracks);
    
  event.vtxPosx.resize(nGFTracks);
  event.vtxPosy.resize(nGFTracks);
  event.vtxPosz.resize(nGFTracks);
    
  event.htofPx.resize(nGFTracks);
  event.htofPy.resize(nGFTracks);
  event.htofPz.resize(nGFTracks);
    
  event.htofPosx.resize(nGFTracks);
  event.htofPosy.resize(nGFTracks);
  event.htofPosz.resize(nGFTracks);

  event.GFTrackLength.resize(nGFTracks);
  event.GFTrackLengthHTOF.resize(nGFTracks);
  event.GFTrackTOF.resize(nGFTracks);
  event.GFTrackTOFHTOF.resize(nGFTracks);

  for (int igf = 0 ; igf < nGFTracks ; igf++) {
    auto tpcTrack = DCAna.GetTrackTPCGenFit(igf);
    auto gfTrack = tpcTrack->GetGenFitTrack();
    Int_t gfStatus = tpcTrack->GetStatus();
    event.GFStatus[igf] = gfStatus;
    if (gfTrack == NULL || gfStatus == 0) continue;

    Int_t helixID = tpcTrack->GetHelixTrackID();
    Int_t gfTrackID = tpcTrack->GetGFTrackID();
    Int_t gfCharge = tpcTrack->GetGFCharge();
    Int_t gfVtxID = tpcTrack->GetVertexID();

    event.HelixID[igf] = helixID;
    event.GenFitID[igf] = gfTrackID;
    event.VertexID[igf] = gfVtxID;
    event.GFCharge[igf] = gfCharge;
    
    Double_t chisqr = tpcTrack->GetChisqr();
    Double_t pValue = tpcTrack->GetPvalue();

    event.GFChisqr[igf] = chisqr;
    event.GFpValue[igf] = pValue;
    
    TVector3 vtxMom = tpcTrack->GetMomentum();
    TVector3 vtxPos = tpcTrack->GetVertex();

    event.vtxPx[igf] = vtxMom.X();
    event.vtxPy[igf] = vtxMom.Y();
    event.vtxPz[igf] = vtxMom.Z();

    event.vtxPosx[igf] = vtxPos.X();
    event.vtxPosy[igf] = vtxPos.Y();
    event.vtxPosz[igf] = vtxPos.Z();

    TVector3 htofMom = tpcTrack->GetMomentumHTOF();
    TVector3 htofPos = tpcTrack->GetPosAtHTOF();

    event.htofPx[igf] = htofMom.X();
    event.htofPy[igf] = htofMom.Y();
    event.htofPz[igf] = htofMom.Z();

    event.htofPosx[igf] = htofPos.X();
    event.htofPosy[igf] = htofPos.Y();
    event.htofPosz[igf] = htofPos.Z();

    Double_t gfTrackLength = tpcTrack->GetTrackLength();
    Double_t gfTrackLengthHTOF = tpcTrack->GetTrackLengthHTOF();

    event.GFTrackLength[igf] = gfTrackLength;
    event.GFTrackLengthHTOF[igf] = gfTrackLengthHTOF;

    Double_t gfTrackTOF = tpcTrack->GetTrackTOF();
    Double_t gfTrackTOFHTOF = tpcTrack->GetTrackTOFHTOF();
    
    event.GFTrackTOF[igf] = gfTrackTOF;
    event.GFTrackTOFHTOF[igf] = gfTrackTOFHTOF;

  }

  Int_t nVertices = DCAna.GetNVerticesTPCRave();
  event.numVertices = nVertices;
  
  event.RaveVtxID.resize(nVertices);
  event.VtxNumTrack.resize(nVertices);
  event.RaveChisqr.resize(nVertices);
  event.RaveVtxPosx.resize(nVertices);
  event.RaveVtxPosy.resize(nVertices);
  event.RaveVtxPosz.resize(nVertices);
  event.RaveVtxCovx.resize(nVertices);
  event.RaveVtxCovy.resize(nVertices);
  event.RaveVtxCovz.resize(nVertices);

  for (int iVtx = 0 ; iVtx < nVertices ; iVtx++) {
    auto vertex = DCAna.GetVertexTPCRave(iVtx);

    Int_t vertexID = vertex->GetVertexID();
    Int_t vtxNumTrack = vertex->GetNumTrack();
    TVector3 vtxPos = vertex->GetPosition();
    Double_t chisqr = vertex->GetChisqr();
    TMatrixDSym vtxCov = vertex->GetCovariance();

    event.RaveVtxID[iVtx] = vertexID;
    event.VtxNumTrack[iVtx] = vtxNumTrack;
    event.RaveChisqr[iVtx] = chisqr;
    event.RaveVtxPosx[iVtx] = vtxPos.X();
    event.RaveVtxPosy[iVtx] = vtxPos.Y();
    event.RaveVtxPosz[iVtx] = vtxPos.Z();
    event.RaveVtxCovx[iVtx] = vtxCov[0][0];
    event.RaveVtxCovy[iVtx] = vtxCov[1][1];
    event.RaveVtxCovz[iVtx] = vtxCov[2][2];
    
  }
    
#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[nhittpc]: "<<std::setw(2)<<src.nhittpc<<" "<<std::endl;
#endif
  
  //  delete DCAna;
  
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

  tree->Branch("nhPrm",&event.nhPrm,"nhPrm/I");
  tree->Branch("xPrm",event.xPrm,"xPrm[nhPrm]/D");
  tree->Branch("yPrm",event.yPrm,"yPrm[nhPrm]/D");
  tree->Branch("zPrm",event.zPrm,"zPrm[nhPrm]/D");
  tree->Branch("pxPrm",event.pxPrm,"pxPrm[nhPrm]/D");
  tree->Branch("pyPrm",event.pyPrm,"pyPrm[nhPrm]/D");
  tree->Branch("pzPrm",event.pzPrm,"pzPrm[nhPrm]/D");
  tree->Branch("ppPrm",event.ppPrm,"ppPrm[nhPrm]/D");
  tree->Branch("thetaPrm",event.thetaPrm,"thetaPrm[nhPrm]/D");
  tree->Branch("phiPrm",event.phiPrm,"phiPrm[nhPrm]/D");

  tree->Branch("nhittpc",&event.nhittpc,"nhittpc/I");
  tree->Branch("max_ititpc",&event.max_ititpc,"max_ititpc/I");  
  tree->Branch("ititpc",event.ititpc,"ititpc[nhittpc]/I");
  tree->Branch("nhittpc_iti",event.nhittpc_iti,"nhittpc_iti[max_ititpc]/I");

  tree->Branch("xtpc",event.xtpc,"xtpc[nhittpc]/D");
  tree->Branch("ytpc",event.ytpc,"ytpc[nhittpc]/D");
  tree->Branch("ztpc",event.ztpc,"ztpc[nhittpc]/D");
  tree->Branch("x0tpc",event.x0tpc,"x0tpc[nhittpc]/D");
  tree->Branch("y0tpc",event.y0tpc,"y0tpc[nhittpc]/D");
  tree->Branch("z0tpc",event.z0tpc,"z0tpc[nhittpc]/D");
  tree->Branch("pxtpc",event.xtpc,"pxtpc[nhittpc]/D");
  tree->Branch("pytpc",event.ytpc,"pytpc[nhittpc]/D");
  tree->Branch("pztpc",event.ztpc,"pztpc[nhittpc]/D");

  tree->Branch("edeptpc",event.edeptpc,"edeptpc[nhittpc]/D");
  tree->Branch("laytpc",event.laytpc,"laytpc[nhittpc]/D");
  tree->Branch("rowtpc",event.rowtpc,"rowtpc[nhittpc]/D");
  
  tree->Branch("ntTpc", &event.ntTpc );
  tree->Branch("nhit",&event.nhit,"nhit[ntTpc]/I");
  tree->Branch("ncluster",&event.ncluster,"ncluster[ntTpc]/I");
  tree->Branch("radius",&event.radius,"radius[ntTpc]/D");
  tree->Branch("xcenter",event.xcenter,"xcenter[ntTpc]/D");
  tree->Branch("zcenter",event.zcenter,"zcenter[ntTpc]/D");
  tree->Branch("slope",event.slope,"event.slope[ntTpc]/D");
  tree->Branch("yoffset",event.yoffset,"event.offset[ntTpc]/D");
  tree->Branch("dipangle",event.dipangle,"event.dipangle[ntTpc]/D");

  tree->Branch("alpha",event.alpha,"event.alpha[ntTpc]/D");
  tree->Branch("alphaHead",event.alphaHead,"event.alphaHead[ntTpc]/D");
  tree->Branch("alphaTail",event.alphaTail,"event.alphaTail[ntTpc]/D");

  tree->Branch("helixmom",event.helixmom,"event.helixmom[ntTpc]/D");
  tree->Branch("tracklength",event.tracklength,"tracklength[ntTpc]/D");

  tree->Branch("rmsW",event.rmsW,"rmsW[ntTpc]/D");
  tree->Branch("rmsH",event.rmsH,"rmsH[ntTpc]/D");

  tree->Branch("charge",event.charge,"charge[ntTpc]/D");
  tree->Branch("xPOCA",event.xPOCA,"xPOCA[ntTpc]/D");
  tree->Branch("yPOCA",event.yPOCA,"yPOCA[ntTpc]/D");
  tree->Branch("zPOCA",event.zPOCA,"zPOCA[ntTpc]/D");
  tree->Branch("closestDist",event.closestDist,"closestDist[ntTpc]/D");

  tree->Branch("xtgtDir",event.xtgtDir,"xtgtDir[ntTpc]/D");
  tree->Branch("ytgtDir",event.ytgtDir,"ytgtDir[ntTpc]/D");
  tree->Branch("ztgtDir",event.ztgtDir,"ztgtDir[ntTpc]/D");

  tree->Branch("xhitpos",&event.xhitpos);
  tree->Branch("yhitpos",&event.yhitpos);
  tree->Branch("zhitpos",&event.zhitpos);
  tree->Branch("hittime",&event.hittime);
  tree->Branch("layerhit",&event.layerhit);
  tree->Branch("rowhit",&event.rowhit);
  
  tree->Branch("isStable",&event.isStable);
  tree->Branch("xclusterpos",&event.xclusterpos);
  tree->Branch("yclusterpos",&event.yclusterpos);
  tree->Branch("zclusterpos",&event.zclusterpos);
  tree->Branch("layercluster",&event.layercluster);
  tree->Branch("nhitcluster",&event.nhitcluster);
  tree->Branch("clusterxdir",&event.clusterxdir);
  tree->Branch("clusterydir",&event.clusterydir);
  tree->Branch("clusterzdir",&event.clusterzdir);

  tree->Branch("dEdx",&event.dEdx,"dEdx[ntTpc]/D");
  tree->Branch("cdEdx",&event.cdEdx,"cdEdx[ntTpc]/D");
  tree->Branch("dEdx_truncated",&event.dEdx_truncated,"dEdx_truncated[ntTpc]/D");
  tree->Branch("cdEdx_truncated",&event.cdEdx_truncated,"cdEdx_truncated[ntTpc]/D");

    //GenFit fit results
  tree->Branch("nGFTracks",&event.nGFTracks,"nGFTracks/I");
  tree->Branch("GFStatus",&event.GFStatus);
  tree->Branch("HelixID",&event.HelixID);
  tree->Branch("GenFitID",&event.GenFitID);
  tree->Branch("GFCharge",&event.GFCharge);

  tree->Branch("GFChisqr",&event.GFChisqr);
  tree->Branch("GFpValue",&event.GFpValue);

  tree->Branch("vtxPx",&event.vtxPx);
  tree->Branch("vtxPy",&event.vtxPy);
  tree->Branch("vtxPz",&event.vtxPz);

  tree->Branch("vtxPosx",&event.vtxPosx);
  tree->Branch("vtxPosy",&event.vtxPosy);
  tree->Branch("vtxPosz",&event.vtxPosz);

  tree->Branch("htofPx",&event.htofPx);
  tree->Branch("htofPy",&event.htofPy);
  tree->Branch("htofPz",&event.htofPz);

  tree->Branch("htofPosx",&event.htofPosx);
  tree->Branch("htofPosy",&event.htofPosy);
  tree->Branch("htofPosz",&event.htofPosz);

  tree->Branch("GFTrackLength",&event.GFTrackLength);
  tree->Branch("GFTrackLengthHTOF",&event.GFTrackLengthHTOF);

  tree->Branch("GFTrackTOF",&event.GFTrackTOF);
  tree->Branch("GFTrackTOFHTOF",&event.GFTrackTOFHTOF);

  tree->Branch("numVertices",&event.numVertices);
  tree->Branch("RaveVtxID",&event.RaveVtxID);
  tree->Branch("VtxNumTrack",&event.VtxNumTrack);
  tree->Branch("RaveChisqr",&event.RaveChisqr);
  
  tree->Branch("RaveVtxPosx",&event.RaveVtxPosx);
  tree->Branch("RaveVtxPosy",&event.RaveVtxPosy);
  tree->Branch("RaveVtxPosz",&event.RaveVtxPosz);
  tree->Branch("RaveVtxCovx",&event.RaveVtxCovx);
  tree->Branch("RaveVtxCovy",&event.RaveVtxCovy);
  tree->Branch("RaveVtxCovz",&event.RaveVtxCovz);
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
  TTreeCont[kTPCGeant]->SetBranchStatus("thetaPrm",  1);
  TTreeCont[kTPCGeant]->SetBranchStatus("phiPrm",  1);

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
  TTreeCont[kTPCGeant]->SetBranchStatus("laytpc", 1);
  TTreeCont[kTPCGeant]->SetBranchStatus("rowtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("parentID", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("iPadtpc", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("xtpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("ytpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("ztpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dxtpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dytpc_pad", 1);
  // TTreeCont[kTPCGeant]->SetBranchStatus("dztpc_pad", 1);


  TTreeCont[kTPCGeant]->SetBranchAddress("evnum", &src.evnum);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhittpc", &src.nhittpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("nhPrm", &src.nhPrm);

  TTreeCont[kTPCGeant]->SetBranchAddress("xPrm", src.xPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("yPrm", src.yPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("zPrm", src.zPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pxPrm", src.pxPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pyPrm", src.pyPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("pzPrm", src.pzPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("ppPrm", src.ppPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("thetaPrm", src.thetaPrm);
  TTreeCont[kTPCGeant]->SetBranchAddress("phiPrm", src.phiPrm);

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
  TTreeCont[kTPCGeant]->SetBranchAddress("laytpc", src.laytpc);
  TTreeCont[kTPCGeant]->SetBranchAddress("rowtpc", src.rowtpc);
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

