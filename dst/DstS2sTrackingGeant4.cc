// -*- C++ -*-

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "HodoAnalyzer.hh"
#include "HodoHit.hh"
#include "MathTools.hh"
#include "RMAnalyzer.hh"
#include "S2sLib.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"

#include "DstHelper.hh"
#include <TRandom.h>
#include <TParticle.h>
#include "DebugCounter.hh"

#define TotCut   0
#define UseTOF   0
#define LocalTrackSearch 0

namespace
{
using namespace root;
using namespace dst;
const std::string& class_name("DstS2sTrackingGeant4");
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto qnan = TMath::QuietNaN();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& zTOF = gGeom.LocalZ("TOF");
debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();

const Int_t MaxTPCHits = 10000;
const Int_t MaxTPCTracks = 100;
const Int_t MaxTPCnHits = 50;
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kS2sGeant4, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[S2sGeant4]",
  "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "g4s2s", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{

  Int_t evnum;
  Double_t x0;
  Double_t y0;
  Double_t z0;
  Double_t xp0;
  Double_t yp0;
  Double_t phi0;
  Double_t theta0;
  Double_t p0;
  Double_t pB;

  Int_t nhTof;
  Double_t TofSeg[MaxHits];

  Int_t ntSdcIn;
  Int_t nlSdcIn;
  Int_t nhSdcIn[MaxHits];
  Double_t wposSdcIn[NumOfLayersSdcIn];
  Double_t chisqrSdcIn[MaxHits];
  Double_t u0SdcIn[MaxHits];
  Double_t v0SdcIn[MaxHits];
  Double_t x0SdcIn[MaxHits];
  Double_t y0SdcIn[MaxHits];

  Int_t ntSdcOut;
  Int_t nlSdcOut;
  Int_t nhSdcOut[MaxHits];
  Double_t wposSdcOut[NumOfLayersSdcOut];
  Double_t chisqrSdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];

  Int_t ntS2s;
  Int_t nlS2s;
  Int_t nhS2s[MaxHits];
  Double_t chisqrS2s[MaxHits];
  Double_t path[MaxHits];
  Double_t stof[MaxHits];
  Double_t pS2s[MaxHits];
  Double_t qS2s[MaxHits];
  Double_t m2[MaxHits];
  Double_t resP[MaxHits];
  Double_t vpx[NumOfLayersVP];
  Double_t vpy[NumOfLayersVP];
  Double_t vpu[NumOfLayersVP];
  Double_t vpv[NumOfLayersVP];

  Double_t xtgtS2s[MaxHits];
  Double_t ytgtS2s[MaxHits];
  Double_t utgtS2s[MaxHits];
  Double_t vtgtS2s[MaxHits];
  Double_t thetaS2s[MaxHits];
  Double_t phiS2s[MaxHits];

  Double_t xtofS2s[MaxHits];
  Double_t ytofS2s[MaxHits];
  Double_t utofS2s[MaxHits];
  Double_t vtofS2s[MaxHits];

  Double_t lxtofS2s[MaxHits];
  Double_t lytofS2s[MaxHits];
  Double_t lutofS2s[MaxHits];
  Double_t lvtofS2s[MaxHits];
  Double_t tofsegS2s[MaxHits];

  std::vector< std::vector<Double_t> > resL;
  std::vector< std::vector<Double_t> > resG;

  // Calib
  enum eParticle { Pion, Kaon, Proton, nParticle };
  Double_t tTofCalc[nParticle];
  Double_t utTofSeg[NumOfSegTOF];
  Double_t dtTofSeg[NumOfSegTOF];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];
  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];

};

//_____________________________________________________________________
struct Src
{
  TTreeReaderValue<Int_t>* evnum;
  TTreeReaderValue<Double_t>* x0;
  TTreeReaderValue<Double_t>* y0;
  TTreeReaderValue<Double_t>* z0;
  TTreeReaderValue<Double_t>* xp0;
  TTreeReaderValue<Double_t>* yp0;
  TTreeReaderValue<Double_t>* phi0;
  TTreeReaderValue<Double_t>* theta0;
  TTreeReaderValue<Double_t>* p0;
  TTreeReaderValue<Double_t>* pB;

  TTreeReaderArray<TParticle>* PRM;
  TTreeReaderArray<TParticle>* SDC1;
  TTreeReaderArray<TParticle>* SDC2;
  TTreeReaderArray<TParticle>* SDC3;
  TTreeReaderArray<TParticle>* SDC4;
  TTreeReaderArray<TParticle>* SDC5;
  TTreeReaderArray<TParticle>* TOF;
  TTreeReaderArray<TParticle>* AC1;
  TTreeReaderArray<TParticle>* WC;
  TTreeReaderArray<TParticle>* VP;
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
main(int argc, char **argv)
{
  std::vector<std::string> arg(argv, argv+argc);

  if(!CheckArg(arg))
    return EXIT_FAILURE;
  if(!DstOpen(arg))
    return EXIT_FAILURE;
  if(!gConf.Initialize(arg[kConfFile]))
    return EXIT_FAILURE;
  if(!gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  Int_t skip = gUnpacker.get_skip();
  if (skip < 0) skip = 0;
  Int_t max_loop = gUnpacker.get_max_loop();
  Int_t nevent = GetEntries(TTreeCont);
  if (max_loop > 0) nevent = skip + max_loop;

  CatchSignal::Set();

  Int_t ievent = skip;
  for(; ievent<nevent && !CatchSignal::Stop(); ++ievent){
    gCounter.check();
    InitializeEvent();
    if(DstRead(ievent)) tree->Fill();
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

  event.evnum = 0;
  event.x0 = qnan;
  event.y0 = qnan;
  event.z0 = qnan;
  event.xp0 = qnan;
  event.yp0 = qnan;
  event.phi0 = qnan;
  event.theta0 = qnan;
  event.p0 = qnan;
  event.pB = qnan;

  event.nhTof = 0;

  event.ntSdcIn  = 0;
  event.nlSdcIn  = 0;
  event.ntSdcOut = 0;
  event.nlSdcOut = 0;
  event.ntS2s = 0;
  event.nlS2s = 0;


  for(Int_t it=0; it<MaxHits; it++){
    event.TofSeg[it] = -1;
  }

  for(Int_t it=0; it<NumOfLayersSdcIn; ++it){
    event.wposSdcIn[it] = qnan;
  }
  for(Int_t it=0; it<NumOfLayersSdcOut; ++it){
    event.wposSdcOut[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    event.nhSdcIn[it] = 0;
    event.chisqrSdcIn[it] = qnan;
    event.x0SdcIn[it] = qnan;
    event.y0SdcIn[it] = qnan;
    event.u0SdcIn[it] = qnan;
    event.v0SdcIn[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    event.nhSdcOut[it] = 0;
    event.chisqrSdcOut[it] = qnan;
    event.x0SdcOut[it] = qnan;
    event.y0SdcOut[it] = qnan;
    event.u0SdcOut[it] = qnan;
    event.v0SdcOut[it] = qnan;
  }

  for(Int_t it=0; it<MaxHits; it++){
    event.nhS2s[it]     = 0;
    event.chisqrS2s[it] = qnan;
    event.stof[it]         = qnan;
    event.path[it]         = qnan;
    event.pS2s[it]      = qnan;
    event.qS2s[it]      = qnan;
    event.m2[it]           = qnan;
    event.xtgtS2s[it]  = qnan;
    event.ytgtS2s[it]  = qnan;
    event.utgtS2s[it]  = qnan;
    event.vtgtS2s[it]  = qnan;
    event.thetaS2s[it] = qnan;
    event.phiS2s[it]   = qnan;
    event.resP[it]        = qnan;
    event.xtofS2s[it]  = qnan;
    event.ytofS2s[it]  = qnan;
    event.utofS2s[it]  = qnan;
    event.vtofS2s[it]  = qnan;
    event.lxtofS2s[it]  = qnan;
    event.lytofS2s[it]  = qnan;
    event.lutofS2s[it]  = qnan;
    event.lvtofS2s[it]  = qnan;
    event.tofsegS2s[it] = qnan;
  }

  for(Int_t i=0; i<PlMaxTOF; ++i){
    event.resL[i].clear();
    event.resG[i].clear();
  }


  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstOpen(std::vector<std::string> arg)
{
  Int_t open_file = 0;
  Int_t open_tree = 0;
  for(std::size_t i=0; i<nArgc; ++i){
    if(i==kProcess || i==kConfFile || i==kOutFile) continue;
    open_file += OpenFile(TFileCont[i], arg[i]);
    open_tree += OpenTree(TFileCont[i], TTreeCont[i], TreeName[i]);
  }

  if(open_file!=open_tree || open_file!=nArgc-3)
    return false;
  if(!CheckEntries(TTreeCont))
    return false;

  TFileCont[kOutFile] = new TFile(arg[kOutFile].c_str(), "recreate");

  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstRead(Int_t ievent)
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");

  static const auto MinTotSDC1 = gUser.GetParameter("MinTotSDC1");
  static const auto MinTotSDC2 = gUser.GetParameter("MinTotSDC2");
  static const auto MinTotSDC3 = gUser.GetParameter("MinTotSDC3");
  static const auto MinTotSDC4 = gUser.GetParameter("MinTotSDC4");
  static const auto MinTotSDC5 = gUser.GetParameter("MinTotSDC5");
  static const auto MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
  static const auto MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  if(ievent%10000==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1(1, 0);

  event.evnum  = **src.evnum;
  event.x0     = **src.x0;
  event.y0     = **src.y0;
  event.z0     = **src.z0;
  event.xp0    = **src.xp0;
  event.yp0    = **src.yp0;
  event.phi0   = **src.phi0;
  event.theta0 = **src.theta0;
  event.p0     = **src.p0;
  event.pB     = **src.pB;

  event.nhTof = (*src.TOF).GetSize();
  for(Int_t i=0; i<event.nhTof; ++i){
    event.TofSeg[i] = (*src.TOF)[i].GetSecondMother();
  }

  DCAnalyzer DCAna;

  DCAna.DecodeSdcInHitsGeant4(*src.SDC1, *src.SDC2);
  DCAna.DecodeSdcOutHitsGeant4(*src.SDC3, *src.SDC4, *src.SDC5);

#if TotCut
  DCAna.TotCutSDC1(MinTotSDC1);
  DCAna.TotCutSDC2(MinTotSDC2);
  DCAna.TotCutSDC3(MinTotSDC3);
  DCAna.TotCutSDC4(MinTotSDC4);
  DCAna.TotCutSDC5(MinTotSDC5);
#endif

  HF1(1, 1);

  Double_t multi_SdcIn  = 0.;
  ////////////// SdcIn number of hit layer
  {
    Int_t nlSdcIn = 0;
    for(Int_t l=0; l<NumOfLayersSdcIn; ++l){
      const auto& contSdcIn =DCAna.GetSdcInHC(l);
      Int_t nhSdcIn = contSdcIn.size();
      if(nhSdcIn==1){
	auto hit = contSdcIn[0];
	Double_t wpos = hit->GetWirePosition();
	event.wposSdcIn[l] = wpos;
      }
      multi_SdcIn += Double_t(nhSdcIn);
      if(nhSdcIn>0) nlSdcIn++;
    }
    event.nlSdcIn   = nlSdcIn;
    event.nlS2s += nlSdcIn;
  }
  if(multi_SdcIn/Double_t(NumOfLayersSdcIn) > MaxMultiHitSdcIn){
    // return true;
  }

  Double_t multi_SdcOut = 0.;
  ////////////// SdcOut number of hit layer
  {
    Int_t nlSdcOut = 0;
    for(Int_t l=0; l<NumOfLayersSdcOut; ++l){
      const auto& contSdcOut =DCAna.GetSdcOutHC(l);
      Int_t nhSdcOut=contSdcOut.size();
      if(nhSdcOut==1){
	auto hit = contSdcOut[0];
	Double_t wpos = hit->GetWirePosition();
	event.wposSdcOut[l] = wpos;
      }
      multi_SdcOut += Double_t(nhSdcOut);
      if(nhSdcOut>0) nlSdcOut++;
    }
    event.nlSdcOut = nlSdcOut;
    event.nlS2s += nlSdcOut;
  }
  if(multi_SdcOut/Double_t(NumOfLayersSdcOut) > MaxMultiHitSdcOut){
    // return true;
  }

  HF1(1, 11.);

  // std::cout << "==========TrackSearch SdcIn============" << std::endl;
#if LocalTrackSearch
  DCAna.TrackSearchSdcIn();
#else
  DCAna.MakeTrackSdcInGeant4();
#endif
  // DCAna.ChiSqrCutSdcIn(50.);
  Int_t ntSdcIn = DCAna.GetNtracksSdcIn();
  if(MaxHits<ntSdcIn){
    std::cout << "#W too many ntSdcIn " << ntSdcIn << "/" << MaxHits << std::endl;
    ntSdcIn = MaxHits;
  }
  event.ntSdcIn = ntSdcIn;
  HF1(10, Double_t(ntSdcIn));
  for(Int_t it=0; it<ntSdcIn; ++it){
    const auto& track = DCAna.GetTrackSdcIn(it);
    Int_t nh=track->GetNHit();
    Double_t chisqr=track->GetChiSquare();
    Double_t x0=track->GetX0(), y0=track->GetY0();
    Double_t u0=track->GetU0(), v0=track->GetV0();
    HF1(11, Double_t(nh));
    HF1(12, chisqr);
    HF1(14, x0); HF1(15, y0);
    HF1(16, u0); HF1(17, v0);
    HF2(18, x0, u0); HF2(19, y0, v0);
    HF2(20, x0, y0);
    event.nhSdcIn[it] = nh;
    event.chisqrSdcIn[it] = chisqr;
    event.x0SdcIn[it] = x0;
    event.y0SdcIn[it] = y0;
    event.u0SdcIn[it] = u0;
    event.v0SdcIn[it] = v0;
    for(Int_t ih=0; ih<nh; ++ih){
      const auto& hit=track->GetHit(ih);

      Int_t layerId = hit->GetLayer();
      HF1(13, hit->GetLayer());
      Double_t wire = hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+1, wire-0.5);
      HF1(100*layerId+2, dt);
      HF1(100*layerId+3, dl);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1(100*layerId+4, pos);
      HF1(100*layerId+5, res);
      HF2(100*layerId+6, pos, res);
      HF2(100*layerId+7, xcal, ycal);
      event.resL[layerId-1].push_back(res);
    }
  }
  if(ntSdcIn<1) return true;
  //  if(!(ntSdcIn==1)) return true;

  HF1(1, 12.);

#if 0
  //////////////SdcIn vs Tof cut for Proton event
  {
    Int_t ntOk=0;
    for(Int_t it=0; it<ntSdcIn; ++it){
      const auto& track = DCAna.GetTrackSdcIn(it);
      if(!track) continue;

      Int_t nh=track->GetNHit();
      Double_t chisqr=track->GetChiSquare();
      Double_t u0=track->GetU0(), v0=track->GetV0();
      Double_t x0=track->GetX0(), y0=track->GetY0();

      Bool_t condTof=false;
      for(Int_t j=0; j<ncTof; ++j){
	HodoCluster *clTof=hodoAna.GetCluster("TOF", j);
	if(!clTof || !clTof->GoodForAnalysis()) continue;
	Double_t ttof=clTof->CMeanTime()-time0;
	//------------------------Cut
	if(MinModTimeTof< ttof+14.0*u0 && ttof+14.0*u0 <MaxModTimeTof)
	  condTof=true;
      }
      if(condTof){
	++ntOk;
	for(Int_t j=0; j<ncTof; ++j){
	  HodoCluster *clTof=hodoAna.GetCluster("TOF", j);
	  if(!clTof || !clTof->GoodForAnalysis()) continue;
	  Double_t ttof=clTof->CMeanTime()-time0;
	}
	// if(ntOk>0) track->GoodForTracking(false);
      }
      else {
	// track->GoodForTracking(false);
      }
    }
    // if(ntOk<1) return true;
  }
#endif

  HF1(1, 13.);

  //////////////SdcOut tracking
  //std::cout << "==========TrackSearch SdcOut============" << std::endl;

#if LocalTrackSearch

#if UseTOF
  DCAna.TrackSearchSdcOut(TOFCont);
#else
  DCAna.TrackSearchSdcOut();
#endif

#else
  DCAna.MakeTrackSdcOutGeant4();
#endif

  // DCAna.ChiSqrCutSdcOut(50.);
  Int_t ntSdcOut = DCAna.GetNtracksSdcOut();
  if(MaxHits<ntSdcOut){
    std::cout << "#W too many ntSdcOut " << ntSdcOut << "/" << MaxHits << std::endl;
    ntSdcOut = MaxHits;
  }
  event.ntSdcOut=ntSdcOut;
  HF1(40, Double_t(ntSdcOut));
  for(Int_t it=0; it<ntSdcOut; ++it){
    const auto& track = DCAna.GetTrackSdcOut(it);
    Int_t nh=track->GetNHit();
    Double_t chisqr=track->GetChiSquare();
    Double_t u0=track->GetU0(), v0=track->GetV0();
    Double_t x0=track->GetX0(), y0=track->GetY0();
    HF1(41, Double_t(nh));
    HF1(42, chisqr);
    HF1(44, x0); HF1(45, y0);
    HF1(46, u0); HF1(47, v0);
    HF2(48, x0, u0); HF2(49, y0, v0);
    HF2(50, x0, y0);
    event.nhSdcOut[it] = nh;
    event.chisqrSdcOut[it] = chisqr;
    event.x0SdcOut[it] = x0;
    event.y0SdcOut[it] = y0;
    event.u0SdcOut[it] = u0;
    event.v0SdcOut[it] = v0;
    for(Int_t ih=0; ih<nh; ++ih){
      const auto& hit = track->GetHit(ih);
      Int_t layerId = hit->GetLayer();
      HF1(43, hit->GetLayer());
      Double_t wire=hit->GetWire();
      Double_t dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+1, wire-0.5);
      HF1(100*layerId+2, dt);
      HF1(100*layerId+3, dl);
      Double_t xcal=hit->GetXcal(), ycal=hit->GetYcal();
      Double_t pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1(100*layerId+4, pos);
      HF1(100*layerId+5, res);
      HF2(100*layerId+6, pos, res);
      HF2(100*layerId+7, xcal, ycal);
      event.resL[layerId-1].push_back(res);
    }
  }

  if(ntSdcOut<1) return true;

  HF1(1, 14.);

  for(Int_t i1=0; i1<ntSdcIn; ++i1){
    const auto& trSdcIn=DCAna.GetTrackSdcIn(i1);
    Double_t xin=trSdcIn->GetX0(), yin=trSdcIn->GetY0();
    Double_t uin=trSdcIn->GetU0(), vin=trSdcIn->GetV0();
    for(Int_t i2=0; i2<ntSdcOut; ++i2){
      const auto& trSdcOut=DCAna.GetTrackSdcOut(i2);
      Double_t xout=trSdcOut->GetX0(), yout=trSdcOut->GetY0();
      Double_t uout=trSdcOut->GetU0(), vout=trSdcOut->GetV0();
      HF2(20001, xin, xout); HF2(20002, yin, yout);
      HF2(20003, uin, uout); HF2(20004, vin, vout);
    }
  }

  HF1(1, 20.);

  // if(ntSdcIn*ntSdcOut > 4) return true;

  HF1(1, 21.);

  ///// BTOF BH2-Target
  static const auto StofOffset = gUser.GetParameter("StofOffset");

  //////////////S2S Tracking
  if( false && event.nhTof > 0 ){
    Double_t seg = event.TofSeg[0];
    Double_t par[3] = {1.59, -3.07e-2, 4.19e-4};
    Double_t pMag   = par[0]+seg*par[1]+seg*seg*par[2]; // Magnitude
    Double_t scale  = 1.;
    Double_t initial_momentum = pMag*scale;
    DCAna.TrackSearchS2s(initial_momentum);
  }
  else{
    DCAna.TrackSearchS2s(event.p0/1000.); // [GeV/c]
    // DCAna.TrackSearchS2s(); // use PK18 as initial momentum
  }

  Int_t ntS2s = DCAna.GetNTracksS2s();
  if(MaxHits < ntS2s){
    std::cout << "#W too many ntS2s " << ntS2s << "/" << MaxHits << std::endl;
    ntS2s = MaxHits;
  }
  event.ntS2s = ntS2s;
  HF1(70, ntS2s);

  for(Int_t i=0; i<ntS2s; ++i){
    auto track = DCAna.GetS2sTrack(i);
    if(!track) continue;
    // track->Print();
    Int_t nh = track->GetNHits();
    Double_t chisqr = track->ChiSquare();
    const auto& Pos = track->PrimaryPosition();
    const auto& Mom = track->PrimaryMomentum();
    // hddaq::cout << std::fixed
    // 		<< "Pos = " << Pos << std::endl
    // 		<< "Mom = " << Mom << std::endl;
    Double_t path = track->PathLengthToTOF();
    Double_t xt = Pos.x(), yt = Pos.y();
    Double_t p = Mom.Mag();
    Double_t q = track->Polarity();
    Double_t ut = Mom.x()/Mom.z(), vt = Mom.y()/Mom.z();
    Double_t cost = 1./TMath::Sqrt(1.+ut*ut+vt*vt);
    Double_t theta = TMath::ACos(cost)*TMath::RadToDeg();
    Double_t phi = TMath::ATan2(ut, vt);
    Double_t initial_momentum = track->GetInitialMomentum();
    HF1(71, nh);
    HF1(72, chisqr);
    HF1(74, xt); HF1(75, yt); HF1(76, ut); HF1(77,vt);
    HF2(78, xt, ut); HF2(79, yt, vt); HF2(80, xt, yt);
    HF1(91, p); HF1(92, q); HF1(93, path);
    event.nhS2s[i] = nh;
    event.chisqrS2s[i] = chisqr;
    event.path[i] = path;
    event.pS2s[i] = p;
    event.qS2s[i] = q;
    event.xtgtS2s[i] = xt;
    event.ytgtS2s[i] = yt;
    event.utgtS2s[i] = ut;
    event.vtgtS2s[i] = vt;
    event.thetaS2s[i] = theta;
    event.phiS2s[i] = phi;
    event.resP[i] = p - initial_momentum;
    for(Int_t i = 0; i<NumOfLayersVP; ++i){
      Int_t l = PlMinVP + i;
      Double_t vpx, vpy;
      Double_t vpu, vpv;
      track->GetTrajectoryLocalPosition(l, vpx, vpy);
      track->GetTrajectoryLocalDirection(l, vpu, vpv);
      event.vpx[i] = vpx;
      event.vpy[i] = vpy;
      event.vpu[i] = vpu;
      event.vpv[i] = vpv;
      HF2(100*l+1, vpx, vpu); HF2(100*l+2, vpy, vpv); HF2(100*l+3, vpx, vpy);
    }
    Double_t tof_seg = track->TofSeg()+1; // 1-origin
    if( tof_seg > 0 ){
      event.tofsegS2s[i] = tof_seg;
      TVector3 lposTof, lmomTof;
      track->TofLocalPos(lposTof);
      track->TofLocalMom(lmomTof);
      event.lxtofS2s[i] = lposTof.x();
      event.lytofS2s[i] = lposTof.y();
      event.lutofS2s[i] = lmomTof.x()/lmomTof.z();
      event.lvtofS2s[i] = lmomTof.y()/lmomTof.z();
      const auto& posTof = track->TofPos();
      const auto& momTof = track->TofMom();
      event.xtofS2s[i] = posTof.x();
      event.ytofS2s[i] = posTof.y();
      event.utofS2s[i] = momTof.x()/momTof.z();
      event.vtofS2s[i] = momTof.y()/momTof.z();
    }
    /*
#if UseTOF
    Double_t tof_seg = track->TofSeg();
    event.tofsegS2s[i] = tof_seg;
#else
    Double_t tof_x = track->GetLocalTrackOut()->GetX(zTOF);
    Double_t tof_seg = MathTools::Round(tof_x/75. + (NumOfSegTOF + 1)*0.5);
    event.tofsegS2s[i] = tof_seg;
# if 0
    std::cout << "posTof " << posTof << std::endl;
    std::cout << "momTof " << momTof << std::endl;
    std::cout << std::setw(10) << vecTof.X()
	      << std::setw(10) << vecTof.Y()
	      << std::setw(10) << sign*vecTof.Mod()
	      << std::setw(10) << TofSegS2s << std::endl;
# endif
    Double_t minres = 1.0e10;
#endif
    Double_t time = qnan;
    for(const auto& hit: hodoAna.GetHitContainer("TOF")){
   if(!hit) continue;
      Int_t seg = hit->SegmentId() + 1;
#if UseTOF
      if((Int_t)tof_seg == seg){
	time = hit->CMeanTime() - time0 + StofOffset;
      }
#else
      Double_t res = TMath::Abs(tof_seg - seg);
      if(res < minres){
	minres = res;
	time = hit->CMeanTime() - time0 + StofOffset;
      }
#endif
    }
    event.stof[i] = time;
    if(time > 0.){
      Double_t m2 = Kinematics::MassSquare(p, path, time);
      HF1(94, m2);
      event.m2[i] = m2;
# if 0
      std::ios::fmtflags pre_flags     = std::cout.flags();
      std::size_t        pre_precision = std::cout.precision();
      std::cout.setf(std::ios::fixed);
      std::cout.precision(5);
      std::cout << FUNC_NAME << std::endl
		<< "   Mom  = " << p     << std::endl
		<< "   Path = " << path << std::endl
		<< "   Time = " << time  << std::endl
		<< "   m2   = " << m2    << std::endl;
      std::cout.flags(pre_flags);
      std::cout.precision(pre_precision);
# endif
    }
    */

    for(Int_t j=0; j<nh; ++j){
      const auto& hit = track->GetHit(j);
      if(!hit) continue;
      Int_t layerId = hit->GetLayer();
      HF1(73, hit->GetLayer());
      Double_t wire = hit->GetHit()->GetWire();
      Double_t dt   = hit->GetHit()->GetDriftTime();
      Double_t dl   = hit->GetHit()->GetDriftLength();
      Double_t pos  = hit->GetCalLPos();
      Double_t res  = hit->GetResidual();
      const auto& lhit = hit->GetHit();
      Double_t xcal = lhit->GetXcal();
      Double_t ycal = lhit->GetYcal();
      HF1(100*layerId+11, Double_t(wire)-0.5);
      Double_t wp   = lhit->GetWirePosition();
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      HF1(100*layerId+14, pos);

      if(nh>17 && q>=0. && chisqr<200.){
      	HF1(100*layerId+15, res);
	HF2(100*layerId+16, pos, res);
      }
      HF2(100*layerId+17, xcal, ycal);
      if (std::abs(dl-std::abs(xcal-wp))<2.0){
	HF2(100*layerId+22, dt, std::abs(xcal-wp));
      }
      event.resG[layerId-1].push_back(res);
    }

    const auto& trSdcIn  = track->GetLocalTrackIn();
    const auto& trSdcOut = track->GetLocalTrackOut();
    if(trSdcIn){
      Int_t nhSdcIn=trSdcIn->GetNHit();
      Double_t x0in=trSdcIn->GetX0(), y0in=trSdcIn->GetY0();
      Double_t u0in=trSdcIn->GetU0(), v0in=trSdcIn->GetV0();
      Double_t chiin=trSdcIn->GetChiSquare();
      HF1(21, Double_t(nhSdcIn)); HF1(22, chiin);
      HF1(24, x0in); HF1(25, y0in); HF1(26, u0in); HF1(27, v0in);
      HF2(28, x0in, u0in); HF2(29, y0in, v0in); HF2(30, y0in, x0in);
      for(Int_t jin=0; jin<nhSdcIn; ++jin){
	Int_t layer=trSdcIn->GetHit(jin)->GetLayer();
	HF1(23, layer);
      }
    }
    if(trSdcOut){
      Int_t nhSdcOut=trSdcOut->GetNHit();
      Double_t x0out=trSdcOut->GetX(zTOF), y0out=trSdcOut->GetY(zTOF);
      Double_t u0out=trSdcOut->GetU0(), v0out=trSdcOut->GetV0();
      Double_t chiout=trSdcOut->GetChiSquare();
      HF1(51, Double_t(nhSdcOut)); HF1(52, chiout);
      HF1(54, x0out); HF1(55, y0out); HF1(56, u0out); HF1(57, v0out);
      ;;      HF2(58, x0out, u0out); HF2(59, y0out, v0out); HF2(60, y0out, x0out);
      for(Int_t jout=0; jout<nhSdcOut; ++jout){
	Int_t layer=trSdcOut->GetHit(jout)->GetLayer();
	HF1(53, layer);
      }
    }
  } // for in ntS2s

  for(Int_t i=0; i<ntS2s; ++i){
    const auto& track = DCAna.GetS2sTrack(i);
    if(!track) continue;
    const auto& trSdcIn =track->GetLocalTrackIn();
    const auto& trSdcOut=track->GetLocalTrackOut();
    if(!trSdcIn || !trSdcOut) continue;
    Double_t xin=trSdcIn->GetX0(), yin=trSdcIn->GetY0();
    Double_t uin=trSdcIn->GetU0(), vin=trSdcIn->GetV0();
    Double_t xout=trSdcOut->GetX0(), yout=trSdcOut->GetY0();
    Double_t uout=trSdcOut->GetU0(), vout=trSdcOut->GetV0();
    HF2(20011, xin, xout); HF2(20012, yin, yout);
    HF2(20013, uin, uout); HF2(20014, vin, vout);
  }

  if(ntS2s==0) return true;
  const auto& track = DCAna.GetS2sTrack(0);
  Double_t path = track->PathLengthToTOF();
  Double_t p    = track->PrimaryMomentum().Mag();
  //Double_t tTof[Event::nParticle];
  Double_t calt[Event::nParticle] = {
    Kinematics::CalcTimeOfFlight(p, path, pdg::PionMass()),
    Kinematics::CalcTimeOfFlight(p, path, pdg::KaonMass()),
    Kinematics::CalcTimeOfFlight(p, path, pdg::ProtonMass())
  };
  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = calt[i];
  }

#if 0
  // TOF
  {
    Int_t nh = hodoAna.GetNHits("TOF");
    for(Int_t i=0; i<nh; ++i){
      const auto& hit = hodoAna.GetHit("TOF", i);
      if(!hit) continue;
      Int_t seg=hit->SegmentId()+1;
      Double_t tu = hit->GetTUp(), td=hit->GetTDown();
      // Double_t ctu=hit->GetCTUp(), ctd=hit->GetCTDown();
      Double_t cmt=hit->CMeanTime();//, t= cmt-time0+StofOffset;//cmt-time0;
      Double_t ude=hit->GetAUp(), dde=hit->GetADown();
      // Double_t de=hit->DeltaE();
      // Double_t m2 = Kinematics::MassSquare(p, path, t);
      // event.tofmt[seg-1] = hit->MeanTime();
      event.utTofSeg[seg-1] = tu - time0 + StofOffset;
      event.dtTofSeg[seg-1] = td - time0 + StofOffset;
      // event.uctTofSeg[seg-1] = ctu - time0 + offset;
      // event.dctTofSeg[seg-1] = ctd - time0 + offset;
      event.udeTofSeg[seg-1] = ude;
      event.ddeTofSeg[seg-1] = dde;
      // event.ctTofSeg[seg-1]  = t;
      // event.deTofSeg[seg-1]  = de;
      // HF2(30000+100*seg+83, ude, calt[Event::Pion]+time0-StofOffset-cmt);
      // HF2(30000+100*seg+84, dde, calt[Event::Pion]+time0-StofOffset-cmt);
      // HF2(30000+100*seg+83, ude, calt[Event::Pion]+time0-StofOffset-tu);
      // HF2(30000+100*seg+84, dde, calt[Event::Pion]+time0-StofOffset-td);
    }

    const auto& cont = rawData.GetHodoRawHitContainer("TOF");
    Int_t NofHit = cont.size();
    for(Int_t i = 0; i<NofHit; ++i){
      auto hit = cont[i];
      Int_t seg = hit->SegmentId();
      event.tofua[seg] = hit->GetAdcUp();
      event.tofda[seg] = hit->GetAdcDown();
    }
  }
#endif

  HF1(1, 22.);



#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[nhittpc]: "<<std::setw(2)<<src.nhittpc<<" "<<std::endl;
#endif

  return true;
}

//_____________________________________________________________________
Bool_t
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  const std::size_t n = TFileCont.size();
  for(std::size_t i=0; i<n; ++i){
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  const Int_t    NbinSDC1DT =  240;
  const Double_t MinSDC1DT  = -30.;
  const Double_t MaxSDC1DT  = 170.;
  const Int_t    NbinSDC1DL =   55;
  const Double_t MinSDC1DL  = -0.5;
  const Double_t MaxSDC1DL  =  5.0;

  const Int_t    NbinSDC2DT =  360;
  const Double_t MinSDC2DT  = -50.;
  const Double_t MaxSDC2DT  = 250.;
  const Int_t    NbinSDC2DL =   85;
  const Double_t MinSDC2DL  = -0.5;
  const Double_t MaxSDC2DL  =   8.;

  const Int_t    NbinSDC3DT =  360;
  const Double_t MinSDC3DT  = -50.;
  const Double_t MaxSDC3DT  = 250.;
  const Int_t    NbinSDC3DL =   90;
  const Double_t MinSDC3DL  = -2.0;
  const Double_t MaxSDC3DL  =  7.0;

  const Int_t    NbinSDC4DT =  360;
  const Double_t MinSDC4DT  = -50.;
  const Double_t MaxSDC4DT  = 250.;
  const Int_t    NbinSDC4DL =   90;
  const Double_t MinSDC4DL  = -2.0;
  const Double_t MaxSDC4DL  =  7.0;

  const Int_t    NbinSDC5DT =  360;
  const Double_t MinSDC5DT  = -50.;
  const Double_t MaxSDC5DT  = 250.;
  const Int_t    NbinSDC5DL =   90;
  const Double_t MinSDC5DL  = -2.0;
  const Double_t MaxSDC5DL  =  7.0;

  const Int_t NbinSdcInRes   =  200;
  const Double_t MinSdcInRes = -1.;
  const Double_t MaxSdcInRes =  1.;

  const Int_t NbinSdcOutRes   = 1000;
  const Double_t MinSdcOutRes =  -5.;
  const Double_t MaxSdcOutRes =   5.;

  HB1(1, "Status", 30, 0., 30.);

  HB1(10, "#Tracks SdcIn", 10, 0., 10.);
  HB1(11, "#Hits of Track SdcIn", 15, 0., 15.);
  HB1(12, "Chisqr SdcIn", 500, 0., 50.);
  HB1(13, "LayerId SdcIn", 15, 0., 15.);
  HB1(14, "X0 SdcIn", 400, -100., 100.);
  HB1(15, "Y0 SdcIn", 400, -100., 100.);
  HB1(16, "U0 SdcIn", 200, -0.20, 0.20);
  HB1(17, "V0 SdcIn", 200, -0.20, 0.20);
  HB2(18, "U0%X0 SdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(19, "V0%Y0 SdcIn", 100, -100., 100., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 SdcIn", 100, -100., 100., 100, -100., 100.);

  HB1(21, "#Hits of Track SdcIn [S2sTrack]", 15, 0., 15.);
  HB1(22, "Chisqr SdcIn [S2sTrack]", 500, 0., 50.);
  HB1(23, "LayerId SdcIn [S2sTrack]", 15, 0., 15.);
  HB1(24, "X0 SdcIn [S2sTrack]", 400, -100., 100.);
  HB1(25, "Y0 SdcIn [S2sTrack]", 400, -100., 100.);
  HB1(26, "U0 SdcIn [S2sTrack]", 200, -0.20, 0.20);
  HB1(27, "V0 SdcIn [S2sTrack]", 200, -0.20, 0.20);
  HB2(28, "U0%X0 SdcIn [S2sTrack]", 100, -100., 100., 100, -0.20, 0.20);
  HB2(29, "V0%Y0 SdcIn [S2sTrack]", 100, -100., 100., 100, -0.20, 0.20);
  HB2(30, "X0%Y0 SdcIn [S2sTrack]", 100, -100., 100., 100, -100., 100.);

  HB1(40, "#Tracks SdcOut", 10, 0., 10.);
  HB1(41, "#Hits of Track SdcOut", 20, 0., 20.);
  HB1(42, "Chisqr SdcOut", 500, 0., 50.);
  HB1(43, "LayerId SdcOut", 20, 30., 50.);
  HB1(44, "X0 SdcOut", 1400, -1200., 1200.);
  HB1(45, "Y0 SdcOut", 1000, -500., 500.);
  HB1(46, "U0 SdcOut",  700, -0.35, 0.35);
  HB1(47, "V0 SdcOut",  200, -0.20, 0.20);
  HB2(48, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40);
  HB2(49, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20);
  HB2(50, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500., 500.);

  HB1(51, "#Hits of Track SdcOut [S2sTrack]", 20, 0., 20.);
  HB1(52, "Chisqr SdcOut [S2sTrack]", 500, 0., 50.);
  HB1(53, "LayerId SdcOut [S2sTrack]", 20, 30., 50.);
  HB1(54, "X0 SdcOut [S2sTrack]", 1400, -1200., 1200.);
  HB1(55, "Y0 SdcOut [S2sTrack]", 1000, -500., 500.);
  HB1(56, "U0 SdcOut [S2sTrack]",  700, -0.35, 0.35);
  HB1(57, "V0 SdcOut [S2sTrack]",  200, -0.10, 0.10);
  HB2(58, "U0%X0 SdcOut [S2sTrack]", 120, -600., 600., 100, -0.40, 0.40);
  HB2(59, "V0%Y0 SdcOut [S2sTrack]", 100, -500., 500., 100, -0.10, 0.10);
  HB2(60, "X0%Y0 SdcOut [S2sTrack]", 100, -700., 700., 100, -500., 500.);

  HB1(70, "#Tracks S2S", 10, 0., 10.);
  HB1(71, "#Hits of S2sTrack", 50, 0., 50.);
  HB1(72, "Chisqr S2sTrack", 500, 0., 50.);
  HB1(73, "LayerId S2sTrack", 90, 0., 90.);
  HB1(74, "Xtgt S2sTrack", 200, -100., 100.);
  HB1(75, "Ytgt S2sTrack", 200, -100., 100.);
  HB1(76, "Utgt S2sTrack", 300, -0.30, 0.30);
  HB1(77, "Vtgt S2sTrack", 300, -0.20, 0.20);
  HB2(78, "U%Xtgt S2sTrack", 100, -100., 100., 100, -0.25, 0.25);
  HB2(79, "V%Ytgt S2sTrack", 100, -100., 100., 100, -0.10, 0.10);
  HB2(80, "Y%Xtgt S2sTrack", 100, -100., 100., 100, -100., 100.);

  HB1(91, "P S2sTrack", 500, 0.00, 2.50);
  HB1(92, "q S2sTrack", 4, -2., 2.);
  HB1(93, "PathLength S2sTrack", 300, 7000., 10000.);
  HB1(94, "MassSqr", 600, -0.4, 1.4);

  // SdcInTracking
  for( Int_t i = 1; i <= NumOfLayersSdcIn; ++i ){
    TString tag;
    Int_t nwire = 0, nbindt = 1, nbindl = 1;
    Double_t mindt = 0., maxdt = 1., mindl = 0., maxdl = 1.;
    Int_t l = i + PlMinSdcIn - 1;
    if( i <= NumOfLayersSDC1 ){
      tag    = "SDC1";
      nwire  = MaxWireSDC1;
      nbindt = NbinSDC1DT;
      mindt  = MinSDC1DT;
      maxdt  = MaxSDC1DT;
      nbindl = NbinSDC1DL;
      mindl  = MinSDC1DL;
      maxdl  = MaxSDC1DL;
    }else{
      tag    = "SDC2";
      nwire  = MaxWireSDC2;
      nbindt = NbinSDC2DT;
      mindt  = MinSDC2DT;
      maxdt  = MaxSDC2DT;
      nbindl = NbinSDC2DL;
      mindl  = MinSDC2DL;
      maxdl  = MaxSDC2DL;
    }
    TString title1 = Form("HitPat SdcIn%2d", i);
    TString title2 = Form("DriftTime SdcIn%d", i);
    TString title3 = Form("DriftLength SdcIn%2d", i);
    TString title4 = Form("Position SdcIn%2d", i);
    TString title5 = Form("Residual SdcIn%2d", i);
    TString title6 = Form("Resid%%Pos SdcIn%2d", i);
    TString title7 = Form("Y%%Xcal SdcIn%2d", i);
    HB1(100*l+1, title1, nwire, 0., Double_t(nwire));
    HB1(100*l+2, title2, nbindt, mindt, maxdt);
    HB1(100*l+3, title3, nbindl, mindl, maxdl);
    HB1(100*l+4, title4, 500, -100., 100.);
    HB1(100*l+5, title5, NbinSdcInRes, MinSdcInRes, MaxSdcInRes);
    HB2(100*l+6, title6, 100, -100., 100., NbinSdcInRes, MinSdcInRes, MaxSdcInRes);
    HB2(100*l+7, title7, 100, -100., 100., 100, -50., 50.);
    title1 += " [S2sTrack]";
    title2 += " [S2sTrack]";
    title3 += " [S2sTrack]";
    title4 += " [S2sTrack]";
    title5 += " [S2sTrack]";
    title6 += " [S2sTrack]";
    title7 += " [S2sTrack]";
    HB1(100*l+11, title1, nwire, 0., Double_t(nwire));
    HB1(100*l+12, title2, nbindt, mindt, maxdt);
    HB1(100*l+13, title3, nbindl, mindl, maxdl);
    HB1(100*l+14, title4, 500, -100., 100.);
    HB1(100*l+15, title5, NbinSdcInRes, MinSdcInRes, MaxSdcInRes);
    HB2(100*l+16, title6, 100, -100., 100., NbinSdcInRes, MinSdcInRes, MaxSdcInRes);
    HB2(100*l+17, title7, 100, -100., 100., 100, -50., 50.);
    TString title22 = Form("DriftLength%%DriftTime SdcIn%2d [S2sTrack]", i);
    HB2(100*l+22, title22, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
  }

  // SdcOutTracking
  for( Int_t i = 1; i <= NumOfLayersSdcOut; ++i ){
    TString tag;
    Int_t nwire = 0, nbindt = 1, nbindl = 1;
    Double_t mindt = 0., maxdt = 1., mindl = 0., maxdl = 1.;
    Int_t l = i + PlMinSdcOut - 1;
    if( i <= NumOfLayersSDC3 ){
      tag    = "SDC3";
      nwire  = MaxWireSDC3;
      nbindt = NbinSDC3DT;
      mindt  = MinSDC3DT;
      maxdt  = MaxSDC3DT;
      nbindl = NbinSDC3DL;
      mindl  = MinSDC3DL;
      maxdl  = MaxSDC3DL;
    }else if( i <= NumOfLayersSDC3 + NumOfLayersSDC4 ){
      tag    = "SDC4";
      nwire  = MaxWireSDC4;
      nbindt = NbinSDC4DT;
      mindt  = MinSDC4DT;
      maxdt  = MaxSDC4DT;
      nbindl = NbinSDC4DL;
      mindl  = MinSDC4DL;
      maxdl  = MaxSDC4DL;
    }else{
      tag    = "SDC5";
      nwire  = (i==11 || i==12) ? MaxWireSDC5X : MaxWireSDC5Y;
      nbindt = NbinSDC5DT;
      mindt  = MinSDC5DT;
      maxdt  = MaxSDC5DT;
      nbindl = NbinSDC5DL;
      mindl  = MinSDC5DL;
      maxdl  = MaxSDC5DL;
    }
    TString title1 = Form("HitPat SdcOut%2d", i);
    TString title2 = Form("DriftTime SdcOut%d", i);
    TString title3 = Form("DriftLength SdcOut%2d", i);
    TString title4 = Form("Position SdcOut%2d", i);
    TString title5 = Form("Residual SdcOut%2d", i);
    TString title6 = Form("Resid%%Pos SdcOut%2d", i);
    TString title7 = Form("Y%%Xcal SdcOut%2d", i);
    HB1(100*l+1, title1, nwire, 0., Double_t(nwire));
    HB1(100*l+2, title2, nbindt, mindt, maxdt);
    HB1(100*l+3, title3, nbindl, mindl, maxdl);
    HB1(100*l+4, title4, 2000, -1000., 1000.);
    HB1(100*l+5, title5, NbinSdcOutRes, MinSdcOutRes, MaxSdcOutRes);
    HB2(100*l+6, title6, 2000, -1000., 1000., NbinSdcOutRes, MinSdcOutRes, MaxSdcOutRes);
    HB2(100*l+7, title7, 100, -600., 600., 100, -600., 600.);
    title1 += " [S2sTrack]";
    title2 += " [S2sTrack]";
    title3 += " [S2sTrack]";
    title4 += " [S2sTrack]";
    title5 += " [S2sTrack]";
    title6 += " [S2sTrack]";
    title7 += " [S2sTrack]";
    HB1(100*l+11, title1, nwire, 0., Double_t(nwire));
    HB1(100*l+12, title2, nbindt, mindt, maxdt);
    HB1(100*l+13, title3, nbindl, mindl, maxdl);
    HB1(100*l+14, title4, 2000, -1000., 1000.);
    HB1(100*l+15, title5, NbinSdcOutRes, MinSdcOutRes, MaxSdcOutRes);
    HB2(100*l+16, title6, 2000, -1000., 1000., NbinSdcOutRes, MinSdcOutRes, MaxSdcOutRes);
    HB2(100*l+17, title7, 100, -600., 600., 100, -600., 600.);
    TString title22 = Form("DriftLength%%DriftTime SdcOut%2d [S2sTrack]", i);
    HB2(100*l+22, title22, nbindt, mindt, maxdt, nbindl, mindl, maxdl);
  }
  /////////////////////

  // TOF in SdcOut/S2sTracking
  for( Int_t i = 1; i <= NumOfLayersTOF; ++i ){
  const Int_t    NbinSDC4DT =  360;
  const Double_t MinSDC4DT  = -50.;
  const Double_t MaxSDC4DT  = 250.;
  const Int_t    NbinSDC4DL =   90;
  const Double_t MinSDC4DL  = -2.0;
  const Double_t MaxSDC4DL  =  7.0;
    Int_t l = i + PlMinTOF - 1;
    TString title1 = Form("HitPat Tof%d", i);
    TString title2 = Form("DriftTime Tof%d", i);
    TString title3 = Form("DriftLength Tof%2d", i);
    TString title4 = Form("Position Tof%d", i);
    TString title5 = Form("Residual Tof%d", i);
    TString title6 = Form("Resid%%Pos Tof%d", i);
    TString title7 = Form("Y%%Xcal Tof%d", i);
    HB1(100*l+1, title1, 200, 0., 200.);
    HB1(100*l+2, title2, 400, 0., 200.);
    HB1(100*l+3, title3, 100, -2., 8.);
    HB1(100*l+4, title4, 1000, -1000., 1000.);
    HB1(100*l+5, title5, 200, -20., 20.);
    HB2(100*l+6, title6, 100, -1000., 1000., 100, -200., 200.);
    HB2(100*l+7, title6, 100, -1000., 1000., 100, -1000., 1000.);
    title1 += " [S2sTrack]";
    title2 += " [S2sTrack]";
    title3 += " [S2sTrack]";
    title4 += " [S2sTrack]";
    title5 += " [S2sTrack]";
    title6 += " [S2sTrack]";
    title7 += " [S2sTrack]";
    HB1(100*l+11, title1, 200, 0., 200.);
    HB1(100*l+12, title2, 400, 0., 200.);
    HB1(100*l+13, title3, 100, -2., 8.);
    HB1(100*l+14, title4, 1000, -1000., 1000.);
    HB1(100*l+15, title5, 200, -20., 20.);
    HB2(100*l+16, title6, 100, -1000., 1000., 100, -200., 200.);
    HB2(100*l+17, title7, 100, -1000., 1000., 100, -1000., 1000.);
  }

  for( Int_t i = 1; i <= NumOfLayersVP; ++i ){
    Int_t l = i + PlMinVP - 1;
    TString title1 = Form("U%%X VP%d S2sTrack", i);
    TString title2 = Form("V%%Y VP%d S2sTrack", i);
    TString title3 = Form("Y%%X VP%d S2sTrack", i);
    HB2(100*l+1, title1, 500, -500., 500., 1000, -0.5, 0.5);
    HB2(100*l+2, title2, 300, -300., 300., 200, -0.1, 0.1);
    HB2(100*l+3, title3, 500, -500., 500., 300, -300., 300.);
  }

  HB2(20001, "Xout%Xin", 100, -200., 200., 100, -200., 200.);
  HB2(20002, "Yout%Yin", 100, -200., 200., 100, -200., 200.);
  HB2(20003, "Uout%Uin", 100, -0.5,  0.5,  100, -0.5,  0.5);
  HB2(20004, "Vin%Vout", 100, -0.1,  0.1,  100, -0.1,  0.1);

  HB2(20011, "Xout%Xin [S2sTrack]", 100, -200., 200., 100, -200., 200.);
  HB2(20012, "Yout%Yin [S2sTrack]", 100, -200., 200., 100, -200., 200.);
  HB2(20013, "Uout%Uin [S2sTrack]", 100, -0.5,  0.5,  100, -0.5,  0.5);
  HB2(20014, "Vin%Vout [S2sTrack]", 100, -0.1,  0.1,  100, -0.1,  0.1);

  HBTree("s2s", "tree of DstS2sTracking using Geant4 input");
  // Geant4
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("x0",     &event.x0,     "x0/D");
  tree->Branch("y0",     &event.y0,     "y0/D");
  tree->Branch("z0",     &event.z0,     "z0/D");
  tree->Branch("xp0",    &event.xp0,    "xp0/D");
  tree->Branch("yp0",    &event.yp0,    "yp0/D");
  tree->Branch("phi0",   &event.phi0,   "phi0/D");
  tree->Branch("theta0", &event.theta0, "theta0/D");
  tree->Branch("p0",     &event.p0,     "p0/D");
  tree->Branch("pB",     &event.pB,     "pB/D");

  tree->Branch("nhTof",  &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",  event.TofSeg,  "TofSeg[nhTof]/D");

  //Tracking
  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("nlSdcIn",    &event.nlSdcIn,     "nlSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("wposSdcIn",  event.wposSdcIn,  Form("wposSdcIn[%d]/D", NumOfLayersSdcIn));
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  tree->Branch("ntSdcOut",   &event.ntSdcOut,     "ntSdcOut/I");
  tree->Branch("nlSdcOut",   &event.nlSdcOut,     "nlSdcOut/I");
  tree->Branch("nhSdcOut",    event.nhSdcOut,     "nhSdcOut[ntSdcOut]/I");
  tree->Branch("wposSdcOut", event.wposSdcOut, Form("wposSdcOut[%d]/D", NumOfLayersSdcOut));
  tree->Branch("chisqrSdcOut",event.chisqrSdcOut, "chisqrSdcOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut",    event.x0SdcOut,     "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut",    event.y0SdcOut,     "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut",    event.u0SdcOut,     "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",    event.v0SdcOut,     "v0SdcOut[ntSdcOut]/D");

  // S2S Tracking
  tree->Branch("ntS2s",    &event.ntS2s,     "ntS2s/I");
  tree->Branch("nlS2s",    &event.nlS2s,     "nlS2s/I");
  tree->Branch("nhS2s",     event.nhS2s,     "nhS2s[ntS2s]/I");
  tree->Branch("chisqrS2s", event.chisqrS2s, "chisqrS2s[ntS2s]/D");
  tree->Branch("stof",         event.stof,         "stof[ntS2s]/D");
  tree->Branch("path",         event.path,         "path[ntS2s]/D");
  tree->Branch("pS2s",      event.pS2s,      "pS2s[ntS2s]/D");
  tree->Branch("qS2s",      event.qS2s,      "qS2s[ntS2s]/D");
  tree->Branch("m2",           event.m2,           "m2[ntS2s]/D");

  tree->Branch("xtgtS2s",   event.xtgtS2s,   "xtgtS2s[ntS2s]/D");
  tree->Branch("ytgtS2s",   event.ytgtS2s,   "ytgtS2s[ntS2s]/D");
  tree->Branch("utgtS2s",   event.utgtS2s,   "utgtS2s[ntS2s]/D");
  tree->Branch("vtgtS2s",   event.vtgtS2s,   "vtgtS2s[ntS2s]/D");

  tree->Branch("thetaS2s",  event.thetaS2s,  "thetaS2s[ntS2s]/D");
  tree->Branch("phiS2s",    event.phiS2s,    "phiS2s[ntS2s]/D");
  tree->Branch("resP",    event.resP,   "resP[ntS2s]/D");

  tree->Branch("xtofS2s",   event.xtofS2s,   "xtofS2s[ntS2s]/D");
  tree->Branch("ytofS2s",   event.ytofS2s,   "ytofS2s[ntS2s]/D");
  tree->Branch("utofS2s",   event.utofS2s,   "utofS2s[ntS2s]/D");
  tree->Branch("vtofS2s",   event.vtofS2s,   "vtofS2s[ntS2s]/D");

  tree->Branch("lxtofS2s",   event.lxtofS2s,   "lxtofS2s[ntS2s]/D");
  tree->Branch("lytofS2s",   event.lytofS2s,   "lytofS2s[ntS2s]/D");
  tree->Branch("lutofS2s",   event.lutofS2s,   "lutofS2s[ntS2s]/D");
  tree->Branch("lvtofS2s",   event.lvtofS2s,   "lvtofS2s[ntS2s]/D");
  tree->Branch("tofsegS2s", event.tofsegS2s, "tofsegS2s[ntS2s]/D");

  tree->Branch("vpx",          event.vpx,          Form("vpx[%d]/D", NumOfLayersVP));
  tree->Branch("vpy",          event.vpy,          Form("vpy[%d]/D", NumOfLayersVP));
  tree->Branch("vpu",          event.vpu,          Form("vpu[%d]/D", NumOfLayersVP));
  tree->Branch("vpv",          event.vpv,          Form("vpv[%d]/D", NumOfLayersVP));

  event.resL.resize(PlMaxTOF);
  for( Int_t i = PlMinSdcIn;  i<= PlMaxSdcIn;  i++ ) tree->Branch(Form("ResL%d", i), &event.resL[i-1]);
  for( Int_t i = PlMinSdcOut; i<= PlMaxSdcOut; i++ ) tree->Branch(Form("ResL%d", i), &event.resL[i-1]);
  for( Int_t i = PlMinTOF;    i<= PlMaxTOF;    i++ ) tree->Branch(Form("ResL%d", i), &event.resL[i-1]);

  event.resG.resize(PlMaxTOF);
  for( Int_t i = PlMinSdcIn;  i<= PlMaxSdcIn;  i++ ) tree->Branch(Form("ResG%d", i), &event.resG[i-1]);
  for( Int_t i = PlMinSdcOut; i<= PlMaxSdcOut; i++ ) tree->Branch(Form("ResG%d", i), &event.resG[i-1]);
  for( Int_t i = PlMinTOF;    i<= PlMaxTOF;    i++ ) tree->Branch(Form("ResG%d", i), &event.resG[i-1]);

  ////////// Bring Address From Geant4 root file
  TTreeReaderCont[kS2sGeant4] = new TTreeReader("g4s2s", TFileCont[kS2sGeant4]);
  const auto& reader = TTreeReaderCont[kS2sGeant4];
  src.evnum  = new TTreeReaderValue<Int_t>(*reader,     "evnum");
  src.x0     = new TTreeReaderValue<Double_t>(*reader,  "x0");
  src.y0     = new TTreeReaderValue<Double_t>(*reader,  "y0");
  src.z0     = new TTreeReaderValue<Double_t>(*reader,  "z0");
  src.xp0    = new TTreeReaderValue<Double_t>(*reader,  "xp0");
  src.yp0    = new TTreeReaderValue<Double_t>(*reader,  "yp0");
  src.phi0   = new TTreeReaderValue<Double_t>(*reader,  "phi0");
  src.theta0 = new TTreeReaderValue<Double_t>(*reader,  "theta0");
  src.p0     = new TTreeReaderValue<Double_t>(*reader,  "p0");
  src.pB     = new TTreeReaderValue<Double_t>(*reader,  "pB");
  src.PRM    = new TTreeReaderArray<TParticle>(*reader, "PRM");
  src.SDC1   = new TTreeReaderArray<TParticle>(*reader, "SDC1");
  src.SDC2   = new TTreeReaderArray<TParticle>(*reader, "SDC2");
  src.SDC3   = new TTreeReaderArray<TParticle>(*reader, "SDC3");
  src.SDC4   = new TTreeReaderArray<TParticle>(*reader, "SDC4");
  src.SDC5   = new TTreeReaderArray<TParticle>(*reader, "SDC5");
  src.TOF    = new TTreeReaderArray<TParticle>(*reader, "TOF");
  src.AC1    = new TTreeReaderArray<TParticle>(*reader, "AC1");
  src.WC     = new TTreeReaderArray<TParticle>(*reader, "WC");
  src.VP     = new TTreeReaderArray<TParticle>(*reader, "VP");

  return true;
}

//_____________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<FieldMan>("FLDMAP")   &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
