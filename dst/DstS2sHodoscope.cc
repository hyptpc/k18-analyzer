// -*- C++ -*-

#include "DstHelper.hh"

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <TMath.h>

#include <filesystem_util.hh>
#include <UnpackerManager.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "Kinematics.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "HodoPHCMan.hh"

namespace
{
using namespace root;
using namespace dst;
using hddaq::unpacker::GUnpacker;
const auto qnan = TMath::QuietNaN();
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
const auto& gCounter = debug::ObjectCounter::GetInstance();
const Bool_t USE_M2 = false;
const Bool_t USE_XYCut = false;
TString ClassName() { return TString("DstS2sHodoscope"); }
}

namespace dst
{
enum kArgc
{
  kProcess, kConfFile,
  kS2sTracking, kHodoscope, kOutFile, nArgc
};
std::vector<TString> ArgName =
{ "[Process]", "[ConfFile]", "[S2sTracking]",
  "[Hodoscope]", "[OutFile]" };
std::vector<TString> TreeName =
{ "", "", "s2s", "hodo", "" };
std::vector<TFile*> TFileCont;
std::vector<TTree*> TTreeCont;
std::vector<TTreeReader*> TTreeReaderCont;
}

//_____________________________________________________________________
struct Event
{
  Int_t status;

  // SdcOut
  Int_t ntSdcOut;
  Double_t chisqrSdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];

  // S2sTracking
  Int_t ntS2s;
  Double_t path[MaxHits];
  Double_t pS2s[MaxHits];
  Double_t qS2s[MaxHits];
  Double_t chisqrS2s[MaxHits];
  Double_t thetaS2s[MaxHits];
  Double_t xtgtS2s[MaxHits];
  Double_t ytgtS2s[MaxHits];
  Double_t utgtS2s[MaxHits];
  Double_t vtgtS2s[MaxHits];
  Double_t xtofS2s[MaxHits];
  Double_t ytofS2s[MaxHits];
  Double_t utofS2s[MaxHits];
  Double_t vtofS2s[MaxHits];
  Double_t lxtofS2s[MaxHits];
  Double_t lytofS2s[MaxHits];
  Double_t lutofS2s[MaxHits];
  Double_t lvtofS2s[MaxHits];
  Double_t tofsegS2s[MaxHits];
  Double_t vpx[NumOfLayersVP*MaxHits];
  Double_t vpy[NumOfLayersVP*MaxHits];
  Double_t vpu[NumOfLayersVP*MaxHits];
  Double_t vpv[NumOfLayersVP*MaxHits];

  // Hodoscope
  Int_t trigflag[NumOfSegTrig];
  Int_t trigpat[NumOfSegTrig];

  Int_t nhBh1;
  Int_t csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];

  Int_t nhBh2;
  Int_t csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t btof[NumOfSegBH1*MaxDepth];
  Double_t Time0Seg;

  Int_t nhTof;
  Int_t csTof[NumOfSegTOF];
  Double_t TofSeg[NumOfSegTOF];
  Double_t tTof[NumOfSegTOF];
  Double_t dtTof[NumOfSegTOF];
  Double_t deTof[NumOfSegTOF];

  Int_t nhAc1;
  Int_t csAc1[NumOfSegAC1];
  Double_t Ac1Seg[NumOfSegAC1];
  Double_t tAc1[NumOfSegAC1];

  // S2sHodoscope
  Int_t    m2Combi;
  Double_t beta[MaxHits];
  Double_t stof[MaxHits];
  Double_t cstof[MaxHits];
  Double_t m2[MaxHits];

  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];
  enum eParticle { Pion, Kaon, Proton, nParticle };
  Double_t tTofCalc[nParticle];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

};

//_____________________________________________________________________
struct Src
{
  Int_t trigflag[NumOfSegTrig];
  Int_t trigpat[NumOfSegTrig];

  Int_t ntSdcOut;
  Double_t chisqrSdcOut[MaxHits];
  Double_t x0SdcOut[MaxHits];
  Double_t y0SdcOut[MaxHits];
  Double_t u0SdcOut[MaxHits];
  Double_t v0SdcOut[MaxHits];

  Int_t    ntS2s;
  Double_t path[MaxHits];
  Double_t pS2s[MaxHits];
  Double_t qS2s[MaxHits];
  Double_t chisqrS2s[MaxHits];
  Double_t thetaS2s[MaxHits];
  Double_t xtgtS2s[MaxHits];
  Double_t ytgtS2s[MaxHits];
  Double_t utgtS2s[MaxHits];
  Double_t vtgtS2s[MaxHits];
  Double_t xtofS2s[MaxHits];
  Double_t ytofS2s[MaxHits];
  Double_t utofS2s[MaxHits];
  Double_t vtofS2s[MaxHits];
  Double_t lxtofS2s[MaxHits];
  Double_t lytofS2s[MaxHits];
  Double_t lutofS2s[MaxHits];
  Double_t lvtofS2s[MaxHits];
  Double_t tofsegS2s[MaxHits];
  Double_t vpx[NumOfLayersVP];
  Double_t vpy[NumOfLayersVP];
  Double_t vpu[NumOfLayersVP];
  Double_t vpv[NumOfLayersVP];

  Int_t    nhBh1;
  Int_t    csBh1[NumOfSegBH1*MaxDepth];
  Double_t Bh1Seg[NumOfSegBH1*MaxDepth];
  Double_t tBh1[NumOfSegBH1*MaxDepth];
  Double_t dtBh1[NumOfSegBH1*MaxDepth];
  Double_t deBh1[NumOfSegBH1*MaxDepth];

  Int_t    nhBh2;
  Int_t    csBh2[NumOfSegBH2*MaxDepth];
  Double_t Bh2Seg[NumOfSegBH2*MaxDepth];
  Double_t tBh2[NumOfSegBH2*MaxDepth];
  Double_t dtBh2[NumOfSegBH2*MaxDepth];
  Double_t t0Bh2[NumOfSegBH2*MaxDepth];
  Double_t deBh2[NumOfSegBH2*MaxDepth];

  Double_t Time0Seg;
  Double_t deTime0;
  Double_t Time0;
  Double_t CTime0;

  Int_t    nhTof;
  Int_t    csTof[NumOfSegTOF];
  Double_t TofSeg[NumOfSegTOF];
  Double_t tTof[NumOfSegTOF];
  Double_t dtTof[NumOfSegTOF];
  Double_t deTof[NumOfSegTOF];

  Int_t nhAc1;
  Int_t csAc1[NumOfSegAC1];
  Double_t Ac1Seg[NumOfSegAC1];
  Double_t tAc1[NumOfSegAC1];

  ////////// for HodoParam
  Double_t tofua[NumOfSegTOF];
  Double_t tofda[NumOfSegTOF];
  Double_t utTofSeg[NumOfSegTOF][MaxDepth];
  Double_t dtTofSeg[NumOfSegTOF][MaxDepth];
  Double_t udeTofSeg[NumOfSegTOF];
  Double_t ddeTofSeg[NumOfSegTOF];

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
bool
dst::InitializeEvent()
{
  event.status   = 0;
  event.ntSdcOut = 0;
  event.ntS2s = 0;
  event.nhBh1    = 0;
  event.nhBh2    = 0;
  event.nhTof    = 0;
  event.nhAc1    = 0;
  event.m2Combi  = 0;
  event.Time0Seg = qnan;

  for(Int_t i=0; i<NumOfSegTrig; ++i){
    event.trigflag[i] = -1;
    event.trigpat[i]  = -1;
  }

  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = qnan;
  }

  for(Int_t i=0;i<MaxHits;++i){
    event.chisqrSdcOut[i] = qnan;
    event.x0SdcOut[i]     = qnan;
    event.y0SdcOut[i]     = qnan;
    event.u0SdcOut[i]     = qnan;
    event.v0SdcOut[i]     = qnan;
    event.path[i]         = qnan;
    event.pS2s[i]      = qnan;
    event.qS2s[i]      = qnan;
    event.chisqrS2s[i] = qnan;
    event.thetaS2s[i] = qnan;
    event.xtgtS2s[i]  = qnan;
    event.ytgtS2s[i]  = qnan;
    event.utgtS2s[i]  = qnan;
    event.vtgtS2s[i]  = qnan;
    event.xtofS2s[i]  = qnan;
    event.ytofS2s[i]  = qnan;
    event.utofS2s[i]  = qnan;
    event.vtofS2s[i]  = qnan;
    event.lxtofS2s[i]  = qnan;
    event.lytofS2s[i]  = qnan;
    event.lutofS2s[i]  = qnan;
    event.lvtofS2s[i]  = qnan;
    event.tofsegS2s[i]  = qnan;
  }

  for (Int_t l = 0; l < NumOfLayersVP; ++l) {
    event.vpx[l] = qnan;
    event.vpy[l] = qnan;
    event.vpu[l] = qnan;
    event.vpv[l] = qnan;
  }

  for(Int_t i=0;i<NumOfSegBH1*MaxDepth;++i){
    event.Bh1Seg[i] = -1;
    event.csBh1[i]  = 0;
    event.tBh1[i]   = qnan;
    event.dtBh1[i]  = qnan;
    event.deBh1[i]  = qnan;
    event.btof[i]   = qnan;
  }

  for(Int_t i=0;i<NumOfSegBH2*MaxDepth;++i){
    event.Bh2Seg[i] = -1;
    event.csBh2[i]  = 0;
    event.tBh2[i]   = qnan;
    event.dtBh2[i]  = qnan;
    event.t0Bh2[i]  = qnan;
    event.deBh2[i]  = qnan;
  }

  for(Int_t i=0;i<NumOfSegTOF;++i){
    event.TofSeg[i] = -1;
    event.csTof[i]  = 0;
    event.tTof[i]   = qnan;
    event.dtTof[i]  = qnan;
    event.deTof[i]  = qnan;
    event.tofua[i] = qnan;
    event.tofda[i] = qnan;
    for(Int_t m=0; m<MaxDepth; ++m){
      event.utTofSeg[i][m] = qnan;
      event.dtTofSeg[i][m] = qnan;
    }
    event.udeTofSeg[i] = qnan;
    event.ddeTofSeg[i] = qnan;
  }

  for(Int_t i=0;i<NumOfSegAC1;++i){
    event.Ac1Seg[i] = -1;
    event.csAc1[i]  = 0;
    event.tAc1[i]   = qnan;
  }

  for(Int_t i=0; i<Event::nParticle; ++i){
    event.tTofCalc[i] = qnan;
  }

  for(Int_t i=0; i<MaxHits; ++i){
    event.beta[i]  = qnan;
    event.stof[i]  = qnan;
    event.cstof[i] = qnan;
    event.m2[i]    = qnan;
  }
  return true;
}

//_____________________________________________________________________
bool
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
bool
dst::DstRead(Int_t ievent)
{
  static const auto StofOffset = gUser.GetParameter("StofOffset");

  if(ievent%10000==0){
    std::cout << "#D Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1(1, event.status++);

  event.ntSdcOut = src.ntSdcOut;
  event.ntS2s = src.ntS2s;
  event.nhBh1    = src.nhBh1;
  event.nhBh2    = src.nhBh2;
  event.nhTof    = src.nhTof;
  event.nhAc1    = src.nhAc1;

#if 0
  std::cout<<"[event]: "<<std::setw(6)<<ievent<<" ";
  std::cout<<"[ntS2s]: "<<std::setw(2)<<src.ntS2s<<" ";
  std::cout<<"[nhBh1]: "<<std::setw(2)<<src.nhBh1<<" ";
  std::cout<<"[nhBh2]: "<<std::setw(2)<<src.nhBh2<<" ";
  std::cout<<"[nhTof]: "<<std::setw(2)<<src.nhTof<<" "<<std::endl;
  std::cout<<"[nhAc1]: "<<std::setw(2)<<src.nhAc1<<" "<<std::endl;
#endif

  for(Int_t i=0;i<NumOfSegTrig;++i){
    Int_t tdc = src.trigflag[i];
    if(tdc<=0) continue;
    event.trigpat[i]  = i + 1;
    event.trigflag[i] = tdc;
  }

  // if(event.nhBh1<=0) return true;
  HF1(1, event.status++);

  // if(event.nhBh2<=0) return true;
  HF1(1, event.status++);

  // if(event.nhTof<=0) return true;
  HF1(1, event.status++);

  // if(event.ntS2s<=0) return true;
  // if(event.ntS2s>MaxHits)
  //   event.ntS2s = MaxHits;

  HF1(1, event.status++);

  Double_t time0 = src.CTime0;
  for(Int_t i=0; i<src.nhBh2; ++i){
    event.Bh2Seg[i] = src.Bh2Seg[i];
    event.tBh2[i]   = src.tBh2[i];
    event.t0Bh2[i]  = src.t0Bh2[i];
    event.deBh2[i]  = src.deBh2[i];
  }
  event.Time0Seg = src.Time0Seg;

  ////////// for BeamTof
  // Double_t btof = qnan;
  for(Int_t i=0; i<src.nhBh1; ++i){
    event.Bh1Seg[i] = src.Bh1Seg[i];
    event.tBh1[i]   = src.tBh1[i];
    event.deBh1[i]  = src.deBh1[i];
    event.btof[i]   = src.tBh1[i] - time0;
    // if(i==0) btof = src.tBh1[i] - time0;
  }

  for(Int_t it=0; it<src.ntSdcOut; ++it){
    event.chisqrSdcOut[it] = src.chisqrSdcOut[it];
    event.x0SdcOut[it] = src.x0SdcOut[it];
    event.y0SdcOut[it] = src.y0SdcOut[it];
    event.u0SdcOut[it] = src.u0SdcOut[it];
    event.v0SdcOut[it] = src.v0SdcOut[it];
  }

  for(Int_t it=0; it<src.nhTof; ++it){
    event.csTof[it]  = src.csTof[it];
    event.TofSeg[it] = src.TofSeg[it];
    event.tTof[it]   = src.tTof[it];
    event.dtTof[it]  = src.dtTof[it];
    event.deTof[it]  = src.deTof[it];
  }
  for(Int_t it=0; it<NumOfSegTOF; ++it){
    event.tofua[it] = src.tofua[it];
    event.tofda[it] = src.tofda[it];
    for(Int_t m=0; m<MaxDepth; ++m){
      event.utTofSeg[it][m] = src.utTofSeg[it][m];
      event.dtTofSeg[it][m] = src.dtTofSeg[it][m];
    }
    event.udeTofSeg[it] = src.udeTofSeg[it];
    event.ddeTofSeg[it] = src.ddeTofSeg[it];
    // TOF ADC-Pedestal
    Bool_t is_hit = false;
    for(Int_t i=0; i<event.nhTof; ++i){
      if(it == event.TofSeg[i]-1) is_hit = true;
    }
    if(!is_hit){
      HF1(30000+(it+1)*100+1,  event.tofua[it]);
      HF1(30000+(it+1)*100+11, event.tofda[it]);
    }
  }

  for(Int_t it=0; it<src.nhAc1; ++it){
    event.csAc1[it]  = src.csAc1[it];
    event.Ac1Seg[it] = src.Ac1Seg[it];
    event.tAc1[it]   = src.tAc1[it];
  }

  Int_t m2Combi = event.nhTof*event.ntS2s;
  if(m2Combi>MaxHits || m2Combi<0){
    std::cout << FUNC_NAME << " too much m2Combi : " << m2Combi << std::endl;
    return false;
  }
  event.m2Combi = m2Combi;

  HF1(1, event.status++);
  Int_t mm=0;
  for(Int_t it=0; it<src.ntS2s; ++it){
    event.path[it] = src.path[it];
    event.pS2s[it] = src.pS2s[it];
    event.qS2s[it] = src.qS2s[it];
    event.chisqrS2s[it] = src.chisqrS2s[it];
    event.thetaS2s[it]  = src.thetaS2s[it];
    event.xtgtS2s[it] = src.xtgtS2s[it];
    event.ytgtS2s[it] = src.ytgtS2s[it];
    event.utgtS2s[it] = src.utgtS2s[it];
    event.vtgtS2s[it] = src.vtgtS2s[it];
    event.xtofS2s[it] = src.xtofS2s[it];
    event.ytofS2s[it] = src.ytofS2s[it];
    event.utofS2s[it] = src.utofS2s[it];
    event.vtofS2s[it] = src.vtofS2s[it];
    event.lxtofS2s[it] = src.lxtofS2s[it];
    event.lytofS2s[it] = src.lytofS2s[it];
    event.lutofS2s[it] = src.lutofS2s[it];
    event.lvtofS2s[it] = src.lvtofS2s[it];
    event.tofsegS2s[it] = src.tofsegS2s[it];
    Double_t xtgt  = event.xtgtS2s[it];
    Double_t ytgt  = event.ytgtS2s[it];
    Double_t lytof = event.lytofS2s[it];
    Double_t pS2s = event.pS2s[it];
    Double_t qS2s = event.qS2s[it];
    if(event.chisqrS2s[it] < 200){
      HF1(10, pS2s);
      HF1(13, event.path[it]);
      HF1(50, event.chisqrS2s[it]);
      HF1(51, event.xtofS2s[it]);
      HF1(52, event.ytofS2s[it]);
      HF1(53, event.utofS2s[it]);
      HF1(54, event.vtofS2s[it]);
      HF1(55, event.lxtofS2s[it]);
      HF1(56, event.lytofS2s[it]);
      HF1(57, event.lutofS2s[it]);
      HF1(58, event.lvtofS2s[it]);
      HF1(59, event.tofsegS2s[it]);
    }
    if(src.ntS2s == 1){
      for (Int_t l = 0; l < NumOfLayersVP; ++l) {
	event.vpx[l] = src.vpx[l];
	event.vpy[l] = src.vpy[l];
	event.vpu[l] = src.vpu[l];
	event.vpv[l] = src.vpv[l];
      }
    }
    if(it == 0){
      event.tTofCalc[Event::Pion] =
        Kinematics::CalcTimeOfFlight(event.pS2s[it],
                                     event.path[it],
                                     pdg::PionMass());
      event.tTofCalc[Event::Kaon] =
        Kinematics::CalcTimeOfFlight(event.pS2s[it],
                                     event.path[it],
                                     pdg::KaonMass());
      event.tTofCalc[Event::Proton] =
        Kinematics::CalcTimeOfFlight(event.pS2s[it],
                                     event.path[it],
                                     pdg::ProtonMass());
    }

    for(Int_t itof=0; itof<src.nhTof; ++itof){
      Int_t tofseg = (Int_t)event.TofSeg[itof];
      if( tofseg != event.tofsegS2s[0] ) continue;

      ////////// TimeCut
      Double_t stof = event.tTof[itof] - time0 + StofOffset;
      Double_t cstof = stof;
      // gPHC.DoStofCorrection(8, 0, tofseg-1, 2, stof, btof, cstof);
      Double_t beta = event.path[it]/cstof/MathTools::C();
      event.beta[mm] = beta;
      event.stof[mm] = stof;// - event.tTofCalc[0];
      event.cstof[mm] = cstof;
      Double_t m2 = Kinematics::MassSquare(pS2s, event.path[it], cstof);
      event.m2[mm] = m2;
      // Bool_t is_pion = (// qS2s > 0
      //                   // &&
      //                   TMath::Abs(m2-pdg::PionMass()*pdg::PionMass()) < 0.1);
      // Bool_t is_kaon = (qS2s > 0
      //                   && TMath::Abs(m2-pdg::KaonMass()*pdg::KaonMass()) < 0.12);
#if 0
      std::cout << "#D DebugPrint() Event : " << ievent << std::endl
		<< std::setprecision(3) << std::fixed
		<< "   time0   : " << time0 << std::endl
		<< "   offset  : " << StofOffset << std::endl
		<< "   tTof    : " << event.tTof[itof] << std::endl
		<< "   stof    : " << stof << std::endl
		<< "   pS2s    : " << pS2s << std::endl
		<< "   m2      : " << m2 << std::endl;
#endif

      if(event.chisqrS2s[it] < 200.){
        for(Int_t ip=0; ip<Event::nParticle; ++ip){
          if(USE_M2){
            if(ip == Event::Pion && TMath::Abs(m2-0.0194) > 0.1) continue;
            if(ip == Event::Proton && TMath::Abs(m2-0.88) > 0.2) continue;
          }
          HF1(33000+tofseg*100+ip+1,  event.tofua[tofseg-1]);
          HF1(33000+tofseg*100+ip+11, event.tofda[tofseg-1]);
        }
        if(TMath::Abs(lytof) < 5.0){
          HF1(20000+tofseg*100+1, event.dtTof[itof]);
          HF2(20000+1, event.dtTof[itof], tofseg-1);
        }

	Bool_t xy_ok = USE_XYCut
	  ? (TMath::Abs(xtgt) < 25. && TMath::Abs(ytgt) < 20.)
	  : true;
	if(xy_ok){
	  if( it == 0 ){
	    for(Int_t ip=0; ip<Event::nParticle; ++ip){
	      if(USE_M2){
		if(ip == Event::Pion && m2 > 0.2) continue;
		if(ip == Event::Kaon) continue;
		if(ip == Event::Proton && (m2 < 0.5 || !xy_ok)) continue;
	      }
	      else{
		Bool_t flagAc1 = false;
		for( Int_t it=0, nhit=event.nhAc1; it<nhit; it++ ){
		  if( event.Ac1Seg[it] == 21 ){ flagAc1=true; break; }
		}
		if(ip == Event::Pion && !flagAc1) continue; // select Pion
		if(ip == Event::Kaon) continue;
		if(ip == Event::Proton && (m2 < 0.5 || !xy_ok)) continue;
	      }
	      HF2(10000+ip+1, cstof-event.tTofCalc[ip], tofseg-1);
	      HF1(10000+tofseg*100+ip+1, cstof-event.tTofCalc[ip]);
	      if(TMath::Abs(event.dtTof[itof] < 0.1)){
		for(Int_t mh=0; mh<MaxDepth; ++mh){
		  HF2(40000+tofseg*100+ip+1,
		      event.udeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		  HF2(40000+tofseg*100+ip+11,
		      event.ddeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		  HF2(40000+tofseg*100,
		      event.udeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		  HF2(40000+tofseg*100+10,
		      event.ddeTofSeg[tofseg-1],
		      event.tTofCalc[ip] - stof);
		}
	      }
	    }
	  }
	  HF1(11, qS2s*m2);
	  HF1(12, beta);
	  HF1(14, stof);
	  HF2(20, qS2s*m2, pS2s);
	  if(!TMath::IsNaN(event.Time0Seg)){
	    HF1(1100+(Int_t)event.Time0Seg, qS2s*m2);
	    HF2(1200+(Int_t)event.Time0Seg, qS2s*m2, pS2s);
	  }
	  HF1(2100+tofseg, qS2s*m2);
	  HF2(2200+tofseg, qS2s*m2, pS2s);
	}
      }
      ++mm;
    }
  }

  HF1(1, event.status++);

  return true;
}

//_____________________________________________________________________
bool
dst::DstClose()
{
  TFileCont[kOutFile]->Write();
  std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  TFileCont[kOutFile]->Close();

  for(Int_t i=0, n=TFileCont.size(); i<n; ++i){
    if(TTreeCont[i]) delete TTreeCont[i];
    if(TFileCont[i]) delete TFileCont[i];
  }
  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeHistograms()
{
  const Int_t NBinP = 220;
  const Double_t MinP =   0;
  const Double_t MaxP = 2.2;
  const Int_t NBinM2 = 340;
  const Double_t MinM2 = -1.4;
  const Double_t MaxM2 =  2.0;

  HB1(1, "Status", 21, 0., 21.);

  TString name[Event::nParticle] = { "Pion", "Kaon", "Proton" };

  HB1(10, "pS2s",    NBinP,  MinP, MaxP);
  HB1(11, "ChargexMassSquare", NBinM2, MinM2, MaxM2);
  HB1(12, "beta", 500, 0., 1.);
  HB1(13, "path", 500, 3000., 8000.);
  HB1(14, "stof", 500, 0., 100.);
  HB2(20, "pS2s % ChargexMassSquare", NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);

  HB1(50, "Chisqr S2sTrack", 500, 0., 50.);
  HB1(51, "Xtof S2sTrack", 200, -100., 100.);
  HB1(52, "Ytof S2sTrack", 400, -200., 200.);
  HB1(53, "Utof S2sTrack", 300, -0.30, 0.30);
  HB1(54, "Vtof S2sTrack", 300, -0.30, 0.30);
  HB1(55, "Local Xtof S2sTrack", 700, -700., 700.);
  HB1(56, "Local Ytof S2sTrack", 350, -350., 350.);
  HB1(57, "Local Utof S2sTrack", 300, -0.30, 0.30);
  HB1(58, "Local Vtof S2sTrack", 300, -0.30, 0.30);
  HB1(59, "tofseg S2sTrack", 21, 0, 21);

  for(Int_t i=0;i<NumOfSegBH2;++i){
    HB1(1100+i+1, Form("ChargexMassSquare [BH2-%d]", i+1),
        NBinM2, MinM2, MaxM2);
    HB2(1200+i+1, Form("pS2s %% ChargexMassSquare [BH2-%d]", i+1),
        NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);
  }
  for(Int_t i=0;i<NumOfSegTOF;++i){
    HB1(2100+i+1, Form("ChargexMassSquare [TOF-%d]", i+1),
        NBinM2, MinM2, MaxM2);
    HB2(2200+i+1, Form("pS2s %% ChargexMassSquare [TOF-%d]", i+1),
        NBinM2, MinM2, MaxM2, NBinP, MinP, MaxP);
  }
  // for TOF Param
  for(Int_t ip=0; ip<Event::nParticle; ++ip){
    HB2(10000+ip+1,
        Form("TofTime-%sTime %% TofSeg", name[ip].Data()),
        500, -25., 25., NumOfSegTOF, 0., (Double_t)NumOfSegTOF);
  }
  HB2(20001, "Tof TimeDiff U-D",
      1000, -25., 25., NumOfSegTOF, 0., (Double_t)NumOfSegTOF);
  for(Int_t i=0; i<NumOfSegTOF; ++i){
    HB1(20000+(i+1)*100+1,
	Form("Tof TimeDiff U-D %d", i+1),   1000, -25., 25.);
    HB1(30000+(i+1)*100+1,
        Form("TOF ADC-Pedestal %d-U", i+1), 4000, 0., 4000.);
    HB1(30000+(i+1)*100+11,
        Form("TOF ADC-Pedestal %d-D", i+1), 4000, 0., 4000.);
    HB2(40000+(i+1)*100,
        Form("tCalc-Time %% TOF De %d-U [All]", i+1), 100, 0., 4., 100, -3., 3.);
    HB2(40000+(i+1)*100+10,
        Form("tCalc-Time %% TOF De %d-D [All]", i+1), 100, 0., 4., 100, -3., 3.);
    for(Int_t ip=0; ip<Event::nParticle; ++ip){
      HB1(10000+(i+1)*100+ip+1,
          Form("Tof-%d TofTime-%sTime", i+1, name[ip].Data()),
          500, -25., 25.);
      HB1(33000+(i+1)*100+ip+1,
          Form("TOF ADC-Signal %d-U [%s]", i+1, name[ip].Data()),
          4000, 0., 4000.);
      HB1(33000+(i+1)*100+ip+11,
          Form("TOF ADC-Signal %d-D [%s]", i+1, name[ip].Data()),
          4000, 0., 4000.);
      HB2(40000+(i+1)*100+ip+1,
          Form("tCalc-Time %% TOF De %d-U [%s]", i+1, name[ip].Data()),
          100, 0., 4., 100, -3., 3.);
      HB2(40000+(i+1)*100+ip+11,
          Form("tCalc-Time %% TOF De %d-D [%s]", i+1, name[ip].Data()),
          100, 0., 4., 100, -3., 3.);
    }
  }

  HBTree("shodo", "tree of DstS2sHodoscope");
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("status",     &event.status,      "status/I");

  tree->Branch("nhBh1", &event.nhBh1, "nhBh1/I");
  tree->Branch("csBh1",  event.csBh1, "csBh1[nhBh1]/D");
  tree->Branch("Bh1Seg", event.Bh1Seg,"Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",   event.tBh1,  "tBh1[nhBh1]/D");
  tree->Branch("dtBh1",  event.dtBh1, "dtBh1[nhBh1]/D");
  tree->Branch("deBh1",  event.deBh1, "deBh1[nhBh1]/D");

  tree->Branch("nhBh2", &event.nhBh2, "nhBh2/I");
  tree->Branch("csBh2",  event.csBh2, "csBh2[nhBh2]/D");
  tree->Branch("Bh2Seg", event.Bh2Seg,"Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",   event.tBh2,  "tBh2[nhBh2]/D");
  tree->Branch("dtBh2",  event.dtBh2, "dtBh2[nhBh2]/D");
  tree->Branch("t0Bh2",  event.t0Bh2, "t0Bh2[nhBh2]/D");
  tree->Branch("deBh2",  event.deBh2, "deBh2[nhBh2]/D");

  tree->Branch("btof",   event.btof, "btof[nhBh1]/D");
  tree->Branch("Time0Seg", &event.Time0Seg, "Time0Seg/D");

  tree->Branch("nhTof",  &event.nhTof, "nhTof/I");
  tree->Branch("csTof",   event.csTof, "csTof[nhTof]/D");
  tree->Branch("TofSeg",  event.TofSeg,"TofSeg[nhTof]/D");
  tree->Branch("tTof",    event.tTof,  "tTof[nhTof]/D");
  tree->Branch("dtTof",   event.dtTof, "dtTof[nhTof]/D");
  tree->Branch("deTof",   event.deTof, "deTof[nhTof]/D");
  tree->Branch("tofua",   event.tofua, Form("tofua[%d]/D", NumOfSegTOF));
  tree->Branch("tofda",   event.tofda, Form("tofda[%d]/D", NumOfSegTOF));
  tree->Branch("utTofSeg",  event.utTofSeg,  Form("utTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("dtTofSeg",  event.dtTofSeg,  Form("dtTofSeg[%d][%d]/D", NumOfSegTOF, MaxDepth));
  tree->Branch("udeTofSeg", event.udeTofSeg, Form("udeTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("ddeTofSeg", event.ddeTofSeg, Form("ddeTofSeg[%d]/D", NumOfSegTOF));
  tree->Branch("nhAc1",     &event.nhAc1,     "nhAc1/I");
  tree->Branch("csAc1",      event.csAc1,     "csAc1[nhAc1]/I");
  tree->Branch("Ac1Seg",     event.Ac1Seg,    "Ac1Seg[nhAc1]/D");
  tree->Branch("tAc1",       event.tAc1,      "tAc1[nhAc1]/D");

  tree->Branch("ntSdcOut", &event.ntSdcOut,  "ntSdcOut/I");
  tree->Branch("chisqrSdcOut", event.chisqrSdcOut, "chisqrSdcOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut", event.x0SdcOut, "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut", event.y0SdcOut, "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut", event.u0SdcOut, "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut", event.v0SdcOut, "v0SdcOut[ntSdcOut]/D");
  tree->Branch("ntS2s", &event.ntS2s,  "ntS2s/I");
  tree->Branch("path",      event.path,      "path[ntS2s]/D");
  tree->Branch("pS2s",   event.pS2s,   "pS2s[ntS2s]/D");
  tree->Branch("qS2s",   event.qS2s,   "qS2s[ntS2s]/D");
  tree->Branch("chisqrS2s", event.chisqrS2s, "chisqrS2s[ntS2s]/D");
  tree->Branch("thetaS2s",   event.thetaS2s,  "thetaS2s[ntS2s]/D");
  tree->Branch("xtgtS2s",    event.xtgtS2s,   "xtgtS2s[ntS2s]/D");
  tree->Branch("ytgtS2s",    event.ytgtS2s,   "ytgtS2s[ntS2s]/D");
  tree->Branch("utgtS2s",    event.utgtS2s,   "utgtS2s[ntS2s]/D");
  tree->Branch("vtgtS2s",    event.vtgtS2s,   "vtgtS2s[ntS2s]/D");
  tree->Branch("xtofS2s",    event.xtofS2s,   "xtofS2s[ntS2s]/D");
  tree->Branch("ytofS2s",    event.ytofS2s,   "ytofS2s[ntS2s]/D");
  tree->Branch("utofS2s",    event.utofS2s,   "utofS2s[ntS2s]/D");
  tree->Branch("vtofS2s",    event.vtofS2s,   "vtofS2s[ntS2s]/D");
  tree->Branch("lxtofS2s",    event.lxtofS2s,   "lxtofS2s[ntS2s]/D");
  tree->Branch("lytofS2s",    event.lytofS2s,   "lytofS2s[ntS2s]/D");
  tree->Branch("lutofS2s",    event.lutofS2s,   "lutofS2s[ntS2s]/D");
  tree->Branch("lvtofS2s",    event.lvtofS2s,   "lvtofS2s[ntS2s]/D");
  tree->Branch("tofsegS2s",  event.tofsegS2s, "tofsegS2s[ntS2s]/D");
  tree->Branch("vpx", event.vpx, Form("vpx[%d]/D", NumOfLayersVP));
  tree->Branch("vpy", event.vpy, Form("vpy[%d]/D", NumOfLayersVP));
  tree->Branch("vpu", event.vpu, Form("vpu[%d]/D", NumOfLayersVP));
  tree->Branch("vpv", event.vpv, Form("vpv[%d]/D", NumOfLayersVP));

  tree->Branch("tTofCalc",  event.tTofCalc,  Form("tTofCalc[%d]/D", Event::nParticle));

  tree->Branch("m2Combi", &event.m2Combi, "m2Combi/I");
  tree->Branch("beta",     event.beta,    "beta[m2Combi]/D");
  tree->Branch("stof",     event.stof,    "stof[m2Combi]/D");
  tree->Branch("cstof",    event.cstof,   "cstof[m2Combi]/D");
  tree->Branch("m2",       event.m2,      "m2[m2Combi]/D");

  ////////// Bring Address From Dst
  TTreeCont[kHodoscope]->SetBranchStatus("*", 0);
  TTreeCont[kHodoscope]->SetBranchStatus("trigflag",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("trigpat",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh1Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh1",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Bh2Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tBh2",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("t0Bh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deBh2",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0Seg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTime0",   1);
  TTreeCont[kHodoscope]->SetBranchStatus("Time0",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("CTime0",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("TofSeg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tTof",      1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("deTof",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("tofua",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("tofda",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("utTofSeg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("dtTofSeg",  1);
  TTreeCont[kHodoscope]->SetBranchStatus("udeTofSeg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("ddeTofSeg", 1);
  TTreeCont[kHodoscope]->SetBranchStatus("nhAc1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("csAc1",     1);
  TTreeCont[kHodoscope]->SetBranchStatus("Ac1Seg",    1);
  TTreeCont[kHodoscope]->SetBranchStatus("tAc1",      1);

  TTreeCont[kHodoscope]->SetBranchAddress("trigflag", src.trigflag);
  TTreeCont[kHodoscope]->SetBranchAddress("trigpat",  src.trigpat);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh1",  &src.nhBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh1",  src.csBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh1Seg", src.Bh1Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh1",   src.tBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh1",  src.dtBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh1",  src.deBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh2", &src.nhBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh2",  src.csBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh2Seg",src.Bh2Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh2",  src.tBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh2", src.dtBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("t0Bh2", src.t0Bh2);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh2", src.deBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0Seg", &src.Time0Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("deTime0",  &src.deTime0);
  TTreeCont[kHodoscope]->SetBranchAddress("Time0",    &src.Time0);
  TTreeCont[kHodoscope]->SetBranchAddress("CTime0",   &src.CTime0);
  TTreeCont[kHodoscope]->SetBranchAddress("nhTof", &src.nhTof);
  TTreeCont[kHodoscope]->SetBranchAddress("csTof", src.csTof);
  TTreeCont[kHodoscope]->SetBranchAddress("TofSeg",src.TofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tTof",  src.tTof);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTof", src.dtTof);
  TTreeCont[kHodoscope]->SetBranchAddress("deTof", src.deTof);
  TTreeCont[kHodoscope]->SetBranchAddress("tofua", src.tofua);
  TTreeCont[kHodoscope]->SetBranchAddress("tofda", src.tofda);
  TTreeCont[kHodoscope]->SetBranchAddress("utTofSeg", src.utTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTofSeg", src.dtTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("udeTofSeg", src.udeTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("ddeTofSeg", src.ddeTofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("nhAc1", &src.nhAc1);
  TTreeCont[kHodoscope]->SetBranchAddress("csAc1", src.csAc1);
  TTreeCont[kHodoscope]->SetBranchAddress("Ac1Seg",src.Ac1Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tAc1",  src.tAc1);

  TTreeCont[kS2sTracking]->SetBranchStatus("*",      0);
  TTreeCont[kS2sTracking]->SetBranchStatus("ntSdcOut", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("chisqrSdcOut", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("x0SdcOut", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("y0SdcOut", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("u0SdcOut", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("v0SdcOut", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("ntS2s",     1);
  TTreeCont[kS2sTracking]->SetBranchStatus("path",         1);
  TTreeCont[kS2sTracking]->SetBranchStatus("pS2s",      1);
  TTreeCont[kS2sTracking]->SetBranchStatus("qS2s",      1);
  TTreeCont[kS2sTracking]->SetBranchStatus("chisqrS2s", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("thetaS2s",  1);
  TTreeCont[kS2sTracking]->SetBranchStatus("xtgtS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("ytgtS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("utgtS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("vtgtS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("xtofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("ytofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("utofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("vtofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("lxtofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("lytofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("lutofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("lvtofS2s",   1);
  TTreeCont[kS2sTracking]->SetBranchStatus("tofsegS2s", 1);
  TTreeCont[kS2sTracking]->SetBranchStatus("vpx",          1);
  TTreeCont[kS2sTracking]->SetBranchStatus("vpy",          1);
  TTreeCont[kS2sTracking]->SetBranchStatus("vpu",          1);
  TTreeCont[kS2sTracking]->SetBranchStatus("vpv",          1);

  TTreeCont[kS2sTracking]->SetBranchAddress("ntSdcOut", &src.ntSdcOut);
  TTreeCont[kS2sTracking]->SetBranchAddress("chisqrSdcOut", src.chisqrSdcOut);
  TTreeCont[kS2sTracking]->SetBranchAddress("x0SdcOut", src.x0SdcOut);
  TTreeCont[kS2sTracking]->SetBranchAddress("y0SdcOut", src.y0SdcOut);
  TTreeCont[kS2sTracking]->SetBranchAddress("u0SdcOut", src.u0SdcOut);
  TTreeCont[kS2sTracking]->SetBranchAddress("v0SdcOut", src.v0SdcOut);
  TTreeCont[kS2sTracking]->SetBranchAddress("ntS2s", &src.ntS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("path",     src.path);
  TTreeCont[kS2sTracking]->SetBranchAddress("pS2s",  src.pS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("qS2s",  src.qS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("chisqrS2s", src.chisqrS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("thetaS2s", src.thetaS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("xtgtS2s", src.xtgtS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("ytgtS2s", src.ytgtS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("utgtS2s", src.utgtS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("vtgtS2s", src.vtgtS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("xtofS2s",   src.xtofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("ytofS2s",   src.ytofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("utofS2s",   src.utofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("vtofS2s",   src.vtofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("lxtofS2s",   src.lxtofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("lytofS2s",   src.lytofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("lutofS2s",   src.lutofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("lvtofS2s",   src.lvtofS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("tofsegS2s", src.tofsegS2s);
  TTreeCont[kS2sTracking]->SetBranchAddress("vpx",        src.vpx);
  TTreeCont[kS2sTracking]->SetBranchAddress("vpy",        src.vpy);
  TTreeCont[kS2sTracking]->SetBranchAddress("vpu",        src.vpu);
  TTreeCont[kS2sTracking]->SetBranchAddress("vpv",        src.vpv);

  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")   &&
     InitializeParameter<UserParamMan>("USER") &&
     InitializeParameter<HodoPHCMan>("HDPHC"));
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}
