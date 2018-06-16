/**
 *  file: DstXiAna.cc
 *  date: 2017.04.10
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <signal.h>
#include <sstream>
#include <string>

#include <filesystem_util.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "DCGeomMan.hh"
#include "Kinematics.hh"
#include "LorentzVector.hh"
#include "MathTools.hh"
#include "NuclearMass.hh"
#include "RootHelper.hh"
#include "ThreeVector.hh"
#include "UserParamMan.hh"

#include "DstHelper.hh"

namespace
{
  using namespace root;
  using namespace dst;
  const std::string& class_name("DstXiAna");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kKKAna, kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[KKAna]", "[OutFile]" };
  std::vector<TString> TreeName =
    { "", "", "kk", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;

  // HodoCut
  const double MinBH2MeanTime = -1.0;
  const double MaxBH2MeanTime =  1.0;
  const double MinBTOF        = -1.0;
  const double MaxBTOF        =  1.0;
  // Kaon Minus
  const double MaxChisqrK18    =  10.; // chisqr of BcOut
  const double MaxChisqrKurama = 200.;
  // Kaon Plus
  const double MinKuramaP    = 0.90; // for Kaon
  const double MaxKuramaP    = 1.50; // for Kaon
  const double MinMassSquare = 0.05; // for Kaon
  const double MaxMassSquare = 0.45; // for Kaon
  // KK reaction
  const double MinVertexX   =  -50.;
  const double MaxVertexX   =   50.;
  const double MinVertexY   =  -30.;
  const double MaxVertexY   =   30.;
  const double MinVertexZ   =  -80.;
  const double MaxVertexZ   =   80.;
  const double MaxCloseDist =   10.;
  // Xi analysis
  const double MinDeSSDXi     = 15000.;
  const double MinVertexXiX   = -30.;
  const double MaxVertexXiX   =  30.;
  const double MinVertexXiY   = -20.;
  const double MaxVertexXiY   =  20.;
  const double MinVertexXiZ   = -40.;
  const double MaxVertexXiZ   =  40.;
  const double MaxCloseDistXi =  10.;
  const double ResMissU       = 0.20; // 1-sigma
  const double ResMissV       = 0.19; // 1-sigma
  // const double ResMissU       = 0.08; // 1-sigma (CH2)
  // const double ResMissV       = 0.08; // 1-sigma (CH2)
  const double ResVertexX     = 1.12*5.0; // 1-sigma
  const double ResVertexY     = 1.50*5.0; // 1-sigma
  const double ResVertexZ     = 7.41*5.0; // 1-sigma

}

//_____________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  int spill;

  //Trigger
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Hodoscope
  int    nhBh1;
  int    csBh1[NumOfSegBH1];
  double Bh1Seg[NumOfSegBH1];
  double tBh1[NumOfSegBH1];
  double dtBh1[NumOfSegBH1];
  double deBh1[NumOfSegBH1];
  double btof[NumOfSegBH1];

  int    nhBh2;
  int    csBh2[NumOfSegBH2];
  double Bh2Seg[NumOfSegBH2];
  double tBh2[NumOfSegBH2];
  double t0Bh2[NumOfSegBH2];
  double dtBh2[NumOfSegBH2];
  double deBh2[NumOfSegBH2];

  int    nhTof;
  int    csTof[NumOfSegTOF];
  double TofSeg[NumOfSegTOF];
  double tTof[NumOfSegTOF];
  double dtTof[NumOfSegTOF];
  double deTof[NumOfSegTOF];

  int    nhBac;
  double BacSeg[NumOfSegBAC];
  double tBac[NumOfSegBAC];
  double deBac[NumOfSegBAC];

  int    nhPvac;
  double PvacSeg[NumOfSegPVAC];
  double tPvac[NumOfSegPVAC];
  double dePvac[NumOfSegPVAC];

  int    nhFac;
  double FacSeg[NumOfSegFAC];
  double tFac[NumOfSegFAC];
  double deFac[NumOfSegFAC];

  //Fiber
  int    nhBft;
  int    csBft[NumOfSegBFT];
  double tBft[NumOfSegBFT];
  double wBft[NumOfSegBFT];
  double BftPos[NumOfSegBFT];
  double BftSeg[NumOfSegBFT];
  int    nhSch;
  int    csSch[NumOfSegSCH];
  double tSch[NumOfSegSCH];
  double wSch[NumOfSegSCH];
  double SchPos[NumOfSegSCH];
  double SchSeg[NumOfSegSCH];
  int    nhFbh;
  int    csFbh[NumOfSegCFBH];
  double tFbh[NumOfSegCFBH];
  double wFbh[NumOfSegCFBH];
  double FbhPos[NumOfSegCFBH];
  double FbhSeg[NumOfSegCFBH];

  //DC Beam
  int ntBcOut;
  int nlBcOut;
  int nhBcOut[MaxHits];
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  double xtgtBcOut[MaxHits];
  double ytgtBcOut[MaxHits];
  double xbh2BcOut[MaxHits];
  double ybh2BcOut[MaxHits];

  int    ntK18;
  int    nhK18[MaxHits];
  double chisqrK18[MaxHits];
  double pK18[MaxHits];
  double xtgtK18[MaxHits];
  double ytgtK18[MaxHits];
  double utgtK18[MaxHits];
  double vtgtK18[MaxHits];
  double thetaK18[MaxHits];

  //DC KURAMA
  int ntSdcIn;
  int nlSdcIn;
  int nhSdcIn[MaxHits];
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];

  int ntSdcOut;
  int nlSdcOut;
  int nhSdcOut[MaxHits];
  double chisqrSdcOut[MaxHits];
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double x0SdcOut[MaxHits];
  double y0SdcOut[MaxHits];

  int    ntKurama;
  int    nhKurama[MaxHits];
  double chisqrKurama[MaxHits];
  double stof[MaxHits];
  double path[MaxHits];
  double pKurama[MaxHits];
  double qKurama[MaxHits];
  double m2[MaxHits];
  double xtgtKurama[MaxHits];
  double ytgtKurama[MaxHits];
  double utgtKurama[MaxHits];
  double vtgtKurama[MaxHits];
  double thetaKurama[MaxHits];

  int    ntXi;
  double x0Xi[MaxHits];
  double y0Xi[MaxHits];
  double u0Xi[MaxHits];
  double v0Xi[MaxHits];
  double thetaXi[MaxHits];
  double deXi[NumOfLayersSsdIn][MaxHits];
  double deKaon[NumOfLayersSsdIn][MaxHits];

  //Reaction
  int    nKn;
  int    nKp;
  int    nKK;
  double vtx[MaxHits];
  double vty[MaxHits];
  double vtz[MaxHits];
  double closeDist[MaxHits];
  double theta[MaxHits];
  double MissP[MaxHits];
  double MissPu[MaxHits];
  double MissPv[MaxHits];
  double MissMass[MaxHits];
  double MissMassC[MaxHits];
  double thetaCM[MaxHits];
  double costCM[MaxHits];

  double xkn[MaxHits];
  double ykn[MaxHits];
  double ukn[MaxHits];
  double vkn[MaxHits];
  double xkp[MaxHits];
  double ykp[MaxHits];
  double ukp[MaxHits];
  double vkp[MaxHits];

  double pOrg[MaxHits];
  double pCalc[MaxHits];

  // Xi
  int    nMM;
  int    nXi;
  int    nKXi;
  double vtxXi[MaxHits];
  double vtyXi[MaxHits];
  double vtzXi[MaxHits];
  double closeXi[MaxHits];

  // EMC
  double xpos;
  double ypos;
  int    state;
};

//_____________________________________________________________________
struct Src
{
  int runnum;
  int evnum;
  int spill;

  //Trigger
  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  //Hodoscope
  int    nhBh1;
  int    csBh1[NumOfSegBH1];
  double Bh1Seg[NumOfSegBH1];
  double tBh1[NumOfSegBH1];
  double dtBh1[NumOfSegBH1];
  double deBh1[NumOfSegBH1];
  double btof[NumOfSegBH1];

  int    nhBh2;
  int    csBh2[NumOfSegBH2];
  double Bh2Seg[NumOfSegBH2];
  double tBh2[NumOfSegBH2];
  double t0Bh2[NumOfSegBH2];
  double dtBh2[NumOfSegBH2];
  double deBh2[NumOfSegBH2];

  int    nhTof;
  int    csTof[NumOfSegTOF];
  double TofSeg[NumOfSegTOF];
  double tTof[NumOfSegTOF];
  double dtTof[NumOfSegTOF];
  double deTof[NumOfSegTOF];

  int    nhBac;
  double BacSeg[NumOfSegBAC];
  double tBac[NumOfSegBAC];
  double deBac[NumOfSegBAC];

  int    nhPvac;
  double PvacSeg[NumOfSegPVAC];
  double tPvac[NumOfSegPVAC];
  double dePvac[NumOfSegPVAC];

  int    nhFac;
  double FacSeg[NumOfSegFAC];
  double tFac[NumOfSegFAC];
  double deFac[NumOfSegFAC];

  //Fiber
  int    nhBft;
  int    csBft[NumOfSegBFT];
  double tBft[NumOfSegBFT];
  double wBft[NumOfSegBFT];
  double BftPos[NumOfSegBFT];
  double BftSeg[NumOfSegBFT];
  int    nhSch;
  int    csSch[NumOfSegSCH];
  double tSch[NumOfSegSCH];
  double wSch[NumOfSegSCH];
  double SchPos[NumOfSegSCH];
  double SchSeg[NumOfSegSCH];
  int    nhFbh;
  int    csFbh[NumOfSegCFBH];
  double tFbh[NumOfSegCFBH];
  double wFbh[NumOfSegCFBH];
  double FbhPos[NumOfSegCFBH];
  double FbhSeg[NumOfSegCFBH];

  //DC Beam
  int ntBcOut;
  int nlBcOut;
  int nhBcOut[MaxHits];
  double chisqrBcOut[MaxHits];
  double x0BcOut[MaxHits];
  double y0BcOut[MaxHits];
  double u0BcOut[MaxHits];
  double v0BcOut[MaxHits];
  double xtgtBcOut[MaxHits];
  double ytgtBcOut[MaxHits];
  double xbh2BcOut[MaxHits];
  double ybh2BcOut[MaxHits];

  int    ntK18;
  int    nhK18[MaxHits];
  double chisqrK18[MaxHits];
  double pK18[MaxHits];
  double xtgtK18[MaxHits];
  double ytgtK18[MaxHits];
  double utgtK18[MaxHits];
  double vtgtK18[MaxHits];
  double thetaK18[MaxHits];

  //DC KURAMA
  int ntSdcIn;
  int nlSdcIn;
  int nhSdcIn[MaxHits];
  double chisqrSdcIn[MaxHits];
  double x0SdcIn[MaxHits];
  double y0SdcIn[MaxHits];
  double u0SdcIn[MaxHits];
  double v0SdcIn[MaxHits];

  int ntSdcOut;
  int nlSdcOut;
  int nhSdcOut[MaxHits];
  double chisqrSdcOut[MaxHits];
  double u0SdcOut[MaxHits];
  double v0SdcOut[MaxHits];
  double x0SdcOut[MaxHits];
  double y0SdcOut[MaxHits];

  int    ntKurama;
  int    nhKurama[MaxHits];
  double chisqrKurama[MaxHits];
  double stof[MaxHits];
  double path[MaxHits];
  double pKurama[MaxHits];
  double qKurama[MaxHits];
  double m2[MaxHits];
  double xtgtKurama[MaxHits];
  double ytgtKurama[MaxHits];
  double utgtKurama[MaxHits];
  double vtgtKurama[MaxHits];
  double thetaKurama[MaxHits];

  int    ntXi;
  double x0Xi[MaxHits];
  double y0Xi[MaxHits];
  double u0Xi[MaxHits];
  double v0Xi[MaxHits];
  double thetaXi[MaxHits];
  double deXi[NumOfLayersSsdIn][MaxHits];
  double deKaon[NumOfLayersSsdIn][MaxHits];

  //Reaction
  int    nKn;
  int    nKp;
  int    nKK;
  double vtx[MaxHits];
  double vty[MaxHits];
  double vtz[MaxHits];
  double closeDist[MaxHits];
  double theta[MaxHits];
  double MissP[MaxHits];
  double MissPu[MaxHits];
  double MissPv[MaxHits];
  double MissMass[MaxHits];
  double thetaCM[MaxHits];
  double costCM[MaxHits];

  double xkn[MaxHits];
  double ykn[MaxHits];
  double ukn[MaxHits];
  double vkn[MaxHits];
  double xkp[MaxHits];
  double ykp[MaxHits];
  double ukp[MaxHits];
  double vkp[MaxHits];

  double pOrg[MaxHits];
  double pCalc[MaxHits];

  // EMC
  double xpos;
  double ypos;
  int    state;

};

//_____________________________________________________________________
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
  event.runnum = 0;
  event.evnum  = 0;
  event.spill  = 0;
  event.xpos   = -9999.;
  event.ypos   = -9999.;
  event.state  = -1;

  //Trigger
  for( int i=0; i<NumOfSegTrig; ++i ){
    event.trigpat[i]  = -1;
    event.trigflag[i] = -1;
  }

  //Hodoscope
  event.nhBh2  = 0;
  event.nhBh1  = 0;
  event.nhTof  = 0;
  event.nhBac  = 0;
  event.nhPvac = 0;
  event.nhFac  = 0;

  for( int i=0; i<NumOfSegBH1; i++ ){
    event.csBh1[i]  = 0;
    event.Bh1Seg[i] = -1;
    event.tBh1[i]   = -9999.0;
    event.dtBh1[i]  = -9999.0;
    event.deBh1[i]  = -9999.0;
    event.btof[i]   = -9999.0;
  }
  for( int i=0; i<NumOfSegBH2; i++ ){
    event.csBh2[i]  = 0;
    event.Bh2Seg[i] = -1;
    event.tBh2[i]   = -9999.0;
    event.t0Bh2[i]  = -9999.0;
    event.dtBh2[i]  = -9999.0;
    event.deBh2[i]  = -9999.0;
  }
  for( int i=0; i<NumOfSegTOF; i++ ){
    event.csTof[i]  = 0;
    event.TofSeg[i] = -1;
    event.tTof[i]   = -9999.0;
    event.dtTof[i]  = -9999.0;
    event.deTof[i]  = -9999.0;
  }
  for( int i=0; i<NumOfSegBAC; i++ ){
    event.BacSeg[i] = -1;
    event.tBac[i]   = -9999.0;
    event.deBac[i]  = -9999.0;
  }
  for( int i=0; i<NumOfSegPVAC; i++ ){
    event.PvacSeg[i] = -1;
    event.tPvac[i]   = -9999.0;
    event.dePvac[i]  = -9999.0;
  }
  for( int i=0; i<NumOfSegFAC; i++ ){
    event.FacSeg[i] = -1;
    event.tFac[i]   = -9999.0;
    event.deFac[i]  = -9999.0;
  }

  //Fiber
  event.nhBft = 0;
  event.nhSch = 0;
  event.nhFbh = 0;
  for( int i=0; i<NumOfSegBFT; i++ ){
    event.csBft[i]  = 0;
    event.tBft[i]   = -9999.;
    event.wBft[i]   = -9999.;
    event.BftPos[i] = -9999.;
  }
  for( int i=0; i<NumOfSegSCH; i++ ){
    event.csSch[i]   = 0;
    event.tSch[i]    = -9999.;
    event.wSch[i]    = -9999.;
    event.SchPos[i]  = -9999.;
  }
  for( int i=0; i<NumOfSegFBH; i++ ){
    event.csFbh[i]  = 0;
    event.tFbh[i]   = -9999.;
    event.wFbh[i]   = -9999.;
    event.FbhPos[i] = -9999.;
  }

  //DC
  event.nlBcOut  = 0;
  event.nlSdcIn  = 0;
  event.nlSdcOut = 0;
  event.ntBcOut  = 0;
  event.ntSdcIn  = 0;
  event.ntSdcOut = 0;
  event.ntK18    = 0;
  event.ntKurama = 0;
  event.ntXi     = 0;

  //Beam DC
  for( int i=0; i<MaxHits; i++){
    event.nhBcOut[i]     = 0;
    event.chisqrBcOut[i] = -1.0;
    event.x0BcOut[i] = -9999.0;
    event.y0BcOut[i] = -9999.0;
    event.u0BcOut[i] = -9999.0;
    event.v0BcOut[i] = -9999.0;

    event.xtgtBcOut[i] = -9999.0;
    event.ytgtBcOut[i] = -9999.0;
    event.xbh2BcOut[i] = -9999.0;
    event.ybh2BcOut[i] = -9999.0;
  }

  for( int i=0; i<MaxHits; i++){
    event.nhK18[i]     = 0;
    event.chisqrK18[i] = -1.0;
    event.xtgtK18[i]   = -9999.;
    event.ytgtK18[i]   = -9999.;
    event.utgtK18[i]   = -9999.;
    event.vtgtK18[i]   = -9999.;
    event.pK18[i]      = -9999.;
    event.thetaK18[i]  = -9999.;
  }

  //KURAMA DC
  for( int i=0; i<MaxHits; i++){
    event.nhSdcIn[i]     = 0;
    event.chisqrSdcIn[i] = -1.;
    event.x0SdcIn[i]     = -9999.;
    event.y0SdcIn[i]     = -9999.;
    event.u0SdcIn[i]     = -9999.;
    event.v0SdcIn[i]     = -9999.;

    event.nhSdcOut[i]     = 0;
    event.chisqrSdcOut[i] = -1.;
    event.u0SdcOut[i] = -9999.;
    event.v0SdcOut[i] = -9999.;
    event.x0SdcOut[i] = -9999.;
    event.y0SdcOut[i] = -9999.;

    event.nhKurama[i]     = 0;
    event.chisqrKurama[i] = -1.;
    event.xtgtKurama[i] = -9999.;
    event.ytgtKurama[i] = -9999.;
    event.utgtKurama[i] = -9999.;
    event.vtgtKurama[i] = -9999.;

    event.pKurama[i]  = -9999.;
    event.qKurama[i]  = -9999.;
    event.stof[i]     = -9999.;
    event.path[i]     = -9999.;
    event.m2[i]       = -9999.;
    event.thetaKurama[i]   = -9999.;
  }

  for( int i=0; i<MaxHits; ++i ){
    event.x0Xi[i]    = -9999.;
    event.y0Xi[i]    = -9999.;
    event.u0Xi[i]    = -9999.;
    event.v0Xi[i]    = -9999.;
    event.thetaXi[i] = -9999.;
    for( int j=0; j<NumOfLayersSsdIn; ++j ){
      event.deXi[j][i]   = -9999.;
      event.deKaon[j][i] = -9999.;
    }
    event.vtxXi[i] = -9999.;
    event.vtyXi[i] = -9999.;
    event.vtzXi[i] = -9999.;
    event.closeXi[i] = -9999.;
  }

  //Reaction
  event.nKn  = 0;
  event.nKp  = 0;
  event.nKK  = 0;
  event.nMM  = 0;
  event.nXi  = 0;
  event.nKXi = 0;

  for( int i=0; i<MaxHits; ++i ){
    event.vtx[i]       = -9999.;
    event.vty[i]       = -9999.;
    event.vtz[i]       = -9999.;
    event.closeDist[i] = -9999.;
    event.theta[i]     = -9999.;
    event.thetaCM[i]   = -9999.;
    event.costCM[i]    = -9999.;
    event.MissP[i]     = -9999.;
    event.MissPu[i]    = -9999.;
    event.MissPv[i]    = -9999.;
    event.MissMass[i]  = -9999.;
    event.MissMassC[i] = -9999.;

    event.xkn[i] = -9999.0;
    event.ykn[i] = -9999.0;
    event.ukn[i] = -9999.0;
    event.vkn[i] = -9999.0;
    event.xkp[i] = -9999.0;
    event.ykp[i] = -9999.0;
    event.ukp[i] = -9999.0;
    event.vkp[i] = -9999.0;
    event.pOrg[i] = -9999.0;
    event.pCalc[i] = -9999.0;
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

  static const double zTarget =
    gGeom.GetGlobalPosition("Target").z();
  static const double zEmulsion =
    gGeom.GetGlobalPosition("Emulsion").z();

  static const double KaonMass   = pdg::KaonMass();
  static const double ProtonMass = pdg::ProtonMass();
  static const double XiMass     = pdg::XiMass();
  static const double CarbonMass = 12.*NuclearMass::AMU();

  if( ievent%10000==0 ){
    std::cout << "#D " << func_name << " Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  {
    const std::size_t n = TTreeCont.size();
    for( std::size_t i=0; i<n; ++i ){
      if( TTreeCont[i] ) TTreeCont[i]->GetEntry( ievent );
    }
  }

  event.runnum   = src.runnum;
  event.evnum    = src.evnum;
  event.spill    = src.spill;
  event.xpos     = src.xpos;
  event.ypos     = src.ypos;
  event.state    = src.state;
  event.ntBcOut  = src.ntBcOut;
  event.ntSdcIn  = src.ntSdcIn;
  event.ntSdcOut = src.ntSdcOut;
  event.ntKurama = src.ntKurama;
  event.ntK18    = src.ntK18;
  // event.ntXi     = src.ntXi;
  event.nhBft    = src.nhBft;
  event.nhBh1    = src.nhBh1;
  event.nhBh2    = src.nhBh2;
  event.nhBac    = src.nhBac;
  event.nhFbh    = src.nhFbh;
  event.nhPvac   = src.nhPvac;
  event.nhFac    = src.nhFac;
  event.nhSch    = src.nhSch;
  event.nhTof    = src.nhTof;
  int ntBcOut  = event.ntBcOut;
  int ntSdcIn  = event.ntSdcIn;
  int ntSdcOut = event.ntSdcOut;
  int ntKurama = event.ntKurama;
  int ntK18    = event.ntK18;
  int ntXi     = src.ntXi;
  int nhBft    = event.nhBft;
  int nhBh1    = event.nhBh1;
  int nhBh2    = event.nhBh2;
  int nhBac    = event.nhBac;
  int nhFbh    = event.nhFbh;
  int nhPvac   = event.nhPvac;
  int nhFac    = event.nhFac;
  int nhSch    = event.nhSch;
  int nhTof    = event.nhTof;

#if 0
  std::cout << "#D DebugPrint ========================" << std::endl
	    << " event  : " << std::setw(6) << ievent   << std::endl
	    << " BcOut  : " << std::setw(6) << ntBcOut  << std::endl
	    << " SdcIn  : " << std::setw(6) << ntSdcIn  << std::endl
	    << " SdcOut : " << std::setw(6) << ntSdcOut << std::endl
	    << " Kurama : " << std::setw(6) << ntKurama << std::endl
	    << " K18    : " << std::setw(6) << ntK18    << std::endl
	    << " Xi     : " << std::setw(6) << ntXi     << std::endl
	    << " BFT    : " << std::setw(6) << nhBft    << std::endl
	    << " BH1    : " << std::setw(6) << nhBh1    << std::endl
	    << " BH2    : " << std::setw(6) << nhBh2    << std::endl
	    << " BAC    : " << std::setw(6) << nhBac    << std::endl
	    << " FBH    : " << std::setw(6) << nhFbh    << std::endl
	    << " SCH    : " << std::setw(6) << nhSch    << std::endl
	    << " PVAC   : " << std::setw(6) << nhPvac   << std::endl
	    << " FAC    : " << std::setw(6) << nhFac    << std::endl
	    << " TOF    : " << std::setw(6) << nhTof    << std::endl
	    << std::endl;
#endif

  // TFlag
  for( int i=0; i<NumOfSegTrig; ++i ){
    int seg = src.trigpat[i];
    if( seg<=0 ) continue;
    int tdc = src.trigflag[seg-1];
    if( tdc<=0 ) continue;
    event.trigpat[i] = seg;
    event.trigflag[seg-1] = tdc;
  }

  HF1( 1, 0. );
  if( event.nhBh1==0 ) return false;
  HF1( 1, 1. );
  if( event.nhBh2==0 ) return false;
  HF1( 1, 2. );
  if( event.nhTof==0 ) return false;
  HF1( 1, 3. );
  // if( event.nhSch==0 ) return false;
  HF1( 1, 4. );
  // if( event.nhFbh==0 ) return false;
  HF1( 1, 5. );
  if( event.ntK18==0 ) return false;
  HF1( 1, 6. );
  if( event.ntKurama==0 ) return false;
  HF1( 1, 7. );

  // BH2
  HF1( 101, nhBh2 );
  // double time0 = -9999.;
  for( int i=0; i<nhBh2; ++i ){
    event.csBh2[i]  = src.csBh2[i];
    event.Bh2Seg[i] = src.Bh2Seg[i];
    event.tBh2[i]   = src.tBh2[i];
    event.t0Bh2[i]  = src.t0Bh2[i];
    event.dtBh2[i]  = src.dtBh2[i];
    event.deBh2[i]  = src.deBh2[i];
    // time0 = event.tBh2[i];
    HF1( 102, event.csBh2[i]  ); HF1( 103, event.Bh2Seg[i] );
    HF1( 104, event.tBh2[i]   ); HF1( 105, event.deBh2[i]  );
    HF1( 106, event.t0Bh2[i]  );
  }

  // BFT
  for( int i=0; i<nhBft; ++i ){
    event.csBft[i]  = src.csBft[i];
    event.tBft[i]   = src.tBft[i];
    event.wBft[i]   = src.wBft[i];
    event.BftPos[i] = src.BftPos[i];
  }

  // BH1
  HF1( 201, nhBh1 );
  double btof = -9999.;
  for( int i=0; i<nhBh1; ++i ){
    event.csBh1[i]  = src.csBh1[i];
    event.Bh1Seg[i] = src.Bh1Seg[i];
    event.tBh1[i]   = src.tBh1[i];
    event.dtBh1[i]  = src.dtBh1[i];
    event.deBh1[i]  = src.deBh1[i];
    event.btof[i]   = src.btof[i];
    if( std::abs( event.btof[i] )< std::abs( btof ) ){
      btof = event.btof[i];
    }
    HF1( 202, event.csBh1[i]  ); HF1( 203, event.Bh1Seg[i] );
    HF1( 204, event.tBh1[i]   ); HF1( 205, event.deBh1[i]  );
  }
  HF1( 206, btof );

  // BAC
  for( int i=0; i<nhBac; ++i ){
    event.BacSeg[i] = src.BacSeg[i];
    event.tBac[i]   = src.tBac[i];
    event.deBac[i]  = src.deBac[i];
  }

  // PVAC
  for( int i=0; i<nhPvac; ++i ){
    event.PvacSeg[i] = src.PvacSeg[i];
    event.tPvac[i]   = src.tPvac[i];
    event.dePvac[i]  = src.dePvac[i];
  }

  // BAC
  for( int i=0; i<nhFac; ++i ){
    event.FacSeg[i] = src.FacSeg[i];
    event.tFac[i]   = src.tFac[i];
    event.deFac[i]  = src.deFac[i];
  }

  // FBH
  for( int i=0; i<nhFbh; ++i ){
    event.csFbh[i]  = src.csFbh[i];
    event.tFbh[i]   = src.tFbh[i];
    event.wFbh[i]   = src.wFbh[i];
    event.FbhPos[i] = src.FbhPos[i];
  }

  // PVAC
  for( int i=0; i<nhPvac; ++i ){
    event.PvacSeg[i] = src.PvacSeg[i];
    event.tPvac[i]   = src.tPvac[i];
    event.dePvac[i]  = src.dePvac[i];
  }

  // FAC
  for( int i=0; i<nhFac; ++i ){
    event.FacSeg[i] = src.FacSeg[i];
    event.tFac[i]   = src.tFac[i];
    event.deFac[i]  = src.deFac[i];
  }

  // SCH
  for( int i=0; i<nhSch; ++i ){
    event.csSch[i]  = src.csSch[i];
    event.tSch[i]   = src.tSch[i];
    event.wSch[i]   = src.wSch[i];
    event.SchPos[i] = src.SchPos[i];
  }

  // TOF
  HF1( 301, nhTof );
  for( int i=0; i<nhTof; ++i ){
    event.csTof[i]  = src.csTof[i];
    event.TofSeg[i] = src.TofSeg[i];
    event.tTof[i]   = src.tTof[i];
    event.dtTof[i]  = src.dtTof[i];
    event.deTof[i]  = src.deTof[i];
    HF1( 302, event.csTof[i] ); HF1( 303, event.TofSeg[i] );
    HF1( 304, event.tTof[i]  ); HF1( 305, event.deTof[i]  );
  }

  // HodoCut
  // if( time0<MinBH2MeanTime || MaxBH2MeanTime<time0 ) return false;
  // if( btof<MinBTOF || MaxBTOF<btof ) return false;

  HF1( 1, 8. );

  HF1( 111, nhBh2 );
  for( int i=0; i<nhBh2; ++i ){
    HF1( 112, event.csBh2[i]  ); HF1( 113, event.Bh2Seg[i] );
    HF1( 114, event.tBh2[i]   ); HF1( 115, event.deBh2[i]  );
    HF1( 116, event.t0Bh2[i]  );
  }
  HF1( 211, nhBh1 ); HF1( 216, btof );
  for( int i=0; i<nhBh1; ++i ){
    HF1( 212, event.csBh1[i]  ); HF1( 213, event.Bh1Seg[i] );
    HF1( 214, event.tBh1[i]   ); HF1( 215, event.deBh1[i]  );
  }
  HF1( 311, nhTof );
  for( int i=0; i<nhTof; ++i ){
    HF1( 312, event.csTof[i] ); HF1( 313, event.TofSeg[i] );
    HF1( 314, event.tTof[i]  ); HF1( 315, event.deTof[i]  );
  }

  ////////// BcOut
  HF1( 1001, double(ntBcOut) );
  event.nlBcOut = src.nlBcOut;
  for( int i=0; i<ntBcOut; ++i ){
    event.nhBcOut[i]     = src.nhBcOut[i];
    event.chisqrBcOut[i] = src.chisqrBcOut[i];
    event.x0BcOut[i]     = src.x0BcOut[i];
    event.y0BcOut[i]     = src.y0BcOut[i];
    event.u0BcOut[i]     = src.u0BcOut[i];
    event.v0BcOut[i]     = src.v0BcOut[i];
    HF1( 1002, event.nhBcOut[i] );
    HF1( 1003, event.chisqrBcOut[i] );
    HF1( 1004, event.x0BcOut[i] );
    HF1( 1005, event.y0BcOut[i] );
    HF1( 1006, event.u0BcOut[i] );
    HF1( 1007, event.v0BcOut[i] );
    HF2( 1008, event.x0BcOut[i], event.u0BcOut[i] );
    HF2( 1009, event.y0BcOut[i], event.v0BcOut[i] );
    HF2( 1010, event.x0BcOut[i], event.y0BcOut[i] );
  }

  ////////// SdcIn
  HF1( 1101, double(ntSdcIn) );
  event.nlSdcIn = src.nlSdcIn;
  for( int i=0; i<ntSdcIn; ++i ){
    event.nhSdcIn[i]     = src.nhSdcIn[i];
    event.chisqrSdcIn[i] = src.chisqrSdcIn[i];
    event.x0SdcIn[i]     = src.x0SdcIn[i];
    event.y0SdcIn[i]     = src.y0SdcIn[i];
    event.u0SdcIn[i]     = src.u0SdcIn[i];
    event.v0SdcIn[i]     = src.v0SdcIn[i];
    HF1( 1102, event.nhSdcIn[i] );
    HF1( 1103, event.chisqrSdcIn[i] );
    HF1( 1104, event.x0SdcIn[i] );
    HF1( 1105, event.y0SdcIn[i] );
    HF1( 1106, event.u0SdcIn[i] );
    HF1( 1107, event.v0SdcIn[i] );
    HF2( 1108, event.x0SdcIn[i], event.u0SdcIn[i] );
    HF2( 1109, event.y0SdcIn[i], event.v0SdcIn[i] );
    HF2( 1110, event.x0SdcIn[i], event.y0SdcIn[i] );
 }

  ////////// SdcOut
  HF1( 1201, double(ntSdcIn) );
  event.nlSdcOut = src.nlSdcOut;
  for( int i=0; i<ntSdcOut; ++i ){
    event.nhSdcOut[i]     = src.nhSdcOut[i];
    event.chisqrSdcOut[i] = src.chisqrSdcOut[i];
    event.x0SdcOut[i]     = src.x0SdcOut[i];
    event.y0SdcOut[i]     = src.y0SdcOut[i];
    event.u0SdcOut[i]     = src.u0SdcOut[i];
    event.v0SdcOut[i]     = src.v0SdcOut[i];
    HF1( 1202, event.nhSdcOut[i] );
    HF1( 1203, event.chisqrSdcOut[i] );
    HF1( 1204, event.x0SdcOut[i] );
    HF1( 1205, event.y0SdcOut[i] );
    HF1( 1206, event.u0SdcOut[i] );
    HF1( 1207, event.v0SdcOut[i] );
    HF2( 1208, event.x0SdcOut[i], event.u0SdcOut[i] );
    HF2( 1209, event.y0SdcOut[i], event.v0SdcOut[i] );
    HF2( 1210, event.x0SdcOut[i], event.y0SdcOut[i] );
  }

  ////////// Kaon
  std::vector< std::vector<double> > KaonDeCont;
  for( int i=0; i<ntKurama; ++i ){
    std::vector<double> dECont;
    for( int j=0; j<NumOfLayersSsdIn; ++j ){
      double de = src.deKaon[j][i];
      event.deKaon[j][i] = de;
      HF1( 10000+j+1, de );
      if( de>0. ) dECont.push_back( de );
    }
    KaonDeCont.push_back( dECont );
  }

  HF1( 1, 9. );

  ////////// Kaon Minus
  std::vector <ThreeVector> KnPCont, KnXCont;
  HF1( 2001, double(ntK18) );
  for( int itK18=0; itK18<ntK18; ++itK18 ){
    double nh     = src.nhK18[itK18];
    double chisqr = src.chisqrK18[itK18];
    double p = src.pK18[itK18];
    double x = src.xtgtK18[itK18], y = src.ytgtK18[itK18];
    double u = src.utgtK18[itK18], v = src.vtgtK18[itK18];
    event.nhK18[itK18]     = nh;
    event.chisqrK18[itK18] = chisqr;
    event.pK18[itK18]      = p;
    event.xtgtK18[itK18] = x; event.ytgtK18[itK18] = y;
    event.utgtK18[itK18] = u; event.vtgtK18[itK18] = v;
    double pt = p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    HF1( 2002, nh ); HF1( 2003, chisqr );
    HF1( 2004, x ); HF1( 2005, y ); HF1( 2006, u ); HF1( 2007, v );
    HF2( 2008, x, u ); HF2( 2009, y, v ); HF2( 2010, x, y );
    HF1( 2011, p );
    if( chisqr>MaxChisqrK18 ) continue;
    HF1( 2102, nh ); HF1( 2103, chisqr );
    HF1( 2104, x ); HF1( 2105, y ); HF1( 2106, u ); HF1( 2107, v );
    HF2( 2108, x, u ); HF2( 2109, y, v ); HF2( 2110, x, y );
    HF1( 2111, p );
    KnPCont.push_back(Mom); KnXCont.push_back(Pos);
  }

  HF1( 2101, double(KnPCont.size()) );
  if( KnPCont.size()==0 ) return false;

  HF1( 1, 10. );

  ////////// Kaon Plus
  std::vector<double> chiKpCont, m2Cont;
  std::vector<ThreeVector> KpPCont, KpXCont;
  for( int itKurama=0; itKurama<ntKurama; ++itKurama ){
    int nh= src.nhKurama[itKurama];
    double chisqr = src.chisqrKurama[itKurama];
    double p = src.pKurama[itKurama];
    double q = src.qKurama[itKurama];
    double x = src.xtgtKurama[itKurama];
    double y = src.ytgtKurama[itKurama];
    double u = src.utgtKurama[itKurama];
    double v = src.vtgtKurama[itKurama];
    double path = src.path[itKurama];
    double stof = src.stof[itKurama];
    double theta = src.thetaKurama[itKurama];
    double pt = p/std::sqrt(1.+u*u+v*v);
    double m2 = src.m2[itKurama];
    ThreeVector PosCorr( x, y, 0. );
    ThreeVector MomCorr( pt*u, pt*v, pt );
    double mean_de = math::Average( KaonDeCont[itKurama] );
    event.chisqrKurama[itKurama] = chisqr;
    event.pKurama[itKurama] = p;
    event.qKurama[itKurama] = q;
    event.xtgtKurama[itKurama] = x;
    event.ytgtKurama[itKurama] = y;
    event.utgtKurama[itKurama] = u;
    event.vtgtKurama[itKurama] = v;
    event.thetaKurama[itKurama] = theta;
    event.stof[itKurama] = stof;
    event.path[itKurama] = path;
    event.m2[itKurama] = m2;
    HF1( 3002, double(nh) ); HF1( 3003, chisqr );
    HF1( 3004, x ); HF1( 3005, y ); HF1( 3006, u ); HF1( 3007, v );
    HF2( 3008, x, u ); HF2( 3009, y, v ); HF2( 3010, x, y );
    HF1( 3011, p ); HF1( 3012, path ); HF1( 3013, stof );
    HF1( 3014, m2 ); HF2( 3015, m2, p );
    HF1( 3016, mean_de ); HF2( 3017, mean_de, p );
    if( q<0 ) continue;
    HF1( 3102, double(nh) ); HF1( 3103, chisqr );
    HF1( 3104, x ); HF1( 3105, y ); HF1( 3106, u ); HF1( 3107, v );
    HF2( 3108, x, u ); HF2( 3109, y, v ); HF2( 3110, x, y );
    HF1( 3111, p ); HF1( 3112, path ); HF1( 3113, stof );
    HF1( 3114, m2 ); HF2( 3115, m2, p );
    HF1( 3116, mean_de ); HF2( 3117, mean_de, p );
    if( chisqr>MaxChisqrKurama ) continue;
    HF1( 3202, double(nh) ); HF1( 3203, chisqr );
    HF1( 3204, x ); HF1( 3205, y ); HF1( 3206, u ); HF1( 3207, v );
    HF2( 3208, x, u ); HF2( 3209, y, v ); HF2( 3210, x, y );
    HF1( 3211, p ); HF1( 3212, path ); HF1( 3213, stof );
    HF1( 3214, m2 ); HF2( 3215, m2, p );
    HF1( 3216, mean_de ); HF2( 3217, mean_de, p );
    if( p<MinKuramaP || MaxKuramaP<p ) continue;
    HF1( 3302, double(nh) ); HF1( 3303, chisqr );
    HF1( 3304, x ); HF1( 3305, y ); HF1( 3306, u ); HF1( 3307, v );
    HF2( 3308, x, u ); HF2( 3309, y, v ); HF2( 3310, x, y );
    HF1( 3311, p ); HF1( 3312, path ); HF1( 3313, stof );
    HF1( 3314, m2 ); HF2( 3315, m2, p );
    HF1( 3316, mean_de ); HF2( 3317, mean_de, p );
    if( m2<MinMassSquare || MaxMassSquare<m2 ) continue;
    HF1( 3402, double(nh) ); HF1( 3403, chisqr );
    HF1( 3404, x ); HF1( 3405, y ); HF1( 3406, u ); HF1( 3407, v );
    HF2( 3408, x, u ); HF2( 3409, y, v ); HF2( 3410, x, y );
    HF1( 3411, p ); HF1( 3412, path ); HF1( 3413, stof );
    HF1( 3414, m2 ); HF2( 3415, m2, p );
    HF1( 3416, mean_de ); HF2( 3417, mean_de, p );
    KpXCont.push_back( PosCorr );
    KpPCont.push_back( MomCorr );
    chiKpCont.push_back( chisqr);
    m2Cont.push_back( m2 );
  }

  if( KpPCont.size()==0 ) return false;

  HF1( 1, 11. );

  //MissingMass
  std::vector<double> MissMassCont;
  std::vector<double> CchiKpCont;
  std::vector<double> Cm2Cont;
  std::vector<ThreeVector> MissPCont, MissXCont;
  std::vector<ThreeVector> CKpPCont, CKpXCont; // good for Xi
  std::vector<ThreeVector> CKnPCont, CKnXCont; // good for Xi
  const int nKn = KnPCont.size();
  const int nKp = KpPCont.size();
  const int nKK = nKn*nKp;
  event.nKn   = nKn;
  event.nKp   = nKp;
  event.nKK   = nKK;
  HF1( 10, double(nKn) );
  HF1( 20, double(nKp) );
  HF1( 30, double(nKK) );
  int nkk = 0;
  for( int ikp=0; ikp<nKp; ++ikp ){
    if( nkk==MaxHits ) break;
    ThreeVector pkp = KpPCont[ikp], xkp = KpXCont[ikp];
    double chisqr = chiKpCont[ikp];
    double m2     = m2Cont[ikp];
    for( int ikn=0; ikn<nKn; ++ikn ){
      if( nkk==MaxHits ) break;
      ThreeVector pkn  = KnPCont[ikn], xkn = KnXCont[ikn];
      ThreeVector vert = Kinematics::VertexPoint( xkn, xkp, pkn, pkp );
      double close = Kinematics::closeDist( xkn, xkp, pkn, pkp );
      double us = pkp.x()/pkp.z(), vs = pkp.y()/pkp.z();
      double ub = pkn.x()/pkn.z(), vb = pkn.y()/pkn.z();
      double cost = pkn*pkp/(pkn.Mag()*pkp.Mag());
      double theta = std::acos(cost)*math::Rad2Deg();
      double pk0  = pkp.Mag();
      double tkn = std::sqrt( KaonMass*KaonMass+pkn.Mag2() );
      double tkp = std::sqrt( KaonMass*KaonMass+pkp.Mag2() );
      LorentzVector LvKn( pkn, tkn );
      LorentzVector LvKp( pkp, tkp );
      LorentzVector LvP( 0., 0., 0., ProtonMass );
      LorentzVector LvC( 0., 0., 0., CarbonMass );
      LorentzVector LvRp = LvKn+LvP-LvKp;
      LorentzVector LvRc = LvKn+LvC-LvKp;
      ThreeVector MissP = LvRp.Vect();
      double MisPp = std::sqrt( LvRp.E()*LvRp.E() - LvRp.Mag2() );
      double MisPu = LvRp.Px()/LvRp.Pz();
      double MisPv = LvRp.Py()/LvRp.Pz();
      double MisMass = LvRp.Mag();
      double MisMassC = LvRc.Mag();
      //Primary frame
      LorentzVector PrimaryLv     = LvKn+LvC;
      double        TotalEnergyCM = PrimaryLv.Mag();
      ThreeVector   beta          = 1./PrimaryLv.E()*PrimaryLv.Vect();
      //CM
      double TotalMomCM
	= 0.5*std::sqrt( ( TotalEnergyCM*TotalEnergyCM
			   -( KaonMass+XiMass )*( KaonMass+XiMass ) ) *
			 ( TotalEnergyCM*TotalEnergyCM
			   -( KaonMass-XiMass )*( KaonMass-XiMass ) ) )
	/ TotalEnergyCM;
      double costLab = cost;
      double cottLab = costLab/std::sqrt(1.-costLab*costLab);
      double bt = beta.Mag(), gamma=1./std::sqrt(1.-bt*bt);
      double gbep=gamma*bt*std::sqrt(TotalMomCM*TotalMomCM+KaonMass*KaonMass)/TotalMomCM;
      double a  = gamma*gamma+cottLab*cottLab;
      double bp = gamma*gbep;
      double c  = gbep*gbep-cottLab*cottLab;
      double dd = bp*bp-a*c;
      if( dd<0. ){
	std::cerr << "dd<0." << std::endl;
	dd = 0.;
      }
      double costCM = (std::sqrt(dd)-bp)/a;
      if( costCM>1. || costCM<-1. ){
	std::cerr << "costCM>1. || costCM<-1." << std::endl;
	costCM=-1.;
      }
      double sintCM  = std::sqrt(1.-costCM*costCM);
      double KaonMom = TotalMomCM*sintCM/std::sqrt(1.-costLab*costLab);
      HF1( 5001, vert.x() ); HF1( 5002, vert.y() );
      HF1( 5003, vert.z() ); HF1( 5004, close );
      HF1( 5005, MisMass );  HF1( 5006, theta );
      HF2( 5007, theta, pkp.Mag() );
      HF1( 5008, MisMassC ); HF2( 5009, MisMass, MisMassC );
      HF2( 5011, MisMass, us ); HF2( 5012, MisMass, vs );
      HF2( 5013, MisMass, ub ); HF2( 5014, MisMass, vb );

      // Vertex Cut
      if( vert.x()<MinVertexX || MaxVertexX<vert.x() ) continue;
      if( vert.y()<MinVertexY || MaxVertexY<vert.y() ) continue;
      if( vert.z()<MinVertexZ || MaxVertexZ<vert.z() ) continue;
      if( MaxCloseDist<close ) continue;

      HF1( 5101, vert.x() ); HF1( 5102, vert.y() );
      HF1( 5103, vert.z() ); HF1( 5104, close );
      HF1( 5105, MisMass );  HF1( 5106, theta );
      HF2( 5107, theta, pkp.Mag() );
      HF1( 5108, MisMassC ); HF2( 5109, MisMass, MisMassC );
      HF2( 5111, MisMass, us ); HF2( 5112, MisMass, vs );
      HF2( 5113, MisMass, ub ); HF2( 5114, MisMass, vb );
      // for Kaon Plus
      HF1( 3503, chisqr );
      HF1( 3504, xkp.x() ); HF1( 3505, xkp.y() );
      HF1( 3506, us ); HF1( 3507, vs );
      HF2( 3508, xkp.x(), us ); HF2( 3509, xkp.y(), vs );
      HF2( 3510, xkp.x(), xkp.y() );
      HF1( 3511, pkp.Mag() ); HF1( 3514, m2 ); HF2( 3515, m2, pkp.Mag() );

      event.vtx[nkk] = vert.x(); event.vty[nkk] = vert.y();
      event.vtz[nkk] = vert.z();
      event.closeDist[nkk] = close;
      event.theta[nkk]     = theta;
      event.thetaCM[nkk]   = std::acos(costCM)*math::Rad2Deg();
      event.costCM[nkk]   = costCM;
      event.MissMass[nkk] = MisMass;
      event.MissP[nkk]  = MisPp;
      event.MissPu[nkk] = MisPu; event.MissPv[nkk] = MisPv;
      event.xkp[nkk] = xkp.x(); event.ykp[nkk] = xkp.y();
      event.ukp[nkk] = us;      event.vkp[nkk] = vs;
      event.xkn[nkk] = xkn.x(); event.ykn[nkk] = xkn.y();
      event.ukn[nkk] = ub;      event.vkn[nkk] = vb;
      event.pOrg[nkk]  = pk0;
      event.pCalc[nkk] = KaonMom;
      MissMassCont.push_back( MisMass );
      MissPCont.push_back( MissP );
      MissXCont.push_back( vert  );
      CKpPCont.push_back( pkp ); CKpXCont.push_back( xkp );
      CKnPCont.push_back( pkn ); CKnXCont.push_back( xkn );
      CchiKpCont.push_back( chisqr );
      Cm2Cont.push_back( m2 );
      nkk++;
    }
  }

  HF1( 31, double(nkk) );
  if( MissPCont.size()==0 ) return false;

  HF1( 1, 12. );

  HF1( 121, nhBh2 );
  for( int i=0; i<nhBh2; ++i ){
    HF1( 122, event.csBh2[i] );  HF1( 123, event.Bh2Seg[i] );
    HF1( 124, event.tBh2[i]  );  HF1( 125, event.deBh2[i]  );
    HF1( 126, event.t0Bh2[i] );
  }
  HF1( 221, nhBh1 ); HF1( 226, btof );
  for( int i=0; i<nhBh1; ++i ){
    HF1( 222, event.csBh1[i] ); HF1( 223, event.Bh1Seg[i] );
    HF1( 224, event.tBh1[i]  ); HF1( 225, event.deBh1[i]  );
  }
  HF1( 321, nhTof );
  for( int i=0; i<nhTof; ++i ){
    HF1( 322, event.csTof[i] ); HF1( 323, event.TofSeg[i] );
    HF1( 324, event.tTof[i]  ); HF1( 325, event.deTof[i]  );
  }

  ////////// Xi
  std::vector<ThreeVector> XiPCont, XiXCont;
  std::vector< std::vector<double> > XiDeCont;
  HF1( 50, (double)ntXi );
  for( int itXi=0; itXi<ntXi; ++itXi ){
    double x0 = src.x0Xi[itXi];
    double y0 = src.y0Xi[itXi];
    double u0 = src.u0Xi[itXi];
    double v0 = src.v0Xi[itXi];
    std::vector<double> DeCont;
    for( int j=0; j<NumOfLayersSsdIn; ++j ){
      double de = src.deXi[j][itXi];
      DeCont.push_back( de );
    }
    double xtgt = x0 + u0*zTarget;
    double ytgt = y0 + v0*zTarget;
    double pt   = 1./std::sqrt(1.+u0*u0+v0*v0);
    XiPCont.push_back( ThreeVector( pt*u0, pt*v0, pt ) );
    XiXCont.push_back( ThreeVector( xtgt,  ytgt,  0. ) );
    XiDeCont.push_back( DeCont );
  }

  HF1( 51, (double)XiPCont.size() );
  if( XiPCont.size()==0 ) return false;

  HF1( 1, 13. );

  // Xi Analysis
  std::vector<double>      FinalChiKpCont;
  std::vector<double>      FinalM2Cont;
  std::vector<double>      FinalMissMassCont;
  std::vector<ThreeVector> FinalKKVertexCont;
  std::vector<ThreeVector> FinalKXiVertexCont;
  std::vector<ThreeVector> FinalKnPCont;
  std::vector<ThreeVector> FinalKpPCont;
  std::vector<ThreeVector> FinalMMPCont;
  std::vector<ThreeVector> FinalXiPCont;
  std::vector<ThreeVector> FinalXiXCont;
  std::vector< std::vector<double> > FinalSSDdECont;
  const int nMM = MissPCont.size();
  const int nXi = XiPCont.size();
  event.nMM  = nMM;
  event.nXi  = nXi;
  HF1( 40, double(nMM) );
  HF1( 52, double(nXi) );
  int nkxi = 0;
  for( int imm=0; imm<nMM; ++imm ){
    if( nkxi==MaxHits ) break;
    ThreeVector pmm = MissPCont[imm], xmm = MissXCont[imm];
    ThreeVector pkp = CKpPCont[imm],   xkp = CKpXCont[imm];
    ThreeVector pkn = CKnPCont[imm],   xkn = CKnXCont[imm];
    double MissMass = MissMassCont[imm];
    double chisqr   = CchiKpCont[imm];
    double m2       = Cm2Cont[imm];
    for( int ixi=0; ixi<nXi; ++ixi ){
      if( nkxi==MaxHits ) break;
      double mean_de = math::Average( XiDeCont[ixi] );
      ThreeVector pxi    = XiPCont[ixi], xxi = XiXCont[ixi];
      ThreeVector vertXi = Kinematics::VertexPoint( xxi, xkp, pxi, pkp );
      double closeXi = Kinematics::closeDist( xxi, xkp, pxi, pkp );
      double vtxXi = vertXi.x();
      double vtyXi = vertXi.y();
      double vtzXi = vertXi.z();
      double costKpXi  = pkp*pxi/(pkp.Mag()*pxi.Mag());
      double thetaKpXi = std::acos(costKpXi)*math::Rad2Deg();
      double costKpMM  = pkp*pmm/(pkp.Mag()*pmm.Mag());
      double thetaKpMM = std::acos(costKpMM)*math::Rad2Deg();
      double costKnXi  = pkn*pxi/(pkn.Mag()*pxi.Mag());
      double thetaKnXi = std::acos(costKnXi)*math::Rad2Deg();
      double costKnMM  = pkn*pmm/(pkn.Mag()*pmm.Mag());
      double thetaKnMM = std::acos(costKnMM)*math::Rad2Deg();
      double costKnKp  = pkn*pkp/(pkn.Mag()*pkp.Mag());
      double thetaKnKp = std::acos(costKnKp)*math::Rad2Deg();
      double umm = pmm.x()/pmm.z();
      double vmm = pmm.y()/pmm.z();
      double uxi = pxi.x()/pxi.z();
      double vxi = pxi.y()/pxi.z();
      // if( std::abs( uxi ) > 1.0 ) continue;
      // if( std::abs( vxi ) > 1.0 ) continue;
      double resu = uxi - umm;
      double resv = vxi - vmm;
      // x^2/a^2 + y^2/b^2 < 1
      double x2a2 = (resu*resu)/(ResMissU*ResMissU);
      double y2b2 = (resv*resv)/(ResMissV*ResMissV);

      HF1( 6001, vtxXi ); HF1( 6002, vtyXi ); HF1( 6003, vtzXi );
      HF1( 6004, closeXi ); HF1( 6005, MissMass ); HF1( 6006, mean_de );
      HF1( 6007, pmm.Mag() ); HF2( 6008, mean_de, pmm.Mag() );
      HF2( 6011, MissMass, umm ); HF2( 6012, MissMass, vmm );
      HF2( 6013, MissMass, uxi ); HF2( 6014, MissMass, vxi );
      HF2( 6021, uxi, umm );  HF2( 6022, vxi, vmm );
      HF1( 6023, uxi - umm ); HF1( 6024, vxi - vmm );
      HF2( 6025, uxi - umm, vxi - vmm );
      HF1( 6026, x2a2 + y2b2 );
      HF1( 6027, thetaKpXi ); HF1( 6028, thetaKpMM );
      HF1( 6029, thetaKnXi ); HF1( 6030, thetaKnMM ); HF1( 6031, thetaKnKp );
      HF2( 6032, thetaKpXi, thetaKpMM ); HF1( 6033, thetaKpXi - thetaKpMM );
      HF2( 6034, thetaKnXi, thetaKnMM ); HF1( 6035, thetaKnXi - thetaKnMM );
      HF2( 6041, vtxXi, xmm.x() ); HF2( 6042, vtyXi, xmm.y() );
      HF2( 6043, vtzXi, xmm.z() ); HF1( 6044, vtxXi - xmm.x() );
      HF1( 6045, vtyXi - xmm.y() ); HF1( 6046, vtzXi - xmm.z() );
      for( int j=0; j<NumOfLayersSsdIn; ++j ){
	HF2( 6051+j, XiDeCont[ixi][j], pmm.Mag() );
      }

      // SsdDeCut
      double min_de = *std::min_element( XiDeCont[ixi].begin(),
					 XiDeCont[ixi].end() );
      if( min_de<MinDeSSDXi ) continue;

      HF1( 6101, vtxXi ); HF1( 6102, vtyXi ); HF1( 6103, vtzXi );
      HF1( 6104, closeXi ); HF1( 6105, MissMass ); HF1( 6106, mean_de );
      HF1( 6107, pmm.Mag() ); HF2( 6108, mean_de, pmm.Mag() );
      HF2( 6111, MissMass, umm ); HF2( 6112, MissMass, vmm );
      HF2( 6113, MissMass, uxi ); HF2( 6114, MissMass, vxi );
      HF2( 6121, uxi, umm );  HF2( 6122, vxi, vmm );
      HF1( 6123, uxi - umm ); HF1( 6124, vxi - vmm );
      HF2( 6125, uxi - umm, vxi - vmm );
      HF1( 6126, x2a2 + y2b2 );
      HF1( 6127, thetaKpXi ); HF1( 6128, thetaKpMM );
      HF1( 6129, thetaKnXi ); HF1( 6130, thetaKnMM ); HF1( 6131, thetaKnKp );
      HF2( 6132, thetaKpXi, thetaKpMM ); HF1( 6133, thetaKpXi - thetaKpMM );
      HF2( 6134, thetaKnXi, thetaKnMM ); HF1( 6135, thetaKnXi - thetaKnMM );
      HF2( 6141, vtxXi, xmm.x() ); HF2( 6142, vtyXi, xmm.y() );
      HF2( 6143, vtzXi, xmm.z() ); HF1( 6144, vtxXi - xmm.x() );
      HF1( 6145, vtyXi - xmm.y() ); HF1( 6146, vtzXi - xmm.z() );
      for( int j=0; j<NumOfLayersSsdIn; ++j ){
	HF2( 6151+j, XiDeCont[ixi][j], pmm.Mag() );
      }

      // Vertex Cut
      if( vertXi.x()<MinVertexXiX || MaxVertexXiX<vertXi.x() ) continue;
      if( vertXi.y()<MinVertexXiY || MaxVertexXiY<vertXi.y() ) continue;
      if( vertXi.z()<MinVertexXiZ || MaxVertexXiZ<vertXi.z() ) continue;
      if( MaxCloseDistXi<closeXi ) continue;

      HF1( 6201, vtxXi ); HF1( 6202, vtyXi ); HF1( 6203, vtzXi );
      HF1( 6204, closeXi ); HF1( 6205, MissMass ); HF1( 6206, mean_de );
      HF1( 6207, pmm.Mag() ); HF2( 6208, mean_de, pmm.Mag() );
      HF2( 6211, MissMass, umm ); HF2( 6212, MissMass, vmm );
      HF2( 6213, MissMass, uxi ); HF2( 6214, MissMass, vxi );
      HF2( 6221, uxi, umm );  HF2( 6222, vxi, vmm );
      HF1( 6223, uxi - umm ); HF1( 6224, vxi - vmm );
      HF2( 6225, uxi - umm, vxi - vmm );
      HF1( 6226, x2a2 + y2b2 );
      HF1( 6227, thetaKpXi ); HF1( 6228, thetaKpMM );
      HF1( 6229, thetaKnXi ); HF1( 6230, thetaKnMM ); HF1( 6231, thetaKnKp );
      HF2( 6232, thetaKpXi, thetaKpMM ); HF1( 6233, thetaKpXi - thetaKpMM );
      HF2( 6234, thetaKnXi, thetaKnMM ); HF1( 6235, thetaKnXi - thetaKnMM );
      HF2( 6241, vtxXi, xmm.x() ); HF2( 6242, vtyXi, xmm.y() );
      HF2( 6243, vtzXi, xmm.z() ); HF1( 6244, vtxXi - xmm.x() );
      HF1( 6245, vtyXi - xmm.y() ); HF1( 6246, vtzXi - xmm.z() );
      for( int j=0; j<NumOfLayersSsdIn; ++j ){
	HF2( 6251+j, XiDeCont[ixi][j], pmm.Mag() );
      }

      // Angle Cut w/3-sigma
      if( x2a2 + y2b2 > 3.*3. ) continue;

      HF1( 6301, vtxXi ); HF1( 6302, vtyXi ); HF1( 6303, vtzXi );
      HF1( 6304, closeXi ); HF1( 6305, MissMass ); HF1( 6306, mean_de );
      HF1( 6307, pmm.Mag() ); HF2( 6308, mean_de, pmm.Mag() );
      HF2( 6311, MissMass, umm ); HF2( 6312, MissMass, vmm );
      HF2( 6313, MissMass, uxi ); HF2( 6314, MissMass, vxi );
      HF2( 6321, uxi, umm );  HF2( 6322, vxi, vmm );
      HF1( 6323, uxi - umm ); HF1( 6324, vxi - vmm );
      HF2( 6325, uxi - umm, vxi - vmm );
      HF1( 6326, x2a2 + y2b2 );
      HF1( 6327, thetaKpXi ); HF1( 6328, thetaKpMM );
      HF1( 6329, thetaKnXi ); HF1( 6330, thetaKnMM ); HF1( 6331, thetaKnKp );
      HF2( 6332, thetaKpXi, thetaKpMM ); HF1( 6333, thetaKpXi - thetaKpMM );
      HF2( 6334, thetaKnXi, thetaKnMM ); HF1( 6335, thetaKnXi - thetaKnMM );
      HF2( 6341, vtxXi, xmm.x() ); HF2( 6342, vtyXi, xmm.y() );
      HF2( 6343, vtzXi, xmm.z() ); HF1( 6344, vtxXi - xmm.x() );
      HF1( 6345, vtyXi - xmm.y() ); HF1( 6346, vtzXi - xmm.z() );
      for( int j=0; j<NumOfLayersSsdIn; ++j ){
	HF2( 6351+j, XiDeCont[ixi][j], pmm.Mag() );
      }

      // Vertex Residual Cut
      if( std::abs( vtxXi - xmm.x() )>ResVertexX ) continue;
      if( std::abs( vtyXi - xmm.y() )>ResVertexY ) continue;
      if( std::abs( vtzXi - xmm.z() )>ResVertexZ ) continue;

      HF1( 6401, vtxXi ); HF1( 6402, vtyXi ); HF1( 6403, vtzXi );
      HF1( 6404, closeXi ); HF1( 6405, MissMass ); HF1( 6406, mean_de );
      HF1( 6407, pmm.Mag() ); HF2( 6408, mean_de, pmm.Mag() );
      HF2( 6411, MissMass, umm ); HF2( 6412, MissMass, vmm );
      HF2( 6413, MissMass, uxi ); HF2( 6414, MissMass, vxi );
      HF2( 6421, uxi, umm );  HF2( 6422, vxi, vmm );
      HF1( 6423, uxi - umm ); HF1( 6424, vxi - vmm );
      HF2( 6425, uxi - umm, vxi - vmm );
      HF1( 6426, x2a2 + y2b2 );
      HF1( 6427, thetaKpXi ); HF1( 6428, thetaKpMM );
      HF1( 6429, thetaKnXi ); HF1( 6430, thetaKnMM ); HF1( 6431, thetaKnKp );
      HF2( 6432, thetaKpXi, thetaKpMM ); HF1( 6433, thetaKpXi - thetaKpMM );
      HF2( 6434, thetaKnXi, thetaKnMM ); HF1( 6435, thetaKnXi - thetaKnMM );
      HF2( 6441, vtxXi, xmm.x() ); HF2( 6442, vtyXi, xmm.y() );
      HF2( 6443, vtzXi, xmm.z() ); HF1( 6444, vtxXi - xmm.x() );
      HF1( 6445, vtyXi - xmm.y() ); HF1( 6446, vtzXi - xmm.z() );
      for( int j=0; j<NumOfLayersSsdIn; ++j ){
	HF2( 6451+j, XiDeCont[ixi][j], pmm.Mag() );
	HF1( 10011+j, XiDeCont[ixi][j] );
      }

      event.x0Xi[nkxi]    = xxi.x();
      event.y0Xi[nkxi]    = xxi.y();
      event.u0Xi[nkxi]    = uxi;
      event.v0Xi[nkxi]    = vxi;
      event.thetaXi[nkxi] = thetaKpXi;
      event.vtxXi[nkxi]   = vtxXi;
      event.vtyXi[nkxi]   = vtyXi;
      event.vtzXi[nkxi]   = vtzXi;
      event.closeXi[nkxi] = closeXi;
      // Kaon Plus
      HF1( 3603, chisqr );
      HF1( 3604, xkp.x() ); HF1( 3605, xkp.y() );
      HF1( 3606, pkp.x()/pkp.z() ); HF1( 3607, pkp.y()/pkp.z() );
      HF2( 3608, xkp.x(), pkp.x()/pkp.z() );
      HF2( 3609, xkp.y(), pkp.y()/pkp.z() );
      HF2( 3610, xkp.x(), xkp.y() );
      HF1( 3611, pkp.Mag() );
      HF1( 3614, m2 ); HF2( 3615, m2, pkp.Mag() );

      ///// Xi Track on Emulsion
      double xemul = xxi.x() + uxi*( zEmulsion - zTarget );
      double yemul = xxi.y() + vxi*( zEmulsion - zTarget );
      double xemc  = event.xpos + xemul;
      double yemc  = event.ypos + yemul;
      // on Target
      HF1( 7001, xxi.x() ); HF1( 7002, xxi.y() ); HF1( 7003, uxi ); HF1( 7004, vxi );
      HF2( 7005, xxi.x(), uxi ); HF2( 7006, xxi.y(), vxi ); HF2( 7007, xxi.x(), xxi.y() );
      // on Emulsion
      HF1( 7101, xemul ); HF1( 7102, yemul ); HF1( 7103, uxi ); HF1( 7104, vxi );
      HF2( 7105, xemul, uxi ); HF2( 7106, yemul, vxi ); HF2( 7107, xemul, yemul );
      // on Emulsion w/EM coordinate
      HF1( 7201, xemc ); HF1( 7202, yemc ); HF1( 7203, uxi ); HF1( 7204, vxi );
      HF2( 7205, xemc, uxi ); HF2( 7206, yemc, vxi ); HF2( 7207, xemc, yemc );

      HF1( 2, event.runnum );
      event.x0Xi[nkxi] = xxi.x(); event.y0Xi[nkxi] = xxi.y();
      event.u0Xi[nkxi] = uxi; event.v0Xi[nkxi] = vxi;
      // push back final values
      FinalChiKpCont.push_back( chisqr );
      FinalM2Cont.push_back( m2 );
      FinalMissMassCont.push_back( MissMass );
      FinalKKVertexCont.push_back( xmm );
      FinalKXiVertexCont.push_back( vertXi );
      FinalKnPCont.push_back( pkn );
      FinalKpPCont.push_back( pkp );
      FinalMMPCont.push_back( pmm );
      FinalXiPCont.push_back( pxi );
      FinalXiXCont.push_back( xxi );
      FinalSSDdECont.push_back( XiDeCont[ixi] );
      nkxi++;
    }
  }

  event.nKXi = nkxi;//nMM*nXi;
  event.ntXi = nkxi;

  if( nkxi==0 ) return false;

  HF1( 53, double(nkxi) );

  HF1( 1, 14. );

  HF1( 11, event.xpos ); HF1( 12, event.ypos );
  HF2( 13, event.xpos, event.ypos ); HF1( 14, event.state );

  //Final Hodoscope histograms
  HF1( 131, nhBh2 );
  for( int i=0; i<nhBh2; ++i ){
    HF1( 132, event.csBh2[i] );  HF1( 133, event.Bh2Seg[i] );
    HF1( 134, event.tBh2[i]  );  HF1( 135, event.deBh2[i]  );
    HF1( 136, event.t0Bh2[i] );
  }
  HF1( 231, nhBh1 ); HF1( 236, btof );
  for( int i=0; i<nhBh1; ++i ){
    HF1( 232, event.csBh1[i] ); HF1( 233, event.Bh1Seg[i] );
    HF1( 234, event.tBh1[i]  ); HF1( 235, event.deBh1[i]  );
  }
  HF1( 331, nhTof );
  for( int i=0; i<nhTof; ++i ){
    HF1( 332, event.csTof[i] ); HF1( 333, event.TofSeg[i] );
    HF1( 334, event.tTof[i]  ); HF1( 335, event.deTof[i]  );
  }

  HF1( 1, 19. );

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
  HB1( 1, "Status", 20, 0., 20. );
  HB1( 2, "Run Number", 2000, 0., 2000. );
  HB1( 11, "EMC Xpos", 200, -200., 200. );
  HB1( 12, "EMC Ypos", 200, -200., 200. );
  HB2( 13, "EMC Xpos%Ypos", 200, -200., 200., 200, -200., 200. );
  HB1( 14, "EMC State", 10, 0., 10. );

  HB1( 10, "NKn", 20, 0., 20. );
  HB1( 20, "NKp", 20, 0., 20. );
  HB1( 30, "NKK", 20, 0., 20. );
  HB1( 31, "NKK [w/VertexCut]", 20, 0., 20. );
  HB1( 40, "nMM", 20, 0., 20. );
  HB1( 50, "nXi", 20, 0., 20. );
  HB1( 51, "nXi [w/DeCut]", 20, 0., 20. );
  HB1( 52, "nXi [w/DeCut]", 20, 0., 20. );
  HB1( 53, "nXi [Final]", 20, 0., 20. );

  { // BH2
    std::vector<TString> s{ "", " [Hodo]", " [KK]", " [Xi]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 101+i*10, "#Clusters BH2"+s[i], 7, 0., 7. );
      HB1( 102+i*10, "ClusterSize BH2"+s[i], 7, 0., 7. );
      HB1( 103+i*10, "HitPat BH2"+s[i], 10, 0., 10. );
      HB1( 104+i*10, "MeanTime BH2"+s[i], 400, -10., 10. );
      HB1( 105+i*10, "Delta-E BH2"+s[i], 200, -0.5, 4.5 );
      HB1( 106+i*10, "Time0"+s[i], 400, -10., 10. );
    }
  }

  { // BH1
    std::vector<TString> s{ "", " [Hodo]", " [KK]", " [Xi]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 201+i*10, "#Clusters BH1"+s[i], 11, 0., 11. );
      HB1( 202+i*10, "ClusterSize BH1"+s[i], 11, 0., 11. );
      HB1( 203+i*10, "HitPat BH1"+s[i], 11, 0., 11. );
      HB1( 204+i*10, "MeanTime BH1"+s[i], 400, -10., 10. );
      HB1( 205+i*10, "Delta-E BH1"+s[i], 200, -0.5, 4.5 );
      HB1( 206+i*10, "Beam ToF"+s[i], 400, -10., 10. );
    }
  }

  { // TOF
    std::vector<TString> s{ "", " [Hodo]", " [KK]", " [Xi]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 301+i*10, "#Clusters TOF"+s[i], 32, 0., 32. );
      HB1( 302+i*10, "ClusterSize TOF"+s[i], 32, 0., 32. );
      HB1( 303+i*10, "HitPat TOf"+s[i], 32, 0., 32. );
      HB1( 304+i*10, "TimeOfFlight TOF"+s[i], 500, -50., 100. );
      HB1( 305+i*10, "Delta-E TOF"+s[i], 200, -0.5, 4.5 );
    }
  }

  HB1( 1001, "#Tracks BcOut", 10, 0., 10. );
  HB1( 1002, "#Hits BcOut", 20, 0., 20. );
  HB1( 1003, "Chisqr BcOut", 200, 0., 100. );
  HB1( 1004, "Xtgt BcOut", 500, -200., 200. );
  HB1( 1005, "Ytgt BcOut", 500, -200., 200. );
  HB1( 1006, "Utgt BcOut",  700, -0.35, 0.35 );
  HB1( 1007, "Vtgt BcOut",  400, -0.20, 0.20 );
  HB2( 1008, "Xtgt%U BcOut", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1009, "Ytgt%V BcOut", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1010, "Xtgt%Ytgt BcOut", 100, -200., 200., 100, -200., 200. );

  HB1( 1101, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 1102, "#Hits SdcIn", 20, 0., 20. );
  HB1( 1103, "Chisqr SdcIn", 200, 0., 100. );
  HB1( 1104, "X0 SdcIn", 500, -200., 200. );
  HB1( 1105, "Y0 SdcIn", 500, -200., 200. );
  HB1( 1106, "U0 SdcIn",  700, -0.35, 0.35 );
  HB1( 1107, "V0 SdcIn",  400, -0.20, 0.20 );
  HB2( 1108, "X0%U0 SdcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1109, "Y0%V0 SdcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1110, "X0%Y0 SdcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 1201, "#Tracks SdcOut", 10, 0., 10. );
  HB1( 1202, "#Hits SdcOut", 20, 0., 20. );
  HB1( 1203, "Chisqr SdcOut", 200, 0., 100. );
  HB1( 1204, "X0 SdcOut", 600, -1200., 1200. );
  HB1( 1205, "Y0 SdcOut", 600, -600., 600. );
  HB1( 1206, "U0 SdcOut",  700, -0.35, 0.35 );
  HB1( 1207, "V0 SdcOut",  400, -0.20, 0.20 );
  HB2( 1208, "X0%U0 SdcOut", 120, -1200., 1200., 100, -0.35, 0.35 );
  HB2( 1209, "Y0%V0 SdcOut", 100,  -600.,  600., 100, -0.20, 0.20 );
  HB2( 1210, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -600., 600. );

  { // Kaon Minus (K18Tracking)
    std::vector<TString> s{ "", " [chisqrK18]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 2001+i*100, "NTracks K18"+s[i], 10, 0., 10. );
      HB1( 2002+i*100, "NHits K18"+s[i], 30, 0., 30. );
      HB1( 2003+i*100, "Chisqr K18"+s[i], 500, 0., 50. );
      HB1( 2004+i*100, "Xtgt K18"+s[i], 500, -200., 200. );
      HB1( 2005+i*100, "Ytgt K18"+s[i], 500, -100., 100. );
      HB1( 2006+i*100, "Utgt K18"+s[i], 400, -0.35, 0.35 );
      HB1( 2007+i*100, "Vtgt K18"+s[i], 200, -0.20, 0.20 );
      HB2( 2008+i*100, "Xtgt%U K18"+s[i], 100, -200., 200., 100, -0.35, 0.35 );
      HB2( 2009+i*100, "Ytgt%V K18"+s[i], 100, -100., 100., 100, -0.20, 0.20 );
      HB2( 2010+i*100, "Xtgt%Ytgt K18"+s[i], 100, -200., 200., 100, -100., 100. );
      HB1( 2011+i*100, "Momentum K18"+s[i], 500, 1.5, 2.0 );
      // HB1( 2012+i*100, "Xin  K18"+s[i], 500, -100., 100. );
      // HB1( 2013+i*100, "Yin  K18"+s[i], 500, -100., 100. );
      // HB1( 2014+i*100, "Uin  K18"+s[i], 400, -0.35, 0.35 );
      // HB1( 2015+i*100, "Vin  K18"+s[i], 200, -0.20, 0.20 );
    }
  }

  { // Kaon Plus (KuramaTracking)
    std::vector<TString> s{ "", " [qKurama]", " [chisqrKurama]", " [pKurama]", " [m2]", " [KK]", " [Xi]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 3001+i*100, "#Tracks Kurama"+s[i], 10, 0., 10. );
      HB1( 3002+i*100, "#Hits Kurama"+s[i], 30, 0., 30. );
      HB1( 3003+i*100, "Chisqr Kurama"+s[i], 500, 0., 500. );
      HB1( 3004+i*100, "Xtgt Kurama"+s[i], 500, -200., 200. );
      HB1( 3005+i*100, "Ytgt Kurama"+s[i], 500, -100., 100. );
      HB1( 3006+i*100, "Utgt Kurama"+s[i], 400, -0.35, 0.35 );
      HB1( 3007+i*100, "Vtgt Kurama"+s[i], 200, -0.20, 0.20 );
      HB2( 3008+i*100, "Xtgt%U Kurama"+s[i], 100, -200., 200., 100, -0.35, 0.35 );
      HB2( 3009+i*100, "Ytgt%V Kurama"+s[i], 100, -100., 100., 100, -0.20, 0.20 );
      HB2( 3010+i*100, "Xtgt%Ytgt Kurama"+s[i], 100, -200., 200., 100, -100., 100. );
      HB1( 3011+i*100, "Momentum Kurama"+s[i], 200, 0.0, 3.0 );
      HB1( 3012+i*100, "PathLength Kurama"+s[i], 600, 3000., 4000. );
      HB1( 3013+i*100, "Time Of Flight Kurama"+s[i], 200, 0., 100. );
      HB1( 3014+i*100, "Mass Square"+s[i], 400, -0.2, 1.8 );
      HB2( 3015+i*100, "pKurama%Mass Square"+s[i], 200, -0.2, 1.8, 200, 0., 3. );
      HB1( 3016+i*100, "SSD MeanDeltaE Kaon"+s[i], 200, 0., 100000. );
      HB2( 3017+i*100, "SSD MeanDeltaE Kaon%pKurama"+s[i],
	   100, 0., 100000., 100, 0., 3. );
    }
  }

  { // KK Reaction
    std::vector<TString> s{ ""," [VertexKK]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 5001+i*100, "KK Vertex X"+s[i], 200, -100., 100. );
      HB1( 5002+i*100, "KK Vertex Y"+s[i], 200, -100., 100. );
      HB1( 5003+i*100, "KK Vertex Z"+s[i], 200, -200., 200. );
      HB1( 5004+i*100, "KK Closest Distance"+s[i], 200, 0., 20. );
      HB1( 5005+i*100, "KK MissingMass"+s[i], 800, 1.0, 1.8 );
      HB1( 5006+i*100, "KK Theta"+s[i], 180, 0., 90. );
      HB2( 5007+i*100, "pKurama%Theta"+s[i], 90, 0., 45., 100, 0.75, 1.75 );
      HB1( 5008+i*100, "KK MissingMass w/C"+s[i], 800, 10.5, 12.5 );
      HB2( 5009+i*100, "KK MissingMass w/P%w/C"+s[i], 100, 1.0, 1.8, 100, 10.5, 12.5 );
      HB2( 5011+i*100, "MissingMass%Us"+s[i], 200, 0.0, 2.50, 100, -0.40, 0.40 );
      HB2( 5012+i*100, "MissingMass%Vs"+s[i], 200, 0.0, 2.50, 100, -0.20, 0.20 );
      HB2( 5013+i*100, "MissingMass%Ub"+s[i], 200, 0.0, 2.50, 100, -0.30, 0.30 );
      HB2( 5014+i*100, "MissingMass%Vb"+s[i], 200, 0.0, 2.50, 100, -0.10, 0.10 );
    }
  }

  { // Xi analysis
    std::vector<TString> s{ ""," [dE]", " [VertexXi]", " [Angle]", " [VertexRes]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 6001+i*100, "KXi Vertex X"+s[i], 200, -100., 100. );
      HB1( 6002+i*100, "KXi Vertex Y"+s[i], 200, -100., 100. );
      HB1( 6003+i*100, "KXi Vertex Z"+s[i], 200, -200., 200. );
      HB1( 6004+i*100, "KXi ClosestDistance"+s[i], 200, 0., 20. );
      HB1( 6005+i*100, "MissingMass"+s[i], 800, 1.0, 1.8 );
      HB1( 6006+i*100, "SSD MeanDeltaE"+s[i], 200, 0., 150000. );
      HB1( 6007+i*100, "Missing Momentum"+s[i], 200, 0., 1.6 );
      HB2( 6008+i*100, "SSD MeanDeltaE%Missing Momentum"+s[i],
	   100, 0., 150000., 100, 0., 1.6 );
      HB2( 6011+i*100, "MissingMass%MissU"+s[i], 200, 1., 2., 100, -1., 1. );
      HB2( 6012+i*100, "MissingMass%MissV"+s[i], 200, 1., 2., 100, -1., 1. );
      HB2( 6013+i*100, "MissingMass%XiU"+s[i], 200, 1., 2., 100, -1., 1. );
      HB2( 6014+i*100, "MissingMass%XiV"+s[i], 200, 1., 2., 100, -1., 1. );
      HB2( 6021+i*100, "MissU%XiU"+s[i], 100, -1.0, 1.0, 100, -1.0, 1.0 );
      HB2( 6022+i*100, "MissV%XiV"+s[i], 100, -1.0, 1.0, 100, -1.0, 1.0 );
      HB1( 6023+i*100, "Residual MissU-XiU"+s[i], 100, -1.0, 1.0 );
      HB1( 6024+i*100, "Residual MissV-XiV"+s[i], 100, -1.0, 1.0 );
      HB2( 6025+i*100, "Residual MissV-XiV%MissU-XiU"+s[i], 100, -1., 1., 100, -1., 1. );
      HB1( 6026+i*100, "du^2/su^2+dv^2*/sv^2"+s[i], 200, 0., 10. );
      HB1( 6027+i*100, "Theta KpXi"+s[i], 120, 0., 60. );
      HB1( 6028+i*100, "Theta KpMM"+s[i], 120, 0., 60. );
      HB1( 6029+i*100, "Theta KnXi"+s[i], 120, 0., 60. );
      HB1( 6030+i*100, "Theta KnMM"+s[i], 120, 0., 60. );
      HB1( 6031+i*100, "Theta KnKp"+s[i], 120, 0., 60. );
      HB2( 6032+i*100, "Theta KpXi%KpMM"+s[i], 120, 0., 60., 120, 0., 60. );
      HB1( 6033+i*100, "Residual Theta KpXi-KpMM"+s[i], 120, -30., 30. );
      HB2( 6034+i*100, "Theta KnXi%KnMM"+s[i], 120, 0., 60., 120, 0., 60. );
      HB1( 6035+i*100, "Residual Theta KnXi-KnMM"+s[i], 120, -30., 30. );
      HB2( 6041+i*100, "KK%KXi Vertex X"+s[i], 100, -100., 100., 100, -100., 100. );
      HB2( 6042+i*100, "KK%KXi Vertex Y"+s[i], 100, -60., 60., 100, -60., 60. );
      HB2( 6043+i*100, "KK%KXi Vertex Z"+s[i], 100, -100., 100., 100, -100., 100. );
      HB1( 6044+i*100, "Residual KXi-KK Vertex X"+s[i], 200, -20., 20. );
      HB1( 6045+i*100, "Residual KXi-KK Vertex Y"+s[i], 200, -20., 20. );
      HB1( 6046+i*100, "Residual KXi-KK Vertex Z"+s[i], 200, -40., 40. );
      for( int j=0; j<NumOfLayersSsdIn; ++j ){
	HB2( 6051+i*100+j, Form("MissMom%%deXi Ssd%d", j+1)+s[i],
	     100, 0., 100000., 100, 0., 1.6 );
      }
    }
  }

  { // Final Xi Track; on Target, on Emulsion and EM coordinate
    std::vector<TString> s{ " [Target]"," [Emulsion]", " [EMcoord]" };
    for( int i=0,n=s.size(); i<n; ++i ){
      HB1( 7001+i*100, "Xi X"+s[i], 200, -200., 200. );
      HB1( 7002+i*100, "Xi Y"+s[i], 200, -200., 200. );
      HB1( 7003+i*100, "Xi U"+s[i], 200, -2.0, 2.0 );
      HB1( 7004+i*100, "Xi V"+s[i], 200, -2.0, 2.0 );
      HB2( 7005+i*100, "Xi X%U"+s[i], 100, -200., 200., 100, -2.0, 2.0 );
      HB2( 7006+i*100, "Xi Y%V"+s[i], 100, -200., 200., 100, -2.0, 2.0 );
      HB2( 7007+i*100, "Xi X%Y"+s[i], 100, -250., 250., 100, -250., 250. );
    }
  }

  for( int i=0; i<NumOfLayersSsdIn; ++i ){
    HB1( 10000+i+1, Form("Kaon DeltaE Ssd%d", i+1), 200, 0., 100000. );
    HB1( 10010+i+1, Form("Xi DeltaE Ssd%d", i+1), 200, 0., 100000. );
  }

  ////////////////////////////////////////////
  //Tree
  HBTree("xi","tree of XiAna");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  //Trigger
  tree->Branch("trigpat",   event.trigpat,
	       Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",  event.trigflag,
	       Form("trigflag[%d]/I", NumOfSegTrig));

  //Hodoscope
  tree->Branch("nhBh1",   &event.nhBh1,   "nhBh1/I");
  tree->Branch("csBh1",    event.csBh1,   "csBh1[nhBh1]/I");
  tree->Branch("Bh1Seg",   event.Bh1Seg,  "Bh1Seg[nhBh1]/D");
  tree->Branch("tBh1",     event.tBh1,    "tBh1[nhBh1]/D");
  tree->Branch("dtBh1",    event.dtBh1,   "dtBh1[nhBh1]/D");
  tree->Branch("deBh1",    event.deBh1,   "deBh1[nhBh1]/D");
  tree->Branch("btof",     event.btof,    "btof[nhBh1]/D");

  tree->Branch("nhBh2",   &event.nhBh2,   "nhBh2/I");
  tree->Branch("csBh2",    event.csBh2,   "csBh2[nhBh2]/I");
  tree->Branch("Bh2Seg",   event.Bh2Seg,  "Bh2Seg[nhBh2]/D");
  tree->Branch("tBh2",     event.tBh2,    "tBh2[nhBh2]/D");
  tree->Branch("t0Bh2",    event.t0Bh2,   "t0Bh2[nhBh2]/D");
  tree->Branch("dtBh2",    event.dtBh2,   "dtBh2[nhBh2]/D");
  tree->Branch("deBh2",    event.deBh2,   "deBh2[nhBh2]/D");

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("csTof",    event.csTof,   "csTof[nhTof]/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhBac",   &event.nhBac,   "nhBac/I");
  tree->Branch("BacSeg",   event.BacSeg,  "BacSeg[nhBac]/D");
  tree->Branch("tBac",     event.tBac,    "tBac[nhBac]/D");
  tree->Branch("deBac",    event.deBac,   "deBac[nhBac]/D");
  tree->Branch("nhPvac",   &event.nhPvac,   "nhPvac/I");
  tree->Branch("PvacSeg",   event.PvacSeg,  "PvacSeg[nhPvac]/D");
  tree->Branch("tPvac",     event.tPvac,    "tPvac[nhPvac]/D");
  tree->Branch("dePvac",    event.dePvac,   "dePvac[nhPvac]/D");
  tree->Branch("nhFac",   &event.nhFac,   "nhFac/I");
  tree->Branch("FacSeg",   event.FacSeg,  "FacSeg[nhFac]/D");
  tree->Branch("tFac",     event.tFac,    "tFac[nhFac]/D");
  tree->Branch("deFac",    event.deFac,   "deFac[nhFac]/D");

  //Fiber
  tree->Branch("nhBft",  &event.nhBft,  "nhBft/I");
  tree->Branch("csBft",   event.csBft,  "csBft[nhBft]/I");
  tree->Branch("tBft",    event.tBft,   "tBft[nhBft]/D");
  tree->Branch("wBft",    event.wBft,   "wBft[nhBft]/D");
  tree->Branch("BftPos",  event.BftPos, "BftPos[nhBft]/D");
  tree->Branch("nhSch",  &event.nhSch,  "nhSch/I");
  tree->Branch("csSch",   event.csSch,  "csSch[nhSch]/I");
  tree->Branch("tSch",    event.tSch,   "tSch[nhSch]/D");
  tree->Branch("wSch",    event.wSch,   "wSch[nhSch]/D");
  tree->Branch("SchPos",  event.SchPos, "SchPos[nhSch]/D");
  tree->Branch("nhFbh",  &event.nhFbh,  "nhFbh/I");
  tree->Branch("csFbh",   event.csFbh,  "csFbh[nhFbh]/I");
  tree->Branch("tFbh",    event.tFbh,   "tFbh[nhFbh]/D");
  tree->Branch("wFbh",    event.wFbh,   "wFbh[nhFbh]/D");
  tree->Branch("FbhPos",  event.FbhPos, "FbhPos[nhFbh]/D");

  //Beam DC
  tree->Branch("nlBcOut",   &event.nlBcOut,     "nlBcOut/I");
  tree->Branch("ntBcOut",   &event.ntBcOut,     "ntBcOut/I");
  tree->Branch("nhBcOut",    event.nhBcOut,     "nhBcOut[ntBcOut]/I");
  tree->Branch("chisqrBcOut",event.chisqrBcOut, "chisqrBcOut[ntBcOut]/D");
  tree->Branch("x0BcOut",    event.x0BcOut,     "x0BcOut[ntBcOut]/D");
  tree->Branch("y0BcOut",    event.y0BcOut,     "y0BcOut[ntBcOut]/D");
  tree->Branch("u0BcOut",    event.u0BcOut,     "u0BcOut[ntBcOut]/D");
  tree->Branch("v0BcOut",    event.v0BcOut,     "v0BcOut[ntBcOut]/D");
  tree->Branch("xtgtBcOut",  event.xtgtBcOut,   "xtgtBcOut[ntBcOut]/D");
  tree->Branch("ytgtBcOut",  event.ytgtBcOut,   "ytgtBcOut[ntBcOut]/D");
  tree->Branch("xbh2BcOut",  event.xbh2BcOut,   "xbh2BcOut[ntBcOut]/D");
  tree->Branch("ybh2BcOut",  event.ybh2BcOut,   "ybh2BcOut[ntBcOut]/D");

  tree->Branch("ntK18",      &event.ntK18,     "ntK18/I");
  tree->Branch("nhK18",       event.nhK18,     "nhK18[ntK18]/I");
  tree->Branch("chisqrK18",   event.chisqrK18, "chisqrK18[ntK18]/D");
  tree->Branch("pK18",        event.pK18,      "pK18[ntK18]/D");
  tree->Branch("xtgtK18",     event.xtgtK18,   "xtgtK18[ntK18]/D");
  tree->Branch("ytgtK18",     event.ytgtK18,   "ytgtK18[ntK18]/D");
  tree->Branch("utgtK18",     event.utgtK18,   "utgtK18[ntK18]/D");
  tree->Branch("vtgtK18",     event.vtgtK18,   "vtgtK18[ntK18]/D");
  tree->Branch("thetaK18",    event.thetaK18,  "thetaK18[ntK18]/D");

  //KURAMA
  tree->Branch("nlSdcIn",    &event.nlSdcIn,     "nlSdcIn/I");
  tree->Branch("ntSdcIn",    &event.ntSdcIn,     "ntSdcIn/I");
  tree->Branch("nhSdcIn",     event.nhSdcIn,     "nhSdcIn[ntSdcIn]/I");
  tree->Branch("chisqrSdcIn", event.chisqrSdcIn, "chisqrSdcIn[ntSdcIn]/D");
  tree->Branch("x0SdcIn",     event.x0SdcIn,     "x0SdcIn[ntSdcIn]/D");
  tree->Branch("y0SdcIn",     event.y0SdcIn,     "y0SdcIn[ntSdcIn]/D");
  tree->Branch("u0SdcIn",     event.u0SdcIn,     "u0SdcIn[ntSdcIn]/D");
  tree->Branch("v0SdcIn",     event.v0SdcIn,     "v0SdcIn[ntSdcIn]/D");

  tree->Branch("nlSdcOut",   &event.nlSdcOut,     "nlSdcOut/I");
  tree->Branch("ntSdcOut",   &event.ntSdcOut,     "ntSdcOut/I");
  tree->Branch("nhSdcOut",    event.nhSdcOut,     "nhSdcOut[ntSdcOut]/I");
  tree->Branch("chisqrSdcOut",event.chisqrSdcOut, "chisqrOut[ntSdcOut]/D");
  tree->Branch("x0SdcOut",    event.x0SdcOut,     "x0SdcOut[ntSdcOut]/D");
  tree->Branch("y0SdcOut",    event.y0SdcOut,     "y0SdcOut[ntSdcOut]/D");
  tree->Branch("u0SdcOut",    event.u0SdcOut,     "u0SdcOut[ntSdcOut]/D");
  tree->Branch("v0SdcOut",    event.v0SdcOut,     "v0SdcOut[ntSdcOut]/D");

  tree->Branch("ntKurama",     &event.ntKurama,     "ntKurama/I");
  tree->Branch("nhKurama",      event.nhKurama,     "nhKurama[ntKurama]/I");
  tree->Branch("chisqrKurama",  event.chisqrKurama, "chisqrKurama[ntKurama]/D");
  tree->Branch("path",          event.path,         "path[ntKurama]/D");
  tree->Branch("pKurama",       event.pKurama,      "pKurama[ntKurama]/D");
  tree->Branch("qKurama",       event.qKurama,      "qKurama[ntKurama]/D");
  tree->Branch("m2",            event.m2,           "m2[ntKurama]/D");
  tree->Branch("xtgtKurama",    event.xtgtKurama,   "xtgtKurama[ntKurama]/D");
  tree->Branch("ytgtKurama",    event.ytgtKurama,   "ytgtKurama[ntKurama]/D");
  tree->Branch("utgtKurama",    event.utgtKurama,   "utgtKurama[ntKurama]/D");
  tree->Branch("vtgtKurama",    event.vtgtKurama,   "vtgtKurama[ntKurama]/D");
  tree->Branch("thetaKurama",   event.thetaKurama,  "thetaKurama[ntKurama]/D");

  //Reaction
  tree->Branch("nKn",       &event.nKn,       "nKn/I");
  tree->Branch("nKp",       &event.nKp,       "nKp/I");
  tree->Branch("nKK",       &event.nKK,       "nKK/I");
  tree->Branch("vtx",        event.vtx,       "vtx[nKK]/D");
  tree->Branch("vty",        event.vty,       "vty[nKK]/D");
  tree->Branch("vtz",        event.vtz,       "vtz[nKK]/D");
  tree->Branch("closeDist",  event.closeDist, "closeDist[nKK]/D");
  tree->Branch("theta",      event.theta,     "theta[nKK]/D");
  tree->Branch("MissP",      event.MissP,     "MissP[nKK]/D");
  tree->Branch("MissPu",     event.MissPu,    "MissPu[nKK]/D");
  tree->Branch("MissPv",     event.MissPv,    "MissPv[nKK]/D");
  tree->Branch("MissMass",   event.MissMass,  "MissMass[nKK]/D");
  tree->Branch("thetaCM",    event.thetaCM,   "thetaCM[nKK]/D");
  tree->Branch("costCM",     event.costCM,    "costCM[nKK]/D");

  tree->Branch("xkn",   event.xkn,   "xkn[nKK]/D");
  tree->Branch("ykn",   event.ykn,   "ykn[nKK]/D");
  tree->Branch("ukn",   event.ukn,   "ukn[nKK]/D");
  tree->Branch("vkn",   event.vkn,   "vkn[nKK]/D");
  tree->Branch("xkp",   event.xkp,   "xkp[nKK]/D");
  tree->Branch("ykp",   event.ykp,   "ykp[nKK]/D");
  tree->Branch("ukp",   event.ukp,   "ukp[nKK]/D");
  tree->Branch("vkp",   event.vkp,   "vkp[nKK]/D");
  tree->Branch("pOrg",  event.pOrg,  "pOrg[nKK]/D");
  tree->Branch("pCalc", event.pCalc, "pCalc[nKK]/D");

  tree->Branch("nMM",     &event.nMM,     "nMM/I");
  tree->Branch("nXi",     &event.nXi,     "nXi/I");
  tree->Branch("nKXi",    &event.nKXi,    "nKXi/I");
  tree->Branch("vtxXi",    event.vtxXi,   "vtxXi[nKXi]/D");
  tree->Branch("vtyXi",    event.vtyXi,   "vtyXi[nKXi]/D");
  tree->Branch("vtzXi",    event.vtzXi,   "vtzXi[nKXi]/D");
  tree->Branch("closeXi",  event.closeXi, "closeXi[nKXi]/D");

  tree->Branch("ntXi",    &event.ntXi,    "ntXi/I");
  tree->Branch("x0Xi",     event.x0Xi,    "x0Xi[nKXi]/D");
  tree->Branch("y0Xi",     event.y0Xi,    "y0Xi[nKXi]/D");
  tree->Branch("u0Xi",     event.u0Xi,    "u0Xi[nKXi]/D");
  tree->Branch("v0Xi",     event.v0Xi,    "v0Xi[nKXi]/D");
  tree->Branch("thetaXi",  event.thetaXi, "thetaXi[nKXi]/D");
  tree->Branch("deXi",     event.deXi,    Form("deXi[%d][%d]/D",
					       NumOfLayersSsdIn, MaxHits ) );
  tree->Branch("deKaon",   event.deKaon,  Form("deKaon[%d][%d]/D",
					       NumOfLayersSsdIn, MaxHits ) );

  ////////// Bring Address From Dst
  TTreeCont[kKKAna]->SetBranchAddress("trigflag",  src.trigflag);
  TTreeCont[kKKAna]->SetBranchAddress("trigpat",   src.trigpat);
  TTreeCont[kKKAna]->SetBranchAddress("nhBh1",    &src.nhBh1);
  TTreeCont[kKKAna]->SetBranchAddress("csBh1",     src.csBh1);
  TTreeCont[kKKAna]->SetBranchAddress("Bh1Seg",    src.Bh1Seg);
  TTreeCont[kKKAna]->SetBranchAddress("tBh1",      src.tBh1);
  TTreeCont[kKKAna]->SetBranchAddress("dtBh1",     src.dtBh1);
  TTreeCont[kKKAna]->SetBranchAddress("deBh1",     src.deBh1);
  TTreeCont[kKKAna]->SetBranchAddress("btof",      src.btof);
  TTreeCont[kKKAna]->SetBranchAddress("nhBh2",    &src.nhBh2);
  TTreeCont[kKKAna]->SetBranchAddress("csBh2",     src.csBh2);
  TTreeCont[kKKAna]->SetBranchAddress("Bh2Seg",    src.Bh2Seg);
  TTreeCont[kKKAna]->SetBranchAddress("t0Bh2",     src.t0Bh2);
  TTreeCont[kKKAna]->SetBranchAddress("tBh2",      src.tBh2);
  TTreeCont[kKKAna]->SetBranchAddress("dtBh2",     src.dtBh2);
  TTreeCont[kKKAna]->SetBranchAddress("deBh2",     src.deBh2);
  TTreeCont[kKKAna]->SetBranchAddress("nhBac",    &src.nhBac);
  TTreeCont[kKKAna]->SetBranchAddress("BacSeg",    src.BacSeg);
  TTreeCont[kKKAna]->SetBranchAddress("tBac",      src.tBac);
  TTreeCont[kKKAna]->SetBranchAddress("deBac",     src.deBac);
  TTreeCont[kKKAna]->SetBranchAddress("nhPvac",   &src.nhPvac);
  TTreeCont[kKKAna]->SetBranchAddress("PvacSeg",   src.PvacSeg);
  TTreeCont[kKKAna]->SetBranchAddress("tPvac",     src.tPvac);
  TTreeCont[kKKAna]->SetBranchAddress("dePvac",    src.dePvac);
  TTreeCont[kKKAna]->SetBranchAddress("nhFac",    &src.nhFac);
  TTreeCont[kKKAna]->SetBranchAddress("FacSeg",    src.FacSeg);
  TTreeCont[kKKAna]->SetBranchAddress("tFac",      src.tFac);
  TTreeCont[kKKAna]->SetBranchAddress("deFac",     src.deFac);
  TTreeCont[kKKAna]->SetBranchAddress("nhTof",    &src.nhTof);
  TTreeCont[kKKAna]->SetBranchAddress("csTof",     src.csTof);
  TTreeCont[kKKAna]->SetBranchAddress("TofSeg",    src.TofSeg);
  TTreeCont[kKKAna]->SetBranchAddress("tTof",      src.tTof);
  TTreeCont[kKKAna]->SetBranchAddress("dtTof",     src.dtTof);
  TTreeCont[kKKAna]->SetBranchAddress("deTof",     src.deTof);

  TTreeCont[kKKAna]->SetBranchAddress("ntSdcIn",      &src.ntSdcIn      );
  TTreeCont[kKKAna]->SetBranchAddress("nlSdcIn",      &src.nlSdcIn      );
  TTreeCont[kKKAna]->SetBranchAddress("nhSdcIn",       src.nhSdcIn      );
  TTreeCont[kKKAna]->SetBranchAddress("chisqrSdcIn",   src.chisqrSdcIn  );
  TTreeCont[kKKAna]->SetBranchAddress("x0SdcIn",       src.x0SdcIn      );
  TTreeCont[kKKAna]->SetBranchAddress("y0SdcIn",       src.y0SdcIn      );
  TTreeCont[kKKAna]->SetBranchAddress("u0SdcIn",       src.u0SdcIn      );
  TTreeCont[kKKAna]->SetBranchAddress("v0SdcIn",       src.v0SdcIn      );
  TTreeCont[kKKAna]->SetBranchAddress("ntSdcOut",     &src.ntSdcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("nhSdcOut",      src.nhSdcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("chisqrSdcOut",  src.chisqrSdcOut );
  TTreeCont[kKKAna]->SetBranchAddress("x0SdcOut",      src.x0SdcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("y0SdcOut",      src.y0SdcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("u0SdcOut",      src.u0SdcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("v0SdcOut",      src.v0SdcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("ntXi",   &src.ntXi    );
  TTreeCont[kKKAna]->SetBranchAddress("x0Xi",    src.x0Xi    );
  TTreeCont[kKKAna]->SetBranchAddress("y0Xi",    src.y0Xi    );
  TTreeCont[kKKAna]->SetBranchAddress("u0Xi",    src.u0Xi    );
  TTreeCont[kKKAna]->SetBranchAddress("v0Xi",    src.v0Xi    );
  TTreeCont[kKKAna]->SetBranchAddress("deXi",    src.deXi    );
  TTreeCont[kKKAna]->SetBranchAddress("deKaon",  src.deKaon  );
  TTreeCont[kKKAna]->SetBranchAddress("thetaXi", src.thetaXi );
  TTreeCont[kKKAna]->SetBranchAddress("ntKurama", &src.ntKurama);
  TTreeCont[kKKAna]->SetBranchAddress("nhKurama",  src.nhKurama);
  TTreeCont[kKKAna]->SetBranchAddress("path",      src.path);
  TTreeCont[kKKAna]->SetBranchAddress("pKurama",   src.pKurama);
  TTreeCont[kKKAna]->SetBranchAddress("qKurama",   src.qKurama);
  TTreeCont[kKKAna]->SetBranchAddress("chisqrKurama", src.chisqrKurama);
  TTreeCont[kKKAna]->SetBranchAddress("xtgtKurama", src.xtgtKurama);
  TTreeCont[kKKAna]->SetBranchAddress("ytgtKurama", src.ytgtKurama);
  TTreeCont[kKKAna]->SetBranchAddress("utgtKurama", src.utgtKurama);
  TTreeCont[kKKAna]->SetBranchAddress("vtgtKurama", src.vtgtKurama);
  TTreeCont[kKKAna]->SetBranchAddress("m2",         src.m2);
  TTreeCont[kKKAna]->SetBranchAddress("ntBcOut",     &src.ntBcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("nlBcOut",     &src.nlBcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("nhBcOut",      src.nhBcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("chisqrBcOut",  src.chisqrBcOut );
  TTreeCont[kKKAna]->SetBranchAddress("x0BcOut",      src.x0BcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("y0BcOut",      src.y0BcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("u0BcOut",      src.u0BcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("v0BcOut",      src.v0BcOut     );
  TTreeCont[kKKAna]->SetBranchAddress("ntK18",     &src.ntK18);
  TTreeCont[kKKAna]->SetBranchAddress("nhK18",      src.nhK18);
  TTreeCont[kKKAna]->SetBranchAddress("chisqrK18",  src.chisqrK18);
  TTreeCont[kKKAna]->SetBranchAddress("pK18",    &src.pK18);
  TTreeCont[kKKAna]->SetBranchAddress("xtgtK18",  src.xtgtK18);
  TTreeCont[kKKAna]->SetBranchAddress("ytgtK18",  src.ytgtK18);
  TTreeCont[kKKAna]->SetBranchAddress("utgtK18",  src.utgtK18);
  TTreeCont[kKKAna]->SetBranchAddress("vtgtK18",  src.vtgtK18);
  TTreeCont[kKKAna]->SetBranchAddress("nhBft", &src.nhBft);
  TTreeCont[kKKAna]->SetBranchAddress("csBft",  src.csBft);
  TTreeCont[kKKAna]->SetBranchAddress("tBft",   src.tBft);
  TTreeCont[kKKAna]->SetBranchAddress("wBft",   src.wBft);
  TTreeCont[kKKAna]->SetBranchAddress("BftPos",  src.BftPos);
  TTreeCont[kKKAna]->SetBranchAddress("nhSch", &src.nhSch);
  TTreeCont[kKKAna]->SetBranchAddress("csSch",  src.csSch);
  TTreeCont[kKKAna]->SetBranchAddress("tSch",   src.tSch);
  TTreeCont[kKKAna]->SetBranchAddress("wSch",   src.wSch);
  TTreeCont[kKKAna]->SetBranchAddress("SchPos",  src.SchPos);
  TTreeCont[kKKAna]->SetBranchAddress("nhFbh", &src.nhFbh);
  TTreeCont[kKKAna]->SetBranchAddress("csFbh",  src.csFbh);
  TTreeCont[kKKAna]->SetBranchAddress("tFbh",   src.tFbh);
  TTreeCont[kKKAna]->SetBranchAddress("wFbh",   src.wFbh);
  TTreeCont[kKKAna]->SetBranchAddress("FbhPos",  src.FbhPos);
  TTreeCont[kKKAna]->SetBranchAddress("runnum", &src.runnum);
  TTreeCont[kKKAna]->SetBranchAddress("evnum",  &src.evnum);
  TTreeCont[kKKAna]->SetBranchAddress("spill",  &src.spill);
  TTreeCont[kKKAna]->SetBranchAddress("xpos",   &src.xpos);
  TTreeCont[kKKAna]->SetBranchAddress("ypos",   &src.ypos);
  TTreeCont[kKKAna]->SetBranchAddress("state",  &src.state);

  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<UserParamMan>("USER") );
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
