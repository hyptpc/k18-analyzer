/**
 *  file: DstKKAna.cc
 *  date: 2017.04.10
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include <filesystem_util.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "RootHelper.hh"
#include "DatabasePDG.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "KuramaLib.hh"
#include "K18TrackD2U.hh"
#include "FiberCluster.hh"
#include "NuclearMass.hh"
#include "MathTools.hh"
#include "MsTParamMan.hh"
#include "UserParamMan.hh"

#include "DstHelper.hh"

namespace
{
  using namespace root;
  using namespace dst;
  const std::string& class_name("DstKKAna");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  // const MsTParamMan&  gMsT  = MsTParamMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
}

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kKuramaTracking, kK18Tracking, kHodoscope, kEasiroc, kEMC,
      kOutFile, nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]", "[KuramaTracking]",
      "[K18Tracking]", "[Hodoscope]", "[Easiroc]", "[EMC]",
      "[OutFile]" };
  std::vector<TString> TreeName =
    { "", "", "kurama", "k18track", "hodo", "ea0c", "emc", "" };
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;

  const double pB_offset   = 1.000;
  const double pS_offset   = 1.000;
  const double pK18_offset = 0.000;
  const double pKURAMA_offset = 0.000;
  const double x_off = 0.000;
  const double y_off = 0.000;
  const double u_off = 0.000;
  const double v_off = 0.000;

  // const double CarbonMass  = 12.*NuclearMass::AMU();
  const double DeutronMass = 2.*NuclearMass::AMU()+0.01313672;
  const double DeltaC11    = 0.0106502;
  const double CoreMass    = 11.*NuclearMass::AMU()+DeltaC11;
  const double ThetaMass   = 1.530;
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
  int much; // debug
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
  double xtofKurama[MaxHits];
  double ytofKurama[MaxHits];
  double utofKurama[MaxHits];
  double vtofKurama[MaxHits];
  double tofsegKurama[MaxHits];

  double deKaon[NumOfLayersSsdIn][MaxHits];
  int    ntXi;
  double x0Xi[MaxHits];
  double y0Xi[MaxHits];
  double u0Xi[MaxHits];
  double v0Xi[MaxHits];
  double thetaXi[MaxHits];
  double deXi[NumOfLayersSsdIn][MaxHits];

  // Mass Trigger
  int rm_event;
  int rm_spill;
  int rm_accept;
  int rm_clear;
  int mst_accept;
  int mst_tdc[NumOfSegTOF];
  int mst_tof[NumOfSegTOF];
  int mst_sch[NumOfSegSCH];

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
  double MissMassCorr[MaxHits];
  double MissMassCorrDE[MaxHits];
  double BE[MaxHits];
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
  double pCorr[MaxHits];
  double pCorrDE[MaxHits];

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
  int    nhSch;
  // int    sch_hitpat[NumOfSegSCH];
  int    csSch[NumOfSegSCH];
  double tSch[NumOfSegSCH];
  double wSch[NumOfSegSCH];
  double SchPos[NumOfSegSCH];
  int    nhFbh;
  // int    fbh_hitpat[NumOfSegSCH];
  int    csFbh[NumOfSegCFBH];
  double tFbh[NumOfSegCFBH];
  double wFbh[NumOfSegCFBH];
  double FbhPos[NumOfSegCFBH];

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
  int much;
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
  double xtofKurama[MaxHits];
  double ytofKurama[MaxHits];
  double utofKurama[MaxHits];
  double vtofKurama[MaxHits];
  double tofsegKurama[MaxHits];

  double deKaon[NumOfLayersSsdIn][MaxHits];
  int    ntXi;
  double x0Xi[MaxHits];
  double y0Xi[MaxHits];
  double u0Xi[MaxHits];
  double v0Xi[MaxHits];
  double thetaXi[MaxHits];
  double deXi[NumOfLayersSsdIn][MaxHits];

  // Mass Trigger
  int rm_event;
  int rm_spill;
  int rm_accept;
  int rm_clear;
  int mst_accept;
  int mst_tdc[NumOfSegTOF];
  int mst_tof[NumOfSegTOF];
  int mst_sch[NumOfSegSCH];

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
  event.runnum     = 0;
  event.evnum      = 0;
  event.spill      = 0;
  event.rm_event   = -1;
  event.rm_spill   = -1;
  event.rm_accept  = -1;
  event.rm_clear   = -1;
  event.mst_accept = -1;
  for( int i=0; i<NumOfSegTOF; ++i ){
    event.mst_tof[i] = -1;
    event.mst_tdc[i] = -1;
  }
  for( int i=0; i<NumOfSegSCH; ++i ){
    event.mst_sch[i] = -1;
  }

  //Trigger
  for( int it=0; it<NumOfSegTrig; ++it ){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }

  //Hodoscope
  event.nhBh2  = 0;
  event.nhBh1  = 0;
  event.nhTof  = 0;
  event.nhBac  = 0;
  event.nhPvac = 0;
  event.nhFac  = 0;

  for( int it=0; it<NumOfSegBH1; it++ ){
    event.csBh1[it]  = 0;
    event.Bh1Seg[it] = -1;
    event.tBh1[it]   = -9999.;
    event.dtBh1[it]  = -9999.;
    event.deBh1[it]  = -9999.;
    event.btof[it]   = -9999.;
  }
  for( int it=0; it<NumOfSegBH2; it++ ){
    event.csBh2[it]  = 0;
    event.Bh2Seg[it] = -1;
    event.tBh2[it]   = -9999.;
    event.t0Bh2[it]  = -9999.;
    event.dtBh2[it]  = -9999.;
    event.deBh2[it]  = -9999.;
  }
  for( int it=0; it<NumOfSegTOF; it++ ){
    event.csTof[it]  = 0;
    event.TofSeg[it] = -1;
    event.tTof[it]   = -9999.;
    event.dtTof[it]  = -9999.;
    event.deTof[it]  = -9999.;
  }
  for( int it=0; it<NumOfSegBAC; it++ ){
    event.BacSeg[it] = -1;
    event.tBac[it]   = -9999.;
    event.deBac[it]  = -9999.;
  }
  for( int it=0; it<NumOfSegPVAC; it++ ){
    event.PvacSeg[it] = -1;
    event.tPvac[it]   = -9999.;
    event.dePvac[it]  = -9999.;
  }
  for( int it=0; it<NumOfSegFAC; it++ ){
    event.FacSeg[it] = -1;
    event.tFac[it]   = -9999.;
    event.deFac[it]  = -9999.;
  }

  //Fiber
  event.nhBft = 0;
  event.nhSch = 0;
  event.nhFbh = 0;
  for( int it=0; it<NumOfSegBFT; it++ ){
    event.csBft[it] = -999;
    event.tBft[it]  = -999.;
    event.wBft[it]   = -999.;
    event.BftPos[it]  = -999.;
  }
  for( int it=0; it<NumOfSegSCH; it++ ){
    event.csSch[it] = -999;
    event.tSch[it]  = -999.;
    event.wSch[it]   = -999.;
    event.SchPos[it]  = -999.;
  }
  for( int it=0; it<NumOfSegFBH; it++ ){
    event.csFbh[it] = -999;
    event.tFbh[it]  = -999.;
    event.wFbh[it]   = -999.;
    event.FbhPos[it]  = -999.;
  }

  //DC
  event.nlBcOut  = 0;
  event.nlSdcIn  = 0;
  event.nlSdcOut = 0;
  event.ntBcOut  = 0;
  event.much     = 0;
  event.ntSdcIn  = 0;
  event.ntSdcOut = 0;
  event.ntK18    = 0;
  event.ntKurama = 0;
  event.ntXi     = 0;

  //Beam DC
  for( int it=0; it<MaxHits; it++){
    event.nhBcOut[it]     = 0;
    event.chisqrBcOut[it] = -1.0;
    event.x0BcOut[it] = -9999.0;
    event.y0BcOut[it] = -9999.0;
    event.u0BcOut[it] = -9999.0;
    event.v0BcOut[it] = -9999.0;

    event.xtgtBcOut[it] = -9999.0;
    event.ytgtBcOut[it] = -9999.0;
    event.xbh2BcOut[it] = -9999.0;
    event.ybh2BcOut[it] = -9999.0;
  }

  for( int it=0; it<MaxHits; it++){
    event.nhK18[it]     = 0;
    event.chisqrK18[it] = -1.0;
    event.xtgtK18[it] = -9999.;
    event.ytgtK18[it] = -9999.;
    event.utgtK18[it] = -9999.;
    event.vtgtK18[it] = -9999.;
    event.pK18[it]    = -9999.;
    event.thetaK18[it] = -9999.;
  }

  //KURAMA DC
  for( int it=0; it<MaxHits; it++){
    event.nhSdcIn[it]     = 0;
    event.chisqrSdcIn[it] = -1.;
    event.x0SdcIn[it] = -9999.;
    event.y0SdcIn[it] = -9999.;
    event.u0SdcIn[it] = -9999.;
    event.v0SdcIn[it] = -9999.;

    event.nhSdcOut[it]     = 0;
    event.chisqrSdcOut[it] = -1.;
    event.u0SdcOut[it] = -9999.;
    event.v0SdcOut[it] = -9999.;
    event.x0SdcOut[it] = -9999.;
    event.y0SdcOut[it] = -9999.;

    event.nhKurama[it]     = 0;
    event.chisqrKurama[it] = -1.;
    event.xtgtKurama[it]   = -9999.;
    event.ytgtKurama[it]   = -9999.;
    event.utgtKurama[it]   = -9999.;
    event.vtgtKurama[it]   = -9999.;
    event.pKurama[it]      = -9999.;
    event.qKurama[it]      = -9999.;
    event.stof[it]         = -9999.;
    event.path[it]         = -9999.;
    event.m2[it]           = -9999.;
    event.thetaKurama[it]  = -9999.;
    event.xtofKurama[it]   = -9999.;
    event.ytofKurama[it]   = -9999.;
    event.utofKurama[it]   = -9999.;
    event.vtofKurama[it]   = -9999.;
    event.tofsegKurama[it] = -9999.;
  }

  for( int it=0; it<MaxHits; ++it ){
    event.x0Xi[it]    = -9999.;
    event.y0Xi[it]    = -9999.;
    event.u0Xi[it]    = -9999.;
    event.v0Xi[it]    = -9999.;
    event.thetaXi[it] = -9999.;
    for( int that=0; that<NumOfLayersSsdIn; ++that ){
      event.deKaon[that][it] = -9999.;
      event.deXi[that][it]   = -9999.;
    }
  }

  //Reaction
  event.nKn = 0;
  event.nKp = 0;
  event.nKK = 0;

  for( int it=0; it<MaxHits; ++it ){
    event.vtx[it]       = -9999.;
    event.vty[it]       = -9999.;
    event.vtz[it]       = -9999.;
    event.closeDist[it] = -9999.;
    event.theta[it]     = -9999.;
    event.thetaCM[it]   = -9999.;
    event.costCM[it]    = -9999.;
    event.MissP[it]     = -9999.;
    event.MissPu[it]    = -9999.;
    event.MissPv[it]    = -9999.;
    event.MissMass[it]  = -9999.;
    event.MissMassCorr[it]  = -9999.;
    event.MissMassCorrDE[it]  = -9999.;
    event.BE[it]        = -9999.0;

    event.xkn[it] = -9999.0;
    event.ykn[it] = -9999.0;
    event.ukn[it] = -9999.0;
    event.vkn[it] = -9999.0;
    event.xkp[it] = -9999.0;
    event.ykp[it] = -9999.0;
    event.ukp[it] = -9999.0;
    event.vkp[it] = -9999.0;
    event.pOrg[it] = -9999.0;
    event.pCalc[it] = -9999.0;
    event.pCorr[it] = -9999.0;
    event.pCorrDE[it] = -9999.0;
  }

  // EMC
  event.xpos  = -9999.;
  event.ypos  = -9999.;
  event.state = -1;

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

  static const double OffsetToF  = gUser.GetParameter("OffsetToF");

  static const double KaonMass    = pdg::KaonMass();
  static const double ProtonMass  = pdg::ProtonMass();
  static const double XiMass      = pdg::XiMass();
  static const double LambdaMass  = pdg::LambdaMass();

  if( ievent%10000==0 ){
    std::cout << "#D " << func_name << " Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  event.runnum   = src.runnum;
  event.evnum    = src.evnum;
  event.spill    = src.spill;

  event.xpos  = src.xpos;
  event.ypos  = src.ypos;
  event.state = src.state;

  event.ntBcOut  = src.ntBcOut;
  event.much     = src.much;
  event.ntSdcIn  = src.ntSdcIn;
  event.ntSdcOut = src.ntSdcOut;
  event.ntKurama = src.ntKurama;
  event.ntK18    = src.ntK18;
  event.ntXi     = src.ntXi;
  event.nhBft    = src.nhBft;
  event.nhBh1    = src.nhBh1;
  event.nhBh2    = src.nhBh2;
  event.nhBac    = src.nhBac;
  event.nhFbh    = src.nhFbh;
  event.nhPvac   = src.nhPvac;
  event.nhFac    = src.nhFac;
  event.nhSch    = src.nhSch;
  event.nhTof    = src.nhTof;

  const int ntBcOut  = event.ntBcOut;
  const int ntSdcIn  = event.ntSdcIn;
  const int ntSdcOut = event.ntSdcOut;
  const int ntKurama = event.ntKurama;
  const int ntK18    = event.ntK18;
  const int ntXi     = event.ntXi;
  const int nhBft    = event.nhBft;
  const int nhBh1    = event.nhBh1;
  const int nhBh2    = event.nhBh2;
  const int nhBac    = event.nhBac;
  const int nhFbh    = event.nhFbh;
  const int nhPvac   = event.nhPvac;
  const int nhFac    = event.nhFac;
  const int nhSch    = event.nhSch;
  const int nhTof    = event.nhTof;

#if 0
  std::cout << "#D DebugPrint" << std::endl
	    << " event  : " << std::setw(6) << ievent   << std::endl
	    << " BcOut  : " << std::setw(6) << ntBcOut  << std::endl
	    << " SdcIn  : " << std::setw(6) << ntSdcIn  << std::endl
	    << " SdcOut : " << std::setw(6) << ntSdcOut << std::endl
	    << " Kurama : " << std::setw(6) << ntKurama << std::endl
	    << " K18    : " << std::setw(6) << ntK18    << std::endl
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
    if( seg<0 ) continue;
    int tdc = src.trigflag[seg-1];
    if( tdc<=0 ) continue;
    event.trigpat[i] = seg;
    event.trigflag[seg-1] = tdc;
  }

  HF1( 1, 0. );

  // if( event.ntKurama==0 ) return true;
  // HF1( 1, 1. );
  // if( event.ntK18==0 ) return true;
  // HF1( 1, 2. );
  // if( event.nhBh1==0 ) return true;
  // HF1( 1, 3. );
  // if( event.nhBh2==0 ) return true;
  // HF1( 1, 4. );
  // if( event.nhTof==0 ) return true;
  // HF1( 1, 5. );
  // if( event.nhSch==0 ) return true;
  // HF1( 1, 6. );
  // if( event.nhFbh==0 ) return true;
  // HF1( 1, 7. );

  std::vector <ThreeVector> KnPCont, KnXCont;
  std::vector <ThreeVector> KpPCont, KpXCont;

  // BFT
  for( int i=0; i<nhBft; ++i ){
    event.csBft[i]  = src.csBft[i];
    event.tBft[i]   = src.tBft[i];
    event.wBft[i]   = src.wBft[i];
    event.BftPos[i] = src.BftPos[i];
  }

  // BH1
  for( int i=0; i<nhBh1; ++i ){
    event.csBh1[i]  = src.csBh1[i];
    event.Bh1Seg[i] = src.Bh1Seg[i];
    event.tBh1[i]   = src.tBh1[i];
    event.dtBh1[i]  = src.dtBh1[i];
    event.deBh1[i]  = src.deBh1[i];
    event.btof[i]   = src.btof[i];
  }

  // BH2
  double time0 = -9999.;
  for( int i=0; i<nhBh2; ++i ){
    event.csBh2[i]  = src.csBh2[i];
    event.Bh2Seg[i] = src.Bh2Seg[i];
    event.tBh2[i]   = src.tBh2[i];
    event.t0Bh2[i]  = src.t0Bh2[i];
    event.dtBh2[i]  = src.dtBh2[i];
    event.deBh2[i]  = src.deBh2[i];
    time0 = event.t0Bh2[i];
  }

  // BAC
  for( int i=0; i<nhBac; ++i ){
    event.BacSeg[i] = src.BacSeg[i];
    event.tBac[i]   = src.tBac[i];
    event.deBac[i]  = src.deBac[i];
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
  for( int i=0; i<nhTof; ++i ){
    event.csTof[i]  = src.csTof[i];
    event.TofSeg[i] = src.TofSeg[i];
    event.tTof[i]   = src.tTof[i];
    event.dtTof[i]  = src.dtTof[i];
    event.deTof[i]  = src.deTof[i];
  }

  ////////// BcOut
  event.nlBcOut = src.nlBcOut;
  for( int it=0; it<ntBcOut; ++it ){
    event.nhBcOut[it]     = src.nhBcOut[it];
    event.chisqrBcOut[it] = src.chisqrBcOut[it];
    event.x0BcOut[it]     = src.x0BcOut[it];
    event.y0BcOut[it]     = src.y0BcOut[it];
    event.u0BcOut[it]     = src.u0BcOut[it];
    event.v0BcOut[it]     = src.v0BcOut[it];
  }

  ////////// SdcIn
  event.nlSdcIn = src.nlSdcIn;
  for( int it=0; it<ntSdcIn; ++it ){
    event.nhSdcIn[it]     = src.nhSdcIn[it];
    event.chisqrSdcIn[it] = src.chisqrSdcIn[it];
    event.x0SdcIn[it]     = src.x0SdcIn[it];
    event.y0SdcIn[it]     = src.y0SdcIn[it];
    event.u0SdcIn[it]     = src.u0SdcIn[it];
    event.v0SdcIn[it]     = src.v0SdcIn[it];
  }

  ////////// SdcOut
  event.nlSdcOut = src.nlSdcOut;
  for( int it=0; it<ntSdcOut; ++it ){
    event.nhSdcOut[it]     = src.nhSdcOut[it];
    event.chisqrSdcOut[it] = src.chisqrSdcOut[it];
    event.x0SdcOut[it]     = src.x0SdcOut[it];
    event.y0SdcOut[it]     = src.y0SdcOut[it];
    event.u0SdcOut[it]     = src.u0SdcOut[it];
    event.v0SdcOut[it]     = src.v0SdcOut[it];
  }

  ////////// Kaon
  for( int it=0; it<ntKurama; ++it ){
    for( int il=0; il<NumOfLayersSsdIn; ++il ){
      event.deKaon[il][it] = src.deKaon[il][it];
    }
  }

  ////////// Xi
  for( int it=0; it<ntXi; ++it ){
    event.x0Xi[it] = src.x0Xi[it];
    event.y0Xi[it] = src.y0Xi[it];
    event.u0Xi[it] = src.u0Xi[it];
    event.v0Xi[it] = src.v0Xi[it];
    event.thetaXi[it] = src.thetaXi[it];
    for( int il=0; il<NumOfLayersSsdIn; ++il ){
      event.deXi[il][it] = src.deXi[il][it];
    }
  }

  ////////// Kaon Plus
  for( int itKurama=0; itKurama<ntKurama; ++itKurama ){
    int nh= src.nhKurama[itKurama];
    double chisqr = src.chisqrKurama[itKurama];
    double p = src.pKurama[itKurama];
    double x = src.xtgtKurama[itKurama];
    double y = src.ytgtKurama[itKurama];
    double u = src.utgtKurama[itKurama];
    double v = src.vtgtKurama[itKurama];
    double theta = src.thetaKurama[itKurama];
    double pt = p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    if( std::isnan( Pos.Mag() ) ) continue;
    //Calibration
    ThreeVector PosCorr( x+x_off, y+y_off, 0. );
    ThreeVector MomCorr( pt*(u+u_off), pt*(v+v_off), pt );
    double path = src.path[itKurama];
    double xt = PosCorr.x();
    double yt = PosCorr.y();
    double zt = PosCorr.z();
    double pCorr = MomCorr.Mag();
    double q  = src.qKurama[itKurama];
    double ut = xt/zt;
    double vt = yt/zt;

    // SdcOut vs TOF
    double TofSegKurama = src.tofsegKurama[itKurama];
    double stof = -9999.;
    double m2   = -9999.;
    // w/  TOF
    for( int j=0; j<nhTof; ++j ){
      double seg = src.TofSeg[j];
      if( seg == TofSegKurama ){
	// event.tTof[i]  = src.tTof[j];
	// event.dtTof[i] = src.dtTof[j];
	// event.deTof[i] = src.deTof[j];
	stof  = src.tTof[j] - time0 + OffsetToF;
	m2 = Kinematics::MassSquare( pCorr, path, stof );
      }
    }

    // w/o TOF
    // double minres = 1.0e10;
    // for( int j=0; j<nhTof; ++j ){
    //   double seg = src.TofSeg[j];
    //   double res = TofSegKurama - seg;
    //   if( std::abs( res ) < minres && std::abs( res ) < 1. ){
    // 	minres = res;
    // 	stof = src.tTof[j] - time0 + OffsetToF;
    //   }
    // }
    event.nhKurama[itKurama]   = nh;
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
    HF1( 3202, double(nh) );
    HF1( 3203, chisqr );
    HF1( 3204, xt ); HF1( 3205, yt );
    HF1( 3206, ut ); HF1( 3207, vt );
    HF2( 3208, xt, ut ); HF2( 3209, yt, vt );
    HF2( 3210, xt, yt );
    HF1( 3211, pCorr );
    HF1( 3212, path );

    //HF1( 3213, m2 );

    //     double xTof=(clTof->MeanSeg()-7.5)*70.;
    //     double yTof=(clTof->TimeDif())*800./12.;
    //     HF2( 4011, clTof->MeanSeg()-0.5, xtof );
    //     HF2( 4013, clTof->TimeDif(), ytof );
    //     HF1( 4015, xtof-xTof ); HF1( 4017, ytof-yTof );
    //     HF2( 4019, xtof-xTof, ytof-yTof );

    //       double ttof=clTof->CMeanTime()-time0,
    //       HF1( 322, clTof->ClusterSize() );
    //       HF1( 323, clTof->MeanSeg()-0.5 );
    //       HF1( 324, ttof );
    //       HF1( 325, clTof->DeltaE() );
    //       double u0in=trIn->GetU0();
    //       HF2( 4001, u0in, ttof ); HF2( 4003, u0in, ttof+12.5*u0in );
    HF1( 4202, double(nh) );
    HF1( 4203, chisqr );
    HF1( 4204, xt ); HF1( 4205, yt );
    HF1( 4206, ut ); HF1( 4207, vt );
    HF2( 4208, xt, ut ); HF2( 4209, yt, vt );
    HF2( 4210, xt, yt );
    HF1( 4211, pCorr );
    HF1( 4212, path );
    KpXCont.push_back( PosCorr );
    KpPCont.push_back( MomCorr );
  }

  // if( KpPCont.size()==0 ) return true;

  HF1( 1, 8. );

  ////////// Kaon Minus
  for( int itK18=0; itK18<ntK18; ++itK18 ){
    double nh = src.nhK18[itK18];
    double chisqr = src.chisqrK18[itK18];
    //Calibration
    double p = src.pK18[itK18]*pB_offset+pK18_offset;
    double x = src.xtgtK18[itK18];
    double y = src.ytgtK18[itK18];
    double u = src.utgtK18[itK18];
    double v = src.vtgtK18[itK18];
    event.nhK18[itK18]     = nh;
    event.chisqrK18[itK18] = chisqr;
    event.pK18[itK18]      = p;
    event.xtgtK18[itK18]   = x;
    event.ytgtK18[itK18]   = y;
    event.utgtK18[itK18]   = u;
    event.vtgtK18[itK18]   = v;
    // double loss_bh2 = 1.09392e-3;
    // p = p - loss_bh2;
    double pt=p/std::sqrt(1.+u*u+v*v);
    ThreeVector Pos( x, y, 0. );
    ThreeVector Mom( pt*u, pt*v, pt );
    //double xo=trOut->GetX0(), yo=trOut->GetY0();
    HF1( 4104, p );
    HF1( 4105, x ); HF1( 4106, y );
    //HF1( 4107, xo ); HF1( 4108, yo ); HF1( 4109, u ); HF1( 4110, v );

    KnPCont.push_back(Mom); KnXCont.push_back(Pos);
  }

  // if( KnPCont.size()==0 ) return true;

  HF1( 1, 9. );

  //MissingMass
  int nKn = KnPCont.size();
  int nKp = KpPCont.size();
  event.nKn = nKn;
  event.nKp = nKp;
  event.nKK = nKn*nKp;
  HF1( 4101, double(nKn));
  HF1( 4201, double(nKp) );
  int nkk=0;
  for( int ikp=0; ikp<nKp; ++ikp ){
    ThreeVector pkp = KpPCont[ikp], xkp = KpXCont[ikp];
    for( int ikn=0; ikn<nKn; ++ikn ){
      ThreeVector pkn  = KnPCont[ikn], xkn = KnXCont[ikn];
      ThreeVector vert = Kinematics::VertexPoint( xkn, xkp, pkn, pkp );
      // std::cout << "vertex : " << vert << " " << vert.Mag() << std::endl;
      double closedist = Kinematics::closeDist( xkn, xkp, pkn, pkp );

      double us = pkp.x()/pkp.z(), vs = pkp.y()/pkp.z();
      double ub = pkn.x()/pkn.z(), vb = pkn.y()/pkn.z();
      double cost = pkn*pkp/(pkn.Mag()*pkp.Mag());

      double pk0   = pkp.Mag();
      double pCorr = pk0;

      ThreeVector pkpCorr( pCorr*pkp.x()/pkp.Mag(),
			   pCorr*pkp.y()/pkp.Mag(),
			   pCorr*pkp.z()/pkp.Mag() );

      ThreeVector pknCorrDE = Kinematics::CorrElossIn( pkn, xkn, vert, KaonMass );
      ThreeVector pkpCorrDE = Kinematics::CorrElossOut( pkpCorr, xkp, vert, KaonMass );

      LorentzVector LvKn( pkn, std::sqrt( KaonMass*KaonMass+pkn.Mag2() ) );
      LorentzVector LvKnCorrDE( pknCorrDE, sqrt( KaonMass*KaonMass+pknCorrDE.Mag2() ) );

      LorentzVector LvKp( pkp, std::sqrt( KaonMass*KaonMass+pkp.Mag2() ) );
      LorentzVector LvKpCorr( pkpCorr, std::sqrt( KaonMass*KaonMass+pkpCorr.Mag2() ) );
      LorentzVector LvKpCorrDE( pkpCorrDE, std::sqrt( KaonMass*KaonMass+pkpCorrDE.Mag2() ) );

      LorentzVector LvC( 0., 0., 0., ProtonMass );
      LorentzVector LvCore( 0., 0., 0., 0. );

      LorentzVector LvRc       = LvKn+LvC-LvKp;
      LorentzVector LvRcCorr   = LvKn+LvC-LvKpCorr;
      LorentzVector LvRcCorrDE = LvKnCorrDE+LvC-LvKpCorrDE;
      double MisP  = std::sqrt( LvRc.E()*LvRc.E() - LvRc.Mag2() );
      double MisPu = LvRc.Px()/LvRc.Pz();
      double MisPv = LvRc.Py()/LvRc.Pz();
      double MisMass       = LvRc.Mag();//-LvC.Mag();
      double MisMassCorr   = LvRcCorr.Mag();//-LvC.Mag();
      double MisMassCorrDE = LvRcCorrDE.Mag();//-LvC.Mag();

      double BE       = LvRc.Mag()       - ( CoreMass+LambdaMass );
      // double BECorr   = LvRcCorr.Mag()   - ( CoreMass+LambdaMass );
      // double BECorrDE = LvRcCorrDE.Mag() - ( CoreMass+LambdaMass );

      //Primary frame
      LorentzVector PrimaryLv = LvKn+LvC;
      double TotalEnergyCM = PrimaryLv.Mag();
      ThreeVector beta( 1/PrimaryLv.E()*PrimaryLv.Vect() );

      //CM
      double TotalMomCM
	= 0.5*std::sqrt(( TotalEnergyCM*TotalEnergyCM
			  -( KaonMass+XiMass )*( KaonMass+XiMass ))
			*( TotalEnergyCM*TotalEnergyCM
			   -( KaonMass-XiMass )*( KaonMass-XiMass )))/TotalEnergyCM;

      double costLab = cost;
      double cottLab = costLab/std::sqrt(1.-costLab*costLab);
      double bt=beta.Mag(), gamma=1./std::sqrt(1.-bt*bt);
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

      if (nkk<MaxHits) {
	event.vtx[nkk]=vert.x();
	event.vty[nkk]=vert.y();
	event.vtz[nkk]=vert.z();
	event.closeDist[nkk] = closedist;
	event.theta[nkk]     = std::acos(cost)*math::Rad2Deg();
	event.thetaCM[nkk]   = std::acos(costCM)*math::Rad2Deg();
	event.costCM[nkk]    = costCM;

	event.MissP[nkk]          = MisP;
	event.MissPu[nkk]         = MisPu;
	event.MissPv[nkk]         = MisPv;
	event.MissMass[nkk]       = MisMass;
	event.MissMassCorr[nkk]   = MisMassCorr;
	event.MissMassCorrDE[nkk] = MisMassCorrDE;

	event.BE[nkk]=BE;

	event.xkp[nkk] = xkp.x();
	event.ykp[nkk] = xkp.y();
	event.ukp[nkk] = us;
	event.vkp[nkk] = vs;

	event.xkn[nkk] = xkn.x();
	event.ykn[nkk] = xkn.y();
	event.ukn[nkk] = ub;
	event.vkn[nkk] = vb;
	event.pOrg[nkk] = pk0;
	event.pCalc[nkk] = KaonMom;
	event.pCorr[nkk] = pCorr;
	event.pCorrDE[nkk] = pkpCorrDE.Mag();
	nkk++;
      }

      HF1( 5001, vert.z() );

      HF1( 5002, MisMass );
      HF2( 5011, MisMass, us );
      HF2( 5012, MisMass, vs );
      HF2( 5013, MisMass, ub );
      HF2( 5014, MisMass, vb );
    }
  }

  HF1( 1, 10. );

  //Final Hodoscope histograms
  for( int i=0; i<nhBh2; ++i ){
    HF1( 152, event.csBh2[i] );
    HF1( 153, event.Bh2Seg[i] );
    HF1( 154, event.tBh2[i] );
    HF1( 155, event.deBh2[i] );
  }

  for( int i=0; i<nhBh1; ++i ){
    HF1( 252, event.csBh1[i] );
    HF1( 253, event.Bh1Seg[i] );
    HF1( 254, event.tBh1[i] );
    HF1( 255, event.deBh1[i] );
  }

  for( int i=0; i<nhTof; ++i ){
    HF1( 352, event.csTof[i] );
    HF1( 353, event.TofSeg[i] );
    HF1( 354, event.tTof[i] );
    HF1( 355, event.deTof[i] );
  }

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
  HB1( 1, "Status", 60, 0., 60. );

  HB1( 101, "#Clusters BH2",   7, 0., 7. );
  HB1( 102, "ClusterSize BH2", 7, 0., 7. );
  HB1( 103, "HitPat BH2", 10, 0., 10. );
  HB1( 104, "MeanTime BH2", 200, -10., 10. );
  HB1( 105, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 112, "ClusterSize BH2", 7, 0., 7. );
  HB1( 113, "HitPat BH2", 10, 0., 10. );
  HB1( 114, "MeanTime BH2", 200, -10., 10. );
  HB1( 115, "Delta-E BH2", 200, -0.5, 4.5 );

  HB1( 122, "ClusterSize BH2 [T0]", 7, 0., 7. );
  HB1( 123, "HitPat BH2 [T0]", 10, 0., 10. );
  HB1( 124, "MeanTime BH2 [T0]", 200, -10., 10. );
  HB1( 125, "Delta-E BH2 [T0]", 200, -0.5, 4.5 );

  HB1( 152, "ClusterSize BH2 [kk]", 7, 0., 7. );
  HB1( 153, "HitPat BH2 [kk]", 10, 0., 10. );
  HB1( 154, "MeanTime BH2 [kk]", 200, -10., 10. );
  HB1( 155, "Delta-E BH2 [kk]", 200, -0.5, 4.5 );

  HB1( 201, "#Clusters BH1",  11, 0., 11. );
  HB1( 202, "ClusterSize BH1",11, 0., 11. );
  HB1( 203, "HitPat BH1", 11, 0., 11. );
  HB1( 204, "MeanTime BH1", 200, -10., 10. );
  HB1( 205, "Delta-E BH1", 200, -0.5, 4.5 );
  HB1( 206, "Beam ToF", 200, -10., 10. );

  HB1( 211, "#Clusters BH1 [pi]",  11, 0., 11. );
  HB1( 212, "ClusterSize BH1 [pi]",11, 0., 11. );
  HB1( 213, "HitPat BH1 [pi]", 11, 0., 11. );
  HB1( 214, "MeanTime BH1 [pi]", 200, -10., 10. );
  HB1( 215, "Delta-E BH1 [pi]", 200, -0.5, 4.5 );

  HB1( 252, "ClusterSize BH1 [kk]",11, 0., 11. );
  HB1( 253, "HitPat BH1 [kk]", 11, 0., 11. );
  HB1( 254, "MeanTime BH1 [kk]", 200, -10., 10. );
  HB1( 255, "Delta-E BH1 [kk]", 200, -0.5, 4.5 );
  HB1( 256, "Beam ToF [kk]", 200, -10., 10. );

  HB1( 301, "#Clusters Tof",  32, 0., 32. );
  HB1( 302, "ClusterSize Tof",32, 0., 32. );
  HB1( 303, "HitPat Tof", 32, 0., 32. );
  HB1( 304, "TimeOfFlight Tof", 500, -50., 100. );
  HB1( 305, "Delta-E Tof", 200, -0.5, 4.5 );

  HB1( 311, "#Clusters Tof [Good]",  32, 0., 32. );
  HB1( 312, "ClusterSize Tof [Good]",32, 0., 32. );
  HB1( 313, "HitPat Tof [Good]", 32, 0., 32. );
  HB1( 314, "TimeOfFlight Tof [Good]", 500, -50., 100. );
  HB1( 315, "Delta-E Tof [Good]", 200, -0.5, 4.5 );

  HB1( 352, "ClusterSize Tof [kk]",32, 0., 32. );
  HB1( 353, "HitPat Tof [kk]", 32, 0., 32. );
  HB1( 354, "TimeOfFlight Tof [kk]", 500, -50., 100. );
  HB1( 355, "Delta-E Tof [kk]", 200, -0.5, 4.5 );
  HB2( 501, "SegLC%SegTOF", 32, 0., 32., 28, 0., 28. );

  HB1( 1001, "#Tracks SdcIn", 10, 0., 10. );
  HB1( 1002, "#Hits SdcIn", 20, 0., 20. );
  HB1( 1003, "Chisqr SdcIn", 200, 0., 100. );
  HB1( 1004, "X0 SdcIn", 500, -200., 200. );
  HB1( 1005, "Y0 SdcIn", 500, -200., 200. );
  HB1( 1006, "U0 SdcIn",  700, -0.35, 0.35 );
  HB1( 1007, "V0 SdcIn",  400, -0.20, 0.20 );
  HB2( 1008, "X0%U0 SdcIn", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1009, "Y0%V0 SdcIn", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1010, "X0%Y0 SdcIn", 100, -200., 200., 100, -200., 200. );

  HB1( 1101, "#Tracks BcOut", 10, 0., 10. );
  HB1( 1102, "#Hits BcOut", 20, 0., 20. );
  HB1( 1103, "Chisqr BcOut", 200, 0., 100. );
  HB1( 1104, "X0 BcOut", 500, -200., 200. );
  HB1( 1105, "Y0 BcOut", 500, -200., 200. );
  HB1( 1106, "U0 BcOut",  700, -0.35, 0.35 );
  HB1( 1107, "V0 BcOut",  400, -0.20, 0.20 );
  HB2( 1108, "X0%U0 BcOut", 120, -200., 200., 100, -0.35, 0.35 );
  HB2( 1109, "Y0%V0 BcOut", 100, -200., 200., 100, -0.20, 0.20 );
  HB2( 1110, "X0%Y0 BcOut", 100, -200., 200., 100, -200., 200. );
  HB1( 1111, "Xtgt BcOut", 500, -200., 200. );
  HB1( 1112, "Ytgt BcOut", 500, -200., 200. );
  HB2( 1113, "Xtgt%Ytgt BcOut", 100, -200., 200., 100, -200., 200. );

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

  HB1( 2201, "#Tracks K18", 10, 0., 10. );
  HB1( 2202, "#Hits K18", 30, 0., 30. );
  HB1( 2203, "Chisqr K18", 500, 0., 50. );
  HB1( 2204, "P K18", 1000, 0.5, 2.0 );
  HB1( 2251, "#Tracks K18 [Good]", 10, 0., 10. );

  HB1( 3001, "#Tracks Kurama", 10, 0., 10. );
  HB1( 3002, "#Hits Kurama", 30, 0., 30. );
  HB1( 3003, "Chisqr Kurama", 500, 0., 500. );
  HB1( 3004, "Xtgt Kurama", 500, -200., 200. );
  HB1( 3005, "Ytgt Kurama", 500, -100., 100. );
  HB1( 3006, "Utgt Kurama", 400, -0.35, 0.35 );
  HB1( 3007, "Vtgt Kurama", 200, -0.20, 0.20 );
  HB2( 3008, "Xtgt%U Kurama", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3009, "Ytgt%V Kurama", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3010, "Xtgt%Ytgt Kurama", 100, -200., 200., 100, -100., 100. );
  HB1( 3011, "P Kurama", 200, 0.0, 1.0 );
  HB1( 3012, "PathLength Kurama", 600, 3000., 6000. );
  HB1( 3013, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3101, "#Tracks Kurama [Good]", 10, 0., 10. );
  HB1( 3102, "#Hits Kurama [Good]", 30, 0., 30. );
  HB1( 3103, "Chisqr Kurama [Good]", 500, 0., 500. );
  HB1( 3104, "Xtgt Kurama [Good]", 500, -200., 200. );
  HB1( 3105, "Ytgt Kurama [Good]", 500, -100., 100. );
  HB1( 3106, "Utgt Kurama [Good]", 700, -0.35, 0.35 );
  HB1( 3107, "Vtgt Kurama [Good]", 400, -0.20, 0.20 );
  HB2( 3108, "Xtgt%U Kurama [Good]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3109, "Ytgt%V Kurama [Good]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3110, "Xtgt%Ytgt Kurama [Good]", 100, -200., 200., 100, -100., 100. );
  HB1( 3111, "P Kurama [Good]", 200, 0.0, 1.0 );
  HB1( 3112, "PathLength Kurama [Good]", 600, 3000., 6000. );
  HB1( 3113, "MassSqr", 600, -1.2, 1.2 );

  HB1( 3201, "#Tracks Kurama [Good2]", 10, 0., 10. );
  HB1( 3202, "#Hits Kurama [Good2]", 30, 0., 30. );
  HB1( 3203, "Chisqr Kurama [Good2]", 500, 0., 500. );
  HB1( 3204, "Xtgt Kurama [Good2]", 500, -200., 200. );
  HB1( 3205, "Ytgt Kurama [Good2]", 500, -100., 100. );
  HB1( 3206, "Utgt Kurama [Good2]", 700, -0.35, 0.35 );
  HB1( 3207, "Vtgt Kurama [Good2]", 400, -0.20, 0.20 );
  HB2( 3208, "Xtgt%U Kurama [Good2]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 3209, "Ytgt%V Kurama [Good2]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 3210, "Xtgt%Ytgt Kurama [Good2]", 100, -200., 200., 100, -100., 100. );
  HB1( 3211, "P Kurama [Good2]", 200, 0.0, 1.0 );
  HB1( 3212, "PathLength Kurama [Good2]", 600, 3000., 6000. );
  HB1( 3213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 4101, "#Tracks K18 [KuramaP]", 10, 0., 10. );
  HB1( 4102, "#Hits K18 [KuramaP]", 30, 0., 30. );
  HB1( 4103, "Chisqr K18 [KuramaP]", 500, 0., 50. );
  HB1( 4104, "P K18 [KuramaP]", 500, 0.5, 2.0 );
  HB1( 4105, "Xtgt K18 [KuramaP]", 500, -200., 200. );
  HB1( 4106, "Ytgt K18 [KuramaP]", 500, -100., 100. );
  HB1( 4107, "Xout K18 [KuramaP]", 500, -200., 200. );
  HB1( 4108, "Yout K18 [KuramaP]", 500, -100., 100. );
  HB1( 4109, "Uout K18 [KuramaP]", 400, -0.35, 0.35 );
  HB1( 4110, "Vout K18 [KuramaP]", 200, -0.20, 0.20 );
  HB1( 4111, "Xin  K18 [KuramaP]", 500, -100., 100. );
  HB1( 4112, "Yin  K18 [KuramaP]", 500, -100., 100. );
  HB1( 4113, "Uin  K18 [KuramaP]", 400, -0.35, 0.35 );
  HB1( 4114, "Vin  K18 [KuramaP]", 200, -0.20, 0.20 );

  HB1( 4201, "#Tracks Kurama [Proton]", 10, 0., 10. );
  HB1( 4202, "#Hits Kurama [Proton]", 30, 0., 30. );
  HB1( 4203, "Chisqr Kurama [Proton]", 500, 0., 500. );
  HB1( 4204, "Xtgt Kurama [Proton]", 500, -200., 200. );
  HB1( 4205, "Ytgt Kurama [Proton]", 500, -100., 100. );
  HB1( 4206, "Utgt Kurama [Proton]", 700, -0.35, 0.35 );
  HB1( 4207, "Vtgt Kurama [Proton]", 400, -0.20, 0.20 );
  HB2( 4208, "Xtgt%U Kurama [Proton]", 100, -200., 200., 100, -0.35, 0.35 );
  HB2( 4209, "Ytgt%V Kurama [Proton]", 100, -100., 100., 100, -0.20, 0.20 );
  HB2( 4210, "Xtgt%Ytgt Kurama [Proton]", 100, -200., 200., 100, -100., 100. );
  HB1( 4211, "P Kurama [Proton]", 200, 0.0, 1.0 );
  HB1( 4212, "PathLength Kurama [Proton]", 600, 3000., 6000. );
  HB1( 4213, "MassSqr", 600, -1.2, 1.2 );

  HB1( 5001, "Zvert [KK]", 1000, -1000., 1000. );
  HB1( 5002, "MissingMass [KK]", 1000, 0.0, 2.0 );

  HB2( 5011, "MissingMass%Us", 200, 0.0, 2.50, 100, -0.40, 0.40 );
  HB2( 5012, "MissingMass%Vs", 200, 0.0, 2.50, 100, -0.20, 0.20 );
  HB2( 5013, "MissingMass%Ub", 200, 0.0, 2.50, 100, -0.30, 0.30 );
  HB2( 5014, "MissingMass%Vb", 200, 0.0, 2.50, 100, -0.10, 0.10 );

  ////////////////////////////////////////////
  //Tree
  HBTree("kk","tree of KKAna");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  //Trigger
  tree->Branch("trigpat",   event.trigpat,  Form( "trigpat[%d]/I", NumOfSegTrig ) );
  tree->Branch("trigflag",  event.trigflag, Form( "trigflag[%d]/I", NumOfSegTrig ) );

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
  tree->Branch("much",       &event.much,        "much/I");
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

  tree->Branch("deKaon",   event.deKaon,  Form("deKaon[%d][%d]/D",
					       NumOfLayersSsdIn, MaxHits ) );
  tree->Branch("ntXi",    &event.ntXi,    "ntXi/I");
  tree->Branch("x0Xi",     event.x0Xi,    "x0Xi[ntXi]/D");
  tree->Branch("y0Xi",     event.y0Xi,    "y0Xi[ntXi]/D");
  tree->Branch("u0Xi",     event.u0Xi,    "u0Xi[ntXi]/D");
  tree->Branch("v0Xi",     event.v0Xi,    "v0Xi[ntXi]/D");
  tree->Branch("thetaXi",  event.thetaXi, "thetaXi[ntXi]/D");
  tree->Branch("deXi",     event.deXi,    Form("deXi[%d][%d]/D",
					       NumOfLayersSsdIn, MaxHits ) );

  tree->Branch("ntKurama",     &event.ntKurama,     "ntKurama/I");
  tree->Branch("nhKurama",      event.nhKurama,     "nhKurama[ntKurama]/I");
  tree->Branch("chisqrKurama",  event.chisqrKurama, "chisqrKurama[ntKurama]/D");
  tree->Branch("stof",          event.stof,         "stof[ntKurama]/D");
  tree->Branch("path",          event.path,         "path[ntKurama]/D");
  tree->Branch("pKurama",       event.pKurama,      "pKurama[ntKurama]/D");
  tree->Branch("qKurama",       event.qKurama,      "qKurama[ntKurama]/D");
  tree->Branch("m2",            event.m2,           "m2[ntKurama]/D");
  tree->Branch("xtgtKurama",    event.xtgtKurama,   "xtgtKurama[ntKurama]/D");
  tree->Branch("ytgtKurama",    event.ytgtKurama,   "ytgtKurama[ntKurama]/D");
  tree->Branch("utgtKurama",    event.utgtKurama,   "utgtKurama[ntKurama]/D");
  tree->Branch("vtgtKurama",    event.vtgtKurama,   "vtgtKurama[ntKurama]/D");
  tree->Branch("thetaKurama",   event.thetaKurama,  "thetaKurama[ntKurama]/D");
  tree->Branch("xtofKurama",    event.xtofKurama,   "xtofKurama[ntKurama]/D");
  tree->Branch("ytofKurama",    event.ytofKurama,   "ytofKurama[ntKurama]/D");
  tree->Branch("utofKurama",    event.utofKurama,   "utofKurama[ntKurama]/D");
  tree->Branch("vtofKurama",    event.vtofKurama,   "vtofKurama[ntKurama]/D");
  tree->Branch("tofsegKurama",  event.tofsegKurama, "tofsegKurama[ntKurama]/D");

  tree->Branch("rm_event",   &event.rm_event,   "rm_event/I");
  tree->Branch("rm_spill",   &event.rm_spill,   "rm_spill/I");
  tree->Branch("rm_accept",  &event.rm_accept,  "rm_accept/I");
  tree->Branch("rm_clear",   &event.rm_clear,   "rm_clear/I");
  tree->Branch("mst_accept", &event.mst_accept, "mst_accept/I");
  tree->Branch("mst_tdc",     event.mst_tdc,    Form( "mst_tdc[%d]/I", NumOfSegTOF ) );
  tree->Branch("mst_tof",     event.mst_tof,    Form( "mst_tof[%d]/I", NumOfSegTOF ) );
  tree->Branch("mst_sch",     event.mst_sch,    Form( "mst_sch[%d]/I", NumOfSegSCH ) );

  //Reaction
  tree->Branch("nKn",           &event.nKn,            "nKn/I");
  tree->Branch("nKp",           &event.nKp,            "nKp/I");
  tree->Branch("nKK",           &event.nKK,            "nKK/I");
  tree->Branch("vtx",            event.vtx,            "vtx[nKK]/D");
  tree->Branch("vty",            event.vty,            "vty[nKK]/D");
  tree->Branch("vtz",            event.vtz,            "vtz[nKK]/D");
  tree->Branch("closeDist",      event.closeDist,      "closeDist[nKK]/D");
  tree->Branch("theta",          event.theta,          "theta[nKK]/D");
  tree->Branch("MissP",          event.MissP,          "MissP[nKK]/D");
  tree->Branch("MissPu",         event.MissPu,         "MissPu[nKK]/D");
  tree->Branch("MissPv",         event.MissPv,         "MissPv[nKK]/D");
  tree->Branch("MissMass",       event.MissMass,       "MissMass[nKK]/D");
  tree->Branch("MissMassCorr",   event.MissMassCorr,   "MissMassCorr[nKK]/D");
  tree->Branch("MissMassCorrDE", event.MissMassCorrDE, "MissMassCorrDE[nKK]/D");
  tree->Branch("BE",      event.BE,       "BE[nKK]/D");
  tree->Branch("thetaCM", event.thetaCM,  "thetaCM[nKK]/D");
  tree->Branch("costCM",  event.costCM,   "costCM[nKK]/D");

  tree->Branch("xkn",        event.xkn,      "xkn[nKK]/D");
  tree->Branch("ykn",        event.ykn,      "ykn[nKK]/D");
  tree->Branch("ukn",        event.ukn,      "ukn[nKK]/D");
  tree->Branch("vkn",        event.vkn,      "vkn[nKK]/D");
  tree->Branch("xkp",        event.xkp,      "xkp[nKK]/D");
  tree->Branch("ykp",        event.ykp,      "ykp[nKK]/D");
  tree->Branch("ukp",        event.ukp,      "ukp[nKK]/D");
  tree->Branch("vkp",        event.vkp,      "vkp[nKK]/D");
  tree->Branch("pOrg",       event.pOrg,      "pOrg[nKK]/D");
  tree->Branch("pCalc",      event.pCalc,     "pCalc[nKK]/D");
  tree->Branch("pCorr",      event.pCorr,     "pCorr[nKK]/D");
  tree->Branch("pCorrDE",    event.pCorrDE,   "pCorrDE[nKK]/D");

  tree->Branch("xpos",  &event.xpos,  "xpos/D");
  tree->Branch("ypos",  &event.ypos,  "ypos/D");
  tree->Branch("state", &event.state, "state/I");

  ////////// Bring Address From Dst
  TTreeCont[kHodoscope]->SetBranchAddress("trigflag",  src.trigflag);
  TTreeCont[kHodoscope]->SetBranchAddress("trigpat",   src.trigpat);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh1",    &src.nhBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh1",     src.csBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh1Seg",    src.Bh1Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh1",      src.tBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh1",     src.dtBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh1",     src.deBh1);
  TTreeCont[kHodoscope]->SetBranchAddress("btof",      src.btof);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBh2",    &src.nhBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("csBh2",     src.csBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("Bh2Seg",    src.Bh2Seg);
  TTreeCont[kHodoscope]->SetBranchAddress("t0Bh2",     src.t0Bh2);
  TTreeCont[kHodoscope]->SetBranchAddress("tBh2",      src.tBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("dtBh2",     src.dtBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("deBh2",     src.deBh2);
  TTreeCont[kHodoscope]->SetBranchAddress("nhBac",    &src.nhBac);
  TTreeCont[kHodoscope]->SetBranchAddress("BacSeg",    src.BacSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tBac",      src.tBac);
  TTreeCont[kHodoscope]->SetBranchAddress("deBac",     src.deBac);
  TTreeCont[kHodoscope]->SetBranchAddress("nhPvac",   &src.nhPvac);
  TTreeCont[kHodoscope]->SetBranchAddress("PvacSeg",   src.PvacSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tPvac",     src.tPvac);
  TTreeCont[kHodoscope]->SetBranchAddress("dePvac",    src.dePvac);
  TTreeCont[kHodoscope]->SetBranchAddress("nhFac",    &src.nhFac);
  TTreeCont[kHodoscope]->SetBranchAddress("FacSeg",    src.FacSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tFac",      src.tFac);
  TTreeCont[kHodoscope]->SetBranchAddress("deFac",     src.deFac);
  TTreeCont[kHodoscope]->SetBranchAddress("nhTof",    &src.nhTof);
  TTreeCont[kHodoscope]->SetBranchAddress("csTof",     src.csTof);
  TTreeCont[kHodoscope]->SetBranchAddress("TofSeg",    src.TofSeg);
  TTreeCont[kHodoscope]->SetBranchAddress("tTof",      src.tTof);
  TTreeCont[kHodoscope]->SetBranchAddress("dtTof",     src.dtTof);
  TTreeCont[kHodoscope]->SetBranchAddress("deTof",     src.deTof);

  TTreeCont[kKuramaTracking]->SetBranchAddress("much",         &src.much         );
  TTreeCont[kKuramaTracking]->SetBranchAddress("ntSdcIn",      &src.ntSdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nlSdcIn",      &src.nlSdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhSdcIn",       src.nhSdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrSdcIn",   src.chisqrSdcIn  );
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0SdcIn",       src.x0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0SdcIn",       src.y0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0SdcIn",       src.u0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0SdcIn",       src.v0SdcIn      );
  TTreeCont[kKuramaTracking]->SetBranchAddress("ntSdcOut",     &src.ntSdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhSdcOut",      src.nhSdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrSdcOut",  src.chisqrSdcOut );
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0SdcOut",      src.x0SdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0SdcOut",      src.y0SdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0SdcOut",      src.u0SdcOut     );
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0SdcOut",      src.v0SdcOut     );

  TTreeCont[kKuramaTracking]->SetBranchAddress("deKaon",  src.deKaon  );
  TTreeCont[kKuramaTracking]->SetBranchAddress("ntXi",   &src.ntXi    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("x0Xi",    src.x0Xi    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("y0Xi",    src.y0Xi    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("u0Xi",    src.u0Xi    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("v0Xi",    src.v0Xi    );
  TTreeCont[kKuramaTracking]->SetBranchAddress("thetaXi", src.thetaXi );
  TTreeCont[kKuramaTracking]->SetBranchAddress("deXi",    src.deXi    );

  TTreeCont[kKuramaTracking]->SetBranchAddress("ntKurama",    &src.ntKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("nhKurama",     src.nhKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("stof",         src.stof);
  TTreeCont[kKuramaTracking]->SetBranchAddress("path",         src.path);
  TTreeCont[kKuramaTracking]->SetBranchAddress("pKurama",      src.pKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("qKurama",      src.qKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("chisqrKurama", src.chisqrKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtgtKurama",   src.xtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytgtKurama",   src.ytgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("utgtKurama",   src.utgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vtgtKurama",   src.vtgtKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("thetaKurama",  src.thetaKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("xtofKurama",   src.xtofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("ytofKurama",   src.ytofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("utofKurama",   src.utofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("vtofKurama",   src.vtofKurama);
  TTreeCont[kKuramaTracking]->SetBranchAddress("tofsegKurama", src.tofsegKurama);

  TTreeCont[kK18Tracking]->SetBranchAddress("ntBcOut",     &src.ntBcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("nlBcOut",     &src.nlBcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("nhBcOut",      src.nhBcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrBcOut",  src.chisqrBcOut );
  TTreeCont[kK18Tracking]->SetBranchAddress("x0BcOut",      src.x0BcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("y0BcOut",      src.y0BcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("u0BcOut",      src.u0BcOut     );
  TTreeCont[kK18Tracking]->SetBranchAddress("v0BcOut",      src.v0BcOut     );

  TTreeCont[kK18Tracking]->SetBranchAddress("ntK18", &src.ntK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("nhK18", &src.nhK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("chisqrK18", &src.chisqrK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("p_3rd", &src.pK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("xtgtK18", &src.xtgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("ytgtK18", &src.ytgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("utgtK18", &src.utgtK18);
  TTreeCont[kK18Tracking]->SetBranchAddress("vtgtK18", &src.vtgtK18);

  TTreeCont[kEasiroc]->SetBranchAddress("bft_ncl",    &src.nhBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_clsize", &src.csBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_ctime",  &src.tBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_ctot",   &src.wBft);
  TTreeCont[kEasiroc]->SetBranchAddress("bft_clpos",  &src.BftPos);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ncl",    &src.nhSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_clsize", &src.csSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ctime",  &src.tSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_ctot",   &src.wSch);
  TTreeCont[kEasiroc]->SetBranchAddress("sch_clpos",  &src.SchPos);
  TTreeCont[kEasiroc]->SetBranchAddress("fbh_ncl",    &src.nhFbh);
  TTreeCont[kEasiroc]->SetBranchAddress("fbh_clsize", &src.csFbh);
  TTreeCont[kEasiroc]->SetBranchAddress("fbh_ctime",  &src.tFbh);
  TTreeCont[kEasiroc]->SetBranchAddress("fbh_ctot",   &src.wFbh);
  TTreeCont[kEasiroc]->SetBranchAddress("fbh_clpos",  &src.FbhPos);

  TTreeCont[kEMC]->SetBranchAddress("runnum", &src.runnum);
  TTreeCont[kEMC]->SetBranchAddress("evnum",  &src.evnum);
  TTreeCont[kEMC]->SetBranchAddress("spill",  &src.spill);
  TTreeCont[kEMC]->SetBranchAddress("xpos",   &src.xpos);
  TTreeCont[kEMC]->SetBranchAddress("ypos",   &src.ypos);
  TTreeCont[kEMC]->SetBranchAddress("state",  &src.state);

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
