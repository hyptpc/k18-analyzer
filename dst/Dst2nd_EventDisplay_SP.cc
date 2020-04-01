/**
 *  file: DstSkeleton.cc
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
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "EventDisplayDST.hh"
#include "DstHelper.hh"

namespace
{
  using namespace root;
  using namespace dst;
  const std::string& class_name("Dst2nd_EventDisplay");
  ConfMan&            gConf = ConfMan::GetInstance();
  const DCGeomMan&    gGeom = DCGeomMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
  EventDisplayDST&        gEvDisp = EventDisplayDST::GetInstance();
}

//functions for calculation
ThreeVector calcStartPoint(ThreeVector pos, ThreeVector dir, ThreeVector vert);
ThreeVector calcEndPoint(ThreeVector pos, ThreeVector dir, int mode);

namespace dst
{
  enum kArgc
    {
      kProcess, kConfFile,
      kDstPiKana,
      nArgc
    };
  std::vector<TString> ArgName =
    { "[Process]", "[ConfFile]",
      "[DSTPiKana]"};
  std::vector<TString> TreeName =
    { "", "", "pik"};
  std::vector<TFile*> TFileCont;
  std::vector<TTree*> TTreeCont;
}

//_____________________________________________________________________
struct Event
{
  int runnum;
  int evnum;
  int spill;
  std::string runno;
};

//_____________________________________________________________________
struct Src
{
  int runnum;
  int evnum;
  int spill;

  //K18 track
  int bft_ncl;
  int ntK18;
  double pK18[MaxHits];
  double xtgtK18[MaxHits];
  double ytgtK18[MaxHits];
  double utgtK18[MaxHits];
  double vtgtK18[MaxHits];
  int priority_K18[MaxHits];

  //KURAMA
  int    ntKurama;
  double chisqrKurama[MaxHits];
  double pKurama[MaxHits];
  double m2[MaxHits];
  double cm2[MaxHits];
  double xtgtKurama[MaxHits];
  double ytgtKurama[MaxHits];
  double utgtKurama[MaxHits];
  double vtgtKurama[MaxHits];

  //CFT
  int ntCFT;
  int ntProton;
  int ntOther;
  int tracknoproton[MaxHits];
  int tracknoother[MaxHits];
  double theta_cft[MaxHits];
  double PosT0_x[MaxHits];
  double PosT0_y[MaxHits];
  double PosT0_z[MaxHits];
  double Dir_x[MaxHits];
  double Dir_y[MaxHits];
  double Dir_z[MaxHits];
  int segBGOt[MaxHits];
  double energyBGOt[MaxHits];
  double Total_dE[MaxHits];
  int segPiIDt[MaxHits];
  int protonflagt[MaxHits];
  int BGOnohitt[MaxHits];
  int simKuramat[MaxHits];
  int simKuramapartner[MaxHits];
  int pipelasticflag[MaxHits];
  int PiIDflagt[MaxHits];

  double energyBGO[NumOfSegBGO];
  int PiIDflag[NumOfSegPiID];

  //Reaction
  int    nPi;
  int    nK;
  int    nPiK;
  double vtx_KURAMA[MaxHits];
  double vty_KURAMA[MaxHits];
  double vtz_KURAMA[MaxHits];
  double ccm2[MaxHits];
  double chisqrKuramapik[MaxHits];
  double pKuramapik[MaxHits];
  double closedist_KURAMA[MaxHits];
  int priority_K18pik[MaxHits];
  double theta[MaxHits];
  double MissMass[MaxHits];
  double MissMassCorr[MaxHits];
  double MissMassCorrDE[MaxHits];
  double MissMassPiPi[MaxHits];
  double MissMassPiP[MaxHits];
  double thetaCM[MaxHits];
  double costCM[MaxHits];

  double MissMom[MaxHits];
  double MissMomx[MaxHits];
  double MissMomy[MaxHits];
  double MissMomz[MaxHits];
  double MissMomcal[MaxHits];
  double MissMomxcal[MaxHits];
  double MissMomycal[MaxHits];
  double MissMomzcal[MaxHits];
  double SigmaBeta[MaxHits];

  double xpi[MaxHits];
  double ypi[MaxHits];
  double upi[MaxHits];
  double vpi[MaxHits];
  double xk[MaxHits];
  double yk[MaxHits];
  double uk[MaxHits];
  double vk[MaxHits];
  double uc[MaxHits];
  double vc[MaxHits];

  double pOrg[MaxHits];
  double pCalc[MaxHits];
  double pCorr[MaxHits];
  double pCorrDE[MaxHits];

  int KURAMAPID[MaxHits];//1:K+ 2:proton 3:pi+ 4:pi-

  //combination
  int nCatch;
  int nPiKCatch;
  double vtx_K18Catch[MaxHits];
  double vty_K18Catch[MaxHits];
  double vtz_K18Catch[MaxHits];
  
  //NPScat
  double DeltaE_NPScat;
  double ProtonMom_NPScat;
  double DecayNeutronMom;
  double DecayPionMom;
  double MissMassSigmaP_NPScat;
  double vtx_Decay2np;
  double vty_Decay2np;
  double vtz_Decay2np;
  double vtx_NPScat;
  double vty_NPScat;
  double vtz_NPScat;
  double cdistDecay2np;
  double cdistNPScat;
  double theta_NPScat;
  double thetaCM_NPScat;
  int ptrno_NPScat;
  int pitrno_NPScat;
  int pikno_NPScat;
  int priority_NPScat;

  //decayp assumption
  double DeltaE_DecayP;
  double ProtonMom_DecayP;
  double MissMassSigmaP_DecayP;
  double cdistDecayP;
  double SigmaMomcor_DecayP;
  double SigmaLength;
  double SigmaLengthcor;
  double theta_DecayP;
  int ptrno_DecayP;
  int pikno_DecayP;
  int priority_DecayP;

  //PiPScat assumption
  double DeltaP_PiPScat;
  double DeltaPD_PiPScat;
  double DeltaPD_PiPScatM;
  double DecayNeutronMom2;
  double DecayPionMom2;
  double vtx_Decay2pip;
  double vty_Decay2pip;
  double vtz_Decay2pip;
  double vtx_PiPScat;
  double vty_PiPScat;
  double vtz_PiPScat;
  double cdistDecay2pip;
  double cdistPiPScat;
  int ptrno_PiPScat;
  int pitrno_PiPScat;
  int pikno_PiPScat;
  int priority_PiPScat;

  //PPScat assumption
  double DeltaP_PPScat;
  double theta_PPScat;
  double vtx_Decay2pp;
  double vty_Decay2pp;
  double vtz_Decay2pp;
  double vtx_PPScat;
  double vty_PPScat;
  double vtz_PPScat;
  double cdistDecay2pp;
  double cdistPPScat;
  int ptrno_PPScat;
  int p2trno_PPScat;
  int pikno_PPScat;
  int priority_PPScat;

  //sigma p scattering assumption
  double DeltaE_SigmaPScat;
  double MissMassSigmaP_SigmaPScat;
  double vtx_SigmaPScat;
  double vty_SigmaPScat;
  double vtz_SigmaPScat;
  double cdistSigmaPScat;
  int pikno_SigmaPScat;
  int ptrno_SigmaPScat;

  //npi+ decay mode
  double DeltaE_SigmaPScat2npi;
  double DeltaE2_SigmaPScat2npi;
  double ProtonMom_SigmaPScat2npi;
  double MissMassSigmaP_SigmaPScat2npi;
  double vtx_ScatSigmaDecay2npi;
  double vty_ScatSigmaDecay2npi;
  double vtz_ScatSigmaDecay2npi;
  double cdistScatSigmaDecay2npi;
  double vdistance_SigmaPScat2npi;
  double theta_SigmaPScat2npi;
  double thetaCM_SigmaPScat2npi;
  double dcostheta_SigmaL_SigmaPScat2npi;
  double dcostheta_SigmaL2_SigmaPScat2npi;
  int pikno_SigmaPScat2npi;
  int ptrno_SigmaPScat2npi;
  int pitrno_SigmaPScat2npi;
  int priority_SigmaPScat2npi;

  //ppi decay mode
  double DeltaE_SigmaPScat2ppi;
  double DeltaE2_SigmaPScat2ppi;
  double ProtonMom_SigmaPScat2ppi;
  double theta2proton_SigmaPScat2ppi;
  double MissMassSigmaP_SigmaPScat2ppi;
  double MissMassSigmaP2_SigmaPScat2ppi;
  double vtx_ScatSigmaDecay2ppi;
  double vty_ScatSigmaDecay2ppi;
  double vtz_ScatSigmaDecay2ppi;
  double cdistScatSigmaDecay2ppi;
  double vdistance_SigmaPScat2ppi;
  double theta_SigmaPScat2ppi;
  double thetaCM_SigmaPScat2ppi;
  int pikno_SigmaPScat2ppi;
  int ptrno_SigmaPScat2ppi;
  int p2trno_SigmaPScat2ppi;
  int priority_SigmaPScat2ppi;
  double DeltaP_PPScatSP;
  double theta_PPScatSP;

};

//_____________________________________________________________________
namespace root
{
  Event  event;
  Src    src;
  TH1   *h[MaxHist];
  //TTree *tree;
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

  //std::string runno=arg[kDstPiKana].substr(arg[kDstPiKana].length()-26,5);
  std::string runno="hoge";
  std::cout<<runno<<std::endl;

  event.runno=runno;

  int ievent = 0;
  for( ; ievent<nevent && !CatchSignal::Stop(); ++ievent ){
    InitializeEvent();
    if( DstRead( ievent )){
      gEvDisp.EndOfEvent();//tree->Fill();
    } 
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
  return true;
}

//_____________________________________________________________________
bool
dst::DstOpen( std::vector<std::string> arg )
{
  int open_file = 0;
  int open_tree = 0;
  for( std::size_t i=0; i<nArgc; ++i ){
    if( i==kProcess || i==kConfFile ) continue;
    open_file += OpenFile( TFileCont[i], arg[i] );
    open_tree += OpenTree( TFileCont[i], TTreeCont[i], TreeName[i] );
  }

  if( open_file!=open_tree || open_file!=nArgc-2 )
    return false;
  if( !CheckEntries( TTreeCont ) )
    return false;

  //  TFileCont[kOutFile] = new TFile( arg[kOutFile].c_str(), "recreate" );

  return true;
}

//_____________________________________________________________________
bool
dst::DstRead( int ievent )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");

  if( ievent%10000==0 ){
    std::cout << "#D " << func_name << " Event Number: "
	      << std::setw(6) << ievent << std::endl;
  }

  GetEntry(ievent);

  HF1( 1, 0. );

  event.runnum = src.runnum;
  event.evnum  = src.evnum;
  event.spill  = src.spill;

  bool cutcondition = false;

  //sigmapscat2npi mode
  if(src.pikno_SigmaPScat2npi>-1 && src.MissMass[src.pikno_SigmaPScat2npi]<1.25 && src.MissMass[src.pikno_SigmaPScat2npi]>1.12 && src.ntProton==1 && src.ntOther==1 && src.cdistSigmaPScat<25 && abs(src.vtx_SigmaPScat)<30 && abs(src.vty_SigmaPScat)<30 && abs(src.vtz_SigmaPScat)<200 && abs(src.vtx_KURAMA[src.pikno_SigmaPScat2npi])<25 && abs(src.vty_KURAMA[src.pikno_SigmaPScat2npi])<25 && abs(src.vtz_KURAMA[src.pikno_SigmaPScat2npi])<200 && src.cdistScatSigmaDecay2npi<25 && abs(src.vtx_ScatSigmaDecay2npi)<30 && abs(src.vty_ScatSigmaDecay2npi)<30 && abs(src.vtz_ScatSigmaDecay2npi)<160 && pow(src.vtx_KURAMA[src.pikno_SigmaPScat2npi],2)+pow(src.vty_KURAMA[src.pikno_SigmaPScat2npi],2)<18*18 && src.vtx_KURAMA[src.pikno_SigmaPScat2npi]<14 && (src.vtz_ScatSigmaDecay2npi-src.vtz_SigmaPScat<30 && src.vtz_ScatSigmaDecay2npi-src.vtz_SigmaPScat>-30 && src.vtz_ScatSigmaDecay2npi-src.vtz_KURAMA[src.pikno_SigmaPScat2npi]<50 && src.vtz_ScatSigmaDecay2npi-src.vtz_KURAMA[src.pikno_SigmaPScat2npi]>-50) && !(src.MissMassSigmaP_SigmaPScat2npi<0.035 && src.MissMassSigmaP_SigmaPScat2npi>0.005) && !(src.DeltaE_NPScat<10 && src.DeltaE_NPScat>-15) && -src.DeltaP_PiPScat>0.015 && src.pipelasticflag[src.ptrno_SigmaPScat2npi]==0 && cos(src.thetaCM_SigmaPScat2npi*3.141592/180)<0.5 && cos(src.thetaCM_SigmaPScat2npi*3.141592/180)>-0.6 && src.dcostheta_SigmaL_SigmaPScat2npi>-10 && src.dcostheta_SigmaL2_SigmaPScat2npi>-10){
    cutcondition=true;
  }

  if(cutcondition ==false){
    return true;
  }

  
  gEvDisp.DrawText( 0.1, 0.2, Form("Run# %5s Event# %6d",
				   event.runno.c_str(),
				   event.evnum) );
  
  HF1( 1, 1. );

  for(int it=0;it<NumOfSegBGO;it++){
    if(src.energyBGO[it]<0) continue;
    gEvDisp.ShowHitBGO(it,src.energyBGO[it]);
  }
  /*
  for(int it=0;it<NumOfSegPiID;it++){
    if(src.PiIDflag[it]!=it) continue;
    gEvDisp.ShowHitPiID(it);
  }
  */

  for(int it=0;it<src.ntCFT;it++){
    if(src.PiIDflagt[it]>0)  gEvDisp.ShowHitPiID(src.segPiIDt[it]);
  }

  //for sigmapscat2npi
  if(src.pikno_SigmaPScat2npi>-1){
    int selectedbno=src.pikno_SigmaPScat2npi/src.nPi;
    int selectedkno=src.pikno_SigmaPScat2npi%src.nPi;
    ThreeVector vert_KURAMA(src.vtx_KURAMA[src.pikno_SigmaPScat2npi],src.vty_KURAMA[src.pikno_SigmaPScat2npi],src.vtz_KURAMA[src.pikno_SigmaPScat2npi]);
    gEvDisp.SetVertex(vert_KURAMA);
    ThreeVector vert_Scat(src.vtx_SigmaPScat,src.vty_SigmaPScat,src.vtz_SigmaPScat);
    gEvDisp.SetVertex(vert_Scat);
    ThreeVector vert_Decay(src.vtx_ScatSigmaDecay2npi,src.vty_ScatSigmaDecay2npi, src.vtz_ScatSigmaDecay2npi);
    gEvDisp.SetVertex(vert_Decay);

    for(int it=0;it<src.ntK18;it++){
      if(it!=selectedkno){
	ThreeVector bX0(src.xtgtK18[it]-220*src.utgtK18[it],src.ytgtK18[it]-220*src.vtgtK18[it],-220);
	ThreeVector bX1(src.xtgtK18[it]-200*src.utgtK18[it],src.ytgtK18[it]-200*src.vtgtK18[it],-200);
	gEvDisp.SetK18Track(bX0,bX1,0);
      }else{
	ThreeVector bX0(src.xtgtK18[it]-220*src.utgtK18[it],src.ytgtK18[it]-220*src.vtgtK18[it],-220);
	ThreeVector input1(src.xtgtK18[it],src.ytgtK18[it],0);
	ThreeVector input2(src.utgtK18[it],src.vtgtK18[it],1);
	ThreeVector bX1=calcStartPoint(input1,input2,vert_Scat);
	gEvDisp.SetK18Track(bX0,bX1,1);
      }
    }

    for(int it=0;it<src.ntKurama;it++){
      if(it!=selectedbno){
	ThreeVector bX0(src.xtgtKurama[it]+150*src.utgtKurama[it],src.ytgtKurama[it]+150*src.vtgtKurama[it],150);
	ThreeVector bX1(src.xtgtKurama[it]+420*src.utgtKurama[it],src.ytgtKurama[it]+420*src.vtgtKurama[it],420);
	gEvDisp.SetKuramaTrack(bX0,bX1,0);
      }else{
	ThreeVector bX1(src.xtgtKurama[it]+420*src.utgtKurama[it],src.ytgtKurama[it]+420*src.vtgtKurama[it],420);
	ThreeVector input1(src.xtgtKurama[it],src.ytgtKurama[it],0);
	ThreeVector input2(src.utgtKurama[it],src.vtgtKurama[it],1);
	ThreeVector bX0=calcStartPoint(input1,input2,vert_Scat);
	gEvDisp.SetKuramaTrack(bX0,bX1,1);
      }
    }

    ThreeVector sigmavec(src.MissMomxcal[src.pikno_SigmaPScat2npi],src.MissMomycal[src.pikno_SigmaPScat2npi],src.MissMomzcal[src.pikno_SigmaPScat2npi]);
    ThreeVector sigmaX1=calcStartPoint(vert_KURAMA,sigmavec,vert_Scat);
    gEvDisp.SetSigmaTrack(vert_KURAMA,sigmaX1,0);
    gEvDisp.SetSigmaTrack(vert_Scat,vert_Decay,1);

    ThreeVector ppos(src.PosT0_x[src.ptrno_SigmaPScat2npi],src.PosT0_y[src.ptrno_SigmaPScat2npi],src.PosT0_z[src.ptrno_SigmaPScat2npi]);
    ThreeVector pdir(src.Dir_x[src.ptrno_SigmaPScat2npi],src.Dir_y[src.ptrno_SigmaPScat2npi],src.Dir_z[src.ptrno_SigmaPScat2npi]);

    ThreeVector pX0=calcStartPoint(ppos,pdir,vert_Scat);
    int pmode=1;
    if(src.energyBGOt[src.ptrno_SigmaPScat2npi]==0){
      pmode=0;
      if(src.BGOnohitt[src.ptrno_SigmaPScat2npi]==1){
	pmode=2;
      }
    }
    ThreeVector pX1=calcEndPoint(ppos,pdir,pmode);

    gEvDisp.SetCFTTrack(pX0,pX1, 1, 1);
    std::cout<<"SetProton : energy="<<src.energyBGOt[src.ptrno_SigmaPScat2npi]+src.Total_dE[src.ptrno_SigmaPScat2npi]<<std::endl;

    ThreeVector pipos(src.PosT0_x[src.pitrno_SigmaPScat2npi],src.PosT0_y[src.pitrno_SigmaPScat2npi],src.PosT0_z[src.pitrno_SigmaPScat2npi]);
    ThreeVector pidir(src.Dir_x[src.pitrno_SigmaPScat2npi],src.Dir_y[src.pitrno_SigmaPScat2npi],src.Dir_z[src.pitrno_SigmaPScat2npi]);
    ThreeVector piX0=calcStartPoint(pipos,pidir,vert_Decay);
    int pimode=2;
    if(src.energyBGOt[src.pitrno_SigmaPScat2npi]==0){
      pimode=0;
      if(src.BGOnohitt[src.pitrno_SigmaPScat2npi]==1){
	pimode=2;
      }
    }else if(src.PiIDflagt[src.pitrno_SigmaPScat2npi]==0){
      pimode=1;
    }
    ThreeVector piX1=calcEndPoint(pipos,pidir,pimode);

    gEvDisp.SetCFTTrack(piX0,piX1, 0,1);

    //    std::cout<<"DeltaE_SigmaPScat2npi="<<src.DeltaE_SigmaPScat2npi<<std::endl;
    gEvDisp.DrawText( 0.6, 0.2, Form("#DeltaE(#Sigmap Scattering to n#pi^{+})=%.2f",
				     src.DeltaE_SigmaPScat2npi) );

  }

    gEvDisp.UpdateCATCH();

  if( 1 )
    gEvDisp.GetCommand();

  return true;
}

//_____________________________________________________________________
ThreeVector calcStartPoint(ThreeVector pos, ThreeVector dir, ThreeVector vert)
{
  //vert:
  ThreeVector dpos=pos-vert;
  double t=-1*(dpos*dir)/dir.Mag2();

  return pos+t*dir;

}
//_____________________________________________________________________
ThreeVector calcEndPoint(ThreeVector pos, ThreeVector dir, int mode)
{
  //mode 0:stopped in CFT, 1:stopped in BGO 2:penetlated PIID(other)
  double r;
  
  switch(mode){
  case 0:
  r=80.0;
  break;
  case 1:
  r=120.0;
  break;
  case 2:
  r=200;
  }

  double px=pos.x();
  double py=pos.y();
  double dx=dir.x();
  double dy=dir.y();

  double D=pow(dx*px+dy*py,2)-(dx*dx+dy*dy)*(px*px+py*py-r*r);
  
  if(D<0){
    return pos+3*dir;
  }

  double t=(-(dx*px+dy*py)+sqrt(D))/(dx*dx+dy*dy);

  return pos+t*dir;

}


bool
dst::DstClose( void )
{
  //TFileCont[kOutFile]->Write();
  //std::cout << "#D Close : " << TFileCont[kOutFile]->GetName() << std::endl;
  // TFileCont[kOutFile]->Close();

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

  ////////////////////////////////////////////
  //Tree
  /*
  HBTree("dst","tree of DstSkeleton");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  */
  ////////// Bring Address From Dst
  TTreeCont[kDstPiKana]->SetBranchStatus("*", 1);

  TTreeCont[kDstPiKana]->SetBranchAddress("runnum", &src.runnum);
  TTreeCont[kDstPiKana]->SetBranchAddress("evnum",  &src.evnum);
  TTreeCont[kDstPiKana]->SetBranchAddress("spill",  &src.spill);
  TTreeCont[kDstPiKana]->SetBranchAddress("bft_ncl", &src.bft_ncl);
  TTreeCont[kDstPiKana]->SetBranchAddress("ntK18", &src.ntK18);
  TTreeCont[kDstPiKana]->SetBranchAddress("pK18", &src.pK18);
  TTreeCont[kDstPiKana]->SetBranchAddress("xtgtK18", &src.xtgtK18);
  TTreeCont[kDstPiKana]->SetBranchAddress("ytgtK18", &src.ytgtK18);
  TTreeCont[kDstPiKana]->SetBranchAddress("utgtK18", &src.utgtK18);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtgtK18", &src.vtgtK18);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_K18", &src.priority_K18);

  TTreeCont[kDstPiKana]->SetBranchAddress("ntKurama", &src.ntKurama);
  TTreeCont[kDstPiKana]->SetBranchAddress("chisqrKurama", &src.chisqrKurama);
  TTreeCont[kDstPiKana]->SetBranchAddress("pKurama", &src.pKurama);
  TTreeCont[kDstPiKana]->SetBranchAddress("m2", &src.m2);
  TTreeCont[kDstPiKana]->SetBranchAddress("cm2", &src.cm2);
  TTreeCont[kDstPiKana]->SetBranchAddress("xtgtKurama", &src.xtgtKurama);
  TTreeCont[kDstPiKana]->SetBranchAddress("ytgtKurama", &src.ytgtKurama);
  TTreeCont[kDstPiKana]->SetBranchAddress("utgtKurama", &src.utgtKurama);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtgtKurama", &src.vtgtKurama);
  //CFT
  TTreeCont[kDstPiKana]->SetBranchAddress("ntCFT", &src.ntCFT);
  TTreeCont[kDstPiKana]->SetBranchAddress("ntProton", &src.ntProton);
  TTreeCont[kDstPiKana]->SetBranchAddress("ntOther", &src.ntOther);
  TTreeCont[kDstPiKana]->SetBranchAddress("tracknoproton", &src.tracknoproton);
  TTreeCont[kDstPiKana]->SetBranchAddress("tracknoother", &src.tracknoother);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta_cft", &src.theta_cft);
  TTreeCont[kDstPiKana]->SetBranchAddress("PosT0_x", &src.PosT0_x);
  TTreeCont[kDstPiKana]->SetBranchAddress("PosT0_y", &src.PosT0_y);
  TTreeCont[kDstPiKana]->SetBranchAddress("PosT0_z", &src.PosT0_z);
  TTreeCont[kDstPiKana]->SetBranchAddress("Dir_x", &src.Dir_x);
  TTreeCont[kDstPiKana]->SetBranchAddress("Dir_y", &src.Dir_y);
  TTreeCont[kDstPiKana]->SetBranchAddress("Dir_z", &src.Dir_z);  
  TTreeCont[kDstPiKana]->SetBranchAddress("segBGOt", &src.segBGOt);
  TTreeCont[kDstPiKana]->SetBranchAddress("energyBGOt", &src.energyBGOt);
  TTreeCont[kDstPiKana]->SetBranchAddress("Total_dE", &src.Total_dE);
  TTreeCont[kDstPiKana]->SetBranchAddress("segPiIDt", &src.segPiIDt);
  TTreeCont[kDstPiKana]->SetBranchAddress("protonflagt", &src.protonflagt);
  TTreeCont[kDstPiKana]->SetBranchAddress("BGOnohitt", &src.BGOnohitt);
  TTreeCont[kDstPiKana]->SetBranchAddress("simKuramat", &src.simKuramat);
  TTreeCont[kDstPiKana]->SetBranchAddress("simKuramapartner", &src.simKuramapartner);
  TTreeCont[kDstPiKana]->SetBranchAddress("pipelasticflag", &src.pipelasticflag);
  TTreeCont[kDstPiKana]->SetBranchAddress("energyBGO", &src.energyBGO);
  TTreeCont[kDstPiKana]->SetBranchAddress("PiIDflag", &src.PiIDflag);
  TTreeCont[kDstPiKana]->SetBranchAddress("PiIDflagt", &src.PiIDflagt);

  TTreeCont[kDstPiKana]->SetBranchAddress("nPi", &src.nPi);
  TTreeCont[kDstPiKana]->SetBranchAddress("nK", &src.nK);
  TTreeCont[kDstPiKana]->SetBranchAddress("nPiK", &src.nPiK);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_KURAMA", &src.vtx_KURAMA);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_KURAMA", &src.vty_KURAMA);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_KURAMA", &src.vtz_KURAMA);
  TTreeCont[kDstPiKana]->SetBranchAddress("ccm2", &src.ccm2);
  TTreeCont[kDstPiKana]->SetBranchAddress("chisqrKuramapik", &src.chisqrKuramapik);
  TTreeCont[kDstPiKana]->SetBranchAddress("pKuramapik", &src.pKuramapik);
  TTreeCont[kDstPiKana]->SetBranchAddress("closedist_KURAMA", &src.closedist_KURAMA);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_K18pik", &src.priority_K18pik);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta", &src.theta);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMass", &src.MissMass);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassCorr", &src.MissMassCorr);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassCorrDE", &src.MissMassCorrDE);

  TTreeCont[kDstPiKana]->SetBranchAddress("MissMomcal", &src.MissMomcal);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMomxcal", &src.MissMomxcal);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMomycal", &src.MissMomycal);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMomzcal", &src.MissMomzcal);

  TTreeCont[kDstPiKana]->SetBranchAddress("xpi", &src.xpi);
  TTreeCont[kDstPiKana]->SetBranchAddress("ypi", &src.ypi);
  TTreeCont[kDstPiKana]->SetBranchAddress("upi", &src.upi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vpi", &src.vpi);
  TTreeCont[kDstPiKana]->SetBranchAddress("xk", &src.xk);
  TTreeCont[kDstPiKana]->SetBranchAddress("yk", &src.yk);
  TTreeCont[kDstPiKana]->SetBranchAddress("uk", &src.uk);
  TTreeCont[kDstPiKana]->SetBranchAddress("vk", &src.vk);
  TTreeCont[kDstPiKana]->SetBranchAddress("uc", &src.uc);
  TTreeCont[kDstPiKana]->SetBranchAddress("vc", &src.vc);
  TTreeCont[kDstPiKana]->SetBranchAddress("pOrg", &src.pOrg);
  TTreeCont[kDstPiKana]->SetBranchAddress("pCalc", &src.pCalc);
  TTreeCont[kDstPiKana]->SetBranchAddress("pCorr", &src.pCorr);
  TTreeCont[kDstPiKana]->SetBranchAddress("pCorrDE", &src.pCorrDE);

  TTreeCont[kDstPiKana]->SetBranchAddress("KURAMAPID", &src.KURAMAPID);

  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaE_NPScat", &src.DeltaE_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("ProtonMom_NPScat", &src.ProtonMom_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("DecayNeutronMom", &src.DecayNeutronMom);
  TTreeCont[kDstPiKana]->SetBranchAddress("DecayPionMom", &src.DecayPionMom);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassSigmaP_NPScat", &src.MissMassSigmaP_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_Decay2np", &src.vtx_Decay2np);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_Decay2np", &src.vty_Decay2np);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_Decay2np", &src.vtz_Decay2np);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_NPScat", &src.vtx_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_NPScat", &src.vty_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_NPScat", &src.vtz_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistDecay2np", &src.cdistDecay2np);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistNPScat", &src.cdistNPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta_NPScat", &src.theta_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("thetaCM_NPScat", &src.thetaCM_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("ptrno_NPScat", &src.ptrno_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("pitrno_NPScat", &src.pitrno_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("pikno_NPScat", &src.pikno_NPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_NPScat", &src.priority_NPScat);

  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaE_DecayP", &src.DeltaE_DecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("ProtonMom_DecayP", &src.ProtonMom_DecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassSigmaP_DecayP", &src.MissMassSigmaP_DecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistDecayP", &src.cdistDecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("SigmaMomcor_DecayP", &src.SigmaMomcor_DecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("SigmaLength", &src.SigmaLength);
  TTreeCont[kDstPiKana]->SetBranchAddress("SigmaLengthcor", &src.SigmaLengthcor);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta_DecayP", &src.theta_DecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("ptrno_DecayP", &src.ptrno_DecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("pikno_DecayP", &src.pikno_DecayP);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_DecayP", &src.priority_DecayP);

  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaP_PiPScat", &src.DeltaP_PiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("DecayNeutronMom2", &src.DecayNeutronMom2);
  TTreeCont[kDstPiKana]->SetBranchAddress("DecayPionMom2", &src.DecayPionMom2);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_Decay2pip", &src.vtx_Decay2pip);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_Decay2pip", &src.vty_Decay2pip);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_Decay2pip", &src.vtz_Decay2pip);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_PiPScat", &src.vtx_PiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_PiPScat", &src.vty_PiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_PiPScat", &src.vtz_PiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistDecay2pip", &src.cdistDecay2pip);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistPiPScat", &src.cdistPiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("ptrno_PiPScat", &src.ptrno_PiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("pitrno_PiPScat", &src.pitrno_PiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("pikno_PiPScat", &src.pikno_PiPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_PiPScat", &src.priority_PiPScat);

  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaP_PPScat", &src.DeltaP_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta_PPScat", &src.theta_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_Decay2pp", &src.vtx_Decay2pp);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_Decay2pp", &src.vty_Decay2pp);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_Decay2pp", &src.vtz_Decay2pp);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_PPScat", &src.vtx_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_PPScat", &src.vty_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_PPScat", &src.vtz_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistDecay2pp", &src.cdistDecay2pp);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistPPScat", &src.cdistPPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("ptrno_PPScat", &src.ptrno_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("p2trno_PPScat", &src.p2trno_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("pikno_PPScat", &src.pikno_PPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_PPScat", &src.priority_PPScat);

  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaE_SigmaPScat", &src.DeltaE_SigmaPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassSigmaP_SigmaPScat", &src.MissMassSigmaP_SigmaPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_SigmaPScat", &src.vtx_SigmaPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_SigmaPScat", &src.vty_SigmaPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_SigmaPScat", &src.vtz_SigmaPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistSigmaPScat", &src.cdistSigmaPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("pikno_SigmaPScat", &src.pikno_SigmaPScat);
  TTreeCont[kDstPiKana]->SetBranchAddress("ptrno_SigmaPScat", &src.ptrno_SigmaPScat);

  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaE_SigmaPScat2npi", &src.DeltaE_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaE2_SigmaPScat2npi", &src.DeltaE2_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("ProtonMom_SigmaPScat2npi", &src.ProtonMom_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassSigmaP_SigmaPScat2npi", &src.MissMassSigmaP_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_ScatSigmaDecay2npi", &src.vtx_ScatSigmaDecay2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_ScatSigmaDecay2npi", &src.vty_ScatSigmaDecay2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_ScatSigmaDecay2npi", &src.vtz_ScatSigmaDecay2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistScatSigmaDecay2npi", &src.cdistScatSigmaDecay2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vdistance_SigmaPScat2npi", &src.vdistance_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta_SigmaPScat2npi", &src.theta_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("thetaCM_SigmaPScat2npi", &src.thetaCM_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("dcostheta_SigmaL_SigmaPScat2npi", &src.dcostheta_SigmaL_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("dcostheta_SigmaL2_SigmaPScat2npi", &src.dcostheta_SigmaL2_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("pikno_SigmaPScat2npi", &src.pikno_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("ptrno_SigmaPScat2npi", &src.ptrno_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("pitrno_SigmaPScat2npi", &src.pitrno_SigmaPScat2npi);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_SigmaPScat2npi", &src.priority_SigmaPScat2npi);

  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaE_SigmaPScat2ppi", &src.DeltaE_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaE2_SigmaPScat2ppi", &src.DeltaE2_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("ProtonMom_SigmaPScat2ppi", &src.ProtonMom_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta2proton_SigmaPScat2ppi", &src.theta2proton_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassSigmaP_SigmaPScat2ppi", &src.MissMassSigmaP_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("MissMassSigmaP2_SigmaPScat2ppi", &src.MissMassSigmaP2_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtx_ScatSigmaDecay2ppi", &src.vtx_ScatSigmaDecay2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vty_ScatSigmaDecay2ppi", &src.vty_ScatSigmaDecay2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vtz_ScatSigmaDecay2ppi", &src.vtz_ScatSigmaDecay2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("cdistScatSigmaDecay2ppi", &src.cdistScatSigmaDecay2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("vdistance_SigmaPScat2ppi", &src.vdistance_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta_SigmaPScat2ppi", &src.theta_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("thetaCM_SigmaPScat2ppi", &src.thetaCM_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("pikno_SigmaPScat2ppi", &src.pikno_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("ptrno_SigmaPScat2ppi", &src.ptrno_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("p2trno_SigmaPScat2ppi", &src.p2trno_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("priority_SigmaPScat2ppi", &src.priority_SigmaPScat2ppi);
  TTreeCont[kDstPiKana]->SetBranchAddress("DeltaP_PPScatSP", &src.DeltaP_PPScatSP);
  TTreeCont[kDstPiKana]->SetBranchAddress("theta_PPScatSP", &src.theta_PPScatSP);


  return true;
}

//_____________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")   &&
      InitializeParameter<UserParamMan>("USER") &&
      InitializeParameter<EventDisplayDST>());
}

//_____________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
