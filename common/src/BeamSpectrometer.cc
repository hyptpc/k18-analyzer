// BeamSpectrometer.cpp
#include <iostream>
#include "BeamSpectrometer.hh"
#include "TransferMatrixMan.hh"
#include "BLDCWireMapMan.hh"

#define DEBUG 0
#define DEBUG2 0

// parameters for TMinuit
static const Double_t  FitStep[5] = { 1e-3, 1e-3, 1e-3,1e-3,1e-3 };
static const Double_t LowLimit[5] = { -200, -200, -200,-200, -10 };
static const Double_t  UpLimit[5] = { 200,  200, 200, 200, 10 };

// global variables for TMinuit
static TVector3 gBHDHitPos;
static Double_t gBHDWeight;

static Int_t gBLC1NumOfHits;
static TVector3 gBLC1HitPos[MAX_NUM_OF_HITS_BEAM]; //local position
static Double_t gBLC1Weight[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC1Rotation[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC1GX,gBLC1GY,gBLC1GZ,gBLC1GdX,gBLC1GdY;

static Int_t gBLC2NumOfHits;
static TVector3 gBLC2HitPos[MAX_NUM_OF_HITS_BEAM]; //local position
static Double_t gBLC2Weight[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC2Rotation[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC2GX,gBLC2GY,gBLC2GZ,gBLC2GdX,gBLC2GdY;

static TMatrixD gD5Matrix1st;
static TMatrixD gD5Matrix2nd[5];

static TMatrixD gBLC1VIMatrix1st;

namespace
{
  const std::string& class_name("BeamSpectrometer");
  const TransferMatrixMan& gTM = TransferMatrixMan::GetInstance();
  const BLDCWireMapMan& gBLDC = BLDCWireMapMan::GetInstance();
  const double mm=0.1;
  const double cm=10*mm;
  const double unit=cm;
}

// --------------------------------------------------------------//
// functions for TMinuit
static void BLC1GtoL( Double_t *in, Double_t *out ){
  double tmp1=TMath::Tan(in[1]/1000.);
  double tmp3=TMath::Tan(in[3]/1000.);

  TVector3 pos(in[0]-tmp1*gBLC1GZ-gBLC1GX,in[2]-tmp3*gBLC1GZ-gBLC1GY,0);
  pos.RotateY(-gBLC1GdY);
  pos.RotateX(-gBLC1GdX);
  TVector3 dir(tmp1,tmp3,1);
  dir.RotateY(-gBLC1GdY);
  dir.RotateX(-gBLC1GdX);

  out[1]=dir.X()/dir.Z();              //dx
  out[0]=pos.X()-out[1]*pos.Z(); //x
  out[3]=dir.Y()/dir.Z();              //dy
  out[2]=pos.Y()-out[3]*pos.Z(); //y
  out[4]=0.;
  out[5]=in[4];

  out[1]=TMath::ATan(out[1])*1000.;
  out[3]=TMath::ATan(out[3])*1000.;
  //  std::cout<<"blc1gtol:\t"<<in[4]<<"\t"<<out[5]<<std::endl;
}
static void ParBLC1toBLC2( Double_t *parblc1, Double_t *parblc2){
  // mm, mrad
  double parin[6]={parblc1[0],parblc1[1],parblc1[2],parblc1[3],0.,parblc1[4]}; //global
  parin[0]/=unit;
  parin[2]/=unit;
  TMatrixD in;
  in.Use(6,1,parin);
  //  std::cout<<"blc1toblc:\t"<<parblc1[4]<<"\t"<<in[5][0]<<std::endl;
  TMatrixD out1st;
  out1st.ResizeTo(6,1);
  //  gD5Matrix1st.Print();
  //  in.Print();
  out1st.Mult(gD5Matrix1st,in);

  for(int i=0;i<6;i++){
    double out2nd=0;
    if(i<5)
      for(int j=0;j<6;j++)
	for(int k=0;k<6;k++)
	  out2nd+=gD5Matrix2nd[i][j][k]*in[j][0]*in[k][0];
    parblc2[i]=out1st[i][0]+out2nd;
  }
  parblc2[1]=TMath::Tan(parblc2[1]/1000.); // mrad to dx/dz
  parblc2[3]=TMath::Tan(parblc2[3]/1000.);
  parblc2[0]*=cm; 
  parblc2[2]*=cm;

  // global -> local
  TVector3 pos(parblc2[0]-parblc2[1]*gBLC2GZ-gBLC2GX,parblc2[2]-parblc2[3]*gBLC2GZ-gBLC2GY,0);
  pos.RotateY(-gBLC2GdY);
  pos.RotateX(-gBLC2GdX);
  TVector3 dir(parblc2[1],parblc2[3],1);
  dir.RotateY(-gBLC2GdY);
  dir.RotateX(-gBLC2GdX);

  parblc2[1]=dir.X()/dir.Z();              //dx
  parblc2[0]=pos.X()-parblc2[1]*pos.Z(); //x
  parblc2[3]=dir.Y()/dir.Z();              //dy
  parblc2[2]=pos.Y()-parblc2[3]*pos.Z(); //y

  parblc2[1]=TMath::ATan(parblc2[1])*1000; // dx/dz -> mrad
  parblc2[3]=TMath::ATan(parblc2[3])*1000;
}

static void ParBLC1toVI( Double_t *parblc1, Double_t *parvi){
  // mm, mrad
  double parin[6]={parblc1[0],parblc1[1],parblc1[2],parblc1[3],0.,parblc1[4]}; //global
  parin[0]/=unit;
  parin[2]/=unit;
  TMatrixD in;
  in.Use(6,1,parin);
  TMatrixD out1st;
  out1st.ResizeTo(6,1);
  //  gD5Matrix1st.Print();
  //  in.Print();
  out1st.Mult(gBLC1VIMatrix1st,in);
  for(int i=0;i<6;i++){
    double out2nd=0;
#if 0
    if(i<5)
      for(int j=0;j<6;j++)
	for(int k=0;k<6;k++)
	  out2nd+=gD5Matrix2nd[i][j][k]*in[j][0]*in[k][0];
#endif
    parvi[i]=out1st[i][0]+out2nd;
  }
  parvi[0]*=cm; //cm->mm
  parvi[2]*=cm;
}

static void fcn( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  Double_t chisq=0.;  Int_t dof = 0;
  double parblc1[6];
  BLC1GtoL(par,parblc1);

  for( Int_t i=0; i<gBLC1NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = gBLC1HitPos[i].x();
    fitx = parblc1[0]+TMath::Tan(parblc1[1]/1000.)*gBLC1HitPos[i].z();
    fity = parblc1[2]+TMath::Tan(parblc1[3]/1000.)*gBLC1HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-gBLC1Rotation[i]);
    chisq +=( (hitx-vec2.X())*(hitx-vec2.X()) )/gBLC1Weight[i]/gBLC1Weight[i];
    dof++;
  }

  Double_t parblc2[6];
  ParBLC1toBLC2(par,parblc2);
#if DEBUG2
  std::cout<<"parblc2=";
  for(int i=0;i<6;i++) std::cout<<"  "<<parblc2[i];
  std::cout<<std::endl;
#endif
  for( Int_t i=0; i<gBLC2NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = gBLC2HitPos[i].x();
    fitx = parblc2[0]+TMath::Tan(parblc2[1]/1000.)*gBLC2HitPos[i].z();
    fity = parblc2[2]+TMath::Tan(parblc2[3]/1000.)*gBLC2HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-gBLC2Rotation[i]);
    chisq +=( (hitx-vec2.X())*(hitx-vec2.X()) )/gBLC2Weight[i]/gBLC2Weight[i];
#if DEBUG2
    std::cout<<"BLC2: "<<i<<"  "<<hitx<<"  "<<vec2.X()<<std::endl;
#endif
    dof++;
  }
  f = chisq/(dof-5);
#if DEBUG
  std::cout<<f<<std::endl;
#endif
}

static void fcn2( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  Double_t chisq=0.;  Int_t dof = 0;
  double parblc1[6];
  BLC1GtoL(par,parblc1);

  for( Int_t i=0; i<gBLC1NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = gBLC1HitPos[i].x();
    fitx = parblc1[0]+TMath::Tan(parblc1[1]/1000.)*gBLC1HitPos[i].z();
    fity = parblc1[2]+TMath::Tan(parblc1[3]/1000.)*gBLC1HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-gBLC1Rotation[i]);
    chisq +=( (hitx-vec2.X())*(hitx-vec2.X()) )/gBLC1Weight[i]/gBLC1Weight[i];
    dof++;
  }

  Double_t parbhd[6];
  ParBLC1toVI(par,parbhd);
  double hitx = gBHDHitPos.x();
  double fitx = parbhd[0];
  chisq +=( (hitx-fitx)*(hitx-fitx) )/gBHDWeight/gBHDWeight;
  dof++;

  f = chisq/(dof-5);
#if DEBUG
  std::cout<<f<<std::endl;
#endif
}

BeamSpectrometer::BeamSpectrometer()
{
  Clear();  

  if( !gTM.IsReady() ){
    std::cerr << "Error in BeamSpectrometer!  Set InitializeParameter<TransferMatrixMan>." << std::endl;
    return;
  }
  if( !gBLDC.IsReady() ){
    std::cerr << "Error in BeamSpectrometer!  Set InitializeParameter <BLDCWireMapMan>." << std::endl;
    return;
  }

  minuit = new TMinuit(5);
  TROOT minexam("BeamSpectrometer","Beam fit using TMinuit");  
  MomCenter=gTM.GetCentralMomentum();
  D5Matrix1st.ResizeTo(6,6,-1);
  const double *ele = gTM.GetD5Matrix();
  D5Matrix1st.Use(6,6,ele);

  BLC1VIMatrix1st.ResizeTo(6,6,-1);
  const double *ele2 = gTM.GetBLC1VIMatrix();
  BLC1VIMatrix1st.Use(6,6,ele2);
  BLC1VIMatrix1st.Invert();

  for(int i=0; i<5 ;i++){
    D5Matrix2nd[i].ResizeTo(6,6,-1);
    ele = gTM.GetD5Matrix2nd(i);
    D5Matrix2nd[i].Use(6,6,ele);
  }
}

BeamSpectrometer::~BeamSpectrometer()
{
  delete minuit;
}

void BeamSpectrometer::SetLocalTrack(const LocalTrack &blc1,const LocalTrack &blc2)
{
  BLC1Track=LocalTrack(blc1);
  BLC2Track=LocalTrack(blc2);
  CalcInitPar();
}

void BeamSpectrometer::SetTrackParams()
{
  TVector3 gpos,grot;
  gBLDC.GetGParam(DetIdBLC1a,gpos,grot);
  BLC1GX=gpos.X();
  BLC1GY=gpos.Y();
  BLC1GZ=gpos.Z();
  BLC1GdX=grot.X()*TMath::DegToRad(); // deg -> rad
  BLC1GdY=grot.Y()*TMath::DegToRad();
  int ihit=0;
  for(int xy=0;xy<2;xy++){
    for( int i=0; i<blc1track()->nhit(xy); i++ ) {
      if(ihit>=MAX_NUM_OF_HITS_BEAM) break;
      TrackHit hit=blc1track()->hit(xy,i);
      if(xy==0){	
	BLC1HitPos[ihit]=hit.hitpos;
      }else{
	BLC1HitPos[ihit].SetXYZ(hit.hitpos.Y(),hit.hitpos.X(),hit.hitpos.Z());
      }
      BLC1Weight[ihit] = 0.2*mm;
      //conf->GetReslMapManager()->GetParam(hit->cid(),hit->layer(),0,BLC1Weight[i],tmp);
      //conf->GetReslMapManager()->GetResolution(CID_BLC,Weight)
      BLC1Rotation[ihit]= (grot.Z()+hit.rotation)*TMath::DegToRad();
#if DEBUG
      std::cout<<"layer: "	<<std::setw(5)<<hit.layer
	       <<", wire: "	<<std::setw(5)<<hit.wire
	       <<", hitpos: "	<<std::setw(10)<<hit.hitpos.X()
	       <<", "		<<std::setw(10)<<hit.hitpos.Y()
	       <<", "		<<std::setw(10)<<hit.hitpos.Z()
	       <<", wpos: "	<<std::setw(10)<<hit.wpos.X()
	       <<", "		<<std::setw(10)<<hit.wpos.Y()
	       <<", "		<<std::setw(10)<<hit.wpos.Z()
	       <<", rot: "      <<std::setw(10)<<BLC1Rotation[ihit]
	       <<", grotz: "    <<std::setw(10)<<grot.Z()
	       <<", rot: "	<<std::setw(10)<<hit.rotation
	       <<std::endl;
#endif
      ihit++;
    }
  }
  BLC1NumOfHits=ihit;

  gBLDC.GetGParam(DetIdBLC2a,gpos,grot);
  BLC2GX=gpos.X();
  BLC2GY=gpos.Y();
  BLC2GZ=gpos.Z()+130.*cm;
  BLC2GdX=grot.X()*TMath::DegToRad(); // deg -> rad
  BLC2GdY=grot.Y()*TMath::DegToRad();

  ihit=0;
  for(int xy=0;xy<2;xy++){
    for( int i=0; i<blc2track()->nhit(xy); i++ ) {
      if(ihit>=MAX_NUM_OF_HITS_BEAM) break;
      TrackHit hit=blc2track()->hit(xy,i);
      if(xy==0){	
	BLC2HitPos[ihit]=hit.hitpos;//+TVector3(0,0,1300);
      }else{
	BLC2HitPos[ihit].SetXYZ(hit.hitpos.Y(),hit.hitpos.X(),hit.hitpos.Z());//+130*cm);
      }
      BLC2Weight[ihit] = 0.2*mm;
      //conf->GetReslMapManager()->GetParam(hit->cid(),hit->layer(),0,BLC2Weight[i],tmp);
      //conf->GetReslMapManager()->GetResolution(CID_BLC,Weight)
      BLC2Rotation[ihit]= (grot.Z()+hit.rotation)*TMath::DegToRad();
      ihit++;
    }
  }
  BLC2NumOfHits=ihit;
  return;
}
void BeamSpectrometer::Clear()
{
  BLC1Track.Clear();
  BLC2Track.Clear();
  for( int i=0; i<5; i++ ){
    Par[i] = Err[i] = -999.;
  }

  BLC1NumOfHits = BLC2NumOfHits =-999;
  FitDof = FitStat = -999;
  FitChi2 = -999.;
  //  MomCenter=-999.;

  BLC1GX=BLC1GY=BLC1GZ=BLC1GdX=BLC1GdY=-999.;  
  BLC2GX=BLC2GY=BLC2GZ=BLC2GdX=BLC2GdY=-999.;  
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ){
    BLC1HitPos[i].SetXYZ(-999,-999,-999);
    BLC1Weight[i] = -999.;
    BLC1Rotation[i]=-999.;
    BLC2HitPos[i].SetXYZ(-999,-999,-999);
    BLC2Weight[i] = -999.;
    BLC2Rotation[i]=-999.;
  }
}

void BeamSpectrometer::GetBLC2fromBLC1(LocalTrack *blc1,double *parblc,const double &momentum)
{
}

bool BeamSpectrometer::CalcInitPar()
{
  double a,b,c,d,tmp;
  blc1track()->gabc(tmp,b,a);
  blc1track()->gdef(tmp,d,c);
  double a2,b2,c2,d2;
  blc2track()->gabc(tmp,b2,a2);
  blc2track()->gdef(tmp,d2,c2);

  double x,y;
  blc1track()->XYPosatZ(0,x,y);
  double x2,y2;
  blc2track()->XYPosatZ(-130*cm,x2,y2);

  Par[0]=x;
  Par[1]=TMath::ATan(-b)*1000; //mrad
  Par[2]=y;
  Par[3]=TMath::ATan(-d)*1000; //mrad
  double cang=-1./12.5;
  double cplus=3.75/125.;
  double cminus=-0.11/125.;  
  double angle=-TMath::ATan(-b)*1000.+TMath::ATan(-b2)*1000.;
  double cmom=cang*angle+cplus*(x+x2)/unit+cminus*(x2-x)/unit;
  cmom *= -1;
  //  std::cout<<"cmom:\t"<<cmom<<std::endl;
  Par[4]=cmom;
  return true;
}
void BeamSpectrometer::SetParameters( const double *param )
{
  for( int i=0; i<5; i++ ) Par[i] = (Double_t)param[i];
}

void BeamSpectrometer::GetParameters( double *param )
{
  for( int i=0; i<5; i++ ) param[i] = (double)Par[i];
}

void BeamSpectrometer::SetGlobalVariables()
{
  gBLC1NumOfHits = BLC1NumOfHits;
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ){
    gBLC1HitPos[i] = BLC1HitPos[i];
    gBLC1Weight[i] = BLC1Weight[i];
    gBLC1Rotation[i]= BLC1Rotation[i];
  }
  gBLC2NumOfHits = BLC2NumOfHits;
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ){
    gBLC2HitPos[i] = BLC2HitPos[i];
    gBLC2Weight[i] = BLC2Weight[i];
    gBLC2Rotation[i]= BLC2Rotation[i];
  }

  gBLC1GX=BLC1GX;
  gBLC1GY=BLC1GY;
  gBLC1GZ=BLC1GZ;
  gBLC1GdX=BLC1GdX;
  gBLC1GdY=BLC1GdY;
  
  gBLC2GX=BLC2GX;
  gBLC2GY=BLC2GY;
  gBLC2GZ=BLC2GZ;
  gBLC2GdX=BLC2GdX;
  gBLC2GdY=BLC2GdY;

  gBLC1VIMatrix1st.Use(6,6,BLC1VIMatrix1st.GetMatrixArray());
  gD5Matrix1st.Use(6,6,D5Matrix1st.GetMatrixArray());
  for(int i=0;i<5;i++)
    gD5Matrix2nd[i].Use(6,6,D5Matrix2nd[i].GetMatrixArray());
}

void BeamSpectrometer::TMinuitFit(LocalTrack *blc1,LocalTrack *blc2)
{
  SetLocalTrack(*blc1,*blc2);
  fit();
}

void BeamSpectrometer::TMinuitFitBHD(TVector3 bhdpos,double bhdw)
{
  BHDHitPos=bhdpos;
  BHDWeight=bhdw;
  gBHDHitPos=BHDHitPos;
  gBHDWeight=BHDWeight;
  fit2();
}

void BeamSpectrometer::fit()
{
#if 0
  std::cout << "!!! BeamSpectrometer::fit() !!!" << std::endl;
#endif
  Int_t plevel=-1;
  //  Int_t plevel=1;
  SetTrackParams();
  SetGlobalVariables();
  minuit->SetPrintLevel( plevel );
  minuit->SetFCN( fcn );
  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  //minuit->mnexcm("SET ERR", arglist,1,ierflg);
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  // Set starting values and step sizes for parameters
  TString name[5] = {"x_blc1", "theta_blc1", "y_blc1","phi_blc1","dp"};
  for( Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i],Par[i],FitStep[i],LowLimit[i],UpLimit[i],ierflg);
  }

  minuit->Command("SET STRategy 0");
  // Now ready for minimization step
  arglist[0] = 1000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

#if DEBUG
  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm,  errdef, nvpar, nparx, icstat);
  FitStat = icstat;
  minuit->mnprin(5,amin);
#endif
  
  Int_t err;
  Double_t bnd1, bnd2;
  for( Int_t i=0; i<5; i++ ){
    minuit->mnpout(i, name[i], Par[i], Err[i], bnd1, bnd2, err);
  }  
  CalcChi2();
  //  std::cout<<<<"  "<<FitChi2<<std::endl;
}

void BeamSpectrometer::fit2()
{
#if 0
  std::cout << "!!! BeamSpectrometer::fit() !!!" << std::endl;
#endif
  Int_t plevel=-1;
  minuit->SetPrintLevel( plevel );
  minuit->SetFCN( fcn2 );
  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  //minuit->mnexcm("SET ERR", arglist,1,ierflg);
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  // Set starting values and step sizes for parameters
  TString name[5] = {"x_blc1", "theta_blc1", "y_blc1","phi_blc1","dp"};
  for( Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i],Par[i],FitStep[i],LowLimit[i],UpLimit[i],ierflg);
  }

  minuit->Command("SET STRategy 0");
  // Now ready for minimization step
  arglist[0] = 1000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

#if DEBUG
  //#if 1
  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm,  errdef, nvpar, nparx, icstat);
  FitStat = icstat;
  minuit->mnprin(5,amin);
#endif
  
  Int_t err;
  Double_t bnd1, bnd2;
  for( Int_t i=0; i<5; i++ ){
    minuit->mnpout(i, name[i], Par[i], Err[i], bnd1, bnd2, err);
  }  
  FitChi2=minuit->fAmin;
  //  std::cout<<minuit->fAmin<<"  "<<FitChi2<<std::endl;
}

void BeamSpectrometer::CalcParBLC1toBLC2( double *parblc1, double *parblc2, const bool &GTOL){
  double parin[6]={parblc1[0],parblc1[1],parblc1[2],parblc1[3],parblc1[4],parblc1[5]};
  parin[0]/=unit;
  parin[2]/=unit;
  TMatrixD in;
  in.Use(6,1,parin);
  TMatrixD out1st;
  out1st.ResizeTo(6,1);
  out1st.Mult(D5Matrix1st,in);
  for(int i=0;i<6;i++){
    double out2nd=0;
    if(i<5)
      for(int j=0;j<6;j++)
	for(int k=0;k<6;k++)
	  out2nd+=D5Matrix2nd[i][j][k]*in[j][0]*in[k][0];
    parblc2[i]=out1st[i][0]+out2nd;
  }
  parblc2[0]*=cm;
  parblc2[2]*=cm;
  parblc2[4]*=cm;
  if(GTOL){
    parblc2[1]=TMath::Tan(parblc2[1]/1000.); // mrad -> dx/dz
    parblc2[3]=TMath::Tan(parblc2[3]/1000.);

    TVector3 pos(parblc2[0]-parblc2[1]*BLC2GZ-BLC2GX,parblc2[2]-parblc2[3]*BLC2GZ-BLC2GY,0);
    pos.RotateY(-BLC2GdY);
    pos.RotateX(-BLC2GdX);
    TVector3 dir(parblc2[1],parblc2[3],1);
    dir.RotateY(-BLC2GdY);
    dir.RotateX(-BLC2GdX);
    
    parblc2[1]=dir.X()/dir.Z();              //dx
    parblc2[0]=pos.X()-parblc2[1]*pos.Z(); //x
    parblc2[3]=dir.Y()/dir.Z();              //dy
    parblc2[2]=pos.Y()-parblc2[3]*pos.Z(); //y

    parblc2[1]=TMath::ATan(parblc2[1])*1000; //dx/dz ->mmrad
    parblc2[3]=TMath::ATan(parblc2[3])*1000;
  }
}
void BeamSpectrometer::CalcBLC1GtoL( Double_t *in, Double_t *out ){

  double tmp[6];
  tmp[0]=in[0];
  tmp[1]=TMath::Tan(in[1]/1000.);
  tmp[2]=in[2];
  tmp[3]=TMath::Tan(in[3]/1000.);
  TVector3 pos(in[0]-tmp[1]*BLC1GZ-BLC1GX,tmp[2]-tmp[3]*BLC1GZ-BLC1GY,0);
  pos.RotateY(-BLC1GdY);
  pos.RotateX(-BLC1GdX);
  TVector3 dir(tmp[1],tmp[3],1);
  dir.RotateY(-BLC1GdY);
  dir.RotateX(-BLC1GdX);

  out[1]=dir.X()/dir.Z();              //dx
  out[0]=pos.X()-out[1]*pos.Z(); //x
  out[3]=dir.Y()/dir.Z();              //dy
  out[2]=pos.Y()-out[3]*pos.Z(); //y

  out[4]=in[4];
  out[5]=in[5];
  //  std::cout<<"blc1gtol:\t"<<in[4]<<"\t"<<out[5]<<std::endl;

  out[1]=TMath::ATan(out[1])*1000.;
  out[3]=TMath::ATan(out[3])*1000.;
  //  for(int i=0;i<6;i++) std::cout<<in[i]<<"  "; std::cout<<std::endl;
  //  for(int i=0;i<6;i++) std::cout<<out[i]<<"  "; std::cout<<std::endl;
}

void BeamSpectrometer::CalcChi2()
{  
  double chisq=0.;
  int dof = 0;  
  double parorg[6]={Par[0],Par[1],Par[2],Par[3],0.,Par[4]};
  double parblc2[6];
  double parblc1[6];
  CalcBLC1GtoL(parorg,parblc1);
  CalcParBLC1toBLC2(parorg,parblc2,1);
#if DEBUG
  std::cout<<"Par0,1"<<parorg[0]<<" , "<<parorg[1]<<"  , "<<parorg[2]<<std::endl;
  std::cout<<"BLC1X: "<<parblc1[0]<<std::endl;
  std::cout<<"BLC1dX:  "<<parblc1[1]<<" , "<<TMath::Tan(parblc1[1]/1000.)<<std::endl;
#endif
  for( Int_t i=0; i<BLC1NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = BLC1HitPos[i].x();
    fitx = parblc1[0]+TMath::Tan(parblc1[1]/1000.)*BLC1HitPos[i].z();
    fity = parblc1[2]+TMath::Tan(parblc1[3]/1000.)*BLC1HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-BLC1Rotation[i]);    
    double tmp=( (hitx-vec2.X())*(hitx-vec2.X()) )/BLC1Weight[i]/BLC1Weight[i];
    chisq +=tmp;
#if DEBUG
    std::cout<<"BLC1 "<<i<<"  "<<hitx<<"  "<<vec2.X()<<"  "<<tmp<<"  "<<BLC1Rotation[i]<<std::endl;
#endif
    dof++;
  }
#if DEBUG
  std::cout<<"parblc2=";
  for(int i=0;i<6;i++) std::cout<<"  "<<parblc2[i];
  std::cout<<std::endl;
#endif
  for( Int_t i=0; i<BLC2NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = BLC2HitPos[i].x();
    fitx = parblc2[0]+TMath::Tan(parblc2[1]/1000.)*BLC2HitPos[i].z();
    fity = parblc2[2]+TMath::Tan(parblc2[3]/1000.)*BLC2HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-BLC2Rotation[i]);
    double tmp=( (hitx-vec2.X())*(hitx-vec2.X()) )/BLC2Weight[i]/BLC2Weight[i];
    chisq +=tmp;
#if DEBUG
    std::cout<<"BLC2 "<<i<<"  "<<hitx<<"  "<<BLC2HitPos[i].y()<<"  "<<vec2.X()<<"  "<<tmp<<std::endl;
#endif
    dof++;
  }

  FitDof = dof - 5;
  FitChi2 = chisq / FitDof;
}

