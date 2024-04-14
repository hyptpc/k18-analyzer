// Helixfit2.cpp
#include <iostream>
#include "HelixFit.hh"
#include "MathTools.hh"
#include "TMath.h"

#define DEBUG 0

namespace{
  // parameters for TMinuit  
  static const Double_t  FitStep[5] = { 1e-20, 1e-21, 1e-22,1e-20,1e-20 };
  static const Double_t LowLimit[5] = { -500, -2*TMath::Pi(), -0.5,-1000,-10 };
  static const Double_t  UpLimit[5] = { 500, 2*TMath::Pi(), 0.5,1000,10 };
  
  // global variables for TMinuit
  static Int_t gNumOfHits;
  static TVector3 gWirePos[MAX_NUM_OF_HITS_H2];
  static TVector3 gWireDir[MAX_NUM_OF_HITS_H2];
  static Double_t gWeight[MAX_NUM_OF_HITS_H2];
  static Double_t gDrift[MAX_NUM_OF_HITS_H2];

  // --------------------------------------------------------------//
  // functions for TMinuit
  static bool LineToHelix(const TVector3 &a, const  TVector3 &dline ,
			  const double *par, double &dis, TVector3 &lnest)
  {
    // a: origin of the line
    // dline: direction of the line
    // par: helix parameter
    // dis: closest distance from the helix to the line
    //
    TVector3 dlineu=dline.Unit();
    double phi=0,phi_b=0,phi_a=0;
    phi=math::CalcHelixPhi(a.x(),a.y(),par);
  
    double philen=1./par[2]*sqrt(1+par[4]*par[4]);
    double dist=1.;
    int trial=14; //sigma_position<1.2 micron 

    // find LTH initial param
    while(dist<1024)
      {
	phi_b=phi-dist/philen;
	phi_a=phi+dist/philen;
	double dlen_b=math::dfunc_LTH(a,dlineu,phi_b,par);
	double dlen_a=math::dfunc_LTH(a,dlineu,phi_a,par);
	if(dlen_b*dlen_a<=0) break;
	else {dist*=2;trial++;}
      }
  
    if(dist>=1024) 
      {
#if DEBUG
	std::cout<<"Can not find LTH inital param!!"<<std::endl;
#endif
	return false;
      }
 
    // Bisection Method
    for(int i=0;i<trial;i++)
      {
	phi_b=phi-dist/philen;
	phi_a=phi+dist/philen;
	double dlen_b=math::dfunc_LTH(a,dlineu,phi_b,par);
	//      double dlen_a=dfunc_LTH(a,dlineu,phi_a,par);
	double dlen=math::dfunc_LTH(a,dlineu,phi,par);
	if(dlen*dlen_b<=0 ) {phi=(phi_b+phi)/2.0;dist=dist/2.0;}
	else  {phi=(phi_a+phi)/2.0;dist=dist/2.0;}      
      }
    //
    // helix nearest point
    TVector3 hnest=math::GetPosition(phi,par);
    double k=(hnest-a)*dlineu;
    // line nearest point
    lnest=a+k*dlineu;
    TVector3 tmp;
    tmp=hnest-lnest;
    dis=tmp.Mag();
    return true;
  }

  static void fcn2( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
  {
    Double_t chisq=0.;
    Int_t dof = 0;
    Double_t hitx, hity, fitx, fity, dist=0.;
    TVector3 fittmp;
  
    for( Int_t i=0; i<gNumOfHits; i++ ){
      if(gWirePos[i].z()<-1000)
	{
	  hitx = gWirePos[i].x();
	  hity = gWirePos[i].y();
	  //	NearestPos(par, hitx, hity, fitx, fity);
	  math::PointToCircle(par, hitx, hity, fitx, fity,dist);
	  chisq +=pow( ( dist - gDrift[i] )/gWeight[i] , 2 );
	}
      else
	{
	  if(!LineToHelix(gWirePos[i], gWireDir[i],par, dist, fittmp) ) 
	    //	if(!math::PointToHelix(gHitPos[i],fittmp,par) ) 
	    {	
	      //    std::cout<<"Miss Newton Method in fcn!!"<<std::endl;
	      chisq+=9999;
	      continue;
	    }
	  chisq +=pow( ( dist - gDrift[i] )/gWeight[i] , 2 );
	}
      dof++;
    }  
    f = chisq/(dof-5);
  }
} // namespace
// --------------------------------------------------------------//
HelixFit::HelixFit()
{
  for( int i=0; i<5; i++ ){
    Par[i] = Err[i] = -9999.;
  }

  NumOfHits = FitDof = FitStat = -9999;
  FitChi2 = -9999.;

  for( int i=0; i<MAX_NUM_OF_HITS_H2; i++ ){
    WirePos[i].SetXYZ(-9999,-9999,-9999);
    WireDir[i].SetXYZ(-9999,-9999,-9999);
    Weight[i] = -9999.;
    Drift[i] = -9999.;
  }
  
  minuit = new TMinuit(5);
  TROOT minexam("HelixFit","helix fit using TMinuit");
}

HelixFit::HelixFit( const double *initPar, const TVector3 *wirepos,
		    const TVector3 *wiredir,const double *weight,
		    const double *drift, const int &numofhit)
{
  NumOfHits = (Int_t)numofhit;

  for( int i=0; i<5; i++ ) {
    Par[i] = (Double_t)initPar[i];
    Err[i] = -9999.;
  }

  for( int i=0; i<MAX_NUM_OF_HITS_H2; i++ ) {
    if( i<NumOfHits ){
      WirePos[i]= wirepos[i];
      WireDir[i]= wiredir[i];
      Weight[i] = (Double_t)weight[i];      
      Drift[i]  = (Double_t)drift[i];      
    }
    else {
      WirePos[i].SetXYZ(-9999,-9999,-9999);
      WireDir[i].SetXYZ(-9999,-9999,-9999);
      Weight[i] = -9999.;
      Drift[i]  = -9999.;
    }
  }

  FitDof = FitStat = -9999;
  FitChi2 = -9999.;

  minuit = new TMinuit(5);
  TROOT minexam("HelixFit","Helix fit using TMinuit");

  fit();
}

HelixFit::~HelixFit()
{
  delete minuit;
}
void HelixFit::SetParameters( const double *param )
{
  for( int i=0; i<5; i++ ) Par[i] = (Double_t)param[i];
}
void HelixFit::SetWirePos( const TVector3 hitpos,  const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    WirePos[i]=hitpos;
}
void HelixFit::SetWireDir( const TVector3 hitpos,  const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    WireDir[i]=hitpos;
}
void HelixFit::SetWeight( const double *weight, const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    Weight[i] = (Double_t)weight[i];
}
void HelixFit::SetDrift( const double *weight, const int &numofhit )
{
  for( int i=0; i<numofhit; i++ )    Drift[i] = (Double_t)weight[i];
}
void HelixFit::GetParameters( double *param )
{
  for( int i=0; i<5; i++ ) param[i] = (double)Par[i];
}

void HelixFit::SetGlobalVariables()
{
  gNumOfHits = NumOfHits;
  for( int i=0; i<MAX_NUM_OF_HITS_H2; i++ ){
    gWirePos[i] = WirePos[i];
    gWireDir[i] = WireDir[i];
    gWeight[i] = Weight[i];
    gDrift[i] = Drift[i];
  }
}

void HelixFit::fit()
{
#if 0
  std::cout << "!!! HelixFit::fit() !!!" << std::endl;
#endif

#if DEBUG
  Int_t plevel=1;
#else
  Int_t plevel=-1;
#endif

  SetGlobalVariables();

  minuit->SetPrintLevel( plevel );
  minuit->SetFCN( fcn2 );

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  //  minuit->mnexcm("SET ERR", arglist,1,ierflg);
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  // Set starting values and step sizes for parameters
  TString name[5] = {"d_rho", "phi_0", "rho","d_z","tanL"};
  for( Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i],Par[i],FitStep[i],LowLimit[i],UpLimit[i],ierflg);
  }
  minuit->Command("SET STRategy 0");
  // Now ready for minimization step
  arglist[0] = 1000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  //minuit->Command("TMProve 100");

  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm,  errdef, nvpar, nparx, icstat);
#if DEBUG
  minuit->mnprin(5,amin);
#endif
  
  Int_t err;
  Double_t bnd1, bnd2;
  for( Int_t i=0; i<5; i++ ){
    minuit->mnpout(i, name[i], Par[i], Err[i], bnd1, bnd2, err);
  }
#if DEBUG
  for( Int_t i=0; i<5; i++ )
    std::cout<<Par[i]<<"  ";
  std::cout<<std::endl;
#endif
  FitStat = icstat;
  CalcChi2();
}

void HelixFit::CalcChi2()
{  
  Double_t chisq=0.;
  Int_t dof = 0;
  TVector3 fittmp,fittmp2;
  Double_t hitx, hity, fitx, fity, dist;
  for( Int_t i=0; i<NumOfHits; i++ ){    
    if( WirePos[i].z()<-100. )
      {
	hitx = WirePos[i].x();
	hity = WirePos[i].y();
	math::PointToCircle(Par, hitx, hity, fitx, fity,dist);
	chisq +=pow((dist-Drift[i])/Weight[i], 2);
      }
    else if( fabs( WirePos[i].z() )<50. )
      {
	math::LineToHelix( WirePos[i], WireDir[i], Par, fittmp, fittmp2, dist);
	chisq += pow( (dist-Drift[i])/Weight[i] , 2 );
	//	std::cout<<i<<"  dist  "<<dist<<std::endl;
      }
    dof++;
  }
  FitDof = dof - 5;
  FitChi2 = chisq;
}

