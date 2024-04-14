/* HelixFit.h */

// ======================================================== //
// Reference for circle(helix) fit :                        //   
//   http://www-jlc.kek.jp/subg/offl/lib/docs/helix_manip/  //
// ======================================================== //

#ifndef HelixFit_h
#define HelixFit_h 1

#include "TROOT.h"
#include "TVector3.h"
#include "TMinuit.h"
const int MAX_NUM_OF_HITS_H2 = 60;

class HelixFit
{
  // ------------------------------------------- //
  // Parameters for circle fit :                 //
  //   par[3] = { d_rho, phi_0, alpha/k(=rho), dZ, tan_theta }  //
  // ------------------------------------------- //
 private:
  Double_t Par[5];
  Double_t Err[5];
  
  Int_t NumOfHits;

  TVector3 WirePos[MAX_NUM_OF_HITS_H2];
  TVector3 WireDir[MAX_NUM_OF_HITS_H2];
  Double_t Drift[MAX_NUM_OF_HITS_H2];
  Double_t Weight[MAX_NUM_OF_HITS_H2];
  
  Double_t FitChi2;
  Int_t FitDof;
  Int_t FitStat;

 private:
  TMinuit *minuit;

 public:
  HelixFit();
  HelixFit(const double *initPar, const TVector3 *wirepos,
	   const TVector3 *wiredir,const double *weight,
	   const double *drift, const int &numofhit);
  ~HelixFit();

 public:
  void SetParameters( const double *param );
  void SetParameter( const int &i, const double &param){ Par[i]=(Double_t)param; }
  void SetWirePos( const TVector3 hitpos, const int &numofhit );
  void SetWireDir( const TVector3 hitpos, const int &numofhit );
  void SetWeight( const double *weight, const int &numofhit );
  void SetDrift( const double *weight, const int &numofhit );
  void SetNumOfHit( const int &i ){ NumOfHits=(Int_t)i; }

  void GetParameters( double *param );
  double param( const int &i ){ return (double)Par[i]; }
  double chisquare() { return (double)FitChi2; }
  int dof() { return (int)FitDof; }
  int stat() { return (int)FitStat; }

 public:
  void SetGlobalVariables();
  void fit();
  void CalcChi2();
};

#endif 
