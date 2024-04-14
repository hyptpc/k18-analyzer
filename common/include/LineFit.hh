/* LineFit.h */

// ======================================================== //
// Reference for circle(line) fit :                        //   
//   http://www-jlc.kek.jp/subg/offl/lib/docs/line_manip/  //
// ======================================================== //

#ifndef LineFit_h
#define LineFit_h 1

#include "TROOT.h"
#include "TVector3.h"
#include "TMinuit.h"
//#include "MathTools.h"
const int MAX_NUM_OF_HITS_L = 60;

class LineFit
{
  // ------------------------------------------- //
  // Parameters for circle fit :                 //
  //   par[4] = { a,b,c,d }  //
  //   ax + by + cz + d = 0
  // ------------------------------------------- //
 private:
  Double_t Par[6];
  Double_t Err[6];
  
  Int_t NumOfHits;

  TVector3 WirePos[MAX_NUM_OF_HITS_L];
  TVector3 WireDir[MAX_NUM_OF_HITS_L];
  Double_t Drift[MAX_NUM_OF_HITS_L];
  Double_t Weight[MAX_NUM_OF_HITS_L];
  
  Double_t FitChi2;
  Int_t FitDof;
  Int_t FitStat;

 private:
  TMinuit *minuit;

 public:
  LineFit();
  LineFit(const double *initPar, const TVector3 *wirepos,
	   const TVector3 *wiredir,const double *weight,
	   const double *drift, const int &numofhit);
  ~LineFit();

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
