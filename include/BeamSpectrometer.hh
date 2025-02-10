/* BeamSpectrometer.h */

// ======================================================== //
// Reference for circle(helix) fit :                        //   
//   http://www-jlc.kek.jp/subg/offl/lib/docs/helix_manip/  //
// ======================================================== //

#ifndef BeamSpectrometer_h
#define BeamSpectrometer_h 1

#include "TROOT.h"
#include "TVector3.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TransferMatrixMan.hh"
#include "LocalTrack.hh"

const int MAX_NUM_OF_HITS_BEAM = 32;

class BeamSpectrometer
{
  // ---------------------------------------------------- //
  // Parameters for beam mom fit :                        //
  //   par[5] = { x_blc1 theta_blc1, y_blc1, phi_blc1, mom }  //
  // ---------------------------------------------------- //
 private:
  Double_t Par[5];
  Double_t Err[5];

  LocalTrack BLC1Track;
  LocalTrack BLC2Track;
  
  Int_t BLC1NumOfHits;
  Int_t BLC2NumOfHits;

  TVector3 BLC1HitPos[MAX_NUM_OF_HITS_BEAM];
  Double_t BLC1Weight[MAX_NUM_OF_HITS_BEAM];
  Double_t BLC1Rotation[MAX_NUM_OF_HITS_BEAM];
  Double_t BLC1GX,BLC1GY,BLC1GZ,BLC1GdX,BLC1GdY;

  TVector3 BLC2HitPos[MAX_NUM_OF_HITS_BEAM];
  Double_t BLC2Weight[MAX_NUM_OF_HITS_BEAM];
  Double_t BLC2Rotation[MAX_NUM_OF_HITS_BEAM];
  Double_t BLC2GX,BLC2GY,BLC2GZ,BLC2GdX,BLC2GdY;

  TVector3 BHDHitPos;
  Double_t BHDWeight;
 
  Double_t FitChi2;
  Int_t FitDof;
  Int_t FitStat;
  TMatrixD BLC1VIMatrix1st;
  TMatrixD D5Matrix1st;
  TMatrixD D5Matrix2nd[5];
  Double_t MomCenter;

 private:
  TMinuit *minuit;

 public:
  BeamSpectrometer();
  ~BeamSpectrometer();
  
 public:
  void SetParameters( const double *param );
  void SetParameter( const int &i, const double &param){ Par[i]=(Double_t)param; }
  void GetParameters( double *param );
  double param( const int &i ){ return (double)Par[i]; }
  double mom(){ return (1.+Par[4]/100.)*MomCenter; }
  double momcenter(){ return MomCenter; }

  double chisquare() { return (double)FitChi2; }
  void GetBLC2fromBLC1 (LocalTrack *blc1,double *blcpar,const double &momentum);
  int dof() { return (int)FitDof; }
  int stat() { return (int)FitStat; }
  LocalTrack *blc1track(){ return &BLC1Track;}
  LocalTrack *blc2track(){ return &BLC2Track;}
  void SetMomentum(const double &momentum) { Par[4]=momentum; }
  void SetLocalTrack(const LocalTrack &blc1,const LocalTrack &blc2);
  void SetTrackParams();

  void TMinuitFitBHD(TVector3 bhdpos,double bhdw);
  void TMinuitFit(LocalTrack *blc1,LocalTrack *blc2);
  void Clear();
  void fit();
  void fit2();
  void CalcParBLC1toBLC2(double *blc1param,double *blc2param,const bool &GTOL=false);
  void CalcBLC1GtoL(double *in,double *out);

 private:
  void SetGlobalVariables();
  void CalcChi2();
  bool CalcInitPar();
};

#endif 
