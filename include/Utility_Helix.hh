#ifndef Utillity_Helix_h
#define Utillity_Helix_h 1

#include "TVector3.h"
#include "TLorentzVector.h"

namespace Utility_Helix
{
  //Temporary
  const double FieldPara = -1.;// -1T
  const double TargetPositionZ = -143.0; //Target z position
  const double TargetRadius = 27.0; //Target Radius 
  const double TargetHeight = 50.0; //Target Height +-50 mm 
  const double TargetXC = 30.0/2.0; //Diamond target length in x
  const double TargetYC = 20.0/2.0; //Diamond target length in y
  const double TargetZC = 20.0/2.0; //Diamond target length in z
  const double KuramaXshift = 87.0; //Kurama shift in x-direction 

  //const double TargetDensity = 0.0709; //liquid H2, g/cm3  
  const double TargetDensity = 3.51; //diamond, g/cm3    
  //const double TargetDensity = 0.936353; //diamond, g/cm3           
  //Residual
  //  TVector3 GaussPosition( TVector3 tv );
  
  
  //Helix Fit
  double CalcRad( const double x, const double y );
  double CalcHelixPhi(const double &x, const double &y, const double *par);
  bool HelixToHelix(const double *par1, const double  *par2,TVector3 &fitpos1,TVector3 &fitpos2, double &dis );
  TVector3 GetPosition(const double &helixphi, const double *par);
  bool PointToHelix(const TVector3 &hitpos, TVector3 &fitpos, const double *par);
  double dfunc_PTH( const TVector3 &pos, const double &helixphi, const double *par);
  bool PointToHelix(const TVector3 &hitpos, const double *par, TVector3 &fitpos,double &dis);
  TVector3 CalcHelixMom(const double par[5], double z);

  //Energy Loss
  bool EnergyLossBethe( const double distance, const TLorentzVector lv, TLorentzVector &lv_corr);
  double calc_dE_dx(double beta);

  //Energy Loss
  double MyGamma(double beta);
  double MyBeta(double energy,double mormentum);
  double CalcHelixArc(const double par[5], const TVector3 &pos1,const TVector3 &pos2);
  bool CheckInTarget(TVector3 pos);
  double CalcFlightLengthInTarget( const double par[5], const TVector3 vertex );

  bool CheckInTargetC(TVector3 pos);
  double CalcFlightLengthInTargetC( const double par[5], const TVector3 vertex );

  double CalcHelixPhi(const TVector3 &pos, const double *par);

  TVector3 GetVertex(TVector3 v1, TVector3 v2);
  void LineToLine(const TVector3 &x1, const TVector3 &p1, const TVector3 &x2, const TVector3 &p2, TVector3 &vp, double &dist);

  bool LineToHelix(const TVector3 &a, const TVector3 &dline, const double *par, TVector3 &lnest, TVector3 &hnest, double &dis);

  double dfunc_LTH(const TVector3 &lpos,const TVector3 &dline,const double &helixphi,const double *par);

  
  //Beam modification
  TVector3 BeamPMod( TVector3 tv );

  //Tof energy modification
  double TofEMod( double edep );

}



#endif

