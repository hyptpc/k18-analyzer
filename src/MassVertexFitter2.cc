//#include "KinFit.cc"
#include "MassVertexFitter2.hh"
#include "TString.h"
#ifndef MassVertexFitter2_cc
#define MassVertexFitter2_cc
#define Debug 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit


MassVertexFitter2::MassVertexFitter2(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_
                                 ,TVector3 V_P, TVector3 V_Q){ 
  //Kinematic fitting for R -> P+Q decay, with vertex constraint
  P=P_;
  Q=Q_;
  R=R_;
  VP=V_P;
  VQ=V_Q;
  VR=0.5*(V_P+V_Q);
  ScaleParams = 0;//Set to false.
  Initialize();
};
double
MassVertexFitter2::CalcHelixDistance2(std::vector<double> pars){
  double p1 = pars[0], th1 = pars[1], ph1 = pars[2], cx1 = pars[3], cy1 = pars[4], z01 = pars[5];
  double p2 = pars[6], th2 = pars[7], ph2 = pars[8], cx2 = pars[9], cy2 = pars[10], z02 = pars[11];
  int sign = charge_config==0? 1 : -1;//config 0 is the configuration where the particle 1 trajectory is on the positive side of the particle 2 trajectory.

  double px1 = p1*sin(th1)*cos(ph1);
  double py1 = p1*sin(th1)*sin(ph1);
  double pz1 = p1*cos(th1);
  double pt1 = p1*sin(th1);//note: th < pi -> sin(th)>0 gaurenteed for polar angle
  double r1 = pt1/RadToMom;
  double dz1 = 1./tan(th1)*sign;

  double px2 = p2*sin(th2)*cos(ph2 );
  double py2 = p2*sin(th2)*sin(ph2 );
  double pz2 = p2*cos(th2);
  double pt2 = p2*sin(th2);
  double r2 = pt2/RadToMom;
  double dz2 = -1./tan(th2)*sign;
  
  double d_cx = cx1 - cx2;
  double d_cy = cy1 - cy2;
  double d_cz = z01 - z02;

  double ph1_h = ph1 - sign*M_PI/2;//For positive particle, ph_helix(= helix_t) = ph_mom - pi/2.  
  double ph2_h = ph2 + sign*M_PI/2;

  double Dx = d_cx + pt1/RadToMom*cos(ph1_h) - pt2/RadToMom*cos(ph2_h);
  double Dy = d_cy + pt1/RadToMom*sin(ph1_h) - pt2/RadToMom*sin(ph2_h);
  double Dz = d_cz + pt1/RadToMom*dz1*(ph1_h) - pt2/RadToMom*dz2*(ph2_h);

  return Dx * Dx + Dy * Dy + Dz * Dz;
}
double
MassVertexFitter2::DerivativeHelixDistance2(int index, std::vector<double> pars){
  double p1 = pars[0], th1 = pars[1], ph1 = pars[2], cx1 = pars[3], cy1 = pars[4], z01 = pars[5];
  double p2 = pars[6], th2 = pars[7], ph2 = pars[8], cx2 = pars[9], cy2 = pars[10], z02 = pars[11];
  int sign = charge_config==0? 1 : -1;//config 0 is the configuration where the particle 1 trajectory is on the positive side of the particle 2 trajectory.
  
  double px1 = p1*sin(th1)*cos(ph1);
  double py1 = p1*sin(th1)*sin(ph1);
  double pz1 = p1*cos(th1);
  double pt1 = p1*sin(th1);
  double r1 = pt1/RadToMom;
  double dz1 = 1./tan(th1) *sign;// pz / pt = cos /sin = cot.

  double px2 = p2*sin(th2)*cos(ph2 );
  double py2 = p2*sin(th2)*sin(ph2 );
  double pz2 = p2*cos(th2);
  double pt2 = p2*sin(th2);
  double r2 = pt2/RadToMom;
  double dz2 = -1./tan(th2) *sign;// charge is opposite, so the sign is flipped.
  
  double d_cx = cx1 - cx2;
  double d_cy = cy1 - cy2;
  double d_cz = z01 - z02;
  
  double ph1_h = ph1 - sign*M_PI/2;//For positive particle, ph_helix(= helix_t) = ph_mom - pi/2.  
  double ph2_h = ph2 + sign*M_PI/2;

  // Dx = d_cx + r1*cos(ph1_h) - r2*cos(ph2_h);
  // r = pt/RadToMom = p*sin(th)/RadToMom, r*dz =p*sin(th)/RadToMom*cot(th) = p*cos(th)/RadToMom.
  // z = z0 + r*dz*helix_t = z0 + r*cos(th) * helix_t* charge
  double Dx = d_cx + p1*sin(th1)/RadToMom*cos(ph1_h) - p2*sin(th2)/RadToMom*cos(ph2_h);
  double Dy = d_cy + p1*sin(th1)/RadToMom*sin(ph1_h) - p2*sin(th2)/RadToMom*sin(ph2_h);
  double Dz = d_cz + p1*cos(th1)/RadToMom*(ph1_h)*sign + p2*cos(th2)/RadToMom*(ph2_h)*sign;// + p2(...) due to charge
  //spatial components are written in terms of kinematic components, to emphasize the dependencies.

  std::vector<std::vector<double>> dD_dpars;
  
  double dDx_dp1 = sin(th1)/RadToMom*cos(ph1_h);
  double dDy_dp1 = sin(th1)/RadToMom*sin(ph1_h);
  double dDz_dp1 = cos(th1)/RadToMom*sign*(ph1_h);
  std::vector<double> dD_dp1 = {dDx_dp1, dDy_dp1, dDz_dp1};
  dD_dpars.push_back(dD_dp1);
  
  double dDx_dth1 = p1*cos(th1)/RadToMom*cos(ph1_h);
  double dDy_dth1 = p1*cos(th1)/RadToMom*sin(ph1_h);
  double dDz_dth1 = -p1*sin(th1)/RadToMom*sign*(ph1_h);
  std::vector<double> dD_dth1 = {dDx_dth1, dDy_dth1, dDz_dth1};
  dD_dpars.push_back(dD_dth1);

  // ph_h = ph -charge*pi/2, -> d ph = d ph_h
  double dDx_dph1 = - p1*sin(th1)/RadToMom*sin(ph1_h);
  double dDy_dph1 = + p1*sin(th1)/RadToMom*cos(ph1_h);
  double dDz_dph1 = p1*cos(th1)/RadToMom*sign;
  std::vector<double> dD_dph1 = {dDx_dph1, dDy_dph1, dDz_dph1};
  dD_dpars.push_back(dD_dph1);

  double dDx_dcx1 = 1;
  double dDy_dcx1 = 0;
  double dDz_dcx1 = 0;
  std::vector<double> dD_dcx1 = {dDx_dcx1, dDy_dcx1, dDz_dcx1};
  dD_dpars.push_back(dD_dcx1);

  double dDx_dcy1 = 0;
  double dDy_dcy1 = 1;
  double dDz_dcy1 = 0;
  std::vector<double> dD_dcy1 = {dDx_dcy1, dDy_dcy1, dDz_dcy1};
  dD_dpars.push_back(dD_dcy1);

  double dDx_dz01 = 0;
  double dDy_dz01 = 0;
  double dDz_dz01 = 1;
  std::vector<double> dD_dz01 = {dDx_dz01, dDy_dz01, dDz_dz01};
  dD_dpars.push_back(dD_dz01);

  //Dx = d_cx + p1*sin(th1)/RadToMom*cos(ph1_h) - p2*sin(th2)/RadToMom*cos(ph2_h);
  //Dy = d_cy + p1*sin(th1)/RadToMom*sin(ph1_h) - p2*sin(th2)/RadToMom*sin(ph2_h);
  //Dz = d_cz + p1*cos(th1)/RadToMom*(ph1_h)*sign + p2*cos(th2)/RadToMom*(ph2_h)*sign;
  double dDx_dp2 = -sin(th2)/RadToMom*cos(ph2_h);
  double dDy_dp2 = -sin(th2)/RadToMom*sin(ph2_h);
  double dDz_dp2 = cos(th2)/RadToMom*sign*(ph2_h);
  std::vector<double> dD_dp2 = {dDx_dp2, dDy_dp2, dDz_dp2};
  dD_dpars.push_back(dD_dp2);

  double dDx_dth2 = -p2*cos(th2)/RadToMom*cos(ph2_h);
  double dDy_dth2 = -p2*cos(th2)/RadToMom*sin(ph2_h);
  double dDz_dth2 = -p2*sin(th2)/RadToMom*sign*(ph2_h);
  std::vector<double> dD_dth2 = {dDx_dth2, dDy_dth2, dDz_dth2};
  dD_dpars.push_back(dD_dth2);

  double dDx_dph2 =  + p2*sin(th2)/RadToMom*sin(ph2_h);
  double dDy_dph2 =  - p2*sin(th2)/RadToMom*cos(ph2_h);
  double dDz_dph2 =  p2*cos(th2)/RadToMom*sign;
  std::vector<double> dD_dph2 = {dDx_dph2, dDy_dph2, dDz_dph2};
  dD_dpars.push_back(dD_dph2);

  double dDx_dcx2 = -1;
  double dDy_dcx2 = 0;
  double dDz_dcx2 = 0;
  std::vector<double> dD_dcx2 = {dDx_dcx2, dDy_dcx2, dDz_dcx2};
  dD_dpars.push_back(dD_dcx2);
  
  double dDx_dcy2 = 0;
  double dDy_dcy2 = -1;
  double dDz_dcy2 = 0;
  std::vector<double> dD_dcy2 = {dDx_dcy2, dDy_dcy2, dDz_dcy2};
  
  dD_dpars.push_back(dD_dcy2);
  double dDx_dz02 = 0;
  double dDy_dz02 = 0;
  double dDz_dz02 = -1;
  std::vector<double> dD_dz02 = {dDx_dz02, dDy_dz02, dDz_dz02};
  dD_dpars.push_back(dD_dz02);

  auto dD_dparam = dD_dpars[index];
  double derivative = 2 * (Dx * dD_dparam[0] + Dy * dD_dparam[1] + Dz * dD_dparam[2]);
  return derivative;
}
void
MassVertexFitter2::GetHelixParameters(double p, double th, double ph, TVector3 position, int charge, double* pars){
  /*
  ph_h: Helix_t in Helix fit. 
  Helix = (cx,cy,z0) + r *(cos(ph_h),sin(ph_h), dz * ph_h)
  For positive charge, momentum is
   p = (pT * cos(ph_h + pi/2), pT * sin(ph_h + pi/2), pz)
  For negative charge, momentum is
   p =  -(pT * cos(ph_h + pi/2), pT * sin(ph_h + pi/2), pz)
   = (pT * cos(ph_h - pi/2), pT * sin(ph_h - pi/2), -pz)
   HelixPars = {cx, cy, z0, r, dz}
  */
  TVector3 mom(p*sin(th)*cos(ph),p*sin(th)*sin(ph),p*cos(th));
  double th_h = charge > 0 ? th : M_PI - th;
  double ph_h = charge > 0 ? ph - M_PI/2 : ph + M_PI/2;

  double pt = hypot(mom.x(),mom.y());
  double r = pt/RadToMom;
  double dz = charge > 0 ? mom.z() / pt : -mom.z() / pt;

  double cx = position.x() - r*cos(ph_h);
  double cy = position.y() - r*sin(ph_h);
  double z0 = position.z() - r*dz*ph_h;
  pars[0] = cx;
  pars[1] = cy;
  pars[2] = z0;
  pars[3] = r;
  pars[4] = dz;
}
TVector3
MassVertexFitter2::CalcHelixPosition(double p, double th, double ph, TVector3 center, int sign){
  double pt = p*sin(th);//Ordinary, sin(th) >> 0, since it is defined as polar angle.
  double r = pt/RadToMom;
  double dz = sign*1./tan(th);

  //center = (cx,cy,z0)
  double ph_h = sign > 0 ? ph - M_PI/2 : ph + M_PI/2;//For positive particle, ph_helix(= helix_t) = ph_mom - pi/2. For negative particle, ph_helix = ph_mom + pi/2.
  double x = center.x() + r*cos(ph_h);
  double y = center.y() + r*sin(ph_h);
  double z = center.z() + r*dz*(ph_h);
  return TVector3(x,y,z);
}
TVector3
MassVertexFitter2::DerivativeHelixPosition(int index, double p, double th, double ph, TVector3 center, int sign){
  double pt = p*sin(th);
  double r = pt/RadToMom;
  double dz = sign*1./tan(th);

  //center = (cx,cy,z0)
  double ph_h = sign > 0 ? ph - M_PI/2 : ph + M_PI/2;//For positive particle, ph_helix(= helix_t) = ph_mom - pi/2. For negative particle, ph_helix = ph_mom + pi/2.
  double x = center.x() + p*sin(th)/RadToMom*cos(ph_h);
  double y = center.y() + p*sin(th)/RadToMom*sin(ph_h);
  double z = center.z() + sign*p*cos(th)/RadToMom*(ph_h);

  std::vector<TVector3> derivatives;
  double dx_dp = sin(th)/RadToMom*cos(ph_h);
  double dy_dp = sin(th)/RadToMom*sin(ph_h);
  double dz_dp = sign*cos(th)/RadToMom*(ph_h); 
  derivatives.push_back(TVector3(dx_dp, dy_dp, dz_dp));

  double dx_dth = p*cos(th)/RadToMom*cos(ph_h);
  double dy_dth = p*cos(th)/RadToMom*sin(ph_h);
  double dz_dth = -sign*p*sin(th)/RadToMom*(ph_h);
  derivatives.push_back(TVector3(dx_dth, dy_dth, dz_dth));

  double dx_dph = -p*sin(th)/RadToMom*sin(ph_h);
  double dy_dph = p*sin(th)/RadToMom*cos(ph_h);
  double dz_dph = sign*p*cos(th)/RadToMom;
  derivatives.push_back(TVector3(dx_dph, dy_dph, dz_dph));

  double dx_dcx = 1;
  double dy_dcx = 0;
  double dz_dcx = 0;
  derivatives.push_back(TVector3(dx_dcx, dy_dcx, dz_dcx));

  double dx_dcy = 0;
  double dy_dcy = 1;
  double dz_dcy = 0;
  derivatives.push_back(TVector3(dx_dcy, dy_dcy, dz_dcy));

  double dx_dz0 = 0;
  double dy_dz0 = 0;
  double dz_dz0 = 1;
  derivatives.push_back(TVector3(dx_dz0, dy_dz0, dz_dz0));

  return derivatives[index];
}

void MassVertexFitter2::Initialize(){
  nMeas = 12;
  // p, th, ph, cx, cy, z0 for both particles.
  nUnkn = 3;// p, th, ph
  nConst = 5;// px, py, pz conservations, E conservation(Mass constraint), Vertex 
  //ndf = nMeas - nUnkn - nConst;
  mP = P.Mag();
  TVector3 TV_P = P.Vect();
  double p_P = TV_P.Mag();  
  double th_P = TV_P.Theta();
  double ph_P = TV_P.Phi();
  double parsP[5];
  GetHelixParameters(p_P, th_P, ph_P, VP, charge_config==0 ? 1 : -1, parsP);
  double cx_P = parsP[0];
  double cy_P = parsP[1];
  double z0_P = parsP[2];


  

  mQ = Q.Mag();
  TVector3 TV_Q = Q.Vect(); 
  double p_Q = TV_Q.Mag();  
  double th_Q = TV_Q.Theta();
  double ph_Q = TV_Q.Phi(); 
  double parsQ[5];
  GetHelixParameters(p_Q, th_Q, ph_Q, VQ, charge_config==0 ? -1 : 1, parsQ);
  double cx_Q = parsQ[0];
  double cy_Q = parsQ[1];
  double z0_Q = parsQ[2];

  TVector3 TV_R = R.Vect(); 
  double p_R = TV_R.Mag();  
  double th_R = TV_R.Theta();
  double ph_R = TV_R.Phi();
  VR = 0.5 * (VP + VQ);//Initial vertex position is set to the middle point of the two trajectory centers. 

  double meas[12];
  double unkn[6];
  double temp[] = {p_P,th_P,ph_P,cx_P,cy_P,z0_P,p_Q,th_Q,ph_Q,cx_Q,cy_Q,z0_Q};
  for(int i=0;i<nMeas;++i)meas[i]=temp[i];
  double temp2[] = {p_R,th_R,ph_R,VR.x(),VR.y(),VR.z()};
  for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
  TMatrixD Meas0(nMeas,1,meas);  
  TMatrixD Unkn0(nUnkn,1,unkn);
  vector<double>Pull;
  Pull.resize(nMeas);
  vector<double>UPull;
  UPull.resize(nUnkn);
  Measurements.push_back(Meas0);
  Unknowns.push_back(Unkn0);
  Pulls.push_back(Pull);
  UPulls.push_back(UPull);
  Chi2s.push_back(-1);
  MassDiffs.push_back(1e9);
  VtxCov.ResizeTo(3,3);
  StepVtxCov.push_back(TMatrixD(3,3));
}
void MassVertexFitter2::SetConstraints(){
// Loading Variables...
  auto Meas = Measurements.at(step); 
  auto Unkn = Unknowns.at(step);
  double p_R=  Unkn(0,0); 
  double th_R= Unkn(1,0); 
  double ph_R= Unkn(2,0);

  double p_P=  Meas(0,0); 
  double th_P= Meas(1,0); 
  double ph_P= Meas(2,0);
  double cx_P=  Meas(3,0);
  double cy_P=  Meas(4,0);
  double z0_P=  Meas(5,0);

  double p_Q=  Meas(6,0); 
  double th_Q= Meas(7,0); 
  double ph_Q= Meas(8,0); 
  double cx_Q=  Meas(9,0);
  double cy_Q=  Meas(10,0);
  double z0_Q=  Meas(11,0);

  vector<double> meas_vect= {p_P,th_P,ph_P,cx_P,cy_P,z0_P,p_Q,th_Q,ph_Q,cx_Q,cy_Q,z0_Q};
  
  // Constraints
  double f1 = 
    -p_R*sin(th_R)*cos(ph_R) 
    +p_P*sin(th_P)*cos(ph_P) 
    +p_Q*sin(th_Q)*cos(ph_Q) ;//Constraint on x momentum
  double f2 = 
    -p_R*sin(th_R)*sin(ph_R) 
    +p_P*sin(th_P)*sin(ph_P) 
    +p_Q*sin(th_Q)*sin(ph_Q) ;//Constraint on y momentum
  double f3 =  
    -p_R*cos(th_R) 
    +p_P*cos(th_P)
    +p_Q*cos(th_Q);//Constraint on z momentum 
  double f4  =
    - sqrt(p_R*p_R+mR*mR)
    + sqrt(p_P*p_P+mP*mP)
    + sqrt(p_Q*p_Q+mQ*mQ);//Constraint on Energy
  double f5 = CalcHelixDistance2(meas_vect);//Constraint on vertex (helix distance)


  // Jacobians

  //f1 derivatives//
  double df1dp_R = -sin(th_R)*cos(ph_R);//df1 / d(P_R)
  double df1dth_R = -p_R*cos(th_R)*cos(ph_R);//It could be df1/dm1 in case of 3-C fit.However, I didnt want to change the token... Mathematically it should be df1 / d (Th_R)
  double df1dph_R = p_R*sin(th_R)*sin(ph_R);//df1 / d(Ph_R)
	// Vertex distance is not explictly related to kinematics
  
  double df1dp_P = sin(th_P)*cos(ph_P);//...
  double df1dth_P = p_P*cos(th_P)*cos(ph_P);// d/ dth_P
  double df1dph_P = -p_P*sin(th_P)*sin(ph_P);
	double df1dcx_P = 0;
	double df1dcy_P = 0;
	double df1dz0_P = 0;
 
  double df1dp_Q = sin(th_Q)*cos(ph_Q);
  double df1dth_Q = p_Q*cos(th_Q)*cos(ph_Q);// d/ dth_Q
  double df1dph_Q = -p_Q*sin(th_Q)*sin(ph_Q);
  double df1dcx_Q = 0;
  double df1dcy_Q = 0;
  double df1dz0_Q = 0;
  
  //f2 derivatives//
  double df2dp_R = -sin(th_R)*sin(ph_R);
  double df2dth_R = -p_R*cos(th_R)*sin(ph_R);// d/ dth_R
  double df2dph_R = -p_R*sin(th_R)*cos(ph_R);

  double df2dp_P = sin(th_P)*sin(ph_P);
  double df2dth_P = p_P*cos(th_P)*sin(ph_P);// d/ dth_P
  double df2dph_P = p_P*sin(th_P)*cos(ph_P);
  double df2dcx_P = 0;
	double df2dcy_P = 0;
	double df2dz0_P = 0;

  double df2dp_Q = sin(th_Q)*sin(ph_Q);
  double df2dth_Q = p_Q*cos(th_Q)*sin(ph_Q);// d/ dth_Q
  double df2dph_Q = p_Q*sin(th_Q)*cos(ph_Q);
	double df2dcx_Q = 0;
	double df2dcy_Q = 0;
	double df2dz0_Q = 0;

  //f3 derivatives//
  double df3dp_R = -cos(th_R);
  double df3dth_R = p_R*sin(th_R);
  double df3dph_R = 0;

  double df3dp_P = cos(th_P);
  double df3dth_P = -p_P*sin(th_P);
  double df3dph_P = 0;
  double df3dcx_P = 0;
	double df3dcy_P = 0;
	double df3dz0_P = 0;

  double df3dp_Q = cos(th_Q);
  double df3dth_Q = -p_Q*sin(th_Q);
  double df3dph_Q = 0;
	double df3dcx_Q = 0;
	double df3dcy_Q = 0;
	double df3dz0_Q = 0;

  //f4 derivatives//
  double ER = sqrt(p_R*p_R+mR*mR);
  double df4dp_R = -p_R/ER;
  double df4dth_R = 0;
  double df4dph_R = 0;

  double df4dp_P = p_P/sqrt(p_P*p_P+mP*mP);
  double df4dth_P = 0;
  double df4dph_P = 0;
  double df4dcx_P = 0;
  double df4dcy_P = 0;
  double df4dz0_P = 0;
  
  double df4dp_Q = p_Q/sqrt(p_Q*p_Q+mQ*mQ);
  double df4dth_Q = 0;
  double df4dph_Q = 0;
  double df4dcx_Q = 0;
  double df4dcy_Q = 0;
  double df4dz0_Q = 0;




  //f5 derivatives//
  double df5dp_R = 0; 
  double df5dth_R = 0; 
  double df5dph_R = 0;

  double df5dp_P = DerivativeHelixDistance2(0, meas_vect);
  double df5dth_P = DerivativeHelixDistance2(1, meas_vect);
  double df5dph_P = DerivativeHelixDistance2(2, meas_vect);
  double df5dcx_P = DerivativeHelixDistance2(3, meas_vect);
  double df5dcy_P = DerivativeHelixDistance2(4, meas_vect);
  double df5dz0_P = DerivativeHelixDistance2(5, meas_vect);
  double df5dp_Q = DerivativeHelixDistance2(6, meas_vect);
  double df5dth_Q = DerivativeHelixDistance2(7, meas_vect);
  double df5dph_Q = DerivativeHelixDistance2(8, meas_vect);
  double df5dcx_Q = DerivativeHelixDistance2(9, meas_vect);
  double df5dcy_Q = DerivativeHelixDistance2(10, meas_vect);
  double df5dz0_Q = DerivativeHelixDistance2(11, meas_vect);





  double df1du1du1 = 0;
  double df1du1du2 = -cos(th_R)*cos(ph_R);
  double df1du1du3 = sin(th_R)*sin(ph_R);
  
  double df1du2du1 = df1du1du2;
  double df1du2du2 = p_R*sin(th_R)*cos(ph_R);
  double df1du2du3 = p_R*cos(th_R)*sin(ph_R);
  
  double df1du3du1 = df1du1du3;
  double df1du3du2 = df1du2du3;
  double df1du3du3 = p_R*sin(th_R)*cos(ph_R);

  double df2du1du1 = 0;
  double df2du1du2 = -cos(th_R)*sin(ph_R);
  double df2du1du3 = -sin(th_R)*cos(ph_R);

  double df2du2du1 = df2du1du2;
  double df2du2du2 = p_R*sin(th_R)*sin(ph_R);
  double df2du2du3 = -p_R*cos(th_R)*cos(ph_R);
  
  double df2du3du1 = df2du1du3;
  double df2du3du2 = df2du2du3;
  double df2du3du3 = p_R*sin(th_R)*sin(ph_R);
  double df3du1du1 = 0;
  double df3du1du2 = sin(th_R);
  double df3du1du3 = 0;
  
  double df3du2du1 = df3du1du2;
  double df3du2du2 = p_R*cos(th_R);
  double df3du2du3 = 0;

  double df3du3du1 = 0;
  double df3du3du2 = 0;
  double df3du3du3 = 0;
  
  double df4du1du1 = -mR*mR/ER/ER/ER;
  double df4du1du2 = 0;
  double df4du1du3 = 0;
  
  double df4du2du1 = 0;
  double df4du2du2 = 0;
  double df4du2du3 = 0;

  double df4du3du1 = 0;
  double df4du3du2 = 0;
  double df4du3du3 = 0;



  double fs[]={f1,f2,f3,f4,f5};
  double dfdms[200] ;
  double dfdus[200] ;
  double temp[] = {// 12 meas, 5 const -> 12 X 5 matrix
    df1dp_P, df1dth_P, df1dph_P, df1dcx_P, df1dcy_P, df1dz0_P, df1dp_Q, df1dth_Q, df1dph_Q, df1dcx_Q, df1dcy_Q, df1dz0_Q,
    df2dp_P, df2dth_P, df2dph_P, df2dcx_P, df2dcy_P, df2dz0_P, df2dp_Q, df2dth_Q, df2dph_Q, df2dcx_Q, df2dcy_Q, df2dz0_Q,
    df3dp_P, df3dth_P, df3dph_P, df3dcx_P, df3dcy_P, df3dz0_P, df3dp_Q, df3dth_Q, df3dph_Q, df3dcx_Q, df3dcy_Q, df3dz0_Q,
    df4dp_P, df4dth_P, df4dph_P, df4dcx_P, df4dcy_P, df4dz0_P, df4dp_Q, df4dth_Q, df4dph_Q, df4dcx_Q, df4dcy_Q, df4dz0_Q,
    df5dp_P, df5dth_P, df5dph_P, df5dcx_P, df5dcy_P, df5dz0_P, df5dp_Q, df5dth_Q, df5dph_Q, df5dcx_Q, df5dcy_Q, df5dz0_Q 
  };
  for(int i=0;i<nMeas*nConst;++i){
    dfdms[i]=temp[i];
  };
  double tempu[] = {
    df1dp_R, df1dth_R, df1dph_R,
    df2dp_R, df2dth_R, df2dph_R,
    df3dp_R, df3dth_R, df3dph_R,
    df4dp_R, df4dth_R, df4dph_R,
    df5dp_R, df5dth_R, df5dph_R
  };
  for(int i=0;i<nUnkn*nConst;++i){
    dfdus[i]=tempu[i];
  };

  //Hessian : Not supported yet
  double temp1[] = {
    df1du1du1,df1du1du2,df1du1du3,
    df1du2du1,df1du2du2,df1du2du3,
    df1du3du1,df1du3du2,df1du3du3
  };
  double temp2[] = {
    df2du1du1,df2du1du2,df2du1du3,
    df2du2du1,df2du2du2,df2du2du3,
    df2du3du1,df2du3du2,df2du3du3
  };
  double temp3[] = {
    df3du1du1,df3du1du2,df3du1du3,
    df3du2du1,df3du2du2,df3du2du3,
    df3du3du1,df3du3du2,df3du3du3
  };
  double temp4[] = {
    df4du1du1,df4du1du2,df4du1du3,
    df4du2du1,df4du2du2,df4du2du3,
    df4du3du1,df4du3du2,df4du3du3
  };
  TMatrixD d2F1dU(nUnkn,nUnkn,temp1);
  TMatrixD d2F2dU(nUnkn,nUnkn,temp2);
  TMatrixD d2F3dU(nUnkn,nUnkn,temp3);
  TMatrixD d2F4dU(nUnkn,nUnkn,temp4);
  TMatrixD d2F5dU(nUnkn,nUnkn,temp4);
  vector<TMatrixD> d2FdU = {
    d2F1dU,
    d2F2dU,
    d2F3dU,
    d2F4dU,
    d2F5dU
  };
  d2Fd2Us.push_back(d2FdU);
  // Hessian //


  



  TMatrixD FMat(nConst,1,fs);
#if Debug
  cout<<"Constraint";
  FMat.Print();
#endif
  TMatrixD dFdM(nConst,nMeas,dfdms);
  TMatrixD dFdU(nConst,nUnkn,dfdus);
  FMats.push_back(FMat);//Constraint Matrices for Each step
  dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
  dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void MassVertexFitter2::SampleStepPoint(int steps){
  auto Meas = Measurements.at(steps); 
  auto Unkn = Unknowns.at(steps); 
  double p_R= Unkn(0,0); 
  double th_R= Unkn(1,0); 
  double ph_R= Unkn(2,0);
  
  double p_P=  Meas(0,0); 
  double th_P= Meas(1,0); 
  double ph_P= Meas(2,0);
  double cx_P=  Meas(3,0);
  double cy_P=  Meas(4,0);
  double z0_P=  Meas(5,0);
  double p_Q=  Meas(6,0); 
  double th_Q= Meas(7,0); 
  double ph_Q= Meas(8,0);
  double cx_Q=  Meas(9,0);
  double cy_Q=  Meas(10,0);
  double z0_Q=  Meas(11,0);

  double px_P = p_P*sin(th_P)*cos(ph_P);
  double py_P = p_P*sin(th_P)*sin(ph_P);
  double pz_P = p_P*cos(th_P);
  double px_Q = p_Q*sin(th_Q)*cos(ph_Q);
  double py_Q = p_Q*sin(th_Q)*sin(ph_Q);
  double pz_Q = p_Q*cos(th_Q);
  double px_R = p_R*sin(th_R)*cos(ph_R);
  double py_R = p_R*sin(th_R)*sin(ph_R);
  double pz_R = p_R*cos(th_R);
  TLorentzVector PP(px_P,py_P,pz_P,hypot(mP,p_P));
  TLorentzVector QQ(px_Q,py_Q,pz_Q,hypot(mQ,p_Q));
  TLorentzVector RR(px_R,py_R,pz_R,hypot(mR,p_R));
  auto PPQQ = PP + QQ;
  double MassDiff = PPQQ.Mag()-mR;
  PCor = PP;
  QCor = QQ;
  RCor = RR;
  TVector3 C_P(cx_P,cy_P,z0_P);
  TVector3 C_Q(cx_Q,cy_Q,z0_Q);
  TVector3 VPCor,VQCor;

  VPCor = CalcHelixPosition(p_P, th_P, ph_P, C_P, charge_config==0 ? 1 : -1);
  VQCor = CalcHelixPosition(p_Q, th_Q, ph_Q, C_Q, charge_config==0 ? -1 : 1);
  VRCor = 0.5*(VPCor + VQCor);
  MassDiffs.push_back(MassDiff);


  /*
  Vtx Covariance estimation:
  From the helix positons xp,yp,zp, xq,yq,zq, the resulting vtx is (xp+xq)/2, (yp+yq)/2, (zp+zq)/2.
  Since the position X is a function of the measured parameters m, the variance of X can be calculated as
  V_X = dX_dm * V * dX_dm^T.
  since m is a 12 dim vector, but since there is no explicict dependency on
   Q variables to P vertex position, dV_p/dm_6 .. 11 == 0 and vise versa. 
  */
  TMatrixD dX_dm(nUnkn, nMeas);
  for(int i=0;i<12;++i){
    double dp = 1e-5;
    double meas_up[12];
    // 
    dX_dm(0,i) = 0.5*(i<6?DerivativeHelixPosition(i%6,p_P, th_P, ph_P, C_P, charge_config==0? 1:-1 ).X()
                          :DerivativeHelixPosition(i%6,p_Q, th_Q, ph_Q, C_Q, charge_config==0? -1:1 ).X());
    dX_dm(1,i) = 0.5*(i<6?DerivativeHelixPosition(i%6,p_P, th_P, ph_P, C_P, charge_config==0? 1:-1 ).Y()
                          :DerivativeHelixPosition(i%6,p_Q, th_Q, ph_Q, C_Q, charge_config==0? -1:1 ).Y());
    dX_dm(2,i) = 0.5*(i<6?DerivativeHelixPosition(i%6,p_P, th_P, ph_P, C_P, charge_config==0? 1:-1 ).Z()
                          :DerivativeHelixPosition(i%6,p_Q, th_Q, ph_Q, C_Q, charge_config==0? -1:1 ).Z());
  }
  TMatrixD dX_dmT = TransposeMatrix(dX_dm);
  TMatrixD V_m = Variancies.at(0);
  TMatrixD V_X = dX_dm * V_m * dX_dmT;
  VtxCov = V_X;
  StepVtxCov.push_back(V_X);

}
TMatrixD
MassVertexFitter2::JacobianSphToCart(double p, double th, double ph){
  // x = p sin(th) cos(ph)
  // y = p sin(th) sin(ph)
  // z = p cos(th)
  //V_c = J^T V J|->
  //    dxdp, dxdth,dxdph
  //J  =  dydp, dydth,dydph
  //    dzdp, dzdth,dzdph

  double dxdp = sin(th)*cos(ph);
  double dydp = sin(th)*sin(ph);
  double dzdp = cos(th);

  double dxdth = p*cos(th)*cos(ph);
  double dydth = p*cos(th)*sin(ph);
  double dzdth = -p*sin(th);

  double dxdph = -p*sin(th)*sin(ph);
  double dydph = p*sin(th)*cos(ph);
  double dzdph = 0;
  double mat[9] = 
  { dxdp, dxdth, dxdph,
    dydp, dydth, dydph,
    dzdp, dzdth, dzdph
  };
  /*
  double mat[9] = 
  { dxdp, dydp, dzdp,
    dxdth, dydth, dzdth,
    dxdph, dydph, dzdph
  };
  */
  return TMatrixD(3,3,mat);

}
void
MassVertexFitter2::CalcVariance(int istep){
  //Not supproted yet
  /*
  auto Meas = Measurements.at(istep); 
  auto Unkn = Unknowns.at(istep);
  double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q;
  p_R= Unkn(0,0); 
  th_R= Unkn(1,0); 
  ph_R= Unkn(2,0); 
  p_P= Meas(0,0); 
  th_P= Meas(1,0); 
  ph_P= Meas(2,0);
  p_Q= Meas(3,0); 
  th_Q= Meas(4,0); 
  ph_Q= Meas(5,0); 
  TMatrixD Jsc_P = JacobianSphToCart(p_P,th_P,ph_P);
  TMatrixD Jsc_Q = JacobianSphToCart(p_Q,th_Q,ph_Q);
  
  double El_Jsc_PQ[6*6]= {0};
  for(int ic =0;ic<3;++ic){
  for(int ir =0;ir<3;++ir){
    int col_P = ic, row_P = ir;
    int col_Q = ic+3, row_Q = ir+3;
    El_Jsc_PQ[row_P+6*col_P] = Jsc_P(ic,ir);
    El_Jsc_PQ[row_Q+6*col_Q] = Jsc_Q(ic,ir);
  }
  }
  TMatrixD Jsc_PQ = TMatrixD(6,6,El_Jsc_PQ);
#if Debug > 1
  Jsc_P.Print();
  Jsc_Q.Print();
  Jsc_PQ.Print();
  cin.ignore();
#endif
  TMatrixD Jsc_PQ_T = TransposeMatrix(Jsc_PQ);
  TMatrixD VMat = Variancies.at(istep);
  TMatrixD dV = dVMats.at(istep);
  TMatrixD VMat_C = Jsc_PQ_T*(VMat-dV)*Jsc_PQ;
//  TMatrixD VMat_C = Jsc_PQ_T*(VMat)*Jsc_PQ;
  double El_contract[18]={//reduce matrix dimension
    1,0,0,1,0,0,
    0,1,0,0,1,0,
    0,0,1,0,0,1
  };
  TMatrixD ContT(3,6,El_contract);
  TMatrixD Cont = TransposeMatrix(ContT);
  TMatrixD UVMat_C = ContT*VMat_C*Cont;
  TMatrixD Jcs_R = JacobianSphToCart(p_R,th_R,ph_R);
  Jcs_R.Invert();
  auto Jcs_RT = TransposeMatrix(Jcs_R);
  TMatrixD UVMat = Jcs_RT*UVMat_C*Jcs_R;
//  VarianciesU.push_back(UVMat);
  */
}

void
MassVertexFitter2::Rotate(){
  auto VMat = Variancies.at(0);
  Initialize();
  Variancies.push_back(VMat);
  TMatrixD J;
  RotateVariance(J);
}
void
MassVertexFitter2::ToDecayPlane(){
  auto Zaxis =(P + Q).Vect();
  auto vP = P.Vect();
  auto vQ = Q.Vect();
  auto Yaxis = vP.Cross(vQ);
//  double YNorm = 1./(Yaxis.Mag());
//  Yaxis = YNorm * Yaxis;
  double Th_F = Zaxis.Theta();
  double Ph_F = Zaxis.Phi();
  double RotZ[9] ={
    cos(-Ph_F),  -sin(-Ph_F),  0,  
    sin(-Ph_F),  cos(-Ph_F),    0,
    0,          0,            1
  };
  double RotY[9] ={
    cos(Th_F),  0,          -sin(Th_F),
    0,          -1,          0,
    sin(Th_F),  0,          cos(Th_F)
  };
  TMatrixD RZ(3,3,RotZ);
  TMatrixD RY(3,3,RotY);
  TMatrixD R_F = RY * RZ;
  Yaxis = R_F * Yaxis;
  double Th_Y = Yaxis.Theta();
  double Ph_Y = Yaxis.Phi();
  double RotX[9] ={
    1,        0,          0,
    0,        cos(-Ph_Y),  -sin(-Ph_Y),
    0,        sin(-Ph_Y),  cos(-Ph_Y)
  };
  TMatrixD RX(3,3,RotX);
  Yaxis = RX * Yaxis;

}
#endif
