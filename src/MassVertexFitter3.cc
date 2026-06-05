//#include "KinFit.cc"
#include "MassVertexFitter3.hh"
#include "TString.h"
#ifndef MassVertexFitter3_cc
#define MassVertexFitter3_cc
#define Debug 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit


MassVertexFitter3::MassVertexFitter3(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_
                                 ,TVector3 V_P, TVector3 V_Q){ 
  /*
  Kinematic fitting for R -> P+Q decay, with vertex constraint
   Note that, the relationship between the kinematic parameter 
  and vertex position should be defined in the off-diagonal term
  in the covariance matrix. These off-diagonal terms should 'drag'
  the vertex position according to the change of the kinematic parameters.
  Without these terms, the fitter will simply return the 1-C fit result,
  with P Q vertices just moved to the averaged vertex.
  */
  P=P_;
  Q=Q_;
  R=R_;
  VP=V_P;
  VQ=V_Q;
  VR=0.5*(V_P+V_Q);
  ScaleParams = 0;//Set to false.
  Initialize();
};

void MassVertexFitter3::Initialize(){
  nMeas = 12;
  // p, th, ph, vx, vy, vz for both particles.
  nUnkn = 6;// p, th, ph, vx, vy, vz
  nConst = 10;// px, py, pz conservations, E conservation(Mass constraint), dx, dy, dz = 0 (vertex constraint) 
  //ndf = nConst - nUnkn;
  mP = P.Mag();
  TVector3 TV_P = P.Vect();
  double p_P = TV_P.Mag();  
  double th_P = TV_P.Theta();
  double ph_P = TV_P.Phi();
  double vx_P = VP.x();
  double vy_P = VP.y();
  double vz_P = VP.z();



  mQ = Q.Mag();
  TVector3 TV_Q = Q.Vect(); 
  double p_Q = TV_Q.Mag();  
  double th_Q = TV_Q.Theta();
  double ph_Q = TV_Q.Phi();
  double vx_Q = VQ.x();
  double vy_Q = VQ.y();
  double vz_Q = VQ.z();

  TVector3 TV_R = R.Vect(); 
  double p_R = TV_R.Mag();  
  double th_R = TV_R.Theta();
  double ph_R = TV_R.Phi();
  double vx_R = VR.x();
  double vy_R = VR.y();
  double vz_R = VR.z();

  double meas[12];
  double unkn[6];
  double temp[] = {p_P,th_P,ph_P,vx_P,vy_P,vz_P,p_Q,th_Q,ph_Q,vx_Q,vy_Q,vz_Q};
  for(int i=0;i<nMeas;++i)meas[i]=temp[i];
  double temp2[] = {p_R,th_R,ph_R,vx_R,vy_R,vz_R};
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
void MassVertexFitter3::SetConstraints(){
// Loading Variables...
  auto Meas = Measurements.at(step); 
  auto Unkn = Unknowns.at(step);
  double p_R=  Unkn(0,0); 
  double th_R= Unkn(1,0); 
  double ph_R= Unkn(2,0);
  double vx_R= Unkn(3,0);
  double vy_R= Unkn(4,0);
  double vz_R= Unkn(5,0);

  double p_P=  Meas(0,0); 
  double th_P= Meas(1,0); 
  double ph_P= Meas(2,0);
  double vx_P=  Meas(3,0);
  double vy_P=  Meas(4,0);
  double vz_P=  Meas(5,0);

  double p_Q=  Meas(6,0); 
  double th_Q= Meas(7,0); 
  double ph_Q= Meas(8,0); 
  double vx_Q=  Meas(9,0);
  double vy_Q=  Meas(10,0);
  double vz_Q=  Meas(11,0);

  double dist_PQ = hypot(vx_P-vx_Q,vy_P-vy_Q,vz_P-vz_Q);
  vector<double> meas_vect= {p_P,th_P,ph_P,vx_P,vy_P,vz_P,p_Q,th_Q,ph_Q,vx_Q,vy_Q,vz_Q};
  
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
  double f5 = -vx_R + 0.5*(vx_P + vx_Q);
  double f6 = -vy_R + 0.5*(vy_P + vy_Q);
  double f7 = -vz_R + 0.5*(vz_P + vz_Q);
	double f8 = vx_P - vx_Q;
	double f9 = vy_P - vy_Q;
	double f10= vz_P - vz_Q;

  // Jacobians

  //f1 derivatives//
  double df1dp_R = -sin(th_R)*cos(ph_R);//df1 / d(P_R)
  double df1dth_R = -p_R*cos(th_R)*cos(ph_R);//It could be df1/dm1 in case of 3-C fit.However, I didnt want to change the token... Mathematically it should be df1 / d (Th_R)
  double df1dph_R = p_R*sin(th_R)*sin(ph_R);//df1 / d(Ph_R)
  double df1dvx_R = 0;double df1dvy_R = 0;double df1dvz_R = 0;
	// Vertex distance is not explictly related to kinematics
  
  double df1dp_P = sin(th_P)*cos(ph_P);//...
  double df1dth_P = p_P*cos(th_P)*cos(ph_P);// d/ dth_P
  double df1dph_P = -p_P*sin(th_P)*sin(ph_P);
	double df1dvx_P = 0;double df1dvy_P = 0;double df1dvz_P = 0;
 
  double df1dp_Q = sin(th_Q)*cos(ph_Q);
  double df1dth_Q = p_Q*cos(th_Q)*cos(ph_Q);// d/ dth_Q
  double df1dph_Q = -p_Q*sin(th_Q)*sin(ph_Q);
  double df1dvx_Q = 0;double df1dvy_Q = 0;double df1dvz_Q = 0;
  
  //f2 derivatives//
  double df2dp_R = -sin(th_R)*sin(ph_R);
  double df2dth_R = -p_R*cos(th_R)*sin(ph_R);// d/ dth_R
  double df2dph_R = -p_R*sin(th_R)*cos(ph_R);
  double df2dvx_R = 0;double df2dvy_R = 0;double df2dvz_R = 0;

  double df2dp_P = sin(th_P)*sin(ph_P);
  double df2dth_P = p_P*cos(th_P)*sin(ph_P);// d/ dth_P
  double df2dph_P = p_P*sin(th_P)*cos(ph_P);
  double df2dvx_P = 0;double df2dvy_P = 0;double df2dvz_P = 0;

  double df2dp_Q = sin(th_Q)*sin(ph_Q);
  double df2dth_Q = p_Q*cos(th_Q)*sin(ph_Q);// d/ dth_Q
  double df2dph_Q = p_Q*sin(th_Q)*cos(ph_Q);
	double df2dvx_Q = 0;double df2dvy_Q = 0;double df2dvz_Q = 0;

  //f3 derivatives//
  double df3dp_R = -cos(th_R);
  double df3dth_R = p_R*sin(th_R);
  double df3dph_R = 0;
  double df3dvx_R = 0;
  double df3dvy_R = 0;
  double df3dvz_R = 0;

  double df3dp_P = cos(th_P);
  double df3dth_P = -p_P*sin(th_P);
  double df3dph_P = 0;
  double df3dvx_P = 0;double df3dvy_P = 0;double df3dvz_P = 0;

  double df3dp_Q = cos(th_Q);
  double df3dth_Q = -p_Q*sin(th_Q);
  double df3dph_Q = 0;
	double df3dvx_Q = 0;double df3dvy_Q = 0;double df3dvz_Q = 0;

  //f4 derivatives//
  double ER = sqrt(p_R*p_R+mR*mR);
  double df4dp_R = -p_R/ER;
  double df4dth_R = 0;
  double df4dph_R = 0;
  double df4dvx_R = 0; double df4dvy_R = 0; double df4dvz_R = 0;

  double df4dp_P = p_P/sqrt(p_P*p_P+mP*mP);
  double df4dth_P = 0;
  double df4dph_P = 0;
  double df4dvx_P = 0;double df4dvy_P = 0;double df4dvz_P = 0;
  
  double df4dp_Q = p_Q/sqrt(p_Q*p_Q+mQ*mQ);
  double df4dth_Q = 0;
  double df4dph_Q = 0;
  double df4dvx_Q = 0;double df4dvy_Q = 0;double df4dvz_Q = 0;

  //f5 derivatives//
  double df5dp_R = 0; double df5dth_R = 0; double df5dph_R = 0;
  double df5dvx_R = -1;double df5dvy_R = 0;double df5dvz_R = 0;

  double df5dp_P = 0;  double df5dth_P = 0;double df5dph_P = 0;
  double df5dvx_P =0.5;double df5dvy_P = 0;double df5dvz_P = 0;

  double df5dp_Q = 0;   double df5dth_Q = 0;double df5dph_Q = 0;
  double df5dvx_Q = 0.5;double df5dvy_Q = 0;double df5dvz_Q = 0;
  
  //f6 derivatives//
  double df6dp_R = 0; double df6dth_R = 0; double df6dph_R = 0;
  double df6dvx_R = 0;double df6dvy_R = 1;double df6dvz_R = 0;

  double df6dp_P = 0;double df6dth_P = 0;  double df6dph_P = 0;
  double df6dvx_P = 0;double df6dvy_P =0.5;double df6dvz_P = 0;

  double df6dp_Q = 0;double df6dth_Q = 0;  double df6dph_Q = 0;
  double df6dvx_Q = 0;double df6dvy_Q =0.5;double df6dvz_Q = 0;

  //f7 derivatives//
  double df7dp_R = 0; double df7dth_R = 0; double df7dph_R = 0;
  double df7dvx_R = 0;double df7dvy_R = 0;double df7dvz_R = -1;

  double df7dp_P = 0;double df7dth_P = 0;double df7dph_P = 0;
  double df7dvx_P = 0;double df7dvy_P = 0;double df7dvz_P = 0.5;
  double df7dp_Q = 0;double df7dth_Q = 0;double df7dph_Q = 0;
  double df7dvx_Q = 0;double df7dvy_Q = 0;double df7dvz_Q = 0.5;


  double df8dp_R = 0;double df8dth_R = 0;double df8dph_R = 0;
  double df8dvx_R = 0;double df8dvy_R = 0;double df8dvz_R = 0;

  double df8dp_P = 0; double df8dth_P = 0;double df8dph_P = 0;
  double df8dvx_P = 1;double df8dvy_P = 0;double df8dvz_P = 0;
  double df8dp_Q = 0; double df8dth_Q = 0;double df8dph_Q = 0;
  double df8dvx_Q =-1;double df8dvy_Q = 0;double df8dvz_Q = 0;


  double df9dp_R = 0;double df9dth_R = 0;double df9dph_R = 0;
  double df9dvx_R = 0;double df9dvy_R = 0;double df9dvz_R = 0;

  double df9dp_P = 0; double df9dth_P = 0;double df9dph_P = 0;
  double df9dvx_P = 0;double df9dvy_P = 1;double df9dvz_P = 0;
  double df9dp_Q = 0; double df9dth_Q = 0;double df9dph_Q = 0;
  double df9dvx_Q = 0;double df9dvy_Q =-1;double df9dvz_Q = 0;


  double df10dp_R = 0;double df10dth_R = 0;double df10dph_R = 0;
  double df10dvx_R = 0;double df10dvy_R = 0;double df10dvz_R = 0;
  double df10dp_P = 0; double df10dth_P = 0;double df10dph_P = 0;
  double df10dvx_P = 0;double df10dvy_P = 0;double df10dvz_P = 1;
  double df10dp_Q = 0; double df10dth_Q = 0;double df10dph_Q = 0;
  double df10dvx_Q = 0;double df10dvy_Q = 0;double df10dvz_Q =-1;



  double fs[]={f1,f2,f3,f4,f5,f6,f7,f8};
  double dfdms[200];
  double dfdus[200];
  double temp[] = {// 12 meas, 10 const -> 12 X 10 matrix
    df1dp_P,  df1dth_P,  df1dph_P,  df1dvx_P,  df1dvy_P,  df1dvz_P,  df1dp_Q,  df1dth_Q,  df1dph_Q,  df1dvx_Q,  df1dvy_Q,  df1dvz_Q,
    df2dp_P,  df2dth_P,  df2dph_P,  df2dvx_P,  df2dvy_P,  df2dvz_P,  df2dp_Q,  df2dth_Q,  df2dph_Q,  df2dvx_Q,  df2dvy_Q,  df2dvz_Q,
    df3dp_P,  df3dth_P,  df3dph_P,  df3dvx_P,  df3dvy_P,  df3dvz_P,  df3dp_Q,  df3dth_Q,  df3dph_Q,  df3dvx_Q,  df3dvy_Q,  df3dvz_Q,
    df4dp_P,  df4dth_P,  df4dph_P,  df4dvx_P,  df4dvy_P,  df4dvz_P,  df4dp_Q,  df4dth_Q,  df4dph_Q,  df4dvx_Q,  df4dvy_Q,  df4dvz_Q,
    df5dp_P,  df5dth_P,  df5dph_P,  df5dvx_P,  df5dvy_P,  df5dvz_P,  df5dp_Q,  df5dth_Q,  df5dph_Q,  df5dvx_Q,  df5dvy_Q,  df5dvz_Q, 
    df6dp_P,  df6dth_P,  df6dph_P,  df6dvx_P,  df6dvy_P,  df6dvz_P,  df6dp_Q,  df6dth_Q,  df6dph_Q,  df6dvx_Q,  df6dvy_Q,  df6dvz_Q,
    df7dp_P,  df7dth_P,  df7dph_P,  df7dvx_P,  df7dvy_P,  df7dvz_P,  df7dp_Q,  df7dth_Q,  df7dph_Q,  df7dvx_Q,  df7dvy_Q,  df7dvz_Q,
    df8dp_P,  df8dth_P,  df8dph_P,  df8dvx_P,  df8dvy_P,  df8dvz_P,  df8dp_Q,  df8dth_Q,  df8dph_Q,  df8dvx_Q,  df8dvy_Q,  df8dvz_Q,
    df9dp_P,  df9dth_P,  df9dph_P,  df9dvx_P,  df9dvy_P,  df9dvz_P,  df9dp_Q,  df9dth_Q,  df9dph_Q,  df9dvx_Q,  df9dvy_Q,  df9dvz_Q,
    df10dp_P, df10dth_P, df10dph_P, df10dvx_P, df10dvy_P, df10dvz_P, df10dp_Q, df10dth_Q, df10dph_Q, df10dvx_Q, df10dvy_Q, df10dvz_Q 
  };
  for(int i=0;i<nMeas*nConst;++i){
    dfdms[i]=temp[i];
  };
  double tempu[] = {
    df1dp_R,  df1dth_R,  df1dph_R,  df1dvx_R,  df1dvy_R,  df1dvz_R,
    df2dp_R,  df2dth_R,  df2dph_R,  df2dvx_R,  df2dvy_R,  df2dvz_R,
    df3dp_R,  df3dth_R,  df3dph_R,  df3dvx_R,  df3dvy_R,  df3dvz_R,
    df4dp_R,  df4dth_R,  df4dph_R,  df4dvx_R,  df4dvy_R,  df4dvz_R,
    df5dp_R,  df5dth_R,  df5dph_R,  df5dvx_R,  df5dvy_R,  df5dvz_R,
    df6dp_R,  df6dth_R,  df6dph_R,  df6dvx_R,  df6dvy_R,  df6dvz_R,
    df7dp_R,  df7dth_R,  df7dph_R,  df7dvx_R,  df7dvy_R,  df7dvz_R,
    df8dp_R,  df8dth_R,  df8dph_R,  df8dvx_R,  df8dvy_R,  df8dvz_R,
    df9dp_R,  df9dth_R,  df9dph_R,  df9dvx_R,  df9dvy_R,  df9dvz_R,
    df10dp_R, df10dth_R, df10dph_R, df10dvx_R, df10dvy_R, df10dvz_R
  };
  for(int i=0;i<nUnkn*nConst;++i){
    dfdus[i]=tempu[i];
  };

  //Hessian : Not supported yet//
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
void MassVertexFitter3::SampleStepPoint(int steps){
  auto Meas = Measurements.at(steps); 
  auto Unkn = Unknowns.at(steps); 
  double p_R= Unkn(0,0); 
  double th_R= Unkn(1,0); 
  double ph_R= Unkn(2,0);
  double vx_R= Unkn(3,0);
  double vy_R= Unkn(4,0);
  double vz_R= Unkn(5,0);
  
  double p_P=  Meas(0,0); 
  double th_P= Meas(1,0); 
  double ph_P= Meas(2,0);
  double vx_P=  Meas(3,0);
  double vy_P=  Meas(4,0);
  double vz_P=  Meas(5,0);
  double p_Q=  Meas(6,0); 
  double th_Q= Meas(7,0); 
  double ph_Q= Meas(8,0);
  double vx_Q=  Meas(9,0);
  double vy_Q=  Meas(10,0);
  double vz_Q=  Meas(11,0);

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

  VPCor = TVector3(vx_P,vy_P,vz_P);
  VQCor = TVector3(vx_Q,vy_Q,vz_Q);
  VRCor = TVector3(vx_R,vy_R,vz_R);
  MassDiffs.push_back(MassDiff);

}
TMatrixD
MassVertexFitter3::JacobianSphToCart(double p, double th, double ph){
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
MassVertexFitter3::CalcVariance(int istep){
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
MassVertexFitter3::Rotate(){
  auto VMat = Variancies.at(0);
  Initialize();
  Variancies.push_back(VMat);
  TMatrixD J;
  RotateVariance(J);
}
void
MassVertexFitter3::ToDecayPlane(){
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
