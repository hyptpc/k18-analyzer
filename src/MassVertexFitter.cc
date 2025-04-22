//#include "KinFit.cc"
#include "MassVertexFitter.hh"
#include "TString.h"
#ifndef MassVertexFitter_cc
#define MassVertexFitter_cc
#define Debug 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit


MassVertexFitter::MassVertexFitter(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_
                                  ,TVector3 V_P, TVector3 V_Q){ 
  //Kinematic fitting for R -> P+Q decay, with vertex constraint
  P=P_;
  Q=Q_;
  R=R_;
  TVector3 norm = P.Vect().Cross(Q.Vect());
  TVector3 perp = P.Vect().Cross(norm);//Shift the position reference, to avoid crossing-problem at cylindrical coordinate.
  V0 = 0.5*(V_P + V_Q) + 20 * perp;// Position reference of every verticies are from this point. 
  VP = V_P - V0;
  VQ = V_Q - V0;
  Initialize();
};
void
MassVertexFitter::GetRZparameters(TVector3 Vert, TVector3 Dir, double& R, double& Z){
  TVector3 close_point = Vert - (Vert * Dir) * Dir ;
  //track = close_point + t * Dir;. Vert = close_point + t0 * Dir;
  //close_point* Dir = 0, because close point should be perpendicular to the direction of the track.
  //t0 = Vert * Dir;
  //close_point = Vert - t0 * Dir = Vert - (Vert * Dir) * Dir;
  double phi = Dir.Phi();
  R = close_point.x()*cos(phi + M_PI/2) + close_point.y()*sin(phi + M_PI/2);
  Z = close_point.Z();
#if Debug
  cout<<Form("Step %d, Vertex (%g,%g,%g)",step,Vert.X(),Vert.Y(),Vert.Z())<<endl;
  cout<<Form("Direction (%g,%g,%g)",Dir.X(),Dir.Y(),Dir.Z())<<endl;
  cout<<Form("RZ (%g,%g)",R,Z)<<endl;
#endif

}
double
MassVertexFitter::CalcVertexDistance(double th1, double ph1, double r1, double z1,
			double th2, double ph2, double r2, double z2){
  TVector3 dir1(
    sin(th1) * cos(ph1),
    sin(th1) * sin(ph1),
    cos(th1)
  );// Direction of the particle 1 trajectory
  TVector3 dir2(
    sin(th2) * cos(ph2),
    sin(th2) * sin(ph2),
    cos(th2)
  );
  double dph_1 = 0.5 * M_PI;
  double dph_2 = 0.5 * M_PI;

  TVector3 base1(
    r1 * cos(ph1+dph_1),
    r1 * sin(ph1+dph_1),
    z1
  );// Closest point of the particle 1 trajectory, to the origin(V0) with linear assumption.
  TVector3 base2(
    r2 * cos(ph2+dph_2),
    r2 * sin(ph2+dph_2),
    z2
  );// Closest point of the trajectory 2 from the vertex should be at the opposite side of the trajectory 1. 
  TVector3 diff = base1 - base2;
  TVector3 norm = dir1.Cross(dir2);
  double d = norm * diff;
#if Debug
  cout<<Form("Step %d, VertexDistance = %g",step,d)<<endl;
#endif
  return d;
}
void
MassVertexFitter::CalcClosePoint(double th1, double ph1, double r1, double z1,
			double th2, double ph2, double r2, double z2,
      TVector3& vert1, TVector3& vert2){
  double dph_1 = 0.5 * M_PI;
  double dph_2 = 0.5 * M_PI;
  TVector3 base1(
    r1 * cos(ph1+dph_1),
    r1 * sin(ph1+dph_1),
    z1
  );// Closest point of the particle 1 trajectory, to the origin(V0) with linear assumption.
  TVector3 base2(
    r2 * cos(ph2+dph_2),
    r2 * sin(ph2+dph_2),
    z2
  );// Closest point of the trajectory 2 from the vertex should be at the opposite side of the trajectory 1. 
  TVector3 dir1(
    sin(th1) * cos(ph1),
    sin(th1) * sin(ph1),
    cos(th1)
  );// Direction of the particle 1 trajectory
  TVector3 dir2(
    sin(th2) * cos(ph2),
    sin(th2) * sin(ph2),
    cos(th2)
  );

  TVector3 diff = base2 - base1;
  
  double dir_cos  = dir1 * dir2;
  double dir_sin2 = 1 - dir_cos * dir_cos;
  double proj_1 = dir1 * diff;
  double proj_2 = dir2 * diff;
  double t1 = (dir1 * diff - dir2 * diff *dir_cos)/ dir_sin2;
  vert1 = base1 + t1 * dir1;
  double t2 = (-dir2 * diff + dir1*diff*dir_cos)/dir_sin2;
  vert2 = base2 + t2 * dir2;
#if Debug
  cout<<Form("Direction1 (%g,%g,%g)",dir1.X(),dir1.Y(),dir1.Z())<<endl;
  cout<<Form("Vertex1 (%g,%g,%g)",base1.X(),base1.Y(),base1.Z())<<endl;
#endif

}	


void MassVertexFitter::Initialize(){
  nMeas = 10;
  nUnkn = 3;
  nConst = 5;// px, py, pz conservations, E conservation(Mass constraint), Vertex distance  
  //ndf = nMeas - nUnkn - nConst;
  mP = P.Mag();
  TVector3 TV_P = P.Vect();
  double p_P = TV_P.Mag();  
  double th_P = TV_P.Theta();
  double ph_P = TV_P.Phi();
  double r_P,z_P;
  GetRZparameters(VP,TV_P.Unit(),r_P,z_P);
  

  mQ = Q.Mag();
  TVector3 TV_Q = Q.Vect(); 
  double p_Q = TV_Q.Mag();  
  double th_Q = TV_Q.Theta();
  double ph_Q = TV_Q.Phi(); 
  double r_Q,z_Q;
  GetRZparameters(VQ,TV_Q.Unit(),r_Q,z_Q);


  TVector3 TV_R = R.Vect(); 
  double p_R = TV_R.Mag();  
  double th_R = TV_R.Theta();
  double ph_R = TV_R.Phi();

  double meas[10];
  double unkn[3];
  double temp[] = {p_P,th_P,ph_P,r_P,z_P,p_Q,th_Q,ph_Q,r_Q,z_Q};
  for(int i=0;i<nMeas;++i)meas[i]=temp[i];
  double temp2[] = {p_R,th_R,ph_R};
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
}
void MassVertexFitter::SetConstraints(){
// Loading Variables...
  auto Meas = Measurements.at(step); 
  auto Unkn = Unknowns.at(step);
  double p_R=  Unkn(0,0); 
  double th_R= Unkn(1,0); 
  double ph_R= Unkn(2,0); 
  double p_P=  Meas(0,0); 
  double th_P= Meas(1,0); 
  double ph_P= Meas(2,0);
  double r_P=  Meas(3,0);
  double z_P=  Meas(4,0);
  double p_Q=  Meas(5,0); 
  double th_Q= Meas(6,0); 
  double ph_Q= Meas(7,0); 
  double r_Q=  Meas(8,0);
  double z_Q=  Meas(9,0);

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
  double f5 = CalcVertexDistance(th_P,ph_P,r_P,z_P,
    th_Q,ph_Q,r_Q,z_Q);//Vertex constraint


  // Jacobians
  double df1dp_R = -sin(th_R)*cos(ph_R);//df1 / d(P_R)
  double df1dth_R = -p_R*cos(th_R)*cos(ph_R);//It could be df1/dm1 in case of 3-C fit.However, I didnt want to change the token... Mathematically it should be df1 / d (Th_R)
  double df1dph_R = p_R*sin(th_R)*sin(ph_R);//df1 / d(Ph_R)
  double df1dp_P = sin(th_P)*cos(ph_P);//...
  double df1dth_P = p_P*cos(th_P)*cos(ph_P);// d/ dth_P
  double df1dph_P = -p_P*sin(th_P)*sin(ph_P);
  double df1dp_Q = sin(th_Q)*cos(ph_Q);
  double df1dth_Q = p_Q*cos(th_Q)*cos(ph_Q);// d/ dth_Q
  double df1dph_Q = -p_Q*sin(th_Q)*sin(ph_Q);
	// Vertex distance is not explictly related to kinematics
	double df1dr_P = 0;
	double df1dz_P = 0;
	double df1dr_Q = 0;
	double df1dz_Q = 0;
	//
  

  double df2dp_R = -sin(th_R)*sin(ph_R);
  double df2dth_R = -p_R*cos(th_R)*sin(ph_R);// d/ dth_R
  double df2dph_R = -p_R*sin(th_R)*cos(ph_R);
  double df2dp_P = sin(th_P)*sin(ph_P);
  double df2dth_P = p_P*cos(th_P)*sin(ph_P);// d/ dth_P
  double df2dph_P = p_P*sin(th_P)*cos(ph_P);
  double df2dp_Q = sin(th_Q)*sin(ph_Q);
  double df2dth_Q = p_Q*cos(th_Q)*sin(ph_Q);// d/ dth_Q
  double df2dph_Q = p_Q*sin(th_Q)*cos(ph_Q);
  double df2dr_P = 0;
	double df2dz_P = 0;
	double df2dr_Q = 0;
	double df2dz_Q = 0;


  double df3dp_R = -cos(th_R);
  double df3dth_R = p_R*sin(th_R);
  double df3dph_R = 0;
  double df3dp_P = cos(th_P);
  double df3dth_P = -p_Q*sin(th_P);
  double df3dph_P = 0;
  double df3dp_Q = cos(th_Q);
  double df3dth_Q = -p_Q*sin(th_Q);
  double df3dph_Q = 0;
  double df3dr_P = 0;
	double df3dz_P = 0;
	double df3dr_Q = 0;
	double df3dz_Q = 0;


  double ER = sqrt(p_R*p_R+mR*mR);
  double df4dp_R = -p_R/ER;
  double df4dth_R = 0;
  double df4dph_R = 0;
  double df4dp_P = p_P/sqrt(p_P*p_P+mP*mP);
  double df4dth_P = 0;
  double df4dph_P = 0;
  double df4dp_Q = p_Q/sqrt(p_Q*p_Q+mQ*mQ);
  double df4dth_Q = 0;
  double df4dph_Q = 0;
  double df4dr_P = 0;
	double df4dz_P = 0;
	double df4dr_Q = 0;
	double df4dz_Q = 0;


  TVector3 base_P(-r_P * sin(ph_P),r_P * cos(ph_P),z_P);
  TVector3 base_Q(-r_Q * sin(ph_Q),r_Q * cos(ph_Q),z_Q);

  TVector3 dbase_Pdph_P(-r_P * cos(ph_P),-r_P * sin(ph_P),0);
  TVector3 dbase_Pdr_P(-sin(ph_P),cos(ph_P),0);
  TVector3 dbase_Pdz_P(0,0,1);
  TVector3 dbase_Qdph_Q(-r_Q * cos(ph_Q),-r_Q * sin(ph_Q),0);
  TVector3 dbase_Qdr_Q(-sin(ph_Q),cos(ph_Q),0);
  TVector3 dbase_Qdz_Q(0,0,1);

  TVector3 dir_P(sin(th_P) * cos(ph_P),sin(th_P) * sin(ph_P),cos(th_P));
  TVector3 dir_Q(sin(th_Q) * cos(ph_Q),sin(th_Q) * sin(ph_Q),cos(th_Q));

  TVector3 ddir_Pdth_P(cos(th_P) * cos(ph_P),cos(th_P) * sin(ph_P),-sin(th_P));
  TVector3 ddir_Pdph_P(-sin(th_P) * sin(ph_P),sin(th_P) * cos(ph_P),0);
  TVector3 ddir_Qdth_Q(cos(th_Q) * cos(ph_Q),cos(th_Q) * sin(ph_Q),-sin(th_Q));
  TVector3 ddir_Qdph_Q(-sin(th_Q) * sin(ph_Q),sin(th_Q) * cos(ph_Q),0);


  TVector3 norm_PQ = dir_P.Cross(dir_Q);
  //reminder: f5 = (base_P - base_Q) * (dir_P.Cross(dir_Q))
  double df5dp_R  = 0;
	double df5dth_R = 0;
	double df5dph_R = 0;
  double df5dp_P = 0;
	double df5dp_Q = 0; // 
  double df5dth_P = (base_P - base_Q)* (ddir_Pdth_P.Cross(dir_Q));
  double df5dph_P = dbase_Pdph_P*(dir_P.Cross(dir_Q))
                  + (base_P - base_Q)* ddir_Pdph_P.Cross(dir_Q);
  double df5dth_Q = (base_P - base_Q)* (dir_P.Cross(ddir_Qdth_Q));
  double df5dph_Q = dbase_Qdph_Q*(dir_P.Cross(dir_Q))
                  + (base_P - base_Q)* dir_P.Cross(ddir_Qdph_Q);
  double df5dr_P = dbase_Pdr_P*(dir_P.Cross(dir_Q));
  double df5dz_P = dbase_Pdz_P*(dir_P.Cross(dir_Q));
  double df5dr_Q = dbase_Qdr_Q*(dir_P.Cross(dir_Q));
  double df5dz_Q = dbase_Qdz_Q*(dir_P.Cross(dir_Q));



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
  double temp[] = {
    df1dp_P , df1dth_P , df1dph_P , df1dr_P , df1dz_P , df1dp_Q , df1dth_Q , df1dph_Q , df1dr_Q , df1dz_Q,
    df2dp_P , df2dth_P , df2dph_P , df2dr_P , df2dz_P , df2dp_Q , df2dth_Q , df2dph_Q , df2dr_Q , df2dz_Q,
    df3dp_P , df3dth_P , df3dph_P , df3dr_P , df3dz_P , df3dp_Q , df3dth_Q , df3dph_Q , df3dr_Q , df3dz_Q,
    df4dp_P , df4dth_P , df4dph_P , df4dr_P , df4dz_P , df4dp_Q , df4dth_Q , df4dph_Q , df4dr_Q , df4dz_Q,
    df5dp_P , df5dth_P , df5dph_P , df5dr_P , df5dz_P , df5dp_Q , df5dth_Q , df5dph_Q , df5dr_Q , df5dz_Q
  };
  for(int i=0;i<nMeas*nConst;++i){
    dfdms[i]=temp[i];
  };
  double tempu[] = {
    df1dp_R,df1dth_R,df1dph_R,
    df2dp_R,df2dth_R,df2dph_R,
    df3dp_R,df3dth_R,df3dph_R,
    df4dp_R,df4dth_R,df4dph_R,
    df5dp_R,df5dth_R,df5dph_R,
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
void MassVertexFitter::SampleStepPoint(int steps){
  auto Meas = Measurements.at(steps); 
  auto Unkn = Unknowns.at(steps); 
  double p_R= Unkn(0,0); 
  double th_R= Unkn(1,0); 
  double ph_R= Unkn(2,0); 
  double p_P=  Meas(0,0); 
  double th_P= Meas(1,0); 
  double ph_P= Meas(2,0);
  double r_P=  Meas(3,0);
  double z_P=  Meas(4,0);
  double p_Q=  Meas(5,0); 
  double th_Q= Meas(6,0); 
  double ph_Q= Meas(7,0);
  double r_Q=  Meas(8,0);
  double z_Q=  Meas(9,0);

  double Ppx = p_P*sin(th_P)*cos(ph_P);
  double Ppy = p_P*sin(th_P)*sin(ph_P);
  double Ppz = p_P*cos(th_P);
  double Qpx = p_Q*sin(th_Q)*cos(ph_Q);
  double Qpy = p_Q*sin(th_Q)*sin(ph_Q);
  double Qpz = p_Q*cos(th_Q);
  double Rpx = p_R*sin(th_R)*cos(ph_R);
  double Rpy = p_R*sin(th_R)*sin(ph_R);
  double Rpz = p_R*cos(th_R);
  TLorentzVector PP(Ppx,Ppy,Ppz,hypot(mP,p_P));
  TLorentzVector QQ(Qpx,Qpy,Qpz,hypot(mQ,p_Q));
  TLorentzVector RR(Rpx,Rpy,Rpz,hypot(mR,p_R));
  auto V = PP + QQ;
  double MassDiff = V.Mag()-mR;
  PCor = PP;
  QCor = QQ;
  RCor = RR;
  TVector3 VPCor,VQCor;
  CalcClosePoint(th_P,ph_P,r_P,z_P,
    th_Q,ph_Q,r_Q,z_Q,
    VPCor,VQCor);
  MassDiffs.push_back(MassDiff);
}
TMatrixD
MassVertexFitter::JacobianSphToCart(double p, double th, double ph){
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
MassVertexFitter::CalcVariance(int istep){
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
MassVertexFitter::Rotate(){
  auto VMat = Variancies.at(0);
  Initialize();
  Variancies.push_back(VMat);
  TMatrixD J;
  RotateVariance(J);
}
void
MassVertexFitter::ToDecayPlane(){
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
void
MassVertexFitter::GetVertexResolution(double& dr, double& dz ){
  auto VMat = Variancies.at(0);
  dr = sqrt(VMat(3,3) + VMat(8,8)) /2;
  dz = sqrt(VMat(4,4) + VMat(9,9)) /2;
}
#endif
