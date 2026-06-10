#include "MassVertexFitter.hh"
#ifndef MassVertexFitter_cc
#define MassVertexFitter_cc
#define DebugCKF 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void MassVertexFitter::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	UseVertexFlag = status;
	Clear();
	Initialize();
}
MassVertexFitter::MassVertexFitter(TLorentzVector P_,TVector3 VP_,
	TLorentzVector Q_, TVector3 VQ_){
	P=P_;
	Q=Q_;
	VP=VP_;
	VQ=VQ_;
	L=P_+Q_;
	VL=0.5*(VP_+VQ_);
	ScaleParams = 0;
	Initialize();
};
void MassVertexFitter::Initialize(){
	Initialized = 1;
#if DebugCKF
	cout<<"Initializing..."<<endl;
#endif
	Clear();
	if(UseVertexFlag){
		nMeas = 11;nUnkn = 2; nConst = 9; 
	}
	else{
		nMeas = 12;nUnkn = 6; nConst = 10;
		//Meas : px,py,pz,vx,vy,vz for P,Q
		//Unkn : px,py,pz,vx,vy,vz for L
		//Const: px,py,pz conservation(3), E conservation(Mass constraint) for L
		// vertex constraints for P - Q = 0
		// vertex constraints for L - 0.5(P + Q) = 0
		//ndf = nConst - nUnkn = 4
	}
	mP = P.Mag();
	TVector3 TV_P = P.Vect();
	double px_P = TV_P.x();	
	double py_P = TV_P.y();
	double pz_P = TV_P.z();
	double vx_P = VP.x();
	double vy_P = VP.y();
	double vz_P = VP.z();
	
	mQ = Q.Mag();
	TVector3 TV_Q = Q.Vect(); 
	double px_Q = TV_Q.x();
	double py_Q = TV_Q.y();
	double pz_Q = TV_Q.z();
	double vx_Q = VQ.x();
	double vy_Q = VQ.y();
	double vz_Q = VQ.z();


	TVector3 TV_L = L.Vect();
	VL = (VP + VQ)*0.5;//Initial value of Lambda vertex is set to the averaged vertex of P and Q. This is not a constraint, but just an initial value. The fitter will move the vertex according to the constraints and the covariance matrix.
	double px_L = TV_L.x();
	double py_L = TV_L.y();
	double pz_L = TV_L.z();
	double vx_L = VL.x();
	double vy_L = VL.y();
	double vz_L = VL.z();
	TVector3 D_L = TV_L.Unit();

	std::vector<double> MV = {
	px_P,py_P,pz_P,vx_P,vy_P,vz_P,
	px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q
	};
	std::vector<double> UV = {
	px_L,py_L,pz_L,
	vx_L,vy_L,vz_L
	};

	double meas[20];
	double unkn[20];
	if(UseVertexFlag){
	}
	else{
		double temp[] = {
			px_P,py_P,pz_P,vx_P,vy_P,vz_P,
			px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {px_L,py_L,pz_L,vx_L,vy_L,vz_L};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	TMatrixD Meas0(nMeas,1,meas);	
	TMatrixD Unkn0(nUnkn,1,unkn);
#if DebugCKF
	cout<<"Meas0 : ";
	Meas0.Print();
	cout<<"Unkn0 : ";
	Unkn0.Print();
#endif
	vector<double>Pull;
	Pull.resize(nMeas);
	vector<double>UPull;
	UPull.resize(nUnkn);
	Measurements.push_back(Meas0);
	Unknowns.push_back(Unkn0);
	Pulls.push_back(Pull);
	UPulls.push_back(UPull);
	Chi2s.push_back(-1);
	MassDiffsL.push_back(1e9);
	std::cout<<"MassVertexFitter::Initialize() done"<<std::endl;
}
void MassVertexFitter::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double px_P,py_P,pz_P,px_Q,py_Q,pz_Q,px_L,py_L,pz_L;
	double vx_P,vy_P,vz_P,vx_Q,vy_Q,vz_Q,vx_L,vy_L,vz_L;
	if(UseVertexFlag){
	}
	else{
		px_P= Meas(0,0); 
		py_P= Meas(1,0); 
		pz_P= Meas(2,0);
		vx_P= Meas(3,0);
		vy_P= Meas(4,0);
		vz_P= Meas(5,0);

		px_Q= Meas(6,0); 
		py_Q= Meas(7,0); 
		pz_Q= Meas(8,0);
		vx_Q= Meas(9,0);
		vy_Q= Meas(10,0);
		vz_Q= Meas(11,0);

		px_L= Unkn(0,0);
		py_L= Unkn(1,0);
		pz_L= Unkn(2,0);
		vx_L= Unkn(3,0);
		vy_L= Unkn(4,0);
		vz_L= Unkn(5,0);
	}

	TVector3 TV_P(px_P,py_P,pz_P);
	double p_P = TV_P.Mag();
	double E_P = hypot(p_P,mP);
	double dE_Pdpx_P = px_P/E_P;
	double dE_Pdpy_P = py_P/E_P;
	double dE_Pdpz_P = pz_P/E_P;

	TVector3 TV_Q(px_Q,py_Q,pz_Q);
	double p_Q = TV_Q.Mag();
	double E_Q = hypot(p_Q,mQ);
	double dE_Qdpx_Q = px_Q/E_Q;
	double dE_Qdpy_Q = py_Q/E_Q;
	double dE_Qdpz_Q = pz_Q/E_Q;

	TVector3 TV_L(px_L,py_L,pz_L);
	double p_L = TV_L.Mag();
	double E_L = hypot(p_L,mL);
	double dE_Ldpx_L = px_L/E_L;
	double dE_Ldpy_L = py_L/E_L;
	double dE_Ldpz_L = pz_L/E_L;
	
	
	double f1 = -px_L + px_P + px_Q;
	double f2 = -py_L + py_P + py_Q;
	double f3 = -pz_L + pz_P + pz_Q;
	double f4 = -E_L + E_P + E_Q;//Constraint on Lambda Energy
	double f5 = vx_P - vx_Q;
	double f6 = vy_P - vy_Q;
	double f7 = vz_P - vz_Q;
	double f8 = vx_L - 0.5*(vx_P + vx_Q);
	double f9 = vy_L - 0.5*(vy_P + vy_Q);
	double f10= vz_L - 0.5*(vz_P + vz_Q);
	
	std::vector<double> MV = {
		px_P,py_P,pz_P,vx_P,vy_P,vz_P,
		px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q};
	std::vector<double> UV = {
		px_L,py_L,pz_L,
		vx_L,vy_L,vz_L};
	//f1 - f5: Kinematic Constraints//
	//f1: -px_L + px_P + px_Q = 0
	double df1du1 =-1, df1du2 = 0, df1du3 = 0;
	double df1du4 = 0, df1du5 = 0, df1du6 = 0;
	double df1dm1 = 1, df1dm2 = 0, df1dm3 = 0;
	double df1dm4 = 0, df1dm5 = 0, df1dm6 = 0;
	double df1dm7 = 1, df1dm8 = 0, df1dm9 = 0;
	double df1dm10= 0, df1dm11= 0, df1dm12= 0;
	//f2: -py_L + py_P + py_Q = 0
	double df2du1 = 0, df2du2 =-1, df2du3 = 0;
	double df2du4 = 0, df2du5 = 0, df2du6 = 0;
	double df2dm1 = 0, df2dm2 = 1, df2dm3 = 0;
	double df2dm4 = 0, df2dm5 = 0, df2dm6 = 0;
	double df2dm7 = 0, df2dm8 = 1, df2dm9 = 0;
	double df2dm10= 0, df2dm11= 0, df2dm12= 0;
	//f3: -pz_L + pz_P + pz_Q = 0
	double df3du1 = 0, df3du2 = 0, df3du3 =-1;
	double df3du4 = 0, df3du5 = 0, df3du6 = 0;
	double df3dm1 = 0, df3dm2 = 0, df3dm3 = 1;
	double df3dm4 = 0, df3dm5 = 0, df3dm6 = 0;
	double df3dm7 = 0, df3dm8 = 0, df3dm9 = 1;
	double df3dm10= 0, df3dm11= 0, df3dm12= 0;
	//f4: -E_L + E_P + E_Q = 0
	double df4du1 = -dE_Ldpx_L;
	double df4du2 = -dE_Ldpy_L;
	double df4du3 = -dE_Ldpz_L;
	double df4du4 = 0, df4du5 = 0, df4du6 = 0;
	double df4dm1 = dE_Pdpx_P;
	double df4dm2 = dE_Pdpy_P;
	double df4dm3 = dE_Pdpz_P;
	double df4dm4 = 0, df4dm5 = 0, df4dm6 = 0;
	double df4dm7 = dE_Qdpx_Q;
	double df4dm8 = dE_Qdpy_Q;
	double df4dm9 = dE_Qdpz_Q;
	double df4dm10 = 0, df4dm11 = 0, df4dm12 = 0;
	//f5 - f7: Vertex Constraints for P Q//
	//f5: vx_P - vx_Q = 0
	double df5du1 = 0, df5du2 = 0, df5du3 = 0;
	double df5du4 = 1, df5du5 = 0, df5du6 = 0;
	double df5dm1 = 0, df5dm2 = 0, df5dm3 = 0;
	double df5dm4 = 1, df5dm5 = 0, df5dm6 = 0;
	double df5dm7 = 0, df5dm8 = 0, df5dm9 = 0;
	double df5dm10=-1, df5dm11 = 0, df5dm12 = 0;
	//f6: vy_P - vy_Q = 0
	double df6du1 = 0, df6du2 = 0, df6du3 = 0;
	double df6du4 = 0, df6du5 = 0, df6du6 = 0;
	double df6dm1 = 0, df6dm2 = 0, df6dm3 = 0;
	double df6dm4 = 0, df6dm5 = 1, df6dm6 = 0;
	double df6dm7 = 0, df6dm8 = 0, df6dm9 = 0;
	double df6dm10= 0, df6dm11=-1, df6dm12= 0;
	//f7: vz_P - vz_Q = 0
	double df7du1 = 0, df7du2 = 0, df7du3 = 0;
	double df7du4 = 0, df7du5 = 0, df7du6 = 0;
	double df7dm1 = 0, df7dm2 = 0, df7dm3 = 0;
	double df7dm4 = 0, df7dm5 = 0, df7dm6 = 1;
	double df7dm7 = 0, df7dm8 = 0, df7dm9 = 0;
	double df7dm10= 0, df7dm11= 0, df7dm12=-1;
	//f8 -f10: L vertex determination.//
	//f8: -vx_L + 0.5 *(vx_P + vx_Q) = 0
	double df8du1 = 0, df8du2 = 0, df8du3 = 0;
	double df8du4 =-1, df8du5 = 0, df8du6 = 0;
	double df8dm1 = 0, df8dm2 = 0, df8dm3 = 0;
	double df8dm4 =0.5,df8dm5 = 0, df8dm6 = 0;
	double df8dm7 = 0, df8dm8 = 0, df8dm9 = 0;
	double df8dm10=0.5,df8dm11= 0, df8dm12= 0;
	//f9: -vy_L + 0.5 *(vy_P + vy_Q) = 0
	double df9du1 = 0, df9du2 = 0, df9du3 = 0;
	double df9du4 = 0, df9du5 =-1, df9du6 = 0;
	double df9dm1 = 0, df9dm2 = 0, df9dm3 = 0;
	double df9dm4 = 0, df9dm5 =0.5,df9dm6 = 0;
	double df9dm7 = 0, df9dm8 = 0, df9dm9 = 0;
	double df9dm10= 0, df9dm11=0.5,df9dm12= 0;
	//f10:-vz_L + 0.5 *(vz_P + vz_Q) = 0
	double df10du1 = 0, df10du2 = 0, df10du3 = 0;
	double df10du4 = 0, df10du5 = 0, df10du6 =-1;
	double df10dm1 = 0, df10dm2 = 0, df10dm3 = 0;
	double df10dm4 = 0, df10dm5 = 0, df10dm6 =0.5;
	double df10dm7 = 0, df10dm8 = 0, df10dm9 = 0;
	double df10dm10= 0, df10dm11= 0, df10dm12=0.5;

	double fs[10] = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10};
	double dfdms[400] ;
	double dfdus[400] ;
	if(UseVertexFlag){
	}
	else{
		double temp[] = {// 10 constraints, 12 measurement params. =  120 elements
			df1dm1, df1dm2, df1dm3, df1dm4	,df1dm5, df1dm6, df1dm7, df1dm8, df1dm9, df1dm10, df1dm11, df1dm12, 
			df2dm1, df2dm2, df2dm3, df2dm4	,df2dm5, df2dm6, df2dm7, df2dm8, df2dm9, df2dm10, df2dm11, df2dm12, 
			df3dm1, df3dm2, df3dm3, df3dm4	,df3dm5, df3dm6, df3dm7, df3dm8, df3dm9, df3dm10, df3dm11, df3dm12, 
			df4dm1, df4dm2, df4dm3, df4dm4	,df4dm5, df4dm6, df4dm7, df4dm8, df4dm9, df4dm10, df4dm11, df4dm12, 
			df5dm1, df5dm2, df5dm3, df5dm4	,df5dm5, df5dm6, df5dm7, df5dm8, df5dm9, df5dm10, df5dm11, df5dm12, 
			df6dm1, df6dm2, df6dm3, df6dm4	,df6dm5, df6dm6, df6dm7, df6dm8, df6dm9, df6dm10, df6dm11, df6dm12, 
			df7dm1, df7dm2, df7dm3, df7dm4	,df7dm5, df7dm6, df7dm7, df7dm8, df7dm9, df7dm10, df7dm11, df7dm12, 
			df8dm1, df8dm2, df8dm3, df8dm4	,df8dm5, df8dm6, df8dm7, df8dm8, df8dm9, df8dm10, df8dm11, df8dm12, 
			df9dm1, df9dm2, df9dm3, df9dm4	,df9dm5, df9dm6, df9dm7, df9dm8, df9dm9, df9dm10, df9dm11, df9dm12, 
			df10dm1,df10dm2,df10dm3,df10dm4	,df10dm5,df10dm6,df10dm7,df10dm8,df10dm9,df10dm10,df10dm11,df10dm12
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = { // 10 constraints, 6 unknown params. = 60 elements
			df1du1, df1du2, df1du3, df1du4, df1du5, df1du6,
			df2du1, df2du2, df2du3, df2du4, df2du5, df2du6,
			df3du1, df3du2, df3du3, df3du4, df3du5, df3du6,
			df4du1, df4du2, df4du3, df4du4, df4du5, df4du6,
			df5du1, df5du2, df5du3, df5du4, df5du5, df5du6,
			df6du1, df6du2, df6du3, df6du4, df6du5, df6du6,
			df7du1, df7du2, df7du3, df7du4, df7du5, df7du6,
			df8du1, df8du2, df8du3, df8du4, df8du5, df8du6,
			df9du1, df9du2, df9du3, df9du4, df9du5, df9du6,
			df10du1,df10du2,df10du3,df10du4,df10du5,df10du6
		};
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
	}

	TMatrixD FMat(nConst,1,fs);
	TMatrixD dFdM(nConst,nMeas,dfdms);
	TMatrixD dFdU(nConst,nUnkn,dfdus);
	FMats.push_back(FMat);//Constraint Matrices for Each step
	dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
	dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void MassVertexFitter::SampleStepPoint(int steps){
#if DebugCKF
	std::cout<<"MassVertexFitter::SampleStepPoint() step = "<<steps<<std::endl;
#endif
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double px_P,py_P,pz_P,px_Q,py_Q,pz_Q,px_L,py_L,pz_L;
	double vx_P,vy_P,vz_P,vx_Q,vy_Q,vz_Q,vx_L,vy_L,vz_L;
	if(UseVertexFlag){
	}
	else{
		px_P = Meas(0,0); 
		py_P = Meas(1,0); 
		pz_P = Meas(2,0);
		vx_P = Meas(3,0);
		vy_P = Meas(4,0);
		vz_P = Meas(5,0);

		px_Q = Meas(6,0); 
		py_Q = Meas(7,0); 
		pz_Q = Meas(8,0); 
		vx_Q = Meas(9,0);
		vy_Q = Meas(10,0);
		vz_Q = Meas(11,0);

		px_L = Unkn(0,0);
		py_L = Unkn(1,0);
		pz_L = Unkn(2,0);
		vx_L = Unkn(3,0);
		vy_L = Unkn(4,0);
		vz_L = Unkn(5,0);
	}
	TVector3 TV_P(px_P,py_P,pz_P);
	TVector3 TV_Q(px_Q,py_Q,pz_Q);
	TVector3 TV_L(px_L,py_L,pz_L);
	double p_P = TV_P.Mag();double p_Q = TV_Q.Mag();double p_L = TV_L.Mag();

	VPCor = TVector3(vx_P,vy_P,vz_P);
	VQCor = TVector3(vx_Q,vy_Q,vz_Q);
	VLCor = TVector3(vx_L,vy_L,vz_L);
	vector<double> MV = {
	px_P,py_P,pz_P,vx_P,vy_P,vz_P,
	px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q
	};
	vector<double> UV = {
	px_L,py_L,pz_L,
	vx_L,vy_L,vz_L
	};

	TLorentzVector PP(px_P,py_P,pz_P,hypot(mP,p_P));
	TLorentzVector QQ(px_Q,py_Q,pz_Q,hypot(mQ,p_Q));
	TLorentzVector LL(px_L,py_L,pz_L,hypot(mL,p_L));
	auto L_PQ = PP+QQ;

	double MassDiffL = L_PQ.Mag()-mL;
	PCor = PP;
	QCor = QQ;
	LCor = LL;
	MassDiffsL.push_back(MassDiffL);
}
TMatrixD
MassVertexFitter::JacobianSphToCart(double p, double th, double ph){//Legacy
	// x = p sin(th) cos(ph)
	// y = p sin(th) sin(ph)
	// z = p cos(th)
	//V_c = J^T V J|->
	//		dxdp, dxdth,dxdph
	//J	=	dydp, dydth,dydph
	//		dzdp, dzdth,dzdph

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
MassVertexFitter::CalcVariance(int istep){//Legacy
	return;
}








void
MassVertexFitter::Rotate(){//Legacy
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
MassVertexFitter::ToDecayPlane(){//Legacy
	auto Zaxis =(P + Q).Vect();
	auto vP = P.Vect();
	auto vQ = Q.Vect();
	auto Yaxis = vP.Cross(vQ);
//	double YNorm = 1./(Yaxis.Mag());
//	Yaxis = YNorm * Yaxis;
	double Th_F = Zaxis.Theta();
	double Ph_F = Zaxis.Phi();
	double RotZ[9] ={
		cos(-Ph_F),	-sin(-Ph_F),	0,	
		sin(-Ph_F),	cos(-Ph_F),		0,
		0,					0,						1
	};
	double RotY[9] ={
		cos(Th_F),	0,					-sin(Th_F),
		0,					-1,					0,
		sin(Th_F),	0,					cos(Th_F)
	};
	TMatrixD RZ(3,3,RotZ);
	TMatrixD RY(3,3,RotY);
	TMatrixD R_F = RY * RZ;
	Yaxis = R_F * Yaxis;
	double Th_Y = Yaxis.Theta();
	double Ph_Y = Yaxis.Phi();
	double RotX[9] ={
		1,				0,					0,
		0,				cos(-Ph_Y),	-sin(-Ph_Y),
		0,				sin(-Ph_Y),	cos(-Ph_Y)
	};
	TMatrixD RX(3,3,RotX);
	Yaxis = RX * Yaxis;

}

double
MassVertexFitter::CalcLambdaExtrapolationParameter(vector<double> Meas_,vector<double> Unkn_){
	/*
	We need to back-propagate Lambda, and evaluate the Xi decay vertex
	 as the midpoint of the Lambda and pi2(R, in this code).
	 The Lambda passes through the Lambda decay vertex VL, and the direction is cos(th_L), phi_L.
	 The Lambda trajectory is expresses as:
	 L = VL + t*DL
	 where DL = (sin(th_L)*cos(phi_L), sin(th_L)*sin(phi_L), cos(th_L)) the direction of Lambda.
	 For the closest point of L to VR, the following condition should be satisfied:
	 (VR - VL - t*DL) . DL = 0
	 since the direction from VR to the closest point should be vertical to the Lambda trajectory. 
	 Hence, t = (VR - VL) . DL / (DL . DL) = (VR - VL) . DL, since DL is a unit vector
	*/
	double px_P  = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double vx_P = Meas_[3];
	double vy_P = Meas_[4];
	double vz_P = Meas_[5];
	double px_Q  = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];
	double vx_Q = Meas_[9];
	double vy_Q = Meas_[10];
	double vz_Q = Meas_[11];
	double vx_R_= Meas_[15];
	double vy_R_= Meas_[16];
	double vz_R_= Meas_[17];

	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();
	TVector3 D_L(u_L, v_L, w_L);

	TVector3 V_P(vx_P,vy_P,vz_P);
	TVector3 V_Q(vx_Q,vy_Q,vz_Q);
	TVector3 V_L = 0.5*(V_P + V_Q);
	
	TVector3 VR(vx_R_,vy_R_,vz_R_);
	return (VR - V_L).Dot(D_L);
}
TVector3
MassVertexFitter::CalcLambdaDirectionalDerivativesM(int idx, vector<double> Meas_){
	double px_P = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double px_Q = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];
	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double p_L = TV_L.Mag();
	double u_L = TV_L.x()/p_L;
	double v_L = TV_L.y()/p_L;
	double w_L = TV_L.z()/p_L;
	TVector3 D_L(u_L, v_L, w_L);
	
	double du_Ldpx_P = (1 - u_L*u_L)/p_L;
	double dv_Ldpx_P = -v_L*u_L/p_L;
	double dw_Ldpx_P = -w_L*u_L/p_L;
	
	double du_Ldpy_P = -u_L*v_L/p_L;
	double dv_Ldpy_P = (1 - v_L*v_L)/p_L;
	double dw_Ldpy_P = -w_L*v_L/p_L;
	
	double du_Ldpz_P = -u_L*w_L/p_L;
	double dv_Ldpz_P = -v_L*w_L/p_L;
	double dw_Ldpz_P = (1 - w_L*w_L)/p_L;
	
	double du_Ldpx_Q = (1 - u_L*u_L)/p_L;
	double dv_Ldpx_Q = -v_L*u_L/p_L;
	double dw_Ldpx_Q = -w_L*u_L/p_L;
	
	double du_Ldpy_Q = -u_L*v_L/p_L;
	double dv_Ldpy_Q = (1 - v_L*v_L)/p_L;
	double dw_Ldpy_Q = -w_L*v_L/p_L;
	
	double du_Ldpz_Q = -u_L*w_L/p_L;
	double dv_Ldpz_Q = -v_L*w_L/p_L;
	double dw_Ldpz_Q = (1 - w_L*w_L)/p_L;

	TVector3 dD_Ldpx_P(du_Ldpx_P, dv_Ldpx_P, dw_Ldpx_P);
	TVector3 dD_Ldpy_P(du_Ldpy_P, dv_Ldpy_P, dw_Ldpy_P);
	TVector3 dD_Ldpz_P(du_Ldpz_P, dv_Ldpz_P, dw_Ldpz_P);
	TVector3 dD_Ldpx_Q(du_Ldpx_Q, dv_Ldpx_Q, dw_Ldpx_Q);
	TVector3 dD_Ldpy_Q(du_Ldpy_Q, dv_Ldpy_Q, dw_Ldpy_Q);
	TVector3 dD_Ldpz_Q(du_Ldpz_Q, dv_Ldpz_Q, dw_Ldpz_Q);

	switch(idx){
		case 0:
			return dD_Ldpx_P;
		case 1:
			return dD_Ldpy_P;
		case 2:
			return dD_Ldpz_P;
		case 6:
			return dD_Ldpx_Q;
		case 7:
			return dD_Ldpy_Q;
		case 8:
			return dD_Ldpz_Q;
		default:
			return TVector3(0,0,0);
	}
}
double 
MassVertexFitter::CalcLambdaExtrapolationParameterDerivativesM(int idx, vector<double> Meas_, vector<double> Unkn_){
	//A function to calculate the derivatives of the Lambda extrapolation parameter t by the measurement parameters. 
	double px_P = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double vx_P = Meas_[3];
	double vy_P = Meas_[4];
	double vz_P = Meas_[5];
	double px_Q = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];
	double vx_Q = Meas_[9];
	double vy_Q = Meas_[10];
	double vz_Q = Meas_[11];
	double vx_R = Meas_[15];
	double vy_R = Meas_[16];
	double vz_R = Meas_[17];

	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();
	double p_L = TV_L.Mag();
	TVector3 D_L(u_L, v_L, w_L);
	

	TVector3 V_P(vx_P,vy_P,vz_P);
	TVector3 V_Q(vx_Q,vy_Q,vz_Q);
	TVector3 V_L = 0.5*(V_P + V_Q);
	TVector3 V_R(vx_R,vy_R,vz_R);

	double extrap = (V_R - V_L).Dot(D_L);

	double derivative = 0;
	switch(idx){
		case 0:
			derivative = (V_R -V_L).Dot(Calc_dDdM(0,Meas_));
			break;
		case 1:
			derivative = (V_R - V_L).Dot(Calc_dDdM(1,Meas_));
			break;
		case 2:
			derivative = (V_R -V_L).Dot(Calc_dDdM(2,Meas_));
			break;
		case 3:
			derivative = -0.5*D_L.x();
			break;
		case 4:
			derivative = -0.5*D_L.y();
			break;
		case 5:
			derivative = -0.5*D_L.z();
			break;
		case 6:
			derivative = (V_R - V_L).Dot(Calc_dDdM(6,Meas_));
			break;
		case 7:
			derivative = (V_R - V_L).Dot(Calc_dDdM(7,Meas_));
			break;
		case 8:
			derivative = (V_R - V_L).Dot(Calc_dDdM(8,Meas_));
			break;
		case 9:
			derivative = -0.5*D_L.x();
			break;
		case 10:
			derivative = -0.5*D_L.y();
			break;
		case 11:
			derivative = -0.5*D_L.z();
			break;
		case 15:
			derivative = D_L.x();
			break;
		case 16:
			derivative = D_L.y();
			break;
		case 17:
			derivative = D_L.z();
			break;
		default:
			derivative = 0;
			break;
	}
	return derivative;
}
double
MassVertexFitter::CalcLambdaExtrapolationParameterDerivativesU(int idx, vector<double> Meas_, vector<double> Unkn_){
	//A function to calculate the derivatives of the Lambda extrapolation parameter t by the Unkn parameters.
	//Not used here, since L vtx are not unkn parameters in this code. 
	return 0;
}
TVector3
MassVertexFitter::CalcV_LX(vector<double> Meas_){
	double px_P  = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double px_Q  = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];
	double vx_P = Meas_[3];
	double vy_P = Meas_[4];
	double vz_P = Meas_[5];
	double vx_Q = Meas_[9];
	double vy_Q = Meas_[10];
	double vz_Q = Meas_[11];
	
	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();

	double t_cor_ = CalcLambdaExtrapolationParameter(Meas_, {});//Unmeas are left blank.
	double vx_LX = 0.5*(vx_P + vx_Q) + t_cor_ * u_L; 
	double vy_LX = 0.5*(vy_P + vy_Q) + t_cor_ * v_L; 
	double vz_LX = 0.5*(vz_P + vz_Q) + t_cor_ * w_L; 

	return TVector3(vx_LX, vy_LX, vz_LX);
}
TVector3
MassVertexFitter::CalcV_LXDerivativesM(int idx, vector<double> Meas_){
	double px_P  = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double px_Q  = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];
	double vx_P = Meas_[3];
	double vy_P = Meas_[4];
	double vz_P = Meas_[5];
	double vx_Q = Meas_[9];
	double vy_Q = Meas_[10];
	double vz_Q = Meas_[11];

	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();
	double t_cor_ = CalcLambdaExtrapolationParameter(Meas_, {});//Unmeas are left blank.
	TVector3 D_L(u_L, v_L, w_L);
	double vx_LX = 0.5*(vx_P + vx_Q) + t_cor_ * u_L; 
	double vy_LX = 0.5*(vy_P + vy_Q) + t_cor_ * v_L; 
	double vz_LX = 0.5*(vz_P + vz_Q) + t_cor_ * w_L; 

	double dvx_LXdM, dvy_LXdM, dvz_LXdM;
	switch(idx){
		case 0:
			dvx_LXdM = u_L * (Calc_dtdM(0,Meas_,{}))
			 + t_cor_ * Calc_dDdM(0,Meas_).x();
			dvy_LXdM = v_L * (Calc_dtdM(0,Meas_,{}))
			 + t_cor_ * Calc_dDdM(0,Meas_).y();
			dvz_LXdM = w_L * (Calc_dtdM(0,Meas_,{}))
			 + t_cor_ * Calc_dDdM(0,Meas_).z();
			break;
		case 1:
			dvx_LXdM = u_L * (Calc_dtdM(1,Meas_,{}))
			 + t_cor_ * Calc_dDdM(1,Meas_).x();
			dvy_LXdM = v_L * (Calc_dtdM(1,Meas_,{}))
			 + t_cor_ * Calc_dDdM(1,Meas_).y();
			dvz_LXdM = w_L * (Calc_dtdM(1,Meas_,{}))
			 + t_cor_ * Calc_dDdM(1,Meas_).z();
			break;
		case 2:
			dvx_LXdM = u_L * (Calc_dtdM(2,Meas_,{}))
			 + t_cor_ * Calc_dDdM(2,Meas_).x();
			dvy_LXdM = v_L * (Calc_dtdM(2,Meas_,{}))
			 + t_cor_ * Calc_dDdM(2,Meas_).y();
			dvz_LXdM = w_L * (Calc_dtdM(2,Meas_,{}))
			 + t_cor_ * Calc_dDdM(2,Meas_).z();
			break;
		case 3:
			dvx_LXdM = 0.5 + u_L * (Calc_dtdM(3,Meas_,{}));
			dvy_LXdM = v_L * (Calc_dtdM(3,Meas_,{}));
			dvz_LXdM = w_L * (Calc_dtdM(3,Meas_,{}));
			break;
		case 4:
			dvx_LXdM = u_L * (Calc_dtdM(4,Meas_,{}));
			dvy_LXdM = 0.5 + v_L * (Calc_dtdM(4,Meas_,{}));
			dvz_LXdM = w_L * (Calc_dtdM(4,Meas_,{}));
			break;
		case 5:
			dvx_LXdM = u_L * (Calc_dtdM(5,Meas_,{}));
			dvy_LXdM = v_L * (Calc_dtdM(5,Meas_,{}));
			dvz_LXdM = 0.5 + w_L * (Calc_dtdM(5,Meas_,{}));
			break;
		case 6:
			dvx_LXdM = u_L * (Calc_dtdM(6,Meas_,{}))
			 + t_cor_ * Calc_dDdM(6,Meas_).x();
			dvy_LXdM = v_L * (Calc_dtdM(6,Meas_,{}))
			 + t_cor_ * Calc_dDdM(6,Meas_).y();
			dvz_LXdM = w_L * (Calc_dtdM(6,Meas_,{}))
			 + t_cor_ * Calc_dDdM(6,Meas_).z();
			break;
		case 7:
			dvx_LXdM = u_L * (Calc_dtdM(7,Meas_,{}))
			 + t_cor_ * Calc_dDdM(7,Meas_).x();
			dvy_LXdM = v_L * (Calc_dtdM(7,Meas_,{}))
			 + t_cor_ * Calc_dDdM(7,Meas_).y();
			dvz_LXdM = w_L * (Calc_dtdM(7,Meas_,{}))
			 + t_cor_ * Calc_dDdM(7,Meas_).z();
			break;
		case 8:
			dvx_LXdM = u_L * (Calc_dtdM(8,Meas_,{}))
			 + t_cor_ * Calc_dDdM(8,Meas_).x();
			dvy_LXdM = v_L * (Calc_dtdM(8,Meas_,{}))
			 + t_cor_ * Calc_dDdM(8,Meas_).y();
			dvz_LXdM = w_L * (Calc_dtdM(8,Meas_,{}))
			 + t_cor_ * Calc_dDdM(8,Meas_).z();
			break;
		case 9:
			dvx_LXdM = 0.5 + u_L * (Calc_dtdM(9,Meas_,{}));
			dvy_LXdM = v_L * (Calc_dtdM(9,Meas_,{}));
			dvz_LXdM = w_L * (Calc_dtdM(9,Meas_,{}));
			break;
		case 10:
			dvx_LXdM = u_L * (Calc_dtdM(10,Meas_,{}));
			dvy_LXdM = 0.5 + v_L * (Calc_dtdM(10,Meas_,{}));
			dvz_LXdM = w_L * (Calc_dtdM(10,Meas_,{}));
			break;
		case 11:
			dvx_LXdM = u_L * (Calc_dtdM(11,Meas_,{}));
			dvy_LXdM = v_L * (Calc_dtdM(11,Meas_,{}));
			dvz_LXdM = 0.5 + w_L * (Calc_dtdM(11,Meas_,{}));
			break;
		case 15:
			dvx_LXdM = u_L * (Calc_dtdM(15,Meas_,{}));
			dvy_LXdM = v_L * (Calc_dtdM(15,Meas_,{}));
			dvz_LXdM = w_L * (Calc_dtdM(15,Meas_,{}));
			break;
		case 16:
			dvx_LXdM = u_L * (Calc_dtdM(16,Meas_,{}));
			dvy_LXdM = v_L * (Calc_dtdM(16,Meas_,{}));
			dvz_LXdM = w_L * (Calc_dtdM(16,Meas_,{}));
			break;
		case 17:
			dvx_LXdM = u_L * (Calc_dtdM(17,Meas_,{}));
			dvy_LXdM = v_L * (Calc_dtdM(17,Meas_,{}));
			dvz_LXdM = w_L * (Calc_dtdM(17,Meas_,{}));
			break;
		default:
			dvx_LXdM = 0;
			dvy_LXdM = 0;
			dvz_LXdM = 0;
			break;
	}
	return TVector3(dvx_LXdM, dvy_LXdM, dvz_LXdM);
}

TVector3
MassVertexFitter::CalcLambdaE1Vector(vector<double> Meas_){
	double px_P  = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double px_Q  = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];

	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();
	double norm = hypot(v_L,u_L);
	TVector3 E1(-v_L/norm, u_L/norm, 0);
	return E1;
}
TVector3
MassVertexFitter::CalcLambdaE1VectorDerivativesM(int idx, vector<double> Meas_){
	double px_P  = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double px_Q  = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];

	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();
	//E1 = z X D_L / |z X D_L|, where z = (0,0,1)
	// norm = 1./|z X D_L| = 1./hypot(v_L,u_L)
	double norm = 1./hypot(v_L,u_L);
	TVector3 E1_dir = TVector3(-v_L, u_L, 0);
	TVector3 E1 = norm * TVector3(-v_L, u_L, 0);
	double dnormdM = -(v_L*(Calc_dDdM(idx,Meas_).y()) + u_L*(Calc_dDdM(idx,Meas_).x()))*pow(norm,3);
	TVector3 Derivatives = dnormdM * E1_dir
	 + norm * TVector3(0,0,1).Cross(Calc_dDdM(idx,Meas_));
	return Derivatives;
}
TVector3
MassVertexFitter::CalcLambdaE2Vector(vector<double> Meas_){
	double px_P  = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double px_Q  = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];
	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();
	TVector3 D_L(u_L, v_L, w_L);

	TVector3 E1 = CalcLambdaE1Vector(Meas_);
	TVector3 E2 = D_L.Cross(E1);
	return E2;
}
TVector3
MassVertexFitter::CalcLambdaE2VectorDerivativesM(int idx, vector<double> Meas_){
	double px_P  = Meas_[0];
	double py_P = Meas_[1];
	double pz_P = Meas_[2];
	double px_Q  = Meas_[6];
	double py_Q = Meas_[7];
	double pz_Q = Meas_[8];
	TVector3 TV_P(px_P, py_P, pz_P);
	TVector3 TV_Q(px_Q, py_Q, pz_Q);
	TVector3 TV_L = TV_P + TV_Q;
	double u_L = TV_L.x()/TV_L.Mag();
	double v_L = TV_L.y()/TV_L.Mag();
	double w_L = TV_L.z()/TV_L.Mag();
	TVector3 D_L(u_L, v_L, w_L);
	TVector3 E1 = CalcLambdaE1Vector(Meas_);

	TVector3 dD_LdM = Calc_dDdM(idx,Meas_);
	TVector3 dE1dM = CalcLambdaE1VectorDerivativesM(idx,Meas_);
	TVector3 Derivatives = dD_LdM.Cross(E1) + D_L.Cross(dE1dM);
	return Derivatives;
}
#endif
