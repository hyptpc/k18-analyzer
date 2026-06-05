#include "CascadeFitter2.hh"
#ifndef CascadeFitter2_cc
#define CascadeFitter2_cc
#define DebugCKF 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void CascadeFitter2::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	UseVertexFlag = status;
	Clear();
	if(UseVertexFlag){//Production and Decay vertex of R is known, hence the direction is known. However, the momentum magnitude is not directly measured. I used the sum of P and Q as an initial value.
	}
	Initialize();
}
void CascadeFitter2::Initialize(){
	Initialized = 1;
#if DebugCKF
	cout<<"Initializing..."<<endl;
#endif
	Clear();
	if(UseVertexFlag){
		nMeas = 11;nUnkn = 2; nConst = 9; 
	}
	else{
		nMeas = 9;nUnkn = 3; nConst = 5;
	}
	mP = P.Mag();
	TVector3 TV_P = P.Vect();
	double px_P = TV_P.x();	
	double py_P = TV_P.y();
	double pz_P = TV_P.z();
	
	mQ = Q.Mag();
	TVector3 TV_Q = Q.Vect(); 
	double px_Q = TV_Q.x();
	double py_Q = TV_Q.y();
	double pz_Q = TV_Q.z();

	mR = R.Mag();
	TVector3 TV_R = R.Vect(); 
	double px_R = TV_R.x();
	double py_R = TV_R.y();
	double pz_R = TV_R.z();


	TVector3 TV_X = (P+Q+R).Vect();
	double px_X = TV_X.x();
	double py_X = TV_X.y();
	double pz_X = TV_X.z();
	double meas[20];
	double unkn[20];
	if(UseVertexFlag){
	}
	else{
		double temp[] = {px_P,py_P,pz_P,px_Q,py_Q,pz_Q,px_R,py_R,pz_R};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {px_X,py_X,pz_X};
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
	MassDiffsX.push_back(1e9);
	std::cout<<"CascadeFitter2::Initialize() done"<<std::endl;
}
void CascadeFitter2::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double px_R,py_R,pz_R,px_P,py_P,pz_P,px_Q,py_Q,pz_Q,px_X,py_X,pz_X;	
	if(UseVertexFlag){
	}
	else{
		px_P= Meas(0,0); 
		py_P= Meas(1,0); 
		pz_P= Meas(2,0);

		px_Q= Meas(3,0); 
		py_Q= Meas(4,0); 
		pz_Q= Meas(5,0); 

		px_R= Meas(6,0);
		py_R= Meas(7,0);
		pz_R= Meas(8,0);

		px_X= Unkn(0,0);
		py_X= Unkn(1,0);
		pz_X= Unkn(2,0);
	}
	TVector3 TV_P(px_P,py_P,pz_P);
	TVector3 TV_Q(px_Q,py_Q,pz_Q);
	TVector3 TV_R(px_R,py_R,pz_R);
	TVector3 TV_L = TV_P + TV_Q;
	TVector3 TV_X(px_X,py_X,pz_X);
	double p_P = TV_P.Mag();double p_Q = TV_Q.Mag();double p_R = TV_R.Mag();double p_L = TV_L.Mag();double p_X = TV_X.Mag();

	double E_P = hypot(p_P,mP);
	double dE_Pdpx_P = px_P/E_P;
	double dE_Pdpy_P = py_P/E_P;
	double dE_Pdpz_P = pz_P/E_P;
	double dEdp_P = p_P/E_P;
	
	double E_Q = hypot(p_Q,mQ);
	double dE_Qdpx_Q = px_Q/E_Q;
	double dE_Qdpy_Q = py_Q/E_Q;
	double dE_Qdpz_Q = pz_Q/E_Q;

	double E_R = hypot(p_R,mR);
	double dEdp_R = p_R/E_R;
	double dE_Rdpx_R = px_R/E_R;
	double dE_Rdpy_R = py_R/E_R;
	double dE_Rdpz_R = pz_R/E_R;

	double E_L = hypot(p_L,mL);
	double dE_Ldpx_P =(px_P + px_Q)/E_L;
	double dE_Ldpy_P =(py_P + py_Q)/E_L;
	double dE_Ldpz_P =(pz_P + pz_Q)/E_L;
	double dE_Ldpx_Q =(px_P + px_Q)/E_L;
	double dE_Ldpy_Q =(py_P + py_Q)/E_L;
	double dE_Ldpz_Q =(pz_P + pz_Q)/E_L;

	double E_X = hypot(p_X,mX);
	double dE_Xdpx_X = px_X/E_X;
	double dE_Xdpy_X = py_X/E_X;
	double dE_Xdpz_X = pz_X/E_X;

	double f1 = -px_X + px_P + px_Q + px_R;
	double f2 = -py_X + py_P + py_Q + py_R;
	double f3 = -pz_X + pz_P + pz_Q + pz_R;
	double f4 = -E_L + E_P + E_Q;//Constraint on Lambda Energy
	double f5 = -E_X + E_P +E_Q + E_R;//Constraint on Xi Energy
	//f1 - f5: Kinematic Constraints//
	//f1: -px_X + px_P + px_Q + px_R = 0
	double df1du1 =-1, df1du2 = 0, df1du3 = 0;
	double df1dm1 = 1, df1dm2 = 0, df1dm3 = 0;
	double df1dm4 = 1, df1dm5 = 0, df1dm6 = 0;
	double df1dm7 = 1, df1dm8 = 0, df1dm9 = 0;
	
	//f2: -py_X + py_P + py_Q + py_R = 0
	double df2du1 = 0, df2du2 =-1, df2du3 = 0;
	double df2dm1 = 0, df2dm2 = 1, df2dm3 = 0;
	double df2dm4 = 0, df2dm5 = 1, df2dm6 = 0;
	double df2dm7 = 0, df2dm8 = 1, df2dm9 = 0;
	//f3: -pz_X + pz_P + pz_Q + pz_R = 0
	double df3du1 = 0, df3du2 = 0, df3du3 =-1;
	double df3dm1 = 0, df3dm2 = 0, df3dm3 = 1;
	double df3dm4 = 0, df3dm5 = 0, df3dm6 = 1;
	double df3dm7 = 0, df3dm8 = 0, df3dm9 = 1;
	//f4: -E_L + E_P + E_Q = 0
	double df4du1 = 0, df4du2 = 0, df4du3 = 0;
	double df4dm1 = -dE_Ldpx_P + dE_Pdpx_P;
	double df4dm2 = -dE_Ldpy_P + dE_Pdpy_P;
	double df4dm3 = -dE_Ldpz_P + dE_Pdpz_P;
	double df4dm4 = -dE_Ldpx_Q + dE_Qdpx_Q;
	double df4dm5 = -dE_Ldpy_Q + dE_Qdpy_Q;
	double df4dm6 = -dE_Ldpz_Q + dE_Qdpz_Q;
	double df4dm7 = 0, df4dm8 = 0, df4dm9 = 0;
	//f5: -E_X + E_P + E_Q + E_R = 0
	double df5du1 = -dE_Xdpx_X;
	double df5du2 = -dE_Xdpy_X;
	double df5du3 = -dE_Xdpz_X;
	double df5dm1 = dE_Pdpx_P;
	double df5dm2 = dE_Pdpy_P;
	double df5dm3 = dE_Pdpz_P;
	double df5dm4 = dE_Qdpx_Q;
	double df5dm5 = dE_Qdpy_Q;
	double df5dm6 = dE_Qdpz_Q;
	double df5dm7 = dE_Rdpx_R;
	double df5dm8 = dE_Rdpy_R;
	double df5dm9 = dE_Rdpz_R;

	double fs[]={f1,f2,f3,f4,f5};
	double dfdms[200] ;
	double dfdus[200] ;
	if(UseVertexFlag){
	}
	else{
		double temp[] = {
			df1dm1,df1dm2,df1dm3,df1dm4	,df1dm5,df1dm6,df1dm7,df1dm8,df1dm9,
			df2dm1,df2dm2,df2dm3,df2dm4	,df2dm5,df2dm6,df2dm7,df2dm8,df2dm9,
			df3dm1,df3dm2,df3dm3,df3dm4	,df3dm5,df3dm6,df3dm7,df3dm8,df3dm9,
			df4dm1,df4dm2,df4dm3,df4dm4	,df4dm5,df4dm6,df4dm7,df4dm8,df4dm9,
			df5dm1,df5dm2,df5dm3,df5dm4	,df5dm5,df5dm6,df5dm7,df5dm8,df5dm9,
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = {
			df1du1,df1du2,df1du3,
			df2du1,df2du2,df2du3,
			df3du1,df3du2,df3du3,
			df4du1,df4du2,df4du3,
			df5du1,df5du2,df5du3
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
void CascadeFitter2::SampleStepPoint(int steps){
#if DebugCKF
	std::cout<<"CascadeFitter2::SampleStepPoint() step = "<<steps<<std::endl;
#endif
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double px_R,py_R,pz_R,px_P,py_P,pz_P,px_Q,py_Q,pz_Q,px_X,py_X,pz_X;	
	if(UseVertexFlag){
	}
	else{
		px_P= Meas(0,0); 
		py_P= Meas(1,0); 
		pz_P= Meas(2,0);

		px_Q= Meas(3,0); 
		py_Q= Meas(4,0); 
		pz_Q= Meas(5,0); 

		px_R= Meas(6,0);
		py_R= Meas(7,0);
		pz_R= Meas(8,0);

		px_X= Unkn(0,0);
		py_X= Unkn(1,0);
		pz_X= Unkn(2,0);
	}
	TVector3 TV_P(px_P,py_P,pz_P);
	TVector3 TV_Q(px_Q,py_Q,pz_Q);
	TVector3 TV_R(px_R,py_R,pz_R);
	TVector3 TV_L = TV_P + TV_Q;
	TVector3 TV_X(px_X,py_X,pz_X);
	double px_L = TV_L.X();double py_L = TV_L.Y();double pz_L = TV_L.Z();
	double p_P = TV_P.Mag();double p_Q = TV_Q.Mag();double p_R = TV_R.Mag();double p_L = TV_L.Mag();double p_X = TV_X.Mag();

	TLorentzVector PP(px_P,py_P,pz_P,hypot(mP,p_P));
	TLorentzVector QQ(px_Q,py_Q,pz_Q,hypot(mQ,p_Q));
	TLorentzVector RR(px_R,py_R,pz_R,hypot(mR,p_R));
	TLorentzVector LL(px_L,py_L,pz_L,hypot(mL,p_L));
	TLorentzVector XX(px_X,py_X,pz_X,hypot(mX,p_X));
	auto L_PQ = PP+QQ;
	auto X_PQR = PP+QQ+RR;

	double MassDiffL = L_PQ.Mag()-mL;
	double MassDiffX = X_PQR.Mag()-mX;
	PCor = PP;
	QCor = QQ;
	RCor = RR;
	LCor = LL;
	XCor = XX;
	MassDiffsL.push_back(MassDiffL);
	MassDiffsX.push_back(MassDiffX);
}
CascadeFitter2::CascadeFitter2(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_){
	//Does Kinematic fitting for R -> P+Q decay.
	P=P_;
	Q=Q_;
	R=R_;
	Initialize();
};
TMatrixD
CascadeFitter2::JacobianSphToCart(double p, double th, double ph){
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
CascadeFitter2::CalcVariance(int istep){//Legacy
	return;
	auto Meas = Measurements.at(istep); 
	auto Unkn = Unknowns.at(istep);
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_L,th_L,ph_L,p_X,th_X,ph_X;	
	if(UseVertexFlag){
	}
	else{
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);

		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 

		p_R= Meas(6,0);
		th_R= Meas(7,0);
		ph_R= Meas(8,0);

		p_L= Unkn(0,0);
		th_L= Unkn(1,0);
		ph_L= Unkn(2,0);

		p_X= Unkn(3,0);
		th_X= Unkn(4,0);
		ph_X= Unkn(5,0);
	}
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
	/*
	Jsc_P.Print();
	Jsc_Q.Print();
	Jsc_PQ.Print();
	cin.ignore();
	*/
	TMatrixD Jsc_PQ_T = TransposeMatrix(Jsc_PQ);
	TMatrixD VMat = Variancies.at(istep);
	TMatrixD dV = dVMats.at(istep);
	TMatrixD VMat_C = Jsc_PQ_T*(VMat-dV)*Jsc_PQ;
//	TMatrixD VMat_C = Jsc_PQ_T*(VMat)*Jsc_PQ;
	double El_contract[18]={
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
//	VarianciesU.push_back(UVMat);
}















void
CascadeFitter2::Rotate(){
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
CascadeFitter2::ToDecayPlane(){
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
#endif
