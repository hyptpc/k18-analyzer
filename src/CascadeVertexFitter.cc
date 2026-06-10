#include "CascadeVertexFitter.hh"
#ifndef CascadeVertexFitter_cc
#define CascadeVertexFitter_cc
#define DebugCKF 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void CascadeVertexFitter::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	UseVertexFlag = status;
	Clear();
	Initialize();
}
CascadeVertexFitter::CascadeVertexFitter(TLorentzVector P_,TVector3 VP_,
	TLorentzVector Q_, TVector3 VQ_,
	 TLorentzVector R_ ,TVector3 VR_){
	//Does simultaneous Kinematic fitting for L -> P + Q and X -> L + R decays.

	P=P_;
	Q=Q_;
	R=R_;
	VP=VP_;
	VQ=VQ_;
	VR=VR_;
	ScaleParams = 0;
	Initialize();
};
void CascadeVertexFitter::Initialize(){
	Initialized = 1;
#if DebugCKF
	cout<<"Initializing..."<<endl;
#endif
	Clear();
	if(UseVertexFlag){
		nMeas = 11;nUnkn = 2; nConst = 9; 
	}
	else{
		nMeas = 18;nUnkn = 6; nConst = 13;
		//Meas : px,py,pz,vx,vy,vz for P,Q,R
		//Unkn : px,py,pz,vx,vy,vz for X
		//Const: px,py,pz conservation(3), E conservation(Mass constraint)(2) for L and X,
		// vertex constraint for L and X (6)
		//ndf = nConst - nUnkn = 7
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

	mR = R.Mag();
	TVector3 TV_R = R.Vect(); 
	double px_R = TV_R.x();
	double py_R = TV_R.y();
	double pz_R = TV_R.z();
	double vx_R = VR.x();
	double vy_R = VR.y();
	double vz_R = VR.z();

	TVector3 TV_L = (P+Q).Vect();
	double p_L = TV_L.Mag();
	double th_L = TV_L.Theta();
	double ph_L = TV_L.Phi();
	VL = (VP + VQ)*0.5;//Initial value of Lambda vertex is set to the averaged vertex of P and Q. This is not a constraint, but just an initial value. The fitter will move the vertex according to the constraints and the covariance matrix.
	double vx_L = VL.x();
	double vy_L = VL.y();
	double vz_L = VL.z();
	TVector3 D_L = TV_L.Unit();

	TVector3 TV_X = (P+Q+R).Vect();
	double px_X = TV_X.x();
	double py_X = TV_X.y();
	double pz_X = TV_X.z();
	double vx_X = 0;
	double vy_X = 0;
	double vz_X = 0;
	std::vector<double> MV = {
	px_P,py_P,pz_P,vx_P,vy_P,vz_P,
	px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q,
	px_R,py_R,pz_R,vx_R,vy_R,vz_R
	};
	std::vector<double> UV = {
	px_X,py_X,pz_X,
	vx_X,vy_X,vz_X//t should be evaluated before V_X derermination.
	};
	t = CalcLambdaExtrapolationParameter(MV,UV);
	VX = (VL + t*D_L + VR)*0.5;
	vx_X = VX.x();
	vy_X = VX.y();
	vz_X = VX.z();

	double meas[20];
	double unkn[20];
	if(UseVertexFlag){
	}
	else{
		double temp[] = {px_P,py_P,pz_P,vx_P,vy_P,vz_P,px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q,px_R,py_R,pz_R,vx_R,vy_R,vz_R};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {px_X,py_X,pz_X,vx_X,vy_X,vz_X};
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
	std::cout<<"CascadeVertexFitter::Initialize() done"<<std::endl;
}
void CascadeVertexFitter::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double px_P,py_P,pz_P,px_Q,py_Q,pz_Q,px_R,py_R,pz_R,px_X,py_X,pz_X;
	double vx_R,vy_R,vz_R,vx_P,vy_P,vz_P,vx_Q,vy_Q,vz_Q,vx_L,vy_L,vz_L,vx_X,vy_X,vz_X;
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

		px_R= Meas(12,0);
		py_R= Meas(13,0);
		pz_R= Meas(14,0);
		vx_R= Meas(15,0);
		vy_R= Meas(16,0);
		vz_R= Meas(17,0);

		px_X= Unkn(0,0);
		py_X= Unkn(1,0);
		pz_X= Unkn(2,0);
		vx_X= Unkn(3,0);
		vy_X= Unkn(4,0);
		vz_X= Unkn(5,0);
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

	TVector3 TV_R(px_R,py_R,pz_R);
	double p_R = TV_R.Mag();
	double E_R = hypot(p_R,mR);
	double dE_Rdpx_R = px_R/E_R;
	double dE_Rdpy_R = py_R/E_R;
	double dE_Rdpz_R = pz_R/E_R;

	double px_L = px_P + px_Q;
	double py_L = py_P + py_Q;
	double pz_L = pz_P + pz_Q;
	TVector3 TV_L(px_L,py_L,pz_L);
	double p_L = TV_L.Mag();
	double E_L = hypot(p_L,mL);
	double dEdp_L = p_L/E_L;
	
	double dE_Ldpx_P =(px_P + px_Q)/E_L;
	double dE_Ldpy_P =(py_P + py_Q)/E_L;
	double dE_Ldpz_P =(pz_P + pz_Q)/E_L;
	double dE_Ldpx_Q =(px_P + px_Q)/E_L;
	double dE_Ldpy_Q =(py_P + py_Q)/E_L;
	double dE_Ldpz_Q =(pz_P + pz_Q)/E_L;
	
	double u_L = TV_L.Unit().x(), v_L = TV_L.Unit().y(), w_L = TV_L.Unit().z();

	TVector3 TV_X(px_X,py_X,pz_X);
	double p_X = TV_X.Mag();
	double E_X = hypot(p_X,mX);
	double dE_Xdpx_X = px_X/E_X;
	double dE_Xdpy_X = py_X/E_X;
	double dE_Xdpz_X = pz_X/E_X;
	

	double f1 = -px_X + px_P + px_Q + px_R;
	double f2 = -py_X + py_P + py_Q + py_R;
	double f3 = -pz_X + pz_P + pz_Q + pz_R;
	double f4 = -E_L + E_P + E_Q;//Constraint on Lambda Energy
	double f5 = -E_X + E_P +E_Q + E_R;//Constraint on Xi Energy
	double f6 = vx_P - vx_Q;
	double f7 = vy_P - vy_Q;
	double f8 = vz_P - vz_Q;
	
	std::vector<double> MV = {
		px_P,py_P,pz_P,vx_P,vy_P,vz_P,
		px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q,
		px_R,py_R,pz_R,vx_R,vy_R,vz_R};
	std::vector<double> UV = {
		px_X,py_X,pz_X,
		vx_X,vy_X,vz_X};
	/*
	 Reconstructed Laambda should propagate back to the Xi decay vertex.
	 Note that this extrapolation parameter t_cor is also a function of the kinematic parameters,
	 so its derivatives should be considered in the calculation of the Jacobian matrix.
	*/
	t_cor = CalcLambdaExtrapolationParameter(MV,UV);
	TVector3 V_LX = CalcV_LX(MV);
	double f9 = -vx_X + 0.5 *(vx_R + V_LX.x()); 
	double f10= -vy_X + 0.5 *(vy_R + V_LX.y());
	double f11= -vz_X + 0.5 *(vz_R + V_LX.z());
	TVector3 V_R(vx_R,vy_R,vz_R);
	/*
		Note! Vertex constraint for Extrapolated Lambda and Pi2 should be 2-dimensional.
		From the definition of t_cor, (V_LX - V_R)*(D_L) = 0 is already constrained.
		The vertex constraint should only have two degree of freedom, so following
		cartesian representation, which has three degree of freedom, is invalid. 
		|f12 = vx_R - 0.5*(vx_P + vx_Q) - t_cor*u_L|
		|f13 = vy_R - 0.5*(vy_P + vy_Q) - t_cor*v_L|
		|f14 = vz_R - 0.5*(vz_P + vz_Q) - t_cor*w_L|

		Instead, two unit vectors, normal to D_L, should be defined.
		E_1 = TVector3(0,0,1).Cross(D_L).Unit();
		E_2 = D_L.Cross(E_1).Unit();
		Then, the constraint should be defined as:
		f12 := (V_R - V_LX)*E_1 = 0;
		f13 := (V_R - V_LX)*E_2 = 0;

		Here, the Z axis (0,0,1) is Y axis in E42 frame.
		No Lambda goes close to this direction, so this definition is numerically stable.
		Keep in mind that, two unit vectors, E_1 and E_2, are also functions of the kinematic parameters.
		Their derivatives should also be considered in the Jacobian matrix calculation.
	*/
	TVector3 E_1 = CalcLambdaE1Vector(MV); //(0,0,1) X (u,v,w)
	TVector3 E_2 = CalcLambdaE2Vector(MV); //E2 = (u,v,w) X E1
	double f12 = (V_R - V_LX)*E_1;
	double f13 = (V_R - V_LX)*E_2;
	//f1 - f5: Kinematic Constraints//
	//f1: -px_X + px_P + px_Q + px_R = 0
	double df1du1 =-1, df1du2 = 0, df1du3 = 0;
	double df1du4 = 0, df1du5 = 0, df1du6 = 0;
	double df1dm1 = 1, df1dm2 = 0, df1dm3 = 0;
	double df1dm4 = 0, df1dm5 = 0, df1dm6 = 0;
	double df1dm7 = 1, df1dm8 = 0, df1dm9 = 0;
	double df1dm10= 0, df1dm11= 0, df1dm12= 0;
	double df1dm13= 1, df1dm14= 0, df1dm15= 0;
	double df1dm16= 0, df1dm17= 0, df1dm18= 0;
	//f2: -py_X + py_P + py_Q + py_R = 0
	double df2du1 = 0, df2du2 =-1, df2du3 = 0;
	double df2du4 = 0, df2du5 = 0, df2du6 = 0;
	double df2dm1 = 0, df2dm2 = 1, df2dm3 = 0;
	double df2dm4 = 0, df2dm5 = 0, df2dm6 = 0;
	double df2dm7 = 0, df2dm8 = 1, df2dm9 = 0;
	double df2dm10= 0, df2dm11= 0, df2dm12= 0;
	double df2dm13= 0, df2dm14= 1, df2dm15= 0;
	double df2dm16= 0, df2dm17= 0, df2dm18= 0;
	//f3: -pz_X + pz_P + pz_Q + pz_R = 0
	double df3du1 = 0, df3du2 = 0, df3du3 =-1;
	double df3du4 = 0, df3du5 = 0, df3du6 = 0;
	double df3dm1 = 0, df3dm2 = 0, df3dm3 = 1;
	double df3dm4 = 0, df3dm5 = 0, df3dm6 = 0;
	double df3dm7 = 0, df3dm8 = 0, df3dm9 = 1;
	double df3dm10= 0, df3dm11= 0, df3dm12= 0;
	double df3dm13= 0, df3dm14= 0, df3dm15= 1;
	double df3dm16= 0, df3dm17= 0, df3dm18= 0;
	//f4: -E_L + E_P + E_Q = 0
	double df4du1 = 0, df4du2 = 0, df4du3 = 0;
	double df4du4 = 0, df4du5 = 0, df4du6 = 0;
	double df4dm1 = -dE_Ldpx_P + dE_Pdpx_P;
	double df4dm2 = -dE_Ldpy_P + dE_Pdpy_P;
	double df4dm3 = -dE_Ldpz_P + dE_Pdpz_P;
	double df4dm4 = 0, df4dm5 = 0, df4dm6 = 0;
	double df4dm7 = -dE_Ldpx_Q + dE_Qdpx_Q;
	double df4dm8 = -dE_Ldpy_Q + dE_Qdpy_Q;
	double df4dm9 = -dE_Ldpz_Q + dE_Qdpz_Q;
	double df4dm10= 0, df4dm11= 0, df4dm12= 0;
	double df4dm13= 0, df4dm14= 0, df4dm15= 0;
	double df4dm16= 0, df4dm17= 0, df4dm18= 0;
	//f5: -E_X + E_P + E_Q + E_R = 0
	double df5du1 = -dE_Xdpx_X;
	double df5du2 = -dE_Xdpy_X;
	double df5du3 = -dE_Xdpz_X;
	double df5du4 = 0, df5du5 = 0, df5du6 = 0;
	double df5dm1 = dE_Pdpx_P;
	double df5dm2 = dE_Pdpy_P;
	double df5dm3 = dE_Pdpz_P;
	double df5dm4 = 0, df5dm5 = 0, df5dm6 = 0;
	double df5dm7 = dE_Qdpx_Q;
	double df5dm8 = dE_Qdpy_Q;
	double df5dm9 = dE_Qdpz_Q;
	double df5dm10= 0, df5dm11= 0, df5dm12= 0;
	double df5dm13= dE_Rdpx_R;
	double df5dm14= dE_Rdpy_R;
	double df5dm15= dE_Rdpz_R;
	double df5dm16= 0, df5dm17= 0, df5dm18= 0;
	//f6 - f8: Vertex Constraints for p pi//
	//f6: vx_P - vx_Q = 0
	double df6du1 = 0, df6du2 = 0, df6du3 = 0;
	double df6du4 = 0, df6du5 = 0, df6du6 = 0;
	double df6dm1 = 0, df6dm2 = 0, df6dm3 = 0;
	double df6dm4 = 1, df6dm5 = 0, df6dm6 = 0;
	double df6dm7 = 0, df6dm8 = 0, df6dm9 = 0;
	double df6dm10=-1, df6dm11= 0, df6dm12= 0;
	double df6dm13= 0, df6dm14= 0, df6dm15= 0;
	double df6dm16= 0, df6dm17= 0, df6dm18= 0;
	//f7: vy_P - vy_Q = 0
	double df7du1 = 0, df7du2 = 0, df7du3 = 0;
	double df7du4 = 0, df7du5 = 0, df7du6 = 0;
	double df7dm1 = 0, df7dm2 = 0, df7dm3 = 0;
	double df7dm4 = 0, df7dm5 = 1, df7dm6 = 0;
	double df7dm7 = 0, df7dm8 = 0, df7dm9 = 0;
	double df7dm10= 0, df7dm11=-1, df7dm12= 0;
	double df7dm13= 0, df7dm14= 0, df7dm15= 0;
	double df7dm16= 0, df7dm17= 0, df7dm18= 0;
	//f8: vz_P - vz_Q = 0
	double df8du1 = 0, df8du2 = 0, df8du3 = 0;
	double df8du4 = 0, df8du5 = 0, df8du6 = 0;
	double df8dm1 = 0, df8dm2 = 0, df8dm3 = 0;
	double df8dm4 = 0, df8dm5 = 0, df8dm6 = 1;
	double df8dm7 = 0, df8dm8 = 0, df8dm9 = 0;
	double df8dm10= 0, df8dm11= 0, df8dm12=-1;
	double df8dm13= 0, df8dm14= 0, df8dm15= 0;
	double df8dm16= 0, df8dm17= 0, df8dm18= 0;
	//f9 -f11: Vertex Constraints for Xi decay//
	//f9: -vx_X + 0.5 *(vx_R + V_LX.x()) = 0
	double df9du1 = 0, df9du2 = 0, df9du3 = 0;
	double df9du4 =-1, df9du5 = 0, df9du6 = 0;//m1 = px_X,m4 = vx_P,m7 = px_Q,m10= vx_R...
	double df9dm1 = 0.5 * Calc_dV_LXdM(0,MV).x();
	double df9dm2 = 0.5 * Calc_dV_LXdM(1,MV).x();
	double df9dm3 = 0.5 * Calc_dV_LXdM(2,MV).x();
	double df9dm4 = 0.5 * Calc_dV_LXdM(3,MV).x();
	double df9dm5 = 0.5 * Calc_dV_LXdM(4,MV).x();
	double df9dm6 = 0.5 * Calc_dV_LXdM(5,MV).x();
	double df9dm7 = 0.5 * Calc_dV_LXdM(6,MV).x();
	double df9dm8 = 0.5 * Calc_dV_LXdM(7,MV).x();
	double df9dm9 = 0.5 * Calc_dV_LXdM(8,MV).x();
	double df9dm10= 0.5 * Calc_dV_LXdM(9,MV).x();
	double df9dm11= 0.5 * Calc_dV_LXdM(10,MV).x();
	double df9dm12= 0.5 * Calc_dV_LXdM(11,MV).x();
	double df9dm13= 0, df9dm14= 0, df9dm15= 0;
	double df9dm16= 0.5 *(1 + Calc_dV_LXdM(15,MV).x());
	double df9dm17= 0.5 * Calc_dV_LXdM(16,MV).x();
	double df9dm18= 0.5 * Calc_dV_LXdM(17,MV).x();
	//f10: -vy_X + 0.5 *(vy_R + V_LX.y()) = 0
	double df10du1 = 0, df10du2 = 0, df10du3 = 0;
	double df10du4 = 0, df10du5 =-1, df10du6 = 0;
	double df10dm1 = 0.5 * Calc_dV_LXdM(0,MV).y();
	double df10dm2 = 0.5 * Calc_dV_LXdM(1,MV).y();
	double df10dm3 = 0.5 * Calc_dV_LXdM(2,MV).y();
	double df10dm4 = 0.5 * Calc_dV_LXdM(3,MV).y();
	double df10dm5 = 0.5 * Calc_dV_LXdM(4,MV).y();
	double df10dm6 = 0.5 * Calc_dV_LXdM(5,MV).y();
	double df10dm7 = 0.5 * Calc_dV_LXdM(6,MV).y();
	double df10dm8 = 0.5 * Calc_dV_LXdM(7,MV).y();
	double df10dm9 = 0.5 * Calc_dV_LXdM(8,MV).y();
	double df10dm10= 0.5 * Calc_dV_LXdM(9,MV).y();
	double df10dm11= 0.5 * Calc_dV_LXdM(10,MV).y();
	double df10dm12= 0.5 * Calc_dV_LXdM(11,MV).y();
	double df10dm13= 0, df10dm14= 0, df10dm15= 0;
	double df10dm16= 0.5 * Calc_dV_LXdM(15,MV).y();
	double df10dm17= 0.5 * (1 + Calc_dV_LXdM(16,MV).y());
	double df10dm18= 0.5 * Calc_dV_LXdM(17,MV).y();
	//f11: -vz_X + 0.5 *(vz_R + V_LX.z()) = 0
	double df11du1 = 0, df11du2 = 0, df11du3 = 0;
	double df11du4 = 0, df11du5 = 0, df11du6 =-1;
	double df11dm1 = 0.5 * Calc_dV_LXdM(0,MV).z();
	double df11dm2 = 0.5 * Calc_dV_LXdM(1,MV).z();
	double df11dm3 = 0.5 * Calc_dV_LXdM(2,MV).z();
	double df11dm4 = 0.5 * Calc_dV_LXdM(3,MV).z();
	double df11dm5 = 0.5 * Calc_dV_LXdM(4,MV).z();
	double df11dm6 = 0.5 * Calc_dV_LXdM(5,MV).z();
	double df11dm7 = 0.5 * Calc_dV_LXdM(6,MV).z();
	double df11dm8 = 0.5 * Calc_dV_LXdM(7,MV).z();
	double df11dm9 = 0.5 * Calc_dV_LXdM(8,MV).z();
	double df11dm10= 0.5 * Calc_dV_LXdM(9,MV).z();
	double df11dm11= 0.5 * Calc_dV_LXdM(10,MV).z();
	double df11dm12= 0.5 * Calc_dV_LXdM(11,MV).z();
	double df11dm13= 0, df11dm14= 0, df11dm15= 0;
	double df11dm16= 0.5 * Calc_dV_LXdM(15,MV).z();
	double df11dm17= 0.5 * Calc_dV_LXdM(16,MV).z();
	double df11dm18= 0.5 * (1 + Calc_dV_LXdM(17,MV).z());
	//f12-f14 Vertex Constraints for L pi
	//f12: vx_R - 0.5*(vx_P + vx_Q) - t_cor*u_L = 0
	
	#if 0
	double du_Ldpx_P = Calc_dDdM(0, MV).x();
	double du_Ldpy_P = Calc_dDdM(1, MV).x();
	double du_Ldpz_P = Calc_dDdM(2, MV).x();
	double dv_Ldpx_P = Calc_dDdM(0, MV).y();
	double dv_Ldpy_P = Calc_dDdM(1, MV).y();
	double dv_Ldpz_P = Calc_dDdM(2, MV).y();
	double dw_Ldpx_P = Calc_dDdM(0, MV).z();
	double dw_Ldpy_P = Calc_dDdM(1, MV).z();
	double dw_Ldpz_P = Calc_dDdM(2, MV).z();
	double du_Ldpx_Q = Calc_dDdM(0, MV).x();
	double du_Ldpy_Q = Calc_dDdM(1, MV).x();
	double du_Ldpz_Q = Calc_dDdM(2, MV).x();
	double dv_Ldpx_Q = Calc_dDdM(0, MV).y();
	double dv_Ldpy_Q = Calc_dDdM(1, MV).y();
	double dv_Ldpz_Q = Calc_dDdM(2, MV).y();
	double dw_Ldpx_Q = Calc_dDdM(0, MV).z();
	double dw_Ldpy_Q = Calc_dDdM(1, MV).z();
	double dw_Ldpz_Q = Calc_dDdM(2, MV).z();
	double df12du1 = 0, df12du2 = 0, df12du3 = 0;
	double df12du4 = 0, df12du5 = 0, df12du6 = 0;
	double df12dm1 = -(t_cor * du_Ldpx_P + u_L * Calc_dtdM(0,MV,UV));
	double df12dm2 = -(t_cor * du_Ldpy_P + u_L * Calc_dtdM(1,MV,UV));
	double df12dm3 = -(t_cor * du_Ldpz_P + u_L * Calc_dtdM(2,MV,UV));
	double df12dm4 = -0.5 - u_L * Calc_dtdM(3,MV,UV);
	double df12dm5 = -u_L * Calc_dtdM(4,MV,UV);
	double df12dm6 = -u_L * Calc_dtdM(5,MV,UV);
	double df12dm7 = -(t_cor * du_Ldpx_Q + u_L * Calc_dtdM(6,MV,UV));
	double df12dm8 = -(t_cor * du_Ldpy_Q + u_L * Calc_dtdM(7,MV,UV));
	double df12dm9 = -(t_cor * du_Ldpz_Q + u_L * Calc_dtdM(8,MV,UV));
	double df12dm10= -0.5 - u_L * Calc_dtdM(9,MV,UV);
	double df12dm11= -u_L * Calc_dtdM(10,MV,UV);
	double df12dm12= -u_L * Calc_dtdM(11,MV,UV);
	double df12dm13= 0, df12dm14= 0, df12dm15= 0;
	double df12dm16= 1 - u_L * Calc_dtdM(15,MV,UV);
	double df12dm17= -u_L * Calc_dtdM(16,MV,UV);
	double df12dm18= -u_L * Calc_dtdM(17,MV,UV);
	//f13: vy_R - 0.5*(vy_P + vy_Q) - t_cor*v_L
	double df13du1 = 0, df13du2 = 0, df13du3 = 0;
	double df13du4 = 0, df13du5 = 0, df13du6 = 0;
	double df13dm1 = -(t_cor * dv_Ldpx_P + v_L * Calc_dtdM(0,MV,UV));
	double df13dm2 = -(t_cor * dv_Ldpy_P + v_L * Calc_dtdM(1,MV,UV));
	double df13dm3 = -(t_cor * dv_Ldpz_P + v_L * Calc_dtdM(2,MV,UV));
	double df13dm4 = -v_L * Calc_dtdM(3,MV,UV);
	double df13dm5 = -0.5 - v_L * Calc_dtdM(4,MV,UV);
	double df13dm6 = -v_L * Calc_dtdM(5,MV,UV);
	double df13dm7 = -(t_cor * dv_Ldpx_Q + v_L * Calc_dtdM(6,MV,UV));
	double df13dm8 = -(t_cor * dv_Ldpy_Q + v_L * Calc_dtdM(7,MV,UV));
	double df13dm9 = -(t_cor * dv_Ldpz_Q + v_L * Calc_dtdM(8,MV,UV));
	double df13dm10= -v_L * Calc_dtdM(9,MV,UV);
	double df13dm11= -0.5 - v_L * Calc_dtdM(10,MV,UV);
	double df13dm12= -v_L * Calc_dtdM(11,MV,UV);
	double df13dm13= 0, df13dm14= 0, df13dm15= 0;
	double df13dm16= -v_L * Calc_dtdM(15,MV,UV);
	double df13dm17= 1 - v_L * Calc_dtdM(16,MV,UV);
	double df13dm18= -v_L * Calc_dtdM(17,MV,UV);
	//f14: vz_R - 0.5*(vz_P + vz_Q) - t_cor*w_L
	double df14du1 = 0, df14du2 = 0, df14du3 = 0;
	double df14du4 = 0, df14du5 = 0, df14du6 = 0;
	double df14dm1 = -(t_cor * dw_Ldpx_P + w_L * Calc_dtdM(0,MV,UV));
	double df14dm2 = -(t_cor * dw_Ldpy_P + w_L * Calc_dtdM(1,MV,UV));
	double df14dm3 = -(t_cor * dw_Ldpz_P + w_L * Calc_dtdM(2,MV,UV));
	double df14dm4 = -w_L * Calc_dtdM(3,MV,UV);
	double df14dm5 = -w_L * Calc_dtdM(4,MV,UV);
	double df14dm6 = -0.5 - w_L * Calc_dtdM(5,MV,UV);
	double df14dm7 = -(t_cor * dw_Ldpx_Q + w_L * Calc_dtdM(6,MV,UV));
	double df14dm8 = -(t_cor * dw_Ldpy_Q + w_L * Calc_dtdM(7,MV,UV));
	double df14dm9 = -(t_cor * dw_Ldpz_Q + w_L * Calc_dtdM(8,MV,UV));
	double df14dm10= -w_L * Calc_dtdM(9,MV,UV);
	double df14dm11= -w_L * Calc_dtdM(10,MV,UV);
	double df14dm12= -0.5 - w_L * Calc_dtdM(11,MV,UV);
	double df14dm13= 0, df14dm14= 0, df14dm15= 0;
	double df14dm16= -w_L * Calc_dtdM(15,MV,UV);
	double df14dm17= -w_L * Calc_dtdM(16,MV,UV);
	double df14dm18= 1 - w_L * Calc_dtdM(17,MV,UV);
	#endif
	//f12 = (V_R - V_LX)*E_1
	double df12du1 = 0, df12du2 = 0, df12du3 = 0;
	double df12du4 = 0, df12du5 = 0, df12du6 = 0;
	double df12dm1 = -Calc_dV_LXdM(0,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(0,MV);
	double df12dm2 = -Calc_dV_LXdM(1,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(1,MV);
	double df12dm3 = -Calc_dV_LXdM(2,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(2,MV);
	double df12dm4 = -Calc_dV_LXdM(3,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(3,MV);
	double df12dm5 = -Calc_dV_LXdM(4,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(4,MV);
	double df12dm6 = -Calc_dV_LXdM(5,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(5,MV);
	double df12dm7 = -Calc_dV_LXdM(6,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(6,MV);
	double df12dm8 = -Calc_dV_LXdM(7,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(7,MV);
	double df12dm9 = -Calc_dV_LXdM(8,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(8,MV);
	double df12dm10= -Calc_dV_LXdM(9,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(9,MV);
	double df12dm11= -Calc_dV_LXdM(10,MV)*E_1+ (V_R - V_LX) *Calc_dE1dM(10,MV);
	double df12dm12= -Calc_dV_LXdM(11,MV)*E_1+ (V_R - V_LX) *Calc_dE1dM(11,MV);
	double df12dm13= 0, df12dm14= 0, df12dm15= 0;
	double df12dm16= E_1.x() -Calc_dV_LXdM(15,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(15,MV);
	double df12dm17= E_1.y() -Calc_dV_LXdM(16,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(16,MV);
	double df12dm18= E_1.z() -Calc_dV_LXdM(17,MV)*E_1 + (V_R - V_LX) * Calc_dE1dM(17,MV);
	//f13 = (V_R - V_LX)*E_2
	double df13du1 = 0, df13du2 = 0, df13du3 = 0;
	double df13du4 = 0, df13du5 = 0, df13du6 = 0;
	double df13dm1 = -Calc_dV_LXdM(0,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(0,MV);
	double df13dm2 = -Calc_dV_LXdM(1,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(1,MV);
	double df13dm3 = -Calc_dV_LXdM(2,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(2,MV);
	double df13dm4 = -Calc_dV_LXdM(3,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(3,MV);
	double df13dm5 = -Calc_dV_LXdM(4,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(4,MV);
	double df13dm6 = -Calc_dV_LXdM(5,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(5,MV);
	double df13dm7 = -Calc_dV_LXdM(6,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(6,MV);
	double df13dm8 = -Calc_dV_LXdM(7,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(7,MV);
	double df13dm9 = -Calc_dV_LXdM(8,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(8,MV);
	double df13dm10= -Calc_dV_LXdM(9,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(9,MV);
	double df13dm11= -Calc_dV_LXdM(10,MV)*E_2+ (V_R - V_LX) *Calc_dE2dM(10,MV);
	double df13dm12= -Calc_dV_LXdM(11,MV)*E_2+ (V_R - V_LX) *Calc_dE2dM(11,MV);
	double df13dm13= 0, df13dm14= 0, df13dm15= 0;
	double df13dm16= E_2.x() -Calc_dV_LXdM(15,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(15,MV);
	double df13dm17= E_2.y() -Calc_dV_LXdM(16,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(16,MV);
	double df13dm18= E_2.z() -Calc_dV_LXdM(17,MV)*E_2 + (V_R - V_LX) * Calc_dE2dM(17,MV);
	
	//double fs[]={f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14};
	double fs[]={f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13};
	double dfdms[400] ;
	double dfdus[400] ;
	if(UseVertexFlag){
	}
	else{
		double temp[] = {// 14 constraints, 18 measurement params. =  252 elements
			df1dm1, df1dm2, df1dm3, df1dm4	,df1dm5, df1dm6, df1dm7, df1dm8, df1dm9, df1dm10, df1dm11, df1dm12, df1dm13, df1dm14, df1dm15, df1dm16, df1dm17, df1dm18,
			df2dm1, df2dm2, df2dm3, df2dm4	,df2dm5, df2dm6, df2dm7, df2dm8, df2dm9, df2dm10, df2dm11, df2dm12, df2dm13, df2dm14, df2dm15, df2dm16, df2dm17, df2dm18,
			df3dm1, df3dm2, df3dm3, df3dm4	,df3dm5, df3dm6, df3dm7, df3dm8, df3dm9, df3dm10, df3dm11, df3dm12, df3dm13, df3dm14, df3dm15, df3dm16, df3dm17, df3dm18,
			df4dm1, df4dm2, df4dm3, df4dm4	,df4dm5, df4dm6, df4dm7, df4dm8, df4dm9, df4dm10, df4dm11, df4dm12, df4dm13, df4dm14, df4dm15, df4dm16, df4dm17, df4dm18,
			df5dm1, df5dm2, df5dm3, df5dm4	,df5dm5, df5dm6, df5dm7, df5dm8, df5dm9, df5dm10, df5dm11, df5dm12, df5dm13, df5dm14, df5dm15, df5dm16, df5dm17, df5dm18,
			df6dm1, df6dm2, df6dm3, df6dm4	,df6dm5, df6dm6, df6dm7, df6dm8, df6dm9, df6dm10, df6dm11, df6dm12, df6dm13, df6dm14, df6dm15, df6dm16, df6dm17, df6dm18,
			df7dm1, df7dm2, df7dm3, df7dm4	,df7dm5, df7dm6, df7dm7, df7dm8, df7dm9, df7dm10, df7dm11, df7dm12, df7dm13, df7dm14, df7dm15, df7dm16, df7dm17, df7dm18,
			df8dm1, df8dm2, df8dm3, df8dm4	,df8dm5, df8dm6, df8dm7, df8dm8, df8dm9, df8dm10, df8dm11, df8dm12, df8dm13, df8dm14, df8dm15, df8dm16, df8dm17, df8dm18,
			df9dm1, df9dm2, df9dm3, df9dm4	,df9dm5, df9dm6, df9dm7, df9dm8, df9dm9, df9dm10, df9dm11, df9dm12, df9dm13, df9dm14, df9dm15, df9dm16, df9dm17, df9dm18,
			df10dm1,df10dm2,df10dm3,df10dm4	,df10dm5,df10dm6,df10dm7,df10dm8,df10dm9,df10dm10,df10dm11,df10dm12,df10dm13,df10dm14,df10dm15,df10dm16,df10dm17,df10dm18,
			df11dm1,df11dm2,df11dm3,df11dm4	,df11dm5,df11dm6,df11dm7,df11dm8,df11dm9,df11dm10,df11dm11,df11dm12,df11dm13,df11dm14,df11dm15,df11dm16,df11dm17,df11dm18,
			df12dm1,df12dm2,df12dm3,df12dm4	,df12dm5,df12dm6,df12dm7,df12dm8,df12dm9,df12dm10,df12dm11,df12dm12,df12dm13,df12dm14,df12dm15,df12dm16,df12dm17,df12dm18,
			df13dm1,df13dm2,df13dm3,df13dm4	,df13dm5,df13dm6,df13dm7,df13dm8,df13dm9,df13dm10,df13dm11,df13dm12,df13dm13,df13dm14,df13dm15,df13dm16,df13dm17,df13dm18
			//df14dm1,df14dm2,df14dm3,df14dm4	,df14dm5,df14dm6,df14dm7,df14dm8,df14dm9,df14dm10,df14dm11,df14dm12,df14dm13,df14dm14,df14dm15,df14dm16,df14dm17,df14dm18
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = { // 14 constraints, 6 unknown params. = 84 elements
			df1du1, df1du2, df1du3, df1du4, df1du5, df1du6,
			df2du1, df2du2, df2du3, df2du4, df2du5, df2du6,
			df3du1, df3du2, df3du3, df3du4, df3du5, df3du6,
			df4du1, df4du2, df4du3, df4du4, df4du5, df4du6,
			df5du1, df5du2, df5du3, df5du4, df5du5, df5du6,
			df6du1, df6du2, df6du3, df6du4, df6du5, df6du6,
			df7du1, df7du2, df7du3, df7du4, df7du5, df7du6,
			df8du1, df8du2, df8du3, df8du4, df8du5, df8du6,
			df9du1, df9du2, df9du3, df9du4, df9du5, df9du6,
			df10du1,df10du2,df10du3,df10du4,df10du5,df10du6,
			df11du1,df11du2,df11du3,df11du4,df11du5,df11du6,
			df12du1,df12du2,df12du3,df12du4,df12du5,df12du6,
			df13du1,df13du2,df13du3,df13du4,df13du5,df13du6
			//df14du1,df14du2,df14du3,df14du4,df14du5,df14du6
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
void CascadeVertexFitter::SampleStepPoint(int steps){
#if DebugCKF
	std::cout<<"CascadeVertexFitter::SampleStepPoint() step = "<<steps<<std::endl;
#endif
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double px_P,py_P,pz_P,px_Q,py_Q,pz_Q,px_R,py_R,pz_R,px_X,py_X,pz_X;
	double vx_R,vy_R,vz_R,vx_P,vy_P,vz_P,vx_Q,vy_Q,vz_Q,vx_X,vy_X,vz_X;
	if(UseVertexFlag){
	}
	else{
		px_P = Meas(0,0); 
		py_P = Meas(1,0); 
		pz_P = Meas(2,0);
		vx_P= Meas(3,0);
		vy_P= Meas(4,0);
		vz_P= Meas(5,0);

		px_Q = Meas(6,0); 
		py_Q= Meas(7,0); 
		pz_Q= Meas(8,0); 
		vx_Q= Meas(9,0);
		vy_Q= Meas(10,0);
		vz_Q= Meas(11,0);

		px_R = Meas(12,0);
		py_R= Meas(13,0);
		pz_R= Meas(14,0);
		vx_R= Meas(15,0);
		vy_R= Meas(16,0);
		vz_R= Meas(17,0);

		px_X = Unkn(0,0);
		py_X= Unkn(1,0);
		pz_X= Unkn(2,0);
		vx_X= Unkn(3,0);
		vy_X= Unkn(4,0);
		vz_X= Unkn(5,0);
	}
	TVector3 TV_P(px_P,py_P,pz_P);
	TVector3 TV_Q(px_Q,py_Q,pz_Q);
	TVector3 TV_R(px_R,py_R,pz_R);
	TVector3 TV_L = TV_P + TV_Q;
	TVector3 TV_X(px_X,py_X,pz_X);
	double px_L = TV_L.X();double py_L = TV_L.Y();double pz_L = TV_L.Z();
	double p_P = TV_P.Mag();double p_Q = TV_Q.Mag();double p_R = TV_R.Mag();double p_L = TV_L.Mag();double p_X = TV_X.Mag();

	double vx_L= 0.5*(vx_P + vx_Q);double vy_L= 0.5*(vy_P + vy_Q);double vz_L= 0.5*(vz_P + vz_Q);
	VPCor = TVector3(vx_P,vy_P,vz_P);
	VQCor = TVector3(vx_Q,vy_Q,vz_Q);
	VRCor = TVector3(vx_R,vy_R,vz_R);
	VLCor = TVector3(vx_L,vy_L,vz_L);
	VXCor = TVector3(vx_X,vy_X,vz_X);
	vector<double> MV = {
	px_P,py_P,pz_P,vx_P,vy_P,vz_P,
	px_Q,py_Q,pz_Q,vx_Q,vy_Q,vz_Q,
	px_R,py_R,pz_R,vx_R,vy_R,vz_R
	};
	vector<double> UV = {
	px_X,py_X,pz_X,
	vx_X,vy_X,vz_X
	};
	t_cor = CalcLambdaExtrapolationParameter(MV,UV);


	
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
TMatrixD
CascadeVertexFitter::JacobianSphToCart(double p, double th, double ph){//Legacy
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
CascadeVertexFitter::CalcVariance(int istep){//Legacy
	return;
}








void
CascadeVertexFitter::Rotate(){//Legacy
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
CascadeVertexFitter::ToDecayPlane(){//Legacy
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
CascadeVertexFitter::CalcLambdaExtrapolationParameter(vector<double> Meas_,vector<double> Unkn_){
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
CascadeVertexFitter::CalcLambdaDirectionalDerivativesM(int idx, vector<double> Meas_){
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
CascadeVertexFitter::CalcLambdaExtrapolationParameterDerivativesM(int idx, vector<double> Meas_, vector<double> Unkn_){
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
CascadeVertexFitter::CalcLambdaExtrapolationParameterDerivativesU(int idx, vector<double> Meas_, vector<double> Unkn_){
	//A function to calculate the derivatives of the Lambda extrapolation parameter t by the Unkn parameters.
	//Not used here, since L vtx are not unkn parameters in this code. 
	return 0;
}
TVector3
CascadeVertexFitter::CalcV_LX(vector<double> Meas_){
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
CascadeVertexFitter::CalcV_LXDerivativesM(int idx, vector<double> Meas_){
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
CascadeVertexFitter::CalcLambdaE1Vector(vector<double> Meas_){
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
CascadeVertexFitter::CalcLambdaE1VectorDerivativesM(int idx, vector<double> Meas_){
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
CascadeVertexFitter::CalcLambdaE2Vector(vector<double> Meas_){
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
CascadeVertexFitter::CalcLambdaE2VectorDerivativesM(int idx, vector<double> Meas_){
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
