#ifndef KEKKinematicFit_c
#define KEKKinematicFit_c
#include "KEKKinematicFit.hh"
#include "UserParamMan.hh"
#include "ConfMan.hh"
namespace
{
  const auto& gUser = UserParamMan::GetInstance();
}
KEKFourVectorFitter::KEKFourVectorFitter(
	TLorentzVector LV1, TMatrixD Cov1,
	TLorentzVector LV2, TMatrixD Cov2,double Mass
){
	double Diagonals[6] = {
		Cov1[0][0],Cov1[1][1],Cov1[2][2],Cov2[0][0],Cov2[1][1],Cov2[2][2]
	};
	auto Offdiagonals = MathTools::MergeOffdiagonals(Cov1,Cov2);
	//In KF framework, Y and Z should be swapped
	//260529: KF coordinate is now (-x,z,y)
	auto HLV1 = TLorentzVector(-LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(-LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
	auto HLV3 = HLV1 + HLV2;
	Fitter = new FourVectorFitter(HLV1,HLV2,HLV3);
	Fitter->ScaleParameters(ScaleParams);
	ThisFitter()->SetInvMass(Mass);
	Fitter->SetMaximumStep(5);
	Fitter->SetVariance(Diagonals);
	Fitter->AddOffdiagonals(Offdiagonals);
}
KEKCascadeFitter::KEKCascadeFitter(
	TLorentzVector LV1, TMatrixD Cov1,
	TLorentzVector LV2, TMatrixD Cov2,
	TLorentzVector LV3, TMatrixD Cov3,
	double Mass1, double Mass2
){
	double Diagonals[9] = {
		Cov1[0][0],Cov1[1][1],Cov1[2][2],
		Cov2[0][0],Cov2[1][1],Cov2[2][2],
		Cov3[0][0],Cov3[1][1],Cov3[2][2]
	};
	auto Cov12 = MathTools::MergeOffdiagonals(Cov1,Cov2);
	auto Offdiagonals = MathTools::MergeOffdiagonals(Cov12,Cov3);
	auto HLV1 = TLorentzVector(-LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(-LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
	auto HLV3 = TLorentzVector(-LV3.Px(), LV3.Pz(), LV3.Py(), LV3.E());
	Fitter = new CascadeFitter(HLV1,HLV2,HLV3);
	Fitter->ScaleParameters(ScaleParams);
	ThisFitter()->SetInvMass(Mass1,Mass2);
	Fitter->SetMaximumStep(5);
	Fitter->SetVariance(Diagonals);
	Fitter->AddOffdiagonals(Offdiagonals);
}
KEKCascadeFitter2::KEKCascadeFitter2(
	TLorentzVector LV1, TMatrixD Cov1,
	TLorentzVector LV2, TMatrixD Cov2,
	TLorentzVector LV3, TMatrixD Cov3,
	double Mass1, double Mass2
){
	#if 0
	static const Double_t MomResScale = gUser.GetParameter("MomResScale") ;
	static const Double_t dZResScale = gUser.GetParameter("dZResScale") ;
	#else
	static const Double_t MomResScale  = 2.4;
	static const Double_t dZResScale   = 1;
	#endif
	TMatrixD MomCov1(3,3);
	TMatrixD MomCov2(3,3);
	TMatrixD MomCov3(3,3);
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			MomCov1(i,j) = Cov1(i,j);
			MomCov2(i,j) = Cov2(i,j);
			MomCov3(i,j) = Cov3(i,j);
		}
	}
	TMatrixD MomScaleMat(3,3);
	MomScaleMat.Zero();
	MomScaleMat(0,0) = MomResScale;//px
	MomScaleMat(1,1) = MomResScale*dZResScale;//py
	MomScaleMat(2,2) = MomResScale;//pz
	MomCov1 = MomScaleMat*MomCov1*MomScaleMat;
	MomCov2 = MomScaleMat*MomCov2*MomScaleMat;
	MomCov3 = MomScaleMat*MomCov3*MomScaleMat;
	double Diagonals[9] = {
		MomCov1[0][0],MomCov1[1][1],MomCov1[2][2],
		MomCov2[0][0],MomCov2[1][1],MomCov2[2][2],
		MomCov3[0][0],MomCov3[1][1],MomCov3[2][2]
	};
	auto MomCov12 = MathTools::MergeOffdiagonals(MomCov1,MomCov2);
	auto Offdiagonals = MathTools::MergeOffdiagonals(MomCov12,MomCov3);
	double scale_offdiagonals = 0.9;//To avoid singular matrix, the off-diagonal elements are scaled down by this factor.
	for(int i=0;i<Offdiagonals.GetNrows();++i){
		for(int j=0;j<Offdiagonals.GetNcols();++j){
			Offdiagonals(i,j) *= scale_offdiagonals;
		}
	}

	auto HLV1 = TLorentzVector(-LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(-LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
	auto HLV3 = TLorentzVector(-LV3.Px(), LV3.Pz(), LV3.Py(), LV3.E());
	Fitter = new CascadeFitter2(HLV1,HLV2,HLV3);
	Fitter->ScaleParameters(ScaleParams);
	ThisFitter()->SetInvMass(Mass1,Mass2);
	Fitter->SetMaximumStep(5);
	Fitter->SetVariance(Diagonals);
	Fitter->AddOffdiagonals(Offdiagonals);
}
KEKCascadeVertexFitter::KEKCascadeVertexFitter(
	TLorentzVector LV1, TVector3 Vtx1, TMatrixD Cov1,
	TLorentzVector LV2, TVector3 Vtx2, TMatrixD Cov2,
	TLorentzVector LV3, TVector3 Vtx3, TMatrixD Cov3,
	double Mass1, double Mass2
){

	auto HLV1 = TLorentzVector(-LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(-LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
	auto HLV3 = TLorentzVector(-LV3.Px(), LV3.Pz(), LV3.Py(), LV3.E());
	auto HVtx1= TVector3(-Vtx1.X(), Vtx1.Z(), Vtx1.Y());
	auto HVtx2= TVector3(-Vtx2.X(), Vtx2.Z(), Vtx2.Y());
	auto HVtx3= TVector3(-Vtx3.X(), Vtx3.Z(), Vtx3.Y());


	#if 0
	static const Double_t MomResScale = gUser.GetParameter("MomResScale") ;
	static const Double_t dZResScale = gUser.GetParameter("dZResScale") ;
	static const Double_t PhiResScale = gUser.GetParameter("PhiResScale") ;
	static const Double_t TransResScale = gUser.GetParameter("TransResScale") ;
	static const Double_t VertResScale = gUser.GetParameter("VertResScale") ;
	#else
	static const Double_t MomResScale   = 2.4;
	static const Double_t dZResScale    = 1;
	/*
	static const Double_t ResXScaleP    = 1.5;
	static const Double_t ResYScaleP    = 1.5;
	static const Double_t ResZScaleP    = 0.3;
	static const Double_t ResXScalePi1  = 1.1;
	static const Double_t ResYScalePi1  = 1.2;
	static const Double_t ResZScalePi1  = 0.4;
	static const Double_t ResXScalePi2  = 2.1;
	static const Double_t ResYScalePi2  = 2.2;
	static const Double_t ResZScalePi2  = 2;
	*/
	static const Double_t ResXScaleP    = 1;
	static const Double_t ResYScaleP    = 1;
	static const Double_t ResZScaleP    = 1;
	static const Double_t ResXScalePi1  = 1;
	static const Double_t ResYScalePi1  = 1;
	static const Double_t ResZScalePi1  = 1;
	static const Double_t ResXScalePi2  = 1;
	static const Double_t ResYScalePi2  = 1;
	static const Double_t ResZScalePi2  = 1;
	#endif
	/*
	Each point in the track has 6 informations, (x,y,z, px, py, pz).
	However, if we 'select' a point, we lose one degree of freedom. 
	Hence the covariance matrix of the selected point has a rank of 5,
	leading to a singularity in (x,y,z,px,py,pz) representation.
	In the covariance at the Vertex point, there is another degree of freedom
	in the 'selection' of the point: We do not 'know' the exact position of the vertex.
	Let the 'true' information be X0, 'measured' be X, the track parameter at the measured be t0.
	If we take a linear approximation, the extrapolated information at the vertex can be written as:
	X' = X + U (t - t0)
	where U is the momentum direction.
	If we assume some linear approximation, the covariance of X' can be approximated as:
	Cov(X') = Cov(X) + dt^2 *U*U^T
	where the cross term X * U dt are negleced.
	The scale parameter, dt, can be estimated from the vertex counterpart, since the vertex resolution
	is dominated by the extrapolation uncertainty. If the position part of the counterpart tarck is V,
	dt can be defined as:
	dt^2 = U * V * U^T 
	* 20260607 Modification:
	dt^2 = U * V * U^T * 1./sin(opening angle)^2.
	A variance in U will induce a variance proportional tu U* 1/sin(opening angle) in the counterpart plane.
	
	For pi2, it has no counterpart.
	We take the average of p pi position covariance.
	*/
	TVector3 U1 = HLV1.Vect().Unit();
	TVector3 U2 = HLV2.Vect().Unit();
	TVector3 U3 = HLV3.Vect().Unit();

	double cth_opening_12 = U1.Dot(U2);
	double sth_opening_12 = sqrt(1 - cth_opening_12*cth_opening_12);
	//The sign does not matter. We will take the squared value.
	//Maybe a numerical cap may be required?
	double cth_opening_3 = U3.Dot((U1+U2).Unit());
	double sth_opening_3 = sqrt(1 - cth_opening_3*cth_opening_3);

	TMatrixD UMat1(3,1); TMatrixD UMatT1(1,3);
	TMatrixD UMat2(3,1); TMatrixD UMatT2(1,3);
	TMatrixD UMat3(3,1); TMatrixD UMatT3(1,3);
	for(int i=0;i<3;++i){
		UMat1(i,0) = U1(i);UMatT1(0,i) = U1(i);
		UMat2(i,0) = U2(i);UMatT2(0,i) = U2(i);
		UMat3(i,0) = U3(i);UMatT3(0,i) = U3(i);
	}

	TMatrixD V1(3,3);TMatrixD V2(3,3);TMatrixD V3(3,3);
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			V1(i,j) = Cov1(i+3,j+3);
			V2(i,j) = Cov2(i+3,j+3);
			V3(i,j) = 0.5 * (Cov1(i+3,j+3) + Cov2(i+3,j+3));
		}
	}
	double dt1 = sqrt((UMatT1*V2*UMat1)(0,0))*sth_opening_12;
	double dt2 = sqrt((UMatT2*V1*UMat2)(0,0))*sth_opening_12;
	double dt3 = sqrt((UMatT3*V3*UMat3)(0,0))*sth_opening_3;
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			Cov1(i+3,j+3) += dt1*dt1*U1(i)*U1(j);
			Cov2(i+3,j+3) += dt2*dt2*U2(i)*U2(j);
			Cov3(i+3,j+3) += dt3*dt3*U3(i)*U3(j);
		}
	}

	TMatrixD MomPosScaleMatP(6,6);
	MomPosScaleMatP.Zero();
	MomPosScaleMatP(0,0) = MomResScale;//px
	MomPosScaleMatP(1,1) = MomResScale;//pz = py in helix coordinate
	MomPosScaleMatP(2,2) = MomResScale*dZResScale;//
	MomPosScaleMatP(3,3) = ResXScaleP;//x
	MomPosScaleMatP(4,4) = ResZScaleP;//z
	MomPosScaleMatP(5,5) = ResYScaleP;//y
	TMatrixD MomPosScaleMatPi1(6,6);
	MomPosScaleMatPi1.Zero();
	MomPosScaleMatPi1(0,0) = MomResScale;//px
	MomPosScaleMatPi1(1,1) = MomResScale;//pz
	MomPosScaleMatPi1(2,2) = MomResScale*dZResScale;//py
	MomPosScaleMatPi1(3,3) = ResXScalePi1;
	MomPosScaleMatPi1(4,4) = ResZScalePi1;//z
	MomPosScaleMatPi1(5,5) = ResYScalePi1;//y
	TMatrixD MomPosScaleMatPi2(6,6);
	MomPosScaleMatPi2.Zero();
	MomPosScaleMatPi2(0,0) = MomResScale;//px
	MomPosScaleMatPi2(1,1) = MomResScale;//pz
	MomPosScaleMatPi2(2,2) = MomResScale*dZResScale;//py
	MomPosScaleMatPi2(3,3) = ResXScalePi2;//x
	MomPosScaleMatPi2(4,4) = ResZScalePi2;//z
	MomPosScaleMatPi2(5,5) = ResYScalePi2;//y

	Cov1 = MomPosScaleMatP*Cov1*MomPosScaleMatP;
	Cov2 = MomPosScaleMatPi1*Cov2*MomPosScaleMatPi1;
	Cov3 = MomPosScaleMatPi2*Cov3*MomPosScaleMatPi2;
	double Diagonals[18] = {
		Cov1[0][0],Cov1[1][1],Cov1[2][2],Cov1[3][3],Cov1[4][4],Cov1[5][5],
		Cov2[0][0],Cov2[1][1],Cov2[2][2],Cov2[3][3],Cov2[4][4],Cov2[5][5],
		Cov3[0][0],Cov3[1][1],Cov3[2][2],Cov3[3][3],Cov3[4][4],Cov3[5][5]
	};
	auto Cov12 = MathTools::MergeOffdiagonals(Cov1,Cov2);
	auto Offdiagonals = MathTools::MergeOffdiagonals(Cov12,Cov3);
	double scale_offdiagonals = 0.9;//
	for(int i=0;i<Offdiagonals.GetNrows();++i){
		for(int j=0;j<Offdiagonals.GetNcols();++j){
	//		Offdiagonals(i,j) *= scale_offdiagonals;
		}
	}
	Fitter = new CascadeVertexFitter(HLV1,HVtx1,HLV2,HVtx2,HLV3,HVtx3);
	//Fitter->ScaleParameters(ScaleParams);
	Fitter->ScaleParameters(0);
	ThisFitter()->SetInvMass(Mass1,Mass2);
	Fitter->SetMaximumStep(5);
	Fitter->SetVariance(Diagonals);
	Fitter->AddOffdiagonals(Offdiagonals);
}
KEKMassVertexFitter::KEKMassVertexFitter(
	TLorentzVector LV1, TVector3 Vtx1, TMatrixD Cov1,
	TLorentzVector LV2, TVector3 Vtx2, TMatrixD Cov2,double Mass
){
	double Diagonals[10] = {
		Cov1[0][0],Cov1[1][1],Cov1[2][2],Cov1[3][3],Cov1[4][4],
		Cov2[0][0],Cov2[1][1],Cov2[2][2],Cov2[3][3],Cov2[4][4]
	};
	auto Offdiagonals = MathTools::MergeOffdiagonals(Cov1,Cov2);
	//In KF framework, Y and Z should be swapped
	auto HLV1 = TLorentzVector(-LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(-LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
	auto HLV3 = HLV1 + HLV2;
	Fitter = new MassVertexFitter(HLV1,HLV2,HLV3,Vtx1,Vtx2);
	Fitter->ScaleParameters(ScaleParams);
	ThisFitter()->SetInvMass(Mass);
	Fitter->SetMaximumStep(5);
	Fitter->SetVariance(Diagonals);
	Fitter->AddOffdiagonals(Offdiagonals);
}
KEKMassVertexFitter3::KEKMassVertexFitter3(
	TLorentzVector LV1, TVector3 Vtx1, TMatrixD Cov1,
	TLorentzVector LV2, TVector3 Vtx2, TMatrixD Cov2,
	double Mass
){
	double Diagonals[15] = {
		Cov1[0][0],Cov1[1][1],Cov1[2][2],Cov1[3][3],Cov1[4][4],
		Cov2[0][0],Cov2[1][1],Cov2[2][2],Cov2[3][3],Cov2[4][4]
	};
	auto Offdiagonals = MathTools::MergeOffdiagonals(Cov1,Cov2);
	auto HLV1 = TLorentzVector(-LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(-LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
	auto HLV3 = HLV1 + HLV2; 
	Fitter = new MassVertexFitter3(HLV1,HLV2,HLV3,Vtx1,Vtx2);
	Fitter->ScaleParameters(ScaleParams);
	ThisFitter()->SetInvMass(Mass);
	Fitter->SetMaximumStep(5);
	Fitter->SetVariance(Diagonals);
	Fitter->AddOffdiagonals(Offdiagonals);
}
#endif