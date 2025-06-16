#ifndef KEKKinematicFit_c
#define KEKKinematicFit_c
#include "KEKKinematicFit.hh"

KEKFourVectorFitter::KEKFourVectorFitter(
	TLorentzVector LV1, TMatrixD Cov1,
	TLorentzVector LV2, TMatrixD Cov2,double Mass
){
	double Diagonals[6] = {
		Cov1[0][0],Cov1[1][1],Cov1[2][2],Cov2[0][0],Cov2[1][1],Cov2[2][2]
	};
	auto Offdiagonals = MathTools::MergeOffdiagonals(Cov1,Cov2);
	//In KF framework, Y and Z should be swapped
	auto HLV1 = TLorentzVector(LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
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
	cout<<"Dims of Cov1: "<<Cov1.GetNrows()<<"x"<<Cov1.GetNcols()<<endl;
	cout<<"Dims of Cov2: "<<Cov2.GetNrows()<<"x"<<Cov2.GetNcols()<<endl;
	cout<<"Dims of Cov3: "<<Cov3.GetNrows()<<"x"<<Cov3.GetNcols()<<endl;
	cout<<"Dims of Cov12: "<<Cov12.GetNrows()<<"x"<<Cov12.GetNcols()<<endl;
	cout<<"Dims of Offdiagonals: "<<Offdiagonals.GetNrows()<<"x"<<Offdiagonals.GetNcols()<<endl;
	auto HLV1 = TLorentzVector(LV1.Px(), LV1.Pz(), LV1.Py(), LV1.E());
	auto HLV2 = TLorentzVector(LV2.Px(), LV2.Pz(), LV2.Py(), LV2.E());
	auto HLV3 = TLorentzVector(LV3.Px(), LV3.Pz(), LV3.Py(), LV3.E());
	Fitter = new CascadeFitter(HLV1,HLV2,HLV3);
	Fitter->ScaleParameters(ScaleParams);
	ThisFitter()->SetInvMass(Mass1,Mass2);
	Fitter->SetMaximumStep(5);
	Fitter->SetVariance(Diagonals);
	Fitter->AddOffdiagonals(Offdiagonals);
}
#endif