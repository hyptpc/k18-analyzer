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

#endif