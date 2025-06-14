#ifndef KEKKinematicFit_h
#define KEKKinematicFit_h
#include "FourVectorFitter.hh"
#include "MassVertexFitter.hh"
//#include "CascadeFitter.hh"
#include "TPCLocalTrackHelix.hh"
#include "MathTools.hh"
class KEKKinematicFitter{
	protected:
		int nSteps = 5;
		int nDoF = -1;
		double Chi2 = -1;
		double Pvalue = -1;
		bool ScaleParams = 1;
		KinematicFitter* Fitter = nullptr;
		vector<double> Pulls;
	public:
		KEKKinematicFitter(){}
		~KEKKinematicFitter(){
		}
		void ScaleVaraince(bool status = true){
			ScaleParams = status;
			Fitter->ScaleParameters(ScaleParams);
		}
		void SetMaximumStep(int nsteps){
			nSteps = nsteps;
			Fitter->SetMaximumStep(nSteps);
		}
		int GetNSteps(){
			return nSteps;
		}
		int GetNDoF(){
			return nDoF;
		}
		TMatrixD GetVariance(int i = 0){
			return Fitter->GetVariance(0);
		}
		double GetChi2(){
			return Chi2;
		}
		double GetPValue(){
			return Pvalue;
		}
		vector<double> GetPull(){
			return Pulls;
		}
		TMatrixD GetUnmeasuredCovariance(){
			return Fitter->GetUnmeasuredCovariance();
		};
		double DoKinematicFit(){
			Chi2 = Fitter->DoKinematicFit();
			Pvalue = Fitter->GetPValue();
			Pulls = Fitter->GetPull();
			nDoF = Fitter->GetNDF();
			return Chi2;
		}
};
class KEKFourVectorFitter: virtual public KEKKinematicFitter{
	private:
		TLorentzVector LV1Cor;
		TLorentzVector LV2Cor;
		TLorentzVector LV3Cor;
	public:
		KEKFourVectorFitter(){}
		FourVectorFitter* ThisFitter(){//Necessary to call functions in FourVectorFitter. ex: SetInvMass.
			return dynamic_cast<FourVectorFitter*>(Fitter);
		}
		KEKFourVectorFitter(
		TLorentzVector LV1, TMatrixD Cov1,
		TLorentzVector LV2, TMatrixD Cov2,double Mass);
		std::vector<TLorentzVector> GetFittedLV(){
			auto HLVs = ThisFitter()->GetFittedLV();
			auto HLV1Cor = HLVs.at(0);
			auto HLV2Cor = HLVs.at(1);
			auto HLV3Cor = HLVs.at(2);
			LV1Cor = TLorentzVector(HLV1Cor.Px(), HLV1Cor.Pz(), HLV1Cor.Py(), HLV1Cor.E());
			LV2Cor = TLorentzVector(HLVs.at(1).Px(), HLVs.at(1).Pz(), HLVs.at(1).Py(), HLVs.at(1).E());
			LV3Cor = TLorentzVector(HLVs.at(2).Px(), HLVs.at(2).Pz(), HLVs.at(2).Py(), HLVs.at(2).E());
			return std::vector<TLorentzVector>{LV1Cor, LV2Cor, LV3Cor};
		}
		~KEKFourVectorFitter(){
		}
};
#endif
