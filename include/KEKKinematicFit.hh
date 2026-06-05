#ifndef KEKKinematicFit_h
#define KEKKinematicFit_h
#include "FourVectorFitter.hh"
#include "MassVertexFitter.hh"
#include "CascadeFitter.hh"
#include "CascadeFitter2.hh"
#include "CascadeVertexFitter.hh"
#include "MassVertexFitter3.hh"
#include "TPCLocalTrackHelix.hh"
#include "MathTools.hh"
class KEKKinematicFitter{
	protected:
		int nSteps = 5;
		int nDoF = -1;
		double Chi2 = -1;
		double Pvalue = -1;
		bool ScaleParams = 0;
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
			return Fitter->GetVariance(i);
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
			LV1Cor = TLorentzVector(-HLV1Cor.Px(), HLV1Cor.Pz(), HLV1Cor.Py(), HLV1Cor.E());
			LV2Cor = TLorentzVector(-HLVs.at(1).Px(), HLVs.at(1).Pz(), HLVs.at(1).Py(), HLVs.at(1).E());
			LV3Cor = TLorentzVector(-HLVs.at(2).Px(), HLVs.at(2).Pz(), HLVs.at(2).Py(), HLVs.at(2).E());
			return std::vector<TLorentzVector>{LV1Cor, LV2Cor, LV3Cor};
		}
		~KEKFourVectorFitter(){
		}
};
class KEKCascadeFitter: virtual public KEKKinematicFitter{
	private:
		TLorentzVector LV1Cor;
		TLorentzVector LV2Cor;
		TLorentzVector LV3Cor;//Decay products, eg: p,pi_L, pi_Xi
		TLorentzVector LV4Cor;//Intermediate products, eg: Lambda
		TLorentzVector LV5Cor;//Initial state, eg: Xi
	public:
		KEKCascadeFitter(){}
		CascadeFitter* ThisFitter(){
			return dynamic_cast<CascadeFitter*>(Fitter);
		}
		KEKCascadeFitter(
		TLorentzVector LV1, TMatrixD Cov1,
		TLorentzVector LV2, TMatrixD Cov2,
		TLorentzVector LV3, TMatrixD Cov3,
		double Mass1, double Mass2);
		std::vector<TLorentzVector> GetFittedLV(){
			auto HLVs = ThisFitter()->GetFittedLV();
			LV1Cor = TLorentzVector(-HLVs.at(0).Px(), HLVs.at(0).Pz(), HLVs.at(0).Py(), HLVs.at(0).E());
			LV2Cor = TLorentzVector(-HLVs.at(1).Px(), HLVs.at(1).Pz(), HLVs.at(1).Py(), HLVs.at(1).E());
			LV3Cor = TLorentzVector(-HLVs.at(2).Px(), HLVs.at(2).Pz(), HLVs.at(2).Py(), HLVs.at(2).E());
			LV4Cor = TLorentzVector(-HLVs.at(3).Px(), HLVs.at(3).Pz(), HLVs.at(3).Py(), HLVs.at(3).E());
			LV5Cor = TLorentzVector(-HLVs.at(4).Px(), HLVs.at(4).Pz(), HLVs.at(4).Py(), HLVs.at(4).E());
			return std::vector<TLorentzVector>{LV1Cor, LV2Cor, LV3Cor, LV4Cor,LV5Cor};
		}
		~KEKCascadeFitter(){
		}
};
class KEKCascadeFitter2: virtual public KEKKinematicFitter{
	private:
		TLorentzVector LV1Cor;
		TLorentzVector LV2Cor;
		TLorentzVector LV3Cor;//Decay products, eg: p,pi_L, pi_Xi
		TLorentzVector LV4Cor;//Intermediate products, eg: Lambda
		TLorentzVector LV5Cor;//Initial state, eg: Xi
	public:
		KEKCascadeFitter2(){}
		CascadeFitter2* ThisFitter(){
			return dynamic_cast<CascadeFitter2*>(Fitter);
		}
		KEKCascadeFitter2(
		TLorentzVector LV1, TMatrixD Cov1,
		TLorentzVector LV2, TMatrixD Cov2,
		TLorentzVector LV3, TMatrixD Cov3,
		double Mass1, double Mass2);
		std::vector<TLorentzVector> GetFittedLV(){
			auto HLVs = ThisFitter()->GetFittedLV();
			LV1Cor = TLorentzVector(-HLVs.at(0).Px(), HLVs.at(0).Pz(), HLVs.at(0).Py(), HLVs.at(0).E());
			LV2Cor = TLorentzVector(-HLVs.at(1).Px(), HLVs.at(1).Pz(), HLVs.at(1).Py(), HLVs.at(1).E());
			LV3Cor = TLorentzVector(-HLVs.at(2).Px(), HLVs.at(2).Pz(), HLVs.at(2).Py(), HLVs.at(2).E());
			LV4Cor = TLorentzVector(-HLVs.at(3).Px(), HLVs.at(3).Pz(), HLVs.at(3).Py(), HLVs.at(3).E());
			LV5Cor = TLorentzVector(-HLVs.at(4).Px(), HLVs.at(4).Pz(), HLVs.at(4).Py(), HLVs.at(4).E());
			return std::vector<TLorentzVector>{LV1Cor, LV2Cor, LV3Cor, LV4Cor,LV5Cor};
		}
		~KEKCascadeFitter2(){
		}
};
class KEKCascadeVertexFitter: virtual public KEKKinematicFitter{
	private:
		TLorentzVector LV1Cor;
		TLorentzVector LV2Cor;
		TLorentzVector LV3Cor;//Decay products, eg: p,pi_L, pi_Xi
		TLorentzVector LV4Cor;//Intermediate products, eg: Lambda
		TLorentzVector LV5Cor;//Initial state, eg: Xi
		TVector3 Vtx1Cor;
		TVector3 Vtx2Cor;
		TVector3 Vtx3Cor;
		TVector3 Vtx4Cor;
		TVector3 Vtx5Cor;

	public:
		KEKCascadeVertexFitter(){}
		CascadeVertexFitter* ThisFitter(){
			return dynamic_cast<CascadeVertexFitter*>(Fitter);
		}
		KEKCascadeVertexFitter	(
		TLorentzVector LV1, TVector3 Vtx1, TMatrixD Cov1,
		TLorentzVector LV2, TVector3 Vtx2, TMatrixD Cov2,
		TLorentzVector LV3, TVector3 Vtx3, TMatrixD Cov3,
		double Mass1, double Mass2);
		std::vector<TLorentzVector> GetFittedLV(){
			auto HLVs = ThisFitter()->GetFittedLV();
			LV1Cor = TLorentzVector(-HLVs.at(0).Px(), HLVs.at(0).Pz(), HLVs.at(0).Py(), HLVs.at(0).E());
			LV2Cor = TLorentzVector(-HLVs.at(1).Px(), HLVs.at(1).Pz(), HLVs.at(1).Py(), HLVs.at(1).E());
			LV3Cor = TLorentzVector(-HLVs.at(2).Px(), HLVs.at(2).Pz(), HLVs.at(2).Py(), HLVs.at(2).E());
			LV4Cor = TLorentzVector(-HLVs.at(3).Px(), HLVs.at(3).Pz(), HLVs.at(3).Py(), HLVs.at(3).E());
			LV5Cor = TLorentzVector(-HLVs.at(4).Px(), HLVs.at(4).Pz(), HLVs.at(4).Py(), HLVs.at(4).E());
			return std::vector<TLorentzVector>{LV1Cor, LV2Cor, LV3Cor, LV4Cor,LV5Cor};
		}
		std::vector<TVector3> GetFittedVerticies(){
			auto Verts = ThisFitter()->GetFittedVerticies();
			Vtx1Cor = TVector3(-Verts.at(0).X(), Verts.at(0).Z(), Verts.at(0).Y());
			Vtx2Cor = TVector3(-Verts.at(1).X(), Verts.at(1).Z(), Verts.at(1).Y());
			Vtx3Cor = TVector3(-Verts.at(2).X(), Verts.at(2).Z(), Verts.at(2).Y());
			Vtx4Cor = TVector3(-Verts.at(3).X(), Verts.at(3).Z(), Verts.at(3).Y());
			Vtx5Cor = TVector3(-Verts.at(4).X(), Verts.at(4).Z(), Verts.at(4).Y());
			return std::vector<TVector3>{Vtx1Cor, Vtx2Cor, Vtx3Cor, Vtx4Cor, Vtx5Cor};
		}
		~KEKCascadeVertexFitter(){
		}
		TVector3 GetVX(){
			TVector3 VX = TVector3(-ThisFitter()->GetVX().X(), ThisFitter()->GetVX().Z(), ThisFitter()->GetVX().Y());
			return VX;
		}
		double Gett(){
			return ThisFitter()->Gett();
		}
		double Gettcor(){
			return ThisFitter()->Gettcor();
		}
};
class KEKMassVertexFitter: virtual public KEKKinematicFitter{
	private:
		TLorentzVector LV1Cor;
		TLorentzVector LV2Cor;
		TLorentzVector LV3Cor;
	public:
		KEKMassVertexFitter(){}
		MassVertexFitter* ThisFitter(){//Necessary to call functions in FourVectorFitter. ex: SetInvMass.
			return dynamic_cast<MassVertexFitter*>(Fitter);
		}
		KEKMassVertexFitter(
		TLorentzVector LV1, TVector3 Vtx1, TMatrixD Cov1,
		TLorentzVector LV2, TVector3 Vtx2, TMatrixD Cov2, double Mass);
		std::vector<TLorentzVector> GetFittedLV(){
			auto HLVs = ThisFitter()->GetFittedLV();
			auto HLV1Cor = HLVs.at(0);
			auto HLV2Cor = HLVs.at(1);
			auto HLV3Cor = HLVs.at(2);
			LV1Cor = TLorentzVector(-HLV1Cor.Px(), HLV1Cor.Pz(), HLV1Cor.Py(), HLV1Cor.E());
			LV2Cor = TLorentzVector(-HLVs.at(1).Px(), HLVs.at(1).Pz(), HLVs.at(1).Py(), HLVs.at(1).E());
			LV3Cor = TLorentzVector(-HLVs.at(2).Px(), HLVs.at(2).Pz(), HLVs.at(2).Py(), HLVs.at(2).E());
			return std::vector<TLorentzVector>{LV1Cor, LV2Cor, LV3Cor};
		}
		std::vector<TVector3> GetFittedVerticies(){
			auto Verts = ThisFitter()->GetFittedVerticies();
			auto VPCor = Verts.at(0);
			auto VQCor = Verts.at(1);
			auto VRCor = Verts.at(2);
			return std::vector<TVector3>{VPCor, VQCor, VRCor};
		}
		~KEKMassVertexFitter(){
		}
};
class KEKMassVertexFitter3: virtual public KEKKinematicFitter{
	private:
		TLorentzVector LV1Cor;
		TLorentzVector LV2Cor;
		TLorentzVector LV3Cor;
	public:
		KEKMassVertexFitter3(){}
		MassVertexFitter3* ThisFitter(){//Necessary to call functions in FourVectorFitter. ex: SetInvMass.
			return dynamic_cast<MassVertexFitter3*>(Fitter);
		}
		KEKMassVertexFitter3(
		TLorentzVector LV1, TVector3 Vtx1, TMatrixD Cov1,
		TLorentzVector LV2, TVector3 Vtx2, TMatrixD Cov2, double Mass);
		std::vector<TLorentzVector> GetFittedLV(){
			auto HLVs = ThisFitter()->GetFittedLV();
			auto HLV1Cor = HLVs.at(0);
			auto HLV2Cor = HLVs.at(1);
			auto HLV3Cor = HLVs.at(2);
			LV1Cor = TLorentzVector(-HLV1Cor.Px(), HLV1Cor.Pz(), HLV1Cor.Py(), HLV1Cor.E());
			LV2Cor = TLorentzVector(-HLVs.at(1).Px(), HLVs.at(1).Pz(), HLVs.at(1).Py(), HLVs.at(1).E());
			LV3Cor = TLorentzVector(-HLVs.at(2).Px(), HLVs.at(2).Pz(), HLVs.at(2).Py(), HLVs.at(2).E());
			return std::vector<TLorentzVector>{LV1Cor, LV2Cor, LV3Cor};
		}
		std::vector<TVector3> GetFittedVerticies(){
			auto Verts = ThisFitter()->GetFittedVerticies();
			auto VPCor = Verts.at(0);
			auto VQCor = Verts.at(1);
			auto VRCor = Verts.at(2);
			return std::vector<TVector3>{VPCor, VQCor, VRCor};
		}
		~KEKMassVertexFitter3(){
		}
};
#endif