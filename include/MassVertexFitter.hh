#ifndef MassVertexFitter_h
#define MassVertexFitter_h
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
#include "KinFit.hh"
#include <TVector3.h>
#include <TLorentzVector.h>
class MassVertexFitter: virtual public KinematicFitter{
	// P: proton, Q: pion, R: pion, L: Lambda, X: Xi
	//L -> P + Q;
	//X -> L + R;
	protected:	
		/*
		Note that, the relationship between the kinematic parameter 
		and vertex position should be defined in the off-diagonal term
		in the covariance matrix. These off-diagonal terms should 'drag'
		the vertex position according to the change of the kinematic parameters.
		Without these terms, the fitter will simply return the 1-C fit result,
		with P Q vertices just moved to the averaged vertex.
		*/
		TLorentzVector P;
		TVector3 Pres;
		TLorentzVector PCor;
		TVector3 VP;
		TVector3 VPCor;
		double mP;

		TLorentzVector Q;
		TVector3 Qres;
		TLorentzVector QCor;
		TVector3 VQ;
		TVector3 VQCor;
		double mQ;
		
		TLorentzVector L;
		TVector3 Lres;
		TLorentzVector LCor;
		TVector3 VL;
		TVector3 VLCor;
		vector<double> MassDiffsL;
		double mL;
		
		double t;//Lambda extrapolation parameter. left for future...
		double t_cor;

		bool UseVertexFlag = 0;
	public:
		MassVertexFitter(){}
		MassVertexFitter(TLorentzVector P_,TVector3 V_P,
			TLorentzVector Q_,TVector3 V_Q);
		void SetInvMass(double ML){
			mL = ML;
		}
		vector<TLorentzVector> GetFittedLV(){
			vector<TLorentzVector> ret = {PCor,QCor,LCor};	
			return ret;
		}
		vector<TVector3> GetFittedVerticies(){
			vector<TVector3> ret = {VPCor,VQCor,VLCor};
			return ret;
		}
		void UseVertex(bool status,TVector3 Vert1,TVector3 Vert2);
		void ToDecayPlane();
		double CalcLambdaExtrapolationParameter(vector<double> Meas_, vector<double> Unkn_);
		double CalcLambdaExtrapolationParameterDerivativesM(int idx, vector<double> Meas_, vector<double> Unkn_);
		double CalcLambdaExtrapolationParameterDerivativesU(int idx, vector<double> Meas_, vector<double> Unkn_);
		double Calc_dtdM(int idx, vector<double> Meas_, vector<double> Unkn_){
			return CalcLambdaExtrapolationParameterDerivativesM(idx,Meas_,Unkn_);
		};
		double Calc_dtdU(int idx, vector<double> Meas_, vector<double> Unkn_){
			return CalcLambdaExtrapolationParameterDerivativesU(idx,Meas_,Unkn_);
		};
		TVector3 CalcLambdaDirectionalDerivativesM(int idx, vector<double> Meas_);
		TVector3 Calc_dDdM(int idx, vector<double> Meas_){
			return CalcLambdaDirectionalDerivativesM(idx,Meas_);
		};
		TVector3 CalcLambdaE1Vector(vector<double> Meas_);
		TVector3 CalcLambdaE1VectorDerivativesM(int idx, vector<double> Meas_);
		TVector3 Calc_dE1dM(int idx, vector<double> Meas_){
			return CalcLambdaE1VectorDerivativesM(idx,Meas_);
		};
		TVector3 CalcLambdaE2Vector(vector<double> Meas_);
		TVector3 CalcLambdaE2VectorDerivativesM(int idx, vector<double> Meas_);
		TVector3 Calc_dE2dM(int idx, vector<double> Meas_){
			return CalcLambdaE2VectorDerivativesM(idx,Meas_);
		};
		TVector3 CalcV_LX(vector<double> Meas_);
		TVector3 CalcV_LXDerivativesM(int idx, vector<double> Meas_);
		TVector3 Calc_dV_LXdM(int idx, vector<double> Meas_){
			return CalcV_LXDerivativesM(idx,Meas_);
		};

		double Gett(){
			return t;
		}
		double Gettcor(){
			return t_cor;
		}
	protected:
		virtual void Initialize();
		virtual void SampleStepPoint(int steps);
		virtual void SetConstraints();
		virtual void Rotate();
		TMatrixD JacobianSphToCart(double p, double th, double ph);		
		virtual void CalcVariance(int istep);
};
#endif