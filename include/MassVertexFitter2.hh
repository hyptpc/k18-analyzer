#ifndef MassVertexFitter2_h
#define MassVertexFitter2_h
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
// The coordinate to work on will be helix coordinate.
// Note that the covariance sign of phi should be the opposite.
#include "KinFit.hh"
#include <TVector3.h>
#include <TLorentzVector.h>
class MassVertexFitter2: virtual public KinematicFitter{
	//R -> P + Q;
	protected:	
		TLorentzVector P;
		TVector3 VP;//put vertex here.
		//Note that, the 'position covariance'
		//is defined for helix parameters, cx,cy,z0.
		TLorentzVector PCor;
		TVector3 VPCor;
		double mP;

		TLorentzVector Q;
		TVector3 VQ;
		TLorentzVector QCor;
		TVector3 VQCor;
		double mQ;

		TLorentzVector R;
		TVector3 VR;
		TLorentzVector RCor;
		TVector3 VRCor;
		vector<double> MassDiffs;
		double mR;

		bool MeasDir = false;
		double RadToMom = 0;
		int charge_config = 0;//0: particle 1 trajectory is on the positive side of the particle 2 trajectory, 1: opposite. This is for the definition of R in cylindrical coordinate.

		TMatrixD VtxCov;//Since the resulting vertex, i.e. V_Lambda, is not a direct parameter of the fit, its covariance is not directly calculated. We need an explicit calculation using error propagation. This matrix is for that purpose. It is a 3x3 matrix for the x,y,z coordinates of the vertex.
		vector<TMatrixD> StepVtxCov;//Since the resulting vertex, i.e. V_Lambda, is not a direct parameter of the fit, its covariance is not directly calculated. We need an explicit calculation using error propagation. This matrix is for that purpose. It is a 3x3 matrix for the x,y,z coordinates of the vertex.
	public:
		MassVertexFitter2(){}
		MassVertexFitter2(TLorentzVector P_,TLorentzVector Q_,TLorentzVector R_
			,TVector3 V_P, TVector3 V_Q);
		void SetInvMass(double IM){
			mR = IM;
		}
		vector<TLorentzVector> GetFittedLV(){
			vector<TLorentzVector> ret = {PCor,QCor,RCor};
			return ret;
		}
		vector<TVector3> GetFittedVerticies(){
			vector<TVector3> ret = {VPCor,VQCor,0.5*(VPCor+VQCor)};
			return ret;
		}
		void ToDecayPlane();
		void SetBField(double B){
			Double_t ConstC = 0.299792458; //=c/10^9
			RadToMom = ConstC * B * 1e-3; //1e-3: mm -> m
		}
		double CalcHelixDistance2(std::vector<double> pars);
		double DerivativeHelixDistance2(int index, std::vector<double> pars);
		void GetHelixParameters(double p, double th, double ph, TVector3 center, int sign, double* pars);
		TVector3 CalcHelixPosition(double p, double th, double ph, TVector3 center, int sign);
		TVector3 DerivativeHelixPosition(int index, double p, double th, double ph, TVector3 center, int sign);
	protected:
		virtual void Initialize();
		virtual void SampleStepPoint(int steps);
		virtual void SetConstraints();
		virtual void Rotate();
		TMatrixD JacobianSphToCart(double p, double th, double ph);		
		virtual void CalcVariance(int istep);
};
#endif
