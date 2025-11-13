#ifndef MassVertexFitter_h
#define MassVertexFitter_h
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
#include "KinFit.hh"
#include <TVector3.h>
#include <TLorentzVector.h>
class MassVertexFitter: virtual public KinematicFitter{
	//R -> P + Q;
	protected:	
		TLorentzVector P;
		TVector3 VP;
		TVector3 Pres;
		TLorentzVector PCor;
		TVector3 VPCor;
		double mP;

		TLorentzVector Q;
		TVector3 VQ;
		TVector3 Qres;
		TLorentzVector QCor;
		TVector3 VQCor;
		double mQ;

		TLorentzVector R;
		TVector3 VR;
		TVector3 Rres;
		TLorentzVector RCor;
		vector<double> MassDiffs;
		double mR;

		bool MeasDir = false;
		TVector3 V0;
	public:
		MassVertexFitter(){}
		MassVertexFitter(TLorentzVector P_,TLorentzVector Q_,TLorentzVector R_
			,TVector3 V_P, TVector3 V_Q);
		void SetInvMass(double IM){
			mR = IM;
		}
		vector<TLorentzVector> GetFittedLV(){
			vector<TLorentzVector> ret = {PCor,QCor,RCor};
			return ret;
		}
		vector<TVector3> GetFittedVerticies(){
			vector<TVector3> ret = {VPCor+V0,VQCor+V0,0.5*(VPCor+VQCor)+V0};
			return ret;
		}
		void GetRZparameters(TVector3 Vert, TVector3 Dir, double &R, double &Z);
		double CalcVertexDistance(double th1, double ph1, double r1, double z1,
			double th2, double ph2, double r2, double z2);
		void CalcClosePoint(double th1, double ph1, double r1, double z1,
			double th2, double ph2, double r2, double z2, TVector3& vert1, TVector3& vert2);//Assigns the closest point of the particle 1 trajectory and particle 2 trajectory to each other.
		void GetVertexResolution(double& dr, double& dz );
		void ToDecayPlane();
	protected:
		virtual void Initialize();
		virtual void SampleStepPoint(int steps);
		virtual void SetConstraints();
		virtual void Rotate();
		TMatrixD JacobianSphToCart(double p, double th, double ph);		
		virtual void CalcVariance(int istep);
};
#endif
