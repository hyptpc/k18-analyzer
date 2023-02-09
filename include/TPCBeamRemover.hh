#ifndef TPCBeam_Remover_HH
#define TPCBeam_Remover_HH
#include <TSpectrum.h>
#include <vector>
#include <TString.h>
#include <TVector3.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include "DCAnalyzer.hh"
class TPCCluster;
class TPCHit;
class TPCBeamRemover{
	public:
		static TString ClassName();
		TPCBeamRemover(std::vector<TPCClusterContainer>ClCont);
		~TPCBeamRemover();	
	private:
		int MaxNBeam = 3;
		int MinBeamHit = 10;
		double MinZCut = -180;
		double MaxZCut = 150;
		double MaxHoughWindow = 15;
		double MaxHoughWindowY = 5;

		const Int_t nBin_rdiff = 110;
 		const Double_t rdiff_min = -110.;
  	const Double_t rdiff_max = 110.;
  	const Int_t nBin_theta = 180;
  	const Double_t theta_min = (1.5)*acos(-1);
  	const Double_t theta_max = (2.5)*acos(-1);//Charge < -1 region.
  	const Int_t nBin_p = 100;
  	const Double_t p_min = 1600.;//MeV/c
  	const Double_t p_max = 2000.;//MeV/c
  	const int    thetaY_ndiv =  180;
  	const double thetaY_min  =  60.;
  	const double thetaY_max  = 120.;
  	const int    r_ndiv =  2000;
  	const double r_min  = -5000.;
  	const double r_max  =  5000.;




	private:
		std::vector<TPCClusterContainer> m_ClCont_array;
		TPCClusterContainer m_Cl_array;
		std::vector<TPCClusterContainer> m_PeakCl_array;
		std::vector<TPCClusterContainer> m_BeamClCont_array;
		TPCClusterContainer m_BeamCl_array;
		std::vector<double>peaks;
		std::vector<int> m_layer_info;
		std::vector<std::vector<int>> m_Peaklayer;
		std::vector<bool> m_Accidental;
		std::vector<double> beam_p0;
		std::vector<double> beam_p1;
		std::vector<double> beam_p2;
		std::vector<double> beam_y;
		std::vector<double> beam_v;
		double x_min,x_max,y_min,y_max;
		double Xcut,Ycut;
		double Ywidth; 
		bool enable;
		bool enableHough;
		int np,nh;

		std::vector<std::vector<double>> h_cx;
		std::vector<std::vector<double>> h_cy;
		std::vector<std::vector<double>> h_z0;
		std::vector<std::vector<double>> h_r;
		std::vector<std::vector<double>> h_dz;
		std::vector<std::vector<double>> acc_tparam;
		std::vector<std::vector<int>> h_flag;
		std::vector<double> Allh_cx;
		std::vector<double> Allh_cy;
		std::vector<double> Allh_z0;
		std::vector<double> Allh_r;
		std::vector<double> Allh_dz;
		
		
		
		
		double Const = 0.299792458; // =c/10^9
		double HS_field_0 = 0.9860;
		double HS_field_Hall_calc;
		double HS_field_Hall;
		double dMagneticField;

		double quadratic(double z,double p0,double p1,double p2){
			return p0+p1*z+p2*z*z;
		}
		double linear(double z,double p0,double p1){
			return p0+p1*z;
		}
		TF1* Quadratic = new TF1("quadratic","pol2",-250,250);
		TF1* Linear = new TF1("linear","pol1",-250,250);
		TF1* Arc = new TF1("Arc","-TMath::Sqrt([2]*[2]-(x-[1])*(x-[1]))+[0]",-150,400);
		
		TH1D* hist_y = new TH1D("hist_y","hist_y",140,-350,350);
		TH2D* hist_beam = new TH2D("hist_beam","hist_beam",100,-250,250,120,-120,120);
		TH2D* hist_beamZY = new TH2D("hist_beamZY","hist_beamZY",100,-250,250,700,-350,350);
		
		TGraph* ArcGraph;
		TGraph* YThetaGraph;
		TH3D* hist_Ci;
		TH2D* hist_YTheta;
		
	private:
		int SearchPeaks(TH1D* hist,std::vector<double> &peaks);
		
		void DoQuadraticSearch(int i);
		void DoHoughSearch(int i);
		
		void DoCircleHough(int i);
		void DoYThetaHough(int i);
		void DoYThetaFit(int i);
		void CircleFit(std::vector<TVector3> pos,double* param);
		void LinearFit(std::vector<TVector3> pos,double* param);

		int CompareHough(TVector3 pos, std::vector<double> hcx,std::vector<double>hcy,std::vector<double> hr); 
		bool IsThisBeam(int hflag, int ib);	
	
	public:
		std::vector<double>GetParameter(int i){
			if(enableHough){
				switch(i){
					case 0:
						return Allh_cx; 
						break;
					case 1:
						return Allh_cy; 
						break;
					case 2:
						return Allh_z0; 
						break;
					case 3:
						return Allh_r; 
						break;
					case 4:
						return Allh_dz; 
						break;
				}
			}
			else{
				switch(i){
					case 0:
						return beam_p0;
						break;
					case 1:
						return beam_p1;
						break;
					case 2:
						return beam_p2;
						break;
					case 3:
						return beam_y;
						break;
					case 4:
						return beam_v;
						break;
				}
			}
			std::vector<double>Dummy;
			return Dummy;
		}
		std::vector<TPCClusterContainer> GetClCont(){
			return m_ClCont_array;
		};
		std::vector<TPCClusterContainer> GetBeamClCont(){
			return m_BeamClCont_array;//ClContainor of accidental beams.
		};

		void 	SearchAccidentalBeam(double xmin,double xmax,double ymin,double ymax);//(y_min,y_max) = exception region
		int IsAccidental(TPCHit* hit);
		int IsAccidental(TVector3 pos, TVector3 res, double& Adist);
		void Enable(bool status){
			enable = status;
		}
		void EnableHough(bool status){
			enableHough = status;
		}
};

inline TString
TPCBeamRemover::ClassName()
{
	static TString s_name("TPCBeamRemover");
	return s_name;
};







#endif
