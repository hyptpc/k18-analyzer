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
#include "TPCLocalTrackHelix.hh"
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
		double MaxHoughWindowY = 2.5;




	private:
		std::vector<TPCClusterContainer> m_ClCont_array;
		std::vector<TPCLocalTrackHelix*> m_Track_array;
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
		double Ywindow;//=10 
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
		std::vector<double> AllPeaks;
		
		
		
		int Acc_flag_base = 400;	
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
		
		
	private:
		int SearchPeaks(TH1D* hist,std::vector<double> &peaks);
		
		void DoHoughSearch(int i);
			
		void DoCircleHough(int i);
		void DoYThetaHough(int i);
		void DoYThetaFit(int i);

		void DoHelixFit(TPCLocalTrackHelix* Track,const std::vector<TPCClusterContainer>& ClCont,int MinNumOfHits);

		int CompareHoughDist(TVector3 pos, std::vector<double> hcx,std::vector<double>hcy,std::vector<double> hr); 
		bool IsThisBeam(int hflag, int ib);	
	
	public:
		std::vector<double>GetParameter(int i){
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
			std::vector<double>Dummy;
			return Dummy;
		}
		int GetNAccBeam(){
			return Allh_cx.size();
		}
		std::vector<TPCClusterContainer> GetClCont(){
			return m_ClCont_array;
		};
		std::vector<TPCClusterContainer> GetBeamClCont(){
			return m_BeamClCont_array;//ClContainor of accidental beams.
		};
		TPCLocalTrackHelix* GetAccTrack(int iacc){
			return m_Track_array.at(iacc);
		}


		void 	SearchAccidentalBeam(double xmin,double xmax,double ymin,double ymax);//(y_min,y_max) = exception region
		void ConstructAccidentalTracks();	
		
		int IsAccidental(TPCHit* hit);
		int IsAccidental(TVector3 pos, TVector3 res, double& PullDist);
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
