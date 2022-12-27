#ifndef TPC_Beam_Remover_HH
#define TPC_Beam_Remover_HH
#include <TSpectrum.h>
#include <vector>
#include <TString.h>
#include <TVector3.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include "DCAnalyzer.hh"
class TPCCluster;
class TPCHit;

class TPCBeamRemover{
	public:
		static TString ClassName();
		TPCBeamRemover(std::vector<TPCClusterContainer>ClCont);
		~TPCBeamRemover();	
	private:
		std::vector<TPCClusterContainer> m_ClCont_array;
		TPCClusterContainer m_Cl_array;
		std::vector<TPCClusterContainer> m_BeamClCont_array;
		TPCClusterContainer m_BeamCl_array;
		std::vector<int> m_layer_info;
		std::vector<bool> m_Accidental;
		std::vector<double> beam_p0;
		std::vector<double> beam_p1;
		std::vector<double> beam_p2;
		std::vector<double> beam_y;
		std::vector<double> beam_v;
		double y_min,y_max;
		double Ywidth; 
		double Ycut;
		double Xcut;
		bool enable;
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
		TH1D* hist_y = new TH1D("hist_y","hist_y",140,-350,350);
		TH2D* hist_beam = new TH2D("hist_beam","hist_beam",100,-250,250,120,-120,120);
		TH2D* hist_beamZY = new TH2D("hist_beamZY","hist_beamZY",100,-250,250,700,-350,350);
	public:
		std::vector<double>GetParameter(int i){
			if(i==0){
				return beam_y;
			}else if(i==1){
				return beam_p0;
			}else if(i==2){
				return beam_p1;
			}else if(i==3){
				return beam_p2;
			}else if (i==4){
				return beam_v;
			}
		}
		std::vector<TPCClusterContainer> GetClCont(){
			return m_ClCont_array;
		};
		std::vector<TPCClusterContainer> GetBeamClCont(){
			return m_BeamClCont_array;//ClContainor of accidental beams.
		};

		int SearchPeaks(TH1D* hist,std::vector<double> &peaks);
		void 	SearchAccidentalBeam(double ymin,double ymax);//(y_min,y_max) = exception region
		void 	SortAccidentalBeam();
		bool IsAccidental(TPCHit* hit);
		bool IsAccidental(TVector3 pos);
		void Enable(bool status){
			enable = status;
		}
};

	inline TString
TPCBeamRemover::ClassName()
{
	static TString s_name("TPCBeamRemover");
	return s_name;
}







#endif
