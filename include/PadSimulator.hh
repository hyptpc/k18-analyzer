// -*- C++ -*-

#ifndef PADSIMULATOR_HH
#define PADSIMULATOR_HH

#include <vector>
#include <TString.h>

#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCPadHelper.hh"
#include "TRandom.h"
class RawData;

class TPCHit;
class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrackHelix;


typedef std::vector<TPCHit*>        TPCHitContainer;
typedef std::vector<TPCCluster*>    TPCClusterContainer;
typedef std::vector<TPCLocalTrack*> TPCLocalTrackContainer;
typedef std::vector<TPCLocalTrackHelix*> TPCLocalTrackHelixContainer;

static const int max_pad =  5768;
//_____________________________________________________________________________
class PadSimulator
{
public:
  static const TString& ClassName();
  PadSimulator();
  ~PadSimulator();
private:
  PadSimulator(const PadSimulator&);
  PadSimulator& operator =(const PadSimulator&);
private:
  enum e_type
  { kBcIn, kBcOut, kSdcIn, kSdcOut, kTPC, kTOF, n_type };
	std::vector<double>* DriftLength = new std::vector<double>(max_pad,0);
  std::vector<TPCHitContainer>       m_TPCHitCont;
  std::vector<TPCHitContainer>       m_TempTPCHitCont;
  std::vector<TPCClusterContainer>   m_TPCClCont;
  TPCLocalTrackContainer             m_TPCTC;
  TPCLocalTrackHelixContainer        m_TPCTC_Helix;
  TPCLocalTrackHelixContainer        m_TPCTC_HelixFailed;
	double zK18HS,TPC_Dt;
	int np;
public:


//	void InitializePadHist(int i);
  const TPCHitContainer& GetTPCHC(Int_t l) const { return m_TPCHitCont.at(l); }
  const TPCClusterContainer& GetTPCClCont(Int_t l) const
    { return m_TPCClCont.at(l); }

  Bool_t TrackSearchTPC();
  Bool_t TrackSearchTPCHelix();

  Int_t GetNTracksTPC() const { return m_TPCTC.size(); }
  Int_t GetNTracksTPCHelix() const { return m_TPCTC_Helix.size(); }
  // Exclusive Tracks

  TPCLocalTrack* GetTrackTPC(Int_t l) const { return m_TPCTC.at(l); }
  TPCLocalTrackHelix* GetTrackTPCHelix(Int_t l) const
    { return m_TPCTC_Helix.at(l); }
  // Exclusive Tracks


	Bool_t Simulate(std::vector<Double_t> chisqr,
			std::vector<Double_t> x0,
			std::vector<Double_t> y0,
			std::vector<Double_t> u0,
			std::vector<Double_t> v0,
			TH2Poly* padHist);

  void HoughYCut(Double_t min_y, Double_t max_y){};
	void SetzK18HS(double zK18HS_){zK18HS = zK18HS_;}
	void SetNumberOfPoints(int np_){np = np_;}
protected:

//  void ClearDCHits();
  void ClearTPCHits();
  void ClearTPCClusters();

  void ClearTracksTPC();
  static Bool_t MakeUpTPCClusters(const TPCHitContainer& HitCont,TPCClusterContainer& ClCont, Double_t maxdy,Double_t de_cut=0);
public:
//	void KillMe(){~PadSimulator();}
};

//_____________________________________________________________________________
inline const TString&
PadSimulator::ClassName()
{
  static TString s_name("PadSimulator");
  return s_name;
}

class SimTrack{
	private:
  	static const TString& ClassName();
		double m_chisqr,m_x0,m_y0,m_u0,m_v0;
	public:
		SimTrack(){
			debug::ObjectCounter::increase(ClassName());
		};
		SimTrack(double chisqr,double x0,double y0,double u0,double v0){
			m_chisqr = chisqr;m_x0=x0;m_y0 = y0;m_u0 = u0;m_v0 = v0;
			debug::ObjectCounter::increase(ClassName());
		}
		~SimTrack();
		TVector2 RandomGaus(double sig);
		TVector3 GetDiffusedPosition(double t, double* par){
			double z = t;
			double zK18HS = par[0];
			double TPC_Dt = par[1];
			t+=zK18HS;
			double x = m_x0+m_u0*z;
			double y = m_y0+m_v0*z;
			double sig = TMath::Sqrt(2)*TPC_Dt*TMath::Sqrt(y*0.1+30);
			TVector2 Diffusion = RandomGaus(sig);
			x+=Diffusion.X();
			z+=Diffusion.Y();


			return TVector3(x,y,z);
		}

};
inline const TString&
SimTrack::ClassName()
{
  static TString s_name("SimTrack");
  return s_name;
}

#endif
