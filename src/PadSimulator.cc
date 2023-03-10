// -*- C++ -*-

#include "PadSimulator.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <TH2D.h>

#include "ConfMan.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DCTrackSearch.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FiberCluster.hh"
#include "FuncName.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "K18Parameters.hh"
//#include "K18TrackU2D.hh"
#include "K18TrackD2U.hh"
#include "KuramaTrack.hh"
#include "MathTools.hh"
#include "MWPCCluster.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCTrackSearch.hh"

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

// Tracking routine selection __________________________________________________
/* TPCTracking */
#define UseTpcCluster     1

namespace
{
	using namespace K18Parameter;
	const auto& gConf   = ConfMan::GetInstance();
	const auto& gGeom   = DCGeomMan::GetInstance();
	const auto& gTPC  = TPCParamMan::GetInstance();
	const auto& gTPCPos = TPCPositionCorrector::GetInstance();
	const auto& gUser   = UserParamMan::GetInstance();

	//_____________________________________________________________________________
	const Double_t& pK18 = ConfMan::Get<Double_t>("PK18");


	const Int_t    MaxNumOfTrackTPC = 20;
	const Int_t    MaxRowDifTPC = 2; // for cluster

	//_____________________________________________________________________________

	//_____________________________________________________________________________
	inline void
		printConnectionFlag(const std::vector<std::deque<Bool_t> >& flag)
		{
			for(Int_t i=0, n=flag.size(); i<n; ++i){
				hddaq::cout << std::endl;
				for(Int_t j=0, m=flag[i].size(); j<m; ++j){
					hddaq::cout << " " << flag[i][j];
				}
			}
			hddaq::cout << std::endl;
		}
}

//_____________________________________________________________________________
PadSimulator::PadSimulator()
	: m_TPCHitCont(NumOfLayersTPC),
	//	m_TempTPCHitCont(NumOfLayersTPC),
	m_TPCClCont(NumOfLayersTPC)
{

	zK18HS=gGeom.LocalZ("K18HS");
	TPC_Dt=gUser.GetParameter("TPC_Dt");
	debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
PadSimulator::~PadSimulator()
{
	delete DriftLength;
	ClearTracksTPC();
	ClearTPCClusters();
	ClearTPCHits();
	debug::ObjectCounter::decrease(ClassName());
}
SimTrack::~SimTrack()
{
	debug::ObjectCounter::decrease(ClassName());
}
	void
PadSimulator::ClearTracksTPC()
{
	del::ClearContainer(m_TPCTC);
	del::ClearContainer(m_TPCTC_Helix);
}
//_____________________________________________________________________________

//_____________________________________________________________________________
	Bool_t
PadSimulator::TrackSearchTPC()
{
	static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

	tpc::LocalTrackSearch(m_TPCClCont, m_TPCTC, m_TPCTC_Failed, MinLayer);
	return true;
}

//_____________________________________________________________________________
	Bool_t
PadSimulator::TrackSearchTPCHelix()
{
	static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

	tpc::HelixTrackSearch(0, 0, m_TPCClCont, m_TPCTC_Helix, m_TPCTC_HelixFailed, MinLayer);
	return true;
}

//_____________________________________________________________________________
	void
PadSimulator::ClearTPCHits()
{
	del::ClearContainerAll(m_TPCHitCont);
	//	del::ClearContainerAll(m_TempTPCHitCont);
}

//_____________________________________________________________________________
	void
PadSimulator::ClearTPCClusters()
{
	del::ClearContainerAll(m_TPCClCont);
}


	Bool_t
PadSimulator::Simulate(std::vector<Double_t> chisqr,
		std::vector<Double_t> x0,
		std::vector<Double_t> y0,
		std::vector<Double_t> u0,
		std::vector<Double_t> v0,
		TH2Poly* padHist)
{
	double param[2] = {zK18HS,TPC_Dt};
	int ntrk = chisqr.size();
	SimTrack Trk[20];
	if(ntrk==0){
		return false;
	}
	for(int layer=0;layer<32;++layer){
		if(m_TPCHitCont[layer].size()>0) std::cout<<"Hit?"<<std::endl;
		m_TPCHitCont[layer].clear();
		m_TPCClCont[layer].clear();
	}
	static const Double_t MaxYDif = gUser.GetParameter("MaxYDifClusterTPC");

	for(int i=0;i<1;++i){
		Trk[i]=SimTrack(chisqr[0],x0[0],y0[0],u0[0],v0[0]);
		for(int j=0;j<np;++j){
			double t = 500.*(j+0.5)/np-250;
			TVector3 pos = Trk[i].GetDiffusedPosition(t,param);
			double x =pos.X(),y=pos.Y(),z=pos.Z();
			int padID = tpc::findPadID(z,x);
			if(padID<1 or padID>max_pad) continue;
			if(tpc::Dead(padID)) continue;
			padHist ->Fill(z,x);

			if(padID>0)DriftLength->at(padID-1)+=y;
		}
	}
	for(int ipad = 1;ipad<max_pad+1;++ipad){
		if(tpc::Dead(ipad)) continue;
		double de = padHist->GetBinContent(ipad);
		if(!de) continue;
		TVector3 pos = tpc::getPosition(ipad);
		double y = DriftLength->at(ipad-1)/de;
		pos.SetY(y);
		Int_t layer = tpc::getLayerID(ipad);
		Int_t row = tpc::getRowID(ipad);
#if 1
		auto hit = new TPCHit(layer,row);
		hit->AddHit(TMath::QuietNaN(),TMath::QuietNaN());
		hit->SetPad(ipad);
		hit->SetLayer(layer);
		hit->SetRow(row);
		hit->SetDe(de);
		hit->SetPosition(pos);
		if(layer>31) std::cout<<"Warning , layer = "<<layer<<std::endl;
		m_TPCHitCont[layer].push_back(hit);
#endif

	}
	for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
		//    if(layer==ExlayerID) continue; //exclusive layer
		MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], MaxYDif);
	}
	return true;
}


	Bool_t
PadSimulator::MakeUpTPCClusters(const TPCHitContainer& HitCont,
		TPCClusterContainer& ClCont,
		Double_t maxdy, Double_t de_cut)
{
	//static const Double_t ClusterDeCut = gUser.GetParameter("MinClusterDeTPC");
	static const Double_t ClusterDeCut = de_cut;
	const auto nh = HitCont.size();
	if(nh==0) return false;

	std::vector<Int_t> joined(nh, 0);
	for(Int_t i=0; i<nh; ++i){
		if(joined[i] > 0) continue;
		TPCHitContainer CandCont;
		TPCHit* hit = HitCont[i];
		if(!hit) continue;
		Int_t layer = hit->GetLayer();
		CandCont.push_back(hit);
		joined[i]++;
		for(Int_t j=0; j<nh; ++j){
			if(i==j || joined[j]>0) continue;
			TPCHit* thit = HitCont[j];
			if(!thit ) continue;
			Int_t rowID = thit->GetRow();
			for(const auto& c_hit: CandCont){
				Int_t c_rowID = c_hit->GetRow();
				if(tpc::IsClusterable(layer, rowID, c_rowID, MaxRowDifTPC)
						&& TMath::Abs(thit->GetY() - c_hit->GetY()) < maxdy){
					CandCont.push_back(thit);
					joined[j]++;
					break;
				}
			}
		}
		TPCCluster* cluster = new TPCCluster(layer, CandCont);
		if(!cluster) continue;
		if(cluster->Calculate() && cluster->GetDe()>ClusterDeCut){
			ClCont.push_back(cluster);
		}else{
			delete cluster;
		}
	}

	return true;
}
TVector2
SimTrack::RandomGaus(double sig){
	double x = gRandom->Gaus(0,sig);
	double y = gRandom->Gaus(0,sig);
	return TVector2(x,y);
}
