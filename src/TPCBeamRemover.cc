#ifndef TPC_Beam_Remover_CC
#define TPC_Beam_Remover_CC
#include "TPCBeamRemover.hh"
#include "DCGeomMan.hh"
#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"
#include "TPCCluster.hh"
#include "RootHelper.hh"
namespace
{
	const auto& gConf = ConfMan::GetInstance();
	const auto& gGeom = DCGeomMan::GetInstance();
	const auto& gUser = UserParamMan::GetInstance();
	const auto& zTarget    = gGeom.LocalZ("Target");
	const auto& zK18Target = gGeom.LocalZ("K18Target");

	// TPC Tracking

}
TPCBeamRemover::TPCBeamRemover(std::vector<TPCClusterContainer>ClCont){
	m_ClCont_array.clear();
	m_Cl_array.clear();
	m_BeamClCont_array.clear();
	m_BeamCl_array.clear();
	m_layer_info.clear();
	m_Accidental.clear();
	beam_p0.clear();
	beam_p1.clear();
	beam_p2.clear();
	beam_y.clear();
	beam_v.clear();
	enable = true;
	for(int layer=0;layer<NumOfLayersTPC;layer++){
		TPCClusterContainer temp_cl_arr;
		temp_cl_arr.clear();
		for(auto hitcl:ClCont[layer]){
			temp_cl_arr.push_back(hitcl);
			m_Cl_array.push_back(hitcl);
			m_layer_info.push_back(layer);
		}
		m_ClCont_array.push_back(temp_cl_arr);
	}
	Ywidth = 10; 
	Ycut = 5; 
	Xcut = 10; 
	hist_y->Reset("ICES");
	hist_beam->Reset("ICES");
	hist_beamZY->Reset("ICES");
	HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
	HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
	dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

	debug::ObjectCounter::increase(ClassName());
}

TPCBeamRemover::~TPCBeamRemover()
{
	delete Quadratic;
	delete Linear;
	delete hist_beam;
	delete hist_beamZY;
	delete hist_y;
	debug::ObjectCounter::decrease(ClassName());
}


int 
TPCBeamRemover::SearchPeaks(TH1D* hist,std::vector<double> &peaks){
	TSpectrum spec(30);
	double sig=1,th=0.15;
	int npeaks = spec.Search(hist,sig,"goff",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	return npeaks;
}


void 
TPCBeamRemover::SearchAccidentalBeam(double ymin = -50.,double ymax = 50.){
	y_min=ymin;y_max=ymax;
	if(!enable) return;
	for(auto hitcl:m_Cl_array){
		TVector3 pos = hitcl->GetPosition();
		double y = pos.Y();
		if(y>y_min and y<y_max)continue;
		else hist_y->Fill(y);
	}
	std::vector<double>peaks;
	int np = SearchPeaks(hist_y,peaks);
	int nh = m_Cl_array.size();
	m_Accidental.resize(nh,false);
	Linear->SetParameter(0,0);
	Linear->SetParLimits(0,-400,400);

	Linear->SetParameter(1,0);
	Linear->SetParLimits(1,-0.1,0.1);

	Quadratic->SetParameter(0,0);
	Quadratic->SetParLimits(0,-120,120);

	Quadratic->SetParameter(1,-0.05);
	Quadratic->SetParLimits(1,-0.1,-0.01);

	Quadratic->SetParameter(2,-7.6e-5);
	Quadratic->SetParLimits(2,-3e-4,-3e-5);

	for(int i = 0;i<np;++i){
		double peak = peaks.at(i);
		hist_beam->Reset("ICES");
		hist_beamZY->Reset("ICES");
		for(auto hitp:m_Cl_array){
			TVector3 pos = hitp->GetPosition();
			double x=pos.X(),y=pos.Y(),z=pos.Z();
			if(y>y_min and y<y_max)continue;
			if(abs(peak-y)>Ywidth) continue;
			hist_beam->Fill(z,x);	
			hist_beamZY->Fill(z,y);	
		}
		Linear->SetParameter(0,peak);
		bool isAccidental = false;
		double min_z = 255;
		double dist_cut[3]={20,10,5};
		double dist_cutY[3]={10,7.5,5};
		for(int nitr=0;nitr<3;++nitr){
			min_z = 255;
			Quadratic->SetParameter(1,-0.05);
			if(hist_beam->GetEffectiveEntries()<5) continue;
			hist_beam->Fit("quadratic","QRNB");
			hist_beamZY->Fit("linear","QRNB");
			double p0,p1,p2,v=0;
			hist_beam->Reset("ICES");
			hist_beamZY->Reset("ICES");
			p0=Quadratic->GetParameter(0);
			p1=Quadratic->GetParameter(1);
			p2=Quadratic->GetParameter(2);
			peak = Linear->GetParameter(0);
			v = Linear->GetParameter(1);
			for(auto hitp:m_Cl_array){
				TVector3 pos = hitp->GetPosition();
				double x = pos.X(),y = pos.Y(),z = pos.Z();
				if(y>y_min and y<y_max)continue;
				if(abs(peak+v*z-y)>dist_cutY[nitr] ) continue;
				double val = quadratic(z,p0,p1,p2);
				if(abs(val-x)>dist_cut[nitr]) continue;
				hist_beam->Fill(pos.Z(),pos.X());	
				hist_beamZY->Fill(pos.Z(),pos.Y());	
				min_z = std::min(pos.Z(),min_z);
			}
		}//nitr;
		Quadratic->SetParameter(2,-7.6e-5);//~1.8GeV;
		int ent_cut = 10;
		if(hist_beam->GetEffectiveEntries()<ent_cut)continue;
		hist_beam->Fit("quadratic","QRNB");
		hist_beamZY->Fit("linear","QRNB");
		peak = Linear->GetParameter(0);
		double par0 = Quadratic->GetParameter(0);
		double par1 = Quadratic->GetParameter(1);
		double par2 = Quadratic->GetParameter(2);
		double beamv = Linear->GetParameter(1);
		double rad = 1/(2*par2);// x = p0 + p1 z + p2 z^2 -> BendingAngle = (dx/dz)_zOut-(dx/dz)_zIn 2*p2*(zOut-zIn). rad = arclength/angle ~ (zOut-zIn)/BendingAngle = 1/(2*p2)
		double mom = rad*(Const*dMagneticField);
		if(hist_beam->GetEntries()>ent_cut  and abs(mom)>500 and min_z<-150 )isAccidental=true;	
		if(isAccidental){
			beam_p0.push_back(par0);
			beam_p1.push_back(par1);
			beam_p2.push_back(par2);
			beam_y.push_back(peak);
			beam_v.push_back(beamv);
		}
	}//np = npeaks
	//	std::cout<<"NBeams = "<<beam_y.size()<<std::endl;
}


void  
TPCBeamRemover::SortAccidentalBeam(){ std::vector<TPCCluster*> triggered_hit_array;std::vector<int> triggered_hit_layer;
	triggered_hit_array.clear();triggered_hit_layer.clear();
	std::vector<TPCCluster*> accidental_hit_array;std::vector<int> accidental_hit_layer;
	accidental_hit_array.clear();accidental_hit_layer.clear();
	for(int ih = 0; ih<m_Cl_array.size();++ih){
		auto mhit = m_Cl_array[ih];
		auto pos = mhit->GetPosition();
		m_Accidental.push_back(IsAccidental(pos));
	}
	for(int ih = 0; ih<m_Cl_array.size();++ih){
		if(m_Accidental[ih]){
			accidental_hit_array.push_back(m_Cl_array[ih]);
			accidental_hit_layer.push_back(m_layer_info[ih]);
		}
		else {
			triggered_hit_array.push_back(m_Cl_array[ih]);
			triggered_hit_layer.push_back(m_layer_info[ih]);
		} 
	}
	m_Cl_array.clear();
	m_ClCont_array.clear();
	m_ClCont_array.resize(NumOfLayersTPC);
	m_BeamClCont_array.clear();
	m_BeamClCont_array.resize(NumOfLayersTPC);
	for(int ih = 0;ih<triggered_hit_array.size();++ih){
		int layer = triggered_hit_layer[ih];
		m_ClCont_array[layer].push_back(triggered_hit_array[ih]);
	}
	for(int ih = 0;ih<accidental_hit_array.size();++ih){
		int layer = accidental_hit_layer[ih];
		m_BeamClCont_array[layer].push_back(accidental_hit_array[ih]);
	}
}
bool
TPCBeamRemover::IsAccidental(TPCHit* hit){
	auto pos = hit->GetPosition();
	return IsAccidental(pos);
}
bool
TPCBeamRemover::IsAccidental(TVector3 pos){
	double x = pos.X(),y=pos.Y(),z=pos.Z();
	if(y>y_min and y<y_max) return false;
	for(int i=0;i<beam_y.size();++i){
		double p0= beam_p0[i],p1=beam_p1[i],p2=beam_p2[i],by=beam_y[i],bv=beam_v[i];
		double beamX = p0+p1*z+p2*z*z;
		if( (abs(by+z*bv-y)<Ycut)and(abs(beamX-x)<Xcut) ){
			return true;
		}
	}
	return false;
}



#endif
