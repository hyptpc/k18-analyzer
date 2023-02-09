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
	double ZTarget = -143;
  static std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
  static TF1 fint("fint",s_tmp.c_str(),-4.,4.);

	// TPC Tracking

}
TPCBeamRemover::TPCBeamRemover(std::vector<TPCClusterContainer>ClCont){
	enable = true;
	enableHough = true;
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
	hist_y->Reset();
	hist_beam->Reset();
	hist_beamZY->Reset();
	HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
	HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
	dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
	hist_Ci = new TH3D("hist_circle_acc","rd(mm);theta(rad);p(MeV/c)",
			nBin_rdiff,rdiff_min,rdiff_max,
			nBin_theta,theta_min,theta_max,
			nBin_p,p_min,p_max);
	hist_YTheta = new TH2D("histY_acc","theta(deg),r(mm)",
			thetaY_ndiv,thetaY_min,thetaY_max,
			r_ndiv,r_min,r_max);
	nh = m_Cl_array.size();

	debug::ObjectCounter::increase(ClassName());
}

TPCBeamRemover::~TPCBeamRemover()
{
	delete Quadratic;
	delete Linear;
	delete Arc;
	delete hist_beam;
	delete hist_beamZY;
	delete hist_y;
	delete hist_Ci;
	delete hist_YTheta;
	debug::ObjectCounter::decrease(ClassName());
}


int 
TPCBeamRemover::SearchPeaks(TH1D* hist,std::vector<double> &peaks){
	TSpectrum spec(30);
	double sig=1,th=0.2;
	int npeaks = spec.Search(hist,sig,"goff",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	m_PeakCl_array.resize(npeaks);
	h_cx.resize(npeaks);
	h_cy.resize(npeaks);
	h_z0.resize(npeaks);
	h_r.resize(npeaks);
	h_dz.resize(npeaks);
	h_flag.resize(npeaks);
	for(auto hit : m_Cl_array){
		auto pos = hit->GetPosition();
		double x = pos.X(),y=pos.Y();
		if(x>299792458)x+=0;//To avoid warning message.
		for(int ip = 0;ip<npeaks;++ip){
			double peak = peaks.at(ip);
			if(x_min<x and x<x_max and y_min<y and y<y_max)continue;
			if(abs(peak - y)<Ywidth){
				m_PeakCl_array[ip].push_back(hit);
				break;
			}
		}
	}
	return npeaks;
}


void 
TPCBeamRemover::SearchAccidentalBeam(double xmin = -30., double xmax = 30.,double ymin = -50.,double ymax = 50.){
	x_min=xmin;x_max=xmax;
	y_min=ymin;y_max=ymax;
	if(!enable) return;
	for(auto hitcl:m_Cl_array){
		TVector3 pos = hitcl->GetPosition();
		double x = pos.X();
		double y = pos.Y();
		if(x_min<x and x<x_max and y_min<y and y<y_max)continue;
		else hist_y->Fill(y);
	}
	np = SearchPeaks(hist_y,peaks);
	m_Accidental.resize(nh,false);

	for(int i = 0;i<np;++i){
		hist_beam->Reset();
		hist_beamZY->Reset();
		for(auto hitp:m_PeakCl_array[i]){
			TVector3 pos = hitp->GetPosition();
			double x=pos.X(),y=pos.Y(),z=pos.Z();
			if(x_min<x and x<x_max and y_min<y and y<y_max)continue;
			hist_beam->Fill(z,x);	
			hist_beamZY->Fill(z,y);	
		}
		if(enableHough) DoHoughSearch( i);
		else DoQuadraticSearch( i);
	}//np = npeaks
	int cnt=0;
	if(enableHough){
		for(int ip=0;ip < np;++ip){
			int nb = h_cx[ip].size();
			for(int ib = 0; ib< nb; ++ib){
				Allh_cx.push_back(h_cx[ip][ib]);
				Allh_cy.push_back(h_cy[ip][ib]);
				Allh_z0.push_back(h_z0[ip][ib]);
				Allh_r.push_back(h_r[ip][ib]);
				Allh_dz.push_back(h_dz[ip][ib]);
				cnt++;
			}
		}
	}
}


int
TPCBeamRemover::IsAccidental(TPCHit* hit){
	auto pos = hit->GetPosition();
	auto res = hit->GetResolutionVect();
	double Adist=0;
	int flag =  IsAccidental(pos,res,Adist);
	if(Adist< 5000) hit->SetHoughDist(Adist);
	return flag;
}

int
TPCBeamRemover::IsAccidental(TVector3 pos,TVector3 res, double& Adist){
	double x = pos.X(),y=pos.Y(),z=pos.Z();
	double sx = res.X(),sy=res.Y(),sz=res.Z();
	//	if(x_min<x and x<x_max and y_min<y and y<y_max) return false;
	if(enableHough){
		int nb = Allh_cx.size();
		double min_delta = 1e9;
		int beamID = -1;
		for(int ib = 0; ib < nb; ++ib){
			double cx = Allh_cx[ib];	
			double cy = Allh_cy[ib];	
			double z0 = Allh_z0[ib];	
			double r = Allh_r[ib];	
			double dz = Allh_dz[ib];	
			TVector3 pos_(-x,z-ZTarget,y);
			double fpar[8] = {cx,cy,z0,r,dz,pos_.X(),pos_.Y(),pos_.Z()};
  		fint.SetParameters(fpar);
  		double t = fint.GetMinimumX();
			TVector3 HelixPos(fpar[0]+fpar[3]*cos(t),fpar[1]+fpar[3]*sin(t),fpar[2]+(fpar[4]*fpar[3]*t));
			TVector3 dX = HelixPos-pos_;
			double delx = dX.X(),dely = dX.Y(),delz=dX.Z();
			double delta = sqrt(delx*delx/sx/sx + dely*dely/sy/sy+delz*delz/sz/sz);
			if( delta<min_delta){
				min_delta = delta;
				beamID = ib;
			}
		}
		Adist = min_delta;

		if(min_delta < 10){
			return 1000+beamID;
		}else{
			return 0;
		}
	}else{//Quadratic mode (Legacy)
		for(int ib=0;ib<beam_y.size();++ib){
			double p0= beam_p0[ib],p1=beam_p1[ib],p2=beam_p2[ib],by=beam_y[ib],bv=beam_v[ib];
			double beamX = p0+p1*z+p2*z*z;
			if( (abs(by+z*bv-y)<Ycut)and(abs(beamX-x)<Xcut) ){
				return 1000+ib;
			}
		}
	}
	return 0;
}
void 
TPCBeamRemover::DoHoughSearch(int i){
	DoCircleHough(i);
	DoYThetaFit(i);
//	DoYThetaHough(i);
}
void 
TPCBeamRemover::DoCircleHough(int i){
	int nhits = m_PeakCl_array[i].size();
	MaxNBeam = 3;
	h_flag[i].resize(nhits);
	std::vector<double> hough_x;
	std::vector<double> hough_y;
	std::vector<double> hough_z;
	std::vector<double> hough_p;
	std::vector<double> hough_rd;
	std::vector<double> hough_theta;
	for(int ib = 0 ;ib<MaxNBeam;++ib){
		hist_Ci->Reset();
		for(int ih = 0; ih<nhits;++ih){
			if(h_flag[i].at(ih)>0) continue;
			auto pos = m_PeakCl_array[i].at(ih)->GetPosition();
			for(int ird = 0;ird<nBin_rdiff;++ird){
				double rd = hist_Ci->GetXaxis()->GetBinCenter(ird+1);
				for(int ip=0;ip<nBin_p;++ip){
					double x = -pos.x();
					double y = pos.z()-ZTarget;
					double p = hist_Ci->GetZaxis()->GetBinCenter(ip+1);
					double r = p/(Const*dMagneticField);
					double a = 2.*(r+rd)*y;	
					double b = 2.*(r+rd)*x;	
					double c = -1.*(rd*rd+2.*r*rd+x*x+y*y);
					double r0 = sqrt(a*a+b*b);
					if(fabs(-1*c/r0)>1){
						continue;
					}
					double theta1_alpha = asin(-1.*c/r0);
					double theta2_alpha;
					if(theta1_alpha>0.)
						theta2_alpha = acos(-1.) - theta1_alpha;
					else
						theta2_alpha = -1.*acos(-1.) - theta1_alpha;
					double theta_alpha = atan2(b, a);
					double xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
					double ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
					double xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
					double ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
					double theta1 = atan2(ycenter1, xcenter1)+2*acos(-1);
					double theta2 = atan2(ycenter2, xcenter2)+2*acos(-1);
					hist_Ci->Fill(rd, theta1, p);
					hist_Ci->Fill(rd, theta2, p);
				}// ip
			}//ird
		}//ih
		int maxbin = hist_Ci->GetMaximumBin();
		int mx,my,mz;
		hist_Ci->GetBinXYZ(maxbin, mx, my, mz);
		hough_x.push_back(mx);
		hough_y.push_back(my);
		hough_z.push_back(mz);
		hough_rd.push_back(hist_Ci->GetXaxis()->GetBinCenter(mx));
		hough_theta.push_back( hist_Ci->GetYaxis()->GetBinCenter(my));
		hough_p.push_back(hist_Ci->GetZaxis()->GetBinCenter(mz));
		double hough_r = hough_p[ib]/(Const*dMagneticField);
		double hough_cx = (hough_r + hough_rd[ib])*cos(hough_theta[ib]);
		double hough_cy = (hough_r + hough_rd[ib])*sin(hough_theta[ib]);
		int hough_count = 0;
		double minz = 1115.683;
		double maxz = -1115.683;
		std::vector<double>zdist;
		for( int ih = 0; ih < nhits; ++ih){	
			auto pos = m_PeakCl_array[i].at(ih)->GetPosition();
			double x = -pos.X();
			double y = pos.Z() - ZTarget;
			double r_cal = sqrt(pow(x-hough_cx,2)+pow(y-hough_cy,2));
			double dist = abs(r_cal-hough_r);
			if(dist < MaxHoughWindow){
				if(pos.Z()<minz) minz = pos.Z();
				if(pos.Z()>maxz) maxz = pos.Z();
				h_flag[i].at(ih) += pow(2,ib);
				hough_count++;
				zdist.push_back(pos.Z());
			}
		}
		std::sort(zdist.begin(),zdist.end());
		for(int iz = 0; iz<hough_count-1;++iz){
			double dist = abs(zdist.at(iz)-zdist.at(iz+1));
			if(dist>125) hough_count = 0;
		}
		if(hough_count > MinBeamHit and minz < MinZCut and maxz > MaxZCut){
			h_cx[i].push_back(hough_cx);
			h_cy[i].push_back(hough_cy);
			h_r[i].push_back(hough_r);
		}
	}//ib
	int nbeam = h_cx[i].size();
	for(int ib = 0;ib<nbeam;++ib){
		double hcx = h_cx[i].at(ib);	
		double hcy = h_cy[i].at(ib);	
		double hr = h_r[i].at(ib);	
		std::vector<TVector3>AccArr(0);
		for(int ih = 0; ih<nhits; ++ih){
			if(IsThisBeam(h_flag[i].at(ih),ib)){
				auto pos = m_PeakCl_array[i].at(ih)->GetPosition();
				TVector3 HelixPos(-pos.X(),pos.Z()-ZTarget,pos.Y());
				AccArr.push_back(pos);
			}
		}
		double params[3] = {0,0,0};
		CircleFit(AccArr,params);	
		h_cx[i].at(ib)= params[0];
		h_cy[i].at(ib)= params[1];
		h_r[i].at(ib)= params[2];
	}
	
}
void 
TPCBeamRemover::DoYThetaHough(int i){
	int nb = h_cx[i].size();
	int nhits = m_PeakCl_array[i].size();
	for(int ib = 0; ib<nb;++ib){
		double hcx = h_cx[i].at(ib);
		double hcy = h_cy[i].at(ib);
		double hr = h_r[i].at(ib);
		hist_YTheta->Reset();
		for(int ih = 0;ih<nhits;++ih){
			if(h_flag[i].at(ih) != pow(2,ib)) continue;
			auto pos = m_PeakCl_array[i].at(ih)->GetPosition();
			double x = -pos.X();
			double y = pos.Z()-ZTarget;
			double z = pos.Y();
			for(int ti=0; ti<thetaY_ndiv; ti++){
				double theta = thetaY_min+ti*(thetaY_max-thetaY_min)/thetaY_ndiv;
				double tmp_t = atan2(y - hcy,x-hcx);
				double tmp_xval = hr * tmp_t;
				hist_YTheta->Fill(theta, cos(theta*acos(-1.)/180.)*tmp_xval+sin(theta*acos(-1.)/180.)*z);
			}
		}//ih
    int maxbinY = hist_YTheta->GetMaximumBin();
    int mxY,myY,mzY;
    hist_YTheta->GetBinXYZ(maxbinY, mxY, myY, mzY);
    double mtheta = hist_YTheta->GetXaxis()->GetBinCenter(mxY)*acos(-1.)/180.;
    double mr = hist_YTheta->GetYaxis()->GetBinCenter(myY);
		double p0 = mr/sin(mtheta);
    double p1 = -cos(mtheta)/sin(mtheta);
		h_z0[i].push_back(p0);
		h_dz[i].push_back(p1);
		for(int ih = 0;ih<nhits;++ih){
			if(h_flag[i].at(ih) != pow(2,ib)) continue;
			auto pos = m_PeakCl_array[i].at(ih)->GetPosition();
	  	double tmpx = -pos.x();
	  	double tmpy = pos.z() - ZTarget;
	  	double tmpz = pos.y();
	  	double tmp_t = atan2(tmpy - hcy,tmpx - hcx);
	  	double tmp_xval = hr * tmp_t;
	  	double distY = fabs(p1*tmp_xval - tmpz + p0)/sqrt(pow(p1, 2) + 1);
			if(distY > MaxHoughWindowY) h_flag[i].at(i) = 0;
		}
	}
}
void 
TPCBeamRemover::DoYThetaFit(int i){
	int nb = h_cx[i].size();
	int nhits = m_PeakCl_array[i].size();
	for(int ib = 0; ib<nb;++ib){
		double hcx = h_cx[i].at(ib);
		double hcy = h_cy[i].at(ib);
		double hr = h_r[i].at(ib);
		std::vector<double> ttmp;
		std::vector<TVector3> AccArr(0);
		for(int ih = 0;ih<nhits;++ih){
			if(IsThisBeam(h_flag[i].at(ih),ib)){
				auto pos = m_PeakCl_array[i].at(ih)->GetPosition();
				double x = -pos.X();
				double y = pos.Z()-ZTarget;
				double z = pos.Y();
				double t = atan2(y - hcy,x-hcx);
				double rt = t * hr;
				ttmp.push_back(t);
				AccArr.push_back(TVector3(rt,z,y));			
			}
		}//ih
		double params[2]={0,0};
		LinearFit(AccArr,params);
		std::sort(ttmp.begin(),ttmp.end());
		double p0 = params[0];
		double p1 = params[1];
		h_z0[i].push_back(p0);
		h_dz[i].push_back(p1);
	}
}
int
TPCBeamRemover::CompareHough(TVector3 pos, std::vector<double> hcx,std::vector<double>hcy,std::vector<double> hr){ 
			
	int beamid = -1;
	double x = -pos.X();
	double y = pos.Z() - ZTarget;
	int nbeam = hcx.size();
	double min_dist = MaxHoughWindow;
	for(int ib=0;ib<nbeam;++ib){
		double hough_cx = hcx.at(ib),hough_cy = hcy.at(ib),hough_r = hr.at(ib);
		double r_cal = sqrt(pow(x-hough_cx,2)+pow(y-hough_cy,2));
		double dist = abs(r_cal-hough_r);
		if(dist<min_dist){
			min_dist = dist;
			beamid = ib;
		}
	}
	return beamid;
}
bool 
TPCBeamRemover::IsThisBeam(int hflag, int ib){
	int mask = pow(2,ib+1);
	int flag= hflag%mask;
	return flag/int(pow(2,ib));
}

void 
TPCBeamRemover::DoQuadraticSearch(int i){
	double peak = peaks.at(i);
	Linear->SetParameter(0,peak);
	Linear->SetParLimits(0,-400,400);

	Linear->SetParameter(1,0);
	Linear->SetParLimits(1,-0.1,0.1);

	Quadratic->SetParameter(0,0);
	Quadratic->SetParLimits(0,-120,120);

	Quadratic->SetParameter(1,-0.05);
	Quadratic->SetParLimits(1,-0.1,-0.01);

	Quadratic->SetParameter(2,-7.6e-5);
	Quadratic->SetParLimits(2,-3e-4,-3e-5);
	double min_z = MinZCut;
	double dist_cut[3]={20,10,5};
	double dist_cutY[3]={10,7.5,5};
	bool isAccidental = false;
	for(int nitr=0;nitr<3;++nitr){
		min_z = 255;
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
		for(auto hitp:m_PeakCl_array[i]){
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
	int ent_cut = 15;
	if(hist_beam->GetEffectiveEntries()<ent_cut)return;
	hist_beam->Fit("quadratic","QRNB");
	hist_beamZY->Fit("linear","QRNB");
	peak = Linear->GetParameter(0);
	double par0 = Quadratic->GetParameter(0);
	double par1 = Quadratic->GetParameter(1);
	double par2 = Quadratic->GetParameter(2);
	double beamv = Linear->GetParameter(1);
	double rad = 1/(2*par2);// x = p0 + p1 z + p2 z^2 -> BendingAngle = (dx/dz)_zOut-(dx/dz)_zIn 2*p2*(zOut-zIn). rad = arclength/angle ~ (zOut-zIn)/BendingAngle = 1/(2*p2)
	double mom = rad*(Const*dMagneticField);
	if(hist_beam->GetEntries()>ent_cut  and abs(mom)>500 and min_z<MinZCut )isAccidental=true;	
	if(isAccidental){
		beam_p0.push_back(par0);
		beam_p1.push_back(par1);
		beam_p2.push_back(par2);
		beam_y.push_back(peak);
		beam_v.push_back(beamv);
	}
}
void
TPCBeamRemover::CircleFit(std::vector<TVector3> posarr,double* param){
	int n =posarr.size();
	double Sumx=0,Sumy=0;
	double Sumx2=0,Sumy2=0,Sumxy=0;
	double Sumx3=0,Sumy3=0,Sumx2y=0,Sumy2x=0;
	for (Int_t i=0;i<n;i++) {
		double x = posarr.at(i).X();
		double y = posarr.at(i).Y();
		Sumx+=x;Sumy+=y;
		Sumx2+=x*x;Sumy2+=y*y;Sumxy+=x*y;
		Sumx3+=x*x*x;Sumy3+=y*y*y;Sumx2y+=x*x*y;Sumy2x+=y*y*x;
	}
	double a_1 = Sumx3 + Sumy2x, a_2 = Sumx2, a_3 = Sumx;
	double b_1 = Sumy3 + Sumx2y, b_2 = Sumy2, b_3 = Sumy;
	double c_1 = a_2+b_2,c_2 = Sumx,c_3 = Sumy;
	double A = (a_1*(b_2*n-b_3*c_3)+a_3*(b_1*c_3-b_2*c_1)+Sumxy*(b_3*c_1-b_1*n))
		/ ( Sumxy*(n*Sumxy-a_3*c_3-b_3*c_2)+a_2*b_3*c_3+a_3*b_2*c_2-a_2*b_2*n);
	double B = (b_1*(a_2*n-a_3*c_2)+b_3*(a_1*c_2-a_2*c_1)+Sumxy*(a_3*c_1-a_1*n))      
		/ ( Sumxy*(n*Sumxy-a_3*c_3-b_3*c_2)+a_2*b_3*c_3+a_3*b_2*c_2-a_2*b_2*n);
	double C = (c_1*(a_2*b_2-Sumxy*Sumxy)+c_2*(Sumxy*b_1-a_1*b_2)+c_3*(Sumxy*a_1-b_1*a_2))
		/ ( Sumxy*(n*Sumxy-a_3*c_3-b_3*c_2)+a_2*b_3*c_3+a_3*b_2*c_2-a_2*b_2*n);
	double cx = -0.5 * A;
	double cy = -0.5 * B;
	double rad = sqrt(cx*cx+cy*cy-C);
	param[0] = cx;
	param[1] = cy;
	param[2] = rad;
}
void
TPCBeamRemover::LinearFit(std::vector<TVector3> posarr,double* param){
	// y = ax + b;
	int n =posarr.size();
	double Sumx=0,Sumy=0;
	double Sumx2=0,Sumy2=0,Sumxy=0;
	for (Int_t i=0;i<n;i++) {
		double x = posarr.at(i).X();
		double y = posarr.at(i).Y();
		Sumx+=x;Sumy+=y;
		Sumx2+=x*x;Sumy2+=y*y;Sumxy+=x*y;
	}
	double a = (Sumxy-Sumx*Sumy)/( Sumx*Sumxy+Sumx2);	
	double b = (Sumx2*Sumy-Sumxy*Sumxy)/( Sumx*Sumxy+Sumx2);	
	param[0]=b;
	param[1]=a;
}


#endif
