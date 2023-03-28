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
#define DebugDisplay 0
namespace {
	const auto& gConf = ConfMan::GetInstance();
	const auto& gGeom = DCGeomMan::GetInstance();
	const auto& gUser = UserParamMan::GetInstance();
	const auto& zTarget    = gGeom.LocalZ("Target");
	const auto& zK18Target = gGeom.LocalZ("K18Target");
	double ZTarget = -143;
  static std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
  static TF1 fint("fint",s_tmp.c_str(),-4.,4.);


	const Int_t nBin_rdiff = 220;
	const Double_t rdiff_min = -220.;
	const Double_t rdiff_max = 220.;
	const Int_t nBin_theta = 180;
	const Double_t theta_min = -0.3*acos(-1);
	const Double_t theta_max = 0.3*acos(-1);//Charge < -1
	const Int_t nBin_p = 100;
	const Double_t p_min = 1600.;//MeV/c
	const Double_t p_max = 2000.;//MeV/c


	double Y_shift = 500;//To avoid Flip on hough angle.
	const int    thetaY_ndiv =  180;
	const double thetaY_min  =  60.;
	const double thetaY_max  =  120.;
	const int    r_ndiv =  1000;
	const double r_min  = 000.;
	const double r_max  = 1000.;
	
	const int ZYtheta_ndiv = 120;
	const double ZYtheta_min = 0.45*acos(-1);
	const double ZYtheta_max = 0.55*acos(-1);//-0.4*acos(-1);
	const int rho_ndiv = 500;
	const double rho_min = -500;
	const double rho_max = 500;


	TH3D* hist_Ci = new TH3D("hist_circle_acc","rd(mm);theta(rad);p(MeV/c)",
			nBin_rdiff,rdiff_min,rdiff_max,
			nBin_theta,theta_min,theta_max,
			nBin_p,p_min,p_max);
	TH2D* hist_YTheta = new TH2D("histY_acc","theta(deg),r(mm)",
			thetaY_ndiv,thetaY_min,thetaY_max,
			r_ndiv,r_min,r_max);

	TH1D* hist_y = new TH1D("peakhist_y","peakhist_y",200,-500,500);
	TH2D* hist_ZY = new TH2D("histZY","theta[rad] : rho[mm]"
			,ZYtheta_ndiv, ZYtheta_min,ZYtheta_max
			,rho_ndiv, rho_min,rho_max);


	
	
	
	void CircleFit(std::vector<TVector3> posarr,double* param){
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
#if DebugDisplay >1
		std::cout<<Form("(cx,cy,rad) = (%f,%f,%f)",cx,cy,rad)<<std::endl;
#endif
	}
	void
		LinearFit(std::vector<TVector3> posarr,double* param){
			// y = ax + b;
			int n =posarr.size();
			double Sumx=0,Sumy=0;
			double Sumx2=0,Sumxy=0;
			for (Int_t i=0;i<n;i++) {
				double x = posarr.at(i).X();
				double y = posarr.at(i).Y();
				Sumx+=x;Sumy+=y;
				Sumx2+=x*x;Sumxy+=x*y;
			}
			double b = (Sumx*Sumxy-Sumx2*Sumy)/( Sumx*Sumx-n*Sumx2);	
			double a = (Sumx*Sumy-n*Sumxy)/( Sumx*Sumx-n*Sumx2);	
			param[0]=b;//b = z0;
			param[1]=a;//a = dz;
#if DebugDisplay >1
			std::cout<<Form("(z0,dz) = (%f,%f)",b,a)<<std::endl;
#endif
		}
		bool SortY(TVector3 a, TVector3 b){
			return a.Y()<b.Y();
		}
}
TPCBeamRemover::TPCBeamRemover(std::vector<TPCClusterContainer>ClCont){
	enable = true;
	enableHough = true;
	for(int layer=0;layer<NumOfLayersTPC;layer++){
		TPCClusterContainer layer_cl_array;
		layer_cl_array.clear();
		for(auto hitcl:ClCont[layer]){
			layer_cl_array.push_back(hitcl);
			m_Cl_array.push_back(hitcl);
		}
		m_ClCont_array.push_back(layer_cl_array);
	}
	hist_y->Reset();
	HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
	HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");
	dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);
	nh = m_Cl_array.size();
	m_Cl_Flag.resize(nh,0);
	m_ZYHough_Flag.resize(nh,0);
	debug::ObjectCounter::increase(ClassName());
}

TPCBeamRemover::~TPCBeamRemover()
{
	m_Track_array.clear();	
	debug::ObjectCounter::decrease(ClassName());
}


int 
TPCBeamRemover::SearchPeaks(TH1D* hist,std::vector<double> &peaks){
	TSpectrum spec(30);
	double sig=2,th=0.1;
	int npeaks = spec.Search(hist,sig,"goff",th);
	double* peaks_ = spec.GetPositionX();
	peaks.clear();
	for(int i=0;i<npeaks;++i){
		peaks.push_back(peaks_[i]);
	}
	m_PeakCl_array.resize(npeaks);
	m_PeakCl_ID.resize(npeaks);
	h_cx.resize(npeaks);
	h_cy.resize(npeaks);
	h_z0.resize(npeaks);
	h_r.resize(npeaks);
	h_dz.resize(npeaks);
	h_flag.resize(npeaks);
	for(int ih = 0; ih < nh; ++ih){
		auto hit = m_Cl_array.at(ih);
		auto pos = hit->GetPosition();
		double x = pos.X(),y=pos.Y();
		if(x>299792458)x+=0;//Dummy statement. But you will get warning message without this.
		for(int ip = 0;ip<npeaks;++ip){
			double peak = peaks.at(ip);
//			if(x_min<x and x<x_max and y_min<y and y<y_max)continue;
			if(abs(peak - y)<Ywindow){
				m_PeakCl_array[ip].push_back(hit);
				m_PeakCl_ID[ip].push_back(ih);
//				break;
			}
		}
	}
#if DebugDisplay 
	std::cout<<FUNC_NAME<<"Number of peaks = "<<npeaks<<std::endl;
#endif
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
	for(int ip = 0;ip<np;++ip){
		for(auto hitp:m_PeakCl_array[ip]){
			TVector3 pos = hitp->GetPosition();
			double x=pos.X(),y=pos.Y(),z=pos.Z();
			if(x_min<x and x<x_max and y_min<y and y<y_max)continue;
		}
		DoPeakHoughSearch( ip);
	}//np = npeaks
	DoZYHough();	
	int nZYHough = m_ZYHoughCl_array.size();
	for(int inz=0;inz<nZYHough;++inz){
		DoZYHoughSearch( inz);
	}
	int cnt=0;
	for(int ip=0;ip < peaks.size();++ip){
		int nb = h_cx[ip].size();
		for(int ib = 0; ib< nb; ++ib){
			Allh_cx.push_back(h_cx[ip][ib]);
			Allh_cy.push_back(h_cy[ip][ib]);
			Allh_z0.push_back(h_z0[ip][ib]);
			Allh_r.push_back(h_r[ip][ib]);
			Allh_dz.push_back(h_dz[ip][ib]);
			AllPeaks.push_back(peaks.at(ip));
			cnt++;
		}
	}
}


int
TPCBeamRemover::IsAccidental(TPCHit* hit){
	auto pos = hit->GetPosition();
	auto res = hit->GetResolutionVect();
	double PullDist=0;
	int flag =  IsAccidental(pos,res,PullDist);
	if(PullDist< 1000) hit->SetPull(PullDist);
	return flag;
}

int
TPCBeamRemover::IsAccidental(TVector3 pos,TVector3 res, double& PullDist){
	double x = pos.X(),y=pos.Y(),z=pos.Z();
	double sx = res.X(),sy=res.Y(),sz=res.Z();
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
		double t = fint.GetMinimumX(0*acos(-1),2*acos(-1));
		TVector3 HelixPos(fpar[0]+fpar[3]*cos(t),fpar[1]+fpar[3]*sin(t),fpar[2]+(fpar[4]*fpar[3]*t));
		TVector3 dX = HelixPos-pos_;
		double delx = dX.X(),dely = dX.Y(),delz=dX.Z();
		double delta = sqrt(delx*delx/sx/sx + dely*dely/sy/sy+delz*delz/sz/sz);
		if( delta<min_delta){
			min_delta = delta;
			beamID = ib;
		}
		PullDist = min_delta;
	}
	if(min_delta < MaxPullDist){
		return Acc_flag_base+beamID;
	}else{
		return 0;
	}
	return 0;
}

void
TPCBeamRemover::ConstructAccidentalTracks(){
  int nacc = GetNAccBeam();
	std::vector<std::vector<TPCClusterContainer>>ClContVect;
	ClContVect.resize(nacc);
	for(int iacc=0; iacc<nacc;++iacc){
		ClContVect.at(iacc).resize(NumOfLayersTPC);
		auto Track = new TPCLocalTrackHelix();
		Track->SetAcx(GetParameter(0).at(iacc));
		Track->SetAcy(GetParameter(1).at(iacc));
		Track->SetAz0(GetParameter(2).at(iacc));
		Track->SetAr(GetParameter(3).at(iacc));
		Track->SetAdz(GetParameter(4).at(iacc));
		Track->SetFlag(8);
		Track->SetHoughFlag(Acc_flag_base+iacc);
		m_Track_array.push_back(Track);
	}
	
	for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(auto cl:m_ClCont_array[layer]){
			auto hit = cl -> GetMeanHit();
			auto pos = hit->GetPosition();
			int accf = IsAccidental(hit);
			int iacc = accf - Acc_flag_base;	
			if(hit->GetHoughFlag() !=0) continue;
			if((-1 > iacc or iacc > 100  )){
				continue;
			}
			else{
				hit->SetHoughFlag(accf);
				m_Track_array.at(iacc)->AddTPCHit(new TPCLTrackHit(hit));
			}
			for(auto peak : AllPeaks){
				if(abs(peak - pos.Y())<Ywindow+10){
					ClContVect.at(iacc)[layer].push_back(cl);
				}
			}
		}
	}
	
	int nt = nacc;	
	for(int it = 0; it < nt; ++it){
		auto Track = m_Track_array.at(it);
		DoHelixFit(Track,ClContVect.at(it),8);
	}
}






void 
TPCBeamRemover::DoPeakHoughSearch(int i){
	auto ClCont = m_PeakCl_array[i];
	DoCircleHough(ClCont, i);
	DoYThetaHough(ClCont, i);
}
void 
TPCBeamRemover::DoZYHoughSearch(int i){
	auto ClCont = m_ZYHoughCl_array[i];
//	int np = m_PeakCl_array.size();
	DoCircleHough(ClCont, np+i);
	DoYThetaHough(ClCont, np+i);
}


void 
TPCBeamRemover::DoCircleHough(TPCClusterContainer ClCont, int i){
	int nhits = ClCont.size();
#if DebugDisplay 
	std::cout<<FUNC_NAME<<"np / Number of peaks = "<<np<<" / "<<peaks.size()<<std::endl;
	if(i<np){
		std::cout<<FUNC_NAME<<"(peak# , peak, nhits )=  ( "<<i<<" , "<<peaks.at(i)<<" , "<<nhits<<" )"<<std::endl;
	}else{
		std::cout<<FUNC_NAME<<"(ZYHough# ,peak, nhits )=  ( "<<i-np<<" , "<<peaks.at(i)<<" , "<<nhits<<" )"<<std::endl;
	}
#endif
	MaxNBeam = 5;
	h_flag[i].resize(nhits);
	std::vector<double> hough_x;
	std::vector<double> hough_y;
	std::vector<double> hough_z;
	std::vector<double> hough_p;
	std::vector<double> hough_rd;
	std::vector<double> hough_theta;
	for(int ib = 0 ;ib<MaxNBeam;++ib){
		hist_Ci->Reset();
		int nh = 0;
		for(int ih = 0; ih<nhits;++ih){
			if(h_flag[i].at(ih)>0) continue;
			nh ++;
		}
		if(nh < 8) continue;
		for(int ih = 0; ih<nhits;++ih){
			if(h_flag[i].at(ih)>0) continue;
			auto pos = ClCont.at(ih)->GetPosition();
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
					double theta1 = atan2(ycenter1, xcenter1);
					double theta2 = atan2(ycenter2, xcenter2);
//					if(theta1<0) theta1+= 2*acos(-1);
//					if(theta2<0) theta2+= 2*acos(-1);
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
		double mrd = (hist_Ci->GetXaxis()->GetBinCenter(mx));
		double mtheta = ( hist_Ci->GetYaxis()->GetBinCenter(my));
		double mp = ( hist_Ci->GetZaxis()->GetBinCenter(mz));
		hough_rd.push_back(hist_Ci->GetXaxis()->GetBinCenter(mx));
		hough_theta.push_back( hist_Ci->GetYaxis()->GetBinCenter(my));
		hough_p.push_back(hist_Ci->GetZaxis()->GetBinCenter(mz));
		double hough_r = hough_p[ib]/(Const*dMagneticField);
		double hough_cx = (hough_r + hough_rd[ib])*cos(hough_theta[ib]);
		double hough_cy = (hough_r + hough_rd[ib])*sin(hough_theta[ib]);
#if DebugDisplay
		std::cout<<FUNC_NAME<<"(rd, theta,p, p0,p1  )= "<<Form("(%f, %f, %f)",mrd,mtheta,mp)<<std::endl;
#endif
		int hough_count = 0;
		double minz = 1115.683;
		double maxz = -1115.683;
		std::vector<double>zdist;
		for( int ih = 0; ih < nhits; ++ih){	
				//No cuts for h_flag in this loop.Double counting is intentionally allowed here. Note that this is not Hough transformation.
			auto pos = ClCont.at(ih)->GetPosition();
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
			double gap = abs(zdist.at(iz)-zdist.at(iz+1));
			if(gap>125) hough_count = 0;
		}
		if(hough_count > MinBeamHit and minz < MinZCut and maxz > MaxZCut){
#if DebugDisplay
		std::cout<<FUNC_NAME<<"param size : "<<h_cx.size()<<" i = "<<i<<std::endl;
#endif
			h_cx[i].push_back(hough_cx);
			h_cy[i].push_back(hough_cy);
			h_r[i].push_back(hough_r);
		if(i < np){
		for( int ih = 0; ih < nhits; ++ih){	
				auto pos = ClCont.at(ih)->GetPosition();
				double x = -pos.X();
				double y = pos.Z() - ZTarget;
				double r_cal = sqrt(pow(x-hough_cx,2)+pow(y-hough_cy,2));
				double dist = abs(r_cal-hough_r);
				if(dist < MaxHoughWindow){
					int id = m_PeakCl_ID[i].at(ih);
					m_Cl_Flag[id]+=1;
					//Sort Hits In PeakSearch
				}
			}
		}
#if DebugDisplay
			std::cout<<FUNC_NAME<<"track accepted."<<std::endl;
#endif
		}
		else{
#if DebugDisplay
			std::cout<<FUNC_NAME<<"track rejected. (minz ,maxz ) = ( "<<minz<<" , "<<maxz<<" )"<<std::endl;
#endif
		}
#if DebugDisplay
		std::cout<<FUNC_NAME<<"( ib, hough_count ) = ( "<<ib<<" , "<<hough_count<<" )"<<std::endl;
#endif
	}//ib
	int nbeam = h_cx[i].size();
	for(int ib = 0;ib<nbeam;++ib){
		std::vector<TVector3>AccArr;
		for(int ih = 0; ih<nhits; ++ih){
			auto pos = ClCont.at(ih)->GetPosition();
			int beamID = CompareHoughDist(pos,h_cx[i],h_cy[i],h_r[i]);
			if(IsThisBeam(beamID,ib)){
				TVector3 HelixPos(-pos.X(),pos.Z()-ZTarget,pos.Y());
				AccArr.push_back(HelixPos);
			}
		}
		double params[3] = {0,0,0};
		std::sort(AccArr.begin(),AccArr.end(),SortY);		
		std::vector<TVector3>AccArrFit;
		double x0 = AccArr.at(0).X();
		for(auto pv : AccArr){
//			if(abs(x0 - pv.X())>10) continue;
//			if(x0 - pv.X()>0) continue;
//			x0 = pv.X();
			AccArrFit.push_back(pv);
		}
		CircleFit(AccArrFit,params);	
		if(params[2]>3000 and params[0]>1000){//Update param only if radius>3000
			h_cx[i].at(ib)= params[0];
			h_cy[i].at(ib)= params[1];
			h_r[i].at(ib)= params[2];
		}
		else{
#if DebugDisplay > 1 
			std::cout<<FUNC_NAME<<"CircleFitFail! nhits = "<<AccArrFit.size()<<std::endl;
			std::cout<<Form("Params(%f,%f,%f)",params[0],params[1],params[2])<<std::endl;
#endif
			for(auto pv : AccArrFit){
#if DebugDisplay > 1
				if(params[2]<500)std::cout<<Form("Pos(%f,%f,%f)",-pv.X(),pv.Z(),pv.Y()+ZTarget)<<std::endl;
#endif
			}
		}
	}
	
}
void 
TPCBeamRemover::DoYThetaHough(TPCClusterContainer ClCont,int i){
	int nb = h_cx[i].size();
	int nhits = ClCont.size();
#if DebugDisplay 
	std::cout<<FUNC_NAME<<"(peak# , peak, nb, nhits )=  ( "<<i<<" , "<<peaks.at(i)<<" , "<<nb<<" , "<<nhits<<" )"<<std::endl;
#endif
	for(int ib = 0; ib<nb;++ib){
		double hcx = h_cx[i].at(ib);
		double hcy = h_cy[i].at(ib);
		double hr = h_r[i].at(ib);
		hist_YTheta->Reset();
		for(int ih = 0;ih<nhits;++ih){
//			if(h_flag[i].at(ih) != pow(2,ib)) continue;
			if(!IsThisBeam(h_flag[i].at(ih),ib)) continue;
			auto pos = ClCont.at(ih)->GetPosition();
			double x = -pos.X();
			double y = pos.Z()-ZTarget;
			double z = pos.Y()+Y_shift;
			for(int ti=0; ti<thetaY_ndiv; ti++){
				double theta = thetaY_min+ti*(thetaY_max-thetaY_min)/thetaY_ndiv;
				double tmp_t = atan2(y - hcy,x-hcx);
				if(tmp_t < 0) tmp_t += 2*acos(-1);
				tmp_t-=acos(-1);
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
#if DebugDisplay
		std::cout<<FUNC_NAME<<"(mr, mtheta , p0,p1  )= "<<Form("(%f, %f, %f, %f)",mr,mtheta,p0,p1)<<std::endl;
#endif
		int hough_count = 0;	
		for(int ih = 0;ih<nhits;++ih){
			if(!IsThisBeam(h_flag[i].at(ih),ib)) continue;
			auto pos = ClCont.at(ih)->GetPosition();
			double tmpx = -pos.x();
			double tmpy = pos.z() - ZTarget;
			double tmpz = pos.y()+Y_shift;
			double tmp_t = atan2(tmpy - hcy,tmpx - hcx);
			if(tmp_t < 0) tmp_t += 2*acos(-1);
			tmp_t-=acos(-1);
			double tmp_xval = hr * tmp_t;
			double distY = fabs(p1*tmp_xval - tmpz + p0)/sqrt(pow(p1, 2) + 1);
#if DebugDisplay > 1
				std::cout<<FUNC_NAME<<Form("(tx, z)= (%f,%f)  DistY= %f",tmp_xval,tmpz,distY)<<std::endl;
#endif
			if(distY > MaxHoughWindowY){
				h_flag[i].at(ih) -= pow(2,ib);
			}
			else{
				hough_count++;
			}
		}
		if(hough_count <MinBeamHit+1){
#if DebugDisplay
		std::cout<<FUNC_NAME<<Form("Deleting Track %d: %d/%d: nhits = %d",i,ib,nb,hough_count)<<std::endl;
#endif
			h_cx[i].erase(h_cx[i].begin()+ib);
			h_cy[i].erase(h_cy[i].begin()+ib);
			h_r[i].erase(h_r[i].begin()+ib);
			nb-=1;
			ib-=1;
		}
	}
	DoYThetaFit(ClCont, i);
}


void 
TPCBeamRemover::DoYThetaFit(TPCClusterContainer ClCont,int i){
#if DebugDisplay
	std::cout<<FUNC_NAME<<std::endl;
#endif
	int nb = h_cx[i].size();
	int nhits = ClCont.size();
#if DebugDisplay
	std::cout<<FUNC_NAME<<Form("(nb,nh) = (%d,%d)",nb,nhits)<<std::endl;
#endif
	for(int ib = 0; ib<nb;++ib){
		double hcx = h_cx[i].at(ib);
		double hcy = h_cy[i].at(ib);
		double hr = h_r[i].at(ib);
		std::vector<double> ttmp;
		std::vector<double> ztmp;
		std::vector<TVector3> AccArr(0);
		for(int ih = 0;ih<nhits;++ih){
#if DebugDisplay >1
	std::cout<<FUNC_NAME<<Form("h_flag[%d].at(%d) = %d",i,ih,h_flag[i].at(ih))<<std::endl;
#endif
			if(IsThisBeam(h_flag[i].at(ih),ib)){
				auto pos = ClCont.at(ih)->GetPosition();
				double x = -pos.X();
				double y = pos.Z()-ZTarget;
				double z = pos.Y();
				double t = atan2(y - hcy,x-hcx);
				if(t < 0) t += 2*acos(-1);
				double rt = t * hr;
				ttmp.push_back(t);
				ztmp.push_back(z);
				AccArr.push_back(TVector3(rt,z,y));			
			}
		}//ih
		std::sort(ttmp.begin(),ttmp.end());
		std::sort(ztmp.begin(),ztmp.end());
		int nz = ztmp.size();
		double z_med = ztmp.at(nz/2);
		double t_med = ttmp.at(nz/2);
		double params[2]={0,0};
		LinearFit(AccArr,params);

		double p0 = params[0];
		double p1 = params[1];
		double dist = abs( hr*t_med*p1 + p0 - z_med);
		if(dist>MaxHoughWindowY*2){
#if DebugDisplay
			std::cout<<"TPCBeamRemover::DoYThetaFit Z0 Fail, hr*dz ,t_med, z_med, p0, p1 ,dist= ( "<<hr*p1<<" , "<<t_med<<" , "<<z_med<<" , "<<p0<<" , "<<p1<<" , "<<dist<<" )"<<std::endl;
#endif
			p0 = z_med - hr*t_med*p1;
		}
		
		
		h_z0[i].push_back(p0);
		h_dz[i].push_back(p1);
#if DebugDisplay> 1 
		if(abs(params[1])>0.1){
			for(auto v :AccArr){
				std::cout<<Form("pos(%f,%f,%f)",v.X(),v.Y(),v.Z())<<std::endl;
			}
		}
#endif
	}
}

void
TPCBeamRemover::DoZYHough(){
	int MaxNt = 5;
	for(int it=0;it< MaxNt; ++it){
		hist_ZY->Reset();
		for(int ih=0;ih<nh;++ih){
			if(m_Cl_Flag.at(ih)) continue;
			if(m_ZYHough_Flag.at(ih)) continue;
			auto hit = m_Cl_array.at(ih);
			auto pos = hit->GetPosition();
			double x = pos.X(),y=pos.Y(),z=pos.Z();
			for(int ti = 0; ti<ZYtheta_ndiv;++ti){
				double theta = hist_ZY->GetXaxis()->GetBinCenter(ti+1);
				double rho = cos(theta)*z + sin(theta)*y;
				hist_ZY->Fill(theta,rho);
			}//ti
		}//ih
		int maxbin = hist_ZY->GetMaximumBin();
		int mx,my,mz;
		hist_ZY->GetBinXYZ(maxbin,mx,my,mz);
		double mtheta = hist_ZY->GetXaxis()->GetBinCenter(mx);
		double mrho = hist_ZY->GetYaxis()->GetBinCenter(my);
		double y0 = mrho/sin(mtheta);
		double v = -cos(mtheta)/sin(mtheta);
		std::vector<int> index;
		int nc = 0;
		for(int ih = 0; ih < nh; ++ih){
			if(m_Cl_Flag.at(ih)) continue;
			if(m_ZYHough_Flag.at(ih)) continue;
			auto hit = m_Cl_array.at(ih);
			auto pos = hit->GetPosition();
			double x = pos.X(),y=pos.Y(),z=pos.Z();
			double Ydist = sqrt( (y0 + v*z  - y)/(v*v+1));
			if(Ydist < 20){
				index.push_back(ih);
				nc++;
			}
		}
		TPCClusterContainer ZYClCont;
		if(nc> 8){
			for(auto ih : index){
				m_ZYHough_Flag.at(ih)=true;
				auto hit = m_Cl_array.at(ih);
				ZYClCont.push_back(hit);
			}
			m_ZYHoughCl_array.push_back(ZYClCont);
			peaks.push_back(y0);
			std::vector<int> hough_flag;	
			std::vector<double> hough_cx;	
			std::vector<double> hough_cy;	
			std::vector<double> hough_z0;	
			std::vector<double> hough_r;	
			std::vector<double> hough_dz;
			h_flag.push_back(hough_flag);
			h_cx.push_back(hough_cx);
			h_cy.push_back(hough_cy);
			h_z0.push_back(hough_z0);
			h_r.push_back(hough_r);
			h_dz.push_back(hough_dz);
#if DebugDisplay
			std::cout<<Form("Track %d: hits = %d , y = %f +  %f z ",it,nc,y0,v)<<std::endl;
			std::cout<<Form("Params: (rho , theta ) = ( %f , %f) ",mrho,mtheta)<<std::endl;
#endif
		}
	}//it
	int nZYHough = m_ZYHoughCl_array.size();
#if DebugDisplay
	std::cout<<ClassName()<<"::DoZYHough() "<<"Number of ZYHoughs : "<<nZYHough<<std::endl;
#endif
}






void
TPCBeamRemover::DoHelixFit(TPCLocalTrackHelix* Track,const std::vector<TPCClusterContainer>& ClCont,int MinNumOfHits = 8){
	int hf = Track->GetHoughFlag();
#if DebugDisplay  
	std::cout<<FUNC_NAME<<"HoughFlag= "<<hf<<std::endl;
#endif
	int it = hf - Acc_flag_base; 
	if(Track->DoFit(MinNumOfHits)){
    TPCLocalTrackHelix *copied_track = new TPCLocalTrackHelix(Track);
    bool Add_rescheck = false;
    for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
      for(auto cl:ClCont[layer]){
				TPCHit* hit = cl->GetMeanHit();
				if(hit->GetHoughFlag() != 0) continue;
				auto pos = hit->GetPosition();
				auto res = hit->GetResolutionVect();
				double residual = 0;
				if(copied_track->ResidualCheck(pos,res,residual,25)){
	  			copied_track->AddTPCHit(new TPCLTrackHit(hit));
	  			Add_rescheck = true;
				}
			}//cl
		}//layer
#if DebugDisplay  
	std::cout<<FUNC_NAME<<"ResCheck= "<<Add_rescheck<<std::endl;
#endif
		if(!Add_rescheck){
			delete copied_track;
		}
		else{
			if(copied_track->DoFit(MinNumOfHits)){
				copied_track->SetHoughFlag(hf);
				copied_track->SetFitFlag(2);
				m_Track_array.at(it) = copied_track;
				delete Track;
			}
		}
	}//if(DoFit)
	else{
		Track->SetFlag(8);
		Track->SetHoughFlag(hf);
		Track->SetParamUsingHoughParam();
	}
}

int
TPCBeamRemover::CompareHoughDist(TVector3 pos, std::vector<double> hcx,std::vector<double>hcy,std::vector<double> hr){ 
	int beamid = -1;
	double x = -pos.X();
	double y = pos.Z() - ZTarget;
	int nbeam = hcx.size();
	double min_dist = MaxHoughWindow/2;
	for(int ib=0;ib<nbeam;++ib){
		double hough_cx = hcx.at(ib),hough_cy = hcy.at(ib),hough_r = hr.at(ib);
		double r_cal = sqrt(pow(x-hough_cx,2)+pow(y-hough_cy,2));
		double dist = abs(r_cal-hough_r);
		if(dist<min_dist){
			min_dist = dist;
			beamid = pow(2,ib);
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


#endif
