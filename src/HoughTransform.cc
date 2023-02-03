// -*- C++ -*-

#include "TPCTrackSearch.hh"

#include <chrono>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TH2D.h>
#include <TH3D.h>

#include "DebugTimer.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "DeleteUtility.hh"
#include "ConfMan.hh"
#include "TPCPadHelper.hh"
#include "TPCCluster.hh"

#include "RootHelper.hh"

#define DebugDisp      0
#define DebugEvDisp    0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
static Int_t cannum = 0;
#endif
namespace
{
  const auto qnan = TMath::QuietNaN();
  const double& HS_field_0 = 0.9860;
  const double& valueHSCalc = ConfMan::Get<Double_t>("HSFLDCALC");
  const double& valueHSHall = ConfMan::Get<Double_t>("HSFLDHALL");
  const double Const = 0.299792458;

  // HoughXZ binning
  static const Int_t    Li_theta_ndiv = 180;
  static const Double_t Li_theta_min  =   0;
  static const Double_t Li_theta_max  = 180;
  static const Int_t    Li_rho_ndiv =  180;
  static const Double_t Li_rho_min  = -720;
  static const Double_t Li_rho_max  =  720;

  // HoughCircle binning
  const Int_t    nBin_rdiff = 11;
  const Double_t rdiff_min = -110.;
  const Double_t rdiff_max = 110.;
  const Int_t    nBin_theta = 720;
  const Double_t theta_min = -1.*acos(-1);
  const Double_t theta_max = acos(-1);
  const Int_t    nBin_p = 900;
  const Double_t pmin = 50.;
  const Double_t pmax = 1550.; //MeV/c

  // HoughY binning
  //const Ynt_t    thetaY_ndiv = 360; //previous value
  const Int_t    thetaY_ndiv = 1080;
  const Double_t thetaY_min  =   0.;
  const Double_t thetaY_max  = 180.;
  const Int_t    r_ndiv =  2000;
  const Double_t r_min  = -5000.;
  const Double_t r_max  =  5000.;

  Double_t GetMagneticField()
  {
    return HS_field_0*(valueHSHall/valueHSCalc);
  }
}

//_____________________________________________________________________________
namespace tpc
{

//_____________________________________________________________________________
Bool_t
HoughTransformLineXZ(std::vector<TVector3> gHitPos, Int_t *MaxBin,
		     Double_t *LinearPar, Int_t MinNumOfHits /*=8*/)
{

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif
  Bool_t status = true;
  TH2D *histXZ = new TH2D("hist_xz", "; #theta (deg.); #rho (mm)",
			   Li_theta_ndiv, Li_theta_min, Li_theta_max,
			   Li_rho_ndiv, Li_rho_min, Li_rho_max);

  //Hough space
  //for TPC linear track
  //rho = z*cos(theta) + x*sin(theta)
  //x = p[0] + p[1]*z

  for(Int_t i=0; i<gHitPos.size(); ++i){
    for(int ti=0; ti<thetaY_ndiv; ti++){
      Double_t theta = histXZ->GetXaxis()->GetBinCenter(ti+1);
      Double_t rho = TMath::Cos(theta*TMath::DegToRad())*gHitPos[i].Z()
	+ TMath::Sin(theta*TMath::DegToRad())*gHitPos[i].X();
      histXZ->Fill(theta, rho);
    } //ti
  } //i

#if DebugEvDisp
  histXZ->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = histXZ->GetMaximumBin();
  Int_t mx, my, mz;
  histXZ->GetBinXYZ(maxbin, mx, my, mz);
  MaxBin[0] = mx; MaxBin[1] = my; MaxBin[2] = mz;

  Double_t mtheta = histXZ->GetXaxis()->GetBinCenter(mx)*TMath::DegToRad();
  Double_t mr = histXZ->GetYaxis()->GetBinCenter(my);

  LinearPar[0] = mr/TMath::Sin(mtheta);
  LinearPar[1]  = -TMath::Cos(mtheta)/TMath::Sin(mtheta);

#if DebugDisp
  std::cout<<FUNC_NAME+" (Z, X) hough vote : "<<histXZ->GetBinContent(maxbin)<<std::endl;
  std::cout<<FUNC_NAME+" Hough (theta, r)=("<<mtheta<<", "<<mr<<")"<<std::endl;
#endif

  if(histXZ->GetMaximum() < 0.5*MinNumOfHits) status = false;

  delete histXZ;
  return status;

}

//_____________________________________________________________________________
Bool_t
HoughTransformLineXZ(const std::vector<TPCClusterContainer>& ClCont,
		     Int_t *MaxBin, Double_t *LinearPar,
		     Int_t MinNumOfHits /*=8*/)
{

  std::vector<TVector3> gHitPos;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      gHitPos.push_back(pos);
    } //ci
  } //layer

  return HoughTransformLineXZ(gHitPos, MaxBin, LinearPar, MinNumOfHits);
}

//_____________________________________________________________________________
void
HoughTransformLineYPhi(std::vector<TVector3> gHitPos,
		       Double_t *HelixPar, Double_t MaxHoughWindow)
{

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  TH2D *histY = new TH2D("histY",";theta (deg.);r (mm)",
			 thetaY_ndiv, thetaY_min, thetaY_max,
			 r_ndiv, r_min, r_max);

  //Hough space
  //mr = tmp_xval*cos(mtheta) + tmpz*sin(mtheta)
  //tmpz = p[2] + p[4]*tmp_xval
  for(Int_t i=0; i<gHitPos.size(); ++i){
    Double_t tmpx = -gHitPos[i].X();
    Double_t tmpy = gHitPos[i].Z() - tpc::ZTarget;
    Double_t tmpz = gHitPos[i].Y();
    Double_t r_cal = sqrt(pow(tmpx - HelixPar[0], 2) + pow(tmpy - HelixPar[1], 2));
    Double_t dist = fabs(r_cal - HelixPar[3]);
    if(dist < MaxHoughWindow){
      for(int ti=0; ti<thetaY_ndiv; ti++){
	Double_t theta = histY->GetXaxis()->GetBinCenter(ti+1);
	Double_t tmpt = atan2(tmpy - HelixPar[1], tmpx - HelixPar[0]);
	Double_t tmp_xval = HelixPar[3]*tmpt;
	histY->Fill(theta, cos(theta*acos(-1.)/180.)*tmp_xval
		    +sin(theta*acos(-1.)/180.)*tmpz);
      } //ti
    } //dist
  }

#if DebugEvDisp
  histY->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = histY->GetMaximumBin();
  Int_t mx, my, mz;
  histY->GetBinXYZ(maxbin, mx, my, mz);

  double mtheta = histY->GetXaxis()->GetBinCenter(mx)*acos(-1.)/180.;
  Double_t mr = histY->GetYaxis()->GetBinCenter(my);
  HelixPar[2] = mr/sin(mtheta);
  HelixPar[4] = -cos(mtheta)/sin(mtheta);

#if DebugDisp
  std::cout<<FUNC_NAME+" (Y, phi) hough vote : "<<histY->GetBinContent(maxbin)<<std::endl;
  std::cout<<FUNC_NAME+" Hough (theta, r)=("<<mtheta<<", "<<mr<<")"<<std::endl;
#endif
  delete histY;

}

//_____________________________________________________________________________
void
HoughTransformLineYPhi(const std::vector<TPCClusterContainer>& ClCont,
		       Double_t *HelixPar, Double_t MaxHoughWindow)
{

  std::vector<TVector3> gHitPos;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      gHitPos.push_back(pos);
    } //ci
  } //layer

  HoughTransformLineYPhi(gHitPos, HelixPar, MaxHoughWindow);
}

//_____________________________________________________________________________
Bool_t
HoughTransformCircleXZ(std::vector<TVector3> gHitPos,
		       Int_t *MaxBin, Double_t *HelixPar,
		       Int_t MinNumOfHits /*=8*/)
{

  // Equation
  // (x - (r + rd)*cos(theta))^2 + (y - (r + rd)*sin(theta))^2 = r^2
  // p = r * Const * dMagneticField;
  Double_t dMagneticField = GetMagneticField();
  Bool_t status = true;

#if DebugEvDisp
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 900, 900);
  c1.cd();
#endif

  //for TPC circle track
  //Hough-transform
  TH3D *hist_circle = new TH3D("hist_circle",";rd (mm); theta (rad); p(MeV/c)",
			       nBin_rdiff, rdiff_min,  rdiff_max,
			       nBin_theta, theta_min, theta_max,
			       nBin_p, pmin, pmax);

  for(Int_t i=0; i<gHitPos.size(); ++i){
    Double_t x = -gHitPos[i].X();
    Double_t y = gHitPos[i].Z() - tpc::ZTarget;
    for(Int_t ird=0; ird<nBin_rdiff; ++ird){
      Double_t rd = hist_circle->GetXaxis()->GetBinCenter(ird+1);
      for(Int_t ip=0; ip<nBin_p; ++ip){
	Double_t p = hist_circle->GetZaxis()->GetBinCenter(ip+1);
	Double_t r = p/(Const*dMagneticField);

	//a*sin(theta) + b*cos(theta) + c = 0
	Double_t a = 2.*(r+rd)*y;
	Double_t b = 2.*(r+rd)*x;
	Double_t c = -1.*(rd*rd + 2.*r*rd + x*x + y*y);

	Double_t r0 = sqrt(a*a + b*b);
	if(fabs(-1.*c/r0)>1.) continue;
	Double_t theta1_alpha = asin(-1.*c/r0);

	Double_t theta2_alpha;
	if(theta1_alpha>0.) theta2_alpha = acos(-1.) - theta1_alpha;
	else theta2_alpha = -1.*acos(-1.) - theta1_alpha;

	Double_t theta_alpha = atan2(b, a);

	Double_t xcenter1 = (r+rd)*cos(theta1_alpha - theta_alpha);
	Double_t ycenter1 = (r+rd)*sin(theta1_alpha - theta_alpha);
	Double_t r_re1 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

	Double_t xcenter2 = (r+rd)*cos(theta2_alpha - theta_alpha);
	Double_t ycenter2 = (r+rd)*sin(theta2_alpha - theta_alpha);
	Double_t r_re2 = sqrt(pow(x-xcenter1,2) + pow(y-ycenter1,2));

	Double_t theta1 = atan2(ycenter1, xcenter1);
	Double_t theta2 = atan2(ycenter2, xcenter2);
	if(TMath::IsNaN(theta1)){
	  std::cout<<"theta1="<<theta1<<", x="<<x<<", y="<<y
		   <<"rd="<<rd<<", r"<<r<<std::endl;
	}
	if(fabs(r-r_re1)>0.01||fabs(r-r_re2)>0.01){
	  std::cout<<"r="<<r<<", r_re1="<<r_re1<<", r_re1="<<r_re2<<std::endl;
	  std::cout<<"x:"<<x<<", y:"<<y
		   <<", theta1:"<<theta1<<", theta2:"<<theta2
		   <<", theta_alpha:"<<theta_alpha<<std::endl;
	}
	hist_circle->Fill(rd, theta1, p);
	hist_circle->Fill(rd, theta2, p);
      } // p bins
    } // rdiff bins
  } //i

#if DebugEvDisp
  hist_circle->Draw("colz");
  gPad->Modified();
  gPad->Update();
  c1.Print(Form("c%d.pdf",cannum));
  getchar();
  cannum++;
#endif

  Int_t maxbin = hist_circle->GetMaximumBin();
  Int_t mx, my, mz;
  hist_circle->GetBinXYZ(maxbin, mx, my, mz);
  MaxBin[0] = mx; MaxBin[1] = my; MaxBin[2] = mz;

  Double_t hough_rd = hist_circle->GetXaxis()->GetBinCenter(mx);
  Double_t hough_theta = hist_circle->GetYaxis()->GetBinCenter(my);
  Double_t hough_p = hist_circle->GetZaxis()->GetBinCenter(mz);
  HelixPar[3] = hough_p/(Const*dMagneticField); //helix r
  HelixPar[0] = (HelixPar[3] + hough_rd)*cos(hough_theta); //helix cx
  HelixPar[1] = (HelixPar[3] + hough_rd)*sin(hough_theta); //helix cy

#if DebugDisp
  std::cout<<FUNC_NAME+" XZ hough vote : "<<hist_circle->GetMaximum()<<std::endl;
  std::cout<<FUNC_NAME+" Circle hough maxbin : "<<hist_circle->GetBinContent(maxbin)<<std::endl;
  std::cout<<FUNC_NAME+" Hough(rd, theta, p)=("<<hough_rd<<", "
	   <<hough_theta<<", "<<hough_p<<")"<<std::endl;
#endif

  if(hist_circle->GetMaximum() < 0.5*MinNumOfHits) status = false;

  delete hist_circle;
  return status;
}

//_____________________________________________________________________________
Bool_t
HoughTransformCircleXZ(const std::vector<TPCClusterContainer>& ClCont,
		       Int_t *MaxBin, Double_t *HelixPar,
		       Int_t MinNumOfHits /*=8*/)
{

  std::vector<TVector3> gHitPos;
  for(Int_t layer=0; layer<NumOfLayersTPC; layer++){
    for(Int_t ci=0, n=ClCont[layer].size(); ci<n; ci++){
      auto cl = ClCont[layer][ci];
      TPCHit* hit = cl->GetMeanHit();
      if(hit->GetHoughFlag()>0) continue;
      TVector3 pos = cl->GetPosition();
      gHitPos.push_back(pos);
    } //ci
  } //layer

  return HoughTransformCircleXZ(gHitPos, MaxBin, HelixPar, MinNumOfHits);
}

} //namespace tpc
