// -*- C++ -*-

#include "TPCHit.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TMath.h>
#include <TSpectrum.h>

#include <std_ostream.hh>

//#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
//#include "DCParameters.hh"
//#include "DCTdcCalibMan.hh"
#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "DebugCounter.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "TPCLTrackHit.hh"
#include "TPCPadHelper.hh"
#include "TPCParamMan.hh"
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

//#define QuickAnalysis  1 // User EventSelectionTPCHits in RawData.cc
//#define FitPedestal    1
#define FitPedestal    0
//#define UseGaussian    1 //for gate noise fit
#define SingleFit      0
#define UseGaussian    0
#define UseLandau      0
#define UseExponential 0
#define UseGumbel      0
#define MultiFit       1
#define UseMultiGumbel 1
#define DebugEvDisp    0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
#endif

namespace
{
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gTPC  = TPCParamMan::GetInstance();
const Int_t MaxADC = 4096;
const Int_t MaxIteration = 3;
const Int_t MaxPeaks = 20;
const Double_t MaxChisqr = 1000.;
const Double_t zTgtTPC = -143.;
}

//_____________________________________________________________________________
TPCHit::TPCHit(TPCRawHit* rhit)
  : DCHit(rhit->LayerId(), rhit->RowId()),
    m_rhit(rhit),
    m_layer(),
    m_row(),
    m_pad(),
    m_pedestal(TMath::QuietNaN()),
    m_rms(TMath::QuietNaN()),
    m_de(),
    m_time(),
    m_chisqr(),
    m_cde(),
    m_ctime(),
    m_drift_length(),
    m_charge(),
    m_pos(),
    m_is_good(false),
    m_is_calculated(false),
    m_mrow(),
    m_tpc_flag(),
    m_hough_flag(),
    m_resx(),
    m_resy(),
    m_resz(),
    m_belong_track(false),
    m_hit_xz(),
    m_hit_yz()
{
  m_houghY_num.clear();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCHit::TPCHit(Int_t layer, Double_t mrow)
  : DCHit(layer, mrow),
    m_layer(layer),
    m_mrow(mrow)
{
  if(mrow-(int)mrow<0.5)
    m_row = (int)mrow;
  else
    m_row = 1+(int)mrow;
  //  m_row = (int)mrow;
  //  m_pad = tpc::GetPadId(layer, (int)mrow);
  m_pad = tpc::GetPadId(layer, m_row);
  m_charge = 0.;
  m_charge_center = 0.;
  m_pos=TVector3(0.,0.,0.);
  m_is_good = true;
  m_is_calculated = false;
  m_hough_flag = 0;
  m_houghY_num.clear();
  m_hit_xz = new DCHit(m_layer, m_row);
  m_hit_xz->SetWirePosition(m_pos.x());
  m_hit_xz->SetZ(m_pos.z());
  m_hit_xz->SetTiltAngle(0.);
  m_hit_xz->SetDummyPair();

  m_hit_yz = new DCHit(m_layer, m_row);
  m_hit_yz->SetWirePosition(m_pos.y());
  m_hit_yz->SetZ(m_pos.z());
  m_hit_yz->SetTiltAngle(90.);
  m_hit_yz->SetDummyPair();

  m_clsize = 0;
  debug::ObjectCounter::increase(ClassName());
}


//_____________________________________________________________________________
TPCHit::~TPCHit()
{
  ClearRegisteredHits();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCHit::AddDeTime(Double_t de, Double_t time)
{
  m_de.push_back(de);
  m_time.push_back(time);
}

//_____________________________________________________________________________
void
TPCHit::SetHoughYnum(Int_t houghY_num)
{
  m_houghY_num.push_back(houghY_num); 
}

//_____________________________________________________________________________
Bool_t
TPCHit::Calculate(Double_t clock)
{
  if(m_de.size() != m_time.size()){
    hddaq::cerr << FUNC_NAME << "found size mismatch: "
                << "m_de.size()=" << m_de.size() << ", "
                << "m_time.size()=" << m_time.size() << std::endl;
    return false;
  }

  if(!gTPC.IsReady()){
    hddaq::cerr << FUNC_NAME << " TPCParamMan must be initialized" << std::endl;
    return false;
  }

  for(Int_t i=0, n=m_de.size(); i<n; ++i){
    Double_t cde, ctime, dl;
    if(!gTPC.GetCDe(m_layer, m_row, m_de[i], cde)){
      hddaq::cerr << FUNC_NAME << " something is wrong at GetCDe("
                  << m_layer << ", " << m_row << ", " << m_de[i]
                  << ", " << cde << ")" << std::endl;
    }
    if(!gTPC.GetCTime(m_layer, m_row, m_time[i], ctime)){
      hddaq::cerr << FUNC_NAME << " something is wrong at GetCTime("
                  << m_layer << ", " << m_row << ", " << m_time[i]
                  << ", " << ctime << ")" << std::endl;
    }
    
    // Time Correction
    double c_clock = gTPC.GetC_Clock(m_layer, m_row, clock);
    //    ctime -= clock;
    //ctime -= c_clock;
    ctime += c_clock;
    //
    if(!gTPC.GetDriftLength(m_layer, m_row, ctime, dl)){
      hddaq::cerr << FUNC_NAME << " something is wrong at GetDriftLength("
                  << m_layer << ", " << m_row << ", " << ctime
                  << ", " << dl << ")" << std::endl;
    }
    m_cde.push_back(cde);
    m_ctime.push_back(ctime);
    m_drift_length.push_back(dl);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCHit::DoFit()
{
  static const Double_t MinDe = gUser.GetParameter("MinDeTPC");
  static const Double_t MinRms = gUser.GetParameter("MinRmsTPC");
  static const Int_t MinTimeBucket = gUser.GetParameter("TimeBucketTPC", 0);
  static const Int_t MaxTimeBucket = gUser.GetParameter("TimeBucketTPC", 1);

  if(m_is_calculated){
    hddaq::cerr << FUNC_NAME << " already calculated" << std::endl;
    return false;
  }

  if(!m_rhit){
    hddaq::cerr << FUNC_NAME << " m_rhit is nullptr" << std::endl;
    return false;
  }

#if DebugEvDisp
  //  gStyle->SetOptStat(1110);
  //gStyle->SetOptFit(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat(0);
  //  Option_t* option = "";
  Option_t* option = "W";//to avoid overflow events
#else
  //  Option_t* option = "Q";
  Option_t* option = "Q+W";//to avoid overflow events
#endif

  const auto level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  {
    Double_t mean = m_rhit->Mean();
    Double_t rms = m_rhit->RMS();
    m_pedestal = mean;
    m_rms = rms;
  }
// #if QuickAnalysis
//   {
//     Double_t max_adc = m_rhit->MaxAdc();
//     //std::cout<<"max_dE = "<<max_adc - m_pedestal<<std::endl;
//     if (max_adc-m_pedestal < MinDe) {
//       return false;
//     }
//   }
// #endif

  m_layer = m_rhit->LayerId();
  m_row = m_rhit->RowId();
  m_pad = tpc::GetPadId(m_layer, m_row);
  static TCanvas c1("c1", "c1", 800, 600);
  c1.cd();
  TH1D h1(FUNC_NAME+"-h1", "Pedestal",
           MaxADC, 0, MaxADC);
  TH1D h2(FUNC_NAME+"-h2",
           Form("Layer#%d Row#%d;Time Bucket;ADC",
                m_layer, m_row),
           NumOfTimeBucket, 0, NumOfTimeBucket);
  for(Int_t i=0, n=m_rhit->Fadc().size(); i<n; ++i){
    Double_t tb = i + 1;
    Double_t adc = m_rhit->Fadc().at(i);
    if(i < MinTimeBucket || MaxTimeBucket < i)
      continue;
    h1.Fill(adc);
    //avoid overflow
    if(adc<4075){
      h2.SetBinContent(tb, adc);
      //h2.SetBinError(tb, 1.);
    }
  }
  h2.GetXaxis()->SetRangeUser(MinTimeBucket, MaxTimeBucket);
  h2.GetYaxis()->SetRangeUser(m_pedestal-100., m_rhit->MaxAdc()+200);
  //std::cout<<"pedestal:"<<m_pedestal<<", max_adc:"<<m_rhit->MaxAdc()<<std::endl;

#if FitPedestal
  { //___ Pedestal Fitting
    Double_t mean = h1.GetXaxis()->GetBinCenter(h1.GetMaximumBin());
    Double_t rms = m_rhit->RMS();
    if(mean > 1000.) mean = 400.;
    TF1 f1("f1", "gaus", 0, MaxADC);
    f1.SetParameter(0, NumOfTimeBucket/rms);
    f1.SetParLimits(0, 1, NumOfTimeBucket);
    f1.SetParameter(1, mean);
    f1.SetParLimits(1, mean-3*rms, mean+3*rms);
    f1.SetParameter(2, rms);
    f1.SetParLimits(2, 0.2*rms, 1.*rms);
    for(UInt_t i=0; i<MaxIteration; ++i){
      h1.Fit("f1", option, "", mean - 2*rms, mean + 2*rms);
      mean = f1.GetParameter(1);
      rms = f1.GetParameter(2);
    }
    Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
    if(chisqr < 100.){
      m_pedestal = mean;
      // m_rms = rms;
    } else {
      // hddaq::cerr << FUNC_NAME << " failed to fit pedestal" << std::endl;
      // Print();
    }
#if DebugEvDisp
    hddaq::cout << FUNC_NAME << " (mean,rms,chisqr)=("
		<< mean << "," << rms << "," << chisqr << ")" << std::endl;
    h1.GetXaxis()->SetRangeUser(mean-10*rms, mean+10*rms);
    h1.Draw();
    c1.Modified();
    c1.Update();
    // getchar();
    Double_t max_adc = m_rhit->MaxAdc();
    if(max_adc - mean < MinDe){
      return false;
    }
#endif
  }
#endif

  if(m_rms < MinRms){
    return false;
  }


#if 1
  { //___ Peak Search
    Double_t sigma = 3;
    Double_t threshold = 0.4;
    TSpectrum spec(MaxPeaks);
    const Int_t n_peaks = spec.Search(&h2, sigma, "", threshold);
    if(n_peaks == 0)
      return false;
    Double_t* x_peaks = spec.GetPositionX();
    std::vector<Double_t> sorted_x_peaks(n_peaks);
    sorted_x_peaks.assign(x_peaks, x_peaks+n_peaks);
    std::stable_sort(sorted_x_peaks.begin(), sorted_x_peaks.end());

#if SingleFit
    for(Int_t i=0; i<n_peaks; ++i){
      Double_t time = sorted_x_peaks[i];
      Double_t de = h2.GetBinContent(h2.FindBin(time)) - m_pedestal;
      // hddaq::cout << FUNC_NAME << " time,de=" << time << "," << de << std::endl;
#if UseGaussian
      TF1 f1("f1", "gaus(0)+pol1(3)", 0, NumOfTimeBucket);
      sigma = 3.;
      f1.SetParameter(0, de);
      f1.SetParLimits(0, de*0.8, de*1.2);
      f1.SetParameter(1, time);
      f1.SetParLimits(1, time-2*sigma, time+2*sigma);
      f1.SetParameter(2, sigma);
      f1.SetParLimits(2, 2., 5.);
      // f1.FixParameter(2, 1.0);
      f1.SetParameter(3, m_pedestal);
      f1.SetParLimits(3, m_pedestal-3*m_rms, m_pedestal+3*m_rms);
      f1.FixParameter(4, 0);
      // f1.SetParLimits(4, -0.01, 0.01);
      for(UInt_t j=0; j<MaxIteration; ++j){
        h2.Fit("f1", option, "", time-2.5*sigma, time+2.5*sigma);
        time = f1.GetParameter(1);
        sigma = f1.GetParameter(2);
      }
      auto p = f1.GetParameters();
      time = f1.GetMaximumX();
      de = f1.GetMaximum() - p[3];
      // de = p[0]*p[2]*TMath::Sqrt(2*TMath::Pi());
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseLandau
      TF1 f1("f1", "landau(0)+pol1(3)", 0, NumOfTimeBucket);
      sigma = 1.;
      f1.SetParameter(0, de/0.18);
      f1.SetParLimits(0, de/0.18*0.8, de/0.18*1.2);
      f1.SetParameter(1, time);
      f1.SetParLimits(1, time-2*sigma, time+2*sigma);
      f1.SetParameter(2, sigma);
      f1.SetParLimits(2, 0.5, 1.5);
      // f1.FixParameter(2, 1.0);
      f1.SetParameter(3, m_pedestal);
      f1.SetParLimits(3, m_pedestal-3*m_rms, m_pedestal+3*m_rms);
      f1.FixParameter(4, 0);
      // f1.SetParLimits(4, -0.01, 0.01);
      for(UInt_t j=0; j<MaxIteration; ++j){
        h2.Fit("f1", option, "", time-3*sigma, time+7*sigma);
        time = f1.GetParameter(1);
        sigma = f1.GetParameter(2);
      }
      auto p = f1.GetParameters();
      time = f1.GetMaximumX();
      de = f1.GetMaximum() - p[3]; //p[0]*p[2]*TMath::Sqrt(2*TMath::Pi());
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseExponential
      TF1 f1("f1", "[0]*(x-[1])*exp(-(x-[1])/[2])+[3]", 0, NumOfTimeBucket);
      Double_t c = 3.; // decay constant
      f1.SetParameter(0, de*TMath::Exp(1)/c);
      f1.SetParLimits(0, TMath::Exp(1)/c, m_rhit->MaxAdc()*TMath::Exp(1)/c);
      f1.SetParameter(1, time+c);
      f1.SetParLimits(1, time+c-12, time+c+12);
      f1.SetParameter(2, c);
      f1.SetParLimits(2, 1.5, 7.5);
      f1.FixParameter(3, m_pedestal);
      // f1.SetParameter(3, m_pedestal);
      // f1.SetParLimits(3, m_pedestal - m_rms, m_pedestal + m_rms);
      h2.Fit("f1", option, "", time-5, time+10);
      auto p = f1.GetParameters();
      time = p[1] + p[2]; // peak time
      // c = p[2];
      de = p[0]*p[2]*TMath::Exp(-1); // amplitude
      sigma = p[2];
      // de = p[0]*p[2]*p[2]; // integral
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseGumbel
      TF1 f1("f1", "[0]*exp(-(x-[1])/[2])*exp(-exp(-(x-[1])/[2]))+[3]",
             0, NumOfTimeBucket);
      f1.SetParameter(0, de);
      //f1.SetParLimits(0, 0, m_rhit->MaxAdc()+m_pedestal+m_rms);
      //f1.SetParLimits(0, 0, de+m_rms);
      f1.SetParameter(1, time);
      f1.SetParLimits(1, time-10, time+10);
      f1.SetParameter(2, 3);
      f1.SetParLimits(2, 1., 10.);
      f1.SetParameter(3, m_pedestal);
      f1.SetParLimits(3, m_pedestal - 0.5*m_rms, m_pedestal + 0.5*m_rms);
      //h2.Fit("f1", option, "", time-6, time+10);
      //h2.Fit("f1", option, "",MinTimeBucket , MaxTimeBucket);
      h2.Fit("f1", option, "", time-6.-de*0.07, time+10.+de*0.11);
      auto p = f1.GetParameters();
      time = p[1]; // peak time
      de = f1.Eval(time) - p[3]; // amplitude
      sigma = p[2];
      Double_t chisqr = f1.GetChisquare()/f1.GetNDF();
      //std::cout<<"chisqr:"<<chisqr<<std::endl;
#endif
      if (chisqr > MaxChisqr || de < MinDe)
        continue;
      m_time.push_back(time);
      m_de.push_back(de);
      m_chisqr.push_back(chisqr);
      m_sigma.push_back(sigma);
#if DebugEvDisp
      hddaq::cout << FUNC_NAME << " (time,de,chisqr, rms)=("
		  << time << "," << de << "," << chisqr << "," << m_rms <<")" << std::endl;
#endif
    }
#endif

#if MultiFit
#if UseMultiGumbel
    Double_t time[n_peaks];
    Double_t de[n_peaks];
    TString func;
    for(Int_t i=0; i<n_peaks; ++i){
      time[i] = sorted_x_peaks[i];
      de[i] = h2.GetBinContent(h2.FindBin(time[i])) - m_pedestal;
      if(i!=0)
	func += "+";
      func += "["+std::to_string(i*3)+"]";
      func += "*exp(-(x-["+std::to_string(i*3+1)+"])";
      func += "/["+std::to_string(i*3+2)+"])";
      func += "*exp(-exp(-(x-["+std::to_string(i*3+1)+"])/["+std::to_string(i*3+2)+"]))";
    }
    func += "+["+std::to_string(n_peaks*3)+"]";
    TF1 f1("f1", func,
	   0, NumOfTimeBucket);
    for(Int_t i=0; i<n_peaks; ++i){
      f1.SetParameter(i*3, de[i]);
      f1.SetParameter(i*3+1, time[i]);
      f1.SetParLimits(i*3+1, time[i]-10, time[i]+10);
      f1.SetParameter(i*3+2, 3);
      f1.SetParLimits(i*3+2, 1., 10.);
      //f1.SetParLimits(i*3+2, 2.5, 10.);
    }
    f1.SetParameter(n_peaks*3, m_pedestal);
    f1.SetParLimits(n_peaks*3, m_pedestal-0.5*m_rms, m_pedestal+0.5*m_rms);
    double FitRangeMin = MinTimeBucket;
    double FitRangeMax = MaxTimeBucket;
    if(MinTimeBucket<time[0]-6.-de[0]*0.07)
      FitRangeMin = time[0]-6.-de[0]*0.07 ;
    if(MaxTimeBucket>time[n_peaks-1]+10.+de[n_peaks-1]*0.11)
      FitRangeMax = time[n_peaks-1]+10.+de[n_peaks-1]*0.11;
    h2.Fit("f1", option, "", FitRangeMin, FitRangeMax);
    auto p = f1.GetParameters();
    Double_t chisqr = f1.GetChisquare()/f1.GetNDF();
    if(chisqr < MaxChisqr){
      for(Int_t i=0; i<n_peaks; ++i){
	time[i] = p[i*3+1]; // peak time
	de[i] = f1.Eval(time[i]) - p[n_peaks*3]; // amplitude
	sigma = p[i*3+2];
	if(de[i] > MinDe && sigma > 2.5){
	  m_time.push_back(time[i]);
	  m_de.push_back(de[i]);
	  m_chisqr.push_back(chisqr);
	  m_sigma.push_back(sigma);
	}
      }
    }
#endif
#endif
  }
#endif

#if DebugEvDisp
  //h2.Draw("L");
  // double maxde = 0.;
  // if(m_de.size()>0)
  //   maxde = TMath::MaxElement(m_de.size(), m_de.data());
  // if(m_time.size()>3&&maxde>150.){
  //h2.Draw("");
  c1.Modified();
  c1.Update();
  gSystem->ProcessEvents();
  Print();
  c1.Print("c1.pdf");
  getchar();
  //  }
#endif

  gErrorIgnoreLevel = level;
  m_is_good = true;
  return true;
}

//_____________________________________________________________________________
void
TPCHit::ClearRegisteredHits()
{
  del::ClearContainer(m_register_container);
  if(m_hit_xz) delete m_hit_xz;
  if(m_hit_yz) delete m_hit_yz;
}

//_____________________________________________________________________________
double
TPCHit::GetResolutionX()
{
  //calculated by using NIM paper
  double y_pos= m_pos.Y();
//double s0 = 0.204;// mm HIMAC result //To do parameter
  double s0 = gUser.GetParameter("TPC_sigma0");
  //s0 is considered for common resolution
//  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double Dt = gUser.GetParameter("TPC_Dt");
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1
  double e_ALD = exp(-1.*A*L_D);
  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  //double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT_r = sqrt(sT2);
  if(m_clsize ==1){
    double Rpad = tpc::padParameter[m_layer][2];
    double padsize = Rpad*2.*acos(-1)/tpc::padParameter[m_layer-1][3];
    sT_r = padsize/sqrt(12.);
  }
  double sT_padlen = tpc::padParameter[m_layer][5]/sqrt(12.);

  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);

  double res_x = sqrt(pow(sT_r*cos(alpha),2)+pow(sT_padlen*sin(alpha),2));

  return res_x;
}

//_____________________________________________________________________________
double
TPCHit::GetResolutionZ()
{
  //calculated by using NIM paper
  double y_pos= m_pos.Y();
//double s0 = 0.204;// mm HIMAC result //To do parameter
  double s0 = gUser.GetParameter("TPC_sigma0");
  //s0 is considered for common resolution
//  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double Dt = gUser.GetParameter("TPC_Dt");
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1
  double e_ALD = exp(-1.*A*L_D);
  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  //double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT_r = sqrt(sT2);
  if(m_clsize ==1){
    double Rpad = tpc::padParameter[m_layer][2];
    double padsize = Rpad*2.*acos(-1)/tpc::padParameter[m_layer-1][3];
    sT_r = padsize/sqrt(12.);
  }
  double sT_padlen = tpc::padParameter[m_layer][5]/sqrt(12.);

  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);

  double res_z = sqrt(pow(sT_r*sin(alpha),2)+pow(sT_padlen*cos(alpha),2));

  return res_z;
}

//_____________________________________________________________________________
double
TPCHit::GetResolutionY()
{
  // temporary
  return 0.5;
  //  return 1.5;
  //return 2.;
}

//_____________________________________________________________________________
double
TPCHit::GetResolution()
{
  double res_x = GetResolutionX();
  double res_y = GetResolutionY();
  double res_z = GetResolutionZ();

  double tot_res = sqrt(res_x*res_x + res_y*res_y + res_z*res_z);
  return tot_res;
}

//_____________________________________________________________________________
void
TPCHit::Print(const std::string& arg, std::ostream& ost) const
{
  const int w = 16;
  ost << "#D " << FUNC_NAME << " " << arg << std::endl
      << std::setw(w) << std::left << "layer" << m_layer << std::endl
      << std::setw(w) << std::left << "row"  << m_row  << std::endl
      << std::setw(w) << std::left << "pad"  << m_pad  << std::endl
      << std::setw(w) << std::left << "charge" << m_charge << std::endl
      << std::setw(w) << std::left << "posx" << m_pos.x() << std::endl
      << std::setw(w) << std::left << "posy" << m_pos.y() << std::endl
      << std::setw(w) << std::left << "posz" << m_pos.z() << std::endl;
  for(Int_t i=0, n=GetNHits(); i<n; ++i){
    ost << std::setw(3) << std::right << i << " "
	<< "(time,de,chisqr,ctime,cde,dl)=("
        << std::fixed << std::setprecision(3) << std::setw(9) << m_time[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_de[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_chisqr[i]
      //<< std::fixed << std::setprecision(3) << std::setw(9) << m_ctime[i]
      //<< std::fixed << std::setprecision(3) << std::setw(9) << m_cde[i]
        // << std::fixed << std::setprecision(3) << std::setw(9) << m_drift_length[i]
        << ")" << std::endl;
  }
}
