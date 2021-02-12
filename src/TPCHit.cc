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
#include "TPCRawHit.hh"
#include "UserParamMan.hh"

#define QuickAnalysis  1
//#define FitPedestal    1
#define FitPedestal    0
#define UseGaussian    0
#define UseLandau      0
#define UseExponential 1
#define DebugEvDisp    0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
namespace { TApplication app( "DebugApp", nullptr, nullptr ); }
#endif

namespace
{
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const Int_t MaxADC = 4096;
const Int_t MaxIteration = 3;
const Int_t MaxPeaks = 20;
const Double_t MaxChisqr = 1000.;
const Double_t zTgtTPC = -143.;
}

//_____________________________________________________________________________
TPCHit::TPCHit( TPCRawHit* rhit )
  : DCHit( rhit->LayerId(), rhit->RowId() ),
    m_rhit( rhit ),
    m_row(),
    m_pad(),
    m_pedestal( TMath::QuietNaN() ),
    m_rms( TMath::QuietNaN() ),
    m_de(),
    m_time(),
    m_chisqr(),
    m_charge(),
    m_pos(),
    m_is_good( false ),
    m_is_calculated( false ),
    m_mrow(),
    m_tpc_flag(),
    m_resx(),
    m_resy(),
    m_resz(),
    m_belong_track( false ),
    m_hit_xz(),
    m_hit_yz()
{
  // m_pos = tpc::getPosition(padid);
  // m_pos.SetY(y);
  // m_layer = tpc::getLayerID(padid);
  // m_row   = tpc::getRowID(padid);
  // m_hit_xz = new DCHit(m_layer, m_row);
  // m_hit_xz->SetWirePosition(m_pos.x());
  // m_hit_xz->SetZ(m_pos.z());
  // m_hit_xz->SetTiltAngle(0.);
  // m_hit_xz->SetDummyPair();

  // m_hit_yz = new DCHit(m_layer, m_row);
  // m_hit_yz->SetWirePosition(m_pos.y());
  // m_hit_yz->SetZ(m_pos.z());
  // m_hit_yz->SetTiltAngle(90.);
  // m_hit_yz->SetDummyPair();

  debug::ObjectCounter::increase( ClassName() );
}

//_____________________________________________________________________________
TPCHit::TPCHit( int padid, TVector3 pos, double charge )
  : m_pad(padid),
    m_charge(charge)
{
  m_pos=pos;
  m_layer = tpc::getLayerID(padid);
  m_row   = tpc::getRowID(padid);
  m_is_good = true;
  m_is_calculated = false;
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

  debug::ObjectCounter::increase( ClassName() );
}


//______________________________________________________________________________
TPCHit::~TPCHit( void )
{
  ClearRegisteredHits();
  debug::ObjectCounter::decrease( ClassName() );
}

//______________________________________________________________________________
void
TPCHit::ClearRegisteredHits( void )
{
  del::ClearContainer( m_register_container );
  if( m_hit_xz ) delete m_hit_xz;
  if( m_hit_yz ) delete m_hit_yz;
}

//______________________________________________________________________________
Bool_t
TPCHit::CalcTPCObservables( void )
{
  static const Double_t MinDeTPC = gUser.GetParameter( "MinDeTPC" );
  static const Double_t MinRmsTPC = gUser.GetParameter( "MinRmsTPC" );
  static const Double_t MinTimeBucketTPC
    = gUser.GetParameter( "TimeBucketTPC", 0 );
  static const Double_t MaxTimeBucketTPC
    = gUser.GetParameter( "TimeBucketTPC", 1 );

  if( m_is_calculated ){
    hddaq::cerr << FUNC_NAME << " already calculated" << std::endl;
    return false;
  }

#if DebugEvDisp
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1);
  Option_t* option = "";
#else
  Option_t* option = "Q";
#endif

  const auto level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  {
    Double_t mean = m_rhit->Mean();
    Double_t rms = m_rhit->RMS();
    m_pedestal = mean;
    m_rms = rms;
  }
#if QuickAnalysis
  {
    Double_t max_adc = m_rhit->MaxAdc();
    if(max_adc-m_pedestal<MinDeTPC)
      return false;
  }
#endif

  static TCanvas c1;
  c1.cd();
  TH1D h1( ClassName()+"::h1", "h1", MaxADC, 0, MaxADC );
  TH1D h2( ClassName()+"::h2", "h2", NumOfTimeBucket, 0, NumOfTimeBucket );
  for( Int_t i=0, n=m_rhit->Fadc().size(); i<n; ++i ){
    Double_t tb = i + 1;
    Double_t adc = m_rhit->Fadc().at( i );
    if( i < MinTimeBucketTPC || MaxTimeBucketTPC < i )
      continue;
    h1.Fill( adc );
    h2.SetBinContent( tb, adc );
  }

#if FitPedestal
  { //___ Pedestal Fitting
    Double_t mean = h1.GetXaxis()->GetBinCenter( h1.GetMaximumBin() );
    Double_t rms = m_rhit->RMS();
    if( mean > 1000. ) mean = 400.;
    TF1 f1("f1", "gaus", 0, MaxADC);
    f1.SetParameter( 0, NumOfTimeBucket/rms );
    f1.SetParLimits( 0, 1, NumOfTimeBucket );
    f1.SetParameter( 1, mean );
    f1.SetParLimits( 1, mean-3*rms, mean+3*rms );
    f1.SetParameter( 2, rms );
    f1.SetParLimits( 2, 0.2*rms, 1.*rms );
    for( UInt_t i=0; i<MaxIteration; ++i ){
      h1.Fit("f1", option, "", mean - 2*rms, mean + 2*rms );
      mean = f1.GetParameter(1);
      rms = f1.GetParameter(2);
    }
    Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
    if( chisqr < 100. ){
      m_pedestal = mean;
      // m_rms = rms;
    } else {
      // hddaq::cerr << FUNC_NAME << " failed to fit pedestal" << std::endl;
      // Print();
    }
#if DebugEvDisp
    hddaq::cout << FUNC_NAME << " (mean,rms,chisqr)=("
		<< mean << "," << rms << "," << chisqr << ")" << std::endl;
    h1.GetXaxis()->SetRangeUser( mean-10*rms, mean+10*rms );
    h1.Draw();
    c1.Modified();
    c1.Update();
    gSystem->ProcessEvents();
    // getchar();
    Double_t max_adc = m_rhit->MaxAdc();
    if( max_adc - mean < 40 )
      return false;
#endif
  }
#endif

  if( m_rms < MinRmsTPC ){
    return false;
  }

#if 1
  { //___ Peak Search
    Double_t sigma = 3;
    Double_t threshold = 0.4;
    TSpectrum spec( MaxPeaks );
    Int_t n_peaks = spec.Search( &h2, sigma, "", threshold );
    if( n_peaks == 0 )
      return false;
    Double_t* x_peaks = spec.GetPositionX();
    std::vector<Double_t> sorted_x_peaks( n_peaks );
    sorted_x_peaks.assign( x_peaks, x_peaks+n_peaks );
    std::stable_sort( sorted_x_peaks.begin(), sorted_x_peaks.end() );
    for( Int_t i=0; i<n_peaks; ++i ){
      Double_t time = sorted_x_peaks[i];
      Double_t de = h2.GetBinContent( h2.FindBin( time ) ) - m_pedestal;
      // hddaq::cout << FUNC_NAME << " time,de=" << time << "," << de << std::endl;
#if UseGaussian
      TF1 f1( "f1", "gaus(0)+pol1(3)", 0, NumOfTimeBucket );
      sigma = 3.;
      f1.SetParameter( 0, de );
      f1.SetParLimits( 0, de*0.8, de*1.2 );
      f1.SetParameter( 1, time );
      f1.SetParLimits( 1, time-2*sigma, time+2*sigma );
      f1.SetParameter( 2, sigma );
      f1.SetParLimits( 2, 2., 5. );
      // f1.FixParameter( 2, 1.0 );
      f1.SetParameter( 3, m_pedestal );
      f1.SetParLimits( 3, m_pedestal-3*m_rms, m_pedestal+3*m_rms );
      f1.FixParameter( 4, 0 );
      // f1.SetParLimits( 4, -0.01, 0.01 );
      for( UInt_t j=0; j<MaxIteration; ++j ){
        h2.Fit( "f1", option, "", time-2.5*sigma, time+2.5*sigma );
        time = f1.GetParameter(1);
        sigma = f1.GetParameter(2);
      }
      auto p = f1.GetParameters();
      time = f1.GetMaximumX();
      de = f1.GetMaximum() - p[3];
      // de = p[0]*p[2]*TMath::Sqrt( 2*TMath::Pi() );
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseLandau
      TF1 f1( "f1", "landau(0)+pol1(3)", 0, NumOfTimeBucket );
      sigma = 1.;
      f1.SetParameter( 0, de/0.18 );
      f1.SetParLimits( 0, de/0.18*0.8, de/0.18*1.2 );
      f1.SetParameter( 1, time );
      f1.SetParLimits( 1, time-2*sigma, time+2*sigma );
      f1.SetParameter( 2, sigma );
      f1.SetParLimits( 2, 0.5, 1.5 );
      // f1.FixParameter( 2, 1.0 );
      f1.SetParameter( 3, m_pedestal );
      f1.SetParLimits( 3, m_pedestal-3*m_rms, m_pedestal+3*m_rms );
      f1.FixParameter( 4, 0 );
      // f1.SetParLimits( 4, -0.01, 0.01 );
      for( UInt_t j=0; j<MaxIteration; ++j ){
        h2.Fit( "f1", option, "", time-3*sigma, time+7*sigma );
        time = f1.GetParameter(1);
        sigma = f1.GetParameter(2);
      }
      auto p = f1.GetParameters();
      time = f1.GetMaximumX();
      de = f1.GetMaximum() - p[3]; //p[0]*p[2]*TMath::Sqrt( 2*TMath::Pi() );
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#elif UseExponential
      TF1 f1( "f1", "[0]*(x-[1])*exp(-(x-[1])/[2])+[3]", 0, NumOfTimeBucket );
      Double_t c = 3.; // decay constant
      f1.SetParameter( 0, de*TMath::Exp(1)/c );
      f1.SetParLimits( 0, TMath::Exp(1)/c, m_rhit->MaxAdc()*TMath::Exp(1)/c );
      f1.SetParameter( 1, time+c );
      f1.SetParLimits( 1, time-2*c, time+4*c);
      f1.SetParameter( 2, c );
      f1.SetParLimits( 2, 0.5*c, 1.5*c );
      f1.FixParameter( 3, m_pedestal );
      h2.Fit( "f1", option, "", time-2*c, time+4*c );
      auto p = f1.GetParameters();
      time = p[1] + p[2]; // peak time
      // c = p[2];
      de = p[0]*p[2]*TMath::Exp(-1); // amplitude
      // de = p[0]*p[2]*p[2]; // integral
      Double_t chisqr = f1.GetChisquare() / f1.GetNDF();
#endif
      if( chisqr > MaxChisqr ||
	  de < MinDeTPC )
	continue;
      m_time.push_back( time );
      m_de.push_back( de );
      m_chisqr.push_back( chisqr );
#if DebugEvDisp
      hddaq::cout << FUNC_NAME << " (time,de,chisqr)=("
                << time << "," << de << "," << chisqr << ")" << std::endl;
#endif
    }

  }
#endif

#if DebugEvDisp
  h2.Draw("L");
  c1.Modified();
  c1.Update();
  gSystem->ProcessEvents();
  Print();
  getchar();
#endif

  gErrorIgnoreLevel = level;
  m_is_good = true;
  return true;
}


//______________________________________________________________________________
double
TPCHit::GetResolutionX( void )
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
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
  //double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);

  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);
  double rho =  sqrt(pow(z_pos,2)+pow(x_pos,2));
  double dalpha =sT/rho;
  double smear_alpha = alpha + dalpha;
  double res_xdiff = fabs(rho*(sin(smear_alpha)-sin(alpha)));
  double res_x = sqrt(s0*s0+res_xdiff*res_xdiff);

  //hddaq::cout<<"res_x: "<<res_x<<std::endl;

  return res_x;
  //return 0.2;
}

//______________________________________________________________________________
double
TPCHit::GetResolutionZ( void )
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  double y_pos= m_pos.Y();
  double s0 = gUser.GetParameter("TPC_sigma0");
  //s0 is considered for common resolution
//  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double Dt = gUser.GetParameter("TPC_Dt");
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1
  double e_ALD = exp(-1.*A*L_D);
  //  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);


  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);
  double rho =  sqrt(pow(z_pos,2)+pow(x_pos,2));
  double dalpha =sT/rho;
  double smear_alpha = alpha + dalpha;
  double res_zdiff = fabs(rho*(cos(smear_alpha)-cos(alpha)));
  double res_z = sqrt(s0*s0 + res_zdiff*res_zdiff);

  //hddaq::cout<<"res_z: "<<res_z<<std::endl;

  return res_z;
  //return 0.2;
}

//______________________________________________________________________________
double
TPCHit::GetResolutionY( void )
{
  // temporary
  return 0.5;
}

//______________________________________________________________________________
double
TPCHit::GetResolution( void )
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  double y_pos= m_pos.Y();
  double s0 = gUser.GetParameter("TPC_sigma0");
  //s0 is considered for common resolution
//  double Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  double Dt = gUser.GetParameter("TPC_Dt");
  double L_D = 30.+(y_pos*0.1);//cm
  double N_eff = 42.8;
  double A = 0.0582*0.01;//m-1 -> cm-1
  double e_ALD = exp(-1.*A*L_D);
  //  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);

  double x_pos= m_pos.X();
  double z_pos= m_pos.Z() - zTgtTPC;
  double alpha =  atan2(x_pos,z_pos);
  double rho =  sqrt(pow(z_pos,2)+pow(x_pos,2));
  double dalpha =sT/rho;
  double smear_alpha = alpha + dalpha;
  double res_x_diff = fabs(rho*(sin(smear_alpha)-sin(alpha)));
  double res_z_diff = fabs(rho*(cos(smear_alpha)-cos(alpha)));

  double res_x = sqrt(s0*s0+res_x_diff*res_x_diff);
  double res_y = 0.5;
  double res_z = sqrt(s0*s0+res_z_diff*res_z_diff);

  double tot_res = sqrt(res_x*res_x + res_y*res_y + res_z*res_z);
  //double tot_res = sqrt(0.5*0.5 + 0.5*0.5 + 0.5*0.5);
  //double tot_res = sqrt(0.2*0.2 + 0.5*0.5 + 0.2*0.2);

  return tot_res;
}


//______________________________________________________________________________
void
TPCHit::Print( const std::string& arg, std::ostream& ost ) const
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
  for( Int_t i=0, n=GetNHits(); i<n; ++i ){
    ost << std::setw(3) << std::right << i << " "
        << "(time,de,chisqr)=("
        << std::fixed << std::setprecision(1) << std::setw(9) << m_time[i]
        << std::fixed << std::setprecision(1) << std::setw(9) << m_de[i]
        << std::fixed << std::setprecision(3) << std::setw(9) << m_chisqr[i]
        << ")" << std::endl;
  }
}
