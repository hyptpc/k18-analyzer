/**
 *  file: DCHit.cc
 *  date: 2017.04.10
 *
 */

#include "DCHit.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCParameters.hh"
#include "DCTdcCalibMan.hh"
#include "DCLTrackHit.hh"
#include "DebugCounter.hh"
#include "MathTools.hh"
#include "RootHelper.hh"
#include "SsdParamMan.hh"

namespace
{
  const std::string& class_name("DCHit");
  const DCGeomMan&       gGeom  = DCGeomMan::GetInstance();
  const DCTdcCalibMan&   gTdc   = DCTdcCalibMan::GetInstance();
  const DCDriftParamMan& gDrift = DCDriftParamMan::GetInstance();
  // const SsdParamMan&     gSsd   = SsdParamMan::GetInstance();
  const bool SelectTDC1st  = false;
}

//______________________________________________________________________________
DCHit::DCHit( void )
  : m_layer(-1), m_wire(-1),
    m_wpos(-9999.), m_angle(0.),
    m_cluster_size(0.),
    m_mwpc_flag(false),
    m_is_ssd(false),
    m_zero_suppressed(false),
    m_time_corrected(false),
    m_good_waveform(false),
    m_pedestal(-999),
    m_peak_height(-999),
    m_peak_position(-999),
    m_deviation(-999.),
    m_amplitude(-999.),
    m_peak_time(-999.),
    m_adc_sum(-999.),
    m_de(-999.),
    m_rms(-999.),
    m_chisqr(9999.),
    m_belong_kaon(false)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer )
  : m_layer( layer ), m_wire(-1),
    m_wpos(-9999.), m_angle(0.),
    m_cluster_size(0.),
    m_mwpc_flag(false),
    m_is_ssd(false),
    m_zero_suppressed(false),
    m_time_corrected(false),
    m_good_waveform(false),
    m_pedestal(-999),
    m_peak_height(-999),
    m_peak_position(-999),
    m_deviation(-999.),
    m_amplitude(-999.),
    m_peak_time(-999.),
    m_adc_sum(-999.),
    m_de(-999.),
    m_rms(-999.),
    m_chisqr(-999.),
    m_belong_kaon(false)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer, double wire )
  : m_layer(layer), m_wire(wire),
    m_wpos(-9999.), m_angle(0.),
    m_cluster_size(0.),
    m_mwpc_flag(false),
    m_is_ssd(false),
    m_zero_suppressed(false),
    m_time_corrected(false),
    m_good_waveform(false),
    m_pedestal(-999),
    m_peak_height(-999),
    m_peak_position(-999),
    m_deviation(-999.),
    m_amplitude(-999.),
    m_peak_time(-999.),
    m_adc_sum(-999.),
    m_de(-999.),
    m_rms(-999.),
    m_chisqr(-999.),
    m_belong_kaon(false)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::~DCHit( void )
{
  ClearRegisteredHits();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
DCHit::SetDummyPair()
{
  data_pair a_pair = {0.,
		      0.,
		      std::numeric_limits<double>::quiet_NaN(),
		      std::numeric_limits<double>::quiet_NaN(),
		      -1,
		      false,
		      true};
  m_pair_cont.push_back(a_pair);
}

//______________________________________________________________________________
void
DCHit::ClearRegisteredHits( void )
{
  int n = m_register_container.size();
  for(int i=0; i<n; ++i){
    delete m_register_container[i];
  }
}

//______________________________________________________________________________
bool
DCHit::CalcDCObservables( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( !gGeom.IsReady() ||
      !gTdc.IsReady()  ||
      !gDrift.IsReady() ){
    return false;
  }

  m_wpos  = gGeom.CalcWirePosition( m_layer,m_wire );
  m_angle = gGeom.GetTiltAngle( m_layer );
  m_z     = gGeom.GetLocalZ( m_layer );

  bool status = true;
  int  nh_tdc      = m_tdc.size();
  int  nh_trailing = m_trailing.size();

  IntVec leading_cont, trailing_cont;

  // Prepare
  {
    for ( int m = 0; m < nh_tdc; ++m ) {
      leading_cont.push_back( m_tdc.at( m ) );
    }
    for ( int m = 0; m < nh_trailing; ++m ) {
      trailing_cont.push_back( m_trailing.at( m ) );
    }

    std::sort(leading_cont.begin(),  leading_cont.end(),  std::greater<int>());
    std::sort(trailing_cont.begin(), trailing_cont.end(), std::greater<int>());

    int i_t = 0;
    for(int i = 0; i<nh_tdc; ++i){
      data_pair a_pair = {0., 0.,
			  std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN(),
			  -1, false, false};

      int leading  = leading_cont.at(i);
      while(i_t < nh_trailing){
	int trailing = trailing_cont.at(i_t);

	if(leading > trailing){
	  a_pair.index_t = i_t;
	  m_pair_cont.push_back(a_pair);
	  break;
	}else{
	  ++i_t;
	}// Goto next trailing
      }

      if(i_t == nh_trailing){
	a_pair.index_t = -1;
	m_pair_cont.push_back(a_pair);
	continue;
      }// no more trailing data
    }// for(i)
  }

  // Delete duplication index_t
  for(int i = 0; i<nh_tdc-1; ++i){
    if(true
       && m_pair_cont.at(i).index_t != -1
       && m_pair_cont.at(i).index_t == m_pair_cont.at(i+1).index_t)
      {
	m_pair_cont.at(i).index_t = -1;
      }
  }// for(i)


  for ( int i=0; i<nh_tdc; ++i ) {
    double ctime;
    if( !gTdc.GetTime( m_layer, m_wire, leading_cont.at(i), ctime ) ){
      return false;
    }

    double dtime, dlength;
    if( !gDrift.CalcDrift( m_layer, m_wire, ctime, dtime, dlength ) ){
      status = false;
    }

    m_pair_cont.at(i).drift_time   = dtime;
    m_pair_cont.at(i).drift_length = dlength;

    if(m_pair_cont.at(i).index_t != -1){
      double trailing_ctime;
      gTdc.GetTime( m_layer, m_wire, trailing_cont.at(m_pair_cont.at(i).index_t), trailing_ctime );
      m_pair_cont.at(i).trailing_time = trailing_ctime;
      m_pair_cont.at(i).tot           = ctime - trailing_ctime;
    }else{
      m_pair_cont.at(i).trailing_time = std::numeric_limits<double>::quiet_NaN();
      m_pair_cont.at(i).tot           = std::numeric_limits<double>::quiet_NaN();
    }


    switch( m_layer ){
      // BC3,4
    case 113: case 114: case 115: case 116: case 117: case 118:
    case 119: case 120: case 121: case 122: case 123: case 124:
      if( MinDLBc[m_layer-100] < m_pair_cont.at(i).drift_length && m_pair_cont.at(i).drift_length < MaxDLBc[m_layer-100] ){
	m_pair_cont.at(i).dl_range = true;
      }
      break;

      // SDC1,2,3
    case 1: case 2: case 3: case 4: case 5: case 6:
    case 31: case 32: case 33: case 34:
    case 35: case 36: case 37: case 38:
      if( MinDLSdc[m_layer] < m_pair_cont.at(i).drift_length && m_pair_cont.at(i).drift_length < MaxDLSdc[m_layer] ){
	m_pair_cont.at(i).dl_range = true;
      }
      break;
    default:
      hddaq::cout << "#E " << func_name << " "
		  << "invalid layer id : " << m_layer << std::endl;
      status = false;
      break;
    }
  }

  if( SelectTDC1st && status ){
    int       tdc1st   = 0;
    data_pair pair1st;
    for( int i=0; i<nh_tdc; ++i ){
      if( tdc1st < leading_cont[i]
	  /* && m_dl_range[i] */ ){
	tdc1st   = leading_cont[i];
	pair1st  = m_pair_cont.at(i);
      }
    }
    if( tdc1st>0 ){
      m_tdc.clear();
      m_pair_cont.clear();
      m_tdc.push_back( tdc1st );
      m_pair_cont.push_back( pair1st );
    } else {
      m_tdc.clear();
      m_pair_cont.clear();
    }
  }

  return status;
}

//______________________________________________________________________________
bool
DCHit::CalcMWPCObservables( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( !gGeom.IsReady() || !gTdc.IsReady() )
    return false;

  m_angle = gGeom.GetTiltAngle( m_layer );
  m_wpos  = gGeom.CalcWirePosition( m_layer, m_wire );
  m_z     = gGeom.GetLocalZ( m_layer );

  bool status = true;
  int  nh_tdc      = m_tdc.size();
  int  nh_trailing = m_trailing.size();

  IntVec leading_cont, trailing_cont;

  // Prepare
  {
    for ( int m = 0; m < nh_tdc; ++m ) {
      leading_cont.push_back( m_tdc.at( m ) );
    }
    for ( int m = 0; m < nh_trailing; ++m ) {
      trailing_cont.push_back( m_trailing.at( m ) );
    }

    std::sort(leading_cont.begin(),  leading_cont.end(),  std::greater<int>());
    std::sort(trailing_cont.begin(), trailing_cont.end(), std::greater<int>());

    int i_t = 0;
    for(int i = 0; i<nh_tdc; ++i){
      data_pair a_pair = {0., 0.,
			  std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN(),
			  -1, false, false};

      int leading  = leading_cont.at(i);
      while(i_t < nh_trailing){
	int trailing = trailing_cont.at(i_t);

	if(leading > trailing){
	  a_pair.index_t = i_t;
	  m_pair_cont.push_back(a_pair);
	  break;
	}else{
	  ++i_t;
	}// Goto next trailing
      }

      if(i_t == nh_trailing){
	a_pair.index_t = -1;
	m_pair_cont.push_back(a_pair);
	continue;
      }// no more trailing data
    }// for(i)
  }

  // Delete duplication index_t
  for(int i = 0; i<nh_tdc-1; ++i){
    if(true
       && m_pair_cont.at(i).index_t != -1
       && m_pair_cont.at(i).index_t == m_pair_cont.at(i+1).index_t)
      {
	m_pair_cont.at(i).index_t = -1;
      }
  }// for(i)


  for ( int i=0; i<nh_tdc; ++i ) {
    double ctime;
    if( !gTdc.GetTime( m_layer, m_wire, leading_cont.at(i), ctime ) ){
      return false;
    }

    double dtime, dlength;
    if( !gDrift.CalcDrift( m_layer, m_wire, ctime, dtime, dlength ) ){
      status = false;
    }

    m_pair_cont.at(i).drift_time   = dtime;
    m_pair_cont.at(i).drift_length = dlength;

    if(m_pair_cont.at(i).index_t != -1){
      double trailing_ctime;
      gTdc.GetTime( m_layer, m_wire, trailing_cont.at(i), trailing_ctime );
      m_pair_cont.at(i).trailing_time = trailing_ctime;
      m_pair_cont.at(i).tot           = ctime - trailing_ctime;
    }else{
      m_pair_cont.at(i).trailing_time = std::numeric_limits<double>::quiet_NaN();
      m_pair_cont.at(i).tot           = std::numeric_limits<double>::quiet_NaN();
    }

    if( m_pair_cont.at(i).drift_time > MinDLBc[m_layer-100] && m_pair_cont.at(i).drift_time < MaxDLBc[m_layer-100] ){
      m_pair_cont.at(i).dl_range = true;
    }else{
      status = false;
    }
  }

  return status;
}

//______________________________________________________________________________
bool
DCHit::CalcFiberObservables( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( !gGeom.IsReady() ) return false;

  m_angle = gGeom.GetTiltAngle( m_layer );
  m_z     = gGeom.GetLocalZ( m_layer );

  bool status = true;
  std::size_t nh_tdc = m_tdc.size();

  for( std::size_t i=0; i<nh_tdc; i++ ){
    data_pair a_pair = {(double)m_tdc[i], 0.,
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN(),
			-1, false, true};
    m_pair_cont.push_back( a_pair );
  }

  return status;
}

//______________________________________________________________________________
bool
DCHit::CalcSsdObservables( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  // if( !gGeom.IsReady() || !gSsd.IsReady() ) return false;

  m_angle = gGeom.GetTiltAngle( m_layer );
  m_wpos  = gGeom.CalcWirePosition( m_layer, m_wire );
  m_z     = gGeom.GetLocalZ( m_layer );

  std::size_t nhadc = m_adc.size();
  if( nhadc!=NumOfSampleSSD ){
    hddaq::cerr << "#W " << func_name << " the number of sample is wrong " << std::endl
		<< "   layer#" << m_layer << " segment#" << m_wire
		<< " [" << nhadc << "/" << NumOfSampleSSD << "]" << std::endl;
    return false;
  }

  if( m_trailing.size()>0 ) m_zero_suppressed = true;

  double pedestal = m_adc[0];
  std::vector<double> dE(nhadc);
  std::vector<double> rms(nhadc);
  for( std::size_t i=0; i<nhadc; ++i ){
    data_pair a_pair = {(double)m_tdc[i], 0.,
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN(),
			-1, false, true};
    m_pair_cont.push_back( a_pair );

    // if( m_adc[i]<pedestal )
    //   pedestal = m_adc[i];

    // if( !gSsd.GetDe( m_layer, m_wire, i, m_adc[i], dE[i] ) ){
    //   hddaq::cerr << func_name << " : something is wrong at GetDe("
    // 		  << m_layer   << ", " << m_wire << ", " << i << ", "
    // 		  << m_adc[i]  << ", " << dE[i] << ")"
    // 		  << std::endl;
    //   return false;
    // }

    // if( !gSsd.GetRms( m_layer,m_wire, i, rms[i] ) ){
    //   hddaq::cerr << func_name << " : something is wrong at GetRms("
    // 		  << m_layer   << ", " << m_wire << ", " << i << ", "
    // 		  << rms[i]   << ")"
    // 		  <<std::endl;
    //   return false;
    // }

    dE[i] = m_adc[i] - pedestal;
    m_waveform.push_back( dE[i] );
    m_time.push_back( m_tdc[i]*SamplingIntervalSSD );

    if( m_adc[i]>m_peak_height ){
      m_peak_time     = m_time[i];
      m_rms           = rms[i];
      m_peak_height   = m_adc[i];
      m_peak_position = m_tdc[i];
    }
  }

  m_pedestal  = pedestal;
  m_adc_sum   = math::Accumulate(dE);
  m_deviation = math::Deviation(dE);

  /*** SSD Waveform Fitting ***************************
   *
   *  f(x) = a * (x-b) * exp(-(x-b)/c)
   *    a : scale factor
   *    b : start timing
   *    c : decay constant
   *
   *  f'(x)          = a/c * (-x+b+c) * exp(-(x-b)/c)
   *  f'(x)|x=b+c    = 0.
   *  f(b+c)         = a * c * exp(-1)
   *  Sf(x)dx|b->inf = a * c^2
   *
   *  peak time = b + c
   *  amplitude = a * c * exp(-1)
   *  integral  = a * c^2
   *
   ****************************************************/

  // TGraphErrors graph( NumOfSampleSSD, &m_time[0], dE, 0, rms );
  TGraph graph( NumOfSampleSSD, &(m_time[0]), &(dE[0]) );

  double xmin =  40.;
  double xmax = 210.;
  TF1 func( "func", "[0]*(x-[1])*exp(-(x-[1])/[2])", xmin, xmax );
  func.SetParameter( 0, dE[3]*std::exp(1)/60. );
  func.SetParLimits( 0, 0., 50000.*std::exp(-1.) );
  func.SetParameter( 1, 40. );
  func.SetParLimits( 1, 10., 100. );
  func.FixParameter( 2, 50. );

  graph.Fit("func", "RQ");
  double p[3];
  func.GetParameters(p);

  m_peak_time = p[1] + p[2];
  m_amplitude = p[0] * p[2] * std::exp(-1.);
  // m_de        = func.Integral( p[1], math::Infinity() );
  m_de        = p[0]*p[2]*p[2];
  m_chisqr    = func.GetChisquare() / func.GetNDF();

  m_is_ssd      = true;
  m_good_waveform = true;

  if( m_de<0.1 ) m_good_waveform = false;

#if 0
  static bool flag = true;
  if( flag && m_chisqr<5. && m_de>1000. ){
    TCanvas c1;
    graph.SetTitle("");
    graph.GetXaxis()->SetTitle("[ns]");
    graph.GetYaxis()->SetTitle("[ch]");
    graph.SetMarkerStyle(8);
    graph.SetLineWidth(2);
    graph.SetMarkerSize(2);
    graph.Draw("APL");
    c1.Print("tmp.pdf");
    flag = false;
  }
#endif

  return true;
}

// bool DCHit::CalcObservablesSimulation( double dlength)
// {
//   static const std::string func_name="[DCHit::CalcObservablesSimulation]";

//   if( !gGeom.IsReady() ) return false;

//   m_wpos=gGeom.CalcWirePosition(m_layer,m_wire);
//   m_angle=gGeom.GetTiltAngle(m_layer);

//   m_dl = dlength;
//   bool status=true;

//   if(m_layer>=100){
//     if( m_dl>MinDLBc[m_layer-100] && m_dl<MaxDLBc[m_layer-100] )
//       m_dl_range=true;
//   }
//   else {
//     if( m_dl>MinDLSdc[m_layer] && m_dl<MaxDLSdc[m_layer] )
//       m_dl_range=true;
//   }
//   return status;
// }

//______________________________________________________________________________
bool
DCHit::DoTimeCorrection( double offset )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_time_corrected ){
    Print(func_name+" already corrected!");
    return false;
  }

#if 0
  Print("Before Correction");
#endif

  DoubleVec ctime;
  int nh = m_time.size();
  for(int i=0; i<nh; ++i){
    ctime.push_back( m_time[i] + offset );
  }

  m_time = ctime;
  m_peak_time    += offset;
  m_time_corrected = true;

#if 0
  Print("After Correction");
#endif

  return true;
}

//______________________________________________________________________________
double
DCHit::GetResolution( void ) const
{
  return gGeom.GetResolution(m_layer);
}

//______________________________________________________________________________
void
DCHit::TotCut(double min_tot, bool adopt_nan)
{
  auto itr_new_end = std::remove_if(m_pair_cont.begin(), m_pair_cont.end(),
				    [min_tot, adopt_nan](data_pair a_pair)->bool
				    {return (isnan(a_pair.tot) && adopt_nan) ? false : !(a_pair.tot > min_tot);}
				    );
  m_pair_cont.erase(itr_new_end, m_pair_cont.end());
}

//______________________________________________________________________________
void
DCHit::Print( const std::string& arg, std::ostream& ost ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const int w = 16;
  ost << "#D " << func_name << " " << arg << std::endl
      << std::setw(w) << std::left << "layer" << m_layer << std::endl
      << std::setw(w) << std::left << "wire"  << m_wire  << std::endl
      << std::setw(w) << std::left << "wpos"  << m_wpos  << std::endl
      << std::setw(w) << std::left << "angle" << m_angle << std::endl
      << std::setw(w) << std::left << "z"     << m_z     << std::endl
      << std::setw(w) << std::left << "kaon"  << m_belong_kaon << std::endl;

  if(m_is_ssd){
    ost << std::setw(w) << std::left << "zero_suppressed"
	<< m_zero_suppressed << std::endl;
    ost << std::setw(w) << std::left << "time_corrected"
	<< m_time_corrected << std::endl;
    ost << std::setw(w) << std::left << "good_waveform"
	<< m_good_waveform << std::endl;
    ost << std::setw(w) << std::left << "adc" << m_adc.size() << " : ";
    std::copy(m_adc.begin(), m_adc.end(),
	      std::ostream_iterator<int>(ost, " "));
    ost << std::endl;
    ost << std::setw(w) << std::left << "peak_height"
	<< m_peak_height << std::endl;
    ost << std::setw(w) << std::left << "peak_position"
	<< m_peak_position << std::endl;
  }

  ost << std::setw(w) << std::left << "tdc" << m_tdc.size() << " : ";
  std::copy( m_tdc.begin(), m_tdc.end(),
	     std::ostream_iterator<int>(ost, " ") );
  ost << std::endl;

  if(!m_is_ssd){
    std::size_t n_pair = m_pair_cont.size();
    DoubleVec dt(n_pair);
    DoubleVec dl(n_pair);
    DoubleVec ttime(n_pair);
    DoubleVec tot(n_pair);
    BoolVec   belong(n_pair);
    BoolVec   range(n_pair);

    for(int i = 0; i<m_pair_cont.size(); ++i){
      dt[i]     = m_pair_cont.at(i).drift_time;
      dl[i]     = m_pair_cont.at(i).drift_length;
      ttime[i]  = m_pair_cont.at(i).trailing_time;
      tot[i]    = m_pair_cont.at(i).tot;
      belong[i] = m_pair_cont.at(i).belong_track;
      range[i]  = m_pair_cont.at(i).dl_range;
    }

    ost << std::endl << std::setw(w) << std::left
	<< "trailing" << m_trailing.size() << " : ";
    std::copy(m_trailing.begin(), m_trailing.end(),
	      std::ostream_iterator<int>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "drift time" << n_pair << " : ";
    std::copy(dt.begin(), dt.end(),
	      std::ostream_iterator<double>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "drift length" << n_pair << " : ";
    std::copy(dl.begin(), dl.end(),
	      std::ostream_iterator<double>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "trailing time" << n_pair << " : ";
    std::copy(ttime.begin(), ttime.end(),
	      std::ostream_iterator<double>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "tot" << n_pair << " : ";
    std::copy(tot.begin(), tot.end(),
	      std::ostream_iterator<double>(ost, " "));
    ost << std::endl;
    ost << std::setw(w) << std::left
	<< "belongTrack" << n_pair << " : ";
    std::copy(belong.begin(), belong.end(),
	      std::ostream_iterator<bool>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "dlRange" << n_pair << " : ";
    std::copy(range.begin(), range.end(),
	      std::ostream_iterator<bool>(ost, " "));
    ost << std::endl;
  }

  if(m_is_ssd){
    ost << std::setw(w) << std::left << "waveform" << m_waveform.size() << " : ";
    std::copy(m_waveform.begin(), m_waveform.end(),
	      std::ostream_iterator<int>(ost, " "));
    ost << std::endl;
    ost << std::setw(w) << std::left
	<< "time" << m_time.size() << " : ";
    std::copy(m_time.begin(), m_time.end(),
	      std::ostream_iterator<int>(ost, " "));
    ost << std::endl;
    ost << std::setw(w) << std::left << "chisqr"    << m_chisqr    << std::endl;
    ost << std::setw(w) << std::left << "deviation" << m_deviation << std::endl;
    ost << std::setw(w) << std::left << "amplitude" << m_amplitude << std::endl;
    ost << std::setw(w) << std::left << "delta E"   << m_de        << std::endl;
  }

  if(m_mwpc_flag){
    ost << std::endl
	<< std::setw(w) << std::left << "clsize" << m_cluster_size << std::endl
	<< std::setw(w) << std::left << "mean wire" << m_mwpc_wire << std::endl
	<< std::setw(w) << std::left << "mean pos"  << m_mwpc_wpos << std::endl;
  }
}
