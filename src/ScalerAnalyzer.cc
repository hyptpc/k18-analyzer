/**
 *  file: ScalerAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "ScalerAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "DCRawHit.hh"
#include "HodoRawHit.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

namespace
{
  using namespace hddaq::unpacker;
  const std::string& class_name("ScalerAnalyzer");
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
}

//______________________________________________________________________________
ScalerAnalyzer::ScalerAnalyzer( void )
  : m_ost( std::cout ), //m_ost( hddaq::cout ),
    m_info( MaxColumn, std::vector<ScalerInfo>(MaxRow) ),
    m_spill_increment(false),
    m_separate_comma(false)
{
  Initialize();
}

//______________________________________________________________________________
ScalerAnalyzer::~ScalerAnalyzer( void )
{
}

//______________________________________________________________________________
void
ScalerAnalyzer::Clear( void )
{
  for( std::size_t i=0; i<MaxColumn; ++i ){
    for( std::size_t j=0; j<MaxRow; ++j ){
      m_info[i][j].data = 0;
    }
  }
}

//______________________________________________________________________________
bool
ScalerAnalyzer::Decode( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_spill_increment = false;

  //////////////////// Trigger Flag
  {
    static const int tdc_id = gUnpacker.get_data_id("TFlag", "tdc");
    int nhit = gUnpacker.get_entries( DetIdTrig, 0, SpillEndFlag, 0, tdc_id );
    if( nhit>0 ){
      int tdc = gUnpacker.get( DetIdTrig, 0, SpillEndFlag, 0, tdc_id );
      if( tdc>0 ) m_spill_increment = true;
    }
  }

  //////////////////// Scaler Data
  for( std::size_t i=0; i<MaxColumn; ++i ){
    for( std::size_t j=0; j<MaxRow; ++j ){
      if( !m_info[i][j].flag_disp ) continue;
      int module_id = m_info[i][j].module_id;
      int channel   = m_info[i][j].channel;
      int nhit = gUnpacker.get_entries( DetIdScaler, module_id, 0, channel, 0 );
      if( nhit<=0 ) continue;
      Scaler val = gUnpacker.get( DetIdScaler, module_id, 0, channel, 0 );

      if( m_info[i][j].prev > val ){
	m_info[i][j].prev = 0;
      }

      m_info[i][j].curr  = val;
      m_info[i][j].data += val - m_info[i][j].prev;
      m_info[i][j].prev  = m_info[i][j].curr;
    }
  }

  if( m_spill_increment ){
    m_info[1][0].data++; // spill
  }

  return true;
}

//______________________________________________________________________________
Scaler
ScalerAnalyzer::Get( const std::string& name ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  std::vector<ScalerInfo>::const_iterator itr, itr_end;
  for( std::size_t i=0; i<MaxColumn; ++i ){
    itr_end = m_info[i].end();
    for( itr=m_info[i].begin(); itr!=itr_end; ++itr ){
      if( itr->name == name ) return itr->data;
    }
  }

  hddaq::cerr << "#W " << func_name << " "
	      << "no such ScalerInfo : " << name << std::endl;

  return 0;
}

//______________________________________________________________________________
bool
ScalerAnalyzer::Initialize( void )
{
  for( std::size_t i=0; i<MaxColumn; ++i ){
    for( std::size_t j=0; j<MaxRow; ++j ){
      m_info[i][j] = ScalerInfo("n/a", i, j, false );
    }
  }

  // scaler information is defined from here.
  // using a space character is not recommended.

  // left column (counter info)
  int index = 0;
  m_info[kLeft][index++] = ScalerInfo( "K_beam",       0, 19 );
  m_info[kLeft][index++] = ScalerInfo( "Pi_beam",      0, 20 );
  m_info[kLeft][index++] = ScalerInfo( "/p_beam",      0, 21 );
  // m_info[kLeft][index++] = ScalerInfo( "/p_beam(0.6)", 0, 26 );
  m_info[kLeft][index++] = ScalerInfo( "BH1",          0,  0 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-1",        0,  1 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-2",        0,  2 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-3",        0,  3 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-4",        0,  4 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-5",        0,  5 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-6",        0,  6 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-7",        0,  7 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-8",        0,  8 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-9",        0,  9 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-10",       0, 10 );
  m_info[kLeft][index++] = ScalerInfo( "BH1-11",       0, 11 );
  m_info[kLeft][index++] = ScalerInfo( "BH2",          0, 12 );
  m_info[kLeft][index++] = ScalerInfo( "BAC1",         0, 13 );
  m_info[kLeft][index++] = ScalerInfo( "BAC2",         0, 14 );
  m_info[kLeft][index++] = ScalerInfo( "FBH",          0, 10 );
  m_info[kLeft][index++] = ScalerInfo( "PVAC",         0, 15 );
  m_info[kLeft][index++] = ScalerInfo( "FAC",          0, 16 );
  m_info[kLeft][index++] = ScalerInfo( "SCH",          0, 17 );
  m_info[kLeft][index++] = ScalerInfo( "TOF",          0, 18 );

  // right column (DAQ info)
  index = 0;
  m_info[kRight][index++] = ScalerInfo( "Spill",        0, 32 );
  m_info[kRight][index++] = ScalerInfo( "10M_Clock",    0, 33 );
  m_info[kRight][index++] = ScalerInfo( "IM",           0, 34 );
  m_info[kRight][index++] = ScalerInfo( "TM",           0, 35 );
  m_info[kRight][index++] = ScalerInfo( "Real_Time",    0, 36 );
  m_info[kRight][index++] = ScalerInfo( "Live_Time",    0, 37 );
  m_info[kRight][index++] = ScalerInfo( "L1_Req",       0, 38 );
  m_info[kRight][index++] = ScalerInfo( "L1_Acc",       0, 39 );
  m_info[kRight][index++] = ScalerInfo( "MTX_Acc",      0, 40 );
  m_info[kRight][index++] = ScalerInfo( "MTX-1",        0, 41 );
  m_info[kRight][index++] = ScalerInfo( "MsT_Acc",      0, 44 );
  m_info[kRight][index++] = ScalerInfo( "MsT_Clear",    0, 45 );
  m_info[kRight][index++] = ScalerInfo( "MsT_ClearPS",  0, 46 );
  m_info[kRight][index++] = ScalerInfo( "L2_Clear",     0, 47 );
  m_info[kRight][index++] = ScalerInfo( "L2_Req",       0, 48 );
  m_info[kRight][index++] = ScalerInfo( "L2_Acc",       0, 49 );
  m_info[kRight][index++] = ScalerInfo( "(ub)",         0, 50 );
  m_info[kRight][index++] = ScalerInfo( "(ub,ub)",      0, 51 );
  m_info[kRight][index++] = ScalerInfo( "(Pi,TOF)",     0, 52 );
  m_info[kRight][index++] = ScalerInfo( "(K,K)",        0, 53 );
  m_info[kRight][index++] = ScalerInfo( "(ub)PS",       0, 54 );
  m_info[kRight][index++] = ScalerInfo( "(ub,ub)PS",    0, 55 );
  m_info[kRight][index++] = ScalerInfo( "(Pi,TOF)PS",   0, 56 );
  m_info[kRight][index++] = ScalerInfo( "(K,K)PS",      0, 57 );
  m_info[kRight][index++] = ScalerInfo( "K_in",         0, 25 );
  m_info[kRight][index++] = ScalerInfo( "Pi_in",        0, 26 );
  m_info[kRight][index++] = ScalerInfo( "K_out",        0, 27 );
  m_info[kRight][index++] = ScalerInfo( "Pi_out",       0, 28 );
  m_info[kRight][index++] = ScalerInfo( "BH1xBH2",      0, 22 );
  m_info[kRight][index++] = ScalerInfo( "BH2xTOF",      0, 23 );
  m_info[kRight][index++] = ScalerInfo( "PVACx/FAC",    0, 24 );

  return true;
}

//______________________________________________________________________________
void
ScalerAnalyzer::Print( const std::string& arg ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  m_ost << "\033[2J"
	<< func_name << " " << arg << std::endl;

  double real_time = (double)Get("Real_Time");
  double live_time = (double)Get("Live_Time");
  double l1_req    = (double)Get("L1_Req");
  double l1_acc    = (double)Get("L1_Acc");
  double l2_acc    = (double)Get("L2_Acc");
  double real_live = live_time/real_time;
  double daq_eff   = l1_acc/l1_req;
  double l2_eff    = l2_acc/l1_acc;
  // double kbeam     = (double)Get("K_beam");
  // double pibeam    = (double)Get("Pi_beam");
  double duty_factor = daq_eff/(1-daq_eff)*(1/real_live-1);
  if( duty_factor >= 1. ) duty_factor = 1.;

  for( std::size_t i=0; i<MaxRow; ++i ){
    if( m_separate_comma ){
      m_ost << std::left  << std::setw(16) << m_info[kLeft][i].name
	    << std::right << std::setw(16) << SeparateComma( m_info[kLeft][i].data )
	    << " : "
	    << std::left  << std::setw(16) << m_info[kRight][i].name
	    << std::right << std::setw(16) << SeparateComma( m_info[kRight][i].data )
	    <<std::endl;
    }else{
      m_ost << std::left  << std::setw(16) << m_info[kLeft][i].name
	    << std::right << std::setw(16) << m_info[kLeft][i].data
	    << " : "
	    << std::left  << std::setw(16) << m_info[kRight][i].name
	    << std::right << std::setw(16) << m_info[kRight][i].data
	    <<std::endl;
    }
  }
  m_ost << std::endl  << std::setprecision(6) << std::fixed
	<< std::left  << std::setw(16) << "Live/Real"
	<< std::right << std::setw(16) << real_live << " : "
	<< std::left  << std::setw(16) << "DAQ Eff"
	<< std::right << std::setw(16) << daq_eff << std::endl
	<< std::left  << std::setw(16) << "L2 Eff"
	<< std::right << std::setw(16) << l2_eff << " : "
	<< std::left  << std::setw(16) << "Duty Factor"
	<< std::right << std::setw(16) << duty_factor << std::endl
	<< std::endl;
}

//______________________________________________________________________________
std::string
ScalerAnalyzer::SeparateComma( Scaler number ) const
{
  std::vector<Scaler> sep_num;

  while(number/1000){
    sep_num.push_back(number%1000);
    number /= 1000;
  }

  std::stringstream ss;  ss << number;
  std::vector<Scaler>::reverse_iterator
    itr, itr_end = sep_num.rend();
  for( itr=sep_num.rbegin(); itr!=itr_end; ++itr ){
    ss << "," << std::setfill('0') << std::setw(3) << *itr;
  }
  return ss.str();
}
