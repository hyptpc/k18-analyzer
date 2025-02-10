/**
 *  file: ScalerAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "ScalerAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
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
  enum eDisp   { kLeft, kCenter, kRight, MaxColumn };
  static const std::size_t MaxRow = 17;
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
      m_info.at(i).at(j).data = 0;
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
  // {
  //   static const int trig_id = gUnpacker.get_device_id("TriggerFlag");
  //   static const int tdc_id = gUnpacker.get_data_id("TriggerFlag", "leading");
  //   int nhit = gUnpacker.get_entries( trig_id, 0, , 0, tdc_id );
  //   if( nhit>0 ){
  //     int tdc = gUnpacker.get( DetIdTrig, 0, SpillEndFlag, 0, tdc_id );
  //     if( tdc>0 ) m_spill_increment = true;
  //   }
  // }

  //////////////////// Check spill increment
  // if(0){
  //   int i=0;
  //   int j=0;
  //   int module_id = m_info.at(i).at(j).module_id;
  //   int channel   = m_info.at(i).at(j).channel;
  //   int nhit = gUnpacker.get_entries( DetIdScaler, module_id, 0, channel, 0 );
  //   if( nhit>0 ){
  //     Scaler val = gUnpacker.get( DetIdScaler, module_id, 0, channel, 0 );
  //     if( m_info.at(i).at(j).prev > val ){
  // 	//
  //     }
  //   }
  // }
  //////////////////// Scaler Data
  for( std::size_t i=0; i<MaxColumn; ++i ){
    for( std::size_t j=0; j<MaxRow; ++j ){
      if( !m_info.at(i).at(j).flag_disp ) continue;
      int module_id = m_info.at(i).at(j).module_id;
      int channel   = m_info.at(i).at(j).channel;
      int nhit = gUnpacker.get_entries( DetIdScaler, module_id, 0, channel, 0 );
      if( nhit<=0 ) continue;
      Scaler val = gUnpacker.get( DetIdScaler, module_id, 0, channel, 0 );
      if( m_info.at(i).at(j).prev > val ){
	m_spill_increment = true;
	m_info.at(i).at(j).prev = 0;
      }
      m_info.at(i).at(j).curr  = val;
      if(!m_init)  m_info.at(i).at(j).data += val - m_info.at(i).at(j).prev;
      m_info.at(i).at(j).prev  = m_info.at(i).at(j).curr;
    }
  }
  if( m_init ) m_init=false;   
  if( m_spill_increment ){
    m_info[0][0].data++; // spill
    //    std::cout<<"Spill incremented "<<m_info[0][0].data<<"  "<<m_info[2][16].data<<"  "<<m_info[2][16].curr<<std::endl;
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
      m_info.at(i).at(j) = ScalerInfo("n/a", 0, 0, false );
    }
  }
  m_init=true;
  return true;
}
//______________________________________________________________________________
bool
ScalerAnalyzer::Setup( const std::string& filename )
{
  std::cout<<filename<<std::endl;
  std::ifstream ifs(filename);
  if(!ifs){
    std::cout<<"#E "<<filename<<" cannot be opend"<<std::endl;
    exit(0);
  }
  std::cout<<filename<<" opened"<<std::endl; 
  std::string name;
  int plane,channel;
  int column,row;
  while(ifs>>name>>plane>>channel>>column>>row){
#if 1
    std::cout<<std::setw(20)<<name
             <<std::setw(5)<<plane
             <<std::setw(5)<<channel
             <<std::setw(5)<<column
             <<std::setw(5)<<row
             <<std::endl;
#endif
    if(column<0||row<0) continue;
    if(column>=MaxColumn||row>=MaxRow) continue;
    Set( column, row, ScalerInfo(name,plane,channel) );
  }
  std::cout<<""<<__func__<<" finished"<<std::endl;
  ifs.close();
  return true;
}
bool
ScalerAnalyzer::Setup_RunNum( int run_number )
{
  return true;
}

//______________________________________________________________________________
void
ScalerAnalyzer::Print( const std::string& arg ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  std::cout << "\033[2J"
	<< func_name << " " << arg << std::endl;

  double l1_req = (double)Get("Request");
  double l1_acc = (double)Get("Accept");
  double daq_eff= l1_acc/l1_req;

  for( std::size_t i=0; i<MaxRow; ++i ){
    if( m_separate_comma ){
      std::cout << std::left  << std::setw(16) << m_info[kLeft][i].name
	    << std::right << std::setw(16) << SeparateComma( m_info[kLeft][i].data )
	    << " : "
	    << std::left  << std::setw(16) << m_info[kRight][i].name
	    << std::right << std::setw(16) << SeparateComma( m_info[kRight][i].data )
	    <<std::endl;
    }else{
      std::cout << std::left  << std::setw(16) << m_info[kLeft][i].name
	    << std::right << std::setw(16) << m_info[kLeft][i].data
	    << " : "
	    << std::left  << std::setw(16) << m_info[kCenter][i].name
	    << std::right << std::setw(16) << m_info[kCenter][i].data
	    << " : "
	    << std::left  << std::setw(16) << m_info[kRight][i].name
	    << std::right << std::setw(16) << m_info[kRight][i].data
	    <<std::endl;
    }
  }
  std::cout << std::endl  << std::setprecision(6) << std::fixed
	<< std::left  << std::setw(16) << "DAQ Eff"
	<< std::right  << std::setw(16) << daq_eff
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
//______________________________________________________________________________
void
ScalerAnalyzer::Set( Int_t i, Int_t j, const ScalerInfo& info )
{
  if( i >= MaxColumn || j >= MaxRow ){
    return;
  }
  else {
    m_info.at(i).at(j) = info;
  }
}
//______________________________________________________________________________
void
ScalerAnalyzer::WriteToFile(std::ofstream &ofs)
{
  for( std::size_t i=0; i<MaxColumn; ++i ){
    for( std::size_t j=0; j<MaxRow; ++j ){
      if( !FlagDisp(i,j) ) continue;
      // ofs<<std::setw(16)<<gScaler.GetName(i,j)
      // 	 <<std::setw(16)<<gScaler.Get(i,j)
      ofs<<GetName(i,j)<<'\t'
	 <<Get(i,j)
	 <<std::endl;
    }
  }
}
