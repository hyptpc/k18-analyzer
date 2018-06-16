/**
 *  file: UserScaler.cc
 *  date: 2017.04.10
 *
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <filesystem_util.hh>
#include <lexical_cast.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "ScalerAnalyzer.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#define USE_COMMA   0
#define SPILL_RESET 0
#define MAKE_LOG    1

namespace
{
  using namespace root;
  const std::string& class_name("EventScaler");
  RMAnalyzer&     gRM     = RMAnalyzer::GetInstance();
  ScalerAnalyzer& gScaler = ScalerAnalyzer::GetInstance();
}

//______________________________________________________________________________
VEvent::VEvent( void )
{
}

//______________________________________________________________________________
VEvent::~VEvent( void )
{
}

//______________________________________________________________________________
class EventScaler : public VEvent
{
public:
        EventScaler( void );
       ~EventScaler( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );

};

//______________________________________________________________________________
struct Event
{
  int evnum;
};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//______________________________________________________________________________
EventScaler::EventScaler( void )
  : VEvent()
{
}

//______________________________________________________________________________
EventScaler::~EventScaler( void )
{
}

//______________________________________________________________________________
bool
EventScaler::ProcessingBegin( void )
{
  EventScaler::InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventScaler::ProcessingNormal( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");

  gRM.Decode();
  gScaler.Decode();

  event.evnum = gRM.EventNumber();

  static const int runnum = gRM.RunNumber();
  std::stringstream ss;
  ss << "in " << func_name << std::endl
     << std::left  << std::setw(16) << "RUN"
     << std::right << std::setw(16) << runnum << " : "
     << std::left  << std::setw(16) << "Event number"
     << std::right << std::setw(16) << event.evnum
     << std::endl;

  if( event.evnum%400==0 ) gScaler.Print( ss.str() );

#if SPILL_RESET
  if( gScaler.SpillIncrement() ) gScaler.Clear();
#endif

  return true;
}

//______________________________________________________________________________
void
EventScaler::InitializeEvent( void )
{
  event.evnum = -1;
}

//______________________________________________________________________________
bool
EventScaler::ProcessingEnd( void )
{
  return true;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventScaler;
}

//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
#if USE_COMMA
  gScaler.SetComma();
#endif
  HBTree( "scaler", "tree of Scaler" );

  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return true;
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  static const std::string func_name = "[ConfMan::FinalizeProcess()]";

  if( event.evnum==0 ) return true;

  const int run_number = gRM.RunNumber();
  std::stringstream ss;
  ss << "in " << func_name << std::endl
     << std::left  << std::setw(16) << "RUN"
     << std::right << std::setw(16) << run_number << " : "
     << std::left  << std::setw(16) << "Event number"
     << std::right << std::setw(16) << event.evnum
     << std::endl;

  gScaler.Print( ss.str() );

#if MAKE_LOG
  const std::string& bin_dir( hddaq::dirname(hddaq::selfpath()) );
  // assume symbolic link "data" to be used
  const std::string& data_dir( bin_dir+"/../data" );

  std::stringstream run_number_ss; run_number_ss << run_number;
  const std::string& recorder_log( data_dir+"/recorder.log" );
  std::ifstream ifs( recorder_log.c_str() );
  if( !ifs.is_open() ){
    std::cerr << "#E " << func_name << " "
	      << "cannot open recorder.log : "
	      << recorder_log << std::endl;
    return false;
  }

  const std::string& scaler_dir( bin_dir+"/../scaler" );
  const std::string& scaler_txt = Form( "%s/scaler_%05d.txt",
					scaler_dir.c_str(), run_number );

  std::ofstream ofs( scaler_txt.c_str() );
  if( !ofs.is_open() ){
    std::cerr << "#E " << func_name << " "
	      << "cannot open scaler.txt : "
	      << scaler_txt << std::endl;
    return false;
  }

  int recorder_event_number = 0;
  bool found_run_number = false;
  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() ) continue;
    std::istringstream input_line( line );
    std::istream_iterator<std::string> line_begin( input_line );
    std::istream_iterator<std::string> line_end;
    std::vector<std::string> log_column( line_begin, line_end );
    if( log_column[0]!="RUN" ) continue;
    if( log_column[1]!=run_number_ss.str() ) continue;
    recorder_event_number = hddaq::a2i(log_column[15]);
    ofs << line << std::endl;
    found_run_number = true;
  }

  if( !found_run_number ){
    std::cerr << "#E " << func_name << " "
	      << "not found run# " << run_number
	      << " in " << recorder_log << std::endl;
    return false;
  }

  ofs << std::endl;
  ofs << std::left  << std::setw(15) << "" << "\t"
      << std::right << std::setw(15) << "Integral" << std::endl;
  ofs << std::left  << std::setw(15) << "Event"    << "\t"
      << std::right << std::setw(15) << event.evnum << std::endl;

  if( recorder_event_number != event.evnum ){
    std::cerr << "#W " << func_name << " "
	      << "event number mismatch" << std::endl
	      << "   recorder : " << recorder_event_number << std::endl
	      << "   decode   : " << event.evnum << std::endl;
  }

  double spill     = gScaler.Get("Spill");
  double real_time = gScaler.Get("Real_Time");
  double live_time = gScaler.Get("Live_Time");
  double l1_req    = gScaler.Get("L1_Req");
  double l1_acc    = gScaler.Get("L1_Acc");
  // double l2_acc    = gScaler.Get("L2_Acc");
  double tm        = gScaler.Get("TM");
  double kbeam     = gScaler.Get("K_beam");
  double pibeam    = gScaler.Get("Pi_beam");
  double pbeam     = gScaler.Get("/p_beam");
  double bh1bh2    = gScaler.Get("BH1xBH2");
  double kk        = gScaler.Get("(K,K)PS");
  double real_live   = live_time/real_time;
  double daq_eff     = l1_acc/l1_req;
  // double l2_eff      = l2_acc/l1_acc;
  double krate       = kbeam/spill;
  double pirate      = pibeam/spill;
  double prate       = pbeam/spill;
  double bh12rate    = bh1bh2/spill;
  double duty_factor = daq_eff/(1.-daq_eff)*(1./real_live-1.);
  if( duty_factor >= 1. ) duty_factor = 1.;

  // DAQ Info
  for(int i=0; i<NumOfSegScaler; i++){
    std::string name = gScaler.GetName( 1, i );
    if( name=="n/a" ) continue;
    //if(name.substr(0,3)=="MsT") break;
    ofs << std::left  << std::setw(15) << name << "\t"
	<< std::right << std::setw(15) << gScaler.Get( 1, i ) << std::endl;
  }
  ofs << std::endl;

  // Counter Info
  for(int i=0; i<NumOfSegScaler; i++){
    std::string name = gScaler.GetName( 0, i );
    if( name=="n/a" ) continue;
    ofs << std::left  << std::setw(15) << name << "\t"
	<< std::right << std::setw(15) << gScaler.Get( 0, i ) << std::endl;
  }

  ofs << std::fixed << std::setprecision(6) << std::endl
      << std::left  << std::setw(18) << "Live/Real"      << "\t"
      << std::right << std::setw(12) << real_live        << std::endl
      << std::left  << std::setw(18) << "DAQ_Eff"        << "\t"
      << std::right << std::setw(12) << daq_eff          << std::endl
      << std::left  << std::setw(18) << "Duty_Factor"    << "\t"
      << std::right << std::setw(12) << duty_factor      << std::endl
      << std::left  << std::setw(18) << "K_beam/TM"      << "\t"
      << std::right << std::setw(12) << kbeam/tm         << std::endl
      << std::left  << std::setw(18) << "K_beam/BH1xBH2" << "\t"
      << std::right << std::setw(12) << kbeam/bh1bh2     << std::endl
      << std::left  << std::setw(18) << "K_beam/Pi_beam" << "\t"
      << std::right << std::setw(12) << kbeam/pibeam     << std::endl
      << std::left  << std::setw(18) << "(K,K)/K_beam"   << "\t"
      << std::right << std::setw(12) << kk/kbeam         << std::endl
      << std::setprecision(0)
      << std::left  << std::setw(18) << "K_beam/Spill"   << "\t"
      << std::right << std::setw(12) << krate            << std::endl
      << std::left  << std::setw(18) << "Pi_beam/Spill"  << "\t"
      << std::right << std::setw(12) << pirate           << std::endl
      << std::left  << std::setw(18) << "/p_beam/Spill"   << "\t"
      << std::right << std::setw(12) << prate            << std::endl
      << std::left  << std::setw(18) << "BH1xBH2/Spill"  << "\t"
      << std::right << std::setw(12) << bh12rate         << std::endl;

#endif

  return true;
}
