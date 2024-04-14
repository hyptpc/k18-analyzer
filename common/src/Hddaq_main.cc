#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>

#include <filesystem_util.hh>
#include <std_ostream.hh>


#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "DeleteUtility.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

namespace
{
  using namespace hddaq::unpacker;
  const std::string& class_name("");
  ConfMan&              gConf     = ConfMan::GetInstance();
  debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();
  UnpackerManager&      gUnpacker = GUnpacker::get_instance();
}

//______________________________________________________________________________
int
hddaq_main()
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  debug::Timer timer(func_name+" End of Analyzer");
  new TFile( gConf.get_outfilename().c_str(), "recreate" );

  if( !gConf.Initialize() || !gConf.InitializeUnpacker() )
    return EXIT_FAILURE;

  gUnpacker.set_istream( gConf.get_infilename().c_str() );
  gUnpacker.set_parameter("max_loop", std::to_string(gConf.get_maxloop()));
  gUnpacker.set_parameter("skip",     std::to_string(gConf.get_skip()));
  gUnpacker.initialize();
  CatchSignal::Set(SIGINT);

  gConf.BeginRun();

  for( ; !gUnpacker.eof() && !CatchSignal::Stop(); ++gUnpacker ){
    VEvent* event = gConf.EventAllocator();
    bool tmp=event->ProcessingBegin() &&
      event->ProcessingNormal() &&
      event->ProcessingEnd();
    delete event;
    gCounter.check();
    if(!tmp) break;
  }
  gConf.Finalize();

  gFile->Write();
  gFile->Close();

  return EXIT_SUCCESS;
}
