// -*- C++ -*-

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include <std_ostream.hh>

#include <spdlog/spdlog.h>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "DeleteUtility.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

namespace
{
using hddaq::unpacker::GUnpacker;
auto& gConf     = ConfMan::GetInstance();
auto& gCounter  = debug::ObjectCounter::GetInstance();
auto& gUnpacker = GUnpacker::get_instance();
enum EArg
{
  kArgProcess,
  kArgConfFile,
  kArgInFile,
  kArgOutFile,
  kArgc
};
}

TROOT theROOT("k18analyzer", "k18analyzer");

//______________________________________________________________________________
int
main(int argc, char **argv)
{
  spdlog::set_level((spdlog::level::level_enum)SPDLOG_ACTIVE_LEVEL);
  spdlog::set_pattern("%^#%L%$ %v");
  // spdlog::set_pattern("%^%L%$ %v");
  // spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%^%L%$] %v");

  std::vector<TString> arg(argv, argv + argc);
  const TString& process = gSystem->BaseName(arg[kArgProcess]);
  if(argc!=kArgc){
    spdlog::info("Usage: {} [conf-file] [input-stream] [output-file]",
                 process.Data());
    return EXIT_SUCCESS;
  }

  debug::Timer timer("End of Analyzer");

  const TString& conf_file = arg[kArgConfFile];
  const TString& in_file   = arg[kArgInFile];
  const TString& out_file  = arg[kArgOutFile];

  spdlog::info("recreate root file : {}", out_file.Data());
  new TFile(out_file, "recreate");

  if (!gConf.Initialize(conf_file) || !gConf.InitializeUnpacker())
    return EXIT_FAILURE;

  gUnpacker.set_istream(in_file.Data());
  // gUnpacker.enable_istream_bookmark();
  gUnpacker.initialize();
  gConf.InitializeHistograms();
  gConf.WriteParameters();

  CatchSignal::Set(SIGINT);

  for(; !gUnpacker.eof() && !CatchSignal::Stop(); ++gUnpacker){
    ProcessBegin();
    ProcessNormal();
    ProcessEnd();
    gCounter.check();
  }
  gConf.Finalize();

  gFile->Write();
  gFile->Close();

  return EXIT_SUCCESS;
}
