// -*- C++ -*-

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include <filesystem_util.hh>
#include <std_ostream.hh>

#include "CatchSignal.hh"
#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "DeleteUtility.hh"
#include "VEvent.hh"
#include "ReslMapMan.hh"

#include "RootData.hh"

namespace
{
  const std::string& class_name("");
  ConfMan&              gConf     = ConfMan::GetInstance();
  debug::ObjectCounter& gCounter  = debug::ObjectCounter::GetInstance();
}

//______________________________________________________________________________
int
knucl_main()
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  debug::Timer timer(func_name+" End of Analyzer");

  TFile *f = new TFile( gConf.get_infilename().c_str() );
  TTree *evtree = (TTree*)f->Get("tree");

  EventHeaderMC* header = 0;
  DetectorData *detectorData = 0;
  MCData* mcData = 0;
  ReactionData* reactionData=0;

  // evtree->SetBranchAddress( "EventHeaderMC", &header );
  // evtree->SetBranchAddress( "DetectorData", &detectorData );
  // evtree->SetBranchAddress( "MCData", &mcData );
  // evtree->SetBranchAddress( "ReactionData", &reactionData );

  TFile *fout=new TFile( gConf.get_outfilename().c_str(), "recreate" );
  if( !gConf.Initialize() )
    return EXIT_FAILURE;
  if( !gConf.InitializeParameter<ReslMapMan>("RESOL") )
    return EXIT_FAILURE;

  CatchSignal::Set(SIGINT);
  gConf.BeginRun();

  int nev = evtree->GetEntries();
  std::cout << " AllEvent : " << nev << std::endl;
  std::cout << "|0%                  |50%                |100% "
            <<   std::endl;
  int moniter=0;
  for( int iev=0; iev<nev && !CatchSignal::Stop(); iev++ ){
    if( iev==0 )
      std::cout << "|";
    if( (nev>40) && (iev%(nev/40)==0) )
      std::cout << "*"<<std::flush;
    if( iev%(100)==0 )
      {
        if(moniter==0) std::cout << "\\\b"<<std::flush;
        else if(moniter==1) std::cout << "-\b"<<std::flush;
        else if(moniter==2) std::cout << "/\b"<<std::flush;
        else if(moniter==3) std::cout << "|\b"<<std::flush;
        moniter++;
        if(moniter==4) moniter=0;
      }
    if( iev==nev-1 )
      std::cout << "| fin!" << std::endl;

    if(gConf.get_maxloop()>0&&iev>gConf.get_maxloop())
      break;
    if(gConf.get_skip()>0&&iev<gConf.get_skip())
      continue;

    evtree->GetEvent(iev);
    VEvent* event = gConf.EventAllocator();
    bool tmp2=
      event->ProcessingBeginMC(detectorData,mcData,reactionData,iev) &&
      event->ProcessingNormal() &&
      event->ProcessingEnd();
    delete event;
    gCounter.check();
    if(!tmp2) break;
  }
  gConf.Finalize();

  fout->Write();
  fout->Close();
  f->Close();
  return EXIT_SUCCESS;
}
