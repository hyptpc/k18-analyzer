/**
 *  file: UserEMC.cc
 *  date: 2017.04.10
 *
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "EMCAnalyzer.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "KuramaLib.hh"
#include "VEvent.hh"

namespace
{
  using namespace root;
  const std::string& class_name("EMC");
  RMAnalyzer&  gRM  = RMAnalyzer::GetInstance();
  EMCAnalyzer& gEmc = EMCAnalyzer::GetInstance();
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
class EventEMC : public VEvent
{
private:

public:
        EventEMC( void );
       ~EventEMC( void );
  bool ProcessingBegin( void );
  bool ProcessingEnd( void );
  bool ProcessingNormal( void );
  bool InitializeHistograms( void );
  void InitializeEvent( void );
};

//______________________________________________________________________________
EventEMC::EventEMC( void )
  : VEvent()
{
}

//______________________________________________________________________________
EventEMC::~EventEMC( void )
{
}

//______________________________________________________________________________
struct Event
{
  int    runnum;
  int    evnum;
  int    spill;
  int    serial;
  double xpos;
  double ypos;
  int    state;
  unsigned long long time;
};

//______________________________________________________________________________
namespace root
{
  Event  event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//______________________________________________________________________________
bool
EventEMC::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventEMC::ProcessingNormal( void )
{
  gRM.Decode();
  gEmc.Decode();

  event.runnum = gRM.RunNumber();
  event.evnum  = gRM.EventNumber();
  event.spill  = gRM.SpillNumber();
  event.serial = gEmc.Serial();
  event.xpos   = gEmc.X();
  event.ypos   = gEmc.Y();
  event.state  = gEmc.State();
  event.time   = gEmc.Time();

  HF1( 11, event.serial );
  HF1( 12, event.xpos );
  HF1( 13, event.ypos );
  HF1( 14, event.state );
  HF1( 15, event.time );
  HF2( 20, event.xpos, event.ypos );

  if( gRM.LocalEventNumber()==1 ){ // first event of each spill
    HF2( 21, event.spill, event.serial );
    HF2( 22, event.spill, event.xpos );
    HF2( 23, event.spill, event.ypos );
    HF2( 24, event.spill, event.state );
  }

  if( event.state==1 ){
    hddaq::cout << "#D Event Number " << std::setw(6) << event.evnum
		<< "  EMC state : \e[35;1m" << event.state << "\e[m" << std::endl;
  // } else {
  //   hddaq::cout << "#D Event Number " << std::setw(6) << event.evnum
  // 		<< "  EMC state : \e[32;1m" << event.state << "\e[m" << std::endl;
  }

  return true;
}

//______________________________________________________________________________
bool
EventEMC::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
void
EventEMC::InitializeEvent( void )
{
  event.runnum = 0;
  event.evnum  = 0;
  event.spill  = 0;
  event.serial = 0;
  event.xpos   = -999.;
  event.ypos   = -999.;
  event.state  = -1;
  event.time   = -1;
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventEMC;
}

//______________________________________________________________________________
namespace
{
  const int NbinSpill   = 1000;
  const double MinSpill =    0.;
  const double MaxSpill = 1000.;

  const int NbinSerial   = 200;
  const double MinSerial =   0.;
  const double MaxSerial = 1e6;

  const int NbinXpos   =  200;
  const double MinXpos = -200.;
  const double MaxXpos =  200.;

  const int NbinYpos   =  200;
  const double MinYpos = -200.;
  const double MaxYpos =  200.;

  const int NbinState   = 5;
  const double MinState = 0.;
  const double MaxState = 5.;

  const int NbinTime   = 500;
  const double MinTime = 1.5e6;
  const double MaxTime = 2.0e6;
}

//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );
  HB1( 11, "EMC serial", NbinSerial, MinSerial, MaxSerial );
  HB1( 12, "EMC xpos",   NbinXpos,   MinXpos,   MaxXpos   );
  HB1( 13, "EMC ypos",   NbinYpos,   MinYpos,   MaxYpos   );
  HB1( 14, "EMC state",  NbinState,  MinState,  MaxState  );
  HB1( 15, "EMC time",   NbinTime,   MinTime,   MaxTime   );

  HB2( 21, "EMC serial%spill",
       NbinSpill, MinSpill, MaxSpill,
       NbinSerial, MinSerial, MaxSerial );
  HB2( 22, "EMC xpos%spill",
       NbinSpill, MinSpill, MaxSpill,
       NbinXpos,   MinXpos,   MaxXpos   );
  HB2( 23, "EMC ypos%spill",
       NbinSpill, MinSpill, MaxSpill,
       NbinYpos,   MinYpos,   MaxYpos   );
  HB2( 24, "EMC state%spill",
       NbinSpill, MinSpill, MaxSpill,
       NbinState,  MinState,  MaxState  );
  HB2( 25, "EMC time%spill",
       NbinSpill, MinSpill, MaxSpill,
       NbinTime,   MinTime,   MaxTime   );

  // h15 = new TH1F( "h15", "EMC time",   NbinTime,   MinTime,   MaxTime   );
  // h15->GetXaxis()->SetTimeDisplay(1);
  // h15->GetXaxis()->SetLabelOffset(0.02);
  // h15->GetXaxis()->SetTimeFormat("#splitline{%Y/%m/%d}{  %H:%M:%S}");
  // h15->GetXaxis()->SetTimeOffset(0,"jpn");

  HB2( 20, "EMC XY Position",
       NbinXpos, MinXpos, MaxXpos,
       NbinYpos, MinYpos, MaxYpos );

  ////////////////////////////////////////////
  //Tree
  HBTree( "emc","tree of EMC" );
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");
  tree->Branch("serial", &event.serial, "serial/I");
  tree->Branch("xpos",   &event.xpos,   "xpos/D");
  tree->Branch("ypos",   &event.ypos,   "ypos/D");
  tree->Branch("state",  &event.state,  "state/I");
  tree->Branch("time",   &event.time,   "time/l");

  HPrint();
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
  return true;
}
