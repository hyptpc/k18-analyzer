// -*- C++ -*-

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"

namespace
{
using namespace root;
auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
}

//_____________________________________________________________________________
VEvent::VEvent()
{
}

//_____________________________________________________________________________
VEvent::~VEvent()
{
}

//_____________________________________________________________________________
class EventSkeleton : public VEvent
{
private:

public:
  EventSkeleton();
  ~EventSkeleton();
  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal();
  bool InitializeHistograms();
  void InitializeEvent();
};

//_____________________________________________________________________________
EventSkeleton::EventSkeleton()
  : VEvent()
{
}

//_____________________________________________________________________________
EventSkeleton::~EventSkeleton()
{
}

//_____________________________________________________________________________
struct Event
{
  Int_t runnum;
  Int_t evnum;
  Int_t spill;
  void clear(){
    runnum = -1;
    evnum = -1;
    spill = -1;
  }
};

//_____________________________________________________________________________
namespace root
{
Event  event;
TH1   *h[MaxHist];
TTree *tree;
}

//_____________________________________________________________________________
bool
EventSkeleton::ProcessingBegin()
{
  event.clear();
  return true;
}

//_____________________________________________________________________________
bool
EventSkeleton::ProcessingNormal()
{
  event.runnum = gUnpacker.get_run_number();
  event.evnum  = gUnpacker.get_event_number();
  // rawData = new RawData;
  // rawData->DecodeHits();

  // gRM.Decode();

  // event.runnum = gRM.RunNumber();
  // event.evnum  = gRM.EventNumber();
  // event.spill  = gRM.SpillNumber();

  // for(Int_t i=0; i<100; ++i){
  //   HF1(i, (double)i);
  // }

  return true;
}

//_____________________________________________________________________________
bool
EventSkeleton::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
void
EventSkeleton::InitializeEvent()
{
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new EventSkeleton;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeHistograms()
{
  for(Int_t i=0; i<100; ++i){
    HB1(i, Form("hist %d", i), 100, 0., 100.);
  }

  HBTree("skeleton","tree of Skeleton");
  tree->Branch("runnum", &event.runnum, "runnum/I");
  tree->Branch("evnum",  &event.evnum,  "evnum/I");
  tree->Branch("spill",  &event.spill,  "spill/I");

  HPrint();
  // gUnpacker.disable_istream_bookmark();
  return true;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")    &&
     InitializeParameter<HodoParamMan>("HDPRM"));
}

//_____________________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}
