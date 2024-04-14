#include <iostream>
#include <sstream>
#include <cmath>
#include "TString.h"

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RootHelper.hh"
#include "MTDCRawHit.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"
#include "RawData.hh"
#include "CDCWireMapMan.hh"
#include "VEvent.hh"
#include "HistTools.hh"

#include "UnpackerManager.hh"
#include "DAQNode.hh"

#include "setup.hh"

#define DEBUG 0
namespace
{
  enum eType { kTag, kLock, nType };
  enum eData { kEvent, kSpill, kSerial, kTime, nData };
  using namespace e73_2024;
  using namespace root;
  using namespace hddaq::unpacker;
  using namespace hddaq;
  const std::string& classname("EventRaw");
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
class EventBeam : public VEvent
{
private:

public:
  EventBeam( void );
  ~EventBeam( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventBeam::EventBeam( void )
  : VEvent()
{
}

//______________________________________________________________________________
EventBeam::~EventBeam( void )
{
}

//______________________________________________________________________________
namespace root
{
}

//______________________________________________________________________________
bool
EventBeam::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventBeam::ProcessingNormal( void )
{
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  // DAQ
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
  //  gUnpacker.print_tag();
  //  gUnpacker.show_summary();
  int evnum=gUnpacker.get_event_number();
  static const int k_eb      = gUnpacker.get_fe_id("k18eb");
  static const int k_vme     = gUnpacker.get_fe_id("vme_qdc1");
  static const int k_sca     = gUnpacker.get_fe_id("hulscaler-139");
  static const int k_flag    = gUnpacker.get_fe_id("hulmhtdc-138");
  static const int k_hr      = gUnpacker.get_fe_id("hulhrtdc-121");
  static const int k_cdc1    = 1650;
  static const int k_cdc2    = 1660;
  static const int k_bldc    = gUnpacker.get_fe_id("hulmhtdc-101");
  {
    int data_size = gUnpacker.get_node_header( k_eb, DAQNode::k_data_size);
    hist::H1("DataSize",data_size,10000,0,20000);
    hist::H2("DataSize_vs_EvNum",evnum,data_size,300,0,3e6,200,0,20000);
  }
  { // VME node
    const char* name="VME-RM";
    const int k_device = gUnpacker.get_device_id(name);
    for( int i=0; i<2; ++i ){
      int node_id = k_vme+i;
      int data_size = gUnpacker.get_node_header( node_id, DAQNode::k_data_size);
      hist::H1(Form("DataSize_node%04d",node_id),data_size,1000,0,1000);
      hist::H2(Form("DataSize_node%04d_vs_EvNum",node_id),
	       evnum,data_size,300,0,3e6,100,0,2000);
      hist::H2(Form("DataSize_vme"),
	       i,data_size,2,-0.5,1.5,1000,0,2000);
      int tmp_evnum = gUnpacker.get_node_header( node_id, DAQNode::k_event_number);
      // gUnpacker.dump_data_fe(node_id);
      int spill=-1;
      int evnum=-1;
      for( int ch=0; ch<nType; ++ch ){
        for( int data=0; data<nData; ++data ){
          int nhit = gUnpacker.get_entries( k_device, i, 0, ch, data );
          if( nhit<=0 ) continue;
          unsigned int val = gUnpacker.get( k_device, i, 0, ch, data );
          if( ch==kLock && val!=0 ) val = 1;
          if( ch==kTag && data==kSpill ) spill=val;
          if( ch==kTag && data==kEvent ) evnum=val-1;
        }
      }
#if 0
      std::cout<<std::setw(10)<<name
               <<std::setw(10)<<"plane"
               <<std::setw(5)<<i<<" : "
               <<std::setw(10)<<"spill"
               <<std::setw(8)<<spill<<" : "
               <<std::setw(15)<<"Event_number"
               <<std::setw(8)<<evnum
               <<std::setw(8)<<tmp_evnum
               <<std::endl;
#endif
    }
  }

  { // HUL
    for( int i=0; i<2; ++i ){
      int node_id = k_bldc+i;
      int data_size = gUnpacker.get_node_header( node_id, DAQNode::k_data_size);
      hist::H1(Form("DataSize_node%04d",node_id),data_size,1000,0,1000);
      hist::H2(Form("DataSize_node%04d_vs_EvNum",node_id),
	       evnum,data_size,300,0,3e6,100,0,2000);
      hist::H2(Form("DataSize_hul"),
	       i,data_size,50,-0.5,49.5,1000,0,2000);
      int tmp_evnum = gUnpacker.get_node_header( node_id, DAQNode::k_event_number);
      //      gUnpacker.dump_data_fe(node_id);
    }
    {
      const char* name="HUL-RM";
      const int k_device = gUnpacker.get_device_id(name);
      for(int seg=0;seg<2;seg++){
	int spill=-1;
	int evnum=-1;
	int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, 4);
	if(nhit==1){
	  spill = gUnpacker.get(k_device, 0, seg, 0, 4);
	}
	nhit = gUnpacker.get_entries(k_device, 0, seg, 0, 3);
	if(nhit==1){
	  evnum = gUnpacker.get(k_device, 0, seg, 0, 3);
	}
#if 0
	std::cout<<std::setw(10)<<name
		 <<std::setw(10)<<"seg"
		 <<std::setw(5)<<seg<<" : "
		 <<std::setw(10)<<"spill"
		 <<std::setw(8)<<spill<<" : "
		 <<std::setw(15)<<"Event_number"
		 <<std::setw(8)<<evnum
	  //		 <<std::setw(8)<<tmp_evnum
		 <<std::endl;
#endif
      }
    }
  }
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  return true;
}

//______________________________________________________________________________
bool
EventBeam::ProcessingEnd( void )
{
  return true;
}

//______________________________________________________________________________
void
EventBeam::InitializeEvent( void )
{
}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventBeam;
}
//______________________________________________________________________________
bool
ConfMan::InitializeHistograms( void )
{
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return ( InitializeParameter<CDCWireMapMan>("CDCGeom", "CDCASD" ) );
}

//______________________________________________________________________________
bool
ConfMan::BeginRunProcess()
{
  return true;
}
//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
