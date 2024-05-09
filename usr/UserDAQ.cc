// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iostream>
#include <sstream>

#include <UnpackerManager.hh>
#include <Unpacker.hh>
#include <DAQNode.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "ScalerAnalyzer.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UserParamMan.hh"
#include "DCGeomMan.hh"

#include "RMAnalyzer.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"

namespace
{
  using namespace root;
  const auto qnan = TMath::QuietNaN();
  using hddaq::unpacker::DAQNode;
  auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
  auto& gRM = RMAnalyzer::GetInstance();
  auto& gUser = UserParamMan::GetInstance();
  auto& gScaler = ScalerAnalyzer::GetInstance();
}

//_____________________________________________________________________________
struct Event
{
  Int_t evnum;
  Int_t spill;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  // data size
  Int_t eb1_datasize;
  Int_t vme_datasize;
  Int_t ea0c_datasize;
  Int_t hul_datasize;
  Int_t aft_datasize;
  Int_t cobo_datasize;
  // trigger rate
  Int_t  L1req;	 
  Int_t  L1acc;	 
  Int_t  realtime;
  Int_t  livetime;
  // scaler 
  Int_t scaler[NumOfScaler][NumOfSegScaler];
  
  Double_t unixtime;
  Int_t clock10M;

  void clear();
};

//_____________________________________________________________________________
void
Event::clear()
{
  evnum       = 0;
  spill       = 0;
  L1req       = 0;	 
  L1acc       = 0;	 
  realtime    = 0;
  livetime    = 0;
  eb1_datasize = 0;
  vme_datasize = 0;
  ea0c_datasize = 0;
  hul_datasize = 0;
  aft_datasize = 0;
  cobo_datasize = 0;
  unixtime      = -0.999;
  clock10M      = 0;

  for(Int_t i=0; i<NumOfScaler; ++i){
    for(Int_t it=0; it<NumOfSegScaler; ++it){
      scaler[i][it]=0;
    }
  }

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }
}

//_____________________________________________________________________________
struct Dst
{
  Int_t evnum;
  Int_t spill;

  Int_t eb1_datasize;
  Int_t vme_datasize;
  Int_t ea0c_datasize;
  Int_t hul_datasize;
  Int_t aft_datasize;
  Int_t cobo_datasize;

  Int_t L1req;
  Int_t L1acc;
  Int_t realtime;
  Int_t livetime;

  Int_t trigpat[NumOfSegTrig];
  Int_t trigflag[NumOfSegTrig];

  Double_t unixtime;

  // scaler 
  Int_t scaler[NumOfScaler][NumOfSegScaler];

  void clear();
};

//_____________________________________________________________________________
void
Dst::clear()
{
  evnum    = 0;
  spill    = 0;

  L1req       = 0;	 
  L1acc       = 0;	 
  realtime    = 0;
  livetime    = 0;
  eb1_datasize = 0;
  vme_datasize = 0;
  ea0c_datasize = 0;
  hul_datasize = 0;
  aft_datasize = 0;
  cobo_datasize = 0;
  unixtime      = -0.999;

  for(Int_t i=0; i<NumOfScaler; ++i){
    for(Int_t it=0; it<NumOfSegScaler; ++it){
      scaler[i][it]=0;
    }
  }

  for(Int_t it=0; it<NumOfSegTrig; ++it){
    trigpat[it]  = -1;
    trigflag[it] = -1;
  }

}


//_____________________________________________________________________________
namespace root
{
Event  event;
Dst    dst;
TH1*   h[MaxHist];
TTree* tree;
TTree* daq;
enum eDetHid {
  BH1Hid    = 10000,
  BH2Hid    = 20000,
  BACHid    = 30000,
  TOFHid    = 60000,
  AC1Hid    = 70000,
  WCHid     = 80000,

  DAQHid    = 90000
};
}

//_____________________________________________________________________________
Bool_t
ProcessingBegin()
{
  event.clear();
  dst.clear();
  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingNormal()
{
  RawData rawData;
  
  gRM.Decode();
  
  // event.evnum = gRM.EventNumber();
  event.evnum = gUnpacker.get_event_number();
  event.spill = gRM.SpillNumber();
  //  dst.evnum   = gRM.EventNumber();
  dst.evnum   = gUnpacker.get_event_number();
  dst.spill   = gRM.SpillNumber();

  // data size
  ///// DAQ
  //___ node id
  static const Int_t k_eb = gUnpacker.get_fe_id("k18eb");
  std::vector<Int_t> vme_fe_id;
  std::vector<Int_t> hul_fe_id;
  std::vector<Int_t> ea0c_fe_id;
  std::vector<Int_t> cobo_fe_id;
  std::vector<Int_t> aft_fe_id;
  for(auto&& c : gUnpacker.get_root()->get_child_list()) {
    if (!c.second) continue;
    TString n = c.second->get_name();
    auto id = c.second->get_id();
    if (n.Contains("vme"))
      vme_fe_id.push_back(id);
    if (n.Contains("hul"))
      hul_fe_id.push_back(id);
    if (n.Contains("easiroc"))
      ea0c_fe_id.push_back(id);
    if (n.Contains("cobo"))
      cobo_fe_id.push_back(id);
    if (n.Contains("aft"))
      aft_fe_id.push_back(id);
  }
  
  { //___ EB
    auto data_size = gUnpacker.get_node_header(k_eb, DAQNode::k_data_size);
    event.eb1_datasize += data_size;
    HF1(DAQHid+100+1, data_size);
  }

  { //___ VME
    for(Int_t i=0, n=vme_fe_id.size(); i<n; ++i) {
      auto data_size = gUnpacker.get_node_header(vme_fe_id[i], DAQNode::k_data_size);
      event.vme_datasize += data_size;
      //      HF1(DAQHid+200+i+1, data_size);
    }
  }

  { // EASIROC node
    for(Int_t i=0, n=ea0c_fe_id.size(); i<n; ++i) {
      auto data_size = gUnpacker.get_node_header(ea0c_fe_id[i], DAQNode::k_data_size);
      event.ea0c_datasize += data_size;
      //      HF1(DAQHid+300+i+1, data_size);
    }
  }

  { //___ HUL node
    for(Int_t i=0, n=hul_fe_id.size(); i<n; ++i) {
      auto data_size = gUnpacker.get_node_header(hul_fe_id[i], DAQNode::k_data_size);
      event.hul_datasize += data_size;
      //      HF1(DAQHid+400+i+1, data_size);
    }
  }

  { //___ AFT node
    for(Int_t i=0, n=aft_fe_id.size(); i<n; ++i) {
      auto data_size = gUnpacker.get_node_header(aft_fe_id[i], DAQNode::k_data_size);
      event.aft_datasize += data_size;
      //      HF1(DAQHid+500+i+1, data_size);
    }
  }
  {
  //  Int_t c = ScalerAnalyzer::kLeft;
  Int_t c = 0;
  Int_t r = 0;
  gScaler.Decode();
  gScaler.Set(c, r++, ScalerInfo("Spill",       -1,  -1));
  gScaler.Set(c, r++, ScalerInfo("L1-Req",        0,   3));
  gScaler.Set(c, r++, ScalerInfo("L1-Acc",        0,   4));
  gScaler.Set(c, r++, ScalerInfo("Real-Time",     0,   1));
  gScaler.Set(c, r++, ScalerInfo("Live-Time",     0,   2));

  event.L1req = gScaler.Get("L1-Req");
  event.L1acc = gScaler.Get("L1-Acc");
  event.realtime = gScaler.Get("Real-Time");
  event.livetime = gScaler.Get("Live-Time");

  HF1(DAQHid+900+1, event.L1req);
  HF1(DAQHid+900+2, event.L1acc);
  HF1(DAQHid+900+3, event.realtime);
  HF1(DAQHid+900+4, event.livetime);
  }
  HF1(1, 0);

  // Trigger Flag
  rawData.DecodeHits("TFlag");
  std::bitset<NumOfSegTrig> trigger_flag;
  for(const auto& hit: rawData.GetHodoRawHC("TFlag")){
    Int_t seg = hit->SegmentId();
    Int_t tdc = hit->GetTdc();
    if(tdc > 0){
      event.trigpat[trigger_flag.count()] = seg;
      event.trigflag[seg] = tdc;
      dst.trigpat[trigger_flag.count()] = seg;
      dst.trigflag[seg] = tdc;
      trigger_flag.set(seg);
      HF1(10, seg);
      HF1(10+seg, tdc);
    }
  }

  // UnixTime
  {
    int nhit  = 0;
    int plane = 0;
    int seg   = 0;
    nhit = gUnpacker.get_entries( DetIdUnixTime, plane, 0, seg, 0 );
    if( nhit>0 ){
      Long_t lo = gUnpacker.get( DetIdUnixTime, plane, 0, seg, 0);
      Long_t hi = gUnpacker.get( DetIdUnixTime, plane, 0, seg, 1);
      Double_t ut = Double_t(hi*0xfffffff + lo)/1e6;
      event.unixtime = ut;
      dst.unixtime= ut;
    }
  }
  {
    for(int plane=0; plane<NumOfScaler; ++plane){
      for(int seg=0; seg<NumOfSegScaler; ++seg){
	Int_t nhit = gUnpacker.get_entries( DetIdScaler, plane, 0, seg, 0 );
	if( nhit>0 ){
	  Int_t data = gUnpacker.get( DetIdScaler, plane, 0, seg, 0);
	  event.scaler[plane][seg] = data;
	  dst.scaler[plane][seg] = data;
	  if(plane==2 and seg==1){
	    event.clock10M=data;
	  }
	}
      }
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd]){
    // std::cout << "L1-Acc: " << gScaler.Get("L1-Acc") << std::endl;
    return true;
  }
  HF1(1, 1); 

  return true;
}

//_____________________________________________________________________________
Bool_t
ProcessingEnd()
{
  tree->Fill();
  daq->Fill();
  return true;
}

//_____________________________________________________________________________
namespace
{
const Int_t    NbinAdc = 4096;
const Double_t MinAdc  =    0.;
const Double_t MaxAdc  = 4096.;

const Int_t    NbinTdc = 4096;
const Double_t MinTdc  =    0.;
const Double_t MaxTdc  = 4096.;

const Int_t    NbinTdcHr = 1e6/10;
const Double_t MinTdcHr  =  0.;
const Double_t MaxTdcHr  = 1e6;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);
  HB1(10, "Trigger HitPat", NumOfSegTrig, 0., Double_t(NumOfSegTrig));
  for(Int_t i=0; i<NumOfSegTrig; ++i){
    HB1(10+i+1, Form("Trigger Trig %d", i+1), 0x1000, 0, 0x1000);
  }
  HB1(DAQHid+100+1, "eb1 datasize", 10000, 0., 100000.);
  
  HB1(DAQHid+900+1, "L1 request", 10000, 0., 10000.);
  HB1(DAQHid+900+2, "L1 accept", 10000, 0., 10000.);
  HB1(DAQHid+900+3, "realtime", 10000, 0., 10000.);
  HB1(DAQHid+900+4, "livetime", 10000, 0., 10000.);
  HB2(DAQHid+900+5, "data size of vme-easiroc nodes", NumOfSegVMEEASIROC, 0., NumOfSegVMEEASIROC, 370, 90, 460);

  ////////////////////////////////////////////
  //Tree
  HBTree("tree","tree of Counter");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("spill",     &event.spill,     "spill/I");
  tree->Branch("eb1_datasize",     &event.eb1_datasize,     "eb1_datasize/I");
  tree->Branch("vme_datasize",     &event.vme_datasize,     "vme_datasize/I");
  tree->Branch("ea0c_datasize",     &event.ea0c_datasize,     "ea0c_datasize/I");
  tree->Branch("hul_datasize",     &event.hul_datasize,     "hul_datasize/I");
  tree->Branch("aft_datasize",     &event.aft_datasize,     "aft_datasize/I");
  tree->Branch("L1req",     &event.L1req,     "L1req/I");
  tree->Branch("L1acc",     &event.L1acc,    "L1acc/I");
  tree->Branch("realtime",     &event.realtime,     "realtime/I");
  tree->Branch("livetime",     &event.livetime,     "livetime/I");
  tree->Branch("unixtime",     &event.unixtime,     "UnixTime/D");
  tree->Branch("clock10M",     &event.clock10M,     "clock10M/I");

  //Trig
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  //Scaler
  tree->Branch("scaler",   event.scaler,  Form("scaler[%d][%d]/I", NumOfScaler, NumOfSegScaler));

  ////////////////////////////////////////////
  //Dst
  daq = new TTree("daq","Data Summary Table of Hodoscope");
  daq->Branch("evnum",     &dst.evnum,     "evnum/I");
  daq->Branch("spill",     &dst.spill,     "spill/I");
  daq->Branch("eb1_datasize",     &event.eb1_datasize,     "eb1_datasize/I");
  daq->Branch("vme_datasize",     &event.vme_datasize,     "vme_datasize/I");
  daq->Branch("ea0c_datasize",     &event.ea0c_datasize,     "ea0c_datasize/I");
  daq->Branch("hul_datasize",     &event.hul_datasize,     "hul_datasize/I");
  daq->Branch("aft_datasize",     &event.aft_datasize,     "aft_datasize/I");
  daq->Branch("L1req",     &event.L1req,     "L1req/I");
  daq->Branch("L1acc",     &event.L1acc,    "L1acc/I");
  daq->Branch("realtime",     &event.realtime,     "realtime/I");
  daq->Branch("livetime",     &event.livetime,     "livetime/I");
  daq->Branch("unixtime",     &event.unixtime,     "UnixTime/D");

  daq->Branch("trigpat",    dst.trigpat,   Form("trigpat[%d]/I", NumOfSegTrig));
  daq->Branch("trigflag",   dst.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));
  daq->Branch("scaler",   dst.scaler,  Form("scaler[%d][%d]/I", NumOfScaler, NumOfSegScaler));

  HPrint();
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO") &&
     InitializeParameter<HodoParamMan>("HDPRM") &&
     InitializeParameter<HodoPHCMan>("HDPHC")   &&
     InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
Bool_t
ConfMan::FinalizeProcess()
{
  return true;
}
