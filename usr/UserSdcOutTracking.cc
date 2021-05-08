// -*- C++ -*-

#include "VEvent.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <UnpackerManager.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "RMAnalyzer.hh"
#include "RootHelper.hh"
#include "KuramaLib.hh"
#include "RawData.hh"

#define HodoCut     0
#define TotCut      0
#define Chi2Cut     1
#define MaxMultiCut 0
#define UseTOF      1 // use or not TOF for tracking

namespace
{
using namespace root;
const auto& gGeom = DCGeomMan::GetInstance();
auto& gRM = RMAnalyzer::GetInstance();
const auto& gUser = UserParamMan::GetInstance();
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const double& zTOF = gGeom.LocalZ("TOF");
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
class EventSdcOutTracking : public VEvent
{
private:
  RawData*      rawData;
  DCAnalyzer*   DCAna;
  HodoAnalyzer* hodoAna;

public:
  static TString ClassName();
        EventSdcOutTracking();
       ~EventSdcOutTracking();
  bool  ProcessingBegin();
  bool  ProcessingEnd();
  bool  ProcessingNormal();
  bool  InitializeHistograms();
  void  InitializeEvent();
};

//_____________________________________________________________________________
TString
EventSdcOutTracking::ClassName()
{
  static TString s_name("EventSdcOutTracking");
  return s_name;
}

//_____________________________________________________________________________
EventSdcOutTracking::EventSdcOutTracking()
  : VEvent(),
    rawData(new RawData),
    DCAna(new DCAnalyzer),
    hodoAna(new HodoAnalyzer)
{
}

//_____________________________________________________________________________
EventSdcOutTracking::~EventSdcOutTracking()
{
  if(DCAna)   delete DCAna;
  if(hodoAna) delete hodoAna;
  if(rawData) delete rawData;
}

//_____________________________________________________________________________
bool
EventSdcOutTracking::ProcessingBegin()
{
  InitializeEvent();
  return true;
}

//_____________________________________________________________________________
struct Event
{
  int evnum;

  int trignhits;
  int trigpat[MaxHits];
  int trigflag[NumOfSegTrig];

  int nhBh1;
  double tBh1[MaxHits];
  double deBh1[MaxHits];

  int nhBh2;
  double tBh2[MaxHits];
  double deBh2[MaxHits];
  double Bh2Seg[MaxHits];

  double Time0Seg;
  double deTime0;
  double Time0;
  double CTime0;

  double btof;
  double stof[MaxHits];

  int nhTof;
  double TofSeg[MaxHits];
  double tTof[MaxHits];
  double dtTof[MaxHits];
  double deTof[MaxHits];

  int nhit[NumOfLayersSdcOut+2];
  int nlayer;
  double pos[NumOfLayersSdcOut+2][MaxHits];

  int ntrack;
  double chisqr[MaxHits];
  double x0[MaxHits];
  double y0[MaxHits];
  double u0[MaxHits];
  double v0[MaxHits];
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
EventSdcOutTracking::ProcessingNormal()
{
#if HodoCut
  static const double MinDeBH2   = gUser.GetParameter("DeBH2", 0);
  static const double MaxDeBH2   = gUser.GetParameter("DeBH2", 1);
  static const double MinDeBH1   = gUser.GetParameter("DeBH1", 0);
  static const double MaxDeBH1   = gUser.GetParameter("DeBH1", 1);
  static const double MinBeamToF = gUser.GetParameter("BTOF",  1);
  static const double MaxBeamToF = gUser.GetParameter("BTOF",  1);
#endif
  static const double MinDeTOF   = gUser.GetParameter("DeTOF",   0);
  static const double MaxDeTOF   = gUser.GetParameter("DeTOF",   1);
  static const double MinTimeTOF = gUser.GetParameter("TimeTOF", 0);
  static const double MaxTimeTOF = gUser.GetParameter("TimeTOF", 1);
  static const double StopTimeDiffSdcOut      = gUser.GetParameter("StopTimeDiffSdcOut",   0);
  static const double MinStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 0);
  static const double MaxStopTimingSdcOut = gUser.GetParameter("StopTimingSdcOut", 1);
#if TotCut
  static const double MinTotSDC3 = gUser.GetParameter("MinTotSDC3", 0);
  static const double MinTotSDC4 = gUser.GetParameter("MinTotSDC4", 0);
#endif
#if MaxMultiCut
  static const double MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");
#endif

  rawData->DecodeHits();

  gRM.Decode();

  event.evnum = gRM.EventNumber();

  // Trigger Flag
  std::bitset<NumOfSegTrig> trigger_flag;
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int trignhits = 0;
    int nh=cont.size();
    for(int i=0; i<nh; ++i){
      HodoRawHit *hit=cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if(tdc>0){
	trigger_flag.set(seg-1);
	event.trigpat[trignhits++] = seg;
	event.trigflag[seg-1]      = tdc;
      }
    }
    event.trignhits = trignhits;
  }

  HF1(1, 0.);

  if(trigger_flag[trigger::kSpillEnd] ||
      trigger_flag[trigger::kL1SpillOff])
    return true;

  HF1(1, 1.);

 //////////////BH2 time 0
  hodoAna->DecodeBH2Hits(rawData);
  int nhBh2 = hodoAna->GetNHitsBH2();
  event.nhBh2 = nhBh2;
#if HodoCut
  if(nhBh2==0) return true;
#endif
  HF1(1, 2);

  double time0 = -9999.;
  //////////////BH2 Analysis
  for(int i=0; i<nhBh2; ++i){
    BH2Hit* hit = hodoAna->GetHitBH2(i);
    if(!hit) continue;
    double seg = hit->SegmentId()+1;
    double cmt = hit->CMeanTime();
    double dE  = hit->DeltaE();

#if HodoCut
    if(dE<MinDeBH2 || MaxDeBH2<dE) continue;
#endif
    event.tBh2[i]   = cmt;
    event.deBh2[i]  = dE;
    event.Bh2Seg[i] = seg;
  }

  BH2Cluster *cl_time0 = hodoAna->GetTime0BH2Cluster();
  if(cl_time0){
    event.Time0Seg = cl_time0->MeanSeg()+1;
    event.deTime0  = cl_time0->DeltaE();
    event.Time0    = cl_time0->Time0();
    event.CTime0   = cl_time0->CTime0();
    time0          = cl_time0->CTime0();
  } else {
#if HodoCut
    return true;
#endif
  }

  HF1(1, 3.);

  //////////////BH1 Analysis
  hodoAna->DecodeBH1Hits(rawData);
  int nhBh1 = hodoAna->GetNHitsBH1();
  event.nhBh1 = nhBh1;
#if HodoCut
  if(nhBh1==0) return true;
#endif
  HF1(1, 4);

  for(int i=0; i<nhBh1; ++i){
    Hodo2Hit *hit = hodoAna->GetHitBH1(i);
    if(!hit) continue;
    double cmt = hit->CMeanTime();
    double dE  = hit->DeltaE();
#if HodoCut
    if(dE<MinDeBH1 || MaxDeBH1<dE) continue;
    if(btof<MinBeamToF || MaxBeamToF<btof) continue;
#endif
    event.tBh1[i]  = cmt;
    event.deBh1[i] = dE;
  }

  double btof0 = -999.;
  HodoCluster* cl_btof0 = event.Time0Seg > 0 ?
    hodoAna->GetBtof0BH1Cluster(event.CTime0) : nullptr;
  if(cl_btof0) btof0 = cl_btof0->CMeanTime() - time0;
  event.btof = btof0;

  HF1(1, 5.);


  //////////////Tof Analysis
  HodoClusterContainer TOFCont;
  hodoAna->DecodeTOFHits(rawData);
  hodoAna->TimeCutTOF(7, 25);
  int nhTof = hodoAna->GetNClustersTOF();
#if HodoCut
  if(nhTof!=0) return true;
#endif
  event.nhTof = nhTof;
  {
    int nhOk = 0;
    for(int i=0; i<nhTof; ++i){
      HodoCluster *hit = hodoAna->GetClusterTOF(i);
      if(!hit) continue;
      double cmt  = hit->CMeanTime();
      double dt   = hit->TimeDif();
      double de   = hit->DeltaE();
      double stof = cmt-time0;
      event.TofSeg[i] = hit->MeanSeg()+1;
      event.tTof[i]   = cmt;
      event.dtTof[i]  = dt;
      event.deTof[i]  = de;
      event.stof[i]   = stof;
      TOFCont.push_back(hit);
      if(MinDeTOF<de  && de<MaxDeTOF  &&
	  MinTimeTOF<stof && stof<MaxTimeTOF){
	++nhOk;
      }
    }
#if HodoCut
    if(nhOk==0) return true;
#endif
  }

  //Tof flag
  bool flag_tof_stop = false;
  {
    static const int device_id    = gUnpacker.get_device_id("TFlag");
    static const int data_type_id = gUnpacker.get_data_id("TFlag", "tdc");

    int mhit = gUnpacker.get_entries(device_id, 0, trigger::kTofTiming, 0, data_type_id);
    for(int m = 0; m<mhit; ++m){
      int tof_timing = gUnpacker.get(device_id, 0, trigger::kTofTiming, 0, data_type_id, m);
      if(!(MinStopTimingSdcOut < tof_timing && tof_timing < MaxStopTimingSdcOut)) flag_tof_stop = true;
    }// for(m)
  }

  HF1(1, 6.);

  //////SdcIn
#if 0
  static const double MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
  //////////////SdcIn number of hit layer
  DCAna->DecodeSdcInHits(rawData);
  double multi_SdcIn=0.;
  for(int layer=1; layer<=NumOfLayersSdcIn; ++layer){
    int nhIn = DCAna->GetSdcInHC(layer).size();
    multi_SdcIn += double(nhIn);
  }
# if MaxMultiCut
  if(multi_SdcIn/double(NumOfLayersSdcIn) > MaxMultiHitSdcIn)
   return true;
# endif
#endif


  HF1(1, 10.);
  double offset = flag_tof_stop ? 0 : StopTimeDiffSdcOut;
  DCAna->DecodeSdcOutHits(rawData, offset);
#if TotCut
  DCAna->TotCutSDC3(MinTotSDC3);
  DCAna->TotCutSDC4(MinTotSDC4);
#endif
  double multi_SdcOut = 0.;
  {
    for(int layer=1; layer<=NumOfLayersSdcOut; ++layer) {
      const DCHitContainer &contOut =DCAna->GetSdcOutHC(layer);
      int nhOut=contOut.size();
      event.nhit[layer-1] = nhOut;
      if(nhOut>0) event.nlayer++;
      multi_SdcOut += double(nhOut);
      HF1(100*layer, nhOut);
      int plane_eff = (layer-1)*3;
      bool fl_valid_sig = false;

      for(int i=0; i<nhOut; ++i){
	DCHit *hit=contOut[i];
	double wire=hit->GetWire();
	HF1(100*layer+1, wire-0.5);
	int nhtdc = hit->GetTdcSize();
	int tdc1st = -1;
	for(int k=0; k<nhtdc; k++){
	  int tdc = hit->GetTdcVal(k);
	  HF1(100*layer+2, tdc);
	  HF1(10000*layer+int(wire), tdc);
	  //	    HF2(1000*layer, tdc, wire-0.5);
	  if(tdc > tdc1st){
	    tdc1st = tdc;
	    fl_valid_sig = true;
	  }
	}
	HF1(100*layer+6, tdc1st);

	if(i<MaxHits)
	  event.pos[layer-1][i] = hit->GetWirePosition();

	int nhdt = hit->GetDriftTimeSize();
	int tot1st = -1;
	for(int k=0; k<nhdt; k++){
	  double dt = hit->GetDriftTime(k);
	  if(flag_tof_stop) HF1(100*layer+3, dt);
	  else              HF1(100*layer+8, dt);
	  HF1(10000*layer+1000+int(wire), dt);

	  double tot = hit->GetTot(k);
	  HF1(100*layer+5, tot);
	  if(tot > tot1st){
	    tot1st = tot;
	  }
	}
	HF1(100*layer+7, tot1st);
	int nhdl = hit->GetDriftTimeSize();
	for(int k=0; k<nhdl; k++){
	  double dl = hit->GetDriftLength(k);
	  HF1(100*layer+4, dl);
	}
      }
      if(fl_valid_sig) ++plane_eff;
      HF1(38, plane_eff);
    }
  }

#if MaxMultiCut
  if(multi_SdcOut/double(NumOfLayersSdcOut) > MaxMultiHitSdcOut)
    return true;
#endif

  HF1(1, 11.);

  // std::cout << "==========TrackSearch SdcOut============" << std::endl;
  if(flag_tof_stop){
#if UseTOF
    DCAna->TrackSearchSdcOut(TOFCont);
#else
    DCAna->TrackSearchSdcOut();
#endif
  }else{
    DCAna->TrackSearchSdcOut();
  }

#if 1
 #if Chi2Cut
  DCAna->ChiSqrCutSdcOut(30.);
 #endif

  int nt=DCAna->GetNtracksSdcOut();
  if(MaxHits<nt){
    std::cout << "#W " << FUNC_NAME << " "
	      << "too many ntSdcOut " << nt << "/" << MaxHits << std::endl;
    nt = MaxHits;
  }
  event.ntrack=nt;
  HF1(10, double(nt));
  for(int it=0; it<nt; ++it){
    DCLocalTrack *tp=DCAna->GetTrackSdcOut(it);
    if(!tp) continue;
    int nh=tp->GetNHit();
    double chisqr    = tp->GetChiSquare();
    double chisqr1st = tp->GetChiSquare1st();
    double x0=tp->GetX0(), y0=tp->GetY0();
    double u0=tp->GetU0(), v0=tp->GetV0();
    double theta = tp->GetTheta();
    int    nitr  = tp->GetNIteration();
    event.chisqr[it]=chisqr;
    event.x0[it]=x0;
    event.y0[it]=y0;
    event.u0[it]=u0;
    event.v0[it]=v0;

    HF1(11, double(nh));
    HF1(12, chisqr);
    HF1(14, x0); HF1(15, y0);
    HF1(16, u0); HF1(17, v0);
    HF2(18, x0, u0); HF2(19, y0, v0);
    HF2(20, x0, y0);

    HF1(28, nitr);
    HF1(29, chisqr1st);
    HF1(30, chisqr1st-chisqr);
    if(theta<=10.)
      HF1(31, chisqr1st-chisqr);
    if(10.<theta && theta<=20.)
      HF1(32, chisqr1st-chisqr);
    if(20.<theta && theta<=30.)
      HF1(33, chisqr1st-chisqr);
    if(30.<theta && theta<=40.)
      HF1(34, chisqr1st-chisqr);
    if(40.<theta)
      HF1(35, chisqr1st-chisqr);

    double xtof=tp->GetX(zTOF), ytof=tp->GetY(zTOF);
    double utof=u0, vtof=v0;
    HF1(21, xtof); HF1(22, ytof);
    HF1(23, utof); HF1(24, vtof);
    HF2(25, xtof, utof); HF2(26, ytof, vtof);
    HF2(27, xtof, ytof);

    for(int ih=0; ih<nh; ++ih){
      DCLTrackHit *hit=tp->GetHit(ih);
      if(!hit) continue;

      int layerId = 0;
      layerId = hit->GetLayer();
      if(layerId <= PlMaxSdcOut) layerId -= PlOffsSdcOut;
      else                         layerId -= PlOffsTOF - NumOfLayersSdcOut;


      HF1(13, hit->GetLayer());
      HF1(36, double(nh));
      HF1(37, chisqr);

      double wire=hit->GetWire();
      double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
      HF1(100*layerId+11, wire-0.5);
      HF1(100*layerId+12, dt);
      HF1(100*layerId+13, dl);
      double xcal=hit->GetXcal(), ycal=hit->GetYcal();
      double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
      HF1(100*layerId+14, pos);
      HF1(100*layerId+15, res);
      HF2(100*layerId+16, pos, res);
      HF2(100*layerId+17, xcal, ycal);
      //      HF1(100000*layerId+50000+wire, res);
      double wp=hit->GetWirePosition();
      double sign=1.;
      if(pos-wp<0.) sign=-1;
      HF2(100*layerId+18, sign*dl, res);
      double xlcal=hit->GetLocalCalPos();
      HF2(100*layerId+19, dt, xlcal-wp);

      HF2(100*layerId+31, xcal, res);
      HF2(100*layerId+32, ycal, res);
      HF2(100*layerId+33, u0, res);
      HF2(100*layerId+34, v0, res);

      double tot = hit->GetTot();
      HF1(100*layerId+40, tot);

      if (theta>=0 && theta<15)
	HF1(100*layerId+71, res);
      else if (theta>=15 && theta<30)
	HF1(100*layerId+72, res);
      else if (theta>=30 && theta<45)
	HF1(100*layerId+73, res);
      else if (theta>=45)
	HF1(100*layerId+74, res);

      if (std::abs(dl-std::abs(xlcal-wp))<2.0) {
	HFProf(100*layerId+20, dt, std::abs(xlcal-wp));
	HF2(100*layerId+22, dt, std::abs(xlcal-wp));
	HFProf(100000*layerId+3000+int(wire), xlcal-wp,dt);
	HF2(100000*layerId+4000+int(wire), xlcal-wp,dt);
      }
    }
  }
#endif

  HF1(1, 12.);

  return true;
}

//_____________________________________________________________________________
void
EventSdcOutTracking::InitializeEvent()
{
  event.evnum     =  0;
  event.trignhits =  0;
  event.nlayer    =  0;
  event.ntrack    =  0;
  event.nhBh2     =  0;
  event.nhBh1     =  0;
  event.nhTof     =  0;

  event.btof      = -999.;

  event.Time0Seg  = -1;
  event.deTime0   = -1;
  event.Time0     = -999;
  event.CTime0    = -999;

  for(int it=0; it<MaxHits; it++){
    event.tBh1[it]   = -9999.;
    event.deBh1[it]  = -9999.;

    event.Bh2Seg[it] = -1;
    event.tBh2[it]   = -9999.;
    event.deBh2[it]  = -9999.;

    event.stof[it]   = -9999.;

    event.TofSeg[it] = -1;
    event.tTof[it]   = -9999.;
    event.dtTof[it]  = -9999.;
    event.deTof[it]  = -9999.;

    event.chisqr[it] = -1.;
    event.x0[it] = -9999.;
    event.y0[it] = -9999.;
    event.u0[it] = -9999.;
    event.v0[it] = -9999.;

    event.trigpat[it]  = -1;
  }

  for(int it=0; it<NumOfSegTrig; it++){
    event.trigflag[it] = -1;
  }

  for(int it=0; it<NumOfLayersSdcOut; ++it){
    event.nhit[it] = -1;
    for(int that=0; that<MaxHits; ++that){
      event.pos[it][that] = -9999.;
    }
  }
}

//_____________________________________________________________________________
bool
EventSdcOutTracking::ProcessingEnd()
{
  tree->Fill();
  return true;
}

//_____________________________________________________________________________
VEvent*
ConfMan::EventAllocator()
{
  return new EventSdcOutTracking;
}

//_____________________________________________________________________________
namespace
{
  const int    NbinSdcOutTdc = 2000;
  const double MinSdcOutTdc  =    0.;
  const double MaxSdcOutTdc  = 2000.;

  const int    NbinSDC3DT = 240;
  const double MinSDC3DT  = -50.;
  const double MaxSDC3DT  = 150.;
  const int    NbinSDC3DL =  90;
  const double MinSDC3DL  =  -2.;
  const double MaxSDC3DL  =   7.;

  const int    NbinSDC4DT = 480;
  const double MinSDC4DT  = -50.;
  const double MaxSDC4DT  = 350.;
  const int    NbinSDC4DL = 180;
  const double MinSDC4DL  =  -3.;
  const double MaxSDC4DL  =  15.;
}
//_____________________________________________________________________________
bool
ConfMan::InitializeHistograms()
{
  HB1(1, "Status", 20, 0., 20.);

  for(int i=1; i<=NumOfLayersSdcOut+NumOfLayersTOF; ++i){

    std::string tag;
    int nwire = 0, nbindt = 1, nbindl = 1;
    double mindt = 0., maxdt = 1., mindl = 0., maxdl = 1.;
    if(i<=NumOfLayersSDC3){
      tag    = "SDC3";
      nwire  = MaxWireSDC3;
      nbindt = NbinSDC3DT;
      mindt  = MinSDC3DT;
      maxdt  = MaxSDC3DT;
      nbindl = NbinSDC3DL;
      mindl  = MinSDC3DL;
      maxdl  = MaxSDC3DL;
    }else if(i<=NumOfLayersSdcOut){
      tag = "SDC4";
      nwire   = (i==5 || i==6) ? MaxWireSDC4Y : MaxWireSDC4X;
      nbindt = NbinSDC4DT;
      mindt  = MinSDC4DT;
      maxdt  = MaxSDC4DT;
      nbindl = NbinSDC4DL;
      mindl  = MinSDC4DL;
      maxdl  = MaxSDC4DL;
    }else if(i<=NumOfLayersSdcOut+NumOfLayersTOF){
      tag = "TOF";
      nwire = NumOfSegTOF;
      nbindt = NbinSDC4DT;
      mindt  = MinSDC4DT;
      maxdt  = MaxSDC4DT;
      nbindl = NbinSDC4DL;
      mindl  = MinSDC4DL;
      maxdl  = MaxSDC4DL;
    }

    if(i<=NumOfLayersSdcOut){
      TString title0 = Form("#Hits %s#%2d", tag.c_str(), i);
      TString title1 = Form("Hitpat %s#%2d", tag.c_str(), i);
      TString title2 = Form("Tdc %s#%2d", tag.c_str(), i);
      TString title3 = Form("Drift Time %s#%2d", tag.c_str(), i);
      TString title4 = Form("Drift Length %s#%2d", tag.c_str(), i);
      TString title5 = Form("TOT %s#%2d", tag.c_str(), i);
      TString title6 = Form("Tdc 1st %s#%2d", tag.c_str(), i);
      TString title7 = Form("TOT 1st %s#%2d", tag.c_str(), i);
      TString title8 = Form("Drift Time %s#%2d (BH2 timing)", tag.c_str(), i);
      HB1(100*i+0, title0, nwire+1, 0., double(nwire+1));
      HB1(100*i+1, title1, nwire, 0., double(nwire));
      HB1(100*i+2, title2, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
      HB1(100*i+3, title3, nbindt, mindt, maxdt);
      HB1(100*i+4, title4, nbindl, mindl, maxdl);
      HB1(100*i+5, title5, 500,    0, 500);
      HB1(100*i+6, title6, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
      HB1(100*i+7, title7, 500,    0, 500);
      HB1(100*i+8, title8, nbindt, mindt, maxdt);
      for (int wire=1; wire<=nwire; wire++){
	TString title11 = Form("Tdc %s#%2d  Wire#%4d", tag.c_str(), i, wire);
	TString title12 = Form("DriftTime %s#%2d Wire#%4d", tag.c_str(), i, wire);
	TString title13 = Form("DriftLength %s#%2d Wire#%d", tag.c_str(), i, wire);
	TString title14 = Form("DriftTime %s#%2d Wire#%4d [Track]", tag.c_str(), i, wire);
	HB1(10000*i+wire, title11, NbinSdcOutTdc, MinSdcOutTdc, MaxSdcOutTdc);
	HB1(10000*i+1000+wire, title12, nbindt, mindt, maxdt);
	HB1(10000*i+2000+wire, title13, nbindl, mindl, maxdl);
	HB1(10000*i+5000+wire, title14, nbindt, mindt, maxdt);
      }
    }


    // Tracking Histgrams
    TString title11 = Form("HitPat SdcOut%2d [Track]", i);
    TString title12 = Form("DriftTime SdcOut%2d [Track]", i);
    TString title13 = Form("DriftLength SdcOut%2d [Track]", i);
    TString title14 = Form("Position SdcOut%2d", i);
    TString title15 = Form("Residual SdcOut%2d", i);
    TString title16 = Form("Resid%%Pos SdcOut%2d", i);
    TString title17 = Form("Y%%Xcal SdcOut%2d", i);
    TString title18 = Form("Res%%DL SdcOut%2d", i);
    TString title19 = Form("HitPos%%DriftTime SdcOut%2d", i);
    TString title20 = Form("DriftLength%%DriftTime SdcOut%2d", i);
    TString title21 = title15 + " [w/o Self]";
    TString title22 = title20;
    TString title40 = Form("TOT SdcOut%2d [Track]", i);
    TString title71 = Form("Residual SdcOut%2d (0<theta<15)", i);
    TString title72 = Form("Residual SdcOut%2d (15<theta<30)", i);
    TString title73 = Form("Residual SdcOut%2d (30<theta<45)", i);
    TString title74 = Form("Residual SdcOut%2d (45<theta)", i);
    HB1(100*i+11, title11, nwire, 0., nwire);
    HB1(100*i+12, title12, 600, -100, 400);
    HB1(100*i+13, title13, 100, -5, maxdl);
    HB1(100*i+14, title14, 100, -1000., 1000.);
    if(i<=NumOfLayersSdcOut)
      HB1(100*i+15, title15, 1000, -5.0, 5.0);
    else
      HB1(100*i+15, title15, 1000, -1000.0, 1000.0);
    if(i<=NumOfLayersSdcOut)
        HB2(100*i+16, title16, 400, -1000., 1000., 100, -1.0, 1.0);
    else
      HB2(100*i+16, title16, 100, -1000., 1000., 100, -1000.0, 1000.0);
    HB2(100*i+17, title17, 100, -1000., 1000., 100, -1000., 1000.);
    if(i<=NumOfLayersSDC3)
      HB2(100*i+18, title18, 110, -5.5, 5.5, 100, -1.0, 1.0);
    else
      HB2(100*i+18, title18, 110, -11., 11., 100, -1.0, 1.0);
    HB2(100*i+19, title19, 200, -50., maxdt, 200, -maxdl, maxdl);
    HBProf(100*i+20, title20, 100, -50, 300, 0, maxdl);
    HB2(100*i+22, title22, 200, -50, maxdt, 100, -0.5, maxdl);
    HB1(100*i+21, title21, 200, -5.0, 5.0);
    HB2(100*i+31, Form("Resid%%X SdcOut %d", i), 100, -1000., 1000., 100, -2., 2.);
    HB2(100*i+32, Form("Resid%%Y SdcOut %d", i), 100, -1000., 1000., 100, -2., 2.);
    HB2(100*i+33, Form("Resid%%U SdcOut %d", i), 100, -0.5, 0.5, 100, -2., 2.);
    HB2(100*i+34, Form("Resid%%V SdcOut %d", i), 100, -0.5, 0.5, 100, -2., 2.);
    HB1(100*i+40, title40, 360,    0, 300);
    HB1(100*i+71, title71, 200, -5.0, 5.0);
    HB1(100*i+72, title72, 200, -5.0, 5.0);
    HB1(100*i+73, title73, 200, -5.0, 5.0);
    HB1(100*i+74, title74, 200, -5.0, 5.0);

    for (int j=1; j<=nwire; j++) {
      TString title = Form("XT of Layer %2d Wire #%4d", i, j);
      HBProf(100000*i+3000+j, title, 100, -12., 12., -30, 300);
      HB2(100000*i+4000+j, title, 100, -12., 12., 100, -30., 300.);
    }

  }

  // Tracking Histgrams
  HB1(10, "#Tracks SdcOut", 10, 0., 10.);
  HB1(11, "#Hits of Track SdcOut", 20, 0., 20.);
  HB1(12, "Chisqr SdcOut", 500, 0., 50.);
  HB1(13, "LayerId SdcOut", 60, 30., 90.);
  HB1(14, "X0 SdcOut", 1400, -1200., 1200.);
  HB1(15, "Y0 SdcOut", 1000, -500., 500.);
  HB1(16, "U0 SdcOut", 200, -0.35, 0.35);
  HB1(17, "V0 SdcOut", 200, -0.20, 0.20);
  HB2(18, "U0%X0 SdcOut", 120, -1200., 1200., 100, -0.40, 0.40);
  HB2(19, "V0%Y0 SdcOut", 100, -500., 500., 100, -0.20, 0.20);
  HB2(20, "X0%Y0 SdcOut", 100, -1200., 1200., 100, -500, 500);
  HB1(21, "Xtof SdcOut", 1400, -1200., 1200.);
  HB1(22, "Ytof SdcOut", 1000, -500., 500.);
  HB1(23, "Utof SdcOut", 200, -0.35, 0.35);
  HB1(24, "Vtof SdcOut", 200, -0.20, 0.20);
  HB2(25, "Utof%Xtof SdcOut", 100, -1200., 1200., 100, -0.40, 0.40);
  HB2(26, "Vtof%Ytof SdcOut", 100, -500., 500., 100, -0.20, 0.20);
  HB2(27, "Xtof%Ytof SdcOut", 100, -1200., 1200., 100, -500, 500);

  HB1(28, "NIteration SdcOut", 100, 0., 100.);
  HB1(29, "Chisqr1st SdcOut", 500, 0., 50.);
  HB1(30, "Chisqr1st-Chisqr SdcOut", 500, 0., 10.);
  HB1(31, "Chisqr1st-Chisqr SdcOut (0<theta<10)", 500, 0., 10.);
  HB1(32, "Chisqr1st-Chisqr SdcOut (10<theta<20)", 500, 0., 10.);
  HB1(33, "Chisqr1st-Chisqr SdcOut (20<theta<30)", 500, 0., 10.);
  HB1(34, "Chisqr1st-Chisqr SdcOut (30<theta<40)", 500, 0., 10.);
  HB1(35, "Chisqr1st-Chisqr SdcOut (40<theta)", 500, 0., 10.);
  HB1(36, "#Hits of Track SdcOut(SDC)", 20, 0., 20.);
  HB1(37, "Chisqr SdcOut(SDC)", 500, 0., 50.);

  // Plane eff
  HB1(38, "Plane Eff", 30, 0, 30);

  ////////////////////////////////////////////
  //Tree
  HBTree("sdcout","tree of SdcOutTracking");
  tree->Branch("evnum",     &event.evnum,     "evnum/I");
  tree->Branch("trignhits", &event.trignhits, "trignhits/I");
  tree->Branch("trigpat",    event.trigpat,   Form("trigpat[%d]/I", MaxHits));
  tree->Branch("trigflag",   event.trigflag,  Form("trigflag[%d]/I", NumOfSegTrig));

  tree->Branch("nhBh1",    &event.nhBh1,   "nhBh1/I");
  tree->Branch("tBh1",      event.tBh1,    Form("tBh1[%d]/D",   MaxHits));
  tree->Branch("deBh1",     event.deBh1,   Form("deBh1[%d]/D",  MaxHits));

  tree->Branch("nhBh2",    &event.nhBh2,   "nhBh2/I");
  tree->Branch("tBh2",      event.tBh2,    Form("tBh2[%d]/D",   MaxHits));
  tree->Branch("deBh2",     event.deBh2,   Form("deBh2[%d]/D",  MaxHits));
  tree->Branch("Bh2Seg",    event.Bh2Seg,  Form("Bh2Seg[%d]/D", MaxHits));

  tree->Branch("Time0Seg", &event.Time0Seg,  "Time0Seg/D");
  tree->Branch("deTime0",  &event.deTime0,   "deTime0/D");
  tree->Branch("Time0",    &event.Time0,     "Time0/D");
  tree->Branch("CTime0",   &event.CTime0,    "CTime0/D");

  tree->Branch("btof",     &event.btof,     "btof/D");
  tree->Branch("stof",     event.stof,    Form("stof[%d]/D", MaxHits));

  tree->Branch("nhTof",   &event.nhTof,   "nhTof/I");
  tree->Branch("TofSeg",   event.TofSeg,  "TofSeg[nhTof]/D");
  tree->Branch("tTof",     event.tTof,    "tTof[nhTof]/D");
  tree->Branch("dtTof",    event.dtTof,   "dtTof[nhTof]/D");
  tree->Branch("deTof",    event.deTof,   "deTof[nhTof]/D");

  tree->Branch("nhit",     &event.nhit,     Form("nhit[%d]/I", NumOfLayersSdcOut));
  tree->Branch("nlayer",   &event.nlayer,   "nlayer/I");
  tree->Branch("pos",      &event.pos,     Form("pos[%d][%d]/D", NumOfLayersSdcOut, MaxHits));
  tree->Branch("ntrack",   &event.ntrack,   "ntrack/I");
  tree->Branch("chisqr",    event.chisqr,   "chisqr[ntrack]/D");
  tree->Branch("x0",        event.x0,       "x0[ntrack]/D");
  tree->Branch("y0",        event.y0,       "y0[ntrack]/D");
  tree->Branch("u0",        event.u0,       "u0[ntrack]/D");
  tree->Branch("v0",        event.v0,       "v0[ntrack]/D");

  HPrint();
  return true;
}

//_____________________________________________________________________________
bool
ConfMan::InitializeParameterFiles()
{
  return
    (InitializeParameter<DCGeomMan>("DCGEO")        &&
      InitializeParameter<DCDriftParamMan>("DCDRFT") &&
      InitializeParameter<DCTdcCalibMan>("DCTDC")    &&
      InitializeParameter<HodoParamMan>("HDPRM")     &&
      InitializeParameter<HodoPHCMan>("HDPHC")       &&
      InitializeParameter<UserParamMan>("USER"));
}

//_____________________________________________________________________________
bool
ConfMan::FinalizeProcess()
{
  return true;
}
