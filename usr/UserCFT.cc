/**
 *  file: UserCFT.cc
 *  date: 2018.07.20
 *  based on UserBFT.cc
 */

#include <cmath>
#include <iostream>
#include <sstream>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "RMAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "KuramaLib.hh"
#include "RawData.hh"
#include "RootHelper.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#include "Kinematics.hh"

#define HodoCut 0 // with BH1/BH2
#define TimeCut 0 // in cluster analysis

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

namespace
{
  using namespace root;
  const std::string& class_name("EventCFT");
  RMAnalyzer&         gRM   = RMAnalyzer::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();
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
class EventCFT : public VEvent
{
private:
  RawData      *rawData;
  DCAnalyzer   *DCAna;
  HodoAnalyzer *hodoAna;

public:
        EventCFT( void );
       ~EventCFT( void );
  bool  ProcessingBegin( void );
  bool  ProcessingEnd( void );
  bool  ProcessingNormal( void );
  bool  InitializeHistograms( void );
  void  InitializeEvent( void );
};

//______________________________________________________________________________
EventCFT::EventCFT( void )
  : VEvent(),
    rawData(0),
    DCAna( new DCAnalyzer ),
    hodoAna( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventCFT::~EventCFT( void )
{
  if ( hodoAna ){
    delete hodoAna;
    hodoAna = NULL;
  }
  if ( DCAna ){
    delete DCAna;
    DCAna   = NULL;
  }
  if ( rawData ){
    delete rawData;
    rawData = NULL;
  }
}

//______________________________________________________________________________
struct Event
{
  int evnum;

  int trigpat[NumOfSegTrig];
  int trigflag[NumOfSegTrig];

  // Fiber Hit
  
  int    nhits;
  int    unhits;
  int    dnhits;

  // CFT
  int ntCFT;
  int ncl;
  double phi[MaxDepth];
  double theta[MaxDepth];

  double dphi[NumOfPlaneCFT][MaxDepth];
  double phi_ini[NumOfPlaneCFT][MaxDepth];
  double phi_track[NumOfPlaneCFT][MaxDepth];

  double dz[NumOfPlaneCFT][MaxDepth];
  double z_ini[NumOfPlaneCFT][MaxDepth];
  double z_track[NumOfPlaneCFT][MaxDepth];
  ThreeVector Pos[MaxDepth];
  ThreeVector Dir[MaxDepth];
  double vtx_x[MaxDepth], vtx_y[MaxDepth], vtx_z[MaxDepth];
  double vtxAB_x, vtxAB_y, vtxAB_z;

  // for cosmic ray tracking
  double dphi16[NumOfPlaneCFT][2];
  double phi16_ini[NumOfPlaneCFT][2];
  double phi16_track[NumOfPlaneCFT][2];

  double dz16[NumOfPlaneCFT][2];
  double z16_ini[NumOfPlaneCFT][2];
  double z16_track[NumOfPlaneCFT][2];

};

//______________________________________________________________________________
namespace root
{
  Event event;
  TH1   *h[MaxHist];
  TTree *tree;
}

//______________________________________________________________________________
bool
EventCFT::ProcessingBegin( void )
{
  InitializeEvent();
  return true;
}

//______________________________________________________________________________
bool
EventCFT::ProcessingNormal( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"]");


  rawData = new RawData;
  rawData->DecodeHits();

  gRM.Decode();

  int evnum = gRM.EventNumber();
  event.evnum = evnum;

  //**************************************************************************
  //******************RawData

  // CFT   
  {
    for(int p = 0; p<NumOfPlaneCFT; ++p){
      int sumnhit = 0;
      int layer = p;

      double Maxpe = -10.;
      double Maxpeup = -10.;
      double Maxpedown = -10.;      
      int Maxpe_seg = -10;
      int Maxpe_segup = -10;
      int Maxpe_segdown = -10;
      double NofP = 0.;
      double newPeak = 0.;
          
      const HodoRHitContainer &cont = rawData->GetCFTRawHC(p);
      int nhit = cont.size();
      HF1 (1000*(p+1)+1, nhit);
      //event.nHitCFT[p] = nhit;
      
      for(int i = 0; i<nhit; ++i){	
	HodoRawHit *hit = cont.at(i);
	int seg = hit->SegmentId(); 
	int NhitT = hit->GetSizeTdcUp();
	int NhitT_tr = hit->GetSizeTdcDown();
	int NhitAH = hit->GetSizeAdcUp();	  
	int NhitAL = hit->GetSizeAdcDown();
	
	//TDC
	bool fl_tdc_ok = false;
	
	bool thitflag = false;
	bool cosmic_tflag = false;
	bool cosmic_wflag = false;
	bool cosmic_tflag2 = false;
	bool cosmic_wflag2 = false;
	
	for(int m = 0; m<NhitT; ++m){	    
	  int bufT = hit->GetTdcUp(m);
	  HF2 (1000*(p+1)+100, seg, bufT);  //TDC Nhits 
	  
	  int width = -999, bufTtr = -999;	      
	  int flag_w = -1;
	  if(m<NhitT_tr){
	    bufTtr = hit->GetTdcDown(m);	      
	    width = bufT - bufTtr;
	    HF2 (1000*(p+1)+101, seg, bufTtr);
	    HF2 (1000*(p+1)+104, seg, width);
	    if(width>30){
	      HF2 (1000*(p+1) +122, seg, bufT);
	    }
	  }
	}
	
	//ADC Hi
	for(int m = 0; m<NhitAH; ++m){
	  int bufAH = hit->GetAdcUp();	
	  HF2 (1000*(p+1)+200, seg, bufAH);
	}
	
	//ADC Low
	for(int m = 0; m<NhitAL; ++m){	    
	  int bufAL = hit->GetAdcDown();	
	  HF2 (1000*(p+1)+201, seg, bufAL);
	}
      }
    }
  }
  

  // FiberHitCFT    
  hodoAna->DecodeCFTHits(rawData);
  for(int p = 0; p<NumOfPlaneCFT; ++p){

    int nhit = hodoAna->GetNHitsCFT(p);          
    int nhit_t = 0;
    for(int i = 0; i<nhit; ++i){
      const FiberHit* hit = hodoAna->GetHitCFT(p, i);
      int mhit = hit->GetNumOfHit();
      int seg_id = hit->PairId();
      HF1(1000*(p+1) +20, seg_id);

      bool fl_m = false;
      for(int m = 0; m<mhit; ++m){

	double leading = hit->GetLeading(m);
	
	double ctime = hit->GetCTime(m);
	double time = hit->GetTime(m);
	double width = hit->GetWidth(m);	
	double adcHi = hit->GetAdcHi();
	double adcLow = hit->GetAdcLow();  

	double CFT_r       = hit->GetPositionR();
	double CFT_phi   = hit->GetPositionPhi();	  

	HF2 (1000*(p+1) +102, seg_id, leading);
	HF2 (1000*(p+1) +103, seg_id, ctime);
	HF2 (1000*(p+1) +202, seg_id, adcHi);
	HF2 (1000*(p+1) +203, seg_id, adcLow);
	
	if(-30 < ctime && ctime < 30){
	  if(fl_m==false){
	    fl_m=true;
	    nhit_t++;
	    HF1(1000*(p+1)  +21, seg_id);
	    HF2( 3,CFT_r*cos(CFT_phi*Deg2Rad),CFT_r*sin(CFT_phi*Deg2Rad));
	    HF2(1000*(p+1) +222, seg_id, adcHi);
	    HF2(1000*(p+1) +233, seg_id, adcLow);
	  }
	}

      }// mhit      
    }//nhit
    HF1 (1000*(p+1)+11, nhit_t);      
  }
  

#if 1
  // Fiber Cluster
  for(int p = 0; p<NumOfPlaneCFT; ++p){
    hodoAna->TimeCutCFT(p, -30, 30); // CATCH@J-PARC  

    int ncl = hodoAna->GetNClustersCFT(p);
    HF1 (1000*(p+1)+12, ncl);      
    
    for(int i=0; i<ncl; ++i){
      FiberCluster *cl = hodoAna->GetClusterCFT(p,i);
      if(!cl) continue;
      double size  = cl->ClusterSize();
      double seg = cl->MeanSeg();
      double ctime = cl->CMeanTime();
      double r   = cl->MeanPositionR();
      double phi   = cl->MeanPositionPhi();

      HF2(4, r*cos(phi*Deg2Rad),r*sin(phi*Deg2Rad));      
    }
    
  }
#endif


  // CFT tracking
#if 1
  DCAna->DecodeCFTHits( rawData );

  DCAna->TrackSearchCFT();

  int ntCFT=DCAna->GetNtracksCFT();// vtx limit ver.
  event.ntCFT = ntCFT;
  DCLocalTrack *tpp[2];
  for( int i=0; i<ntCFT; ++i ){

    DCLocalTrack *tp=DCAna->GetTrackCFT(i);
    if(i==0){tpp[0]=tp;}
    else if(i==1){tpp[1]=tp;}

    int nh   = tp->GetNHit();
    int nhUV = tp->GetNHitUV();
    double chisqrXY=tp->GetChiSquareXY();
    double chisqrXYZ=tp->GetChiSquareZ();
    double vtx_z =tp->GetVtxZ();
    double theta =tp->GetThetaCFT();
    int xyFlag = tp->GetCFTxyFlag();
    int zFlag  = tp->GetCFTzFlag() ;
       
    ThreeVector Pos0 = tp->GetPos0();
    ThreeVector Dir = tp->GetDir();
    double A=(Dir.x()*Dir.x()+Dir.y()*Dir.y());

    // aka 
    double D=(Dir.x()*Dir.x()+Dir.y()*Dir.y()+Dir.z()*Dir.z());

    double phi = -999.;
    if(Dir.x()>=0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*Rad2Deg;
    }//0~90
    else if (Dir.x()<0 && Dir.y()>=0){
      phi = acos(Dir.x()/sqrt(A))*Rad2Deg;
    }//90~180
    else if (Dir.x()<0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*Rad2Deg;
    }//180~270
    else if (Dir.x()>=0 && Dir.y()<0){
      phi = 360. - acos(Dir.x()/sqrt(A))*Rad2Deg; ;
    }//270~360
    else{
    }


    HF1(5, vtx_z); 
    HF1(6, theta);
    event.theta[i] = theta;
    event.phi[i] = phi;

    event.Pos[i] = Pos0;
    event.Dir[i] = Dir;

    // vertex
    ThreeVector  bPos(0., 0., 0.);
    ThreeVector  bdir(0., 0., 1.);
    double dist = 1000.;
    if(Pos0.x()>-500.){
      ThreeVector vtx = Kinematics::VertexPoint3D(bPos, Pos0, bdir, Dir, dist);
      event.vtx_x[i] = vtx.x();
      event.vtx_y[i] = vtx.y();
      event.vtx_z[i] = vtx.z();
    }

    // straight layer
    for(int ip=0; ip<nh; ip++){
      DCLTrackHit *hit = tp->GetHit(ip);
      int layer = hit->GetLayer();
      int seg = (int)hit->GetMeanSeg();

      double phi_ini   = tp->GetPhiIni(layer);      
      double phi_track   = tp->GetPhiTrack(layer);      
      double z_track = tp->GetZTrack(layer);
      double dphi  = tp->GetdPhi(layer);

      event.phi_ini[layer][i] = phi_ini;
      event.phi_track[layer][i] = phi_track;
      event.dphi[layer][i] = dphi;
      event.z_track[layer][i] = z_track;
      HF2(1000*(layer+1)+300, z_track, dphi);
    }

    // spiral layer
    for(int ip=0; ip<nhUV; ip++){
      DCLTrackHit *hit = tp->GetHitUV(ip);
      int layer = hit->GetLayer();
      int seg = (int)hit->GetMeanSeg();
      
      double phi_track   = tp->GetPhiTrack(layer);      
      double z_track = tp->GetZTrack(layer);
      double z_ini   = tp->GetZIni(layer);      
      double dz    = tp->GetdZ(layer);
      event.phi_track[layer][i] = phi_track;
      event.z_ini[layer][i] = z_ini;
      event.z_track[layer][i] = z_track;
      event.dz[layer][i] = dz;
      HF2(1000*(layer+1)+310, phi_track, dz);
    }    
  }
 
  if(ntCFT>1){
    // vertex of 2 tracks
    ThreeVector PosA = event.Pos[0];
    ThreeVector PosB = event.Pos[1];
    ThreeVector DirA = event.Dir[0];
    ThreeVector DirB = event.Dir[1];
    double distAB = -999.;

    ThreeVector vtxAB = Kinematics::VertexPoint3D(PosA, PosB, DirA, DirB, distAB);
    event.vtxAB_x = vtxAB.x();
    event.vtxAB_y = vtxAB.y();
    event.vtxAB_z = vtxAB.z();

#endif

#if 0    
    // conbine to 16 layers tracking    
    double d_phi = fabs(event.phi[0]-event.phi[1]);
    if(d_phi>0){
      DCAna->DecodeCFT16Hits(rawData,tpp[0],0); // 1st track
      DCAna->DecodeCFT16Hits(rawData,tpp[1],1); // 2nd track
      
      DCAna->TrackSearchCFT16();
      int ntCFT16=DCAna->GetNtracksCFT16();
      for( int it=0; it<ntCFT16; ++it ){

	DCLocalTrack *tp=DCAna->GetTrackCFT16(it);
	int nh   = tp->GetNHit();
	int nhUV = tp->GetNHitUV();
	
	// straight layer
	int i=0;
	for(int ip=0; ip<nh; ip++){
	  DCLTrackHit *hit = tp->GetHit(ip);
	  int layer = hit->GetLayer();
	  int seg = (int)hit->GetMeanSeg();	  
	  double phi_ini   = tp->GetPhiIni(layer);      
	  double phi_track = tp->GetPhiTrack(layer);      
	  double z_track   = tp->GetZTrack(layer);
	  double dphi      = tp->GetdPhi(layer);
	  
	  if(layer<8){i=0;}
	  else if(layer>=8){layer-=8;i=1;}

	  event.phi16_ini[layer][i]   = phi_ini;
	  event.phi16_track[layer][i] = phi_track;
	  event.dphi16[layer][i]      = dphi;
	  event.z16_track[layer][i]   = z_track;

	}	
	// spiral layer
	for(int ip=0; ip<nhUV; ip++){
	  DCLTrackHit *hit = tp->GetHitUV(ip);
	  int layer = hit->GetLayer();
	  int seg = (int)hit->GetMeanSeg();	  
	  double phi_track = tp->GetPhiTrack(layer);      
	  double z_track   = tp->GetZTrack(layer);
	  double z_ini     = tp->GetZIni(layer);      
	  double dz        = tp->GetdZ(layer);

	  if(layer<8){i=0;}
	  else if(layer>=8){layer-=8;i=1;}

	  event.phi16_track[layer][i] = phi_track;
	  event.z16_ini[layer][i] = z_ini;
	  event.z16_track[layer][i] = z_track;
	  event.dz16[layer][i] = dz;

	}	     
      }
      
    }


#endif
  }



  // Trigger Flag
  {
    const HodoRHitContainer &cont=rawData->GetTrigRawHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int seg = hit->SegmentId()+1;
      int tdc = hit->GetTdc1();
      if( tdc ){
	event.trigpat[i]      = seg;
	event.trigflag[seg-1] = tdc;
      }
    }
  }


  return true;
}

//______________________________________________________________________________
bool
EventCFT::ProcessingEnd( void )
{
  tree->Fill();
  return true;
}

//______________________________________________________________________________
void
EventCFT::InitializeEvent( void )
{
  event.evnum  = 0;
  event.nhits  = 0;
  event.unhits = 0;
  event.dnhits = 0;
  event.ncl    = 0;

  for( int it=0; it<NumOfSegTrig; ++it ){
    event.trigpat[it]  = -1;
    event.trigflag[it] = -1;
  }


  event.ntCFT  = -1;  
  for(int i = 0; i<MaxDepth; ++i){
    event.phi[i]  = -999.;
    event.theta[i]  = -999.;
    event.vtx_x[i]  = -999.;
    event.vtx_y[i]  = -999.;
    event.vtx_z[i]  = -999.;
    for(int j = 0; j<3; ++j){
      event.Pos[i][j]  = -999.;
      event.Dir[i][j]  = -999.;
    }
    
    for( int p=0; p<NumOfPlaneCFT; ++p ){      
      event.dphi[p][i]       = -999.;
      event.phi_ini[p][i]    = -999.;
      event.phi_track[p][i]  = -999.;
      event.dz[p][i]       = -999.;
      event.z_ini[p][i]    = -999.;
      event.z_track[p][i]  = -999.;
    }
    
  }
  event.vtxAB_x = -999.;
  event.vtxAB_y = -999.; 
  event.vtxAB_z = -999.;

  for( int p=0; p<NumOfPlaneCFT; ++p ){      
    for(int i = 0; i<2; ++i){
      event.dphi16[p][i]       = -999.;
      event.phi16_ini[p][i]    = -999.;
      event.phi16_track[p][i]  = -999.;
      event.dz16[p][i]       = -999.;
      event.z16_ini[p][i]    = -999.;
      event.z16_track[p][i]  = -999.;
    }
  }


}

//______________________________________________________________________________
VEvent*
ConfMan::EventAllocator( void )
{
  return new EventCFT;
}

//______________________________________________________________________________
namespace
{
  const int    NbinTdc = 1000;
  const double MinTdc  =    0.;
  const double MaxTdc  = 1000.;

  const int    NbinTot =  136;
  const double MinTot  =   -8.;
  const double MaxTot  =  128.;

  const int    NbinTime = 1000;
  const double MinTime  = -500.;
  const double MaxTime  =  500.;
}
//______________________________________________________________________________
bool
ConfMan:: InitializeHistograms( void )
{
  HB1( 1, "Status", 20, 0., 20. );
  HB1( 2, "analys flag", 20, 0., 20. );

  //CFT
  HB2( 3, "Fiber Position", 400,-100,100,400,-100,100 );
  HB2( 4, "Fiber Mean Position", 400,-100,100,400,-100,100 );
  HB1( 5, "vertex z (CFT tracking)", 1000,-500,500 );
  HB1( 6, "theta (CFT tracking)", 180, 0 ,180 );
  for(int i=0; i<NumOfPlaneCFT; i++){
    std::ostringstream title[20];
   
    if(i%2 == 0){// spiral layer
      int layer = (int)i/2 +1;
      title[0] << "CFT UV%d  N hit" << layer;
      title[1] << "CFT UV%d  N hit w/ Tdc cut"     << layer;
      title[2] << "CFT UV%d  NCluster hit"         << layer;
      title[3] << "CFT UV%d  Hit Pattern"          << layer;
      title[4] << "CFT UV%d  Hit Pattern w/ Tdc cut" << layer;
      // TDC
      title[5] << "CFT UV%d  Tdc(Leading) vs seg"  << layer;     
      title[6] << "CFT UV%d  Tdc(Trailing) vs seg" << layer;
      // TDC Fiber Hit
      title[7] << " CFT UV%d Time(Leading) vs seg" << layer;  
      title[8] << " CFT UV%d Time(Leading) vs seg" << layer;  
      title[9] << " CFT UV%d CTime vs seg"         << layer;  
      title[10]<< " CFT UV%d width vs seg"         << layer; 
      // ADC
      title[11]<< "CFT UV%d Adc(High) vs seg"      << layer;
      title[12]<< "CFT UV%d Adc(Low) vs seg"       << layer;     
      // ADC Fiber Hit
      title[13] << "CFT UV%d Adc(High)-pedestal vs seg"            << layer;  
      title[14] << "CFT UV%d Adc(Low)-pedestal vs seg"             << layer;  
      title[15] << "CFT UV%d Adc(High)-pedestal vs seg w/ Tdc cut" << layer;  
      title[16] << "CFT UV%d Adc(Low)-pedestal vs seg w/ Tdc cut"  << layer; 
      // tracking
      title[17]<< "CFT UV%d  dphi vs z"          << layer;
      title[18]<< "CFT UV%d  dz vs phi"          << layer;     
    }else if(i%2 == 1){// straight layer
      int layer = (int)i/2 +1;
      title[0] << "CFT Phi%d  N hit" << layer;
      title[1] << "CFT Phi%d  N hit w/ Tdc cut" << layer;
      title[2] << "CFT Phi%d  NCluster hit" << layer;
      title[3] << "CFT Phi%d  Hit Pattern" << layer;
      title[4] << "CFT Phi%d  Hit Pattern w/ Tdc cut" << layer;
      // TDC
      title[5] << "CFT Phi%d  Tdc(Leading) vs seg" << layer;     
      title[6] << "CFT Phi%d  Tdc(Trailing) vs seg" << layer;
      // TDC Fiber Hit
      title[7] << " CFT Phi%d Time(Leading) vs seg" << layer;  
      title[8] << " CFT Phi%d Time(Leading) vs seg" << layer;  
      title[9] << " CFT Phi%d CTime vs seg"         << layer;  
      title[10]<< " CFT Phi%d width vs seg"         << layer; 
      // ADC
      title[11]<< "CFT Phi%d Adc(High) vs seg"      << layer;
      title[12]<< "CFT Phi%d Adc(Low) vs seg"       << layer;     
      // ADC Fiber Hit
      title[13] << "CFT Phi%d Adc(High)-pedestal vs seg"            << layer;  
      title[14] << "CFT Phi%d Adc(Low)-pedestal vs seg"             << layer;  
      title[15] << "CFT Phi%d Adc(High)-pedestal vs seg w/ Tdc cut" << layer;  
      title[16] << "CFT Phi%d Adc(Low)-pedestal vs seg w/ Tdc cut"  << layer; 
      // tracking
      title[17]<< "CFT Phi%d  dphi vs z"          << layer;
      title[18]<< "CFT Phi%d  dz vs phi"          << layer;     
    }

    HB1( 1000*(i+1)+1,  title[0].str().c_str(),50,0,50);
    HB1( 1000*(i+1)+11, title[1].str().c_str(),50,0,50);
    HB1( 1000*(i+1)+12, title[2].str().c_str(),50,0,50);
    HB1( 1000*(i+1)+20, title[3].str().c_str(), NumOfSegCFT[i]+20, -10, NumOfSegCFT[i]+10);
    HB1( 1000*(i+1)+21, title[4].str().c_str(), NumOfSegCFT[i]+20, -10, NumOfSegCFT[i]+10);

    //TDC
    HB2( 1000*(i+1)+100, title[5].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);
    HB2( 1000*(i+1)+101, title[6].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 1024,0,1024);    
    // Fiber Hit
    HB2( 1000*(i+1)+102, title[7].str().c_str(),  NumOfSegCFT[i], 0, NumOfSegCFT[i],1000,0,1000);
    HB2( 1000*(i+1)+122, title[8].str().c_str(),  NumOfSegCFT[i], 0, NumOfSegCFT[i],1000,0,1000);
    HB2( 1000*(i+1)+103, title[9].str().c_str(),  NumOfSegCFT[i], 0, NumOfSegCFT[i],1000,-500,500);
    HB2( 1000*(i+1)+104, title[10].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i],100,0,100);

    //ADC
    HB2( 1000*(i+1)+200, title[11].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4096,0,4096);
    HB2( 1000*(i+1)+201, title[12].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4096,0,4096);
    
    // Fiber Hit
    HB2( 1000*(i+1)+202, title[13].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4010,-10,4000);
    HB2( 1000*(i+1)+203, title[14].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4010,-10,4000);
    HB2( 1000*(i+1)+222, title[15].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4010,-10,4000);
    HB2( 1000*(i+1)+233, title[16].str().c_str(), NumOfSegCFT[i], 0, NumOfSegCFT[i], 4010,-10,4000);
    
    // tracking
    HB2( 1000*(i+1)+300, title[17].str().c_str(), 200,0,400,200,-5,5);
    HB2( 1000*(i+1)+310, title[18].str().c_str(), 180,0,360,200,-10,10);
    
  }

  HB1( 10000, "nhvalid", 10,0,10);

  //Tree
  HBTree( "tree","tree of Counter" );
  tree->Branch("evnum",    &event.evnum,    "evnum/I");
  tree->Branch("trigpat",   event.trigpat,  Form("trigpat[%d]/I", NumOfSegTrig));
  tree->Branch("trigflag",  event.trigflag, Form("trigflag[%d]/I", NumOfSegTrig));
  tree->Branch("nhits",     &event.nhits,        "nhits/I");
  
  tree->Branch("ntCFT",     &event.ntCFT,    "ntCFT/I");
  tree->Branch("theta",     event.theta,    "theta[ntCFT]/D");
  tree->Branch("phi",       event.phi,      "phi[ntCFT]/D");
  tree->Branch("vtx_x",     event.vtx_x,    "vtx_x[ntCFT]/D");
  tree->Branch("vtx_y",     event.vtx_y,    "vtx_y[ntCFT]/D");
  tree->Branch("vtx_z",     event.vtx_z,    "vtx_z[ntCFT]/D");
  tree->Branch("vtxAB_x",   &event.vtxAB_x,  "vtxAB_x/D");
  tree->Branch("vtxAB_y",   &event.vtxAB_y,  "vtxAB_y/D");
  tree->Branch("vtxAB_z",   &event.vtxAB_z,  "vtxAB_z/D");
  
  tree->Branch("dphi",     event.dphi,     Form("dphi[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("phi_ini",  event.phi_ini,  Form("phi_ini[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("phi_track",event.phi_track,Form("phi_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );

  tree->Branch("dz",      event.dz,      Form("dz[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("z_ini",   event.z_ini,   Form("z_ini[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );
  tree->Branch("z_track", event.z_track, Form("z_track[%d][%d]/D", NumOfPlaneCFT, MaxDepth ) );


  tree->Branch("dphi16",     event.dphi16,     Form("dphi16[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("phi16_ini",  event.phi16_ini,  Form("phi16_ini[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("phi16_track",event.phi16_track,Form("phi16_track[%d][%d]/D", NumOfPlaneCFT, 2 ) );

  tree->Branch("dz16",      event.dz16,      Form("dz16[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("z16_ini",   event.z16_ini,   Form("z16_ini[%d][%d]/D", NumOfPlaneCFT, 2 ) );
  tree->Branch("z16_track", event.z16_track, Form("z16_track[%d][%d]/D", NumOfPlaneCFT, 2 ) );


  HPrint();
  return true;
}

//______________________________________________________________________________
bool
ConfMan::InitializeParameterFiles( void )
{
  return
    ( InitializeParameter<DCGeomMan>("DCGEO")    &&
      InitializeParameter<HodoParamMan>("HDPRM") &&
      InitializeParameter<HodoPHCMan>("HDPHC")   &&
      InitializeParameter<UserParamMan>("USER")  );
}

//______________________________________________________________________________
bool
ConfMan::FinalizeProcess( void )
{
  return true;
}
