/**
 *  file: RawData.cc
 *  date: 2017.04.10
 *
 */

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "MTDCRawHit.hh"
#include "HodoRawHit.hh"
#include "SDDRawHit.hh"
#include "DCRawHit.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

#define OscillationCut 0
#define MsTPhase1      1

namespace
{
  using namespace hddaq::unpacker;
  const std::string& class_name("RawData");
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
  //  const UserParamMan&    gUser     = UserParamMan::GetInstance();
  enum EUorD { kOneSide=1, kBothSide=2 };
  enum EDCDataType { kLeading, kTrailing, kNDCDataType };
#if OscillationCut
  const int  MaxMultiHitDC  = 16;
#endif
  
  //______________________________________________________________________________
  inline bool
  AddHodoRawHit( HodoRHitContainer& cont,
		 int id, int plane, int seg, int UorD, int AorT, int data, bool DEBUG=false )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");
    
    HodoRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      HodoRawHit *q = cont[i];
      if( q->DetectorId()==id &&
	  q->PlaneId()==plane &&
	  q->SegmentId()==seg ){
	if(DEBUG) std::cout<<"add data"<<std::endl;
	p=q; break;
      }
    }
    if( !p ){
      p = new HodoRawHit( id, plane, seg );
      cont.push_back(p);
      if(DEBUG) std::cout<<"========= new hodorawhit ======="<<std::endl;
    }
    if(DEBUG){
      std::cout   << "DetectorId = " << id    << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "SegId      = " << seg   << std::endl
		  << "AorT       = " << AorT  << std::endl
		  << "UorD       = " << UorD  << std::endl;      
    }
    switch(AorT){
    case 0:
      if( UorD==0 ) p->SetAdcUp(data);
      else          p->SetAdcDown(data);
      break;
    case 1:
      if( UorD==0 ) p->SetTdcUp(data);
      else          p->SetTdcDown(data);
      break;
    default:
      hddaq::cerr << func_name << " wrong AorT " << std::endl
		  << "DetectorId = " << id    << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "SegId      = " << seg   << std::endl
		  << "AorT       = " << AorT  << std::endl
		  << "UorD       = " << UorD  << std::endl;
      return false;
    }
    return true;
  }
  //______________________________________________________________________________
  inline bool
  AddSDDRawHit( SDDRHitContainer& cont,
		int id, int port, int unit, int type, int data )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");    
    SDDRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      SDDRawHit *q = cont[i];
      if(q->PortId()==port && q->UnitId()==unit ){
	p=q; break;
      }
    }
    if( !p ){
      p = new SDDRawHit( id, port, unit );
      cont.push_back(p);
    }
    switch(type){
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
      p->SetAdc(type,data);
      break;
    case 8:
      p->SetLeading(data);
      break;
    case 9:
      p->SetTrailing(data);
      break;
    case 11:
      p->SetResetLeading(data);
      break;
    case 12:
      p->SetResetTrailing(data);
      break;
    default:
      return false;
    }
    return true;
  }
  //______________________________________________________________________________
  inline bool
  AddACRawHit( HodoRHitContainer& cont,
	       int id, int plane, int seg, int type, int data, bool DEBUG=false )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");
    
    HodoRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      HodoRawHit *q = cont[i];
      if( q->DetectorId()==id &&
	  q->PlaneId()==plane &&
	  q->SegmentId()==seg ){
	p=q; break;
      }
    }
    if( !p ){
      p = new HodoRawHit( id, plane, seg );
      cont.push_back(p);
    }
    switch(type){
    case 0:
    case 1:
    case 2:
    case 3:
      p->SetAdcUp(data);
      break;
    case 4:
      p->SetTdcUp(data);
      break;
    default:
      hddaq::cerr << func_name << " wrong type " << std::endl
		  << "DetectorId = " << id    << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "SegId      = " << seg   << std::endl
		  << "type       = " << type  << std::endl;
      return false;
    }
    return true;
  }
  //______________________________________________________________________________
  inline bool
  AddDCRawHit( DCRHitContainer& cont,
	       int plane, int wire, int tdc, int type=kLeading )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");
    DCRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      DCRawHit *q = cont[i];
      if( q->PlaneId()==plane &&
	  q->WireId()==wire ){
	p=q; break;
      }
    }
    if( !p ){
      p = new DCRawHit( plane, wire );
      cont.push_back(p);
    }

    switch(type){
    case kLeading:
      p->SetTdc(tdc);
      break;
    case kTrailing:
      p->SetTrailing(tdc);
      break;
    default:
      hddaq::cerr << func_name << " wrong data type " << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "WireId     = " << wire  << std::endl
		  << "DataType   = " << type  << std::endl;
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  inline bool
  AddMTDCRawHit( MTDCRHitContainer& cont,
		 int id, int seg, int type, int data )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");
    MTDCRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      MTDCRawHit *q = cont[i];
      if(q->SegmentId()==seg ){
	p=q; break;
      }
    }
    if( !p ){
      p = new MTDCRawHit( id, 0, seg );
      cont.push_back(p);
    }
    
    switch(type){
    case kLeading:
      p->SetLeading(data);
      break;
    case kTrailing:
      p->SetTrailing(data);
      break;
    default:
      hddaq::cerr << func_name << " wrong data type " << std::endl
		  << "SegmentId  = " << seg  << std::endl
		  << "DataType   = " << type  << std::endl;
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  inline void
  DecodeHodo( int id, int plane, int nseg, int nch, HodoRHitContainer& cont, bool DEBUG=false )
  {
    for( int seg=0; seg<nseg; ++seg ){
      for( int UorD=0; UorD<nch; ++UorD ){
	for( int AorT=0; AorT<2; ++AorT ){
	  int nhit = gUnpacker.get_entries( id, plane, seg, UorD, AorT );
	  if( nhit<=0 ) continue;
	  for(int ihit=0;ihit<nhit;++ihit){
	    int data = gUnpacker.get( id, plane, seg, UorD, AorT, ihit );
	    AddHodoRawHit( cont, id, plane, seg, UorD, AorT, data, DEBUG );
	  }
	}
      }
    }
  }
  //______________________________________________________________________________
  inline void
  DecodeBHT( int id, int plane, int nseg, int nch, HodoRHitContainer& cont, bool DEBUG=false )
  {
    for( int seg=0; seg<nseg; ++seg ){
      for( int UorD=0; UorD<nch; ++UorD ){
	int nhit = gUnpacker.get_entries( id, plane, seg, UorD, kLeading ); //leading
	if( nhit<=0 ) continue;
	for(int ihit=0;ihit<nhit;++ihit){
	  int data = gUnpacker.get( id, plane, seg, UorD, kLeading, ihit );
	  AddHodoRawHit( cont, id, plane, seg, UorD, 1, data, DEBUG );
	}
	nhit = gUnpacker.get_entries( id, plane, seg, UorD, kTrailing ); //trailing
	if( nhit<=0 ) continue;
	for(int ihit=0;ihit<nhit;++ihit){
	  int data = gUnpacker.get( id, plane, seg, UorD, kTrailing, ihit );
	  AddHodoRawHit( cont, id, plane, seg, UorD, 0, data, DEBUG );
	}
      }
    }
  }
  //______________________________________________________________________________
  inline void
  DecodeSDD( int id, int nport, SDDRHitContainer& cont )
  {
    for( int port=0; port<nport; ++port ){
      for( int unit=0; unit<4; ++unit ){
	for( int type=0; type<14; ++type ){
	  int nhit = gUnpacker.get_entries( id, port, 0, unit, type );
	  if( nhit<=0 ) continue;
	  for(int ihit=0;ihit<nhit;++ihit){
	    int data = gUnpacker.get(id, port, 0, unit,  type, ihit );
	    AddSDDRawHit( cont, id, port, unit ,type, data );
	  }
	}
      }
    }
  }
  inline void
  DecodeAC( int id, int nadc, HodoRHitContainer& cont )
  {
    for( int type=0; type<nadc+1; ++type ){
      int nhit = gUnpacker.get_entries( id, 0, 0, 0, type );
      if( nhit<=0 ) continue;
      for(int ihit=0;ihit<nhit;++ihit){
	int data = gUnpacker.get(id, 0, 0, 0,  type, ihit );
	AddACRawHit( cont, id, 0, 0 ,type, data );
      }
    }
  }
  //______________________________________________________________________________
  inline void
  DecodeMTDC( int id, int nseg, MTDCRHitContainer& cont )
  {
    for( int seg=0; seg<nseg; ++seg ){
      for( int type=0; type<2; ++type ){
	int nhit = gUnpacker.get_entries( id, 0, seg, 0, type );
	//	std::cout<<"DecodeMTDC  "<<id<<"  "<<seg<<"  "<<type<<"  "<<nhit<<std::endl;
	if( nhit<=0 ) continue;
	for(int ihit=0;ihit<nhit;++ihit){
	  int data = gUnpacker.get(id, 0, seg, 0, type, ihit );
	  AddMTDCRawHit( cont, id, seg,type, data );
	}
      }
    }
  }
  //______________________________________________________________________________
  inline void
  DecodeDC( int id, int layer, int nwire, DCRHitContainer& cont)
  {
    for(int wire=0; wire<nwire; ++wire){
#ifndef E62 // large BPC
      //      if(id==DetIdBPC&&(wire<6||wire>25)) continue;
#endif
      int nhit = gUnpacker.get_entries( id, layer, 0, wire, 0 );
#if OscillationCut
      if( nhit>MaxMultiHitDC ) continue;
#endif
      for(int i=0; i<nhit; i++ ){
	int data = gUnpacker.get( id, layer, 0, wire, 0, i);
	//	if(id==DetIdVFT) std::cout<<id<<"  "<<layer<<"  "<<wire<<"  "<<i<<"  "<<data<<std::endl;
	//	if( data<tdcmin || tdcmax<data ) continue;
	AddDCRawHit( cont, layer, wire, data );
      }
      // trailing edge
      nhit = gUnpacker.get_entries( id, layer, 0, wire, 1 );
      for(int i=0; i<nhit; i++ ){
	int data = gUnpacker.get( id, layer, 0, wire, 1, i);
	//	if(id==DetIdFDC) std::cout<<id<<"  "<<layer<<"  "<<wire<<"  "<<i<<"  "<<data<<std::endl;
	//	if( data<tdcmin || tdcmax<data ) continue;
	AddDCRawHit( cont, layer, wire, data, kTrailing );
      }
    }
  }
  
  //______________________________________________________________________________
  inline void
  DecodeHodo( int id, int nseg, int nch, HodoRHitContainer& cont )
  {
    DecodeHodo( id, 0, nseg, nch, cont );
  }
  
}

//______________________________________________________________________________
RawData::RawData( void )
  : m_is_decoded(false),
    m_BHDRawHC(),
    m_T0RawHC(),
    m_T1RawHC(),
    m_E0RawHC(),
    m_DEFRawHC(),
    m_ACRawHC(),
    m_SDDRawHC(),
    m_TrigRawHC(),
    m_BLC1aRawHC(8),
    m_BLC1bRawHC(8),
    m_BLC2aRawHC(8),
    m_BLC2bRawHC(8),
    m_BPC2RawHC(8),
#if CDS
    m_CDHRawHC(),
    m_CDCRawHC(118),
#endif
#ifdef E15
    m_FDCRawHC(6),
    m_BVCRawHC(),
    m_CVCRawHC(),
    m_NCRawHC(),
    m_PCRawHC(),
    m_LBRawHC(),
    m_WVCRawHC(),
    m_BDRawHC(),
    m_BPDRawHC(),
    m_IHRawHC(),
#elif E62
    m_StartRawHC(),
    m_StopRawHC(),
    m_SDCRawHC(8),
    m_FDCRawHC(6),
#elif E73
    m_PbF2RawHC(),
    m_Veto0RawHC(),
    m_VetoRawHC(),
    m_BTCRawHC(),
    m_LeakRawHC(),
#elif T98
    m_PbGRawHC(),
    m_PbF2RawHC(),
    m_VetoRawHC(),
    m_BTCRawHC(),
    m_RCRawHC(),
    m_VFTRawHC(14),
#elif E73_2024
    m_PbGRawHC(),
    m_PbF2RawHC(),
    m_VetoRawHC(),
    m_BTCRawHC(),
    m_CVCRawHC(),
    m_NCRawHC(),
    m_BPC1RawHC(8),
    m_VFTRawHC(14),
#endif
    m_VmeCalibRawHC()
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
RawData::~RawData( void )
{
  ClearAll();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
RawData::ClearAll( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::ClearContainer( m_SDDRawHC );
  del::ClearContainer( m_TrigRawHC );
  del::ClearContainer( m_BHDRawHC );
  del::ClearContainer( m_T0RawHC );
  del::ClearContainer( m_E0RawHC );
  del::ClearContainer( m_DEFRawHC );
  del::ClearContainer( m_ACRawHC );
  del::ClearContainer( m_T1RawHC );
  del::ClearContainerAll( m_BLC1aRawHC );
  del::ClearContainerAll( m_BLC1bRawHC );
  del::ClearContainerAll( m_BLC2aRawHC );
  del::ClearContainerAll( m_BLC2bRawHC );
  del::ClearContainerAll( m_BPC2RawHC );
#ifdef E15
  del::ClearContainer( m_BVCRawHC );
  del::ClearContainer( m_CVCRawHC );
  del::ClearContainer( m_NCRawHC  );
  del::ClearContainer( m_PCRawHC  );
  del::ClearContainer( m_LBRawHC  );
  del::ClearContainer( m_WVCRawHC );
  del::ClearContainer( m_BDRawHC  );
  del::ClearContainer( m_BPDRawHC );
  del::ClearContainer( m_IHRawHC  );
  del::ClearContainerAll( m_FDCRawHC );
#elif E62
  del::ClearContainer( m_StartRawHC );
  del::ClearContainer( m_StopRawHC );
  del::ClearContainerAll( m_SDCRawHC );
  del::ClearContainerAll( m_FDCRawHC );
#elif E73
  del::ClearContainer( m_PbF2RawHC );
  del::ClearContainer( m_Veto0RawHC );
  del::ClearContainer( m_VetoRawHC );
  del::ClearContainer( m_BTCRawHC );  
  del::ClearContainer( m_LeakRawHC );
#elif T98
  del::ClearContainer( m_PbGRawHC );
  del::ClearContainer( m_PbF2RawHC );
  del::ClearContainer( m_VetoRawHC );
  del::ClearContainer( m_BTCRawHC );  
  del::ClearContainer( m_RCRawHC );
  del::ClearContainerAll( m_VFTRawHC );
#elif E73_2024
  del::ClearContainer( m_PbGRawHC );
  del::ClearContainer( m_PbF2RawHC );
  del::ClearContainer( m_VetoRawHC );
  del::ClearContainer( m_BTCRawHC );  
  del::ClearContainer( m_CVCRawHC );
  del::ClearContainer( m_NCRawHC );
  del::ClearContainerAll( m_BPC1RawHC );
  del::ClearContainerAll( m_VFTRawHC );
#endif
#ifdef CDS
  del::ClearContainer( m_CDHRawHC );
  del::ClearContainerAll( m_CDCRawHC );
#endif
}

//______________________________________________________________________________
bool
RawData::DecodeMTDCHits( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  del::ClearContainer( m_TrigRawHC );
  DecodeMTDC( DetIdTrigFlag, 32 , m_TrigRawHC);
  return true;
}
//______________________________________________________________________________
bool
RawData::DecodeHits( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_decoded ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded!" << std::endl;
    return false;
  }
  ClearAll();
  DecodeHodoHits();
  DecodeDCHits();
  DecodeMTDCHits();
  DecodeCDCHits();
  m_is_decoded = true;
  return true;
}
//______________________________________________________________________________
bool
RawData::DecodeHodoHits( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
#if T98
  DecodeBHT(	DetIdBHD	,  0, 63, kBothSide, m_BHDRawHC );
#elif E73_2024
  DecodeBHT(	DetIdBHD	,  0, 63, kBothSide, m_BHDRawHC );
#else
  DecodeHodo(	DetIdBHD	,  0, 16, kBothSide, m_BHDRawHC );
#endif
  DecodeHodo(	DetIdT0		,  0,  5, kBothSide, m_T0RawHC );
  DecodeHodo(	DetIdT1	        ,  0,  1, kBothSide, m_T1RawHC );
  DecodeHodo(	DetIdE0		,  0,  3, kBothSide, m_E0RawHC );
  DecodeAC(	DetIdAC		,  4, m_ACRawHC );
#ifdef E62
  DecodeHodo(	DetIdDEF	,  0,  3, kBothSide, m_DEFRawHC );
  DecodeSDD(	DetIdSDD	,  1, m_SDDRawHC);
  DecodeHodo(   DetIdStart      ,  0,  8, kBothSide, m_StartRawHC );
  DecodeHodo(   DetIdStop       ,  0,  5, kBothSide, m_StopRawHC );
#elif E57
  DecodeHodo(	DetIdDEF	,  0,  5, kBothSide, m_DEFRawHC );
  DecodeSDD(	DetIdSDD	, 12, m_SDDRawHC);
#elif E73
  DecodeHodo(	DetIdDEF	,  0,  5, kBothSide, m_DEFRawHC );
  DecodeHodo(	DetIdPbF2	,  0, 40, kOneSide, m_PbF2RawHC );
  DecodeHodo(	DetIdVeto 	,  0,  2, kBothSide,m_VetoRawHC );
  DecodeHodo(	DetIdVeto0	,  0,  1, kOneSide, m_Veto0RawHC );
  DecodeHodo(	DetIdBTC	,  0,  2, kOneSide, m_BTCRawHC );
  DecodeHodo(	DetIdLeak	,  0,  6, kOneSide, m_LeakRawHC );
  //  DecodeHodo(	DetIdFinger	,  0,  2, kOneSide, m_FingerRawHC );
#elif T98
  DecodeHodo(	DetIdDEF	,  0,  5, kBothSide, m_DEFRawHC );
  DecodeHodo(	DetIdPbG	,  0, 40, kOneSide,  m_PbGRawHC );
  DecodeHodo(	DetIdPbF2	,  0, 36, kOneSide,  m_PbF2RawHC );
  DecodeHodo(	DetIdVeto 	,  0,  4, kBothSide, m_VetoRawHC );
  DecodeHodo(	DetIdBTC	,  0,  2, kBothSide, m_BTCRawHC );
  DecodeHodo(	DetIdRC	        ,  0,  8, kBothSide, m_RCRawHC );
#elif E73_2024
  DecodeHodo(	DetIdDEF	,  0,  4, kBothSide, m_DEFRawHC );
  DecodeHodo(	DetIdPbG	,  0, 40, kOneSide,  m_PbGRawHC );
  DecodeHodo(	DetIdPbF2	,  0, 36, kOneSide,  m_PbF2RawHC );
  DecodeHodo(	DetIdVeto 	,  0,  4, kBothSide, m_VetoRawHC );
  DecodeHodo(	DetIdBTC	,  0,  2, kBothSide, m_BTCRawHC );
  DecodeHodo(	DetIdCVC	,  0, 10, kBothSide, m_CVCRawHC );
  DecodeHodo(	DetIdNC	        ,  0,  6, kBothSide, m_NCRawHC );
#endif
#ifdef CDS
  DecodeHodo(	DetIdCDH	,  0, 36, kBothSide, m_CDHRawHC );
#ifdef E15
  DecodeHodo(	DetIdIH 	,  0, 24, kOneSide, m_IHRawHC );
#endif
#endif
  return true;
}
//______________________________________________________________________________
bool
RawData::DecodeDCHits( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  for(int i=0;i<8;i++){
    DecodeDC( DetIdBLC1a, i, 32, m_BLC1aRawHC[i] );
    DecodeDC( DetIdBLC1b, i, 32, m_BLC1bRawHC[i] );
    DecodeDC( DetIdBLC2a, i, 32, m_BLC2aRawHC[i] );
    DecodeDC( DetIdBLC2b, i, 32, m_BLC2bRawHC[i] );
#ifdef E62
    DecodeDC( DetIdBPC2  , i, 16, m_BPC2RawHC[i] );
    DecodeDC( DetIdSDC   , i, 16, m_SDCRawHC[i] );
    if(i<6){
      DecodeDC( DetIdFDC, i, 64, m_FDCRawHC[i] );
    }
#elif E57 
    DecodeDC( DetIdBPC2  , i, 32, m_BPC2RawHC[i] );
#elif E73 
    DecodeDC( DetIdBPC2  , i, 32, m_BPC2RawHC[i] );
#elif T98
    DecodeDC( DetIdBPC2  , i, 32, m_BPC2RawHC[i] );
    //    DecodeDC( DetIdBPC2 , i, 15, m_BPC2RawHC[i] );
  }
  for(int i=0;i<14;i++){
    DecodeDC( DetIdVFT  , i, 64, m_VFTRawHC[i] ); 
#elif E73_2024
    DecodeDC( DetIdBPC1 , i, 15, m_BPC1RawHC[i] );
    DecodeDC( DetIdBPC2 , i, 32, m_BPC2RawHC[i] );
  }
  for(int i=0;i<14;i++){
    DecodeDC( DetIdVFT  , i, 64, m_VFTRawHC[i] ); 
#endif    
  }
  return true;
}
//______________________________________________________________________________
bool
RawData::DecodeCDCHits( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
#ifdef CDS
  for(int i=0;i<118;i++){
    DecodeDC( DetIdCDC  , i, 16, m_CDCRawHC[i] );
  }
#endif
  return true;
}
//______________________________________________________________________________
bool
RawData::DecodeCalibHits( void )
{
  return true;
}

//______________________________________________________________________________
void
RawData::PrintHodo( const int &detid )
{
  std::string name="";
  switch(detid){
  case DetIdBHD:   name="BHD";   break;
  case DetIdT0:    name="T0";    break;
  case DetIdT1:    name="T1";    break;
  case DetIdE0:    name="E0";    break;
  case DetIdDEF:   name="DEF";   break;
  case DetIdCDH:   name="CDH";   break;
  default: name="ERROR";
  }
  const HodoRHitContainer &cont = GetHodoRawHC(detid);
  int nh = cont.size();
  std::cout<<"------------------------"<<std::endl;
  std::cout<<"[ "<<name<<" Raw Hit ] size: "<<nh<<std::endl;
  for( int i=0; i<nh; ++i ){
    HodoRawHit *raw = cont[i];
    if(!raw) continue;
    int seg = raw->SegmentId();
    double au = raw->GetAdcUp();
    double ad = raw->GetAdcDown();  
    int ntu = raw->GetSizeTdcUp();
    int ntd = raw->GetSizeTdcDown();
    std::cout<<"[Segment "<<std::setw(2)<<seg<<"]"<<std::endl;
    std::cout<<"   ADCud  :"<<std::setw(8)<<au<<std::setw(8)<<ad<<std::endl;
    std::cout<<"   TDCup  :["<<std::setw(2)<<ntu<<"]";      
    for(int it=0;it<ntu;it++){   
      int tu  = raw->GetTdcUp(it);
      std::cout<<std::setw(15)<<tu;
    }
    std::cout<<std::endl;
    std::cout<<"   TDCdown:["<<std::setw(2)<<ntd<<"]";      
    for(int it=0;it<ntd;it++){   
      int td  = raw->GetTdcDown(it);
      std::cout<<std::setw(15)<<td;
    }
    std::cout<<std::endl;
  }
  std::cout<<"------------------------"<<std::endl;
  return;
}
