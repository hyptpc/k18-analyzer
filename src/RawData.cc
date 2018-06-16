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
#include "DCRawHit.hh"
#include "HodoRawHit.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

#define OscillationCut 0
#define E07CntIs 0
#define E07SsdIs 0

namespace
{
  using namespace hddaq::unpacker;
  const std::string& class_name("RawData");
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
  const UserParamMan&    gUser     = UserParamMan::GetInstance();
  enum EUorD { kOneSide=1, kBothSide=2 };
  enum EHodoDataType { kHodoAdc, kHodoLeading, kHodoTrailing, kHodoOverflow, kHodoNDataType };
  enum EDCDataType   { kDcLeading, kDcTrailing, kDcOverflow, kDcNDataType };
#if OscillationCut
  const int  MaxMultiHitDC  = 16;
#endif

  //______________________________________________________________________________
  inline bool
  AddHodoRawHit( HodoRHitContainer& cont,
		 int id, int plane, int seg, int UorD, int type, int data )
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
    case kHodoAdc:
      if( UorD==0 ) p->SetAdcUp(data);
      else          p->SetAdcDown(data);
      break;
    case kHodoLeading:
      if( UorD==0 ) p->SetTdcUp(data);
      else          p->SetTdcDown(data);
      break;
    case kHodoTrailing:
      if( UorD==0 ) p->SetTdcTUp(data);
      else          p->SetTdcTDown(data);
      break;
    case kHodoOverflow:
      p->SetTdcOverflow(data);
      break;
    default:
      hddaq::cerr << func_name << " wrong data type " << std::endl
		  << "DetectorId = " << id    << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "SegId      = " << seg   << std::endl
		  << "AorT       = " << type  << std::endl
		  << "UorD       = " << UorD  << std::endl;
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  inline bool
  AddDCRawHit( DCRHitContainer& cont,
	       int plane, int wire, int data, int type=kDcLeading )
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
    case kDcLeading:
      p->SetTdc(data);
      break;
    case kDcTrailing:
      p->SetTrailing(data);
      break;
    case kDcOverflow:
      p->SetTdcOverflow(data);
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
  inline void
  DecodeHodo( int id, int plane, int nseg, int nch, HodoRHitContainer& cont )
  {
    for( int seg=0; seg<nseg; ++seg ){
      for( int UorD=0; UorD<nch; ++UorD ){
	for( int AorT=0; AorT<2; ++AorT ){
	  int nhit = gUnpacker.get_entries( id, plane, seg, UorD, AorT );
	  if( nhit<=0 ) continue;
	  int data = gUnpacker.get( id, plane, seg, UorD, AorT );
	  AddHodoRawHit( cont, id, plane, seg, UorD, AorT, data );
	}
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
    m_T1RawHC(),
    m_T2RawHC(),
    m_T3RawHC(),
    m_T4RawHC(),
    m_S1RawHC(),
    m_S2RawHC(),
    m_BH1RawHC(),
    m_BH2RawHC(),
    m_BACRawHC(),
    m_PVACRawHC(),
    m_FACRawHC(),
    m_SACRawHC(),
    m_TOFRawHC(),
    m_HtTOFRawHC(),
    m_LCRawHC(),
    m_BFTRawHC(NumOfPlaneBFT),
    m_SFTRawHC(NumOfPlaneSFT),
    m_SCHRawHC(),
    m_FBHRawHC(),
    m_SSDTRawHC(),
    m_BcInRawHC(),
    m_BcOutRawHC(NumOfLayersBcOut+1),
    m_SdcInRawHC(NumOfLayersSdcIn+1),
    m_SdcOutRawHC(NumOfLayersSdcOut+1),
    m_SsdInRawHC(NumOfLayersSsdIn+1),
    m_SsdOutRawHC(NumOfLayersSsdOut+1),
    m_ScalerRawHC(),
    m_TrigRawHC(),
    m_VmeCalibRawHC(),
    m_FpgaBH1RawHC(),
    m_FpgaBH2RawHC(),
    m_FpgaBH2MtRawHC()
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

  del::ClearContainer( m_T1RawHC );
  del::ClearContainer( m_T2RawHC );
  del::ClearContainer( m_T3RawHC );
  del::ClearContainer( m_T4RawHC );
  del::ClearContainer( m_S1RawHC );
  del::ClearContainer( m_S2RawHC );
  del::ClearContainer( m_BH1RawHC );
  del::ClearContainer( m_BH2RawHC );
  del::ClearContainer( m_BACRawHC );
  del::ClearContainer( m_PVACRawHC );
  del::ClearContainer( m_FACRawHC );
  del::ClearContainer( m_SACRawHC );
  del::ClearContainer( m_TOFRawHC );
  del::ClearContainer( m_HtTOFRawHC );
  del::ClearContainer( m_LCRawHC );

  del::ClearContainerAll( m_BFTRawHC );
  del::ClearContainerAll( m_SFTRawHC );

  del::ClearContainer( m_SCHRawHC );
  del::ClearContainer( m_FBHRawHC );
  del::ClearContainer( m_SSDTRawHC );

  // del::ClearContainerAll( m_BcInRawHC );
  del::ClearContainerAll( m_BcOutRawHC );
  del::ClearContainerAll( m_SdcInRawHC );
  del::ClearContainerAll( m_SdcOutRawHC );
  del::ClearContainerAll( m_SsdInRawHC );
  del::ClearContainerAll( m_SsdOutRawHC );

  del::ClearContainer( m_ScalerRawHC );
  del::ClearContainer( m_TrigRawHC );
  del::ClearContainer( m_VmeCalibRawHC );

  del::ClearContainer( m_FpgaBH1RawHC );
  del::ClearContainer( m_FpgaBH2RawHC );
  del::ClearContainer( m_FpgaBH2MtRawHC );
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

  // T1
  DecodeHodo( DetIdT1, NumOfSegT1, kBothSide, m_T1RawHC );
  // T2
  DecodeHodo( DetIdT2, NumOfSegT2, kBothSide, m_T2RawHC );
  // T3
  DecodeHodo( DetIdT3, NumOfSegT3, kBothSide, m_T3RawHC );
  // T4
  DecodeHodo( DetIdT4, NumOfSegT4, kBothSide, m_T4RawHC );
  // S1
  DecodeHodo( DetIdS1, NumOfSegS1, kBothSide, m_S1RawHC );
  // S2
  DecodeHodo( DetIdS2, NumOfSegS2, kBothSide, m_S2RawHC );

  //SSDT
  DecodeHodo( DetIdSSDT, NumOfSegSSDT, kOneSide, m_SSDTRawHC );
  for( int seg=0; seg<NumOfSegSSDT; ++seg ){
    int nhit = gUnpacker.get_entries( DetIdSSDT, 0, seg*2, 0, 1 );
    if( nhit<=0 ) continue;
    int data = gUnpacker.get( DetIdSSDT, 0, seg*2, 0, 1 );
    AddHodoRawHit( m_SSDTRawHC, DetIdSSDT, 0, seg, 0, 1, data );
  }

  // SsdIn (SSD1)
  for( int plane=0; plane<NumOfLayersSsdIn; ++plane ){
  //   for( int type=0; type<kDcNDataType; ++type ){
    for( int type=0; type<kDcNDataType - 1; ++type ){
      for( int seg=0; seg<NumOfSegSSD1; ++seg ){
  	int nhit = gUnpacker.get_entries( DetIdSSD1, plane, seg, 0, type );
  	for( int i=0; i<nhit; ++i ){
  	  int data = gUnpacker.get( DetIdSSD1, plane, seg, 0, type, i );
  	  AddDCRawHit( m_SsdInRawHC[plane+1], plane+PlMinSsdIn, seg+1, data, type );
  	}//for(i)
      }//for(seg)
    }//for(type)
  }//for(plane)

  // SsdOut (SSD2)
  for( int plane=0; plane<NumOfLayersSsdOut; ++plane ){
  //   for( int type=0; type<kDcNDataType; ++type ){
    for( int type=0; type<kDcNDataType - 1; ++type ){
      for( int seg=0; seg<NumOfSegSSD2; ++seg ){
  	int nhit = gUnpacker.get_entries( DetIdSSD2, plane, seg, 0, type );
  	for(int i=0; i<nhit; ++i){
  	  int data = gUnpacker.get( DetIdSSD2, plane, seg, 0, type, i );
  	  AddDCRawHit( m_SsdOutRawHC[plane+1], plane+PlMinSsdOut, seg+1, data, type );
  	}//for(i)
      }//for(seg)
    }//for(type)
  }//for(plane)

  // Scaler
  for( int l = 0; l<NumOfScaler; ++l){
    for( int seg=0; seg<NumOfSegScaler; ++seg ){
      int nhit = gUnpacker.get_entries( DetIdScaler, l, 0, seg, 0 );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdScaler, 0, 0, seg, 0 );
	AddHodoRawHit( m_ScalerRawHC, DetIdScaler, l, seg, 0, 0, data );
      }
    }
  }// for(l)

  // trigger Flag ---------------------------------------------------------------
  DecodeHodo( DetIdTrig, NumOfSegTrig, kOneSide, m_TrigRawHC );

  m_is_decoded = true;
  return true;
}

//______________________________________________________________________________
bool
RawData::DecodeCalibHits( void )
{
  del::ClearContainer( m_VmeCalibRawHC );

  for( int plane=0; plane<NumOfPlaneVmeCalib; ++plane ){
    DecodeHodo( DetIdVmeCalib, plane, NumOfSegVmeCalib,
  		kOneSide, m_VmeCalibRawHC );
  }

  return true;
}
