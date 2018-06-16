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

  static const double MinTdcBC3  = gUser.GetParameter("TdcBC3", 0);
  static const double MaxTdcBC3  = gUser.GetParameter("TdcBC3", 1);
  static const double MinTdcBC4  = gUser.GetParameter("TdcBC4", 0);
  static const double MaxTdcBC4  = gUser.GetParameter("TdcBC4", 1);
  static const double MinTdcSDC1 = gUser.GetParameter("TdcSDC1", 0);
  static const double MaxTdcSDC1 = gUser.GetParameter("TdcSDC1", 1);
  static const double MinTdcSDC2 = gUser.GetParameter("TdcSDC2", 0);
  static const double MaxTdcSDC2 = gUser.GetParameter("TdcSDC2", 1);
  static const double MinTdcSDC3 = gUser.GetParameter("TdcSDC3", 0);
  static const double MaxTdcSDC3 = gUser.GetParameter("TdcSDC3", 1);

  if( m_is_decoded ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded!" << std::endl;
    return false;
  }

  ClearAll();

  // BH1
  DecodeHodo( DetIdBH1, NumOfSegBH1, kBothSide, m_BH1RawHC );
  // BH2
  DecodeHodo( DetIdBH2, NumOfSegBH2, kBothSide, m_BH2RawHC );
  // SAC
  DecodeHodo( DetIdSAC, NumOfSegSAC, kOneSide,  m_SACRawHC );
  // TOF
  DecodeHodo( DetIdTOF, NumOfSegTOF, kBothSide, m_TOFRawHC );
  // LC
  DecodeHodo( DetIdLC,  NumOfSegLC,  kOneSide,  m_LCRawHC );

#if E07CntIs
  // BAC
  DecodeHodo( DetIdBAC, NumOfSegBAC, kOneSide, m_BACRawHC );
  // PVAC
  DecodeHodo( DetIdPVAC, NumOfSegPVAC, kOneSide, m_PVACRawHC );
  // FAC
  DecodeHodo( DetIdFAC, NumOfSegFAC, kOneSide, m_FACRawHC );
#endif

  //BFT
  for( int plane=0; plane<NumOfPlaneBFT; ++plane ){
    for(int seg = 0; seg<NumOfSegBFT; ++seg){
      int nhit = gUnpacker.get_entries( DetIdBFT, plane, 0, seg, 0 );
      if( nhit>0 ){
	for(int i = 0; i<nhit; ++i){
	  int leading  = gUnpacker.get( DetIdBFT, plane, 0, seg, 0, i )  ;
	  int trailing = gUnpacker.get( DetIdBFT, plane, 0, seg, 1, i )  ;
	  AddHodoRawHit( m_BFTRawHC[plane], DetIdBFT, plane, seg , 0, kHodoLeading,  leading );
	  AddHodoRawHit( m_BFTRawHC[plane], DetIdBFT, plane, seg , 0, kHodoTrailing, trailing );
	}
      }
    }
  }

  //SFT
  for( int plane=0; plane<NumOfPlaneSFT; ++plane ){
    int nseg = 0;
    switch ( plane ) {
    case 0: nseg = NumOfSegSFT_UV;  break;
    case 1: nseg = NumOfSegSFT_UV;  break;
    case 2: nseg = NumOfSegSFT_X; break;
    case 3: nseg = NumOfSegSFT_X; break;
    default: break;
    }
    for ( int seg = 0; seg < nseg; ++seg ) {
      for ( int LorT = 0; LorT < 2; ++LorT ) {
	int nhit = gUnpacker.get_entries( DetIdSFT, plane, 0, seg, LorT );
	if( nhit>0 ){
	  for(int i = 0; i<nhit; ++i){
	    int edge = gUnpacker.get( DetIdSFT, plane, 0, seg, LorT, i )  ;
	    AddHodoRawHit( m_SFTRawHC[plane], DetIdSFT, plane, seg, 0, kHodoLeading+LorT, edge );
	  }
	}
      }
    }
  }

  //SCH
  for(int seg=0; seg<NumOfSegSCH; ++seg){
    int nhit = gUnpacker.get_entries( DetIdSCH, 0, seg, 0, 0 );
    if( nhit>0 ){
      for(int i = 0; i<nhit; ++i){
	int leading  = gUnpacker.get( DetIdSCH, 0, seg, 0, 0, i );
	int trailing = gUnpacker.get( DetIdSCH, 0, seg, 0, 1, i );
	AddHodoRawHit( m_SCHRawHC, DetIdSCH, 0, seg , 0, kHodoLeading,  leading );
	AddHodoRawHit( m_SCHRawHC, DetIdSCH, 0, seg , 0, kHodoTrailing, trailing );
      }
    }
  }

#if 0
  // BC1&BC2 MWPC
  for(int plane=0; plane<NumOfLayersBcIn; ++plane ){
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC1; ++wire){
  	int nhit = gUnpacker.get_entries( DetIdBC1, plane, 0, wire, 0 );
  	if( nhit>0 ){
  	  for(int i=0; i<nhit; i++ ){
  	    int leading = gUnpacker.get( DetIdBC1, plane, 0, wire, 0, i );
  	    int trailing = gUnpacker.get( DetIdBC1, plane, 0, wire, 1, i );
  	    AddDCRawHit( m_BcInRawHC[plane+1], plane+PlMinBcIn, wire+1,
  			 leading, kLeading );
  	    AddDCRawHit( m_BcInRawHC[plane+1], plane+PlMinBcIn, wire+1,
			 trailing, kTrailing );
  	  }
  	}
      }
    }else{
      for(int wire=0; wire<MaxWireBC2; ++wire){
      	int nhit = gUnpacker.get_entries( DetIdBC2, plane-NumOfLayersBc, 0, wire, 0 );
      	if( nhit>0 ){
  	  for(int i=0; i<nhit; i++ ){
  	    int leading = gUnpacker.get( DetIdBC2, plane-NumOfLayersBc, 0, wire, 0, i );
  	    int trailing = gUnpacker.get( DetIdBC2, plane-NumOfLayersBc, 0, wire, 1, i );
  	    AddDCRawHit( m_BcInRawHC[plane+1], plane+PlMinBcIn, wire+1,
  			 leading, kLeading);
  	    AddDCRawHit( m_BcInRawHC[plane+1], plane+PlMinBcIn, wire+1,
  			 trailing, kTrailing);
  	  }
  	}
      }
    }
  }
#endif

  // BC3&BC4 MWDC
  for(int plane=0; plane<NumOfLayersBcOut; ++plane ){
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC3; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC3, plane, 0, wire, 0 );
#if OscillationCut
	if( nhit>MaxMultiHitDC ) continue;
#endif
	for(int i=0; i<nhit; i++ ){
	  int data = gUnpacker.get( DetIdBC3, plane, 0, wire, 0, i);
	  if( data<MinTdcBC3 || MaxTdcBC3<data ) continue;
	  AddDCRawHit( m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data );
	}
      }
    }
    else{
      for(int wire=0; wire<MaxWireBC4; ++wire){
	int nhit = gUnpacker.get_entries( DetIdBC4, plane-NumOfLayersBc, 0, wire, 0 );
#if OscillationCut
	if( nhit>MaxMultiHitDC ) continue;
#endif
	for(int i=0; i<nhit; i++ ){
	  int data =  gUnpacker.get( DetIdBC4, plane-NumOfLayersBc, 0, wire, 0, i );
	  if( data<MinTdcBC4 || MaxTdcBC4<data ) continue;
	  AddDCRawHit( m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data );
	}
      }
    }
  }

#if E07SsdIs
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
#endif

  // SdcIn (SDC1)
  for( int plane=0; plane<NumOfLayersSDC1; ++plane ){
    for( int wire=0; wire<MaxWireSDC1; ++wire ){
      int nhit = gUnpacker.get_entries( DetIdSDC1, plane, 0, wire, 0 );
#if OscillationCut
      if( nhit>MaxMultiHitDC ) continue;
#endif
      for(int i=0; i<nhit; i++ ){
	int data = gUnpacker.get( DetIdSDC1, plane, 0, wire, 0, i ) ;
	if( data<MinTdcSDC1 || MaxTdcSDC1<data ) continue;
	// AddDCRawHit( m_SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data );
	AddDCRawHit( m_SdcInRawHC[plane+1], plane+1, wire+1, data );
      }
    }
  }

  // SdcOut (SDC2&SDC3)
  for( int plane=0; plane<NumOfLayersSdcOut; ++plane ){
    if( plane<NumOfLayersSDC2 ){
      for( int wire=0; wire<MaxWireSDC2; ++wire ){
	for(int lt = 0; lt<2; ++lt){
	  int nhit = gUnpacker.get_entries( DetIdSDC2, plane, 0, wire, lt );
#if OscillationCut
	  if( nhit>MaxMultiHitDC ) continue;
#endif
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC2, plane, 0, wire, lt, i );
	    if( data<MinTdcSDC2 || MaxTdcSDC2<data ) continue;
	    AddDCRawHit( m_SdcOutRawHC[plane+1], plane+PlMinSdcOut, wire+1, data , lt);
	  }// for(i)
	}// for(lt)
      }// for(wire)
    }// if(SDC2/3)
    else{
      int MaxWireSDC3;
      if( plane==NumOfLayersSDC2 || plane==(NumOfLayersSDC2+1) )
	MaxWireSDC3 = MaxWireSDC3Y;
      else
	MaxWireSDC3 = MaxWireSDC3X;
      for( int wire=0; wire<MaxWireSDC3; ++wire ){
	for(int lt = 0; lt<2; ++lt){
	  int nhit = gUnpacker.get_entries( DetIdSDC3, plane-NumOfLayersSDC2, 0, wire, lt );
#if OscillationCut
	  if( nhit>MaxMultiHitDC ) continue;
#endif
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC3, plane-NumOfLayersSDC2, 0, wire, lt ,i );
	    if( data<MinTdcSDC3 || MaxTdcSDC3<data ) continue;
	    AddDCRawHit( m_SdcOutRawHC[plane+1],  plane+PlMinSdcOut, wire+1, data , lt);
	  }// for(i)
	}// for(lt)
      }//for(wire)
    }// if(SDC2/3)
  }// for(plane)

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

  // FPGA HR-TDC ----------------------------------------------------------------
  // BH1
  for(int seg = 0; seg<NumOfSegBH1; ++seg){
    for(int ud = 0; ud<2; ++ud){
      int mhit = gUnpacker.get_entries( DetIdBH1, 0, seg, ud, 2 );
      for(int m = 0; m<mhit; ++m){
	int data = gUnpacker.get( DetIdBH1, 0, seg, ud, 2 , m);
	AddHodoRawHit( m_FpgaBH1RawHC, DetIdFpgaBH1, 0, seg, ud, kHodoLeading, data );
      }// for(m)
    }// for(ud)
  }// for(seg)

  // BH2
  for(int seg = 0; seg<NumOfSegBH2; ++seg){
    for(int ud = 0; ud<2; ++ud){
      int mhit = gUnpacker.get_entries( DetIdBH2, 0, seg, ud, 2 );
      for(int m = 0; m<mhit; ++m){
	int data = gUnpacker.get( DetIdBH2, 0, seg, ud, 2 , m);
	AddHodoRawHit( m_FpgaBH2RawHC, DetIdFpgaBH2, 0, seg, ud, kHodoLeading, data );
      }// for(m)
    }// for(ud)
  }// for(seg)

#if 1
  // BH2Mt
  for(int seg = 0; seg<NumOfSegBH2; ++seg){
    int mhit = gUnpacker.get_entries( DetIdBH2, 0, seg, 0, 6 );
    for(int m = 0; m<mhit; ++m){
      int data = gUnpacker.get( DetIdBH2, 0, seg, 0, 6 , m);
      AddHodoRawHit( m_FpgaBH2MtRawHC, DetIdFpgaBH2Mt, 0, seg, 0, kHodoLeading, data );
    }// for(m)
  }// for(seg)
#endif

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
