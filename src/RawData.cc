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
#define Run2019 0

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
	  for(int m = 0; m<nhit; ++m){
	    int data = gUnpacker.get( id, plane, seg, UorD, AorT, m );
	    AddHodoRawHit( cont, id, plane, seg, UorD, AorT, data );
	  }
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
    m_SACRawHC(),
    m_TOFRawHC(),
    m_HtTOFRawHC(),
    m_LACRawHC(),
    m_LCRawHC(),
    m_BFTRawHC(NumOfPlaneBFT),
    m_SFTRawHC(NumOfPlaneSFT),
    m_CFTRawHC(NumOfPlaneCFT),
    m_BGORawHC(),
    m_BGOFadcRawHC(NumOfSegBGO),
    m_PiIDRawHC(),
    m_SCHRawHC(),
    m_FBT1RawHC(2*NumOfLayersFBT1),
    m_FBT2RawHC(2*NumOfLayersFBT2),
    m_BcInRawHC(NumOfLayersBcIn+1),
    m_BcOutRawHC(NumOfLayersBcOut+1),
    m_SdcInRawHC(NumOfLayersSdcIn+1),
    m_SdcOutRawHC(NumOfLayersSdcOut+1),
    m_ScalerRawHC(),
    m_TrigRawHC(),
    m_VmeCalibRawHC(),
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
  del::ClearContainer( m_SACRawHC );
  del::ClearContainer( m_TOFRawHC );
  del::ClearContainer( m_HtTOFRawHC );
  del::ClearContainer( m_LACRawHC );
  del::ClearContainer( m_LCRawHC );

  del::ClearContainerAll( m_BFTRawHC );
  del::ClearContainerAll( m_SFTRawHC );
  del::ClearContainerAll( m_CFTRawHC );
  del::ClearContainer( m_BGORawHC );
  //del::ClearContainerAll( m_BGOFadcRawHC );
  for (int i=0; i<NumOfSegBGO; i++) {
    m_BGOFadcRawHC[i].clear();
    m_BGOFadcRawHC[i].shrink_to_fit();
  }



  del::ClearContainer( m_PiIDRawHC );
  del::ClearContainerAll( m_FBT1RawHC );
  del::ClearContainerAll( m_FBT2RawHC );

  del::ClearContainer( m_SCHRawHC );

  del::ClearContainerAll( m_BcInRawHC );
  del::ClearContainerAll( m_BcOutRawHC );
  del::ClearContainerAll( m_SdcInRawHC );
  del::ClearContainerAll( m_SdcOutRawHC );

  del::ClearContainer( m_ScalerRawHC );
  del::ClearContainer( m_TrigRawHC );
  del::ClearContainer( m_VmeCalibRawHC );

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
  // TOF-HT
  DecodeHodo( DetIdHtTOF,NumOfSegHtTOF,kOneSide,  m_HtTOFRawHC );
#if Run2019
  // LAC
  DecodeHodo( DetIdLAC,  NumOfSegLAC,  kOneSide,  m_LACRawHC );
#endif
  // LC
  DecodeHodo( DetIdLC,  NumOfSegLC,  kOneSide,  m_LCRawHC );

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

  //CFT
  for( int plane=0; plane<NumOfPlaneCFT; ++plane ){    
    //CFT TDC 
    for(int seg = 0; seg<NumOfSegCFT[plane]; ++seg){
      int nhit_tdc      = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 0 );
      int nhit_trailing = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 1 );
      if( nhit_tdc>0 ){
	for(int i = 0; i<nhit_tdc; ++i){
	  int leading  = gUnpacker.get( DetIdCFT, plane, seg, 0, 0, i );
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 0, kHodoLeading, leading );
	}
      }	
      if(nhit_trailing>0){
	for(int i = 0; i<nhit_trailing; ++i){
	  int trailing = gUnpacker.get( DetIdCFT, plane, seg, 0, 1, i )  ;
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 0, kHodoTrailing, trailing );
	}	  	
      }
      //if( nhit_tdc==0 )continue; // w/ or w/o TDC,   comment out => pedestal              

      //CFT ADC HIGH
      int nhit_adc_hi = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 2 );
      if( nhit_adc_hi>0 ){
	for(int i = 0; i<nhit_adc_hi; ++i){
	  int adc_hi  = gUnpacker.get( DetIdCFT, plane, seg, 0, 2, i )  ;
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 0, 0, adc_hi );//ADC Hi	  
	}
      }
      else continue;
      //CFT ADC LOW
      int nhit_adc_low = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 3 );
      if( nhit_adc_low>0 ){
	for(int i = 0; i<nhit_adc_low; ++i){
	  int adc_low = gUnpacker.get( DetIdCFT, plane, seg, 0, 3, i );
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 1, 0, adc_low );//ADC Low
	}
      }
      else continue;      
    }
  }
  
  //BGO
  for(int seg=0; seg<NumOfSegBGO; ++seg){
    int ped      = 0;
    int integral = 0;

    int nhit_a = gUnpacker.get_entries( DetIdBGO, 0, seg, 0, 0 );
    if( nhit_a>0 ){
      for(int i = 0; i<nhit_a; ++i){
	unsigned int fadc = gUnpacker.get(DetIdBGO, 0, seg, 0, 0 ,i);	

	m_BGOFadcRawHC[seg].push_back(fadc);

	if( fadc == 0xffff ){
	}else{
	  if( ped == 0 ){ped = fadc;} // 1st value is seemed to be pedestal
	  int delta = ped - fadc;
	  /*
	  std::cout << "BGO : seg=" << seg  << ", i=" << i << "/nhit=" << nhit
		    << ", fadc=" << fadc << ", pede=" << ped 
		    << ", integral=" << integral << ", +=" << delta // (ped - fadc)
		    << std::endl;
	  */
	  integral += delta;
	}	
      }
      AddHodoRawHit( m_BGORawHC, DetIdBGO, 0, seg , 0, kHodoAdc, integral );
    }

    //TDC
    unsigned int nhit_t = gUnpacker.get_entries(DetIdBGO, 0, seg, 0, 1);
    for(unsigned int m = 0; m<nhit_t; ++m){
      int tdc = gUnpacker.get(DetIdBGO, 0, seg, 0, kHodoLeading, m);
      AddHodoRawHit( m_BGORawHC, DetIdBGO, 0, seg , 0, kHodoLeading, tdc );
      //std::cout << "BGO TDC : seg=" << seg  << ", tdc=" << tdc	<< std::endl;
    }

  }

  //PiID counter
  for(int seg=0; seg<NumOfSegPiID; ++seg){
    int nhit_l = gUnpacker.get_entries( DetIdPiID, 0, seg, 0, 0 );
    int nhit_t = gUnpacker.get_entries( DetIdPiID, 0, seg, 0, 1 );
    if( nhit_l>0 ){
      for(int i = 0; i<nhit_l; ++i){
	int leading  = gUnpacker.get( DetIdPiID, 0, seg, 0, 0, i );
	AddHodoRawHit( m_PiIDRawHC, DetIdPiID, 0, seg , 0, kHodoLeading,  leading );
	//std::cout << "PiID TDC : seg=" << seg  << ", tdc=" << leading << std::endl;
      }
    }
    if( nhit_t>0 ){
      for(int j = 0; j<nhit_t; ++j){	    
	int trailing = gUnpacker.get( DetIdPiID, 0, seg, 0, 1, j );
	AddHodoRawHit( m_PiIDRawHC, DetIdPiID, 0, seg , 0, kHodoTrailing, trailing );
	//std::cout << "PiID TDC trailing : seg=" << seg  << ", tdc=" << trailing << std::endl;
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

  //FBT1
  for( int layer=0; layer<NumOfLayersFBT1; ++layer ){
    for(int seg = 0; seg<MaxSegFBT1; ++seg){
      for(int UorD = 0; UorD<2; ++UorD){
	for(int LorT = 0; LorT<2; ++LorT){
	  int nhit = gUnpacker.get_entries( DetIdFBT1, layer, seg, UorD, LorT);
	  for(int i = 0; i<nhit; ++i){
	    int time  = gUnpacker.get( DetIdFBT1, layer, seg, UorD, LorT, i )  ;
	    AddHodoRawHit( m_FBT1RawHC[2*layer + UorD], DetIdFBT1, layer, seg , UorD, kHodoLeading+LorT, time );
	  }// multihit
	}// LorT
      }// UorD
    }// seg
  }// layer

  //FBT2
  for( int layer=0; layer<NumOfLayersFBT2; ++layer ){
    for(int seg = 0; seg<MaxSegFBT2; ++seg){
      for(int UorD = 0; UorD<2; ++UorD){
	for(int LorT = 0; LorT<2; ++LorT){
	  int nhit = gUnpacker.get_entries( DetIdFBT2, layer, seg, UorD, LorT);
	  for(int i = 0; i<nhit; ++i){
	    int time  = gUnpacker.get( DetIdFBT2, layer, seg, UorD, LorT, i )  ;
	    AddHodoRawHit( m_FBT2RawHC[2*layer + UorD], DetIdFBT2, layer, seg , UorD, kHodoLeading+LorT, time );
	  }// multihit
	}// LorT
      }// UorD
    }// seg
  }// layer

  // BC3&BC4 MWDC
  for(int plane=0; plane<NumOfLayersBcOut; ++plane ){
    int nlt = 1;
#if Run2019
    nlt = 2;
#endif
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC3; ++wire){
        for ( int lt = 0; lt < nlt; ++lt ) {
            int nhit = gUnpacker.get_entries( DetIdBC3, plane, 0, wire, lt );
#if OscillationCut
            if( nhit>MaxMultiHitDC ) continue;
#endif
            for(int i=0; i<nhit; i++ ){
              int data = gUnpacker.get( DetIdBC3, plane, 0, wire, lt, i);
              if ( lt == 0 && ( data<MinTdcBC3 || MaxTdcBC3<data ) ) continue;
              if ( lt == 1 && data<MinTdcBC3 ) continue;
              AddDCRawHit( m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data, lt );
            }   // for ( i )
          } // for ( lt )
        }   // for ( wire )
      } // if ( plane )
    else{
      for(int wire=0; wire<MaxWireBC4; ++wire){
        for ( int lt = 0; lt < nlt; ++lt ) {
            int nhit = gUnpacker.get_entries( DetIdBC4, plane-NumOfLayersBc, 0, wire, lt );
#if OscillationCut
            if( nhit>MaxMultiHitDC ) continue;
#endif
            for(int i=0; i<nhit; i++ ){
              int data =  gUnpacker.get( DetIdBC4, plane-NumOfLayersBc, 0, wire, lt, i );
              if ( lt == 0 && ( data<MinTdcBC4 || MaxTdcBC4<data ) ) continue;
              if ( lt == 1 && data<MinTdcBC4 ) continue;
              AddDCRawHit( m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data, lt );
            }   // for ( i )
        }   // for ( lt )
      } // for ( wire )
    }   // if ( plane )
  } // for ( plane )

  // SdcIn (SDC1)
  for( int plane=0; plane<NumOfLayersSDC1; ++plane ){
    int nlt = 1;
#if Run2019
    nlt = 2;
#endif
    for( int wire=0; wire<MaxWireSDC1; ++wire ){
      for ( int lt = 0; lt < nlt; ++lt ) {
          int nhit = gUnpacker.get_entries( DetIdSDC1, plane, 0, wire, lt );
#if OscillationCut
          if( nhit>MaxMultiHitDC ) continue;
#endif
          for(int i=0; i<nhit; i++ ){
              int data = gUnpacker.get( DetIdSDC1, plane, 0, wire, lt, i ) ;
              if ( lt == 0 && ( data<MinTdcSDC1 || MaxTdcSDC1<data ) ) continue;
              if ( lt == 1 && data<MinTdcSDC1 ) continue;
	// AddDCRawHit( m_SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data );
              AddDCRawHit( m_SdcInRawHC[plane+1], plane+1, wire+1, data, lt );
          } // for ( i )
      } // for ( lt )
    } // for ( wire )
  } // for ( plane )

  // SdcOut (SDC2&SDC3)
  for( int plane=0; plane<NumOfLayersSDC2+NumOfLayersSDC3; ++plane ){
    if( plane<NumOfLayersSDC2 ){
      for( int wire=0; wire<MaxWireSDC2; ++wire ){
	for(int lt = 0; lt<2; ++lt){
	  int nhit = gUnpacker.get_entries( DetIdSDC2, plane, 0, wire, lt );
#if OscillationCut
	  if( nhit>MaxMultiHitDC ) continue;
#endif
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC2, plane, 0, wire, lt, i );
	    if( lt == 0 && ( data<MinTdcSDC2 || MaxTdcSDC2<data ) ) continue;
	    if( lt == 1 && data<MinTdcSDC2 ) continue;
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
	    if( lt == 0 && ( data<MinTdcSDC3 || MaxTdcSDC3<data ) ) continue;
	    if( lt == 1 && data<MinTdcSDC3 ) continue;
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
	int data = gUnpacker.get( DetIdScaler, l, 0, seg, 0 );
	AddHodoRawHit( m_ScalerRawHC, DetIdScaler, l, seg, 0, 0, data );
      }
    }
  }// for(l)

  // trigger Flag ---------------------------------------------------------------
  DecodeHodo( DetIdTrig, NumOfSegTrig, kOneSide, m_TrigRawHC );

  {
    // BH2Mt
    static const int type_mt = gUnpacker.get_data_id("BH2", "fpga_meantime");
    for(int seg = 0; seg<NumOfSegBH2; ++seg){
      int mhit = gUnpacker.get_entries( DetIdBH2, 0, seg, 0, type_mt );
      for(int m = 0; m<mhit; ++m){
	int data = gUnpacker.get( DetIdBH2, 0, seg, 0, type_mt , m);
	AddHodoRawHit( m_FpgaBH2MtRawHC, DetIdFpgaBH2Mt, 0, seg, 0, kHodoLeading, data );
      }// for(m)
    }// for(seg)
  }

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
