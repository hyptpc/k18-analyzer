// -*- C++ -*-

#ifndef DC_ANALYZER_HH
#define DC_ANALYZER_HH

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>
#include <map>
#include <string>

#include <TSystem.h>

#include "FuncName.hh"
#include "XTMapMan.hh"
#include "BLDCWireMapMan.hh"
#include "DCTdcCalibMan.hh"
#include "LocalTrack.hh"
#include "DCCluster.hh"
#include "DCHit.hh"

// class DCHit;
// class DCCluster;
 class RawData;
// class LocalTrack;

typedef std::vector<DCHit*>        DCHitContainer;
typedef std::vector<DCCluster*>    DCClusterContainer;
typedef std::vector<LocalTrack*>   LocalTrackContainer;
typedef std::vector<DCClusterContainer> DCClusterList;
typedef std::map<int, DCClusterList > DCClusterListContainer;

//_____________________________________________________________________________
class DCAnalyzer
{
public:
  DCAnalyzer(const RawData& raw_data);
  ~DCAnalyzer();
  enum e_type { k_BLC1a, k_BLC1b, k_BLC2a, k_BLC2b, k_SDC, k_BPC1=k_SDC, k_BPC, k_BPC2=k_BPC, k_FDC,  n_type };

private:
  DCAnalyzer( const DCAnalyzer& );
  DCAnalyzer& operator =( const DCAnalyzer& );

private:
  const RawData*              m_raw_data;
  std::vector<bool>           m_is_decoded;
  std::vector<int>            m_much_combi;
  std::vector<DCHitContainer> m_BLC1aHC;
  std::vector<DCHitContainer> m_BLC1bHC;
  std::vector<DCHitContainer> m_BLC2aHC;
  std::vector<DCHitContainer> m_BLC2bHC;
  std::vector<DCHitContainer> m_SDCHC;
  std::vector<DCHitContainer> m_BPCHC;
  std::vector<DCHitContainer> m_FDCHC;

  DCClusterListContainer m_DCCC;

  LocalTrackContainer m_BLC1aTC;
  LocalTrackContainer m_BLC1bTC;
  LocalTrackContainer m_BLC2aTC;
  LocalTrackContainer m_BLC2bTC;
  LocalTrackContainer m_BLC1TC;
  LocalTrackContainer m_BLC2TC;
  LocalTrackContainer m_SDCTC;
  LocalTrackContainer m_BPCTC;
  LocalTrackContainer m_FDCTC;

  inline int MakeKey(int cid,int xy){ return cid<<3 | xy; }

public:
  bool DecodeRawHits(double retiming_t0=0, double retiming_def=0 );
  bool DecodeRawHits(e_type k_type,const int &detid, double retiming=0 );
  bool DecodeDCHits(const int &detid);

  DCHitContainer& GetDCHC( const int &detid, int layer )
  { return const_cast<DCHitContainer&>(std::as_const(*this).GetDCHC(detid, layer)); }
  const DCHitContainer& GetDCHC( const int &detid, int layer ) const;

  inline DCClusterList& GetDCCL( const int &detid, const int &xy );
  inline LocalTrackContainer& GetTC( const int &detid ) ;
  inline int GetNClusters( const int &detid, const int &xy, const int &i );
  inline int GetNClusterContainers( const int &detid, const int &xy );
  inline const int GetNTracks( const int &detid );
  inline LocalTrack* GetTrack( const int &detid,int i );
  inline DCCluster* GetCluster( const int &detid, const int &xy, const int &i, const int &j );
  //  inline DCClusterList& GetClusterList( const int &detid, const int &xy) const;

  bool TrackSearchAll( );
  //  int TrackSearch( const int &cid, const int &debug=0 );
  int TrackSearch( const int &cid, const bool &TIMING, const int &debug=0 );
  bool MakePairs( const int &detid, const int &layer1, const int &layer2, const bool &isMC, const double &maxsub);
  bool MakePairsAll( const bool &isMC=false, const double &maxsub=999 );
  bool MakePairs( const int &detid, const bool &isMC=false, const double &maxsub=999  );
  int DeleteBadClusters( const int &detid, const int &xy, const int &ith);
  bool ReCalcDCHits( std::vector<DCHitContainer>& cont,
		     bool applyRecursively=false );
  bool ReCalcDCHits( bool applyRecursively=false );

  bool ReCalcAll();
  void Print(const int &detid);
protected:
  void ClearDCHits();
  void ClearDCHits( const int &detid );
  void ClearDCTracks();
  void ClearDCTracks( const int &detid );
  void ClearDCClusters();
};

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetDCHC( const int &detid,int layer ) const
{
  static DCHitContainer null_container;
  // std::cout<<"DCAnalyzer::GetDCHC() "<<detid<<"  "<<layer<<std::endl;
  if( layer>8 ) layer=0;
  switch(detid){
  case DetIdSDC:    return m_SDCHC.at(layer);
  case DetIdBLC1a:  return m_BLC1aHC.at(layer);
  case DetIdBLC1b:  return m_BLC1bHC.at(layer);
  case DetIdBLC2a:  return m_BLC2aHC.at(layer);
  case DetIdBLC2b:  return m_BLC2bHC.at(layer);
  case DetIdBPC:    return m_BPCHC.at(layer);
  case DetIdFDC:
    if( layer>5 ) layer=0;
    //    std::cout<<"DCAnalyzer::GetDCHC() "<<detid<<"  "<<layer<<std::endl;
    return m_FDCHC.at(layer);
  case DetIdVFT:
    // ignore VFT
    return null_container;
  default:
    std::cout << "E# invalid detector id "<< detid << std::endl;
    // gSystem->Exit(1);
    return null_container;
  }
}

inline const int
DCAnalyzer::GetNTracks( const int &detid )
{
  return GetTC(detid).size();
}

inline LocalTrack*
DCAnalyzer::GetTrack( const int &detid, int i )
{
  if(i< GetTC(detid).size() ) return GetTC(detid)[i];
  else return 0;
}
inline LocalTrackContainer&
DCAnalyzer::GetTC( const int &detid )
{
  static const std::string func_name(std::string("[DCAnalyzer::")+__func__+"()]");
  switch(detid){
  case DetIdSDC:    return m_SDCTC;
  case DetIdBLC1:   return m_BLC1TC;
  case DetIdBLC1a:  return m_BLC1aTC;
  case DetIdBLC1b:  return m_BLC1bTC;
  case DetIdBLC2:   return m_BLC2TC;
  case DetIdBLC2a:  return m_BLC2aTC;
  case DetIdBLC2b:  return m_BLC2bTC;
  case DetIdBPC:    return m_BPCTC;
  case DetIdFDC:    return m_FDCTC;
  default:
    std::cout<<"E# "<<func_name<<" invalid detector id "<< detid<<std::endl;
    exit(0);
  }
}
inline DCClusterList&
DCAnalyzer::GetDCCL( const int &detid, const int &xy )
{
  return m_DCCC[MakeKey(detid,xy)];
}
inline int
DCAnalyzer::GetNClusterContainers( const int &detid, const int &xy )
{
  return GetDCCL(detid,xy).size();
}
inline int
DCAnalyzer::GetNClusters( const int &detid, const int &xy, const int &i )
{
  if(GetNClusterContainers(detid,xy)>i){
    return m_DCCC[MakeKey(detid,xy)][i].size();
  }
  return 0;
}
inline DCCluster*
DCAnalyzer::GetCluster( const int &detid, const int &xy, const int &i, const int &j )
{
  if(GetNClusters(detid,xy,i)>j){
    return GetDCCL(detid,xy)[i][j];
  }
  return 0;
}
#endif
