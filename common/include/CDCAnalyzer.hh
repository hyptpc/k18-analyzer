/**
 *  file: CDCAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef CDC_ANALYZER_HH
#define CDC_ANALYZER_HH

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>
#include <map>
#include <string>

#include "Event.hh"
#include "DCHit.hh"
class DCCluster;
class CDSTrack;
class RawData;

typedef std::vector<DCHit*>        DCHitContainer;
typedef std::vector<DCCluster*>    DCClusterContainer;
typedef std::vector<DCClusterContainer> DCClusterList;
typedef std::vector<CDSTrack*>   CDSTrackContainer;
//______________________________________________________________________________
class CDCAnalyzer
{
public:
  CDCAnalyzer( void );
  ~CDCAnalyzer( void );
  enum e_type { k_CDC, n_type };
  
private:
  CDCAnalyzer( const CDCAnalyzer& );
  CDCAnalyzer& operator =( const CDCAnalyzer& );
  
private:
  std::vector<bool>     m_is_decoded;
  std::vector<int>      m_much_combi;
  std::vector<DCHitContainer>       m_CDCHC;

  DCClusterList m_DCCC;
  CDSTrackContainer m_CDCTC;

public:
  bool DecodeRawHits( RawData* rawData, double retiming=0 );
  bool DecodeRawHits( RawData* rawData, e_type k_type,const int &detid, double retiming=0 );

  bool TrackSearch( );
  bool MakeClusters( const int &slayer, const double &maxsub);
  bool FindCircleTrackCandidates(const bool &LINE=false);
  bool CircleFitting();
  bool LineAxialFitting();
  bool FindStereoHits();
  bool HelixFitting();
  bool LineFitting();

  void SelectSharedHit();
  void ClearDCHits( const int &detid );
  void ClearDCHits( void );
  inline DCHitContainer& GetDCHC( const int &detid, int layer );
  inline DCHitContainer& GetDCHC( const int &layer ) { return GetDCHC(DetIdCDC,layer); }
  inline DCHit* GetDCHit(const int &layer, const int &wire);
  inline int GetNDCHit(const int &layer){ return GetDCHC(layer).size(); }

  inline CDSTrackContainer& GetTC(  ) { return m_CDCTC; }
  inline int GetNClusters( const int &slayer);
  inline const int GetNTracks( ) { return m_CDCTC.size(); }
  inline CDSTrack* GetTrack( int i ) { return m_CDCTC[i]; }
  inline DCCluster* GetCluster( const int &slayer, const int &i );
  void ClearDCTracks( void );
  void ClearDCClusters( void );

  void SetCDCTracksFromTree( CDCTrackContainer* cont );

  void Print(const int &detid);
};

//______________________________________________________________________________
inline DCHit*
CDCAnalyzer::GetDCHit( const int &layer,const int &wire )
{
  const DCHitContainer &cont = GetDCHC(layer);
  int mul=cont.size();
  for(int ihit=0;ihit<mul;ihit++){
    DCHit* hit=cont[ihit];
    if(hit->GetWire()==wire) return hit;
  }
  return 0;
}
//______________________________________________________________________________
inline DCHitContainer&
CDCAnalyzer::GetDCHC(const int &detid,int layer )
{
  const std::string& class_name("CDCAnalyzer");
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( layer>15 ) layer=0;
  switch(detid){
  case DetIdCDC:
    return m_CDCHC[layer];
  default:
    std::cout<<"E# "<<func_name<<" invalid detector id "<< detid<<std::endl;
    exit(0);
  }
}
inline int
CDCAnalyzer::GetNClusters( const int &slayer) 
{
  if(7>slayer){
    return m_DCCC[slayer].size();
  }
  return 0;
}
inline DCCluster*
CDCAnalyzer::GetCluster( const int &slayer, const int &i) 
{
  if(GetNClusters(slayer)>i){
    return m_DCCC[slayer][i];
  }
  return 0;
}

#endif
