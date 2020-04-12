/**
 *  file: TPCCluster.hh
 *  date: 2020.04.11
 *
 */

#ifndef TPC_CLUSTER_HH
#define TPC_CLUSTER_HH

#include <cstddef>
#include <string>
#include <vector>

#include <TVector3.h>

#include "TPCHit.hh"

typedef std::vector<TPCHit*>       TPCHitContainer;

//______________________________________________________________________________
class TPCCluster
{
public:
  TPCCluster( int layer, TPCHitContainer& HitCont );
  ~TPCCluster( void );

private:
  int      m_layer_id;
  double   m_charge;
  TVector3 m_pos;
  std::vector<TPCHit*> m_tpchits;
  bool m_pos_calculated;

  void CalculateWeightedMean( void );

public:
  int  LayerId( void )		const { return m_layer_id; }
  double Charge( void )		const { return m_charge; }
  
  void Print( const std::string& arg="" ) const;
  
  int  GetClusterSize()		const { return m_tpchits.size(); }
  void AddTPCHit(TPCHit* hit);
  TVector3  Position( void );
  double X( void );
  double Y( void );
  double Z( void );
  std::vector<TPCHit*> GetTPCHits() const { return m_tpchits; }
  void ClearTPCHits();

};

#endif
