//  Authors: Wooseung Jung

#ifndef HYPTPCHit_HH
#define HYPTPCHit_HH

//ROOT
#include <TVector3.h>
#include <TObject.h>

//k18-analyzer
#include "TPCLTrackHit.hh"

class HypTPCHit : public TObject
{
public:

  HypTPCHit(TPCLTrackHit& right);
  ~HypTPCHit();

  //return the TPCLTrackit
  TPCLTrackHit &GetHit() const;

private:

  TPCLTrackHit &m_hit;

  ClassDef(HypTPCHit, 1)
}; //class HypTPCHit

#endif // HypTPCHit_hh
