// -*- C++ -*-

#ifndef BHT_HIT_HH
#define BHT_HIT_HH

#include <string>

#include <TObject.h>
#include <TString.h>

#include <std_ostream.hh>

#include "Hodo2Hit.hh"

//______________________________________________________________________________
class BHTHit : public Hodo2Hit
{
public:
  explicit  BHTHit( HodoRawHit *object, int index=0 );
  virtual  ~BHTHit( void );

private:
  BHTHit( const BHTHit&  );
  BHTHit& operator =( const BHTHit&);

public:
  virtual Bool_t Calculate( void ) override;
  virtual Bool_t ReCalc( Bool_t allpyRecursively=false )
  { return BHTHit::Calculate(); }

  //  ClassDef(BHTHit,0);
};

#endif
