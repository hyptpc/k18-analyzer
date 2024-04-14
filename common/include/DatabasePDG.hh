/**
 *  file: DatabasePDG.hh
 *  date: 2017.04.10
 *
 */

#ifndef DATABASE_PDG_HH
#define DATABASE_PDG_HH

#include <TPDGCode.h>

namespace pdg
{
  const int kDeuteron= 1000010020; //deuteron
  const int kTriton  = 1000010030; //tirton
  const int kHe3     = 1000020030; //He3
  const int kHe4     = 1000020040; //He4
  const int kOther   = 9999;
  // Mass [GeV/c2]
  double Mass( int pdg_code );
  double KaonMass( void );
  double PionMass( void );
  double ProtonMass( void );
  double NeutronMass( void );
  double LambdaMass( void );
  double SigmaNMass( void );
  double SigmaPMass( void );
  double XiMass( void );
  void   Print( int pdg_code );
  void   Print( void );
}

#endif
