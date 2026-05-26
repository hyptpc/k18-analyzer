/**
 *  file: DatabasePDG.cc
 *  date: 2017.04.10
 *  note: PDG code is defined in ROOT/include/TPDGCode.h
 *        Mass unit [GeV/c2]
 *
 */

#include "DatabasePDG.hh"

#include <iostream>
#include <string>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

namespace
{
  const std::string& name("DatabasePDG");
}

//______________________________________________________________________________
namespace pdg
{
  //______________________________________________________________________________
  Double_t
  Mass(Int_t pdg_code)
  {
    if( pdg_code==kDeuteron ) return 1.875613;
    if( pdg_code==kTriton   ) return 2.808921;
    if( pdg_code==kHe3      ) return 2.808391;
    if( pdg_code==kHe4      ) return 3.727379;
    TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
    return ( particle ? particle->Mass() : -1. );
  }

  //______________________________________________________________________________
  Double_t
  KaonMass()
  {
    return Mass(kKMinus);
  }

  //______________________________________________________________________________
  Double_t
  PionMass()
  {
    return Mass(kPiMinus);
  }

  //______________________________________________________________________________
  Double_t
  ProtonMass()
  {
    return Mass(kProton);
  }

  //______________________________________________________________________________
  Double_t
  NeutronMass()
  {
    return Mass(kNeutron);
  }

  //______________________________________________________________________________
  Double_t
  ElectronMass()
  {
    return Mass(kElectron);
  }

  //______________________________________________________________________________
  Double_t
  DeuteronMass()
  {
    return Mass(kDeuteron);
  }

  //______________________________________________________________________________
  Double_t
  TritonMass()
  {
    return Mass(kTriton);
  }

  //______________________________________________________________________________
  Double_t
  LambdaMass()
  {
    return Mass(kLambda0);
  }

  //______________________________________________________________________________
  Double_t
  SigmaNMass()
  {
    return Mass(kSigmaMinus);
  }

  //______________________________________________________________________________
  Double_t
  SigmaPMass()
  {
    return Mass(kSigmaPlus);
  }

  //______________________________________________________________________________
  Double_t
  XiMass()
  {
    return Mass(kXiMinus);
  }

  //______________________________________________________________________________
  void
  Print(Int_t pdg_code)
  {
    TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg_code);
    if( particle ) particle->Print();
  }

  //______________________________________________________________________________
  void
  Print()
  {
    TDatabasePDG::Instance()->Print();
  }

}
