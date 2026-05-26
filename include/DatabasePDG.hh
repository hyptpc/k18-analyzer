/**
 *  file: DatabasePDG.hh
 *  date: 2017.04.10
 *
 */

#ifndef DATABASE_PDG_HH
#define DATABASE_PDG_HH

#include <Rtypes.h>
#include <TPDGCode.h>

namespace pdg
{
  const Int_t kDeuteron = 1000010020; //deuteron
  const Int_t kTriton   = 1000010030; //tirton
  const Int_t kHe3      = 1000020030; //He3
  const Int_t kHe4      = 1000020040; //He4
  const Int_t kOther    = 9999;

  /// Returns particle mass [GeV/c^2] (ROOT TDatabasePDG; d,t,He3,He4 use fixed values).
  Double_t Mass(Int_t pdg_code);
  Double_t PionMass();     // [GeV/c^2]
  Double_t KaonMass();     // [GeV/c^2]
  Double_t ProtonMass();   // [GeV/c^2]
  Double_t NeutronMass();  // [GeV/c^2]
  Double_t ElectronMass(); // [GeV/c^2]
  Double_t DeuteronMass(); // [GeV/c^2]
  Double_t TritonMass();   // [GeV/c^2]
  Double_t LambdaMass();   // [GeV/c^2]
  Double_t SigmaNMass();   // [GeV/c^2]
  Double_t SigmaPMass();   // [GeV/c^2]
  Double_t XiMass();       // [GeV/c^2]
  void Print(Int_t pdg_code);
  void Print();
}

#endif
