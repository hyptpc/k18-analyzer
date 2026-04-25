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
  // Mass [GeV/c2]
  Double_t Mass(Int_t pdg_code);
  Double_t KaonMass();
  Double_t PionMass();
  Double_t ProtonMass();
  Double_t NeutronMass();
  Double_t LambdaMass();
  Double_t SigmaNMass();
  Double_t SigmaPMass();
  Double_t XiMass();
  void Print(Int_t pdg_code);
  void Print();
}

#endif
