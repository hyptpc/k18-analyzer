// ====================================================================
//    ComCrossSectionTable.hh
//
//
// ====================================================================
#ifndef ComCrossSectionTable_h
#define ComCrossSectionTable_h 1

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip> 
#include <vector>
#include <list>

#include "TObject.h"

//-----------------------------------------------//
// class CorssSection
//-----------------------------------------------//
class CrossSection : public TObject
{
private:
  int    fId;       // reaction ID, based on CERN-HERA-83-02
  double fCs;       // cross-section [mb]
  double fCsFactor; // factor for input CS
  double fInitMom;  // momentum of the beam particle, fInitPdg[0]
  std::vector <int>    fInitPdg; // pdg codes of initial particles 0:beam, 1:target, 2:reaction
  std::vector <int>    fFinlPdg; // pdg codes of finlhters
  std::vector <int>    fSpecPdg; // pdg codes of spectators
  std::vector <double> fPol;     // coefficients A(n) for Legendre polynomial;
  double fPolMax; // maximum value for random-number generation

public:
  CrossSection(){} // CAUTION :: not initialze here
  ~CrossSection(){}

  void Init(double initMom);
  void SetId(int id) { fId = id; }
  void SetCs(double cs) { fCs = cs; }
  void SetCsFactor(double csfactor) { fCsFactor = csfactor; }
  void SetInitMom( double mom ) { fInitMom = mom; }
  void SetInitPdg( const std::vector <int>& pdg ) { fInitPdg = pdg; }
  void SetFinlPdg( const std::vector <int>& pdg ) { fFinlPdg = pdg; }
  void SetSpecPdg( const std::vector <int>& pdg ) { fSpecPdg = pdg; }
  void SetPol( const std::vector <double>& pol ) { fPol = pol; }
  void SetPolMax(double val) { fPolMax = val; }

  int    Id() const { return fId; }
  double Cs() const { return fCs; }
  double CsFactor() const { return fCsFactor; }
  double InitMom() const { return fInitMom; }
  int InitPdgSize() const { return fInitPdg.size(); }
  int FinlPdgSize() const { return fFinlPdg.size(); }
  int SpecPdgSize() const { return fSpecPdg.size(); }
  const std::vector <int>& InitPdg() { return fInitPdg; }
  const std::vector <int>& FinlPdg() { return fFinlPdg; }
  const std::vector <int>& SpecPdg() { return fSpecPdg; }
  int    InitPdg(int i) const { return fInitPdg[i]; }
  int    FinlPdg(int i) const { return fFinlPdg[i]; }
  int    SpecPdg(int i) const { return fSpecPdg[i]; }
  int PolSize() const { return fPol.size(); }
  const std::vector <double>& Pol() { return fPol; }
  double Pol(int i) const { return fPol[i]; }
  double PolMax() const { return fPolMax; }

  void PrintCS();

  ClassDef(CrossSection,1)
};

inline void CrossSection::Init(double initMom)
{
  fId = 0; fCs = 0; fCsFactor = 0;
  fInitMom = initMom;
  fInitPdg.clear();
  fFinlPdg.clear();
  fSpecPdg.clear();
  fPol.clear();
  fPolMax = 0;
}



//-----------------------------------------------//
// class CorssSectionTable
//-----------------------------------------------//
class CrossSectionTable : public TObject
{
private:
  std::vector <CrossSection> fCS;
  double fTotalCS;

public:
  CrossSectionTable(){}
  CrossSectionTable(const char *CSFileName, double initMom);
  ~CrossSectionTable(){}

  void SetCS(const CrossSection& val) { fCS.push_back(val); }

  int CSSize() const { return fCS.size(); }
  const std::vector <CrossSection>& CS() { return fCS; }
  const CrossSection& CS(int i) { return fCS[i]; }
  double TotalCS() const { return fTotalCS; }

  void Init();
  double CheckCSList();
  double CalcTotalCS();
  void PrintAllCS();

  ClassDef(CrossSectionTable,1)
};

inline void CrossSectionTable::Init()
{
  fTotalCS = 0;
  fCS.clear();
}

#endif
