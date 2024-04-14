// TransferMatrixMan.h

#ifndef TransferMatrixMan_h
#define TransferMatrixMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TObject.h"
#include "TMatrixD.h"
#include "TString.h"

class TransferMatrixMan
{
 public:
  static TransferMatrixMan& GetInstance(void);
  static const std::string& ClassName( void );
  ~TransferMatrixMan();

 private:
   TransferMatrixMan();
   TransferMatrixMan( const TString & filename );
   TransferMatrixMan( const TransferMatrixMan &right );

 public:
   void SetFileName( const TString & filename ){  FileName = filename; }
   bool Initialize();
   bool Initialize( const TString & filename );

 private:
  TString FileName;
  double CentralMomentum;
  double BLC1VIele[36];
  TMatrixD BLC1VIMatrix;
  TMatrixD BLC2VIMatrix;
  TMatrixD D5Matrix;
  double D5Matrix1st[36];
  double D5Matrix2nd[6][36];
  TMatrixD BLC2BLC1Matrix1st;
  TMatrixD BLC2BLC1Matrix2nd[6];
  bool m_isready;

 public:
  bool IsReady( void ) const { return m_isready; }
  TString GetFileName() { return FileName; }
  double GetCentralMomentum() const { return CentralMomentum; } 
  const double* GetBLC1VIMatrix() const { return BLC1VIele; }
  const double* GetD5Matrix() const { return D5Matrix1st; }
  const double* GetD5Matrix2nd(const int &i) const { return (i<6) ? D5Matrix2nd[i] : 0 ; }
  void CalcBLC1toBLC2(double *parblc1, double *parblc2,const int &order=1);
  void CalcBLC2toBLC1(double *parblc2, double *parblc1, const int &order=1);
  void CalcBLC1toVI(double *parblc1, double *parvi,const int &order=1);
  void CalcBLC2toVI(double *parblc2, double *parvi,const int &order=1);
  void Clear();
};
//______________________________________________________________________________
inline TransferMatrixMan&
TransferMatrixMan::GetInstance( void )
{
  static TransferMatrixMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
TransferMatrixMan::ClassName( void )
{
  static std::string g_name("TransferMatrixMan");
  return g_name;
}

#endif
