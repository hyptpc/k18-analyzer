// -*- C++ -*-

#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include <TString.h>
#include <TVector3.h>

//_____________________________________________________________________________
namespace Kinematics
{
const TString& ClassName();
Double_t MassSquare(Double_t p, Double_t path, Double_t time);
Double_t CalcTimeOfFlight(Double_t p, Double_t path, Double_t mass);
TVector3 VertexPoint(const TVector3& Xin, const TVector3& Xout,
                     const TVector3& Pin, const TVector3& Pout);
TVector3 VertexPointByHonly(const TVector3& Xin,
                            const TVector3& Xout,
                            const TVector3& Pin,
                            const TVector3& Pout);
TVector3 VertexPoint(const TVector3& Xin, const TVector3& Xout,
                     const TVector3& Pin, const TVector3& Pout,
                     Double_t& dist);
TVector3 VertexPoint3D(const TVector3& Xin, const TVector3& Xout,
                       const TVector3& Pin, const TVector3& Pout,
                       Double_t& dist);
TVector3 VertexPointReal(const TVector3& Xin, const TVector3& Xout,
                         const TVector3& Pin, const TVector3& Pout,
                         Double_t& dist);
TVector3 VertexPointTF2(const TVector3& Xin, const TVector3& Xout,
                        const TVector3& Pin, const TVector3& Pout,
                        Double_t& dist);
TVector3 VertexPoint_Helix(const Double_t par1[5], const Double_t par2[5],
                           Double_t& dist, Double_t& t1, Double_t& t2);
Double_t CloseDist(const TVector3& Xin, const TVector3& Xout,
                   const TVector3& Pin, const TVector3& Pout);
TVector3 CorrElossIn(const TVector3& Pin,
                     const TVector3& Xin,
                     const TVector3& vtx, Double_t mass);
Double_t CalcLengthBeam(const TVector3& Pin,
                        const TVector3& Xin,
                        const TVector3& vtx);
TVector3 CorrElossOut(const TVector3& Pout,
                      const TVector3& Xout,
                      const TVector3& vtx, Double_t mass);
Double_t CalcLengthScat(const TVector3& Pout,
                        const TVector3& Xout,
                        const TVector3& vtx);
TVector3 CorrElossOutCheck(const TVector3& Pout,
                           const TVector3& Xout,
                           const TVector3& vtx, Double_t mass);
Bool_t   IsInsideTarget(const TVector3&point);
Bool_t   CalcCrossingPoint(Double_t u, Double_t v, const TVector3& Point,
                           Double_t *z1, Double_t *z2);
Double_t DiffE(Double_t mass, Double_t E, Double_t length, Double_t Elast);
Int_t    CalcDe(Double_t momentum, Double_t mass, Double_t distance,
                Double_t *momentum_cor, Double_t *energy_cor);
Double_t CalcDedx(Double_t beta);
Double_t Gamma(Double_t beta);
Double_t Beta(Double_t energy, Double_t mormentum);
}

//_____________________________________________________________________________
inline const TString&
Kinematics::ClassName()
{
  static TString s_name("Kinematic");
  return s_name;
}

#endif
