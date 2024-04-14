#ifndef GeomTools_hh
#define GeomTools_hh 1

#include <string>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TPolyLine.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoPgon.h"

namespace geom
{
  void MakeGeometry();
  void ConstructHadronHall();
  void ConstructTarget();
  void ConstructHodoscopes(); 
  void ConstructChambers();
  void ConstructMan();
  void ConstructMaterial();

  bool   CrossBoundary(const TVector3 &pos1,const TVector3 &pos2,const double &step, TString &mat, double &tmpstep, int &newid);
  double CrossBoundary(const TVector3 &pos1,const TVector3 &pos2,const double &step);
  bool GetParam(const int &cid,const int &seg, double *param);

  inline int GetID(const TVector3 &pos);
  inline int GetID2(const TVector3 &pos);
  inline int GetIDMat(const TVector3 &pos,TString &mat);
  inline TString GetMaterial(const TVector3 &pos);
  inline TString GetVolName(const TVector3 &pos);

  bool StepToNextVolume(const TVector3 &in, const TVector3 &out,double &length, TString &newmat);
  bool HelixStepToNextVolume(const double param[5],const TVector3 &pos1,TVector3 &pos2, 
			     double &length, TString &mat, int &id);
  bool IsSameVolume(const TVector3 &pos1, const TVector3 &pos2,double margin=0.);
  bool IsSameVolumeHelix(const double param[5],const TVector3 &pos1, const TVector3 &pos2,double margin=0.);
  double CalcLengthinFiducial(const TVector3 &pos1, const TVector3 &dir2);

  TGeoVolume* ConstructShape(Int_t CID,TGeoVolume* mother=0);
  TGeoVolume* MakeShape(double *param, const TString &name, TGeoMedium* medium);
  TGeoCombiTrans *MakeTrans( double *param);
};

inline int geom::GetID(const TVector3 &pos)
{
  //  if(pos.Perp()>200) return 0;
  TGeoNode *node=gGeoManager->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return 0;
  int id=node->GetNumber();
  id/=1000;
  return id;
}

inline int geom::GetID2(const TVector3 &pos)
{
  //  if(pos.Perp()>200) return 0;
  TGeoNode *node=gGeoManager->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return 0;
  int id=node->GetNumber();
  return id;
}

inline int geom::GetIDMat(const TVector3 &pos,TString &mat)
{
  //  if(pos.Mag()>100) return 0;
  TGeoNode *node=gGeoManager->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return 0;
  int id=node->GetNumber();
  id/=1000;
  mat = node->GetMedium()->GetName();
  if(mat.IsNull()){
    std::cout<<node->GetNumber()<<"  "<<node->GetName()<<"  ";
    pos.Print();
  }
  return id;
}

inline TString geom::GetMaterial(const TVector3 &pos)
{
  //  if(pos.Mag()>100) return 0;
  TGeoNode *node=gGeoManager->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return "None";
  return node->GetMedium()->GetName();
}

inline TString geom::GetVolName(const TVector3 &pos)
{
  //  if(pos.Mag()>100) return 0;
  TGeoNode *node=gGeoManager->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return "None";
  return node->GetName();
}
#endif
