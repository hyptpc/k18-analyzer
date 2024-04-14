#include "GeomTools.hh"
#include "MathTools.hh"
#include "DetectorID.hh"
#include "DetectorList.hh"
#include "CDCWireMapMan.hh"
#include "BLDCWireMapMan.hh"
#include "GeomMapMan.hh"

#include "TMath.h"

#define DEBUG 0
// #if !defined(R__ALPHA) && !defined(R__SOLARIS) && !defined(R__ACC) 
// NamespaceImp(GeomTools) 
// #endif 

namespace{
  const DetectorList&   dlist = DetectorList::GetInstance();
  const CDCWireMapMan&  gCDC  = CDCWireMapMan::GetInstance();
  const BLDCWireMapMan& gBLDC = BLDCWireMapMan::GetInstance();
  const GeomMapMan&     gGeom = GeomMapMan::GetInstance();
  enum gGeoParam { kPosX=0,
		   kPosR=0,
		   kPosY=1,
		   kPosPhi=1,	       
		   kPosZ=2,
		   kRotX=3,
		   kRotY=4,
		   kRotZ=5,
		   kSizeX=6,
		   kRmin=6,
		   kSizeY=7,
		   kRmax=7,
		   kSizeZ=8,
		   kdPhi=8,
		   kFlag=9,
		   kCylZ=9,
		   kMother=10
  };
  //  enum gHelixParam{ };
  const double mm=0.1; 
  const double m=1000.*mm; 
  const double cm=10.*mm; 
  const double um=1.e-3*mm; 
  const double degree=1.; 
  const double deg=1.; 
  const double unit=cm; // unit in parameter file
}
void geom::MakeGeometry()
{
  std::cout<<"============================="<<std::endl;
  //  std::cout<<"start Makeing TGeo"<<std::endl;
  new TGeoManager("K1.8BR spectrometer", "geometry for K1.8BR spectrometer");
  ConstructMaterial();
  ConstructHadronHall();
  ConstructHodoscopes();
  ConstructChambers();
  ConstructTarget();
  //--- close the geometry
  gGeoManager->CloseGeometry();
  gGeoManager->CheckOverlaps(0.00001, "d"); // unit??
  gGeoManager->PrintOverlaps();
  std::cout<<"============================="<<std::endl;
}

bool geom::IsSameVolume(const TVector3 &pos1, const TVector3 &pos2,double margin)
{
  TVector3 dir=(pos2-pos1).Unit();
  gGeoManager->InitTrack(pos1.X(),pos1.Y(),pos1.Z(),dir.X(),dir.Y(),dir.Z());
  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep();
  double length=gGeoManager->GetStep();
  //  std::cout<<length<<std::endl;
  if( (pos2-pos1).Mag()<(length+margin) ) return true;
  return false;
}

bool geom::IsSameVolumeHelix(const double param[5],const TVector3 &pos1, const TVector3 &pos2,double margin)
{
  //  std::cout<<"[geom::IsSameVolumeHelix] start!!!"<<std::endl;
  TVector3 tmppos;
  double tmpl;
  TString mat;
  int id;
  if(!geom::HelixStepToNextVolume(param,pos1,tmppos,tmpl,mat,id)) return false;  
  double tmpl2=math::CalcHelixArc(param,pos1,pos2);
  if(tmpl2<(tmpl+margin)) return true;
  return false;
}

bool geom::HelixStepToNextVolume(const double param[5],const TVector3 &in,TVector3 &pos2,double &length, TString &mat, int &id)
{
#if DEBUG
  std::cout<<"[geom::HelixStepToNextVolume] start!!!"<<std::endl;
#endif
  double defaultstep=1*cm;
  double step=defaultstep;  
  length=0;
  TVector3 pos1=in;
  pos2=math::CalcHelixStep(param,pos1,step);
  while(step>9*um){
#if DEBUG
    std::cout<<geom::GetMaterial(pos2)<<"  "<<step<<"  "<<length<<"  ";pos2.Print();
#endif
    if(geom::IsSameVolume(pos1,pos2)){
      length+=step;
      pos1=pos2;
      pos2=math::CalcHelixStep(param,pos1,step);
    }else{
      step/=10.;
      pos2=math::CalcHelixStep(param,pos1,step);
    }    
    if(length>100*cm){
      std::cout<<"[geom::HelixStepToNextVolume] too long"<<std::endl;
      in.Print();
      pos2.Print();
      std::cout<<geom::GetMaterial(in)<<"  "<<geom::GetMaterial(pos2)<<"  "<<step<<"  "<<length<<std::endl;
      return false;
    }
  }
  step*=10;
  pos2=math::CalcHelixStep(param,pos1,step);
  double tmpstep;
  geom::CrossBoundary(pos1,pos2,step,mat,tmpstep,id);
  id/=1000;
  pos2=math::CalcHelixStep(param,pos1,tmpstep+1*um);
  length+=tmpstep;
  //  std::cout<<length<<std::endl;
  return true;
}
bool geom::StepToNextVolume(const TVector3 &pos1,const TVector3 &pos2,double &length, TString &newmat)
{
  // almost the same as cross boundary. to be checked later.
  TVector3 dir=(pos2-pos1).Unit();
  gGeoManager->InitTrack(pos1.X(),pos1.Y(),pos1.Z(),dir.X(),dir.Y(),dir.Z());
  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep();
  length=gGeoManager->GetStep()+1*um;
  if(length>(pos2-pos1).Mag()) return false;
#if 0
  const double *pos=gGeoManager->GetCurrentPoint();
  for(int i=0;i<3;i++) std::cout<<pos1[i]<<"  ";
  for(int i=0;i<3;i++) std::cout<<pos2[i]<<"  ";
  std::cout<<length<<std::endl;
  //  std::cout<<std::endl;
#endif
  //  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep();
  if(!node){
#if DEBUG
    std::cout<<"newnode cannot be defined"<<std::endl;
#endif
    return false;
  }
  gGeoManager->CdNext(); // needed?

  if(node->GetMedium())
    newmat = node->GetMedium()->GetName();
  else
    std::cout<<"no medium in "<<node->GetNumber()<< "  "<<node->GetName()<<std::endl;
  if(newmat=="Concrete"||newmat=="Fe"){
    pos1.Print();
    pos2.Print();
    std::cout<<newmat<<" !!!"<<std::endl;
    return false;
  }else  if(newmat.IsNull()){
    std::cout<<node->GetNumber()<<"  "<<node->GetName()<<"  ";
    pos1.Print();
  }
  return true;
}

bool geom::CrossBoundary(const TVector3 &pos1,const TVector3 &pos2,const double &step, TString &mat, double &tmpstep, int &newid){
#if DEBUG
  for(int i=0;i<3;i++) std::cout<<pos1[i]<<"  ";
  for(int i=0;i<3;i++) std::cout<<pos2[i]<<"  ";
  std::cout<<step<<std::endl;
#endif
  TVector3 dir=(pos2-pos1).Unit();
  gGeoManager->InitTrack(pos1.X(),pos1.Y(),pos1.Z(),dir.X(),dir.Y(),dir.Z());
  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep(step,true);
  newid=node->GetNumber();
  mat = node->GetMedium()->GetName();
  tmpstep=gGeoManager->GetStep();
  if(mat.IsNull()){
    std::cout<<newid<<"  "<<node->GetName()<<"  ";
    pos1.Print();
  }
  return true;
}

double geom::CrossBoundary(const TVector3 &pos1,const TVector3 &pos2,const double &step){
  TVector3 dir=(pos2-pos1).Unit();
  gGeoManager->InitTrack(pos1.X(),pos1.Y(),pos1.Z(),dir.X(),dir.Y(),dir.Z());
  gGeoManager->FindNextBoundaryAndStep(step,true);
  return gGeoManager->GetStep();
}

double geom::CalcLengthinFiducial(const TVector3 &pos,const TVector3 &dir){
  TVector3 dir2=(dir-pos).Unit();
  TVector3 tmppos=pos;
  TString newmat;
  double length;
  bool FIDUCIAL=false;
  TVector3 in,out;
  while(geom::StepToNextVolume(tmppos,dir,length,newmat)){
    tmppos+=dir2*length;
    int id=GetID(tmppos);
    if(FIDUCIAL&&id!=DetIdFiducial){
      out=tmppos;
      break;
    }
    if(id==DetIdFiducial){
      in=tmppos;
      FIDUCIAL=true;
    }
    if(tmppos.Z()>10) return 0;
  }
  return (out-in).Mag();
}

void geom::ConstructHadronHall()
{
  std::cout<<"["<<__func__<<"]"<<std::endl;
  // We need to construct hadronhall first 
  ConstructShape(DetIdHall);
  ConstructShape(DetIdFloor);
  ConstructShape(DetIdDoraemon);
  //  ConstructMan();
}
void geom::ConstructTarget()
{
  std::cout<<"["<<__func__<<"]"<<std::endl;
  ConstructShape(DetIdCDCCFRP);
  ConstructShape(DetIdCDCMylar);

  ConstructShape(DetIdTarSys);
  ConstructShape(DetIdRadS);
  ConstructShape(DetIdTarCFRP);
  ConstructShape(DetIdTarCap);
  ConstructShape(DetIdTarCell);  
  ConstructShape(DetIdTarRing);
  //ConstructShape(DetIdTarget);
  ConstructShape(DetIdCellTube);
  ConstructShape(DetIdCellFlange);
  //ConstructShape(DetIdCellRing); // comment out 20150422
  // ConstructShape(DetIdBShield);
  // ConstructShape(DetIdBFrange);
  ConstructShape(DetIdFiducial);
}
void geom::ConstructHodoscopes(){
  std::cout<<"["<<__func__<<"]"<<std::endl;
  //  ConstructShape(DetIdBHD);
  ConstructShape(DetIdT0);
  ConstructShape(DetIdAC);
  ConstructShape(DetIdT0new);
  ConstructShape(DetIdDEF);
  //ConstructShape(DetIdVeto0);
#ifdef E73
  ConstructShape(DetIdVeto1);
  ConstructShape(DetIdPbF2);
  ConstructShape(DetIdBTC);
#endif
  ConstructShape(DetIdCDH);
#ifdef E15
  ConstructShape(DetIdIH);
#endif
}
void geom::ConstructChambers(){
  std::cout<<"["<<__func__<<"]"<<std::endl;
  ConstructShape( DetIdBLC2a );
  ConstructShape( DetIdBLC2b );
  ConstructShape( DetIdBPC   );  
  ConstructShape( DetIdCDC   );
}

bool geom::GetParam(const int &cid,const int &seg, double *param)
{
  bool status=false;
  if(cid==DetIdCDC)
    status=gCDC.GetGParam(param);
  else if(dlist.IsChamber(cid))
    status=gBLDC.GetGParam(cid,0,param);
  else
    status=gGeom.GetParam(cid,seg,param);
  return status;
}

TGeoVolume* geom::ConstructShape(Int_t CID,TGeoVolume *mother){
  //====================================================
  //--- define detector based on CID
  //====================================================
  Int_t nsegments=dlist.GetNsegs(CID);
  TString name=dlist.GetName(CID);
  std::cout<<__func__<<"  "<<CID<<"  "<<name<<std::endl;
  if(name=="NULL"){
    std::cout<<"No counter registerd for CID = "<<CID<<std::endl;
    return 0;
    exit(1);
  }
  //  std::cout<<nsegments<<"  "<<name<<"  "<<dlist.GetMaterial(CID)<<std::endl;
  TGeoMedium *medium1= gGeoManager->GetMedium("Air");
  TGeoMedium *medium2= gGeoManager->GetMedium(dlist.GetMaterial(CID).Data());
  if(nsegments==0) medium1=medium2;  

  double param[20];
  TGeoVolume* assembly=0;
  if(geom::GetParam(CID,-1,param)){
    assembly=geom::MakeShape(param,name,medium1);
#if 0
    std::cout<<CID<<"  "<<name<<"  "<<dlist.GetMaterial(CID).Data();
    for(int i=0;i<11;i++)
      std::cout<<"  "<<param[i];
    std::cout<<std::endl;
#endif
    if(param[kMother]==-1){
      gGeoManager->SetTopVolume(assembly); return assembly;
    }else{
      if(!mother)
	mother=gGeoManager->GetVolume(dlist.GetName((int)param[kMother]).Data());
      if(!mother)
	mother=gGeoManager->GetTopVolume();
      //    std::cout<<"mother= "<<CID<<"  "<<param[10]<<"  "<<mother->GetName()<<std::endl;
      if(mother&&assembly){
	TGeoCombiTrans* assembly_trans=geom::MakeTrans(param);
	assembly->SetLineColor(kBlack);
	assembly->SetTransparency(1);
	mother->AddNode(assembly, CID*1000, assembly_trans);
	if(nsegments==0) return assembly;
      }
    }
  }else{
    //  if(!assembly){
    std::cout<<"#W No assembly volume is defined."<<std::endl;
    std::cout<<"#W Each segment will be directly put in the mother volume."<<std::endl;
    if(!geom::GetParam(CID,0,param)) return 0;
    assembly=gGeoManager->GetVolume(dlist.GetName((int)param[kMother]).Data());
  }
  if(!assembly)
    return 0;
  
  TGeoVolume *segment;
  for (Int_t i=0; i<nsegments; i++){
    TString segname=Form("%s_segment%d",name.Data(),i);
    //    std::cout<<segname<<std::endl;
    if(!geom::GetParam(CID,i,param)) continue;
    segment = geom::MakeShape(param,segname,medium2); 
    TGeoCombiTrans *segment_trans = geom::MakeTrans(param);
    segment->SetLineColor(kBlue);
    segment->SetTransparency(1);
    assembly->AddNode(segment, CID*1000+i, segment_trans);
  }
  return assembly;
}

TGeoVolume* geom::MakeShape(double *param,const TString &name, TGeoMedium* medium ){
  Double_t z    = param[kCylZ]*unit/2.0;
  TGeoVolume *shape=0;
  if(z>0){
    Double_t rmin = param[kRmin]*unit;
    Double_t rmax = param[kRmax]*unit;
    Double_t phi  = param[kdPhi]*deg;
    //    std::cout<<"MakeTubs( "<<name<<" ) "<<rmin<<"  "<<rmax<<"  "<<z<<"  "<<phi<<std::endl;
    shape= gGeoManager->MakeTubs(name, medium, rmin, rmax, z, 0, phi);
  }else{
    Double_t x = param[kSizeX]*unit/2.0;
    Double_t y = param[kSizeY]*unit/2.0;
    z    = param[kSizeZ]*unit/2.0;
    if(z>0){
      //      std::cout<<"MakeBox( "<<name<<" ) "<<x<<"  "<<y<<"  "<<z<<std::endl;
      shape= gGeoManager->MakeBox(name, medium,x,y,z);
    }else{
      return 0;
    }
  }    
  return shape;
}

TGeoCombiTrans* geom::MakeTrans(double *param){
  Double_t pos_x = param[kPosX]*unit;
  Double_t pos_y = param[kPosY]*unit;
  Double_t pos_z = param[kPosZ]*unit;
  if(param[kFlag]<0&&param[kdPhi]<90){
    Double_t pos_r   = param[kPosR]*unit;
    Double_t pos_phi = param[kPosPhi]*TMath::DegToRad();
    Double_t pos_z   = param[kPosZ]*unit;
    TVector3 pos(pos_r,0,pos_z);
    pos.RotateZ(pos_phi);
    pos_x = pos.X();
    pos_y = pos.Y();
    pos_z = pos.Z();
  }
  TGeoRotation *rot   = new TGeoRotation();
  rot->RotateX(param[kRotX]);
  rot->RotateY(param[kRotY]);
  rot->RotateZ(param[kRotZ]);
  // std::cout<<"MakeTrans"<<pos_x<<"  "<<pos_y<<"  "<<pos_z<<"  "
  // 	   <<param[kRotX]<<"  "<<param[kRotY]<<"  "<<param[kRotZ]<<std::endl;
  TGeoCombiTrans *trans = new TGeoCombiTrans(pos_x,pos_y,pos_z,rot);
  return trans;
}

void geom::ConstructMaterial(){
  std::cout<<"["<<__func__<<"]"<<std::endl;
  //====================================================
  //--- define materials
  //====================================================
  TGeoMaterial *mat[30];
  for(int i=0;i<30;i++)
    mat[i]= new TGeoMaterial(Form("mat%d",i),    1.008,  1,  1.e-25);
  int i=0;
  new TGeoMedium("Vacuum"		, i, mat[i++]);
  new TGeoMedium("LHelium-3"		, i, mat[i++]);
  new TGeoMedium("LHelium-4"		, i, mat[i++]);
  new TGeoMedium("Beryllium"		, i, mat[i++]);
  new TGeoMedium("Air"			, i, mat[i++]);
  new TGeoMedium("CFRP"			, i, mat[i++]);
  new TGeoMedium("Mylar"		, i, mat[i++]);
  new TGeoMedium("Scinti"		, i, mat[i++]);
  new TGeoMedium("Plastic"		, i, mat[i++]);
  new TGeoMedium("Aluminum"		, i, mat[i++]);
  new TGeoMedium("Concrete"		, i, mat[i++]);
  new TGeoMedium("Iron"			, i, mat[i++]);
  new TGeoMedium("Tungsten"		, i, mat[i++]);
  new TGeoMedium("CDCGas"		, i, mat[i++]);
  new TGeoMedium("BLDCGas"		, i, mat[i++]);
  new TGeoMedium("Aerogel"		, i, mat[i++]);
  new TGeoMedium("AlBeMet"		, i, mat[i++]);
  new TGeoMedium("STAINLESS-STEEL"	, i, mat[i++]);
  new TGeoMedium("Graphite"		, i, mat[i++]);
  new TGeoMedium("Bi"			, i, mat[i++]);
  new TGeoMedium("LHydrogen"		, i, mat[i++]);
  new TGeoMedium("LDeuterium"		, i, mat[i++]);
  new TGeoMedium("G10"			, i, mat[i++]);
  new TGeoMedium("Li"			, i, mat[i++]);
  new TGeoMedium("PbF2"			, i, mat[i++]);
}

void geom::ConstructMan(){
  TGeoVolume *mother=gGeoManager->GetTopVolume();
  TGeoMedium* Air=gGeoManager->GetMedium("Air");
  //====================================================
  //--- define operator
  //====================================================
  Double_t beam_line_height = 2.0*m;
  Double_t operator_assembly_x = 60*cm/2.0;
  Double_t operator_assembly_y = 180*cm/2.0;
  Double_t operator_assembly_z = 30.0*cm/2.0;
  TGeoVolume *operator_assembly = gGeoManager->MakeBox("operator_assembly", 
						      Air, 
						      operator_assembly_x,
						      operator_assembly_y,
						      operator_assembly_z);
  Double_t operator_assembly_pos_x =  3.0*m;
  Double_t operator_assembly_pos_y = -beam_line_height + operator_assembly_y;
  Double_t operator_assembly_pos_z =  7.0*m;
  Double_t operator_assembly_angle =  90.0*degree; 
  TGeoRotation *operator_assembly_rot = new TGeoRotation();
  operator_assembly_rot->RotateY(operator_assembly_angle);
  TGeoCombiTrans *operator_assembly_trans = new TGeoCombiTrans(operator_assembly_pos_x,
							       operator_assembly_pos_y,
							       operator_assembly_pos_z,
							       operator_assembly_rot);
  operator_assembly->SetLineColor(kBlue);
  mother->AddNode(operator_assembly, 0, operator_assembly_trans);
  Double_t operator_head_rmin     = 0.0*cm;
  Double_t operator_head_rmax     = 15.0*cm;
  Double_t operator_head_thetamin = 0.0*degree;
  Double_t operator_head_thetamax = 180.0*degree;
  Double_t operator_head_phimin   = 0.0*degree;
  Double_t operator_head_phimax   = 360.0*degree;
  TGeoVolume *operator_head = gGeoManager->MakeSphere("operator_head", 
						     Air,
						     operator_head_rmin,
						     operator_head_rmax,
						     operator_head_thetamin,
						     operator_head_thetamax,
						     operator_head_phimin,
						     operator_head_phimax);
  Double_t operator_head_pos_x =  0.0*m;
  Double_t operator_head_pos_y =  0.75*m;
  Double_t operator_head_pos_z =  0.0*m;
  Double_t operator_head_angle =  0.0*degree; 
  TGeoRotation *operator_head_rot = new TGeoRotation();
  operator_head_rot->RotateY(operator_head_angle);
  TGeoCombiTrans *operator_head_trans = new TGeoCombiTrans(operator_head_pos_x,
							   operator_head_pos_y,
							   operator_head_pos_z,
							   operator_head_rot);
  operator_head->SetLineColor(kBlue);
  operator_assembly->AddNode(operator_head, 0, operator_head_trans);
  Double_t operator_body_x = 40.0*cm/2.0;
  Double_t operator_body_y = 50.0*cm/2.0;
  Double_t operator_body_z = 10.0*cm/2.0;
  TGeoVolume *operator_body = gGeoManager->MakeBox("operator_body", 
						  Air,
						  operator_body_x,
						  operator_body_y,
						  operator_body_z);
  Double_t operator_body_pos_x =  0.0*m;
  Double_t operator_body_pos_y =  0.35*m;
  Double_t operator_body_pos_z =  0.0*m;
  Double_t operator_body_angle =  0.0*degree; 
  TGeoRotation *operator_body_rot = new TGeoRotation();
  operator_body_rot->RotateY(operator_body_angle);
  TGeoCombiTrans *operator_body_trans = new TGeoCombiTrans(operator_body_pos_x,
							   operator_body_pos_y,
							   operator_body_pos_z,
							   operator_body_rot);
  operator_body->SetLineColor(kBlue);
  operator_assembly->AddNode(operator_body, 0, operator_body_trans);
  Double_t operator_arm_rmin = 0.0*cm;
  Double_t operator_arm_rmax = 5.0*cm;
  Double_t operator_arm_z    = 70.0*cm/2.0;
  TGeoVolume *operator_arm = gGeoManager->MakeTube("operator_arm", 
						  Air,
						  operator_arm_rmin,
						  operator_arm_rmax,
						  operator_arm_z);
  for(Int_t i=0; i<2; i++){
    Double_t operator_arm_pos_x;
    if( i==0 )
      operator_arm_pos_x =  20.0*cm + operator_arm_rmax;
    if( i==1 )
      operator_arm_pos_x = -20.0*cm - operator_arm_rmax;
    Double_t operator_arm_pos_y =  0.25*m;
    Double_t operator_arm_pos_z =  0.0*m;
    Double_t operator_arm_angle =  90.0*degree; 
    TGeoRotation *operator_arm_rot = new TGeoRotation();
    operator_arm_rot->RotateX(operator_arm_angle);
    TGeoCombiTrans *operator_arm_trans = new TGeoCombiTrans(operator_arm_pos_x,
							    operator_arm_pos_y,
							    operator_arm_pos_z,
							    operator_arm_rot);
    operator_arm->SetLineColor(kBlue);
    operator_assembly->AddNode(operator_arm, i, operator_arm_trans);
  }
  Double_t operator_leg_rmin = 0.0*cm;
  Double_t operator_leg_rmax = 10.0*cm;
  Double_t operator_leg_z    = 100.0*cm/2.0;
  TGeoVolume *operator_leg = gGeoManager->MakeTube("operator_leg", 
						  Air,
						  operator_leg_rmin,
						  operator_leg_rmax,
						  operator_leg_z);
  for(Int_t i=0; i<2; i++){
    Double_t operator_leg_pos_x;
    if( i==0 )
      operator_leg_pos_x =  operator_leg_rmax;
    if( i==1 )
      operator_leg_pos_x = -operator_leg_rmax;
    Double_t operator_leg_pos_y =  -40*cm;
    Double_t operator_leg_pos_z =   0.0*m;
    Double_t operator_leg_angle =   90.0*degree; 
    TGeoRotation *operator_leg_rot = new TGeoRotation();
    operator_leg_rot->RotateX(operator_leg_angle);
    TGeoCombiTrans *operator_leg_trans = new TGeoCombiTrans(operator_leg_pos_x,
							    operator_leg_pos_y,
							    operator_leg_pos_z,
							    operator_leg_rot);
    operator_leg->SetLineColor(kBlue);
    operator_assembly->AddNode(operator_leg, i, operator_leg_trans);
  }
}

