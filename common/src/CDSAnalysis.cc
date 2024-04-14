// CDSAnalysis.cpp
#include "CDSAnalysis.hh"
#include "TF1.h"
#include "DatabasePDG.hh"
#include "GeomTools.hh"
#include "MathTools.hh"
#include "ELossTools.hh"
#include "GeomMapMan.hh"
#include "UserParamMan.hh"
#define DEBUG 0

namespace{
  const GeomMapMan&     gGeom = GeomMapMan::GetInstance();
  const UserParamMan&   gUser = UserParamMan::GetInstance();
  TVector3 DEFVECT(-999,-999,-999);
  double degree=TMath::DegToRad();
  const double mm=0.1;
  const double cm=10*mm;
  double cdc_rin=154.*mm; 
  double cdh_rin=544.*mm; 
  double cdh_rout=574.*mm;
  double cdh_r=cdh_rin;
  //  double cdh_r=(cdh_rin+cdh_rout)*0.5; // from 20221122 to 20230113; does not fit with MC analysis
  double cdh_dphi=7.*degree;
  double mybeta(double mass, double mom){
    return fabs(mom)/sqrt(mass*mass+mom*mom);
  }
  double mymass2(double mom, double beta){
    return mom*mom*(1/(beta*beta)-1);
  }
}
double cds::Mass2Resol(int pid, double mom, double fl)
{
  double tofresol=0.160;//ns
  double lv=TMath::C()*1e-6*mm;//mm/ns  
  double mass=pdg::Mass(pid);
  double momres=MomResol(pid,mom);
  double tmp2=4*pow(mass,4)*pow(momres,2)+4*pow(mom,2)*(pow(mom,2)+pow(mass,2))*pow(lv/fl*tofresol,2);
  if(tmp2>0)
    return sqrt(tmp2);
  return 999.;
}

double cds::MomResol(int pid, double mom)
{
  int id=-1;
  if(pid==kPiPlus||pid==kPiMinus) id=0;
  else if(pid==kProton) id=1;
  else if(pid==kKMinus)  id=2;
  else return 999;
  // for mom resol. [0]+[1]*x+[2]/(x-[3])
#if 0
  // parameter for de corrected
  double param[3][4]={{1.03626e-05, 0.062906, 0.000538243, 0.0421856},
		      {-5.6182e-05, 0.0632374, 0.0002, 0.195},
		      {0.00216906, 0.0596655, 0.000187301, 0.13}
  };
#else  
  //parameter at CDC 
  double param[3][4]={{-0.000154926, 0.0628548, 0.000544216, 0.0203668}, 
		      {-0.0119525, 0.0699056, 0.0058263, 0.00385919}, 
		      {-0.00514371, 0.0653424, 0.00293869, -0.0220217}
  };
#endif
  double tmp=param[id][0]+param[id][1]*mom+param[id][2]/(mom-param[id][3]);
  if(tmp>0)
    return tmp;
  return 999;
}

bool cds::FindMass2(CDSTrack *track, LocalTrack* bpc,double tof, double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof)
{
  // --- INPUT
  // CDSTrack track: cdc track
  // LocalTrack bpc: beam track
  // double tof: time of flight between t0 and cdh
  // double beammom: beam momentum at D5
  // --- OUTPUT
  // double beta_calc: calculated beta of CDC track
  // double mass2: CDCtrack
  // double tmptof: time of flight from vertex to CDH?

  TVector3 vtxb1,vtxb2;
  double dis;
  cds::CalcLineHelixVertex(bpc,track,vtxb1,vtxb2,dis);
  TVector3 t0pos=bpc->GetPosatZ(-110*cm-0.5*cm);
  bool tmp=cds::FindMass2(track,t0pos,vtxb1,vtxb2,tof,beammom,pdg::Mass(pid_beam),beta_calc,mass2,tmptof);
  return tmp;
}

bool cds::FindMass2(CDSTrack *track, pBeam *beam,double tof ,double &beta_calc,double &mass2,double &tmptof)
{
  TVector3 vtxb1,vtxb2;
  double dis;
  cds::CalcLineHelixVertex(beam,track,vtxb1,vtxb2,dis);
  TVector3 t0pos=beam->t0pos();
  double beammom=beam->mom();
  int pid_beam=beam->pid();
  bool tmp=cds::FindMass2(track,t0pos,vtxb1,vtxb2,tof,beammom,pdg::Mass(pid_beam),beta_calc,mass2,tmptof);
  return tmp;
}

bool cds::FindMass2(CDSTrack *track, TVector3 t0pos,TVector3 vtxbpc,TVector3 vtxcdc,double tof_t0_cdh, double beammom, double beammass,double &beta_calc,double &mass2,double &tof_vtx_cdcin)
{
  double param[5];
  TVector3 tmp;
  track->GetParameters(DetIdCDC,param,tmp);   
  TVector3 cdhvtx=math::CalcHelixPosatR(param,cdh_r);
  TVector3 pos_cdcin=math::CalcHelixPosatR(param,cdc_rin);
  if(fabs(cdhvtx.Z())>1000 || fabs(pos_cdcin.Z())>1000 ){
    return false;
  }
  double cdc_dis=math::CalcHelixArc(param,cdhvtx,pos_cdcin);
  
  double mom,time_beam;
  eloss::CalcElossBeam(t0pos,vtxbpc,beammom,beammass,mom,time_beam);
  double tof_vtx_cdh=tof_t0_cdh-time_beam;

  mom = track->mom();
  int count =0;
  double tmpl,tmpmass2;
  do{
    tmpmass2=mass2;
    double tmpmass=sqrt(tmpmass2);
    if(tmpmass2<1e-8) tmpmass=1e-4; // 100keV
    if( !eloss::CalcHelixElossToVertex(param,vtxcdc,mom,tmpmass,tmpl,tof_vtx_cdcin,cdc_rin)){
      std::cout<<"failed in eloss::calchelixelosstovertex()"<<std::endl;
      // track->Print();
      //    pos_cdcin.Print();
      // vtxcdc.Print();
      // tmp.Print();
      // for(int i=0;i<5;i++) std::cout<<param[i]<<"  ";
      // std::cout<<std::endl;
      // exit(0);
      return false;
    }
    beta_calc=(cdc_dis)/(tof_vtx_cdh-tof_vtx_cdcin)/(TMath::C()*1e-6*mm); //mm/ns
    mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
    count++;
    if(count>50)   return false;
  } while( TMath::Abs(tmpmass2-mass2)>0.0001);
  return true;
}

bool cds::CalcVertex2HelixdE(CDSTrack *cds1, CDSTrack *cds2, const double &mass1, const double &mass2, TVector3 &vtx1, TVector3 &vtx2)
{
  // calc vertex of 2 cdc tracks considering curavature change by dE loss 
  //  std::cout<<"[ "<<__FUNCTION__<<" ]"<<std::endl;
  double par1[5];
  double par2[5];
  TVector3 tmppos1,tmppos2;
//  double margin=10.*mm;
  double dis, tmpdis=999;
  double tmp1,tmp2;
  cds1->GetParameters(DetIdCDC,par1,tmppos1,tmp1,tmp2);
  cds2->GetParameters(DetIdCDC,par2,tmppos2,tmp1,tmp2);
  math::HelixToHelix(par1,par2,vtx1,vtx2,dis);
  if(dis>100*cm) {
#if DEBUG
    std::cout<<"#W too much distance between the two helix"<<std::endl;
    vtx1.Print();
    vtx2.Print();
#endif
    return false;
  }
  int id1=geom::GetID(vtx1);
  int id2=geom::GetID(vtx2);
  int id1_prev=-1;
  int id2_prev=-1;
  int count=0;
  while(id1!=id1_prev||id2!=id2_prev){
    id1_prev=id1;
    id2_prev=id2;
    if(!cds1->GetParameters(id1,par1,tmppos1,tmp1,tmp2)){
      cds1->CalcEnergyLoss(vtx1,mass1,par1,tmppos1);
    }
    if(!cds2->GetParameters(id2,par2,tmppos2,tmp1,tmp2)){
      cds2->CalcEnergyLoss(vtx2,mass2,par2,tmppos2);
    }
    math::HelixToHelix(par1,par2,vtx1,vtx2,dis);
    id1=geom::GetID(vtx1);
    id2=geom::GetID(vtx2);
    count++;
    if(count>10){
      std::cout<<"#E too many loops. something is wrong."<<std::endl;
      return false;
    }
  }
  return true;
}

bool cds::CalcVertex2Helix2(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2)
{
  // copy from knana
  double par1[5];
  double par2[5];
  TVector3 tmppos1,tmppos2;
  TVector3 tmp1,tmp2;
  int id,id2;
  double margin=1.*mm;
  double dis, tmpdis=999;
  for(int i=0;i<cds1->nParamSets();i++){
    cds1->GetNthParameters(i,id,par1,tmppos1);
    for(int i2=0;i2<cds2->nParamSets();i2++){
      cds2->GetNthParameters(i2,id2,par2,tmppos2);
      math::HelixToHelix(par1,par2,tmp1,tmp2,dis);
      //      if((GeomTools::GetID(tmp1)==id||GeomTools::GetID(tmp2)==id2)&&dis<tmpdis){
      //      std::cout<<GeomTools::GetID(tmp1)<<"  "<<id<<"   "<<tmp1.Perp()<<"  "<<tmp1.Z()<<std::endl;
      if((geom::GetID(tmp1)==id||geom::IsSameVolumeHelix(par1,tmppos1,tmp1,margin))){	
	//	std::cout<<dis<<"  "<<tmpdis<<std::endl;
	//	std::cout<<GeomTools::GetID(tmp2)<<"  "<<id2<<"   "<<tmp2.Perp()<<"  "<<tmp2.Z()<<std::endl;
	if((geom::GetID(tmp2)==id2||geom::IsSameVolumeHelix(par2,tmppos2,tmp2,margin))
	   &&dis<tmpdis){
	  vtx1=tmp1;
	  vtx2=tmp2;
	  tmpdis=dis;
	  //	  std::cout<<"!!!"<<std::endl;
	}
      }
    }
  }
  //  std::cout<<"-----------!!!"<<std::endl;
  if(tmpdis<100*cm) return true;
  else{
#if DEBUG
    std::cout<<"++++++++++++++++++++++++++++++++"<<tmpdis<<std::endl;
    cds1->Print();
    cds2->Print();
    exit(0);
#endif
    return false;
  }
}

bool cds::CalcVertex2Helix(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2)
{
  // simple vertexing of two helices without energy loss
  double par1[5];
  double par2[5];
  TVector3 tmpvtx;
  double dis,tmp1,tmp2;
  cds1->GetParameters(DetIdCDC,par1,tmpvtx,tmp1,tmp2);;
  cds2->GetParameters(DetIdCDC,par2,tmpvtx,tmp1,tmp2);;
  math::HelixToHelix(par1,par2,vtx1,vtx2,dis);
  if(dis<100*mm) return true; //mm
  else return false;
}

bool cds::CalcLineHelixVertex(LocalTrack *track,CDSTrack *cds, TVector3 &vtxline, TVector3 &vtxhelix, double &dis)
{
  double par[5];
  TVector3 tmpvtx;
  double tmp1,tmp2;
  cds->GetParameters(DetIdCDC,par,tmpvtx,tmp1,tmp2);
  TVector3 lpos=track->GetPosatZ(0);
  TVector3 ldir=track->GetMomDir();
  if(!math::LineToHelix(lpos,ldir,par,vtxline,vtxhelix,dis)) return false;
  if(dis<100*mm) return true; //mm
  else return false;
}

bool cds::CalcLineHelixVertex(pBeam *beam,CDSTrack *cds,TVector3 &vtx1,TVector3 &vtx2,double &dis, bool ELOSS)
{
  if(ELOSS){
    cds->GetVertex(beam->bpcpos(),beam->bpcdir(),vtx2,vtx1);
    dis=(vtx1-vtx2).Mag();
  }else{
    double par[5];
    TVector3 tmpvtx;
    double tmp1,tmp2;
    cds->GetParameters(DetIdCDC,par,tmpvtx,tmp1,tmp2);
    TVector3 lpos=beam->bpcpos();
    TVector3 ldir=beam->bpcdir();
    if(!math::LineToHelix(lpos,ldir,par,vtx1,vtx2,dis)) return false;
  }
  if(dis<100*mm) return true; //mm
  else return false;
}

void cds::CalcBetaMass(TVector3 vertex,LocalTrack *beam, CDSTrack* cdc, 
		       int beam_pid,double tof, double &beta,double &mass2)
{
  double param[5];
  cdc->GetParameters(param);
  if(param[2]==0 ||param[4]==0   ) return;
  double mom = cdc->mom();
  TVector3 gpos;
  //  conf->GetGeomMapManager()->GetPos(DetIdT0,0,gpos);
  double z_t0=gpos.Z(),z_vtxb=vertex.Z();
  TVector3 t0pos=beam->GetPosatZ(z_t0);
  TVector3 vpos=beam->GetPosatZ(z_vtxb);
  double beam_dis=(t0pos-vpos).Mag();

  double beta_b=mybeta(pdg::Mass(beam_pid),1.0); // assume 1.0 GeV/c beam momentum
  double time_beam=beam_dis/beta_b/(TMath::C()*1e-6*mm);
  
  //#####CDC dis#######		  	
  TVector3 cdhvtx=math::CalcHelixPosatR(param,cdh_r);
  double cdc_dis=math::CalcHelixArc(param,cdhvtx,vertex);

  double beta_calc=cdc_dis/(tof-time_beam)/(TMath::C()*1e-6*mm);
  double calc_mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
  
  beta = beta_calc;
  mass2 = calc_mass2;
#if 0
  std::cout<<"vertex.x,y,z"<<vertex.X()<<"\t"<<vertex.Y()<<"\t"<<vertex.Z()<<std::endl;
  std::cout<<"cdh x,y,z"<<cdhvtx.X()<<"\t"<<cdhvtx.Y()<<"\t"<<cdhvtx.Z()<<std::endl;
  std::cout<<"tof, beta, mass"<<tof<<"\t"<<beta<<"\t"<<mass<<std::endl;
#endif
}

int cds::PID1d(double mom,double mass2)
{
  int ptype=pdg::kOther;
  if(mom>0) {
    if(mass2<0.4) ptype=kPiPlus;
    else if(mass2<2.2) ptype=kProton;
    else if(mass2<5.) ptype=pdg::kDeuteron;
    else if(mass2<9.) ptype=pdg::kTriton;
    else ptype=pdg::kOther;
  }
  else{
    if(mass2<0.14) ptype=kPiMinus;
    else if(mass2<0.5) ptype=kKMinus;
    else ptype=pdg::kOther;
  }
  return ptype;
}

int cds::PID2d(double mom,double mass2)
{
  int ptype=pdg::kOther;
  //  double fmom=fabs(mom);
  /* 2016.07.05 -----> */
  double pi_mass2 = 0.139570*0.139570;
  double k_mass2  = 0.497234*0.497234;
  double p_mass2  = 0.918312*0.918312;
  double d_mass2  = 1.829990*1.829990;
  /* 2016.07.05 -----> */
  double p[4][4]={{0.00084734,0.00011146,pi_mass2,0.00674811},
		  {0.00084734,0.00011146,k_mass2 ,0.00674811},
		  {0.00084734,0.00011146,p_mass2 ,0.00674811},
		  {0.00084734,0.00011146,d_mass2 ,0.00674811}};
  /* <----- 2016.07.05 */
  
  double pi_sigma = sqrt(4.0*p[0][2]*p[0][0]*mom*mom+
			 4.0*p[0][2]*p[0][2]*p[0][1]*(1.0+p[0][2]/(mom*mom))+
			 4.0*p[0][3]*mom*mom*(p[0][2]+mom*mom));
  double k_sigma  = sqrt(4.0*p[1][2]*p[1][0]*mom*mom+
			 4.0*p[1][2]*p[1][2]*p[1][1]*(1.0+p[1][2]/(mom*mom))+
			 4.0*p[1][3]*mom*mom*(p[1][2]+mom*mom));
  double p_sigma  = sqrt(4.0*p[2][2]*p[2][0]*mom*mom+
			 4.0*p[2][2]*p[2][2]*p[2][1]*(1.0+p[2][2]/(mom*mom))+
			 4.0*p[2][3]*mom*mom*(p[2][2]+mom*mom));
  double d_sigma  = sqrt(4.0*p[3][2]*p[3][0]*mom*mom+
			 4.0*p[3][2]*p[3][2]*p[3][1]*(1.0+p[3][2]/(mom*mom))+
			 4.0*p[3][3]*mom*mom*(p[3][2]+mom*mom));
  /*=== positive ===*/
  double n_sigma = 2.5;
  double pi_ll = pi_mass2 - n_sigma*pi_sigma;
  double pi_ul = pi_mass2 + n_sigma*pi_sigma;
  double k_ll = k_mass2 - n_sigma*k_sigma;
  double k_ul = k_mass2 + n_sigma*k_sigma;
  double p_ll = p_mass2 - n_sigma*p_sigma;
  double p_ul = p_mass2 + n_sigma*p_sigma;
  double d_ll = d_mass2 - n_sigma*d_sigma;
  double d_ul = d_mass2 + n_sigma*d_sigma;
  if(mom>0){
    bool pi_flag = false;
    bool p_flag = false;
    bool d_flag = false;
    if(pi_ll<mass2&&mass2<pi_ul){
      pi_flag = true;
    }
    if(p_ll<mass2&&mass2<p_ul){
      if(mom>0.1){
        p_flag = true;
      }
    }
    if(d_ll<mass2&&mass2<d_ul){
      if(mom>0.2){
        d_flag = true;
      }
    }
    /*=== pi+ ===*/
    if(pi_flag){
      if(!p_flag){
        ptype = kPiPlus; return ptype;
      }
    }
    /*=== proton ===*/
    if(p_flag){
      if(!pi_flag&&!d_flag){
        ptype = kProton; return ptype;
      }
    }
    /*=== deuteron ===*/
    if(d_flag){
      if(!p_flag){
        ptype = pdg::kDeuteron; return ptype;
      }
    }
  }
  /*=== negative ===*/
  else {
    bool pi_flag = false;
    bool k_flag = false;
    if(pi_ll<mass2&&mass2<pi_ul){
      pi_flag = true;
    }
    if(k_ll<mass2&&mass2<k_ul){
      if(mom<-0.05){
        k_flag = true;
      }
    }
    /*=== pi- ===*/
    if(pi_flag){
      if(!k_flag){
        ptype = kPiMinus; return ptype;
      }
    }
    /*=== k- ===*/
    if(k_flag){
      if(!pi_flag){
        ptype = kKMinus; return ptype;
      }
    }
  }

  ptype = pdg::kOther; return ptype;
}

pCDS* cds::CalcSingleAll(pBeam *beam,CDSTrack* cdc,HodoAnalyzer *hodo,bool ELOSS,bool REFIT)
{
  //  std::cout<<"[cds::CalcSingleAll] start!!!"<<std::endl;
  TVector3 vtxcdc,vtxbpc;
  double par[5];
  cdc->GetParameters(DetIdCDC,par,vtxcdc);
  double mom=cdc->mom();
  double dis;
  cds::CalcLineHelixVertex(beam,cdc,vtxbpc,vtxcdc,dis);
  if(vtxcdc.Mag()<1e-10){
    std::cout<<" calc vertex failed"<<std::endl;
    return 0;
  }
  double time_beam=beam->CalcVertexTime(vtxbpc);

  //  std::cout<<"# new pCDS"<<std::endl;
  pCDS* cdstrack=new pCDS();
  cdstrack->SetParameters(par);
  cdstrack->SetTrackID(cdc->trackid());
  int pid=cdc->pid_tot();

  //  std::cout<<"# CDH matching"<<std::endl;
  TVector3 cdhvtx=math::CalcHelixPosatR(par,cdh_r);
  TVector3 cdhout=math::CalcHelixPosatR(par,cdh_rout);
  double cdc_dis=math::CalcHelixArc(par,cdhvtx,vtxcdc);
  double tofvtxcdc,cdhtime=999,cdhctime=999,diff=999;
  int cdhseg=-1;
  int nh = hodo->GetNHits(DetIdCDH);
  //  cdhvtx.Print();
  for( int i=0; i<nh; ++i ){
    Hodo2Hit *hit = hodo->GetHit(DetIdCDH,i);
    if(!hit) continue;
    int nind= hit->GetIndex();
    for(int ii=0;ii<nind;ii++){      
      double cmt  = hit->CMeanTime(ii);
      int seg=hit->SegmentId();
      if(!gUser.Check("CDHTOF",cmt-beam->t0ctime())) continue;
      TVector3 tmppos;
      if(!gGeom.GetGPos(DetIdCDH,seg,tmppos)) continue;
      //  tmppos.Print();
      double dphi1=TMath::Abs(tmppos.DeltaPhi(cdhvtx));
      double dphi2=TMath::Abs(tmppos.DeltaPhi(cdhout));
      if(dphi1>cdh_dphi&&dphi2>cdh_dphi) continue;
      cdstrack->SetCDHSeg(seg);
      cdstrack->AddDeCDH(hit->DeltaE());
      //      if(cmt<cdhtime){
      if(diff>fabs(0.5*(dphi1+dphi2))){
	diff=fabs(0.5*(dphi1+dphi2));
	cdhseg=seg;
	cdhtime=hit->MeanTime(ii);
	cdhctime=cmt;
      }
    }
  }//for(ihit)
  //  std::cout<<"# CDH analysis "<<cdstrack->ncdh()<<std::endl;
  double mass2=-1;
  if(cdstrack->ncdh()>0){
    double tof=cdhctime-beam->t0ctime();
    double beta=cdc_dis/(cdhctime-time_beam)/(TMath::C()*1e-6*mm);
    mass2=mom*mom*(1/(beta*beta)-1);
    if(ELOSS){
      cds::FindMass2(cdc,beam,tof,beta,mass2,tofvtxcdc);
    }
    pid=cds::PID1d(mom,mass2);
    if(REFIT){
      double cdc_mass=pdg::Mass(pid);      
      double beta_calc=fabs(mom)/sqrt(mom*mom+cdc_mass*cdc_mass);    
      double time_calc2=cdc_dis/(beta_calc*TMath::C()*1e-6*mm);      
      double cdhdt=tof-time_calc2-time_beam;
      cdc->Retiming(cdhdt,beta_calc,true);
      mom = cdc->mom();   
      cds::FindMass2(cdc,beam,tof,beta,mass2,tofvtxcdc);
    }
    cdstrack->SetTOF(cdhctime);
    cdstrack->SetBeta(beta);
    cdstrack->SetCDHSeg(cdhseg);
    cdstrack->SetMass(sqrt(mass2));
    cdstrack->SetMass2(mass2);
  }
  //  
  //  std::cout<<"# vertex with dE"<<std::endl;
  double mass=pdg::Mass(pid);
  if(mass<0) ELOSS=false;
  double tmpl;
  if(ELOSS){    
    //    std::cout<<"calcvertextimelength"<<std::endl;
    bool flag=cdc->CalcVertexTimeLength(beam->bpcpos(),beam->bpcdir(),mass,vtxbpc,vtxcdc,tofvtxcdc,tmpl,true);
    if(!flag){
      std::cout<<"[cds::CalcSingleAll] failed in CalcVertexTimeLength"<<std::endl;
      std::cout<<mass<<"  "<<cdc->mom()<<"  "<<tofvtxcdc<<"  "<<tmpl<<" ";
      vtxcdc.Print();
      delete cdstrack;
      return 0;
    }
  }  
  dis=(vtxbpc-vtxcdc).Mag();
  //  std::cout<<"# momentum with dE"<<std::endl;
  TVector3 Pd;
  double tof,length;
  if( !cdc->GetMomentum(vtxcdc, mass, Pd , tof,length, ELOSS) ){    
    std::cout<<"[cds::CalcSingleAll] failed in GetMomentum"<<std::endl;
    delete cdstrack;
    return 0;
  }      

  //  pid=PID2d(mom,mass2);

  cdstrack->SetVertexCDC(vtxcdc);
  cdstrack->SetVertexBeam(vtxbpc);
  cdstrack->SetVDis(dis);
  cdstrack->SetMomDir(Pd);
  cdstrack->SetPID(pid);
  cdstrack->SetRawMomentum(mom);
  cdstrack->SetMomentum(mom/TMath::Abs(mom)*Pd.Mag());
  cdstrack->SetAngleLab(beam->bpcdir().Angle(Pd));    
  cdstrack->SetChi2(cdc->chi2());
  cdstrack->SetFL(cdc_dis);
  cdstrack->SetDeCDC(cdc->totave2());

  TVector3 pos1;
  cdc->GetParameters(DetIdCDCCFRP,par,pos1);
  cdc->GetParameters(DetIdCDC,par,vtxcdc);
  cdc_dis=math::CalcHelixArc(par,cdhvtx,pos1);

  if(cdstrack->ncdh()>0){
    double ctime_beam=beam->CalcVertexTime(vtxbpc);
    time_beam=ctime_beam-beam->t0ctime()+beam->t0time();
    double beta=mybeta(mass,mom);
    double time_cdc=cdc_dis/beta/(TMath::C()*1e-6*mm);
    double time_cdc2=cdstrack->fl()/(TMath::C()*1e-6*mm);
    double cdcdt=(cdhtime-time_beam-tofvtxcdc)-time_cdc;
    double cdcdt2=(cdhtime-time_beam)-time_cdc2;
    double cdccdt=(cdhctime-ctime_beam-tofvtxcdc)-time_cdc;
    cdstrack->SetTOF(cdcdt2);
    cdstrack->SetFT(time_cdc);
    cdstrack->SetDt(cdcdt);
    cdstrack->SetCDt(cdccdt);
  }
  //  std::cout<<"[cds::CalcSingleAll] finished!!!"<<std::endl;
  return cdstrack;
}

pCDS *cds::Calc2HelixAll(pBeam *beam, pCDS* cds1, pCDS* cds2,CDCAnalyzer *cdcAna, bool ELOSS,int debug)
{
#if DEBUG
  std::cout<<"[cds::Calc2HelixAll] start!!!"<<std::endl;
#endif
  if(!cds1||!cds2) return 0;
  if(cds1->ncdh()==0||cds2->ncdh()==0){
    if(debug>0) std::cout<<"#W [cds::Calc2HelixAll] no cdh hit"<<std::endl;
    return 0;
  }
  CDSTrack *track1=cdcAna->GetTrack(cds1->id());	
  CDSTrack *track2=cdcAna->GetTrack(cds2->id());	
  TVector3 vtx1=DEFVECT,vtx2=DEFVECT;  
  //  if( ELOSS && !cds::CalcVertex2HelixdE(track1,track2,cds1->pdgmass(),cds2->pdgmass(),vtx1,vtx2) ){
  if( ELOSS && !cds::CalcVertex2Helix2(track1,track2,vtx1,vtx2) ){
    if(debug>0) std::cout<<"#W [cds::Calc2HelixAll] failed in CalcVertex2Helix2"<<std::endl;
    return 0;
  }else if( !ELOSS && !cds::CalcVertex2Helix(track1,track2,vtx1,vtx2) ){
    if(debug>0) std::cout<<"#W [cds::Calc2HelixAll] failed in CalcVertex2Helix"<<std::endl;
    return 0;
  }
  TVector3 Pp1,Pp2;
  double tof1,tof2,fl1,fl2;
  //  std::cout<<"===== GetMomentum"<<std::endl;
  if( !track1->GetMomentum(vtx1,cds1->pdgmass(),Pp1,tof1,fl1,ELOSS) ||
      !track2->GetMomentum(vtx2,cds2->pdgmass(),Pp2,tof2,fl2,ELOSS) ){
    if(debug>0){
      std::cout<<"#W [cds::Calc2HelixAll] failed in GetMomentum"<<std::endl;
      vtx1.Print();
      vtx2.Print();
    }
    return 0;
  }
  TVector3 Pp= Pp1+Pp2;
  TLorentzVector L1; L1.SetVectM( Pp1, cds1->pdgmass() );
  TLorentzVector L2; L2.SetVectM( Pp2, cds2->pdgmass() );
  double im = (L1+L2).M();	
  
  double dist,dltmp=0;
  TVector3 xest,nest;
  TVector3 vtx=(vtx1+vtx2)*0.5;
  math::LineToLine( vtx,Pp.Unit(),beam->bpcpos(), beam->bpcdir(),dltmp,dist,xest,nest );

  double time_beam=beam->CalcVertexTime(nest);
  double beta_pro=mybeta(im,Pp.Mag());
  double time_to_decay=(xest-vtx).Mag()/beta_pro/(TMath::C()*1e-6*mm);
  
  TVector3 pos1,pos2;
  double param[5];
  track1->GetParameters(DetIdCDCCFRP,param,pos1);
  track1->GetParameters(DetIdCDC,param,pos2);
  pos2=math::CalcHelixPosatR(param,cdh_r);
  double cdc_dis1=math::CalcHelixArc(param,pos1,pos2);
  track2->GetParameters(DetIdCDCCFRP,param,pos1);
  track2->GetParameters(DetIdCDC,param,pos2);
  pos2=math::CalcHelixPosatR(param,cdh_r);
  double cdc_dis2=math::CalcHelixArc(param,pos1,pos2);

  double time_cdc1=cdc_dis1/mybeta(cds1->pdgmass(),cds1->rawmom())/(TMath::C()*1e-6*mm);
  double time_cdc2=cdc_dis2/mybeta(cds2->pdgmass(),cds2->rawmom())/(TMath::C()*1e-6*mm);
  time_beam=0;
  double dt1=cds1->tof()-time_beam-time_to_decay-tof1-time_cdc1;
  double dt2=cds2->tof()-time_beam-time_to_decay-tof2-time_cdc2;
  //  std::cout<<time_beam<<"  "<<beam->t0ctime()<<std::endl;
  pCDS* pro=new pCDS();    
  pro->SetDaughterID1(cds1->id());
  pro->SetDaughterID2(cds2->id());
  pro->SetDaughterLMom1(L1);
  pro->SetDaughterLMom2(L2);
  pro->SetCombID((int)(pow(2,cds1->pid())+pow(2,cds2->pid())));
  pro->SetMomentum(Pp.Mag());
  pro->SetMass(im);
  pro->SetBeta(beta_pro);
  pro->SetGamma(1./sqrt(1-beta_pro*beta_pro));
  pro->SetVertex(vtx);
  pro->SetVertexBeam(nest);
  pro->SetVertexCDC(xest);
  pro->SetMomDir(Pp.Unit());
  pro->SetVDis((vtx1-vtx2).Mag()); // 2 helix
  pro->SetVBDis((xest-vtx).Mag()); // displaced vertex 
  pro->SetPBDCA(dist); // DCA btw beam and product 
  pro->SetOA(Pp1.Angle(Pp2));
  pro->SetAngleLab(beam->bpcdir().Angle(Pp));
  pro->SetDt1(dt1);
  pro->SetDt2(dt2);
  //  std::cout<<"[cds::Calc2HelixAll] finished!!!"<<std::endl;
  return pro;
}
void cds::PDFLambda(double* per, double* pdf,bool YAMAGA){
  /* Input : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  /*Output : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  double c0,c1,c2,c3,c4,c5;
  double tmppdf = 0;

  if(!YAMAGA){
    double cutparam[5][6]=
      {  1.0,   1.1155,    0.002, 0.0, 0.0, 0.0, //mass
	 1.0,-1.578198,-0.000075, 0.0, 0.0, 0.0, //dca pi-p
	 1.0,-3.078679,-0.000329, 0.0, 0.0, 0.0, //dca Lambda-K
	 1.0,-3.403185,-0.000388, 0.0, 0.0, 0.0, //dca p-K
	 2.967292, -2.445, -7.8313095, -6.826, -0.2258,  1.298 //dca4 Lambda-p
      };
    
    /* mass */
    c0=cutparam[0][0];
    c1=cutparam[0][1];
    c2=cutparam[0][2];
    tmppdf = c0 * TMath::Exp(-TMath::Power((per[0]-c1)/c2,2)/2.0);
    pdf[0] = tmppdf;
    
    for( int i=1; i<4; i++){
      c0=cutparam[i][0];
      c1=cutparam[i][1];
      c2=cutparam[i][2];
      tmppdf = c0 *TMath::Exp(c1*per[i]+c2*per[i]*per[i]);
      pdf[i] = tmppdf;
    }
    
    c0=cutparam[4][0];
    c1=cutparam[4][1];
    c2=cutparam[4][2];
    c3=cutparam[4][3];
    c4=cutparam[4][4];
    c5=cutparam[4][5];  
    tmppdf = c0*TMath::Exp(c1*per[4]) + c2*TMath::Exp(c3*TMath::Power((per[4]-c4),c5));
    pdf[4] = tmppdf;  
    return;
  }else{
    // double cutparam[5][6]=
    //   {  1.0,    1.1155,      0.002, 0.0, 0.0, 0.0, //mass
    // 	 1.0, -1.449435,  -0.772811, 0.0, 0.0, 0.0, //dca pi-p
    // 	 1.0, -2.600850, -10.915050, 0.0, 0.0, 0.0, //dca Lambda-K
    // 	 1.0, -2.721580, -15.428124, 0.0, 0.0, 0.0, //dca p-K
    // 	 1.0, -1.222023, 0.007465, 0.0, 0.0, 0.0, //dca Lambda-p
    //   };      
    double cutparam[5][6]={1.0,   1.1156,      0.002, 0.0, 0.0, 0.0,  //mass
			   1.0, -1.75127,  -0.229124, 0.0, 0.0, 0.0,  //dca pi-p
			   1.0, -4.01229,   -5.12492, 0.0, 0.0, 0.0,  //dca Lambda-K
			   1.0, -2.42777,   -14.5197, 0.0, 0.0, 0.0,  //dca p-K
			   1.0, -0.642945, -0.492661, 0.0, 0.0, 0.0};  //dca Lambda-p

    /* mass */
    c0=cutparam[0][0];
    c1=cutparam[0][1];
    c2=cutparam[0][2];
    tmppdf = c0 * TMath::Exp(-TMath::Power((per[0]-c1)/c2,2)/2.0);
    pdf[0] = tmppdf;
    
    for( int i=1; i<5; i++){
      c0=cutparam[i][0];
      c1=cutparam[i][1];
      c2=cutparam[i][2];
      tmppdf = c0 * TMath::Exp(c1*per[i]+c2*per[i]*per[i]);
      pdf[i] = tmppdf;
    }
    return;
  }
  //std::cout << "<------- PDFLambda ###" << std::endl;
}
