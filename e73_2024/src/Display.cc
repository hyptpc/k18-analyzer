#include <cmath>
#include "Display.hh"
#include <TF1.h>
#include <TH2F.h>
#include <TPolyLine.h>
#include <TArc.h>
#include <TLine.h>
#include <TMarker.h>
#include <TBox.h>

#include "DetectorID.hh"
#include "CDCWireMapMan.hh"
#include "GeomMapMan.hh"
#include "UserParamMan.hh"
#include "BLDCWireMapMan.hh"
#include "DetectorList.hh"
#include "DCCluster.hh"

#define YOKE 0

namespace{
  const DetectorList&   dList = DetectorList::GetInstance();
  const CDCWireMapMan&  gCDC  = CDCWireMapMan::GetInstance();
  GeomMapMan&     gGeom = GeomMapMan::GetInstance();
  const UserParamMan& gUser = UserParamMan::GetInstance();

  TArc arc       = TArc();
  TLine line     = TLine();
  TBox box       = TBox();
  TMarker marker = TMarker();
  TPolyLine pline= TPolyLine();

  void SetLine(int color,double width,int style=1){
    line.SetLineColor(color);
    line.SetLineStyle(style);
    line.SetLineWidth(width);
  }
  void SetMarker(int color,int style,double size=1){
    marker.SetMarkerColor(color);
    marker.SetMarkerStyle(style);
    marker.SetMarkerSize(size);
  }
  void SetArc(int color,double width=1,int style=1,int fill=0){
    arc.SetLineColor(color);
    arc.SetLineWidth(width);
    arc.SetLineStyle(style);
    arc.SetFillStyle(fill);
  }
}

bool disp::Wait()
{
  char dispflag;
  char dispin[100]="";
  std::cout << " Input any word or return" << std::endl;
  std::cout << " (q:quit)" << std::endl;
  fgets(dispin,100,stdin);
  if(sscanf(dispin,"%c",&dispflag)==1){
    if( dispflag=='q' ){
      return false;
    }
  }
  return true;
}
// ====================================================
bool disp::DrawPLine(const double &xc,const double &yc, const double &dx, const double &dy, const double &rot)
{
  double x[5],y[5];
  TVector3 tmppos[5];
  tmppos[0]=TVector3( dx/2, dy/2,0);
  tmppos[1]=TVector3(-dx/2, dy/2,0);
  tmppos[2]=TVector3(-dx/2,-dy/2,0);
  tmppos[3]=TVector3( dx/2,-dy/2,0);
  tmppos[4]=TVector3( dx/2, dy/2,0);
  for(int i=0;i<5;i++){
    tmppos[i].RotateZ(rot);
    x[i]=tmppos[i].X()+xc;
    y[i]=tmppos[i].Y()+yc;
  }
  pline.DrawPolyLine(5,x,y);
  return true;
}
// ====================================================
bool disp::DrawSegmentsXY( TVirtualPad *pad, int cid )
{
  if( cid!=DetIdCDH && cid != DetIdIH ) return false;
  pad->cd();
  int nseg=dList.GetNsegs(cid);
  for( int seg=0; seg<nseg; seg++ ){
    gGeom.GetXYCDS(cid,seg,pline);
  }
  return true;
}
// ====================================================
bool disp::DrawSegmentsZY( TVirtualPad *pad, int cid )
{
  pad->cd();
  if(cid==DetIdCDH){
    double xc=55.9,wid=3.,gxc=0.;
    double zc=0.,th=79.0,gzc=0.;
    double rot=0.;
    disp::DrawPLine(zc+gzc, xc+gxc,th,wid,rot);
    disp::DrawPLine(zc+gzc,-xc+gxc,th,wid,rot);
  }else if(cid==DetIdIH){
    double xc=13.9,wid=0.3,gxc=0.;
    double zc=0.,th=60.0,gzc=0.;
    double rot=0.;
    disp::DrawPLine(zc+gzc, xc+gxc,th,wid,rot);
    disp::DrawPLine(zc+gzc,-xc+gxc,th,wid,rot);
  }
  return true;
}
// ====================================================
bool disp::DrawSegmentsZR( TVirtualPad *pad, int cid )
{
  pad->cd();
  if(cid==DetIdCDH){
    double xc=55.9,wid=3.,gxc=0.;
    double zc=0.,th=79.0,gzc=0.;
    double rot=0.;
    disp::DrawPLine(zc+gzc, xc+gxc,th,wid,rot);
  }else if(cid==DetIdIH){
    double xc=13.9,wid=0.3,gxc=0.;
    double zc=0.,th=60.0,gzc=0.;
    double rot=0.;
    disp::DrawPLine(zc+gzc, xc+gxc,th,wid,rot);
  }
  return true;
}
// ====================================================
bool disp::DrawSegmentsZX( TVirtualPad *pad, int cid, bool GLOBAL )
{
  TPolyLine pl;
  pad->cd();
  if( cid==DetIdCVC || cid==DetIdT0 || cid==DetIdBHD
      //      || cid==DetIdPC || cid==DetIdNC
      ){
    int nsegs=dList.GetNsegs(cid);
    for( int seg=0; seg<nsegs; seg++ ){
      gGeom.GetZX(cid,seg,pline,GLOBAL);
    }
  }else if(cid==DetIdCDH || cid==DetIdIH){
    DrawSegmentsZY(pad,cid);
  }
  return true;
}
// ====================================================
bool disp::DrawCDCLayersXY( TVirtualPad *pad )
{
  pad->cd();
  Double_t rad,phi0,dphi,tilt;
  SetArc(1);
  for( Int_t lay=0; lay<15; lay++ ){
    gCDC.GetGeom( lay, rad, phi0, dphi, tilt );
    if( 0<tilt )      arc.SetLineColor(4);
    else if( tilt<0 ) arc.SetLineColor(6);
    else              arc.SetLineColor(3);
    //    std::cout<<lay<<"  "<<rad<<"  "<<tilt<<std::endl;
    arc.DrawArc( 0., 0., rad, 0., 360. );
  }
  SetArc(1,2);
  Double_t rin,rout,zlen;
  gCDC.GetFrame(zlen,rin,rout);
  arc.DrawArc( 0., 0., rin , 0., 360. );
  arc.DrawArc( 0., 0., rout, 0., 360. );
  return true;
}
// ====================================================
bool disp::DrawCDCLayersZY( TVirtualPad *pad, bool ZR )
{
  pad->cd();
  Double_t rin,rout,zlen;
  gCDC.GetFrame(zlen,rin,rout);
  line.SetLineWidth(1);
  Double_t rad,phi0,dphi,tilt;
  for( Int_t lay=0; lay<15; lay++ ){
    gCDC.GetGeom( lay, rad, phi0, dphi, tilt );
    if( 0<tilt ) line.SetLineColor(4);
    else if( tilt<0 ) line.SetLineColor(6);
    else line.SetLineColor(3);
    //    std::cout<<lay<<"  "<<rad<<"  "<<zlen<<std::endl;
    line.DrawLine(-zlen/2., rad, zlen/2., rad);
    if(!ZR) line.DrawLine(-zlen/2.,-rad, zlen/2.,-rad);
  }
  box.SetFillStyle(0);
  box.SetLineWidth(1);
  box.SetLineColor(1);
  box.DrawBox(-zlen/2., rin, zlen/2., rout);
  if(!ZR)  box.DrawBox(-zlen/2.,-rin, zlen/2.,-rout);
  return true;
}
// ====================================================
bool disp::DrawCDCLayersZX( TVirtualPad *pad )
{
  return DrawCDCLayersZY(pad);
}
// ====================================================
bool disp::DrawCDCLayersZR( TVirtualPad *pad )
{
  return DrawCDCLayersZY(pad,true);
}
// ====================================================
bool disp::DrawCDSHits( TVirtualPad *pad, HodoAnalyzer* hodo, int cid, enum XYZ xyz )
{
  if( cid!=DetIdCDH && cid!=DetIdIH) return false;
  pad->cd();
  SetMarker(1,20);
  int nh = hodo->GetNHits(cid);
  for( int i=0; i<nh; ++i ){
    Hodo2Hit *hit = hodo->GetHit(cid,i);
    if(!hit) continue;
    int nind= hit->GetIndex();
    for(int ii=0;ii<nind;ii++){
      double mt  = hit->MeanTime(ii);
      double tsub= hit->TimeDiff(ii);
      if(cid==DetIdCDH){
	if(!gUser.IsInRange("CDHTOF",mt)) continue;
      }else{
	if(!gUser.IsInRange("HodoGATE",mt)) continue;
      }
      int seg=hit->SegmentId();
      TVector3 pos,pos2;
      gGeom.GetGPos(cid,seg,tsub,pos);
      gGeom.GetGPos(cid,seg,pos2);
      if(xyz==kXY) marker.DrawMarker(pos.X(),pos.Y());
      if(xyz==kZX) marker.DrawMarker(pos.Z(),pos.Perp());
      if(xyz==kZY) marker.DrawMarker(pos.Z(),pos.Perp());
      if(xyz==kZR) marker.DrawMarker(pos.Z(),pos.Perp());
      if(xyz==kZPhi){
	marker.DrawMarker(pos.Z(),pos.Phi());
	if(pos.Phi()>3)  marker.DrawMarker(pos.Z(),pos.Phi()-2*TMath::Pi());
	if(pos.Phi()<-3) marker.DrawMarker(pos.Z(),pos.Phi()+2*TMath::Pi());
      }
      break;
    }
  }
  return true;
}
// ====================================================
bool disp::DrawCDCHits( TVirtualPad *pad, CDCAnalyzer* cdc, enum XYZ xyz )
{
  pad->cd();
  for( int layer=0; layer<15; layer++ ){
    const DCHitContainer &cont = cdc->GetDCHC(DetIdCDC, layer);
    int mul=cont.size();
    for(int ihit=0;ihit<mul;ihit++){
      DCHit* hit=cont[ihit];
      TVector3 pos=hit->GetWirePosition();
      TVector3 dir=hit->GetWireDirection();
      pos-=pos.Z()/dir.Z()*dir;
      SetMarker(1,1,1);
      if(xyz==kXY)   marker.DrawMarker(pos.X(),pos.Y());

      bool FLAG=false;
      int multi=hit->GetDriftLengthSize();
      for ( int m1=0; m1<multi; ++m1 ) {
	if( !hit->IsWithinTotRange(m1) ) continue;
	if( !hit->IsWithinDtRange(m1) )  continue;
	FLAG=true;
	break;
      }
      if(!FLAG) continue;
      SetMarker(1,24,0.5);
      SetLine(1,1,2);
      if(xyz==kXY)     marker.DrawMarker(pos.X(),pos.Y());
      if(xyz==kZPhi){
	if(dir.Perp()<1) continue;
	TVector3 pos1=pos+0.5*dir;
	TVector3 pos2=pos-0.5*dir;
	line.DrawLine(pos1.Z(),pos1.Phi(),pos2.Z(),pos2.Phi());
      }
    }
  }
  return true;
}
// ====================================================
bool disp::DrawCDCClusterHits( TVirtualPad *pad, CDCAnalyzer* cdc, enum XYZ xyz )
{
  pad->cd();
  for( int layer=0; layer<7; layer++ ){
    int ncl=cdc->GetNClusters(layer);
    for(int icl=0;icl<ncl;icl++){
      DCCluster* cl=cdc->GetCluster(layer,icl);
      int mul=cl->nhit();
      for(int ihit=0;ihit<mul;ihit++){
	DCHit* hit=cl->hit(ihit);
	TVector3 pos=hit->GetWirePosition();
	TVector3 dir=hit->GetWireDirection();
	pos-=pos.Z()/dir.Z()*dir;
	SetMarker(1,1,1);
	if(xyz==kXY)   marker.DrawMarker(pos.X(),pos.Y());

	bool FLAG=false;
	int multi=hit->GetDriftLengthSize();
	for ( int m1=0; m1<multi; ++m1 ) {
	  if( !hit->IsWithinTotRange(m1) ) continue;
	  if( !hit->IsWithinDtRange(m1) )  continue;
	  FLAG=true;
	  break;
	}
	if(!FLAG) continue;

	SetMarker(1,20,0.5);
	SetLine(1,1,2);
	if(xyz==kXY)     marker.DrawMarker(pos.X(),pos.Y());
	if(xyz==kZPhi){
	  if(dir.Perp()<1) continue;
	  TVector3 pos1=pos+0.5*dir;
	  TVector3 pos2=pos-0.5*dir;
	  line.DrawLine(pos1.Z(),pos1.Phi(),pos2.Z(),pos2.Phi());
	}
      }
    }
  }
  return true;
}
// ====================================================
bool disp::DrawCDCTracks( TVirtualPad *pad, CDCAnalyzer* cdc, enum XYZ xyz )
{
  pad->cd();
  for(int i=0;i<cdc->GetNTracks();i++){
    SetMarker(2+i,20,0.5);
    CDSTrack *tr=cdc->GetTrack(i);
    // hits
    for(int layer=0;layer<15;layer++){
      for(int nhit=0;nhit<tr->nhits(layer);nhit++){
	TVector3 pos=tr->hit(layer,nhit)->wpos;
	TVector3 dir=tr->hit(layer,nhit)->wdir;
	pos-=pos.Z()/dir.Z()*dir;
	if(xyz==kXY)  marker.DrawMarker(pos.X(),pos.Y());
	pos=tr->hit(layer,nhit)->hitpos;
	if(xyz==kZX)  marker.DrawMarker(pos.Z(),pos.X());
	if(xyz==kZY)  marker.DrawMarker(pos.Z(),pos.Y());
	if(xyz==kZR)  marker.DrawMarker(pos.Z(),pos.Perp());
	if(xyz==kZPhi){
	  double phi=pos.Phi();
	  marker.DrawMarker(pos.Z(),phi);
	  if(phi>3)  marker.DrawMarker(pos.Z(),phi-2*TMath::Pi());
	  if(phi<-3) marker.DrawMarker(pos.Z(),phi+2*TMath::Pi());
	}
      }
    }
    //
    if(tr->IsLine()){ // straight
      if(xyz==kXY){
	SetLine(2+i,1);
	TVector3 p0=tr->GetPos0();
	TVector3 p1=tr->GetPositionatR(60);
	line.DrawLine(p0.X(),p0.Y(),p1.X(),p1.Y());
      }
    }else{
      TVector3 circle=TVector3(tr->CircleX(),tr->CircleY(),0);
      if(xyz==kXY){
	SetArc(2+i);
	arc.DrawArc(tr->CircleX(),tr->CircleY(),tr->CircleR(),0,360);
      }
      {
	SetLine(2+i,1);
	SetMarker(2+i,24,1);
	TVector3 posin =tr->GetPositionatR(10);
	TVector3 posout=tr->GetPositionatR(60);
	if(xyz==kXY){
	  marker.DrawMarker(posin.X() ,posin.Y());
	  marker.DrawMarker(posout.X(),posout.Y());
	}
	if(xyz==kZR)   line.DrawLine(posin.Z(),posin.Perp(),posout.Z(),posout.Perp());
	if(xyz==kZPhi) line.DrawLine(posin.Z(),posin.Phi(),posout.Z(),posout.Phi());
      }
    }
    //CDH
    // TVector3 cdhpos;
    // for(int n=0;n<track->nCDHHit();n++)
    //   {
    // 	cdhpos.SetXYZ(track->CDHHit(cdsMan,n)->x(),
    // 		      track->CDHHit(cdsMan,n)->y(),
    // 		      track->CDHHit(cdsMan,n)->z() );
    // 	SetMarker(2+i,20,0.7);
    // 	marker.DrawMarker(cdhpos.x(),cdhpos.y());
    //   }
  }
  return true;
}
#if 0
// ====================================================
bool disp::DrawCDCTrackYZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  for(int i=0;i<trackMan->nTrack();i++)
    {
      CDSTrack *track=trackMan->Track(i);
      double aparam[5];
      track->GetParameters(aparam);

      if(aparam[2]==0)
	{
	  TLine lzy;
	  double a,b,c;
	  a=track->A(); b=track->B(); c=track->C();
	  double x1,y1,z1,x2,y2,z2;

	  x1=aparam[0]*cos(aparam[1])-(0)*sin(aparam[1]);
	  y1=aparam[0]*sin(aparam[1])+(0)*cos(aparam[1]);
	  z1=aparam[3]-(0)*aparam[4];
	  x2=aparam[0]*cos(aparam[1])-60*sin(aparam[1]);
	  y2=aparam[0]*sin(aparam[1])+60*cos(aparam[1]);
	  z2=aparam[3]-60*aparam[4];

	  lzy.SetLineColor(2+i);
	  lzy.DrawLine(z1,y1,z2,y2);
	}
      else
	{

	  TF1 *func_y=new TF1("func_y","([0])*sin([1])+1./[2]*(sin([1])-sin( [1]+([3]-x)/([4]*1./[2]) ) )",aparam[3]-1./aparam[2]*aparam[4]*(-TMath::Pi()/2. ),aparam[3]-1./aparam[2]*aparam[4]*(TMath::Pi()/2. ) );
	  func_y->SetParameters(aparam[0],aparam[1],aparam[2],aparam[3],aparam[4]);
	  func_y->SetLineColor(2+i);
	  func_y->SetLineWidth(1);

	  func_y->Draw("LPSAME");
	}
    }
  return true;
}
// ====================================================
bool disp::DrawCDCTrackXZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan )
{
  pad->cd();
  double zlen = conf->GetCDCWireMapManager()->zlen();
  //  double rin = conf->GetCDCWireMapManager()->rin();
  double rout = conf->GetCDCWireMapManager()->rout();
  double dz = zlen/500.;
  double xx[1000], yy[1000], zz[1000];
  TPolyLine pline;
  pline.SetLineColor(2);
  for( int itr=0; itr<trackMan->nTrack(); itr++ ){
    int count=0;
    for( int i=0; i<=500; i++ ){
      double z = -zlen/2. + dz*i;
      double x,y;
      trackMan->Track(itr)->XYatZ(x,y,z);
      double r = sqrt(x*x+y*y);
      if( r<rout ){
	xx[count]=x;
	yy[count]=y;
	zz[count]=z;
	count++;
      }
    }
    pline.DrawPolyLine(count,zz,xx);
  }
  return true;
}
// ====================================================
bool disp::DrawBLHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid )
{

  if( cid==DetIdCVC || cid==DetIdT0 || cid==DetIdBHD  ||
      cid==DetIdPC || cid==DetIdNC
      ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(2);
    int numseg=DetectorList::GetInstance()->GetNsegs(cid);

    for( int seg=0; seg<numseg; seg++ ){
      HodoscopeLikeHit *hod = bl->Hodo( cid, seg );
      if( hod==0 ) continue;
      if( !hod->CheckRange() ) continue;
      TVector3 pos=hod->pos();
      mark.DrawMarker(pos.Z(),pos.X());
    }
    return true;
  }

  return false;
}
// ====================================================
bool disp::DrawBLHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid )
{

  if( cid==DetIdCVC || cid==DetIdT0 || cid==DetIdBHD
      ){
    pad->cd();
    TMarker mark;
    mark.SetMarkerStyle(20);
    mark.SetMarkerColor(1);
    int numseg=DetectorList::GetInstance()->GetNsegs(cid);
    for( int seg=0; seg<numseg; seg++ ){
      HodoscopeLikeHit *hod = bl->Hodo( cid, seg );
      if( hod==0 ) continue;
      if( !hod->CheckRange() ) continue;
      TVector3 pos=hod->pos();
      mark.DrawMarker(pos.Z(),pos.Y());
    }
    return true;
  }

  return false;
}
//#######for BLDC  ########
// ====================================================
bool disp::DrawBLDCLayersXZ( TVirtualPad *pad, ConfMan *conf,int cid )
{
  pad->cd();
  TMarker mark;
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;
  for( int layer=0; layer<8; layer++ ){
    if(!conf->GetBLDCWireMapManager()->GetParam( cid, layer,
						 nw, z, xy, xy0, dxy, wl, tilt, ra )) continue;
    // std::cout<<"conf "<<lay<<rad<<tilt<<std::endl;
    if( xy!=0 ) continue;
    mark.SetMarkerStyle(7);

    for(int wire=0;wire<nw;wire++)
      {
	mark.DrawMarker( z,xy0+wire*dxy );
      }
  }
  return true;
}
// ====================================================
bool disp::DrawBLDCLayersYZ( TVirtualPad *pad, ConfMan *conf,int cid )
{
  pad->cd();
  TMarker mark;
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;
  for( int layer=1; layer<=8; layer++ ){
    conf->GetBLDCWireMapManager()->GetParam( cid, layer,
		 nw, z, xy, xy0, dxy, wl, tilt, ra );
    // std::cout<<"conf "<<lay<<rad<<tilt<<std::endl;
    if( xy!=1 ) continue;
    mark.SetMarkerStyle(7);
    for(int wire=0;wire<nw;wire++)
      {
	mark.DrawMarker( z,xy0+wire*dxy );
      }
  }
  return true;
}
// ====================================================
bool disp::DrawBLDCHit( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *blMan, int cid,int xy )
{
  //xz plane :xy=0
  //yz plane :xy=1
  if( !(xy==0 || xy==1) ) return false;
  if( !( cid==DetIdBLC1a || cid==DetIdBLC1b || cid==DetIdBLC1 || cid==DetIdBLC2 ||cid==DetIdBLC2a || cid==DetIdBLC2b || cid==DetIdBPC ) ) return false;

  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerSize(0.5);
  for( int layer=0; layer<8; layer++ )
    {
      int nw,xytmp;
      double z,xy0,dxy,wl,tilt,ra;
      conf->GetBLDCWireMapManager()->GetParam( cid, layer,
					       nw, z, xytmp, xy0,
					       dxy, wl, tilt, ra );

      if(xytmp!=xy)	  continue;
      for( int i=0; i<blMan->nBLDC(cid,layer); i++ )
	{
	  double x,z;
	  ChamberLikeHit *hit=blMan->BLDC(cid,layer,i);
	  //	  std::cout<<"lay,wire,dl= "<<layer<<"\t"<<hit->wire()<<"\t"<<hit->dl()<<"\t"<<hit->CheckRange()<<std::endl;
	  if(hit->CheckRange())   mark.SetMarkerColor(2);
	  else   mark.SetMarkerColor(3);
	  if(xy==0)   x=hit->wx();
	  else    x=hit->wy();
	  z=hit->wz();
	  mark.DrawMarker(z,x);
	}
    }

  return true;
}
// ====================================================
bool disp::DrawBLDCTrack( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan* blMan, BeamLineTrackMan *blTrackMan, int cid,int xy )
{
  //xz plane :xy=0
  //yz plane :xy=1
  if( !(xy==0 || xy==1) ) return false;
  if( !( cid==DetIdBLC1a || cid==DetIdBLC1b || cid==DetIdBLC1 || cid==DetIdBLC2 ||cid==DetIdBLC2a || cid==DetIdBLC2b || cid==DetIdBPC ) ) return false;

  pad->cd();
  for(int itr=0;itr<blTrackMan->ntrackBLDC(cid);itr++)
    {
      LocalTrack *bldc=blTrackMan->trackBLDC(cid,itr);
      // BLDCWireMapMan *BLwireman=conf->GetBLDCWireMapManager();
      // TVector3 gpos,gdir;
      // BLwireman->GetGParam(DetIdBLC2a,gpos,gdir);
      double z1,z2,x1,x2,y1,y2;
      z1=-25;z2=25;
      //      bldc->XYLocalPosatZ(z1+gpos.z(),x1,y1);
      //      bldc->XYLocalPosatZ(z2+gpos.z(),x2,y2);
      bldc->XYLocalPosatZ(z1,x1,y1);
      bldc->XYLocalPosatZ(z2,x2,y2);

      TLine line;
      line.SetLineStyle(2);
      line.SetLineColor(itr+2);
      if(xy==0)       line.DrawLine(z1,x1,z2,x2);
      else if(xy==1)  line.DrawLine(z1,y1,z2,y2);
      double a,b,c;
      bldc->abc(a,b,c);

      for(int ih=0;ih<bldc->nhit();ih++)
   	{
   	  ChamberLikeHit *hit=bldc->hit(blMan,ih);
   	  TMarker blt_m;
   	  double hpos,wpos,dltrack;
   	  if(hit->xy()==0)
   	    {
   	      hpos=hit->x(); wpos=hit->wx();dltrack=fabs(hpos-wpos);
	      if(xy==0)
		{
		  blt_m.SetMarkerStyle(20);
		  blt_m.SetMarkerSize(0.2);
		  blt_m.SetMarkerColor(3+itr );
		  blt_m.DrawMarker(hit->wz(),hit->wx() );
		}
   	    }
   	  else if(hit->xy()==1)
   	    {
   	      hpos=hit->y(); wpos=hit->wy();dltrack=fabs(hpos-wpos);
	      if(xy==1)
		{
		  blt_m.SetMarkerStyle(20);
		  blt_m.SetMarkerSize(0.2);
		  blt_m.SetMarkerColor(3+itr );
		  blt_m.DrawMarker(hit->wz(),hit->wy() );
		}
   	    }
   	}
     }
  return true;
}
// ====================================================
bool disp::DrawBLDCTrackHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id)
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(4);
  mark.SetMarkerSize(0.5);

  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in X plane"<<std::endl;
  for( int itr=0; itr<ntr; itr++ ){
    LocalTrack *track=trackMan->trackBLDC(id,itr);
    std::cout<<"track "<<itr<<"\t chi2xz= "<<track->chi2xz()
	     <<"\t chi2= "<<track->chi2all()
	     <<"\t tracktimeX= "<<track->GetTrackTime(0)
	     <<"\t tracktime= "<<track->GetTrackTime()
	     <<"\t RMS= "<<track->GetTrackTimeRMS()<<std::endl;
    for( int ih=0; ih<track->nhit(0); ih++ ){
      double x,y,z;
      ChamberLikeHit *hit=track->hit(0,ih);
      if(hit->leftright())
	x = hit->wx()-hit->dl();
      else
	x = hit->wx()+hit->dl();
      y = hit->y();
      z = hit->z();
      mark.DrawMarker(z,x);
    }
  }
  return true;
}
// ====================================================
bool disp::DrawBLDCTrackXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in X plane"<<std::endl;

  pline.SetLineColor(col);
  for( int itr=0; itr<ntr; itr++ ){
    if(col==-1)
      pline.SetLineColor(itr+2);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
    pline.DrawPolyLine(2,z,x);
  }
  return true;
}
// ====================================================
bool disp::DrawBLDCTrackHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan* blMan, BeamLineTrackMan *trackMan , const int &id)
{
  pad->cd();
  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(4);
  mark.SetMarkerSize(0.5);

  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in Y plane"<<std::endl;

  for( int itr=0; itr<ntr; itr++ ){
    LocalTrack *track=trackMan->trackBLDC(id,itr);
    std::cout<<"track "<<itr<<"\t chi2yz= "<<track->chi2yz()
	     <<"\t chi2= "<<track->chi2all()
	     <<"\t tracktimeY= "<<track->GetTrackTime(1)
	     <<"\t tracktime= "<<track->GetTrackTime()<<std::endl;
    for( int ih=0; ih<track->nhit(1); ih++ ){
      double x,y,z;
      ChamberLikeHit *hit=track->hit(blMan,1,ih);
      x = hit->x();
      if(hit->leftright())
	y = hit->wy()-hit->dl();
      else
	y = hit->wy()+hit->dl();
      z = hit->z();
      mark.DrawMarker(z,y);
    }
  }
  return true;
}
// ====================================================
bool disp::DrawBLDCTrackYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  int ntr=trackMan->ntrackBLDC(id);
  std::cout<<id<<"\ttracking status: "<<trackMan->status(id)<<"\t"<<ntr<<" track in Y plane"<<std::endl;

  pline.SetLineColor(col);
  for( int itr=0; itr<ntr; itr++ ){
    if(col==-1)
      pline.SetLineColor(itr+2);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[0],x[0],y[0]);
    trackMan->trackBLDC(id,itr)->XYLocalPosatZ(z[1],x[1],y[1]);
    pline.DrawPolyLine(2,z,y);
  }
  return true;
}
// ====================================================
bool disp::DrawClusterTime( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *trackMan , const int &id, const int &col)
{
  pad->cd();
  double x,y;
  int ntr=trackMan->ntrackBLDC(id);

  TMarker mark;
  mark.SetMarkerSize(0.5);
  mark.SetMarkerStyle(20);

  for( int itr=0; itr<ntr; itr++ ){
    LocalTrack *track=trackMan->trackBLDC(id,itr);
    mark.SetMarkerColor(itr+2);
    for(int xy=0;xy<2;xy++)
      for(int i=0;i<track->ncluster(xy);i++){
	x=track->cluster(xy,i)->hit(blMan,0)->z();
	y=track->cluster(xy,i)->GetCTime();
	mark.DrawMarker(x,y);
      }
  }
  return true;
}
// ====================================================
bool disp::DrawBLC2TrackfromBLC1YZ( TVirtualPad *pad, ConfMan *conf, LocalTrack *track, const double &mom,const int &col)
{
  /*
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  double a,b,c,d,e,f;
  track->gabc(a,b,c);
  track->gdef(d,e,f);

  double parblc1[6];
  parblc1[0]=c; //cm
  parblc1[1]=TMath::ATan(b)*1000;//mrad
  parblc1[2]=f; //cm
  parblc1[3]=TMath::ATan(e)*1000;//mrad
  parblc1[4]=0.;
  parblc1[5]=(mom-1000)/1000.*100.;
  double mat[36];

  conf->GetTransferMatrixManager()->GetD5Matrix(mat);
  double parblc[6];
  for(int i=0;i<6;i++){
    parblc[i]=0;
    for(int j=0;j<6;j++){
      parblc[i]+=parblc1[j]*mat[6*i+j];
    }
  }
  TVector2 pos(parblc[0],parblc[2]);
  TVector2 dir(parblc[1],parblc[3]);

  TVector2 pos2=pos.Rotate(-135./180.*TMath::Pi());
  TVector2 dir2=dir.Rotate(-135./180.*TMath::Pi());

  y[0]=pos2.Y()+z[0]*TMath::Tan(dir2.Y()/1000.);
  y[1]=pos2.Y()+z[1]*TMath::Tan(dir2.Y()/1000.);
  TPolyLine pline;
  pline.SetLineColor(col);
  pline.DrawPolyLine(2,z,y);
  */
  return true;
}
// ====================================================
bool disp::DrawBLC2TrackfromBLC1XZ( TVirtualPad *pad, ConfMan *conf, LocalTrack *track, const double &mom,const int &col)
{
  /*
  pad->cd();
  double x[2],y[2];
  double z[2]={-25,25};
  double a,b,c,d,e,f;
  track->gabc(a,b,c);
  track->gdef(d,e,f);

  double parblc1[6];
  parblc1[0]=c; //cm
  parblc1[1]=TMath::ATan(b)*1000;//mrad
  parblc1[2]=f; //cm
  parblc1[3]=TMath::ATan(e)*1000;//mrad
  parblc1[4]=0.;
  parblc1[5]=(mom-1000)/1000.*100.;
  double mat[36];

  conf->GetTransferMatrixManager()->GetD5Matrix(mat);
  double parblc[6];
  for(int i=0;i<6;i++){
    parblc[i]=0;
    for(int j=0;j<6;j++){
      parblc[i]+=parblc1[j]*mat[6*i+j];
    }
  }
  TVector2 pos(parblc[0],parblc[2]);
  TVector2 dir(parblc[1],parblc[3]);

  TVector2 pos2=pos.Rotate(-135./180.*TMath::Pi());
  TVector2 dir2=dir.Rotate(-135./180.*TMath::Pi());

  x[0]=pos2.X()+z[0]*TMath::Tan(dir2.X()/1000.);
  x[1]=pos2.X()+z[1]*TMath::Tan(dir2.X()/1000.);
  pline.SetLineColor(col);
  pline.DrawPolyLine(2,z,x);
  */
  return true;
}
// ====================================================
bool disp::DrawBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCLayersXZ(pad,conf,cid);
  DrawBLDCHit(pad,conf,bl,cid,0);
  DrawBLDCTrackHitXZ(pad,conf,track,cid);
  DrawBLDCTrackXZ(pad,conf,track,cid);
  return true;
}
// ====================================================
bool disp::DrawTrackBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCTrackHitXZ(pad,conf,track,cid);
  DrawBLDCTrackXZ(pad,conf,track,cid,3);
  return true;
}
// ====================================================
bool disp::DrawBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCLayersYZ(pad,conf,cid);
  DrawBLDCHit(pad,conf,bl,cid,1);
  DrawBLDCTrackHitYZ(pad,conf,bl,track,cid);
  DrawBLDCTrackYZ(pad,conf,track,cid);
  return true;
}
// ====================================================
bool disp::DrawTrackBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
			  BeamLineTrackMan *track,int cid)
{
  DrawBLDCTrackHitYZ(pad,conf,bl,track,cid);
  DrawBLDCTrackYZ(pad,conf,track,cid,3);
  return true;
}
// ====================================================
bool disp::DrawFDC(TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl )
{
  DrawFDCWire( pad, conf );
  DrawFDCHitWire( pad, conf ,bl);
  return true;
}
// ====================================================
bool disp::DrawFDCWire( TVirtualPad *pad, ConfMan *conf )
{
  pad->cd();
  BLDCWireMapMan *wireman=conf->GetBLDCWireMapManager();

  TLine line;
  line.SetLineStyle(2);
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;

  for(int layer=0;layer<3;layer++){
    wireman->GetParam(DetIdFDC1,layer*2+1,
		      nw,z,xy,xy0,dxy,wl,tilt,ra);
    TVector2 tmp(0,wl/2);
    TVector2 rot=tmp.Rotate(TMath::Pi()*ra/180.);
    for(int wire=0;wire<nw;wire++){
      TVector2 zero(xy0+dxy*(wire-1),0);
      TVector2 rot2=zero.Rotate(TMath::Pi()*ra/180.);

#if 0
      line.DrawLine( (zero+rot).X(), (zero+rot).Y(),
		     (zero-rot).X(), (zero-rot).Y() );
#else
      line.DrawLine( (rot2+rot).X(), (rot2+rot).Y(),
		     (rot2-rot).X(), (rot2-rot).Y() );
#endif
    }
  }
  return true;
}
// ====================================================
bool disp::DrawFDCHitWire( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl )
{
  pad->cd();
  BLDCWireMapMan *wireman=conf->GetBLDCWireMapManager();

  TLine line;
  line.SetLineWidth(2);
  int nw,xy;
  double z,xy0,dxy,wl,tilt,ra;

  for(int layer=0;layer<6;layer++){
    line.SetLineColor(layer+1);
    wireman->GetParam(DetIdFDC1,layer,
		      nw,z,xy,xy0,dxy,wl,tilt,ra);
    TVector2 tmp(0,wl/2);
    TVector2 rot=tmp.Rotate(TMath::Pi()*ra/180.);
    for(int i=0;i<bl->nBLDC(DetIdFDC1,layer);i++){
      int wire=bl->BLDC(DetIdFDC1,layer,i)->wire();
      TVector2 zero(xy0+dxy*(wire-1),0);
      line.DrawLine( (zero+rot).X(), (zero+rot).Y(),
		     (zero-rot).X(), (zero-rot).Y() );
    }
  }
  return true;
}
// ====================================================
#endif
